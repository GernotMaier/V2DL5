""""
Main analysis class.
"""

import logging
from pathlib import Path

import astropy.units as u
import numpy as np
import yaml
from gammapy.datasets import Datasets, SpectrumDataset
from gammapy.estimators import FluxPointsEstimator, LightCurveEstimator
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    PowerLawSpectralModel,
    SkyModel,
)
from IPython.display import display

import v2dl5.plot as v2dl5_plot
import v2dl5.time as v2dl5_time


class Analysis:
    """
    Analysis class.

    Parameters
    ----------
    configuration : str
        Analysis configuration
    output_dir : str
        Output directory
    sky_regions: SkyRegions
        Sky regions
    data : Data
        Datastore object

    """

    def __init__(self, args_dict=None, sky_regions=None, v2dl5_data=None):
        self._logger = logging.getLogger(__name__)

        self.args_dict = args_dict
        self._output_dir = args_dict.get("output_dir", None)
        self.sky_regions = sky_regions
        self.v2dl5_data = v2dl5_data

        self.datasets = None
        self.spectral_model = None
        self.flux_points = None
        self.light_curve_per_obs = None
        self.light_curve_per_night = None

    def run(self):
        """
        Run analysis.

        """

        self._data_reduction()
        if self.args_dict["datasets"]["stack"]:
            _data_sets = Datasets(self.datasets).stack_reduce()
        else:
            _data_sets = self.datasets
        self._define_spectral_models(model=self.args_dict["fit"]["model"])
        self._spectral_fits(datasets=_data_sets)
        self._flux_points(_data_sets)
        self.light_curve_per_obs = self._light_curve(_data_sets, None)
        self.light_curve_per_night = self._nightly_light_curve(_data_sets)

    def plot(self):
        """
        Plot all results.

        """

        _plot_dir = Path(f"{self._output_dir}/plots")

        plotter = v2dl5_plot.Plot(
            v2dl5_data=self.v2dl5_data,
            data_set=self.datasets,
            on_region=self.sky_regions.on_region,
            output_dir=_plot_dir,
        )
        plotter.plot_event_histograms()
        plotter.plot_irfs()
        plotter.plot_maps(exclusion_mask=self.sky_regions.exclusion_mask)
        plotter.plot_source_statistics()
        plotter.plot_spectra(
            flux_points=self.flux_points,
            model=self.spectral_model,
        )
        plotter.plot_light_curves(
            light_curve_per_obs=self.light_curve_per_obs,
            light_curve_per_night=self.light_curve_per_night,
        )

    def write(self):
        """
        Write all results.

        """

        for dataset in self.datasets:
            self._write_datasets(dataset, f"{dataset.name}.fits.gz")
        self._write_datasets(self.flux_points, "flux_points.ecsv", "gadf-sed")
        self._write_datasets(
            self.light_curve_per_obs, "light_curve_per_obs.ecsv", "lightcurve", "flux"
        )
        self._write_datasets(
            self.light_curve_per_night, "light_curve_per_night.ecsv", "lightcurve", "flux"
        )
        if self.spectral_model:
            self._write_yaml(self.spectral_model.to_dict(), "spectral_model.yaml")

    def _write_datasets(self, datasets, filename, file_format=None, sed_type=None):
        """
        Write datasets to disk.

        Parameters
        ----------
        datasets : Datasets
            Datasets
        filename : str
            Filename
        file_format: str
            Format specification (gammapy)

        """

        if datasets is None:
            return

        _out_file = f"{self._output_dir}/data/{filename}"
        self._logger.info("Writing datasets to %s", _out_file)
        if file_format is not None:
            if sed_type is not None:
                datasets.write(_out_file, overwrite=True, format=file_format, sed_type=sed_type)
            else:
                datasets.write(_out_file, overwrite=True, format=file_format)
        else:
            datasets.write(_out_file, overwrite=True)

    def _write_yaml(self, data_dict, filename):
        """
        Write model to disk.

        Parameters
        ----------
        data_dict : dict
            Data dictionary
        filename : str
            Filename

        """

        _data_dir = Path(f"{self._output_dir}/data")
        _data_dir.mkdir(parents=True, exist_ok=True)
        _out_file = f"{_data_dir}/{filename}"
        self._logger.info("Writing dataset to %s", _out_file)
        with open(_out_file, "w", encoding="utf-8") as outfile:
            yaml.dump(data_dict, outfile, default_flow_style=False)

    def _data_reduction(self):
        """
        Reduce data using the reflected region maker.

        """

        energy_axis = self._get_energy_axis(name="energy")
        energy_axis_true = self._get_energy_axis(name="energy_true")

        geom = RegionGeom.create(region=self.sky_regions.on_region, axes=[energy_axis])

        dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

        dataset_maker = SpectrumDatasetMaker(
            containment_correction=False, selection=["counts", "exposure", "edisp"]
        )
        bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=self.sky_regions.exclusion_mask)
        safe_mask_masker = SafeMaskMaker(
            methods=self.args_dict["datasets"]["safe_mask"]["methods"],
            aeff_percent=self.args_dict["datasets"]["safe_mask"]["parameters"]["aeff_percent"],
        )
        self._logger.info(
            "Mask applied: %s, aeff_percent = %d",
            safe_mask_masker.methods,
            safe_mask_masker.aeff_percent,
        )

        self.datasets = Datasets()

        for obs_id, observation in zip(self.v2dl5_data.runs, self.v2dl5_data.get_observations()):
            dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
            dataset_on_off = bkg_maker.run(dataset, observation)
            dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
            self.datasets.append(dataset_on_off)

        print("Run-wise results:")
        self._print_results(self.datasets.info_table(cumulative=False))
        print("Cumulative results:")
        self._print_results(self.datasets.info_table(cumulative=True))

    def _print_results(self, info_table):
        """
        Print results per run.

        Parameters
        ----------
        info_table : Table
            Table with results per run

        """

        for col in info_table.colnames:
            if info_table[col].dtype.kind == "f":
                info_table[col].format = "{:.2f}"

        print(
            info_table[
                "name",
                "livetime",
                "ontime",
                "counts",
                "counts_off",
                "background",
                "alpha",
                "excess",
                "sqrt_ts",
            ]
        )
        print()

    def _spectral_fits(self, datasets=None):
        """
        Spectral fitting.

        """

        datasets.models = self.spectral_model

        _fit = Fit()
        result_joint = _fit.run(datasets=datasets)
        print(result_joint)
        display(result_joint.models.to_parameters_table())

    def _define_spectral_models(self, model=None):
        """ "
        Spectral models

        """

        _spectral_model = None
        if model == "pl":
            _spectral_model = PowerLawSpectralModel(
                amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
                index=2,
                reference=1 * u.TeV,
            )
        elif model == "ecpl":
            _spectral_model = ExpCutoffPowerLawSpectralModel(
                amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
                index=2,
                lambda_=0.1 * u.Unit("TeV-1"),
                reference=1 * u.TeV,
            )

        self.spectral_model = SkyModel(spectral_model=_spectral_model, name=model)

    def _flux_points(self, datasets):
        """
        Calculate flux points.

        """

        energy_edges = (
            np.geomspace(
                u.Quantity(self.args_dict["flux_points"]["energy"]["min"]).value,
                u.Quantity(self.args_dict["flux_points"]["energy"]["max"]).value,
                self.args_dict["flux_points"]["energy"]["nbins"],
            )
            * u.TeV
        )

        self._logger.info("Estimating flux points")
        fpe = FluxPointsEstimator(energy_edges=energy_edges, selection_optional="all")
        self.flux_points = fpe.run(datasets=datasets)

        display(self.flux_points.to_table(sed_type="dnde", formatted=True))
        display(self.flux_points.to_table(sed_type="e2dnde", formatted=True))

    def _light_curve(self, datasets, time_intervals=None):
        """
        Calculate light curve.

        Parameters
        ----------
        datasets : Datasets
            Datasets
        time_intervals : `~astropy.time.Time`
            Time intervals

        Returns
        -------
        _light_curve : `~gammapy.estimators.LightCurve`
            Light curve

        """

        self._logger.info(f"Estimating light curve {time_intervals}")

        lc_maker_1d = LightCurveEstimator(
            energy_edges=[
                u.Quantity(self.args_dict["light_curve"]["energy"]["min"]).to("TeV").value,
                u.Quantity(self.args_dict["light_curve"]["energy"]["max"]).to("TeV").value,
            ]
            * u.TeV,
            time_intervals=time_intervals,
            reoptimize=False,
            selection_optional="all",
        )
        _light_curve = lc_maker_1d.run(datasets)

        display(_light_curve.to_table(sed_type="flux", format="lightcurve"))

        return _light_curve

    def _nightly_light_curve(self, _data_sets):
        """
        Combine observations per night and calculate light curve.
        Not applicable for stacked analysis.

        Parameters
        ----------
        _data_sets : Datasets
            Datasets

        """

        if self.args_dict["datasets"]["stack"]:
            return None

        self._logger.info("Estimating daily light curve")

        time_intervals = v2dl5_time.get_list_of_nights(
            _data_sets, time_zone=self.args_dict["light_curve"]["time_zone"]
        )

        print(time_intervals, type(time_intervals))

        return self._light_curve(_data_sets, time_intervals=time_intervals)

    def _get_energy_axis(self, name="energy"):
        """
        Get energy axis.

        """

        _axes_dict = self.args_dict["datasets"]["geom"]["axes"][name]

        return MapAxis.from_energy_bounds(
            u.Quantity(_axes_dict["min"]).to("TeV").value,
            u.Quantity(_axes_dict["max"]).to("TeV").value,
            nbin=_axes_dict["nbins"],
            per_decade=True,
            unit="TeV",
            name=name,
        )
