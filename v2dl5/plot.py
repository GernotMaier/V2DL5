"""
Plotting

"""

import logging
import warnings

import matplotlib.pyplot as plt
from astropy import units as u
from gammapy.datasets import FluxPointsDataset
from gammapy.makers.utils import make_theta_squared_table
from gammapy.maps import MapAxis
from gammapy.visualization import (
    plot_spectrum_datasets_off_regions,
    plot_theta_squared_table,
)

warnings.filterwarnings("ignore", category=UserWarning)


class Plot:
    """
    Plotting functions

    """

    def __init__(self, v2dl5_data, data_set, on_region=None, output_dir=None):
        self._logger = logging.getLogger(__name__)

        self.v2dl5_data = v2dl5_data
        self.data_set = data_set
        self.on_region = on_region
        self.output_dir = output_dir

    def default_offsets(self):
        """
        List of default offsets

        """

        _offsets = [0.5, 0.7, 1.0, 1.5] * u.deg
        self._logger.info(f"Default offsets for plotting: {_offsets}")
        return _offsets

    def default_energy_true(self):
        """
        List of default true energies

        """

        _energy_true = [0.2, 0.3, 1.0, 3.0, 10.0, 20.0] * u.TeV
        self._logger.info(f"Default true energies for plotting: {_energy_true}")
        return _energy_true

    def plot_maps(self, exclusion_mask=None):
        """
        Map and geometry related plots

        """

        self.plot_regions(exclusion_mask=exclusion_mask)
        self.plot_theta2()

    def plot_spectra(self, flux_points=None, model=None):
        """
        Spectrum related plots

        """

        for dataset in self.data_set:
            self.plot_fit(dataset)

        self.plot_flux_points(flux_points)
        self.plot_sed(
            FluxPointsDataset(data=flux_points, models=model.copy()),
        )

    def plot_light_curves(self, light_curve_per_obs=None, light_curve_per_night=None):
        """
        Light curve related plots

        """

        self.plot_light_curve(light_curve_per_obs, "per observation")
        self.plot_light_curve(light_curve_per_night, "per night")

    def plot_event_histograms(self):
        """
        Plot event histograms per observation

        """

        for obs in self.v2dl5_data.get_observations():
            obs.events.select_offset([0, 2.5] * u.deg).peek()
            try:
                self._plot(
                    plot_name=f"events_obs_{obs.obs_id}",
                    output_dir=self.output_dir / "events",
                )
            except TypeError:
                pass

    def plot_source_statistics(self):
        """
        Plot significance vs observation time

        """

        info_table = self.data_set.info_table(cumulative=True)
        _, (ax_excess, ax_sqrt_ts) = plt.subplots(figsize=(10, 4), ncols=2, nrows=1)
        ax_excess.plot(
            info_table["livetime"].to("h"),
            info_table["excess"],
            marker="o",
            ls="none",
        )

        ax_excess.set_title("Excess")
        ax_excess.set_xlabel("Livetime [h]")
        ax_excess.set_ylabel("Excess events")

        ax_sqrt_ts.plot(
            info_table["livetime"].to("h"),
            info_table["sqrt_ts"],
            marker="o",
            ls="none",
        )

        ax_sqrt_ts.set_title("Sqrt(TS)")
        ax_sqrt_ts.set_xlabel("Livetime [h]")
        ax_sqrt_ts.set_ylabel("Sqrt(TS)")
        try:
            self._plot(
                plot_name="source_statistics",
                output_dir=self.output_dir,
            )
        except TypeError:
            pass

    def plot_irfs(self):
        """
        Plot instrument response functions per observation

        """

        for obs in self.v2dl5_data.get_observations():
            self._plot_effective_area(obs)
            self._plot_energy_dispersion(obs)

    def plot_fit(self, data_set):
        """
        Plot successful fit results and residuals.

        """

        try:
            ax_spectrum, _ = data_set.plot_fit()
        except ValueError:
            return
        ax_spectrum.set_ylim(0.1, 40)
        data_set.plot_masks(ax=ax_spectrum)
        try:
            self._plot(
                plot_name=f"{data_set.name}_{data_set.models[0].name}_fit",
                output_dir=self.output_dir / "fit",
            )
        except TypeError:
            pass

    def plot_flux_points(self, flux_point_dataset):
        """
        Plot flux points

        """

        _, ax = plt.subplots()
        flux_point_dataset.plot(ax=ax, sed_type="dnde", color="darkorange")
        flux_point_dataset.plot_ts_profiles(ax=ax, sed_type="dnde")
        self._plot(plot_name="flux_points", output_dir=self.output_dir)

    def plot_sed(self, flux_point_dataset):
        """
        Plot spectral energy distribution

        """

        kwargs_model = {"color": "grey", "ls": "--", "sed_type": "dnde"}
        kwargs_fp = {"color": "black", "marker": "o", "sed_type": "dnde"}
        flux_point_dataset.plot_spectrum(kwargs_fp=kwargs_fp, kwargs_model=kwargs_model)
        self._plot(plot_name="spectrum", output_dir=self.output_dir)
        try:
            flux_point_dataset.plot_residuals(method="diff/model")
            self._plot(plot_name="residuals", output_dir=self.output_dir)
        except ValueError:
            pass

    def plot_light_curve(self, light_curve, plot_name):
        """
        Plot light curve

        """

        _, ax = plt.subplots(
            figsize=(8, 6),
            gridspec_kw={"left": 0.16, "bottom": 0.2, "top": 0.98, "right": 0.98},
        )

        try:
            light_curve.plot(ax=ax, marker="o", label=plot_name, sed_type="flux", time_format="mjd")
            ax.set_yscale("linear")
            self._plot(
                plot_name="light_curve_" + plot_name.replace(" ", "_"), output_dir=self.output_dir
            )
        except AttributeError:
            pass

    def plot_regions(self, exclusion_mask):
        """
        Plot on and off regions, exclusion mask

        """

        ax = exclusion_mask.plot()
        self.on_region.to_pixel(ax.wcs).plot(ax=ax, edgecolor="k")
        try:
            plot_spectrum_datasets_off_regions(ax=ax, datasets=self.data_set)
            self._plot(plot_name="regions", output_dir=self.output_dir)
        except AttributeError:
            pass

    def plot_theta2(self):
        """
        Plot theta2 distribution

        """

        theta2_axis = MapAxis.from_bounds(0, 0.2, nbin=20, interp="lin", unit="deg2")

        theta2_table = make_theta_squared_table(
            observations=self.v2dl5_data.get_observations(),
            position=self.on_region.center,
            theta_squared_axis=theta2_axis,
        )
        self._logger.info(f"Theta2 table {theta2_table}")

        try:
            plt.figure(figsize=(10, 5))
            plot_theta_squared_table(theta2_table)
            self._plot(plot_name="theta2", output_dir=self.output_dir)
        except AttributeError:
            pass

    def _plot(self, plot_name=None, output_dir=None):
        """
        Plotting helper function

        """
        if output_dir is not None:
            output_dir.mkdir(parents=True, exist_ok=True)
            _ofile = f"{output_dir}/{plot_name}.png"
            self._logger.info("Plotting %s", _ofile)
            plt.savefig(_ofile)
        else:
            plt.show()
        plt.close()

    def _plot_effective_area(self, obs):
        """
        Plot effective area

        """

        _, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
        obs.aeff.plot(ax=axes[2])
        obs.aeff.plot_energy_dependence(ax=axes[0], offset=self.default_offsets())
        obs.aeff.plot_offset_dependence(ax=axes[1], energy=self.default_energy_true())
        plt.tight_layout()

        try:
            self._plot(
                plot_name=f"aeff_obs_{obs.obs_id}",
                output_dir=self.output_dir / "irfs",
            )
        except TypeError:
            pass

    def _plot_energy_dispersion(self, obs):
        """
        Plot energy dispersion

        """

        _, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
        obs.edisp.plot_bias(
            ax=axes[0],
            offset=self.default_offsets()[0],
        )
        _mig_ax = obs.edisp.plot_migration(
            ax=axes[1],
            offset=self.default_offsets()[0],
            energy_true=self.default_energy_true(),
        )
        _mig_ax.legend(loc="upper right")
        _edisp = obs.edisp.to_edisp_kernel(offset=self.default_offsets()[0])
        _edisp.plot_matrix(ax=axes[2])

        plt.tight_layout()

        try:
            self._plot(
                plot_name=f"edisp_obs_{obs.obs_id}",
                output_dir=self.output_dir / "irfs",
            )
        except TypeError:
            pass
