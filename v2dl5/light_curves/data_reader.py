"""Data reader for light-curve analysis."""

import logging
import os

import astropy.units as u
import numpy as np
import yaml
from astropy.table import Table

import v2dl5.orbital_phase as orbit
from v2dl5 import binaries


class LightCurveDataReader:
    """
    Data reader for light-curve analysis.

    Allows to read gamma-ray, optical or x-ray data. Enrich data
    with information required for orbital analysis.

    Parameters
    ----------
    configuration_file: str
        Configuration file with data description.
    binary: dict
        Binary parameters.

    """

    def __init__(self, configuration_file, binary=None):

        self._logger = logging.getLogger(__name__)

        self.data_dict = {}
        try:
            with open(configuration_file) as file:
                _yml_in = yaml.safe_load(file)
                self.config = _yml_in["data"]
        except FileNotFoundError as err:
            self._logger.error("Configuration file not found: %s", configuration_file)
            raise err
        except KeyError as err:
            self._logger.error("Data key not found in configuration file %s", configuration_file)
            raise err

        self.binary = binary
        self._logger.info("Binary properties: %s", self.binary)

    def read_data(self):
        """
        Read data from files as described in configuration file.

        Returns
        -------
        dict
            Data dictionary with data for each instrument.

        """
        for data_config in self.config:
            self.data_dict[data_config["instrument"]] = self._read_fluxes_from_file(data_config)
            self._add_orbital_parameters(self.data_dict[data_config["instrument"]], data_config)

    def _read_fluxes_from_file(self, data_config):
        """Read flux from file."""
        try:
            data_config["file_name"] = os.path.expandvars(data_config["file_name"])
        except KeyError:
            self._logger.error(f"File name not found in configuration {data_config}")
            raise KeyError

        if data_config["file_name"].endswith((".csv", ".ecsv")):
            self._logger.info("Reading data from %s", data_config["file_name"])
            return self._read_fluxes_from_ecsv_file(
                file_name=data_config["file_name"],
                mjd_min=data_config.get("mjd_min", -1.0),
                mjd_max=data_config.get("mjd_max", -1.0),
            )

        return None

    def _apply_phase_mask(self, data, data_config):
        """
        Apply phase mask to data.

        Requires 'phase_cut_binary' to be set in data_config.
        """
        if not data_config.get("phase_cut_binary"):
            return data

        try:
            _binary = binaries.binary_properties()[data_config["phase_cut_binary"]]
        except KeyError:
            raise KeyError(f"Binary {data_config['phase_cut_binary']} not found in binaries.py")

        phase_min = data_config.get("phase_min", 0.0)
        phase_max = data_config.get("phase_max", 1.0)

        phases = [
            orbit.get_orbital_phase(
                mjd,
                orbital_period=_binary["orbital_period"],
                mjd_0=_binary["mjd_0"],
                phase_reduce=True
            )
            for mjd in data["MJD"]
        ]
        phase_mask = [
            (p >= phase_min) & (p <= phase_max) for p in phases
        ]
        for key in data.keys():
            data[key] = [val for val, mask in zip(data[key], phase_mask) if mask]

        return data

    def _add_orbital_parameters(self, data, data_config):
        """
        Add orbital phase to light-curve data and filter by phase range.

        Parameters
        ----------
        data: dict
            Data dictionary with light-curve data.
        data_config: dict
            Data configuration dictionary.
        """
        orbital_period = self.binary["orbital_period"]
        mjd_0 = self.binary["mjd_0"]

        self._apply_phase_mask(data, data_config)

        data["phase"] = [
            orbit.get_orbital_phase(
                mjd, orbital_period=orbital_period, mjd_0=mjd_0, phase_reduce=True
            )
            for mjd in data["MJD"]
        ]
        data["phaseN"] = [
            orbit.get_orbital_phase(
                mjd, orbital_period=orbital_period, mjd_0=mjd_0, phase_reduce=False
            )
            for mjd in data["MJD"]
        ]
        data["phase_err_low"] = [
            orbit.get_orbital_phase_range(
                a, b, c, upper_error=False, orbital_period=orbital_period, mjd_0=mjd_0
            )
            for a, b, c in zip(data["time_min"], data["time_max"], data["phase"])
        ]
        data["phase_err_hig"] = [
            orbit.get_orbital_phase_range(
                a, b, c, upper_error=True, orbital_period=orbital_period, mjd_0=mjd_0
            )
            for a, b, c in zip(data["time_min"], data["time_max"], data["phase"])
        ]
        data["phase_err"] = [data["phase_err_low"], data["phase_err_hig"]]

        data["orbit_number"] = [
            orbit.get_orbit_number(
                mjd=mjd,
                orbital_period=orbital_period,
                mjd_orbit_count=float(self.binary.get("mjd_orbit_count", self.binary.get("mjd_0")))
            )
            for mjd in data["MJD"]
        ]

    def convert_photon_to_energy_flux(self, c_e, e_0, gamma):
        """Convert photon to energy flux."""
        f = (-1.0 * gamma + 1) / (-1.0 * gamma + 2)
        # conversion to erg
        f = f * (e_0.to(u.erg)).value
        return [v * f for v in self], [e * f for e in c_e]

    def _read_fluxes_from_ecsv_file(
        self, file_name, time_min_max=True, mjd_min=-1.0, mjd_max=-1.0
    ):
        """
        Read gamma-ray fluxes from ecsv file (open gamma-ray format).

        Parameters
        ----------
        file_name: str, Path
            Name of the file to read.
        time_min_max: bool
            flag to read time_min and time_max
        mjd_min: float
            MJD min value for MJD cut
        mjd_max: float
            MJD max value for MJD cut

        """
        table = Table.read(file_name)
        f = {}

        if "time" not in table.colnames and "MJD" in table.colnames:
            table.rename_column("MJD", "time")
            time_min_max = False

        if not time_min_max:
            table["time_min"] = table["time"].data
            table["time_max"] = table["time"].data + 0.1

        # MJD filter
        condition = np.ones(len(table), dtype=bool)
        if mjd_min > -1:
            condition &= table["time_min"] > mjd_min
        if mjd_max > -1:
            condition &= table["time_max"] < mjd_max
        table = table[condition]

        f = {}
        f["time_min"] = table["time_min"].data.tolist()
        f["time_max"] = table["time_max"].data.tolist()
        f["flux"] = table["flux"].data.flatten().tolist()
        if "flux_err" in table.colnames:
            f["flux_err"] = table["flux_err"].data.flatten().tolist()
        else:
            up = table["flux_up"].data.flatten().tolist()
            down = table["flux_down"].data.flatten().tolist()
            f["flux_err"] = [0.5 * abs(u - d) for u, d in zip(up, down)]
        f["MJD"] = [0.5 * (a + b) for a, b in zip(f["time_min"], f["time_max"])]
        f["MJD_err"] = [0.5 * (b - a) for a, b in zip(f["time_min"], f["time_max"])]
        if "flux_ul" in table.colnames:
            flux_ul = table["flux_ul"].data.flatten().tolist()
            is_ul = table["is_ul"].data.flatten().tolist()
            f["flux_ul"] = [flux if is_ul else -1.0 for flux, is_ul in zip(flux_ul, is_ul)]
        else:
            f["flux_ul"] = [-1.0 if fe > 0 else fl for fe, fl in zip(f["flux_err"], f["flux"])]
        for column_name in ["significance", "index", "index_err", "live_time"]:
            if column_name in table.colnames:
                f[column_name] = table[column_name].data.flatten().tolist()

        return f
