import logging

import astropy.units as u
import numpy as np
import yaml
from astropy.table import Table

import v2dl5.light_curves.orbital_period as orbit


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
            with open(configuration_file, "r") as file:
                _yml_in = yaml.safe_load(file)
                self.config = _yml_in["data"]
        except FileNotFoundError as err:
            self._logger.error("Configuration file not found: %s", configuration_file)
            raise err
        except KeyError as err:
            self._logger.error("Data key not found in configuration file %s", configuration_file)
            raise err

        self.binary = binary

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
            self._add_orbital_parameters(self.data_dict[data_config["instrument"]])

    def _read_fluxes_from_file(self, data_config):
        """
        Read flux from file.

        """

        try:
            if data_config["file_name"].endswith(".csv") or data_config["file_name"].endswith(
                ".ecsv"
            ):
                return self._read_fluxes_from_ecsv_file(data_config)
        except KeyError:
            self._logger.error("File name not found in configuration file")
            raise KeyError

    def _add_orbital_parameters(self, data):
        """
        Add orbital phase to light-curve data data.

        Parameters
        ----------
        data: dict
            Data dictionary with light-curve data.

        """

        orbital_period = self.binary["orbital_period"]
        mjd_0 = self.binary["mjd_0"]

        data["phase"] = [
            orbit.get_orbital_phase(mjd, orbital_period, mjd_0, True) for mjd in data["MJD"]
        ]
        data["phaseN"] = [
            orbit.get_orbital_phase(mjd, orbital_period, mjd_0, False) for mjd in data["MJD"]
        ]
        data["phase_err_low"] = [
            orbit.get_orbital_phase_range(a, b, c, False, orbital_period, mjd_0)
            for a, b, c in zip(data["time_min"], data["time_max"], data["phase"])
        ]
        data["phase_err_hig"] = [
            orbit.get_orbital_phase_range(a, b, c, True, orbital_period, mjd_0)
            for a, b, c in zip(data["time_min"], data["time_max"], data["phase"])
        ]
        data["phase_err"] = [data["phase_err_low"], data["phase_err_hig"]]

    def convert_photon_to_energy_flux(C_v, C_e, E_0, gamma):
        """
        Convert photon to energy flux

        """

        f = (-1.0 * gamma + 1) / (-1.0 * gamma + 2)
        # conversion to erg
        f = f * (E_0.to(u.erg)).value
        return [v * f for v in C_v], [e * f for e in C_e]

    def _read_fluxes_from_ecsv_file(self, data_config, TimeMinMax=True, MJD_min=-1.0, MJD_max=-1.0):
        """
        Read gamma-ray fluxes from ecsv file (open gamma-ray format)

        Parameters:
        -----------
        data_config: dict
            configuration dictionary
        TimeMinMax: bool
            flag to read time_min and time_max
        MJD_min: float
            MJD min value for MJD cut
        MJD_max: float
            MJD max value for MJD cut

        """
        table = Table.read(data_config["file_name"])
        f = {}
        if not TimeMinMax:
            table["time_min"] = table["time"].data
            table["time_max"] = table["time"].data + 0.1

        # MJD filter
        condition = np.ones(len(table), dtype=bool)
        if MJD_min > -1:
            condition &= table["time_min"] > MJD_min
        if MJD_max > -1:
            condition &= table["time_max"] < MJD_max
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
            f["flux_ul"] = [-1.0 for _ in f["flux"]]
        return f
