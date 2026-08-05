"""Tests for light-curve input and plotting contracts."""

import astropy.units as u
import numpy as np
from astropy.table import Table

from v2dl5.light_curves.binary_plotting import BinaryLightCurvePlotter
from v2dl5.light_curves.data_reader import LightCurveDataReader


def test_reader_reads_ecsv_and_creates_phase_errors(tmp_path):
    table_file = tmp_path / "light_curve.ecsv"
    Table(
        {
            "time_min": [55000.0],
            "time_max": [55001.0],
            "flux": [1.0],
            "flux_err": [0.1],
        }
    ).write(table_file, format="ascii.ecsv")
    config_file = tmp_path / "config.yml"
    config_file.write_text(
        "data:\n  - instrument: test\n    file_name: "
        f"{table_file}\n",
        encoding="utf-8",
    )

    reader = LightCurveDataReader(
        config_file,
        binary={"orbital_period": 10.0, "mjd_0": 55000.0},
    )
    reader.read_data()

    data = reader.data_dict["test"]
    assert "phase_err_high" in data
    assert len(data["phase_err"]) == 1
    assert len(data["phase_err"][0]) == 2


def test_photon_to_energy_flux_uses_input_values():
    reader = object.__new__(LightCurveDataReader)
    values, errors = reader.convert_photon_to_energy_flux(
        [1.0, 2.0], 1.0 * u.TeV, 2.5, [0.1, 0.2]
    )

    factor = 3.0 * 1.602176634
    assert np.allclose(values, np.array([1.0, 2.0]) * factor)
    assert np.allclose(errors, np.array([0.1, 0.2]) * factor)


def test_plotter_uses_instrument_names_and_handles_missing_significance():
    data = {
        "second": {
            "MJD": [1.0],
            "phase": [0.1],
            "orbit_number": [1],
            "flux": [2.0],
            "flux_err": [0.1],
        },
        "first": {
            "MJD": [1.0],
            "phase": [0.1],
            "orbit_number": [1],
            "flux": [1.0],
            "flux_err": [0.1],
        },
    }
    config = [
        {"instrument": "first", "plot_variable": "flux"},
        {"instrument": "second", "plot_axis": ["flux"]},
    ]
    plotter = BinaryLightCurvePlotter(data, config, {"name": "test"})

    x, y, error, _, _ = plotter._get_light_curve_in_mjd_limits(
        data["second"], "flux", "MJD", min_significance=5
    )

    assert (x, y, error) == ([1.0], [2.0], [0.1])
    assert plotter.plot_this_instrument(config[0], "flux")
