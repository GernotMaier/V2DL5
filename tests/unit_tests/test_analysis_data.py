"""Tests for analysis and data-selection contracts."""

from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

import v2dl5.analysis as analysis_module
import v2dl5.data as data_module
from v2dl5.analysis import Analysis
from v2dl5.data import Data


def test_flux_point_edges_match_requested_bin_count(monkeypatch):
    captured = {}

    class FakeTable:
        def pprint(self):
            pass

        def pprint_all(self):
            pass

    class FakeFluxPoints:
        def to_table(self, **_kwargs):
            return FakeTable()

    class FakeEstimator:
        def __init__(self, energy_edges, **_kwargs):
            captured["energy_edges"] = energy_edges

        def run(self, datasets):
            assert datasets is not None
            return FakeFluxPoints()

    monkeypatch.setattr(analysis_module, "FluxPointsEstimator", FakeEstimator)
    analysis = Analysis(
        args_dict={
            "output_dir": ".",
            "flux_points": {
                "energy": {"min": "100 GeV", "max": "10 TeV", "nbins": 3}
            },
        }
    )

    analysis._flux_points(object())

    assert len(captured["energy_edges"]) - 1 == 3
    assert captured["energy_edges"].unit == u.TeV


def test_default_star_catalogue_is_packaged():
    from importlib.resources import files

    catalogue = files("v2dl5").joinpath("data", "hip_mag9.fits.gz")
    with fits.open(catalogue) as hdul:
        assert len(hdul) > 1


def test_bti_is_applied_to_each_observation_cache(monkeypatch):
    observation = SimpleNamespace(obs_id=7, gti="original")

    class Store:
        obs_table = Table({"OBS_ID": [7]})

        def get_observations(self, *_args, **_kwargs):
            return [observation]

    class FakeBTI:
        def __init__(self, obs):
            assert obs is observation

        def update_gti(self, pairs):
            assert pairs == [(1.0, 3.0)]
            return "updated"

    data = Data.__new__(Data)
    data._logger = data_module.logging.getLogger(__name__)
    data._observation_cache = {}
    data._data_store = Store()
    data.runs = [7]
    data._bti = [{"run": 7, "bti_start": 1.0, "bti_length": 2.0}]
    monkeypatch.setattr(data_module.BTI, "BTI", FakeBTI)

    assert data.get_observations(reflected_region=True)[0].gti == "updated"
    assert data.get_observations(reflected_region=False)[0].gti == "updated"


def test_target_selection_accepts_quantity_strings_and_rejects_empty():
    class ObservationTable:
        pointing_radec = SkyCoord([83.6], [22.0], unit="deg")

        def __getitem__(self, mask):
            return Table({"OBS_ID": np.array([1])[mask]})

    data = Data.__new__(Data)
    data._logger = data_module.logging.getLogger(__name__)
    data.target = SkyCoord(83.6, 22.0, unit="deg")
    data._data_store = SimpleNamespace(obs_table=ObservationTable())

    assert data._from_target("1 deg").tolist() == [1]
    data.target = SkyCoord(90.0, 22.0, unit="deg")
    with pytest.raises(ValueError, match="No observations"):
        data._from_target("0.001 deg")
