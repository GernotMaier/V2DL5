"""Tests for run-list utilities."""

from astropy.table import Table

from v2dl5.run_lists import split_binary_run_list


def test_split_binary_run_list_matches_numeric_ids(tmp_path):
    run_list = tmp_path / "runs.txt"
    run_list.write_text("1\n2\n", encoding="utf-8")
    observation_table = tmp_path / "observations.ecsv"
    Table(
        {
            "OBS_ID": [1, 2],
            "DATE-OBS": ["2020-01-01T00:00:00", "2020-01-02T00:00:00"],
            "LIVETIME": [10.0, 20.0],
        }
    ).write(observation_table, format="ascii.ecsv")

    split_binary_run_list(run_list, observation_table, "LS I +61 303", 2)

    output = tmp_path / "runs_orbital_bin_00.txt"
    assert output.exists()
