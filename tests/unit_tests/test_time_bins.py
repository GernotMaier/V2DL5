"""Tests for time-bin producers and consumers."""

import numpy as np
from astropy.table import Table

from v2dl5.scripts.grouped_time_bins import find_groups, write_time_bins
from v2dl5.scripts.nightly_time_bins import get_unique_nights
from v2dl5.time import get_time_bins_from_file


def test_find_groups_honors_gap_and_group_size():
    data = Table(
        {
            "time_min": [0.0, 0.2, 0.4, 2.0],
            "time_max": [0.1, 0.3, 0.5, 2.1],
        }
    )

    assert find_groups(data, max_gap=0.5, max_group=2) == [[0], [1, 2], [3]]


def test_grouped_time_bins_round_trip_as_ecsv(tmp_path):
    data = Table(
        {
            "time_min": [0.25, 0.75],
            "time_max": [0.5, 1.0],
        }
    )
    output = tmp_path / "nested" / "time_bins.ecsv"

    write_time_bins(output, data, [[0, 1]])

    intervals = get_time_bins_from_file(output)
    assert intervals[0].mjd.tolist() == [0.25, 1.0]


def test_unique_nights_includes_intermediate_days():
    nights = get_unique_nights(np.array([10.9]), np.array([12.1]))

    assert nights.tolist() == [10, 11, 12]
