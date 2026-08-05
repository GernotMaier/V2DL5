#!/usr/bin/python

"""
Light curve binning based on nightly observations.

Groups observations in order to obtainer longer bins based on the following rules:

- Never group observations with gaps > 1 day
- Prefer groups of 3 observations when possible
- Avoid single observations at end of sequences
- Split longer sequences into optimal groups

Example usage:

    v2dl5-grouped-time-bins --light_curve_file light_curve.ecsv \
        --time_bins_file time_bins.ecsv
    v2dl5-grouped-time-bins --light_curve_file light_curve.ecsv \
        --time_bins_file time_bins.ecsv --max_gap 2.0 --max_group 4
"""

import argparse
from pathlib import Path

import numpy as np
from astropy.table import Table


def should_group(times, max_gap=1.0, max_group=4):
    """
    Determine if a sequence of times should be grouped together.

    Parameters
    ----------
    times : array-like
        Array of time values
    max_gap : float
        Maximum allowed gap between observations in days
    max_group : int
        Maximum allowed group size (default 2 to prefer pairs)
    """
    if len(times) < 2:
        return True

    gaps = np.diff(times)
    return all(gap <= max_gap for gap in gaps) and len(times) <= max_group


def _split_sequence(sequence, max_group):
    """Split a contiguous sequence without leaving a final singleton."""
    groups = []
    sequence = list(sequence)
    while len(sequence) > max_group:
        take = max_group
        if len(sequence) - take == 1:
            take -= 1
        groups.append(sequence[:take])
        sequence = sequence[take:]
    if sequence:
        groups.append(sequence)
    return groups


def find_groups(data, max_gap=1.0, max_group=4):
    """
    Find groups of observations based on the specified rules.

    In detail, the rules are:

    - Never group observations with gaps > 1 day
    - Prefer groups of 3 observations when possible
    - Avoid single observations at end of sequences
    - Split longer sequences into optimal groups
    """
    if max_group < 2:
        raise ValueError("max_group must be at least 2")

    groups = []
    sequence = []
    for index in range(len(data)):
        if sequence:
            gap = data["time_min"][index] - data["time_max"][sequence[-1]]
            if gap > max_gap:
                groups.extend(_split_sequence(sequence, max_group))
                sequence = []
        sequence.append(index)
    groups.extend(_split_sequence(sequence, max_group))
    return groups


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Group light curve observations into broader time bins."
    )
    parser.add_argument(
        "--light_curve_file",
        type=str,
        required=True,
        help="Input light curve file in ECSV format"
    )
    parser.add_argument(
        "--time_bins_file",
        type=str,
        required=True,
        help="Output file for time bins (ECSV format)"
    )
    parser.add_argument(
        "--max_gap",
        type=float,
        default=1.0,
        help="Maximum gap between observations in days (default: 1.0)"
    )
    parser.add_argument(
        "--max_group",
        type=int,
        default=4,
        help="Maximum group size (default: 4)"
    )
    return parser.parse_args()


def write_time_bins(filename, data, groups):
    """Write time bins as an ECSV table."""
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    output = Table(
        {
            "time_min": [data["time_min"][group[0]] for group in groups],
            "time_max": [data["time_max"][group[-1]] for group in groups],
        }
    )
    output.write(filename, format="ascii.ecsv", overwrite=True)


def main():
    """Group light curve observations into broader time bins."""
    args = parse_args()

    data = Table.read(filename=args.light_curve_file, format="ascii.ecsv")

    groups = find_groups(data, max_gap=args.max_gap, max_group=args.max_group)

    write_time_bins(args.time_bins_file, data, groups)

    print(f"Found {len(groups)} groups")
    for i, group in enumerate(groups):
        times = data['time_min'][group]
        print(f"\nGroup {i+1}:")
        print(f"Start MJD: {times[0]:.3f}")
        print(f"End MJD: {data['time_max'][group[-1]]:.3f}")
        print(f"Number of observations: {len(group)}")
        print(f"Timestamps: {', '.join(f'{t:.3f}' for t in times)}")

    print(f"\nTime bins written to: {args.time_bins_file}")


if __name__ == "__main__":
    main()
