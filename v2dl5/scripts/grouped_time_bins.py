#!/usr/bin/python

"""
Light curve binning based on nightly observations.

Groups observations in order to obtainer longer bins based on the following rules:

- Never group observations with gaps > 1 day
- Prefer groups of 3 observations when possible
- Avoid single observations at end of sequences
- Split longer sequences into optimal groups

Example usage:

    python group.py
    python group.py --max_gap 2.0 --max_group 4
"""

import argparse

import numpy as np
from astropy.table import Table


def should_group(times, max_gap=1.0, max_group=2):
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


def find_groups(data):
    """
    Find groups of observations based on the specified rules.

    In detail, the rules are:

    - Never group observations with gaps > 1 day
    - Prefer groups of 3 observations when possible
    - Avoid single observations at end of sequences
    - Split longer sequences into optimal groups
    """
    groups = []
    used_indices = set()
    i = 0

    while i < len(data):
        if i in used_indices:
            i += 1
            continue

        # Look ahead for consecutive observations
        sequence = [i]
        j = i + 1
        while j < len(data) and j not in used_indices:
            gap = data['time_min'][j] - data['time_max'][sequence[-1]]
            if gap <= 1.0:
                sequence.append(j)
            else:
                break
            j += 1

        # Handle different sequence lengths
        if len(sequence) >= 7:
            # Split into triplets plus remainder
            num_groups = len(sequence) // 3
            remainder = len(sequence) % 3

            if remainder == 1:  # Avoid single observation at end
                num_groups -= 1
                remainder = 4  # Make last group a quartet

            for k in range(num_groups):
                groups.append(sequence[k*3:(k+1)*3])
            if remainder:
                groups.append(sequence[num_groups*3:])

            used_indices.update(sequence)
            i = sequence[-1] + 1

        elif len(sequence) == 6:
            # Split into two triplets
            groups.append(sequence[:3])
            groups.append(sequence[3:])
            used_indices.update(sequence)
            i = sequence[-1] + 1

        elif len(sequence) == 5:
            # Split 3+2 or 2+3 based on context
            if i + 5 < len(data) and data['time_min'][i+5] - data['time_max'][sequence[-1]] <= 1.0:
                # If next observation is close, do 2+3 to avoid single observation
                groups.append(sequence[:2])
                groups.append(sequence[2:])
            else:
                # Otherwise do 3+2
                groups.append(sequence[:3])
                groups.append(sequence[3:])
            used_indices.update(sequence)
            i = sequence[-1] + 1

        elif len(sequence) == 4:
            # Keep as quartet to avoid single observation
            groups.append(sequence)
            used_indices.update(sequence)
            i = sequence[-1] + 1

        elif len(sequence) == 3 or len(sequence) == 2:
            groups.append(sequence)
            used_indices.update(sequence)
            i = sequence[-1] + 1

        else:
            groups.append([i])
            used_indices.add(i)
            i += 1

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
        help="Output file for time bins (ASCII format)"
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
    """Write time bins to ASCII file."""
    with open(filename, 'w') as f:
        for group in groups:
            start_time = data['time_min'][group[0]]
            end_time = data['time_max'][group[-1]]
            f.write(f"{int(start_time)} {int(np.ceil(end_time))}\n")


def main():
    """Group light curve observations into broader time bins."""
    args = parse_args()

    data = Table.read(filename=args.light_curve_file, format='ascii.ecsv')

    groups = find_groups(data)

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
