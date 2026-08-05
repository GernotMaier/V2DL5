#!/usr/bin/python

"""
Split a run list into run lists per orbital phase bins.

Used as input for gammapy or anasum analysis.

"""

import argparse
import logging

from v2dl5 import run_lists


def _parse():
    """
    Parse command line arguments.

    Returns
    -------
    dict
        Command line parameters.

    """
    parser = argparse.ArgumentParser(
        description="Split a run list into run list per orbital phase bins."
    )

    parser.add_argument(
        "--run_list",
        type=str,
        required=True,
        help="Path to the run list.",
    )
    parser.add_argument(
        "--obs_table",
        type=str,
        required=True,
        help="Path to observation table.",
    )
    parser.add_argument(
        "--binary_name",
        type=str,
        required=True,
        default=None,
        help="Binary name (e.g., LS I +61 303; see v2dl5.binaries for definition).",
    )
    parser.add_argument(
        "--orbital_bins",
        type=int,
        required=False,
        default=10,
        help="Number of bins in orbital period for averaging.",
    )

    return parser.parse_args()


def main():
    """Split a run list into run list per orbital phase bins."""
    logging.root.setLevel(logging.INFO)
    args = _parse()

    run_lists.split_binary_run_list(
        run_list_file=args.run_list,
        obs_table=args.obs_table,
        binary_name=args.binary_name,
        orbital_bins=args.orbital_bins,
    )


if __name__ == "__main__":
    main()
