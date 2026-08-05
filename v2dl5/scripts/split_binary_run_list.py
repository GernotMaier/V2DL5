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
        help="Binary name (e.g., LS I +61 303; see v2dl5.binaries for definition).",
    )
    parser.add_argument(
        "--equal_phase_bin_count",
        type=int,
        required=False,
        default=None,
        help="Number of equal-width orbital-phase bins (default: 10).",
    )
    parser.add_argument(
        "--phase_bin_edges",
        type=float,
        nargs="+",
        default=None,
        help="Explicit orbital-phase bin edges, e.g. 0.0 0.2 0.3 ... 1.0.",
    )

    args = parser.parse_args()
    if args.equal_phase_bin_count is not None and args.phase_bin_edges is not None:
        parser.error("--equal_phase_bin_count and --phase_bin_edges are mutually exclusive")
    if args.equal_phase_bin_count is None and args.phase_bin_edges is None:
        args.equal_phase_bin_count = 10
    return args


def main():
    """Split a run list into run list per orbital phase bins."""
    logging.root.setLevel(logging.INFO)
    args = _parse()

    run_lists.split_binary_run_list(
        run_list_file=args.run_list,
        obs_table=args.obs_table,
        binary_name=args.binary_name,
        equal_phase_bin_count=args.equal_phase_bin_count,
        phase_bin_edges=args.phase_bin_edges,
    )


if __name__ == "__main__":
    main()
