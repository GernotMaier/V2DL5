#!/usr/bin/python3
"""
Generate a run list for a given sky region or a given target.

Reads in a configuration yaml file (see example `examples/run_selection.yaml`)
with criteria for selecting runs.

Generates three output files:

1. A list of runs that satisfy the criteria.
2. A detailed printout of the observation table for all runs which satisfy the criteria.
3. A detailed printout of the observation table for all runs which do not satisfy the criteria.

Applies a simple outlier detection based on mean, median, standard, and absolute deviation of L3Rate
and FIR1. This should be used as indications for further investigation, not as hard cuts to remove
runs.

Example:

.. code-block:: console

        $ python generate_run_list.py \
            --obs_table /path/to/ob_table.fits \
            --config examples/run_selection.yaml \
            --output_dir run_lists

"""

import argparse
import logging
from pathlib import Path

import v2dl5.configuration
import v2dl5.run_lists
import v2dl5.sky_regions

_logger = logging.getLogger(__name__)


def _parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate a run list.")

    parser.add_argument(
        "--obs_table",
        type=str,
        required=True,
        help="Path to observation table.",
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to configuration file.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory.",
    )

    return parser.parse_args()


def main():
    """Generate a run list."""
    logging.root.setLevel(logging.INFO)
    args_dict = v2dl5.configuration.configuration(args=_parse(), generate_dqm_run_list=True)
    output_path = Path(args_dict["output_dir"])
    output_path.mkdir(parents=True, exist_ok=True)

    sky_regions = v2dl5.sky_regions.SkyRegions(args_dict=args_dict)
    v2dl5.run_lists.generate_run_list(args_dict=args_dict, target=sky_regions.target)


if __name__ == "__main__":
    main()
