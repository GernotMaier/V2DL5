#!/usr/bin/python3
"""
Generate run lists for list of targets.

Reads in a configuration yaml file (see example `examples/run_selection.yaml`)
with criteria for selecting runs (ignore 'on_region:target' in this example file).

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

from astropy.table import Table

import v2dl5.configuration
import v2dl5.run_lists
import v2dl5.sky_regions

logger = logging.getLogger(__name__)


def _parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate a run list.")

    parser.add_argument(
        "--target_list",
        type=str,
        required=True,
        help="List of targets to be processed.",
    )
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


def read_list_of_targets(target_list_file):
    """Read list of targets from file."""
    if Path(target_list_file).suffix == ".ecsv":
        object_table = Table.read(target_list_file, format="ascii.ecsv")
        target_list = [
            typed_id if typed_id and len(str(typed_id).strip()) > 0 else main_id
            for typed_id, main_id in zip(object_table["TYPED_ID"], object_table["MAIN_ID"])
        ]
    else:
        with open(target_list_file, encoding="utf-8") as file:
            target_list = [line.strip() for line in file]

    logger.info(f"Found {len(target_list)} targets in {target_list_file}")
    return target_list


def main():
    """Generate run lists for a list of targets."""
    logging.root.setLevel(logging.INFO)

    args_dict = v2dl5.configuration.configuration(args=_parse(), generate_dqm_run_list=True)

    logger.info("args_dict: %s", args_dict)

    target_list = read_list_of_targets(args_dict["target_list"])

    output_path_main = Path(args_dict["output_dir"])

    for target in target_list:
        target_dir_name = (
            target.replace(" ", "_")
            .replace(",", "")
            .replace("+", "p")
            .replace("-", "m")
            .replace("[", "")
            .replace("]", "")
        )

        args_dict["output_dir"] = output_path_main / target_dir_name

        output_path = Path(args_dict["output_dir"])
        output_path.mkdir(parents=True, exist_ok=True)

        args_dict["on_region"]["target"] = target

        sky_regions = v2dl5.sky_regions.SkyRegions(args_dict=args_dict)
        v2dl5.run_lists.generate_run_list(args_dict=args_dict, target=sky_regions.target)


if __name__ == "__main__":
    main()
