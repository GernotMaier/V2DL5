#!/usr/bin/python3
"""
Perform reflected region analysis for a given run list.

Example:

.. code-block:: console

    $ python reflected_region_analysis.py \
        --run_list region_run_list.txt \
        --output_dir 1D_reflected_region_output \
        --config examples/1D_reflected_region_config.yaml

"""

import argparse
import logging

import v2dl5.analysis
import v2dl5.configuration
import v2dl5.data
import v2dl5.sky_regions


def _parse():
    """
    Parse command line arguments.

    """
    parser = argparse.ArgumentParser(description="Perform reflected region analysis.")

    parser.add_argument(
        "--run_list",
        type=str,
        required=False,
        help="List of run IDs to be processed.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to configuration file.",
    )

    # Get command line arguments
    return parser.parse_args()


def main():
    """
    Reflected region analysis main function.

    """
    logging.basicConfig(level=logging.INFO)

    args_dict = v2dl5.configuration.configuration(args=_parse())

    sky_regions = v2dl5.sky_regions.SkyRegions(args_dict=args_dict)
    data = v2dl5.data.Data(args_dict=args_dict, target=sky_regions.target)
    sky_regions.update_regions(
        args_dict=args_dict,
        on_region_radius=data.get_on_region_radius(),
        max_wobble_distance=data.get_max_wobble_distance(),
    )

    analysis = v2dl5.analysis.Analysis(
        args_dict=args_dict,
        sky_regions=sky_regions,
        data=data,
    )
    analysis.run()
    analysis.plot()
    analysis.write()


if __name__ == "__main__":
    main()