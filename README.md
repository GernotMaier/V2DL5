# V2DL5 - high-level analysis for VERITAS with gammapy

[![DOI](https://zenodo.org/badge/673002313.svg)](https://zenodo.org/badge/latestdoi/673002313)
[![gammapy](https://img.shields.io/badge/powered%20by-gammapy-orange.svg?style=flat)](https://www.gammapy.org/)
[![GitHub Super-Linter](https://github.com/GernotMaier/V2DL5/actions/workflows/linter.yml/badge.svg)](https://github.com/marketplace/actions/super-linter)

This is a collection of simple scripts for the high-level analysis of VERITAS data with gammapy.
The focus is on binary light curve analysis, meaning this is mostly a reflected region analysis with some binary specific tools.

Allows to run analysis scripts for a given list of runs or for a cone search around the given on\_region direction.

- source detection analysis including integral flux (or flux upper limits), reflection region model
- spectral analysis, reflected region model

Reminder on data levels (as defined e.g., by CTAO):

- DL3: event lists for gamma-ray like events
- DL4: binned data (in time, space, energy)
- DL5: science data products (e.g., sky maps, energy spectra, light curves)

This repository provides scripts to generate DL5 science data products.

## Acknowledgement

This work relies heavily on the [gammapy](https://gammapy.org/) development and especially on the excellent [tutorials](https://docs.gammapy.org/1.1/tutorials/index.html) provided by the gammapy team.

## Installation

To install the required python packages, run:

```bash
conda env create -f environment.yml
```

Activate the environment to start the analysis:

```bash
conda activate v2dl5
pip install -e .
```

## Run list generator

The tool `runlist.py` allows to generate a list of runs for a given observation time and zenith angle.
It generates a lot of printout which should be used to fine tune the run selection.

Example:

```console
python v2dl5/scripts/generate_runlist.py \
    --obs_table ../../../VTS/DL3/v490/dl3_pointlike_moderate2tel/obs-index.fits.gz \
    --config examples/run_selection.yml \
   --output_dir my_output_dir
```

## Reflected region analysis

```console
python v2dl5/scripts/reflected_region_analysis.py \
    --obs_table ../../../VTS/DL3/v490/dl3_pointlike_moderate2tel/obs-index.fits.gz \
    --runlist my_output_dir/runlist.txt \
    --config examples/reflected_region_analysis.yml \
    --output_dir my_output_dir
```

## Binary light curve plotting

```console
python v2dl5/scripts/plot_binary_light_curves.py \
    --instrument VERITAS \
    --configuration examples/binary_lightcurve_plotting.yml
```
