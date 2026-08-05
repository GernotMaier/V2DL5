# V2DL5 - High-level analysis for binary light curves for VERITAS

[![LICENSE](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://github.com/GernotMaier/V2DL5-Binary-Analysis/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/673002313.svg)](https://zenodo.org/badge/latestdoi/673002313)

[![gammapy](https://img.shields.io/badge/powered%20by-gammapy-orange.svg?style=flat)](https://www.gammapy.org/)

This is a collection of scripts for the high-level analysis of binary data for VERITAS.
Focus is on the usage of gammapy and anasum for binary light curve analysis, meaning this is mostly a reflected region analysis with some binary specific tools.

Allows to run analysis scripts for a given list of runs or for a cone search around the given on\_region direction.

- source detection analysis including integral flux (or flux upper limits), reflection region model
- spectral analysis, reflected region model
- run-list generation
- light-curve analysis and plotting

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

## Analysis

### Run list generator

The tool `generate_run_list.py` allows to generate a list of runs for a given observation time and zenith angle.
It generates a lot of printout which should be used to fine tune the run selection.

Example:

```console
python v2dl5/scripts/generate_run_list.py \
    --obs_table ../../../VTS/DL3/v490/dl3_pointlike_moderate2tel/obs-index.fits.gz \
    --config examples/run_selection.yml \
   --output_dir my_output_dir
```

### Reflected region analysis

Reflected region analysis using gammapy:

```console
python v2dl5/scripts/reflected_region_analysis.py \
    --run_list my_output_dir/runlist.txt \
    --config examples/reflected_region.yml \
    --output_dir my_output_dir
```

### Binary light curve plotting

Light-curve plotting is independent of the analysis code.
Input required is a table (ecsv format) with the light-curve data.

```console
python v2dl5/scripts/plot_binary_light_curves.py \
    --instrument VERITAS \
    --configuration examples/binary_lightcurve_plotting.yml
```

### Auxiliary data

Auxiliary data is stored in v2dl5/data and available at runtime. This includes:

#### Star catalogues

Hippargos catalog for stars with magnitude < 9: [v2dl5/data/hip_mag9.fits.gz](v2dl5/data/hip_mag9.fits.gz).
Star catalogs are listed in the configuration files as

```yaml
datasets:
    exclusion_region:
        on_radius: 0.5 deg
        magnitude_B: 7
        star_exclusion_radius: 0.30 deg
        fov: 3.5 deg
        star_file: hip_mag9.fits.gz
```

Star catalogs are expected to be in the [v2dl5/data](v2dl5/data) directory and of FITS format.

## External Eventdisplay and anasum analysis

High-level analysis is also possible with the [Eventdisplay](https://github.com/VERITAS-Observatory/EventDisplay_v4) package and its high-level tool anasum.

The following commands belong to the external Eventdisplay/anasum container workflow and are
not shipped by this Python package. The Python entry points listed above are the supported
commands in this repository.

### Light-curve analysis

```console
./light_curve_analysis.sh <anasum file> <time cuts file>
```

Use `RUNWISE` for the `time cuts file` to generate a light curve for each run.

### Spectral analysis

```console
./spectral_analysis.sh <anasum file> <output file>
```

### Time cuts file

The time cuts file is a simple text file with the following format:

```text
MJD_start MJD_stop
```

To generate a time cuts file for nightly binning:

```console
v2dl5-nightly-time-bins --anasum_file <anasum file> --output_file nights.ecsv
```

The generated ECSV file can be passed to the light-curve analysis as a time-bin
file; its columns are `time_min` and `time_max` in MJD.

### Validate time cuts and generate run lists

Check that every run from the anasum file is inside one of the defined time cuts:

```console
./runlist_from_time_bins.sh <anasum file> <time cuts file> false
```

Additionally write run lists for each of the MJD bins:

```console
./runlist_from_time_bins.sh <anasum file> <time cuts file> true
```

The tool will also write a figure with the length-of-time bins vs time bin.

## Acknowledgement

This work relies heavily on the [gammapy](https://gammapy.org/) development and especially on the excellent [tutorials](https://docs.gammapy.org/1.1/tutorials/index.html) provided by the gammapy team.
