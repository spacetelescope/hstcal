# HSTCAL

Calibration software for HST/WFC3, HST/ACS, and HST/STIS.

## Install using Conda (Complete Calibration Environment)

HSTCAL can be obtained as part of the
[Space Telescope Environment (stenv)](https://stenv.readthedocs.io/en/latest/), 
a conda environment.  The instructions found at this URL describe 
choosing the ``stenv`` release, building the ``stenv``
environment from the YAML file, and activating your new environment. Once your environment
is activated, you can invoke any of the HSTCAL executables.

```bash
# Use HSTCAL
$ calacs.e [...]
$ calwf3.e [...]
$ cs[...].e
```

## Install using Conda (HSTCAL only)

The [HSTCAL package is resident on conda-forge](https://anaconda.org/conda-forge/hstcal)
and you can use the following command to perform the installation:

```bash
$ conda install -c conda-forge hstcal==X.Y.Z
```
The `X.Y.Z` is the desired version number.

## Source Installation

See [detailed installation instructions](INSTALL.md).

## Dev notes

Nightly regression test results are available from
[RegressionTests hstcal workflow](https://github.com/spacetelescope/RegressionTests/actions/workflows/hstcal.yml).

HSTCAL conda-forge recipe is hosted at [hstcal-feedstock](https://github.com/conda-forge/hstcal-feedstock/).

DMS deliveries are managed by [stasis](https://github.com/spacetelescope/stasis).

## Help Desk

If you need further assistance, please contact [HST Help Desk](https://hsthelp.stsci.edu).
