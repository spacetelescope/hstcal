# HSTCAL

Calibration software for HST/WFC3, HST/ACS, and HST/STIS.

## Install using Conda (Complete Calibration Environment)

[![release](https://img.shields.io/github/v/release/spacetelescope/stenv)](https://github.com/spacetelescope/stenv/releases)

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

![Anaconda ver badge](https://anaconda.org/conda-forge/hstcal/badges/version.svg)
![Anaconda platform badge](https://anaconda.org/conda-forge/hstcal/badges/platforms.svg)

The [HSTCAL package is resident on conda-forge](https://anaconda.org/conda-forge/hstcal)
and you can use the following command to perform the installation:

```bash
$ conda install -c conda-forge hstcal==X.Y.Z
```
The `X.Y.Z` is the desired version number.

## Source Installation

See [detailed installation instructions](INSTALL.md).

## Environment Variables

The following environment variables must be defined for the respective `cal[...].e` pipelines
to be able to find calibration reference files. For more information, see the
[CRDS user guide](https://hst-crds.stsci.edu/static/users_guide/index.html).

### ACS

`jref` must point to a path containing ACS reference files. For example:

```bash
export jref="/grp/crds/cache/references/hst/"
```

For more information, please see the [ACS Data Handbook](https://hst-docs.stsci.edu/acsdhb). 

### WFC3

`iref` must point to a path containing WFC3 reference files. For example:

```bash
export iref="/grp/crds/cache/references/hst/"
```

For more information, please see the [WFC3 Data Handbook](https://hst-docs.stsci.edu/wfc3dhb).

### STIS

`oref` must point to a path containing STIS reference files. For example:

```bash
export oref="/grp/crds/cache/references/hst/"
```

For more information, please see the [STIS Data Handbook](https://hst-docs.stsci.edu/stisdhb). 

## Dev notes

[![CMake](https://github.com/spacetelescope/hstcal/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/spacetelescope/hstcal/actions/workflows/cmake.yml)
[![HST Calibration Pipeline](https://github.com/spacetelescope/RegressionTests/actions/workflows/hstcal.yml/badge.svg?branch=main)](https://github.com/spacetelescope/RegressionTests/actions/workflows/hstcal.yml)

Nightly regression test results are available from
[RegressionTests hstcal workflow](https://github.com/spacetelescope/RegressionTests/actions/workflows/hstcal.yml).

HSTCAL conda-forge recipe is hosted at [hstcal-feedstock](https://github.com/conda-forge/hstcal-feedstock/).

DMS deliveries are managed by [stasis](https://github.com/spacetelescope/stasis).

## Help Desk

If you need further assistance, please contact [HST Help Desk](https://hsthelp.stsci.edu).
