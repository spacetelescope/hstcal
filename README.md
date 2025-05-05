# HSTCAL

Calibration software for HST/WFC3, HST/ACS, and HST/STIS.

HSTCAL is a C-based package which is comprised of the science calibration
software in support of the Advanced Camera for Surveys (ACS), Space Telescope
Imaging Spectrograph (STIS), and Wide Field Camera 3 pipelines. Initially, the
pipelines were written using C, but were encapsulated within the IRAF/STSDAS
environment, relying on IRAF to perform I/O and other basic interface
functions. HSTCAL replaces all the IRAF-based functionality with routines
based on the third-party package CFITSIO. This allows all the pipeline
software to be compiled and run without any dependence on IRAF. Not only
can HST data be processed using the C code directly via the C executables,
but the pipelines can also be run by using a high-level Python interface
(via [subprocess](https://docs.python.org/3/library/subprocess.html))
that serves as thin wrappers for the C executables.

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

## Instrument-specific Instructions

The following environment variables must be defined for the respective `cal[...].e` pipelines
to be able to find calibration reference files. For more information, see the
[CRDS user guide](https://hst-crds.stsci.edu/static/users_guide/index.html).

### ACS

`jref` must point to a path containing ACS reference files. For example:

```bash
export jref="/grp/crds/cache/references/hst/"
```

For more information, please see the [ACS Data Handbook](https://hst-docs.stsci.edu/acsdhb)
and the [acstools documentation](https://acstools.readthedocs.io/).

### WFC3

`iref` must point to a path containing WFC3 reference files. For example:

```bash
export iref="/grp/crds/cache/references/hst/"
```

For more information, please see the [WFC3 Data Handbook](https://hst-docs.stsci.edu/wfc3dhb)
and the [wfc3tools documentation](https://wfc3tools.readthedocs.io/).

### STIS

`oref` must point to a path containing STIS reference files. For example:

```bash
export oref="/grp/crds/cache/references/hst/"
```

For more information, please see the [STIS Data Handbook](https://hst-docs.stsci.edu/stisdhb)
and the [stistools documentation](https://stistools.readthedocs.io/).

## Other Software

The following are also part of HST pipeline but not in this package:

* COS pipeline: See [COS Data Handbook](https://www.stsci.edu/hst/cos/documents/handbooks/datahandbook/COS_cover.html) and [costools documentation](https://costools.readthedocs.io/)
* [Drizzlepac](https://www.stsci.edu/scientific-community/software/drizzlepac.html)

## Dev notes

[![CMake](https://github.com/spacetelescope/hstcal/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/spacetelescope/hstcal/actions/workflows/cmake.yml)
[![HST Calibration Pipeline](https://github.com/spacetelescope/RegressionTests/actions/workflows/hstcal.yml/badge.svg?branch=main)](https://github.com/spacetelescope/RegressionTests/actions/workflows/hstcal.yml)

Nightly regression test results are available from
[RegressionTests hstcal workflow](https://github.com/spacetelescope/RegressionTests/actions/workflows/hstcal.yml).

HSTCAL conda-forge recipe is hosted at [hstcal-feedstock](https://github.com/conda-forge/hstcal-feedstock/).

DMS deliveries are managed by [stasis](https://github.com/spacetelescope/stasis).

## Help Desk

If you need further assistance, please contact [HST Help Desk](https://hsthelp.stsci.edu).
