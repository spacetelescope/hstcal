# HSTCAL

[![Jenkins CI](https://ssbjenkins.stsci.edu/job/STScI/job/hstcal/job/master/badge/icon)](https://ssbjenkins.stsci.edu/job/STScI/job/hstcal/job/master/)

Calibration software for HST/WFC3, HST/ACS, and HST/STIS.

Nightly regression test results are available only from within the STScI network at this time.
https://plwishmaster.stsci.edu:8081/job/RT/job/hstcal/test_results_analyzer/

## Install using Conda (Complete Calibration Environment)

HSTCAL can be obtained as part of the
**Space Telescope Environment (stenv)** at <https://stenv.readthedocs.io/en/latest/>, 
a conda environment.  The instructions found at this URL describe choosing Miniconda 
or Anaconda (if not already installed), choosing the ``stenv`` release, building the ``stenv``
environment from the YAML file, and activating your new environment. Once your environment
is activated, you can invoke any of the HSTCAL executables.

```bash
# Use HSTCAL
$ calacs.e [...]
$ calwf3.e [...]
$ cs[...].e
```

## Install using Conda (HSTCAL only)

The HSTCAL package is resident on conda-forge
<https://anaconda.org/conda-forge/hstcal>, and you can use the following
command to perform the installation.

```bash
$ conda install -c conda-forge hstcal==X.Y.Z
```
The `X.Y.Z` is the desired version number.

## Source Installation

**Note:
This involves compilation of C code which requires that a C compiler and all dependent libraries be
installed on your system prior to installing this package; specifically,**
  - cfitsio >=3.430 (https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
  - gcc >=4.4.7 (https://gcc.gnu.org)
  - openmp (http://www.openmp.org/)
  - pkg-config (https://www.freedesktop.org/wiki/Software/pkg-config)
  - python >=2.7 (https://www.python.org)
  - GCC may be supplemented by Clang under the following conditions:
    * Clang must be compiled with OpenMP support
    * A Fortran compiler must be available in $PATH

Instructions:

1. Clone the package from github onto your local system using:

  `git clone https://github.com/spacetelescope/hstcal.git`

2. Compile the code using the instructions provided in `INSTALL.md`.
