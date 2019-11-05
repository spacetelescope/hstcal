# HSTCAL

[![Travis CI](https://travis-ci.org/spacetelescope/hstcal.svg?branch=master)](https://travis-ci.org/spacetelescope/hstcal)
[![Jenkins CI](https://ssbjenkins.stsci.edu/job/STScI/job/hstcal/job/master/badge/icon)](https://ssbjenkins.stsci.edu/job/STScI/job/hstcal/job/master/)

Calibration software for HST/WFC3, HST/ACS, and HST/STIS.

Nightly regression test results are available only from within the STScI network at this time.
https://plwishmaster.stsci.edu:8081/job/RT/job/hstcal/test_results_analyzer/

## Install using Conda

HSTCAL can be obtained as part of the
**Astroconda Pipeline Builds** at <https://astroconda.readthedocs.io/en/latest/releases.html>.  
The simplest way to get the most recent environment would be:

1.  Install Miniconda or Anaconda (if not already installed) using the instructions at <https://astroconda.readthedocs.io/en/latest/getting_started.html>.  

2.  Create a conda environment with the HSTCAL included:

```bash
# Add the AstroConda channel to ~/.condarc
$ conda channel --add http://ssb.stsci.edu/astroconda

# Create a dedicated calibration environment
$ conda create -n calib python=3 stsci-hst

# Activate the environment
$ source activate calib

# Use HSTCAL
$ calacs.e [...]
$ calwf3.e [...]
$ cs[...].e
```

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
