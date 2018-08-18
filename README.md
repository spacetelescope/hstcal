# HSTCAL

[![Travis CI](https://travis-ci.org/spacetelescope/hstcal.svg?branch=master)](https://travis-ci.org/spacetelescope/hstcal)

[![Jenkins CI](https://ssbjenkins.stsci.edu/job/STScI/job/hstcal/job/master/badge/icon)](https://ssbjenkins.stsci.edu/job/STScI/job/hstcal/job/master/)

Calibration software for HST/WFC3, HST/ACS, and HST/STIS.

## Install using Conda

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
$ cal[...].e
```

## Source Installation

To install HSTCAL from source please follow the instructions provided in `INSTALL.md`.
