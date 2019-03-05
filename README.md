# HSTCAL

[![Travis CI](https://travis-ci.org/spacetelescope/hstcal.svg?branch=master)](https://travis-ci.org/spacetelescope/hstcal)
[![Jenkins CI](https://ssbjenkins.stsci.edu/job/STScI/job/hstcal/job/master/badge/icon)](https://ssbjenkins.stsci.edu/job/STScI/job/hstcal/job/master/)

Calibration software for HST/WFC3, HST/ACS, and HST/STIS.

Nightly regression test results are available only from within the STScI network at this time.

https://plwishmaster.stsci.edu:8081/job/RT/job/hstcal/test_results_analyzer/


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
$ cs[...].e
```

## Source Installation

To install HSTCAL from source please follow the instructions provided in
`INSTALL.md`.


## JupyterHub Access

*Note:** This is currently still in research-and-development stage and is subject to change.

To run a pre-installed pipeline in JupyterHub:

* Click on https://dev.science.stsci.edu/hub/spawn?image=793754315137.dkr.ecr.us-east-1.amazonaws.com/datb-tc-pipeline-nb:hstdp-snapshot and sign in.
* Click "Terminal" to:
    * Do a `which calacs.e` to see if CALACS is installed.
      You can repeat this for other HSTCAL executables, as desired.
    * Do a `calacs.e --version` to see which CALACS version is installed.
      You can repeat this for other HSTCAL executables, if applicable, as desired.
    * Run `pip freeze` to see what Python packages are installed.
    * Install any optional Python packages using `pip install`.
    * You can download the necessary data files using HTTP/HTTPS protocol.
    * Set up your `jref`, `iref, etc. as desired.
    * Run the pipeline.

Latest release of any packages is not guaranteed in this environment. Amazon Web Services charges may apply.