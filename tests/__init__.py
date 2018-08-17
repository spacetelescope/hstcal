"""
This directory contains HSTCAL regression tests, which are written in
Python but call the C executables.

These environment variables must be set up prior::

    export BIG_DATA=/path/to/bigdata/scsb_hstcal
    export jref=/grp/hst/cdbs/jref/
    export iref=/grp/hst/cdbs/iref/
    export oref=/grp/hst/cdbs/oref/

Example::

    pytest tests --big-data [--slow] [-sv] [--basetemp=/my/tmp]

"""

import sys

__minimum_python_version__ = '3.5'


class UnsupportedPythonError(Exception):
    pass


# This is the same check as the one at the top of setup.py
if sys.version_info < tuple((int(val) for val in __minimum_python_version__.split('.'))):  # noqa
    raise UnsupportedPythonError("HSTCAL tests do not support Python < {}".format(__minimum_python_version__))  # noqa
