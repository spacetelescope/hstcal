"""HSTCAL regression test helpers."""
import os
import shutil
import subprocess
import sys

import pytest
import requests
from astropy.io import fits

__all__ = ['remote_data', 'use_calacs', 'use_calwf3', 'use_calstis',
           'HAS_CALXXX', 'download_file_cgi', 'download_ref_ftp']

PREFIX = '/tmp/hstcal'  # As set in .travis.yml
HAS_CALXXX = {}   # Set by set_exe_marker()


def set_exe_marker(instrument):
    """Set pytest marker for given instrument calibration executable."""
    global HAS_CALXXX
    instrument = instrument.lower()

    if instrument == 'acs':
        calxxx = 'calacs.e'
    elif instrument == 'wfc3':
        calxxx = 'calwf3.e'
    elif instrument == 'stis':
        calxxx = 'cs0.e'
    else:
        raise ValueError('{} is not supported'.format(instrument))

    HAS_CALXXX[instrument] = False
    travis_path = os.path.join(PREFIX, calxxx)

    # This is for Travis CI installation check.
    if os.path.isfile(travis_path):
        HAS_CALXXX[instrument] = travis_path

    # This is for local installation check for Python 3.3 or later.
    elif sys.version_info >= (3, 3):
        if os.path.isfile(shutil.which(calxxx)):
            HAS_CALXXX[instrument] = calxxx

    # This is for local installation check before Python 3.3.
    # NOTE: This will fail on Windows!
    else:
        try:
            subprocess.check_call(['which', calxxx])
        except Exception:
            pass
        else:
            HAS_CALXXX[instrument] = calxxx


# pytest marker to mark tests which get data from the web
remote_data = pytest.mark.remote_data

# pytest markers to mark tests that use CALXXX

set_exe_marker('acs')
use_calacs = pytest.mark.skipif(not HAS_CALXXX['acs'], reason='no CALACS')

set_exe_marker('wfc3')
use_calwf3 = pytest.mark.skipif(not HAS_CALXXX['wfc3'], reason='no CALWF3')

set_exe_marker('stis')
use_calstis = pytest.mark.skipif(not HAS_CALXXX['stis'], reason='no CALSTIS')


def _download_file(url, filename, filemode='wb', timeout=None):
    """Generic remote data download."""
    r = requests.get(url, timeout=timeout)
    with open(filename, filemode) as fout:
        fout.write(r.content)


def _download_crds(hdr, key, sep='$', timeout=None):
    """Download a CRDS file from FTP to current directory."""
    refname = hdr[key].strip()

    # Local file: Do nothing.
    if sep not in refname:
        return

    s = refname.split(sep)
    filename = s[1]

    # CRDS file for given name never changes, so no need to re-download.
    if os.path.exists(filename):
        return

    url = 'ftp://ftp.stsci.edu/cdbs/{}/{}'.format(s[0], filename)
    _download_file(url, filename, timeout=timeout)


def download_file_cgi(tree, project, filename, filemode='wb', timeout=30):
    """
    Download remote data to current directory using CGI script interface.

    Parameters
    ----------
    tree : {'rt', 'rtx', 'null'}
        Test tree:

        * rt = dev
        * rtx = public
        * null = neither, grab first match (e.g., for STAK)

    project : str
        Path containing data in the test tree.
        For example, ``hstcal/acs/calacs_e``.

    filename : str
        Filename of the data.

    filemode : str
        Mode of the output file writer. By default, it saves as binary
        and overwrites existing file.

    timeout : number or `None`
        This is not a time limit on the entire response download;
        rather, an exception is raised if the server has not issued
        a response for timeout seconds (more precisely, if no bytes
        have been received on the underlying socket for timeout seconds).
        If no timeout is specified explicitly, requests do not time out.

    """
    url = ('http://ssb.stsci.edu/cgi-bin/remote_testing.cgi?'
           'tree={}&project={}&name={}'.format(tree, project, filename))
    _download_file(url, filename, filemode=filemode, timeout=timeout)


def download_ref_ftp(input_image, timeout=30):
    """
    Download necessary CRDS reference files for given input image.
    The image primary header will be scanned for ``*CORR``.
    Reference file for any ``CORR`` set to ``PERFORM`` will be
    downloaded to current directory from STScI FTP site.

    Parameters
    ----------
    input_image : str
        Input image filename.

    timeout : number or `None`
        See :func:`download_file_cgi`.

    """

    # TODO: Add additional mapping as needed.
    # Map mandatory reference file for instrument/detector combo.
    det_lookup = {
        ('ACS', 'WFC'): ['CCDTAB'],
        ('ACS', 'HRC'): ['CCDTAB'],
        ('ACS', 'SBC'): [''],
        ('WFC3', 'UVIS1'): [''],
        ('WFC3', 'UVIS2'): [''],
        ('WFC3', 'IR'): [''],
        ('STIS', 'CCD'): [''],
        ('STIS', 'FUVMAMA'): [''],
        ('STIS', 'NUVMAMA'): ['']}

    # TODO: Add additional mapping as needed.
    # Map *CORR to associated CRDS reference file.
    corr_lookup = {
        'DQICORR': ['BPIXTAB'],
        'ATODCORR': ['ATODTAB'],
        'BLEVCORR': ['OSCNTAB'],
        'SINKCORR': ['SNKCFILE'],
        'BIASCORR': ['BIASFILE'],
        'PCTECORR': ['PCTETAB', 'DRKCFILE'],
        'FLSHCORR': ['FLSHFILE'],
        'CRCORR': ['CRREJTAB'],
        'SHADCORR': ['SHADFILE'],
        'DARKCORR': ['DARKFILE'],
        'FLATCORR': ['PFLTFILE', 'DFLTFILE', 'LFLTFILE'],
        'PHOTCORR': ['IMPHTTAB'],
        'GLINCORR': ['MLINTAB']}

    hdr = fits.getheader(input_image, ext=0)

    for key in det_lookup[(hdr['INSTRUME'], hdr['DETECTOR'])]:
        _download_crds(hdr, key, timeout=timeout)

    for step in corr_lookup:
        # Not all images have the CORR step and it is not always on.
        if (step not in hdr) or (hdr[step].strip().upper() != 'PERFORM'):
            continue

        for key in corr_lookup[step]:
            _download_crds(hdr, key, timeout=timeout)
