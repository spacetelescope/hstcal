"""Test the test helpers."""
from __future__ import absolute_import, division, print_function

import os
import subprocess

import pytest
from astropy.io import fits
from astropy.table import Table

from ..helpers import (HAS_CALXXX, use_calacs, use_calwf3, use_calstis,
                       ref_from_image, raw_from_asn, remote_data,
                       download_file_cgi, download_crds)


@use_calacs
def test_use_calacs():
    """Call CALACS without need for file access."""
    retval = subprocess.call([HAS_CALXXX['acs'], '--version'])
    assert retval == 0


@use_calwf3
def test_use_calwf3():
    """Call CALWF3 without need for file access."""
    retval = subprocess.call([HAS_CALXXX['wfc3'], '--version'])
    assert retval == 0


@use_calstis
def test_use_calstis():
    """Call CALSTIS without need for file access."""
    retval = subprocess.call([HAS_CALXXX['stis'], '--version'])
    assert retval == 0


def test_ref_from_image(tmpdir):
    # Make a dummy input file (to avoid package data headache)
    hdu = fits.PrimaryHDU()
    hdu.header['INSTRUME'] = 'ACS'
    hdu.header['DETECTOR'] = 'WFC'
    hdu.header['CCDTAB'] = 'jref$xa81715gj_ccd.fits   '  # Compulsory
    hdu.header['DQICORR'] = 'PERFORM'
    hdu.header['BPIXTAB'] = 'jref$t3n1116nj_bpx.fits'  # For DQICORR
    hdu.header['DARKCORR'] = 'perform'
    hdu.header['DARKFILE'] = 'dummy_file_1.fits'  # Some local file
    hdu.header['SINKCORR'] = 'PERFORM'
    hdu.header['SNKCFILE'] = 'N/A'
    hdu.header['PCTECORR'] = 'OMIT'
    hdu.header['PCTETAB'] = 'jref$dummy_file_2.fits'  # Should not matter
    datafile = tmpdir.join('dummy_raw.fits').strpath
    hdu.writeto(datafile, overwrite=True)

    ref_files = ref_from_image(datafile)
    assert sorted(ref_files) == ['dummy_file_1.fits',
                                 'jref$t3n1116nj_bpx.fits',
                                 'jref$xa81715gj_ccd.fits']


def test_raw_from_asn(tmpdir):
    # Make a dummy input file (to avoid package data headache)
    tab = Table()
    tab['MEMNAME'] = ['J6LQ01NAQ', 'J6LQ01NDQ', 'J6LQ01011']
    tab['MEMTYPE'] = ['EXP-CRJ', 'EXP-CRJ', 'PROD-CRJ']
    tab['MEMPRSNT'] = [True, True, True]
    datafile = tmpdir.join('dummy_asn.fits').strpath
    tab.write(datafile, format='fits', overwrite=True)

    raw_files = raw_from_asn(datafile)
    assert raw_files == ['j6lq01naq_raw.fits', 'j6lq01ndq_raw.fits']

    # Make sure do not download existing file.
    # This will fail if download is attemped.
    download_crds('jref', datafile)


@remote_data
class TestRemoteData(object):
    """Test remote data acess."""
    prevdir = os.getcwd()  # Need to be here for teardown access

    @pytest.fixture(autouse=True)
    def setup_class(self, tmpdir):
        """Run test in own dir to mimic how CALXXX works."""
        p = tmpdir.mkdir('utils').strpath
        os.chdir(p)

    def teardown_class(self):
        """Change back to test directory."""
        os.chdir(self.prevdir)

    def test_marker(self, pytestconfig):
        """Check if test not marked for remote data somehow runs this."""
        if pytestconfig.getoption('remote_data') is not True:
            pytest.fail('@remote_data was not skipped')

    def test_acs_input_data(self):
        """
        Test CALACS input data access.
        This tests central storage access on internal machine and
        CGI (HTTP) access for Travis CI or external machine.
        """
        # Get a small file that we know exists
        filename = 'newhrc_osc.fits'
        download_file_cgi('rt', 'hstcal/acs/calacs_e', filename)
        assert os.path.isfile(filename)

        # Get the path only; do not download
        filename = 'jbdf10ykq_flt_ref.fits'
        filepath = download_file_cgi(
            'rt', 'hstcal/acs/calacs_e/ref', filename, allow_remote_ref=True)
        assert (filepath.startswith(('http', '/eng/ssb2')) and
                filepath.endswith(filename))
        assert not os.path.exists(filename)

    def test_ftp_jref(self):
        """
        Test CRDS FTP access for ACS (JREF).
        This mechanism is used by Travis CI or undefined JREF.
        For JREF pointing to local dir or central storage, HSTCAL already
        knows how to handle that; hence not tested here.
        """
        filename = 't3n1116nj_bpx.fits'
        download_crds('jref', filename)
        assert os.path.isfile(filename)
