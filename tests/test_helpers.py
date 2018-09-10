"""Test the test helpers."""
import subprocess

from astropy.io import fits

from .helpers import (HAS_CALXXX, use_calacs, use_calwf3, use_calstis,
                      calref_from_image)


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


def test_calref_from_image(_jail):
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
    hdu.header['PCTETAB'] = 'jref$dummy_file_2.fits'
    datafile = 'dummy_raw.fits'
    hdu.writeto(datafile, overwrite=True)

    # Get all the references files at this time regardless of 'CORR' setting
    ref_files = calref_from_image(datafile)
    assert sorted(ref_files) == ['dummy_file_1.fits',
                                 'jref$dummy_file_2.fits',
                                 'jref$t3n1116nj_bpx.fits',
                                 'jref$xa81715gj_ccd.fits']
