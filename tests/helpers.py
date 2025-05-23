"""HSTCAL regression test helpers."""

import os
import shutil
from functools import partial

import pytest
from ci_watson.artifactory_helpers import (
    generate_upload_params, generate_upload_schema)
from ci_watson.artifactory_helpers import get_bigdata as _get_bigdata
from ci_watson.hst_helpers import raw_from_asn, ref_from_image, download_crds

from astropy.io import fits
from astropy.io.fits import FITSDiff

__all__ = ['use_calacs', 'use_calwf3', 'use_calstis', 'calref_from_image',
           'BaseACS', 'BaseWFC3', 'BaseSTIS']

HAS_CALXXX = {}   # Set by set_exe_marker()

# Overload generic get_bigdata to include repo root dir.
# This is to accomodate developers who have to run big data tests across
# several repositories using the same TEST_BIGDATA env var.
get_bigdata = partial(_get_bigdata, 'scsb-hstcal')


# NOTE: This is because HSTCAL allows partial installation via --targets
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
    else:  # pragma: no cover
        raise ValueError('{} is not supported'.format(instrument))

    cal = shutil.which(calxxx)

    # This is for installation check.
    if os.path.isfile(cal):
        HAS_CALXXX[instrument] = cal
    else:  # pragma: no cover
        HAS_CALXXX[instrument] = False


# pytest markers to mark tests that use CALXXX

set_exe_marker('acs')
use_calacs = pytest.mark.skipif(not HAS_CALXXX['acs'], reason='no CALACS')

set_exe_marker('wfc3')
use_calwf3 = pytest.mark.skipif(not HAS_CALXXX['wfc3'], reason='no CALWF3')

set_exe_marker('stis')
use_calstis = pytest.mark.skipif(not HAS_CALXXX['stis'], reason='no CALSTIS')


def calref_from_image(input_image):
    """
    Return a list of reference filenames, as defined in the primary
    header of the given input image, necessary for calibration.
    """
    # NOTE: Add additional mapping as needed.
    # Map mandatory CRDS reference file for instrument/detector combo.
    # This is for file not tied to any particular *CORR or used throughout.
    det_lookup = {
        ('ACS', 'WFC'): ['CCDTAB'],
        ('ACS', 'HRC'): ['CCDTAB'],
        ('ACS', 'SBC'): ['CCDTAB'],
        ('WFC3', 'UVIS'): ['CCDTAB'],
        ('WFC3', 'IR'): ['CCDTAB'],
        ('STIS', 'CCD'): ['CCDTAB'],
        ('STIS', 'FUV-MAMA'): ['CCDTAB', 'DISPTAB', 'INANGTAB', 'APDESTAB',
                               'SPTRCTAB'],
        ('STIS', 'NUV-MAMA'): ['CCDTAB', 'DISPTAB', 'INANGTAB', 'APDESTAB',
                               'SPTRCTAB']}

    # NOTE: Add additional mapping as needed.
    # Map *CORR to associated CRDS reference file.
    corr_lookup = {
        'DQICORR': ['BPIXTAB', 'SNKCFILE'],
        'ATODCORR': ['ATODTAB'],
        'BLEVCORR': ['OSCNTAB'],
        'SINKCORR': ['SNKCFILE'],
        'BIASCORR': ['BIASFILE', 'SATUFILE'],
        'PCTECORR': ['PCTETAB', 'DRKCFILE', 'BIACFILE'],
        'FLSHCORR': ['FLSHFILE'],
        'CRCORR': ['CRREJTAB'],
        'SHADCORR': ['SHADFILE'],
        'DARKCORR': ['DARKFILE', 'TDCTAB'],
        'FLATCORR': ['PFLTFILE', 'DFLTFILE', 'LFLTFILE'],
        'PHOTCORR': ['IMPHTTAB'],
        'LFLGCORR': ['MLINTAB'],
        'GLINCORR': ['MLINTAB'],
        'NLINCORR': ['NLINFILE'],
        'ZSIGCORR': ['DARKFILE', 'NLINFILE'],
        'WAVECORR': ['LAMPTAB', 'WCPTAB', 'SDCTAB'],
        'SGEOCORR': ['SDSTFILE'],
        'X1DCORR': ['XTRACTAB', 'SDCTAB'],
        'SC2DCORR': ['CDSTAB', 'ECHSCTAB', 'EXSTAB', 'RIPTAB', 'HALOTAB',
                     'TELTAB', 'SRWTAB'],
        'BACKCORR': ['XTRACTAB'],
        'FLUXCORR': ['APERTAB', 'PHOTTAB', 'PCTAB', 'TDSTAB']}

    hdr = fits.getheader(input_image, ext=0)
    ref_files = ref_from_image(
        input_image, det_lookup[(hdr['INSTRUME'], hdr['DETECTOR'])])

    for step in corr_lookup:
        # Not all images have the CORR step and it is not always on.
        # Download ALL reference files associated with a calibration
        # step present in the header
        if step not in hdr:
            continue

        single_step_files = ref_from_image(input_image, corr_lookup[step])
        if single_step_files:
            ref_files += single_step_files

    return ref_files


# Base classes for actual tests.
# NOTE: Named in a way so pytest will not pick them up here.
# NOTE: bigdata marker requires TEST_BIGDATA environment variable to
#       point to a valid big data directory, whether locally or on Artifactory.
# NOTE: envopt would point tests to "dev" or "stable".
# NOTE: _jail fixture ensures each test runs in a clean tmpdir.
@pytest.mark.bigdata
@pytest.mark.usefixtures('_jail', 'envopt')
class BaseCal:
    # Timeout in seconds for file downloads.
    timeout = 30

    # To be defined by instrument.
    instrument = ''
    ignore_keywords = []  # List of FITS keywords to ignore during fitsdiff

    # To be defined by test class in actual test modules.
    detector = ''

    # Where artifacts go in Artifactory.
    # Set to None if you do not want this functionality.
    results_root = 'scsb-hstcal-results'

    @pytest.fixture(autouse=True)
    def setup_class(self, envopt):
        """
        Class-level setup that is done at the beginning of the test.

        Parameters
        ----------
        envopt : {'dev', 'stable'}
            This is a ``pytest`` fixture that defines the test
            environment in which input and truth files reside.

        """
        self.env = envopt

    def get_input_file(self, filename):
        """
        Copy input file (ASN, RAW, etc) into the working directory.
        If ASN is given, RAW files in the ASN table is also copied.
        The associated CRDS reference files are also copied or
        downloaded, if necessary.

        Data directory layout for HSTCAL::

            instrument/
                detector/
                    input/
                    truth/

        .. note::

            If given filename has ".orig" suffix, the suffix is
            also automatically removed.

        Parameters
        ----------
        filename : str
            Filename of the ASN/RAW/etc to copy over, along with its
            associated files.

        """
        # Copy over main input file.
        get_bigdata(self.env, self.instrument, self.detector,
                    'input', filename)

        # For historical reason, need to remove ".orig" suffix if it exists.
        if filename.endswith('.orig'):
            newfilename = filename.rstrip('.orig')
            os.rename(filename, newfilename)
            filename = newfilename

        if filename.endswith('_asn.fits'):
            all_raws = raw_from_asn(filename)
            for raw in all_raws:  # Download RAWs in ASN.
                get_bigdata(self.env, self.instrument, self.detector,
                            'input', raw)
        else:
            all_raws = [filename]

        for raw in all_raws:
            ref_files = calref_from_image(raw)

            for ref_file in ref_files:
                # Special reference files that live with inputs.
                if ('$' not in ref_file and
                        os.path.basename(ref_file) == ref_file):
                    get_bigdata(self.env, self.instrument, self.detector,
                                'input', ref_file)
                    continue

                # Download reference files, if needed only.
                download_crds(ref_file)

    def compare_outputs(self, outputs, atol=0, rtol=1e-7, raise_error=True,
                        ignore_keywords_overwrite=None, verbose=True):
        """
        Compare CALXXX output with "truth" using ``fitsdiff``.

        Parameters
        ----------
        outputs : list of tuple
            A list of tuples, each containing filename (without path)
            of CALXXX output and truth, in that order. Example::

                [('output1.fits', 'truth1.fits'),
                 ('output2.fits', 'truth2.fits'),
                 ...]

        atol, rtol : float
            Absolute and relative tolerance for data comparison.

        raise_error : bool
            Raise ``AssertionError`` if difference is found.

        ignore_keywords_overwrite : list of str or `None`
            If not `None`, these will overwrite
            ``self.ignore_keywords`` for the calling test.

        verbose : bool
            Print extra info to screen.

        Returns
        -------
        report : str
            Report from ``fitsdiff``.
            This is part of error message if ``raise_error=True``.

        """
        all_okay = True
        creature_report = ''
        updated_outputs = []  # To track outputs for Artifactory JSON schema

        if ignore_keywords_overwrite is None:
            ignore_keywords = self.ignore_keywords
        else:
            ignore_keywords = ignore_keywords_overwrite

        for actual, desired in outputs:
            desired = get_bigdata(self.env, self.instrument, self.detector,
                                  'truth', desired)
            fdiff = FITSDiff(actual, desired, rtol=rtol, atol=atol,
                             ignore_keywords=ignore_keywords)
            creature_report += fdiff.report()

            if not fdiff.identical:
                all_okay = False
                # Only keep track of failed results which need to
                # be used to replace the truth files (if OK).
                updated_outputs.append((actual, desired))

        if not all_okay:
            if self.results_root is not None:  # pragma: no cover
                schema_pattern, tree, testname = generate_upload_params(
                    self.results_root, updated_outputs, verbose=verbose)
                generate_upload_schema(schema_pattern, tree, testname)

            if raise_error:
                raise AssertionError(os.linesep + creature_report)

        return creature_report


@use_calacs
class BaseACS(BaseCal):
    """Base class for all ACS tests."""
    instrument = 'acs'
    ignore_keywords = ['filename', 'date', 'iraf-tlm', 'fitsdate',
                       'opus_ver', 'cal_ver', 'proctime', 'history']


@use_calwf3
class BaseWFC3(BaseCal):
    """Base class for all WFC3 tests."""
    instrument = 'wfc3'
    ignore_keywords = ['filename', 'date', 'iraf-tlm', 'fitsdate',
                       'extname', 'rootname', 'cal_ver', 'history']


@use_calstis
class BaseSTIS(BaseCal):
    """Base class for all STIS tests."""
    instrument = 'stis'
    ignore_keywords = ['filename', 'date', 'cal_ver', 'history']
