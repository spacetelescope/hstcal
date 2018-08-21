"""HSTCAL regression test helpers."""

import os
from distutils.spawn import find_executable

import pytest
from ci_watson.artifactory_helpers import get_bigdata
from ci_watson.hst_helpers import raw_from_asn, ref_from_image, download_crds

from astropy.io import fits
from astropy.io.fits import FITSDiff

__all__ = ['use_calacs', 'use_calwf3', 'use_calstis', 'calref_from_image',
           'BaseACS', 'BaseWFC3', 'BaseSTIS', 'fix_keywords']

HAS_CALXXX = {}   # Set by set_exe_marker()


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

    cal = find_executable(calxxx)

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
    header of the given input image, necessary for calibration; i.e.,
    only those associated with ``*CORR`` set to ``PERFORM`` will be
    considered.
    """
    # NOTE: Add additional mapping as needed.
    # Map mandatory CRDS reference file for instrument/detector combo.
    # This is for file not tied to any particular *CORR or used throughout.
    det_lookup = {
        ('ACS', 'WFC'): ['CCDTAB'],
        ('ACS', 'HRC'): ['CCDTAB'],
        ('ACS', 'SBC'): ['CCDTAB'],
        ('WFC3', 'UVIS1'): ['CCDTAB'],
        ('WFC3', 'UVIS2'): ['CCDTAB'],
        ('WFC3', 'IR'): ['CCDTAB'],
        ('STIS', 'CCD'): ['CCDTAB'],
        ('STIS', 'FUV-MAMA'): ['CCDTAB', 'DISPTAB', 'INANGTAB', 'APDESTAB',
                               'SPTRCTAB'],
        ('STIS', 'NUV-MAMA'): ['CCDTAB', 'DISPTAB', 'INANGTAB', 'APDESTAB',
                               'SPTRCTAB']}

    # NOTE: Add additional mapping as needed.
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
        if (step not in hdr) or (hdr[step].strip().upper() != 'PERFORM'):
            continue

        ref_files += ref_from_image(input_image, corr_lookup[step])

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
        dest = get_bigdata(self.env, self.instrument, self.detector,
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

        first_pass = ('JENKINS_URL' in os.environ and
                      'ssbjenkins' in os.environ['JENKINS_URL'])

        for raw in all_raws:
            ref_files = calref_from_image(raw)

            for ref_file in ref_files:
                # Special reference files that live with inputs.
                if ('$' not in ref_file and
                        os.path.basename(ref_file) == ref_file):
                    get_bigdata(self.env, self.instrument, self.detector,
                                'input', ref_file)
                    continue

                # Jenkins cannot see Central Storage on push event,
                # and somehow setting, say, jref to "." does not work anymore.
                # So, we need this hack.
                if '$' in ref_file and first_pass:
                    first_pass = False
                    if not os.path.isdir('/grp/hst/cdbs'):
                        ref_path = os.path.dirname(dest) + os.sep
                        var = ref_file.split('$')[0]
                        os.environ[var] = ref_path  # hacky hack hack

                # Download reference files, if needed only.
                download_crds(ref_file, timeout=self.timeout)

    def compare_outputs(self, outputs, atol=0, rtol=1e-7, raise_error=True,
                        ignore_keywords_overwrite=None):
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

        Returns
        -------
        report : str
            Report from ``fitsdiff``.
            This is part of error message if ``raise_error=True``.

        """
        all_okay = True
        creature_report = ''

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

            if not fdiff.identical and all_okay:
                all_okay = False

        if not all_okay and raise_error:
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


def fix_keywords(key_fix_directive, file_pattern='*.fits', upload_login=None,
                 upload_env='dev', upload_type=None):
    """
    Fix keywords that trip up Jenkins CI.
    For example, reference files hardcoded to
    ``/grp/hst/...`` that is not always accessible.

    This is not used directly by CI, but rather is
    used to prepare files for CI runs.

    Parameters
    ----------
    key_fix_directive : dict
        Dictionary defining the directive on how to
        fix up keyword(s) in the selected file(s)
        in the form of ``{key: (ext, wrong, correct)}``.
        Keyword is case-insensitive but value is case-sensitive.
        Example::

            {'PCTETAB': (0, '/grp/hst/cdbs/jref/filename.fits',
                         'jref$filename.fits'),
             'FLSHFILE': (0, 'jref$localfile.fits',
                          'localfile.fits')}

    file_pattern : str
        File selection pattern used by :py:mod:`glob`.

    upload_login : tuple of str or `None`
        ``(username, token)`` is needed for uploading
        corrected file to Artifactory. If `None` or
        invalid, upload is skipped.

    upload_env : {'dev', 'stable'}
        Testing environment of upload destination.
        Default is 'dev'.

    upload_type : str or `None`
        A valid type (e.g., 'input' or 'truth') is
        needed for upload. Otherwise, upload is skipped.

    """
    import glob
    import subprocess

    do_upload = (upload_login is not None and upload_type is not None and
                 len(upload_login) == 2)
    base_url = ('https://bytesalad.stsci.edu/artifactory/'
                'scsb-hstcal/{}'.format(upload_env))
    status_msgs = []

    for fname in glob.iglob(file_pattern):
        updated = []

        with fits.open(fname, mode='update') as pf:
            ins = pf[0].header['INSTRUME'].strip().lower()
            det = pf[0].header['DETECTOR'].strip().lower()
            for key, (ext, wrong, correct) in key_fix_directive.items():
                hdr = pf[ext].header
                if key in hdr and hdr[key].strip() == wrong:
                    pf[ext].header[key] = correct
                    updated.append(key)

        if len(updated) > 0:
            status_msgs.append('{} updated: {}'.format(fname, updated))
            if do_upload:
                url = '{}/{}/{}/{}/'.format(base_url, ins, det, upload_type)
                retcode = subprocess.call([
                    'curl', '-u{}:{}'.format(*upload_login), '-T', fname, url])
                if retcode == 0:
                    status_msgs.append('\tUploaded to {}'.format(url))
                else:
                    status_msgs.append(
                        '\tUpload to {} failed with {}'.format(url, retcode))

    # Print status messages at the very end so they won't get drown out
    # by Artifactory JSON output.
    print()
    print(os.linesep.join(status_msgs))
