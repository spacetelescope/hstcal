import subprocess
from ci_watson.artifactory_helpers import get_bigdata

from ..helpers import BaseSTIS


class TestLev1FUVSpec(BaseSTIS):
    """Level 1 CCD spectroscopy."""
    detector = 'ccd'

    def test_lev1(self):
        raw_file = 'o55108010_raw.fits'
        wav_file = 'o55108010_wav.fits'

        # Prepare input files.
        self.get_input_file(raw_file)
        get_bigdata('scsb-hstcal', self.env, 'stis', 'ccd', 'input',
                    wav_file)

        # Run CALSTIS (equivalent to stistools.calstis.calstis)
        subprocess.call(['cs0.e', raw_file, '-v'])

        # Compare results
        outputs = [('o55108010_flt.fits', 'calstis_ccdspec_flt.fits'),
                   ('o55108010_sx1.fits', 'calstis_ccdspec_sx1.fits'),
                   ('o55108010_crj.fits', 'calstis_ccdspec_crj.fits'),
                   ('o55108010_sx2.fits', 'calstis_ccdspec_sx2.fits')]
        self.compare_outputs(outputs)
