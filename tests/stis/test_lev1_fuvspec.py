import subprocess
from ci_watson.artifactory_helpers import get_bigdata

from ..helpers import BaseSTIS


class TestLev1FUVSpec(BaseSTIS):
    """Level 1 FUV-MAMA spectroscopy."""
    detector = 'fuv-mama'

    def test_lev1(self):
        """This was ``calstis_lev1_FUVspec``."""
        raw_file = 'o5cl02040_raw.fits'
        wav_file = 'o5cl02040_wav.fits'

        # Prepare input files.
        self.get_input_file(raw_file)
        get_bigdata(self.env, 'stis', 'fuv-mama', 'input', wav_file)

        # Run CALSTIS (equivalent to stistools.calstis.calstis)
        subprocess.call(['cs0.e', raw_file, '-v'])

        # Compare results
        outputs = [('o5cl02040_flt.fits', 'calstis_lev1_FUVspec_flt.fits'),
                   ('o5cl02040_x1d.fits', 'calstis_lev1_FUVspec_x1d.fits')]
        self.compare_outputs(outputs)
