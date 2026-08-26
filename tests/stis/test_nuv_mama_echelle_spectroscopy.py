import subprocess
from ci_watson.artifactory_helpers import get_bigdata

from ..helpers import BaseSTIS


class TestNUVMAMAEchelleSpectroscopy(BaseSTIS):
    """NUV-MAMA echelle spectroscopy."""
    detector = 'nuv-mama'

    def test_nuv_mama_echelle_spectroscopy(self):
        """Test NUV-MAMA echelle spectroscopy."""
        raw_file = 'o61l01030_raw.fits'
        wav_file = 'o61l01030_wav.fits'

        # Prepare input files.
        self.get_input_file(raw_file)
        get_bigdata('scsb-hstcal', self.env, 'stis', 'nuv-mama', 'input',
                    wav_file)

        # Run CALSTIS (equivalent to stistools.calstis.calstis)
        subprocess.call(['cs0.e', raw_file, '-v'])

        # Compare results
        outputs = [('o61l01030_flt.fits', 'o61l01030_flt.fits'),
                   ('o61l01030_sfl.fits', 'o61l01030_sfl.fits'),
                   ('o61l01030_x1d.fits', 'o61l01030_x1d.fits')]
        self.compare_outputs(outputs)
