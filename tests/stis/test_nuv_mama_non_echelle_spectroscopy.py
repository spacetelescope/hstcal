import subprocess
from ci_watson.artifactory_helpers import get_bigdata

from ..helpers import BaseSTIS


class TestNUVMAMANonEchelleSpectroscopy(BaseSTIS):
    """NUV-MAMA non-echelle spectroscopy."""
    detector = 'nuv-mama'

    def test_nuv_mama_non_echelle_spectroscopy(self):
        """Test NUV-MAMA non-echelle spectroscopy."""
        raw_file = 'obq002030_raw.fits'
        wav_file = 'obq002030_wav.fits'

        # Prepare input files.
        self.get_input_file(raw_file)
        get_bigdata('scsb-hstcal', self.env, 'stis', 'nuv-mama', 'input',
                    wav_file)

        # Run CALSTIS (equivalent to stistools.calstis.calstis)
        subprocess.call(['cs0.e', raw_file, '-v'])

        # Compare results
        outputs = [('obq002030_flt.fits', 'obq002030_flt.fits'),
                   ('obq002030_x1d.fits', 'obq002030_x1d.fits')]
        self.compare_outputs(outputs)
