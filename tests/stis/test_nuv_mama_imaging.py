import subprocess
from ci_watson.artifactory_helpers import get_bigdata

from ..helpers import BaseSTIS


class TestNUVMAMAImage(BaseSTIS):
    """Level 1 NUV-MAMA imaging."""
    detector = 'nuv-mama'

    def test_nuv_mama_imaging(self):
        """Test NUV-MAMA imaging."""
        raw_file = 'odnv02ahq_raw.fits'
        spt_file = 'odnv02ahq_spt.fits'

        # Prepare input files.
        self.get_input_file(raw_file)
        get_bigdata('scsb-hstcal', self.env, 'stis', 'nuv-mama', 'input',
                    spt_file)

        # Run CALSTIS (equivalent to stistools.calstis.calstis)
        subprocess.call(['cs0.e', raw_file, '-v'])

        # Compare results
        outputs = [('odnv02ahq_flt.fits', 'odnv02ahq_flt.fits'),
                   ('odnv02ahq_x2d.fits', 'odnv02ahq_x2d.fits')]
        self.compare_outputs(outputs)
