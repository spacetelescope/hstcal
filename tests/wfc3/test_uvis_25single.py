import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS25Single(BaseWFC3):
    """
    Test pos UVIS2 NGC7318B data
    """
    
    detector = 'uvis'

    def _single_raw_calib(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', raw_file, '-vt'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname),
                    '{}_flt_ref.fits'.format(rootname))]
        self.compare_outputs(outputs)

    # Ported from ``calwf3_uv_25``.
    @pytest.mark.parametrize(
        'rootname', ['iacr52vjq', 
                     'iacr52vlq', 
                     'iacr52voq', 
                     'iacr52vqq'])
    def test_uvis_25single(self, rootname):
        self._single_raw_calib(rootname)
