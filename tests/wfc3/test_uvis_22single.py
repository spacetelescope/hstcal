import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS22Single(BaseWFC3):
    """
    Test pos UVIS2 NGC4639 data
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

    # Ported from ``calwf3_uv_22``.
    @pytest.mark.parametrize(
        'rootname', ['ib1f23qyq'])
    def test_uvis_22single(self, rootname):
        self._single_raw_calib(rootname)
