import subprocess
import pytest

from ..helpers import BaseWFC3


class TestIR13Single(BaseWFC3):
    """Tests for WFC3/IR."""
    detector = 'ir'

    def _single_raw_calib(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', raw_file, '-v'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname),
                    '{}_flt_ref.fits'.format(rootname)),
                   ('{}_ima.fits'.format(rootname),
                    '{}_ima_ref.fits'.format(rootname))]
        self.compare_outputs(outputs)

    # iabg21awq = Ported from calwf3_ir_13f
    # iabg21axq = Ported from calwf3_ir_13g
    @pytest.mark.parametrize(
        'rootname', ['iabg21awq', 'iabg21axq'])
    def test_ir_13single(self, rootname):
        self._single_raw_calib(rootname)
