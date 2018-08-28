import subprocess
import pytest

from ..helpers import BaseWFC3


class TestIR14Single(BaseWFC3):
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

    # ibbt02atq = Ported from calwf3_ir_14a
    # ibbt02auq = Ported from calwf3_ir_14b
    # ibbt02avq = Ported from calwf3_ir_14c
    # ibbt02awq = Ported from calwf3_ir_14d
    # ibbt02ayq = Ported from calwf3_ir_14e
    # ibbt02azq = Ported from calwf3_ir_14f
    @pytest.mark.parametrize(
        'rootname', ['ibbt02atq', 'ibbt02auq', 'ibbt02avq',
                     'ibbt02awq', 'ibbt02ayq', 'ibbt02azq'])
    def test_ir_14single(self, rootname):
        self._single_raw_calib(rootname)
