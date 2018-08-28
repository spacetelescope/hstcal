import subprocess
import pytest

from ..helpers import BaseWFC3


class TestIR10Single(BaseWFC3):
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

    # ib2k03dlq = Ported from calwf3_ir_10
    @pytest.mark.parametrize(
        'rootname', ['ib2k03dlq'])
    def test_ir_10single(self, rootname):
        self._single_raw_calib(rootname)
