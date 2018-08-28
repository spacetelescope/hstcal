import subprocess
import pytest

from ..helpers import BaseWFC3


class TestIR06Single(BaseWFC3):
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

    # ib3m55wxq = Ported from ``calwf3_ir_06``.
    @pytest.mark.parametrize(
        'rootname', ['ib3m55wxq'])
    def test_ir_06single(self, rootname):
        self._single_raw_calib(rootname)
