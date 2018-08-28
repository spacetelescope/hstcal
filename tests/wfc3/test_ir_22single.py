import subprocess
import pytest

from ..helpers import BaseWFC3


class TestIR22Single(BaseWFC3):
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

    # Ported from calwf3_ir_22 (all single images to be processed)
    @pytest.mark.parametrize(
        'rootname', ['ibh714a1q', 'ibh714a3q', 'ibh714a4q',
                     'ibh714a6q', 'ibh714a7q', 'ibh714a9q'])
    def test_ir_22single(self, rootname):
        self._single_raw_calib(rootname)
