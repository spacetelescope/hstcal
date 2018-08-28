import subprocess
import pytest

from ..helpers import BaseWFC3


class TestIR19Single(BaseWFC3):
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

    # Ported from calwf3_ir_19 (all single images to be processed)
    @pytest.mark.parametrize(
        'rootname', ['iabg31eeq', 'iabg31egq', 'iabg31ehq', 'iabg31eiq', 
                 'iabg31ejq', 'iabg31f3q', 'iabg31f4q', 'iabg31fdq', 
                 'iabg31feq', 'iabg31ffq', 'iabg31fgq', 'iabg31fhq', 
                 'iabg31fiq', 'ibm901gpq', 'ibm901gqq', 'iabg31efq'])
    def test_ir_19single(self, rootname):
        self._single_raw_calib(rootname)
