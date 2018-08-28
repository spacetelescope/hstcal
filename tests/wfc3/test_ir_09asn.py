import subprocess

from ..helpers import BaseWFC3


class TestIR09ASN(BaseWFC3):
    """Tests for WFC3/IR."""
    detector = 'ir'

    def test_ir_09asn(self):
        """Ported from ``calwf3_ir_09``."""
        asn_file = 'ibvd11020_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        flist = ['ibvd11bwq', 'ibvd11bsq', 'ibvd11bgq', 'ibvd11bcq', 'ibvd11b8q',
                 'ibvd11baq', 'ibvd11beq', 'ibvd11bqq', 'ibvd11buq', 'ibvd11byq']

        for rn in flist:
            outputs = [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                       (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
