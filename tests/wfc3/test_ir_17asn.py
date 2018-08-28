import subprocess

from ..helpers import BaseWFC3


class TestIR17ASN(BaseWFC3):
    """Tests for WFC3/IR."""
    detector = 'ir'

    def test_ir_17asn(self):
        """Ported from ``calwf3_ir_17``."""
        asn_file = 'iacr61010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('iacr61010_spt.fits', 'iacr61010_spt_ref.fits')
        flist = ['iacr61ivq', 'iacr61iwq', 'iacr61iyq', 'iacr61izq']

        for rn in flist:
            outputs = [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                       (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
