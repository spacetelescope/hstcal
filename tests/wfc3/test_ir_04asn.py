import subprocess

from ..helpers import BaseWFC3


class TestIR04ASN(BaseWFC3):
    """Tests for WFC3/IR."""
    detector = 'ir'

    def test_ir_04asn(self):
        """Ported from ``calwf3_ir_04``."""
        asn_file = 'ib8t01010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('ib8t01010_spt.fits', 'ib8t01010_spt_ref.fits')
        flist = ['ib8t01hwq', 'ib8t01htq']

        for rn in flist:
            outputs = [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                       (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
