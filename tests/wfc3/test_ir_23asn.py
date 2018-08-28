import subprocess

from ..helpers import BaseWFC3


class TestIR23ASN(BaseWFC3):
    """Tests for WFC3/IR."""
    detector = 'ir'

    def test_ir_23asn(self):
        """Ported from ``calwf3_ir_23``."""
        asn_file = 'iabg21060_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('iabg21061_spt.fits', 'iabg21061_spt_ref.fits')
        flist = ['iabg21arq', 'iabg21asq', 'iabg21atq', 'iabg21auq', 'iabg21avq']
        outputs = [('iabg21061_crj.fits', 'iabg21061_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                       (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
