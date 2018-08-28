import subprocess

from ..helpers import BaseWFC3


class TestIR13CASN(BaseWFC3):
    """Tests for WFC3/IR - Subarray dark associations of small sizes."""
    detector = 'ir'

    def test_ir_13casn(self):
        """Ported from ``calwf3_ir_13c``."""
        asn_file = 'iabg21030_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('iabg21031_spt.fits', 'iabg21031_spt_ref.fits')
        flist = ['iabg21abq', 'iabg21acq',
                 'iabg21adq', 'iabg21aeq',
                 'iabg21afq']
        outputs = [('iabg21031_crj.fits', 'iabg21031_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                        (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
