import subprocess

from ..helpers import BaseWFC3


class TestIR13AASN(BaseWFC3):
    """Tests for WFC3/IR - Subarray dark associations of small sizes."""
    detector = 'ir'

    def test_ir_13aasn(self):
        """Ported from ``calwf3_ir_13a``."""
        asn_file = 'iabg21010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('iabg21011_spt.fits', 'iabg21011_spt_ref.fits')
        flist = ['iabg21a1q', 'iabg21a2q',
                 'iabg21a3q', 'iabg21a4q',
                 'iabg21a5q']
        outputs = [('iabg21011_crj.fits', 'iabg21011_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                        (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
