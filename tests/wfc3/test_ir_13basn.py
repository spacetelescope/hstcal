import subprocess

from ..helpers import BaseWFC3


class TestIR13BASN(BaseWFC3):
    """Tests for WFC3/IR - Subarray dark associations of small sizes."""
    detector = 'ir'

    def test_ir_13basn(self):
        """Ported from ``calwf3_ir_13b``."""
        asn_file = 'iabg21020_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('iabg21021_spt.fits', 'iabg21021_spt_ref.fits')
        flist = ['iabg21a6q', 'iabg21a7q',
                 'iabg21a8q', 'iabg21a9q',
                 'iabg21aaq']
        outputs = [('iabg21021_crj.fits', 'iabg21021_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                        (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
