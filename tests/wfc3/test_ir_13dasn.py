import subprocess

from ..helpers import BaseWFC3


class TestIR13DASN(BaseWFC3):
    """Tests for WFC3/IR - Subarray dark associations of small sizes."""
    detector = 'ir'

    def test_ir_13dasn(self):
        """Ported from ``calwf3_ir_13d``."""
        asn_file = 'iabg21040_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('iabg21041_spt.fits', 'iabg21041_spt_ref.fits')
        flist = ['iabg21agq', 'iabg21aiq',
                 'iabg21ajq', 'iabg21akq',
                 'iabg21alq']
        outputs = [('iabg21041_crj.fits', 'iabg21041_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                        (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
