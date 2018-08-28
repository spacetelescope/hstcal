import subprocess

from ..helpers import BaseWFC3


class TestIR13EASN(BaseWFC3):
    """Tests for WFC3/IR - Subarray dark associations of small sizes."""
    detector = 'ir'

    def test_ir_13easn(self):
        """Ported from ``calwf3_ir_13e``."""
        asn_file = 'iabg21050_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('iabg21051_spt.fits', 'iabg21051_spt_ref.fits')
        flist = ['iabg21amq', 'iabg21anq',
                 'iabg21aoq', 'iabg21apq',
                 'iabg21aqq']
        outputs = [('iabg21051_crj.fits', 'iabg21051_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                        (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
