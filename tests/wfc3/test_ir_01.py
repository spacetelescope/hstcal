import subprocess

from ..helpers import BaseWFC3


class TestIR01(BaseWFC3):
    """Tests for WFC3/IR."""
    detector = 'ir'

    def test_ir_01(self):
        """Ported from ``calwf3_ir_01``."""
        asn_file = 'iabg31030_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('iabg31031_spt.fits', 'iabg31031_spt_ref.fits')
        flist = ['iabg31f2q', 'iabg31f1q', 'iabg31f0q',
                 'iabg31ezq', 'iabg31exq', 'iabg31ewq', 'iabg31euq',
                 'iabg31evq', 'iabg31esq', 'iabg31etq']
        outputs = [('iabg31031_crj.fits', 'iabg31031_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                        (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
