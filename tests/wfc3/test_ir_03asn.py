import subprocess

from ..helpers import BaseWFC3


class TestIR03ASN(BaseWFC3):
    """Tests for WFC3/IR."""
    detector = 'ir'

    def test_ir_03asn(self):
        """Ported from ``calwf3_ir_03``."""
        asn_file = 'ib1f23010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('ib1f23010_spt.fits', 'ib1f23010_spt_ref.fits')
        flist = ['ib1f23qmq', 'ib1f23qlq', 'ib1f23qjq', 'ib1f23qhq']

        for rn in flist:
            outputs = [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                       (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
