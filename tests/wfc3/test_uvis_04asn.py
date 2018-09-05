import subprocess

from ..helpers import BaseWFC3


class TestUVIS04ASN(BaseWFC3):
    """
    Tests for WFC3/UVIS.
    """

    detector = 'uvis'

    def test_uvis_04asn(self):
        """Ported from ``calwf3_uvis_04``."""
        asn_file = 'ib1f23040_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-vt'])

        # Compare results
        # Excluded ('ib1f23040_spt.fits', 'ib1f23040_spt_ref.fits')
        flist = ['ib1f23r0q', 'ib1f23qyq']

        for rn in flist:
            outputs = [(rn + '_flt.fits', rn + '_flt_ref.fits')]

        self.compare_outputs(outputs)
