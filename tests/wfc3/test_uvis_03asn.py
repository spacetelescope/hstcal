import subprocess

from ..helpers import BaseWFC3


class TestUVIS03ASN(BaseWFC3):
    """
    Tests for WFC3/UVIS.
    """

    detector = 'uvis'

    def test_uvis_03asn(self):
        """Ported from ``calwf3_uvis_03``."""
        asn_file = 'ibbq01010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-vt'])

        # Compare results
        # Excluded ('ibbq01010_spt.fits', 'ibbq01010_spt_ref.fits')
        flist = ['ibbq01msq', 'ibbq01mrq', 'ibbq01mqq']

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits')]

        self.compare_outputs(outputs)
