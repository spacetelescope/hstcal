import subprocess

from ..helpers import BaseWFC3


class TestUVIS08ASN(BaseWFC3):
    """
    Tests for WFC3/UVIS.
    Post UVIS2
    """

    detector = 'uvis'

    def test_uvis_08asn(self):
        """Ported from ``calwf3_uvis_08``."""
        asn_file = 'ibbsa5010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-vt'])

        # Compare results
        # Excluded ('ibbsa5010_spt.fits', 'ibbsa5010_spt_ref.fits')
        flist = ['ibbsa5ivq', 'ibbsa5iuq', 'ibbsa5itq', 'ibbsa5isq']

        outputs = []
        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits')]

        self.compare_outputs(outputs)
