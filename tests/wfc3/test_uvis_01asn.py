import subprocess

from ..helpers import BaseWFC3


class TestUVIS01ASN(BaseWFC3):
    """
    Tests for WFC3/UVIS.
    Post UVIS2 1995DW2
    """

    detector = 'uvis'

    def test_uvis_01asn(self):
        """Ported from ``calwf3_uvis_01``."""
        asn_file = 'iaao01010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-vt'])

        # Compare results
        # Excluded ('iaao01011_spt.fits', 'iaao01011_spt_ref.fits')
        flist = ['iaao01k9q', 'iaao01k8q']

        outputs = [('iaao01011_crj.fits', 'iaao01011_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                        (rn + '_flc.fits', rn + '_flc_ref.fits')]

        self.compare_outputs(outputs)
