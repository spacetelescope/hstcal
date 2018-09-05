import subprocess

from ..helpers import BaseWFC3


class TestUVIS02ASN(BaseWFC3):
    """
    Tests for WFC3/UVIS.
    Post UVIS2 1995DW2
    """

    detector = 'uvis'

    def test_uvis_02asn(self):
        """Ported from ``calwf3_uvis_02``."""
        asn_file = 'iaal01010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-vt'])

        # Compare results
        # Excluded ('iaal01011_spt.fits', 'iaal01011_spt_ref.fits')
        flist = ['iaal01hyq', 'iaal01hxq']

        outputs = [('iaal01011_crj.fits', 'iaal01011_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits')]

        self.compare_outputs(outputs)
