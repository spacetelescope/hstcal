import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS05ASN(BaseWFC3):
    """
    Tests for WFC3/UVIS.
    """

    detector = 'uvis'

    def test_uvis_05asn(self):
        """Ported from ``calwf3_uvis_05``."""
        asn_file = 'iacr51010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-vt'])

        # Compare results
        # Excluded ('iacr51010_spt.fits', 'iacr51010_spt_ref.fits')
        # Excluded ('iacr51012_spt.fits', 'iacr51012_spt_ref.fits')
        flist = ['iacr51omq', 'iacr51okq', 'iacr51ojq', 'iacr51ohq']

        outputs = [('iacr51011_crj.fits', 'iacr51011_crj_ref.fits'),
                   ('iacr51011_crc.fits', 'iacr51011_crc_ref.fits'),
                   ('iacr51012_crc.fits', 'iacr51012_crc_ref.fits'),
                   ('iacr51012_crj.fits', 'iacr51012_crj_ref.fits')]
                 

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                        (rn + '_flc.fits', rn + '_flc_ref.fits')]

        self.compare_outputs(outputs)
