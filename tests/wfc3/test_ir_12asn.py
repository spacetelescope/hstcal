import subprocess
import pytest

from ..helpers import BaseWFC3

#@pytest.mark.slow
class TestIR12ASN(BaseWFC3):
    """Tests for WFC3/IR."""
    detector = 'ir'

    def test_ir_12asn(self):
        """Ported from ``calwf3_ir_12``."""
        asn_file = 'ibh719010_asn.fits'

        # Prepare input file.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        # Excluded ('ibh719011_spt.fits', 'ibh719011_spt_ref.fits')
        flist=['ibh719gmq',
               'ibh719gnq',
               'ibh719goq',
               'ibh719gpq',
               'ibh719gqq',
               'ibh719grq',
               'ibh719gsq',
               'ibh719gtq',
               'ibh719guq',
               'ibh719gvq',
               'ibh719gwq',
               'ibh719gxq',
               'ibh719gzq',
               'ibh719h0q',
               'ibh719h1q',
               'ibh719h2q',
               'ibh719h3q',
               'ibh719h4q',
               'ibh719h5q',
               'ibh719h6q',
               'ibh719h7q',
               'ibh719h8q',
               'ibh719h9q',
               'ibh719i2q',
               'ibh719i4q',
               'ibh719i5q',
               'ibh719i6q',
               'ibh719i7q',
               'ibh719i8q',
               'ibh719i9q',
               'ibh719iaq',
               'ibh719ibq',
               'ibh719icq',
               'ibh719idq',
               'ibh719ieq',
               'ibh719ifq',
               'ibh719ihq',
               'ibh719iiq',
               'ibh719ijq',
               'ibh719ikq',
               'ibh719ilq',
               'ibh719imq',
               'ibh719inq',
               'ibh719ioq',
               'ibh719ipq',
               'ibh719iqq',
               'ibh719irq',
               'ibh719isq',
               'ibh719iuq',
               'ibh719ivq',
               'ibh719iwq',
               'ibh719ixq',
               'ibh719iyq',
               'ibh719izq',
               'ibh719j0q',
               'ibh719j1q',
               'ibh719j2q',
               'ibh719j3q',
               'ibh719j4q',
               'ibh719j5q',
               'ibh719j7q',
               'ibh719j8q',
               'ibh719j9q',
               'ibh719jaq',
               'ibh719jbq',
               'ibh719jcq',
               'ibh719jdq',
               'ibh719jeq',
               'ibh719jfq',
               'ibh719jgq',
               'ibh719jhq',
               'ibh719jiq',
               'ibh719jkq',
               'ibh719jlq',
               'ibh719jmq',
               'ibh719jnq',
               'ibh719joq',
               'ibh719jpq',
               'ibh719jqq',
               'ibh719jrq',
               'ibh719jsq',
               'ibh719jtq',
               'ibh719juq',
               'ibh719jvq',
               'ibh719jxq',
               'ibh719jyq',
               'ibh719jzq',
               'ibh719k0q',
               'ibh719k1q',
               'ibh719k2q',
               'ibh719k3q',
               'ibh719k4q',
               'ibh719k5q',
               'ibh719k6q',
               'ibh719k7q',
               'ibh719k8q',
               'ibh719k9q',
               'ibh719kaq']

        outputs = [('ibh719011_crj.fits', 'ibh719011_crj_ref.fits')]

        for rn in flist:
            outputs += [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                        (rn + '_ima.fits', rn + '_ima_ref.fits')]

        self.compare_outputs(outputs)
