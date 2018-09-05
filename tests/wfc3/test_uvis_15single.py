import subprocess
import pytest

from ..helpers import BaseWFC3

@pytest.mark.slow
class TestUVIS15Single(BaseWFC3):
    """
    Test pos UVIS2 GRW+75D5824
    """
    
    detector = 'uvis'

    def _single_raw_calib(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', raw_file, '-vt'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname),
                    '{}_flt_ref.fits'.format(rootname))]
        self.compare_outputs(outputs)

    # Ported from ``calwf3_uv_15``.
    @pytest.mark.parametrize(
        'rootname', ['iaau01k3q',
                     'iaau01k4q',
                     'iaau01k5q', 
                     'iaau01k6q',
                     'iaau01kcq',
                     'iaau01kdq',
                     'iaau01keq',
                     'iaau01kfq',
                     'iaau01kgq',
                     'iaau01khq',
                     'iaau01kiq',
                     'iaau01kjq',
                     'iaau01kkq',
                     'iaau01kmq',
                     'iaau01knq',
                     'iaau01koq',
                     'iaau01kpq',
                     'iaau01kqq',
                     'iaau01krq',
                     'iaau01kxq',
                     'iaau01kyq',
                     'iaau01kzq',
                     'iaau01l0q',
                     'iaau01l1q',
                     'iaau01l2q',
                     'iaau01l3q',
                     'iaau01l4q',
                     'iaau01l5q',
                     'iaau01l9q',
                     'iaau01laq',
                     'iaau01lbq',
                     'iaau01lcq',
                     'iaaua1ldq',
                     'iaaua1leq',
                     'iaaua1lfq',
                     'iaaua1lgq',
                     'iaaua1lhq',
                     'iaaua1liq',
                     'iaaua1lkq',
                     'iaaua1llq',
                     'iaaua1lsq',
                     'iaaua1ltq',
                     'iaaua1luq',
                     'iaaua1lvq',
                     'iaaua1lwq',
                     'iaaua1lxq',
                     'iaaua1lzq',
                     'iaaua1m0q',
                     'iaaua1m1q',
                     'iaaua1m2q',
                     'iaaua1m3q',
                     'iaaua1mtq',
                     'iaaua1muq',
                     'iaaua1mvq',         
                     'iaaua1mwq',
                     'iaaua1mxq',
                     'iaaua1myq',
                     'iaaua1mzq',
                     'iaaua1n0q',
                     'iaaua1n1q',
                     'iaaua1n2q',
                     'iaaua1n3q',
                     'iaaua1n4q',
                     'iaaua1n5q',
                     'ibbs05i6q',
                     'ibbs05i7q',
                     'ibbs05i8q',
                     'ibbs05i9q',
                     'ibbs05iaq',
                     'ibbs05ibq',
                     'ibbs05icq',
                     'ibbs05idq',
                     'ibbs05ieq',
                     'ibbs05ifq',
                     'ibbs05igq',
                     'ibbs05ihq',
                     'ibbs05iiq',
                     'ibbs05ijq',
                     'ibbs05ikq',
                     'ibbs05imq',
                     'ibbs05inq',
                     'ibbs05ioq',
                     'ibbsa5isq',
                     'ibbsa5iuq',
                     'ibbsa5ivq',
                     'ibbsa5iwq',
                     'ibbsa5ixq',
                     'ibbsa5iyq',
                     'ibbsa5izq',
                     'ibbsa5j0q',
                     'ibbsa5j1q',
                     'ibbsa5j2q',
                     'ibbsa5j3q',
                     'ibbsa5j4q',
                     'ibbsa5j5q',
                     'ibbsa5j6q',
                     'ibbsa5j7q',
                     'ibbsa5j8q',
                     'ibbsa5j9q'])

    def test_uvis_15single(self, rootname):
        self._single_raw_calib(rootname)
