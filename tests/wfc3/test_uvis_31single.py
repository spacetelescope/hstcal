import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS31Single(BaseWFC3):
    """
    Test pos UVIS2 WR14 data
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

    # Ported from ``calwf3_uv_31``.
    @pytest.mark.parametrize(
        'rootname', ['ibbr02buq'])
        # 'rootname', ['ibbr02buq', 
        #             'ibbr02bvq', 
        #             'ibbr02bwq', 
        #             'ibbr02bxq', 
        #             'ibbr02byq', 
        #             'ibbr02bzq', 
        #             'ibbr02c0q', 
        #             'ibbr02c1q', 
        #             'ibbr02c2q', 
        #             'ibbr02c3q', 
        #             'ibbr02c5q'])
    def test_uvis_31single(self, rootname):
        self._single_raw_calib(rootname)
