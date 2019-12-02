import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS12Single(BaseWFC3):
    """
    Test pos UVIS2 BIAS data
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

    # Ported from ``calwf3_uv_12``.
    @pytest.mark.parametrize(
        'rootname', ['iaao09l0q'])
        # 'rootname', ['iaao09l0q', 
        #             'iaao09l1q', 
        #             'iaao11odq', 
        #             'iaao11oeq', 
        #             'ibbq01n4q', 
        #             'ibbq01n5q', 
        #             'iblk57bzq', 
        #             'iblk57c0q', 
        #             'iblk57c3q'])
    def test_uvis_12single(self, rootname):
        self._single_raw_calib(rootname)
