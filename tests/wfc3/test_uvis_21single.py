import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS21Single(BaseWFC3):
    """
    Test pos UVIS2 Neptune data
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

    # Ported from ``calwf3_uv_21``.
    @pytest.mark.parametrize(
        'rootname', ['ib2s04jvq'])
        # 'rootname', ['ib2s04jvq',
        #             'ib2s04jwq',
        #             'ib2s04jxq',
        #             'ib2s04jyq',
        #             'ib2s04jzq',
        #             'ib2s04k1q'])
    def test_uvis_21single(self, rootname):
        self._single_raw_calib(rootname)
