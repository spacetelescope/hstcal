import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS17Single(BaseWFC3):
    """
    Test pos UVIS2 LSPM data
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

    # Ported from ``calwf3_uv_17``.
    @pytest.mark.parametrize(
        'rootname', ['ib3w01cnq'])
        """
        'rootname', ['ib3w01cnq', 
                     'ib3w01coq', 
                     'ib3w01csq', 
                     'ib3w01ctq', 
                     'ib3w01cwq', 
                     'ib3w01cxq', 
                     'ib3w01d0q', 
                     'ib3w01d1q', 
                     'ib3w01d4q', 
                     'ib3w01d5q'])
        """
    def test_uvis_17single(self, rootname):
        self._single_raw_calib(rootname)
