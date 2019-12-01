import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS30Single(BaseWFC3):
    """
    Test pos UVIS2 Tungsten data
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

    # Ported from ``calwf3_uv_30``.
    @pytest.mark.parametrize(
        'rootname', ['iaao06nfq'])
        """
        'rootname', ['iaao06nfq', 
                     'iaao06ngq', 
                     'iacs02toq', 
                     'iacs02trq'])
        """
    def test_uvis_30single(self, rootname):
        self._single_raw_calib(rootname)
