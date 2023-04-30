import subprocess
import pytest

from ..helpers import BaseWFC3


@pytest.mark.xfail
class TestUVIS33Single(BaseWFC3):
    """
    Test pos UVIS2 subarray data with on CTE correction
    """
    
    detector = 'uvis'

    def _single_raw_calib(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', raw_file, '-vts'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname), '{}_flt_ref.fits'.format(rootname))]
        self.compare_outputs(outputs)

    # Ported from ``calwf3_uv_33``.
    @pytest.mark.parametrize(
        'rootname', ['ib7i70b7q'])
    def test_uvis_33single(self, rootname):
        self._single_raw_calib(rootname)
