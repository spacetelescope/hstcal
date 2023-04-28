import subprocess
import pytest

from ..helpers import BaseWFC3


@pytest.mark.xfail
class TestUVIS06Single(BaseWFC3):
    """
    Tests for WFC3/UVIS.
    Post UVIS2 
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

    # ibcz21dtq = Ported from ``calwf3_uv_06``.
    @pytest.mark.parametrize(
        'rootname', ['ibcz21dtq'])
    def test_uvis_06single(self, rootname):
        self._single_raw_calib(rootname)
