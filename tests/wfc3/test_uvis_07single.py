import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS07Single(BaseWFC3):
    """
    Tests for WFC3/UVIS.
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

    # iacs01t9q = Ported from ``calwf3_uv_07``.
    @pytest.mark.parametrize(
        'rootname', ['iacs01t9q'])
    def test_uvis_07single(self, rootname):
        self._single_raw_calib(rootname)
