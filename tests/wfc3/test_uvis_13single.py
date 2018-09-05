import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS13Single(BaseWFC3):
    """
    Test pos UVIS2 DARK images
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

    # Ported from ``calwf3_uv_13``.
    @pytest.mark.parametrize(
        'rootname', ['iaao09l2q', 'iaao09l3q', 'iaao11ofq', 'iaao11ogq', 'iblk57c1q'])
    def test_uvis_13single(self, rootname):
        self._single_raw_calib(rootname)
