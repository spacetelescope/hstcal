import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS11Single(BaseWFC3):
    """
    Test pos UVIS2 30 Dor
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

    # Ported from ``calwf3_uv_11``.
    @pytest.mark.parametrize(
        'rootname', ['iacs01t9q'])
        # 'rootname', ['iacs01t9q', 'iacs01tbq', 'iacs02tiq', 'iacs02tlq'])
    def test_uvis_11single(self, rootname):
        self._single_raw_calib(rootname)
