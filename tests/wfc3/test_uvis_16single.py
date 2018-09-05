import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS16Single(BaseWFC3):
    """
    Test pos UVIS2 Jupiter impact site data
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

    # Ported from ``calwf3_uv_16``.
    @pytest.mark.parametrize(
        'rootname', ['ibcz21doq', 
                     'ibcz21dpq', 
                     'ibcz21drq', 
                     'ibcz21dsq', 
                     'ibcz21dtq', 
                     'ibcz21duq', 
                     'ibcz21dvq'])
    def test_uvis_16single(self, rootname):
        self._single_raw_calib(rootname)
