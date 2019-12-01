import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS19Single(BaseWFC3):
    """
    Test pos UVIS2 Messier-82 data
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

    # Ported from ``calwf3_uv_19``.
    @pytest.mark.parametrize(
        'rootname', ['ib6wr9axq'])
        """
        'rootname', ['ib6wr9axq', 
                     'ib6wr9azq', 
                     'ib6wr9baq', 
                     'ib6wr9bcq', 
                     'ib6wr9bpq', 
                     'ib6wr9brq'])
        """
    def test_uvis_19single(self, rootname):
        self._single_raw_calib(rootname)
