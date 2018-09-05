import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS29Single(BaseWFC3):
    """
    Test pos UVIS2 PN-G351.3+07.6 data
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

    # Ported from ``calwf3_uv_29``.
    @pytest.mark.parametrize(
        'rootname', ['ib1b0ko9q', 
                     'ib1b0koaq', 
                     'ib1b0kobq', 
                     'ib1b0kocq', 
                     'ib1b0kodq', 
                     'ib1b0koeq', 
                     'ib1b0kofq'])
    def test_uvis_29single(self, rootname):
        self._single_raw_calib(rootname)
