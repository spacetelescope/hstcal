import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS32Single(BaseWFC3):
    """
    Test pos UVIS2 subarray data with CTE correction
    """
    
    detector = 'uvis'

    def _single_raw_calib(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', raw_file, '-vts'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname), '{}_flt_ref.fits'.format(rootname)),
                   ('{}_flc.fits'.format(rootname), '{}_flc_ref.fits'.format(rootname)),
                   ('{}_rac_tmp.fits'.format(rootname), '{}_rac_tmp_ref.fits'.format(rootname))]
        self.compare_outputs(outputs)

    # Ported from ``calwf3_uv_32``.
    @pytest.mark.parametrize(
        'rootname', ['ib3805v0q'])
        # 'rootname', ['ib3805v0q', 
        #             'ib2kabmaq', 
        #             'ib3503wwq', 
        #             'ibde04msq', 
        #             'icoc14hcq'])
    def test_uvis_32single(self, rootname):
        self._single_raw_calib(rootname)
