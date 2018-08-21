import subprocess
import pytest

from ..helpers import BaseACS


class TestWFCSubarray(BaseACS):
    """
    Process a single WFC subarray with all standard calibration steps on.
    """
    detector = 'wfc'

    # NOTE: This has PCTECORR=OMIT.
    # NOTE:
    # j9j902b6q = was wfc_sub1
    # jb2t11seq = Post-SM4 data with overscan, was wfc_sub2
    # j8c103xaq = Pre-SM4, was wfc_sub3
    # jb2t11se2 = Post-SM4 with overscan and FLSHCORR, was wfc_sub4
    @pytest.mark.parametrize(
        'rootname', ['j9j902b6q', 'jb2t11seq', 'j8c103xaq', 'jb2t11se2'])
    def test_subarray(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname),
                    '{}_flt_ref.fits'.format(rootname))]
        self.compare_outputs(outputs)

    # NOTE: This is slow test due to PCTECORR=PERFORM
    @pytest.mark.slow
    def test_subarray_pctecorr(self):
        """
        This was ``wfc_sub5``, single 2K subarray processed with PCTECORR.
        """
        raw_file = 'jb5s01fnq_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('jb5s01fnq_flt.fits', 'jb5s01fnq_flt_ref.fits'),
                   ('jb5s01fnq_flc.fits', 'jb5s01fnq_flc_ref_gen2cte.fits')]
        self.compare_outputs(outputs)
