import subprocess
import pytest

from ..helpers import BaseACS


class TestDarktimeSingle(BaseACS):
    """
    The updated darktime correction (March 2020) is applied to data
    according to the following categories:
    Full-frame, flashed and unflashed, pre-SM4
    Full-frame, flashed and unflashed, post-SM4
    Subarrays, flashed and unflashed, pre-Cycle 24
    Subarrays, flashed and unflashed, post-Cycle 24
    """
    detector = 'wfc'
    ignore_keywords = ['filename', 'date', 'iraf-tlm', 'fitsdate',
                       'opus_ver', 'cal_ver', 'proctime', 'history',
                       'bitpix', 'naxis', 'extend', 'simple']

    def _single_raw_cte(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname),
                    '{}_flt_ref.fits'.format(rootname)),
                   ('{}_flc.fits'.format(rootname),
                    '{}_flc_ref.fits'.format(rootname))]
        self.compare_outputs(outputs, ignore_keywords_overwrite=TestDarktimeSingle.ignore_keywords)

    def _single_raw(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname),
                    '{}_flt_ref.fits'.format(rootname))]
        self.compare_outputs(outputs, ignore_keywords_overwrite=TestDarktimeSingle.ignore_keywords)

    # FULL-FRAME
    # NOTE: This is slow test due to PCTECORR=PERFORM
    # 
    # j8ux05d3q, pre-SM4, postflashed, full-frame
    # j59l54gjq, pre-SM4, unflashed, full-frame
    # jch001zdq, post-SM4, postflashed, full-frame
    # jdw002gdq, post-SM4, unflashed, full-frame
    @pytest.mark.skip(reason="CTE is too slow to run under nominal conditions.")
    @pytest.mark.bigdata
    @pytest.mark.slow
    @pytest.mark.parametrize(
        'rootname', ['j8ux05d3q', 'j59l54gjq', 'jch001zdq', 'jdw002gdq'])
    def test_fullframe_darktime_single(self, rootname):
        self._single_raw_cte(rootname)

    # SUBARRAYS
    # no data,  pre-Cycle24, postflashed, subarray
    # j8mv02r4q, pre-Cycle24, unflashed, subarray
    # jdtr02oqq, post-Cycle24, postflashed, subarray
    # jdd701acq, post-Cycle24, unflashed, subarray
    @pytest.mark.bigdata
    @pytest.mark.parametrize(
        'rootname', ['j8mv02r4q', 'jdtr02oqq','jdd701acq'])
    def test_subarray_darktime_single(self, rootname):
        self._single_raw(rootname)
