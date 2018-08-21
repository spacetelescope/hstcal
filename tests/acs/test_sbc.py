"""Tests for ACS/SBC."""
import subprocess
import pytest

from ..helpers import BaseACS


class TestSingleFrame(BaseACS):
    """Process a single SBC dataset."""
    detector = 'sbc'

    # NOTE:
    # j9ic01vpq = pre-SM4, was sbc_single1
    # jbdf10ykq = post-SM4, was sbc_single2
    @pytest.mark.parametrize(
        'rootname',
        ['j9ic01vpq',
         'jbdf10ykq'])
    def test_single_frame(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname),
                    '{}_flt_ref.fits'.format(rootname))]
        self.compare_outputs(outputs)
