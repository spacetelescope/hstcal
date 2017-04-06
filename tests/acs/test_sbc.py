from __future__ import absolute_import, division, print_function

import subprocess

from ..helpers import BaseACS


class TestSingleFrame(BaseACS):
    """Process a single SBC dataset."""
    subdir = 'acs_sbc_single'

    def test_pre_sm4(self):
        """This was ``sbc_single1``."""
        raw_file = 'j9ic01vpq_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('j9ic01vpq_flt.fits', 'j9ic01vpq_flt_ref.fits')]
        self.compare_outputs(outputs)

    def test_post_sm4(self):
        """This was ``sbc_single2``."""
        raw_file = 'jbdf10ykq_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('jbdf10ykq_flt.fits', 'jbdf10ykq_flt_ref.fits')]
        self.compare_outputs(outputs)
