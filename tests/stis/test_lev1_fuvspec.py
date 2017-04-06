from __future__ import absolute_import, division, print_function

import subprocess

from ..helpers import BaseSTIS, download_file_cgi


class TestLev1FUVSpec(BaseSTIS):
    """
    Level 1 FUV spec.

    .. note:: This was ``calstis_lev1_FUVspec``.

    """
    subdir = 'stis_fuvspec_lev1'

    def test_lev1(self):
        rootname = 'o5cl02040'
        raw_file = rootname + '_raw.fits'
        wav_file = rootname + '_wav.fits'

        # Prepare input files.
        self.get_input_file(raw_file)
        download_file_cgi(self.tree, self.input_loc, wav_file,
                          timeout=self.timeout)

        # Run CALSTIS (equivalent to stistools.calstis.calstis)
        subprocess.call(['cs0.e', raw_file, '-v'])

        # Compare results
        outputs = [('o5cl02040_flt.fits', 'calstis_lev1_FUVspec_flt.fits'),
                   ('o5cl02040_x1d.fits', 'calstis_lev1_FUVspec_x1d.fits')]
        self.compare_outputs(outputs)
