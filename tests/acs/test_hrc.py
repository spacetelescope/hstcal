from __future__ import absolute_import, division, print_function

import os
import subprocess

from ..helpers import BaseACS, download_file_cgi, raw_from_asn


class TestMosaicBox(BaseACS):
    """
    Process an HRC mosaic dataset using the standard HRC-MOSAIC-BOX pattern
    with CR-SPLIT=1 at each of the 4 dither positions.

    .. note:: This was ``hrc_asn1``.

    """
    subdir = 'acs_hrc_asn1'

    def test_4point_mosaic(self):
        rootname = 'j6m901020'
        asn_file = rootname + '_asn.fits'
        asn0_file = asn_file + '.orig'

        # Prepare input files.
        download_file_cgi(self.tree, self.input_loc, asn0_file,
                          timeout=self.timeout)
        os.rename(asn0_file, asn_file)
        for raw_file in raw_from_asn(asn_file):
            self.get_input_file(raw_file)

        # Run CALACS
        subprocess.call(['calacs.e', asn_file, '-v'])

        # Compare results
        outputs = [('j6m901bzq_flt.fits', 'j6m901bzq_flt_ref.fits'),
                   ('j6m901c3q_flt.fits', 'j6m901c3q_flt_ref.fits'),
                   ('j6m901d9q_flt.fits', 'j6m901d9q_flt_ref.fits'),
                   ('j6m901deq_flt.fits', 'j6m901deq_flt_ref.fits')]
        self.compare_outputs(outputs)
