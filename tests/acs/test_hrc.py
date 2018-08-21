"""Tests for ACS/HRC."""

import subprocess

from ..helpers import BaseACS


class TestMosaicBox(BaseACS):
    """
    Process an HRC mosaic dataset using the standard HRC-MOSAIC-BOX pattern
    with CR-SPLIT=1 at each of the 4 dither positions.

    .. note:: This was ``hrc_asn1``.

    """
    detector = 'hrc'

    def test_4point_mosaic(self):
        rootname = 'j6m901020'
        asn_file = rootname + '_asn.fits'

        # Prepare input files.
        self.get_input_file(asn_file)

        # Run CALACS
        subprocess.call(['calacs.e', asn_file, '-v'])

        # Compare results
        outputs = [('j6m901bzq_flt.fits', 'j6m901bzq_flt_ref.fits'),
                   ('j6m901c3q_flt.fits', 'j6m901c3q_flt_ref.fits'),
                   ('j6m901d9q_flt.fits', 'j6m901d9q_flt_ref.fits'),
                   ('j6m901deq_flt.fits', 'j6m901deq_flt_ref.fits')]
        self.compare_outputs(outputs)
