"""Tests for ACS/HRC."""

import subprocess
import pytest

from ..helpers import BaseACS


class TestMosaic(BaseACS):
    """
    Process HRC mosaic datasets.
    """
    detector = 'hrc'

    def test_4point_mosaic(self):
        """
        Process an HRC mosaic dataset using the standard HRC-MOSAIC-BOX 
        pattern with CR-SPLIT=1 at each of the 4 dither positions.

        .. note:: This was ``hrc_asn1``.

        """
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

    @pytest.mark.xfail
    def test_3point_mosaic(self):
        """ 
        Process an HRC mossaic dataset using 3 dither positions with
        NUM_EXP=2 (RPT-OBS) at each of the dither positions and only
        DQICORR, BLEVCORR, BIASCORR, DARKCORR, FLATCORR, PHOTCORR,
        RPTCORR, and DITHCORR are set to PERFORM.
   
        .. note:: This was ``hrc_asn2``.
        """

        rootname = 'j8cd02011'
        asn_file = rootname + '_asn.fits'

        # Prepare input files.
        self.get_input_file(asn_file)

        # Run CALACS
        subprocess.call(['calacs.e', asn_file, '-v'])

        # Compare results
        outputs = [('j8cd02011_crj.fits', 'j8cd02011_sfl_ref.fits'),
                   ('j8cd020a1_crj.fits', 'j8cd020a1_sfl_ref.fits'),
                   ('j8cd020k1_crj.fits', 'j8cd020k1_sfl_ref.fits')]
        self.compare_outputs(outputs)
