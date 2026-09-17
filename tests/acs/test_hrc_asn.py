import subprocess
import pytest

from ..helpers import BaseACS


class TestSpotTabASN(BaseACS):
    """
    Process single HRC ASN using all standard calibration steps.

    This test covers the route where ACS2D is run after ACSREJ
    and flatfield correction applies SPOTTAB.
    """
    detector = 'hrc'

    ignore_keywords = ['filename', 'date', 'iraf-tlm', 'fitsdate',
                       'opus_ver', 'cal_ver', 'proctime', 'history',
                       'bitpix', 'naxis', 'extend', 'simple']

    def test_asn_01(self):
        #   j8cw08nwq_raw.fits  j8cw08nyq_raw.fits
        asn_file ="j8cw08070_asn.fits"

        # Prepare input files.
        self.get_input_file(asn_file)

        # Run CALACS
        subprocess.call(['calacs.e', asn_file, '-v'])

        # Compare results.
        # The first outroot is the output from whole ASN,
        # the rest are individual members.
        outputs = [("j8cw08071_crj.fits", "j8cw08071_crj_ref.fits"),
                   ("j8cw08nwq_flt.fits", "j8cw08nwq_flt_ref.fits"),
                   ("j8cw08nyq_flt.fits", "j8cw08nyq_flt_ref.fits")]
        self.compare_outputs(outputs, ignore_keywords_overwrite=TestSpotTabASN.ignore_keywords)
