"""Tests for ACS/WFC ASN files."""
import subprocess
import pytest

from ..helpers import BaseACS


# NOTE: This is slow test due to PCTECORR=PERFORM
@pytest.mark.slow
class TestFullFrameASN(BaseACS):
    """
    Process pre- and post-SM4 fullframe WFC datasets using CR-SPLIT=2 with all
    standard calibration steps turned on.
    """
    detector = 'wfc'

    # NOTE:
    # j6lq01010 = pre-SM4, was wfc_asn1
    # jbdf08010 = post-SM4, was wfc_asn2
    @pytest.mark.parametrize(
        ('rootname', 'outroots'),
        [('j6lq01010', ['j6lq01011', 'j6lq01naq', 'j6lq01ndq']),
         ('jbdf08010', ['jbdf08011', 'jbdf08ufq', 'jbdf08uhq'])])
    def test_fullframe(self, rootname, outroots):
        asn_file = rootname + '_asn.fits'

        # Prepare input files.
        self.get_input_file(asn_file)

        # Run CALACS
        subprocess.call(['calacs.e', asn_file, '-v'])

        # Compare results.
        # The first outroot is the output from whole ASN,
        # the rest are individual members.
        outputs = [('{}_crj.fits'.format(outroots[0]),
                    '{}_crj_ref.fits'.format(outroots[0])),
                   ('{}_crc.fits'.format(outroots[0]),
                    '{}_crc_ref_gen2cte.fits'.format(outroots[0]))]
        for outroot in outroots[1:]:
            outputs += [('{}_flt.fits'.format(outroot),
                         '{}_flt_ref.fits'.format(outroot)),
                        ('{}_flc.fits'.format(outroot),
                         '{}_flc_ref_gen2cte.fits'.format(outroot))]
        self.compare_outputs(outputs)
