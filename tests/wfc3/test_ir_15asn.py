"""Tests for WFC3/IR ASN files."""
import subprocess
import pytest

from ..helpers import BaseWFC3


class TestIR15ASN(BaseWFC3):
    """Tests for WFC3/IR."""
    detector = 'ir'

    # NOTE:
    # ib1f23010 = Ported from calwf3_ir_15a
    # ib1f23020 = Ported from calwf3_ir_02 (15b internally)
    # NOTE:
    # https://github.com/spacetelescope/hstcal/pull/71#issuecomment-414132998
    @pytest.mark.parametrize(
        ('rootname', 'outroots'),
        [('ib1f23010', ['ib1f23qhq', 'ib1f23qjq', 'ib1f23qlq', 'ib1f23qmq']),
         ('ib1f23020', ['ib1f23qoq', 'ib1f23qwq'])])
    def test_fullframe(self, rootname, outroots):
        asn_file = rootname + '_asn.fits'

        # Prepare input files.
        self.get_input_file(asn_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results.
        # The first outroot is the output from whole ASN,
        # the rest are individual members.
        # Excluded ('iabg23010_spt.fits', 'iabg23010_spt_ref.fits')
        # Excluded ('iabg23020_spt.fits', 'iabg23020_spt_ref.fits')
        #outputs = [('{}_spt.fits'.format(outroots[0]),
        #            '{}_spt_ref.fits'.format(outroots[0]))]
        for outroot in outroots[1:]:
            #outputs += [('{}_ima.fits'.format(outroot),
            outputs  = [('{}_ima.fits'.format(outroot),
                         '{}_ima_ref.fits'.format(outroot)),
                        ('{}_flt.fits'.format(outroot),
                         '{}_flt_ref.fits'.format(outroot))]
        self.compare_outputs(outputs)
