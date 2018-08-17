import subprocess
import pytest

from ..helpers import BaseACS


class TestFullFrameSingle(BaseACS):
    """
    Process pre- and post-SM4 single fullframe WFC datasets
    using all standard calibration steps and the
    Generation 2 CTE correction.

    For pre-SM4, apply the original BLEVCORR algorithm,
    but no bias drift correction.

    For post-SM4, apply the 'new' BLEVCORR algorithm,
    which includes bias shift, cross talk, and destripe corrections.
    """
    detector = 'wfc'

    def _single_raw_calib(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname),
                    '{}_flt_ref.fits'.format(rootname)),
                   ('{}_flc.fits'.format(rootname),
                    '{}_flc_ref_gen2cte.fits'.format(rootname))]
        self.compare_outputs(outputs)

    # NOTE: This is slow test due to PCTECORR=PERFORM
    # NOTE:
    # j6lq01naq = pre-SM4, was wfc_single1
    # jbdf08ufq = post-SM4, was wfc_single2
    @pytest.mark.slow
    @pytest.mark.parametrize(
        'rootname', ['j6lq01naq', 'jbdf08ufq'])
    def test_fullframe_single(self, rootname):
        self._single_raw_calib(rootname)

    # NOTE: This is not marked slow to run one PCTECORR=PERFORM
    #       for a push event on GitHub. This alone takes about 8 mins.
    # NOTE:
    # jbdf08uf2 = post-SM4 with FLSHCORR, was wfc_single3
    def test_fullframe_single_postsm4_flshcorr(self):
        self._single_raw_calib('jbdf08uf2')
