import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS09Single(BaseWFC3):
    """
    Check postflash correction for WFC3 UVIS (Post UVIS2) full and subarrays.
    """
    
    detector = 'uvis'

    def _single_raw_calib(self, rootname):
        raw_file = '{}_raw.fits'.format(rootname)

        # Prepare input file.
        self.get_input_file(raw_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', raw_file, '-vt'])

        # Compare results
        outputs = [('{}_flt.fits'.format(rootname),
                    '{}_flt_ref.fits'.format(rootname))]
        self.compare_outputs(outputs)

    # Ported from ``calwf3_uv_09``.
    # ic5p02e6q: subarray which straddles the overscan region in the postflash reference file
    #   ccdamp=A, straddle=1, offset=(0,0), ampx,ampy=(512,512), x0,y0=(1610,1557)
    # ibly01g2q: full frame image with postflash
    #   ccdamp=A, straddle=0, offset=(0,0), ampx,ampy=(2070,2050), x0,y0=(2,19)
    # ic5p02e0q: subarray only in Amp A
    #   ccdamp=A, straddle=0, offset=(0,0), ampx,ampy=(512,512), x0,y0=(1086,1381)
    # Sort out reference file issue for ic5p02eoq
    #   'rootname', ['ic5p02e6q', 'ibly01g2q', 'ic5p02e0q'])
    @pytest.mark.parametrize(
        'rootname', ['ic5p02e6q', 'ibly01g2q'])
    def test_uvis_09single(self, rootname):
        self._single_raw_calib(rootname)
