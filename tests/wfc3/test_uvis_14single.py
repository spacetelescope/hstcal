import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS14Single(BaseWFC3):
    """
    Test pos UVIS2 GD-71
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

    # Ported from ``calwf3_uv_14``.
    @pytest.mark.parametrize(
        'rootname', ['ibbq01moq'])
        """
        'rootname', ['ibbq01moq', 
                     'ibbq01mpq', 
                     'ibbq01mqq', 
                     'ibbq01mrq', 
                     'ibbq01msq', 
                     'ibbq01mtq', 
                     'ibbq01muq', 
                     'ibbq01mvq', 
                     'ibbq01mwq', 
                     'ibbq01mxq', 
                     'ibbq01myq', 
                     'ibbq01n0q', 
                     'ibbq01n1q', 
                     'ibbq01n2q', 
                     'ibbq01n3q'])
        """
    def test_uvis_14single(self, rootname):
        self._single_raw_calib(rootname)
