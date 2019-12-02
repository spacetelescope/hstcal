import subprocess
import pytest

from ..helpers import BaseWFC3


class TestUVIS20Single(BaseWFC3):
    """
    Test pos UVIS2 Messier-83 data
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

    # Ported from ``calwf3_uv_20``.
    @pytest.mark.parametrize(
        'rootname', ['ib6w62toq'])
        # 'rootname', ['ib6w62toq',
        #             'ib6w62tvq',
        #             'ib6w62txq',
        #             'ib6w62tzq',
        #             'ib6w63b3q',
        #             'ib6w63bcq',
        #             'ib6w63bhq'])
    def test_uvis_20single(self, rootname):
        self._single_raw_calib(rootname)
