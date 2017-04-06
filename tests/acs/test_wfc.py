from __future__ import absolute_import, division, print_function

import os
import subprocess

from ..helpers import BaseACS, download_file_cgi, raw_from_asn


class TestFullframePreSM4(BaseACS):
    """
    Process a pre-SM4 fullframe WFC dataset using CR-SPLIT=2 with all
    standard calibration steps turned on.

    .. note:: This was ``wfc_asn1``.

    """
    subdir = 'acs_wfc_asn1'

    def test_fullframe_presm4(self):
        rootname = 'j6lq01010'
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
        # j6lq01naq_fl?.fits tested in wfc_single1
        outputs = [('j6lq01011_crj.fits', 'j6lq01011_crj_ref.fits'),
                   ('j6lq01011_crc.fits', 'j6lq01011_crc_ref.fits'),
                   ('j6lq01ndq_flt.fits', 'j6lq01ndq_flt_ref.fits'),
                   ('j6lq01ndq_flc.fits', 'j6lq01ndq_flc_ref.fits')]
        self.compare_outputs(outputs)
