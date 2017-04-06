from __future__ import absolute_import, division, print_function

import os
import subprocess

from ..helpers import BaseWFC3IR, download_file_cgi, raw_from_asn


class TestIR01(BaseWFC3IR):
    """Ported from ``calwf3_ir_01``."""
    subdir = 'wfc3_ir_01'

    def test_ir_01(self):
        rootname = 'iabg31030'
        asn_file = rootname + '_asn.fits'
        asn0_file = rootname + '_asn_orig.fits'

        # Prepare input files.
        download_file_cgi(self.tree, self.input_loc, asn0_file,
                          timeout=self.timeout)
        os.rename(asn0_file, asn_file)
        for raw_file in raw_from_asn(asn_file):
            self.get_input_file(raw_file)

        # Run CALWF3
        subprocess.call(['calwf3.e', asn_file, '-v'])

        # Compare results
        flist = ['iabg31f2q', 'iabg31f1q', 'iabg31f0q',
                 'iabg31ezq', 'iabg31exq', 'iabg31ewq', 'iabg31euq',
                 'iabg31evq', 'iabg31esq', 'iabg31etq']
        outputs = [('iabg31031_crj.fits', 'iabg31031_crj_ref.fits')]

        # TODO: Ask Megan why this file is missing.
        #outputs.append(('iabg31031_spt.fits', 'iabg31031_spt_ref.fits'))

        for rn in flist:
            subouts = [(rn + '_flt.fits', rn + '_flt_ref.fits'),
                       (rn + '_ima.fits', rn + '_ima_ref.fits')]
            outputs += subouts

        self.compare_outputs(outputs)
