from __future__ import absolute_import, division, print_function

import os
import subprocess

from astropy.io import fits

from ..helpers import BaseACS, download_file_cgi, raw_from_asn


class TestFullframeAsn(BaseACS):
    """
    Process pre- and post-SM4 fullframe WFC datasets using CR-SPLIT=2 with all
    standard calibration steps turned on.
    """
    subdir = 'acs_wfc_asn'

    def test_fullframe_presm4(self):
        """This was ``wfc_asn1.``"""
        rootname = 'j6lq01010'
        asn_file = rootname + '_asn.fits'
        asn0_file = asn_file + '.orig'

        # Prepare input files.
        download_file_cgi(self.tree, self.input_loc, asn0_file,
                          timeout=self.timeout)
        os.rename(asn0_file, asn_file)
        for raw_file in raw_from_asn(asn_file):
            self.get_input_file(raw_file)

            # Disable PCTECORR for now until we can handle long
            # execution time without timeout from CI provider.
            with fits.open(raw_file, mode='update') as pf:
                pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', asn_file, '-v'])

        # Compare results
        # j6lq01naq_fl?.fits tested in wfc_single1
        #
        # TAKEN OUT FOR NOW:
        # ('j6lq01011_crc.fits', 'j6lq01011_crc_ref_gen2cte.fits')
        # ('j6lq01ndq_flc.fits', 'j6lq01ndq_flc_ref_gen2cte.fits')
        outputs = [('j6lq01011_crj.fits', 'j6lq01011_crj_ref.fits'),
                   ('j6lq01ndq_flt.fits', 'j6lq01ndq_flt_ref.fits')]
        self.compare_outputs(outputs)

    def test_fullframe_postsm4(self):
        """This was ``wfc_asn2.``"""
        rootname = 'jbdf08010'
        asn_file = rootname + '_asn.fits'
        asn0_file = asn_file + '.orig'

        # Prepare input files.
        download_file_cgi(self.tree, self.input_loc, asn0_file,
                          timeout=self.timeout)
        os.rename(asn0_file, asn_file)
        for raw_file in raw_from_asn(asn_file):
            self.get_input_file(raw_file)

            # Disable PCTECORR for now until we can handle long
            # execution time without timeout from CI provider.
            with fits.open(raw_file, mode='update') as pf:
                pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', asn_file, '-v'])
        #subprocess.call(['calacs.e', asn_file, '-v', '--ctegen', '2', '--pctetab',
        #        '/grp/hst/cdbs/jref/16k1747tj_cte.fits'])

        # Compare results
        # jbdf08ufq_fl?.fits tested in wfc_single1
        #
        # TAKEN OUT FOR NOW:
        # ('jbdf08011_crc.fits', 'jbdf08011_crc_ref_gen2cte.fits')
        # ('jbdf08uhq_flc.fits', 'jbdf08uhq_flc_ref_gen2cte.fits')
        outputs = [('jbdf08011_crj.fits', 'jbdf08011_crj_ref.fits'),
                   ('jbdf08uhq_flt.fits', 'jbdf08uhq_flt_ref.fits')]
        self.compare_outputs(outputs)
