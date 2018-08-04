from __future__ import absolute_import, division, print_function

import subprocess

from astropy.io import fits

from ..helpers import BaseACS


class TestWFCSubarray(BaseACS):
    """Process a single WFC subarray with all standard calibration steps on."""
    subdir = 'acs_wfc_subarray'

    def test_subarray(self):
        """This was ``wfc_sub1, all standard calibration steps.``"""
        raw_file = 'j9j902b6q_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Disable PCTECORR for now until we can handle long
        # execution time without timeout from CI provider.
        with fits.open(raw_file, mode='update') as pf:
            pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('j9j902b6q_flt.fits', 'j9j902b6q_flt_ref.fits')]
        self.compare_outputs(outputs)


    def test_subarray_postsm4(self):
        """This was ``wfc_sub2, Post-SM4 data with overscan, all standard calibration steps.``"""
        raw_file = 'jb2t11seq_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Disable PCTECORR for now until we can handle long
        # execution time without timeout from CI provider.
        with fits.open(raw_file, mode='update') as pf:
            pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('jb2t11seq_flt.fits', 'jb2t11seq_flt_ref.fits')]
        self.compare_outputs(outputs)


    def test_subarray_presm4(self):
        """This was ``wfc_sub3, Pre-SM4, all standard calibration steps.``"""
        raw_file = 'j8c103xaq_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Disable PCTECORR for now until we can handle long
        # execution time without timeout from CI provider.
        with fits.open(raw_file, mode='update') as pf:
            pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('j8c103xaq_flt.fits', 'j8c103xaq_flt_ref.fits')]
        self.compare_outputs(outputs)


    def test_subarray_postsm4_flshcorr(self):
        """This was ``wfc_sub4, Post-SM4 with overscan and FLSHCORR, all standard calibration steps.``"""
        raw_file = 'jb2t11se2_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Disable PCTECORR for now until we can handle long
        # execution time without timeout from CI provider.
        with fits.open(raw_file, mode='update') as pf:
            pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])

        # Compare results
        outputs = [('jb2t11se2_flt.fits', 'jb2t11se2_flt_ref.fits')]
        self.compare_outputs(outputs)


    def test_subarray_pctecorr(self):
        """This was ``wfc_sub5, single 2K subarray processed with PCTECORR, all standard calibration steps.``"""
        raw_file = 'jb5s01fnq_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Disable PCTECORR for now until we can handle long
        # execution time without timeout from CI provider.
        with fits.open(raw_file, mode='update') as pf:
            pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])
        #subprocess.call(['calacs.e', raw_file, '-v', '--ctegen', '2',
        #                 '--pctetab', '/grp/hst/cdbs/jref/16k1747tj_cte.fits'])

        # Compare results
        outputs = [('jb5s01fnq_flt.fits', 'jb5s01fnq_flt_ref.fits'), 
                   ('jb5s01fnq_flc.fits', 'jb5s01fnq_flc_ref.fits')]
        self.compare_outputs(outputs)
