from __future__ import absolute_import, division, print_function

import subprocess

from astropy.io import fits

from ..helpers import BaseACS


class TestFullframeSingle(BaseACS):
    """
    Process pre- and post-SM4 single fullframe WFC datasets.
    usage: pytest -r a -s --remote-data --basetemp=/processingDirectory test_wfc_single.py

    This class is a subclass of BaseCal->BaseACS and extends the functionality of these superclasses 
    to invoke very specific regression test cases.
    """

    subdir = 'acs_wfc_single'

    def test_fullframe_single_presm4(self):
        """
        Process a single WFC file of pre-SM4 data using all standard calibration steps.

        Apply the original BLEVCORR algorithm, but no bias drift correction, and use
        the Generation 2 CTE correction. (wfc_single1)
        """
        raw_file = 'j6lq01naq_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Disable PCTECORR for now until we can handle long
        # execution time without timeout from CI provider.
        with fits.open(raw_file, mode='update') as pf:
            pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])
        #subprocess.call(['calacs.e', raw_file, '-v', '--ctegen', '2',
        #                        '--pctetab', '/grp/hst/cdbs/jref/16k1747tj_cte.fits'])

        # Compare results
        outputs = [('j6lq01naq_flt.fits', 'j6lq01naq_flt_ref_gen2cte.fits')]
        #outputs = [('j6lq01naq_flt.fits', 'j6lq01naq_flt_ref_gen2cte.fits'),
        #           ('j6lq01naq_flc.fits', 'j6lq01naq_flc_ref_gen2cte.fits')]
        self.compare_outputs(outputs)

    def test_fullframe_single_postsm4(self):
        """
        Process a single WFC file of post-SM4 data using all standard calibration steps.

        Apply the ``new`` BLEVCORR algorithm which includes bias shift, cross talk, 
        and destripe corrections, and use the Generation 2 CTE correction. (wfc_single2)
        """
        raw_file = 'jbdf08ufq_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Disable PCTECORR for now until we can handle long
        # execution time without timeout from CI provider.
        with fits.open(raw_file, mode='update') as pf:
            pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])
        #subprocess.call(['calacs.e', raw_file, '-v', '--ctegen', '2',
        #                        '--pctetab', '/grp/hst/cdbs/jref/16k1747tj_cte.fits'])

        # Compare results
        outputs = [('jbdf08ufq_flt.fits', 'jbdf08ufq_flt_ref_gen2cte.fits')]
        #outputs = [('jbdf08ufq_flt.fits', 'jbdf08ufq_flt_ref_gen2cte.fits'),
        #           ('jbdf08ufq_flc.fits', 'jbdf08ufq_flc_ref_gen2cte.fits')]
        self.compare_outputs(outputs)

    def test_fullframe_single_postsm4_flshcorr(self):
        """
        Process a single WFC file of post-SM4 data using all standard calibration steps, 

        plus FLSHCORR.  Apply the ``new`` BLEVCORR algorithm which includes bias shift, 
        cross talk, and destripe corrections, and use the Generation 2 CTE correction. (wfc_single3)
        """
        raw_file = 'jbdf08uf2_raw.fits'

        # Prepare input file.
        self.get_input_file(raw_file)

        # Disable PCTECORR for now until we can handle long
        # execution time without timeout from CI provider.
        with fits.open(raw_file, mode='update') as pf:
            pf[0].header['PCTECORR'] = 'OMIT'

        # Run CALACS
        subprocess.call(['calacs.e', raw_file, '-v'])
        #subprocess.call(['calacs.e', raw_file, '-v', '--ctegen', '2',
        #                        '--pctetab', '/grp/hst/cdbs/jref/16k1747tj_cte.fits'])

        # Compare results
        outputs = [('jbdf08uf2_flt.fits', 'jbdf08uf2_flt_ref_gen2cte.fits')]
        #outputs = [('jbdf08uf2_flt.fits', 'jbdf08uf2_flt_ref_gen2cte.fits'),
        #           ('jbdf08uf2_flc.fits', 'jbdf08uf2_flc_ref_gen2cte.fits')]
        self.compare_outputs(outputs)
