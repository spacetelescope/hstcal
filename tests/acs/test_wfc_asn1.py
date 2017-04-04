from __future__ import absolute_import, division, print_function

import os
import shutil
import subprocess

import pytest
from astropy.io.fits import FITSDiff
from astropy.table import Table

from ..helpers import (HAS_CALXXX, remote_data, use_calacs,
                       download_file_cgi, download_ref_ftp)


@remote_data
@use_calacs
class TestFullframePreSM4(object):
    """
    Process a pre-SM4 fullframe WFC dataset using CR-SPLIT=2 with all
    standard calibration steps turned on.
    """
    prevdir = os.getcwd()  # Need to be here for teardown access
    prevjref = os.environ.get('jref')

    @pytest.fixture(autouse=True)
    def setup_class(self, tmpdir):
        """
        Run test in own dir so we can keep results separate from
        other tests.
        """
        self.tree = 'rt'  # Use dev for now
        self.input_loc = 'hstcal/acs/calacs_e'
        self.ref_loc = 'hstcal/acs/calacs_e/ref'

        p = tmpdir.mkdir('wfc_asn1').strpath
        os.chdir(p)

        # HSTCAL cannot open remote JREF
        os.environ['jref'] = '.'

    def teardown_class(self):
        """Change back to test directory."""
        os.chdir(self.prevdir)
        os.environ['jref'] = self.prevjref

    def test_fullframe_presm4(self):
        """Meat of the test."""
        ignore_keywords = ['filename', 'date', 'iraf-tlm', 'fitsdate',
                           'opus_ver', 'cal_ver', 'proctime']

        rootname = 'j6lq01010'
        asn_file = rootname + '_asn.fits'
        asn0_file = asn_file + '.orig'

        # j6lq01naq_flt.fits is tested in wfc_single1_noniraf
        # j6lq01naq_flc.fits is tested in wfc_single1_noniraf
        outputs = [('j6lq01011_crj.fits', 'j6lq01011_crj_ref.fits'),
                   ('j6lq01011_crc.fits', 'j6lq01011_crc_ref.fits'),
                   ('j6lq01ndq_flt.fits', 'j6lq01ndq_flt_ref.fits'),
                   ('j6lq01ndq_flc.fits', 'j6lq01ndq_flc_ref.fits')]

        # Prepare input files.
        download_file_cgi(self.tree, self.input_loc, asn0_file)
        shutil.copy(asn0_file, asn_file)
        tab = Table.read(asn_file, format='fits')
        for row in tab:
            if row['MEMTYPE'].startswith('PROD'):
                continue
            raw_file = row['MEMNAME'].lower().strip() + '_raw.fits'
            download_file_cgi(self.tree, self.input_loc, raw_file)

            # Prepare CRDS JREF files
            download_ref_ftp(raw_file)

        # Prepare local JREF files
        for fname in ('newwfc_osc.fits', 'pctetab_pcte.fits'):
            download_file_cgi(self.tree, self.input_loc, fname)

        # Run CALACS
        exe_to_use = HAS_CALXXX['acs']
        subprocess.call([exe_to_use, asn_file, '-v'])

        # Compare results
        all_okay = True
        creature_report = ''
        for actual, desired in outputs:

            # Download "truth" image
            download_file_cgi(self.tree, self.ref_loc, desired)

            fdiff = FITSDiff(actual, desired, ignore_keywords=ignore_keywords)
            creature_report += fdiff.report() + '\n\n'
            if not fdiff.identical and all_okay:
                all_okay = False

        if not all_okay:
            raise AssertionError(creature_report)
