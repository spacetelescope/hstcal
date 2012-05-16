# include <stdlib.h>
# include "xtables.h"

# include "stis.h"
# include "calstis6.h"

/* Initialize the calstis6 structure.  This includes information about the
   input and output images, calibration files, and flags to specify which
   calibration steps are to be done.




   Revision history:
   ----------------
   20 Feb 97  -  Adapted from similar routine in calstis7 (I.Busko)
   09 Apr 97  -  Changes after code review (IB):
                 - replaced literal by ALL_BITS constant.
   25 Apr 97  -  Removed esptab (IB)
   08 May 97  -  Added obsmode and info flags (IB)
   14 May 97  -  Flag signals calstis6 is being run from the pipeline (IB)
   23 May 97  -  Add DefSwitch function call (IB)
   03 Sep 97  -  Add ATODGAIN correction (IB)
   26 Jan 98  -  Add PCTAB support (IB)
   05 Feb 98  -  Remove pipeline flag / do not initialize calib. switches (IB)
   10 Apr 98  -  Replace debug switch by extrloc, remove debug file (IB)
   13 Apr 98  -  Add extraction info columns to output table (IB)
   22 Jun 98  -  Global croscor fit mode (IB)
   23 Sep 98  -  Output weights image name (IB)
   18 Nov 98  -  Background value and error from command line (IB)
   06 Dec 99  -  Profile DQ array, sigma-clip thresholds (IB)
   05 Jan 00  -  Lee filter window size (IB)
   17 Feb 00  -  Default IMSET selection (IB)
   04 Dec 00  -  Default subsampling factor (IB)
   08 Jan 02  -  Time-dependent sensitivity (IB)
   16 Apr 02  -  Blaze shift (IB)
   24 Jul 02  -  Background smoothing (IB)
   20 Jul 04  -  Call InitRefTab for ccdtab (PEH)
   27 Dec 04  -  Remove InitRefTab from this file; set detector_temp (PEH)
   08 Apr 05  -  Call InitRefTab for gactab (PEH)
   20 Apr 05  -  Initialize xoffset to 0. (PEH)
   11 Jan 06  -  Remove mofftab (PEH)
*/

void StisInit6 (StisInfo6 *sts) {

	/* Assign default values. */

	sts->input[0]       = '\0';
	sts->output[0]      = '\0';
	sts->rootname[0]    = '\0';
	sts->outw[0]        = '\0';
	sts->cl_a2center    = NO_POSITION;	/* command-line parameters */
	sts->xoffset        = 0.;
	sts->maxsearch      = NO_RANGE;
	sts->extrsize       = NO_SIZE;
	sts->bksize[0]      = NO_SIZE;
	sts->bksize[1]      = NO_SIZE;
	sts->bkoffset[0]    = NO_SIZE;
	sts->bkoffset[1]    = NO_SIZE;
	sts->bktilt         = NO_TILT;
	sts->bkord          = NO_ORDER;
	sts->sporder        = NO_ORDER;
	sts->obsmode[0]     = '\0';
	sts->xtracalg[0]    = '\0';
	sts->aperture[0]    = '\0';
	sts->cc_global      = 0;
	sts->wavecal        = 0;
	sts->detector       = UNKNOWN_DETECTOR;
	sts->opt_elem[0]    = '\0';
	sts->atodgain       = 1.0;
	sts->nimages        = 1;
	sts->sdqflags       = ALL_BITS;		/* all bits set */
	sts->sdqflags_orig  = ALL_BITS;		/* all bits set */
	sts->exptime        = 1.;
	sts->detector_temp  = -1.;		/* invalid temperature */
	sts->hfactor        = 1.;
	sts->ltm[0]         = 1.;
	sts->ltm[1]         = 1.;
	sts->ltv[0]         = 0.;
	sts->ltv[1]         = 0.;
	sts->crscroff       = 0.;
	sts->nm_a2center    = 0.;
	sts->bck[0]         = 0.;
	sts->bck[1]         = 0.;
	sts->vbck[0]        = 0.;
	sts->vbck[1]        = 0.;
	sts->ebck           = 0.;
	sts->dqbck          = 0;
	sts->cenwave        = 0;
	sts->dispaxis       = 1;
	sts->verbose        = 0;
	sts->printtime      = 0;
	sts->extrloc        = 1;
	sts->dispinfo       = 0;
	sts->fluxinfo       = 0;
	sts->sgeoinfo       = 0;
        sts->ccdinfo        = 0;
	sts->cc_off         = NULL;		/* croscor offsets */
	sts->cc_spord       = NULL;
	sts->cc_size        = 0;
	sts->backval        = NO_VALUE;		/* background */
	sts->backerr        = 0.0;
	sts->optimal        = 0;
	sts->lfilter        = 17;
	sts->imset          = 0;
	sts->blazeshift     = NO_VALUE;
	sts->bks_mode       = BKS_OFF;
	sts->bks_order      = 3;
	sts->avoid1a        = 0;
	sts->avoid2a        = 0;
	sts->avoid1b        = 0;
	sts->avoid2b        = 0;

	sts->dither_offset[0] = 0.;
	sts->dither_offset[1] = 0.;
	sts->msm_offset[0]    = 0.;
	sts->msm_offset[1]    = 0.;
	sts->ap_offset[0]     = 0.;
	sts->ap_offset[1]     = 0.;

	sts->profile       = NULL;		/* profile generator */
	sts->profile_dq    = NULL;
	sts->profile_rejf  = NULL;
	sts->profile_minw  = NULL;
	sts->profile_maxw  = NULL;
	sts->profile_minp  = NULL;
	sts->profile_maxp  = NULL;
	sts->profile_sn    = NULL;
	sts->profile_msize = 0;
	sts->sclip         = 6.0;
	sts->psclip        = 5.0;
	sts->subscale      = 0.0;

	/* If calstis6 was called, it is presumed that extraction
           has to be performed anyway.
        */
	sts->x1d       = PERFORM;
	sts->x1d_o     = PERFORM;

	/* No info yet about whether reference files exist, and
	   therefore no values for pedigree and descrip.
	*/

	sts->sdstfile.name[0]      = '\0';		/* image */
	sts->sdstfile.pedigree[0]  = '\0';
	sts->sdstfile.descrip[0]   = '\0';
	sts->sdstfile.exists       = EXISTS_UNKNOWN;
	sts->sdstfile.goodPedigree = PEDIGREE_UNKNOWN;

	/* tables */
	InitRefTab (&sts->apdestab);
	InitRefTab (&sts->apertab);
	InitRefTab (&sts->phottab);
	InitRefTab (&sts->ccdtab);
	InitRefTab (&sts->tdstab);
	InitRefTab (&sts->disptab);
	InitRefTab (&sts->inangtab);
	InitRefTab (&sts->sptrctab);
	InitRefTab (&sts->xtrctab);
	InitRefTab (&sts->pctab);
	InitRefTab (&sts->gactab);
	InitRefTab (&sts->pftab);
	InitRefTab (&sts->pxtab);

	InitRefTab (&sts->echsctab);
	InitRefTab (&sts->exstab);
	InitRefTab (&sts->cdstab);
	InitRefTab (&sts->riptab);
	InitRefTab (&sts->srwtab);
	InitRefTab (&sts->halotab);
	InitRefTab (&sts->psftab);

	/* Initialize trace rotation angle. */
	sts->trace_rotation = 0.;
}



/* Initialize the elements of a TblDesc structure. */

void InitTblDesc (TblDesc *table) {

	table->tp         = 0;   /* to avoid complaining by the Alpha C */
	table->nrows      = 0;
	table->array_size = 0;
	table->sporder    = 0;
	table->npts       = 0;
	table->wave       = 0;
	table->gross      = 0;
	table->back       = 0;
	table->net        = 0;
	table->flux       = 0;
	table->error      = 0;
	table->dq         = 0;
        table->a2center   = 0;
	table->extrsize   = 0;
	table->bk1size    = 0;
	table->bk2size    = 0;
	table->bk1offset  = 0;	
	table->bk2offset  = 0;	
	table->maxsearch  = 0;
	table->extrlocy   = 0;
	table->cc_offset  = 0;
}



/* Initialize the elements of a ProfTblDesc structure. */

void InitProfTblDesc (ProfTblDesc *table) {

	table->tp             = 0;   /* to avoid complaining by the Alpha C */
	table->nrows          = 0;
	table->array_size     = 0;
	table->array_size_off = 0;
	table->sporder        = 0;
	table->npts           = 0;
	table->minwave        = 0;
	table->maxwave        = 0;
	table->s_n            = 0;
	table->profoff        = 0;
	table->prof           = 0;
}




/* Initialize the elements of a RowContents structure. */

void InitRowContents (RowContents *row) {

	row->sporder   = 0;
	row->npts      = 0;
        row->a2center  = 0.0F;
	row->extrsize  = 0.0F;
	row->bk1size   = 0.0F;
	row->bk2size   = 0.0F;
	row->bk1offset = 0.0F;	
	row->bk2offset = 0.0F;	
	row->maxsearch = 0;
	row->wave      = NULL;
	row->gross     = NULL;
	row->back      = NULL;
	row->net       = NULL;
	row->flux      = NULL;
	row->error     = NULL;
	row->dq        = NULL;
	row->extrlocy  = NULL;
	row->cc_offset = 0.0F;	
}
