# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis4.h"
# include "err.h"
# include "stisdq.h"
# include "stisdef.h"

static void GetSDC (StisInfo4 *);

/* This routine gets keyword values from an extension header.

   Phil Hodge, 1998 Oct 5:
	Change status value 1040 to GENERIC_ERROR_CODE.

   Phil Hodge, 1999 Sept 23:
	Get EXPTIME.

   Phil Hodge, 2000 Jan 5:
	Call GetLT0, because now we need both ltm and ltv for echelle data.

   Phil Hodge, 2000 Oct 18:
	Reset several bits in sdqflags, so that these flags will not be
	considered "serious" for any function in calstis4.  The bits that
	are reset are DATAMASKED (behind an occulting bar), HOTPIX (hot but
	probably corrected), and SMALLBLEM.

   Phil Hodge, 2001 Mar 7:
	For prism data, the input coordinate parameters will be for RA & Dec,
	so replace them with appropriate values.  Add function GetSDC, and
	call it for prism to get a2center and cdelt2.

   Nadia Dencheva, 2006 Apr 10:
	Get EXPSTART (needed for trace rotation).

   Phil Hodge, 2008 Dec 12:
	Get IMSET_OK.
*/

int GetGrpInfo4 (StisInfo4 *sts, Hdr *hdr) {

/* arguments:
StisInfo4 *sts  io: calibration switches and info
Hdr *hdr        i: header of current imset
*/

	int status;

	int sdqflags;			/* serious data quality flags */
	int use_def = 1;		/* use default if missing keyword */
	Bool value;			/* value of IMSET_OK keyword */

	/* Get IMSET_OK, which will be false if the current imset was
	   flagged as having zero exptime or constant pixel values.
	*/
	if ((status = Get_KeyB (hdr, "IMSET_OK", use_def, True, &value)) != 0)
	    return status;
	if (value)
	    sts->imset_ok = 1;
	else
	    sts->imset_ok = 0;

	/* Get the dispersion axis. */
	if ((status = Get_KeyI (hdr, "DISPAXIS", use_def, 1, &sts->dispaxis)))
	    return (status);
	if (sts->dispaxis < 1 || sts->dispaxis > 2) {
	    printf (
	"ERROR    DISPAXIS = %d is invalid for spectroscopic data.\n",
		sts->dispaxis);
	    return (GENERIC_ERROR_CODE);
	}

	/* Find out which data quality bits (all, by default) are considered
	   serious.  Note that we reset some bits, so that we can include
	   pixels that are flagged with these conditions.
	*/
	if ((status = Get_KeyI (hdr, "SDQFLAGS", use_def, 32767, &sdqflags)))
	    return (status);
	sdqflags &= ~DATAMASKED;	/* behind occulting bar */
	sdqflags &= ~HOTPIX;		/* hot (but probably just warm) pixel */
	sdqflags &= ~SMALLBLEM;		/* too small to be concerned about */
	sts->sdqflags = (short) sdqflags;

	/* Get the LTMi_i and LTVi keywords, converted to zero indexing. */
	if ((status = GetLT0 (hdr, sts->ltm, sts->ltv)))
	    return (status);
	sts->scale[0] = 1. / sts->ltm[0];
	sts->scale[1] = 1. / sts->ltm[1];

	/* Get the coordinate parameters. */
	if (sts->disp_type == PRISM_DISP) {

	    /* For prism data prior to 2-D rectification, the coordinate
		parameters in the header are for RA & Dec, but we need
		cross-dispersion spatial coordinates.
	    */
	    GetSDC (sts);
	    sts->crpix[1] = sts->crpix[1] * sts->ltm[1] + sts->ltv[1];
	    sts->cdelt[1] /= sts->ltm[1];

	} else {

	    if ((status = Get_KeyD (hdr, "CRPIX1", use_def, 0., &sts->crpix[0])))
		return (status);
	    if ((status = Get_KeyD (hdr, "CRPIX2", use_def, 0., &sts->crpix[1])))
		return (status);
	    if ((status = Get_KeyD (hdr, "CRVAL1", use_def, 0., &sts->crval[0])))
		return (status);
	    if ((status = Get_KeyD (hdr, "CRVAL2", use_def, 0., &sts->crval[1])))
		return (status);
	    if ((status = Get_KeyD (hdr, "CD1_1", use_def, 1., &sts->cdelt[0])))
		return (status);
	    if ((status = Get_KeyD (hdr, "CD2_2", use_def, 1., &sts->cdelt[1])))
		return (status);
	    /* Convert CRPIX to zero index. */
	    sts->crpix[0]--;
	    sts->crpix[1]--;
	}

	/* We only need the exposure time to verify that it's greater
	   than zero.  That's why the default is zero.
	*/
	if ((status = Get_KeyD (hdr, "EXPTIME", use_def, 0., &sts->exptime)))
	    return (status);

	/* Exposure start time is needed for trace rotation */
	if ((status = Get_KeyD (hdr, "EXPSTART", 0, 0., &sts->expstart)))
	  return (status);

	return (0);
}

/* This routine opens the SDC table, finds the row for PRISM, and reads
   A2CENTER into crpix[1] and CDELT2 into cdelt[1].  A2CENTER is taken
   instead of CRPIX2 because the PRISM data have not been 2-D rectified.
*/

static void GetSDC (StisInfo4 *sts) {

	IRAFPointer tp;
	IRAFPointer cp_opt_elem;
	IRAFPointer cp_a2center, cp_cdelt2;
	IRAFPointer cp_pedigree, cp_descrip;
	char opt_elem[STIS_CBUF+1];		/* grating or prism name */
	int nrows, row;
	int foundit;

	/* default values */
	sts->crpix[1] = PRISM_CRPIX2 - 1.;	/* zero indexed */
	sts->cdelt[1] = PRISM_CDELT2;		/* degrees per pixel */
	sts->crpix[0] = 0.;		/* the rest are not used */
	sts->cdelt[0] = 0.;
	sts->crval[0] = 0.;
	sts->crval[1] = 0.;

	tp = c_tbtopn (sts->sdctab.name, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("Warning  SDCTAB `%s' not found; default values used.\n",
		sts->sdctab.name);
	    clear_cvoserr();
	    return;
	}
	nrows = c_tbpsta (tp, TBL_NROWS);

	c_tbcfnd1 (tp, "OPT_ELEM", &cp_opt_elem);
	c_tbcfnd1 (tp, "A2CENTER", &cp_a2center);
	c_tbcfnd1 (tp, "CDELT2", &cp_cdelt2);
	if (cp_opt_elem == 0 || cp_a2center == 0 || cp_cdelt2 == 0) {
	    printf ("Warning  Column(s) not found in SDCTAB; defaults used.\n");
	    c_tbtclo (tp);
	    return;
	}
	c_tbcfnd1 (tp, "PEDIGREE", &cp_pedigree);
	c_tbcfnd1 (tp, "DESCRIP", &cp_descrip);

	foundit = 0;
	for (row = 1;  row <= nrows;  row++) {

	    c_tbegtt (tp, cp_opt_elem, row, opt_elem, STIS_CBUF);
	    if (c_iraferr()) {
		strcpy (opt_elem, "dummy");
		clear_cvoserr();
	    }

	    if (SameString (opt_elem, sts->opt_elem)) {

		foundit = 1;
		c_tbegtd (tp, cp_a2center, row, &sts->crpix[1]);
		c_tbegtd (tp, cp_cdelt2, row, &sts->cdelt[1]);
		sts->crpix[1] -= 1.;		/* zero indexed */
		sts->cdelt[1] /= 3600.;		/* degrees per pixel */

		RowPedigree (&sts->sdctab, row, tp, cp_pedigree, cp_descrip);

		break;
	    }
	}

	c_tbtclo (tp);

	if (!foundit)
	    printf ("Warning  PRISM not found in SDCTAB; defaults used.\n");
}
