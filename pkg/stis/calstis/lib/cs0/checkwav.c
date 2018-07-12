# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis0.h"
# include "hstcalerr.h"
# include "stisdef.h"

static int goodImage (SingleGroup *);

/* This routine checks whether the wavecal file has exptime <= 0 or
   all zero pixel values in every image set.  If so, WAVECORR, X1DCORR
   and X2DCORR will all be set to OMIT.

   Phil Hodge, 2007 May 2:
	Function created.
*/

int checkWav (StisInfo *sts) {

/* arguments:
StisInfo *sts         i: calibration flags and other info
*/

	int status;

	IODescPtr im;		/* descriptor for an image */
	SingleGroup x;		/* one image set */
	Hdr phdr;		/* primary header */
	Hdr *hdr;		/* pointer to extension header */
	double texptime;	/* sum of exposure times in all imsets */
	double exptime;		/* exposure time */
	int nextend;		/* number of FITS extensions in wavecal file */
	int nimsets;		/* number of image sets */
	int extver;		/* loop index over image sets */
	int all_bad;		/* true if no imset has good data (boolean) */
	int use_default = 1;	/* use default if keyword is missing */

	/* Read primary header of wavfile into phdr. */
	initHdr (&phdr);
	im = openInputImage (sts->wavfile, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (im, &phdr);		/* get primary header */
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);

	/* First check texptime. */
	if ((status = Get_KeyD (&phdr, "TEXPTIME", use_default, 0., &texptime)))
	    return (status);
	if (texptime <= 0.) {
	    printf (
"Warning  Wavecal file has zero total exposure time,\\\n");
	    printf (
"         so WAVECORR, X1DCORR, X2DCORR will be set to OMIT.\n");
	    sts->sci_wavecorr = OMIT;
	    sts->sci_2d_rect = OMIT;
	    sts->sci_1d_extract = OMIT;
	    freeHdr (&phdr);
	    return (0);
	}

	/* Get the number of extensions. */
	if ((status = Get_KeyI (&phdr, "NEXTEND",
				use_default, EXT_PER_GROUP, &nextend)))
	    return (status);
	nimsets = nextend / EXT_PER_GROUP;

	freeHdr (&phdr);

	/* Loop over image sets.  If the exposure time is greater than zero
	   and there is at least one non-zero pixel value, return without
	   doing anything further.  If no imset has any non-zero data,
	   turn off wavecorr, x1dcorr, x2dcorr.
	*/

	all_bad = 1;			/* pessimistic initial value */
	for (extver = 1;  extver <= nimsets;  extver++) {
	    /* Get the current image set. */
	    initSingleGroup (&x);
	    getSingleGroup (sts->wavfile, extver, &x);
	    if (hstio_err())
		return (OPEN_FAILED);
	    hdr = &x.sci.hdr;
	    if ((status = Get_KeyD (hdr, "EXPTIME", use_default, 0., &exptime)))
		return (status);
	    if (exptime > 0. && goodImage (&x)) {
		all_bad = 0;
		freeSingleGroup (&x);
		break;
	    }
	    freeSingleGroup (&x);
	}
	if (all_bad) {
	    printf (
"Warning  Wavecal file has zero exposure time or zero data values, \\\n");
	    printf (
"         so WAVECORR, X1DCORR, X2DCORR will be set to OMIT.\n");
	    sts->sci_wavecorr = OMIT;
	    sts->sci_2d_rect = OMIT;
	    sts->sci_1d_extract = OMIT;
	}

	return (0);
}

static int goodImage (SingleGroup *x) {

	int nx, ny;
	int i, j;
	int good;

	nx = x->sci.data.nx;
	ny = x->sci.data.ny;

	good = 0;
	for (i = 0;  i < nx;  i++) {
	    for (j = 0;  j < ny;  j++) {
		if (Pix (x->sci.data, i, j) > 0.) {
		    good = 1;
		    break;
		}
	    }
	}

	return (good);
}
