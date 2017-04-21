# include <xtables.h>
# include "hstio.h"	/* defines HST I/O functions */

# include "wf3.h"
# include "wf3info.h"
# include "err.h"
# include "wf3dq.h"

typedef struct {
        IRAFPointer tp;                 /* pointer to table descriptor */
        IRAFPointer cp_xstart, cp_ystart; /* starting pixel for bad values */
        IRAFPointer cp_length;          /* this many bad values */
        IRAFPointer cp_axis;            /* repeat along X or Y axis (1 or 2) */
        IRAFPointer cp_flag;            /* this is the data quality value */
	IRAFPointer cp_amp;		/* selection columns */
	IRAFPointer cp_ccdchip;
	IRAFPointer cp_ccdgain;
	int axlen1, axlen2;		/* DQ array size specified in header */
        int nrows;                      /* number of rows in table */
} TblInfo;

typedef struct {
        /* values read from a table row */
        int xstart, ystart;
        int length;
        int axis;
        short flag;
	/* values used for selecting row to be used */
	char ccdamp[SZ_CBUF+1];
	float ccdgain;
	int ccdchip;
} TblRow;

extern int status;

static int dqicorr (WF3Info *, SingleNicmosGroup *, SingleNicmosGroup *);


/* DODQIIR: Load the data quality initialization table (BPIXTAB) and
**          call dqicorr for each readout of a MultiAccum.
**
**	The BPIXTAB should contain the following integer header
**	parameters:
**		SIZAXIS1, SIZAXIS2: full size of data quality array;
**			these are expected to be 1024x1024 for WFC3 IR data
**	and the following integer columns:
**		PIX1, PIX2: starting pixel (one indexed) of a range of
**			pixels to be assigned an initial value
**		LENGTH: number of pixels to be assigned value
**		AXIS: 1 --> X axis, 2 --> Y axis
**		VALUE: value to be ORed with data quality array
**	Each row of the table specifies a set of LENGTH pixels to be flagged
**	in either the X or Y direction with the value VALUE. The value
**	assigned for each pixel will be the OR of VALUE and any previous
**	value at that pixel.
**
** Revision history:
** H.Bushouse	24-May-2001	Revised from calnica maskcorr for WFC3.
** H.Bushouse	19-Mar-2002	Modifications to handle subarrays.
** H.Bushouse	21-Jun-2002	Revised to select BPIXTAB rows based on
**				amp, chip, and gain values.
** H.Bushouse	16-Oct-2003	Revised to use floating-point gain values.
** H.Bushouse	18-Nov_2008	Updated to check for missing CCDGAIN and CCDAMP
**				columns in BPIXTAB and default to a match with
**				the science data (same logic as in lib/dodqi.c).
**				(PR 61436)
*/

int doDQIIR (WF3Info *wf3, MultiNicmosGroup *input) {

/* Arguments:
**	wf3	 i: WF3 info structure
**	input	io: image to be masked
*/

	/* Local variables */
	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow  tabrow;		/* values read from a table row */

	/* mappings from one coordinate system to another */
	double ri_m[2], ri_v[2];	/* reference to image */
	int row;			/* loop index for row number */
	int nrows;			/* number of rows applied */
	int npix[2];			/* size of current image */
	int sameamp, samegain, samechip;

	SingleNicmosGroup mask;		/* temporary DQ mask image */

	/* Function definitions */
	int  GetLT0 (Hdr *, double *, double *);
	int  SameInt (int, int);
	int  SameFlt (float, float);
	int  SameString (char *, char *);

	int  OpenBpixTab (char *, TblInfo *);
	int  ReadBpixTab (TblInfo *, int, TblRow *);
	int  CloseBpixTab (TblInfo *);

	void DQINormal (DQHdrData *, double *, TblRow *);
	void PrSwitch (char *, int);

	/* Do the DQ initialization */
	if (wf3->dqicorr == PERFORM) {

	    /* Get the linear transformations */
	    if (GetLT0 (&(input->group[0].sci.hdr), ri_m, ri_v))
		return (status);

	    /* Open the DQ ref table, find columns, etc. */
	    if (OpenBpixTab (wf3->bpix.name, &tabinfo))
		return (status);

	    /* Size of current (science) image */
	    npix[0] = input->group[0].dq.data.nx;
	    npix[1] = input->group[0].dq.data.ny;

	    /* Create scratch image */
	    initSingleNicmosGroup (&mask);
	    if (allocSingleNicmosGroup (&mask, npix[0], npix[1]) == -1) {
		trlerror ("doDQI couldn't allocate data quality array.");
		return (status = OUT_OF_MEMORY);
	    }

	    /* Read each row of the table and load DQ values into mask image */
	    nrows = 0;
	    for (row = 1; row <= tabinfo.nrows; row++) {
		 if (ReadBpixTab (&tabinfo, row, &tabrow)) {
		     trlerror ("Error reading BPIXTAB.");
		     return (status);
		 }

		 /* If CCDAMP column does not exist or it has a wildcard value,
		 ** always return a match. */
		 if (tabinfo.cp_amp == 0 || SameString(tabrow.ccdamp, "N/A")) {
		     sameamp = 1;
		 } else {
		     sameamp = SameString (tabrow.ccdamp, wf3->ccdamp);
		 }

		 /* If CCDGAIN column does not exist or it has a wildcard value,
		 ** always return a match. */
		 if (tabinfo.cp_ccdgain == 0 || tabrow.ccdgain == -999) {
		     samegain = 1;
		 } else {
		     samegain = SameFlt (tabrow.ccdgain, wf3->ccdgain);
		 }

		 /* If CCDCHIP column does not exist or it has a wildcard value,
		 ** always return a match. */
		 if (tabinfo.cp_ccdchip == 0 || tabrow.ccdchip == -999) {
		     samechip = 1;
		 } else {
		     samechip = SameInt (tabrow.ccdchip, wf3->chip);
		 }

		 /* Check for a match with selection criteria */
		 if (sameamp && samegain && samechip) {

		     /* Assign the flag values to all relevant pixels */
		     DQINormal (&(mask.dq), ri_v, &tabrow);
		     nrows += 1;
		 }
	    }

	    if (CloseBpixTab (&tabinfo))	/* done with the table */
		return (status);

	    if (nrows == 0) {
	        sprintf (MsgText, "No rows from BPIXTAB applied to DQ array.");
	        trlwarn (MsgText);
	    }

	    /* Loop over all MultiAccum groups, applying mask to each one */
	    for (wf3->group=wf3->ngroups; wf3->group >= 1; wf3->group--) {
		 if (dqicorr (wf3, &(input->group[wf3->group-1]), &mask))
		     return (status);
	    }

	    /* Free the mask image */
	    freeSingleNicmosGroup (&mask);

	    PrSwitch ("dqicorr", COMPLETE);
	}

	/* Successful return */
	return (status = 0);
}

/* DQICORR: Combine static bad pixel mask with science image DQ array.
** The input DQ array is logically "OR'd" with the mask DQ array.
** Also set all DQ values to DETECTORPROB if there was a TDF transition.
** The input SCI, ERR, SAMP, and TIME arrays are unchanged.
**
** Revision history:
** H.Bushouse	Oct. 2000	Copied from calnica maskcorr for WFC3.
** H.Bushouse	08-May-2002	Changed DQ flag macro "BADPIX" to
**				"DETECTORPROB". Added use of wf3dq.h.
** M. Sosey     11-Feb-2013, the DQ array no longer reflects TDFTRANS
*/

static int dqicorr (WF3Info *wf3, SingleNicmosGroup *input,
		    SingleNicmosGroup *mask) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: image to be masked
**	mask	 i: mask image
*/

	/* Function definitions */
	void aor (SingleNicmosGroup *, SingleNicmosGroup *);

	/* Combine the DQ mask with the input DQ */
	aor (input, mask);

	/* Successful return */
	return (status = 0);
}

