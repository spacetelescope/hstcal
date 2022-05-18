# include <stdio.h>
# include <string.h>

#include "hstcal.h"
# include "hstio.h"
# include "xtables.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

# define NUMCOLS 26		/* Defined locally to be specific to CCDTAB */

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_amp;		/* column descriptors */
	IRAFPointer cp_ccdchip;	
	IRAFPointer cp_ccdgain;
	IRAFPointer cp_bin1;
	IRAFPointer cp_bin2;
	IRAFPointer cp_ccdoffset[NAMPS];
	IRAFPointer cp_bias[NAMPS];
	IRAFPointer cp_atodgain[NAMPS];
	IRAFPointer cp_readnoise[NAMPS];
	IRAFPointer cp_ampx;
	IRAFPointer cp_ampy;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char ccdamp[SZ_CBUF+1];
	float ccdgain;
	int ccdchip;
	int bin[2];
	int ccdoffset[NAMPS];
	float bias[NAMPS];
	float atodgain[NAMPS];
	float readnoise[NAMPS];
	int ampx;
	int ampy;
} TblRow;

static int OpenCCDTab (char *, TblInfo *);
static int ReadCCDTab (TblInfo *, int, TblRow *);
static int CloseCCDTab (TblInfo *);

/* This routine gets the gain, bias and readnoise for the CCD from
   the CCD parameters table (keyword name CCDTAB).

   The CCD parameters table should contain the following:
	header parameters:
		none needed
	columns:
	CCDAMP:  identifies which amp was used (string A-D)
        CCDCHIP: identifies which chip corresponds to the extension      
	CCDGAIN:  commanded gain of the CCD (float)
	BINAXIS1, BINAXIS2:  commanded bin sizes (int)
	CCDOFST[A,B,C,D]:  commanded bias for amp 1,2,3,4 of CCD (int)
	CCDBIAS[A,B,C,D]:  calibrated bias offset for each amp (double)
	ATODGN[A,B,C,D]:  actual gain of amps 1-4 used for readouts(double)
	READNSE[A,B,C,D]:  typical value of readout noise for each amp(double)
	SATURATE:  CCD saturation threshold
	AMPX:	first column affected by second amp on multiamp readout (int)
	AMPY: 	first line affected by second amp(s) on multiamp readout (int)

   The table is read to find the row for which the CCDAMP matches the
   expected amplifier string (from the image header) and for which the
   commanded CCD gain, commanded CCD bias, CCDCHIP, and bin sizes match the
   values read from the image header.  The CCDAMP value from the image header
   must match exactly the value read from the row, including the order of
   the AMPs listed in CCDAMP.  The AMPs should be listed in increasing
   value, from A to D, for whichever AMPs are used.  The matching row is
   then read to get the actual gain, bias, readnoise, and saturation
   level.
   
** This was only modified to read ACS specific columns from CCDTAB.  May
** still need to be modified to support multiple AMP readout from a single
** chip by adding AMPX (and AMPY) column(s), which specify the column (and row)
** boundaries for each AMP configuration for that chip. 
** AMPX will be set to 0 (zero) for 1 AMP readout of the chip.
** In addition, AMPY gets set to 0 (zero) for 1 AMP readout.
    
    19 Feb 1999 WJH:
        Revised to handle more columns correctly, increased NUMCOLS to 23
        and added CCDOFST[C,D] to the list of columns in used for printing
        missing columns.

    16 Nov 2001 H.A.Bushouse:
	Revised to replace CCDBIAS column with CCDBIAS[A,B,C,D] in table.
	Also use CCDOFST[A,B,C,D] to select appropriate row as well.

    24 Jul 2002: H.Bushouse:
	Upgraded to use BINAXIS1,BINAXIS2 as additional reference table
	row selectors, along with amp, gain, chip, and offset. This had
	been removed from original calstis code for use in calacs, which
	doesn't use binning.

    16 Oct 2002: H.Bushouse:
	Upgraded to use floating-point commanded gain values for WFC3.

    20 Feb 2004: H.Bushouse:
	Changed use of wf3->binaxis to wf3->bin to get proper handling of
	binned science images.

    27 Aug 2008: H.Bushouse:
	Modified to reset ampx to <= dimx only for WFC3/UVIS images; leave
	it alone for WFC3/IR images.

    21 Oct 2009: H.Bushouse:
	Added computation of mean_gain to GetCCDTab. mean_gain is now used
	in flatcorr steps for doing gain conversion.

    16 Feb 2022: M. De La Pena:
    The "saturate" variable became obsolete once the full-well saturation
    map was implemented. Removed "saturate" as a cleanup operation.
*/

int GetCCDTab (WF3Info *wf3, int dimx, int dimy) {

/* arguments:
WF3Info *wf3     io: calibration switches, etc
int     dimx      i: number of columns in exposure
int     dimy      i: number of lines in exposure
*/

	extern int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int i;
	char amp[5] = "ABCD";
	
	int foundit;		/* has the correct row been found? */
	int RowPedigree (RefTab *, int, IRAFPointer, IRAFPointer, IRAFPointer);
	int SameInt (int, int);
	int SameFlt (float, float);
	int SameString (char *, char *);

	/* Open the CCD parameters table and find columns. */
	if (OpenCCDTab (wf3->ccdpar.name, &tabinfo))
	    return (status);

	/* Check each row for a match with ccdamp, ccdgain, ccdoffst,
	   binaxis1, and binaxis2, and get info from the matching row.
	*/

	foundit = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    /* Read the current row into tabrow. */
	    if (ReadCCDTab (&tabinfo, row, &tabrow))
		return (status);

	    if (SameString (tabrow.ccdamp,    wf3->ccdamp) &&
		SameFlt (tabrow.ccdgain,      wf3->ccdgain) &&
		SameInt (tabrow.ccdchip,      wf3->chip) &&
		SameInt (tabrow.ccdoffset[0], wf3->ccdoffset[0]) &&
		SameInt (tabrow.ccdoffset[1], wf3->ccdoffset[1]) &&
		SameInt (tabrow.ccdoffset[2], wf3->ccdoffset[2]) &&
		SameInt (tabrow.ccdoffset[3], wf3->ccdoffset[3]) &&
		SameInt (tabrow.bin[0],       wf3->bin[0]) &&
		SameInt (tabrow.bin[1],       wf3->bin[1])) {

		foundit = 1;
		if (RowPedigree (&wf3->ccdpar, row,
			tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip))
		    return (status);
		if (wf3->ccdpar.goodPedigree == DUMMY_PEDIGREE) {
		    sprintf (MsgText, "Row %d of CCDTAB is DUMMY.", row);
			trlwarn (MsgText);
		}
		wf3->mean_gain = 0;
		for (i = 0; i < NAMPS; i++){
			wf3->mean_gain += tabrow.atodgain[i];
			/* If the amp is used, keep the value,
			** otherwise set to zero*/
			if (strchr(wf3->ccdamp, amp[i]) != NULL) {
				wf3->atodgain[i] = tabrow.atodgain[i];
				wf3->readnoise[i] = tabrow.readnoise[i];
				wf3->ccdbias[i] = tabrow.bias[i];
				
			} else {
				wf3->atodgain[i] = 0.;
				wf3->readnoise[i] = 0.;			
				wf3->ccdbias[i] = 0.;
			}
		}
		wf3->mean_gain /= NAMPS;

        /* For WFC3/UVIS exposures, correct ampx to match the actual
           size of the image for more seamless processing of subarrays.
           Leave ampx as is for WFC3/IR exposures. */
		if (wf3->detector == CCD_DETECTOR)
		    wf3->ampx = (tabrow.ampx > dimx) ? dimx : tabrow.ampx;
		else
		    wf3->ampx = tabrow.ampx;
		wf3->ampy = tabrow.ampy;		

		break;
	    }
	}

	if (!foundit) {
	    sprintf (MsgText, "Matching row not found in CCDTAB `%s'.", wf3->ccdpar.name);
	    trlerror (MsgText);
	    sprintf (MsgText,
		    "CCDCHIP %d, CCDAMP %s, CCDGAIN %g, CCDOFFST %d,%d,%d,%d.",
		     wf3->chip, wf3->ccdamp, wf3->ccdgain, wf3->ccdoffset[0],
		     wf3->ccdoffset[1], wf3->ccdoffset[2], wf3->ccdoffset[3]);
	    trlerror (MsgText);
	    sprintf (MsgText, "BINAXIS %d,%d.",
		     wf3->bin[0], wf3->bin[1]);
	    trlerror (MsgText);
		
	    CloseCCDTab (&tabinfo);
	    return (status = TABLE_ERROR);
	}

	if (CloseCCDTab (&tabinfo))		/* close the table */
	    return (status);

	return (status);
}

/* This routine opens the CCD parameters table, finds the columns that we
   need, and gets the total number of rows in the table.  The columns are 
   CCDAMP, CCDGAIN, CCDBIAS, and READNSE.
*/

static int OpenCCDTab (char *tname, TblInfo *tabinfo) {

	extern int status;

	int nocol[NUMCOLS];
	int i, j, missing;
	
	char *colnames[NUMCOLS] ={"CCDAMP", "CCDCHIP", "CCDGAIN", "BINAXIS1", "BINAXIS2", "CCDOFSTA", "CCDOFSTB", "CCDOFSTC", "CCDOFSTD","CCDBIASA", "CCDBIASB", "CCDBIASC", "CCDBIASD", "ATODGNA", "ATODGNB", "ATODGNC", "ATODGND", "READNSEA", "READNSEB", "READNSEC", "READNSED", "AMPX", "AMPY", "SATURATE"};

	int PrintMissingCols (int, int, int *, char **, char *, IRAFPointer);

	for (j = 0; j < NUMCOLS; j++)
		nocol[j] = NO;

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    sprintf (MsgText, "CCDTAB `%s' not found.", tname);
	    trlerror (MsgText);
		return (status = OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "CCDAMP", &tabinfo->cp_amp);
	c_tbcfnd1 (tabinfo->tp, "CCDCHIP", &tabinfo->cp_ccdchip);
	c_tbcfnd1 (tabinfo->tp, "CCDGAIN", &tabinfo->cp_ccdgain);
	c_tbcfnd1 (tabinfo->tp, "BINAXIS1", &tabinfo->cp_bin1);
	c_tbcfnd1 (tabinfo->tp, "BINAXIS2", &tabinfo->cp_bin2);
	c_tbcfnd1 (tabinfo->tp, "CCDOFSTA", &tabinfo->cp_ccdoffset[0]);
	c_tbcfnd1 (tabinfo->tp, "CCDOFSTB", &tabinfo->cp_ccdoffset[1]);
	c_tbcfnd1 (tabinfo->tp, "CCDOFSTC", &tabinfo->cp_ccdoffset[2]);
	c_tbcfnd1 (tabinfo->tp, "CCDOFSTD", &tabinfo->cp_ccdoffset[3]);
	c_tbcfnd1 (tabinfo->tp, "CCDBIASA", &tabinfo->cp_bias[0]);
	c_tbcfnd1 (tabinfo->tp, "CCDBIASB", &tabinfo->cp_bias[1]);
	c_tbcfnd1 (tabinfo->tp, "CCDBIASC", &tabinfo->cp_bias[2]);
	c_tbcfnd1 (tabinfo->tp, "CCDBIASD", &tabinfo->cp_bias[3]);
	c_tbcfnd1 (tabinfo->tp, "ATODGNA", &tabinfo->cp_atodgain[0]);
	c_tbcfnd1 (tabinfo->tp, "ATODGNB", &tabinfo->cp_atodgain[1]);
	c_tbcfnd1 (tabinfo->tp, "ATODGNC", &tabinfo->cp_atodgain[2]);
	c_tbcfnd1 (tabinfo->tp, "ATODGND", &tabinfo->cp_atodgain[3]);
	c_tbcfnd1 (tabinfo->tp, "READNSEA", &tabinfo->cp_readnoise[0]);
	c_tbcfnd1 (tabinfo->tp, "READNSEB", &tabinfo->cp_readnoise[1]);
	c_tbcfnd1 (tabinfo->tp, "READNSEC", &tabinfo->cp_readnoise[2]);
	c_tbcfnd1 (tabinfo->tp, "READNSED", &tabinfo->cp_readnoise[3]);
	c_tbcfnd1 (tabinfo->tp, "AMPX", &tabinfo->cp_ampx);
	c_tbcfnd1 (tabinfo->tp, "AMPY", &tabinfo->cp_ampy);
	
	/* Initialize counters here... */
	missing = 0;
	i=0;
		
    /* Increment i for every column, mark only missing columns in
        nocol as YES.  WJH 27 July 1999
    */
	if (tabinfo->cp_amp == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_ccdchip == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_ccdgain == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_bin1 == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_bin2 == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_ccdoffset[0] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_ccdoffset[1] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_ccdoffset[2] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_ccdoffset[3] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_bias[0] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_bias[1] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_bias[2] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_bias[3] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_atodgain[0] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_atodgain[1] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_atodgain[2] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_atodgain[3] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_readnoise[0] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_readnoise[1] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_readnoise[2] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_readnoise[3] == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_ampx == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_ampy == 0 ) { missing++; nocol[i] = YES;} i++;
	
	if (PrintMissingCols (missing, NUMCOLS, nocol, colnames, "CCDTAB", tabinfo->tp) )
		return(status);
		
	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (status);
}

int PrintMissingCols (int missing, int numcols, int *nocol, char *colnames[], char *tabname, IRAFPointer tp) {

/* Parameters:
int missing			i: number of missing columns
int numcols			i: number of columns expected in table
int *nocol			i: array of YES/NO for each column, YES means missing
char **colnames		i: array of colnames
char *tabname		i: name of table columns were read from
IRAFPointer tp		i: pointer to table, close it if necessary
*/

	extern int status;
	int j;
	/* If any columns are missing... */
	if (missing) {
 	    sprintf (MsgText,"%d columns not found in %s.", missing, tabname);
		trlerror (MsgText);
       
		for (j=0; j< numcols; j++) {
			/* Recall which ones were marked missing... */
			if (nocol[j]) {
				/*... and print out that column's name */
	    		sprintf (MsgText,"Column %s not found in %s.", colnames[j], tabname);
				trlerror (MsgText);
			}
		}
	    c_tbtclo (tp);
	    return (status = COLUMN_NOT_FOUND);
	}
	return(status);
}
/* This routine reads the relevant data from one row.  The amplifier
   number, CCD gain, bias, and readnoise are gotten.
*/

static int ReadCCDTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	extern int status;

	c_tbegtt (tabinfo->tp, tabinfo->cp_amp, row,
			tabrow->ccdamp, SZ_CBUF);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_ccdchip, row, &tabrow->ccdchip);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_ccdgain, row, &tabrow->ccdgain);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_bin1, row, &tabrow->bin[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_bin2, row, &tabrow->bin[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_ccdoffset[0], row, &tabrow->ccdoffset[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);		
	c_tbegti (tabinfo->tp, tabinfo->cp_ccdoffset[1], row, &tabrow->ccdoffset[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_ccdoffset[2], row, &tabrow->ccdoffset[2]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);		
	c_tbegti (tabinfo->tp, tabinfo->cp_ccdoffset[3], row, &tabrow->ccdoffset[3]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_bias[0], row, &tabrow->bias[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_bias[1], row, &tabrow->bias[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_bias[2], row, &tabrow->bias[2]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_bias[3], row, &tabrow->bias[3]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
		
	c_tbegtr (tabinfo->tp, tabinfo->cp_atodgain[0], row, &tabrow->atodgain[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_atodgain[1], row, &tabrow->atodgain[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_atodgain[2], row, &tabrow->atodgain[2]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_atodgain[3], row, &tabrow->atodgain[3]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_readnoise[0], row, &tabrow->readnoise[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_readnoise[1], row, &tabrow->readnoise[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_readnoise[2], row, &tabrow->readnoise[2]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegtr (tabinfo->tp, tabinfo->cp_readnoise[3], row, &tabrow->readnoise[3]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_ampx, row, &tabrow->ampx);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_ampy, row, &tabrow->ampy);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
		
	return (status);
}

/* This routine closes the ccdtab table. */

static int CloseCCDTab (TblInfo *tabinfo) {

	extern int status;

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	return (status);
}
