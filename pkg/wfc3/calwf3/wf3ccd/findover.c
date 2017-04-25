# include <stdio.h>
# include "hstio.h"

# include "xtables.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

# define FULL_FRAME_READOUT  1
# define SUBARRAY_READOUT    2
# define BINNED_READOUT      3

# define NUMCOLS 28

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_amp;		/* column descriptors */
	IRAFPointer cp_chip;
	IRAFPointer cp_binx;
	IRAFPointer cp_biny;	
	IRAFPointer cp_nx;
	IRAFPointer cp_ny;
	IRAFPointer cp_trimx1;
	IRAFPointer cp_trimx2;
	IRAFPointer cp_trimx3;
	IRAFPointer cp_trimx4;
	IRAFPointer cp_trimy1;
	IRAFPointer cp_trimy2;
	IRAFPointer cp_vx1;
	IRAFPointer cp_vx2;
	IRAFPointer cp_vx3;
	IRAFPointer cp_vx4;
	IRAFPointer cp_vy1;
	IRAFPointer cp_vy2;
	IRAFPointer cp_vy3;
	IRAFPointer cp_vy4;
	IRAFPointer cp_biassecta1;
	IRAFPointer cp_biassecta2;		
	IRAFPointer cp_biassectb1;
	IRAFPointer cp_biassectb2;		
	IRAFPointer cp_biassectc1;
	IRAFPointer cp_biassectc2;		
	IRAFPointer cp_biassectd1;
	IRAFPointer cp_biassectd2;		
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char ccdamp[SZ_CBUF+1];
	int chip;
	int nx;
	int ny;
	int trimx[4];
	int trimy[2];
	int vx[4];
	int vy[4];
	int biassecta[2];
	int biassectb[2];	
	int biassectc[2];	
	int biassectd[2];	
	int bin[2];
} TblRow;


static int OpenOverTab (char *, TblInfo *);
static int ReadOverTab (TblInfo *, int, TblRow *);
static int CloseOverTab (TblInfo *);
	
/* This routine assigns the widths of the overscan region on each edge
   of the image by assuming particular region sizes, depending on the
   binning and which amplifier was used for readout.
  
**   
**   This routine reads in the values from the OSCNTAB reference table.
**
**   The OSCNTAB refers to all regions in image coordinates, even when 
**	   the image is binned.


	    parallel readout direction, amp A
	    |
	    |
	    v
	(VX1,VX2)     (VX2,VY2)
  A       /         \     B
    -----/-----------+--- 
   |   |/     |      |   | } TRIMY2
   |   +------+------    |
   |   |      |      |   | <--- serial readout direction, amp A
   |   |      |      |   | 
   | - | -----+------|---|<-- AMPY
   |   |      |      |   |
   |   |      |      |   |
   |   |      |      |   |
   |    ------+------    |
   |   |      |      |   | } TRIMY1
    ---------------------
  C /  \      ^       / \ D
   A1  A2     |      B1  B2 
   			AMPX

	A,B,C,D   - Amps
	AMPX 	  - First column affected by second AMP
	AMPY 	  - First line affected by second set of AMPS 	 
	(VX1,VY1) - image coordinates of virtual overscan region origin
	(VX2,VY2) - image coordinates of top corner of virtual overscan region
	A1,A2 	  - beginning and ending columns for leading bias section
					(BIASSECTA1,BIASSECTA2 from OSCNTAB)
	B1,B2 	  - beginning and ending columns for trailing bias section
					(BIASSECTB1,BIASSECTB2 from OSCNTAB)
	TRIMX1    - Number of columns to trim off beginning of each line,
					contains A1,A2
	TRIMX2    - Number of columns to trim off end of each line,
					contains B1,B2
	TRIMY1    - Number of lines to trim off beginning of each column
	TRIMY2    - Number of line to trim off end of each column
		

   Warren Hack, 1998 June 2:
   	Initial version based on GetCCDTab...
	
   Warren Hack, 1998 Oct 20:
	Second version.  Condensed all input/output through ACSInfo 
	structure, and included both trailing and leading overscan
	region designations.  
   Warren Hack, 1999 Feb 19:
        Revised to read CCDCHIP column to match CCDCHIP keyword.
   Howard Bushouse, 2000 Aug 29:
	Revised for WFC3 use.
   H.Bushouse, 2001 Nov 16:
	Updated to track CALACS changes - Fixed a bug in calculating the second
	trim region in a subarray; Removed another bug in determining trim
	regions for subarrays; modified checks to be relative to trimmed size
	of a full frame; check to see if overlap beginning trim region first
	(and do both), then check for second region overlap.
   H.Bushouse, 2002 Jan 28:
	Added trimx3, trimx4, biassectc, biassectd, vx3, vx4, vy3, vy4 to
	accomodate extra serial virtual overscan regions in WFC3 UVIS CCD's.
   H.Bushouse, 2002 Jun 17:
	Modified to subtract 1 from biassect values read in from OSCNTAB in
	order to conform with C style zero-indexing (following CALACS changes).
   H.Bushouse, 2003 Oct 23:
	Modified for WFC3 to zero-out both serial and parallel virtual
	biassect and trim values for subarray images. Also fixed a bug in
	which one of the biassect values was not being converted from
	1-indexed to 0-indexed for subarray images.
   H.Bushouse, 2005 Feb 14:
	Modified FindOverscan routine for WFC3 IR channel to select oscntab
        row based on image size (nx,ny) instead of binning.
   H.Bushouse, 2009 Apr 24:
	Fixed a bug in the computation of biassect values for subarrays that
	include the physical overscan on the amp B/D edge of the image. In
	addition to converting the column values from 1-indexed to 0-indexed
	values, I also upgraded the logic to be equivalent to that used for
	amp A/C subarrays, which uses the biassect values from the OSCNTAB
	to set the limits. The previous logic had simply been using all
	overscan columns present in the image.
*/

int FindOverscan (WF3Info *wf3, int nx, int ny, int *overscan) {

/* arguments:
WF3Info *wf3		i: structure with all values from OSCNTAB 
int nx, ny		i: lengths of first and second image axes
int offsetx, offsety	i: Trim values from LTV1,2
*/

	extern int status;
	int row;
	TblInfo tabinfo;
	TblRow tabrow;
	int foundit;

	int cx0, cx1, tx1, tx2;

	int SameInt (int, int);
	int SameString (char *, char *);

	/* Open the CCD parameters table and find columns. */
	if (OpenOverTab (wf3->oscn.name, &tabinfo))
	    return (status);

	foundit = NO;
	*overscan = NO;
	
	/* Check each row for a match with ccdamp, chip, nx, ny,
	** binx, and biny, and get info from the matching row.  */
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    /* Read the current row into tabrow. */
	    if (ReadOverTab (&tabinfo, row, &tabrow))
		return (status);

	    if (wf3->detector == CCD_DETECTOR &&
		SameString (tabrow.ccdamp, wf3->ccdamp) &&
		SameInt (tabrow.chip, wf3->chip) &&
		SameInt (tabrow.bin[0], wf3->bin[0]) &&
		SameInt (tabrow.bin[1], wf3->bin[1])) {

		foundit = YES;
		*overscan = YES;             

		if (wf3->subarray == YES) {

		    /* We are working with a subarray... */
		    /* There is never any virtual overscan. */
		    wf3->trimy[0] = 0; wf3->trimy[1] = 0;
		    wf3->vx[0] = 0; wf3->vx[1] = 0;
		    wf3->vx[2] = 0; wf3->vx[3] = 0;
		    wf3->vy[0] = 0; wf3->vy[1] = 0;
		    wf3->vy[2] = 0; wf3->vy[3] = 0;
		    wf3->trimx[2] = 0; wf3->trimx[3] = 0;
		    wf3->biassectc[0] = 0; wf3->biassectc[1] = 0;
		    wf3->biassectd[0] = 0; wf3->biassectd[1] = 0;

		    /* Determine whether the subarray extends into the
		    ** physical overscan regions on either side of the chip */
		    tx1 = (int)(nx - wf3->offsetx);
		    cx1 = tabrow.nx - tabrow.trimx[1] - tabrow.trimx[0];
		    cx0 = (int)(tabrow.trimx[0] - wf3->offsetx);

		    if (wf3->offsetx > 0) {
		        /* Subarray starts in the first overscan region */
		        wf3->trimx[0] = (wf3->offsetx < tabrow.trimx[0]) ?
					 wf3->offsetx : tabrow.trimx[0];
			/* Check to see if it extends into second
			** overscan region */
			tx2 = (int)(nx - wf3->trimx[0]) - cx1;
			wf3->trimx[1] = (tx2 < 0) ? 0 : tx2;
			wf3->biassecta[0] = (tabrow.biassecta[0]-1 - cx0 > 0) ?
					    (tabrow.biassecta[0]-1 - cx0) : 0;
			wf3->biassecta[1] = tabrow.biassecta[1]-1 - cx0;
			wf3->biassectb[0] = 0;
			wf3->biassectb[1] = 0;
		    } else if (tx1 > cx1 && tx1 <= tabrow.nx) {

			/* then the subarray overlaps biassectb... */
			wf3->trimx[0] = 0;
			wf3->trimx[1] = tx1 - cx1;
			wf3->biassecta[0] = 0;
			wf3->biassecta[1] = 0;
			wf3->biassectb[0] = tabrow.biassecta[0]-1 - cx0;
			wf3->biassectb[1] = tabrow.biassecta[1]-1 - cx0;
			if (wf3->biassectb[0] >= nx) {
			    wf3->biassectb[0] = 0;
			    wf3->biassectb[1] = 0;
			}
			if (wf3->biassectb[1] >= nx) wf3->biassectb[1] = nx-1;

		    } else {

			/* Subarray doesn't overlap either physical
			** overscan region */
			wf3->trimx[0] = 0;
			wf3->trimx[1] = 0;
			wf3->biassecta[0] = 0;
			wf3->biassecta[1] = 0;
			wf3->biassectb[0] = 0;
			wf3->biassectb[1] = 0;                               
			*overscan = NO;             
		    }

		} else {

		    /* We are working with a full chip image... */		
		    wf3->trimx[0] = tabrow.trimx[0];
		    wf3->trimx[1] = tabrow.trimx[1];
		    wf3->trimx[2] = tabrow.trimx[2];
		    wf3->trimx[3] = tabrow.trimx[3];
		    wf3->trimy[0] = tabrow.trimy[0];
		    wf3->trimy[1] = tabrow.trimy[1];

		    /* Subtract one from table values to
		    ** conform to C array indexing... */
		    wf3->vx[0] = tabrow.vx[0]-1;
		    wf3->vx[1] = tabrow.vx[1]-1;
		    wf3->vx[2] = tabrow.vx[2]-1;
		    wf3->vx[3] = tabrow.vx[3]-1;
		    wf3->vy[0] = tabrow.vy[0]-1;
		    wf3->vy[1] = tabrow.vy[1]-1;
		    wf3->vy[2] = tabrow.vy[2]-1;
		    wf3->vy[3] = tabrow.vy[3]-1;
		    wf3->biassecta[0] = tabrow.biassecta[0] - 1;
		    wf3->biassecta[1] = tabrow.biassecta[1] - 1;
		    wf3->biassectb[0] = tabrow.biassectb[0] - 1;
		    wf3->biassectb[1] = tabrow.biassectb[1] - 1;
		    wf3->biassectc[0] = tabrow.biassectc[0] - 1;
		    wf3->biassectc[1] = tabrow.biassectc[1] - 1;
		    wf3->biassectd[0] = tabrow.biassectd[0] - 1;
		    wf3->biassectd[1] = tabrow.biassectd[1] - 1;
		}

		break;

	    } else if (wf3->detector == IR_DETECTOR &&
		SameString (tabrow.ccdamp, wf3->ccdamp) &&
		SameInt (tabrow.chip, wf3->chip) &&
		SameInt (tabrow.nx, nx) &&
		SameInt (tabrow.ny, ny)) {

		foundit = YES;
		*overscan = YES;             

		/* Copy table values to local variables */
		wf3->trimx[0] = tabrow.trimx[0];
		wf3->trimx[1] = tabrow.trimx[1];
		wf3->trimx[2] = tabrow.trimx[2];
		wf3->trimx[3] = tabrow.trimx[3];
		wf3->trimy[0] = tabrow.trimy[0];
		wf3->trimy[1] = tabrow.trimy[1];

		/* Subtract one from table values to
		** conform to C array indexing... */
		wf3->vx[0] = tabrow.vx[0]-1;
		wf3->vx[1] = tabrow.vx[1]-1;
		wf3->vx[2] = tabrow.vx[2]-1;
		wf3->vx[3] = tabrow.vx[3]-1;
		wf3->vy[0] = tabrow.vy[0]-1;
		wf3->vy[1] = tabrow.vy[1]-1;
		wf3->vy[2] = tabrow.vy[2]-1;
		wf3->vy[3] = tabrow.vy[3]-1;
		wf3->biassecta[0] = tabrow.biassecta[0] - 1;
		wf3->biassecta[1] = tabrow.biassecta[1] - 1;
		wf3->biassectb[0] = tabrow.biassectb[0] - 1;
		wf3->biassectb[1] = tabrow.biassectb[1] - 1;
		wf3->biassectc[0] = tabrow.biassectc[0] - 1;
		wf3->biassectc[1] = tabrow.biassectc[1] - 1;
		wf3->biassectd[0] = tabrow.biassectd[0] - 1;
		wf3->biassectd[1] = tabrow.biassectd[1] - 1;
	    }
	}

	if (foundit == NO) {
	    status = ROW_NOT_FOUND;
	    sprintf(MsgText, "Could not find appropriate row from OSCNTAB. ");
	    trlwarn(MsgText);
	}	

	if (CloseOverTab (&tabinfo))
	    return(status);
			
	if (wf3->verbose == YES) {
	    sprintf (MsgText, "Found trim values of: x(%d,%d,%d,%d) y(%d,%d)",
		     wf3->trimx[0], wf3->trimx[1], wf3->trimx[2], wf3->trimx[3],
		     wf3->trimy[0], wf3->trimy[1]);
	    trlmessage (MsgText);
	}

	return (status);
}


/* This routine opens the overscan parameters table, finds the columns that we
   need, and gets the total number of rows in the table. 
*/

static int OpenOverTab (char *tname, TblInfo *tabinfo) {

	extern int status;
	
	int nocol[NUMCOLS];
	int i, j, missing;
	
	char *colnames[NUMCOLS] =
	      {"CCDAMP","CCDCHIP","BINX","BINY","NX","NY","TRIMX1","TRIMX2",
		"TRIMX3","TRIMX4","TRIMY1","TRIMY2","VX1","VX2","VY1","VY2",
		"VX3","VX4","VY3","VY4",
		"BIASSECTA1","BIASSECTA2","BIASSECTB1","BIASSECTB2",
		"BIASSECTC1","BIASSECTC2","BIASSECTD1","BIASSECTD2"};

	int PrintMissingCols (int, int, int *, char **, char *, IRAFPointer);

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    sprintf (MsgText, "OSCNTAB `%s' not found.", tname);
	    trlerror (MsgText);
	    return (status = OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	for (j = 0; j < NUMCOLS; j++)
	     nocol[j] = NO;
	
	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "CCDAMP", &tabinfo->cp_amp);
	c_tbcfnd1 (tabinfo->tp, "CCDCHIP", &tabinfo->cp_chip);
	c_tbcfnd1 (tabinfo->tp, "BINX", &tabinfo->cp_binx);
	c_tbcfnd1 (tabinfo->tp, "BINY", &tabinfo->cp_biny);	
	c_tbcfnd1 (tabinfo->tp, "NX", &tabinfo->cp_nx);
	c_tbcfnd1 (tabinfo->tp, "NY", &tabinfo->cp_ny);	
	c_tbcfnd1 (tabinfo->tp, "TRIMX1", &tabinfo->cp_trimx1);
	c_tbcfnd1 (tabinfo->tp, "TRIMX2", &tabinfo->cp_trimx2);
	c_tbcfnd1 (tabinfo->tp, "TRIMX3", &tabinfo->cp_trimx3);
	c_tbcfnd1 (tabinfo->tp, "TRIMX4", &tabinfo->cp_trimx4);
	c_tbcfnd1 (tabinfo->tp, "TRIMY1", &tabinfo->cp_trimy1);
	c_tbcfnd1 (tabinfo->tp, "TRIMY2", &tabinfo->cp_trimy2);	
	c_tbcfnd1 (tabinfo->tp, "VX1", &tabinfo->cp_vx1);
	c_tbcfnd1 (tabinfo->tp, "VX2", &tabinfo->cp_vx2);
	c_tbcfnd1 (tabinfo->tp, "VX3", &tabinfo->cp_vx3);
	c_tbcfnd1 (tabinfo->tp, "VX4", &tabinfo->cp_vx4);
	c_tbcfnd1 (tabinfo->tp, "VY1", &tabinfo->cp_vy1);	
	c_tbcfnd1 (tabinfo->tp, "VY2", &tabinfo->cp_vy2);
	c_tbcfnd1 (tabinfo->tp, "VY3", &tabinfo->cp_vy3);	
	c_tbcfnd1 (tabinfo->tp, "VY4", &tabinfo->cp_vy4);
	c_tbcfnd1 (tabinfo->tp, "BIASSECTA1", &tabinfo->cp_biassecta1);
	c_tbcfnd1 (tabinfo->tp, "BIASSECTA2", &tabinfo->cp_biassecta2);
	c_tbcfnd1 (tabinfo->tp, "BIASSECTB1", &tabinfo->cp_biassectb1);
	c_tbcfnd1 (tabinfo->tp, "BIASSECTB2", &tabinfo->cp_biassectb2);
	c_tbcfnd1 (tabinfo->tp, "BIASSECTC1", &tabinfo->cp_biassectc1);
	c_tbcfnd1 (tabinfo->tp, "BIASSECTC2", &tabinfo->cp_biassectc2);
	c_tbcfnd1 (tabinfo->tp, "BIASSECTD1", &tabinfo->cp_biassectd1);
	c_tbcfnd1 (tabinfo->tp, "BIASSECTD2", &tabinfo->cp_biassectd2);

	/* Check which columns are missing */
	i=0; missing = 0;
	if (tabinfo->cp_amp == 0)        { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_chip == 0)       { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_binx == 0)       { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_biny == 0)       { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_nx == 0)         { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_ny == 0)         { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_vx1 == 0)        { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_vx2 == 0)        { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_vx3 == 0)        { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_vx4 == 0)        { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_vy1 == 0)        { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_vy2 == 0)        { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_vy3 == 0)        { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_vy4 == 0)        { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_trimx1 == 0)     { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_trimx2 == 0)     { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_trimx3 == 0)     { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_trimx4 == 0)     { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_trimy1 == 0)     { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_trimy2 == 0)     { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_biassecta1 == 0) { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_biassecta2 == 0) { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_biassectb1 == 0) { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_biassectb2 == 0) { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_biassectc1 == 0) { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_biassectc2 == 0) { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_biassectd1 == 0) { missing++; nocol[i] = YES; i++;}
	if (tabinfo->cp_biassectd2 == 0) { missing++; nocol[i] = YES; i++;}

	if (PrintMissingCols (missing, NUMCOLS, nocol, colnames, "OSCNTAB",
			      tabinfo->tp))
	    return(status);

	return (status);
}

/* This routine reads the relevant data from one row. 
	It reads the amp configuration, chip number, binning factors,
	and x/y size of the image for selecting the row appropriate
	for the data.  It also reads the trim regions, virtual overscan
	corner coordinates, and the bias section columns for both the
	trailing and leading overscan regions.
*/

static int ReadOverTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	extern int status;

	c_tbegtt (tabinfo->tp, tabinfo->cp_amp, row, tabrow->ccdamp, SZ_CBUF);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_chip, row, &tabrow->chip);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
		
	c_tbegti (tabinfo->tp, tabinfo->cp_binx, row, &tabrow->bin[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);		
	c_tbegti (tabinfo->tp, tabinfo->cp_biny, row, &tabrow->bin[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
		
	c_tbegti (tabinfo->tp, tabinfo->cp_nx, row, &tabrow->nx);
	if (c_iraferr())
	    return (status = TABLE_ERROR);	
	c_tbegti (tabinfo->tp, tabinfo->cp_ny, row, &tabrow->ny);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_trimx1, row, &tabrow->trimx[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_trimx2, row, &tabrow->trimx[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_trimx3, row, &tabrow->trimx[2]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_trimx4, row, &tabrow->trimx[3]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_trimy1, row, &tabrow->trimy[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_trimy2, row, &tabrow->trimy[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_vx1, row, &tabrow->vx[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_vx2, row, &tabrow->vx[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_vx3, row, &tabrow->vx[2]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_vx4, row, &tabrow->vx[3]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_vy1, row, &tabrow->vy[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_vy2, row, &tabrow->vy[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_vy3, row, &tabrow->vy[2]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_vy4, row, &tabrow->vy[3]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_biassecta1, row,
		  &tabrow->biassecta[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_biassecta2, row,
		  &tabrow->biassecta[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_biassectb1, row,
		  &tabrow->biassectb[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_biassectb2, row,
		  &tabrow->biassectb[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_biassectc1, row,
		  &tabrow->biassectc[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_biassectc2, row,
		  &tabrow->biassectc[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_biassectd1, row,
		  &tabrow->biassectd[0]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_biassectd2, row,
		  &tabrow->biassectd[1]);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	return (status);
}

/* This routine closes the oscn table. */

static int CloseOverTab (TblInfo *tabinfo) {

	extern int status;

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	return (status);
}
