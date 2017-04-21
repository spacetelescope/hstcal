# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "hstcalerr.h"
# include "stispht.h"


typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_sporder;		/* column descriptors */
	IRAFPointer cp_nelem;
	IRAFPointer cp_wave;
	IRAFPointer cp_intens;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	int sporder;
} TblRow;


static int OpenIntensTab (char *, TblInfo *, double *, int);
static int ReadIntensTab (TblInfo *, int, TblRow *);
static int ReadIntensArray (StisInfo6 *, TblInfo *, int, IntensArray *, int);
static int CloseIntensTab (TblInfo *);
static int FluxToNet (StisInfo6 *, IntensArray *, int);



/* This routine gets intensity/wavelength information from OSPECTAB
   (_x1d) and saves the info in the intensity information structure.
   This information is used by the optimal extraction algorithm.

   The intensity info table should contain the following:
	header parameters:
		none needed
	columns:
		SPORDER:      spectral order number
		NELEM:        actual number of elements in arrays (int)
		WAVELENGTH:   wavelengths in Angstroms (double)
		NET:          net intensity (float)
			      (corresponding to WAVELENGTH)
                or
		FLUX:         flux (float)
			      (also corresponding to WAVELENGTH)

   Rows are selected on SPORDER. If a matching row is found, the sizes
   of the wavelength and throughput arrays are gotten from column NELEM,
   memory is allocated, and WAVELENGTH and NET/FLUX are read.

   The EXPTIME keyword is also read from the extension header. It is used
   to convert the NET values back to raw counts. These raw counts are
   used by the optimal extraction algorithm to compute the Poisson
   contribution to the total variance.

   When done, the memory should be freed by calling FreeIntensity.




   Revision history:
   ----------------
   16 Sep 98  -  Implemented (I.Busko)
   29 Oct 98  -  Read EXPTIME and convert NET to counts (IB)
   10 May 99  -  Choice between FLUX or NET reference spectrum (IB)
   27 Dec 04  -  Remove declaration of InitRefTab, moved to stis.h (PEH)
   01 Jun 05  -  Initialize slit.gac_allocated to 0 (PEH)
*/

int GetIntens (StisInfo6 *sts, int sporder, IntensArray *inta) {

/* arguments:
StisInfo6 *sts        i: calibration switches and info
IntensArray *inta     o: description of inta
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit = 0;	/* true if parameters found in table */
	double exptime;		/* exposure time */
	int i;

	/* Open the intensity/wavelength table. */
	exptime = 1.0;
	if ((status = OpenIntensTab (sts->pxtab.name,
                                     &tabinfo, &exptime, sts->fflux)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadIntensTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (!SameInt (tabrow.sporder, sporder))
	        continue;

	    foundit = 1;

	    /* Get pedigree & descrip from the row. */
	    if ((status = RowPedigree (&sts->pxtab, row,
                    tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
	        return (status);

	    /* Read wavelengths and intensities into structure. */
	    if ((status = ReadIntensArray (sts, &tabinfo, row, inta, sporder)))
	        return (status);
	}

            if ((status = CloseIntensTab (&tabinfo)))
	    return (status);

	if (!foundit) {
	    printf ("Warning  Matching row not found in OSPECTAB %s\n",
			sts->pxtab.name);
	    return (ROW_NOT_FOUND);
	}

	/* Convert to counts. */
	for (i = 0;  i < inta->nelem; inta->intens[i++] *= exptime);

	return (0);
}




/* This routine opens the intensity table, finds the columns that we need,
   and gets the total number of rows in the table. It also gets the
   EXPTIME keyword from the table header.
*/

static int OpenIntensTab (char *tname, TblInfo *tabinfo, double *exptime,
                          int fflux) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    OSPECTAB `%s' not found\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, SPORDER, &tabinfo->cp_sporder);
	c_tbcfnd1 (tabinfo->tp, NELEM, &tabinfo->cp_nelem);
	c_tbcfnd1 (tabinfo->tp, WAVELENGTH, &tabinfo->cp_wave);
	if (fflux)
	    c_tbcfnd1 (tabinfo->tp, "FLUX", &tabinfo->cp_intens);
	else
	    c_tbcfnd1 (tabinfo->tp, NET, &tabinfo->cp_intens);

	if (tabinfo->cp_sporder == 0 ||
            tabinfo->cp_nelem   == 0 ||
	    tabinfo->cp_wave    == 0 ||
	    tabinfo->cp_intens  == 0) {
	    printf ("ERROR    Column not found in OSPECTAB\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	/* Get EXPTIME */
	*exptime = c_tbhgtd (tabinfo->tp, "EXPTIME");

	return (0);
}




/* This routine reads the column (sporder) used to select the correct row.
*/

static int ReadIntensTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegti (tabinfo->tp, tabinfo->cp_sporder, row, &tabrow->sporder);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}



/* This routine reads the array data from one row. The number of elements
   in the arrays is gotten, the arrays are allocated, and the wavelengths
   and intensities are read into the arrays.
*/

static int ReadIntensArray (StisInfo6 *sts, TblInfo *tabinfo, int row,
                            IntensArray *inta, int sporder) {

	int status;
	int nwave, nintens;	/* number of elements actually read */

	/* Find out how many elements there are in the intensity arrays,
	   and allocate space for the arrays to be read from the table.
	*/
	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &inta->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Allocate memory. */
	inta->wave   = (double *) calloc (inta->nelem, sizeof(double));
	inta->intens = (double *) calloc (inta->nelem, sizeof(double));
	if (inta->wave == NULL || inta->intens == NULL) {
	    CloseIntensTab (tabinfo);
	    return (OUT_OF_MEMORY);
	}
	inta->allocated = 1;

	nwave = c_tbagtd (tabinfo->tp, tabinfo->cp_wave, row,
                          inta->wave, 1, inta->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	nintens = c_tbagtd (tabinfo->tp, tabinfo->cp_intens, row,
                            inta->intens, 1, inta->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Go get data to transform flux into net counts. */
	if (sts->fflux) {
	    status = FluxToNet (sts, inta, sporder);
	    if (status != STIS_OK)
	        return (status);
	}

	if (nwave < inta->nelem || nintens < inta->nelem) {
	    c_tbtclo (tabinfo->tp);
	    free (inta->wave);
	    free (inta->intens);
	    printf ("ERROR    Not all coefficients were read from OSPECTAB\n");
	    return (TABLE_ERROR);
	}

	return (0);
}



/* This routine closes the OSPECTAB table. */

static int CloseIntensTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}



/*  Converts the FLUX array to net counts.

    This function performs the inverse operation in AbsFlux6: given the
    flux values, convert then back to net counts/s.

    This function DOES NOT include the time-dependent sensitivity
    correction, since small errors in the derived net counts have
    negligible impact in the optimal extraction algorithm.
*/

static int FluxToNet (StisInfo6 *sts, IntensArray *inta, int sporder) {

	/* This is used to store information from the fflux file in a
           form suitable for the reference file input routines.
        */
	StisInfo6 fsts;
	ApInfo slit;
	PhotInfo phot;

	IODescPtr im;
	Hdr phdr;
	double photfactor, throughput, response, dispersion;
	double atodgain, readnoise;
	float correction;
	int i, dispc, helc, status;
	int abs_starti, thr_starti;
	int dummy;

	void FreePhot6 (PhotInfo *);
	void FreeThroughput6 (ApInfo *);
	int GetAbsPhot6 (StisInfo6 *, int, PhotInfo *, int, int *);
	int GetApDes6 (StisInfo6 *, ApInfo *);
	int GetApThr6 (StisInfo6 *, ApInfo *);
	int Get_KeyD (Hdr *, char *, int, double, double *);
	int Get_KeyS (Hdr *, char *, int, char *, char *, int);
	int GetSwitch (Hdr *, char *, int *);
	double interp1d (double, double *, double *, int, int *);
	void StisInit6 (StisInfo6 *);

	photfactor = H_PLANCK * C_LIGHT / HST_AREA;

	/* Initialize local data structures. */
	StisInit6 (&fsts);
        InitRefTab (&fsts.phottab);
        InitRefTab (&fsts.apertab);
        InitRefTab (&fsts.apdestab);
	slit.allocated  = 0;
	slit.gac_allocated  = 0;
	phot.allocated  = 0;
	phot.pcorr      = NULL;

	/* Handling the primary header here is not efficient. But keeps
           this new code manageable since everything new is added at a
           single point. In the future we may move this to outside the
           main loop and pass the necessary values as part of the sts
           structure.
        */
	initHdr (&phdr);
	im = openInputImage (sts->pxtab.name, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (im, &phdr);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);

	/* Abort if both helcorr and dispcorr weren't performed.
           The criterion is: if a keyword is set to either COMPLETE
           or PERFORM, we assume that the operation was performed.
           This is because UpdHdrSwitch in Do1Dx only updates the
           keywords to COMPLETE if they are set to PERFORM in the
           input file.  (note:  UpdHdrSwitch is no longer used)
        */
	if ((status = GetSwitch (&phdr, "DISPCORR", &dispc)))
	    return (status);
	if ((status = GetSwitch (&phdr, "HELCORR", &helc)))
	    return (status);
	if (!((dispc == PERFORM || dispc == COMPLETE) &&
              (helc  == PERFORM || helc  == COMPLETE))) {
	    printf ("ERROR    No DISPCORR/HELCORR in fflux file.\n");
	    return (ERROR_RETURN);
	}

	/* Read header keywords. */
	if ((status = Get_KeyD (&phdr, "READNSE", 1, 0., &readnoise)))
	    return (status);
	if ((status = Get_KeyD (&phdr, "ATODGAIN", 1, 1., &atodgain)))
	    return (status);
	if ((status = Get_KeyS (&phdr, "PHOTTAB", FATAL, "",
                                fsts.phottab.name, STIS_LINE)))
	    return (status);
	if ((status = Get_KeyS (&phdr, "APDESTAB", FATAL, "",
                                fsts.apdestab.name, STIS_LINE)))
	    return (status);
	if ((status = Get_KeyS (&phdr, "APERTAB", FATAL, "",
                                fsts.apertab.name, STIS_LINE)))
	    return (status);

	/* Copy stuff from primary data structure into local one. */
	fsts.x1d_o    = sts->x1d_o;
	fsts.dispcorr = sts->dispcorr;
	fsts.fluxcorr = sts->fluxcorr;
	fsts.pctcorr  = sts->pctcorr;
	fsts.cenwave  = sts->cenwave;
	strcpy (fsts.opt_elem, sts->opt_elem);
	strcpy (fsts.aperture, sts->aperture);

	/* Read the required reference info. */
	dummy = 0;
	if ((status = GetAbsPhot6 (&fsts, sporder, &phot, 0, &dummy)))
	    return (status);
	if ((status = GetApDes6 (&fsts, &slit)))
	    return (status);
        if ((status = GetApThr6 (&fsts, &slit)))
	    return (status);

	abs_starti = 1;				/* initial values */
	thr_starti = 1;

	/* Loop over flux array. */
	for (i = 0;  i < inta->nelem;  i++) {
	    response   = interp1d (inta->wave[i], phot.wl, phot.thru,
                                   phot.nelem, &abs_starti);
	    throughput = interp1d (inta->wave[i], slit.wl, slit.thr,
                                   slit.nelem, &thr_starti);
	    if (i > 0)
	        dispersion = inta->wave[i] - inta->wave[i-1];
	    else
	        dispersion = inta->wave[1] - inta->wave[0];

	    /* This check is provisional; final version awaits IS's words. */
	    if (response   <= 0.0 ||
	        dispersion <= 0.0 ||
	        throughput <= 0.0) {
	        printf ("ERROR    Error in fflux file contents.\n");
	        return (ERROR_RETURN);
	    }

	    correction = (float) (photfactor / (response * throughput *
                         inta->wave[i] * dispersion * atodgain *
                         CM_PER_ANGSTROM));

	    inta->intens[i] = inta->intens[i] / correction;
	}

	FreeThroughput6 (&slit);
	FreePhot6 (&phot);

	freeHdr (&phdr);
	return STIS_OK;
}
