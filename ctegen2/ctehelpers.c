#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "hstcal_memory.h"
#include "hstcal.h"
#include "wf3.h" //need to remove this dependency
#include "ctegen2.h"

void initCTEParamsFast(CTEParamsFast * pars, const unsigned _nTraps,
        const unsigned _nRows, const unsigned _nColumns,
        const unsigned _nScaleTableColumns, const unsigned _maxThreads)
{
    pars->maxThreads = _maxThreads;
    pars->verbose = False;
    pars->nRows = _nRows;
    pars->nColumns = _nColumns;
    pars->nTraps = _nTraps;
    pars->nScaleTableColumns = _nScaleTableColumns;

    pars->noise_mit=0;
    pars->thresh=0;
    pars->cte_date0=0;
    pars->cte_date1=0;
    pars->cte_traps=0;
    pars->cte_len=0;
    pars->rn_amp=0;
    pars->n_forward=0;
    pars->n_par=0;
    pars->scale_frac=0; /*will be updated during routine run*/
    pars->fix_rocr = 0;

    pars->nRowsPerFullFrame = 0;
    pars->nColumnsPerFullFrame = 0;
    pars->nColumnsPerChip = 0;
    pars->nRowsPerChip = 0;
    pars->nColumnsPerQuad = 0;
    pars->nRowsPerQuad = 0;

    pars->isSubarray = False;
    pars->chip = 0;
    pars->refAndIamgeBinsIdenticle = True;
    pars->rowOffset = 0; //from begining of chip (for WFC3 this includes overscan, for ACS this does not, i.e. it has already been trimmed)
    pars->columnOffset = 0; //from begining of chip (for WFC3 this includes overscan, for ACS this does not, i.e. it has already been trimmed)
    pars->imageRowsStart = 0;
    pars->imageRowsEnd = 0;
    pars->postscanWidth = 0;
    pars->prescanWidth = 0;
    pars->parallelOverscanWidth = 0;

    //alignment relative to old RAZ format, both chips side by side ordered CDAB
    pars->razColumnOffset = 0;

    {unsigned i;
    for (i = 0; i < 2; ++i)
    {
        pars->imageColumnsStart[i] = 0;
        pars->imageColumnsEnd[i] = 0;
        pars->hasPrescan[i] = False;
        pars->hasPostscan[i] = False;
        pars->quadExists[i] = False;
    }}

    pars->iz_data = NULL;
    pars->wcol_data = NULL;
    pars->scale512 = NULL;
    pars->scale1024 = NULL;
    pars->scale1536 = NULL;
    pars->scale2048 = NULL;
    pars->qlevq_data = NULL;
    pars->dpdew_data = NULL;

    pars->rprof = NULL; /* differential trail profile as image */
    pars->cprof = NULL; /* cummulative trail profile as image */

    *pars->cte_name='\0';
    *pars->cte_ver='\0';
}

int allocateCTEParamsFast(CTEParamsFast * pars)
{
    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    void * tmp = NULL;
    tmp = newAndZero((void*)&pars->iz_data, pars->nScaleTableColumns, sizeof(*pars->iz_data));
    addPtr(&ptrReg, tmp, &free);
    if (!tmp)
    {
        freeOnExit(&ptrReg);
        trlerror ("Out of memory.\n");
        return OUT_OF_MEMORY;
    }
    tmp = newAndZero((void*)&pars->scale512, pars->nScaleTableColumns, sizeof(*pars->scale512));
    addPtr(&ptrReg, tmp, &free);
    if (!tmp)
    {
        freeOnExit(&ptrReg);
        trlerror ("Out of memory.\n");
        return OUT_OF_MEMORY;
    }
    tmp = newAndZero((void*)&pars->scale1024, pars->nScaleTableColumns, sizeof(*pars->scale1024));
    addPtr(&ptrReg, tmp, &free);
    if (!tmp)
    {
        freeOnExit(&ptrReg);
        trlerror ("Out of memory.\n");
        return OUT_OF_MEMORY;
    }
    tmp = newAndZero((void*)&pars->scale1536, pars->nScaleTableColumns, sizeof(*pars->scale1536));
    addPtr(&ptrReg, tmp, &free);
    if (!tmp)
    {
        freeOnExit(&ptrReg);
        trlerror ("Out of memory.\n");
        return OUT_OF_MEMORY;
    }
    tmp = newAndZero((void*)&pars->scale2048, pars->nScaleTableColumns, sizeof(*pars->scale2048));
    addPtr(&ptrReg, tmp, &free);
    if (!tmp)
    {
        freeOnExit(&ptrReg);
        trlerror ("Out of memory.\n");
        return OUT_OF_MEMORY;
    }
    tmp = newAndZero((void*)&pars->wcol_data, pars->nTraps, sizeof(*pars->wcol_data));
    addPtr(&ptrReg, tmp, &free);
    if (!tmp)
    {
        freeOnExit(&ptrReg);
        trlerror ("Out of memory.\n");
        return OUT_OF_MEMORY;
    }
    tmp = newAndZero((void*)&pars->qlevq_data, pars->nTraps, sizeof(*pars->qlevq_data));
    addPtr(&ptrReg, tmp, &free);
    if (!tmp)
    {
        freeOnExit(&ptrReg);
        trlerror ("Out of memory.\n");
        return OUT_OF_MEMORY;
    }
    tmp = newAndZero((void*)& pars->dpdew_data, pars->nTraps, sizeof(*pars->dpdew_data));
    addPtr(&ptrReg, tmp, &free);
    if (!tmp)
    {
        freeOnExit(&ptrReg);
        trlerror ("Out of memory.\n");
        return OUT_OF_MEMORY;
    }

    freeReg(&ptrReg);
    return 0;
}

void freeCTEParamsFast(CTEParamsFast * pars)
{
    delete((void*)&pars->iz_data);
    delete((void*)&pars->wcol_data);
    delete((void*)&pars->scale512);
    delete((void*)&pars->scale1024);
    delete((void*)&pars->scale1536);
    delete((void*)&pars->scale2048);
    delete((void*)&pars->qlevq_data);
    delete((void*)&pars->dpdew_data);

    freeFloatHdrData(pars->rprof);
    freeFloatHdrData(pars->cprof);
    delete((void*)&pars->rprof);
    delete((void*)&pars->cprof);
}

/************ HELPER SUBROUTINES ****************************/
/*
MLS 2015: read in the CTE parameters from the PCTETAB file

    Jan 19, 2016: MLS  Updated to check for existence of PCTENSMD in
                  raw science header

MDD Sept 2023: Implementation to support a PCTETAB which is an
    update to the Generation 2 CTE correction.  This file will
    now contain both parallel(y) and serial(x) CTE information.
    Extensions 1-4 contain the parallel CTE correction data.  The
    remaining extensions contain the serial CTE correction, one
    full set of corrections (QPROF, SCLBYCOL, RPROF, and CPROF)
    per amplifier.
*/

int loadPCTETAB (char *filename, CTEParamsFast *pars, int extn, Bool skipLoadPrimary) {
/* Read the cte parameters from the reference table PCTETAB

   These are taken from the PCTETAB global header (aka primary header):
   CTE_NAME - name of cte algorithm
   CTE_VER  - version number of cte algorithm
   PCTERNOI - readnoise amplitude and clipping level
   PCTENSMD - readnoise mitigation algorithm
   PCTETRSH - over-subtraction threshold
   FIXROCR  - account for cosmic ray over-subtraction?

   These values are now amp-dependent for the *serial* CTE correction.  This
   means that they are read from the QPROF extension header of the amp being 
   processed.  In contrast, there is only one set for the *parallel* CTE 
   correction which applies to all amps.  The parallel values are also read 
   from the QPROF extension header where the parallel correction is defined
   as discussed below.  
   CTEDATE0 - date of instrument installation in HST, in fractional years (MJD)
   CTEDATE1 - reference date of CTE model pinning, in fractional years (MJD)
   PCTENFOR - number of iterations used in CTE forward modeling
   PCTENPAR - number of iterations used in the parallel transfer
   PCTETLEN - max length of CTE trail

   The parallel/serial PCTETAB which supports the latest CTE algorithm
      CTE_VER = '3.0     '
      CTE_NAME= 'Par/Serial PixelCTE 2023'

   The updated PCTETAB has a primary HDU and 20 extensions.  Extensions 1-4 
   apply to the parallel CTE and the remaining extensions apply to the serial CTE.

No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU      75   ()
  1  QPROF         1 BinTableHDU     21   9999R x 4C   [I, J, E, 20A]
  2  SCLBYCOL      1 BinTableHDU     25   8192R x 6C   [I, E, E, E, E, 20A]
  3  RPROF         1 ImageHDU        21   (9999, 100)   float32
  4  CPROF         1 ImageHDU        21   (9999, 100)   float32
  5  QPROF         2 BinTableHDU     22   9999R x 4C   [I, J, E, 20A]
  6  SCLBYCOL      2 BinTableHDU     26   8192R x 6C   [I, E, E, E, E, 20A]
  7  RPROF         2 ImageHDU        22   (9999, 100)   float32
  8  CPROF         2 ImageHDU        22   (9999, 100)   float32
  9  QPROF         3 BinTableHDU     22   9999R x 4C   [I, J, E, 20A]
 10  SCLBYCOL      3 BinTableHDU     26   8192R x 6C   [I, E, E, E, E, 20A]
 11  RPROF         3 ImageHDU        22   (9999, 100)   float32
 12  CPROF         3 ImageHDU        22   (9999, 100)   float32
 13  QPROF         4 BinTableHDU     22   9999R x 4C   [I, J, E, 20A]
 14  SCLBYCOL      4 BinTableHDU     26   8192R x 6C   [I, E, E, E, E, 20A]
 15  RPROF         4 ImageHDU        22   (9999, 100)   float32
 16  CPROF         4 ImageHDU        22   (9999, 100)   float32
 17  QPROF         5 BinTableHDU     22   9999R x 4C   [I, J, E, 20A]
 18  SCLBYCOL      5 BinTableHDU     26   8192R x 6C   [I, E, E, E, E, 20A]
 19  RPROF         5 ImageHDU        22   (9999, 100)   float32
 20  CPROF         5 ImageHDU        22   (9999, 100)   float32

   The "extn" parameter indicates the starting extension for the QPROF table to be read
   as this routine will be called multiple times - first for the parallel CTE correction
   and then for the multiple sets (one set per amp) of serial CTE corrections. The serial
   CTE is performed first.  A "set" here means QPROF, SCLBYCOL, RPROF, and CPROF extensions.
      Parallel: starting extension 1 (extensions 1 - 4)
      Serial Amp A: starting extension 5 (extensions 5 - 8)
      Serial Amp B: starting extension 9 (extensions 9 - 12)
      Serial Amp C: starting extension 13 (extensions 13 - 16)
      Serial Amp D: starting extension 17 (extensions 17 - 20)
*/

    extern int status; /* variable for return status */

    /* VARIABLE FOR FILENAME + EXTENSION NUMBER. */
    char filename_wext[strlen(filename) + 4];

    /* NAMES OF DATA COLUMNS WE WANT FROM THE FILE, DATA WILL BE STORED IN THE PARS STRUCTURE */
    /* QPROF table - the last column is a description string */
    const char wcol[] = "W";
    const char qlevq[] = "QLEV_Q";
    const char dpdew[] = "DPDE_W";

    /* SCLBYCOL table - the last column is a description string */
    const char iz[]    = "IZ";
    const char sens512[]= "SENS_0512";
    const char sens1024[] = "SENS_1024";
    const char sens1536[] = "SENS_1536";
    const char sens2048[] = "SENS_2048";

    /* Read in the primary header keywords */
    if (!skipLoadPrimary)
    {
        sprintf(MsgText, "(ctehelpers) Reading PRIMARY.  cte_name: %s skipLoadPrimary: %c\n", pars->cte_name, skipLoadPrimary);  
        /* HSTIO VARIABLES */
        Hdr hdr_ptr;
        initHdr(&hdr_ptr);
        /* LOAD PRIMARY HEADER */
        if (LoadHdr(filename, &hdr_ptr)) {
            sprintf(MsgText,"(pctecorr) Error loading header from %s",filename);
            cteerror(MsgText);
            status = OPEN_FAILED;
            return status;
        }

		/* GET CTE_NAME KEYWORD */
		if (GetKeyStr (&hdr_ptr, "CTE_NAME", NO_DEFAULT, "", pars->cte_name, SZ_CBUF)) {
			cteerror("(pctecorr) Error reading CTE_NAME keyword from PCTETAB");
			status = KEYWORD_MISSING;
			return status;
		}
		sprintf(MsgText,"\nCTE_NAME: %s",pars->cte_name);
		trlmessage(MsgText);

		/* GET VERSION NUMBER  */
		if (GetKeyStr(&hdr_ptr, "CTE_VER", NO_DEFAULT, "", pars->cte_ver, SZ_CBUF)) {
			cteerror("(pctecorr) Error reading CTE_VER keyword from PCTETAB");
			status = KEYWORD_MISSING;
			return status;
		}
		sprintf(MsgText,"CTE_VER: %s",pars->cte_ver);
		trlmessage(MsgText);

		/* GET READ NOISE CLIPPING LEVEL */
		if (GetKeyDbl(&hdr_ptr, "PCTERNOI", NO_DEFAULT, -999, &pars->rn_amp)) {
			cteerror("(pctecorr) Error reading PCTERNOI keyword from PCTETAB");
			status = KEYWORD_MISSING;
			return status;
		}
		sprintf(MsgText,"PCTERNOI: %f",pars->rn_amp);
		trlmessage(MsgText);

		/* GET READ NOISE MITIGATION ALGORITHM*/
		if (GetKeyInt(&hdr_ptr, "PCTENSMD", NO_DEFAULT, -999, &pars->noise_mit)) {
			cteerror("(pctecorr) Error reading PCTENSMD keyword from PCTETAB");
			status = KEYWORD_MISSING;
			return status;
		}
		sprintf(MsgText,"PCTENSMD: %d",pars->noise_mit);
		trlmessage(MsgText);

		/* GET OVER SUBTRACTION THRESHOLD */
		if (GetKeyDbl(&hdr_ptr, "PCTETRSH", NO_DEFAULT, -999, &pars->thresh)) {
			cteerror("(pctecorr) Error readsng PCTETRSH keyword from PCTETAB");
			status = KEYWORD_MISSING;
			return status;
		}
		sprintf(MsgText,"PCTETRSH: %g",pars->thresh);
		trlmessage(MsgText);

		/* FIX THE READOUT CR'S? */
		if (GetKeyInt(&hdr_ptr, "FIXROCR", NO_DEFAULT, -999, &pars->fix_rocr)){
			cteerror("(pctecorr) Error reading FIXROCR keyword from PCTETAB");
			status = KEYWORD_MISSING;
			return status;
		}
		sprintf(MsgText,"FIXROCR: %d",pars->fix_rocr);
		trlmessage(MsgText);

        /* DONE READING STUFF FROM THE PRIMARY HEADER */
        freeHdr(&hdr_ptr);
    }

    /*
     Read in the remaining keywords necessary for proper processing and are
     amp-dependent from the associated "extn". 

     The variable extn contains the number of the FITS extension which is the starting
     extension for the "set" of qprof/sclbycol/rprof/cprof data.  The values for the
     starting extension of a new set are: 5, 9, 13, 17, and 1.  For example, the starting
     serial CTE correct for Amp A as seen from the description of the PCTETAB earlier
     in this file:
     Extension EXTNAME     EXTVER
     5         QPROF         2
     6         SCLBYCOL      2
     7         RPROF         2
     8         CPROF         2
    */

    /*****************************************************************************/
    /* READ DATA FROM THE SPECIFIED QPROF TABLE - EXTENSIONS 5, 9, 13, 17, and 1 */
    sprintf(filename_wext, "%s[%i]", filename, extn);

    sprintf(MsgText,"Opening %s to read QPROF table.",filename_wext);
    ctemessage(MsgText);
    /* OPEN PARAMETERS FILE TO QPROF EXTENSION */
    IRAFPointer tbl_ptr = c_tbtopn(filename_wext, IRAF_READ_ONLY, 0); // xtables table pointer
    if (c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error opening %s with xtables",filename_wext);
        cteerror(MsgText);
        status = OPEN_FAILED;
        c_tbtclo(tbl_ptr);
        return status;
    }

    /* 
       First get the keywords necessary for processing from the extension header
       All the extension headers pertaining to an amp correction have the same 
       values for these keywords, so they only have to be read from the QPROF 
       extension.
    */

	/* GET DATE SERIAL CTE BECAME RELEVANT */
    pars->cte_date0 = c_tbhgtd(tbl_ptr, "CTEDATE0");
    if (status = c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error reading CTEDATE0 from extension %d.", extn);
        cteerror(MsgText);
        c_tbtclo(tbl_ptr);
        return status;
    }
	sprintf(MsgText,"CTEDATE0: %g",pars->cte_date0);
	trlmessage(MsgText);

	/* GET REFRENCE DATE OF CTE MODEL PINNING */
    pars->cte_date1 = c_tbhgtd(tbl_ptr, "CTEDATE1");
    if (status = c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error reading CTEDATE1 from extension %d.", extn);
        cteerror(MsgText);
        c_tbtclo(tbl_ptr);
        return status;
    }
	sprintf(MsgText,"CTEDATE1: %g",pars->cte_date1);
	trlmessage(MsgText);

	/* READ MAX LENGTH OF CTE TRAIL */
    pars->cte_len = c_tbhgti(tbl_ptr, "PCTETLEN");
    if (status = c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error reading PCTETLEN from extension %d.", extn);
        cteerror(MsgText);
        c_tbtclo(tbl_ptr);
        return status;
    }
	sprintf(MsgText,"PCTETLEN: %d",pars->cte_len);
	trlmessage(MsgText);

	/* GET NUMBER OF ITERATIONS USED IN FORWARD MODEL */
    pars->n_forward = c_tbhgti(tbl_ptr, "PCTENFOR");
    if (status = c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error reading PCTENFOR from extension %d.", extn);
        cteerror(MsgText);
        c_tbtclo(tbl_ptr);
        return status;
    }
	sprintf(MsgText,"PCTENFOR: %d",pars->n_forward);
	trlmessage(MsgText);

	/* GET NUMBER OF ITERATIONS USED IN TRANSFER*/
    pars->n_par = c_tbhgti(tbl_ptr, "PCTENPAR");
    if (status = c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error reading PCTENPAR from extension %d.", extn);
        cteerror(MsgText);
        c_tbtclo(tbl_ptr);
        return status;
    }
	sprintf(MsgText,"PCTENPAR: %d",pars->n_par);
	trlmessage(MsgText);

    /* READ DATA FROM TABLE */

    /* get column pointer for w */
    IRAFPointer w_ptr = c_tbcfnd1_retPtr(tbl_ptr, wcol);
    if (c_iraferr() || !w_ptr) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",wcol);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    /* GET COLUMN POINTER FOR QLEVQ */
    IRAFPointer qlevq_ptr = c_tbcfnd1_retPtr(tbl_ptr, qlevq);
    if (c_iraferr() || !qlevq_ptr) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",qlevq);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    /* GET COLUMN POINTER FOR DPDEW */
    IRAFPointer dpdew_ptr = c_tbcfnd1_retPtr(tbl_ptr, dpdew);
    if (c_iraferr() || !dpdew_ptr) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",dpdew);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    // LOOP OVER TABLE ROWS UP TO SIZE TRAPS
    int ctraps = 0; // actual usable traps, i.e. see if more traps were added to reference file
    {unsigned j;
    for (j = 0; j < pars->nTraps; ++j) {

        /* GET W FROM THIS ROW */
    	pars->wcol_data[j] = c_tbeGetInt(tbl_ptr, w_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, wcol);
            cteerror(MsgText);
            status = TABLE_ERROR;
            return status;
        }

        /* GET QLEVQ FROM THIS ROW */
        pars->qlevq_data[j] = c_tbeGetDouble(tbl_ptr, qlevq_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, qlevq);
            cteerror(MsgText);
            status = TABLE_ERROR;
            return status;
        }

        if (pars->qlevq_data[j] < 999999.)
            ctraps+=1;

        /* GET DPDEW FROM THIS ROW */
        pars->dpdew_data[j] = c_tbeGetDouble(tbl_ptr, dpdew_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, dpdew);
            cteerror(MsgText);
            status = TABLE_ERROR;
            return status;
        }
        if (ctraps > pars->nTraps){
            sprintf(MsgText,"More TRAPS in reference file than available, update TRAPS: %i -> %i",pars->nTraps,(int)ctraps);
            trlmessage(MsgText);
        }
    }}

    /*IF CTRAPS EVER OVERFLOWS INT THIS NEEDS TO BE CHANGED*/
    pars->cte_traps = ctraps;

    /*
    sprintf(MsgText,"(pctecorr) data check for PCTETAB QPROF, row %i, %i\t%g\t%g\ttraps=%i\n",20,
            pars->wcol_data[19],pars->qlevq_data[19], pars->dpdew_data[19], pars->cte_traps);
    trlmessage(MsgText);
    */

    /* CLOSE CTE PARAMETERS FILE FOR QPROF EXTENSION */
    c_tbtClose((void*)&tbl_ptr);
    assert(!tbl_ptr);

    /****************************************************************************/
    /****************************************************************************/
    /* READ CTE SCALING DATA FROM THE SPECIFIED SCLBYCOL TABLE - EXTENSIONS 6, 10, 14, 18, and 2 */
    sprintf(filename_wext, "%s[%i]", filename, extn + 1);

    sprintf(MsgText,"Opening %s to read SCLBYCOL table.",filename_wext);
    ctemessage(MsgText);
    tbl_ptr = c_tbtopn(filename_wext, IRAF_READ_ONLY, 0);
    if (c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error opening %s with xtables",filename_wext);
        cteerror(MsgText);
        status = OPEN_FAILED;
        c_tbtclo(tbl_ptr);
        return status;
    }

    /*get column pointer for iz column*/
    IRAFPointer iz_ptr = c_tbcfnd1_retPtr(tbl_ptr, iz);
    if (c_iraferr() || iz_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",iz);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    /* get column pointer for sens512 */
    IRAFPointer sens512_ptr = c_tbcfnd1_retPtr(tbl_ptr, sens512);
    if (c_iraferr() || w_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",sens512);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    /* get column pointer for sens1024 */
    IRAFPointer sens1024_ptr = c_tbcfnd1_retPtr(tbl_ptr, sens1024);
    if (c_iraferr() || w_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",sens1024);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }
    /* get column pointer for sens1536 */
    IRAFPointer sens1536_ptr = c_tbcfnd1_retPtr(tbl_ptr, sens1536);
    if (c_iraferr() || w_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",sens1536);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }
    /* get column pointer for sens2048 */
    IRAFPointer sens2048_ptr = c_tbcfnd1_retPtr(tbl_ptr, sens2048);
    if (c_iraferr() || w_ptr == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",sens2048);
        cteerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
    }

    /* read data from table */
    /* loop over table rows */
    {unsigned j;
    for (j = 0; j < pars->nScaleTableColumns; ++j)
    {
        /* get trap from this row */
        pars->iz_data[j] = c_tbeGetInt(tbl_ptr, iz_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, iz);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }
        pars->scale512[j] = c_tbeGetDouble(tbl_ptr, sens512_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, sens512);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }

        pars->scale1024[j] = c_tbeGetDouble(tbl_ptr, sens1024_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, sens1024);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }
        pars->scale1536[j] = c_tbeGetDouble(tbl_ptr, sens1536_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, sens1536);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }
        pars->scale2048[j] = c_tbeGetDouble(tbl_ptr, sens2048_ptr, j+1);
        if (c_iraferr()) {
            sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, sens2048);
            cteerror(MsgText);
            return (status = TABLE_ERROR);
        }
    }}
    // for testing
    /*{
    unsigned j = pars->nColumns;
    sprintf(MsgText,"(pctecorr) data check for PCTETAB SCLBYCOL row %d, %d %g\t%g\t%g\t%g\ntotal traps = %i",
            j,pars->iz_data[j-1],pars->scale512[j-1],pars->scale1024[j-1],pars->scale1536[j-1],pars->scale2048[j-1],pars->cte_traps);
    trlmessage(MsgText);
    }
    */

    /* close CTE parameters file for SCLBYCOL extension */
    c_tbtClose((void*)&tbl_ptr);
    assert(!tbl_ptr);

    /*********************************************************************************/
    /* extensions 7, 11, 15, 19, and 3 - RPROF : differential trail profile as image */
    /* Images are read by EXTNAM and EXTVER and not extension number                 */

    /* Determine the extension version number for the RPROF/CPROF images */
    int extver = (extn / 4) + 1;

    sprintf(MsgText,"Reading in image from RPROF EXTVER %d.", extver);
    ctemessage(MsgText);

    /* Get the coefficient images from the PCTETAB */
    pars->rprof = malloc(sizeof(*pars->rprof));
    if (pars->rprof == NULL){
        sprintf (MsgText, "Can't allocate memory for RPROF ref data");
        trlerror (MsgText);
        return (status = 1);
    }
   
    initFloatHdrData(pars->rprof);
    pars->rprof->data.storageOrder = COLUMNMAJOR;
    if (getFloatHD (filename, "RPROF", extver, pars->rprof)){
        return (status=1);
    }


    /*********************************************************************************/
    /* extensions 8, 12, 16, 20, and 4 -  CPROF : cummulative trail profile as image */
    /* This image has the same EXTVER as the RPROF.                                  */

    sprintf(MsgText,"Reading in image from CPROF EXTVER %d.", extver);
    ctemessage(MsgText);

    pars->cprof  = malloc(sizeof(*pars->cprof));
    if (pars->cprof == NULL){
        sprintf (MsgText, "Can't allocate memory for CPROF ref data");
        trlerror (MsgText);
        return (status = 1);
    }

    /* Get the coefficient images from the PCTETAB */
    initFloatHdrData (pars->cprof);
    pars->cprof->data.storageOrder = COLUMNMAJOR;
    if (getFloatHD (filename, "CPROF", extver, pars->cprof)){
        return (status=1);
    }

    return(status);
}

/*
  Some CTE parameters  will likely end up in the image headers so users can tune them.
  This will first check whether these parameters are set in the image header. If they
  are, then the values found there will be used instead of the values read from
  the PCTETAB. If the values are not set they will be populated with the values
  from the PCTETAB, namely these:

        'CTE_NAME':'Par/Serial pixelCTE 2023', #name of cte algorithm
        'CTE_VER':'3.0' ,  #version number of algorithm
        'CTEDATE0':nnnnn.0, #date of instrument installation in HST in MJD
        'CTEDATE1':nnnnn.0, #reference date of cte model pinning in MJd
        'PCTETLEN':100, #max length of CTE trail
        'PCTERNOI':n.nn, #read noise amplitude, clipping limit
        'PCTENFOR':n ,#number of iterations used in cte forward modeling
        'PCTENPAR':n ,#number of iterations used in parallel transfer
        'PCTENSMD':n ,#read noise mitigation algorithm
        'PCTETRSH':-10.0 ,#over subtraction threshold, always use reference value
        'FIXROCR' : 1, #set to 1 for true, fix the readout cr's

 */

int populateImageFileWithCTEKeywordValues(SingleGroup *group, CTEParamsFast *pars)
{
    extern int status;

    if ((status = updateKeyOrAddAsHistKeyStr(group->globalhdr,"CTE_NAME", pars->cte_name, "CTE algorithm name")))
    {
        trlerror("(pctecorr) failed to update (or add as history) CTE_NAME keyword in image header");
        return status;
    }

    if ((status = updateKeyOrAddAsHistKeyStr(group->globalhdr,"CTE_VER", pars->cte_ver, "CTE algorithm version")))
    {
        trlerror("(pctecorr) failed to update (or add as history) CTE_VER keyword in image header");
        return status;
    }

    if ((status = updateKeyOrAddAsHistKeyDouble(group->globalhdr, "CTEDATE0", pars->cte_date0,"Date of instrument installation")))
    {
        trlerror("(pctecorr) failed to update (or add as history) CTEDATE0 keyword in header");
        return status;
    }

    if ((status = updateKeyOrAddAsHistKeyDouble(group->globalhdr, "PCTETRSH", pars->thresh,"CTE over subtraction threshold")))
    {
        trlerror("(pctecorr) failed to update (or add as history) PCTETRSH keyword in header");
        return status;
    }


    if ((status = updateKeyOrAddAsHistKeyDouble(group->globalhdr, "CTEDATE1", pars->cte_date1, "Date of CTE model pinning")))
    {
        trlerror("(pctecorr) failed to update (or add as history) CTEDATE1 keyword in header");
        return status;
    }

    if ((status = updateKeyOrAddAsHistKeyDouble(group->globalhdr, "PCTEFRAC", pars->scale_frac, "CTE scaling factor")))
    {
        trlerror("(pctecorr) failed to update (or add as history) PCTEFRAC to image header");
        return status;
    }

    if ((status = updateKeyOrAddAsHistKeyInt(group->globalhdr,"PCTETLEN",pars->cte_len,"max length of CTE trail")))
    {
        trlerror("(pctecorr) failed to update (or add as history) PCTETLEN in header");
        return status;
    }

    // The ACS CTE correction doesn't use this amp independent value, it uses the amp dependent values read from the CCDTAB.
    // The values used are written the the resultant image file as history keywords (pcteHistory()).
    /*if ((status = updateKeyOrAddAsHistKeyDouble(group->globalhdr, "PCTERNOI", pars->rn_amp,"read noise amp clip limit")))
    {
        trlerror("(pctecorr) failed to update (or add as history) PCTERNOI in header");
        return status;
    }*/

    if ((status = updateKeyOrAddAsHistKeyInt(group->globalhdr, "PCTENFOR",pars->n_forward,"Number of iter in forward model")))
    {
        trlerror("(pctecorr) failed to update (or add as history) PCTENFOR in header");
        return status;
    }

    if ((status = updateKeyOrAddAsHistKeyInt(group->globalhdr, "PCTENPAR",pars->n_par,"Number of iter in parallel transfer")))
    {
        trlerror("(pctecorr) failed to update (or add as history) PCTENPAR in header");
        return status;
    }

    if ((status = updateKeyOrAddAsHistKeyInt(group->globalhdr,"FIXROCR",pars->fix_rocr,"fix readout cosmic rays")))
    {
        trlerror("(pctecorr) failed to update (or add as history) FIXROCR keyword in header");
        return status;
    }

    return HSTCAL_OK;
}

int getCTEParsFromImageHeader(SingleGroup *group, CTEParamsFast *pars) {

    extern int status;

    double rn_amp = 0;
    int cte_len = 0;
    int n_forward = 0;
    int n_par = 0;
    int fix_rocr = 0;
    int noise_mit = 0;

    trlmessage("Begin attempt to read CTE keywords from input image header. Missing CTE keywords");
    trlmessage("from input image header are not a fatal error as the values from PCTETAB will be used.");
    int tempStatus = HSTCAL_OK;
    /*check the PCTENSMD keyword in the header*/
    tempStatus = GetKeyInt(group->globalhdr, "PCTENSMD", NO_DEFAULT, -999, &noise_mit);
    if (tempStatus == HSTCAL_OK)
        pars->noise_mit = noise_mit;
    else if (tempStatus != KEYWORD_MISSING)
    {
        trlerror("(pctecorr) Error reading PCTENSMD keyword from header");
        return (status = tempStatus);
    }
    else if (tempStatus == KEYWORD_MISSING)
        status = HSTCAL_OK;

    /*check the PCTEDIM keyword in header*/
    tempStatus = GetKeyInt(group->globalhdr, "PCTETLEN", NO_DEFAULT, -999, &cte_len);
    if (tempStatus == HSTCAL_OK)
    {
        if (cte_len > 1)
            pars->cte_len = cte_len;
    }
    else if (tempStatus != KEYWORD_MISSING)
    {
        trlerror("(pctecorr) Error reading PCTETLEN keyword from header");
        return (status = tempStatus);
    }
    else if (tempStatus == KEYWORD_MISSING)
        status = HSTCAL_OK;

    /*check the PCTERNOI keyword in header*/
    tempStatus = GetKeyDbl(group->globalhdr, "PCTERNOI", NO_DEFAULT, -999, &rn_amp);
    if (tempStatus == HSTCAL_OK)
    {
        if (rn_amp > 1.)
            pars->rn_amp = rn_amp;
    }
    else if (tempStatus != KEYWORD_MISSING)
    {
        trlerror("(pctecorr) Error reading PCTERNOI keyword from header");
        return (status = tempStatus);
    }
    else if (tempStatus == KEYWORD_MISSING)
        status = HSTCAL_OK;

    /* get number of iterations used in forward model */
    tempStatus = GetKeyInt(group->globalhdr, "PCTENFOR", NO_DEFAULT, -999, &n_forward);
    if (tempStatus == HSTCAL_OK)
    {
        if (n_forward > 1)
            pars->n_forward = n_forward;
    }
    else if (tempStatus != KEYWORD_MISSING)
    {
        trlerror("(pctecorr) Error reading PCTENFOR keyword from header");
        return (status = tempStatus);
    }
    else if (tempStatus == KEYWORD_MISSING)
        status = HSTCAL_OK;

    /* get number of iterations used in parallel transfer*/
    tempStatus = GetKeyInt(group->globalhdr, "PCTENPAR", NO_DEFAULT, -999, &n_par);
    if (tempStatus == HSTCAL_OK)
    {
        if (n_par > 1)
            pars->n_par = n_par;
    }
    else if (tempStatus != KEYWORD_MISSING)
    {
        trlerror("(pctecorr) Error reading PCTENPAR keyword from header");
        return (status = tempStatus);
    }
    else if (tempStatus == KEYWORD_MISSING)
        status = HSTCAL_OK;

    /*fix the readout Cr's? */
    tempStatus = GetKeyInt(group->globalhdr, "FIXROCR", NO_DEFAULT, -999, &fix_rocr);
    if (tempStatus == HSTCAL_OK)
    {
        if (fix_rocr == 0 || fix_rocr == 1)
            pars->fix_rocr = fix_rocr;
        else
        {
            char msgBuffer[256];
            *msgBuffer = '\0';
            sprintf(msgBuffer, "(pctecorr) FIXROCR keyword from header has invalid value, '%d'. Only 0 or 1 accepted.", fix_rocr);
            trlerror(msgBuffer);
            return (status = INVALID_VALUE);
        }
    }
    else if (tempStatus != KEYWORD_MISSING)
    {
        trlerror("(pctecorr) Error reading FIXROCR keyword from header");
        return (status = tempStatus);
    }
    else if (tempStatus == KEYWORD_MISSING)
        status = HSTCAL_OK;

    trlmessage("End attempt to read CTE keywords from input image header.");

    return HSTCAL_OK;
}

void ctewarn (char *message) {

    char line[CHAR_LINE_LENGTH+1];

    // Use macro for beginning of Warning message
    sprintf(line,"%s",WARN_PREFIX);
    strcat (line,message);

    ctemessage(line);
}

void cteerror (char *message) {

    char line[CHAR_LINE_LENGTH+1];

    // Use macro for beginning of Warning message
    sprintf(line,"%s",ERR_PREFIX);
    strcat (line,message);

    ctemessage(line);
}

void ctemessage (char *message) {
    printf ("%s\n", message);
        fflush(stdout);
}

