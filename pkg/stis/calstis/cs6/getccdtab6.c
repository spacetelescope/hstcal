# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"
# include "stis.h"
# include "calstis6.h"
# include "err.h"
# include "stisdef.h"

typedef struct {
    IRAFPointer tp;			/* pointer to table descriptor */
    IRAFPointer cp_amp;		        /* column descriptors */
    IRAFPointer cp_ccdgain;
    IRAFPointer cp_ccdoffset;
    IRAFPointer cp_bin1;
    IRAFPointer cp_bin2;
    IRAFPointer cp_atodgain;
    IRAFPointer cp_ctinorm;
    IRAFPointer cp_ctigpower;
    IRAFPointer cp_ctibgfac;
    IRAFPointer cp_ctibgpower;
    IRAFPointer cp_ctitimfac;
    IRAFPointer cp_ctirefmjd;
    IRAFPointer cp_spurcharge;
    IRAFPointer cp_halofac;
    IRAFPointer cp_halominfrac;
    IRAFPointer cp_pedigree;
    IRAFPointer cp_descrip;
    int nrows;			        /* number of rows in table */
} TblInfo;

typedef struct {
    char ccdamp[STIS_CBUF+1];
    int ccdgain;
    int bin[2];
    int ccdoffset;
    float atodgain;
    float ctinorm;
    float ctigpower;
    float ctibgfac;
    float ctibgpower;
    float ctitimfac;
    float ctirefmjd;
    float spurcharge;
    float halofac;
    float halominfrac;
} TblRow;

static int OpenCCDTab6 (char *, TblInfo *);
static int ReadCCDTab6 (TblInfo *, int, TblRow *);
static int CloseCCDTab6 (TblInfo *);

/* This routine gets the CCD charge transfer inefficiency (CTI)
   parameters from the CCD parameters table (keyword name CCDTAB).

   The CCD parameters table should contain the following:
   header parameters:
       none needed
   columns:
       CCDAMP:      identifies which amp was used (string A-D)
       CCDGAIN:     commanded gain of the CCD (int)
       CCDOFFST:    commanded bias (int)
       BINAXIS1,    BINAXIS2:  commanded bin sizes (int)
       ATODGAIN:    actual gain of the CCD (double)
       CCDBIAS:     typical value of bias (double)
       READNSE:     typical value of readout noise (double)
       SATURATE:    CCD saturation threshold

       CTINORM:     Normalization for zero background (double)
       CTIGPOWER:   Power value in CTI dependence on gross counts (double)
       CTIBGFAC:    Scaling factor in exponential roll-off (double)
       CTIBGPOWER:  Power value in exponential roll-off (double)
       CTIREFMJD:   Reference MJD for CTI correction (double)
       CTITIMFAC:   Scaling factor in CTI time dependence (double)
       SPURCHARGE:  Spurious Charge per pixel (double)
       HALOFAC:     Factor multiplying halopar (double)
       HALOMINFRAC: Minimum allowed value of frac-halo (double)

   The table is read to find the row for which the CCDAMP matches the
   expected amplifier name (a single letter, from the image header)
   and for which the commanded CCD gain, commanded CCD bias, and bin
   sizes match the values read from the image header.  Then that row
   is read to get the actual CTI norm, gpow, bgfac, bgpow, mjd, tmfac,
   and spurious charge.

   Paul Barrett, 2003 Jun 16:
      Copied from cs1/getccdtab.c and modified.

   Phil Hodge, 2006 Feb 20
      Also get columns HALOFAC and HALOMINFRAC.
*/

int GetCCDTab6(StisInfo6 *sts, CTICorrInfo *cti) {

/* arguments:
StisInfo1 *sts     io: calibration switches, etc
CTICorrInfo *cti   io: CTI correction constants
*/

    int status;

    TblInfo tabinfo;	        /* pointer to table descriptor, etc */
    TblRow tabrow;		/* values read from a table row */

    int row;		        /* loop index */
    int foundit;		/* has the correct row been found? */

    /* Open the CCD parameters table and find columns. */
    if ((status = OpenCCDTab6(sts->ccdtab.name, &tabinfo)))
        return (status);

    /* Check each row for a match with ccdamp, ccdgain, ccdoffst,
       binaxis1, and binaxis2, and get info from the matching row.
    */

    foundit = 0;
    for (row = 1; row <= tabinfo.nrows; row++) {

        /* Read the current row into tabrow. */
        if ((status = ReadCCDTab6(&tabinfo, row, &tabrow)))
            return (status);

        /* Check for matching row.  If found, copy data to relavant
           data structures.
        */
        if (SameString(tabrow.ccdamp, sts->ccdamp) &&
            SameInt(tabrow.ccdgain, sts->ccdgain) &&
            SameInt(tabrow.ccdoffset, sts->ccdoffset) &&
            SameInt(tabrow.bin[0], sts->binaxis[0]) &&
            SameInt(tabrow.bin[1], sts->binaxis[1])) {

            foundit = 1;
            if ((status = RowPedigree(&sts->ccdtab, row, tabinfo.tp,
                                      tabinfo.cp_pedigree, tabinfo.cp_descrip)))
                return (status);
            if (sts->ccdtab.goodPedigree == DUMMY_PEDIGREE)
                printf ("Warning  Row %d of CCDTAB is DUMMY.\n", row);
            sts->atodgain    = tabrow.atodgain;
            cti->ctinorm     = tabrow.ctinorm;
            cti->ctigpower   = tabrow.ctigpower;
            cti->ctibgfac    = tabrow.ctibgfac;
            cti->ctibgpower  = tabrow.ctibgpower;
            cti->ctirefmjd   = tabrow.ctirefmjd;
            cti->ctitimfac   = tabrow.ctitimfac;
            cti->spurcharge  = tabrow.spurcharge;
            cti->halofac     = tabrow.halofac;
            cti->halominfrac = tabrow.halominfrac;
            break;
        }
    }

    if (!foundit) {
        printf("ERROR    Matching row not found in CCDTAB `%s'.\n",
               sts->ccdtab.name);
        printf("ERROR    CCDAMP %s, CCDGAIN %d, CCDOFFST %d, BINAXIS %d,%d.\n",
               sts->ccdamp, sts->ccdgain, sts->ccdoffset,
               sts->binaxis[0], sts->binaxis[1]);
        CloseCCDTab6(&tabinfo);
        return (TABLE_ERROR);
    }

    if ((status = CloseCCDTab6(&tabinfo)))
        return (status);

    return (0);
}

/* This routine opens the CCD parameters table, finds the columns that
   we need, and gets the total number of rows in the table.  The
   columns are CCDAMP, CCDGAIN, CCDOFFSET, BINAXIS1, BINAXIS2,
   ATODGAIN, CTINORM, CTIGPOWER, CTIBGFAC, CTIBGPOWER, CTITIMFAC,
   CTIREFMJD, and SPURCHARGE.
*/

static int OpenCCDTab6(char *tname, TblInfo *tabinfo)
{
    tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
    if (c_iraferr()) {
        printf ("ERROR    CCDTAB `%s' not found.\n", tname);
        return (OPEN_FAILED);
    }

    tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

    /* Find the columns. */
    c_tbcfnd1(tabinfo->tp, "CCDAMP",     &tabinfo->cp_amp);
    c_tbcfnd1(tabinfo->tp, "CCDGAIN",    &tabinfo->cp_ccdgain);
    c_tbcfnd1(tabinfo->tp, "CCDOFFST",   &tabinfo->cp_ccdoffset);
    c_tbcfnd1(tabinfo->tp, "BINAXIS1",   &tabinfo->cp_bin1);
    c_tbcfnd1(tabinfo->tp, "BINAXIS2",   &tabinfo->cp_bin2);
    c_tbcfnd1(tabinfo->tp, "ATODGAIN",   &tabinfo->cp_atodgain);
    c_tbcfnd1(tabinfo->tp, "CTINORM",    &tabinfo->cp_ctinorm);
    c_tbcfnd1(tabinfo->tp, "CTIGPOWER",  &tabinfo->cp_ctigpower);
    c_tbcfnd1(tabinfo->tp, "CTIBGFAC",   &tabinfo->cp_ctibgfac);
    c_tbcfnd1(tabinfo->tp, "CTIBGPOWER", &tabinfo->cp_ctibgpower);
    c_tbcfnd1(tabinfo->tp, "CTITIMFAC",  &tabinfo->cp_ctitimfac);
    c_tbcfnd1(tabinfo->tp, "CTIREFMJD",  &tabinfo->cp_ctirefmjd);
    c_tbcfnd1(tabinfo->tp, "SPURCHARGE", &tabinfo->cp_spurcharge);
    if (tabinfo->cp_amp == 0 ||
        tabinfo->cp_ccdgain == 0 ||
        tabinfo->cp_bin1 == 0 ||
        tabinfo->cp_bin2 == 0 ||
        tabinfo->cp_atodgain == 0 ||
        tabinfo->cp_ctinorm == 0 ||
        tabinfo->cp_ctigpower == 0 ||
        tabinfo->cp_ctibgfac == 0 ||
        tabinfo->cp_ctibgpower == 0 ||
        tabinfo->cp_ctirefmjd == 0 ||
        tabinfo->cp_ctitimfac == 0 ||
        tabinfo->cp_spurcharge == 0) {
        c_tbtclo (tabinfo->tp);
        return (COLUMN_NOT_FOUND);
    }

    /* These are optional, only needed for G750L and G750M. */
    c_tbcfnd1(tabinfo->tp, "HALOFAC", &tabinfo->cp_halofac);
    c_tbcfnd1(tabinfo->tp, "HALOMINFRAC", &tabinfo->cp_halominfrac);

    /* Pedigree and descrip are optional columns. */
    c_tbcfnd1(tabinfo->tp, "PEDIGREE",   &tabinfo->cp_pedigree);
    c_tbcfnd1(tabinfo->tp, "DESCRIP",    &tabinfo->cp_descrip);

    return (0);
}

/* This routine reads the relevant data from one row.  The amplifier
   number, CCD gain, and CTI values are gotten.
*/

static int ReadCCDTab6(TblInfo *tabinfo, int row, TblRow *tabrow)
{
    c_tbegtt(tabinfo->tp, tabinfo->cp_amp, row, tabrow->ccdamp, STIS_CBUF);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegti(tabinfo->tp, tabinfo->cp_ccdgain, row, &tabrow->ccdgain);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegti(tabinfo->tp, tabinfo->cp_ccdoffset, row, &tabrow->ccdoffset);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegti(tabinfo->tp, tabinfo->cp_bin1, row, &tabrow->bin[0]);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegti(tabinfo->tp, tabinfo->cp_bin2, row, &tabrow->bin[1]);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegtr(tabinfo->tp, tabinfo->cp_atodgain, row, &tabrow->atodgain);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegtr(tabinfo->tp, tabinfo->cp_ctinorm, row, &tabrow->ctinorm);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegtr(tabinfo->tp, tabinfo->cp_ctigpower, row, &tabrow->ctigpower);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegtr(tabinfo->tp, tabinfo->cp_ctibgfac, row, &tabrow->ctibgfac);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegtr(tabinfo->tp, tabinfo->cp_ctibgpower, row, &tabrow->ctibgpower);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegtr(tabinfo->tp, tabinfo->cp_ctitimfac, row, &tabrow->ctitimfac);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegtr(tabinfo->tp, tabinfo->cp_ctirefmjd, row, &tabrow->ctirefmjd);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegtr(tabinfo->tp, tabinfo->cp_spurcharge, row, &tabrow->spurcharge);
    if (c_iraferr())
        return (TABLE_ERROR);

    if (tabinfo->cp_halofac > 0) {
        c_tbegtr(tabinfo->tp, tabinfo->cp_halofac, row, &tabrow->halofac);
        if (c_iraferr())
            return (TABLE_ERROR);
    } else {
        tabrow->halofac = 0.;
    }

    if (tabinfo->cp_halominfrac > 0) {
        c_tbegtr(tabinfo->tp, tabinfo->cp_halominfrac, row,
			&tabrow->halominfrac);
        if (c_iraferr())
            return (TABLE_ERROR);
    } else {
        tabrow->halominfrac = 0.;
    }

    return (0);
}

/* This routine closes the ccdtab table. */

static int CloseCCDTab6(TblInfo *tabinfo)
{
    c_tbtclo (tabinfo->tp);
    if (c_iraferr())
        return (TABLE_ERROR);

    return (0);
}
