# include <stdio.h>
# include <hstio.h>

# include "xtables.h"
# include "acs.h"
# include "acsinfo.h"
# include "acserr.h"

# define FULL_FRAME_READOUT  1
# define SUBARRAY_READOUT    2
# define BINNED_READOUT      3

# define NUMCOLS 18

typedef struct {
    IRAFPointer tp;            /* pointer to table descriptor */
    IRAFPointer cp_amp;        /* column descriptors */
    IRAFPointer cp_chip;
    IRAFPointer cp_binx;
    IRAFPointer cp_biny;
    IRAFPointer cp_nx;
    IRAFPointer cp_ny;
    IRAFPointer cp_trimx1;
    IRAFPointer cp_trimx2;
    IRAFPointer cp_trimy1;
    IRAFPointer cp_trimy2;
    IRAFPointer cp_vx1;
    IRAFPointer cp_vx2;
    IRAFPointer cp_vy1;
    IRAFPointer cp_vy2;
    IRAFPointer cp_biassecta1;
    IRAFPointer cp_biassecta2;
    IRAFPointer cp_biassectb1;
    IRAFPointer cp_biassectb2;
    int nrows;            /* number of rows in table */
} TblInfo;

typedef struct {
    char ccdamp[ACS_CBUF];
    int chip;
    int nx;
    int ny;
    int trimx[2];
    int trimy[2];
    int vx[2];
    int vy[2];
    int biassecta[2];
    int biassectb[2];
    int bin[2];
} TblRow;

static int OpenOverTab (char *, TblInfo *);
static int ReadOverTab (TblInfo *, int, TblRow *);
static int CloseOverTab (TblInfo *);


/* This routine assigns the widths of the overscan region on each edge
   of the image by assuming particular region sizes, depending on the
   binning and which amplifier was used for readout.

   This routine reads in the values from the OSCNTAB reference table.

   The OSCNTAB refers to all regions in image coordinates, even when
   the image is binned.


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
    AMPX      - First column affected by second AMP
    AMPY      - First line affected by second set of AMPS
    (VX1,VY1) - image coordinates of virtual overscan region origin
    (VX2,VY2) - image coordinates of top corner of virtual overscan region
    A1,A2     - beginning and ending columns for leading bias section
                (BIASSECTA1,BIASSECTA2 from OSCNTAB)
    B1,B2     - beginning and ending columns for trailing bias section
                (BIASSECTB1,BIASSECTB2 from OSCNTAB)
    TRIMX1    - Number of columns to trim off beginning of each line,
                contains A1,A2
    TRIMX2    - Number of columns to trim off end of each line,
                contains B1,B2
    TRIMY1    - Number of lines to trim off beginning of each column
    TRIMY2    - Number of line to trim off end of each column


   HISTORY:

   1998-06-02 Warren Hack  Initial version based on GetCCDTab.
   1998-10-20 Warren Hack  Second version. Condensed all input/output through
                           ACSInfo structure, and included both trailing and
                           leading overscan region designations.
   1999-02-19 Warren Hack  Revised to read CCDCHIP column to match CCDCHIP
                           keyword.
   2001-07-24 Warren Hack  Fixed a bug in calculating the second trim region
                           in a subarray.
   2001-09-24 Warren Hack  Removed another bug in determining trim regions for
                           subarrays, specifically, modified checks to be
                           relative to trimmed size of a full frame.
                           Also, check to see if overlap beginning trim
                           region first (and do both), then check for second
                           region overlap.
   2016-07-07 P. L. Lim    Removed no virtual overscan assumptions for new
                           subarrays added by FSW change in May 2016.
*/
int FindOverscan (ACSInfo *acs, int nx, int ny, int *overscan) {
    /* arguments:
       ACSInfo *acs          i: structure with all values from OSCNTAB
       int nx, ny            i: lengths of first and second image axes
       int overscan          o: indicates whether there is overscan or not
    */

    extern int status;
    int row;
    TblInfo tabinfo;
    TblRow tabrow;
    int foundit;
    int is_newsub;

    int cx0, cx1, tx1, tx2;
    int cy0, cy1, ty1, ty2;
    int full_nx;
    int full_ny;

    int SameInt (int, int);
    int SameString (char *, char *);

    /* Open the OSCNTAB parameters table and find columns. */
    if (OpenOverTab (acs->oscn.name, &tabinfo))
        return (status);

    /* Check each row for a match with ccdamp, chip,
       binx, and biny, and get info from the matching row. */

    foundit = NO;
    is_newsub = NO;
    *overscan = NO;

    for (row = 1; row <= tabinfo.nrows; row++) {

        /* Read the current row into tabrow. */
        if (ReadOverTab (&tabinfo, row, &tabrow))
            return (status);

        /* The overscan regions only depend on the full chip size. */
        if (SameString (tabrow.ccdamp, acs->ccdamp) &&
                SameInt (tabrow.chip, acs->chip) &&
                SameInt (tabrow.bin[0], acs->bin[0]) &&
                SameInt (tabrow.bin[1], acs->bin[1])) {

            foundit = YES;
            *overscan = YES;

            /* We are working with a subarray.
               There is never any virtual overscan EXCEPT for new subarrays
               added by FSW change in May 2016. */
            if (acs->subarray == YES) {
                if (SameString(acs->aperture, "WFC1A-512") ||
                        SameString(acs->aperture, "WFC1A-1K") ||
                        SameString(acs->aperture, "WFC1A-2K") ||
                        SameString(acs->aperture, "WFC1B-512") ||
                        SameString(acs->aperture, "WFC1B-1K") ||
                        SameString(acs->aperture, "WFC1B-2K") ||
                        SameString(acs->aperture, "WFC2C-512") ||
                        SameString(acs->aperture, "WFC2C-1K") ||
                        SameString(acs->aperture, "WFC2C-2K") ||
                        SameString(acs->aperture, "WFC2D-512") ||
                        SameString(acs->aperture, "WFC2D-1K") ||
                        SameString(acs->aperture, "WFC2D-2K")) {
                    is_newsub = YES;
                } else {
                    acs->trimy[0] = 0;
                    acs->trimy[1] = 0;
                }

                acs->vx[0] = 0;
                acs->vx[1] = 0;
                acs->vy[0] = 0;
                acs->vy[1] = 0;

                /* Virtual overscan processing is similar to physical
                   overscan processing below. */
                if (is_newsub == YES) {
                    /* Determine whether the subarray extends into the
                       virtual overscan regions on either side of the chip */
                    ty1 = (int)(ny - acs->offsety);
                    cy1 = tabrow.ny - tabrow.trimy[1] - tabrow.trimy[0];
                    cy0 = (int)(tabrow.trimy[0] - acs->offsety);

                    /* Subarray starts in the first overscan region... */
                    if (acs->offsety > 0) {
                        acs->trimy[0] = (acs->offsety < tabrow.trimy[0]) ?
                            acs->offsety : tabrow.trimy[0];
                        /* Check to see if it extends into second overscan
                           region. */
                        full_ny = tabrow.ny - (tabrow.trimy[0] +
                                               tabrow.trimy[1]);
                        ty2 = (int)(ny - acs->trimy[0]) - full_ny;
                        acs->trimy[1] = (ty2 < 0) ? 0 : ty2;

                    /* Subarray overlaps second overscan region */
                    } else if ( ty1 > cy1 && ty1 <= tabrow.ny) {
                        acs->trimy[0] = 0;
                        acs->trimy[1] = ty1 - cy1;

                    /* Subarray does not have virtual overscan */
                    } else {
                        acs->trimy[0] = 0;
                        acs->trimy[1] = 0;
                    }
                }

                /* Determine whether the subarray extends into the
                   physical overscan regions on either side of the chip */
                tx1 = (int)(nx - acs->offsetx);
                cx1 = tabrow.nx - tabrow.trimx[1] - tabrow.trimx[0];
                cx0 = (int)(tabrow.trimx[0] - acs->offsetx);

                /* Subarray starts in the first overscan region... */
                if (acs->offsetx > 0) {
                    acs->trimx[0] = (acs->offsetx < tabrow.trimx[0]) ?
                        acs->offsetx : tabrow.trimx[0];
                    /* Check to see if it extends into second overscan region.
                       Fixed 24 July 2001 WJH. */
                    full_nx = tabrow.nx - (tabrow.trimx[0] + tabrow.trimx[1]);
                    tx2 = (int)(nx - acs->trimx[0]) - full_nx;
                    acs->trimx[1] = (tx2 < 0) ? 0 : tx2;
                    acs->biassecta[0] = (tabrow.biassecta[0] - 1 - cx0 > 0) ?
                        (tabrow.biassecta[0] - 1 - cx0) : 0;
                    acs->biassecta[1] = tabrow.biassecta[1] - 1 - cx0;
                    acs->biassectb[0] = 0;
                    acs->biassectb[1] = 0;

                /* ... then the subarray overlaps biassectb... */
                } else if ( tx1 > cx1 && tx1 <= tabrow.nx) {
                    acs->trimx[0] = 0;
                    acs->trimx[1] = tx1 - cx1;
                    acs->biassecta[0] = 0;
                    acs->biassecta[1] = 0;
                    acs->biassectb[0] = nx - (tx1 - cx1);
                    acs->biassectb[1] = nx - 1;

                /* Subarray doesn't overlap either physical overscan region */
                } else {
                    acs->trimx[0] = 0;
                    acs->trimx[1] = 0;
                    acs->biassecta[0] = 0;
                    acs->biassecta[1] = 0;
                    acs->biassectb[0] = 0;
                    acs->biassectb[1] = 0;
                    *overscan = NO;
                }

            /* We are working with a full chip image... */
            } else {
                acs->trimx[0] = tabrow.trimx[0];
                acs->trimx[1] = tabrow.trimx[1];
                acs->trimy[0] = tabrow.trimy[0];
                acs->trimy[1] = tabrow.trimy[1];
                /* Subtract one from table values to
                   conform to C array indexing... */
                acs->vx[0] = tabrow.vx[0]-1;
                acs->vx[1] = tabrow.vx[1]-1;
                acs->vy[0] = tabrow.vy[0]-1;
                acs->vy[1] = tabrow.vy[1]-1;
                acs->biassecta[0] = tabrow.biassecta[0] - 1;
                acs->biassecta[1] = tabrow.biassecta[1] - 1;
                acs->biassectb[0] = tabrow.biassectb[0] - 1;
                acs->biassectb[1] = tabrow.biassectb[1] - 1;
            }

            break;
        }
    }

    if (foundit == NO) {
        status = ROW_NOT_FOUND;
        sprintf(MsgText, "Could not find appropriate row from OSCNTAB. ");
        trlwarn(MsgText);
    }

    if(CloseOverTab (&tabinfo))
        return(status);

    if (acs->verbose == YES) {
        sprintf(MsgText, "Found trim values of: x(%d,%d) y(%d,%d)",
                acs->trimx[0], acs->trimx[1], acs->trimy[0], acs->trimy[1]);
        trlmessage(MsgText);
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

    char *colnames[NUMCOLS] = {"CCDAMP", "CCDCHIP", "BINX", "BINY", "NX", "NY",
                               "TRIMX1", "TRIMX2", "TRIMY1", "TRIMY2",
                               "VX1", "VX2", "VY1", "VY2",
                               "BIASSECTA1", "BIASSECTA2",
                               "BIASSECTB1", "BIASSECTB2"};

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
    c_tbcfnd1 (tabinfo->tp, "TRIMY1", &tabinfo->cp_trimy1);
    c_tbcfnd1 (tabinfo->tp, "TRIMY2", &tabinfo->cp_trimy2);
    c_tbcfnd1 (tabinfo->tp, "VX1", &tabinfo->cp_vx1);
    c_tbcfnd1 (tabinfo->tp, "VX2", &tabinfo->cp_vx2);
    c_tbcfnd1 (tabinfo->tp, "VY1", &tabinfo->cp_vy1);
    c_tbcfnd1 (tabinfo->tp, "VY2", &tabinfo->cp_vy2);
    c_tbcfnd1 (tabinfo->tp, "BIASSECTA1", &tabinfo->cp_biassecta1);
    c_tbcfnd1 (tabinfo->tp, "BIASSECTA2", &tabinfo->cp_biassecta2);
    c_tbcfnd1 (tabinfo->tp, "BIASSECTB1", &tabinfo->cp_biassectb1);
    c_tbcfnd1 (tabinfo->tp, "BIASSECTB2", &tabinfo->cp_biassectb2);

    /* Initialize counters here... */
    missing = 0;
    i=0;
    /* Check which columns are missing */
    if (tabinfo->cp_amp == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_chip == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_binx == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_biny == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_nx == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_ny == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_vx1 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_vx2 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_vy1 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_vy2 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_trimx1 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_trimx2 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_trimy1 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_trimy2 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_biassecta1 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_biassecta2 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_biassectb1 == 0 ) { missing++; nocol[i] = YES; i++;}
    if (tabinfo->cp_biassectb2 == 0 ){ missing++; nocol[i] = YES; i++;}

    if (PrintMissingCols (missing, NUMCOLS, nocol, colnames, "OSCNTAB",
                          tabinfo->tp) )
        return(status);

    return (status);
}


/* This routine reads the relevant data from one row.
   It reads the amp configuration, chip number, binning factors,
   and x/y size of the image for selecting the row appropriate
   for the data. It also reads the trim regions, virtual overscan
   corner coordinates, and the bias section columns for both the
   trailing and leading overscan regions.
*/
static int ReadOverTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

    extern int status;

    c_tbegtt (tabinfo->tp, tabinfo->cp_amp, row,
            tabrow->ccdamp, ACS_CBUF-1);
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
    c_tbegti (tabinfo->tp, tabinfo->cp_vy1, row, &tabrow->vy[0]);
    if (c_iraferr())
        return (status = TABLE_ERROR);
    c_tbegti (tabinfo->tp, tabinfo->cp_vy2, row, &tabrow->vy[1]);
    if (c_iraferr())
        return (status = TABLE_ERROR);

    c_tbegti (tabinfo->tp, tabinfo->cp_biassecta1, row, &tabrow->biassecta[0]);
    if (c_iraferr())
        return (status = TABLE_ERROR);

    c_tbegti (tabinfo->tp, tabinfo->cp_biassecta2, row, &tabrow->biassecta[1]);
    if (c_iraferr())
        return (status = TABLE_ERROR);

    c_tbegti (tabinfo->tp, tabinfo->cp_biassectb1, row, &tabrow->biassectb[0]);
    if (c_iraferr())
        return (status = TABLE_ERROR);

    c_tbegti (tabinfo->tp, tabinfo->cp_biassectb2, row, &tabrow->biassectb[1]);
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
