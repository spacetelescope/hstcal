# include <stdio.h>
# include <string.h>  /* memcpy */
# include <stdlib.h>  /* malloc */

# include <xtables.h>
# include <hstio.h>

# include "../stis.h"
# include "../stisdef.h"
# include "../stiserr.h"
# include "../stistds.h"
# include "calstis1.h"

typedef struct {
        IRAFPointer tp;                 /* pointer to table descriptor */
        IRAFPointer cp_wl;              /* column descriptors */
        IRAFPointer cp_time;
        IRAFPointer cp_slope;
        IRAFPointer cp_nwl;
        IRAFPointer cp_nt;
        IRAFPointer cp_optelem;
        IRAFPointer cp_pedigree;
        IRAFPointer cp_descrip;
        IRAFPointer cp_reftemp;
        IRAFPointer cp_tempsens;
        int nrows;                      /* number of rows in table */
} TblInfo;

typedef struct {
        char opt_elem[STIS_CBUF];
} TblRow;

static int OpenTdsTab(char *, TblInfo *);
static int ReadTdsTab(TblInfo *, int, TblRow *);
static int ReadTdsArray(TblInfo *, int, TdsInfo *);
static int CloseTdsTab(TblInfo *);
static void FreeTds(TdsInfo *);

/*
  This routine gets the wavelength, time, and slopes from the
  time-dependent sensitivity (TDS) parameter table (keyword name
  TDSTAB).

  The TDS parameters table should contain the following data:
        header parameters:
                None
        columns:
                OPT_ELEM    optical element name
                WAVELENGTH  wavelength array
                TIME        time array
                SLOPE       slope at time and wavelength
                NWL         number of wavelengths in wavelength array
                NT          number of times in time array
                REFTEMP     reference temperature
                TEMPSENS    array (NWL) of temperature sensitivity factors
                PEDIGREE    row pedigree
                DESCRIP     row descriptor

  Rows are selected on OPT_ELEM.  If a matching row is found, the sizes
  of the wavelength and time arrays are read from NWL and NT, memory is
  allocated for the WAVELENGTH, TIME, and SLOPE arrays, and then they
  are read.  A linear interpolation is then done based on the exposure
  start date and the resulting TDS sensitivity values are copied into
  PhotInfo structure so they can be multiply by the throughput values.

  The wavelength, time, and slope arrays are double precision.

  Paul Barrett, 2003 Sep 25:
       Initial development, based on GetApThr1.

  Phil Hodge, 2004 Dec 27:
       Also get temperature info, if present.

  Phil Hodge, 2007 Mar 19:
       Don't include <c_iraf_priv.h>, and remove the declaration of c_tbciga.
*/

int GetTds1(StisInfo1 *sts, PhotInfo *phot) {
    /*
      StisInfo1 *sts   i: calibration switches and info
      PhotInfo  *phot  o: throughput values
    */

    void TdsCorrection (TdsInfo *, double, double, double *);

    TblInfo tabinfo;         /* pointer to table descriptor */
    TblRow tabrow;           /* value read from table row */

    int status = 0;
    int foundit = 0;
    int row;

    /* Open the TDS table. */
    if ((status = OpenTdsTab (sts->tdstab.name, &tabinfo)))
        return (status);

    for (row = 1; row <= tabinfo.nrows; row++) {
        if ((status = ReadTdsTab(&tabinfo, row, &tabrow)))
            return (status);

        if (SameString(tabrow.opt_elem, sts->opt_elem)) {
            TdsInfo tds;                /* TDS correction info struct */
            tds.allocated = 0;          /* memory allocation flag */
            foundit = 1;

            /* Get pedigree & descrip from the row. */
            if ((status = RowPedigree(&sts->tdstab, row, tabinfo.tp,
                                      tabinfo.cp_pedigree,
                                      tabinfo.cp_descrip)))
                return (status);
            if (sts->tdstab.goodPedigree == DUMMY_PEDIGREE) {
                printf("Warning  DUMMY pedigree in row %d of %s.\n",
                       row, sts->tdstab.name);
                sts->tdscorr = OMIT;
                return CloseTdsTab(&tabinfo);
            }
            /* Read time-dependent sensitivity info into TDS structure. */
            if ((status = ReadTdsArray(&tabinfo, row, &tds)))
                return (status);

            /* Set up PhotInfo structures.  */
            if (sts->filtcorr == PERFORM) {
                free(phot->f_wl);
                free(phot->f_thru);
            }
            phot->f_nelem = tds.nwl;
            phot->f_wl   = (double *) malloc(tds.nwl*sizeof(double));
            phot->f_thru = (double *) malloc(tds.nwl*sizeof(double));
            if (phot->f_wl == NULL || phot->f_thru == NULL)
                return (OUT_OF_MEMORY);

            /* Copy TDS row info into PhotInfo structure. */
            memcpy(phot->f_wl, tds.wl, tds.nwl*sizeof(double));
            TdsCorrection(&tds, sts->expstart, sts->detector_temp,
                          phot->f_thru);

            FreeTds(&tds);
        }
    }
    if ((status = CloseTdsTab(&tabinfo)))
        return (status);

    if (!foundit) {
        printf("Warning  Matching row not found in TDSTAB %s\n",
               sts->tdstab.name);
        printf("Warning  OPT_ELEM %s.\n", sts->opt_elem);
        printf("Warning  Skipping TDS correction.\n");
        sts->tdscorr = OMIT;
    }
    return (status);
}

/* This routine opens the time-dependent sensitivity table, finds the
   columns that we need, and gets the total number of rows in the table.
*/

static int OpenTdsTab(char *tname, TblInfo *tabinfo) {

    tabinfo->tp = c_tbtopn(tname, IRAF_READ_ONLY, 0);
    if (c_iraferr()) {
        printf("ERROR    TDSTAB `%s' not found\n", tname);
        return (OPEN_FAILED);
    }

    tabinfo->nrows = c_tbpsta(tabinfo->tp, TBL_NROWS);

    /* Find the columns. */

    c_tbcfnd1(tabinfo->tp, "OPT_ELEM",   &tabinfo->cp_optelem);
    c_tbcfnd1(tabinfo->tp, "WAVELENGTH", &tabinfo->cp_wl);
    c_tbcfnd1(tabinfo->tp, "TIME",       &tabinfo->cp_time);
    c_tbcfnd1(tabinfo->tp, "SLOPE",      &tabinfo->cp_slope);
    c_tbcfnd1(tabinfo->tp, "NWL",        &tabinfo->cp_nwl);
    c_tbcfnd1(tabinfo->tp, "NT",         &tabinfo->cp_nt);
    if (tabinfo->cp_wl    == 0 ||
        tabinfo->cp_time  == 0 ||
        tabinfo->cp_slope == 0 ||
        tabinfo->cp_nwl   == 0 ||
        tabinfo->cp_nt    == 0) {
        printf("ERROR    Column not found in TDSTAB\n");
        c_tbtclo(tabinfo->tp);
        return (COLUMN_NOT_FOUND);
    }
    /* Pedigree and descrip are optional columns. */
    c_tbcfnd1(tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
    c_tbcfnd1(tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

    /* Find the columns (if present) that are used for correcting
       for temperature dependence.
    */
    c_tbcfnd1 (tabinfo->tp, "REFTEMP",  &tabinfo->cp_reftemp);
    c_tbcfnd1 (tabinfo->tp, "TEMPSENS", &tabinfo->cp_tempsens);
    if (tabinfo->cp_reftemp == 0 ||tabinfo->cp_tempsens == 0) {
        printf ("Warning  Column REFTEMP or TEMPSENS not found in %s;\n",
                tname);
        printf ("Warning  no temperature correction applied to sensitivity\n");
    }

    return (0);
}

/* This routine reads the column name used to select the correct row. */

static int ReadTdsTab(TblInfo *tabinfo, int row, TblRow *tabrow) {

    c_tbegtt(tabinfo->tp, tabinfo->cp_optelem, row, tabrow->opt_elem,
             STIS_CBUF-1);
    if (c_iraferr())
        return (TABLE_ERROR);

    return (0);
}

/* This routine reads array data from one row. */

static int ReadTdsArray(TblInfo *tabinfo, int row, TdsInfo *tds) {

    int nwl, nt, ns, ntemp, dim[2], ini, ndim, i;

    /* Find out how many elements there are in the arrays. */

    c_tbegti(tabinfo->tp, tabinfo->cp_nwl, row, &tds->nwl);
    if (c_iraferr())
        return (TABLE_ERROR);
    c_tbegti(tabinfo->tp, tabinfo->cp_nt, row, &tds->nt);
    if (c_iraferr())
        return (TABLE_ERROR);

    /* Allocate memory. */

    tds->wl        = (double *) calloc(tds->nwl, sizeof(double));
    tds->temp_sens = (double *) calloc(tds->nwl, sizeof(double));
    tds->time      = (double *) calloc(tds->nt,  sizeof(double));
    tds->slope     = (double **) calloc(tds->nt, sizeof(double *));
    if (tds->temp_sens == NULL || tds->wl == NULL || tds->time == NULL ||
        tds->slope == NULL) {
        CloseTdsTab(tabinfo);
        return (OUT_OF_MEMORY);
    }
    for (i = 0; i < tds->nt; i++) {
        tds->slope[i] = (double *) calloc(tds->nwl, sizeof(double));
        if (tds->slope[i] == NULL) {
            c_tbtclo(tabinfo->tp);
            return (OUT_OF_MEMORY);
        }
    }
    tds->allocated = 1;

    /* Read 1-D arrays first. */

    nwl = c_tbagtd(tabinfo->tp, tabinfo->cp_wl, row, tds->wl, 1, tds->nwl);
    if (c_iraferr())
        return (TABLE_ERROR);
    nt = c_tbagtd(tabinfo->tp, tabinfo->cp_time, row, tds->time, 1, tds->nt);
    if (c_iraferr())
        return (TABLE_ERROR);

    /* The slope array is 2-dimensional, so we must get the stride first. */

    c_tbciga(tabinfo->tp, tabinfo->cp_slope, &ndim, dim, 2);
    ini = 1;    /* Arrays are 1-indexed in spp ! */
    for (i = 0; i < tds->nt; i++) {
        ns = c_tbagtd(tabinfo->tp, tabinfo->cp_slope, row, tds->slope[i],
                      ini, tds->nwl);
        if (c_iraferr())
            return (TABLE_ERROR);
        ini += dim[0];
    }
    /* Check if ingestion was properly satisfied. */

    if (nwl < tds->nwl || nt < tds->nt) {
        c_tbtclo(tabinfo->tp);
        printf("ERROR    Not all values were read from TDSTAB\n");
        return (TABLE_ERROR);
    }

    /* Get the reference temperature and temperature sensitivity,
       if they're present in this TDS table.
    */
    if (tabinfo->cp_reftemp == 0 || tabinfo->cp_tempsens == 0) {
        /* set the temperature sensitivity to zero */
        tds->ref_temp = -1.;
        for (i = 0; i < tds->nwl; i++)
            tds->temp_sens[i] = 0.;
    } else {
        c_tbegtd (tabinfo->tp, tabinfo->cp_reftemp, row, &tds->ref_temp);
        if (c_iraferr())
            return (TABLE_ERROR);
        ntemp = c_tbagtd (tabinfo->tp, tabinfo->cp_tempsens, row,
                    tds->temp_sens, 1, tds->nwl);
        if (c_iraferr())
            return (TABLE_ERROR);
        if (ntemp < tds->nwl) {
            c_tbtclo (tabinfo->tp);
            printf (
            "ERROR    Not all TEMPSENS values were read from TDSTAB\n");
            return (TABLE_ERROR);
        }
    }

    return (0);
}

/* This routine closes the TDSTAB table. */

static int CloseTdsTab(TblInfo *tabinfo) {

    c_tbtclo(tabinfo->tp);
    if (c_iraferr())
        return (TABLE_ERROR);

    return (0);
}

/* And this routine frees memory. */

void FreeTds(TdsInfo *tds) {

    int i;

    if (tds->allocated) {
        free(tds->wl);
        free(tds->temp_sens);
        free(tds->time);
        for (i = 0; i < tds->nt; i++)
            free(tds->slope[i]);
        free(tds->slope);
        tds->allocated = 0;
    }
}
