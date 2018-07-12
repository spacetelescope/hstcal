# include <stdio.h>
# include <string.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "hstcalerr.h"
# include "stisdq.h"

# define SZ_KWD 80      /* size of keyword string buffer */

static int InheritHeader (StisInfo6 *, TblDesc *, int);
static int isInvalid (char *);


/*
   Routines that handle the output tables. There is support for both
   the standard calstis6 output table, and the profile generator table.




   Revision history:
   ----------------
   20 Feb 97  -  Implemented (I.Busko)
   10 Apr 97  -  Changes after code review (IB):
                 - replaced literals by SZ_KWD constant.
                 - renamed output column "Angstroms" to comply to ISR spec.
   22 Apr 97  -  Changed column names: NELEM, WAVELENGTH, BACKGROUND (IB)
   22 Apr 97  -  Do not write empty spectrum (IB)
   08 May 97  -  Conform to new _trl standard (IB)
   29 May 97  -  Changed arguments in WriteRow (IB)
   17 Jul 97  -  Background is in counts/sec (IB)
   17 Jul 97  -  Add HistoryAddHeader routine (IB)
   13 Apr 98  -  Add extraction info columns to output table (IB)
   23 Jun 98  -  Add column with croscorr offset in each sporder (IB)
   24 Jun 98  -  Profile generator (IB)
   02 Oct 98  -  Write keywords from all input IMSETs into output (IB)
   28 Oct 98  -  Add strUpdateHeader function (IB)
   17 Dec 98  -  Changed column name from MAXSEARCH to MAXSRCH (IB)
   23 Feb 00  -  Add doubleUpdateHeader function (IB)
   11 Apr 00  -  Interpolated profile is one pixel smaller (IB)
   16 Jun 00  -  Write uncompressed profile arrays (IB)
   01 Nov 00  -  Profile offset (IB)
   01 Dec 00  -  Output subsampled profile (IB)
   30 Jan 01  -  Centroid output (IB)
   09 May 01  -  Remove doubleUpdateHeader function (IB)
   17 Apr 02  -  Add doubleUpdateHeader function (IB)
   13 Nov 02  -  Remove invalid kwds from output table header (OPR #46956) (IB)
   10 Apr 03  -  Remove additional invalid kwds from header (OPR #48303) (PB)
   10 Apr 03  -  Rewrote isInvalid function to use a string array. (PB)
    3 Nov 08  -  Change the description for extnum in CreaTable; extnum is
                 the image set in the input file, not the extension in the
                 output table (PEH)
   21 Oct 11  -  New column NET_ERROR was added (PEH)
*/



/* Open the cs6 output table extension, create column descriptors and
   physically create the table. The table extension must be opened
   with the raw output file name, so the table routines can properly
   append it as a FITS extension into the existing file.
*/

int CreaTable (StisInfo6 *sts, int extnum, TblDesc *table) {

/* arguments:
StisInfo6 *sts      i: calibration switches and info
int       extnum    i: current image set in input file
TblDesc   *table    o: table descriptor
*/

        int status;

        if ((table->tp = c_tbtopn (sts->output, IRAF_NEW_FILE,0)) == 0)
            return (OPEN_FAILED);
        c_tbcdef1 (table->tp, &(table->sporder), SPORDER, "","",
                   IRAF_SHORT, 1);
        c_tbcdef1 (table->tp, &(table->npts),    NELEM, "","",
                   IRAF_SHORT, 1);
        c_tbcdef1 (table->tp, &(table->wave),    WAVELENGTH, "Angstroms","",
                   IRAF_DOUBLE, table->array_size);
        c_tbcdef1 (table->tp, &(table->gross),   GROSS, "Counts/s","",
                   IRAF_REAL, table->array_size);
        c_tbcdef1 (table->tp, &(table->back),    BACKGROUND,"Counts/s",
                   "", IRAF_REAL, table->array_size);
        c_tbcdef1 (table->tp, &(table->net),     NET, "Counts/s","",
                   IRAF_REAL, table->array_size);
        c_tbcdef1 (table->tp, &(table->flux),    "FLUX",
                   "erg/s/cm**2/Angstrom","",IRAF_REAL, table->array_size);

        if (sts->fluxcorr == PERFORM)
            c_tbcdef1 (table->tp, &(table->error), "ERROR",
                      "erg/s/cm**2/Angstrom","",IRAF_REAL, table->array_size);
        else
            c_tbcdef1 (table->tp, &(table->error), "ERROR",
                      "Counts/s","",IRAF_REAL, table->array_size);
	c_tbcdef1 (table->tp, &(table->net_error), "NET_ERROR",
		"Counts/s", "", IRAF_REAL, table->array_size);

        c_tbcdef1 (table->tp, &(table->dq),      "DQ", "","",
                   IRAF_SHORT, table->array_size);

        if (sts->extrloc) {
            c_tbcdef1 (table->tp, &(table->a2center),  "A2CENTER",
                       "pixel", "", IRAF_REAL, 1);
            c_tbcdef1 (table->tp, &(table->extrsize),  "EXTRSIZE",
                       "pixel", "", IRAF_REAL, 1);
            c_tbcdef1 (table->tp, &(table->maxsearch), "MAXSRCH",
                       "pixel", "", IRAF_SHORT, 1);
            c_tbcdef1 (table->tp, &(table->bk1size),   "BK1SIZE",
                       "pixel", "", IRAF_REAL, 1);
            c_tbcdef1 (table->tp, &(table->bk2size),   "BK2SIZE",
                       "pixel", "", IRAF_REAL, 1);
            c_tbcdef1 (table->tp, &(table->bk1offset), "BK1OFFST",
                       "pixel", "", IRAF_REAL, 1);
            c_tbcdef1 (table->tp, &(table->bk2offset), "BK2OFFST",
                       "pixel", "", IRAF_REAL, 1);
            c_tbcdef1 (table->tp, &(table->extrlocy),  EXTRLOCY,
                      "pixel", "", IRAF_REAL, table->array_size);
        }

        c_tbcdef1 (table->tp, &(table->cc_offset),  "OFFSET",
                   "pixel", "", IRAF_REAL, 1);
        c_tbtcre (table->tp);

        /* The output table extension header inherits all keywords found
           in the [SCI] extension of the current IMSET in the input image.
        */
        if ( (status = InheritHeader (sts, table, extnum)) )
            return (status);

        return (STIS_OK);

}




/* Open the profile output table extension, create column descriptors
   and physically create the table. The table extension must be opened
   with the raw output file name, so the table routines can properly
   append it as a FITS extension into the existing file.
*/

int CreaProfTable (StisInfo6 *sts, ProfTblDesc *table) {

/* arguments:
StisInfo6 *sts             i: calibration switches and info
ProfTblDesc   *table    o: table descriptor
*/

        if ((table->tp = c_tbtopn (sts->output, IRAF_NEW_FILE,0)) == 0)
            return (OPEN_FAILED);

        c_tbcdef1 (table->tp, &(table->sporder), "SPORDER", "","",
                   IRAF_SHORT, 1);
        c_tbcdef1 (table->tp, &(table->npts), "NPTS", "","",
                   IRAF_SHORT, 1);
        c_tbcdef1 (table->tp, &(table->nptsoff), "NPTS_OFFSET", "","",
                   IRAF_SHORT, 1);
        c_tbcdef1 (table->tp, &(table->subscale), "SUB_FACTOR", "","",
                   IRAF_DOUBLE, 1);
        c_tbcdef1 (table->tp, &(table->minwave), "MIN_WAVE", "","",
                   IRAF_DOUBLE, 1);
        c_tbcdef1 (table->tp, &(table->maxwave), "MAX_WAVE", "","",
                   IRAF_DOUBLE, 1);
        c_tbcdef1 (table->tp, &(table->minpix), "MIN_PIX", "","",
                   IRAF_SHORT, 1);
        c_tbcdef1 (table->tp, &(table->maxpix), "MAX_PIX", "","",
                   IRAF_SHORT, 1);
        c_tbcdef1 (table->tp, &(table->s_n), "S_N", "","",
                   IRAF_REAL, 1);
        c_tbcdef1 (table->tp, &(table->profoff), "PROF_OFFSET", "", "",
                   IRAF_REAL, table->array_size_off);
        c_tbcdef1 (table->tp, &(table->profcent), "PROF_CENTROID", "", "",
                   IRAF_REAL, table->array_size_off);
        c_tbcdef1 (table->tp, &(table->prof), "PROF", "", "",
                   IRAF_REAL, table->array_size);

        c_tbtcre (table->tp);

        return (STIS_OK);
}



/* Write an output row with the results of 1-D extraction. The row number
   is incremented before writing the data. If an empty spectrum is found,
   no output takes place and the row number is not incremented (NOT ACTIVE).
*/

int WriteRow (StisInfo6 *sts, TblDesc *table, RowContents *row,
              int *row_number) {

/* arguments:
StisInfo6 *sts          i:  calibration switches and info
TblDesc table           i:  table descriptor
RowContents row         i:  row contents
int *row_number         io: row where to write in output table
*/

        int i, count;

        /* If spectrum is empty, do not write it. Note that the specific
           values used to detect an empty pixel must be the same used by
           routine Interp2D.
        */
        count = 0;
        for (i = 0; i < row->npts; i++) {
            if (row->gross[i] != 0.0F        ||
                row->error[i] != 1.0F        ||
                row->dq[i]    != DETECTORPROB)
                count++;
        }
        if (count == 0) {
/*          printf (
            "Warning  Empty spectrum, do not write row %d\n",(*row_number)+1);
            return (1);
        }
*/
            printf ("Warning  Empty spectrum.\n");
        }

        /* Increment row number. */
        (*row_number)++;

        c_tbapts (table->tp, table->sporder, *row_number, &(row->sporder),1,1);
        c_tbapts (table->tp, table->npts, *row_number, &(row->npts),   1,1);
        c_tbaptd (table->tp, table->wave,  *row_number, row->wave,  1,
                  (int)row->npts);
        c_tbaptr (table->tp, table->gross, *row_number, row->gross, 1,
                  (int)row->npts);
        c_tbaptr (table->tp, table->back,  *row_number, row->back,  1,
                  (int)row->npts);
        c_tbaptr (table->tp, table->net,   *row_number, row->net,   1,
                  (int)row->npts);
        c_tbaptr (table->tp, table->flux,  *row_number, row->flux,  1,
                  (int)row->npts);
        c_tbaptr (table->tp, table->error, *row_number, row->error, 1,
                  (int)row->npts);
	c_tbaptr (table->tp, table->net_error, *row_number, row->net_error, 1,
                  (int)row->npts);
        c_tbapts (table->tp, table->dq, *row_number, row->dq,       1,
                  (int)row->npts);
        if (sts->extrloc) {
            c_tbaptr (table->tp, table->a2center, *row_number,
                      &(row->a2center), 1, 1);
            c_tbaptr (table->tp, table->extrsize, *row_number,
                      &(row->extrsize), 1, 1);
            c_tbaptr (table->tp, table->bk1size, *row_number,
                      &(row->bk1size), 1, 1);
            c_tbaptr (table->tp, table->bk2size, *row_number,
                      &(row->bk2size), 1, 1);
            c_tbaptr (table->tp, table->bk1offset, *row_number,
                      &(row->bk1offset), 1, 1);
            c_tbaptr (table->tp, table->bk2offset, *row_number,
                      &(row->bk2offset), 1, 1);
            c_tbapts (table->tp, table->maxsearch, *row_number,
                      &(row->maxsearch), 1, 1);
            c_tbaptr (table->tp, table->extrlocy, *row_number,
                      row->extrlocy, 1, (int)row->npts);
        }
        c_tbaptr (table->tp, table->cc_offset, *row_number,
                  &(row->cc_offset), 1, 1);

        if (sts->verbose == 1 || sts->verbose == 2)
            printf ("         Row %d written to disk.\n", *row_number);

        return (0);
}




/* Write rows into output profile file. One or more rows are written
   depending on how the spectrum is broken up in pieces by the
   wavelength/pixel range selector.
*/

int WriteProfRow (StisInfo6 *sts, ProfTblDesc *table, RowContents *row,
                  int *row_number) {

/* arguments:
StisInfo6 *sts          i:  calibration switches and info
ProfTblDesc table       i:  table descriptor
RowContents row         i:  row contents structure
int *row_number         io: row where to start writing in output table.
                            This parameter is updated so when this function
                             returns it points to the last written row.
*/
        int i, ii, j, asize, k, k_off;
        short shold, asize_off;
        float *profile, *profile_off, *profile_cent;

        void BuildOutputProfile (StisInfo6 *, RowContents *);
        void FreeOutputProfile (StisInfo6 *);

        /* Build compressed profile. */
        BuildOutputProfile (sts, row);

        /* Create local buffers. These are used to convert double arrays
           to float. The conversion is necessary since internally these
           arrays have to be kept as doubles since they are reused by the
           1-D extraction code. But in the file they have to be kept as
           floats to decrease file size. Notice that the array size is
           taken from the first range bin. This presumes that bin sizes
           are constant except eventually for the last one. Notice also
           that bin ranges are only available after BuildOutputProfile
           executes.
        */
        asize = sts->subprof_size;
        asize_off = sts->profile_maxp[0] - sts->profile_minp[0] + 1;

        profile      = (float *) malloc (asize     * sizeof (float));
        profile_off  = (float *) malloc (asize_off * sizeof (float));
        profile_cent = (float *) malloc (asize_off * sizeof (float));
        if (profile == NULL || profile_off == NULL || profile_cent == NULL)
            return (OUT_OF_MEMORY);

        for (i = 0; i < sts->profile_msize; i++) {

            /* Increment row number. */
            (*row_number)++;

            /* Write a single row. */

            c_tbapts (table->tp, table->sporder, *row_number,
                      &(row->sporder),1,1);
            shold = sts->subprof_size;
            c_tbapts (table->tp, table->npts, *row_number,
                      &shold, 1,1);
            c_tbapts (table->tp, table->nptsoff, *row_number,
                      &asize_off, 1,1);
            c_tbaptd (table->tp, table->subscale, *row_number,
                      &(sts->subscale), 1,1);
            c_tbaptd (table->tp, table->minwave, *row_number,
                      &(sts->profile_minw[i]), 1,1);
            c_tbaptd (table->tp, table->maxwave, *row_number,
                      &(sts->profile_maxw[i]), 1,1);
            c_tbapts (table->tp, table->minpix, *row_number,
                      &(sts->profile_minp[i]), 1,1);
            c_tbapts (table->tp, table->maxpix, *row_number,
                      &(sts->profile_maxp[i]), 1,1);
            c_tbaptd (table->tp, table->s_n, *row_number,
                      &(sts->profile_sn[i]), 1,1);

            /* Move profile array data to output float arrays so they
               can be written with a single call to the table routines.
            */
            k     = 0;
            k_off = 0;
            for (ii =  sts->profile_minp[i] - 1;
                 ii <= sts->profile_maxp[i] - 1; ii++) {
                profile_off[k_off]  = (float)sts->profile_offset[ii];
                profile_cent[k_off] = (float)sts->profile_centroid[ii];
                k_off++;
                for (j =  0; j < sts->subprof_size; j++)
                    profile[j]  = (float)sts->subprof[i][j];
            }

            c_tbaptr (table->tp, table->prof, *row_number, profile, 1, asize);
            c_tbaptr (table->tp, table->profoff, *row_number, profile_off, 1,
                      asize_off);
            c_tbaptr (table->tp, table->profcent, *row_number, profile_cent, 1,
                      asize_off);

            printf ("         Row %d written to disk.\n", *row_number);
        }

        free (profile);
        free (profile_off);
        free (profile_cent);

        /* Free compressed profile. */
        FreeOutputProfile (sts);

        return (0);
}


/*  Utility functions to update (or write) header keywords. */

int intUpdateHeader (char *name, char *keywname, int keyval, char *comment,
                     int extension) {

        int status;
        Hdr hdr;                /* header structure */
        IODescPtr tb;           /* image or table descriptor */

        initHdr (&hdr);
        tb = openUpdateImage (name, "", extension, &hdr);
        if (hstio_err())
            return (OPEN_FAILED);
        if ( (status = Put_KeyI (&hdr, keywname, keyval, comment)) )
            return (status);
        putHeader (tb);
        if (hstio_err())
            return (OPEN_FAILED);
        closeImage (tb);
        freeHdr (&hdr);
        return (0);
}

int doubleUpdateHeader (char *name, char *keywname, double keyval, char *comment,
                        int extension) {

        int status;
        Hdr hdr;                /* header structure */
        IODescPtr tb;           /* image or table descriptor */

        initHdr (&hdr);
        tb = openUpdateImage (name, "", extension, &hdr);
        if (hstio_err())
            return (OPEN_FAILED);
        if ( (status = Put_KeyD (&hdr, keywname, keyval, comment)) )
            return (status);
        putHeader (tb);
        if (hstio_err())
            return (OPEN_FAILED);
        closeImage (tb);
        freeHdr (&hdr);
        return (0);
}

int strUpdateHeader (char *name, char *keywname, char *keyval, char *comment,
                     int extension) {

        int status;

        Hdr hdr;                /* header structure */
        IODescPtr tb;           /* image or table descriptor */

        initHdr (&hdr);
        tb = openUpdateImage (name, "", extension, &hdr);
        if (hstio_err())
            return (OPEN_FAILED);
        if ( (status = Put_KeyS (&hdr, keywname, keyval, comment)) )
            return (status);
        putHeader (tb);
        if (hstio_err())
            return (OPEN_FAILED);
        closeImage (tb);
        freeHdr (&hdr);
        return (0);
}


/*  Write HISTORY keywords in primary header. This is used to
    add HISTORY keywords in the primary header of the output table,
    which is handled as the primary header of an image.
*/

int HistoryAddHeader (char *name, char *text) {

        int status;

        Hdr hdr;                /* header structure */
        IODescPtr tb;           /* image or table descriptor */

        int addHistoryKw (Hdr *, char *);

        initHdr (&hdr);
        tb = openUpdateImage (name, "", 0, &hdr);
        if (hstio_err())
            return (OPEN_FAILED);
        if ( (status = addHistoryKw (&hdr, text)) )
            return (status);
        putHeader (tb);
        if (hstio_err())
            return (OPEN_FAILED);
        closeImage (tb);
        freeHdr (&hdr);
        return (0);
}




/*  Copies all keywords found in the SCI extension of the curent
    IMSET of the input image into the output table extension header.

    This function filters out keywords that are not allowed in a table.
    At this point, only the keywords explicitly listed in OPR #46956 are
    included in the filter. We still need a throughout test pass with
    fverify to generate a complete keyword list.
*/
static int InheritHeader (StisInfo6 *sts, TblDesc *table, int extnum) {

        IODescPtr im;           /* input extension header descriptor */
        Hdr hdr;                /* header structure */
        FitsKw kw;              /* keyword structure */
        char *kw_name;          /* keyword name */
        char *kw_comm;          /* keyword comment */
        FitsDataType kw_type;   /* keyword type */
        Bool kw_bool;           /* boolean keyword value */
        int kw_int;             /* integer keyword value */
        double kw_double;       /* double keyword value */
        char kw_string[SZ_KWD]; /* string keyword value */

        im = openInputImage (sts->input, "SCI", extnum);
        if (hstio_err())
            return (OPEN_FAILED);
        initHdr (&hdr);
        getHeader (im, &hdr);
        if (hstio_err())
            return (OPEN_FAILED);
        closeImage (im);
        kw = insertfirst (&hdr);

        while ((kw = next (kw)) != NotFound) {
            kw_name = getKwName (kw);
            if (isInvalid (kw_name))
                continue;
            kw_type = getKwType (kw);
            kw_comm = getKwComm (kw);
            switch (kw_type) {
            case FITSLOGICAL:
                kw_bool = getBoolKw (kw);
                c_tbhadb (table->tp, kw_name, kw_bool);
                break;
            case FITSLONG:
                kw_int = getIntKw (kw);
                c_tbhadi (table->tp, kw_name, kw_int);
                break;
            case FITSDOUBLE:
                kw_double = getDoubleKw (kw);
                c_tbhadd (table->tp, kw_name, kw_double);
                break;
            case FITSCHAR:
                getStringKw (kw, kw_string, SZ_KWD);
                c_tbhadt (table->tp, kw_name, kw_string);
                break;
            case FITSNOVALUE:
            case FITSBIT:
            case FITSBYTE:
            case FITSSHORT:
            case FITSFLOAT:
            case FITSCOMPLEX:
            case FITSICOMPLEX:
            case FITSDCOMPLEX:
            case FITSVADESC:
                /* Not handled */
                break;
            }
            c_tbhpcm (table->tp, kw_name, kw_comm);
        }
        freeHdr (&hdr);

        return (0);
}

static int isInvalid (char *kw) {

        /*
         *  When adding additional invalid keywords, make sure they are
         *  added before the NULL string.
         */
        static char *invalid_kwds[] = {
            "BUNIT", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2",
            "CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
            "CUNIT1", "CUNIT2", NULL};

        int j;

        for (j=0; invalid_kwds[j] != NULL; j++) {
            if (strcmp(kw, invalid_kwds[j]) == 0)
                return (1);
        }
        return (0);
}
