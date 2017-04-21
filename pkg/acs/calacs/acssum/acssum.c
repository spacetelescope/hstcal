/* acssum -- Sum repeatobs data

   This file contains:
    AcsSum
    AcsInit
    GetSumKeyInfo
    SumGrps
    PutSumHdrInfo
    SquareErr
    SquareErrLine
    SqrtErr
    RptSumLine
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

# include "hstio.h"

# include "acs.h"
# include "acssum.h"
# include "err.h"

static int GetSumKeyInfo (AcsSumInfo *, Hdr *);
static int PutSumHdrInfo (SingleGroup *, double, double, int, int);
static int RptSumLine (SingleGroup *, int, SingleGroupLine *);
static void SqrtErr (SingleGroup *);
static void SquareErr (SingleGroup *);
static void SquareErrLine (SingleGroupLine *);
static void AcsInit (AcsSumInfo *, int);
static int SumGrps (AcsSumInfo *, char *mtype);
static void FreeAcsInput (char **, int);

/* This routine sums together the IMSETs chips of Repeat-Obs exposures.

    Warren Hack, 1998 Oct 12:
       Initial version.

    WJH, 1999 Apr 19:
        Revised SumGrps and GetSumKeyInfo to read EXPTIMEs and EXPEND
            from PHDR headers.
    WJH, 1999 Nov 11:
        Added creation of SPT file for output product.

    WJH 19 April 2001
        Removed 'save_tmp' as a parameter to ACSSUM.  This parameter will
        only be controlled by CALACS, not by individual tasks.
    WJH 17 April 2002
        Removed all references to 'STATFLAG' and 'statcorr'.
        Always perform 'doStat'.
    WJH 6 June 2002
        Corrected an error in allocating space for list of image names in
        InitSumTrl.
*/

int AcsSum (char *input, char *output, char *mtype, int printtime, int verbose)
{

    extern int status;

    IRAFPointer tpin;

    AcsSumInfo acs;
    IODescPtr im;        /* descriptor for input image */
    Hdr phdr;        /* primary header for input image */
    int nimgs;
    int i;
    char acs_input[ACS_FNAME];

    int          FileExists (char *);
    void         TimeStamp (char *, char *);
    void         PrBegin (char *);
    void         PrEnd (char *);
    void         PrFileName (char *, char *);
    void         PrHdrInfo (char *, char *, char *, char *);
    int          MkName (char *, char *, char *, char *, char *, int);
    void         WhichError (int);
    void         InitSumTrl (char *input, char *output);
    void         FindAsnRoot (char *, char *);
    int          mkNewSpt (char *, char *, char *);

    /* Determine input and output trailer files, then initialize
        output file by combining inputs into output file */
    InitSumTrl (input, output);

    PrBegin ("ACSSUM");
    nimgs = 0;

    if (printtime)
        TimeStamp ("ACSSUM started", "");

    /* open the input file template */
    tpin = c_imtopen (input);

    nimgs = c_imtlen(tpin);

    /* Initialize structure containing acssum information. */
    AcsInit (&acs,nimgs);

    /* Copy command-line arguments into acs. */
    for (i = 0; i < nimgs; i++) {
        c_imtgetim (tpin, acs.input[i], ACS_FNAME);
        PrFileName ("input", acs.input[i]);
    }

    /* close file template */
    c_imtclose (tpin);

    strcpy (acs.output, output);
    acs.printtime = printtime;
    acs.verbose = verbose;

    PrFileName ("output", acs.output);
    FindAsnRoot (output, acs.rootname);

    initHdr (&phdr);

    /* Check whether the output file already exists. */
    if (FileExists (acs.output)) {
        FreeAcsInput (acs.input, nimgs);
        return (status);
    }
    strcpy(acs_input,acs.input[0]);

    /* Open input image in order to read its primary header. */
    im = openInputImage (acs_input, "", 0);

    if (hstio_err()) {
        FreeAcsInput (acs.input, nimgs);
        return (status = OPEN_FAILED);
    }

    getHeader (im, &phdr);        /* get primary header */
    if (hstio_err()) {
        FreeAcsInput (acs.input, nimgs);
        return (status = OPEN_FAILED);
    }
    closeImage (im);

    /* Get keyword values from primary header. */
    if (GetSumKeyInfo (&acs, &phdr)) {
        FreeAcsInput (acs.input, nimgs);
        return (status);
    }
    freeHdr (&phdr);

    /* Print information about this image. */
    PrHdrInfo (acs.aperture, acs.filter1, acs.filter2, acs.det);

    if (acs.printtime)
        TimeStamp ("Begin processing", acs.rootname);

    /* Sum all imsets. */
    if (SumGrps (&acs,mtype)){
        FreeAcsInput (acs.input, nimgs);
        return (status);
    }

    /* create new SPT file for output product */
    if (mkNewSpt (input, mtype, output)) {
        return(status);
    }

    /* Done... */
    trlmessage ("\n");
    PrEnd ("ACSSUM");

    if (acs.printtime)
        TimeStamp ("ACSSUM completed", acs.rootname);

    /* Write out temp trailer file to final file */
    WriteTrlFile ();

    FreeAcsInput (acs.input, nimgs);
    return (status);
}

static void FreeAcsInput (char **input, int nimgs) {
    int i;

    for (i = 0; i < nimgs; i++)
        free((char *)input[i]);

    free((char *)input);

}

/* Initialize the AcsSumInfo structure.  This includes information about the
   input and output images, calibration files, and flags to specify which
   calibration steps are to be done.
*/

static void AcsInit (AcsSumInfo *acs, int nimages) {

    int i;

    acs->input = (char **) calloc (nimages, sizeof(char *));
    for (i=0;i<nimages; i++) {
        acs->input[i] = (char *) calloc (ACS_LINE, sizeof(char));
        /* Initialize here ... */
        acs->input[i][0] = '\0';
    }

    /* Assign default values. */
    acs->output[0] = '\0';
    acs->rootname[0] = '\0';
    acs->obsmode[0] = '\0';
    acs->aperture[0] = '\0';
    acs->filter1[0] = '\0';
    acs->filter2[0] = '\0';
    acs->det[0] = '\0';
    acs->detector = UNKNOWN_DETECTOR;
    acs->nimages = nimages;
    acs->nimsets = 0;
    acs->exptime = 0;
    acs->expend = 0;
    acs->sdqflags = MAX_DQ;
}

/* This routine gets keyword values from the primary header.
*/

static int GetSumKeyInfo (AcsSumInfo *acs, Hdr *phdr) {

/* arguments:
AcsSumInfo *acs  io: calibration switches and info
Hdr *phdr       i: primary header
*/

    extern int status;
    int sdqflags;            /* "serious" data quality flags */
    int nextend;            /* number of FITS extensions */
    int nrptexp;            /* number of exposures */
    int no_def = 0;            /* missing keyword is fatal error */

    int GetKeyInt (Hdr *, char *, int, int, int *);
    int GetKeyStr (Hdr *, char *, int, char *, char *, int);
    int GetKeyDbl (Hdr *, char *, int, double, double *);

    /* Have the data been converted from counts to absolute flux? */

    if (GetKeyStr (phdr, "APERTURE", no_def, "", acs->aperture, ACS_CBUF))
        return (status);

    if (GetKeyStr (phdr, "DETECTOR", no_def, "", acs->det, ACS_CBUF))
        return (status);
    if (strcmp (acs->det, "SBC") == 0) {
        acs->detector = MAMA_DETECTOR;
    } else if (strcmp (acs->det, "HRC") == 0) {
        acs->detector = HRC_CCD_DETECTOR;
    } else if (strcmp (acs->det, "WFC") == 0) {
        acs->detector = WFC_CCD_DETECTOR;
    } else {
        sprintf (MsgText, "DETECTOR = %s is invalid", acs->det);
        trlerror (MsgText);
        return (status = HEADER_PROBLEM);
    }

    /* Filter names. */
    if (GetKeyStr (phdr, "FILTER1", USE_DEFAULT, "", acs->filter1, ACS_CBUF))
        return (status);
    if (GetKeyStr (phdr, "FILTER2", USE_DEFAULT, "", acs->filter2, ACS_CBUF))
        return (status);

    if (GetKeyDbl (phdr, "EXPTIME", NO_DEFAULT, 0., &acs->exptime))
        return (status);

    if (GetKeyInt (phdr, "SDQFLAGS", USE_DEFAULT, MAX_DQ, &sdqflags))
        return (status);
    acs->sdqflags = sdqflags;

    /* Find out how many extensions there are in this file. */
    if (GetKeyInt (phdr, "NEXTEND", USE_DEFAULT, EXT_PER_GROUP, &nextend))
        return (status);

    /* Convert number of extensions to number of SingleGroups. */
    acs->nimsets = nextend / EXT_PER_GROUP;

    /* Get NRPTEXP and compare with nimages. */
    if (GetKeyInt (phdr, "NRPTEXP", no_def, 0, &nrptexp))
        return (status);

    if (nrptexp != acs->nimages) {
        sprintf (MsgText, "%d exposures will be summed, however, NRPTEXP was %d.",
            acs->nimages, nrptexp);
        trlwarn (MsgText);
    }

    if (acs->nimages <= 1) {
        trlwarn ("ACSSUM: Can not continue. There is only one input image!");
        return (status = NOTHING_TO_DO);
    }

    return (status);
}

static int SumGrps (AcsSumInfo *acs, char *mtype) {

    extern int status;
    SingleGroup x;                /* first imset */
    SingleGroupLine y;            /* line from Nth imset */
    double exptime;                /* exposure time of current image */
    double sumexptime = 0.;        /* accumulated exposure time */
    char *message;                 /* for printtime info */
    int extver;                    /* imset number */
    int i;                    /* counter for current image */
    int chip, ychip;            /*Chip being summed */
    int extchip;            /* Extension of chip being summed */
    int line;                /* Line of chip being summed */
    char        uroot[ACS_FNAME];   /* Upper case version of rootname */

    int doStat (SingleGroup *, short);
    void TimeStamp (char *, char *);
    void PrGrpBegin (char *, int);
    void PrGrpEnd (char *, int);
    void PrSwitch (char *, int);
    void UCalVer (Hdr *);
    void UFilename (char *, Hdr *);
    void UMemType (char *, Hdr *);
    void UExpname (char *, Hdr *);
    int DetCCDChip (char *, int, int, int *);
    void    UpperAll (char *, char *, int);

    int GetKeyInt (Hdr *, char *, int, int, int *);
    int GetKeyDbl (Hdr *, char *, int, double, double *);
    int PutKeyStr (Hdr *, char *, char *, char *);

    initSingleGroup (&x);
    initSingleGroupLine (&y);

    if (acs->printtime) {
        if ((message = calloc (ACS_LINE+1, sizeof (char))) == NULL)
            return (status = OUT_OF_MEMORY);
    }

    for (extver = 1;  extver <= acs->nimsets;  extver++) {

        PrGrpBegin ("imset", extver);

        getSingleGroup (acs->input[0], extver, &x);
        if (hstio_err())
            return (status = OPEN_FAILED);
        if (acs->printtime)
            TimeStamp ("first imset read", acs->input[0]);

        /* get from x */
        if (GetKeyInt (&x.sci.hdr, "CCDCHIP", USE_DEFAULT, 1, &chip))
            return (status);
        if (GetKeyDbl (x.globalhdr, "EXPEND", NO_DEFAULT, 0., &acs->expend))
            return (status);

        sumexptime = acs->exptime;

        /* Square the errors to convert to variance. */
        SquareErr (&x);                /* operate on x */

        /* For each imset/extver, loop over all images */
        for (i = 1; i < acs->nimages; i++) {

            /* Determine which extension corresponds to desired chip
            ** for the remainder of the images.
            */
            extchip = 0;

            if (DetCCDChip(acs->input[i], chip, acs->nimsets, &extchip) ) {
                return (status);
            }

            /* Get the first line of bias image data. */
            openSingleGroupLine (acs->input[i], extchip, &y);
            if (hstio_err())
                return (status = OPEN_FAILED);

            /* Update exposure time info. */
            /* get from y */

            if (GetKeyInt (&y.sci.hdr, "CCDCHIP", USE_DEFAULT, 1, &ychip))
                return (status);
            if (GetKeyDbl (y.globalhdr, "EXPTIME", NO_DEFAULT, 0., &exptime))
                return (status);
            if (GetKeyDbl (y.globalhdr, "EXPEND", NO_DEFAULT, 0., &acs->expend))
                return (status);

            sumexptime += exptime;

            /*Loop over lines in each subsequent image */
            for (line = 0; line < x.sci.data.ny; line++) {
                status = getSingleGroupLine (acs->input[i], line, &y);
                if (status) {
                    sprintf(MsgText,"Could not read line %d from image %d.",line+1,i+1);
                    trlerror(MsgText);
                    return (status = OPEN_FAILED);
                }

                SquareErrLine (&y);                /* operate on y */

                /* Add current imset to sum (i.e. add y to x).  This differs
                    from add2d in that RptSum adds variances, rather than
                    adding errors in quadrature.
                */
                if (RptSumLine (&x, line, &y))
                    return (status);

            } /*End loop over lines */

            if (acs->printtime) {
                if (i == 1)
                    strcpy (message, "1st imset added");
                else if (i == 2)
                    strcpy (message, "2nd imset added");
                else if (i == 3)
                    strcpy (message, "3rd imset added");
                else
                    sprintf (message, "%dth imset added", i);

                TimeStamp (message, acs->input[i]);
            }

            closeSingleGroupLine (&y);
        } /* End loop over images */

        freeSingleGroupLine (&y);

        /* Take the square root of variance to convert back to errors. */
        SqrtErr (&x);

        /* Compute statistics and update keywords in output headers. */
        trlmessage ("\n");
        if (doStat (&x, acs->sdqflags))
        return (status);
        if (acs->printtime)
            TimeStamp ("Image statistics computed", acs->rootname);
        /* Update header info in output. */
        if (PutSumHdrInfo (&x, sumexptime, acs->expend, acs->nimages, acs->nimsets))
            return (status);

        /* Update CAL_VER and FILENAME, then write output file.
            EXPNAME values modified for all extensions in a SingleGroup.
                    WJH 7 July 1999
        */
        UCalVer (x.globalhdr);
        UFilename (acs->output, x.globalhdr);
        UMemType (mtype, x.globalhdr);
        UExpname (acs->rootname, &x.sci.hdr);
        UExpname (acs->rootname, &x.err.hdr);
        UExpname (acs->rootname, &x.dq.hdr);
        UpperAll (acs->rootname, uroot, strlen(acs->rootname)+1 );

        PutKeyStr (x.globalhdr, "ROOTNAME", uroot,"Rootname of the observation set");

        putSingleGroup (acs->output, extver, &x, 0);
        if (hstio_err())
                return (status = 1001);
        freeSingleGroup (&x);

        PrGrpEnd ("imset", extver);

        if (acs->printtime)
            TimeStamp ("Output written to disk", acs->rootname);
    } /* End loop over imsets */

    if (acs->printtime)
        free (message);
    return (status);
}

/* This routine adds history info and updates RPTCORR, NEXTEND, and
   NCOMBINE in the primary header.  It also updates EXPTIME and EXPEND
   in the Global extension header.
*/

static int PutSumHdrInfo (SingleGroup *out, double sumexptime, double expend, int nimages, int nimsets) {

/* arguments:
SingleGroup *out  i: current imset
double *exptime   i: accumulated exposure time
double *expend    i: last exposure end time read from input file
int nimages       i: the number of imsets that were combined
*/

    extern int status;
    int PutKeyDbl (Hdr *, char *, double, char *);
    int PutKeyInt (Hdr *, char *, int, char *);
    int PutKeyStr (Hdr *, char *, char *, char *);

    /* Set the switch to indicate that rptcorr has been done. */
    if (PutKeyStr (out->globalhdr, "RPTCORR", "COMPLETE",
            "add individual repeat observations"))
        return (status);

    /* Write history records. */

    addHistoryKw (out->globalhdr, "RPTCORR complete");
    if (hstio_err())
        return (status = HEADER_PROBLEM);

    addHistoryKw (out->globalhdr, "Statistics computed after rptcorr.");
    if (hstio_err())
        return (status = HEADER_PROBLEM);

    /* Update NEXTEND in primary header, to indicate only one imset. */
    if (PutKeyInt (out->globalhdr, "NEXTEND", EXT_PER_GROUP * nimsets,
            "number of extensions"))
        return (status);

    /* Update NCOMBINE in primary header, to tell how many imsets were
       combined into this one output imset.
    */
    if (PutKeyInt (out->globalhdr, "NCOMBINE", nimages,
            "number of imsets combined"))
        return (status);

    /* Update exposure time info in SCI extension header. */
    if (PutKeyDbl (out->globalhdr, "EXPTIME", sumexptime, "exposure time"))
        return (status);
    if (PutKeyDbl (out->globalhdr, "EXPEND", expend, "exposure end time"))
        return (status);

    return (status);
}

/* convert error to variance */

static void SquareErr (SingleGroup *x) {

    int i, j;

    for (j = 0;  j < x->err.data.ny;  j++) {
        for (i = 0;  i < x->err.data.nx;  i++) {
        Pix (x->err.data, i, j) =
            Pix (x->err.data, i, j) * Pix (x->err.data, i, j);
        }
    }
}
static void SquareErrLine (SingleGroupLine *y) {

    int i;

    for (i = 0;  i < y->err.tot_nx;  i++) {
        y->err.line[i] = y->err.line[i] * y->err.line[i];
    }

}

/* convert variance to error */

static void SqrtErr (SingleGroup *x) {

    int i, j;

    for (j = 0;  j < x->err.data.ny;  j++) {
        for (i = 0;  i < x->err.data.nx;  i++) {
        Pix (x->err.data, i, j) = sqrt (Pix (x->err.data, i, j));
        }
    }
}

/* Add one SingleGroup triplet with a SingleGroupLine,
    leaving the result in the first.

   (*a) += (*b[i])

   The science data arrays are added together; the error arrays are
   added; the data quality arrays are ORed.

   This differs from add2d in that RptSum[Line] assumes the error arrays
   contain variance rather than standard deviations, so those values
   will simply be added.
*/

static int RptSumLine (SingleGroup *a, int line, SingleGroupLine *b) {

/* arguments:
SingleGroup *a      io: input data; output sum
int line             i: line from input/output to be summed
SingleGroupLine *b   i: second input data line
*/

    extern int status;

    int i;
    short dqa, dqb, dqab;    /* data quality for a, b, combined */

    if (a->sci.data.nx != b->sci.tot_nx)
        return (status = SIZE_MISMATCH);

    /* science data */
    for (i = 0;  i < a->sci.data.nx;  i++) {
        Pix (a->sci.data, i, line) = Pix(a->sci.data, i, line) + b->sci.line[i];
    }

    /* error array (actually contains variance) */
        for (i = 0;  i < a->err.data.nx;  i++) {
        Pix (a->err.data, i, line) = Pix (a->err.data, i, line) + b->err.line[i];
        }

    /* data quality */
        for (i = 0;  i < a->dq.data.nx;  i++) {
        dqa = DQPix (a->dq.data, i, line);
        dqb = b->dq.line[i];
        dqab = dqa | dqb;
        DQSetPix (a->dq.data, i, line, dqab);
        }

    return (status);
}



void InitSumTrl (char *input, char *output) {

    extern int status;
    IRAFPointer    tpin;
    int n;

    char *trl_in;     /* trailer filename for input */
    char trl_out[ACS_LINE+1];     /* output trailer filename */
    char in_name[ACS_FNAME+1];
    char out_name[ACS_FNAME+1];

    int trl_len;

    char *isuffix[] = {"_crj", "_flt"};
    char *osuffix[] = {"_sfl", "_sfl"};
    char *trlsuffix[] = {"", ""};
    int nsuffix = 2;

    int MkOutName (char *, char **, char **, int, char *, int);
    int MkNewExtn (char *, char *);
    void WhichError (int);

    trl_in = realloc (NULL, (ACS_LINE));
    trl_len = ACS_LINE;

    if (trl_in == NULL) {
        trlerror ("Out of memory: Couldn't allocate for CRJ_TMP trailer file.");
         status = OUT_OF_MEMORY;
        trl_len = 0;
     }
    /* Initialize TRL filenames */
    trl_in[0] = '\0';
    trl_out[0] = '\0';

    /* open the input file template */
    tpin = c_imtopen (input);

     for (n = 0; n < c_imtlen(tpin); ++n) {
        c_imtgetim (tpin, in_name, ACS_FNAME);
        out_name[0] = '\0';


        /* Start by stripping off suffix from input/output filenames */
        if (MkOutName (in_name, isuffix, trlsuffix, nsuffix, out_name, ACS_LINE)) {
            WhichError (status);
            sprintf (MsgText, "Couldn't create trailer filename for %s", in_name);
            trlerror (MsgText);
            continue;
        }

        /* Now, convert trailer filename extensions from '.fits' to '.trl' */
        if (MkNewExtn (out_name, TRL_EXTN) ) {
            sprintf(MsgText, "Error with input trailer filename %s", out_name);
            trlerror (MsgText);
            WhichError (status);
        }

        if ( (strlen(out_name) + strlen(trl_in) + 1) >= trl_len) {
            /*
                Add 1 to out_name to account for comma to be appended.
                WJH  4 June 2002
            */
            trl_len += strlen(out_name) + 1;
            trl_in = realloc (trl_in, trl_len);
        }

        /* Append each filename to create list of input trailer files */
        strcat(trl_in, out_name);

        /* But don't put a comma after the last filename */
        if (n < (c_imtlen(tpin)-1)) strcat (trl_in, ",");

    }

    if (MkOutName (output, osuffix, trlsuffix, nsuffix, trl_out, ACS_LINE)) {
        WhichError (status);
        sprintf (MsgText, "Couldn't create trailer filename for %s", output);
        trlerror (MsgText);
    }

    /* Now, convert trailer filename extensions from '.fits' to '.trl' */
    if (MkNewExtn (trl_out, TRL_EXTN) ) {
        sprintf(MsgText, "Error with input trailer filename %s", trl_out);
        trlerror (MsgText);
        WhichError (status);
    }
    /* Sets up temp trailer file for output and copies input
        trailer file into it.
    */
    InitTrlFile (trl_in, trl_out);

    /* Deallocate memory */
    free(trl_in);
    c_imtclose (tpin);
}
