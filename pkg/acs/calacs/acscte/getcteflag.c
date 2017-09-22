# include <stdio.h>
# include <string.h>        /* for strncmp, strcmp */

#include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"        /* defines error codes */
#include "trlbuf.h"

static int checkPCTE (Hdr *, ACSInfo *, int *, int *);
static int checkCCD (Hdr *, ACSInfo *, int *);


/* This routine gets the names of reference images and tables from the
 primary header and checks for dummy pedigree.

 Pey Lian Lim, 2013 Aug 12:
 Separated PCTECORR from ACSCCD.
 **

*/
int GetCTEFlags (ACSInfo *acs, Hdr *phdr) {

    extern int status;

    int missing = 0;    /* true if any calibration file is missing */
    int nsteps = 0;     /* number of calibration steps to perform */

    int GetcteSw (ACSInfo *, Hdr *);

    /* Get the values for the Calibration Switches from the
    **    header for processing.
    */
    if (GetcteSw (acs, phdr) )
        return(status);

    /* Although this should never be called with MAMA data, we
       just want to be safe...
    */
    if (acs->detector != MAMA_DETECTOR) {
        if (checkCCD (phdr, acs, &missing))
            return (status);
    }

    /* Check each reference file that we need. */

    if (checkPCTE (phdr, acs, &missing, &nsteps))
        return (status);

    if (missing) {
        return (status = CAL_FILE_MISSING);
    } else if (nsteps < 1) {
        trlwarn ("No calibration switch was set to PERFORM, ");
        trlwarn ("            or all reference files had PEDIGREE = DUMMY.");
        return (status = NOTHING_TO_DO);
    } else {
        return (status);
    }
}


/* If this step is to be performed, check for the existence of the
   CTE parameters file.  If it exists, get the pedigree and descrip
   keyword values.
*/
static int checkPCTE (Hdr *phdr, ACSInfo *acs, int *missing, int *nsteps) {

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs   i: switches, file names, etc
       int *missing     io: incremented if the file is missing
       int *nsteps      io: incremented if this step can be performed
    */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
    void MissingFile (char *, char *, int *);
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    int checkTabRefPedigree (char * filename, RefTab *, int *);

    /* Are we supposed to do this step? */
    if (acs->pctecorr == PERFORM) {

        if (GetSwitch (phdr, "PCTECORR", &calswitch))
            return (status);

        if (calswitch == COMPLETE) {
            acs->pctecorr = OMIT;
            return (status);
        }

        if (acs->pcteTabNameFromCmd && *acs->pcteTabNameFromCmd != '\0')
        {
            char msgBuffer[CHAR_LINE_LENGTH];
            *msgBuffer = '\0';
            sprintf(msgBuffer, "(pctecorr) USING PCTETAB SPECIFIED BY '--pctetab %s' AND NOT THAT FROM IMAGE HEADER!!!", acs->pcteTabNameFromCmd);
            trlwarn(msgBuffer);
            if ((status = checkTabRefPedigree(acs->pcteTabNameFromCmd, &acs->pcte, &acs->pctecorr)))
                return (status);
        }
        else
        {
            if (GetTabRef (acs->refnames, phdr,
                           "PCTETAB", &acs->pcte, &acs->pctecorr))
                return (status);
        }

        if (acs->pcte.exists != EXISTS_YES)
            MissingFile ("PCTETAB", acs->pcte.name, missing);

        if (acs->pctecorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}


/* If the detector is the CCD, we need the CCD parameters table for
   BIASCORR, DARKCORR, and PHOTCORR.  This routine checks that the table
   exists.

   We also need the table for initializing the error array, but we
   don't have a flag for that step.  That's why we need this table
   regardless of which steps are to be performed.
*/
static int checkCCD (Hdr *phdr, ACSInfo *acs, int *missing) {

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs   i: switches, file names, etc
       int *missing     io: incremented if the table is missing
    */

    extern int status;
    int calswitch;            /* returned by GetTabRef and ignored */
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    void MissingFile (char *, char *, int *);
    int GotFileName (char *);

    if (acs->detector == MAMA_DETECTOR)
        return (status);

    if (GetTabRef (acs->refnames, phdr,
                   "CCDTAB", &acs->ccdpar, &calswitch))
        return (status);

    if (acs->ccdpar.exists != EXISTS_YES) {

        MissingFile ("CCDTAB", acs->ccdpar.name, missing);

    } else if (acs->ccdpar.goodPedigree != GOOD_PEDIGREE) {

        (*missing)++;
        sprintf (MsgText, "CCDTAB `%s' is a dummy table.", acs->ccdpar.name);
        trlerror (MsgText);
    }

    /* Get OSCNTAB here as it applies to all CCD processing as well.
       WJH 6 May 1999
    */
    if (GetTabRef (acs->refnames, phdr,
                   "OSCNTAB", &acs->oscn, &acs->blevcorr))
        return (status);

    if (acs->oscn.exists != EXISTS_YES) {
        if (GotFileName (acs->oscn.name)) {
            MissingFile ("OSCNTAB", acs->oscn.name, missing);
        }
    }

    return (status);
}
