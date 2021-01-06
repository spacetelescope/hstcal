/* These functions initialize values for the ACSInfo structure.
   The values set here are NOT default values, just initialization.

   10-Jun-1998 WJH - Initial ACS version.
   24-Sep-1998 WJH - Removed BIAS_REJ as a keyword.
   07-Oct-1998 WJH - Added routines for getting RefImage and RefTab info.
   11-Sep-2000 WJH - Added support for post-FLASH processing step.
   29-Oct-2001 WJH - replaced 'phot' and 'aper' with 'graph' and 'comp',
                     and 'ccdbias' with 'ccdbias[4]'. removed 'filtcorr'.
   04-Dec-2001 WJH - added expend, along with expstart.
   17-Apr-2002 WJH - removed all references to 'statcorr'.
   12-Dec-2012 PLL - added CTE corrected flash reference file.
   12-Aug-2013 PLL - Tidied up code layout.
   05-Dec-2019 MDD - Added overhead for post-flash, unflashed, and DARKTIME values.
   29-Apr-2020 MDD - Added overhead for atod_saturate value.
   11-May-2020 MDD - Added variable, satmap - the reference image for full-well saturation.
                     Use of this image renders acs->saturate variable OBSOLETE.
*/
#include <string.h>

#include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "acscorr.h"  /* calibration switch names */
# include "hstcal.h"


void ACSInit (ACSInfo *acs) {
    void InitRefTab (RefTab *tab);
    void InitRefImg (RefImage *img);

    int i;

    /* Assign default values. */

    acs->input[0] = '\0';
    acs->output[0] = '\0';
    acs->rootname[0] = '\0';
    acs->printtime = 0;
    acs->verbose = 0;
    acs->nThreads = 1;
    acs->cteAlgorithmGen = 0;
    acs->pcteTabNameFromCmd[0] = '\0';

    acs->det[0] = '\0';
    acs->aperture[0] = '\0';
    acs->filter1[0] = '\0';
    acs->filter2[0] = '\0';
    acs->obstype[0] = '\0';
    acs->jwrotype[0] = '\0';
    acs->detector = UNKNOWN_DETECTOR;
    acs->chip = 1;
    acs->ncombine = 1;
    acs->nimsets = 1;
    acs->members = 1;
    acs->mtype[0] = '\0';

    acs->exptime = 1.;
    acs->expstart = 0.;
    acs->expend = 0.;

    acs->sdqflags = MAX_DQ;  /* 16 bits set */

    acs->subarray = 0;
    acs->bin[0] = 1;
    acs->bin[1] = 1;
    acs->offsetx = 0.;
    acs->offsety = 0.;

    /* MAMA-specific info */
    acs->global_limit = 0.;
    acs->tau = 0.;
    acs->local_limit = 0.;
    acs->expand = 0.;
    acs->globrate = 0.;

    /* CCD-specific info */
    acs->ccdamp[0] = '\0';
    acs->ccdgain = 1.;
    acs->binaxis[0] = 1;
    acs->binaxis[1] = 1;
    for (i=0; i < NAMPS; i++) {
        acs->ccdoffset[i] = 0;
        acs->ccdbias[i] = 0.;
        acs->atodgain[i] = 0.;
        acs->readnoise[i] = 0.;
        acs->blev[i] = 0.;
    }
    acs->ampx = 0;
    acs->ampy = 0;
    acs->atod_saturate = 0;
    acs->saturate = 0.;
    acs->trimx[0] = 0;
    acs->trimx[1] = 0;
    acs->trimy[0] = 0;
    acs->trimy[1] = 0;
    acs->vx[0] = 0;
    acs->vx[1] = 0;
    acs->vy[0] = 0;
    acs->vy[1] = 0;
    acs->biassecta[0] = 0;
    acs->biassecta[1] = 0;
    acs->biassectb[0] = 0;
    acs->biassectb[1] = 0;
    acs->flashdur = 0;
    acs->flashstatus[0] = '\0';
    acs->overhead_postflashed = 0.;
    acs->overhead_unflashed = 0.;
    acs->darktime = 0.;

    /* Initialize Calibration switches */
    acs->dqicorr = OMIT;
    acs->atodcorr = OMIT;
    acs->blevcorr = OMIT;
    acs->biascorr = OMIT;
    acs->noisecorr = OMIT;
    acs->glincorr = OMIT;
    acs->lflgcorr = OMIT;
    acs->darkcorr = OMIT;
    acs->flatcorr = OMIT;
    acs->flashcorr = OMIT;
    acs->pctecorr = OMIT;
    acs->pfltcorr = OMIT;
    acs->dfltcorr = OMIT;
    acs->lfltcorr = OMIT;
    acs->shadcorr = OMIT;
    acs->photcorr = OMIT;

    /* Initialize reference images and tables for ACSCCD */
    InitRefImg (&(acs->bias));
    InitRefTab (&(acs->bpix));
    InitRefTab (&(acs->ccdpar));
    InitRefTab (&(acs->oscn));
    InitRefTab (&(acs->atod));
    InitRefTab (&(acs->pcte));
    InitRefImg (&(acs->satmap));

    /* Initialize reference images and tables for ACS2D */
    InitRefImg (&(acs->dark));
    InitRefImg (&(acs->darkcte));
    InitRefImg (&(acs->flash));
    InitRefImg (&(acs->flashcte));
    InitRefImg (&(acs->pflt));
    InitRefImg (&(acs->dflt));
    InitRefImg (&(acs->lflt));
    InitRefImg (&(acs->cflt));
    InitRefImg (&(acs->shad));
    InitRefTab (&(acs->mlin));
    InitRefTab (&(acs->phot));
}


void InitRefTab (RefTab *tab) {
    tab->name[0] = '\0';
    tab->pedigree[0] = '\0';
    tab->descrip[0] = '\0';
    tab->descrip2[0] = '\0';
    tab->exists = EXISTS_UNKNOWN;
    tab->goodPedigree = PEDIGREE_UNKNOWN;
}


void InitRefImg (RefImage *img) {
    img->name[0] = '\0';
    img->pedigree[0] = '\0';
    img->descrip[0] = '\0';
    img->exists = EXISTS_UNKNOWN;
    img->goodPedigree = PEDIGREE_UNKNOWN;
}


/* This routine gets the name of a reference image, checks whether it
   exists, and gets pedigree and descrip if they are present.  The image
   name will be null if the keyword is not present in the header; this is
   not an error.
*/

int GetImageRef (RefFileInfo *refnames, Hdr *phdr,
                 char *keyword, RefImage *image, int *calswitch) {
    extern int status;

    int GetRefName (RefFileInfo *, Hdr *, char *, char *);
    int ImgPedigree (RefImage *);

    /* Get the reference image name. */
    if (GetRefName (refnames, phdr, keyword, image->name))
        return (status);

    /* ImgPedigree opens the image to verify that it exists, and if so,
       gets pedigree & descrip.
    */
    if (ImgPedigree (image))
        return (status);
    if (image->exists == EXISTS_YES) {
        if (image->goodPedigree != GOOD_PEDIGREE)
            *calswitch = DUMMY;
    }

    return (status);
}


/* This routine gets the name of a reference table, checks whether it
   exists, and gets pedigree and descrip if they are present.  The table
   name will be null if the keyword is not present in the header; this is
   not an error.
*/
int GetTabRef (RefFileInfo *refnames, Hdr *phdr,
               char *keyword, RefTab *table, int *calswitch) {
    extern int status;

    int GetRefName (RefFileInfo *, Hdr *, char *, char *);
    int TabPedigree (RefTab *);

    /* Get the reference table name. */
    if (GetRefName (refnames, phdr, keyword, table->name))
        return (status);

    /* TabPedigree opens the table to verify that it exists, and if so,
       gets pedigree & descrip.
    */
    if (TabPedigree (table))
        return (status);
    if (table->exists == EXISTS_YES) {
        if (table->goodPedigree != GOOD_PEDIGREE)
            *calswitch = DUMMY;
    }

    return (status);
}

int checkTabRefPedigree (char *filename, RefTab *table, int *calswitch)
{
    extern int status;

    int TabPedigree (RefTab *);

    strcpy(table->name, filename);

    /* TabPedigree opens the table to verify that it exists, and if so,
       gets pedigree & descrip.
    */
    if ((status = TabPedigree (table)))
        return (status);
    if (table->exists == EXISTS_YES) {
        if (table->goodPedigree != GOOD_PEDIGREE)
            *calswitch = DUMMY;
    }

    return (status);
}

void MissingFile (char *keyword, char *filename, int *missing) {
    sprintf (MsgText, "%s `%s' not found or can't open.", keyword, filename);
    trlerror (MsgText);
    (*missing)++;
}


void initSwitch (CalSwitch *sw) {
    sw->atodcorr = OMIT;
    sw->biascorr = OMIT;
    sw->blevcorr = OMIT;
    sw->crcorr = OMIT;
    sw->darkcorr = OMIT;
    sw->dqicorr = OMIT;
    sw->flatcorr = OMIT;
    sw->flashcorr = OMIT;
    sw->pctecorr = OMIT;
    sw->glincorr = OMIT;
    sw->lflgcorr = OMIT;
    sw->photcorr = OMIT;
    sw->rptcorr = OMIT;
    sw->shadcorr = OMIT;
}
