# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
#include "hstcal.h"
# include "ctables.h"

/* Image file name template descriptor.
   NOTE:  This version handles only explicit file names, no wildcards.
   More than one file name may be given, separated by a blank or comma
   (or both).  If a file name includes an expression in brackets at the
   end of the name, only the portion before the '[' will be used when
   checking whether the file exists.

   example:

    r_imt = c_imtopen (rawlist);
    n_raw = c_imtlen (r_imt);
    j = c_imtgetim (r_imt, rawfile, STIS_LINE);
    c_imtrew (r_imt);
    c_imtclose (r_imt);
*/

typedef struct {
        char *pattern;          /* file name as specified by the user */
        char **files;           /* array of file names matching the template */
        int nfiles;             /* number of file names in 'files' */
        int alloc_nfiles;       /* allocated length of 'files' array */
        int current_index;      /* index of current file in 'files' */
} ImtDescr;

static void findFiles (ImtDescr *imt_descr);

IRAFPointer c_imtopen (char *pattern) {

/* Create a file name template object.
argument:
char *pattern           i: name of file; wildcard characters are supposed to be
                           allowed, but currently only one explicit file name
                           may be given
function value          o: file name template descriptor
*/

        ImtDescr *imt_descr;
        IRAFPointer imt;

        imt_descr = (ImtDescr *)calloc (1, sizeof(ImtDescr));
        imt_descr->pattern = (char *)calloc (strlen(pattern)+1, sizeof(char));
        strcpy (imt_descr->pattern, pattern);

        /* allocate and populate the list of file names */
        findFiles (imt_descr);

        imt_descr->current_index = 0;

        imt = (void *)imt_descr;

        return imt;
}

int c_imtlen (IRAFPointer imt) {

/* Return the number of file names that match the template.
argument:
IRAFPointer imt         i: file name template descriptor

function value          o: number of file names in the list
*/

        ImtDescr *imt_descr;
        int nfiles;

        if (imt == NULL) {
            nfiles = 0;
        } else {
            imt_descr = (ImtDescr *)imt;
            nfiles = imt_descr->nfiles;
        }
        return nfiles;
}

void c_imtrew (IRAFPointer imt) {

/* "Rewind" the list, i.e. reset the file name index to 0.
argument:
IRAFPointer imt         i: file name template descriptor
*/

        ImtDescr *imt_descr;

        imt_descr = (ImtDescr *)imt;
        imt_descr->current_index = 0;
}

int c_imtgetim (IRAFPointer imt, char *outstr, int maxch) {

/* Get the next file name in the list.
arguments:
IRAFPointer imt         i: file name template descriptor
char *outstr            o: file name
int maxch               i: maximum length of the string (not incl NULL)

function value          o: length of the file name, or 0 if there are no
                           more names in the list
*/

        ImtDescr *imt_descr;
        int i;                  /* current_index, but not incremented */

        if (imt == NULL) {
            outstr[0] = '\0';
            return 0;
        }

        imt_descr = (ImtDescr *)imt;

        if (imt_descr->current_index + 1 > imt_descr->nfiles) {
            outstr[0] = '\0';
            return 0;
        }

        i = imt_descr->current_index;

        if (strlen (imt_descr->files[i]) <= maxch) {
            strcpy (outstr, imt_descr->files[i]);
        } else {
            setError (ERR_STRING_TOO_LONG,
                        "c_imtgetim:  file name is too long");
            strncpy (outstr, imt_descr->files[i], maxch);
            outstr[maxch] = '\0';
        }

        /* point to the next name in the list */
        imt_descr->current_index += 1;

        return (strlen (imt_descr->files[i]));
}

void c_imtclose (IRAFPointer imt) {

/* Deallocate memory for the imt descriptor.
argument:
IRAFPointer imt         i: file name template descriptor
*/

        if (!imt)
            return;

        ImtDescr * imt_descr = (ImtDescr *)imt;
        if (imt_descr->pattern)
            free (imt_descr->pattern);
        {unsigned i;
        for (i = 0;  i < imt_descr->nfiles;  i++)
        {
            if (imt_descr->files[i])
                free (imt_descr->files[i]);
        }}
        if (imt_descr->files)
            free (imt_descr->files);
        if (imt_descr)
            free (imt_descr);
}

static void findFiles (ImtDescr *imt_descr) {

/* Allocate and populate the list of file names. */

        char *filename;         /* a name from the list */
        char *fullname;         /* name with environment variable resolved */
        int nfiles;             /* allocated size of file list */
        int n;                  /* number of files that exist */
        int i, j;
        int done, done_2, done_3;
        int status = 0;

        /* Count the number of names in a comma- or blank-separated list.
           This can give one more than the actual number of names, e.g. if
           the string begins with a blank or ends with a blank, or if no name
           is given, but the extra (blank) name will be rejected later.
           Also, we want alloc_nfiles to be at least one because we'll use
           this number for the length when allocating the 'files' array.
        */
        nfiles = 0;
        i = 0;
        done = 0;
        while (!done) {
            if (imt_descr->pattern[i] == ',' || imt_descr->pattern[i] == ' ')
                i++;
            else
                done = 1;
        }
        done = 0;
        while (!done) {
            if (imt_descr->pattern[i] == '\0') {
                nfiles++;
                done = 1;
            } else if (imt_descr->pattern[i] == ',' ||
                       imt_descr->pattern[i] == ' ') {
                nfiles++;
                i++;
                /* check for and skip over repeated commas or blanks */
                done_2 = 0;
                while (!done_2) {
                    if (imt_descr->pattern[i] == ',' ||
                        imt_descr->pattern[i] == ' ') {
                        i++;
                    } else {
                        done_2 = 1;
                    }
                }
            } else {
                i++;
            }
        }
        imt_descr->alloc_nfiles = nfiles;

        filename = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));

        imt_descr->files = (char **)calloc (imt_descr->alloc_nfiles,
                                sizeof(char *));

        imt_descr->nfiles = 0;          /* set in loop */

        /* preserve i as an index in pattern through all the following loops */
        i = 0;
        /* skip over leading blanks or commas */
        done = 0;
        while (!done) {
            if (imt_descr->pattern[i] == ',' || imt_descr->pattern[i] == ' ')
                i++;
            else
                done = 1;
        }

        n = 0;
        done = 0;
        while (!done) {

            j = 0;
            done_2 = 0;
            while (!done_2) {
                if (imt_descr->pattern[i] == '\0' ||
                    imt_descr->pattern[i] == ',' ||
                    imt_descr->pattern[i] == ' ') {
                    done_2 = 1;
                    filename[j] = '\0';
                    if (imt_descr->pattern[i] != '\0') {
                        i++;
                        done_3 = 0;
                        while (!done_3) {
                            if (imt_descr->pattern[i] == '\0') {
                                done_3 = 1;
                            } else if (imt_descr->pattern[i] == ',' ||
                                       imt_descr->pattern[i] == ' ') {
                                i++;
                            } else {        /* next file name in list */
                                done_3 = 1;
                            }
                        }
                    }
                } else {
                    filename[j] = imt_descr->pattern[i];
                    i++;
                    j++;
                }
            }
            if (imt_descr->pattern[i] == '\0')
                done = 1;

            fullname = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
            status = c_vfn2osfn (filename, fullname);
            if (status != 0) {
                setError (status, "c_imtopen:  error from c_vfn2osfn");
                free (fullname);
                return;
            }

            /* if a file name was specified, add it to the list */
            if (fullname[0] != '\0' && fullname[0] != ' ') {
                imt_descr->files[n] = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
                strcpy (imt_descr->files[n], fullname);
                n++;
                imt_descr->nfiles = n;
            }
            free (fullname);
        }
        free (filename);
}
