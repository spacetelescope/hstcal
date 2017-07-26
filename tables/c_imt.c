#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#  include <glob.h>
#  define WXYZ_ZYXW_Unix
#endif
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

# define DELTA_MAXFILES 100

/* Expand a file name template.
   More than one file name may be given, separated by a blank or comma
   (or both).  For Unix-like systems, the file names may include wildcard
   characters.  If a file name includes an expression in brackets at the
   end of the name (e.g. [sci,1], for specifying a FITS extension), that
   expression will be stripped off before expanding wildcards, and then
   it will be appended back on to each matching file name.

   example:

    IRAFPointer r_imt;
    int n_raw, nchar;

    r_imt = c_imtopen(rawlist);
    n_raw = c_imtlen(r_imt);
    nchar = c_imtgetim(r_imt, rawfile, STIS_LINE);
    c_imtrew(r_imt);
    c_imtclose(r_imt);
*/

typedef struct {
        char *pattern;          /* file name as specified by the user */
        char **files;           /* array of file names matching the template */
        int nfiles;             /* number of file names in 'files' */
        int alloc_nfiles;       /* allocated length of 'files' array */
        int current_index;      /* index of current file in 'files' */
} ImtDescr;

static int findFiles(ImtDescr *imt_descr);
static int add_filename(ImtDescr *imt_descr, char *filename);
static int more_names(ImtDescr *imt_descr, int new_maxfiles);

IRAFPointer c_imtopen(char *pattern) {

/* Create a file name template object.
argument:
char *pattern           i: name of file; wildcard characters are supposed to be
                           allowed, but currently only one explicit file name
                           may be given
function value          o: file name template descriptor
*/

        ImtDescr *imt_descr;
        IRAFPointer imt;
        int status;

        imt_descr = (ImtDescr *)calloc(1, sizeof(ImtDescr));
        if (imt_descr == NULL) {
            return NULL;
        }
        imt_descr->pattern = (char *)calloc(strlen(pattern)+1, sizeof(char));
        if (imt_descr->pattern == NULL) {
            free(imt_descr);
            return NULL;
        }
        imt_descr->files = NULL;
        imt_descr->nfiles = 0;
        imt_descr->alloc_nfiles = 0;
        imt_descr->current_index = 0;

        strcpy(imt_descr->pattern, pattern);

        /* Populate the list of file names. */
        status = findFiles(imt_descr);
        if (status != 0) {
            return NULL;
        }

        imt = (void *)imt_descr;

        return imt;
}

int c_imtlen(IRAFPointer imt) {

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

void c_imtrew(IRAFPointer imt) {

/* "Rewind" the list, i.e. reset the file name index to 0.
argument:
IRAFPointer imt         i: file name template descriptor
*/

        ImtDescr *imt_descr;

        imt_descr = (ImtDescr *)imt;
        imt_descr->current_index = 0;
}

int c_imtgetim(IRAFPointer imt, char *outstr, int maxch) {

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

        if (strlen(imt_descr->files[i]) <= maxch) {
            strcpy(outstr, imt_descr->files[i]);
        } else {
            setError(ERR_STRING_TOO_LONG,
                        "c_imtgetim:  file name is too long");
            strncpy(outstr, imt_descr->files[i], maxch);
            outstr[maxch] = '\0';
        }

        /* point to the next name in the list */
        imt_descr->current_index += 1;

        return strlen(imt_descr->files[i]);
}

void c_imtclose(IRAFPointer imt) {

/* Deallocate memory for the imt descriptor.
argument:
IRAFPointer imt         i: file name template descriptor
*/

        ImtDescr *imt_descr;
        int i;

        if (imt != NULL) {
            imt_descr = (ImtDescr *)imt;
            free(imt_descr->pattern);
            for (i = 0;  i < imt_descr->nfiles;  i++)
                free(imt_descr->files[i]);
            free(imt_descr->files);
            free(imt_descr);
        }
}

static int findFiles(ImtDescr *imt_descr) {

/* Allocate and populate the list of file names. */

#ifdef WXYZ_ZYXW_Unix
        glob_t globbuf;
#endif

        char *fullname;         /* name with environment variable resolved */
        char *filename;         /* a file name from the pattern */
        char *ext_spec;         /* expression in brackets at end of name */
        int len_ext_spec;       /* length of expression with brackets */
        int len_fname;          /* length of a file name */
        int max_strlen;         /* maximum length of filename pattern */
        int flags;              /* for the glob function */
        int nfiles;             /* allocated size of file list */
        int i, j;
        int k;
        int brackets;           /* incr. with '[', decremented with ']' */
        int done, done_2, done_3;       /* boolean */
        int status = 0;

        /* Count the number of names in a comma- or blank-separated list,
           and get the length of the longest file name in the list.
        */
        nfiles = 0;
        max_strlen = 0;
        i = 0;
        done = 0;
        while (!done) {
            if (imt_descr->pattern[i] == ',' || imt_descr->pattern[i] == ' ')
                i++;
            else
                done = 1;
        }
        done = 0;
        brackets = 0;
        j = 0;
        while (!done) {
            if (imt_descr->pattern[i] == '\0') {
                nfiles++;
                max_strlen = (i - j > max_strlen) ? (i - j) : max_strlen;
                j = i;
                done = 1;
            } else if ((imt_descr->pattern[i] == ',' ||
                        imt_descr->pattern[i] == ' ') && !brackets) {
                nfiles++;
                max_strlen = (i - j > max_strlen) ? (i - j) : max_strlen;
                j = i;
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
                if (imt_descr->pattern[i] == '[')
                    brackets++;
                else if (imt_descr->pattern[i] == ']')
                    brackets--;
                i++;
            }
        }
        if (brackets != 0) {
            printf("Warning  Unmatched bracket in pattern '%s'\n",
                   imt_descr->pattern);
        }
        if (max_strlen < 1)
            max_strlen = SZ_FNAME;

        filename = (char *)calloc(max_strlen+1, sizeof(char));

        /* Allocate space for nfiles files.  This may be increased later,
           if one or more names in the list include wildcard characters.
        */
        status = more_names(imt_descr, nfiles);
        if (status != 0) {
            setError(status, "c_imtopen:  out of memory");
            return status;
        }

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

        brackets = 0;
        done = 0;
        while (!done) {

            j = 0;
            done_2 = 0;
            while (!done_2) {
                if (imt_descr->pattern[i] == '\0' ||
                    ((imt_descr->pattern[i] == ',' ||
                      imt_descr->pattern[i] == ' ') &&
                     !brackets)) {
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
                    if (imt_descr->pattern[i] == '[')
                        brackets++;
                    else if (imt_descr->pattern[i] == ']')
                        brackets--;
                    i++;
                    j++;
                }
            }
            if (imt_descr->pattern[i] == '\0')
                done = 1;

            fullname = (char *)calloc(SZ_FNAME+1, sizeof(char));
            status = c_vfn2osfn(filename, fullname);
            if (status != 0) {
                setError(status, "c_imtopen:  error from c_vfn2osfn");
                free(fullname);
                return status;
            }

#ifdef WXYZ_ZYXW_Unix
            /* If fullname ends with an expression in brackets (e.g. an
               extension name or number), strip it off before expanding
               wildcards, then append it to each file name.
            */
            len_fname = strlen(fullname);
            ext_spec = NULL;
            if (len_fname > 2 && fullname[len_fname-1] == ']') {
                for (k = len_fname - 1;  k > 0;  k--) {
                    if (fullname[k] == '[') {
                        len_ext_spec = len_fname - k;
                        ext_spec = malloc((len_ext_spec + 10) * sizeof(char));
                        if (ext_spec == NULL) {
                            setError(status, "c_imtopen:  out of memory");
                            return 1;
                        }
                        strcpy(ext_spec, fullname+k);
                        fullname[k] = '\0';     /* chop off "[whatever]" */
                        break;
                    }
                }
            }
            globbuf.gl_offs = 0;
            /* GLOB_MARK --> Append a slash to each path which corresponds
               to a directory.
               GLOB_NOCHECK --> If no pattern matches, return the original
               pattern.  This is needed to allow specifying output files,
               which should not exist yet.
            */
            flags = GLOB_MARK | GLOB_NOCHECK;
            status = glob(fullname, flags, NULL, &globbuf);
            if (status != 0) {
                char msg[SZ_FNAME];
                sprintf(msg, "c_imtopen:  error %d from glob", status);
                setError(status, msg);
                free(fullname);
                free(filename);
                if (ext_spec != NULL)
                    free(ext_spec);
                return status;
            }
            if (ext_spec == NULL) {
                for (k = 0;  k < globbuf.gl_pathc;  k++) {
                    status = add_filename(imt_descr, globbuf.gl_pathv[k]);
                    if (status != 0) {
                        setError(status, "c_imtopen:  out of memory");
                        return status;
                    }
                }
            } else {
                /* re-use fullname */
                for (k = 0;  k < globbuf.gl_pathc;  k++) {
                    strcpy(fullname, globbuf.gl_pathv[k]);
                    strcat(fullname, ext_spec);
                    status = add_filename(imt_descr, fullname);
                    if (status != 0) {
                        setError(status, "c_imtopen:  out of memory");
                        return status;
                    }
                }
                free(ext_spec);
            }
            globfree(&globbuf);
#else

            /* if a file name was specified, add it to the list */
            if (fullname[0] != '\0' && fullname[0] != ' ') {
                status = add_filename(imt_descr, fullname);
            }
#endif

            free(fullname);
        }
        free(filename);

        return status;
}

static int add_filename(ImtDescr *imt_descr, char *filename) {

/* Add one file name to the list in imt_descr.  Memory will be reallocated
   if necessary to allow more room for file names.

arguments:
ImtDescr *imt_descr     i: Pointer to an ImtDescr struct.
char *filename          i: Pointer to a string containing a file name.
function value          o: status; 0 is OK, 1 means out of memory.
*/

        int len;                /* length of the file name */
        int n;                  /* number of file names in imt */
        int status = 0;

        if (imt_descr->nfiles >= imt_descr->alloc_nfiles) {
            status = more_names(imt_descr, -1);
            if (status != 0) {
                return status;
            }
        }

        len = strlen(filename) + 1;
        n = imt_descr->nfiles + 1;
        imt_descr->files[n-1] = malloc(len * sizeof(char));
        if (imt_descr->files[n-1] == NULL) {
            status = 1;
        } else {
            strcpy(imt_descr->files[n-1], filename);
        }
        imt_descr->nfiles = n;

        return status;
}

static int more_names(ImtDescr *imt_descr, int new_maxfiles) {

/* Allocate or reallocate the `files` member of imt_descr.
   If new_maxfiles is less than or equal to zero, it will be set to the
   current maximum number of files plus an increment that is defined in
   this file as a macro.  If the new maximum number of files is less than
   or equal to the current maximum, this function will return without
   changing anything; otherwise, the `files` member of imt_descr will be
   reallocated to increase its length.

arguments:
ImtDescr *imt_descr     i: Pointer to an ImtDescr struct.  Memory may be
                           allocated or reallocated, and the alloc_nfiles
                           member may be updated.
int new_maxfiles        i: The new maximum number of files; the `files` array
                           of pointers will be reallocated (if necessary) to
                           be this length.  If new_maxfiles is negative or
                           zero, the current array length will be increased
                           by a default amount.
function value          o: status; 0 is OK, 1 means out of memory.
*/

        void *x;
        int i;
        int status = 0;

        if (new_maxfiles <= 0) {
            new_maxfiles = imt_descr->alloc_nfiles + DELTA_MAXFILES;
        }
        if (new_maxfiles > imt_descr->alloc_nfiles) {
            x = realloc(imt_descr->files, new_maxfiles * sizeof(char *));

            if (x == NULL) {
                status = 1;
            } else {
                imt_descr->files = x;
                for (i = imt_descr->alloc_nfiles;  i < new_maxfiles;  i++) {
                    imt_descr->files[i] = NULL;
                }
                imt_descr->alloc_nfiles = new_maxfiles;
                status = 0;
            }
        } else {
            status = 0;         /* don't need to reallocate */
        }

        return status;
}
