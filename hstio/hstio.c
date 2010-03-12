/*
** Allen Farris - Original Implementation.
**
** M.D. De La Pena, 25 February 1998: Modified use of "long" to "int" to
** enforce compatibility with SPP/IRAF.  Removed numeric.h - not used.
**
** Version 2.0, 23 March 1998:
** M.D. De La Pena: Added structures SingleGroupLine, FloatHdrLine, and
**      ShortHdrLine.  Also, new functions to support the acquisition of
**      obtaining single lines from a SingleGroup.
**
** 07 April 1998:
** M.D. De La Pena: Added functionality to output a subsection of each image
** of an imset to a file where the subsection is the full size (naxis1/naxis2)
** of the output image.
**
** 26 August 1998:
** M.D. De La Pena: Modified getSingleGroupLine to increment the line number
** by 1.  The public interface should have the line numbers from zero through
** (size of data - 1).  However, the IRAF routines to acquire lines require
** the line numbers to range from one through (size of data).
**
** 30 September 1998:
** M.D. De La Pena: Removed EXTVER parameter from getSingleGroupLine since
** it was unnecessary.
**
** 07 October 1998
** M.D. De La Pena: Modified putSect[Float/Short]HD to update the WCS keywords
** (LTV[1/2] and CRPIX[1/2], if these keywords are present in the header.
** A new supporting routine, updateWCS, performs the update.
**
** 09 November 1998:
** M.D. De La Pena: The following routine names have been updated to append
** "Sect" to the end of the name rather than using it as a prefix - consistency
** request: put[Sci/Err/DQ]Sect, put[Float/Short]HDSect, putSingleGroupSect,
** and put[Float/Short]Sect.
**
** 12 November 1998:
** M.D. De La Pena: Added use of the HSTIO_VERSION macro so that a search can
** be done on the HSTIO library or the executable using the HSTIO library to
** determine the version.  Corrected an omitted call to "initHdr" in
** allocSingleGroup, allocSingleGroupLine, allocMultiGroup, and
** allocMultiNicmosGroup.
**
** 28 December 1998:
** M.D. De La Pena: Modified getSingleGroup and openSingleGroupLine to avoid
** allocating pointers to space of zero length; now agrees with code from
** getSingleNicmosGroup.
**
** 20 May 1999:
** R.L. Williamson: Modified type for updateWCS from int to void.
**
** 19 August 1999: M.D. De La Pena - Added void to functions which have no
** parameters.
**
** 13 January 2000:
** M.D. De La Pena - Modified freeMultiGroup and freeMultiNicmosGroup to
** handle structures which have been initialized, but never allocated before
** they are freed.
**
** 24 February 2000:
** M.D. De La Pena - Modified getSingleGroupLine such that the input line
** number (which should be zero-based) is no longer incremented by one in
** this routine.  The increment by one is done in get[Float/Short]Line.
** Routine initSingleGroupLine() now initializes line_num attribute to -1.
**
** 16 February 2007:
** H.A. Bushouse - Added putSingleNicmosGroupSect routine, patterned after
** the existing putSingleGroupSect, and the supporting putSmplSect and
** putIntgSect routines.
**
** M. Droettboom, January 2010:
** Change to use CFITSIO rather than IRAF IMIO routines.
**
** Table of Contents
**
** Section 1.
**      Defines that isolate data type dependencies between Sci, DQ, and Err.
** Section 2.
**      Declarations and functions related to error handling.
** Section 3.
**      Functions that initialize, allocate, and free storage.
** Section 4 .
**      Low-level I/O functions.
** Section 5.
**      High-level I/O Functions.
** Section 6.
**      Functions to manipulate the header array.
*/
# include "hstio.h"
# include <stdio.h>
# include <string.h>
# include <ctype.h>
# if defined(VMS)
# include <stat.h>
# else
# include <sys/types.h>
# include <sys/stat.h>
# endif

#include <fitsio.h>

/*
** String defined to allow determination of the HSTIO library version
** from the library file (*.a) or the executable using the library.
*/
const char *hstio_version = HSTIO_VERSION;

/*
** Section 1.
** Defines that isolate data type dependencies between Sci, DQ, and Err
** and private data declarations.
*/

typedef struct {
  /* CFITSIO TODO: Remove redundant members here */
        fitsfile *ff;           /* Ptr to cfitsio fitsfile object       */
        char *filename;         /* File name.                           */
        char *extname;          /* FITS EXTNAME value.                  */
        int extver;             /* FITS EXTVER value.                   */
        int hflag;              /* Flag indicating header update.       */
        Hdr *hdr;               /* Address of header lines.             */
        long dims[2];           /* FITS NAXIS values.                   */
        int type;               /* FITS data type.                      */
        unsigned int options;   /* I/O options.                         */
} IODesc;

/* The default allocation unit for Header arrays */
# define HdrUnit 36

/*
** Section 2.
** Declarations and functions related to error handling.
*/
# define ERRLINEWIDTH 2048
static HSTIOError error_status;
static char error_msg[ERRLINEWIDTH];
static HSTIOErrHandler errhandler[32];
static int max_err_handlers = 32;
static int errtop = -1;

HSTIOError hstio_err(void) {
        return error_status;
}

char *hstio_errmsg(void) {
        return error_msg;
}

int push_hstioerr(HSTIOErrHandler x) {
        if (errtop == (max_err_handlers - 1)) return -1;
        ++errtop;
        errhandler[errtop] = x;
        return errtop + 1;
}

int pop_hstioerr(void) {
        if (errtop == -1) return -1;
        --errtop;
        return errtop + 1;
}

static void clear_err(void) { error_status = HSTOK; error_msg[0] = '\0'; }

void clear_hstioerr(void) { error_status = HSTOK; error_msg[0] = '\0'; }

# if defined(__cplusplus)
extern "C" {
void error(HSTIOError, char *);
}
# endif

void error(HSTIOError e, char *str) {
        int n;

        error_status = e;
        if (str != 0) {
            n = strlen(str);
            strncpy(error_msg,str,(n > ERRLINEWIDTH ? ERRLINEWIDTH : n));
            error_msg[n] = '\0';
        }
        switch(error_status) {
            /* Do not make these messages longer than 80 chars. */
            case HSTOK:
                error_msg[0] = '\0';
                break;
            case NOMEM:
                strcat(error_msg,"\nNo memory left to allocate data.");
                break;
            case BADOPEN:
                strcat(error_msg,"\nError opening image array.");
                break;
            case BADCLOSE:
                strcat(error_msg,"\nError closing image array.");
                break;
            case BADREAD:
                strcat(error_msg,"\nError reading image array.");
                break;
            case BADWRITE:
                strcat(error_msg,"\nError writing image array.");
                break;
            case BADEXTNAME:
                strcat(error_msg,"\nInvalid EXTNAME name");
                break;
            case BADHSIZE:
                strcat(error_msg,"\nInvalid size for header array.");
                break;
            case NOGET:
                strcat(error_msg,"\nIncorrect I/O mode for get operation.");
                break;
            case NOPUT:
                strcat(error_msg,"\nIncorrect I/O mode for put operation.");
                break;
            case BADDIMS:
                strcat(error_msg,"\nImage has wrong number of dimensions.");
                break;
            case BADTYPE:
                strcat(error_msg,"\nImage has wrong data type.");
                break;
            case NOSCI:
                strcat(error_msg,"\nNo Sci array corresponding to DQ or Err arrays");
                break;
            case BADSCIDIMS:
                strcat(error_msg,"\nSci array has wrong number of dimensions.");
                break;
            case BADGROUP:
                strcat(error_msg,"\nGroup number is out of range.");
                break;
            case BADGET:
                strcat(error_msg,"\nKeyword specified in get_Kw function was not found.");
                break;
            case BADFITSEQ:
                strcat(error_msg,"\nFITS card has no value indicator.");
                break;
            case BADFITSQUOTE:
                        strcat(error_msg,"\nFITS card has no ending quote.");
                break;
            case BADFITSNUMERIC:
                strcat(error_msg,"\nFITS card has invalid numeric field.");
                break;
            case BADFITSTYPE:
                strcat(error_msg,"\nWrong data type specified in get_Kw function.");
                break;
            case BADPUT:
                strcat(error_msg,"\nKeyword specified in put_Kw function was not found.");
                break;
            case BADNAME:
                strcat(error_msg,"\nKeyword name specified in add_Kw function is too long.");
                break;
            case BADBITPIX:
                strcat(error_msg,"\nWrong data type specified in making primary array or image extension.");
                break;
            case BADNDIM:
                strcat(error_msg,"\nWrong number of dimensions in making primary array or image extension.");
                break;
            case BADEXIST:
                strcat(error_msg,"\nFile already exists.  Operation would overwrite existing file.");
                break;
            case BADREMOVE:
                strcat(error_msg,"\nError removing file.");
                break;
        }
        if (errtop > -1 && errhandler[errtop] != 0)
            errhandler[errtop]();
}

/*
** Section 3.
** Functions that initialize, allocate, and free storage in data structures.
*/
void initFloatData(FloatTwoDArray *x) {
        x->buffer = NULL;
        x->buffer_size = 0;
        x->tot_nx = 0;
        x->tot_ny = 0;
        x->nx = 0;
        x->ny = 0;
        x->data = NULL;
# if defined (DEBUG)
        printf("initFloatData: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
}

int allocFloatData(FloatTwoDArray *x, int i, int j) {
# if defined (DEBUG)
        printf("allocFloatData-1: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
        if (x->buffer == NULL || x->buffer_size != (i * j)) {
            if (x->buffer != NULL)
                free(x->buffer);
            x->buffer_size = i * j;
            x->buffer = (float *)calloc(x->buffer_size, sizeof(float));
            if (x->buffer == NULL) {
                initFloatData(x);
                error(NOMEM,"Allocating SciData");
                return -1;
            }
        }
        x->tot_nx = i;
        x->tot_ny = j;
        x->nx = x->tot_nx;
        x->ny = x->tot_ny;
        x->data = x->buffer;
# if defined (DEBUG)
        printf("allocFloatData-2: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
        return 0;
}

void freeFloatData(FloatTwoDArray *x) {
# if defined (DEBUG)
        printf("freeFloatData: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
        if (x->buffer != NULL)
            free(x->buffer);
        initFloatData(x);
}

void initShortData(ShortTwoDArray *x) {
        x->buffer = NULL;
        x->buffer_size = 0;
        x->tot_nx = 0;
        x->tot_ny = 0;
        x->nx = 0;
        x->ny = 0;
        x->data = NULL;
# if defined (DEBUG)
        printf("initShortData: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
}

int allocShortData(ShortTwoDArray *x, int i, int j) {
# if defined (DEBUG)
        printf("allocShortData-1: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
        if (x->buffer == NULL || x->buffer_size != (i * j)) {
            if (x->buffer != NULL)
                free(x->buffer);
            x->buffer_size = i * j;
            x->buffer = (short *)calloc(x->buffer_size, sizeof(short));
            if (x->buffer == NULL) {
                initShortData(x);
                error(NOMEM,"Allocating DQData");
                return -1;
            }
        }
        x->tot_nx = i;
        x->tot_ny = j;
        x->nx = x->tot_nx;
        x->ny = x->tot_ny;
        x->data = x->buffer;
# if defined (DEBUG)
        printf("allocShortData-2: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
        return 0;
}

void freeShortData(ShortTwoDArray *x) {
# if defined (DEBUG)
        printf("freeShortData: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
        if (x->buffer != NULL)
            free(x->buffer);
        initShortData(x);
}

void initFloatLine (FloatHdrLine *x) {
        x->line   = NULL;
        x->tot_nx = 0;
}

int allocFloatLine (FloatHdrLine *x, int i) {
# if defined (DEBUG)
        printf("allocFloatLine-1: %x %d\n",
                (int)(x->line),i);
# endif
        if (x->line == NULL || x->tot_nx != i) {
            if (x->line != NULL)
                free (x->line);
            x->tot_nx = i;
            x->line = (float *) calloc (x->tot_nx, sizeof(float));
            if (x->line == NULL) {
                initFloatLine (x);
                error (NOMEM,"Allocating float line");
                return (-1);
            }
        }
# if defined (DEBUG)
        printf("allocFloatLine-2: %x %d\n",
                (int)(x->line),x->tot_nx);
# endif
        return (0);
}

void freeFloatLine (FloatHdrLine *x) {
# if defined (DEBUG)
        printf("freeFloatLine: %x %x %d\n",
                (int)x,(int)(x->line),x->tot_nx);
# endif
        if (x->line != NULL)
            free (x->line);
        initFloatLine (x);
}

void initShortLine (ShortHdrLine *x) {
        x->line   = NULL;
        x->tot_nx = 0;
}

int allocShortLine (ShortHdrLine *x, int i) {
# if defined (DEBUG)
        printf("allocShortLine-1: %x %d\n",
                (int)(x->line),x->tot_nx);
# endif
        if (x->line == NULL || x->tot_nx != i) {
            if (x->line != NULL)
                free (x->line);
            x->tot_nx = i;
            x->line = (short *) calloc (x->tot_nx, sizeof(short));
            if (x->line == NULL) {
                initShortLine (x);
                error (NOMEM,"Allocating short line");
                return (-1);
            }
        }
# if defined (DEBUG)
        printf("allocShortLine-2: %x %d\n",
                (int)(x->line),x->tot_nx);
# endif
        return (0);
}

void freeShortLine (ShortHdrLine *x) {
# if defined (DEBUG)
        printf("freeShortLine: %x %x %d\n",
                (int)x,(int)(x->line),x->tot_nx);
# endif
        if (x->line != NULL)
            free (x->line);
        initShortLine (x);
}

void initHdr(Hdr *h) {
# if defined (DEBUG)
        printf("initHdr: %x %d %d %x\n",
                (int)h,h->nlines,h->nalloc,(int)(h->array));
# endif
        h->nlines = 0;
        h->nalloc = 0;
        h->array = NULL;
}

int allocHdr(Hdr *h, int n) {
# if defined (DEBUG)
        printf("allocHdr-1: %x %d %d %x %d\n",
                (int)h,h->nlines,h->nalloc,(int)(h->array),n);
# endif
        h->nlines = 0;
        if (h->array == NULL || h->nalloc != n) {
            if (h->array != NULL)
                free(h->array);
            h->nalloc = n;
            h->array = (HdrArray *)calloc(n,sizeof(HdrArray));
            if (h->array == NULL) {
                h->nalloc = 0;
                error(NOMEM,"Allocating Hdr");
                return -1;
            }
        }
# if defined (DEBUG)
        printf("allocHdr-2: %x %d %d %x %d\n",
                (int)h,h->nlines,h->nalloc,(int)(h->array),n);
# endif
        return 0;
}

int reallocHdr(Hdr *h, int n) {
        int i;
        HdrArray *tmp;
# if defined (DEBUG)
        printf("reallocHdr-1: %x %d %d %x %d\n",
                (int)h,h->nlines,h->nalloc,(int)(h->array),n);
# endif
        if (h->array == NULL || n <= h->nalloc) return -1;
        tmp = (HdrArray *)calloc(n,sizeof(HdrArray));
        if (tmp == NULL) {
            error(NOMEM,"Re-llocating Hdr");
            return -1;
        }
        h->nalloc = n;
        for (i = 0; i < h->nlines; ++i)
            strcpy(tmp[i],h->array[i]);
        free(h->array);
        h->array = tmp;
# if defined (DEBUG)
        printf("reallocHdr-2: %x %d %d %x %d\n",
                (int)h,h->nlines,h->nalloc,(int)(h->array),n);
# endif
        return 0;
}

void freeHdr(Hdr *h) {
# if defined (DEBUG)
        printf("freeHdr: %x %d %d %x\n",
                (int)h,h->nlines,h->nalloc,(int)(h->array));
# endif
        if (h != NULL && h->array != NULL)
            free(h->array);
        if (h != NULL)
            initHdr(h);
}

int copyHdr(Hdr *to, Hdr *from) {
        int i;
        if (allocHdr(to,from->nalloc)) return -1;
        for (i = 0; i < from->nlines; ++i)
            strcpy(to->array[i],from->array[i]);
        to->nlines = from->nlines;
        return 0;
}

/*
** The above are the basic cases, now for the composite cases.
*/
void initFloatHdrData(FloatHdrData *x) {
        x->iodesc = NULL;
        x->section.x_beg = 0;
        x->section.y_beg = 0;
        x->section.sx = 0;
        x->section.sy = 0;
        initHdr(&(x->hdr));
        initFloatData(&(x->data));
}

int allocFloatHdrData(FloatHdrData *x, int i, int j) {
        if (allocFloatData(&(x->data),i,j)) return -1;
        if (allocHdr(&(x->hdr),HdrUnit)) return -1;
        x->section.x_beg = 0;
        x->section.y_beg = 0;
        x->section.sx = i;
        x->section.sy = j;
        return 0;
}

void freeFloatHdrData(FloatHdrData *x) {
        freeFloatData(&(x->data));
        freeHdr(&(x->hdr));
        initFloatHdrData(x);
}

void initShortHdrData(ShortHdrData *x) {
        x->iodesc = NULL;
        x->section.x_beg = 0;
        x->section.y_beg = 0;
        x->section.sx = 0;
        x->section.sy = 0;
        initHdr(&(x->hdr));
        initShortData(&(x->data));
}

int allocShortHdrData(ShortHdrData *x, int i, int j) {
        if (allocShortData(&(x->data),i,j)) return -1;
        if (allocHdr(&(x->hdr),HdrUnit)) return -1;
        x->section.x_beg = 0;
        x->section.y_beg = 0;
        x->section.sx = i;
        x->section.sy = j;
        return 0;
}

void freeShortHdrData(ShortHdrData *x) {
        freeShortData(&(x->data));
        freeHdr(&(x->hdr));
        initShortHdrData(x);
}

void initFloatHdrLine (FloatHdrLine *x) {
        x->iodesc = NULL;
        initHdr (&(x->hdr));
        x->ehdr_loaded = False;
        x->tot_nx = 0;
        x->line   = NULL;
}

int allocFloatHdrLine (FloatHdrLine *x, int i) {
        if (allocFloatLine (x, i)) return (-1);
        if (allocHdr (&(x->hdr),HdrUnit)) return (-1);
        return (0);
}

void freeFloatHdrLine (FloatHdrLine *x) {
        if (x->line != NULL)
            free(x->line);
        freeHdr (&(x->hdr));
        initFloatHdrLine (x);
}

void initShortHdrLine (ShortHdrLine *x) {
        x->iodesc = NULL;
        initHdr (&(x->hdr));
        x->ehdr_loaded = False;
        x->tot_nx = 0;
        x->line   = NULL;
}

int allocShortHdrLine (ShortHdrLine *x, int i) {
        if (allocShortLine (x, i)) return (-1);
        if (allocHdr (&(x->hdr),HdrUnit)) return (-1);
        return (0);
}

void freeShortHdrLine (ShortHdrLine *x) {
        if (x->line != NULL)
            free(x->line);
        freeHdr (&(x->hdr));
        initShortHdrLine (x);
}

void initSingleGroup(SingleGroup *x) {
        x->filename = NULL;
        x->group_num = 0;
        x->globalhdr = NULL;
        initFloatHdrData(&(x->sci));
        initShortHdrData(&(x->dq));
        initFloatHdrData(&(x->err));
}

int allocSingleGroup(SingleGroup *x, int i, int j) {
        if (x->globalhdr == NULL) {
            x->globalhdr = (Hdr *)calloc(1,sizeof(Hdr));
            if (x->globalhdr == NULL) return -1;
            initHdr(x->globalhdr);
        }
        if (allocFloatHdrData(&(x->sci),i,j)) return -1;
        if (allocShortHdrData(&(x->dq),i,j)) return -1;
        if (allocFloatHdrData(&(x->err),i,j)) return -1;
        return 0;
}

void freeSingleGroup(SingleGroup *x) {
        freeFloatHdrData(&(x->err));
        freeShortHdrData(&(x->dq));
        freeFloatHdrData(&(x->sci));
        freeHdr(x->globalhdr);
        if (x->globalhdr != NULL)
                free(x->globalhdr);
        if (x->filename != NULL)
                free(x->filename);
        initSingleGroup(x);
}

void initMultiGroup(MultiGroup *x) {
        x->ngroups = 0;
        x->group = NULL;
}

int allocMultiGroup(MultiGroup *x, int n) {
        int i;
        if (x->group != NULL)
                freeMultiGroup(x);
        x->ngroups = n;
        x->group = (SingleGroup *)calloc(n,sizeof(SingleGroup));
        if (x->group == NULL) {
            x->ngroups = 0;
            error(NOMEM,"Allocating MultiGroup");
            return -1;
        }
        for (i = 0; i < x->ngroups; ++i)
            initSingleGroup(&(x->group[i]));
        x->group[0].globalhdr = (Hdr *)calloc(1,sizeof(Hdr));
        if (x->group[0].globalhdr == NULL) return -1;
        initHdr(x->group[0].globalhdr);
        for (i = 1; i < x->ngroups; ++i)
            x->group[i].globalhdr = x->group[0].globalhdr;
        return 0;
}

void freeMultiGroup(MultiGroup *x) {
        int i;
        if (x->group != NULL) {
            freeSingleGroup(&(x->group[0]));
            for (i = 1; i < x->ngroups; ++i) {
                x->group[i].globalhdr = NULL;
                x->group[i].filename = NULL;
                freeSingleGroup(&(x->group[i]));
            }
        }
        initMultiGroup(x);
}

void initSingleNicmosGroup(SingleNicmosGroup *x) {
        x->filename = NULL;
        x->group_num = 0;
        x->globalhdr = NULL;
        initFloatHdrData(&(x->sci));
        initFloatHdrData(&(x->err));
        initShortHdrData(&(x->dq));
        initShortHdrData(&(x->smpl));
        initFloatHdrData(&(x->intg));
}

int allocSingleNicmosGroup(SingleNicmosGroup *x, int i, int j) {
        if (x->globalhdr == NULL) {
            x->globalhdr = (Hdr *)calloc(1,sizeof(Hdr));
            if (x->globalhdr == NULL) return -1;
            initHdr(x->globalhdr);
        }
        if (allocFloatHdrData(&(x->sci),i,j)) return -1;
        if (allocFloatHdrData(&(x->err),i,j)) return -1;
        if (allocShortHdrData(&(x->dq),i,j)) return -1;
        if (allocShortHdrData(&(x->smpl),i,j)) return -1;
        if (allocFloatHdrData(&(x->intg),i,j)) return -1;
        return 0;
}

void freeSingleNicmosGroup(SingleNicmosGroup *x) {
        freeFloatHdrData(&(x->intg));
        freeShortHdrData(&(x->smpl));
        freeShortHdrData(&(x->dq));
        freeFloatHdrData(&(x->err));
        freeFloatHdrData(&(x->sci));
        freeHdr(x->globalhdr);
        if (x->globalhdr != NULL)
                free(x->globalhdr);
        if (x->filename != NULL)
                free(x->filename);
        initSingleNicmosGroup(x);
}

void initMultiNicmosGroup(MultiNicmosGroup *x) {
        x->ngroups = 0;
        x->group = NULL;
}

int allocMultiNicmosGroup(MultiNicmosGroup *x, int n) {
        int i;
        if (x->group != NULL)
                freeMultiNicmosGroup(x);
        x->ngroups = n;
        x->group = (SingleNicmosGroup *)calloc(n,sizeof(SingleNicmosGroup));
        if (x->group == NULL) {
            x->ngroups = 0;
            error(NOMEM,"Allocating MultiNicmosGroup");
            return -1;
        }
        for (i = 0; i < x->ngroups; ++i)
            initSingleNicmosGroup(&(x->group[i]));
        x->group[0].globalhdr = (Hdr *)calloc(1,sizeof(Hdr));
        if (x->group[0].globalhdr == NULL) return -1;
        initHdr(x->group[0].globalhdr);
        for (i = 1; i < x->ngroups; ++i)
            x->group[i].globalhdr = x->group[0].globalhdr;
        return 0;
}

void freeMultiNicmosGroup(MultiNicmosGroup *x) {
        int i;
        if (x->group != NULL) {
            freeSingleNicmosGroup(&(x->group[0]));
            for (i = 1; i < x->ngroups; ++i) {
                x->group[i].globalhdr = NULL;
                x->group[i].filename = NULL;
                freeSingleNicmosGroup(&(x->group[i]));
            }
        }
        initMultiNicmosGroup(x);
}

void initSingleGroupLine (SingleGroupLine *x) {
        x->filename    = NULL;
        x->group_num   = 0;
        x->line_num    = -1;
        x->phdr_loaded = False;
        x->globalhdr   = NULL;
        initFloatHdrLine (&(x->sci));
        initFloatHdrLine (&(x->err));
        initShortHdrLine (&(x->dq));
}

int allocSingleGroupLine (SingleGroupLine *x, int i) {
        if (x->globalhdr == NULL) {
            x->globalhdr = (Hdr *) calloc (1,sizeof(Hdr));
            if (x->globalhdr == NULL) return (-1);
            initHdr(x->globalhdr);
        }
        if (allocFloatHdrLine (&(x->sci),i)) return (-1);
        if (allocFloatHdrLine (&(x->err),i)) return (-1);
        if (allocShortHdrLine (&(x->dq),i))  return (-1);
        return (0);
}

void freeSingleGroupLine (SingleGroupLine *x) {
        freeFloatHdrLine (&(x->sci));
        freeFloatHdrLine (&(x->err));
        freeShortHdrLine (&(x->dq));
        freeHdr (x->globalhdr);
        if (x->globalhdr != NULL)
                free (x->globalhdr);
        if (x->filename != NULL)
                free (x->filename);
        initSingleGroupLine (x);
}

/*                                                                      **
** Allocate space for the lines of data from each extension of a Single **
** Group                                                                */
int allocSciLine (SingleGroupLine *x) {
    IODesc *xio;
    xio = (IODesc *)(x->sci.iodesc);
    if (allocFloatLine (&(x->sci), xio->dims[0])) return (-1);
    return (0);
}

int allocErrLine (SingleGroupLine *x) {
    IODesc *xio;
    xio = (IODesc *)(x->err.iodesc);
    if (allocFloatLine (&(x->err), xio->dims[0])) return (-1);
    return (0);
}

int allocDQLine (SingleGroupLine *x) {
    IODesc *xio;
    xio = (IODesc *)(x->dq.iodesc);
    if (allocShortLine (&(x->dq), xio->dims[0])) return (-1);
    return (0);
}

/*
** Section 4.
** Low-level I/O functions.
*/

char *getFilename(IODescPtr p) { return ((IODesc *)p)->filename; }
char *getExtname(IODescPtr p) { return ((IODesc *)p)->extname; }
int getExtver(IODescPtr p) { return ((IODesc *)p)->extver; }
int getNaxis1(IODescPtr p) { return ((IODesc *)p)->dims[0]; }
int getNaxis2(IODescPtr p) { return ((IODesc *)p)->dims[1]; }
int getType(IODescPtr p) { return ((IODesc *)p)->type; }

# include "c_iraf.h"

# include "hstioirf.c"

/*
** Section 5.
** High-level I/O Functions.
*/

/*
**
** The function ckNewFile() checks to see whether a file exists or not.  It
** takes a single argument, the full path name of the file.  If there is a
** system environment variable "imclobber" that is set to "yes", the
** function then takes action to remove that file.  In the case of VMS, if
** "imclobber" is set to "yes", the function removes all previous versions
** of the file.
**
** The function returns the following integer values:
**      -1 = The file did exist but was removed.
**       0 = The file does not exist.
**       1 = The file exists and was not removed.
**       2 = The file exists and an error occured attempting to remove it.
**
** The function is used in the following manner.
**      if (ckNewFile("filename.ext") > 0) {
**              error("File already exists");
**      }
**
** The ckNewFile() function runs under both UNIX and VMS.
*/
int ckNewFile(char *fname) {
        char *value;
        FILE *x = fopen(fname,"r");
        if (x == NULL)
            return 0; /* file does not exist */
        /* file exists */
        fclose(x);
        value = getenv("imclobber");
        if (value == NULL)
            return 1; /* file exists and was not removed */
        if ((strcmp(value,"yes") != 0) && (strcmp(value,"YES") != 0))
            return 1; /* file exists and was not removed */
        /* file exists and imclobber is yes */
        if (remove(fname) != 0)
            return 2;
        while (remove(fname) == 0); /* The while loop is for VMS */
        return -1;
}

int openFitsFile(char *filename, unsigned int option) {
        return 0;
}

int closeFitsFile(char *filename) {
        return 0;
}

/*
** Routine to open the input file, read in the primary header information, *
** acquire file pointers to the SingleGroup extensions, read the headers   *
** of the extensions, and allocate space of the appropriate length         *
** for the respective lines of data.  Access to the SingleGroup extensions *
** remains open.                                                           *
**                                                                         */
int openSingleGroupLine (char *fname, int ever, SingleGroupLine *x) {
        IODescPtr in;
        in = openInputImage(fname,"",0); if (hstio_err()) return (-1);
        if (x->globalhdr != NULL)
            free(x->globalhdr);
        if (x->filename != NULL)
            free(x->filename);
        x->filename = (char *) calloc ((strlen(fname) + 1),sizeof(char));
        strcpy (x->filename,fname);
        x->globalhdr = (Hdr *)calloc(1,sizeof(Hdr));
        if (x->globalhdr == NULL) return -1;
        initHdr(x->globalhdr);
        getHeader (in,x->globalhdr); if (hstio_err()) return (-1);
        x->phdr_loaded = True;
        closeImage (in);
        x->group_num = ever;

        /* obtain the file pointers to the individual SingleGroup     *
         * extensions, read the headers, and allocate the proper size *
         * storage for the line arrays.                               */
        getSciHdr (fname,ever,&(x->sci)); if (hstio_err()) return (-1);
        x->sci.ehdr_loaded = True;
        getErrHdr (fname,ever,&(x->err)); if (hstio_err()) return (-1);
        x->err.ehdr_loaded = True;
        getDQHdr  (fname,ever,&(x->dq)); if (hstio_err()) return (-1);
        x->dq.ehdr_loaded = True;
        allocSciLine (x);
        allocErrLine (x);
        allocDQLine  (x);
        clear_err();
        return (0);
}

void closeSingleGroupLine (SingleGroupLine *x) {
        closeImage (x->sci.iodesc);
        closeImage (x->err.iodesc);
        closeImage (x->dq.iodesc);
}

int getFloatHD(char *fname, char *ename, int ever, FloatHdrData *x) {
        IODesc *xio;
        x->iodesc = openInputImage(fname,ename,ever);
        xio = (IODesc *)(x->iodesc);
        if (hstio_err()) return -1;
        x->section.sx = xio->dims[0];
        x->section.sy = xio->dims[1];
        getHeader(x->iodesc,&(x->hdr));
        if (hstio_err()) return -1;
        getFloatData(x->iodesc,&(x->data));
        if (hstio_err()) return -1;
        closeImage(x->iodesc);
        clear_err();
        return 0;
}

int putFloatHD(char *fname, char *ename, int ever, FloatHdrData *x, int option) {
        if (option == 0)
            x->iodesc = openOutputImage(fname, ename, ever, &(x->hdr),
                x->data.tot_nx, x->data.tot_ny, FITSFLOAT);
        else if (option & Overwrite)
            x->iodesc = openUpdateImage(fname, ename, ever, &(x->hdr));
        if (hstio_err()) return -1;
        putFloatData(x->iodesc,&(x->data));
        if (hstio_err()) return -1;
        closeImage(x->iodesc);
        clear_err();
        return 0;
}

int getShortHD(char *fname, char *ename, int ever, ShortHdrData *x) {
        IODesc *xio;
        x->iodesc = openInputImage(fname,ename,ever);
        xio = (IODesc *)(x->iodesc);
        if (hstio_err()) return -1;
        x->section.sx = xio->dims[0];
        x->section.sy = xio->dims[1];
        getHeader(x->iodesc,&(x->hdr));
        if (hstio_err()) return -1;
        getShortData(x->iodesc,&(x->data));
        if (hstio_err()) return -1;
        closeImage(x->iodesc);
        clear_err();
        return 0;
}

int putShortHD(char *fname, char *ename, int ever, ShortHdrData *x, int option) {
        if (option == 0)
            x->iodesc = openOutputImage(fname, ename, ever, &(x->hdr),
                x->data.tot_nx, x->data.tot_ny, FITSSHORT);
        else if (option & Overwrite)
            x->iodesc = openUpdateImage(fname, ename, ever, &(x->hdr));
        if (hstio_err()) return -1;
        putShortData(x->iodesc,&(x->data));
        if (hstio_err()) return -1;
        closeImage(x->iodesc);
        clear_err();
        return 0;
}

/* Routine to support the routines which write out a subsection of data in *
 * memory to output files.  XBEG and YBEG are zero-indexed values.  The    *
 * coordinate values in the headers are one-indexed.                       */
void updateWCS (Hdr *hdr, int xbeg, int ybeg) {
        FitsKw kw;
        float  old_LTV, new_LTV;
        double old_CRPIX, new_CRPIX;

        kw = findKw(hdr,"LTV1");
        if (kw != 0) {
            old_LTV = getFloatKw (kw);
            new_LTV = old_LTV - (float)xbeg;
            putFloatKw (kw, new_LTV);
        }
        kw = findKw(hdr,"LTV2");
        if (kw != 0) {
            old_LTV = getFloatKw (kw);
            new_LTV = old_LTV - (float)ybeg;
            putFloatKw (kw, new_LTV);
        }

        kw = findKw(hdr,"CRPIX1");
        if (kw != 0) {
            old_CRPIX = getDoubleKw (kw);
            new_CRPIX = old_CRPIX - (double)xbeg;
            putDoubleKw (kw, new_CRPIX);
        }
        kw = findKw(hdr,"CRPIX2");
        if (kw != 0) {
            old_CRPIX = getDoubleKw (kw);
            new_CRPIX = old_CRPIX - (double)ybeg;
            putDoubleKw (kw, new_CRPIX);
        }
}

int putFloatHDSect(char *fname, char *ename, int ever, FloatHdrData *x, int xbeg, int ybeg, int xsize, int ysize, int option) {

        /* Update the LTV keywords */
        updateWCS (&(x->hdr), xbeg, ybeg);

        if (option == 0)
            x->iodesc = openOutputImage(fname, ename, ever, &(x->hdr),
                xsize, ysize, FITSFLOAT);
        else if (option & Overwrite)
            x->iodesc = openUpdateImage(fname, ename, ever, &(x->hdr));
        if (hstio_err()) return -1;
        putFloatSect(x->iodesc,&(x->data),xbeg,ybeg,xsize,ysize);
        if (hstio_err()) return -1;
        closeImage(x->iodesc);
        clear_err();
        return 0;
}

int putShortHDSect(char *fname, char *ename, int ever, ShortHdrData *x, int xbeg, int ybeg, int xsize, int ysize, int option) {

        /* Update the LTV keywords */
        updateWCS (&(x->hdr), xbeg, ybeg);

        if (option == 0)
            x->iodesc = openOutputImage(fname, ename, ever, &(x->hdr),
                xsize, ysize, FITSSHORT);
        else if (option & Overwrite)
            x->iodesc = openUpdateImage(fname, ename, ever, &(x->hdr));
        if (hstio_err()) return -1;
        putShortSect(x->iodesc,&(x->data),xbeg,ybeg,xsize,ysize);
        if (hstio_err()) return -1;
        closeImage(x->iodesc);
        clear_err();
        return 0;
}

int getFloatHdr (char *fname, char *ename, int ever, FloatHdrLine *x) {
        IODesc *xio;
        FitsKw kw;
        int dim1, dim2, no_dims;
        int status = 0;

        x->iodesc = openInputImage (fname,ename,ever);
        if (hstio_err()) return (-1);
        getHeader (x->iodesc,&(x->hdr));

        /* determine dimensions for images which contain a constant value */
        xio = (IODesc *)(x->iodesc);
        if (fits_get_img_dim(xio->ff, &no_dims, &status)) {
            ioerr(BADDIMS, xio, status);
            return -1;
        }
        if (no_dims == 0) {
            kw = findKw(xio->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADDIMS,xio,0); return -1; }
            dim1 = getIntKw(kw);
            kw = findKw(xio->hdr,"NPIX2");
            if (kw == 0) { ioerr(BADDIMS,xio,0); return -1; }
            dim2 = getIntKw(kw);
            xio->dims[0] = dim1;
            xio->dims[1] = dim2;
        }

        if (hstio_err()) return (-1);
        clear_err();
        return (0);
}

int getShortHdr (char *fname, char *ename, int ever, ShortHdrLine *x) {
        IODesc *xio;
        FitsKw kw;
        int dim1, dim2, no_dims;
        int status = 0;

        x->iodesc = openInputImage (fname,ename,ever);
        if (hstio_err()) return (-1);
        getHeader (x->iodesc,&(x->hdr));

        /* determine dimensions for images which contain a constant value */
        xio = (IODesc *)(x->iodesc);
        if (fits_get_img_dim(xio->ff, &no_dims, &status)) {
            ioerr(BADDIMS, xio, status);
            return -1;
        }
        if (no_dims == 0) {
            kw = findKw(xio->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADDIMS,xio,0); return -1; }
            dim1 = getIntKw(kw);
            kw = findKw(xio->hdr,"NPIX2");
            if (kw == 0) { ioerr(BADDIMS,xio,0); return -1; }
            dim2 = getIntKw(kw);
            xio->dims[0] = dim1;
            xio->dims[1] = dim2;
        }

        if (hstio_err()) return (-1);
        clear_err();
        return (0);
}

int getSci(char *fname, int ever, SciHdrData *x) {
        return getFloatHD(fname,"SCI",ever,x); }
int putSci(char *fname, int ever, SciHdrData *x, int option) {
        return putFloatHD(fname,"SCI",ever,x,option); }
int getErr(char *fname, int ever, ErrHdrData *x) {
        return getFloatHD(fname,"ERR",ever,x); }
int putErr(char *fname, int ever, ErrHdrData *x, int option) {
        return putFloatHD(fname,"ERR",ever,x,option); }
int getDQ(char *fname, int ever, DQHdrData *x) {
        return getShortHD(fname,"DQ",ever,x); }
int putDQ(char *fname, int ever, DQHdrData *x, int option) {
        return putShortHD(fname,"DQ",ever,x,option); }
int getSmpl(char *fname, int ever, SmplHdrData *x) {
        return getShortHD(fname,"SAMP",ever,x); }
int putSmpl(char *fname, int ever, SmplHdrData *x, int option) {
        return putShortHD(fname,"SAMP",ever,x,option); }
int getIntg(char *fname, int ever, IntgHdrData *x) {
        return getFloatHD(fname,"TIME",ever,x); }
int putIntg(char *fname, int ever, IntgHdrData *x, int option) {
        return putFloatHD(fname,"TIME",ever,x,option); }

/*                                                                      **
** Routines to output a subsection of an image in memory to a disk file **
** where the subsection is the full size (NAXIS1/NAXIS2) of the output  **
** image.                                                               **
**                                                                      */
int putSciSect(char *fname, int ever, SciHdrData *x, int xbeg, int ybeg,
               int xsize, int ysize, int option) {
    return (putFloatHDSect(fname,"SCI",ever,x,xbeg,ybeg,xsize,ysize,option)); }

int putErrSect(char *fname, int ever, ErrHdrData *x, int xbeg, int ybeg,
               int xsize, int ysize, int option) {
    return (putFloatHDSect(fname,"ERR",ever,x,xbeg,ybeg,xsize,ysize,option)); }

int putDQSect(char *fname, int ever, DQHdrData *x, int xbeg, int ybeg,
               int xsize, int ysize, int option) {
    return (putShortHDSect(fname,"DQ",ever,x,xbeg,ybeg,xsize,ysize,option)); }

int putSmplSect(char *fname, int ever, SmplHdrData *x, int xbeg, int ybeg,
               int xsize, int ysize, int option) {
    return (putShortHDSect(fname,"SAMP",ever,x,xbeg,ybeg,xsize,ysize,option)); }

int putIntgSect(char *fname, int ever, IntgHdrData *x, int xbeg, int ybeg,
               int xsize, int ysize, int option) {
    return (putFloatHDSect(fname,"TIME",ever,x,xbeg,ybeg,xsize,ysize,option)); }

/* Get just the header for the extension */
int getSciHdr (char *fname, int ever, SciHdrLine *x) {
        return getFloatHdr (fname,"SCI",ever,x); }
int getErrHdr (char *fname, int ever, ErrHdrLine *x) {
        return getFloatHdr (fname,"ERR",ever,x); }
int getDQHdr(char *fname, int ever, DQHdrLine *x) {
        return getShortHdr (fname,"DQ",ever,x); }

/* Get just the data line for the extension */
int getSciLine (SciHdrLine *x, int line_num) {
        return (getFloatLine (x->iodesc, line_num, x->line));
}
int getErrLine (ErrHdrLine *x, int line_num) {
        return (getFloatLine (x->iodesc, line_num, x->line));
}
int getDQLine (DQHdrLine *x, int line_num) {
        return (getShortLine (x->iodesc, line_num, x->line));
}

int getSingleGroup(char *fname, int ever, SingleGroup *x) {
        IODescPtr in;
        in = openInputImage(fname,"",0); if (hstio_err()) return -1;
        if (x->globalhdr != NULL)
            free(x->globalhdr);
        if (x->filename != NULL)
            free(x->filename);
        x->filename = (char *)calloc((strlen(fname) + 1),sizeof(char));
        strcpy(x->filename,fname);
        x->globalhdr = (Hdr *)calloc(1,sizeof(Hdr));
        if (x->globalhdr == NULL) return -1;
        initHdr(x->globalhdr);
        getHeader(in,x->globalhdr); if (hstio_err()) return -1;
        closeImage(in);
        x->group_num = ever;
        getSci(fname,ever,&(x->sci)); if (hstio_err()) return -1;
        getErr(fname,ever,&(x->err)); if (hstio_err()) return -1;
        getDQ(fname,ever,&(x->dq)); if (hstio_err()) return -1;
        clear_err();
        return 0;
}

int getSingleGroupLine (char *fname, int line, SingleGroupLine  *x) {
        x->line_num = line;
        getSciLine(&(x->sci), line);
        if (hstio_err()) return (-1);
        getErrLine(&(x->err), line);
        if (hstio_err()) return (-1);
        getDQLine(&(x->dq), line);
        if (hstio_err()) return (-1);
        clear_err();
        return (0);
}

int putSingleGroupHdr(char *fname, SingleGroup *x, int option) {
        IODescPtr out = NULL;
        if (option == 0)
            out = openOutputImage(fname,"",0,x->globalhdr,0,0,FITSBYTE);
        else if (option & Overwrite)
            out = openUpdateImage(fname,"",0,x->globalhdr);
        if (hstio_err()) return -1;
        closeImage(out);
        clear_err();
        return 0;
}

int putSingleGroup(char *fname, int ever, SingleGroup *x, int option) {
        struct stat buf;
        if (option == 0) {
            if (stat(fname,&buf) == -1)
                putSingleGroupHdr(fname,x,0);
        }
        putSci(fname,ever,&(x->sci),option); if (hstio_err()) return -1;
        putErr(fname,ever,&(x->err),option); if (hstio_err()) return -1;
        putDQ(fname,ever,&(x->dq),option); if (hstio_err()) return -1;
        clear_err();
        return 0;
}

/*                                                                           **
** Routine to output a subsection of an imset in memory to a disk file where **
** the subsection is the full size (NAXIS1/NAXIS2) of the output image.      **
**                                                                           */
int putSingleGroupSect(char *fname, int ever, SingleGroup *x, int xbeg,
    int ybeg, int xsize, int ysize, int option) {
        struct stat buf;

        if (option == 0) {
            if (stat(fname,&buf) == -1)
                putSingleGroupHdr(fname,x,0);
        }

        putSciSect(fname,ever,&(x->sci),xbeg,ybeg,xsize,ysize,option);
        if (hstio_err()) return -1;

        putErrSect(fname,ever,&(x->err),xbeg,ybeg,xsize,ysize,option);
        if (hstio_err()) return -1;

        putDQSect (fname,ever,&(x->dq),xbeg,ybeg,xsize,ysize,option);
        if (hstio_err()) return -1;

        clear_err ();
        return 0;
}

int getSingleNicmosGroup(char *fname, int ever, SingleNicmosGroup *x) {
        IODescPtr in;
        in = openInputImage(fname,"",0); if (hstio_err()) return -1;
        if (x->globalhdr != NULL)
            free(x->globalhdr);
        if (x->filename != NULL)
            free(x->filename);
        x->filename = (char *)calloc((strlen(fname) + 1),sizeof(char));
        strcpy(x->filename,fname);
        x->globalhdr = (Hdr *)calloc(1,sizeof(Hdr));
        if (x->globalhdr == NULL) return -1;
        initHdr(x->globalhdr);
        getHeader(in,x->globalhdr); if (hstio_err()) return -1;
        closeImage(in);
        x->group_num = ever;
        getSci(fname,ever,&(x->sci)); if (hstio_err()) return -1;
        getErr(fname,ever,&(x->err)); if (hstio_err()) return -1;
        getDQ(fname,ever,&(x->dq)); if (hstio_err()) return -1;
        getSmpl(fname,ever,&(x->smpl)); if (hstio_err()) return -1;
        getIntg(fname,ever,&(x->intg)); if (hstio_err()) return -1;
        clear_err();
        return 0;
}

int putSingleNicmosGroupHdr(char *fname, SingleNicmosGroup *x, int option) {
        IODescPtr out = NULL;
        if (option == 0)
            out = openOutputImage(fname,"",0,x->globalhdr,0,0,FITSBYTE);
        else if (option & Overwrite)
            out = openUpdateImage(fname,"",0,x->globalhdr);
        if (hstio_err()) return -1;
        closeImage(out);
        clear_err();
        return 0;
}

int putSingleNicmosGroup(char *fname, int ever, SingleNicmosGroup *x,
        int option) {
        struct stat buf;
        if (option == 0) {
            if (stat(fname,&buf) == -1)
                putSingleNicmosGroupHdr(fname,x,0);
        }
        putSci(fname,ever,&(x->sci),option); if (hstio_err()) return -1;
        putErr(fname,ever,&(x->err),option); if (hstio_err()) return -1;
        putDQ(fname,ever,&(x->dq),option); if (hstio_err()) return -1;
        putSmpl(fname,ever,&(x->smpl),option); if (hstio_err()) return -1;
        putIntg(fname,ever,&(x->intg),option); if (hstio_err()) return -1;
        clear_err();
        return 0;
}

/*                                                                           **
** Routine to output a subsection of an imset in memory to a disk file where **
** the subsection is the full size (NAXIS1/NAXIS2) of the output image.      **
**                                                                           */
int putSingleNicmosGroupSect(char *fname, int ever, SingleNicmosGroup *x,
    int xbeg, int ybeg, int xsize, int ysize, int option) {
        struct stat buf;

        if (option == 0) {
            if (stat(fname,&buf) == -1)
                putSingleNicmosGroupHdr(fname,x,0);
        }

        putSciSect(fname,ever,&(x->sci),xbeg,ybeg,xsize,ysize,option);
        if (hstio_err()) return -1;

        putErrSect(fname,ever,&(x->err),xbeg,ybeg,xsize,ysize,option);
        if (hstio_err()) return -1;

        putDQSect (fname,ever,&(x->dq),xbeg,ybeg,xsize,ysize,option);
        if (hstio_err()) return -1;

        putSmplSect(fname,ever,&(x->smpl),xbeg,ybeg,xsize,ysize,option);
        if (hstio_err()) return -1;

        putIntgSect(fname,ever,&(x->intg),xbeg,ybeg,xsize,ysize,option);
        if (hstio_err()) return -1;

        clear_err ();
        return 0;
}

int getMultiGroupHdr(char *fname, MultiGroup *x) {
        IODescPtr in;
        int i;
        in = openInputImage(fname,"",0); if (hstio_err()) return -1;
        getHeader(in,x->group[0].globalhdr); if (hstio_err()) return -1;
        closeImage(in);
        if (x->group[0].filename != NULL)
            free(x->group[0].filename);
        x->group[0].filename = (char *)calloc((strlen(fname) + 1),sizeof(char));
        strcpy(x->group[0].filename,fname);
        for (i = 1; i < x->ngroups; ++i) {
            x->group[i].filename = x->group[0].filename;
            x->group[i].globalhdr = x->group[0].globalhdr;
        }
        clear_err();
        return 0;
}

int getMultiGroup(MultiGroup *x, int ng, int ever) {
        if (ng < 0 || ng > x->ngroups) { error(BADGROUP,""); return -1; }
        x->group[ng].group_num = ever;
        getSci(x->group[ng].filename,ever,&(x->group[ng].sci));
        if (hstio_err()) return -1;
        getErr(x->group[ng].filename,ever,&(x->group[ng].err));
        if (hstio_err()) return -1;
        getDQ(x->group[ng].filename,ever,&(x->group[ng].dq));
        if (hstio_err()) return -1;
        clear_err();
        return 0;
}

int putMultiGroupHdr(char *fname, MultiGroup *x, int option) {
        IODescPtr out = NULL;
        if (option == 0)
            out = openOutputImage(fname,"",0,x->group[0].globalhdr,0,0,FITSBYTE);
        else if (option & Overwrite)
            out = openUpdateImage(fname,"",0,x->group[0].globalhdr);
        if (hstio_err()) return -1;
        putHeader(out); if (hstio_err()) return -1;
        closeImage(out);
        clear_err();
        return 0;
}

int putMultiGroup(char *fname, int ever, MultiGroup *x, int ng, int option) {
        struct stat buf;
        if (ng < 0 || ng > x->ngroups) { error(BADGROUP,""); return -1; }
        if (option == 0) {
            if (stat(fname,&buf) == -1)
                putMultiGroupHdr(fname,x,0);
        }
        putSci(fname,ever,&(x->group[ng].sci),option); if (hstio_err()) return -1;
        putErr(fname,ever,&(x->group[ng].err),option); if (hstio_err()) return -1;
        putDQ(fname,ever,&(x->group[ng].dq),option); if (hstio_err()) return -1;
        clear_err();
        return 0;
}

int getMultiNicmosGroupHdr(char *fname, MultiNicmosGroup *x) {
        IODescPtr in;
        int i;
        in = openInputImage(fname,"",0); if (hstio_err()) return -1;
        getHeader(in,x->group[0].globalhdr); if (hstio_err()) return -1;
        closeImage(in);
        if (x->group[0].filename != NULL)
            free(x->group[0].filename);
        x->group[0].filename = (char *)calloc((strlen(fname) + 1),sizeof(char));
        strcpy(x->group[0].filename,fname);
        for (i = 1; i < x->ngroups; ++i) {
            x->group[i].filename = x->group[0].filename;
            x->group[i].globalhdr = x->group[0].globalhdr;
        }
        clear_err();
        return 0;
}

int getMultiNicmosGroup(MultiNicmosGroup *x, int ng, int ever) {
        if (ng < 0 || ng > x->ngroups) { error(BADGROUP,""); return -1; }
        x->group[ng].group_num = ever;
        getSci(x->group[ng].filename,ever,&(x->group[ng].sci));
        if (hstio_err()) return -1;
        getErr(x->group[ng].filename,ever,&(x->group[ng].err));
        if (hstio_err()) return -1;
        getDQ(x->group[ng].filename,ever,&(x->group[ng].dq));
        if (hstio_err()) return -1;
        getSmpl(x->group[ng].filename,ever,&(x->group[ng].smpl));
        if (hstio_err()) return -1;
        getIntg(x->group[ng].filename,ever,&(x->group[ng].intg));
        if (hstio_err()) return -1;
        clear_err();
        return 0;
}

int putMultiNicmosGroupHdr(char *fname, MultiNicmosGroup *x, int option) {
        IODescPtr out = NULL;
        if (option == 0)
            out = openOutputImage(fname,"",0,x->group[0].globalhdr,0,0,FITSBYTE);
        else if (option & Overwrite)
            out = openUpdateImage(fname,"",0,x->group[0].globalhdr);
        if (hstio_err()) return -1;
        putHeader(out); if (hstio_err()) return -1;
        closeImage(out);
        clear_err();
        return 0;
}

int putMultiNicmosGroup(char *fname, int ever, MultiNicmosGroup *x, int ng,
        int option) {
        struct stat buf;
        if (ng < 0 || ng > x->ngroups) { error(BADGROUP,""); return -1; }
        if (option == 0) {
            if (stat(fname,&buf) == -1)
                putMultiNicmosGroupHdr(fname,x,0);
        }
        putSci(fname,ever,&(x->group[ng].sci),option); if (hstio_err()) return -1;
        putErr(fname,ever,&(x->group[ng].err),option); if (hstio_err()) return -1;
        putDQ(fname,ever,&(x->group[ng].dq),option); if (hstio_err()) return -1;
        putSmpl(fname,ever,&(x->group[ng].smpl),option); if (hstio_err()) return -1;
        putIntg(fname,ever,&(x->group[ng].intg),option); if (hstio_err()) return -1;
        clear_err();
        return 0;
}

/*
** Section 6.
** Functions to manipulate the header array.
**
** See the file keyword.c
*/
