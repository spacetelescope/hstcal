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
** Phil Hodge, May 2011:
** In putHeader, check status after calling fits_update_key or fits_delete_key.
** For fits_update_key, if status is non-zero return -1.  For fits_delete_key,
** however, if status is KEY_NO_EXIST, clear the error messages and set status
** to 0.
**
**
** Sara Ogaz, April 2017
** Update putHeader: fits_update_keyword should update the keyword and
** comment values if they exist or append new values, and it should be
** a NULL pointer to leave the comment field alone. putHeader was supplying
** the empty string directly, this has been changed to NULL.
**
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
# include <fitsio.h>
# include <ctype.h>
# include <stdio.h>
# include <string.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <time.h>
# include <unistd.h>
# include <assert.h>
# include <stdlib.h>

# include "hstio.h"
# include "hstcalerr.h"

/*
** String defined to allow determination of the HSTIO library version
** from the library file (*.a) or the executable using the library.
*/
const char *hstio_version = HSTIO_VERSION;

void initPtrRegister(PtrRegister * reg)
{
    reg->cursor = 0; //points to last ptr NOT next slot
    reg->length = PTR_REGISTER_LENGTH_INC;
    reg->ptrs = calloc(reg->length+1, sizeof(*reg->ptrs));
    assert(reg->ptrs);
    reg->freeFunctions = calloc(reg->length+1, sizeof(*reg->freeFunctions));
    if (!reg->freeFunctions)
    {
        free(reg->ptrs);
        assert(0);
    }
    reg->ptrs[0] = reg; //this ptr
    reg->freeFunctions[0] = &free;
}
void addPtr(PtrRegister * reg, void * ptr, void * freeFunc)
{
    if (!reg || !ptr || !freeFunc)
        return;

    //check ptr isn't already registered? - go on then.
    {int i;
    for (i = reg->cursor; i >= 0 ; --i)// i >= 0 prevents adding self again
    {
        if (reg->ptrs[i] == ptr)
            return;
    }}

    if (++reg->cursor >= reg->length)
    {
        reg->length += PTR_REGISTER_LENGTH_INC;
        assert(reg->ptrs = realloc(reg->ptrs, reg->length*sizeof(*reg->ptrs)));
        assert(reg->freeFunctions = realloc(reg->freeFunctions, reg->length*sizeof(*reg->freeFunctions)));
    }
    reg->ptrs[reg->cursor] = ptr;
    reg->freeFunctions[reg->cursor] = freeFunc;
}
void freePtr(PtrRegister * reg, void * ptr)
{
    //Can't be used to free itself, use freeReg(), use of i > 0 in below for is reason.
    if (!reg || !ptr)
        return;

    int i;
    for (i = reg->cursor; i > 0 ; --i)
    {
        if (reg->ptrs[i] == ptr)
            break;
    }

    //call function to free ptr
    reg->freeFunctions[i](ptr);

    if (i == reg->cursor)
    {
        reg->ptrs[i] = NULL;
        reg->freeFunctions[i] = NULL;
    }
    else
    {
        //move last one into gap to close - not a stack so who cares
        reg->ptrs[i] = reg->ptrs[reg->cursor];
        reg->ptrs[reg->cursor] = NULL;
        reg->freeFunctions[i] = reg->freeFunctions[reg->cursor];
        reg->freeFunctions[reg->cursor] = NULL;
    }
    --reg->cursor;
}
void freeAll(PtrRegister * reg)
{
    if (!reg || reg->length == 0 || reg->cursor == 0)
        return;

    {unsigned i;
    for (i = 1; i <= reg->cursor; ++i)
    {
        if (reg->freeFunctions[i] && reg->ptrs[i])
        {
            reg->freeFunctions[i](reg->ptrs[i]);
            reg->ptrs[i] = NULL;
            reg->freeFunctions[i] = NULL;
        }
    }}
    reg->cursor = 0;
}
void freeReg(PtrRegister * reg)
{
    /* THIS SHOULD NEVER CALL freeALL()
     * This is designed to be used when allocating multiple persistent memory allocations,
     * registering each allocation in turn. If one allocation fails freeOnExit() can be
     * called to free prior successful allocations and if all allocations are successful
     * this function can be called to free the registers without freeing the actual
     * pointers just allocated.
     */

    if (!reg || reg->length == 0)
        return;

    reg->cursor = 0;
    reg->length = 0;
    // free 'itself'
    free(reg->ptrs);
    reg->ptrs = NULL;
    free(reg->freeFunctions);
    reg->freeFunctions = NULL;
}
void freeOnExit(PtrRegister * reg)
{
    //Free everything registered
    freeAll(reg);
    //Free itself
    freeReg(reg);
}

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

static void ioerr(HSTIOError e, IODescPtr x_, int status) {
        IODesc *x;
        char cfitsio_errmsg[81];
        x = (IODesc *)x_;
        sprintf(&error_msg[strlen(error_msg)],
                "Filename %s EXTNAME %s EXTVER %d CFITSIO status %d\n",
                x->filename, x->extname, x->extver, status);
        while (fits_read_errmsg(cfitsio_errmsg)) {
            strncat(error_msg, cfitsio_errmsg, 80);
        }
        error(e,0);
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
        x->storageOrder = ROWMAJOR;
        x->data = NULL;
# if defined (DEBUG)
        printf("initFloatData: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
}

int allocFloatData(FloatTwoDArray *x, int i, int j, Bool zeroInitialize) {
    //WARNING: target (x) must be initialized by caller
# if defined (DEBUG)
        printf("allocFloatData-1: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
        if (x->buffer == NULL || x->buffer_size != (i * j)) {
            if (x->buffer != NULL)
            {
                free(x->buffer);
                x->buffer = NULL;
            }
            x->buffer_size = i * j;
            if (zeroInitialize)
                x->buffer = calloc(x->buffer_size, sizeof(*x->buffer));
            else
                x->buffer = malloc(x->buffer_size * sizeof(*x->buffer));
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
        if (!x)
            return;
# if defined (DEBUG)
        printf("freeFloatData: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
        if (x->buffer != NULL)
            free(x->buffer);
        initFloatData(x);
}

int copyFloatData(FloatTwoDArray * target, const FloatTwoDArray * source, enum StorageOrder targetStorageOrder)
{
    if (!target || !source)
        return -1;

    //should this check be raise higher up the call stack also?
    if (targetStorageOrder != source->storageOrder)
    {
        //assumes target initialized
        if (!target->buffer)
        {
            if (allocFloatData(target, source->nx, source->ny, False))
                return ALLOCATION_PROBLEM; //allocFloatData() initializes before returning
        }
        return swapFloatStorageOrder(target, source, targetStorageOrder);
        //fall through and copy normally
    }

    if (allocFloatData(target, source->nx, source->ny, False))
        return ALLOCATION_PROBLEM; //allocFloatData() initializes before returning

    //allocFloatData() correctly initializes all other members leaving only buffer (data points to buffer)
    memcpy(target->buffer, source->buffer, source->nx*source->ny*sizeof(*source->buffer));
    return 0;
}

int swapFloatStorageOrder(FloatTwoDArray * target, const FloatTwoDArray * source, enum StorageOrder targetStorageOrder)
{
    //this probably breaks use of Pix on target? Do we need to swap nx & ny?
    if (!target || !source)
        return -1;

    target->storageOrder = targetStorageOrder;
    if (targetStorageOrder == source->storageOrder)
        return 0;

    const unsigned nRows = target->ny;
    const unsigned nCols = target->nx;

    {unsigned j;
    for (j = 0; j < nCols; ++j)
    {
        {unsigned i;
        for (i = 0; i < nRows; ++i)
        {
            if (targetStorageOrder == COLUMNMAJOR)
                target->data[j*nRows + i] = source->data[i*nCols + j];
            else
                target->data[i*nCols + j] = source->data[j*nRows + i];
        }}
    }}
    return 0;
}

int swapShortStorageOrder(ShortTwoDArray * target, const ShortTwoDArray * source, enum StorageOrder targetStorageOrder)
{
    //this probably breaks use of Pix on target? Do we need to swap nx & ny?
    if (!target || !source)
        return -1;

    target->storageOrder = targetStorageOrder;
    if (targetStorageOrder == source->storageOrder)
        return 0;

    const unsigned nRows = target->ny;
    const unsigned nCols = target->nx;

    {unsigned j;
    for (j = 0; j < nCols; ++j)
    {
        {unsigned i;
        for (i = 0; i < nRows; ++i)
        {
            if (targetStorageOrder == COLUMNMAJOR)
                target->data[j*nRows + i] = source->data[i*nCols + j];
            else
                target->data[i*nCols + j] = source->data[j*nRows + i];
        }}
    }}
    return 0;
}

void initShortData(ShortTwoDArray *x) {
        x->buffer = NULL;
        x->buffer_size = 0;
        x->tot_nx = 0;
        x->tot_ny = 0;
        x->nx = 0;
        x->ny = 0;
        x->storageOrder = ROWMAJOR;
        x->data = NULL;
# if defined (DEBUG)
        printf("initShortData: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
}

int allocShortData(ShortTwoDArray *x, int i, int j, Bool zeroInitialize) {
# if defined (DEBUG)
        printf("allocShortData-1: %x %x %d\n",
                (int)x,(int)(x->buffer),x->buffer_size);
# endif
        if (x->buffer == NULL || x->buffer_size != (i * j)) {
            if (x->buffer != NULL)
                free(x->buffer);
            x->buffer_size = i * j;
            if (zeroInitialize)
                x->buffer = calloc(x->buffer_size, sizeof(*x->buffer));
            else
                x->buffer = malloc(x->buffer_size * sizeof(*x->buffer));
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
int copyShortData(ShortTwoDArray * target, const ShortTwoDArray * source, enum StorageOrder targetStorageOrder)
{
    if (!target || !source)
        return -1;

    //should this check be raise higher up the call stack also?
    if (targetStorageOrder != source->storageOrder)
    {
        //assumes target initialized
        if (!target->buffer)
        {
            if (allocShortData(target, source->nx, source->ny, False))
                return ALLOCATION_PROBLEM; //allocShortData() initializes before returning
        }
        return swapShortStorageOrder(target, source, targetStorageOrder);
        //fall through and copy normally
    }

    if (allocShortData(target, source->nx, source->ny, False))
        return ALLOCATION_PROBLEM; //allocShortData() initializes before returning

    //allocShortData() correctly initializes all other members leaving only buffer (data points to buffer)
    memcpy(target->buffer, source->buffer, source->nx*source->ny*sizeof(*source->buffer));
    return 0;
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

int allocHdr(Hdr *h, int n, Bool zeroInitialize) {
# if defined (DEBUG)
        printf("allocHdr-1: %x %d %d %x %d\n",
                (int)h,h->nlines,h->nalloc,(int)(h->array),n);
# endif
        h->nlines = 0;
        if (h->array == NULL || h->nalloc != n) {
            if (h->array != NULL)
                free(h->array);
            h->nalloc = n;
            if (zeroInitialize)
                h->array = calloc(n,sizeof(*h->array));
            else
                h->array = malloc(n * sizeof(*h->array));

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
        if (!h)
        	return;
        if (h->array)
            free(h->array);
        initHdr(h);
}

int copyHdr(Hdr *to, const Hdr *from) {
        if (!to || !from)
            return -1;
        //allcoHdr only allocates if to->array == NULL or sizes differ
        if (allocHdr(to,from->nalloc, False)) return -1;
        memcpy(to->array, from->array, to->nalloc*sizeof(*to->array));
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

int allocFloatHdrData(FloatHdrData *x, int i, int j, Bool zeroInitialize) {
        if (allocFloatData(&(x->data),i,j, zeroInitialize)) return -1;
        if (allocHdr(&(x->hdr),HdrUnit, zeroInitialize)) return -1;
        x->section.x_beg = 0;
        x->section.y_beg = 0;
        x->section.sx = i;
        x->section.sy = j;
        return 0;
}

int copyFloatHdrData(FloatHdrData * target, const FloatHdrData * src, enum StorageOrder targetStorageOrder)
{
    if (!target || !src)
        return -1;

    target->iodesc = src->iodesc;

    //Since DataSection section refers to image IO, keep as source (I think?).
    copyDataSection(&target->section, &src->section);//No allocations

    if (copyHdr(&target->hdr, &src->hdr))//This allocates
        return ALLOCATION_PROBLEM;

    return copyFloatData(&target->data, &src->data, targetStorageOrder);
}

void freeFloatHdrData(FloatHdrData *x) {
        freeFloatData(&(x->data));
        freeHdr(&(x->hdr));
        initFloatHdrData(x);
}

void copyDataSection(DataSection * dest, const DataSection * src)
{
    dest->x_beg = src->x_beg;
    dest->y_beg = src->y_beg;
    dest->sx = src->sx;
    dest->sy = src->sy;
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

int allocShortHdrData(ShortHdrData *x, int i, int j, Bool zeroInitialize) {
        if (allocShortData(&(x->data),i,j, zeroInitialize)) return -1;
        if (allocHdr(&(x->hdr),HdrUnit, zeroInitialize)) return -1;
        x->section.x_beg = 0;
        x->section.y_beg = 0;
        x->section.sx = i;
        x->section.sy = j;
        return 0;
}

int copyShortHdrData(ShortHdrData * target, const ShortHdrData * src, enum StorageOrder targetStorageOrder)
{
    if (!target || !src)
        return -1;

    target->iodesc = src->iodesc;

    //Since DataSection section refers to image IO, keep as source (I think?).
    copyDataSection(&target->section, &src->section);//No allocations

    if (copyHdr(&target->hdr, &src->hdr))//This allocates
        return ALLOCATION_PROBLEM;

    return copyShortData(&target->data, &src->data, targetStorageOrder);
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
        if (allocHdr (&(x->hdr),HdrUnit, True)) return (-1);
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
        if (allocHdr (&(x->hdr),HdrUnit, True)) return (-1);
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

int allocSingleGroup(SingleGroup *x, int i, int j, Bool zeroInitialize)
{
    if (allocSingleGroupHeader(&x->globalhdr, zeroInitialize) ||
            allocFloatHdrData(&(x->sci),i,j, zeroInitialize)  ||
            allocShortHdrData(&(x->dq),i,j, zeroInitialize)   ||
            allocFloatHdrData(&(x->err),i,j, zeroInitialize))
        return ALLOCATION_PROBLEM;
    return 0;
}

int allocSingleGroupHeader(Hdr ** hdr, Bool zeroInitialize)
{
    if (!hdr)
        return ALLOCATION_PROBLEM;

    if (*hdr)
        return 0; //Already allocated

    if (zeroInitialize)
        *hdr = calloc(1,sizeof(*hdr));
    else
        *hdr = malloc(sizeof(*hdr));

    if (!*hdr)
        return ALLOCATION_PROBLEM;

    initHdr(*hdr);
    return HSTCAL_OK;
}

int allocSingleGroupExts(SingleGroup *x, int i, int j, unsigned extension, Bool zeroInitialize)
{
    if (allocSingleGroupHeader(&x->globalhdr, zeroInitialize))
        return ALLOCATION_PROBLEM;

    if (extension & SCIEXT)
    {
        if (allocFloatHdrData(&(x->sci),i,j, zeroInitialize))
            return ALLOCATION_PROBLEM;
    }
    if (extension & ERREXT)
    {
        if (allocFloatHdrData(&(x->err),i,j, zeroInitialize))
            return ALLOCATION_PROBLEM;
    }
    if (extension & DQEXT)
    {
        if (allocShortHdrData(&(x->dq),i,j, zeroInitialize))
            return ALLOCATION_PROBLEM;
    }
    return 0;
}

void setStorageOrder(SingleGroup * group, enum StorageOrder storageOrder)
{
    if (!group)
        return;

    group->sci.data.storageOrder = storageOrder;
    group->err.data.storageOrder = storageOrder;
    group->dq.data.storageOrder = storageOrder;

}

void copyOffsetFloatData(float * output, const float * input,
        unsigned nRows, unsigned nColumns,
        unsigned outputOffset, unsigned inputOffset,
        unsigned outputSkipLength, unsigned inputSkipLength)
{
    //WARNING - assumes row major storage
    {unsigned ithRow;
#ifdef _OPENMP
    #pragma omp parallel for shared(output, input) private(ithRow) schedule(static)
#endif
    for (ithRow = 0; ithRow < nRows; ++ithRow)
        memcpy(output + outputOffset + ithRow*outputSkipLength, input + inputOffset + ithRow*inputSkipLength, nColumns*sizeof(*output));
    }
}
void copyOffsetShortData(short * output, const short * input,
        unsigned nRows, unsigned nColumns,
        unsigned outputOffset, unsigned inputOffset,
        unsigned outputSkipLength, unsigned inputSkipLength)
{
    //WARNING - assumes row major storage
    {unsigned ithRow;
#ifdef _OPENMP
    #pragma omp parallel for shared(output, input) private(ithRow) schedule(static)
#endif
    for (ithRow = 0; ithRow < nRows; ++ithRow)
        memcpy(output + outputOffset + ithRow*outputSkipLength, input + inputOffset + ithRow*inputSkipLength, nColumns*sizeof(*output));
    }
}

void copyOffsetSingleGroup(SingleGroup * output, const SingleGroup * input,
        unsigned nRows, unsigned nColumns,
		unsigned outputOffset, unsigned inputOffset,
		unsigned outputSkipLength, unsigned inputSkipLength)
{
    if (!output || !input)
        return;
    //WARNING - assumes row major storage
    assert(output->sci.data.storageOrder == ROWMAJOR && output->sci.data.storageOrder == ROWMAJOR);

    //sci data
    if (output->sci.data.data && input->sci.data.data)
        copyOffsetFloatData(output->sci.data.data, input->sci.data.data, nRows, nColumns, outputOffset, inputOffset, outputSkipLength, inputSkipLength);
    //err data
    if (output->err.data.data && input->err.data.data)
        copyOffsetFloatData(output->err.data.data, input->err.data.data, nRows, nColumns, outputOffset, inputOffset, outputSkipLength, inputSkipLength);
    //dq data
    if (output->dq.data.data && input->dq.data.data)
        copyOffsetShortData(output->dq.data.data, input->dq.data.data, nRows, nColumns, outputOffset, inputOffset, outputSkipLength, inputSkipLength);
}


int copySingleGroup(SingleGroup * target, const SingleGroup * source, enum StorageOrder targetStorageOrder)
{
    //WARNING assumes target pre allocated and initialized (entire tree). This way data can be copied to pre
    //allocated target, i.e. copy(a, b) .. do something .. copy(b, a)
    //NOTE: If structs contained total size we could just use malloc & memcpy and be done with it.

    if (!target || !source)
        return ALLOCATION_PROBLEM;

    setStorageOrder(target, targetStorageOrder);

    if (source->filename)
    {
        size_t filenameLength = strlen(source->filename)+1;
        if (!target->filename || (target->filename && strlen(target->filename) != filenameLength))
        {
            if (target->filename)
                free(target->filename);
            target->filename = malloc(filenameLength*sizeof(*source->filename));
        }
        if (!target->filename)
        {
            initSingleGroup(target);
            return ALLOCATION_PROBLEM;
        }
        memcpy(target->filename, source->filename, filenameLength);
    }

    target->group_num = source->group_num;

    copyHdr(target->globalhdr, source->globalhdr); //This allocates

    if (source->sci.data.data)
    {
        if (copyFloatHdrData(&target->sci, &source->sci, targetStorageOrder))
        {
            initSingleGroup(target);
            return ALLOCATION_PROBLEM;
        }
    }
    if (source->err.data.data)
    {
        if (copyFloatHdrData(&target->err, &source->err, targetStorageOrder))
        {
            initSingleGroup(target);
            return ALLOCATION_PROBLEM;
        }
    }
    if (source->dq.data.data)
    {
        if (copyShortHdrData(&target->dq, &source->dq, targetStorageOrder))
        {
            initSingleGroup(target);
            return ALLOCATION_PROBLEM;
        }
    }
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
        if (allocFloatHdrData(&(x->sci),i,j, True)) return -1;
        if (allocFloatHdrData(&(x->err),i,j, True)) return -1;
        if (allocShortHdrData(&(x->dq),i,j, True)) return -1;
        if (allocShortHdrData(&(x->smpl),i,j, True)) return -1;
        if (allocFloatHdrData(&(x->intg),i,j, True)) return -1;
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
        assert(x);
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

/*
** Section 7.
**
** High-level functions formerly in hstioirf.c
 */

/* CFITSIO TODO: store axes in IODesc object as array to more
   conveniently interface with CFITSIO */

static void detect_iraferr(void) {
        sprintf(error_msg,"\nIRAF error %d: %s\n",c_iraferr(),
                c_iraferrmsg());
}

/*
** Make_iodesc takes a filename, extname, and extver and creates and
** initializes an IODesc structure.  In the process, it builds a
** correct filename to be used in the open statement to IRAF.  This
** constructed filename is returned.
*/
static char *make_iodesc(IODesc **x, char *fname, char *ename, int ever) {
        int i, n, flen;
        char *tmp;
        IODesc *iodesc;
        char xname[9];

        iodesc = (IODesc *)calloc(1,sizeof(IODesc));
        if (iodesc == NULL) {
            error(NOMEM,"Allocating I/O descriptor");
            return NULL;
        }
        iodesc->ff = NULL;
        iodesc->filename = NULL;
        iodesc->extname = NULL;
        iodesc->extver = 0;
        iodesc->hflag = 0;
        iodesc->hdr = NULL;
        iodesc->dims[0] = 0;
        iodesc->dims[1] = 0;
        iodesc->type = 0;
        if (fname == 0) fname = "";
        if (ename == 0) ename = "";
        iodesc->filename = (char *)calloc(((flen = strlen(fname)) + 1),
                sizeof(char));
        if (iodesc->filename == NULL) {
            free(iodesc);
            error(NOMEM,"Allocating I/O descriptor");
            return NULL;
        }
        n = strlen(ename);
        if (n > 8) { ioerr(BADEXTNAME,iodesc,0); return NULL; }
        for (i = 0; i < n; ++i)
            xname[i] = toupper(ename[i]);
        for (--i; i >= 0 && xname[i] == ' '; --i) ;
        ++i;
        xname[i] = '\0';
        iodesc->extname = (char *)calloc((strlen(xname) + 1),sizeof(char));
        if (iodesc->extname == NULL) {
            free(iodesc->filename);
            free(iodesc);
            error(NOMEM,"Allocating I/O descriptor");
            return NULL;
        }
        strcpy(iodesc->filename,fname);
        strcpy(iodesc->extname,xname);
        iodesc->extver = ever;

        /* make up the proper filename */
        /* check for a request for the primary HDU */
        tmp = (char *)calloc((flen + 80),sizeof(char));
        if (tmp == NULL) { error(NOMEM,"Allocating I/O descriptor"); return NULL; }
        strcpy(tmp,fname);
        if (ever == 0 || ename == 0 || ename[0] == '\0' || ename[0] == ' ')
            strcat(tmp,"[0]");
        else
            sprintf(&tmp[flen],"[%s,%d]",xname,ever);

        *x = iodesc;
        return tmp;
}

IODescPtr openInputImage(char *fname, char *ename, int ever) {
        IODesc *iodesc;
        int no_dims;
        char *tmp;
        char ospath[SZ_PATHNAME];
        int open_mode;
        int status = 0;

        /* CFITSIO: Error handling */
        c_pusherr(detect_iraferr);

        tmp = make_iodesc(&iodesc, fname, ename, ever);
        if (tmp == NULL) return NULL;
        iodesc->options = ReadOnly;

        /* CFITSIO: Resolve this inheritance stuff */
        /* p = strstr(tmp,"[0]"); */
        /* if (p == NULL) { */
        /*     tmp[strlen(tmp) - 1] = '\0'; */
        /*     strcat(tmp, ",NOINHERIT]"); */
        /* } */

        /* open the file using CFITSIO */
        if (c_vfn2osfn(fname, ospath)) {
            free(tmp);
            return NULL;
        }

        open_mode = READONLY;

        if (c_vfn2osfn(tmp, ospath)) {
            free(tmp);
            return NULL;
        }
        free(tmp);

        if (fits_open_file(&iodesc->ff, ospath, open_mode, &status)) {
            ioerr(BADOPEN, iodesc, status);
            free(iodesc->extname);
            free(iodesc->filename);
            free(iodesc);
            return NULL;
        }

        /* get the dimensions and type */
        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        fits_get_img_equivtype(iodesc->ff, &iodesc->type, &status);
        fits_get_img_size(iodesc->ff, 2, iodesc->dims, &status);
        if (status) {
            ioerr(BADDIMS, iodesc, status);
            return NULL;
        }
        if (no_dims == 2) {
            /* Nothing */
        } else if (no_dims == 1) {
            iodesc->dims[1] = 0;
        } else if (no_dims == 0) {
            iodesc->dims[0] = 0;
            iodesc->dims[1] = 0;
        } else {
            ioerr(BADDIMS, iodesc, 0);
            return NULL;
        }

        clear_err();

        return iodesc;
}

IODescPtr openOutputImage(char *fname, char *ename, int ever, Hdr *hd,
        int d1, int d2, FitsDataType typ) {
        IODesc *iodesc;
        char *tmp;
        char date[12];
        char date_card[81];
        time_t t;
        struct tm *time_tmp;
        FitsKw kw;
        char ename_val[9];
        int ever_val;
        char ospath[SZ_PATHNAME];
        int status = 0;

        /* CFITSIO: Error handling */
        c_pusherr(detect_iraferr);

        tmp = make_iodesc(&iodesc, fname, ename, ever);
        if (tmp == NULL) return NULL;
        iodesc->options = WriteOnly;
        if (ever == 0 || ename == 0 || ename[0] == '\0' || ename[0] == ' ') {
            int rtn = ckNewFile(fname);
            if (rtn == 1) {
                ioerr(BADEXIST, iodesc, 0);
                return NULL;
            } else if (rtn == 2) {
                ioerr(BADREMOVE, iodesc, 0);
                return NULL;
            }
        }

        /* CFITSIO: Check this INHERIT, APPEND nonsense works */
        /* p = strstr(tmp, "[0]"); */
        /* if (p == NULL) { */
        /*     tmp[strlen(tmp) - 1] = '\0'; */
        /*     strcat(tmp, ",INHERIT,APPEND]"); */
        /* } */
        /* else { */
        /*     tmp[strlen(tmp) - 3] = '\0'; /\* eliminate the "[0]" *\/ */
        /* } */

        /* make sure ename and ever are in the header array */
        kw = findKw(hd, "EXTNAME");
        if (kw == NotFound) {
            if (ever != 0 && ename != 0 &&
                ename[0] != '\0' && ename[0] != ' ') {
                kw = insertfirst(hd);
                kw = insertStringKw(kw, "EXTNAME", ename, "Name of the extension");
            }
        } else {
            /* Make sure it has the right value */
            getStringKw(kw,ename_val,8);
            if (strncpy(ename_val, ename, strlen(ename)) != 0)
                putStringKw(kw,ename);
        }
        kw = findKw(hd,"EXTVER");
        if (kw == NotFound) {
            if (ever != 0 && ename != 0 &&
                ename[0] != '\0' && ename[0] != ' ') {
                kw = findKw(hd, "EXTNAME");
                kw = insertIntKw(kw, "EXTVER", ever, "Extension version");
            }
        } else {
            /* Make sure it has the right value */
            ever_val = getIntKw(kw);
            if (ever != ever_val)
                putIntKw(kw, ever);
        }

        /* open or create the file using CFITSIO */
        if (ever == 0 || ename == 0 || ename[0] == '\0' || ename[0] == ' ') {
            c_vfn2osfn(fname, ospath);
            fits_create_file(&iodesc->ff, ospath, &status);
        } else {
            c_vfn2osfn(fname, ospath);
            fits_open_file(&iodesc->ff, ospath, READWRITE, &status);
        }
        if (status) {
            ioerr(BADOPEN, iodesc, status);
            free(iodesc->extname);
            free(iodesc->filename);
            free(iodesc);
            return NULL;
        }
        free(tmp);

        iodesc->dims[0] = d1;
        iodesc->dims[1] = d2;
        /* IMIO would always set bitpix to 16 when the dimensions are
           naught */
        if (d1 == 0 && d2 == 0) {
            iodesc->type = SHORT_IMG;
        } else {
            switch (typ) {
            case FITSBYTE:
                iodesc->type = BYTE_IMG;
                break;
            case FITSSHORT:
                iodesc->type = SHORT_IMG;
                break;
            case FITSLONG:
                iodesc->type = LONG_IMG;
                break;
            case FITSFLOAT:
                iodesc->type = FLOAT_IMG;
                break;
            case FITSDOUBLE:
                iodesc->type = DOUBLE_IMG;
                break;
            default:
                iodesc->type = SHORT_IMG;
                break;
            }
        }
        iodesc->hdr = hd;

        if (fits_create_img(iodesc->ff, iodesc->type, 2, iodesc->dims, &status)) {
            ioerr(BADOPEN, iodesc, status);
            return NULL;
        }

        if (fits_write_record(iodesc->ff, "ORIGIN  = 'HSTIO/CFITSIO March 2010' / FITS file originator", &status)) {
            ioerr(BADWRITE, iodesc, status);
            return NULL;
        }

        t = time(NULL);
        time_tmp = localtime(&t);
        strftime(date, 12, "%Y-%m-%d", time_tmp);
        snprintf(date_card, 80,
                 "DATE    = '%s' / date this file was written (yyyy-mm-dd)", date);
        if (fits_write_record(iodesc->ff, date_card, &status)) {
            ioerr(BADWRITE, iodesc, status);
            return NULL;
        }

        iodesc->hflag = 1; /* mark to write header */
        if (iodesc->dims[0] == 0) {
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        clear_err();

        return iodesc;
}

IODescPtr openUpdateImage(char *fname, char *ename, int ever, Hdr *hd) {
        IODesc *iodesc;
        int no_dims;
        char *tmp;
        char ospath[SZ_PATHNAME];
        int status = 0;

        /* CFITSIO: Error handling */
        c_pusherr(detect_iraferr);

        tmp = make_iodesc(&iodesc, fname, ename, ever);
        if (tmp == NULL) return NULL;
        iodesc->options = ReadWrite;

        /* CFITSIO: Resolve this keyword inheritance stuff */
        /* p = strstr(tmp,"[0]"); */
        /* if (p == NULL) { */
        /*     tmp[strlen(tmp) - 1] = '\0'; */
        /*     strcat(tmp,",NOINHERIT]"); */
        /* } */

        /* open the file using CFITSIO */
        c_vfn2osfn(tmp, ospath);

        if (fits_open_file(&iodesc->ff, ospath, READWRITE, &status)) {
            ioerr(BADOPEN, iodesc, status);
            free(tmp);
            free(iodesc->extname);
            free(iodesc->filename);
            free(iodesc);
            return NULL;
        }
        free(tmp);

        /* get the dimensions and type */
        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        fits_get_img_equivtype(iodesc->ff, &iodesc->type, &status);
        fits_get_img_size(iodesc->ff, 2, iodesc->dims, &status);
        if (status) {
            ioerr(BADDIMS, iodesc, status);
            return NULL;
        }
        if (no_dims == 2) {
            /* Nothing */
        } else if (no_dims == 1) {
            iodesc->dims[1] = 0;
        } else if (no_dims == 0) {
            iodesc->dims[0] = 0;
            iodesc->dims[1] = 0;
        } else {
            ioerr(BADDIMS, iodesc, 0);
            return NULL;
        }

        /* read the user area into the header array */
        getHeader(iodesc, hd);

        clear_err();
        return iodesc;
}

void closeImage(IODescPtr iodesc_) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int status = 0;

        if (iodesc->options != ReadOnly && iodesc->dims[0] != 0)
            putHeader(iodesc);

        if (fits_close_file(iodesc->ff, &status)) {
            /* TODO: Raise error */
        }

        /* This is a handy check to use pyfits to validate the file upon every close */
        /* c_vfn2osfn(iodesc->filename, ospath); */
        /* sprintf(system_string, "python -c \"import pyfits; pyfits.open('%s')\"", ospath); */
        /* if (system(system_string)) { */
        /*   printf("LOG: pyfits corruption!!!\n"); */
        /*   exit(1); */
        /* } */

        /* if (there is an IRAF error) */
        /*      ioerr(IRAF_CLOSE,iodesc); */
        free(iodesc->extname);
        free(iodesc->filename);
        free(iodesc);

        /* CFITSIO: Error handling */
        c_poperr();
}


/* According to the imio documentation, the following reserved
   keywords are recognized:

   SIMPLE BITPIX DATATYPE NAXIS* GROUPS GCOUNT PCOUNT PSIZE
   PTYPE* PDTYPE* PSIZE* XTENSION
*/

static char* reservedKwds[] = {
    "BITPIX  ",
    "BSCALE  ",
    "BZERO   ",
    "DATAMAX ",
    "DATAMIN ",
    "DATATYPE",
    "DATE    ",
    "EXTEND  ",
    "GCOUNT  ",
    "GROUPS  ",
    "NAXIS   ",
    "NAXIS*  ",
    "ORIGIN  ",
    "PCOUNT  ",
    "PDTYPE* ",
    "PSIZE   ",
    "PSIZE*  ",
    "PTYPE*  ",
    "SIMPLE  ",
    "XTENSION",
    NULL
};

/* Whole cards to not allow in the user area.

   Must remain alphabetized.
*/
static char* reservedCards[] = {
    "COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy",
    "COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H ",
    NULL
};

int isReservedKwd(const char* card) {
        /* CFITSIO: Should this be made case-insensitive? */
        /* TODO: Maybe use a binary search? */

        /* Returns 1 if the card matches one of the reserved keywords */
        int i;
        int match;
        char** kwd = reservedKwds;
        int cmp;

        for (kwd = reservedCards; *kwd != NULL; ++kwd) {
            cmp = strncmp(*kwd, card, 79);
            if (cmp == 0) {
                return 1;
            }
        }

        for (kwd = reservedKwds; *kwd != NULL; ++kwd) {
            /* Short-circuit if we're certain not to find the kwd
               later in the list */
            if ((*kwd)[0] > card[0]) {
                return 0;
            }

            match = 1;
            for (i = 0; i < 8; ++i) {
                /* '*' indicates a digit */
                if ((*kwd)[i] == '*') {
                    if (card[i] < '0' || card[i] > '9') {
                        match = 0;
                        break;
                    }
                } else if ((*kwd)[i] != card[i]) {
                    match = 0;
                    break;
                }
            }
            if (match) {
                return 1;
            }
        }

        return 0;
}

int getHeader(IODescPtr iodesc_, Hdr *hd) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int ncards, i, j;
        char source[HDRSize];
        char *target;
        int status = 0;

        if (iodesc->options == WriteOnly) {
            ioerr(NOGET, iodesc, 0);
            return -1;
        }

        /* get the number of cards in the header */
        if (fits_get_hdrspace(iodesc->ff, &ncards, NULL, &status)) {
            ioerr(BADHSIZE, iodesc, status);
            return -1;
        }

        /* allocate space for the header cards */
        if (allocHdr(hd, ncards, True) == -1) return -1;

        /* translate the data */
        hd->nlines = 0;
        for (i = 0; i < ncards; ++i) {
            if (fits_read_record(iodesc->ff, i+1, source, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
            if (!isReservedKwd(source)) {
                target = hd->array[hd->nlines];
                for (j = 0; j < (HDRSize -1); ++j) {
                    *target++ = source[j];
                }
                *target++ = '\0';
                hd->nlines++;
            }
        }
        iodesc->hdr = hd;

        clear_err();
        return 0;
}

int putHeader(IODescPtr iodesc_) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j, tmp;
        int numkeys;
        int found_non_space;
        char *source;
        char card[81];
        int status = 0;

        if (iodesc->options == ReadOnly) {
            ioerr(NOPUT, iodesc, status);
            return -1;
        }

        if (iodesc->hflag) {
            /* CFITSIO: We probably need to move this in front of all
               calls to fits_create_img */

            /* If the image is actually 1-dimensional, modify the naxis2
             * value so the output header is written with only NAXIS and
             * NAXIS1 keywords, where NAXIS=1, and NAXIS1=number.
             */
            if (iodesc->dims[0] != 0 && iodesc->dims[1] == 1)
                iodesc->dims[1] = 0;

            /* set the pixel type */
            fits_update_key(iodesc->ff, TINT, "BITPIX", &(iodesc->type), NULL, &status);
            if (status) {
                ioerr(BADWRITE, iodesc, status);
                return -1;
            }
            if (iodesc->dims[0] == 0 && iodesc->dims[1] == 0) {
                tmp = 0;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &tmp, NULL, &status);
                if (status) {
                    ioerr(BADWRITE, iodesc, status);
                    return -1;
                }
                fits_delete_key(iodesc->ff, "NAXIS1", &status);
                if (status == KEY_NO_EXIST) {
                    fits_clear_errmsg();
                    status = 0;
                }
                fits_delete_key(iodesc->ff, "NAXIS2", &status);
                if (status == KEY_NO_EXIST) {
                    fits_clear_errmsg();
                    status = 0;
                }
            } else if (iodesc->dims[0] != 0 && iodesc->dims[1] == 0) {
                /* set the number of dimensions */
                tmp = 1;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &tmp, NULL, &status);
                if (status) {
                    ioerr(BADWRITE, iodesc, status);
                    return -1;
                }
                /* set dim1 */
                fits_update_key(iodesc->ff, TINT, "NAXIS1", &iodesc->dims[0], NULL, &status);
                if (status) {
                    ioerr(BADWRITE, iodesc, status);
                    return -1;
                }
                fits_delete_key(iodesc->ff, "NAXIS2", &status);
                if (status == KEY_NO_EXIST) {
                    fits_clear_errmsg();
                    status = 0;
                }
            } else {
                /* set the number of dimensions */
                tmp = 2;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &tmp, NULL, &status);
                /* set dim1 and dim2 */
                fits_update_key(iodesc->ff, TINT, "NAXIS1", &iodesc->dims[0], NULL, &status);
                fits_update_key(iodesc->ff, TINT, "NAXIS2", &iodesc->dims[1], NULL, &status);
            }

            if (status) {
                ioerr(BADWRITE, iodesc, status);
            }
        }

        /* Verify the size of the user area */
        /* The original code just memcopies the cards into the "user
           area" of the header.  CFITSIO doesn't have the concept of a
           "user area", so we need to carefully only copy the cards
           that are not "reserved".
        */
        if (fits_get_hdrspace(iodesc->ff, &numkeys, NULL, &status)) {
            ioerr(BADWRITE, iodesc, status);
            return -1;
        }

        for (i = 0, j = numkeys; i < numkeys; ++i, --j) {
            if (fits_read_record(iodesc->ff, j, card, &status)) {
                ioerr(BADWRITE, iodesc, status);
                return -1;
            }
            if (!isReservedKwd(card)) {
                if (fits_delete_record(iodesc->ff, j, &status)) {
                    ioerr(BADWRITE, iodesc, status);
                    return -1;
                }
            } else {
                ++j;
            }
        }

        /* translate the data */

        /* Skip blank cards at the beginning */
        found_non_space = 0;
        for (i = 0; i < iodesc->hdr->nlines; ++i) {
            source = iodesc->hdr->array[i];
            for (j = 0; j < 80; ++j) {
                if (source[j] != ' ' &&
                    source[j] != '\n' &&
                    source[j] != 0) {
                    found_non_space = 1;
                    break;
                }
            }
            if (found_non_space) {
                break;
            }
        }

        for (/* i from above */; i < iodesc->hdr->nlines; ++i) {
            source = iodesc->hdr->array[i];
            if (!isReservedKwd(source)) {
                if (fits_write_record(iodesc->ff, source, &status)) {
                    ioerr(BADWRITE, iodesc, status);
                    return -1;
                }
            }
        }

        /* If we don't explicitly set BSCALE and BZERO to 1.0 and 0.0
           here, their values could be inadvertently brought over from
           the source image.  This was the source of a very
           hard-to-find bug. */
        if (iodesc->type == TFLOAT || iodesc->type == TDOUBLE) {
            fits_set_bscale(iodesc->ff, 1.0, 0.0, &status);
        }

        clear_err();
        return 0;
}

int getFloatData(IODescPtr iodesc_, FloatTwoDArray *da) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int no_dims, i, j;
        long fpixel[2];
        int anynul;
        int type;
        FitsKw kw;
        float val;
        int status = 0;

        if (iodesc->options == WriteOnly) { ioerr(NOGET,iodesc, 0); return -1; }

        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        fits_get_img_size(iodesc->ff, 2, iodesc->dims, &status);
        if (status) {
            ioerr(BADDIMS, iodesc, status);
            return -1;
        }

        /*
           If the number  of dimensions of the image is zero, need to
           determine how many dimensions the image is supposed to have
           according to the NPIX[1/2] keyword(s).
        */
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc, 0); return -1; }
            val = getFloatKw(kw);

            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc, 0); return -1; }
            iodesc->dims[0] = getIntKw(kw);

            /* If NPIX2 is not found, the image should be 1D; dim2 = 1 and *
             * not 0 for purposes of memory allocation.                    */
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw == 0)  {
                iodesc->dims[1] = 1;
            } else {
                iodesc->dims[1] = getIntKw(kw);
            }

            if (allocFloatData(da, iodesc->dims[0], iodesc->dims[1], False)) return -1;
            for (j = 0; j < iodesc->dims[1]; ++j) {
                for (i = 0; i < iodesc->dims[0]; ++i) {
                    PPix(da, i, j) = val;
                }
            }
        } else if (no_dims == 1) {
            iodesc->dims[1] = 1;
            fits_get_img_equivtype(iodesc->ff, &type, &status);
            /* CFITSIO TODO: Should we verify the type is correct
               here?  Original code gets type, but then does nothing
               with it. */
            if (allocFloatData(da, iodesc->dims[0], iodesc->dims[1], True)) return -1;
            fpixel[0] = 1;
            fpixel[1] = 1;
            if (fits_read_pix(iodesc->ff, TFLOAT, fpixel, iodesc->dims[0], 0,
                              (float *)&(PPix(da, 0, 0)), &anynul, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
        } else if (no_dims == 2) {
            fits_get_img_equivtype(iodesc->ff, &type, &status);
            /* CFITSIO TODO: Should we verify the type is correct
               here?  Original code gets type, but then does nothing
               with it. */
            if (allocFloatData(da, iodesc->dims[0], iodesc->dims[1], True)) return -1;

            fpixel[0] = 1;
            if (da->storageOrder == ROWMAJOR)
            {
                for (i = 0; i < iodesc->dims[1]; ++i) {
                    fpixel[1] = i + 1;
                    if (fits_read_pix(iodesc->ff, TFLOAT, fpixel, iodesc->dims[0], 0,
                            &(PPix(da, 0, i)), &anynul, &status)) {
                        ioerr(BADREAD,iodesc, status);
                        return -1;
                    }
                }
            }
            else
            {
                unsigned nColumns = iodesc->dims[0];
                float * row = malloc(nColumns*sizeof(float));
                if (!row)
                    return OUT_OF_MEMORY;
                for (i = 0; i < iodesc->dims[1]; ++i)
                {
                    fpixel[1] = i + 1;
                    if (fits_read_pix(iodesc->ff, TFLOAT, fpixel, nColumns, 0,
                            row, &anynul, &status)) {
                        ioerr(BADREAD,iodesc, status);
                        return -1;
                    }
                    {unsigned j;
                    for (j = 0; j < nColumns; ++j)
                        PPixColumnMajor(da, i, j) = row[j];
                    }
                }
                if (row)
                    free(row);
            }
        } else {
            ioerr(BADDIMS,iodesc,0);
            return -1;
        }

        clear_err();
        return 0;
}

int putFloatData(IODescPtr iodesc_, FloatTwoDArray *da) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j;
        float tmp;
        long fpixel[2];
        FitsKw kw;
        int is_eq;
        int naxis;
        long dims[2];
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }

        /* check for a constant array, if not SCI data */
        if (strcmp(iodesc->extname,"SCI") != 0
            && da->tot_nx != 0 && da->tot_ny != 0) {
            tmp = PPix(da,0,0);
            for (i = 0, is_eq = 1; (i < da->tot_nx) && is_eq; ++i) {
                for (j = 0; (j < da->tot_ny); ++j) {
                    if (PPix(da,i,j) != tmp) {
                        is_eq = 0;
                        break;
                    }
                }
            }
            if (is_eq) {
                /* This is a constant array. */
                /* add NPIX1, NPIX2 (if necessary), and PIXVALUE keywords */
                kw = findKw(iodesc->hdr,"PIXVALUE");
                if (kw == 0) /* add it */
                    addFloatKw(iodesc->hdr,"PIXVALUE",tmp,
                        "values of pixels in constant array");
                else
                    putFloatKw(kw,tmp);

                kw = findKw(iodesc->hdr,"NPIX1");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"NPIX1",iodesc->dims[0],
                        "length of constant array axis 1");
                else
                    putIntKw(kw,iodesc->dims[0]);

                /* NPIX2 should only be added if the y-dimension is > 1. */
                if (da->tot_ny > 1) {
                    kw = findKw(iodesc->hdr,"NPIX2");
                    if (kw == 0) /* add it */
                        addIntKw(iodesc->hdr,"NPIX2",iodesc->dims[1],
                            "length of constant array axis 2");
                    else
                        putIntKw(kw,iodesc->dims[1]);
                }

                naxis = 0;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &naxis, NULL, &status);
                iodesc->dims[0] = 0;
                iodesc->dims[1] = 0;

                if (fits_resize_img(iodesc->ff, FLOAT_IMG, 0, iodesc->dims, &status)) {
                    ioerr(BADWRITE, iodesc, status); return -1;
                }

                /* update the header, etc. */
                if (iodesc->hflag) {
                    iodesc->type = FLOAT_IMG;
                    putHeader(iodesc);
                    iodesc->hflag = 0;
                }
                fits_flush_file(iodesc->ff, &status);

                return 0;
            }
        }

        /* If not a constant array, make sure NPIX1, NPIX2, and PIXVALUE *
         * are NOT present in the header to be written out.              */
        kw = findKw(iodesc->hdr,"NPIX1");
        if (kw != 0) /* remove it */
            delKw(kw);

        if (da->tot_ny > 1) {
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw != 0) /* remove it */
                delKw(kw);
        }

        kw = findKw(iodesc->hdr,"PIXVALUE");
        if (kw != 0) /* remove it */
            delKw(kw);

        /* Get the current CFITSIO size, and if it's different, resize it */
        fits_get_img_size(iodesc->ff, 2, dims, &status);
        if (dims[0] != da->nx || dims[1] != da->ny) {
            iodesc->dims[0] = da->nx;
            iodesc->dims[1] = da->ny;
            if (fits_resize_img(iodesc->ff, FLOAT_IMG, 2, iodesc->dims, &status)) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }
        }

        /* update the header area */
        if (iodesc->hflag) {
            iodesc->type = FLOAT_IMG;
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        fpixel[0] = 1;
        for (i = 0; i < da->ny; ++i) {
            fpixel[1] = i + 1;
            if (fits_write_pix(iodesc->ff, TFLOAT, fpixel, da->nx,
                               (float *)&(PPix(da, 0, i)), &status)) {
                ioerr(BADWRITE, iodesc, status);
                return -1;
            }
        }

        fits_flush_file(iodesc->ff, &status);

        clear_err();
        return 0;
}

/*                                                                     **
** Write output a subsection of an image in memory to a file where the **
** subsection is the full size of the output data.                     **
**                                                                     */
int putFloatSect(IODescPtr iodesc_, FloatTwoDArray *da, int xbeg,
                 int ybeg, int xsize, int ysize) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j, xend, yend;
        float tmp;
        FitsKw kw;
        long fpixel[2];
        long dims[2];
        int is_eq;
        int naxis;
        int status = 0;

        /* CFITSIO: Verify that the section is within range? */

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc, 0); return -1; }

        xend = xbeg + xsize;
        yend = ybeg + ysize;
        /* check for a constant array, if not SCI data */
        if (strcmp(iodesc->extname,"SCI") != 0
            && da->tot_nx != 0 && da->tot_ny != 0) {
            tmp = PPix(da, 0, 0);
            for (i = xbeg, is_eq = 1; (i < xend) && is_eq; ++i) {
                for (j = ybeg; (j < yend); ++j) {
                    if (PPix(da,i,j) != tmp) {
                        is_eq = 0;
                        break;
                    }
                }
            }
            if (is_eq) {
                /* This is a constant array. */
                /* add NPIX1, NPIX2 (if necessary), and PIXVALUE keywords */
                kw = findKw(iodesc->hdr,"PIXVALUE");
                if (kw == 0) /* add it */
                    addFloatKw(iodesc->hdr,"PIXVALUE",tmp,
                        "values of pixels in constant array");
                else
                    putFloatKw(kw,tmp);

                kw = findKw(iodesc->hdr,"NPIX1");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"NPIX1",iodesc->dims[0],
                        "length of constant array axis 1");
                else
                    putIntKw(kw,iodesc->dims[0]);

                /* NPIX2 should only be added if the y-dimension is > 1. */
                if (da->tot_ny > 1) {
                    kw = findKw(iodesc->hdr,"NPIX2");
                    if (kw == 0) /* add it */
                        addIntKw(iodesc->hdr,"NPIX2",iodesc->dims[1],
                            "length of constant array axis 2");
                    else
                        putIntKw(kw,iodesc->dims[1]);
                }

                naxis = 0;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &naxis, NULL, &status);
                iodesc->dims[0] = 0;
                iodesc->dims[1] = 0;
                /* update the header, etc. */
                if (iodesc->hflag) {
                    iodesc->type = FLOAT_IMG;
                    putHeader(iodesc);
                    iodesc->hflag = 0;
                }

                fits_flush_file(iodesc->ff, &status);

                clear_err();
                return 0;
            }
        }

        /* If not a constant array, make sure NPIX1, NPIX2, and PIXVALUE *
         * are NOT present in the header to be written out.              */
        kw = findKw(iodesc->hdr,"PIXVALUE");
        if (kw != 0) /* remove it */
            delKw(kw);

        kw = findKw(iodesc->hdr,"NPIX1");
        if (kw != 0) /* remove it */
            delKw(kw);

        if (da->tot_ny > 1) {
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw != 0) /* remove it */
                delKw(kw);
        }

        /* Get the current CFITSIO size, and if it's different, resize it */
        fits_get_img_size(iodesc->ff, 2, dims, &status);
        if (dims[0] != xend - xbeg || dims[1] != yend - ybeg) {
            iodesc->dims[0] = xend - xbeg;
            iodesc->dims[1] = yend - ybeg;
            if (fits_resize_img(iodesc->ff, FLOAT_IMG, 2, iodesc->dims, &status)) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }

            /* Note, we don't need to fill the image with the constant value, since
               the image will be entirely over-written by the passed in array da */
        }

        /* update the header area */
        if (iodesc->hflag) {
            iodesc->type = FLOAT_IMG;
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        fpixel[0] = 1;
        for (i = ybeg; i < yend; ++i) {
            fpixel[1] = i - ybeg + 1;
            if (fits_write_pix(iodesc->ff, TFLOAT, fpixel, xsize,
                               (float*)&(PPix(da, xbeg, i)), &status)) {
                ioerr(BADWRITE, iodesc, status);
                return -1;
            }
        }

        fflush(stdout);

        clear_err();
        return 0;
}

int getShortData(IODescPtr iodesc_, ShortTwoDArray *da) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int no_dims, i, j;
        FitsKw kw;
        short val;
        long fpixel[2];
        int anynul = 0;
        int status = 0;

        if (iodesc->options == WriteOnly) { ioerr(NOGET,iodesc, 0); return -1; }

        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        fits_get_img_size(iodesc->ff, 2, iodesc->dims, &status);
        if (status) {
            ioerr(BADDIMS, iodesc, status);
            return -1;
        }

        /*
           If the number  of dimensions of the image is zero, need to
           determine how many dimensions the image is supposed to have
           according to the NPIX[1/2] keyword(s).
        */
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            val = getIntKw(kw);

            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            iodesc->dims[0] = getIntKw(kw);

            /* If NPIX2 is not found, the image should be 1D; dim2 = 1 and *
             * not 0 for purposes of memory allocation.                    */
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw == 0)  {
                iodesc->dims[1] = 1;
            } else {
                iodesc->dims[1] = getIntKw(kw);
            }

            if (allocShortData(da, iodesc->dims[0], iodesc->dims[1], True)) return -1;
            for (j = 0; j < iodesc->dims[1]; ++j)
                for (i = 0; i < iodesc->dims[0]; ++i)
                    PPix(da, i, j) = val;
        } else if (no_dims == 1) {
            iodesc->dims[1] = 1;
            /* CFITSIO TODO: Should we verify the type is correct
               here?  Original code gets type, but then does nothing
               with it. */
            if (allocShortData(da, iodesc->dims[0], iodesc->dims[1], True)) return -1;
            fpixel[0] = 1;
            fpixel[1] = 1;
            if (fits_read_pix(iodesc->ff, TSHORT, fpixel, iodesc->dims[0], NULL,
                              (short *)&(PPix(da, 0, 0)), &anynul, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
        } else if (no_dims == 2) {
            /* CFITSIO TODO: Should we verify the type is correct
               here?  Original code gets type, but then does nothing
               with it. */
            if (allocShortData(da, iodesc->dims[0], iodesc->dims[1], True)) return -1;
            fpixel[0] = 1;
            for (i = 0; i < iodesc->dims[1]; ++i) {
                fpixel[1] = i + 1;
                if (fits_read_pix(iodesc->ff, TSHORT, fpixel, iodesc->dims[0], NULL,
                                  (short *)&(PPix(da, 0, i)), &anynul, &status)) {
                    ioerr(BADREAD, iodesc, status);
                    return -1;
                }
            }
        } else {
            ioerr(BADDIMS, iodesc, 0);
            return -1;
        }

        clear_err();
        return 0;
}

int putShortData(IODescPtr iodesc_, ShortTwoDArray *da) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j;
        short tmp;
        long fpixel[2];
        FitsKw kw;
        int is_eq;
        int naxis;
        long dims[2];
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }

        /* check for a constant array, if not SCI data */
        if (strcmp(iodesc->extname,"SCI") != 0
            && da->tot_nx != 0 && da->tot_ny != 0) {
            tmp = PPix(da,0,0);
            for (i = 0, is_eq = 1; (i < da->tot_nx) && is_eq; ++i) {
                for (j = 0; (j < da->tot_ny); ++j) {
                    if (PPix(da,i,j) != tmp) {
                        is_eq = 0;
                        break;
                    }
                }
            }
            if (is_eq) {
                /* This is a constant array. */
                /* add NPIX1, NPIX2 (if necessary), and PIXVALUE keywords */
                kw = findKw(iodesc->hdr,"PIXVALUE");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"PIXVALUE",(int)tmp,
                        "values of pixels in constant array");
                else
                    putIntKw(kw,(int)tmp);

                kw = findKw(iodesc->hdr,"NPIX1");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"NPIX1",iodesc->dims[0],
                        "length of constant array axis 1");
                else
                    putIntKw(kw,iodesc->dims[0]);

                /* NPIX2 should only be added if the y-dimension is > 1. */
                if (da->tot_ny > 1) {
                    kw = findKw(iodesc->hdr,"NPIX2");
                    if (kw == 0) /* add it */
                        addIntKw(iodesc->hdr,"NPIX2",iodesc->dims[1],
                            "length of constant array axis 2");
                    else
                        putIntKw(kw,iodesc->dims[1]);
                }

                naxis = 0;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &naxis, NULL, &status);
                iodesc->dims[0] = 0;
                iodesc->dims[1] = 0;

                if (fits_resize_img(iodesc->ff, SHORT_IMG, 0, iodesc->dims, &status)) {
                    ioerr(BADWRITE, iodesc, status); return -1;
                }

                /* update the header, etc. */
                if (iodesc->hflag) {
                    iodesc->type = SHORT_IMG;
                    putHeader(iodesc);
                    iodesc->hflag = 0;
                }

                fits_flush_file(iodesc->ff, &status);

                clear_err();
                return 0;
            }
        }

        /* If not a constant array, make sure NPIX1, NPIX2, and PIXVALUE *
         * are NOT present in the header to be written out.              */
        kw = findKw(iodesc->hdr,"NPIX1");
        if (kw != 0) /* remove it */
            delKw(kw);

        if (da->tot_ny > 1) {
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw != 0) /* remove it */
                delKw(kw);
        }

        kw = findKw(iodesc->hdr,"PIXVALUE");
        if (kw != 0) /* remove it */
            delKw(kw);

        /* Get the current CFITSIO size, and if it's different, resize it */
        fits_get_img_size(iodesc->ff, 2, dims, &status);
        if (dims[0] != da->nx || dims[1] != da->ny) {
            iodesc->dims[0] = da->nx;
            iodesc->dims[1] = da->ny;
            if (fits_resize_img(iodesc->ff, SHORT_IMG, 2, iodesc->dims, &status)) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }
        }

        /* update the header area */
        if (iodesc->hflag) {
            iodesc->type = SHORT_IMG;
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        fpixel[0] = 1;
        for (i = 0; i < da->ny; ++i) {
            fpixel[1] = i + 1;
            if (fits_write_pix(iodesc->ff, TSHORT, fpixel, da->nx,
                               (short *)&(PPix(da, 0, i)), &status)) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }
        }

        fits_flush_file(iodesc->ff, &status);

        clear_err();
        return 0;
}

/*                                                                     **
** Write output a subsection of an image in memory to a file where the **
** subsection is the full size of the output data.                     **
**                                                                     */
int putShortSect(IODescPtr iodesc_, ShortTwoDArray *da, int xbeg, int ybeg,
                     int xsize, int ysize) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j, xend, yend;
        short tmp;
        FitsKw kw;
        int naxis;
        int is_eq;
        long fpixel[2];
        long dims[2];
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }

        xend = xbeg + xsize;
        yend = ybeg + ysize;
        /* check for a constant array, if not SCI data */
        if (strcmp(iodesc->extname,"SCI") != 0
            && da->tot_nx != 0 && da->tot_ny != 0) {
            tmp = PPix(da,0,0);
            for (i = xbeg, is_eq = 1; (i < xend) && is_eq; ++i) {
                for (j = ybeg; (j < yend); ++j) {
                    if (PPix(da,i,j) != tmp) {
                        is_eq = 0;
                        break;
                    }
                }
            }
            if (is_eq) {
                /* This is a constant array. */
                /* add NPIX1, NPIX2 (if necessary), and PIXVALUE keywords */
                kw = findKw(iodesc->hdr,"PIXVALUE");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"PIXVALUE",(int)tmp,
                        "values of pixels in constant array");
                else
                    putIntKw(kw,(int)tmp);

                kw = findKw(iodesc->hdr,"NPIX1");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"NPIX1",iodesc->dims[0],
                        "length of constant array axis 1");
                else
                    putIntKw(kw,iodesc->dims[0]);

                /* NPIX2 should only be added if the y-dimension is > 1. */
                if (da->tot_ny > 1) {
                    kw = findKw(iodesc->hdr,"NPIX2");
                    if (kw == 0) /* add it */
                        addIntKw(iodesc->hdr,"NPIX2",iodesc->dims[1],
                            "length of constant array axis 2");
                    else
                        putIntKw(kw,iodesc->dims[1]);
                }

                naxis = 0;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &naxis, NULL, &status);
                iodesc->dims[0] = 0;
                iodesc->dims[1] = 0;
                /* update the header, etc. */
                if (iodesc->hflag) {
                    iodesc->type = SHORT_IMG;
                    putHeader(iodesc);
                    iodesc->hflag = 0;
                }

                fits_flush_file(iodesc->ff, &status);

                clear_err();
                return 0;
            }
        }

        /* If not a constant array, make sure NPIX1, NPIX2, and PIXVALUE *
         * are NOT present in the header to be written out.              */
        kw = findKw(iodesc->hdr,"PIXVALUE");
        if (kw != 0) /* remove it */
            delKw(kw);

        kw = findKw(iodesc->hdr,"NPIX1");
        if (kw != 0) /* remove it */
            delKw(kw);

        if (da->tot_ny > 1) {
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw != 0) /* remove it */
                delKw(kw);
        }

        /* Get the current CFITSIO size, and if it's different, resize it */
        fits_get_img_size(iodesc->ff, 2, dims, &status);
        if (dims[0] != xend - xbeg || dims[1] != yend - ybeg) {
            iodesc->dims[0] = xend - xbeg;
            iodesc->dims[1] = yend - ybeg;
            if (fits_resize_img(iodesc->ff, SHORT_IMG, 2, iodesc->dims, &status)) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }

            /* Note, we don't need to fill the image with the constant value, since
               the image will be entirely over-written by the passed in array da */
        }

        /* update the header area */
        if (iodesc->hflag) {
            iodesc->type = SHORT_IMG;
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        fpixel[0] = 1;
        for (i = ybeg; i < yend; ++i) {
            fpixel[1] = i - ybeg + 1;
            if (fits_write_pix(iodesc->ff, TSHORT, fpixel, xsize,
                               (short *)&(PPix(da, xbeg, i)), &status)) {
                ioerr(BADWRITE,iodesc, status);
                return -1;
            }
        }

        fits_flush_file(iodesc->ff, &status);

        clear_err();
        return 0;
}

int getFloatLine(IODescPtr iodesc_, int line, float *ptr) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int no_dims, i, dim1;
        long dims[2];
        FitsKw kw;
        float val;
        long fpixel[2];
        int anynul;
        int status = 0;

        if (iodesc->options == WriteOnly) { ioerr(NOGET,iodesc,0); return -1; }

        if (fits_get_img_dim(iodesc->ff, &no_dims, &status)) {
            ioerr(BADDIMS, iodesc, status);
            return -1;
        }
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            dim1 = getIntKw(kw);
            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            val = getFloatKw(kw);
            for (i = 0; i < dim1; ++i) {
                ptr[i] = val;
            }
        } else {
            if (fits_get_img_size(iodesc->ff, 2, dims, &status)) {
                ioerr(BADDIMS, iodesc, status);
                return -1;
            }
            fpixel[0] = 1;
            fpixel[1] = line + 1;
            if (fits_read_pix(iodesc->ff, TFLOAT, fpixel, dims[0], NULL,
                              ptr, &anynul, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
        }
        clear_err();
        return 0;
}

int putFloatLine(IODescPtr iodesc_, int line, float *ptr) {
        IODesc *iodesc = (IODesc *)iodesc_;
        long fpixel[2];
        int no_dims;
        long dims[2];
        FitsKw kw;
        float* buffer;
        float val;
        long i, j;
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }
        if (iodesc->hflag) { iodesc->hflag = 0; putHeader(iodesc); }

        /* If a constant array, convert to a non-constant array */
        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) {
                ioerr(BADSCIDIMS, iodesc, 0);
                return -1;
            } else {
                dims[0] = getIntKw(kw);
                delKw(kw);
            }

            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw == 0) {
                ioerr(BADSCIDIMS, iodesc, 0);
                return -1;
            } else {
                dims[1] = getIntKw(kw);
                delKw(kw);
            }

            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) {
                ioerr(BADSCIDIMS,iodesc,0);
                return -1;
            } else {
                val = getFloatKw(kw);
                delKw(kw);
            }

            iodesc->dims[0] = dims[0];
            iodesc->dims[1] = dims[1];
            if (fits_resize_img(iodesc->ff, FLOAT_IMG, 2, dims, &status)) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }

            buffer = malloc(dims[0] * sizeof(float));
            if (buffer == NULL) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }

            for (i = 0; i < dims[0]; ++i) {
                buffer[i] = val;
            }

            /* Write the constant value into CFITSIO's array */
            fpixel[0] = 1;
            for (j = 0; j < dims[1]; ++j) {
                fpixel[1] = j + 1;
                if (fits_write_pix(iodesc->ff, TFLOAT, fpixel, dims[0], buffer, &status)) {
                    ioerr(BADWRITE, iodesc, status);
                    free(buffer);
                    return -1;
                }
            }

            free(buffer);
        }

        fpixel[0] = 1;
        fpixel[1] = line + 1;
        if (fits_write_pix(iodesc->ff, TFLOAT, fpixel, iodesc->dims[0],
                           ptr, &status)) {
            ioerr(BADWRITE, iodesc, status);
            return -1;
        }

        clear_err();
        return 0;
}

int getShortLine(IODescPtr iodesc_, int line, short *ptr) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int no_dims, dim1, i;
        long dims[2];
        FitsKw kw;
        short val;
        long fpixel[2];
        int anynul;
        int status = 0;

        if (iodesc->options == WriteOnly) { ioerr(NOGET,iodesc,0); return -1; }

        if (fits_get_img_dim(iodesc->ff, &no_dims, &status)) {
            ioerr(BADDIMS, iodesc, status);
            return -1;
        }
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            dim1 = getIntKw(kw);
            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            val = getIntKw(kw);
            for (i = 0; i < dim1; ++i)
                ptr[i] = val;
        } else {
            if (fits_get_img_size(iodesc->ff, 2, dims, &status)) {
                ioerr(BADDIMS, iodesc, status);
                return -1;
            }
            fpixel[0] = 1;
            fpixel[1] = line + 1;
            if (fits_read_pix(iodesc->ff, TSHORT, fpixel, dims[0], NULL,
                              ptr, &anynul, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
        }
        clear_err();
        return 0;
}

int putShortLine(IODescPtr iodesc_, int line, short *ptr) {
        IODesc *iodesc = (IODesc *)iodesc_;
        long fpixel[2];
        long dims[2];
        int no_dims;
        FitsKw kw;
        short val;
        short* buffer;
        long i, j;
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }
        if (iodesc->hflag) { iodesc->hflag = 0; putHeader(iodesc); }

        /* If a constant array, convert to a non-constant array */
        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) {
                ioerr(BADSCIDIMS, iodesc, 0);
                return -1;
            } else {
                dims[0] = getIntKw(kw);
                delKw(kw);
            }

            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw == 0) {
                ioerr(BADSCIDIMS, iodesc, 0);
                return -1;
            } else {
                dims[1] = getIntKw(kw);
                delKw(kw);
            }

            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) {
                ioerr(BADSCIDIMS,iodesc,0);
                return -1;
            } else {
                val = getIntKw(kw);
                delKw(kw);
            }

            iodesc->dims[0] = dims[0];
            iodesc->dims[1] = dims[1];
            if (fits_resize_img(iodesc->ff, SHORT_IMG, 2, dims, &status)) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }

            buffer = malloc(dims[0] * sizeof(short));
            if (buffer == NULL) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }

            for (i = 0; i < dims[0]; ++i) {
                buffer[i] = val;
            }

            /* Write the constant value into CFITSIO's array */

            fpixel[0] = 1;
            for (j = 0; j < dims[1]; ++j) {
                fpixel[1] = j + 1;
                if (fits_write_pix(iodesc->ff, TSHORT, fpixel, dims[0], buffer, &status)) {
                    ioerr(BADWRITE, iodesc, status);
                    free(buffer);
                    return -1;
                }
            }

            free(buffer);
        }

        fpixel[0] = 1;
        fpixel[1] = line + 1;
        if (fits_write_pix(iodesc->ff, TSHORT, fpixel, iodesc->dims[0],
                           ptr, &status)) {
            ioerr(BADWRITE, iodesc, status);
            return -1;
        }

        clear_err();
        return 0;
}
