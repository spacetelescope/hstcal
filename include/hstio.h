# if !defined(HSTIO_)
# define HSTIO_

#if defined(__cplusplus)
extern "C" {
# endif

# include <stdlib.h>

# define HSTIO_VERSION "HSTIO Version 2.6 (11-Mar-2010)"

/*
** Data Structures and I/O Function Declarations for
** STScI Calibration Pipeline Software for STIS and NICMOS
**
** Version 2.1
**
** 23 March 1998
**
** Table of Contents
**      Basic structures to represent two-dimensional arrays.
**      Definitions of Science, Data Quality and Error data.
**      Macros to implement 2-d array indexing.
**      Definition of data sections.
**      Definition of Header "card" image arrays.
**      Definition of I/O descriptor.
**      Definitions of Single and Multi-Group data structures.
**      Declarations of error handling functions.
**      Declarations of high-level I/O functions.
**      Declarations of header manipulation functions.
**      Declarations of low-level I/O functions.
**      Declarations of functions to initialize, allocate and free
**              data storage contained within data structures.
**
** History
** Allen Farris - Original Implementation.
** Version 1.00 (06/21/95).
** Revisions to version 1.00 to produce version 1.01 (08/03/95).
**      Corrected bug in definition of IHdr.
**      Changed declaration of putMultiGroup function.
**      Changed the ordering within structures to Sci, Err, Dq.
** Revisions to version 1.01 to produce version 1.10 (09/28/95).
**      Corrected bug in definition of OUT_OPTION.
**      Added ifdef extern C to accomodate C++.
**      Added support for NICMOS samples and integration time data.
**      Made IODesc structure private by changing to void pointer
**          and changed the definitions of the IO_OPTIONS and added
**          functions to retrieve data from the IODesc structure.
**      Added open and close FITS file functions.
**      Major revision to functions to manipulate header arrays.
**      Added additional error handling functions
**      Revised definition of single and multi-group structures
**          to achieve greater compatibility between the two
** Revisions to version 1.10 to produce version 1.11 (09/29/95).
**      Added an option parameter to the putGroupHdr functions.
** Revisions to version 1.11 to produce version 1.20 (11/29/95).
**      Major revisions to the low-level I/O function and to the
**          low-level data structures.  References to Sci, Err, and
**          DQ were removed from the low-level structure and function
**          names in favor of generic names such as float and short.
** M.D. De La Pena  25 February 1998: Modified use of "long" to "int" to
** enforce compatibility with SPP/IRAF.
**
** Version 2.0, 23 March 1998:
** M.D. De La Pena: Added structures SingleGroupLine, FloatHdrLine, and
**      ShortHdrLine.  Also, new functions to support the acquisition of
**      obtaining single lines from a SingleGroup.  Moved the definition
**      of Bool_.
**
** 07 April 1998:
** M.D. De La Pena: Added functionality to output a subsection of each image
** of an imset to a file where the subsection is the full size (naxis1/naxis2)
** of the output image.
**
** 30 September 1998:
** M.D. De La Pena: Removed EXTVER parameter from getSingleGroupLine since
** it was unnecessary.
**
** 07 October 1998
** M.D. De La Pena: Added the function prototype for updateWCS.  This is a
** low-level support routine for putSect[Float/Short]HD.
**
** 09 November 1998:
** M.D. De La Pena: Retired the old get/put[Float/Short]Sect routines; these
** routines were meant to be used with the [Float/Short]TwoDArray data
** structure, but were never implemented as intended.  The following routine
** names have been updated to append "Sect" to the end of the name rather
** than using it as a prefix - consistency request: put[Sci/Err/DQ]Sect,
** put[Float/Short]HDSect, putSingleGroupSect, and put[Float/Short]Sect.
**
** 12 November 1998:
** M.D. De La Pena: Incorporated high-level keyword access routines written
** by H. Bushouse.   Modified the HSTIO_VERSION macro.
**
** 20 May 1999:
** R.L. Williamson: Modify prototype of updateWCS from type int to type void.
**
** 27 May 1999:
** M.D. De La Pena: Updated library version to 2.2 for RLW change.
**
** 18 August 1999:
** M.D. De La Pena: Updated library version to 2.3.  Added void to functions
** which have no parameters.
**
** 24 February 2000:
** Added comment to SingleGroupLine structure, line_num attribute is a
** zero-based reference system value.  Updated library version to 2.4.
**
** 07 May 2000:
** Added full support to accommodate reading/writing one-dimensional images.
** Updated library to 2.5. Note: Version 2.4 was implemented, but called 2.3.
**
** 16 Feb 2007:
** H.A. Bushouse: Added putSingleNicmosGroupSect, putSmplSect, and
** putIntgSect definitions.
**
** M. Droettboom, January 2010:
** Change to use CFITSIO rather than IRAF IMIO routines.
**
*/

/* Below is a general use pointer register (book keeping) to allow easy cleanup upon exit (pseudo C++ destructor)
 * Use example:
 *
 * PtrRegister ptrReg;
 * PtrRegister initPtrRegister(&ptrReg); // must always initialize before further use
 *
 * void * array = malloc(size);
 * addPtr(&ptrReg, array, &free); // returns if array == NULL so don't need to pre check - this is also self expanding
 *
 * Bool cancel = False;
 * ...
 * do something
 * cancel = True; //uh oh, something went wrong
 * ...
 * if (cancel)
 * {
 *     freeOnExit(&ptrReg); // This frees itself, i.e. the register array holding the ptrs
 *     return;
 * }
 * ...
 *
 *NOTE: This pattern is considered integral to all use and as such internal failed allocations are asserted
 */

typedef void (*FreeFunction)(void*); // Only trivial functions accepted
#define PTR_REGISTER_LENGTH_INC 10
typedef struct {
    unsigned cursor;
    unsigned length;
    void ** ptrs;
    FreeFunction * freeFunctions;
} PtrRegister;
void initPtrRegister(PtrRegister * reg);
void addPtr(PtrRegister * reg, void * ptr, void * freeFunc); // ptr list is self expanding
void freePtr(PtrRegister * reg, void * ptr);
void freeOnExit(PtrRegister * reg); //only calls freeAll() followed by freeReg()
void freeAll(PtrRegister * reg); // frees all ptrs registered (excluding itself)
void freeReg(PtrRegister * reg); // frees itself i.e. the register array holding the ptrs

# define SZ_PATHNAME 511

# if !defined(BOOL_)
# define BOOL_
enum Bool_ { False = 0, True = 1 };
typedef enum Bool_ Bool;
# endif

enum StorageOrder
{
    ROWMAJOR,
    COLUMNMAJOR
};

#define SCIEXT 0x01
#define ERREXT 0x02
#define DQEXT 0x04

typedef struct {
        short *buffer;          /* the pointer to the beg. of the buffer */
        int buffer_size;        /* the size of the full 2-d array in the */
        int tot_nx;             /* buffer.                               */
        int tot_ny;             /*     buffer_size = tot_nx*tot_ny       */
        int nx;                 /* The size of the current "view" of the */
        int ny;                 /* full 2-d array.                       */
        enum StorageOrder storageOrder;
        short *data;            /* The pointer to the beginning of the   */
                                /* subsection of the full 2-d array.     */
} ShortTwoDArray;

typedef struct {
        float *buffer;
        int buffer_size;
        int tot_nx;
        int tot_ny;
        int nx;
        int ny;
        enum StorageOrder storageOrder;
        float *data;
} FloatTwoDArray;


/*
** The Following macros are used to represent 2-d indexing.
** Two dimensional arrays are stored in FITS order.
**
**        ny
**        ^
**      N | a05   a15   a25   a35
**      A | a04   a14   a24   a34
**      X | a03   a13   a23   a33
**      I | a02   a12   a22   a32
**      S | a01   a11   a21   a31
**      2 | a00   a10   a20   a30
**         ---------------------------> nx
**            NAXIS1
**
**      NAXIS1 is 4 and NAXIS2 is 6
**      PIX(a,1,4) accesses a14
*/

# define Pix(a,i,j)      (a).data[(j)*(a).tot_nx + (i)]
# define PixColumnMajor(a,i,j) (a).data[(j)*(a).tot_ny + (i)]
# define DQPix(a,i,j)    (a).data[(j)*(a).tot_nx + (i)]
# define DQSetPix(a,i,j,value) (a).data[(j)*(a).tot_nx + (i)] = (value)

# define PPix(a,i,j)      (a)->data[(j)*(a)->tot_nx + (i)]
# define PPixColumnMajor(a,i,j) (a)->data[(j)*(a)->tot_ny + (i)]
# define PDQPix(a,i,j)    (a)->data[(j)*(a)->tot_nx + (i)]
# define PDQSetPix(a,i,j,value) (a)->data[(j)*(a)->tot_nx + (i)] = (value)

/*
**              Examples of using macros.
**
** Suppose we have the following declarations and that these data
** structures have been properly initialized and allocated.
**
**      SciData picture;
**      DQData dataqual_in1, dataqual_in2, dataqual_out;
**      short val;
**
** Now set all the diagonal pixels in the picture array to 1.
**
**      for (i = 0; i < min(picture.nx,picture.ny); i++)
**              Pix(picture,i,i) = 1.0F;
**
** Now combine the two data quality arrays.
**
**      for ( j = 0; j < dataqual_out.ny; j++)
**          for (i = 0; i < dataqual_out.nx; i++) {
**              val = DQPix(dataqual_in1,i,j) | DQPix(dataqual_in2,i,j);
**              DQSetPix(dataqual_out,i,j,val);
**          }
**
*/

/*
** The data section structure can be used to read sections of an image
** from disk into memory.  It is intended to be used only for high-level
** I/O operations.  At present, the entire array is read into memory.
*/
typedef struct {
        int x_beg;              /* The beginning coordinates of the     */
        int y_beg;              /* section.                             */
        int sx;                 /* The sizes of the X and Y dimensions  */
        int sy;                 /* of the section.                      */
} DataSection;
void copyDataSection(DataSection * dest, const DataSection * src);


# define HDRSize 81
typedef char HdrArray[HDRSize]; /* Headers are simply an array of fixed */
                                /* length, null terminated strings that */
                                /* represent FITS card images.          */

typedef struct {
        int nlines;             /* The number of lines actually used.   */
        int nalloc;             /* Number of lines currently allocated. */
        HdrArray *array;        /* The buffer of card images.           */
} Hdr;

/*
** I/O Options and the I/O Descriptor Pointer
**
** The I/O descriptor pointer points to an I/O structure that is
** used internally within the I/O functions.
*/
# if !defined(IO_OPTIONS_)
# define IO_OPTIONS_
enum IO_OPTION_ {
        Sequential = 0x1,
        Random = 0x2,
        ReadOnly = 0x4,
        WriteOnly = 0x8,
        ReadWrite = 0x10,
        New = 0x20,
        Overwrite = 0x40,
        Force = 0x80,
        Inherit = 0x100,
        Dupname = 0x200,
        NoDupname = 0x400,
        RemainOpen = 0x800,
        Geis = 0x1000,
        Tape9 = 0x2000,
        Std = 0x4000,
        Disk = 0x8000
};
typedef enum IO_OPTION_ OP_OPTION;
# endif
typedef void *IODescPtr;

/*
** The following data structures define combinations of both
** headers and data.
*/
typedef struct {
        IODescPtr iodesc;
        DataSection section;
        Hdr hdr;
        FloatTwoDArray data;
} FloatHdrData;

typedef struct {
        IODescPtr iodesc;
        DataSection section;
        Hdr hdr;
        ShortTwoDArray data;
} ShortHdrData;

typedef FloatHdrData SciHdrData;
typedef FloatHdrData ErrHdrData;
typedef ShortHdrData DQHdrData;
typedef ShortHdrData SmplHdrData;
typedef FloatHdrData IntgHdrData;

/*
** Data structures which define combinations of both headers and
** a single line of data.
*/

typedef struct {
        IODescPtr iodesc;         /* file descriptor                   */
        Hdr       hdr;            /* header structure for extension    */
        Bool      ehdr_loaded;    /* flag, is extension header loaded? */
        int       tot_nx;         /* number of values in line          */
        float     *line;          /* pointer to line buffer            */
} FloatHdrLine;

typedef struct {
        IODescPtr iodesc;
        Hdr       hdr;
        Bool      ehdr_loaded;
        int       tot_nx;
        short     *line;
} ShortHdrLine;

typedef FloatHdrLine SciHdrLine;
typedef FloatHdrLine ErrHdrLine;
typedef ShortHdrLine DQHdrLine;

/*
** The SingleGroup data structure is used for single images or single
** groups from a multi-group file.  The disk file name together with
** the group number uniquely identifies the source of the data array
** in the file.  The group number is the FITS EXTVER keyword value.
*/
/*\group(smgroup)*/
typedef struct {
        char *filename;
        int group_num;
        Hdr *globalhdr;
        SciHdrData sci;
        ErrHdrData err;
        DQHdrData dq;
} SingleGroup;

typedef struct {
        char *filename;
        int group_num;
        Hdr *globalhdr;
        SciHdrData sci;
        ErrHdrData err;
        DQHdrData dq;
        SmplHdrData smpl;
        IntgHdrData intg;
} SingleNicmosGroup;
/*\endgroup*/

/*
** The SingleGroupLine data structure is used to access lines from a single
** group from a multi-group file.
*/
typedef struct {
        char *filename;         /* filename                               */
        int  group_num;         /* EXTVER or group number                 */
        int  line_num;          /* current line number loaded, zero-based */
        Bool phdr_loaded;       /* flag, is primary header loaded?        */
        Hdr  *globalhdr;        /* header structure for primary           */
        SciHdrLine sci;         /* science line data structure            */
        ErrHdrLine err;         /* error line data structure              */
        DQHdrLine dq;           /* dq line data structure                 */
} SingleGroupLine;

/*
** The MultiGroup data structure is used to read members of a multi-group
** file.
*/
typedef struct {
        int ngroups;            /* number of groups contained in the struct  */
        SingleGroup *group;
} MultiGroup;

typedef struct {
        int ngroups;            /* number of groups contained in the struct  */
        SingleNicmosGroup *group;
} MultiNicmosGroup;

/*
** Error Handling Function Declarations
*/
enum HSTIOERROR_ {
        HSTOK, NOMEM, BADOPEN, BADCLOSE, BADREAD, BADWRITE, BADEXTNAME,
        BADHSIZE, NOGET, NOPUT, BADDIMS, BADTYPE, NOSCI, BADSCIDIMS,
        BADGROUP, BADGET, BADFITSEQ, BADFITSQUOTE, BADFITSNUMERIC,
        BADFITSTYPE, BADPUT, BADNAME, BADBITPIX, BADNDIM, BADEXIST,
        BADREMOVE
};
typedef enum HSTIOERROR_ HSTIOError;
HSTIOError hstio_err(void);
char *hstio_errmsg(void);
typedef void (*HSTIOErrHandler)(void);
int push_hstioerr(HSTIOErrHandler);
int pop_hstioerr(void);
void clear_hstioerr(void);

/*
** High-level I/O Function Declarations
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
int  openSingleGroupLine  (char *filename, int extver, SingleGroupLine *);
void closeSingleGroupLine (SingleGroupLine *);

int ckNewFile(char *fname);
int openFitsFile(char *filename, unsigned int option);
int closeFitsFile(char *filename);
int getSci(char *filename, int extver, SciHdrData *);
int putSci(char *filename, int extver, SciHdrData *, int option);
int getErr(char *filename, int extver, ErrHdrData *);
int putErr(char *filename, int extver, ErrHdrData *, int option);
int getDQ(char *filename, int extver, DQHdrData *);
int putDQ(char *filename, int extver, DQHdrData *, int option);
int getSmpl(char *filename, int extver, SmplHdrData *);
int putSmpl(char *filename, int extver, SmplHdrData *, int option);
int getIntg(char *filename, int extver, IntgHdrData *);
int putIntg(char *filename, int extver, IntgHdrData *, int option);

int putSciSect(char *filename, int extver, SciHdrData *, int xbeg, int ybeg,
               int xsize, int ysize, int option);
int putErrSect(char *filename, int extver, ErrHdrData *, int xbeg, int ybeg,
               int xsize, int ysize, int option);
int putDQSect(char *filename, int extver, DQHdrData *, int xbeg, int ybeg,
               int xsize, int ysize, int option);
int putSmplSect(char *filename, int extver, SmplHdrData *, int xbeg, int ybeg,
               int xsize, int ysize, int option);
int putIntgSect(char *filename, int extver, IntgHdrData *, int xbeg, int ybeg,
               int xsize, int ysize, int option);

int getSciHdr(char *filename, int extver, SciHdrLine *);
int getErrHdr(char *filename, int extver, ErrHdrLine *);
int getDQHdr(char *filename, int extver, DQHdrLine *);
int getSciLine (SciHdrLine *x, int);
int getErrLine (ErrHdrLine *x, int);
int getDQLine  (DQHdrLine *x, int);
int getFloatHD(char *filename, char *extname, int extver, FloatHdrData *);
int putFloatHD(char *filename, char *extname, int extver, FloatHdrData *, int option);
int getShortHD(char *filename, char *extname, int extver, ShortHdrData *);
int putShortHD(char *filename, char *extname, int extver, ShortHdrData *, int option);

int putFloatHDSect(char *filename, char *extname, int extver, FloatHdrData *,
                   int xbeg, int ybeg, int xsize, int ysize, int option);
int putShortHDSect(char *filename, char *extname, int extver, ShortHdrData *,
                   int xbeg, int ybeg, int xsize, int ysize, int option);

int getFloatHdr(char *filename, char *extname, int extver, FloatHdrLine *);
int getShortHdr(char *filename, char *extname, int extver, ShortHdrLine *);
int getSingleGroup(char *filename, int extver, SingleGroup *);
int getSingleGroupLine(char *filename, int line, SingleGroupLine *);
int putSingleGroupHdr(char *filename, SingleGroup *, int option);
int putSingleGroup(char *filename, int extver, SingleGroup *, int option);

int putSingleGroupSect(char *filename, int extver, SingleGroup *,
                       int xbeg, int ybeg, int xsize, int ysize, int option);

int getSingleNicmosGroup(char *filename, int extver, SingleNicmosGroup *);
int putSingleNicmosGroupHdr(char *filename, SingleNicmosGroup *, int option);
int putSingleNicmosGroup(char *filename, int extver, SingleNicmosGroup *,
        int option);
int putSingleNicmosGroupSect(char *filename, int extver, SingleNicmosGroup *,
        int xbeg, int ybeg, int xsize, int ysize, int option);
int getMultiGroupHdr(char *filename, MultiGroup *);
int getMultiGroup(MultiGroup *, int ngroup, int extver);
int putMultiGroupHdr(char *filename, MultiGroup *, int option);
int putMultiGroup(char *filename, int extver, MultiGroup *, int ngroup,
        int option);
int getMultiNicmosGroupHdr(char *filename, MultiNicmosGroup *);
int getMultiNicmosGroup(MultiNicmosGroup *, int ngroup, int extver);
int putMultiNicmosGroupHdr(char *filename, MultiNicmosGroup *, int option);
int putMultiNicmosGroup(char *filename, int extver, MultiNicmosGroup *,
        int ngroup, int option);

/*
** Functions for Accessing and Manipulating Keywords in a FITS Keyword List
*/

# define NotFound NULL
typedef void *FitsKw;
enum FitsDataType_ { FITSNOVALUE = 0, FITSLOGICAL, FITSBIT, FITSCHAR,
        FITSBYTE, FITSSHORT, FITSLONG, FITSFLOAT, FITSDOUBLE, FITSCOMPLEX,
        FITSICOMPLEX, FITSDCOMPLEX, FITSVADESC
};
typedef enum FitsDataType_ FitsDataType;

/*
** Making a FITS Header
**
** These return 0 is successful and -1 if an error occurred.
*/
/*\group(makehdr)*/
int makePrimaryArrayHdr(Hdr *, FitsDataType, int dims, int *ndim);
int makeImageExtHdr(Hdr *, FitsDataType, int dims, int *ndim,
        char *extname, int extver);
/*\endgroup*/

/*
** Finding keywords in a FITS keyword list.
**
** The "find" function starts at the beginning of the list and finds the
** first occurence of the desired keyword.  The "findnext" function
** begins at the current position and returns the next occurence of the
** keyword.  NotFound is returned if no such keyword is found.  The "first"
** function sets the current position to the beginning of the list.
** The "next" function sets the position to the next keyword in the list.
** "Next" returns NotFound after the last keyword in the list.  "GetKw"
** returns the n-th keyword.  The special function "insertfirst" is used to
** position the cursor prior to the first item in the list so that a new
** keyword may be inserted at the beginning of the list.  (Insert functions
** always insert keywords at a position after the cursor.)
*/
/*\group(find)*/
FitsKw findKw(Hdr *, char *name);
FitsKw findnextKw(Hdr *, char *name);
FitsKw first(Hdr *);
FitsKw next(FitsKw);
FitsKw getKw(Hdr *, int n);
FitsKw insertfirst(Hdr *);
/*\endgroup*/
/*
** Accessing and changing name, value, or comment in a FITS keyword.
**
** These functions must be preceded by a "find" or "findnext" operation.
** Data conversions in the "get" routines:
**              FITS Type       Possible Requested Type
**              ---------       -----------------------
**              Logical         Bool, Int
**              Int             Int, Double
**              Float           Float, Double
**              Double          Double
**              String          String
** If the requested type conversion is not allowed, an error code is set.
*/
/*\group(get)*/
/* High-level keyword access routines */
int getKeyB (Hdr *, char *keyword, Bool *value);
int getKeyI (Hdr *, char *keyword, int *value);
int getKeyF (Hdr *, char *keyword, float *value);
int getKeyD (Hdr *, char *keyword, double *value);
int getKeyS (Hdr *, char *keyword, char *value);
int putKeyB (Hdr *, char *keyword, Bool value, char *comment);
int putKeyI (Hdr *, char *keyword, int value, char *comment);
int putKeyF (Hdr *, char *keyword, float value, char *comment);
int putKeyD (Hdr *, char *keyword, double value, char *comment);
int putKeyS (Hdr *, char *keyword, char *value, char *comment);
/* End high-level keyword access routines */

char *getKwName(FitsKw);
char *getKwComm(FitsKw);
FitsDataType getKwType(FitsKw);
Bool getBoolKw(FitsKw);
int  getIntKw(FitsKw);
float getFloatKw(FitsKw);
double getDoubleKw(FitsKw);
int getStringKw(FitsKw, char *str, int maxch);
int putKwName(FitsKw, char *name);
void putKwComm(FitsKw, char *comment);
int putBoolKw(FitsKw, Bool value);
int putIntKw(FitsKw, int value);
int putFloatKw(FitsKw, float value);
int putDoubleKw(FitsKw, double value);
int putStringKw(FitsKw, char *value);
/*\endgroup*/

/*
** Adding keywords to a FITS Header.
**
** The new keywords are added to the end of the keyword list.  There
** is no "END" keyword maintained in the list.  The software writing
** the FITS header automatically adds an "END" keyword before writing
** the header.
**
*/
/*\group(add)*/
int addBoolKw(Hdr *, char *name, Bool value, char *comment);
int addIntKw(Hdr *, char *name, int value, char *comment);
int addFloatKw(Hdr *, char *name, float value, char *comment);
int addDoubleKw(Hdr *, char *name, double value, char *comment);
int addStringKw(Hdr *, char *name, char * value, char *comment);
int addSpacesKw(Hdr *, char *comment);
int addCommentKw(Hdr *, char *comment);
int addHistoryKw(Hdr *, char *comment);
/*\endgroup*/

/*
** Inserting keywords in a FITS Header.
**
** The keywords are inserted after the current position.  The position
** is set with the "find", "first" and "next" functions.
**
*/
/*\group(insert)*/
FitsKw insertBoolKw(FitsKw, char *name, Bool value, char *comment);
FitsKw insertIntKw(FitsKw, char *name, int value, char *comment);
FitsKw insertFloatKw(FitsKw, char *name, float value, char *comment);
FitsKw insertDoubleKw(FitsKw, char *name, double value, char *comment);
FitsKw insertStringKw(FitsKw, char *name, char * value, char *comment);
FitsKw insertSpacesKw(FitsKw, char *comment);
FitsKw insertCommentKw(FitsKw, char *comment);
FitsKw insertHistoryKw(FitsKw, char *comment);
/*\endgroup*/

/*
** Adding and inserting already formatted FITS cards to the list
**
** These functions take an already formatted FITS "card" and either adds
** the string to the list or inserts it after the current position.  If
** the string is less than 80 bytes it is padded with spaces to 80 bytes.
** If the string is more than 80 bytes it is truncated.
*/
/*\group(allformat)*/
int addFitsCard(Hdr *, char *card);
FitsKw insertFitsCard(FitsKw, char *card);
/*\endgroup*/

/*
** Deleting keywords in a FITS Header.
**
** The "del" function must be preceeded by a "find" or "findnext" operation.
** The "delAllKw" function deletes all keywords.
*/
/*\group(del)*/
void delKw(FitsKw);
void delAllKw(Hdr *);
/*\endgroup*/

/*
** Low-level I/O Function Declarations
*/
IODescPtr openInputImage(char *filename, char *extname, int extver);
IODescPtr openOutputImage(char *filename, char *extname, int extver, Hdr *hdr,
        int dim1, int dim2, FitsDataType type);
IODescPtr openUpdateImage(char *filename, char *extname, int extver, Hdr *hdr);
void closeImage(IODescPtr );

char *getFilename(IODescPtr);
char *getExtname(IODescPtr);
int getExtver(IODescPtr);
int getNaxis1(IODescPtr);
int getNaxis2(IODescPtr);
int getType(IODescPtr);

int getHeader(IODescPtr, Hdr *);
int putHeader(IODescPtr);

int getFloatData(IODescPtr, FloatTwoDArray *);
int putFloatData(IODescPtr, FloatTwoDArray *);
int getShortData(IODescPtr, ShortTwoDArray *);
int putShortData(IODescPtr, ShortTwoDArray *);

int putFloatSect(IODescPtr, FloatTwoDArray *, int, int, int, int);
int putShortSect(IODescPtr, ShortTwoDArray *, int, int, int, int);

int getFloatLine(IODescPtr, int line, float *);
int putFloatLine(IODescPtr, int line, float *);
int getShortLine(IODescPtr, int line, short *);
int putShortLine(IODescPtr, int line, short *);
/*
** Low-level Support Function Declarations
*/
void updateWCS (Hdr *hdr, int xbeg, int ybeg);
/*
** Initialization, Allocation, and Freeing Storage Function Declarations
*/
# define IFloatData { NULL, 0, 0, 0, 0, 0, NULL }
void initFloatData(FloatTwoDArray *);
int allocFloatData(FloatTwoDArray *, int, int, Bool zeroInitialize);
void freeFloatData(FloatTwoDArray *);
int swapFloatStorageOrder(FloatTwoDArray * target, const FloatTwoDArray * source, enum StorageOrder targetStorageOrder);
# define IShortData { NULL, 0, 0, 0, 0, 0, NULL }
void initShortData(ShortTwoDArray *);
int allocShortData(ShortTwoDArray *, int, int, Bool zeroInitialize);
void freeShortData(ShortTwoDArray *);
int swapShortStorageOrder(ShortTwoDArray * target, const ShortTwoDArray * source, enum StorageOrder targetStorageOrder);

# define IFloatLine { NULL, 0 }
void initFloatLine  (FloatHdrLine *);
int  allocFloatLine (FloatHdrLine *, int );
void freeFloatLine  (FloatHdrLine *);
# define IShortLine { NULL, 0 }
void initShortLine  (ShortHdrLine *);
int  allocShortLine (ShortHdrLine *, int );
void freeShortLine  (ShortHdrLine *);

int  allocSciLine (SingleGroupLine *);
int  allocErrLine (SingleGroupLine *);
int  allocDQLine  (SingleGroupLine *);

# define IHdr { 0, 0, NULL }
void initHdr(Hdr *);
int allocHdr(Hdr *, int, Bool zeroInitialize);
int reallocHdr(Hdr *, int);
void freeHdr(Hdr *);
int copyHdr(Hdr *to, const Hdr *from);

# define IFloatHdrData { NULL, { 0, 0, 0, 0 }, IHdr, IFloatData }
void initFloatHdrData(FloatHdrData *);
int allocFloatHdrData(FloatHdrData *, int, int, Bool zeroInitialize);
int copyFloatHdrData(FloatHdrData * target, const FloatHdrData * src, enum StorageOrder targetStorageOrder);
void freeFloatHdrData(FloatHdrData *);
# define IShortHdrData { NULL, { 0, 0, 0, 0 }, IHdr, IShortData }
void initShortHdrData(ShortHdrData *);
int allocShortHdrData(ShortHdrData *, int, int, Bool zeroInitialize);
int copyShortHdrData(ShortHdrData * target, const ShortHdrData * src, enum StorageOrder targetStorageOrder);
void freeShortHdrData(ShortHdrData *);

# define IFloatHdrLine { NULL, IHdr, False, 0, NULL }
void initFloatHdrLine  (FloatHdrLine *);
int  allocFloatHdrLine (FloatHdrLine *, int);
void freeFloatHdrLine  (FloatHdrLine *);
# define IShortHdrLine { NULL, IHdr, False, 0, NULL }
void initShortHdrLine  (ShortHdrLine *);
int  allocShortHdrLine (ShortHdrLine *, int);
void freeShortHdrLine  (ShortHdrLine *);

# define ISingleGroup { NULL, 0, NULL, IFloatHdrData, IFloatHdrData, \
IShortHdrData }
void initSingleGroup(SingleGroup *);
int allocSingleGroup(SingleGroup *, int, int, Bool zeroInitialize);
int allocSingleGroupHeader(Hdr ** hdr, Bool zeroInitialize);
int allocSingleGroupExts(SingleGroup *x, int i, int j, unsigned extension, Bool zeroInitialize);
void freeSingleGroup(SingleGroup *);
void setStorageOrder(SingleGroup * group, enum StorageOrder storageOrder);
int copySingleGroup(SingleGroup * target, const SingleGroup * source, enum StorageOrder targetStorageOrder);
void copyOffsetSingleGroup(SingleGroup * output, const SingleGroup * input, unsigned nRows, unsigned nColumns, unsigned outputOffset, unsigned inputOffset, unsigned outputSkipLength, unsigned inputSkipLength);
# define IMultiGroup { 0, NULL }
void initMultiGroup(MultiGroup *);
int allocMultiGroup(MultiGroup *, int);
void freeMultiGroup(MultiGroup *);

int copyFloatData(FloatTwoDArray * target, const FloatTwoDArray * source, enum StorageOrder targetStorageOrder);
void copyOffsetFloatData(float * output, const float * input, unsigned nRows, unsigned nColumns, unsigned outputOffset, unsigned inputOffset, unsigned outputSkipLength, unsigned inputSkipLength);
int copyShortData(ShortTwoDArray * target, const ShortTwoDArray * source, enum StorageOrder targetStorageOrder);
void copyOffsetShortData(short * output, const short * input, unsigned nRows, unsigned nColumns, unsigned outputOffset, unsigned inputOffset, unsigned outputSkipLength, unsigned inputSkipLength);

# define ISingleNicmosGroup { NULL, 0, NULL, IFloatHdrData, IFloatHdrData, \
IShortHdrData, IShortHdrData, IFloatHdrData }
void initSingleNicmosGroup(SingleNicmosGroup *);
int allocSingleNicmosGroup(SingleNicmosGroup *, int, int);
void freeSingleNicmosGroup(SingleNicmosGroup *);
# define IMultiNicmosGroup { 0, NULL }
void initMultiNicmosGroup(MultiNicmosGroup *);
int allocMultiNicmosGroup(MultiNicmosGroup *, int);
void freeMultiNicmosGroup(MultiNicmosGroup *);

# define ISingleGroupLine { NULL, 0, 0, False, IHdr, IFloatHdrLine, IFloatHdrLine, \
IShortHdrLine }
void initSingleGroupLine  (SingleGroupLine *);
int  allocSingleGroupLine (SingleGroupLine *, int);
void freeSingleGroupLine  (SingleGroupLine *);

#if defined(__cplusplus)
}
# endif

# endif
