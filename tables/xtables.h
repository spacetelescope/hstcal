#ifndef INCL_XTABLES_H
#define INCL_XTABLES_H

# include <c_iraf.h>

/* These constants do NOT include space for the NULL character */
# define SZ_KEYWORD 8
# define SZ_COLNAME 79
# define SZ_COLUNITS 79
# define SZ_COLFMT 79
# define SZ_PARREC 80

/* Indefinite valued numbers. Taken from IRAF's $hlib/iraf.h */
# define IRAF_INDEFS (-32767)
# define IRAF_INDEFL (0x80000001)
# define IRAF_INDEFI IRAF_INDEFL
# define IRAF_INDEFR 1.6e38
# if defined(__VMS)
# define IRAF_INDEFD 1.6e38
# else
# define IRAF_INDEFD 1.6e308
# endif

/* Error codes used in tables. */
# define ERR_OUT_OF_MEMORY             3
# define ERR_STRING_TOO_LONG           5
# define ERR_ARRAY_TOO_LARGE           7
# define ERR_PARAMETER_UNKNOWN         9
# define ERR_DATATYPE_UNKNOWN         11
# define ERR_NOT_A_TABLE              13
# define ERR_COLUMN_ALREADY_EXISTS    15

/* These may be read by c_tbpsta but may not be set: */
# define TBL_WHTYPE         5   /* which type of table?  (TBL_TYPE_FITS) */
# define TBL_NROWS         21   /* number of rows written to */
# define TBL_NCOLS         22   /* number of columns # defined */
# define TBL_NPAR          24   /* number of keywords */

/* these are table types (see TBL_WHTYPE) that may be relevant */
# define TBL_TYPE_TEXT     13   /* type is text file */
# define TBL_TYPE_FITS     14   /* type is FITS table */

/* These are for information about a column. */
# define TBL_COL_NAME      41   /* column name */
# define TBL_COL_UNITS     42   /* units for column */
# define TBL_COL_FMT       43   /* print format for displaying values */
# define TBL_COL_DATATYPE  44   /* data type (-n for char string) */
# define TBL_COL_NUMBER    45   /* column number */
# define TBL_COL_FMTLEN    46   /* length for printing using print fmt */
# define TBL_COL_LENDATA   47   /* number of elements if it's an array */
# define TBL_COL_DIMENSION 48   /* dimension of array */

/* These five are the image file name template functions. */
IRAFPointer c_imtopen (char *pattern);
int c_imtlen (IRAFPointer imt);
void c_imtrew (IRAFPointer imt);
int c_imtgetim (IRAFPointer imt, char *outstr, int maxch);
void c_imtclose (IRAFPointer imt);

IRAFPointer c_tbtopn (char *tablename, int iomode, IRAFPointer template);
void c_tbtcre (IRAFPointer tp);
void c_tbtclo (IRAFPointer tp);
int c_tbtacc (char *tablename);
void c_tbtnam (IRAFPointer tp, char *tablename, int maxch);
void c_tbfpri (char *intable, char *outtable, int *copied);
int c_tbparse (char *tablename, char *fname, char *extname, int maxch,
                int *hdu);

int c_tbpsta (IRAFPointer tp, int param);

void c_tbcdef1 (IRAFPointer tp, IRAFPointer *cp,
        char *colname, char *colunits, char *colfmt, int datatype, int nelem);
void c_tbcfnd1 (IRAFPointer tp, const char *colname, IRAFPointer *cp);
IRAFPointer c_tbcnum (IRAFPointer tp, int colnum);

int c_tbcigi (IRAFPointer cp, int param);
void c_tbcigt (IRAFPointer cp, int param, char *outstr, int maxch);
void c_tbciga (IRAFPointer tp, IRAFPointer cp, int *ndim, int *axlen,
                int maxdim);
void c_tbcinf (IRAFPointer cp, int *colnum, char *colname, char *colunits,
                char *colfmt, int *datatype, int *nelem, int *lenfmt);

int c_tbfres (char *keyword);
void c_tbhcal (IRAFPointer itp, IRAFPointer otp);
void c_tbhgnp (IRAFPointer tp, int parnum,
                char *keyword, int *dtype, char *str);
void c_tbhgcm (IRAFPointer tp, char *keyword, char *comment, int maxch);
void c_tbhpcm (IRAFPointer tp, char *keyword, char *comment);
Bool c_tbhgtb (IRAFPointer tp, char *keyword);
double c_tbhgtd (IRAFPointer tp, char *keyword);
float c_tbhgtr (IRAFPointer tp, char *keyword);
int c_tbhgti (IRAFPointer tp, char *keyword);
void c_tbhgtt (IRAFPointer tp, char *keyword, char *text, int maxch);
void c_tbhadb (IRAFPointer tp, char *keyword, Bool value);
void c_tbhadd (IRAFPointer tp, char *keyword, double value);
void c_tbhadr (IRAFPointer tp, char *keyword, float value);
void c_tbhadi (IRAFPointer tp, char *keyword, int value);
void c_tbhadt (IRAFPointer tp, char *keyword, char *text);
void c_tbhptt (IRAFPointer tp, char *keyword, char *text);

void c_tbrgtr (IRAFPointer tp, IRAFPointer *cp, float *buffer, Bool *nullflag,
                int numcols, int row);
void c_tbrcsc (IRAFPointer itp, IRAFPointer otp,
                IRAFPointer *icp, IRAFPointer *ocp,
                int irow, int orow, int ncols);
void c_tbrudf (IRAFPointer tp, IRAFPointer *cp, int numcols, int row);

void c_tbegtb (IRAFPointer tp, IRAFPointer cp, int rownum, Bool *buffer);
void c_tbegtd (IRAFPointer tp, IRAFPointer cp, int rownum, double *buffer);
void c_tbegtr (IRAFPointer tp, IRAFPointer cp, int rownum, float *buffer);
void c_tbegti (IRAFPointer tp, IRAFPointer cp, int rownum, int *buffer);
void c_tbegts (IRAFPointer tp, IRAFPointer cp, int rownum, short *buffer);
void c_tbegtt (IRAFPointer tp, IRAFPointer cp, int rownum, char *buffer,
                int maxch);
void c_tbeptb (IRAFPointer tp, IRAFPointer cp, int rownum, Bool buffer);
void c_tbeptd (IRAFPointer tp, IRAFPointer cp, int rownum, double buffer);
void c_tbeptr (IRAFPointer tp, IRAFPointer cp, int rownum, float buffer);
void c_tbepti (IRAFPointer tp, IRAFPointer cp, int rownum, int buffer);
void c_tbepts (IRAFPointer tp, IRAFPointer cp, int rownum, short buffer);
void c_tbeptt (IRAFPointer tp, IRAFPointer cp, int rownum, char *buffer);

int c_tbagtb (IRAFPointer tp, IRAFPointer cp, int row, Bool *buffer,
                int first, int nelem);
int c_tbagtd (IRAFPointer tp, IRAFPointer cp, int row, double *buffer,
                int first, int nelem);
int c_tbagtr (IRAFPointer tp, IRAFPointer cp, int row, float *buffer,
                int first, int nelem);
int c_tbagti (IRAFPointer tp, IRAFPointer cp, int row, int *buffer,
                int first, int nelem);
int c_tbagts (IRAFPointer tp, IRAFPointer cp, int row, short *buffer,
                int first, int nelem);
int c_tbagtt (IRAFPointer tp, IRAFPointer cp, int row, char **buffer,
                int first, int nelem, int maxch);
void c_tbaptb (IRAFPointer tp, IRAFPointer cp, int row, Bool *buffer,
                int first, int nelem);
void c_tbaptd (IRAFPointer tp, IRAFPointer cp, int row, double *buffer,
                int first, int nelem);
void c_tbaptr (IRAFPointer tp, IRAFPointer cp, int row, float *buffer,
                int first, int nelem);
void c_tbapti (IRAFPointer tp, IRAFPointer cp, int row, int *buffer,
                int first, int nelem);
void c_tbapts (IRAFPointer tp, IRAFPointer cp, int row, short *buffer,
                int first, int nelem);
void c_tbaptt (IRAFPointer tp, IRAFPointer cp, int row, char **cbuf,
                int maxch, int first, int nelem);

#endif /* INCL_XTABLES_H */
