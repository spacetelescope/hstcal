#ifndef INCL_CTABLES_H
#define INCL_CTABLES_H

/* NOTE:  this file (ctables.h) depends on fitsio.h */

# include "xtables.h"

# define SZ_ERRMESS  1024
# define SZ_FITS_STR   80

/* Descriptive information for a table. */
typedef struct {
        char *tablename;        /* table name as specified by the user */
        char *fullname;         /* OS file name (may still include brackets) */
        char *filename;         /* OS file name (without brackets) */
        char *brackets;         /* expression (if any) in brackets */
        fitsfile *fptr;         /* CFITSIO descriptor */
        IRAFPointer template;   /* tp for template, or NULL */
        int table_exists;       /* 1 if the table has been created */
        int iomode;             /* IRAF code for I/O mode */
        int hdunum;             /* HDU number (CFITSIO convention) */
        int hdutype;            /* type of HDU (CFITSIO convention) */
        long nrows;             /* number of rows in the table */
        int ncols;              /* number of columns */
        int alloc_cols;         /* current size of 'columns' array */
        IRAFPointer *columns;   /* array of descriptors for columns */
} TableDescr;

/* Descriptive information for a column.
   For a column containing strings, repeat = nelem * width;
   for columns of other data types, repeat = nelem.
*/
typedef struct {
        IRAFPointer tp;         /* tp for table that contains this column */
        int colnum;             /* column number (one indexed) */
        char *name;             /* column name */
        char *tform;            /* tform (data type as a string) */
        char *tunit;            /* units for column */
        char *tdisp;            /* display format */
        int typecode;           /* data type of column (CFITSIO typecode) */
        int datatype;           /* data type of column (IRAF code) */
        int var_length;         /* 1 if variable-length array */
        long repeat;            /* repeat count (r in rAw for strings) */
        long nelem;             /* number of elements in array */
        int width;              /* for a string, size of one element */
} ColumnDescr;

IRAFPointer init_tp (void);
IRAFPointer init_cp (IRAFPointer tp);
void free_tp (IRAFPointer tp);
void columnSpace (TableDescr *tbl_descr, int newcols);
void initCol (IRAFPointer tp, IRAFPointer *cp,
        char *colname, char *colunits, char *colfmt, int datatype, int nelem);
void addCol (IRAFPointer tp, IRAFPointer cp, char *colname, char *colunits);
void tbSaveInfo (IRAFPointer tp, int *status);
void tbCopyPrimary (fitsfile *template_fptr, fitsfile *fptr, int *status);
void tbCopyHeader (fitsfile *template_fptr, fitsfile *fptr, int *status);
void tbCopyTmpl (IRAFPointer tp);

/* in tbl_util.c */
char *expandfn (char *filename);
int checkExists (char *filename);
void cToFortran (char *c_fmt, char *ftn_fmt);
void trimString (char *value);
void copyString (char *output, char *input, int maxch);
void str_lower (char lc_name[], const char name[]);

/* in cerror.c */
void setError (int status, char *msg);
void clearError (void);
int checkError (void);

/* in c_vfn2osfn.c */
int c_vfn2osfn(const char *const vfn, char *const osfn);

#endif /* INCL_CTABLES_H */
