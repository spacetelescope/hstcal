# include <fitsio.h>
# include "xtables.h"

# define SZ_FNAME    1025
# define SZ_ERRMESS  1025
# define SZ_FITS_STR   81

# define ERR_OUT_OF_MEMORY     5
# define ERR_DATATYPE_UNKNOWN  7
# define ERR_NOT_A_TABLE       9

# define TABLE_DESCR  1
# define COLUMN_DESCR 2

/* Descriptive information for a table. */
typedef struct {
	char *tablename;	/* table name as specified by the user */
	char *filename;		/* OS file name */
	char *brackets;		/* expression (if any) in brackets */
	fitsfile *fptr;		/* CFITSIO descriptor */
	int hdunum;		/* HDU number (CFITSIO convention) */
	int hdutype;		/* type of HDU */
	long nrows;		/* number of rows in the table */
	int ncols;		/* number of columns */
	int alloc_cols;		/* number of elements allocated for columns */
	IRAFPointer *columns;	/* array of descriptors for columns */
	int table_exists;	/* 1 if the table has been created */
} TableDescr;

/* Descriptive information for a column. */
typedef struct {
	int colnum;		/* column number, one indexed */
	char *name;		/* column name */
	char *tunit;		/* units for column */
	char *tdisp;		/* display format */
	int dtype;		/* data type of column */
	long repeat;		/* repeat count (number of elements) */
} ColumnDescr;

/* This is an element of the array indexed by IRAFPointer. */
typedef struct {
	int ptype;	/* indicates whether table or column is specified */
	TableDescr *t_descr;
	ColumnDescr *c_descr;
} TblInfoPtr;

IRAFPointer init_tp (void);
IRAFPointer init_cp (void);
void free_tp (IRAFPointer tp);
void columnSpace (TableDescr *t_descr, int newcols);
TblInfoPtr *makeTableDescr (TableDescr *t_descr);
TblInfoPtr *makeColumnDescr (ColumnDescr *c_descr);
IRAFPointer saveIrafP (TblInfoPtr *tptr);
TableDescr *getTableDescr (IRAFPointer tp);
ColumnDescr *getColumnDescr (IRAFPointer cp);
char *expandfn (char *filename);

void setError (int status, char *msg);
void clearError (void);
int checkError (void);
int hstio_err (void);
int c_iraferr (void);
char *c_iraferrmsg (void);
