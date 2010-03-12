/* should get these from c_iraf.h */

/* iomode */
# define IRAF_READ_ONLY   1
# define IRAF_READ_WRITE  2
# define IRAF_NEW_FILE    5
# define IRAF_TEMP_FILE   6
# define IRAF_NEW_COPY    7

/* data types for tables */
# define IRAF_BOOL        1
# define IRAF_CHAR        2
# define IRAF_SHORT       3
# define IRAF_INT         4
# define IRAF_REAL        6
# define IRAF_DOUBLE      7

typedef int IRAFPointer;

IRAFPointer c_tbtopn (char *tablename, int iomode, IRAFPointer template);
