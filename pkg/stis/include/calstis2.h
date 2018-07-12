#ifndef INCL_CALSTIS2_H
#define INCL_CALSTIS2_H

# include <c_iraf.h>
# include <xclio.h>
# include <ximio.h>
 
# define	MAX_ITER	20

# define	MAX_NEXP	3
# define	NONEXIST	-10000.

# define	IM_MAXDIM	7

# define	PIX(v,i,j,nx)	v[(i) + (j) * (nx)]
# define SQ(x)  ((x) * (x))

#endif /* INCL_CALSTIS2_H */
