/* This routine uses the LTV and LTM keyword values to compute the corner
   location and pixel size.  The corner location
   and pixel size are in units of unbinned pixels. 

   ltm[0] and ltm[1] are assumed to be greater than zero.
*/

# include "trlbuf.h"
# include "wf3.h"

# define NINT(x)  ((x >= 0.) ? (int) (x + 0.5) : (int) (x - 0.5))

int FromLT (int rsize, double *ltm, double *ltv, int *bin, int *corner) {

/* arguments:
int rsize        i: reference pixel size
double ltm[2]    i: diagonal elements of MWCS matrix
double ltv[2]    i: MWCS linear transformation vector
int bin[2]       o: pixel size in X and Y
int corner[2]    o: corner of subarray in X and Y
*/

	extern int status;
	double dbinx, dbiny, dxcorner, dycorner;

	dbinx = (double)rsize / ltm[0];
	dbiny = (double)rsize / ltm[1];

	dxcorner = (dbinx - rsize) - dbinx * ltv[0];
	dycorner = (dbiny - rsize) - dbiny * ltv[1];

	/* Round off to the nearest integer. */
	corner[0] = NINT (dxcorner);
	corner[1] = NINT (dycorner);
	bin[0] = NINT (dbinx);
	bin[1] = NINT (dbiny);
    
	return (status);
}
