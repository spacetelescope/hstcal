# include <stdio.h>

# include "stis.h"
# include "stispht.h"

/* This routine applies the MSM/blaze correction to the photometry
   wavelength array.



   Revision history:
   ----------------
   06 Feb 2002  -  Implemented (I. Busko)
   16 Apr 2002  -  Parameter from command line (IB)
   17 Apr 2002  -  Output computed blazeshift value (IB)
   27 Jun 2006  -  Include zero point offset (PEH)

*/

void BlazeCorr (PhotInfo *phot, double start, double end , int sporder,
                double blazeshift, double *bs, double dispersion) {
/* arguments:
PhotInfo *phot		io: photometry data structure
double start;		i:  exposure start time, in MJD
double end;		i:  exposure end time, in MJD
int sporder;		i:  spectral order number
double blazeshift;	i:  override value from command line (in pixels)
double *bs;		o:  blaze shift actually computed/used
double dispersion;	i:  dispersion (used only when no reference values
                            are available
*/
	double date, deltamw, dw, dx, dy, dt, bshift;
	int i;

	/* This is the case where we have valid table values and 
	   no command-line override. 
	*/
	if (phot->blazecorr == PERFORM && blazeshift == NO_VALUE) {

	    date    = (end - start) / 2.0 + start;
	    deltamw = phot->disp * phot->mref;

	    dw = phot->wref - phot->wpos;
	    dx = dw / phot->disp;
	    dy = phot->ypos - phot->yref;
	    dt = date - phot->mjd;

	    bshift = dx * phot->mx + dy * phot->my + dt * phot->mt + phot->m0;

	    for (i = 0; i < phot->nelem; 
	        phot->wl[i++] += bshift * deltamw / sporder);

	    *bs = bshift;

	    printf (
	    "         Blaze shift of %g pixels applied to sensitivity curve.\n",
	     bshift);

	} else if (blazeshift != NO_VALUE) {

	    /* This is the case where we have a command-line override.
	       It is applied no matter what the table values are.
	    */

	    for (i = 0; i < phot->nelem; 
	        phot->wl[i++] += blazeshift * dispersion);

	    *bs = blazeshift;

	    printf (
	    "         Blaze shift of %g pixels applied to sensitivity curve.\n",
	     blazeshift);
	} 
}
