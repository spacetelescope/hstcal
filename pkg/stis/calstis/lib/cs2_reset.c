# include "c_iraf.h"
# include "stis.h"
# include "cs2.h"


void cs2_reset (clpar *par, int newpar[])
{

/* arguments
clpar *par;		o: user specified parameters
int newpar;		o: parameter been set by the user?
*/

	/* First, set parameters specifiable in the command line to be
			non-existent. */
	par->tbname[0] = '\0';
	par->verbose = 0;
	par->printtime = 0;

	newpar[TOTAL] = 0;
	newpar[SCALENSE] = 0;
	newpar[INITGUES] = 0;
	newpar[SKYSUB] = 0;
	newpar[CRSIGMAS] = 0;
	newpar[CRRADIUS] = 0;
	newpar[CRTHRESH] = 0;
	newpar[BADINPDQ] = 0;
	newpar[CRMASK] = 0;
}
