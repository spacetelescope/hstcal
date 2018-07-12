# include	<stdio.h>
# include	<string.h>
# include 	<stdlib.h>
# include	"hstio.h"

# include	"stis.h"
# include	"calstis2.h"
# include	"cs2.h"

/*  crrej_init -- get the initial average pixel values

  Description:
  ------------
  Get the initial average according to the specified scheme
  
  Date		Author		Description
  ----		------		-----------
  26-Aug-1993  J.-C. Hsu	design and coding
  10-Feb-2000  Phil Hodge	Check for out of memory; return int.
   8-May-2000  Phil Hodge	Pass zero-indexed line number to getFloatLine.
   7-Oct-2002  Phil Hodge	Copy John Blakeslee's changes to acsrej_init.c
				to this file, to handle scalense correctly;
				i.e. subtract the sky before applying scalense.
*/

int crrej_init (IODescPtr ipsci[], clpar *par, int nimgs, int dim_x, int dim_y,
		float noise[], float gain[], float efac[], float skyval[], 
		FloatTwoDArray *ave, FloatTwoDArray *avevar,   float *work)
{

	float		rog2, exp2, scale, val, raw, raw0, signal0;
	float		*buf;
	int		i, j, n;
	int		dum;
	int		*npts;

	void		piksrt (float [], int);

/* -------------------------------- begin ---------------------------------- */

	scale = par->scalenoise / 100.;

	npts = calloc (dim_x, sizeof(int));
	buf = calloc (dim_x, sizeof(float));
	if (npts == NULL || buf == NULL) {
	    printf ("ERROR    out of memory in crrej_init\n");
	    return (2);
	}

	/* use the stack median to construct the initial average */
	if (strncmp(par->initial, "median", 3) == 0) {
	    for (j = 0; j < dim_y; j++) {
		for (i = 0; i < dim_x; i++)
		    npts[i] = 0;

		for (n = 0; n < nimgs; n++) {
		    getFloatLine (ipsci[n], j, buf);
		    for (i = 0; i < dim_x; i++) {
			PIX(work,npts[i],i,nimgs) = (buf[i] - skyval[n]) / 
							efac[n];
			npts[i] += 1;
		    }
		}
		for (i = 0; i < dim_x; i++) {
		    dum = npts[i];
		    if (dum == 0)
			PPix(ave,i,j) = 0.0F;
		    else {
			piksrt (&PIX(work,0,i,nimgs), dum);
			if ((dum/2)*2 == dum) 
			    PPix(ave,i,j) = (PIX(work,dum/2-1,i,nimgs) + 
						  PIX(work,dum/2,i,nimgs))/2.;
			else
			    PPix(ave,i,j) = PIX(work,dum/2,i,nimgs);
		    }
		}
	    }

	/* use the minimum to construct the initial average */
	} else {
	    if (strncmp(par->initial, "minimum", 3) != 0) {
		printf ("Invalid INITGUES value %s, reset it to 'minimum'\n",
			par->initial);
		strcpy (par->initial, "minimum");
	    }
	    for (j = 0; j < dim_y; j++) {
		for (n = 0; n < nimgs; n++) {
		    exp2 = SQ(efac[n]);
		    rog2 = SQ(noise[n]);
		    getFloatLine (ipsci[n], j, buf);
		    for (i = 0; i < dim_x; i++) {
			raw = buf[i];
			raw0 = (raw > 0.)? raw : 0.;
			signal0 = ((raw - skyval[n]) > 0.) ?
				   (raw - skyval[n]) : 0.;
			val = (raw - skyval[n]) / efac[n];
			if (n == 0) {
			    PPix(ave,i,j) = val;
			    PPix(avevar,i,j) = (rog2 + raw0/gain[n] + 
						SQ(scale*signal0)) / exp2;
			} else if (val < PPix(ave,i,j)) {
			    PPix(ave,i,j) = val;
			    PPix(avevar,i,j) = (rog2 + raw0/gain[n] + 
						SQ(scale*signal0)) / exp2;
			}
		    }
		}
	    }
	}

	/* free the memory */
	free (npts);
	free (buf);

	return (0);
}
