# include 	<stdio.h>
# include 	"hstio.h"

# include	"stis.h"

/*  cr_scaling -- Determine the scaling factors according to exposure times or
		 other user specified scheme. */

/*
Description:
------------
If using the exposure time, the scaling factors are normalized to ratios
relative to the max exposure.

  Date		Author		Description
  ----		------		-----------
  02-May-1996  J.-C. Hsu	Adapt from the SPP code cr_scaling.x
  10-Feb-2000  Phil Hodge	Add exp_range to the calling sequence,
				and compute its values; replace get_key_f
				with getKeyF; replace 'exit' with 'return',
				and return int instead of void.
  22-May-2012  Phil Hodge	Change the declaration of imgname.
*/

int cr_scaling (char *expname, IODescPtr ipsci[], char *imgname[],
		int nimgs, float efac[], float tfac[], double exp_range[])
{
	Hdr	scihdr;
	int	nzero, k;
	double	expstart, expend;	/* MJD of exposure start and stop */

/* -------------------------------- begin ---------------------------------- */

        initHdr (&scihdr);

	/* if the parameter scaling is null, all images have equal weight. */
	if (expname[0] == '\0') {
	    for (k = 0; k < nimgs; ++k)
		efac[k] = 1.;
	    return (0);
	}

	/* Use exposure time as scaling factor */
	nzero = 0;
	for (k = 0; k < nimgs; ++k) {
	    getHeader (ipsci[k], &scihdr);

	    if (getKeyF(&scihdr, expname, &efac[k]) != 0) {
		printf ("cannot read the keyword '%s' from the SCI "
			"extension of '%s'\n", expname, imgname[k]);
		return (2);
	    }
	    if (efac[k] < 0.) {
		printf ("exposure time of file '%s' is negative\n",
		    	imgname[k]);
		return (2);
	    }
	    if (efac[k] == 0.) {
		nzero++;
	    }
	    if (getKeyF(&scihdr, "OCCDHTAV", &tfac[k]) != 0) {
		tfac[k] = -1.;
	    }

	    /* Update min exposure start and max exposure stop times. */
	    if (getKeyD (&scihdr, "EXPSTART", &expstart))
		expstart = 0.;
	    if (getKeyD (&scihdr, "EXPEND", &expend))
		expend = 0.;
	    if (k == 0) {
		exp_range[0] = expstart;
		exp_range[1] = expend;
	    } else {
		if (expstart < exp_range[0])
		    exp_range[0] = expstart;
		if (expend > exp_range[1])
		    exp_range[1] = expend;
	    }

	    freeHdr (&scihdr);
	}
	if (nzero > 0 && nzero < nimgs) {
	    printf ("some (but not all) input imsets have zero exposure time\n");
	    return (2);
	}

	return (0);
}
