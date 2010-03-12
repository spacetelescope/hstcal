# include <stdio.h>
# include "hstio.h"

# include "stis.h"
# include "cs2.h"
# include "calstis2.h"

/*  calstis2 -- cosmic ray rejection

  Description:
  ------------

  Date		Author		Description
  ----		------		-----------
  05-06-1996   J.-C. Hsu	Adapt from the SPP code t_crrej.x
  03-Sep-1999  Phil Hodge	Remove niter from calling sequence of crrej_do.
  10-Feb-2000  Phil Hodge	Replace exit with return;
				crrej_do is now int instead of void;
				print input (if only one) and output file names.
*/

int	CalStis2 (char *input, char *fout, clpar *par, int newpar[]) 
{

	IRAFPointer	tpin;
	char		*fdata;
	float		sigma[MAX_ITER];
	int		flag;

	int		crrej_do (IRAFPointer, char *, clpar *, int [],
				 float []);

	/*  announce start of the task */
        PrBegin (2);

	if (par->printtime)
            TimeStamp ("CALSTIS-2 started", "");

	/* make sure the output does not exist */
        flag = ckNewFile (fout);
        if (flag > 0) {
            if (flag == 1) {
                printf ("ERROR    Output file `%s' already exists.\n", fout);
                return (2);
            } else {
                printf ("ERROR    Can't clobber `%s'.\n", fout);
                return (2);
            }
        }

	/* open the input file template */
	tpin = c_imtopen (input);

	/* if there's only one input file, print its name */
	if (c_imtlen (tpin) == 1) {
	    fdata = malloc (STIS_FNAME * sizeof (char));
	    c_imtgetim (tpin, fdata, STIS_FNAME);
	    PrFileName ("input", fdata);
	    free (fdata);
	}
	PrFileName ("output", fout);

        printf ("\n");
        PrSwitch ("crcorr", PERFORM);

	/* perform the calculation */
	if (crrej_do (tpin, fout, par, newpar, sigma))
	    return (2);

        PrSwitch ("crcorr", COMPLETE);
	/*
	printf ("Created output file '%s'\n", fout);
	*/

	/* close file template */
	c_imtclose (tpin);

        printf ("\n");
        PrEnd (2);

	if (par->printtime)
            TimeStamp ("CALSTIS-2 complete", "");

	return (0);
}
