# include <stdio.h>
# include <stdlib.h>	/* atoi, atof */
# include <string.h>

# include "xtables.h"

# include "stis.h"
# include "calstis6.h"
# include "err.h"

static void RemoveBlanks (char *, char *);
static double WaveToPix (double, RowContents *);

/* 
   Translates a range string into a rejection flag array.

   This is used by the profile generator to reject wavelength ranges
   in the input image. This function presumes the necessary arrays
   were previously allocated and initialized elsewhere.




   Revision history:
   ----------------
   29 Jun 98  -  Implemented (I.Busko)
   28 Feb 01  -  Warn if invalid range is found (IB)
*/

void SetRanges (StisInfo6 *sts, RowContents *row) {

	char str[STIS_LINE], *colon, *range;
	double wave;
	int i, index1, index2;

	/* Clean up string from any embedded blanks. */
	RemoveBlanks (sts->rejranges, str);

	/* Scan comma-separated ranges. */
	range = strtok (str, ",");
	while (range != NULL) {

	    /* See if there is a colon separator. */
	    colon = strchr (range, ':');

	    /* Found colon, assume it's a range. */
	    if (colon != NULL) {

	        /* Cannot use strtok again... */
	        *colon = '\0';
	        if ((wave = atof (range)) == 0.0)
	            continue;
	        index1 = WaveToPix (wave, row);
	        if ((wave = atof (colon+1)) == 0.0)
	            continue;
	        index2 = WaveToPix (wave, row);

	        /* Indexes might be inverted. */
	        if (index1 > index2) {
	            i = index2;
	            index2 = index1;
	            index1 = i;
	        }

	        /* Check for valid range. */
                if (index1 < 0 || index2 >= row->npts) {
	            printf ("Warning  Invalid wavelength range.\n");
	            continue;
	        }

	        /* Set all flags in range. */
	        for (i = index1-1; i < index2; i++) {
	        /*  if (i >= 1 && i <= row->npts)  not required anymore */
	                sts->profile_rej[i-1] = 1;
	        }

	    /* No colon detected, assume it is a single pixel index. */
	    } else {
	        if ((wave = atof (range)) == 0.0)
	            continue;
	        index1 = WaveToPix (wave, row);

	        /* Set flag in case of valid index. */
	        if (index1 > 0 && index1 <= row->npts)
	            sts->profile_rej[index1-1] = 1;
	    }

	    /* Get next range. */
	    range = strtok (NULL, ",");
	}
}



/*  Removes embedded blanks. */

static void RemoveBlanks (char *in, char *out) {

	int i, k;

	k = 0;
	for (i = 0; i < strlen (in); i++) {
	    if (in[i] != ' ')
	        out[k++] = in[i];
	}
	out[k] = '\0';
}



/*  Converts wavelength to pixel. */

static double WaveToPix (double wave, RowContents *row) {

	return ((wave - row->wave[0]) / 
                (row->wave[row->npts-1] - row->wave[0]) * row->npts);
}

