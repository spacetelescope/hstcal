static int BinarySearch (double, double [], int);

/* This routine finds the index (*index) in the wl array such that
   wl[*index] < wavelength <= wl[*index+1]
   Searching starts with index *index, so that value must have been
   assigned.  Before calling FindWL for the first time, *index should
   be initialized to -1.

   An output value of -1 for *index means the wavelength is smaller than
   the first element of wl (i.e. wl[0]); while *index = nelem+1 means the
   wavelength is greater than the last element of wl.

   Phil Hodge, 2001 Feb 23:
	Add a check for the case that the starting index is larger
	than the matching index.
*/

void FindWL (double wavelength, double wl[], int nelem, int *index) {

/* arguments:
double wavelength   i: wavelength to be found within the wl array
double wl[nelem+1]  i: wavelengths for the template spectrum (edges of pixels)
int nelem           i: nelem+1 is size of wl array
int *index         io: starting index; location of wavelength in wl array
*/

	int i;
	int start;

	if (wavelength < wl[0]) {

	    *index = -1;

	} else if (wavelength > wl[nelem]) {

	    *index = nelem + 1;

	} else if (nelem == 1 || wavelength == wl[0]) {

	    *index = 0;

	} else {

	    if (*index < 0 || *index > nelem) {

		*index = BinarySearch (wavelength, wl, nelem+1);

	    } else {

		start = *index;

		if (wavelength <= wl[start]) {
		    *index = BinarySearch (wavelength, wl, nelem+1);
		} else {
		    for (i = start;  i < nelem;  i++) {
			if (wavelength > wl[i] && wavelength <= wl[i+1]) {
			    break;
			}
		    }
		    *index = i;
		}
	    }
	}
}

/* This function does a binary search and returns the index i such that
   x is between wl[i] and wl[i+1].  The returned value will be out of
   range (-1 or n) if x is outside the wl array.
*/

static int BinarySearch (double x, double wl[], int n) {

/* arguments:
double x        i: value to be found within the wl array
double wl[n+1]  i: wavelengths at edges of pixels
int n           i: n+1 is size of wl array
*/

	int low, high;		/* range of elements to consider */
	int k;			/* middle element between low and high */

	if (x < wl[0])
	    return (-1);
	if (x >= wl[n-1])
	    return (n);

	low = 0;
	high = n - 1;

	while (high - low > 1) {

	    k = (low + high) / 2;
	    if (x < wl[k]) {
		high = k;
	    } else if (x > wl[k+1]) {
		low = k;
	    } else {
		high = k;
		low = k;
	    }
	}

	if (low < high) {
	    if (wl[high] <= x)
		return (high);
	    else
		return (low);
	} else {
	    return (low);
	}
}
