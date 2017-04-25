# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"
# include "stisdq.h"
# include "stisdef.h"

/* This routine computes the minimum, mean, and maximum of the science
   data values that are flagged as good in the data quality array.  These
   three quantities are updated in the science image header.  These values
   are also determined and saved for the error array.

   NOTE:  This version differs from the one in calstis1 in that sdqflags
   is a short int rather than unsigned short.  Also, the local variable
   flagval is short, rather than short *.

   Phil Hodge, 1997 Sept 12:
	The message about "all output pixels flagged as bad" was printed
	in cases where it wasn't really correct; print other messages if
	appropriate.

   Phil Hodge, 1998 Oct 7:
	Allow the error array value to be zero if the science array value
	is less than or equal to zero; in this case the signal-to-noise
	will be set to zero.
*/

int doStat (SingleGroup *out, short sdqflags) {

/* arguments:
SingleGroup *out  io: image to be calibrated; the headers are modified
short sdqflags    i: "serious" data quality flags
*/

	int status;

	double value;			/* current data value */
	double valsum, valmin, valmax;
	double stddev;			/* current error estimate */
	double errsum, errmin, errmax;
	double snr;			/* current signal-to-noise ratio */
	double snrsum, snrmin, snrmax;
	int numgood;			/* number of good pixels */
	int num_bad_stddev;		/* number of pixels with err = 0 */
	int area;			/* total number of pixels */
	int i, j;
	short flagval;			/* data quality flag value */

	/* Statistics for the science data. */

	numgood = 0;
	num_bad_stddev = 0;
	valsum = 0.;
	errsum = 0.;
	snrsum = 0.;
	for (j = 0;  j < out->sci.data.ny;  j++) {
	    for (i = 0;  i < out->sci.data.nx;  i++) {
		flagval = DQPix (out->dq.data, i, j);
		if ( ! (sdqflags & flagval) ) {
		    /* no serious flag bit set */
		    value = Pix (out->sci.data, i, j);
		    stddev = Pix (out->err.data, i, j);
		    if (stddev <= 0.) {
			if (value <= 0.) {
			    snr = 0.;		/* S/N of 0/0 = 0 */
			} else {
			    num_bad_stddev++;
			    continue;		/* bad error value */
			}
		    } else {
			snr = value / stddev;
		    }
		    if (numgood < 1) {
			valsum = value;
			valmin = value;
			valmax = value;
			errsum = stddev;
			errmin = stddev;
			errmax = stddev;
			snrsum = snr;
			snrmin = snr;
			snrmax = snr;
			numgood = 1;
		    } else {
			valsum += value;
			errsum += stddev;
			snrsum += snr;
			if (value < valmin)
			    valmin = value;
			if (value > valmax)
			    valmax = value;
			if (stddev < errmin)
			    errmin = stddev;
			if (stddev > errmax)
			    errmax = stddev;
			if (snr < snrmin)
			    snrmin = snr;
			if (snr > snrmax)
			    snrmax = snr;
			numgood++;
		    }
		}
	    }
	}

	if (numgood > 0) {
	    valsum /= (double) numgood;
	    errsum /= (double) numgood;
	    snrsum /= (double) numgood;
	} else {
	    area = out->sci.data.ny * out->sci.data.nx;
	    if (area == 0) {
		printf ("Warning  Output image size is zero.\n");
	    } else if (num_bad_stddev > 0) {
		if (num_bad_stddev == area) {
		    printf ("Warning  No ERR values > 0.\n");
		} else {
		    printf (
		"Warning  All output pixels either flagged as bad \\\n");
		    printf ("Warning  or ERR <= 0.\n");
		}
	    } else {
		printf ("Warning  All output pixels flagged as bad.\n");
	    }
	    status = Put_KeyI (&out->sci.hdr, "NGOODPIX", numgood, "");
	    status = Put_KeyI (&out->err.hdr, "NGOODPIX", numgood, "");
	    return (0);
	}

	/* Update header values for the science array. */

	if ((status = Put_KeyI (&out->sci.hdr, "NGOODPIX", numgood,
                                "number of good pixels")))
	    return (status);

	if ((status = Put_KeyF (&out->sci.hdr, "GOODMIN", (float) valmin,
                                "minimum good data value")))
	    return (status);

	if ((status = Put_KeyF (&out->sci.hdr, "GOODMAX", (float) valmax,
                                "maximum good data value")))
	    return (status);

	if ((status = Put_KeyF (&out->sci.hdr, "GOODMEAN", (float) valsum,
                                "average of good data values")))
	    return (status);

	if ((status = Put_KeyF (&out->sci.hdr, "SNRMIN", (float) snrmin,
                                "minimum S/N of good data values")))
	    return (status);

	if ((status = Put_KeyF (&out->sci.hdr, "SNRMAX", (float) snrmax,
                                "maximum S/N of good data values")))
	    return (status);

	if ((status = Put_KeyF (&out->sci.hdr, "SNRMEAN", (float) snrsum,
                                "mean S/N of good data values")))
	    return (status);

	/* Update header values for the error array. */

	if ((status = Put_KeyI (&out->err.hdr, "NGOODPIX", numgood,
                                "number of good pixels")))
	    return (status);

	if ((status = Put_KeyF (&out->err.hdr, "GOODMIN", (float) errmin,
                                "minimum sigma for good data")))
	    return (status);

	if ((status = Put_KeyF (&out->err.hdr, "GOODMAX", (float) errmax,
                                "maximum sigma for good data")))
	    return (status);

	if ((status = Put_KeyF (&out->err.hdr, "GOODMEAN", (float) errsum,
                                "average of sigma for good data")))
	    return (status);

	return (0);
}
