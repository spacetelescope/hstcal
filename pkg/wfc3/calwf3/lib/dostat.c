# include <stdio.h>

# include "hstio.h"
# include "trl.h"		/* For trlwarn */
# include "err.h"
# include "wf3dq.h"

/* This routine computes the minimum, mean, and maximum of the science
   data values that are flagged as good in the data quality array.  These
   three quantities are updated in the science image header.  These values
   are also determined and saved for the error array.

   NOTE:  This version differs from the one in calstis1 in that sdqflags
   is a short int rather than unsigned short.  Also, the local variable
   flagval is short, rather than short *.

   Warren Hack, 1998 June 15:
   	Initial ACS version.
*/

int doStat (SingleGroup *out, short sdqflags) {

/* arguments:
SingleGroup *out  io: image to be calibrated; the headers are modified
short sdqflags     i: "serious" data quality flags
*/

	extern int status;

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
	int dimx, dimy;
    
	short flagval;			/* data quality flag value */
	int PutKeyFlt (Hdr *, char *, float, char *);
	int PutKeyInt (Hdr *, char *, int, char *);

	/* Statistics for the science data. */
	numgood = 0;
    valmin=0.0f;
    valmax=0.0f;
    errmin=0.0f;
    errmax=0.0f;
    snrmin=0.0f;
    snrmax=0.0f;
    
	num_bad_stddev = 0;
	valsum = 0.;
	errsum = 0.;
	snrsum = 0.;
	dimx = out->sci.data.nx;
	dimy = out->sci.data.ny;

	for (j = 0;  j < dimy;  j++) {
	     for (i = 0;  i < dimx;  i++) {
		  flagval = DQPix (out->dq.data, i, j);

		  if (!(sdqflags & flagval)) {

		      /* no serious flag bit set */
		      value = Pix (out->sci.data, i, j);
		      stddev = Pix (out->err.data, i, j);
		      if (stddev <= 0.) {
			  num_bad_stddev++;
			  continue;		/* bad error value */
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
	    area = dimy * dimx;
	    if (area == 0) {
		trlwarn ("Output image size is zero.");
	    } else if (num_bad_stddev > 0) {
		if (num_bad_stddev == area) {
		    trlwarn ("No ERR values > 0.");
		} else {
		    trlwarn 
		       ("All output pixels either flagged as bad or ERR <= 0.");
		}
	    } else {
		trlwarn ("All output pixels flagged as bad.");
	    }
	    PutKeyInt (&out->sci.hdr, "NGOODPIX", numgood, "");
	    PutKeyInt (&out->err.hdr, "NGOODPIX", numgood, "");
	    return (status);
	}

	/* Update header values for the science array. */

	if (PutKeyInt (&out->sci.hdr, "NGOODPIX", numgood,
			"number of good pixels"))
	    return (status);

	if (PutKeyFlt (&out->sci.hdr, "GOODMIN", (float) valmin,
			"minimum good data value"))
	    return (status);

	if (PutKeyFlt (&out->sci.hdr, "GOODMAX", (float) valmax,
			"maximum good data value"))
	    return (status);

	if (PutKeyFlt (&out->sci.hdr, "GOODMEAN", (float) valsum,
			"average of good data values"))
	    return (status);

	if (PutKeyFlt (&out->sci.hdr, "SNRMIN", (float) snrmin,
			"minimum S/N of good data values"))
	    return (status);

	if (PutKeyFlt (&out->sci.hdr, "SNRMAX", (float) snrmax,
			"maximum S/N of good data values"))
	    return (status);

	if (PutKeyFlt (&out->sci.hdr, "SNRMEAN", (float) snrsum,
			"mean S/N of good data values"))
	    return (status);

	/* Update header values for the error array. */

	if (PutKeyInt (&out->err.hdr, "NGOODPIX", numgood,
			"number of good pixels"))
	    return (status);

	if (PutKeyFlt (&out->err.hdr, "GOODMIN", (float) errmin,
			"minimum sigma for good data"))
	    return (status);

	if (PutKeyFlt (&out->err.hdr, "GOODMAX", (float) errmax,
			"maximum sigma for good data"))
	    return (status);

	if (PutKeyFlt (&out->err.hdr, "GOODMEAN", (float) errsum,
			"average of sigma for good data"))
	    return (status);

	return (status);
}
