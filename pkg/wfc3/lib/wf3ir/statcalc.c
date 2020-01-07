# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>

# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"

extern int status;

/* STATCALC: Compute statistics for science images. The min, max,
** mean and standard deviation of "good" (unflagged) pixels are computed
** for the entire image, as well as each quadrant. A tally of the
** number of pixels with a each DQ flag is also computed. All results are
** written to keywords in the science image extension header. This
** routine does NOT modify the input image data in any way.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	20-Mar-2002	Deactivated computation of NQUAL*, MEDN,
**				and STDV values.
** H.Bushouse	02-Apr-2002	Modified to compute STIS, ACS-style
**				statistics. Also confirmed readiness to
**				support WFC3/IR subarrays.
** H.Bushouse	10-Apr-2002	Skip IR reference pixels by basing limits
**				on trim values from OSCNTAB.
*/

int statcalc (WF3Info *wf3, SingleNicmosGroup *input, short sdqflags) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: input image
**	sdqflags i: serious data quality flags
*/

	/* Local variables */
	int   ibeg, iend, jbeg, jend;	/* image limits */
	int   npix;			/* number of good pixel values */
	float mean,  min,  max;		/* statistics values */
	float emean, emin, emax;
	float smean, smin, smax;

	/* Function definitions */
	int qstats (SingleNicmosGroup *, int, int, int, int, short, int *,
		    float *, float *, float *, float *, float *, float *,
		    float *, float *, float *);

	/* Set boundaries for statistics computation */
	ibeg = wf3->trimx[0]; iend = input->sci.data.nx - wf3->trimx[1] - 1;
	jbeg = wf3->trimy[0]; jend = input->sci.data.ny - wf3->trimy[1] - 1;

	/* Calculate statistics for the SCI and ERR images */
	if (qstats (input, ibeg, iend, jbeg, jend, sdqflags, 
		    &npix, &mean, &min, &max, &emean, &emin, &emax,
		    &smean, &smin, &smax))
	    return (status = 1);

	/* Store results in header keywords */
	if (npix < 1) {
	    if (putKeyI (&input->sci.hdr, "NGOODPIX", npix, ""))
		return (status = 1);
	    if (putKeyI (&input->err.hdr, "NGOODPIX", npix, ""))
		return (status = 1);
	    return (status = 0);
	}

	if (putKeyI (&input->sci.hdr, "NGOODPIX", npix, ""))
	    return (status = 1);
	if (putKeyF (&input->sci.hdr, "GOODMEAN", mean, ""))
	    return (status = 1);
	if (putKeyF (&input->sci.hdr, "GOODMIN",  min, ""))
	    return (status = 1);
	if (putKeyF (&input->sci.hdr, "GOODMAX",  max, ""))
	    return (status = 1);
	if (putKeyF (&input->sci.hdr, "SNRMEAN",  smean, ""))
	    return (status = 1);
	if (putKeyF (&input->sci.hdr, "SNRMIN",   smin, ""))
	    return (status = 1);
	if (putKeyF (&input->sci.hdr, "SNRMAX",   smax, ""))
	    return (status = 1);

	if (putKeyI (&input->err.hdr, "NGOODPIX", npix, ""))
	    return (status = 1);
	if (putKeyF (&input->err.hdr, "GOODMEAN", emean, ""))
	    return (status = 1);
	if (putKeyF (&input->err.hdr, "GOODMIN",  emin, ""))
	    return (status = 1);
	if (putKeyF (&input->err.hdr, "GOODMAX",  emax, ""))
	    return (status = 1);

	/* Successful return */
	return (status = 0);
}

