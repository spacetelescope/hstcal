#ifndef INCL_STISPHT_H
#define INCL_STISPHT_H

/* This contains factors for converting from count rate to absolute flux,
   read from the photometry table _pht.

   This was moved out from calstis6.h and calstis7.h to include the
   blaze shift correction parameters (IB, 04Feb02).
*/

typedef struct {
	char bunit[STIS_FITS_REC];	/* units for flux-calibrated data */
	int allocated;		/* true if memory has been allocated */
	int nelem;		/* size of wavelength and flux corr arrays */
	double *wl;		/* array of wavelengths */
	double *thru;		/* array of throughputs (QE) */
	double *error;		/* array of throughput errors */
	double *pcorr;		/* array of correction factors */

	int blazecorr;		/* calibration switch */

	int mref;		/* reference order # of sens. function */
	double wref;		/* central wavel. of mref order in sens. func.*/
	double yref;		/* A2 position of mref order in sens. func.*/
	double mjd;		/* observ. date of sens. function */

	double mx;		/* x shift scale factor */
	double my;		/* y shift scale factor */
	double mt;		/* time shift scale factor (pixel/day)*/
	double m0;		/* blaze shift zero point offset */

	double wpos;		/* wavel. of pixel 512 in order mref in data */
	double ypos;		/* A2 position of mref order in data */
	double disp;		/* dispersion for order mref in data */
} PhotInfo;

#endif /* INCL_STISPHT_H */
