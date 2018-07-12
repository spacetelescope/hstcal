#ifndef INCL_IDTALG_H
#define INCL_IDTALG_H

# include <ximio.h>

# define	NITER		3
# define	NSBIG		3000
# define	NREFWAVE	3

/* Singles out a pixel in a 2-D Image object. */

# define PIX(a,i,j)	(a).pix[(j)*((a).nx) + (i)]
# define PPIX(a,i,j)	(a)->pix[(j)*((a)->nx) + (i)]


/* a single scattering function */

typedef struct {		
	int	sporder;	/* spectral order */
	int	nelem;		/* # of elements in function */
	double	*values;	/* function values */
} ScFunc;


/* a single ripple function */

typedef struct {
	int	sporder;	/* spectral order */
	int	nelem;		/* # of elements in function */
	double	*wavelengths;	/* wavelength values */
	double	*values;	/* function values */
} RippleFunc;


/* Scattering functions */

typedef struct {
	char opt_elem[200];	/* grating */ 
	int	    nsc;	/* # of elements in scattering array */
	ScFunc	*scfunc;	/* scattering array elements */
	int	    nrp;	/* # of elements in ripple array */
	RippleFunc *rpfunc;	/* ripple array elements */
	float     *spsf;	/* echelle x-disp function */
	int       nspsf;	/* and its size */
	float    *xdisp;	/* x-disp function USED IN MAKE_FFT */
	int      nxdisp;	/* and its size */
	double kernw[NREFWAVE+1];  /* wavelengths of halo functions */
	double  psfw[NREFWAVE+1];  /* wavelengths of telesc. PSF functions */
	int       nwave;	/* # of wavelengths in kernw and psfw */
	CmplxArray  ft1;	/* combined halo/PSF Fourier transforms */
	CmplxArray  ft2;
	CmplxArray  ft3;
	CmplxArray fto1;
	CmplxArray fto2;
	CmplxArray fto3;
} ScatterFunctions;


/* image array */

typedef struct {
	int       imio;		/* was buffer allocated by imio ? */
	IRAFPointer im;		/* image descriptor */
	int	 naxis;		/* # of axis */
	int	    nx;		/* 1st axis size */
	int	    ny;		/* 2nd axis size */
	int	  npix;		/* # of pixels */
	float     *pix;		/* pixel buffer */
} Image;


/* spliced spectrum */

typedef struct {
	int	npts;		/* # of bins in spectrum */
	double	*wmerge;	/* wavelength values */
	double	*fmerge;	/* flux (net) values */
} Spliced;


#endif /* INCL_IDTALG_H */
