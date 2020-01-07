#ifndef INCL_STISDQ_H
#define INCL_STISDQ_H

/* Data quality flags for STIS. */

# define GOODPIXEL             0  /* OK */
# define SOFTERR               1  /* Reed-Solomon decoding error */
# define DATALOST              2  /* data replaced by fill value */
# define DETECTORPROB          4  /* bad detector pixel or beyond aperture */
# define DATAMASKED            8  /* masked by occulting bar */
# define HOTPIX               16  /* hot pixel */
# define LARGEBLEM            32  /* large blemish */
/* # define TBD                  64  */
# define OVERSCAN            128  /* can be used for finding bias level */
# define SATPIXEL            256  /* saturated pixel */
# define CALIBDEFECT         512  /* bad pixel in reference file */
# define SMALLBLEM          1024  /* small blemish */
# define X1D_BAD_BACKGROUND 2048  /* bad pixels in background region */
# define X1D_DISCARDED      4096  /* pixel discarded from extraction region */
# define DATAREJECT         8192  /* rejected during image combination */
# define NOT_CTI_CORR      16384  /* pixel not CTI corrected */

#endif /* INCL_STISDQ_H */
