#ifndef INCL_STISSIZES_H
#define INCL_STISSIZES_H

/* Sizes of detectors. */

/* The size of the MAMA detectors. */
# define MAMA_NPIX_X     1024
# define MAMA_NPIX_Y     1024

/* The size of the CCD detector excluding overscan. */
# define CCD_NPIX_X      1024
# define CCD_NPIX_Y      1024

/* The size of the CCD overscan region on each edge of chip. */
# define OVERSCAN_LEFT     19	/* leading overscan in serial direction */
# define OVERSCAN_RIGHT    19	/* trailing overscan in serial direction */
# define OVERSCAN_BOTTOM    0
# define OVERSCAN_TOP      20	/* trailing overscan in parallel direction */

/* In CCD subarray mode, one pixel is lost all around the edge. */
# define SUBARRAY_BORDER    1

/* Use this range of pixels (zero indexed) for finding the bias level. */
# define BIASSECT_1L     1046
# define BIASSECT_1R     1060
# define BIASSECT_1L_B   1046
# define BIASSECT_1R_B   1052
# define BIASSECT_2L      524
# define BIASSECT_2R      530
# define BIASSECT_4L      263
# define BIASSECT_4R      269

#endif /* INCL_STISSIZES_H */
