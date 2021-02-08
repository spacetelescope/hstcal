#ifndef INCL_ACSDQ_H
#define INCL_ACSDQ_H

/* Data quality flags for CALACS. */

# define GOODPIXEL         0    /* OK */
# define SOFTERR           1    /* Reed-Solomon decoding error */
# define DATALOST          2    /* data replaced by fill value */
# define DETECTORPROB      4    /* bad detector pixel or beyond aperture */
# define DATAMASKED        8    /* masked by occulting bar */
# define HOTPIX           16    /* hot pixel */
# define LARGEBLEM        32    /* large blemish */
/* # define TBD           64    */
# define OVERSCAN        128    /* can be used for finding bias level */
# define SATPIXEL        256    /* full-well saturated pixel */
# define CALIBDEFECT     512    /* bad pixel in reference file */
# define SMALLBLEM      1024    /* small blemish */
# define ATODSAT        2048    /* a-to-d saturated pixel */
/* # define TBD         4096    */
# define DATAREJECT     8192    /* rejected during image combination */
/* # define TBD        16384    */
/* # define TBD        32768    */

#endif /* INCL_ACSDQ_H */
