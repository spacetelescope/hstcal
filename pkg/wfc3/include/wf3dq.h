#ifndef INCL_WF3DQ_H
#define INCL_WF3DQ_H

/* Data quality flags for CALWF3. */
/* 2009-01-28 	H.Bushouse	New WFC3 DQ assignments (PR 61741). 
   2015  M. Sosey               Sink pixels added to the mask, see #1093
*/

# define GOODPIXEL         0    /* OK */
# define SOFTERR           1    /* Reed-Solomon decoding error */
# define DATALOST          2    /* data replaced by fill value */
# define DETECTORPROB      4    /* bad detector pixel or beyond aperture */
# define DATAMASKED        8    /* masked by occulting bar */
# define BADZERO           8    /* deviant IR zero-read pixel */
# define HOTPIX           16    /* hot pixel */
# define CTETAIL          32    /* UVIS CTE tail */
# define UNSTABLE         32    /* IR unstable pixel */
# define WARMPIX          64    /* warm pixel */
# define BADBIAS         128    /* bad bias value */
# define SATPIXEL        256    /* full-well saturated pixel */
# define BADFLAT         512    /* bad flatfield value */
# define TRAP           1024    /* UVIS charge trap, SINK pixel */
# define SPIKE          1024    /* CR spike detected during cridcalc IR */ 
# define ATODSAT        2048    /* a-to-d saturated pixel */
# define ZEROSIG        2048    /* IR zero-read signal correction */
/* # define TBD         4096       reserved for Multidrizzle CR rej */
# define DATAREJECT     8192    /* rejected during image combination UVIS, IR CR rejection* */
# define HIGH_CURVATURE 16384   /* pixel has more than max CR's  */
# define CROSSTALK     16384   /* ghost or crosstalk  */
# define RESERVED2     32768    /* can't use */

#endif /* INCL_WF3DQ_H */
