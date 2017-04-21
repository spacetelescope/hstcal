# include   <stdio.h>
# include   <string.h>
# include   <stdlib.h>
# include   "hstio.h"
# include   <limits.h>

# include   "acs.h"    /* for message output */
# include   "err.h"

# define    MINVAL      -15000
# define    BIN_WIDTH   1
# define    MIN_BINS    1000

/* crrej_sky -- Calculate the sky for an image. */

void acsrej_sky (char *sky, IODescPtr ipsci[], IODescPtr ipdq[], 
        int nimgs, short badinpdq, float skyval[])
{
    extern int status;

    int         *histgrm;   /* pointer to the histogram */
    int         nbins;      /* number of bins */
    int         min, max;
    float       data_min;
    float       hwidth;     /* histogram resolution */
    float       hmin;       /* minimum histogram value */
    int         i, k, h;
    float       *a;
    short       *b;
    int         line, npt;
    int         dimx, dimy;
    float       sum, mean;
    Hdr         dqhdr;

    float   cr_mode (int *, int, float, float);

    /* -------------------------------- begin ---------------------------------- */
    /* decide what to do according to sky */
    if (strcmp (sky, "none") == 0) {
        for (k = 0; k < nimgs; ++k) {
            skyval[k] = 0.;
        }
        return;
    }
    if (strcmp (sky, "mode") != 0) {
        trlerror ("illegal sky value");
        status = INVALID_VALUE;
        return;
    }

    dimx = getNaxis1(ipsci[0]);
    dimy = getNaxis2(ipsci[0]);


    a = (float *) calloc (dimx, sizeof(float));
    b = (short *) calloc (dimx, sizeof(short));
    if (a == NULL || b == NULL) {
        trlerror ("Couldn't allocate memory for sky arrays");
        status = OUT_OF_MEMORY;
        return; 
    }


    /* compute MIN and MAX values for data */
    /* use the minimum and twice of the mean to determine the data 
    range */
    data_min = INT_MAX;
    sum = 0.;
    npt = 0;
    
    initHdr(&dqhdr);
    getHeader(ipdq[0],&dqhdr);

    for (line = 0; line < dimy; line++) {

        /* read the data in */
        getFloatLine (ipsci[0], line, a);
        getShortLine (ipdq[0], line, b);

        for (i = 0; i < dimx; ++i) {
            if ( (b[i] & badinpdq) == ACS_OK) {
                data_min = (a[i] < data_min) ? a[i] : data_min;
                sum += a[i];
                npt++;
            }
        }
    } /* End of loop over lines */

    freeHdr(&dqhdr);
    
    /* Compute min and max for histogram.
        MIN is min of good data or MINVAL, which ever is greater
        DELTA is difference between mean of data and MIN
        MAX is mean plus DELTA
        This insures that the mean falls in the center of the range
            between min and max.
    */
    if (npt == 0) npt = 1;
    min = (data_min < MINVAL) ? MINVAL : data_min;
    mean = (sum > 0.) ? (int) ( (sum / (float)npt) + 1) : 1;
    max = 2 * mean - min;
    
    /* use the mode as sky value, use the bin size of 1 (DN) */
    nbins = max - min + 1; 

    
    /* Insure that there are at least MIN_BINS in the histogram
        and reset the width accordingly.
    */ 
    if (nbins < MIN_BINS) {
        nbins = MIN_BINS;
        hwidth = (float) nbins / (float)MIN_BINS;
    } else {
        hwidth = 1.;
    }
    hmin = (float) min;

    /*
    sprintf(MsgText,"sky computed using min %d, max %d, mean %g, and bins %d",min,max,mean,nbins);
    trlmessage(MsgText);
    */

    for (k = 0; k < nimgs; ++k) {
        /*  
            set up a new histogram array for each image
        */
        histgrm = calloc (nbins, sizeof(int));
        if (histgrm == NULL){
            trlerror ("Couldn't allocate memory for sky histogram array");
            status = OUT_OF_MEMORY;
            return;   
        }

        initHdr(&dqhdr);
        getHeader(ipdq[k],&dqhdr);
        for (line = 0; line < dimy; line++) {

            /* read the data in */
            getFloatLine (ipsci[k], line, a);
            getShortLine (ipdq[k], line, b);

            /* construct the histogram */
            for (i = 0; i < dimx; ++i) {
                if (b[i] == ACS_OK) {
                    /* adjust the bin position by the width of each bin */
                    h = (int) ((a[i] - min) / hwidth);

	                if (h >= 0 && h < nbins)
                        histgrm[h]++;
	            }
            }
        } /* End of loop over lines */
        freeHdr(&dqhdr);
        
        /* calculate the mode from the histogram */
        skyval[k] = cr_mode (histgrm, nbins, hwidth, hmin);

        free (histgrm);

    }
    free(a);
    free(b);
}
