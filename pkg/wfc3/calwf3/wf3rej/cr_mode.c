# define    ABS(x)      ((x>0.) ? x : -x)

# define    EPSILON     1.e-20

/* CR_MODE -- Procedure to compute the mode from the histogram. 

    Adapted from images$iminfo/t_imstat.x
    May-03-1996	J.-C. Hsu
*/

float cr_mode (int *histgrm, int nbins, float hwidth, float hmin) {

    int     i, bpeak;
    float   mode, hpeak, dh1, dh2, denom;

    /* If there is a single bin return the midpoint of that bin. */
    if (nbins == 1)
        return (hmin + 0.5 * hwidth);

    /*  If there are two bins return the midpoint of the greater bin. */
    if (nbins == 2) {
        if (histgrm[0] > histgrm[1])
            mode = hmin + 0.5 * hwidth;
        else if (histgrm[1] > histgrm[0])
            mode = hmin + 1.5 * hwidth;
        else
            mode = hmin + hwidth;
        return (mode);
    }

    /* Find the bin containing the histogram maximum. */
    hpeak = histgrm[0];
    bpeak = 0;
    for (i = 1; i < nbins; ++i) {
        if (histgrm[i] > hpeak) {
            hpeak = histgrm[i];
            bpeak = i;
        }
    }

    /* If the maximum is in the first bin return the midpoint of the bin. */
    if (bpeak == 0)
        return (hmin + 0.5 * hwidth);

    /* If the maximum is in the last bin return the midpoint of the bin. */
    if (bpeak == (nbins-1))
        return (hmin + (nbins - 0.5) * hwidth);

    /* Do a parabolic interpolation to find the peak. */
    dh1 = histgrm[bpeak] - histgrm[bpeak-1];
    dh2 = histgrm[bpeak] - histgrm[bpeak+1];
    denom = dh1 + dh2;
    if (ABS(denom) < EPSILON) {
        mode = hmin + (bpeak + 0.5) * hwidth;
    } else {
        mode = bpeak + 0.5 + 0.5 * (dh1 - dh2) / denom;
        mode = hmin + mode * hwidth;
    }
    return (mode);
}
