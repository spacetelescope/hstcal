#ifndef INCL_ACSCORR_H
#define INCL_ACSCORR_H

/* calibration switches for 'acsccd', 'acscte', and 'acs2d' */

typedef struct {
    int dqicorr;
    int atodcorr;
    int blevcorr;
    int biascorr;
} acsccd_switch;

typedef struct {
    int pctecorr;
} acscte_switch;

typedef struct {
    int dqicorr;
    int glincorr;
    int lflgcorr;
    int darkcorr;
    int flashcorr;
    int flatcorr;
    int shadcorr;
    int photcorr;
} acs2d_switch;

/* All the calibration switches. */

typedef struct {
    int atodcorr;
    int biascorr;
    int blevcorr;
    int crcorr;
    int darkcorr;
    int flashcorr;
    int dqicorr;
    int flatcorr;
    int pctecorr;
    int glincorr;
    int lflgcorr;
    int photcorr;
    int rptcorr;
    int shadcorr;
    int expscorr;
} CalSwitch;

#endif /* INCL_ACSCORR_H */
