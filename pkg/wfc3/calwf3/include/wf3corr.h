#ifndef INCL_WF3CORR_H
#define INCL_WF3CORR_H

/* calibration switches for 'wf3ccd' and 'wf3ir' */

/* Revision history:
**
** H.Bushouse	17 June 2002: Removed "statcorr" from CCD_Switch and
**		and IR_Switch (tracking CALACS changes).
** M. Sosey 4 Dec 2013: added fluxcorr for UVIS
**
** M. Sosey August 2014: added pctecorr for UVIS
*/

typedef struct {

	int atodcorr;
	int biascorr;
	int blevcorr;
	int crcorr;
	int darkcorr;
	int dqicorr;
	int flatcorr;
	int flashcorr;
    int fluxcorr;
    int pctecorr;
	int photcorr;
	int rptcorr;
	int shadcorr;
	int expscorr;

} CCD_Switch;

typedef struct {

	int zsigcorr;
	int zoffcorr;
	int dqicorr;
	int blevcorr;
	int noiscorr;
	int darkcorr;
	int nlincorr;
	int flatcorr;
	int unitcorr;
	int photcorr;
	int crcorr;
	int rptcorr;

} IR_Switch;


#endif /* INCL_WF3CORR_H */
