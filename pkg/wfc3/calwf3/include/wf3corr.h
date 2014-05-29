/* calibration switches for 'wf3ccd' and 'wf3ir' */

/* Revision history:
**
** H.Bushouse	17 June 2002: Removed "statcorr" from CCD_Switch and
**		and IR_Switch (tracking CALACS changes).
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

