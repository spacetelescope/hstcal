/* calibration switches for 'acsccd' and 'acs2d' */

typedef struct {
	int dqicorr;
	int atodcorr;
	int blevcorr;
	int biascorr;
	int flashcorr;
  int pctecorr;
} acsccd_switch;

typedef struct {
	int dqicorr;
	int glincorr;
	int lflgcorr;
	int darkcorr;
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
	int dqicorr;
	int flatcorr;
	int flashcorr;
  int pctecorr;
	int glincorr;
	int lflgcorr;
	int photcorr;
	int rptcorr;
	int shadcorr;
	int expscorr;

} CalSwitch;
