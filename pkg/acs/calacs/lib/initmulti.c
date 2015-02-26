# include <string.h>
# include	"acs.h"

/* Initialize multiamp structure to ZERO */
void initmulti (multiamp *amp) {
	int i;
	
	for (i=0; i < NAMPS; i++) 
		amp->val[i] = 0.;
		
	amp->colx = 0;
	amp->coly = 0;
	amp->chip = 0;
	amp->detector = UNKNOWN_DETECTOR;
}
