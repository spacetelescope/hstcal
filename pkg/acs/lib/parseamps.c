
# include <string.h>
# include "acs.h"

/* This function parses the CCDAMP string from the image header and
    returns a new string which only specifies those amps used by
    the CCDCHIP being processed.

24 Jul 01 (WJH): Fixed a problem with parsing the CCDAMP string.
10 Oct 01 (WJH): CCDAMP string was 'strcat'ing entire input instead of
                just 1 character, so used 'strncat' to fix.
*/

void parseWFCamps (char *acsamps, int chip, char *ccdamp) {

/* Parameters:
char *acsamps           i: keyword from image header
int chip                i: value of CCDCHIP from image header
char *ccdamp            o: string with amps used by chip
*/

	int i,j,k;
    int out_max;
    char wfcamps[3];

		i = 0;
		k = 0;
        out_max = 2;
        wfcamps[0] = '\n';
        
        /* Set up string of possible amps used with the chip */
        if (chip == 2)
            strcpy (wfcamps,AMPSTR1);
        else 
            strcpy (wfcamps,AMPSTR2);
        
        /* Pick out only those amps actually used... */
		for (j = 0; j < out_max; j++) {
            
			if (strchr (acsamps, *(wfcamps+j)) != NULL) {
				strncat(ccdamp,&wfcamps[j],1);
            }
 		}
}       
