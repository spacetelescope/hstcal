# include <string.h>
# include "acs.h"

/* This function parses the CCDAMP string from the image header and
    returns a new string which only specifies those amps used by
    the CCDCHIP being processed.

24 Jul 01 (WJH): Fixed a problem with parsing the CCDAMP string.
10 Oct 01 (WJH): CCDAMP string was 'strcat'ing entire input instead of
                 just 1 character, so used 'strncat' to fix.
04 Sep 26 (PLL): Refactored function to avoid compiler warning
                 about strncat possible truncation.
*/

void parseWFCamps (char *acsamps, int chip, char *ccdamp) {

/* Parameters:
char *acsamps           i: keyword from image header
int chip                i: value of CCDCHIP from image header
char *ccdamp            o: string with amps used by chip
*/

    char wfcamps[3], curamp[2];

    /* Set up string of possible amps used with the chip */
    if (chip == 2)
        strcpy(wfcamps, AMPSTR1);
    else
        strcpy(wfcamps, AMPSTR2);

    /* Pick out only those amps actually used... */
    for (size_t j = 0; j < strlen(wfcamps); j++) {
        if (strchr(acsamps, wfcamps[j]) != NULL) {
            curamp[0] = wfcamps[j];
            curamp[1] = '\0';
            strcat(ccdamp, curamp);
        }
    }
}
