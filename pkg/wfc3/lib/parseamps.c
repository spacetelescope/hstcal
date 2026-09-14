# include <string.h>
# include "wf3.h"

/* This function parses the CCDAMP string from the image header and
    returns a new string which only specifies those amps used by
    the CCDCHIP being processed.

   Howard Bushouse, 2001 Nov 16:
    Fixed a problem with parsing the CCDAMP string.
    CCDAMP string was 'strcat'ing entire input instead of just 1
    character, so used 'strncat' to fix.

   P.L. Lim, 2026 Sep 11:
    Refactored function to avoid compiler warning
    about strncat possible truncation.
*/

void parseWFCamps (char *wf3amps, int chip, char *ccdamp) {

/* Parameters:
char *wf3amps           i: keyword from image header
int chip                i: value of CCDCHIP from image header
char *ccdamp            o: string with amps used by chip
*/

    char wfcamps[3], curamp[2];

    /* Set up string of possible amps used with the chip */
    if (chip == 2)
        strcpy (wfcamps, AMPSTR1);
    else
        strcpy (wfcamps, AMPSTR2);

    /* Pick out only those amps actually used... */
    for (size_t j = 0; j < strlen(wfcamps) ; j++) {
        if (strchr(wf3amps, wfcamps[j]) != NULL) {
            curamp[0] = wfcamps[j];
            curamp[1] = '\0';
            strcat(ccdamp, curamp);
        }
    }
}
