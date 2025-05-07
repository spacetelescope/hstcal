# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string.h>
# include <ctype.h>

#include "str_util.h"
#include "trlbuf.h"
#include "timestamp.h"

/* This routine prints the current date and time into a string
   in the format preferred by OPUS, date & time = YYYYDDDHHMMSS.
   Example:

1996212152529-I--------------- Generic Conversion Completed: N2RUSG54B ---------
*/

void TimeStamp(char *message, const char *rootname) {
    /* arguments:
    char *message   i: string to append to date & time info
    char *rootname  i: root name to include in printed string
    */

    char *timestring = NULL; /* string to be printed */
    char *uc_rootname = NULL; /* rootname converted to upper case */
    size_t lenrootname = 0; /* length of rootname string */
    time_t *tp = NULL;
    time_t now;

    if ((timestring = malloc(SZ_TIMESTRING * sizeof(char))) == NULL) {
        trlerror("(TimeStamp) Out of memory.");
        return;
    }

    if (rootname) {
        lenrootname = strlen(rootname);
        if (lenrootname > 0) {
            uc_rootname = malloc((lenrootname + 1) * sizeof(char));
            if (uc_rootname == NULL) {
                trlerror("(TimeStamp) Out of memory.");
                free(timestring);
                return;
            }
            strcpy(uc_rootname, rootname);
            for (size_t i = 0; i < strlen(rootname); i++) {
                if (islower(uc_rootname[i])) {
                    uc_rootname[i] = (char) toupper(uc_rootname[i]);
                }
            }
        }
    }

    tp = NULL;
    now = time(tp);
    strftime(timestring, SZ_TIMESTRING, "%Y%j%H%M%S", localtime(&now));
    snprintf(timestring + strlen(timestring),
        SZ_TIMESTRING - strlen(timestring),
        "-I--------------- %s%s%s ", message, lenrootname > 0 ? ": " : "", lenrootname > 0 ? uc_rootname : "");

    /* Fill out to 80 bytes with dashes. */
    repchar_s('-', SZ_EIGHTY, timestring + strlen(timestring), SZ_TIMESTRING - strlen(timestring));

    trlmessage("%s", timestring);

    free(uc_rootname);
    free(timestring);
}
