# include <stdio.h>
# include <ctype.h>
# include <string.h>
# include "hstio.h"
# include <time.h>

# include "acs.h"

int parseObsDate (Hdr *phdr, time_t *date) {
    /* This function reads in the DATE-OBS keyword from the header
        object phdr, parses it out into a string structure, then
        converts it to a time in seconds since 1970.
    */

	extern int status;
    char dateobs[ACS_CBUF];
    struct tm stime;

    const char delimiters[] = " /-";
    char *token, *cp;

    int GetKeyStr (Hdr *, char *, int, char *, char *, int);

    if (GetKeyStr (phdr, "DATE-OBS", USE_DEFAULT, "", dateobs, ACS_CBUF))
        	    return (status);

    /* Initialize entries in stime for values not set by DATE-OBS */
    stime.tm_min = 0;
    stime.tm_hour = 0;
    stime.tm_sec = 0;

    /* Parse out values for year, month, and day from DATE-OBS string */
    cp = strdup (dateobs);                /* Make writable copy.  */
    token = strtok (cp, delimiters);      /* token => year */
    stime.tm_year = atoi(token) - 1900;
    token = strtok (NULL, delimiters);    /* token => month */
    stime.tm_mon = atoi(token);
    token = strtok (NULL, delimiters);    /* token => day */
    stime.tm_mday = atoi(token);

    /* Convert the stime structure into a real number */
    *date = mktime(&stime);

    /* Clean up */
    free(cp);

    return (status);
}

int parseObsDateVal (char *dateobs, time_t *date) {
    /* This function reads in the DATE-OBS keyword from the header
        object phdr, parses it out into a string structure, then
        converts it to a time in seconds since 1970.
    */

	extern int status;
    struct tm stime;

    const char delimiters[] = " /-";
    char *token, *cp;


    /* Initialize entries in stime for values not set by DATE-OBS */
    stime.tm_min = 0;
    stime.tm_hour = 0;
    stime.tm_sec = 0;
    stime.tm_isdst = -1;

    /* Parse out values for year, month, and day from DATE-OBS string */
    cp = strdup (dateobs);                /* Make writable copy.  */
    token = strtok (cp, delimiters);      /* token => year */
    stime.tm_year = atoi(token) - 1900;
    token = strtok (NULL, delimiters);    /* token => month */
    stime.tm_mon = atoi(token);
    token = strtok (NULL, delimiters);    /* token => day */
    stime.tm_mday = atoi(token);

    /* Convert the stime structure into a real number */
    *date = mktime(&stime);

    /* Free memory allocated by strdup() */
    free(cp);

    return (status);
}

int parseTabDate (char *date, time_t *dtime) {
    /* This function reads in the date string from the SPOTTAB
        reference table, parses it out into a string structure, then
        converts it to a time in seconds since 1970.
    */

	extern int status;
    struct tm stime;

    const char delimiters[] = " /-";
    char *token, *cp;

    /* Initialize entries in stime for values not set by DATE-OBS */
    stime.tm_min = 0;
    stime.tm_hour = 12;
    stime.tm_sec = 0;

    /* Parse out values for year, month, and day from DATE-OBS string */
    cp = strdup(date);                /* Make writable copy.  */
    token = strtok (cp, delimiters);      /* token => day */
    stime.tm_mday = atoi(token);
    token = strtok (NULL, delimiters);    /* token => month */
    stime.tm_mon = atoi(token);
    token = strtok (NULL, delimiters);    /* token => year */
    stime.tm_year = atoi(token) + 100;

    /* Convert the stime structure into a real number */
    *dtime = mktime(&stime);

    /* Clean up */
    free(cp);

    return (status);
}
