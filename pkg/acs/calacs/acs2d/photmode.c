# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

/*
	This function builds the PHOTMODE string for the image header.
	The string will be used as input to 'synphot' for determining
	the values of the photometric keywords.

*/

int PhotMode (ACSInfo *acs2d, SingleGroup *x) {

/* arguments:
ACSInfo *acs    i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; primary header is modified
*/

    extern int status;

    char *photstr;		/* the photmode string */
    char *scratch;
    float cenwave;		/* central wavelength */

    int PutKeyStr (Hdr *, char *, char *, char *);
    int GetKeyFlt (Hdr *, char *, int, float, float *);

	photstr = calloc (ACS_LINE+1, sizeof (char));
	scratch = calloc (ACS_FITS_REC+1, sizeof (char));
	if (photstr == NULL || scratch == NULL)
	    return (status = OUT_OF_MEMORY);

	strcpy (photstr, "ACS");

	/* Add detector type and chip number for WFC*/
	if (acs2d->detector == MAMA_DETECTOR)
	    strcat (photstr, " SBC");
	else if (acs2d->detector == HRC_CCD_DETECTOR)
	    strcat (photstr, " HRC");
	else if (acs2d->detector == WFC_CCD_DETECTOR) {
	    strcat (photstr, " WFC1");
    }

    /* Check to see if HRC Coronographic mode was used */
    if (strncmp ("CORON",acs2d->obstype,5) == 0) {
	    strcat (photstr, " CORON");
    }

	if (GetKeyFlt (x->globalhdr, "LRFWAVE", USE_DEFAULT, -1.0, &cenwave))
	    return (status);

    if (strncmp("N/A",acs2d->filter1,3) != 0) {
	    if (strncmp ("CLEAR", acs2d->filter1,5) != 0) {

            if (strncmp("POL", acs2d->filter1, 3) == 0){
                if (strstr(acs2d->filter1,"UV") != NULL) {
                    sprintf(scratch,"POL_UV");
                } else {
                    sprintf(scratch,"POL_V");
                }
            } else {
                sprintf(scratch, "%s", acs2d->filter1);
            }
	        strcat (photstr, " ");
	        strcat (photstr, scratch);

            /* If filter is a ramp filter, append value of LRFWAVE
            for use in synphot. */
            if (strncmp("FR", acs2d->filter1, 2) == 0){
	            sprintf (scratch, "#%0.4g", cenwave);
	            strcat (photstr, scratch);
            }
	    }
    }
    if (strncmp("N/A",acs2d->filter2,3) != 0) {
	    if (strncmp ("CLEAR", acs2d->filter2,5) != 0) {
            if (strncmp("POL", acs2d->filter2, 3) == 0){
                if (strstr(acs2d->filter2,"UV") != NULL ){
                    sprintf(scratch,"POL_UV");
                } else {
                    sprintf(scratch,"POL_V");
                }
            } else {
                sprintf(scratch, "%s", acs2d->filter2);
            }
	        strcat (photstr, " ");
	        strcat (photstr, scratch);

            /* If filter is a ramp filter, append value of LRFWAVE
            for use in synphot. */
            if (strncmp("FR", acs2d->filter2, 2) == 0){
	            sprintf (scratch, "#%0.4g", cenwave);
	            strcat (photstr, scratch);
            }
	    }
    }
    /* Add 'mjd#' keyword to PHOTMODE string, but only for WFC and HRC */
	if (acs2d->detector != MAMA_DETECTOR) {
        sprintf (scratch, " MJD#%0.4f", acs2d->expstart);
        strcat (photstr,scratch);
    }

    if (acs2d->verbose){
        sprintf(MsgText,"Keyword PHOTMODE built as: %s",photstr);
        trlmessage(MsgText);
    }
	/* Update PHOTMODE in the extension header. */
	if (PutKeyStr (&x->sci.hdr, "PHOTMODE", photstr, "list of components"))
	    return (status);

	free (scratch);
	free (photstr);

	return (status);
}
