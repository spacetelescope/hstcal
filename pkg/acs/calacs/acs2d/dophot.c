# include <stdio.h>
# include <stdlib.h>	/* malloc */
# include <string.h>    /* strchr */

# include "xsynphot.h"		/* for c_phopar */

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "acserr.h"

# define NPHOT      4		/* size of phot returned by c_phopar */
# define ARR_SIZE    10000    /* size of arrays used for throughputs */


/*static void MultFilter (PhotInfo *); */
static void Phot2Obs (char *, char *);

/* This routine gets the absolute flux conversion from PHOTTAB and the
   aperture throughput (filter throughput, really) from APERTAB.  These
   are multiplied together, and then the synphot routine phopar is called
   to determine the inverse sensitivity, reference magnitude (actually a
   constant), pivot wavelength, and RMS bandwidth.  These are written to
   keywords in the primary header.

   Warren Hack, 1998 June 12:
   	Initial ACS version.  
		Changes were made only to GetPhot, except for include file and 
		variable names (sts -> acs2d, StisInfo1 -> ACSInfo).
   Warren Hack, 2002 Jan 30:
    Removed gain normalization of PHOTFLAM keyword value to be consistent
        with image datatype after flat-field calibration. Also, increased
        ARR_SIZE to match the default value used by SYNPHOT.
*/

int doPhot (ACSInfo *acs2d, SingleGroup *x) {

/* arguments:
ACSInfo *acs2d     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; primary header is modified
*/

    extern int status;

    PhotInfo phot;		/* QE and filter throughput, vs wavelength */
    float atodgain;

    float photkey[NPHOT];	/* values of photometry keywords */
    char photmode[ACS_LINE],obsmode[ACS_LINE];
    char graphtbl[ACS_LINE],comptbl[ACS_LINE];
    char grafpath[1025];

    int GetPhot (ACSInfo *, PhotInfo *);
    int GetApThr (ACSInfo *, PhotInfo *);
    int GetKeyStr (Hdr *, char *, int, char *, char *, int);
    int PutKeyFlt (Hdr *, char *, float, char *);

    if (acs2d->photcorr == DUMMY)
        return (status);
        
    /* Extract photmode from sci extension header*/ 
    if (GetKeyStr(&x->sci.hdr,"PHOTMODE",USE_DEFAULT,"",photmode, ACS_LINE))
        return (status);

    /* Add commas to PHOTMODE string and change all strings to
        lower case for use in synphot. */
    Phot2Obs(photmode,obsmode);
    if (acs2d->verbose){
        sprintf(MsgText,"Created SYNPHOT obsmode of: %s",obsmode);
        trlmessage(MsgText);
    }

    /* Initialize photometry keyword structures */
    phot.nelem = ARR_SIZE;
    phot.f_wl = calloc (phot.nelem, sizeof(float));
    phot.f_thru = calloc (phot.nelem, sizeof(float));
    phot.f_err = calloc (phot.nelem, sizeof(float));

    if (phot.f_wl == NULL || phot.f_thru == NULL || phot.f_err == NULL)
        return (status = OUT_OF_MEMORY);

    /* Comment out synphot calls so that the linker does not 
       try to look for these functions at all.
        
      Pass obsmode string to synphot to generate throughput array *
    c_getbandx (obsmode, acs2d->graph.name, acs2d->comp.name, True, 
        phot.nelem, phot.f_wl, phot.f_thru, phot.f_err);
    if (c_iraferr()) {
        trlerror ("Error return from c_getbandx");
        return (status = ERROR_RETURN);
    }

    * Compute photometry values. *
    c_phopar (phot.nelem, phot.f_wl, phot.f_thru, photkey);
    if (c_iraferr()) {
        trlerror ("Error return from c_phopar");
        return (status = ERROR_RETURN);
    }
    c_listpath(obsmode,acs2d->graph.name,acs2d->comp.name,grafpath,1025);
    trlmessage(grafpath);
    */
    /* Add this information as a HISTORY comment */
    addHistoryKw(&x->sci.hdr,grafpath);
    if (hstio_err())
        return (status = HEADER_PROBLEM);

    if (acs2d->verbose){
        sprintf(MsgText,"Computed PHOTFLAM value of %g",photkey[0]);
        trlmessage(MsgText);
    }

    /* Update the photometry keyword values in the SCI header. 
        Revised 10 March 1999 WJH
    */

    if (PutKeyFlt (&x->sci.hdr, "PHOTFLAM", photkey[0], "inverse sensitivity"))
        return (status);

    if (PutKeyFlt (&x->sci.hdr, "PHOTZPT", photkey[1], "zero point"))
        return (status);

    if (PutKeyFlt (&x->sci.hdr, "PHOTPLAM", photkey[2], "pivot wavelength"))
        return (status);

    if (PutKeyFlt (&x->sci.hdr, "PHOTBW", photkey[3], "RMS bandwidth"))
        return (status);

    free (phot.f_wl);
    free (phot.f_thru);
    free (phot.f_err);

    return (status);
}

/* This routine interpolates the aperture throughput (APERTAB) to the
   same wavelengths as the QE (PHOTTAB), and then multiplies the
   QE throughput in-place by the aperture throughput.
*/

static void MultFilter (PhotInfo *phot) {

	double wl;			/* a wavelength in QE array */
	double filt_throughput;		/* interpolated filter throughput */
	int filt_starti;		/* index to begin search in interp1d */
	int i;
	double interp1d (double, double *, double *, int, int *);

	filt_starti = 1;			/* initial value */

/*
	for (i = 0;  i < phot->p_nelem;  i++) {

	    wl = phot->p_wl[i];

	    filt_throughput = interp1d (wl, phot->f_wl, phot->f_thru,
			phot->f_nelem, &filt_starti);

	    phot->p_thru[i] *= filt_throughput;
	}
*/
}

/*
This function converts the PHOTMODE string into an OBSMODE
string suitable for use with synphot functions.
    PHOTMODE - all upper case component names separated by blanks
    OBSMODE - all lower case names separated by commas
NOTE: The operation occurs **in-place** on the photmode string.
*/
static void Phot2Obs(char *photmode, char *obsmode)
{
    char blank = 32, comma = 44;
    int i, len, c1;
        
    len = strlen(photmode);
    
    for (i=0; i < len; i++) {
        c1 = photmode[i];
        if (c1 == blank) {
            obsmode[i] = comma;
        } else 
            obsmode[i] = tolower(c1);
    }
    obsmode[len] = '\0';
}
