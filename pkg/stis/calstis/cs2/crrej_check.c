# include 	<stdio.h>
# include 	<string.h>
# include 	"hstio.h"

# include	"stis.h"
# include	"cs2.h"
# include	"calstis2.h"

/*  crrej_check -- check input files of crrej

  Description:
  ------------
  Open and check input images/masks to have consistent dimensions and number
  of groups.

  Date		Author		Description
  ----		------		-----------
  05-May-1996  J.-C. Hsu	Adapt from the SPP code crrej_check.x
  10-Feb-2000  Phil Hodge	Replace get_key with getKey[];
				remove iperr from calling sequence;
				move the call to closeImage to inside the loop;
				return int instead of void,
				and replace exit with return.
  13-Aug-2001  Phil Hodge	Check for error after calls to openInputImage.
  24-Oct-2008  Phil Hodge	Add a call to checkImsetOK.
  22-May-2012  Phil Hodge	Change the declaration of imgname;
				initialize ipsci[0] and ipdq[0] to NULL;
				move the call to closeImage to just after the
				call to getHeader.
*/

int crrej_check (IRAFPointer tpin, clpar *par, int newpar[],
			char *imgname[], int grp[],
			IODescPtr ipsci[], IODescPtr ipdq[],
			float noise[], float gain[], int *dim_x, int *dim_y,
			int *nimgs)
{
	IODescPtr	ip;
	Hdr		prihdr;			/* primary header structure */
	char		fdata[STIS_FNAME];
	char		ccdamp[STIS_FITS_REC+1], ccdamp0[STIS_FITS_REC+1];
	float		gn, ron;
	int		k, n, nk;
	int		status, imset_ok;

/* -------------------------------- begin ---------------------------------- */

	initHdr (&prihdr);

	/* Initial values, so we can check (in crrej_do) that images
	   were actually opened.
	*/
	ipsci[0] = NULL;
	ipdq[0] = NULL;

	/* rewind the image template */
	c_imtrew (tpin);
	
	/* loop all input files */
	*nimgs = 0;
	ccdamp0[0] = '\0';
	for (n = 0; n < c_imtlen(tpin); ++n) {

	    /* read the next input image name in the template list */
	    c_imtgetim (tpin, fdata, STIS_FNAME);

	    /* open the primary header */
	    ip = openInputImage (fdata, "", 0);
	    if (hstio_err()) {
		printf ("ERROR    HSTIO error %s\n", hstio_errmsg());
		return (2);
	    }
	
	    /* how many groups in each file */
	    getHeader (ip, &prihdr);
	    closeImage (ip);

	    if (getKeyI (&prihdr, "NEXTEND", &nk) != 0)
		nk = 0;
	    nk /= EXT_PER_GROUP;

	    /* read the CRREJ reference table name from the first file,
		if necessary */
	    if (n == 0) {
		if (newpar[0] < MAX_PAR && par->tbname[0] == '\0') {
	    	    if (getKeyS (&prihdr, "CRREJTAB", par->tbname)) {
			printf ("Keyword CRREJTAB not found in file %s.\n",
				fdata);
			return (2);
		    }
		}

		/* also read the keyword OBSTYPE */
	    	if (getKeyS (&prihdr, "OBSTYPE", par->obstype)) {
		    printf ("Keyword OBSTYPE not found in file %s.\n", fdata);
		    return (2);
		}
	    }

	    /* get the readout noise and A-to-D gain */
	    if (getKeyF (&prihdr, "READNSE", &ron) != 0) {
		printf ("Keyword READNSE not found in file %s.\n", fdata);
		return (2);
	    }
	    if (getKeyF (&prihdr, "ATODGAIN", &gn) != 0) {
		printf ("Keyword ATODGAIN not found in file %s.\n", fdata);
		return (2);
	    }
	    if (gn == 0.) {
		printf ("Keyword ATODGAIN in file %s is 0.\n", fdata);
		return (2);
	    }

	    /* make sure the same CCDAMP is used */
	    if (getKeyS (&prihdr, "CCDAMP", ccdamp) != 0)
		ccdamp[0] = '\0';
	    if (n == 0) {
		strcpy (ccdamp0, ccdamp);
	    } else {
		if (strcmp (ccdamp0, ccdamp) != 0) {
		    printf ("%s uses different CCDAMP.\n", fdata);
		    return (2);
		}
	    }

	    /* loop through the groups */
	    for (k = 1; k <= nk; ++k) {
		/* ignore if IMSET_OK is F */
		if ((status = checkImsetOK (fdata, k, &imset_ok)) != 0) {
		    printf ("ERROR    HSTIO error %s\n", hstio_errmsg());
		    return (2);
		}
		if (!imset_ok)
		    continue;
	    	ipsci[*nimgs] = openInputImage (fdata, "SCI", k);
		if (hstio_err()) {
		    printf ("ERROR    HSTIO error %s\n", hstio_errmsg());
		    return (2);
		}
	    	ipdq[*nimgs]  = openInputImage (fdata, "DQ",  k);
		if (hstio_err()) {
		    printf ("ERROR    HSTIO error %s\n", hstio_errmsg());
		    return (2);
		}
		noise[*nimgs] = ron/gn;
		gain[*nimgs]  = gn;
		sprintf (imgname[*nimgs], "%s", fdata);
		grp[*nimgs] = k;
		(*nimgs)++;
	    }

	    freeHdr (&prihdr);
	}

	/* make sure there is more than one image */
	if (*nimgs < 2) {
	    printf ("Needs more than one input images\n");
	    return (2);
	}
		
	for (k = 0; k < *nimgs; ++k) {

	    /* use the first image's attributes to compare with the rest of
	     the files */
	    if (k == 0) {
		*dim_x = getNaxis1(ipsci[k]);
		*dim_y = getNaxis2(ipsci[k]);
	    }

	    /* verify the image size to be the same as the first image */
	    if (getNaxis1(ipsci[k]) != *dim_x || getNaxis2(ipsci[k]) != *dim_y){
		printf ("file '%s[%d][sci]' does not have the same size as the first image\n", imgname[k], grp[k]);
		return (2);
	    }
	}

	return (0);
}
