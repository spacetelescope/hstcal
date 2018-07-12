# include <ctype.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#include "hstcal.h"
# include "hstio.h"	/* defines HST I/O functions */
# include "acs.h"	/* defines ACS data structures */
# include "acsasn.h"	/* defines ACS Association data structures */
# include "hstcalerr.h"
# include "calacs.h"	/* defines ACS observation data structures */

/* GETASNMEMBER: Copy information from association table structure
**	to a single image structure for use by the remainder of the 
**	processing tasks.
*/
int GetAsnMember (AsnInfo *asn, int prodid, int posid, int expid, ACSInfo *acs) {

/* arguments:
AsnInfo *asn      	i: calibration flags and other info
int prodid			i: product id for exposure
int posid			i: sub-product id for exposure 
int expid			i: id of exposure within sub-product
ACSInfo *acs		o: exposure specific flags and info
*/
	extern int status;
	
	char rootname[CHAR_FNAME_LENGTH+1];
	char mtype[SZ_STRKWVAL+1];
    int mlen;
	void FindAsnRoot (const char *, char *);
    void UpperAll (char *, char *, int);
    
	/* find out if the member we want exists... */	
	if (asn->product[prodid].subprod[posid].exp[expid].name[0] == '\0') {
		sprintf(MsgText,"Couldn't find exposure %d of sub-product %d for product %d. ", expid, posid, prodid);
		trlerror (MsgText);
		return (status = NO_GOOD_DATA);
	}

	/* Initialize local strings... */
	rootname[0] = '\0';
    mtype[0] = '\0';

    const char * outroot = asn->product[prodid].subprod[posid].name;
	
	/* Make sure we are only passing a rootname, and not a full filename.*/
	FindAsnRoot (outroot, rootname);
	strcpy (acs->outroot, rootname);
					
	strcpy (acs->rootname, asn->product[prodid].subprod[posid].exp[expid].name);
	strcpy (acs->crj_root, asn->product[prodid].subprod[posid].spname);

	/* Make sure we are only passing a rootname, and not a full filename.*/
	FindAsnRoot (acs->rootname, rootname);
	strcpy (acs->rootname, rootname);
	
	if (asn->debug){
		sprintf (MsgText, "GetAsnMember: Rootname: %s, Output rootname: %s",rootname, outroot);
		trlmessage (MsgText);
	}
	/* Check to see that this value of rootname is what we really need... */
	strcpy (acs->asn_table, asn->asn_table);	
	strcpy (acs->rawfile, asn->product[prodid].subprod[posid].exp[expid].expname);
	
	/* Set sci_* flags for acs */	
	acs->sci_crcorr = asn->crcorr;
	acs->sci_dthcorr = asn->dthcorr;
	acs->sci_rptcorr = asn->rptcorr;
	acs->detector = asn->detector;
	acs->nimages = asn->spmems[posid];

    /* Set MemType appropriate for output */
    if (!strncmp (rootname, acs->asn_table, 8) ) {
        mlen = strlen(acs->mtype);
        UpperAll (acs->mtype, mtype, mlen);
        strcpy(mtype, asn->product[prodid].subprod[posid].mtype);
    }
    
	return (status);
}


int GetSingle (AsnInfo *asn, ACSInfo *acs) {

/* arguments:
AsnInfo *asn      	i: calibration flags and other info
ACSInfo *acs		o: exposure specific flags and info
*/
	extern int status;
	char rootname[CHAR_FNAME_LENGTH];
	*rootname = '\0';
	void FindAsnRoot (const char *, char *);
	
	const char * outroot = asn->filename;
	
	/* Make sure we are only passing a rootname, and not a full filename.*/
	FindAsnRoot (outroot, rootname);
	strcpy (acs->outroot, rootname);
	strcpy (acs->rootname, rootname);
	
	if (asn->debug) {
		sprintf (MsgText, "GetSingle: Rootname: %s, Output rootname: %s",rootname, outroot);
		trlmessage (MsgText);
	}

	/* Check to see that this value of rootname is what we really need... */
	strcpy (acs->asn_table, asn->asn_table);	
	strcpy (acs->rawfile, asn->filename);
	
	/* Set sci_* flags for acs */	
	acs->sci_crcorr = asn->crcorr;
	acs->sci_dthcorr = asn->dthcorr;
	acs->sci_rptcorr = asn->rptcorr;
	acs->detector = asn->detector;
	acs->nimages = 1;

    /* Set MemType appropriate for output */
    acs->mtype[0] = '\0';

	return (status);
}
