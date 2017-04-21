# include <ctype.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "hstio.h"	/* defines HST I/O functions */

# include "wf3.h"	/* defines WF3 data structures */
# include "wf3asn.h"	/* defines WF3 Association data structures */
# include "hstcalerr.h"
# include "calwf3.h"	/* defines WF3 observation data structures */

/* GETASNMEMBER: Copy information from association table structure
**	to a single image structure for use by the remainder of the 
**	processing tasks.
*/
int GetAsnMember (AsnInfo *asn, int prodid, int posid, int expid, WF3Info *wf3){

/* arguments:
AsnInfo *asn      	i: calibration flags and other info
int prodid		i: product id for exposure
int posid		i: sub-product id for exposure 
int expid		i: id of exposure within sub-product
WF3Info *wf3		o: exposure specific flags and info
*/
	extern int status;
	
	char rootname[SZ_FNAME+1];
	char outroot[SZ_CBUF+1];
	char mtype[SZ_FITS_VAL+1];
	int  mlen;
	void FindAsnRoot (char *, char *);
	void UpperAll (char *, char *, int);
    int MkName(char *, char *, char *, char *, char *, int);
    
	/* FIND OUT IF THE MEMBER WE WANT EXISTS */	
	if (asn->product[prodid].subprod[posid].exp[expid].name[0] == '\0') {
	    sprintf (MsgText,
		 "Couldn't find exposure %d of sub-product %d for product %d. ",		 expid, posid, prodid);
	    trlerror (MsgText);
	    return (status = NO_GOOD_DATA);
	}

	/* Initialize local strings */
	rootname[0] = '\0';
	outroot[0]  = '\0';
	mtype[0]    = '\0';

	strcpy (outroot, asn->product[prodid].subprod[posid].name);
	
	/* MAKE SURE WE ARE ONLY PASSING A ROOTNAME, AND NOT A FULL FILENAME.*/
	FindAsnRoot (outroot, rootname);
	strcpy (wf3->outroot, rootname);
					
	strcpy (wf3->rootname,
		asn->product[prodid].subprod[posid].exp[expid].name);
	
    strcpy (wf3->crj_root, asn->product[prodid].subprod[posid].spname);
    strcpy (wf3->crc_root, asn->product[prodid].subprod[posid].spname_cte);
        

	/* Make sure we are only passing a rootname, and not a full filename.*/
	FindAsnRoot (wf3->rootname, rootname);
	strcpy (wf3->rootname, rootname);
	
	if (asn->debug) {
	    sprintf (MsgText, "GetAsnMember: Rootname: %s, Output rootname: %s",
		     rootname, outroot);
	    trlmessage (MsgText);
	}

	/* Check to see that this value of rootname is what we really need */
	strcpy (wf3->asn_table, asn->asn_table);	
	strcpy (wf3->rawfile,
		asn->product[prodid].subprod[posid].exp[expid].expname);
	
	/* Set sci_* flags for wf3 */	
	wf3->sci_crcorr  = asn->crcorr;
	wf3->sci_dthcorr = asn->dthcorr;
	wf3->sci_rptcorr = asn->rptcorr;
	wf3->detector    = asn->detector;
	wf3->nimages     = asn->spmems[posid];

	/* Set MemType appropriate for output */
	if (!strncmp (rootname, wf3->asn_table, 8) ) {
	    mlen = strlen(wf3->mtype);
	    UpperAll (wf3->mtype, mtype, mlen);
	    strcpy (mtype, asn->product[prodid].subprod[posid].mtype);
	}
    
	return (status);
}


int GetSingle (AsnInfo *asn, WF3Info *wf3) {

/* arguments:
AsnInfo *asn      	i: calibration flags and other info
WF3Info *wf3		o: exposure specific flags and info
*/
	extern int status;
	char rootname[SZ_FNAME+1];
	char outroot[SZ_CBUF+1];
	void FindAsnRoot (char *, char *);
	
	strcpy(outroot, asn->filename);
	
	/* Make sure we are only passing a rootname, and not a full filename.*/
	FindAsnRoot (outroot, rootname);
	strcpy (wf3->outroot, rootname);
	strcpy (wf3->rootname, rootname);
	
	if (asn->debug) {
	    sprintf (MsgText, "GetSingle: Rootname: %s, Output rootname: %s",
		     rootname, outroot);
	    trlmessage (MsgText);
	}

	/* Check to see that this value of rootname is what we really need. */
	strcpy (wf3->asn_table, asn->asn_table);	
	strcpy (wf3->rawfile, asn->filename);
	
	/* Set sci_* flags for wf3 */	
	wf3->sci_crcorr  = asn->crcorr;
	wf3->sci_dthcorr = asn->dthcorr;
	wf3->sci_rptcorr = asn->rptcorr;
	wf3->detector    = asn->detector;
	wf3->nimages     = 1;

	/* Set MemType appropriate for output */
	wf3->mtype[0] = '\0';

	return (status);
}

