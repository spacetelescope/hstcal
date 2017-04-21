# include <stdio.h>
# include <string.h>		/* for strncmp, strcmp */

# include "hstio.h"
# include "wf3.h"
# include "msg.h"
# include "wf3info.h"
# include "err.h"		/* defines error codes */
# include "cte.h"

static int checkCCD (Hdr *, WF3Info *, int *, int *);
static int checkBiac (Hdr *, WF3Info *, int *);

/* This routine gets the names of reference images and tables from the
   primary header and checks for dummy pedigree for the CTE reference files
   which are used during the CTE correction. The DARKCFILE is used during regular
   dark correction on data which has already been corrected.

   M.Sosey, 2014 August
   Updated to include the new PCTETAB for the CTE uvis correction
   
   M.Sose, 2015 August
   Updated to check BIASC filename based on PCTECORR and added check
   for BIASCORR complete with clean exit, #1216
 */

int GetCTEFlags (WF3Info *wf3, Hdr *phdr) {

	extern int status;

	int missing = 0;	/* true if any calibration file is missing */
	int nsteps = 0;		/* number of calibration steps to perform */

	/* Get the values for the Calibration Switches from the
	 **	header for processing.  */
 
    
	if (GetCTESwitch (wf3, phdr) )
		return(status);
    
	/* Check each reference file that we need. */

	if (checkCCD(phdr, wf3, &missing, &nsteps))
		return (status);

	if (checkBiac (phdr, wf3, &missing))
	    return (status);
                        
	if (missing) {
		return (status = CAL_FILE_MISSING);

	} else if (nsteps < 1) {
		trlwarn("No calibration switch was set to PERFORM, ");
		trlwarn("            or all reference files had PEDIGREE = DUMMY.");
		return (status = NOTHING_TO_DO);
	} else {
		return (status);
	}
}


/* If the detector is the CCD, we need the CTE parameters table for
   PCTECORR.  This routine checks that the table exists.
 */

static int checkCCD (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

	/* arguments:
	   Hdr *phdr         i: primary header
	   WF3Info *wf3      i: switches, file names, etc
	   int *missing     io: incremented if the table is missing
	   int *nsteps      io: incremented if this step can be performed
	 */

	extern int status;

	int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	void MissingFile (char *, char *, int *);
	void CheckTabType (RefTab *, char *, char *, int *);
    int GotFileName(char *);
    
	if (wf3->pctecorr == PERFORM) {        
		if (GetTabRef (wf3->refnames, phdr, "PCTETAB", &wf3->pctetab, &wf3->pctecorr))
			return (status);

		if (wf3->pctetab.exists != EXISTS_YES) {
			if (GotFileName (wf3->pctetab.name)) {
				MissingFile ("PCTETAB", wf3->pctetab.name, missing);
			}

		} else {

			/* Is the FILETYPE appropriate for a CTE table? */
			CheckTabType (&wf3->pctetab, "PIXCTE", "PCTETAB", missing);
		}

		if (wf3->pctecorr == PERFORM)
			(*nsteps)++;
	}

	return (status);
}
static int checkBiac (Hdr *phdr, WF3Info *wf3, int *missing) {

/* arguments:
Hdr *phdr         i: primary header
WF3Info *wf3      i: switches, file names, etc
int *missing     io: incremented if the file is missing
*/

	extern int status;

	int calswitch;
	int GetSwitch (Hdr *, char *, int *);
	int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
	void MissingFile (char *, char *, int *);
	void CheckImgType (RefImage *, char *, char *, int *);
	int CheckGain (char *, float, char *, int *);
    int GotFileName(char *);
    
	/* Are we supposed to do this step? */
	if (wf3->pctecorr == PERFORM) {

	    if (GetImageRef (wf3->refnames, phdr, "BIACFILE", &wf3->biac, 
			     &wf3->biascorr))
    		return (status=CAL_FILE_MISSING);
	            
        if (wf3->biac.exists != EXISTS_YES) {
    		if (GotFileName (wf3->biac.name)) {
                MissingFile ("BIACFILE", wf3->biac.name, missing);
                return (status=CAL_FILE_MISSING);
            }
	    } else {

		    /* Is the FILETYPE appropriate for a BIAS file? */
		    CheckImgType (&wf3->biac, "CTEBIAS", "BIACFILE", missing);
        }

        /*check if the user has already performed biascorr on the data*/
	    if (GetSwitch (phdr, "BIASCORR", &calswitch))
		    return (status);
	    if (calswitch == COMPLETE) {
            sprintf(MsgText,"\n *** BIASCORR already performed, stopping, cannot run PCTECORR with BIASCORR == COMPLETE");
            trlmessage(MsgText);
		    return (status=ERROR_RETURN);
	    }

    }
	return (status);
}

