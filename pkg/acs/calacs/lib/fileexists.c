# include <stdio.h>
#include "hstcal.h"
# include "hstio.h"	/* for ckNewFile */
# include "acs.h"	/* for message output */

/* This routine takes action if the output file already exists.  If not,
   zero is returned.  If it does, the environment variable imclobber is
   checked.  If imclobber is defined and its value is "yes", the file will
   be deleted.  If imclobber is not defined or has some other value, an
   error message will be printed, and an error code will be returned.
*/

int FileExists (char *fname) {

	extern int status;
    int exists;

	int flag;

    /* Initialize exists to 0 */
    exists = EXISTS_NO;
    
	flag = ckNewFile (fname);
	if (flag > 0) {
	    if (flag == 1) {
			sprintf (MsgText, "Output file `%s' already exists.", fname);
			trlwarn (MsgText);
            status = 1021;
            exists = EXISTS_YES;
	    } else {
			sprintf (MsgText, "Can't clobber `%s'.", fname);
			trlerror (MsgText);
            status = 1023;
            exists = EXISTS_YES;
	    }
	}
	return (exists);
}


/* 
	This function explicitly tries to open the file 'trlname' 
	to test whether it already exists and whether it can be 
	opened for adding new comments.

	It returns:
		0 			new trailer file
		1			revising existing file
		1023		couldn't open a file at all
*/

int TrlExists (char *trlname) {

	extern int status;
    int exists;
	
	FILE *fp;

	/* Initialize exists to 0 */
    exists = EXISTS_NO;

	if ( (fp = fopen(trlname,"r") ) == NULL)  {
		/* File does NOT already exist, try to create a new one... */
		if ( (fp = fopen (trlname, "a+") ) == NULL ) {
			/* Can't create new file! */
			sprintf(MsgText,"Can't update/overwrite trailer file '%s'.",trlname);
			trlerror (MsgText);
			
			/* Return error condition, can't create TRL file... */
			status = 1023;	
			return (exists = EXISTS_YES);			
			
		} else {
			/* New Trailer File */
			sprintf (MsgText, "Creating new trailer file `%s'.", trlname);
			trlmessage (MsgText);
			(void)fcloseWithStatus(&fp);
			return (exists);
		}
		
	} else {
		/* File exists already */
		sprintf (MsgText, "Revising existing trailer file `%s'.", trlname);
		trlmessage (MsgText);

		/* This flag is used to set OverwriteMode in TRL files */
		exists = EXISTS_YES;
		(void)fcloseWithStatus(&fp);
		return (exists);		
	}
}
