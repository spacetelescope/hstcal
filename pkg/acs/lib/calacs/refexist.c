# include <stdio.h>
# include <stdlib.h>

#include "hstcal.h"
# include "hstio.h"
# include "ximio.h"
# include "xtables.h"
# include "c_iraf.h"	/* clear_cvoserr */

# include "acs.h"
# include "calacs.h"

static int SpecialCheck (char *, char *, int, int *, int *);
static void MissingFile (char *, char *, int *);
int DoesFileExist (char *);

/* This function checks for the existence of the reference files that
   we need.  There are a couple of cases where a reference file need not
   be specified (dqitab, if CCD, and flat fields), but if a file name
   was given, we require that the file exist.

   This function may be called more than once, e.g. for science and
   wavecal RefFileInfo structs.  Note that InitNames must have been
   called prior to calling this routine, as SaveName is called.

   The keywords that are given special treatment are WBIAFILE, BPIXTAB,
   PFLTFILE, DFLTFILE, and LFLTFILE.

   Phil Hodge, 1997 Dec 10:
	File created.

   Phil Hodge, 1998 Jan 13:
	Include PCTAB and APERTAB as "special" files.

   Warren Hack, 1998 May 26:
    Revised to work with ACS.
    
   Warren Hack, 2001 July 3:
        Revised to report more informative error messages for
        situations where a NULL string is given as a ref file name.
*/

void RefExist (RefFileInfo *ref, int detector, int *missing) {

/* arguments:
RefFileInfo *ref   i: the list of reference info
int detector       i: detector code, to check for CCD
int *missing       o: a count of missing (or blank) reference files
*/

	RefFileInfo *current;
	int exist;			/* EXISTS_YES or EXISTS_NO */
	int flat_not_specified = 0;
	int streq_ic (char *, char *);	/* strings equal? (case insensitive) */

	current = ref->next;		/* skip dummy first record */
	while (current != NULL) {
		/* Handle missing filenames, rather than searching for
			the file 'N/A'...
		*/        
		if (streq_ic (current->filename, "N/A") || streq_ic (current->filename, "") ){ 
			exist = EXISTS_NO;
            
		} else {
	    	/* Check if file exists. */
	    	exist = DoesFileExist (current->filename);
		}

	    if (exist != EXISTS_YES) {

		/*
		   If this is a "special" keyword, treat it individually;
		   otherwise, log it as missing.
		*/
		if (!SpecialCheck (current->keyword, current->filename,
				detector, &flat_not_specified, missing)) {
            MissingFile (current->keyword, current->filename, missing);
            }
	    }

	    current = current->next;
	}

	/* All three flats not specified?  That's not OK. */
	if (flat_not_specified == 3) {
	    (*missing)++;
	    trlwarn ("PFLTFILE, DFLTFILE, LFLTFILE are all blank");
	}
}

static void MissingFile (char *keyword, char *filename, int *missing) {

	sprintf (MsgText, "%s `%s' not found or can't open.", keyword, filename);
	trlwarn (MsgText);
	(*missing)++;
}

/* This routine checks whether a file exists by trying to open it
   read-only.  If it does exist, the macro EXISTS_YES will be returned
   as the function value; if not, EXISTS_NO will be returned.  
*/

int DoesFileExist (char *filename) {

	int flag;
	IODescPtr im;		/* descriptor for reference image */
	IRAFPointer tp;		/* for checking if file exists */

		im = openInputImage (filename, "", 0);
		if (hstio_err()) {
		    clear_cvoserr();
		    clear_hstioerr();
		    tp = c_tbtopn (filename, IRAF_READ_ONLY, 0);
		    if (c_iraferr()) {
				flag = EXISTS_NO;
				clear_cvoserr();
		    } else {
				flag = EXISTS_YES;
				c_tbtclo (tp);
		    }
		} else {
		    flag = EXISTS_YES;
		    closeImage (im);
		}

	return (flag);
}

/* This routine checks whether the current keyword is "special".  If it is,
   it is handled individually.

   The function value will be one if it is a special keyword, and it will
   be zero otherwise.
*/

static int SpecialCheck (char *keyword, char *filename,
		int detector,
		int *flat_not_specified, int *missing) {

	int flag = 0;		/* to be returned as the function value */
	int streq_ic (char *, char *);	/* strings equal? (case insensitive) */

    if (streq_ic (filename, "")){
        /* 
            If a blank filename is specified, always flag an error. 
        */
		MissingFile (keyword, filename, missing); 
       
	} else if (streq_ic (keyword, "BPIXTAB")) {

	    /* We could do DQICORR just to check saturation,
		if the detector is the CCD.
	    */
	    if (detector == MAMA_DETECTOR )
		    MissingFile (keyword, filename, missing);
	    flag = 1;

	} else if (streq_ic (keyword, "PFLTFILE") ||
		   streq_ic (keyword, "DFLTFILE") ||
		   streq_ic (keyword, "LFLTFILE")) {

	    /* It's OK for some flats to not be specified, but keep count. */
	    (*flat_not_specified)++;
	    flag = 1;

	} else if (streq_ic (keyword, "APERTAB")) {

	    /* It's OK for this to not be specified only for PHOTCORR. */
	    flag = 1;
	}

	return (flag);
}
