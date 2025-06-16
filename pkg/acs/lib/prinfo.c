/* This file contains:
	PrBegin (label)
	PrEnd (label)
	PrFileName (label, filename)
	PrHdrInfo (obsmode, aperture, opt_elem, detector)
	PrSwitch (keyword, value)
	PrRefInfo (keyword, filename, pedigree, descrip, descrip2)
	PrGrpBegin (label, n)
	PrGrpEnd (label, n)

arguments:

char *label       e.g. "input", "output", "imset", "order"
char *filename    name of input or output file name
char *obsmode     value of OBSMODE keyword:  "ACCUM", etc
char *aperture    value of APERTURE keyword
char *opt_elem    value of OPT_ELEM keyword
char *detector    value of DETECTOR keyword
char *keyword     keyword name of calibration switch
char *pedigree    value of PEDIGREE read from header or from table row
char *descrip     value of DESCRIP read from file or table header
char *descrip2    value of DESCRIP read from table row
int   n           imset number (ignored if <= 0)


	17-Nov-1998 WJH - Revised to support trailer file comments
*/

# include <stdio.h>
# include <ctype.h>
# include <string.h>
# include <time.h>
#include "hstcal.h"
# include "acs.h"
# include "acsversion.h"		/* ACS_CAL_VER */
# include "trlbuf.h"

/* The beginning string will be padded to this many characters, plus one
   to ensure that there's at least one separator.
*/
# define PAD_SIZE       8

# define SCRATCH_SIZE  20	/* length of scratch string */
# define TIME_SIZE     82	/* length of string for date & time */

static char buf[SCRATCH_SIZE];
static char timebuf[TIME_SIZE];		/* string for date & time */

/* Print a message at the beginning of calacs task. */

void PrBegin (char *label) {

	char *GetDateTime (void);

	trlmessage("\n"
                   "%s*** %s -- Version %s ***\n"
                   "Begin    %s",
                   TRL_PREFIX, label, ACS_CAL_VER, GetDateTime());
}

/* Print a message at the end of calacs task. */

void PrEnd (char *label) {

	char *GetDateTime (void);

	trlmessage("End      %s\n\n"
                   "*** %s complete ***", GetDateTime(), label);
}

/* Print input or output file name. */

void PrFileName (char *label, char *filename) {

/* arguments:
char *label     i: e.g. "Input" or "Output"
char *filename  i: name of input or output file
*/
	void Upper1 (char *, char *, int);

	Upper1 (label, buf, SCRATCH_SIZE);

	trlmessage("%s %s", buf, filename);
}

/* Print info from primary header. */

void PrHdrInfo (char *aperture, char *filter1, char *filter2, char *detector) {

	trlmessage("APERTURE %s\n"
                   "FILTER1 %s\n"
                   "FILTER2 %s\n"
                   "DETECTOR %s",
                   aperture, filter1, filter2, detector);
}

/* Print a calibration switch keyword and value. */

void PrSwitch (char *keyword, int value) {

/* arguments:
char *keyword  i: keyword name of calibration switch
int value      i: value of switch (OMIT, PERFORM, etc)
*/

	void UpperAll (char *, char *, int);

	UpperAll (keyword, buf, SCRATCH_SIZE);

	const char *value_s = NULL;
	if (value == OMIT)
	    value_s = "OMIT";
	else if (value == PERFORM)
	    value_s = "PERFORM";
	else if (value == DUMMY)
	    value_s = "SKIPPED";
	else if (value == SKIPPED)
	    value_s = "SKIPPED";
	else if (value == COMPLETE)
	    value_s = "COMPLETE";
	else
	    value_s = "unknown";

	trlmessage ("%s %s", buf, value_s);
}

/* Print a message at the beginning of an imset or spectral order. */

void PrGrpBegin (char *label, int n) {

/* arguments:
char *label  i: to be printed at beginning of line (e.g. "Imset" or "Order")
int n        i: number to be printed following label, if n > 0
*/
	char *GetTime (void);
	void Upper1 (char *, char *, int);

	if (strlen (label) > PAD_SIZE){
		/* The label was longer than the buffer! */
	    strcpy (buf, "*****");
   	} else if (n > 0)
	    snprintf (buf, sizeof(buf), "%s %d", label, n);
	else
	    snprintf (buf, sizeof(buf), "%s", label);

	Upper1 (buf, buf, SCRATCH_SIZE);

	/* Is PAD_SIZE too small for the label and number? */
	if (strlen (buf) > PAD_SIZE && n > 0) {
	    snprintf (buf, sizeof(buf), "%s%d", label, n);	/* delete the space */
	    Upper1 (buf, buf, SCRATCH_SIZE);
	}

	trlmessage("%s Begin %s", buf, GetTime());
}

/* Print a message at the end of an imset or spectral order. */

void PrGrpEnd (char *label, int n) {

/* arguments:
char *label  i: to be printed at beginning of line (e.g. "Imset" or "Order")
int n        i: number to be printed following label, if n > 0
*/
	char *GetTime (void);
	void Upper1 (char *, char *, int);

	if (strlen (label) > PAD_SIZE)
	    strcpy (buf, "*****");
	else if (n > 0)
	    snprintf (buf, sizeof(buf), "%s %d", label, n);
	else
	    snprintf (buf, sizeof(buf), "%s", label);

	Upper1 (buf, buf, SCRATCH_SIZE);

	if (strlen (buf) > PAD_SIZE && n > 0) {
	    snprintf (buf, sizeof(buf), "%s%d", label, n);	/* too long; delete the space */
	    Upper1 (buf, buf, SCRATCH_SIZE);
	}

	trlmessage("%s End %s", buf, GetTime());
}

/* Print reference file keyword and name, plus any pedigree and
   descrip strings that were found.
*/

void PrRefInfo (char *keyword, char *filename,
	char *pedigree, char *descrip, char *descrip2) {

/* arguments:
char *keyword   i: keyword name for reference file
char *filename  i: name of reference file (image or table)
char *pedigree  i: pedigree for image file, table header, or table row
char *descrip   i: first descrip, from image or table header
char *descrip2  i: second descrip, if any, from table row
*/
	void UpperAll (char *, char *, int);

	UpperAll (keyword, buf, SCRATCH_SIZE);

	trlmessage("%s %s", buf, filename);

	if (pedigree[0] != '\0') {
	    trlmessage("%s PEDIGREE=%s", buf, pedigree);
	}

	if (descrip[0] != '\0') {
	    trlmessage("%s DESCRIP =%s", buf, descrip);
	}

	if (descrip2[0] != '\0') {
	    trlmessage("%s DESCRIP =%s", buf, descrip2);
	}
}

/* The remaining routines are used internally. */

/* Format the date and time into the timebuf string, and return that string. */

char *GetDateTime (void) {

	time_t *tp, now;

	tp = NULL;
	now = time (tp);
	strftime (timebuf, TIME_SIZE, "%d-%b-%Y %H:%M:%S %Z", localtime (&now));
	return (timebuf);
}

/* Format the time into the timebuf string, and return that string. */

char *GetTime (void) {

	time_t *tp, now;

	tp = NULL;
	now = time (tp);
	strftime (timebuf, TIME_SIZE, "%H:%M:%S %Z", localtime (&now));
	return (timebuf);
}

/* Convert a string to upper case and pad with blanks. */

void UpperAll (char *instr, char *outstr, int maxch) {

/* arguments:
char *instr   i: input string
char *outstr  o: string converted to upper case
int maxch     i: allocated length of outstr (including EOS)
*/

	int i;

	/* maxch should be larger than PAD_SIZE */
	if (PAD_SIZE >= maxch) {
	    strcpy (outstr, "*****");
	    return;
	}

	for (i = 0;  i < maxch-1;  i++) {
	    if (instr[i] == '\0')
			break;
	    else if (islower (instr[i]))
			outstr[i] = toupper (instr[i]);
	    else
			outstr[i] = instr[i];
	}
	while (i < PAD_SIZE) {
	    outstr[i] = ' ';
	    i++;
	}
	outstr[i] = '\0';
}

/* Convert the first character of a string to upper case, convert the rest
   to lower case, and pad with blanks.
*/

void Upper1 (char *instr, char *outstr, int maxch) {

/* arguments:
char *instr   i: input string
char *outstr  o: string with only the first letter capitalized
int maxch     i: allocated length of outstr (including EOS)
*/

	int i = 0;

	if (PAD_SIZE >= maxch) {
	    strcpy (outstr, "*****");
	    return;
	}

	if (instr[0] != '\0') {
	    if (islower (instr[0]))
			outstr[0] = toupper (instr[0]);
	    else
			outstr[0] = instr[0];

		for (i = 1;  i < maxch-1;  i++) {
			if (instr[i] == '\0')
		    	break;
			else if (isupper (instr[i]))
		    	outstr[i] = tolower (instr[i]);
			else
		    	outstr[i] = instr[i];
	    }
	}
	while (i < PAD_SIZE) {
	    outstr[i] = ' ';
	    i++;
	}
	outstr[i] = '\0';
}
