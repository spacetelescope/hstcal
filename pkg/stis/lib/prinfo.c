/* This file contains:
	PrBegin (csnumber)
	PrEnd (csnumber)
	PrFullVersion()
	PrVersion()
	PrFileName (label, filename)
	PrHdrInfo (obsmode, aperture, opt_elem, detector)
	PrSwitch (keyword, value)
	PrRefInfo (keyword, filename, pedigree, descrip, descrip2)
	PrGrpBegin (label, n)
	PrGrpEnd (label, n)

arguments:

int   csnumber    calstis number (i.e. 0, 1, 2, 4, 6, 7, 8, 11, 12)
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
int   n           imset or spectral order number (ignored if <= 0)
*/

# include <stdio.h>
# include <stdlib.h>
# include <ctype.h>
# include <string.h>
# include <time.h>

# include "stis.h"
# include "hstcalerr.h"
# include "stisversion.h"		/* STIS_CAL_VER */

/* The beginning string will be padded to this many characters, plus one
   to ensure that there's at least one separator.
*/
# define PAD_SIZE       8

# define SCRATCH_SIZE  20	/* length of scratch string */
# define TIME_SIZE     82	/* length of string for date & time */

static char buf[SCRATCH_SIZE];
static char timebuf[TIME_SIZE];		/* string for date & time */

static char *GetDateTime (void);	/* these write into and return buf */
static char *GetTime (void);
static void UpperAll (char *, char *, int);
static void Upper1 (char *, char *, int);

/* Print a message at the beginning of calstis-N.

   Phil Hodge, 2000 Jan 6:
	Call fflush.
*/

void PrBegin (int csnumber) {

	printf ("\n");
	printf ("*** CALSTIS-%d -- Version %s ***\n", csnumber, STIS_CAL_VER);
	printf ("Begin    %s\n", GetDateTime());

	fflush (stdout);
}

/* Print a message at the end of calstis-N. */

void PrEnd (int csnumber) {

	printf ("End      %s\n", GetDateTime());
	printf ("\n");
	printf ("*** CALSTIS-%d complete ***\n", csnumber);

	fflush (stdout);
}

void PrFullVersion (void) {

	printf ("%s\n", STIS_CAL_VER);
	fflush (stdout);
}

void PrVersion (void) {

	int i;
	int len;
	char *vstring;

	len = strlen (STIS_CAL_VER);
	vstring = (char *) calloc (len + 1, sizeof(char));
	if (vstring == NULL)
	    exit (OUT_OF_MEMORY);
	for (i = 0;  i <= len;  i++) {
	    if (STIS_CAL_VER[i] == ' ') {
		vstring[i] = '\0';
		break;
	    }
	    vstring[i] = STIS_CAL_VER[i];
	}
	printf ("%s\n", vstring);
	fflush (stdout);
	free (vstring);
}

/* Print input or output file name. */

void PrFileName (char *label, char *filename) {

/* arguments:
char *label     i: e.g. "Input" or "Output"
char *filename  i: name of input or output file
*/

	Upper1 (label, buf, SCRATCH_SIZE);

	printf ("%s %s\n", buf, filename);

	fflush (stdout);
}

/* Print info from primary header. */

void PrHdrInfo (char *obsmode, char *aperture, char *opt_elem, char *detector) {

	printf ("OBSMODE  %s\n", obsmode);
	printf ("APERTURE %s\n", aperture);
	printf ("OPT_ELEM %s\n", opt_elem);
	printf ("DETECTOR %s\n", detector);

	fflush (stdout);
}

/* Print a calibration switch keyword and value. */

void PrSwitch (char *keyword, int value) {

/* arguments:
char *keyword  i: keyword name of calibration switch
int value      i: value of switch (OMIT, PERFORM, etc)
*/

	UpperAll (keyword, buf, SCRATCH_SIZE);

	printf ("%s", buf);

	if (value == OMIT)
	    printf (" OMIT\n");
	else if (value == PERFORM)
	    printf (" PERFORM\n");
	else if (value == DUMMY)
	    printf (" SKIPPED\n");
	else if (value == SKIPPED)
	    printf (" SKIPPED\n");
	else if (value == COMPLETE)
	    printf (" COMPLETE\n");
	else
	    printf (" unknown\n");

	fflush (stdout);
}

/* Print a message at the beginning of an imset or spectral order. */

void PrGrpBegin (char *label, int n) {

/* arguments:
char *label  i: to be printed at beginning of line (e.g. "Imset" or "Order")
int n        i: number to be printed following label, if n > 0
*/

	if (strlen (label) > PAD_SIZE)
	    strcpy (buf, "*****");
	else if (n > 0)
	    sprintf (buf, "%s %d", label, n);
	else
	    sprintf (buf, "%s", label);

	Upper1 (buf, buf, SCRATCH_SIZE);

	/* Is PAD_SIZE too small for the label and number? */
	if (strlen (buf) > PAD_SIZE && n > 0) {
	    sprintf (buf, "%s%d", label, n);	/* delete the space */
	    Upper1 (buf, buf, SCRATCH_SIZE);
	}

	printf ("%s Begin %s\n", buf, GetTime());

	fflush (stdout);
}

/* Print a message at the end of an imset or spectral order. */

void PrGrpEnd (char *label, int n) {

/* arguments:
char *label  i: to be printed at beginning of line (e.g. "Imset" or "Order")
int n        i: number to be printed following label, if n > 0
*/

	if (strlen (label) > PAD_SIZE)
	    strcpy (buf, "*****");
	else if (n > 0)
	    sprintf (buf, "%s %d", label, n);
	else
	    sprintf (buf, "%s", label);

	Upper1 (buf, buf, SCRATCH_SIZE);

	if (strlen (buf) > PAD_SIZE && n > 0) {
	    sprintf (buf, "%s%d", label, n);	/* too long; delete the space */
	    Upper1 (buf, buf, SCRATCH_SIZE);
	}

	printf ("%s End %s\n", buf, GetTime());

	fflush (stdout);
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

	UpperAll (keyword, buf, SCRATCH_SIZE);

	printf ("%s %s\n", buf, filename);

	if (pedigree[0] != '\0')
	    printf ("%s PEDIGREE=%s\n", buf, pedigree);

	if (descrip[0] != '\0')
	    printf ("%s DESCRIP =%s\n", buf, descrip);

	if (descrip2[0] != '\0')
	    printf ("%s DESCRIP =%s\n", buf, descrip2);

	fflush (stdout);
}

/* The remaining routines are used internally. */

/* Format the date and time into the timebuf string, and return that string. */

static char *GetDateTime (void) {

	time_t *tp, now;

	tp = NULL;
	now = time (tp);
	strftime (timebuf, TIME_SIZE, "%d-%b-%Y %H:%M:%S %Z", localtime (&now));
	return (timebuf);
}

/* Format the time into the timebuf string, and return that string. */

static char *GetTime (void) {

	time_t *tp, now;

	tp = NULL;
	now = time (tp);
	strftime (timebuf, TIME_SIZE, "%H:%M:%S %Z", localtime (&now));
	return (timebuf);
}

/* Convert a string to upper case and pad with blanks. */

static void UpperAll (char *instr, char *outstr, int maxch) {

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

static void Upper1 (char *instr, char *outstr, int maxch) {

/* arguments:
char *instr   i: input string
char *outstr  o: string with only the first letter capitalized
int maxch     i: allocated length of outstr (including EOS)
*/

	int i;

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
