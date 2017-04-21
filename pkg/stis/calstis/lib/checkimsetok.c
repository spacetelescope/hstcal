# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "stisdef.h"
# include "err.h"

/* This routine gets keyword IMSET_OK from an extension header and
   returns the value of this keyword (default is True).  A value of
   1 means the specified imset is OK, 0 means the imset either has
   zero exposure time or the raw data values are constant.

   The function value will be 0 if the image set exists and there was
   no error opening it or reading the keyword (keyword missing is not
   considered an error in this case).  The function value will be -1
   if the file exists (primary header can be opened) but the image set
   (identified by extver) does not.  If there was another error, the
   hstio status will be returned.

   Phil Hodge, 2008 Dec 12:
	Function created.

   Phil Hodge, 2010 Mar 9:
	Calls to c_tbtclo were added.

   Phil Hodge, 2011 Apr 12:
	Delete one of the calls to c_tbtclo, the one in the section for
	when c_iraferr() was non-zero after calling c_tbtopn().
*/

int checkImsetOK (char *input, int extver, int *imset_ok) {

/* arguments:
char *input      i: name of input file
int extver       i: imset number
int *imset_ok    o: value of IMSET_OK keyword
*/

	int status = 0;
	IODescPtr fd;
	Hdr hdr;
	Bool value;		/* value of IMSET_OK keyword */
	int use_def = 1;	/* use default if keyword is missing */
	int len_input;		/* length of input string */
	char *tablename;	/* name including "SCI" and extver */
	IRAFPointer tp;		/* to open the image as a table */

	*imset_ok = 0;		/* default */

	len_input = strlen (input);
	len_input += 10;	/* allow for [SCI,<number>] */
	if ((tablename = (char *) calloc (len_input+1, sizeof(char))) == NULL)
	    return (OUT_OF_MEMORY);
	sprintf (tablename, "%s[SCI,%d]", input, extver);
	tp = c_tbtopn (tablename, IRAF_READ_ONLY, 0);
	free (tablename);
	if (c_iraferr()) {
	    *imset_ok = 0;	/* image set is missing, so not OK */
	    clear_cvoserr();
	    return -1;
	}
	c_tbtclo (tp);

	initHdr (&hdr);
	fd = openInputImage (input, "SCI", extver);
	if ((status = hstio_err()) != 0)
	    return status;
	getHeader (fd, &hdr);
	if ((status = Get_KeyB (&hdr, "IMSET_OK", use_def, True, &value)) != 0)
	    return status;
	if (value)
	    *imset_ok = 1;
	else
	    *imset_ok = 0;
	freeHdr (&hdr);
	closeImage (fd);

	return status;
}
