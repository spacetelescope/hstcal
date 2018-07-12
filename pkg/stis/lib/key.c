# include <stdio.h>
# include <ctype.h>		/* for isspace */
# include <string.h>		/* for strcpy, strlen */
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"
# include "stisdef.h"

static void KeyMissing (char *);	/* prints error message */
static void AddKeywordMessage (char *);	/* prints warning message */

/* This file contains a set of routines to get keyword values from
   a header and routines to put keyword values into a header.
   The get routines allow specifying a default value which will be
   returned if the keyword is not found.  The put routines will add
   a keyword and its comment if the keyword doesn't already exist.
   In all cases the returned value is status.

	double:   Get_KeyD   Put_KeyD
	float:    Get_KeyF   Put_KeyF
	int:      Get_KeyI   Put_KeyI
	boolean:  Get_KeyB   Put_KeyB
	string:   Get_KeyS   Put_KeyS

   Phil Hodge, 1998 Jan 15:
	Print a warning when adding a new keyword.

   Phil Hodge, 1998 Jan 23:
	Add Get_KeyB and Put_KeyB.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to HEADER_PROBLEM.

   Phil Hodge, 2000 Feb 9:
	Put_KeyB was erroneously defined to be of type Bool; change it to int.
*/

/* Get routines. */

/* double */

int Get_KeyD (Hdr *hd, char *keyword,
	int use_def, double def, double *value) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
double def        i: value to be returned if keyword not found
double *value     o: value gotten
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    if (use_def) {
		*value = def;
	    } else {
		KeyMissing (keyword);
		return (KEYWORD_MISSING);
	    }
	} else {
	    *value = getDoubleKw (key);
	    if (hstio_err()) {
		KeyMissing (keyword);
		return (HEADER_PROBLEM);
	    }
	}

	return (0);
}

/* float */

int Get_KeyF (Hdr *hd, char *keyword,
	int use_def, float def, float *value) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
float def         i: value to be returned if keyword not found
float *value      o: value gotten
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    if (use_def) {
		*value = def;
	    } else {
		KeyMissing (keyword);
		return (KEYWORD_MISSING);
	    }
	} else {
	    *value = getFloatKw (key);
	    if (hstio_err()) {
		KeyMissing (keyword);
		return (HEADER_PROBLEM);
	    }
	}

	return (0);
}

/* int */

int Get_KeyI (Hdr *hd, char *keyword,
	int use_def, int def, int *value) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
int def           i: value to be returned if keyword not found
int *value        o: value gotten
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    if (use_def) {
		*value = def;
	    } else {
		KeyMissing (keyword);
		return (KEYWORD_MISSING);
	    }
	} else {
	    *value = getIntKw (key);
	    if (hstio_err()) {
		KeyMissing (keyword);
		return (HEADER_PROBLEM);
	    }
	}

	return (0);
}

/* boolean */

int Get_KeyB (Hdr *hd, char *keyword,
	int use_def, Bool def, Bool *value) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
Bool def          i: value to be returned if keyword not found
Bool *value       o: value gotten
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    if (use_def) {
		*value = def;
	    } else {
		KeyMissing (keyword);
		return (KEYWORD_MISSING);
	    }
	} else {
	    *value = getBoolKw (key);
	    if (hstio_err()) {
		KeyMissing (keyword);
		return (HEADER_PROBLEM);
	    }
	}

	return (0);
}

/* string */

int Get_KeyS (Hdr *hd, char *keyword,
	int use_def, char *def, char *value, int maxch) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
char *def         i: value to be returned if keyword not found
char *value       o: value gotten
int maxch         i: allocated size of string
*/

	FitsKw key;		/* location of keyword in header */
	int i;

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    if (use_def) {
		strcpy (value, def);
	    } else {
		KeyMissing (keyword);
		return (KEYWORD_MISSING);
	    }
	} else {
	    getStringKw (key, value, maxch);
	    if (hstio_err()) {
		KeyMissing (keyword);
		return (HEADER_PROBLEM);
	    }
	}

	for (i = strlen (value) - 1;  i >= 0;  i--) {
	    if (isspace (value[i]))
		value[i] = '\0';
	    else
		break;
	}

	return (0);
}

/* Put routines. */

/* double */

int Put_KeyD (Hdr *hd, char *keyword, double value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
double value      i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    AddKeywordMessage (keyword);	/* print a warning */
	    addDoubleKw (hd, keyword, value, comment);
	} else {
	    putDoubleKw (key, value);
	}

	if (hstio_err())
	    return (HEADER_PROBLEM);

	return (0);
}

/* float */

int Put_KeyF (Hdr *hd, char *keyword, float value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
float value       i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    AddKeywordMessage (keyword);
	    addFloatKw (hd, keyword, value, comment);
	} else {
	    putFloatKw (key, value);
	}

	if (hstio_err())
	    return (HEADER_PROBLEM);

	return (0);
}

/* int */

int Put_KeyI (Hdr *hd, char *keyword, int value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int value         i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    AddKeywordMessage (keyword);
	    addIntKw (hd, keyword, value, comment);
	} else {
	    putIntKw (key, value);
	}

	if (hstio_err())
	    return (HEADER_PROBLEM);

	return (0);
}

/* boolean */

int Put_KeyB (Hdr *hd, char *keyword, Bool value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
Bool value        i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    AddKeywordMessage (keyword);
	    addBoolKw (hd, keyword, value, comment);
	} else {
	    putBoolKw (key, value);
	}

	if (hstio_err())
	    return (HEADER_PROBLEM);

	return (0);
}

/* string */

int Put_KeyS (Hdr *hd, char *keyword, char *value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
char *value       i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    AddKeywordMessage (keyword);
	    addStringKw (hd, keyword, value, comment);
	} else {
	    putStringKw (key, value);
	}

	if (hstio_err())
	    return (HEADER_PROBLEM);

	return (0);
}

static void KeyMissing (char *keyword) {

	printf ("ERROR    Keyword = `%s'.\n", keyword);
}

static void AddKeywordMessage (char *keyword) {

	printf ("Warning  Keyword `%s' is being added to header.\n", keyword);
}
