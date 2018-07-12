# include <stdio.h>
# include <ctype.h>                 /* for isspace */
# include <string.h>                /* for strcpy, strlen */
#include "hstcal.h"
# include "hstio.h"
# include "hstcalerr.h"
# include "acs.h"                /* for message output */

static void KeyMissing (char *);    /* prints error message */

/* This file contains a set of routines to get keyword values from
   a header and routines to put keyword values into a header.
   The get routines allow specifying a default value which will be
   returned if the keyword is not found.  The put routines will add
   a keyword and its comment if the keyword doesn't already exist.
   In all cases the returned value is status.

    double:   GetKeyDbl   PutKeyDbl
    float:    GetKeyFlt   PutKeyFlt
    int:      GetKeyInt   PutKeyInt
    boolean:  GetKeyBool   PutKeyBool
    string:   GetKeyStr   PutKeyStr

    Warren Hack, 1998 May 26:
    Revised for ACS

    Warren Hack, 1999 Jan 4:
    Revised the PutKey functions to simply call the appropriate 
        HSTIO functions.  The GetKey functions were unchanged in order
        to retain support for error messages within the trailer files,
        as well as the use of default values for the operation.
        
    Warren Hack, 2000 Jan 12:
        Renamed all the functions to avoid name conflicts with HSTIO 
        functions under VMS.
*/

/* Get routines. */

/* ------------------------------------------------------------------*/
/*                          GetKeyDbl                                  */
/* ------------------------------------------------------------------*/


/* double */

int GetKeyDbl (Hdr *hd, char *keyword,
    int use_def, double def, double *value) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
double def        i: value to be returned if keyword not found
double *value     o: value gotten
*/

    extern int status;
    FitsKw key;		/* location of keyword in header */

    key = findKw (hd, keyword);
    if (key == NotFound) {
        if (use_def) {
            *value = def;
        } else {
            KeyMissing (keyword);
            return (status = KEYWORD_MISSING);
        }
    } else {
        *value = getDoubleKw (key);
        if (hstio_err()) {
            KeyMissing (keyword);
            return (status = HEADER_PROBLEM);
	    }
    }

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          GetKeyFlt                                  */
/* ------------------------------------------------------------------*/


/* float */

int GetKeyFlt (Hdr *hd, char *keyword, int use_def, float def, float *value) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
float def         i: value to be returned if keyword not found
float *value      o: value gotten
*/

    extern int status;
    FitsKw key;		/* location of keyword in header */

    key = findKw (hd, keyword);
    if (key == NotFound) {
        if (use_def) {
            *value = def;
        } else {
            KeyMissing (keyword);
            return (status = KEYWORD_MISSING);
        }
    } else {
        *value = getFloatKw (key);
        if (hstio_err()) {
            KeyMissing (keyword);
            return (status = HEADER_PROBLEM);
        }
    }

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          GetKeyInt                                */
/* ------------------------------------------------------------------*/


/* int */

int GetKeyInt (Hdr *hd, char *keyword,
    int use_def, int def, int *value) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
int def           i: value to be returned if keyword not found
int *value        o: value gotten
*/

    extern int status;
    FitsKw key;		/* location of keyword in header */

    key = findKw (hd, keyword);
    if (key == NotFound) {
        if (use_def) {
            *value = def;
        } else {
            KeyMissing (keyword);
            return (status = KEYWORD_MISSING);
        }
    } else {
        *value = getIntKw (key);
        if (hstio_err()) {
            KeyMissing (keyword);
            return (status = HEADER_PROBLEM);
        }
    }

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          GetKeyBool                                  */
/* ------------------------------------------------------------------*/


/* boolean */

int GetKeyBool (Hdr *hd, char *keyword,
    int use_def, Bool def, Bool *value) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
Bool def          i: value to be returned if keyword not found
Bool *value       o: value gotten
*/

    extern int status;
    FitsKw key;		/* location of keyword in header */

    key = findKw (hd, keyword);
    if (key == NotFound) {
        if (use_def) {
            *value = def;
        } else {
            KeyMissing (keyword);
            return (status = KEYWORD_MISSING);
        }
    } else {
        *value = getBoolKw (key);
        if (hstio_err()) {
            KeyMissing (keyword);
            return (status = HEADER_PROBLEM);
        }
    }

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          GetKeyStr                                  */
/* ------------------------------------------------------------------*/


/* string */

int GetKeyStr (Hdr *hd, char *keyword,
    int use_def, char *def, char *value, int maxch) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int use_def       i: if true and keyword not found, return default
char *def         i: value to be returned if keyword not found
char *value       o: value gotten
int maxch         i: allocated size of string
*/

    extern int status;
    FitsKw key;		/* location of keyword in header */
    int i;

    key = findKw (hd, keyword);
    if (key == NotFound) {
        if (use_def) {
            strcpy (value, def);
        } else {
            KeyMissing (keyword);
            return (status = KEYWORD_MISSING);
        }
    } else {
        getStringKw (key, value, maxch);
        if (hstio_err()) {
            KeyMissing (keyword);
            return (status = HEADER_PROBLEM);
        }
    }

    for (i = strlen (value) - 1;  i >= 0;  i--) {
        if (isspace (value[i]))
            value[i] = '\0';
        else
            break;
    }

    return (status);
}

/* Put routines. */

/* These routines call the HSTIO 'putKey' routines initially added to
	HSTIO version 2.1 and later, which were derived from CALNIC's
	functions.  They had identical arguments but a slightly different
	name, so the changes were confined to this file...
	
	WJH 4Jan99
*/ 

/* ------------------------------------------------------------------*/
/*                          PutKeyDbl                                  */
/* ------------------------------------------------------------------*/


/* double */

int PutKeyDbl (Hdr *hd, char *keyword, double value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
double value      i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

    extern int status;

    status = putKeyD (hd, keyword, value, comment);

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          PutKeyFlt                                  */
/* ------------------------------------------------------------------*/


/* float */

int PutKeyFlt (Hdr *hd, char *keyword, float value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
float value       i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

    extern int status;

    status = putKeyF (hd, keyword, value, comment);

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          PutKeyInt                                  */
/* ------------------------------------------------------------------*/


/* int */

int PutKeyInt (Hdr *hd, char *keyword, int value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
int value         i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

    extern int status;

    status = putKeyI (hd, keyword, value, comment);

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          PutKeyBool                                  */
/* ------------------------------------------------------------------*/


/* boolean */

Bool PutKeyBool (Hdr *hd, char *keyword, Bool value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
Bool value        i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

    extern int status;

    status = putKeyB (hd, keyword, value, comment);

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          PutKeyStr                                  */
/* ------------------------------------------------------------------*/


/* string */

int PutKeyStr (Hdr *hd, char *keyword, char *value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
char *value       i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

    extern int status;

    status = putKeyS (hd, keyword, value, comment);

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          KeyMissing                               */
/* ------------------------------------------------------------------*/


static void KeyMissing (char *keyword) {

    sprintf (MsgText, "Keyword = `%s'.", keyword);
    trlerror (MsgText);
}
