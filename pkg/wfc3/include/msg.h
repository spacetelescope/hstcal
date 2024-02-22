#ifndef INCL_MSG_H
#define INCL_MSG_H

# include "hstcal.h"

# define SZ_CBUF           24  /* small buffer for e.g. rootname */
# define SZ_FITS_VAL       68
# define SZ_KEYWORD	    8

/* HSTIO error check function */
void errchk ();

/* The following function definitions handle the messages created
	during operations.  These will have counterparts which
	send output both the STDOUT and a TRL file. 
*/
void asnwarn (char *message);
void asnerror (char *message);
void ctemessage (char *message);
void cteerror (char *message);
void ctewarn (char *message);

#endif /* INCL_MSG_H */
