# define SZ_CBUF           24  /* small buffer for e.g. rootname */
# define SZ_FNAME          255
# define SZ_LINE           255
# define SZ_FITS_VAL       68
# define SZ_KEYWORD	    8

/* Standard string buffer for use in messages */
char MsgText[SZ_LINE+1];

/* HSTIO error check function */
void errchk ();

/* The following function definitions handle the messages created
	during operations.  These will have counterparts which
	send output both the STDOUT and a TRL file. 
*/
void asnwarn (char *message);
void asnerror (char *message);
void asnmessage (char *message);
void ctemessage (char *message);
void cteerror (char *message);
void ctewarn (char *message);

# define WARN_PREFIX    "Warning    "
# define ERR_PREFIX     "ERROR:    "

