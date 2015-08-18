/* mainrej -- Main program for controlling WF3REJ */

# include <stdio.h>
# include <stdlib.h>        /* calloc */
# include <string.h>

# include "c_iraf.h"        /* for c_irafinit */
# include "hstio.h"

# include "wf3.h"
# include "wf3err.h"
# include "wf3rej.h"

int status = 0;             /* zero is OK */

int main (int argc, char **argv) {

    char    *input, output[SZ_LINE+1];      /* file names */
    clpar   par;                                    /* parameters used */
    int     newpar[MAX_PAR+1];          /* user specifiable parameters */
    char    mtype[SZ_FITS_VAL+1];      /* Role of exposure in association */

    int rej_command (int, char **, char **, char *, clpar *, int []);
    int Wf3Rej (char *, char *, char *, clpar *, int []);
    void WhichError (int);

    /* Initialize the IRAF interface */
    c_irafinit (argc, argv);

    /* Post HSTIO error handler */
    push_hstioerr (errchk);

    /* Initialize mtype to NULL to signal no change in ASN_MTYP for output*/
    mtype[0] = '\0';
    input = NULL;

    /* Initialize the structure for managing trailer file comments */
    InitTrlBuf ();

    
    /* Get input and output file names and switches in the command line. */
    if (rej_command (argc, argv, &input, output, &par, newpar)){
        exit (ERROR_RETURN);
    }


    /* Reject cosmic rays. */
    if (Wf3Rej (input, output, mtype, &par, newpar)) {
        sprintf (MsgText,"Error processing %s.", input);
        asnerror (MsgText);
    }

    CloseTrlBuf ();

    if (status) {
        WhichError (status);
        exit (ERROR_RETURN);
    } else{
        exit (status);
    }
}

