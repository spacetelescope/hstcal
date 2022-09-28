/* mainrej -- Main program for controlling WF3REJ */

# include <stdio.h>
# include <stdlib.h>        /* calloc */
# include <string.h>

#include "hstcal_memory.h"
#include "hstcal.h"
# include "c_iraf.h"        /* for c_irafinit */
# include "hstio.h"

# include "wf3.h"
# include "hstcalerr.h"
# include "wf3rej.h"
# include "hstcalversion.h"
# include "trlbuf.h"

extern int status;
struct TrlBuf trlbuf = { 0 };

/* Standard string buffer for use in messages */
char MsgText[MSG_BUFF_LENGTH]; // Global char auto initialized to '\0'

int main (int argc, char **argv) {

    char    *input, output[CHAR_LINE_LENGTH+1];      /* file names */
    clpar   par;                                    /* parameters used */
    int     newpar[MAX_PAR+1];          /* user specifiable parameters */
    char    mtype[SZ_FITS_VAL+1];      /* Role of exposure in association */

    int rej_command (int, char **, char **, char *, clpar *, int []);
    int Wf3Rej (char *, char *, char *, clpar *, int []);
    void WhichError (int);

    status = 0;

    /* Initialize the IRAF interface */
    c_irafinit (argc, argv);

    /* Post HSTIO error handler */
    push_hstioerr (errchk);

    /* Initialize mtype to NULL to signal no change in ASN_MTYP for output*/
    mtype[0] = '\0';
    input = NULL;

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
    /* Initialize the structure for managing trailer file comments */
    InitTrlBuf ();
    addPtr(&ptrReg, &trlbuf, &CloseTrlBuf);
    
    /* Get input and output file names and switches in the command line. */
    if (rej_command (argc, argv, &input, output, &par, newpar)){
        if (input)
            free(input);
        freeOnExit(&ptrReg);
        exit (ERROR_RETURN);
    }
    addPtr(&ptrReg, input, &free);

    /* Reject cosmic rays. */
    if (Wf3Rej (input, output, mtype, &par, newpar)) {
        sprintf (MsgText,"Error processing %s.", input);
        asnerror (MsgText);
    }


    if (status) {
        WhichError (status);
        freeOnExit(&ptrReg);
        exit (ERROR_RETURN);
    } else{
        freeOnExit(&ptrReg);
        exit (status);
    }
}

