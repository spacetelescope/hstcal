#ifndef TRLBUF_INCL
#define TRLBUF_INCL

#include <stdio.h>
#include "hstcal.h"

# define TRL_PREFIX "CALBEG"
# define TRL_EXTN   ".tra" // default extension for Trailer files

struct TrlBuf
{
    int overwrite;          // overwrite previous comments or append
    int quiet;              // Suppress STDOUT output?
    int usepref;            // Switch to specify whether preface is used
    int init;
    char *buffer;
    char *preface;          // comments common to all inputs
    char trlfile[CHAR_FNAME_LENGTH+1]; // name of output trailer file
    FILE *fp;               // pointer to open trailer file
};

int InitTrlFile (char *inlist, char *output);
int WriteTrlFile (void);
int InitTrlBuf (void);
void SetTrlPrefaceMode (int use);
void SetTrlOverwriteMode (int owrite);
void SetTrlQuietMode (int quiet);
void InitTrlPreface (void);
void ResetTrlPreface (void);
void CloseTrlBuf (void);
void trlmessage (char *message);
void trlwarn (char *message);
void trlerror (char *message);
void trlopenerr (char *filename);
void trlreaderr (char *filename);
void trlkwerr (char *keyword, char *filename);
void trlfilerr (char *filename);
void printfAndFlush (char *message);

#endif
