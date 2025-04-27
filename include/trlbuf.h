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
void CloseTrlBuf (struct TrlBuf * ptr);
void trlmessage (const char *fmt, ...);
void trlwarn (const char *fmt, ...);
void trlerror (const char *fmt, ...);
void trlopenerr (const char *filename);
void trlreaderr (const char *filename);
void trlkwerr (const char *keyword, const char *filename);
void trlfilerr (const char *filename);
void printfAndFlush (const char *message);
void trlGitInfo(void);

#endif
