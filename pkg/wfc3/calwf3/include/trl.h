
# define TRL_EXTN   ".tra"      /* default extension for Trailer files */
# define FITS_EXTN  ".fits"     /* default extension */

/* Trailer file management routines */
int InitTrlBuf (void);
void SetTrlQuietMode (int quiet);
void SetTrlPrefaceMode (int use);
void CloseTrlBuf (void);
void InitTrlPreface (void);
void ResetTrlPreface (void);
int InitTrlFile (char *input, char *output);
int WriteTrlFile (void);

/* Trailer file comment output routines */
void trlmessage (char *message);
void trlwarn (char *message);
void trlerror (char *message);
void trlopenerr (char *filename);
void trlreaderr (char *name);
void trlkwerr (char *keyword, char *file);
void trlfilerr (char *name);
