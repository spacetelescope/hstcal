#ifndef INCL_CTE_H
#define INCL_CTE_H

#define NUM_SCALE 4 /*number of scaling points, this is the 4 columns in the second table extension*/
#define TRAPS 999 /*max number of traps per column = rows in pctetab[1], valid traps are < 999999 in qlev*/
#define CTEFLAG 9999999 /*flag to ignore value in array during cte calculation*/


/* structure to hold CTE parameters from the reference files */
typedef struct {
    double scale512[RAZ_COLS]; /*scaling appropriate at row 512 */
    double scale1024[RAZ_COLS];/*scaling appropriate at row 1024 */
    double scale1536[RAZ_COLS];/*scaling appropriate at row 1536 */
    double scale2048[RAZ_COLS];/*scaling appropriate at row 2048 */
    double qlevq_data[TRAPS];/*charge packet size in electrons*/
    double dpdew_data[TRAPS];/*trap size in electrons*/
    double cte_date0; /*date of uvis install on hst in mjd*/
    double cte_date1; /*date of cte model pinning mjd*/
    double scale_frac; /*scaling of cte model relative to ctedate1*/
    double thresh; /*over subtraction threshold*/
    int cte_len; /*max length of cte trail */
    int n_forward; /* number of forward modeling iterations */
    int n_par; /*number of iterations in parallel transfer */
    int noise_mit; /*read noise mitigation algorithm*/
    int cte_traps; /*number of valid TRAPS in file for reallocation*/
    int wcol_data[TRAPS]; /*trap number, insync with number of traps*/
    int   iz_data[RAZ_COLS]; /*column number in raz format*/
    int fix_rocr; /*make allowance for readout cosmic rays*/
    char descrip2[CHAR_LINE_LENGTH+1]; /*descrip from table row, not read in for cte purposes*/
    char cte_name[CHAR_LINE_LENGTH+1]; /*name of cte algorithm */
    char cte_ver[CHAR_LINE_LENGTH+1]; /*version of algorithm */
    FloatHdrData *rprof; /*differential trail profile as image*/
    FloatHdrData *cprof; /*cummulative trail profile as image*/
} CTEParams;


/*USEFUL LIB FUNCTIONS*/
int GetGrp (WF3Info *, Hdr *);
void PrBegin (char *);
void PrEnd (char *);
void PrFileName (char *, char *);
void PrHdrInfo (char *, char *, char *);
void PrGrpBegin (char *, int);
void PrGrpEnd (char *, int);
int LoadHdr (char *, Hdr *);
void WF3Init (WF3Info *);
int MkName (char *, char *, char *, char *, char *, int);
int GetKeys (WF3Info *, Hdr *);
int GetSwitch (Hdr *, char *, int *);
int resistmean(float *, int, float, float *,float *,float *,float *);
void TimeStamp (char *, char *);
int  FileExists (char *);
void PrRefInfo (char *, char *,char *, char *, char *);
void PrSwitch (char *, int );
void WhichError (int);
int sub1d (SingleGroup *, int, SingleGroupLine *);
int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
int FindLine (SingleGroup *, SingleGroupLine *, int *, int *,int *, int *, int *);
int GetKeyInt (Hdr *, char *, int , int , int *);
int GetKeyDbl (Hdr *, char *, int , double , double *);
int GetKeyFlt (Hdr *, char *, int , float , float *);
int PutKeyInt (Hdr *, char *, int , char *);
int PutKeyFlt (Hdr *, char *, float , char *);
int PutKeyDbl (Hdr *, char *, double , char *);
int PutKeyStr(Hdr *, char *, char *, char *);
int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
int GetKeyStr (Hdr *, char *, int, char *, char *, int);
int streq_ic (char *, char *); /* case insensitive string equal */
int  MkOutName (char *, char **, char **, int, char *, int);
int MkNewExtn (char *, char *);
int Sub2Full(WF3Info *, SingleGroup *, SingleGroup *, int, int, int );
int Full2Sub(WF3Info *, SingleGroup *, SingleGroup *, int, int, int );
int CreateEmptyChip(WF3Info *, SingleGroup *);
int GetCorner (Hdr *, int, int *, int *);

/*FUNCTION SIGNATURES FOR CTE SPECIFIC CODE*/

int GetCTEPars (char *, CTEParams *);
void initCTEParams(CTEParams *);
int doCteBias (WF3Info *, SingleGroup *);
int GetCTEFlags (WF3Info *, Hdr *);
int a2d_raz(WF3Info *);
int raw2raz(WF3Info *, SingleGroup *, SingleGroup *, SingleGroup *);
int raz2rsz(WF3Info *, SingleGroup *, SingleGroup *, double , int );
int findPostScanBias(SingleGroup *, float *, float * );
int findPreScanBias(SingleGroup *, float *, float *);
int find_dadj(int ,int , double [][RAZ_ROWS], double [][RAZ_ROWS], double , double *);
int rsz2rsc(WF3Info *, SingleGroup *, SingleGroup *, CTEParams * );
int inverse_cte_blur(SingleGroup *, SingleGroup *, SingleGroup *, CTEParams *, int, double);
int sim_colreadout_l(double *, double *, double *, CTEParams *);
int CompareCTEParams(SingleGroup *, CTEParams *);
int cteHistory (WF3Info *, Hdr *);
int free_array(float **ptr, int rows, int columns);
int GetCTESwitch (WF3Info *, Hdr *);
int initCTETrl (char *, char *);

int makeRAZ(SingleGroup *, SingleGroup *, SingleGroup *);
int undoRAZ(SingleGroup *, SingleGroup *, SingleGroup *);

#endif /* INCL_CTE_H */
