#define NUM_SCALE 4 /*number of scaling points, this is the 4 columns in the second table extension*/
#define COLUMNS  8412    /*columns for each row */
#define ROWS 2070 /*2070 rows*/
#define TRAPS 999 /*max number of traps per column = rows in pctetab[1], valid traps are < 999999 in qlev*/   
   
/* structure to hold CTE parameters from the reference files */
typedef struct {
    char cte_name; /*name of cte algorithm */
    char cte_ver; /*version of algorithm */
    double cte_date0; /*date of uvis install on hst in mjd*/
    double cte_date1; /*date of cte model pinning mjd*/
    int cte_len; /*max length of cte trail */
    float   rn_amp; /*read noise amplitude for clipping */
    int n_forward; /* number of forward modeling iterations */
    int n_par; /*numver of iterations in parallel transfer */
    float scale_frac; /*scaling of cte model relative to ctedate1*/
    int noise_mit; /*read noise mitigation algorithm*/
    float thresh; /*over subtraction threshold*/        
    int cte_traps; /*number of valid TRAPS in file for reallocation*/
    int wcol_data[TRAPS]; /*trap number, insync with number of traps*/
    int qlevq_data[TRAPS];/*charge packet size in electrons*/
    float dpdew_data[TRAPS];/*trap size in electrons*/  
    int   iz_data[COLUMNS]; /*column number in raz format*/
    double scale512[COLUMNS]; /*scaling appropriate at row 512 */
    double scale1024[COLUMNS];/*scaling appropriate at row 1024 */
    double scale1536[COLUMNS];/*scaling appropriate at row 1536 */
    double scale2048[COLUMNS];/*scaling appropriate at row 2048 */
    char descrip2; /*descrip from table row, not read in for cte purposes*/
    int fix_rocr; /*make allowance for readout cosmic rays*/
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

int sub1d (SingleGroup *, int, SingleGroupLine *);
int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
int FindLine (SingleGroup *, SingleGroupLine *, int *, int *,int *, int *, int *);
int GetKeyInt (Hdr *hd, char *keyword, int use_def, int def, int *value);
int GetKeyDbl (Hdr *hd, char *keyword, int use_def, double def, double *value);
int GetKeyFlt (Hdr *hd, char *keyword, int use_def, float def, float *value);
int PutKeyInt (Hdr *hd, char *keyword, int value, char *comment);
int PutKeyFlt (Hdr *hd, char *keyword, float value, char *comment);
int PutKeyDbl (Hdr *hd, char *keyword, double value, char *comment);
int PutKeyStr(Hdr *hd, char *keyword, char *value, char *comment);
int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
int GetKeyStr (Hdr *, char *, int, char *, char *, int);
int streq_ic (char *, char *); /* case insensitive string equal */

/*FUNCTION SIGNATURES FOR CTE SPECIFIC CODE*/

int GetCTEPars (char *, CTEParams *);
void initCTEParams(CTEParams *);
int doCteBias (WF3Info *, SingleGroup *);
int GetCTEFlags (WF3Info *, Hdr *);
int a2d_raz(WF3Info *);
int raw2raz(WF3Info *, SingleGroup *, SingleGroup *, SingleGroup *);
int raz2rsz(WF3Info *, SingleGroup *, SingleGroup *, float , int );
int findPostScanBias(SingleGroup *, float *, float *);
int findPreScanBias(SingleGroup *, float *, float *);
int find_dadj(int ,int , float [][ROWS], float [][ROWS], float , double *);
int rsz2rsc(WF3Info *, SingleGroup *, SingleGroup *, CTEParams * );
int inverse_cte_blur(SingleGroup *, SingleGroup *, SingleGroup *, CTEParams *, int, double);
int sim_colreadout_l(float *, float *, float *, CTEParams *);
int CompareCTEParams(SingleGroup *, CTEParams *);
int cteHistory (WF3Info *, Hdr *);
float **alloc_array(int rows, int columns);
int free_array(float **ptr, int rows, int columns);
int GetCTESwitch (WF3Info *, Hdr *);
int initCTETrl (char *, char *);

int makesciRAZ(SingleGroup *, SingleGroup *, SingleGroup *);
int undosciRAZ(SingleGroup *, SingleGroup *, SingleGroup *);
