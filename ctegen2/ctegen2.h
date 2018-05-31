#ifndef CTEGEN2_INCL
#define CTEGEN2_INCL

#include "hstio.h"

typedef struct {
    unsigned maxThreads;
    Bool verbose;
    int noise_mit; /*read noise mitigation algorithm*/
    int cte_len; // max length of cte trail
    int fix_rocr; /*make allowance for readout cosmic rays*/
    unsigned nTraps;
    unsigned nScaleTableColumns; //length of the arrays iz_data, scale512, scale1024, scale1536, scale2048 (nColumns both chips)
    unsigned nRows;
    unsigned nColumns;
    int n_forward; // number of forward modeling iterations
    int n_par; // number of iterations in parallel transfer
    unsigned cte_traps; // number of valid TRAPS in file for reallocation
    double thresh; /*over subtraction threshold*/
    double rn_amp; // read noise amplitude for clipping
    double cte_date0; /*date of instrument install on hst in mjd*/
    double cte_date1; /*date of cte model pinning mjd*/
    double scale_frac; /*scaling of cte model relative to ctedate1*/

    unsigned nRowsPerFullFrame; //includes overscan
    unsigned nColumnsPerFullFrame; //includes pre post overscan
    unsigned nColumnsPerChip; //nColumnsPerFullFrame / 2
    unsigned nRowsPerChip; //nRowsPerFullFrame
    unsigned nColumnsPerQuad; //nColumnsPerChip / 2
    unsigned nRowsPerQuad; //nRowsPerFullFrame

    //subarray offset positions within chip
    Bool isSubarray;
    unsigned chip;
    Bool refAndIamgeBinsIdenticle;
    unsigned rowOffset;
    unsigned columnOffset;

    //column alignment relative to old RAZ format, both chips side by side ordered CDAB
    unsigned razColumnOffset;

    //Actual image boundaries, i.e. excluding all overscan. Index is per quad
    //This is in reference to an aligned chip (amps bottom left corners, e.g. C & D or A & B)
    unsigned imageColumnsStart[2];
    unsigned imageColumnsEnd[2];
    unsigned imageRowsStart;
    unsigned imageRowsEnd;

    //prescan & post scan
    Bool hasPrescan[2]; //whether quad has physical serial  overscan
    Bool hasPostscan[2]; //whether quad has virtual serial overscan
    unsigned postscanWidth;
    unsigned prescanWidth;
    unsigned parallelOverscanWidth;


    //Bool indicating whether image exists in 1st and/or 2nd quad
    //0 index => A or C quad
    //1 indeax => B or D quad
    Bool quadExists[2];

    FloatHdrData * rprof; // differential trail profile as image
    FloatHdrData * cprof; // cumulative trail profile as image

    int * wcol_data; //trap number, insync with number of traps
    int    * iz_data; //column number in raz format
    double * scale512;  //scaling appropriate at row 512
    double * scale1024; //scaling appropriate at row 1024
    double * scale1536; //scaling appropriate at row 1536
    double * scale2048; //scaling appropriate at row 2048
    double * qlevq_data; // charge packet size in electrons
    double * dpdew_data; // trap size in electrons


    char descrip2[256]; /*descrip from table row, not read in for cte purposes*/
    char cte_name[256]; /*name of cte algorithm */
    char cte_ver[256]; /*version of algorithm */
} CTEParamsFast;

int inverseCTEBlurWithRowMajorIput(const SingleGroup * rsz, SingleGroup * rsc, const SingleGroup * trapPixelMap, CTEParamsFast * cte);
int inverseCTEBlur(const SingleGroup * rsz, SingleGroup * rsc, SingleGroup * trapPixelMap, CTEParamsFast * cte);
int forwardModel(const SingleGroup * input, SingleGroup * output, SingleGroup * trapPixelMap, CTEParamsFast * ctePars);

int simulatePixelReadout_v1_1(double * const pixelColumn, const float * const traps, const CTEParamsFast * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows);

int simulatePixelReadout_v1_2(double * const pixelColumn, const float * const traps, const CTEParamsFast * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows);

int simulateColumnReadout(double * const pixelColumn, const float * const traps, const CTEParamsFast * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows, const unsigned nPixelShifts);

Bool correctCROverSubtraction(float * const traps, const double * const pix_model, const double * const pix_observed,
        const unsigned nRows, const double threshHold);

int populateTrapPixelMap(SingleGroup * input, CTEParamsFast * params);
int cteSmoothImage(const SingleGroup * input, SingleGroup * output, CTEParamsFast * ctePars, double readNoiseAmp);
double find_dadjFast(const unsigned i ,const unsigned j, const unsigned nRows, const float * obsloc[3], const float * rszloc[3], const double readNoiseAmp);

//helpers
void initCTEParamsFast(CTEParamsFast * pars, const unsigned _nTraps, const unsigned _nRows, const unsigned _nColumns, const unsigned _nScaleTableColumns, const unsigned maxThreads);
int allocateCTEParamsFast(CTEParamsFast * pars);
void freeCTEParamsFast(CTEParamsFast * pars);

int populateImageFileWithCTEKeywordValues(SingleGroup *group, CTEParamsFast *pars);
int getCTEParsFromImageHeader(SingleGroup * input, CTEParamsFast * params);
int loadPCTETAB(char *filename, CTEParamsFast * params);
void ctewarn (char *message);
void cteerror (char *message);
void ctemessage (char *message);

//These shouldn't be here, they belong in their own header, in a central hstcal/inlclude with the code
//needing to be in hstcal/lib
int LoadHdr (char *, Hdr *);
int GetKeyStr (Hdr *, char *, int, char *, char *, int);
int GetKeyInt (Hdr *, char *, int , int , int *);
int GetKeyDbl (Hdr *, char *, int , double , double *);
int PutKeyInt (Hdr *, char *, int , char *);
int PutKeyDbl (Hdr *, char *, double , char *);
int PutKeyStr(Hdr *, char *, char *, char *);

#endif
