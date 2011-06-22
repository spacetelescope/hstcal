#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "xtables.h"
#include "hstio.h"

#include "acs.h"
#include "acserr.h"
#include "pcte.h"

/*
 * PixCteParams reads CTE characterization parameters from a FITS table file.
 * Data are read into a CTEParams struct, which is defined in pcte.h.
 * The reference file will change as Jay Anderson's characterization of CTE
 * develops and as new instruments are added.
 * - MRD 18 Feb. 2011
 */
int PixCteParams (const char *filename, CTEParams *pars) {
  
  extern int status; /* variable for return status */
  
  /* local variables */
  /* hstio variables */
  Hdr hdr_ptr;  
  
  /* xtables variables */
  IRAFPointer tbl_ptr;        /* xtables table pointer */
  IRAFPointer col_ptr_psin;   /* xtables column pointer for psi nodes */
  IRAFPointer col_ptr_logq;   /* xtables column pointer for log q */
  IRAFPointer col_ptr_dtde;   /* xtables column pointer for dtde */
  IRAFPointer col_ptr_qdtde;  /* xtables column pointer for q_dtde */
  IRAFPointer col_ptr_levs;   /* xtables column pointer for levels */
  
  /* names of data columns we want from the file */
  const char dtde[] = "DTDE";
  const char qdtde[] = "Q";
  const char node[] = "NODE";
  const char *logq_keys[NUM_LOGQ] = {"LOG_Q_1","LOG_Q_2","LOG_Q_3","LOG_Q_4"};
  const char level[] = "LEVEL";
  
  /* variable for filename + extension number. 
   * not sure it's cool to use strlen here but it works for now */
  char filename_wext[strlen(filename) + 4];
  
  /* iteration variable */
  int j, k;
  
  /* functions from calacs/lib */
  int LoadHdr (char *input, Hdr *phdr);
  int GetKeyInt (Hdr *hd, char *keyword, int use_def, int def, int *value);
  int GetKeyDbl (Hdr *hd, char *keyword, int use_def, double def, double *value);

  /* load primary header */
  if (LoadHdr(filename, &hdr_ptr)) {
    sprintf(MsgText,"(pctecorr) Error loading header from %s",filename);
    trlerror(MsgText);
    status = OPEN_FAILED;
    return status;
  }
  
  /* read RN2_NIT keyword from primary header */
  if (GetKeyDbl(&hdr_ptr, "RN_CLIP", NO_DEFAULT, -999, &pars->rn_clip)) {
    trlerror("(pctecorr) Error reading RN_CLIP keyword from PCTETAB");
    status = KEYWORD_MISSING;
    return status;
  }
  
  /* read SIM_NIT keyword from primary header */
  if (GetKeyInt(&hdr_ptr, "SIM_NIT", NO_DEFAULT, -999, &pars->sim_nit)) {
    trlerror("(pctecorr) Error reading SIM_NIT keyword from PCTETAB");
    status = KEYWORD_MISSING;
    return status;
  }
  
  /* read SHFT_NIT keyword from primary header */
  if (GetKeyInt(&hdr_ptr, "SHFT_NIT", NO_DEFAULT, -999, &pars->shft_nit)) {
    trlerror("(pctecorr) Error reading SHFT_NIT keyword from PCTETAB");
    status = KEYWORD_MISSING;
    return status;
  }
  
  /* done reading stuff from the primary header */
  freeHdr(&hdr_ptr);
  
  /****************************************************************************/
  /* read DTDE/Q data from first table extensions */
  /* make filename + ext number 1 */
  sprintf(filename_wext, "%s[%i]", filename, 1);
  
  /* open CTE parameters file to extension number 1 */
  tbl_ptr = c_tbtopn(filename_wext, IRAF_READ_ONLY, 0);
  if (c_iraferr()) {
    sprintf(MsgText,"(pctecorr) Error opening %s with xtables",filename_wext);
    trlerror(MsgText);
    status = OPEN_FAILED;
    return status;
  }
  
  /* read data from table */
  /* get column pointer for dtde */
  c_tbcfnd1(tbl_ptr, dtde, &col_ptr_dtde);
  if (c_iraferr() || col_ptr_dtde == 0) {
    sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",dtde);
    trlerror(MsgText);
    status = COLUMN_NOT_FOUND;
    return status;
  }
  
  /* get column pointer for q_dtde */
  c_tbcfnd1(tbl_ptr, qdtde, &col_ptr_qdtde);
  if (c_iraferr() || col_ptr_qdtde == 0) {
    sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",qdtde);
    trlerror(MsgText);
    status = COLUMN_NOT_FOUND;
    return status;
  }
  
  /* loop over table rows */
  for (j = 0; j < NUM_PHI; j++) {
    /* get dtde from this row */
    c_tbegtd(tbl_ptr, col_ptr_dtde, j+1, &pars->dtde_l[j]);
    if (c_iraferr()) {
      sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, dtde);
      trlerror(MsgText);
      status = TABLE_ERROR;
      return status;
    }
    
    /* get q_dtde from this row */
    c_tbegti(tbl_ptr, col_ptr_qdtde, j+1, &pars->q_dtde[j]);
    if (c_iraferr()) {
      sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, qdtde);
      trlerror(MsgText);
      status = TABLE_ERROR;
      return status;
    }
  }
  
  /* close CTE parameters file for extension 1
   * end of reading DTDE/Q data */
  c_tbtclo(tbl_ptr);
  /****************************************************************************/
  
  /****************************************************************************/
  /* read NODE/LOG_Q_# data from second table extension */
  /* make filename + ext number 1 */
  sprintf(filename_wext, "%s[%i]", filename, 2);
  
  /* open CTE parameters file to extension number 2 */
  tbl_ptr = c_tbtopn(filename_wext, IRAF_READ_ONLY, 0);
  if (c_iraferr()) {
    sprintf(MsgText,"(pctecorr) Error opening %s with xtables",filename_wext);
    trlerror(MsgText);
    status = OPEN_FAILED;
    return status;
  }
  
  /* read data from table */
  /* get column pointer for the psi nodes */
  c_tbcfnd1(tbl_ptr, node, &col_ptr_psin);
  if (c_iraferr() || col_ptr_psin == 0) {
    sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",node);
    trlerror(MsgText);
    status = COLUMN_NOT_FOUND;
    return status;
  }
  
  /* iterate over table rows */
  for (j = 0; j < NUM_PSI; j++) {
    /* get the psi node number */
    c_tbegti(tbl_ptr, col_ptr_psin, j+1, &pars->psi_node[j]);
    if (c_iraferr()) {
      sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, node);
      trlerror(MsgText);
      status = TABLE_ERROR;
      return status;
    }
    
    /* loop over table columns to read log q values */
    for (k = 0; k < NUM_LOGQ; k++) {
      /* get column pointer for this log q */
      c_tbcfnd1(tbl_ptr, logq_keys[k], &col_ptr_logq);
      if (c_iraferr() || col_ptr_logq == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",logq_keys[k]);
        trlerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
      }
      
      /* read this log q value */
      c_tbegtd(tbl_ptr, col_ptr_logq, j+1, &pars->chg_leak[j*NUM_LOGQ + k]);
      if (c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",
                j+1, logq_keys[k]);
        trlerror(MsgText);
        status = TABLE_ERROR;
        return status;
      }
    }
  }
  
  /* close CTE parameters file for extension 2
   * end of reading NODE/LOG_Q_# data */
  c_tbtclo(tbl_ptr);
  /****************************************************************************/
  
  /****************************************************************************/
  /* read LEVEL data from third table extensions */
  /* make filename + ext number 3 */
  sprintf(filename_wext, "%s[%i]", filename, 3);
  
  /* open CTE parameters file to extension number 3 */
  tbl_ptr = c_tbtopn(filename_wext, IRAF_READ_ONLY, 0);
  if (c_iraferr()) {
    sprintf(MsgText,"(pctecorr) Error opening %s with xtables",filename_wext);
    trlerror(MsgText);
    status = OPEN_FAILED;
    return status;
  }
  
  /* read data from table */
  /* get column pointer for level */
  c_tbcfnd1(tbl_ptr, level, &col_ptr_levs);
  if (c_iraferr() || col_ptr_levs == 0) {
    sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",level);
    trlerror(MsgText);
    status = COLUMN_NOT_FOUND;
    return status;
  }
  
  /* loop over table rows */
  for (j = 0; j < NUM_LEV; j++) {
    /* get level from this row */
    c_tbegti(tbl_ptr, col_ptr_levs, j+1, &pars->levels[j]);
    if (c_iraferr()) {
      sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, level);
      trlerror(MsgText);
      status = TABLE_ERROR;
      return status;
    }
  }
  
  /* close CTE parameters file for extension 3
   * end of reading LEVEL data */
  c_tbtclo(tbl_ptr);
  /****************************************************************************/
  
  return status;
}

/*
 * Some CTE parameters, especially different numbers of iterations, will likely
 * end up as header keywords in images so that users can tune them. This routine
 * will first check whether these parameters are set in the image header. If they
 * are then the values found there will be used instead of the values read from
 * the PCTETAB. If the values are not set they will be populated with the values
 * from the PCTETAB.
 *
 * For now this doesn't do much because it hasn't been decided which keywords
 * will be user tunable in this manner, this is just the recipe.
 * MRD 17 Mar. 2011
 */
int CompareCteParams(SingleGroup *x, CTEParams *pars) {
  
  extern int status;
  
  /* functions from calacs/lib/key.c */
  int GetKeyInt (Hdr *hd, char *keyword, int use_def, int def, int *value);
  int GetKeyDbl (Hdr *hd, char *keyword, int use_def, double def, double *value);
  int PutKeyInt (Hdr *hd, char *keyword, int value, char *comment);
  int PutKeyDbl (Hdr *hd, char *keyword, double value, char *comment);
  
  int def_value = -99;
  double def_value_dbl = -99.0;
  
  double rn_clip;
  int sim_nit;
  int shft_nit;
  
  /****************************************************************************/
  /* get value and perform check for read noise iterations */
  if (GetKeyDbl(x->globalhdr, "PCTERNCL", USE_DEFAULT, def_value_dbl, &rn_clip)) {
    trlerror("(pctecorr) Error reading PCTERNCL keyword from image header");
    status = KEYWORD_MISSING;
    return status;
  }
  
  /* if the keyword value is not there or set to zero, update it with the value
   * read from the PCTETAB, otherwise use the read value in the CTE correction. */
  if (rn_clip == def_value || rn_clip == 0) {
    if (PutKeyDbl(x->globalhdr, "PCTERNCL", pars->rn_clip, "PCTE readnoise amplitude")) {
      trlerror("(pctecorr) Error updating PCTERNCL keyword in image header");
      return (status = HEADER_PROBLEM);
    }
  } else {
    pars->rn_clip = rn_clip;
  }
  /****************************************************************************/
  
  /****************************************************************************/
  /* get value and perform check for readout simulation iterations */
  if (GetKeyInt(x->globalhdr, "PCTESMIT", USE_DEFAULT, def_value, &sim_nit)) {
    trlerror("(pctecorr) Error reading PCTESMIT keyword from image header");
    status = KEYWORD_MISSING;
    return status;
  }
  
  /* if the keyword value is not there or set to zero, update it with the value
   * read from the PCTETAB, otherwise use the read value in the CTE correction. */
  if (sim_nit == def_value || sim_nit == 0) {
    if (PutKeyInt(x->globalhdr, "PCTESMIT", pars->sim_nit, "PCTE readout simulation iterations")) {
      trlerror("(pctecorr) Error updating PCTESMIT keyword in image header");
      return (status = HEADER_PROBLEM);
    }
  } else {
    pars->sim_nit = sim_nit;
  }
  /****************************************************************************/
  
  /****************************************************************************/
  /* get value and perform check for readout shift iterations */
  if (GetKeyInt(x->globalhdr, "PCTESHFT", USE_DEFAULT, def_value, &shft_nit)) {
    trlerror("(pctecorr) Error reading PCTESHFT keyword from image header");
    status = KEYWORD_MISSING;
    return status;
  }
  
  /* if the keyword value is not there or set to zero, update it with the value
   * read from the PCTETAB, otherwise use the read value in the CTE correction. */
  if (shft_nit == def_value || shft_nit == 0) {
    if (PutKeyInt(x->globalhdr, "PCTESHFT", pars->sim_nit, "PCTE readout number of shifts")) {
      trlerror("(pctecorr) Error updating PCTESHFT keyword in image header");
      return (status = HEADER_PROBLEM);
    }
  } else {
    pars->shft_nit = shft_nit;
  }
  /****************************************************************************/
  
  return status;
}

/*
 * CalcCteFrac calculates the multiplicative factor that accounts for the
 * worsening of CTE over time. This is currently a linear function valid over
 * the whole life of ACS/WFC, but this will change soon. Jay Anderson has
 * discovered that the slope of the CTE scaling is not constant so the function
 * for CTE frac will be different depending on the obs. start time. The plan
 * is to add the CTE frac parameterization to the CTE params reference file
 * once Jay has them nailed down.
 * - MRD 18 Feb. 2011
 *
 * Constants for instrument names are defined in PixCteCorr.h.
 */
double CalcCteFrac(const double mjd, const int instrument) {
  
  /* variables used to calculate the CTE scaling slope */
  double mjd_pt1, mjd_pt2; /* the MJD points at which the scaling is defined */
  double cte_pt1, cte_pt2; /* the CTE frac points at which the scaling is defined */
  
  double cte_frac;         /* return value */
  
  if (instrument == ACSWFC) {
    cte_pt1 = 0.0;
    cte_pt2 = 1.0;
    mjd_pt1 = 52335.0;  /* March 2, 2002 */
    mjd_pt2 = 55263.0;  /* March 8, 2010 */
  } else {
    printf("Instrument not found: %i\n",instrument);
    return -9999.0;
  }
  
  cte_frac = ((cte_pt2 - cte_pt1) / (mjd_pt2 - mjd_pt1)) * (mjd - mjd_pt1);
  
  return cte_frac;
}

/*
 * InterpolatePsi fills in the sparser array describing trail profiles read
 * from the CTE parameters file. Is there any reason the whole profile can't
 * be written to the parameters file? It's only 100 elements long.
 * - MRD 18 Feb. 2011
 *
 * Inputs chg_leak and psi_node are arrays read from the CTE parameters file.
 * Output chg_leak_interp has all the data chg_leak plus interpolated data where
 * chg_leak has none.
 */
int InterpolatePsi(const double chg_leak[NUM_PSI*NUM_LOGQ], const int psi_node[NUM_PSI],
                   double chg_leak_interp[MAX_TAIL_LEN*NUM_LOGQ],
                   double chg_open_interp[MAX_TAIL_LEN*NUM_LOGQ]) {
  
  /* status variable for return */
  extern int status;
  
  /* index variables for tracking where we are in psi_node/chg_leak.
   * these will always be one apart. so we don't really need two, but
   * I like it for cleanliness. */
  int pn_i1 = 0;
  int pn_i2 = 1;
  
  double interp_frac; /* the fraction of the distance between psi_node1 and
                       * psi_node2 are we interpolating at */
  
  double sum_rel;     /* total probability of release */
  double sum_cum;     /* running total probability of release */
  
  /* iteration variables */
  int i, j;
  
  /* loop over all pixels in the trail and calculate the profile at each q
   * if it isn't already in chg_leak */
  for (i = 0; i < MAX_TAIL_LEN; i++) { 
    /* do we match an existing chg_leak row? */
    if (i+1 == psi_node[pn_i1]) {
      /* if so, copy it over */
      for (j = 0; j < NUM_LOGQ; j++) {
        chg_leak_interp[i*NUM_LOGQ + j] = chg_leak[pn_i1*NUM_LOGQ + j];
      }      
    } else {
      /* no match, need to interpolate */
      interp_frac = ((double) (i+1 - psi_node[pn_i1])) / 
      ((double) (psi_node[pn_i2] - psi_node[pn_i1]));
      /* loop over each q column */
      for (j = 0; j < NUM_LOGQ; j++) {
        chg_leak_interp[i*NUM_LOGQ + j] = chg_leak[pn_i1*NUM_LOGQ + j] + 
        (interp_frac * 
         (chg_leak[pn_i2*NUM_LOGQ + j] - chg_leak[pn_i1*NUM_LOGQ + j]));
      }
    }
    /* if the next of row psi_node is an existing row in chg_leak we should
     * increment our indices so we start interpolating between the next pair */
    if ((i+2 == psi_node[pn_i2]) && (i < MAX_TAIL_LEN)) {
      pn_i1++;
      pn_i2++;
    }
  }
  
  /* perform tail normalization and cumulative release probability calculation */
  for (i = 0; i < NUM_LOGQ; i++) {
    sum_rel = 0.0;
    
    /* get total in this Q column */
    for (j = 0; j < MAX_TAIL_LEN; j++) {
      sum_rel += chg_leak_interp[j*NUM_LOGQ + i];
    }
    
    /* normalize chg_leak_interp by total */
    for (j = 0; j < MAX_TAIL_LEN; j++) {
      chg_leak_interp[j*NUM_LOGQ + i] = chg_leak_interp[j*NUM_LOGQ + i]/sum_rel;
    }
    
    /* calculate cumulative probability of release */
    sum_cum = 0.0;
    
    for (j = 0; j < MAX_TAIL_LEN; j++) {
      sum_cum += chg_leak_interp[j*NUM_LOGQ + i];
      chg_open_interp[j*NUM_LOGQ + i] = 1.0 - sum_cum;
    }
  }
  
  return status;
}

/*
 * InterpolatePhi interpolates information read from the CTE parameters file
 * that describes the amount of charge in the CTE tail and where the charge is.
 * -MRD 21 Feb. 2011
 *
 * Input dtde_l is read from the CTE parameters file and cte_frac is calculated
 * from the observation start date by CalcCteFrac.
 * Outputs dtde_q, q_pix_array, and pix_q_array are arrays MAX_PHI long (should
 * be 99999) and ycte_qmax is an integer.
 */
int InterpolatePhi(const double dtde_l[NUM_PHI], const int q_dtde[NUM_PHI],
                   const int shft_nit, double dtde_q[MAX_PHI]) {
  
  /* status variable for return */
  extern int status;
  
  int p; /* iteration variable over phi nodes in reference file */
  int q; /* iteration variable over single phi values between nodes in ref file */
  
  /* interpolation calculation variables */
  double interp_pt;   /* point at which we're interpolating data */
  double interp_dist; /* difference between interp_pt and low_node */
  double interp_val;  /* interpolated value */
  
  /* upper and lower bounds of interpolation range */
  double log_qa, log_qb;
  double log_da, log_db;
  
  /* something for holding intermediate calculation results */
  double qtmp;
  
  for (p = 0; p < NUM_PHI-1; p++) {
    log_qa = log10((double) q_dtde[p]);
    log_qb = log10((double) q_dtde[p+1]);
    log_da = log10(dtde_l[p]);
    log_db = log10(dtde_l[p+1]);
    
    for (q = q_dtde[p]; q < q_dtde[p+1]; q++) {
      interp_pt = log10((double) q);
      interp_dist = (interp_pt - log_qa) / (log_qb - log_qa);
      interp_val = log_da + (interp_dist * (log_db - log_da));
      
      qtmp = pow(10, interp_val)/(double) CTE_REF_ROW;
      qtmp = pow((1.0 - qtmp), (double) CTE_REF_ROW/ (double) shft_nit);
      dtde_q[q-1] = 1.0 - qtmp;
    }
  }
  
  qtmp = pow((1.0 - (dtde_l[NUM_PHI-1]/CTE_REF_ROW)),CTE_REF_ROW/shft_nit);
  dtde_q[MAX_PHI-1] = 1.0 - qtmp;
  
  return status;
}

/* In this function we're interpolating the tail arrays over the Q dimension
 * and reducing the arrays to contain data at only the charge levels
 * specified in the levels array. */
int FillLevelArrays(const double chg_leak_kt[MAX_TAIL_LEN*NUM_LOGQ],
                    const double chg_open_kt[MAX_TAIL_LEN*NUM_LOGQ],
                    const double dtde_q[MAX_PHI], const int levels[NUM_LEV],
                    double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
                    double chg_open_lt[MAX_TAIL_LEN*NUM_LEV],
                    double dpde_l[NUM_LEV],
                    int tail_len[NUM_LEV]) {
  
  /* status variable for return */
  extern int status;
  
  int l,t;  /* iteration variables for tail and levels */
  int q;    /* iteration variable for q levels in between those specified in levels */
  
  /* container for cumulative dtde_q */
  double cpde_l[NUM_LEV];
  
  /* variable for running sum of dtde_q */
  double sum = 0.0;
  
  int logq_ind; /* index of lower logq used for interpolation */
  
  double logq;           /* log of charge level */
  double logq_min = 1.0; /* min value for logq */
  double logq_max = 3.999; /* max value for logq */
  
  double interp_dist; /* difference between logp and lower logq */
  
  dpde_l[0] = 0.0;
  cpde_l[0] = 0.0;
  
  for (t = 0; t < MAX_TAIL_LEN; t++) {
    chg_leak_lt[t*NUM_LEV] = chg_leak_kt[t*NUM_LOGQ];
    chg_open_lt[t*NUM_LEV] = chg_open_kt[t*NUM_LOGQ];
  }
  
  for (l = 1; l < NUM_LEV; l++) {
    for (q = levels[l-1]; q < levels[l]; q++) {
      sum += dtde_q[q];
    }
    
    cpde_l[l] = sum;
    dpde_l[l] = cpde_l[l] - cpde_l[l-1];
    
    /* calculate logq with min/max clipping */
    logq = log10((double) q);
    if (logq < logq_min) {
      logq = logq_min;
    } else if (logq > logq_max) {
      logq = logq_max;
    }
    
    /* set logq_ind for this logq */
    if (logq < 2) {
      logq_ind = 0;
    } else if (logq < 3) {
      logq_ind = 1;
    } else {
      logq_ind = 2;
    }
    
    interp_dist = logq - floor(logq);
    
    for (t = 0; t < MAX_TAIL_LEN; t++) {
      chg_leak_lt[t*NUM_LEV + l] = ((1.0 - interp_dist) * chg_leak_kt[t*NUM_LOGQ + logq_ind]) +
      (interp_dist * chg_leak_kt[t*NUM_LOGQ + logq_ind+1]);
      chg_open_lt[t*NUM_LEV + l] = ((1.0 - interp_dist) * chg_open_kt[t*NUM_LOGQ + logq_ind]) +
      (interp_dist * chg_open_kt[t*NUM_LOGQ + logq_ind+1]);
    }
  }
  
  /* calculate max tail lengths for each level */
  for (l = 0; l < NUM_LEV; l++) {
    tail_len[l] = MAX_TAIL_LEN;
    
    for (t = MAX_TAIL_LEN-1; t >= 0; t--) {
      if (chg_leak_lt[t*NUM_LEV + l] == 0) {
        tail_len[l] = t+1;
      } else {
        break;
      }
    }
  }
  
  return status;
}

/*
 * Attempt to separate readout noise from signal, since CTI happens before
 * readout noise is added to the signal.
 *
 * The clipping parameter pclip controls the maximum amount by which a pixel
 * will be modified, or the maximum amplitude of the read noise.
 */
int DecomposeRN(const int arrx, const int arry, const double data[arrx*arry],
                const double pclip, double sig_arr[arrx*arry], double noise_arr[arrx*arry]) {
  
  /* status variable for return */
  extern int status;
  
  /* iteration variables */
  int i, j;
  
  /* array to hold pixel means from 20 surrounding pixels */
  double * means;
  
  /* array to hold clipped difference between data and means */
  double * diffs;
  double diff;
  
  /* array to hold smoothed diffs */
  double * sm_diffs;
  
  /* index variables */
  int iind, jind;
  
  /* get space for arrays */
  means = (double *) malloc(arrx * arry * sizeof(double));
  diffs = (double *) malloc(arrx * arry * sizeof(double));
  sm_diffs = (double *) malloc(arrx * arry * sizeof(double));
  
  /* calculate means array as average of 20 surrounding pixel, not including
   * central pixel. */
  for (i = 0; i < arrx; i++) {
    for (j = 0; j < arry; j++) {
      /* if this pixel is within 2 rows/columns of the edge get the median
       * from the pixel in the third row/column from the edge. */
      iind = i;
      if (iind <= 1) {
        iind = 2;
      } else if (iind >= arrx-2) {
        iind = arrx - 3;
      }
      jind = j;
      if (jind <= 1) {
        jind = 2;
      } else if (jind >= arry-2) {
        jind = arry - 3;
      }
      
      means[i*arry + j] = (data[(iind+0)*arry + jind+1] + 
                           data[(iind+0)*arry + jind-1] + 
                           data[(iind+1)*arry + jind+0] + 
                           data[(iind-1)*arry + jind-0] + 
                           data[(iind-1)*arry + jind+1] + 
                           data[(iind-1)*arry + jind+1] + 
                           data[(iind-1)*arry + jind-1] + 
                           data[(iind-1)*arry + jind-1] + 
                           data[(iind+0)*arry + jind+2] + 
                           data[(iind-0)*arry + jind-2] + 
                           data[(iind+2)*arry + jind+0] + 
                           data[(iind-2)*arry + jind-0] + 
                           data[(iind+1)*arry + jind+2] + 
                           data[(iind-1)*arry + jind+2] + 
                           data[(iind+1)*arry + jind-2] + 
                           data[(iind-1)*arry + jind-2] + 
                           data[(iind+2)*arry + jind+1] + 
                           data[(iind+2)*arry + jind+1] + 
                           data[(iind-2)*arry + jind-1] + 
                           data[(iind-2)*arry + jind-1])/20.0;
    }
  }
  
  /* calculate clipped differences array */
  for (i = 0; i < arrx; i++) {
    for (j = 0; j < arry; j++) {
      diff = data[i*arry + j] - means[i*arry + j];
      
      if (diff < -pclip) {
        diffs[i*arry + j] = -pclip;
      } else if (diff > pclip) {
        diffs[i*arry + j] = pclip;
      } else {
        diffs[i*arry + j] = diff;
      }
    }
  }
  
  /* to avoid systematic reduction of sources we insist that the average
   * clipping within any 5-pixel vertical window is zero. */
  for (i = 0; i < arrx; i++) {
    for (j = 0; j < arry; j++) {
      iind = i;
      if (iind <= 2) {
        iind = 2;
      } else if (iind >= arrx-2) {
        iind = arrx - 3;
      }
      
      sm_diffs[i*arry + j] = (diffs[(iind-2)*arry + j] + 
                              diffs[(iind-1)*arry + j] + 
                              diffs[(iind-0)*arry + j] + 
                              diffs[(iind+1)*arry + j] + 
                              diffs[(iind+2)*arry + j])/5.0;
      
      diff = diffs[i*arry + j] - sm_diffs[i*arry + j];
      
      if (diff < -pclip) {
        noise_arr[i*arry + j] = -pclip;
      } else if (diff > pclip) {
        noise_arr[i*arry + j] = pclip;
      } else {
        noise_arr[i*arry + j] = diff;
      }
      
      sig_arr[i*arry + j] = data[i*arry + j] - noise_arr[i*arry + j];
    }
  }
  
  free(means);
  free(diffs);
  free(sm_diffs);
  
  return status;
}
