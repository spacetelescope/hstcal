#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "xtables.h"
#include "hstio.h"

#include "acs.h"
#include "hstcalerr.h"
#include "pcte.h"

/*
 * PixCteParams reads CTE characterization parameters from a FITS table file.
 * Data are read into a CTEParams struct, which is defined in pcte.h.
 * The reference file will change as Jay Anderson's characterization of CTE
 * develops and as new instruments are added.
 * - MRD 18 Feb. 2011
 */

int getCTE_NAME(char * filename, char * cteName, int cteNameBufferLength)
{
    extern int status;

    Hdr hdr;
    initHdr(&hdr);//For correctness (this is foolishly initialized in LoadHdr())

    if (LoadHdr(filename, &hdr))
    {
        sprintf(MsgText,"(pctecorr) Error loading header from %s",filename);
        trlerror(MsgText);
        return (status = OPEN_FAILED);
    }

    /* GET CTE_NAME KEYWORD */
    if (GetKeyStr(&hdr, "CTE_NAME", USE_DEFAULT, "", cteName, cteNameBufferLength))
    {
        cteerror("(pctecorr) Error reading CTE_NAME keyword from PCTETAB");
        return KEYWORD_MISSING;//NOTE: DO NOT update status, this is the callers responsibility
    }

    return (status = ACS_OK);
}

int PixCteParams (char *filename, const double expstart, ACSCTEParams *pars) {

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
  /*IRAFPointer col_ptr_mjd;*/    /* xtables column pointer for mjd */
  /*IRAFPointer col_ptr_scale;*/  /* xtables column pointer for scale */
  IRAFPointer col_ptr_amp;    /* xtables column pointer for amp */

  /* names of data columns we want from the file */
  const char dtde[] = "DTDE";
  const char qdtde[] = "Q";
  const char node[] = "NODE";
  const char *logq_keys[NUM_LOGQ] = {"LOG_Q_1","LOG_Q_2","LOG_Q_3","LOG_Q_4"};
  const char level[] = "LEVEL";
  const char mjdstr[] = "MJD";
  const char scalestr[] = "SCALE";
  const char *amp_keys[NAMPS] = {"AMPA","AMPB","AMPC","AMPD"};

  /* variable for filename + extension number. */
  char filename_wext[strlen(filename) + 4];

  /* iteration variable */
  int j, k, l;

  /* number of CHG_LEAK# extensions in the file */
  int nchg_leak;

  /* arrays to hold CTE scaling data */
  double scalemjd[NUM_SCALE];
  double scaleval[NUM_SCALE];

  /* MJD values from CHG_LEAK extension headers */
  double chg_mjd1, chg_mjd2;

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

  /* read RN_CLIP keyword from primary header */
  if (GetKeyDbl(&hdr_ptr, "RN_CLIP", NO_DEFAULT, -999, &pars->rn_clip)) {
    trlerror("(pctecorr) Error reading RN_CLIP keyword from PCTETAB");
    status = KEYWORD_MISSING;
    return status;
  }

  /* read NSEMODEL keyword from primary header */
  if (GetKeyInt(&hdr_ptr, "NSEMODEL", NO_DEFAULT, -999, &pars->noise_model)) {
    trlerror("(pctecorr) Error reading NSEMODEL keyword from PCTETAB");
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

  /* read SUBTHRSH keyword from primary header */
  if (GetKeyDbl(&hdr_ptr, "SUBTHRSH", NO_DEFAULT, -999, &pars->sub_thresh)) {
    trlerror("(pctecorr) Error reading SUBTHRSH keyword from PCTETAB");
    status = KEYWORD_MISSING;
    return status;
  }

  /* read NCHGLEAK keyword from primary header.
   * descripes number of CHG_LEAK# extensions in the file */
  if (GetKeyInt(&hdr_ptr, "NCHGLEAK", NO_DEFAULT, -999, &nchg_leak)) {
    trlerror("(pctecorr) Error reading NCHGLEAK keyword from PCTETAB");
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
  /* read LEVEL data from second table extension */
  /* make filename + ext number 2 */
  sprintf(filename_wext, "%s[%i]", filename, 2);

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

  /* close CTE parameters file for extension 2
   * end of reading LEVEL data */
  c_tbtclo(tbl_ptr);
  /****************************************************************************/

  /****************************************************************************/
  /* read MJD/SCALE data from third table extension */
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
  /* get column pointer for the MJD points */
  c_tbcfnd1(tbl_ptr, mjdstr, &col_ptr_psin);
  if (c_iraferr() || col_ptr_psin == 0) {
    sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",mjdstr);
    trlerror(MsgText);
    status = COLUMN_NOT_FOUND;
    return status;
  }

  /* get column pointer for CTE scale */
  c_tbcfnd1(tbl_ptr, scalestr, &col_ptr_logq);
  if (c_iraferr() || col_ptr_logq == 0) {
    sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",scalestr);
    trlerror(MsgText);
    status = COLUMN_NOT_FOUND;
    return status;
  }

  /* iterate over table rows */
  for (j = 0; j < NUM_SCALE; j++) {
    /* get the MJD value */
    c_tbegtd(tbl_ptr, col_ptr_psin, j+1, &scalemjd[j]);
    if (c_iraferr()) {
      sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",j+1, mjdstr);
      trlerror(MsgText);
      status = TABLE_ERROR;
      return status;
    }

    /* read this scale value */
    c_tbegtd(tbl_ptr, col_ptr_logq, j+1, &scaleval[j]);
    if (c_iraferr()) {
      sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",
              j+1, scalestr);
      trlerror(MsgText);
      status = TABLE_ERROR;
      return status;
    }
  }

  /* close CTE parameters file for extension 3
   * end of reading MJD/SCALE data */
  c_tbtclo(tbl_ptr);

  /* calculate cte_frac */
   if ((status = CalcCteFrac(&pars->cte_frac, expstart, scalemjd, scaleval)))
       return status;
  /****************************************************************************/

  /****************************************************************************/
  /* read column by column scaling from 4th table extension */
  /* make filename + ext number 4 */
  sprintf(filename_wext, "%s[%i]", filename, 4);

  /* open CTE parameters file to extension number 4 */
  tbl_ptr = c_tbtopn(filename_wext, IRAF_READ_ONLY, 0);
  if (c_iraferr()) {
    sprintf(MsgText,"(pctecorr) Error opening %s with xtables",filename_wext);
    trlerror(MsgText);
    status = OPEN_FAILED;
    return status;
  }

  /* iterate over table rows */
  for (j = 0; j < AMP_COLS; j++) {
    /* loop over table columns to read log q values */
    for (k = 0; k < NAMPS; k++) {
      /* get column pointer for this log q */
      c_tbcfnd1(tbl_ptr, amp_keys[k], &col_ptr_amp);
      if (c_iraferr() || col_ptr_amp == 0) {
        sprintf(MsgText,"(pctecorr) Error getting column %s of PCTETAB",amp_keys[k]);
        trlerror(MsgText);
        status = COLUMN_NOT_FOUND;
        return status;
      }

      /* read this scale value */
      c_tbegtd(tbl_ptr, col_ptr_amp, j+1, &pars->col_scale[j*NAMPS + k]);
      if (c_iraferr()) {
        sprintf(MsgText,"(pctecorr) Error reading row %d of column %s in PCTETAB",
                j+1, amp_keys[k]);
        trlerror(MsgText);
        status = TABLE_ERROR;
        return status;
      }
    } /* end loop over columns */
  } /* end iterating over table rows */

  /* close CTE parameters file for extension */
  c_tbtclo(tbl_ptr);
  /* end of reading column by column scaling */
  /****************************************************************************/

  /****************************************************************************/
  /* read NODE/LOG_Q_# data from fifth+ table extension
   * there may be multiple CHG_LEAK extensions based on the time dependence
   * of the CTE trail profiles. we need to open them up until we find the one
   * that matches the time of observation and return that data.
   */

  for (l = 0; l < nchg_leak; l++) {
    /* make filename + ext number 4 */
    sprintf(filename_wext, "%s[%i]", filename, l+5);

    /* open CTE parameters file to extension number l+5 */
    tbl_ptr = c_tbtopn(filename_wext, IRAF_READ_ONLY, 0);
    if (c_iraferr()) {
      sprintf(MsgText,"(pctecorr) Error opening %s with xtables",filename_wext);
      trlerror(MsgText);
      status = OPEN_FAILED;
      return status;
    }

    chg_mjd1 = c_tbhgtr(tbl_ptr, "MJD1");
    chg_mjd2 = c_tbhgtr(tbl_ptr, "MJD2");

    /* check if we're in the date range for this table */
    if (chg_mjd1 <= expstart && chg_mjd2 > expstart) {
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
        } /* end loop over columns */
      } /* end iterating over table rows */

      /* close CTE parameters file for extension */
      c_tbtclo(tbl_ptr);
      break;
    } else {
      /* close CTE parameters file for extension */
      c_tbtclo(tbl_ptr);
    }
  }
  /* end of reading NODE/LOG_Q_# data */
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
int CompareCteParams(SingleGroup *x, ACSCTEParams *pars) {

  extern int status;

  /* functions from calacs/lib/key.c */
  int GetKeyInt (Hdr *hd, char *keyword, int use_def, int def, int *value);
  int GetKeyDbl (Hdr *hd, char *keyword, int use_def, double def, double *value);
  int PutKeyInt (Hdr *hd, char *keyword, int value, char *comment);
  int PutKeyDbl (Hdr *hd, char *keyword, double value, char *comment);

  int def_value = -99;
  double def_value_dbl = -99.0;

  double rn_clip;
  int noise_model;
  double sub_thresh;
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
  /* get value and perform check for read noise smoothing algorithm */
  if (GetKeyInt(x->globalhdr, "PCTENSMD", USE_DEFAULT, def_value, &noise_model)) {
    trlerror("(pctecorr) Error reading PCTENSMD keyword from image header");
    status = KEYWORD_MISSING;
    return status;
  }

  /* if the keyword value is not there or set to zero, update it with the value
   * read from the PCTETAB, otherwise use the read value in the CTE correction. */
  if (noise_model == def_value || noise_model == 0) {
    if (PutKeyInt(x->globalhdr, "PCTENSMD", pars->noise_model, "PCTE read noise smoothing algorithm")) {
      trlerror("(pctecorr) Error updating PCTENSMD keyword in image header");
      return (status = HEADER_PROBLEM);
    }
  } else {
    pars->noise_model = noise_model;
  }
  /****************************************************************************/

  /****************************************************************************/
  /* get value and perform check for over-subtraction threshold */
  if (GetKeyDbl(x->globalhdr, "PCTETRSH", USE_DEFAULT, def_value_dbl, &sub_thresh)) {
    trlerror("(pctecorr) Error reading PCTETRSH keyword from image header");
    status = KEYWORD_MISSING;
    return status;
  }

  /* if the keyword value is not there or set to zero, update it with the value
   * read from the PCTETAB, otherwise use the read value in the CTE correction. */
  if (sub_thresh == def_value || sub_thresh == 0) {
    if (PutKeyDbl(x->globalhdr, "PCTETRSH", pars->sub_thresh, "PCTE over subtraction threshold")) {
      trlerror("(pctecorr) Error updating PCTETRSH keyword in image header");
      return (status = HEADER_PROBLEM);
    }
  } else {
    pars->sub_thresh = sub_thresh;
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
 * worsening of CTE over time.
 */
int CalcCteFrac(double * cte_frac, const double expstart, const double scalemjd[NUM_SCALE], const double scaleval[NUM_SCALE])
{

  extern int status;

  /* iteration variables */
  int i;

  /* variables used to calculate the CTE scaling slope */
  double mjd_pt1 = 0;
  double mjd_pt2 = 0;      /* the MJD points at which the scaling is defined */
  double cte_pt1, cte_pt2; /* the CTE frac points at which the scaling is defined */
  Bool good2Go = False;

  /* find the values that bound this exposure */
  for (i = 0; i < NUM_SCALE-1; i++) {
    if (expstart >= scalemjd[i] && expstart < scalemjd[i+1]) {
      mjd_pt1 = scalemjd[i];
      mjd_pt2 = scalemjd[i+1];
      cte_pt1 = scaleval[i];
      cte_pt2 = scaleval[i+1];
      good2Go = True;
      break;
    }
  }

  /* it's possible this exposure is not bounded by any of defining points,
   * in that case we're extrapolating based on the last two points. */
  if (expstart >= scalemjd[NUM_SCALE-1] && mjd_pt1 == 0 && mjd_pt2 == 0)
  {
    mjd_pt1 = scalemjd[NUM_SCALE-2];
    mjd_pt2 = scalemjd[NUM_SCALE-1];
    cte_pt1 = scaleval[NUM_SCALE-2];
    cte_pt2 = scaleval[NUM_SCALE-1];
    good2Go = True;
  } else if (mjd_pt1 == 0 && mjd_pt2 == 0) {
    trlerror("(pctecorr) No suitable CTE scaling data found in PCTETAB");
    return (status = NO_GOOD_DATA);
  }

  if (!good2Go)
  {
      *cte_frac = 0; //Init anyhow
      trlerror("(pctecorr) No suitable CTE scaling fraction computed");
      return (status = INVALID_VALUE);
  }

  *cte_frac = ((cte_pt2 - cte_pt1) / (mjd_pt2 - mjd_pt1)) * (expstart - mjd_pt1) + cte_pt1;

  return status;
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
                    double dpde_l[NUM_LEV]) {

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

  return status;
}


/*
 * Attempt to separate readout noise from signal, since CTI happens before
 * readout noise is added to the signal.
 *
 * The clipping parameter read_noise controls the maximum amount by which a pixel
 * will be modified in each iteration. The noise_model parameter controls
 * which read noise smoothing algorithm is used, should be 0, 1, or 2.
 *
 * 0 = no RN correction
 * 1 = typical correction
 * 2 = hyper-conservative ; attribute as much as possible to RN
 */
int DecomposeRN(const int arrx, const int arry, const double data[arrx*arry],
                const double read_noise, const int noise_model,
                double sig_arr[arrx*arry], double noise_arr[arrx*arry]) {
  /* arguments:
     const int arrx                i: Number of rows in input array.
     const int arry                i: Number of columns in input array.
     const double data[arrx*arry]  i: Data from one amp.
     const double read_noise       i: RN_CLIP from PCTETAB (electrons).
     const int noise_model         i: NSEMODEL from PCTETAB.
     double sig_arr[arrx*arry]     o: Noiseless data to be corrected.
     double noise_arr[arrx*arry]   o: Noise to be added back after correction.
  */

  /* status variable for return */
  extern int status;

  /* iteration variables */
  int i, j, it_count;

  /* local constants */
  const int max_it1 = 25;
  const int max_it2 = 30;
  const double rms_fac = 1.10;
  const double f_fac = 1.25;

  /* accumulation variables */
  double sum;
  int num;

  double * local_noise;
  double d1, f, rms;

  /* check for valid noise_model */
  if (noise_model != 0 && noise_model != 1 && noise_model != 2) {
    return (status = ERROR_RETURN);
  }

  for (i = 0; i < arrx; i++) {
    for (j = 0; j < arry; j++) {
      sig_arr[i*arry + j] = data[i*arry + j];
      noise_arr[i*arry + j] = 0.0;
    }
  }

  /* no smoothing */
  if (noise_model == 0) {
    return status;
  }

  /* adjust each pixel to be more similar to the three pixels below it */
  it_count = 0;

  do {
    sum = 0.0;
    num = 0;

    for (i = 1; i < arrx; i++) {
      for (j = 0; j < arry; j++) {
        d1 = sig_arr[i*arry + j] - sig_arr[(i-1)*arry + j];

        if (d1 > 0.1*read_noise) {
          d1 = 0.1*read_noise;
        } else if (d1 < -0.1*read_noise) {
          d1 = -0.1*read_noise;
        }

        noise_arr[i*arry + j] += d1;

        sum += pow(noise_arr[i*arry + j], 2);
        num++;
      }
    }

    for (i = 0; i < arrx; i++) {
      for (j = 0; j < arry; j++) {
        sig_arr[i*arry + j] = data[i*arry + j] - noise_arr[i*arry + j];
      }
    }

    rms = sqrt(sum / (double) num);

    it_count++;

  } while (it_count < max_it1 || (it_count < max_it2 && rms < rms_fac*read_noise));

  if (noise_model == 1) {
    /* allocate local copy of noise */
    local_noise = (double *) malloc(arrx * arry * sizeof(double));

    /* remove any excessive local corrections
     */
    for (i = 0; i < arrx; i++) {
      for (j = 0; j < arry; j++) {
        f = 1.0;

        if (i > 1 && i < (arrx-1) && j > 0 && j < (arry-1)) {
          f = pow(noise_arr[(i-1) * arry + (j+1)], 2) +
              pow(noise_arr[ i    * arry + (j+1)], 2) +
	      pow(noise_arr[(i+1) * arry + (j+1)], 2) +
              pow(noise_arr[(i-1) * arry +  j   ], 2) +
              pow(noise_arr[ i    * arry +  j   ], 2) +
	      pow(noise_arr[(i+1) * arry +  j   ], 2) +
              pow(noise_arr[(i-1) * arry + (j-1)], 2) +
              pow(noise_arr[ i    * arry + (j-1)], 2) +
              pow(noise_arr[(i+1) * arry + (j-1)], 2);

          f = sqrt(f/9.0) / (f_fac * read_noise); /* scale-down factor */

          if (f < 1) {
            f = 1.00;
          }
        }

        local_noise[i*arry + j] = noise_arr[i*arry + j] / f;
      } /* end for j */
    } /* end for i */

    /* Copy scaled down noise to output */
    for (i = 0; i < arrx; i++) {
      for (j = 0; j < arry; j++) {
        noise_arr[i*arry + j] = local_noise[i*arry + j];
      }
    }

    free(local_noise);

  } /* end if noise_model==1 */

  /* Calculate noiseless data using final modeled noise.
     If noise_model==2, there is no source preservation.
   */
  for (i = 0; i < arrx; i++) {
    for (j = 0; j < arry; j++) {
      sig_arr[i*arry + j] = data[i*arry + j] - noise_arr[i*arry + j];
    }
  }

  return status;
}
