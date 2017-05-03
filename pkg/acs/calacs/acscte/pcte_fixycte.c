#include <stdio.h>
#include <math.h>

# ifdef _OPENMP
#include <omp.h>
# endif

#include "pcte.h"
#include "acs.h"

/* function prototypes */
int sim_readout(const int arrx, double pix_cur[arrx], double pix_read[arrx],
                const double cte_frac_col[arrx], const int levels[NUM_LEV],
                const double dpde_l[NUM_LEV],
                const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
                const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV]);
int sim_readout_nit(const int arrx, double pix_cur[arrx], double pix_read[arrx],
                    const int shft_nit, const double cte_frac_col[arrx],
                    const int levels[NUM_LEV], const double dpde_l[NUM_LEV],
                    const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
                    const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV]);
int FixYCte(const int arrx, const int arry, const double sig_cte[arrx*arry],
            double sig_cor[arrx*arry], const int sim_nit, const int shft_nit,
            const double too_low, double cte_frac[arrx*arry],
            const int levels[NUM_LEV], const double dpde_l[NUM_LEV],
            const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
            const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV], int onecpu) {

  /* status variable for return */
  extern int status;

  /* iteration variables */
  int i, i2, j, n;

  /* arrays to hold columns of data */
  double pix_obs[arrx];
  double pix_cur[arrx];
  double pix_read[arrx];

  /* a column of the CTE scale array */
  double cte_frac_col[arrx];

  /* recalculated CTE scale, only needed in cases of over subtraction */
  double new_cte_frac, ncf_top, ncf_bot;

  /* flag for whether we've found a pixel with added charge */
  short int high_found;
  int high_location;

  /* track how many times we run the column. */
  short int redo_col = 0;
  int num_redo;

  /* Only use OpenMP, if specified by user and OpenMP was available for compilation */
  if (onecpu == 1) {

    trlmessage("Using single-CPU processing for YCTE correction.\n");

    /* loop over columns. columns are independent of each other. */
    for (j = 0; j < arry; j++) {
      /* copy column data */
      for (i = 0; i < arrx; i++) {
        pix_obs[i] = sig_cte[i*arry + j];//This is row major storage.
        pix_cur[i] = pix_obs[i];
        pix_read[i] = 0.0;

        cte_frac_col[i] = cte_frac[i*arry + j];
      }

      num_redo = 0;

      do {
        for (n = 0; n < sim_nit; n++) {
          status = sim_readout_nit(arrx, pix_cur, pix_read, shft_nit, cte_frac_col,
                                   levels, dpde_l, chg_leak_lt, chg_open_lt);

          if (status == 0) {
            for (i = 0; i < arrx; i++) {
              pix_cur[i] += pix_obs[i] - pix_read[i];
            }
          }
        }

        if (status == 0) {
          /* assume we won't have to redo this column */
          redo_col = 0;

          /* check this column for over subtracted pixels and maybe fix them. */
          for (i = 2; i < arrx-2; i++) {
            if (pix_cur[i] - pix_obs[i] < too_low &&
                pix_cur[i] < too_low && !redo_col) {
              high_found = 0;

              /* search for an upstream pixel with added charge */
              for (i2 = i-1; i2 > 0; i2--) {
                if (pix_cur[i2] - pix_obs[i2-1] < 0) {
                  high_found = 1;
                  high_location = i2 - 1;
                  break;
                }
              }

              /* if no added pixel was found then we can't do anything, move on */
              if (high_found == 0) {
                continue;
              } else {
                /* gonna have to redo this */
                redo_col = 1;
              }

              /* recalculate a new CTE scaling factor */
              ncf_top = fmax(pix_obs[i], 0.0);
              ncf_bot = ncf_top - pix_cur[i];

              if (ncf_top == 0) {
                new_cte_frac = 0.0;
              } else if (ncf_bot == 0) {
                new_cte_frac = 0.0;
              } else {
                new_cte_frac = ncf_top / ncf_bot;
              }

              /* distribute the new scaling factor */
              for (i2 = high_location; i2 <= i; i2++) {
                cte_frac_col[i2] *= new_cte_frac;

                if (cte_frac_col[i2] < 0) {
                  cte_frac_col[i2] = 0.0;
                }
              }

              if (i+1 < arrx) {
                cte_frac_col[i+1] *= 1.0 - 0.8 * (1.0 - new_cte_frac);
              }
              if (i+2 < arrx) {
                cte_frac_col[i+2] *= 1.0 - 0.6 * (1.0 - new_cte_frac);
              }
              if (i+3 < arrx) {
                cte_frac_col[i+3] *= 1.0 - 0.4 * (1.0 - new_cte_frac);
              }
              if (i+4 < arrx) {
                cte_frac_col[i+4] *= 1.0 - 0.2 * (1.0 - new_cte_frac);
              }

              if (redo_col) {
                break;
              }
            }
          }

          num_redo++;
        }
      } while (redo_col && num_redo < 10);

      if (status == 0) {
        /* copy fixed column to output */
        for (i = 0; i < arrx; i++) {
          sig_cor[i*arry + j] = pix_cur[i];
          cte_frac[i*arry + j] = cte_frac_col[i];
        }
      }
    } /* end loop over columns */

  } else {
    /* Use OpenMP for parallel processing, if specified by user */
#   ifdef _OPENMP
      trlmessage("Using parallel processing provided by OpenMP for YCTE correction.\n");
#   endif
#   ifndef _OPENMP
      trlmessage("Parallel processing for YCTE correction not used... OpenMP missing.\n");
#   endif
#   pragma omp parallel for schedule(dynamic) \
      private(i,j,n,status,cte_frac_col,new_cte_frac,ncf_top,ncf_bot,\
              high_found,high_location,redo_col,num_redo,pix_obs,pix_cur,pix_read) \
      shared(sig_cte,sig_cor,cte_frac)
    /* loop over columns. columns are independent of each other. */
    for (j = 0; j < arry; j++) {
      /* copy column data */
      for (i = 0; i < arrx; i++) {
        pix_obs[i] = sig_cte[i*arry + j];
        pix_cur[i] = pix_obs[i];
        pix_read[i] = 0.0;

        cte_frac_col[i] = cte_frac[i*arry + j];
      }

      num_redo = 0;

      do {
        for (n = 0; n < sim_nit; n++) {
          status = sim_readout_nit(arrx, pix_cur, pix_read, shft_nit, cte_frac_col,
                                   levels, dpde_l, chg_leak_lt, chg_open_lt);

          if (status == 0) {
            for (i = 0; i < arrx; i++) {
              pix_cur[i] += pix_obs[i] - pix_read[i];
            }
          }
        }

        if (status == 0) {
          /* assume we won't have to redo this column */
          redo_col = 0;

          /* check this column for over subtracted pixels and maybe fix them. */
          for (i = 2; i < arrx-2; i++) {
            if (pix_cur[i] - pix_obs[i] < too_low &&
                pix_cur[i] < too_low && !redo_col) {
              high_found = 0;

              /* search for an upstream pixel with added charge */
              for (i2 = i-1; i2 > 0; i2--) {
                if (pix_cur[i2] - pix_obs[i2-1] < 0) {
                  high_found = 1;
                  high_location = i2 - 1;
                  break;
                }
              }

              /* if no added pixel was found then we can't do anything, move on */
              if (high_found == 0) {
                continue;
              } else {
                /* gonna have to redo this */
                redo_col = 1;
              }

              /* recalculate a new CTE scaling factor */
              ncf_top = fmax(pix_obs[i], 0.0);
              ncf_bot = ncf_top - pix_cur[i];

              if (ncf_top == 0) {
                new_cte_frac = 0.0;
              } else if (ncf_bot == 0) {
                new_cte_frac = 0.0;
              } else {
                new_cte_frac = ncf_top / ncf_bot;
              }

              /* distribute the new scaling factor */
              for (i2 = high_location; i2 <= i; i2++) {
                cte_frac_col[i2] *= new_cte_frac;

                if (cte_frac_col[i2] < 0) {
                  cte_frac_col[i2] = 0.0;
                }
              }

              if (i+1 < arrx) {
                cte_frac_col[i+1] *= 1.0 - 0.8 * (1.0 - new_cte_frac);
              }
              if (i+2 < arrx) {
                cte_frac_col[i+2] *= 1.0 - 0.6 * (1.0 - new_cte_frac);
              }
              if (i+3 < arrx) {
                cte_frac_col[i+3] *= 1.0 - 0.4 * (1.0 - new_cte_frac);
              }
              if (i+4 < arrx) {
                cte_frac_col[i+4] *= 1.0 - 0.2 * (1.0 - new_cte_frac);
              }

              if (redo_col) {
                break;
              }
            }
          }

          num_redo++;
        }
      } while (redo_col && num_redo < 10);

      if (status == 0) {
        /* copy fixed column to output */
        for (i = 0; i < arrx; i++) {
          sig_cor[i*arry + j] = pix_cur[i];
          cte_frac[i*arry + j] = cte_frac_col[i];
        }
      }
    } /* end loop over columns */
  }

  return status;
}

/* call sim_readout shft_nit times as per Jay's new algorithm */
int sim_readout_nit(const int arrx, double pix_cur[arrx], double pix_read[arrx],
                    const int shft_nit, const double cte_frac_col[arrx],
                    const int levels[NUM_LEV], const double dpde_l[NUM_LEV],
                    const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
                    const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV]) {
  /* status variable for return */
  extern int status;

  /* iteration variables */
  int i,j;

  /* local container of column data */
  double pix_local[arrx];

  for (i = 0; i < arrx; i++) {
    pix_local[i] = pix_cur[i];
  }

  for (j = 0; j < shft_nit; j++) {
    status = sim_readout(arrx, pix_local, pix_read, cte_frac_col, levels,
                         dpde_l, chg_leak_lt, chg_open_lt);
    if (status != 0) {
      return status;
    }

    /* don't need to copy this back the last time through */
    if (j < shft_nit - 1) {
      for (i = 0; i < arrx; i++) {
        pix_local[i] = pix_read[i];
      }
    }
  }

  return status;
}


/* workhorse function that moves simulates CTE by shifting charge down the column
 * and keeping track of charge trapped and released */
int sim_readout(const int arrx, double pix_cur[arrx], double pix_read[arrx],
                const double cte_frac_col[arrx], const int levels[NUM_LEV],
                const double dpde_l[NUM_LEV],
                const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
                const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV]) {

  /* status variable for return */
  extern int status;

  /* iteration variables */
  int i, l;

  /* holds some trap info, I guess */
  double ftrap_l[NUM_LEV];
  int ttrap_l[NUM_LEV];

  double pix0, pix1;    /* current pixel containers */
  double ffil;          /* fraction of this trap that gets filled */

  double add_charge1, add_charge2;
  double rem_charge, rem_charge_temp;

  /* initialize traps with nothing in them */
  for (l = 0; l < NUM_LEV; l++) {
    ftrap_l[l] = 0.0;
    ttrap_l[l] = 999;
  }

  /* iterate over every pixel in the column. each pixel gets changed by
   * charge being added to it by trap releases and losing charge to traps. */
  for (i = 0; i < arrx; i++) {
    pix0 = pix_cur[i];

    if (i > 0) {
      if (cte_frac_col[i] < cte_frac_col[i-1]) {
        for (l = 0; l < NUM_LEV; l++) {
          ftrap_l[l] *= cte_frac_col[i] / cte_frac_col[i-1];
        }
      }
    }

    add_charge1 = 0.0;

    for (l = 0; l < NUM_LEV; l++) {
      if (ttrap_l[l] < MAX_TAIL_LEN) {
        ttrap_l[l]++;

        add_charge1 += chg_leak_lt[(ttrap_l[l]-1)*NUM_LEV + l] * ftrap_l[l];
      }
    }

    pix1 = pix0 + add_charge1;

    add_charge2 = 0.0;
    rem_charge = 0.0;

    for (l = 0; l < NUM_LEV - 1; l++) {
      /* skip the rest of the levels if we don't have enough charge to reach
       * any more traps */
      if (pix1 < levels[l]) {
        break;
      }

      /* can usually fill an entire trap, but if not calculate the fraction */
      ffil = 1.000;

      if (pix1 < levels[l+1]) {
        ffil = (pix1 - levels[l]) / (levels[l+1] - levels[l]);
      }

      rem_charge_temp = ffil * dpde_l[l+1] * cte_frac_col[i];

      if (ttrap_l[l] <= MAX_TAIL_LEN) {
        add_charge2 += chg_open_lt[(ttrap_l[l]-1)*NUM_LEV + l] * ftrap_l[l];
      }

      ttrap_l[l] = 0;

      ftrap_l[l] = rem_charge_temp;

      rem_charge += rem_charge_temp;
    }

    pix_read[i] = pix_cur[i] + add_charge1 + add_charge2 - rem_charge;
  }

  return status;
}
