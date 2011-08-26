#include <stdio.h>
#include <math.h>

# ifdef _OPENMP
#include <omp.h>
# endif

#include "pcte.h"
#include "acs.h"

/* function prototypes */
int sim_readout(const int arrx, double pix_cur[arrx], double pix_read[arrx],
                const double cte_frac, const int levels[NUM_LEV],
                const double dpde_l[NUM_LEV], const int tail_len[NUM_LEV],
                const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
                const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV]);
int sim_readout_nit(const int arrx, double pix_cur[arrx], double pix_read[arrx],
                    const double cte_frac, const int shft_nit, const int levels[NUM_LEV],
                    const double dpde_l[NUM_LEV], const int tail_len[NUM_LEV],
                    const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
                    const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV]);


int FixYCte(const int arrx, const int arry, const double sig_cte[arrx*arry],
            double sig_cor[arrx*arry], const double cte_frac, const int sim_nit,
            const int shft_nit, const int levels[NUM_LEV],
            const double dpde_l[NUM_LEV], const int tail_len[NUM_LEV],
            const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
            const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV], int onecpu) {
  
  /* status variable for return */
  extern int status;
  
  /* iteration variables */
  int i, j, n;
  
  /* arrays to hold columns of data */
  double pix_obs[arrx];
  double pix_cur[arrx];
  double pix_read[arrx];
    
  /* Only use OpenMP, if specified by user and OpenMP was available for compilation */
  if (onecpu == 1) {

      trlmessage("Using single-CPU processing for YCTE correction.\n");
      /* loop over columns. columns are independent of each other. */
      for (i = 0; i < arry; i++) {    
        /* copy column data */
        for (j = 0; j < arrx; j++) {
          pix_obs[j] = sig_cte[j*arry + i];
          pix_cur[j] = pix_obs[j];
          pix_read[j] = 0.0;
        }

        for (n = 0; n < sim_nit; n++) {
          status = sim_readout_nit(arrx, pix_cur, pix_read, cte_frac, shft_nit,
                                   levels, dpde_l, tail_len, chg_leak_lt, chg_open_lt);
          if (status == 0) {

              for (j = 0; j < arrx; j++) {
                pix_cur[j] += pix_obs[j] - pix_read[j];
              }
          } /* Only do this if sim_readout_nit() succeeded */

        }
        if (status == 0){
            /* copy fixed column to output */
            for (j = 0; j < arrx; j++) {
              sig_cor[j*arry + i] = pix_cur[j];
            }
        } /* Only do this if sim_readout_nit() succeeded */
      } /* end loop over columns */
  } else {
    /* Use OpenMP for parallel processing, if specified by user */
#   ifdef _OPENMP
      trlmessage("Using parallel processing provided by OpenMP for YCTE correction.\n");
#   endif
#   ifndef _OPENMP
      trlmessage("Parallel processing for YCTE correction not used... OpenMP missing.\n"); 
#   endif
#   pragma omp parallel for schedule(dynamic) private(i,j,n,status,pix_obs,pix_cur,pix_read) shared(sig_cte,sig_cor)
    /* loop over columns. columns are independent of each other. */
    for (i = 0; i < arry; i++) {
      /* copy column data */
      for (j = 0; j < arrx; j++) {
        pix_obs[j] = sig_cte[j*arry + i];
        pix_cur[j] = pix_obs[j];
        pix_read[j] = 0.0;
      }
      
      for (n = 0; n < sim_nit; n++) {
        status = sim_readout_nit(arrx, pix_cur, pix_read, cte_frac, shft_nit,
                                 levels, dpde_l, tail_len, chg_leak_lt, chg_open_lt);
        if (status == 0) {
          for (j = 0; j < arrx; j++) {
            pix_cur[j] += pix_obs[j] - pix_read[j];
          }
        } /* Only do this if sim_readout_nit() succeeded */
      }
      
      if (status == 0){
        /* copy fixed column to output */
        for (j = 0; j < arrx; j++) {          
          sig_cor[j*arry + i] = pix_cur[j];
        }
      } /* Only do this if sim_readout_nit() succeeded */
    } /* end loop over columns */
  }  
  return status;
}

/* call sim_readout shft_nit times as per Jay's new algorithm */
int sim_readout_nit(const int arrx, double pix_cur[arrx], double pix_read[arrx],
                    const double cte_frac, const int shft_nit, const int levels[NUM_LEV],
                    const double dpde_l[NUM_LEV], const int tail_len[NUM_LEV],
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
    status = sim_readout(arrx, pix_local, pix_read, cte_frac, levels, dpde_l,
                         tail_len, chg_leak_lt, chg_open_lt);
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
                const double cte_frac, const int levels[NUM_LEV],
                const double dpde_l[NUM_LEV], const int tail_len[NUM_LEV],
                const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
                const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV]) {
  
  /* status variable for return */
  extern int status;
  
  /* iteration variables */
  int i,l,t;
  int tmax;
  
  /* holds some trap info, I guess */
  double ftrap_lj[arrx*NUM_LEV];
  
  double pix0;    /* current pixel container */
  double fpix;    /* fraction of this pixel involved */
  double fopn;    /* fraction of this trap that is open */
  double ffil;    /* fraction of this trap that gets filled */
  double dpix;    /* amount of charge that gets transferred */
  
  /* copy input to output */
  for (i = 0; i < arrx; i++) {
    pix_read[i] = pix_cur[i];
    
    /* initialize trap array to zeros */
    for (l = 0; l < NUM_LEV; l++) {
      ftrap_lj[i*NUM_LEV + l] = 0.0;
    }
  }
  
  /* iterate over every pixel in the column */
  for (i = 0; i < arrx; i++) {
    pix0 = pix_read[i];
    
    for (l = 1; l < NUM_LEV; l++) {
      /* skip the rest of the levels if we don't have enough charge to reach
       * any more traps */
      if (pix0 < levels[l-1]) {
        continue;
      }
      
      /* can usually fill an entire trap, but if not calculate the fraction */
      fpix = 1.000;
      if (pix0 < levels[l]) {
        fpix = (pix0 - levels[l-1]) / (levels[l] - levels[l-1]);
      }
      
      /* how much of trap is available for filling */
      fopn = 1.0 - ftrap_lj[i*NUM_LEV + l];
      
      /* what fraction of the trap can I fill given current conditions? */
      ffil = fpix*fopn;
      
      /* how many electrons can this take? */
      dpix = cte_frac * dpde_l[l] * ffil * ((double) (i+1) / (double) CTE_REF_ROW);
      
      /* remove electrons from the pixel */
      pix_read[i] -= dpix;
      
      /* redistribute electrons in the tail */
      if ((i + tail_len[l]) < arrx) {
        tmax = tail_len[l];
      } else {
        tmax = arrx - i - 1;
      }
      
      for (t = 1; t <= tmax; t++) {
        pix_read[i+t] += (dpix * chg_leak_lt[(t-1)*NUM_LEV + l]);
        ftrap_lj[(i+t)*NUM_LEV + l] += (ffil * chg_open_lt[(t-1)*NUM_LEV + l]);
      }
    }
  }
  
  return status;
}
