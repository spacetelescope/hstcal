#include <stdio.h>
#include <string.h>
#include "hstcal.h"
#include "c_iraf.h"
#include "imphttab.h"


/* Function to implement a self-test of this library */
int main(int argc, char **argv){
    void c_irafinit (int, char **);
    void compute_values(char [], char [], double, float, float, FILE *);
    char *outfile;
    FILE *fp;

    outfile = "testphot_out.txt";
    fp = fopen(outfile, "w");

    fprintf(fp, "==== Starting self-test for getphttab.c ====\n");

    /* Initialize IRAF environment */
    c_irafinit(argc, argv);
    fprintf(fp, "==> Initialized IRAF environment for getphttab.\n");

    compute_values(
        "acs,wfc1,fr853n#8344.75,mjd#54718",
        "Half-way points for both pars.",
        1.2695867062843802e-18, 8344.119572140797, 56.793040505098716, fp);

    compute_values(
        "acs,wfc1,mjd#54718.0,fr853n#8344.75",
        "Like above but in different order (results should be same).",
        1.2695867062843802e-18, 8344.119572140797, 56.793040505098716, fp);

    compute_values(
        "acs,wfc1,fr853n#8158,mjd#52334",
        "Both pars at lower boundary.",
        1.1984913607659768e-18, 8158.271171088088, 59.462191932802028, fp);

    compute_values(
        "acs,wfc1,fr853n#8158,mjd#54718",
        "Only FR at lower boundary.",
        1.2230171665655153e-18, 8158.273423487024, 59.320487871226852, fp);

    compute_values(
        "acs,wfc1,fr853n#8344.75,mjd#59214",
        "Only MJD at upper boundary.",
        1.2544879212041671e-18, 8344.116261683303, 56.820520968140734, fp);

    compute_values(
        "acs,wfc1,fr853n#8905,mjd#59214",
        "Both pars at upper boundary.",
        1.522786245291286e-18, 8900.34723883603, 77.771312111848076, fp);

    compute_values(
        "acs,wfc1,f625w,fr462n#4462,mjd#55000",
        "Random 2 pars with another filter.",
        1.0157653348302224e-13, 5873.489318646888, 926.61307046259992, fp);

    compute_values(
        "acs,wfc1,f814w,mjd#55021.0122",
        "Random 1 par with a filter.",
        7.0073064555265768e-20, 8059.349015606429, 652.38059499466738, fp);

    compute_values(
        "acs,wfc1,f814w",
        "No par with a filter.",
        7.0767633114602767e-20, 8058.784799323781, 652.34210160962766, fp);

    /* Extrapolation is not direct equivalent with pysynphot.
       Here, extrapolate from last 2 data points with straight line. */

    compute_values(
        "acs,wfc1,fr853n#6500,mjd#54718",
        "FR extrapolation (should not be allowed).",
        -9999, -9999, -9999, fp);

    compute_values(
        "acs,wfc1,fr853n#8344.75,mjd#48000",
        "MJD extrapolation into the past.",
        1.2339682890817842e-18, 8344.08667342431, 57.038104938309452, fp);

    compute_values(
        "acs,wfc1,fr853n#8344.75,mjd#68000",
        "MJD extrapolation into the future.",
        1.2544879212041671e-18, 8344.116261683303, 56.820520968140734, fp);

    fclose(fp);
    printf("Wrote test results to %s\n", outfile);

    return 0;
}


void compute_values(char photmode[CHAR_LINE_LENGTH], char test_title[CHAR_LINE_LENGTH],
                    double ans_photflam, float ans_photplam, float ans_photbw,
                    FILE *fp) {
    PhotPar obs;
    int status;
    char refname[CHAR_LINE_LENGTH], refpedigree[SZ_FITS_REC];
    double diff_photflam, diff_photplam, diff_photbw;

    /*strcpy(refname,"/grp/hst/cdbs/jref/w3m17171j_imp.fits");*/  /* -9999 */
    strcpy(refname, "w3m17171j_extrap_imp.fits");  /* Extrapolate */
    strcpy(refpedigree, "Inflight calibrations");

    fprintf(fp, "\n==> PHOTMODE: %s\n==> %s\n", photmode, test_title);

    InitPhotPar(&obs, refname, refpedigree);
    fprintf(fp, "==> Finished InitPhotPar()\n");
    status = GetPhotTab(&obs, photmode);
    fprintf(fp, "==> Finished computing photometry keyword values.\n");
    FreePhotPar(&obs);
    fprintf(fp, "==> Cleaned up memory with FreePhotPar()\n");

    if (status == 0){
        diff_photflam = (obs.photflam - ans_photflam) * 100.0 / ans_photflam;
        diff_photplam = (obs.photplam - ans_photplam) * 100.0 / ans_photplam;
        diff_photbw = (obs.photbw - ans_photbw) * 100.0 / ans_photbw;

        fprintf(fp, "Test Results:\n");
        fprintf(fp, "  PHOTFLAM: %0.8g (%0.8g, pct_diff=%g)\n",
                obs.photflam, ans_photflam, diff_photflam);
        fprintf(fp, "  PHOTPLAM: %g (%g, pct_diff=%g)\n",
                obs.photplam, ans_photplam, diff_photplam);
        fprintf(fp, "  PHOTBW: %g (%g, pct_diff=%g)\n",
                obs.photbw, ans_photbw, diff_photbw);
    } else {
        fprintf(fp, "Did not successfully determine PHOT values.\n");
    }
}
