#include <stdio.h>
#include <string.h>
#include "c_iraf.h"
#include "imphttab.h"


/* Function to implement a self-test of this library */
int main(int argc, char **argv){

    char photmode[SZ_LINE];

	void c_irafinit (int, char **);
    void computeValues(char *);
    
    printf("==== Starting self-test for getphttab.c ====\n");

    status = 0;
	/* Initialize IRAF environment */
    c_irafinit(argc,argv);
    printf("==> Initialized IRAF environment for getphttab.\n");

    /* Test cases for interpolation half-way points for 2 par. values */
    strcpy(photmode,"acs,wfc1,fr853n#8344.75,mjd#54718.0");
    computeValues(photmode);
    strcpy(photmode,"acs,wfc1,fr853n#8307.4,mjd#54718.0");
    computeValues(photmode);
    strcpy(photmode,"acs,wfc1,mjd#54718.0,fr853n#8382.1");
    computeValues(photmode);
    strcpy(photmode,"acs,wfc1,fr853n#8344.75,mjd#55516.0");
    computeValues(photmode);
    strcpy(photmode,"acs,wfc1,mjd#53920.0,fr853n#8344.75");
    computeValues(photmode);
    strcpy(photmode,"acs,wfc1,fr853n#8382.1,mjd#53920.0");
    computeValues(photmode);

    /* Random test case of 2 parameterized variables */
    strcpy(photmode,"acs,wfc1,mjd#55021.0122,fr853n#8350.0");
    computeValues(photmode);
    
    strcpy(photmode,"acs,wfc1,f814w,mjd#55021.0122");
    computeValues(photmode);
}

void computeValues(char *photmode){

    PhotPar obs;
    int status;
    char refname[SZ_LINE], refpedigree[SZ_FITS_REC];
    
    strcpy(refname,"test_wfc1_imp.fits");
    strcpy(refpedigree,"Inflight calibrations");

    printf("\n\n==> PHOTMODE: %s\n",photmode);
    
    InitPhotPar(&obs, refname, refpedigree);
    printf("==> Finished InitPhotPar()\n");
    status = GetPhotTab(&obs, photmode);
    printf("==> Finished computing photometry keyword values.\n");
    FreePhotPar(&obs);
    printf("==> Cleaned up memory with FreePhotPar()\n");

    if (status == 0){
        printf("Test Results:\n");
        printf("  PHOTFLAM: %0.8g\n  PHOTPLAM: %g\n  PHOTBW: %g\n", obs.photflam, obs.photplam, obs.photbw);
    } else {
        printf("Did not successfully determine PHOT values.\n");
    }

}
