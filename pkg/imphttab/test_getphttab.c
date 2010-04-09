#include "c_iraf.h"
#include "imphttab.h"


/* Function to implement a self-test of this library */
int main(int argc, char **argv){
    int status;


    PhotPar obs;
    char photmode[SZ_LINE];
    char refname[SZ_LINE], refpedigree[SZ_FITS_REC];

	void c_irafinit (int, char **);
    
    printf("==== Starting self-test for getphttab.c ====\n");

    status = 0;
	/* Initialize IRAF environment */
    c_irafinit(argc,argv);
    printf("==> Initialized IRAF environment for getphttab.\n");

    strcpy(photmode,"acs,wfc1,f814w,mjd#55021.0122");
    strcpy(refname,"new_acs814_imp.fits");
    strcpy(refpedigree,"Inflight calibrations");
    
    printf("==> photmode: %s\n",&photmode);
    
    InitPhotPar(&obs, refname, refpedigree);
    printf("==> Finished InitPhotPar()\n");
    status = GetPhotTab(&obs, photmode);
    printf("==> Finished computing photometry keyword values.\n");
    FreePhotPar(&obs);
    printf("==> Cleaned up memory with FreePhotPar()\n");

    if (status == 0){
        printf("Test Results:\n");
        printf("  PHOTFLAM: %g\n  PHOTPLAM: %g\n  PHOTBW: %g\n", obs.photflam, obs.photplam, obs.photbw);
    } else {
        printf("Did not successfully determine PHOT values.\n");
    }
}
