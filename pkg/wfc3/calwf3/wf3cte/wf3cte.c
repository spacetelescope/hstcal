/* WFC3 -- CTE loss correction for UVIS

   M. sosey  Aug-2014  Adapted for the pipeline from Jay Andersons CTE correction code for wfc3 UVIS
   raw2raz_wfc3uv.F , an edited file was delivered december 2014, and both are different from the
   fortran code currently served on the wfc3 website.

   M. Sosey Aug-2016 Adapted to be used with Subarrays as well as full frame arrays,
   as long as the subarray contains physical overscan pixels, which don't include the science team subarrays
   which can span quads.

   M.D. De La Pena Dec-2019: This routine has been significantly upgraded by Jay Anderson (JA) and delivered
   in November 2019.  As JA is the important resource for this algorithm, I have only cleaned up the comments
   in his original delivered version, fixed up brace placement, and created defines for some of the 
   hard-coded values.  Minimal changes were done explicitly to keep the code in a form familiar to JA for 
   possible future modifications.

   M.D. De La Pena Mar-2020: Further changes to accommodate subarrays - only evaluate valid (non-zero) pixels.
   Updates received from Jay Anderson. Removed deprecated routines: find_dadj, rsz2rsc, inverse_cte_blur, and raz2rsz.
   Small bug found in original subarray code during testing.

   M.D. De La Pena Apr-2021: Fix to address a problem detected when processing a Tungsten flat with a high
   background.  Uninitialized values were used for further computation causing an eventual exception. 
*/

# include <time.h>
# include <string.h>
# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <float.h>

# ifdef _OPENMP
#  include <omp.h>
# endif

# include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "wf3corr.h"
# include "cte.h"
# include "trlbuf.h"

/* 
    These are defined in wf3.h.
    NAMPS       4
    RAZ_COLS 8412
    RAZ_ROWS 2070
*/

# define NITMAX 299          /* Maximum number of iterations */
# define P_OVRSCN 25         /* Physical overscan */
# define V_OVRSCNX2 60       /* Virtual overscan x 2 */
# define XAMP_SCI_DIM 2048   /* X dimension of each AMP of science pixels */
# define YAMP_SCI_DIM 2051   /* Y dimension of each AMP of science pixels */
# define WsMAX 999           /* Maximum number of traps */

/* Used in find_raz2rnoival */
# define SPREAD_FOR_HISTO 4.5
# define LOW_CLIP 3.75
# define HIGH_CLIP 9.75
# define NUM_BINS 1001

int sub_ctecor_v2c(float *, 
                   float *,
                   int,  
                   double *,
                   double *,
                   float *, 
                   float *, 
                   float,
                   float,
                   int,
                   int,
                   float *); 

float find_raz2rnoival(float *, 
                       float *, 
                       float *);

int WF3cte (char *input, char *output, CCD_Switch *cte_sw,
        RefFileInfo *refnames, int printtime, int verbose, int onecpu) {

    /*
    input: filename
    output: filename
    cte_sw: the calibration flags
    refnames: the names of the calibration reference files
    onecpu: use parallel processing?

    The following are new primary header keywords which will be added to the data
    so that they can be updated by the code. They are also specified in the PCTETAB
    reference file.

    These are taken from the PCTETAB
    CTE_NAME - name of cte algorithm
    CTE_VER - version number of cte algorithm
    CTEDATE0 - date of wfc3/uvis installation in HST, in MJD
    CTEDATE1 - reference date of CTE model pinning, in MJD

    PCTETLEN - max length of CTE trail
    PCTESMIT - number of iterations used in CTE forward modeling
    PCTESHFT - number of iterations used in the parallel transfer
    PCTENSMD - readnoise mitigation algorithm
    PCTETRSH - over-subtraction threshold
    PCTEFRAC - cte scaling frac calculated from expstart
    PCTERNOI - the readnoise clipping level to use ***NOTE: This value is no longer used
               from the PCTETAB.  If PCTERNOI keyword value in the raw science image 
               header is non-zero, it will be used for the CTE computations.  Otherwise,
               the value is computed on-the-fly based upon the raw image data. (March 2020)

    #These are taken from getreffiles.c
    DRKCFILE is a new dark reference file used only in the CTE branch *_DRC.fits
    BIACFILE is a new super-bias reference file used only in the CTE branch *_BIC.fits
    PCTETAB is a new reference file FITS table which will contain the software parameter switches for the CTE correction *_CTE.fit

    This is the main workhorse function for removing the CTE from WFC3 UVIS images

    Unfortunately this happens before anything else in wfc3, so there's a lot of reading files
    at the beginning in order to populate needed information. The rest of the pipeline works
    on one chip at a time and the structures are all defined to support that. None of these
    structures are defined until the code enters the single chip loops. This differs from the
    CTE correction in ACS which occurs later in the process after basic structures are defined.
*/

    extern int status;

    WF3Info wf3; /*structure with calibration switches and reference files for passing*/
    Hdr phdr;    /*primary header for input image, all output information saved here*/
    Hdr scihdr;  /*science header in case of subarray image to detect chip*/
    IODescPtr ip = NULL;

    CTEParams cte_pars; /*STRUCTURE HOLDING THE MODEL PARAMETERS*/
    SingleGroup cd; /*SCI 1, chip 2*/
    SingleGroup ab; /*SCI 2, chip 1*/
    SingleGroup subcd; /*subarray chip*/
    SingleGroup subab; /*subarray chip*/
    SingleGroup raz; /* THE LARGE FORMAT COMBINATION OF CDAB*/
    SingleGroup rsz; /* LARGE FORMAT READNOISE CORRECTED IMAGE */
    SingleGroup rsc; /* CTE CORRECTED*/
    SingleGroup rzc; /* FINAL CTE CORRECTED IMAGE */
    SingleGroup chg; /* THE CHANGE DUE TO CTE  */
    SingleGroup raw; /* THE RAW IMAGE IN RAZ FORMAT */

    int i,j; /*loop vars*/
    int max_threads=1;
    clock_t begin;
    double  time_spent;
    float hardset=0.0;

    /* These are used to find subarrays with physical overscan */
    int sci_bin[2];			/* bin size of science image */
    int sci_corner[2];		/* science image corner location */
    int ref_bin[2];
    int ref_corner[2];
    int rsize = 1;          /* reference pixel size */
    int start=0;            /*where the subarray starts*/
    int finish=0;           /*where the subarray ends*/

    /* init header vars */
    initHdr(&phdr);
    initHdr(&scihdr);

    float readNoise = 0.0;

    int ret;

    /*check if this is a subarray image.
      This is necessary because the CTE routine will start with the raw images
      from scratch and read them in so that both chips can be used. CTE is
      outside of the normal processing where one chip goes through the pipeline
      at a time, both chips are used at the same time for the correction.

      For the case of subarrays, a fake second chip needs to be created.
      The subarray is also placed inside the confines of a full size image
      and a mask is created to ignore pixels not associated with the original
      data during the cte correction. This is necessary because the pixel location
      itself is used as part of the correction. A secondary option would be to set
      the looping arrays to variable sizes and make sure all array references were
      consistent with the current data being processed. I decided on masking which
      might allow for other considerations in future updates.

      Only subarrays which were taken with physical overscan pixels are currently valid
      This distinction can be made with the CRDS ruleset for PCTECORR but it
      should also be checked here incase users update the header themselves for
      local runs. In order to check for overscan pixels I'm using the array start
      location instead of the APERTURE keyword information (there are known user
      apertures which do not have overscan pixels, but this gets around string
      comparisons and any future name changes or aperture additions in the future)
     */
    begin = (double)clock();

    /*CONTAIN PARALLEL PROCESSING TO A SINGLE THREAD AS USER OPTION*/
#   ifdef _OPENMP
    trlmessage("Using parallel processing provided by OpenMP inside CTE routine");
    if (onecpu){
        omp_set_dynamic(0);
        max_threads=1;
        sprintf(MsgText,"onecpu == TRUE, Using only %i threads/cpu", max_threads);
    } else {
        omp_set_dynamic(0);
        max_threads = omp_get_num_procs(); /*be nice, use 1 less than avail?*/
        sprintf(MsgText,"Setting max threads to %i of %i cpus",max_threads, omp_get_num_procs());
    }
    omp_set_num_threads(max_threads);
    trlmessage(MsgText);
#   endif


    /* COPY COMMAND-LINE ARGUMENTS INTO WF3. */
    WF3Init (&wf3); /*sets default information*/
    strcpy (wf3.input, input);
    strcpy (wf3.output, output);

    PrBegin ("WFC3CTE");
    if (wf3.printtime)
        TimeStamp("WFC3CTE Started: ",wf3.rootname);

    /* CHECK WHETHER THE OUTPUT FILE ALREADY EXISTS. */
    if (FileExists (wf3.output)){
        WhichError(status);
        return (ERROR_RETURN);
    }

    wf3.pctecorr = cte_sw->pctecorr;
    wf3.darkcorr = cte_sw->darkcorr;
    wf3.biascorr = cte_sw->biascorr;
    wf3.blevcorr = cte_sw->blevcorr;
    wf3.printtime = printtime;
    wf3.verbose = verbose;
    wf3.refnames = refnames;

    PrFileName ("input", wf3.input);
    PrFileName ("output", wf3.output);

    if (wf3.biascorr == COMPLETE){
        trlmessage("BIASCORR complete for input image, CTE can't be performed");
        return(ERROR_RETURN);
    }
    if (wf3.darkcorr == COMPLETE){
        trlmessage("DARKCORR complete for input image, CTE can't be performed");
        return(ERROR_RETURN);
    }
    if (wf3.blevcorr == COMPLETE){
        trlmessage("BLEVCORR complete for input image, CTE can't be performed");
        return(ERROR_RETURN);
    }

    /* DETERMINE THE NAMES OF THE TRAILER FILES BASED ON THE INPUT
       AND OUTPUT FILE NAMES, THEN INITIALIZE THE TRAILER FILE BUFFER
       WITH THOSE NAMES.
       */
    if (initCTETrl (input, output))
        return (status);

    /* OPEN INPUT IMAGE IN ORDER TO READ ITS PRIMARY HEADER. */
    if (LoadHdr (wf3.input, &phdr) ){
        WhichError(status);
        return (ERROR_RETURN);
    }

    /* GET KEYWORD VALUES FROM PRIMARY HEADER. */
    if (GetKeys (&wf3, &phdr)) {
        freeHdr (&phdr);
        return (status);
    }

    if (GetCTEFlags (&wf3, &phdr)) {
        freeHdr(&phdr);
        return (status);
    }


    /*SET UP THE ARRAYS WHICH WILL BE PASSED AROUND*/
    initSingleGroup(&raz);
    allocSingleGroup(&raz, RAZ_COLS, RAZ_ROWS, True);

    initSingleGroup(&rsz);
    allocSingleGroup(&rsz, RAZ_COLS, RAZ_ROWS, True);

    initSingleGroup(&rsc);
    allocSingleGroup(&rsc, RAZ_COLS, RAZ_ROWS, True);

    initSingleGroup(&rzc);
    allocSingleGroup(&rzc, RAZ_COLS, RAZ_ROWS, True);

    initSingleGroup(&raw);
    allocSingleGroup(&raw, RAZ_COLS, RAZ_ROWS, True);

    initSingleGroup(&chg);
    allocSingleGroup(&chg, RAZ_COLS, RAZ_ROWS, True);

    /*hardset the science arrays*/
    for (i=0;i<RAZ_COLS;i++){
        for(j=0;j<RAZ_ROWS;j++){
            Pix(raw.sci.data,i,j)=hardset;
            Pix(raz.sci.data,i,j)=hardset;
            Pix(rsz.sci.data,i,j)=hardset;
            Pix(rsc.sci.data,i,j)=hardset;
            Pix(rzc.sci.data,i,j)=hardset;
            Pix(chg.sci.data,i,j)=hardset;
        }
    }

    /*READ IN THE CTE PARAMETER TABLE*/
    initCTEParams(&cte_pars);
    if (GetCTEPars (wf3.pctetab.name, &cte_pars))
        return (status);

    if (verbose){
        PrRefInfo ("pctetab", wf3.pctetab.name, wf3.pctetab.pedigree,
                wf3.pctetab.descrip, wf3.pctetab.descrip2);
    }

    /* Full frame and subarrays always have group 1
       If it's a subarray, the group can be from either chip
       and will still be labled group 1 because it's the FIRST
       and only group, so look at the ccdchip instead.

       amps ab are in chip1, sci,2
       amps cd are in chip2, sci,1

    */
    if (wf3.subarray) {
        /* OPEN INPUT IMAGE IN ORDER TO READ ITS SCIENCE HEADER. */
        ip = openInputImage (wf3.input, "SCI", 1);
        if (hstio_err()) {
            sprintf (MsgText, "Image: \"%s\" is not present", wf3.input);
            trlerror (MsgText);
            return (status = OPEN_FAILED);
        }
        getHeader (ip, &scihdr);
        if (ip != NULL)
            closeImage (ip);

        /* Get CCD-specific parameters. */
        if (GetKeyInt (&scihdr, "CCDCHIP", USE_DEFAULT, 1, &wf3.chip)){
            freeHdr(&scihdr);
            return (status);
        }
        freeHdr(&scihdr);

        if (wf3.chip == 2){ /*sci1,cd*/
            start=0;
            finish=0;
            /*get CD subarray from first extension*/
            initSingleGroup (&subcd);
            getSingleGroup (wf3.input, 1, &subcd);
            if (hstio_err()){
                freeSingleGroup(&subcd);
                return (status = OPEN_FAILED);
            }

            /*create an empty full size chip for pasting*/
            initSingleGroup(&cd);
            allocSingleGroup(&cd,RAZ_COLS/2,RAZ_ROWS, True);
            cd.group_num=1;
            CreateEmptyChip(&wf3, &cd);

            if (GetCorner(&subcd.sci.hdr, rsize, sci_bin, sci_corner))
                return (status);
            if (GetCorner(&cd.sci.hdr, rsize, ref_bin, ref_corner))
                return (status);

            start = sci_corner[0] - ref_corner[0];
            finish = start + subcd.sci.data.nx;
            if ( start >= P_OVRSCN &&  finish + V_OVRSCNX2 <= (RAZ_COLS/2) - P_OVRSCN){
                sprintf(MsgText,"Subarray not taken with physical overscan (%i %i)\nCan't perform CTE correction\n",start,finish);
                trlmessage(MsgText);
                return(ERROR_RETURN);
            }

            /*SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
              AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
              CTE CODE VARIABLES.
              */
            if (CompareCTEParams(&subcd, &cte_pars))
                return (status);

            /*Put the subarray data into full frame*/
            Sub2Full(&wf3, &subcd, &cd, 0, 1, 1);

            /* now create an empty chip 1*/
            initSingleGroup(&ab);
            allocSingleGroup(&ab,RAZ_COLS/2,RAZ_ROWS, True);
            ab.group_num=2;
            CreateEmptyChip(&wf3, &ab);

            /* SAVE A COPY OF THE RAW IMAGE BEFORE BIAS FOR LATER */
            makeRAZ(&cd,&ab,&raw);

            /* Subtract the BIAC file from the subarray before continuing
               The bias routine will take care of cutting out the correct
               image location for the subarray.*/

            if (doCteBias(&wf3,&subcd)){
                freeSingleGroup(&subcd);
                return(status);
            }

            /*reset the array after bias subtraction*/
            Sub2Full(&wf3, &subcd, &cd, 0, 1, 1);


        } else { /*chip is 1, ab, sci2*/
            start=0;
            finish=0;
            initSingleGroup(&subab);
            getSingleGroup(wf3.input, 1, &subab);
            if (hstio_err()){
                freeSingleGroup(&subab);
                return (status = OPEN_FAILED);
            }

            /*make an empty fullsize chip for pasting*/
            initSingleGroup(&ab);
            allocSingleGroup(&ab,RAZ_COLS/2,RAZ_ROWS, True);
            ab.group_num=2;
            CreateEmptyChip(&wf3, &ab);

            if ( GetCorner(&subab.sci.hdr, rsize, sci_bin, sci_corner))
                return (status);

            if ( GetCorner(&ab.sci.hdr, rsize, ref_bin, ref_corner))
                return (status);

            start = sci_corner[0] - ref_corner[0];
            finish = start + subab.sci.data.nx;
            if ( start >= P_OVRSCN &&  finish + V_OVRSCNX2 <= (RAZ_COLS/2) - P_OVRSCN){
                sprintf(MsgText,"Subarray not taken with physical overscan (%i %i)\nCan't perform CTE correction\n",start,finish);
                trlmessage(MsgText);
                return(ERROR_RETURN);
            }
            /*add subarray to full frame image*/
            Sub2Full(&wf3, &subab, &ab, 0, 1, 1);

            /*SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
              AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
              CTE CODE VARIABLES.
              */
            if (CompareCTEParams(&subab, &cte_pars))
                return (status);

            /* now create an empty chip 2*/
            initSingleGroup(&cd);
            allocSingleGroup(&cd,RAZ_COLS/2,RAZ_ROWS, True);
            cd.group_num=1;
            CreateEmptyChip(&wf3, &cd);

            /* SAVE A COPY OF THE RAW IMAGE FOR LATER */
            makeRAZ(&cd,&ab,&raw);

            /* Subtract the BIAC file from the subarray before continuing*/
            subab.group_num=2;
            if (doCteBias(&wf3,&subab)){
                freeSingleGroup(&subab);
                return(status);
            }

            /*reset the array after bias subtraction*/
            Sub2Full(&wf3, &subab, &ab, 0, 1, 1);
        }

    } else {
        /* Full frame image, just read in the groups
           and init the mask to use all pixels
        */

        initSingleGroup (&cd);
        getSingleGroup (wf3.input, 1, &cd);
        if (hstio_err()){
            return (status = OPEN_FAILED);
        }

        initSingleGroup (&ab);
        getSingleGroup (wf3.input, 2, &ab);
        if (hstio_err()){
            return (status = OPEN_FAILED);
        }

        /*setup the mask*/
        for(i=0; i< ab.dq.data.nx; i++){
            for(j=0; j< ab.dq.data.ny; j++){
                PPix(&ab.dq.data, i, j) = 1;
                PPix(&cd.dq.data, i, j) = 1;
            }
        }

        /* SAVE A COPY OF THE RAW IMAGE FOR LATER */
        makeRAZ(&cd,&ab,&raw);

        /***SUBTRACT THE CTE BIAS FROM BOTH CHIPS IN PLACE***/
        if (doCteBias(&wf3,&cd)){
            freeSingleGroup(&cd);
            return(status);
        }

        if (doCteBias(&wf3,&ab)){
            freeSingleGroup(&ab);
            return(status);
        }
        /*SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
          AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
          CTE CODE VARIABLES.
          */
        if (CompareCTEParams(&cd, &cte_pars))
            return (status);

    }

    /*CONVERT TO RAZ, SUBTRACT BIAS AND CORRECT FOR GAIN*/
    if (raw2raz(&wf3, &cd, &ab, &raz))
        return (status);

    SingleGroup fff;                                              
    initSingleGroup(&fff);
    allocSingleGroup(&fff, RAZ_COLS, RAZ_ROWS, True);

    double cte_ff;

    cte_ff=  (wf3.expstart       - cte_pars.cte_date0)/ 
             (cte_pars.cte_date1 - cte_pars.cte_date0);

    printf("CTE_FF: %8.3f \n",cte_ff);

    cte_pars.scale_frac=cte_ff;                                      

    for(i=0;i<RAZ_COLS;i++) { 
        for(j=0;j<RAZ_ROWS;j++) {  
            Pix(fff.sci.data,i,j) =  cte_ff * (j+1)/((double)XAMP_SCI_DIM);
        }
    }

    /* 
     * If the PCTERNOI value from the primary header of the science image is non-zero, it is
     * used in the CTE algorithm.  Otherwise the read noise must be computed via find_raz2rnoival.
     * FLOAT_RNOIVAL and FLOAT_BKGDVAL are designed to be for diagnostic purposes only.
     */
    float FLOAT_RNOIVAL = 0.;
    float FLOAT_BKGDVAL = 0.;
    readNoise = wf3.pcternoi;
    sprintf(MsgText, "PCTERNOI: %8.4f (source: primary header of science image)\n\n", readNoise);
    trlmessage(MsgText);
    /* Comparison should be OK - read from FITS header and no computation */
    if (readNoise == 0.0) {
        readNoise = find_raz2rnoival(raz.sci.data.data, &FLOAT_RNOIVAL, &FLOAT_BKGDVAL);
        sprintf(MsgText, "RNOIVAL: %8.4f BKGDVAL: %8.4f\n", FLOAT_RNOIVAL, FLOAT_BKGDVAL);
        trlmessage(MsgText);
        sprintf(MsgText, "PCTERNOI: %8.4f (source: computed on-the-fly from science image)", readNoise);
        trlmessage(MsgText);
        sprintf(MsgText, "This computed value supersedes any value obtained from the primary\nheader of the science image.\n\n");
        trlmessage(MsgText);
    }

    /* The PCTERNOI value actually used is written to the PCTERNOI keyword in
     * the output image header when it is updated below for a final time.
     */
      
    /* Invoke the updated CTE correction which does the read noise 
       mitigation in each of the three forward-model iterations.
    */
    trlmessage("CTE: jumping into the routine...");
    ret = sub_ctecor_v2c(raz.sci.data.data,
                         fff.sci.data.data,
                         WsMAX,
                         cte_pars.qlevq_data,
                         cte_pars.dpdew_data,
                         cte_pars.rprof->data.data,
                         cte_pars.cprof->data.data,
                         readNoise,
                         cte_pars.thresh,
                         cte_pars.n_forward,
                         cte_pars.n_par,
                         rzc.sci.data.data);
    trlmessage("CTE: returning from the routine...");

    for (i=0;i<RAZ_COLS;i++){
        for (j=0;j<RAZ_ROWS;j++){
             Pix(chg.sci.data,i,j) = (Pix(rzc.sci.data,i,j) - Pix(raz.sci.data,i,j))/wf3.ccdgain; 
             Pix(rzc.sci.data,i,j) =  Pix(raw.sci.data,i,j) + Pix(chg.sci.data,i,j);
        }
    }

    freeSingleGroup(&fff);

    /*BACK TO NORMAL FORMATTING*/
    /*Copies rzc data to cd->sci.data and ab->sci.data */
    undoRAZ(&cd,&ab,&rzc);

    /* COPY BACK THE SCIENCE SUBARRAYS AND
       SAVE THE NEW RAW FILE WITH UPDATED SCIENCE
       ARRAYS AND PRIMARY HEADER TO RAC
       */
    if (wf3.subarray) {
        if (wf3.chip == 2) {
            /*** SAVE USEFUL HEADER INFORMATION ***/
            if (cteHistory (&wf3, subcd.globalhdr))
                return (status);

            /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
            PutKeyDbl(subcd.globalhdr, "PCTEFRAC", cte_pars.scale_frac,"CTE scaling fraction based on expstart");
            trlmessage("PCTEFRAC saved to header");

            PutKeyFlt(subcd.globalhdr, "PCTERNOI", readNoise,"read noise amp clip limit");
            trlmessage("PCTERNOI saved to header");

            Full2Sub(&wf3, &subcd, &cd, 0, 1, 1);
            putSingleGroup(output, 1, &subcd,0);
            freeSingleGroup(&subcd);
        } else {

            /*** SAVE USEFUL HEADER INFORMATION ***/
            if (cteHistory (&wf3, subab.globalhdr))
                return (status);

            /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
            PutKeyDbl(subab.globalhdr, "PCTEFRAC", cte_pars.scale_frac,"CTE scaling fraction based on expstart");
            trlmessage("PCTEFRAC saved to header");

            PutKeyFlt(subab.globalhdr, "PCTERNOI", readNoise,"read noise amp clip limit");
            trlmessage("PCTERNOI saved to header");

            Full2Sub(&wf3, &subab, &ab, 0, 1, 1);
            putSingleGroup(output, 1, &subab,0);
            freeSingleGroup(&subab);
        }

    } else { /*FUll FRAME*/
        /*** SAVE USEFUL HEADER INFORMATION ***/
        if (cteHistory (&wf3, cd.globalhdr))
            return (status);

        /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
        PutKeyDbl(cd.globalhdr, "PCTEFRAC", cte_pars.scale_frac,"CTE scaling fraction based on expstart");
        trlmessage("PCTEFRAC saved to header");

        PutKeyFlt(cd.globalhdr, "PCTERNOI", readNoise,"read noise amp clip limit");
        trlmessage("PCTERNOI saved to header");

        putSingleGroup(output,cd.group_num, &cd,0);
        putSingleGroup(output,ab.group_num, &ab,0);
    }

    /** CLEAN UP ON AISLE 3 **/
    freeSingleGroup(&rzc);
    freeSingleGroup(&rsc);
    freeSingleGroup(&chg);
    freeSingleGroup(&raz);
    freeSingleGroup(&rsz);
    freeSingleGroup(&raw);
    freeSingleGroup(&cd);
    freeSingleGroup(&ab);

    time_spent = ((double) clock()- begin +0.0) / CLOCKS_PER_SEC;
    if (verbose){
        sprintf(MsgText,"CTE run time: %.2f(s) with %i procs/threads\n",time_spent/max_threads,max_threads);
        trlmessage(MsgText);
    }

    PrSwitch("pctecorr", COMPLETE);
    if(wf3.printtime)
        TimeStamp("PCTECORR Finished",wf3.rootname);

    return (status);
}


/********************* SUPPORTING SUBROUTINES *****************************/

int raw2raz(WF3Info *wf3, SingleGroup *cd, SingleGroup *ab, SingleGroup *raz){
    /*

       convert a raw file to raz file: CDAB longwise amps, save data array
       for comparison with what jay has during testing

       -->do an additional bias correction using the  residual bias level measured for each amplifier from the
       steadiest pixels in the horizontal overscan and subtracted fom the pixels for that amplifier.

       ---> convert into electrons at the end
       ---> add supplemental bias info to the header

       allocate contiguous 2d array on the heap
       with pointers and return the pointer to the head of the array

       The Following macros are used to represent 2-d indexing.
       Two dimensional arrays are stored in FITS order.

       ny
       ^
       N | a05   a15   a25   a35
       A | a04   a14   a24   a34
       X | a03   a13   a23   a33
       I | a02   a12   a22   a32
       S | a01   a11   a21   a31
       2 | a00   a10   a20   a30
       ---------------------------> nx
       NAXIS1

       NAXIS1 is 4 and NAXIS2 is 6
       PIX(a,1,4) accesses a14

       In the raz image, each quadrant has been rotated such that the readout amp is located at the lower left.
       The reoriented four quadrants are then arranged into a single 8412x2070 image (science pixels plus overscan),
       with amps C, D, A, and B, in that order. In the raz image, pixels are all parallel-shifted down,
       then serial-shifted to the left.

*/
    extern int status;

    int i,j,k;              /*loop counters*/
    int subcol = (RAZ_COLS/NAMPS); /* for looping over quads  */
    extern int status;      /* variable for return status */
    float bias_post[NAMPS];
    float bsig_post[NAMPS];
    float bias_pre[NAMPS];
    float bsig_pre[NAMPS];
    float gain;

    /*INIT THE ARRAYS*/
    for(i=0;i<NAMPS;i++){
        bias_post[i]=0.;
        bsig_post[i]=0.;
        bias_pre[i]=0.;
        bsig_pre[i]=0.;
    }

    gain=wf3->ccdgain;

    /*REFORMAT TO RAZ*/
    makeRAZ(cd,ab,raz);


    /*SUBTRACT THE EXTRA BIAS CALCULATED, AND MULTIPLY BY THE GAIN
      Note that for user subarray the image is in only 1 quad, and only
      has prescan bias pixels so the regions are different for full and subarrays
    */
    if (wf3->subarray){
        findPreScanBias(raz, bias_pre, bsig_pre);
        for (k=0;k<NAMPS;k++){
            for (i=0; i<subcol;i++){
                for (j=0;j<RAZ_ROWS; j++){
                    if(Pix(raz->dq.data,i+k*subcol,j)){
                        Pix(raz->sci.data,i+k*subcol,j) -= bias_pre[k];
                        Pix(raz->sci.data,i+k*subcol,j) *= gain;
                    }
                }
            }
        }
    } else {
        findPostScanBias(raz, bias_post, bsig_post);
        for (k=0;k<NAMPS;k++){
            for (i=0; i<subcol;i++){
                for (j=0;j<RAZ_ROWS; j++){
                    Pix(raz->sci.data,i+k*subcol,j) -= bias_post[k];
                    Pix(raz->sci.data,i+k*subcol,j) *= gain;
                }
            }
        }
    }

    return(status);
}

/*calculate the post scan and bias after the biac file has been subtracted
  add some history information to the header

  Jay gave no explanation why plist is limited to 55377 for full arrays, his
  subarray limitation was just 1/4 of this value.  The value 55377 is the number of 
  post-scan pixels in the physical pixel vertical extent (27 x 2051 = 55377).
  Value 2051 is the veritical number of science pixels in an amp, and 
  27 is the 30 post-scan pixels with two pixels stripped from the left boundary
  and one stripped from the right boundary.

  the serial virtual overscan pixels are also called the trailing-edge pixels
  these only exist in full frame images
  */

int findPostScanBias(SingleGroup *raz, float *mean, float *sigma){

    extern int status;
    int arrsize = 55377;
    int i,j,k;              /*Looping variables */
    float plist[arrsize];  /*bias bpixels to measure*/
    float *plistSub;
    float min=0.0;
    float max=0.0;
    float rmean=0.0;
    float rsigma=0.0;
    float sigreg =7.5; /*sigma clip*/


    int subcol = RAZ_COLS/4;
    int npix=0; /*track array size for resistant mean*/

    /*init plist for full size
      We'll allocate heap memory for smaller arrays
      */
    for (i=0;i<arrsize;i++){
        plist[i]=0.;
    }

    for (k=0;k<NAMPS;k++){  /*for each quadrant cdab = 0123*/
        npix=0; /*reset for each quad*/
        rmean=0.;
        rsigma=0.;
        for (i=RAZ_ROWS+5;i<= subcol-1; i++){ /*quad area for post scan bias pixels*/
            for (j=0; j<YAMP_SCI_DIM; j++){
                if (npix < arrsize){
                    if ( Pix(raz->dq.data,i+k*subcol,j)) {
                        plist[npix] = Pix(raz->sci.data,i+k*subcol,j);
                        npix+=1;
                    }
                }
            }
        }
        if (npix > 0 ){
            plistSub = (float *) calloc(npix, sizeof(float));
            if (plistSub == NULL){
                trlerror("out of memory for resistmean entrance in findPostScanBias.");
                free(plistSub);
                return (ERROR_RETURN);
            }
            for(i=0; i<npix; i++){
                plistSub[i]=plist[i];
            }
            resistmean(plistSub, npix, sigreg, &rmean, &rsigma, &min, &max);
            free(plistSub);
        }
        mean[k]= rmean;
        sigma[k] = rsigma;
    }
    return status;
}

/*CALCULATE THE PRE SCAN AND BIAS AFTER THE BIAC FILE HAS BEEN SUBTRACTED

  The serial physical overscan pixels are also known as the serial prescan,
  they are the only pixels available for subarrays. For full frame arrays
  the prescan is not used as part of the correction, instead the virtual
  overscan pixels are used and modeled in findPostScanBias.

*/

int findPreScanBias(SingleGroup *raz, float *mean, float *sigma){
    /** this calls resistmean, which does a better job clipping outlying pixels
      that just a standard stddev clip single pass*/

    extern int status;
    int arrsize = 55377;
    int i,j,k;              /*Looping variables */
    float plist[arrsize];    /*bias pixels to measure*/
    float *plistSub; /*heap allocation for variable size plist array*/
    float min=0.0;
    float max=0.0;
    float rmean;
    float rsigma;
    float sigreg =7.5; /*sigma clip*/
    int subcol = RAZ_COLS/4;
    int npix=0; /*track array size for resistant mean*/


    /*init plist*/
    for (i=0;i<arrsize;i++){
        plist[i]=0.;
    }

    for (k=0;k<NAMPS;k++){  /*for each quadrant, CDAB ordered*/
        npix=0;
        rmean=0.;
        rsigma=0.;
        for (i=5;i<P_OVRSCN; i++){
            for (j=0; j<YAMP_SCI_DIM; j++){ /*all rows*/
                if (npix < arrsize ){
                    if (Pix(raz->dq.data,i+(k*subcol),j)){
                        plist[npix] = Pix(raz->sci.data,i+k*subcol,j);
                        npix+=1;
                    }
                }
            }
         }

        if (0 < npix ){
            plistSub = (float *) calloc(npix, sizeof(float));
            if (plistSub == NULL){
                trlerror("out of memory for resistmean entrance in findPreScanBias.");
                free(plistSub);
                return (ERROR_RETURN);
            }
            for(i=0; i<npix; i++){
                plistSub[i]=plist[i];
            }
            resistmean(plistSub, npix, sigreg, &rmean, &rsigma, &min, &max);
            free(plistSub);
        }

        mean[k]= rmean;
        sigma[k] = rsigma;
        if(npix>0)
            printf("npix=%i\nmean[%i]=%f\nsigma[%i] = %f\n",npix,k+1,rmean,k+1,rsigma);
    }
    return status;
}

/*This is the workhorse subroutine; it simulates the readout
  of one column pixi() and outputs this to pixo() using a single
  iteration.  It can be called successively to do the transfer
  in steps.


  JDIM == RAZ_ROWS
  WDIM == TRAPS  Ws is the input traps number < WsMAX
  NITs == cte_pars->n_par

  These are already in the parameter structure CTEParams
  int     Ws              the number of traps < WsMAX
  float     q_w[TRAPS];     the run of charge with level  == qlevq_data
  float   dpde_w[TRAPS];  the run of charge loss with level == dpdew_data
  float   rprof_wt[TRAPS][100]; the emission probability as fn of downhill pixel == rprof fits image
  float   cprof_wt[TRAPS][100]; the cumulative probability cprof_t( 1)  = 1. - rprof_t(1)  == cprof fits image


  W = wcol_data = trap id

  q_w[TRAP] = qlev_q from QPROF  traps as function of packet size = cte->qlevq_data[TRAP]

  pixi (curr), pixo (read) , pixf(cteff) are passed and are 1d arrays which have values for a particular column

  the ttrap reference to the image array has to be -1 for C
  */

int sim_colreadout_l(double *pixi, double *pixo, double *pixf, CTEParams *cte){

    extern int status;
    int j;
    int ttrap;

    int w;
    double ftrap;
    double pix_1;
    double padd_2;
    double padd_3;
    double prem_3;
    double pmax;
    double fcarry;

    padd_3=0.0;
    prem_3=0.0;
    padd_2=0.0;
    fcarry=0.0;
    pix_1=0.0;
    w=0;
    j=0;
    ftrap=0.0;
    ttrap=0;

    FloatHdrData *rprof;
    FloatHdrData *cprof;

    /*from the reference table*/
    rprof = cte->rprof;
    cprof = cte->cprof;

    /*FIGURE OUT WHICH TRAPS WE DON'T NEED TO WORRY ABOUT IN THIS COLUMN
      PMAX SHOULD ALWAYS BE POSITIVE HERE  */
    pmax=10.;
    for(j=0; j<RAZ_ROWS; j++){
        pixo[j] = pixi[j];
        if (pixo[j] > pmax)
            pmax=pixo[j];
    }

    /*GO THROUGH THE TRAPS ONE AT A TIME, FROM HIGHEST TO LOWEST Q,
      AND SEE WHEN THEY GET FILLED AND EMPTIED, ADJUST THE PIXELS ACCORDINGLY*/
    for (w = cte->cte_traps-1; w>=0; w--){
        if ( cte->qlevq_data[w] <= pmax ) {

            ftrap = 0.0e0;
            ttrap = cte->cte_len; /*for referencing the image at 0*/
            fcarry = 0.0e0;

            /*GO UP THE COLUMN PIXEL BY PIXEL*/
            for(j=0; j<RAZ_ROWS;j++){
                pix_1 = pixo[j];


                if ( (ttrap < cte->cte_len) || ( pix_1 >= cte->qlevq_data[w] - 1. ) ){
                    if (pixo[j] >= 0 ){
                        pix_1 = pixo[j] + fcarry; /*shuffle charge in*/
                        fcarry = pix_1 - floor(pix_1); /*carry the charge remainder*/
                        pix_1 = floor(pix_1); /*reset pixel*/
                    }

                    /*HAPPENS AFTER FIRST PASS*/
                    /*SHUFFLE CHARGE IN*/
                    if ( j> 0  ) {
                        if (pixf[j] < pixf[j-1])
                            ftrap *= (pixf[j] /  pixf[j-1]);
                    }

                    /*RELEASE THE CHARGE*/
                    padd_2=0.0;
                    if (ttrap <cte->cte_len){
                        ttrap += 1;
                        padd_2 = Pix(rprof->data,w,ttrap-1) *ftrap;
                    }

                    padd_3 = 0.0;
                    prem_3 = 0.0;
                    if ( pix_1 >= cte->qlevq_data[w]){
                        prem_3 =  cte->dpdew_data[w] / cte->n_par * pixf[j];  /*dpdew is 1 in file */
                        if (ttrap < cte->cte_len)
                            padd_3 = Pix(cprof->data,w,ttrap-1)*ftrap;
                        ttrap=0;
                        ftrap=prem_3;
                    }

                    pixo[j] += padd_2 + padd_3 - prem_3;
                } /*replaces trap continue*/
            }/*end if j>0*/
        }/* end if qlevq > pmax, replaces continue*/

    }/*end for w*/

    return(status);

}


int initCTETrl (char *input, char *output) {

    extern int status;

    char trl_in[CHAR_LINE_LENGTH+1];     /* trailer filename for input */
    char trl_out[CHAR_LINE_LENGTH+1];    /* output trailer filename */
    int exist;


    int MkName (char *, char *, char *, char *, char *, int);
    int TrlExists (char *);

    /* Initialize internal variables */
    trl_in[0] = '\0';
    trl_out[0] = '\0';
    exist = EXISTS_UNKNOWN;

    /* Input and output suffixes. */
    char *isuffix[] = {"_raw"};
    char *osuffix[] = {"_rac_tmp"};
    char *trlsuffix[] = {""};

    int nsuffix = 1;


    /* Start by stripping off suffix from input/output filenames */
    if (MkOutName (input, isuffix, trlsuffix, nsuffix, trl_in, CHAR_LINE_LENGTH)) {
        WhichError (status);
        sprintf (MsgText, "Couldn't determine trailer filename for %s",
                input);
        trlmessage (MsgText);
    }
    if (MkOutName (output, osuffix, trlsuffix, nsuffix, trl_out, CHAR_LINE_LENGTH)) {
        WhichError (status);
        sprintf (MsgText, "Couldn't create trailer filename for %s",
                output);
        trlmessage (MsgText);
    }

    /* NOW, CONVERT TRAILER FILENAME EXTENSIONS FROM '.FITS' TO '.TRL' */

    if (MkNewExtn (trl_in, TRL_EXTN) ) {
        sprintf (MsgText, "Error with input trailer filename %s", trl_in);
        trlerror (MsgText);
        WhichError (status);
    }
    if (MkNewExtn (trl_out, TRL_EXTN) ) {
        sprintf (MsgText, "Error with output trailer filename %s", trl_out);
        trlerror (MsgText);
        WhichError (status);
    }

    /* If we are working with a RAW file, then see if a TRL file
       needs to be overwritten after the generic conversion comments.  */
    if (strstr(input, isuffix[0]) != NULL) {
        /* Test whether the output file already exists */
        exist = TrlExists(trl_out);
        if (exist == EXISTS_YES) {
            /* The output file exists, so we want to add to them
             ** the new trailer comments.  */
            SetTrlOverwriteMode (NO);
        }
    }

    /* Sets up temp trailer file for output and copies input
     ** trailer file into it.  */
    InitTrlFile (trl_in, trl_out);

    return(status);
}


/* 

   #2 int sim_colreadout_l_uvis_w    --- CTE correction for one column
   #3 int sub_ctecor_v2c             --- reverse CTE correction for image

*/ 




/* ------------------------------------------- */
/*                                             */
/*  Readnoise correction for a single column.  */
/*                                             */
/* ------------------------------------------- */


int rm_rnZ_colj(double *pixj_chg, 
                double *pixj_rnz, 
                double *pixj_rsz, 
                double RNMIT)      {

       int    NIT;
       int    j;

       double   dd;
       double   dtot;
       int      ntot;
       double   ftot;
       double   rtot;
       double   RNN;
       double   RNU;
       double   sqrt();

       int      pixj_ffu[RAZ_ROWS];

       for(j=0;j<RAZ_ROWS;j++) {  
           pixj_rsz[j] = pixj_chg[j];
           pixj_rnz[j] = 0.0;
           pixj_ffu[j] = 0.0;
       } 

       /*
        * Find the upper and lower limits where there are valid
        * pixels - this is done to accommodate subarrays
        */
       int j1 = XAMP_SCI_DIM-1;
       int j2 = 2;

       /* There are no "greater than zero" pixels below j1 */
       for (j=2; j<=XAMP_SCI_DIM-1; j++) {
           if (j1==XAMP_SCI_DIM-1 && pixj_chg[j] > 0) 
               j1 = j;
       }
       /* There are no "greater than zero" pixels above j2*/
       for (j=XAMP_SCI_DIM-1; j>=2; j--) {
           if (j2==2 && pixj_chg[j] > 0) 
               j2 = j;
       }

       /* 
        * For each interation, allow a bit more noise.  This way we can
        * stop when just enough is allowed to go through the column from 2nd
        * to 2nd from the top.
        */
       for(NIT=1;NIT<NITMAX;NIT++) {  
           RNN = RNMIT*(1.00+5.0*NIT/(float)NITMAX);

           /* Bounds of this loop include only the non-zero portion of the data */
           for(j=j1; j<=j2; j++) {

              /* Compare each pixel to the average of its up/down neighbors */
              dd =  pixj_rsz[j]-(pixj_rsz[j+1]+pixj_rsz[j-1])/2.0;
              pixj_ffu[j] = 0.0;

              /* If the pixel is within the readnoise window... */
              if (fabs(dd) < RNN) { 
                  /* Cap the adjustments we are willing to make at any time to 20% */
                  if (dd >  RNMIT/5.0) dd =  RNMIT/5.0;
                  if (dd < -RNMIT/5.0) dd = -RNMIT/5.0;
                  /* Take half of the maximum adjustment from the current pixel... */
                  pixj_rnz[j  ] = pixj_rnz[j  ] + dd*0.50;
                  /* ...and give half of the value to the pixel below and above */
                  pixj_rnz[j-1] = pixj_rnz[j-1] - dd*0.25;
                  pixj_rnz[j+1] = pixj_rnz[j+1] - dd*0.25;
                  /* Flag this pixel to be in the readnoise range so we can track
                   * the total noise in the nominal pixels 
                   */
                  pixj_ffu[j  ] = 1.0; 
              }
          }
          /* Bounds of this loop include only the non-zero portion of the data */
          for(j=j1; j<=j2; j++) { 
              pixj_rsz[j] = pixj_chg[j] - pixj_rnz[j]; 
          }
          dtot = 0.;
          ntot = 0.;
          ftot = 0.;
          rtot = 0.;
          /* Bounds of this loop include only the non-zero portion of the data */
          for(j=j1; j<=j2; j++) { 
             ftot = ftot + pixj_ffu[j];
             dtot = dtot + pixj_rnz[j]*pixj_rnz[j];
             dd =  pixj_rsz[j]-(pixj_rsz[j+1]+pixj_rsz[j-1])/2.0;
             rtot = rtot + dd*dd;
             ntot = ntot + 1; 
          }
          RNU = sqrt(dtot/ftot);

          if (RNU > RNMIT) return(0);
       }

       return(0); 
}



/* ---------------------------------------------------------------*/
/*                                                                */
/*   CTE correction for a single column.                          */ 
/*                                                                */
/*   TDIM is the length of the trails that are being considered.  */
/*                                                                */
/* ---------------------------------------------------------------*/

#define _TDIM_ 60
 
int sim_colreadout_l_uvis_w(double *pixi,     // input column array (JDIM)
                            double *pixo,     // outout column array (JDIM)
                            double *pixf,     // scaling of model for each pixel (JDIM)
                            int    J1,        // bottom and top pixel in column
                            int    J2,        // bottom and top pixel in column
                            int    JDIM,      // number of pixels in column
                            double *q_w,      // the charge level for trap#W (_WDIM_)
                            double *dpde_w,   // the amt of charge this trap grabs (_WDIM_)
                            int    NITs,      // the num of iterations (dpde-->dpde/NITs)
                            float *rprof_wt,  // the trap emission for t=T (_WDIM_,100)
                            float *cprof_wt,  // the amount left in trap after t=T emission (_WDIM_,100)
                            int    Ws)  {

      int    j;       // pixel location up the column

      double ftrap;   // total number of electrons in the trap
      int    ttrap;   // shifts since the trap was last filled

      int    w;       // trap counter

      double pmax;    // max pixel value in the column - tells us the highest relevant trap numbrer
      int    Wf;      // highest relevant trap number

      double prel_1;  // amount in incidental release from trap
      double prel_2;  // amount in flush release of trap
      double pgrb_3;  // amount grabbed by filling trap

      float rprof_t[_TDIM_];              
      float cprof_t[_TDIM_];              

      /* Bounds checking */
      if (Ws>WsMAX) {
          printf("Ws error\n");
          return(1); 
      }

      /* Figure out which traps we do not need 
         to worry about in this column
      */
      pmax = 10;
      for(j=0;j<JDIM;j++) {
         pixo[j] = pixi[j];
         if (pixo[j] > pmax) pmax = pixo[j];
      }

      /* Figure out the highest trap number we need to consider */
      Wf = 1;
      for (w=0;w<Ws;w++) {
         if (pmax >=q_w[w]) Wf = w;
      }

      /* Go thru the traps one at a time (from highest to lowest q)
         and see when they get filled and emptied; adjust the
         pixel values accordingly
      */
      for (w=Wf;w>=0;w--) {   // loop backwards

         for(ttrap=0;ttrap<_TDIM_;ttrap++) { 
             rprof_t[ttrap] = rprof_wt[w+ttrap*WsMAX];
             cprof_t[ttrap] = cprof_wt[w+ttrap*WsMAX];
         }
                                                       
         /* Initialize the flux in the trap to zero */
         ftrap =   0.0;

         /* Initialize the time-since-flush to the max */
         ttrap =  _TDIM_ + 1;

         /* Go up the column, pixel-by-pixel */
         for (j=J1;j<J2;j++) {
            
            /* If we have an inversion of the density (i.e., a readout-cosmic issue),
               then we do not want to flush too much
            */
            if (j>J1) {
                if (pixf[j] < pixf[j-1]) {
                    ftrap = pixf[j]/
                            pixf[j-1]*ftrap;
                }
            }

            /* Set up accounting of pixel value changes */
            prel_1 = 0.;   // charge to be released
            prel_2 = 0.;   // charge to be flushed out
            pgrb_3 = 0.;   // charge to be trapped

            /* Filled/refilled trap#W */
            if (pixo[j] >= q_w[w]) {

                /* Need to flush before filling? */
                if (ttrap < (_TDIM_)) {

                    /* Increment time since filled */
                    ttrap = ttrap + 1;

                    /* Flush out amount for this shift, and ...*/
                    prel_1 = rprof_t[ttrap-1]*ftrap;

                    /* ...flush out the rest */
                    prel_2 = cprof_t[ttrap-1]*ftrap;
                }

                /* Amount to hold in the trap */
                ftrap  = dpde_w[w]/NITs*pixf[j];
 
                /* Subtract that amount held from the pixel and reset 
                   the time-since-filled counter
                */
                pgrb_3 = ftrap;
                ttrap = 0;
            }

            /* trap#W not filled at this time */
            else {

                /* Check if the trap contains charge, and if so, then
                   release the appropriate number of electrons
                */
                if (ttrap < (_TDIM_)) {
                    ttrap = ttrap + 1;
                    prel_1 = rprof_t[ttrap-1]*ftrap;
                }
            }

            /* Make adjustments to the output pixel: add the trail emission, 
               flush the trap, and fill the trap
            */
            pixo[j] = pixo[j] + prel_1 + prel_2 - pgrb_3;
         }
      }
      return(0);
     }


/* --------------------------------- */
/*                                   */ 
/*   CTE correction for one column.  */
/*                                   */ 
/* --------------------------------- */


int sub_ctecor_v2c(float *pixz_raz, 
                   float *pixz_fff,
                   int Ws,
                   double *q_w,
                   double *dpde_w,
                   float *rprof_wt, 
                   float *cprof_wt, 
                   float PCTERNOI,
                   float FIX_ROCR,
                   int PCTENFOR,
                   int PCTENPAR,
                   float *pixz_rzc) 
      { 

      extern int status;

      int i;
      int j;
      int jj;
      int jmax;

      int    NITFOR, NITFORs;
      int    NITPAR, NITPARs;
      double RNOI;
      int ret;

      double *pixj_fff;
      double *pixj_raz;
      double *pixj_mod;
      double *pixj_rnz;
      double *pixj_rsz;
      double *pixj_org;
      double *pixj_obs;
      double *pixj_chg;

      int NCRX;
      int DONE;

      int NDONE = 0;

      RNOI    = PCTERNOI;
      NITFORs = PCTENFOR;
      NITPARs = PCTENPAR;

      printf("                             \n");
      printf("   INSIDE sub_ctecor_v2.f... \n");
      printf("          ---> PCTERNOI: %8.4f \n",PCTERNOI);
      printf("          ---> FIX_ROCR: %8.4f \n",FIX_ROCR);
      printf("          --->  NITFORs: %5d \n",NITFORs);
      printf("          --->  NITPARs: %5d \n",NITPARs);
      printf("                             \n");

      #pragma omp parallel \
       shared(pixz_raz,pixz_fff,pixz_rzc,              \
              NITPARs,NITFORs,                         \
              q_w,dpde_w,                              \
              rprof_wt,cprof_wt,Ws,NDONE)              \
       private(i,j,ret,NCRX, DONE, NITFOR,NITPAR,      \
               pixj_fff, pixj_raz, pixj_mod, pixj_rnz, \
               pixj_rsz, pixj_org, pixj_obs, pixj_chg)

      #pragma omp for


      for(i=0;i<RAZ_COLS;i++) {  

         pixj_fff = malloc(RAZ_ROWS*8);
         pixj_raz = malloc(RAZ_ROWS*8);
         pixj_mod = malloc(RAZ_ROWS*8);
         pixj_rnz = malloc(RAZ_ROWS*8);
         pixj_rsz = malloc(RAZ_ROWS*8);
         pixj_org = malloc(RAZ_ROWS*8);
         pixj_obs = malloc(RAZ_ROWS*8);
         pixj_chg = malloc(RAZ_ROWS*8);

         for(j=0;j<RAZ_ROWS;j++) { 
            pixj_raz[j] = pixz_raz[i+j*RAZ_COLS]; 
            pixj_fff[j] = pixz_fff[i+j*RAZ_COLS]; 
         }

         NCRX = 0;
         DONE = 0;
         while(!DONE) { 
             NCRX = NCRX + 1;
             DONE = 1;
             for (j=0;j<RAZ_ROWS;j++) { 
                pixj_mod[j] = pixj_raz[j];           
                pixj_chg[j] = 0.0;
             }
             for(NITFOR=0;NITFOR<NITFORs;NITFOR++) {
                 ret = rm_rnZ_colj(pixj_mod,pixj_rnz,pixj_rsz,RNOI);
                 for(j=0;j<RAZ_ROWS;j++) { 
                    pixj_org[j] = pixj_rsz[j]; 
                 }
                 for(NITPAR=1;NITPAR<=NITPARs;NITPAR++) { 
                     ret = sim_colreadout_l_uvis_w(pixj_org,
                                                   pixj_obs,
                                                   pixj_fff,
                                                   1,RAZ_ROWS,RAZ_ROWS,
                                                   q_w,dpde_w,NITPARs,
                                                   rprof_wt,cprof_wt,Ws);
                     for (j=0;j<RAZ_ROWS;j++) { 
                        pixj_org[j] = pixj_obs[j];
                     }
                 }
                 for(j=0;j<RAZ_ROWS;j++) { 
                     pixj_chg[j] = pixj_obs[j] - pixj_rsz[j];
                     pixj_mod[j] = pixj_raz[j] - pixj_chg[j];
                 }
             }
             if (FIX_ROCR<0) { 
                 for(j=15;j<=2060;j++) { 
                     if (pixj_mod[j] < FIX_ROCR && 
                         pixj_mod[j]-pixj_raz[j] < FIX_ROCR && 
                         pixj_mod[j] < pixj_mod[j+1] && 
                         pixj_mod[j] < pixj_mod[j-1]) {
                         jmax = j;
                         for(jj=j-2;jj<j;jj++) {    
                             if (pixj_mod[jj  ]-pixj_raz[jj  ] >    
                                 pixj_mod[jmax]-pixj_raz[jmax]) {
                                 jmax = jj;
                             }
                         }
                         if (pixj_mod[jmax]-pixj_raz[jmax] > -2.5*FIX_ROCR) {
                             pixj_fff[jmax] = pixj_fff[jmax]*0.90;        
                             DONE = NCRX >= 10;
                         }
                     }
                 }
             }
         }

         for(j=0;j<RAZ_ROWS;j++) { 
             pixz_rzc[i+j*RAZ_COLS] = pixj_mod[j];
             pixz_fff[i+j*RAZ_COLS] = pixj_fff[j];
         }

         free(pixj_fff);
         free(pixj_raz);
         free(pixj_mod);
         free(pixj_rnz);
         free(pixj_rsz);
         free(pixj_org);
         free(pixj_obs);
         free(pixj_chg);

         /* This variable exists for debuggin purposes. */
         NDONE++;
         /*if (NDONE==(NDONE/100)*100) { printf("  i = %5d   %5d  %5d \n",i,NCRX,NDONE); }*/

      }

      return(status);
}

/*
 * This routine dynamically determines the noise in the input image.
 */
float find_raz2rnoival(float *raz_cdab, float *FLOAT_RNOIVAL, float *FLOAT_BKGDVAL) { 

      /* Return value */
      float RNOIVALo;

      int i, j, ik;

      float     b, d;
      int       ih, iih;
      long      dhist[NUM_BINS], dcum[NUM_BINS], vtot;
      long      vhist[NUM_BINS], vcum[NUM_BINS], dtot;
      int       ivmin, id1, id2;
      int       idmin, iv1, iv2;
      long long vsum;
      long long nsum;

      int       j1, j2;

      float RNOIVAL;
      float RNOIVALu;

      float BKGDVAL;
      float BKGDVALu;

      FLOAT_RNOIVAL[0] = 3.33;
      FLOAT_BKGDVAL[0] = raz_cdab[0];

      iv1 = 1;
      iv2 = 999;
      id1 = 1;
      id2 = 999;

      /*
       * Distill the image variation information and background into quick histograms 
       */
      for (ih=1; ih<=NUM_BINS; ih++) {
          dhist[ih-1] = 0;
          vhist[ih-1] = 0;
      }

      for (i=2; i<=RAZ_COLS; i++) { 
          /*
           * Find the upper and lower limits where there are valid pixels
           * to accommodate subarrays
           */
          j1 = XAMP_SCI_DIM-1;
          j2 = 2;

          /* There are no "greater than zero" pixels below j1 */
          for (j=2; j<=XAMP_SCI_DIM-1; j++) { 
             if (j1==XAMP_SCI_DIM-1 && raz_cdab[i+(j-1)*RAZ_COLS-1] > 0) 
                 j1 = j;
          } 

          /* There are no "greater than zero" pixels above j2 */
          for (j=XAMP_SCI_DIM-1; j>=2; j--) { 
             if (j2==2 && raz_cdab[i+(j-1)*RAZ_COLS-1] > 0) 
                 j2 = j;
          } 

          /*
           * Process the valid pixels. "ik" is the chip horizontal locator.
           * For each pixel in the recorded part of the detector, find the
           * average of the surrounding 8 pixels
           */
          for (j=j1; j<=j2; j++) {
              ik = i - (i-1)/2103*2103;
              if (ik > 25 && ik < XAMP_SCI_DIM+24) {
                  b = (raz_cdab[i-1+(j+1-1)*RAZ_COLS-1]
                      +raz_cdab[i  +(j+1-1)*RAZ_COLS-1]
                      +raz_cdab[i+1+(j+1-1)*RAZ_COLS-1]
                      +raz_cdab[i-1+(j  -1)*RAZ_COLS-1] 
                      +raz_cdab[i+1+(j  -1)*RAZ_COLS-1] 
                      +raz_cdab[i-1+(j-1-1)*RAZ_COLS-1]
                      +raz_cdab[i  +(j-1-1)*RAZ_COLS-1]
                      +raz_cdab[i+1+(j-1-1)*RAZ_COLS-1])/8.00;

                 /* Local residual as proxy for noise */
                 d = raz_cdab[i+(j-1)*RAZ_COLS-1] - b;

                 /* Locate within the histogram bin; 4.5 is to spread the values out.
                  * The value of 4.5 is 3 times the gain, where the gain is 1.5
                  */
                 ih = 501 + d*(SPREAD_FOR_HISTO) + 0.5;
                 if (ih < 1)
                     ih = 1;                                   
                 if (ih > NUM_BINS)
                     ih = NUM_BINS;

                 /* Increment the histogram bin */
                 dhist[ih-1] = dhist[ih-1] + 1;
  
                 /* Locate the pixel value within the histogram bin */
                 ih = 501 + raz_cdab[i+(j-1)*RAZ_COLS-1]*(SPREAD_FOR_HISTO) + 0.5;

                 if (ih < 1)
                     ih = 1;
                 if (ih > NUM_BINS)
                     ih = NUM_BINS;
                
                 /* Increment the histogram bin */
                 vhist[ih-1] = vhist[ih-1] + 1;
              }
          }
      }

      /* Compute the cumulative distribution for the noise and for the background */
      dtot = 0;
      vtot = 0;
      for (ih=1; ih<=NUM_BINS; ih++) {                  
         if (ih > 1 && ih < NUM_BINS) {
             dtot = dtot + dhist[ih-1];
             vtot = vtot + vhist[ih-1];
         }
         dcum[ih-1] = dtot;
         vcum[ih-1] = vtot;
      }

      idmin = 999;
      ivmin = 999;
    
      /*
       * Find the closest 75% of the points and use them to determine the noise
       * and the background
       */
      for (ih=1; ih<=NUM_BINS-1; ih++) {
          for (iih=ih+1; iih<=NUM_BINS-1; iih++) {
              if (dcum[iih-1]-dcum[ih-1] > 0.75*dtot && iih-ih < idmin) {
                  id1   = ih;
                  id2   = iih;
                  idmin = iih-ih;   
              }

              if (vcum[iih-1]-vcum[ih-1] > 0.75*vtot && iih-ih < ivmin) { 
                  iv1   = ih;
                  iv2   = iih;
                  ivmin = iih-ih;
              }      
          }
      }

      nsum = 0;
      vsum = 0;

      for (ih=iv1; ih<=iv2; ih++) { 
          nsum = nsum + vhist[ih-1];
          vsum = vsum + vhist[ih-1]*(ih-501);
      }

      if (vsum==0 || nsum==0) {
         RNOIVALu = 9.75 ;
         BKGDVALu = 999.9 ;
         *FLOAT_RNOIVAL = RNOIVALu;
         *FLOAT_BKGDVAL = BKGDVALu;
         return(RNOIVALu);
      }

      /* For debugging purposes only
      printf("  \n");
      printf("   vsum: %12lld  \n",vsum);
      printf("   nsum: %12lld  \n",nsum);
      printf("  \n");
      printf("   dbar: %12.2f   \n",idmin/2.30*(SPREAD_FOR_HISTO));
      printf("   vbar: %12.2f %12lld %12lld \n",vsum/nsum*(SPREAD_FOR_HISTO),vsum,nsum);
      printf("  \n");
      */

      RNOIVAL  = (int)(idmin/2.30/(SPREAD_FOR_HISTO)/sqrt(1+1/8.0)*4+0.5)/4.00;
      RNOIVAL  =       idmin/2.30/(SPREAD_FOR_HISTO)/sqrt(1+1/8.0);
      RNOIVALu = RNOIVAL;
      if (RNOIVALu > HIGH_CLIP) 
          RNOIVALu = HIGH_CLIP;

      BKGDVAL  = 1.*vsum/nsum/(SPREAD_FOR_HISTO);
      BKGDVALu = BKGDVAL;
      if (BKGDVALu > 999.9) 
          BKGDVALu = 999.9;

      /* Values which can be used for diagnostic analysis */
      *FLOAT_RNOIVAL = RNOIVALu;
      *FLOAT_BKGDVAL = BKGDVALu;

      /* LOW_CLIP and HIGH_CLIP are imposed limits on the computed noise value */
      RNOIVALo = RNOIVALu;
      if (RNOIVALo < LOW_CLIP) 
          RNOIVALo = LOW_CLIP;
      if (RNOIVALo > HIGH_CLIP) 
          RNOIVALo = HIGH_CLIP;

      return(RNOIVALo);
}
