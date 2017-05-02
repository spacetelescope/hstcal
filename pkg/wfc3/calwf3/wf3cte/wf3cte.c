/* WFC3 -- CTE loss correction for UVIS

   M. sosey  Aug-2014  Adapted for the pipeline from Jay Andersons CTE correction code for wfc3 UVIS
   raw2raz_wfc3uv.F , an edited file was delivered december 2014, and both are different from the
   fortran code currently served on the wfc3 website.

   M. Sosey Aug-2016 Adapted to be used with Subarrays as well as full frame arrays,
   as long as the subarray contains physical overscan pixels, which don't include the science team subarrays
   which can span quads.
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

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "wf3corr.h"
# include "cte.h"


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
    PCTERNOI - readnoise amplitude for clipping
    PCTESMIT - number of iterations used in CTE forward modeling
    PCTESHFT - number of iterations used in the parallel transfer
    PCTENSMD - readnoise mitigation algorithm
    PCTETRSH - over-subtraction threshold
    PCTEFRAC - cte scaling frac calculated from expstart
    PCTERNOI - the readnoise clipping level to use

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
    Hdr phdr; /*primary header for input image, all output information saved here*/
    Hdr scihdr; /*science header in case of subarray image to detect chip*/
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
            if ( start >= 25 &&  finish + 60 <= (RAZ_COLS/2) - 25){
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
            if ( start >= 25 &&  finish + 60 <= (RAZ_COLS/2) - 25){
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

    /***CALCULATE THE SMOOTH READNOISE IMAGE***/
    trlmessage("CTE: Calculating smooth readnoise image");


    /***CREATE THE NOISE MITIGATION MODEL ***/
    if (cte_pars.noise_mit == 0) {
        if (raz2rsz(&wf3, &raz, &rsz, cte_pars.rn_amp, max_threads))
            return (status);
    } else {
        trlmessage("Only noise model 0 implemented!");
        return (status=ERROR_RETURN);
    }

    /***CONVERT THE READNOISE SMOOTHED IMAGE TO RSC IMAGE
      THIS IS WHERE THE CTE GETS CALCULATED         ***/
    if (rsz2rsc(&wf3, &rsz, &rsc, &cte_pars))
        return (status);

    /*** CREATE THE FINAL CTE CORRECTED IMAGE, PUT IT BACK INTO ORIGNAL RAW FORMAT***/
    for (i=0;i<RAZ_COLS;i++){
        for(j=0; j<RAZ_ROWS; j++){
           Pix(chg.sci.data,i,j) = (Pix(rsc.sci.data,i,j) - Pix(rsz.sci.data,i,j))/wf3.ccdgain;
           Pix(rzc.sci.data,i,j) =  Pix(raw.sci.data,i,j) + Pix(chg.sci.data,i,j);
        }
    }

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
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    extern int status;      /* variable for return status */
    float bias_post[4];
    float bsig_post[4];
    float bias_pre[4];
    float bsig_pre[4];
    float gain;

    /*INIT THE ARRAYS*/
    for(i=0;i<4;i++){
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
        for (k=0;k<4;k++){
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
        for (k=0;k<4;k++){
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
  subarray limitation was just 1/4 of this value

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

    for (k=0;k<4;k++){  /*for each quadrant cdab = 0123*/
        npix=0; /*reset for each quad*/
        rmean=0.;
        rsigma=0.;
        for (i=RAZ_ROWS+5;i<= subcol-1; i++){ /*quad area for post scan bias pixels*/
            for (j=0; j<2051; j++){
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

    for (k=0;k<4;k++){  /*for each quadrant, CDAB ordered*/
        npix=0;
        rmean=0.;
        rsigma=0.;
        for (i=5;i<25; i++){
            for (j=0; j<2051; j++){ /*all rows*/
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


int raz2rsz(WF3Info *wf3, SingleGroup *raz, SingleGroup *rsz, double rnsig, int max_threads){
    /*
       This routine will read in a RAZ image and will output the smoothest
       image that is consistent with being the observed image plus readnoise. (RSZ image)
       This is necessary because we want the CTE-correction algorithm to produce the smoothest
       possible reconstruction, consistent with the original image and the
       known readnoise.  This algorithm constructs a model that is smooth
       where the pixel-to-pixel variations can be thought of as being related
       to readnoise, but if the variations are too large, then it respects
       the pixel values.  Basically... it uses a 2-sigma threshold.

       This is strategy #1 in a two-pronged strategy to mitigate the readnoise
       amplification.  Strategy #2 will be to not iterate when the deblurring
       is less than the readnoise.

*/

    extern int status;

    int i, j, NIT; /*loop variables*/
    int imid;
    double dptr=0.0;
    double  rms=0.0;
    double  rmsu=0.0;
    double nrms=0.0;
    double nrmsu=0.0;
    float hardset=0.0f;
    double setdbl=0.0;


    /*1D ARRAYS FOR CENTRAL AND NEIGHBORING RAZ_COLS*/
    double obs_loc[3][RAZ_ROWS] ;
    double rsz_loc[3][RAZ_ROWS] ;

    NIT=1;

    /*ALL ELEMENTS TO FLAG*/
    for(i=0;i<3;i++){
        for (j=0; j<RAZ_ROWS; j++){
            obs_loc[i][j]=setdbl;
            rsz_loc[i][j]=setdbl;
        }
    }

    /***INITIALIZE THE LOCAL IMAGE GROUPS***/
    SingleGroup rnz;
    initSingleGroup(&rnz);
    allocSingleGroup(&rnz, RAZ_COLS, RAZ_ROWS, True);

    SingleGroup zadj;
    initSingleGroup(&zadj);
    allocSingleGroup(&zadj, RAZ_COLS, RAZ_ROWS, True);


    /*COPY THE RAZ IMAGE INTO THE RSZ OUTPUT IMAGE
      AND INITIALIZE THE OTHER IMAGES*/
    for(i=0;i<RAZ_COLS;i++){
        for (j=0;j<RAZ_ROWS;j++){
            Pix(rsz->sci.data,i,j) = Pix(raz->sci.data,i,j);
            Pix(rsz->dq.data,i,j) = Pix(raz->dq.data,i,j);
            Pix(rnz.sci.data,i,j) = hardset;
            Pix(zadj.sci.data,i,j) = hardset;
        }
    }


    /*THE RSZ IMAGE JUST GETS UPDATED AS THE RAZ IMAGE IN THIS CASE*/
    if (rnsig < 0.1){
        trlmessage("rnsig < 0.1, No read-noise mitigation needed");
        return(status);
    }

    /*GO THROUGH THE ENTIRE IMAGE AND ADJUST PIXELS TO MAKE THEM
      SMOOTHER, BUT NOT SO MUCH THAT IT IS NOT CONSISTENT WITH
      READNOISE.  DO THIS IN BABY STEPS SO THAT EACH ITERATION
      DOES VERY LITTLE ADJUSTMENT AND INFORMATION CAN GET PROPAGATED
      DOWN THE LINE.
      */

    rms=setdbl;

    for(NIT=1; NIT<=100; NIT++){
        #pragma omp parallel for schedule(dynamic) \
        private(i,j,imid,obs_loc,rsz_loc,dptr)\
        shared(raz, rsz, rnsig,rms,nrms, zadj)
        for(i=0; i<RAZ_COLS; i++){
            imid=i;
            /*RESET TO MIDDLE RAZ_COLS AT ENDPOINTS*/
            if (imid < 1)
                imid=1;
            if (imid == RAZ_COLS-1)
                imid = RAZ_COLS-2;

            /*COPY THE MIDDLE AND NEIGHBORING PIXELS FOR ANALYSIS*/
            for(j=0; j<RAZ_ROWS; j++){
                obs_loc[0][j] = Pix(raz->sci.data,imid-1,j);
                obs_loc[1][j] = Pix(raz->sci.data,imid,j);
                obs_loc[2][j] = Pix(raz->sci.data,imid+1,j);

                rsz_loc[0][j] = Pix(rsz->sci.data,imid-1,j);
                rsz_loc[1][j] = Pix(rsz->sci.data,imid,j);
                rsz_loc[2][j] = Pix(rsz->sci.data,imid+1,j);
            }
            for (j=0; j<RAZ_ROWS; j++){
             if(Pix(raz->dq.data,imid,j)) {
                find_dadj(1+i-imid,j, obs_loc, rsz_loc, rnsig, &dptr);
                Pix(zadj.sci.data,i,j) = dptr;
              }
            }
        } /*end the parallel for*/

        /*NOW GO OVER ALL THE RAZ_COLS AND RAZ_ROWS AGAIN TO SCALE THE PIXELS
        */
        for(i=0; i<RAZ_COLS;i++){
            for(j=0; j<RAZ_ROWS; j++){
                if (Pix(raz->dq.data,i,j)){
                    Pix(rsz->sci.data,i,j) +=  (Pix(zadj.sci.data,i,j)*0.75);
                    Pix(rnz.sci.data,i,j) = (Pix(raz->sci.data,i,j) - Pix(rsz->sci.data,i,j));
                }
            }
        }

        rms=setdbl;
        nrms=setdbl;

        /*This is probably a time sink because the arrays are being
          accessed out of storage order, careful of page faults */
        #pragma omp parallel for schedule(dynamic,1)\
        private(i,j,rmsu,nrmsu) \
        shared(raz,rsz,rms,rnsig,nrms)
        for(j=0; j<RAZ_ROWS; j++){
            nrmsu=setdbl;
            rmsu=setdbl;
            for(i = 0;i<RAZ_COLS; i++){
                if ( (fabs(Pix(raz->sci.data,i,j)) > 0.1) ||
                        (fabs(Pix(rsz->sci.data,i,j)) > 0.1) ){
                    rmsu  +=  ( Pix(rnz.sci.data,i,j) * Pix(rnz.sci.data,i,j) );
                    nrmsu += 1.0;
                }
            }
            #pragma omp critical (rms)
            {   rms  += rmsu;
                nrms += nrmsu;
            }
        }
        rms = sqrt(rms/nrms);

        /*epsilon type comparison*/
        if ( (rnsig-rms) < 0.00001) break; /*this exits the NIT for loop*/
    } /*end NIT*/

    freeSingleGroup(&zadj);
    freeSingleGroup(&rnz);

    return (status);
}


int find_dadj(int i ,int j, double obsloc[][RAZ_ROWS], double rszloc[][RAZ_ROWS], double rnsig, double *d){
    /*
       This function determines for a given pixel how it can
       adjust in a way that is not inconsistent with its being
       readnoise.  To do this, it looks at its upper and lower
       neighbors and sees whether it is consistent with either
       (modulo readnoise).  To the extent that it is consistent
       then move it towards them.  But also bear in mind that
       that we don't want it to be more than 2 RN sigmas away
       from its original value.  This is pretty much a tug of
       war... with readnoise considerations pushing pixels to
       be closer to their neighbors, but the original pixel
       values also pull to keep the pixel where it was.  Some
       accommodation is made for both considerations.
       */

    extern int status;

    double mval=0.0;
    double    dval0, dval0u, w0;
    double    dval9, dval9u, w9;
    double    dmod1, dmod1u, w1;
    double    dmod2, dmod2u, w2;

    dval0=0.;
    dval0u=0.;
    w0=0.;
    dval9=0.;
    dval9u=0.;
    w9=0.;
    dmod1=0.;
    dmod1u=0.;
    w1=0.;
    dmod2=0.;
    dmod2u=0.;
    w2=0.;

    mval = rszloc[i][j];
    dval0  = obsloc[i][j] - mval;
    dval0u = dval0;

    if (dval0u >1.0)
        dval0u =  1.0;
    if (dval0u <-1.0)
        dval0u = -1.0;

    dval9 = 0.;

    /*COMPARE THE SURROUNDING PIXELS*/
    if (i==1 &&  RAZ_ROWS-1>=j  && j>0 ) {

        dval9 = obsloc[i][j-1]  - rszloc[i][j-1] +
            obsloc[i][j]    - rszloc[i][j]  +
            obsloc[i][j+1]  - rszloc[i][j+1] +
            obsloc[i-1][j-1]- rszloc[i-1][j-1] +
            obsloc[i-1][j]  - rszloc[i-1][j] +
            obsloc[i-1][j+1]- rszloc[i-1][j+1] +
            obsloc[i+1][j-1]- rszloc[i+1][j-1] +
            obsloc[i+1][j]  - rszloc[i+1][j] +
            obsloc[i+1][j+1]- rszloc[i+1][j+1];
    }

    dval9 =dval9 / 9.;
    dval9u = dval9;

    if (dval9u > (rnsig*0.33))
        dval9u =  rnsig*0.33;
    if (dval9u <  rnsig*-0.33)
        dval9u = rnsig*-0.33;

    dmod1 = 0.;
    if (j>0)
        dmod1 = rszloc[i][j-1] - mval;

    dmod1u = dmod1;
    if (dmod1u > rnsig*0.33)
        dmod1u =  rnsig*0.33;
    if (dmod1u < rnsig*-0.33)
        dmod1u = rnsig*-0.33;

    dmod2 = 0.;
    if (j < RAZ_ROWS-1)
        dmod2 =  rszloc[i][j+1] - mval;

    dmod2u = dmod2;
    if (dmod2u > rnsig*0.33)
        dmod2u =  rnsig*0.33;
    if (dmod2u < rnsig*-0.33)
        dmod2u = rnsig*-0.33;


    /*
       IF IT'S WITHIN 2 SIGMA OF THE READNOISE, THEN
       TEND TO TREAT AS READNOISE; IF IT'S FARTHER OFF
       THAN THAT, THEN DOWNWEIGHT THE INFLUENCE
       */
    w0 =   (dval0*dval0) / ((dval0*dval0)+ 4.0*(rnsig*rnsig));
    w9 =   (dval9*dval9) / ((dval9*dval9)+ 18.0*(rnsig*rnsig));
    w1 = (4*rnsig*rnsig) / ((dmod1*dmod1)+4.0*(rnsig*rnsig));
    w2 = (4*rnsig*rnsig) / ((dmod2*dmod2)+4.0*(rnsig*rnsig));

    /*(note that with the last two, if a pixel
      is too discordant with its upper or lower
      that neighbor has less of an ability to
      pull it)*/

    *d = ((dval0u * w0 * 0.25f) + /* desire to keep the original pixel value */
            (dval9u*w9*0.25f) + /* desire to keep the original sum over 3x3*/
            (dmod1u*w1*0.25f) + /*desire to get closer to the pixel below*/
            (dmod2u*w2*0.25f)) ; /*desire to get closer to the pixel above*/

    return(status);
}


/*** THIS ROUTINE PERFORMS THE CTE CORRECTIONS
  rsz is the readnoise smoothed image
  rsc is the coorection output image
  rac = raw + ((rsc-rsz) / gain )

 ***/
int rsz2rsc(WF3Info *wf3, SingleGroup *rsz, SingleGroup *rsc, CTEParams *cte) {

    extern int status;

    int i,j;
    double cte_i=0.0;
    double cte_j=0.0;
    double ro=0;
    int io=0;
    double ff_by_col[RAZ_COLS][4];
    float hardset=0.0;

    /*These are already in the parameter structure
      int     Ws              the number of traps < 999999, taken from pctetab read
      int     q_w[TRAPS];     the run of charge with level  cte->qlevq_data[]
      float   dpde_w[TRAPS];  the run of charge loss with level cte->dpdew_data[]

      float   rprof_wt[TRAPS][100]; the emission probability as fn of downhill pixel, TRAPS=999
      float   cprof_wt[TRAPS][100]; the cumulative probability cprof_t( 1)  = 1. - rprof_t(1)

      The rprof array gives the fraction of charge that comes out of every parallel serial-shift
      the cummulative distribution in cprof then tells you what's left

*/

    SingleGroup pixz_fff;
    initSingleGroup(&pixz_fff);
    allocSingleGroup(&pixz_fff, RAZ_COLS, RAZ_ROWS, True);

    /*SCALE BY 1 UNLESS THE PCTETAB SAYS OTHERWISE, I IS THE PACKET NUM
      THIS IS A SAFETY LOOP INCASE NOT ALL THE COLUMNS ARE POPULATED
      IN THE REFERENCE FILE*/

    for(i=0; i<RAZ_COLS;i++){
        ff_by_col[i][0]=1.;
        ff_by_col[i][1]=1.;
        ff_by_col[i][2]=1.;
        ff_by_col[i][3]=1.;
        j= cte->iz_data[i]; /*which column to scale*/
        ff_by_col[j][0]=cte->scale512[i];
        ff_by_col[j][1]=cte->scale1024[i];
        ff_by_col[j][2]=cte->scale1536[i];
        ff_by_col[j][3]=cte->scale2048[i];

        /*CALCULATE THE CTE CORRECTION FOR EVERY PIXEL
          Index is figured on the final size of the image
          not the current size. Moved above
          */

        for(j=0; j<RAZ_ROWS; j++){
            Pix(pixz_fff.sci.data,i,j)=hardset;
            ro = j/512.0; /*ro can be zero, it's an index*/
            if (ro <0 ) ro=0.;
            if (ro > 2.999) ro=2.999; /*only 4 quads, 0 to 3*/
            io = (int) floor(ro); /*force truncation towards 0 for pos numbers*/
            cte_j= (j+1) / 2048.0;
            cte_i= ff_by_col[i][io] + (ff_by_col[i][io+1] -ff_by_col[i][io]) * (ro-io);
            Pix(pixz_fff.sci.data,i,j) =  (cte_i*cte_j);
        }
    }

    /*FOR REFERENCE TO JAYS CODE, FF_BY_COL IS WHAT'S IN THE SCALE BY COLUMN

      int   iz_data[RAZ_ROWS];  column number in raz format
      double scale512[RAZ_ROWS];      scaling appropriate at row 512
      double scale1024[RAZ_ROWS];     scaling appropriate at row 1024
      double scale1536[RAZ_ROWS];     scaling appropriate at row 1536
      double scale2048[RAZ_ROWS];     scaling appropriate at row 2048
      */

    /*THIS IS RAZ2RAC_PAR IN JAYS CODE - MAIN CORRECTION LOOP IN HERE*/
    inverse_cte_blur(rsz, rsc, &pixz_fff, cte, wf3->verbose,wf3->expstart);
    freeSingleGroup(&pixz_fff);
    return(status);
}



/*** this routine does the inverse CTE blurring... it takes an observed
  image and generates the image that would be pushed through the readout
  algorithm to generate the observation

  CTE_FF is found using the observation date of the data
  FIX_ROCRs is cte->fix_rocr
  Ws is the number of TRAPS that are < 999999

  this is sub_wfc3uv_raz2rac_par in jays code

  floor rounds to negative infinity
  ceiling rounds to positive infinity
  truncate rounds up or down to zero
  round goes to the nearest integer

  fff is the input cte scaling array calculated over all pixels
  This is a big old time sink function
 ***/

int inverse_cte_blur(SingleGroup *rsz, SingleGroup *rsc, SingleGroup *fff, CTEParams *cte, int verbose, double expstart){

    extern int status;

    /*looping vars*/
    int NREDO, REDO;
    int NITINV, NITCTE;
    int i;
    int j,jj;
    double dmod;
    int jmax;
    float hardset=0.0f;
    int totflux=0;

    double cte_ff; /*cte scaling based on observation date*/
    double setdbl=0.0;

    /*DEFINE TO MAKE PRIVATE IN PARALLEL RUN*/
    double *pix_obsd=&setdbl;
    double *pix_modl=&setdbl;
    double *pix_curr=&setdbl;
    double *pix_init=&setdbl;
    double *pix_read=&setdbl;
    double *pix_ctef=&setdbl;

    /*STARTING DEFAULTS*/
    NITINV=1;
    NITCTE=1;
    cte_ff=0.0;
    jmax=0;
    dmod=0.0;

    /*LOCAL IMAGES TO PLAY WITH, THEY WILL REPLACE THE INPUTS*/
    SingleGroup rz; /*pixz_raz*/
    initSingleGroup(&rz);
    allocSingleGroup(&rz, RAZ_COLS, RAZ_ROWS, True);

    SingleGroup rc; /*pixz_rac*/
    initSingleGroup(&rc);
    allocSingleGroup(&rc, RAZ_COLS, RAZ_ROWS, True);

    SingleGroup pixz_fff; /*pixz_fff*/
    initSingleGroup(&pixz_fff);
    allocSingleGroup(&pixz_fff, RAZ_COLS, RAZ_ROWS, True);


    /*USE EXPSTART YYYY-MM-DD TO DETERMINE THE CTE SCALING
      APPROPRIATE FOR THE GIVEN DATE. WFC3/UVIS WAS
      INSTALLED AROUND MAY 11,2009 AND THE MODEL WAS
      CONSTRUCTED TO BE VALID AROUND SEP 3, 2012, A LITTLE
      OVER 3 YEARS AFTER INSTALLATION*/

    cte_ff=  (expstart - cte->cte_date0)/ (cte->cte_date1 - cte->cte_date0);
    cte->scale_frac=cte_ff;   /*save to param structure for header update*/

    if(verbose){
        sprintf(MsgText,"CTE_FF (scaling fraction by date) = %g",cte_ff);
        trlmessage(MsgText);
    }

    /*SET UP THE SCALING ARRAY WITH INPUT DATA, hardset arrays for safety*/
    for (i=0;i<RAZ_COLS;i++){
        for(j=0;j<RAZ_ROWS;j++){
            Pix(rc.sci.data,i,j)=hardset;
            Pix(rz.sci.data,i,j)=hardset;
            Pix(pixz_fff.sci.data,i,j)=hardset;
            Pix(rz.sci.data,i,j) = Pix(rsz->sci.data,i,j);
            Pix(rz.dq.data,i,j) = Pix(rsz->dq.data,i,j);
            Pix(pixz_fff.sci.data,i,j) =  cte_ff * Pix(fff->sci.data,i,j);
        }
    }

    #pragma omp parallel for schedule (dynamic,1) \
    private(dmod,i,j,jj,jmax,REDO,NREDO,totflux, \
            pix_obsd,pix_modl,pix_curr,pix_init,\
            pix_read,pix_ctef,NITINV,NITCTE)\
    shared(rc,rz,cte,pixz_fff)

    for (i=0; i< RAZ_COLS; i++){
        pix_obsd = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_modl = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_curr = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_init = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_read = (double *) calloc(RAZ_ROWS, sizeof(double));
        pix_ctef = (double *) calloc(RAZ_ROWS, sizeof(double));

        totflux=0;
        /*HORIZONTAL PRE/POST SCAN POPULATION */
        for (j=0; j< RAZ_ROWS; j++){
            if(Pix(rz.dq.data,i,j)){
                pix_obsd[j] = Pix(rz.sci.data,i,j); /*starts as input RAZ*/
                totflux += 1;
            }
        }

        if (totflux >= 1) {/*make sure the column has flux in it*/
            NREDO=0; /*START OUT NOT NEEDING TO MITIGATE CRS*/
            REDO=0; /*FALSE*/
            do { /*replacing goto 9999*/
                /*STARTING WITH THE OBSERVED IMAGE AS MODEL, ADOPT THE SCALING FOR THIS COLUMN*/
                for (j=0; j<RAZ_ROWS; j++){
                    pix_modl[j] =  Pix(rz.sci.data,i,j);
                    pix_ctef[j] =  Pix(pixz_fff.sci.data,i,j);
                }
                /*START WITH THE INPUT ARRAY BEING THE LAST OUTPUT
                  IF WE'VE CR-RESCALED, THEN IMPLEMENT CTEF*/
                for (NITINV=1; NITINV<=cte->n_forward; NITINV++){
                    for (j=0; j<RAZ_ROWS; j++){
                        pix_curr[j]=pix_modl[j];
                        pix_read[j]=pix_modl[j];
                        pix_ctef[j]=Pix(pixz_fff.sci.data,i,j);
                    }

                    /*TAKE EACH PIXEL DOWN THE DETECTOR IN NCTENPAR=7*/
                    for (NITCTE=1; NITCTE<=cte->n_par; NITCTE++){
                        sim_colreadout_l(pix_curr, pix_read, pix_ctef, cte);

                        /*COPY THE JUST UPDATED READ OUT IMAGE INTO THE INPUT IMAGE*/
                        for (j=0; j< RAZ_ROWS; j++){
                            pix_curr[j]=pix_read[j];
                        }
                    } /* end NITCTE */

                    /*DAMPEN THE ADJUSTMENT IF IT IS CLOSE TO THE READNOISE, THIS IS
                      AN ADDITIONAL AID IN MITIGATING THE IMPACT OF READNOISE*/
                    for (j=0; j< RAZ_ROWS; j++){
                        dmod =  (pix_obsd[j] - pix_read[j]);
                        if (NITINV < cte->n_forward){
                            dmod *= (dmod*dmod) /((dmod*dmod) + (cte->rn_amp * cte->rn_amp));
                        }
                        pix_modl[j] += dmod; /*dampen each pixel as the best is determined*/
                    }
                } /*NITINV end*/

                /*LOOK FOR AND DOWNSCALE THE CTE MODEL IF WE FIND
                  THE TELL-TALE SIGN OF READOUT CRS BEING OVERSUBTRACTED;
                  IF WE FIND ANY THEN GO BACK UP AND RERUN THIS COLUMN


                  THE WFC3 UVIS MODEL SEARCHES FOR OVERSUBTRACTED TRAILS.
                  WHICH ARE  DEFINED AS EITHER:

                  - A SINGLE PIXEL VALUE BELOW -10E-
                  - TWO CONSECUTIVE PIXELS TOTALING -12 E-
                  - THREE TOTALLING -15 E-

                  WHEN WE DETECT SUCH AN OVER-SUBTRACTED TAIL, WE ITERATIVELY REDUCE
                  THE LOCAL CTE SCALING BY 25% UNTIL THE TRAIL IS
                  NO LONGER NEGATIVE  THIS DOES NOT IDENTIFY ALL READOUT-CRS, BUT IT DOES
                  DEAL WITH MANY OF THEM. FOR IMAGES THAT HAVE BACKGROUND GREATER THAN 10 OR SO,
                  THIS WILL STILL END UP OVERSUBTRACTING CRS A BIT, SINCE WE ALLOW
                  THEIR TRAILS TO BE SUBTRACTED DOWN TO -10 RATHER THAN 0.

*/
                if (cte->fix_rocr) {
                    for (j=10; j< RAZ_ROWS-2; j++){
                        if (  (( cte->thresh > pix_modl[j] ) &&
                                    ( cte->thresh > (pix_modl[j] - pix_obsd[j]))) ||

                                (((pix_modl[j] + pix_modl[j+1]) < -12.) &&
                                 (pix_modl[j] + pix_modl[j+1] - pix_obsd[j] - pix_obsd[j+1] < -12.)) ||

                                (((pix_modl[j] + pix_modl[j+1] + pix_modl[j+2]) < -15.) &&
                                 ((pix_modl[j] + pix_modl[j+1] + pix_modl[j+2] -pix_obsd[j] -
                                   pix_obsd[j+1] - pix_obsd[j+2]) <-15.))  ){

                            jmax=j;

                            /*GO DOWNSTREAM AND LOOK FOR THE OFFENDING CR*/
                            for (jj=j-10; jj<=j;jj++){
                                if ( (pix_modl[jj] - pix_obsd[jj]) >
                                        (pix_modl[jmax] - pix_obsd[jmax]) ) {
                                    jmax=jj;
                                }
                            }
                            /* DOWNGRADE THE CR'S SCALING AND ALSO FOR THOSE
                               BETWEEN THE OVERSUBTRACTED PIXEL AND IT*/
                            for (jj=jmax; jj<=j;jj++){
                                Pix(pixz_fff.sci.data,i,jj) *= 0.75;
                            }
                            REDO=1; /*TRUE*/
                        } /*end if*/
                    } /*end for  j*/
                }/*end fix cr*/

                if (REDO) NREDO +=1;
                if (NREDO == 5)  REDO=0; /*stop*/
            } while (REDO); /*replacing goto 9999*/
        } /*totflux > 1, catch for subarrays*/

#pragma omp critical (cte)
        for (j=0; j< RAZ_ROWS; j++){
            if (Pix(rz.dq.data,i,j)){
                Pix(rc.sci.data,i,j)= pix_modl[j];
            }
        }

        free(pix_obsd);
        free(pix_modl);
        free(pix_curr);
        free(pix_init);
        free(pix_read);
        free(pix_ctef);

    } /*end i*/

    for (i=0; i< RAZ_COLS; i++){
        for (j=0; j< RAZ_ROWS; j++){
            if(Pix(rsz->dq.data,i,j)){
                Pix(rsz->sci.data,i,j) = Pix(rz.sci.data,i,j);
                Pix(rsc->sci.data,i,j) = Pix(rc.sci.data,i,j);
                Pix(fff->sci.data,i,j) = Pix(pixz_fff.sci.data,i,j);
            }
        }
    }

    freeSingleGroup(&rz);
    freeSingleGroup(&rc);
    freeSingleGroup(&pixz_fff);

    return(status);
}


/*This is the workhorse subroutine; it simulates the readout
  of one column pixi() and outputs this to pixo() using a single
  iteration.  It can be called successively to do the transfer
  in steps.


  JDIM == RAZ_ROWS
  WDIM == TRAPS  Ws is the input traps number < 999999
  NITs == cte_pars->n_par

  These are already in the parameter structure CTEParams
  int     Ws              the number of traps < 999999
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

    char trl_in[SZ_LINE+1];     /* trailer filename for input */
    char trl_out[SZ_LINE+1];    /* output trailer filename */
    int exist;


    int MkName (char *, char *, char *, char *, char *, int);
    int TrlExists (char *);
    void SetTrlOverwriteMode (int);

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
    if (MkOutName (input, isuffix, trlsuffix, nsuffix, trl_in, SZ_LINE)) {
        WhichError (status);
        sprintf (MsgText, "Couldn't determine trailer filename for %s",
                input);
        trlmessage (MsgText);
    }
    if (MkOutName (output, osuffix, trlsuffix, nsuffix, trl_out, SZ_LINE)) {
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
