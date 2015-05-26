/* WFC3 -- CTE loss correction for UVIS 

   M. sosey  Aug-2014  Adapted for the pipeline from Jay Andersons CTE correction code for wfc3 UVIS
   raw2raz_wfc3uv.F. 

*/

# include <time.h>
# include <string.h>
# include <math.h>
# include <stdlib.h>
# include <stdio.h>

# ifdef _OPENMP
#  include <omp.h>
# endif

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "wf3err.h"
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
CTEDATE0 - date of wfc3/uvis installation in HST, in fractional years
CTEDATE1 - reference date of CTE model pinning, in fractional years

PCTETDIM - max length of CTE trail
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
        
    CTEParams cte_pars; /*STRUCTURE HOLDING THE MODEL PARAMETERS*/
    SingleGroup cd; /*SCI 1*/
    SingleGroup ab; /*SCI 2*/
    SingleGroup raz; /* THE LARGE FORMAT COMBINATION OF CDAB*/
    SingleGroup rsz; /* LARGE FORMAT READNOISE CORRECTED IMAGE */
    SingleGroup rsc; /* CTE CORRECTED*/
    SingleGroup rzc; /*part of the conversion to rac*/
    SingleGroup chg; /*part of the conversion to rac*/
    
    int i,j; /*loop var*/
    int max_threads;
    
    Bool subarray; /* to verify that no subarray is being used, it's not implemented yet*/

    /*CONTAIN PARALLEL PROCESSING TO A SINGLE THREAD*/
#   ifdef _OPENMP
    trlmessage("\nUsing parallel processing provided by OpenMP inside CTE routine\n");
    if (onecpu){
        omp_set_dynamic(0);
        max_threads=1;
        sprintf(MsgText,"onecpu == TRUE, Using only %i threads/cpu", max_threads);
    } else {
        omp_set_dynamic(0);
        max_threads = omp_get_num_procs(); /*be nice, use 1 less than avail?*/
        sprintf(MsgText,"Setting max threads to %i of %i cpus\n",max_threads, omp_get_num_procs()); 
    }
    omp_set_num_threads(max_threads);
    trlmessage(MsgText);
#   endif

    /*JUST FOR VERIFICIATION IMAGES*/
    Hdr junkhdr;
    initHdr(&junkhdr);
    IODescPtr out;
    char tmpout[SZ_LINE+1]; 
   

    /* DETERMINE THE NAMES OF THE TRAILER FILES BASED ON THE INPUT
       AND OUTPUT FILE NAMES, THEN INITIALIZE THE TRAILER FILE BUFFER
       WITH THOSE NAMES.
       */
    if (initCTETrl (input, output))
        return (status);

    /* COPY COMMAND-LINE ARGUMENTS INTO WF3. */ 
    WF3Init (&wf3);
    strcpy (wf3.input, input);
    strcpy (wf3.output, output);

    PrBegin ("WFC3CTE");
    if (wf3.printtime)
        TimeStamp("WFC3CTE Started: ",wf3.rootname);

    wf3.pctecorr = cte_sw->pctecorr;
    wf3.darkcorr = cte_sw->darkcorr;
    wf3.biascorr = cte_sw->biascorr;
    wf3.printtime = printtime;
    wf3.verbose = verbose;
    wf3.refnames = refnames;

    PrFileName ("input", wf3.input);
    PrFileName ("output", wf3.output);

    /* CHECK WHETHER THE OUTPUT FILE ALREADY EXISTS. */
    if (FileExists (wf3.output))
        return (status);

    /* OPEN INPUT IMAGE IN ORDER TO READ ITS PRIMARY HEADER. */
    if (LoadHdr (wf3.input, &phdr) )
        return (status);

    /* GET KEYWORD VALUES FROM PRIMARY HEADER. */
    if (GetKeys (&wf3, &phdr)) {
        freeHdr (&phdr);
        return (status);
    }

    if (GetCTEFlags (&wf3, &phdr)) {
        freeHdr(&phdr);
        return (status);
    }

    /*READ IN THE CTE PARAMETER TABLE*/
    initCTEParams(&cte_pars);
    if (GetCTEPars (wf3.pctetab.name,&cte_pars))
        return (status);

    if (verbose){
        PrRefInfo ("pctetab", wf3.pctetab.name, wf3.pctetab.pedigree,
                wf3.pctetab.descrip, wf3.pctetab.descrip2);
    }
    

    /* OPEN THE INPUT IMAGES AND GET THE  SCIENCE EXTENSIONS  */
    initSingleGroup (&cd);
    initSingleGroup (&ab);

    getSingleGroup (wf3.input, 1, &cd);
    if (hstio_err())
        return (status = OPEN_FAILED);

    getSingleGroup (wf3.input, 2, &ab);
    if (hstio_err())
        return (status = OPEN_FAILED);


    /*** MAKE SURE THIS IS NOT A SUBARRAY ***/
    if (GetKeyBool (cd.globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
        return (status);

    if (subarray) {
        sprintf(MsgText,"SUBARRAY= %i; **SUBARRAY images are not yet supported for CTE**\n",wf3.subarray);
        trlmessage(MsgText);
        status=ERROR_RETURN;
        return(status);
    }

    if (GetKeyBool (ab.globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
        return (status);

    if (subarray) {
        sprintf(MsgText,"SUBARRAY= %i; **SUBARRAY images are not yet supported for CTE**\n",wf3.subarray);
        trlmessage(MsgText);
        status=ERROR_RETURN;
        return(status);
    }

    /*Save the PCTETABLE information to the header of the science image
      after checking to see if the user has specified any changes to the
      CTE code variables.
    */
    if (CompareCTEParams(&cd, &cte_pars)){
        return (status);
    }

    /***SUBTRACT THE BIAS FROM BOTH CHIPS IN PLACE***/
    if (doCteBias(&wf3,&cd)){
        freeSingleGroup(&cd);
        return(status);
    }

    /*SAVE THE OUTPUT BIAS SUBTRACTED IMAGE FOR TEST VERIFICATION*/
    if (verbose) {
        strcpy(tmpout,wf3.rootname);
        strcat(tmpout,"_bsub_cd.fits");
        out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
        putHeader(out);
        putFloatData(out,&cd.sci.data);
        closeImage(out);
        sprintf(MsgText,"Bias subtracted image for amps cd written to %s\n",tmpout);
        trlmessage(MsgText);
    }
    
    if (doCteBias(&wf3,&ab)){
        freeSingleGroup(&ab);
        return(status);
    }

    /*SAVE THE OUTPUT BIAS SUBTRACTED IMAGE FOR TEST VERIFICATION*/
    if (verbose) {
        strcpy(tmpout,wf3.rootname);
        strcat(tmpout,"_bsub_ab.fits");
        out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
        putHeader(out);
        putFloatData(out,&ab.sci.data);
        closeImage(out);
        sprintf(MsgText,"Bias subtracted image for amps ab written to %s\n",tmpout);
        trlmessage(MsgText);
    }
    

    
    /*CONVERT TO RAZ FORMAT AND CORRECT FOR GAIN*/
    initSingleGroup(&raz);
    allocSingleGroup(&raz, COLUMNS, ROWS);

    if (raw2raz(&wf3, &cd, &ab, &raz))
        return (status);

    /*SAVE CORRECTED RAZ AS TEMPORARY FILE FOR USER TO EXAMINE */
    if (verbose) {
        strcpy(tmpout,wf3.rootname);
        strcat(tmpout,"_raz_gaincorr.fits");
        out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
        putHeader(out);
        putFloatData(out,&raz.sci.data);
        closeImage(out);
        sprintf(MsgText,"RAZ format gain (%f) corrected image to check: %s\n",wf3.ccdgain,tmpout);
        trlmessage(MsgText);
    }    

    /***CALCULATE THE SMOOTH READNOISE IMAGE***/
    if (verbose)
        trlmessage("CTE: Calculating smooth readnoise image\n");

    initSingleGroup(&rsz);
    allocSingleGroup(&rsz, COLUMNS, ROWS);
    
    /***CHECK THE NOISE MITIGATION MODEL  cte_pars.rn_amp ***/
    if (cte_pars.noise_mit == 0) {
        if (raz2rsz(&wf3, &raz, &rsz, cte_pars.rn_amp, max_threads))
            return (status);

        if (verbose) {
            strcpy(tmpout,wf3.rootname);
            strcat(tmpout,"_rsz.fits");
            out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
            putHeader(out);
            putFloatData(out,&rsz.sci.data);
            closeImage(out);
            sprintf(MsgText,"RSZ (input for rsc) format image to check: %s\n",tmpout);
            trlmessage(MsgText);
        }
    } else {
        trlmessage("Only noise model 0 implemented!");
        return (status=ERROR_RETURN);
    }


    /***CONVERT the READNOISE SMNOOTHED IMAGE TO RSC IMAGE***/
    if (verbose)
        trlmessage("CTE: Converting RSZ to RSC\n");

    initSingleGroup(&rsc);
    allocSingleGroup(&rsc, COLUMNS, ROWS);

    /* uncomment to remove the rsz2rsc stage to test pipeline for accuracy
    for (i=0;i<COLUMNS;i++){
        for(j=0;j<ROWS;j++){
            memcpy(&Pix(rsc.sci.data,i,j), &Pix(rsz.sci.data,i,j), sizeof(float));
        }
    }*/


    if (rsz2rsc(&wf3, &rsz, &rsc, &cte_pars))
        return (status);
    
    
    if (verbose) {
        strcpy(tmpout,wf3.rootname);
        strcat(tmpout,"_rsc.fits");
        out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
        putHeader(out);
        putFloatData(out,&rsc.sci.data);
        closeImage(out);
        sprintf(MsgText,"CTE: Saved RSC (cte corrected) image to check: %s\n",tmpout);
        trlmessage(MsgText);
    }
    
    /*** CREATE THE FINAL CTE CORRECTED IMAGE, PUT IT BACK INTO ORIGNAL RAW FORMAT***/


    /*** SAVE USEFULL HEADER INFORMATION ***/
    if (cteHistory (&wf3, cd.globalhdr))
        return (status);

    initSingleGroup(&rzc);
    allocSingleGroup(&rzc, COLUMNS, ROWS);

    /*THIS IS THE CHANGE IMAGE FOR VALIDATION*/
    initSingleGroup(&chg);
    allocSingleGroup(&chg, COLUMNS, ROWS);
    
    for (i=0;i<COLUMNS;i++){
        for(j=0; j<ROWS; j++){
            Pix(rzc.sci.data,i,j) = Pix(raz.sci.data,i,j) + (Pix(rsc.sci.data,i,j) - Pix(rsz.sci.data,i,j))/wf3.ccdgain;
            Pix(chg.sci.data,i,j) = Pix(rsc.sci.data,i,j) - Pix(rsz.sci.data,i,j);
        }
    }
    
    /*WRITE OUT THE CHG  and RZC IMAGES FOR CHECKING*/
    if (verbose){
        strcpy(tmpout,wf3.rootname);
        strcat(tmpout,"_chg.fits");
        out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
        putHeader(out);
        putFloatData(out,&chg.sci.data);
        closeImage(out);
        sprintf(MsgText,"CTE: Saved CHANGE image: %s\n",tmpout);
        trlmessage(MsgText);
     
        /*WRITE OUT THE RZC IMAGE FOR CHECKING*/
        strcpy(tmpout,wf3.rootname);
        strcat(tmpout,"_rzc.fits");
        out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
        putHeader(out);
        putFloatData(out,&rzc.sci.data);
        closeImage(out);
        sprintf(MsgText,"CTE: Saved RZC image to check: %s\n",tmpout);
        trlmessage(MsgText);
    }
    
    undosciRAZ(&cd,&ab,&rzc);
    
    /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
    PutKeyFlt (cd.globalhdr, "PCTEFRAC", cte_pars.scale_frac,"cte scaling fraction based on expstart");
    trlmessage("PCTEFRAC saved to header");
    
    /*SAVE THE NEW RAW FILE WITH UPDATED SCIENCE ARRAYS AND PRIMARY HEADER*/
    sprintf(MsgText,"Writing cd[sci,%i] to %s",cd.group_num,output);
    trlmessage(MsgText);
    putSingleGroup(output,cd.group_num, &cd,0);
    

    sprintf(MsgText,"Writing ab[sci,%i] to %s",ab.group_num,output);
    trlmessage(MsgText);
    putSingleGroup(output,ab.group_num, &ab,0);

    /** clean up **/    
    freeSingleGroup(&rzc);
    freeSingleGroup(&rsc);
    freeSingleGroup(&chg);

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
        with pointers and returns the pointer to the head of the array


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
    int subcol = (COLUMNS/4); /* for looping over quads  */
    extern int status;      /* variable for return status */
    float bias_post[4];
    float bsig_post[4];
    float bias_pre[4];
    float bsig_pre[4];
    float gain;
    
    if (wf3->verbose)
        trlmessage("CTE: Converting RAW to RAZ format\n");

    gain=wf3->ccdgain;
    
    /*REFORMAT TO RAZ*/
    makesciRAZ(cd,ab,raz);

    /*SAVE A COPY OF THE RAZ FILE TO DISK FOR VALIDATING*/
    if (wf3->verbose){
        IODescPtr out;
        char razout[SZ_LINE+1];
        Hdr junkhdr;
        initHdr (&junkhdr);
        strcpy(razout,wf3->rootname);
        strcat(razout,"_raz_nogain.fits");
        out = openOutputImage(razout,"",0,&junkhdr,0,0,FITSBYTE);
        putHeader(out);
        putFloatData(out,&raz->sci.data);
        closeImage(out);
        sprintf(MsgText,"RAZ format conversion image to check: %s\n",razout);
        trlmessage(MsgText);
        freeHdr(&junkhdr);
    }
    
    findPostScanBias(raz, bias_post, bsig_post);

    /*SUBTRACT THE BIAS CALCULATED ABOVE, SPANNING THE QUADS WITH K*/
    for (k=1;k<=4;k++){
        for (i=0; i<subcol;i++){
            for (j=0;j<ROWS; j++){
                Pix(raz->sci.data,i+(k-1)*subcol,j) -= bias_post[k];
            }
        }   
    }
    /*done just for reporting*/
    findPreScanBias(raz, bias_pre, bsig_pre);

    /* convert to electrons, multiply by the gain */
    for (i=0; i<COLUMNS;i++){
        for (j=0;j<ROWS; j++){
            Pix(raz->sci.data,i,j) =  Pix(raz->sci.data,i,j) * gain;
        }
    }

    return(status);
}


/*calculate the post scan and bias after the biac file has been subtracted
  add some history information to the header
  
Jay gave no explanation why plist is limited to 55377  
  
*/

int findPostScanBias(SingleGroup *raz, float *mean, float *sigma){

    extern int status;
    
    int i,j,k;              /*Looping variables */
    float plist[55377];  /*bias bpixels to measure*/
    float min=0;
    float max=0;
    float rmean;
    float rsigma;
    float sigreg =7.5; /*sigma clip*/
    int nuse;
    

    int subcol = COLUMNS/4;
    int npix; /*track array size for resistant mean*/
    rmean=0.;
    rsigma=0.;
    npix=0;
    nuse=0;

    /*init plist*/
    for (i=0;i<55377;i++){
        plist[i]=0.;
    }
    
    trlmessage("\nPost scan bias measures:\n");
    
    for (k=1;k<=4;k++){  /*for each quadrant cdab = 0123*/
        npix=0;
        for (i=ROWS+5;i<= subcol-1; i++){ /*specific area for post scan bias pixels*/
            for (j=0; j<2051; j++){
                if (npix < 55377){
                    memcpy(&plist[npix], &Pix(raz->sci.data,i+(k-1)*subcol,j),sizeof(float));
                    npix+=1;
                }
            }
        }
                
        resistmean(plist, npix, sigreg, &rmean, &rsigma, &min, &max);
        mean[k]= rmean;
        sigma[k] = rsigma;
        sprintf(MsgText,"mean=%f\tsigma=%f",rmean,rsigma);
        trlmessage(MsgText);
        
    }
    return status;
}

/*CALCULATE THE PRE SCAN AND BIAS AFTER THE BIAC FILE HAS BEEN SUBTRACTED 
  ADD SOME HISTORY INFORMATION TO THE HEADER
  
  Jay gave no explanation why plist is limited to 55377  

*/

int findPreScanBias(SingleGroup *raz, float *mean, float *sigma){
    /** this calls resistmean, which does a better job clipping outlying pixels
        that just a standard stddev clip single pass*/
        
    extern int status;
    
    int i,j,k;              /*Looping variables */
    float plist[55377];    /*bias pixels to measure*/
    float min=0;
    float max=0;
    float rmean;
    float rsigma;
    float sigreg =7.5; /*sigma clip*/
    int nuse;
    
    int subcol = COLUMNS/4;
    int npix; /*track array size for resistant mean*/
    rmean=0.;
    rsigma=0.;
    npix=0;
    nuse=0;
    
    /*init plist*/
    for (i=0;i<55377;i++){
        plist[i]=0.;
    }
    
    trlmessage("\nPrescan residual bias measures:\n");

    for (k=1;k<=4;k++){  /*for each quadrance cdab*/
        npix=0;
        for (i=5;i<25; i++){ /*specific area fr post scan bias pixels*/
            for (j=0; j<2051; j++){
                if (npix < 55377 ){
                    memcpy(&plist[npix], &Pix(raz->sci.data,i+(k-1)*subcol,j),sizeof(float));
                    npix+=1;
                }
            }
        }
        
        resistmean(plist, npix, sigreg, &rmean, &rsigma, &min, &max);
        mean[k]= rmean;
        sigma[k] = rsigma;
        sprintf(MsgText,"mean=%f\tsigma=%f",rmean,rsigma);
        trlmessage(MsgText);
        
    }
    return status;
}


int raz2rsz(WF3Info *wf3, SingleGroup *raz, SingleGroup *rsz, float rnsig, int max_threads){
    /*
       This routine will read in a RAZ image and will output the smoothest
       image that is consistent with being the observed image plus readnoise.  
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
    char tmpout[SZ_LINE+1];
    

    if (wf3->verbose)
        trlmessage("Starting RAZ -> RSZ\n");
    
    int i, j, NIT; /*loop variables*/
    int imid;
    double d;
    double  rms, rmsu;
    int nrms, nrmsu;
    
    /*1D ARRAYS FOR CENTRAL AND NEIGHBORING COLUMNS*/
    float obs_loc[3][ROWS] = {0}; /*all elements to zero*/
    float rsz_loc[3][ROWS] = {0};
        
    d=0.;
    rms=0.;
    rmsu=0.;
    nrms=0;
    nrmsu=0;
    NIT=1;
        
    /***INITIALIZE THE LOCAL IMAGE GROUPS***/
    SingleGroup rnz;
    initSingleGroup(&rnz);
    allocSingleGroup(&rnz, COLUMNS, ROWS);

        
    SingleGroup zadj;
    initSingleGroup(&zadj);
    allocSingleGroup(&zadj, COLUMNS, ROWS);

    /*COPY THE RAZ IMAGE INTO THE RSZ OUTPUT IMAGE*/
    for(i=0;i<COLUMNS;i++){
        for (j=0;j<ROWS;j++){
            memcpy( &Pix(rsz->sci.data,i,j), &Pix(raz->sci.data,i,j),sizeof(float) );
        }
    }


    /*THE RSZ IMAGE JUST GETS UPDATED AS THE RAZ IMAGE IN THIS CASE*/
    if (rnsig < 0.1){
        trlmessage("rnsig < 0.1, no rnoise mitigation needed");
        return(status); 
    }

    /*GO THROUGH THE ENTIRE IMAGE AND ADJUST PIXELS TO MAKE THEM
      SMOOTHER, BUT NOT SO MUCH THAT IT IS NOT CONSISTENT WITH
      READNOISE.  DO THIS IN BABY STEPS SO THAT EACH ITERATION
      DOES VERY LITTLE ADJUSTMENT AND INFORMATION CAN GET PROPAGATED
      DOWN THE LINE.
    */
             

    rms = 0.;      
    if(wf3->verbose) {
        sprintf(MsgText,"RNSIG=%3.2f\nNIT\tCHGrms\t\ticol\tORIG\t\tDIFF\t\tRSZ\n",rnsig);
        trlmessage(MsgText);
    }

    for(NIT=1; NIT<=100; NIT++){ 

        #pragma omp parallel for schedule(dynamic) \
           private(i,j,imid,obs_loc,rsz_loc)\
           shared(raz, rsz, rnsig)
           
        for(i=0; i<COLUMNS; i++){
            imid=i;

            /*RESET TO MIDDLE COLUMNS AT ENDPOINTS*/
            if (imid == 0)
                imid=1;
            if (imid == COLUMNS-1)
                imid = COLUMNS-2;

            /*COPY THE MIDDLE AND NEIGHBORING PIXELS FOR ANALYSIS*/    
            for(j=0; j<ROWS; j++){
                memcpy(&obs_loc[0][j],&Pix(raz->sci.data,imid-1,j),sizeof(float));
                memcpy(&obs_loc[1][j],&Pix(raz->sci.data,imid,j),sizeof(float));
                memcpy(&obs_loc[2][j],&Pix(raz->sci.data,imid+1,j),sizeof(float));

                memcpy(&rsz_loc[0][j],&Pix(rsz->sci.data,imid-1,j),sizeof(float));
                memcpy(&rsz_loc[1][j],&Pix(rsz->sci.data,imid,j),sizeof(float));
                memcpy(&rsz_loc[2][j],&Pix(rsz->sci.data,imid+1,j),sizeof(float));
            }
            for (j=0; j<ROWS; j++){  
                find_dadj(1+i-imid,j, obs_loc, rsz_loc, rnsig, &d);
                Pix(zadj.sci.data,i,j)=d;
                if (wf3->verbose){
                    if (j==1999 && i==19){
                        sprintf(MsgText,"%2i\t%f\t%4i\t%f\t%f\t%f",NIT,rms,i+1,
                            fminf(Pix(raz->sci.data,i,j),999.9),
                            (Pix(raz->sci.data,i,j) - Pix(rsz->sci.data,i,j)),
                            fminf(Pix(rsz->sci.data,i,j), 999.9));
                        trlmessage(MsgText);
                    }
                }
            }
        } /*end the parallel for*/
    

        /*NOW GO OVER ALL THE COLUMNS AND ROWS AGAIN TO SCALE THE PIXELS */
        for(i=0; i<COLUMNS;i++){
            for(j=0; j<ROWS; j++){
                Pix(rsz->sci.data,i,j) += Pix(zadj.sci.data,i,j)*0.75;
                Pix(rnz.sci.data,i,j) = Pix(raz->sci.data,i,j) - Pix(rsz->sci.data,i,j);
            }                
        }

        nrms = 0;
        rms  = 0.;
       
       /*#pragma omp parllel for default(none) private(i,j,rmsu,nrmsu) shared(raz,rsz,rms)*/
       for(j=0; j<ROWS; j++){
            nrmsu = 0;
            rmsu  = 0.;
            for(i = 0;i<COLUMNS; i++){
                if ( (abs(Pix(raz->sci.data,i,j)) > 0.1) || 
                     (abs(Pix(rsz->sci.data,i,j)) > 0.1) ){

                    rmsu  +=  pow( ((double)Pix(raz->sci.data,i,j) - (double)Pix(rsz->sci.data,i,j)),2) ;
                    nrmsu += 1;
                }
            }
            #pragma omp critical 
            {rms  += rmsu;
            nrms += nrmsu;}
        }
        rms = sqrt(rms/(double)nrms);

        if ( rms > (double) rnsig) break; /*this exits the NIT for loop*/
     } /*end NIT*/


    /*write the rnz image for validation*/
    if (wf3->verbose){
        IODescPtr out;
        Hdr junkhdr;
        initHdr(&junkhdr);
        strcpy(tmpout,wf3->rootname);
        strcat(tmpout,"_rnz.fits");
        out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
        putHeader(out);
        putFloatData(out,&rnz.sci.data);
        closeImage(out);
        sprintf(MsgText,"Wrote RNZ file for verification  %s",tmpout);
        trlmessage(MsgText);
        freeHdr(&junkhdr);
    }            

    freeSingleGroup(&zadj);
    freeSingleGroup(&rnz);
 
    
    return (status);
}

int find_dadj(int i ,int j, float obsloc[][ROWS], float rszloc[][ROWS], float rnsig, double *d){
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

    float mval;
    float    dval0, dval0u, w0;
    float    dval9, dval9u, w9;
    float    dmod1, dmod1u, w1;
    float    dmod2, dmod2u, w2;

    mval = rszloc[i][j];
    dval0  = obsloc[i][j] - mval;
    dval0u = dval0;

    if (dval0u >1.0)
        dval0u =  1.0;
    if (dval0u <-1.0)
        dval0u = -1.0;

    dval9 = 0.;
    
    /*COMPARE THE SURROUNDING PIXELS*/
    if (i==1 &&  ROWS>j  && j>0 ) {
    
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

    dval9 /=9.;
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
    if (j < ROWS) 
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
    w0 =   (dval0*dval0) / ((dval0*dval0)+ 4.*(rnsig*rnsig));
    w9 =   (dval9*dval9) / ((dval9*dval9)+ 18.*(rnsig*rnsig));
    w1 = (4*rnsig*rnsig) / ((dmod1*dmod1)+4.*(rnsig*rnsig));
    w2 = (4*rnsig*rnsig) / ((dmod2*dmod2)+4.*(rnsig*rnsig));

    /*(note that with the last two, if a pixel
    is too discordant with its upper or lower
    that neighbor has less of an ability to
    pull it)*/   
   
    *d = (dval0u * w0 * 0.25) + /* desire to keep the original pixel value */
    (dval9u*w9*0.25) + /* desire to keep the original sum over 3x3*/
    (dmod1u*w1*0.25) + /*desire to get closer to the pixel below*/
    (dmod2u*w2*0.25) ; /*desire to get closer to the pixel above*/

    return(status);
}


/*** THIS ROUTINE PERFORMS THE CTE CORRECTIONS ***/
int rsz2rsc(WF3Info *wf3, SingleGroup *rsz, SingleGroup *rsc, CTEParams *cte) {

    extern int status;
    
    int i,j;
    float cte_i, cte_j;
    float ro;
    int io;
    float ff_by_col[COLUMNS][4];
    char tmpout[SZ_LINE+1]; /*to write validation images*/
        
    /*These are already in the parameter structure
    int     Ws              the number of traps < 999999, taken from pctetab read 
    int     q_w[TRAPS];     the run of charge with level  cte->qlevq_data[]
    float   dpde_w[TRAPS];  the run of charge loss with level cte->dpdew_data[]
    
    float   rprof_wt[TRAPS][100]; the emission probability as fn of downhill pixel
    float   cprof_wt[TRAPS][100]; the cumulative probability cprof_t( 1)  = 1. - rprof_t(1)
    */
    
    if (wf3->verbose){
        trlmessage("Starting RSZ to RSC subroutine");
    }
    
    SingleGroup pixz_fff;
    initSingleGroup(&pixz_fff);
    allocSingleGroup(&pixz_fff, COLUMNS, ROWS);
    
    
    /*For reference to jays code, ff_by_col is what's in the scale by column
    
    int   iz_data[ROWS];  column number in raz format
    double scale512[ROWS];      scaling appropriate at row 512 
    double scale1024[ROWS];     scaling appropriate at row 1024
    double scale1536[ROWS];     scaling appropriate at row 1536
    double scale2048[ROWS];     scaling appropriate at row 2048
    */
    
    /*scale by 1 unless the PCTETAB says otherwise, i is the packet num
      This is a safety loop inase not all the columns are populated
      in the reference file*/
    for (i=0; i<COLUMNS;i++){
        ff_by_col[i][0]=1.;
        ff_by_col[i][1]=1.;
        ff_by_col[i][2]=1.;
        ff_by_col[i][3]=1.;
    }
        
    for (i=0;i<COLUMNS;i++){
        j=(int) cte->iz_data[i]; /*which column to scale*/
        ff_by_col[j][0]=cte->scale512[i];
        ff_by_col[j][1]=cte->scale1024[i];
        ff_by_col[j][2]=cte->scale1536[i];
        ff_by_col[j][3]=cte->scale2048[i];
    }
        
    
    
    /*CALCULATE THE CTE CORRECTION FOR EVERY PIXEL
     Index is figured on the final size of the image
     not the current size
    */
    for (i=0; i<COLUMNS; i++){
        for (j=0;j<ROWS; j++){
            ro = j/512.; /*ro can be zero, it's an index*/
            if (ro <0 ) ro=0;
            if (ro > 2.999) ro=2.999; /*only 4 quads, 0 to 3*/
            io = (int)floor(ro); /*force truncation towards 0*/
            cte_j= (float)(j+1) / 2048.; 
            cte_i= ff_by_col[i][io] + (ff_by_col[i][io+1] -ff_by_col[i][io]) * (ro-io);
            Pix(pixz_fff.sci.data,i,j) = cte_i*cte_j;
        }
    }

    /*write the fff image for validation*/
    IODescPtr out;
    Hdr junkhdr;
    initHdr(&junkhdr);
    strcpy(tmpout,wf3->rootname);
    strcat(tmpout,"_pixz_fff.fits");
    out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
    putHeader(out);
    putFloatData(out,&pixz_fff.sci.data);
    closeImage(out);
    sprintf(MsgText,"Wrote FFF file for verification (CTE scalings)  %s",tmpout);
    trlmessage(MsgText);
    freeHdr (&junkhdr);

    /*THIS IS RAZ2RAC_PAR IN JAYS CODE*/    
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
***/

int inverse_cte_blur(SingleGroup *rsz, SingleGroup *rsc, SingleGroup *fff, CTEParams *cte, int verbose, double expstart){
    
    extern int status;
    
    int NREDO, REDO;
    int NITINV, NITCTE; /*looping vars*/
    int i;
    int j,jj;
    float dmod;
    int jmax;
    float cte_ff; /*cte scaling based on observation date*/
    
    /*STARTING DEFAULTS*/
    NITINV=1;
    NITCTE=1;
    
    /*LOCAL IMAGES TO PLAY WITH, THEY WILL REPLACE THE INPUTS*/
    SingleGroup rz; /*pixz_raz*/
    initSingleGroup(&rz);
    allocSingleGroup(&rz, COLUMNS, ROWS);
    
    SingleGroup rc; /*pixz_rac*/
    initSingleGroup(&rc);
    allocSingleGroup(&rc, COLUMNS, ROWS);
    
    SingleGroup pixz_fff; /*pixz_fff*/
    initSingleGroup(&pixz_fff);
    allocSingleGroup(&pixz_fff, COLUMNS, ROWS);
        
    /*use EXPSTART yyyy-mm-dd to determine the CTE scaling
    appropriate for the given date. WFC3/UVIS was
    installed around May 11,2009 and the model was
    constructed to be valid around Sep 3, 2012, a little
    over 3 years after installation*/

    cte_ff= (float) (expstart - cte->cte_date0)/ (cte->cte_date1 - cte->cte_date0);    
    cte->scale_frac=cte_ff;   /*save to param structure for header update*/ 
    
    if(verbose){
        sprintf(MsgText,"cte_ff = %f",cte_ff);
        trlmessage(MsgText);
    }
    
    /*SET UP THE SCALING ARRAY WITH INPUT DATA*/
    for (i=0; i<COLUMNS; i++){
        for (j=0; j< ROWS; j++){
            memcpy(&Pix(rz.sci.data,i,j), &Pix(rsz->sci.data,i,j),sizeof(float));
            memcpy(&Pix(rc.sci.data,i,j), &Pix(rsc->sci.data,i,j),sizeof(float));
            Pix(pixz_fff.sci.data,i,j) = cte_ff * (float)(j+1)/2048.; /*j+1 so not zero*/
        }          
    }
    
    if (verbose){
        /*WRITE THE FFF IMAGE FOR VALIDATION*/
        IODescPtr out;
        Hdr junkhdr;
        char tmpout[SZ_LINE+1];
        initHdr(&junkhdr);
        strcpy(tmpout,"pixz_fff_local.fits");
        out = openOutputImage(tmpout,"",0,&junkhdr,0,0,FITSBYTE);
        putHeader(out);
        putFloatData(out,&pixz_fff.sci.data);
        closeImage(out);
        sprintf(MsgText,"Wrote PIXZ_FFF file for verification (CTE scalings in blur)  %s",tmpout);
        trlmessage(MsgText);
        freeHdr (&junkhdr);
    }     
    
    if (verbose){
        sprintf(MsgText,"cte->thresh=%f, cte->fix_rocr=%i, n_for=%i, n_par=%i",cte->thresh,cte->fix_rocr,cte->n_forward,cte->n_par);
        trlmessage(MsgText);
        trlmessage("Col[2000]  ori\tcor\tchg");
    }
        



#pragma omp parallel for schedule (dynamic) private(dmod, i,\
                j,jj,jmax,REDO,NREDO) \
                shared(cte,rz,pixz_fff)

    for (i=0; i< COLUMNS; i++){           
    /*These will be farmed out when parallel*/    
        float *pix_obsd = (float *) malloc (ROWS *sizeof(float));
        float *pix_modl = (float *) malloc (ROWS *sizeof(float));
        float *pix_curr = (float *) malloc (ROWS *sizeof(float));
        float *pix_init = (float *) malloc (ROWS *sizeof(float));
        float *pix_read = (float *) malloc (ROWS *sizeof(float));
        float *pix_ctef = (float *) malloc (ROWS *sizeof(float));

        /*HORIZONTAL PRE/POST SCAN POPULATION */
        for (j=0; j< ROWS; j++){
            pix_obsd[j]=Pix(rz.sci.data,i,j); /*starts as input RAZ*/
        }
        
        NREDO=0; /*start out not needing to mitigate CRs*/
        REDO=1; /*true - to get into the loop the first time*/
        while (REDO)  { /*replacing goto 9999*/
            REDO=0; /*false*/
            /*STARTING WITH THE OBSERVED IMAGE AS MODEL, ADOPT THE SCALING FOR THIS COLUMN*/
            for (j=0; j<ROWS; j++){
                pix_modl[j] = Pix(rz.sci.data,i,j);
                pix_ctef[j] = Pix(pixz_fff.sci.data,i,j);     
            }
            /*START WITH THE INPUT ARRAY BEING THE LAST OUTPUT
              IF WE'VE CR-RESCALED, THEN IMPLEMENT CTEF*/
            for (NITINV=1; NITINV<cte->n_forward; NITINV++){
                for (j=0; j<ROWS; j++){
                    pix_curr[j]=pix_modl[j];
                    pix_read[j]=pix_modl[j];
                    pix_ctef[j]=Pix(pixz_fff.sci.data,i,j);
                 }

                /*TAKE EACH PIXEL DOWN THE DETECTOR IN NCTENPAR=7*/
                for (NITCTE=1; NITCTE<cte->n_par; NITCTE++){
                    sim_colreadout_l(pix_curr, pix_read, pix_ctef, cte);


                    if (NITCTE == 1 && NITINV ==1 && i==0 && NREDO==0){
                        FILE *readfile = fopen("pix_read_megan.dat", "w");
                        if (readfile == NULL)
                        {
                            printf("Error opening file!\n");
                            exit(1);
                        }

                        FILE *modlfile = fopen("pix_modl_megan.dat", "w");
                        if (modlfile == NULL)
                        {
                            printf("Error opening file!\n");
                            exit(1);
                        }

                        FILE *obsdfile = fopen("pix_obsd_megan.dat", "w");
                        if (obsdfile == NULL)
                        {
                            printf("Error opening file!\n");
                            exit(1);
                        }
    
                        for(j=0; j<ROWS; j++){
                            fprintf(readfile,"%i %f\n",j+1,pix_read[j]);
                            fprintf(modlfile,"%i %f\n",j+1,pix_modl[j]);
                            fprintf(obsdfile,"%i %f\n",j+1,pix_obsd[j]);
                        }
                    trlmessage("Wrote pix check files");
                    fclose(readfile);
                    fclose(modlfile);
                    fclose(obsdfile);
                    trlmessage("Closed pix check files");
                    }
                    /*COPY THE JUST UPDATED READ OUT IMAGE INTO THE INPUT IMAGE*/
                    for (j=0; j< ROWS; j++){
                        pix_curr[j]=pix_read[j];
                    }
                } /* end NITCTE */
                
                    
                /*DAMPEN THE ADJUSTMENT IF IT IS CLOSE TO THE READNOISE, THIS IS
                  AN ADDITIONAL AID IN MITIGATING THE IMPACT OF READNOISE*/
                for (j=0; j< ROWS; j++){
                    dmod = pix_obsd[j] - pix_read[j];
                    if (NITINV < cte->n_forward){ 
                        dmod *=  (dmod*dmod)/((dmod*dmod) + (3.25*3.25));
                    }
                    pix_modl[j] += dmod; /*dampen each pixel as the best is determined*/
                }
            } /*NITINV end*/

            /*LOOK FOR AND DOWNSCALE THE CTE MODEL IF WE FIND
            THE TELL-TALE SIGN OF READOUT CRS BEING OVERSUBTRACTED;
            IF WE FIND ANY THEN GO BACK UP AND RERUN THIS COLUMN
            */
            if (cte->fix_rocr) {
                for (j=10; j< ROWS-1; j++){                        
                    if (  ((pix_modl[j] < cte->thresh) && 
                           (pix_modl[j] - pix_obsd[j] < cte->thresh)) ||
                          ((pix_modl[j] + pix_modl[j+1] < -12) &&
                           (pix_modl[j] + pix_modl[j+1] - pix_obsd[j] - pix_obsd[j+1] < -12)) ||
                          ((pix_modl[j] + pix_modl[j+1] + pix_modl[j+2] < -15) &&
                           (pix_modl[j] + pix_modl[j+1] + pix_modl[j+2] -pix_obsd[j] - 
                                 pix_obsd[j+1] - pix_obsd[j+2] <-15))) {
                            jmax=j;
                            /*GO DOWNSTREAM AND LOOK FOR THE OFFENDING CR*/
                            for (jj=j-10; jj<j;jj++){
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
        } /*replacing goto 9999*/
        
        #pragma omp critical (modl2rac)                 
        for (j=0; j< ROWS; j++){
            memcpy(&Pix(rc.sci.data,i,j),&pix_modl[j],sizeof(float)); /*copy modl to RSC image*/
        }         
        if (verbose){
            if (i ==0){
                sprintf(MsgText,"%15s AMPLIFIER C %2s","*****","*****");
                trlmessage(MsgText);
            }
            if (i==2103){
                sprintf(MsgText,"%15s AMPLIFIER D %2s","*****","*****");
                trlmessage(MsgText);
            }
            if (i==4206){
                sprintf(MsgText,"%15s AMPLIFIER A %2s","*****","*****");
                trlmessage(MsgText);
            }
            if (i==6309){
                sprintf(MsgText,"%15s AMPLIFIER B %2s","*****","*****");
                trlmessage(MsgText);
            }
            
            if ((i+1)%100 == 0){
                sprintf(MsgText,"%d\t%d\t%d\t%d",i+1, 
                        (int)floor(pix_obsd[1999]), 
                        (int)floor(pix_modl[1999]), 
                        ((int) (floor(pix_modl[1999]) - floor(pix_obsd[1999]))));
                trlmessage(MsgText);
            }
        }
    free(pix_obsd);
    free(pix_modl);
    free(pix_curr);
    free(pix_init);
    free(pix_read);
    free(pix_ctef);
                           
              
    } /*end i*/                 


    for (i=0; i< COLUMNS; i++){
        for (j=0; j< ROWS; j++){
            memcpy(&Pix(rsz->sci.data,i,j), &Pix(rz.sci.data,i,j),sizeof(float));
            memcpy(&Pix(rsc->sci.data,i,j), &Pix(rc.sci.data,i,j),sizeof(float));
            memcpy(&Pix(fff->sci.data,i,j), &Pix(pixz_fff.sci.data,i,j),sizeof(float));
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


JDIM == ROWS
WDIM == TRAPS  Ws is the input traps number < 999999
NITs == cte_pars->n_par

These are already in the parameter structure CTEParams
    int     Ws              the number of traps < 999999
    int     q_w[TRAPS];     the run of charge with level  == qlevq_data
    float   dpde_w[TRAPS];  the run of charge loss with level == dpdew_data
    float   rprof_wt[TRAPS][100]; the emission probability as fn of downhill pixel == rprof fits image
    float   cprof_wt[TRAPS][100]; the cumulative probability cprof_t( 1)  = 1. - rprof_t(1)  == cprof fits image
  

   W = wcol_data = trap id
   
   q_w[TRAP] = qlev_q from QPROF  traps as function of packet size = cte->qlevq_data[TRAP]   
   
   pixi (curr), pixo (read) , pixf(cteff) are passed and are 1d arrays which have values for a particular column
*/

int sim_colreadout_l(float *pixi, float *pixo, float *pixf, CTEParams *cte){

    extern int status;
    int j;
    float ftrap;
    int ttrap;
        
    int w;
    float pix_1;
    float padd_2;
    float padd_3;
    float prem_3;
    int pmax;
    float fcarry, fcarry0;
    
    padd_3=0.;
    prem_3=0.;
    padd_2=0.;
    fcarry0=0.0;

    FloatHdrData *rprof;
    FloatHdrData *cprof;
    
    rprof= cte->rprof;
    cprof = cte->cprof;
    

    /*FIGURE OUT WHICH TRAPS WE DON'T NEED TO WORRY ABOUT IN THIS COLUMN*/
    pmax=10;
    for(j=0; j<ROWS; j++){
        pixo[j] = pixi[j];
        if (pixo[j] > pmax)  
            pmax=pixo[j];
    }
       
    /*GO THROUGH THE TRAPS ONE AT A TIME, FROM HIGHEST TO LOWEST Q,
      AND SEE WHEN THEY GET FILLED AND EMPTIED, ADJUST THE PIXELS ACCORDINGLY*/
    for (w = cte->cte_traps -1; w>=0; w--){  
        if (cte->qlevq_data[w] <= pmax) {
            ftrap = 0.0;
            ttrap = cte->cte_len;
            fcarry = fcarry0;
            
            /*GO UP THE COLUMN PIXEL BY PIXEL*/
            for(j=0; j<ROWS;j++){   
                pix_1 = pixo[j];
                
                if ( (ttrap < cte->cte_len) || (pix_1 >= cte->qlevq_data[w] - 1) ){                     
                    if (pixo[j] >= 0){
                        pix_1 = pixo[j] + fcarry; /*shuffle charge in*/
                        fcarry = pix_1 -  (int)floor(pix_1); /*carry the charge remainder*/
                        pix_1 = (int) floor(pix_1); /*reset pixel*/
                    }
                
                    /*HAPPENS AFTER FIRST PASS*/
                    if (j> 0) {
                        if (pixf[j] < pixf[j-1])
                            ftrap= (pixf[j] /  pixf[j-1]) *ftrap;                             
                    }

                    padd_2=0.;
                    if (ttrap < cte->cte_len){
                        ttrap += 1;
                        padd_2 = (rprof->data,w,ttrap) *ftrap;
                    }

                    padd_3 = 0.;
                    prem_3 = 0.;
                    if (pix_1 >=  cte->qlevq_data[w]){
                        prem_3 = ((float) cte->dpdew_data[w] / (float)cte->n_par) * pixf[j];
                        
                        if (ttrap < cte->cte_len)
                            padd_3 = (cprof->data,w,ttrap)*ftrap;
                        ttrap=0;
                        ftrap=prem_3;
                    }  
                    pixo[j] += padd_2 + padd_3 - prem_3;
                    
                } /*replaces trap continue*/
            }/*end if j>0*/                    
        }       /* end if qlevq > pmax, replaces continue*/
    }/*end for w*/
    
    return(status);
    
}


int initCTETrl (char *input, char *output) {

    extern int status;

    char trl_in[SZ_LINE+1];     /* trailer filename for input */
    char trl_out[SZ_LINE+1];    /* output trailer filename */
    int exist;


    char isuffix[] = "_raw";
    char osuffix[] = "_rac";

    int MkName (char *, char *, char *, char *, char *, int);
    void WhichError (int);
    int TrlExists (char *);
    void SetTrlOverwriteMode (int);

    /* Initialize internal variables */
    trl_in[0] = '\0';
    trl_out[0] = '\0';
    exist = EXISTS_UNKNOWN;


    if (MkName (output, osuffix, "_rac", TRL_EXTN, trl_out, SZ_LINE)) {
        WhichError (status);
        sprintf (MsgText, "Couldn't create trailer filename for %s",
                output);
        trlmessage (MsgText);
    }

    /* Test whether the output file already exists */
    exist = TrlExists(trl_out);
    if (exist == EXISTS_YES) {
        /* The output file exists, so we want to overwrite them with
         ** the new trailer comments.  */
        SetTrlOverwriteMode (YES);  
    }

    /* Sets up temp trailer file for output and copies input
     ** trailer file into it.  */
    InitTrlFile (trl_in, trl_out);
    
    return(status);
}




