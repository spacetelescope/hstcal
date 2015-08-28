/* WFC3 -- CTE loss correction for UVIS 

   M. sosey  Aug-2014  Adapted for the pipeline from Jay Andersons CTE correction code for wfc3 UVIS
   raw2raz_wfc3uv.F , an edited file was delivered december 2014, and both are different from the 
   fortran code currently served on the wfc3 website.

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
        
    CTEParams cte_pars; /*STRUCTURE HOLDING THE MODEL PARAMETERS*/
    SingleGroup cd; /*SCI 1*/
    SingleGroup ab; /*SCI 2*/
    SingleGroup raz; /* THE LARGE FORMAT COMBINATION OF CDAB*/
    SingleGroup rsz; /* LARGE FORMAT READNOISE CORRECTED IMAGE */
    SingleGroup rsc; /* CTE CORRECTED*/
    SingleGroup rzc; /*part of the conversion to rac*/
    SingleGroup chg; /*part of the conversion to rac*/
    
    int i,j; /*loop var*/
    int max_threads=1;
    float hardset=0.0000000;

    Bool subarray; /* to verify that no subarray is being used, it's not implemented yet*/

    /*CONTAIN PARALLEL PROCESSING TO A SINGLE THREAD AS USER OPTION*/
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
    IODescPtr out=0;
    char tmpout[SZ_LINE+1]; 
   

    /* COPY COMMAND-LINE ARGUMENTS INTO WF3. */ 
    WF3Init (&wf3);
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

    /* OPEN THE INPUT IMAGES AND GET THE  SCIENCE EXTENSIONS  */
    initSingleGroup (&cd);
    getSingleGroup (wf3.input, 1, &cd);
    if (hstio_err())
        return (status = OPEN_FAILED);

    /*** MAKE SURE THIS IS NOT A SUBARRAY ***/
    if (GetKeyBool (cd.globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
        return (status=KEYWORD_MISSING);

    if (subarray) {
        sprintf(MsgText,"\n**SUBARRAY FOUND!; SUBARRAY images are not yet supported for CTE**\n");
        trlmessage(MsgText);
        status=ERROR_RETURN;
        return(status);
    }   
    
    initSingleGroup (&ab);
    getSingleGroup (wf3.input, 2, &ab);
    if (hstio_err())
        return (status = OPEN_FAILED);

    if (GetKeyBool (ab.globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
        return (status=KEYWORD_MISSING);

    if (subarray) {
        sprintf(MsgText,"SUBARRAY FOUND; **SUBARRAY images are not yet supported for CTE**\n");
        trlmessage(MsgText);
        status=ERROR_RETURN;
        return(status);
    }

    /*READ IN THE CTE PARAMETER TABLE*/
    initCTEParams(&cte_pars);
    if (GetCTEPars (wf3.pctetab.name,&cte_pars))
        return (status);

    if (verbose){
        PrRefInfo ("pctetab", wf3.pctetab.name, wf3.pctetab.pedigree,
                wf3.pctetab.descrip, wf3.pctetab.descrip2);
    }
    
    /*SAVE THE PCTETABLE INFORMATION TO THE HEADER OF THE SCIENCE IMAGE
      AFTER CHECKING TO SEE IF THE USER HAS SPECIFIED ANY CHANGES TO THE
      CTE CODE VARIABLES.
    */
    if (CompareCTEParams(&cd, &cte_pars)){
        return (status);
    }    

    /***SUBTRACT THE CTE BIAS FROM BOTH CHIPS IN PLACE***/
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
      
    /*SET UP THE ARRAYS WHICH WILL BE PASSED AROUND*/  
    initSingleGroup(&raz);
    allocSingleGroup(&raz, RAZ_COLS, RAZ_ROWS);

    initSingleGroup(&rsz);
    allocSingleGroup(&rsz, RAZ_COLS, RAZ_ROWS);

    initSingleGroup(&rsc);
    allocSingleGroup(&rsc, RAZ_COLS, RAZ_ROWS);

    initSingleGroup(&rzc);
    allocSingleGroup(&rzc, RAZ_COLS, RAZ_ROWS);
    
    /*THIS IS THE CHANGE IMAGE FOR VALIDATION*/
    initSingleGroup(&chg);    
    allocSingleGroup(&chg, RAZ_COLS, RAZ_ROWS);

    /*INITIALIZE ARRAYS*/
    for (i=0;i<RAZ_COLS;i++){
        for(j=0;j<RAZ_ROWS;j++){
            memcpy(&Pix(rsz.sci.data,i,j),&hardset,sizeof(float));
            memcpy(&Pix(raz.sci.data,i,j),&hardset,sizeof(float));            
            memcpy(&Pix(rsc.sci.data,i,j),&hardset,sizeof(float));
            memcpy(&Pix(rzc.sci.data,i,j),&hardset,sizeof(float));
            memcpy(&Pix(chg.sci.data,i,j),&hardset,sizeof(float));
        }
    }


    /*CONVERT TO RAZ FORMAT AND CORRECT FOR GAIN*/
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
    trlmessage("CTE: Calculating smooth readnoise image\n");
    
    /***CREATE THE NOISE MITIGATION MODEL ***/
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
            sprintf(MsgText,"RSZ (input for rsz2rsc) format image to check: %s\n",tmpout);
            trlmessage(MsgText);
        }
    } else {
        trlmessage("Only noise model 0 implemented!");
        return (status=ERROR_RETURN);
    }


    /***CONVERT THE READNOISE SMNOOTHED IMAGE TO RSC IMAGE
        THIS IS WHERE THE CTE GETS CALCULATED         ***/
    if (verbose)
        trlmessage("CTE: Converting RSZ to RSC\n");

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
    

    /*** SAVE USEFULL HEADER INFORMATION ***/
    if (cteHistory (&wf3, cd.globalhdr))
        return (status);

    
    /*** CREATE THE FINAL CTE CORRECTED IMAGE, PUT IT BACK INTO ORIGNAL RAW FORMAT***/
    for (i=0;i<RAZ_COLS;i++){
        for(j=0; j<RAZ_ROWS; j++){
            Pix(rzc.sci.data,i,j) = Pix(raz.sci.data,i,j) + (Pix(rsc.sci.data,i,j) - Pix(rsz.sci.data,i,j));
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
    
    /*REMOVE THE  GAIN*/
   for (i=0;i<RAZ_COLS;i++){
        for(j=0; j<RAZ_ROWS; j++){
            Pix(rzc.sci.data,i,j) /= wf3.ccdgain;
        }
    }
    
    /*BACK TO NORMAL FORMATTING*/
    undosciRAZ(&cd,&ab,&rzc);
    
    /*UPDATE THE OUTPUT HEADER ONE FINAL TIME*/
    PutKeyFlt (cd.globalhdr, "PCTEFRAC", cte_pars.scale_frac,"CTE scaling fraction based on expstart");
    trlmessage("PCTEFRAC saved to header");
    
    /*SAVE THE NEW RAW FILE WITH UPDATED SCIENCE ARRAYS AND PRIMARY HEADER TO RAC*/
    sprintf(MsgText,"Writing cd[sci,%i] to %s",cd.group_num,output);
    trlmessage(MsgText);
    putSingleGroup(output,cd.group_num, &cd,0);
    
    sprintf(MsgText,"Writing ab[sci,%i] to %s",ab.group_num,output);
    trlmessage(MsgText);
    putSingleGroup(output,ab.group_num, &ab,0);


    /** CLEAN UP **/    
    freeSingleGroup(&rzc);
    freeSingleGroup(&rsc);
    freeSingleGroup(&chg);
    freeSingleGroup(&raz);
    freeSingleGroup(&rsz);
    
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
            for (j=0;j<RAZ_ROWS; j++){
                Pix(raz->sci.data,i+(k-1)*subcol,j) -= bias_post[k];
            }
        }   
    }
    /*done just for reporting*/
    findPreScanBias(raz, bias_pre, bsig_pre);

    /* convert to electrons, multiply by the gain */
    for (i=0; i<RAZ_COLS;i++){
        for (j=0;j<RAZ_ROWS; j++){
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
    

    int subcol = RAZ_COLS/4;
    int npix; /*track array size for resistant mean*/
    rmean=0.;
    rsigma=0.;
    npix=0;

    /*init plist*/
    for (i=0;i<55377;i++){
        plist[i]=0.;
    }
    
    trlmessage("\nPost scan bias measures:\n");
    
    for (k=1;k<=4;k++){  /*for each quadrant cdab = 0123*/
        npix=0;
        for (i=RAZ_ROWS+5;i<= subcol-1; i++){ /*specific area for post scan bias pixels*/
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
    
    int subcol = RAZ_COLS/4;
    int npix=0; /*track array size for resistant mean*/
    rmean=0.;
    rsigma=0.;
    
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
    char tmpout[SZ_LINE+1];
    float comparison=999.9;
    
    if (wf3->verbose)
        trlmessage("Starting RAZ -> RSZ\n");
    
    int i, j, NIT; /*loop variables*/
    int imid;
    float dptr=0.;
    double  rms=0.0;
    double  rmsu=0.0;
    double nrms, nrmsu;
    float hardset=0.0000000;
    double dblzero=0.00000000;

    /*1D ARRAYS FOR CENTRAL AND NEIGHBORING RAZ_COLS*/
    float obs_loc[3][RAZ_ROWS] ; 
    float rsz_loc[3][RAZ_ROWS] ;
        
    nrms=0.;
    nrmsu=0.;
    NIT=1;
    
    /*ALL ELEMENTS TO ZERO*/
    for(i=0;i<3;i++){
        for (j=0; j<RAZ_ROWS; j++){
            obs_loc[i][j]=hardset;
            rsz_loc[i][j]=hardset;
        }
    }
        
    /***INITIALIZE THE LOCAL IMAGE GROUPS***/
    SingleGroup rnz;
    initSingleGroup(&rnz);
    allocSingleGroup(&rnz, RAZ_COLS, RAZ_ROWS);
        
    SingleGroup zadj;
    initSingleGroup(&zadj);
    allocSingleGroup(&zadj, RAZ_COLS, RAZ_ROWS);  


    /*COPY THE RAZ IMAGE INTO THE RSZ OUTPUT IMAGE
      AND INITIALIZE THE OTHER IMAGES*/
    for(i=0;i<RAZ_COLS;i++){
        for (j=0;j<RAZ_ROWS;j++){
            memcpy( &Pix(rsz->sci.data,i,j), &Pix(raz->sci.data,i,j),sizeof(float) );
            memcpy( &Pix(rnz.sci.data,i,j), &hardset, sizeof(float));
            memcpy( &Pix(zadj.sci.data,i,j), &hardset, sizeof(float));
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
             

    sprintf(MsgText,"RNSIG=%3.2f\nNIT\tCHGrms\t\ticol\tORIG\t\tDIFF\t\tRSZ\n",rnsig);
    trlmessage(MsgText);
    memcpy(&rms,&dblzero,sizeof(double));
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
                memcpy(&obs_loc[0][j],&Pix(raz->sci.data,imid-1,j),sizeof(float));
                memcpy(&obs_loc[1][j],&Pix(raz->sci.data,imid,j),sizeof(float));
                memcpy(&obs_loc[2][j],&Pix(raz->sci.data,imid+1,j),sizeof(float));

                memcpy(&rsz_loc[0][j],&Pix(rsz->sci.data,imid-1,j),sizeof(float));
                memcpy(&rsz_loc[1][j],&Pix(rsz->sci.data,imid,j),sizeof(float));
                memcpy(&rsz_loc[2][j],&Pix(rsz->sci.data,imid+1,j),sizeof(float));
            }
            for (j=0; j<RAZ_ROWS; j++){  
                find_dadj(1+i-imid,j, obs_loc, rsz_loc, rnsig, &dptr);
                memcpy(&Pix(zadj.sci.data,i,j),&dptr,sizeof(float));
                    if (j==1999 && i==19){
                        sprintf(MsgText,"%2i\t%8.4f\t%4i\t%8.4f\t%8.4f\t%8.4f",NIT,rms,i+1,
                            fminf(Pix(raz->sci.data,i,j),comparison),
                            (Pix(raz->sci.data,i,j) - Pix(rsz->sci.data,i,j)),
                            fminf(Pix(rsz->sci.data,i,j), comparison));
                        trlmessage(MsgText);
                    }
                 
            }
        } /*end the parallel for*/
    

        /*NOW GO OVER ALL THE RAZ_COLS AND RAZ_ROWS AGAIN TO SCALE THE PIXELS */
        for(i=0; i<RAZ_COLS;i++){
            for(j=0; j<RAZ_ROWS; j++){
                Pix(rsz->sci.data,i,j) += Pix(zadj.sci.data,i,j)*0.75;
                Pix(rnz.sci.data,i,j) = (Pix(raz->sci.data,i,j) - Pix(rsz->sci.data,i,j));
            }                
        }

        memcpy(&rms,&dblzero,sizeof(double));
        memcpy(&nrms,&dblzero,sizeof(double));
               
       #pragma omp parallel for schedule(static) \
        private(i,j,rmsu,nrmsu) \
        shared(raz,rsz,rms,rnsig,nrms)
       for(j=0; j<RAZ_ROWS; j++){
            memcpy(&nrmsu,&dblzero,sizeof(double));
            memcpy(&rmsu,&dblzero,sizeof(double));
            for(i = 0;i<RAZ_COLS; i++){
                if ( (fabs(Pix(raz->sci.data,i,j)) > 0.1) || 
                     (fabs(Pix(rsz->sci.data,i,j)) > 0.1) ){
                    rmsu  += (double) (Pix(rnz.sci.data,i,j) *  Pix(rnz.sci.data,i,j));
                    nrmsu += 1.;
                }
            }
            
            #pragma omp critical (rms)
            {rms  += rmsu;
            nrms += nrmsu;}
        }
        
        rms = sqrt(rms/nrms);
        /*epsilon type comparison*/
        if ( (rnsig-rms) < 0.) break; /*this exits the NIT for loop*/
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
        sprintf(MsgText,"Wrote RNZ (readnoise) file for verification  %s",tmpout);
        trlmessage(MsgText);
        freeHdr(&junkhdr);
    }            

    freeSingleGroup(&zadj);
    freeSingleGroup(&rnz);
 
    
    return (status);
}

int find_dadj(int i ,int j, float obsloc[][RAZ_ROWS], float rszloc[][RAZ_ROWS], float rnsig, float *d){
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

    double mval;
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
    
    mval = (double)rszloc[i][j];
    dval0  = (double)obsloc[i][j] - mval;
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

    dval9 =(double)dval9 / 9.;
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
        dmod2 =  (double)rszloc[i][j+1] - mval;

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


/*** THIS ROUTINE PERFORMS THE CTE CORRECTIONS 
     rsz is the readnoise smoothed image
     rsc is the output image
     
     rac = raw + (rsc-rsz) / gain 
     
***/
int rsz2rsc(WF3Info *wf3, SingleGroup *rsz, SingleGroup *rsc, CTEParams *cte) {

    extern int status;
    
    int i,j;
    float cte_i=0.;
    float cte_j=0.;
    float ro=0;
    int io=0;
    float ff_by_col[RAZ_COLS][4];
    char tmpout[SZ_LINE+1]; /*to write validation images*/
    float hardset=0.;
    
    /*These are already in the parameter structure
    int     Ws              the number of traps < 999999, taken from pctetab read 
    int     q_w[TRAPS];     the run of charge with level  cte->qlevq_data[]
    float   dpde_w[TRAPS];  the run of charge loss with level cte->dpdew_data[]
    
    float   rprof_wt[TRAPS][100]; the emission probability as fn of downhill pixel, TRAPS=999
    float   cprof_wt[TRAPS][100]; the cumulative probability cprof_t( 1)  = 1. - rprof_t(1)
    */
    
    if (wf3->verbose){
        trlmessage("Starting RSZ to RSC ...");
    }
    
    SingleGroup pixz_fff;
    initSingleGroup(&pixz_fff);
    allocSingleGroup(&pixz_fff, RAZ_COLS, RAZ_ROWS);
    
    for(i=0; i<RAZ_COLS;i++){
        for(j=0; j<RAZ_ROWS; j++){
            memcpy(&Pix(pixz_fff.sci.data,i,j),&hardset,sizeof(float));
        }
    }
    
    /*FOR REFERENCE TO JAYS CODE, FF_BY_COL IS WHAT'S IN THE SCALE BY COLUMN
    
    int   iz_data[RAZ_ROWS];  column number in raz format
    double scale512[RAZ_ROWS];      scaling appropriate at row 512 
    double scale1024[RAZ_ROWS];     scaling appropriate at row 1024
    double scale1536[RAZ_ROWS];     scaling appropriate at row 1536
    double scale2048[RAZ_ROWS];     scaling appropriate at row 2048
    */
    
    /*SCALE BY 1 UNLESS THE PCTETAB SAYS OTHERWISE, I IS THE PACKET NUM
      THIS IS A SAFETY LOOP INCASE NOT ALL THE COLUMNS ARE POPULATED
      IN THE REFERENCE FILE*/
    for (i=0; i<RAZ_COLS;i++){
        ff_by_col[i][0]=1.;
        ff_by_col[i][1]=1.;
        ff_by_col[i][2]=1.;
        ff_by_col[i][3]=1.;
    }
        
    for (i=0;i<RAZ_COLS;i++){
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
    for (i=0; i<RAZ_COLS; i++){
        for (j=0;j<RAZ_ROWS; j++){
            ro = (float)j/512.; /*ro can be zero, it's an index*/
            if (ro <0 ) ro=0.;
            if (ro > 2.999) ro=2.999; /*only 4 quads, 0 to 3*/
            io = (int) floor(ro); /*force truncation towards 0 for pos numbers*/
            cte_j= (float)(j+1) / 2048.; 
            cte_i= ff_by_col[i][io] + (ff_by_col[i][io+1] -ff_by_col[i][io]) * (ro-(float)io);
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
     
***/

int inverse_cte_blur(SingleGroup *rsz, SingleGroup *rsc, SingleGroup *fff, CTEParams *cte, int verbose, double expstart){
    
    extern int status;
        
    /*looping vars*/
    int NREDO, REDO;
    int NITINV, NITCTE; 
    int i;
    int j,jj;
    float dmod;
    int jmax;
    
    float cte_ff; /*cte scaling based on observation date*/
    
    /*STARTING DEFAULTS*/
    NITINV=1;
    NITCTE=1;
    cte_ff=0.0;
    jmax=0;
    dmod=0.;
    
    /*LOCAL IMAGES TO PLAY WITH, THEY WILL REPLACE THE INPUTS*/
    SingleGroup rz; /*pixz_raz*/
    initSingleGroup(&rz);
    allocSingleGroup(&rz, RAZ_COLS, RAZ_ROWS);
    
    SingleGroup rc; /*pixz_rac*/
    initSingleGroup(&rc);
    allocSingleGroup(&rc, RAZ_COLS, RAZ_ROWS);
    
    SingleGroup pixz_fff; /*pixz_fff*/
    initSingleGroup(&pixz_fff);
    allocSingleGroup(&pixz_fff, RAZ_COLS, RAZ_ROWS);
        
    /*USE EXPSTART YYYY-MM-DD TO DETERMINE THE CTE SCALING
    APPROPRIATE FOR THE GIVEN DATE. WFC3/UVIS WAS
    INSTALLED AROUND MAY 11,2009 AND THE MODEL WAS
    CONSTRUCTED TO BE VALID AROUND SEP 3, 2012, A LITTLE
    OVER 3 YEARS AFTER INSTALLATION*/

    cte_ff= (float) (expstart - cte->cte_date0)/ (cte->cte_date1 - cte->cte_date0);    
    cte->scale_frac=cte_ff;   /*save to param structure for header update*/ 
    
    if(verbose){
        sprintf(MsgText,"cte_ff (scaling fraction by date) = %f",cte_ff);
        trlmessage(MsgText);
    }
    
    /*SET UP THE SCALING ARRAY WITH INPUT DATA*/
    for (i=0; i<RAZ_COLS; i++){
        for (j=0; j< RAZ_ROWS; j++){
            memcpy(&Pix(rz.sci.data,i,j), &Pix(rsz->sci.data,i,j),sizeof(float));
            Pix(pixz_fff.sci.data,i,j) = cte_ff * (float)(j+1)/2048.; /*j+1 so not zero*/
        }          
    }
    
    if (verbose){
        sprintf(MsgText,"cte->thresh=%5.3g, cte->fix_rocr=%i, n_for=%i, n_par=%i",cte->thresh,cte->fix_rocr,cte->n_forward,cte->n_par);
        trlmessage(MsgText);
    }
    
    trlmessage("Col\n[2000]  orig\tcorr\tdiff");

    /*DEFINE TO MAKE PRIVATE IN PARALLEL RUN*/
    float setzero=0.;
    
    float *pix_obsd = &setzero;   
    float *pix_modl = &setzero;   
    float *pix_curr = &setzero;   
    float *pix_init = &setzero;   
    float *pix_read = &setzero;   
    float *pix_ctef = &setzero;   

    #pragma omp parallel for schedule (dynamic,1) \
        private(dmod,i,j,jj,jmax,REDO,NREDO, \
          pix_obsd,pix_modl,pix_curr,pix_init,\
          pix_read,pix_ctef,NITINV,NITCTE)\
        shared(rc,rz,cte,pixz_fff)

    for (i=0; i< RAZ_COLS; i++){           
       pix_obsd = (float *) calloc(RAZ_ROWS, sizeof(float));   
       pix_modl = (float *) calloc(RAZ_ROWS, sizeof(float));   
       pix_curr = (float *) calloc(RAZ_ROWS, sizeof(float));   
       pix_init = (float *) calloc(RAZ_ROWS, sizeof(float));   
       pix_read = (float *) calloc(RAZ_ROWS, sizeof(float));   
       pix_ctef = (float *) calloc(RAZ_ROWS, sizeof(float));   

        /*HORIZONTAL PRE/POST SCAN POPULATION */
        for (j=0; j< RAZ_ROWS; j++){
            pix_obsd[j] = Pix(rz.sci.data,i,j); /*starts as input RAZ*/
        }
        
        NREDO=0; /*start out not needing to mitigate CRs*/
        REDO=0; /*false*/
        do { /*replacing goto 9999*/
            /*STARTING WITH THE OBSERVED IMAGE AS MODEL, ADOPT THE SCALING FOR THIS COLUMN*/
            for (j=0; j<RAZ_ROWS; j++){
                pix_modl[j] = Pix(rz.sci.data,i,j);
                pix_ctef[j] = Pix(pixz_fff.sci.data,i,j);     
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
                        dmod *= (dmod*dmod) /((dmod*dmod) + powf(cte->rn_amp , 2));
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
             DEAL WITH MANY OF THEM. FOR IMAGES THAT HAVE BACKGROUND GREAT THAN 10 OR SO,
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
        
        #pragma omp critical (cte)        
        for (j=0; j< RAZ_ROWS; j++){
             Pix(rc.sci.data,i,j)=pix_modl[j]; 
        }         
        
        if (i ==0){
            sprintf(MsgText,"%15s AMPLIFIER C %2s","*****","*****");
            trlmessage(MsgText);
        }
        if (i==2104){
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
                   (int) roundf(pix_obsd[1999]), 
                    (int) roundf(pix_modl[1999]), 
                    (int)(roundf(pix_modl[1999]) - roundf(pix_obsd[1999])));
            trlmessage(MsgText);
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

int sim_colreadout_l(float *pixi, float *pixo, float *pixf, CTEParams *cte){

    extern int status;
    int j;
    int ttrap;
        
    int w;
    float ftrap;
    float pix_1;
    float padd_2;
    float padd_3;
    float prem_3;
    float pmax;
    float fcarry;
    
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
            
            /*sprintf(MsgText,"QPROF[%i]=%f\t pmax=%f\n",w,cte->qlevq_data[w],pmax);
            trlmessage(MsgText);*/
            
            ftrap = 0.0e0;
            ttrap = cte->cte_len; /*for referencing the image at 0*/
            fcarry = 0.0e0;
            
            /*GO UP THE COLUMN PIXEL BY PIXEL*/
            for(j=0; j<RAZ_ROWS;j++){   
                pix_1 = pixo[j];
                
                
                if ( (ttrap < cte->cte_len) || ( pix_1 >= cte->qlevq_data[w] - 1. ) ){                     
                    if (pixo[j] >= 0){
                        pix_1 = pixo[j] + fcarry; /*shuffle charge in*/
                        fcarry = pix_1 - floor(pix_1); /*carry the charge remainder*/
                        pix_1 = floor(pix_1); /*reset pixel*/
                    }    

                    /*HAPPENS AFTER FIRST PASS*/
                    /*SHUFFLE CHARGE IN*/
                    if (j> 0) {
                        if (pixf[j] < pixf[j-1])
                            ftrap *= (pixf[j] /  pixf[j-1]);                             
                    }
                  
                    /*RELEASE THE CHARGE*/
                    padd_2=0.0;
                    if (ttrap <cte->cte_len){
                        ttrap += 1;
                        padd_2 = Pix(rprof->data,w,ttrap-1) *ftrap;
                        /*sprintf(MsgText,"rprof(w:%i,ttr:%i)=%f\tpadd_2=%f\tftrap=%f\n",w,ttrap-1,Pix(rprof->data,w,ttrap-1),padd_2,ftrap);
                        trlmessage(MsgText);*/
                    }

                    padd_3 = 0.0;
                    prem_3 = 0.0;
                    if ( pix_1 >= cte->qlevq_data[w] ){
                        prem_3 =  cte->dpdew_data[w] / cte->n_par * pixf[j];  /*dpdew is 1 in file */                      
                        if (ttrap < cte->cte_len)
                            padd_3 = Pix(cprof->data,w,ttrap-1)*ftrap;
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


    int MkName (char *, char *, char *, char *, char *, int);
    int TrlExists (char *);
    void SetTrlOverwriteMode (int);

    /* Initialize internal variables */
    trl_in[0] = '\0';
    trl_out[0] = '\0';
    exist = EXISTS_UNKNOWN;


    if (MkName (output, "_rac", "", TRL_EXTN, trl_out, SZ_LINE)) {
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




