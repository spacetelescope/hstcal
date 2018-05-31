#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

# ifdef _OPENMP
#include <omp.h>
# endif

#include "hstcal_memory.h"
#include "ctegen2.h"
#include "hstcalerr.h"
#include "hstcal.h"
#include "trlbuf.h"

static void setAtomicFlag(Bool * atom)
{
    if (!atom)
        return;
#ifdef _OPENMP
    #pragma omp critical(critSecAtomicFlag)
#endif
    {
        *atom = True;
    }
}
static void setAtomicInt(int * atom, const int value)
{
    if (!atom)
        return;
#ifdef _OPENMP
    #pragma omp critical(critSecAtomicInt)
#endif
    {
        *atom = value;
    }
}
int forwardModel(const SingleGroup * input, SingleGroup * output, SingleGroup * trapPixelMap, CTEParamsFast * ctePars)
{
    extern int status;

   if (!input || !output || !trapPixelMap || !ctePars)
       return (status = ALLOCATION_PROBLEM);

   //WARNING - assumes column major storage order
   assert(trapPixelMap->sci.data.storageOrder == COLUMNMAJOR);
   assert(input->sci.data.storageOrder == COLUMNMAJOR);
   output->sci.data.storageOrder = COLUMNMAJOR;

   const unsigned nRows = output->sci.data.ny;
   const unsigned nColumns = output->sci.data.nx;
   const FloatTwoDArray * cteRprof  = &ctePars->rprof->data;
   const FloatTwoDArray * cteCprof = &ctePars->cprof->data;

   Bool allocationFail = False;
   Bool runtimeFail = False;
#ifdef _OPENMP
   #pragma omp parallel shared(input, output, ctePars, cteRprof, cteCprof, trapPixelMap, allocationFail, runtimeFail, status)
#endif
   {
       int localStatus = HSTCAL_OK; //Note: used to set extern int status atomically, note global status takes last set value
       //Thread local pointer register
       PtrRegister localPtrReg;
       initPtrRegister(&localPtrReg);

       double * model = malloc(sizeof(*model)*nRows);
       addPtr(&localPtrReg, model, &free);
       if (!model)
           setAtomicFlag(&allocationFail);

       float * traps = NULL;

       //Allocate all local memory before anyone proceeds
#ifdef _OPENMP
       #pragma omp barrier
#endif
       if (!allocationFail)
       {
           {unsigned j;
#ifdef _OPENMP
           #pragma omp for schedule(dynamic)
#endif
           for (j = 0; j < nColumns; ++j)
           {
               // Can't use memcpy as diff types
               // Do in place (in a distributed context)
               {unsigned i;
               for (i = 0; i < nRows; ++i)
                   model[i] = PixColumnMajor(input->sci.data,i,j);
               }

               traps = &(PixColumnMajor(trapPixelMap->sci.data, 0, j));

               if ((localStatus = simulateColumnReadout(model, traps, ctePars, cteRprof, cteCprof, nRows, ctePars->n_par)))
               {
                   setAtomicFlag(&runtimeFail);
                   setAtomicInt(&status, localStatus);
               }
               // Update source array
               // Can't use memcpy as arrays of diff types
               {unsigned i;
               for (i = 0; i < nRows; ++i)
                   PixColumnMajor(output->sci.data, i, j) = model[i];
               }
           }} //end loop over columns
       }
       freeOnExit(&localPtrReg);
   }// close scope for #pragma omp parallel
   if (allocationFail)
   {
       sprintf(MsgText, "Out of memory in inverseCTEBlur()");
       trlerror(MsgText);
       return (status = OUT_OF_MEMORY);
   }
   if (runtimeFail)
   {
       sprintf(MsgText, "Runtime fail in inverseCTEBlur()");
       trlerror(MsgText);
       return status;
   }
   return (status);
}
int inverseCTEBlur(const SingleGroup * input, SingleGroup * output, SingleGroup * trapPixelMap, CTEParamsFast * ctePars)
{
    extern int status;

    if (!input || !output || !trapPixelMap || !ctePars)
        return (status = ALLOCATION_PROBLEM);

    //WARNING - assumes column major storage order
    assert(trapPixelMap->sci.data.storageOrder == COLUMNMAJOR);
    assert(input->sci.data.storageOrder == COLUMNMAJOR);
    output->sci.data.storageOrder = COLUMNMAJOR;

    const unsigned nRows = output->sci.data.ny;
    const unsigned nColumns = output->sci.data.nx;
    const double rnAmp2 = ctePars->rn_amp * ctePars->rn_amp;
    const FloatTwoDArray * cteRprof  = &ctePars->rprof->data;
    const FloatTwoDArray * cteCprof = &ctePars->cprof->data;

    Bool allocationFail = False;
    Bool runtimeFail = False;
#ifdef _OPENMP
    #pragma omp parallel shared(input, output, ctePars, cteRprof, cteCprof, trapPixelMap, allocationFail, runtimeFail, status)
#endif
    {
        int localStatus = HSTCAL_OK; //Note: used to set extern int status atomically, note global status takes last set value
        //Thread local pointer register
        PtrRegister localPtrReg;
        initPtrRegister(&localPtrReg);

        double * model = malloc(sizeof(*model)*nRows);
        addPtr(&localPtrReg, model, &free);
        if (!model)
            setAtomicFlag(&allocationFail);

        double * tempModel = NULL;
        if (!allocationFail)
            tempModel = malloc(sizeof(*tempModel)*nRows);
        addPtr(&localPtrReg, tempModel, &free);
        if (!tempModel)
            setAtomicFlag(&allocationFail);

        double * observed = NULL;
        if (!allocationFail)
            observed = malloc(sizeof(*observed)*nRows);
        addPtr(&localPtrReg, observed, &free);
        if (!observed)
            setAtomicFlag(&allocationFail);

        float * traps = NULL;

        //Allocate all local memory before anyone proceeds
#ifdef _OPENMP
        #pragma omp barrier
#endif
        if (!allocationFail)
        {
            Bool localOK = True;
            {unsigned j;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic)
#endif
            for (j = 0; j < nColumns; ++j)
            {
                // Can't use memcpy as diff types
                {unsigned i;
                for (i = 0; i < nRows; ++i)
                    observed[i] = PixColumnMajor(input->sci.data,i,j);
                }

                traps = &(PixColumnMajor(trapPixelMap->sci.data, 0, j));
                unsigned NREDO = 0;
                Bool REDO;
                do
                {
                    REDO = False; /*START OUT NOT NEEDING TO MITIGATE CRS*/
                    /*STARTING WITH THE OBSERVED IMAGE AS MODEL, ADOPT THE SCALING FOR THIS COLUMN*/
                    memcpy(model, observed, nRows*sizeof(*observed));

                    /*START WITH THE INPUT ARRAY BEING THE LAST OUTPUT
                      IF WE'VE CR-RESCALED, THEN IMPLEMENT CTEF*/
                    {unsigned NITINV;
                    for (NITINV = 1; NITINV <= ctePars->n_forward - 1; ++NITINV)
                    {
                        memcpy(tempModel, model, nRows*sizeof(*model));
                        if ((localStatus = simulateColumnReadout(model, traps, ctePars, cteRprof, cteCprof, nRows, ctePars->n_par)))
                        {
                            setAtomicFlag(&runtimeFail);
                            setAtomicInt(&status, localStatus);
                            localOK = False;
                            break;
                        }

                        //Now that the updated readout has been simulated, subtract this from the model
                        //to reproduce the actual image, without the CTE trails.
                        //Whilst doing so, DAMPEN THE ADJUSTMENT IF IT IS CLOSE TO THE READNOISE, THIS IS
                        //AN ADDITIONAL AID IN MITIGATING THE IMPACT OF READNOISE
                        {unsigned i;
                        for (i = 0; i < nRows; ++i)
                        {
                            double delta = model[i] - observed[i];
                            double delta2 = delta * delta;

                            //DAMPEN THE ADJUSTMENT IF IT IS CLOSE TO THE READNOISE
                            delta *= delta2 / (delta2 + rnAmp2);

                            //Now subtract the simulated readout
                            model[i] = tempModel[i] - delta;
                        }}
                    }}
                    if (!localOK)
                        break;

                    //Do the last forward iteration but don't dampen... no idea why???
                    memcpy(tempModel, model, sizeof(*model)*nRows);
                    if ((localStatus = simulateColumnReadout(model, traps, ctePars, cteRprof, cteCprof, nRows, ctePars->n_par)))
                    {
                        setAtomicFlag(&runtimeFail);
                        setAtomicInt(&status, localStatus);
                        localOK = False;
                        break;
                    }
                    //Now subtract the simulated readout
                    {unsigned i;
                    for (i = 0; i < nRows; ++i)
                        model[i] = tempModel[i] - (model[i] - observed[i]);
                    }

                    REDO = ctePars->fix_rocr ? correctCROverSubtraction(traps, model, observed, nRows,
                            ctePars->thresh) : False;

                } while (localOK && REDO && ++NREDO < 5); //If really wanting 5 re-runs then use NREDO++

                // Update source array
                // Can't use memcpy as arrays of diff types
                {unsigned i;
                for (i = 0; i < nRows; ++i)
                    PixColumnMajor(output->sci.data, i, j) = model[i];
                }
            }} //end loop over columns
        }
        freeOnExit(&localPtrReg);
    }// close scope for #pragma omp parallel
    if (allocationFail)
    {
        sprintf(MsgText, "Out of memory in inverseCTEBlur()");
        trlerror(MsgText);
        return (status = OUT_OF_MEMORY);
    }
    if (runtimeFail)
    {
        sprintf(MsgText, "Runtime fail in inverseCTEBlur()");
        trlerror(MsgText);
        return status;
    }
    return (status);
}

int simulatePixelReadout_v1_1(double * const pixelColumn, const float * const traps, const CTEParamsFast * const ctePars,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows)
{
    //For performance this does not NULL check passed in ptrs

    double chargeToAdd;
    double extraChargeToAdd;
    double chargeToRemove;
    double pixel;
    double releasedFlux;
    double trappedFlux;
    int nTransfersFromTrap;

    /*FIGURE OUT WHICH TRAPS WE DON'T NEED TO WORRY ABOUT IN THIS COLUMN
      PMAX SHOULD ALWAYS BE POSITIVE HERE*/
    //Look into whether this really has to be computed each iteration?
    //Since this is simulating the readout and thus moving pixels down and out, pmax can only get smaller with
    //each pixel transfer, never greater.
    double maxPixel = 10;
    {unsigned i;
    for (i = 0; i < nRows; ++i)
        maxPixel = pixelColumn[i] > maxPixel ? pixelColumn[i] : maxPixel;
    }

    //Find highest charge trap to not exceed i.e. map pmax to an index
    unsigned maxChargeTrapIndex = ctePars->cte_traps-1;
    {int w;
    for (w = maxChargeTrapIndex; w >= 0; --w)//go up or down? (if swap, change below condition)
    {
        if (ctePars->qlevq_data[w] <= maxPixel)//is any of this even needed or can we just directly map?
        {
            maxChargeTrapIndex = w;
            break;
        }
    }}

    /*GO THROUGH THE TRAPS ONE AT A TIME, FROM HIGHEST TO LOWEST Q,
      AND SEE WHEN THEY GET FILLED AND EMPTIED, ADJUST THE PIXELS ACCORDINGLY*/
    {int w;
    for (w = maxChargeTrapIndex; w >= 0; --w)
    {
        nTransfersFromTrap = ctePars->cte_len; //for referencing the image at 0
        trappedFlux = 0;
        releasedFlux = 0;

        /*GO UP THE COLUMN PIXEL BY PIXEL*/
        {unsigned i;
        for (i = 0; i < nRows; ++i)
        {
            pixel = pixelColumn[i];
            Bool isInsideTrailLength = nTransfersFromTrap < ctePars->cte_len;
            Bool isAboveChargeThreshold = pixel >= ctePars->qlevq_data[w] - 1.;

            if (!isInsideTrailLength && !isAboveChargeThreshold)
                continue;

            if (pixelColumn[i] >= 0 )//seems a shame to check this every iteration
            {
                pixel = pixelColumn[i] + releasedFlux; /*shuffle charge in*/
                double floored = floor(pixel);
                releasedFlux = pixel - floored; /*carry the charge remainder*/
                pixel = floored; /*reset pixel*/
            }

            /*HAPPENS AFTER FIRST PASS*/
            /*SHUFFLE CHARGE IN*/
            //move out of loop to separate instance?
            if (i > 0)
            {
                if ((double)traps[i] < (double)traps[i-1])
                    trappedFlux *= ((double)traps[i] / (double)traps[i-1]);
            }

            /*RELEASE THE CHARGE*/
            chargeToAdd = 0;
            if (isInsideTrailLength)
            {
                ++nTransfersFromTrap;
                chargeToAdd = rprof->data[w*rprof->ny + nTransfersFromTrap-1] * trappedFlux;
            }

            extraChargeToAdd = 0;
            chargeToRemove = 0;
            if (pixel >= ctePars->qlevq_data[w])
            {
                chargeToRemove =  ctePars->dpdew_data[w] / ctePars->n_par * (double)traps[i];  /*dpdew is 1 in file */
                if (nTransfersFromTrap < ctePars->cte_len)
                    extraChargeToAdd = cprof->data[w*cprof->ny + nTransfersFromTrap-1] * trappedFlux; //ttrap-1 may not be the same index as ref'd in rprof???
                nTransfersFromTrap = 0;
                trappedFlux = chargeToRemove;
            }

            pixelColumn[i] += chargeToAdd + extraChargeToAdd - chargeToRemove;
        }} //end for i
    }} //end for w
    return HSTCAL_OK;
}

int simulatePixelReadout_v1_2(double * const pixelColumn, const float * const traps, const CTEParamsFast * const ctePars,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows)
{
    //NOTE: this version of the function, simulatePixelReadout, matches Jay Anderson's update
    //to the algorithm (https://github.com/spacetelescope/hstcal/issues/48).
    //For performance this does not NULL check passed in ptrs

    double chargeToAdd;
    double extraChargeToAdd;
    double chargeToRemove;
    double pixel;
    double trappedFlux;
    int nTransfersFromTrap;

    /*FIGURE OUT WHICH TRAPS WE DON'T NEED TO WORRY ABOUT IN THIS COLUMN
      PMAX SHOULD ALWAYS BE POSITIVE HERE*/
    //Look into whether this really has to be computed each iteration?
    //Since this is simulating the readout and thus moving pixels down and out, pmax can only get smaller with
    //each pixel transfer, never greater.
    double maxPixel = 10;
    {unsigned i;
    for (i = 0; i < nRows; ++i)
        maxPixel = pixelColumn[i] > maxPixel ? pixelColumn[i] : maxPixel;
    }

    //Find highest charge trap to not exceed i.e. map pmax to an index
    unsigned maxChargeTrapIndex = ctePars->cte_traps-1;
    {int w;
    for (w = maxChargeTrapIndex; w >= 0; --w)//go up or down? (if swap, change below condition)
    {
        if (ctePars->qlevq_data[w] <= maxPixel)//is any of this even needed or can we just directly map?
        {
            maxChargeTrapIndex = w;
            break;
        }
    }}

    /*GO THROUGH THE TRAPS ONE AT A TIME, FROM HIGHEST TO LOWEST Q,
      AND SEE WHEN THEY GET FILLED AND EMPTIED, ADJUST THE PIXELS ACCORDINGLY*/
    {int w;
    for (w = maxChargeTrapIndex; w >= 0; --w)
    {
        nTransfersFromTrap = ctePars->cte_len; //for referencing the image at 0
        trappedFlux = 0;

        /*GO UP THE COLUMN PIXEL BY PIXEL*/
        {unsigned i;
        for (i = 0; i < nRows; ++i)
        {
            pixel = pixelColumn[i];
            Bool isInsideTrailLength = nTransfersFromTrap < ctePars->cte_len;

            //check if this are needed
            chargeToAdd = 0;
            extraChargeToAdd = 0;
            chargeToRemove = 0;

            /*HAPPENS AFTER FIRST PASS*/
            /*SHUFFLE CHARGE IN*/
            //move out of loop to separate instance?
            if (i > 0)
            {
                if ((double)traps[i] < (double)traps[i-1])
                    trappedFlux *= ((double)traps[i] / (double)traps[i-1]);
            }

            if (pixel >= ctePars->qlevq_data[w])
            {
                if (isInsideTrailLength)
                {
                    ++nTransfersFromTrap;
                    chargeToAdd = rprof->data[w*rprof->ny + nTransfersFromTrap-1] * trappedFlux;
                    extraChargeToAdd = cprof->data[w*cprof->ny + nTransfersFromTrap-1] * trappedFlux;
                }
                trappedFlux = ctePars->dpdew_data[w] / ctePars->n_par * (double)traps[i];
                chargeToRemove = trappedFlux;
                nTransfersFromTrap = 0;
            }
            else
            {
                if (isInsideTrailLength)
                {
                    ++nTransfersFromTrap;
                    chargeToAdd = rprof->data[w*rprof->ny + nTransfersFromTrap-1] * trappedFlux;
                }
            }

            pixelColumn[i] += chargeToAdd + extraChargeToAdd - chargeToRemove;
        }} //end for i
    }} //end for w
    return HSTCAL_OK;
}

int simulateColumnReadout(double * const pixelColumn, const float * const traps, const CTEParamsFast * const cte,
        const FloatTwoDArray * const rprof, const FloatTwoDArray * const cprof, const unsigned nRows, const unsigned nPixelShifts)
{
    //For performance this does not NULL check passed in ptrs

    int localStatus = HSTCAL_OK;
    //Take each pixel down the detector
    {unsigned shift;
    for (shift = 1; shift <= nPixelShifts; ++shift)
    {
        if ((localStatus = simulatePixelReadout_v1_2(pixelColumn, traps, cte, rprof, cprof, nRows)))
            return localStatus;
    }}

    return localStatus;
}

Bool correctCROverSubtraction(float * const traps, const double * const pix_model, const double * const pix_observed,
        const unsigned nRows, const double threshHold)
{
    /*LOOK FOR AND DOWNSCALE THE CTE MODEL IF WE FIND
      THE TELL-TALE SIGN OF READOUT CRS BEING OVERSUBTRACTED;
      IF WE FIND ANY THEN GO BACK UP AND RERUN THIS COLUMN

      THIS MODEL SEARCHES FOR OVERSUBTRACTED TRAILS.
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

    //For performance this does not NULL check passed in ptrs

    Bool redo = False;
    {unsigned i;
    for (i = 10; i < nRows-2; ++i)
    {
        if ( (( threshHold > pix_model[i] ) &&
                    ( threshHold > (pix_model[i] - pix_observed[i]))) ||

                (((pix_model[i] + pix_model[i+1]) < -12.) &&
                 (pix_model[i] + pix_model[i+1] - pix_observed[i] - pix_observed[i+1] < -12.)) ||

                (((pix_model[i] + pix_model[i+1] + pix_model[i+2]) < -15.) &&
                 ((pix_model[i] + pix_model[i+1] + pix_model[i+2] - pix_observed[i] -
                   pix_observed[i+1] - pix_observed[i+2]) < -15.))  )
        {
            redo = True;
            unsigned iMax = i;

            /*GO DOWNSTREAM AND LOOK FOR THE OFFENDING CR*/
            {unsigned ii;
            for (ii = i-10; ii <= i; ++ii)
            {
                if ( (pix_model[ii] - pix_observed[ii]) > (pix_model[iMax] - pix_observed[iMax]) )
                    iMax = ii;
            }}
            /* DOWNGRADE THE CR'S SCALING AND ALSO FOR THOSE
               BETWEEN THE OVERSUBTRACTED PIXEL AND IT*/
            {unsigned ii;
            for (ii = iMax; ii <= i; ++ii)
                traps[ii] *= 0.75;
            }
        }
    }} /*end for  j*/
    return redo;
}

int populateTrapPixelMap(SingleGroup * trapPixelMap, CTEParamsFast * ctePars)
{
    /*
        int   iz_data[<cte->nScaleTableColumns>];  column number in raz format
        double scale512[<cte->nScaleTableColumns>];      scaling appropriate at row 512
        double scale1024[<cte->nScaleTableColumns>];     scaling appropriate at row 1024
        double scale1536[<cte->nScaleTableColumns>];     scaling appropriate at row 1536
        double scale2048[<cte->nScaleTableColumns>];     scaling appropriate at row 2048
    */

    //For performance this does not NULL check passed in ptrs

    //WARNING - OUTPUTS column major storage order
    trapPixelMap->sci.data.storageOrder = COLUMNMAJOR;

    clock_t begin = clock();

    extern int status;

    const unsigned nRows = trapPixelMap->sci.data.ny;
    const unsigned nColumns = trapPixelMap->sci.data.nx;
    const double cteScale = ctePars->scale_frac;


#ifdef _OPENMP
    #pragma omp parallel shared(trapPixelMap, ctePars)
#endif
    {
    double trapColumnScale[4];
    double cte_i;
    double cte_j;
    double ro;
    int io;

    {unsigned i;
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (i = 0; i < ctePars->nScaleTableColumns; ++i)
    {
        unsigned column = ctePars->iz_data[i] - ctePars->razColumnOffset; //which column to scale
        if (column < 0 || column >= nColumns)//vec blocker
            continue;
        trapColumnScale[0] = ctePars->scale512[i];
        trapColumnScale[1] = ctePars->scale1024[i];
        trapColumnScale[2] = ctePars->scale1536[i];
        trapColumnScale[3] = ctePars->scale2048[i];
        //CALCULATE THE CTE CORRECTION FOR EVERY PIXEL
        //  Index is figured on the final size of the image
        //  not the current size.
        {unsigned  j;
        for (j = 0; j < nRows; ++j)
        {
            ro = j / 512.0; //ro can be zero, it's an index
            if (ro > 2.999)
                ro = 2.999; // only 4 quads, 0 to 3
            else if (ro < 0)
                ro = 0;
            io = (int) floor(ro); //force truncation towards 0 for pos numbers
            cte_j = (j+1) / 2048.0;
            cte_i = trapColumnScale[io] + (trapColumnScale[io+1] - trapColumnScale[io]) * (ro - io);
            PixColumnMajor(trapPixelMap->sci.data, j, column) = cte_i * cte_j * cteScale;
        }}
    }}
    } // end parallel block

    if (ctePars->verbose)
    {
        double timeSpent = ((double)(clock() - begin))/CLOCKS_PER_SEC;
        sprintf(MsgText,"(pctecorr) Time taken to populate pixel trap map image: %.2f(s) with %i threads",timeSpent/ctePars->maxThreads, ctePars->maxThreads);
        trlmessage(MsgText);
    }

    return(status);
}

int cteSmoothImage(const SingleGroup * input, SingleGroup * output, CTEParamsFast * ctePars, double ampReadNoise)
{
    /*
       This routine will output the smoothest
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
    if (!input || !output || !ctePars)
        return (status = ALLOCATION_PROBLEM);

    //WARNING - assumes column major storage order
    assert(input->sci.data.storageOrder == COLUMNMAJOR);
    output->sci.data.storageOrder = COLUMNMAJOR;

    extern int status;

    const unsigned nRows = input->sci.data.ny;
    const unsigned nColumns = input->sci.data.nx;

    double rmsGlobal=0;
    double nrmsGlobal=0;

    clock_t begin = clock();

    copySingleGroup(output, input, input->sci.data.storageOrder);

    //Is the readnoise diff per amp? Current method assumes not.
    if (ampReadNoise < 0.1){
        trlmessage("rnsig < 0.1, No read-noise mitigation needed");
        return(status);
    }

    /*GO THROUGH THE ENTIRE IMAGE AND ADJUST PIXELS TO MAKE THEM
      SMOOTHER, BUT NOT SO MUCH THAT IT IS NOT CONSISTENT WITH
      READNOISE.  DO THIS IN BABY STEPS SO THAT EACH ITERATION
      DOES VERY LITTLE ADJUSTMENT AND INFORMATION CAN GET PROPAGATED
      DOWN THE LINE.
      */

    //To remove the below code adjust in place i.e. using only output:
    //Don't use pointers to output for obs_loc & rsz_loc
    //Copy columns and then just shift these copies (boundary case might be annoying)
    //Use schedule(static) and pre (inner for loop) copy boundary columns to avoid race conditions
    SingleGroup adjustment;
    initSingleGroup(&adjustment);
    allocSingleGroup(&adjustment, nColumns, nRows, False);

    SingleGroup readNoise;
    initSingleGroup(&readNoise);
    allocSingleGroup(&readNoise, nColumns, nRows, False);

#ifdef _OPENMP
    #pragma omp parallel shared(input, output, ampReadNoise, rmsGlobal, nrmsGlobal, readNoise)
#endif
    {
    const float * obs_loc[3];
    const float * rsz_loc[3];

    double rmsLocal;
    double nrmsLocal;
    {unsigned iter;
    for (iter = 0; iter < 100; ++iter)
    {
        rmsLocal = 0;
        nrmsLocal = 0;
        {unsigned i;
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (i = 0; i < nColumns; ++i)
        {
            unsigned imid = i;
            /*RESET TO MIDDLE nColumns AT ENDPOINTS*/
            // This seems odd, the edge columns get accounted for twice?
            if (i == 0)
                imid = 1;
            else if (i == nColumns-1) // NOTE: use of elseif breaks if nColumns = 1
                imid = nColumns-2;

            /*LOCATE THE MIDDLE AND NEIGHBORING PIXELS FOR ANALYSIS*/
            obs_loc[0] = input->sci.data.data + (imid-1)*nRows;
            obs_loc[1] = obs_loc[0] + nRows;
            obs_loc[2] = obs_loc[1] + nRows;

            rsz_loc[0] = output->sci.data.data + (imid-1)*nRows;
            rsz_loc[1] = rsz_loc[0] + nRows;
            rsz_loc[2] = rsz_loc[1] + nRows;

            {unsigned j;
            for (j = 0; j < nRows; ++j)
                PixColumnMajor(adjustment.sci.data, j, i) = find_dadjFast(1+i-imid, j, nRows, obs_loc, rsz_loc, ampReadNoise);
            }
        }} /*end the parallel for*/ //implicit omp barrier

        //NOW GO OVER ALL THE nColumns AND nRows AGAIN TO SCALE THE PIXELS
        {unsigned i;
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (i = 0; i < nColumns; ++i)
        {
            {unsigned j;
            for(j = 0; j < nRows; ++j)
            {
                PixColumnMajor(output->sci.data,j,i) += (PixColumnMajor(adjustment.sci.data,j, i)*0.75);
                PixColumnMajor(readNoise.sci.data,j,i) = (PixColumnMajor(input->sci.data,j,i) - PixColumnMajor(output->sci.data,j,i));
            }}
        }}//implicit omp barrier

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            rmsGlobal=0;
            nrmsGlobal=0;
        }

        {unsigned j;
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (j = 0; j < nColumns; ++j)
        {
            {unsigned i;
            for (i = 0; i < nRows; ++i)
            {
                if ( (fabs(PixColumnMajor(input->sci.data, i, j)) > 0.1 ||
                     fabs(PixColumnMajor(output->sci.data, i, j)) > 0.1))
                {
                    double tmp = PixColumnMajor(readNoise.sci.data, i, j);
                    rmsLocal  +=  tmp*tmp;
                    ++nrmsLocal;
                }
            }}
        }}//implicit omp barrier

#ifdef _OPENMP
        #pragma omp critical (aggregate)
#endif
        {
            rmsGlobal  += rmsLocal;
            nrmsGlobal += nrmsLocal;
        }
#ifdef _OPENMP
        #pragma omp barrier
#endif

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            rmsGlobal = sqrt(rmsGlobal/nrmsGlobal);
        } //implicit barrier

        // if it is true that one breaks then it is true for all
        /*epsilon type comparison*/
        if ((ampReadNoise - rmsGlobal) < 0.00001)
            break; // this exits loop over iter

#ifdef _OPENMP
        #pragma omp barrier
#endif
    }} // end loop over iter
    } // close parallel block
    freeSingleGroup(&adjustment);
    freeSingleGroup(&readNoise);

    if (ctePars->verbose)
    {
        double timeSpent = ((double)(clock() - begin))/CLOCKS_PER_SEC;
        sprintf(MsgText,"(pctecorr) Time taken to smooth image: %.2f(s) with %i threads", timeSpent/ctePars->maxThreads, ctePars->maxThreads);
        trlmessage(MsgText);
    }

    return (status);
}

double find_dadjFast(const unsigned i, const unsigned j, const unsigned nRows, const float * obsloc[3], const float * rszloc[3], const double readNoiseAmp)
{
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

    //For performance this does not NULL check passed in ptrs

    const double mval = (double)*(rszloc[i] + j);
    const double dval0  = (double)*(obsloc[i] + j) - mval;
    double dval0u = dval0;

    if (dval0u > 1)
        dval0u = 1;
    else if (dval0u < -1)
        dval0u = -1;

    /*COMPARE THE SURROUNDING PIXELS*/
    double dval9 = 0.0;
    if (i == 1 &&  j < nRows-1 && j > 0)
    {
        dval9 = (double)*(obsloc[i]   + j-1)   - (double)*(rszloc[i]   + j-1) +
                (double)*(obsloc[i]   + j)     - (double)*(rszloc[i]   + j)   +
                (double)*(obsloc[i]   + j+1)   - (double)*(rszloc[i]   + j+1) +
                (double)*(obsloc[i-1] + j-1)   - (double)*(rszloc[i-1] + j-1) +
                (double)*(obsloc[i-1] + j)     - (double)*(rszloc[i-1] + j)   +
                (double)*(obsloc[i-1] + j+1)   - (double)*(rszloc[i-1] + j+1) +
                (double)*(obsloc[i+1] + j-1)   - (double)*(rszloc[i+1] + j-1) +
                (double)*(obsloc[i+1] + j)     - (double)*(rszloc[i+1] + j)  +
                (double)*(obsloc[i+1] + j+1)   - (double)*(rszloc[i+1] + j+1);

        dval9 = dval9 / 9.0;
    }
    const double readNoiseAmpFraction = 0.33;
    double dval9u = dval9;
    if (dval9u > readNoiseAmp*readNoiseAmpFraction)
        dval9u = readNoiseAmp*readNoiseAmpFraction;
    else if (dval9u < readNoiseAmp*-readNoiseAmpFraction)
        dval9u = readNoiseAmp*-readNoiseAmpFraction;

    const double dmod1 = j > 0 ? (double)*(rszloc[i] + j-1) - mval : 0;
    const double dmod2 = j < nRows-1 ? (double)*(rszloc[i] + j+1) - mval : 0;

    double dmod1u = dmod1;
    if (dmod1u > readNoiseAmp*readNoiseAmpFraction)
        dmod1u = readNoiseAmp*readNoiseAmpFraction;
    else if (dmod1u < readNoiseAmp*-readNoiseAmpFraction)
        dmod1u = readNoiseAmp*-readNoiseAmpFraction;

    double dmod2u = dmod2;
    if (dmod2u > readNoiseAmp*readNoiseAmpFraction)
        dmod2u = readNoiseAmp*readNoiseAmpFraction;
    else if (dmod2u < readNoiseAmp*-readNoiseAmpFraction)
        dmod2u = readNoiseAmp*-readNoiseAmpFraction;

    /*
       IF IT'S WITHIN 2 SIGMA OF THE READNOISE, THEN
       TEND TO TREAT AS READNOISE; IF IT'S FARTHER OFF
       THAN THAT, THEN DOWNWEIGHT THE INFLUENCE
       */
    const double readNoiseAmp2 = readNoiseAmp*readNoiseAmp;
    const double w0 =     dval0 * dval0 / (dval0 * dval0 + 4.0 * readNoiseAmp2);
    const double w9 =     dval9 * dval9 / (dval9 * dval9 + 18.0 * readNoiseAmp2);
    const double w1 = 4 * readNoiseAmp2 / (dmod1 * dmod1 + 4.0 * readNoiseAmp2);
    const double w2 = 4 * readNoiseAmp2 / (dmod2 * dmod2 + 4.0 * readNoiseAmp2);

    /*(note that with the last two, if a pixel
      is too discordant with its upper or lower
      that neighbor has less of an ability to
      pull it)*/

    return  dval0u * w0 * 0.25f + /* desire to keep the original pixel value */
            dval9u * w9 * 0.25f + /* desire to keep the original sum over 3x3*/
            dmod1u * w1 * 0.25f + /* desire to get closer to the pixel below*/
            dmod2u * w2 * 0.25f; /* desire to get closer to the pixel above*/
}
