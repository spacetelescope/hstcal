#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#include "hstcal_memory.h"
#include "hstcal.h"
#include "hstio.h"

#include "acs.h"
#include "acsinfo.h"
#include "hstcalerr.h"
#include "trlbuf.h"

#include "pcte.h"

# ifdef _OPENMP
#  include <omp.h>
# endif
# include "../../../../ctegen2/ctegen2.h"
#include <assert.h>

int get_amp_array_size_acs_cte(const ACSInfo *acs, SingleGroup *amp,
                               char *amploc, char *ccdamp,
                               int *xsize, int *ysize, int *xbeg,
                               int *xend, int *ybeg, int *yend);

static int extractAmp(SingleGroup * amp, const SingleGroup * image, const unsigned ampID, CTEParamsFast * ctePars);
static int insertAmp(SingleGroup * amp, const SingleGroup * image, const unsigned ampID, CTEParamsFast * ctePars);
static int alignAmpData(FloatTwoDArray * amp, const unsigned ampID);
static int alignAmp(SingleGroup * amp, const unsigned ampID);
static int rotateAmp(SingleGroup * amp, const unsigned ampID, bool derotate, char ccdamp);
static int rotateAmpData(FloatTwoDArray * amp, const unsigned ampID);
static int derotateAmpData(FloatTwoDArray * amp, const unsigned ampID);

static void transpose(FloatTwoDArray *amp);
static void side2sideFlip(FloatTwoDArray *amp);
static void top2bottomFlip(FloatTwoDArray *amp);

int doPCTEGen3 (ACSInfo *acs, CTEParamsFast * ctePars, SingleGroup * chipImage, const bool forwardModelOnly, char * corrType, char *ccdamp, int nthAmp, char *amploc, int ampID)

{

    /* arguments:
       acs     i: calibration switches, etc
       x      io: image to be calibrated; written to in-place
    */

    extern int status;

    if (!acs || !ctePars || !chipImage)
        return (status = ALLOCATION_PROBLEM);

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
#ifdef _OPENMP
    if (acs->nThreads > 1)
        trlmessage("(pctecorr) Using multiprocessing provided by OpenMP inside CTE routine.");
#endif

    int nRows, nColumns;    /* int for C array dimension sizes */
    int amp_xsize, amp_ysize;  /* int for amp array size (x/y in CCD coords) */
    int amp_xbeg, amp_xend;    /* int for beg and end of amp arrays on chip */
    int amp_ybeg, amp_yend;

    /* functions from calacs/lib */
    void TimeStamp (char *message, char *rootname);
    int PutKeyDbl (Hdr *hd, char *keyword, double value, char *comment);

    if (acs->printtime) {
        TimeStamp("Starting CTE correction...","");
    }

    sprintf(MsgText, "(pctecorr) Performing CTE correction type %s for amp %c.", corrType, ccdamp[nthAmp]);
    trlmessage(MsgText);

    ctePars->rn_amp = acs->readnoise[ampID];
    sprintf(MsgText, "(pctecorr) Read noise level from CCDTAB: %f.", ctePars->rn_amp);
    trlmessage(MsgText);

    /* get amp array size */
    if (get_amp_array_size_acs_cte(acs, chipImage, amploc, ccdamp,
                                   &amp_xsize, &amp_ysize, &amp_xbeg,
                                   &amp_xend, &amp_ybeg, &amp_yend))
    {
        freeOnExit(&ptrReg);
        return (status);
    }

    nRows = amp_ysize;
    nColumns = amp_xsize;
    ctePars->nRows = nRows;
    ctePars->nColumns = nColumns;
    ctePars->columnOffset = amp_xbeg;
    ctePars->rowOffset = amp_ybeg;
    // razColumnOffset is used to align the image with the column scaling (SCLBYCOL) in the PCTETAB.
    // The SCLBYCOL array is 8192 cols wide which is 2048*4 and therefore does NOT include the physical
    // overscan columns (24 for ACS WFC) like this array does for calwf3.
    // The raz format is all 4 amps side by side in CDAB order.
    ctePars->razColumnOffset = ctePars->columnOffset;
    if (ampID == AMP_A || ampID == AMP_B)
        ctePars->razColumnOffset += ACS_WFC_N_COLUMNS_PER_CHIP_EXCL_OVERSCAN;

    //This is used for the final output
    SingleGroup ampImage;
    initSingleGroup(&ampImage);
    addPtr(&ptrReg, &ampImage, &freeSingleGroup);
    if (allocSingleGroupExts(&ampImage, nColumns, nRows, SCIEXT | ERREXT, False) != 0)
    {
        freeOnExit(&ptrReg);
        return (status = OUT_OF_MEMORY);
    }

    // read data from the SingleGroup into an array containing data from
    // just one amp
    if ((status = extractAmp(&ampImage, chipImage, ampID, ctePars)))
    {
        freeOnExit(&ptrReg);
        return (status);
    }

    // rotate amp in preparation for the serial CTE correction
    // This step is done ONLY to accommodate the serial CTE correction
    // which preceeds the standard parallel CTE correction.
    if (strcmp(corrType, "serial") == 0)
    {
        sprintf(MsgText, "\n(pctecorr) Invoking rotation for CTE Correction Type: %s", corrType);
        trlmessage(MsgText);
        if ((status = rotateAmp(&ampImage, ampID, False, ccdamp[nthAmp])))
        {
            freeOnExit(&ptrReg);
            return (status);
        }
    }

    if ((status = alignAmp(&ampImage, ampID)))
    {
        freeOnExit(&ptrReg);
        return (status);
    }

    //copy to column major storage
    SingleGroup columnMajorImage;
    initSingleGroup(&columnMajorImage);
    addPtr(&ptrReg, &columnMajorImage, &freeSingleGroup);
    if (allocSingleGroupExts(&columnMajorImage, nColumns, nRows, SCIEXT, False) != 0)
    {
        freeOnExit(&ptrReg);
        return (status = OUT_OF_MEMORY);
    }
    if ((status = copySingleGroup(&columnMajorImage, &ampImage, COLUMNMAJOR)))
    {
        freeOnExit(&ptrReg);
        return (status = ALLOCATION_PROBLEM);
    }

    if (forwardModelOnly)
    {
        // Not needed as we won't be smoothing.
    }
    else
    {
        //CALCULATE THE SMOOTH READNOISE IMAGE
        trlmessage("(pctecorr) Calculating smooth readnoise image...");
        if (ctePars->noise_mit != 0)
        {
            trlerror("(pctecorr) Only noise model 0 implemented!");
            freeOnExit(&ptrReg);
            return (status = ERROR_RETURN);
        }
    }
    SingleGroup smoothedImage;
    initSingleGroup(&smoothedImage);
    addPtr(&ptrReg, &smoothedImage, &freeSingleGroup);
    if (allocSingleGroupExts(&smoothedImage, nColumns, nRows, SCIEXT, False) != 0)
    {
        freeOnExit(&ptrReg);
        return (status = OUT_OF_MEMORY);
    }
    setStorageOrder(&smoothedImage, COLUMNMAJOR);

    if (forwardModelOnly)
    {
        // Don't actually smooth but copy to smoothedImage
        if ((status = copySingleGroup(&smoothedImage, &columnMajorImage, COLUMNMAJOR)))
        {
            freeOnExit(&ptrReg);
            return (status = ALLOCATION_PROBLEM);
        }
    }
    else
    {
        // do some smoothing on the data so we don't amplify the read noise.
        if ((status = cteSmoothImage(&columnMajorImage, &smoothedImage, ctePars, ctePars->rn_amp)))
        {
            freeOnExit(&ptrReg);
            return (status);
        }
    }
    trlmessage("(pctecorr) ...complete.");

    trlmessage("(pctecorr) Creating charge trap image...");
    SingleGroup trapPixelMap;
    initSingleGroup(&trapPixelMap);
    addPtr(&ptrReg, &trapPixelMap, &freeSingleGroup);
    if (allocSingleGroupExts(&trapPixelMap, nColumns, nRows, SCIEXT, False) != 0)
    {
        freeOnExit(&ptrReg);
        return (status = OUT_OF_MEMORY);
    }
    setStorageOrder(&trapPixelMap, COLUMNMAJOR);
    if ((status = populateTrapPixelMap(&trapPixelMap, ctePars)))
    {
        freeOnExit(&ptrReg);
        return status;
    }
    trlmessage("(pctecorr) ...complete.");

    SingleGroup * cteCorrectedImage = &columnMajorImage;
    if (forwardModelOnly)
    {
        trlmessage("(pctecorr) Running forward model simulation...");
        //perform CTE correction
        if ((status = forwardModel(&smoothedImage, cteCorrectedImage, &trapPixelMap, ctePars)))
        {
            freeOnExit(&ptrReg);
            return status;
        }
    }
    else
    {
        trlmessage("(pctecorr) Running correction algorithm...");
        //perform CTE correction
        if ((status = inverseCTEBlur(&smoothedImage, cteCorrectedImage, &trapPixelMap, ctePars)))
        {
            freeOnExit(&ptrReg);
            return status;
        }
    }
    freePtr(&ptrReg, &trapPixelMap);
    trlmessage("(pctecorr) ...complete.");

    // add 10% correction to error in quadrature.
    float totalCounts = 0;
    float totalRawCounts = 0;
#ifdef _OPENMP
    #pragma omp parallel shared(nRows, nColumns, cteCorrectedImage, smoothedImage, ampImage)
#endif
    {
        double temp_err;
        float threadCounts = 0;
        float threadRawCounts = 0;
        float delta;
        {   unsigned k;
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            //This loop order is row major as more row major storage is accessed than column.
            //MORE: look into worth splitting out ops, prob needs a order swap (copy) so perhaps not worth it.
            for (k = 0; k < nRows; ++k)
            {
                {   unsigned m;
                    for (m = 0; m < nColumns; ++m)
                    {
                        delta = (PixColumnMajor(cteCorrectedImage->sci.data,k,m) - PixColumnMajor(smoothedImage.sci.data,k,m));
                        threadCounts += delta;
                        threadRawCounts += Pix(ampImage.sci.data, m, k);
                        temp_err = 0.1 * fabs(delta);
                        Pix(ampImage.sci.data, m, k) += delta;

                        float err2 = Pix(ampImage.err.data, m, k);
                        err2 *= err2;
                        Pix(ampImage.err.data, m, k) = sqrt(err2 + temp_err*temp_err);
                    }
                }
            }
        }

#ifdef _OPENMP
        #pragma omp critical(deltaAggregate)
#endif
        {
            totalCounts += threadCounts;
            totalRawCounts += threadRawCounts;
        }
    }//close parallel block
    sprintf(MsgText, "(pctecorr) Total count difference (corrected-raw) incurred from correction: %f (%f%%)", totalCounts, totalCounts/totalRawCounts*100);
    trlmessage(MsgText);

    // put the CTE corrected data back into the SingleGroup structure
    if ((status = alignAmp(&ampImage, ampID)))
    {
        freeOnExit(&ptrReg);
        return (status);
    }

    // Derotate amp so the parallel CTE correction can proceed
    if (strcmp(corrType, "serial") == 0)
    {
        sprintf(MsgText, "\n(pctecorr) Invoking de-rotation for CTE Correction Type: %s", corrType);
        trlmessage(MsgText);
        if ((status = rotateAmp(&ampImage, ampID, True, ccdamp[nthAmp])))
        {
            freeOnExit(&ptrReg);
            return (status);
        }
    }

    if ((status = insertAmp(chipImage, &ampImage, ampID, ctePars)))
    {
        freeOnExit(&ptrReg);
        return (status);
    }

    // free mem used by our amp arrays
    freePtr(&ptrReg, &ampImage);
    freePtr(&ptrReg, &columnMajorImage);
    freePtr(&ptrReg, &smoothedImage);
    if (acs->printtime)
        TimeStamp("CTE corrections complete","");

    freeOnExit(&ptrReg);
    return (status);
}

static int extractAmp(SingleGroup * amp,  const SingleGroup * image, const unsigned ampID, CTEParamsFast * ctePars)
{
    extern int status;

    if (!amp || !amp->sci.data.data || !image || !image->sci.data.data)
        return (status = ALLOCATION_PROBLEM);

    //WARNING - assumes row major storage
    assert(amp->sci.data.storageOrder == ROWMAJOR && image->sci.data.storageOrder == ROWMAJOR);

    if (ampID != AMP_A && ampID != AMP_B && ampID != AMP_C && ampID != AMP_D)
    {
        trlerror("Amp number not recognized, must be 0-3.");
        return (status = ERROR_RETURN);
    }

    const unsigned nRows = amp->sci.data.ny;
    const unsigned nColumns = amp->sci.data.nx;

    unsigned rowSkipLength = image->sci.data.nx;
    unsigned offset = ctePars->columnOffset;

    copyOffsetSingleGroup(amp, image, nRows, nColumns, 0, offset, nColumns, rowSkipLength);
    return status;
}

/*
   The rotation routines support the serial CTE correction such that the data
   is configured to mimic the parallel CTE data.  The idea is to 
   rotate the amp to put the serial trails in the same orientation as 
   those of the parallel trails and then the existing CTE functions are applied.
   Essentially, the CTE correction routines can then be applied in the identical
   manner for BOTH the serial and the parallel CTE correction cases. 

   Once the serial CTE correction has been performed, the amps then need to be
   "de-rotated" so the parallel CTE correction can then be applied.
*/
static int rotateAmp(SingleGroup * amp, const unsigned ampID, bool derotate, char ccdamp)
{
    extern int status;

    if (!amp)
        return (status = ALLOCATION_PROBLEM);

    // The case of rotating the amp into position to mimic parallel processing
    if (!derotate)
    {
        //sci data
        if (amp->sci.data.data)
        {
            sprintf(MsgText, "(pctecorr) Rotation for Amp: %c\n", ccdamp);
            trlmessage(MsgText);
            if ((status = rotateAmpData(&amp->sci.data, ampID)))
                return status;
        }

        //err data
        if (amp->err.data.data)
        {
            if ((status = rotateAmpData(&amp->err.data, ampID)))
                return status;
        }
    }
    // The case of returning the amp to its original orientation
    else
    {
        //sci data
        if (amp->sci.data.data)
        {
            sprintf(MsgText, "(pctecorr) Derotation for Amp: %c\n", ccdamp);
            trlmessage(MsgText);
            if ((status = derotateAmpData(&amp->sci.data, ampID)))
                return status;
        }

        //err data
        if (amp->err.data.data)
        {
            if ((status = derotateAmpData(&amp->err.data, ampID)))
                return status;
        }
    }

    return status;
}

static int rotateAmpData(FloatTwoDArray * amp, const unsigned ampID)
{
    extern int status;
    if (!amp || !amp->data)
        return (status = ALLOCATION_PROBLEM);

    /*
       Rotate the amp to put the serial trails in the same orientation
       as the parallel trails would be. A rotation requires a transpose
       and then a flip. Always transpose the data first.
       
    */
    transpose(amp);

    /*
      To complete the correct rotation, the flip either has to be
      from side-to-side about the central column or top-to-bottom
      about the central row.
    */
    if (ampID == AMP_B || ampID == AMP_C) {
        side2sideFlip(amp);
    } else {
        top2bottomFlip(amp);
    }

    return status;
}

static int derotateAmpData(FloatTwoDArray * amp, const unsigned ampID)
{
    extern int status;
    if (!amp || !amp->data)
        return (status = ALLOCATION_PROBLEM);

    /*
       To derotate the amp, you are reversing the direction of the initial
       rotation.  Always transpose the data first, and then apply the 
       side-to-side or top-to-bottom flip to put the amp back into its 
       original orientation so the parallel CTE correction can proceed.
    */
    transpose(amp);

    /*
      The rotation here is in the opposite direction from the initial
      rotation for the amp in question, so if a side-to-side flip were
      done initially, now do a top-to-bottom flip (for example).
    */
    if (ampID == AMP_B || ampID == AMP_C) {
        top2bottomFlip(amp);
    } else {
        side2sideFlip(amp);
    }

    return status;
}

/*
   Transpose the amp data
*/
static void transpose(FloatTwoDArray * amp)
{
    const unsigned nRows = amp->ny;
    const unsigned nColumns = amp->nx;
    float temp;

    for (unsigned i = 0; i < nRows; i++) {
        for (unsigned j = i; j < nColumns; j++) {
            temp = PPix(amp, j, i);
            PPix(amp, j, i) = PPix(amp, i, j);
            PPix(amp, i, j) = temp;
        }
    }
}

/*
   Flip the amp about the Y-axis central column (i.e., flip from left to right)
*/
static void side2sideFlip(FloatTwoDArray * amp)
{
    const unsigned nRows = amp->ny;
    const unsigned nColumns = amp->nx;
    float temp;

    for (unsigned i = 0; i < nRows; i++) {
        for (unsigned j = 0; j < nColumns/2; j++) {
            temp = PPix(amp, j, i);
            PPix(amp, j, i) = PPix(amp, nColumns - j - 1, i);
            PPix(amp, nColumns - j - 1, i) = temp;
        }
    }
}

static int insertAmp(SingleGroup * image, const SingleGroup * amp, const unsigned ampID, CTEParamsFast * ctePars)
{
    extern int status;

    if (!amp || !amp->sci.data.data || !image || !image->sci.data.data)
        return (status = ALLOCATION_PROBLEM);

    //WARNING - assumes row major storage
    assert(amp->sci.data.storageOrder == ROWMAJOR && image->sci.data.storageOrder == ROWMAJOR);

    if (ampID != AMP_A && ampID != AMP_B && ampID != AMP_C && ampID != AMP_D)
    {
        trlerror("Amp number not recognized, must be 0-3.");
        return (status = ERROR_RETURN);
    }

    const unsigned nRows = amp->sci.data.ny;
    const unsigned nColumns = amp->sci.data.nx;

    unsigned rowSkipLength = image->sci.data.nx;
    unsigned offset = ctePars->columnOffset;

    copyOffsetSingleGroup(image, amp, nRows, nColumns, offset, 0, rowSkipLength, nColumns);
    return status;
}

static int alignAmp(SingleGroup * amp, const unsigned ampID)
{
    extern int status;
    if (!amp)
        return (status = ALLOCATION_PROBLEM);

    //sci data
    if (amp->sci.data.data)
    {
        if ((status = alignAmpData(&amp->sci.data, ampID)))
            return status;
    }
    //err data
    if (amp->err.data.data)
    {
        if ((status = alignAmpData(&amp->err.data, ampID)))
            return status;
    }
    //dq data
    if (amp->dq.data.data)
        assert(0);//unimplemented

    return status;
}

static int alignAmpData(FloatTwoDArray * amp, const unsigned ampID)
{
    //NOTE: There is a similar version of this in wfc3 - code changes should be reflected in both.

    //Align image quadrants such that the amps are at the bottom left, i.e. aligned with amp C.
    extern int status;

    if (!amp || !amp->data)
        return (status = ALLOCATION_PROBLEM);

    //Amp C is already correctly aligned
    if (ampID == AMP_C)
        return status;

    //WARNING - assumes row major storage
    assert(amp->storageOrder == ROWMAJOR);

    //Flip about y axis, i.e. about central column
    if (ampID == AMP_B || ampID == AMP_D)
    {
        side2sideFlip(amp);
    }

    //Only thing left is to flip AB chip upside down
    //Flip about x axis, i.e. central row
    if (ampID == AMP_A || ampID == AMP_B)
    {
        top2bottomFlip(amp);
    }

    return status;
}

/*
   Flip the amp about the X-axis central row (i.e., flip from top to bottom)

   Note: This routine was originally in alignAmpData flow of processing.  It 
   is now encapulated here as it is used multiple times due to the rotations
   needed to accommodate the serial CTE correction.  The original OPENMP 
   specifications have been left intact.
*/
static void top2bottomFlip(FloatTwoDArray * amp)
{

    const unsigned nColumns = amp->nx;
    const unsigned nRows = amp->ny;

    Bool allocationFail = False;
#ifdef _OPENMP
    #pragma omp parallel shared(amp, allocationFail)
#endif
    {
        float * tempRow = NULL;
        size_t rowSize = nColumns*sizeof(*tempRow);
        tempRow = malloc(rowSize);
        if (!tempRow)
        {
#ifdef _OPENMP
            #pragma omp critical(allocationCheck) // This really isn't needed
#endif
            {
                allocationFail = True;
            }
        }

#ifdef _OPENMP
        #pragma omp barrier
#endif
        if (!allocationFail)
        {
            {   unsigned i;
#ifdef _OPENMP
                #pragma omp for schedule(static)
#endif
                for (i = 0; i < nRows/2; ++i)
                {
                    float * topRow = amp->data + i*nColumns;
                    float * bottomRow = amp->data + (nRows-1-i)*nColumns;
                    memcpy(tempRow, topRow, rowSize);
                    memcpy(topRow, bottomRow, rowSize);
                    memcpy(bottomRow, tempRow, rowSize);
                }
            }
        }
        if (tempRow)
            free(tempRow);
    }//close parallel block
    if (allocationFail)
    {
        sprintf(MsgText, "Out of memory for 'tempRow' in 'alignAmpData'");
        trlerror(MsgText);
    }
}
