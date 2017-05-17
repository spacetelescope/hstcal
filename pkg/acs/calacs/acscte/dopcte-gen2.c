#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "hstio.h"

#include "acs.h"
#include "acsinfo.h"
#include "hstcalerr.h"

#include "pcte.h"

# ifdef _OPENMP
#  include <omp.h>
# endif
# include "../../../../ctegen2/ctegen2.h"
#include <assert.h>

int get_amp_array_size_acs_cte(const ACSInfo *acs, SingleGroup *amp,
                              const int ampID, char *amploc, char *ccdamp,
                              int *xsize, int *ysize, int *xbeg,
                              int *xend, int *ybeg, int *yend);

static int extractAmp(SingleGroup * amp, const SingleGroup * image, const unsigned ampID);
static int insertAmp(SingleGroup * amp, const SingleGroup * image, const unsigned ampID);
static int alignAmpData(FloatTwoDArray * amp, const unsigned ampID);
static int alignAmp(SingleGroup * amp, const unsigned ampID);

int doPCTEGen2 (ACSInfo *acs, CTEParamsFast * ctePars, SingleGroup * chipImage)
{

    /* arguments:
       acs     i: calibration switches, etc
       x      io: image to be calibrated; written to in-place
    */

    extern int status;

    if (!acs || !ctePars || !chipImage)
        return (status = ALLOCATION_PROBLEM);

    if (acs->subarray != 0)
    {
        sprintf(MsgText, "UNIMPLEMENTED - CTE correction for fullframe images only");
        trlerror(MsgText);
        return (status = INTERNAL_ERROR);
    }

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
#ifdef _OPENMP
    if (acs->nThreads > 1)
        trlmessage("(pctecorr) Using parallel processing provided by OpenMP inside CTE routine");
#endif

    char ccdamp[strlen(AMPSTR1)+1]; /* string to hold amps on current chip */
    int numamps;               /* number of amps on chip */
    int ampID;                   /* index amp A:0, B:1, etc. */
    char * amploc;             /* pointer to amp character in AMPSORDER */
    int nRows, nColumns;    /* int for C array dimension sizes */
    int amp_xsize, amp_ysize;  /* int for amp array size (x/y in CCD coords) */
    int amp_xbeg, amp_xend;    /* int for beg and end of amp arrays on chip */
    int amp_ybeg, amp_yend;

    /* functions from calacs/lib */
    void parseWFCamps (char *acsamps, int chip, char *ccdamp);
    void TimeStamp (char *message, char *rootname);
    int PutKeyDbl (Hdr *hd, char *keyword, double value, char *comment);

    /* test whether this we've been given ACS WFC data, since the CTE algorithm
       is currently only valid for WFC data. */
    if (acs->detector != WFC_CCD_DETECTOR) {
        trlerror("(pctecorr) only valid for WFC CCD data, PCTECORR should be OMIT for all others.");
        return (status = ERROR_RETURN);
    }

    if (acs->printtime) {
        TimeStamp("Starting CTE correction...","");
    }

    /* need to figure out which amps are on this chip */
    ccdamp[0] = '\0'; /* "reset" the string for reuse */
    parseWFCamps(acs->ccdamp, acs->chip, ccdamp);

    /* loop over amps on this chip and do CTE correction */
    numamps = strlen(ccdamp);

    //Now loop over each amp and compute cte correction
    {unsigned nthAmp;
    for (nthAmp = 0; nthAmp < numamps; ++nthAmp)
    {
        sprintf(MsgText, "(pctecorr) Performing CTE correction for amp %c", ccdamp[nthAmp]);
        trlmessage(MsgText);

        /* get the amp letter and number where A:0, B:1, etc. */
        amploc = strchr(AMPSORDER, ccdamp[nthAmp]);
        ampID = *amploc - AMPSORDER[0];

        /* get amp array size */
        if (get_amp_array_size_acs_cte(acs, chipImage, ampID, amploc, ccdamp,
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
        ctePars->columnOffset = 0;//amp_xbeg;//acs->offsetx;
        ctePars->rowOffset = 0;//amp_ybeg;//acs->offsety;
        ctePars->razColumnOffset = nthAmp*nColumns;

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
        if ((status = extractAmp(&ampImage, chipImage, ampID)))
        {
            freeOnExit(&ptrReg);
            return (status);
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

        //CALCULATE THE SMOOTH READNOISE IMAGE
        trlmessage("(pctecorr) Calculating smooth readnoise image...");
        if (ctePars->noise_mit != 0)
        {
            trlerror("(pctecorr) Only noise model 0 implemented!");
            freeOnExit(&ptrReg);
            return (status = ERROR_RETURN);
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

        // do some smoothing on the data so we don't amplify the read noise.
        if ((status = cteSmoothImage(&columnMajorImage, &smoothedImage, ctePars, ctePars->rn_amp)))
        {
            freeOnExit(&ptrReg);
            return (status);
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

       trlmessage("(pctecorr) Running correction algorithm...");
       //perform CTE correction
       SingleGroup * cteCorrectedImage = &columnMajorImage;
       if ((status = inverseCTEBlur(&smoothedImage, cteCorrectedImage, &trapPixelMap, ctePars)))
       {
           freeOnExit(&ptrReg);
           return status;
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
            {unsigned k;
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            //This loop order is row major as more row major storage is accessed than column.
            //MORE: look into worth splitting out ops, prob needs a order swap (copy) so perhaps not worth it.
            for (k = 0; k < nRows; ++k)
            {
                {unsigned m;
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
                }}
            }}

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

        if ((status = insertAmp(chipImage, &ampImage, ampID)))
        {
            freeOnExit(&ptrReg);
            return (status);
        }

        // free mem used by our amp arrays
        freePtr(&ptrReg, &ampImage);
        freePtr(&ptrReg, &columnMajorImage);
        freePtr(&ptrReg, &smoothedImage);
    }}
    if (acs->printtime)
        TimeStamp("CTE corrections complete","");

    freeOnExit(&ptrReg);
    return (status);
}

static int extractAmp(SingleGroup * amp,  const SingleGroup * image, const unsigned ampID)
{
    extern int status;

    if (!amp || !amp->sci.data.data || !image || !image->sci.data.data)
        return (status = ALLOCATION_PROBLEM);

    //WARNING - assumes row major storage
    assert(amp->sci.data.storageOrder == ROWMAJOR && image->sci.data.storageOrder == ROWMAJOR);

    const unsigned nRows = amp->sci.data.ny;
    const unsigned nColumns = amp->sci.data.nx;

    unsigned rowSkipLength = image->sci.data.nx;
    unsigned offset = 0;
    if (ampID == AMP_B || ampID == AMP_D)
        offset = nColumns;
    else if (ampID != AMP_A && ampID != AMP_C)
    {
        trlerror("Amp number not recognized, must be 0-3.");
        return (status = ERROR_RETURN);
    }

    copyOffsetSingleGroup(amp, image, nRows, nColumns, 0, offset, nColumns, rowSkipLength);
    return status;
}

static int insertAmp(SingleGroup * image, const SingleGroup * amp, const unsigned ampID)
{
    extern int status;

    if (!amp || !amp->sci.data.data || !image || !image->sci.data.data)
        return (status = ALLOCATION_PROBLEM);

    //WARNING - assumes row major storage
    assert(amp->sci.data.storageOrder == ROWMAJOR && image->sci.data.storageOrder == ROWMAJOR);

    const unsigned nRows = amp->sci.data.ny;
    const unsigned nColumns = amp->sci.data.nx;

    unsigned rowSkipLength = image->sci.data.nx;
    unsigned offset = 0;
    if (ampID == AMP_B || ampID == AMP_D)
        offset = nColumns;
    else if (ampID != AMP_A && ampID != AMP_C)
    {
        trlerror("Amp number not recognized, must be 0-3.");
        return (status = ERROR_RETURN);
    }

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

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    //WARNING - assumes row major storage
    assert(amp->storageOrder == ROWMAJOR);

    const unsigned nColumns = amp->nx;
    const unsigned nRows = amp->ny;

    //Flip about y axis, i.e. about central column
    if (ampID == AMP_B || ampID == AMP_D)
    {
        //grab a row, flip it, put it back
        const unsigned rowLength = nColumns;
        {unsigned i;
#ifdef _OPENMP
        #pragma omp parallel for shared(amp) private(i) schedule(static)
#endif
        for (i = 0; i < nRows; ++i)
        {
            //find row
            float * row = amp->data + i*nColumns;
            //flip right to left
            float tempPixel;
            {unsigned j;
            for (j = 0; j < rowLength/2; ++j)
            {
                tempPixel = row[j];
                row[j] = row[rowLength-1-j];
                row[rowLength-1-j] = tempPixel;
            }}
        }}
    }

    //Only thing left is to flip AB chip upside down
    //Flip about x axis, i.e. central row
    if (ampID == AMP_A || ampID == AMP_B)
    {
        //either physically align all or propagate throughout a mechanism to work on the array upside down (flip b quad though)
        //we'll just flip all for now. See if there's info in the header specifying amp location rel to pix in file,
        //i.e. a way to know whether the chip is 'upside down'. Could then reverse cte corr trail walk direction
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
                {unsigned i;
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
                }}
            }
            if (tempRow)
                free(tempRow);
        }//close parallel block
        if (allocationFail)
        {
            sprintf(MsgText, "Out of memory for 'tempRow' in 'alignAmpData'");
            trlerror(MsgText);
            return (status = OUT_OF_MEMORY);
        }
    }

    return status;
}
