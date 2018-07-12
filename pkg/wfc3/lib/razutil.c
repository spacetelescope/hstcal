# include <string.h>
# include "hstio.h"
# include "wf3.h"
# include "hstcalerr.h"

/* These routines facilitate moving between the regular WFC3 image structre
   to the RAZ image structure used in the CTE correction and Sink pixel flagging


   convert a raw file to raz file: CDAB longwise amps, save data array
   for comparison with what jay has during testing

   In the raz image, each quadrant has been rotated such that the readout amp is located at the lower left.
   The reoriented four quadrants are then arranged into a single 8412x2070 image (science pixels plus overscan),
   with amps C, D, A, and B, in that order. In the raz image, pixels are all parallel-shifted down,
   then serial-shifted to the left.

   The code for the DQ arrays and plain float arrays only converts one chip at a time so that it can be run through the regular
   wf3ccd pipeline which operates on 1 chip at a time.

   Megan Sosey, May 2015

*/



/*convert the sci and dq extensions to the long format*/
int makeRAZ(SingleGroup *cd, SingleGroup *ab, SingleGroup *raz){

    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

        for (i=0; i<subcol; i++){
            for (j=0; j<RAZ_ROWS; j++){
                Pix(raz->sci.data,i,j)=Pix(cd->sci.data,i,j);
                Pix(raz->sci.data,i+subcol,j)=Pix(cd->sci.data,subcol*2-i-1,j);
                Pix(raz->sci.data,i+2*subcol,j)=Pix(ab->sci.data,i,RAZ_ROWS-j-1);
                Pix(raz->sci.data,i+3*subcol,j)=Pix(ab->sci.data,subcol*2-i-1,RAZ_ROWS-j-1);

                Pix(raz->dq.data,i,j)=Pix(cd->dq.data,i,j);
                Pix(raz->dq.data,i+subcol,j)=Pix(cd->dq.data,subcol*2-i-1,j);
                Pix(raz->dq.data,i+2*subcol,j)=Pix(ab->dq.data,i,RAZ_ROWS-j-1);
                Pix(raz->dq.data,i+3*subcol,j)=Pix(ab->dq.data,subcol*2-i-1,RAZ_ROWS-j-1);
            }
        }

    return(status);

}


/* Transform a RAZ format image back into the separate input arrays calwf3 likes*/
int undoRAZ(SingleGroup *cd, SingleGroup *ab, SingleGroup *raz){

    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    /*REVERSE THE AMPS TO THE RAW FORMAT*/
    for (i=0; i< subcol; i++){
        for (j=0; j<RAZ_ROWS; j++){
             Pix(cd->sci.data,i,j) = Pix(raz->sci.data,i,j);
             Pix(cd->sci.data,subcol*2-i-1,j) = Pix(raz->sci.data,i+subcol,j);
             Pix(ab->sci.data,i,RAZ_ROWS-j-1) = Pix(raz->sci.data,i+2*subcol,j);
             Pix(ab->sci.data,subcol*2-i-1,RAZ_ROWS-j-1) = Pix(raz->sci.data,i+3*subcol,j);

             Pix(cd->dq.data,i,j) = Pix(raz->dq.data,i,j);
             Pix(cd->dq.data,subcol*2-i-1,j) = Pix(raz->dq.data,i+subcol,j);
             Pix(ab->dq.data,i,RAZ_ROWS-j-1) = Pix(raz->dq.data,i+2*subcol,j);
             Pix(ab->dq.data,subcol*2-i-1,RAZ_ROWS-j-1) = Pix(raz->dq.data,i+3*subcol,j);


        }
    }
    return (status);

}


/***********THE ROUTINES BELOW OPERATE ON A SINGLE GROUP***********/

/*
The DQ versions here are written to work with SINKDETECT
and work on 1 chip at a time, with Half the Columns of the
full size raz image and all the rows
*/

int makedqRAZ(SingleGroup *x, SingleGroup *raz){
    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    if (x->group_num == 1){
        for (i=0; i<subcol; i++){
            for (j=0;j<RAZ_ROWS; j++){
                Pix(raz->dq.data,i,j) = Pix(x->dq.data,i,j);
                Pix(raz->dq.data,i+subcol,j) = Pix(x->dq.data,subcol*2-i-1,j);
            }
        }

    } else {
        if (x->group_num == 2){
            for (i=0; i<subcol; i++){
                for (j=0;j<RAZ_ROWS; j++){
                    Pix(raz->dq.data,i,j) = Pix(x->dq.data,i,RAZ_ROWS-j-1);
                    Pix(raz->dq.data,i+subcol,j) = Pix(x->dq.data,subcol*2-i-1,RAZ_ROWS-j-1);
                }
            }

        } else {
            trlmessage("Invalid group number passed to makedqRAZ");
            return(status=INVALID_VALUE);
        }
    }

   return(status);

}


/* Transform dq in a  RAZ format image back into the separate input arrays calwf3 likes*/
int undodqRAZ(SingleGroup *x, SingleGroup *raz){

    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    if (x->group_num == 1){
        for (i=0; i< subcol; i++){
            for (j=0; j<RAZ_ROWS; j++){
                 Pix(x->dq.data,i,j) = Pix(raz->dq.data,i,j);
                 Pix(x->dq.data,2*subcol-i-1,j) = Pix(raz->dq.data,i+subcol,j);
            }
        }
    } else {
        if (x->group_num == 2){
            for (i=0; i< subcol; i++){
                for (j=0; j<RAZ_ROWS; j++){
                     Pix(x->dq.data,i,RAZ_ROWS-j-1) = Pix(raz->dq.data,i,j);
                     Pix(x->dq.data,subcol*2-i-1,RAZ_ROWS-j-1) = Pix(raz->dq.data,i+subcol,j);
                }
            }
        } else {
            trlmessage("Invalid group number passed to makedqRAZ");
            return(status=INVALID_VALUE);
        }
    }

    return (status);

}

/*convert the science image of a single group to RAZ
 for use in the SINK pixel detection*/

int makeSciSingleRAZ(SingleGroup *x, SingleGroup *raz){
    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    if (x->group_num == 1){
        for (i=0; i<subcol; i++){
            for (j=0;j<RAZ_ROWS; j++){
                Pix(raz->sci.data,i,j) = Pix(x->sci.data,i,j);
                Pix(raz->sci.data,i+subcol,j) = Pix(x->sci.data,subcol*2-i-1,j);
            }
        }

    } else {
        if (x->group_num == 2){
            for (i=0; i<subcol; i++){
                for (j=0;j<RAZ_ROWS; j++){
                    Pix(raz->sci.data,i,j) = Pix(x->sci.data,i,RAZ_ROWS-j-1);
                    Pix(raz->sci.data,i+subcol,j) = Pix(x->sci.data,subcol*2-i-1,RAZ_ROWS-j-1);
                }
            }

        } else {
            trlmessage("Invalid group number passed to makeSciSingleRAZ");
            return(status=INVALID_VALUE);
        }
    }

   return(status);

}

/*convert  floating point arrays into raz format*/
int makeFloatRaz(FloatTwoDArray *x, FloatTwoDArray  *raz, int group){

    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    if (group == 1){
        for (i=0; i<subcol; i++){
            for (j=0;j<RAZ_ROWS; j++){
                PPix(raz,i,j) = PPix(x,i,j);
                PPix(raz,i+subcol,j) = PPix(x,subcol*2-i-1,j);
            }
        }
    } else {
        if (group == 2){
            for (i=0; i<subcol; i++){
                for (j=0;j<RAZ_ROWS; j++){
                    PPix(raz,i,j) = PPix(x,i,RAZ_ROWS-j-1);
                    PPix(raz,i+subcol,j) =PPix(x,subcol*2-i-1,RAZ_ROWS-j-1);
                }
            }
        } else {
            trlmessage("Invalid group number passed to makedqRAZ");
            return(status=INVALID_VALUE);
        }
    }

    return(status);
}
