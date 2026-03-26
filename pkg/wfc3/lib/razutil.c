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
            trlerror("Invalid group number passed to makedqRAZ");
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
            trlerror("Invalid group number passed to makedqRAZ");
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
            trlerror("Invalid group number passed to makeSciSingleRAZ");
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
            trlerror("Invalid group number passed to makedqRAZ");
            return(status=INVALID_VALUE);
        }
    }

    return(status);
}


/* This is Jay's quick-and-easy FITS writer that can be used to
   write out temporary RAZ images for internal testing.
   Do not use for production in pipeline.
*/
int writfits_r4(char *filename, float *data, int nx, int ny) {

    typedef unsigned char Byte;

    char  cbuff[2880];
    Byte  bbuff[2880];
    Byte  bbufr[2880];
    float rbuff[ 720];

    int i;
    int n;
    int NRs, NRo;
    int nb, NBs;

    Byte btemp4[4];
    Byte  *pByte;

    for (i=0;i<2880;i++) {
        cbuff[i] = ' ';
    }

    /* simplest possible header */
    sprintf(&cbuff[0+ 0*80],"SIMPLE  =                    T");
    sprintf(&cbuff[0+ 1*80],"BITPIX  =                  -32");
    sprintf(&cbuff[0+ 2*80],"NAXIS   =                    2");
    sprintf(&cbuff[0+ 3*80],"NAXIS1  = %20d",nx);
    sprintf(&cbuff[0+ 4*80],"NAXIS2  = %20d",ny);
    sprintf(&cbuff[0+ 5*80],"DATATYPE= 'REAL*4'            ");
    sprintf(&cbuff[0+ 6*80],"COMMENT                       ");
    sprintf(&cbuff[0+ 7*80],"COMMENT    Jay's Quickie      ");
    sprintf(&cbuff[0+ 8*80],"COMMENT    C fits writer      ");
    sprintf(&cbuff[0+ 9*80],"COMMENT    writfits_r4()      ");
    sprintf(&cbuff[0+10*80],"COMMENT                       ");
    sprintf(&cbuff[0+11*80],"END                           ");

    for (i=0;i<2880;i++) {
        if (cbuff[i]==0) {
            cbuff[i] = ' ';
        }
    }

    FILE *fp = fopen(filename,"w");
    fwrite(cbuff,2880,1,fp);

    NRs = nx*ny;
    NBs = 1 + ((nx*ny-1)/720);  /* the -1 covers the case where the data have an exact multiple
                                   of 2880 bytes; it can happen! */

    /* distill the data into 2880-byte buffers
       we must flip the bytes, since INTEL machines
       have the opposite Endian-ness to the FITS
       standard
    */
    for (nb=0;nb<NBs;nb++) {
         NRo = nb*720;
         for (n=0;n<720;n++) {
              rbuff[n] = 0;
              if (NRo+n<NRs) {
                  rbuff[n] = data[NRo+n];
              } /* make sure we haven't gone off the end */
         }
         pByte = (Byte *)rbuff;
         for (n=0;n<2880;n++) {
             bbuff[n] = *(pByte+n);
         }
         for (n=0;n<720;n++) {
              btemp4[0] = bbuff[0+n*4];
              btemp4[1] = bbuff[1+n*4];
              btemp4[2] = bbuff[2+n*4];
              btemp4[3] = bbuff[3+n*4];
              bbufr[0+n*4] = btemp4[3];
              bbufr[1+n*4] = btemp4[2];
              bbufr[2+n*4] = btemp4[1];
              bbufr[3+n*4] = btemp4[0];
         }
         fwrite(bbufr,2880,1,fp);
    }
    fclose(fp);

    return(0);
}
