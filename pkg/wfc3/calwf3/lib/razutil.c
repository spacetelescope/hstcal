# include <string.h>
# include "hstio.h"
# include "wf3.h"
# include "wf3err.h"

/* These routines facilitate moving between the regular WFC3 image structre
   to the RAZ image structure used in the CTE correction and Sink pixel flagging
   

   convert a raw file to raz file: CDAB longwise amps, save data array
   for comparison with what jay has during testing

   In the raz image, each quadrant has been rotated such that the readout amp is located at the lower left. 
   The reoriented four quadrants are then arranged into a single 8412x2070 image (science pixels plus overscan), 
   with amps C, D, A, and B, in that order. In the raz image, pixels are all parallel-shifted down, 
   then serial-shifted to the left.

   Megan Sosey, May 2015

*/


/*convert  floating point arrays into raz format*/
int makeFloatRaz(FloatTwoDArray *cd, FloatTwoDArray *ab, FloatTwoDArray  *raz){

    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

    for (i=0; i<subcol; i++){
        for (j=0;j<RAZ_ROWS; j++){
            memcpy( &PPix(raz,i,j), &PPix(cd,i,j), sizeof(float));
            memcpy( &PPix(raz,i+subcol,j), &PPix(cd,subcol*2-i-1,j),sizeof(float));
            memcpy( &PPix(raz,i+2*subcol,j), &PPix(ab,i,RAZ_ROWS-j-1),sizeof(float));
            memcpy( &PPix(raz,i+3*subcol,j), &PPix(ab,subcol*2-i-1,RAZ_ROWS-j-1), sizeof(float));
        }
    }

        
    return(status);      
}

/*convert the sci extension*/
int makesciRAZ(SingleGroup *cd, SingleGroup *ab, SingleGroup *raz){

    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

        for (i=0; i<subcol; i++){
            for (j=0; j<RAZ_ROWS; j++){
                memcpy( &Pix(raz->sci.data,i,j), &Pix(cd->sci.data,i,j), sizeof(float));
                memcpy( &Pix(raz->sci.data,i+subcol,j), &Pix(cd->sci.data,subcol*2-i-1,j),sizeof(float));
                memcpy( &Pix(raz->sci.data,i+2*subcol,j), &Pix(ab->sci.data,i,RAZ_ROWS-j-1),sizeof(float));
                memcpy( &Pix(raz->sci.data,i+3*subcol,j), &Pix(ab->sci.data,subcol*2-i-1,RAZ_ROWS-j-1), sizeof(float));
            }
        }
        
    return(status);      

}


/* Transform a RAZ format image back into the separate input arrays calwf3 likes*/
int undosciRAZ(SingleGroup *cd, SingleGroup *ab, SingleGroup *raz){
    
    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;
        
    /*REVERSE THE AMPS TO THE RAW FORMAT*/
    for (i=0; i< subcol; i++){
        for (j=0; j<RAZ_ROWS; j++){
             memcpy( &Pix(cd->sci.data,i,j), &Pix(raz->sci.data,i,j), sizeof(float));
             memcpy( &Pix(cd->sci.data,subcol*2-i-1,j), &Pix(raz->sci.data,i+subcol,j),sizeof(float));
             memcpy( &Pix(ab->sci.data,i,RAZ_ROWS-j-1), &Pix(raz->sci.data,i+2*subcol,j),sizeof(float));
             memcpy( &Pix(ab->sci.data,subcol*2-i-1,RAZ_ROWS-j-1), &Pix(raz->sci.data,i+3*subcol,j), sizeof(float));
        }
    }
    return (status);

}

/*convert the DQ extension*/
int makedqRAZ(SingleGroup *cd, SingleGroup *ab, SingleGroup *raz){
    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;

        for (i=0; i<subcol; i++){
            for (j=0;j<RAZ_ROWS; j++){
                memcpy( &Pix(raz->dq.data,i,j), &Pix(cd->dq.data,i,j), sizeof(short));
                memcpy( &Pix(raz->dq.data,i+subcol,j), &Pix(cd->dq.data,subcol*2-i-1,j),sizeof(short));
                memcpy( &Pix(raz->dq.data,i+2*subcol,j), &Pix(ab->dq.data,i,RAZ_ROWS-j-1),sizeof(short));
                memcpy( &Pix(raz->dq.data,i+3*subcol,j), &Pix(ab->dq.data,subcol*2-i-1,RAZ_ROWS-j-1), sizeof(short));
            }
        }
        
    return(status);      

}


/* Transform dqin a  RAZ format image back into the separate input arrays calwf3 likes*/
int undodqRAZ(SingleGroup *cd, SingleGroup *ab, SingleGroup *raz){
    
    extern int status;
    int subcol = (RAZ_COLS/4); /* for looping over quads  */
    int i,j;
        
    /*REVERSE THE AMPS TO THE RAW FORMAT*/
    for (i=0; i< subcol; i++){
        for (j=0; j<RAZ_ROWS; j++){
             memcpy( &Pix(cd->dq.data,i,j), &Pix(raz->dq.data,i,j), sizeof(short));
             memcpy( &Pix(cd->dq.data,subcol*2-i-1,j), &Pix(raz->dq.data,i+subcol,j),sizeof(short));
             memcpy( &Pix(ab->dq.data,i,RAZ_ROWS-j-1), &Pix(raz->dq.data,i+2*subcol,j),sizeof(short));
             memcpy( &Pix(ab->dq.data,subcol*2-i-1,RAZ_ROWS-j-1), &Pix(raz->dq.data,i+3*subcol,j), sizeof(short));
        }
    }
    return (status);

}
