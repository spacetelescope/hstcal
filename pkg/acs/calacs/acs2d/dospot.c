# include "ximio.h"
# include "hstio.h"
# include <math.h>

# include "err.h"
# define    PIX(v,i,j,nx)   v[(i) + (j) * (nx)]

int readSpotImage(char *spotname, SingleGroup *inspot, SingleGroupLine *spotline) {
    
    extern int status;
    int dimx;
	
    initSingleGroup (inspot);
	/* Open the input image. */
	getSingleGroup (spotname, 1, inspot);
	if (hstio_err())
	    return (status = OPEN_FAILED);

    dimx = inspot->sci.data.nx;
        
    initSingleGroupLine(spotline);    
    allocSingleGroupLine(spotline, dimx);  
  
    return (status);
}


void copySpotLine(SingleGroup *spot, int line, SingleGroupLine *spotline) {

/* Calling function will need to allocate/deallocate space for 
    SingleGroupLine object.
*/    
    int dimx;
    int j;
    
    dimx = spot->sci.data.nx;
        
    /* Copy data from SingleGroup into SingleGroupLine */ 
    for (j=0;j<dimx;j++) {   
     spotline->sci.line[j] = Pix(spot->sci.data,j,line);
     spotline->err.line[j] = Pix(spot->err.data,j,line);
     spotline->dq.line[j] = Pix(spot->dq.data,j,line);
    }
}

int shiftSpot (SingleGroup *inspot, float dx, float dy, SingleGroup *outspot) {

/* arguments:
char *spotname          i: name of input spot flat
float dx,dy             i: shift in x,y
int nx,ny               i: size of output shifted spot flat
float **outspot          o: output product
*/
    extern int status;

	int i, j;
    int dimx, dimy;
    float xout,yout, xo, yo;
    float value;

    float binterp(float, float, int, int, SingleGroup *);

	/* science data dimensions */
    dimx = inspot->sci.data.nx;
    dimy = inspot->sci.data.ny;
    
    /*For each output line... */    
	for (j = 0;  j < dimy;  j++) {
        /* ... apply linear shift to pixel positions... */
        /* This does not allow for scale or rotation changes */
        xout = -dx;
        yout = j - dy;
        
	    for (i = 0;  i < dimx;  i++) {
            xo = xout + i;
            yo = yout;
            if ( 0 <= xo && xo < dimx && 0 <= yo && yo < dimy){
                value = binterp(xo,yo, dimx,dimy, inspot);
                
                Pix(outspot->sci.data,i,j) = value;
                Pix(outspot->err.data,i,j) = sqrt(value);
            } else {
                Pix(outspot->sci.data,i,j) = 1.0;
                Pix(outspot->err.data,i,j) = 1.0;
            }
        }
    }
	return (status);
}

float binterp (float x, float y, int nxpix, int nypix, SingleGroup *img) {

    float value;
    float sx, sy, tx, ty;
    float hold12, hold21, hold22;
    int nx, ny;
    
    nx = (int)x;
    ny = (int)y;
        
    sx = x - nx;
    tx = 1. - sx;
    sy = y - ny;
    ty = 1. - sy;
    
    if (nx >= nxpix) {      
        hold21 = 2. * Pix(img->sci.data,nx,ny) - Pix(img->sci.data,nx-1,ny);
    } else {     
        hold21 = Pix(img->sci.data,nx+1,ny);
    } 

    if (ny >= nypix) {
        hold12 = 2. * Pix(img->sci.data,nx,ny) - Pix(img->sci.data,nx,ny-1);
    } else {
        hold12 = Pix(img->sci.data,nx,ny+1);
    }    

    if (nx >= nxpix && ny >= nypix) {
        hold22 = 2. * hold21 - (2. * Pix(img->sci.data,nx,ny-1) - Pix(img->sci.data,nx-1,ny-1));
    } else if (nx >= nxpix) {
        hold22 = 2. * hold12 - Pix(img->sci.data,nx-1,ny+1);
    } else if (ny >= nypix) { 
       hold22 = 2. * hold21 - Pix(img->sci.data,nx+1,ny-1);
    } else {     
       hold22 = Pix(img->sci.data,nx+1,ny+1);
    } 

    value = tx * ty * Pix(img->sci.data,nx,ny) + sx * ty * hold21 + 
            sy * tx * hold12 + sx * sy * hold22 ;

    return (value);
}

float ninterp (float x, float y, int nxpix, int nypix, SingleGroup *img) {
    float value;
    
    value = Pix(img->sci.data, (int)(y+0.5),(int)(x+0.5));
    
    return(value);
}
