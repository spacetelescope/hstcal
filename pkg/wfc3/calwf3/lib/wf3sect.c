/* This file contains:
	initWF3sect
	allocWF3sect
	freeWF3sect
	copySectLine
	getWF3sect
*/

/*
# include <xclio.h>
*/

# include <stdio.h>
# include <string.h>
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"


void initWF3sect (WF3sect *x) {
        x->filename     = NULL;
        x->nlines       = 0;
	x->npix         = 0;
	x->line_num	=0;
	x->group_num    =0;
        x->phdr_loaded = False;
        x->globalhdr   = NULL;
        x->sci         = NULL;
	x->err         = NULL;
        x->dq          = NULL;
}
 
int allocWF3sect (WF3sect *x, int npix, int nlines) {

	int i;
		
        if (x->globalhdr == NULL) {     
            x->globalhdr = (Hdr *) calloc (1,sizeof(Hdr));
            if (x->globalhdr == NULL) return (-1);
        }
		
	x->sci = (SciHdrLine *) calloc(nlines, sizeof(SciHdrLine));
	x->err = (ErrHdrLine *) calloc(nlines, sizeof(ErrHdrLine));
	x->dq = (DQHdrLine *) calloc(nlines, sizeof(DQHdrLine));
		
	for (i=0; i < nlines; i++) {
        	if (allocFloatHdrLine (&(x->sci[i]),npix)) return (-1);
        	if (allocFloatHdrLine (&(x->err[i]),npix)) return (-1);
        	if (allocShortHdrLine (&(x->dq[i]),npix))  return (-1);
	} 
	x->nlines = nlines;
	x->npix = npix;

        return (0);
}

void freeWF3sect (WF3sect *x) {

	int i;

	for (i=0;i < x->nlines; i++) {
        	freeFloatHdrLine (&(x->sci[i]));
        	freeFloatHdrLine (&(x->err[i]));
        	freeShortHdrLine (&(x->dq[i]));
	} 
		
	free (x->sci);
	free (x->err);
	free (x->dq);
		
        freeHdr (x->globalhdr);
        if (x->globalhdr != NULL)
                free (x->globalhdr);
        if (x->filename != NULL)
                free (x->filename);
        initWF3sect (x);
}

void closeWF3sect (WF3sect *x) {

	closeImage (x->sci->iodesc);
	closeImage (x->err->iodesc);
	closeImage (x->dq->iodesc);

}

/* 
	This function extracts a SingleGroupLine from an WF3sect.
	Calling function will have to initialize and allocate/free space for 
	output SingleGroupLine.
*/
void copySectLine (WF3sect *x, int line, SingleGroupLine *y) {

/* Arguments:
WF3sect *x			i: Section of input data
int line			i: line from input data to extract
SingleGroupLine *y	o: Output line of data
*/
	y->line_num = x->line_num + line;
	y->group_num = x->group_num;
	y->phdr_loaded = x->phdr_loaded;
        
        /* Copy the acquired image line into the local buffer */
        memcpy (y->sci.line, x->sci[line].line, x->npix * sizeof(float));
        memcpy (y->err.line, x->err[line].line, x->npix * sizeof(float));
        memcpy (y->dq.line, x->dq[line].line, x->npix * sizeof(short));

}

/* This functions reads in an image section line-by-line
    The SingleGroupLine y contains the pointer to the correct extension
    already, and gets managed outside this function. This avoids opening
    this SingleGroup every time we get a section.

*/
void getWF3sect (char *fname, SingleGroupLine *y, int y0, int nlines,
		 WF3sect *x){

/* Arguments:
char *fname			i: input filename
SingleGroupLine *y		i: scratch space for line of data, 
int y0				i: y offset in input image
int nlines			i: number of lines to read
WF3sect *x			o: output section
*/
	int i;
	
	for (i = 0; i < nlines; i++) {
	     getSingleGroupLine (fname, i+y0, y);
		
	     x->line_num = i;
	     x->group_num = y->group_num;
	     x->phdr_loaded = y->phdr_loaded;
		        
	     /* Copy the acquired image line into the local buffer */
	     memcpy (x->sci[i].line, y->sci.line, y->sci.tot_nx*sizeof(float));
	     memcpy (x->err[i].line, y->err.line, y->sci.tot_nx*sizeof(float));
	     memcpy (x->dq[i].line,  y->dq.line,  y->sci.tot_nx*sizeof(short)); 
	}

}

