# include <stdio.h>
# include <string.h>

#include "hstcal.h"
# include "ximio.h"	/* defines IRAF imio functions */
# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "trlbuf.h"

extern int status;

/* IMAGEIO: Contains routines for reading and writing image
** data to be calibrated. These are the routines:
**
** getPriHdr : Reads the primary header of an input file.
**
** getRawData: Reads one group of raw data from input file.
**
** putCalData: Writes one group of calibrated data to output file.
**
** putCalDataSect: Writes a section of one group of calibrated data to an
**                 output file.
**
** putMultiCalData: Writes multiple groups of calibrated data to output file.
**
** copyGroup: Copy the contents of one group to a new one.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	20-Mar-2002	Removed ckImgSize routine (not needed for WFC3).
** H.Bushouse	08-May-2002	Modified to use trlopenerr and trlreaderr.
** H.Bushouse	16-Feb-2007	Enhanced copyGroup to only copy filename if
**				input name is not Null. Added new 
**				putCalDataSect routine.
*/

/* GETPRIHDR: Read primary header from an input file. */

int getPriHdr (WF3Info *wf3, Hdr *PriHdr) {

/* Arguments:
**	wf3	i: WFC3 info structure
**	PriHdr	o: header structure
*/

	/* Local variables */
	IODescPtr im;		/* file descriptor pointer */

	/* Open the image */
	im = openInputImage (wf3->input, "", 0);
	if (hstio_err()) {
	    trlopenerr (wf3->input);
	    return (status = 1);
	}

	/* Initialize the header data structure */
	initHdr (PriHdr);

	/* Read the header */
	if (getHeader (im, PriHdr))
	    status = 1;
	if (hstio_err() || status) {
	    trlreaderr (wf3->input);
	    closeImage (im);
	    freeHdr (PriHdr);
	    return (status = 1);
	}

	/* Close the image */
	closeImage (im);

	/* Successful return */
	return (status = 0);
}

/* GETRAWDATA: Read raw data from input file. One group is read. */

int getRawData (WF3Info *wf3, MultiNicmosGroup *in) {

/* Arguments:
**	wf3	i: WFC3 info structure
**	in	o: input image data
*/

	/* Function definitions */
	int getGroupInfo (WF3Info *, SingleNicmosGroup *);

	/* Initialize the input data structure */
	initMultiNicmosGroup (in);

	/* Allocate the group data structures */
	if (allocMultiNicmosGroup (in, wf3->ngroups))
	    return (status = 1);

	/* Read the global header */
	if (getMultiNicmosGroupHdr (wf3->input, in))
	    status = 1;
	if (hstio_err() || status) {
	    trlreaderr (wf3->input);
	    freeMultiNicmosGroup (in);
	    return (status = 1);
	}

	/* Read the input data groups */
	for (wf3->group = 1; wf3->group <= wf3->ngroups; wf3->group++) {
	     if (getMultiNicmosGroup (in, wf3->group-1, wf3->group))
		 status = 1;
	     if (hstio_err() || status) {
		 trlreaderr (wf3->input);
		 freeMultiNicmosGroup (in);
		 return (status = 1);
	     }

	     /* Get group-specific information from the file headers
	     ** and check for null input data */
	     if (getGroupInfo (wf3, &(in->group[wf3->group-1]))) {
		 freeMultiNicmosGroup (in);
		 return (status);
	     }
	}

	/* Successful return */
	return (status = 0);
}

/* PUTCALDATA: Write calibrated data to a single-group ouput file. */

int putCalData (SingleNicmosGroup *out, char *fname) {

/* Arguments:
**	out	i: image data to be written to file
**	fname	i: output file name
*/

	/* Function definitions */
	int updateHdr (SingleNicmosGroup *, char *);

	/* Update output image header keywords */
	if (updateHdr (out, fname))
	    return (status);

	/* Write a single group header and data */
	if (putSingleNicmosGroup (fname, out->group_num, out, 0))
	    status = 1;
	if (hstio_err() || status) {
	    sprintf (MsgText, "Can't write to output image %s", fname);
	    trlerror (MsgText);
	    return (status = 1);
	}

	/* Reset the DATAMIN/MAX keywords in the output file */

	/* Successful return */
	return (status = 0);
}

/* PUTCALDATASECT: Write calibrated data section to a single-group ouput file.*/

int putCalDataSect (SingleNicmosGroup *out, char *fname, int x1, int y1,
		    int xsize, int ysize) {

/* Arguments:
**	out	i: image data to be written to file
**	fname	i: output file name
**	x1	i: x corner of section
**	y1	i: y corner of section
**	xsize	i: x size of section
**	ysize	i: y size of section
*/

	/* Function definitions */
	int updateHdr (SingleNicmosGroup *, char *);

	/* Update output image header keywords */
	if (updateHdr (out, fname))
	    return (status);

	/* Write a single group header and data section */
	if (putSingleNicmosGroupSect (fname, out->group_num, out, x1, y1,
				      xsize, ysize, 0))
	    status = 1;
	if (hstio_err() || status) {
	    sprintf (MsgText, "Can't write to output image %s", fname);
	    trlerror (MsgText);
	    return (status = 1);
	}

	/* Reset the DATAMIN/MAX keywords in the output file */

	/* Successful return */
	return (status = 0);
}

/* PUTMULTICALDATA: Write multiple groups of calibrated data to ouput file.
** The file is created and the primary header is written when writing the 
** first data group. */

int putMultiCalData (MultiNicmosGroup *out, char *fname) {

/* Arguments:
**	out	i: image data to be written to file
**	fname	i: output file name
*/

	/* Local variables */
	int i;				/* loop index */

	/* Function definitions */
	int updateHdr (SingleNicmosGroup *, char *);

	/* Update output image header keywords */
	for (i = 0; i < out->ngroups; i++) {
	     if (updateHdr (&(out->group[i]), fname))
		 return (status);
	}

	/* Write the group headers and data */
	for (i = 0; i < out->ngroups; i++) {

	     if (putSingleNicmosGroup (fname, i+1, &(out->group[i]), 0))
		 status = 1;
	     if (hstio_err() || status) {
		 sprintf (MsgText, "Can't write to output image %s", fname);
		 trlerror (MsgText);
		 return (status = 1);
	     }

	     /* Reset the DATAMIN/MAX keywords in the output file */
	}

	/* Successful return */
	return (status = 0);
}

/* COPYGROUP: Copy one group structure to another.
** Initializes and allocates the structure for the new group.
*/
 
int copyGroup (SingleNicmosGroup *to, SingleNicmosGroup *from) {
 
/* Arguments:
**      to      o: new group
**      from    i: group to be copied
*/
 
        /* Initialize and allocate new group structure */
        initSingleNicmosGroup (to);
        if (allocSingleNicmosGroup 
		(to, from->sci.data.nx, from->sci.data.ny) == -1) {
            sprintf (MsgText, "in copyGroup; can't allocate new group");
            trlerror (MsgText);
            return (status = 1);
        }
 
	/* Copy the filename */
	if (from->filename != NULL) {
	    to->filename = (char *)calloc((strlen(from->filename)+1),sizeof(char));
	    strcpy (to->filename, from->filename);
	}

	/* Copy the group number */
	to->group_num = from->group_num;

        /* Copy the global and extension headers */
        if (copyHdr (to->globalhdr,   from->globalhdr))   return(status=1);
        if (copyHdr (&(to->sci.hdr),  &(from->sci.hdr)))  return(status=1);
        if (copyHdr (&(to->err.hdr),  &(from->err.hdr)))  return(status=1);
        if (copyHdr (&(to->dq.hdr),   &(from->dq.hdr)))   return(status=1);
        if (copyHdr (&(to->smpl.hdr), &(from->smpl.hdr))) return(status=1);
        if (copyHdr (&(to->intg.hdr), &(from->intg.hdr))) return(status=1);
 
        /* Copy the extension images data arrays */
        memcpy (to->sci.data.data, from->sci.data.data,
                from->sci.data.nx*from->sci.data.ny*sizeof(float));
        memcpy (to->err.data.data, from->err.data.data,
                from->sci.data.nx*from->sci.data.ny*sizeof(float));
        memcpy (to->dq.data.data,  from->dq.data.data,
                from->sci.data.nx*from->sci.data.ny*sizeof(short));
        memcpy (to->smpl.data.data,from->smpl.data.data,
                from->sci.data.nx*from->sci.data.ny*sizeof(short));
        memcpy (to->intg.data.data,from->intg.data.data,
                from->sci.data.nx*from->sci.data.ny*sizeof(float));
 
        /* Successful return */
        return (status = 0);
}

