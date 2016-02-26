# include <stdio.h>

# include "hstio.h"     /* defines HST I/O functions */
# include "wf3.h"

extern int status;

/* UPDATEHDR: Update various keywords in output image header at the
** end of calibration processing. The keywords that get updated are:
**
** NEXTEND: set to 5*wf3.group
** DATAMIN, DATAMAX: updated in the header of all extensions
** FILENAME: set to name of file being written
** CAL_VER: set to WF3_CAL_VER
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
**
*/

int updateHdr (SingleNicmosGroup *input, char *fname) {

/* Arguments:
**	input	io: input image
**	fname	 i: output file name
*/

	/* Local variables */

	/* Function definitions */
	void UCalVer (Hdr *);
	int updMinMaxf (FloatHdrData *);
	int updMinMaxs (ShortHdrData *);

	/* Update the "NEXTEND" keyword */
	if (putKeyI (input->globalhdr, "NEXTEND", 5*input->group_num, ""))
	    return (status = 1);

	/* Update the "FILENAME" keyword */
	if (putKeyS (input->globalhdr, "FILENAME", fname, ""))
	    return (status = 1);

	/* Update the "CAL_VER" keyword */
	UCalVer (input->globalhdr);

	/* Update the datamin/datamax keywords of the SCI image */
	if (updMinMaxf (&input->sci))
	    return (status);

	/* Update the datamin/datamax keywords of the ERR image */
	if (updMinMaxf (&input->err))
	    return (status);

	/* Update the datamin/datamax keywords of the DQ  image */
	if (updMinMaxs (&input->dq))
	    return (status);

	/* Update the datamin/datamax keywords of the SAMP image */
	if (updMinMaxs (&input->smpl))
	    return (status);

	/* Update the datamin/datamax keywords of the TIME image */
	if (updMinMaxf (&input->intg))
	    return (status);

	/* Successful return */
	return (status = 0);
}

int updMinMaxf (FloatHdrData *image) {

	/* Local variables */
	float min, max;		/* data min and max */

	/* Function definitions */
	void compMinMaxf (FloatTwoDArray *, float *, float *);

	/* Compute the min and max of the data array */
	compMinMaxf (&image->data, &min, &max);

	/* Update the DATAMIN keyword */
	if (putKeyF (&image->hdr, "DATAMIN", min, ""))
	    return (status = 1);

	/* Update the DATAMAX keyword */
	if (putKeyF (&image->hdr, "DATAMAX", max, ""))
	    return (status = 1);

	/* Successful return */
	return (status = 0);
}

int updMinMaxs (ShortHdrData *image) {

	/* Local variables */
	float min, max;		/* data min and max */

	/* Function definitions */
	void compMinMaxs (ShortTwoDArray *, float *, float *);

	/* Compute the min and max of the data array */
	compMinMaxs (&image->data, &min, &max);

	/* Update the DATAMIN keyword */
	if (putKeyF (&image->hdr, "DATAMIN", min, ""))
	    return (status = 1);

	/* Update the DATAMAX keyword */
	if (putKeyF (&image->hdr, "DATAMAX", max, ""))
	    return (status = 1);

	/* Successful return */
	return (status = 0);
}

void compMinMaxf (FloatTwoDArray *data, float *min, float *max) {

	/* Local variables */
	int i, j;		/* pixel indexes */
	float val;		/* pixel value */

	/* Compute the min and max of the data array */
	val = PPix(data,0,0);
	*min = val; *max = val;
	for (j=0; j<data->ny; j++) {
	     for (i=0; i<data->nx; i++) {
		  val = PPix(data,i,j);
		  if (val < *min) *min = val;
		  if (val > *max) *max = val;
	     }
	}
}

void compMinMaxs (ShortTwoDArray *data, float *min, float *max) {

	/* Local variables */
	int i, j;		/* pixel indexes */
	float val;		/* pixel value */

	/* Compute the min and max of the data array */
	val = PPix(data,0,0);
	*min = val; *max = val;
	for (j=0; j<data->ny; j++) {
	     for (i=0; i<data->nx; i++) {
		  val = PPix(data,i,j);
		  if (val < *min) *min = val;
		  if (val > *max) *max = val;
	     }
	}
}
