# include <stdio.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "stisvoltages.h"

/* Value to assign if keyword is not found. */
# define DEFAULT_VOLTAGE  (-9999.)

/* This routine gets voltages from the SCI extension header and compares with
   minimum values.  If any voltage is too low, a warning message will be
   printed.

   Phil Hodge, 1998 May 1:
	Function created.

   Phil Hodge, 2001, Aug 27:
	Move this file from cs0 to this directory.  Also, use getKeyD
	instead of Get_KeyD, and print a warning if any of the keywords
	is missing from the header.
*/

int CheckVolt (Hdr *hdr) {

/* argument:
Hdr *hdr         i: SCI extension header
*/

	int status;

	double abav, cdav;	/* amplifier bias voltages */
	double lgcdv;		/* last gate C & D voltage */
	double swalv;		/* summing wells */
	double rcdlv;		/* reset gate C & D low voltage */
	int low_voltage = 0;	/* was any voltage low? */

	/* Get the voltages. */
	status = getKeyD (hdr, "OCBABAV",&abav);
	if (status) {
	    printf ("Warning  Keyword OCBABAV not found.\n");
	    abav = DEFAULT_VOLTAGE;
	}
	status = getKeyD (hdr, "OCBCDAV",&cdav);
	if (status) {
	    printf ("Warning  Keyword OCBCDAV not found.\n");
	    cdav = DEFAULT_VOLTAGE;
	}
	status = getKeyD (hdr, "OCBLGCDV",&lgcdv);
	if (status) {
	    printf ("Warning  Keyword OCBLGCDV not found.\n");
	    lgcdv = DEFAULT_VOLTAGE;
	}
	status = getKeyD (hdr, "OCBSWALV",&swalv);
	if (status) {
	    printf ("Warning  Keyword OCBSWALV not found.\n");
	    swalv = DEFAULT_VOLTAGE;
	}
	status = getKeyD (hdr, "OCBRCDLV",&rcdlv);
	if (status) {
	    printf ("Warning  Keyword OCBRCDLV not found.\n");
	    rcdlv = DEFAULT_VOLTAGE;
	}

	/* Check the voltages against the minimum allowed values. */
	if (abav != DEFAULT_VOLTAGE && abav < MIN_OCBABAV) {
	    low_voltage = 1;
	    printf ("Warning  OCBABAV = %.6g\n", abav);
	}
	if (cdav != DEFAULT_VOLTAGE && cdav < MIN_OCBCDAV) {
	    low_voltage = 1;
	    printf ("Warning  OCBCDAV = %.6g\n", cdav);
	}
	if (lgcdv != DEFAULT_VOLTAGE && lgcdv < MIN_OCBLGCDV) {
	    low_voltage = 1;
	    printf ("Warning  OCBLGCDV = %.6g\n", lgcdv);
	}
	if (swalv != DEFAULT_VOLTAGE && swalv < MIN_OCBSWALV) {
	    low_voltage = 1;
	    printf ("Warning  OCBSWALV = %.6g\n", swalv);
	}
	if (rcdlv != DEFAULT_VOLTAGE && rcdlv < MIN_OCBRCDLV) {
	    low_voltage = 1;
	    printf ("Warning  OCBRCDLV = %.6g\n", rcdlv);
	}

	if (low_voltage) {
	    printf (
	"Warning  Data are likely to be corrupted; CCD settings indicate \\\n");
	    printf (
	"Warning  a cosmic ray induced CCD reset has occured, causing \\\n");
	    printf (
	"Warning  biases and clocking voltages to be incorrectly set.\n");
	}

	return (0);
}
