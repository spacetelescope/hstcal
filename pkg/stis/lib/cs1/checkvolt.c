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
	    trlwarn("Keyword OCBABAV not found.");
	    abav = DEFAULT_VOLTAGE;
	}
	status = getKeyD (hdr, "OCBCDAV",&cdav);
	if (status) {
	    trlwarn("Keyword OCBCDAV not found.");
	    cdav = DEFAULT_VOLTAGE;
	}
	status = getKeyD (hdr, "OCBLGCDV",&lgcdv);
	if (status) {
	    trlwarn("Keyword OCBLGCDV not found.");
	    lgcdv = DEFAULT_VOLTAGE;
	}
	status = getKeyD (hdr, "OCBSWALV",&swalv);
	if (status) {
	    trlwarn("Keyword OCBSWALV not found.");
	    swalv = DEFAULT_VOLTAGE;
	}
	status = getKeyD (hdr, "OCBRCDLV",&rcdlv);
	if (status) {
	    trlwarn("Keyword OCBRCDLV not found.");
	    rcdlv = DEFAULT_VOLTAGE;
	}

	/* Check the voltages against the minimum allowed values. */
	if (abav != DEFAULT_VOLTAGE && abav < MIN_OCBABAV) {
	    low_voltage = 1;
	    trlwarn("OCBABAV = %.6g", abav);
	}
	if (cdav != DEFAULT_VOLTAGE && cdav < MIN_OCBCDAV) {
	    low_voltage = 1;
	    trlwarn("OCBCDAV = %.6g", cdav);
	}
	if (lgcdv != DEFAULT_VOLTAGE && lgcdv < MIN_OCBLGCDV) {
	    low_voltage = 1;
	    trlwarn("OCBLGCDV = %.6g", lgcdv);
	}
	if (swalv != DEFAULT_VOLTAGE && swalv < MIN_OCBSWALV) {
	    low_voltage = 1;
	    trlwarn("OCBSWALV = %.6g", swalv);
	}
	if (rcdlv != DEFAULT_VOLTAGE && rcdlv < MIN_OCBRCDLV) {
	    low_voltage = 1;
	    trlwarn("OCBRCDLV = %.6g", rcdlv);
	}

	if (low_voltage) {
	    trlwarn("Data are likely to be corrupted; CCD settings indicate ");
	    trlwarn("a cosmic ray induced CCD reset has occured, causing ");
	    trlwarn("biases and clocking voltages to be incorrectly set.");
	}

	return (0);
}
