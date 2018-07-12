#ifndef INCL_STISTEMPERATURE_H
#define INCL_STISTEMPERATURE_H

/* This file is used by cs1/ getgrpinfo1.c and dodark.c for scaling
   the NUV dark reference image depending on temperature, as gotten
   from the OM2CAT keyword.

   It also stores the reference CCD housing temperature used for dark 
   scaling of side 2 CCD data.

   The macro BEGIN_SIDE2 gives a date (MJD) after which any data taken
   will be using side-2 electronics.  (This is not a unique date, since
   STIS was down for almost two months.)
*/

# define CELSIUS_TO_KELVIN  273.15	/* add to Celsius to get Kelvin */

# define MIN_OM2CAT          10.	/* 10 deg C is minimim valid value */

/* DARKRATE is the nominal NUV dark count rate per second over the entire
   detector, gotten by summing counts in i9u1759io_drk.fits[1], excluding
   pixels flagged as bad.  This is used in NUVFactor in dodark.c.

   With this value of DARKRATE, the factor for scaling the NUV dark image
   will be one if the temperature is 35.5764 degrees C.
*/
# define DARKRATE 1189.3

# define CCD_REF_TEMP  18.0   /* degrees Celsius */

/* MJD for 2001 July 1
   The last side-1 data were taken before 2001 May 16, and the
   first side-2 data were taken on or around 2001 July 11.  Any
   data taken after 2001 July 1, therefore, would be side-2.
*/
# define BEGIN_SIDE2  (52091.)

#endif /* INCL_STISTEMPERATURE_H */
