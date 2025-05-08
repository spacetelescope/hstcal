#ifndef INCL_ACSVERSION_H
#define INCL_ACSVERSION_H

/* This string is written to the output primary header as CAL_VER. */
#define ACS_CAL_VER "10.5.0 (08-May-2025)"
#define ACS_CAL_VER_NUM "10.5.0"

/* This Generation of the CTE algorithm is obsolete.  These strings
   are only maintained in case a user has an older version of the
   PCTETAB. */
/* name and version number of the CTE correction algorithm */
#define ACS_GEN1_CTE_NAME "PixelCTE 2012"
#define ACS_GEN1_CTE_VER "3.3"

/* This Generation of the CTE algorithm is obsolete.  These strings
   are only maintained in case a user has an older version of the
   PCTETAB. */
/* name and version number of generation 2 CTE correction algorithm */
#define ACS_GEN2_CTE_NAME "PixelCTE 2017"
#define ACS_GEN2_CTE_VER "2.0"

/* name and version number of generation 3 CTE correction algorithm */
#define ACS_GEN3_CTE_NAME "Par/Serial PixelCTE 2023"
#define ACS_GEN3_CTE_VER "3.0"

#endif /* INCL_ACSVERSION_H */
