#ifndef INCL_STISTDS_H
#define INCL_STISTDS_H

# define DAYS_PER_YEAR  365.25

/* Structure that stores time-dependent-sensitivity information. It is
   populated from the time-dependent sensitivity table (_tds).
   ref_temp and temp_sens are used for the temperature-dependent
   sensitivity correction.
*/

# define ORIGINAL_TDS_FORMAT 0
# define COS_TDS_FORMAT      1

typedef struct {
	int allocated;		/* true if memory has been allocated */
	int format;		/* specifies original STIS or COS-like TDS */
	int nwl;		/* number of wavelengths */
	int nt;			/* number of times (linear segments) */
	double ref_time;	/* reference time (MJD), from table header */
	double *wl;		/* array of wavelengths */
	double *time;		/* array of times (MJD) */
	double **slope;		/* array of slopes */
	double **intercept;	/* array of intercepts */
	double ref_temp;	/* reference temperature, degrees C */
	double *temp_sens;	/* array of temperature sensitivity values */
} TdsInfo;

#endif /* INCL_STISTDS_H */
