# define DAYS_PER_YEAR  365.25

/* Structure that stores time-dependent-sensitivity information. It is
   populated from the time-dependent sensitivity table (_tds).
   ref_temp and temp_sens are used for the temperature-dependent
   sensitivity correction.
*/

typedef struct {
	int allocated;		/* true if memory has been allocated */
	int nwl;		/* number of wavelengths */
	int nt;			/* number of times (linear segments) */
	double *wl;		/* array of wavelengths */
	double *time;		/* array of times (MJD) */
	double **slope;		/* array of slopes */
	double ref_temp;	/* reference temperature, degrees C */
	double *temp_sens;	/* array of temperature sensitivity values */
} TdsInfo;
