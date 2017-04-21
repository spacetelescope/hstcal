# include <stdio.h>
# include "stis.h"
# include "calstis0.h"
# include "err.h"

static void StisDefaults (StisInfo *);
static int InsertSuffix (StisInfo *);

/* Initialize the calstis structure.  This includes information about the
   input and output images, calibration files, and flags to specify which
   calibration steps are to be done.  Get switches and reference file
   names from the primary headers, and check that the files exist.

   Phil Hodge, 1997 Dec 10:
	Change calling sequnces of GetSciInfo and GetWavInfo;
	call CompFiles and RefExist.

   Phil Hodge, 2001 Feb 27:
	In StisDefaults, initialize echelle and prism to false.

   Phil Hodge, 2007 May 1:
	Call checkWav to see if the wavecal has exposure time greater than
	zero and some non-zero pixel values; if not, reset wavecorr, x1dcorr,
	and x2dcorr to OMIT.
*/

int StisInit (StisInfo *sts, CalSwitch *sci_sw, CalSwitch *wav_sw,
		RefFileInfo *sciref, RefFileInfo *wavref) {

/* arguments:
StisInfo *sts      i: calibration flags and other info
CalSwitch *sci_sw  o: all calibration switches (0 or 1) for science file
CalSwitch *wav_sw  o: all calibration switches (0 or 1) for wavecal
*/

	int status;

	int missing;		/* number of missing reference files */

	int GetSciInfo (StisInfo *, CalSwitch *, RefFileInfo *);
	int GetWavInfo (StisInfo *, CalSwitch *, RefFileInfo *);
	int checkWav (StisInfo *);
	void CompFiles (RefFileInfo *, RefFileInfo *, int);
	void RefExist (RefFileInfo *, int, int *);
	int InitNames (void);
	void FreeNames (void);

	/* Assign default values in sts. */
	StisDefaults (sts);

	/* Construct output and temporary file names. */
	if ((status = InsertSuffix (sts)) != 0)
	    return (status);

	if ((status = InitNames()) != 0) /* no reference file names yet */
	    return (status);

	/* Get switches and reference file names for science file.
	   Also get the wavecal file name, if that wasn't specified
	   on the command line.
	*/
	if ((status = GetSciInfo (sts, sci_sw, sciref)) != 0)
	    return (status);

	/* Get switches and reference file names for wavecal, if any. */
	if (sts->sci_wavecorr == PERFORM) {
	    if (sts->wavfile[0] == '\0') {
		printf ("ERROR    WAVECORR = PERFORM, but WAVECAL is blank.\n");
		return (OPEN_FAILED);
	    }
	    PrFileName ("wavecal", sts->wavfile);	/* print wavecal name */
	    if ((status = GetWavInfo (sts, wav_sw, wavref)) != 0)
		return (status);
	    /* reset wavecorr, x1dcorr, x2dcorr if wavecal exptime is zero */
	    if ((status = checkWav (sts)) != 0)
		return (status);
	} else if (GotFileName (sts->wavfile)) {
	    printf (
	"Warning  WAVECAL was specified, but WAVECORR is not PERFORM.\n");
	}

	/* Compare the lists of reference files for science and wavecal. */
	if (sts->sci_wavecorr == PERFORM)
	    CompFiles (sciref, wavref, sts->samebin);

	/* Verify that the reference files that we need do exist. */
	missing = 0;			/* initial value */
	RefExist (sciref, sts->detector, &missing);	/* for science file */
	if (sts->sci_wavecorr == PERFORM)
	    RefExist (wavref, sts->detector, &missing);	/* wavecal */

	FreeNames();			/* done with ref file names */

	if (missing > 0) {
	    if (missing == 1)
		printf ("ERROR    One reference file was missing.\n");
	    else
		printf ("ERROR    %d reference files were missing.\n", missing);
	    return (CAL_FILE_MISSING);
	}

	return (0);
}

static void StisDefaults (StisInfo *sts) {

			/* rawfile has been assigned already */
			/* wavfile may have been assigned already */
	sts->crjfile[0] = '\0';
	sts->fltfile[0] = '\0';
	sts->x1dfile[0] = '\0';
	sts->x2dfile[0] = '\0';
	sts->sx2file[0] = '\0';
	sts->sx1file[0] = '\0';
	sts->sflfile[0] = '\0';
	sts->blv_tmp[0] = '\0';
	sts->crj_tmp[0] = '\0';
	sts->fwv_tmp[0] = '\0';
	sts->cwv_tmp[0] = '\0';
	sts->w2d_tmp[0] = '\0';
	sts->rootname[0] = '\0';
	sts->detector = UNKNOWN_DETECTOR;
	sts->obstype = UNKNOWN_TYPE;
	sts->nimages = 1;
	sts->echelle = 0;
	sts->prism = 0;

	/* Initialize flags to not perform the step. */
	sts->sci_basic_2d_a = OMIT;
	sts->sci_basic_2d = OMIT;
	sts->sci_expscorr = OMIT;
	sts->sci_crcorr = OMIT;
	sts->sci_rptcorr = OMIT;
	sts->sci_wavecorr = OMIT;
	sts->sci_2d_rect = OMIT;
	sts->sci_1d_extract = OMIT;
	sts->sci_geocorr = OMIT;
	sts->wav_basic_2d = OMIT;
	sts->wav_subscicorr = OMIT;
}

/* Construct output and temporary file names from outroot.

   Note that we do not construct the wavecal (_wav) file name.  This is
   because it is either given on the command line or is gotten as the
   value of the WAVECAL keyword.
*/

static int InsertSuffix (StisInfo *sts) {

	int status;

	int MkName (char *, char *, char *, char *, int);

	if ((status = MkName (sts->outroot, "_raw", "_crj", sts->crjfile,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_flt", sts->fltfile,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_x1d", sts->x1dfile,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_x2d", sts->x2dfile,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_sx2", sts->sx2file,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_sx1", sts->sx1file,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_sfl", sts->sflfile,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_blv_tmp", sts->blv_tmp,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_crj_tmp", sts->crj_tmp,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_fwv_tmp", sts->fwv_tmp,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_cwv_tmp", sts->cwv_tmp,
		STIS_LINE)) != 0)
	    return (status);

	if ((status = MkName (sts->outroot, "_raw", "_w2d_tmp", sts->w2d_tmp,
		STIS_LINE)) != 0)
	    return (status);

	return (0);
}
