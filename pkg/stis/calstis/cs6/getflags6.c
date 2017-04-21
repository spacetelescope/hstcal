# include <stdio.h>
# include <string.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "err.h"

static int CheckX1D (Hdr *, StisInfo6 *, int *);
static int CheckOptimal (Hdr *, StisInfo6 *);
static int CheckBack (Hdr *, StisInfo6 *);
static int CheckDisp (Hdr *, StisInfo6 *, int *);
static int CheckFlux (Hdr *, StisInfo6 *, int *);
static int CheckHelio (Hdr *, StisInfo6 *);
static int CheckSmGeo (Hdr *, StisInfo6 *, int *);
static int CheckWave (Hdr *, StisInfo6 *);

/*
   Get calibration flag values and names of reference images and tables
   from the primary header. Check the existence of the reference files
   and get pedigree and descrip, and also check for possible inconsistency
   between switches.

   The calibration switches are read from the header only if calstis6
   is being run from the pipeline. In standalone mode the switches must
   be supplied from the command line.

   The optimal extraction-related file names may be not present in the
   header; this does not prevent the program to continue. The file
   names can be provided in the command line, or the program could
   use unweighted extraction instead.

   On 25Apr01 removed all switch reading code. This code is not necessary
   since switches are always input into calstis6 via argument list.
   Reading code here actually prevented switches from being set by the
   pipeline. This caused problems when calstis6std was run from wihin the
   sc2d algorithm.



   Revision history:
   ----------------
   20 Feb 97  -  Adapted from similar routine in calstis7 (I.Busko)
   25 Apr 97  -  Removed ESPTAB (IB)
   08 May 97  -  Conform to new _trl standard (IB)
   14 May 97  -  Flag signals CalStis6 is being run from the pipeline (IB)
   14 Nov 97  -  Get SDCTAB (IB)
   27 Jan 98  -  Get PCTAB (IB)
   13 Apr 98  -  Removed restriction on flux calibration (IB)
   20 Jul 98  -  APDESTAB is always read (A2 offset is always needed, IB)
   15 Sep 98  -  Check optimal extraction-related reference files IB)
   08 Nov 00  -  Do not check for MOFFTAB (IB)
   25 Apr 01  -  Turn off switch input (IB)
   08 Jan 02  -  Check for TDSTAB keyword (IB)
   17 Jun 03  -  Check for CCDTAB keyword (PB)
   08 Apr 05  -  Check for GACTAB (PEH)
*/

int GetFlags6 (StisInfo6 *sts, Hdr *phdr) {

        int status;

        int missing = 0;        /* true if any calibration file is missing */

        /* Check main reference files. */

        if ((status = CheckX1D (phdr, sts, &missing)))
            return (status);

        if ((status = CheckOptimal (phdr, sts)))
            return (status);

        /* Check remaining calibration switches and reference files. */

        if ((status = CheckBack (phdr, sts)))
            return (status);

        if ((status = CheckDisp (phdr, sts, &missing)))
            return (status);

        if ((status = CheckHelio (phdr, sts)))
            return (status);

        if (sts->dispcorr != PERFORM && sts->heliocorr == PERFORM) {
            printf (
        "ERROR    No wavelengths - cannot apply heliocentric correction.\n");
            return (ERROR_RETURN);
        }

        if ((status = CheckFlux (phdr, sts, &missing)))
            return (status);

        if (sts->fluxcorr == PERFORM &&
            sts->dispcorr != PERFORM) {
            printf (
           "ERROR    No wavelengths - cannot flux-calibrate.\n");
            return (ERROR_RETURN);
        }

        if (sts->detector != CCD_DETECTOR) {
            if ((status = CheckSmGeo (phdr, sts, &missing)))
                return (status);
        }

        if ((status = CheckWave (phdr, sts)))             /* just check flag */
            return (status);

        if (missing)
            return (CAL_FILE_MISSING);
        else
            return (0);
}



/* Check X1DCORR reference files. */

static int CheckX1D (Hdr *phdr, StisInfo6 *sts, int *missing) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo6 *sts  i: switches, file names, etc
int *missing    io: incremented if the file is missing
*/

        int status;

        int GetCheckRef (Hdr *, char *, RefTab *, int *, int *, int);

        /* Spectrum trace table. */
        if ((status = GetCheckRef (phdr, "SPTRCTAB", &sts->sptrctab,
                                   &sts->x1d, missing, FATAL)))
            return (status);

        /* SDC table. */
        if ((status = GetCheckRef (phdr, "SDCTAB", &sts->distntab,
                                   &sts->x1d, missing, FATAL)))
            return (status);

        /* 1-D extraction parameters table. */
        if ((status = GetCheckRef (phdr, "XTRACTAB", &sts->xtrctab,
                                   &sts->x1d, missing, FATAL)))
            return (status);

        /* Aperture description table. */
        if ((status = GetCheckRef (phdr, "APDESTAB", &sts->apdestab,
                                   &sts->dispcorr, missing, FATAL)))
            return (status);

        /* CCD description table. */
        if (sts->detector == CCD_DETECTOR) {
            if (sts->ctecorr == PERFORM) {
                if ((status = GetCheckRef (phdr, "CCDTAB", &sts->ccdtab,
                                  &sts->x1d, missing, FATAL)) != 0) {
                    return (status);
                }
            }
        } else {
            sts->ctecorr = OMIT;
        }

        if (sts->x1d != PERFORM) {
            printf (
            "Warning  X1DCORR skipped due to dummy reference file.\n");
            return (NOTHING_TO_DO);
        }
        return (0);
}



/* Check optimal extraction-related reference files. Notice that the
   "missing" global counter is replaced by a local counter. It is no
   error if the optimal extraction files are missing in the header.
*/

static int CheckOptimal (Hdr *phdr, StisInfo6 *sts) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo6 *sts  i: switches, file names, etc
int *missing    io: incremented if the file is missing
*/

        int status;
        int missing;

        int GetCheckRef (Hdr *, char *, RefTab *, int *, int *, int);

        missing = 0;

        /* Profile table. */
        if ((status = GetCheckRef (phdr, "OPROFTAB",
                                   &sts->pftab, &sts->x1d, &missing, NO_FATAL)))
            return (status);

        /* Flux table. */
        if ((status = GetCheckRef (phdr, "OSPECTAB",
                                   &sts->pxtab, &sts->x1d, &missing, NO_FATAL)))
            return (status);

/*
        if (missing != 0 && streq_ic (sts->xtracalg, OPTIMAL))
            printf (
"Warning  No reference file names for optimal extraction exist in header.\n");
*/


        return (0);
}



/* Check whether we should extract and subtract the background. No need
   to check for the spectrum extraction table, since this step can be
   performed only if the spectrum itself is being extracted, and both
   steps share the same reference table.
*/

static int CheckBack (Hdr *phdr, StisInfo6 *sts) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo6 *sts  i: switches, file names, etc
*/
        /*
        if (sts->pipeline) {
            if ((status = GetSwitch (phdr, "BACKCORR", &sts->backcorr)))
                return (status);
        }
        */

        return (0);
}



/* Check whether we should generate the wavelength array.  If so,
   check for the existence of the necessary reference tables.
*/

static int CheckDisp (Hdr *phdr, StisInfo6 *sts, int *missing) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo6 *sts  i: switches, file names, etc
int *missing    io: incremented if the table is missing
*/

        int status;

        int GetCheckRef (Hdr *, char *, RefTab *, int *, int *, int);

        /*
        if (sts->pipeline) {
            if ((status = GetSwitch (phdr, "DISPCORR", &sts->dispcorr)))
                return (status);
        }
        */

        if (sts->dispcorr == PERFORM) {

            /* Dispersion coefficients. */
            if ((status = GetCheckRef (phdr, "DISPTAB", &sts->disptab,
                                       &sts->dispcorr, missing, FATAL)))
                return (status);

            /* Incidence-angle correction table. */
            if ((status = GetCheckRef (phdr, "INANGTAB", &sts->inangtab,
                                       &sts->dispcorr, missing, FATAL)))
                return (status);
        }

        return (0);
}



/* Check whether we should convert to absolute flux units.  If so,
   check for the existence of the absolute sensitivity table and
   aperture throughput table.
   Also check for the photometric correction table (pctab), but
   if that table is missing, the pctcorr switch will be reset to
   indicate that that particular correction should not be included.
*/

static int CheckFlux (Hdr *phdr, StisInfo6 *sts, int *missing) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo6 *sts  i: switches, file names, etc
int *missing    io: incremented if the table is missing
*/

        int status;
        int l_missing;          /* local variable for PCTAB */

        int GetCheckRef (Hdr *, char *, RefTab *, int *, int *, int);

        /*
        if (sts->pipeline) {
            if ((status = GetSwitch (phdr, "FLUXCORR", &sts->fluxcorr)))
                return (status);
        }
        */

        if (sts->fluxcorr == PERFORM) {

            /* Throughput table. */
            if ((status = GetCheckRef (phdr, "PHOTTAB", &sts->phottab,
                                       &sts->fluxcorr, missing, FATAL)))
                return (status);

            /* Relative aperture throughput table. */
            if ((status = GetCheckRef (phdr, "APERTAB", &sts->apertab,
                                       &sts->fluxcorr, missing, FATAL)))
                return (status);

            /* Table for corrections for slit height.  Note that we pass
                pctcorr instead of fluxcorr, and l_missing is local.
            */
            l_missing = 0;
            if ((status = GetCheckRef (phdr, "PCTAB", &sts->pctab,
                                       &sts->pctcorr, &l_missing, NO_FATAL)))
                return (status);
            if (l_missing > 0)
                sts->pctcorr = OMIT;

            /* Table for grating-aperture correction.  Note that we pass
                gaccorr instead of fluxcorr, and l_missing is local.
            */
            l_missing = 0;
            if ((status = GetCheckRef (phdr, "GACTAB", &sts->gactab,
                                       &sts->gaccorr, &l_missing, NO_FATAL)))
                return (status);
            if (l_missing > 0) {
                sts->gaccorr = OMIT;
                printf (
        "Warning  Grating-aperture throughput correction table (GACTAB)"
        " was not found,\n");
                printf ("         and no gac corrections will be applied\n");
            }

            /* Time-dependent sensitivity table. Handled like pct. */
            l_missing = 0;
            if ((status = GetCheckRef (phdr, "TDSTAB", &sts->tdstab,
                                       &sts->tdscorr, &l_missing, NO_FATAL)))
                return (status);
            sts->tdscorr = PERFORM;
            if (l_missing > 0)
                sts->tdscorr = OMIT;
        }

        return (0);
}



/* Check whether this step is to be performed.  There is no associated
   reference file.
*/

static int CheckHelio (Hdr *phdr, StisInfo6 *sts) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo6 *sts  i: switches, file names, etc
*/
    /*
        int status;


        if (sts->pipeline) {
            if ((status = GetSwitch (phdr, "HELCORR", &sts->heliocorr)))
                return (status);
        }
    */

        return (0);
}



/* If this step is to be performed, check for the existence of the
   small-scale distortion file sdstfile, and get pedigree and descrip.
*/

static int CheckSmGeo (Hdr *phdr, StisInfo6 *sts, int *missing) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo6 *sts  i: switches, file names, etc
int *missing    io: incremented if the file is missing
*/

        int status;

        int use_def = 0;                /* missing keyword is fatal error */

        /*
        if (sts->pipeline) {
            if ((status = GetSwitch (phdr, "SGEOCORR", &sts->sgeocorr)))
                return (status);
        }
        */

        /* Are we supposed to do this step? */
        if (sts->sgeocorr == PERFORM) {

            /* Small-scale distortion file. */
            if ((status = Get_KeyS (phdr, "SDSTFILE", use_def, "",
                                    sts->sdstfile.name, STIS_LINE)))
                return (status);

            /* Open the image to verify that it exists, and if it does,
                get pedigree & descrip.
            */
            if ((status = ImgPedigree (&sts->sdstfile)))
                return (status);
            if (sts->sdstfile.exists != EXISTS_YES) {
                (*missing)++;
                printf ("ERROR    SDSTFILE `%s' not found\n",
                        sts->sdstfile.name);
            }
            if (sts->sdstfile.goodPedigree != GOOD_PEDIGREE)
                sts->sgeocorr = DUMMY;
        }

        return (0);
}



/* Check whether the wavecal has been used to update the coordinates. */

static int CheckWave (Hdr *phdr, StisInfo6 *sts) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo6 *sts  i: switches, file names, etc
*/

        /*
        if (sts->pipeline) {
            if ((status = GetSwitch (phdr, "WAVECORR", &sts->wavecorr)))
                return (status);
        }
        */

        return (0);
}
