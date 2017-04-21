# include <stdio.h>
# include <string.h>
# include <math.h>
# include <float.h>

# include "xtables.h"
# include "ximio.h"
# include "hstio.h"

# include "stis.h"
# include "err.h"
# include "stisdef.h"
# include "calstis6.h"
# include "idtalg.h"

static int GetScatter (char *, ScatterFunctions *);
static int GetRefWave (Hdr *, char *, ScatterFunctions *);
static int GetEchelleScatter (char *, ScatterFunctions *);
static int GetHalo (char *, ScatterFunctions *, Image *, Image *, Image *);
static int GetPSF (char *, Hdr *, double, double, ScatterFunctions *,
                   Image *, Image *, Image *);
static int GetXDisp (char *, ScatterFunctions *);
static int GetRipple (Hdr *, char *, ScatterFunctions *);
static int MakeFT (Image *, Image *, ScatterFunctions *, CmplxArray *,
                   CmplxArray *);

/*static int Debug (char *, CmplxArray *);*/ /* Not used */

/*
   Read scattering functions and build kernels.

   Revision history:
   ----------------
   25 Feb 00  -  Implemented from IDL version (I.Busko)
   31 May 00  -  Adapt to new reference file formats (IB)
   17 Jul 03  -  Initialize scf->kernw[k] and psfw through k = NREFWAVE. (PEH)
    8 Jun 06  -  In MakeFT, allow for the possibility that psf is larger
                 than zpsf, i.e. either copy all of psf into the middle
                 of zpsf or copy the middle of psf into all of zpsf. (PEH)
   13 Nov 13  -  Corrected the section in GetPSF for the number of detector
                 pixels illuminated by the aperture. (PEH)
*/

int EchScatRead (Hdr *phdr, double xsize, double ysize,
                 ScatterFunctions *scf, int verbose) {

/* arguments
Hdr *phdr               i: primary header
double xsize, ysize;    i: aperture size in arcsec
ScatterFunctions *scf;  o: data structure with scattering functions
int verbose;            i: verbosity flag
*/
        char fname[STIS_LINE];          /* file names */
        Image halo1;                    /* halo images */
        Image halo2;
        Image halo3;
        Image psf1;                     /* psf images */
        Image psf2;
        Image psf3;

        int no_default = 0;             /* missing keyword is fatal error */
        int status;

        void InitImage (Image *);
        void FreeImage (Image *);

        InitImage (&halo1);
        InitImage (&halo2);
        InitImage (&halo3);
        InitImage (&psf1);
        InitImage (&psf2);
        InitImage (&psf3);

        /* Get scattering functions. */

        if ((status = Get_KeyS (phdr, "ECHSCTAB", no_default, "",
                                fname, STIS_LINE))) {
            printf (
"WARNING  The input file *might* not have the required keywords\n");
            printf (
"         to run the scattered light correction algorithm.\n");
            fflush (stdout);
            return (status);
        }
        if (verbose) {
            printf ("Reading ECHSCTAB: %s\n", fname);
            fflush (stdout);
        }
        if ((status = GetScatter (fname, scf)))
            return (status);

        /* Get echelle x-disp spread function. */

        if ((status = Get_KeyS (phdr, "EXSTAB", no_default, "",
                                fname, STIS_LINE)))
            return (status);
        if (verbose) {
            printf ("Reading EXSTAB: %s\n", fname);
            fflush (stdout);
        }
        if ((status = GetEchelleScatter (fname, scf)))
            return (status);

        /* Get x-disp function. */

        if ((status = Get_KeyS (phdr, "CDSTAB", no_default, "",
                                fname, STIS_LINE)))
            return (status);
        if (verbose) {
            printf ("Reading CDSTAB: %s\n", fname);
            fflush (stdout);
        }
        if ((status = GetXDisp (fname, scf)))
            return (status);

        /* Get echelle ripple functions. */

        if ((status = Get_KeyS (phdr, "RIPTAB", no_default, "",
                                fname, STIS_LINE)))
            return (status);
        if (verbose) {
            printf ("Reading RIPTAB: %s\n", fname);
            fflush (stdout);
        }
        if ((status = GetRipple (phdr, fname, scf)))
            return (status);

        /* Get wavelengths of halo and telescope PSF images. */

        if ((status = Get_KeyS (phdr, "SRWTAB", no_default, "",
                                fname, STIS_LINE)))
            return (status);
        if (verbose) {
            printf ("Reading SRWTAB: %s\n", fname);
            fflush (stdout);
        }

        if ((status = GetRefWave (phdr, fname, scf)))
            return (status);

        /* Get halo and telescope psf images, and build FTs */

        if ((status = Get_KeyS (phdr, "HALOTAB", no_default, "",
                                fname, STIS_LINE)))
            return (status);
        if (verbose) {
            printf ("Reading HALOTAB: %s\n", fname);
            fflush (stdout);
        }
        if ((status = GetHalo (fname, scf, &halo1, &halo2, &halo3)))
            return (status);

        if ((status = Get_KeyS (phdr, "TELTAB", no_default, "",
                                fname, STIS_LINE)))
            return (status);
        if (verbose) {
            printf ("Reading TELTAB: %s\n", fname);
            fflush (stdout);
        }
        if ((status = GetPSF (fname, phdr, xsize, ysize, scf,
                              &psf1, &psf2, &psf3)))
            return (status);

        InitCmplxArray (&(scf->ft1));
        InitCmplxArray (&(scf->ft2));
        InitCmplxArray (&(scf->ft3));
        InitCmplxArray (&(scf->fto1));
        InitCmplxArray (&(scf->fto2));
        InitCmplxArray (&(scf->fto3));

        if (verbose) {
            printf (
            "Computing Fourier Transforms of kernel at wavelength  %g\n",
                     scf->kernw[0]);
            fflush (stdout);
        }
        if ((status = MakeFT (&halo1, &psf1, scf, &(scf->ft1), &(scf->fto1))))
            return (status);

        if (scf->nwave > 1) {
            if (verbose) {
                printf (
                "Computing Fourier Transforms of kernel at wavelength  %g\n",
                         scf->kernw[1]);
                fflush (stdout);
            }
            if ((status = MakeFT (&halo2, &psf2, scf, &(scf->ft2), &(scf->fto2))))
                return (status);
        }

        if (scf->nwave > 2) {
            if (verbose) {
                printf (
                "Computing Fourier Transforms of kernel at wavelength  %g\n",
                         scf->kernw[2]);
                fflush (stdout);
            }
            if ((status = MakeFT (&halo3, &psf3, scf, &(scf->ft3), &(scf->fto3))))
                return (status);
        }

        FreeImage (&halo1);
        FreeImage (&halo2);
        FreeImage (&halo3);
        FreeImage (&psf1);
        FreeImage (&psf2);
        FreeImage (&psf3);

        return (STIS_OK);
}



/*
    Read scattering functions from ECHSCTAB reference file.
*/

static int GetScatter (char *name, ScatterFunctions *scf) {

/* arguments
char *name;             i: name of ECHSCTAB reference file
ScatterFunctions *scf;  o: data structure with scattering functions
*/
        IRAFPointer tp;
        IRAFPointer cp_optelem, cp_sporder, cp_nelem, cp_scat;
        char opt_elem[STIS_CBUF];
        int row, nrows, i;

        tp = c_tbtopn (name, IRAF_READ_ONLY, 0);
        if (c_iraferr()) {
            printf ("ERROR    ECHSCTAB `%s' not found\n", name);
            return (OPEN_FAILED);
        }
        nrows = c_tbpsta (tp, TBL_NROWS);
        c_tbcfnd1 (tp, "OPT_ELEM", &cp_optelem);
        c_tbcfnd1 (tp, "NELEM",    &cp_nelem);
        c_tbcfnd1 (tp, "SPORDER",  &cp_sporder);
        c_tbcfnd1 (tp, "ECHSCAT",  &cp_scat);
        if (cp_optelem == 0 || cp_nelem == 0 ||
            cp_sporder == 0 || cp_scat  == 0) {
            printf( "ERROR    Column not found in ECHSCTAB\n");
            c_tbtclo (tp);
            return (COLUMN_NOT_FOUND);
        }

        /* This creates some small amount of unused memory but the
           cost is negligible.
        */
        scf->scfunc = (ScFunc *) calloc (nrows , sizeof (ScFunc));
        scf->nsc = 0;

        for (row = 1; row <= nrows ; row++) {

            c_tbegtt (tp, cp_optelem, row, opt_elem, STIS_CBUF-1);

            if (streq_ic (opt_elem, scf->opt_elem)) {

                c_tbegti (tp, cp_sporder, row,
                          &(scf->scfunc[(scf->nsc)].sporder));
                if (c_iraferr())
                    return (TABLE_ERROR);
                c_tbegti (tp, cp_nelem, row, &(scf->scfunc[(scf->nsc)].nelem));
                if (c_iraferr())
                    return (TABLE_ERROR);
                scf->scfunc[(scf->nsc)].values = (double *) calloc (
                                                 scf->scfunc[(scf->nsc)].nelem,
                                                 sizeof (double));
                i = c_tbagtd (tp, cp_scat, row,
                                           scf->scfunc[(scf->nsc)].values, 1,
                                           scf->scfunc[(scf->nsc)].nelem);
                if (c_iraferr())
                    return (TABLE_ERROR);

                (scf->nsc)++;
            }
        }

        c_tbtclo (tp);

        if (scf->nsc == 0) {
            printf( "ERROR    Row not found in ECHSCTAB\n");
            return (ROW_NOT_FOUND);
        }

        return (STIS_OK);
}



/*
    Read reference wavelngths from SRWTAB reference file.
*/

static int GetRefWave (Hdr *phdr, char *name, ScatterFunctions *scf) {

/* arguments
Hdr *phdr               i: primary header
char *name;             i: name of SRWTAB reference file
ScatterFunctions *scf;  o: data structure with reference wavelengths
*/
        IRAFPointer tp;
        IRAFPointer cp_optelem, cp_cenwave, cp_nrw, cp_hrwlist, cp_prwlist;
        char opt_elem[STIS_CBUF];
        int row, nrows, i, cenwave, status;
        double holdh[10], holdp[10];   /* holds a maximum of 10 wavelengths */

        /* This is necessary to avoid rui warnings from the debugger. */
        for (i = 0; i < 10; i++) {
            holdh[i] = 0.0;
            holdp[i] = 0.0;
        }
        for (i = 0; i <= NREFWAVE; i++) {
            scf->kernw[i] = 0.0;
            scf->psfw[i] = 0.0;
        }

        /* get reference CENWAVE. */

        if ((status = Get_KeyI (phdr, "CENWAVE", 0, 0, &cenwave)))
            return (status);

        tp = c_tbtopn (name, IRAF_READ_ONLY, 0);
        if (c_iraferr()) {
            printf ("ERROR    SRWTAB `%s' not found\n", name);
            return (OPEN_FAILED);
        }
        nrows = c_tbpsta (tp, TBL_NROWS);
        c_tbcfnd1 (tp, "OPT_ELEM",  &cp_optelem);
        c_tbcfnd1 (tp, "CENWAVE",   &cp_cenwave);
        c_tbcfnd1 (tp, "NRW",       &cp_nrw);
        c_tbcfnd1 (tp, "HALOWAVES", &cp_hrwlist);
        c_tbcfnd1 (tp, "PSFWAVES",  &cp_prwlist);
        if (cp_optelem == 0 || cp_cenwave == 0 ||
            cp_nrw     == 0 || cp_hrwlist == 0 || cp_prwlist == 0) {
            printf( "ERROR    Column not found in SRWTAB\n");
            c_tbtclo (tp);
            return (COLUMN_NOT_FOUND);
        }

        scf->nwave = 0;
        for (row = 1; row <= nrows ; row++) {

            c_tbegti (tp, cp_cenwave, row, &i);
            c_tbegtt (tp, cp_optelem, row, opt_elem, STIS_CBUF-1);
            if (c_iraferr())
                return (TABLE_ERROR);
            if ((i == cenwave) && (streq_ic (opt_elem, scf->opt_elem))) {

                c_tbegti (tp, cp_nrw, row, &(scf->nwave));
                if (c_iraferr())
                    return (TABLE_ERROR);
                i = c_tbagtd (tp, cp_hrwlist, row, holdh, 1, scf->nwave);
                if (c_iraferr())
                    return (TABLE_ERROR);
                i = c_tbagtd (tp, cp_prwlist, row, holdp, 1, scf->nwave);
                if (c_iraferr())
                    return (TABLE_ERROR);
                c_tbtclo (tp);

                scf->kernw[0] = holdh[0];
                scf->psfw[0]  = holdp[0];
                if (scf->nwave > 1) {
                    scf->kernw[1] = holdh[1];
                    scf->psfw[1]  = holdp[1];
                }
                if (scf->nwave > 2) {
                    scf->kernw[2] = holdh[2];
                    scf->psfw[2]  = holdp[2];
                }

                return (STIS_OK);
            }
        }

        printf( "ERROR    Row not found in SRWTAB\n");
        c_tbtclo (tp);
        return (ROW_NOT_FOUND);
}


/*
    Read ripple functions from RIPTAB reference file.
*/

static int GetRipple (Hdr *phdr, char *name, ScatterFunctions *scf) {

/* arguments
Hdr *phdr               i: primary header
char *name;             i: name of RIPTAB  reference file
ScatterFunctions *scf;  o: data structure with ripple functions
*/
        IRAFPointer tp;
        IRAFPointer cp_optelem, cp_cenwave, cp_sporder, cp_wave,
                    cp_nelem, cp_ripple;
        char opt_elem[STIS_CBUF];
        int row, nrows, i, k, cenwave, status;

        /* get reference CENWAVE. */

        if ((status = Get_KeyI (phdr, "CENWAVE", 0, 0, &cenwave)))
            return (status);

        /* Open table and get column pointers. */

        tp = c_tbtopn (name, IRAF_READ_ONLY, 0);
        if (c_iraferr()) {
            printf ("ERROR    RIPTAB `%s' not found\n", name);
            return (OPEN_FAILED);
        }
        c_tbcfnd1 (tp, "OPT_ELEM",   &cp_optelem);
        c_tbcfnd1 (tp, "CENWAVE",    &cp_cenwave);
        c_tbcfnd1 (tp, "SPORDER",    &cp_sporder);
        c_tbcfnd1 (tp, "NELEM",      &cp_nelem);
        c_tbcfnd1 (tp, "WAVELENGTH", &cp_wave);
        c_tbcfnd1 (tp, "RIPPLE",     &cp_ripple);
        if (cp_cenwave == 0 || cp_sporder == 0 ||
            cp_wave    == 0 || cp_ripple  == 0 ||
            cp_nelem   == 0 || cp_optelem == 0) {
            printf( "ERROR    Column not found in RIPTAB\n");
            c_tbtclo (tp);
            return (COLUMN_NOT_FOUND);
        }

        /* Get # of matching rows. */

        nrows = c_tbpsta (tp, TBL_NROWS);
        scf->nrp = 0;
        for (row = 1; row <= nrows; row++) {
            c_tbegti (tp, cp_cenwave, row, &i);
            c_tbegtt (tp, cp_optelem, row, opt_elem, STIS_CBUF-1);
            if (c_iraferr())
                return (TABLE_ERROR);

            if ((i == cenwave) && (streq_ic (opt_elem, scf->opt_elem)))
                (scf->nrp)++;
        }
        if (scf->nrp == 0) {
            printf ("ERROR    No matching rows in RIPTAB `%s'\n", name);
            return (TABLE_ERROR);
        }

        /* Alloc memory. */

        scf->rpfunc = (RippleFunc *) calloc (scf->nrp, sizeof (RippleFunc));
        if (scf->rpfunc == NULL)
            return (OUT_OF_MEMORY);

        /* Ingest data from matching rows. */

        k = 0;
        for (row = 1; row <= nrows ; row++) {

            c_tbegti (tp, cp_cenwave, row, &i);
            c_tbegtt (tp, cp_optelem, row, opt_elem, STIS_CBUF-1);
            if (c_iraferr())
                return (TABLE_ERROR);

            if ((i == cenwave) && (streq_ic (opt_elem, scf->opt_elem))) {

                c_tbegti (tp, cp_sporder, row, &(scf->rpfunc[k].sporder));
                if (c_iraferr())
                    return (TABLE_ERROR);
                scf->rpfunc[k].nelem = c_tbcigi (cp_ripple,
                                                 TBL_COL_LENDATA);
                scf->rpfunc[k].wavelengths = (double *) calloc (
                                             scf->rpfunc[k].nelem,
                                             sizeof (double));
                if (scf->rpfunc[k].wavelengths == NULL)
                    return (OUT_OF_MEMORY);
                scf->rpfunc[k].values = (double *) calloc (
                                         scf->rpfunc[k].nelem,
                                         sizeof (double));
                if (scf->rpfunc[k].values == NULL)
                    return (OUT_OF_MEMORY);
                i = c_tbagtd (tp, cp_wave, row, scf->rpfunc[k].wavelengths,
                              1, scf->rpfunc[k].nelem);
                if (c_iraferr())
                    return (TABLE_ERROR);
                i = c_tbagtd (tp, cp_ripple, row, scf->rpfunc[k].values,
                              1, scf->rpfunc[k].nelem);
                if (c_iraferr())
                    return (TABLE_ERROR);

                k++;
            }
        }

        c_tbtclo (tp);

        return (STIS_OK);
}



/*
   Read echelle scattering x-dispersion spread from EXSTAB file.
*/

static int GetEchelleScatter (char *name, ScatterFunctions *scf) {

/* arguments
char *name;             i: name of EXSTAB reference file
ScatterFunctions *scf;  o: data structure with scattering functions
*/
        IRAFPointer tp;
        IRAFPointer cp_optelem, cp_nelem, cp_exscat;
        char opt_elem[STIS_CBUF];
        int row, nrows, i;

        tp = c_tbtopn (name, IRAF_READ_ONLY, 0);
        if (c_iraferr()) {
            printf ("ERROR    EXSTAB `%s' not found\n", name);
            return (OPEN_FAILED);
        }
        nrows = c_tbpsta (tp, TBL_NROWS);
        c_tbcfnd1 (tp, "OPT_ELEM", &cp_optelem);
        c_tbcfnd1 (tp, "NELEM",    &cp_nelem);
        c_tbcfnd1 (tp, "EXSCAT",   &cp_exscat);
        if (cp_optelem == 0 || cp_nelem == 0 || cp_exscat == 0) {
            printf( "ERROR    Column not found in EXSTAB\n");
            c_tbtclo (tp);
            return (COLUMN_NOT_FOUND);
        }

        for (row = 1; row <= nrows ; row++) {
            c_tbegtt (tp, cp_optelem, row, opt_elem, STIS_CBUF-1);
            if (streq_ic (opt_elem, scf->opt_elem)) {
                c_tbegti (tp, cp_nelem, row, &(scf->nspsf));
                scf->spsf = (float *) calloc (scf->nspsf, sizeof (float));
                if (scf->spsf == NULL )
                    return (OUT_OF_MEMORY);
                i = c_tbagtr (tp, cp_exscat, row, scf->spsf, 1, scf->nspsf);
                if (c_iraferr())
                    return (TABLE_ERROR);
                c_tbtclo (tp);
                return (STIS_OK);
            }
        }

        printf( "ERROR    Row with %s optical element not found in EXSTAB\n",
                 opt_elem);
        c_tbtclo (tp);
        return (ROW_NOT_FOUND);
}



/*
   Read x-disp scattering function from CDSTAB file.
*/

static int GetXDisp (char *name, ScatterFunctions *scf) {

/* arguments
char *name;             i: name of CDSTAB reference file
ScatterFunctions *scf;  o: data structure with scattering functions
*/
        IRAFPointer tp;
        IRAFPointer cp_optelem, cp_nelem, cp_cdscat;
        char opt_elem[STIS_CBUF];
        int row, nrows, i;

        tp = c_tbtopn (name, IRAF_READ_ONLY, 0);
        if (c_iraferr()) {
            printf ("ERROR    CDSTAB `%s' not found\n", name);
            return (OPEN_FAILED);
        }
        nrows = c_tbpsta (tp, TBL_NROWS);
        c_tbcfnd1 (tp, "OPT_ELEM", &cp_optelem);
        c_tbcfnd1 (tp, "NELEM",    &cp_nelem);
        c_tbcfnd1 (tp, "CDSCAT",   &cp_cdscat);
        if (cp_optelem == 0 || cp_nelem == 0 || cp_cdscat == 0) {
            printf( "ERROR    Column not found in CDSTAB\n");
            c_tbtclo (tp);
            return (COLUMN_NOT_FOUND);
        }

        for (row = 1; row <= nrows ; row++) {
            c_tbegtt (tp, cp_optelem, row, opt_elem, STIS_CBUF-1);
            if (streq_ic (opt_elem, scf->opt_elem)) {
                c_tbegti (tp, cp_nelem, row, &(scf->nxdisp));
                scf->xdisp = (float *) calloc (scf->nxdisp, sizeof (float));
                if (scf->xdisp == NULL )
                    return (OUT_OF_MEMORY);
                i = c_tbagtr (tp, cp_cdscat, row, scf->xdisp, 1, scf->nxdisp);
                if (c_iraferr())
                    return (TABLE_ERROR);
                c_tbtclo (tp);
                return (STIS_OK);
            }
        }

        printf( "ERROR    Row with %s optical element not found in CDSTAB\n",
                 opt_elem);
        c_tbtclo (tp);
        return (ROW_NOT_FOUND);
}


/*
   Read halo functions from HALOTAB file.
*/


static int GetHalo (char *name, ScatterFunctions *scf,
                    Image *halo1, Image *halo2, Image *halo3) {

/*
char *name;             i: name of HALOTAB reference file
ScatterFunctions *scf;  o: data structure with scattering functions
Image *halo1,2,3;       o: halo images, previously initialized
*/

        IRAFPointer tp;
        IRAFPointer cp_optelem, cp_refwave, cp_haldim, cp_halo;
        char opt_elem[STIS_CBUF];
        int row, nrows, i, status, k, haldim;
        double rw;

        int Alloc2DImage (Image *, int, int);

        tp = c_tbtopn (name, IRAF_READ_ONLY, 0);
        if (c_iraferr()) {
            printf ("ERROR    HALOTAB `%s' not found\n", name);
            return (OPEN_FAILED);
        }
        nrows = c_tbpsta (tp, TBL_NROWS);
        c_tbcfnd1 (tp, "OPT_ELEM", &cp_optelem);
        c_tbcfnd1 (tp, "HALOWAVE", &cp_refwave);
        c_tbcfnd1 (tp, "HALDIM",   &cp_haldim);
        c_tbcfnd1 (tp, "HALO",     &cp_halo);
        if (cp_optelem == 0 || cp_refwave == 0 ||
            cp_haldim  == 0 || cp_halo    == 0) {
            printf( "ERROR    Column not found in HALOTAB\n");
            c_tbtclo (tp);
            return (COLUMN_NOT_FOUND);
        }

        k = 0;
        for (row = 1; row <= nrows ; row++) {

            c_tbegtt (tp, cp_optelem, row, opt_elem, STIS_CBUF-1);
            c_tbegtd (tp, cp_refwave, row, &rw);
            if (c_iraferr())
                return (TABLE_ERROR);

            if ((rw == scf->kernw[k]) &&
                (streq_ic (opt_elem, scf->opt_elem))) {

                c_tbegti (tp, cp_haldim, row, &haldim);

                /* This awful code is a consequence of the
                   change in reference file format that took place
                   after the original code was completed.
                */
                switch (k) {
                case 0:
                    if ((status = Alloc2DImage (halo1, haldim, haldim)))
                        return (status);
                    i = c_tbagtr (tp, cp_halo, row, halo1->pix, 1, halo1->npix);
                    if (c_iraferr())
                        return (TABLE_ERROR);
                    break;
                case 1:
                    if ((status = Alloc2DImage (halo2, haldim, haldim)))
                        return (status);
                    i = c_tbagtr (tp, cp_halo, row, halo2->pix, 1, halo2->npix);
                    if (c_iraferr())
                        return (TABLE_ERROR);
                    break;
                case 2:
                    if ((status = Alloc2DImage (halo3, haldim, haldim)))
                        return (status);
                    i = c_tbagtr (tp, cp_halo, row, halo3->pix, 1, halo3->npix);
                    if (c_iraferr())
                        return (TABLE_ERROR);
                    break;
                default:
                    break;
                }

                k++;
            }
        }

        if (k == 0) {
            printf ("ERROR    No matching rows in HALOTAB `%s'\n", name);
            return (TABLE_ERROR);
        }

        c_tbtclo (tp);
        return (STIS_OK);
}


/*
   Read telescope PSF from TELTAB reference file.
*/

static int GetPSF (char *name, Hdr *phdr, double xsize, double ysize,
                   ScatterFunctions *scf,
                   Image *psf1, Image *psf2, Image *psf3) {

/*
char *name;             i: name of TELTAB reference file
Hdr *phdr               i: primary header
double xsize, ysize;    i: aperture size in arcsec
ScatterFunctions *scf;  o: data structure with scattering functions
Image *psf1,2,3;        o: PSF images, previously initialized
*/
        IRAFPointer tp;
        IRAFPointer cp_refwave, cp_pscale, cp_psfdim, cp_telepsf;
        int row, nrows, psfdim;
        double rw, psfscale;
        Image ospsf;
        Image imtemp, psf;
        char keyword[STIS_CBUF];
        double xplate, yplate;
        double sbox, lbox;
        double rstart, rstop;
        double frac1, frac2;
        float fhold;
        double sum;
        int i, j, k, kk, ms, ml, ns, nl, s, l;
        int ns2, nl2, s1, s2, l1, l2;
        int istart, istop;
        int status;

        int Alloc1DImage (Image *, int);
        int Alloc2DImage (Image *, int, int);
        void InitImage (Image *);
        void FreeImage (Image *);

        /* Define plate scale. */

        if (streq_ic (scf->opt_elem, "E140M"))
            xplate = 0.036;
        else if (streq_ic (scf->opt_elem, "E140H"))
            xplate = 0.0466;
        else if (streq_ic (scf->opt_elem, "E230M"))
            xplate = 0.033;
        else if (streq_ic (scf->opt_elem, "E230H"))
            xplate = 0.0466;
        else {
            printf ("Non supported grating.\n");
            return (ERROR_RETURN);
        }
        yplate = 0.029;

        /* Truncate the aperture size if it extends beyond the borders
           of an image.
        */

        /* See if uniformly illuminated aperture was used. */

        if ((status = Get_KeyS (phdr, "SCLAMP", 0, "", keyword, STIS_CBUF)))
            return (status);
        if (!streq_ic (keyword, "NONE")) {
            ml = (int)floor (ysize / yplate);
            if (ml % 2) {
                status  = Alloc1DImage (psf1, ml+1);
                status |= Alloc1DImage (psf2, ml+1);
                status |= Alloc1DImage (psf3, ml+1);
                if (status)
                    return (status);
                for (i = 1; i < ml; i++) {
                    psf1->pix[i] = 1.0 / ml;
                    psf2->pix[i] = 1.0 / ml;
                    psf3->pix[i] = 1.0 / ml;
                }
                psf1->pix[0]  = 0.5 / ml;
                psf1->pix[ml] = 0.5 / ml;
                psf2->pix[0]  = 0.5 / ml;
                psf2->pix[ml] = 0.5 / ml;
                psf3->pix[0]  = 0.5 / ml;
                psf3->pix[ml] = 0.5 / ml;
            } else {
                status  = Alloc1DImage (psf1, ml);
                status |= Alloc1DImage (psf2, ml);
                status |= Alloc1DImage (psf3, ml);
                if (status)
                    return (status);
                for (i = 0; i < ml; i++) {
                    psf1->pix[i] = 1.0 / ml;
                    psf2->pix[i] = 1.0 / ml;
                    psf3->pix[i] = 1.0 / ml;
                }
            }
            return (STIS_OK);
        }

        /* Scan reference file and process each PSF. */

        tp = c_tbtopn (name, IRAF_READ_ONLY, 0);
        if (c_iraferr()) {
            printf ("ERROR    TELTAB `%s' not found\n", name);
            return (OPEN_FAILED);
        }
        nrows = c_tbpsta (tp, TBL_NROWS);
        c_tbcfnd1 (tp, "PSFWAVE", &cp_refwave);
        c_tbcfnd1 (tp, "PSCALE",  &cp_pscale);
        c_tbcfnd1 (tp, "PSFDIM",  &cp_psfdim);
        c_tbcfnd1 (tp, "TELEPSF", &cp_telepsf);
        if (cp_refwave == 0 || cp_psfdim  == 0 ||
            cp_pscale  == 0 || cp_telepsf == 0) {
            printf( "ERROR    Column not found in TELTAB\n");
            c_tbtclo (tp);
            return (COLUMN_NOT_FOUND);
        }

        kk = 0;
        for (row = 1; row <= nrows ; row++) {

            c_tbegtd (tp, cp_refwave, row, &rw);
            if (c_iraferr())
                return (TABLE_ERROR);

            /* Matching reference wavelenghts. */

            if (rw == scf->psfw[kk]) {

                /* Read PSF image and associated data. */

                c_tbegtd (tp, cp_pscale, row, &psfscale);
                c_tbegti (tp, cp_psfdim, row, &psfdim);

                InitImage (&imtemp);
                if ((status = Alloc2DImage (&imtemp, psfdim, psfdim)))
                    return (status);
                i = c_tbagtr (tp, cp_telepsf, row, imtemp.pix, 1,
                              imtemp.npix);
                if (c_iraferr())
                    return (TABLE_ERROR);

                /* Find peak. */

                ns = imtemp.nx;
                nl = imtemp.ny;
                fhold = -FLT_MAX;
                for (i = 0; i < ns; i++) {
                    for (j = 0; j < nl; j++) {
                        if (PIX (imtemp, i, j) > fhold) {
                            s = i;
                            l = j;
                            fhold = PIX (imtemp, i, j);
                        }
                    }
                }

                /* Determine portion of PSF that made it through aperture. */

                ns2 = ((int)NINT (xsize / psfscale)) / 2;
                nl2 = ((int)NINT (ysize / psfscale)) / 2;
                s1  = ((s - ns2) >= 0) ? (s - ns2) : 0;
                s2  = ((s + ns2) <= (ns - 1)) ? (s + ns2) : (ns - 1);
                l1  = ((l - nl2) >= 0) ? (l - nl2) : 0;
                l2  = ((l + nl2) <= (nl - 1)) ? (l + nl2) : (nl - 1);

                /* Extract sub-image. */

                InitImage (&ospsf);
                if ((status = Alloc2DImage (&ospsf, s2 - s1 + 1, l2 - l1 + 1)))
                    return (status);
                k = 0;
                /* location of braces looks like a typo, but harmless */
                for (j = l1; j <= l2; j++)
                    for (i = s1; i <= s2; i++) {
                        ospsf.pix[k++] = PIX (imtemp, i, j);
                }
                FreeImage (&imtemp);
                ns = ospsf.nx;
                nl = ospsf.ny;

                /* # of STIS pixels (ms*ml) illuminated by aperture. */

                /* modified on 2013 Nov 13 by PEH
                   xplate / psfscale and yplate / psfscale are the factors
                   by which the PSF in the _tel.fits file are oversampled.
                */
                ms = (NINT (ns / (xplate / psfscale))) / 2;
                ml = (NINT (nl / (yplate / psfscale))) / 2;
                ms = 2 * ms + 1;
                ml = 2 * ml + 1;

                /* Bin oversampled PSF to STIS pixel scale. */

                /* Bin columns. */

                if ((status = Alloc2DImage (&imtemp, ms, nl)))
                    return (status);
                sbox = ns / (double)ms;
                for (i = 0; i < ms; i++) {
                    rstart = i     * sbox;
                    rstop  = (i+1) * sbox;
                    istart = (int) floor (rstart);
                    istop  = (int) floor (rstop);
                    istop = (istop < ns) ? istop : ns-1;
                    frac1 = rstart - istart;
                    frac2 = 1.0 - (rstop - istop);
                    for (j = 0; j < nl; j++) {
                        for (k = istart+1; k < istop;
                            PIX (imtemp, i, j) += PIX (ospsf, k++, j));
                        PIX (imtemp, i, j) += PIX (ospsf, istart, j) *
                                              (1.0 - frac1);
                        PIX (imtemp, i, j) += PIX (ospsf, istop,  j) *
                                              (1.0 - frac2);
                    }
                }
                FreeImage (&ospsf);

                /* Bin rows. */

                InitImage (&psf);
                if ((status = Alloc2DImage (&psf, ms, ml)))
                    return (status);
                lbox = nl / (double)ml;
                for (j = 0; j < ml; j++) {
                    rstart = j     * lbox;
                    rstop  = (j+1) * lbox;
                    istart = (int) floor (rstart);
                    istop  = (int) floor (rstop);
                    istop = (istop < nl) ? istop : nl-1;
                    frac1 = rstart - istart;
                    frac2 = 1.0 - (rstop - istop);
                    for (i = 0; i < ms; i++) {
                        for (k = istart+1; k < istop;
                             PIX (psf, i, j) += PIX (imtemp, i, k++));
                        PIX (psf, i, j) += PIX (imtemp, i, istart) *
                                            (1.0 - frac1);
                        PIX (psf, i, j) += PIX (imtemp, i, istop) *
                                            (1.0 - frac2);
                    }
                }
                FreeImage (&imtemp);

                /* Normalize PSF to unit area. */

                sum = 0.0;
                for (i = 0; i < psf.npix; sum += psf.pix[i++]);
                for (i = 0; i < psf.npix; psf.pix[i++] /= sum);

                /* This awful code is a consequence of the
                   change in reference file format that took place
                   after the original code was completed.
                */
                switch (kk) {
                case 0:
                    if ((status = Alloc2DImage (psf1, psf.nx, psf.ny)))
                        return (status);
                    for (i= 0 ; i < psf1->npix; i++)
                        psf1->pix[i] = psf.pix[i];
                    break;
                case 1:
                    if ((status = Alloc2DImage (psf2, psf.nx, psf.ny)))
                        return (status);
                    for (i= 0 ; i < psf2->npix; i++)
                        psf2->pix[i] = psf.pix[i];
                    break;
                case 2:
                    if ((status = Alloc2DImage (psf3, psf.nx, psf.ny)))
                        return (status);
                    for (i= 0 ; i < psf3->npix; i++)
                        psf3->pix[i] = psf.pix[i];
                    break;
                default:
                    break;
                }

                FreeImage (&psf);
                kk++;
            }
        }

        if (kk == 0) {
            printf ("ERROR    No matching rows in TELTAB `%s'\n", name);
            return (TABLE_ERROR);
        }

        c_tbtclo (tp);
        return (STIS_OK);
}



/*
    Build Fourier Transforms of combined scattering function.
*/

static int MakeFT (Image *iso, Image *psf, ScatterFunctions *scf,
                   CmplxArray *ft, CmplxArray *fto) {

/* arguments
Image *iso;             i: isotropic halo function (1024 x 1024)
Image *psf;             i: telescope PSF for the aperture (normalized)
ScatterFunctions *scf;  i: scattering functions
CmplxArray *ft;         o: FT
CmplxArray *fto;        o: object (central 11 X 11) FT
*/
        CmplxArray ziso, zpsf;
        float *hold;
        int i, j, k, i1, i2, j1, j2, ki, kj, kk;
        int ki1, ki2, kj1, kj2;
        int status, wsize;

        int FFTConvolve (CmplxArray *, CmplxArray *);

        /* Initialize complex arrays. */
        InitCmplxArray (&ziso);
        InitCmplxArray (&zpsf);

        /* Copy real into complex arrays. PSF must be copied into
           larger array so its size matches the halo array size.
        */

        if ((status = AllocCmplxArray (&ziso, iso->nx, iso->ny)))
            return (status);
        for (j = 1; j < iso->ny - 1; j++) {
            for (i = 1; i < iso->nx - 1; i++)
                RPIX2D (&ziso, i, j) = PPIX (iso, i-1, j-1);
        }

        if ((status = AllocCmplxArray (&zpsf, iso->nx, iso->ny)))
            return (status);
        /* The following code only works when both axis are simultaneously
           either even or odd-sized. Half-pixel shifting uses linear interp.
                (Is that comment correct?)
        */

        if (psf->nx <= iso->nx) {
            /* this is the normal case */
            i1 = (iso->nx - psf->nx) / 2;
            i2 = i1 + psf->nx;
            ki1 = 0;
            ki2 = psf->nx;
        } else {
            i1 = 0;
            i2 = iso->nx;
            ki1 = (psf->nx - iso->nx) / 2;
            ki2 = ki1 + iso->nx;
        }

        if (psf->ny <= iso->ny) {
            /* this is the normal case, i.e. a small slit */
            j1 = (iso->ny - psf->ny) / 2;
            j2 = j1 + psf->ny;
            kj1 = 0;
            kj2 = psf->ny;
        } else {
            /* long slit, greater than about 29 arcsec */
            j1 = 0;
            j2 = iso->ny;
            kj1 = (psf->ny - iso->ny) / 2;
            kj2 = kj1 + iso->ny;
        }
/* xxx ...
        if ((i1 % 2) == 0) {
... xxx */

        /* copy psf into zpsf */
        for (j = j1, kj = kj1; j < j2 && kj < kj2; j++, kj++) {
            for (i = i1, ki = ki1; i < i2 && ki < ki2; i++, ki++)
                RPIX2D (&zpsf, i, j) = PPIX (psf, ki, kj);
        }

/* xxx ...
        } else {

            for (j = j1, kj = 1; j < j2 - 1; j++, kj++) {
                for (i = i1, ki = 1; i < i2 - 1; i++, ki++) {
                    RPIX2D (&zpsf, i, j) = (PPIX (psf, ki,   kj)   +
                                            PPIX (psf, ki-1, kj)   +
                                            PPIX (psf, ki,   kj-1) +
                                            PPIX (psf, ki-1, kj-1)) / 4.0;
                }
            }

        }
... xxx */

        /* Convolve halo with telescope PSF. */

        if ((status = FFTConvolve (&zpsf, &ziso)))
            return (status);

        /* Now convolve column-by-column with x-disp scattering function. */

        for (i = 0; i < iso->nx; i++) {
            hold = (float *) calloc (iso->ny, sizeof (float));
            if (hold == NULL)
                return (OUT_OF_MEMORY);
            for (j = 0; j < iso->ny; j++) {
                kk = j - scf->nxdisp / 2;
                for (k = 0; k <scf->nxdisp; k++, kk++) {
                    kk = (kk < 0) ? 0 : kk;
                    kk = (kk >= iso->ny) ? iso->ny - 1 : kk;
                    hold[j] += RPIX2D (&zpsf, i, kk) * scf->xdisp[k];
                }
            }
            for (j = 0; j < iso->ny; j++)
                RPIX2D (&zpsf, i, j) = hold[j];
            free (hold);
        }

        /* Insert array into twice as large array and compute FT.
           No need to normalize !
        */

        if ((status = AllocCmplxArray (ft, 2 * iso->nx, 2 * iso->ny)))
            return (status);
        /* Only works when both are even-sized. */
        i1 = iso->nx / 2;
        j1 = iso->ny / 2;
        i2 = i1 + iso->nx;
        j2 = j1 + iso->ny;
        for (j = j1, kj = 0; j < j2; j++, kj++) {
            for (i = i1, ki = 0; i < i2; i++, ki++)
                RPIX2D (ft, i, j) = RPIX2D (&zpsf, ki, kj);
        }

        if ((status = fft2d (ft)))
            return (status);

        /* Now build "object" FT (central 11 X 11 portion).
           No need to normalize !
        */

        if ((status = AllocCmplxArray (fto, 2 * iso->nx, 2 * iso->ny)))
            return (status);

        wsize = 6;

        i1 = iso->nx - wsize;
        j1 = iso->ny - wsize;
        i2 = iso->nx + wsize;
        j2 = iso->ny + wsize;

        for (j = j1, kj = iso->ny / 2 - wsize; j <= j2; j++, kj++) {
            for (i = i1, ki = iso->nx / 2 - wsize; i <= i2; i++, ki++)
                RPIX2D (fto, i, j) = RPIX2D (&zpsf, ki, kj);
        }

        if ((status = fft2d (fto)))
            return (status);

        FreeCmplxArray (&ziso);
        FreeCmplxArray (&zpsf);

        return (STIS_OK);
}



/***************************************************************************/


/* This function dumps a complex array as a FITS IMSET. (Not used) */
/*
static int Debug (char *name, CmplxArray *z) {

        SingleGroup out;
        int i, j;

        initSingleGroup (&out);
        if (allocSingleGroup (&out, z->nx, z->ny) == -1)
            return (OUT_OF_MEMORY);
        for (j = 0; j < z->ny; j++) {
            for (i = 0; i < z->nx; i++) {
                Pix (out.sci.data, i, j) = RPIX2D (z, i, j);
                Pix (out.err.data, i, j) = IPIX2D (z, i, j);
            }
        }

        printf ("Writing %s image.\n", name);

        if (putSingleGroup (name, 1, &out, 0))
            return (OPEN_FAILED);

        printf ("Done writing.\n");

        return (STIS_OK);
}
*/
