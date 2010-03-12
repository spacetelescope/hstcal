/*
** M. Droettboom, January 2010:
** Change to use CFITSIO rather than IRAF IMIO routines.
*/

#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>

/* CFITSIO TODO: store axes in IODesc object as array to more
   conveniently interface with CFITSIO */

static void detect_iraferr(void) {
        sprintf(error_msg,"\nIRAF error %d: %s\n",c_iraferr(),
                c_iraferrmsg());
}

static void ioerr(HSTIOError e, IODescPtr x_, int status) {
        IODesc *x;
        char cfitsio_errmsg[80];
        x = (IODesc *)x_;
        sprintf(&error_msg[strlen(error_msg)],
                "Filename %s EXTNAME %s EXTVER %d CFITSIO status %d\n",
                x->filename, x->extname, x->extver, status);
        while (fits_read_errmsg(cfitsio_errmsg)) {
            strncat(error_msg, cfitsio_errmsg, 80);
        }
        error(e,0);
}

/*
** Make_iodesc takes a filename, extname, and extver and creates and
** initializes an IODesc structure.  In the process, it builds a
** correct filename to be used in the open statement to IRAF.  This
** constructed filename is returned.
*/
static char *make_iodesc(IODesc **x, char *fname, char *ename, int ever) {
        int i, n, flen;
        char *tmp;
        IODesc *iodesc;
        char xname[9];

        iodesc = (IODesc *)calloc(1,sizeof(IODesc));
        if (iodesc == NULL) {
            error(NOMEM,"Allocating I/O descriptor");
            return NULL;
        }
        iodesc->ff = NULL;
        iodesc->filename = NULL;
        iodesc->extname = NULL;
        iodesc->extver = 0;
        iodesc->hflag = 0;
        iodesc->hdr = NULL;
        iodesc->dims[0] = 0;
        iodesc->dims[1] = 0;
        iodesc->type = 0;
        if (fname == 0) fname = "";
        if (ename == 0) ename = "";
        iodesc->filename = (char *)calloc(((flen = strlen(fname)) + 1),
                sizeof(char));
        if (iodesc->filename == NULL) {
            free(iodesc);
            error(NOMEM,"Allocating I/O descriptor");
            return NULL;
        }
        n = strlen(ename);
        if (n > 8) { ioerr(BADEXTNAME,iodesc,0); return NULL; }
        for (i = 0; i < n; ++i)
            xname[i] = toupper(ename[i]);
        for (--i; i >= 0 && xname[i] == ' '; --i) ;
        ++i;
        xname[i] = '\0';
        iodesc->extname = (char *)calloc((strlen(xname) + 1),sizeof(char));
        if (iodesc->extname == NULL) {
            free(iodesc->filename);
            free(iodesc);
            error(NOMEM,"Allocating I/O descriptor");
            return NULL;
        }
        strcpy(iodesc->filename,fname);
        strcpy(iodesc->extname,xname);
        iodesc->extver = ever;

        /* make up the proper filename */
        /* check for a request for the primary HDU */
        tmp = (char *)calloc((flen + 80),sizeof(char));
        if (tmp == NULL) { error(NOMEM,"Allocating I/O descriptor"); return NULL; }
        strcpy(tmp,fname);
        if (ever == 0 || ename == 0 || ename[0] == '\0' || ename[0] == ' ')
            strcat(tmp,"[0]");
        else
            sprintf(&tmp[flen],"[%s,%d]",xname,ever);

        *x = iodesc;
        return tmp;
}

IODescPtr openInputImage(char *fname, char *ename, int ever) {
        IODesc *iodesc;
        int no_dims;
        char *tmp;
        char ospath[SZ_PATHNAME];
        int open_mode;
        int status = 0;

        /* CFITSIO: Error handling */
        c_pusherr(detect_iraferr);

        tmp = make_iodesc(&iodesc, fname, ename, ever);
        if (tmp == NULL) return NULL;
        iodesc->options = ReadOnly;

        /* CFITSIO: Resolve this inheritance stuff */
        /* p = strstr(tmp,"[0]"); */
        /* if (p == NULL) { */
        /*     tmp[strlen(tmp) - 1] = '\0'; */
        /*     strcat(tmp, ",NOINHERIT]"); */
        /* } */

        /* open the file using CFITSIO */
        if (c_vfn2osfn(fname, ospath)) {
            free(tmp);
            return NULL;
        }

        /*
        This is to get around a safety feature of CFITSIO.  If a file
        is opened READONLY, and then later opened READWRITE (without
        closing it first), CFITSIO raises an error, since it assumes
        that once a file is READONLY, it is always READONLY.  We can't
        just always open READWRITE here, since that will fail if the
        permissions on the file are not writable.  So, we check to see
        if the file is writable, and then open accordingly.
         */
        if (access(ospath, W_OK)) {
            open_mode = READONLY;
        } else {
            open_mode = READWRITE;
        }

        if (c_vfn2osfn(tmp, ospath)) {
            free(tmp);
            return NULL;
        }
        free(tmp);

        if (fits_open_file(&iodesc->ff, ospath, open_mode, &status)) {
            ioerr(BADOPEN, iodesc, status);
            free(iodesc->extname);
            free(iodesc->filename);
            free(iodesc);
            return NULL;
        }

        /* get the dimensions and type */
        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        fits_get_img_equivtype(iodesc->ff, &iodesc->type, &status);
        fits_get_img_size(iodesc->ff, 2, iodesc->dims, &status);
        if (status) {
            ioerr(BADDIMS, iodesc, status);
            return NULL;
        }
        if (no_dims == 2) {
            /* Nothing */
        } else if (no_dims == 1) {
            iodesc->dims[1] = 0;
        } else if (no_dims == 0) {
            iodesc->dims[0] = 0;
            iodesc->dims[1] = 0;
        } else {
            ioerr(BADDIMS, iodesc, 0);
            return NULL;
        }

        clear_err();

        return iodesc;
}

IODescPtr openOutputImage(char *fname, char *ename, int ever, Hdr *hd,
        int d1, int d2, FitsDataType typ) {
        IODesc *iodesc;
        char *tmp;
        char date[12];
        char date_card[80];
        time_t t;
        struct tm *time_tmp;
        FitsKw kw;
        char ename_val[9];
        int ever_val;
        char ospath[SZ_PATHNAME];
        int status = 0;

        /* CFITSIO: Error handling */
        c_pusherr(detect_iraferr);

        tmp = make_iodesc(&iodesc, fname, ename, ever);
        if (tmp == NULL) return NULL;
        iodesc->options = WriteOnly;
        if (ever == 0 || ename == 0 || ename[0] == '\0' || ename[0] == ' ') {
            int rtn = ckNewFile(fname);
            if (rtn == 1) {
                ioerr(BADEXIST, iodesc, 0);
                return NULL;
            } else if (rtn == 2) {
                ioerr(BADREMOVE, iodesc, 0);
                return NULL;
            }
        }

        /* CFITSIO: Check this INHERIT, APPEND nonsense works */
        /* p = strstr(tmp, "[0]"); */
        /* if (p == NULL) { */
        /*     tmp[strlen(tmp) - 1] = '\0'; */
        /*     strcat(tmp, ",INHERIT,APPEND]"); */
        /* } */
        /* else { */
        /*     tmp[strlen(tmp) - 3] = '\0'; /\* eliminate the "[0]" *\/ */
        /* } */

        /* make sure ename and ever are in the header array */
        kw = findKw(hd, "EXTNAME");
        if (kw == NotFound) {
            if (ever != 0 && ename != 0 &&
                ename[0] != '\0' && ename[0] != ' ') {
                kw = insertfirst(hd);
                kw = insertStringKw(kw, "EXTNAME", ename, "Name of the extension");
            }
        } else {
            /* Make sure it has the right value */
            getStringKw(kw,ename_val,8);
            if (strncpy(ename_val, ename, strlen(ename)) != 0)
                putStringKw(kw,ename);
        }
        kw = findKw(hd,"EXTVER");
        if (kw == NotFound) {
            if (ever != 0 && ename != 0 &&
                ename[0] != '\0' && ename[0] != ' ') {
                kw = findKw(hd, "EXTNAME");
                kw = insertIntKw(kw, "EXTVER", ever, "Extension version");
            }
        } else {
            /* Make sure it has the right value */
            ever_val = getIntKw(kw);
            if (ever != ever_val)
                putIntKw(kw, ever);
        }

        /* open or create the file using CFITSIO */
        if (ever == 0 || ename == 0 || ename[0] == '\0' || ename[0] == ' ') {
            c_vfn2osfn(fname, ospath);
            fits_create_file(&iodesc->ff, ospath, &status);
        } else {
            c_vfn2osfn(fname, ospath);
            fits_open_file(&iodesc->ff, ospath, READWRITE, &status);
        }
        if (status) {
            ioerr(BADOPEN, iodesc, status);
            free(iodesc->extname);
            free(iodesc->filename);
            free(iodesc);
            return NULL;
        }
        free(tmp);

        iodesc->dims[0] = d1;
        iodesc->dims[1] = d2;
        /* IMIO would always set bitpix to 16 when the dimensions are
           naught */
        if (d1 == 0 && d2 == 0) {
            iodesc->type = SHORT_IMG;
        } else {
            switch (typ) {
            case FITSBYTE:
                iodesc->type = BYTE_IMG;
                break;
            case FITSSHORT:
                iodesc->type = SHORT_IMG;
                break;
            case FITSLONG:
                iodesc->type = LONG_IMG;
                break;
            case FITSFLOAT:
                iodesc->type = FLOAT_IMG;
                break;
            case FITSDOUBLE:
                iodesc->type = DOUBLE_IMG;
                break;
            default:
                iodesc->type = SHORT_IMG;
                break;
            }
        }
        iodesc->hdr = hd;

        if (fits_create_img(iodesc->ff, iodesc->type, 2, iodesc->dims, &status)) {
            ioerr(BADOPEN, iodesc, status);
            return NULL;
        }

        if (fits_write_record(iodesc->ff, "ORIGIN  = 'HSTIO/CFITSIO March 2010'", &status)) {
            ioerr(BADWRITE, iodesc, status);
            return NULL;
        }

        /* CFITSIO TODO: write a real date */
        t = time(NULL);
        time_tmp = localtime(&t);
        strftime(date, 12, "%Y-%m-%d", time_tmp);
        snprintf(date_card, 80,
                 "DATE    = '%s' / date this file was written (yyyy-mm-dd)", date);
        if (fits_write_record(iodesc->ff, date_card, &status)) {
            ioerr(BADWRITE, iodesc, status);
            return NULL;
        }

        iodesc->hflag = 1; /* mark to write header */
        if (iodesc->dims[0] == 0) {
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        clear_err();

        return iodesc;
}

IODescPtr openUpdateImage(char *fname, char *ename, int ever, Hdr *hd) {
        IODesc *iodesc;
        int no_dims;
        char *tmp;
        char ospath[SZ_PATHNAME];
        int status = 0;

        /* CFITSIO: Error handling */
        c_pusherr(detect_iraferr);

        tmp = make_iodesc(&iodesc, fname, ename, ever);
        if (tmp == NULL) return NULL;
        iodesc->options = ReadWrite;

        /* CFITSIO: Resolve this keyword inheritance stuff */
        /* p = strstr(tmp,"[0]"); */
        /* if (p == NULL) { */
        /*     tmp[strlen(tmp) - 1] = '\0'; */
        /*     strcat(tmp,",NOINHERIT]"); */
        /* } */

        /* open the file using CFITSIO */
        c_vfn2osfn(tmp, ospath);

        if (fits_open_file(&iodesc->ff, ospath, READWRITE, &status)) {
            ioerr(BADOPEN, iodesc, status);
            free(tmp);
            free(iodesc->extname);
            free(iodesc->filename);
            free(iodesc);
            return NULL;
        }
        free(tmp);

        /* get the dimensions and type */
        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        fits_get_img_equivtype(iodesc->ff, &iodesc->type, &status);
        fits_get_img_size(iodesc->ff, 2, iodesc->dims, &status);
        if (status) {
            ioerr(BADDIMS, iodesc, status);
            return NULL;
        }
        if (no_dims == 2) {
            /* Nothing */
        } else if (no_dims == 1) {
            iodesc->dims[1] = 0;
        } else if (no_dims == 0) {
            iodesc->dims[0] = 0;
            iodesc->dims[1] = 0;
        } else {
            ioerr(BADDIMS, iodesc, 0);
            return NULL;
        }

        /* read the user area into the header array */
        getHeader(iodesc, hd);

        clear_err();
        return iodesc;
}

void closeImage(IODescPtr iodesc_) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int status = 0;

        if (iodesc->options != ReadOnly && iodesc->dims[0] != 0)
            putHeader(iodesc);

        if (fits_close_file(iodesc->ff, &status)) {
            /* TODO: Raise error */
        }

        /* This is a handy check to use pyfits to validate the file upon every close */
        /* c_vfn2osfn(iodesc->filename, ospath); */
        /* sprintf(system_string, "python -c \"import pyfits; pyfits.open('%s')\"", ospath); */
        /* if (system(system_string)) { */
        /*   printf("LOG: pyfits corruption!!!\n"); */
        /*   exit(1); */
        /* } */

        /* if (there is an IRAF error) */
        /*      ioerr(IRAF_CLOSE,iodesc); */
        free(iodesc->extname);
        free(iodesc->filename);
        free(iodesc);

        /* CFITSIO: Error handling */
        c_poperr();
}


/* According to the imio documentation, the following reserved
   keywords are recognized:

   SIMPLE BITPIX DATATYPE NAXIS* GROUPS GCOUNT PCOUNT PSIZE
   PTYPE* PDTYPE* PSIZE* XTENSION
*/

static char* reservedKwds[] = {
    "BITPIX  ",
    "BSCALE  ",
    "BZERO   ",
    "DATAMAX ",
    "DATAMIN ",
    "DATATYPE",
    "DATE    ",
    "GCOUNT  ",
    "GROUPS  ",
    "NAXIS   ",
    "NAXIS*  ",
    "ORIGIN  ",
    "PCOUNT  ",
    "PDTYPE* ",
    "PSIZE   ",
    "PSIZE*  ",
    "PTYPE*  ",
    "SIMPLE  ",
    "XTENSION",
    NULL
};

int isReservedKwd(const char* card) {
        /* CFITSIO: Should this be made case-insensitive? */
        /* TODO: Maybe use a binary search? */

        /* Returns 1 if the card matches one of the reserved keywords */
        int i;
        int match;
        char** kwd = reservedKwds;

        while (*kwd != NULL) {
            /* Short-circuit if we're certain not to find the kwd
               later in the list */
            if ((*kwd)[0] > card[0]) {
                return 0;
            }

            match = 1;
            for (i = 0; i < 8; ++i) {
                /* '*' indicates a digit */
                if ((*kwd)[i] == '*') {
                    if (card[i] < '0' || card[i] > '9') {
                        match = 0;
                        break;
                    }
                } else if ((*kwd)[i] != card[i]) {
                    match = 0;
                    break;
                }
            }
            if (match) {
                return 1;
            }

            ++kwd;
        }

        return 0;
}

int getHeader(IODescPtr iodesc_, Hdr *hd) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int ncards, i, j;
        char source[HDRSize];
        char *target;
        int status = 0;

        if (iodesc->options == WriteOnly) {
            ioerr(NOGET, iodesc, 0);
            return -1;
        }

        /* get the number of cards in the header */
        if (fits_get_hdrspace(iodesc->ff, &ncards, NULL, &status)) {
            ioerr(BADHSIZE, iodesc, status);
            return -1;
        }

        /* allocate space for the header cards */
        if (allocHdr(hd, ncards) == -1) return -1;

        /* translate the data */
        hd->nlines = 0;
        for (i = 0; i < ncards; ++i) {
            target = hd->array[i];
            if (fits_read_record(iodesc->ff, i+1, source, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
            if (!isReservedKwd(source)) {
                for (j = 0; j < (HDRSize -1); ++j) {
                    *target++ = source[j];
                }
            }
            *target++ = '\0';
            hd->nlines++;
        }
        iodesc->hdr = hd;

        clear_err();
        return 0;
}

int putHeader(IODescPtr iodesc_) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j, tmp;
        int numkeys;
        char *source;
        char card[80];
        int status = 0;

        if (iodesc->options == ReadOnly) {
            ioerr(NOPUT, iodesc, status);
            return -1;
        }

        if (iodesc->hflag) {
            /* CFITSIO: We probably need to move this in front of all
               calls to fits_create_img */

            /* If the image is actually 1-dimensional, modify the naxis2
             * value so the output header is written with only NAXIS and
             * NAXIS1 keywords, where NAXIS=1, and NAXIS1=number.
             */
            if (iodesc->dims[0] != 0 && iodesc->dims[1] == 1)
                iodesc->dims[1] = 0;

            /* set the pixel type */
            fits_update_key(iodesc->ff, TINT, "BITPIX", &(iodesc->type), "", &status);
            if (iodesc->dims[0] == 0 && iodesc->dims[1] == 0) {
                tmp = 0;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &tmp, "", &status);
                fits_delete_key(iodesc->ff, "NAXIS1", &status);
                fits_delete_key(iodesc->ff, "NAXIS2", &status);
            } else if (iodesc->dims[0] != 0 && iodesc->dims[1] == 0) {
                /* set the number of dimensions */
                tmp = 1;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &tmp, "", &status);
                /* set dim1 */
                fits_update_key(iodesc->ff, TINT, "NAXIS1", &iodesc->dims[0], "", &status);
                fits_delete_key(iodesc->ff, "NAXIS2", &status);
            } else {
                /* set the number of dimensions */
                tmp = 2;
                fits_update_key(iodesc->ff, TINT, "NAXIS", &tmp, "", &status);
                /* set dim1 and dim2 */
                fits_update_key(iodesc->ff, TINT, "NAXIS1", &iodesc->dims[0], "", &status);
                fits_update_key(iodesc->ff, TINT, "NAXIS2", &iodesc->dims[1], "", &status);
            }

            if (status) {
                ioerr(BADWRITE, iodesc, status);
            }
        }

        /* Verify the size of the user area */
        /* The original code just memcopies the cards into the "user
           area" of the header.  CFITSIO doesn't have the concept of a
           "user area", so we need to carefully only copy the cards
           that are not "reserved".
        */
        if (fits_get_hdrspace(iodesc->ff, &numkeys, NULL, &status)) {
            ioerr(BADWRITE, iodesc, status);
            return -1;
        }

        for (i = 0, j = numkeys; i < numkeys; ++i, --j) {
            if (fits_read_record(iodesc->ff, j, card, &status)) {
                ioerr(BADWRITE, iodesc, status);
                return -1;
            }
            if (!isReservedKwd(card)) {
                if (fits_delete_record(iodesc->ff, j, &status)) {
                    ioerr(BADWRITE, iodesc, status);
                    return -1;
                }
            } else {
                ++j;
            }
        }

        /* translate the data */

        /* Skip blank cards at the beginning */
        for (i = 0; i < iodesc->hdr->nlines; ++i) {
            source = iodesc->hdr->array[i];
            for (j = 0; j < 80; ++j) {
                if (source[j] != ' ' && 
                    source[j] != '\n' && 
                    source[j] != 0) {
                    goto found_non_space;
                }
            }
        }

 found_non_space:
        for (/* i from above */; i < iodesc->hdr->nlines; ++i) {
            source = iodesc->hdr->array[i];
            /* CFITSIO TODO: Strictly speaking, we may want to write
               even reserved keywords here... not sure */
            if (!isReservedKwd(source)) {
                if (fits_write_record(iodesc->ff, source, &status)) {
                    ioerr(BADWRITE, iodesc, status);
                    return -1;
                }
            }
        }

        /* If we don't explicitly set BSCALE and BZERO to 1.0 and 0.0
           here, their values could be inadvertently brought over from
           the source image.  This was the source of a very
           hard-to-find bug. */
        if (iodesc->type == TFLOAT || iodesc->type == TDOUBLE) {
          fits_set_bscale(iodesc->ff, 1.0, 0.0, &status);
        }

        clear_err();
        return 0;
}

int getFloatData(IODescPtr iodesc_, FloatTwoDArray *da) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int no_dims, i, j;
        long fpixel[2];
        int anynul;
        int type;
        FitsKw kw;
        float val;
        int status = 0;

        if (iodesc->options == WriteOnly) { ioerr(NOGET,iodesc, 0); return -1; }

        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        fits_get_img_size(iodesc->ff, 2, iodesc->dims, &status);
        if (status) {
            ioerr(BADDIMS, iodesc, status);
            return -1;
        }

        /*
           If the number  of dimensions of the image is zero, need to
           determine how many dimensions the image is supposed to have
           according to the NPIX[1/2] keyword(s).
        */
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc, 0); return -1; }
            val = getFloatKw(kw);

            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc, 0); return -1; }
            iodesc->dims[0] = getIntKw(kw);

            /* If NPIX2 is not found, the image should be 1D; dim2 = 1 and *
             * not 0 for purposes of memory allocation.                    */
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw == 0)  {
                iodesc->dims[1] = 1;
            } else {
                iodesc->dims[1] = getIntKw(kw);
            }

            if (allocFloatData(da, iodesc->dims[0], iodesc->dims[1])) return -1;
            for (j = 0; j < iodesc->dims[1]; ++j) {
                for (i = 0; i < iodesc->dims[0]; ++i) {
                    PPix(da, i, j) = val;
                }
            }
        } else if (no_dims == 1) {
            iodesc->dims[1] = 1;
            fits_get_img_equivtype(iodesc->ff, &type, &status);
            /* CFITSIO TODO: Should we verify the type is correct
               here?  Original code gets type, but then does nothing
               with it. */
            if (allocFloatData(da, iodesc->dims[0], iodesc->dims[1])) return -1;
/*
            if (c_imgnlr(iodesc->fdesc,&x,linevector) == IRAF_EOF) {
                ioerr(BADREAD,iodesc); return -1; }
*/
            fpixel[0] = 1;
            fpixel[1] = 1;
            if (fits_read_pix(iodesc->ff, TFLOAT, fpixel, iodesc->dims[0], 0,
                              (float *)&(PPix(da, 0, 0)), &anynul, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
        } else if (no_dims == 2) {
            fits_get_img_equivtype(iodesc->ff, &type, &status);
            /* CFITSIO TODO: Should we verify the type is correct
               here?  Original code gets type, but then does nothing
               with it. */
            if (allocFloatData(da, iodesc->dims[0], iodesc->dims[1])) return -1;

            fpixel[0] = 1;
            for (i = 0; i < iodesc->dims[1]; ++i) {
                fpixel[1] = i + 1;
                if (fits_read_pix(iodesc->ff, TFLOAT, fpixel, iodesc->dims[0], 0,
                                  (float *)&(PPix(da, 0, i)), &anynul, &status)) {
                    ioerr(BADREAD,iodesc, status);
                    return -1;
                }
            }
        } else {
            ioerr(BADDIMS,iodesc,0);
            return -1;
        }

        clear_err();
        return 0;
}

int putFloatData(IODescPtr iodesc_, FloatTwoDArray *da) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j;
        float tmp;
        long fpixel[2];
        FitsKw kw;
        int is_eq;
        int naxis;
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }

        /* check for a constant array, if not SCI data */
        if (strcmp(iodesc->extname,"SCI") != 0
            && da->tot_nx != 0 && da->tot_ny != 0) {
            tmp = PPix(da,0,0);
            for (i = 0, is_eq = 1; (i < da->tot_nx) && is_eq; ++i) {
                for (j = 0; (j < da->tot_ny); ++j) {
                    if (PPix(da,i,j) != tmp) {
                        is_eq = 0;
                        break;
                    }
                }
            }
            if (is_eq) {
                /* This is a constant array. */
                /* add NPIX1, NPIX2 (if necessary), and PIXVALUE keywords */
                kw = findKw(iodesc->hdr,"PIXVALUE");
                if (kw == 0) /* add it */
                    addFloatKw(iodesc->hdr,"PIXVALUE",tmp,
                        "values of pixels in constant array");
                else
                    putFloatKw(kw,tmp);

                kw = findKw(iodesc->hdr,"NPIX1");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"NPIX1",iodesc->dims[0],
                        "length of constant array axis 1");
                else
                    putIntKw(kw,iodesc->dims[0]);

                /* NPIX2 should only be added if the y-dimension is > 1. */
                if (da->tot_ny > 1) {
                    kw = findKw(iodesc->hdr,"NPIX2");
                    if (kw == 0) /* add it */
                        addIntKw(iodesc->hdr,"NPIX2",iodesc->dims[1],
                            "length of constant array axis 2");
                    else
                        putIntKw(kw,iodesc->dims[1]);
                }

                naxis = 0;
                fits_write_key(iodesc->ff, TINT, "NAXIS", &naxis, NULL, &status);
                iodesc->dims[0] = 0;
                iodesc->dims[1] = 0;
                /* update the header, etc. */
                if (iodesc->hflag) {
                    iodesc->type = FLOAT_IMG;
                    putHeader(iodesc);
                    iodesc->hflag = 0;
                }
                fits_flush_file(iodesc->ff, &status);

                return 0;
            }
        }

        /* If not a constant array, make sure NPIX1, NPIX2, and PIXVALUE *
         * are NOT present in the header to be written out.              */
        kw = findKw(iodesc->hdr,"NPIX1");
        if (kw != 0) /* remove it */
            delKw(kw);

        if (da->tot_ny > 1) {
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw != 0) /* remove it */
                delKw(kw);
        }

        kw = findKw(iodesc->hdr,"PIXVALUE");
        if (kw != 0) /* remove it */
            delKw(kw);

        /* update the header area */
        if (iodesc->hflag) {
            iodesc->type = FLOAT_IMG;
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        fpixel[0] = 1;
        for (i = 0; i < da->ny; ++i) {
            fpixel[1] = i + 1;
            if (fits_write_pix(iodesc->ff, TFLOAT, fpixel, da->nx,
                               (float *)&(PPix(da, 0, i)), &status)) {
                ioerr(BADWRITE, iodesc, status);
                return -1;
            }
        }

        fits_flush_file(iodesc->ff, &status);

        clear_err();
        return 0;
}

/*                                                                     **
** Write output a subsection of an image in memory to a file where the **
** subsection is the full size of the output data.                     **
**                                                                     */
int putFloatSect(IODescPtr iodesc_, FloatTwoDArray *da, int xbeg,
                 int ybeg, int xsize, int ysize) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j, xend, yend;
        float tmp;
        FitsKw kw;
        long fpixel[2];
        int is_eq;
        int naxis;
        int status = 0;

        /* CFITSIO: Verify that the section is within range? */

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc, 0); return -1; }

        xend = xbeg + xsize;
        yend = ybeg + ysize;
        /* check for a constant array, if not SCI data */
        if (strcmp(iodesc->extname,"SCI") != 0
            && da->tot_nx != 0 && da->tot_ny != 0) {
            tmp = PPix(da, 0, 0);
            for (i = xbeg, is_eq = 1; (i < xend) && is_eq; ++i) {
                for (j = ybeg; (j < yend); ++j) {
                    if (PPix(da,i,j) != tmp) {
                        is_eq = 0;
                        break;
                    }
                }
            }
            if (is_eq) {
                /* This is a constant array. */
                /* add NPIX1, NPIX2 (if necessary), and PIXVALUE keywords */
                kw = findKw(iodesc->hdr,"PIXVALUE");
                if (kw == 0) /* add it */
                    addFloatKw(iodesc->hdr,"PIXVALUE",tmp,
                        "values of pixels in constant array");
                else
                    putFloatKw(kw,tmp);

                kw = findKw(iodesc->hdr,"NPIX1");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"NPIX1",iodesc->dims[0],
                        "length of constant array axis 1");
                else
                    putIntKw(kw,iodesc->dims[0]);

                /* NPIX2 should only be added if the y-dimension is > 1. */
                if (da->tot_ny > 1) {
                    kw = findKw(iodesc->hdr,"NPIX2");
                    if (kw == 0) /* add it */
                        addIntKw(iodesc->hdr,"NPIX2",iodesc->dims[1],
                            "length of constant array axis 2");
                    else
                        putIntKw(kw,iodesc->dims[1]);
                }

                naxis = 0;
                fits_write_key(iodesc->ff, TINT, "NAXIS", &naxis, NULL, &status);
                iodesc->dims[0] = 0;
                iodesc->dims[1] = 0;
                /* update the header, etc. */
                if (iodesc->hflag) {
                    iodesc->type = FLOAT_IMG;
                    putHeader(iodesc);
                    iodesc->hflag = 0;
                }

                fits_flush_file(iodesc->ff, &status);

                clear_err();
                return 0;
            }
        }

        /* If not a constant array, make sure NPIX1, NPIX2, and PIXVALUE *
         * are NOT present in the header to be written out.              */
        kw = findKw(iodesc->hdr,"PIXVALUE");
        if (kw != 0) /* remove it */
            delKw(kw);

        kw = findKw(iodesc->hdr,"NPIX1");
        if (kw != 0) /* remove it */
            delKw(kw);

        if (da->tot_ny > 1) {
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw != 0) /* remove it */
                delKw(kw);
        }

        /* update the header area */
        if (iodesc->hflag) {
            iodesc->type = FLOAT_IMG;
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        fpixel[0] = 1;
        for (i = ybeg; i < yend; ++i) {
            fpixel[1] = i - ybeg + 1;
            if (fits_write_pix(iodesc->ff, TFLOAT, fpixel, xsize,
                               (float*)&(PPix(da, xbeg, i)), &status)) {
                ioerr(BADWRITE, iodesc, status);
                return -1;
            }
        }

        fflush(stdout);

        clear_err();
        return 0;
}

int getShortData(IODescPtr iodesc_, ShortTwoDArray *da) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int no_dims, i, j;
        FitsKw kw;
        short val;
        long fpixel[2];
        int anynul = 0;
        int status = 0;

        if (iodesc->options == WriteOnly) { ioerr(NOGET,iodesc, 0); return -1; }

        fits_get_img_dim(iodesc->ff, &no_dims, &status);
        fits_get_img_size(iodesc->ff, 2, iodesc->dims, &status);
        if (status) {
            ioerr(BADDIMS, iodesc, status);
            return -1;
        }

        /*
           If the number  of dimensions of the image is zero, need to
           determine how many dimensions the image is supposed to have
           according to the NPIX[1/2] keyword(s).
        */
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            val = getIntKw(kw);

            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            iodesc->dims[0] = getIntKw(kw);

            /* If NPIX2 is not found, the image should be 1D; dim2 = 1 and *
             * not 0 for purposes of memory allocation.                    */
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw == 0)  {
                iodesc->dims[1] = 1;
            } else {
                iodesc->dims[1] = getIntKw(kw);
            }

            if (allocShortData(da, iodesc->dims[0], iodesc->dims[1])) return -1;
            for (j = 0; j < iodesc->dims[1]; ++j)
                for (i = 0; i < iodesc->dims[0]; ++i)
                    PPix(da, i, j) = val;
        } else if (no_dims == 1) {
            iodesc->dims[1] = 1;
            /* CFITSIO TODO: Should we verify the type is correct
               here?  Original code gets type, but then does nothing
               with it. */
            /* type = (IRAFType)c_imgtypepix(iodesc->fdesc); */
            if (allocShortData(da, iodesc->dims[0], iodesc->dims[1])) return -1;
/*
            if (c_imgnls(iodesc->fdesc,&x,linevector) == IRAF_EOF) {
                ioerr(BADREAD,iodesc); return -1; }
*/
            fpixel[0] = 1;
            fpixel[1] = 1;
            if (fits_read_pix(iodesc->ff, TSHORT, fpixel, iodesc->dims[0], NULL,
                              (short *)&(PPix(da, 0, 0)), &anynul, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
        } else if (no_dims == 2) {
            /* CFITSIO TODO: Should we verify the type is correct
               here?  Original code gets type, but then does nothing
               with it. */
            /* type = (IRAFType)c_imgtypepix(iodesc->fdesc); */
            if (allocShortData(da, iodesc->dims[0], iodesc->dims[1])) return -1;
            fpixel[0] = 1;
            for (i = 0; i < iodesc->dims[1]; ++i) {
                fpixel[1] = i + 1;
                if (fits_read_pix(iodesc->ff, TSHORT, fpixel, iodesc->dims[0], NULL,
                                  (short *)&(PPix(da, 0, i)), &anynul, &status)) {
                    ioerr(BADREAD, iodesc, status);
                    return -1;
                }
            }
        } else {
            ioerr(BADDIMS, iodesc, 0);
            return -1;
        }

        clear_err();
        return 0;
}

int putShortData(IODescPtr iodesc_, ShortTwoDArray *da) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j;
        short tmp;
        long fpixel[2];
        FitsKw kw;
        int is_eq;
        int naxis;
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }

        /* check for a constant array, if not SCI data */
        if (strcmp(iodesc->extname,"SCI") != 0
            && da->tot_nx != 0 && da->tot_ny != 0) {
            tmp = PPix(da,0,0);
            for (i = 0, is_eq = 1; (i < da->tot_nx) && is_eq; ++i) {
                for (j = 0; (j < da->tot_ny); ++j) {
                    if (PPix(da,i,j) != tmp) {
                        is_eq = 0;
                        break;
                    }
                }
            }
            if (is_eq) {
                /* This is a constant array. */
                /* add NPIX1, NPIX2 (if necessary), and PIXVALUE keywords */
                kw = findKw(iodesc->hdr,"PIXVALUE");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"PIXVALUE",(int)tmp,
                        "values of pixels in constant array");
                else
                    putIntKw(kw,(int)tmp);

                kw = findKw(iodesc->hdr,"NPIX1");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"NPIX1",iodesc->dims[0],
                        "length of constant array axis 1");
                else
                    putIntKw(kw,iodesc->dims[0]);

                /* NPIX2 should only be added if the y-dimension is > 1. */
                if (da->tot_ny > 1) {
                    kw = findKw(iodesc->hdr,"NPIX2");
                    if (kw == 0) /* add it */
                        addIntKw(iodesc->hdr,"NPIX2",iodesc->dims[1],
                            "length of constant array axis 2");
                    else
                        putIntKw(kw,iodesc->dims[1]);
                }

                naxis = 0;
                fits_write_key(iodesc->ff, TINT, "NAXIS", &naxis, NULL, &status);
                iodesc->dims[0] = 0;
                iodesc->dims[1] = 0;
                /* update the header, etc. */
                if (iodesc->hflag) {
                    iodesc->type = SHORT_IMG;
                    putHeader(iodesc);
                    iodesc->hflag = 0;
                }

                fits_flush_file(iodesc->ff, &status);

                clear_err();
                return 0;
            }
        }

        /* If not a constant array, make sure NPIX1, NPIX2, and PIXVALUE *
         * are NOT present in the header to be written out.              */
        kw = findKw(iodesc->hdr,"NPIX1");
        if (kw != 0) /* remove it */
            delKw(kw);

        if (da->tot_ny > 1) {
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw != 0) /* remove it */
                delKw(kw);
        }

        kw = findKw(iodesc->hdr,"PIXVALUE");
        if (kw != 0) /* remove it */
            delKw(kw);

        /* update the header area */
        if (iodesc->hflag) {
            iodesc->type = SHORT_IMG;
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        fpixel[0] = 1;
        for (i = 0; i < da->ny; ++i) {
            fpixel[1] = i + 1;
            if (fits_write_pix(iodesc->ff, TSHORT, fpixel, da->nx,
                               (short *)&(PPix(da, 0, i)), &status)) {
                ioerr(BADWRITE, iodesc, status); return -1;
            }
        }

        fits_flush_file(iodesc->ff, &status);

        clear_err();
        return 0;
}

/*                                                                     **
** Write output a subsection of an image in memory to a file where the **
** subsection is the full size of the output data.                     **
**                                                                     */
int putShortSect(IODescPtr iodesc_, ShortTwoDArray *da, int xbeg, int ybeg,
                     int xsize, int ysize) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int i, j, xend, yend;
        short tmp;
        FitsKw kw;
        int naxis;
        int is_eq;
        long fpixel[2];
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }

        xend = xbeg + xsize;
        yend = ybeg + ysize;
        /* check for a constant array, if not SCI data */
        if (strcmp(iodesc->extname,"SCI") != 0
            && da->tot_nx != 0 && da->tot_ny != 0) {
            tmp = PPix(da,0,0);
            for (i = xbeg, is_eq = 1; (i < xend) && is_eq; ++i) {
                for (j = ybeg; (j < yend); ++j) {
                    if (PPix(da,i,j) != tmp) {
                        is_eq = 0;
                        break;
                    }
                }
            }
            if (is_eq) {
                /* This is a constant array. */
                /* add NPIX1, NPIX2 (if necessary), and PIXVALUE keywords */
                kw = findKw(iodesc->hdr,"PIXVALUE");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"PIXVALUE",(int)tmp,
                        "values of pixels in constant array");
                else
                    putIntKw(kw,(int)tmp);

                kw = findKw(iodesc->hdr,"NPIX1");
                if (kw == 0) /* add it */
                    addIntKw(iodesc->hdr,"NPIX1",iodesc->dims[0],
                        "length of constant array axis 1");
                else
                    putIntKw(kw,iodesc->dims[0]);

                /* NPIX2 should only be added if the y-dimension is > 1. */
                if (da->tot_ny > 1) {
                    kw = findKw(iodesc->hdr,"NPIX2");
                    if (kw == 0) /* add it */
                        addIntKw(iodesc->hdr,"NPIX2",iodesc->dims[1],
                            "length of constant array axis 2");
                    else
                        putIntKw(kw,iodesc->dims[1]);
                }

                naxis = 0;
                fits_write_key(iodesc->ff, TINT, "NAXIS", &naxis, NULL, &status);
                iodesc->dims[0] = 0;
                iodesc->dims[1] = 0;
                /* update the header, etc. */
                if (iodesc->hflag) {
                    iodesc->type = SHORT_IMG;
                    putHeader(iodesc);
                    iodesc->hflag = 0;
                }

                fits_flush_file(iodesc->ff, &status);

                clear_err();
                return 0;
            }
        }

        /* If not a constant array, make sure NPIX1, NPIX2, and PIXVALUE *
         * are NOT present in the header to be written out.              */
        kw = findKw(iodesc->hdr,"PIXVALUE");
        if (kw != 0) /* remove it */
            delKw(kw);

        kw = findKw(iodesc->hdr,"NPIX1");
        if (kw != 0) /* remove it */
            delKw(kw);

        if (da->tot_ny > 1) {
            kw = findKw(iodesc->hdr,"NPIX2");
            if (kw != 0) /* remove it */
                delKw(kw);
        }

        /* update the header area */
        if (iodesc->hflag) {
            iodesc->type = SHORT_IMG;
            putHeader(iodesc);
            iodesc->hflag = 0;
        }

        fpixel[0] = 1;
        for (i = ybeg; i < yend; ++i) {
            fpixel[1] = i - ybeg + 1;
            if (fits_write_pix(iodesc->ff, TSHORT, fpixel, xsize,
                               (short *)&(PPix(da, xbeg, i)), &status)) {
                ioerr(BADWRITE,iodesc, status);
                return -1;
            }
        }

        fits_flush_file(iodesc->ff, &status);

        clear_err();
        return 0;
}

int getFloatLine(IODescPtr iodesc_, int line, float *ptr) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int no_dims, i, dim1;
        long dims[2];
        FitsKw kw;
        float val;
        long fpixel[2];
        int anynul;
        int status = 0;

        if (iodesc->options == WriteOnly) { ioerr(NOGET,iodesc,0); return -1; }

        if (fits_get_img_dim(iodesc->ff, &no_dims, &status)) {
            ioerr(BADDIMS, iodesc, status);
            return -1;
        }
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            dim1 = getIntKw(kw);
            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            val = getFloatKw(kw);
            for (i = 0; i < dim1; ++i) {
                ptr[i] = val;
            }
        } else {
            if (fits_get_img_size(iodesc->ff, 2, dims, &status)) {
                ioerr(BADDIMS, iodesc, status);
                return -1;
            }
            fpixel[0] = 1;
            fpixel[1] = line + 1;
            if (fits_read_pix(iodesc->ff, TFLOAT, fpixel, dims[0], NULL,
                              ptr, &anynul, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
        }
        clear_err();
        return 0;
}

int putFloatLine(IODescPtr iodesc_, int line, float *ptr) {
        IODesc *iodesc = (IODesc *)iodesc_;
        long fpixel[2];
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }
        if (iodesc->hflag) { iodesc->hflag = 0; putHeader(iodesc); }

        fpixel[0] = 1;
        fpixel[1] = line + 1;
        if (fits_write_pix(iodesc->ff, TFLOAT, fpixel, iodesc->dims[0],
                           ptr, &status)) {
            ioerr(BADWRITE, iodesc, status);
            return -1;
        }

        clear_err();
        return 0;
}

int getShortLine(IODescPtr iodesc_, int line, short *ptr) {
        IODesc *iodesc = (IODesc *)iodesc_;
        int no_dims, dim1, i;
        long dims[2];
        FitsKw kw;
        short val;
        long fpixel[2];
        int anynul;
        int status = 0;

        if (iodesc->options == WriteOnly) { ioerr(NOGET,iodesc,0); return -1; }

        if (fits_get_img_dim(iodesc->ff, &no_dims, &status)) {
            ioerr(BADDIMS, iodesc, status);
            return -1;
        }
        if (no_dims == 0) {
            kw = findKw(iodesc->hdr,"NPIX1");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            dim1 = getIntKw(kw);
            kw = findKw(iodesc->hdr,"PIXVALUE");
            if (kw == 0) { ioerr(BADSCIDIMS,iodesc,0); return -1; }
            val = getIntKw(kw);
            for (i = 0; i < dim1; ++i)
                ptr[i] = val;
        } else {
            if (fits_get_img_size(iodesc->ff, 2, dims, &status)) {
                ioerr(BADDIMS, iodesc, status);
                return -1;
            }
            fpixel[0] = 1;
            fpixel[1] = line + 1;
            if (fits_read_pix(iodesc->ff, TSHORT, fpixel, dims[0], NULL,
                              ptr, &anynul, &status)) {
                ioerr(BADREAD, iodesc, status);
                return -1;
            }
        }
        clear_err();
        return 0;
}

int putShortLine(IODescPtr iodesc_, int line, short *ptr) {
        IODesc *iodesc = (IODesc *)iodesc_;
        long fpixel[2];
        int status = 0;

        if (iodesc->options == ReadOnly) { ioerr(NOPUT,iodesc,0); return -1; }
        if (iodesc->hflag) { iodesc->hflag = 0; putHeader(iodesc); }

        fpixel[0] = 1;
        fpixel[1] = line + 1;
        if (fits_write_pix(iodesc->ff, TSHORT, fpixel, iodesc->dims[0],
                           ptr, &status)) {
            ioerr(BADWRITE, iodesc, status);
            return -1;
        }

        clear_err();
        return 0;
}
