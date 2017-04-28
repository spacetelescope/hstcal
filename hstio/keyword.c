# include "hstio.h"
# include "numeric.h"
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <ctype.h>

# include <float.h>
# include <limits.h>

/*
 * M.D. De La Pena 28 January 1998 - Addressed problems with improperly
 * overwriting an existing string keyword value and possibly also corrupting
 * the associated comment field.  Source of problem is the public routine
 * putStringKw.  Modified the putStringKw function to be only an interim
 * function to obtain and copy the existing comment field; this routine will
 * also re-write the original comment after the call to a new routine,
 * putString, which puts the new keyword string value.  putString is identical
 * to the original putStringKw function; this was done so that the public
 * interface would not have to change.  The routine insertStringKw was also
 * modified to call putString.
 *
 * M.D. De La Pena 16 April 1998: Modified putKwName to ensure the output
 * keyword is upper case.  Modified putString to check for embedded single
 * quotes in input strings; force the string written to the header to
 * conform to the FITS standard.
 *
 * M.D. De La Pena 28 September 1998 - Created a new routine, searchBlankLine,
 * which searches backwards from the last line in the header to determine if
 * the header contains any trailing blank lines.  Modified addcommentary,
 * insertname, and insertFitsCard to call a new routine, searchBlankLine,
 * which determines if there are blank lines padding the end of the header.
 *
 * M.D. De La Pena 07 October 1998 - Modified getIntKw to accommodate reading
 * FITSDOUBLE values and getFloatKw to read FITSLONG values.  Casts are
 * applied as appropriate.
 *
 * M.D. De La Pena 10 November 1998 - Incorporated high-level keyword access
 * routines written by H. Bushouse with some modification.
 *
 * M.D. De La Pena 10 December 1998 - Added range checking for keyword values
 * being converted into INTs or FLOATs.  A new routine, checkRange, is used
 * by getIntKw and getFloatKw.
 * VMS issues: Modified putStringKw routine to return 0 and not NULL if
 * putString == -1.  searchBlankLine is no longer a static routine.
 *
 * M.D. De La Pena 19 August 1999 - Moved searchBlankLine to top of file and
 * made static.
 *
 * Phil Hodge 31 August 2006 - Modified putKwComm to avoid writing beyond
 * 80 characters.  This occurred while writing " / " after a long
 * character-string value.
 *
 */


typedef struct {
        char *text;             /* the text string */
        Hdr *hdr;               /* the header array */
        char *value_end;        /* position of end of value */
        Bool isparsed;          /* flag telling whether parsed or not */
        int index;              /* position of string within header arrray */
        char name[12];          /* the keyword name */
        FitsDataType type;      /* the FITS data type of the value */
        Bool bresult;           /* the Bool value, if a FITS LOGICAL */
        char cresult[81];       /* the string values, if FITS STRING */
        NumericResult nresult;  /* the numeric value, if numeric */
} FitsKwInfo;
/*printf("FitsKwInfo: isparsed = %d index = %d name = %s type = %d \
nresult = %d\n|%s|\n|%s|\n",kw->isparsed,kw->index,kw->name,kw->type,
kw->nresult.data.l,kw->text,kw->value_end);*/

# define HdrUnit 36

# define BLANK   0x20

# define BLANK_CARD "\
                                                                                \
"

# if defined(__cplusplus)
extern "C" {
# endif
extern void error(HSTIOError, char *);
# if defined(__cplusplus)
}
# endif

enum TargetDataType_ { T_INTEGER = 0, T_FLOAT};
typedef enum TargetDataType_ TargetDataType;

int checkRange (TargetDataType target, NumericResult result);
int   putString(FitsKw kw_, char *txt);

static char *keymsg(char *k) {
        static char tmp[81] = "Searching for keyword ";
        if (strlen(k) > (size_t)58)
                { strncpy(&tmp[22],k,58); tmp[80] = '\0'; }
        else
                strcpy(&tmp[22],k);
        return tmp;
}

static FitsKwInfo findkw = { NULL, NULL, NULL, False, -1, {'\0'},
        FITSNOVALUE, False, {'\0'} };

static int find(Hdr *h, char *nm, FitsKwInfo *kw) {
        int i, n;

        n = strlen(nm);
        for (i = kw->index; i < h->nlines; ++i)
            if (strncmp(h->array[i],nm,n) == 0) {
                if (n == 8) break;
                if (h->array[i][n] == ' ') break;
            }
        if (i == h->nlines) return -1;

        strcpy(kw->name,nm);
        kw->hdr = h;
        kw->index = i;
        kw->text = h->array[i];
        kw->isparsed = False;
        return 0;
}

static int getvalue(FitsKwInfo *kw) {
        char *p, *s, *t;

        kw->type = FITSNOVALUE;
        kw->bresult = False;
        kw->cresult[0] = '\0';
        kw->nresult.data.l = 0;

        if (kw->text[8] != '=') {
            kw->isparsed = True;
            kw->value_end = &(kw->text[7]);
            return 0;
        }

        if (kw->index == -1 || kw->hdr == NULL) {
                error(BADGET,keymsg(kw->name)); return -1; }
        if (kw->text[8] != '=') {
                error(BADFITSEQ,keymsg(kw->name)); return -1; }
        p = &kw->text[10];

        /* look for a Bool */
        while (*p == ' ') ++p;
        if (*p == 'T' || *p == 'F') {
            kw->type = FITSLOGICAL;
            kw->bresult = (*p == 'T' ? True : False);
            kw->isparsed = True;
            kw->value_end = p;
            return 0;
        }
        /* look for a string */
        if (*p == '\'') {
            ++p;
            s = p;
            t = kw->cresult;
            for (;;) {
                if (*s == '\0') {
                        error(BADFITSQUOTE,keymsg(kw->name)); return -1; }
                if (*s == '\'') {
                    if (*(s + 1) == '\'') {
                        ++s;
                        *t = '\'';
                    } else {
                        --t;
                        while (*t == ' ' && t > &kw->cresult[0]) --t;
                        ++t;
                        *t = '\0';
                        kw->type = FITSCHAR;
                        kw->isparsed = True;
                        kw->value_end = s;
                        return 0;
                    }
                }
                *t++ = *s++;
            }
        }
        /* then, it must be numeric */
        get_numeric(p, 70, &kw->nresult);
        switch (kw->nresult.type) {
            case 0: kw->type = FITSNOVALUE;
                        error(BADFITSNUMERIC,keymsg(kw->name)); return -1;
            case 1: kw->type = FITSLONG; break;
            case 2: kw->type = FITSDOUBLE; break;
        }
        kw->isparsed = True;
        kw->value_end = &p[kw->nresult.endpos];
        return 0;
}

/* Routine to search from the last known line backwards to determine if *
 * there are trailing blank lines in the file.                          */
static int searchBlankLine (Hdr *h) {
    int  i, j;
    int  nblankline = 0;
    Bool blankline;
    char *c;

    /* Search backwards from the last header line (not including the END *
     * card) to determine if there exist trailing blank lines which can  *
     * be overwritten.                                                   */
    if (h->nlines != 0) {
        for (i = h->nlines - 1; i >= 0; i--) {
             c = h->array[i];
             blankline = True;
             for (j = 0; j < 80; j++) {
                  if (*c != BLANK) {
                      blankline = False;
                      break;
                  }
                  c++;
             }
             if (blankline == False)
                 break;
             else
                 nblankline++;
        }
    }
    return (nblankline);
}

static int insertname(FitsKwInfo *kw, char *nm) {
        int i;
        int nblankline  = 0;
        int orig_nlines = 0;

        if (kw->hdr->nalloc == 0) {
            if (allocHdr(kw->hdr,HdrUnit, True) != 0)
            return -1;
        }
        if (kw->hdr->nlines == kw->hdr->nalloc) {
            if (reallocHdr(kw->hdr,(kw->hdr->nalloc + HdrUnit)) != 0)
            return -1;
        }

        /* Look for the last true line of text (i.e., search backwards over *
         * any blank lines) not including the END card.                     */
        nblankline = searchBlankLine (kw->hdr);

        /* Save the original number of header lines and reset nlines to the *
         * first trailing blank line.                                       */
        if (nblankline != 0) {
            orig_nlines     = kw->hdr->nlines;
            kw->hdr->nlines = kw->hdr->nlines - nblankline;
            if (kw->index == orig_nlines - 1)
                kw->index = kw->hdr->nlines - 1;
        }

        kw->text = kw->hdr->array[kw->index + 1];
        for (i = kw->hdr->nlines - 1; i > kw->index; --i)
            memcpy(kw->hdr->array[i + 1],kw->hdr->array[i],81);

        /* Only increment the number of nlines if there were no trailing *
         * blank lines.  If there were, the blank lines are overwritten. */
        if (nblankline == 0)
            kw->hdr->nlines++;

        if (putKwName(kw,nm) == -1) return -1;
        kw->text[8] = '=';
        kw->text[9] = ' ';

        /* Restore nlines (if necessary) to preserve the blank lines */
        if (orig_nlines > kw->hdr->nlines)
            kw->hdr->nlines = orig_nlines;

        return 0;
}

static int addcommentary(Hdr *h, char *text, char *type) {
        int  n, l;
        int  orig_nlines = 0;
        int  nblankline  = 0;
        char *t;
        if (h->nalloc == 0) {
            if (allocHdr(h,HdrUnit, True) != 0)
                return -1;
        }
        if (h->nlines == h->nalloc) {
            if (reallocHdr(h,(h->nalloc + HdrUnit)) != 0)
                return -1;
        }

        /* Look for the last true line of text (i.e., search backwards over *
         * any blank lines) not including the END card.                     */
        nblankline  = searchBlankLine (h);

        /* Save the original number of header lines and reset nlines to the *
         * first trailing blank line.                                       */
        orig_nlines = h->nlines;
        h->nlines   = h->nlines - nblankline;

        n = strlen(text);
        do {
            t = h->array[h->nlines];
            memcpy(t,type,8);
            l = n < 72 ? n : 72;
            memcpy(&t[8],text,l);
            if (l < 72) memcpy(&t[l + 8],BLANK_CARD,(72 - l));
            t[80] = '\0';
            text += l;
            n -= l;
            h->nlines++;
        } while (n);

        /* Restore nlines (if necessary) to preserve the blank lines */
        if (orig_nlines > h->nlines)
            h->nlines = orig_nlines;

        return 0;
}

static FitsKw insertcommentary(FitsKwInfo *kw, char *str, char *type) {
        int n;
        n = strlen(str);
        if (n > 72) n = 72;
        if (insertname(kw,type) == -1) {
                error(BADNAME,keymsg(type)); return NULL; }
        memcpy(&(kw->text[8]),str,n);
        if (n < 72) memcpy(&(kw->text[n + 8]),BLANK_CARD,(72 - n));
        kw->text[80] = '\0';
        return (kw->index == kw->hdr->nlines - 1) ? kw : next(kw);
}

int makePrimaryArrayHdr(Hdr *h, FitsDataType t, int dims, int *ndim) {
        int i, j, n;
        char dig[3];
        char naxisn[9];
        char comm[29] = "Number of values in axis 999";
        h->nlines = 0;
        if (addBoolKw(h,"SIMPLE",True,"Standard FITS file") == -1) return -1;
        switch (t) {
            case FITSBYTE:      n = 8; break;
            case FITSSHORT:     n = 16; break;
            case FITSLONG:      n = 32; break;
            case FITSFLOAT:     n = -32; break;
            case FITSDOUBLE:    n = -64; break;
            default:            error(BADBITPIX,""); return -1;
        }
        if (addIntKw(h,"BITPIX",n,"Number of data bits") == -1) return -1;
        if (addIntKw(h,"NAXIS",dims,"Number of axes") == -1) return -1;
        for (i = 0; i < dims; ++i) {
            strcpy(naxisn,"NAXIS");
            n = i + 1;
            if (n < 0 || n > 999) { error(BADNDIM,""); return -1; }
            for (j = 0; n > 0; ++j, n /= 10)
                dig[j] = n % 10 + '0';
            --j;
            for (n = 5; j >= 0; ++n, --j) {
                naxisn[n] = dig[j];
                comm[20 + n] = dig[j];
            }
            naxisn[n] = '\0';
            comm[20 + n] = '\0';
            if (addIntKw(h,naxisn,ndim[i],comm) == -1)
                return -1;
        }
        if (addBoolKw(h,"EXTEND",True,"Extensions are allowed") == -1)
                return -1;
        return 0;
}

int makeImageExtHdr(Hdr *h, FitsDataType t, int dims, int *ndim,
        char *extname, int extver) {
        int i, j, n;
        char dig[3];
        char naxisn[9];
        char comm[29] = "Number of values in axis 999";
        h->nlines = 0;
        if (addStringKw(h,"XTENSION","IMAGE   ","Image extension") == -1) return -1;
        switch (t) {
            case FITSBYTE:      n = 8; break;
            case FITSSHORT:     n = 16; break;
            case FITSLONG:      n = 32; break;
            case FITSFLOAT:     n = -32; break;
            case FITSDOUBLE:    n = -64; break;
            default:            error(BADBITPIX,""); return -1;
        }
        if (addIntKw(h,"BITPIX",n,"Number of data bits") == -1) return -1;
        if (addIntKw(h,"NAXIS",dims,"Number of axes") == -1) return -1;
        for (i = 0; i < dims; ++i) {
            strcpy(naxisn,"NAXIS");
            n = i + 1;
            if (n < 0 || n > 999) { error(BADNDIM,""); return -1; }
            for (j = 0; n > 0; ++j, n /= 10)
                dig[j] = n % 10 + '0';
            --j;
            for (n = 5; j >= 0; ++n, --j) {
                naxisn[n] = dig[j];
                comm[20 + n] = dig[j];
            }
            naxisn[n] = '\0';
            comm[20 + n] = '\0';
            if (addIntKw(h,naxisn,ndim[i],comm) == -1)
                return -1;
        }
        if (addIntKw(h,"PCOUNT",0,"Standard value for image extensions") == -1)
                return -1;
        if (addIntKw(h,"GCOUNT",1,"Standard value for image extensions") == -1)
                return -1;
        if (addStringKw(h,"EXTNAME",extname,"Name of extension") == -1)
                return -1;
        if (addIntKw(h,"EXTVER",extver,"Standard FITS keyword -- extension number") == -1)
                return -1;
        return 0;
}

FitsKw findKw(Hdr *h, char *name) {
        int i, n;
        char tmp[9];

        n = strlen(name);
        if (n > 8) return NotFound;
        for (i = 0; i < n; ++i)
            tmp[i] = toupper(name[i]);
        tmp[i] = '\0';

        findkw.index = 0;
        return find(h,tmp,&findkw) == 0 ? &findkw : NotFound;
}

FitsKw findnextKw(Hdr *h, char *name) {
        int i, n;
        char tmp[9];

        n = strlen(name);
        if (n > 8) return NotFound;
        for (i = 0; i < n; ++i)
            tmp[i] = toupper(name[i]);
        tmp[i] = '\0';

        ++findkw.index;
        return find(h,tmp,&findkw) == 0 ? &findkw : NotFound;
}

FitsKw insertfirst(Hdr *h) {
        findkw.hdr = h;
        findkw.index = -1;
        findkw.isparsed = False;
        findkw.name[0] = '\0';
        findkw.text = (h->nlines > 0 ? h->array[0] : "");
        findkw.value_end = "";
        findkw.type = FITSNOVALUE;
        return &findkw;
}

FitsKw first(Hdr *h) {
        return next(insertfirst(h));
}

FitsKw next(FitsKw kw_) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        kw->index++;
        if (kw->index >= kw->hdr->nlines) return NotFound;
        kw->name[0] = '\0';
        kw->text = kw->hdr->array[kw->index];
        kw->isparsed = False;
        return kw;
}

FitsKw getKw(Hdr *h, int n) {
        if (n < 0 || n >= h->nlines) return NotFound;
        findkw.hdr = h;
        findkw.index = n;
        findkw.isparsed = False;
        findkw.name[0] = '\0';
        findkw.text = h->array[n];
        return &findkw;
}

char *getKwName(FitsKw kw_) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        int i;
        if (kw->name[0] == '\0') {
            for (i = 0; i < 8 && kw->text[i] != ' ' && kw->text[i] != '='; ++i)
                kw->name[i] = kw->text[i];
            kw->name[i] = '\0';
        }
        return kw->name;
}

char *getKwComm(FitsKw kw_) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        char *s;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return NULL;
        }
        s = kw->value_end;
        ++s;
        while (*s == ' ') ++s;
        if (*s == '/' ) {
            ++s;
            while (*s == ' ') ++s;
        }
        return s;
}

FitsDataType getKwType(FitsKw kw_) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return FITSNOVALUE;
        }
        return kw->type;
}

/*
** High-level keyword access routines written by H. Bushouse and incorporated
** into HSTIO.
**
** These routines read and write keyword values in an image header.
** When writing a keyword value and the keyword doesn't already exist,
** the keyword and comment are added to the header. The routines are:
**
** getKeyB: read a Boolean keyword
** getKeyD: read a double precision keyword
** getKeyF: read a (single precision) floating-point keyword
** getKeyI: read an integer keyword
** getKeyS: read a string keyword
**
** putKeyB: write a Boolean keyword
** putKeyD: write a double precision keyword
** putKeyF: write a (single precision) floating-point keyword
** putKeyI: write an integer keyword
** putKeyS: write a string keyword
**
*/

int getKeyB (Hdr *hdr, char *keyword, Bool *value) {

/* Arguments:
**      hdr     i: pointer to header to be read
**      keyword i: name of keyword
**      value   o: value of keyword
*/

        /* Local variables */
        FitsKw kw;              /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound) {
            return (1);
        } else
            *value = getBoolKw (kw);

        if (hstio_err())
            return (1);

        return (0);
}

int getKeyD (Hdr *hdr, char *keyword, double *value) {

/* Arguments:
**      hdr     i: pointer to header to be read
**      keyword i: name of keyword
**      value   o: value of keyword
*/

        /* Local variables */
        FitsKw kw;              /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound) {
            return (1);
        } else
            *value = getDoubleKw (kw);

        if (hstio_err())
            return (1);

        return (0);
}

int getKeyF (Hdr *hdr, char *keyword, float *value) {

/* Arguments:
**      hdr     i: pointer to header to be read
**      keyword i: name of keyword
**      value   o: value of keyword
*/

        /* Local variables */
        FitsKw kw;              /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound) {
            return (1);
        } else
            *value = getFloatKw (kw);

        if (hstio_err())
            return (1);

        return (0);
}

int getKeyI (Hdr *hdr, char *keyword, int *value) {

/* Arguments:
**      hdr     i: pointer to header to be read
**      keyword i: name of keyword
**      value   o: value of keyword
*/

        /* Local variables */
        FitsKw kw;              /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound) {
            return (1);
        } else
            *value = getIntKw (kw);

        if (hstio_err())
            return (1);

        return (0);
}

int getKeyS (Hdr *hdr, char *keyword, char *value) {

/* Arguments:
**      hdr     i: pointer to header to be read
**      keyword i: name of keyword
**      value   o: value of keyword
*/

        /* Local variables */
        FitsKw kw;              /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound) {
            return (1);
        } else
            getStringKw (kw, value, 70);

        if (hstio_err())
            return (1);

        return (0);
}

int putKeyB (Hdr *hdr, char *keyword, Bool value, char *comment) {

/* Arguments:
**      hdr     i: pointer to header to be updated
**      keyword i: name of keyword
**      value   i: value of keyword
**      comment i: comment to add with keyword if keyword doesn't exist
*/

        FitsKw kw;      /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound)
            addBoolKw (hdr, keyword, value, comment);
        else
            putBoolKw (kw, value);

        if (hstio_err())
            return (1);

        return (0);
}

int putKeyD (Hdr *hdr, char *keyword, double value, char *comment) {

/* Arguments:
**      hdr     i: pointer to header to be updated
**      keyword i: name of keyword
**      value   i: value of keyword
**      comment i: comment to add with keyword if keyword doesn't exist
*/

        FitsKw kw;      /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound)
            addDoubleKw (hdr, keyword, value, comment);
        else
            putDoubleKw (kw, value);

        if (hstio_err())
            return (1);

        return (0);
}

int putKeyF (Hdr *hdr, char *keyword, float value, char *comment) {

/* Arguments:
**      hdr     i: pointer to header to be updated
**      keyword i: name of keyword
**      value   i: value of keyword
**      comment i: comment to add with keyword if keyword doesn't exist
*/

        FitsKw kw;      /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound)
            addFloatKw (hdr, keyword, value, comment);
        else
            putFloatKw (kw, value);

        if (hstio_err())
            return (1);

        return (0);
}

int putKeyI (Hdr *hdr, char *keyword, int value, char *comment) {

/* Arguments:
**      hdr     i: pointer to header to be updated
**      keyword i: name of keyword
**      value   i: value of keyword
**      comment i: comment to add with keyword if keyword doesn't exist
*/

        FitsKw kw;              /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound)
            addIntKw (hdr, keyword, value, comment);
        else
            putIntKw (kw, value);

        if (hstio_err())
            return (1);

        return (0);
}

int putKeyS (Hdr *hdr, char *keyword, char *value, char *comment) {

/* Arguments:
**      hdr     i: pointer to header to be updated
**      keyword i: name of keyword
**      value   i: value of keyword
**      comment i: comment to add with keyword if keyword doesn't exist
*/

        FitsKw kw;              /* keyword pointer */

        kw = findKw (hdr, keyword);
        if (kw == NotFound)
            addStringKw (hdr, keyword, value, comment);
        else
            putStringKw (kw, value);

        if (hstio_err())
            return (1);

        return (0);
}
/*
**
** End high-level keyword access routines.
**
*/

Bool   getBoolKw(FitsKw kw_) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return False;
        }
        if (kw->type == FITSLOGICAL)
            return kw->bresult;
        error(BADFITSTYPE,keymsg(kw->name));
        return False;
}

int    getIntKw(FitsKw kw_) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        double value;

        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return 0;
        }
        if (kw->type == FITSLOGICAL)
            return (int)(kw->bresult);
        else if (kw->type == FITSLONG) {
            /* TODO: Implicit declaration */
            if (checkRange (T_INTEGER, kw->nresult)) {
                error(BADFITSTYPE,keymsg(kw->name));
                return 0;
            }
            return (int)(kw->nresult.data.l);
        }
        else if (kw->type == FITSDOUBLE) {
            value = fabs(kw->nresult.data.d);
            if ((value - (int)value) > 0.0) {
                error(BADFITSTYPE,keymsg(kw->name));
                return 0;
            } else {
                if (checkRange (T_INTEGER, kw->nresult)) {
                   error(BADFITSTYPE,keymsg(kw->name));
                   return 0;
                }
                return ((int)kw->nresult.data.d);
              }
        }
        error(BADFITSTYPE,keymsg(kw->name));
        return 0;
}

float  getFloatKw(FitsKw kw_) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return 0;
        }
        if (kw->type == FITSDOUBLE) {
            if (checkRange (T_FLOAT, kw->nresult)) {
                error(BADFITSTYPE,keymsg(kw->name));
                return 0.0F;
            }
            return (float)kw->nresult.data.d;
        }
        else if (kw->type == FITSLONG) {
            if (checkRange (T_FLOAT, kw->nresult)) {
                error(BADFITSTYPE,keymsg(kw->name));
                return 0.0F;
            }
            return (float)kw->nresult.data.l;
        }
        error(BADFITSTYPE,keymsg(kw->name));
        return 0.0F;
}

double getDoubleKw(FitsKw kw_) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return 0;
        }
        if (kw->type == FITSDOUBLE)
            return kw->nresult.data.d;
        else if (kw->type == FITSLONG)
            return (double)kw->nresult.data.l;
        error(BADFITSTYPE,keymsg(kw->name));
        return 0.0;
}

int   getStringKw(FitsKw kw_, char *str, int maxch) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return 0;
        }
        if (kw->type == FITSCHAR) {
            if ((int)strlen(kw->cresult) < maxch)
                strcpy(str,kw->cresult);
            else
                { strncpy(str,kw->cresult,maxch); str[maxch] = '\0'; }
            return 0;
        }
        error(BADFITSTYPE,keymsg(kw->name));
        return -1;
}

int putKwName(FitsKw kw_, char *nm) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        int i;
        char tmp[9];

        int n = strlen(nm);
        if (n > 8 ) return -1;

        /* Ensure that the keyword is upper case */
        for (i = 0; i < n; ++i)
            tmp[i] = toupper(nm[i]);
        tmp[i] = '\0';

        memcpy(kw->text,tmp,n);
        if (n < 8) memcpy(&(kw->text[n]),BLANK_CARD,(8 - n));
        return 0;
}

void putKwComm(FitsKw kw_, char *comm) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        int n, len, pos;
        if (comm == NULL)
            len = 0;
        else
            len = strlen(comm);
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return;
        }
        for (pos = 0; pos < 80; ++pos) /* find position of end of value */
            if (kw->value_end == &(kw->text[pos])) break;
        ++pos;
        if (pos < 30) {
            memcpy(&(kw->text[pos]), BLANK_CARD, (30 - pos - 1));
            pos = 30;
        }
        if (pos > (80-3)) {
            memcpy(&(kw->text[pos]), BLANK_CARD, 80-pos);
            return;
        }
        kw->text[pos++] = ' ';
        kw->text[pos++] = '/';
        kw->text[pos++] = ' ';
        n = 80 - pos;
        if (len > n) len = n;
        memcpy(&(kw->text[pos]), comm, len);
        if (len < n)
            memcpy(&(kw->text[len + pos]), BLANK_CARD, (n - len));
}

int   putBoolKw(FitsKw kw_, Bool bval){
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return -1;
        }
        memcpy(&(kw->text[10]),BLANK_CARD,20);
        kw->text[29] = (bval == True ? 'T' : 'F');
        kw->value_end = &(kw->text[29]);
        return 0;
}

int   putIntKw(FitsKw kw_, int lval) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return -1;
        }
        memcpy(&(kw->text[10]),BLANK_CARD,20);
        sprintf(&(kw->text[18]),"%12d",lval);
        kw->text[30] = ' ';
        kw->value_end = &(kw->text[29]);
        return 0;
}

int   putFloatKw(FitsKw kw_, float fval) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return -1;
        }
        memcpy(&(kw->text[10]),BLANK_CARD,20);
        sprintf(&kw->text[16],"%#14.7E",fval);
        kw->text[30] = ' ';
        kw->value_end = &(kw->text[29]);
        return 0;
}

int   putDoubleKw(FitsKw kw_, double dval) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        int i;
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return -1;
        }
        memcpy(&(kw->text[10]),BLANK_CARD,20);
        sprintf(&kw->text[10],"%#20.12E",dval); /* optimum %23.15E */
        for (i = 29; i < 10; --i) /* change the E to a D */
            if (kw->text[i] == 'E') { kw->text[i] = 'D'; break; }
        kw->text[30] = ' ';
        kw->value_end = &(kw->text[29]);
        return 0;
}

int   putStringKw(FitsKw kw_, char *txt) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        char *s;
        char tmp[81];
        tmp[0] = '\0';
        s = getKwComm(kw_);
        strcpy (tmp, s);
        /* TODO: Implicit declaration */
        if (putString(kw_,txt) == -1) return 0;
        putKwComm(kw_, tmp);
        kw->text[80] = '\0';
        return 0;
}

int   putString(FitsKw kw_, char *txt) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        int i, j, length;
        char *p;
        char newtxt[81];
        if (kw->isparsed == False) {
            if (getvalue(kw) == -1) return -1;
        }
        kw->text[10] = '\'';
        p = &kw->text[11];

        /* Determine if the input text contains embedded single quotes. *
         * Represent a single quote within the text as two successive   *
         * single quotes to conform to the FITS standard.               */
        length = strlen (txt);
        for (i = 0, j = 0; i < length; i++, j++) {
             newtxt[j] = txt[i];
             if (txt[i] == '\'') {
                 j++;
                 newtxt[j] = '\'';
             }
        }
        newtxt[j] = '\0';

        for (i = 0; i < 68 && newtxt[i] != '\0'; ++i)
            *p++ = newtxt[i];
        for (; i < 8; ++i)
            *p++ = ' ';
        *p++ = '\'';
        for (i += 2; i < 20; ++i, ++p)
            *p = ' ';
        --p;
        kw->value_end = p;
        return 0;
}

int   addBoolKw(Hdr *h, char *name, Bool value, char *comment) {
        FitsKwInfo kw;
        if (h->nlines != 0 &&
            strncmp(h->array[h->nlines - 1],"END     ",8) == 0) --(h->nlines);
        kw.index = h->nlines - 1;
        kw.hdr = h;
        return insertBoolKw(&kw,name,value,comment) == NULL ? -1 : 0;
}

int   addIntKw(Hdr *h, char *name, int value, char *comment) {
        FitsKwInfo kw;
        if (h->nlines != 0 &&
            strncmp(h->array[h->nlines - 1],"END     ",8) == 0) --(h->nlines);
        kw.index = h->nlines - 1;
        kw.hdr = h;
        return insertIntKw(&kw,name,value,comment) == NULL ? -1 : 0;
}

int   addFloatKw(Hdr *h, char *name, float value, char *comment) {
        FitsKwInfo kw;
        if (h->nlines != 0 &&
            strncmp(h->array[h->nlines - 1],"END     ",8) == 0) --(h->nlines);
        kw.index = h->nlines - 1;
        kw.hdr = h;
        return insertFloatKw(&kw,name,value,comment) == NULL ? -1 : 0;
}

int   addDoubleKw(Hdr *h, char *name, double value, char *comment) {
        FitsKwInfo kw;
        if (h->nlines != 0 &&
            strncmp(h->array[h->nlines - 1],"END     ",8) == 0) --(h->nlines);
        kw.index = h->nlines - 1;
        kw.hdr = h;
        return insertDoubleKw(&kw,name,value,comment) == NULL ? -1 : 0;
}

int   addStringKw(Hdr *h, char *name, char *value, char *comment) {
        FitsKwInfo kw;
        if (h->nlines != 0 &&
            strncmp(h->array[h->nlines - 1],"END     ",8) == 0) --(h->nlines);
        kw.index = h->nlines - 1;
        kw.hdr = h;
        return insertStringKw(&kw,name,value,comment) == NULL ? -1 : 0;
}

int   addSpacesKw(Hdr *h, char *text) {
        if (h->nlines != 0 &&
            strncmp(h->array[h->nlines - 1],"END     ",8) == 0) --(h->nlines);
        return addcommentary(h,text,"        ");
}

int   addCommentKw(Hdr *h, char *text) {
        if (h->nlines != 0 &&
            strncmp(h->array[h->nlines - 1],"END     ",8) == 0) --(h->nlines);
        return addcommentary(h,text,"COMMENT ");
}

int   addHistoryKw(Hdr *h, char *text) {
        if (h->nlines != 0 &&
            strncmp(h->array[h->nlines - 1],"END     ",8) == 0) --(h->nlines);
        return addcommentary(h,text,"HISTORY ");
}

FitsKw  insertBoolKw(FitsKw kw_, char *name, Bool value, char *comment) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        kw->type = FITSLOGICAL;
        kw->isparsed = True;
        if (insertname(kw,name) == -1) {
                error(BADNAME,keymsg(name)); return NULL; }
        if (putBoolKw(kw,value) == -1) return NULL;
        kw->text[80] = '\0';
        putKwComm(kw,comment);
        return (kw->index >= kw->hdr->nlines - 1) ? kw : next(kw);
}

FitsKw  insertIntKw(FitsKw kw_, char *name, int value, char *comment) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        kw->type = FITSLONG;
        kw->isparsed = True;
        if (insertname(kw,name) == -1) {
                error(BADNAME,keymsg(name)); return NULL; }
        if (putIntKw(kw,value) == -1) return NULL;
        putKwComm(kw,comment);
        kw->text[80] = '\0';
        return (kw->index == kw->hdr->nlines - 1) ? kw : next(kw);
}

FitsKw  insertFloatKw(FitsKw kw_, char *name, float value, char *comment) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        kw->type = FITSFLOAT;
        kw->isparsed = True;
        if (insertname(kw,name) == -1) {
                error(BADNAME,keymsg(name)); return NULL; }
        if (putFloatKw(kw,value) == -1) return NULL;
        putKwComm(kw,comment);
        kw->text[80] = '\0';
        return (kw->index == kw->hdr->nlines - 1) ? kw : next(kw);
}

FitsKw  insertDoubleKw(FitsKw kw_, char *name, double value, char *comment) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        kw->type = FITSDOUBLE;
        kw->isparsed = True;
        if (insertname(kw,name) == -1) {
                error(BADNAME,keymsg(name)); return NULL; }
        if (putDoubleKw(kw,value) == -1) return NULL;
        putKwComm(kw,comment);
        kw->text[80] = '\0';
        return (kw->index == kw->hdr->nlines - 1) ? kw : next(kw);
}

FitsKw  insertStringKw(FitsKw kw_, char *name, char *value, char *comment) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        kw->type = FITSCHAR;
        kw->isparsed = True;
        if (insertname(kw,name) == -1) {
                error(BADNAME,keymsg(name)); return NULL; }
        if (putString(kw,value) == -1) return NULL;
        putKwComm(kw,comment);
        kw->text[80] = '\0';
        return (kw->index == kw->hdr->nlines - 1) ? kw : next(kw);
}

FitsKw  insertSpacesKw(FitsKw kw_, char *text) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        return insertcommentary(kw,text,"        ");
}

FitsKw  insertCommentKw(FitsKw kw_, char *text) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        return insertcommentary(kw,text,"COMMENT ");
}

FitsKw  insertHistoryKw(FitsKw kw_, char *text) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        return insertcommentary(kw,text,"HISTORY ");
}

int     addFitsCard(Hdr *h, char *card) {
        FitsKwInfo kw;
        if (h->nlines != 0 &&
            strncmp(h->array[h->nlines - 1],"END     ",8) == 0) --(h->nlines);
        kw.index = h->nlines - 1;
        kw.hdr = h;
        return insertFitsCard(&kw,card) == NULL ? -1 : 0;
}

FitsKw  insertFitsCard(FitsKw kw_, char *card) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        int len, i;
        int nblankline  = 0;
        int orig_nlines = 0;

        kw->isparsed = False;

        /* make sure we have room */
        if (kw->hdr->nalloc == 0) {
            if (allocHdr(kw->hdr,HdrUnit, True) != 0)
                { error(NOMEM,""); return NULL; }
        }
        if (kw->hdr->nlines == kw->hdr->nalloc) {
            if (reallocHdr(kw->hdr,(kw->hdr->nalloc + HdrUnit)) != 0)
                { error(NOMEM,""); return NULL; }
        }

        /* Look for the last true line of text (i.e., search backwards over *
         * any blank lines) not including the END card.                     */
        nblankline = searchBlankLine (kw->hdr);

        /* Save the original number of header lines and reset nlines to the *
         * first trailing blank line.                                       */
        if (nblankline != 0) {
            orig_nlines     = kw->hdr->nlines;
            kw->hdr->nlines = kw->hdr->nlines - nblankline;
            if (kw->index == orig_nlines - 1)
                kw->index = kw->hdr->nlines - 1;
        }

        /* make room in the array for the inserted card */
        kw->text = kw->hdr->array[kw->index + 1];
        for (i = kw->hdr->nlines - 1; i > kw->index; --i)
            memcpy(kw->hdr->array[i + 1],kw->hdr->array[i],81);

        /* Only increment the number of nlines if there were no trailing *
         * blank lines.  If there were, the blank lines are overwritten. */
        if (nblankline == 0)
            kw->hdr->nlines++;

        /* move the card into the header array */
        len = strlen(card);
        if (len >= 80)
            strncpy(kw->text,card,80);
        else {
            strcpy(kw->text,card);
            for (i = len; i < 80; ++i)
                kw->text[i] = ' ';
        }
        kw->text[80] = '\0';

        /* Restore nlines (if necessary) to preserve the blank lines */
        if (orig_nlines > kw->hdr->nlines)
            kw->hdr->nlines = orig_nlines;

        return (kw->index == kw->hdr->nlines - 1) ? kw : next(kw);
}

void delKw(FitsKw kw_) {
        FitsKwInfo *kw = (FitsKwInfo *)kw_;
        int i;
        for (i = kw->index + 1; i < kw->hdr->nlines; ++i)
            memcpy(kw->hdr->array[i - 1],kw->hdr->array[i],81);
        kw->isparsed = False;
        kw->name[0] = '\0';
        kw->hdr->nlines--;
}

void delAllKw(Hdr *h) {
        h->nlines = 0;
}

# define LONG   1
# define DOUBLE 2

/**********************************************************************
**
** checkRange - Compare the valid range of the requested datatype
** to the data value read to make sure the value can be accommodated.
** It is only necessary to check for casting the data which has been
** read as a LONG or a DOUBLE and being demoted to an INT or a FLOAT.
** This routine checks the range BEFORE the conversion would be done.
**
** Michele De La Pena  10 December 1998
**
**********************************************************************/

int
checkRange (TargetDataType target, NumericResult result) {

    switch (target) {

         case T_INTEGER:
             if (result.type == LONG) {
                if ((result.data.l > INT_MAX) || (result.data.l < INT_MIN))
                   return (1);
             } else if (result.type == DOUBLE) {
                if ((result.data.d > INT_MAX) || (result.data.d < INT_MIN))
                   return (1);
             }
             break;

         case T_FLOAT:
             if (result.type == LONG) {
                if (fabs((double)result.data.l) > FLT_MAX)
                   return (1);
             } else if (result.type == DOUBLE) {
                if (fabs(result.data.d) > FLT_MAX)
                   return (1);
             }
             break;

         default:
             fprintf (stderr, "Invalid target datatype.\n");
             return (1);
             break;
    }

    return (0);
}
