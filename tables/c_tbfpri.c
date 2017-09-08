# include <string.h>
# include <fitsio.h>
#include "hstcal.h"
# include "ctables.h"

void c_tbfpri (char *intable, char *outtable, int *copied) {

/* Copy the primary header from intable to outtable.  If outtable already
   exists, nothing will be done.  If the primary HDU of intable contains
   a data block (i.e. it's not just a header), nothing will be done.
arguments:
char *intable           i: name of the input file
char *outtable          i: name of the output file
int *copied             o: set to true (1) if the primary header was actually
                           copied, else false (0); it is considered an error
                           if the header was not copied due to one of the
                           conditions mentioned above
*/

        fitsfile *infptr, *outfptr;
        char *infile, *outfile;         /* file names */
        char *filename;                 /* for updating keyword in output */
        int naxis;
        int i;
        int status = 0;

        *copied = 0;            /* initial value */

        outfile = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
        status = c_vfn2osfn (outtable, outfile);
        if (status != 0) {
            setError (status, "c_tbfpri:  error from c_vfn2osfn");
            free (outfile);
            return;
        }
        for (i = strlen (outfile) - 1;  i >= 0;  i--) { /* remove brackets */
            if (outfile[i] == '[') {
                outfile[i] = '\0';
                break;
            }
        }

        /* if the output file already exists, there's nothing to do */
        if (checkExists (outtable)) {
            free (outfile);
            return;
        }

        infile = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
        status = c_vfn2osfn (intable, infile);
        if (status != 0) {
            setError (status, "c_tbfpri:  error from c_vfn2osfn");
            free (outfile);
            free (infile);
            return;
        }
        for (i = strlen (infile) - 1;  i >= 0;  i--) {  /* remove brackets */
            if (infile[i] == '[') {
                infile[i] = '\0';
                break;
            }
        }

        /* fits_open_file = ffopen */
        fits_open_file (&infptr, infile, READONLY, &status);
        free (infile);
        if (status != 0) {
            setError (status, "c_tbfpri:  couldn't open input file");
            free (outfile);
            return;
        }
        /* fits_read_key = ffgky */
        fits_read_key (infptr, TINT, "NAXIS", &naxis, NULL, &status);
        if (status != 0) {
            setError (status, "c_tbfpri:  couldn't read NAXIS keyword");
            status = 0;
            /* fits_close_file = ffclos */
            fits_close_file (infptr, &status);
            free (outfile);
            return;
        }
        /* the input primary HDU must be just a header */
        if (naxis > 0) {
            fits_close_file (infptr, &status);
            free (outfile);
            return;
        }

        /* fits_create_file = ffinit */
        fits_create_file (&outfptr, outfile, &status);
        if (status != 0) {
            setError (status, "c_tbfpri:  couldn't create output file");
            status = 0;
            fits_close_file (infptr, &status);
            free (outfile);
            return;
        }

        /* fits_copy_header = ffcphd */
        fits_copy_header (infptr, outfptr, &status);
        if (status == 0) {
            *copied = 1;
        } else {
            setError (status, "c_tbfpri:  couldn't copy primary header");
            status = 0;
        }

        /* add or update FILENAME in the output primary header */
        filename = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
        for (i = strlen (outfile) - 1;  i >= 0;  i--) {
            if (outfile[i] == '/')
                break;
        }
        if (i < 0)
            i = 0;
        strcpy (filename, outfile+i);
        fits_update_key (outfptr, TSTRING, "FILENAME", filename,
                            "File name", &status);
        free (filename);

        /* fits_write_date = ffpdat */
        fits_write_date (outfptr, &status);     /* write the DATE keyword */

        fits_close_file (outfptr, &status);
        fits_close_file (infptr, &status);
        free (outfile);
        if (status != 0)
            setError (status, "c_tbfpri:  error closing file");
}
