/*
** This version incorporates both UNIX and VMS versions.
*/

# include <stdlib.h>
# include <string.h>
# include <c_iraf.h>

void c_irafinit (int argc, char **argv) {
        ;
}

/* Convert a virtual filename (path) to an OS-specific filename

This handles both the IRAF and UNIX forms of variable substitution:

  - name$str
  - $name/str

  vfn: input virtual filename
  osfn: output OS filename (must be a buffer of size SZ_PATHNAME)
*/
#define SZ_PATHNAME 1024
int c_vfn2osfn(const char *const vfn, char *const osfn) {
        size_t vfn_len;
        const char *ip;
        char *op;
        char fname[SZ_PATHNAME+1];
        char *ldir;
        size_t ldir_len;

        vfn_len = strlen(vfn);

        /* Recursively expand logical directories, but don't do
         * anything about subdirectories, extensions, etc.  This is
         * all that is needed for UNIX.
         */
        if (vfn[0] == '$') {
            /* Handle UNIX-style */
          for (ip=vfn+1, op=fname; (*op = *ip++); ++op) {
                if (*op == '/') {
                    *op = 0;
                    if ((ldir = getenv(fname))) {
                        ldir_len = strlen(ldir);
                        if (ldir_len > SZ_PATHNAME) {
                            return -1;
                        }
                        strcpy(fname, ldir);
                    } else {
                        ldir_len = 0;
                    }
                    if (ldir_len + vfn_len - (ip - vfn) > SZ_PATHNAME) {
                        return -1;
                    }
                    strcat(fname, ip - 1);
                    return c_vfn2osfn(fname, osfn);
                }
            }
        } else {
            /* Handle IRAF-style */
          for (ip=vfn, op=fname; (*op = *ip++); ++op) {
                if (*op == '$') {
                    *op = 0;
                    if ((ldir = getenv(fname))) {
                        ldir_len = strlen(ldir);
                        if (ldir_len > SZ_PATHNAME) {
                            return -1;
                        }
                        strcpy(fname, ldir);
                    } else {
                        ldir_len = 0;
                    }
                    if (ldir_len + vfn_len - (ip - vfn) > SZ_PATHNAME) {
                        return -1;
                    }
                    strcat(fname, ip);
                    return c_vfn2osfn(fname, osfn);
                }
            }
        }

        /* Copy filename to the output string.  Fix up the "//"
         * sequences that occur because IRAF likes the / at the end of
         * logical directory names.
         */
        for (ip = fname, op = osfn; (*op = *ip++); ++op) {
            if (*op == '/' && op > osfn && *(op - 1) == '/') {
                --op;
            }
        }

        return 0;
}
