/* This contains functions to process IMPHTTAB and calculate PHOTXXXX
   values for PHOTCORR step in CALXXX.

History:

MLS, 12/2013 Updated to account for new WFC3 UVIS PHTFLAM1 and PHTFLAM2,
             which were added to the PhotPars setup to be used with
             WFC3's new IMPHTTAB, which has those 2 added extensions.
             This should not break use of the other instruments.

PLL, 12/2013 Allow MJD extrapolation using simple straight line from
             last 2 data points, if and only if EXTRAP exists in
             IMPHTTAB primary header and is set to True.

MLS: 07/2015 Cleand up for unused variables and warning
MLS: 08/2015 Added some initializations the clang complained about

*/
# include <stdio.h>
# include <string.h>
# include <stdlib.h>    /* for strtol function */
# include <ctype.h>
# include <math.h>

# include "hstio.h"
# include "xtables.h"
# include "imphttab.h"
# include "c_iraf.h"  /* For Bool type */

/* Internal functions to be used to interpret IMPHTTAB ref tables */
static int OpenPhotTab (char *, char *, PhtCols *);
static int ReadPhotTab (PhtCols *, int, PhtRow *);
static int ReadPhotArray (PhtCols *, int, PhtRow *);
static double ComputeValue (PhtRow *, PhotPar *);
static int ClosePhotTab (PhtCols *);
static int InterpretPhotmode(char *, PhotPar *);
static int CompareObsModes(char *obsmode1, char *obsmode2);
static int PhotRowPedigree (PhotPar *, int , IRAFPointer, IRAFPointer, IRAFPointer);
static void ClosePhotRow (PhtRow *);

static char *photnames[6] = {"PHOTZPT","PHOTFLAM","PHOTPLAM","PHOTBW","PHTFLAM1","PHTFLAM2"};

/* This routine gets the values from IMPHTTAB.

   The table should contain the following:

       header parameters:
           EXTRAP:   set to True to allow MJD extrapolation
           DESCRIP:  description of the source of the values
           PEDIGREE: description of the pedigree of the values

       columns:
           OBSMODE:      observation mode for exposure
           DATACOL:      name of column which tneeds to be read to get result
           NELEM<n>:     number of elements for parameter <n>
           PAR<n>VALUES: values for parameterized variable
           <result>:     column with single result value
           <result><n>:  column with result based on parameterized variable <n>

   The table is read to find the row for which the OBSMODE matches the
   generated PHOTMODE string (from the image header). There will be one
   table extension for each photometry <result>:

       PHOTFLAM
       PHOTPLAM
       PHOTBW
       PHTFLAM<n> (WFC3 only)
*/
int GetPhotTab (PhotPar *obs, char *photmode) {
    /* arguments:
       RefTab  *obs      io: name and pedigree values
       char    *photmode  i: photmode from science header for row selection

       Calling routine will need to initialize and free PhotPar object using:

           RefTab *ref;
           PhotPar *obs;
           InitPhotPar(obs, ref->name, ref->pedigree);

           FreePhotPar(obs);

       where ref is the information about the reference file from the header
    */

    extern int status;

    PhtCols tabinfo;    /* pointer to table descriptor, etc */
    PhtRow tabrow;        /* values read from a table row */

    int row;        /* loop index */
    int extn;
    char phdrname[SZ_FNAME];
    IODescPtr im;        /* descriptor for primary header unit */
    Hdr tphdr;        /* primary header */
    FitsKw key;        /* location of keyword in header */
    char pname[SZ_COLNAME];
    double value;
    int numkeys;

    int foundit;        /* has the correct row been found? */
    int PhotRowPedigree (PhotPar *, int, IRAFPointer, IRAFPointer, IRAFPointer);
    int SameInt (int, int);
    int SameFlt (float, float);
    int streq_ic (char *, char *);
    double ComputeValue(PhtRow *, PhotPar *);

    extern char *photnames[];

    /* initialize status for this run */
    status = 0;

    /* Interpret OBSMODE string from science file header for
       comparison with obsmode values in reference table
       primarily by stripping out any parameterized values */
    /* Find the length of the string to be interpreted */
    status = InterpretPhotmode(photmode, obs);
    if (status != PHOT_OK) {
        printf("*** Error in InterpretPhotmode\n");
        return(status=ERROR_RETURN);
    }

    /* Create name of table to access PRIMARY header keywords */
    strcpy(phdrname, obs->name);
    /*strcat(phdrname,"[0]"); */

    /* Open PRIMARY header for reading keywords */
    initHdr (&tphdr);
    /* Open the primary header of the reference file. */
    im = openInputImage (phdrname, "", 0);
    if (hstio_err()) {
        printf ("IMPHTTAB `%s' not found.\n", obs->name);
        clear_hstioerr();
        status = OPEN_FAILED;
        return (status);
    }
    getHeader (im, &tphdr);
    if (hstio_err()){
        printf ("\n==>ERROR: IMPHTTAB `%s' not found.\n", obs->name);
        return (status = HEADER_PROBLEM);
    }
    /* Read in keywords from PRIMARY header for use in interpreting
       the tables */
    key = findKw (&tphdr, "PARNUM");
    if (key == NotFound) {
        printf ("\n==>ERROR: Trying to get PARNUM...\n");
        closeImage(im);
        return (status = KEYWORD_MISSING);
    } else {
        tabinfo.parnum = getIntKw(key);
    }


    /* Read in PHOTZPT keyword value from Primary header */
    key = findKw (&tphdr, "PHOTZPT");
    obs->photzpt = getDoubleKw (key);
    if (hstio_err()){
        printf("\n==>ERROR: Keyword `PHOTZPT` not found in PRIMARY header.\n");
        closeImage(im);
        return (status = KEYWORD_MISSING);
    }
    key = findKw (&tphdr, "NEXTEND");
    numkeys = getIntKw (key);
    if (hstio_err()){
        printf("\n==>ERROR: Keyword `NEXTEND` not found in PRIMARY header.\n");
        closeImage(im);
        return (status = KEYWORD_MISSING);
    }

    /* Read in EXTRAP keyword value from Primary header */
    key = findKw (&tphdr, "EXTRAP");
    if (key == NotFound) {
        obs->extrap = False;
    } else {
        obs->extrap = getBoolKw (key);
        if (hstio_err()){
            printf("\n==>ERROR: Keyword `EXTRAP` found in PRIMARY header but unable to be extracted.\n");
            closeImage(im);
            return (status = HEADER_PROBLEM);
        }
    }

    /* Close PRIMARY header of reference table */
    closeImage(im);
    freeHdr (&tphdr);

    /* Now step through each of the extensions to compute the
       different photometry keyword values, one from each
       extension.
     */
    for (extn = 1; extn <= numkeys; extn++){
        strcpy(pname,photnames[extn]);
        /* Open the photometry parameters table and find columns. */
        if (OpenPhotTab (obs->name, pname, &tabinfo)) {
            printf("*** Error in OpenPhotTab %d\n",status);
            return (status);
        }
        /* Check each row for a match with obsmode,
           and get info from the matching row.
         */

        foundit = 0;
        for (row = 1;  row <= tabinfo.nrows;  row++) {

            /* Read the current row into tabrow. */
            if (ReadPhotTab (&tabinfo, row, &tabrow)) {
                printf("*** Error in ReadPhotTab\n");
                return (status);
            }

            if (CompareObsModes(tabrow.obsmode, obs->obsmode) == PHOT_OK) {
                foundit = 1;
                if (PhotRowPedigree (obs, row,
                            tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip))
                    return (status);

                if (obs->goodPedigree == DUMMY_PEDIGREE) {
                    printf ("==>Warning: Row %d of IMPHTTAB is DUMMY.\n", row);
                    /* trlwarn (MsgText); */
                }

                /* Read in photometry values from table row */
                if ( (status = ReadPhotArray(&tabinfo, row, &tabrow))) {
                    printf("*** Error in ReadPhotArray\n");
                    return (status);
                }

                /* Interpret values to compute return value.
                   value = 0.0 is allowed.
                 */
                value = ComputeValue (&tabrow, obs);
                if ((value != 0.0) && (value == '\0')){
                    return(status=INTERNAL_ERROR);
                }
                printf("==> Value of %s = %0.8g\n",photnames[extn], value);

                /* Free memory used to read in this row */
                ClosePhotRow(&tabrow);

                break;
            } /* End of if CompareObsModes() */
        } /* End of loop over table rows */

        if (foundit == 0) {
            printf ("\n==>ERROR: Matching row not found in IMPHTTAB `%s'.\n", obs->obsmode);
            /*trlerror (MsgText); */
            printf ("\n==>ERROR: OBSMODE %s\n",    obs->obsmode);
            /*trlerror (MsgText); */

            ClosePhotTab (&tabinfo);
            return (status = TABLE_ERROR);
        }

        if (ClosePhotTab (&tabinfo))        /* close the table */
            return (status);

        /* Store computed value for keyword as appropriate member in
           PhotPar struct */
        if (strncmp(pname,"PHOTFLAM",8) == 0) {
            obs->photflam = value;
        } else if (strncmp(pname,"PHOTPLAM",8) == 0) {
            obs->photplam = value;
        } else if (strncmp (pname,"PHOTBW",6) == 0){
            obs->photbw = value;
        } else if (strncmp (pname,"PHTFLAM1",8) == 0) {
            obs->phtflam1 = value;
        } else if (strncmp (pname,"PHTFLAM2",8) == 0) {
            obs->phtflam2 = value;
        } else {
            return (status = TABLE_ERROR);
        }

    } /* End of loop over photometry keyword names */

    return (status);
}


/* This routine opens the CCD parameters table, finds the columns that we
   need, and gets the total number of rows in the table.  The columns are
   CCDAMP, CCDGAIN, CCDBIAS, and READNSE.
 */
static int OpenPhotTab (char *tabname, char *photvar, PhtCols *tabinfo) {

    extern int status;

    char tname[SZ_FNAME];
    char **colnames, **ecolnames, **pncolnames, **pvcolnames;

    int *nocol;
    int i, j, missing;
    int parnum;
    
    int PrintMissingCols_IMPHTTAB (int, int, int *, char **, char *, IRAFPointer);
    int buildTabName (char *, char *, char *);

    /* allocate space for column names to be read from this table */
    colnames = (char **)calloc(tabinfo->parnum+1, sizeof(char *));
    ecolnames = (char **)calloc(tabinfo->parnum+1, sizeof(char *));
    pncolnames = (char **)calloc(tabinfo->parnum+1, sizeof(char *));
    pvcolnames = (char **)calloc(tabinfo->parnum+1, sizeof(char *));
    for (i=0;i <= tabinfo->parnum; i++){
        colnames[i] = (char *)calloc(SZ_COLNAME, sizeof(char));
        ecolnames[i] = (char *)calloc(SZ_COLNAME, sizeof(char));
        pncolnames[i] = (char *)calloc(SZ_COLNAME, sizeof(char));
        pvcolnames[i] = (char *)calloc(SZ_COLNAME, sizeof(char));
    }
    /* Copy in the column names now */
    strcpy(colnames[0],photvar);
    for (parnum=1;parnum <= tabinfo->parnum; parnum++){
        sprintf(colnames[parnum],"%s%i",photvar,parnum);
        sprintf(ecolnames[parnum],"NELEM%i",parnum);
        sprintf(pncolnames[parnum],"PAR%iNAMES",parnum);
        sprintf(pvcolnames[parnum],"PAR%iVALUES",parnum);
    }

    /* Create name of table with extension to be opened */
    sprintf(tname,"%s[%s]",tabname,photvar);

    /* keep track of what extension we are processing here */
    strcpy(tabinfo->photvar, photvar);

    tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
    if (c_iraferr()) {
        printf ("\n==>ERROR: IMPHTTAB extension`%s' not found.\n", tname);
        /*trlerror (MsgText); */
        return (status = OPEN_FAILED);
    }

    /* Start getting info on this table */
    tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);
    tabinfo->ncols = c_tbpsta (tabinfo->tp, TBL_NCOLS);

    /* Initialize missing column array */
    nocol = (int *)calloc(tabinfo->ncols+1,sizeof(int));
    for (j = 0; j < tabinfo->ncols; j++){
        nocol[j] = NO;
    }

    /* Find the columns. */
    c_tbcfnd1 (tabinfo->tp, "OBSMODE", &tabinfo->cp_obsmode);
    c_tbcfnd1 (tabinfo->tp, "DATACOL", &tabinfo->cp_datacol);
    c_tbcfnd1 (tabinfo->tp, photvar,   &tabinfo->cp_result[0]);
    missing = 0;
    i=0;
    for (parnum = 1;parnum <= tabinfo->parnum; parnum++){
        /*
           Read in NELEM<parnum> columns
         */
        c_tbcfnd1 (tabinfo->tp, ecolnames[parnum], &tabinfo->cp_nelem[parnum]);
        /*
           Read in PAR<parnum>NAMES columns
         */
        c_tbcfnd1 (tabinfo->tp, pncolnames[parnum], &tabinfo->cp_parnames[parnum]);
        /*
           Read in PAR<parnum>VALUES columns
         */
        c_tbcfnd1 (tabinfo->tp, pvcolnames[parnum], &tabinfo->cp_parvalues[parnum]);
        /*
           Read in the PHOT*<parnum> columns
         */
        c_tbcfnd1 (tabinfo->tp, colnames[parnum], &tabinfo->cp_result[parnum]);
        if (tabinfo->cp_result[parnum] == 0 ) { missing++; nocol[i] = YES;} i++;
    }
    /* Initialize counters here... */

    /* Increment i for every column, mark only missing columns in
       nocol as YES.  WJH 27 July 1999
     */
    if (tabinfo->cp_obsmode == 0 ) { missing++; nocol[i] = YES;} i++;
    if (tabinfo->cp_datacol == 0 ) { missing++; nocol[i] = YES;} i++;

    if (PrintMissingCols_IMPHTTAB(missing, tabinfo->ncols, nocol, colnames,
                "IMPHTTAB", tabinfo->tp) )
        return(status);

    /* Pedigree and descrip are optional columns. */
    c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
    c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

    /* Free memory */
    for (i=0; i <= tabinfo->parnum; i++){
        free(colnames[i]);
        free(ecolnames[i]);
        free(pncolnames[i]);
        free(pvcolnames[i]);
    }
    free(colnames);
    free(ecolnames);
    free(pncolnames);
    free(pvcolnames);
    free(nocol);

    return (status);
}


/* Interprets PHOTMODE string (as a comma-separated, all-lowercase string)
   from the science image header into obsmode for comparison with
   obsmode column values in IMPHTTAB while also recording the parameterized
   variable values from the PHOTMODE string as outputs parnames, parvalues.
*/
static int InterpretPhotmode(char *photmode, PhotPar *obs){
    int i;
    int numpar, n;
    int tok;
    int indx;
    char *loc;
    char *chtok,*chval,*chtmp;
    char **obsnames;
    int obselems;

    extern int status;

    /* Temporary string for parameterized variable value from photmode */
    char tempmode[SZ_FITS_REC];

    int AllocPhotPar(PhotPar *, int);

    numpar = 0;
    n=0;
    
    /* scan entire photmode string and count how many # symbols are found */
    obselems = 0;
    strcpy(tempmode,photmode);
    /* get first token from the string */
    chtok = strtok(tempmode,",");
    while (chtok != NULL) {
        obselems++;
        if (strstr(chtok,"#") != NULL) {
            numpar++;
            printf("Found parameterized variable %d.\n",numpar);
        }
        chtok = strtok(NULL,",");
    }
    /* Define number of pointers to allocate memory for */
    if (numpar == 0) {
        n = 1;
    } else {
        n = numpar;
    }
    printf("NUMPAR=%d, N=%d\n", numpar, n);
    /* set up output arrays for all the parameterized variables
       found in input photmode string
     */
    status = AllocPhotPar(obs,n);

    if (status > PHOT_OK){
        printf("\n==>ERROR: Problems allocating memory for ObsmodeVars.\n");
        return(status);
    }
    strcpy(obs->photmode,photmode);

    if (numpar > 0){
        /* Start extracting values from photmode string */
        tok = 0;
        /* keep track of the number of terms in the obsmode string */
        strcpy(tempmode,photmode);
        /* get first token from the string */
        chtok = strtok(tempmode,",");
        while (chtok != NULL) {
            if (strstr(chtok,"#") != NULL) {
                printf("Adding parameter %s as parnames[%d]\n",chtok,tok);
                strcpy(obs->parnames[tok],chtok);
                tok++;
            }
            chtok = strtok(NULL,",");
        }

        /* Now, parse each of the stored parameter strings from the obsmode
           into parnames and  parvalues. */
        for (i=0;i<numpar;i++) {
            chtmp = (char *) calloc(strlen(obs->parnames[i]) + 1, sizeof(char));
            strcpy(chtmp, obs->parnames[i]);

            chtok = strtok(chtmp,"#");
            chval = strtok(NULL,"#");
            /* Replace parnames[i] with just the par name */
            strcpy(obs->parnames[i],chtok);
            /* save the value in parvalues array */
            obs->parvalues[i] = atof(chval);
            /* concatentate the '#'  back to the end of the parname */
            strcat(obs->parnames[i],"#");

            free(chtmp);
        }

        /* Now build up the new obsmode string for comparison with the
           reference IMPHTTAB table obsmode */
        /*
           Initialize intermediate variable
         */
        obsnames = (char **)calloc(obselems,sizeof(char *));
        for (i=0;i<obselems;i++)
            obsnames[i] = (char *) calloc(SZ_COLNAME,sizeof(char));
        /* Start adding elements to the intermediate variable */
        strcpy(tempmode,photmode);
        chtok=strtok(tempmode,",");
        tok = 0;
        obselems = 0;
        while (chtok != NULL) {
            loc = strchr(chtok,'#');
            if (loc == NULL) {
                strcpy(obsnames[tok],chtok);
            } else {
                indx = loc - chtok + 1;
                strncpy(obsnames[tok],chtok,indx);
            }
            obselems++;
            chtok = strtok(NULL,",");
            tok++;
        }
        /* build up new obsmode string */
        for (i=0;i<obselems;i++){
            strcat(obs->obsmode,obsnames[i]);
            if (i < obselems-1)
                strcat(obs->obsmode,",");
        }

        /* Free memory */
        for (i=0;i<obselems;i++)
            free(obsnames[i]);
        free(obsnames);

    } else {
        strcpy(obs->obsmode, photmode);
    }

    return(PHOT_OK);
}


/* Function to compare whether two obsmode strings are the same ignoring
 * case and order. Returns 0 if they are the same, non-zero otherwise.
 * MRD 4 Apr. 2011
 */
static int CompareObsModes(char *obsmode1, char *obsmode2) {

    char temp1[SZ_FITS_REC], temp2[SZ_FITS_REC];

    char * chtok;

    int i;

    /* regardless of order or case they should be the same length */
    if (strlen(obsmode1) != strlen(obsmode2))
        return 1;

    /* convert them both to all lowercase */
    strcpy(temp1, obsmode1);
    strcpy(temp2, obsmode2);

    for (i = 0; obsmode1[i]; i++) {
        temp1[i] = tolower(obsmode1[i]);
        temp2[i] = tolower(obsmode2[i]);
    }

    chtok = strtok(temp1, ",");
    while (chtok != NULL) {
        if (strstr(temp2, chtok) == NULL)
            return 1;

        chtok = strtok(NULL, ",");
    }

    return PHOT_OK;
}


int PrintMissingCols_IMPHTTAB (int missing, int numcols, int *nocol,
        char *cnames[], char *tabname, IRAFPointer tp) {

    /* Parameters:
       int missing       i: number of missing columns
       int numcols       i: number of columns expected in table
       int *nocol        i: array of YES/NO for each column, YES means missing
       char **cnames     i: array of colnames
       char *tabname     i: name of table columns were read from
       IRAFPointer tp    i: pointer to table, close it if necessary
    */

    extern int status;
    int j;
    /* If any columns are missing... */
    if (missing) {
        printf ("\n==>ERROR: %d columns not found in %s.\n", missing, tabname);
        /*trlerror (MsgText); */

        for (j=0; j< numcols; j++) {
            /* Recall which ones were marked missing... */
            if (nocol[j]) {
                /*... and print out that column's name */
                printf ("\n==>ERROR: Column %s not found in %s.\n", cnames[j], tabname);
                /*trlerror (MsgText); */
            }
        }
        c_tbtclo (tp);
        return (status = COLUMN_NOT_FOUND);
    }
    return(status);
}


/* This routine reads the relevant data from one row. */
static int ReadPhotTab (PhtCols *tabinfo, int row, PhtRow *tabrow) {

    extern int status;

    c_tbegtt (tabinfo->tp, tabinfo->cp_obsmode, row,
            tabrow->obsmode, SZ_FITS_REC-1);
    if (c_iraferr())
        return (status = TABLE_ERROR);

    return (status);
}


/* This routine reads the array data from one row.  The number of elements
   in the arrays is gotten, the array is allocated, and the wavelengths
   and throughputs are read into the arrays.
 */
static int ReadPhotArray (PhtCols *tabinfo, int row, PhtRow *tabrow) {

    int len_datacol, len_photvar;
    char col_nelem[SZ_COLNAME];
    char col_parval[SZ_COLNAME]="PAR";
    int nret;
    int n, col, i;
    
    n=0;
    col=0;
    i=0;

    extern int status;

    c_tbegtt (tabinfo->tp, tabinfo->cp_datacol, row,
            tabrow->datacol, SZ_COLNAME-1);
    if (c_iraferr())
        return (status = TABLE_ERROR);

    /* interpret the last character of 'datacol' as int for parameterized row */
    /* Find out how many elements there are in the throughput arrays,
       and allocate space for the arrays to be read from the table.
     */
    len_datacol = strlen(tabrow->datacol);
    len_photvar = strlen(tabinfo->photvar);
    if (len_datacol > len_photvar) {
        tabrow->parnum = (int)strtol(&tabrow->datacol[len_photvar],(char **)NULL, 10);
    } else {
        tabrow->parnum = 0;
    }

    tabrow->telem = 1;

    if (tabrow->parnum == 0){
        tabrow->parvals = (double **) malloc(sizeof(double *));
        tabrow->parnames = (char **) malloc(sizeof(char *));
        tabrow->nelem = (int *) malloc(sizeof(int));
        tabrow->results = (double *) malloc(sizeof(double));
        c_tbegtd (tabinfo->tp, tabinfo->cp_result[0], row, &tabrow->results[0]);
        if (c_iraferr())
            return (status = TABLE_ERROR);
    } else {
        /*
           We are working with a paramterized obsmode/row
         */
        /* Initialize NELEM arrays based on number of parameterized values
           that need to be interpreted for this obsmode */
        tabrow->nelem = (int *) malloc((tabrow->parnum+1) * sizeof(int));
        tabrow->parvals = (double **) malloc(tabrow->parnum * sizeof(double *));
        tabrow->parnames = (char **) malloc(tabrow->parnum * sizeof(char *));
        for (i=0,n=1; i<tabrow->parnum; i++,n++){
            tabrow->parnames[i] = (char *) malloc (SZ_COLNAME * sizeof(char));
            /*
               Copy name of parameterized value as defined in PRIMARY header
               into this row's structure for matching with obsmode later
               tabrow->parnames indexing starts at 0 to match obsmode parsing
               tabinfo->parnames indexing starts at 1 to match index in name
            */
            /*strcpy(tabrow->parnames[i],tabinfo->parnames[n]);*/
            c_tbegtt(tabinfo->tp, tabinfo->cp_parnames[n], row,
                     tabrow->parnames[i], SZ_COLNAME-1);
            if (c_iraferr())
                return (status = TABLE_ERROR);
        }
        if (tabrow->nelem == NULL || tabrow->parvals == NULL || tabrow->parnames == NULL)
            return (OUT_OF_MEMORY);

        for (n=1,col=0; n<=tabrow->parnum; n++,col++){
            /* Build up names of columns for parameterized values
               which need to be read in */
            sprintf(col_nelem,"NELEM%d",n);    /* col name: NELEMn */
            sprintf(col_parval,"PAR%dVALUES",n); /* col name: PARnVALUES */

            /* Read in number of elements for this parameterized variable */
            c_tbegti (tabinfo->tp, tabinfo->cp_nelem[n], row, &tabrow->nelem[n]);
            if (c_iraferr())
                return (status = TABLE_ERROR);
            /* Now that we have found how many elements are in this paramaters
               array, we can use that to allocate memory for the values */
            tabrow->parvals[col] = (double *) malloc (tabrow->nelem[n] * sizeof(double));

            /* Read in parameterized values for this variable */
            c_tbagtd (tabinfo->tp, tabinfo->cp_parvalues[n], row, tabrow->parvals[col], 1, tabrow->nelem[n]);
            if (c_iraferr())
                return (status = TABLE_ERROR);
            /* Keep track of the number of elements in the results column */
            tabrow->telem *= tabrow->nelem[n];

        } /* End of loop over parameterized variables columns: n[parnum] */


        /* Now that we have read in all the parameterized variable values,
           we know how many elements will be read in for the results */
        tabrow->results = (double *) malloc(tabrow->telem * sizeof(double));

        /* Start by reading in the full set of results from the
           specified results column */
        nret = c_tbagtd (tabinfo->tp, tabinfo->cp_result[tabrow->parnum], row, tabrow->results,1,tabrow->telem);
        if (c_iraferr() || nret < tabrow->telem)
            return (status = TABLE_ERROR);
    } /* End else statement */

    return (0);
}


/* Perform interpolation, if necessary, to derive the final
   output value for this obsmode from the appropriate row.
 */
static double ComputeValue(PhtRow *tabrow, PhotPar *obs) {
    /* Parameters:
       PhtRow  *tabrow:    values read from matching row of table
       char   *obsmode:    full obsmode string from header

       obsmode string needs to contain values of parameterized values
       appropriate for this observation.
     */

    double value;
    int n,p, nx, ndim;
    double *out;
    double *obsindx; /* index of each obsmode par value in tabrow->parvals */
    double *obsvals; /* value for each obsmode par in the same order as tabrow */

    int *ndpos;
    int **bounds; /* [ndim,2] array for bounds around obsvals values */
    int indx,pdim,ppos,xdim,xpos;
    int tabparlen;
   
    xdim=0;
    /*
       intermediate products used in iterating over dims
     */
    int iter, x;
    int dimpow,iterpow;
    double *ndposd;
    int b0,b1,pindx;
    int deltadim;            /* index of varying dimension */
    double bindx[2],bvals[2]; /* indices into results array for bounding values */
    double rinterp;          /* placeholder for interpolated result */
    BoundingPoint **points;   /* array of bounding points defining the area of interpolation */

    /* Define functions called in this functions */
    double linterp(double *, int, double *, double);
    void byteconvert(int, int *, int);
    int computedeltadim(BoundingPoint *, BoundingPoint *);
    long computeindex(int *, double *, int);
    void computebounds(double *, int, double, int *, int*);
    int strneq_ic(char *, char*, int);
    BoundingPoint **InitBoundingPointArray(int, int);
    void FreeBoundingPointArray(BoundingPoint **, int);

    /* Initialize variables and allocate memory */
    value = 0.0;
    ndim = tabrow->parnum;
    if (ndim == 0) {
        /* No parameterized values, so simply return value
           stored in 1-element array results */
        return(tabrow->results[0]);
    }
    dimpow = pow(2,ndim);

    obsindx = (double *)calloc(ndim, sizeof(double));
    obsvals = (double *)calloc(ndim, sizeof(double)); /* Final answer */
    ndpos = (int *)calloc(ndim, sizeof(int));
    ndposd = (double *)calloc(ndim, sizeof(double));
    bounds = (int **) calloc(ndim, sizeof(int *));
    for (p=0; p < ndim; p++)
        bounds[p] = (int *) calloc(2, sizeof(int));

    /* We have parameterized values, so linear interpolation
       will be needed in all dimensions to get correct value.
       Start by getting the floating-point indices for each
       parameterized value
     */
    /*
       Start by matching up the obsmode parameters with those in the table row
       These are the values along each dimension of the interpolation
     */
    for (p=0;p<ndim;p++){
        tabparlen = strlen(tabrow->parnames[p]);
        for(n=0;n<obs->npar;n++){
            if (strneq_ic(tabrow->parnames[p],obs->parnames[n],tabparlen)){
                obsvals[p] = obs->parvalues[n];
                break;
            }
        }

        if (obsvals[p] == 0.0) {
            printf("ERROR: No obsmode value found for %s\n",tabrow->parnames[p]);

            free(obsindx);
            free(obsvals);
            free(ndpos);
            free(ndposd);
            for (p=0; p < ndim; p++)
                free(bounds[p]);
            free(bounds);

            return ('\0');
        }

        /* Check whether we're going beyond the data in the table
           (extrapolation). If we are, return -9999 */
        nx = tabrow->nelem[p+1];
        if ((obsvals[p] < tabrow->parvals[p][0]) ||
            (obsvals[p] > tabrow->parvals[p][nx-1])) {
            if ((obs->extrap == False) ||
                    (strcasecmp(obs->parnames[p], "mjd#") != 0)) {
                printf("WARNING: Parameter value %s%g is outside table data bounds, returning -9999.\n", tabrow->parnames[p], obsvals[p]);

                free(obsindx);
                free(obsvals);
                free(ndpos);
                free(ndposd);
                for (p=0; p < ndim; p++)
                    free(bounds[p]);
                free(bounds);
                return -9999.0;
            } else {
                printf("WARNING: Parameter value %s%g is outside table data bounds and will be extrapolated.\n", tabrow->parnames[p], obsvals[p]);
            }
        }
    }

    /* Set up array of BoundingPoint objects to keep track of information needed
       for the interpolation
     */
    points = InitBoundingPointArray(dimpow, ndim);

    /* Now find the index of the obsmode parameterized value
       into each parameterized value array
       Equivalent to getting positions in each dimension (x,y,...).
     */
    for (p=0; p < ndim; p++){
        nx = tabrow->nelem[p+1];

        out = (double *) calloc(nx, sizeof(double));

        for (n=0; n < nx; n++)
            out[n] = n;

        /* Special index if extrapolating, otherwise interpolate as usual. */
        if (obsvals[p] < tabrow->parvals[p][0]) {
            value = 0;
        } else if (obsvals[p] > tabrow->parvals[p][nx-1]) {
            value = nx - 1;
        } else {
            value = linterp(tabrow->parvals[p], nx, out, obsvals[p]);
        }

        if (value == -99) {
            free(obsindx);
            free(obsvals);
            free(ndpos);
            free(ndposd);
            for (p=0; p < ndim; p++)
                free(bounds[p]);
            free(bounds);
            free(out);
            return('\0');
        }

        obsindx[p] = value;  /* Index into dimension p */
        computebounds(out, nx, (double)floor(value), &b0, &b1);

        bounds[p][0] = b0;
        bounds[p][1] = b1;
        /* Free memory so we can use this array for the next variable*/
        free(out);
    } /* End loop over each parameterized value */


    /* Loop over each axis and perform interpolation to find final result.
       For each axis, interpolate between all pairs of positions along the
       same axis.

       An example with 3 parameters/dimensions for a point with array index
       (2.2, 4.7, 1.3):

       Iteration 1: for all z and y positions, interpolate between pairs in x
           (2, 4, 1) vs (3, 4, 1),
           (2, 5, 1) vs (3, 5, 1),
           (2, 4, 2) vs (3, 4, 2),
           (2, 5, 2) vs (3, 5, 2)

       Iteration 2: for all z positions, interpolate between pairs from
                    iteration 1 in y
           (2.2, 4, 1) vs (2.2, 5, 1),
           (2.2, 4, 2) vs (2.2, 5, 2)

       Iteration 3: interpolate between pairs from iteration 2 in z
           (2.2, 4.7, 1) vs (2.2, 4.7, 2) ==> final answer
    */
    for (iter=ndim; iter > 0; iter--) {
        iterpow = pow(2, iter);
        for (p=0; p < iterpow; p++) {
            pdim = floor(p / 2);
            ppos = p % 2;

            if (iter == ndim) {
                /* Initialize all intermediate products and perform first
                   set of interpolations over the first dimension.

                   Create a bitmask for each dimension for each position. */
                byteconvert(p, ndpos, ndim);
                for (n=0; n < ndim; n++) {
                    pindx = bounds[n][ndpos[n]];
                    points[pdim][ppos].index[n] = (double)pindx;
                    points[pdim][ppos].pos[n] = tabrow->parvals[n][pindx];
                }

                /* Determine values from tables which correspond to
                   bounding positions to be interpolated */
                indx = computeindex(tabrow->nelem, points[pdim][ppos].index, ndim);
                points[pdim][ppos].value = tabrow->results[indx];
            } /* End if(iter==ndim) */

            if (ppos == 1) {
                /* Determine which axis is varying, so we know which
                   input value from the obsmode string
                   we need to use for the interpolation */
                deltadim = computedeltadim(&points[pdim][0],&points[pdim][1]);
                if (deltadim < 0 || deltadim >= ndim) {
                    printf("ERROR: Deltadim out of range: %i\n",deltadim);
                    free(obsindx);
                    free(obsvals);
                    free (ndpos);
                    free (ndposd);
                    for (p=0;p<ndim;p++)
                        free(bounds[p]);
                    free(bounds);

                    return('\0');
                }
                bindx[0] = points[pdim][0].pos[deltadim];
                bindx[1] = points[pdim][1].pos[deltadim];
                bvals[0] = points[pdim][0].value;
                bvals[1] = points[pdim][1].value;
                /* Perform interpolation now and record the results */
                rinterp = linterp(bindx, 2, bvals, obsvals[deltadim]);

                /* Update intermediate arrays with results in
                   preparation for the next iteration
                 */

                if (rinterp == -99)
                    return('\0');

                /* Determine where the result of this interpolation should go */
                x = floor((p-1)/2);
                xdim = floor(x/2);
                xpos = x%2;
                /* update bpos and bindx for iteration over next dimension
                 */
                points[xdim][xpos].value = rinterp;
                for (n=0;n<ndim;n++) {
                    points[xdim][xpos].index[n] = points[pdim][0].index[n];
                    points[xdim][xpos].pos[n] = points[pdim][0].pos[n];
                }
                points[xdim][xpos].index[deltadim] = obsindx[deltadim];
                points[xdim][xpos].pos[deltadim] = obsvals[deltadim];

            } /* Finished with this pair of positions (end if(ppos==1)) */

        } /* End loop over p, data stored for interpolation in changing dimension */

    } /* End loop over axes(iterations), iter, for interpolation */

    /* Record result */
    value = points[0][0].value;

    /* clean up memory allocated within this function */
    free(obsindx);
    free(obsvals);
    free (ndpos);
    free (ndposd);
    for (p=0;p<tabrow->parnum;p++)
        free(bounds[p]);
    free(bounds);

    FreeBoundingPointArray(points,dimpow);
    return (value);
}


/* This routine implements 1-D linear interpolation
   It returns the interpolated value from f(x) that corresponds
   to the input position xpos, where f(x) is sampled at positions x.
 */
double linterp(double *x, int nx, double *fx, double xpos) {

    int i0, i1;  /* x values that straddle xpos */

    double value;

    void computebounds (double *, int, double , int *, int *);

    /* interpolation calculated as:
       yi + (xpos - xi)*(yi1 - yi)/(xi1-xi)
     */
    /* Start by finding which elements in x array straddle xpos value */
    computebounds(x, nx, xpos, &i0, &i1);
    if ((x[i1] - x[i0]) == 0){
        printf("==>ERROR: Linear interpolation reached singularity...\n");
        return(-99);
    }
    /* Now, compute interpolated value. */
    /* This is just a straight line formula, can also extrapolate. */
    value = fx[i0] + (xpos - x[i0]) * (fx[i1] - fx[i0]) / (x[i1] - x[i0]);

    return(value);
}


/* Compute index into x array of val and returns
   the indices for the bounding values, taking into account
   boundary conditions.
 */
void computebounds (double *x, int nx, double val, int *i0, int *i1) {
    int n;

    /* first test for whether we've got an end case here */
    if (x[nx-1] <= val) {
        *i0 = nx - 2;
        *i1 = nx - 1;
    } else {
        for(n=0; n < nx; n++){
            if (x[n] <= val) {
                *i0 = n;
            } else {
                if (n > 0 && n < (nx - 1)){
                    *i1 = n;
                } else if (n == 0) {
                    *i0 = 0;
                    *i1 = 1;
                } else {
                    *i0 = nx - 2;
                    *i1 = nx - 1;
                }
                break;
            }
        }
    }
}


/* Compute the 1-D index of a n-D (ndim) position given by the array
   of values in pos[]
 */
long computeindex(int *nelem, double *pos, int ndim) {
    int n;
    int szaxis = 1;
    long indx = 0;

    for (n=0; n < ndim; n++) {
        indx += szaxis * pos[n];
        szaxis *= nelem[n+1];
    }
    return(indx);
}


/* Convert an int value into an array of 0,1 values to represent
   which bytes are 0 or 1 in the integer.  The result array must
   already be allocated for the number of bytes to be checked in the
   integer value.
 */
void byteconvert(int val, int *result, int ndim) {
    int i;
    int bval = 1;

    for (i=0; i < ndim; i++) {
        if ((val & bval) > 0){
            result[i] = 1;
        } else {
            result[i] = 0;
        }
        bval = bval << 1;
    }
}


/* Given 2 N dimensional sets of array indices, determine which dimension
   changes from one set to the other.

   NOTE:
   This assumes that the positions only change in 1 dimension at a time.
*/
int computedeltadim(BoundingPoint *pos1, BoundingPoint *pos2){
    int p;
    int xdim=0;
    double diff;

    for (p=0; p < pos1->ndim; p++) {
        diff = pos2->index[p] - pos1->index[p];
        if ( diff != 0) {
            xdim = p;
            break;
        }
    }
    return(xdim);
}


/* This routine frees memory allocated to a row's entries,
   so that the next call to 'ReadPhotTab' will have an empty
   structure to use for the storing the rows values.
 */
static void ClosePhotRow (PhtRow *tabrow) {

    int i;

    for (i=0; i<tabrow->parnum; i++){
        free(tabrow->parvals[i]);
        free(tabrow->parnames[i]);
    }
    free(tabrow->parvals);
    free(tabrow->parnames);

    free(tabrow->nelem);
    free(tabrow->results);
}


/* This routine closes the imphttab table. */
static int ClosePhotTab (PhtCols *tabinfo){

    extern int status;

    c_tbtclo (tabinfo->tp);
    if (c_iraferr())
        return (status = TABLE_ERROR);

    return (status);
}


/* Initialize the array of BoundingPoint objects */
BoundingPoint **InitBoundingPointArray(int npoints, int ndim){
    int i;
    int pdim;
    void InitBoundingPoint(BoundingPoint *, int);
    BoundingPoint **points;

    pdim = npoints/2;
    points = (BoundingPoint **)calloc(pdim,sizeof(BoundingPoint *));
    for (i=0;i<pdim;i++) {
        points[i] = (BoundingPoint *)calloc(2,sizeof(BoundingPoint));
        InitBoundingPoint(&points[i][0],ndim);
        InitBoundingPoint(&points[i][1],ndim);
    }
    return(points);
}


void InitBoundingPoint(BoundingPoint *point, int ndim){
    point->index = (double *)calloc(ndim, sizeof(double));
    point->pos = (double *)calloc(ndim, sizeof(double));
    point->ndim = ndim;
    point->value=0.0;
}


/* Free the memory allocated to an array of BoundingPoint objects */
void FreeBoundingPointArray(BoundingPoint **points, int npoints){
    int i;
    int pdim;
    void FreeBoundingPoint(BoundingPoint *);
    pdim = npoints/2;

    for (i=0;i<pdim;i++) {
        FreeBoundingPoint(&points[i][0]);
        FreeBoundingPoint(&points[i][1]);
        free(points[i]);
    }
    free(points);
}


void FreeBoundingPoint(BoundingPoint *point){
    free(point->index);
    free(point->pos);
}


/* This routine gets pedigree and descrip from the current table row.

   The pedigree and descrip columns need not be defined.  If not, this
   is flagged by their column pointers being zero, in which case no change
   will be made to the ref struct for this table.

   If the pedigree column is present, the value in the current row
   will replace any value previously gotten from the header.  As with
   the value from the header, if the first five letters of pedigree are
   DUMMY, then goodPedigree will be set to DUMMY_PEDIGREE; if not,
   goodPedigree will not be changed from its previous value (either as
   initialized or as gotten from the header).

   If the descrip column is present, the value is assigned to a second
   string, descrip2, rather than overwriting the one from the header.
   If the column is not found, descrip2 will be set to null.
 */
static int PhotRowPedigree (PhotPar *obs, int row,
        IRAFPointer tp, IRAFPointer cp_pedigree, IRAFPointer cp_descrip) {

    extern int status;

    /* Get pedigree and descrip.  If either or both are missing,
       that's not an error in this case.
     */
    if (cp_pedigree > 0) {
        c_tbegtt (tp, cp_pedigree, row, obs->pedigree, SZ_FITS_REC);
        if (c_iraferr())
            return (status = TABLE_ERROR);
        /* Is this row flagged as dummy? */
        if (strncmp (obs->pedigree, "DUMMY", 5) == 0)
            obs->goodPedigree = DUMMY_PEDIGREE;
        else
            obs->goodPedigree = GOOD_PEDIGREE;
    }

    if (cp_descrip > 0) {
        c_tbegtt (tp, cp_descrip, row, obs->descrip2, SZ_FITS_REC);
        if (c_iraferr())
            return (status = PHOTTABLE_ERROR);
    } else {
        obs->descrip2[0] = '\0';
    }


    return (status);
}


/* Initializes the PhotPar structure for use in this routine.

   This routine should be called by the user's code, and is not
   explicitly called within this library.
   The parameter's 'name' and 'pedigree' should come from RefTab.
*/
void InitPhotPar(PhotPar *obs, char *name, char *pedigree) {
    obs->name[0] = '\0';
    obs->pedigree[0] = '\0';

    /* Start by copying in required values from input table RefTab */
    strcpy(obs->name, name);
    strcpy(obs->pedigree, pedigree);

    /* Initialize remainder of fields to be used as output in this code */
    obs->descrip2[0] = '\0';
    obs->goodPedigree = PEDIGREE_UNKNOWN;

    obs->obsmode[0] = '\0';     /* obsmode of science data */
    obs->photmode[0] = '\0'; /* obsmode used for comparison with IMPHTTAB */

    /* parsed out value of any parameterized values */
    /* tab->obspars=NULL; */
    obs->npar = 0;

    /* Output values derived from table */
    obs->photflam = 0;
    obs->photplam = 0;
    obs->photbw = 0;
    obs->photzpt = 0;
    obs->phtflam1 = 0;
    obs->phtflam2 = 0;
}


int AllocPhotPar(PhotPar *obs, int npar){
    extern int status;
    int i;

    obs->npar = npar;

    obs->parnames = (char **)calloc(npar, sizeof(char *));
    for (i=0;i<npar;i++) {
        obs->parnames[i] = (char *)calloc(SZ_FITS_REC, sizeof(char));
        obs->parnames[i][0] = '\0';
    }
    obs->parvalues = (double *)calloc(npar, sizeof(double));
    if (obs->parnames == NULL || obs->parvalues == NULL) {
        return(status=OUT_OF_MEMORY);
    }

    printf("Allocated %d parnames\n", npar);
    return(status);
}


void FreePhotPar(PhotPar *obs){
    int n;

    for (n=0;n<obs->npar;n++){
        free(obs->parnames[n]);
    }
    free(obs->parnames);
    free(obs->parvalues);
}


/* This function compares two strings without regard to case, returning
   one if the strings are equal.

   Phil Hodge, 1997 Dec 12:
   Function created.
 */
int streq_ic_IMPHTTAB (char *s1, char *s2) {
    int c1, c2;
    int i;

    c1 = 1;
    for (i = 0;  c1 != 0;  i++) {

        c1 = s1[i];
        c2 = s2[i];
        if (isupper(c1))
            c1 = tolower (c1);
        if (isupper(c2))
            c2 = tolower (c2);
        if (c1 != c2)
            return (0);
    }
    return (1);
}


int strneq_ic (char *s1, char *s2, int n) {
    int c1, c2;
    int i;

    if (n == 0)
        return 0;

    c1 = 1;
    for (i = 0;  i < n;  i++) {

        c1 = s1[i];
        c2 = s2[i];
        if (isupper(c1))
            c1 = tolower (c1);
        if (isupper(c2))
            c2 = tolower (c2);
        if (c1 != c2)
            return (0);
    }
    return (1);
}
