# include <stdio.h>
# include <string.h>
# include <stdlib.h>    /* for strtol function */
# include <ctype.h>

# include "c_iraf.h"

# include "hstio.h"
# include "xtables.h"
# include "imphttab.h"


typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_obsmode;		/* column descriptors */
	IRAFPointer cp_datacol;	
	IRAFPointer cp_result[MAXPARS+1];
	IRAFPointer cp_nelem[MAXPARS];
	IRAFPointer cp_parvalues[MAXPARS];
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
    char photvar[SZ_COLNAME]; /* photometric parameter in this table */
	int nrows;			/* number of rows in table */
    int ncols;          /* number of columns in table */
    int parnum;         /* number of parameterized variables in tables */
    char **parnames; /* par. variable names from table hdr keywords*/
} PhtCols;

typedef struct {
	char obsmode[SZ_FITS_REC];  /* obsmode string read from table row */
	char datacol[SZ_FITS_REC];
    char **parnames; /* record the par. variable names for comparison with obsmode string */
    int parnum;      /* number of parameterized variables */
	double *results;  /* results[telem] or results[nelem1*nelem2*...] */
    int telem;     /* total number of parameterized variable values */
    int *nelem;    /* multiple paramerized variables will each N values */
	double **parvals; /* need to support multiple parameterized variables */
} PhtRow;


/* Internal functions to be used to interpret IMPHTTAB ref tables */
static int OpenPhotTab (char *, char *, PhtCols *);
static int ReadPhotTab (PhtCols *, int, PhtRow *);
static int ReadPhotArray (PhtCols *, int, PhtRow *);
static double ComputeValue (PhtRow *, PhotPar *);
static int ClosePhotTab (PhtCols *);
static int InterpretPhotmode(char *, PhotPar *);
static int PhotRowPedigree (PhotPar *, int , IRAFPointer, IRAFPointer, IRAFPointer);
static void ClosePhotRow (PhtRow *);

/* This routine gets the gain, bias and readnoise for the CCD from
   the CCD parameters table (keyword name CCDTAB).

   The CCD parameters table should contain the following:
	header parameters:
		none needed
	columns:
        OBSMODE:  observation mode for exposure
        DATACOL:  name of column which tneeds to be read to get result
        NELEM<n>: number of elements for parameter <n>
        PAR<n>VALUES: values for parameterized variable
        <result>: column with single result value 
        <result><n>: column with result based on parameterized variable <n>
        DESCRIP:  description of the source of the values
        PEDIGREE: description of the pedigree of the values
        
   The table is read to find the row for which the OBSMODE matches the
   generated PHOTMODE string (from the image header). There will be one
   table extension for each photometry <result>: 
       PHOTFLAM, PHOTPLAM, and PHOTBW.
*/

int GetPhotTab (PhotPar *obs, char *photmode) {

/* arguments:
RefTab  *obs     io: name and pedigree values
char    *photmode  i: photmode from science header for row selection

Calling routine will need to initialize and free PhotPar object using:

    RefTab *ref;
    PhotPar *obs;
    InitPhotPar(obs, ref->name, ref->pedigree);
    
    FreePhotPar(obs);

    where ref is the information about the reference file from the header
*/

	extern int status;

	PhtCols tabinfo;	/* pointer to table descriptor, etc */
	PhtRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int i;
    int extn;
    char phdrname[SZ_FNAME];
    IODescPtr im;		/* descriptor for primary header unit */
    Hdr tphdr;		/* primary header */
    FitsKw key;		/* location of keyword in header */
    char pname[SZ_COLNAME];
    int parnum,numpar;
    double value;
    int npar;
    int plen, n;
    int numkeys;
	
    int foundit;		/* has the correct row been found? */
    int PhotRowPedigree (PhotPar *, int, IRAFPointer, IRAFPointer, IRAFPointer);
    int SameInt (int, int);
    int SameFlt (float, float);
    int streq_ic (char *, char *);
    double ComputeValue(PhtRow *, PhotPar *);
    
    extern int status;
    extern char **photnames;

    
    /* Interpret OBSMODE string from science file header for
        comparison with obsmode values in reference table
        primarily by stripping out any parameterized values */
    /* Find the length of the string to be interpreted */
    status = InterpretPhotmode(photmode, obs);
    if (status != PHOT_OK)
        return(status=ERROR_RETURN);

    /* Create name of table to access PRIMARY header keywords */
    strcpy(phdrname, obs->name);
    /*strcat(phdrname,"[0]"); */
    
    /* Open PRIMARY header for reading keywords */
    initHdr (&tphdr);
    /* Open the primary header of the reference file. */
    im = openInputImage (phdrname, "", 0);
    if (hstio_err()) {
	    printf ("IMPHTTAB `%s' not found.", obs->name);
	    clear_hstioerr();
        status = OPEN_FAILED;
	    return (status);
    }
    getHeader (im, &tphdr);
    if (hstio_err()){
	    printf ("==>ERROR: IMPHTTAB `%s' not found.", obs->name);
	    return (status = HEADER_PROBLEM);
    }
    /* Read in keywords from PRIMARY header for use in interpreting
        the tables */
	key = findKw (&tphdr, "PARNUM");
	if (key == NotFound) {
		printf ("==>ERROR: Trying to get PARNUM...");
        closeImage(im);
		return (status = KEYWORD_MISSING);
	} else {
	    tabinfo.parnum = getIntKw(key);
	}
    
    tabinfo.parnames = (char **)calloc(tabinfo.parnum+1, sizeof(char *));
    for (i=0; i<=tabinfo.parnum; i++){
        tabinfo.parnames[i] = (char *) calloc (SZ_COLNAME, sizeof(char));
    }
	if (tabinfo.parnames == NULL)
	    return (OUT_OF_MEMORY);
    
    /* For each parameterized value that went into this obsmode,
        read in the column name for that value */
    for (parnum = 1; parnum <= tabinfo.parnum; parnum++){
        sprintf(pname,"PAR%dNAME",parnum);
    	key = findKw (&tphdr, pname);
	    getStringKw (key, tabinfo.parnames[parnum], SZ_COLNAME);
    }
    /* Read in PHOTZPT keyword value from Primary header */
	key = findKw (&tphdr, "PHOTZPT");
    obs->photzpt = getDoubleKw (key);
    if (hstio_err()){
        printf("==>ERROR: Keyword `PHOTZPT` not found in PRIMARY header.");
        closeImage(im);
        return (status = KEYWORD_MISSING);
    }
	key = findKw (&tphdr, "NEXTEND");
    numkeys = getIntKw (key);
    if (hstio_err()){
        printf("==>ERROR: Keyword `NEXTEND` not found in PRIMARY header.");
        closeImage(im);
        return (status = KEYWORD_MISSING);
    }
    
    /* Close PRIMARY header of reference table */
    closeImage(im);
	freeHdr (&tphdr);
    
    /* Now step through each of the extensions to compute the 
        difference photometry keyword values, one from each
        extension.
    */
    for (extn = 1; extn <= numkeys; extn++){ 
        strcpy(pname,&photnames[extn]);
	    /* Open the photometry parameters table and find columns. */
	    if (OpenPhotTab (obs->name, pname, &tabinfo))
	        return (status);

	    /* Check each row for a match with obsmode, 
            and get info from the matching row.
	    */

	    foundit = 0;
	    for (row = 1;  row <= tabinfo.nrows;  row++) {

	        /* Read the current row into tabrow. */
	        if (ReadPhotTab (&tabinfo, row, &tabrow))
		    return (status);

	        if (streq_ic (tabrow.obsmode, obs->obsmode)) {

		        foundit = 1;
		        if (PhotRowPedigree (obs, row,
			        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip))
		            return (status);
		        if (obs->goodPedigree == DUMMY_PEDIGREE) {
		            printf ("==>Warning: Row %d of IMPHTTAB is DUMMY.", row);
			        /* trlwarn (MsgText); */
		        }

                /* Read in photometry values from table row */
                if (status = ReadPhotArray(&tabinfo, row, &tabrow))
                    return (status);

                /* interpret values to compute return value */
                value = ComputeValue (&tabrow, obs);
                if (value == '\0'){
                    return(status=INTERNAL_ERROR);
                }
                printf("==> Value of %s = %g",photnames[extn], &value);

                /* Free memory used to read in this row */
                ClosePhotRow(&tabrow);

                break;
            } /* End of if SameString() */
	    } /* End of loop over table rows */

	    if (!foundit) {
	        printf ("==>ERROR: Matching row not found in IMPHTTAB `%s'.", obs->obsmode);
	        /*trlerror (MsgText); */
		    printf ("==>ERROR: OBSMODE %s",	obs->obsmode);
		    /*trlerror (MsgText); */

	        ClosePhotTab (&tabinfo);
	        return (status = TABLE_ERROR);
	    }
        
	    if (ClosePhotTab (&tabinfo))		/* close the table */
	        return (status);
            
        /* Store computed value for keyword as appropriate member in 
            PhotPar struct */
        if (strncmp(photnames[extn],"PHOTFLAM",6)) {
            obs->photflam = value;
        } else if (strncmp(photnames[extn],"PHOTPLAM",6)) {
            obs->photplam = value;
        } else {
            obs->photbw = value;
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

    int colnum, datatype, lendata, lenfmt;
    char tname[SZ_FNAME];    
    char **colnames, **ecolnames, **pcolnames;

	int *nocol;
	int i, j, missing;
    int parnum;
    	
	int PrintMissingCols (int, int, int *, char **, char *, IRAFPointer);
    int buildTabName (char *, char *, char *);

    /* allocate space for column names to be read from this table */
    colnames = (char **)calloc(tabinfo->parnum+1, sizeof(char *));
    ecolnames = (char **)calloc(tabinfo->parnum+1, sizeof(char *));
    pcolnames = (char **)calloc(tabinfo->parnum+1, sizeof(char *));
    for (i=0;i <= tabinfo->parnum; i++){
        colnames[i] = (char *)calloc(SZ_COLNAME, sizeof(char));
        ecolnames[i] = (char *)calloc(SZ_COLNAME, sizeof(char));
        pcolnames[i] = (char *)calloc(SZ_COLNAME, sizeof(char));
    }
    /* Copy in the column names now */
    strcpy(colnames[0],photvar);
    for (parnum=1;parnum <= tabinfo->parnum; parnum++){
        sprintf(colnames[parnum],"%s%d",photvar,parnum);
        sprintf(ecolnames[parnum],"NELEM%d",parnum);
        sprintf(pcolnames[parnum],"PAR%dVALUES",parnum);
    }    
    
    /* Create name of table with extension to be opened */
    sprintf(tname,"%s[%s]",tabname,photvar);
    
    /* keep track of what extension we are processing here */
    strcpy(tabinfo->photvar, photvar);

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("==>ERROR: IMPHTTAB extension`%s' not found.", tname);
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
    for (parnum = 1;parnum <= tabinfo->parnum; parnum++){
        /* 
            Read in NELEM<parnum> columns 
        */
	    c_tbcfnd1 (tabinfo->tp, ecolnames[parnum], &tabinfo->cp_nelem[parnum]);
	    /* 
            Read in PAR<parnum>VALUES columns 
        */
        c_tbcfnd1 (tabinfo->tp, pcolnames[parnum], &tabinfo->cp_parvalues[parnum]);
        /* 
            Read in the PHOT*<parnum> columns 
        */
	    c_tbcfnd1 (tabinfo->tp, colnames[parnum], &tabinfo->cp_result[parnum]);
    }	
	/* Initialize counters here... */
	missing = 0;
	i=0;
		
    /* Increment i for every column, mark only missing columns in
        nocol as YES.  WJH 27 July 1999
    */
	if (tabinfo->cp_obsmode == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_datacol == 0 ) { missing++; nocol[i] = YES;} i++;
	if (tabinfo->cp_result[0] == 0 ) { missing++; nocol[i] = YES;} i++;
	
	if (PrintMissingCols (missing, tabinfo->ncols, nocol, colnames, "IMPHTTAB", tabinfo->tp) )
		return(status);
		
	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);


	return (status);
}

static int InterpretPhotmode(char *photmode, PhotPar *obs){
/* Interprets PHOTMODE string (as a comma-separated, all-lowercase string)
    from the science image header into obsmode for comparison with 
    obsmode column values in IMPHTTAB while also recording the parameterized
    variable values from the PHOTMODE string as outputs parnames, parvalues.
*/
    
    int i,j;
    int numpar, n;
    int tok;
    char *chtok,*chval;
    char **obsnames;
    int obselems;

	extern int status;
    
    /* Temporary string for parameterized variable value from photmode */
    char tempmode[SZ_FITS_REC];
    
    int AllocPhotPar(PhotPar *, int);

    numpar = 0;
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
    printf("NUMPAR = %d, N=%d\n",numpar,n);
    /* set up output arrays for all the parameterized variables 
        found in input photmode string
    */
    status = AllocPhotPar(obs,n);
    
    if (status > PHOT_OK){
        printf("==>ERROR: Problems allocating memory for ObsmodeVars.");
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
            chtok = strtok(obs->parnames[i],"#");
            chval = strtok(NULL,"#");
            /* Replace parnames[i] with just the par name */
            strcpy(obs->parnames[i],chtok);
            /* save the value in parvalues array */
            obs->parvalues[i] = atof(chval);
        }
            
        /* Now build up the new obsmode string for comparison with the 
        reference IMPHTTAB table obsmode */
        /* 
            Initialize intermediate variable 
        */
        obsnames = (char **)calloc(obselems,sizeof(char *));
        for (i=0;i<obselems;i++) obsnames[i]=(char *) calloc(SZ_COLNAME,sizeof(char));
        /* Start adding elements to the intermediate variable */
        strcpy(tempmode,photmode);
        chtok=strtok(tempmode,",");
        tok = 0;
        obselems = 0;
        while (chtok != NULL) {
            if (strstr(chtok,"#") == NULL) {
                strcpy(obsnames[tok],chtok);
                obselems++;
            }
            chtok = strtok(NULL,",");
            tok++;
        }
        /* build up now obsmode string */
        for (i=0;i<obselems;i++){
            strcat(obs->obsmode,obsnames[i]);
            if (i < obselems-1)strcat(obs->obsmode,",");
        }
    } else {
        strcpy(obs->obsmode, photmode);
    }   
    /* Free memory */ 
    for (i=0;i<obselems;i++) free(obsnames[i]);
    free(obsnames);
    
    return(PHOT_OK);
}

int PrintMissingCols (int missing, int numcols, int *nocol, char *cnames[], char *tabname, IRAFPointer tp) {

/* Parameters:
int missing			i: number of missing columns
int numcols			i: number of columns expected in table
int *nocol			i: array of YES/NO for each column, YES means missing
char **cnames		i: array of colnames
char *tabname		i: name of table columns were read from
IRAFPointer tp		i: pointer to table, close it if necessary
*/

	extern int status;
	int j;
	/* If any columns are missing... */
	if (missing) {
 	    printf ("==>ERROR: %d columns not found in %s.", missing, tabname);
		/*trlerror (MsgText); */
       
		for (j=0; j< numcols; j++) {
			/* Recall which ones were marked missing... */
			if (nocol[j]) {
				/*... and print out that column's name */
	    		printf ("==>ERROR: Column %s not found in %s.", cnames[j], tabname);
				/*trlerror (MsgText); */
			}
		}
	    c_tbtclo (tp);
	    return (status = COLUMN_NOT_FOUND);
	}
	return(status);
}
/* This routine reads the relevant data from one row.  
*/

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

	int nwl, nthru;		/* number of elements actually read */
    int len_datacol, len_photvar;
    char col_nelem[SZ_COLNAME];
    char col_parval[SZ_COLNAME]="PAR";
    int nret;
    int parnum;
    int n, col, i;

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
        tabrow->parvals = (double **)calloc(1,sizeof(double *));
        tabrow->parvals[0] = (double *) calloc(1,sizeof(double));
        tabrow->nelem = calloc(1,sizeof(float));
        tabrow->results = calloc(1,sizeof(double));
	    c_tbegtd (tabinfo->tp, tabinfo->cp_result[0], row, &tabrow->results[0]);
	    if (c_iraferr())
	        return (status = TABLE_ERROR);
    } else {
        /* 
            We are working with a paramterized obsmode/row 
        */
        /* Initialize NELEM arrays based on number of parameterized values 
            that need to be interpreted for this obsmode */
        tabrow->nelem = calloc(tabrow->parnum+1, sizeof(int));
        tabrow->parvals = (double **)calloc(tabrow->parnum, sizeof(double *));
        tabrow->parnames = (char **)calloc(tabrow->parnum, sizeof(char *));
        for (i=0,n=1; i<tabrow->parnum; i++,n++){
            tabrow->parnames[i] = (char *) calloc (SZ_COLNAME, sizeof(char));
            /* 
                Copy name of parameterized value as defined in PRIMARY header 
                into this row's structure for matching with obsmode later
                tabrow->parnames indexing starts at 0 to match obsmode parsing
                tabinfo->parnames indexing starts at 1 to match index in name 
            */
            strcpy(tabrow->parnames[i],tabinfo->parnames[n]);
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
            tabrow->parvals[col] = (double *) calloc (tabrow->nelem[n], sizeof(double));

            /* Read in parameterized values for this variable */
	        c_tbagtd (tabinfo->tp, tabinfo->cp_parvalues[n], row, tabrow->parvals[col], 1, tabrow->nelem[n]);
	        if (c_iraferr())
	            return (status = TABLE_ERROR);        
            /* Keep track of the number of elements in the results column */
            tabrow->telem *= tabrow->nelem[n];
            
        } /* End of loop over parameterized variables columns: n[parnum] */
            
        
        /* Now that we have read in all the parameterized variable values,
            we know how many elements will be read in for the results */
        tabrow->results = calloc(tabrow->telem,sizeof(double));

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
    int parnum;
    int n,p, nx, ndim;
    double *out;
    double *obsindx; /* index of each obsmode par value in tabrow->parvals */ 
    double *obsvals; /* value for each obsmode par in the same order as tabrow */

    int *ndpos;
    int **bounds; /* [ndim,2] array for bounds around obsvals values */
    double resvals[2]; /* Values from tabrow->results for interpolation */
    double resindx[2]; /* 1-D indices into tabrow->results for bounding positions*/
    int posindx;      /* N dimensional array of indices for a position */
    int indx,pdim,ppos;
    /* 
       intermediate products used in iterating over dims 
    */
    int iter, x;
    int dimpow,iterpow;
    double *ndposd;
    int b0,b1;
    int deltadim;            /* index of varying dimension */
    double **bvals, ***bpos; /* values and indices of results into arrays */ 
    double **bindx;             /* indices into results array for bounding values */
    double rinterp;          /* placeholder for interpolated result */
        
    /* Define functions called in this functions */ 
    double linterp(double *, int, double *, double);
    void byteconvert(int, int *, int);
    int computedeltadim(double *, double *);  
    long computeindex(int *, int *, int);  
    void computebounds(double *, int, double, int *, int*);  
    int streq_ic(char *, char*);

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
    for (p=0;p<ndim;p++) bounds[p] = (int *) calloc(2, sizeof(int));
    
    bpos = (double ***)calloc((dimpow)/2, sizeof(double **));
    bvals = (double **)calloc((dimpow)/2, sizeof(double *));
    bindx = (double **)calloc((dimpow)/2, sizeof(int *));
    for (p=0;p<ndim;p++){
        bvals[p] = (double *) calloc(2, sizeof(double));
        bindx[p] = (double *) calloc(2, sizeof(double));
        bpos[p] = (double **) calloc(2, sizeof(double *));
        for (n=0;n<2;n++) bpos[p][n] = (double *) calloc(ndim, sizeof(double));
    }
    
    
    /* We have parameterized values, so linear interpolation 
        will be needed in all dimensions to get correct value.
        Start by getting the floating-point indices for each 
        parameterized value
    */
    /* 
        Start by matching up the obsmode parameters with those in the table row 
    */
    for (p=0;p<ndim;p++){
        for(n=0;n<obs->npar;n++){
            if (streq_ic(tabrow->parnames[p],obs->parnames[n])){
                obsvals[p] = obs->parvalues[n];
                break;
            }
        }
        if (obsvals[p] == 0.0) {
            printf("ERROR: No obsmode value found for %s\n",tabrow->parnames[p]);
            free(obsindx);
            free(obsvals);
            free (ndpos);
            free (ndposd);
            for (p=0;p<ndim;p++)free(bounds[p]);
            free(bounds);
            for (p=0;p<ndim;p++){
                for (n=0;n<2;n++) free(bpos[p][n]);
                free(bvals[p]);
                free(bindx[p]);
                free(bpos[p]);
            }
            free(bpos);
            free(bvals);
            free(bindx);
            
            return('\0');
        }
    }
    /* Now find the index of the obsmode parameterized value
        into each parameterized value array
        Equivalent to getting positions in each dimension (x,y,...).
    */
    for (p=0;p<ndim;p++){    
        nx = tabrow->nelem[p+1];
        out = (double *)calloc(nx, sizeof(double));
        for (n=0; n<nx;n++) out[n] = n;
        value = linterp(tabrow->parvals[p], nx, out, obsvals[p]);
        if (value == '\0') {
            return('\0');
        }
        obsindx[p] = value;
        /* remember the bounding values for this position */
        computebounds(out, nx, (double)floor(obsindx[p]), &b0, &b1);
        bounds[p][0] = tabrow->parvals[p][b0];
        bounds[p][1] = tabrow->parvals[p][b1];
        
        /* Free memory so we can use this array for the next variable*/
        free(out);
    } /* End loop over each parameterized value */
    /* 
        Loop over each axis and perform interpolation to find final result 
    */
    for (iter=ndim; iter >0; iter--) {
        iterpow = pow(2,iter);
        for (p=0;p < iterpow;p++){    
            pdim = floor(p/2);
            ppos = p%2;

            if (iter == ndim) {
                /* Initialize all intermediate products and perform first 
                    set of interpolations over the first dimension 
                */
                /* Create a bitmask for each dimension for each position */
                byteconvert(p,ndpos,ndim);
                for (n=0;n<ndim;n++) bpos[pdim][ppos][n] = (double)ndpos[n];

                /* Determine values from tables which correspond to 
                    bounding positions to be interpolated */
                indx = computeindex(tabrow->nelem, ndpos, ndim);
                bvals[pdim][ppos] = tabrow->results[indx];

            } /* End if(iter==ndim) */ 
            if (ppos == 1) {
                /* Determine which axis is varying, so we know which 
                    input value from the obsmode string 
                    we need to use for the interpolation */
                deltadim = computedeltadim(bpos[pdim][1],bpos[pdim][0]);
                bindx[pdim][0] = bounds[deltadim][0];
                bindx[pdim][1] = bounds[deltadim][1];      
                /*Perform interpolation now and record the results */
                rinterp = linterp(bindx[pdim], 2, bvals[pdim],obsindx[deltadim]);

                /* 
                    Update intermediate arrays with results in 
                    preparation for the next iteration 
                */
                if (rinterp == '\0') return(rinterp);
                /* Determine where the result of this interpolation should go */
                x = (p-1)/2;
                bvals[x][x%2] = rinterp;
                /* update bpos and bindx for iteration over next dimension */
                for (n=0;n<ndim;n++) bpos[x][x%2][n] = bpos[pdim][ppos][n];
                bpos[x][x%2][deltadim] = obsvals[deltadim];
                
            } /* Finished working out what dimension is being interpolated */
            
        } /* End loop over p, data stored for interpolation in changing dimension */
        
    } /* End loop over axes(iterations), iter, for interpolation */

    /* Record result */
    value = bvals[0][0];
    
    /* clean up memory allocated within this function */
    free(obsindx);
    free(obsvals);
    free (ndpos);
    free (ndposd);
    for (p=0;p<tabrow->parnum;p++)free(bounds[p]);
    free(bounds);

    for (p=0;p<ndim;p++){
        for (n=0;n<2;n++) free(bpos[p][n]);
        free(bvals[p]);
        free(bindx[p]);
        free(bpos[p]);
    }
    free(bpos);
    free(bvals);
    free(bindx);

    return (value);
}


/* This routine implements 1-D linear interpolation
    It returns the interpolated value from f(x) that corresponds
    to the input position xpos, where f(x) is sampled at positions x.
*/
double linterp(double *x, int nx, double *fx, double xpos) {

    int i0, i1;  /* x values that straddle xpos */

    int n;
    double value;

    void computebounds (double *, int, double , int *, int *);
    
    /* interpolation calculated as: 
        yi + (xpos - xi)*(yi1 - yi)/(xi1-xi)
    */
    /* Start by finding which elements in x array straddle xpos value */
    computebounds(x, nx, xpos, &i0, &i1);
    
    if ((x[i1] - x[i0]) == 0){
        printf("==>ERROR: Linear interpolation reached singularity...\n");
        return('\0');
    }
    /* Now, compute interpolated value */
    value = fx[i0] + (xpos - x[i0])*(fx[i1] - fx[i0])/(x[i1]-x[i0]);

    return(value);
} 

/* Compute index into x array of val and returns 
    the indices for the bounding values, taking into account 
    boundary conditions. 
*/
void computebounds (double *x, int nx, double val, int *i0, int *i1) {
    int n;
    
    
    for(n=0;n < nx; n++){
        if (x[n] <= val) {
            *i0 = n;
        } else {
            if (n > 0 && n < nx -1 ){
                *i1 = n;
            } else if (n == 0) {
                *i0 = 0;
                *i1 = 1;
            } else {
                *i0 = nx-2;
                *i1 = nx-1;
            }
            break;
        }
    }
}

/* Compute the 1-D index of a n-D (ndim) position given by the array 
    of values in pos[] 
*/
long computeindex(int *nelem, int *pos, int ndim) {
    int n, szaxis;    
    long indx;
    
    indx = 0;
    szaxis = 1;
    for (n=0;n<ndim;n++) {
        indx += szaxis*pos[n];
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
    int i,bval;
    
    bval = 1;
    for (i=0;i<ndim;i++){
        if ((val & bval) > 0){
            result[i] = 1;
        } else {
            result[i] = 0;
        }
        bval = bval << 1;
    }
}

/*
    Given 2 N dimensional sets of array indices, determine which dimension 
    changes from one set to the other.  
    
    NOTE: 
        This assumes that the positions only change in 1 dimension at a time. 
*/
int computedeltadim(double *ndpos1, double *ndpos2){
    int p,ndim;
    int xdim;
    
    ndim = sizeof(ndpos1)/sizeof(*ndpos1);
    for (p=0;p<ndim;p++){
        if (abs(ndpos2[p] - ndpos1[p]) > 0) {
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

void InitPhotPar(PhotPar *obs, char *name, char *pedigree) {
/* Initializes the PhotPar structure for use in this routine.

    This routine should be called by the user's code, and is not 
    explicitly called within this library.
    The parameter's 'name' and 'pedigree' should come from RefTab.
*/
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
    
}

int AllocPhotPar(PhotPar *obs, int npar){
    extern int status;

    int i;

    obs->npar = npar;
    
    obs->parnames = (char **)calloc(npar, sizeof(char *));
    printf("Allocated %d parnames\n",npar);
    for (i=0;i<npar;i++) {
        obs->parnames[i] = (char *)calloc(SZ_FITS_REC, sizeof(char));
        obs->parnames[i][0] = '\0';
    }
    obs->parvalues = (double *)calloc(npar, sizeof(double));
    if (obs->parnames == NULL || obs->parvalues == NULL) {
        return(status=OUT_OF_MEMORY);
    }
       
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

int streq_ic (char *s1, char *s2) {

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
