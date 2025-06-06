/* These routines provide the interface to the WF3 association tables.
 **	
 **	Howard Bushouse, 2000 Aug 23
 **		Initial version.
 **		[Adapted from acstable.c by W. Hack]
 **	H.Bushouse, 2001 May 8
 **		Modified in sync with changes to acstable:
 **		Corrected bugs dealing with dither associations that have
 **		single images at each pointing.
 **	H.Bushouse, 2002 June 20
 **		Modified to keep in sync with CALACS acstable: Always set
 **		dthcorr to PERFORM.
 **	H.Bushouse, 2002 Nov 26
 **		Fixed bug with interpreting RPn MEMTYPEs. RPTLEN reset to 2,
 **		and searching for '-rp' instead of '-rpt' (following CALACS
 **		changes).
 **	H.Bushouse, 2006 June 20
 **		Modified to keep in sync with CALACS acstable: Revised
 **		'getAsnTable' to only populate sub-products if at least
 **		one input for that product is present.
 **	H.Bushouse, 2008 Aug 20
 **		Changed all occurences of dither product file name suffix
 **		from "dth" to "drz".
 **	H.Bushouse, 2010 Oct 20
 **		Modified GetAsnTable to turn off CRCORR/RPTCORR processing if
 **		there aren't any sub-products with > 1 member. (PR 66366)
 **  M. Sosey 2015
 **      Updates for CTE calibration code, see #1193 
 */

# include <ctype.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#include "hstcal.h"
# include "xtables.h"	/* defines TABLE I/O functions */
# include "hstio.h"	/* defines HST I/O functions */

# include "wf3.h"	/* defines WF3 data structures */
# include "wf3asn.h"	/* defines WF3 Association data structures */
# include "hstcalerr.h"
# include "calwf3.h"
# include "trlbuf.h"

# define CRLEN	2
# define RPTLEN 2
# define DTHLEN 3
# define NSUF	14

# define MEMABSENT -1
# define DASH_CHAR '-'
# define UNDERLINE_CHAR '_'

typedef struct {
    /* char memname[CHAR_FNAME_LENGTH+1];
       char mtype[CHAR_FNAME_LENGTH+1];
       Bool prsnt;
       int posid;
       char type[CHAR_FNAME_LENGTH+1]; */
    char memname[20+1];
    char mtype[10+1];
    Bool prsnt;
    int posid;
    char type[10+1];
} RowInfo;

static int IsProduct (char *);
static int UpdateHdr (char *);
int checkGlobalInfo (AsnInfo *);

# define NCOLS	3	/* number of columns in ASNTAB */

/* LoadAsn: Parse the input string, determine if it is an assoc. table
 **	name or a single image.  Then call GetAsnTable accordingly.
 */

int LoadAsn (AsnInfo *asn) {

    /* Arguments:
     **  input		 i: Name of input file or image
     **  asn			io: Association info structure
     */

    extern int status;
    void printInfo (AsnInfo *);
    int SetInput (AsnInfo *);
    int SetAsnSingle (AsnInfo *);
    int GetAsnTable (AsnInfo *);
    int GetGlobalInfo (AsnInfo *);

    /* Determine whether input is a single file, an association table,
     ** or an entry from an association table. */
    if (SetInput (asn))
        return (status);

    if (asn->process == FULL) {
        trlmessage("LoadAsn:  Processing FULL Association");
    } else if (asn->process == PARTIAL) {
        trlmessage("LoadAsn:  Processing PART of Association");
    } else {
        trlmessage("LoadAsn:  Processing SINGLE exposure");
    }

    /* Read in global info from ASN table's primary header */	
    if (GetGlobalInfo (asn)) {
        trlerror(" Problem getting primary header information.");
        return (status);
    }


    /* Read in ASN table, and load appropriate members info into memory */
    if (asn->process == SINGLE) {
        /* Set ASN structure values to process a single exposure */
        if (SetAsnSingle (asn))
            return (status);
    } else {
        if (GetAsnTable (asn))
            return (status);
    }

    if (asn->debug) { 
        trlmessage("LoadAsn:  Read in ASN table %s ", asn->asn_table);
    }


    /* Print a summary of information about the association */
    if (asn->verbose)
        printInfo (asn);

    return (status);
}


/* SETASNSINGLE: Set the values in the structure ASN for single
 **	exposure processing.
 */

int SetAsnSingle (AsnInfo *asn){

    extern int status;
    int numsp = 1;
    int i;

    /* Allocate the member structures */
    asn->spmems = (int *)calloc(numsp+1,sizeof(int));		
    for (i=0; i <=numsp; i++)
        asn->spmems[i] = 1;

    asn->process = SINGLE;
    asn->crcorr  = OMIT;
    asn->rptcorr = OMIT;
    asn->dthcorr = OMIT;
    asn->numprod = 1;
    asn->numsp   = 1;		
    asn->numasn  = 1;
    asn->product = NULL;		

    return(status);
}

/* SETINPUT: Determine whether input is a single image, an association
 **	table, or an product from an association table.  It will then make sure
 **	the specified input exists.
 */

int SetInput (AsnInfo *asn) {

    /* Arguments:
     **  asn			io: Association info structure
     */	
    extern int status;

    /* Local Variables */
    char filename[CHAR_FNAME_LENGTH+1];	
    int exist;			/* EXISTS_YES or EXISTS_NO */
    int in_dot;
    char linput[CHAR_FNAME_LENGTH+1];	/* Lower case version of input */
    int incase;			/* What kind of input do we have? */

    int  DoesFileExist (char *);
    char *lowcase (char *, char *);
    int  GetAsnName (char *, char *);
    void FindAsnRoot (const char *, char *);

    /* Initialize internal variables here... */
    filename[0] = '\0';

    /* convert input to lowercase for use only in looking for extensions */
    lowcase (linput, asn->input);

    /* Determine what kind of input: root only, 
     **				 root + suffix, or 
     **				 root + suffix + extension
     ** We will look for '.fit' or '.fits' extensions.
     */

    /* First, let's complete the input filename with a FITS extension.*/
    if (strstr(linput, ".fit") == NULL) {
        strcpy (filename, asn->input);
        strcat (filename, ".fits");
    } else {
        strcpy (filename, asn->input);
        strcat (filename, "\0");
    } 		

    /* Initialize local variable */
    incase = 0;

    /* Now, we find out what kind of input name was provided... */
    if (strstr(linput, "_asn") != NULL){
        incase = 0;
    } else if (strstr(linput, "_raw") != NULL ) {
        incase = 1;
    } else if (!IsProduct (linput)) {
        /* Not a sub-product/intermediate product... */
        incase = 2;
    } else {
        /* treat as if it was only a rootname (plus '.fits') */
        incase = 3;
    }

    if (asn->debug) {
        trlmessage("GetAsnTable: incase = %d",incase);
    }

    /* Given this full filename for a file which exists, 
     ** copy out just the rootname. */
    FindAsnRoot(filename, asn->rootname);

    /* Can we find the file as it is input? */
    exist = DoesFileExist(filename);

    /* Treat filename according to what was input... */
    switch (incase) {

        case 0:  /* We have an ASN file explicitly specified */

            if (exist == EXISTS_YES) {

                /* Found an ASN file:  Set values in ASN structure */
                asn->process = FULL;
                strcpy(asn->asn_table, filename);
                strcpy(asn->filename, filename);

                if (asn->verbose) {
                    trlmessage("Processing FULL ASN table...");
                }			

            } else {

                /* Couldn't find specified ASN file... */
                trlerror("File %s not found for processing", filename);
                return (status = OPEN_FAILED);
            }
            break;

        case 1:  /* We have a RAW file explicitly specified */

            if (exist == EXISTS_YES) {

                /* Found a RAW file: Set values in ASN structure */
                strcpy (asn->filename, filename);
                strcpy (asn->asn_table, filename);
                strcat (asn->asn_table, "\0");

                if (asn->verbose) {
                    trlmessage("Processing SINGLE image %s...", asn->filename);
                }
                asn->process = SINGLE;

            } else {

                /* Couldn't find specified RAW file... */
                trlerror("File %s not found for processing", filename);
                return (status = OPEN_FAILED);
            }
            break;

        case 2:  /* We have a sub-product/intermediate file specified
                  ** for re-processing. They should really just run the
                  ** stand-alone tasks separately by hand, but...  */

            if (exist == EXISTS_YES) {
                strcpy (asn->filename, filename);

                /* Look for ASN_TAB in file's header, and copy to
                 ** ASN->asn_table if it is found. */
                if (!GetAsnName (filename, asn->asn_table) ) {

                    /* No ASN table listed, process as a
                     ** SINGLE exposure */
                    asn->process = SINGLE;				
                    if (asn->verbose) {
                        trlmessage("Re-Processing a SINGLE image from ASN table...");
                    }
                } else {

                    /* ASN table given in file, PARTIALLY process
                     ** table */ 	
                    asn->process = PARTIAL;	
                    if (asn->verbose) {
                        trlmessage("Re-Processing PART of ASN table...");
                    }
                }

            } else {

                /* Couldn't find specified file... */
                trlerror("File %s not found for re-processing", filename);
                return (status = OPEN_FAILED);
            }
            break;

        case 3:  /* We only have a rootname (with .fits extension)
                  ** specified */
        default:		

            if (exist == EXISTS_YES) {
                strcpy (asn->filename, filename);

                /* Look for ASN_TAB in file's header, and copy to
                 ** ASN->asn_table if it is found. */
                if (!GetAsnName (filename, asn->asn_table)) {

                    /* No ASN table listed, process as a SINGLE
                     ** exposure */
                    asn->process = SINGLE;				
                    if (asn->verbose) {
                        trlmessage("Processing a SINGLE image from ASN table...");
                    }
                } else {

                    /* ASN table given in file, PARTIALLY process
                     ** table */ 	
                    asn->process = PARTIAL;	
                    if (asn->verbose) {
                        trlmessage("Processing PART of ASN table...");
                    }
                }			

            } else {

                /* Let's start fresh, shall we?  At this point we only
                 ** have as rootname and (maybe) a .fits extension... 

                 in_dot = strcspn (filename, ".fit");
                 */	
                in_dot = strlen(filename) -
                    strlen(strstr (filename, ".fit"));

                if (asn->debug){
                    trlmessage("For file %s, in_dot = %d", filename, in_dot);
                }

                /* Truncate extension off of filename */
                filename[in_dot] = '\0';

                /* Create full ASN table name from input rootname */
                strcat (filename, "_asn.fits");
                exist = DoesFileExist (filename);

                if (exist == EXISTS_YES) {
                    strcat (filename, "\0");
                    /* If a name was specified, the file must exist. */
                    strcpy (asn->asn_table, filename);
                    if (asn->verbose) {
                        trlmessage("Found and processing FULL ASN table...");
                    }
                    asn->process = FULL;
                } else {

                    /* Couldn't find ASN file, so
                     ** trim off suffix previously tested. */
                    in_dot = strlen (filename) - 
                        strlen (strstr(filename, "_asn")) ;
                    filename[in_dot] = '\0';

                    /* ... and look for RAW file */
                    strcat (filename, "_raw.fits");
                    exist = DoesFileExist (filename);

                    if (exist == EXISTS_YES) {

                        /* Found a RAW file;
                         ** Set values in ASN structure */
                        strcat (filename, "\0");
                        strcpy (asn->asn_table, filename);
                        if (asn->verbose) {
                            trlmessage("Processing a SINGLE image from ASN table...");
                        }
                        asn->process = SINGLE;

                    } else {

                        /* Is INPUT the rootname of a product which
                         ** hasn't been created yet? */
                        /* Extract just the rootname plus last character
                         ** of rootname (which may be the subproduct ID*/
                        in_dot = strlen (filename) - 
                            strlen (strstr (filename, "_raw"));
                        filename[in_dot -1] = '\0';

                        /* Replace last character of filename with
                         ** 0_asn.fits	*/
                        strcat (filename, "0_asn.fits");
                        exist = DoesFileExist (filename);

                        if (exist == EXISTS_YES) {
                            strcpy (asn->asn_table, filename);
                            if (asn->verbose) {
                                trlmessage("PARTIAL processing of ASN table...");
                            }
                            asn->process = PARTIAL;
                        } else {
                            /* We can't find a matching file or ASN
                             ** table to process */
                            trlerror("File %s not found for processing", filename);
                            return (status = OPEN_FAILED);
                        }
                    }
                }			
            } /* End of Exists ELSE statement for this case */

            break; /* End of Case 3 and DEFAULT case */
    } /* End of switch statement */

    if (asn->debug) {
        trlmessage("SetInput: Determined what ASN table to process and how.");
    }

    return (status);
} 

/* GETASNTABLE: Read the Association (ASN) table and load the
 ** list of member dataset names. Also check to see if any of the
 ** members are listed as missing. 
 */

int GetAsnTable (AsnInfo *asn) {

    /* Arguments:
     **	asn		o: Association info structure
     */
    extern int status;

    /* Local variables */
    int i;				/* loop index */
    int nrows;			/* number of rows in ASNTAB */
    int col, row;			/* loop indexes */
    IRAFPointer tp;			/* ASNTAB table pointer */
    IRAFPointer colptr[NCOLS];	/* ASNTAB column pointers */
    int numsp;			/* number of sub-products in asn */
    int	posid;			/* id of sub-product */
    RowInfo *exp;		/* Internal structure for table information */
    int *spmems, *expmem;	/* number of EXP for each sub-product */
    int poslen;		/* Length of memtype string minus POSID */
    int prodid;		/* id of product 	*/
    int expid;
    int defid;		/* default position id for exposures */
    int procprod;
    int maxspmems;		/* max number of sub-product members */
    char *word;
    Bool proddth;
    char spname_ext[SZ_CBUF+1];	/* Extension for sub-product */
    char spname_ext_cte[SZ_CBUF+1];	/* Extension for sub-product */
    char memtype[SZ_COLNAME+1];
    char memsubtype[SZ_COLNAME+1];

    /* ASNTAB column names */
    char colname[NCOLS][SZ_COLNAME+1] = {
        "MEMNAME",
        "MEMTYPE",
        "MEMPRSNT"
    };

    /* Function definitions */
    void freeAsnInfo (AsnInfo *);
    int  allocAsnInfo (AsnInfo *, int, int *);
    char *lowcase (char *, char *);
    int  MkName (char *, char *, char *, char *, char *, int);
    void initRowInfo (RowInfo *);
    int  streq_ic (char *, char *);  /* strings equal? (case insensitive) */

    if (asn->debug) {
        trlmessage("GetAsnTable: ASN_TABLE is %s",asn->asn_table);
    }

    /* Open the ASN table */
    tp = c_tbtopn (asn->asn_table, IRAF_READ_ONLY, 0);
    if (c_iraferr()) {
        trlopenerr (asn->asn_table);
        return (status = TABLE_ERROR);
    }

    /* Get pointers to columns in ASNTAB */
    for (col=0; col < NCOLS; col++) {
        c_tbcfnd1 (tp, colname[col], &(colptr[col]));
        if (c_iraferr() || colptr[col] == 0) {
            trlerror("Can't find column %s in %s", colname[col], asn->asn_table);
            c_tbtclo (tp);
            return (status = COLUMN_NOT_FOUND);
        }
    }

    /* Find out how many rows are in ASNTAB */
    nrows = 0;
    nrows = c_tbpsta (tp, TBL_NROWS);
    if (nrows <= 0) {
        trlerror("Invalid number of rows in %s", asn->asn_table);
        c_tbtclo (tp);
        return (status = TABLE_ERROR);
    }

    /* Initialize total number of members based on number of
     ** rows in table */
    asn->numasn = nrows;
    poslen = 0;

    /* Allocate space for internal variables */
    exp = NULL;
    exp = (RowInfo *)calloc(nrows, sizeof(RowInfo));
    for (row = 0; row < nrows; row++)
        initRowInfo (&(exp[row]));

    /* Read each row of ASNTAB into a local structure */
    for (row = 0; row < nrows; row++) {

        /* Get the MEMBER NAME in this row */
        c_tbegtt (tp, colptr[0], row+1, exp[row].memname, SZ_CBUF);
        if (c_iraferr()) {
            trlerror("Can't read %s in row %d in %s", colname[0], row+1, asn->asn_table);
            c_tbtclo (tp);
            free (exp);
            return (status = ELEMENT_NOT_FOUND);
        }

        /* Convert to lowercase for use as a file name */
        for (i = 0; i < strlen(exp[row].memname); i++)
            exp[row].memname[i] = tolower(exp[row].memname[i]);

        /* Get the TYPE in this row */
        c_tbegtt (tp, colptr[1], row+1, exp[row].mtype, SZ_CBUF);
        if (c_iraferr()) {
            trlerror("Can't read %s in row %d in %s", colname[1], row+1, asn->asn_table);
            c_tbtclo (tp);
            free (exp);
            return (status = ELEMENT_NOT_FOUND);
        }

        /* Convert to lowercase for use later in routine.
         ** Also, if a value of MEMTYPE contains a UNDERLINE, then 
         ** record conversion to DASH in trailer file and warn
         ** user to correct the value. */
        lowcase (exp[row].type, exp[row].mtype);
        for (i = 0; i < strlen(exp[row].type); i++) {
            if (exp[row].type[i] == UNDERLINE_CHAR) {       
                exp[row].type[i] = DASH_CHAR;
                trlwarn("MEMTYPE %s in row %d was INVALID and needs to be corrected.", exp[row].mtype, row+1);
            }
        }

        /* Get the STATUS in this row */
        c_tbegtb (tp, colptr[2], row+1, &(exp[row].prsnt));
        if (c_iraferr()) {
            trlerror("Can't read %s in row %d in %s", colname[2], row+1, asn->asn_table);
            c_tbtclo (tp);
            free (exp);
            return (status = ELEMENT_NOT_FOUND);
        }
        if (asn->debug) {
            trlmessage("GetAsnTable: Read in row %d from ASN table", row);
        }
    }

    /*  
     ** Determine whether CRCORR or RPTOBS processing will be required
     ** by searching for MEMTYPE of EXP_CR* or EXP_RPT*, respectively.
     ** Once it is determined, go on to next step...
     */
    for (row = 0; row < nrows; row++) {

        /* Check to see if we need to do CR-SPLIT combination ... */
        if (strstr (exp[row].type, "-cr") != NULL) {
            asn->crcorr = PERFORM;
            asn->rptcorr = OMIT;
            poslen = CRLEN;
            break;

            /* ... or REPEAT-OBS combination ... 	*/
        } else if (strstr (exp[row].type, "-rp") != NULL) { 
            asn->rptcorr = PERFORM;
            asn->crcorr = OMIT;
            poslen = RPTLEN;
            break;

            /* ... or neither at all	*/
        } else {
            asn->rptcorr = OMIT;
            asn->crcorr = OMIT;
            poslen = DTHLEN;
        }	
    }

    /* Default to always perform DRIZCORR step */
    asn->dthcorr = PERFORM;

    if (asn->debug) {
        trlmessage("GetAsnTable: CRCORR = %d, RPTCORR = %d", asn->crcorr,asn->rptcorr);
    }

    /* Sort through the list figuring out which are input vs. output
     ** files, and see if any input files are missing. */
    /* Initialize local variables */
    numsp   = 0;
    posid   = 0;
    prodid  = 0;
    proddth = False;
    defid   = 1;

    /* Find out how many products/sub-products are in the association
     ** and determine the posid for each member. */	
    for (row = 0; row < nrows; row++) {
        memtype[0] = '\0';
        memsubtype[0] = '\0';

        /* As long as this is not the final product,
         ** count number of members in each product/sub-product...*/
        if (strstr(exp[row].type,"prod-dth") != NULL) {
            exp[row].posid = 0;	
            posid = 0;	

            /* If we have a dither product listed, we want to eventually
             ** perform dither combining step... */	
            prodid++;
            proddth = True;

            /* We always want to produce a product, but if
             ** not set to PERFORM, then create an empty product... */
            if (asn->dthcorr == OMIT) 
                asn->dthcorr = DUMMY;

        } else {

            strcpy(memtype,exp[row].type);	
            /* Let's start by breaking up the MEMTYPE string */
            word = strtok(memtype,"-");
            strcpy(memsubtype,word);
            /* The end of this second part of MEMTYPE has the POSID */
            word = strtok(NULL,"\n");

            /* If the last character of the memtype is a number or letter, 
             ** convert to an integer value...  */
            if (streq_ic(word,"crj") || streq_ic(word,"rpt") || streq_ic(word,"crc") ) {
                posid = 1;
            } else {
                if (exp[row].prsnt || streq_ic(memsubtype,"prod")) {

                    if (isalnum(word[poslen])) {
                        /* Interpret the POSID from the second part
                         ** of MEMTYPE */
                        posid=(int)strtol(&word[poslen],(char **)NULL,10);
                    } else {
                        /* Otherwise, assume it is a different pointing
                         ** and assign an incremented position ID. After
                         ** all, it is NOT a CR-SPLIT or RPT-OBS exposure.*/
                        posid = defid;
                        defid++;
                    }
                }
            }

            /* Keep track of how many sub-products there are based on
             ** POSID from MEMTYPE values */

            if (posid > numsp) {
                if ((exp[row].prsnt && (strstr(memtype,"exp") != NULL)) ||
                        strstr(memtype,"prod") != NULL) {
                    numsp++;
                }
            }

            exp[row].posid = posid;		
        }

        if (asn->debug) {
            trlmessage("GetAsnTable: Posid = %d for row %d",posid, row);
        }

        /* If the member is missing, give a warning */
        /* This implies that MEMPRSNT must be set to YES for EXP_*
         ** This is not fatal for WF3.  Simply decrement tmembers. */
        if (!exp[row].prsnt && strncmp(exp[row].type, "prod-",5) != 0) {
            trlwarn("Member \"%s\" is not present", exp[row].memname);
            asn->numasn--;
            /* Now flag row as being absent so it doesn't get passed
             ** along for further processing... */
            exp[row].posid = MEMABSENT;
        }
    }

    if (asn->debug) {
        trlmessage("GetAsnTable: NUMSP = %d, PRODID = %d",numsp, prodid);
    }

    /* Check for existence of enough data to process */
    if (asn->numasn < 1) {
        trlerror("No data available for given assoc. table");
        freeAsnInfo (asn);
        return (status = ERROR_RETURN);
    }

    /* Record the number of products found in association */
    if (proddth) {
        asn->numprod = prodid;
    } else {
        /* If there are no PROD-DTH entries in ASN table, set
         ** numprod to 1 */
        asn->numprod = prodid + 1;
    }

    /* Determine what elements should be processed based on
     ** initial input, either FULL or PARTIAL processing.  */
    procprod = 0;
    if (asn->process != FULL) {
        for (row=0; row < nrows; row++){
            /* Look for entry with same name as input */
            if (streq_ic(exp[row].memname, asn->rootname)) {
                /* We only want to process this product	*/
                procprod = exp[row].posid;
                numsp = 1;
            }	
        }
    }

    /* Allocate space to count number of exposures per sub-product
     ** with spmems[0] being reserved for final output product
     ** since POSID starts indexing at 1. */
    spmems = (int *)calloc(numsp+1,sizeof(int));
    expmem = (int *)calloc(numsp+1,sizeof(int));

    /* Initialize arrays */
    for (i=0; i <= numsp; i++) {
        spmems[i] = 0;
        expmem[i] = 0;
    }

    /* For each sub-product,
     ** identify each EXP that belongs to that posid and 
     ** is present to be processed. */
    for (row=0; row < nrows; row++) {
        if (strstr(exp[row].type, "exp-") && exp[row].posid > MEMABSENT) {

            if ((asn->process != FULL && exp[row].posid == procprod) ||
                    asn->process == FULL) {

                spmems[exp[row].posid]++;
                /* Exposure IDs will start at 1 to be consistent 
                 ** with POSID numbering. Initialize here, count later. */
                expmem[exp[row].posid] = 1;
            }
        }
    }

    /* Allocate slots for all members in ASN info structure */
    if (allocAsnInfo (asn, numsp, spmems)) {
        return (status = TABLE_ERROR);
    }

    asn->product[prodid].prodid = prodid;

    /* Reset prodid for filling ASN table */
    prodid = 0;

    /* Copy summary information about ASN relationships to ASN structure. */
    maxspmems = 0;
    asn->numsp = numsp;
    for (i=0; i <= numsp; i++) {
        asn->spmems[i] = spmems[i];
        if (spmems[i] > maxspmems)
            maxspmems = spmems[i];
    }

    /* If there aren't any sub-products with more than 1 member, then
     ** omit crcorr/rptcorr processing */
    if ((maxspmems < 2) && asn->crcorr == PERFORM) {
        trlwarn("No sub-product with more than 1 member; CRCORR will be skipped");
        asn->crcorr = DUMMY;
    } else if ((maxspmems < 2) && asn->rptcorr == PERFORM) {
        trlwarn("No sub-product with more than 1 member; RPTCORR will be skipped");
        asn->rptcorr = DUMMY;
    }

    /* Copy information read-in into ASN structure now... */
    for (row = 0; row < nrows; row++) {
        if (exp[row].posid != MEMABSENT) {
            if ((asn->process != FULL && exp[row].posid == procprod) ||
                    asn->process == FULL) {	  
                posid = exp[row].posid;

                /* Is this row the final product entry? */
                if (strstr(exp[row].type, "prod-dth") != NULL) {
                    strcpy (asn->product[prodid].name, exp[row].memname);
                    strcpy (asn->product[prodid].mtype, exp[row].type);
                    asn->product[prodid].prsnt = exp[row].prsnt;
                    asn->product[prodid].numsp = numsp;
                    asn->product[prodid].asnrow = row+1;

                    /* Create full file name for this image */
                    if (MkName (exp[row].memname, "_raw", "_drz", "",
                                asn->product[prodid].prodname, CHAR_LINE_LENGTH)){
                        strcpy(asn->product[prodid].prodname,
                                exp[row].memname);
                        strcat (asn->product[prodid].prodname,
                                "_drz.fits");
                    }
                    
                    /* Create full file name for this CTE image */
                    if (MkName (exp[row].memname, "_rac_tmp", "_drc", "",
                                asn->product[prodid].prodname_cte, CHAR_LINE_LENGTH)){
                        strcpy(asn->product[prodid].prodname_cte,
                                exp[row].memname);
                        strcat (asn->product[prodid].prodname_cte,
                                "_drc.fits");
                    }
                    

                    /* Or, is this row an input exposure? */
                } else if (strstr(exp[row].type, "exp-") != NULL) {

                    if (exp[row].posid > MEMABSENT) {

                        /* Internal counter for which exposure we want for
                         ** a position */
                        expid = expmem[posid];
                        strcpy(asn->product[prodid].subprod[posid].exp[expid].name,
                                exp[row].memname);		
                        strcpy(asn->product[prodid].subprod[posid].exp[expid].mtype,
                                exp[row].type);
                        asn->product[prodid].subprod[posid].exp[expid].prsnt =
                            exp[row].prsnt;
                        asn->product[prodid].subprod[posid].exp[expid].asnrow=
                            row+1;

                        /* Create full file name for this image */
                        if (MkName (exp[row].memname, "_raw", "_raw", "",
                                    asn->product[prodid].subprod[posid].exp[expid].expname,
                                    CHAR_LINE_LENGTH)) {

                            strcpy(asn->product[prodid].subprod[posid].exp[expid].expname,
                                    exp[row].memname);
                            strcat(asn->product[prodid].subprod[posid].exp[expid].expname,
                                    "_raw.fits");
                        }

                        /* Fill-in sub-product information for EXP-DTH
                         ** exposures which don't create sub-products */
                        if (strstr(exp[row].type, "exp-dth") != NULL) {
                            if (!MkName (exp[row].memname, "_raw", "_flt", "",
                                        asn->product[prodid].subprod[posid].spname,
                                        CHAR_LINE_LENGTH) ) {

                                strcpy(asn->product[prodid].subprod[posid].name,
                                        exp[row].memname);
                                strcpy(asn->product[prodid].subprod[posid].mtype,
                                        exp[row].type);
                                asn->product[prodid].subprod[posid].posid =
                                    exp[row].posid;
                            }
                        }

                        /* Increment that counter for next exposure's id
                         ** for this posid */
                        expmem[posid]++;
                    }

                    /* If neither, it must be a sub-product */
                } else {

                    if (spmems[posid] > 0) {

                        strcpy(asn->product[prodid].subprod[posid].name,
                                exp[row].memname);
                        strcpy(asn->product[prodid].subprod[posid].mtype,
                                exp[row].type);
                        asn->product[prodid].subprod[posid].prsnt =
                            exp[row].prsnt;
                        asn->product[prodid].subprod[posid].asnrow = row+1;


                        /* Create full file name for this image for 
                         ** DTHCORR input */
                        spname_ext[0] = '\0';
                        spname_ext_cte[0] = '\0';
                        if (asn->crcorr || asn->rptcorr) {
                            strcpy (spname_ext, "_crj");
                            strcpy (spname_ext_cte,"_crc");
                        } else {
                            strcpy (spname_ext, "_sfl");
                            strcpy (spname_ext_cte, "_sfl");
                        }

                        if (MkName (exp[row].memname, "_raw", spname_ext, "",
                                    asn->product[prodid].subprod[posid].spname,
                                    CHAR_LINE_LENGTH)) {

                            strcpy(asn->product[prodid].subprod[posid].spname,
                                    exp[row].memname);
                            strcat(asn->product[prodid].subprod[posid].spname,
                                    spname_ext);		
                            strcat(asn->product[prodid].subprod[posid].spname,
                                    ".fits");		
                                    
                        }
                        if (MkName (exp[row].memname, "_raw", spname_ext_cte, "",
                                    asn->product[prodid].subprod[posid].spname_cte,
                                    CHAR_LINE_LENGTH)) {

                            strcpy(asn->product[prodid].subprod[posid].spname_cte,
                                    exp[row].memname);
                            strcat(asn->product[prodid].subprod[posid].spname_cte,
                                    spname_ext_cte);		
                            strcat(asn->product[prodid].subprod[posid].spname_cte,
                                    ".fits");		
                                    
                        }

                        asn->product[prodid].subprod[posid].numexp =
                            spmems[posid];
                        asn->product[prodid].subprod[posid].posid = posid;

                    }
                }
            }
        } /* Process only those exposures where MEMPRSNT == YES */
    }

    /* Close the ASN table.  We are done reading it in. */
    c_tbtclo (tp);



    /* Clean up memory usage as well. */
    free (spmems);
    free (expmem);
    free (exp);

    if (asn->debug) {
        trlmessage("GetAsnTable: Info from ASN read into memory.");
    }

    /* Successful return */
    return (status);
}


void initRowInfo (RowInfo *exp) {

    exp->memname[0] = '\0';
    exp->mtype[0]   = '\0';
    exp->prsnt      = False;
    exp->posid      = 0;
    exp->type[0]    = '\0';

}

void initAsnInfo (AsnInfo *asn) {

    asn->input[0]     = '\0';
    asn->asn_table[0] = '\0';
    asn->process      = -1;
    asn->crcorr       = DUMMY;
    asn->rptcorr      = DUMMY;
    /* for dthcorr, OMIT == never do, DUMMY == produce dummy product */
    asn->dthcorr      = OMIT; 
    asn->numprod      = 0;
    asn->numsp        = 0;
    asn->spmems       = NULL;
    asn->numasn       = 0;
    asn->product      = NULL;
    asn->instr[0]     = '\0';
    asn->detector     = UNKNOWN_DETECTOR;
    asn->verbose      = 0;
    asn->debug        = 0;

}

int allocAsnInfo (AsnInfo *asn, int numsp, int *spmems) {

    extern int status;

    /* Local variables */
    int i;		/* loop index */
    int prodid = 0;	/* product ID: for looping over products later */
    int j;		/* loop index for EXPs */
    int numexp;	

    /* Function definitions */
    void initAsnProduct (ProdInfo *, int);
    void initAsnSubProd (SubProdInfo *, int);
    void initAsnExp (ExpInfo *);

    void freeAsnInfo (AsnInfo *);
    void initAsnExp (ExpInfo *);

    /* Free the structure if necessary */
    if (asn->product != NULL)
        freeAsnInfo(asn);

    /* Allocate the member structures */
    asn->spmems  = (int *)calloc(numsp+1,sizeof(int));
    asn->product = (ProdInfo *)calloc(1,sizeof(ProdInfo ));

    /* Initialize each member structure */
    for (i=0; i < asn->numprod; i++) {
        initAsnProduct (&(asn->product[i]), numsp+1);
    }

    asn->product[prodid].subprod =
        (SubProdInfo *)calloc(numsp+1, sizeof(SubProdInfo));

    for (i=0; i <= numsp; i++) {
        if (spmems[i] > 0) 
            numexp = spmems[i];
        else
            numexp = 1;

        initAsnSubProd(&(asn->product[prodid].subprod[i]), numexp);
    }

    for (i=0; i <= numsp; i++) {
        asn->product[0].subprod[i].exp =
            (ExpInfo *)calloc(spmems[i]+1, sizeof(ExpInfo));
        for (j=0; j <= spmems[i]; j++)
            initAsnExp (&(asn->product[prodid].subprod[i].exp[j]));	
    }

    /* Check for error during allocation */
    if (asn->product == NULL) {
        asn->numprod = 0;
        asn->numasn = 0;
        trlerror("Insufficient memory to allocate ASN structure");
        return (status = 1);
    }

    /* Succesful return */
    return (status);
}

void freeAsnInfo (AsnInfo *asn) {

    /* Local variables */
    int i, j, p;			/* loop index */
    int numexp;

    /* Function definitions */
    void initAsnProduct (ProdInfo *, int);
    void initAsnInfo (AsnInfo *);
    void initAsnSubProd (SubProdInfo *, int);
    void initAsnExp (ExpInfo *);

    /* If ASN table was read in, free it. */
    if (asn->product != NULL) {
        for (p = 0; p < asn->numprod; p++) {
            for (i=0; i <= asn->numsp; i++) {
                for (j=0; j <= asn->spmems[i]; j++) {
                    initAsnExp (&(asn->product[p].subprod[i].exp[j]));
                }
                free (asn->product[p].subprod[i].exp);
                if (asn->spmems[i] > 0) 
                    numexp = asn->spmems[i];
                else
                    numexp = 1;
                initAsnSubProd (&(asn->product[p].subprod[i]), numexp);
            }
            free (asn->product[p].subprod);
            initAsnProduct (&(asn->product[p]), asn->numsp+1);
        }
        free (asn->product);
    }
    free (asn->spmems);

    initAsnInfo (asn);

}

void initAsnProduct (ProdInfo *product, int numsp) {

    product->name[0]     = '\0';
    product->mtype[0]    = '\0';
    product->prsnt       = False;
    product->asnrow      = 0;
    product->numsp       = numsp;
    product->prodid      = 1;
    product->prodname[0] = '\0';
    product->prodname_cte[0] = '\0';
    product->subprod     = NULL;

}

void initAsnSubProd (SubProdInfo *subprod, int numexp) {

    subprod->name[0]    = '\0';
    subprod->mtype[0]   = '\0';
    subprod->prsnt      = False;
    subprod->asnrow     = 0;
    subprod->spname[0]  = '\0';
    subprod->numexp     = numexp;
    subprod->spname_cte[0]  = '\0';
    subprod->crj_tmp[0] = '\0';
    subprod->crc_tmp[0] = '\0';
    subprod->exp        = NULL;
}

void initAsnExp (ExpInfo *exp) {

    exp->name[0]    = '\0';
    exp->mtype[0]   = '\0';
    exp->prsnt      = False;
    exp->asnrow     = 0;
    exp->expname[0] = '\0';
    exp->blv_tmp[0] = '\0';
    exp->blc_tmp[0] = '\0';
    exp->rac_tmp[0] = '\0';
    exp->dx = 0;
    exp->dy = 0;
    exp->xi = 0;
    exp->yi = 0;
}

void initWCS (WCS *wcs) {

    wcs->crpix[0] = 0;
    wcs->crpix[1] = 0;
    wcs->crval[0] = 0;
    wcs->crval[1] = 0;
    wcs->cd[0][0]  = 0;
    wcs->cd[0][1]  = 0;
    wcs->cd[1][0]  = 0;
    wcs->cd[1][1]  = 0;
    wcs->ctype[0][0] = '\0';
    wcs->ctype[1][0] = '\0';

}

/* GETGLOBALINFO: Reads observation flags and indicator keyword values 
 ** from association table primary header.
 **	This task will read the following keywords from the ASN table's primary
 **	header:  INSTRUME, DETECTOR
 */

int GetGlobalInfo (AsnInfo *asn) {

    /* Arguments:
     **	asn	io: association info structure
     */
    extern int status;
    Hdr phdr;               /* primary header */

    char detector[SZ_FITS_REC+1];

    /* Function definitions */
    int GetKeyStr (Hdr *, char *, int, char *, char *, int);
    int LoadHdr (char *, Hdr *);
    /*	int GetSwitch (Hdr *, char *, int *); */

    if (asn->debug) {
        trlmessage("GetGlobalInfo: Ready to open primary header... ");
    } 

    if (asn->debug) {
        trlmessage("GetGlobalInfo: asn_table is %s",asn->asn_table);
    }

    /* Read primary header of ASN file into phdr. */
    if (LoadHdr (asn->asn_table, &phdr)) {
        trlerror("Could not load header from table %s", asn->asn_table);
        return (status);
    }
    if (asn->debug) {
        trlmessage("GetGlobalInfo: Read in header from Image");
    }	

    /* Get the observing mode keyword values from header */
    asn->instr[0] = '\0';
    if (GetKeyStr (&phdr, "INSTRUME", 0, "", asn->instr, SZ_FITS_REC)) {
        trlkwerr ("INSTRUME", asn->asn_table);
        return (status = KEYWORD_MISSING);
    }

    asn->detector = 0;
    detector[0]   = '\0';
    if (GetKeyStr (&phdr, "DETECTOR", 0, "", detector, SZ_FITS_REC)) {
        trlkwerr ("DETECTOR", asn->asn_table);
        return (status = KEYWORD_MISSING);
    }

    /* Convert detector string to usable value */
    if (strncmp (detector, "UVIS", 4) == 0) {
        asn->detector = CCD_DETECTOR;
    } else if (strncmp (detector, "IR", 2) == 0) {
        asn->detector = IR_DETECTOR;
    } else {
        asn->detector = UNKNOWN_DETECTOR;
        return (status = HEADER_PROBLEM);
    }

    checkGlobalInfo(asn);

    /* You can NOT create a summed image with only 1 input */
    if (asn->process == SINGLE) {
        asn->rptcorr = OMIT;
    }

    /* If we are not processing an entire association, then
     ** we will not have the inputs necessary for a DTHCORR. */
    /*	if (dthcorr == PERFORM) { */
    if (asn->process == SINGLE) {
        asn->dthcorr = DUMMY;
    }
    /*	} */

    /* Otherwise, leave asn->dthcorr as set by reading ASN table itself */

    /* Close the ASN table's primary header here. */
    freeHdr (&phdr);



    if (asn->debug) {
        trlmessage("GetGlobalInfo: Detector and Instrument determined");
    }	

    /* Successful return */
    return (status);
}

/* CHECKGLOBALINFO: Checks flags and indicators for validity. */

int checkGlobalInfo (AsnInfo *asn) {

    /* Arguments:
     **	asn	io: Association info structure
     */
    extern int status;

    /* Check instrument = WF3 */
    if (strncmp (asn->instr, "WFC3", 4) != 0 ) {
        trlerror("INSTRUME keyword value \"%s\" not valid in %s", asn->instr, asn->filename);
        status = 1;
    }

    /* Check for valid camera number */
    if (asn->detector < 1 || asn->detector > 2) {
        trlerror("CAMERA keyword value \"%d\" not valid in %s", asn->detector, asn->filename);
        status = 1;
    }

    return (status);
}

void printInfo (AsnInfo *asn) {

    int i, j, k;
    int numprod;

    trlmessage("");

    if (asn->dthcorr == DUMMY) {
        numprod = 0;
    } else {
        numprod = asn->numprod;
    }
    trlmessage("NUMBER of MEMBERS in TABLE: %d  PRODUCTS: %d  SUB-PRODUCTS: %d", asn->numasn, numprod, asn->numsp);

    if (asn->process != SINGLE) {
        for (i=0; i < asn->numprod; i++) {
            if (asn->dthcorr != DUMMY || asn->dthcorr != OMIT) {
                trlmessage("Product-- Member %3d: %s  Product: %2d  Type: %s", i+1,
                    asn->product[i].name, asn->product[i].prodid, asn->product[i].mtype);
            }

            for (j = 1; j <= asn->numsp; j++) {
                trlmessage("Sub-Product-- Member %3d: %s  Posn: %2d  Type: %s",
                        j, asn->product[i].subprod[j].name, 
                        asn->product[i].subprod[j].posid,
                        asn->product[i].subprod[j].mtype);

                for (k = 1; k <= asn->spmems[j]; k++) {
                    trlmessage("Exposure-- Member %3d: %s  Type: %s",
                            k, asn->product[i].subprod[j].exp[k].name, 
                            asn->product[i].subprod[j].exp[k].mtype);
                }
            }
        }

    } else {

        /* Print out info for SINGLE exposure ... */
        trlmessage("Exposure-- Processing SINGLE Exposure %s ", asn->filename);
    }

    trlmessage("");
}

static int IsProduct (char *input) {

    int i;
    /* names of intermediate product suffixes */
    char prodsuf[NSUF][SZ_CBUF+1] = {
        "_crj", "_crj_tmp",
        "_blv", "_blv_tmp",
        "_blc", "_blc_tmp",
        "_crc", "_crc_tmp",
        "_flc",
        "_rac_tmp",
        "_flt", "_sfl",
        "_drz",
        "_drc",
    };

    for (i=0; i < NSUF; i++) {
        if (strstr(input, prodsuf[i]) != NULL)
            return(0);
    }
    return (1);
}


int GetAsnName (char *filename, char *asn_name) {

    extern int status;
    IODescPtr im;           /* descriptor for an image */
    Hdr phdr;               /* primary header */

    /* Function definitions */
    int GetKeyStr (Hdr *, char *, int, char *, char *, int);

    /* Read primary header of ASN file into hdr. */
    initHdr (&phdr);
    im = openInputImage (filename, "", 0);
    if (hstio_err())
        return (status = OPEN_FAILED);

    getHeader (im, &phdr);          /* get primary header */
    if (hstio_err())
        return (status = OPEN_FAILED);

    closeImage (im);

    asn_name[0] = '\0';
    if (GetKeyStr (&phdr, "ASN_TAB", 0, "", asn_name, SZ_FITS_REC)) {
        trlkwerr ("ASN_TAB", asn_name);
        return (status = KEYWORD_MISSING);
    }

    /* Close the file's primary header. */
    freeHdr (&phdr);

    /* Successful return */
    return (status);

}

/* updateAsnTable: Write the Association (ASN) table, updating
 **	the particular sub-product/product MEMPRSNT information.
 */

int updateAsnTable (AsnInfo *asn, int prodid, int posid) {

    /* Arguments:
     **	asn		    i: association info structure
     **	prodid		i: product id for member that needs to be updated
     **	posid		i: position id for member that needs to be updated
     */
    extern int status;

    /* Local variables */
    int col;			/* loop index */
    IRAFPointer asn_tp;		/* table pointers */
    IRAFPointer colptr[NCOLS];	/* ASN table column pointers */
    Bool prsnt;

    /* ASNTAB column names */
    char colname[NCOLS][SZ_COLNAME+1] = {
        "MEMNAME",
        "MEMTYPE",
        "MEMPRSNT"
    };

    /* Initialize the table column pointers */
    for (col = 0; col < NCOLS; col++)
        colptr[col] = 0;

    /* Open the ASN table */
    asn_tp = c_tbtopn (asn->asn_table, IRAF_READ_WRITE, 0);
    if (c_iraferr()) {
        trlopenerr (asn->asn_table);
        return (status = TABLE_ERROR);
    }

    /* Find the columns in the ASN table */
    for (col=0; col < NCOLS; col++) {
        if (colptr[col] == 0) {
            c_tbcfnd1 (asn_tp, colname[col], &(colptr[col]));
            if (c_iraferr() || colptr[col] == 0) {
                trlerror("Can't find column %s in %s", colname[col], asn->asn_table);
                c_tbtclo (asn_tp);
                return (status = 1);
            }
        }
    }

    /* Value MEMPRSNT will be updated to... */
    prsnt = True;

    /* Are we working with a Product or a subproduct... */
    if (posid == 0) {

        /* Write the updated info for PRODUCT to the ASN table */
        c_tbeptb (asn_tp, colptr[2], asn->product[prodid].asnrow, prsnt);

    } else {

        /* Write the updated info for SUB-PRODUCT to the ASN table */
        c_tbeptb (asn_tp, colptr[2],
                asn->product[prodid].subprod[posid].asnrow, prsnt);
    }

    /* Close the ASN table */
    c_tbtclo (asn_tp);

    /* Update ASN_PROD keyword to TRUE to signal successful
     ** creation of product/sub-product. */
    if (UpdateHdr (asn->asn_table) ) {    
        trlerror("Couldn't update ASN table header");
        return(status = KEYWORD_MISSING);
    }

    /* Successful return */
    return (status);
}


static int UpdateHdr (char *output) {

    extern int status;

    Hdr phdr;               /* primary header */
    IODescPtr im;		/* descriptor for output image */

    int PutKeyBool (Hdr *, char *, Bool, char *);

    trlmessage("Trying to open %s...",output);
    initHdr (&phdr);

    /* Open input image in order to read its primary header. */
    im = openUpdateImage (output, "", 0, &phdr);				
    if (hstio_err()) {
        trlopenerr (output);
        closeImage(im);
        return (status = OPEN_FAILED);
    }

    if (PutKeyBool (&phdr, "ASN_PROD", True, "") ) {
        freeHdr (&phdr);
        trlerror("Couldn't update ASN_PROD keyword in ASN table header");
        return(status = KEYWORD_MISSING);
    }

    /* write out primary header */
    if (putHeader (im))
        status = HEADER_PROBLEM;	
    if (hstio_err() || status) {
        trlreaderr (output);
        closeImage (im);
        return (status = OPEN_FAILED);
    }

    closeImage (im);
    /* Close the ASN table's primary header. */
    freeHdr (&phdr);

    trlmessage("Updated Global Header for %s...",output);

    return (status);

}

