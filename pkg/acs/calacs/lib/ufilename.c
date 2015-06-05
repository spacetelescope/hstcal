# include <stdio.h>
# include <ctype.h>
# include <string.h>
# include "hstio.h"

/* This routine updates the FILENAME primary header keyword, or adds it
   to the header if it's not already present.  If the input file name
   begins with an environment variable or an explicit directory, that will
   be skipped over before writing the name to the header.
   The section for finding a directory prefix is based on the iraf zfnbrk
   function.
*/

void UFilename (char *filename, Hdr *phdr) {

/* arguments:
char *filename  i: file name, but may begin with directory specification
Hdr *phdr       io: pointer to header; FILENAME will be updated
*/

	int ch;			/* one character in filename */
	int namelen;		/* length of filename */
	int start = 0;		/* start of file name */
	int i;
	int PutKeyStr(Hdr *, char *, char *, char *);

	namelen = strlen (filename);

	/* If there's a directory prefix, skip over it. */
	for (i = 0;  i < namelen;  i++) {
	    ch = filename[i];
	    if (isalnum(ch) || ch == '_' || ch == '.')
			continue;			/* ordinary character */
	    else if (ch == '\\' && i+1 < namelen)
			i++;			/* skip \ and next character */
	    else
			start = i + 1;	/* not ordinary, so part of directory prefix */
	}

	/* Were we unable to find the filename root? */
	if (start >= namelen - 1)
	    start = 0;

	/* Update the FILENAME keyword, or add it if it doesn't exist. */
	PutKeyStr (phdr, "FILENAME", &filename[start], "name of file");
}

/* This function updates the ASN_MTYPE keyword to reflect the different
    roles each image has in the association.  This primarily affects only
    output products are created from ACSREJ, ACSSUM, and ACSDTH, where
    the ASN_MTYP is different from the input file's value.
    
    WJH 6 May 1999 
*/
void UMemType (char *mtype, Hdr *phdr) {

/* arguments:
char *mtype     i: memtype from association table
Hdr *phdr       io: pointer to header; ASN_MTYP will be updated
*/

    int len;
    char u_mtype[10];
    
    void UpperAll (char *, char *, int);
    int PutKeyStr (Hdr *, char *, char *, char *);
    
    len = 0;
    u_mtype[0] = '\0';
    
    if (mtype[0] != '\0') {
        len = strlen(mtype);
        UpperAll (mtype, u_mtype, len+8);
        PutKeyStr (phdr, "ASN_MTYP", u_mtype, "Role of the Exposure in the Association");
    }
}

/*
    This function updates the EXPNAME keyword in an image header.
	WJH 10May 1999
*/
void UExpname (char *expname, Hdr *phdr) {

/* arguments:
char *expname     i:expname for observation (ROOTNAME)
Hdr *phdr       io: pointer to header; EXPNAME will be updated
*/

   
    int PutKeyStr (Hdr *, char *, char *, char *);
    
    PutKeyStr (phdr, "EXPNAME", expname, "9 character exposure identifier");

}
