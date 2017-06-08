# include <stdio.h>
# include <string.h>

# include "msg.h"

# define NSUF	16

void FindAsnRoot (const char *input, char *root) {

	int i;
	int in_len;
	char filename[SZ_FNAME+1];
	
	/* names of WF3 product suffixes */
	char prodsuf[NSUF][SZ_CBUF+1] = {
		"_raw", "_asn", 
		"_crj", "_crj_tmp",
		"_blv", "_blv_tmp",
        "_blc", "_blc_tmp",
        "_crc", "_crc_tmp",
        "_flc",
		"_flt", "_sfl",
		"_drz",
        "_drc", 
		".fit"
	};

	/* Initialize local variable... */
	filename[0] = '\0';
	
	/* We do not want to alter the input... */
	strcpy (filename, input);
	
	/* Look for each suffix recognized by WF3 in filename... */
	for (i = 0; i < NSUF ; i++){
		/* Is this suffix in filename? 
		**	If so, set in_len to the length of filename up to that point
		**	If not, set in_len to the length of filename
		*/
		in_len = 0;
				
		if (strstr(filename, prodsuf[i]) == NULL) {
			in_len = strlen (filename);
		} else {
			in_len = strlen(filename) - strlen(strstr(filename, prodsuf[i]));
		}
			
		/* If in_len indicates the suffix was found, copy out just
		** that portion of filename up to that suffix...
		** otherwise, set rootname to be the original filename
		**	as the original apparently didn't have one.
		*/
		if (in_len < strlen(filename)) {
			filename[in_len] = '\0';		
			strcpy (root, filename);
			break;
		} else {
		/* If it doesn't find a suffix, filename must not have one,
		**	so, return original filename unchanged...
		*/
			strcpy (root, filename);
		}
		
	} /* Go on to the next suffix, checking for '.fits' last. */
}

