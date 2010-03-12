# include <stdlib.h>
# include <string.h>
# include "ctables.h"

char *expandfn (char *filename) {

	int i, j;
	int doller, start;
	int done;
	char *var;
	char *value, *l_value;

	var = (char *)calloc (SZ_FNAME, sizeof(char));
	l_value = (char *)calloc (SZ_FNAME, sizeof(char));

	i = 0;
	start = 0;
	doller = -1;
	done = 0;
	/* look for a '$', indicating an environment variable */
	while (!done) {
	    if (filename[i] == '\0') {
		done = 1;
	    } else if (filename[i] == '$') {
		done = 1;
		doller = i;
		if (doller == 0) {
		    /* Unix-style file name */
		    for (j = 1;  filename[j] != '\0';  j++) {
			if (j >= SZ_FNAME-1) {
			    var[j] = '\0';
			    start = j + 1;
			    break;
			}
			var[j-1] = filename[j];
			if (filename[j] == '/') {
			    var[j-1] = '\0';
			    start = j + 1;
			    break;
			}
		    }
		} else {
		    /* IRAF-style file name */
		    start = doller + 1;
		    if (doller >= SZ_FNAME-1)
			doller = SZ_FNAME-1;
		    for (j = 0;  j < doller;  j++)
			var[j] = filename[j];
		    var[doller] = '\0';
		}
	    }
	    i++;
	}
	if (doller >= 0) {
	    value = getenv (var);
	    if (value == NULL) {
		strcpy (l_value, filename);
	    } else {
		strcpy (l_value, value);
		i = strlen (l_value) - 1;
		if (l_value[i] != '/')
		    strcat (l_value, "/");
		strcat (l_value, filename+start);
	    }
	    return l_value;
	} else {
	    return filename;
	}
}
