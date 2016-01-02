# include <stdio.h>
# include <stdlib.h>
# include <ctype.h>

/*  strtor -- convert a string of real numbers into a real array 

  Description:
  ------------
  If the input string is blank(s), this routine will return 0.  If there are
  characters other than digits, decimal point, comma, semi-colon, slash or
  blank in the input string, this routine will print an error message and
  return -1.
  A null (or blank) string is valid, and the function value will be zero
  (i.e. no value specified).  However, a null substring (e.g. repeated
  separators) is an error.

  Date		Author		Description
  ----		------		-----------
  09-May-1996  J.-C. Hsu	Adapt from the SPP code strtor.x
  10-Feb-2000  Phil Hodge	Replace exit (2) with return (0).
  20-Mar-2002  Ivo Busko	Move to library, return -1 on error;
				support multiple separators.
  24-Jan-2005  Phil Hodge	Fix indexing error.
*/
int strtor (char *str, float arr[])
{
	/* indexes in the string; the substring to be copied to tmp (from
	   which a numerical value will be read) runs from ip0 to ip
	*/
	int	ip0, ip;
	int	n, i;
	int	done;
	char	tmp[100];

	n = 0;
	ip0 = 0;
	ip = 0;

	while (str[ip0] == ' ') {
	    ip0++;
	    ip++;
        }
	if (str[ip0] == '\0')
	    return (0);

	done = 0;
	while (!done) {
	    if (str[ip] == ',' || str[ip] == ';' || str[ip] == '/' || 
		str[ip] == ' ' || str[ip] == '\0') {
		for (i = 0; i < (ip-ip0); ++i)
	    	    tmp[i] = str[ip0+i];
		tmp[ip-ip0] = '\0';
		if (!(isdigit (tmp[0]) || tmp[0] == '-' || tmp[0] == '.')) {
		    printf ("illegal input string '%s'\n", str);
		    return (-1);
		}
	    	arr[n] = (float) atof(tmp);
		++n;
	        if (str[ip] == '\0')
		    break;
		ip++;			/* increment past the separator */
		ip0 = ip;

		while (str[ip0] == ' ') {
		    ip0++;
		    ip++;
	        }
	    } else {
		++ip;
	    }
	    if (str[ip0] == '\0')
		break;
	}

	return (n);
}
