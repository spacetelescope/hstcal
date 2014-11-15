# if !defined(NUMERIC_LOADED)
# define NUMERIC_LOADED

/*
** The routine get_numeric converts numeric data in the form of
** ASCII strings to binary form.  Results are returned in the structure
** NumericResult.  Input is in the form of a character string and the
** number of characters to search within the string.  (Therefore, the 
** string need not be null terminated.)  The numeric result may be either
** type long or double.
**
** The strings recognized as numeric conform to the usual way of
** writing scientific numbers and NOT to the ANSI C standard.  Leading
** zeros are discarded as insignificant and the exponent letter can be
** either 'E', 'e', 'D', or 'd'.
**
** This module also defines a number of useful constants:
**	minLong, maxLong, minFloat, maxFloat, minDouble, maxDouble,
** The functions:
**	FloatNAN() and DoubleNAN()
** return IEEE values for Not-A-Number.  They must be used with care.
**
** The NumericResult structure returns the data and its type, as well as
** other potentially useful information: the number of significant digits,
** the beginning and ending position of the value within the string and the
** exponent letter. 
**
** If there is no identifiable numeric value, the returned type is 0 (no 
** value) and the returned data has no significance.  If the value is a 
** floating point number it is returned as a double.  
**
** If an error occurrs, errmsg returns a text string error message.
** If no error occurs, errmsg is null.  If a value was intended to be 
** an integer, but would have produced an overflow or underflow, the
** type is marked as 1 (long) and maxLong or minLong is returned.  If 
** a floating point value would have produced an overflow or underflow, 
** the type is marked as 2 (double) and  DoubleNAN() is returned.  In 
** the case of overflows or underflows, errmsg is still set.
**
** 19 August 1999: M.D. De La Pena - Added void to functions which have no
** parameters.
**
*/

typedef struct {
	union {
    	    long   l;
    	    double d;
	} data;
    	int type;           /* 0 = no value, 1 is a long, 2 is a double */
	int sig_digit;	    /* the number of significant digits */
    	int begpos;	    /* beginning position of value within string */
	int endpos;	    /* ending position of value within string */
	char exp;	    /* the exponent letter, if any */
			    /*  'E', 'e', 'D', or 'd' are allowed */
    	const char *errmsg; /* error message, if any */
} NumericResult;

# if defined(__cplusplus)
extern "C" {
# endif

extern const int minLong;
extern const int maxLong;
extern const float minFloat;
extern const float maxFloat;
extern const double minDouble;
extern const double maxDouble;
float FloatNAN(void);
double DoubleNAN(void);

void get_numeric(const char *s, int len, NumericResult *result);

# if defined(__cplusplus)
}
# endif


# endif

