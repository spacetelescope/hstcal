/*
 * A. Farris: Original Implementation.
 *
 * M.D. De La Pena 06 July 1998: Changed the definition of min/maxFloat and
 * min/maxDouble to the ANSI standard of FLT_MIN/MAX and DBL_MIN/MAX as 
 * defined in <float.h>.
 *
 * M.D. De La Pena 01 September 1998: Modified the initialization of 
 * values from 0xffffffff to 0x7fffffff to eliminate warning messages when
 * compiling.
 *
 * M.D. De La Pena 12 November 1998: Removed the redefinition of FLT_MIN/MAX
 * DBL_MIN/MAX for VMS systems.
 *
 * M.D. De La Pena 29 December 1998: Since in the IRAF system an INT and a
 * LONG are the same, a special conditional has been imposed for the ALPHA 
 * running Digital Unix (DU) where maxLong = INT_MAX and minLong = INT_MIN. 
 * This was necessary as the LONG is 8 bytes (versus the INT of 4 bytes). 
 * This is only important in rare cases when input FITS keywords are at the
 * extreme possible values.
 *
 * 19 August 1999: M.D. De La Pena - Added void to functions which have no
 * parameters.
 *
 * 25 April 2002: P.E. Hodge - Removed '# include <values.h>'.
 *
 */

typedef struct {
	union {
    	    long   l;
    	    double d;
	} data;
    	int type;           /* 0 = no value, 1 is a long, 
			       2 is a double */
	int sig_digit;	    /* the number of significant digits */
    	int begpos;	    /* beginning position of value within string */
	int endpos;	    /* ending position of value within string */
	char exp;	    /* the exponent letter, if any */
			    /*  'E', 'e', 'D', or 'd' are allowed */
    	const char *errmsg; /* error message, if any */
} NumericResult;

# include <ctype.h>
# include <limits.h>
# include <float.h>

/* Initialize the static data */
static double tenpower[309] = { 1.0,
	  1.0E1,   1.0E2,   1.0E3,   1.0E4,   1.0E5,   1.0E6,   1.0E7,   1.0E8,
	  1.0E9,  1.0E10,  1.0E11,  1.0E12,  1.0E13,  1.0E14,  1.0E15,  1.0E16,
	 1.0E17,  1.0E18,  1.0E19,  1.0E20,  1.0E21,  1.0E22,  1.0E23,  1.0E24,
	 1.0E25,  1.0E26,  1.0E27,  1.0E28,  1.0E29,  1.0E30,  1.0E31,  1.0E32,
	 1.0E33,  1.0E34,  1.0E35,  1.0E36,  1.0E37,  1.0E38,  1.0E39,  1.0E40,
	 1.0E41,  1.0E42,  1.0E43,  1.0E44,  1.0E45,  1.0E46,  1.0E47,  1.0E48,
	 1.0E49,  1.0E50,  1.0E51,  1.0E52,  1.0E53,  1.0E54,  1.0E55,  1.0E56,
	 1.0E57,  1.0E58,  1.0E59,  1.0E60,  1.0E61,  1.0E62,  1.0E63,  1.0E64,
	 1.0E65,  1.0E66,  1.0E67,  1.0E68,  1.0E69,  1.0E70,  1.0E71,  1.0E72,
	 1.0E73,  1.0E74,  1.0E75,  1.0E76,  1.0E77,  1.0E78,  1.0E79,  1.0E80,
	 1.0E81,  1.0E82,  1.0E83,  1.0E84,  1.0E85,  1.0E86,  1.0E87,  1.0E88,
	 1.0E89,  1.0E90,  1.0E91,  1.0E92,  1.0E93,  1.0E94,  1.0E95,  1.0E96,
	 1.0E97,  1.0E98,  1.0E99, 1.0E100, 1.0E101, 1.0E102, 1.0E103, 1.0E104,
	1.0E105, 1.0E106, 1.0E107, 1.0E108, 1.0E109, 1.0E110, 1.0E111, 1.0E112,
	1.0E113, 1.0E114, 1.0E115, 1.0E116, 1.0E117, 1.0E118, 1.0E119, 1.0E120,
	1.0E121, 1.0E122, 1.0E123, 1.0E124, 1.0E125, 1.0E126, 1.0E127, 1.0E128, 
	1.0E129, 1.0E130, 1.0E131, 1.0E132, 1.0E133, 1.0E134, 1.0E135, 1.0E136,
	1.0E137, 1.0E138, 1.0E139, 1.0E140, 1.0E141, 1.0E142, 1.0E143, 1.0E144,
	1.0E145, 1.0E146, 1.0E147, 1.0E148, 1.0E149, 1.0E150, 1.0E151, 1.0E152, 
	1.0E153, 1.0E154, 1.0E155, 1.0E156, 1.0E157, 1.0E158, 1.0E159, 1.0E160, 
	1.0E161, 1.0E162, 1.0E163, 1.0E164, 1.0E165, 1.0E166, 1.0E167, 1.0E168,
	1.0E169, 1.0E170, 1.0E171, 1.0E172, 1.0E173, 1.0E174, 1.0E175, 1.0E176, 
	1.0E177, 1.0E178, 1.0E179, 1.0E180, 1.0E181, 1.0E182, 1.0E183, 1.0E184, 
	1.0E185, 1.0E186, 1.0E187, 1.0E188, 1.0E189, 1.0E190, 1.0E191, 1.0E192,
	1.0E193, 1.0E194, 1.0E195, 1.0E196, 1.0E197, 1.0E198, 1.0E199, 1.0E200,
	1.0E201, 1.0E202, 1.0E203, 1.0E204, 1.0E205, 1.0E206, 1.0E207, 1.0E208, 
	1.0E209, 1.0E210, 1.0E211, 1.0E212, 1.0E213, 1.0E214, 1.0E215, 1.0E216,
	1.0E217, 1.0E218, 1.0E219, 1.0E220, 1.0E221, 1.0E222, 1.0E223, 1.0E224,
	1.0E225, 1.0E226, 1.0E227, 1.0E228, 1.0E229, 1.0E230, 1.0E231, 1.0E232,
	1.0E233, 1.0E234, 1.0E235, 1.0E236, 1.0E237, 1.0E238, 1.0E239, 1.0E240,
	1.0E241, 1.0E242, 1.0E243, 1.0E244, 1.0E245, 1.0E246, 1.0E247, 1.0E248, 
	1.0E249, 1.0E250, 1.0E251, 1.0E252, 1.0E253, 1.0E254, 1.0E255, 1.0E256, 
	1.0E257, 1.0E258, 1.0E259, 1.0E260, 1.0E261, 1.0E262, 1.0E263, 1.0E264,
	1.0E265, 1.0E266, 1.0E267, 1.0E268, 1.0E269, 1.0E270, 1.0E271, 1.0E272,
	1.0E273, 1.0E274, 1.0E275, 1.0E276, 1.0E277, 1.0E278, 1.0E279, 1.0E280,
	1.0E281, 1.0E282, 1.0E283, 1.0E284, 1.0E285, 1.0E286, 1.0E287, 1.0E288, 
	1.0E289, 1.0E290, 1.0E291, 1.0E292, 1.0E293, 1.0E294, 1.0E295, 1.0E296,
	1.0E297, 1.0E298, 1.0E299, 1.0E300, 1.0E301, 1.0E302, 1.0E303, 1.0E304,
	1.0E305, 1.0E306, 1.0E307, 1.0E308 };
static const int minfltexp = -38;
static const int maxfltexp = 38;
static const int mindblexp = -308;
static const int maxdblexp = 308;
static const int maxsigdigits = 17;
static const int maxdigl = 9; /* max digits in a long */
static const int maxexpdig = 3; /* max digits in an exponent */

# if defined(__cplusplus)
extern "C" {
# endif

# if SIZEOF_INT == 4
     const int minLong = INT_MIN;
     const int maxLong = INT_MAX;
# else
     const int minLong = LONG_MIN;
     const int maxLong = LONG_MAX;
# endif
const float minFloat = FLT_MIN;
const float maxFloat = FLT_MAX;
const double minDouble = DBL_MIN;
const double maxDouble = DBL_MAX;

float FloatNAN(void);
double DoubleNAN(void);
void get_numeric(const char *s, int len, NumericResult *result);

# if defined(__cplusplus)
}
# endif


/* Return NAN values for floats and doubles */
float FloatNAN(void) {
	static long l_nanFloat = 0x7fffffff;
	static float *f_nanFloat = (float *)&l_nanFloat;
	return *f_nanFloat;
}
double DoubleNAN(void) {
	typedef struct { double x; long l1; long l2; } L_nanDouble;
	static L_nanDouble l_nanDouble = { 0.0, 0x7fffffff, 0x7fffffff }; 
	static double *d_nanDouble = (double *)&l_nanDouble.l1;
	return *d_nanDouble;
}

#define digit2bin(c) ((c) - '0')
#define ten(numb,pow) \
(((pow) > 0) ? (((double)(numb)) * tenpower[(pow)]) : \
(((pow) < 0) ? (((double)(numb)) / tenpower[-(pow)]) : ((double)(numb))))


static int ckaccum(double *d, int numb, int pow) {
	/* compute d += numb * 10**pow checking for over/underflow */
	double tmp = (double)numb;
	if (pow > 0) {
	    if (pow > 2 * maxdblexp)
		return 1;
	    if (tmp <= (maxDouble / tenpower[pow]))
		tmp *= tenpower[pow];
	    else
		return 1;
	} else if (pow < 0) {
	    pow = -pow;
	    if (pow > 2 * maxdblexp)
		return -1;
	    if ((minDouble * tenpower[pow]) <= tmp)
		tmp /= tenpower[pow];
	    else
		return -1;
	}
	if ((maxDouble - tmp) < *d)
		return 1;
	*d += tmp;
	return 0;
}

#define NOVALUE 0
#define LONG 1
#define DOUBLE 2

/*
**
	This is the algorithm used to convert ASCII to binary.

	1. Skip any spaces at the beginning
	2. Mark position as the start of the field
	3. Check for any sign that may be present
	4. Skip any leading 0s as insignificant
	5. Get integer part of number.  Get digits, store in longs, 
		and count significant digits 
	6. If the end of the string has been reached, or if the ending
		char is NOT one of `.', `E', `e', `D', `d'
		the number is an integer and we are done.
	Otherwise,
	7. If valid, the number is a floating point number.
           If there is a fraction part,
	       8. Get the fraction part 
	       9. Check next char after fraction part -- if it is NOT
		    one of `E', `e', `D', `d', then goto #14.
	   endif
	10. If there is an exponent,
	      11. Check for any sign that may be present 
	      12. Skip any leading 0s as insignificant 
	      13. Get, at most, the max digits in an exponent 
	    endif
	14. Compute real value 
**
*/

void get_numeric(const char *s, int len, NumericResult *result) {
	int n;			/* the number of chars processed -- */
				/*   n == len signals `at end' */
	int i, j;		/* counter, confined to local context */
	const char *p;		/* utility variable, confined local context */

	int negsign    = 0;	
	long intpart1  = 0;	/* part 1 of digits of integer part */
	long intpart2  = 0;	/* part 2 of digits of integer part */
	int  sigint    = 0;	/* number of significant digits */

	int pointpos   = 0;	/* the position of the decimal point */

	long fracpart1 = 0;	/* part 1 of digits of fraction part */
	long fracpart2 = 0;	/* part 2 of digits of fraction part */
	int  sigfrac   = 0;	/* number of significant digits in fraction */
	int  fracpos   = 0;	/* position of first digit relative to point */
	int  exp       = 0;	/* exponent */
	int  sigexp    = 0;	/* number of significant digits in exponent */
	char exp_type  = ' ';	/* the exponent letter */

	result->type   	  = NOVALUE; /* Initialize result */
	result->sig_digit = 0;
	result->errmsg    = 0;
	result->data.l    = 0;
	result->begpos    = 0;
	result->endpos    = 0;
	result->exp       = '\0';

	/* 1. Skip any spaces at the beginning */
	for (n = 0; (*s == ' ' || *s == '\t') && (n < len); ++n, ++s) ;
	if (n == len) {
	    result->errmsg = "Value field is all blanks";
	    return;
	}

	/* 2. Mark position as the start of the field */
	result->begpos = n;

	/* 3. Check for any sign that may be present */
	negsign = (*s == '-' ? (++s, ++n, 1) :
			(*s == '+' ? (++s, ++n, 0) : 0));
	if (n == len || !(isdigit(*s) || *s == '.')) {
	    result->errmsg = "Not a number";
	    return;
	}

	/* 4. Skip any leading 0s as insignificant */
	for (; *s == '0' && n < len; ++n, ++s) ;
	if (n == len) {
	    result->type = LONG;
	    result->data.l = 0;
	    result->endpos = n - 1;
	    return;
	}

	/* 5. Get integer part of number.  Get digits, store in longs, */
	/*	and count significant digits */
	if (isdigit(*s)) {
	    intpart1 = digit2bin(*s);
	    ++sigint;
	    for (++n, ++s; isdigit(*s) && (n < len); ++n, ++s) {
		++sigint;
		intpart1 = intpart1 * 10 + digit2bin(*s);
		if (sigint == maxdigl) /* Longs can always hold maxdigl digits */
		    break;
	    }
	    if (sigint == maxdigl && n < len) {
		++n;
		++s;
		if (isdigit(*s) && n < len) {
		    intpart2 = digit2bin(*s);
		    ++sigint;
		    for (++n, ++s; isdigit(*s) && (n < len); ++n, ++s) {
			++sigint;
			intpart2 = intpart2 * 10 + digit2bin(*s);
			if (sigint == maxsigdigits)
			    break;
		    }
		}
		if (sigint == maxsigdigits && n < len) { /* Discard digits */
		    for (++n, ++s; isdigit(*s) && (n < len); ++n, ++s)
			++sigint; /* But count them */
		}
	    }
	}

	/* 6. If the end of the string has been reached, or if the ending */
	/*	char is NOT one of `.', `E', `e', `D', `d'		  */
	/*	the number is an integer and we are done.		  */
	if (n == len || (!(*s == '.' || *s == 'E' || *s == 'e' || 
					*s == 'D' || *s == 'd'))) {
	    result->endpos = n - 1;
	    result->type = LONG;
	    result->sig_digit = sigint;
	    if (sigint < 10) {
		result->data.l = negsign ? (-intpart1) : intpart1;
		return;
	    }
	    if (sigint > 10) {
		if (negsign) {
		    result->errmsg = "Integer underflow";
		    result->data.l = minLong;
		} else {
		    result->errmsg = "Integer overflow";
		    result->data.l = maxLong;
		}
		return;
	    }
	    if (negsign) {
		intpart1 = -intpart1;
		if (((minLong + intpart2) / 10) <= intpart1) {
		    result->data.l = intpart1 * 10 - intpart2;
		    return;
		} else {
		    result->errmsg = "Integer underflow";
		    result->data.l = minLong;
		    return;
		}
	    } else {
		if (((maxLong - intpart2) / 10) >= intpart1) {
		    result->data.l = intpart1 * 10 + intpart2;
		    return;
		} else {
		    result->errmsg = "Integer overflow";
		    result->data.l = maxLong;
		    return;
		}
	    }
	}

	/* 7. If valid, the number is float or double.  Get the fraction */
	/*	part, if any. */
	if (*s == '.') {
	    pointpos = n;

	    /* 8. Get the fraction part */
	    ++s;
	    ++n;
	    if (n == len)
		goto real_value;
	    if (*s == 'E' || *s == 'e' || *s == 'D' || *s == 'd')
		goto get_exp;
	    if (!(isdigit(*s)))
		goto real_value;
	    if (sigint == 0) { /* i. e., there are no sig digs in int part */
		for (; *s == '0' && n < len; ++n, ++s) ;
		if (n == len)
		    goto real_value;
		fracpos = n - pointpos - 1;
	    }
	    p = s; /* save start of fraction part */
	    for (; isdigit(*s) && (n < len); ++s, ++n, ++sigfrac) ;
	    --sigfrac;
	    while (p[sigfrac] == '0' && sigfrac >= 0) /* eliminate trailing 0s */
		--sigfrac;
	    ++sigfrac;
	    /* adjust sigfrac to only process the maximum number of sig digits */
	    if (sigint >= maxsigdigits)
		sigfrac = 0;
	    else if ((sigint + sigfrac) > maxsigdigits)
		    sigfrac = maxsigdigits - sigint;
	    /* get sigfrac digits starting at p */
	    if (sigfrac > 0) {
		fracpart1 = digit2bin(*p);
		for (i = 1, ++p; i < sigfrac && i < maxdigl; ++i, ++p)
		    fracpart1 = fracpart1 * 10 + digit2bin(*p);
		if (i == maxdigl && i < sigfrac) {
		    if (i < sigfrac) {
			fracpart2 = digit2bin(*p);
			for (++i, ++p; i < sigfrac; ++i, ++p)
			    fracpart2 = fracpart2 * 10 + digit2bin(*p);
		    }
		}
	    }

	    /* 9. Check next char after fraction part */
	    if (!(*s == 'E' || *s == 'e' || *s == 'D' || *s == 'd'))
		goto real_value;
	}

	/* 10. Get the exponent, if any */
    get_exp:
	exp_type = *s++;
	++n;
	if (n == len)
	    goto real_value;

	/* 11. Check for any sign that may be present */
	i = 0; /* i is the sign */
	if (*s == '+' || *s == '-') {
	    if (*s == '-')
		i = 1;
	    ++n;
	    ++s;
	    if (n == len) {	   /* This is a strange condition, the field */
		result->endpos = n; /* ends with a sign.  Mark it as an error. */
		result->errmsg = "Value field is not a valid number";
		return;
	    }
	}

	/* 12. Skip any leading 0s as insignificant */
	for (; *s == '0' && n < len; ++n, ++s) ;
	if (n == len)
	    goto real_value;

	/* 13. Get, at most, the max digits in an exponent */
	if (isdigit(*s)) {
	    exp = digit2bin(*s);
	    ++sigexp;
	    for (++n, ++s; isdigit(*s) && (n < len); ++n, ++s) {
		++sigexp;
		exp = exp * 10 + digit2bin(*s);
		if (sigexp == maxexpdig) {
		    ++n;
		    if (n == len)
			break;
		    if (isdigit(s[1])) {
			result->endpos = n - 1;
			if (i)
			    result->errmsg = "Exponent underflow";
			else
			    result->errmsg = "Exponent overflow";
			result->type = DOUBLE;
			result->data.d = DoubleNAN();
			return;
		    }
		    break;
		}
	    }
	    if (i)
		exp = -exp;
	}

    real_value: /* 14. Compute real value */
	/* It's real because either there is a point or exp_type != ' '. */
	result->endpos = n - 1;
	result->sig_digit = sigint + sigfrac;
	result->exp = exp_type;
	result->type = DOUBLE;

	if (intpart1 == 0) {
	    if (fracpart1 == 0) {
	        result->data.d = 0.0;
	        return;
	    } else if ((exp - fracpos - 1) > mindblexp &&
	    	(exp - fracpos) < maxdblexp) {
	        result->data.d = (fracpart2 == 0) ? 0.0 :
	    	   ten(fracpart2,(exp - sigfrac - fracpos));
	        i = (sigfrac < maxdigl) ? sigfrac : maxdigl;
	        i = exp - fracpos - i;
	        result->data.d += ten(fracpart1,i);
	        if (negsign)
	    	  result->data.d = -result->data.d;
	        return;
	    }
	} else if ((sigint - 1 + exp) > mindblexp &&
	    	(sigint + exp) < maxdblexp) {
	    result->data.d = (fracpart2 == 0) ? 0.0 :
	           ten(fracpart2,(exp - sigfrac - fracpos));
	    if (fracpart1) {
	        i = (sigfrac < maxdigl) ? sigfrac : maxdigl;
	        i = exp - fracpos - i;
	        result->data.d += ten(fracpart1,i);
	    }
	    if (intpart2) {
	        i = sigint - maxsigdigits;
	        if (i < 0)
	    	i = 0;
	        i += exp;
	        result->data.d += ten(intpart2,i);
	    }
	    i = sigint - maxdigl;
	    if ( i < 0)
	        i = 0;
	    i += exp;
	    result->data.d += ten(intpart1,i);
	    if (negsign)
	        result->data.d = -result->data.d;
	    return;
	}

	/* There might be an overflow. */
	result->data.d = 0.0;
	j = ckaccum(&result->data.d,fracpart2,(exp - sigfrac - fracpos));
	if (j) {
	    result->errmsg = (j > 0) ? "Double overflow" : "Double underflow";
	    result->data.d = DoubleNAN();
	    return;
	}
	if (fracpart1) {
	    i = (sigfrac < maxdigl) ? sigfrac : maxdigl;
	    i = exp - fracpos - i;
	    j = ckaccum(&result->data.d,fracpart1,i);
	    if (j) {
	      result->errmsg = (j > 0) ? "Double overflow" : "Double underflow";
	      result->data.d = DoubleNAN();
	      return;
	    }
	}
	if (intpart2) {
	    i = sigint - maxsigdigits;
	    if (i < 0)
		i = 0;
	    i += exp;
	    j = ckaccum(&result->data.d,intpart2,i);
	    if (j) {
	      result->errmsg = (j > 0) ? "Double overflow" : "Double underflow";
	      result->data.d = DoubleNAN();
	      return;
	    }
	}
	i = sigint - maxdigl;
	if ( i < 0)
	    i = 0;
	i += exp;
	j = ckaccum(&result->data.d,intpart1,i);
	if (j) {
	    result->errmsg = (j > 0) ? "Double overflow" : "Double underflow";
	    result->data.d = DoubleNAN();
	    return;
	}
	return;
}


