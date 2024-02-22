#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

static int  Slfit (double *, double *, double *, int, double *, int *, int,
                   double **, void (*funcs)(double, double [], int));
static int Sgaussj (double **, int, double **, int);
static void Scovsrt (double **covar, int ma, int ia[], int mfit);
static int *intVector(long nLow, long nHigh);
static void free_intVector(int *v, long nLow);
static void free_dblVector(double *v, long nLow);
static double *dblVector(long nLow, long nHigh);
static double **dmatrix(long lowRow, long highRow, long lowCol, long highCol);
static void free_dmatrix(double **m, long lowRow, long lowCol);
static void nrerror(char *);

void swap(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

/*
   Polynomial least squares fit.

   Modified/updated original routines to move away from copyright
   code.  These routines work in double precision, and handle a 
   origin-offseted data set, in order to improve accuracy.  Improved
   comments and variable naming.

   Revision history:
   ----------------
   10 Jan 00  -  Implemented (I.Busko)
   12 Jan 05  -  FitPoly now returns -1 if ndat is too small (Phil Hodge)

*/

/*
  Compute basis functions.
*/

void poly (double x, double values[], int n) {

	for (int i = 1; i <= n; i++)
	    values[i] = pow ((double)x, (double)(i-1));
}


/*
  Compute polynomial coefficients. The coefficient array must be allocated
  with two additional elements besides the ones necessary to hold the
  polynomial coefficients themselves. The first two elements in the array
  will hold the average of the input x and y arrays, respectively. The
  computed coefficients will describe a polynomial relative to these
  two averages. Function ComputePoly will properly take care of restoring
  the offsets.

  The function value (status return) will be -1 if the number of input
  data points is less than 1 + the degree of the polynomial.
*/

int FitPoly (double ix[], double iy[], double iw[], int ndat,
             int degree, double coeff[]) {

/* arguments:
double ix[];		i:  input independent variable
double iy[];		i:  input dependent variable
double iw[];		i:  input weights, always positive or zero
int ndat;		i:  number of data points
int degree;		i:  degree of polynomial to fit
double coeff[];		o:  output coefficients
*/
	int nterm, nav, *ia, status;
	double *a, *x, *y, *sig, **covar, avx, avy;

	nterm = degree + 1;

	if (ndat < nterm)
	    return -1;

	/* Compute averages. */
	avx = 0.0;
	avy = 0.0;
	nav = 0;
	for (int i = 0; i < ndat; i++) {
	    if (iw[i] > 0.0) {
	        avx += ix[i];
	        avy += iy[i];
	        nav++;
	    }
	}
	if (nav > 0) {
	    avx /= nav;
	    avy /= nav;
	} else {
	    avx = 0.0;
	    avy = 0.0;
	}

	/* Alloc memory. */
	ia    = intVector ((long)1, (long)nterm);
	a     = dblVector ((long)1, (long)nterm);
	x     = dblVector ((long)1, (long)ndat);
	y     = dblVector ((long)1, (long)ndat);
	sig   = dblVector ((long)1, (long)ndat);
	covar = dmatrix ((long)1, (long)nterm, (long)1, (long)nterm);

	/* Initialize data arrays for fitting routine. */
	for (int i = 1; i <= ndat; i++) {
	    x[i]   = ix[i-1] - avx;
	    y[i]   = iy[i-1] - avy;
	    if (iw[i-1] > 0.0)
	        sig[i] = 1.0 / iw[i-1];
	    else
	        sig[i] = 1.0;
	}
	for (int i = 1; i <= nterm; ia[i++] = 1);

	/* Solve. */
    printf("******** Calling Slfit ***********");
	if ((status = Slfit (x, y, sig, ndat, a, ia, nterm, covar, poly)))
	    return (status);

	/* Store coefficients. */
	coeff[0] = avx;
	coeff[1] = avy;
	for (int i = 1; i <= nterm; i++)
	    coeff[i+1] = a[i];

	/* Free memory. */
	free_dmatrix (covar, (long)1, (long)1);
	free_dblVector (sig, (long)1);
	free_dblVector (y, (long)1);
	free_dblVector (x, (long)1);
	free_dblVector (a, (long)1);
	free_intVector (ia, (long)1);

	return (0);
}


/*
  Compute fitted values.
*/

void ComputePoly (double ix[], int ndat, double coeff[], int degree,
                 double oy[]) {

/* arguments:
double ix[];		i:  input independent variable
int ndat;		i:  number of data points
double coeff[];		i:  coefficients
int degree;		i:  degree of polynomial
double oy[];		o:  output values
*/
	for (int i = 0; i < ndat; i++) {
	    oy[i] = 0.0;
	    for (int j = 0; j <= degree; j++)
	        oy[i] += coeff[j+2] * pow (ix[i] - coeff[0], (double)j);
	    oy[i] += coeff[1];
	}
}

/*

   The least squares fit:
   - returns an error code
   - use double precision
   - has no chisq computation
*/


static int Slfit (double xdat[], double ydat[], double sigma[], int ndat, double coeff[],
                  int coeffMask[], int ndim, double **covar,
                  void (*funcs)(double, double [], int))
{
	int mfit = 0, status;
	double ym, weight, sig2i, **beta, *basisFunc;

	beta = dmatrix((long)1, (long)ndim, (long)1, (long)1);
	basisFunc = dblVector((long)1, (long)ndim);

	for (int j = 1; j <= ndim; j++)
		if (coeffMask[j]) 
            mfit++;

	if (mfit == 0) {
	    printf ("lfit: There are no parameters to be fit.\n");
	    return (1);
	}

	for (int j = 1; j <= mfit; j++) {
		for (int k = 1; k <= mfit; k++) 
            covar[j][k] = 0.0;
		beta[j][1] = 0.0;
	}

    /* Declare this for the remainder of the routine to avoid
     * confusion for the reader.
     */
    int j, k, l, m;
	for (int i = 1; i <= ndat; i++) {
		(*funcs)(xdat[i], basisFunc, ndim);
		ym = ydat[i];
		if (mfit < ndim) {
			for (j = 1; j <= ndim; j++)
				if (!coeffMask[j]) 
                    ym -= coeff[j] * basisFunc[j];
		}
		sig2i = 1.0 / sqrt(sigma[i]);
		for (j = 0, l = 1; l <= ndim; l++) {
			if (coeffMask[l]) {
				weight = basisFunc[l] * sig2i;
				for (j++, k = 0, m = 1; m <= l; m++)
					if (coeffMask[m]) 
                        covar[j][++k] += weight * basisFunc[m];
				beta[j][1] += ym * weight;
			}
		}
	}

	for (j = 2; j <= mfit; j++)
		for (k = 1; k < j; k++)
			covar[k][j] = covar[j][k];

	if ((status = Sgaussj(covar, mfit, beta, 1)))
	    return (status);

	for (j = 0, l = 1; l <= ndim; l++)
		if (coeffMask[l]) 
            coeff[l] = beta[++j][1];

	Scovsrt (covar, ndim, coeffMask, mfit);

	free_dblVector(basisFunc, (long)1);
	free_dmatrix(beta, (long)1, (long)1);

	return (0);
}

/*
 * Redistribute the covariance values into the full covariance
 * matrix 
 */
static void Scovsrt (double **covar, int ndim, int data[], int mfit)
{
	int k = mfit;

	for (int i = mfit + 1; i <= ndim; i++) {
		for (int j = 1; j <= i; j++) {
            covar[i][j] = 0.0;
            covar[j][i] = 0.0;
        }
    }

	for (int j = ndim; j >= 1; j--) {
		if (data[j]) {
			for (int i = 1; i <= ndim; i++) 
                swap(&covar[i][k], &covar[i][j]);
			for (int i = 1; i <= ndim; i++) 
                swap(&covar[k][i], &covar[j][i]);
			k--;
		}
	}
}

/*
 * Gauss Jordan elimination with full pivot.  The full pivot means
 * there is a row exchange to move the fabs(element) to the pivotal
 * position.  In addition, there is a column exchange to maximize
 * the absolute value of the pivot.
 */
static int Sgaussj (double **a, int n, double **b, int m)
{
	int *index_col, *index_row, *index_pivot;
	int col, row;
	double maxElement, temp, pivot_element;

	index_col = intVector((long)1, (long)n);
	index_row = intVector((long)1, (long)n);
	index_pivot = intVector((long)1, (long)n);

	for (int j = 1; j <= n; j++) 
        index_pivot[j] = 0;

    /*
     * Choose the largest element as the pivot which is overall a good choice.
     */
	for (int i = 1; i <= n; i++) {
		maxElement = 0.0;
		for (int j = 1; j <= n; j++) {
			if (index_pivot[j] != 1)
				for (int k = 1; k <= n; k++) {
					if (index_pivot[k] == 0) {
						if (fabs(a[j][k]) >= maxElement) {
							maxElement = fabs(a[j][k]);
							row = j;
							col = k;
						}
					} else if (index_pivot[k] > 1) {
	                    printf
                        ("gaussj: Singular Matrix when choosing largest element.\n");
	                    return (1);
	                }
				}
        }
		++(index_pivot[col]);

        /*
         * Exchange the rows to put the pivot on the diagonal
         */
		if (row != col) {
			for (int l = 1; l <= n; l++) 
                swap(&a[row][l], &a[col][l]);
			for (int l = 1; l <= m; l++) 
                swap(&b[row][l], &b[col][l]);
		}
 
        /*
         * Divide the pivot row by the pivot element, but make
         * sure there is no divide-by-zero.
         */
		index_row[i] = row;
		index_col[i] = col;
		if (a[col][col] == 0.0) {
            printf("gaussj: Singular Matrix - pivot element is zero.\n");
	        return (1);
	    }

		pivot_element = 1.0 / a[col][col];
		a[col][col] = 1.0;

		for (int l = 1; l <= n; l++) 
            a[col][l] *= pivot_element;

		for (int l = 1; l <= m; l++) 
            b[col][l] *= pivot_element;

        /*
         * Reduce the rows, except for the pivot row.
         */
		for (int j = 1; j <= n; j++) {
			if (j != col) {
				temp = a[j][col];
				a[j][col] = 0.0;
				for (int l = 1; l <= n; l++) 
                   a[j][l] -= a[col][l] * temp;
				for (int l = 1; l <= m; l++) 
                   b[j][l] -= b[col][l] * temp;
			}
        }
	}

    /*
     * Fix up the inverse matrix by re-arranging the columns in the 
     * reverse order the initial matrix was built.
     */
	for (int l = n; l >= 1; l--) {
		if (index_row[l] != index_col[l])
			for (int k = 1; k <= n; k++)
				swap(&a[k][index_row[l]], &a[k][index_col[l]]);
	}

	free_intVector(index_pivot, 1);
	free_intVector(index_row, 1);
	free_intVector(index_col, 1);

	return (0);
}

static int *intVector(long nLow, long nHigh) {
	int *v;

	v = (int *)malloc((unsigned int) ((nHigh - nLow + 1 + 1) * sizeof(int)));
	if (!v) 
        nrerror("Memory allocation failure in intVector().");
	return v - nLow + 1;
}

static double *dblVector(long nLow, long nHigh) {
	double *v;

	v = (double *)malloc((size_t) ((nHigh - nLow + 1 + 1) * sizeof(double)));
	if (!v) 
        nrerror("Memory allocation failure in dblVector().");
	return v - nLow + 1;
}

static double **dmatrix(long lowRow, long highRow, long lowCol, long highCol) {
    long nrow = highRow-lowRow+1; 
    long ncol = highCol-lowCol+1;
	double **m;

	m = (double **) malloc((unsigned int) ((nrow + 1) * sizeof(double*)));
	if (!m) nrerror("Memory allocation failure 1 in matrix().");
	m += 1;
	m -= lowRow;

	m[lowRow]=(double *) malloc((unsigned int) ((nrow*ncol + 1) * sizeof(double)));
	if (!m[lowRow]) nrerror("Memory allocation failure 2 in matrix().");
	m[lowRow] += 1;
	m[lowRow] -= lowCol;

	for (long i = lowRow + 1; i <= highRow; i++) 
        m[i] = m[i-1] + ncol;

	return m;
}

static void free_intVector(int *v, long nLow) {
	free((char*) (v+nLow-1));
}

static void free_dblVector(double *v, long nLow) {
	free((char*) (v + nLow - 1));
}

static void free_dmatrix(double **m, long lowRow, long lowCol) {
	free((char*) (m[lowRow] + lowCol - 1));
	free((char*) (m + lowRow - 1));
}

static void nrerror(char error_text[]) {
	fprintf(stderr, "Run-time error...\n");
	fprintf(stderr, "%s\n",error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}
