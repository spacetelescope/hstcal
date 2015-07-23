#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

static int  Slfit (double *, double *, double *, int, double *, int *, int,
                   double **, void (*funcs)(double, double [], int));
static int Sgaussj (double **, int, double **, int);
static void Scovsrt (double **covar, int ma, int ia[], int mfit);
static int *ivector(long nl, long nh);
static void free_ivector(int *v, long nl, long nh);
static void free_dvector(double *v, long nl, long nh);
static double *dvector(long nl, long nh);
static double **dmatrix(long nrl, long nrh, long ncl, long nch);
static void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
static void nrerror(char *);


/*
   Polynomial least squares fit.

   We use for now Numerical Recipes routines modified to work in
   double precision, and to handle a origin-offseted data set, in
   order to improve accuracy.



   Revision history:
   ----------------
   10 Jan 00  -  Implemented (I.Busko)
   12 Jan 05  -  FitPoly now returns -1 if ndat is too small (Phil Hodge)

*/

/*
  Compute basis functions.
*/

void poly (double x, double values[], int n) {

	int i;

	for (i = 1; i <= n; i++)
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
	int i, nterm, nav, *ia, status;
	double *a, *x, *y, *sig, **covar, avx, avy;

	nterm = degree + 1;

	if (ndat < nterm)
	    return -1;

	/* Compute averages. */
	avx = 0.0;
	avy = 0.0;
	nav = 0;
	for (i = 0; i < ndat; i++) {
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
	ia    = ivector (1, nterm);
	a     = dvector (1, nterm);
	x     = dvector (1, ndat);
	y     = dvector (1, ndat);
	sig   = dvector (1, ndat);
	covar = dmatrix (1, nterm, 1, nterm);

	/* Initialize data arrays for fitting routine. */
	for (i = 1; i <= ndat; i++) {
	    x[i]   = ix[i-1] - avx;
	    y[i]   = iy[i-1] - avy;
	    if (iw[i-1] > 0.0)
	        sig[i] = 1.0 / iw[i-1];
	    else
	        sig[i] = 1.0;
	}
	for (i = 1; i <= nterm; ia[i++] = 1);

	/* Solve. */
	if ((status = Slfit (x, y, sig, ndat, a, ia, nterm, covar, poly)))
	    return (status);

	/* Store coefficients. */
	coeff[0] = avx;
	coeff[1] = avy;
	for (i = 1; i <= nterm; i++)
	    coeff[i+1] = a[i];

	/* Free memory. */
	free_dmatrix (covar, 1, nterm, 1, nterm);
	free_dvector (sig,   1, ndat);
	free_dvector (y,     1, ndat);
	free_dvector (x,     1, ndat);
	free_dvector (a,     1, nterm);
	free_ivector (ia,    1, nterm);

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
	int i, j;

	for (i = 0; i < ndat; i++) {
	    oy[i] = 0.0;
	    for (j = 0; j <= degree; j++)
	        oy[i] += coeff[j+2] * pow (ix[i] - coeff[0], (double)j);
	    oy[i] += coeff[1];
	}
}


/*
   Numerical Recipes code.

   This code was modified from its original form to:

   - return an error code instead of aborting plain and simple.
   - use double precision throughout.
   - chisq computation removed.
*/


static int Slfit (double x[], double y[], double sig[], int ndat, double a[],
                  int ia[], int ma, double **covar,
                  void (*funcs)(double, double [], int))
{
	int i,j,k,l,m,mfit=0,status;
	double ym,wt,sig2i,**beta,*afunc;

	beta=dmatrix(1,ma,1,1);
	afunc=dvector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	if (mfit == 0) {
	    printf ("lfit: no parameters to be fitted\n");
	    return (1);
	}
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}
	for (i=1;i<=ndat;i++) {
		(*funcs)(x[i],afunc,ma);
		ym=y[i];
		if (mfit < ma) {
			for (j=1;j<=ma;j++)
				if (!ia[j]) ym -= a[j]*afunc[j];
		}
		sig2i=1.0/sqrt(sig[i]);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=afunc[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) covar[j][++k] += wt*afunc[m];
				beta[j][1] += ym*wt;
			}
		}
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++)
			covar[k][j]=covar[j][k];
	if ((status = Sgaussj(covar,mfit,beta,1)))
	    return (status);
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) a[l]=beta[++j][1];
	Scovsrt (covar,ma,ia,mfit);
	free_dvector(afunc,1,ma);
	free_dmatrix(beta,1,ma,1,1);
	return (0);
}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

static void Scovsrt (double **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	double swap;

	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) {
		if (ia[j]) {
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}
#undef SWAP

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

static int Sgaussj (double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) {
	                                    printf
                                           ("gaussj: Singular Matrix-1\n");
	                                    return (1);
	                                }
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) {
                    printf("gaussj: Singular Matrix-2\n");
	            return (1);
	        }
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
	return (0);
}
#undef SWAP

static int *ivector(nl,nh)
long nh,nl;
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+1)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+1;
}

static double *dvector(long nl, long nh)
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+1;
}

static double **dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	m=(double **) malloc((unsigned int)((nrow+1)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += 1;
	m -= nrl;

	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+1)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

static void free_ivector(int *v, long nl, long nh)
{
	free((char*) (v+nl-1));
}

static void free_dvector(v,nl,nh)
double *v;
long nh,nl;
{
	free((char*) (v+nl-1));
}

static void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
	free((char*) (m[nrl]+ncl-1));
	free((char*) (m+nrl-1));
}

static void nrerror(char error_text[])
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
