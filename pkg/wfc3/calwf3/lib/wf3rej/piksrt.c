/*
    PIKSRT  -   Sorting by straight insertion. 

    Sorts an array arr[0...n-1] of float type into ascending numerical 
    order of its .value elements. arr is replaced on output by its sorted 
    rearrangement.

    Order: n^2

    From "Numerical Recipes - The Art of Scientific Computing",
    Press, W.H., Teukolsky, S.A., Vetterling, W.T. and Flannery, B.P., 
    2nd edition, Cambridge University Press, 1995.  

*/

/*
void piksrt (float arr[], int n) {

    int         i, j;
    float       a;

    for (j = 1; j < n; j++) {
        a = arr[j];
        i = j - 1;
        while ((i >= 0) && (arr[i] > a)) {
            arr[i+1] = arr[i];
            i--;
        }
        arr[i+1] = a;
    }
}
*/

void ipiksrt (float arr[], int n, int brr[]) {

    int         i, j;
    float       a;
    int		b;

    for (j = 1; j < n; j++) {
        a = arr[j];
	b = brr[j];
        i = j - 1;
        while ((i >= 0) && (arr[i] > a)) {
            arr[i+1] = arr[i];
	    brr[i+1] = brr[i];
            i--;
        }
        arr[i+1] = a;
        brr[i+1] = b;
    }
}

