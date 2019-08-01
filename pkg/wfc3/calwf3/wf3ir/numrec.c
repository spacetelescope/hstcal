#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>

/* This file has been updated to contain only a sort routine, 
** as well as utilities for allocating and freeing arrays.
** All other routines were determined to be obsolete.
**
**	sort 
**	vector
**	free_vector
**	ivector
**	free_ivector
*/

float *vector(int nl, int nh);
void free_vector(float *v, int nl, int nh);
int *ivector(int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void randomizedQSort(float array[], int p, int r);
int randomizedPartition(float array[], int p, int r);
int partition(float array[], int p, int r);

/* This is a recursive implementation of quick sort with a randomized 
 * partition which accepts a floating point array and the number of 
 * points in the array.  The implementation is based on discussion 
 * found in "Introduction to Algorithms" by Cormen, Leiserson, and Rivest.
 * This routine replaces the former proprietary sort code.
 * 
 * In order to maintain the function prototype of the original
 * routine, the new sort() is a thin wrapper around the actual 
 * recursive sorting routine.  The wrapper only accepts two parameters, 
 * and the array is expected to be called as array-1 which is 
 * what the original sort expected. This "off by 1" and then
 * correctd in the wrapper before the recursive routine is called. 
 */ 
int sort (float array[], int n) {
    int p = 0;
    int r = n - 1;

    srand(time(NULL));
    randomizedQSort(array+1, p, r);
 
    return(0);
}

void swap(float *a, float *b) {
  float c = *a;
  *a = *b;
  *b = c;
}

/* Partition the array into two non-empty subarrays which
 * are sorted by recursive calls. Each element of the
 * array[p..pivot-1] <= each element of array[pivot+1..r].
 */
void randomizedQSort(float array[], int p, int r) {
    if (p < r) {
        int pivot = randomizedPartition(array, p, r);
        randomizedQSort(array, p, pivot - 1);
        randomizedQSort(array, pivot + 1, r);
    }
}

/* The randomized partition tries to ensure a good 
 * choice for the pivot on average.
 */
int randomizedPartition(float array[], int p, int r) {

    /* Generate a random number in the range from p...r */
    int tmp = p + (rand() % (r - p + 1));

    swap(&array[tmp], &array[r]);
    return (partition(array, p, r));
}

int partition(float array[], int p, int r) {
    int i = p - 1;
    for (int j = p; j <= r; j++) {
        if (array[j] <= array[r]) {
            swap(&array[++i], &array[j]);
        }
    }
    return(i);
}

#define NR_END 1
#define FREE_ARG char*

float *vector(int nl, int nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	return v-nl+NR_END;
}

void free_vector(float *v, int nl, int nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

int *ivector(int nl, int nh) {
/* allocate an int vector with subscript range v[nl..nh] */
        int *v;
 
	v = NULL;
        v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}
 
void free_ivector(int *v, int nl, int nh) {
/* free an int vector allocated with ivector */
        free((FREE_ARG) (v+nl-NR_END));
}
 
#undef NR_END
#undef FREE_ARG
