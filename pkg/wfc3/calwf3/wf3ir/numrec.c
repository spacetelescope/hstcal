#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#define INSERT_SORT_LIMIT 7
#define AUX_STORAGE 50

/* This file has been updated to contain only a sort routine, 
** as well as utilities for allocating and freeing arrays.
** All other routines were determined to be obsolete.
*/

/* Allocate an int vector with subscript range data[nlow..nhigh] */
int *ivector(int nlow, int nhigh) {
    int *data = NULL;
 
    data = (int *)calloc((size_t)(nhigh - nlow + 1 + 1), (size_t)sizeof(int));
    return data - nlow + 1;
}
 
/* Free an int vector allocated with ivector */
void free_ivector(int *data, int nlow) {
    free((char *) (data + nlow - 1));
}

/* Simple swap routine */
void swap(float *a, float *b) {
  float c = *a;
  *a = *b;
  *b = c;
}

/*
 * Sort array of length n (1-n) in ascending order.
 * ATTN: Note the beginning index of "1" and not "0".
 */
int sort(float array[], int n) {

    int l = 1;
	int m = n;
	int jstack = 0, *auxStack;
	float aDataValue;

    /*
     * Quick sort needs some auxiliary storage as a stack data structure
     */
	auxStack = ivector(1, AUX_STORAGE);
	if (auxStack == NULL) 
        return (1);

    /*
     * Continuous loop
     */
	for (;;) {
        /* 
         * When there are few enough data values, switch to insertion sort.
         */
		if (m - l < INSERT_SORT_LIMIT) {
            int i;
			for (int j = l + 1; j <= m; j++) {
				aDataValue = array[j];
				for (i = j - 1; i >= l; i--) {
					if (array[i] <= aDataValue) 
                        break;
					array[i + 1] = array[i];
				}
				array[i + 1] = aDataValue;
			}
			if (jstack == 0) 
                break;

            /* 
             * Pop the data off the stack in preparation for the new partitioning cycle
             */
			m = auxStack[jstack--];
			l = auxStack[jstack--];

        /*
         * quick sort
         */
		} else {
            /* Median of left, center, and right elements is the pivot */
            int k = (l + m) / 2;
			swap(&array[k], &array[l + 1]);
			if (array[l] > array[m]) {
				swap(&array[l], &array[m]);
			}
			if (array[l + 1] > array[m]) {
				swap(&array[l + 1], &array[m]);
			}
			if (array[l] > array[l + 1]) {
				swap(&array[l], &array[l + 1]);
			}

            /* 
             * Setup the indices for the partitioning
             */
			int i = l + 1;
			int j = m;
			aDataValue = array[l + 1];
			for (;;) {
                /* Find element > aDataValue */
				do {
                    i++; 
                } while (array[i] < aDataValue);

                /* Find element < aDataValue */
				do {
                    j--; 
                } while (array[j] > aDataValue);

                /* The search indices have crossed, so the partition is done. */
				if (j < i) 
                    break;

                /* Swap the elements so they are in the proper subarray according to their magnitude. */
				swap(&array[i], &array[j]);
			}
			array[l + 1] = array[j];
			array[j] = aDataValue;
			jstack += 2;
			if (jstack > AUX_STORAGE) {
			    return (1);
			}
			if (m - i + 1 >= j - l) {
				auxStack[jstack] = m;
				auxStack[jstack - 1] = i;
				m = j - 1;
			} else {
				auxStack[jstack] = j - 1;
				auxStack[jstack - 1] = l;
				l = i;
			}
		}
	}

	free_ivector(auxStack, 1);
	return (0);
}

