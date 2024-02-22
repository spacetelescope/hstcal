#include <stdio.h>
/*
    PIKSRT  -   Sorting by straight insertion. 
    Sorts an array arrayOfFloats[0...n-1] of float type into ascending numerical 
    order of its value elements. arrayOfFloats is replaced on output by its sorted 
    rearrangement.

    IPIKSRT -   Sorting by straight insertion. 
    Same algorithm as PIKSRT, but now an arrayOfInts is along for the ride. It is
    sorted according to the ranking for the arrayOfFloats elements simultaneously.  
    There is an assumption the two input arrays are the same size.

    Order: n^2
*/

void piksrt (float arrayOfFloats[], int arrayLength) {
    float value;
    int position;

    /* Loop over the full input array */
    for (int j = 1; j < arrayLength; j++) {

        /* Get the value to be sorted and the candidate position */
        value = arrayOfFloats[j];
        position = j;

        /* Is preceding value larger than the new value? 
           If so, move the preceding value to the higher ranking position...
        */
        while (position > 0 && arrayOfFloats[position-1] > value) {
            arrayOfFloats[position] = arrayOfFloats[position-1];
            position--;
        }

        /* ...and put the new value in the lower ranking position. */
        if (position != j) {
           arrayOfFloats[position] = value;
        }
    }
}

void ipiksrt (float arrayOfFloats[], int arrayLength, int arrayOfInts[]) {
    float fValue;
    int iValue;
    int position;

    /* Loop over the input array */
    for (int j = 1; j < arrayLength; j++) {

        /* Get the value to be sorted and the candidate position */
        fValue = arrayOfFloats[j];
        iValue = arrayOfInts[j];
        position = j;

        /* Is preceding value larger than the new value? 
           If so, move the preceding value to the higher ranking position...
        */
        while (position > 0 && arrayOfFloats[position-1] > fValue) {
            arrayOfFloats[position] = arrayOfFloats[position-1];
            arrayOfInts[position] = arrayOfInts[position-1];
            position--;
        }

        /* ...and put the new value in the lower ranking position. */
        if (position != j) {
           arrayOfFloats[position] = fValue;
           arrayOfInts[position] = iValue;
        }
    }
}

