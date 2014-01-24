# include "stis.h"

/*

   Tests the printtime flag and calls TimeStamp.


   Revision history:
   ----------------

   11 Apr 97  -  Implemented after code review (I. Busko)

*/

void TimeStamp6 (int prittime, char* msg, char *rootname) {

	if (prittime)
	    TimeStamp (msg, rootname);
} 
