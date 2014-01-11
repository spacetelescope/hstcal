# include "stis.h"

/* Initialize the elements of a RefTab structure.

   Phil Hodge, 2004 Dec 27:
        Duplicate code extracted from several functions to this one.
*/

void InitRefTab (RefTab *table) {

	table->name[0]      = '\0';
	table->pedigree[0]  = '\0';
	table->descrip[0]   = '\0';
	table->descrip2[0]  = '\0';
	table->exists       = EXISTS_UNKNOWN;
	table->goodPedigree = PEDIGREE_UNKNOWN;
}
