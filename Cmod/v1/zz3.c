# include <stdlib.h>
# include "ctables.h"

# define INCREMENT 100
static TblInfoPtr **dataIRAFPointer=NULL;
static int nelem=0;	/* current length of dataIRAFPointer array */

/* public functions (prototypes are in ctables.h):

    TableDescr *tp;
    ColumnDescr *cp;
    TblInfoPtr *tptr;
    IRAFPointer n;

    tptr = makeTableDescr (tp);
    tptr = makeColumnDescr (cp);
    n = saveIrafP (tptr);
    tp = getTableDescr (n);
    cp = getColumnDescr (n);
*/

static IRAFPointer findEmptyElement (void);
static int reallocBuffer (int increment);

/*
saveIrafP        save in list
getTableDescr    get from list
getColumnDescr   get from list
makeTableDescr   allocate table descriptor
makeColumnDescr  allocate column descriptor
reallocBuffer    reallocate
xxx              deallocate (delete from list)
*/

TblInfoPtr *makeTableDescr (TableDescr *tp) {
	TblInfoPtr *tptr;
	tptr = (TblInfoPtr *)calloc (1, sizeof(TblInfoPtr));
	tptr->ptype = TABLE_DESCR;
	tptr->tp = tp;
	tptr->cp = NULL;
	return tptr;
}

TblInfoPtr *makeColumnDescr (ColumnDescr *cp) {
	TblInfoPtr *tptr;
	tptr = (TblInfoPtr *)calloc (1, sizeof(TblInfoPtr));
	tptr->ptype = COLUMN_DESCR;
	tptr->tp = NULL;
	tptr->cp = cp;
	return tptr;
}

IRAFPointer saveIrafP (TblInfoPtr *tptr) {

	IRAFPointer n;

	n = findEmptyElement();
	if (n >= 0)
	    dataIRAFPointer[n] = tptr;

	return n;
}

TableDescr *getTableDescr (IRAFPointer n) {
	/* should check whether the type is correct */
	return (dataIRAFPointer[n]->tp);
}

ColumnDescr *getColumnDescr (IRAFPointer n) {
	/* should check whether the type is correct */
	return (dataIRAFPointer[n]->cp);
}

static IRAFPointer findEmptyElement (void) {

	IRAFPointer n=-1;
	int i;
	int status = 0;

	if (nelem < 1) {
	    if ((status = reallocBuffer (INCREMENT)) > 0)
		return n;
	}
	for (i = 0;  i < nelem;  i++) {
	    if (dataIRAFPointer[i] == NULL) {
		n = i;
		break;
	    }
	}
	if (n < 0) {
	    if ((status = reallocBuffer (INCREMENT)) > 0)
		return n;
	    for (i = 0;  i < nelem;  i++) {
		if (dataIRAFPointer[i] == NULL) {
		    n = i;
		    break;
		}
	    }
	}
	return n;
}

static int reallocBuffer (int increment) {
	int new_nelem;
	int i;

	new_nelem = nelem + increment;
	dataIRAFPointer = (TblInfoPtr **)realloc (dataIRAFPointer,
			increment * sizeof(TblInfoPtr));
	if (dataIRAFPointer == NULL)
	    return 1;
	for (i = nelem;  i < new_nelem;  i++)
	    dataIRAFPointer[i] = NULL;
	nelem = new_nelem;
	return 0;
}
