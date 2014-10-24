# include <stdio.h>
# include <string.h>

# include "xtables.h"

# include "stis.h"
# include "calstis6.h"

/* 
   Print message at the beginning of processing with names of reference 
   files used, and their pedigree and descrip keywords.



   Revision history:
   ----------------
   20 Feb 97  -  Implemented, with some code adapted from similar routines
                 in calstis7 (I.Busko)
   25 Apr 97  -  Removed ESPTAB (IB)
   09 May 97  -  New _trl standard (IB)
   15 Dec 97  -  Removed sdctab reference (IB) 
   26 Jan 98  -  PCTAB support (IB)
   08 Apr 05  -  Print GACTAB info (PEH)
*/


void Message6 (StisInfo6 *sts, int type) {

	switch (type) {

	case XTRAC_INFO:
	    if (sts->x1d == PERFORM) {
	        PrRefInfo ("xtractab", sts->xtrctab.name,
                                       sts->xtrctab.pedigree,
                                       sts->xtrctab.descrip,
                                       sts->xtrctab.descrip2);
	        PrRefInfo ("sptrctab", sts->sptrctab.name,
                                       sts->sptrctab.pedigree,
                                       sts->sptrctab.descrip,
                                       sts->sptrctab.descrip2);
	    }
	    break;

	case DISP_INFO:
	    if (sts->dispcorr == PERFORM && !(sts->dispinfo)) {
	        sts->dispinfo = 1;
	        PrRefInfo ("disptab",  sts->disptab.name,
                                       sts->disptab.pedigree,
                                       sts->disptab.descrip,
                                       sts->disptab.descrip2);
	        PrRefInfo ("apdestab", sts->apdestab.name,
                                       sts->apdestab.pedigree,
                                       sts->apdestab.descrip,
                                       sts->apdestab.descrip2);
	        PrRefInfo ("inangtab", sts->inangtab.name,
                                       sts->inangtab.pedigree,
                                       sts->inangtab.descrip,
                                       sts->inangtab.descrip2);
	    }
	    break;

	case FLUX_INFO:
	    if (sts->fluxcorr == PERFORM && !(sts->fluxinfo)) {
	        sts->fluxinfo = 1;
	        PrRefInfo ("phottab", sts->phottab.name,
                                      sts->phottab.pedigree,
                                      sts->phottab.descrip,
                                      sts->phottab.descrip2);
	        PrRefInfo ("apertab", sts->apertab.name,
                                      sts->apertab.pedigree,
                                      sts->apertab.descrip,
                                      sts->apertab.descrip2);
	        if (sts->gaccorr == PERFORM) {
	            PrRefInfo ("gactab", sts->gactab.name,
                                         sts->gactab.pedigree,
                                         sts->gactab.descrip,
                                         sts->gactab.descrip2);
	        }
	        /* Only print if the table exists ! */
	        if (strlen (sts->pctab.name) > 0) {
	            PrRefInfo ("pctab",   sts->pctab.name,
                                          sts->pctab.pedigree,
                                          sts->pctab.descrip,
                                          sts->pctab.descrip2);
	        }
	        if (sts->tdscorr == PERFORM) {
	            PrRefInfo ("tdstab", sts->tdstab.name,
                                         sts->tdstab.pedigree,
                                         sts->tdstab.descrip,
                                         sts->tdstab.descrip2);
	        } else
	            printf ("Warning  TDSTAB correction not performed.\n");

	    }
	    break;

	case SGEO_INFO:
	    if (sts->sgeocorr == PERFORM && !(sts->sgeoinfo)) {
	        sts->sgeoinfo = 1;
	        printf ("SDSTFILE %s\n", sts->sdstfile.name);
	        printf ("SDSTFILE PEDIGREE=%s\n", sts->sdstfile.pedigree);
	        printf ("SDSTFILE DESCRIP =%s\n", sts->sdstfile.descrip);
	    }
	    break;

        case CCD_INFO:
            if (sts->ctecorr == PERFORM && !(sts->ccdinfo)) {
                sts->ccdinfo = 1;
	        PrRefInfo ("ccdtab", sts->ccdtab.name,
                           sts->ccdtab.pedigree,
                           sts->ccdtab.descrip,
                           sts->ccdtab.descrip2);
            }
           
	}
}
