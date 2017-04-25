/*
 MLS 2015 This code  pulls keyword switches for the CTE routines

*/

# include <stdio.h>
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

int GetSw (Hdr *, char *, int *);


/* Read the value fo PCTECORR. OMIT and PERFORM are the only possible values
 returned by this function. 
*/

int GetCTESwitch (WF3Info *wf3, Hdr *phdr){
    
    extern int status;
    FitsKw key;       
     
    key = findKw (phdr, "PCTECORR");
    if (key == NotFound) {
        sprintf(MsgText, "PCTECORR keyword not found...");
        trlwarn(MsgText);
    }
    if (GetSw(phdr, "PCTECORR", &wf3->pctecorr))
        return (status);
        
    
    key = findKw (phdr, "BIASCORR");
    if (key == NotFound) {
        sprintf(MsgText, "BIASCORR keyword not found...");
        trlwarn(MsgText);
    }
    if (GetSw (phdr, "BIASCORR", &wf3->biascorr))
        return (status);

    return (status);
}


/* This routine just calls GetSwitch to get the value of the switch,
   but then the value returned is limited to OMIT (0) and PERFORM (1).
*/
int GetSw (Hdr *phdr, char *calswitch, int *flag) {

    /* arguments:
       Hdr *phdr        i: primary header
       char *calswitch  i: name of keyword (e.g. FLATCORR)
       int *flag        o: value (0 or 1) of calibration switch
    */

    extern int status;

    int GetSwitch (Hdr *, char *, int *);

    if (GetSwitch (phdr, calswitch, flag))
        return (status);

    if (*flag != PERFORM)
        *flag = OMIT;

    return (status);
}

