#ifndef GETACSKEYS_INCL
#define GETACSKEYS_INCL

# include "hstio.h"
# include "acsinfo.h"

int getACSKeys (ACSInfo *acs, Hdr *phdr);
int getAndCheckACSKeys (ACSInfo *acs, Hdr *phdr);

#endif
