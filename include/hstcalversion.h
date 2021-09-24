#ifndef HSTCALVERSION_INCL
#define HSTCALVERSION_INCL

#include "version.h"

char * getVersionInfo(char ** buffer);
char * sprintfGitInfo(char ** buffer);
void printGitInfo(void);

#endif
