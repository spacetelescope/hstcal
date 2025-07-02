#ifndef HSTCALVERSION_INCL
#define HSTCALVERSION_INCL

#include "config.h"

char * getVersionInfo(char ** buffer);
char * sprintfGitInfo(char ** buffer);
void printGitInfo(void);

#endif
