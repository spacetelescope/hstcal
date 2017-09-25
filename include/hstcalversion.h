#ifndef HSTCALVERSION_INCL
#define HSTCALVERSION_INCL

#ifndef VERSION
#define VERSION "UNKNOWN"
#endif

#ifndef BRANCH
#define BRANCH "UNKNOWN"
#endif

#ifndef COMMIT
#define COMMIT "UNKNOWN"
#endif

char * getVersionInfo(char ** buffer);
char * sprintfGitInfo(char ** buffer);
void printGitInfo(void);

#endif
