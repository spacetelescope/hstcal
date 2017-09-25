#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "hstcalversion.h"

char * sprintfGitInfo(char ** buffer)
{
    if (!buffer)
        return NULL;
    if (*buffer)
        assert(0); // Incorrect usage - NULL ptr must be passed in.

    const char * format = "git tag: %s\ngit branch: %s\nHEAD @: %s";
    size_t length = strlen(format) + strlen(VERSION) + strlen(BRANCH) + strlen(COMMIT); // NOTE: this doesn't require an explicit +1 for '\0' as it is larger than needed due to '%s' & '\n'.
    *buffer = malloc(length*sizeof(char));
    if (!*buffer)
        return NULL;
    sprintf(*buffer, format, VERSION, BRANCH, COMMIT);
    return *buffer;
}

void printGitInfo(void)
{
    char * gitInfo = NULL;
    sprintfGitInfo(&gitInfo);
    printf("%s\n", gitInfo);
    if (gitInfo)
    {
        free(gitInfo);
        gitInfo = NULL;
    }
}
