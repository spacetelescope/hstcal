#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "hstcalversion.h"

char * sprintfGitInfo(char ** buffer)
{
    if (!buffer)
        return NULL;
    if (*buffer)
        return NULL; // Incorrect usage - NULL ptr must be passed in.

    const char * format = "git tag: %s\ngit branch: %s\nHEAD @: %s";
    size_t length = strlen(format)
        + strlen(HSTCAL_VERSION)
        + strlen(HSTCAL_VERSION_BRANCH)
        + strlen(HSTCAL_VERSION_COMMIT); // NOTE: this doesn't require an explicit +1 for '\0' as it is larger than needed due to '%s' & '\n'.
    *buffer = malloc((length + 1) * sizeof(char));
    if (!*buffer)
        return NULL;
    snprintf(*buffer, length, format, HSTCAL_VERSION, HSTCAL_VERSION_BRANCH, HSTCAL_VERSION_COMMIT);
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
