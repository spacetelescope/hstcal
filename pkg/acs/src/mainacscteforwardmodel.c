
#include <stdlib.h>
#include <string.h>
#include <assert.h>

extern int doMainCTE (int argc, char **argv);

int main (int argc, char **argv)
{
    // alloc a new argv + 1 for the above new option
    char ** newArgv = malloc((argc + 1)*sizeof(char*));
    assert(newArgv);

    // Copy the old pointers to the new one
    memcpy(newArgv, argv, argc*sizeof(char*));

    // Add the new opt and inc argc
    newArgv[argc++] = "--forwardModelOnly";

    const int status = doMainCTE(argc, newArgv);
    free(newArgv);
    return status;
}
