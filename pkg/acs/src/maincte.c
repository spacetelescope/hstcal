#include <stdlib.h>
#include <string.h>

extern int doMainCTE (int argc, char **argv);

int main (int argc, char **argv)
{
    const int status = doMainCTE(argc, argv);
    return status;
}
