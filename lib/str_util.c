#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <fitsio.h>

#include "str_util.h"

void upperCase(char * str)
{
    if (!str || *str=='\0')
        return;
    fits_uppercase(str);
}

bool isStrInLanguage(const char * str, const char * alphabet)
{
    assert(str); // Not an empty string, use "" or '\0' instead
    assert(alphabet && *alphabet!='\0'); // alphabet := non-empty finite set

    // For any set A the empty set is a subset of A
    if (*str=='\0')
        return true;

    {unsigned i;
    for (i = 0; i < strlen(str) ; ++i)
    {
        if (!strchr(alphabet, str[i]))
            return false;
    }}
    return true;
}

