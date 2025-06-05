#include <assert.h>
#include <string.h>
#include <fitsio.h>

#include "str_util.h"

/**
 * Fills a buffer with `ch` up to `max_ch` times.
 * The `dest_size` must be large enough to properly terminate `dest`
 *
 * @param ch character to repeat
 * @param max_ch maximum number of times to repeat `ch`
 * @param dest destination buffer to write characters (must be pre-allocated)
 * @param dest_size size of destination buffer
 * @return  0 on success
 *          -1 if `dest` is `NULL`, or `dest_size` is `0`
 */
int repchar_s(const char ch, size_t max_ch, char * restrict dest, size_t dest_size) {
    if (!dest || !dest_size) {
        return -1;
    }
    const size_t limit = (max_ch < dest_size - 1) ? max_ch : dest_size - 1;
    memset(dest, ch, limit);
    dest[limit] = '\0';
    return 0;
}

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

