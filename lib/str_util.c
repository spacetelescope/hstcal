#include <string.h>
#include <fitsio.h>

#include "str_util.h"

/**
 * Checks whether a character is a delimiter based on
 * an accept and reject string.
 *
 *  @code .c
 *  #include "str_util.h"
 *  int main(int argc, char *argv[]) {
 *      char test[] = "filename.txt";
 *      for (size_t i = 0; i < strlen(test); i++) {
 *          int delim_status = delim_check(filename[i]);
*           if (DELIM_NOT_FOUND) {
*               printf("%c was neither accepted or rejected\n", filename[i]);
*               continue;
*            else if (DELIM_REJECTED) {
*               printf("%c was rejected\n", filename[i]);
*               continue;
*           }
*           // DELIM_FOUND
*           printf("%c was accepted\n", filename[i]);
 *      }
 *  }
 *  @endcode
 *
 * @param ch the character to check
 * @param accept a string of accepted delimiter characters
 * @param reject a string of rejected delimiter characters
 * @return a status code:
 *     - `DELIM_FOUND` if `ch` is in `accept`
 *     - `DELIM_REJECTED` if `ch` is in `reject`
 *     - `DELIM_NOT_FOUND`
 */
int delim_check(const int ch, const char *accept, const char *reject) {
    if (accept && strchr(accept, ch)) {
        return DELIM_FOUND;
    }
    if (reject && strchr(reject, ch)) {
        return DELIM_REJECTED;
    }
    return DELIM_NOT_FOUND;
}

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
    // Not an empty str, use "" or '\0' instead.
    // alphabet := non-empty finite set.
    if (!str || !alphabet || *alphabet=='\0')
        return false;

    // For any set A the empty set is a subset of A
    if (*str=='\0')
        return true;

    for (size_t i = 0; i < strlen(str) ; ++i) {
        if (!strchr(alphabet, str[i]))
            return false;
    }
    return true;
}
