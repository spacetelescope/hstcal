# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string.h>
# include <ctype.h>

#include "str_util.h"
#include "trlbuf.h"
#include "timestamp.h"


/**
 * This routine prints the current date and time into a string
 * in the format preferred by OPUS, date & time = YYYYDDDHHMMSS.
 *
 * Example:
 * @code{.c}
 * #include "timestamp.h"
 * int main(int argc, char *argv[]) {
 *     InitTrlBuf();
 *     TimeStamp("Generic Conversion Complete", "N2RUSG54B");
 *     CloseTrlBuf(&trlbuf);
 * }
 * @endcode
 *
 * Output:
 * @code
 * 1996212152529-I--------------- Generic Conversion Completed: N2RUSG54B ---------
 * @endcode
 *
 * @param message string to append to date & time info
 * @param rootname root name to include in printed string
*/

void TimeStamp(char *message, const char *rootname) {
    // Separates time string and the message[: rootname]
    const char *message_sep = "-I---------------";

    // Configure timestamp
    // Y - Year (YYYY)
    // j - Day of the year (001-360)
    // H - Hour (00-23)
    // M - Minute (0-59)
    // S - seconds (0-60)
    const char *time_fmt = "%Y%j%H%M%S";
    char *time_str = calloc(SZ_TIMESTRING, sizeof(*time_str));
    if (!time_str) {
        trlerror("(TimeStamp) Out of memory.");
        return;
    }

    // Populate time string
    const time_t tp = time(NULL);
    const struct tm *now = localtime(&tp);
    strftime(time_str, SZ_TIMESTRING, time_fmt, now);

    // Optional root name
    char *uc_rootname = NULL;
    // Length of root name
    size_t uc_rootname_len = 0;

    if (rootname) {
        uc_rootname = strdup(rootname);
        if (uc_rootname == NULL) {
            trlerror("(TimeStamp) Out of memory.");
            return;
        }
        // Convert root name to upper case
        upperCase(uc_rootname);
        uc_rootname_len = strlen(uc_rootname);
    }

    // Populate the output buffer
    const size_t output_max_len = SZ_EIGHTY;
    char *output = calloc(output_max_len + 1, sizeof(*output));
    snprintf(output, output_max_len, "%s%s %s%s%s ",
        time_str,
        message_sep,
        message,
        uc_rootname_len ? ": " : "",
        uc_rootname_len ? uc_rootname: "");

    // Pad remaining bytes in the output buffer with dashes
    repchar_s('-',
        output_max_len - strlen(output),
        &output[strlen(output)],
        output_max_len);

    // Display output
    trlmessage("%s", output);

    // Clean up
    free(uc_rootname);
    free(time_str);
    free(output);
}