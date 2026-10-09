#ifndef HSTCAL_UNITTEST_H
#define HSTCAL_UNITTEST_H

#include "str_util.h"


#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>

/**
 * Terminal colors
 */
#define TEST_TERM_COLOR_RESET "\033[0m"
#define TEST_TERM_COLOR_BLACK "\033[30m"
#define TEST_TERM_COLOR_RED "\033[31m"
#define TEST_TERM_COLOR_GREEN "\033[32m"
#define TEST_TERM_COLOR_YELLOW "\033[33m"
#define TEST_TERM_COLOR_BLUE "\033[34m"
#define TEST_TERM_COLOR_MAGENTA "\033[35m"
#define TEST_TERM_COLOR_CYAN "\033[36m"
#define TEST_TERM_COLOR_WHITE "\033[37m"
#define TEST_TERM_COLOR_BOLD "\033[1m"
#define TEST_TERM_COLOR_BRIGHT_BLACK "\033[90m"
#define TEST_TERM_COLOR_BRIGHT_RED "\033[91m"
#define TEST_TERM_COLOR_BRIGHT_GREEN "\033[92m"
#define TEST_TERM_COLOR_BRIGHT_YELLOW "\033[93m"
#define TEST_TERM_COLOR_BRIGHT_BLUE "\033[94m"
#define TEST_TERM_COLOR_BRIGHT_MAGENTA "\033[95m"
#define TEST_TERM_COLOR_BRIGHT_CYAN "\033[96m"
#define TEST_TERM_COLOR_BRIGHT_WHITE "\033[97m"
#define TEST_TERM_COLOR_BG_BLACK "\033[40m"
#define TEST_TERM_COLOR_BG_RED "\033[41m"
#define TEST_TERM_COLOR_BG_GREEN "\033[42m"
#define TEST_TERM_COLOR_BG_YELLOW "\033[43m"
#define TEST_TERM_COLOR_BG_BLUE "\033[44m"
#define TEST_TERM_COLOR_BG_MAGENTA "\033[45m"
#define TEST_TERM_COLOR_BG_CYAN "\033[46m"
#define TEST_TERM_COLOR_BG_WHITE "\033[47m"

// Valid unit_test_fixture array indicies
enum {
    TEST_FUNC_SETUP=0,
    TEST_FUNC_TEARDOWN,
    TEST_FUNC_ARRAY_MAX,
};

/**
 * Valid unit_test return codes
 */
enum {
    TEST_T_PASS=0,
    TEST_T_FAIL,
    TEST_T_ERROR,
    TEST_T_SKIP,
    TEST_STATS_ARRAY_MAX,
};

/**
 * Disable line buffering (useful for CI)
 */
#define TEST_DISABLE_BUFFERING() \
    do { \
        setvbuf(stdout, NULL, _IONBF, 0); \
        setvbuf(stderr, NULL, _IONBF, 0); \
    } while (0)

#define TEST_SUITE(NAME) \
    int main(int argc, char *argv[]) { \
        (void) argc; \
        (void) argv; \
        const char *TEST_LOCAL_SUITE_NAME = (NAME); \
        int TEST_LOCAL_STATS[TEST_STATS_ARRAY_MAX] = {0}; \
        TEST_MSG(stdout, NULL, TEST_TERM_COLOR_BRIGHT_BLUE "SUITE", " %s%s...%s", TEST_TERM_COLOR_BRIGHT_WHITE, \
            TEST_LOCAL_SUITE_NAME, TEST_TERM_COLOR_RESET);

#define TEST_SUITE_RETURN \
    TEST_STATS_SHOW(); \
    if (TEST_LOCAL_STATS[TEST_T_FAIL] || TEST_LOCAL_STATS[TEST_T_ERROR]) { \
        return 1; \
    } \
    return 0; \
    } // DO NOT REMOVE THIS BRACE

#define TEST_SUITE_RUN(TESTFUNC_ARRAY) \
    TEST_DISABLE_BUFFERING();          \
    const size_t TEST_LOCAL_TESTFUNC_COUNT = sizeof((TESTFUNC_ARRAY)) / sizeof((TESTFUNC_ARRAY)[0]); \
    for (size_t i = 0; i < TEST_LOCAL_TESTFUNC_COUNT; i++) { \
        const unit_test TEST_LOCAL_TESTFUNC = (TESTFUNC_ARRAY)[i]; \
        const int TEST_LOCAL_RESULT = TEST_LOCAL_TESTFUNC(); \
        TEST_STATS_UPDATE(TEST_LOCAL_RESULT); \
    }

#define TEST_SUITE_SET_FIXTURE_SETUP(func) \
    do { \
        TEST_SUITE_FIXTURES[TEST_FUNC_SETUP] = (func); \
    } while (0)
#define TEST_SUITE_SET_FIXTURE_TEARDOWN(func) \
    do { \
        TEST_SUITE_FIXTURES[TEST_FUNC_TEARDOWN] = (func); \
    } while (0)

#define TEST_SUITE_FIXTURE_RUN(FIXTURE) \
    do { \
        if (TEST_SUITE_FIXTURES[FIXTURE]) { \
            TEST_LOCAL_FIXTURE_RESULT = TEST_SUITE_FIXTURES[FIXTURE](); \
            if (TEST_LOCAL_FIXTURE_RESULT) { \
                TEST_THROW_ERROR("An error occurred in fixture %d", FIXTURE); \
            } \
        } \
    } while (0)

/* Generate a function signature suffixed with FUNC_NAME and configure local variables
 * required by other unit test functions
 *
 * TEST_RETURN must be called to properly close the signature
 *
 * // Basic test definition
 * TEST_UNIT(mytest) {
 *     // To skip a test that is broken or not implemented
 *     // TEST_THROW_SKIP("Reason here")
 *
 *     // To fail a test with an error (i.e. can't continue)
 *     // if (condition) {
 *     //     TEST_THROW_ERROR("Because condition");
 *     // }
 *
 *     // Asserts
 *     TEST_ASSERT(value == 1234, "expected %d, got %d", value);
 *     TEST_ASSERT(stuff == 4321, "%s", "nothing to format here");
 *
 *     // Clean up
 *     TEST_RETURN
 * }
 *
 *
 * TEST_SUITE() {
 *     unit_test tests[] = {
 *         TEST_UNIT_REPR(mytest),
 *     };
 *     TEST_SUITE_RUN(tests)
 *     TEST_SUITE_RETURN
 * }
 *
 * @param FUNC_NAME function name
 */
#define TEST_UNIT(FUNC_NAME) \
    static int test__##FUNC_NAME() { \
        int TEST_LOCAL_FIXTURE_RESULT = 0; \
        int TEST_LOCAL_ERROR_COUNT = 0; \
        struct TestRedirect TEST_LOCAL_REDIRECT = {0}; \
        printf(TEST_TERM_COLOR_BRIGHT_BLUE "  UNIT" TEST_TERM_COLOR_RESET TEST_TERM_COLOR_BRIGHT_WHITE \
            " %s... " TEST_TERM_COLOR_RESET, \
            __func__); \
        TEST_REDIRECT_OUTPUT(&TEST_LOCAL_REDIRECT); \
        TEST_SUITE_FIXTURE_RUN(TEST_FUNC_SETUP);

#define TEST_UNIT_REPR(FUNC_NAME) test__##FUNC_NAME
#define REPR(FUNC_NAME) test__##FUNC_NAME

#define TEST_FAILED \
    do { \
        TEST_LOCAL_ERROR_COUNT++; \
    } while (0);

#define TEST_RETURN \
        do { \
            if (TEST_SUITE_FIXTURES[TEST_FUNC_TEARDOWN]) { \
                TEST_SUITE_FIXTURE_RUN(TEST_FUNC_TEARDOWN); \
            } \
            if (TEST_LOCAL_ERROR_COUNT != 0) { \
                return TEST_THROW(TEST_T_FAIL); \
            } \
            return TEST_THROW(TEST_T_PASS); \
        } while (0); \
    }

#define TEST_MARK(DESC) \
    do { \
        printf(TEST_TERM_COLOR_BRIGHT_CYAN "MARK" TEST_TERM_COLOR_RESET " %s\n" TEST_TERM_COLOR_RESET, DESC); \
    } while (0)

#define TEST_ASSERT(COND, REASON, ...) \
    do { \
        if (!(COND)) { \
            TEST_MSG(stderr, TEST_TERM_COLOR_BRIGHT_RED "ASSERTION FAILED", "", \
                TEST_TERM_COLOR_MAGENTA " " #COND TEST_TERM_COLOR_RESET TEST_TERM_COLOR_RED \
                " BECAUSE " TEST_TERM_COLOR_RESET REASON, \
                ##__VA_ARGS__); \
            TEST_FAILED \
        } \
    } while (0)

struct TestRedirect {
    int fd_stdout;
    int fd_stderr;
    char filename[4096];
};

/**
 * Redirect stdout/stderr to a temporary file. Save file descriptors and file
 * name in TestRedirect structure.
 * @param r TestRedirect structure
 * @return 0 on success
 */
static inline int TEST_REDIRECT_OUTPUT(struct TestRedirect *r) {
    r->fd_stdout = dup(STDOUT_FILENO);
    if (r->fd_stdout < 0) {
        fprintf(stderr, "Unable to dup() stdout\n");
        return -1;
    }

    r->fd_stderr = dup(STDERR_FILENO);
    if (r->fd_stderr < 0) {
        fprintf(stderr, "Unable to dup() stderr\n");
        return -1;
    }

    char template[] = "test_output_XXXXXX";
    const int fd = mkstemp(template);
    if (fd < 0) {
        fprintf(stderr, "Unable to create temporary file\n");
        return -1;
    }
    strncpy(r->filename, template, sizeof(r->filename) - 1);
    r->filename[sizeof(r->filename) - 1] = '\0';

    fflush(stdout);
    fflush(stderr);

    if (dup2(fd, STDOUT_FILENO) < 0) {
        fprintf(stderr, "Unable to dup2() stdout to %s\n", template);
        return -1;
    }
    if (dup2(fd, STDERR_FILENO) < 0) {
        fprintf(stderr, "Unable to dup2() stderr to %s\n", template);
        return -1;
    }
    return close(fd);
}

/**
 * Restores stdout/stderr file descriptors in TestRedirect structure
 * @param r pointer to TestRedirect
 * @return 0 on success
 */
static inline int TEST_REDIRECT_OUTPUT_RESTORE(struct TestRedirect *r) {
    fflush(stdout);
    fflush(stderr);

    if (dup2(r->fd_stdout, STDOUT_FILENO) < 0) {
        fprintf(stderr, "Unable to restore stdout\n");
        return -1;
    }
    close(r->fd_stdout);
    r->fd_stdout = -1;

    if (dup2(r->fd_stderr, STDERR_FILENO) < 0) {
        fprintf(stderr, "Unable to restore stderr\n");
        return -1;
    }
    close(r->fd_stderr);
    r->fd_stderr = -1;

    return 0;
}

// ANSI Control Sequence Introducer (CSI)
#define TEST_ANSI_CSI_ESC_1 0x1b        // ESC
#define TEST_ANSI_CSI_ESC_2 0x5b        // '['

#define TEST_ANSI_CSI_PARAM_LOW 0x30    // '0'
#define TEST_ANSI_CSI_PARAM_HIGH 0x3f   // '?'

#define TEST_ANSI_CSI_INTR_LOW 0x20     // ' '
#define TEST_ANSI_CSI_INTR_HIGH 0x2F    // '/'

#define TEST_ANSI_CSI_FINAL_LOW 0x40      // '@'
#define TEST_ANSI_CSI_FINAL_HIGH 0x7e     // '~'


static inline void strip_ansi_codes(char *s) {
    char *src = s;
    char *dst = s;

    while (*src != '\0') {
        if (src[0] == TEST_ANSI_CSI_ESC_1 && src[1] == TEST_ANSI_CSI_ESC_2) {
            // skip preamble
            char *p = src + 2;

            // skip parameter byte range
            while (*p >= TEST_ANSI_CSI_PARAM_LOW && *p <= TEST_ANSI_CSI_PARAM_HIGH) {
                p++;
            }

            // skip intermediate byte range
            while (*p >= TEST_ANSI_CSI_INTR_LOW && *p <= TEST_ANSI_CSI_INTR_HIGH) {
                p++;
            }

            // check final command byte
            if (*p >= TEST_ANSI_CSI_FINAL_LOW && *p <= TEST_ANSI_CSI_FINAL_HIGH) {
                // jump src to after the final byte
                src = p + 1;
                continue;
            }
        }
        // copy in place
        *dst++ = *src++;
    }
    *dst = '\0';
}

/**
 * Returns the length of a character array, ignoring ANSI terminal control codes
 * @param s character array to tally
 * @return total characters
 */
static size_t strlen_sans_ansi_codes(const char *s) {
    size_t count = 0;

    while (*s != '\0') {
        if (s[0] == TEST_ANSI_CSI_ESC_1 && s[1] == TEST_ANSI_CSI_ESC_2) {
            const char *p = s;
            s += 2;

            while (*s >= TEST_ANSI_CSI_PARAM_LOW && *s <= TEST_ANSI_CSI_PARAM_HIGH) {
                s++;
            }
            while (*s >= TEST_ANSI_CSI_INTR_LOW && *s <= TEST_ANSI_CSI_INTR_HIGH) {
                s++;
            }

            if (*s >= TEST_ANSI_CSI_FINAL_LOW && *s <= TEST_ANSI_CSI_FINAL_HIGH) {
                s++;
                continue;
            }

            s = p;
        }

        count++;
        s++;
    }

    return count;
}

static inline int is_ansi_and_empty(const char *s) {
    if (strrchr(s, '\n') == NULL) {
        const size_t logical_len = strlen_sans_ansi_codes(s);
        if (logical_len == 0) {
            return 1;
        }
    }
    return 0;
}

static inline void hexdump(const void *data, const size_t maxlen) {
    const unsigned char *bytes = data;
    const size_t bytes_per_line = 16;

    for (size_t i = 0; i < maxlen; i += bytes_per_line) {
        size_t bytes_count;
        if (maxlen - i < bytes_per_line) {
            bytes_count = maxlen - i;
        } else {
            bytes_count = bytes_per_line;
        }

        // print offset
        printf("%08zx  ", i);
        for (size_t j = 0; j < bytes_per_line; ++j) {
            if (j < bytes_count) {
                // print byte
                printf("%02x ", bytes[i + j]);
            } else {
                // print padding
                printf("   ");
            }

            if (j == (bytes_per_line - 1) / 2) {
                // pad center of output
                printf(" ");
            }
        }

        printf(" ");
        for (size_t j = 0; j < bytes_count; ++j) {
            const unsigned char c = bytes[i + j];
            putchar(isprint(c) ? c : '.');
        }
        printf("\n");
    }
}

static inline int TEST_REDIRECT_OUTPUT_DUMP(struct TestRedirect *r) {
    fflush(stdout);
    fflush(stderr);
    FILE *fp = fopen(r->filename, "rb");
    if (!fp) {
        fprintf(stderr, "Unable to open log file: %s\n", r->filename);
        return -1;
    }

    char line[0x1000] = {0};
    while (fgets(line, sizeof(line) - 1, fp) != NULL) {
        if (is_ansi_and_empty(line)) {
            fprintf(stdout, "%s", line);
            continue;
        }
        fprintf(stdout, TEST_TERM_COLOR_BRIGHT_BLUE "      " TEST_TERM_COLOR_RESET " %s", line);
    }
    fclose(fp);

    if (strlen(r->filename)) {
        remove(r->filename);
    }
    *r->filename = '\0';

    return 0;
}

static inline int TEST_REDIRECT_SAVE(struct TestRedirect *r, const char *filename) {
    TEST_REDIRECT_OUTPUT_RESTORE(r);
    FILE *outstream = fopen(filename, "w+");
    if (!outstream) {
        fprintf(stderr, "Unable to open output file: %s\n", filename);
        return -1;
    }

    FILE *fp = fopen(r->filename, "rb");
    if (!fp) {
        fprintf(stderr, "Unable to open log file: %s\n", r->filename);
        fclose(outstream);
        return -1;
    }

    char line[0x1000] = {0};
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (is_ansi_and_empty(line)) {
            fprintf(outstream, "%s", line);
            continue;
        }
        fprintf(outstream, "%s", line);
    }

    fclose(fp);
    fclose(outstream);
    TEST_REDIRECT_OUTPUT(r);
    return 0;
}

#define TEST_THROW(ACTION) TEST_THROW_(ACTION, &TEST_LOCAL_REDIRECT)

static inline int TEST_THROW_(const int action, struct TestRedirect *r) {
    const char *color;
    const char *action_msg;
    switch (action) {
        case TEST_T_PASS:
            action_msg = "PASSED";
            color = TEST_TERM_COLOR_GREEN;
            break;
        case TEST_T_FAIL:
            action_msg = "FAILED";
            color = TEST_TERM_COLOR_RED;
            break;
        case TEST_T_SKIP:
            action_msg = "SKIPPED";
            color = TEST_TERM_COLOR_YELLOW;
            break;
        case TEST_T_ERROR:
            action_msg = "ERROR";
            color = TEST_TERM_COLOR_BOLD TEST_TERM_COLOR_RED;
            break;
        default:
            action_msg = "UNHANDLED ACTION";
            color = TEST_TERM_COLOR_RED;
            break;
    }
    TEST_REDIRECT_OUTPUT_RESTORE(r);
    printf("%s%s\n" TEST_TERM_COLOR_RESET, color, action_msg);
    TEST_REDIRECT_OUTPUT_DUMP(r);
    return action;
}

static inline int TEST_VMSG(FILE *stream, const char *color, const char *prefix, const char *fmt, va_list ap) {
    va_list ap_copy;
    va_copy(ap_copy, ap);
    char *output = NULL;
    int len = 0;
    if ((len = vasprintf(&output, fmt, ap_copy)) < 0) {
        fprintf(stream, "%s format encoding error\n", __func__);
        return -1;
    }
    va_end(ap_copy);

    char *tokens = output;
    char *token = NULL;
    while ((token = strsep(&tokens, "\n")) != NULL) {
        if (strlen(token)) {
            fprintf(stream, TEST_TERM_COLOR_BRIGHT_WHITE "%s" TEST_TERM_COLOR_RESET "%s%s%s\n" TEST_TERM_COLOR_RESET,
                prefix && strlen(prefix) ? " " : "",
                prefix && strlen(prefix) ? prefix : "",
                color ? color : "", token);
            fprintf(stream, TEST_TERM_COLOR_RESET);
        }
    }
    free(output);
    return len;
}

static inline int TEST_MSG(FILE *stream, const char *color, const char *prefix, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    const int len = TEST_VMSG(stream, color, prefix, fmt, ap);
    va_end(ap);
    return len;
}

#define TEST_THROW_ERROR(MESSAGE, ...) \
    do { \
        if (MESSAGE && strlen(MESSAGE)) { \
            TEST_MSG(stderr, TEST_TERM_COLOR_BRIGHT_CYAN "EXCEPTION " TEST_TERM_COLOR_BOLD TEST_TERM_COLOR_BRIGHT_RED, "", MESSAGE, ##__VA_ARGS__); \
        } \
        return TEST_THROW(TEST_T_ERROR); \
    } while (0)

#define TEST_THROW_FAIL(MESSAGE, ...) \
    do { \
        if (MESSAGE && strlen(MESSAGE)) { \
            TEST_MSG(stderr, TEST_TERM_COLOR_BRIGHT_CYAN "EXCEPTION " TEST_TERM_COLOR_RED, "", MESSAGE, ##__VA_ARGS__); \
        } \
        return TEST_THROW(TEST_T_FAIL); \
    } while (0)

#define TEST_THROW_SKIP(MESSAGE, ...) \
    do { \
        if (MESSAGE && strlen(MESSAGE)) { \
            TEST_MSG(stderr, TEST_TERM_COLOR_BRIGHT_CYAN "REASON " TEST_TERM_COLOR_YELLOW, "", MESSAGE, ##__VA_ARGS__); \
        } \
        TEST_SUITE_FIXTURE_RUN(TEST_FUNC_TEARDOWN); \
        return TEST_THROW(TEST_T_SKIP); \
    } while (0)

#define TEST_STATS_UPDATE(RES) \
    do { \
        TEST_STATS_UPDATE_(TEST_LOCAL_STATS, (RES)); \
    } while (0);

static inline void TEST_STATS_UPDATE_(int *stats, const int result) {
    switch (result) {
        case TEST_T_FAIL:
            stats[TEST_T_FAIL]++;
            break;
        case TEST_T_ERROR:
            stats[TEST_T_ERROR]++;
            break;
        case TEST_T_PASS:
            stats[TEST_T_PASS]++;
            break;
        case TEST_T_SKIP:
            stats[TEST_T_SKIP]++;
            break;
        default:
            fprintf(stderr, "invalid result\n");
            break;
    }
}

#define TEST_STATS_SHOW() \
    do { \
        TEST_STATS_SHOW_(TEST_LOCAL_STATS, TEST_LOCAL_SUITE_NAME); \
    } while (0);

static inline void TEST_STATS_SHOW_(const int *stats, const char *name) {
    size_t total_tests = 0;
    for (size_t i = 0; i < TEST_STATS_ARRAY_MAX; i++) {
        total_tests += stats[i];
    }
    printf("\n");
    printf("%sREPORT%s %s%s%s\n\n", TEST_TERM_COLOR_BRIGHT_BLUE, TEST_TERM_COLOR_RESET, TEST_TERM_COLOR_BRIGHT_WHITE,
        name, TEST_TERM_COLOR_RESET);
    printf("%s%-6s%s... %-8zu\n", TEST_TERM_COLOR_BRIGHT_WHITE, "Tests", TEST_TERM_COLOR_RESET, total_tests);
    printf("%s%-6s%s... %-8d\n", TEST_TERM_COLOR_GREEN, "Pass", TEST_TERM_COLOR_RESET, stats[TEST_T_PASS]);
    printf("%s%-6s%s... %-8d\n", TEST_TERM_COLOR_RED, "Fail", TEST_TERM_COLOR_RESET, stats[TEST_T_FAIL]);
    printf("%s%-6s%s... %-8d\n", TEST_TERM_COLOR_BOLD TEST_TERM_COLOR_BRIGHT_RED, "Error", TEST_TERM_COLOR_RESET,
        stats[TEST_T_ERROR]);
    printf("%s%-6s%s... %-8d\n", TEST_TERM_COLOR_YELLOW, "Skip", TEST_TERM_COLOR_RESET, stats[TEST_T_SKIP]);
}

static inline char **TEST_FILE_AS_ARRAY(const char *filename, size_t *lines_count) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        return NULL;
    }

    size_t used = 0;
    size_t alloc = 1024;
    char **arr = calloc(alloc + 1,sizeof(*arr));

    size_t line_len = 0;
    char *line = NULL;
    while (getline(&line, &line_len, fp) != -1) {
        arr[used] = line;
        used++;
        line = NULL;

        if (used >= alloc) {
            alloc *= 2;
            char **tmp = realloc(arr, alloc + 1 * sizeof(*arr));
            if (!tmp) {
                return NULL;
            }
            arr = tmp;
        }
        arr[used] = NULL;
    }
    free(line); // final pointer from getline needs to be freed
    line = NULL;

    *lines_count = used;
    fclose(fp);

    return arr;
}

#define TEST_ARRAY_FREE(PTR, COUNT) \
    do { \
        for (size_t i = 0; i < (COUNT); i++) { \
            free((PTR[i])); \
            (PTR[i]) = NULL; \
        } \
        free(PTR); \
        (PTR) = NULL; \
    } while (0)

/**
 * Find a substring in a file
 *
 * @param filename path to file
 * @param pattern a substring to match
 * @param range_start ignore lines before this line (zero-index, -1 for all)
 * @param range_end ignore lines after this line (zero-index, -1 for all)
 * @param result return the string containing the pattern (must be freed by caller)
 * @param result_lineno the line number where the pattern matched (zero-index)
 * @return 0 = not found, 1 = found, -1 = error
 */
static inline int TEST_FILE_CONTAINS(const char *filename, const char *pattern, const ssize_t range_start, const ssize_t range_end, char **result, size_t *result_lineno) {
    size_t line_count = 0;
    char **lines = TEST_FILE_AS_ARRAY(filename, &line_count);
    for (ssize_t i = 0; i < (ssize_t) line_count; i++) {
        if (range_end >= 0 && i > range_end) {
            break;
        }
        if (range_start >= 0 && i < range_start) {
            continue;
        }
        const char *match = strstr(lines[i], pattern);
        if (match) {
            if (result) {
                *result = strdup(lines[i]);
                if (!*result) {
                    return -1;
                }
                if (result_lineno) {
                    *result_lineno = i;
                }
            }
            goto found;
        }
    }
    goto not_found;

    found:
    TEST_ARRAY_FREE(lines, line_count);
    return 1;

    not_found:
    if (result) {
        *result = NULL;
    }
    TEST_ARRAY_FREE(lines, line_count);
    return 0;
}

typedef int (*unit_test)(void);
typedef int (*unit_test_fixture)(void);

unit_test_fixture TEST_SUITE_FIXTURES[TEST_FUNC_ARRAY_MAX] = {NULL, NULL};

#endif // HSTCAL_UNITTEST_H
