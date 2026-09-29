#ifndef HSTCAL_UNITTEST_H
#define HSTCAL_UNITTEST_H

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

/**
 * Valid testfunc return codes
 */
#define TEST_T_PASS 0
#define TEST_T_FAIL 1
#define TEST_T_ERROR 2
#define TEST_T_SKIP 3

#define TEST_STATS_ARRAY_MAX 4

/**
 * Disable line buffering (useful for CI)
 */
#define TEST_DISABLE_BUFFERING() \
    do { \
        setvbuf(stdout, NULL, _IONBF, 0); \
        setvbuf(stderr, NULL, _IONBF, 0); \
    } while (0)

#define TEST_SUITE_BEGIN(NAME) \
    int main(int argc, char *argv[]) { \
        (void) argc; \
        (void) argv; \
        const char *TEST_LOCAL_SUITE_NAME = NAME; \
        int TEST_LOCAL_STATS[TEST_STATS_ARRAY_MAX] = {0}; \
        TEST_MSG(stdout, NULL, TEST_TERM_COLOR_BRIGHT_BLUE " SUITE", "%s%s...%s", TEST_TERM_COLOR_BRIGHT_WHITE, \
            TEST_LOCAL_SUITE_NAME, TEST_TERM_COLOR_RESET);

#define TEST_SUITE_RETURN \
    TEST_STATS_SHOW(); \
    if (TEST_LOCAL_STATS[TEST_T_FAIL] || TEST_LOCAL_STATS[TEST_T_ERROR]) { \
        return 1; \
    } \
    return 0; \
    }

#define TEST_SUITE_RUN(TESTFUNC_ARRAY) \
    TEST_DISABLE_BUFFERING(); \
    const size_t TEST_LOCAL_TESTFUNC_COUNT = sizeof((TESTFUNC_ARRAY)) / sizeof((TESTFUNC_ARRAY)[0]); \
    for (size_t i = 0; i < TEST_LOCAL_TESTFUNC_COUNT; i++) { \
        const testfunc TEST_LOCAL_TESTFUNC = (TESTFUNC_ARRAY)[i]; \
        const int TEST_LOCAL_RESULT = TEST_LOCAL_TESTFUNC(); \
        TEST_STATS_UPDATE(TEST_LOCAL_RESULT); \
    }

/* Generate a function signature suffixed with FN_NAME and configure local variables
 * required by other unit test functions
 *
 * TEST_RETURN must be called to properly close the signature
 *
 * // Basic test definition
 * TEST_BEGIN(mytest) {
 *     // To skip a test that is broken or not implemented
 *     // TEST_FORCE_SKIP
 *
 *     // To fail a test with an error (i.e. can't continue)
 *     // if (condition) {
 *     //     TEST_FORCE_ERROR
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
 * TEST_SUITE_BEGIN() {
 *     testfunc tests[] = {
 *         test_mytest,
 *     };
 *     TEST_SUITE_RUN(tests)
 *     TEST_SUITE_RETURN
 * }
 *
 * @param FN_NAME function name
 */
#define TEST_BEGIN(FN_NAME) \
    static int test_##FN_NAME() { \
        int TEST_LOCAL_ERROR_COUNT = 0; \
        struct TestRedirect TEST_LOCAL_REDIRECT = {0}; \
        printf(TEST_TERM_COLOR_BRIGHT_BLUE "  UNIT" TEST_TERM_COLOR_RESET TEST_TERM_COLOR_BRIGHT_WHITE \
                                           " %s..." TEST_TERM_COLOR_RESET, \
            __func__); \
        TEST_REDIRECT_OUTPUT(&TEST_LOCAL_REDIRECT);

#define TEST_FAILED \
    do { \
        TEST_LOCAL_ERROR_COUNT++; \
    } while (0);

#define TEST_FORCE_ERROR \
    do { \
        TEST_REDIRECT_OUTPUT_RESTORE(&TEST_LOCAL_REDIRECT); \
        printf(TEST_TERM_COLOR_BOLD TEST_TERM_COLOR_BRIGHT_RED " ERROR\n" TEST_TERM_COLOR_RESET); \
        TEST_REDIRECT_OUTPUT_DUMP(&TEST_LOCAL_REDIRECT); \
        return TEST_T_ERROR; \
    } while (0);

#define TEST_FORCE_SKIP \
    do { \
        TEST_REDIRECT_OUTPUT_RESTORE(&TEST_LOCAL_REDIRECT); \
        printf(TEST_TERM_COLOR_YELLOW " SKIPPED\n" TEST_TERM_COLOR_RESET); \
        TEST_REDIRECT_OUTPUT_DUMP(&TEST_LOCAL_REDIRECT); \
        return TEST_T_SKIP; \
    } while (0);

#define TEST_RETURN \
    TEST_REDIRECT_OUTPUT_RESTORE(&TEST_LOCAL_REDIRECT); \
    printf(TEST_TERM_COLOR_RESET " %s\n" TEST_TERM_COLOR_RESET, \
        TEST_LOCAL_ERROR_COUNT ? TEST_TERM_COLOR_RED "FAILED" TEST_TERM_COLOR_RESET \
                               : TEST_TERM_COLOR_GREEN "PASSED" TEST_TERM_COLOR_RESET); \
    TEST_REDIRECT_OUTPUT_DUMP(&TEST_LOCAL_REDIRECT); \
    return (TEST_LOCAL_ERROR_COUNT) ? TEST_T_FAIL : TEST_T_PASS; \
    }

#define TEST_MARK(DESC) \
    do { \
        printf(TEST_TERM_COLOR_BRIGHT_CYAN "MARK" TEST_TERM_COLOR_RESET " %s\n" TEST_TERM_COLOR_RESET, DESC); \
    } while (0)

#define TEST_ASSERT(COND, REASON, ...) \
    do { \
        if (!(COND)) { \
            TEST_MSG(stderr, TEST_TERM_COLOR_RED, TEST_TERM_COLOR_BRIGHT_RED "ASSERTION FAILED", \
                TEST_TERM_COLOR_MAGENTA #COND TEST_TERM_COLOR_RESET TEST_TERM_COLOR_RED \
                " BECAUSE " TEST_TERM_COLOR_RESET REASON, \
                __VA_ARGS__); \
            TEST_FAILED; \
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

/**
 * Returns the length of a character array, ignoring ANSI terminal control codes
 * @param s character array to tally
 * @return total characters
 */
static size_t strlen_sans_ansi_codes(const char *s) {
    const char *p = s;
    size_t count = 0;
    while (*p != '\0') {
        if (p[0] == 0x1b && p[1] == '[') {
            p = &p[2];
            while (*p >= '0' && *p <= '9') {
                p++;
            }
            p++;
            continue;
        }
        count++;
        p++;
    }
    return count;
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
        if (strrchr(line, '\n') == NULL) {
            const size_t logical_len = strlen_sans_ansi_codes(line);
            if (logical_len == 0) {
                // line consists of only control characters and no line feed
                // emit control codes without line information
                printf("%s", line);
                continue;
            }
        }
        printf(TEST_TERM_COLOR_BRIGHT_BLUE "      " TEST_TERM_COLOR_RESET " %s", line);
    }
    if (strlen(r->filename)) {
        remove(r->filename);
    }
    *r->filename = '\0';

    return 0;
}

static inline int TEST_MSG(FILE *stream, const char *color, const char *prefix, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    char *output = NULL;
    int len = 0;
    if ((len = vasprintf(&output, fmt, ap)) < 0) {
        fprintf(stream, "%s format encoding error\n", __func__);
        return -1;
    }
    va_end(ap);

    char *tokens = output;
    char *token = NULL;
    while ((token = strsep(&tokens, "\n")) != NULL) {
        if (strlen(token)) {
            fprintf(stream, TEST_TERM_COLOR_BRIGHT_WHITE "%s" TEST_TERM_COLOR_RESET " %s%s\n" TEST_TERM_COLOR_RESET,
                prefix, color ? color : "", token);
            fprintf(stream, TEST_TERM_COLOR_RESET);
        }
    }
    free(output);
    return len;
}

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

typedef int (*testfunc)(void);
#endif // HSTCAL_UNITTEST_H
