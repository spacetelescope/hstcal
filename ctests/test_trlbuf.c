#include "unittest.h"
#include "trlbuf.h"

// Carve out a chunk of zeroes to compare init/deinit states against
static const struct TrlBuf trlbuf_uninitialized = {0};

// All test functions start off with a fresh GLOBAL trlbuf
static int setup_trl() {
    if (InitTrlBuf()) {
        TEST_MSG(stderr, TEST_TERM_COLOR_BOLD TEST_TERM_COLOR_BRIGHT_RED, "\n", "Unable to initialize trlbuf\n");
        return 1;
    }
    const int init_actually_failed = memcmp(&trlbuf, &trlbuf_uninitialized, sizeof(trlbuf)) == 0;
    if (init_actually_failed) {
        TEST_MSG(stderr, TEST_TERM_COLOR_BOLD TEST_TERM_COLOR_BRIGHT_RED, "\n", "InitTrlBuf returned 0, but trlbuf is not properly initialized\n");
        return 1;
    }

    return 0;
}

// All test functions end with destroying the GLOBAL trlbuf
static int teardown_trl() {
    CloseTrlBuf(&trlbuf);
    const int close_actually_failed = memcmp(&trlbuf, &trlbuf_uninitialized, sizeof(trlbuf)) == 0;
    if (close_actually_failed) {
        TEST_MSG(stderr, TEST_TERM_COLOR_BOLD TEST_TERM_COLOR_BRIGHT_RED, "\n", "CloseTrlBuf called, but trlbuf is not properly deinitialized\n");
        return 1;
    }

    return 0;
}

TEST_UNIT(InitTrlBuf) {
    TEST_ASSERT(trlbuf.init != 0, "trlbuf.init not initialized");
    TEST_ASSERT(trlbuf.fp == NULL, "trlbuf.fp is not NULL");
    TEST_ASSERT(trlbuf.overwrite == 0, "trlbuf.overwrite mode should not be active");
    TEST_ASSERT(trlbuf.quiet == 0, "trlbuf.quiet mode should not be active");
    TEST_ASSERT(strlen(trlbuf.trlfile) == 0, "trlbuf.trlfile should be an empty string");
    TEST_ASSERT(trlbuf.usepref == 1, "trlbuf.usepref should be enabled by default");
    TEST_ASSERT(trlbuf.preface != NULL, "trlbuf.preface should be initialized");
    TEST_ASSERT(strlen(trlbuf.preface) == 0, "trlbuf.preface should be an empty string");
    TEST_ASSERT(trlbuf.buffer != NULL, "trlbuf.buffer should be initialized");
    TEST_ASSERT(strlen(trlbuf.buffer) == 0, "trlbuf.buffer should be an empty string");
    TEST_RETURN
}

TEST_UNIT(InitTrlFile) {
    const char *filename = "output.txt";
    remove(filename);
    const char *test_strings[] = {
        "one",
        "two",
        "three",
    };
    const size_t test_strings_count = sizeof(test_strings)/sizeof(test_strings[0]);

    InitTrlPreface();
    SetTrlPrefaceMode(1);

    // Commit messages to file
    TEST_MARK("Set output file");
    InitTrlFile("", (char *) filename);

    TEST_MARK("Write messages");
    for (size_t i = 0; i < test_strings_count; i++) {
        trlmessage("message %s", test_strings[i]);
    }

    TEST_MARK("check strings written to trailer");
    for (ssize_t i = 0; i < (ssize_t) test_strings_count; i++) {
        char *match = NULL;
        size_t match_lineno = 0;

        const char *pattern = test_strings[i];
        const int have_message = TEST_FILE_CONTAINS(filename, pattern, i, i, &match, &match_lineno) > 0;
        TEST_ASSERT(have_message == 1, "line %zi does not contain expected pattern: '%s' in '%s'", &match_lineno, pattern, match ? match : "NO MATCH");
        free(match);
    }
    remove(filename);

    TEST_RETURN
}

TEST_UNIT(WriteTrlFile) {
    const char *filename = "output.txt";
    TEST_ASSERT(InitTrlFile("", (char *) filename) == 0, "output file creation failed");
    trlmessage("This is a test message");
    TEST_ASSERT(TEST_FILE_CONTAINS(filename, "This is a test message", -1, -1, NULL, NULL) > 0 , "pattern not found in message");
    TEST_ASSERT(trlbuf.fp != NULL, "trlbuf.fp should not be NULL");
    TEST_ASSERT(WriteTrlFile() == 0, "file should close without error");
    TEST_ASSERT(trlbuf.fp == NULL, "trlbuf.fp should be NULL");
    remove(filename);
    TEST_RETURN
}

TEST_UNIT(SetTrlPrefaceMode) {
    TEST_ASSERT(trlbuf.usepref == 1, "preface mode is enabled by default, but isn't");
    SetTrlPrefaceMode(0);
    TEST_ASSERT(trlbuf.usepref == 0, "preface mode is not disabled");
    TEST_RETURN
}

TEST_UNIT(SetTrlOverwriteMode) {
    TEST_ASSERT(trlbuf.overwrite == 0, "overwrite mode should not enabled by default");
    SetTrlOverwriteMode(1);
    TEST_ASSERT(trlbuf.overwrite == 1, "overwrite mode should be enabled");
    TEST_RETURN
}

TEST_UNIT(SetTrlQuietMode) {
    TEST_ASSERT(trlbuf.quiet == 0, "quiet mode should not be enabled");
    SetTrlQuietMode(1);
    TEST_ASSERT(trlbuf.quiet == 1, "quiet mode should enabled");
    TEST_RETURN
}

TEST_UNIT(InitTrlPreface) {
    TEST_ASSERT(trlbuf.preface != NULL, "preface should not be NULL");
    TEST_RETURN
}

TEST_UNIT(ResetTrlPreface) {
    ResetTrlPreface();
    TEST_ASSERT(trlbuf.preface != NULL, "trlbuf.preface should be initialized");
    TEST_ASSERT(strlen(trlbuf.preface) == 0, "trlbuf.preface should be empty");
    TEST_RETURN
}

TEST_UNIT(CloseTrlBuf) {
    CloseTrlBuf(&trlbuf);
    TEST_ASSERT(trlbuf.buffer == NULL, "trlbuf.buffer should be NULL after close");
    TEST_ASSERT(trlbuf.preface == NULL, "trlbuf.preface should be NULL after close");
    TEST_RETURN
}

static const size_t huge_message_size = 512;

TEST_UNIT(trlmessage) {
    const char *filename = "output.txt";
    trlmessage("message");
    trlmessage("message %s", "variadic");
    char *huge_message = malloc(huge_message_size * sizeof(*huge_message));
    if (!huge_message) {
        TEST_THROW_ERROR("unable to allocate memory for huge message");
    }
    memset(huge_message, '?', huge_message_size);
    huge_message[huge_message_size - 1] = '\0';
    trlmessage("%s", huge_message);
    InitTrlFile("", (char *) filename);
    TEST_ASSERT(true, "string too large and truncation not trapped"); // you will never see this if it fails
    free(huge_message);

    remove(filename);
    TEST_RETURN
}

TEST_UNIT(trlwarn) {
    const char *filename = "output.txt";
    trlwarn("message");
    trlwarn("message %s", "variadic");
    char *huge_message = malloc(huge_message_size * sizeof(*huge_message));
    if (!huge_message) {
        TEST_THROW_ERROR("unable to allocate memory for huge message");
    }
    memset(huge_message, '?', huge_message_size);
    huge_message[huge_message_size - 1] = '\0';
    trlwarn("%s", huge_message);
    InitTrlFile("", (char *) filename);
    TEST_ASSERT(true, "string too large and truncation not trapped"); // you will never see this if it fails
    free(huge_message);

    remove(filename);
    TEST_RETURN
}

TEST_UNIT(trlerror) {
    const char *filename = "output.txt";
    trlerror("message");
    trlerror("message %s", "variadic");
    char *huge_message = malloc(huge_message_size * sizeof(*huge_message));
    if (!huge_message) {
        TEST_THROW_ERROR("unable to allocate memory for huge message");
    }
    memset(huge_message, '?', huge_message_size);
    huge_message[huge_message_size - 1] = '\0';
    trlerror("%s", huge_message);
    InitTrlFile("", (char *) filename);
    TEST_ASSERT(true, "string too large and truncation not trapped"); // you will never see this if it fails
    free(huge_message);

    remove(filename);
    TEST_RETURN
}

TEST_UNIT(trlopenerr) {
    // trlopenerr does not check the file, it only prints an error message using trlerror
    const char *filename = "output.txt";
    InitTrlFile("", (char *) filename);
    trlopenerr(filename);

    char *match = NULL;
    size_t match_lineno = 0;
    const char *pattern = "ERROR:    Can't open file output.txt";
    TEST_FILE_CONTAINS(filename, pattern, -1, -1, &match, &match_lineno);
    TEST_ASSERT(match != NULL, "'%s' not found in file '%s'", pattern, filename);
    free(match);
    remove(filename);
    TEST_RETURN
}

TEST_UNIT(trlreaderr) {
    // trlreaderr does not check the file, it only prints an error message using trlerror
    const char *filename = "output.txt";
    InitTrlFile("", (char *) filename);
    trlreaderr(filename);

    char *match = NULL;
    size_t match_lineno = 0;
    const char *pattern = "ERROR:    Can't read file output.txt";
    TEST_FILE_CONTAINS(filename, pattern, -1, -1, &match, &match_lineno);
    TEST_ASSERT(match != NULL, "'%s' not found in file '%s'", pattern, filename);
    free(match);
    remove(filename);
    TEST_RETURN
}

TEST_UNIT(trlkwerr) {
    // trlkwerr does not check the file, it only prints an error message using trlerror
    const char *filename = "output.txt";
    const char *kw = "TEST";
    InitTrlFile("", (char *) filename);
    trlkwerr(kw, filename);

    char *match = NULL;
    size_t match_lineno = 0;
    const char *pattern = "ERROR:    Keyword \"TEST\" not found in output.txt";
    TEST_FILE_CONTAINS(filename, pattern, -1, -1, &match, &match_lineno);
    TEST_ASSERT(match != NULL, "'%s' not found in file '%s'", pattern, filename);
    free(match);
    remove(filename);
    TEST_RETURN
}

TEST_UNIT(trlfilerr) {
    // trlfilerr does not check the file, it only prints an error message using trlerror
    const char *filename = "output.txt";
    InitTrlFile("", (char *) filename);
    trlfilerr(filename);

    char *match = NULL;
    size_t match_lineno = 0;
    const char *pattern = "ERROR:    while trying to read file output.txt";
    TEST_FILE_CONTAINS(filename, pattern, -1, -1, &match, &match_lineno);
    TEST_ASSERT(match != NULL, "'%s' not found in file '%s'", pattern, filename);
    free(match);

    remove(filename);
    TEST_RETURN
}

TEST_UNIT(printfAndFlush) {
    // More of a kernel behavior test than anything. The strings should be written to stdout in order.
    printfAndFlush("flushing first");
    printfAndFlush("flushing second");
    printfAndFlush("flushing third");
    fsync(TEST_LOCAL_REDIRECT.fd_stdout);
    fsync(TEST_LOCAL_REDIRECT.fd_stderr);

    const char *patterns[] = {
        "flushing first\n",
        "flushing second\n",
        "flushing third\n",
    };
    char *match = NULL;
    size_t match_lineno = 0;
    for (size_t i = 0; i < sizeof(patterns)/sizeof(patterns[0]); i++) {
        const char *pattern = patterns[i];
        TEST_FILE_CONTAINS(TEST_LOCAL_REDIRECT.filename, pattern, (ssize_t) i, (ssize_t) i, &match, &match_lineno);
        TEST_ASSERT(match != NULL, "'%s' not found in file '%s'", pattern, TEST_LOCAL_REDIRECT.filename);
        free(match);
    }
    TEST_RETURN
}

TEST_UNIT(trlGitInfo) {
    // NOQA
    trlGitInfo();
    TEST_RETURN
}

TEST_SUITE(__FILE__) {
    TEST_SUITE_SET_FIXTURE_SETUP(setup_trl);
    TEST_SUITE_SET_FIXTURE_TEARDOWN(teardown_trl);

    const unit_test tests[] = {
        TEST_UNIT_REPR(InitTrlBuf),
        TEST_UNIT_REPR(InitTrlFile),
        TEST_UNIT_REPR(WriteTrlFile),
        TEST_UNIT_REPR(SetTrlPrefaceMode),
        TEST_UNIT_REPR(SetTrlOverwriteMode),
        TEST_UNIT_REPR(SetTrlQuietMode),
        TEST_UNIT_REPR(InitTrlPreface),
        TEST_UNIT_REPR(ResetTrlPreface),
        TEST_UNIT_REPR(CloseTrlBuf),
        TEST_UNIT_REPR(trlmessage),
        TEST_UNIT_REPR(trlwarn),
        TEST_UNIT_REPR(trlerror),
        TEST_UNIT_REPR(trlopenerr),
        TEST_UNIT_REPR(trlreaderr),
        TEST_UNIT_REPR(trlkwerr),
        TEST_UNIT_REPR(trlfilerr),
        TEST_UNIT_REPR(printfAndFlush),
        TEST_UNIT_REPR(trlGitInfo),
    };
    TEST_SUITE_RUN(tests);
    TEST_SUITE_RETURN
}
