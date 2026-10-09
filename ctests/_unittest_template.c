#include "unittest.h"

static inline int setup() {
    // return 0 on success, or non-zero on error
    return 0;
}

static inline int teardown() {
    // return 0 on success, or non-zero on error
    return 0;
}

TEST_UNIT(its_true) {
    TEST_ASSERT(true == true, "Message when not true");
    TEST_RETURN
}

TEST_UNIT(its_false) {
    TEST_ASSERT(false == false, "Message when not false");
    TEST_RETURN
}

TEST_UNIT(its_always_false) {
    TEST_ASSERT(false == true, "This cannot be");
    TEST_RETURN
}

TEST_UNIT(its_a_failure) {
    if (true) {
        TEST_THROW_FAIL("Something failed and continuing may segfault");
    }
    // Success
    TEST_RETURN
}

TEST_UNIT(its_an_error) {
    if (true) {
        TEST_THROW_ERROR("Something went wrong (allocation error, etc)");
    }
    // No error
    TEST_RETURN
}

TEST_UNIT(its_a_skip) {
    if (true) {
        TEST_THROW_SKIP("Stub");
    }
    // Not skipped
    TEST_RETURN
}

TEST_UNIT(its_checking_on_screen_output) {
    const char *mesg = "hello world\nwhat a fine day it is\ntoday\n";
    printf("%s", mesg);

    // Flush the stream
    fflush(stdout);

    // Check each line, or a range of lines for substrings
    char *match = NULL;
    size_t match_lineno = 0;
    TEST_ASSERT(TEST_FILE_CONTAINS(TEST_LOCAL_REDIRECT.filename, "world", 0, 0, &match, &match_lineno) != 0, "Message not present in output");
    free(match);
    match = NULL;

    TEST_ASSERT(TEST_FILE_CONTAINS(TEST_LOCAL_REDIRECT.filename, "fine", 1, 1, &match, &match_lineno) != 0, "Message not present in output");
    free(match);
    match = NULL;

    TEST_ASSERT(TEST_FILE_CONTAINS(TEST_LOCAL_REDIRECT.filename, "day", 2, 2, &match, &match_lineno) != 0, "Message not present in output");
    free(match);
    match = NULL;
    TEST_RETURN
}

TEST_UNIT(its_checking_on_screen_output_manually) {
    const char *outfile = "capture.txt";
    const char *mesg = "hello world\nwhat a fine day it is\ntoday\n";

    // Write a message to the output log via stdout
    printf("%s", mesg);
    // Save output log
    TEST_REDIRECT_SAVE(&TEST_LOCAL_REDIRECT, outfile);

    // Parse output log
    FILE *fp = fopen(outfile, "rb");
    if (!fp) {
        TEST_THROW_ERROR("Unable to open IO capture log");
    }

    char buf[0x1000] = {0};
    for (size_t i = 0; fgets(buf, sizeof(buf), fp) != NULL; i++) {
        const char *test_cases[] = {
            "world",
            "fine",
            "day",
        };
        strip_ansi_codes(buf);
        if (buf[strlen(buf) - 1] == '\n') {
            buf[strlen(buf) - 1] = '\0';
        }
        TEST_ASSERT(strstr(buf, test_cases[i]) != NULL, "Message '%s' is not present in '%s' on line %zu of %s", test_cases[i], buf, i, outfile);
    }
    fclose(fp);
    remove(outfile);

    TEST_RETURN
}

TEST_SUITE(__FILE__) {
    const unit_test tests[] = {
        TEST_UNIT_REPR(its_true),
        TEST_UNIT_REPR(its_false),
        TEST_UNIT_REPR(its_always_false),
        TEST_UNIT_REPR(its_an_error),
        TEST_UNIT_REPR(its_a_failure),
        TEST_UNIT_REPR(its_a_skip),
        TEST_UNIT_REPR(its_checking_on_screen_output),
        TEST_UNIT_REPR(its_checking_on_screen_output_manually),
    };

    // Configure setup/teardown fixtures to execute before and after each test
    // Not required
    TEST_SUITE_SET_FIXTURE_SETUP(setup);
    TEST_SUITE_SET_FIXTURE_TEARDOWN(teardown);

    TEST_SUITE_RUN(tests)
    TEST_SUITE_RETURN
}