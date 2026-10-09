#include "unittest.h"
#include "str_util.h"

#include <ctype.h>

TEST_UNIT(repchar_s) {
    char test_string[10];
    const char ch = '.';

    // Begin with a raw string and no terminator
    memset(test_string, '?', sizeof(test_string));

    TEST_MARK("fill string with characters");
    const int result = repchar_s(ch, sizeof(test_string), test_string, sizeof(test_string));
    TEST_ASSERT(result == 0, "function returned non-zero, got %d", result);

    TEST_MARK("check all characters are identical");
    int ch_fail = 0;
    for (size_t i = 0; test_string[i] != '\0'; i++) {
        if (test_string[i] != ch) {
            ch_fail = 1;
        }
    }
    TEST_ASSERT(ch_fail == 0, "string contains characters other than '%c', got '%s'", ch, test_string);

    TEST_MARK("check string is nul-terminated");
    TEST_ASSERT(test_string[sizeof(test_string) - 1] == '\0', "did not terminate string");

    TEST_MARK("Fill the string up to the maximum size using loop");
    memset(test_string, '?', sizeof(test_string));
    for (size_t i = 0; i < sizeof(test_string) - 1; i++) {
        const size_t offset = i + 1;
        repchar_s(ch, offset, test_string, sizeof(test_string));
        TEST_ASSERT(test_string[i] == ch, "expected '%c' character, got '%c'", ch, test_string[i]);
        TEST_ASSERT(strlen(test_string) == offset, "unexpected string length, got %zu", strlen(test_string));
        TEST_ASSERT(test_string[offset] == '\0', "function did not terminate string");
        if (TEST_LOCAL_ERROR_COUNT) {
            printf("(max_ch=%zu, strlen=%zu, dest_size=%zu) %s\n", offset, strlen(test_string), sizeof(test_string), test_string);
        }
    }

    TEST_MARK("Attempt to write more bytes than the maximum allowed size");
    memset(test_string, 0, sizeof(test_string));
    repchar_s(ch, sizeof(test_string) + 1024, test_string, sizeof(test_string));
    TEST_ASSERT(strlen(test_string) == sizeof(test_string) - 1, "max_ch safeguard did not prevent buffer overflow");

    TEST_RETURN
}

static size_t count_upper_chars(const char *s) {
    size_t result = 0;
    for (size_t i = 0; s[i] != '\0'; i++) {
        if (isupper((int) s[i])) {
            result++;
        }
    }
    return result;
}

TEST_UNIT(upperCase) {
    char *test_strings[] = {
        ".....a.....",
        " a b c d e ",
        "_\t.\t.\t.\t_",
        "hello world",
        "Hello World",
        "HeLlO wOrLd",
    };
    const size_t test_strings_size = sizeof(test_strings) / sizeof(test_strings[0]);

    TEST_MARK("convert strings to upper-case");
    for (size_t i = 0; i < test_strings_size; i++) {
        char *before = strdup(test_strings[i]);
        if (!before) {
            TEST_THROW_ERROR("strdup failed");
        }
        char *after = strdup(test_strings[i]);
        if (!after) {
            TEST_THROW_ERROR("strdup failed");
        }
        const size_t before_upper_count = count_upper_chars(before);

        // Convert here. upperCase has no return value
        upperCase(after);

        const size_t after_upper_count = count_upper_chars(after);
        TEST_ASSERT(after_upper_count >= before_upper_count, "%zu: string not converted to upper-case\n    before='%s'\n    after ='%s'\n", i, before, after);

        free(before);
        free(after);
    }

    TEST_RETURN
}

TEST_UNIT(isStrInLanguage) {
    char *test_strings[] = {
        "I am HSTCAL",
        "1 will fail",
        "HSTCAL I am",
        "50 w1ll 1",
    };
    const size_t test_strings_size = sizeof(test_strings) / sizeof(test_strings[0]);
    const char *alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ ";

    TEST_MARK("check characters in string exist in alphabet");
    for (size_t i = 0; i < test_strings_size; i++) {
        const int test_strings_expected[] = {
            1,
            0,
            1,
            0,
        };
        const char *test_string = test_strings[i];
        const int expected = test_strings_expected[i];
        TEST_ASSERT(isStrInLanguage(test_string, alphabet) == expected, "'%s' should have returned %d", test_string, expected);
    }

    TEST_MARK("check default behavior");
    TEST_ASSERT(isStrInLanguage("", alphabet) == 1, "should return true on zero-length input string");
    TEST_ASSERT(isStrInLanguage(NULL, alphabet) == 0, "should not return true on NULL input string");
    TEST_ASSERT(isStrInLanguage(NULL, NULL) == 0, "should return false on NULL input string and NULL alphabet");
    TEST_ASSERT(isStrInLanguage("", "") == 0, "should return false on empty input string and alphabet");

    TEST_RETURN
}

TEST_SUITE(__FILE__) {
    const unit_test tests[] = {
        TEST_UNIT_REPR(repchar_s),
        TEST_UNIT_REPR(upperCase),
        TEST_UNIT_REPR(isStrInLanguage),
    };
    TEST_SUITE_RUN(tests);
    TEST_SUITE_RETURN
}