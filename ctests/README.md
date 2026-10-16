# UNITTEST

A lightweight single-header unit testing framework for C.

## Features

- Single-header library
- Test suites and test cases
- Optional setup and teardown fixtures
- Assertions
- Skip, fail, and error handling
- Capture and examine test output

## Quickstart

```c
#include "unittest.h"

TEST_UNIT(addition) {
    TEST_ASSERT(1 + 1 == 2, "expected 2");
    TEST_RETURN
}

TEST_UNIT(subtraction) {
    TEST_ASSERT(5 - 2 == 3, "expected 3");
    TEST_RETURN
}

TEST_SUITE("Math Operations") {
    unit_test tests[] = {
        TEST_UNIT_REPR(addition),
        TEST_UNIT_REPR(subtraction),
    };
    TEST_SUITE_RUN(tests)
    TEST_SUITE_RETURN
}
```

### Build

```c
gcc -Wall -Wextra -o test_math_operations tests_math_operations.c
```

### Output

```c
./test_math_operations
UNITTEST v1.0.0
Copyright (C) 2026 Association of Universities for Research in Astronomy (AURA)

 SUITE Math Operations...
  UNIT test__addition... PASSED
  UNIT test__subtraction... PASSED

REPORT Math Operations

Tests ... 2
Pass  ... 2
Fail  ... 0
Error ... 0
Skip  ... 0
```

# Assertions

## Basic Assertions

```c
TEST_ASSERT(x != y, "x should never be equal to y")
```

## Formatted Assertions

```c
int expected = 5;
int result = 6;
TEST_ASSERT(result == expected, "expected %d, got %d", expected, result);
```

### Output

```c
ASSERTION FAILED count == expected BECAUSE expected 5, got 6
```

# Test Status Handling

`TEST_ASSERT` increments an internal error counter whenever the tested condition returns `false`.
`TEST_RETURN` returns `TEST_T_PASS` when the internal counter is zero, or `TEST_T_FAIL` otherwise.

## Pass

```c
TEST_UNIT(zero) {
    TEST_ASSERT(0 == 0, "expected zero!");
    TEST_RETURN  // pass
}
```

## Fail

```c
TEST_UNIT(not_zero) {
    TEST_ASSERT(0 == 1, "expected zero!");
    TEST_RETURN  // fail
}
```

## Error

Use `TEST_THROW_ERROR` when a test cannot continue due to an unrecoverable error.

```c
const char *filename = "datafile";
FILE *fp = fopen(datafile, "rb");
if (!fp) {
    TEST_THROW_ERROR("unable to open %s for reading: %s", filename, strerror(errno));
}
```

## Skip

```c
TEST_THROW_SKIP("Not implemented");
```

To conditionally skip a test, use an `if` statement

```c
if (!use_large_dataset) {
    TEST_THROW_SKIP("Test run not configured for large datasets");
}
```

# Fixtures

Setup and teardown fixtures are defined per-suite and run before, and after each unit test.

## Setup/Teardown

```c
static int setup(void) {
    puts("Inside setup");
    return 0;
}

static int teardown(void) {
    puts("Inside teardown");
    return 0;
}
```

## Register Fixtures

```c
TEST_SUITE("Fixtures") {
    TEST_SUITE_SET_FIXTURE_SETUP(setup);
    TEST_SUITE_SET_FIXTURE_TEARDOWN(teardown);
    unit_test tests[] = {
        TEST_UNIT_REPR(mytest)
    }
    TEST_SUITE_RUN(tests);
    TEST_SUITE_RETURN
}
```

### Output

```
 SUITE Fixtures...
  UNIT test__mytest... PASSED
  Inside setup
  Inside teardown
```

# Marking

Use a `TEST_MARK` to document test execution. This aids debugging unexpected crashes (`SIGSEGV`, etc.).

```c
TEST_UNIT(mytest) {
    TEST_MARK("check strstr");
    TEST_ASSERT(strstr("this substring is present", "substring"), "substring should match");
    TEST_MARK("check strncmp");
    TEST_ASSERT(strncmp(badptr + 3, "string", 6) == 0, "always segfaults")
    TEST_MARK("check strpbrk");
    TEST_ASSERT(strpbrk("this has spaces", " ") != NULL, "string should contain spaces");
    TEST_RETURN
}
```

## Output

```
SUITE mysuite...
  UNIT test__mytest... 
  MARK check strstr
  MARK check strncmp
segmentation fault
```

Breadcrumbs from the crash may be available in the unit test's output log.

```
less -R unitest_log/mysuite/test_mytest.c/test__mytest.log
```

# Checkpoints

Each test has a local error counter. Use `TEST_CHECKPOINT` to detect problems in a given block of `TEST_ASSERT` calls.

## Example

```c
TEST_UNIT(multiple_asserts_checkpoint) {
    TEST_CHECKPOINT("testing closely related variables");
    TEST_ASSERT('a' == 'a' "should be 'a'");
    TEST_ASSERT('b' == 'z' "should be 'b'");  // will fail
    TEST_ASSERT('c' == 'c', "should be 'c'");
    
    if (TEST_CHECKPOINT_FAILED()) {
        // assuming the failed assertion introduces an unrecoverable state
        TEST_THROW_ERROR("cannot continue to next step");
    }

    TEST_ASSERT(important_call('b') == 0, "this crashes the program otherwise");
    
    TEST_RETURN
}
```

## Output

```
SUITE mysuite...
  UNIT test__multiple_asserts_checkpoint... ERROR
  CHECKPOINT #0 - testing closely related variables 
  ASSERTION_FAILED 'b' == 'z' BECAUSE should be 'b'
  ERROR cannot continue to next step
```

# Testing Captured Output

Use `TEST_REDIRECT_SAVE` to save the test log to a file.

```c
static void a_print_only_function(void) {
    printf("Did it work?\n");
    printf("Yes.");
}

TEST_UNIT(a_print_only_function) {
    a_print_only_function();

    const char *filename = "capture.txt";
    if (TEST_REDIRECT_SAVE(&TEST_LOCAL_REDIRECT, filename)) {
        TEST_THROW_ERROR("unable to save captured output to %s", filename);
    }

    size_t line_count = 0;
    char **lines = TEST_FILE_AS_ARRAY(filename, &line_count);
    if (!lines) {
        TEST_THROW_ERROR("unable to read %s", filename);
    }

    const int expected_lines = 2;
    TEST_CHECKPOINT_ERROR_COUNT_SET();
    TEST_ASSERT(line_count == expected_lines, "expected %zu lines, got %zu", expected_lines, line_count);
    if (TEST_CHECKPOINT_ERROR_COUNT_DELTA) {
        TEST_MARK("debug output");
        for (size_t i = 0; i < line_count; i++) {
            printf("line %zu:\n", i + 1);
            hexdump(lines[i], strlen(lines[i]) + 1);
        }
    } else {
        // Tests that can only run if we have enough data
        TEST_CHECKPOINT_ERROR_COUNT_CLEAR();
        TEST_ASSERT(strcmp(lines[0], "Did it work?\n") == 0, "should have worked");
        TEST_ASSERT(strcmp(lines[1], "Yes.") == 0, "should have worked");
    }

    TEST_ARRAY_FREE(lines, line_count);
}
```

## Output

### Pass

```
 SUITE test_mytest...
  UNIT test__a_print_only_function... PASSED
```

### Fail

And if a failure is forced by commenting out `printf("Yes.");` in `a_print_only_function`...

```
SUITE test_mytest...
  UNIT test__a_print_only_function... FAILED
       ASSERTION FAILED line_count == expected_lines BECAUSE not enough data in output
       MARK debug output
       line 1:
       00000000  44 69 64 20 69 74 20 77  6f 72 6b 3f 0a 00        Did it work?..
```

# Helper Functions

## hexdump

Prints bytes to stdout in a traditional hexdump format.

```
[ offset ] [               16 hex values                 ] [ 16 characters  ]
 00000000   00 01 02 03 04 05 06  07 08 09 0a 0b 0c 0e 0f   ................
```

### Signature

```c
void hexdump(const char *s, size_t maxlen);
```

### Example

```c
#include "unittest.h"

TEST_UNIT(my_test){
    printf("# Character arrays\n");
    const char *s = "This is useful for debugging, I promise.";
    hexdump(s, strlen(s) + 1);
    
    struct data {
        char key[9];
        char value[71];
    };
    const struct data x = {
        .key = "kw_here", 
        .value = "this is a value / and this is a rather short comment"
    };
    printf("# Structures\n");
    hexdump(&x, sizeof(x));
    return 0;
}
```

### Output

```
# Character arrays
00000000  54 68 69 73 20 69 73 20  75 73 65 66 75 6c 20 66  This is useful f
00000010  6f 72 20 64 65 62 75 67  67 69 6e 67 2c 20 49 20  or debugging, I 
00000020  70 72 6f 6d 69 73 65 2e  00                       promise..
# Structures
00000000  6b 77 5f 68 65 72 65 00  00 74 68 69 73 20 69 73  kw_here..this is
00000010  20 61 20 76 61 6c 75 65  20 2f 20 61 6e 64 20 74   a value / and t
00000020  68 69 73 20 69 73 20 61  20 72 61 74 68 65 72 20  his is a rather 
00000030  73 68 6f 72 74 20 63 6f  6d 6d 65 6e 74 00 00 00  short comment...
00000040  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  ................
```

## TEST_FILE_AS_ARRAY

Read a file, storing each line in an array. The number of lines in the files is returned
in `lines_count`.

### Signature

```c
static inline char **TEST_FILE_AS_ARRAY(const char *filename, size_t *lines_count);
```

### Example

**file.txt**

```
I am real data.
I promise.
Really.
```

**example.c**

```c
#include "unittest.h"

TEST_UNIT(file_as_array) {
    const char *filename = "file.txt";
    
    printf("reading %s\n", filename);
    size_t num_lines = 0;
    char **lines = TEST_FILE_AS_ARRAY(filename, &num_lines);
    if (!lines) {
        TEST_THROW_ERROR("unable to read %s: %s", filename, strerror(errno));
    }
    if (!num_lines) {
        TEST_FREE_ARRAY(lines, num_lines);
        TEST_THROW_FAIL("%s is an empty file", filename);
    }

    TEST_CHECKPOINT("checking contents of %s", filename)
    //
    // TEST_ASSERT(s) here
    //
    
    // If any assertions fail after the checkpoint, display the lines
    if (TEST_CHECKPOINT_FAILED()) {
        for (size_t i = 0; i < num_lines; i++) {
            char *line = lines[i];
            printf("(line %zu) %s", i + 1, line);
        }
    }

    TEST_FREE_ARRAY(lines, num_lines);
    TEST_RETURN
}
```

### Output

```
(line 1) I am real data.

(line 2) I promise.

(line 3) Really.
```

## TEST_FILE_CONTAINS

## Signature

```c
static inline int TEST_FILE_CONTAINS(const char *filename, const char *pattern, const ssize_t range_start, const ssize_t range_end, char **result, size_t *result_lineno) {
```

## Example

```
TEST_UNIT(file_contains) {
    const char *filename = "file.txt";
    const char *pattern = "promise";
    TEST_ASSERT(
        TEST_FILE_CONTAINS(filename, pattern, -1, -1, &match, &match_lineno),
        "String not present in file"
    );
    if (match) {
        printf("match: %s\n", match ? "NO MATCH");
        printf("match line number: %zu\n", match_lineno);
    }
    TEST_RETURN
}
```