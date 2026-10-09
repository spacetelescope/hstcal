#include <stddef.h>
#include "unittest.h"
#include "trlbuf.h"
#include "timestamp.h"

TEST_UNIT(timestamp) {
    const char *filename = "capture.txt";
    InitTrlBuf();
    TimeStamp("Generic Conversion Complete", "N2RUSG54B");
    TEST_REDIRECT_SAVE(&TEST_LOCAL_REDIRECT, filename);
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        TEST_THROW_ERROR("Unable to open log capture file: %s", filename);
    }
    char line[255] = {0};
    if (!fgets(line, sizeof(line), fp)) {
        TEST_THROW_ERROR("fgets() failed");
    }
    // remove linefeed from output
    if (line[strlen(line) - 1] == '\n') {
        line[strlen(line) - 1] = '\0';
    }
    char datestamp[SZ_TIMESTRING] = {0};
    char message[SZ_TIMESTRING] = {0};

    const char *fmt = "%13s%[^*]s";
    const int len = sscanf(line, fmt, datestamp, message);
    if (len < 0) {
        TEST_THROW_FAIL("message '%s' does not match format: %s", fmt);
    }
    TEST_ASSERT(len <= SZ_TIMESTRING, "Message too long, got %zu", len);
    TEST_ASSERT(strlen(datestamp) == 13, "Datestamp length mismatch, should be 13, got %zu", strlen(datestamp));
    TEST_ASSERT(strlen(message) == 67, "Message length mismatch, should be 67, got %zu", strlen(message));
    TEST_ASSERT(strlen(line) == 80, "Message should be padded to maximum length of %zu, got %zu", SZ_TIMESTRING, strlen(line));
    remove(filename);

    CloseTrlBuf(&trlbuf);
    TEST_RETURN
}

TEST_SUITE(__FILE__) {
    const unit_test tests[] = {
        TEST_UNIT_REPR(timestamp),
    };

    TEST_SUITE_RUN(tests)
    TEST_SUITE_RETURN
}