#include "unittest.h"
#include "hstio.h"

// Load symbols that should be exported but aren't
void error(HSTIOError e, char *str);
void ioerr(HSTIOError e, IODescPtr x_, int status);

static void hstio_error_handler() {
#ifdef DEBUG
    TEST_MSG(stderr, NULL, __func__, "called with status %d:\nmessage: '%s'\nmessage length: %zu\n", hstio_err(),
        hstio_errmsg(), strlen(hstio_errmsg()));
#endif
}
static const int handlers_max = 2; // up to 32 handlers, see hstio.h

TEST_UNIT(macro_Pix) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(macro_PixColumnMajor) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(macro_PPixColumnMajor) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(macro_DQPix) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(macro_DQSetPix) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copyDataSection) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(hstio_err) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(hstio_errmsg) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(push_hstioerr) {
    TEST_MARK("check state after pushing handler onto stack");
    for (int i = 0; i < handlers_max; i++) {
        const int stack_level = push_hstioerr(hstio_error_handler);
        TEST_ASSERT(stack_level == i + 1, "push failed, got %d", stack_level);
    }
    TEST_RETURN
}

TEST_UNIT(pop_hstioerr) {
    TEST_MARK("check state after poppping handler from stack");
    for (int i = handlers_max; i >= 0; i--) {
        const int stack_level = pop_hstioerr();
        TEST_ASSERT(stack_level == i - 1, "pop failed, got %d", stack_level);
    }
    TEST_RETURN
}

TEST_UNIT(clear_hstioerr) {
    TEST_MARK("check error state after clear");
    clear_hstioerr();
    const int err = hstio_err();
    const char *msg = hstio_errmsg();
    TEST_MARK("check error code is zero");
    TEST_ASSERT(err == 0, "expected error code to be zero, got %d\n", err);
    TEST_MARK("check error message is zero-length");
    TEST_ASSERT(strlen(msg) == 0, "expected message length to be zero, got %d\n", strlen(msg));
    TEST_RETURN
}

TEST_UNIT(error) {
    size_t error_len = 0;
    push_hstioerr(hstio_error_handler);

    TEST_MARK("check contents of error_msg buffer");
    for (size_t i = 0; i <= BADREMOVE; i++) {
        error(i, "canary");
        const char *error_msg = hstio_errmsg();
        TEST_ASSERT(strstr(error_msg, "canary") != NULL, "Message does not contain expected value! error_msg = '%s'\n",
            error_msg);
    }

    TEST_MARK("check unhandled error");
    error(9999, "this is an unhandled error");
    const char *error_msg = hstio_errmsg();
    TEST_ASSERT(strstr(error_msg, "HSTIOError") != NULL, "Failed to trap unhandled error code! error_msg = '%s'\n",
        error_msg);

    TEST_MARK("check very long message truncation");
    char huge_error_message[4096] = {0};
    memset(huge_error_message, '?', sizeof(huge_error_message) - 1);
    error(BADNAME, huge_error_message);
    error_msg = hstio_errmsg();
    error_len = strlen(error_msg);
    const char *expected = "Keyword name";
    TEST_ASSERT(strstr(error_msg, expected) != NULL,
        "Message does not contain expected value, '%s'! error_msg = '%s'\n", expected, error_msg);

    TEST_MARK("HSTOK with an empty error string produces an empty string");
    error(HSTOK, "");
    error_msg = hstio_errmsg();
    error_len = strlen(error_msg);
    TEST_ASSERT(error_len == 0, "Message should be empty! error_msg = '%s'\n", error_msg);

    pop_hstioerr();
    TEST_RETURN
}

TEST_UNIT(openSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(closeSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(fcloseNull) {
    FILE *fp = fopen("/dev/null", "w");
    if (!fp) {
        TEST_THROW_ERROR("failed to open /dev/null for writing");
    }
    TEST_MARK("check closing a file handle");
    TEST_ASSERT(fcloseNull(fp) == 0, "%s", "stream did not close\n");

    TEST_MARK("check passing NULL file handle returns zero");
    TEST_ASSERT(fcloseNull(NULL) == 0, "%s", "function returned non-zero on NULL file handle\n");
    TEST_RETURN
}

TEST_UNIT(fcloseWithStatus) {
    FILE *fp = fopen("/dev/null", "w");
    if (!fp) {
        TEST_THROW_ERROR("%s", "Error opening /dev/null");
    }
    TEST_MARK("check closing a file handle");
    TEST_ASSERT(fcloseWithStatus(&fp) == 0, "%s", "stream did not close\n");

    TEST_MARK("check file handle is NULL");
    TEST_ASSERT(fp == NULL, "%s", "stream is not NULL\n");

    // fcloseWithStatus returns IO_ERROR when fclose() fails. Intentionally
    // forcing fclose() to fail is difficult. Accessing a stream after you've
    // called fclose() on it is purely undefined behavior.
    TEST_RETURN
}

TEST_UNIT(ckNewFile) {
    TEST_MARK("generating test file");
    const char *filename = "testnewfile";
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        TEST_THROW_ERROR("Error opening %s\n", filename);
    }
    fprintf(fp, "hello world from %s\n", filename);
    fclose(fp);

    int exists = 0;
    exists = ckNewFile((char *) filename);
    TEST_MARK("check non-zero return value if file exists");
    TEST_ASSERT(exists == 1, "existence check failed, returned %d\n", exists);

    TEST_MARK("removing test file");
    remove(filename);
    exists = ckNewFile((char *) filename);
    TEST_MARK("check zero return value if file does not exist");
    TEST_ASSERT(exists == 0, "existence check failed, returned %d\n", exists);
    TEST_RETURN
}

TEST_UNIT(getSci) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSci) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getErr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putErr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getDQ) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putDQ) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getSmpl) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSmpl) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getIntg) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putIntg) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSciSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putErrSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putDQSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSmplSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putIntgSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getSciHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getErrHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getDQHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getSciLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getErrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getDQLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getFloatHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putFloatHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getShortHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putShortHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putFloatHDSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putShortHDSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getFloatHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getShortHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSingleGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSingleGroupSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSingleNicmosGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putSingleNicmosGroupSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getMultiGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putMultiGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getMultiNicmosGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putMultiNicmosGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(makePrimaryArrayHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(makeImageExtHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getKeyB) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getKeyI) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getKeyF) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getKeyD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getKeyS) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putKeyB) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putKeyI) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putKeyF) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putKeyD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putKeyS) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyB) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyI) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyF) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyS) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyOrAddAsHistKeyBool) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyOrAddAsHistKeyInt) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyOrAddAsHistKeyFloat) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyOrAddAsHistKeyDouble) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateKeyOrAddAsHistKeyStr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getKwName) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getKwComm) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getKwType) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getBoolKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getIntKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getFloatKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getDoubleKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getStringKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putKwName) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putKwComm) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putBoolKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putIntKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putFloatKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putDoubleKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putStringKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(addBoolKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(addIntKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(addFloatKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(addDoubleKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(addStringKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(addSpacesKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(addCommentKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(addHistoryKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(insertBoolKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(insertIntKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(insertFloatKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(insertDoubleKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(insertStringKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(insertSpacesKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(insertCommentKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(insertHistoryKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(addFitsCard) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(insertFitsCard) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(delKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(delAllKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(openInputImage) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(openOutputImage) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(openUpdateImage) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(closeImage) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getFilename) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getExtname) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getExtver) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getNaxis1) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getNaxis2) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getType) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getHeader) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putHeader) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putFloatSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putShortSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(putShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(updateWCS) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(swapFloatStorageOrder) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(swapShortStorageOrder) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocSciLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocErrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocDQLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(reallocHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copyHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initFloatHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocFloatHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copyFloatHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeFloatHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initShortHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocShortHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copyShortHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeShortHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initFloatHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocFloatHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeFloatHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initShortHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocShortHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeShortHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocSingleGroupHeader) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocSingleGroupExts) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(setStorageOrder) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copySingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copyOffsetSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copyFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copyOffsetFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copyShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(copyOffsetShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(initSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(allocSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(freeSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(getNumHDUs) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(findTotalNumberOfImsets) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(findTotalNumberOfHDUSets) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}
TEST_UNIT(FloatNAN) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(DoubleNAN) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_UNIT(get_numeric) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_SUITE(__FILE__) {
    const unit_test tests[] = {
        TEST_UNIT_REPR(macro_Pix),
        TEST_UNIT_REPR(macro_PixColumnMajor),
        TEST_UNIT_REPR(macro_PPixColumnMajor),
        TEST_UNIT_REPR(macro_DQPix),
        TEST_UNIT_REPR(macro_DQSetPix),
        TEST_UNIT_REPR(copyDataSection),
        TEST_UNIT_REPR(hstio_errmsg),
        TEST_UNIT_REPR(hstio_err),
        TEST_UNIT_REPR(push_hstioerr),
        TEST_UNIT_REPR(pop_hstioerr),
        TEST_UNIT_REPR(clear_hstioerr),
        TEST_UNIT_REPR(error),
        TEST_UNIT_REPR(openSingleGroupLine),
        TEST_UNIT_REPR(closeSingleGroupLine),
        TEST_UNIT_REPR(fcloseNull),
        TEST_UNIT_REPR(fcloseWithStatus),
        TEST_UNIT_REPR(ckNewFile),
        TEST_UNIT_REPR(getSci),
        TEST_UNIT_REPR(putSci),
        TEST_UNIT_REPR(getErr),
        TEST_UNIT_REPR(putErr),
        TEST_UNIT_REPR(getDQ),
        TEST_UNIT_REPR(putDQ),
        TEST_UNIT_REPR(getSmpl),
        TEST_UNIT_REPR(putSmpl),
        TEST_UNIT_REPR(getIntg),
        TEST_UNIT_REPR(putIntg),
        TEST_UNIT_REPR(putSciSect),
        TEST_UNIT_REPR(putErrSect),
        TEST_UNIT_REPR(putDQSect),
        TEST_UNIT_REPR(putSmplSect),
        TEST_UNIT_REPR(putIntgSect),
        TEST_UNIT_REPR(getSciHdr),
        TEST_UNIT_REPR(getErrHdr),
        TEST_UNIT_REPR(getDQHdr),
        TEST_UNIT_REPR(getSciLine),
        TEST_UNIT_REPR(getErrLine),
        TEST_UNIT_REPR(getDQLine),
        TEST_UNIT_REPR(getFloatHD),
        TEST_UNIT_REPR(putFloatHD),
        TEST_UNIT_REPR(getShortHD),
        TEST_UNIT_REPR(putShortHD),
        TEST_UNIT_REPR(putFloatHDSect),
        TEST_UNIT_REPR(putShortHDSect),
        TEST_UNIT_REPR(getFloatHdr),
        TEST_UNIT_REPR(getShortHdr),
        TEST_UNIT_REPR(getSingleGroup),
        TEST_UNIT_REPR(getSingleGroupLine),
        TEST_UNIT_REPR(putSingleGroupHdr),
        TEST_UNIT_REPR(putSingleGroup),
        TEST_UNIT_REPR(putSingleGroupSect),
        TEST_UNIT_REPR(getSingleNicmosGroup),
        TEST_UNIT_REPR(putSingleNicmosGroupHdr),
        TEST_UNIT_REPR(putSingleNicmosGroup),
        TEST_UNIT_REPR(putSingleNicmosGroupSect),
        TEST_UNIT_REPR(getMultiGroupHdr),
        TEST_UNIT_REPR(getMultiGroup),
        TEST_UNIT_REPR(putMultiGroupHdr),
        TEST_UNIT_REPR(putMultiGroup),
        TEST_UNIT_REPR(getMultiNicmosGroupHdr),
        TEST_UNIT_REPR(getMultiNicmosGroup),
        TEST_UNIT_REPR(putMultiNicmosGroupHdr),
        TEST_UNIT_REPR(putMultiNicmosGroup),
        TEST_UNIT_REPR(makePrimaryArrayHdr),
        TEST_UNIT_REPR(makeImageExtHdr),
        TEST_UNIT_REPR(getKeyB),
        TEST_UNIT_REPR(getKeyI),
        TEST_UNIT_REPR(getKeyF),
        TEST_UNIT_REPR(getKeyD),
        TEST_UNIT_REPR(getKeyS),
        TEST_UNIT_REPR(putKeyB),
        TEST_UNIT_REPR(putKeyI),
        TEST_UNIT_REPR(putKeyF),
        TEST_UNIT_REPR(putKeyD),
        TEST_UNIT_REPR(putKeyS),
        TEST_UNIT_REPR(updateKeyB),
        TEST_UNIT_REPR(updateKeyI),
        TEST_UNIT_REPR(updateKeyF),
        TEST_UNIT_REPR(updateKeyD),
        TEST_UNIT_REPR(updateKeyS),
        TEST_UNIT_REPR(updateKeyOrAddAsHistKeyBool),
        TEST_UNIT_REPR(updateKeyOrAddAsHistKeyInt),
        TEST_UNIT_REPR(updateKeyOrAddAsHistKeyFloat),
        TEST_UNIT_REPR(updateKeyOrAddAsHistKeyDouble),
        TEST_UNIT_REPR(updateKeyOrAddAsHistKeyStr),
        TEST_UNIT_REPR(getKwName),
        TEST_UNIT_REPR(getKwComm),
        TEST_UNIT_REPR(getKwType),
        TEST_UNIT_REPR(getBoolKw),
        TEST_UNIT_REPR(getIntKw),
        TEST_UNIT_REPR(getFloatKw),
        TEST_UNIT_REPR(getDoubleKw),
        TEST_UNIT_REPR(getStringKw),
        TEST_UNIT_REPR(putKwName),
        TEST_UNIT_REPR(putKwComm),
        TEST_UNIT_REPR(putBoolKw),
        TEST_UNIT_REPR(putIntKw),
        TEST_UNIT_REPR(putFloatKw),
        TEST_UNIT_REPR(putDoubleKw),
        TEST_UNIT_REPR(putStringKw),
        TEST_UNIT_REPR(addBoolKw),
        TEST_UNIT_REPR(addIntKw),
        TEST_UNIT_REPR(addFloatKw),
        TEST_UNIT_REPR(addDoubleKw),
        TEST_UNIT_REPR(addStringKw),
        TEST_UNIT_REPR(addSpacesKw),
        TEST_UNIT_REPR(addCommentKw),
        TEST_UNIT_REPR(addHistoryKw),
        TEST_UNIT_REPR(insertBoolKw),
        TEST_UNIT_REPR(insertIntKw),
        TEST_UNIT_REPR(insertFloatKw),
        TEST_UNIT_REPR(insertDoubleKw),
        TEST_UNIT_REPR(insertStringKw),
        TEST_UNIT_REPR(insertSpacesKw),
        TEST_UNIT_REPR(insertCommentKw),
        TEST_UNIT_REPR(insertHistoryKw),
        TEST_UNIT_REPR(addFitsCard),
        TEST_UNIT_REPR(insertFitsCard),
        TEST_UNIT_REPR(delKw),
        TEST_UNIT_REPR(delAllKw),
        TEST_UNIT_REPR(openInputImage),
        TEST_UNIT_REPR(openOutputImage),
        TEST_UNIT_REPR(openUpdateImage),
        TEST_UNIT_REPR(closeImage),
        TEST_UNIT_REPR(getFilename),
        TEST_UNIT_REPR(getExtname),
        TEST_UNIT_REPR(getExtver),
        TEST_UNIT_REPR(getNaxis1),
        TEST_UNIT_REPR(getNaxis2),
        TEST_UNIT_REPR(getType),
        TEST_UNIT_REPR(getHeader),
        TEST_UNIT_REPR(putHeader),
        TEST_UNIT_REPR(getFloatData),
        TEST_UNIT_REPR(putFloatData),
        TEST_UNIT_REPR(getShortData),
        TEST_UNIT_REPR(putShortData),
        TEST_UNIT_REPR(putFloatSect),
        TEST_UNIT_REPR(putShortSect),
        TEST_UNIT_REPR(getFloatLine),
        TEST_UNIT_REPR(putFloatLine),
        TEST_UNIT_REPR(getShortLine),
        TEST_UNIT_REPR(putShortLine),
        TEST_UNIT_REPR(updateWCS),
        TEST_UNIT_REPR(initFloatData),
        TEST_UNIT_REPR(allocFloatData),
        TEST_UNIT_REPR(freeFloatData),
        TEST_UNIT_REPR(swapFloatStorageOrder),
        TEST_UNIT_REPR(initShortData),
        TEST_UNIT_REPR(allocShortData),
        TEST_UNIT_REPR(freeShortData),
        TEST_UNIT_REPR(swapShortStorageOrder),
        TEST_UNIT_REPR(initFloatLine),
        TEST_UNIT_REPR(allocFloatLine),
        TEST_UNIT_REPR(freeFloatLine),
        TEST_UNIT_REPR(initShortLine),
        TEST_UNIT_REPR(allocShortLine),
        TEST_UNIT_REPR(freeShortLine),
        TEST_UNIT_REPR(allocSciLine),
        TEST_UNIT_REPR(allocErrLine),
        TEST_UNIT_REPR(allocDQLine),
        TEST_UNIT_REPR(initHdr),
        TEST_UNIT_REPR(allocHdr),
        TEST_UNIT_REPR(reallocHdr),
        TEST_UNIT_REPR(freeHdr),
        TEST_UNIT_REPR(copyHdr),
        TEST_UNIT_REPR(initFloatHdrData),
        TEST_UNIT_REPR(allocFloatHdrData),
        TEST_UNIT_REPR(copyFloatHdrData),
        TEST_UNIT_REPR(freeFloatHdrData),
        TEST_UNIT_REPR(initShortHdrData),
        TEST_UNIT_REPR(allocShortHdrData),
        TEST_UNIT_REPR(copyShortHdrData),
        TEST_UNIT_REPR(freeShortHdrData),
        TEST_UNIT_REPR(initFloatHdrLine),
        TEST_UNIT_REPR(allocFloatHdrLine),
        TEST_UNIT_REPR(freeFloatHdrLine),
        TEST_UNIT_REPR(initShortHdrLine),
        TEST_UNIT_REPR(allocShortHdrLine),
        TEST_UNIT_REPR(freeShortHdrLine),
        TEST_UNIT_REPR(initSingleGroup),
        TEST_UNIT_REPR(allocSingleGroup),
        TEST_UNIT_REPR(allocSingleGroupHeader),
        TEST_UNIT_REPR(allocSingleGroupExts),
        TEST_UNIT_REPR(freeSingleGroup),
        TEST_UNIT_REPR(setStorageOrder),
        TEST_UNIT_REPR(copySingleGroup),
        TEST_UNIT_REPR(copyOffsetSingleGroup),
        TEST_UNIT_REPR(initMultiGroup),
        TEST_UNIT_REPR(allocMultiGroup),
        TEST_UNIT_REPR(freeMultiGroup),
        TEST_UNIT_REPR(copyFloatData),
        TEST_UNIT_REPR(copyOffsetFloatData),
        TEST_UNIT_REPR(copyShortData),
        TEST_UNIT_REPR(copyOffsetShortData),
        TEST_UNIT_REPR(initSingleNicmosGroup),
        TEST_UNIT_REPR(allocSingleNicmosGroup),
        TEST_UNIT_REPR(freeSingleNicmosGroup),
        TEST_UNIT_REPR(initMultiNicmosGroup),
        TEST_UNIT_REPR(allocMultiNicmosGroup),
        TEST_UNIT_REPR(freeMultiNicmosGroup),
        TEST_UNIT_REPR(initSingleGroupLine),
        TEST_UNIT_REPR(allocSingleGroupLine),
        TEST_UNIT_REPR(freeSingleGroupLine),
        TEST_UNIT_REPR(getNumHDUs),
        TEST_UNIT_REPR(findTotalNumberOfImsets),
        TEST_UNIT_REPR(findTotalNumberOfHDUSets),
        TEST_UNIT_REPR(FloatNAN),
        TEST_UNIT_REPR(DoubleNAN),
        TEST_UNIT_REPR(get_numeric),
    };
    TEST_SUITE_RUN(tests);
    TEST_SUITE_RETURN
}
