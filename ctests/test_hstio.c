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

TEST_BEGIN(macro_Pix) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(macro_PixColumnMajor) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(macro_PPixColumnMajor) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(macro_DQPix) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(macro_DQSetPix) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyDataSection) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_hstio_err) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_hstio_errmsg) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_push_hstioerr) {
    TEST_MARK("check state after pushing handler onto stack");
    for (int i = 0; i < handlers_max; i++) {
        const int stack_level = push_hstioerr(hstio_error_handler);
        TEST_ASSERT(stack_level == i + 1, "push failed, got %d", stack_level);
    }
    TEST_RETURN
}

TEST_BEGIN(fn_pop_hstioerr) {
    TEST_MARK("check state after poppping handler from stack");
    for (int i = handlers_max; i >= 0; i--) {
        const int stack_level = pop_hstioerr();
        TEST_ASSERT(stack_level == i - 1, "pop failed, got %d", stack_level);
    }
    TEST_RETURN
}

TEST_BEGIN(fn_clear_hstioerr) {
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

TEST_BEGIN(fn_error) {
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

TEST_BEGIN(fn_openSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_closeSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_fcloseNull) {
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

TEST_BEGIN(fn_fcloseWithStatus) {
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

TEST_BEGIN(fn_ckNewFile) {
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

TEST_BEGIN(fn_getSci) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSci) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getErr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putErr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getDQ) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putDQ) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSmpl) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSmpl) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getIntg) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putIntg) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSciSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putErrSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putDQSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSmplSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putIntgSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSciHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getErrHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getDQHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSciLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getErrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getDQLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getFloatHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putFloatHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getShortHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putShortHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putFloatHDSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putShortHDSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getFloatHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getShortHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleGroupSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleNicmosGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleNicmosGroupSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getMultiGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putMultiGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getMultiNicmosGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putMultiNicmosGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_makePrimaryArrayHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_makeImageExtHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getKeyB) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getKeyI) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getKeyF) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getKeyD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getKeyS) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putKeyB) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putKeyI) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putKeyF) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putKeyD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putKeyS) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyB) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyI) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyF) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyS) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyOrAddAsHistKeyBool) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyOrAddAsHistKeyInt) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyOrAddAsHistKeyFloat) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyOrAddAsHistKeyDouble) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateKeyOrAddAsHistKeyStr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getKwName) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getKwComm) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getKwType) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getBoolKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getIntKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getFloatKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getDoubleKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getStringKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putKwName) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putKwComm) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putBoolKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putIntKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putFloatKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putDoubleKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putStringKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_addBoolKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_addIntKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_addFloatKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_addDoubleKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_addStringKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_addSpacesKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_addCommentKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_addHistoryKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_insertBoolKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_insertIntKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_insertFloatKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_insertDoubleKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_insertStringKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_insertSpacesKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_insertCommentKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_insertHistoryKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_addFitsCard) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_insertFitsCard) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_delKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_delAllKw) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_openInputImage) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_openOutputImage) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_openUpdateImage) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_closeImage) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getFilename) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getExtname) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getExtver) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getNaxis1) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getNaxis2) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getType) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getHeader) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putHeader) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putFloatSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putShortSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_updateWCS) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_swapFloatStorageOrder) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_swapShortStorageOrder) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeFloatLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeShortLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocSciLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocErrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocDQLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_reallocHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initFloatHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocFloatHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyFloatHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeFloatHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initShortHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocShortHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyShortHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeShortHdrData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initFloatHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocFloatHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeFloatHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initShortHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocShortHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeShortHdrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocSingleGroupHeader) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocSingleGroupExts) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_setStorageOrder) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copySingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyOffsetSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyOffsetFloatData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyOffsetShortData) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_initSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_allocSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_freeSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getNumHDUs) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_findTotalNumberOfImsets) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_findTotalNumberOfHDUSets) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}
TEST_BEGIN(fn_FloatNAN) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_DoubleNAN) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_get_numeric) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_SUITE_BEGIN(__FILE__) {
    const testfunc tests[] = {
        test_macro_Pix,
        test_macro_PixColumnMajor,
        test_macro_PPixColumnMajor,
        test_macro_DQPix,
        test_macro_DQSetPix,
        test_fn_copyDataSection,
        test_fn_hstio_errmsg,
        test_fn_hstio_err,
        test_fn_push_hstioerr,
        test_fn_pop_hstioerr,
        test_fn_clear_hstioerr,
        test_fn_error,
        test_fn_openSingleGroupLine,
        test_fn_closeSingleGroupLine,
        test_fn_fcloseNull,
        test_fn_fcloseWithStatus,
        test_fn_ckNewFile,
        test_fn_getSci,
        test_fn_putSci,
        test_fn_getErr,
        test_fn_putErr,
        test_fn_getDQ,
        test_fn_putDQ,
        test_fn_getSmpl,
        test_fn_putSmpl,
        test_fn_getIntg,
        test_fn_putIntg,
        test_fn_putSciSect,
        test_fn_putErrSect,
        test_fn_putDQSect,
        test_fn_putSmplSect,
        test_fn_putIntgSect,
        test_fn_getSciHdr,
        test_fn_getErrHdr,
        test_fn_getDQHdr,
        test_fn_getSciLine,
        test_fn_getErrLine,
        test_fn_getDQLine,
        test_fn_getFloatHD,
        test_fn_putFloatHD,
        test_fn_getShortHD,
        test_fn_putShortHD,
        test_fn_putFloatHDSect,
        test_fn_putShortHDSect,
        test_fn_getFloatHdr,
        test_fn_getShortHdr,
        test_fn_getSingleGroup,
        test_fn_getSingleGroupLine,
        test_fn_putSingleGroupHdr,
        test_fn_putSingleGroup,
        test_fn_putSingleGroupSect,
        test_fn_getSingleNicmosGroup,
        test_fn_putSingleNicmosGroupHdr,
        test_fn_putSingleNicmosGroup,
        test_fn_putSingleNicmosGroupSect,
        test_fn_getMultiGroupHdr,
        test_fn_getMultiGroup,
        test_fn_putMultiGroupHdr,
        test_fn_putMultiGroup,
        test_fn_getMultiNicmosGroupHdr,
        test_fn_getMultiNicmosGroup,
        test_fn_putMultiNicmosGroupHdr,
        test_fn_putMultiNicmosGroup,
        test_fn_makePrimaryArrayHdr,
        test_fn_makeImageExtHdr,
        test_fn_getKeyB,
        test_fn_getKeyI,
        test_fn_getKeyF,
        test_fn_getKeyD,
        test_fn_getKeyS,
        test_fn_putKeyB,
        test_fn_putKeyI,
        test_fn_putKeyF,
        test_fn_putKeyD,
        test_fn_putKeyS,
        test_fn_updateKeyB,
        test_fn_updateKeyI,
        test_fn_updateKeyF,
        test_fn_updateKeyD,
        test_fn_updateKeyS,
        test_fn_updateKeyOrAddAsHistKeyBool,
        test_fn_updateKeyOrAddAsHistKeyInt,
        test_fn_updateKeyOrAddAsHistKeyFloat,
        test_fn_updateKeyOrAddAsHistKeyDouble,
        test_fn_updateKeyOrAddAsHistKeyStr,
        test_fn_getKwName,
        test_fn_getKwComm,
        test_fn_getKwType,
        test_fn_getBoolKw,
        test_fn_getIntKw,
        test_fn_getFloatKw,
        test_fn_getDoubleKw,
        test_fn_getStringKw,
        test_fn_putKwName,
        test_fn_putKwComm,
        test_fn_putBoolKw,
        test_fn_putIntKw,
        test_fn_putFloatKw,
        test_fn_putDoubleKw,
        test_fn_putStringKw,
        test_fn_addBoolKw,
        test_fn_addIntKw,
        test_fn_addFloatKw,
        test_fn_addDoubleKw,
        test_fn_addStringKw,
        test_fn_addSpacesKw,
        test_fn_addCommentKw,
        test_fn_addHistoryKw,
        test_fn_insertBoolKw,
        test_fn_insertIntKw,
        test_fn_insertFloatKw,
        test_fn_insertDoubleKw,
        test_fn_insertStringKw,
        test_fn_insertSpacesKw,
        test_fn_insertCommentKw,
        test_fn_insertHistoryKw,
        test_fn_addFitsCard,
        test_fn_insertFitsCard,
        test_fn_delKw,
        test_fn_delAllKw,
        test_fn_openInputImage,
        test_fn_openOutputImage,
        test_fn_openUpdateImage,
        test_fn_closeImage,
        test_fn_getFilename,
        test_fn_getExtname,
        test_fn_getExtver,
        test_fn_getNaxis1,
        test_fn_getNaxis2,
        test_fn_getType,
        test_fn_getHeader,
        test_fn_putHeader,
        test_fn_getFloatData,
        test_fn_putFloatData,
        test_fn_getShortData,
        test_fn_putShortData,
        test_fn_putFloatSect,
        test_fn_putShortSect,
        test_fn_getFloatLine,
        test_fn_putFloatLine,
        test_fn_getShortLine,
        test_fn_putShortLine,
        test_fn_updateWCS,
        test_fn_initFloatData,
        test_fn_allocFloatData,
        test_fn_freeFloatData,
        test_fn_swapFloatStorageOrder,
        test_fn_initShortData,
        test_fn_allocShortData,
        test_fn_freeShortData,
        test_fn_swapShortStorageOrder,
        test_fn_initFloatLine,
        test_fn_allocFloatLine,
        test_fn_freeFloatLine,
        test_fn_initShortLine,
        test_fn_allocShortLine,
        test_fn_freeShortLine,
        test_fn_allocSciLine,
        test_fn_allocErrLine,
        test_fn_allocDQLine,
        test_fn_initHdr,
        test_fn_allocHdr,
        test_fn_reallocHdr,
        test_fn_freeHdr,
        test_fn_copyHdr,
        test_fn_initFloatHdrData,
        test_fn_allocFloatHdrData,
        test_fn_copyFloatHdrData,
        test_fn_freeFloatHdrData,
        test_fn_initShortHdrData,
        test_fn_allocShortHdrData,
        test_fn_copyShortHdrData,
        test_fn_freeShortHdrData,
        test_fn_initFloatHdrLine,
        test_fn_allocFloatHdrLine,
        test_fn_freeFloatHdrLine,
        test_fn_initShortHdrLine,
        test_fn_allocShortHdrLine,
        test_fn_freeShortHdrLine,
        test_fn_initSingleGroup,
        test_fn_allocSingleGroup,
        test_fn_allocSingleGroupHeader,
        test_fn_allocSingleGroupExts,
        test_fn_freeSingleGroup,
        test_fn_setStorageOrder,
        test_fn_copySingleGroup,
        test_fn_copyOffsetSingleGroup,
        test_fn_initMultiGroup,
        test_fn_allocMultiGroup,
        test_fn_freeMultiGroup,
        test_fn_copyFloatData,
        test_fn_copyOffsetFloatData,
        test_fn_copyShortData,
        test_fn_copyOffsetShortData,
        test_fn_initSingleNicmosGroup,
        test_fn_allocSingleNicmosGroup,
        test_fn_freeSingleNicmosGroup,
        test_fn_initMultiNicmosGroup,
        test_fn_allocMultiNicmosGroup,
        test_fn_freeMultiNicmosGroup,
        test_fn_initSingleGroupLine,
        test_fn_allocSingleGroupLine,
        test_fn_freeSingleGroupLine,
        test_fn_getNumHDUs,
        test_fn_findTotalNumberOfImsets,
        test_fn_findTotalNumberOfHDUSets,
        test_fn_FloatNAN,
        test_fn_DoubleNAN,
        test_fn_get_numeric,
    };
    TEST_SUITE_RUN(tests);
    TEST_SUITE_RETURN
}
