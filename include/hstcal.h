#ifndef HSTCAL_INCL
#define HSTCAL_INCL

// Note: Use of CHAR_FNAME_LENGTH & CHAR_LINE_LENGTH are interchanged throughout and should therefore be identical (#164).
// However, for now, keep them separated such that they can be altered back (replace all) to their original IRAF equivalents set in include/c_iraf.h:
// # define IRAF_SZ_LINE    1023
// # define IRAF_SZ_FNAME    255
//
// Note: IRAF_SZ_LINE = 1023 relates to the maximum filename size accepted by cfitsio.
#define CHAR_FNAME_LENGTH 255
#define CHAR_LINE_LENGTH 255
#define MSG_BUFF_LENGTH CHAR_LINE_LENGTH + 1

/* Standard string buffer for use in messages */
extern char MsgText[MSG_BUFF_LENGTH];

#endif
