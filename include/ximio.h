#ifndef INCL_XIMIO_H
#define INCL_XIMIO_H

# include <c_iraf.h>

IRAFPointer c_imtopen(char *pattern);

int c_imtlen(IRAFPointer imt);

void c_imtrew(IRAFPointer imt);

int c_imtgetim(IRAFPointer imt, char *outstr, int maxch);

void c_imtclose(IRAFPointer imt);

#endif /* INCL_XIMIO_H */
