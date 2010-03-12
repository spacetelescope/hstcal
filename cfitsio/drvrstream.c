/*  This file, drvrstream.c contains driver routines for */
/*  writing FITS files to the stdin or stdout streams.         */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

#include <stdlib.h>
#include "fitsio2.h"

/*--------------------------------------------------------------------------*/
int stream_open(char *filename, int rwmode, int *handle)
{
    /*
        read from stdin
    */
    rwmode = (int) filename;  /* suppress unused parameter compiler warning */
    *handle = 1;     /*  1 = stdin */   

    return(0);
}
/*--------------------------------------------------------------------------*/
int stream_create(char *filename, int *handle)
{
    /*
        write to stdout
    */

    *handle = (int) filename;  /* suppress unused parameter compiler warning */
    *handle = 2;         /*  2 = stdout */       
    return(0);
}
/*--------------------------------------------------------------------------*/
int stream_size(int handle, LONGLONG *filesize)
/*
  return the size of the file in bytes
*/
{
    handle = 0;  /* suppress unused parameter compiler warning */
    
    /* this operation is not supported in a stream; return large value */
    *filesize = LONG_MAX;
    return(0);
}
/*--------------------------------------------------------------------------*/
int stream_close(int handle)
/*
     don't have to close stdin or stdout 
*/
{
    handle = 0;  /* suppress unused parameter compiler warning */
    
    return(0);
}
/*--------------------------------------------------------------------------*/
int stream_flush(int handle)
/*
  flush the file
*/
{
    if (handle == 2)
       fflush(stdout);  

    return(0);
}
/*--------------------------------------------------------------------------*/
int stream_seek(int handle, LONGLONG offset)
   /* 
      seeking is not allowed in a stream
   */
{
    offset = handle;  /* suppress unused parameter compiler warning */
    return(1);
}
/*--------------------------------------------------------------------------*/
int stream_read(int hdl, void *buffer, long nbytes)
/*
     reading from stdin stream 
*/

{
    long nread;
    
    if (hdl != 1)
       return(1);  /* can only read from stdin */

    nread = (long) fread(buffer, 1, nbytes, stdin);

    if (nread != nbytes)
    {
/*        return(READ_ERROR); */
        return(END_OF_FILE);
    }

    return(0);
}
/*--------------------------------------------------------------------------*/
int stream_write(int hdl, void *buffer, long nbytes)
/*
  write bytes at the current position in the file
*/
{
    if (hdl != 2)
       return(1);  /* can only write to stdout */

    if((long) fwrite(buffer, 1, nbytes, stdout) != nbytes)
        return(WRITE_ERROR);

    return(0);
}


