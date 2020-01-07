/*
    The following function and structure definitions manage the
    output of CALACS trailer file comments, both to STDOUT and to
    the trailer file(s).  Code for using a temp file as an output
    trailer file have been commented out, but left in to insure that
    no unsuspecting problems arise from this change.  They should be
    deleted upon the next release assuming everything works as desired.

    The overwrite switch will be used to specify whether the entire
    input trailer file should be kept or whether old comments starting
    at TRL_PREFIX should be overwritten.  This would be set to YES for
    initial SINGLE exposure trailer files.

    It includes the following functions, provided with their calling
    sequences.  For the trailer files, we use the following terms:
        preface - comments common to more than 1 input file, copied into
            the beginning of the trailer file (after the GENERIC conversion
            comments).  This usually only applies to _RAW files.
        buffer  - trailer file comments for use in a single trailer file
        prefix  - string which marks the beginning of CALACS processing
            comments in the trailer files. (TRL_PREFIX)

    These control the buffer which contains the messages to be written out:
    int InitTrlBuf (void);
        - MUST BE CALLED PRIOR TO ANY OTHER TRAILER FILE FUNCTION!
        - initializes trailer buffer structure 'trlbuf'
    void SetTrlPrefaceMode (int use);
        - This function sets the switch to specify whether the
            preface should be output to the trailer file or
            not.  This is used when no CRREJ is done on a list
            of images, and you don't want the preface in the
            trailer files of the sub-products again.
    void SetTrlOverwriteMode (int owrite);
        - This function sets the overwrite switch.  Once set, all
            comments after the first CALACSBEG string will be
            overwritten.
    void SetTrlQuietMode (int quiet);
        - This function records the value of the command-line
            parameter QUIET into the structure.
    void InitTrlPreface (void);
        - This function will copy contents of the buffer into the preface.
    void ResetTrlPreface (void);
        - This function will reset the preface to a null string, so
            nothing more gets written out to any other trailer files.
    void CloseTrlBuf (void);
        - moves buffer from temp file to final output trailer file,
            then closes temp and output trailer files.

    These control the actual trailer files themselves:
    int InitTrlFile (char *inlist, char *output);
        - opens temporary trailer file based on the
            input trailer file(s) given in the file list 'inlist'
    int WriteTrlFile (void);
        - writes out temp trailer file to final trailer file, closes both
            and deletes the temp file (if all went OK).

    These create the messages to be written out:
    void trlmessage (char *message);
        - calls WrtTrlBuf and printfAndFlush
    void trlwarn (char *message);
        - calls trlmessage after building appropriate message
    void trlerror (char *message);
        - calls trlmessage after building appropriate message
    void trlopenerr (char *filename);
        - calls trlerror after starting to build appropriate message
    void trlreaderr (char *filename);
        - calls trlerror after building appropriate message
    void trlkwerr (char *keyword, char *filename);
        - calls trlerror after building appropriate message
    void trlfilerr (char *filename);
        - calls trlerror after building appropriate message

    The remainder of the functions are NOT USED outside this file:
    static void CatTrlFile(FILE *ip, FILE *op);
        - appends ENTIRE input trailer file (ip) to output trailer file (op)
    static int AppendTrlFile();
        - sets up trailer file so that all new comments start after preface
    static void AddTrlBuf (char *message);
        - adds a new message (line) to trailer comments buffer
    static void ResetTrlBuf (void);
        - clears out comments (sets to blank) already written to file
    static void WriteTrlBuf (char *message);
        - writes out comments in buffer to open trailer file
          then calls ResetTrlBuf to clear buffer
          -OR-
          calls AddTrlBuf to append comment to buffer.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "hstcal_memory.h"
#include "ximio.h"
#include "trlbuf.h"
#include "hstcalerr.h"
#include "hstcal.h"
#include "hstcalversion.h"
#include "hstio.h"

char MsgText[MSG_BUFF_LENGTH];
int status;
const unsigned initLength = 2; // Should always be > 0 to prevent issues with use of realloc().

/* Internal trailer file buffer routines */
static void ResetTrlBuf (void);
static void AddTrlBuf (const char *message);
static void WriteTrlBuf (const char *message);
static void CatTrlFile (FILE *ip, FILE *op);
static void CatTrlFile_NoEOF (FILE *ip, FILE *op);
static int AppendTrlFile(void);

struct TrlBuf trlbuf;

int InitTrlFile (char *inlist, char *output)
{
    /*
        This initialization function sets up the trailer file to be used
        as output for all the messages from the task it is called from.
        It will concatenate multiple input trailer files
        or copy intial portion of input trailer file into one output
        file, creating the single trailer file for the messages to be
        written to...
    */

/* Parameters:
char *inlist        i: list of input trailer filenames
char *output        i: full filename of output (final) trailer file
*/
    extern int status;

    IRAFPointer tpin = NULL;
    FILE * ip = NULL;
    FILE * tp = NULL;
    int n, td;
    char trldata[CHAR_LINE_LENGTH+1];
    static char uniq_outtemplate[] = "tmp_hstcal_XXXXXX";
    char uniq_outname[CHAR_FNAME_LENGTH+1];

    int unlink(const char *);

    trldata[0] = '\0';

    /* Copy name of output file to trlbuf */
    strcpy (trlbuf.trlfile, output);

    /*  Open the output file, for concatenating
        exposure trailers or sub-product trailers, or simply for
        appending new comments to end of existing trailer file.
    */
    if ((trlbuf.fp = fopen(trlbuf.trlfile,"a+")) == NULL) {
        trlopenerr(trlbuf.trlfile);
        return(status=INVALID_TEMP_FILE);
    }
    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
    /* open the input file template */
    tpin = c_imtopen (inlist);
    addPtr(&ptrReg, tpin, &c_imtclose);
    /* Only if we are building a product/sub-product trailer
        file, do we need to concatenate input trailer files
        or append to an input file.
    */
    if (strcmp(inlist, output) != 0) {
        /* Generate temporary output file name.
           This is to avoid infinite loop when output is
           accidentally the same as one of the inputs.
           Not using tmpfile() because it opens binary stream.
           Not using tmpnam() because compiler complains about danger.
        */
        strcpy(uniq_outname, uniq_outtemplate);
        if ( (td = mkstemp(uniq_outname)) < 0 ) {
            trlerror("Failed to create temporary file name.");
            freeOnExit(&ptrReg);
            return(status=INVALID_TEMP_FILE);
        }
        if ( unlink(uniq_outname) < 0) {
            trlerror("Failed to unlink temporary file name.");
            freeOnExit(&ptrReg);
            return(status=INVALID_TEMP_FILE);
        }
        if ( (tp = fdopen(td, "w+") ) == NULL ) {
            trlopenerr(uniq_outname);
            freeOnExit(&ptrReg);
            return(status=INVALID_TEMP_FILE);
        }
        addPtr(&ptrReg, tp, &fclose);
        /* Since we have determined that we are indeed combining
            input trailer files to create a new trailer file, ...
        */
        /* loop all input files */
        for (n = 0; n < c_imtlen(tpin); n++) {

            /* read the next input image name in the template list */
            c_imtgetim (tpin, trldata, CHAR_FNAME_LENGTH);

            /* open the file (read-only) and add to temp output file... */
            if ((ip = fopen (trldata, "r")) == NULL) {
                /* Report the error, but the processing can continue anyway */
                trlopenerr(trldata);
            }

            /* Do we have an input file to start with? */
            if (ip != NULL) {
                /* If so, append the input to output */
                CatTrlFile(ip, tp);
                /* Done with this input file... */
                int tmpStat = HSTCAL_OK;
                if ((tmpStat = fcloseWithStatus(&ip)))
                {
                    freeOnExit(&ptrReg); // ip isn't registered
                    return (status=tmpStat);
                }
            }

        }    /* End loop over all input files */

        fflush(tp);

        /* Copy temporary file content to output trailer file */
        CatTrlFile_NoEOF(tp, trlbuf.fp);
        freePtr(&ptrReg, tp);  /* Also delete the file because of unlink() */
        tp = NULL;
        /* Reset overwrite now to NO */
        SetTrlOverwriteMode(NO);

        /* End if (strcmp) */
    } else if (trlbuf.overwrite == YES) {

        /* We are revising/creating a single trailer file
           We want to overwrite old comments with new,
           so set the file pointer to just before the first line
           which starts with TRL_PREFIX...
        */
        if (AppendTrlFile())
        {
            /* Clean-up... */
            if (trlbuf.buffer)
            {
                free (trlbuf.buffer);
                trlbuf.buffer = NULL;
            }
            if (trlbuf.preface)
            {
                free (trlbuf.preface);
                trlbuf.preface = NULL;
            }
            status = fcloseWithStatus(&trlbuf.fp);
            freeOnExit(&ptrReg);
            return(status);
        }

        /* We are finished overwriting old comments, so reset this mode. */
        SetTrlOverwriteMode(NO);

    }

    freeOnExit(&ptrReg);
    trlbuf.init = 1;

    return (status);
}
static void CatTrlFile(FILE *ip, FILE *op)
{
    /* This function appends the ENTIRE contents of the input file (ip)
        to those of the output file (op)...
    */

    char buf[CHAR_LINE_LENGTH+1];

    buf[0] = '\0';

    /* We need to insure that we start at the beginning of
        the input file (ip)...
    */
    rewind(ip);

    /* Now copy the file into the output trailer file */
    while ( !feof(ip) ) {
        fgets(buf, CHAR_LINE_LENGTH+1, ip);
        fprintf(op,"%s",buf);
    }

    /* Flush the pipe to make sure that nothing gets left behind
        in whatever cache may be in operation.
    */
    fflush (op);
}
static void CatTrlFile_NoEOF(FILE *ip, FILE *op)
{
    // Like CatTrlFile() but ip does not have EOF

    char buf[CHAR_LINE_LENGTH+1];

    buf[0] = '\0';

    /* We need to insure that we start at the beginning of
        the input file (ip)...
    */
    rewind(ip);

    /* Now copy the file into the output trailer file */
    while ( fgets(buf, CHAR_LINE_LENGTH+1, ip) != NULL ) {
        fprintf(op,"%s",buf);
    }

    /* Flush the pipe to make sure that nothing gets left behind
        in whatever cache may be in operation.
    */
    fflush (op);
}
static int AppendTrlFile(void)
{
    /* This function sets the file pointer for the input file (ip)
        up to the line just before encountering TRL_PREFIX ...
        It reads in all the comments prior to TRL_PREFIX into memory,
        then opens the trailer file in overwrite('w+' instead of 'a+') mode
        and prints those comments to it.
        WJH 14 Apr 2000

        18 Jan 2002 (WJH) & 21 June 2002 (HAB): Removed extraneous, bug-ridden
        'tmpptr' buffer used during reallocation of oprefix buffer size.
    */

    extern int status;

    char * oprefix = malloc(initLength*sizeof(*oprefix));
    if (!oprefix){
        trlerror ("Out of memory for trailer file preface.");
        return (status = OUT_OF_MEMORY);
    }
    char buf[CHAR_LINE_LENGTH+1];
    buf[0] = '\0';
    oprefix[0] = '\0';

    /* Make sure we start searching from the beginning of the file */
    rewind (trlbuf.fp);

    while ( !feof(trlbuf.fp) )
    {
        /* Read in a line */
        fgets(buf, CHAR_LINE_LENGTH+1, trlbuf.fp);

        /* If we find the prefix, stop searching */
        if (strstr(buf,TRL_PREFIX))
            break;

        /* Store this line in a buffer to be written out when
            the old file is overwritten...
        */
        size_t newSize = (strlen(oprefix) + strlen(buf) +2 )*sizeof(*oprefix);
        assert(newSize);
        void * ptr = realloc (oprefix, newSize);
        if (!ptr)
        {
            status = OUT_OF_MEMORY;
            printfAndFlush ("Out of memory: Couldn't store trailer file comment.");
            if (oprefix)
            {
                free (oprefix);
                oprefix = NULL;
            }
            return (status);
        }
        oprefix = ptr;
        ptr = NULL;
        strcat (oprefix, buf);
    }
    /* Now we know what needs to be kept, let's close the file... */
    (void)fcloseWithStatus(&trlbuf.fp);

    /* ...reopen it with 'w+' instead of 'a+'... */
    if ( (trlbuf.fp = fopen(trlbuf.trlfile,"w+")) == NULL) {
        trlopenerr(trlbuf.trlfile);
        return(status=INVALID_TEMP_FILE);
    }
    /* ... and write out the preface information we wanted to keep. */
    fprintf (trlbuf.fp, "%s\n", oprefix);
    fflush (trlbuf.fp);

    /* Now, clean up the memory used by the temp buffer. */
    if (oprefix)
    {
        free (oprefix);
        oprefix = NULL;
    }

    return (status);
}
int WriteTrlFile (void)
{
    /*
        This function closes the trailer file.
    */

    extern int status;

    /* Now that we have copied the information to the final
        trailer file, we can close it and the temp file...
    */
    status = fcloseWithStatus(&trlbuf.fp);

    return (status);
}
int InitTrlBuf (void)
{
    /* This initialization routine must be called before any others in this
        file.
    */
    extern int status;

    trlbuf.trlfile[0] = '\0';
    trlbuf.fp = NULL;
    trlbuf.overwrite = NO;  /* Initialize with default of append */
    trlbuf.quiet = NO;      /* Initialize to produce STDOUT messages */
    trlbuf.usepref = YES;      /* Switch to specify whether to output preface */
    trlbuf.buffer = NULL;
    trlbuf.preface = NULL;

    if (!(trlbuf.buffer = malloc(initLength*sizeof(*trlbuf.buffer)))){
        trlerror ("Out of memory for trailer file buffer.");
        return (status = OUT_OF_MEMORY);
    }
    if (!(trlbuf.preface = malloc(initLength*sizeof(*trlbuf.preface)))){
        trlerror ("Out of memory for trailer file preface.");
        if (trlbuf.buffer)
        {
            free(trlbuf.buffer);
            trlbuf.buffer = NULL;
        }
        return (status = OUT_OF_MEMORY);
    }

    trlbuf.buffer[0] = '\0';
    trlbuf.preface[0] = '\0';
    trlbuf.init = 1;
    return(status);
}
void SetTrlPrefaceMode (int use)
{
    /*
        This function sets the switch to specify whether the preface should
        be output to the trailer file or not.  This is used when no CRREJ is
        done on a list of images, and you don't want the preface in the trailer
        files of the sub-products again.
    */

    trlbuf.usepref = use;

}
void SetTrlOverwriteMode (int owrite)
{
    /*
        This function sets the overwrite switch.  Once set, all comments
        after the first CALXXXBEG string will be overwritten.
        This will be used strictly for the RAW file trailer files,
        where the generic conversion comments should be kept, but
        the CALXXX comments should be overwritten with the latest results.
    */

    trlbuf.overwrite = owrite;
}
void SetTrlQuietMode (int quiet)
{
    /* This function records the value of the command-line parameter QUIET
        into the structure...
    */
    trlbuf.quiet = quiet;
}
static void AddTrlBuf (const char *message)
{
    /* Add a new message line to the buffer.
        Re-allocate space for the buffer and append line to new buffer.
    */

    /* arguments:
    char *message         i: new trailer file line to add to buffer
    */
    extern int status;

    if ( ! trlbuf.init )
        assert(0); //TRLBUF NOT INIT, YOU MAY HAVE PROBLEMS

    size_t newSize = (strlen(trlbuf.buffer) + strlen(message) +2)*sizeof(*trlbuf.buffer);
    assert(newSize);
    void * ptr = realloc (trlbuf.buffer, newSize);
    if (!ptr)
    {
        status = OUT_OF_MEMORY;
        printfAndFlush ("Out of memory: Couldn't store trailer file comment.");
        // Don't free anything of trlbuf, it should not be this func's responsibility
        return;
    }
    trlbuf.buffer = ptr;
    ptr = NULL;
    strcat (trlbuf.buffer, message);
    strcat (trlbuf.buffer, "\n"); // Append a newline at the end of every output message
}
void InitTrlPreface (void)
{
    /*
        This function will copy contents of the buffer into the preface
    */
    extern int status;

    size_t newSize = (strlen(trlbuf.buffer) +2)*sizeof(*trlbuf.preface);
    assert(newSize);
    void * ptr = realloc (trlbuf.preface, newSize);
    if (!ptr)
    {
        status = OUT_OF_MEMORY;
        printfAndFlush ("Out of memory: Couldn't store trailer file preface.");
        // Don't free anything of trlbuf, it should not be this func's responsibility
        return;
    }
    trlbuf.preface = ptr;
    ptr = NULL;
    strcpy (trlbuf.preface, trlbuf.buffer);
}
void ResetTrlPreface (void)
{
    extern int status;
    size_t newSize = initLength*sizeof(*trlbuf.preface);
    assert(newSize);
    void * ptr = realloc(trlbuf.preface, newSize);
    if (!ptr)
    {
        status = OUT_OF_MEMORY;
        printfAndFlush ("Out of memory: Couldn't resize preface for trailer file.");
        // Don't free anything of trlbuf, it should not be this func's responsibility
        return;
    }
    trlbuf.preface = ptr;
    ptr = NULL;
    *trlbuf.preface = '\0';
}
static void ResetTrlBuf (void)
{
    extern int status;
    size_t newSize = initLength*sizeof(*trlbuf.buffer);
    assert(newSize);
    void * ptr = realloc(trlbuf.buffer, newSize);
    if (!ptr)
    {
        status = OUT_OF_MEMORY;
        printfAndFlush("Out of memory: Couldn't resize buffer for trailer file.");
        // Don't free anything of trlbuf, it should not be this func's responsibility
        return;
    }
    trlbuf.buffer = ptr;
    ptr = NULL;
    *trlbuf.buffer = '\0';
}
static void WriteTrlBuf (const char *message)
{
    /*
        Write out message to trailer file, if any exist.
        If no file exists, then append message to buffer.
        Once message has been written out to a file,
        reset the trailer buffer to a single newline.
    */

    /* Check to see if trailer file is open for writing... */
    if (trlbuf.fp != NULL) {

        /* Trailer file is open, so write out buffer...*/
        if (trlbuf.preface[0] != '\0' && trlbuf.usepref == YES) {
            fprintf(trlbuf.fp,"%s\n",trlbuf.preface);
            fflush (trlbuf.fp);
            trlbuf.usepref = NO;
        }
        fprintf(trlbuf.fp,"%s\n",message);
        fflush (trlbuf.fp);

        /* ...then reset buffer, so we don't write old messages out
            next time.
        */
        if (trlbuf.buffer[0] != '\0')
            ResetTrlBuf ();

        /* The user must explicitly use 'ResetTrlPreface' in
            calling function in order to reset the preface.
            This allows the preface to be written out to
            multiple trailer files...
        */
    } else {
        /* If not, then we must append message to end of buffer. */
        AddTrlBuf (message);
    }
}
void CloseTrlBuf (struct TrlBuf * buf)
{
    /* Free memory allocated for this list, after writing out any final
        messages to the last used trailer file...
        It will then close any open trailer files as well.
    */

    if (!buf)
        return;

    extern int status;
    FILE *ofp;

    /* Do we have any messages which need to be written out? */
    if (buf->buffer[0] != '\0') {
        /* We do, so open last known trailer file and
            append these messages to that file...
        */

        if ( (ofp = fopen(buf->trlfile,"a+")) == NULL) {
            trlopenerr(buf->trlfile);
            status = INVALID_TEMP_FILE;
            goto cleanup;
        }
        fprintf (ofp,"%s",buf->buffer);

        /* Now that we have copied the information to the final
            trailer file, we can close it and the temp file...
        */
        status = fcloseWithStatus(&ofp);
    }

    cleanup: ;

        if (buf->buffer)
        {
            free (buf->buffer);
            buf->buffer = NULL;
        }
        if (buf->preface)
        {
            free (buf->preface);
            buf->preface = NULL;
        }

        status = fcloseWithStatus(&buf->fp);

}
void trlmessage (const char *message) {

    if (!message || !*message)
        return;

    /* Send output to STDOUT and explicitly flush STDOUT, if desired */
    if (trlbuf.quiet == NO) {
        printfAndFlush (message);
    }

    /* Send output to (temp) trailer file */
    WriteTrlBuf (message);

}
void trlwarn (const char *message) {

    char line[CHAR_LINE_LENGTH+1];

    /* Create full warning message, like that output in ASNWARN */
    /* Use macro to add prefix to beginning of Warning message */
    sprintf(line,"%s",WARN_PREFIX);
    strcat (line,message);

    /* Send output to (temp) trailer file */
    trlmessage (line);

}
void trlerror (const char *message) {

    char line[CHAR_LINE_LENGTH+1];

    /* Create full warning message, like that output in ASNWARN */
    /* Use macro to add prefix to beginning of ERROR message */
    sprintf(line,"%s",ERR_PREFIX);
    strcat (line,message);

    /* Send output to (temp) trailer file */
    trlmessage (line);

}
void trlopenerr (const char *filename) {
    sprintf (MsgText, "Can't open file %s", filename);
    trlerror (MsgText);
}
void trlreaderr (const char *filename) {
    sprintf (MsgText, "Can't read file %s", filename);
    trlerror (MsgText);
}
void trlkwerr (const char *keyword, const char *filename) {
    sprintf (MsgText, "Keyword \"%s\" not found in %s", keyword, filename);
    trlerror (MsgText);
}
void trlfilerr (const char *filename) {
    sprintf (MsgText, "while trying to read file %s", filename);
    trlerror (MsgText);
}
void printfAndFlush (const char *message) {
    printf ("%s\n", message);
        fflush(stdout);
}

void trlGitInfo(void)
{
    char * gitInfo = NULL;
    sprintfGitInfo(&gitInfo);
    trlmessage(gitInfo);
    if (gitInfo)
    {
        free(gitInfo);
        gitInfo = NULL;
    }
}
