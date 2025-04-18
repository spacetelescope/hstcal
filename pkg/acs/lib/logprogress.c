# include <stdio.h>
# include <time.h>
#include "hstcal.h"
# include "acs.h"	/* for message output */

void LogProgress (char *mess, int n) {

	char t[82];
	time_t *tp, now;

	tp = NULL;
	now = time (tp);
	strftime (t, 82, "%a %H:%M:%S %Z %d-%b-%Y", localtime (&now));

	if (n > 0)
	    snprintf(MsgText, sizeof(MsgText), "%s %s %d.", t, mess, n);
	else
	    snprintf(MsgText, sizeof(MsgText), "%s %s", t, mess);

	trlmessage (MsgText);
}
