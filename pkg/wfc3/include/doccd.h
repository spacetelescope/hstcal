#ifndef INCL_DOCCD_H
#define INCL_DOCCD_H

#include "hstio.h"
#include "wf3info.h"

/*USEFUL LIB FUNCTIONS*/
static void AtoDMsg (WF3InfoRef *, int);
static void BiasMsg (WF3InfoRef *, int);
static void SatMsg (WF3InfoRef *, int);
static void FlashMsg (WF3InfoRef *, int);
static void BlevMsg (WF3InfoRef *, int);
static void dqiMsg  (WF3InfoRef *, int);

int checkBinned (SingleGroup *);
int GetCorner (Hdr *, int , int *, int *);
int doAtoD (WF3InfoRef *, SingleGroup *);
int atodHistory (WF3InfoRef *, Hdr *);
int doBias (WF3InfoRef *, SingleGroup *);
int doFullWellSat(WF3InfoRef *, SingleGroup *);
int biasHistory (WF3InfoRef *, Hdr *);
int doFlash (WF3InfoRef *, SingleGroup *, float *);
int flashHistory (WF3InfoRef *, Hdr *);
int doBlev (WF3InfoRef *, SingleGroup *, int, float *, int *, int *);
int blevHistory (WF3InfoRef *, Hdr *, int, int);
int CCDHistory (WF3InfoRef *, Hdr *);
int doDQI (WF3InfoRef *, SingleGroup *, int overscan);
int dqiHistory (WF3InfoRef *, Hdr *);
int doNoise (WF3InfoRef *, SingleGroup *, int *);
int noiseHistory (Hdr *);
int GetGrp (WF3InfoRef *, Hdr *);
int OmitStep (int);
int PutKeyFlt (Hdr *, char *, float, char *);
int PutKeyInt (Hdr *, char *, int, char *);
void PrSwitch (char *, int);
void PrRefInfo (char *, char *, char *, char *, char *);
void UCalVer (Hdr *);
void UFilename (char *, Hdr *);
int FindOverscan (WF3InfoRef *, int, int, int *);
int GetCCDTab (WF3InfoRef *, int, int);
int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
int SinkDetect (WF3InfoRef *, SingleGroup *);
int FindLine (SingleGroup *, SingleGroupLine *, int *, int *,int *,int *, int *);

int Full2Sub(WF3InfoRef *, SingleGroup *, SingleGroup *, int, int, int);
int Sub2Full(WF3InfoRef *, SingleGroup *, SingleGroup *, int, int, int );
int CreateEmptyChip(WF3InfoRef *, SingleGroup *);

#endif /* INCL_DOCCD_H */
