#ifndef INCL_DOCCD_H
#define INCL_DOCCD_H

/*USEFUL LIB FUNCTIONS*/
static void AtoDMsg (WF3Info *, int);
static void BiasMsg (WF3Info *, int);
static void SatMsg (WF3Info *, int);
static void FlashMsg (WF3Info *, int);
static void BlevMsg (WF3Info *, int);
static void dqiMsg  (WF3Info *, int);

int checkBinned (SingleGroup *);
int GetCorner (Hdr *, int , int *, int *);
int doAtoD (WF3Info *, SingleGroup *);
int atodHistory (WF3Info *, Hdr *);
int doBias (WF3Info *, SingleGroup *);
int doFullWellSat(WF3Info *, SingleGroup *);
int biasHistory (WF3Info *, Hdr *);
int doFlash (WF3Info *, SingleGroup *, float *);
int flashHistory (WF3Info *, Hdr *);
int doBlev (WF3Info *, SingleGroup *, int, float *, int *, int *);
int blevHistory (WF3Info *, Hdr *, int, int);
int CCDHistory (WF3Info *, Hdr *);
int doDQI (WF3Info *, SingleGroup *, int overscan);
int dqiHistory (WF3Info *, Hdr *);
int doNoise (WF3Info *, SingleGroup *, int *);
int noiseHistory (Hdr *);
int GetGrp (WF3Info *, Hdr *);
int OmitStep (int);
int PutKeyFlt (Hdr *, char *, float, char *);
int PutKeyInt (Hdr *, char *, int, char *);
void PrSwitch (char *, int);
void PrRefInfo (char *, char *, char *, char *, char *);
void TimeStamp (char *, char *);
void UCalVer (Hdr *);
void UFilename (char *, Hdr *);
int FindOverscan (WF3Info *, int, int, int *);
int GetCCDTab (WF3Info *, int, int);
int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
int SinkDetect (WF3Info *, SingleGroup *);
int FindLine (SingleGroup *, SingleGroupLine *, int *, int *,int *,int *, int *);

int Full2Sub(WF3Info *, SingleGroup *, SingleGroup *, int, int, int);
int Sub2Full(WF3Info *, SingleGroup *, SingleGroup *, int, int, int );
int CreateEmptyChip(WF3Info *, SingleGroup *);

#endif /* INCL_DOCCD_H */
