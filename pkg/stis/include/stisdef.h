#ifndef INCL_STISDEF_H
#define INCL_STISDEF_H

/* Prototypes for public functions that need definitions from hstio,
   and for related functions.
*/

int add2d (SingleGroup *a, SingleGroup *b);

int addk2d (SingleGroup *a, float k);

int bin2d (SingleGroup *a, int xcorner, int ycorner, int binx, int biny,
	int avg, SingleGroup *b);

int BinCoords (Hdr *inhdr, double *block, double *offset,
		Hdr *scihdr, Hdr *errhdr, Hdr *dqhdr);

int div2d (SingleGroup *a, SingleGroup *b);

int doStat (SingleGroup *out, short sdqflags);

int GetCorner (Hdr *hdr, int rsize, int *bin, int *corner);

int GetLT  (Hdr *hdr, double *ltm, double *ltv);
int GetLT0 (Hdr *hdr, double *ltm, double *ltv);

int GetMOC (RefTab *mofftab, char *opt_elem, int cenwave,
		int *mref, double *yref, double *a4corr);

int GetRefName (RefFileInfo *ref, Hdr *phdr, char *keyword, char *refname);

int GetSwitch (Hdr *phdr, char *calswitch, int *flag);

int ImgHistory (RefImage *ref, Hdr *phdr);
int TabHistory (RefTab *ref, Hdr *phdr);

int ImgPedigree (RefImage *ref);
int TabPedigree (RefTab *ref);
int RowPedigree (RefTab *ref, int row,
	IRAFPointer tp, IRAFPointer cp_pedigree, IRAFPointer cp_descrip);

void Interp2D (SingleGroup *in, short sdqflags,
	double ix, double iy, double jacobian,
	int err_algorithm,
	float *oSci, float *oErr, short *oDQ);

int Get_KeyD (Hdr *hd, char *keyword,
	int use_def, double def, double *value);
int Get_KeyF (Hdr *hd, char *keyword,
	int use_def, float def, float *value);
int Get_KeyI (Hdr *hd, char *keyword,
	int use_def, int def, int *value);
int Get_KeyB (Hdr *hd, char *keyword,
	int use_def, Bool def, Bool *value);
int Get_KeyS (Hdr *hd, char *keyword,
	int use_def, char *def, char *value, int maxch);
int Put_KeyD (Hdr *hd, char *keyword, double value, char *comment);
int Put_KeyF (Hdr *hd, char *keyword, float value, char *comment);
int Put_KeyI (Hdr *hd, char *keyword, int value, char *comment);
int Put_KeyB (Hdr *hd, char *keyword, Bool value, char *comment);
int Put_KeyS (Hdr *hd, char *keyword, char *value, char *comment);

int mult2d (SingleGroup *a, SingleGroup *b);

int multk2d (SingleGroup *a, float k);

int sub2d (SingleGroup *a, SingleGroup *b);

void UCalVer (Hdr *phdr);

void UFilename (char *filename, Hdr *phdr);

int unbin2d (SingleGroup *a, SingleGroup *b);

#endif /* INCL_STISDEF_H */
