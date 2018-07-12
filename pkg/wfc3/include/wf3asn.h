#ifndef INCL_WF3ASN_H
#define INCL_WF3ASN_H

# define	BADVAL	-99999

/* Image types */
enum imageTypes_ {EXP, CRJ, RPT, DTH, CRC, CDTH};
typedef enum imageTypes_ imageTypes;

/* Observation type */
enum obsTypes_ {EXT, DARK, IFLAT, EFLAT}; /* Expand as needed */
typedef enum obsTypes_ obsTypes; 

/* Process type: Single image, single product or full table */
enum asnType_ {SINGLE, PARTIAL, FULL};
typedef enum asnType_ asnType;

/* Observation and calibration information */
typedef struct {
	char imagetype[SZ_FITS_VAL+1];
	obsTypes obs_type;
} ObsInfo;

typedef struct {
	float  crpix[2];
	double crval[2];
	double cd[2][2];
	char   ctype[2][SZ_FITS_VAL+1];
} WCS; /* Not needed for CALWF3*/

typedef struct {
	char name[CHAR_FNAME_LENGTH+1];
	char mtype[SZ_FITS_VAL+1];
	Bool status;
	/* WCS wcs;  */
	imageTypes type;
	int posid;
	int prodid;
	int dithid;
	float dx, dy;
	float xi, yi;
} MemberInfo;


typedef struct {
	char name[CHAR_FNAME_LENGTH+1];
	char mtype[SZ_FITS_VAL+1];
	Bool prsnt;
	int asnrow;		      /* row from ASN table */
	char expname[CHAR_FNAME_LENGTH+1];    /* Full filenames for member EXP images */
    char blv_tmp[CHAR_FNAME_LENGTH+1];    /* BLV_TMP files to be deleted */
    char blc_tmp[CHAR_FNAME_LENGTH+1];  /*BLC_TMP file to be deleted*/
    char rac_tmp[CHAR_FNAME_LENGTH+1];  /*RAC file to be deleted*/
	float dx, dy;
	float xi, yi;
} ExpInfo;

typedef struct {
	char name[CHAR_FNAME_LENGTH+1];
	char mtype[SZ_FITS_VAL+1];
	Bool prsnt;
	int asnrow;			/* row from ASN table */
	int numexp;			/* # of EXP making this sub-product */
	int posid;			/* Numbering starts at 1 */
	char spname[CHAR_FNAME_LENGTH+1];     	/* Full filename for sub-product */
	char spname_cte[CHAR_FNAME_LENGTH+1];     	/* Full filename for sub-product */
	char crj_tmp[CHAR_FNAME_LENGTH+1];	/* CRJ_TMP file(s) to be deleted */
    char crc_tmp[CHAR_FNAME_LENGTH+1];  /*CRC_TMP file to be deleted*/
	ExpInfo *exp;			/* List of member EXP information */
} SubProdInfo;

typedef struct {
	char name[CHAR_FNAME_LENGTH+1];
	char mtype[SZ_FITS_VAL+1];
	Bool prsnt;
	int asnrow;		    /* row from ASN table */
	int numsp;		    /* # of sub-products making this product */
	int prodid;		    /* Numbering starts at 0 */	
	char prodname[CHAR_FNAME_LENGTH+1]; /* Full filename for final product */
    char prodname_cte[CHAR_FNAME_LENGTH+1]; /*Full filename for final product with CTE corrected data*/
	SubProdInfo *subprod;	    /* List of member sub-product information */
} ProdInfo;

typedef struct {
	char input[CHAR_FNAME_LENGTH+1];	/* initial input  */
	char filename[CHAR_FNAME_LENGTH+1]; 	/* full filename of input */
	char rootname[CHAR_FNAME_LENGTH+1];  	/* Rootname derived from input */	
	char asn_table[CHAR_FNAME_LENGTH+1];
	asnType process;      /* single image, partial or full ASN processing */
	int crcorr;		/* do cosmic-ray rejection for science files? */
	int rptcorr;		/* combine repeatobs science data? */
	int dthcorr;  	/* Dither combine sub-products into final product? */
	int numprod;	/* Number of dither-combined/polarizer products (>1) */
	int numsp;		/* Number of sub-products in association */
	int *spmems;		/* Number of EXP for each sub-product */
	int numasn;		/* number of rows in association table	*/
	ProdInfo *product;	/* Set to NULL for SINGLE image processing */
	char instr[SZ_FITS_VAL+1];
	int  detector;
	int verbose;		/* Print out messages during processing? */
	int debug;			/* Print out DEBUGGING messages */
} AsnInfo;


#endif /* INCL_WF3ASN_H */
