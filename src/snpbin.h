#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


#define NEARZERO 0.0000000001
#define TRUE 1
#define FALSE 0

typedef short bool;



/*
   =========================
   === CLASS DEFINITIONS ===
   =========================
*/

/* 'bytevecnb' arrays of bytes concatenated into a single array */
/* of dim 'byteveclength' x 'bytevecnb' */
/* nloc is the number of SNPs - used for recoding to integers */
/* naposi indicates the positions of NAs */
/* nanb is the length of naposi */
struct snpbin{
	unsigned char *bytevec;
	int *byteveclength, *bytevecnb, *nloc, *nanb, *naposi, *ploidy; /* all but naposi have length 1 */
};



struct genlightC{
	struct snpbin *x;
	int *nind;
};



/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/


void byteToBinInt(unsigned char in, int *out);
void byteToBinDouble(unsigned char in, double *out);
void bytesToBinInt(unsigned char *vecbytes, int *vecsize, int *vecres);
struct snpbin makesnpbin(unsigned char *bytevec, int *byteveclength, int *bytevecnb, int *nloc, int *nanb, int *naposi, int *ploidy);





/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

void bytesToInt(unsigned char *vecbytes, int *veclength, int *nbvec, int *vecres, int *reslength);
void bytesToDouble(unsigned char *vecbytes, int *veclength, int *nbvec, double *vecres, int *reslength);
void binIntToBytes(int *vecsnp, int *vecsize, unsigned char *vecres, int *ressize);





/*
   =====================
   === CLASS METHODS ===
   =====================
*/

int nLoc(struct snpbin *x);
int ploidy(struct snpbin *x);
void snpbin2intvec(struct snpbin *x, int *out);
void snpbin2freq(struct snpbin *x, double *out);
void printsnpbin(struct snpbin *x);
short int snpbin_isna(struct snpbin *x, int i);
double snpbin_dotprod_int(struct snpbin *x, struct snpbin *y, double *mean, double *sd);
double snpbin_dotprod_freq(struct snpbin *x, struct snpbin *y, double *mean, double *sd);
struct genlightC genlightTogenlightC(unsigned char *gen, int *nbvecperind, int *byteveclength, int *nbnaperind, int *naposi, int *nind, int *nloc, int *ploidy);








/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/


void testRaw(unsigned char *a, int *n);
void testSizePointer(int *sizePointer, int *sizeFirstElement, int *nbElements);
