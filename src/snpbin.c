/*
  Coded by Thibaut Jombart (tjombart@imperial.ac.uk), December 2010.
  Distributed with the adephylo package for the R software.
  Licence: GPL >=2.

  These functions are designed to recode genotypes given as binary integers into new integers which
  map them to unique bytes. One genotype of 8 binary SNPs is mapped uniquely (bijectively) to a
  value between 0 and 255. This is achieved by considering the genotype 'x' in the basis 2^0
  ... 2^7, and summing the values of the vector in this basis. That is, we use the function:

  {0,1}^8 |-> {0,...,255}
  x -> x_1 * 2^0 + ... + x_8 * 2^7 = \sum_i x_i * 2^(i-1)


  # Function named as 'SNPbin...' or 'GL...' are to be called directly from R.
  # The structure 'snpbin' is a C representation of the class 'SNPbin'.
  # Function named as 'snpbin...' are made to be called internally.
*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <R_ext/Print.h>

#include "adesub.h"
#include "snpbin.h"

/* #define NEARZERO 0.0000000001 */
/* #define TRUE 1 */
/* #define FALSE 0 */

/* typedef short bool; */

/* /\* */
/*    ========================= */
/*    === CLASS DEFINITIONS === */
/*    ========================= */
/* *\/ */

/* /\* 'bytevecnb' arrays of bytes concatenated into a single array *\/ */
/* /\* of dim 'byteveclength' x 'bytevecnb' *\/ */
/* /\* nloc is the number of SNPs - used for recoding to integers *\/ */
/* /\* naposi indicates the positions of NAs *\/ */
/* /\* nanb is the length of naposi *\/ */

/* struct snpbin{ */
/* 	unsigned char *bytevec; */
/* 	int *byteveclength, *bytevecnb, *nloc, *nanb, *naposi, *ploidy; /\* all but naposi have length 1 *\/ */
/* }; */




/* struct genlightC{ */
/* 	struct snpbin *x; */
/* 	int *nind; */
/* }; */





/*
  ===========================
  === AUXILIARY FUNCTIONS ===
  ===========================
*/


struct snpbin makesnpbin(unsigned char *bytevec, int *byteveclength, int *bytevecnb, int *nloc, int *nanb, int *naposi, int *ploidy) {
    struct snpbin out;
    int i;

    out.bytevec = bytevec;
    out.byteveclength = byteveclength;
    out.bytevecnb = bytevecnb;
    out.nloc = nloc;
    out.nanb = nanb;
    /* need to decrease the indices of NAs by 1, e.g. [1-10]->[0-9] */
    out.naposi = naposi;
    if(*nanb > 0){
	for(i=0;i< *nanb; i++){
	    out.naposi[i] = out.naposi[i] - 1;
	}
    }
    out.ploidy = ploidy;
    return out;
}


/* Maps one byte from 0-255 to sequences of 8 (binary) integers values */
void byteToBinInt(unsigned char in, int *out){
    short int rest, i, temp;

    rest = (int)in;

    /* initialize all values to 0*/
    for(i=0;i<=7;i++)
	out[i]=0;

    for(i=7;i>=0;i--){
	temp = pow(2, i);
	if(rest >= temp) {
	    out[i] = 1;
	    rest = rest- temp;
	    if(rest == 0) break;
	}
    }
}




/* Maps one byte from 0-255 to sequences of 8 (binary) double values */
void byteToBinDouble(unsigned char in, double *out){
    short rest, i, temp;

    rest = (int) in;

    /* initialize all values to 0*/
    for(i=0;i<=7;i++)
	out[i]=0.0;

    for(i=7;i>=0;i--){
	temp = pow(2, i);
	if(rest >= temp) {
	    out[i] = 1.0;
	    rest = rest- temp;
	    if(rest == 0) break;
	}
    }
}





/* Maps an array of values from 0-255 to sequences of 8 binary values */
/* Input are unsigned char (hexadecimal), outputs are integers */
void bytesToBinInt(unsigned char *vecbytes, int *vecsize, int *vecres){
    int i, j, idres=0, *temp; /* idres: index in vecres*/

    temp = (int *) calloc(8, sizeof(int));

    for(i=0;i<*vecsize;i++){
	byteToBinInt(vecbytes[i], temp);
	for(j=0;j<=7;j++){
	    vecres[j+idres] = temp[j];
	}
	idres = idres + 8;
    }

    free(temp);
} /* end binIntToBytes*/










/*
  ===============================
  === MAIN EXTERNAL FUNCTIONS ===
  ===============================
*/



/* Maps an array of values from 0-255 to integers representing counts of alleles */
/* This is done by adding arrays of 0-1 for indiv with ploidy > 1*/
/* Input are unsigned char (hexadecimal), outputs are integers */
/* veclength is the length of one vector of bytes */
/* nbvec is the nb of input vectors*/
/* input 'vecbytes' is actually concatenated, ie of size veclength * nbvec */
void bytesToInt(unsigned char *vecbytes, int *veclength, int *nbvec, int *vecres, int *reslength){
    int i, j, k, idres=0, *temp; /* idres: index in vecres*/

    temp = (int *) calloc(8, sizeof(int));

    /* initialize result vector to 0 */
    for(i=0; i < *reslength; i++){
	vecres[i]=0;
    }

    /* build output */
    for(k=0;k<*nbvec;k++){ /* for all input vector */
	idres = 0;
	for(i=0;i<*veclength;i++){ /* for one input vector */
	    byteToBinInt(vecbytes[i+ k* *veclength], temp); /* byte -> 8 int (0/1)*/
	    for(j=0;j<=7;j++){ /* fill in the result*/
		vecres[j+idres] += temp[j];
	    }
	    idres = idres + 8;
	}
    }
    free(temp);
} /* end bytesToInt */





void bytesToDouble(unsigned char *vecbytes, int *veclength, int *nbvec, double *vecres, int *reslength){
    int i, j, k, idres=0; /* idres: index in vecres*/
    double *temp;
    temp = (double *) calloc(8, sizeof(double));

    /* initialize result vector to 0 */
    for(i=0; i < *reslength; i++){
	vecres[i]=0.0;
    }

    for(k=0;k<*nbvec;k++){ /* for all input vector */
	idres = 0;
	for(i=0;i<*veclength;i++){ /* for one input vector */
	    byteToBinDouble(vecbytes[i+ k* *veclength], temp); /* byte -> 8 double (0/1)*/
	    for(j=0;j<=7;j++){ /* fill in the result*/
		vecres[j+idres] += temp[j];
	    }
	    idres = idres + 8;
	}
    }
    free(temp);
} /* end bytesToInt */





/* 
   === MAP BINARY SNPS TO 1->256 SCALE ===
   - vecsnp: vector of integers (0/1)
   - vesize: length of vecsnp
   - res: vector of integers valued on 0:255
   - ressize: length of res
*/
void binIntToBytes(int *vecsnp, int *vecsize, unsigned char *vecres, int *ressize){
    /* declarations */
    int i, j, idres, *binBasis; /* must use dynamic allocation */

    /* allocate memory for local variables */
    vecintalloc(&binBasis, 8);

    /* define binary basis */
    for(i=1; i<=8; i++){
	binBasis[i] = pow(2, i-1);
    }

    /* set all values of vecres to 0 */
    for(i=0;i < *ressize;i++){
	vecres[i] = 0x00;
    }



    /* INDICES */
    /* i: idx of snp */
    /* j: idx of binBasis (1:8) */
    /* idres: idx in vector of results */

    idres = 0;
    j = 1;
    for(i=0;i< *vecsize;i++){
	vecres[idres] = vecres[idres] + (unsigned char)(binBasis[j] * vecsnp[i]);
	if(j == 8){
	    idres++;
	    j = 1;
	} else {
	    j++;
	}
    }


    /* free memory */
    freeintvec(binBasis);

} /* end binIntToBytes */








/*
  =====================
  === CLASS METHODS ===
  =====================
*/

int nLoc(struct snpbin *x){
    return *(x->nloc);
}


int ploidy(struct snpbin *x){
    return *(x->ploidy);
}





/* transform a snpbin into a vector of integers */
void snpbin2intvec(struct snpbin *x, int *out){
    int *temp;
    temp= (int *) calloc(1, sizeof(int));
    *temp=nLoc(x);
    bytesToInt(x->bytevec, x->byteveclength, x->bytevecnb, out, temp);
    free(temp);
/*reminders:
  - void bytesToInt(unsigned char *vecbytes, int *veclength, int *nbvec, int *vecres, int reslength){
  - snpbin: unsigned char *bytevec; int *byteveclength, *bytevecnb, *nloc, *nanb, *naposi; */
}



/* transform a snpbin into a vector of frequencies (double) */
void snpbin2freq(struct snpbin *x, double *out){
    double ploid = (double) ploidy(x);
    int *temp;
    temp= (int *) calloc(1, sizeof(int));
    *temp=nLoc(x);
    bytesToDouble(x->bytevec, x->byteveclength, x->bytevecnb, out, temp);
    int i;

    for(i=0; i < nLoc(x); i++){
	out[i] = out[i] / ploid;
    }
    free(temp);
/*reminders:
  - void bytesToInt(unsigned char *vecbytes, int *veclength, int *nbvec, int *vecres, int reslength){
  - snpbin: unsigned char *bytevec; int *byteveclength, *bytevecnb, *nloc, *nanb, *naposi; */
}



/* print a snpbin object - used for debugging */
void printsnpbin(struct snpbin *x){
    int i, *temp;
    temp = (int *) calloc(nLoc(x), sizeof(int));
    snpbin2intvec(x, temp);


    for(i=0;i< *(x->byteveclength);i++){
	Rprintf("%i ", (int) (x->bytevec)[i]);
	/* printf("%i ", (int) (x->bytevec)[i]); */
    }
    Rprintf("   ");
    for(i=0;i<nLoc(x);i++){
	Rprintf("%i ", temp[i]);
	/* printf("%i ", temp[i]); */
    }

    Rprintf("NA posi: ");
    for(i=0;i< *(x->nanb);i++){
	Rprintf("%i ", (x->naposi)[i]);
	/* printf("%i ", (x->naposi)[i]); */
    }

    free(temp);
}



short int snpbin_isna(struct snpbin *x, int i){
    int j = 0;
    if(*(x->nanb) < 1 || i > nLoc(x)) return 0;

    while(j < *(x->nanb)){
	if( i == (x->naposi)[j]) return 1;
	j++;
    }

    return 0;
}





/* Function to compute one dot products between two individuals */
/* centring and scaling is always used */
/* but need to pass vectors of 0 and 1*/
double snpbin_dotprod_int(struct snpbin *x, struct snpbin *y, double *mean, double *sd){
    /* define variables, allocate memory */
    int P = nLoc(x), i;
    double res = 0.0;
    int *vecx, *vecy;
    vecx = (int *) calloc(P, sizeof(int));
    vecy = (int *) calloc(P, sizeof(int));

    /* conversion to integers */
    snpbin2intvec(x, (int *) vecx);
    snpbin2intvec(y, (int *) vecy);


    /* printf("\nvector x: \n"); */
    /* for(i=0;i<P;i++){ */
    /* 	printf("%i", vecx[i]); */
    /* } */

    /* printf("\nvector y: \n"); */
    /* for(i=0;i<P;i++){ */
    /* 	printf("%i", vecy[i]); */
    /* } */

    /* compute dot product */
    for(i=0;i<P;i++){
	if(snpbin_isna(x,i) == 0 && snpbin_isna(y,i) == 0){
	    /* res += ((vecx[i]-mean[i])/sd[i]) * ((vecy[i]-mean[i])/sd[i]); */
	    res += ((vecx[i]-mean[i])/sd[i]) * ((vecy[i]-mean[i])/sd[i]);
	    /* printf("\ntemp value of increment: %f", ((vecx[i]-mean[i])/sd[i]) * ((vecy[i]-mean[i])/sd[i])); */
	    /* printf("\ntemp value of result: %f", res); */
	}
    }

    /* free memory */
    free(vecx);
    free(vecy);

    return res;
}






double snpbin_dotprod_freq(struct snpbin *x, struct snpbin *y, double *mean, double *sd){
    /* define variables, allocate memory */
    int P = nLoc(x), i;
    double res = 0.0;
    double *vecx, *vecy;
    vecx = (double *) calloc(P, sizeof(double));
    vecy = (double *) calloc(P, sizeof(double));
	
    /* conversion to integers or frequencies*/
    snpbin2freq(x, vecx);
    snpbin2freq(y, vecy);
	

    /* printf("\nvector x: \n"); */
    /* for(i=0;i<P;i++){ */
    /* 	printf("%i", vecx[i]); */
    /* } */

    /* printf("\nvector y: \n"); */
    /* for(i=0;i<P;i++){ */
    /* 	printf("%i", vecy[i]); */
    /* } */

    /* compute dot product */
    for(i=0;i<P;i++){
	if(snpbin_isna(x,i) == 0 && snpbin_isna(y,i) == 0){
	    /* res += ((vecx[i]-mean[i])/sd[i]) * ((vecy[i]-mean[i])/sd[i]); */
	    res += ((vecx[i]-mean[i])/sd[i]) * ((vecy[i]-mean[i])/sd[i]);
	    /* printf("\ntemp value of increment: %f", ((vecx[i]-mean[i])/sd[i]) * ((vecy[i]-mean[i])/sd[i])); */
	    /* printf("\ntemp value of result: %f", res); */
	}
    }

    /* free memory */
    free(vecx);
    free(vecy);

    return res;
}






/* Function to convert a 'genlight' object (R side) into an array of 'snpbin' (C side) */
/* Each component of the genlight is concatenated into a single vector */
/* and then used to create different 'snpbin' on the C side */
struct genlightC genlightTogenlightC(unsigned char *gen, int *nbvecperind, int *byteveclength, int *nbnaperind, int *naposi, int *nind, int *nloc, int *ploidy){
    /* declare variables and allocate memory */
    int i, idxByteVec=0, idxNAVec=0;
    struct genlightC out;
    out.x = (struct snpbin *) calloc(*nind, sizeof(struct snpbin));

    /* create the list of snpbin */
    /* printf("\n nind: %d\n", *nind); */
    for(i=0; i < *nind; i++){
	out.x[i] = makesnpbin(&gen[idxByteVec], byteveclength, &nbvecperind[i], nloc, &nbnaperind[i], &naposi[idxNAVec], &ploidy[i]);
	idxByteVec += *byteveclength * nbvecperind[i]; /* update index in byte array */
	idxNAVec +=  nbnaperind[i]; /* update index in byte array */
	/* printf("\nimported genotype %i: ", i+1); */
	/* printsnpbin(&out.x[i]); */
    }

    /* printf("step 3"); */

    out.nind = nind;

    /* printf("step 4"); */
    return out;
}







/*
  =========================
  === TESTING FUNCTIONS ===
  =========================
*/

/* Simple test function */
/* Test: increases for a raw (unsigned char) vector */
void testRaw(unsigned char *a, int *n){
    int i;
    for(i=0; i< *n; i++){
	a[i] = (unsigned char)(i);
    }
}



/* Test: increases for a raw (unsigned char) vector */
void testSizePointer(int *sizePointer, int *sizeFirstElement, int *nbElements){
    double *a;
    a = (double *) calloc(5, sizeof(double));
    *sizePointer = sizeof(a);
    *sizeFirstElement = sizeof(a[0]);
    *nbElements = sizeof(a) / sizeof(a[0]);
    free(a);
}


/* TESTING in R */

/*
  ## test raw conversion
  .C("testRaw", raw(256), 256L, PACKAGE="adegenet")
  .C("testSizePointer", integer(1), integer(1), integer(1), PACKAGE="adegenet")

  ## test raw->int conversion
  x <- sample(0:1,800,replace=TRUE)
  toto <- .bin2raw(x)$snp
  all(.C("bytesToBinInt", toto, length(toto), integer(length(toto)*8))[[3]]==x)

  ## test raw vec -> binary integers
  .C("bytesToBinInt",as.raw(c(12,11)), 2L, integer(16), PACKAGE="adegenet")

  ## test several raw vec -> int (allele counts, any ploidy)
  .C("bytesToInt",as.raw(c(12,11)), 1L, 2L, integer(8), integer(16), PACKAGE="adegenet")


*/

