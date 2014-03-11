/*
  Coded by Thibaut Jombart (tjombart@imperial.ac.uk), December 2010.
  Distributed with the adephylo package for the R software.
  Licence: GPL >=2.

  Functions based on snpbin and genlightC classes, which mirror the R classes SNPbin and genlight on the C side.
*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <R_ext/Print.h>

#include "snpbin.h"


/* Function to compute all dot products between individuals */
/* centring and scaling is always used */
/* but need to pass vectors of 0 and 1*/
void GLdotProd(unsigned char *gen, int *nbvecperind, int *byteveclength, int *nbnaperind, int *naposi, int *nind, int *nloc, int *ploidy, double *mean, double *sd, bool *freq, double *res){
	struct genlightC dat;
	int i, j, k=0;

	/* Check variance vector: do not divide by 0 */
	for(i=0;i< *nloc;i++){
		if(sd[i] < NEARZERO){
			sd[i] = 1;
		}
	}

	dat = genlightTogenlightC(gen, nbvecperind, byteveclength, nbnaperind, naposi, nind, nloc, ploidy);

	if(*freq){
		/* === working on frequencies === */
		/* Lower triangle - without the diagonal */
		for(i=0; i< (*nind-1); i++){
			for(j=i+1; j< *nind; j++){
				/* printf("\n == pair %i-%i ==\n", i+1,j+1); */
				res[k] = snpbin_dotprod_freq(&dat.x[i], &dat.x[j], mean, sd);
				++k;
			}
		}

		/* add the diagonal to the end of the array */
		for(i=0; i< *nind; i++){
			/* printf("\n == pair %i-%i == \n", i+1,i+1); */
			res[k] = snpbin_dotprod_freq(&dat.x[i], &dat.x[i], mean, sd);
			++k;
		}
	} else {
		/* === working on frequencies === */
		/* Lower triangle - without the diagonal */
		for(i=0; i< (*nind-1); i++){
			for(j=i+1; j< *nind; j++){
				/* printf("\n == pair %i-%i ==\n", i+1,j+1); */
				res[k] = snpbin_dotprod_int(&dat.x[i], &dat.x[j], mean, sd);
				++k;
			}
		}

		/* add the diagonal to the end of the array */
		for(i=0; i< *nind; i++){
			/* printf("\n == pair %i-%i == \n", i+1,i+1); */
			res[k] = snpbin_dotprod_int(&dat.x[i], &dat.x[i], mean, sd);
			++k;
		}
	}
}









void GLsumInt(unsigned char *gen, int *nbvecperind, int *byteveclength, int *nbnaperind, int *naposi, int *nind, int *nloc, int *ploidy, int *res){
	struct genlightC dat;
	int i, j;
	int *vecIntTemp;
	vecIntTemp = (int *) calloc(*nloc, sizeof(int));

	/* set res to zeros */
	/* for(j=0;j< *nloc;j++){ */
	/* 	res[j] = 0; */
	/* } */

	/* Internal C representation of the genlight object */
	dat = genlightTogenlightC(gen, nbvecperind, byteveclength, nbnaperind, naposi, nind, nloc, ploidy);


	/* === working on frequencies === */
	/* Lower triangle - without the diagonal */
	for(i=0; i < (*nind); i++){ /* for all individuals*/
		/* conversion to integers of current indiv */
		snpbin2intvec(&(dat.x[i]), vecIntTemp);

		for(j=0; j < *nloc; j++){ /* for all loci */
			if(!snpbin_isna(&(dat.x[i]), j)) res[j] += vecIntTemp[j];
		}
	}
}





void GLsumFreq(unsigned char *gen, int *nbvecperind, int *byteveclength, int *nbnaperind, int *naposi, int *nind, int *nloc, int *ploidy, double *res){
	struct genlightC dat;
	int i, j;
	double *vecFreqTemp;
	vecFreqTemp = (double *) calloc(*nloc, sizeof(double));

	/* set res to zeros */
	/* for(j=0;j< *nloc;j++){ */
	/* 	res[j] = 0.0; */
	/* } */

	/* Internal C representation of the genlight object */
	dat = genlightTogenlightC(gen, nbvecperind, byteveclength, nbnaperind, naposi, nind, nloc, ploidy);

	/* === working on frequencies === */
	/* Lower triangle - without the diagonal */
	for(i=0; i < (*nind); i++){ /* for all individuals*/
		/* conversion to frequencies of current indiv */
		snpbin2freq(&(dat.x[i]), vecFreqTemp);

		for(j=0; j < *nloc; j++){ /* for all loci */
			if(!snpbin_isna(&(dat.x[i]), j)) res[j] += vecFreqTemp[j];
		}
	}
}




/* TESTING in R */

/*

## === DOT PRODUCTS ALLELE COUNTS === ##

library(adegenet)
dat <- rbind("a"=c(1,0,0), "b"=c(1,2,1), "c"=c(1,0,1))
x <- new("genlight",dat)


## RANDOM DATA
dat <- matrix(sample(0:1, 5*1000, replace=TRUE), nrow=5)
x <- new("genlight",dat)
res1 <- glDotProd(x, alle=TRUE)
res2 <- as.matrix(x) %*% t(as.matrix(x))
all(res1==res2)


## CENTRED, NOT SCALED
res1 <- glDotProd(x, cent=TRUE, alle=TRUE)
temp <- as.matrix(x) / ploidy(x)
temp <- scalewt(temp, cent=TRUE, scale=FALSE)
res2 <- temp %*% t(temp)
res2
all(abs(res1-res2)<1e-10)


## CENTRED, SCALED
res1 <- glDotProd(x, cent=TRUE, scale=TRUE, alle=TRUE)
temp <- as.matrix(x) / ploidy(x)
temp <- scalewt(temp, cent=TRUE, scale=TRUE)
res2 <- temp %*% t(temp)
res2
all(abs(res1-res2)<1e-10)


## TEST WITH NAs
library(adegenet)
dat <- list(a=c(1,NA,0,0,2), b=c(1,2,3,4,0), c=c(NA,0,1,NA,2))
x <- new("genlight", dat) # conversion
x
res1 <- glDotProd(x, alle=TRUE)
t(data.frame(dat))
res1







## === DOT PRODUCTS ALLELE FREQUENCIES === ##

library(adegenet)


## RANDOM DATA
dat <- rbind(matrix(sample(0:1, 3*1000, replace=TRUE), nrow=3), 
             matrix(sample(0:2, 2*1000, replace=TRUE), nrow=2))
x <- new("genlight",dat)
res1 <- glDotProd(x)
temp <- as.matrix(x) / ploidy(x)
res2 <- temp %*% t(temp)
all(res1==res2)


## CENTRED, NOT SCALED
res1 <- glDotProd(x, cent=TRUE, alle=FALSE)
temp <- scalewt(temp, cent=TRUE, scale=FALSE)
res2 <- temp %*% t(temp)
res2
all(abs(res1-res2)<1e-10)


## CENTRED, SCALED
res1 <- glDotProd(x, cent=TRUE, scale=TRUE, alle=FALSE)
temp <- as.matrix(x) / ploidy(x)
temp <- scalewt(temp, cent=TRUE, scale=TRUE)
res2 <- temp %*% t(temp)
res2
all(abs(res1-res2)<1e-10)


## TEST WITH NAs
library(adegenet)
dat <- list(a=c(1,NA,0,0,2), b=c(1,2,3,4,0), c=c(NA,0,1,NA,2))
x <- new("genlight", dat) # conversion
x
res1 <- glDotProd(x, alle=FALSE)
temp <- as.matrix(x)/ploidy(x)
temp
res1

*/

