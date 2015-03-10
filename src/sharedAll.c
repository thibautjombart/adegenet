/* 
*** THIS CODE IS PART OF THE *adegenet* PACKAGE FOR R. ***

Code used to compute the proportion of alleles shared among a set of individuals.
The arguments are: 
- matAll a matrix containing alleles coded by integers, with genotypes in rows and loci in columns. Results from the binding by columns of A1 and A2, where A1 stores the one allele and A2 the other allele.
- nRow: the number of rows of matAll, i.e. number of genotypes
- nCol: the number of columns of matAll, i.e. twice the number of loci
- resVec: a vector of length (n(n-1)/2) storing the proportion of shared alleles for each couple of genotypes.

Thibaut Jombart (t.jombart@imperial.ac.uk), 2008.
*/ 

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R_ext/Print.h>

#include "adesub.h"


/* THIS FUNCTION IS DEPRECATED */
void sharedAll(int *matAll, int *nRow, int *nCol, double *resVec)
{
  /* Declare local C variables */
  int i, i1, i2, j, k, n, p, nbAll, **mat;
  n = *nRow;
  p = *nCol;

  int nLoc=p/2;

  /* Memory allocation for local C variables */

  tabintalloc(&mat, n, p); /* function from ade4 */

  /* Local reconstruction of the matrix of alleles
     ! beware: this matrix is to be used from 1 to n and 1 to p included,
     and not from zero to n and p excluded as it is common in C */
  k = 0;
  for (j=1; j<=p; j++) {
    for (i=1; i<=n; i++) {
      mat[i][j] = matAll[k];
      k = k + 1;
    }
  }

  /* == Main Computations: ==
     - i1, i2: indices of genotypes
     - j: index of allele
     - n: number of genotypes
     - p number of columns in mat (i.e. twice the number of loci)
     - each term in mat codes an allele (NAs are coded by 0)
  */

  k=0; /* counter used to browse resVec */
  for(i1=1; i1<=(n-1); i1++){
    for(i2=(i1+1); i2<=n; i2++){
      /* Used for debugging
	 printf("\n\n debug: ## %d-%d ##",i1,i2);
      */

      resVec[k] = 0.0; /* Initialization of the result */
      nbAll = 0; /* counts the number of types alleles */
      for(j=1; j<=nLoc; j++){
	/* Used for debugging
	   printf("\n debug: j=%d",j);
	   printf("\n debug: mat[i1,j]=%d",mat[i1][j]);
	   printf("\n debug: mat[i1,j]=%d",mat[i1][j+nLoc]);
	   printf("\n debug: mat[i2,j]=%d",mat[i2][j]);
	   printf("\n debug: mat[i2,j+nLoc]=%d",mat[i2][j+nLoc]);
	*/
	if(mat[i1][j] != 0 && mat[i1][j+nLoc] != 0 && 
	   mat[i2][j] != 0 && mat[i2][j+nLoc] != 0){
	  /* Used for debugging 
	     printf("\n debug: alleles are typed");
	  */
	  nbAll+=2;
	  /* Used for debugging
	     printf("\n debug: nbAll=%d", nbAll);
	  */
	  /* Compare alleles: 
	     -> either both alleles are in common, 
	     -> or no allele are common, 
	     -> or there is one common allele */
	  /* both alleles common */
	  if((mat[i1][j] == mat[i2][j] 
	      && mat[i1][j+nLoc] == mat[i2][j+nLoc])
	     || (mat[i1][j] == mat[i2][j+nLoc]
		 && mat[i1][j+nLoc] == mat[i2][j])){
	    resVec[k] += 2.0;
	  } else if(!( /* if not 'all alleles differe' */
		      mat[i1][j] != mat[i2][j] 
		      && mat[i1][j] != mat[i2][j+nLoc]
		      && mat[i1][j+nLoc] != mat[i2][j] 
		      && mat[i1][j+nLoc] != mat[i2][j+nLoc])
		    ) resVec[k]++;

	} /* end if */
      } /* end for j in 1:(nLoc) */

      /* Divide the number of shared alleles by the number of typed alleles */
      if(nbAll > 0) resVec[k] = resVec[k]/nbAll;
      /*printf("\n debug: resVec[i1,i2]/nbAll (%d,%d)=# %f #", i1,i2,resVec[k]);*/
      k++;

    } /* end for i2 */
  } /* end for i1*/

  /* Free allocated memory */
  freeinttab(mat);

} /* end sharedAll */



/* SMALL FUNCTION TO RETURN THE SMALLEST OF 2 INTEGERS */
int min_int(int a, int b){
  if(a<b) return a;
  return b;
}




/* THIS IS THE FUNCTION TO USE */
void nb_shared_all(int *in, int *out, int *nind, int *ncol){
  int i, j, k, counter=0, **mat, n = *nind, p = *ncol;

  /* allocate memory for table of allele nb */
  tabintalloc(&mat, n, p);


  /* reconstruct table of allele nb */
  for(j=1;j<=p;j++){
    for(i=1;i<=n;i++){
      mat[i][j] = in[counter++];
    }
  }


  /* perform computations */
  counter = 0;
  for(i=1;i<=(n-1);i++){
    for(j=i+1;j<=n;j++){
      out[counter] = 0; /* initialize result to zero */
      for(k=1;k<=p;k++){
	out[counter] = out[counter] + min_int(mat[i][k], mat[j][k]);
      }
      counter++;
    }
  }


  /* Free allocated memory */
  freeinttab(mat);

} /* end nb_shared_all */
