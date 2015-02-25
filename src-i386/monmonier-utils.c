/* 
A part of this code (intersection between segments) is based on Joseph O'Rourke code to identify intersection between two segments (originaly called 'segseg.c').
The original code was modified in order to handle double coordinates instead of integer.
A part of the code is new and devoted to some monmonier algorithm computations
To compile : R CMD SHLIB monmonier-utils.c

Thibaut Jombart (t.jombart@imperial.ac.uk), 2006, to fit his egocentric needs.

The original copyright follows.
*/ 


/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 7.  It is not written to be comprehensible without the 
explanation in that book.

Written by Joseph O'Rourke.
Last modified: November 1997
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1998 by Joseph O'Rourke.  It may be freely 
redistributed in its entirety provided that this copyright notice is 
not removed.
--------------------------------------------------------------------
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "adesub.h"



#define	EXIT_FAILURE 1
#define	X 0
#define	Y 1
typedef	enum { FALSE, TRUE } bool;

#define	DIM 2               /* Dimension of points */

typedef	double tPointd[DIM];   /* Type double point */
const double NEARZERO=10e-15; /* THRESHOLD for zero. */
/*---------------------------------------------------------------------
Function prototypes.
---------------------------------------------------------------------*/
int SegSeg( tPointd a, tPointd b, tPointd c, tPointd d);
int Parallel( tPointd a, tPointd b, tPointd c, tPointd d, tPointd p );
bool Between( tPointd a, tPointd b, tPointd c );
void Assignpx( tPointd p, tPointd a );
int Collinear( tPointd a, tPointd b, tPointd c );
int AreaSign( tPointd a, tPointd b, tPointd c );
void CheckAllSeg(int *nrow, int *ncol, double *tab, tPointd a, tPointd b, int *answer);
double dAbs(double a);
int dEqual(double a, double b);
/*-------------------------------------------------------------------*/

/* dAbs returns the absolute value of a double */
double dAbs(double a)
{
  if(a>=0.0)  return a;
  else return -a;
}

/* dEqual returns 1 if two doubles are equal, and 0 otherwise */
int dEqual(double a, double b)
{
  if(dAbs(a-b) < NEARZERO) return 1;
  else return 0;
}


/*-------------------------------------------------------------------
CheckAllSeg: performs the intersection test of a segment with 
a given set of segments. Calls SegSeg to perform 2-segments tests.
answer is returned back to R, and is:
0: no intersection
1: at least one intersection
-------------------------------------------------------------------*/
void CheckAllSeg(int *nrow, int *ncol, double *tab, tPointd a, tPointd b, int *answer)
{
/* Declarations de variables C locales */
/* X est le tableau des segments ; chaque ligne est un segment (xP,yP,xQ,yQ) */
int i,j,n,p,k,temp;
double **mat;
tPointd c,d;


/* Memory allocation for local C variables */
n = *nrow;
p = *ncol;

 taballoc(&mat, n, p); /* function from ade4 */

/* Reconstruction of the matrix of segments */
k = 0;
for (j=1; j<=p; j++) {
  for (i=1; i<=n; i++) {
            mat[i][j] = tab[k];
            k = k + 1;
  }
}

/* The segment of interest (ab) is checked for crossing against all other segments (cd)
We stop as soon as one segment is crossed */
 temp = 0;
 i = 1;
 while(temp ==0 && i<=n){
   c[X] = mat[i][1];
   c[Y] = mat[i][2];
   d[X] = mat[i][3];
   d[Y] = mat[i][4];
   temp = SegSeg(a,b,c,d);
   i++;
   }
 *answer = temp;

 /* Free allocated memory */
 freetab(mat);
}


/*---------------------------------------------------------------------
SegSeg: Tests an intersection between two closed
segments ab and cd. Returned values :
  3 : The segments collinearly overlap, sharing at least a point.
  2 : An endpoint (vertex) of one segment is on the other segment,
        but segments aren't collinear.
  1 : The segments intersect properly (i.e. not case 2 or 3)
  0 : The segments do not intersect.
  10 : initial value, i.e. failure.
---------------------------------------------------------------------*/
/* int	SegSeg( tPointd a, tPointd b, tPointd c, tPointd d, tPointd p) */
int	SegSeg( tPointd a, tPointd b, tPointd c, tPointd d)
{
   double  s, t;       /* The two parameters of the parametric eqns. */
   double num, denom;  /* Numerator and denoninator of equations. */
   int code = 10; /* returned value, default 10 is a failure */

   /* For debugging */
   /*printf("\n!!! SegSeg: code initialized at %d\n",code);*/

   /* Initialization of the intersection point 'p' */
   tPointd p;

   p[X] = -1;
   p[Y] = -1;

   denom = a[X] * (double)( d[Y] - c[Y] ) +
           b[X] * (double)( c[Y] - d[Y] ) +
           d[X] * (double)( b[Y] - a[Y] ) +
           c[X] * (double)( a[Y] - b[Y] );
   /* If denom is zero, then segments are parallel: handle separately. 
      Beware to avoid ... == 0 with doubles, 
      as well as ...==... */
   if (dAbs(denom) < NEARZERO) {
     code =  Parallel(a, b, c, d, p);
     /* For debugging */
     /*printf("\n!!! SegSeg: call to Parallel (denom=%f)\n",denom);*/
   }
   else{
     num =    a[X] * (double)( d[Y] - c[Y] ) +
       c[X] * (double)( a[Y] - d[Y] ) +
       d[X] * (double)( c[Y] - a[Y] );
     /* code 2 handled here */
     /*if ( ((num < NEARZERO) && (num > -NEARZERO)) || (num == denom) ) code = 2;*/
     if ( (dAbs(num) < NEARZERO) || (dEqual(num,denom)) ) code = 2;
     
     /* Debugging step 1*/
     /*printf("\n!!! SegSeg step1: dAbs(num)=%f, dEqual(num,denom)=%d), code=%d\n",dAbs(num),dEqual(num,denom),code);
     printf("\nNEARZERO=%f\n",NEARZERO);*/

     s = num / denom;
     
     num = -( a[X] * (double)( c[Y] - b[Y] ) +
	      b[X] * (double)( a[Y] - c[Y] ) +
	      c[X] * (double)( b[Y] - a[Y] ) );
     
     t = num / denom;
     
     /* if ( ((num < NEARZERO) && (num > -NEARZERO)) || (num == denom) ) code = 2;*/
     if ( (dAbs(num) < NEARZERO) || (dEqual(num,denom)) ) code = 2;
     /* Debugging step 2*/
     /*printf("\n!!! SegSeg step2: dAbs(num)=%f, dEqual(num,denom)=%d), code=%d\n",dAbs(num),dEqual(num,denom),code);
     printf("\nNEARZERO=%f\n",NEARZERO);*/
     
     if ( (NEARZERO < s) && (s < 1.0) &&
	  (NEARZERO < t) && (t < 1.0) )
       code = 1;
     else if ( (-NEARZERO > s) || (s > 1.0) ||
	       (-NEARZERO > t) || (t > 1.0) )
       code = 0;
     
     p[X] = a[X] + s * ( b[X] - a[X] );
     p[Y] = a[Y] + s * ( b[Y] - a[Y] );
   }
   
   /* Debugging step 3*/
   /*printf("\n!!! SegSeg step3: final value of code=%d\n",code);*/
   return code;
}

int	Parallel( tPointd a, tPointd b, tPointd c, tPointd d, tPointd p )
{
/* Avoid to consider segments as parallel whenever two points are the same. */
/*if ( (a[X]==b[X] && a[Y]==b[Y]) || (c[X]==d[X] && c[Y]==d[Y]) ) return 0;*/
if ( (dEqual(a[X],b[X]) && dEqual(a[Y],b[Y])) || (dEqual(c[X],d[X]) && dEqual(c[Y],d[Y])) ) return 0;
if ( Collinear( a, b, c)==0 ) return 0;

   if ( Between( a, b, c ) ) {
      Assignpx( p, c );
      return 3;
   }
   if ( Between( a, b, d ) ) {
      Assignpx( p, d );
      return 3;
   }
   if ( Between( c, d, a ) ) {
      Assignpx( p, a );
      return 3;
   }
   if ( Between( c, d, b ) ) {
      Assignpx( p, b );
      return 3;
   }
   return 0;
}
void	Assignpx( tPointd p, tPointd a )
{
   int i;
   for ( i = 0; i < DIM; i++ )
      p[i] = a[i];
}
/*---------------------------------------------------------------------
Returns TRUE iff point c lies on the closed segement ab.
Assumes it is already known that abc are collinear.
---------------------------------------------------------------------*/
bool    Between( tPointd a, tPointd b, tPointd c )
{
   /* If ab not vertical, check betweenness on x; else on y. */
   if ( a[X] != b[X] )
     return ((a[X] <= c[X]) && (c[X] <= b[X])) ||
       ((a[X] >= c[X]) && (c[X] >= b[X]));
   else
     return ((a[Y] <= c[Y]) && (c[Y] <= b[Y])) ||
       ((a[Y] >= c[Y]) && (c[Y] >= b[Y]));
}

int Collinear( tPointd a, tPointd b, tPointd c )
{
	if (AreaSign(a, b, c) ==0) return 1;
	else return 0;
}
int     AreaSign( tPointd a, tPointd b, tPointd c )
{
    double area2;

    area2 = ( b[X] - a[X] ) * (double)( c[Y] - a[Y] ) -
            ( c[X] - a[X] ) * (double)( b[Y] - a[Y] );

    /* The area is not an integer. */
    if      ( area2 >  NEARZERO ) return  1;
    else if ( area2 < -NEARZERO ) return -1;
    else                     return  0;
}
