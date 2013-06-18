#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

int dtodelta (double **data, double *pl);
void initvec (double *v1, double r);
double alea (void);
void aleapermutvec (double *a);
void aleapermutmat (double **a);
void aleapermutmat (double **a);
void aleapermutvec (double *a);
void DiagobgComp (int n0, double **w, double *d, int *rang);
void freeinttab (int **tab);
void freeintvec (int *vec);
void freetab (double **tab);
void freevec (double *vec);
void getpermutation (int *numero, int repet);
void matcentrage (double **A, double *poili, char *typ);
void matmodifcm (double **tab, double *poili);
void matmodifcn (double **tab, double *poili);
void matmodifcp (double **tab, double *poili);
void matmodifcs (double **tab, double *poili);
void matmodiffc (double **tab, double *poili);
void matpermut (double **A, int *num, double **B);
double maxvec (double *vec);
void prodmatAAtB (double **a, double **b);
void prodmatABC (double **a, double **b, double **c);
void prodmatAtAB (double **a, double **b);
void prodmatAtBC (double **a, double **b, double **c);
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut);
double traceXtdLXq (double **X, double **L, double *d, double *q);
void sqrvec (double *v1);
void taballoc (double ***tab, int l1, int c1);
void tabintalloc (int ***tab, int l1, int c1);
void trild (double *x , int *num, int gauche, int droite);
void trildintswap (int *v, int i, int j);
void trildswap (double *v, int i, int j);
void trirap (double *x , int *num);
void trirapideint (int *x , int *num, int gauche, int droite);
void trirapideintswap (int *v, int i, int j);
void vecalloc (double **vec, int n);
void vecintalloc (int **vec, int n);
void vecpermut (double *A, int *num, double *B);
