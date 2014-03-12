#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "adesub.h"

/***********************************************************************/
double traceXtdLXq (double **X, double **L, double *d, double *q)
/*  Produit matriciel XtDLXQ avec LX comme lag.matrix   */
{
    /* Declarations de variables C locales */
    int j, i, lig, col;
    double **auxi, **A, trace;
    
    
    
    /* Allocation memoire pour les variables C locales */
    lig = X[0][0];
    col = X[1][0];
    taballoc(&auxi, lig, col);
    taballoc(&A, col, col);
    
    
    /* Calcul de LX */
    prodmatABC(L, X, auxi);
    
    /* Calcul de DLX */
    for (i=1;i<=lig;i++) {
        for (j=1;j<=col;j++) {
            auxi[i][j] = auxi[i][j] * d[i];
        }       
    }
    
    /* Calcul de XtDLX */
    prodmatAtBC(X,auxi,A);
    
    /* Calcul de trace(XtDLXQ) */
    trace=0;
    for (i=1;i<=col;i++) {
        trace = trace + A[i][i] * q[i];
    }
    
    /* Libération des réservations locales */
    freetab (auxi);
    freetab (A);
    return(trace);
}

/***********************************************************************/
void tabintalloc (int ***tab, int l1, int c1)
/*--------------------------------------------------
* Allocation de memoire dynamique pour un tableau
* d'entiers (l1, c1)
--------------------------------------------------*/
{
    int     i, j;
    
    *tab = (int **) calloc(l1+1, sizeof(int *));

    if ( *tab != NULL) {
        for (i=0;i<=l1;i++) {
            
            *(*tab+i)=(int *) calloc(c1+1, sizeof(int));        
            if ( *(*tab+i) == NULL ) {
                for (j=0;j<i;j++) {
                    free(*(*tab+j));
                }
                return;
            }
        }
    } else return;
    **(*tab) = l1;
    **(*tab+1) = c1;
    for (i=1;i<=l1;i++) {
        for (j=1;j<=c1;j++) {
            (*tab)[i][j] = 0;
        }
    }
}

/***********************************************************************/
void freeinttab (int **tab)
/*--------------------------------------------------
* Allocation de memoire dynamique pour un tableau
--------------------------------------------------*/
{
    int     i, n;
    
    n = *(*(tab));
    
    for (i=0;i<=n;i++) {
            free((char *) *(tab+i) );
    }
    
    free((char *) tab);
}


/*********************/
int dtodelta (double **data, double *pl)
{
    /* la matrice de distances d2ij dans data est associee aux poids pl
    Elle est transformee par aij - ai. -a.j + a..
    aij = -d2ij/2);*/

    int lig, i, j;
    double *moy, a0, moytot;
    
    lig=data[0][0];
    vecalloc(&moy, lig);
    
    for (i=1; i<=lig; i++) {
        for (j=1; j<=lig; j++) data[i][j] = 0.0 - data[i][j] * data[i][j] / 2.0;
    } 

    for (i=1; i<=lig; i++) {
        a0=0;
        for (j=1; j<=lig; j++) a0 = a0 + pl[j]*data[i][j];
        moy[i] = a0;
    }
    moytot=0;
    for (i=1; i<=lig; i++) {
        moytot = moytot+pl[i]*moy[i];
    }
    for (i=1; i<=lig; i++) {
        for (j=1; j<=lig; j++) data[i][j] = data[i][j] - moy[i] - moy[j] + moytot;
    }
    freevec (moy);
    return 1;
}
/***************************/
void initvec (double *v1, double r)
/*--------------------------------------------------
* Initialisation des elements d'un vecteur
--------------------------------------------------*/
{
    int i, c1;
    
    c1 = v1[0];
    for (i=1;i<=c1;i++) {
        v1[i] = r;
    }
}
/**************************/
double alea (void)
{
    double w;
    w = ((double) rand())/ (double)RAND_MAX;
    return (w);
}
/*************************/
void aleapermutmat (double **a)
{
    /* permute au hasard les lignes du tableau a
    Manly p. 42 le tableau est modifié */
    int lig, i,j, col, n, k;
    double z;

    lig = a[0][0];
    col = a[1][0];
    for (i=1; i<=lig-1; i++) {
        j=lig-i+1;
        k = (int) (j*alea ()+1);
        /*k = (int) (j*genrand()+1);*/
        if (k>j) k=j;
        for (n=1; n<=col; n++) {
            z = a[j][n];
            a[j][n]=a[k][n];
            a[k][n] = z;
        }
    }
}
/*************************/
void aleapermutvec (double *a)
{
    /* permute au hasard les ÚlÚments du vecteur a
    Manly p. 42 Le vecteur est modifié
    from Knuth 1981 p. 139*/
    int lig, i,j, k;
    double z;
    
    lig = a[0];
    for (i=1; i<=lig-1; i++) {
        j=lig-i+1;
        k = (int) (j*alea()+1);
        /*k = (int) (j*genrand()+1);*/
        if (k>j) k=j;
        z = a[j];
        a[j]=a[k];
        a[k] = z;
    }
}
/***********************************************************************/
void DiagobgComp (int n0, double **w, double *d, int *rang)
/*--------------------------------------------------
* Diagonalisation
* T. FOUCART Analyse factorielle de tableaux multiples,
* Masson, Paris 1984,185p., p. 62. D'apr?s VPROP et TRIDI,
* de LEBART et coll.
--------------------------------------------------*/
{
    double          *s, epsilon;
    double          a, b, c, x, xp, q, bp, ab, ep, h, t, u , v;
    double          dble;
    int             ni, i, i2, j, k, jk, ijk, ij, l, ix, m, m1, isnou;
    
    vecalloc(&s, n0);
    a = 0.000000001;
    epsilon = 0.0000001;
    ni = 100;
    if (n0 == 1) {
        d[1] = w[1][1];
        w[1][1] = 1.0;
        *rang = 1;
        freevec (s);
        return;
    }
    
    for (i2=2;i2<=n0;i2++) {
        
        b=0.0;
        c=0.0;
        i=n0-i2+2;
        k=i-1;
        if (k < 2) goto Et1;
        for (l=1;l<=k;l++) {
            c = c + fabs((double) w[i][l]);
        }
        if (c != 0.0) goto Et2;
        
Et1:    s[i] = w[i][k];
        goto Etc;
        
Et2:    for (l=1;l<=k;l++) {
            x = w[i][l] / c;
            w[i][l] = x;
            b = b + x * x;
        }
        xp = w[i][k];
        ix = 1;
        if (xp < 0.0) ix = -1;
        
/*      q = -sqrt(b) * ix; */
        dble = b;
        dble = -sqrt(dble);
        q = dble * ix;

        s[i] = c * q;
        b = b - xp * q;
        w[i][k] = xp - q;
        xp = 0;
        for (m=1;m<=k;m++) {
            w[m][i] = w[i][m] / b / c;
            q = 0;
            for (l=1;l<=m;l++) {
                q = q + w[m][l] * w[i][l];
            }
            m1 = m + 1;
            if (k < m1) goto Et3;
            for (l=m1;l<=k;l++) {
                q = q + w[l][m] * w[i][l];
            }
            
Et3:        s[m] = q / b;
            xp = xp + s[m] * w[i][m];
        }
        bp = xp * 0.5 / b;
        for (m=1;m<=k;m++) {
            xp = w[i][m];
            q = s[m] - bp * xp;
            s[m] = q;
            for (l=1;l<=m;l++) {
                w[m][l] = w[m][l] - xp * s[l] - q * w[i][l];
            }
        }
        for (l=1;l<=k;l++) {
            w[i][l] = c * w[i][l];
        }
        
Etc:    d[i] = b;
    } /* for (i2=2;i2<n0;i2++) */
    
    s[1] = 0.0;
    d[1] = 0.0;
    
    for (i=1;i<=n0;i++) {
        
        k = i - 1;
        if (d[i] == 0.0) goto Et4;
        for (m=1;m<=k;m++) {
            q = 0.0;
            for (l=1;l<=k;l++) {
                q = q + w[i][l] * w[l][m];
            }
            for (l=1;l<=k;l++) {
                w[l][m] = w[l][m] - q * w[l][i];
            }
        }
        
Et4:    d[i] = w[i][i];
        w[i][i] = 1.0;
        if (k < 1) goto Et5;
        for (m=1;m<=k;m++) {
            w[i][m] = 0.0;
            w[m][i] = 0.0;
        }

Et5:;
    }
    
    for (i=2;i<=n0;i++) {
        s[i-1] = s[i];
    }
    s[n0] = 0.0;
    
    for (k=1;k<=n0;k++) {

        m = 0;

Et6:    for (j=k;j<=n0;j++) {
            if (j == n0) goto Et7;
            ab = fabs((double) s[j]);
            ep = a * (fabs((double) d[j]) + fabs((double) d[j+1]));
            if (ab < ep) goto Et7;
        }
    
Et7:    isnou = 1;
        h = d[k];
        if (j == k) goto Eta;
        if (m < ni) goto Etd;
        
        /*err_message("Error: can't compute matrix eigenvalues");*/
        
Etd:    m = m + 1;
        q = (d[k+1]-h) * 0.5 / s[k];
        
/*      t = sqrt(q * q + 1.0); */
        dble = q * q + 1.0;
        dble = sqrt(dble);
        t = dble;
        
        if (q < 0.0) isnou = -1;
        q = d[j] - h + s[k] / (q + t * isnou);
        u = 1.0;
        v = 1.0;
        h = 0.0;
        jk = j-k;
        for (ijk=1;ijk<=jk;ijk++) {
            i = j - ijk;
            xp = u * s[i];
            b = v * s[i];
            if (fabs((double) xp) < fabs((double) q)) goto Et8;
            u = xp / q;
            
/*          t = sqrt(u * u + 1); */
            dble = u * u + 1.0;
            dble = sqrt(dble);
            t = dble;
            
            s[i+1] = q * t;
            v = 1 / t;
            u = u * v;
            goto Et9;

Et8:        v = q / xp;

/*          t = sqrt(1 + v * v); */
            dble = 1.0 + v * v;
            dble = sqrt(dble);
            t = dble;
            
            s[i+1] = t * xp;
            u = 1 / t;
            v = v * u;

Et9:
            q = d[i+1] - h;
            t = (d[i] - q) * u + 2.0 * v * b;
            h = u * t;
            d[i+1] = q + h;
            q = v * t - b;
            for (l=1;l<=n0;l++) {
                xp = w[l][i+1];
                w[l][i+1] = u * w[l][i] + v * xp;
                w[l][i] = v * w[l][i] - u * xp;
            }
        }
        d[k] = d[k] - h;
        s[k] = q;
        s[j] = 0.0;
        
        goto Et6;

Eta:;
    } /* for (k=1;k<=n0;k++) */
    
    for (ij=2;ij<=n0;ij++) {
        
        i = ij - 1;
        l = i;
        h = d[i];
        for (m=ij;m<=n0;m++) {
            if (d[m] >= h) {
                l = m;
                h = d[m];
            }
        }
        if (l == i) {
            goto Etb;
        } else {
            d[l] = d[i];
            d[i] = h;
        }
        for (m=1;m<=n0;m++) {
            h = w[m][i];
            w[m][i] = w[m][l];
            w[m][l] = h;
        }

Etb:;
    } /* for (ij=2;ij<=n0;ij++) */

    *rang = 0;
    for (i=1;i<=n0;i++) {
        if (d[i] / d[1] < epsilon) d[i] = 0.0;
        if (d[i] != 0.0) *rang = *rang + 1;
    }
    freevec(s);
} /* DiagoCompbg */
/***********************************************************************/
void freeintvec (int *vec)
/*--------------------------------------------------
* liberation de memoire pour un vecteur
--------------------------------------------------*/
{
    
    free((char *) vec);
    
}
/***********************************************************************/
void freetab (double **tab)
/*--------------------------------------------------
* Allocation de memoire dynamique pour un tableau (l1, c1)
--------------------------------------------------*/
{
    int     i, n;
    
    n = *(*(tab));
    for (i=0;i<=n;i++) {
            free((char *) *(tab+i) );
    }
    free((char *) tab);
}
/***********************************************************************/
void freevec (double *vec)
/*--------------------------------------------------
* liberation de memoire pour un vecteur
--------------------------------------------------*/
{
    free((char *) vec); 
}
/***********************************************************************/
void getpermutation (int *numero, int repet)
/*----------------------
* affectation d'une permutation alÚatoire des n premiers entiers 
* dans dans un vecteur d'entiers de dimension n
* vecintalloc prÚalable exigÚ
* *numero est un vecteur d'entier
* repet est un entier qui peut prendre une valeur arbitraire
* utilise dans le germe du generateur de nb pseudo-aleatoires
* si on l'incremente dans des appels repetes (e.g. simulation) garantit
* que deux appels donnent deux resultats distincts (seed=clock+repet)
------------------------*/
{
    int i, n, seed;
    int *alea;
    
    n=numero[0];
    vecintalloc (&alea,n);
    
    /*-------------
    * numerotation dans numero
    -----------*/
    for (i=1;i<=n;i++) {
        numero[i]=i;
    }
    
    /*-------------
    * affectation de nombres aleatoires dans alea
    ----------------*/
    seed = clock();
    seed = seed + repet;
    srand(seed);
    for (i=1;i<=n;i++) {
        alea[i]=rand();
    }
    
    trirapideint (alea , numero, 1, n);
    freeintvec (alea);
}
/***********************************************************************/
void matcentrage (double **A, double *poili, char *typ)
{
    
    if (strcmp (typ,"nc") == 0) {
        return;
    } else if (strcmp (typ,"cm") == 0) {
        matmodifcm (A, poili);
        return;
    } else if (strcmp (typ,"cn") == 0) {
        matmodifcn (A, poili);
        return;
    } else if (strcmp (typ,"cp") == 0) {
        matmodifcp (A, poili);
        return;
    } else if (strcmp (typ,"cs") == 0) {
        matmodifcs (A, poili);
        return;
    } else if (strcmp (typ,"fc") == 0) {
        matmodiffc (A, poili);
        return;
    } else if (strcmp (typ,"fl") == 0) {
        matmodifcm (A, poili);
        return;
    }
}
/***********************************************************************/
void matmodifcm (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, m colonnes
* disjonctif complet
* poili est un vecteur n composantes
* la procedure retourne tab centre par colonne 
* pour la ponderation poili (somme=1)
* centrage type correspondances multiples
--------------------------------------------------*/
{
    double      poid;
    int             i, j, l1, m1;
    double      *poimoda;
    double      x, z;

    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);


    for (i=1;i<=l1;i++) {
        poid = poili[i];
        for (j=1;j<=m1;j++) {
            poimoda[j] = poimoda[j] + tab[i][j] * poid;
        }
    }
    
    for (j=1;j<=m1;j++) {
        x = poimoda[j];
        if (x==0) {
            for (i=1;i<=l1;i++) tab[i][j] = 0;
        } else {
        
            for (i=1;i<=l1;i++) {
                z = tab[i][j]/x - 1.0;
                tab[i][j] = z;
            }
        }
    }
    freevec (poimoda);
}
/***********************************************************************/
void matmodifcn (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, p colonnes
* poili est un vecteur n composantes
* la procedure retourne tab norme par colonne 
* pour la ponderation poili (somme=1)
--------------------------------------------------*/
{
    double      poid, x, z, y, v2;
    int             i, j, l1, c1;
    double      *moy, *var;

    l1 = tab[0][0];
    c1 = tab[1][0];

    vecalloc(&moy, c1);
    vecalloc(&var, c1);


/*--------------------------------------------------
* calcul du tableau centre/norme
--------------------------------------------------*/

    for (i=1;i<=l1;i++) {
        poid = poili[i];
        for (j=1;j<=c1;j++) {
            moy[j] = moy[j] + tab[i][j] * poid;
        }
    }
    
    for (i=1;i<=l1;i++) {
        poid=poili[i];
        for (j=1;j<=c1;j++) {
            x = tab[i][j] - moy[j];
            var[j] = var[j] + poid * x * x;
        }
    }
    
    for (j=1;j<=c1;j++) {
        v2 = var[j];
        if (v2<=0) v2 = 1;
        v2 = sqrt(v2);
        var[j] = v2;
    }
    
    for (i=1;i<=c1;i++) {
        x = moy[i];
        y = var[i];
        for (j=1;j<=l1;j++) {
            z = tab[j][i] - x;
            z = z / y;
            tab[j][i] = z;
        }
    }
    
    freevec(moy);
    freevec(var);
    
}
/***********************************************************************/
void matmodifcs (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, p colonnes
* poili est un vecteur n composantes
* la procedure retourne tab standardise par colonne 
* pour la ponderation poili (somme=1)
--------------------------------------------------*/
{
    double      poid, x, z, y, v2;
    int         i, j, l1, c1;
    double      *var;

    l1 = tab[0][0];
    c1 = tab[1][0];

    vecalloc(&var, c1);


/*--------------------------------------------------
* calcul du tableau standardise
--------------------------------------------------*/

    for (i=1;i<=l1;i++) {
        poid=poili[i];
        for (j=1;j<=c1;j++) {
            x = tab[i][j];
            var[j] = var[j] + poid * x * x;
        }
    }
    
    for (j=1;j<=c1;j++) {
        v2 = var[j];
        if (v2<=0) v2 = 1;
        v2 = sqrt(v2);
        var[j] = v2;
    }
    
    for (i=1;i<=c1;i++) {
        y = var[i];
        for (j=1;j<=l1;j++) {
            z = tab[j][i];
            z = z / y;
            tab[j][i] = z;
        }
    }
    freevec(var);
}
/***********************************************************************/
void matmodifcp (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, p colonnes
* poili est un vecteur n composantes
* la procedure retourne tab centre par colonne 
* pour la ponderation poili (somme=1)
--------------------------------------------------*/
{
    double      poid;
    int             i, j, l1, c1;
    double      *moy, x, z;

    l1 = tab[0][0];
    c1 = tab[1][0];
    vecalloc(&moy, c1);


/*--------------------------------------------------
* calcul du tableau centre
--------------------------------------------------*/

    for (i=1;i<=l1;i++) {
        poid = poili[i];
        for (j=1;j<=c1;j++) {
            moy[j] = moy[j] + tab[i][j] * poid;
        }
    }
    
    
    for (i=1;i<=c1;i++) {
        x = moy[i];
        for (j=1;j<=l1;j++) {
            z = tab[j][i] - x;
            tab[j][i] = z;
        }
    }
    freevec(moy);
}
/***********************************************************************/
void matmodiffc (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, m colonnes
* de nombres positifs ou nuls
* poili est un vecteur n composantes
* la procedure retourne tab centre doublement 
* pour la ponderation poili (somme=1)
* centrage type correspondances simples
--------------------------------------------------*/
{
    double      poid;
    int             i, j, l1, m1;
    double      *poimoda;
    double      x, z;

    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);


    for (i=1;i<=l1;i++) {
        x = 0;
        for (j=1;j<=m1;j++) {
            x = x + tab[i][j];
        }
        if (x!=0) {
            for (j=1;j<=m1;j++) {
                tab[i][j] = tab[i][j]/x;
            }
        }   
    }

    for (i=1;i<=l1;i++) {
        poid = poili[i];
        for (j=1;j<=m1;j++) {
            poimoda[j] = poimoda[j] + tab[i][j] * poid;
        }
    }
    
    for (j=1;j<=m1;j++) {
        x = poimoda[j];
        if (x==0) {
            /*err_message("column has a nul weight (matmodiffc)");*/
        }
        
        for (i=1;i<=l1;i++) {
            z = tab[i][j]/x - 1.0;
            tab[i][j] = z;
        }
    }
    freevec (poimoda);
}
/***********************************************************************/
void matpermut (double **A, int *num, double **B)
{
/*---------------------------------------
* A est une matrice n-p
* B est une matrice n-p
* num est une permutation aléatoire des n premiers entiers
* B contient en sortie les lignes de A permutÚes
* ---------------------------------------*/

    int lig, col,lig1, col1, lig2, i, j, k;
    
    lig = A[0][0];
    col = A[1][0];
    lig1 = B[0][0];
    col1 = B[1][0];
    lig2 = num[0];
    
    
    if ( (lig!=lig1) || (col!=col1) || (lig!=lig2) ) {
        return;
    }
    
    for (i=1; i<=lig; i++) {
        k=num[i];
        for (j=1; j<=col; j++) {
            B[i][j] = A[k][j];
        }
    }
}
/***********************************************************************/
void prodmatABC (double **a, double **b, double **c)
/*--------------------------------------------------
* Produit matriciel AB
--------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];

    for (i=1;i<=lig;i++) {
        for (k=1;k<=col2;k++) {
            s = 0;
            for (j=1;j<=col;j++) {
                s = s + a[i][j] * b[j][k];
            }
        c[i][k] = s;
        }       
    }
}

/***********************************************************************/
void prodmatAtAB (double **a, double **b)
/*--------------------------------------------------
* Produit matriciel AtA
--------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];

    for (j=1;j<=col;j++) {
        for (k=j;k<=col;k++) {
            s = 0;
            for (i=1;i<=lig;i++) {
                s = s + a[i][k] * a[i][j];
            }
        b[j][k] = s;
        b[k][j] = s;
        }       
    }
}
/***********************************************************************/
void prodmatAtBC (double **a, double **b, double **c)
/*--------------------------------------------------
* Produit matriciel AtB
--------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];

    for (j=1;j<=col;j++) {
        for (k=1;k<=col2;k++) {
            s = 0;
            for (i=1;i<=lig;i++) {
                s = s + a[i][j] * b[i][k];
            }
        c[j][k] = s;
        }       
    }
}
/***********************************************************************/
double maxvec (double *vec)
/*--------------------------------------------------
* calcul le max d'un vecteur
--------------------------------------------------*/
{
    int     i, len;
    double x;
    
    x = vec[1];
    len = vec[0];
    for (i=1;i<=len;i++) {
        if (vec[i] > x) x = vec[i];
    }
    return(x);
}
/***********************************************************************/
void prodmatAAtB (double **a, double **b)
/*--------------------------------------------------
* Produit matriciel B = AAt
--------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];

    for (j=1;j<=lig;j++) {
        for (k=j;k<=lig;k++) {
            s = 0;
            for (i=1;i<=col;i++) {
                s = s + a[j][i] * a[k][i];
            }
        b[j][k] = s;
        b[k][j] = s;
        }       
    }
}
/***********************************************************************/
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut)
/*--------------------------------------------------
* Produit matriciel AtB
* les lignes de B sont permutÚes par la permutation permut
--------------------------------------------------*/
{
    int j, k, i, i0, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];

    for (j=1;j<=col;j++) {
        for (k=1;k<=col2;k++) {
            s = 0;
            for (i=1;i<=lig;i++) {
                i0 = permut[i];
                s = s + a[i][j] * b[i0][k];
            }
        c[j][k] = s;
        }       
    }
}
/***********************************************************************/
void sqrvec (double *v1)
/*--------------------------------------------------
* Racine carree des elements d'un vecteur
--------------------------------------------------*/
{
    int i, c1;
    double v2;
    
    c1 = v1[0];
    
    for (i=1;i<=c1;i++) {
        v2 = v1[i];
        /* if (v2 < 0.0) err_message("Error: Square root of negative number (sqrvec)");*/
        v2 = sqrt(v2);
        v1[i] = v2;
    }
}
/***********************************************************************/
void taballoc (double ***tab, int l1, int c1)
/*--------------------------------------------------
* Allocation de memoire dynamique pour un tableau (l1, c1)
--------------------------------------------------*/
{
    int i, j;
    
    if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != 0) {
        for (i=0;i<=l1;i++) {
            if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == 0 ) {
                return;
                for (j=0;j<i;j++) {
                    free(*(*tab+j));
                }
            }
        }
    }

    **(*tab) = l1;
    **(*tab+1) = c1;
}
/***********************************************************************/
void trild (double *x , int *num, int gauche, int droite)
/*--------------------------------------------------
* Tri d'un tableau de double avec conservation du rang
* dans un tableau entier.
--------------------------------------------------*/
{
    int j, dernier, milieu;
    double  t;

    
    if ( (droite-gauche)<=0) return;
    milieu = (gauche+droite)/2;
    trildswap (x, gauche, milieu);
    
    trildintswap (num, gauche, milieu);
    t=x[gauche];
    dernier=gauche;
    for (j = gauche+1; j<=droite; j++) {
        if (x[j] > t) {
            dernier = dernier + 1;
            trildswap (x, dernier, j);  
            trildintswap (num, dernier, j);
        }
    }
    trildswap (x, gauche, dernier);
    trildintswap (num, gauche, dernier);
    trild (x, num, gauche, dernier-1);
    trild (x, num, dernier+1, droite);
}
/**************************/
void trildintswap (int *v, int i, int j)
{
    int provi;
    
    provi=v[i];
    v[i]=v[j];
    v[j]=provi;
}
/***********************************************************************/
void trildswap (double *v, int i, int j)
/*--------------------------------------------------
* Echange les valeurs de deux double
--------------------------------------------------*/
{
    double provi;
    
    provi=v[i];
    v[i]=v[j];
    v[j]=provi;
}

/***********************************************************************/
void trirap (double *x , int *num)
/*--------------------------------------------------
* Tri d'un tableau de double par ordre croissant
* avec conservation du rang dans un tableau entier.
--------------------------------------------------*/
{
    int             i, n, *num2, gauche, droite;
    double      *x2;
    
    n = x[0];
    gauche = 1;
    droite = n;
    vecalloc(&x2, n);
    vecintalloc(&num2, n);
    for (i=1;i<=n;i++) num[i] = i;
    trild(x, num, gauche, droite);
    for (i=1;i<=n;i++) {
        x2[i] = x[n - i + 1];
        num2[i] = num[n - i + 1];
    }
    for (i=1;i<=n;i++) {
        x[i] = x2[i];
        num[i] = num2[i];
    }
    freevec(x2);
    freeintvec(num2);
}
/***********************************************************************/
void trirapideint (int *x , int *num, int gauche, int droite)
{
    int j, dernier, milieu, t;
    
    if ( (droite-gauche)<=0) return;
    
    milieu = (gauche+droite)/2;
    trirapideintswap (x, gauche, milieu);
    trirapideintswap (num, gauche, milieu);
    
    t=x[gauche];
    dernier=gauche;
    for (j = gauche+1; j<=droite; j++) {
        if (x[j] < t) {
            dernier = dernier + 1;
            trirapideintswap (x, dernier, j);   
            trirapideintswap (num, dernier, j);
        }
    }
    trirapideintswap (x, gauche, dernier);
    trirapideintswap (num, gauche, dernier);
    
    trirapideint (x, num, gauche, dernier-1);
    trirapideint (x, num, dernier+1, droite);
        
}
/***********************************************************************/
void trirapideintswap (int *v, int i, int j)
{
    int provi;
    
    provi=v[i];
    v[i]=v[j];
    v[j]=provi;
}
/***********************************************************************/
void vecalloc (double **vec, int n)
/*--------------------------------------------------
* Allocation de memoire pour un vecteur de longueur n
--------------------------------------------------*/
{
    if ( (*vec = (double *) calloc(n+1, sizeof(double))) != 0) {
        **vec = n;
        return;
    } else {
        return;
    }
}
/***********************************************************************/
void vecintalloc (int **vec, int n)
/*--------------------------------------------------
* Allocation de memoire pour un vecteur d'entiers de longueur n
--------------------------------------------------*/
{
    if ( (*vec = (int *) calloc(n+1, sizeof(int))) != NULL) {
        **vec = n;
        return;
    } else {
        return;
    }
}



/***********************************************************************/
void vecpermut (double *A, int *num, double *B)
{
/*---------------------------------------
* A est un vecteur n elements
* B est une vecteur n elements
* num est une permutation alŽatoire des n premiers entiers
* B contient en sortie les elements de A permutŽes
* ---------------------------------------*/

    int lig, lig1, lig2, i, k;
    
    lig = A[0];
    lig1 = B[0];
    lig2 = num[0];
    
    
    if ( (lig!=lig1) || (lig!=lig2) ) {
        /*err_message ("Illegal parameters (vecpermut)");
        closelisting();*/
    }
    
    for (i=1; i<=lig; i++) {
        k=num[i];
        B[i] = A[k];
    }
}
