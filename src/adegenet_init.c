/* #include <stdlib.h> // for NULL */
/* #include <R_ext/Rdynload.h> */

/* /\* FIXME:  */
/*    Check these declarations against the C/Fortran source code. */
/* *\/ */

/* /\* .C calls *\/ */
/* extern void binIntToBytes(void *, void *, void *, void *); */
/* extern void bytesToBinInt(void *, void *, void *); */
/* extern void bytesToInt(void *, void *, void *, void *, void *); */
/* extern void CheckAllSeg(void *, void *, void *, void *, void *, void *); */
/* extern void GLdotProd(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *); */
/* extern void GLsumFreq(void *, void *, void *, void *, void *, void *, void *, void *, void *); */
/* extern void GLsumInt(void *, void *, void *, void *, void *, void *, void *, void *, void *); */
/* extern void nb_shared_all(void *, void *, void *, void *); */

/* static const R_CMethodDef CEntries[] = { */
/*     {"binIntToBytes", (DL_FUNC) &binIntToBytes,  4}, */
/*     {"bytesToBinInt", (DL_FUNC) &bytesToBinInt,  3}, */
/*     {"bytesToInt",    (DL_FUNC) &bytesToInt,     5}, */
/*     {"CheckAllSeg",   (DL_FUNC) &CheckAllSeg,    6}, */
/*     {"GLdotProd",     (DL_FUNC) &GLdotProd,     12}, */
/*     {"GLsumFreq",     (DL_FUNC) &GLsumFreq,      9}, */
/*     {"GLsumInt",      (DL_FUNC) &GLsumInt,       9}, */
/*     {"nb_shared_all", (DL_FUNC) &nb_shared_all,  4}, */
/*     {NULL, NULL, 0} */
/* }; */

/* void R_init_adegenet(DllInfo *dll) */
/* { */
/*     R_registerRoutines(dll, CEntries, NULL, NULL, NULL); */
/*     R_useDynamicSymbols(dll, FALSE); */
/* } */
