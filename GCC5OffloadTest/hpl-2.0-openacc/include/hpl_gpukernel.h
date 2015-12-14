#ifndef HPL_GPUKERNEL_H
#define HPL_GPUKERNEL_H

void HPL_accdtrsm
(
   const enum HPL_ORDER             ORDER,
   const enum HPL_SIDE              SIDE,
   const enum HPL_UPLO              UPLO,
   const enum HPL_TRANS             TRANS,
   const enum HPL_DIAG              DIAG,
   const int                        M,
   const int                        N,
   const double                     ALPHA,
   const double *                   A,
   const int                        LDA,
   double *                         B,
   const int                        LDB
);

void HPL_accdgemm
(
   const enum HPL_ORDER             ORDER,
   const enum HPL_TRANS             TRANSA,
   const enum HPL_TRANS             TRANSB,
   const int                        M,
   const int                        N,
   const int                        K,
   const double                     ALPHA,
   const double *                   A,
   const int                        LDA,
   const double *                   B,
   const int                        LDB,
   const double                     BETA,
   double *                         C,
   const int                        LDC
);

#endif // HPL_GPUKERNEL_H
