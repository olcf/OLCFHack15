#include "hpl.h"
#include "hpl_gpukernel.h"


static void HPL_accdtrsmLLNN
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, iaik, ibij, ibkj, j, jak, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, ibij= jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = 0, jak  = 0, ibkj = jbj; k < M; k++, jak += LDA, ibkj += 1 )
      {
         B[ibkj] /= A[k+jak];
         for( i = k+1,    iaik  = k+1+jak, ibij  = k+1+jbj;
              i < M; i++, iaik +=1,        ibij += 1 )
         { B[ibij] -= B[ibkj] * A[iaik]; }
      }
   }
}

static void HPL_accdtrsmLLNU
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, iaik, ibij, ibkj, j, jak, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, ibij= jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = 0, jak  = 0, ibkj = jbj; k < M; k++, jak += LDA, ibkj += 1 )
      {
         for( i = k+1,    iaik  = k+1+jak, ibij  = k+1+jbj;
              i < M; i++, iaik +=1,        ibij += 1 )
         { B[ibij] -= B[ibkj] * A[iaik]; }
      }
   }
}

static void HPL_accdtrsmLLTN
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   register double            t0;
   int                        i, iaki, ibij, ibkj, j, jai, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = M-1,     jai  = (M-1)*LDA, ibij  = M-1+jbj;
           i >= 0; i--, jai -= LDA,       ibij -= 1 )
      {
         t0 = ALPHA * B[ibij];
         for( k = i+1,    iaki  = i+1+jai, ibkj  = i+1+jbj;
              k < M; k++, iaki += 1,       ibkj += 1 )
         { t0 -= A[iaki] * B[ibkj]; }
         t0 /= A[i+jai];
         B[ibij] = t0;
      }
   }
}

static void HPL_accdtrsmLLTU
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   register double            t0;
   int                        i, iaki, ibij, ibkj, j, jai, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = M-1,     jai  = (M-1)*LDA, ibij  = M-1+jbj;
           i >= 0; i--, jai -= LDA,       ibij -= 1 )
      {
         t0 = ALPHA * B[ibij];
         for( k = i+1,    iaki  = i+1+jai, ibkj  = i+1+jbj;
              k < M; k++, iaki += 1,       ibkj += 1 )
         { t0 -= A[iaki] * B[ibkj]; }
         B[ibij] = t0;
      }
   }
}

static void HPL_accdtrsmLUNN
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, iaik, ibij, ibkj, j, jak, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = M-1,     jak  = (M-1)*LDA, ibkj  = M-1+jbj;
           k >= 0; k--, jak -= LDA,       ibkj -= 1 )
      {
         B[ibkj] /= A[k+jak];
         for( i = 0,      iaik  = jak, ibij  = jbj;
              i < k; i++, iaik += 1,   ibij += 1 )
         { B[ibij] -= B[ibkj] * A[iaik]; }
      }
   }
}

static void HPL_accdtrsmLUNU
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, iaik, ibij, ibkj, j, jak, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = M-1,     jak  = (M-1)*LDA, ibkj  = M-1+jbj;
           k >= 0; k--, jak -= LDA,       ibkj -= 1 )
      {
         for( i = 0,      iaik  = jak, ibij  = jbj;
              i < k; i++, iaik += 1,   ibij += 1 )
         { B[ibij] -= B[ibkj] * A[iaik]; }
      }
   }
}

static void HPL_accdtrsmLUTN
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, iaki, ibij, ibkj, j, jai, jbj, k;
   register double            t0;

   for( j = 0, jbj  = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, jai  = 0, ibij = jbj; i < M; i++, jai += LDA, ibij += 1 )
      {
         t0 = ALPHA * B[ibij];
         for( k = 0, iaki = jai, ibkj = jbj; k < i; k++, iaki += 1, ibkj += 1 )
         { t0 -= A[iaki] * B[ibkj]; }
         t0 /= A[i+jai];
         B[ibij] = t0;
      }
   }
}

static void HPL_accdtrsmLUTU
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   register double            t0;
   int                        i, iaki, ibij, ibkj, j, jai, jbj, k;

   for( j = 0, jbj  = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, jai  = 0, ibij = jbj; i < M; i++, jai += LDA, ibij += 1 )
      {
         t0 = ALPHA * B[ibij];
         for( k = 0, iaki = jai, ibkj = jbj; k < i; k++, iaki += 1, ibkj += 1 )
         { t0 -= A[iaki] * B[ibkj]; }
         B[ibij] = t0;
      }
   }
}

static void HPL_accdtrsmRLNN
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

   for( j = N-1,      jaj  = (N-1)*LDA, jbj  = (N-1)*LDB;
        j >= 0;  j--, jaj -= LDA,       jbj -= LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = j+1,    iakj  = j+1+jaj, jbk  = (j+1)*LDB;
           k < N; k++, iakj += 1,       jbk += LDB )
      {
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= A[iakj] * B[ibik]; }
      }
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] /= A[j+jaj]; }
   }
}

static void HPL_accdtrsmRLNU
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

   for( j = N-1,      jaj  = (N-1)*LDA, jbj  = (N-1)*LDB;
        j >= 0;  j--, jaj -= LDA,       jbj -= LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = j+1,    iakj  = j+1+jaj, jbk  = (j+1)*LDB;
           k < N; k++, iakj += 1,       jbk += LDB )
      {
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= A[iakj] * B[ibik]; }
      }
   }
}

static void HPL_accdtrsmRLTN
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   register double            t0;
   int                        i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

   for( k = 0, jak = 0, jbk = 0; k < N; k++, jak += LDA, jbk += LDB )
   {
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] /= A[k+jak]; }
      for( j = k+1,    iajk  = (k+1)+jak, jbj  = (k+1)*LDB;
           j < N; j++, iajk += 1,         jbj += LDB )
      {
         t0 = A[iajk];
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= t0 * B[ibik]; }
      }
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] *= ALPHA; }
   }
}

static void HPL_accdtrsmRLTU
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   register double            t0;
   int                        i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

   for( k = 0, jak = 0, jbk = 0; k < N; k++, jak += LDA, jbk += LDB )
   {
      for( j = k+1,    iajk  = (k+1)+jak, jbj  = (k+1)*LDB;
           j < N; j++, iajk += 1,         jbj += LDB )
      {
         t0 = A[iajk];
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= t0 * B[ibik]; }
      }
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] *= ALPHA; }
   }
}

static void HPL_accdtrsmRUNN
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

   for( j = 0, jaj = 0, jbj = 0; j < N; j++, jaj += LDA, jbj += LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = 0, iakj = jaj, jbk = 0; k < j; k++, iakj += 1, jbk += LDB )
      {
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= A[iakj] * B[ibik]; }
      }
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] /= A[j+jaj]; }
   }
}

static void HPL_accdtrsmRUNU
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

   for( j = 0, jaj = 0, jbj = 0; j < N; j++, jaj += LDA, jbj += LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = 0, iakj = jaj, jbk = 0; k < j; k++, iakj += 1, jbk += LDB )
      {
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= A[iakj] * B[ibik]; }
      }
   }
}

static void HPL_dtrsmRUTN
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   register double            t0;
   int                        i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

   for( k = N-1,     jak  = (N-1)*LDA, jbk  = (N-1)*LDB;
        k >= 0; k--, jak -= LDA,       jbk -= LDB )
   {
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] /= A[k+jak]; }
      for( j = 0, iajk = jak, jbj = 0; j < k; j++, iajk += 1, jbj += LDB )
      {
         t0 = A[iajk];
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= t0 * B[ibik]; }
      }
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] *= ALPHA; }
   }
}

static void HPL_accdtrsmRUTU
(
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   register double            t0;
   int                        i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

   for( k = N-1,     jak  = (N-1)*LDA, jbk  = (N-1)*LDB;
        k >= 0; k--, jak -= LDA,       jbk -= LDB )
   {
      for( j = 0, iajk = jak, jbj = 0; j < k; j++, iajk += 1, jbj += LDB )
      {
         t0 = A[iajk];
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= t0 * B[ibik]; }
      }
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] *= ALPHA; }
   }
}

static void HPL_accdtrsm0
(
   const enum HPL_SIDE        SIDE,
   const enum HPL_UPLO        UPLO,
   const enum HPL_TRANS       TRANS,
   const enum HPL_DIAG        DIAG,
   const int                  M,
   const int                  N,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   double                     * B,
   const int                  LDB
)
{
   int                        i, j;

   if( ( M == 0 ) || ( N == 0 ) ) return;

   if( ALPHA == HPL_rzero )
   {
      for( j = 0; j < N; j++ )
      {  for( i = 0; i < M; i++ ) *(B+i+j*LDB) = HPL_rzero; }
      return;
   }

   if( SIDE == HplLeft )
   {
      if( UPLO == HplUpper )
      {
         if( TRANS == HplNoTrans )
         {
            if( DIAG == HplNonUnit )
            {      HPL_accdtrsmLUNN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPL_accdtrsmLUNU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
         else
         {
            if( DIAG == HplNonUnit )
            {      HPL_accdtrsmLUTN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPL_accdtrsmLUTU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
      }
      else
      {
         if( TRANS == HplNoTrans )
         {
            if( DIAG == HplNonUnit )
            {      HPL_accdtrsmLLNN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPL_accdtrsmLLNU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
         else
         {
            if( DIAG == HplNonUnit )
            {      HPL_accdtrsmLLTN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPL_accdtrsmLLTU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
      }
   }
   else
   {
      if( UPLO == HplUpper )
      {
         if( TRANS == HplNoTrans )
         {
            if( DIAG == HplNonUnit )
            {      HPL_accdtrsmRUNN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPL_accdtrsmRUNU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
         else
         {
            if( DIAG == HplNonUnit )
            {      HPL_dtrsmRUTN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPL_accdtrsmRUTU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
      }
      else
      {
         if( TRANS == HplNoTrans )
         {
            if( DIAG == HplNonUnit )
            {      HPL_accdtrsmRLNN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPL_accdtrsmRLNU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
         else
         {
            if( DIAG == HplNonUnit )
            {      HPL_accdtrsmRLTN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPL_accdtrsmRLTU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
      }
   }
}


static void HPL_accdscal
(
   const int                        N,
   const double                     ALPHA,
   double *                         X,
   const int                        INCX
)
{
   register double           x0, x1, x2, x3, x4, x5, x6, x7;
   register const double     alpha = ALPHA;
   const double              * StX;
   register int              i;
   int                       nu;
   const int                 incX2 = 2 * INCX, incX3 = 3 * INCX,
                             incX4 = 4 * INCX, incX5 = 5 * INCX,
                             incX6 = 6 * INCX, incX7 = 7 * INCX,
                             incX8 = 8 * INCX;

   if( ( N > 0 ) && ( alpha != HPL_rone ) )
   {
      if( alpha == HPL_rzero )
      {
         if( ( nu = ( N >> 3 ) << 3 ) != 0 )
         {
            StX = (double *)X + nu * INCX;

            do
            {
               (*X)     = HPL_rzero; X[incX4] = HPL_rzero;
               X[INCX ] = HPL_rzero; X[incX5] = HPL_rzero;
               X[incX2] = HPL_rzero; X[incX6] = HPL_rzero;
               X[incX3] = HPL_rzero; X[incX7] = HPL_rzero; X += incX8;

            } while( X != StX );
         }

         for( i = N - nu; i != 0; i-- ) { *X = HPL_rzero; X += INCX; }
      }
      else
      {
         if( ( nu = ( N >> 3 ) << 3 ) != 0 )
         {
            StX = X + nu * INCX;

            do
            {
               x0 = (*X);     x4 = X[incX4]; x1 = X[INCX ]; x5 = X[incX5];
               x2 = X[incX2]; x6 = X[incX6]; x3 = X[incX3]; x7 = X[incX7];

               x0 *= alpha;   x4 *= alpha;   x1 *= alpha;   x5 *= alpha;
               x2 *= alpha;   x6 *= alpha;   x3 *= alpha;   x7 *= alpha;

               (*X)     = x0; X[incX4] = x4; X[INCX ] = x1; X[incX5] = x5;
               X[incX2] = x2; X[incX6] = x6; X[incX3] = x3; X[incX7] = x7;

               X  += incX8;

            } while( X != StX );
         }

         for( i = N - nu; i != 0; i-- )
         { x0 = (*X); x0 *= alpha; *X = x0; X += INCX; }
      }
   }
}

#define DSCAL(N,ALPHA,X,INCX) \
{ \
   register double           x0, x1, x2, x3, x4, x5, x6, x7; \
   register const double     alpha = ALPHA; \
   const double              * StX; \
   register int              i; \
   int                       nu; \
   const int                 incX2 = 2 * INCX, incX3 = 3 * INCX, \
                             incX4 = 4 * INCX, incX5 = 5 * INCX, \
                             incX6 = 6 * INCX, incX7 = 7 * INCX, \
                             incX8 = 8 * INCX; \
 \
   if( ( N > 0 ) && ( alpha != HPL_rone ) ) \
   { \
      if( alpha == HPL_rzero ) \
      { \
         if( ( nu = ( N >> 3 ) << 3 ) != 0 ) \
         { \
            StX = (double *)(X) + nu * INCX; \
 \
            do \
            { \
               (*X)       = HPL_rzero; *(X+incX4) = HPL_rzero; \
               *(X+INCX)  = HPL_rzero; *(X+incX5) = HPL_rzero; \
               *(X+incX2) = HPL_rzero; *(X+incX6) = HPL_rzero; \
               *(X+incX3) = HPL_rzero; *(X+incX7) = HPL_rzero; X += incX8; \
 \
            } while( X != StX ); \
         } \
 \
         for( i = N - nu; i != 0; i-- ) { *X = HPL_rzero; X += INCX; } \
      } \
      else \
      { \
         if( ( nu = ( N >> 3 ) << 3 ) != 0 ) \
         { \
            StX = X + nu * INCX; \
 \
            do \
            { \
               x0 = (*X);       x4 = *(X+incX4); x1 = *(X+INCX);  x5 = *(X+incX5); \
               x2 = *(X+incX2); x6 = *(X+incX6); x3 = *(X+incX3); x7 = *(X+incX7); \
 \
               x0 *= alpha;   x4 *= alpha;   x1 *= alpha;   x5 *= alpha; \
               x2 *= alpha;   x6 *= alpha;   x3 *= alpha;   x7 *= alpha; \
 \
               (*X)       = x0; *(X+incX4) = x4; *(X+INCX)  = x1; *(X+incX5) = x5; \
               *(X+incX2) = x2; *(X+incX6) = x6; *(X+incX3) = x3; *(X+incX7) = x7; \
 \
               X  += incX8; \
 \
            } while( X != StX ); \
         } \
 \
         for( i = N - nu; i != 0; i-- ) \
         { x0 = (*X); x0 *= alpha; *X = x0; X += INCX; } \
      } \
   } \
}

static void HPL_accdgemmNN
(
   const int                  M,
   const int                  N,
   const int                  K,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   const double               * B,
   const int                  LDB,
   const double               BETA,
   double                     * C,
   const int                  LDC
)
{
   register double            t0;
   int                        i, iail, iblj, icij, j, jal, jbj, jcj, l;
   double                     * Cjcj;

   #pragma acc kernels copyin(A[0:N],B[0:N]) copyout(C[0:N])
   for( j = 0, jbj = 0, jcj  = 0; j < N; j++, jbj += LDB, jcj += LDC )
   {
      Cjcj = C+jcj;
      DSCAL( M, BETA, Cjcj, 1 );
//#pragma acc loop
      for( l = 0, jal = 0, iblj = jbj; l < K; l++, jal += LDA, iblj += 1 )
      {
         t0 = ALPHA * B[iblj];
	iail=jal; icij=jcj;
#pragma acc loop
         for( i = 0; i < M; i++ )
         { C[icij] += A[iail] * t0; 
	  iail += 1;  icij += 1;
	 }
      }

   }
#pragma acc end kernels
}

static void HPL_accdgemmNT
(
   const int                  M,
   const int                  N,
   const int                  K,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   const double               * B,
   const int                  LDB,
   const double               BETA,
   double                     * C,
   const int                  LDC
)
{
   register double            t0;
   int                        i, iail, ibj, ibjl, icij, j, jal, jcj, l;

   for( j = 0, ibj  = 0, jcj  = 0; j < N; j++, ibj += 1, jcj += LDC )
   {
      HPL_accdscal( M, BETA, C+jcj, 1 );
      for( l = 0, jal = 0, ibjl = ibj; l < K; l++, jal += LDA, ibjl += LDB )
      {
         t0 = ALPHA * B[ibjl];
         for( i = 0, iail = jal, icij = jcj; i < M; i++, iail += 1, icij += 1 )
         { C[icij] += A[iail] * t0; }
      }
   }
}

static void HPL_accdgemmTN
(
   const int                  M,
   const int                  N,
   const int                  K,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   const double               * B,
   const int                  LDB,
   const double               BETA,
   double                     * C,
   const int                  LDC
)
{
   register double            t0;
   int                        i, iai, iail, iblj, icij, j, jbj, jcj, l;

   for( j = 0, jbj = 0, jcj = 0; j < N; j++, jbj += LDB, jcj += LDC )
   {
      for( i = 0, icij = jcj, iai = 0; i < M; i++, icij += 1, iai += LDA )
      {
         t0 = HPL_rzero;
         for( l = 0, iail = iai, iblj = jbj; l < K; l++, iail += 1, iblj += 1 )
         { t0 += A[iail] * B[iblj]; }
         if( BETA == HPL_rzero ) C[icij]  = HPL_rzero;
         else                    C[icij] *= BETA;
         C[icij] += ALPHA * t0;
      }
   }
}

static void HPL_accdgemmTT
(
   const int                  M,
   const int                  N,
   const int                  K,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   const double               * B,
   const int                  LDB,
   const double               BETA,
   double                     * C,
   const int                  LDC
)
{
   register double            t0;
   int                        i, iali, ibj, ibjl, icij, j, jai, jcj, l;

   for( j = 0, ibj = 0, jcj  = 0; j < N; j++, ibj += 1, jcj += LDC )
   {
      for( i = 0, icij = jcj, jai = 0; i < M; i++, icij += 1, jai += LDA )
      {
         t0 = HPL_rzero;
         for( l = 0,      iali  = jai, ibjl  = ibj;
              l < K; l++, iali += 1,   ibjl += LDB ) t0 += A[iali] * B[ibjl];
         if( BETA == HPL_rzero ) C[icij]  = HPL_rzero;
         else                    C[icij] *= BETA;
         C[icij] += ALPHA * t0;
      }
   }
}

static void HPL_accdgemm0
(
   const enum HPL_TRANS       TRANSA,
   const enum HPL_TRANS       TRANSB,
   const int                  M,
   const int                  N,
   const int                  K,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   const double               * B,
   const int                  LDB,
   const double               BETA,
   double                     * C,
   const int                  LDC
)
{
   int                        i, j;

   if( ( M == 0 ) || ( N == 0 ) ||
       ( ( ( ALPHA == HPL_rzero ) || ( K == 0 ) ) &&
         ( BETA == HPL_rone ) ) ) return;

   if( ALPHA == HPL_rzero )
   {
      for( j = 0; j < N; j++ )
      {  for( i = 0; i < M; i++ ) *(C+i+j*LDC) = HPL_rzero; }
      return;
   }

   if( TRANSB == HplNoTrans )
   {
      if( TRANSA == HplNoTrans )
      { HPL_accdgemmNN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
      else
      { HPL_accdgemmTN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
   }
   else
   {
      if( TRANSA == HplNoTrans )
      { HPL_accdgemmNT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
      else
      { HPL_accdgemmTT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
   }
}


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
)
{
   if( ORDER == HplColumnMajor )
   {
      HPL_accdtrsm0( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB );
   }
   else
   {
      HPL_accdtrsm0( ( SIDE == HplRight ? HplLeft  : HplRight ),
                  ( UPLO == HplLower ? HplUpper : HplLower ),
                  TRANS, DIAG, N, M, ALPHA, A, LDA, B, LDB );
   }
}

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
)
{
    if( ORDER == HplColumnMajor )
    {
       HPL_accdgemm0( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA,
                   C, LDC );
    }
    else
    {
       HPL_accdgemm0( TRANSB, TRANSA, N, M, K, ALPHA, B, LDB, A, LDA, BETA,
                   C, LDC );
    }
}
