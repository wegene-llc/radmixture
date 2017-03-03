// function for updating F under supervised learning

# include <R.h>
# include <Rinternals.h>
# include <Rmath.h>
# include <stdlib.h>
# include <matrix.h>


SEXP updatef (SEXP n1, SEXP p1, SEXP K1, SEXP Rq0, SEXP Rf0, SEXP Rg) {
    /*
        Update F matrix.
    */
    int n, p, K;
    int i, j, k, m;
    double a, b, temp1, temp2;
    double *ff0;
    Matrix *q0, *f0, *g0, *f, *mtemp1, *mtemp2;
    SEXP ff;
    n = asInteger(n1);
    p = asInteger(p1);
    K = asInteger(K1);
    Rg = coerceVector(Rg, REALSXP);
    Rq0 = coerceVector(Rq0, REALSXP);
    Rf0 = coerceVector(Rf0, REALSXP);
    f0 = init_matrix(K, p, REAL(Rf0));
    q0 = init_matrix(n, K, REAL(Rq0));
    g0 = init_matrix(n, p, REAL(Rg));
    f = init_matrix(K, p, NULL);
    mtemp1 = init_matrix(n, p, NULL);
    mtemp2 = init_matrix(n, p, NULL);

    for (j = 0; j < p; j++)
    {
        for (i = 0; i < n; i++)
        {
            temp1 = temp2 = 0;
            for (m = 0; m < K; m++)
            {
                temp1 += GET(q0, i, m) * GET(f0, m, j);
                temp2 += GET(q0, i, m) * (1 - GET(f0, m, j));
            }
            SET(mtemp1, i, j, temp1);
            SET(mtemp2, i, j, temp2);
        }
    }
    for (k = 0; k < K; k++)
    {
        for (j = 0; j < p; j++)
        {
            a = b = 0;
            for (i = 0; i < n; i++)
            {
                a += (GET(g0, i, j) * GET(q0, i, k) * GET(f0, k, j)) / GET(mtemp1, i, j);
                b += ((2 - GET(g0, i, j)) * GET(q0, i, k) * (1 - GET(f0, k, j))) / GET(mtemp2, i, j);
            }
            SET(f, k, j, a / (a + b));
        }
    }
    PROTECT(ff = allocMatrix(REALSXP, K, p));
    ff0 = REAL(ff);
    for (k = 0; k < K; k++)
    {
        for (j = 0; j < p; j++)
        {
            ff0[k + K * j] = GET(f, k, j);
        }
    }
    UNPROTECT(1);
    free_matrix(f0);
    free_matrix(q0);
    free_matrix(g0);
    free_matrix(f);
    free_matrix(mtemp1);
    free_matrix(mtemp2);
    return (ff);
}
