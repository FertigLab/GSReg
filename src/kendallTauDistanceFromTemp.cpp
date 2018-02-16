#include <Rcpp.h>

#define epsilon 1e-5
#define iszero(x) (((x)<epsilon) && ((-x)<epsilon))

// [[Rcpp::export]]
Rcpp::NumericVector kendalltaudistFromTemp(const Rcpp::NumericMatrix &V, const Rcpp::NumericMatrix &T)
{
    Rcpp::NumericVector D(V.ncol());
    for (unsigned k = 0; k < V.ncol(); ++k)
    {
        for (unsigned i = 0; i < V.nrow(); ++i)
        {
            for (unsigned j = 0; j < V.nrow(); ++j)
            {
                if (i != j)
                {
                    double d_ijk = V(i,k) - V(j,k);
                    int temp = static_cast<int>(T(i,j));
                    if (iszero(d_ijk))
                    {
                        D(k) += 0.5;
                    }
                    else if ((d_ijk > 0.0 && temp > 0) || (d_ijk < 0.0 && temp < 0))
                    {
                        D(k) += 1.0;
                    }
                }
            }
        }
    }
    return D;
}