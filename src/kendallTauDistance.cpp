#include <Rcpp.h>

#define epsilon 1e-5
#define iszero(x) (((x)<epsilon) && ((-x)<epsilon))

// [[Rcpp::export]]
Rcpp::NumericMatrix kendalltaudist(const Rcpp::NumericMatrix &V)
{
    Rcpp::NumericMatrix D(V.ncol(), V.ncol());
    for (unsigned k = 0; k < V.ncol(); ++k)
    {
        for (unsigned l = k + 1; l < V.ncol(); ++l)
        {
            for (unsigned i = 0; i < V.nrow(); ++i)
            {
                for (unsigned j = i + 1; j < V.nrow(); ++j)
                {
                    double d_ijk = V(i,k) - V(j,k);
                    double d_ijl = V(i,l) - V(j,l);
                
                    if (!iszero(d_ijk) && !iszero(d_ijl))
                    {
                        if ((d_ijk > 0.0 && d_ijl < 0.0) || (d_ijk < 0.0 && d_ijl > 0.0))
                        {
                            D(k,l) += 1.0;
                            D(l,k) += 1.0;
                        }
                    }
                    else if ((iszero(d_ijk) && !iszero(d_ijl)) || (!iszero(d_ijk) && iszero(d_ijl)))
                    {
                        D(k,l) += 0.5;
                        D(l,k) += 0.5;
                    }
                }
            }
        }
    }
    return D;
}