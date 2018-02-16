#include <Rcpp.h>

#define epsilon 1e-5
#define iszero(x) (((x)<epsilon) && ((-x)<epsilon))

// [[Rcpp::export]]
Rcpp::NumericMatrix kendalltaudistRestricted(const Rcpp::NumericMatrix &V, const Rcpp::NumericMatrix &R)
{
    Rcpp::NumericMatrix D(V.ncol(), V.ncol());
    for (unsigned k = 0; k < V.ncol(); ++k)
    {
        for (unsigned l = k + 1; l < V.ncol(); ++l)
        {
            for (unsigned i = 0; i < V.nrow(); ++i)
            {
                double sum = 0.0;
                for (unsigned j = 0; j < V.nrow(); ++j)
                {
                    if (static_cast<int>(R(i,j)) > 0)
                    {
                        double d_ijk = V(i,k) - V(j,k);
                        double d_ijl = V(i,l) - V(j,l);
                
                        if (!iszero(d_ijk) && !iszero(d_ijl))
                        {
                            if ((d_ijk > 0.0 && d_ijl < 0.0) || (d_ijk < 0.0 && d_ijl > 0.0))
                            {
                                sum += 1;
                            }
                        }
                        else if ((iszero(d_ijk) && !iszero(d_ijl)) || (!iszero(d_ijk) && iszero(d_ijl)))
                        {
                            sum += 0.5;
                        }
                    }
                }
                D(k,l) += sum;
                D(l,k) += sum;
            }
        }
    }
    return D;
}
