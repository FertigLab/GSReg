#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix Nij(const Rcpp::NumericMatrix &V)
{
    Rcpp::NumericMatrix N(V.nrow(), V.nrow());
    for (unsigned k = 0; k < V.ncol(); ++k)
    {
        for (unsigned i = 0; i < V.nrow(); ++i)
        {
            for (unsigned j = i + 1; j < V.nrow(); ++j)
            {
                if (V(i,k) < V(j,k))
                {
                    N(i,j) += 1.0;
                }
                else if (V(i,k) > V(j,k))
                {
                    N(j,i) += 1.0;
                }
            }
        }
    }
    return N;
}