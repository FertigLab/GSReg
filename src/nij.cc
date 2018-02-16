#include <R.h>
#include "math.h"

#define max(x,y)  ((x)<(y))? (y):(x)
// 
// [N,M] = nij(vect)
//
// vect is a matrix. The columns are samples and the rows are genes.
// N[i,j]= #(vect[,i]<vect[,j])
// M[i,j]= #(vect[,i]>vect[,j])

 extern "C"
{
void Nij( double *vect, int *np,int *mp, double *N)
{
   //I/O variables
    //n number of elements in vect and m the number of vectors
    //vect to be processed
    // comp    should be of size dimN^2*dimM
    int n = *np;
    int m = *mp;
    //Internal variables
    int i,j,k;//indices
    
    for(k=0;k<m;k++)
    {
        for(i=0;i<n;i++)
            for(j=i+1;j<n;j++)
            {  
                if( vect[i+k*n] < vect[j+k*n])
                {
                    N[i+j*n]++;                  
                }else if( vect[i+k*n] > vect[j+k*n]){
                    N[j+i*n]++;
                }
            }   
    }
}
}
