#include <R.h>
#include "math.h"

#define epsilon 1e-5
#define iszero(x) (((x)<epsilon) && ((-x)<epsilon))
#define max(x,y)  ((x)<(y))? (y):(x)
// dist = kendalltaudist(vect);
// dist(i,j)=kendall-tau-dist(vect[,i],vect[,j])
//
extern "C"
{
void kendalltaudistFromTemp( double *vect, int *np,int *mp, int *temp, double *dist)
{
    //I/O variables
    //n number of elements in vect and m the number of vectors
    //vect to be processed
    // comp    should be of size dimN^2*dimM
    int n = *np;
    int m = *mp;
    
    //Internal variables
    int i,j,k;//indices
    double d_ijk;


    for(k=0;k<m;k++)
      {
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
            {  
                  if(i == j)
                    continue;
              
                  d_ijk = vect[ i + k * n ] - vect[ j + k * n ];
                              
                  if( iszero(d_ijk))
                    dist[k]+= 0.5;
                  else if((d_ijk>0 && temp[i + n *j]>0) || (d_ijk<0 && temp[i + n *j]<0))
                      dist[k]++;                                          
                            
            }   
      }

}
}
