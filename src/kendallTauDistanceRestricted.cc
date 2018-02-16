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
void kendalltaudistRestricted( double *vect, int *np,int *mp, int *restriction, double *dist)
{
    //I/O variables
    //n number of elements in vect and m the number of vectors
    //vect to be processed
    // comp    should be of size dimN^2*dimM
    int n = *np;
    int m = *mp;
    bool iseqk,iseql;
    //Internal variables
    int i,j,k,l;//indices
	 
	double sum;
    double d_ijk, d_ijl;

    for( k = 0 ; k < m; k++)
		for( l = k+1; l < m; l++)
		{
			
			for( i = 0 ; i < n ; i++)
			{
				for( sum = 0, j=0; j<n ; j++)
				{  
			

					if( restriction[ i + j * n ] > 0 )
					{	
			
						d_ijk = vect[ i + k * n ] - vect[ j + k * n ];
						d_ijl = vect[ i + l * n ] - vect[ j + l * n ];
            

            
						iseqk = iszero(d_ijk);
						iseql = iszero(d_ijl);
          
						if(!iseqk && !iseql)
						{
							if((d_ijk>0 && d_ijl<0) || (d_ijk<0 && d_ijl>0))
							sum ++;
							/*else continue;*/
						}else if( (iseqk && !iseql) || (!iseqk && iseql))
							sum += 0.5;
						/*else if(iseqk && iseql)
						continue;*/
            			                    
											
					}
          
				}
			
				dist[ k + l * m ] += sum ;	
				dist[ l + k * m ] += sum ;					
			
			
			}

		}
}
}
