#include <R.h>
#include "math.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>
using namespace std;

#define epsilon 1e-5
#define iszero(x) (((x)<epsilon) && ((-x)<epsilon))
#define max(x,y)  ((x)<(y))? (y):(x)
// dist = kendalltaudist(vect);
// dist(i,j)=kendall-tau-dist(vect[,i],vect[,j])
//
extern "C"
{
void kendalltaudistWrap( double *vect, int *np,int *mp, double *incVec, double *dist)
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
    double d_ijk, d_ijl;
    double x[n];
    double y[n];
    double xc[n];
    double yc[n];
  	int tempTauDist;
    for(k=0;k<m;k++){
      for (i=0;i<n;i++){
      	x[i]=vect[i+k*n];
      	xc[i]=incVec[i+k*n];
      }
      for(l=k+1;l<m;l++){
      	for(i=0;i<n;i++){
      		y[i]=vect[i+k*n];
      		yc[i]=incVec[i+l*n];
      	}
      	tempTauDist=ktau_p(x,xc,y,yc);      //computes cenken for each pair of vectors
      	dist[k+l*m]=tempTauDist;
      	dist[l+k*m]=tempTauDist;
        /*for(i=0;i<n;i++){      //Old censoring of kendall tau math
          for(j=i+1;j<n;j++){  
            if((incVec[i+k*n]==0) || (incVec[j+k*n]==0)){
              d_ijk=0;
            }else{
              d_ijk = vect[ i + k * n ] - vect[ j + k * n ];
            }

            if((incVec[i+l*n]==0) || (incVec[j+l*n]==0)){
              d_ijl=0;
            }else{
              d_ijl = vect[ i + l * n ] - vect[ j + l * n ];
            }
            
            iseqk = iszero(d_ijk);
            iseql = iszero(d_ijl);
          
            if(!iseqk && !iseql){
              if((d_ijk>0 && d_ijl<0) || (d_ijk<0 && d_ijl>0)){
                dist[k+l*m]++;
                dist[l+k*m]++;
              }
                      
            }else if( (iseqk && !iseql) || (!iseqk && iseql)){
              dist[k+l*m]+= 0.5;
              dist[l+k*m]+= 0.5;
            }                  
                          
          }
        }   */
      }
    }
}
}
double ktau_p(double x, double xcen, double y, double ycen){
	double tau; //We only need the Tau output, ignore the p-value output
	double tau_b;
	double ymin=Minimum(y);
	for(int i=0;i<Length(*y);i++){
		y[i]=y[i]-ymin/1000*ycen[i];
	}
	vector<double> xx (x,x+Length(x));
	vector<double> cx (xcen,xcen+Length(xcen));
	vector<double> yy (y,y+Length(y));
	vector<double> cy (ycen,ycen+Length(ycen));
	double n=xx.size();
	vector<double> sx=xx;
	vector<double> sy=xy;
	sort(sx.begin(),sx.end());
	sort(sy.begin(),sy.end());
	sx.erase(unique(sx.begin(),sx.end()),sx.end());
	sy.erase(unique(sy.begin(),sy.end()),sy.end());
	double delx=min_element(vector_diff(sx))/1000;
	double dely=min_element(vector_diff(sy))/1000;
	vector<double> dupx=vector_sub(xx,vector_scaMul(cx,delx));
	vector<double> diffx=outer_Sub(dupx,dupx);
	vector<double> diffcx=outer_Sub(cx,cx);
	vector<double> xplus=outer_sub(cx,vector_scaMul(cx,-1));
	vector<double> dupy=vector_sub(yy,vector_scaMul(cy,dely));
	vector<double> diffy=outer_Sub(dupy,dupy);
	vector<double> diffcy=outer_Sub(cy,cy);
	vector<double> yplus=outer_sub(cy,vector_scaMul(cy,-1));
	vector<int> signyx=sign(vector_Prod(diffy,diffx));
	double tt=(vector_sum(vector_scaAdd(vector_scaMul(vector_abs(sign(diffx)),-1),1))-n)/2;
	double uu=(vector_sum(vector_scaAdd(vector_scaMul(vector_abs(sign(diffy)),-1),1))-n)/2;
	vector<double> cix=vector_Prod(sign(diffcx),diffx);
	for(int i=0;i<cix.size();i++){
		if(cix[i]<=0){cix[i]=0;} 
		else{cix[i]=1;} 
	}
	tt=tt+vector_sum(cix)/2;
	signyx=vector_Prod(signyx,vector_scaAdd(vector_scaMul(cix,-1),1));
	vector<double> ciy=vector_Prod(sign(diffcy),diffy);
	for(int i=0;i<ciy.size();i++){
		if(ciy[i]<=0){ciy[i]=0;} 
		else{ciy[i]=1;} 
	}
	uu=uu+vector_sum(ciy)/2;
	signyx=vector_Prod(signyx,vector_scaAdd(vector_scaMul(ciy,-1),1));
	for(int i=0;i<xplus.size();i++){
		if(xplus[i]<=0){xplus[i]=0;} 
		else{xplus[i]=1;} 
	}
	for(int i=0;i<yplus.size();i++){
		if(yplus[i]<=0){yplus[i]=0;} 
		else{yplus[i]=1;} 
	}
	diffx=vector_abs(sign(diffx));
	diffy=vector_abs(sign(diffy));
	vector<double> tplus=vector_Prod(xplus,diffx);
	vector<double> uplus=vector_Prod(yplus,diffy);
	tt=tt+vector_sum(tplus)/2;
	uu=uu+vector_sum(uplus)/2;
	double itot=vector_sum(vector_Prod(signyx,vector_Prod(vector_scaAdd(vector_scaMul(xplus,-1),1),vector_scaAdd(vector_scaMul(yplus,-1),1))));
	double kenS=itot/2;
	tau=tau/(n*(n-1));
	double J=n*(n-1)/2;
	tau_b=kenS/sqrt((J-tt)*(J-uu));
	return tau;
}
vector<double> vector_abs(vector<double> vec){
	vector<double> out;
	for(int i=0;i<vec.size();i++){
		out+=abs(vec[i]);
	}
	return out;
}
double vector_sum(vector<double> vec){
	double out=0;
	for(int i=0;,i<vec.size();i++){
		out=out+vec[i];
	}
	return out;
}
vector<double> vector_scaAdd(vector<double> vec, double scal){
	vector<double>out;
	for(int i=0;i<vec.size();i++){
		out+=vec[i]+scal;
	}
	return out;
}
vector<double> vector_diff(vector<double> vec){
	vector<double> out;
	for(int i=0;i<(vec.size()-1);i++){
		out+=vec[i+1]-vec[i];
	}
	return out;
}
vector<double> vector_sub(vector<double> vec1,vector<double> vect2){
	vector<double> out;
	for(int i=0;i<vec1.size();i++){
		out+=vec1[i]-vec2[i];
	}
	return out;
}
vector<double> vector_scaMul(vector<double> vec, double scal){
	vector<double> out;
	for(int i=0;i<vec.size();i++){
		out+=vec[i]*scal;
	}
	return out;
}
vector<double> vector_Prod(vector<double> vec1, vector<double> vec1){
	vector<double> out;
	for(int i=0;i<vec1.size();i++){
		out+=vec1[i]*vec2[i];
	}
	return out;
}
vector<double> outer_Sub(vector<double> vec1, vector<double> vec2){
	vector<double> out;
	for(int i=0;i<Length(vec1);i++){
		for(int j=0;j<Length(vec2);j++){
			out+=vec1[i]-vec2[j];
		}
	}
	return out;
}
vector<int> sign(vector<double> vec){
	vector<int> out;
	for(int i=0;i<vec.size();i++){
		if(vec[i]<0){out+=-1;}
		else if(vec[i]>0){out+=1;}
		else {out+=0;}
	}
	return out;
}
double Diff(double vec){
	tArrLen=Length(vec);
	double out[tArrLen-1];
	for(int i=0;i<(tArrLen-1);i++){
		out[i]=vec[i+1]-vec[i];
	}
	return out;
}
double Minimum(double vec){
	double out=10^7;
	for(int i=0;i<Length(vec);i++){
		if(out>vec[i]){
			out=vec[i];
		}
	}
	return out;
}
int Length(double vec){
	return (sizeof(vec)/sizeof(vec[0]));
}
/*
void kendallATS(double *y, double *ycen, double *x, double *xcen, double tol=10^(-7), double iter= 10^3){
	double ymin=Minimum(*y);
	for(int i=0;i<Length(*y);i++){
		y[i]=y[i]-ymin/1000*ycen[i];
	}
	double k[]=ktau_b(*x,*xcen,*y,*ycen);
	double bs[]=iter_s(k[0],k[1],*x,*xcen,*y,*ycen,tol=tol,iter=iter);
	double ubs[]=iter_s(bs[1],bs[2],*x,*xcen,*y,*ycen,tol=tol,iter=iter);
	double lbs[]=iter_s(bs[0],bs[1],*x,*xcen,*y,*ycen,tol=tol,iter=iter);
	double slope=0.5*(ubs[1]+lbs[1]);
	double inter=turnbull(*y,*ycen,*x,*xcen,slope,tol);
	double p_tau[]=ktau_p(*x,*xcen,*y,*ycen);
	double out[4];
	out[0]=slope;
	out[1]=inter;
	out[2]=p_tau[0];
	out[3]=p_tau[1];
}
double step_s(double *lb, double *ub, double *x, double *xcen, double *y, double *ycen){
	double out[6];
	out[0]=*lb;
	out[1]=(*lb + *ub)/2;
	out[2]=*ub;
	int tArrLen=Length(*y);
	double res[tArrLen*3];
	for(int i=0;i<tArrLen*3;i++){
		res[i]=y[i%tArrLen]-out[i/tArrLen]*x[i%tArrLen];
	}
	for(int i=0;i<3;i++){
		double resTemp[tArrLen];
		for(int j=0;j<tArrLen;j++){
			resTemp[tArrLen]=res[i*tArrLen+j];
		}
		out[i+3]=ktau_s(*x,*xcen,resTemp,*ycen);
	}
	return out;
}
double iter_s(double *lb, double *ub, double *x, double *xcen, double *y, double *ycen, double tol, double iter){
	double bs[]=step_s(*lb,*ub,*x,*xcen,*y,*ycen);
	for(i=0;i<iter;i++){
		if((bs[3]*bs[4])<=0){
			bs[5]=bs[4];
			bs[2]=bs[1];
		} else{
			bs[3]=bs[4];
			bs[0]=bs[1];
		}
		if((bs[4]==0)||(abs(bs[2]-out[0])<=tol)){
			break;
		}
		bs=step_s(bs[0],bs[2],*x,*xcen,*y,*ycen);
	}
	return bs;
}
double ktau_s(double *x, double *xcen, double *y, double *ycen){

}
double ktau_b(double *x, double *xcen, double *y, double *ycen){

}
double turnbull(double *y, double *ycen, double *x, double *xcen, double slope, double tol){

}
*/

