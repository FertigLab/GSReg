#include <Rcpp.h>
#include "math.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>
using namespace std;

#define epsilon 1e-5
#define iszero(x) (((x)<epsilon) && ((-x)<epsilon))

int Length(double *vec){
    return (sizeof(vec)/sizeof(vec[0]));
}
double Minimum(double *vec){
    double out=10^7;
    for(int i=0;i<Length(vec);i++){
        if(out>vec[i]){
            out=vec[i];
        }
    }
    return out;
}
double vector_min(vector<double> vec){
    double out=10^7;
    for (int i=0;i<vec.size();i++){
        if(out>vec[i]){
            out=vec[i];
        }
    }
    return out;
}
vector<double> vector_abs(vector<double> vec){
    vector<double> out(vec.size());
    for(int i=0;i<vec.size();i++){
        out[i]=abs(vec[i]);
    }
    return out;
}
double vector_sum(vector<double> vec){
    double out=0;
    for(int i=0;i<vec.size();i++){
        out=out+vec[i];
    }
    return out;
}
vector<double> vector_scaAdd(vector<double> vec, double scal){
    vector<double> out(vec.size());
    for(int i=0;i<vec.size();i++){
        out[i]=vec[i]+scal;
    }
    return out;
}
vector<double> vector_diff(vector<double> vec){
    vector<double> out(vec.size());
    for(int i=0;i<(vec.size()-1);i++){
        out[i]=vec[i+1]-vec[i];
    }
    return out;
}
vector<double> vector_sub(vector<double> vec1,vector<double> vec2){
    vector<double> out(vec1.size());
    for(int i=0;i<vec1.size();i++){
        out[i]=vec1[i]-vec2[i];
    }
    return out;
}
vector<double> vector_scaMul(vector<double> vec, double scal){
    vector<double> out(vec.size());
    for(int i=0;i<vec.size();i++){
        out[i]=vec[i]*scal;
    }
    return out;
}
vector<double> vector_Prod(vector<double> vec1, vector<double> vec2){
    vector<double> out(vec1.size());
    for(int i=0;i<vec1.size();i++){
        out[i]=(vec1[i]*vec2[i]);
    }
    return out;
}
vector<double> outer_Sub(vector<double> vec1, vector<double> vec2){
    vector<double> out(vec1.size());
    for(int i=0;i<vec1.size();i++){
        for(int j=0;j<vec2.size();j++){
            out[i]=(vec1[i]-vec2[j]);
        }
    }
    return out;
}
vector<double> sign(vector<double> vec){
    vector<double> out(vec.size());
    for(int i=0;i<vec.size();i++){
        if(vec[i]<0){out[i]=(-1);}
        else if(vec[i]>0){out[i]=(1);}
        else {out[i]=(0);}
    }
    return out;
}

double ktau_p(vector<double> x, vector<double> xcen, vector<double> y, vector<double> ycen){
    double tau; //We only need the Tau output, ignore the p-value output
    double ymin=vector_min(y);
    for(int i=0;i<y.size();i++){
        y[i]=y[i]-ymin/1000*ycen[i];
    }
    vector<double> xx =x;
    vector<double> cx =xcen;
    vector<double> yy =y;
    vector<double> cy =ycen;
    double n=xx.size();
    vector<double> sx=xx;
    vector<double> sy=yy;
    sort(sx.begin(),sx.end());
    sort(sy.begin(),sy.end());
    sx.erase(unique(sx.begin(),sx.end()),sx.end());
    sy.erase(unique(sy.begin(),sy.end()),sy.end());
    double delx=vector_min(vector_diff(sx))/1000;
    double dely=vector_min(vector_diff(sy))/1000;
    vector<double> dupx=vector_sub(xx,vector_scaMul(cx,delx));
    vector<double> diffx=outer_Sub(dupx,dupx);
    vector<double> diffcx=outer_Sub(cx,cx);
    vector<double> xplus=outer_Sub(cx,vector_scaMul(cx,-1));
    vector<double> dupy=vector_sub(yy,vector_scaMul(cy,dely));
    vector<double> diffy=outer_Sub(dupy,dupy);
    vector<double> diffcy=outer_Sub(cy,cy);
    vector<double> yplus=outer_Sub(cy,vector_scaMul(cy,-1));
    vector<double> signyx=sign(vector_Prod(diffy,diffx));
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
    //double kenS=itot/2;
    tau=itot/(n*(n-1));
    double out=0.5-tau/2;
    //double J=n*(n-1)/2;
    //tau_b=kenS/sqrt((J-tt)*(J-uu));
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector kendalltaudistWrap(const Rcpp::NumericVector vect,const  int np, int mp,const Rcpp::NumericVector incVec)
{
    Rcpp::NumericVector dist(mp*mp);
    //I/O variables
    //n number of elements in vect and m the number of vectors
    //vect to be processed
    // comp    should be of size dimN^2*dimM
    int n = np;
    int m = mp;
    //Internal variables
    int i,k,l;//indices
    vector<double> x(n);
    vector<double> y(n);
    vector<double> xc(n);
    vector<double> yc(n);
    double tempTauDist;
    for(k=0;k<m;k++){
      for (i=0;i<n;i++){
        x[i]=(vect[i+k*n]);
        xc[i]=(incVec[i+k*n]);
      }
      for(l=k+1;l<m;l++){
        for(i=0;i<n;i++){
            y[i]=(vect[i+k*n]);
            yc[i]=(incVec[i+l*n]);
        }
        tempTauDist=ktau_p(x,xc,y,yc);      //computes cenken for each pair of vectors
        dist[k+l*m]=tempTauDist;
        dist[l+k*m]=tempTauDist;
      }
    }
    return dist;
}