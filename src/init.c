/* Registration of C routines */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void Nij( double *vect, int *np,int *mp, double *N);
void kendalltaudist( double *vect, int *np,int *mp, double *dist);
void kendalltaudistRestricted( double *vect, int *np,int *mp, int *restriction, double *dist);
void kendalltaudistFromTemp( double *vect, int *np,int *mp, int *temp, double *dist);
void vect2compC( double *vect, int *np,int *mp, double *comp);
void kendalltaudistWrap( double *vect, int *np, int *mp, double *incVec, double *dist);



static const R_CMethodDef R_CDef[] = {
   {"Nij", (DL_FUNC)&Nij,4},
  {"kendalltaudist", (DL_FUNC)&kendalltaudist,4},
  {"kendalltaudistWrap", (DL_FUNC)&kendalltaudistWrap,5},
  {"kendalltaudistRestricted", (DL_FUNC)&kendalltaudistRestricted,5},
  {"vect2compC", (DL_FUNC)&vect2compC,4},
  {"kendalltaudistFromTemp", (DL_FUNC)&kendalltaudistFromTemp,5},
  {NULL, NULL, 0},
};

void R_init_GSReg(DllInfo *info)
{
  R_registerRoutines(info,R_CDef,NULL,NULL,NULL);
}
