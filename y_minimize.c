#include "y_minimize.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>


/******************************** G R A D I E N T A R     M E T H O D S ******************************************/


/********************************* L I N E    S E A R C H *********************************************/


//Inverse quadratic interpolation of extreme point (http://fedc.wiwi.hu-berlin.de/xplore/tutorials/xegbohtmlnode62.html)
inline double _inverse_quadratic_interpolation(register double a,register double b,register double c,register double fa,register double fb,register double fc)
{
return b-.5*( (b-a)*(b-a)*(fb-fc)-(b-c)*(b-c)*(fb-fa) )/( (b-a)*(fb-fc)-(b-c)*(fb-fa) );
}

//Line-search method with square_approximation but without curvatures
//NOTE. *lambda is the lenght of the previous step, *ff is initial function value
unsigned int line_search_square_fapproximation(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int n,double *x0,double *x1,double *g1,double **G,double *p,unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list stack)
{
double a[4], f[4], max_p, pg;
unsigned int _i;

//Stage 0. Find the biggest p component
max_p=fabs(*p), _i=n; while (_i--) if (max_p<fabs(p[_i])) max_p=fabs(p[_i]);
//Stage I. Get upper limit
f[0]=*ff, a[0]=0., a[2]=*lambda;
if (a[2]*max_p>max_variable_change) a[2]=0.9*max_variable_change/max_p;  
vect_mult_summ_vect(n,x1,a[2],p,x0);
f[2]=0.; //Init to avoid occasionaly NAN which is a flag for optimization with exclusions
if ((_i=funct(&f[2],n,x1,g1,G,stack))!=YERROR_OK) return _i;
if (f[2]<f[0])
  {
  if (calc_vect_vect_scalar_product(n,p,g1)>0.)
    { 
    STEP_IN_DEEP: f[1]=f[2], a[1]=a[2];
    if (max_p*a[2]*GOLDEN_RATIO>max_variable_change) { *lambda=a[2]; return YERROR_OK; } //maximal extension reached  
    else 
      {
      a[2]*=GOLDEN_RATIO; //Step deeper
      vect_mult_summ_vect(n,x1,a[2],p,x0);
      if ((_i=funct(&f[2],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
      if ( (f[2]<f[1])&&(calc_vect_vect_scalar_product(n,p,g1)) ) { f[0]=f[1], a[0]=a[1]; goto STEP_IN_DEEP; }
      //Else apply Brant method to find the (exact) minima
      }
    }
  else
    {//Do left-hand truncation |---|--|
    STEP_IN_MIDDLE: a[1]=a[2]-GOLDEN_SECTION*(a[2]-a[0]);
    vect_mult_summ_vect(n,x1,a[2],p,x0);
    if ((_i=funct(&f[1],n,x1,g1,G,stack))!=YERROR_OK) return _i;
    if ( (f[2]<f[1])||(f[0]<f[1]) )
      {
      if (calc_vect_vect_scalar_product(n,p,g1)>0.) { f[0]=f[1], a[0]=a[1]; } else { f[2]=f[1], a[2]=a[1]; }
      goto STEP_IN_MIDDLE;
      }     
    //Else apply Brant method to find the (exact) minima
    }
  }
else
  {
  STEP_ON_SHALLOW: a[1]=a[2]*GOLDEN_SECTION; //Step upper
  vect_mult_summ_vect(n,x1,a[1],p,x0);
  if ((_i=funct(&f[1],n,x1,g1,G,stack))!=YERROR_OK) return _i;
  if (f[1]>=f[0]) 
    {
    f[2]=f[1], a[2]=a[1]; 
    if ((a[2]-a[0])>itol/max_p) goto STEP_ON_SHALLOW;
    else return YERROR_NCONVERGED; 
    }
  //Else apply Brant method to find the (exact) minima 
  }
//Stage II. Quadratic interpolation: the points are A-->B<-C 
while ( ((a[2]-a[0])*max_p>itol)&&( ((f[2]-f[1])>itol)||((f[0]-f[1])>itol) )&&(calc_vect_norm(n,g1)>itol) )
  {//do quadratic interpolation serch
  a[3]=.5*(a[0]+a[1]-(f[1]-f[0])*(a[2]-a[1])/(a[1]-a[0])/((f[2]-f[0])/(a[2]-a[0])-(f[1]-f[0])/(a[1]-a[0])));
       if ( (pg>0.)&&(a[3]>a[1]) ) { if ((a[1]-a[0])*max_p>.5*itol) a[3]=(a[1]+a[0])/2.; else break; } //Contradiction between derivative and a values
  else if ( (pg<0.)&&(a[3]<a[1]) ) { if ((a[2]-a[1])*max_p>.5*itol) a[3]=(a[2]+a[1])/2.; else break; } //Contradiction between derivative and a values
  else   { 
              if (a[3]<=a[0]+(a[1]-a[0])/2.) a[3]=a[0]+(a[1]-a[0])/2.;
         else if (a[3]>=a[2]-(a[2]-a[1])/2.) a[3]=a[2]-(a[2]-a[1])/2.;
         }
  *lambda=a[3], vect_mult_summ_vect(n,x1,a[3],p,x0);
  if ((_i=funct(&f[3],n,x1,g1,G,stack))!=YERROR_OK) return _i;
  if (a[3]>a[1])
    {
    if (f[3]>f[1])             { a[2]=a[3], f[2]=f[3]; }
    else { a[0]=a[1], f[0]=f[1], a[1]=a[3], f[1]=f[3], pg=calc_vect_vect_scalar_product(n,g1,p); }
    }
  else
    {
    if (f[3]>f[1])             { a[0]=a[3], f[0]=f[3]; }
    else { a[2]=a[1], f[2]=f[1], a[1]=a[3], f[1]=f[3], pg=calc_vect_vect_scalar_product(n,g1,p); }
    }
  }
//Stage III. Solution found, exit
if (*lambda!=a[1])
  {//Recalculate previous min if necessary
  *lambda=a[1], vect_mult_summ_vect(n,x1,a[1],p,x0);
  return funct(ff,n,x1,g1,G,stack);
  }
else { *ff=f[1]; return YERROR_OK; }
}

//The same as previous but the search is done with functions evaluations and only the very last call (at return statements) re-calculates the gradient
//NB! The trigger of gradient calculation is to have pointers g!=0x0 and optionally G!=0
unsigned int line_search_square_fapproximation_only(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int n,double *x0,double *x1,double *g1,double **G,double *p,unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list stack)
{
double a[4], f[4], max_p;
unsigned int _i;

//Stage 0. Find the biggest p component
max_p=fabs(*p), _i=n; while (_i--) if (max_p<fabs(p[_i])) max_p=fabs(p[_i]);
//Stage I. Define 3-points problem
f[0]=*ff, a[0]=0., a[2]=*lambda;
if (a[2]*max_p>max_variable_change) a[2]=0.9*max_variable_change/max_p;  
vect_mult_summ_vect(n,x1,a[2],p,x0);
f[2]=0.; //Init to avoid occasionaly NAN which is a flag for optimization with exclusions
if ((_i=funct(&f[2],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
if (f[2]<f[0])
  {
  STEP_IN_DEEP: f[1]=f[2], a[1]=a[2];
  if (max_p*a[2]*GOLDEN_RATIO>max_variable_change) { *lambda=a[2]; return funct(ff,n,x1,g1,G,stack); } //maximal extension reached  
  else 
    {
    a[2]*=GOLDEN_RATIO; //Step deeper
    vect_mult_summ_vect(n,x1,a[2],p,x0);
    if ((_i=funct(&f[2],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
    if (f[2]<f[1]) { f[0]=f[1], a[0]=a[1]; goto STEP_IN_DEEP; }
    //Else apply Brant method to find the (exact) minima
    }
  }
else
  {
  STEP_ON_SHALLOW: a[1]=a[2]*GOLDEN_SECTION; //Step upper
  vect_mult_summ_vect(n,x1,a[1],p,x0);
  if ((_i=funct(&f[1],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
  if (f[1]>=f[0]) 
    {
    f[2]=f[1], a[2]=a[1]; 
    if ((a[2]-a[0])>itol/max_p) goto STEP_ON_SHALLOW;
    else return YERROR_NCONVERGED; 
    }
  //Else apply Brant method to find the (exact) minima 
  }
//Stage II. Quadratic interpolation: the points are A-->B<-C 
while ( ((a[2]-a[0])*max_p>itol)&&( (fabs(f[2]-f[1])>tol)||(fabs(f[0]-f[1])>tol) ) )
  {
  //do quadratic interpolation
  a[3]=.5*(a[0]+a[1]-(f[1]-f[0])*(a[2]-a[1])/(a[1]-a[0])/((f[2]-f[0])/(a[2]-a[0])-(f[1]-f[0])/(a[1]-a[0])));
       if (a[3]<=a[0]+(a[1]-a[0])/2.) a[3]=a[0]+(a[1]-a[0])/2.;
  else if (a[3]>=a[2]-(a[2]-a[1])/2.) a[3]=a[2]-(a[2]-a[1])/2.;
  vect_mult_summ_vect(n,x1,a[3],p,x0);
  if ((_i=funct(&f[3],n,x1,0X0,0X0,stack))!=YERROR_OK) return _i;
  if (a[3]>a[1])
    {
    if (f[3]>f[1])             { a[2]=a[3], f[2]=f[3]; }
    else { a[0]=a[1], f[0]=f[1], a[1]=a[3], f[1]=f[3]; }
    }
  else
    {
    if (f[3]>f[1])             { a[0]=a[3], f[0]=f[3]; }
    else { a[2]=a[1], f[2]=f[1], a[1]=a[3], f[1]=f[3]; }
    }
  }
//Stage III. Solution found, recalculate gradient at a[1]
*lambda=a[1], vect_mult_summ_vect(n,x1,a[1],p,x0);
return funct(ff,n,x1,g1,G,stack);
}

//Doing 4-point schema: generating bad-good-bad situation and then appliing one step of Brent quadratic interpolation
//NB! The trigger of gradient calculation is to have pointers g!=0x0 and optionally G!=0
unsigned int line_qsearch_square_fapproximation_only(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int n,double *x0,double *x1,double *g1,double **G,double *p,unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list stack)
{
double a[3], f[3], max_p;
unsigned int _i;

//Stage 0. Find the biggest p component
max_p=fabs(*p), _i=n; while (_i--) if (max_p<fabs(p[_i])) max_p=fabs(p[_i]);
//Stage I. Define 3-points problem
f[0]=*ff, a[0]=0., a[2]=*lambda;
if (a[2]*max_p>max_variable_change) a[2]=0.9*max_variable_change/max_p;  
vect_mult_summ_vect(n,x1,a[2],p,x0);
f[2]=0.; //Init to avoid occasionaly NAN which is a flag for optimization with exclusions
if ((_i=funct(&f[2],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
if (f[2]<f[0])
  {
  STEP_IN_DEEP: f[1]=f[2], a[1]=a[2];
  if (max_p*a[2]*GOLDEN_RATIO>max_variable_change) { *lambda=a[2]; return funct(ff,n,x1,g1,G,stack); } //maximal extension reached  
  else 
    {
    a[2]*=GOLDEN_RATIO; //Step deeper
    vect_mult_summ_vect(n,x1,a[2],p,x0);
    if ((_i=funct(&f[2],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
    if (f[2]<f[1]) { f[0]=f[1], a[0]=a[1]; goto STEP_IN_DEEP; }
    //Else apply Brant method to find the (exact) minima
    }
  }
else
  {
  STEP_ON_SHALLOW: a[1]=a[2]*GOLDEN_SECTION; //Step upper
  vect_mult_summ_vect(n,x1,a[1],p,x0);
  if ((_i=funct(&f[1],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
  if (f[1]>=f[0]) 
    {
    f[2]=f[1], a[2]=a[1]; 
    if ((a[2]-a[0])>itol/max_p) goto STEP_ON_SHALLOW;
    else return YERROR_NCONVERGED; 
    }
  //Else apply Brant method to find the (exact) minima 
  }
//Stage II. Do 1 step of quadratic interpolation
*lambda=.5*(a[0]+a[1]-(f[1]-f[0])*(a[2]-a[1])/(a[1]-a[0])/((f[2]-f[0])/(a[2]-a[0])-(f[1]-f[0])/(a[1]-a[0]))); //Can't be NAN or INF
vect_mult_summ_vect(n,x1,*lambda,p,x0);
if ((_i=funct(ff,n,x1,g1,G,stack))!=YERROR_OK) return _i;
if (*ff>f[1])
  {//Disaster! Recalculate the previous best
  *lambda=a[1], vect_mult_summ_vect(n,x1,*lambda,p,x0);
  return funct(ff,n,x1,g1,G,stack);
  }
else return YERROR_OK;
}
//The same as previous but it modifies only a subset of variables (from - to) keeping the rest frozed
unsigned int line_qssearch_square_fapproximation_only(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int from,unsigned int to,unsigned int n,double *x0,double *x1,double *g1,double **G,double *p,unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list stack)
{
double a[3], f[3], max_p;
register unsigned int _i;

//Stage 0. Find the biggest p component
max_p=fabs(p[from]), _i=to; while (--_i!=from) if (max_p<fabs(p[_i])) max_p=fabs(p[_i]);
//Stage I. Define 3-points problem
f[0]=*ff, a[0]=0., a[2]=*lambda;
if (a[2]*max_p>max_variable_change) a[2]=0.9*max_variable_change/max_p;  
{ if ( (from)) memcpy(x1,x0,sizeof(double)*from); vect_mult_summ_vect(to-from,&x1[from],a[2],&p[from],&x0[from]); if (to!=n) memcpy(&x1[to],&x0[to],sizeof(double)*(n-to)); } // Unit shift
f[2]=0.; //Init to avoid occasionaly NAN which is a flag for optimization with exclusions
if ((_i=funct(&f[2],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
if (f[2]<f[0])
  {
  STEP_IN_DEEP: f[1]=f[2], a[1]=a[2];
  if (max_p*a[2]*GOLDEN_RATIO>max_variable_change) { *lambda=a[2]; return funct(ff,n,x1,g1,G,stack); } //maximal extension reached  
  else 
    {
    a[2]*=GOLDEN_RATIO; //Step deeper
    { if ( (from)) memcpy(x1,x0,sizeof(double)*from); vect_mult_summ_vect(to-from,&x1[from],a[2],&p[from],&x0[from]); if (to!=n) memcpy(&x1[to],&x0[to],sizeof(double)*(n-to)); } // Unit shift
    if ((_i=funct(&f[2],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
    if (f[2]<f[1]) { f[0]=f[1], a[0]=a[1]; goto STEP_IN_DEEP; }
    //Else apply Brant method to find the (exact) minima
    }
  }
else
  {
  STEP_ON_SHALLOW: a[1]=a[2]*GOLDEN_SECTION; //Step upper
  { if ( (from)) memcpy(x1,x0,sizeof(double)*from); vect_mult_summ_vect(to-from,&x1[from],a[1],&p[from],&x0[from]); if (to!=n) memcpy(&x1[to],&x0[to],sizeof(double)*(n-to)); } // Unit shift
  if ((_i=funct(&f[1],n,x1,0x0,0x0,stack))!=YERROR_OK) return _i;
  if (f[1]>=f[0]) 
    {
    f[2]=f[1], a[2]=a[1]; 
    if ((a[2]-a[0])>itol/max_p) goto STEP_ON_SHALLOW; else return YERROR_NCONVERGED; 
    }
  //Else apply Brant method to find the (exact) minima 
  }
//Stage II. Do 1 step of quadratic interpolation
*lambda=.5*(a[0]+a[1]-(f[1]-f[0])*(a[2]-a[1])/(a[1]-a[0])/((f[2]-f[0])/(a[2]-a[0])-(f[1]-f[0])/(a[1]-a[0]))); //Can't be NAN or INF
{ if ( (from)) memcpy(x1,x0,sizeof(double)*from); vect_mult_summ_vect(to-from,&x1[from],*lambda,&p[from],&x0[from]); if (to!=n) memcpy(&x1[to],&x0[to],sizeof(double)*(n-to)); } // Unit shift
if ((_i=funct(ff,n,x1,g1,G,stack))!=YERROR_OK) return _i;
if (*ff>f[1])
  {//Disaster! Recalculate the previous best
  *lambda=a[1];
  { if ( (from)) memcpy(x1,x0,sizeof(double)*from); vect_mult_summ_vect(to-from,&x1[from],*lambda,&p[from],&x0[from]); if (to!=n) memcpy(&x1[to],&x0[to],sizeof(double)*(n-to)); } // Unit shift
  return funct(ff,n,x1,g1,G,stack);
  }
else return YERROR_OK;
}

//The same as above 'line_qsearch_square_fapproximation_only()' but with separate error vector (for routines like GDIIS)
unsigned int line_qsearch_e_square_fapproximation_only(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int n,double *x0,double *x1,double *g1,double *e1,double **G,double *p,unsigned int (*funct)(double *,unsigned int ,double *,double *,double *,double **,va_list ),va_list stack)
{
double a[3], f[3], max_p;
unsigned int _i;

//Stage 0. Find the biggest p component
max_p=fabs(*p), _i=n; while (_i--) if (max_p<fabs(p[_i])) max_p=fabs(p[_i]);
//Stage I. Define 3-points problem
f[0]=*ff, a[0]=0., a[2]=*lambda;
if (a[2]*max_p>max_variable_change) a[2]=0.9*max_variable_change/max_p;  
vect_mult_summ_vect(n,x1,a[2],p,x0);
f[2]=0.; //Init to avoid occasionaly NAN which is a flag for optimization with exclusions
if ((_i=funct(&f[2],n,x1,0x0,0x0,0x0,stack))!=YERROR_OK) return _i;
if (f[2]<f[0])
  {
  STEP_IN_DEEP: f[1]=f[2], a[1]=a[2];
  if (max_p*a[2]*GOLDEN_RATIO>max_variable_change) { *lambda=a[2]; return funct(ff,n,x1,g1,e1,G,stack); } //maximal extension reached  
  else 
    {
    a[2]*=GOLDEN_RATIO; //Step deeper
    vect_mult_summ_vect(n,x1,a[2],p,x0);
    if ((_i=funct(&f[2],n,x1,0x0,0x0,0x0,stack))!=YERROR_OK) return _i;
    if (f[2]<f[1]) { f[0]=f[1], a[0]=a[1]; goto STEP_IN_DEEP; }
    //Else apply Brant method to find the (exact) minima
    }
  }
else
  {
  STEP_ON_SHALLOW: a[1]=a[2]*GOLDEN_SECTION; //Step upper
  vect_mult_summ_vect(n,x1,a[1],p,x0);
  if ((_i=funct(&f[1],n,x1,0x0,0x0,0x0,stack))!=YERROR_OK) return _i;
  if (f[1]>=f[0]) 
    {
    f[2]=f[1], a[2]=a[1]; 
    if ((a[2]-a[0])>itol/max_p) goto STEP_ON_SHALLOW;
    else return YERROR_NCONVERGED; 
    }
  //Else apply Brant method to find the (exact) minima 
  }
//Stage II. Do one step of quadratic interpolation
*lambda=.5*(a[0]+a[1]-(f[1]-f[0])*(a[2]-a[1])/(a[1]-a[0])/((f[2]-f[0])/(a[2]-a[0])-(f[1]-f[0])/(a[1]-a[0]))); //Can't be NAN or INF
vect_mult_summ_vect(n,x1,*lambda,p,x0);
if ((_i=funct(ff,n,x1,g1,e1,G,stack))!=YERROR_OK) return _i;
if (*ff>f[1])
  {//Disaster! Recalculate the previous best
  *lambda=a[1], vect_mult_summ_vect(n,x1,*lambda,p,x0);
  return funct(ff,n,x1,g1,e1,G,stack);
  }
else return YERROR_OK;
}



/*********************** P O L A K _ R I B I E R E   (  C o n j u g a t e    G r a d i e n t)    P A R T ************************/

//This function performs optimization with Polak-Ribiere (conjugate gradient) method:
// p(k)=-g(k)+b(k-1)*p(k-1),    b(k-1)=max[0,((gk-gk-1).gk)/(gk-1.gk-1)].
unsigned int polak_ribiere(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **g,double **G,double *p,
                           unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),
                           unsigned int (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double *,double *,double *,double **,double *,
                                                       unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... )
{
unsigned int _i, step;
register double _beta, *_dp;
double lambda=1., gg, tol2;
va_list stack;

//#define DEBUG_N_BIGGEST 5
//unsigned int j, t, debug_biggest[DEBUG_N_BIGGEST];

//get stack access
va_start(stack,label);
//Stage I. Do one steepest-descent step
step=0, tol2=tol*tol; 
if ((_i=funct(fx,n,x[0x0],g[0x0],G,stack))!=YERROR_OK) { LABEL_EXTERNAL_FAILURE: if (label) yprintf(YPRINTF_ERROR,"%s encountered failure %s in the objective function.\n",label,get_yerrno(_i)); va_end(stack); return _i; }
else gg=calc_vect_norm(n,g[0x0]);
if (gg<tol2) { if (label) yprintf(YPRINTF_INFO,"%s initial gradient is beyond asked tolerance %16.8f (f=%32.24f |gg|=%32.24f).\n",label,tol,*fx,gg); va_end(stack); return YERROR_OK; }
else if (label) yprintf(YPRINTF_INFO,"%s n=%10.1d f=%32.24f |gg|=%32.24f lambda=%16.8f\n",label,step,*fx,gg,lambda);
LABEL_RELAUNCH: lambda=1.;
vect_inverse_sign(n,p,g[0x0]);
if ((_i=line_search(tol,itol,max_variable_change,&lambda,fx,n,x[0x0],x[0x1],g[0x1],G,p,funct,stack))!=YERROR_OK)
  {
//  //Debugin code
//  for (_i=0; _i<DEBUG_N_BIGGEST; _i++) debug_biggest[_i]=_i;
//  j=DEBUG_N_BIGGEST; while (--j) if (fabs(g[0][*debug_biggest])>fabs(g[0][debug_biggest[j]])) { t=*debug_biggest, *debug_biggest=debug_biggest[j], debug_biggest[j]=t; }
//  for (_i=DEBUG_N_BIGGEST; _i<n; _i++)
//    if (fabs(g[0][*debug_biggest])<fabs(g[0][_i]))
//      {
//      *debug_biggest=_i;
//      j=DEBUG_N_BIGGEST; while (--j) if (fabs(g[0][*debug_biggest])>fabs(g[0][debug_biggest[j]])) { t=*debug_biggest, *debug_biggest=debug_biggest[j], debug_biggest[j]=t; }
//      }
  LABEL_LINE_SEARCH_FAILURE: if (_i!=YERROR_NCONVERGED) goto LABEL_EXTERNAL_FAILURE;
  if (label) yprintf(YPRINTF_ERROR,"%s can't step along anti-gradient vector.\n",label); va_end(stack); return YERROR_SUSPICIOUS;
  }
else { _beta=1./gg, gg=calc_vect_norm(n,g[0x1]); _dp=x[0x0], x[0x0]=x[0x1], x[0x1]=_dp; _dp=g[0x0], g[0x0]=g[0x1], g[0x1]=_dp; } //Swap x[0] <-> x[1] and g[0] <-> g[1]
if (gg<tol2) { LABEL_EXIT_CONVERGENCE: if (label) yprintf(YPRINTF_INFO,"%s has converged to asked tol %16.8f in %10.1d steps (f(x)=%32.24f |ee|=%32.24f)!\n",label,tol,step,*fx,sqrt(gg)); va_end(stack); return YERROR_OK; }
else { if (label) yprintf(YPRINTF_INFO,"%s n=%10.1d f=%32.24f |g|=%32.24f lambda=%16.8f\n",label,step,*fx,sqrt(gg),lambda); }
//Stage II. Do CG opimization
while (++step<nsteps)
  {
  //Stgae II.1. Try CG step
  _beta*=gg-calc_vect_vect_scalar_product(n,g[0x0],g[0x1]);
  if (_beta>0.)
    {
    multiple_self_vect_scalar(n,p,_beta); subt_self_vect(n,p,g[0x0]); // Construct new dirrection
    if (calc_vect_vect_scalar_product(n,p,g[0x0])<=0.0) goto LABEL_BZERO;
    if ((_i=line_search(tol,itol,max_variable_change,&lambda,fx,n,x[0x0],x[0x1],g[0x1],G,p,funct,stack))!=YERROR_OK)
      { if (_i==YERROR_NCONVERGED) goto LABEL_RELAUNCH; else goto LABEL_EXTERNAL_FAILURE; }
    } 
  else 
    {
    LABEL_BZERO: vect_inverse_sign(n,p,g[0x0]); //Reset beta
    if ((_i=line_search(tol,itol,max_variable_change,&lambda,fx,n,x[0x0],x[0x1],g[0x1],G,p,funct,stack))!=YERROR_OK) goto LABEL_LINE_SEARCH_FAILURE;
    }
  //Stage II.2. Update beta on success
  _beta=1./gg, gg=calc_vect_norm(n,g[0x1]); _dp=x[0x0], x[0x0]=x[0x1], x[0x1]=_dp; _dp=g[0x0], g[0x0]=g[0x1], g[0x1]=_dp; //Swap x[0] <-> x[1] and g[0] <-> g[1]
  if (gg<tol2) goto LABEL_EXIT_CONVERGENCE;
  else { if (label) yprintf(YPRINTF_INFO,"%s n=%10.1d f=%32.24f |g|=%32.24f lambda=%16.8f\n",label,step,*fx,sqrt(gg),lambda); }
  }
if (label) yprintf(YPRINTF_WARNING,"Polak-Ribiere exhausted iterations count but haven't converged to the asked tolerance %32.24f :\n%s n=%10.1d f=%32.24f |g|=%32.24f\n",tol,label,step,*fx,sqrt(gg));
va_end(stack); return YERROR_NCONVERGED;  
}   



/************************************* L B F G S   P A R T *******************************************/



//This function updates l-bfgs structure: S, Y, IR and YY
inline unsigned int _update_lbfgs_append(unsigned int n,double *x0,double *x1,double *g0,double *g1,unsigned int size_m,double **S,double **Y,double **IR,double **YY)
{
register unsigned int _i;
register double _d;
subt_vect(n,S[size_m],x1,x0);
subt_vect(n,Y[size_m],g1,g0);
_d=calc_vect_vect_scalar_product(n,S[size_m],Y[size_m]); if (!_d) return YERROR_LEGAL; else IR[size_m][0]=1./_d; //rho :: 1/D
_i=size_m; while (_i--) YY[size_m][_i]=calc_vect_vect_scalar_product(n,S[_i],Y[size_m]);
_i=size_m; while (_i--) IR[_i][size_m-_i]=-IR[size_m][0]*calc_vect_vect_scalar_product(size_m-_i,IR[_i],&YY[size_m][_i]); 
_i=size_m; while (_i--) YY[size_m][_i]=YY[_i][size_m]=calc_vect_vect_scalar_product(n,Y[size_m],Y[_i]);  YY[size_m][size_m]=calc_vect_norm(n,Y[size_m]);
return YERROR_OK;
}
inline unsigned int _update_lbfgs_replace(unsigned int n,double *x0,double *x1,double *g0,double *g1,unsigned int size_m,double **S,double **Y,double **IR,double **YY)
{
register unsigned int _i, _j;
register double *_dp;
--size_m;
for (_dp=*S,  _i=0; _i<size_m; _i++)  S[_i]= S[_i+1];  S[size_m]=_dp;
for (_dp=*Y,  _i=0; _i<size_m; _i++)  Y[_i]= Y[_i+1];  Y[size_m]=_dp;
for (_i=0; _i<size_m; _i++) { for (_j=0; _j<size_m-_i; _j++) IR[_i][_j]=IR[_i+1][_j]; }
for (_dp=*YY, _i=0; _i<size_m; _i++) { for (YY[_i]=YY[_i+1], _j=0; _j<size_m; _j++) YY[_i][_j]=YY[_i][_j+1]; } YY[size_m]=_dp;
return _update_lbfgs_append(n,x0,x1,g0,g1,size_m,S,Y,IR,YY);
}

//This function (separately) updates l-bfgs structure: S, Y, IR and YY
unsigned int update_lbfgs(unsigned int n,double *x0,double *x1,double *g0,double *g1,unsigned int m,unsigned int size_m,double **S,double **Y,double **IR,double **YY)
{
if (m==size_m) return _update_lbfgs_replace(n,x0,x1,g0,g1,size_m,S,Y,IR,YY);
else           return _update_lbfgs_append(n,x0,x1,g0,g1,     m,S,Y,IR,YY); 
}

//This function calculates H-1.g with l-BFGS data
inline void calc_Hg_lbfgs(unsigned int n,double *g,double *p,unsigned int size_m,double gamma,double **S,double **Y,double **IR,double **YY,double *c,double *temp)
{
register unsigned int _i, _j;
//Stage I. Compute weights
_i=size_m; while (_i--) { c[_i]=calc_vect_vect_scalar_product(n,S[_i],g); } //S^T.g
_i=size_m; while (_i--) { c[size_m+_i]=-calc_vect_vect_scalar_product(size_m-_i,IR[_i],&c[_i]); } //-R^-1.S^T.g  { i.e. c1 }
_i=size_m; while (_i--) { temp[_i]=-c[size_m+_i]/IR[_i][0]-gamma*(calc_vect_vect_scalar_product(size_m,YY[_i],&c[size_m])+calc_vect_vect_scalar_product(n,Y[_i],g)); } //D.R^-1.S^T.g + gamma * (YY.R^-1.S^T.g-Y^T.g)
_i=size_m; while (_i--) { c[_i]=0., _j=_i; do { c[_i]+=IR[_j][_i-_j]*temp[_j]; } while (_j--); } //R^-T . ( D.R^-1.S^T.g + gamma * (YY.R^-1.S^T.g-Y^T.g) )  { i.e. c0 }
//Stage II. Synthetise the direction vector -H.g = [ S | gamma*Y ] . [ c0 | c1 ]
multiple_vect_scalar(n,p,S[0],-c[0]),  self_mult_subt_vect(n,p,gamma*c[size_m+0],Y[0]);
_i=size_m; while (--_i) { self_mult_subt_vect(n,p,c[_i],S[_i]), self_mult_subt_vect(n,p,gamma*c[size_m+_i],Y[_i]); }
}

//This function is reverse call of l-BFGS routine
unsigned int rlbfgs(unsigned int n,double *x0,double *x1,double *g0,double *g1,double *p,unsigned int *m,unsigned int size_m,double gamma,double **S,double **Y,double **IR,double **YY,double *c,double *temp)
{
register unsigned int _ylib_errno; 
if (*m==size_m) 
  {
  if ((_ylib_errno=_update_lbfgs_replace(n,x0,x1,g0,g1,size_m,S,Y,IR,YY))!=YERROR_OK) return _ylib_errno; //Stage I. Update S, Y, IR and YY
  else calc_Hg_lbfgs(n,g1,p,size_m,gamma,S,Y,IR,YY,c,temp);                                               //Stage II. Compute l-BFGS direction
  }
else
  {
  if ((_ylib_errno=_update_lbfgs_append(n,x0,x1,g0,g1,*m,S,Y,IR,YY))!=YERROR_OK) return _ylib_errno; //Stage I. Update S, Y, IR and YY
  else calc_Hg_lbfgs(n,g1,p,++(*m),gamma,S,Y,IR,YY,c,temp);                                         //Stage II. Compute l-BFGS direction
  } 
return YERROR_OK;
}

//This function performs unconstrained l-bfgs minimization (see Richard H. Byrd, Jorgr Nocedal and Robert B. Schnabel "Representation of quasi-newton matrices and their use in limited memory methods" )
unsigned int lbfgs_new(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **g,double **G,double *p,
                  unsigned int size_m,double **S,double **Y,double **IR,double **YY,double *c,double *temp,
                  unsigned int funct(double *fx,unsigned int n,double *x,double *g,double **G,va_list stack),
                  unsigned int (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double *,double *,double *,double **,double *,
                                              unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... )
{
register unsigned int _i;
register double *_dp;
unsigned int m, step;
double gg, gamma, val, lambda;
va_list stack;
//Stage 0. General sanity checks
if ( (!size_m)||(size_m==1) ) { if (label) yprintf(YPRINTF_ERROR,"%s encountered too few (<2, while optimal is 6-10) history vectors put into l-BFGS solver.\n",label); return YERROR_EXTERNAL_CODE; }
tol*=tol;
va_start(stack,label); //Get stack access
//Stage I.a. Init l-BFGS with one steepest-descent step
if ( (_i=funct(fx,n,x[0x0],g[0x0],G,stack))!=YERROR_OK) { LABEL_EXTERNAL_FAILURE: if (label) yprintf(YPRINTF_ERROR,"%s encountered failure %s in the objective function.\n",label,get_yerrno(_i)); va_end(stack); return _i; }
else { val=*fx, gg=calc_vect_norm(n,g[0x0]); }
if (gg<tol) { if (label) yprintf(YPRINTF_INFO,"%s found grdient below the desired tolerance %16.8f at the initial step (f=%32.24f |gg|=%32.24f).\n",label,tol,*fx,gg); va_end(stack); return YERROR_OK; } 
else { if (label) yprintf(YPRINTF_INFO,"%s init l-BFGS f=%32.24f |gg|=%32.24f\n",label,*fx,gg); }
//Stage I.b. Do the initial linear search
step=0;
LABEL_RESTART_LBFGS: m=0, lambda=1., vect_inverse_sign(n,p,g[0x0]); //gamma is used as lambda here
if ((_i=line_search(tol,itol,max_variable_change,&lambda,&val,n,x[0x0],x[0x1],g[0x1],G,p,funct,stack))!=YERROR_OK) 
  {
  if (_i!=YERROR_NCONVERGED) goto LABEL_EXTERNAL_FAILURE;
  if (label) yprintf(YPRINTF_ERROR,"%s can't step along errors vector.\n",label); va_end(stack); return YERROR_SUSPICIOUS;
  }
else { *fx=val, gamma=gg; gg=calc_vect_norm(n,g[0x1]); }
if (gg<tol) { LABEL_CONVERGED: if (label) yprintf(YPRINTF_INFO,"%s converged to tolerance %16.8f at %10.1d step (f=%32.24f |gg|=%32.24f).\n",label,tol,step,*fx,gg); va_end(stack); return YERROR_OK; } 
else { if (label) yprintf(YPRINTF_INFO,"%s step %10.1d f=%32.24f |gg|=%32.24f lambda=%16.12f\n",label,step,*fx,gg,lambda); }
//Stage II. l-BFGS itself
while (nsteps!=step++)
  {
  //Stage II.1. Update search dirrection
  gamma=(calc_vect_vect_scalar_product(n,x[0x1],g[0x1])-calc_vect_vect_scalar_product(n,x[0x0],g[0x1])-calc_vect_vect_scalar_product(n,x[0x1],g[0x0])+calc_vect_vect_scalar_product(n,x[0x0],g[0x0]))/(gg-2.*calc_vect_vect_scalar_product(n,g[0x0],g[0x1])+gamma);
  _i=rlbfgs(n,x[0x0],x[0x1],g[0x0],g[0x1],p,&m,size_m,gamma,S,Y,IR,YY,c,temp);  //Stage II.2. Conjugate directions vectors
  self_mult_subt_vect(n,p,gamma,g[0x1]);
  { _dp=x[0x0], x[0x0]=x[0x1], x[0x1]=_dp; _dp=g[0x0], g[0x0]=g[0x1], g[0x1]=_dp; }
  //Stage II.3. Try the suggested point
  lambda=1.;
  if ((_i=line_search(tol,itol,max_variable_change,&lambda,&val,n,x[0x0],x[0x1],g[0x1],G,p,funct,stack))!=YERROR_OK) 
    {
    if (_i!=YERROR_NCONVERGED) goto LABEL_EXTERNAL_FAILURE;
    if (label) yprintf(YPRINTF_ERROR,"%s can't step along l-bfgs errors vector, trying to restart.\n",label);
    goto LABEL_RESTART_LBFGS;
    }
  else
    {//Stage II.4. Accept the new point and check for convergence
    *fx=val; gamma=gg, gg=calc_vect_norm(n,g[0x1]);
    if (gg<tol) { _dp=x[0x0], x[0x0]=x[0x1], x[0x1]=_dp; goto LABEL_CONVERGED; } 
    else { if (label) yprintf(YPRINTF_INFO,"%s step %10.1d f=%32.24f |gg|=%32.24f lambda=%16.12f\n",label,step,*fx,gg,lambda); }
    }
  }
//Stage III. Exit on NCONVERGED
if (label) yprintf(YPRINTF_INFO,"%s excided allowed amount of iterations (%10.1d) but hasn't converged to desired tolerance %32.24f (f=%32.24f |gg|=%32.24f)\n",label,nsteps,tol,*fx,gg);
_dp=x[0x0], x[0x0]=x[0x1], x[0x1]=_dp; va_end(stack); return YERROR_NCONVERGED;
}



//DEPRECATED!!!
//This function performs unconstrained l-bfgs minimization (see Richard H. Byrd, Jorgr Nocedal and Robert B. Schnabel "Representation of quasi-newton matrices and their use in limited memory methods" )
// n - number of variables, m - number of corrections (3...7 recomended), x - two sets of variables-coordinates, g - two sets of gradient, funct - function and its derivative, tol - tolerance of dfunct
//Note. The S and Y matrices are trasposed to use processors cache more effectively
char lbfgs(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **g,double **G,double *z,unsigned int m,char (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),
           char (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double **,double **,double **,double *,char (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... )
{
double **S, **Y, **YY, **R, *d,*sg[0x2], *yg[0x2], *p, gg[2], lambda, gamma;
unsigned int _i, _j, _m, step=0;
double _d;
va_list stack;
void *vp1, *vp2;

//get stack access
va_start(stack,label);

//Stage A. Init memory and perform one quasi steepest descent steps 
if (!(funct(fx,n,x[0],g[0],G,stack)))
  {
  EXTERNAL_FAIL: if (label) printf("Sorry, lbfgs failed due to failure in the function under scope.\n");
  return FALSE;
  }
tol*=tol;
if ((gg[0]=calc_vect_norm(n,g[0x0]))<tol) 
  {
  RETURN_SUCCESS: ;
  if ((label)) printf("The system solved with l-bfgs is already converged to requested gradient norm %32.24f (now |g|=%32.24f f=%32.24f)\n",tol,gg[0],*fx);
  ylib_errno=YERR_OK;
  return TRUE;
  }
else if ((label)) printf("Launching l-bfgs from init point:\n%s n=%10.1d f=%32.24f |g|=%32.24f\n",label,step,*fx,sqrt(gg[0]));

//Stage 0. Prepare all
//Stage B. Map memory
if (!(p=(void*)malloc( (sizeof(double*)+sizeof(double)*n)*m*2+sizeof(double)*m*2+
                       (sizeof(double*)+sizeof(double)*m)*m+sizeof(double*)*m+sizeof(double)*m*(m+1)/2+
                        sizeof(double)*m+sizeof(double)*m*4+sizeof(double)*m) )) { ylib_errno=YERROR_MEMORY; return FALSE; }
S=(void*)p+sizeof(double)*m*2;
Y=(void*)S+(sizeof(double*)+sizeof(double)*n)*m;
S[0]=(void*)S+sizeof(double*)*m;
Y[0]=(void*)Y+sizeof(double*)*m;
for (_i=1;_i<m;_i++) { S[_i]=(void*)S[_i-1]+sizeof(double)*n; Y[_i]=(void*)Y[_i-1]+sizeof(double)*n; }
YY=(void*)Y+(sizeof(double*)+sizeof(double)*n)*m;
R=(void*)YY+(sizeof(double*)+sizeof(double)*m)*m;
YY[0]=(void*)YY+sizeof(double*)*m;
R[0]=(void*)R+sizeof(double*)*m;
for (_i=1;_i<m;_i++) { YY[_i]=(void*)YY[_i-1]+sizeof(double)*m; R[_i]=(void*)R[_i-1]+sizeof(double)*(m-_i+1); }
d=(void*)R+sizeof(double*)*m+sizeof(double)*m*(m+1)/2;
sg[0]=(void*)d+sizeof(double)*m;
sg[1]=(void*)sg[0]+sizeof(double)*m;
yg[0]=(void*)sg[1]+sizeof(double)*m;
yg[1]=(void*)yg[0]+sizeof(double)*m;

//Stage C. optimize with l-bfgs (with H==I) initially
RESTART_LBFGS:
//Fill first S, Y and gg 
_i=n, z+=n, g[0]+=n, gg[0]=sqrt(gg[0]); while(_i--) { z--, g[0]--, *z=-*g[0]/gg[0]; }
lambda=.5;
if (!(line_search(tol,itol,max_variable_change,&lambda,fx,n,x,g,G,z,funct,stack)))
  {
  if (ylib_errno!=YERROR_NCONVERGED) goto EXTERNAL_FAIL;
  if (label) printf("Something strange here: lbfgs cant step even along gradient vector...\n");
  ylib_errno=YERROR_SUSPICIOUS;
  return FALSE;
  }
else step++;
gg[1]=gg[0], gg[0]=calc_vect_norm(n,g[0x0]);
if (gg[0]<tol) goto RETURN_SUCCESS;
if ((label)) { _j=0, _i=n; while (--_i) if (fabs(g[0][_i])>fabs(g[0][_j])) _j=_i; printf("Launching l-bfgs from init point:\n%s n=%10.1d f=%32.24f |g|=%32.24f max_g=%32.24f\n",label,step,*fx,sqrt(gg[0]),g[0][_j]); }

//Create S and Y + calculate g(k-1).g(k)
_d=0., S[0]+=n, x[0]+=n, x[1]+=n, Y[0]+=n, g[0]+=n, g[1]+=n;
_j=n; while (_j--) { S[0]--, x[0]--, x[1]--, Y[0]--, g[0]--, g[1]--, *S[0]=*x[0]-*x[1], *Y[0]=*g[0]-*g[1], _d+=*g[0]**g[1]; }
//Create S.g and Y.g
*sg[0]=0., *yg[0]=0., S[0]+=n, Y[0]+=n, g[0]+=n;
_j=n; while (_j--) { S[0]--, Y[0]--, g[0]--, *sg[0]+=*S[0]**g[0], *yg[0]+=*Y[0]**g[0]; }
//Create d[0], R[0][0], YY[0][0]
*d=0., S[0]+=n, Y[0]+=n, _j=n; while (_j--) { S[0]--, Y[0]--, *d+=*S[0]**Y[0]; } R[0][0]=*d;
YY[0][0]=gg[1]+gg[0]-2.*_d; //YY=g(k-1).g(k-1)+g(k).g(k)-2*g(k-1).g(k)

//Do enering step
lambda=.5, _i=n, z+=n, g[0]+=n, gg[0]=sqrt(gg[0]); while(_i--) { z--, g[0]--, *z=-*g[0]/gg[0]; }
if (!(line_search(tol,itol,max_variable_change,&lambda,fx,n,x,g,G,z,funct,stack)))
  {
  if (ylib_errno!=YERROR_NCONVERGED) goto EXTERNAL_FAIL;
  if (label) printf("Something strange: lbfgs can't step even along gradient vector...\n");
  ylib_errno=YERROR_SUSPICIOUS;
  return FALSE;
  }
else step++;
gg[1]=gg[0], gg[0]=calc_vect_norm(n,g[0x0]);

//Do l-BFGS calculatios
_m=0; 
while (step<=nsteps)
  {
  if (gg[0]<tol)
    {
    if (label) { _j=0, _i=n; while (--_i) if (fabs(g[0][_i])>fabs(g[0][_j])) _j=_i; printf("l-bfgs converged to requested gradient norm %32.24f :\n%s n=%10.1d f=%32.24f |g|=%32.24f max_g=%32.24f\n",sqrt(tol),label,step,*fx,sqrt(gg[0]),g[0][_j]); }
    break;
    }
  else if ((label)) { _j=0, _i=n; while (--_i) if (fabs(g[0][_i])>fabs(g[0][_j])) _j=_i; printf("%s n=%10.1d f=%32.24f |g|=%32.24f max_g=%32.24f lambda=%16.8f\n",label,step,*fx,sqrt(gg[0]),g[0][_j],lambda); }
  //Prepare updating
  if (_m==m-1)
    {
    //Rebuilt S, Y, YY, R, sg, yg 
    vp1=*S, vp2=*Y;
    _i=m;
    while (--_i)
      {
      for (_j=_i;_j<m;_j++) { R[_i-1][_j-_i]=R[_i][_j-_i]; YY[_i-1][_j-1]=YY[_j-1][_i-1]=YY[_i][_j]; }
      S[_i-1]=S[_i];
      Y[_i-1]=Y[_i];
      d[_i-1]=d[_i];
      sg[0][_i-1]=sg[0][_i];
      yg[0][_i-1]=yg[0][_i];
      }
    S[_m]=vp1, Y[_m]=vp2;
    }
  else _m++;
  //Stage 1. Update Sk and Yk
  S[_m]+=n, x[0]+=n, x[1]+=n, Y[_m]+=n, g[0]+=n, g[1]+=n;
  _j=n; while (_j--) { S[_m]--, x[0]--, x[1]--, Y[_m]--, g[0]--, g[1]--; *S[_m]=*x[0]-*x[1]; *Y[_m]=*g[0]-*g[1]; }

  //Stage 2. Compute g(k)^T.g(k), S(k)^T.g(k), Y(k)^T.g(k)
  vp1=sg[0], vp2=yg[0], sg[0]=sg[1], yg[0]=yg[1], sg[1]=vp1, yg[1]=vp2;
  _i=_m+1, sg[0]+=_i, yg[0]+=_i;
  while (_i--)
    {
    sg[0]--, yg[0]--, *sg[0]=0.,  *yg[0]=0., S[_i]+=n, Y[_i]+=n, g[0]+=n;
    _j=n; while (_j--) { S[_i]--, Y[_i]--, g[0]--, *sg[0]+=*S[_i]**g[0], *yg[0]+=*Y[_i]**g[0]; }
    }

  //Stage 3. Calculate s(k-1)^T.g(k-1) and y(k-1)^T.g(k-1) normally
  sg[1][_m]=0., yg[1][_m]=0., S[_m]+=n, Y[_m]+=n, g[1]+=n; 
  _j=n; while (_j--) { S[_m]--, Y[_m]--, g[1]--, sg[1][_m]+=*S[_m]**g[1], yg[1][_m]+=*Y[_m]**g[1]; }

  //Stage 4. Update R, YY and D
  _i=_m;
  while (_i--)
    {
     R[_i][_m-_i]=sg[0][_i]-sg[1][_i];            //keep R triangle
    YY[_i][_m]=YY[_m][_i]=yg[0][_i]-yg[1][_i];
    }
  R[_m][0]=d[_m]=sg[0][_m]-sg[1][_m]; //keep R triangle
  g[0]+=n; g[1]+=n; _d=0.; _i=n; while(_i--) { g[0]--; g[1]--; _d+=(*g[0]-*g[1])**g[0]; }
  YY[_m][_m]=gg[1]-gg[0]+2.*_d;

  //Stage 5. Calculate positive gamma: gamma=|y(k-1).s(k-1)/y(k-1).y(k-1)|
  gamma=fabs(d[_m]/YY[_m][_m]);

  //Stage 6. Calculate p
  //Strage A. Calculate R^-1.(S^T.g) vector
  for (_i=0; _i<_m+1; _i++) p[_m+1+_i]=sg[0][_i]; 
  bsubstitute_tdU(_m+1,R,&p[_m+1]);
  //Stage B. Calculate lambda*((D/lambda+YY).p[m+1]-Yg)
  for (_i=0;_i<=_m;_i++)
    {
    p[_i]=p[_m+1+_i]*d[_i]-gamma*yg[0][_i];  //Use (probably) g[1] as temporary storage
    for (_j=0;_j<=_m;_j++) p[_i]+=gamma*p[_m+1+_j]*YY[_i][_j];
    }
  bsubstitute_tdUT(_m+1,R,p);

  //Stage 7. Compute z=H.g
  for (_j=0;_j<=_m;_j++) { S[_j]+=n; Y[_j]+=n; }
  z+=n; g[0]+=n;
  _i=n;
  while (_i--)
    {
    z--, g[0]--;
    *z=gamma**g[0];
    for (_j=0;_j<=_m;_j++) { S[_j]--, Y[_j]--, (*z)+=*(S[_j])*p[_j]-gamma**(Y[_j])*p[_m+1+_j]; }
    }
  if (sqrd(calc_vect_vect_scalar_product(n,g[0],z))<0.)
    {
    if (label) printf("Enhances direction failed, step along gradien vector has done instead\n");
    vp1=S[0], vp2=Y[0], S[0]=S[_m], Y[0]=Y[_m], S[_m]=vp1, Y[_m]=vp2;
    sg[0][0]=sg[0][_m], yg[0][0]=yg[0][_m], d[0]=d[_m], R[0][0]=R[_m][0], YY[0][0]=YY[_m][_m];
    goto RESTART_LBFGS; 
    }
  //Stage 8. Do line search
  _i=n, z+=n; while(_i--) { z--, *z=-*z; }
  if (!(line_search(tol,itol,max_variable_change,&lambda,fx,n,x,g,G,z,funct,stack)))
    {//Try to make a step along gradient
    if (label) printf("Enhances direction failed, triyng step along gradien vector instead\n");
    goto RESTART_LBFGS; //forget all previous directions
    }
  else step++;
  gg[1]=gg[0], gg[0]=calc_vect_norm(n,g[0x0]);
  }

//Free memory and exit
va_end(stack);
free(p);
ylib_errno=YERR_OK;
return TRUE;
}



/*************************************   D I I S   P A R T   *******************************************/


//This function updates diis structure: x, e and B
inline void _update_diis_append(unsigned int n,unsigned int      m,double **e,double **B,double *val)
{
register unsigned int _i;
B[m][m]=calc_vect_norm(n,e[m]), _i=m; while (_i--) B[_i][m]=B[m][_i]=calc_vect_vect_scalar_product(n,e[_i],e[m]);
}
inline void _update_diis_replace(unsigned int n,unsigned int size_e,double **x,double **e,double **B,double *val)
{
register unsigned int _i, _j;
register double *_dp;
_i=size_e-1, _j=0; while (--_i) if (val[_i]>val[_j]) _j=_i;
val[_j]      =val[size_e-1]; _dp=x[_j      ], x[_j      ]=x[size_e-1], x[size_e-1]=_dp; _dp=e[_j      ], e[_j      ]=e[size_e-1], e[size_e-1]=_dp;
B[_j][_j]=B[size_e-1][size_e-1], _i=size_e-1; while (_i--) if (_j!=_i) B[_i][_j]=B[_j][_i]=B[size_e-1][_i];
val[size_e-1]=val[size_e  ]; _dp=x[size_e-1], x[size_e-1]=x[size_e  ], x[size_e  ]=_dp; _dp=e[size_e-1], e[size_e-1]=e[size_e  ], e[size_e  ]=_dp;
_i=size_e; while (_i--) B[_i][_j]=B[_j][_i]=calc_vect_vect_scalar_product(n,e[_i],e[_j]);
}

//This function (separately) updates l-bfgs structure: S, Y, IR and YY
void update_diis(unsigned int n,unsigned int m,unsigned int size_e,double **x,double **e,double **B,double *val)
{
if (m==size_e) _update_diis_replace(n,   m,x,e,B,val);
else           _update_diis_append(n,size_e,e,B,val); 
}

//The reversibly-calling diis - it computes DIIS interpolation xi for given set of coordinate vectors x and error vectors e
inline unsigned int _weight_diis(unsigned int size_e,double **B,double *c)
{
register unsigned int _i;
//Stage I. Construct constants vector
_i=size_e; while (_i--) { c[_i]=0., B[_i][size_e]=B[size_e][_i]=-1.; }  c[size_e]=-1., B[size_e][size_e]=0.;  
//Stage II. Solve DIIS interpolation
return gauss_solve_dmatrix(size_e+1,B,c,.1*SMALL2);
}

//The reversibly-calling diis - it computes DIIS interpolation weigths for given set of error vectors e
unsigned int rdiis(unsigned int n,unsigned int *m,unsigned int size_e,double **x,double **e,double **B,double *val,double *c)
{
if (*m+1==size_e) 
  {
  _update_diis_replace(n,size_e,x,e,B,val); //Stage I. Update B, e and x
  return _weight_diis(size_e,B,c);           //Stage II. Compute diis weights
  }
else
  {
  _update_diis_append(n, (*m)++,x,B,val); //Stage I. Update B, e and x
  return _weight_diis(    *m,B,c);         //Stage II. Compute diis weights
  } 
}

//This routine create diis vector
void synthesize_diis_vector(unsigned int n,unsigned int size_e,double *vi,double **v,double *c)
{
multiple_vect_scalar(n,vi,*v,*c); while (--size_e) self_mult_summ_vect(n,vi,c[size_e],v[size_e]); 
}

//This function minimizes function using Direct Inversion in the Iterative (sub-)Space method [DIIS]
unsigned int diis(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **G,
                  unsigned int size_e,double **e,double **B,double *c,double *val,
                  unsigned int funct(double *fx,unsigned int n,double *x,double *e,double **G,va_list stack),
                  unsigned int (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double *,double *,double *,double **,double *,
                                              unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... )
{
register unsigned int _i;
unsigned int m, step;
register double *_dp;
double tol2, lambda, ee;
va_list stack;

//Stage 0. General sanity checks
if ( (!size_e)||(size_e==1) ) { if (label) yprintf(YPRINTF_ERROR,"%s encountered too few (<2) errors vectors put into DIIS solver.\n",label); return YERROR_EXTERNAL_CODE; }
if (nsteps<size_e+2) { if (label) yprintf(YPRINTF_ERROR,"Too few (<%10.1d) steps are requested for %s - I can't even prime the DIIS.\n",size_e+2,label); return YERROR_EXTERNAL_CODE; }
tol2=tol*tol; va_start(stack,label); //Get stack access
//Stage I.a. Init DIIS with 1 steepest-descent direction
if ( (_i=funct(fx,n,x[0x0],e[0x0],0x0,stack))==YERROR_OK) { val[0x0]=*fx, B[0x0][0x0]=ee=calc_vect_norm(n,e[0x0]); self_vect_inverse_sign(n,e[m-1]); }
else { LABEL_EXTERNAL_FAILURE: if (label) yprintf(YPRINTF_ERROR,"%s encountered failure %s in the objective function.\n",label,get_yerrno(_i)); va_end(stack); return _i; }
if (ee<tol2) { if (label) yprintf(YPRINTF_INFO,"%s init f=%32.24f |ee|=%32.24f\n",label,*fx,ee); }
else { if (label) yprintf(YPRINTF_INFO,"%s found error vector beyond desired tolerance %16.8f at the initial step (f=%32.24f |ee|=%32.24f).\n",label,tol,*fx,ee); va_end(stack); return YERROR_OK; } 
//Stage I.b. Do the line search
step=0, m=1, lambda=1.;
if ((_i=line_search(tol,itol,max_variable_change,&lambda,fx,n,x[m-1],x[m],e[m],0x0,e[m-1],funct,stack))==YERROR_OK) { val[0x1]=*fx, ee=calc_vect_norm(n,e[m]); self_vect_inverse_sign(n,e[m]); }
else { LABEL_LINEAR_SEARCH_FAILURE:
     if (_i!=YERROR_NCONVERGED) goto LABEL_EXTERNAL_FAILURE; 
     if (label) yprintf(YPRINTF_ERROR,"%s can't step along errors vector.\n",label); va_end(stack); return YERROR_SUSPICIOUS;
     }
{ _dp=x[0x0], x[0x0]=x[0x1], x[0x1]=_dp; _dp=e[0x0], e[0x0]=e[0x1], e[0x1]=_dp; }
if (ee<tol2) { if (label) yprintf(YPRINTF_INFO,"%s converged to desired tolerance %16.8f at the initial iteration (f=%32.24f |ee|=%32.24f).\n",label,tol,*fx,ee); va_end(stack); return YERROR_OK; }
else { if (label) yprintf(YPRINTF_INFO,"%s iteration %10.1d %sf=%32.24f |ee|=%32.24f\n",label,step,(m<size_e) ? "(priming)" : "",*fx,ee); }
//Stage II. DIIS loop
while (nsteps!=step++)
  {
  //Stage II.A. Try DIIS extrapolation 
  if (rdiis(n,&m,size_e,x,e,B,val,c)==YERROR_OK) //it also updates m-counter
    {
    //II.A.1. Evaluate the trial point
    synthesize_diis_vector(n,m,x[m+1],x,c), synthesize_diis_vector(n,m,e[m+1],e,c), subt_self_vect(n,x[m+1],e[m+1]); //??????????????????????????????????????????????????????
    if ((_i=funct(&val[m+1],n,x[m+1],e[m+1],0x0,stack))!=YERROR_OK) goto LABEL_EXTERNAL_FAILURE;
    //III.A.2.a Update
    if (val[m+1]<*fx) *fx=val[m+1]; else { if (label) yprintf(YPRINTF_INFO,"%s iteration (failure) %10.1d f=%32.24f |ee|=%32.24f\n",label,step,*fx,ee); goto LABEL_LINEAR_SEARCH; } 
    }
  else   //Stage II.B. Failure, do linear search instead
    {
    LABEL_LINEAR_SEARCH: lambda=1.;
    if ((_i=line_search(tol,itol,max_variable_change,&lambda,fx,n,x[m],x[m+1],e[m+1],0x0,e[m],funct,stack))!=YERROR_OK) val[m+1]=*fx; else goto LABEL_LINEAR_SEARCH_FAILURE; 
    }
  //Stage III. Check convergence  
  ee=calc_vect_norm(n,e[m+1]); self_vect_inverse_sign(n,e[m+1]);
  if (ee>tol2) { if (label) yprintf(YPRINTF_INFO,"%s iteration %10.1d %sf=%32.24f |ee|=%32.24f\n",label,step,(m+1<size_e) ? "(priming) " : "",*fx,ee); }
  else 
    {
    if (label) yprintf(YPRINTF_INFO,"%s converged to desired tolerance %16.8f in %10.1d iterations (f=%32.24f |ee|=%32.24f).\n",label,tol,step,*fx,ee);
    { _dp=x[0x0], x[0x0]=x[m], x[m]=_dp; } va_end(stack); return YERROR_OK;
    }
  }
//Stage IV. Exit on NCONVERGED
if (label) yprintf(YPRINTF_INFO,"%s excided allowed amount of iterations (%10.1d) but hasn't converged to desired tolerance %32.24f (f=%32.24f |ee|=%32.24f)\n",label,nsteps,tol,*fx,ee);
_i=m; while (_i--) if (val[_i]<val[m]) m=_i;
if ( (m)) { _dp=x[0x0], x[0x0]=x[m], x[m]=_dp; } va_end(stack); return YERROR_NCONVERGED;
}

inline void _update_gdiis_replace(unsigned int n,unsigned int size_e,double **x,double **g,double **e,double **B,double *val)
{
register unsigned int _i, _j;
register double *_dp;
_i=size_e-1, _j=0; while (--_i) if (val[_i]>val[_j]) _j=_i;
val[_j]      =val[size_e-1]; _dp=x[_j      ], x[_j      ]=x[size_e-1], x[size_e-1]=_dp; _dp=g[_j      ], g[_j      ]=g[size_e-1], g[size_e-1]=_dp; _dp=e[_j      ], e[_j      ]=e[size_e-1], e[size_e-1]=_dp;
B[_j][_j]=B[size_e-1][size_e-1], _i=size_e-1; while (_i--) if (_j!=_i) B[_i][_j]=B[_j][_i]=B[size_e-1][_i];
val[size_e-1]=val[size_e  ]; _dp=x[size_e-1], x[size_e-1]=x[size_e  ], x[size_e  ]=_dp; _dp=g[size_e-1], g[size_e-1]=g[size_e  ], g[size_e  ]=_dp; _dp=e[size_e-1], e[size_e-1]=e[size_e  ], e[size_e  ]=_dp;
_i=size_e; while (_i--) B[_i][_j]=B[_j][_i]=calc_vect_vect_scalar_product(n,e[_i],e[_j]);
}

//The reversibly-calling gdiis - it computes DIIS interpolation weigths for given set of error vectors e
unsigned int rgdiis(unsigned int n,unsigned int *m,unsigned int size_e,double **x,double **g,double **e,double **B,double *val,double *c)
{
if (*m+1==size_e) 
  {
  _update_gdiis_replace(n,size_e,x,g,e,B,val); //Stage I. Update B, e and x
  return _weight_diis(size_e,B,c);             //Stage II. Compute diis weights
  }
else
  {
  _update_diis_append(n, (*m)++,e,B,val); //Stage I. Update B, e and x
  return _weight_diis(    *m,B,c);         //Stage II. Compute diis weights
  } 
}

//This routine minimizes function using GDIIS approach
unsigned int gdiis(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **g,double **G,double *p,
                   unsigned int size_e,double **e,double **B,double *c_diis,double *val,
                   unsigned int size_m,double **S,double **Y,double **IR,double **YY,double *c_lbfgs,double *temp_lbfgs,
                   unsigned int funct(double *fx,unsigned int n,double *x,double *g,double **G,va_list stack),
                   unsigned int (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double *,double *,double *,double **,double *,
                                               unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... )
 {
register unsigned int _i;
unsigned int step, m_diis, m_lbfgs, lbfgs_error;
register double *_dp;
double gg, gamma, lambda, tol2;
va_list stack;

//Stage 0. General sanity checks
if ( (!size_e)||(size_e==1) ) { if (label) yprintf(YPRINTF_ERROR,"%s encountered too few (<2) errors vectors put into GDIIS solver.\n",label); return YERROR_EXTERNAL_CODE; }
if ( (!size_m)||(size_m==1) ) { if (label) yprintf(YPRINTF_ERROR,"%s encountered too few (<2) l-bfgs vectors put into GDIIS solver.\n",label); return YERROR_EXTERNAL_CODE; }
if (nsteps<size_e+2) { if (label) yprintf(YPRINTF_ERROR,"Too few (<%10.1d) steps are requested for %s - I can't even prime the DIIS.\n",size_e+2,label); return YERROR_EXTERNAL_CODE; }
tol2=tol*tol; 
va_start(stack,label); //Get stack access
//Stage I.a. Init DIIS with 1 steepest-descent direction
if ( (_i=funct(fx,n,x[0x0],g[0x0],0x0,stack))==YERROR_OK) { val[0x0]=*fx, gg=calc_vect_norm(n,g[0x0]); self_vect_inverse_sign(n,e[m_diis-1]); }
else { LABEL_EXTERNAL_FAILURE: if (label) yprintf(YPRINTF_ERROR,"%s encountered failure %s in the objective function.\n",label,get_yerrno(_i)); va_end(stack); return _i; }
if (gg<tol2) { if (label) yprintf(YPRINTF_INFO,"%s init f=%32.24f |gg|=%32.24f\n",label,*fx,gg); }
else { if (label) yprintf(YPRINTF_INFO,"%s found gradient vector beyond desired tolerance %16.8f at the initial step (f=%32.24f |gg|=%32.24f).\n",label,tol,*fx,gg); va_end(stack); return YERROR_OK; } 
//Stage I.b. Do the line search
step=0, m_lbfgs=0, m_diis=0, lambda=1.; vect_inverse_sign(n,p,g[0x0]);
if ((_i=line_search(tol,itol,max_variable_change,&lambda,fx,n,x[m_diis],x[m_diis+1],g[m_diis+1],0x0,g[m_diis],funct,stack))==YERROR_OK) { val[m_diis+1]=*fx, gamma=gg, gg=calc_vect_norm(n,g[m_diis+1]); subt_vect(n,e[m_diis],x[m_diis+1],x[m_diis]); }
else { LABEL_LINEAR_SEARCH_FAILURE:
     if (_i!=YERROR_NCONVERGED) goto LABEL_EXTERNAL_FAILURE; 
     if (label) yprintf(YPRINTF_ERROR,"%s can't step along errors vector.\n",label); va_end(stack); return YERROR_SUSPICIOUS;
     }
if (gg<tol2) { if (label) yprintf(YPRINTF_INFO,"%s converged to desired tolerance %16.8f at the initial iteration (f=%32.24f |gg|=%32.24f).\n",label,tol,*fx,gg); va_end(stack); return YERROR_OK; }
else { if (label) yprintf(YPRINTF_INFO,"%s iteration %10.1d %sf=%32.24f |gg|=%32.24f\n",label,step,(m_diis<size_e) ? "(priming)" : "",*fx,gg); }
//Stage II. GDIIS loop
lbfgs_error=YERROR_LEGAL; //Not yet initialized
while (nsteps!=step++)
  {
  //Stage II.1. Updates
  //Stage II.1.a. Update l-bfgs matrices
  if ((lbfgs_error=update_lbfgs(n,x[m_diis],x[m_diis+1],g[m_diis],g[m_diis+1],m_lbfgs,size_m,S,Y,IR,YY))==YERROR_OK) { if (m_lbfgs<size_m) m_lbfgs++; } 
  //Stage II.1.b. Compute new e-vector estimation
  gamma=(calc_vect_vect_scalar_product(n,x[m_diis+1],g[m_diis])-calc_vect_vect_scalar_product(n,x[m_diis],g[m_diis+1])-calc_vect_vect_scalar_product(n,x[m_diis+1],g[m_diis])+calc_vect_vect_scalar_product(n,x[m_diis],g[m_diis]))/(gg-2.*calc_vect_vect_scalar_product(n,g[m_diis],g[m_diis+1])+gamma);
  if (lbfgs_error!=YERROR_OK) multiple_vect_scalar(n,e[m_diis+1],g[m_diis+1],-gamma); else { calc_Hg_lbfgs(n,g[m_diis+1],e[m_diis+1],m_lbfgs,gamma,S,Y,IR,YY,c_lbfgs,temp_lbfgs); self_mult_subt_vect(n,e[m_diis+1],gamma,g[m_diis+1]); }
  //Stage II.1.c. Update DIIS 
  if (rgdiis(n,&m_diis,size_e,x,g,e,B,val,c_diis)==YERROR_OK) //it also updates m-counter
    {
    //Stage II.2. Evaluate the trial point
    //Stage II.2.a. Compute DIIS extrapolations
    synthesize_diis_vector(n,m_diis,g[m_diis+1],g,c_diis), synthesize_diis_vector(n,m_diis,x[m_diis+1],x,c_diis), synthesize_diis_vector(n,m_diis,e[m_diis+1],e,c_diis), subt_self_vect(n,x[m_diis+1],e[m_diis+1]); //?????????????????????????
    //Stage II.2.b. Correct it with l-bfgs hessian
    if (lbfgs_error!=YERROR_OK)  multiple_vect_scalar(n,p,g[m_diis+1],-gamma);
    else { calc_Hg_lbfgs(n,g[m_diis+1],e[m_diis+1],m_lbfgs,gamma,S,Y,IR,YY,c_lbfgs,temp_lbfgs); self_mult_subt_vect(n,e[m_diis+1],gamma,g[m_diis+1]); }
    //Stage II.2.c. Compute function value
    if ((_i=funct(&val[m_diis+1],n,x[m_diis+1],e[m_diis+1],0x0,stack))!=YERROR_OK) goto LABEL_EXTERNAL_FAILURE;
    //III.A.2.a Update if possible
    if (val[m_diis+1]<*fx) *fx=val[m_diis+1]; else { if (label) yprintf(YPRINTF_INFO,"%s iteration (failure) %10.1d f=%32.24f |gg|=%32.24f\n",label,step,*fx,gg); goto LABEL_LINEAR_SEARCH; } 
    }
  else
    { //Stage II.2.b. Do linesearch along (l-bfgs) error vector
    LABEL_LINEAR_SEARCH: lambda=1.;
    if ((_i=line_search(tol,itol,max_variable_change,&lambda,fx,n,x[m_diis],x[m_diis+1],e[m_diis+1],0x0,e[m_diis],funct,stack))!=YERROR_OK) val[m_diis+1]=*fx; else goto LABEL_LINEAR_SEARCH_FAILURE; 
    }
  //Stage III. Check convergence  
  gamma=gg, gg=calc_vect_norm(n,g[m_diis+1]); subt_vect(n,e[m_diis],x[m_diis+1],x[m_diis]);
  if (gg>tol2) { if (label) yprintf(YPRINTF_INFO,"%s iteration %10.1d %sf=%32.24f |gg|=%32.24f\n",label,step,(m_diis+1<size_e) ? "(priming) " : "",*fx,gg); }
  else 
    {
    if (label) yprintf(YPRINTF_INFO,"%s converged to desired tolerance %16.8f in %10.1d iterations (f=%32.24f |gg|=%32.24f).\n",label,tol,step,*fx,gg);
    { _dp=x[0x0], x[0x0]=x[m_diis], x[m_diis]=_dp; } va_end(stack); return YERROR_OK; 
    }
  }
//Stage IV. Exit on NCONVERGED
if (label) yprintf(YPRINTF_INFO,"%s excided allowed amount of iterations (%10.1d) but hasn't converged to desired tolerance %32.24f (f=%32.24f |gg|=%32.24f)\n",label,nsteps,tol,*fx,gg);
_i=m_diis; while (_i--) if (val[_i]<val[m_diis]) m_diis=_i;
if ( (m_diis)) { _dp=x[0x0], x[0x0]=x[m_diis], x[m_diis]=_dp; } va_end(stack); return YERROR_NCONVERGED;
}






