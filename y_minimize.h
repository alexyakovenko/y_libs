#define Y_MINIMIZE 0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_LIST
#include "y_list.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif
#ifndef Y_VECTOR
#include "y_vector.h"
#endif
#ifndef Y_MATRIX
#include "y_matrix.h"
#endif
#ifndef Y_INTERPOLATE
#include "y_interpolate.h"
#endif

#define LINE_SEARCH_C1 0.0001
#define LINE_SEARCH_C2 0.8
#define GOLDEN_SECTION 0.381966
#define GOLDEN_RATIO   1.618034


/******************************** G R A D I E N T A R     M E T H O D S ******************************************/

/********************************* L I N E    S E A R C H *********************************************/

//Line-search method with square_approximation but without curvatures
//NOTE. *lambda is the lenght of the previous step, *ff is initial function value
unsigned int line_search_square_fapproximation(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int n,double *x0,double *x1,double *g1,double **G,double *p,unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list stack);
//The same as previous but the search is done with functions evaluations and only the very last call (at return statements) re-calculates the gradient
//NB! The trigger of gradient calculation is to have pointers g!=0x0 and optionally G!=0
unsigned int line_search_square_fapproximation_only(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int n,double *x0,double *x1,double *g1,double **G,double *p,unsigned int  (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list stack);

//Doing 4-point schema: generating bad-good-bad situation and then appliing one step of Brent quadratic interpolation
//NB! The trigger of gradient calculation is to have pointers g!=0x0 and optionally G!=0
unsigned int line_qsearch_square_fapproximation_only(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int n,double *x0,double *x1,double *g1,double **G,double *p,unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list stack);
//The same as previous but it modifies only a subset of variables (from - to) keeping the rest frozed
unsigned int line_qssearch_square_fapproximation_only(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int from,unsigned int to,unsigned int n,double *x0,double *x1,double *g1,double **G,double *p,unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list stack);

//The same as above 'line_qsearch_square_fapproximation_only()' but with separate error vector (for routines like GDIIS)
unsigned int line_qsearch_e_square_fapproximation_only(double tol,double itol,double max_variable_change,double *lambda,double *ff,unsigned int n,double *x0,double *x1,double *g1,double *e1,double **G,double *p,unsigned int (*funct)(double *,unsigned int ,double *,double *,double *,double **,va_list ),va_list stack);



/*********************** P O L A K _ R I B I E R E   (  C o n j u g a t e    g r a d i e n t)    P A R T ************************/



//This function performs optimization with Polak-Ribiere (conjugate gradient) method:
// p(k)=-g(k)+b(k-1)*p(k-1),    b(k-1)=((gk-gk-1)^T.gk)/((gk-gk-1)^T.pk).
unsigned int polak_ribiere(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **g,double **G,double *p,
                           unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),
                           unsigned int (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double *,double *,double *,double **,double *,
                                                       unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... );



/************************************* L B F G S   P A R T *******************************************/



//This function (separately) updates l-bfgs structure: S, Y, IR and YY
unsigned int update_lbfgs(unsigned int n,double *x0,double *x1,double *g0,double *g1,unsigned int m,unsigned int size_m,double **S,double **Y,double **IR,double **YY);

//This function calculates H-1.g with l-BFGS data
inline void calc_Hg_lbfgs(unsigned int n,double *g,double *p,unsigned int size_m,double gamma,double **S,double **Y,double **IR,double **YY,double *c,double *temp);

//This function is reverse call of l-BFGS update routine
unsigned int rlbfgs(unsigned int n,double *x0,double *x1,double *g0,double *g1,double *p,unsigned int *m,unsigned int size_m,double gamma,double **S,double **Y,double **IR,double **YY,double *c,double *temp);

//This function performs unconstrained l-bfgs minimization (see Richard H. Byrd, Jorgr Nocedal and Robert B. Schnabel "Representation of quasi-newton matrices and their use in limited memory methods" )
unsigned int lbfgs_new(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **g,double **G,double *p,
                  unsigned int size_m,double **S,double **Y,double **R,double **YY,double *c,double *temp,
                  unsigned int funct(double *fx,unsigned int n,double *x,double *g,double **G,va_list stack),
                  unsigned int (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double *,double *,double *,double **,double *,
                                              unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... );



//DEPRECATED!!!
//This function performs unconstrained l-bfgs minimization (see Richard H. Byrd, Jorgr Nocedal and Robert B. Schnabel "Representation of quasi-newton matrices and their use in limited memory methods" )
// n - number of variables, m - number of corrections (3...7 recomended), x - two sets of variables-coordinates, g - two sets of gradient, funct - function and its derivative, tol - tolerance of dfunct
//Note. The S and Y matrices are trasposed to use processors cache more effectively
char lbfgs(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **g,double **G,double *z,unsigned int m,char (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),
           char (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double **,double **,double **,double *,char (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... );




/*************************************   D I I S   P A R T   *******************************************/



//This function (separately) updates l-bfgs structure: S, Y, IR and YY
void update_diis(unsigned int n,unsigned int m,unsigned int size_e,double **x,double **e,double **B,double *val);

//The reversibly-calling diis - it computes DIIS interpolation weigths for given set of error vectors e
unsigned int rdiis(unsigned int n,unsigned int *m,unsigned int size_e,double **x,double **e,double **B,double *val,double *c);

//This routine create diis vector
inline void _synthesize_diis_vector(unsigned int n,unsigned int size_e,double *vi,double **v,double *c);

//This function minimizes function using Direct Inversion in the Iterative (sub-)Space method [DIIS]
unsigned int diis(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **G,
                  unsigned int size_e,double **e,double **B,double *c,double *val,
                  unsigned int funct(double *fx,unsigned int n,double *x,double *e,double **G,va_list stack),
                  unsigned int (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double *,double *,double *,double **,double *,
                                              unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... );

//The reversibly-calling gdiis - it computes DIIS interpolation weigths for given set of error vectors e
unsigned int rgdiis(unsigned int n,unsigned int *m,unsigned int size_e,double **x,double **g,double **e,double **B,double *val,double *c);

//This routine minimizes function using GDIIS approach
unsigned int gdiis(double *fx,unsigned int nsteps,double tol,double itol,double max_variable_change,unsigned int n,double **x,double **g,double **G,double *p,
                   unsigned int size_e,double **e,double **B,double *c_diis,double *val,
                   unsigned int size_m,double **S,double **Y,double **IR,double **YY,double *c_lbfgs,double *temp_lbfgs,
                   unsigned int funct(double *fx,unsigned int n,double *x,double *g,double **G,va_list stack),
                   unsigned int (*line_search)(double ,double ,double ,double *,double *,unsigned int ,double *,double *,double *,double **,double *,
                                               unsigned int (*funct)(double *,unsigned int ,double *,double *,double **,va_list ),va_list ),char *label, ... );





