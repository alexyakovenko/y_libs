#define Y_INTERPOLATE 1

#ifndef Y_SYSTEM
#include "y_system.h"
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
#ifndef Y_GEOMETRY
#include "y_geometry.h"
#endif
#ifndef Y_MINIMIZE
#include "y_minimize.h"
#endif


//Inverse quadratic interpolation
inline char inverse_quadratic_interpolation(register double *x,register double xa,register double xb,register double xc,register double fa,register double fb,register double fc);

//------------------------- M O N O T O N E   C U B I C   I N T E R P O L A T I O N --------------------------------

//The interpolation polynome derivative model on the left
inline double _calc_mc_interpolation_derivative_left(register double x[4],register double f[4]);
//The interpolation polynome derivative model
inline double _calc_mc_interpolation_derivative(register double s0,register double s1,register double dx0,register double dx1); 
//The interpolation polynome derivative model on the right
inline double _calc_mc_interpolation_derivative_right(register double x[4],register double f[4]);
//The interpolation polynome
inline void _calc_mc_interpolation_polynome(double c[4],register double f0,register double g0,register double g1,register double s,register double dx0);
//The interpolation itself
inline double _calc_mc_interpolation(register double c[4],register double x0,register double x);
			  
//The calculator. It generates aligned massive of values for the common t0+n*dt checkpoints using arbitrary t,v sets of datapoints.
//To get missing values, monotonic three-cubic interpolation is used (gl is the 'left monotonic' gradient).
//NOTE. The uncovered area in w is marked as NAN;
char cubic_align_dataset_in_time(register double t0,register double dt,register unsigned int n,register double *w,register unsigned int size,register double *t,register double *v,double *g);

//--------------------- I N T E R P O L A T I O N     O N     G R I D -----------------------------------------


typedef struct{
              t_len len;
              t_vec ori;
              double sp;
              char ***c;              
              }t_cgrid;           //three-dimensional char grid

typedef struct{
              t_len len;
              t_vec ori;
              double sp;
              double ***i;              
              }t_igrid;           //three-dimensional int grid

typedef struct{
              t_len len;
              t_vec ori;
              double sp;
              double ***d;              
              }t_dgrid;           //three-dimensional double grid

typedef struct{
              t_len len;
              t_vec ori;
              double sp;
              double (***d)[64];
              }t_tcgrid;          //three-cubic grid

//This function creates grid for 3D char grid
t_cgrid* alloc_cgrid(unsigned int ni,unsigned int nj,unsigned int nk);
//This function creates grid for 3D char grid
t_igrid* alloc_igrid(unsigned int ni,unsigned int nj,unsigned int nk);
//This function creates grid for 3D double interpolation
t_dgrid* alloc_dgrid(unsigned int ni,unsigned int nj,unsigned int nk);
//This function creates grid fot tricubic interpolation
t_tcgrid* alloc_tcgrid(unsigned int ni,unsigned int nj,unsigned int nk);

//This function read 3D double grid file
t_cgrid *read_cgrid(FILE *in);
//This function read 3D double grid file
t_igrid *read_igrid(FILE *in);
//This function read 3D double grid file
t_dgrid *read_dgrid(FILE *in);
//This function read grid file
t_tcgrid *read_tcgrid(FILE *in);

//This function writes 3D char grid to file
char write_cgrid(FILE *out,t_cgrid *cgrid);
//This function writes 3D int grid to file
char write_igrid(FILE *out,t_igrid *igrid);
//This function writes 3D double grid to file
char write_dgrid(FILE *out,t_dgrid *dgrid);
//This function writes tricubic grid to file
char write_tcgrid(FILE *out,t_tcgrid *tcgrid);

//   T H R E E - L I N E A R     I N T E R P O L A T I O N     P A R T 

//--------------------- T R I L I N E A R     P A R T ------------------------------//

//This function calculates values and its derivaives over grid
inline double calc_trilinear_interpolation(double i, double j, double k,
                                           double f000, double f100, double f010, double f001,
                                           double f110, double f101, double f011, double f111);
//This function calculates value of trilinear interpolation on grid
char calc_linear_interpolation_value_dfindif(double *f,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map);
//This function calculates values and its derivaives over grid
inline double calc_dtrilinear_interpolation(double i, double j, double k,
                                            double f000, double f100, double f010, double f001,
                                            double f110, double f101, double f011, double f111,t_vec *df);
//This function calculates gradients of trilinear interpolation on grid
char calc_linear_interpolation_grad_dfindif(double *f,t_vec *g,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map);

//   T H R E E - C U B I C     I N T E R P O L A T I O N     P A R T 
// Lekien F. and Marsden J. Tricubic interpolation in three dimensions Int.J.Numer.Meth.Engng 2005, 63:455-471

//This function calculates threecubic proto-map
char calc_threecubic_protomap(unsigned int ni,unsigned int nj,unsigned int nk,t_vec *origin,double spacering,double (***p_map)[8],char (*funct)(double *f,double *dfdx,double *dfdy,double *dfdz,double *d2fdxdy,double *d2fdxdz,double *d2fdydz,double *d3fdxdydz,double rr,t_vec *r,char label,va_list stack),char label, ... );

//This function transforms interpolation proto-map [nx][ny][nz][8] of {fx,dx,dy,dz,dxdy,dydz,dzdx,dxdydz} into threecubic cooficients [nx-1][ny-1][nz-1][64]    
void transform_threecubic_protomap(unsigned int ni,unsigned int nj,unsigned int nk,double (***a_map)[64],double (***p_map)[8]);

//This function calculates threecubic iterpolated value of function
char calc_threecubic_interpolation_value(double *f,t_vec *origin,double spacering,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double (***a_map)[64]);
//This function calculates threecubic iterpolated value of function and its first derivative
char calc_threecubic_interpolation_derivative(double *f,t_vec *g,t_vec *origin,double spacering,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double (***a_map)[64]);
//This function calculates threecubic iterpolated value of function, its first and second derivatives
char calc_threecubic_interpolation_derivatives(double *f,t_vec *g,t_tensor *G,t_vec *origin,double spacering,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double (***a_map)[64]);

//--------- N O N D E R I V A T I V E    I N T E R P O L A T I O N    P A R T ------------------//

//This function calculates tricubic interpolation value using finite differences instead of derivatives
char calc_tricubic_interpolation_value_dfindif(double *f,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map);
//This function calculates tricubic interpolation value and its derivative using finite differences instead of derivatives
char calc_tricubic_interpolation_derivative_dfindif(double *f,t_vec *g,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map);

//Hyman J.M. Accurate monotonicity preserving cubic interpolation SIAM J.Sci.Stat.Comput., 4, 645-654, 1983
//Fritsch F.N. and Carlson R.E. Monotonicity preserving bicubic interpolation: A progress report Comp.Aid.Geom.Des., 2, 117-121, 1985 
//Carlson R.E. and Fritsch F.N. Monotone piecewise bicubic interpolation SIAM J.Numer.Anal., 22, 386-400, 1985 

//This function calculates tricubic polinome values with satisfing monotonicity inequalities
char calc_tricubic_interpolation_wp_monotonicity(double *f,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map);
//This function calculates tricubic polinome values with satisfing monotonicity inequalities
char calc_tricubic_interpolation_derivative_wp_monotonicity(double *f,t_vec *g,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map);


//------------------------------------------------------ O B J E C T S     R O U N D I N G -------------------------------------------------//

//This function numerically perturbates a tensor so its determinant become |Det(T)^2-1.| < tol
void round_udet_tensor(t_tensor *R,double tol);



