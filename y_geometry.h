//This module contain geometry routine operatins

//Note. this module deals with 3D operation. If you are interested in somthing more complex check y_vec.c/y_vec.h
#define Y_GEOMETRY 0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_LIST
#include "y_list.h"
#endif
#ifndef Y_VECTOR
#include "y_vector.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif
#ifndef Y_LINTAR
#include "y_lintar.h"
#endif
#ifndef Y_MATRIX
#include "y_matrix.h"
#endif


#define MAX_GEOMETRY_VARIABLE_CHANGE 0.1

typedef struct{double x,y;}t_coord;


//--------------------------------------------------- P R I M I T I V E S ----------------------------------------------------------

//This function calculates distance between two points
//Note. This function returns square of distance
inline double calc_distance(t_vec *a,t_vec *b);
//The same as previous but n dimensional
double calc_ndistance(unsigned int n,double *a,double *b);

//This function calculates angle value [-pi...+pi] from its trigonometric functions
double calc_trig_angle(double csA,double snA);

//This function calculates cosine of angle between three points A-B-C
double calc_cos(t_vec *A,t_vec *B,t_vec *C);

//This function calculates dihedrals C-A-B-D value from four points
//Note this function returns cosine of angle instead of real angle value
double calc_dih_cos(t_vec *C,t_vec *A,t_vec *B,t_vec *D);

//This function calculates signed dihedral angle
void calc_dih_angle(double *csA,double *snA,t_vec *vec_i,t_vec *vec_j,t_vec *vec_k,t_vec *vec_l);

//This function converts the output of previous into real angle
double calc_dih_angle_value(t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl);

//--------------------------------------- D I S T A N C E     P A R T -------------------------------------------------------------

//This function solves the system of two euclidian distance equations:
// (x-x0)^2+(y-y0)^2=R0,
// (x-x1)^2+(y-y0)^2=R1.  it returns amount of found solutions pairs (in real numbers only!)
char solve_distances_equations_system_2D(double (*x)[2],double (*y)[2],double x0,double y0,double x1,double y1,double R0,double R1);

//-------------------------------------- 2 D    F I G U R E S    P A R T ----------------------------------------------------------

//This function calculates radii of a circle inscribed in triangle 
//        _______________________
//ir=S/p=V ((p-a)*(p-b)*(p-c))/p  
inline double calc_triangle_iradii(t_vec *a,t_vec *b,t_vec *c);

//This function calculates center of inscribed circle in triangle if the iradii is known
inline void calc_triangle_icenter(t_vec *x,double ir,t_vec *a,t_vec *b,t_vec *c);
 

//---------------------------------- T R A N S F O R M A T I O N   P A R T ----------------------------------------------------------

//This are translation functions A=B+S*_len, A and B might be the same
inline void translate_along_vector(t_vec *A,t_vec *B,register t_vec *s,register double _len);
inline void translate_along_vector_n(register unsigned int n,register double *a,register double *b,register double *s,register double _len);

/*
 * This function return the matrix for rotation around of given unit vector on angle alpha
 * For unit vector only:
 *     | csA+(1-csA)*i*i    | (1-csA)*j*i-snA*k  | (1-csA)*k*i+snA*j  |
 *     | (1-csA)*i*j+snA*k  | csA+(1-csA)*j*j    | (1-csA)*k*j-snA*i  |
 *     | (1-csA)*i*k-snA*j  | (1-csA)*j*k+snA*i  | csA+(1-csA)*k*k    |
 */
//Note. This function is probably worth than lower one as it probably require two square roots (for sin and vector lenght) and three addition multiplication for vector normalization
void rotate_around_uvector(t_tensor *R,t_vec *n,double csA,double snA);

/*
 * This function rotates system around given vector to angle alpha:
 *           | k*k+(i*i+j*j)*csA    k*i-k*i*csA-j*l*snA  k*j-k*j*csA+i*l*snA |
 *   1/L/L * | k*i-k*i*csA+j*l*snA  i*i+(k*k+j*j)*csA    i*j-i*j*csA-k*l*snA |
 *           | k*j-k*j*csA-j*l*snA  i*j-i*j*csA+k*l*snA  j*j+(k*k+i*i)*csA   |
 *
 *   L=sqrt(i*i+j*j+k*k)
 */
//Note. The snlA is  sinA/l that can be computed in a very efficient way as form of sqrt(sinA*sinA/_l2)=sqrt((1-cosA*cosA)/_l2). One square root is only required!
void rotate_around_vector(t_tensor *R,t_vec *n,double l2,double csA, double snlA);

//------------------------------------ E U L E R   A N G L E   P A R T --------------------------------------------

/*
 * This function fill euler angles rotational tensor
 *     | cs[a]cs[b]  cs[a]sn[b]sn[g]-sn[a]*cs[g] cs[a]sn[b]cs[g]+sn[a]sn[g] |
 * R = | sn[a]cs[b]  sn[a]sn[b]sn[g]+cs[a]*cs[g] sn[a]sn[b]cs[g]-cs[a]sn[g] |
 *     | -sn[b]      cs[b]sn[g]                  cs[b]cs[g]                 |
 */
void ueler(double csA,double snA,double csB,double snB,double csG,double snG,t_tensor *R);

/*
 * This function rotates system around axises in yulers angles:
 *    |  cos[b]cos[g]                     cos[b]sin[g]                    -sin[b]       |
 *  R=|  sin[a]sin[b]cos[g]-cos[a]sin[g]  sin[a]sin[b]sin[g]+cos[a]cos[g]  sin[a]cos[b] |
 *    |  cos[a]sin[b]cos[g]+sin[a]sin[g]  cos[a]sin[b]sin[g]-sin[a]cos[g]  cos[a]cos[b] |
 *
 *  A=R.B
 */
inline void calc_R_axis(double csA,double csB,double csG,double snA,double snB,double snG,t_tensor *R);
void rotate_yuler(t_vec *A,t_vec *B,double a,double b,double g);

//---------------------------------------     Q U A T E R N I O N S      P  A R T   -------------------------------------------

//This function calculates quaternion via exponential mapping of R^3 vector
inline void map_exp_uquaternion(t_quaternion *q,t_vec *v);

//This function calculates vector whose exponential mapping from R^3 results in quaternion
inline void imap_exp_uquaternion(t_vec *v,t_quaternion *q);

//This function calculates dq/dv
inline void dimap_exp_quaternion(t_qtensor *dq,t_vec *v);

//This function calculates rotation matrix from given quaternion
inline void calc_R_from_unit_quaternion(t_tensor *R,t_quaternion *q);

//This function calculates derivative of rotational tensor over quaternion components dR/dp
inline void calc_dR_dquaternion(t_tensor *dRi,t_tensor *dRj,t_tensor *dRk,t_qtensor *dq,t_quaternion *q);

//This function numerically perturbates a tensor so its determinant become Det(T)^2-1. < +/-tol
void round_udet_tensor(t_tensor *R,double tol);

//-------------------------------------------- D E R I V A T I V E S   P A R T -------------------------------------------------------

//This functin calculates first derivative of two interaction points in three dimensions (normal bonds)
inline void calc_bond_derivative(double *r,t_vec *ri,t_vec *rj,double *dr);

//This functin calculates first derivatives of three interaction points in three dimensions (normal angles)
inline void calc_angle_derivative(double *csA,double *snA,t_vec *ri,t_vec *rj,t_vec *rk,double *dr);

//This functin calculates first derivatives of four interaction points in three dimensions (dihedral angles)
inline void calc_dih_derivative(double *csA,double *snA,t_vec  *ri,t_vec *rj,t_vec *rk,t_vec *rl,double *dr);

//This functin calculates two first derivatives of two interaction points in three dimensions (normal bonds)
inline void calc_bond_derivatives(double *r,t_vec *ri,t_vec *rj,double *dr,double *ddr);

//This functin calculates two first derivatives of three interaction points in three dimensions (normal angles)
inline void calc_angle_derivatives(double *csA,double *snA,t_vec *ri,t_vec *rj,t_vec *rk,double *dr,double *ddr);

//This functin calculates two first derivatives of four interaction points in three dimensions (dihedral angles)
inline void calc_dih_angle_derivatives(double *csA,double *snA,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl,double *dr,double *ddr);

//------------------------------------- 3 D   P A R T --------------------------------------------------------//

//Taken from http://atheist4.narod.ru/mw/distance.htm the equations obtained as making two pairs vectors to be orthogonal - line and distance.
//This function calculates scale factors for both vectors of  minimal distance between two lines in 3D from given 2 points and 2 vectors
//NOTE. It returns  FALSE if lines are parallel (<SMALL2). 
char _calc_line_line_scale_distance_3D(double *t1,double *t2,t_vec *X1,t_vec *n1,t_vec *X2,t_vec *n2);
//This function calculates the distance between two lines in 3D from given 4 points
//NOTE. It returns distance and coordinates of orthogonal points that are NAN if lines are parallel. In addition it sets ylib_errno=YERROR_LEGAL if lines are parallel.
double calc_line_line_distance_3D(t_vec *M,t_vec *N,t_vec *A,t_vec *B,t_vec *C,t_vec *D);


//------------------------------------- S o S    P A R T -------------------------------------------

//The symbolic math part. See Herbert Edelsbrunner and Ernst Peter Mucke Simulation of Somplicity: a technique to cope with degenerate cases in geometrical algorithms.

//This function resorts size of indexes stored in stack as previous but do not update stack of sorting. For manual call only!!!!

inline unsigned int bjsort2(void *p0,void *p1,void **_p0,void **_p1); //1-1 comparisons
inline unsigned int bjsort3(void *p0,void *p1,void *p2,void **_p0,void **_p1,void **_p2); // 2-3 comparisons
inline unsigned int bjsort4(void *p0,void *p1,void *p2,void *p3,void **_p0,void **_p1,void **_p2,void **_p3); // 3-6 comparisons
inline unsigned int bjsort5(void *p0,void *p1,void *p2,void *p3,void *p4,void **_p0,void **_p1,void **_p2,void **_p3,void **_p4); // 7-16 comparisons
//isort function are the same as jsort but they replaces the address 
inline unsigned int isort2(void **_p0,void **_p1);
inline unsigned int isort3(void **_p0,void **_p1,void **_p2);
inline unsigned int isort4(void **_p0,void **_p1,void **_p2,void **_p3);
inline unsigned int isort5(void **_p0,void **_p1,void **_p2,void **_p3,void **_p4);

//This function resorts size of indexes stored in stack. For manual call only!!!!
unsigned int isort(unsigned int size, ... );
//This function resorts size of indexes stored in stack as previous but do not update stack of sorting. For manual call only!!!!
unsigned int jsort(unsigned int size, ...);

//This function performs 5d matix determinant sign calculations
char sign_ldet4(double **p);

//This function performs 4d matix determinant sign calculations
char sign_det4(double **p);

//This function performs 4d matix determinant sign calculations
char sign_ldet5(double **p);

//----------------------------- 3 D   F I G U R E S    P A R T ----------------------------------------------------//

//This function calculates the rmsd of two structures as difference in summ of corresponding atoms distances (sqares)
double calc_str_str_rmsd(unsigned int size,t_vec *A,t_vec *B);

//This function finds the coordinates of rectlinear tetrahedon that contain all points of given set in it. Function returns the edge lenght of tetrahedron.
double inscribe_tetrahedron(double scale,t_lvec *p0,t_lvec *p1,t_lvec *p2,t_lvec *p3,unsigned int size,t_lvec *lvec);

//This function performs insphere test for 3D case: Is point 'p' located inside the sphere spaned by points 'a', 'b', 'c' and 'd'
char in_sphere(t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d,t_lvec *p);

//This function is like above one but you can save the de4 abcd computations if it is known
char in_sphere_d4(char d4,t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d,t_lvec *p);

//This function performs intetrahedron test for 3D case: Is point p inside tetrahedron spaned by points 'a','b','c' and 'd'
char in_tetrahedron(t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d,t_lvec *p);




