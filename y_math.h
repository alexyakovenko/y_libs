#define Y_MATH 0x1

#ifndef Y_SYSTEM
#include "../y_libs/y_system.h"
#endif

//#define __USE_ISOC99 I_want_ISO
#include <math.h>

#define CALC_INVERT(A) (double)(1.00/((double)A))

#define CRON(i,j) (i==j) ? 1 : 0
#define ETTA(i,j) (i!=j) ? 1 : 0

#define SQRT_2  1.4142135623730950488
#define SQRT_3  1.7320508075688772935
#define SQRT_10 3.162277660168379332
#define LN_2    0.6931471805599453094
#define LN_3    1.0986122886681096914
#define LN_10   2.302585092994045684
#define PI      3.141592653589793238462643383279502884197
#define SQRT_PI 1.7724538509055160273
#define EXP_E   2.7182818284590452353603
#define EXP_G  10.0000000000000000000
//Massive of integger powers of exponent

/* "small number" to avoid overflow in some cases */
#define EPSILON   1.e-72       //The machine precission
#define TINY      1.e-18       //Reasonable precission for mathematical algorithms to work 
#define SQRT_TINY 1.e-8        //Qsuare root of TINY
#define SMALL     1.e-3        //Resonable accuracy for molecular modeling geometry purposes 
#define SMALL2    1.e-6        //For sure resonable accuracy for molecular modeling geometry purposes


//----------- S I G N   M A N I P U L A T I O N   P A R T --------------//

inline double c_signum(char   value);
inline double i_signum(int    value);
inline double f_signum(float  value);
inline double d_signum(double value);

//----------- L O W   P O W E R   R O U T I N E S --------------//

//This function calculate square of given value
inline double sqrd(register double value);

//This functionn calculate cube of give value
inline double cubed(register double value);

//This function get -1 power of the value
inline double calc_invert(register double value);

//This function define root from custom order
inline double get_root(register double order,register double x);

//This function calculates the inverse square root by Newton method
inline double isqrt(register double x);

//This function calculates the dihedral angle from its sin and cos pair
inline double arc_angle(register double csA,register double snA);

//This functions calculates log2 and lg
inline float log2_float(float value);
inline double log2_double(double value);
inline float lg_float(float value);
inline double lg_double(double value);

//This function computes inverse of erf (using fromula from wiki: 0.5*sqrt(P)*{x+Px^3/12+7P^2x^5/480+127P^3x^7/40320+4369P^4x^9/5806080+34807P^5x^11/182476800+...} )
inline double ierf(register double x);

//----------- U N C E R T A I N E S S    P A R T --------------//

//This function calculates sin(alpha)/(alpha) around zero using equation sinT(T)=0.5+T*T/48.
inline double sinT(double theta2);

//----------- E Q U A T I O N S    P A R T -------------//

//This function solves square equations in real roots a*x^2+b*x+c=0
inline unsigned char solve_square_equation(double (*x)[2],double a,double b,double c);

//----------- R A N D O M   N U M B E R   G E N E R A T I O N ------------------//

//Pierre L'ecuyer "Efficient and portable combined random number generators" ACM 1988 Vol 31 N6
//This routine init seed
void push_seed(int32_t _seed1, int32_t _seed2);
//This function do pseudo-random number generation 
double yrnd();

//Another, more primitive generator
//There are fast but robust random number generator
void push_qseed(uint32_t _qseed);
inline double qyrnd();
//The thread-safe fast and robust random number generator
inline double pthread_qyrnd(uint32_t *qseed);


//------------------------------------- S o S    P A R T -------------------------------------------

//The symbolic math part. See Herbert Edelsbrunner and Ernst Peter Mucke Simulation of Somplicity: a technique to cope with degenerate cases in geometrical algorithms.

//This function resorts size of indexes stored in stack as previous but do not update stack of sorting. For manual call only!!!!

//inline unsigned int bjsort2(void *p0,void *p1,void **_p0,void **_p1); //1-1 comparisons
//inline unsigned int bjsort3(void *p0,void *p1,void *p2,void **_p0,void **_p1,void **_p2); // 2-3 comparisons
//inline unsigned int bjsort4(void *p0,void *p1,void *p2,void *p3,void **_p0,void **_p1,void **_p2,void **_p3); // 3-6 comparisons
//inline unsigned int bjsort5(void *p0,void *p1,void *p2,void *p3,void *p4,void **_p0,void **_p1,void **_p2,void **_p3,void **_p4); // 7-16 comparisons
//isort function are the same as jsort but they replaces the address 
//inline unsigned int isort2(void **_p0,void **_p1);
//inline unsigned int isort3(void **_p0,void **_p1,void **_p2);
//inline unsigned int isort4(void **_p0,void **_p1,void **_p2,void **_p3);
//inline unsigned int isort5(void **_p0,void **_p1,void **_p2,void **_p3,void **_p4);

//This function performs 5d matix determinant sign calculations
char sign_ldet4(double **p);

//This function performs 4d matix determinant sign calculations
char sign_ldet5(double **p);

//This function performs 4d matix determinant sign calculations
char sign_det4(double **p);


