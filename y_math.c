#include "y_math.h"
#include "y_system.h"

static  int32_t seed1, seed2;
static uint32_t qseed;


//----------- S I G N   M A N I P U L A T I O N   P A R T --------------//

inline double c_signum(char   value) { return   (char)( (0  < value) - (value < 0 ) ); }
inline double i_signum(int    value) { return    (int)( (0  < value) - (value < 0 ) ); }
inline double f_signum(float  value) { return  (float)( (0. < value) - (value < 0.) ); }
inline double d_signum(double value) { return (double)( (0. < value) - (value < 0.) ); }

//----------- L O W   P O W E R   R O U T I N E S --------------//

//This function calculate square of given value
inline double sqrd(register double value)
{
return value*value;
}

//This functionn calculate cube of give value
inline double cubed(register double value)
{
return value*value*value;
}

//This function get -1 power of the value
inline double calc_invert(register double value)
{
return (double)(1./value);
}

//This function calculates the inverse square root by Newton method.
//To edit!!!!!!!!!!!!!!!!!!!!!!!
inline double isqrt(register double x)
{
return 1./sqrt(x);
}

//This function calculates the dihedral angle from its sin and cos pair
inline double arc_angle(register double csA,register double snA)
{
if (csA<0.)
  {
  if (csA<=-1.) return -PI; 
  if (snA<0.)   return -acos(csA); //it's the 3-rd quater
  else          return +acos(csA); //it's the 2-nd quater
  }
else 
  {
  if (csA>=+1.) return +PI;
  if (snA<0.)   return -acos(csA); //it's the 4-th quater
  else          return +acos(csA); //it's the 1-st quater
  }
}


//This function define root of custom order
inline double get_root(register double x,register double order)
{
return exp(log(x)/order);
}

//This functions calculates lg and log2
inline float log2_float(float value)
{
return logf(value)/logf(2.00);
}
inline double log2_double(double value)
{
return log(value)/log(2.00);
}
inline float lg_float(float value)
{
return logf(value)/logf(10.00);
}
inline double lg_double(double value)
{
return log(value)/log(10.00);
}


//This function computes inverse of erf (using fromula from wiki: 0.5*sqrt(P)*{x+Px^3/12+7P^2x^5/480+127P^3x^7/40320+4369P^4x^9/5806080+34807P^5x^11/182476800+...} )
inline double ierf(register double x)
{
register double _d=PI*x*x;
return 0.5*sqrt(PI)*x*(1.+_d*(1./12.+_d*(7./480.+_d*(127./40320.+_d*(4369./5806080.+_d*(34807./182476800.+0.))))));
}


//----------- U N C E R T A I N E S S    P A R T --------------//


//This function calculates sin(alpha)/(alpha) around zero using equation sinT(T)=0.5+T*T/48.
inline double sinT(double theta2)
{
return .5+theta2/48.;
}


//----------- E Q U A T I O N S    P A R T -------------//

//This function solves square equations in real roots A*x^2+B*x+C=0. Amount of real roots is returned.
inline unsigned char solve_square_equation(double (*x)[2],double A,double B,double C)
{
double D;
//Switch equation type
if (!A) 
  {//It is not square equation actually
       if (!B) return 0;
  else if (!C) (*x)[0]=0.;    
  else         (*x)[0]=C/B;
  return 1; 
  }
else 
  {//Process exceptions
  if (!B)
    {
    if (!C) { (*x)[0]=0.; return 1; }
    else if (C<0.) return 0;
    else { (*x)[0]=sqrt(C); (*x)[1]=-(*x)[0]; return 2; }	
    }
  else if (!C) { (*x)[0]=0; (*x)[1]=-B/A; return 2; }
  else
    {//It is normal solver
    D=B*B-4.*A*C;
    if (fabs(D)<TINY) { (*x)[0]=-B/(2.*A); return 1; }
    if (D<0.) return 0;
    if (B>0.) D=-B/2.+sqrt(D)/2.;	//Awoid subtraction of huge values if b^2 is big
    else      D=-B/2.-sqrt(D)/2.;
	(*x)[0]=D/A, (*x)[1]=C/D;
    return 2;	
    }
  }
}
/*
//This function solves cubic equation for real roots A*x^3+B*x^2+C*x+D=0. Amount of real roots is returned.
inline unsigned char cubic_solve(double x[4],double A,double B,double C,double D)
{
}

//This function solves quartic equation for real roots A*x^4+B*x^3+C*x^2+D*x+E=0. Amount of real roots is returned.
inline unsigned char quartic_solve(double x[4],double A,double B,double C,double D,double E)
{
nroots=0;
double a, b, g, P, Q, R;
if (A) return cubic_solve((double[3])x,B,C,D,E);
B/=A, C/=A, D/=A, E/=A;
P=B*B, a=-3./8.*B+C, b=P*B/8.-B*C/2.+D, g=-3.*P*P/256.+C*P/16.-B*D/4.+E;
if (!B) 
  {
  if ((D=a*a-4.*g)<0.) return FALSE;
  else D=sqrt(D);
  if ((A=-a+D)>0.) { A=sqrt(A); x[0]=-B/4.+A, x[1]=-B/4.-A, nroots=2; }
  if ((A=-a-D)>0.) { A=sqrt(A); x[nroots++]=-B/4.+A, x[nroots++]=-B/4.-A; }
  return nroots;
  }
P=-a*a/12.-g, Q=-a*a/108.+a*g/3.-b*b/8.;
if (Q>0.) R=-Q/2.-sqrt(Q*Q/4.+P*P*P/27.);
else      R=-Q/2.+sqrt(Q*Q/4.+P*P*P/27.);
if (!R)                    P=-5./6.*a-get_cube_root(Q);
else { R=get_cube_root(R); P=-5./6.*a+R-P/(3.*R);       }
if ((Q=a+2.*P)<0.) return FALSE;
else Q=sqrt(Q);
if ((g=-(3.*a+2.*P+2.*b/Q))>0.)
  {
  g=sqrt(g);
  x[0]=-B/4.+(Q-g)/2., x[1]=B/4.+(Q+g)/2., nroots=2;
  }
if ((g=-(3.*a+2.*P-2.*b/Q))>0.)
  {
  g=sqrt(g);
  x[nroots++]=-B/4.-(Q-g)/2., x[nroots++]=B/4.-(Q+g)/2.;
  }
return nroots;
}
*/

//----------- R A N D O M   N U M B E R   G E N E R A T I O N ------------------//

//This routine init seed
void push_seed(int32_t _seed1,int32_t _seed2)
{
seed1=_seed1, seed2=_seed2;
}
//This function do pseudo-random number generation 
//Pierre L'ecuyer "Efficient and portable combined random number generators" ACM 1988 Vol 31 N6
double yrnd()
{
int32_t _k, Z;
_k=seed1/53668;
seed1=40041*(seed1-_k*53668)-_k*12211;
if (seed1<0) seed1+=2147483563;
_k=seed2/52774;
seed2=40692*(seed2-_k*52774)-_k*3791;
if (seed2<0) seed2+=2147483399;
Z=seed1-seed2;
if(Z<1) Z+=2147483562;
return (double)Z*4.656613E-10;
}

//Another, more primitive generator
//There are fast but robust random number generator
void push_qseed(uint32_t _qseed)
{
qseed=_qseed;
}
inline double qyrnd()
{
uint32_t temp;
qseed=1664525L*qseed+1013904223L;
temp=0x3f800000 | ( 0x007fffff & qseed );
return (double)((*(float*)&temp)-1.F);
}
//The thread-safe fast and robust random number generator
inline double pthread_qyrnd(uint32_t *qseed)
{
uint32_t temp;
(*qseed)=1664525L*(*qseed)+1013904223L;
temp=0x3f800000 | ( 0x007fffff & (*qseed) );
return (double)((*(float*)&temp)-1.F);
}

