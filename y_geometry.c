#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include "y_geometry.h"
//This module contain geometry routine operatins

//--------------------------------------------------- P R I M I T I V E S ----------------------------------------------------------

//This function calculates distance between two points
//Note. This function returns square of distance
inline double calc_distance(t_vec *a,t_vec *b)
{
return sqrd(a->i-b->i)+sqrd(a->j-b->j)+sqrd(a->k-b->k);
}
//The same as previous but n dimensional
double calc_ndistance(unsigned int n,double *a,double *b)
{
register double _d=0.00;
while(n--)
  _d+=sqrd(a[n]-b[n]);
return _d;
}

//This function calculates angle value [-pi...+pi] from its trigonometric functions
double calc_trig_angle(double csA,double snA)
{
return (snA>0.) ? +acos(csA) : -acos(csA);
}

//This function calculates cosine of angle between three points A-B-C 
double calc_cos(t_vec *A,t_vec *B,t_vec *C)
{
register double _ab,_bc,_ca,_cos_abc;
_ab=calc_distance(A,B);
_bc=calc_distance(B,C);
_ca=calc_distance(C,A);
_cos_abc=(_ab+_bc-_ca)/(2.00*sqrt(_ab*_bc));
     if (_cos_abc>1.00)
       _cos_abc=1.00;
else if (_cos_abc<-1.00)
       _cos_abc=-1.00;
return _cos_abc;
}

//This function calculates dihedrals C-A-B-D value from four points
//Note this function returns cosine of angle instead of real angle value
double calc_dih_cos(t_vec *C,t_vec *A,t_vec *B,t_vec *D)
{
t_vec r[3],p[2];
register double _d;
r[0].i=B->i-A->i, r[0].j=B->j-A->j, r[0].k=B->k-A->k;  // 3 flops
r[1].i=C->i-A->i, r[1].j=C->j-A->j, r[1].k=C->k-A->k;  // 3 flops    
r[2].i=D->i-B->i, r[2].j=D->j-B->j, r[2].k=D->k-B->k;  // 3 flops 
vec_vec_vmult(&p[0],&r[1],&r[0]);                // 9 flops
vec_vec_vmult(&p[1],&r[2],&r[0]);                // 9 flops
_d=calc_vec_vec_scalar_product(&p[0],&p[1])/sqrt(calc_vec_norm(&p[0])*calc_vec_norm(&p[1])); // 17 flops + 1 sqrt
     if (_d>+1.00) _d=+1.00;
else if (_d<-1.00) _d=-1.00;
return _d;
//total 44 flps + 1 root
}

//This function calculates signed dihedral angle
void calc_dih_angle(double *csA,double *snA,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl)
{
t_vec f, g, h, a, b;
double aa, bb, sqrtgg;
f.i=ri->i-rj->i, f.j=ri->j-rj->j, f.k=ri->k-rj->k; // F                        3 flops
g.i=rj->i-rk->i, g.j=rj->j-rk->j, g.k=rj->k-rk->k; // G                        3 flops
h.i=rl->i-rk->i, h.j=rl->j-rk->j, h.k=rl->k-rk->k; // H                        3 flops
vec_vec_vmult(&a,&f,&g);                           // A                        9 flops
vec_vec_vmult(&b,&h,&g);                           // B                        9 flops
//Calculate cos(alpha) and sin(alpha)
if ((aa=calc_vec_norm(&a))<SMALL2) aa=1./SMALL2; else aa=1./aa;  //              6 flops      
if ((bb=calc_vec_norm(&b))<SMALL2) bb=1./SMALL2; else bb=1./bb;  //              6 flops
sqrtgg=sqrt(aa*bb);                                //                            1 flops + 1 root
*csA=calc_vec_vec_scalar_product(&a,&b)*sqrtgg;    //                            6 flops
*snA=calc_vec_vec_scalar_product(&a,&h)*sqrtgg;    //                            6 flops
sqrtgg=sqrt(calc_vec_norm(&g));                    //                            5 flops + 1 root
*snA*=sqrtgg;                                      //                            1 flops
     if (*csA>+1.) *csA=+1.;
else if (*csA<-1.) *csA=-1.;
else if (*snA>+1.) *snA=+1.;
else if (*snA<-1.) *snA=-1.;
// TOTAL: 60 flops + 2 sqrt
}

//This function converts the output of previous into real angle
double calc_dih_angle_value(t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl)
{
double csF, snF;
calc_dih_angle(&csF,&snF,ri,rj,rk,rl);
return calc_trig_angle(csF,snF);
}

//--------------------------------------- D I S T A N C E     P A R T -------------------------------------------------------------

//This function solves the system of two euclidian distance equations:
// (x-x0)^2+(y-y0)^2=R0,
// (x-x1)^2+(y-y0)^2=R1.  it returns amount of found solutions pairs (in real numbers only!)
char solve_distances_equations_system_2D(double (*x)[2],double (*y)[2],double x0,double y0,double x1,double y1,double R0,double R1)
{
double _d;
char exchange;

//Step 1. Chose leading variable and check exceptions.
(*x)[0]=x1-x0, (*y)[0]=y1-y0;
if (fabs((*x)[0])>fabs((*y)[0]))
  { // y - is leading variable
  if (fabs((*x)[0])<TINY) return 0; //Equations are indeterminable due to singularity
  exchange=TRUE;
  _d=(*x)[0], (*x)[0]=(*y)[0], (*y)[0]=_d, _d=x0, x0=y0, y0=_d, _d=x1, x1=y1, y1=_d;
  }
else 
  {// x -is leading variable
  if (fabs((*y)[0])<TINY) return 0; //Equations are indeterminable due to singularity
  exchange=FALSE;
  }
//Step 2. Calculate y[1]=b and y[0]=k;
(*y)[1]=-(*x)[0]/(*y)[0];                            //k
(*y)[0]=-.5*(R1-x1*x1-y1*y1-R0+x0*x0+y0*y0)/(*y)[0]; //b
//Step 3. Solve square equation
switch (solve_square_equation(x,1.+(*y)[1]*(*y)[1],2.*((*y)[1]*(*y)[0]-y0*(*y)[1]-x0),x0*x0+y0*y0-R0+(*y)[0]*(*y)[0]-2.*y0*(*y)[0]))
  {
  case 1  : { //Only one root exists
            (*y)[0]+=(*y)[1]*(*x)[0];
            if (exchange)
              {//Restore correct order of variables
              _d=(*y)[0], (*y)[0]=(*x)[0], (*x)[0]=_d;
              } 
            return 1; }
  case 2  : { //Two roots exists 
            (*y)[0]+=(*y)[1]*(*x)[0], (*y)[1]=(*y)[0]+(*y)[1]*((*x)[1]-(*x)[0]);
            if (exchange)
              {//Restore correct order of variables
              _d=(*y)[0], (*y)[0]=(*x)[0], (*x)[0]=_d;
              _d=(*y)[1], (*y)[1]=(*x)[1], (*x)[1]=_d;
              } 
            return 2; }
  default : return 0; //No solution exists in real numbers
  }
}


//-------------------------------------- 2 D    F I G U R E S    P A R T ----------------------------------------------------------

//This function calculates radii of a circle inscribed in triangle 
//        _______________________
//ir=S/p=V ((p-a)*(p-b)*(p-c))/p  
inline double calc_triangle_iradii(t_vec *a,t_vec *b,t_vec *c)
{
register double _a, _b, _c, _p;
_a=sqrt(calc_distance(b,c)), _b=sqrt(calc_distance(a,c)), _c=sqrt(calc_distance(a,b)), _p=(_a+_b+_c)/2.;
return sqrt(((_p-_a)*(_p-_b)*(_p-_c))/_p);
}

//This function calculates center of inscribed circle in triangle if the iradii is known
inline void calc_triangle_icenter(t_vec *x,double ir,t_vec *a,t_vec *b,t_vec *c)
{
t_vec _ab, _ac;

_ab.i=b->i-a->i, _ab.j=b->j-a->j, _ab.k=b->k-a->k, multiple_vec_scalar(&_ab,&_ab,1./sqrt(calc_vec_norm(&_ab)));
_ac.i=c->i-a->i, _ac.j=c->j-a->j, _ac.k=c->k-a->k, multiple_vec_scalar(&_ac,&_ac,1./sqrt(calc_vec_norm(&_ac)));
_ac.i=(_ab.i+_ac.i)/2., _ac.j=(_ab.j+_ac.j)/2., _ac.k=(_ab.k+_ac.k)/2.; //_ac is a bisect now
multiple_vec_scalar(&_ac,&_ac,ir/sqrt((1.-sqrd(calc_vec_vec_scalar_product(&_ab,&_ac))/calc_vec_norm(&_ac))*calc_vec_norm(&_ac)));
x->i=a->i+_ac.i, x->j=a->j+_ac.j, x->k=a->k+_ac.k;
}


//---------------------------------- T R A N S F O R M A T I O N   P A R T ----------------------------------------------------------

//This are translation functions A=B+s*len
inline void translate_along_vector(t_vec *A,t_vec *B,register t_vec *s,double _len)
{
if (A==B) { A->i+=s->i*_len,     A->j+=s->j*_len,     A->k+=s->k*_len;     }
else      { A->i=B->i+s->i*_len, A->j=B->j+s->j*_len, A->k=B->k+s->k*_len; }
}
inline void translate_along_vector_n(register unsigned int n,register double *a,register double *b,register double *s,register double _len)
{
if (a==b) while (n--) { *a+=_len**s,   a++,      s++; }
else      while (n--) { *a=*b+_len**s, a++, b++, s++; } 
}

/*
 * This function return the matrix for rotation around of given unit vector on angle alpha
 * For unit vector only:
 *     | csA+(1-csA)*i*i    | (1-csA)*j*i+snA*k  | (1-csA)*k*i-snA*j  |
 *     | (1-csA)*i*j-snA*k  | csA+(1-csA)*j*j    | (1-csA)*k*j+snA*i  |
 *     | (1-csA)*i*k+snA*j  | (1-csA)*j*k-snA*i  | csA+(1-csA)*k*k    |
 */
//Note. This function is probably worth than lower one as it probably require two square roots (for sin and vector lenght) and three addition multiplication for vector normalization
inline void rotate_around_uvector(t_tensor *R,t_vec *n,double csA,double snA)
{
register double ucsA=1.00-csA;           // 1  flops

(*R)[0][0]=csA+ucsA*n->i*n->i;          // 3  flops
(*R)[1][1]=csA+ucsA*n->j*n->j;          // 3  flops
(*R)[2][2]=csA+ucsA*n->k*n->k;          // 3  flops
(*R)[1][0]=(*R)[0][1]=ucsA*n->j*n->i;   // 3  flops
(*R)[2][0]=(*R)[0][2]=ucsA*n->k*n->i;   // 3  flops
(*R)[2][1]=(*R)[1][2]=ucsA*n->j*n->k;   // 3  flops
(*R)[0][1]+=snA*n->k;                   // 2  flops
(*R)[0][2]-=snA*n->j;                   // 2  flops
(*R)[1][0]-=snA*n->k;                   // 2  flops
(*R)[1][2]+=snA*n->i;                   // 2  flops
(*R)[2][0]+=snA*n->j;                   // 2  flops
(*R)[2][1]-=snA*n->i;                   // 2  flops
                                        // 31 flops total
}

/*
 * This function rotates system around given vector to angle alpha:
 *           | k*k+(i*i+j*j)*csA    k*i-k*i*csA-j*l*snA  k*j-k*j*csA+i*l*snA |
 *   1/L/L * | k*i-k*i*csA+j*l*snA  i*i+(k*k+j*j)*csA    i*j-i*j*csA-k*l*snA |
 *           | k*j-k*j*csA-j*l*snA  i*j-i*j*csA+k*l*snA  j*j+(k*k+i*i)*csA   |
 *
 *   L=sqrt(i*i+j*j+k*k)
 */
//Note. The snlA=sinA/l that can be computed in a very efficient way as snlA=(+/-)sqrt(sinA*sinA/_l2)=(+/-)sqrt((1-cosA*cosA)/_l2). One square root is only required!
inline void rotate_around_vector(t_tensor *R,t_vec *n,double l2,double csA, double snlA)
{
double ucsA=(1.00-csA)/l2;              // 2  flops

(*R)[0][0]=csA+ucsA*n->i*n->i;          // 3  flops
(*R)[1][1]=csA+ucsA*n->j*n->j;          // 3  flops
(*R)[2][2]=csA+ucsA*n->k*n->k;          // 3  flops
(*R)[1][0]=(*R)[0][1]=ucsA*n->k*n->i;   // 3  flops
(*R)[2][0]=(*R)[0][2]=ucsA*n->j*n->k;   // 3  flops
(*R)[2][1]=(*R)[1][2]=ucsA*n->i*n->j;   // 3  flops
(*R)[0][1]-=snlA*n->j;                  // 2  flops
(*R)[0][2]+=snlA*n->i;                  // 2  flops
(*R)[1][0]+=snlA*n->j;                  // 2  flops
(*R)[1][2]-=snlA*n->k;                  // 2  flops
(*R)[2][0]-=snlA*n->j;                  // 2  flops
(*R)[2][1]+=snlA*n->k;                  // 2  flops
                                        // 32 flops total
}

//--------------------------------  E U L E R    A N G L E S    P A R T -----------------------------------

/*
 * This function fill euler angles rotational tensor
 *     | cs[a]cs[b]  cs[a]sn[b]sn[g]-sn[a]*cs[g] cs[a]sn[b]cs[g]+sn[a]sn[g] |
 * R = | sn[a]cs[b]  sn[a]sn[b]sn[g]+cs[a]*cs[g] sn[a]sn[b]cs[g]-cs[a]sn[g] |
 *     | -sn[b]      cs[b]sn[g]                  cs[b]cs[g]                 |
 */
void ueler(double csA,double snA,double csB,double snB,double csG,double snG,t_tensor *R)
{
(*R)[0][0]= csA*csB, (*R)[0][1]= csA*snB*snG-snA*csG, (*R)[0][2]= csA*snB*csG+snA*snG;
(*R)[1][0]= snA*csB, (*R)[1][1]= snA*snB*snG+csA*csG, (*R)[1][2]= snA*snB*csG-csA*snG;
(*R)[2][0]=-snB,     (*R)[2][1]= csB*snG,             (*R)[2][2]= csB*csG;
}




/*
 * This function rotates system around axises in yulers angles:
 *    |  cos[b]cos[g]                     cos[b]sin[g]                    -sin[b]       |
 *  R=|  sin[a]sin[b]cos[g]-cos[a]sin[g]  sin[a]sin[b]sin[g]+cos[a]cos[g]  sin[a]cos[b] |
 *    |  cos[a]sin[b]cos[g]+sin[a]sin[g]  cos[a]sin[b]sin[g]-sin[a]cos[g]  cos[a]cos[b] |
 *
 *  A=R.B
 */
inline void calc_R_axis(double csA,double csB,double csG,double snA,double snB,double snG,t_tensor *R)
{
(*R)[0][0]= csB*csG,             (*R)[0][1]= csB*snG,             (*R)[0][2]=-snB;
(*R)[1][0]= snA*snB*csG-csA*snG, (*R)[1][1]= snA*snB*snG+csA*csG, (*R)[1][2]= snA*csB;
(*R)[2][0]= csA*snB*csG+snA*snG, (*R)[2][1]= csA*snB*snG-snA*csG, (*R)[2][2]= csA*csB;
}
void rotate_yuler(t_vec *A,t_vec *B,double a,double b,double g)
{
double cs[3],sn[3];
t_tensor R;

cs[0]=cos(a),cs[1]=cos(b),cs[2]=cos(g);
sn[0]=sin(a),sn[1]=sin(b),sn[2]=sin(g);
calc_R_axis(cs[0],cs[1],cs[2],sn[0],sn[1],sn[2],&R);
multiple_origin_tensor_origin_vec((t_vec*)A,&R,(t_vec*)B);
}


//---------------------------------------     Q U A T E R N I O N S      P  A R T   -------------------------------------------

//This function calculates quaternion via exponential mapping of R^3 vector
inline void map_exp_uquaternion(t_quaternion *q,t_vec *v)
{
double vv, snT;
if ((vv=calc_vec_norm(v))<TINY) snT=sinT(vv);
else
  {
  vv=sqrt(vv);
  snT=sin(0.5*vv)/vv;
  }
q->i=snT*v->i;
q->j=snT*v->j;
q->k=snT*v->k;
q->w=cos(0.5*vv);
}

//This function calculates vector whose exponential mapping from R^3 results in quaternion
inline void imap_exp_uquaternion(t_vec *v,t_quaternion *q)
{
double theta, snT;
//Anihilate numerical errors
     if (q->w>+1.) q->w=+1.;
else if (q->w<-1.) q->w=-1.;
//Calc theta
if ((theta=2.*acos(q->w))<TINY) { v->i=v->j=v->k=0.; }
else
  {
  snT=sin(.5*theta)/theta;
  v->i=q->i/snT;
  v->j=q->j/snT;
  v->k=q->k/snT;
  }
}

//This function calculates dq/dv
//Note. It normalizes input quaternion
void dimap_exp_quaternion(t_qtensor *dq,t_vec *v)
{
double theta, snT, csT;

theta=calc_vec_norm(v);
if (theta>TINY*TINY)
  {//If theta away from zero
  theta=sqrt(theta);
  snT=sin(.5*theta)/theta;
  csT=cos(.5*theta);
  (*dq)[0][0]=v->i*v->i*(.5*csT-snT)/theta/theta+snT;
  (*dq)[0][1]=(*dq)[1][0]=v->i*v->j*(.5*csT-snT)/theta/theta;
  (*dq)[0][2]=(*dq)[2][0]=v->i*v->k*(.5*csT-snT)/theta/theta;
  (*dq)[1][1]=v->j*v->j*(.5*csT-snT)/theta/theta+snT;
  (*dq)[1][2]=(*dq)[2][1]=v->j*v->k*(.5*csT-snT)/theta/theta;
  (*dq)[2][2]=v->k*v->k*(.5*csT-snT)/theta/theta+snT;
  (*dq)[0][3]=-.5*v->i*snT;
  (*dq)[1][3]=-.5*v->j*snT;
  (*dq)[2][3]=-.5*v->k*snT;
  }
else
  {//If theta close to zero
  snT=sinT(theta);
  (*dq)[0][0]=v->i*v->i*(theta/40.-1.)/24.+snT;
  (*dq)[0][1]=(*dq)[1][0]=v->i*v->j*(theta/40.-1.)/24.;
  (*dq)[0][2]=(*dq)[2][0]=v->i*v->k*(theta/40.-1.)/24.;
  (*dq)[1][1]=v->j*v->j*(theta/40.-1.)/24.+snT;
  (*dq)[1][2]=(*dq)[2][1]=v->j*v->k*(theta/40.-1.)/24.;
  (*dq)[2][2]=v->k*v->k*(theta/40.-1.)/24.+snT;
  (*dq)[0][3]=-.5*v->i*snT;
  (*dq)[1][3]=-.5*v->j*snT;
  (*dq)[2][3]=-.5*v->k*snT;
  }
}

//This function calculates rotation matrix from given quaternion
inline void calc_R_from_unit_quaternion(t_tensor *R,t_quaternion *q)
{
(*R)[0][0]=1.-2.*(q->j*q->j+q->k*q->k);
(*R)[0][1]=2.*(q->i*q->j-q->k*q->w);
(*R)[0][2]=2.*(q->i*q->k+q->j*q->w);
(*R)[1][0]=2.*(q->i*q->j+q->k*q->w);
(*R)[1][1]=1.-2.*(q->i*q->i+q->k*q->k);
(*R)[1][2]=2.*(q->j*q->k-q->i*q->w);
(*R)[2][0]=2.*(q->i*q->k-q->j*q->w);
(*R)[2][1]=2.*(q->j*q->k+q->i*q->w);
(*R)[2][2]=1.-2.*(q->i*q->i+q->j*q->j);
}

//This function calculates derivative of rotational tensor over quaternion components dr/dq=R(q)
inline void calc_dR_dquaternion(t_tensor *dRi,t_tensor *dRj,t_tensor *dRk,register t_qtensor *dq,t_quaternion *q)
{
t_quaternion _q;
_q.i=2.*q->i, _q.j=2.*q->j, _q.k=2.*q->k, _q.w=2.*q->w;
//dRi
(*dRi)[0][0]= _q.i*(*dq)[0][0]-_q.j*(*dq)[0][1]-_q.k*(*dq)[0][2]+_q.w*(*dq)[0][3];
(*dRi)[0][1]= _q.j*(*dq)[0][0]+_q.i*(*dq)[0][1]-_q.w*(*dq)[0][2]-_q.k*(*dq)[0][3];
(*dRi)[0][2]= _q.k*(*dq)[0][0]+_q.w*(*dq)[0][1]+_q.i*(*dq)[0][2]+_q.j*(*dq)[0][3];
(*dRi)[1][0]= _q.j*(*dq)[0][0]+_q.i*(*dq)[0][1]+_q.w*(*dq)[0][2]+_q.k*(*dq)[0][3];
(*dRi)[1][1]=-_q.i*(*dq)[0][0]+_q.j*(*dq)[0][1]-_q.k*(*dq)[0][2]+_q.w*(*dq)[0][3];
(*dRi)[1][2]=-_q.w*(*dq)[0][0]+_q.k*(*dq)[0][1]+_q.j*(*dq)[0][2]-_q.i*(*dq)[0][3];
(*dRi)[2][0]= _q.k*(*dq)[0][0]-_q.w*(*dq)[0][1]+_q.i*(*dq)[0][2]-_q.j*(*dq)[0][3];
(*dRi)[2][1]= _q.w*(*dq)[0][0]+_q.k*(*dq)[0][1]+_q.j*(*dq)[0][2]+_q.i*(*dq)[0][3];
(*dRi)[2][2]=-_q.i*(*dq)[0][0]-_q.j*(*dq)[0][1]+_q.k*(*dq)[0][2]+_q.w*(*dq)[0][3];
//dRj
(*dRj)[0][0]= _q.i*(*dq)[1][0]-_q.j*(*dq)[1][1]-_q.k*(*dq)[1][2]+_q.w*(*dq)[1][3];
(*dRj)[0][1]= _q.j*(*dq)[1][0]+_q.i*(*dq)[1][1]-_q.w*(*dq)[1][2]-_q.k*(*dq)[1][3];
(*dRj)[0][2]= _q.k*(*dq)[1][0]+_q.w*(*dq)[1][1]+_q.i*(*dq)[1][2]+_q.j*(*dq)[1][3];
(*dRj)[1][0]= _q.j*(*dq)[1][0]+_q.i*(*dq)[1][1]+_q.w*(*dq)[1][2]+_q.k*(*dq)[1][3];
(*dRj)[1][1]=-_q.i*(*dq)[1][0]+_q.j*(*dq)[1][1]-_q.k*(*dq)[1][2]+_q.w*(*dq)[1][3];
(*dRj)[1][2]=-_q.w*(*dq)[1][0]+_q.k*(*dq)[1][1]+_q.j*(*dq)[1][2]-_q.i*(*dq)[1][3];
(*dRj)[2][0]= _q.k*(*dq)[1][0]-_q.w*(*dq)[1][1]+_q.i*(*dq)[1][2]-_q.j*(*dq)[1][3];
(*dRj)[2][1]= _q.w*(*dq)[1][0]+_q.k*(*dq)[1][1]+_q.j*(*dq)[1][2]+_q.i*(*dq)[1][3];
(*dRj)[2][2]=-_q.i*(*dq)[1][0]-_q.j*(*dq)[1][1]+_q.k*(*dq)[1][2]+_q.w*(*dq)[1][3];
//dRk
(*dRk)[0][0]= _q.i*(*dq)[2][0]-_q.j*(*dq)[2][1]-_q.k*(*dq)[2][2]+_q.w*(*dq)[2][3];
(*dRk)[0][1]= _q.j*(*dq)[2][0]+_q.i*(*dq)[2][1]-_q.w*(*dq)[2][2]-_q.k*(*dq)[2][3];
(*dRk)[0][2]= _q.k*(*dq)[2][0]+_q.w*(*dq)[2][1]+_q.i*(*dq)[2][2]+_q.j*(*dq)[2][3];
(*dRk)[1][0]= _q.j*(*dq)[2][0]+_q.i*(*dq)[2][1]+_q.w*(*dq)[2][2]+_q.k*(*dq)[2][3];
(*dRk)[1][1]=-_q.i*(*dq)[2][0]+_q.j*(*dq)[2][1]-_q.k*(*dq)[2][2]+_q.w*(*dq)[2][3];
(*dRk)[1][2]=-_q.w*(*dq)[2][0]+_q.k*(*dq)[2][1]+_q.j*(*dq)[2][2]-_q.i*(*dq)[2][3];
(*dRk)[2][0]= _q.k*(*dq)[2][0]-_q.w*(*dq)[2][1]+_q.i*(*dq)[2][2]-_q.j*(*dq)[2][3];
(*dRk)[2][1]= _q.w*(*dq)[2][0]+_q.k*(*dq)[2][1]+_q.j*(*dq)[2][2]+_q.i*(*dq)[2][3];
(*dRk)[2][2]=-_q.i*(*dq)[2][0]-_q.j*(*dq)[2][1]+_q.k*(*dq)[2][2]+_q.w*(*dq)[2][3];
}

//----------------------------------   L I N E A R     D E R I V A T I V E S    P A R T  ---------------------------------------

//This functin calculates first derivative of two interaction points in three dimensions (normal bonds)
inline void calc_bond_derivative(double *r,t_vec *ri,t_vec *rj,double *dr)
{
dr[0]=ri->i-rj->i;                                       //  1 flops
dr[1]=ri->j-rj->j;                                       //  1 flops
dr[2]=ri->k-rj->k;                                       //  1 flops
*r=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);            //  5 flops + 1 sqrt
if (*r<SMALL2)
  {
  dr[3]=-(dr[0]/=SMALL2);                                      //  1 flops  
  dr[4]=-(dr[1]/=SMALL2);                                      //  1 flops
  dr[5]=-(dr[2]/=SMALL2);                                      //  1 flops
  }
else
  {
  dr[3]=-(dr[0]/=*r);                                      //  1 flops  
  dr[4]=-(dr[1]/=*r);                                      //  1 flops
  dr[5]=-(dr[2]/=*r);                                      //  1 flops
  }
//TOTAL:                                                    11 flops + 1 sqrt
}

//This functin calculates first derivatives of three interaction points in three dimensions (normal angles)
inline void calc_angle_derivative(double *csA,double *snA,t_vec *ri,t_vec *rj,t_vec *rk,double *dr)
{
t_vec a, b;
double aa, bb, ab;

a.i=ri->i-rj->i, a.j=ri->j-rj->j, a.k=ri->k-rj->k; //  3 flops
b.i=rk->i-rj->i, b.j=rk->j-rj->j, b.k=rk->k-rj->k; //  3 flops
if ((aa=calc_vec_norm(&a))<SMALL2) aa=SMALL2;      //  5 flops
if ((bb=calc_vec_norm(&b))<SMALL2) bb=SMALL2;      //  5 flops
ab=sqrt(aa*bb); //Just one square root required        1 flops + 1 sqrt
*csA=calc_vec_vec_scalar_product(&a,&b)/ab;        //  6 flops
//Check the accuracy of calculations
     if (*csA>=+1.-SMALL2) 
       {
       *csA=+1., *snA=0.;
       *csA=-1., *snA=0.;
       ab*=SMALL2;                                         //  1 flops
       aa*=+SMALL2;                                        //  2 flops
       bb*=+SMALL2;                                        //  2 flops
//Calculate the first derivative of angle over the radius-vectors
       dr[0]=a.i/aa-b.i/ab, dr[1]=a.j/aa-b.j/ab, dr[2]=a.k/aa-b.k/ab; //  9 flops
       dr[6]=b.i/bb-a.i/ab, dr[7]=b.j/bb-a.j/ab, dr[8]=b.k/bb-a.k/ab; //  9 flops
       }
else if (*csA<=-1.+SMALL2) 
       {
       *csA=-1., *snA=0.;
       ab*=SMALL2;                                         //  1 flops
       aa*=-SMALL2;                                        //  2 flops
       bb*=-SMALL2;                                        //  2 flops
//Calculate the first derivative of angle over the radius-vectors
       dr[0]=a.i/aa-b.i/ab, dr[1]=a.j/aa-b.j/ab, dr[2]=a.k/aa-b.k/ab; //  9 flops
       dr[6]=b.i/bb-a.i/ab, dr[7]=b.j/bb-a.j/ab, dr[8]=b.k/bb-a.k/ab; //  9 flops
       }
else if ( (*csA>0.-SMALL2)&&(*csA<0+SMALL2) )
       {
       *snA=sqrt(1.-*csA**csA);
       ab*=*snA;                                    //  1 flops
       dr[0]=-b.i/ab, dr[1]=-b.j/ab, dr[2]=-b.k/ab; //  9 flops
       dr[6]=-a.i/ab, dr[7]=-a.j/ab, dr[8]=-a.k/ab; //  9 flops
       }
else   {
       *snA=sqrt(1.-*csA**csA);                   //  2 flops + 1 sqrt
       ab*=*snA;                                           //  1 flops
       aa*=*snA/(*csA);                                    //  2 flops
       bb*=*snA/(*csA);                                    //  2 flops
//Calculate the first derivative of angle over the radius-vectors
       dr[0]=a.i/aa-b.i/ab, dr[1]=a.j/aa-b.j/ab, dr[2]=a.k/aa-b.k/ab; //  9 flops
       dr[6]=b.i/bb-a.i/ab, dr[7]=b.j/bb-a.j/ab, dr[8]=b.k/bb-a.k/ab; //  9 flops
       }
dr[3]=-dr[0]-dr[6],  dr[4]=-dr[1]-dr[7],  dr[5]=-dr[2]-dr[8];  //  3 flops
//TOTAL:                                              51 flops + 2 sqrt
}

//This functin calculates first derivatives of four interaction points in three dimensions (dihedral angles)
inline void calc_dih_derivative(double *csA,double *snA,t_vec  *ri,t_vec *rj,t_vec *rk,t_vec *rl,double *dr)
{
t_vec f,g,h,a,b;
double sqrtgg,aa,bb;

//Calc init vectors
f.i=ri->i-rj->i, f.j=ri->j-rj->j, f.k=ri->k-rj->k; // F                          3 flops
g.i=rj->i-rk->i, g.j=rj->j-rk->j, g.k=rj->k-rk->k; // G                          3 flops
h.i=rl->i-rk->i, h.j=rl->j-rk->j, h.k=rl->k-rk->k; // H                          3 flops
vec_vec_vmult(&a,&f,&g);                           // A                          9 flops
vec_vec_vmult(&b,&h,&g);                           // B                          9 flops
//Calculate cos(alpha) and sin(alpha)
if ((aa=calc_vec_norm(&a))<SMALL2) aa=1./SMALL2; else aa=1./aa;  //              6 flops      
if ((bb=calc_vec_norm(&b))<SMALL2) bb=1./SMALL2; else bb=1./bb;  //              6 flops
sqrtgg=sqrt(aa*bb);                                //                            1 flops + 1 root
*csA=calc_vec_vec_scalar_product(&a,&b)*sqrtgg;    //                            6 flops
*snA=calc_vec_vec_scalar_product(&a,&h)*sqrtgg;    //                            6 flops
sqrtgg=sqrt(calc_vec_norm(&g));                    //                            5 flops + 1 root
*snA*=sqrtgg;                                      //                            1 flops
     if (*csA>+1.) *csA=+1.;
else if (*csA<-1.) *csA=-1.;
else if (*snA>+1.) *snA=+1.;
else if (*snA<-1.) *snA=-1.;
//Calculate derivatives
multiple_vec_scalar(&a,&a,aa); //                                                3 flops
multiple_vec_scalar(&b,&b,bb); //                                                3 flops
multiple_vec_scalar((t_vec*)&dr[0x0],&a,-sqrtgg);                             // 3 flops  
multiple_vec_scalar((t_vec*)&dr[0x9],&b,+sqrtgg);                             // 3 flops  
multiple_vec_scalar(&a,&a,+calc_vec_vec_scalar_product(&f,&g)/sqrtgg);        // 9 flops  
multiple_vec_scalar(&b,&b,-calc_vec_vec_scalar_product(&h,&g)/sqrtgg);        // 9 flops  
dr[0x3]=a.i+b.i, dr[0x4]=a.j+b.j, dr[0x5]=a.k+b.k;                            // 3 flops
dr[0x6]=-dr[0x9]-dr[0x3], dr[0x7]=-dr[0xA]-dr[0x4], dr[0x8]=-dr[0xB]-dr[0x5]; // 3 flops
dr[0x3]-=dr[0x0],         dr[0x4]-=dr[0x1],         dr[0x5]-=dr[0x2];         // 3 flops
//TOTAL:                                                                      97 flops + 2 sqrt
}

//This functin calculates two first derivatives of two interaction points in three dimensions (normal bonds).Hardly imagine it can be usefull for something, but...
inline void calc_bond_derivatives(double *r,t_vec *ri,t_vec *rj,double *dr,double *ddr)
{
register double _d;
//Calculate bond value
dr[0]=ri->i-rj->i, dr[1]=ri->j-rj->j, dr[2]=ri->k-rj->k; //     3 flops
(*r)=sqrt(calc_vec_norm((t_vec*)&dr[0])), _d=1./(*r);    //     6 flops + 1 root
//Calculate first derivative
dr[3]=dr[0]*=_d, dr[4]=dr[1]*=_d, dr[5]=dr[2]*=_d;   //         3 flps     
dr[3]*=-_d,      dr[4]*=-_d,      dr[5]*=-_d;        //         3 flops                
//Calculate the second derivative
ddr[ 0]=dr[0]*dr[3],
ddr[27]=ddr[ 6]=dr[1]*dr[3], ddr[ 7]=dr[1]*dr[4],
ddr[33]=ddr[12]=dr[2]*dr[3], ddr[34]=ddr[13]=dr[2]*dr[4], ddr[14]=dr[2]*dr[5]; // 6 flops 
ddr[21]=ddr[00]+=_d, ddr[28]=ddr[07]+=_d, ddr[35]=ddr[14]+=_d;  // 3 flops
ddr[18]=-ddr[ 0], 
ddr[24]=-ddr[ 6], ddr[25]=-ddr[ 7], 
ddr[30]=-ddr[12], ddr[31]=-ddr[13], ddr[32]=-ddr[14]; 
ddr[19]=+ddr[24], ddr[20]=+ddr[30], ddr[26]=+ddr[31];
dr[3]=-dr[0], dr[4]=-dr[1], dr[5]=-dr[2];   // Restore first derivatives     
//TOTAL                                                        30 flops + 1 root
}

///The above code is more efficient for calculating of first derivative but this one utilizes items for second derivative calculations
//This functin calculates two first derivatives of three interaction points in three dimensions (normal angles)
inline void calc_angle_derivatives(double *csA,double *snA,t_vec *ri,t_vec *rj,t_vec *rk,double *dr,double *ddr)
{
t_vec *a, *b, *_t;
double aa, bb, ab, _d, _aa, _bb, _ab;
t_tensor *AXB;
register unsigned int _i,_j;

a=(t_vec*)&ddr[0x4*0x9+0x0], b=(t_vec*)&ddr[0x4*0x9+0x3], _t=(t_vec*)&ddr[0x4*0x9+0x6];
a->i=ri->i-rj->i, a->j=ri->j-rj->j, a->k=ri->k-rj->k; //       3 flops 
b->i=rk->i-rj->i, b->j=rk->j-rj->j, b->k=rk->k-rj->k; //       3 flops
aa=calc_vec_norm(a);  //                                       5 flops
bb=calc_vec_norm(b);  //                                       5 flops
_ab=sqrt(aa*bb); //Just one square root required               1 flops + 1 sqrt
ab=calc_vec_vec_scalar_product(a,b); //                        5 flops
*csA=ab/_ab;   //                                              1 flops    
//Check the accuracy of calculations
     if (*csA>+1.) { *csA=+1., *snA=0.; }
else if (*csA<-1.) { *csA=-1., *snA=0.; }
else     *snA=sqrt(1.-*csA**csA);  //                          2 flops + 1 sqrt
_ab=*snA*_ab;      //                                           1 flops
_aa=aa**snA/(*csA); //                                          2 flops
_bb=bb**snA/(*csA); //                                          2 flops
//Calculate the first derivative of angle over the radius-vectors                   
dr[0]=a->i/_aa-b->i/_ab, dr[1]=a->j/_aa-b->j/_ab, dr[2]=a->k/_aa-b->k/_ab; //  9 flops
dr[6]=b->i/_bb-a->i/_ab, dr[7]=b->j/_bb-a->j/_ab, dr[8]=b->k/_bb-a->k/_ab; //  9 flops
dr[3]=-dr[0]-dr[6],  dr[4]=-dr[1]-dr[7],  dr[5]=-dr[2]-dr[8];        //  3 flops
//Calculate second derivatives
AXB=(t_tensor*)&ddr[0x3*0x9];
(*AXB)[0][0]=a->i*b->i, (*AXB)[0][1]=a->i*b->j, (*AXB)[0][2]=a->i*b->k, 
(*AXB)[1][0]=a->j*b->i, (*AXB)[1][1]=a->j*b->j, (*AXB)[1][2]=a->j*b->k,
(*AXB)[2][0]=a->k*b->i, (*AXB)[2][1]=a->k*b->j, (*AXB)[2][2]=a->k*b->k;    //  9 flops
//Calculate d^2a/dri^2=d^2a/dA^2
_bb=bb/_ab, _aa=2./aa, _d=1./_ab/aa; //                                          4 flops
_t->i=_bb*dr[6]-_aa*a->i, _t->j=_bb*dr[7]-_aa*a->j, _t->k=_bb*dr[8]-_aa*a->k; //  9 flops
ddr[ 0]=_t->i*dr[0]-_d*((*AXB)[0][0]-ab),
ddr[ 9]=_t->j*dr[0]-_d*(2.*(*AXB)[1][0]-(*AXB)[0][1]), ddr[10]=_t->j*dr[1]-_d*((*AXB)[1][1]-ab),
ddr[18]=_t->k*dr[0]-_d*(2.*(*AXB)[2][0]-(*AXB)[0][2]), ddr[19]=_t->k*dr[1]-_d*(2.*(*AXB)[2][1]-(*AXB)[1][2]), ddr[20]=_t->k*dr[2]-_d*((*AXB)[2][2]-ab); 
//ddr[ 1]=ddr[ 9], ddr[ 2]=ddr[18], ddr[11]=ddr[19];  //                           27 flops   
//Calculate d^2a/drk^2=d^2a/dB^2
_aa=aa/_ab, _bb=2./bb, _d=1./_ab/bb; //                                          4 flops
_t->i=_aa*dr[0]-_bb*b->i, _t->j=_aa*dr[1]-_bb*b->j, _t->k=_aa*dr[2]-_bb*b->k; //  9 flops
ddr[60]=_t->i*dr[6]-_d*((*AXB)[0][0]-ab),
ddr[69]=_t->j*dr[6]-_d*(2.*(*AXB)[0][1]-(*AXB)[1][0]), ddr[70]=_t->j*dr[7]-_d*((*AXB)[1][1]-ab),
ddr[78]=_t->k*dr[6]-_d*(2.*(*AXB)[0][2]-(*AXB)[2][0]), ddr[79]=_t->k*dr[7]-_d*(2.*(*AXB)[1][2]-(*AXB)[2][1]), ddr[80]=_t->k*dr[8]-_d*((*AXB)[2][2]-ab); 
//ddr[61]=ddr[69], ddr[62]=ddr[78], ddr[71]=ddr[79];  //                           27 flops   
//Calculate d^2a/dri/drk=d^2a/dA/dB
_bb=bb/_ab;          //                                                             1 flops
_t->i=_bb*dr[6], _t->j=_bb*dr[7], _t->k=_bb*dr[8]; //                               3 flops
ddr[54]=_t->i*dr[6]+_d*(b->i*b->i-bb),
ddr[63]=_t->j*dr[6]+_d*(b->j*b->i),    ddr[64]=_t->j*dr[7]+_d*(b->j*b->j-bb),
ddr[72]=_t->k*dr[6]+_d*(b->k*b->i),    ddr[73]=_t->k*dr[7]+_d*(b->k*b->j),    ddr[74]=_t->k*dr[8]+_d*(b->k*b->k-bb); 
ddr[55]=ddr[63], ddr[56]=ddr[72], ddr[65]=ddr[73];  //                           27 flops
//Calculate d^2a/dri/drj=-d^2a/dA^2-d^2a/dA/dB
ddr[27]=-ddr[ 0]-ddr[54],
ddr[36]=-ddr[ 9]-ddr[63], ddr[37]=-ddr[10]-ddr[64],
ddr[45]=-ddr[18]-ddr[72], ddr[46]=-ddr[19]-ddr[73], ddr[47]=-ddr[20]-ddr[74];  
ddr[28]=ddr[36], ddr[29]=ddr[45], ddr[38]=ddr[46];  //                            6 flops
//Calculate d^2/drj/drk=-d^2a/dB^2-d^2a/dA/dB
ddr[57]=-ddr[60]-ddr[54],
ddr[66]=-ddr[69]-ddr[63], ddr[67]=-ddr[70]-ddr[64],
ddr[75]=-ddr[78]-ddr[72], ddr[76]=-ddr[79]-ddr[73], ddr[77]=-ddr[80]-ddr[74];
ddr[58]=ddr[66], ddr[59]=ddr[75], ddr[68]=ddr[76];  //                            6 flops
//Calculate d^2a/drj^2=d^2a/dA^2+2*d^2a/dA/dB+d^2a/dB^2=-d^2a/dridrj-d^2a/drjdrk
ddr[30]=-ddr[27]-ddr[57],
ddr[39]=-ddr[36]-ddr[66], ddr[40]=-ddr[37]-ddr[67],
ddr[48]=-ddr[45]-ddr[75], ddr[49]=-ddr[46]-ddr[76], ddr[50]=-ddr[47]-ddr[77];
//ddr[31]=ddr[39], ddr[32]=ddr[48], ddr[41]=ddr[49];  //                          6 flops
_i=9;
while(_i--)
  {
  _j=_i;
  while (_j--)
    {
    ddr[_j*9+_i]=ddr[_i*9+_j]; 
    }
  }
//TOTAL                                                                         174 flops + 2 root  
}

//This functin calculates two first derivatives of four interaction points in three dimensions (dihedral angles)
inline void calc_dih_angle_derivatives(double *csA,double *snA,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl,double *dr,double *ddr)
{
t_vec f,g,h,a,b,v;
double sqrtgg,aa,bb,fgg,hgg,_d;
t_tensor *aga,*afa,*bgb,*bhb;
register unsigned int _i,_j;

//Calc init vectors
f.i=ri->i-rj->i, f.j=ri->j-rj->j, f.k=ri->k-rj->k; // F                        3 flops
g.i=rj->i-rk->i, g.j=rj->j-rk->j, g.k=rj->k-rk->k; // G                        3 flops
h.i=rl->i-rk->i, h.j=rl->j-rk->j, h.k=rl->k-rk->k; // H                        3 flops
vec_vec_vmult(&a,&f,&g);                           // A                        9 flops
vec_vec_vmult(&b,&h,&g);                           // B                        9 flops
//Calculate cos(alpha) and sin(alpha)
aa=1./calc_vec_norm(&a);                           //                          6 flops      
bb=1./calc_vec_norm(&b);                           //                          6 flops
sqrtgg=sqrt(aa*bb);                                //                          1 flops + 1 root
*csA=calc_vec_vec_scalar_product(&a,&b)*sqrtgg;    //                          6 flops
*snA=calc_vec_vec_scalar_product(&a,&h)*sqrtgg;    //                          6 flops
sqrtgg=sqrt(calc_vec_norm(&g));                    //                          5 flops + 1 root
*snA*=sqrtgg;                                      //                          1 flops
     if (*csA>+1.) *csA=+1.;
else if (*csA<-1.) *csA=-1.;
else if (*snA>+1.) *snA=+1.;
else if (*snA<-1.) *snA=-1.;
//Calculate second derivatives
fgg=calc_vec_vec_scalar_product(&f,&g)/sqrtgg;     //                          6 flops
hgg=calc_vec_vec_scalar_product(&h,&g)/sqrtgg;     //                          6 flops
aga=(t_tensor*)&ddr[ 03], afa=(t_tensor*)&ddr[ 15], bgb=(t_tensor*)&ddr[ 27], bhb=(t_tensor*)&ddr[ 42]; //Setup upper row of coordinates 82
//fill tensors                                                                72 flops
vec_vec_vmult(&v,&g,&a);
(*aga)[0][0]=v.i*a.i, (*aga)[0][1]=v.i*a.j, (*aga)[0][2]=v.i*a.k, (*aga)[1][0]=v.j*a.i, (*aga)[1][1]=v.j*a.j, (*aga)[1][2]=v.j*a.k, (*aga)[2][0]=v.k*a.i, (*aga)[2][1]=v.k*a.j, (*aga)[2][2]=v.k*a.k;
vec_vec_vmult(&v,&f,&a);
(*afa)[0][0]=v.i*a.i, (*afa)[0][1]=v.i*a.j, (*afa)[0][2]=v.i*a.k, (*afa)[1][0]=v.j*a.i, (*afa)[1][1]=v.j*a.j, (*afa)[1][2]=v.j*a.k, (*afa)[2][0]=v.k*a.i, (*afa)[2][1]=v.k*a.j, (*afa)[2][2]=v.k*a.k;
vec_vec_vmult(&v,&g,&b);
(*bgb)[0][0]=v.i*b.i, (*bgb)[0][1]=v.i*b.j, (*bgb)[0][2]=v.i*b.k, (*bgb)[1][0]=v.j*b.i, (*bgb)[1][1]=v.j*b.j, (*bgb)[1][2]=v.j*b.k, (*bgb)[2][0]=v.k*b.i, (*bgb)[2][1]=v.k*b.j, (*bgb)[2][2]=v.k*b.k;
vec_vec_vmult(&v,&h,&b);
(*bhb)[0][0]=v.i*b.i, (*bhb)[0][1]=v.i*b.j, (*bhb)[0][2]=v.i*b.k, (*bhb)[1][0]=v.j*b.i, (*bhb)[1][1]=v.j*b.j, (*bhb)[1][2]=v.j*b.k, (*bhb)[2][0]=v.k*b.i, (*bhb)[2][1]=v.k*b.j, (*bhb)[2][2]=v.k*b.k;
//Here we use v as a storage of temporary variables
//Fill derivatives ddr[l][l]=d^2A/dH^2                                        20 flops 
_d=-sqrtgg*bb*bb;
ddr[117]=_d*((*bgb)[0][0]+(*bgb)[0][0]),
ddr[129]=_d*((*bgb)[0][1]+(*bgb)[1][0]), ddr[130]=_d*((*bgb)[1][1]+(*bgb)[1][1]),
ddr[141]=_d*((*bgb)[0][2]+(*bgb)[2][0]), ddr[142]=_d*((*bgb)[1][2]+(*bgb)[2][1]), ddr[143]=_d*((*bgb)[2][2]+(*bgb)[2][2]);
//Fill derivatives part ddr[i][i]=+d^2A/dF^2 ... and overwrite tensors!!!!    20 flops
_d=+sqrtgg*aa*aa;
ddr[  0]=_d*((*aga)[0][0]+(*aga)[0][0]), 
ddr[ 12]=_d*((*aga)[1][0]+(*aga)[0][1]), ddr[ 13]=_d*((*aga)[1][1]+(*aga)[1][1]), 
ddr[ 24]=_d*((*aga)[2][0]+(*aga)[0][2]), ddr[ 25]=_d*((*aga)[2][1]+(*aga)[1][2]), ddr[ 26]=_d*((*aga)[2][2]+(*aga)[2][2]);
//Fill derivatives part ddr[l][j]=+d^2A/dH/dG                                 36 flops
_d=+bb*bb;
ddr[111]=_d*(sqrtgg*(*bhb)[0][0]+hgg*(*bgb)[0][0]), ddr[112]=_d*(sqrtgg*(*bhb)[0][1]+hgg*(*bgb)[1][0]), ddr[113]=_d*(sqrtgg*(*bhb)[0][2]+hgg*(*bgb)[2][0]),
ddr[123]=_d*(sqrtgg*(*bhb)[1][0]+hgg*(*bgb)[0][1]), ddr[124]=_d*(sqrtgg*(*bhb)[1][1]+hgg*(*bgb)[1][1]), ddr[125]=_d*(sqrtgg*(*bhb)[1][2]+hgg*(*bgb)[2][1]),
ddr[135]=_d*(sqrtgg*(*bhb)[2][0]+hgg*(*bgb)[0][2]), ddr[136]=_d*(sqrtgg*(*bhb)[2][1]+hgg*(*bgb)[1][2]), ddr[137]=_d*(sqrtgg*(*bhb)[2][2]+hgg*(*bgb)[2][2]);
//Fill derivatives part ddr[j][k]=-d^2A/dF/dG                                 36 flops
_d=-aa*aa;
ddr[ 72]=_d*(sqrtgg*(*afa)[0][0]+fgg*(*aga)[0][0]), ddr[ 73]=_d*(sqrtgg*(*afa)[1][0]+fgg*(*aga)[0][1]), ddr[ 74]=_d*(sqrtgg*(*afa)[2][0]+fgg*(*aga)[0][2]),
ddr[ 84]=_d*(sqrtgg*(*afa)[0][1]+fgg*(*aga)[1][0]), ddr[ 85]=_d*(sqrtgg*(*afa)[1][1]+fgg*(*aga)[1][1]), ddr[ 86]=_d*(sqrtgg*(*afa)[2][1]+fgg*(*aga)[1][2]),
ddr[ 96]=_d*(sqrtgg*(*afa)[0][2]+fgg*(*aga)[2][0]), ddr[ 97]=_d*(sqrtgg*(*afa)[1][2]+fgg*(*aga)[2][1]), ddr[ 98]=_d*(sqrtgg*(*afa)[2][2]+fgg*(*aga)[2][2]);
//Fill derivatives part ddr[k][j]=-d^2A/dG^2                                 111 flops
//Note. Fill it particular only to save on future calculations
v.i=1./(2.*sqrd(sqrd(sqrtgg)))/aa, v.j=1./(2.*sqrd(sqrd(sqrtgg)))/bb, fgg*=aa*aa, hgg*=bb*bb; //hgg and fgg are updated now!
ddr[ 39]=ddr[ 78]=v.i*ddr[ 00]+v.j*ddr[117]+fgg*((*afa)[0][0]+(*afa)[0][0])-hgg*((*bhb)[0][0]+(*bhb)[0][0]),
ddr[ 51]=ddr[ 90]=v.i*ddr[ 12]+v.j*ddr[129]+fgg*((*afa)[1][0]+(*afa)[0][1])-hgg*((*bhb)[1][0]+(*bhb)[0][1]), ddr[ 52]=ddr[ 91]=v.i*ddr[ 13]+v.j*ddr[130]+fgg*((*afa)[1][1]+(*afa)[1][1])-hgg*((*bhb)[1][1]+(*bhb)[1][1]),
ddr[ 63]=ddr[102]=v.i*ddr[ 24]+v.j*ddr[141]+fgg*((*afa)[2][0]+(*afa)[0][2])-hgg*((*bhb)[2][0]+(*bhb)[0][2]), ddr[ 64]=ddr[103]=v.i*ddr[ 25]+v.j*ddr[142]+fgg*((*afa)[2][1]+(*afa)[1][2])-hgg*((*bhb)[2][1]+(*bhb)[1][2]), ddr[ 65]=ddr[104]=v.i*ddr[ 26]+v.j*ddr[143]+fgg*((*afa)[2][2]+(*afa)[2][2])-hgg*((*bhb)[2][2]+(*bhb)[2][2]);
ddr[ 75]=-ddr[ 39], ddr[ 76]=-ddr[ 51], ddr[ 77]=-ddr[ 63],
ddr[ 87]=-ddr[ 51], ddr[ 88]=-ddr[ 52], ddr[ 89]=-ddr[ 64],
ddr[ 99]=-ddr[ 63], ddr[100]=-ddr[ 64], ddr[101]=-ddr[ 65];
//Fill derivatives part ddr[l][k]=-d^2A/dH^2+d^2A/dH/dG                        9 flops  
ddr[114]=-ddr[117]-ddr[111], ddr[115]=-ddr[129]-ddr[112], ddr[116]=-ddr[141]-ddr[113],
ddr[126]=-ddr[129]-ddr[123], ddr[127]=-ddr[130]-ddr[124], ddr[128]=-ddr[142]-ddr[125],
ddr[138]=-ddr[141]-ddr[135], ddr[139]=-ddr[142]-ddr[136], ddr[140]=-ddr[143]-ddr[137];
//Fill derivatives part ddr[j][i]=-d^2A/dF^2+d^2A/dF/dG                        9 flops
ddr[ 36]=-ddr[  0]+ddr[ 72], ddr[ 37]=-ddr[ 12]+ddr[ 73], ddr[ 38]=-ddr[ 24]+ddr[ 74],
ddr[ 48]=-ddr[ 12]+ddr[ 84], ddr[ 49]=-ddr[ 13]+ddr[ 85], ddr[ 50]=-ddr[ 25]+ddr[ 86],
ddr[ 60]=-ddr[ 24]+ddr[ 96], ddr[ 61]=-ddr[ 25]+ddr[ 97], ddr[ 62]=-ddr[ 26]+ddr[ 98];
//Fill derivatives part ddr[k][k]=d^2A/dG^2+d^2A/dH^2+d^2A/dG/dH              18 flops
ddr[ 78]+=ddr[117]+ddr[111]+ddr[111], 
ddr[ 90]+=ddr[129]+ddr[123]+ddr[112], ddr[ 91]+=ddr[130]+ddr[124]+ddr[124],
ddr[102]+=ddr[141]+ddr[135]+ddr[113], ddr[103]+=ddr[142]+ddr[136]+ddr[125], ddr[104]+=ddr[143]+ddr[137]+ddr[137];
//Fill derivatives part ddr[j][j]=d^2A/dG^2+d^2A/dF^2+d^2A/dG/dF              18 flops
ddr[ 39]+=ddr[  0]-ddr[ 72]-ddr[ 72], 
ddr[ 51]+=ddr[ 12]-ddr[ 84]-ddr[ 73], ddr[ 52]+=ddr[ 13]-ddr[ 85]-ddr[ 85],
ddr[ 63]+=ddr[ 24]-ddr[ 96]-ddr[ 74], ddr[ 64]+=ddr[ 25]-ddr[ 97]-ddr[ 86], ddr[ 65]+=ddr[ 26]-ddr[ 98]-ddr[ 98];
//Fill derivatives part ddr[k][j]=-d^2A/dG^2+d^2A/dG/dF+d^2A/dG/dH            18 flops
ddr[ 75]+=ddr[ 72]-ddr[111], ddr[ 76]+=ddr[ 73]-ddr[112], ddr[ 77]+=ddr[ 74]-ddr[113],
ddr[ 87]+=ddr[ 84]-ddr[123], ddr[ 88]+=ddr[ 85]-ddr[124], ddr[ 89]+=ddr[ 86]-ddr[125],
ddr[ 99]+=ddr[ 96]-ddr[135], ddr[100]+=ddr[ 97]-ddr[136], ddr[101]+=ddr[ 98]-ddr[137];
//Fill derivatives part ddr[k][i]=0
ddr[108]=ddr[109]=ddr[110]=0., ddr[120]=ddr[121]=ddr[122]=0., ddr[132]=ddr[133]=ddr[134]=0.;
//End of d^2A/dr/dr part

//Calculate first derivatives now 417
multiple_vec_scalar(&a,&a,aa); //                                                3 flops
multiple_vec_scalar(&b,&b,bb); //                                                3 flops
multiple_vec_scalar((t_vec*)&dr[0x0],&a,-sqrtgg);                             // 3 flops  
multiple_vec_scalar((t_vec*)&dr[0x9],&b,+sqrtgg);                             // 3 flops  
multiple_vec_scalar(&a,&a,+calc_vec_vec_scalar_product(&f,&g)/sqrtgg);        // 9 flops  
multiple_vec_scalar(&b,&b,-calc_vec_vec_scalar_product(&h,&g)/sqrtgg);        // 9 flops  
dr[0x3]=a.i+b.i, dr[0x4]=a.j+b.j, dr[0x5]=a.k+b.k;                            // 3 flops
dr[0x6]=-dr[0x9]-dr[0x3], dr[0x7]=-dr[0xA]-dr[0x4], dr[0x8]=-dr[0xB]-dr[0x5]; // 3 flops
dr[0x3]-=dr[0x0],         dr[0x4]-=dr[0x1],         dr[0x5]-=dr[0x2];         // 3 flops
//End of d^A/dr part
_i=12;
while(_i--)
  {
  _j=_i;
  while (_j--)
    {
    ddr[_j*12+_i]=ddr[_i*12+_j]; 
    }
  }
//TOTAL:                                                                     456 flops + 2 sqrt
}

//------------------------------------- 3 D   P A R T --------------------------------------------------------

//Taken from http://atheist4.narod.ru/mw/distance.htm the equations obtained as making two pairs vectors to be orthogonal - line and distance.
//This function calculates scale factors for both vectors of  minimal distance between two lines in 3D from given 2 points and 2 vectors
//NOTE. It returns  FALSE if lines are parallel (<SMALL2)
char _calc_line_line_scale_distance_3D(double *t1,double *t2,t_vec *X1,t_vec *n1,t_vec *X2,t_vec *n2)
{
register double c1,c2,c0,r1,r2;
c0=n1->i*n2->i+n1->j*n2->j+n1->k*n2->k;
c1=n1->i*n1->i+n1->j*n1->j+n1->k*n1->k;
c2=n2->i*n2->i+n2->j*n2->j+n2->k*n2->k;
if (fabs(1.-c0/c1/c2)<SMALL2) return FALSE; //Lines are parallel
r1=(X1->i-X2->i)*n1->i+(X1->j-X2->j)*n1->j+(X1->k-X2->k)*n1->k;
r2=(X2->i-X1->i)*n2->i+(X2->j-X1->j)*n2->j+(X2->k-X1->k)*n2->k;
*t1=r1*c2-r2*c0, *t2=c1*r2-c0*r1, c0*=c0, c0-=c1*c2;
*t1/=-c0, *t2/=-c0;
return TRUE; 
}
//This function calculates the distance between two lines in 3D from given 4 points
//NOTE. It returns distance and coordinates of orthogonal points that are NAN if lines are parallel.
inline double calc_line_line_distance_3D(t_vec *M,t_vec *N,t_vec *A,t_vec *B,t_vec *C,t_vec *D)
{
t_vec n1, n2;
double t1,t2;
n1.i=B->i-A->i, n1.j=B->j-A->j, n1.k=B->k-A->k;
n2.i=D->i-C->i, n2.j=D->j-C->j, n2.k=D->k-C->k;
if (!(_calc_line_line_scale_distance_3D(&t1,&t2,&n1,A,&n2,C)))
  {//Lines are colinear
  n2.i=C->i-A->i, n2.j=C->j-A->j, n2.k=C->k-A->k;
  if (M) M->i=M->j=M->k=(double)NAN;  
  if (N) N->i=N->j=N->k=(double)NAN;  
  ylib_errno=YERROR_LEGAL; 
  return fabsf(calc_vec_norm(&n2)-sqrd(calc_vec_vec_scalar_product(&n1,&n2))/calc_vec_norm(&n1));  
  }
else
  {
  if (M) { M->i=A->i+n1.i*t1, M->j=A->j+n1.j*t1, M->k=A->k+n1.k*t1; }
  if (N) { N->i=C->i+n2.i*t2, N->j=C->j+n2.j*t2, N->k=C->k+n2.k*t2; }
  return sqrd(A->i+n1.i*t1-C->i-n2.i*t2)+sqrd(A->j+n1.j*t1-C->j-n2.j*t2)+sqrd(A->k+n1.k*t1-C->k-n2.k*t2);
  }
}

//----------------------------- S O S   P A R T ------------------------------------------------------------

//------------------------------------- S o S    P A R T -------------------------------------------

//The symbolic math part. See Herbert Edelsbrunner and Ernst Peter Mucke Simulation of Somplicity: a technique to cope with degenerate cases in geometrical algorithms.

//1-1 comparisons
inline unsigned int bjsort2(void *p0,void *p1,void **_p0,void **_p1)
{
if (p0<p1) { *_p0=p0, *_p1=p1; return 0; }
else       { *_p0=p1, *_p1=p0; return 1; }
}
//2-3 comparisons
inline unsigned int _bjsort3(void *p0,void *p1,void *p2,void **_p0,void **_p1,void **_p2)
{
if (p0<p1)
  {
  if (p0<p2)
    {
    if (p1<p2)
      {// p0 < p1 < p2
      *_p0=p0, *_p1=p1, *_p2=p2; return 0;
      }
    else
      {// p0 < p2 < p1
      *_p0=p0, *_p1=p2, *_p2=p1; return 1;
      }  
    } 
  else
    {// p2 < p0 < p1
    *_p0=p2, *_p1=p0, *_p2=p1; return 2;   
    }
  }
else
  {
  if (p1<p2)
    {
    if (p0<p2)
      {// p1 < p0 < p2
      *_p0=p1, *_p1=p0, *_p2=p2; return 1;
      }
    else
      {// p1 < p2 < p0
      *_p0=p1, *_p1=p2, *_p2=p0; return 2;
      }  
    } 
  else
    {// p2 < p1 < p0
    *_p0=p2, *_p1=p1, *_p2=p0; return 1;
    }
  }
}
//3-6 comparisons
inline unsigned int bjsort3(void *p0,void *p1,void *p2,void **_p0,void **_p1,void **_p2)
{
unsigned int nexch=_bjsort3(p0,p1,p2,_p0,_p1,_p2);
if ( (*_p0>*_p1)||(*_p1>*_p2) ) 
  error_exit("Failure in bjsort3()");   
return nexch;
}
//This function resorts size of indexes stored in stack as previous but do not update stack of sorting. For manual call only!!!!
inline unsigned int _bjsort4(void *p0,void *p1,void *p2,void *p3,void **_p0,void **_p1,void **_p2,void **_p3)
{
if (p0<p1)
  {
  if (p0<p2)
    {
    if (p0<p3)
      { // p0 < { p1, p2, p3 }^3
      *_p0=p0; return bjsort3(p1,p2,p3,_p1,_p2,_p3);
      }
    else
      {
      if (p1<p2)
        {// p3 < p0 < p1 < p2
        *_p0=p3, *_p1=p0, *_p2=p1, *_p3=p2; return 3;
        }
      else
        {// p3 < p0 < p2 < p1
        *_p0=p3, *_p1=p0, *_p2=p2, *_p3=p1; return 2;
        }
      }  
    } 
  else
    {
    if (p2<p3)
      {
      if (p1<p3)
        {// p2 < p0 < p1 < p3
        *_p0=p2, *_p1=p0, *_p2=p1, *_p3=p3; return 2;
        }
      else
        {
        if (p0<p3)
          {// p2 < p0 < p3 < p1
          *_p0=p2, *_p1=p0, *_p2=p3, *_p3=p1; return 3; 
          }
        else
          {// p2 < p3 < p0 < p1
          *_p0=p2, *_p1=p3, *_p2=p0, *_p3=p1; return 2;
          }
        }
      }
    else
      {// p3 < p2 < p0 < p1
      *_p0=p3, *_p1=p2, *_p2=p0, *_p3=p1; return 3;
      }
    }
  }
else
  {
  if (p1<p2)
    {
    if (p1<p3)
      {// p1 < { p0, p2, p3 }^3
      *_p0=p1; return 1+bjsort3(p0,p2,p3,_p1,_p2,_p3);
      }
    else
      {
      if (p0<p2)
        {// p3 < p1 < p0 < p2
        *_p0=p3, *_p1=p1, *_p2=p0, *_p3=p2; return 2;
        }
      else
        {// p3 < p1 < p2 < p0
        *_p0=p3, *_p1=p1, *_p2=p2, *_p3=p0; return 1;
        }
      }  
    } 
  else
    {
    if (p2<p3)
      {
      if (p0<p3)
        {// p2 < p1 < p0 < p3 
        *_p0=p2, *_p1=p1, *_p2=p0, *_p3=p3; return 1;
        }
      else
        {
        if (p1<p3)
          {// p2 < p1 < p3 < p0
          *_p0=p2, *_p1=p1, *_p2=p3, *_p3=p0; return 2;
          }
        else
          {// p2 < p3 < p1 < p0
          *_p0=p2, *_p1=p3, *_p2=p1, *_p3=p0; return 3;
          }
        }
      }
    else
      {// p3 < p2 < p1 < p0
      *_p0=p3, *_p1=p2, *_p2=p1, *_p3=p0; return 2;
      }
    }
  }
}
inline unsigned int bjsort4(void *p0,void *p1,void *p2,void *p3,void **_p0,void **_p1,void **_p2,void **_p3)
{
unsigned int nexch=_bjsort4(p0,p1,p2,p3,_p0,_p1,_p2,_p3);
if ( (*_p0>*_p1)||(*_p1>*_p2)||(*_p2>*_p3) ) 
  error_exit("Failure in bjsort4()");   
return nexch;
}
// 7-16 comparisons
inline unsigned int bjsort5(void *p0,void *p1,void *p2,void *p3,void *p4,void **_p0,void **_p1,void **_p2,void **_p3,void **_p4)
{
unsigned int nexch;
     if ( (p0<p1)&&(p0<p2)&&(p0<p3)&&(p0<p4) ) { *_p0=p0; nexch=  bjsort4(p1,p2,p3,p4,_p1,_p2,_p3,_p4); }
else if (          (p1<p2)&&(p1<p3)&&(p1<p4) ) { *_p0=p1; nexch=1+bjsort4(p0,p2,p3,p4,_p1,_p2,_p3,_p4); }
else if (                   (p2<p3)&&(p2<p4) ) { *_p0=p2; nexch=1+bjsort4(p1,p0,p3,p4,_p1,_p2,_p3,_p4); }
else if (                            (p3<p4) ) { *_p0=p3; nexch=1+bjsort4(p1,p2,p0,p4,_p1,_p2,_p3,_p4); }
else                                           { *_p0=p4; nexch=1+bjsort4(p1,p2,p3,p0,_p1,_p2,_p3,_p4); }
if ( (*_p0>*_p1)||(*_p1>*_p2)||(*_p2>*_p3)||(*_p3>*_p4) ) 
  error_exit("Failure in bjsort5()");   
return nexch;
}

//isort function are the same as jsort but they replaces the address 
inline unsigned int isort2(void **_p0,void **_p1)
{
return bjsort2(*_p0,*_p1,_p0,_p1);
}
inline unsigned int isort3(void **_p0,void **_p1,void **_p2)
{
return bjsort3(*_p0,*_p1,*_p2,_p0,_p1,_p2);
}
inline unsigned int isort4(void **_p0,void **_p1,void **_p2,void **_p3)
{
return bjsort4(*_p0,*_p1,*_p2,*_p3,_p0,_p1,_p2,_p3);
}
inline unsigned int isort5(void **_p0,void **_p1,void **_p2,void **_p3,void **_p4)
{
return bjsort5(*_p0,*_p1,*_p2,*_p3,*_p4,_p0,_p1,_p2,_p3,_p4);
}

//This function performs 5d matix determinant sign calculations
char sign_ldet4(double **p)
{
double det_res;
det_res = calc_det4x4(p[0][0], p[0][1], p[0][2], 1.,
                      p[1][0], p[1][1], p[1][2], 1.,
                      p[2][0], p[2][1], p[2][2], 1.,
                      p[3][0], p[3][1], p[3][2], 1.);

//Checking if floating point calculations provide sufficient accuracy.
if (fabs(det_res) > PRECISION)
  return det_res > 0.;
//Floating point calculations unreliable, using SoS. lint_SoS returns 1 and -1,
//so these values should be converted to 1 and 0.
return lint_SoS_udet4(p[0][0], p[0][1], p[0][2],
                      p[1][0], p[1][1], p[1][2],
                      p[2][0], p[2][1], p[2][2],
                      p[3][0], p[3][1], p[3][2], 2) > 0;
}

//This function performs 4d matix determinant sign calculations
char sign_det4(double **p)
{
double det_res;
det_res = calc_det4x4(p[0][0], p[0][1], p[0][2], p[0][3],
                      p[1][0], p[1][1], p[1][2], p[1][3],
                      p[2][0], p[2][1], p[2][2], p[2][3],
                      p[3][0], p[3][1], p[3][2], p[3][3]);

//Checking if floating point calculations provide sufficient accuracy.
if (fabs(det_res) > PRECISION)
  return det_res > 0.;
//Floating point calculations unreliable, using SoS. lint_SoS returns 1 and -1,
//so these values should be converted to 1 and 0.
return lint_SoS_det4(p[0][0], p[0][1], p[0][2], p[0][3],
                     p[1][0], p[1][1], p[1][2], p[1][3],
                     p[2][0], p[2][1], p[2][2], p[2][3],
                     p[3][0], p[3][1], p[3][2], p[3][3], 2) > 0;
}

//This function performs 4d matix determinant sign calculations
char sign_ldet5(double **p)
{
double det_res;
det_res = calc_det5x5(p[0][0], p[0][1], p[0][2], p[0][3], 1.,
                      p[1][0], p[1][1], p[1][2], p[1][3], 1.,
                      p[2][0], p[2][1], p[2][2], p[2][3], 1.,
                      p[3][0], p[3][1], p[3][2], p[3][3], 1.,
                      p[4][0], p[4][1], p[4][2], p[4][3], 1.);

//Checking if floating point calculations provide sufficient accuracy.
if (fabs(det_res) > PRECISION)
  return det_res > 0.;
//Floating point calculations unreliable, using SoS. lint_SoS returns 1 and -1,
//so these values should be converted to 1 and 0.
return lint_SoS_udet5(p[0][0], p[0][1], p[0][2], p[0][3],
                      p[1][0], p[1][1], p[1][2], p[1][3],
                      p[2][0], p[2][1], p[2][2], p[2][3],
                      p[3][0], p[3][1], p[3][2], p[3][3],
                      p[4][0], p[4][1], p[4][2], p[4][3], 2) > 0;
}



//----------------------------- F I G U R E S    P A R T ----------------------------------------------------

//This function calculates the rmsd of two structures as difference in summ of corresponding atoms distances (sqares)
double calc_str_str_rmsd(unsigned int size,t_vec *A,t_vec *B)
{
register unsigned int _i;
double rmsd=0.00;

for (_i=0;_i<size;_i++)
  rmsd+=calc_distance(&A[_i],&B[_i]);
return sqrt(rmsd/(double)size);
}

//This function finds the coordinates of rectlinear tetrahedon that contain all points of given set in it. Function returns the edge lenght of tetrahedron.
double inscribe_tetrahedron(double scale,t_lvec *p0,t_lvec *p1,t_lvec *p2,t_lvec *p3,unsigned int size,t_lvec *lvec)
{
register double _r;
unsigned int _i;
//Calc circumsphere ceneter
p3->i=p3->j=p3->k=0.00;
for (_i=0;_i<size;_i++)
  {
  p3->i+=lvec[_i].i;
  p3->j+=lvec[_i].j;
  p3->k+=lvec[_i].k;
  }
//Calc circumsphere radius
multiple_vec_scalar((t_vec*)p3,(t_vec*)p3,1.00/(double)size);
_r=(sqrd(p3->i-lvec[0].i)+sqrd(p3->j-lvec[0].j)+sqrd(p3->k-lvec[0].k));
for (_i=1;_i<size;_i++)
  if ((sqrd(p3->i-lvec[_i].i)+sqrd(p3->j-lvec[_i].j)+sqrd(p3->k-lvec[_i].k))>_r)
    _r=sqrd(p3->i-lvec[_i].i)+sqrd(p3->j-lvec[_i].j)+sqrd(p3->k-lvec[_i].k);
_r*=scale;
//Calc the side lenght
_r=sqrt(_r*24.00);
//Set the circumscribed tetrahedron.
p0->i=p3->i+_r/2.00;                  //  The description:
p0->j=p3->j-_r/SQRT_3/2.00;           //  the origin of coordinate system is in the center of p0,p1,p2 facet and p3 lay on the positive side of applicate axis
p0->k=p3->k-_r/SQRT_2/SQRT_3/2.00;    //  the abscise axis is oriented parallel with p0,p1 line, the p1 lies on the positive side of abscise axis
p1->i=p3->i-_r/2.00;                  //  the ordinate oxis path throught p1 point that lay on the positive side of ordinate axis
p1->j=p3->j-_r/SQRT_3/2.00;           //  after defining the init sets of coordinates the system should be moved along center vector
p1->k=p3->k-_r/SQRT_2/SQRT_3/2.00;    //  the side lenght of rectlinear tetrahedron is calculated from the circumsphere radius
p2->i=p3->i+0.00;                     //  the circumsphere is builded on the set of given points
p2->j=p3->j+_r/SQRT_3;                //
p2->k=p3->k-_r/SQRT_2/SQRT_3/2.00;    //
p3->i=p3->i+0.00;                     //
p3->j=p3->j+0.00;                     //
p3->k=p3->k+_r*SQRT_3/SQRT_2/2.00;    //
return _r;
}

//This function performs insphere test for 3D case: Is point 'p' located inside the sphere spaned by points 'a', 'b', 'c' and 'd'
char in_sphere(t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d,t_lvec *p)
{
char d4;
void *_p[5];
//Calculate d_abcd
bjsort4((void*)a,(void*)b,(void*)c,(void*)d,&_p[0],&_p[1],&_p[2],&_p[3]); //We do not need to check the sign as we have a common starting configuration
d4=sign_ldet4((double**)_p);
if (bjsort5((void*)&a,(void*)&b,(void*)&c,(void*)&d,(void*)&p,&_p[0],&_p[1],&_p[2],&_p[3],&_p[4])%2)
  return d4==sign_ldet5((double**)_p);
else
  return d4!=sign_ldet5((double**)_p);
}
//This function is like above one but you can save the de4 abcd computations if it is known
char in_sphere_d4(char d4,t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d,t_lvec *p)
{
void *_p[5];
if (bjsort5(a,b,c,d,p,&_p[0],&_p[1],&_p[2],&_p[3],&_p[4])%2)
  return d4==sign_ldet5((double**)_p);
else
  return d4!=sign_ldet5((double**)_p);
}

//This function performs intetrahedron test for 3D case: Is point p inside tetrahedron spaned by points 'a','b','c' and 'd'
char in_tetrahedron(t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d,t_lvec *p)
{
char d4;
double *_p[5];
//Calculate d_abcd
if (bjsort4(a,b,c,d,(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void*)&_p[3])%2)
  d4=-sign_ldet4((double**)_p);
else
  d4=+sign_ldet4((double**)_p);
//Chck all facets keeping clockwise rule
//Check pbdc
if (bjsort4(p,b,c,d,(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)
  {
  if (d4==sign_ldet4((double**)_p))
    return FALSE;
  }
else
  {
  if (d4!=sign_ldet4((double**)_p))
    return FALSE;
  }
//Check pacd
if (bjsort4(p,a,c,d,(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)
  {
  if (d4==sign_ldet4((double**)_p))
    return FALSE;
  }
else
  {
  if (d4!=sign_ldet4((double**)_p))
    return FALSE;
  }
//Check padb
if (bjsort4(p,b,b,d,(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)
  {
  if (d4==sign_ldet4((double**)_p))
    return FALSE;
  }
else
  {
  if (d4!=sign_ldet4((double**)_p))
    return FALSE;
  }
//Check pabc
if (bjsort4(p,a,b,c,(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)
  {
  if (d4==sign_ldet4((double**)_p))
    return FALSE;
  }
else
  {
  if (d4!=sign_ldet4((double**)_p))
    return FALSE;
  }
//OK, the point is inside of all facets
return TRUE;
}



