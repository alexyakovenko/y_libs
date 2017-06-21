//This file contain the routines for vectors manipulations
#include "y_vector.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


//This function print vectro data
void show_vector(register unsigned int size,register double *vec)
{
unsigned int _i;
printf("\n\nvector.size=%d\n",size);
for (_i=0;_i<size;_i++)
  printf("%4.6f ",vec[_i]);
printf("\n");
}

//This function allocate memory for vector
t_vector *alloc_vector(register unsigned int size)
{
register t_vector *vector;
if (!(vector=(t_vector*)malloc(sizeof(t_vector)+sizeof(double)*size)))
  return FALSE;
vector->size=size;
vector->vector=(double*)((void*)vector+sizeof(t_vector));
return vector;
}

//This function resize vector
char resize_vector(t_vector **vector,unsigned int size)
{
if (!((*vector)=(t_vector*)realloc((*vector),sizeof(t_vector)+sizeof(double)*size)))
  return FALSE;
(*vector)->vector=(void*)(*vector)+sizeof(t_vector);
return TRUE;
}

//This function set all vector elements to given value
inline void set_vec(register t_vec *a,register double value)
{
a->i=value, a->j=value, a->k=value;
}
inline void set_vect(register unsigned int n,register double *a,register double value)
{
while (n--) { *a=value, a++; }
}

//This function copy one vector into another. Vectors should not overlap
inline void copy_vec(register t_vec *a,register t_vec *b)
{
a->i=b->i, a->j=b->j, a->k=b->k;
}
inline void copy_vect(register unsigned int n,register double *a,register double *b)
{
a+=n, b+=n; while (n--) { a--, b--, *a=*b; }
}
inline void copy_ivect(register unsigned int n,register double *a,register double *b)
{
a+=n, b+=n; while (n--) { a--, b--, *a=-*b; }
}

//This function read vector from file
char read_vector(register FILE *in,register unsigned int *size,register double **vect)
{
if ( ( (size))&&(fread(size,sizeof(unsigned int),0x1,in)==0x1)&&( (*vect)=(double*)malloc(sizeof(double)*(*size))) )
  {
  if (fread((*vect),sizeof(double),*size,in)==(*size)) return TRUE;
  else free(*vect);
  }
return FALSE;
}

//This function write vector into the hdd
char write_vector(register FILE *out,unsigned int size,register double *vect)
{
return ( (fwrite(&size,sizeof(unsigned int),0x1,out)==0x1)&&(fwrite(vect,sizeof(double),size,out)==size) );
}

//This function imports vector from text file
char import_vector(FILE *in,unsigned int *size,double **vect)
{
unsigned int _i;
extern char buffer[];
if ( (!(size))||(!(vect)) )
  return FALSE;
while ((fgets(buffer,0xFE,in)))
  if ( ( (check_lexem(0x1,buffer)))&&((check_int_type(get_lex(0x1,buffer)))) )
    {
    if (!(*size=(unsigned int)atoi(get_lex(0x1,buffer))))
      return FALSE;
    else
      break;
    }
if ( ( (!(*vect))&&(!((*vect)=(double*)malloc(sizeof(double)*(*size)))) )||( ((*vect))&&(!( (*vect)=(double*)realloc(*vect,sizeof(double)*(*size))) ) ) )
  return FALSE;
_i=0;
while ( (_i<*size)&&((fgets(buffer,0xFE,in))) )
  if ( ( (check_lexem(0x1,buffer)))&&((check_double_type(get_lex(0x1,buffer)))) )
    (*vect)[_i++]=atof(get_lex(0x1,buffer));
if (_i!=*size)
  {
  free(*vect);
  return FALSE;
  }
return TRUE;
}

//This function calc summ of vectors components
inline double summ_vect_components(register unsigned int n,register double *a)
{
register double _summ=0.;
while (n--) { _summ+=*a, a++; }
return _summ;
}

//This function claculate the lenght of vector
inline double calc_vec_norm(register t_vec *a)
{
return a->i*a->i+a->j*a->j+a->k*a->k;
}

//This function claculate the lenght of vector
inline double calc_vect_norm(register unsigned int n,register double *a)
{
register double _norm=0.;
while (n--) { _norm+=*a**a, a++; }
return _norm;
}

//This function multiple all vector coordinates on scalar
inline void multiple_self_vec_scalar(register t_vec *a,register double scalar)
{
a->i*=scalar, a->j*=scalar, a->k*=scalar;
}
inline void multiple_vec_scalar(register t_vec *a,register t_vec *b,register double scalar)
{
a->i=b->i*scalar, a->j=b->j*scalar, a->k=b->k*scalar;
}
inline void multiple_self_vect_scalar(register unsigned int n,register double *a,register double scalar)
{
while (n--) { *a*=scalar, a++; }
}
inline void multiple_vect_scalar(register unsigned int n,register double *a,register double *b,register double scalar)
{
while (n--) { *a=*b*scalar, a++, b++; }
}
//This function invert sign of vectors elements a=-b (like *(-1))
inline void self_vec_inverse_sign(register t_vec *a)
{
a->i=-a->i, a->j=-a->j, a->k=-a->k;
}
inline void self_vect_inverse_sign(register unsigned int n,register double *a)
{
while (n--) { *a=-*a, a++; }
}
inline void vec_inverse_sign(register t_vec *a,register t_vec *b)
{
a->i=-b->i, a->j=-b->j, a->k=-b->k;
}
inline void vect_inverse_sign(register unsigned int n,register double *a,register double *b)
{
while (n--) { *a=-*b, a++, b++; }
}


//This function summ two vectors a+=b.
inline void summ_self_vec(register t_vec *a,register t_vec *b)
{
a->i+=b->i, a->j+=b->j, a->k+=b->k;
}
inline void summ_self_vect(register unsigned int n,register double *a,register double *b)
{
while (n--) { *a+=*b, a++, b++; }
}
inline void summ_vec(register t_vec *a,register t_vec *b,register t_vec *c)
{
a->i=b->i+c->i, a->j=b->j+c->j, a->k=b->k+c->k;
}
inline void summ_vect(register unsigned int n,register double *a,register double *b,register double *c)
{
while (n--) { *a=*b+*c, a++, b++, c++; }
}

//The identical to above
inline void subt_self_vec(register t_vec *a,register t_vec *b)
{
a->i-=b->i, a->j-=b->j, a->k-=b->k;
}
inline void subt_self_vect(register unsigned int n, register double *a, register double *b)
{
while (n--) { *a-=*b, a++, b++; }
}
inline void subt_vec(register t_vec *a,register t_vec *b,register t_vec *c)
{
a->i=b->i-c->i, a->j=b->j-c->j, a->k=b->k-c->k;
}
inline void subt_vect(register unsigned int n, register double *a, register double *b, register double *c)
{
while (n--) { *a=*b-*c, a++, b++, c++; }
}

//This function calculates a=beta*b+c, where a, b and c are vectors and beta is a scalar. It is useful for minimization routines.
inline void vec_mult_summ_vec(register t_vec *a,register double beta,register t_vec *b, register t_vec *c)
{
a->i=beta*b->i+c->i, a->j=beta*b->j+c->j, a->k=beta*b->k+c->k;
}
inline void vect_mult_summ_vect(register unsigned int n,register double *a,register double beta,register double *b, register double *c)
{
while (n--) { *a=beta**b+*c, a++, b++, c++; }
}
//This function calculates a=beta*b-c, where a, b and c are vectors and beta is a scalar. It is useful for minimization routines.
inline void vec_mult_subt_vec(register t_vec *a,register double beta,register t_vec *b, register t_vec *c)
{
a->i=beta*b->i-c->i, a->j=beta*b->j-c->j, a->k=beta*b->k-c->k;
}
inline void vect_mult_subt_vect(register unsigned int n,register double *a,register double beta,register double *b, register double *c)
{
while (n--) { *a=beta**b-*c, a++, b++, c++; }
}

//This function calculates a+=beta*b, where a and b are vectors and beta is a scalar. It is useful for minimization routines.
inline void self_mult_summ_vec(register t_vec *a,register double beta,register t_vec *b)
{
a->i+=beta*b->i, a->j+=beta*b->j, a->k+=beta*b->k;
}
inline void self_mult_summ_vect(register unsigned int n,register double *a,register double beta,register double *b)
{
while (n--) { *a+=beta**b, a++, b++; }
}
//This function calculates a-=beta*b, where a and b are vectors and beta is a scalar. It is useful for minimization routines.
inline void self_mult_subt_vec(register t_vec *a,register double beta,register t_vec *b)
{
a->i-=beta*b->i, a->j-=beta*b->j, a->k-=beta*b->k;
}
inline void self_mult_subt_vect(register unsigned int n,register double *a,register double beta,register double *b)
{
while (n--) { *a-=beta**b, a++, b++; }
}

//This function calculate scalar product of two vectors
inline double calc_vec_vec_scalar_product(register t_vec *a,register t_vec *b)
{
return a->i*b->i+a->j*b->j+a->k*b->k;                                           // 5 flops
}
inline double calc_vect_vect_scalar_product(register unsigned int n,register double *a,register double *b)
{
register double _product=0.;
while (n--) { _product+=*a**b, a++, b++; }
return _product;
}

//This function calc cross product of two vecs axb
inline double calc_vec_vec_cross_product(register t_vec *a,register t_vec *b)
{
return sqrd(a->j*b->k-a->k*b->j)+sqrd(a->k*b->i-a->i*b->k)+sqrd(a->i*b->j-a->j*b->i);  // 15 flops
}
//This function calc cross product of two vecs a=bxc
inline void vec_vec_vmult(register t_vec *a,register t_vec *b,register t_vec *c)
{
a->i=(b->j*c->k-b->k*c->j);                    // 3 flops
a->j=(b->k*c->i-b->i*c->k);                    // 3 flops
a->k=(b->i*c->j-b->j*c->i);                    // 3 flops -> 9 flops total
}

//This function calculates vec vec cosine cos=a.b/(|a||b|)
inline double calc_vec_vec_cos(register t_vec *a,register t_vec *b)
{
return calc_vec_vec_scalar_product(a,b)/sqrt(calc_vec_norm(a)*calc_vec_norm(b));
}
inline double calc_vect_vect_cos(register unsigned int n,register double *a,register double *b)
{
return calc_vect_vect_scalar_product(n,a,b)/sqrt(calc_vect_norm(n,a)*calc_vect_norm(n,b));
}

//---------------------------------      Q U A T E R N I O N      P A R T      ----------------------------------------

//Thi function summs couple quaternions qa=qb+qc
//Note qa and qb and qc might be the same
inline void quaternion_summ(t_quaternion *qa,t_quaternion *qb,t_quaternion *qc)
{
qa->i=qb->i+qc->i;
qa->j=qb->j+qc->j;
qa->k=qb->k+qc->k;
qa->w=qb->w+qc->w;
}

//this function multiply two quaternions qa=qb*qc
inline void quaternion_mult(t_quaternion *qa,t_quaternion *qb,t_quaternion *qc)
{
double e, f, g, h;

e=(qb->i+qb->k)*(qc->i+qc->j); // 3 flops
f=(qb->i-qb->k)*(qc->i-qc->j);
g=(qb->w+qb->j)*(qc->w-qc->k);
h=(qb->w-qb->j)*(qc->w+qc->k);

qa->w=(qb->k-qb->j)*(qc->j-qc->k)-(e+f-g-h)*.5; // 8 flops
qa->i=(qb->w+qb->i)*(qc->w+qc->i)-(e+f+g+h)*.5;
qa->j=(qb->w-qb->i)*(qc->j+qc->k)+(e-f+g-h)*.5;
qa->k=(qb->j+qb->k)*(qc->w-qc->i)+(e-f-g+h)*.5;
//total 44 flops
}

//This function normalizes quaternion
inline void normalize_quaternion(t_quaternion *q)
{
double qq;
if ((qq=q->i*q->i+q->j*q->j+q->k*q->k+q->w*q->w)<TINY) { q->i=q->j=q->k=0., q->w=1.; }
else                                      { qq=sqrt(qq); q->i/=qq, q->j/=qq, q->k/=qq, q->w/=qq; }
}

//This function calculates quaternion x vector x conjugate_quaternion product r=q.r.-q
inline void calc_uquaterion_vec_cuquaternion_product(t_vec *r,t_quaternion *q,t_vec *ro)
{
t_quaternion qt1, qt2;
qt1.i=ro->i, qt1.j=ro->j, qt1.k=ro->k, qt1.w=0.;
quaternion_mult(&qt2,q,&qt1);
q->i=-q->i, q->j=-q->j, q->k=-q->k;
quaternion_mult(&qt1,&qt2,q);
q->i=-q->i, q->j=-q->j, q->k=-q->k;
r->i=qt1.i, r->j=qt1.j, r->k=qt1.k;
}



















