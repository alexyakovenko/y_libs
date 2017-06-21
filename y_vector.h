#define Y_VECTOR 0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_FILE
#include "y_file.h"
#endif
#ifndef Y_TXT
#include "y_txt.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif


//This structure defines 3D coordinate structure for geometry
typedef struct{unsigned int i,j,k;}t_len; //discrette lenght structure
typedef struct{double i,j,k;}t_vec;    //the normal points radius-vector
typedef struct{double i,j,k,w;}t_lvec; //the lifted points radius-vector
typedef t_lvec t_quaternion; //The quaternion

//This structure defines X-dimenssion structure for linear algebra
typedef struct{
               unsigned int size;
               double *vector;
              }t_vector;

typedef union { t_vec v; double d[3]; }t_vecd; //Union version of vector massive

///This function print vectro data
#define show_vec(vec)    show_vector(0x3,(double*)vec)
#define show_vect(vect)  show_vector(vect->size,vect->vector)
void show_vector(unsigned int size,double *vect);

//This function allocate memory for vector
t_vector *alloc_vector(register unsigned int size);

//This function resize vector
char resize_vector(t_vector **vector,unsigned int size);

//This function set all vector elements to given value
inline void set_vec(t_vec *a,double value);
inline void set_vect(register unsigned int size,register double *vect,register double value);

//This function copy one vector into another. Vectors should not overlap
inline void copy_vec(t_vec *a,t_vec *b);
inline void copy_vect(unsigned int n,double *a,double *b);
inline void copy_ivect(unsigned int n,double *a,double *b);

//This function read vector from file
char read_vector(register FILE *in,register unsigned int *size,register double **vect);

//This function write vector into the hdd
#define write_vec(out,vec)   write_vector(out,0x3,(double*)vec)
#define write_vect(out,vect) write_vector(out,vect->size,vect->vector)
char write_vector(register FILE *out,unsigned int size,register double *vect);

//This function imports vector from text file
char import_vector(FILE *in,unsigned int *size,double **vect);

//This function calc summ of vectors components
inline double summ_vect_components(register unsigned int n,register double *a);

//This function claculate the lenght of vector
inline double calc_vec_norm(register t_vec *a);

//This function claculate the lenght of vector
inline double calc_vect_norm(register unsigned int n,register double *a);

//This function multiple all vector coordinates on scalar
inline void multiple_self_vec_scalar(register t_vec *a,register double scalar);
inline void multiple_vec_scalar(register t_vec *a,register t_vec *b,register double scalar);
inline void multiple_self_vect_scalar(register unsigned int n,register double *a,register double scalar);
inline void multiple_vect_scalar(register unsigned int n,register double *a,register double *b,register double scalar);
//This function invert sign of vectors elements a=-b (like *(-1))
inline void self_vec_inverse_sign(register t_vec *a);
inline void self_vect_inverse_sign(register unsigned int n,register double *a);
inline void vec_inverse_sign(register t_vec *a,register t_vec *b);
inline void vect_inverse_sign(register unsigned int n,register double *a,register double *b);

//This function summ two vectors a+=b.
inline void summ_self_vec(register t_vec *a,register t_vec *b);
inline void summ_self_vect(register unsigned int n,register double *a,register double *b);
//This function summ two vectors a=b+c.
inline void summ_vec(register t_vec *a,register t_vec *b,register t_vec *c);
inline void summ_vect(register unsigned int n,register double *a,register double *b,register double *c);

//The identical to above
inline void subt_self_vec(register t_vec *a,register t_vec *b);
inline void subt_self_vect(register unsigned int n, register double *a, register double *b);
inline void subt_vec(register t_vec *a,register t_vec *b,register t_vec *c);
inline void subt_vect(register unsigned int n, register double *a, register double *b, register double *c);

//This function calculates a=beta*b+c, where a, b and c are vectors and beta is a scalar. It is useful for minimization routines.
inline void vec_mult_summ_vec(register t_vec *a,register double beta,register t_vec *b, register t_vec *c);
inline void vect_mult_summ_vect(register unsigned int n,register double *a,register double beta,register double *b, register double *c);
//This function calculates a=beta*b-c, where a, b and c are vectors and beta is a scalar. It is useful for minimization routines.
inline void vec_mult_subt_vec(register t_vec *a,register double beta,register t_vec *b, register t_vec *c);
inline void vect_mult_subt_vect(register unsigned int n,register double *a,register double beta,register double *b, register double *c);

//This function calculates a+=beta*b, where a and b are vectors and beta is a scalar. It is useful for minimization routines.
inline void self_mult_summ_vec(register t_vec *a,register double beta,register t_vec *b);
inline void self_mult_summ_vect(register unsigned int n,register double *a,register double beta,register double *b);
//This function calculates a-=beta*b, where a and b are vectors and beta is a scalar. It is useful for minimization routines.
inline void self_mult_subt_vec(register t_vec *a,register double beta,register t_vec *b);
inline void self_mult_subt_vect(register unsigned int n,register double *a,register double beta,register double *b);

//This function calculate product of two vectors
inline double calc_vec_vec_scalar_product(register t_vec *a,register t_vec *b);
inline double calc_vect_vect_scalar_product(register unsigned int n,register double *a,register double *b);

//This function calc cross product of two vecs axb
inline double calc_vec_vec_cross_product(register t_vec *a,register t_vec *b);

//This function calc cross product of two vecs a=bxc
inline void vec_vec_vmult(t_vec *a,t_vec *b,t_vec *c);

//This function calculates vec vec cosine cos=a.b/(|a||b|)
inline double calc_vec_vec_cos(register t_vec *a,register t_vec *b);
inline double calc_vect_vect_cos(register unsigned int n,register double *a,register double *b);

//---------------------------------      Q U A T E R N I O N      P A R T      ----------------------------------------

//Thi function summs couple quaternions qa=qb+qc
//Note qa and qb and qc might be the same
inline void quaternion_summ(t_quaternion *qa,t_quaternion *qb,t_quaternion *qc);

//this function multiply two quaternions qa=qb*qc
inline void  quaternion_mult(t_quaternion *qa,t_quaternion *qb,t_quaternion *qc);

//This function normalizes quaternion
inline void normalize_quaternion(t_quaternion *q);

//This function calculates quaternion x vector x conjugate_quaternion product r=q.r.-q
inline void calc_uquaterion_vec_cuquaternion_product(t_vec *r,t_quaternion *q,t_vec *ro);



