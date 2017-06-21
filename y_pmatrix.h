#define Y_PMATRIX 0x1
#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_PLIBS
#include "y_plibs.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif

#define MAX_PQANT_pcalc_vect_norm 0x1000
#define MIN_PQANT_pcalc_vect_norm 0x400
#define MIN_PQANT_pcalc_vect_scalar_product 0x100

//This function calculates vector lenght in parallel
inline char pcalc_vect_norm(double *norm,unsigned int size,double *vect,unsigned int pjob_id,t_ypsystem *ypsys);

//This function calculates scaled vector lenght in parallel
inline char pcalc_vect_scaled_norm(double *norm,double scale,unsigned int size,double *vect,unsigned int pjob_id,t_ypsystem *ypsys);
void get_stack_pcalc_vect_scaled_norm(unsigned int size,unsigned int n_threads,unsigned int *n_hpcards,size_t *hpcards_size,unsigned int *n_pcards,size_t *pcards_size);

//This function calculates summ of vector's unsigned elements in parallel
inline char pcalc_vect_unsigned_summ(double *summ,unsigned int size,double *vect,unsigned int pjob_id,t_ypsystem *ypsys);
void get_stack_pcalc_vect_unsigned_summ(unsigned int size,unsigned int n_threads,unsigned int *n_hpcards,size_t *hpcards_size,unsigned int *n_pcards,size_t *pcards_size);

//This function multiples vector on scalar in parallel
unsigned int _pcalc_vect_scalar_product(char *pcard,void *vp);
inline char pcalc_vect_scalar_product(double scalar,unsigned int size,double *vect,unsigned int pjob_id,t_ypsystem *ypsys);

//This function calculates vectors dot product in parallel
unsigned int _pcalc_vect_vect_dot_product(char *pcard,void *vp);
inline char pcalc_vect_vect_dot_product(double *product,unsigned int size,double *vect_a,double *vect_b,unsigned int pjob_id,t_ypsystem *ypsys);

#define MAX_PQUANT_NI_pQL_factorization 0xF
#define MAX_PQUANT_NJ_pQL_factorization 0xFF
#define MIN_PQUANT_NI_pQL_factorization 0x1
#define MIN_PQUANT_NJ_pQL_factorization 0xF


