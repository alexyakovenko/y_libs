#define Y_MATRIX 0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_TXT
#include "y_txt.h"
#endif
#ifndef Y_FILE
#include "y_file.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif
#ifndef Y_LIST   
#include "y_list.h"
#endif
#ifndef Y_VECTOR
#include "y_vector.h"
#endif

//NOTE. We assume to mark matrices with upercase letters and vectors with lowercase letters

//The matrix structure
typedef struct{
               unsigned int ni,nj;
               double **d;
              }t_dmatrix;

//Tensors
typedef double t_tensor [0x3][0x3];
typedef double t_qtensor[0x4][0x4];

typedef union {  t_tensor R; double d[9];  } t_tensord;
typedef union { t_qtensor Q; double d[16]; } t_qtensord;


//kronekker symbol
#define KRON(i,j) (unsigned int)(i==j)

#define DOMINANT_DIAG_ITERATIONS 0xF 

//This function print the matrix
void show_matrix(unsigned int ni,unsigned int nj,double **A);

//These functions allocate memory for double-pointer array
inline double **alloc_marray(register unsigned int ni,register unsigned int nj);
inline double **calloc_marray(register unsigned int ni,register unsigned int nj);
//This function allocates memory for the matrix
t_dmatrix *alloc_dmatrix(register unsigned int ni,register unsigned int nj);
//This function allocates memory for the dense triangle matrix, type means 'L' or 'U' matrix to alloc
t_dmatrix *alloc_tdmatrix(char type,register unsigned int ni);
//This function reads dense matrix from the input file
t_dmatrix *read_dmatrix(FILE *in);
//This function reads dense triangle matrix from the input file
t_dmatrix *read_tdmatrix(FILE *in,char *type);
//This function write dense matrix to hdd
char write_dmatrix(FILE *out,t_dmatrix *dA);
//This function write dense triangle matrix to hdd
char write_tdmatrix(FILE *out,char type,t_dmatrix *tdA);

//This function setup dense matrix
inline void set_dmatrix(t_dmatrix *A,register double value);
//This function sets matrix to be identity
void set_identity_dmatrix(unsigned int n,unsigned int m,double **A);

//This function transpose matrix
//NOTE. It is possiable to use the same matrix as both arguments. 
inline void transpose_dmatrix(unsigned int n,unsigned int m,double **dA,double **dB);
//This function summarise two matrices
//NOTE. A, B and C can be the same
inline void summ_matrix(unsigned int n,unsigned int m,double **dA,double **dB,double **dC);

//This function calculates DU.b=c.
//Note. U is symmetrical matrix stored in upper triangle matrix for memory saving purposes.
void multiple_origin_tdmatrixU_origin_vector(unsigned int n,register  double *c,double **tdU,register double *b);

//
//This function multiple matrix and vector: c = Ab
inline void multiple_origin_matrix_origin_vector(unsigned int ni, unsigned int nj,double **A,double *b,double *c);
//                                              T
//This function get vector matrix product: c = A b
inline void multiple_transp_matrix_origin_vector(unsigned int ni,unsigned int nj,double **A,double *b,register double *c);
//                                          T   T
//This function get vector matrix product: c = b A
inline void multiple_transp_vector_origin_matrix(unsigned int ni,unsigned int nj,double **A,double *b,double *c);
//                                          T   T T
//This function get vector matrix product: c = b A
inline void multiple_transp_vector_transp_matrix(unsigned int ni,unsigned int nj,double **A,double *b,double *c);
//                                                   T
//This function  get vector matrix product: C = a * b . Result expands into matrix
inline void multiple_origin_vector_transp_vector(unsigned int ni,unsigned int nj,double *a,double *b,double **C);

//This function perform multiple matrix C = AB
inline void multiple_origin_matrix_origin_matrix(unsigned int ni,unsigned int nj,unsigned int nk,double **A,double **B,double **C);
//                                           T
//This function perform multiple matrix C = A B
inline void multiple_transp_matrix_origin_matrix(unsigned int ni,unsigned int nj,unsigned int nk,double **A,double **B,double **C);
//                                            T
//This function perform multiple matrix C = AB
inline void multiple_origin_matrix_transp_matrix(unsigned int ni,unsigned int nj,unsigned int nk,double **A,double **B,double **C);
//                                           T T
//This function perform multiple matrix C = A B
inline void multiple_transp_matrix_transp_matrix(unsigned int ni,unsigned int nj,unsigned int nk,double **A,double **B,double **C);

//This function multiple matrix scalar
//Note A and B can be the same
char multiple_matrix_scalar(t_dmatrix *A,t_dmatrix *B,register double _f);
//This function compute matrix row norm
double get_row_norm(unsigned int _i,unsigned int m,double **dA);

/******************************************** L U   d e c o m p o s i t i o n ****************************************************/
/*
  LU-decomposition according to Crout's algorithm with pivoting.
Description:
   char LU_decomposition(int n,double **matrix,int *p_order,char *twoness,double *temp);
Parameters:
   matrix - source matrix (n x n) on input, destination on output;
   n - matrix size;
   p_order - integer array (size n) to remember permutations;
   twoness - on output, contains +1 or -1 for even or odd permutations number.
   temp - temporary array (size n).
   Returns:
   0 - the source matrix is singular (invalid for decomposition),
   1 - if OK.
  Back substitution, using LU decomposed matrix.
Description:
  void LU_back_substitution(double **matrix,int n,int *p_order,double *result_vector);
Parameters:
  a - the matrix decomposed by Crout;
  n - the matrix size;
  indx - permutation order obtained by decomposition algorithm;
  result_vector - the vector (size n) to be substituted on input, the result
                  of the substitution on output.
  Note: matrix and p_order are not modified by this routine and could be
  used in multiple calls.
  Invertation of matrix, using LU decomposed matrix.
Description:
  void LU_invert_dmatrix(double **matrix,int n,int *p_order,double **invert_dmatrix,double *temp);
Parameters:
  matrix - the matrix decomposed by Crout;
  n - the matrix size;
  p_order - permutation order obtained by decomposition algorithm;
  invert_dmatrix - the destination matrix;
  temp - temporary array (size n).
  Note: test for singularity has been already obtained on the
  matrix decomposition, a and indx are not modified by this routine,
  the routine uses multiple backsubstitutions (previous algorithm).
  Obtaining the matrix determinant, using LU-decomposed matrix
Description:
  double LU_determinant(double **matrix,int n,int *p_order,char *twoness);

Parameters:
  matrix - the matrix decomposed by Crout;
  n - the matrix size;
  p_order - permutation order obtained by decomposition algorithm;
  twoness - the parity sign (+1 or -1) obtained at decomposition.

Returns:
  the determinant value. Note: non-zero (the matrix cannot be
  singular, if decomposed properly); matrix, p_order and twoness are not modified
  by this routine.
*/
/* the decomposition itself */
char LU_decomposition (unsigned int n, double **matrix, unsigned int *p_order, char *twoness, double *temp);
// the back substitution
void LU_back_substitution (unsigned int n,double **A,unsigned int *p_order,double *y);
// the matrix invertation
void LU_invert_dmatrix (unsigned int n,double **A,unsigned int *p_order,double **IA,double *_t);
// calculating determinant
double LU_determinant (unsigned int n,double **A,unsigned int *p_order,char twoness);

/*************************************** G A U S S   E l i m i n a t i o n ******************************************/

//This function solves linear system of equation with gauss-elimination method. It is preffered if the matrix will be used only once.
unsigned int gauss_solve_dmatrix(unsigned int n, double **A,double *b,double badly_tol);

/******************************************* S V D ******************************************************************/

// Bidiagonalization. It produces either super_diag or sub_diag
double Hausholder_bidiagonalization (double *diag,double *sdiag,unsigned int ni,unsigned int nj,unsigned int nr,double **A,double **U,double **V);
//This function diagonalize bidiagonal matrix
void diagonalize(unsigned int ni,unsigned int nj,unsigned int nr,double *super_diag,double *diag,double **U,double **V,double epsilon);

//This function compute the SVD. super_diag is the help massive for Hausholder Transformation
//A[m:n], U[m:m], V[[n:n], diag[m], super_diag[m]; m<=n
void singular_value_decomposition(unsigned int ni,unsigned int nj,unsigned int nr,double **A,double **B,double **U,double *diag,double *super_diag,double **V,double epsilon);

//This function calculates svd using previouw one and some other twics:
void svd(char symm,char reorder,unsigned int ni,unsigned int nj,unsigned int nr,double **A,double **B,double **U,double *s,double **V,double *temp);

//This is deprecated algorithm of Jacobi rotations of symmetriñ matrix. Please use svd instead of it anywhere you can!
//Note. This method performs the diagonalization only. No any eigen vectors are calculated. Diagonalized A is S. A can be restored from untouched lower triangle.
char diagonalization_Jacobi(unsigned int n_iter,t_dmatrix *A);

//This function performs Grahamm-Smidt ortonormalization
inline void gram_schmidt_ortonormalization(unsigned int ni,unsigned int nj,double **AT,double *_t);

//------------------- T H E   T R I A N G L E   M A T R I C E S   P A R T ---------------------------------

//This function generates triangle matrix form square self-multiplicaion

//This function calculates C=A.D.AT
inline void multiple_nXDXT(unsigned int ni,unsigned int nj,double **C,register double **A,double *d);
inline void multiple_tXDXT(unsigned int ni,unsigned int nj,double **U,register double **A,double *d);
//This function calculates C=AT.D.A
inline void multiple_nXTDX(unsigned int ni,unsigned int nj,double **C,register double **A,double *d);
inline void multiple_tXTDX(unsigned int ni,unsigned int nj,double **U,register double **A,double *d);
//This function calculates C=A.AT
inline void multiple_nXXT(unsigned int ni,unsigned int nj,double **C,register double **A);
inline void multiple_tXXT(unsigned int ni,unsigned int nj,double **U,register double **A);
//This function calculates C=AT.A
inline void multiple_nXTX(unsigned int ni,unsigned int nj,double **C,register double **A);
inline void multiple_tXTX(unsigned int ni,unsigned int nj,double **U,register double **A);

//This function perform modified Choletski factorization (Square root method)
//Note. The matrix is given here in form of lower triangle of simmetric A. Result define as L.D.LT=P.A.P so to minimize ||E|| in A-E=L.D.LT, where A-E -s positively-defined well-conditioned matrix.
// mu - relaxation factor - recomended to set as 0.1
//Algorithm taken from Haw-ren Fang, Dianne P. O'Leary Modified cholesky algorithms: a catalog with new approaches 2006
//Note this function returns FALSE if exact value was calculated (i.e. ||E||=0. ) and TRUE otherwise 
char modified_cholesky_factorization_dmatrix(double mu,unsigned int n,unsigned int *p,double **A);

//This function perform Choletski transformation (Square root method)
//Note. The matrix is given here in form of upper triangle of simmetric A. 
//IF type=='U' the result define as UT.D.U ELSE IF type=='L' the result is L.D.LT ELSE YERROR_NIMPLEMENT
void cholesky_decomposition_tdmatrix(char type,unsigned int n,double **A);

//This function calculates L.b=c.
//Note. It suggests that L has unit diagonal.
void multiple_tdILmatrix_vector(unsigned int n,register  double *c,double **IL,register double *b);
//This function calculates U.b=c.
//Note. U is symmetrical matrix stored in upper triangle matrix for memory saving purposes.
void multiple_tdmatrixU_vector(unsigned int n,register  double *c,double **tU,register double *b);
//This function calculates U.b=c.
//Note. It suggests that U has unit diagonal.
void multiple_tdIUmatrix_vector(unsigned int n,register  double *c,double **IU,register double *b);
void multiple_tdIUTmatrix_vector(unsigned int n,register double *c,double **IU,register double *b);
//This function invert triangle matrix U
inline void invert_tdIUmatrix(unsigned int n,double **IIU,double **IU);
//This function performs invertion of lower triangle matrix Calc c=L^-1.b

//This function performs invertion of upper triangle matrix Calc y=U^-1.y
inline void bsubstitute_tdU(unsigned int n,double **U,double *y);
//This function performs invertion of transposed upper triangle matrix Calc y=U^-T.y
inline void bsubstitute_tdUT(unsigned int n,double **U,double *y);


//This function performs solvation of lower triangle linear equations system y=L^-1.y
inline void bsubstitute_tdIL(unsigned int n,double **IL,double *y);
//This function performs solvation of upper triangle linear equations system y=U^-1.y
//Note. The matrix U diagonal suggested to be unit
inline void bsubstitute_tdIU(unsigned int n,double **IU,double *y);
//This function performs solvation of transposed lower triangle linear equations system y=L^-1.y
inline void bsubstitute_tdILT(unsigned int n,double **IL,double *y);
//This function performs solvation of transposed upper triangle matrix y=U^-T.y
//Note. The matrix U diagonal suggested to be unit
inline void bsubstitute_tdIUT(unsigned int n,double **IU,register double *y);
//This function solves equation (L.D.LT).y=y
inline void bsubstitute_tdLDLT(unsigned int n,double **tdLD,double *y);
//This function solves equation (UT.D.U).y=y
inline void bsubstitute_tdUTDU(unsigned int n,double **tdDU,double *y);

//This function bsubstitude tdL matrix with given pivot 
inline void bsubstitute_pivoted_tdLDLT(unsigned int n,unsigned int *p,double **LD,double *y,double *_t);
//This function bsubstitude tdU matrix with given pivot 
inline void bsubstitute_pivoted_tdUTDU(unsigned int n,unsigned int *p,double **UD,double *y,double *_t);

//This function inverts cholesky decomposed factor marix
void invert_tdUTDU(unsigned int n,double **IA,double **A,double *_t);

//------------------- T H E   3 D   P A R T   ( t e n s o r s   p a r t ) --------------------------------

//Calculate tensor determinant
#define TENSOR_DET(T) calc_det3x3(T[0][0],T[0][1],T[0][2], \
                                  T[1][0],T[1][1],T[1][2], \
                                  T[2][0],T[2][1],T[2][2])

//This function transposes tensor A=B^T (useful for calculation of inverse rotation matrix)
//Note A and B can be the same
inline void transpose_tensor(t_tensor *A, t_tensor *B);
//This function calculates product of b^T.R.a
inline double calc_bTRa(register t_vec *b,register t_tensor *R,register t_vec *a);
//This function calculates tensor-vcector product a=T.b
inline void multiple_origin_tensor_origin_vec(register t_vec *a,register t_tensor *T,register t_vec *b);
//This function calculates tensor-vcector product a=T^T.b
inline void multiple_transp_tensor_origin_vec(register t_vec *a,register t_tensor *T,register t_vec *b);
//This function calculates tensor-vcector product aT=bT.R
inline void multiple_transp_vec_transp_tensor(register t_vec *a,register t_tensor *T,register t_vec *b);
//This function calculates tensor.tensor product A=B.C
inline void multiple_origin_tensor_origin_tensor(register t_tensor *A,register t_tensor *B,register t_tensor *C);
//This function calculates tensor.tensor^T product A=B.C^T
inline void multiple_origin_tensor_transp_tensor(register t_tensor *A,register t_tensor *B,register t_tensor *C);

//This function solves 3x3 linear equations system with double pivoting
char gauss_solve_tensor(t_tensor *A,t_vec *b);

//This is univewrsializetion cover of jacobi Singular Value Decomposition
//NOTE. T sould have the lowest direction. T can be the same with A, but only if A symetric square matrix. In this case A will be destroed!
//It uses the rule that A^TA = V . S^2.V^T  and AA^T= U.S^2.U^T
unsigned int tensor_svd(char reorder,char eigenvectors,char symmetry_check,double tol,t_tensor *A,t_tensor *U,t_vec *S,t_tensor *V);
//This function makes eighten-values positive
inline void positive_svd(t_tensor *U,t_vec *S,t_tensor *V);

//------------------- T H E   4 D   P A R T   ( q t e n s o r s   p a r t ) --------------------------------

//This function multiple two qutensors into 3D tensor
inline void multiple_qtensor_tqtensor(t_tensor *t,t_qtensor *qt,t_qtensor *qtt);

//This function solves 4x4 linear equations system with iterative Gauss-Seindel method
char gauss_seindel_solve_qtensor(t_lvec *x,t_qtensor *Q,t_lvec *b,unsigned int niter);

//---------------------- S O M E    A C C E L E R A T I O N    P A R T --------------------------------------

//This function calculates 2x2 determinant
inline double calc_det2x2(double v00,double v01,
                          double v10,double v11);
//This function calculates 3x3 determinant
inline double calc_det3x3(double v00,double v01,double v02,
                          double v10,double v11,double v12,
                          double v20,double v21,double v22 );
//This function calculates 4x4 determinant
inline double calc_det4x4(double v00,double v01,double v02,double v03,
                          double v10,double v11,double v12,double v13,
                          double v20,double v21,double v22,double v23,
                          double v30,double v31,double v32,double v33 );
//This function calculates 5x5 determinant
inline double calc_det5x5(double v00,double v01,double v02,double v03,double v04,
                          double v10,double v11,double v12,double v13,double v14,
                          double v20,double v21,double v22,double v23,double v24,
                          double v30,double v31,double v32,double v33,double v34,
                          double v40,double v41,double v42,double v43,double v44);
//This function calculates 6x6 determinant
inline double calc_det6x6(double v00,double v01,double v02,double v03,double v04,double v05,
                          double v10,double v11,double v12,double v13,double v14,double v15,
                          double v20,double v21,double v22,double v23,double v24,double v25,
                          double v30,double v31,double v32,double v33,double v34,double v35,
                          double v40,double v41,double v42,double v43,double v44,double v45,
                          double v50,double v51,double v52,double v53,double v54,double v55);
//We need 6x6 determinant to calculate 3D RT, but all future should be processed with more efficient techniques.
