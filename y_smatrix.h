//This file contains routines to deal with sparse matrices

#define Y_SMATRIX 0x1

#ifndef Y_SYS 
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
#ifndef Y_MATRIX
#include "y_matrix.h"
#endif


//Column-compressed sparse matrix structure
typedef struct{                        //                                      FALSE  ,  TRUE
              unsigned char ctype;     //Type of sparce matrix compression { CCOMPRESS, RCOMPRESS }
              unsigned int ni,nj;      //i-th and j-th dimension of the sparse matrix
              unsigned int nnz;        //Number of nonzero elements in matrix
              unsigned int  *i;        //Pointers to j-column
              unsigned int  *j;        //i-th coordinate of nonzeroelements
              double        *d;        //Nonzero storage
              }t_smatrix;              //Sparse matrix structure

//This function show sparse matrix
void show_smatrix(t_smatrix *sA);

//Allocate memory to smatrix structure
t_smatrix *alloc_smatrix(unsigned char ctype,unsigned int ni,unsigned int nj,unsigned int nnz);

//This function resize sparse matrix structure
char realloc_smatrix(unsigned int new_p,unsigned int new_nnz,t_smatrix *smatrix);

//This function free memory allocated for sparse matrix
void free_smatrix(t_smatrix *sA);

//This function read sparse matrix from file matrix
t_smatrix *read_smatrix(FILE *in);

//This function wtrites sparse matrix into the file
char write_smatrix(FILE *out,t_smatrix *smatrix);

//This function convert dense matrix into sparse matrix
t_smatrix *dmatrix_to_smatrix(register double zero_cutoff,t_dmatrix *dA);
//This function convert sparse matrix into dense matrix
t_dmatrix *smatrix_to_dmatrix(t_smatrix *sA);
//This function convert sparse matrix into triangular dense matrix
t_dmatrix *smatrix_to_tdmatrix(t_smatrix *sA);

//This function sorts entries in a sparse matrix
//Note. This function switch its strategy depending subproblem size
char order_sparse_matrix(t_smatrix *A);

//This function computes transpose of sparse matrix: A^T=Trasnsp(A);
//Note. If 'ordered' flag TRUE the transposed matrix is ordered with 2*o(N^2) complexity else it is unordered with 1*o(N^2) complexity.
t_smatrix *transpose_smatrix(t_smatrix *A);
t_smatrix *transpose_smatrix_potrain(t_smatrix *sA);

//This block of functions do sparse matrix - vector multiplication

//This function calc sparse matrix - dense vector multiplication. Matrix should have row-wise compression
//c=A.b
void multiple_origin_smatrix_dvector(double *c,t_smatrix *sA,double *b);
void multiple_origin_stL_dvector(double *c,t_smatrix *stL,double *b); //The as previous but for L part of symmetric matrix 
//This function calculates c=L.b
void multiple_origin_sL_dvector(double *c,t_smatrix *sL,double *b);
//This function calc sparse matrix - dense vector multiplication. Matrix should have row-wise compression
//c=AT.b
void multiple_transp_smatrix_dvector(double *c,t_smatrix *sAT,double *b);
//This function calculates c=LT.b
void multiple_transp_sL_dvector(double *c,t_smatrix *sL,double *b);

//This function calculates c=L.LT.b product
//NOTE. b vector is destroed during calculations
void multiple_sLLT_dvector(double *c,t_smatrix *sL,double *b);

//This function calculates c=LI.D.LT.b product
//NOTE. b vector is destroed during calculations
void multiple_sLDLT_dvector(double *c,t_smatrix *sLD,double *b);

//These couple functions do linear solvation of diagonal sparse problem

//This function solves row-wise upper triangle sparse problem with unit main diagonal
//NOTE. Resulting vector is replaced input vector
void solve_sIU(t_smatrix *sIU,double *y);
//This function solves row-wise lower triangle sparse problem with unit main diagonal
//NOTE. Resulting vector is replaced input vector
void solve_sIUT(t_smatrix *sIL,double *y);

//This block of functions do sparse matrix multiplication

//This function calculate C=A.A^T
//char *multiple_sparse_AAT(t_smatrix *A,t_smatrix *C);

//This function calculate C=A^T.A
//char *multiple_sparse_ATA(t_smatrix *A,t_smatrix *C);

//This function multiple two sparse matrices C=A.B
t_smatrix *multiple_origin_smatrix_origin_smatrix(t_smatrix *A,t_smatrix *B);

//This function multiple sparse matrices itself as C=A.A^T
t_smatrix *multiple_origin_smatrix_transp_smatrix(t_smatrix *A);

//This function multiple sparse matrices itself as C=A^T.A
t_smatrix *multiple_transp_smatrix_origin_smatrix(t_smatrix *A);

//This block of functions performs Cholesky decomposition of sparse matrix A.

//This function construct L portrain of cholesky decomposed matrix on the base of elimination tree of A
//Jeroen van Grondelle "Symbolic Sparse Cholesky Factorisation Using Elimination Trees"
//Note. it expects sA->j to be row-wise sorted so sA->j[_i]<sA->j[_i+1] if i and i+1 are in th same row
t_smatrix *build_cholesky_portain(t_smatrix *sA);

//Calc sparse cholesky decomposition for portrain
//This function seems build L matrix A=L.LT, might fail due to negative eigenvalues but more precisse then LI.D.LIT decomposition.
char calc_sparse_cholesky_L(t_smatrix *sA,t_smatrix *sL);
//This function seems build LI.D matrix A=LI.D.LIT
void calc_sparse_cholesky_LD(t_smatrix *sA,t_smatrix *sL);

//Cholesky decompositon root function.
//This function require row compessed sparse matrix. It switches on type flag 'l' or 'L' to A=L.LT and 'd' or 'D' to A=LI.D.LIT and NIMPLEMENTED othewise.
t_smatrix *sparse_cholesky_decomposition(char nftype,t_smatrix *sA);


/**************************************** B A C K S U B S T I T U T I O N    B L O C K ****************************************/

//This function do backsubstitution into L.LT matrix
void bsubstitute_sLLT(t_smatrix *sL,double *b);

//This function do backsubstitution into LDLT matrix
void bsubstitute_sLDLT(t_smatrix *sLD,double *b);

//This function solves sparse linear problem from given cholesky decomposition.
//Note. if ((y==x) && (x==b)) than no addition vector storage required, but vector b will be substituted with x data (example: scdsp(L,LT,b,b,temp))
inline char solve_cholesky_decomposed_sparse_problem(t_smatrix *L,t_smatrix *LT,t_vector *b,t_vector *x,t_vector *y);

//This function solve sparse linear system with aid of cholesky decompositin technique. SOLVE[A.x=b] where A sparse matrix, x, b - dense vectors
//IF A is symmetrical matrix (symmetry=TRUE) do: solve[A.x=b] A=LL^T, Lx=L^-T.b, x=L^-1.y (where y=L^-T.b)
//ELSE do: AA=A^T.A, AA=LL^T, c=A^T.b, Lx=L^-T.c, x=L^-1.y (where y=L^-T.c)
char solve_sparse_linear_system_via_cholesky_decomposition(t_smatrix *A,t_vector *x,t_vector *b);

