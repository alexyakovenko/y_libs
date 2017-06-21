#include "y_matrix.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


//This function print the matrix
void show_matrix(unsigned int ni,unsigned int nj,double **A)
{
unsigned int i,j;
printf("\n\nt_dmatrix:\n matrix_ni=%d, matrix_nj=%d;\n",ni,nj);
for (i=0;i<ni;i++)
  {
  for (j=0;j<nj;j++)
    printf("%8.5f ",A[i][j]);
  printf("\n");
  }
}


//These functions allocate memory for double-pointer array
inline double **alloc_marray(register unsigned int ni,register unsigned int nj)
{
extern unsigned int ylib_errno;
register unsigned int _i;
register double **_d;
if (!(_d=(double**)malloc(ni*(sizeof(double*)+sizeof(double)*nj)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
for (*_d=(double*)((void*)_d+sizeof(double*)*ni), _i=1; _i<ni; _i++) _d[_i]=_d[_i-1]+nj;
return _d;
}
inline double **calloc_marray(register unsigned int ni,register unsigned int nj)
{
extern unsigned int ylib_errno;
register unsigned int _i;
register double **_d;
if (!(_d=(double**)calloc(0x1,ni*(sizeof(double*)+sizeof(double)*nj)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
for (*_d=(double*)((void*)_d+sizeof(double*)*ni), _i=1; _i<ni; _i++) _d[_i]=_d[_i-1]+nj;
return _d;
}
//This function allocate memory for the matrix
t_dmatrix *alloc_dmatrix(register unsigned int ni,register unsigned int nj)
{
t_dmatrix *A;
if (!(A=(t_dmatrix*)malloc(sizeof(double)*ni*nj+sizeof(double*)*ni+sizeof(t_dmatrix))))
  return FALSE;
A->ni=ni;
A->nj=nj;
A->d=(double**)((void*)A+sizeof(t_dmatrix));
*A->d=(void*)A->d+sizeof(double*)*A->ni;
for (ni=0x1;ni<A->ni;ni++)
  A->d[ni]=A->d[ni-0x1]+nj;
return A;
}
//This function allocate memory for the dense triangle matrix
t_dmatrix *alloc_tdmatrix(char type,register unsigned int ni)
{
t_dmatrix *A;

if ( (type!='L')&&(type!='l')&&(type!='U')&&(type!='u') ) { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
if (!(A=(t_dmatrix*)malloc(sizeof(double)*ni*(ni+0x1)/0x2+sizeof(double*)*ni+sizeof(t_dmatrix)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
A->ni=A->nj=ni;
A->d=(double**)((void*)A+sizeof(t_dmatrix));
*A->d=(double*)((void*)A->d+sizeof(double*)*A->ni);
if ( (type=='L')||(type=='l') )
  for (ni=1;ni<A->ni;ni++)
    A->d[ni]=(double*)((void*)A->d[ni-1]+sizeof(double)*ni);
else 
  for (ni=1;ni<A->ni;ni++)
    A->d[ni]=(double*)((void*)A->d[ni-1]+sizeof(double)*(A->ni-ni));
return A;
}

//This function reads dense matrix from the input file
t_dmatrix *read_dmatrix(FILE *in)
{
register t_dmatrix *dA;
unsigned int dim[0x2];
if ( (fread(dim,sizeof(unsigned int),0x2,in)==0x2)&&((dA=(t_dmatrix*)alloc_dmatrix(dim[0x0],dim[0x1]))) )
  {
  if  (fread(*dA->d,sizeof(double),dim[0x0]*dim[0x1],in)==dim[0x0]*dim[0x1]) return dA;
  else free(dA);
  }
return FALSE;
}
//This function reads dense triangle matrix from the input file
t_dmatrix *read_tdmatrix(FILE *in,char *type)
{
register t_dmatrix *tdA;
unsigned int dim[2];
if ( (fread(dim,sizeof(unsigned int),0x2,in)==0x2)&&(fread(type,sizeof(char),0x1,in)==0x1) ) 
  {
  if ( (dim[0]==dim[1])&&( (*type=='u')||(*type=='U')||(*type=='l')||(*type=='L') ) )
    {
    if ( (tdA=(t_dmatrix*)alloc_tdmatrix(*type,dim[0])))
      {
      if  (fread(*tdA->d,sizeof(double),dim[0]*(dim[1]+1)/2,in)==dim[0]*(dim[1]+1)/2) return tdA;
      else { free(tdA); goto ERROR_IO; }
      }
    else ylib_errno=YERROR_MEMORY;
    }
  else ylib_errno=YERROR_DATA_CONSISTMENT;
  }
else { ERROR_IO: ylib_errno=YERROR_IO; }
return FALSE;
}
//This function write dense matrix to hdd
char write_dmatrix(FILE *out,t_dmatrix *dA)
{
return ( (fwrite(&dA->ni,sizeof(unsigned int),0x1,out)==0x1)&&(fwrite(&dA->nj,sizeof(unsigned int),0x1,out)==0x1)&&(fwrite(*dA->d,sizeof(double),dA->ni*dA->nj,out)==dA->ni*dA->nj) );
}
//This function write dense triangle matrix to hdd
char write_tdmatrix(FILE *out,char type,t_dmatrix *tdA)
{
if ( (fwrite(&tdA->ni,sizeof(unsigned int),0x1,out)==0x1)&&(fwrite(&tdA->nj,sizeof(unsigned int),0x1,out)==0x1)&&
     (fwrite(&type,sizeof(char),0x1,out)==0x1)&&(fwrite(tdA->d,sizeof(double),tdA->ni*(tdA->ni+1)/2,out)==tdA->ni*(tdA->ni+1)/2) ) return TRUE;
ylib_errno=YERROR_IO; return FALSE;
}

//This function setup dense matrix
inline void set_dmatrix(t_dmatrix *A,register double value)
{
register unsigned int _i=A->ni*A->nj;
register double *_d=*A->d;
while (_i--)
  {
  *_d=value;
   _d++;
  }
}

//This function sets matrix to be identity
void set_identity_dmatrix(unsigned int n,unsigned int m,double **A)
{
register unsigned int _i, _j;
register double *row;
_i=n; while (_i--) { _j=m, row=A[_i]+m; while (_j--) { row--, *row=0.; } }
_i=( n > m ) ? n : m; while (_i--) A[_i][_i]=1.;
}


//This function transpose matrix
//NOTE. It is possiable to use the same matrix as both arguments. 
inline void transpose_dmatrix(unsigned int n,unsigned int m,double **dA,double **dB)
{
register unsigned int _j;
register double _d;

if (dA!=dB) while (n--) { _j=m; while (_j--) dB[_j][n]=dA[n][_j]; } //Just copy it
else 
  {
  if (n!=m)
    {
    if (n>m) do{ _j=m; while (_j--) dA[_j][n]=dA[n][_j]; }while(--n!=m);
    else     do{ _j=m; while (_j--) dA[_j][n]=dA[n][_j]; }while(n!=--m);    
    }
  while (n--) { _j=n; while (--_j) { _d=dB[_j][n], dB[_j][n]=dA[n][_j], dA[n][_j]=_d; } } //Skip diagonal
  }
}

//This function summarise two matrices
//NOTE. A, B and C can be the same
inline void summ_matrix(unsigned int n,unsigned int m,double **dA,double **dB,double **dC)
{
register unsigned int _j;
while (n--) { _j=m; while (_j--) { dC[n][_j]=dA[n][_j]+dB[n][_j]; } }
}



//This function calculates DU.b=c.
//Note. U is symmetrical matrix stored in upper triangle matrix for memory saving purposes.
void multiple_origin_tdmatrixU_origin_vector(unsigned int n,register  double *c,double **tdU,register double *b)
{
unsigned int _i, _j;
double _d, __d;

b+=n-1;
c+=n-1;
_i=n;
while(_i--)
  {
  __d=*b;
  _d=*tdU[_i]*__d;
  for (_j=_i+1;_j<n;_j++)
    {
    c++;
    b++;
    _d+=tdU[_i][_j-_i]**b;
    *c+=tdU[_i][_j-_i]*__d;
    }
  c-=_j-_i;
  b-=_j-_i;
  c[1]=_d;
  }
}


//
//This function multiple matrix and vector: c = Ab
inline void multiple_origin_matrix_origin_vector(unsigned int ni, unsigned int nj,double **A,double *b,double *c)
{
register unsigned int _j;
register double *a;
c+=ni, A+=ni;
while (ni--)
  {
  A--; 
  c--, *c=0.00;
  _j=nj;
  b+=nj;
  a=*A+nj;
  while (_j--)
    {
    b--;
    a--;
    *c+=*a**b;
    }
  }
}
//                                              T
//This function get vector matrix product: c = A b
inline void multiple_transp_matrix_origin_vector(unsigned int ni,unsigned int nj,double **A,double *b,double *c)
{
register unsigned int _j;
register double *a;
_j=nj, c+=nj;
while (_j--) 
  { 
  c--, *c=0.00;
  }
b+=ni, A+=ni;
while (ni--) 
  {
  b--;
  A--;
  _j=nj; 
  c+=nj;
  a=*A+nj; 
  while (_j--) 
    {
    c--;
    a--;
    *c+=*a**b; 
    }
  }
}
//                                          T   T
//This function get vector matrix product: c = b A
inline void multiple_transp_vector_origin_matrix(unsigned int ni,unsigned int nj,double **A,double *b,double *c)
{
multiple_transp_matrix_origin_vector(ni,nj,A,b,c);
}
//                                          T   T T
//This function get vector matrix product: c = b A
inline void multiple_transp_vector_transp_matrix(unsigned int ni,unsigned int nj,double **A,double *b,double *c)
{
multiple_origin_matrix_origin_vector(ni,nj,A,b,c);
}
//                                                   T
//This function  get vector matrix product: C = a * b . Result expands into matrix
inline void multiple_origin_vector_transp_vector(unsigned int ni,unsigned int nj,double *a,double *b,double **C)
{
register unsigned int _i, _j;
register double *c;
b+=ni;
C+=ni;
_i=ni; 
while(_i--)
  {
  C--;
  b--;
  a+=nj;
  c=*C+nj, _j=nj;
  while (_j--)
    {
    a--;
    c--;
    *c=*a**b; 
    }
  }
}
//This function perform multiple matrix C = AB
inline void multiple_origin_matrix_origin_matrix(unsigned int ni,unsigned int nj,unsigned int nk,double **A,double **B,double **C)
{
register unsigned int _i, _j, _k;
for (_i=0;_i<ni;_i++)
  for (_j=0;_j<nj;_j++)
    for (C[_i][_j]=0., _k=0;_k<nk;_k++)
      C[_i][_j]+=A[_i][_k]*B[_k][_j];
}
//                                           T
//This function perform multiple matrix C = A B
inline void multiple_transp_matrix_origin_matrix(unsigned int ni,unsigned int nj,unsigned int nk,double **A,double **B,double **C)
{
register unsigned int _i, _j, _k;
for (_i=0;_i<ni;_i++)
  for (_j=0;_j<nj;_j++)
    for (C[_i][_j]=0., _k=0;_k<nk;_k++)
      C[_i][_j]+=A[_k][_i]*B[_k][_j];
}
//                                            T
//This function perform multiple matrix C = AB
inline void multiple_origin_matrix_transp_matrix(unsigned int ni,unsigned int nj,unsigned int nk,double **A,double **B,double **C)
{
register unsigned int _i, _j, _k;
for (_i=0;_i<ni;_i++)
  for (_j=0;_j<ni;_j++)
    for (C[_i][_j]=0., _k=0;_k<nk;_k++)
      C[_i][_j]+=A[_i][_k]*B[_j][_k];
}
//                                           T T
//This function perform multiple matrix C = A B
inline void multiple_transp_matrix_transp_matrix(unsigned int ni,unsigned int nj,unsigned int nk,double **A,double **B,double **C)
{
register unsigned int _i, _j, _k;
for (_i=0;_i<nj;_i++)
  for (_j=0;_j<ni;_j++)
    for (C[_i][_j]=0., _k=0;_k<nk;_k++)
      C[_i][_j]+=A[_k][_i]*B[_j][_k];
}

//This function multiple matrix scalar
//Note A and B can be the same
char multiple_matrix_scalar(t_dmatrix *A,t_dmatrix *B,register double _f)
{
register unsigned int  _i,_j;
if ((A->ni!=B->ni)||(A->nj!=B->nj))
  return FALSE;
for (_i=0;_i<B->ni;_i++)
  for (_j=0;_j<B->nj;_j++)
    A->d[_i][_j]=B->d[_i][_j]*_f;
return TRUE;
}


//This function compute matrix row norm
double get_row_norm(unsigned int _i,unsigned int m,double **dA)
{
register double norm=0.00;
register double *row;
row=dA[_i];
while (m--) { norm+=sqrd(*row), row++; }
return norm;
}
//This function compute matrix column norm
double get_column_norm(unsigned int n,unsigned int _j,register double **dA)
{
register double norm=0.00;
while (n--) { norm+=sqrd(dA[n][_j]); }
return norm;
}

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
char LU_decomposition (unsigned int n, double **matrix, unsigned int *p_order, char *twoness, double *temp)
{
register unsigned int i,i_max,j,k;
double big,sum,f_temp;
// set 1 for initially even (0) transmutations number
*twoness=1;
// search for the largest element in each row; save the scaling in the f_temporary array temp and return zero if the matrix is singular
for (i=0;i<n;i++)
  {
  big=FALSE;
  for (j=0;j<n;j++)
    if ((f_temp=fabs(matrix[i][j]))>big)
      big=f_temp;
  if (!big)
    return FALSE;
  temp[i]=big;
  }
// the main loop for the Crout's algorithm
for (j=0;j<n;j++)
  {
  // this is the part matrix) of the algorithm except for i==j
  for (i=0;i<j;i++)
    {
    sum=matrix[i][j];
    for (k=0;k<i;k++)
      sum-=matrix[i][k]*matrix[k][j];

    matrix[i][j]=sum;
    }
  // initialize for the search for the largest pivot element
  big=FALSE;
  i_max=j;
  // this is the part matrix) for i==j and part b) for i>j plus pivot search
  for (i=j;i<n;i++)
    {
    sum=matrix[i][j];
    for (k=0;k<j;k++)
      sum-=matrix[i][k]*matrix[k][j];
    matrix[i][j]=sum;
    // is the figure of merit for the pivot better than the best so far?
    if ((f_temp=temp[i]*fabs(sum))>=big)
      {
      big=f_temp;
      i_max=i;
      }
    }
  // interchange rows, if needed, change parity and the scale factor
  if (i_max!=j)
    {
    for (k=0;k<n;k++)
      {
      f_temp=matrix[i_max][k];
      matrix[i_max][k]=matrix[j][k];
      matrix[j][k]=f_temp;
      }
    *twoness=-(*twoness);
    temp[i_max]=temp[j];
    }
  // store the index
  p_order[j]=i_max;
  // if the pivot element is zero, the matrix is singular but for some applications matrix tiny number is desirable instead
  if (!matrix[j][j])
    matrix[j][j]=TINY;
  // finally, divide by the pivot element
  if (j<n-1)
    {
    f_temp=TRUE/matrix[j][j];
    for (i=j+1;i<n;i++)
      matrix[i][j]*=f_temp;
    }
  }
return TRUE;
}

// the back substitution
void LU_back_substitution (unsigned int n,double **A,unsigned int *p_order,double *y)
{
register int _i,_j,ip,ii=-1;
register double sum;
// First step of backsubstitution; the only wrinkle is to unscramble the permutation order. Note: the algorithm is optimized for matrix possibility of large amount of zeroes in result_vector
for (_i=0;_i<n;_i++)
  {
  ip=p_order[_i];
  sum=y[ip];
  y[ip]=y[_i];
  if (ii>=0)
    for (_j=ii;_j<_i;_j++)
      sum-=A[_i][_j]*y[_j];
  else
    if (sum)
      ii=_i; // matrix nonzero element encounted
  y[_i]=sum;
  }
// the second step
for (_i=n-1;_i>=0;_i--)
  {
  sum=y[_i];
  for (_j=_i+1;_j<n;_j++)
    sum-=A[_i][_j]*y[_j];
  y[_i]=sum/A[_i][_i];
  }
}

// the matrix invertation
void LU_invert_dmatrix (unsigned int n,double **A,unsigned int *p_order,double **IA,double *_t)
{
register unsigned int _i,_j;
_j=n;
while (_j--)
  {
  _i=n, _t+=n; while (_i--) { _t--, *_t=0.; }
  _t[_j]=1.;
  LU_back_substitution(n,A,p_order,_t);
  _i=n, _t+=n; while (_i--) { _t--, IA[_i][_j]=*_t; }
  }
}

// calculating determinant
double LU_determinant (unsigned int n,double **A,unsigned int *p_order,char twoness)
{
register unsigned int _i;
double det=(double)(twoness);
_i=n; while (_i--) det*=A[_i][_i];
return det;
}

/*************************************** G A U S S   E l i m i n a t i o n ******************************************/

//This function solves linear system of equation with gauss-elimination method. It is preffered when matrix will be used once only.
unsigned int gauss_solve_dmatrix(unsigned int n, double **A,double *b,double badly_tol)
{
register unsigned int _i, _j, _l;
register double _d, *_dp, *__dp, max_d=0.;
//Stage 1. Build L matrix with eliminations
_i=n;
while (--_i)
  {
  //Stage a. choose the biggest item in the column and reorder rows if necessary
  _l=_i, _d=fabs(A[_i][_i]), _j=_i; while (_j--) if (_d<fabs(A[_j][_i])) { _l=_j, _d=fabs(A[_j][_i]); }
  //Stage b. Check validity of maximal value
  if (!_d) return YERROR_DATA_CONSISTMENT;
  else if (max_d<_d) max_d=_d;
  //Stage c. Do exchanging if necessary
  if (_l!=_i)
    {
    _d=b[_i], b[_i]=b[_l], b[_l]=_d;
    _dp=A[_i], A[_i]=A[_l], A[_l]=_dp;
    }
  //Stage d. do row-wise elimination
  _j=_i; 
  while (_j--) 
    { _d=A[_j][_i]/A[_i][_i], A[_j][_i]=0., _dp=A[_i], __dp=A[_j], _l=_i; while (_l--) { (*__dp)-=_d*(*_dp), _dp++, __dp++; } b[_j]-=_d*b[_i]; }
  }
//Stage 2. Do substitution
for (_i=0;_i<n;_i++)
  { 
  _d=0., _dp=A[_i], __dp=b, _j=_i;
  while (_j--)
    { _d+=(*_dp)*(*__dp), _dp++, __dp++; }
  if (fabs(*_dp)<badly_tol*max_d)
    {
    if (!*_dp) return YERROR_DATA_CONSISTMENT;
    else       return YERROR_LEGAL;  //Badly conditioned (can be solvable with smaller badly_tol)
    }
  *__dp=(*__dp-_d)/(*_dp);
  }
//Exits
return YERROR_OK;
}

/******************************************* S V D ******************************************************************/


/*
 *------------------------------------------------------------------------
 *                               Bidiagonalization
 */
 /// Left Hausholder Transformations
double sdiag__left_Hausholder_transformation(unsigned int ni,unsigned int nj,double **A,double **U,unsigned int i)
{
unsigned int j,k;
double scale, b, _d, new_Aii;  

scale=0.; for (k=i+2;k<ni;k++) scale+=fabs(A[k][i]);
if (!scale) return A[i+1][i]; 
else scale+=fabs(A[i+1][i]);  
b=0.; for (k=i+1;k<ni;k++) b+=sqrd(A[k][i]/=scale);  
new_Aii=sqrt(b);  
b=-1./(b-A[i+1][i]*new_Aii);  
A[i+1][i]-=new_Aii;
for (j=i+1;j<nj;j++)
  {
  _d=0.; for (k=i+1;k<ni;k++) _d+=A[k][i]*A[k][j]; 
  _d*=b; for (k=i+1;k<ni;k++) A[k][j]+=A[k][i]*_d;
  }
for (j=0;j<ni;j++) 
  {
  _d=0.; for (k=i+1;k<ni;k++) _d+=A[k][i]*U[j][k]; 
  _d*=b; for (k=i+1;k<ni;k++) U[j][k]+=A[k][i]*_d;
  }
return new_Aii*scale; 
}
double diag__left_Hausholder_transformation(unsigned int ni,unsigned int nj,double **A,double **U,unsigned int i)
{
unsigned int j,k;
double scale, b, _d, new_Aii;  

scale=0.; for (k=i+1;k<ni;k++) scale+=fabs(A[k][i]);
if (!scale) return A[i][i];
else scale+=fabs(A[i][i]);  
b=0.; for (k=i;k<ni;k++) b+=sqrd(A[k][i]/=scale);  
new_Aii=sqrt(b);  
b=-1./(b-A[i][i]*new_Aii);  
A[i][i]-=new_Aii;
for (j=i+1;j<nj;j++) 
  {
  _d=0.; for (k=i;k<ni;k++) _d+=A[k][i]*A[k][j]; 
  _d*=b; for (k=i;k<ni;k++) A[k][j]+=A[k][i]*_d;
  }
for (j=0;j<ni;j++) 
  {
  _d=0.; for (k=i;k<ni;k++) _d+=A[k][i]*U[j][k]; 
  _d*=b; for (k=i;k<ni;k++) U[j][k]+=A[k][i]*_d;
  }
return new_Aii*scale; 
}

/// Right Hausholder Transformations
double sdiag_right_Hausholder_transformation (unsigned int ni,unsigned int nj,double **A,double **V,unsigned int i)
{
unsigned int j,k;
double scale, b, _d, new_Aii;

scale=0.; for (k=i+2;k<nj;k++) scale+=fabs(A[i][k]);
if (!scale) return A[i][i+1];
else scale+=fabs(A[i][i+1]); 
b=0.; for (k=i+1;k<nj;k++) b+=sqrd(A[i][k]/=scale);
new_Aii=sqrt(b); 
b=-1./(b-A[i][i+1]*new_Aii);  
A[i][i+1]-=new_Aii;  
for (j=i+1;j<ni;j++) 
  {
  _d=0.; for (k=i+1;k<nj;k++) _d+=A[i][k]*A[j][k];
  _d*=b; for (k=i+1;k<nj;k++) A[j][k]+=A[i][k]*_d;
  }
for (j=0;j<nj;j++) 
  {
  _d=0.; for (k=i+1;k<nj;k++) _d+=A[i][k]*V[j][k];
  _d*=b; for (k=i+1;k<nj;k++) V[j][k]+=A[i][k]*_d;
  }
return new_Aii*scale;
}
double diag_right_Hausholder_transformation (unsigned int ni,unsigned int nj,double **A,double **V,unsigned int i)
{
unsigned int j,k;
double scale, b, _d, new_Aii;  

scale=0.; for (k=i+1;k<nj;k++) scale+=fabs(A[i][k]);
if (!scale) return A[i][i];
else scale+=fabs(A[i][i]); 
b=0.; for (k=i;k<nj;k++) b+=sqrd(A[i][k]/=scale); 
new_Aii=sqrt(b); 
b=-1./(b-A[i][i]*new_Aii);  
A[i][i]-=new_Aii;  
for (j=i+1;j<ni;j++) 
  {
  _d=0.; for (k=i;k<nj;k++) _d+=A[i][k]*A[j][k]; 
  _d*=b; for (k=i;k<nj;k++) A[j][k]+=A[i][k]*_d;
  }
for (j=0;j<nj;j++) 
  {
  _d=0.; for (k=i;k<nj;k++) _d+=A[i][k]*V[j][k]; 
  _d*=b; for (k=i;k<nj;k++) V[j][k]+=A[i][k]*_d;
  }
return new_Aii*scale; 
}
/// Bidiagonalization. It produces either super_diag or sub_diag
double Hausholder_bidiagonalization (double *diag,double *sdiag,unsigned int ni,unsigned int nj,unsigned int nr,double **A,double **U,double **V)
{
double norm;
unsigned int i;
*sdiag=norm=0.0; 
if (ni<=nj)  
  {//generate diag/subdiag pair
  for (i=0;i<ni-2;i++)
    {
     diag[i  ]= diag_right_Hausholder_transformation(ni,nj,A,V,i);
    sdiag[i+1]=sdiag__left_Hausholder_transformation(ni,nj,A,U,i);
    if (norm<(fabs(diag[i])+fabs(sdiag[i]))) norm=fabs(diag[i])+fabs(sdiag[i]);
    }
   diag[i  ]=diag_right_Hausholder_transformation(ni,nj,A,V,i);
  sdiag[i+1]=A[i+1][i],  diag[i+1]=diag_right_Hausholder_transformation(ni,nj,A,V,i+1);
  if (norm<(fabs(diag[i])+fabs(sdiag[i])+fabs(diag[i+1]))) norm=fabs(diag[i])+fabs(sdiag[i])+fabs(diag[i+1]);
  }
else
  {//generate diag/superdiag pair
  for (i=0;i<nj-2;i++)
    {
     diag[i  ]= diag__left_Hausholder_transformation(ni,nj,A,U,i);
    sdiag[i+1]=sdiag_right_Hausholder_transformation(ni,nj,A,V,i);
    if (norm<(fabs(diag[i])+fabs(sdiag[i]))) norm=fabs(diag[i])+fabs(sdiag[i]);
    }
   diag[i  ]=diag__left_Hausholder_transformation(ni,nj,A,U,i);
  sdiag[i+1]=A[i][i+1],  diag[i+1]=diag__left_Hausholder_transformation(ni,nj,A,U,i+1);
  if (norm<(fabs(diag[i])+fabs(sdiag[i])+fabs(diag[i+1]))) norm=fabs(diag[i])+fabs(sdiag[i])+fabs(diag[i+1]);
  }

return norm;
}

/**************************  QR-diagonalization of a bidiagonal matrix ************************************************/
 /*
  *
  * After bidiagonalization we get a bidiagonal matrix J:
  *      (1)  J = U' * A * V
  * The present method turns J into a matrix JJ by applying a set of
  * orthogonal transforms
  *      (2)  JJ = S' * J * T
  * Orthogonal matrices S and T are chosen so that JJ were also a
  * bidiagonal matrix, but with superdiag elements smaller than those of J.
  * We repeat (2) until non-diag elements of JJ become smaller than EPS
  * and can be disregarded.
  * Matrices S and T are constructed as
  *      (3)  S = S1 * S2 * S3 ... Sn, and similarly T
  * where Sk and Tk are matrices of simple rotations
  *      (4)  Sk[i,j] = i==j ? 1 : 0 for all i>k or i<k-1
  *               Sk[k-1,k-1] = cos(Phk),  Sk[k-1,k] = -sin(Phk),
  *               SK[k,k-1] = sin(Phk),  Sk[k,k] = cos(Phk), k=2..N
  * Matrix Tk is constructed similarly, only with angle Thk rather than Phk.
  *
  * Thus left multiplication of J by SK' can be spelled out as
  *      (5)  (Sk' * J)[i,j] = J[i,j] when i>k or i<k-1,
  *                                [k-1,j] = cos(Phk)*J[k-1,j] + sin(Phk)*J[k,j]
  *                                [k,j] =  -sin(Phk)*J[k-1,j] + cos(Phk)*J[k,j]

  * That is, k-1 and k rows of J are replaced by their linear combinations;
  * the rest of J is unaffected. Right multiplication of J by Tk similarly
  * changes only k-1 and k columns of J.
  * Matrix T2 is chosen the way that T2'J'JT2 were a QR-transform with a
  * shift. Note that multiplying J by T2 gives rise to a J[2,1] element of
  * the product J (which is below the main diagonal). Angle Ph2 is then
  * chosen so that multiplication by S2' (which combines 1 and 2 rows of J)
  * gets rid of that elemnent. But this will create a [1,3] non-zero element.
  * T3 is made to make it disappear, but this leads to [3,2], etc.
  * In the end, Sn removes a [N,N-1] element of J and matrix S'JT becomes
  * bidiagonal again. However, because of a special choice
  * of T2 (QR-algorithm), its non-diag elements are smaller than those of J.
  *
  * All this process in more detail is described in
  *      J.H. Wilkinson, C. Reinsch. Linear algebra - Springer-Verlag, 1971
  *
  * If during transforms (1), JJ[N-1,N] turns 0, then JJ[N,N] is a singular
  * number (possibly with a wrong (that is, negative) sign). This is a
  * consequence of Frantsis' Theorem, see the book above. In that case, we can
  * eliminate the N-th row and column of JJ and carry out further transforms
  * with a smaller matrix. If any other superdiag element of JJ turns 0,
  * then JJ effectively falls into two independent matrices. We will process
  * them independently (the bottom one first).
  *
  * Since matrix J is a bidiagonal, it can be stored efficiently. As a matter
  * of fact, its N diagonal elements are in array Sig, and superdiag elements
  * are stored in array super_diag.
  */
 // Carry out U * S with a rotation matrix
 // S (which combines i-th and j-th columns
 // of U, i>j)
inline void rotate_matrix(unsigned int ni,double **A,unsigned int i,unsigned int j,double cos_phy,double sin_phy)
{
register unsigned int k;
register double akj;
k=ni;
while(k--)
  {
  akj=A[k][j];
  A[k][j]= cos_phy*akj+sin_phy*A[k][i];
  A[k][i]=-sin_phy*akj+cos_phy*A[k][i];
  }
}
/*
 * A diagonal element J[l-1,l-1] turns out 0 at the k-th step of the
 * algorithm. That means that one of the original matrix' singular numbers

 * shall be zero. In that case, we multiply J by specially selected
 * matrices S' of simple rotations to eliminate a superdiag element J[l-1,l].
 * After that, matrix J falls into two pieces, which can be dealt with
 * in a regular way (the bottom piece first).
 *
 * These special S transformations are accumulated into matrix U: since J
 * is left-multiplied by S', U would be right-multiplied by S. Transform
 * formulas for doing these rotations are similar to (5) above. See the
 * book cited above for more details.
 */
void RipThrough (unsigned int ni,double *super_diag,double *diag,double **U,unsigned int i,unsigned int j,double epsilon)
{
double cos_phy=0.0,sin_phy=1.0,f; // Accumulate cos,sin of Ph
unsigned int k;

// The first step of the loop below
// (when i==l) would eliminate J[l-1,l],
// which is stored in super_diag(l)
// However, it gives rise to J[l-1,l+1]
// and J[l,l+2]
// The following steps eliminate these
// until they fall below
// significance
for (k=i;k<=j;k++)
  {
  f=sin_phy*super_diag[k];
  super_diag[k]*=cos_phy;
  if ( fabs(f) <= epsilon )
    break ; // Current J[l-1,l] became unsignificant
  cos_phy=diag[k];
  sin_phy=-f; // unnormalized sin/cos
  diag[k]=hypot(cos_phy, sin_phy); // sqrt(sin^2+cos^2)
  cos_phy/=diag[k];
  sin_phy/=diag[k]; // Normalize sin/cos
  rotate_matrix(ni,U,k,i-1,cos_phy,sin_phy);
  }
}
// We're about to proceed doing QR-transforms
// on a (bidiag) matrix J[1:k,1:k]. It may happen
// though that the matrix splits (or can be
// split) into two independent pieces. This function
// checks for splitting and returns the lowerbound
// index l of the bottom piece, J[l:k,l:k]
unsigned int getSubmatrixToWorkOn (unsigned int ni,double *super_diag,double *diag,double **U,unsigned int i,double epsilon)
{
unsigned int j;	
for (j=i;j;j--) 
  if (fabs(super_diag[j])<=epsilon)
    return j; // The breaking point: zero J[l-1,l]
  else if (fabs(diag[j-1])<=epsilon) // Diagonal J[l,l] turns out 0
    {
    // meaning J[l-1,l] _can_ be made
    RipThrough (ni,super_diag,diag,U,i,j,epsilon); // zero after some rotations
//    printf("Ripping\n");
    return j;
    }
return FALSE; // Deal with J[1:k,1:k] as a whole
}
//Terrible function!
/// Diagonalize root module
void diagonalize(unsigned int ni,unsigned int nj,unsigned int nr,double *super_diag,double *diag,double **U,double **V,double epsilon)
{
int k, l, i, j;
double shift,Jk2k1,Jk1k,Jk1k1,Jkk,Jll,cos_th,sin_th,Ji1i1,Ji1i,Jii,norm_f,Jii1,cos_ph,sin_ph;
unsigned int sweep=0;
k=nr;
while (k--) // QR-iterate upon J[l:k,l:k]
  {
  while (l=getSubmatrixToWorkOn(ni,super_diag,diag,U,k,epsilon),fabs(super_diag[k])>epsilon) // until superdiag J[k-1,k] becomes 0
    {
    sweep++;
    // Compute a QR-shift from a bottom
    // corner minor of J[l:k,l:k] order 2
    Jk2k1=super_diag[k-1], // J[k-2,k-1]
    Jk1k=super_diag[k],
    Jk1k1=diag[k-1], // J[k-1,k-1]
    Jkk=diag[k], 
    Jll=diag[l]; // J[l,l]

    shift=(Jk1k1-Jkk)*(Jk1k1+Jkk)+(Jk2k1-Jk1k)*(Jk2k1+Jk1k);
    shift/=2*Jk1k*Jk1k1;
    shift+=(shift<0 ? -1 : 1)*sqrt(shift*shift+1);
    shift=((Jll-Jkk)*(Jll+Jkk)+Jk1k*(Jk1k1/shift-Jk1k))/Jll ;
    //end of corner minor           
    // Carry on multiplications by T2, S2, T3...
    cos_th=1;
    sin_th=1;
    Ji1i1=diag[l];  // J[i-1,i-1] at i=l+1...k
    for (i=l+1;i<=k;i++) 
      {
      Ji1i=super_diag[i],Jii=diag[i]; // J[i-1,i] and J[i,i]
      sin_th*=Ji1i;
      Ji1i*=cos_th;
      cos_th=shift;
      norm_f=super_diag[i-1]=hypot(cos_th, sin_th);
      cos_th/=norm_f;
      sin_th/=norm_f;

      // Rotate J[i-1:i,i-1:i] by Ti
      shift=cos_th*Ji1i1+sin_th*Ji1i ; // new J[i-1,i-1]
      Ji1i=-sin_th*Ji1i1+cos_th*Ji1i ; // J[i-1,i] after rotation
      Jii1=Jii*sin_th; // Emerged J[i,i-1]
      Jii*=cos_th; // new J[i,i]
      rotate_matrix(nj,V,i,i-1,cos_th,sin_th); // Accumulate T rotations in V
      //show_matrix(nj,nj,V);
      cos_ph=shift;
      sin_ph=Jii1; // Make Si to get rid of J[i,i-1]
      norm_f=hypot(cos_ph,sin_ph); // New J[i-1,i-1]
      diag[i-1]=norm_f;
      if (norm_f==0)
        {
        // If norm =0, rotation angle
        cos_ph=cos_th;
        sin_ph=sin_th; // can be anything now
        }
      else
        {
        cos_ph/=norm_f;
        sin_ph/=norm_f;
        }
	  
      // Rotate J[i-1:i,i-1:i] by Si
      shift=cos_ph* Ji1i+sin_ph*Jii; // New J[i-1,i]
      Ji1i1=-sin_ph*Ji1i+cos_ph*Jii; // New Jii, would carry over

      // as J[i-1,i-1] for next i
      rotate_matrix(ni,U,i,i-1,cos_ph,sin_ph); // Accumulate S rotations in U
      //show_matrix(ni,ni,U);
      // Jii1 disappears, sin_th would
      cos_th=cos_ph; // carry over a (scaled) J[i-1,i+1]
      sin_th=sin_ph;
      // to eliminate on the next i, cos_th
      // would carry over a scaled J[i,i+1]
      }
    super_diag[l]=0; // Supposed to be eliminated by now
    super_diag[k]=shift;
    diag[k]=Ji1i1;
    } // --- end-of-QR-iterations
  if (diag[k]<0)  { diag[k]=-diag[k], j=nj; while (j--) V[j][k]=-V[j][k]; } // Correct the sign of the sing number
  }
//printf("Done in %1d sweeps with tol==%f\n",sweep,epsilon);
}


//This function compute the SVD. super_diag is the help massive for Hausholder Transformation
//A[m:n], U[m:m], V[[n:n], diag[m], super_diag[m]; m<=n
void singular_value_decomposition(unsigned int ni,unsigned int nj,unsigned int nr,double **A,double **B,double **U,double *diag,double *sdiag,double **V,double epsilon)
{
register unsigned int i, j, k;
register double _d;
i=ni; while (i--) { j=ni; while (j--) U[i][j]=0.; U[i][i]=1.; } 
i=nj; while (i--) { j=nj; while (j--) V[i][j]=0.; V[i][i]=1.; } 
//printf("A:\n");  show_matrix(ni,nj,A);
epsilon*=Hausholder_bidiagonalization(diag,sdiag,ni,nj,nr,A,U,V); // Significance threshold !!!
//for (i=0;i<ni;i++) { for (j=nj;j>0;j--) U[i][j]=U[i][j-1]*sdiag[j]+U[i][j]*diag[j]; U[i][0]*=diag[0]; }
//for (i=0;i<ni;i++) for (j=0;j<nj;j++) for (A[i][j]=0, k=0;k<nj;k++) A[i][j]+=U[i][k]*V[j][k];
if (ni<nj) 
  {//Transposition take some time but do not reduce an accuracy
  i=ni; while (--i) { j=i; while (j--) { _d=U[i][j], U[i][j]=U[j][i], U[j][i]=_d; } }
  diagonalize(nj,ni,ni,sdiag,diag,V,U,epsilon);
  }
else diagonalize(ni,nj,nj,sdiag,diag,U,V,epsilon);
if ( (B)) 
  {//Accuracy check is asked
         i=ni; while (i--) { j=nj; while (j--) { k=nr, A[i][j]=0.; while (k--) A[i][j]+=U[i][k]*diag[k]*V[j][k]; } }
  _d=0., i=ni; while (i--) { j=nj; while (j--) _d+=sqrd(A[i][j]-B[i][j]); } printf("\nInfo. Norm of E matrix after svd is %16.8f\n",_d);
  }
//printf("A`=U.S.V:\n");  show_matrix(ni,nj,A);
//printf("U:\n");         show_matrix(ni,ni,U);
//printf("V:\n");         show_matrix(nj,nj,V);
//printf("diga and sdiag:\n"); for (i=0;i<nr;i++) printf("%1d: d=%f\tsd=%f\n",i,diag[i],sdiag[i]);
}

//This function calculates svd using previous one and some other twics:
//symm - optimizing perfomance by symmetric algorithm
//reorder - sort eigenvalues after svd
void svd(char symm,char reorder,unsigned int ni,unsigned int nj,unsigned int nr,double **A,double **B,double **U,double *s,double **V,double *temp)
{
int *id=(int*)temp;
double _d;
unsigned int i, j, k;

//General assymetrical case [not implemented]
//Calculate SVD of A
singular_value_decomposition(ni,nj,nr,A,B,U,s,temp,V,TINY);
//Order eignenvalues if required
if (reorder)
  {
  i=nr; while (i--) id[i]=i;
  di_qsort(nr,s,id);
  i=nr;
  while (--i)
    if (id[i]!=i)
      {//Place i-th value onto its position
      j=i; while (--j) if (id[j]==i) break; //Find appropriate position
      id[j]=id[i], j=id[i], id[i]=i;
      k=ni; while (k--) { _d=U[k][i], U[k][i]=U[k][j], U[k][j]=_d;}
      k=nj; while (k--) { _d=V[k][i], V[k][i]=V[k][j], V[k][j]=_d; } //Swap vectors
      }
  //Check reordering
  if ( (B)) 
    {//Accuracy check is asked
           i=ni; while (i--) { j=nj; while (j--) { k=nr, A[i][j]=0.; while (k--) A[i][j]+=U[i][k]*s[k]*V[j][k]; } }
    _d=0., i=ni; while (i--) { j=nj; while (j--) _d+=sqrd(A[i][j]-B[i][j]); } printf("\nInfo. Norm of E matrix after reordering of svd is %16.8f\n",_d);
    }
  }
}


#define ROTATE(a_ij,a_kl,_e,_f,_s,tau)  _f=a_ij,_e=a_kl,a_ij=_f-_s*(_e+_f*tau),a_kl=_e-_s*(_f+_e*tau)
//This is deprecated algorithm of Jacobi rotations of symmetriñ matrix. Please use svd instead of it anywhere you can!
//Note. This method performs the diagonalization only. No any eigen vectors are calculated. Diagonalized A is S. A can be restored from untouched lower triangle.
char diagonalization_Jacobi(unsigned int n_iter,t_dmatrix *A)
{
int _i,_j,step;
register int _k;
double tau,_s,_c;
register double _e,_f;

if (A->ni!=A->nj)
  return FALSE;
//Perform rotation cycles
for (step=0;step<n_iter;step++)
  {
  //Summ all offdiagonal elemnts
  for (_f=0.0,_i=0;_i<A->ni;_i++)
    for (_j=_i+1;_j<A->nj;_j++)
      _f+=fabs(A->d[_i][_j]);
  if (!_f) return TRUE; //exit on converge show_matrix(S) show_matrix(A)
  //Choose step depended strategy of elimination
  if (step<4)
    tau=0.2*_f/(A->ni*A->ni);
  else
    tau=0.0; //No limits
  for (_i=0;_i<A->ni-1;_i++)
    for (_j=_i+1;_j<A->nj;_j++)
      {
      _f=1.e+2*fabs(A->d[_i][_j]);
      if ((step>4)&&(_f<TINY))
        {
        A->d[_i][_j]=0.0;
        A->d[_j][_i]=0.0;
        }
      else
        if (fabs(A->d[_i][_j])>tau)
          {
          _e=A->d[_j][_j]-A->d[_i][_i];
          if (((double)fabs(_e)+_f)==(double)fabs(_e))
            _f=A->d[_i][_j]/_e;
          else
            {
            _e*=0.5/A->d[_i][_j];
            _f=1.0/(fabs(_e)+sqrtf(1.0+_e*_e));
            if (_e<0.0)
              _f=-_f;
            _c=1.00/sqrtf(1.00+_f*_f);
            _s=_f*_c;
            tau=_s/(1.00+_c);
            _e=_f*A->d[_i][_j];
            A->d[_i][_i]-=_e;
            A->d[_j][_j]+=_e;
            A->d[_i][_j]=0.0;
            for (_k=0   ;_k<=_i-1;_k++) //Case  1 <= _k <  _i
              ROTATE(A->d[_k][_i],A->d[_k][_j],_e,_f,_s,tau);
            for (_k=_i+1;_k<=_j-1;_k++) //Case _i <  _k <=  _j
              ROTATE(A->d[_i][_k],A->d[_k][_j],_e,_f,_s,tau);
            for (_k=_j+1;_k< A->ni;_k++) //Case _j <  _k < A->ni
              ROTATE(A->d[_i][_k],A->d[_j][_k],_e,_f,_s,tau);
            }
          }
      }
  }
return FALSE; //still not converged, sorry
}
#undef ROTATE


//This function performs Grahamm-Smidt ortonormalization
inline void gram_schmidt_ortonormalization(unsigned int ni,unsigned int nj,double **AT,double *_t)
{
register unsigned int _i, _j, _l;
register double *dp, *dp2, _d;
*_t=calc_vect_norm(nj,*AT);
for (_i=1; _i<ni; _i++)
  for ( _j=0; _j<_i; _j++)
    {//calc proj [u] (v) = v-(<v.u>/<u.u>)*u
    _d=calc_vect_vect_scalar_product(nj,AT[_i],AT[_j])/_t[_j];
    _l=nj, dp=AT[_i], dp2=AT[_j]; while (_l--) { *dp-=*dp2*_d, dp++, dp2++; }
    _t[_i]=calc_vect_norm(nj,AT[_i]);
    }
_i=ni; while (_i--) { _d=1./sqrt(_t[_i]), dp=AT[_i], _l=nj; while (_l--) { *dp*=_d, dp++; } }
}

//------------------- T H E   T R I A N G L E   M A T R I C E S   P A R T ---------------------------------

//This function generates triangle matrix form square self-multiplicaion

//This function calculates C=A.D.AT
inline void multiple_nXDXT(unsigned int ni,unsigned int nj,double **C,register double **A,double *d)
{
register double _d;
register unsigned int _k,_j,_i;
for (_i=0;_i<ni;_i++)
  {
  _j=_i+0x1;
  while (_j--)
    {
    _k=nj;
    _d=0.00;
    while (_k--) _d+=d[_k]*A[_i][_k]*A[_j][_k];
    C[_j][_i]=C[_i][_j]=_d;
    }
  }
}
inline void multiple_tXDXT(unsigned int ni,unsigned int nj,double **U,register double **A,double *d)
{
register double _d;
register unsigned int _k,_j,_i;
for (_i=0;_i<ni;_i++)
  {
  _j=_i+0x1;
  while (_j--)
    {
    _k=nj;
    _d=0.00;
    while (_k--) _d+=d[_k]*A[_i][_k]*A[_j][_k];
    U[_j][_i-_j]=_d;
    }
  }
}

//This function calculates C=AT.D.A
inline void multiple_nXTDX(unsigned int ni,unsigned int nj,double **C,register double **A,double *d)
{
register double _d;
register unsigned int _k,_j,_i;
for (_i=0;_i<nj;_i++)
  {
  _j=_i+0x1;
  while (_j--)
    {
    _k=ni;
    _d=0.00;
    while (_k--) _d+=d[_k]*A[_k][_i]*A[_k][_j];
    C[_j][_i]=C[_i][_j]=_d;
    }
  }
}
inline void multiple_tXTDX(unsigned int ni,unsigned int nj,double **U,register double **A,double *d)
{
register double _d;
register unsigned int _k,_j,_i;
for (_i=0;_i<nj;_i++)
  {
  _j=_i+0x1;
  while (_j--)
    {
    _k=ni;
    _d=0.00;
    while (_k--) _d+=d[_k]*A[_k][_i]*A[_k][_j];
    U[_j][_i-_j]=_d;
    }
  }
}

//This function calculates C=A.AT
inline void multiple_nXXT(unsigned int ni,unsigned int nj,double **C,register double **A)
{
register double _d;
register unsigned int _k,_j,_i;
for (_i=0;_i<ni;_i++)
  {
  _j=_i+0x1;
  while (_j--)
    {
    _k=nj;
    _d=0.00;
    while (_k--) _d+=A[_i][_k]*A[_j][_k];
    C[_j][_i]=C[_i][_j]=_d;
    }
  }
}
inline void multiple_tXXT(unsigned int ni,unsigned int nj,double **U,register double **A)
{
register double _d;
register unsigned int _k,_j,_i;
for (_i=0;_i<ni;_i++)
  {
  _j=+0x1;
  while (_j--)
    {
    _k=nj;
    _d=0.00;
    while (_k--) _d+=A[_i][_k]*A[_j][_k];
    U[_j][_i-_j]=_d;
    }
  }
}

//This function calculates C=AT.A
inline void multiple_nXTX(unsigned int ni,unsigned int nj,double **C,register double **A)
{
register double _d;
register unsigned int _k,_j,_i;
for (_i=0;_i<nj;_i++)
  {
  _j=_i+0x1;
  while (_j--)
    {
    _k=ni;
    _d=0.00;
    while (_k--) _d+=A[_k][_i]*A[_k][_j];
    C[_j][_i]=C[_i][_j]=_d;
    }
  }
}
inline void multiple_tXTX(unsigned int ni,unsigned int nj,double **U,register double **A)
{
register double _d;
register unsigned int _k,_j,_i;
for (_i=0;_i<nj;_i++)
  {
  _j=_i+0x1;
  while (_j--)
    {
    _k=ni;
    _d=0.00;
    while (_k--) _d+=A[_k][_i-_k]*A[_k][_j-_k];
    U[_j][_i-_j]=_d;
    }
  }
}

//This function perform modified Choletski factorization (Square root method)
//Note. The matrix is given here in form of lower triangle of simmetric A. Result define as L.D.LT=P.A.P so to minimize ||E|| in A-E=L.D.L, where A-E -s positively-defined well-conditioned matrix.
// mu - relaxation factor - recomended to set as 0.1
//Algorithm taken from Haw-ren Fang, Dianne P. O'Leary Modified cholesky algorithms: a catalog with new approaches 2006
//Note this function returns FALSE if exact value was calculated (i.e. ||E||=0. ) and TRUE otherwise 
char modified_cholesky_factorization_dmatrix(double mu,unsigned int n,unsigned int *p,double **A)
{
unsigned int _i,_j,_k;
double nu, beta2, sigma, _d;

for (_i=0;_i<n;_i++) p[_i]=_i;
//Stage 1 of relaxed factorization
for (_k=0;_k<n;_k++)
  {
  //Pivot on max diagonal
  for (_d=A[_k][_k], _j=_k, _i=_k+1;_i<n;_i++) if (A[_i][_i]>_d) { _j=_i, _d=A[_i][_i]; }
  if (_k!=_j) 
    {
    _i=p[_k], p[_k]=p[_j], p[_j]=_i;
    for (_i=0;   _i<_k;_i++) { _d=A[_k][_i], A[_k][_i]=A[_j][_i], A[_j][_i]=_d; }
    for (_i=_k+1;_i<_j;_i++) { _d=A[_i][_k], A[_i][_k]=A[_j][_i], A[_j][_i]=_d; }
    for (_i=_j+1;_i<n; _i++) { _d=A[_i][_k], A[_i][_k]=A[_i][_j], A[_i][_j]=_d; }
    _d=A[_k][_k], A[_k][_k]=A[_j][_j], A[_j][_j]=_d;
    }
  //Check the diagonal
  if ((_d=A[_k][_k])<EPSILON) goto STAGE2; 
  for ( _i=_k+1;_i<n;_i++) if ( (A[_i][_i]<-mu*_d)||(A[_i][_i]-A[_i][_k]*A[_i][_k]/_d<-mu*_d) ) goto STAGE2;
  //Do new row/column factorization
  for (_d=A[_k][_k],_i=_k+1;_i<n;_i++) 
    {
    A[_i][_k]/=_d; //L(k+1:n,k:k)=ck/ak
    for (_j=_k+1;_j<=_i;_j++) A[_i][_j]-=A[_i][_k]*A[_j][_k]*_d;    //Ak+1=Ak-ck*ck^T/ak  
    } 
  }
return FALSE;
//Stage2 - approximate calculations of positive-defined well-conditioned A`

STAGE2:
sigma=SMALL2;
//Rip throught
for (;_k<n;_k++)
  {
  //Search for biggest nu 
  if (_k!=n-1) 
    {
    nu=A[_k+1][_k]; for (_i=_k;_i<n;_i++) for (_j=_k;_j<_i;_j++) if (nu<fabs(A[_i][_j])) nu=fabs(A[_i][_j]);
    if (_k!=n-1) beta2=nu/sqrt((double)((n-_k)*(n-_k-1)));
    if (SMALL2>beta2) beta2=SMALL2;  //beta^2
    //Update sigma
    for (_d=0., _i=_k+1;_i<n;_i++) _d+=A[_i][_k]*A[_i][_k]; 
    _d= (_d/beta2 > SMALL2) ? _d/beta2 : SMALL2;
    }
  else _d=SMALL2;
  _d= ((A[_k][_k]+sigma)<_d) ? _d : A[_k][_k]+sigma;
  sigma=_d-A[_k][_k];
  //Do factorization
  for (A[_k][_k]=_d,_i=_k+1;_i<n;_i++) 
    {
    A[_i][_k]/=_d; //L(k+1:n,k:k)=ck/ak
    for (_j=_k+1;_j<=_i;_j++) A[_i][_j]-=A[_i][_k]*A[_j][_k]*_d;    //Ak+1=Ak-ck*ck^T/ak  
    } 
  }
return TRUE;
}

//This function perform Choletski transformation (Square root method)
//Note. The matrix is given here in form of upper triangle of simmetric A. 
//IF type=='U' the result define as UT.D.U ELSE IF type=='L' the result is L.D.LT ELSE YERROR_NIMPLEMENT
void cholesky_decomposition_tdmatrix(char type,unsigned int n,double **A)
{
register double _d,__d;
register unsigned int _k,_j,_i;
     if ( (type=='U')||(type=='u') )
       for (_i=0x0;_i<n;_i++)
         {
         //Get diagonal element
         _d=0.00, _k=_i; while (_k--) _d+=*A[_k]*A[_k][_i-_k]*A[_k][_i-_k];
         __d=(*A[_i]-=_d);
         if (!__d) __d=*A[_i]=TINY;
         //Get off-diagonal elements
         for (_j=_i+1;_j<n;_j++)
           {
           _d=0.00, _k=_i; while (_k--) _d+=*A[_k]*A[_k][_i-_k]*A[_k][_j-_k];
           A[_i][_j-_i]=(A[_i][_j-_i]-_d)/__d;
           }
         }
else if ( (type=='l')||(type=='L') )
       for (_i=1;_i<n;_i++)
         {
         __d=A[_i][0];
         __d*=(A[_i][0]/=A[0][0]); //Calculate new L
         //Calculate Aij-SUMM[k=l...j-1](Lik*Ljk*Dk)
         for (_j=1;_j<_i;_j++)
           {//Calculate summ
           _d=0., _k=_j; while (_k--) _d+=A[_i][_k]*A[_j][_k]*A[_k][_k];
           _d=A[_i][_j]-_d;              // L[_i][_j]*A[_j][_j] 
           A[_i][_j]=_d/A[_j][_j];       //update Aij     
           __d+=A[_i][_j]*_d;   
           }
         A[_i][_i]-=__d;  //Update Aii
         }
else ylib_errno=YERROR_NIMPLEMENTED;
}

//This function calculates L.b=c.
//Note. It suggests that L has unit diagonal.
void multiple_tdILmatrix_vector(unsigned int n,register  double *c,double **IL,register double *b)
{
register double *d;
register unsigned int _i,_j;
_i=n;
while (_i--)
  {
  d=IL[_i];
  c[_i]=b[_i];
  _j=_i;
  while (--_j)
    {
    *c+=*b**d;
    b++, c++, d++;
    }
  b-=_i;
  c-=_i;
  }
}

//This function calculates U.b=c.
//Note. U is symmetrical matrix stored in upper triangle matrix for memory saving purposes.
void multiple_tdmatrixU_vector(unsigned int n,register  double *c,double **tU,register double *b)
{
register unsigned int _i,_j;
register double _d;
_i=n;
while(_i--)
  {
  *c=0.00;
  c++;
  }
_i=n;
while(_i--)
  {
  c--;
  _j=n-_i; //It means j=n in case of uncompressed matrix
  b+=n;
  while (--_j)
    {
    b--;
    _d=tU[_i][_j]**b;
    *c+=_d;
    c[_j]+=_d;
    }
  *c+=*tU[_i]**b;
  }
}

//This function calculates U.b=c.
//Note. It suggests that U has unit diagonal.
void multiple_tdIUmatrix_vector(unsigned int n,register  double *c,double **IU,register double *b)
{
register double *d;
register unsigned int _i,_j;
_i=n;
while (_i--)
  {
  d=IU[_i];
  c[_i]=b[_i];
  _j=n-_i;
  while (--_j)
    {
    *c+=*b**d;
    b++, c++, d++;
    }
  b-=n-_i;
  c-=n-_i;
  }
}
void multiple_tdIUTmatrix_vector(unsigned int n,register double *c,double **IU,register double *b)
{
register unsigned int _i,_j;
*c=*b;
_i=n;
while (_i--) { *c=*b, c++, b++; }
_i=n-0x1;
while (_i--)
  {
  _j=n-_i-0x1;
  c-=_j;
  b-=_j;
  while (_j--)
    {
    *c+=IU[_i][_j]**b;
    c++, b++;
    }
  }
}


//This function invert triangle matrix U
inline void invert_tdIUmatrix(unsigned int n,double **IIU,double **IU)
{
register double _d;
register unsigned int _k,_j,_i;

IIU[0x0][0x1]=-IU[0x0][0x1];
for (_j=0x2;_j<n;_j++)
  {
  _i=_j;
  while (_i--)
    {
    _k=_j-_i;
    _d=-IU[_i][_j-_i];
    while (--_k)
      _d-=IIU[_i+_k][_j-_i-_k]*IU[_i][_k];
    IIU[_i][_j-_i]=_d;
    }
  }
}


//This function performs invertion of upper triangle matrix Calc y=U^-1.y
inline void bsubstitute_tdU(unsigned int n,double **U,double *y)
{
register double *_d,__d;
register unsigned int _i,_k;

_i=--n;
y+=n;
*y/=*U[n];
while (_i--)
  {
  y--;
  __d=*y;
  _k=n-_i;
  _d=U[_i];
  while (_k--)
    {
    y++, _d++;
    __d-=*y**_d;
    }
  y-=n-_i;
  *y=__d/(*U[_i]);
  }
}
//This function performs invertion of transposed upper triangle matrix Calc y=U^-T.y
inline void bsubstitute_tdUT(unsigned int n,double **U,double *y)
{
register double *_d,__d;
register unsigned int _i,_k;

for (_i=0;_i<n;y-=n-(++_i))
  {
  _d=U[_i];
  __d=*y/=*_d;
  _k=n-_i;
  while(--_k)
    {
    y++, _d++;
    *y-=__d**_d;
    }
  y++;
  }
}


//This function performs solvation of lower triangle linear equations system y=L^-1.y
inline void bsubstitute_tdIL(unsigned int n,double **IL,double *y)
{
register double *_d,__d;
register unsigned int _i,_k;

for (++y, _i=1; _i<n; _i++, y++)
  {
  _k=_i, __d=*y, y-=_i, _d=IL[_i];
  while (_k--) { __d-=*y**_d, _d++, y++; }
  *y=__d;
  }
}
//This function performs solvation of upper triangle linear equations system y=U^-1.y
//Note. The matrix U diagonal suggested to be unit
inline void bsubstitute_tdIU(unsigned int n,double **IU,double *y)
{
register double *_d,__d;
register unsigned int _i,_k;

_i=--n;
y+=n;
while (_i--)
  {
  y--;
  __d=*y;
  _k=n-_i;
  _d=IU[_i];
  while (_k--)
    {
    y++, _d++;
    __d-=*y**_d;
    }
  y-=n-_i;
  *y=__d;
  }
}
//This function performs solvation of transposed lower triangle linear equations system y=L^-1.y
inline void bsubstitute_tdILT(unsigned int n,double **IL,double *y)
{
register double *_d,__d;
register unsigned int _i,_k;
_i=n, y+=n;
while (_i--)
  {
  y--;
  __d=*y, _d=&IL[_i][_i];
  _k=_i; while (_k--) { y--, _d--, *y-=__d**_d; }
  y+=_i;
  }
}
//This function performs solvation of transposed upper triangle matrix y=U^-T.y
//Note. The matrix U diagonal suggested to be unit
inline void bsubstitute_tdIUT(unsigned int n,double **IU,register double *y)
{
register double *_d,__d;
register unsigned int _i,_k;

for (_i=0;_i<n;y-=n-(++_i))
  {
  _d=IU[_i], __d=*y;
  _k=n-_i; while(--_k) { y++, _d++, *y-=__d**_d; }
  y++;
  }
}

//This function solves equation (L.D.LT).y=y
inline void bsubstitute_tdLDLT(unsigned int n,double **tdLD,double *y)
{
register unsigned int _i;
bsubstitute_tdIL(n,tdLD,y);
_i=n, y+=n; while(_i--) { y--, *y/=tdLD[_i][_i]; }
bsubstitute_tdILT(n,tdLD,y);
}

//This function solves equation (UT.D.U).y=y
inline void bsubstitute_tdUTDU(unsigned int n,double **tdDU,double *y)
{
register unsigned int _i;
bsubstitute_tdIU(n,tdDU,y);
_i=n, y+=n; while (_i--) { y--, *y/=*tdDU[_i]; }
bsubstitute_tdIUT(n,tdDU,y);
}

//This function bsubstitude tdL matrix with given pivot 
inline void bsubstitute_pivoted_tdLDLT(unsigned int n,unsigned int *p,double **LD,double *y,double *_t)
{
register unsigned int _i;
_i=n; while (_i--) _t[_i]=y[p[_i]];
bsubstitute_tdLDLT(n,LD,_t);
_i=n; while (_i--) y[p[_i]]=_t[_i];
}
//This function bsubstitude tdU matrix with given pivot 
inline void bsubstitute_pivoted_tdUTDU(unsigned int n,unsigned int *p,double **UD,double *y,double *_t)
{
register unsigned int _i;
_i=n; while (_i--) _t[_i]=y[p[_i]];
bsubstitute_tdUTDU(n,UD,y);
_i=n; while (_i--) y[p[_i]]=_t[_i];
}



//This function inverts cholesky decomposed factor marix
void invert_tdUTDU(unsigned int n,double **IA,double **A,double *_t)
{
register unsigned int _i,_j;
_i=n;
while (_i--)
  {
 // if (!(_i%100)) printf("processing column %d of %d columns of given matrix\n",_i,A->ni);
  _j=n, _t+=n; while (_j--) { _t--, *_t=0.; }
  _t[_i]=1.;
  bsubstitute_tdUTDU(n,A,_t);
  _j=_i;
  do{ IA[_j][_i-_j]=_t[_j]; }while (_j--);
  }
}


//------------------- T H E   3 D   P A R T   ( t e n s o r s   p a r t ) --------------------------------

//This function transposes tensor A=B^T
//Note A and B can be the same
inline void transpose_tensor(t_tensor *A, t_tensor *B)
{
register double _d;
if (A==B)
  {
  _d=(*A)[0][1], (*A)[0][1]=(*A)[1][0], (*A)[1][0]=_d;  
  _d=(*A)[0][2], (*A)[0][2]=(*A)[2][0], (*A)[2][0]=_d;  
  _d=(*A)[1][2], (*A)[1][2]=(*A)[2][1], (*A)[2][1]=_d;  
  }
else
  {
  (*A)[0][0]=(*B)[0][0], (*A)[0][1]=(*B)[1][0], (*A)[0][2]=(*B)[2][0];
  (*A)[1][0]=(*B)[0][1], (*A)[1][1]=(*B)[1][1], (*A)[1][2]=(*B)[2][1];
  (*A)[2][0]=(*B)[0][2], (*A)[2][1]=(*B)[1][2], (*A)[2][2]=(*B)[2][2];
  }
}

//This function calculates product of b^T.R.a
inline double calc_bTRa(register t_vec *b,register t_tensor *R,register t_vec *a)
{
return b->i*(a->i*(*R)[0][0]+a->j*(*R)[0][1]+a->k*(*R)[0][2])+
       b->j*(a->i*(*R)[1][0]+a->j*(*R)[1][1]+a->k*(*R)[1][2])+
       b->k*(a->i*(*R)[2][0]+a->j*(*R)[2][1]+a->k*(*R)[2][2]);
}

//This function calculates tensor-vcector product a=T.b
inline void multiple_origin_tensor_origin_vec(register t_vec *c,register t_tensor *T,register t_vec *b)
{
c->i=(*T)[0x0][0x0]*b->i+(*T)[0x0][0x1]*b->j+(*T)[0x0][0x2]*b->k;
c->j=(*T)[0x1][0x0]*b->i+(*T)[0x1][0x1]*b->j+(*T)[0x1][0x2]*b->k;
c->k=(*T)[0x2][0x0]*b->i+(*T)[0x2][0x1]*b->j+(*T)[0x2][0x2]*b->k;
}
//This function calculates tensor-vcector product a=T^T.b
inline void multiple_transp_tensor_origin_vec(register t_vec *c,register t_tensor *T,register t_vec *b)
{
c->i=(*T)[0][0]*b->i+(*T)[1][0]*b->j+(*T)[2][0]*b->k;
c->j=(*T)[0][1]*b->i+(*T)[1][1]*b->j+(*T)[2][1]*b->k;
c->k=(*T)[0][2]*b->i+(*T)[1][2]*b->j+(*T)[2][2]*b->k;
}
//This function calculates tensor-vcector product aT=bT.R
inline void multiple_transp_vec_transp_tensor(register t_vec *c,register t_tensor *T,register t_vec *b)
{
c->i=(*T)[0][0]*b->i+(*T)[1][0]*b->j+(*T)[2][0]*b->k;
c->j=(*T)[0][1]*b->i+(*T)[1][1]*b->j+(*T)[2][1]*b->k;
c->k=(*T)[0][2]*b->i+(*T)[1][2]*b->j+(*T)[2][2]*b->k;
}
//This function calculates tensor.tensor product A=B.C
inline void multiple_origin_tensor_origin_tensor(register t_tensor *A,register t_tensor *B,register t_tensor *C)
{
(*A)[0][0]=(*B)[0][0]*(*C)[0][0]+(*B)[0][1]*(*C)[1][0]+(*B)[0][2]*(*C)[2][0];
(*A)[0][1]=(*B)[0][0]*(*C)[0][1]+(*B)[0][1]*(*C)[1][1]+(*B)[0][2]*(*C)[2][1];
(*A)[0][2]=(*B)[0][0]*(*C)[0][2]+(*B)[0][1]*(*C)[1][2]+(*B)[0][2]*(*C)[2][2];
(*A)[1][0]=(*B)[1][0]*(*C)[0][0]+(*B)[1][1]*(*C)[1][0]+(*B)[1][2]*(*C)[2][0];
(*A)[1][1]=(*B)[1][0]*(*C)[0][1]+(*B)[1][1]*(*C)[1][1]+(*B)[1][2]*(*C)[2][1];
(*A)[1][2]=(*B)[1][0]*(*C)[0][2]+(*B)[1][1]*(*C)[1][2]+(*B)[1][2]*(*C)[2][2];
(*A)[2][0]=(*B)[2][0]*(*C)[0][0]+(*B)[2][1]*(*C)[1][0]+(*B)[2][2]*(*C)[2][0];
(*A)[2][1]=(*B)[2][0]*(*C)[0][1]+(*B)[2][1]*(*C)[1][1]+(*B)[2][2]*(*C)[2][1];
(*A)[2][2]=(*B)[2][0]*(*C)[0][2]+(*B)[2][1]*(*C)[1][2]+(*B)[2][2]*(*C)[2][2];
}
//This function calculates tensor.tensor^T product A=B.C^T
inline void multiple_origin_tensor_transp_tensor(register t_tensor *A,register t_tensor *B,register t_tensor *C)
{
(*A)[0][0]=(*B)[0][0]*(*C)[0][0]+(*B)[0][1]*(*C)[0][1]+(*B)[0][2]*(*C)[0][2];
(*A)[0][1]=(*B)[0][0]*(*C)[1][0]+(*B)[0][1]*(*C)[1][1]+(*B)[0][2]*(*C)[1][2];
(*A)[0][2]=(*B)[0][0]*(*C)[2][0]+(*B)[0][1]*(*C)[2][1]+(*B)[0][2]*(*C)[2][2];
(*A)[1][0]=(*B)[1][0]*(*C)[0][0]+(*B)[1][1]*(*C)[0][1]+(*B)[1][2]*(*C)[0][2];
(*A)[1][1]=(*B)[1][0]*(*C)[1][0]+(*B)[1][1]*(*C)[1][1]+(*B)[1][2]*(*C)[1][2];
(*A)[1][2]=(*B)[1][0]*(*C)[2][0]+(*B)[1][1]*(*C)[2][1]+(*B)[1][2]*(*C)[2][2];
(*A)[2][0]=(*B)[2][0]*(*C)[0][0]+(*B)[2][1]*(*C)[0][1]+(*B)[2][2]*(*C)[0][2];
(*A)[2][1]=(*B)[2][0]*(*C)[1][0]+(*B)[2][1]*(*C)[1][1]+(*B)[2][2]*(*C)[1][2];
(*A)[2][2]=(*B)[2][0]*(*C)[2][0]+(*B)[2][1]*(*C)[2][1]+(*B)[2][2]*(*C)[2][2];
}

//This function solves 3x3 linear equations system with double pivoting
char gauss_solve_tensor(t_tensor *T,t_vec *b)
{
double _f;
unsigned int n, _i, _j, max_i, max_j, order[3];

order[0]=0, order[1]=1, order[2]=2;
//Stage 1. Do elimination
n=3;
while (--n)
  {
  //1. Find biggest value
  _i=n+1;
  max_j=max_i=n-1;
  while (_i--)
    {
    _j=n+1;
    while (_j--)
      if (fabs((*T)[_i][_j])>fabs((*T)[max_i][max_j])) 
        max_i=_i, max_j=_j; 	
    }
  if (fabs((*T)[max_i][max_j])<TINY) { ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }	
  //2. Exchange all if necessary
  if (max_i!=n) 
    {
    _j=n+1; while (_j--) { _f=(*T)[max_i][_j], (*T)[max_i][_j]=(*T)[n][_j], (*T)[n][_j]=_f; } 
    _f=((double*)b)[max_i], ((double*)b)[max_i]=((double*)b)[n], ((double*)b)[n]=_f;
    }
  if (max_j!=n)
    {
    _j=order[n], order[n]=max_j, order[max_j]=_j;
    _i=3; while (_i--) { _f=(*T)[_i][max_j], (*T)[_i][max_j]=(*T)[_i][n], (*T)[_i][n]=_f; }
    }
  //3. Normalize equation
  _f=(*T)[n][n], (*T)[n][n]=1.;
  if (fabs(_f)>TINY) { ((double*)b)[n]/=_f, _j=n; while (_j--) (*T)[n][_j]/=_f; }
  else return FALSE; //No solution exists
  //4. Update equations
  _i=n;
  while(_i--)
    {
    _f=(*T)[_i][n], (*T)[_i][n]=0.;
    if (fabs(_f)>TINY)
      {
      _j=n; while (_j--) (*T)[_i][_j]=(*T)[_i][_j]-_f*(*T)[n][_j];
      ((double*)b)[_i]=((double*)b)[_i]-_f*((double*)b)[n];
      }
    }
  }
//Stage 2. Do backsubstitution
((double*)b)[0]/=(*T)[0][0], ((double*)b)[1]-=((double*)b)[0]*(*T)[1][0], ((double*)b)[2]-=((double*)b)[0]*(*T)[2][0]+((double*)b)[1]*(*T)[2][1];
//Stage 3. Do reordering
if (order[0])
  {
  if (order[0]==1)
    {
    if (order[1])
      {// 1,2,0 - exchange 2-nd and 3-rd values then 1-st and 2-nd
      _f=((double*)b)[1], ((double*)b)[1]=((double*)b)[2], ((double*)b)[2]=_f;
      _f=((double*)b)[0], ((double*)b)[0]=((double*)b)[1], ((double*)b)[1]=_f;
      }
    else
      {// 1,0,2 - exchange 1-st and 2-nd values
      _f=((double*)b)[0], ((double*)b)[0]=((double*)b)[1], ((double*)b)[1]=_f;
      }
    }
  else
    {
    if (order[1])
      {// 2,1,0 - exchange 1-st and 3-rd values
      _f=((double*)b)[0], ((double*)b)[0]=((double*)b)[2], ((double*)b)[2]=_f;
      }
    else
      {// 2,0,1 - eschange 1-st and 2-nd values then 2-nd and 3-rd
      _f=((double*)b)[0], ((double*)b)[0]=((double*)b)[1], ((double*)b)[1]=_f;
      _f=((double*)b)[1], ((double*)b)[1]=((double*)b)[2], ((double*)b)[2]=_f;
      }
    }
  }
else
  {
  if (order[1]==1) ; // 0,1,2 - nothing to do
  else 
    {// 0,2,1 - exchange 2-nd and 3-rd values
    _f=((double*)b)[1], ((double*)b)[1]=((double*)b)[2], ((double*)b)[2]=_f;
    }
  }

return TRUE;
}

#define ROTATE(T1,T2)  g=T1, h=T2, T1=g-sn*(h+g*tau), T2=h+sn*(g-h*tau)

//This function calculates eigenvectors and eigenvalues of symmetric tensor
//Note. it uses Jacobi rotations until tol atchieved in off-diagonal elements
//Note2. See http://www.cyberguru.ru/programming/programming-theory/matrix-vectors-values-page6.html for details
unsigned int tensor_svd_weigenv(double tol,t_tensor *T,t_tensor *U,t_vec *S)
{
unsigned int niterations;
double g,h,tau,t,cs,sn,summ,tresh,theta;

//Init U and S
(*U)[0][0]=(*U)[1][1]=(*U)[2][2]=1.00;
(*U)[0][1]=(*U)[0][2]=(*U)[1][0]=(*U)[1][2]=(*U)[2][0]=(*U)[2][1]=0.00;
S->i=(*T)[0][0], S->j=(*T)[1][1], S->k=(*T)[2][2];
//Do rotations
niterations=0;
do{
  summ=fabs((*T)[0][1])+fabs((*T)[0][2])+fabs((*T)[1][2]); //Calc offdiagonal abs summ
  tresh=(niterations<3) ? 0.02*summ : 0.;  //Set scaling treshold for first 3 iterations
  //Dot [0][1]-[1][0] fragment
  g=100.*fabs((*T)[0][1]);
  if( (niterations>3)&&(g<tol) )  (*T)[0][1]=0.;
  else if (fabs((*T)[0][1])>tresh)
         {
         h=S->j-S->i;
         //calc value of t=s/c
         if ( fabs(h)+g==fabs(h) ) t=(*T)[0][1]/h;
         else
           {
           if ((theta=.5*h/(*T)[0][1])<0.) t=-1./(-theta+sqrt(1.+theta*theta));
           else                            t=+1./( theta+sqrt(1.+theta*theta));
           }
         // calc cos, sin, tau, and others, update diagonal.  Reduce (*T)[0][1]
         cs=1./sqrt(1.+t*t), sn=t*cs, tau=sn/(1.+cs), h=t*(*T)[0][1];
         S->i-=h, S->j+=h;
         (*T)[0][1]=0.;
         //rotation  iq<j<=n
         ROTATE((*T)[0][2],(*T)[1][2]);
         //update eigenvectors 1<=n
         ROTATE((*U)[0][0],(*U)[0][1]);
         ROTATE((*U)[1][0],(*U)[1][1]);
         ROTATE((*U)[2][0],(*U)[2][1]);
         }
  //Dot [1][2]-[2][1] fragment
  g=100.*fabs((*T)[1][2]);
  if( (niterations>3)&&(g<tol) )  (*T)[1][2]=0.;
  else if (fabs((*T)[1][2])>tresh)
         {
         h=S->k-S->j;
         //calc value of t=s/c
         if ( fabs(h)+g==fabs(h) ) t=(*T)[1][2]/h;
         else
           {
           if ((theta=.5*h/(*T)[1][2])<0.) t=-1./(-theta+sqrt(1.+theta*theta));
           else                            t=+1./( theta+sqrt(1.+theta*theta));
           }
         // calc cos, sin, tau, and others, update diagonal.  Reduce (*T)[1][2]
         cs=1./sqrt(1.+t*t), sn=t*cs, tau=sn/(1.+cs), h=t*(*T)[1][2];
         S->j-=h, S->k+=h;
         (*T)[1][2]=0.;
         //rotation  1<=j<ip
         ROTATE((*T)[0][1],(*T)[0][2]);
         //update eigenvectors 1<=n
         ROTATE((*U)[0][1],(*U)[0][2]);
         ROTATE((*U)[1][1],(*U)[1][2]);
         ROTATE((*U)[2][1],(*U)[2][2]);
         }
  //Dot [2][0]-[0][2] fragment
  g=100.*fabs((*T)[0][2]);
  if( (niterations>3)&&(g<tol) )  (*T)[0][2]=0.;
  else if (fabs((*T)[0][2])>tresh)
         {
         h=S->i-S->k;
         //calc value of t=s/c
         if ( fabs(h)+g==fabs(h) ) t=(*T)[0][2]/h;
         else
           {
           if ((theta=.5*h/(*T)[0][2])<0.) t=-1./(-theta+sqrt(1.+theta*theta));
           else                            t=+1./( theta+sqrt(1.+theta*theta));
           }
         // calc cos, sin, tau, and others, update diagonal.  Reduce (*T)[2][0]
         cs=1./sqrt(1.+t*t), sn=t*cs, tau=sn/(1.+cs), h=t*(*T)[0][2];
         S->k-=h, S->i+=h;
         (*T)[0][2]=0.;
         //rotation ip<j<iq
         ROTATE((*T)[1][2],(*T)[0][1]);
         //update eigenvectors 1<=n
         ROTATE((*U)[0][2],(*U)[0][0]);
         ROTATE((*U)[1][2],(*U)[1][0]);
         ROTATE((*U)[2][2],(*U)[2][0]);
         }
  //update diagonal and reinit z
  niterations++;
  }while (summ>tol);
//Exit
return niterations;
}
//The same as previous but eigenvectors are omitted
unsigned int tensor_svd_neigenv(double tol,t_tensor *T,t_vec *S)
{
unsigned int niterations;
double g,h,tau,t,cs,sn,summ,tresh,theta;

//Init U and S
S->i=(*T)[0][0], S->j=(*T)[1][1], S->k=(*T)[2][2];
//Do rotations
niterations=0;
while ((summ=fabs((*T)[0][1])+fabs((*T)[0][2])+fabs((*T)[1][2]))>tol) //Calc offdiagonal abs summ
  {
  tresh=(niterations<3) ? 0.02*summ : 0.00;  //Set scaling treshold for first 3 iterations
  //Dot [0][1]-[1][0] fragment
  g=100.*fabs((*T)[0][1]);
  if( (niterations>3)&&(g<tol) )  (*T)[0][1]=0.00;
  else if (fabs((*T)[0][1])>tresh)
         {
         h=S->j-S->i;
         //calc value of t=s/c
         if ( fabs(h)+g==fabs(h) ) t=(*T)[0][1]/h;
         else
           {
           if ((theta=0.50*h/(*T)[0][1])<0.00) t=-1./(-theta+sqrt(1.00+theta*theta));
           else                                t=+1./( theta+sqrt(1.00+theta*theta));
           }
         // calc cos, sin, tau, and others, update diagonal.  Reduce (*T)[0][1]
         cs=1./sqrt(1+t*t), sn=t*cs, tau=sn/(1.+cs), h=t*(*T)[0][1];
         S->i-=h, S->j+=h;
         (*T)[0][1]=0.00;
         //rotation  iq<j<=n
         ROTATE((*T)[0][2],(*T)[1][2]);
         }
  //Dot [1][2]-[2][1] fragment
  g=100.*fabs((*T)[1][2]);
  if( (niterations>3)&&(g<tol) )  (*T)[1][2]=0.00;
  else if (fabs((*T)[1][2])>tresh)
         {
         h=S->k-S->j;
         //calc value of t=s/c
         if ( fabs(h)+g==fabs(h) ) t=(*T)[1][2]/h;
         else
           {
           if ((theta=0.50*h/(*T)[1][2])<0.00) t=-1./(-theta+sqrt(1.00+theta*theta));
           else                                t=+1./( theta+sqrt(1.00+theta*theta));
           }
         // calc cos, sin, tau, and others, update diagonal.  Reduce (*T)[1][2]
         cs=1.00/sqrt(1.00+t*t), sn=t*cs, tau=sn/(1.00+cs), h=t*(*T)[1][2];
         S->j-=h, S->k+=h;
         (*T)[1][2]=0.00;
         //rotation  1<=j<ip
         ROTATE((*T)[0][1],(*T)[0][2]);
         }
  //Dot [2][0]-[0][2] fragment
  g=100.*fabs((*T)[0][2]);
  if( (niterations>3)&&(g<tol) )  (*T)[0][2]=0.00;
  else if (fabs((*T)[0][2])>tresh)
         {
         h=S->i-S->k;
         //calc value of t=s/c
         if ( fabs(h)+g==fabs(h) ) t=(*T)[0][2]/h;
         else
           {
           if ((theta=0.50*h/(*T)[0][2])<0.00) t=-1./(-theta+sqrt(1.00+theta*theta));
           else                                t=+1./( theta+sqrt(1.00+theta*theta));
           }
         // calc cos, sin, tau, and others, update diagonal.  Reduce (*T)[2][0]
         cs=1.00/sqrt(1.00+t*t), sn=t*cs, tau=sn/(1.00+cs), h=t*(*T)[0][2];
         S->k-=h, S->i+=h;
         (*T)[0][2]=0.00;
         //rotation ip<j<iq
         ROTATE((*T)[1][2],(*T)[0][1]);
         }
  //update diagonal and reinit z
  niterations++;
  }

//Exit
return niterations;
}

#undef ROTATE

//This function sort eightenvalues and corresponding eightenectors
void tensor_eigenv_sort(t_tensor *U,t_vec *S)
{
register double _d;
     if (S->i<S->j)
       {
       if (S->i<S->k)
         {//i<k
         if (S->j<S->k)
           {//i<j<k    i<->k
           _d=S->i, S->i=S->k, S->k=_d;
           if (*U)
             {
             _d=(*U)[0][0], (*U)[0][0]=(*U)[0][2], (*U)[0][2]=_d;
             _d=(*U)[1][0], (*U)[1][0]=(*U)[1][2], (*U)[1][2]=_d;
             _d=(*U)[2][0], (*U)[2][0]=(*U)[2][2], (*U)[2][2]=_d;
             }
           }
         else
           {//i<k<j    j<->k & k<->i
           _d=S->j, S->j=S->k, S->k=_d;
           _d=S->i, S->i=S->k, S->k=_d;
           if (*U)
             {
             _d=(*U)[0][1], (*U)[0][1]=(*U)[0][2], (*U)[0][2]=_d;
             _d=(*U)[1][1], (*U)[1][1]=(*U)[1][2], (*U)[1][2]=_d;
             _d=(*U)[2][1], (*U)[2][1]=(*U)[2][2], (*U)[2][2]=_d;
             _d=(*U)[0][0], (*U)[0][0]=(*U)[0][2], (*U)[0][2]=_d;
             _d=(*U)[1][0], (*U)[1][0]=(*U)[1][2], (*U)[1][2]=_d;
             _d=(*U)[2][0], (*U)[2][0]=(*U)[2][2], (*U)[2][2]=_d;
             }
           }
         }
       else
         {//k<i<j  i<->j
         _d=S->i, S->i=S->j, S->j=_d;
         if (*U)
           {
           _d=(*U)[0][0], (*U)[0][0]=(*U)[0][1], (*U)[0][1]=_d;
           _d=(*U)[1][0], (*U)[1][0]=(*U)[1][1], (*U)[1][1]=_d;
           _d=(*U)[2][0], (*U)[2][0]=(*U)[2][1], (*U)[2][1]=_d;
           }
         }
       }
else if (S->j<S->k)
       {// j<k
       if (S->i<S->k)
         {//j<i<k    i<->j  & i<->k
         _d=S->i, S->i=S->j, S->j=_d;
         _d=S->i, S->i=S->k, S->k=_d;
         if (*U)
           {
           _d=(*U)[0][0], (*U)[0][0]=(*U)[0][1], (*U)[0][1]=_d;
           _d=(*U)[1][0], (*U)[1][0]=(*U)[1][1], (*U)[1][1]=_d;
           _d=(*U)[2][0], (*U)[2][0]=(*U)[2][1], (*U)[2][1]=_d;
           _d=(*U)[0][0], (*U)[0][0]=(*U)[0][2], (*U)[0][2]=_d;
           _d=(*U)[1][0], (*U)[1][0]=(*U)[1][2], (*U)[1][2]=_d;
           _d=(*U)[2][0], (*U)[2][0]=(*U)[2][2], (*U)[2][2]=_d;
           }
         }
       else
         {//j<k<i  j<->k
         _d=S->j, S->j=S->k, S->k=_d;
         if (*U)
           {
           _d=(*U)[0][1], (*U)[0][1]=(*U)[0][2], (*U)[0][2]=_d;
           _d=(*U)[1][1], (*U)[1][1]=(*U)[1][2], (*U)[1][2]=_d;
           _d=(*U)[2][1], (*U)[2][1]=(*U)[2][2], (*U)[2][2]=_d;
           }
         }
       }
      else
       {//k<j<i   Nothing to do
       ;
       }
}

//This is univewrsializetion cover of jacobi Singular Value Decomposition
//NOTE. T sould have the lowest direction. T can be the same with A, but only if A symetric square matrix. In this case A will be destroed!
//It uses the rule that A^TA = V . S^2.V^T  and AA^T= U.S^2.U^T
unsigned int tensor_svd(char reorder,char eigenvectors,char symmetry_check,double tol,t_tensor *A,t_tensor *U,t_vec *S,t_tensor *V)
{
unsigned int _i;
t_tensor T;

if ( (symmetry_check)&&( ((*A)[0][1]!=(*A)[1][0])||((*A)[0][2]!=(*A)[2][0])||((*A)[1][2]!=(*A)[2][1]) ) )
  {
  //Nonsymmetrical case. Do symmetrization before evaluation T=A.A^T  svd(T)=U.S^2.U.
  T[0][0]=(*A)[0][0]*(*A)[0][0]+(*A)[0][1]*(*A)[0][1]+(*A)[0][2]*(*A)[0][2];
  T[0][1]=(*A)[0][0]*(*A)[1][0]+(*A)[0][1]*(*A)[1][1]+(*A)[0][2]*(*A)[1][2];
  T[0][2]=(*A)[0][0]*(*A)[2][0]+(*A)[0][1]*(*A)[2][1]+(*A)[0][2]*(*A)[2][2];
  T[1][0]=T[0][1];
  T[1][1]=(*A)[1][0]*(*A)[1][0]+(*A)[1][1]*(*A)[1][1]+(*A)[1][2]*(*A)[1][2];
  T[1][2]=(*A)[1][0]*(*A)[2][0]+(*A)[1][1]*(*A)[2][1]+(*A)[1][2]*(*A)[2][2];
  T[2][0]=T[0][2];
  T[2][1]=T[1][2];
  T[2][2]=(*A)[2][0]*(*A)[2][0]+(*A)[2][1]*(*A)[2][1]+(*A)[2][2]*(*A)[2][2];
  if (eigenvectors)
    {
    _i=tensor_svd_weigenv(tol,&T,U,S);
    if (reorder) tensor_eigenv_sort(U,S);
    }
  else
    {
    _i=tensor_svd_neigenv(tol,&T,S);
    if (reorder) tensor_eigenv_sort(0x0,S);
    }
  S->i= (S->i<0.00) ? +sqrt(-S->i) : +sqrt(+S->i);
  S->j= (S->j<0.00) ? +sqrt(-S->j) : +sqrt(+S->j);
  S->k= (S->k<0.00) ? +sqrt(-S->k) : +sqrt(+S->k);
  if (eigenvectors)
    {
    if (!*V) { ylib_errno=YERROR_USER; return (unsigned int)-1; }
    //NOTE. Constructing of V can be made in two ways, choose fastest if it is numerically stable
    if (S->k>tol) 
      {//The fastest way. It uses equation vi=(1/si)*A*ui
      (*V)[0][0]=((*A)[0][0]*(*U)[0][0]+(*A)[1][0]*(*U)[1][0]+(*A)[2][0]*(*U)[2][0])/S->i;
      (*V)[0][1]=((*A)[0][1]*(*U)[0][0]+(*A)[1][1]*(*U)[1][0]+(*A)[2][1]*(*U)[2][0])/S->i;
      (*V)[0][2]=((*A)[0][2]*(*U)[0][0]+(*A)[1][2]*(*U)[1][0]+(*A)[2][2]*(*U)[2][0])/S->i;
      (*V)[1][0]=((*A)[0][0]*(*U)[0][1]+(*A)[1][0]*(*U)[1][1]+(*A)[2][0]*(*U)[2][1])/S->j;
      (*V)[1][1]=((*A)[0][1]*(*U)[0][1]+(*A)[1][1]*(*U)[1][1]+(*A)[2][1]*(*U)[2][1])/S->j;
      (*V)[1][2]=((*A)[0][2]*(*U)[0][1]+(*A)[1][2]*(*U)[1][1]+(*A)[2][2]*(*U)[2][1])/S->j;
      (*V)[2][0]=((*A)[0][0]*(*U)[0][2]+(*A)[1][0]*(*U)[1][2]+(*A)[2][0]*(*U)[2][2])/S->k;
      (*V)[2][1]=((*A)[0][1]*(*U)[0][2]+(*A)[1][1]*(*U)[1][2]+(*A)[2][1]*(*U)[2][2])/S->k;
      (*V)[2][2]=((*A)[0][2]*(*U)[0][2]+(*A)[1][2]*(*U)[1][2]+(*A)[2][2]*(*U)[2][2])/S->k;
      }
    else
      {//Get numerically stable svd by reaiting with A^T.A
      if (!reorder) tensor_eigenv_sort(U,S);
      T[0][0]=(*A)[0][0]*(*A)[0][0]+(*A)[1][0]*(*A)[1][0]+(*A)[2][0]*(*A)[2][0];
      T[0][1]=(*A)[0][0]*(*A)[0][1]+(*A)[1][0]*(*A)[1][1]+(*A)[2][0]*(*A)[2][1];
      T[0][2]=(*A)[0][0]*(*A)[0][2]+(*A)[1][0]*(*A)[1][2]+(*A)[2][0]*(*A)[2][2];
      T[1][0]=T[0][1];
      T[1][1]=(*A)[0][1]*(*A)[0][1]+(*A)[1][1]*(*A)[1][1]+(*A)[2][1]*(*A)[2][1];
      T[1][2]=(*A)[0][1]*(*A)[0][2]+(*A)[1][1]*(*A)[1][2]+(*A)[2][1]*(*A)[2][2];
      T[2][0]=T[0][2];
      T[2][1]=T[1][2];
      T[2][2]=(*A)[0][2]*(*A)[0][2]+(*A)[1][2]*(*A)[1][2]+(*A)[2][2]*(*A)[2][2];
      _i=tensor_svd_weigenv(tol,&T,V,S);
      S->i= (S->i<0.00) ? +sqrt(-S->i) : +sqrt(+S->i);
      S->j= (S->j<0.00) ? +sqrt(-S->j) : +sqrt(+S->j);
      S->k= (S->k<0.00) ? +sqrt(-S->k) : +sqrt(+S->k);
      tensor_eigenv_sort(V,S);
      transpose_tensor(V,V);
      //Resolve unclearity with the sign: sign(U^T.S^-1.A)==sign(V), but S is always positive
      if (((*U)[0][0]*(*A)[0][0]+(*U)[1][0]*(*A)[1][0]+(*U)[2][0]*(*A)[2][0]>0.)!=((*V)[0][0]>0.)) (*U)[0][0]=-(*U)[0][0], (*U)[1][0]=-(*U)[1][0], (*U)[2][0]=-(*U)[2][0]; 
      if (((*U)[0][1]*(*A)[0][1]+(*U)[1][1]*(*A)[1][1]+(*U)[2][1]*(*A)[2][1]>0.)!=((*V)[1][1]>0.)) (*U)[0][1]=-(*U)[0][1], (*U)[1][1]=-(*U)[1][1], (*U)[2][1]=-(*U)[2][1]; 
      if (((*U)[0][2]*(*A)[0][2]+(*U)[1][2]*(*A)[1][2]+(*U)[2][2]*(*A)[2][2]>0.)!=((*V)[2][2]>0.)) (*U)[0][2]=-(*U)[0][2], (*U)[1][2]=-(*U)[1][2], (*U)[2][2]=-(*U)[2][2];  
      }
    }
  }
else
  {
  //Symmetrical case
  T[0][0]=(*A)[0][0];
  T[0][1]=(*A)[0][1];
  T[0][2]=(*A)[0][2];
  T[1][0]=T[0][1];
  T[1][1]=(*A)[1][1];
  T[1][2]=(*A)[1][2];
  T[2][0]=T[0][2];
  T[2][1]=T[1][2];
  T[2][2]=(*A)[2][2];
  if (eigenvectors)
    {
    _i=tensor_svd_weigenv(tol,&T,U,S);
    if (reorder) tensor_eigenv_sort(U,S);
    if ( (*V)&&(*U!=*V) )
      {
      (*V)[0][0]=(*U)[0][0];
      (*V)[0][1]=(*U)[1][0];
      (*V)[0][2]=(*U)[2][0];
      (*V)[1][0]=(*U)[0][1];
      (*V)[1][1]=(*U)[1][1];
      (*V)[1][2]=(*U)[2][1];
      (*V)[2][0]=(*U)[0][2];
      (*V)[2][1]=(*U)[1][2];
      (*V)[2][2]=(*U)[2][2];
      }
    }
  else
    {
    _i=tensor_svd_neigenv(tol,&T,S);
    if (reorder) tensor_eigenv_sort(0x0,S);
    }
  }
//Return amount of steps required to atchieve requested tolerance
return _i;
}

//This function makes eighten-values positive
inline void positive_svd(t_tensor *U,t_vec *S,t_tensor *V)
{
if (S->i<0.) { S->i=-S->i; (*V)[0][0]=-(*V)[0][0]; (*V)[0][1]=-(*V)[0][1]; (*V)[0][2]=-(*V)[0][2]; }
if (S->j<0.) { S->j=-S->j; (*V)[1][0]=-(*V)[1][0]; (*V)[1][1]=-(*V)[1][1]; (*V)[1][2]=-(*V)[1][2]; }
if (S->k<0.) { S->k=-S->k; (*V)[2][0]=-(*V)[2][0]; (*V)[2][1]=-(*V)[2][1]; (*V)[2][2]=-(*V)[2][2]; }
}


//------------------- T H E   4 D   P A R T   ( q t e n s o r s   p a r t ) --------------------------------

//This function multiple two qutensors into 3D tensor
inline void multiple_qtensor_tqtensor(t_tensor *t,t_qtensor *qt,t_qtensor *qtt)
{
(*t)[0][0]=(*qt)[0][0]*(*qtt)[0][0]+(*qt)[0][1]*(*qtt)[0][1]+(*qt)[0][2]*(*qtt)[0][2]+(*qt)[0][3]*(*qtt)[0][3];
(*t)[0][1]=(*qt)[0][0]*(*qtt)[1][0]+(*qt)[0][1]*(*qtt)[1][1]+(*qt)[0][2]*(*qtt)[1][2]+(*qt)[0][3]*(*qtt)[1][3];
(*t)[0][2]=(*qt)[0][0]*(*qtt)[2][0]+(*qt)[0][1]*(*qtt)[2][1]+(*qt)[0][2]*(*qtt)[2][2]+(*qt)[0][3]*(*qtt)[2][3];
(*t)[1][0]=(*qt)[1][0]*(*qtt)[0][0]+(*qt)[1][1]*(*qtt)[0][1]+(*qt)[1][2]*(*qtt)[0][2]+(*qt)[1][3]*(*qtt)[0][3];
(*t)[1][1]=(*qt)[1][0]*(*qtt)[1][0]+(*qt)[1][1]*(*qtt)[1][1]+(*qt)[1][2]*(*qtt)[1][2]+(*qt)[1][3]*(*qtt)[1][3];
(*t)[1][2]=(*qt)[1][0]*(*qtt)[2][0]+(*qt)[1][1]*(*qtt)[2][1]+(*qt)[1][2]*(*qtt)[2][2]+(*qt)[1][3]*(*qtt)[2][3];
(*t)[2][0]=(*qt)[2][0]*(*qtt)[0][0]+(*qt)[2][1]*(*qtt)[0][1]+(*qt)[2][2]*(*qtt)[0][2]+(*qt)[2][3]*(*qtt)[0][3];
(*t)[2][1]=(*qt)[2][0]*(*qtt)[1][0]+(*qt)[2][1]*(*qtt)[1][1]+(*qt)[2][2]*(*qtt)[1][2]+(*qt)[2][3]*(*qtt)[1][3];
(*t)[2][2]=(*qt)[2][0]*(*qtt)[2][0]+(*qt)[2][1]*(*qtt)[2][1]+(*qt)[2][2]*(*qtt)[2][2]+(*qt)[2][3]*(*qtt)[2][3];
}

//This function solves 4x4 linear equations system with iterative Gauss-Seindel method
char gauss_seindel_solve_qtensor(t_lvec *x,t_qtensor *Q,t_lvec *b,unsigned int niter)
{
double _d, dx, norm;
unsigned int _i, _j, _k, best_j;

//Optimizing diagonale entrance wih iterative swapping
norm=sqrt((*Q)[0][0]*(*Q)[0][0]+(*Q)[0][1]*(*Q)[0][1]+(*Q)[0][2]*(*Q)[0][2]+(*Q)[0][3]*(*Q)[0][3]+
          (*Q)[1][0]*(*Q)[1][0]+(*Q)[1][1]*(*Q)[1][1]+(*Q)[1][2]*(*Q)[1][2]+(*Q)[1][3]*(*Q)[1][3]+
          (*Q)[2][0]*(*Q)[2][0]+(*Q)[2][1]*(*Q)[2][1]+(*Q)[2][2]*(*Q)[2][2]+(*Q)[2][3]*(*Q)[2][3]+
          (*Q)[3][0]*(*Q)[3][0]+(*Q)[3][1]*(*Q)[3][1]+(*Q)[3][2]*(*Q)[3][2]+(*Q)[3][3]*(*Q)[3][3]);
_k=DOMINANT_DIAG_ITERATIONS;
do{
  _i=4;
  while (--_i)
    {
    dx=0., _j=_i;
    while (_j--)
      if ((fabs((*Q)[_i][_j])>norm*TINY)&&(fabs((*Q)[_j][_i])>norm*TINY))
        if ((_d=1./fabs((*Q)[_i][_i]+SMALL2)+1./fabs((*Q)[_j][_j]+SMALL2)-1./fabs((*Q)[_i][_j]+SMALL2)-1./fabs((*Q)[_j][_i])+SMALL2)>dx)
          {
          dx=_d, best_j=_j;
          }
    if (!(dx)) break;
    _d=(*Q)[_i][0], (*Q)[_i][0]=(*Q)[best_j][0], (*Q)[best_j][0]=_d;  
    _d=(*Q)[_i][1], (*Q)[_i][1]=(*Q)[best_j][1], (*Q)[best_j][1]=_d;  
    _d=(*Q)[_i][2], (*Q)[_i][2]=(*Q)[best_j][2], (*Q)[best_j][2]=_d;  
    _d=(*Q)[_i][3], (*Q)[_i][3]=(*Q)[best_j][3], (*Q)[best_j][3]=_d;  
    _d=((double*)b)[_i], ((double*)b)[_i]=((double*)b)[best_j], ((double*)b)[best_j]=_d;
    }
  }while (--_k);
if ( (fabs((*Q)[0][0])<norm*SMALL2)||(fabs((*Q)[1][1])<norm*SMALL2)||(fabs((*Q)[2][2])<norm*SMALL2) ) return FALSE; //No numerically feasiable evaluation is possible
//Solve the system
((double*)x)[0]=((double*)b)[0], ((double*)x)[1]=((double*)b)[1], ((double*)x)[2]=((double*)b)[2], ((double*)x)[3]=((double*)b)[3];
do{
  for (dx=0., _i=0;_i<4;_i++)
    {
    _d=0., _j=4; while (--_j!=_i) { _d+=(*Q)[_i][_j]*((double*)x)[_j]; }
                 while (_j--)     { _d+=(*Q)[_i][_j]*((double*)x)[_j]; }
    _d=(((double*)b)[_i]-_d)/(*Q)[_i][_i];
    dx+=fabs(((double*)x)[_i]-_d), ((double*)x)[_i]=_d;
    }
  if (dx<norm*TINY) return TRUE; //Converged
  }while (--niter);
return NTNF; //Hasn't converged but results exists
}


//---------------------- S O M E    A C C E L E R A T I O N    P A R T --------------------------------------

//This function calculates 2x2 determinant
inline double calc_det2x2(double v00,double v01,
                          double v10,double v11)
{
return v00*v11-v10*v01;  // 3 flops
}

//This function calculates 3x3 determinant
inline double calc_det3x3(double v00,double v01,double v02,
                          double v10,double v11,double v12,
                          double v20,double v21,double v22 )
{
return v00*(v11*v22-v21*v12)-v01*(v10*v22-v12*v20)+v02*(v10*v21-v11*v20); // 14 flops
}

//This function calculates 4x4 determinant
inline double calc_det4x4(double v00,double v01,double v02,double v03,
                          double v10,double v11,double v12,double v13,
                          double v20,double v21,double v22,double v23,
                          double v30,double v31,double v32,double v33 )
{
return v00*(v11*(v22*v33-v32*v23)-v12*(v21*v33-v31*v23)+v13*(v21*v32-v31*v22))-  // 16 flops
       v01*(v10*(v22*v33-v32*v23)-v12*(v20*v33-v30*v23)+v13*(v20*v32-v30*v22))+
       v02*(v10*(v21*v33-v31*v23)-v11*(v20*v33-v30*v23)+v13*(v20*v31-v30*v21))-
       v03*(v10*(v21*v32-v31*v22)-v11*(v20*v32-v30*v22)+v12*(v20*v31-v30*v21));  //total 63 flops
}

//This function calculates 5x5 determinant
inline double calc_det5x5(double v00,double v01,double v02,double v03,double v04,
                          double v10,double v11,double v12,double v13,double v14,
                          double v20,double v21,double v22,double v23,double v24,
                          double v30,double v31,double v32,double v33,double v34,
                          double v40,double v41,double v42,double v43,double v44)
{
return v00*calc_det4x4(v11,v12,v13,v14, v21,v22,v23,v24, v31,v32,v33,v34, v41,v42,v43,v44 )- //total 64 flops
       v01*calc_det4x4(v10,v12,v13,v14, v20,v22,v23,v24, v30,v32,v33,v34, v40,v42,v43,v44 )+
       v02*calc_det4x4(v10,v11,v13,v14, v20,v21,v23,v24, v30,v31,v33,v34, v40,v41,v43,v44 )-
       v03*calc_det4x4(v10,v11,v12,v14, v20,v21,v22,v24, v30,v31,v32,v34, v40,v41,v42,v44 )+
       v04*calc_det4x4(v10,v11,v12,v13, v20,v21,v22,v23, v30,v31,v32,v33, v40,v41,v42,v43 ); //total 320 flps
}

//This function calculates 6x6 determinant
inline double calc_det6x6(double v00,double v01,double v02,double v03,double v04,double v05,
                          double v10,double v11,double v12,double v13,double v14,double v15,
                          double v20,double v21,double v22,double v23,double v24,double v25,
                          double v30,double v31,double v32,double v33,double v34,double v35,
                          double v40,double v41,double v42,double v43,double v44,double v45,
                          double v50,double v51,double v52,double v53,double v54,double v55)
{
return v00*calc_det5x5(v11,v12,v13,v14,v15, v21,v22,v23,v24,v25, v31,v32,v33,v34,v35, v41,v42,v43,v44,v45, v51,v52,v53,v54,v55 )- //total 321 flops-
       v01*calc_det5x5(v10,v12,v13,v14,v15, v20,v22,v23,v24,v25, v30,v32,v33,v34,v35, v40,v42,v43,v44,v45, v50,v52,v53,v54,v55 )+
       v02*calc_det5x5(v10,v11,v13,v14,v15, v20,v21,v23,v24,v25, v30,v31,v33,v34,v35, v40,v41,v43,v44,v45, v50,v51,v53,v54,v55 )-
       v03*calc_det5x5(v10,v11,v12,v14,v15, v20,v21,v22,v24,v25, v30,v31,v32,v34,v35, v40,v41,v42,v44,v45, v50,v51,v52,v54,v55 )+
       v04*calc_det5x5(v10,v11,v12,v13,v15, v20,v21,v22,v23,v25, v30,v31,v32,v33,v35, v40,v41,v42,v43,v45, v50,v51,v52,v53,v55 )-
       v05*calc_det5x5(v10,v11,v12,v13,v14, v20,v21,v22,v23,v24, v30,v31,v32,v33,v34, v40,v41,v42,v43,v44, v50,v51,v52,v53,v54 ); //total 1926 flops
}
//We need 6x6 determinant to calculate 3D RT, but all future should be processed with more efficient techniques.


