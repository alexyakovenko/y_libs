#include "y_smatrix.h"

#define SMATRIX_REALLOC_QUANT 32768

extern unsigned int ylib_errno;

//This function show sparse matrix
void show_smatrix(t_smatrix *sA)
{
register unsigned int _i,_j;
printf("ni=%d, nj=%d, nnz=%d\ncontent is:\n",sA->ni,sA->nj,sA->nnz);
for (_i=0;_i<sA->ni;_i++)
  {
  printf("\n i=%d -> ",_i);
  for (_j=sA->i[_i];_j<sA->i[_i+0x1];_j++)
    printf(" %d:%f ",sA->j[_j],sA->d[_j]);
  }
printf("\n");
}

//Allocate memory to smatrix structure
t_smatrix *alloc_smatrix(unsigned char ctype,unsigned int ni,unsigned int nj,unsigned int nnz)
{
t_smatrix *sA=0x0;

if (!(sA=(t_smatrix*)calloc(sizeof(t_smatrix),0x1)))
  return FALSE;
sA->nnz=nnz;
if ((sA->ctype=ctype))
  {//Do column-wise allocation
  sA->nj=ni;
  sA->ni=nj;
  }
else
  {//Do row-wise allocation
  sA->ni=ni;
  sA->nj=nj;
  }
if (!(sA->i=(unsigned int*)malloc(sizeof(unsigned int)*(sA->ni+0x1)))) { free(sA); return FALSE; }
if (sA->nnz)
  {
  if ( (!(sA->j=(unsigned int*)malloc(sizeof(unsigned int)*sA->nnz)))   ||
       (!(sA->d=(double*)malloc(sizeof(double)*sA->nnz)))                )
    {
    free_smatrix(sA);
    return FALSE;
    }
  }
return sA;
}

//This function resize sparse matrix structure
char realloc_smatrix(unsigned int new_i,unsigned int new_nnz,t_smatrix *smatrix)
{
register unsigned int *_i=0x0, *_j=0x0;
register double *_d=0x0;
//realloc memory first
if (new_i!=(unsigned int)-1)
  {
  if (!(_i=(unsigned int*)realloc(smatrix->i,sizeof(unsigned int)*(new_i+0x1)))) goto ERROR;
  }
if (new_nnz!=(unsigned int)-1)
  {
  if ( (!(_j=(unsigned int*)realloc(smatrix->j,sizeof(unsigned int)*new_nnz)))||(!(_d=(double*)realloc(smatrix->d,sizeof(double)*new_nnz))) ) goto ERROR;
  }
//Mount new smatrix
if (new_i!=(unsigned int)-1)
  smatrix->i=_i;
if (new_nnz!=(unsigned int)-1)
  {
  smatrix->j=_j;
  smatrix->d=_d;
  }
return TRUE;
ERROR:
if (_i) smatrix->i=_i;
if (_j) smatrix->j=_j;
if (_d) smatrix->d=_d;
return FALSE;
}

//This function free memory allocated for sparse matrix
void free_smatrix(t_smatrix *sA)
{
if (sA)
  {
  if (sA->i) free(sA->i);
  if (sA->j) free(sA->j);
  if (sA->d) free(sA->d);
  free(sA);
  }
}

//This function read sparse matrix from file matrix
t_smatrix *read_smatrix(FILE *in)
{
t_smatrix *smatrix;
if (!(smatrix=(t_smatrix*)calloc(sizeof(t_smatrix),0x1))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
if ( (fread(&smatrix->ctype,sizeof(unsigned char),0x1,in)!=0x1)||(fread(&smatrix->ni,sizeof(unsigned int),0x1,in)!=0x1)||
     (fread(&smatrix->nj,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&smatrix->nnz,sizeof(unsigned int),0x1,in)!=0x1)    )
  { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if ( (!(smatrix->i=(unsigned int*)malloc(sizeof(unsigned int)*smatrix->ni+0x1)))||
     (!(smatrix->j=(unsigned int*)malloc(sizeof(unsigned int)*smatrix->nnz)))   ||
     (!(smatrix->d=(double*)malloc(sizeof(double)*smatrix->nnz)))                ) { free(smatrix); goto LABEL_MEMORY_ERROR; }
if ( (fread(smatrix->i,sizeof(unsigned int),smatrix->ni+0x1,in)!=smatrix->ni+0x1)||
     (fread(smatrix->j,sizeof(unsigned int),smatrix->nnz,in)!=smatrix->nnz)      ||
     (fread(smatrix->d,sizeof(double),smatrix->nnz,in)!=smatrix->nnz)             ) { free_smatrix(smatrix); goto LABEL_IO_ERROR; }
return smatrix;
}

//This function wtrites sparse matrix into the file
char write_smatrix(FILE *out,t_smatrix *smatrix)
{
if ( (fwrite(&smatrix->ctype,sizeof(unsigned char),0x1,out)!=0x1)                  ||
     (fwrite(&smatrix->ni,sizeof(unsigned int),0x1,out)!=0x1)                      ||
     (fwrite(&smatrix->nj,sizeof(unsigned int),0x1,out)!=0x1)                      ||
     (fwrite(&smatrix->nnz,sizeof(unsigned int),0x1,out)!=0x1)                     ||
     (fwrite(smatrix->i,sizeof(unsigned int),smatrix->ni+0x1,out)!=smatrix->ni+0x1)||
     (fwrite(smatrix->j,sizeof(unsigned int),smatrix->nnz,out)!=smatrix->nnz)      ||
     (fwrite(smatrix->d,sizeof(double),smatrix->nnz,out)!=smatrix->nnz)            ) { ylib_errno=YERROR_IO; return FALSE; }
else return TRUE;
}


//This function convert dense matrix into sparse matrix
t_smatrix *dmatrix_to_smatrix(register double zero_cutoff,t_dmatrix *dA)
{
register unsigned int _i,_j;
t_smatrix *sA;

if (!(sA=(t_smatrix*)malloc(sizeof(t_smatrix)))) return FALSE;
sA->ni=dA->ni;
sA->nj=dA->nj;
if (!(sA->i=(unsigned int*)malloc(sizeof(unsigned int)*sA->ni))) { free(sA); return FALSE; }
sA->i[sA->nnz=0x0]=0x0;
for (_i=0x0;_i<dA->ni;sA->i[++_i]=sA->nnz)
  for (_j=0x0;_j<dA->nj;_j++)
    if (fabs(dA->d[_i][_j])>zero_cutoff)
      sA->nnz++;
if (!(sA->j=(unsigned int*)malloc(sizeof(unsigned int)*sA->nnz)))  { free(sA->i); free(sA); return FALSE; }
if (!(sA->d=(double*)malloc(sizeof(double)*sA->nnz))) { free(sA->j); free(sA->i); free(sA); return FALSE; }
sA->nnz=0x0;
for (_i=0x0;_i<dA->ni;_i++)
  for (_j=0x0;_j<dA->nj;_j++)
    if (fabs(dA->d[_i][_j])>zero_cutoff)
      {
      sA->j[sA->nnz]=_j;
      sA->d[sA->nnz]=dA->d[_i][_j];
      sA->nnz++;
      }
return sA;
}
//This function convert sparse matrix into dense matrix
t_dmatrix *smatrix_to_dmatrix(t_smatrix *sA)
{
register unsigned int _i,_j,_k;
t_dmatrix *dA;

if (!(dA=(t_dmatrix*)alloc_dmatrix(sA->ni,sA->nj))) return FALSE;
for (_i=0;_i<sA->ni;_i++)
  {
  for (_k=0,_j=sA->i[_i];_j<sA->i[_i+0x1];_j++)
    {
    while (_k<sA->j[_j]) dA->d[_i][_k++]=0.00;
    dA->d[_i][_k++]=sA->d[_j];
//    if (sA->d[_j]>TINY)  dA->d[_i][_k++]=sA->d[_j];
//    else                 dA->d[_i][_k++]=0.00;
    }
  while (_k<dA->nj) dA->d[_i][_k++]=0.00;
  }
return dA;
}
//This function convert sparse matrix into triangular dense matrix
t_dmatrix *smatrix_to_tdmatrix(t_smatrix *sA)
{
register unsigned int _i,_j;
t_dmatrix *tdU;

if ( (sA->ni!=sA->nj)||(!(tdU=(t_dmatrix*)alloc_tdmatrix('U',sA->ni))) ) return FALSE;
for (_i=0x0;_i<sA->ni;_i++)
  {
  for (_j=_i;_j<tdU->ni;_j++)
    tdU->d[_i][_j-_i]=0.00;
  for (_j=sA->i[_i];_j<sA->i[_i+0x1];_j++)
    if (sA->j[_j]>=_i)
      tdU->d[_i][sA->j[_j]-_i]=sA->d[_j];
  }
return tdU;
}


//This function sorts entries in a sparse matrix
//Note. This function switch its strategy depending subproblem size
inline char order_sparse_matrix(t_smatrix *sA)
{
register unsigned int _i,_j;
_i=sA->ni;
sA->i+=_i;
while (_i--)
  {
  _j=*sA->i;
  sA->i--;
  _j-=*sA->i;
  if (!(ud_qsort(_j,&sA->j[*sA->i],&sA->d[*sA->i]))) return FALSE;
  }
return TRUE;
}

//This function exchanges row and columns compression and vice versa. Technically it looks alike transposition.
void row_column_permutate(t_smatrix *sA,t_smatrix *sAT)
{
register unsigned int _i,_j,_k;
//Gets matrix statistics 2(N)
sAT->i++;
_i=sA->nnz;
while (_i--)
  sAT->i[sA->j[_i]]++;
//Make AT map 2(M)
_i=sAT->ni-0x1;
sAT->i[_i]=sAT->nnz-sAT->i[_i];
while (--_i)
  sAT->i[_i]=sAT->i[_i+0x1]-sAT->i[_i];
*sAT->i=0x0;
//Fill AT accordingly with its map 1(N)
for (_j=0x0,_i=0x0;_i<sA->ni;_i++)
  for (;_j<sA->i[_i+0x1];_j++)
    {
    _k=sAT->i[sA->j[_j]]++;
    sAT->j[_k]=_i;
    sAT->d[_k]=sA->d[_j];
    }
sAT->i--;
// Transposition completed with total complexity 2N+2M
}

//This function exchanges row and columns compression and vice versa. Technically it looks alike transposition.
void row_column_portrain_permutate(t_smatrix *sA,t_smatrix *sAT)
{
register unsigned int _i,_j;
//Gets matrix statistics 2(N)
sAT->i++;
_i=sA->nnz;
while (_i--)
  sAT->i[sA->j[_i]]++;
//Make AT map 2(M)
_i=sAT->ni-0x1;
sAT->i[_i]=sAT->nnz-sAT->i[_i];
while (--_i)
  sAT->i[_i]=sAT->i[_i+0x1]-sAT->i[_i];
*sAT->i=0x0;
//Fill AT accordingly with its map 1(N)
for (_j=0x0,_i=0x0;_i<sA->ni;_i++)
  for (;_j<sA->i[_i+0x1];_j++)
    sAT->j[sAT->i[sA->j[_j]]++]=_i;
sAT->i--;
// Transposition completed with total complexity 2N+2M
}

//This function computes transpose of sparse matrix: A^T=Trasnsp(A);
//Note. If 'ordered' flag TRUE the transposed matrix is ordered with 2*o(N^2) complexity else it is unordered with 1*o(N^2) complexity.
t_smatrix *transpose_smatrix(t_smatrix *sA)
{
t_smatrix *sAT=0x0;
if (!(sAT=(t_smatrix*)malloc(sizeof(t_smatrix)))) return FALSE;
sAT->ni=sA->nj;
sAT->nj=sA->ni;
sAT->nnz=sA->nnz;
if (!(sAT->i=(unsigned int*)calloc(sizeof(unsigned int),sAT->ni+0x1)))              { free(sAT); return FALSE; }
if (!(sAT->j=(unsigned int*)malloc(sizeof(unsigned int)*sAT->nnz)))   { free(sAT->i); free(sAT); return FALSE; }
if (!(sAT->d=(double*)malloc(sizeof(double)*sAT->nnz))) { free(sAT->j); free(sAT->i); free(sAT); return FALSE; }
row_column_permutate(sA,sAT);
sAT->ctype=sA->ctype;
return sAT;
}

//This function computes transpose of sparse matrix: A^T=Trasnsp(A);
//Note. If 'ordered' flag TRUE the transposed matrix is ordered with 2*o(N^2) complexity else it is unordered with 1*o(N^2) complexity.
t_smatrix *transpose_smatrix_potrain(t_smatrix *sA)
{
t_smatrix *sAT=0x0;
if (!(sAT=(t_smatrix*)malloc(sizeof(t_smatrix)))) return FALSE;
sAT->ni=sA->nj;
sAT->nj=sA->ni;
sAT->nnz=sA->nnz;
if (!(sAT->i=(unsigned int*)calloc(sizeof(unsigned int),sAT->ni+0x1)))            { free(sAT); return FALSE; }
if (!(sAT->j=(unsigned int*)malloc(sizeof(unsigned int)*sAT->nnz))) { free(sAT->i); free(sAT); return FALSE; }
row_column_portrain_permutate(sA,sAT);
sAT->d=0x0;
sAT->ctype=sA->ctype;
return sAT;
}

//This block of functions do sparse matrix - vector multiplication

//This function calc sparse matrix - dense vector multiplication. Matrix should have row-wise compression
//c=A.b
void multiple_origin_smatrix_dvector(double *c,t_smatrix *sA,double *b)
{
register unsigned int _i,_j;
_i=sA->ni;
c+=_i;
sA->i+=_i;
sA->j+=sA->nnz;
sA->d+=sA->nnz;
while (_i--)
  {
  _j=*sA->i;
  sA->i--;
  *(--c)=0.00; //init inline
  while (*sA->i!=_j--)
    {
    sA->d--;
    sA->j--;
    *c+=*sA->d*b[*sA->j];
    }
  }
}
void multiple_origin_stL_dvector(double *c,t_smatrix *stL,double *b)
{
register unsigned int _i,_j;
set_vect(stL->ni,c,0.);
if (stL->nnz)
  {
  for (_i=0;_i<stL->ni;_i++)
    {
    if (!(_j=stL->i[_i+1]-stL->i[_i])) continue;
    while (--_j)
      {
      c[_i]+=*stL->d*b[*stL->j], c[*stL->j]+=*stL->d*b[_i], stL->d++, stL->j++;
      }
    //process the last element especially as it can be diagonal
    if (*stL->j==_i) c[_i]+=*stL->d*b[_i];
    else { c[_i]+=*stL->d*b[*stL->j], c[*stL->j]+=*stL->d*b[_i]; }
    stL->d++, stL->j++;
    }
  stL->j-=stL->nnz, stL->d-=stL->nnz;
  }
}
//This function calc sparse matrix - dense vector multiplication. Matrix should have row-wise compression
//c=AT.b
void multiple_transp_smatrix_dvector(double *c,t_smatrix *sAT,double *b)
{
register unsigned int _i,_j;
set_vect(sAT->ni,c,0.);
_i=sAT->ni;
b+=_i;
sAT->i+=_i;
sAT->j+=sAT->nnz;
sAT->d+=sAT->nnz;
while (_i--)
  {
  _j=*sAT->i;
  sAT->i--;
  b--;
  while (*sAT->i!=_j--)
    {
    sAT->d--;
    sAT->j--;
    c[*sAT->j]+=*sAT->d**b;
    }
  }
}

//This function calculates c=LI.D.LIT.b product
//NOTE. b vector is destroed during calculations
void multiple_sLDLT_dvector(double *c,t_smatrix *sLD,double *b)
{
register unsigned int _i,_j;
set_vect(sLD->ni,c,0.);
if (sLD->nnz)
  {
  //Calculate c=LIT.b
  _i=sLD->ni, sLD->j+=sLD->nnz, sLD->d+=sLD->nnz, b+=sLD->ni;
  while (_i--)
    {
    b--;
    if (!(_j=sLD->i[_i+1]-sLD->i[_i])) { *b=0.; continue; }
    //Process the first element especially as it can be the diagonal
    sLD->d--, sLD->j--;
    if (*sLD->j!=_i) { c[_i]+=*sLD->d**b; }
    while (--_j)
      {
      sLD->d--, sLD->j--, c[*sLD->j]+=*sLD->d**b;
      }
    }
  //Calculate b=D.c
  _i=sLD->ni, b+=sLD->ni, c+=sLD->ni; while (_i--) { b--, c--, *b+=*c; if (sLD->j[sLD->i[_i+1]-1]==_i) *b*=sLD->d[sLD->i[_i+1]-1]; *c=*b; }
  //Calculate c=LI.b
  for (_i=0;_i<sLD->ni;_i++)
    {
    if (!(_j=sLD->i[_i+1]-sLD->i[_i])) continue;
    while (--_j)
      {
      *c+=*sLD->d*b[*sLD->j], sLD->d++, sLD->j++;
      }
    //Process the last element especially as it can be the diagonal
    if (*sLD->j!=_i) { *c+=*sLD->d*b[*sLD->j]; }
    sLD->d++, sLD->j++, c++;
    }
  sLD->d-=sLD->nnz, sLD->d-=sLD->nnz; 
  }
}

//This function calculates c=L.LT.b product
//NOTE. b vector is destroed during calculations
void multiple_sLLT_dvector(double *c,t_smatrix *sL,double *b)
{
register unsigned int _i,_j;
set_vect(sL->ni,c,0.);
if (sL->nnz)
  {
  //Calculate c=LT.b
  _i=sL->ni, sL->j+=sL->nnz, sL->d+=sL->nnz, b+=sL->ni;
  while (_i--)
    {
    b--;
    if (!(_j=sL->i[_i+1]-sL->i[_i])) { *b=0.; continue; }
    while (_j--)
      {
      sL->d--, sL->j--, c[*sL->j]+=*sL->d**b;
      }
    }
  //Reconfigure vectors
  _i=sL->ni, b+=sL->ni, c+=sL->ni; while (_i--) { b--, c--, *b=*c, *c=0.; }
  //Calculate c=L.b
  for (_i=0;_i<sL->ni;_i++)
    {
    if (!(_j=sL->i[_i+1]-sL->i[_i])) continue;
    while (_j--)
      {
      *c+=*sL->d*b[*sL->j], sL->d++, sL->j++;
      }
    c++;
    }
  sL->d-=sL->nnz, sL->d-=sL->nnz; 
  }
}


//These couple functions do linear solvation of diagonal sparse problem

//This function solves row-wise upper triangle sparse problem with unit main diagonal
//NOTE. Resulting vector is replaced input vector
void solve_sIU(t_smatrix *sIU,double *y)
{
register unsigned int _i,_j;
register double _d;

_i=sIU->ni;
_j=sIU->nnz;
while (_i--)
  {
  _d=y[_i];
  while (--_j>sIU->i[_i]) _d-=sIU->d[_j]*y[sIU->j[_j]];
  y[_i]=_d;
  _j--;
  }
}

//This function solves row-wise lower triangle sparse problem with unit main diagonal
//NOTE. Resulting vector is replaced input vector
void solve_sIUT(t_smatrix *sIL,double *y)
{
register unsigned int _i,_j;
register double _d;

_i=0x0;
_j=0x0;
while (++_i<sIL->ni)
  {
  _d=y[_i];
  while (++_j<sIL->i[_i]) y[sIL->j[_j]]-=sIL->d[_j]*_d;
  _j++;
  }
}

//This block of functions do sparse matrix multiplication


//This function builds portrain for C[s,s]=AT[s,x].A[x,s] : _tr[s], _ts[s]
//Note. Matrix should be ordered
t_smatrix *build_sparse_AAT_portrain(unsigned char*_tr,unsigned int*_ts,t_smatrix *A)
{
register unsigned int _i,_j,_k,_l;
t_smatrix *C=0x0;
unsigned int mem=0x0, *vp,count;
//Alloc some memory
if ( (!(C=(t_smatrix*)alloc_smatrix(TRUE,A->ni,A->ni,0x0)))                            ||
     (!(A->j=(unsigned int*)malloc(sizeof(unsigned int)*(mem=SMATRIX_REALLOC_QUANT)))) )
  {
  ERROR_EXIT: free_smatrix(C);
  return FALSE;
  }
//Do row by row matrix product generation. Use autocompression into L as AAT is symmetric
count=0x0;
for (_k=0;_k<A->ni;C->i[_k++]=C->nnz)
  {//Generate string
  for (_l=0;_l<=_k;_l++)
    {
    _i=A->i[_k-0x1];
    _j=A->i[_l-0x1];
    while ( (_i<A->i[_k])&&(_j<A->i[_l]) )
           if (A->j[_j]==A->j[_i])
             {
             _ts[_l-0x1]=TRUE;
             _tr[count++]=_l-0x1;
             break;
             }
      else if (A->j[_j]<A->j[_i]) _j++;
           else                   _i++;
    }
  if (mem<C->nnz+count)
    {//Realloc memory for matrix
    if (!(vp=(unsigned int*)realloc(C->j,sizeof(unsigned int)*(mem+=SMATRIX_REALLOC_QUANT)))) goto ERROR_EXIT;
    else C->j=vp;
    }
  //Copy string into smatrix
  for (_i=0;_i<count;_i++)
    C->j[C->nnz++]=_tr[_i];
  }
return C;
}

//This function builds portrain for C[r,s]=A[r,x].B[x,s] : _tr[s], _ts[s]
t_smatrix *build_sparse_AB_portrain(unsigned char *_tr,unsigned int *_ts,t_smatrix *A,t_smatrix *B)
{
register unsigned int _i,_j,_k,_l;
t_smatrix *C=0x0;
unsigned int mem=0x0;
unsigned int *vp=0x0;
//Alloc some memory
if ( (!(C=(t_smatrix*)alloc_smatrix(TRUE,A->ni,B->nj,0x0)))||(!(C->j=(unsigned int*)malloc(sizeof(unsigned int)*(mem=SMATRIX_REALLOC_QUANT)))) )
  {
  ERROR_EXIT: free_smatrix(C); C=0x0;
  return FALSE;
  }
//Init tr
_j=C->nj;
while(_j--) _tr[_j]=FALSE;
//Do row by row C portrain evaluation
C->i[C->nnz=0x0]=0x0;
for (_i=0x0;_i<A->ni;C->i[++_i]=C->nnz)
  {
  for (_l=0x0,_k=A->i[_i];_k<A->i[_i+0x1];_k++)
    for (_j=B->i[A->j[_k]];_j<B->i[A->j[_k]+0x1];_j++)
      if (!(_tr[B->j[_j]]))
        _tr[_ts[_l++]=B->j[_j]]=TRUE;
  //Merge data
  if ((C->nnz+_l)>mem)
    {
    if (!(vp=(unsigned int*)realloc(C->j,sizeof(unsigned int)*(mem=C->nnz+SMATRIX_REALLOC_QUANT)))) { free(vp); goto ERROR_EXIT; }
    else C->j=vp;
    }
  //Copy string into smatrix and reset flags
  if (!(u_qsort(_l,_ts))) { free_smatrix(C); C=0x0; return FALSE; }
  while (_l--) _tr[C->j[C->nnz++]=_ts[_l]]=FALSE;
  }
//Sync memory
if (!(vp=(unsigned int*)realloc(C->j,sizeof(unsigned int)*C->nnz))) { free_smatrix(C); goto ERROR_EXIT; }
else C->j=vp;
if (!(C->d=(double*)malloc(sizeof(double)*C->nnz)))                 { free_smatrix(C); goto ERROR_EXIT; }
return C; //Fast, but a lot of reallocs required
}


//This function calculate C[r,s]=A[r,x].B[x,s] : _td[s]
void multiple_sparse_AB(double *_td,t_smatrix *A,t_smatrix *B,t_smatrix *C)
{
register unsigned int _i,_j,_k;
register double _d;

//Do row by row C numerical evaluation
//Init x
_i=C->nj;
while (_i--) _td[_i]=0.00;
for (_i=0x0;_i<A->ni;)
  {//Calc A(i,?).B(A(i,?),?)
  for (_j=A->i[_i++];_j<A->i[_i];_j++)
    for (_d=A->d[_j],_k=B->i[A->j[_j]];_k<B->i[A->j[_j]+0x1];_k++)
      _td[B->j[_k]]+=_d*B->d[_k];
  //Save data in C matrix
  for (_k=C->i[_i-0x1];_k<C->i[_i];_k++)
    {
    C->d[_k]=_td[C->j[_k]];
    _td[C->j[_k]]=0.00; //reinit _td
    }
  }
}

//This function multiple two sparse matrices C=A.B
t_smatrix *multiple_origin_smatrix_origin_smatrix(t_smatrix *A,t_smatrix *B)
{
double *_t=0x0;
t_smatrix *C=0x0;
if (!(_t=(double*)malloc(sizeof(double)*B->nj))) return FALSE;
if ( (C=build_sparse_AB_portrain((unsigned char*)&((unsigned int*)_t)[B->nj],(unsigned int*)_t,A,B)))
  {
  printf("Multipling portrain buit, performing numerical evaluations...\n");
  multiple_sparse_AB(_t,A,B,C);
  }
free(_t);
return C;
}

//Unfortunately there is no efficient way to calculate AT.A product :(
//This function calculates A.AT product
t_smatrix *multiple_sparse_AAT(unsigned char *_tr,unsigned int *_ts,t_smatrix *A)
{
double *_t=0x0;
t_smatrix *C=0x0;
if (!(_t=(double*)malloc(sizeof(double*)*A->nj))) return FALSE;
//if ( (C=build_sparse_AAT_portrain((unsigned char*)&((unsigned int*)_t)[A->nj],(unsigned int*)_t,A))) multiple_sparse_AAT(_t,A,C);
free(_t);
return C;
}




//This block of functions performs Cholesky decomposition of sparse matrix A.

//This function construct L portrain of cholesky decomposed matrix on the base of elimination tree of A
//Jeroen van Grondelle "Symbolic Sparse Cholesky Factorisation Using Elimination Trees"
//Note. it expects sA->j to be row-wise sorted so sA->j[_i]<sA->j[_i+1] if i and i+1 are in th same row
t_smatrix *build_cholesky_portain(t_smatrix *sA)
{
register unsigned int _i, _j, _k, _l, _p, _t;
unsigned int *pr=0x0, *vr=0x0; //parent and virtual columns
t_smatrix *sL=0x0;
unsigned int *vp;

//Stage 0.0. Adequacy test
if (sA->ni!=sA->nj) { ylib_errno=YERROR_IMPOSSIBLE; return FALSE; }
//Stage 0.1. Memory allocations
//Create pr
if (!(pr=(unsigned int*)malloc(sizeof(unsigned int)*2*sA->ni)))      { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return FALSE; }
else vr=pr+sA->ni;
//Create L
if (!(sL=(t_smatrix*)calloc(sizeof(t_smatrix),0x1)))                 { LABEL_MEMORY_ERROR_1: free(pr);    pr=0x0;    goto LABEL_MEMORY_ERROR_0; }
else { sL->ni=sL->nj=sA->ni; }
if (!(sL->i=(unsigned int*)malloc(sizeof(unsigned int)*(sA->ni+1)))) { LABEL_MEMORY_ERROR_2: free(sL);    sL=0x0;    goto LABEL_MEMORY_ERROR_1; }
else sL->i[0]=0;
if (!(sL->j=(unsigned int*)malloc(sizeof(unsigned int)*sA->nnz)))    { LABEL_MEMORY_ERROR_3: free(sL->i); sL->i=0x0; goto LABEL_MEMORY_ERROR_2; }
//Stage I. Build elimination tree of sA
//Do Lui algorithm of elimination tree construction with path compression
for (_i=0;_i<sA->ni;_i++)
  {
  pr[_i]=(unsigned int)-1, vr[_i]=(unsigned int)-1;
  for (_j=sA->i[_i];_j<sA->i[_i+1]-1;_j++)
    {
    _k=sA->j[_j]; while (vr[_k]<_i) { _t=vr[_k], vr[_k]=_i, _k=_t; }
    if (vr[_k]==(unsigned int)-1) { pr[_k]=_i, vr[_k]=_i; }
    }
  }
//Stage II. Map each row using elimination tree
for (sL->nnz=_i=0; _i<sA->ni; sL->i[++_i]=sL->nnz)
  {
  for (_k=sL->i[_i], _j=sA->i[_i]; _j!=sA->i[_i+1]-1; _j++) 
    {//Merge sL->j[] and a new _t, yet keeping sL->j sorted (sA->j is expected o be already sorted)
    _t=sA->j[_j]; while ( (_k!=sL->nnz)&&(sL->j[_k]<_t) ) _k++;
    if ( (_k==sL->nnz)||(sL->j[_k]!=_t) ) 
      {//Insert a new _t at position _k 
      if ( (sL->nnz>=sA->nnz)&&(!((sL->nnz-sA->nnz)%SMATRIX_REALLOC_QUANT)) )
        {//Memory management insertion
        if (!(vp=(unsigned int*)realloc(sL->j,sizeof(unsigned int)*(sL->nnz+SMATRIX_REALLOC_QUANT))))
          { free(sL->j); sL->j=0x0; goto LABEL_MEMORY_ERROR_3; }
        else sL->j=vp;
        }
      _l=++sL->nnz; while (--_l!=_k) sL->j[_l]=sL->j[_l-1];
      sL->j[_k]=_t, _p=++_k; //We use a new variable _p because parent's chain can be longer than sA->j[j+1] so sL->j would become unsorted
      while ((_t=pr[_t])<_i)
        {//Merge the new _t's parents chain
        while ( (_p!=sL->nnz)&&(sL->j[_p]<_t) ) _p++;
        if ( (_p==sL->nnz)||(sL->j[_p]!=_t) ) 
          {//Insert a new _t at position _p 
          if ( (sL->nnz>=sA->nnz)&&(!((sL->nnz-sA->nnz)%SMATRIX_REALLOC_QUANT)) )
            {//Memory management insertion
            if (!(vp=(unsigned int*)realloc(sL->j,sizeof(unsigned int)*(sL->nnz+SMATRIX_REALLOC_QUANT))))
              { free(sL->j); sL->j=0x0; goto LABEL_MEMORY_ERROR_3; }
            else sL->j=vp;
            }
          _l=++sL->nnz; while (--_l!=_p) sL->j[_l]=sL->j[_l-1];
          sL->j[_p]=_t;
          }
        else break; //Stop merging if the element is already in sL->j
        } 
      }
    } 
  //Add diagonal optionally
  if ( (!sL->nnz)||(sL->j[sL->nnz-1]!=_i) )
    {
    if ( (sL->nnz>=sA->nnz)&&(!((sL->nnz-sA->nnz)%SMATRIX_REALLOC_QUANT)) )
      {//Memory management insertion
      if (!(vp=(unsigned int*)realloc(sL->j,sizeof(unsigned int)*(sL->nnz+SMATRIX_REALLOC_QUANT))))
        { free(sL->j); sL->j=0x0; goto LABEL_MEMORY_ERROR_3; }
      else sL->j=vp;
      }
    sL->j[sL->nnz++]=_i;
    }
  }

//Stage III. Sync memory and exit
if (!(vp=(unsigned int*)realloc(sL->j,sizeof(unsigned int)*sL->nnz)))
  { free(sL->j); sL->j=0x0; goto LABEL_MEMORY_ERROR_3; }
else sL->j=vp;
if (!(sL->d=(double*)malloc(sizeof(double)*sL->nnz)))
  { free(sL->j); sL->j=0x0; goto LABEL_MEMORY_ERROR_3; }
free(pr); pr=0x0;
return sL;
}

//Calc sparse cholesky decomposition for portrain
//This function seems build L matrix A=L.LT, might fail due to negative eigenvalues but more precisse then LI.D.LIT decomposition.
char calc_sparse_cholesky_L(t_smatrix *sA,t_smatrix *sL)
{
register double _d, __d;
register unsigned int row_i, row_j, _i, _j, _k;

for (_i=0;_i<sL->ni;_i++)
  {
  //Do under-diagonal elements Lij=(Aij-SUMM[k=1...j-1](Lik*Ljk*Dk))/Dj and prepare diagonal element Di=Aii-SUMM[k=1...j-1](Lik*Lik*Dk)
  for (__d=0., _k=sA->i[_i], _j=sL->i[_i];_j<sL->i[_i+1]-1;_j++)
    {
    //Calc summ    
    if (sA->j[_k]==sL->j[_j]) { _d=sA->d[_k]; _k++; }
    else                        _d=0.;
    row_j=sL->i[sL->j[_j]], row_i=sL->i[_i];
    while ( (row_i<_j)&&(row_j<sL->i[sL->j[_j]+1]-1) )
      {
           if (sL->j[row_j]< sL->j[row_i]) row_j++;
      else if (sL->j[row_j]==sL->j[row_i]) { _d-=sL->d[row_i]*sL->d[row_j], row_j++, row_i++; }
      else if (sL->j[row_i]< sL->j[   _j]) row_i++;
      else break;
      }
    sL->d[_j]=_d/sL->d[sL->i[sL->j[_j]+0x1]-0x1];
    __d-=sqrd(sL->d[_j]);
    }
  //Calc diagonal element
  while (sA->j[_k]!=_i) _k++;
  __d+=sA->d[_k];
  if (__d<=EPSILON) { ylib_errno=YERROR_LEGAL; return FALSE; } //the matrix is semidefined, fail, try next function instead
  sL->d[sL->i[_i+1]-1]=__d;
  sL->d[sL->i[_i+1]-1]=sqrt(sL->d[sL->i[_i+1]-1]);
  }
return TRUE;
}
//This function seems build LI.D matrix A=LI.D.LIT
void calc_sparse_cholesky_LD(t_smatrix *sA,t_smatrix *sLD)
{
register double _d, __d;
register unsigned int row_i, row_j, _i, _j, _k;

for (_i=0;_i<sLD->ni;_i++)
  {
  //Do under-diagonal elements Lij=(Aij-SUMM[k=1...j-1](Lik*Ljk*Dk))/Dj and prepare diagonal element Di=Aii-SUMM[k=1...j-1](Lik*Lik*Dk)
  for (__d=0., _k=sA->i[_i], _j=sLD->i[_i];_j<sLD->i[_i+1]-1;_j++)
    {
    //Calc summ    
    if (sA->j[_k]==sLD->j[_j]) { _d=sA->d[_k]; _k++; }
    else                        _d=0.;
    row_j=sLD->i[sLD->j[_j]], row_i=sLD->i[_i];
    while ( (row_i<sLD->i[_i+1]-1)&&(row_j<sLD->i[sLD->j[_j]+1]-1) )
      {
           if (sLD->j[row_j]< sLD->j[row_i]) row_j++;
      else if (sLD->j[row_j]==sLD->j[row_i]) { _d-=sLD->d[row_i]*sLD->d[row_j]*sLD->d[sLD->i[sLD->j[row_j]+1]-1], row_j++, row_i++; }
      else if (sLD->j[row_i]< sLD->j[   _j]) row_i++;
      else break;
      }
    sLD->d[_j]=_d/sLD->d[sLD->i[sLD->j[_j]+0x1]-0x1];
    __d-=_d*sLD->d[_j];
    }
  //Calc diagonal element
  while (sA->j[_k]!=_i) _k++;
  __d+=sA->d[_k];
  if (fabs(__d)<=EPSILON)
    {
    if (__d<0.) __d=-EPSILON;
    else        __d=+EPSILON;
    }
  sLD->d[sLD->i[_i+1]-1]=__d;
  }
}

//Cholesky decompositon root function.
//This function require row compessed sparse matrix. It switches on type flag 'l' or 'L' to A=L.LT and 'd' or 'D' to A=LI.D.LIT and NIMPLEMENTED othewise.
t_smatrix *sparse_cholesky_decomposition(char nftype,t_smatrix *sA)
{
t_smatrix *sL=0x0;

//Adequacy test
if (sA->ni!=sA->nj) { ylib_errno=YERROR_IMPOSSIBLE; return FALSE; }
if ( (nftype!='L')&&(nftype!='l')&&(nftype!='d')&&(nftype!='D') ) { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
//Build elimination tree of L first
if (!(sL=build_cholesky_portain(sA))) return FALSE;
//Perform numerical factorization
if ( (nftype=='l')||(nftype=='L') ) calc_sparse_cholesky_L(sA,sL);
else                                calc_sparse_cholesky_LD(sA,sL);
return sL;
}



/**************************************** B A C K S U B S T I T U T I O N    B L O C K ****************************************/

//This function do backsubstitution into L.LT matrix
void bsubstitute_sLLT(t_smatrix *sL,double *b)
{
register unsigned int _i,_j;
register double _d;
if (sL->nnz)  
  {
  //Do second pass L^-1.b
  for (_i=0;_i<sL->ni;_i++)
    {
    if (!(_j=sL->i[_i+1]-sL->i[_i])) continue;
    _d=0.;
    while(--_j)
      {
      _d-=*sL->d*b[*sL->j], sL->j++, sL->d++;
      }
    if (*sL->j!=_i) b[_i]=0.;
    else { b[_i]+=_d, b[_i]/=*sL->d; }
    sL->j++, sL->d++; 
    }
  sL->j--, sL->d--;  
  //Do first pass L^-T.b
  _i=sL->ni, _j=sL->nnz;
  while (_i--)
    {
    if ( (!(_j=sL->i[_i+1]-sL->i[_i]))||(*sL->j!=_i) ) { b[_i]=0.; continue; }
    _d=b[_i]/=(*sL->d), sL->j--, sL->d--;
    while (--_j) 
      {
      b[*sL->j]-=*sL->d*_d, sL->j--, sL->d--;
      }
    }
  sL->j++, sL->d++;
  }
else set_vect(sL->ni,b,0.);
}

//This function do backsubstitution into LDLT matrix
void bsubstitute_sLDLT(t_smatrix *sLD,double *b)
{
register unsigned int _i,_j;
register double _d;
if (sLD->nnz)  
  {
  //Do first pass LI^-1.b
  for (_i=0;_i<sLD->ni;_i++)
    {
    if (!(_j=sLD->i[_i+1]-sLD->i[_i])) continue;
    _d=0.;
    while(--_j)
      {
      _d-=*sLD->d*b[*sLD->j], sLD->j++, sLD->d++;
      }
    if (*sLD->j!=_i) b[_i]=_d-*sLD->d*b[*sLD->j];
    else b[_i]+=_d; 
    sLD->j++, sLD->d++; 
    }
  sLD->j-=sLD->nnz, sLD->d-=sLD->nnz;  
  //Do diagonal 
  _i=sLD->ni, b+=sLD->ni; while (_i--) { b--; if (sLD->j[sLD->i[_i+1]-1]==_i) *b/=sLD->d[sLD->i[_i+1]-1]; }
  //Do second pass LI^-T.b
  _i=sLD->ni, _j=sLD->nnz, sLD->j+=sLD->nnz-1, sLD->d+=sLD->nnz-1;
  while (_i--)
    {
    if (!(_j=sLD->i[_i+1]-sLD->i[_i])) { b[_i]=0.; continue; }
    if (*sLD->j==_i) sLD->j--, sLD->d--; 
    _d=b[_i];
    while (--_j) 
      {
      b[*sLD->j]-=*sLD->d*_d, sLD->j--, sLD->d--;
      }
    }
  sLD->j++, sLD->d++;
  }
else set_vect(sLD->ni,b,0.);
}




//This function solves sparse linear problem from given cholesky decomposition.
//Note. if ((y==x) && (x==b)) than no addition vector storage required, but vector b will be substituted with x data (example: scdsp(L,LT,b,b,temp))
/*inline char solve_cholesky_decomposed_sparse_problem(t_smatrix *L,t_smatrix *LT,t_vector *b,t_vector *x,t_vector *y)
{
return  ( ( (solve_Lsparse_triangle(L,b,y)))&&( (solve_Usparse_triangle(LT,y,x))) );
}

//This function solve sparse linear system with aid of cholesky decompositin technique. SOLVE[A.x=b] where A sparse matrix, x, b - dense vectors
//IF A is symmetrical matrix (symmetry=TRUE) do: solve[A.x=b] A=LL^T, Lx=L^-T.b, x=L^-1.y (where y=L^-T.b)
//ELSE do: AA=A^T.A, AA=LL^T, c=A^T.b, Lx=L^-T.c, x=L^-1.y (where y=L^-T.c)
char solve_sparse_linear_system_via_cholesky_decomposition(t_smatrix *A,t_vector *x,t_vector *b)
{
t_smatrix *AT,*C,*L,*LT;

if (!(AT=transpose_smatrix(A)))
  return FALSE;
if (!(C=multiple_origin_smatrix_origin_smatrix(AT,A)))
  {
  free_smatrix(AT);
  return FALSE;
  }
else multiple_origin_smatrix_dvector(b,AT,x);
free_smatrix(AT);
if (!(L=sparse_cholesky_decomposition(C)))
  {
  free_smatrix(C);
  return FALSE;
  }
free_smatrix(C);
if (!(LT=transpose_smatrix(L)))
  {
  free_smatrix(L);
  return FALSE;
  }
if (!(solve_cholesky_decomposed_sparse_problem(L,LT,x,x,b)))
  {
  free_smatrix(L);
  free_smatrix(LT);
  return FALSE;
  }
free_smatrix(L);
free_smatrix(LT);
return TRUE;
}
*/
/*
//This is a central function that solve sparse problems
char solve_sparse_problem(char apprimate,char symmetry,t_smatrix *A,t_vec *b,t_vec *x)
{
//approximate or full problem solution
if (approximate)
  {//LU on not to LU :)
  if (!symmetry)
    approximate_solve_sparse_linear_system_via_LU_decomposition(A,b,x);
  else
    {
    approximate_solve_sparse_linear_system_via_cholesky_decomposition(A,b,x);
    }
  }
else //Fuck! full sollution requested... Are you really shure to spend such a lot affords???
  {//LU on not to LU :)
  if ( (!symmetry)&&(is_LU_better(A)) )
    solve_sparse_linear_system_via_LU_decomposition(A,b,x);
  else
    {
    if (!symmetry)
      {
      C=multiple_sparse_transp_matrix_sparse_origin_matrix(A,AT);
      multiple_transp_sparse_matrix_dense_vector(C,b,x);
      free_smatrix(A);
      A=C;
      memexchange(b,x,sizeof(t_vector*));
	  }
    solve_sparse_linear_system_via_cholesky_decomposition(A,b,x);
    }
  }
}

//All the following functions works with column compressed matrix form if other is not stated.




//This function perform triangular sparse system solvation L.x=b [x=L^-1.b]
//Note triangular matrix should has its diagonal filled
char solve_ltriangle(t_smatrix *L,t_vector *x,t_vector *b);



*/

/********************************** S - M A T R I X    D - V E C T O R     M U L T I P L I C A T I O N S ********************************/

/*
//This function performs sparse matrix . dense vecrot multiplication: y=A.x
//This function calculates c=A.b
void multiple_origin_smatrix_dvector(double *c,t_smatrix *sA,double *b)
{
register unsigned int _i,_j;
set_vect(sL->ni,c,0.);
if (sL->nnz)
  {
  for (_i=0;_i<sL->ni;_i++, c++)
    {
    if (!(_j=sL->i[_i+1]-sL->i[_i])) continue;
    while (_j--)
      {
      *c+=*sL->d*b[*sL->j], sL->d++, sL->j++;
      }
    sL->d++, sL->j++;
    }
  sL->j-=sL->nnz, sL->d-=sL->nnz, c-=sL->ni;
  }
}
//This function performs transpose sparse matrix . dense vector multiplication: c=A^T.b
void multiple_transp_smatrix_dvector(double *c,t_smatrix *sA,double *b)
{
register unsigned int _i,_j;
set_vect(sA->ni,c,0.);
if (sA->nnz)
  {
  sA->j+=sA->nnz, sA->d+=sL->nnz;
  _i=sA->ni;
  while (_i--)
    {
    if (!(_j=sA->i[_i+1]-sA->i[_i])) continue;
    //process the first element especially as it can be diagonal
    sA->j--, sA->d--; 
    if (*sA->j==_i) c[_i]+=*sA->d*b[_i];
    else { c[*sA->j]+=*sA->d*b[_i]; }
    while (--_j)
      {
      sA->d--, sA->j--, c[*sA->j]+=*sA->d*b[_i];
      }
    }
 }
}
*/




