#include "y_pls.h"

extern int ylib_errno;

//This function performs scaled compressed Q.L matrix decomposition: A=L.Q=L.(I+Z^T.T^T.Z)
void calc_scompressed_QL_decomposition_dmatrix(unsigned int ni,unsigned int nj,unsigned int nr,double **A,double *diag,double **T)
{
register unsigned int _i, _j;
register double _d, *dp, *dp2, scale;
unsigned int i;
for (i=0; i!=nr; i++)
  {
  dp=&A[i][i], _d=0., _j=nj-i; while (_j--) { _d+=fabs(*dp), dp++; } //Get scale
  if (_d>TINY)
    {
    scale=1./_d, diag[i]=_d, _d=0., dp=&A[i][i], _j=nj-i; while (_j--) { *dp*=scale, _d+=*dp**dp, dp++; } //Get norm
    _d=sqrt(_d); 
    if (A[i][i]>0.) { scale=-1./(_d*(_d-A[i][i])), A[i][i]-=_d, diag[i]*=+_d; }
    else            { scale=-1./(_d*(_d+A[i][i])), A[i][i]+=_d, diag[i]*=-_d; }
    for (_i=i+1;_i<ni;_i++) 
      {
      dp=&A[i][i], dp2=&A[_i][i], _d=0.,  _j=nj-i; while (_j--) { _d+=*dp**dp2, dp++, dp2++; }
      dp=&A[i][i], dp2=&A[_i][i], _d*=scale, _j=nj-i; while (_j--) { *dp2+=*dp*_d, dp++, dp2++; }
      }
    //Store transformations in T[i]=-2*T-1.Y-1.v
    for (T[i][i]=scale, _i=0;_i<i;_i++)
      {
      _d=0., dp=&A[i][i], dp2=&A[_i][i], _j=nj-i; while (_j--) { _d+=*dp**dp2, dp++, dp2++; }  
      _d*=scale, T[_i][i]=T[_i][_i]*_d, _j=_i; while (_j--) T[_j][i]+=T[_j][_i]*_d;
      }
    }
  else { diag[i]=0., T[i][i]=0., _i=i; while (_i--) T[_i][i]=0.; }
  }
}

//This function calculates p=Q.y from previous decomposition p=Q.y=Z.T.Z)^T.y=(I+Z^T.T^T.Z).y
void calc_compressed_QL_vector_product(unsigned int ni,unsigned int nj,unsigned int nr,double **L,double **T,double *y,double *p,double *t)
{
unsigned int _i, _j;
double _d, *dp, *dp2;
//Calc p1=Z.y
_i=nr; while (_i--) { dp=&L[_i][_i], dp2=&y[_i], _d=0., _j=nj-_i; while (_j--) { _d+=*dp**dp2, dp++, dp2++; } p[_i]=_d; }
//Calc p2=T^T.p1
dp=t, dp2=T[0], _d=*p, _j=nr; while (_j--) { *dp=*dp2*_d, dp++, dp2++; } 
for (_i=1;_i<nr;_i++) { _d=p[_i], dp=&t[_i], dp2=&T[_i][_i], _j=nr-_i; while (_j--) { *dp+=*dp2*_d, dp++, dp2++; } }
//Calc p3=Z^T.p2 
dp=p, dp2=L[0], _d=*t, _j=nj; while (_j--) { *dp=*dp2*_d, dp++, dp2++; } 
for (_i=1;_i<nr;_i++) { _d=t[_i], dp=&p[_i], dp2=&L[_i][_i], _j=nj-_i; while (_j--) { *dp+=*dp2*_d, dp++, dp2++; } }
//Calc p=y+p3
_i=ni; while (_i--) { *p+=*y, p++, y++; }
}

//This function calculates p=Q.y from previous decomposition p=Q.y=(I+Z^T.T^T.Z)^T.y=(I+Z^T.T.Z).y
void calc_compressed_QLT_vector_product(unsigned int ni,unsigned int nj,unsigned int nr,double **L,double **T,double *y,double *p,double *t)
{
unsigned int _i, _j;
double _d, *dp, *dp2;
//Calc p1=Z.y
_i=nr; while (_i--) { dp=&L[_i][_i], dp2=&y[_i], _d=0., _j=nj-_i; while (_j--) { _d+=*dp**dp2, dp++, dp2++; } p[_i]=_d; }
//Calc p2=T.p1
_i=nr; while (_i--) { _d=0., dp=&p[_i], dp2=&T[_i][_i], _j=nr-_i; while (_j--) { _d+=*dp**dp2, dp++, dp2++; } t[_i]=_d; }
//Calc p3=Z^T.p2 && p=y+p3
dp=p, dp2=L[0], _d=*t, _j=nj; while (_j--) { *dp=*dp2*_d, dp++, dp2++; } 
for (_i=1;_i<nr;_i++) { _d=t[_i], dp=&p[_i], dp2=&L[_i][_i], _j=nj-_i; while (_j--) { *dp+=*dp2*_d, dp++, dp2++; } }
//Calc p=y+p3
_i=ni; while (_i--) { *p+=*y, p++, y++; }
}





//This function calculates PLS model serially. flag_T determines A or A^T is submitted.
char calc_serial_pls_model(unsigned int power,unsigned int ni,unsigned int nj,double **A,char flag_T,double *x,double *y)
{
void *vp;
double **U, **V, *sdiag, *diag, **L, **T, **K, *p, *b, *_t;
unsigned int _i, _j, nr, nk;
double _d, yy;

//Switch startegy and alloc memory accvordingly
if (flag_T) 
  {
  if (ni>nj) { LABEL_NIMPLEMENT: ylib_errno=YERROR_NIMPLEMENTED; return FALSE; } 
  else 
    {
    nr=ni, nk=nj;
    if (5*ni<3*nj)
      {
      if (!(vp=malloc( sizeof(double)*nr*2                      +  //alloc diag and sdiag
//                       sizeof(double)*nk                        +  //alloc b
                       sizeof(double)*nk*2                      +  //alloc p and _t
                      (sizeof(double*)+sizeof(double)*nr)*nr    +  //alloc K
                      (sizeof(double*)+sizeof(double)*nr)*nr*2  +  //alloc T and L
                      (sizeof(double*)+sizeof(double)*nr)*nr*2 ))) //alloc U and V
        { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
      diag=(double*)vp, sdiag=diag+nr, p=sdiag+nr, _t=p+nk;
//      diag=(double*)vp, sdiag=diag+nr, b=sdiag+nr, p=b+nk, _t=p+nk;
      K=(void*)_t+sizeof(double)*nk;                        for (K[0]=(void*)K+sizeof(double*)*nr, _i=1;_i<nr;_i++) K[_i]=K[_i-1]+nr; //using the fact that nr>=power   
      L=(void*)K+(sizeof(double*)+sizeof(double)*nr)*nr;    for (L[0]=(void*)L+sizeof(double*)*nr, _i=1;_i<nr;_i++) L[_i]=L[_i-1]+nr;
      U=(void*)L+(sizeof(double*)+sizeof(double)*nr)*nr;    for (U[0]=(void*)U+sizeof(double*)*nr, _i=1;_i<nr;_i++) U[_i]=U[_i-1]+nr;
      V=(void*)U+(sizeof(double*)+sizeof(double)*nr)*nr;    for (V[0]=(void*)V+sizeof(double*)*nr, _j=1;_j<nr;_j++) V[_j]=V[_j-1]+nr;
      //Validate decomposition part #1
//      _i=ni; while (_i--) x[_i]=yrnd();
//      multiple_transp_matrix_origin_vector(ni,nj,A,x,b);
      //A=(A^T)^T=(Q.L)^T=Q^T.L^T=Q^T.U.S.V^T, L^T=U.S.V^T
      calc_scompressed_QL_decomposition_dmatrix(ni,nj,nr,A,diag,K);
      //svd
      _i=nr; while (_i--) { _j=nr; while (_i!=--_j) L[_j][_i]=0.; L[_i][_i]=diag[_i]; while (_j--) L[_j][_i]=A[_i][_j]; }
      memset(U[0],FALSE,sizeof(double)*nr*nr); _i=nr; while (_i--) U[_i][_i]=1.;
      memset(V[0],FALSE,sizeof(double)*nr*nr); _j=nr; while (_j--) V[_j][_j]=1.;
      _d=Hausholder_bidiagonalization(diag,sdiag,nr,nr,nr,L,U,V);
      diagonalize(nr,nr,nr,sdiag,diag,V,U,_d*SMALL2*SMALL2*SMALL2);
      //Validate decomposition part #2
//      multiple_transp_matrix_origin_vector(nr,nr,V,x,p);
//      _i=nr; while (_i--) p[_i]*=diag[_i];
//      multiple_origin_matrix_origin_vector(nr,nr,U,p,_t);
//      _j=nj; while (ni!=_j--) _t[_j]=0.;
//      calc_compressed_QLT_vector_product(ni,nj,nr,A,K,_t,p,x);
//      while (_j--) if (fabs(b[_j]-p[_j])>SMALL2) printf("SVD failed: vectors differs significantly (%d   A.x=%f   Q^T.U.S.V^T.x=%f)\n",_i,b[_j],p[_j]);
      //Calculate p=(Q^T.U)^T.y=U^T.Q.y 
      calc_compressed_QL_vector_product(ni,nj,nr,A,K,y,_t,p);
      multiple_transp_matrix_origin_vector(nr,nr,U,_t,p); 
      //Krylov
//      yy=sqrt(calc_vect_norm(ni,y));
//      _j=nr; while(_j--) K[0][_j]=diag[_j]*p[_j]/yy;
//      for (_i=1;_i<power;_i++) { _j=nr; while (_j--) K[_i][_j]=diag[_j]*K[_i-1][_j]; }
//      gram_schmidt_ortonormalization(power,nr,K,_t);
      //Final evaluation
//      multiple_origin_matrix_origin_vector(power,nr,K,p,_t);
//      multiple_transp_matrix_origin_vector(power,nr,K,_t,p);
//      _j=nr; while (_j--) p[_j]/=sqrt(diag[_j]);
//      multiple_origin_matrix_origin_vector(nr,nr,V,p,_t);
//      _j=nj; while (ni!=_j) _t[--_j]=0.; while (_j--) _t[_j]*=yy;
      //V-side product
      _j=nr; while (_j--) p[_j]/=diag[_j];
      multiple_origin_matrix_origin_vector(nr,nr,V,p,x);
      free(vp); vp=0x0;
      }  
    else goto LABEL_NIMPLEMENT;
    }
  }
else
  {//Normal A then
  if (ni>nj)
    {
    nr=nj, nk=ni;
    goto LABEL_NIMPLEMENT;
    }
  else
    {
    nr=ni, nk=nj;
    if (5*ni<3*nj) 
      {// svd via L.Q
      if (!(vp=malloc( sizeof(double)*ni                        +  //alloc b
                       sizeof(double)*nr*2                      +  //alloc diag and sdiag
                       sizeof(double)*nj*2                      +  //alloc p and _t
                      (sizeof(double*)+sizeof(double)*nr)*power +  //alloc K
                      (sizeof(double*)+sizeof(double)*nr)*nr*2  +  //alloc T and L
                      (sizeof(double*)+sizeof(double)*nr)*nr*2 ))) //alloc U and V
        goto LABEL_MEMORY_ERROR;
      b=(double*)vp, diag=b+ni, sdiag=diag+nr, p=sdiag+nr, _t=p+nj;
      K=(void*)_t+sizeof(double)*nj;                        for (K[0]=(void*)K+sizeof(double*)*power, _i=1;_i<power;_i++) K[_i]=K[_i-1]+nr;   
      T=(void*)K+(sizeof(double*)+sizeof(double)*nr)*power; for (T[0]=(void*)T+sizeof(double*)*nr, _i=1;_i<nr;_i++) T[_i]=T[_i-1]+nr;
      L=(void*)T+(sizeof(double*)+sizeof(double)*nr)*nr;    for (L[0]=(void*)L+sizeof(double*)*nr, _i=1;_i<nr;_i++) L[_i]=L[_i-1]+nr;
      U=(void*)L+(sizeof(double*)+sizeof(double)*nr)*nr;    for (U[0]=(void*)U+sizeof(double*)*nr, _i=1;_i<nr;_i++) U[_i]=U[_i-1]+nr;
      V=(void*)U+(sizeof(double*)+sizeof(double)*nr)*nr;    for (V[0]=(void*)V+sizeof(double*)*nr, _j=1;_j<nr;_j++) V[_j]=V[_j-1]+nr;
      //L.Q
      push_seed(2011,Y_MAGIC);                             // -.
      _j=nj; while (_j--) x[_j]=yrnd();                    //  |- check correctness
      multiple_origin_matrix_origin_vector(ni,nj,A,x,b);   // -'

      calc_scompressed_QL_decomposition_dmatrix(ni,nj,nr,A,diag,T);
      _i=nr; while (_i--) { _j=nr; while (_i!=--_j) L[_i][_j]=0.; L[_i][_i]=diag[_i]; while (_j--) L[_i][_j]=A[_i][_j]; }
      //svd
      memset(U[0],FALSE,sizeof(double)*nr*nr); _i=nr; while (_i--) U[_i][_i]=1.;
      memset(V[0],FALSE,sizeof(double)*nr*nr); _j=nr; while (_j--) V[_j][_j]=1.;
      _d=Hausholder_bidiagonalization(diag,sdiag,nr,nr,nr,L,U,V);
      diagonalize(nr,nr,nr,sdiag,diag,V,U,_d*SMALL2*SMALL2*SMALL2);

      calc_compressed_QL_vector_product(ni,nj,nr,A,T,x,p,_t);
      multiple_transp_matrix_origin_vector(nr,nr,V,p,_t);
      _i=nr; while (_i--) _t[_i]*=diag[_i];
      multiple_origin_matrix_origin_vector(ni,nr,U,_t,p);
      _i=ni; while (_i--) if (fabs(b[_i]-p[_i])>SMALL2) printf("SVD failed: vectors differs significantly (%d   A.x=%f   U.S.V.Q.x=%f)\n",_i,b[_i],p[_i]);
      //Krylov
//      multiple_transp_matrix_origin_vector(ni,nr,U,y,p);
//      yy=sqrt(calc_vect_norm(ni,y));
//      _j=nr; while(_j--) K[0][_j]=diag[_j]*p[_j]/yy;
//      for (_i=1;_i<power;_i++) { _j=nr; while (_j--) K[_i][_j]=diag[_j]*K[_i-1][_j]; }
//      gram_schmidt_ortonormalization(power,nr,K,_t);
      //Final evaluation
//      multiple_origin_matrix_origin_vector(power,nr,K,p,_t);
//      multiple_transp_matrix_origin_vector(power,nr,K,_t,p);
//      _j=nr; while (_j--) p[_j]/=sqrt(diag[_j]);
//      multiple_origin_matrix_origin_vector(nr,nr,V,p,_t);
//      _j=nj; while (ni!=_j) _t[--_j]=0.; while (_j--) _t[_j]*=yy;
//      calc_compressed_QLT_vector_product(ni,nj,nr,A,T,_t,x,p);

      //Pseudoinverse
      multiple_transp_matrix_origin_vector(ni,nr,U,y,p);
      _j=nr; while (_j--) p[_j]/=diag[_j];
      multiple_origin_matrix_origin_vector(nr,nr,V,p,_t);
      calc_compressed_QLT_vector_product(ni,nj,nr,A,T,_t,x,p);

      free(vp); vp=0x0;
      }
    else 
      {//Normal svd -- NOT DEBUGGED !~!!!!!____
      NORMAL_SVD: 
      if (!(vp=malloc( sizeof(double)*nr*2                      + //alloc diag and sdiag
                       sizeof(double)*nk*2                      + //alloc p and _t
                      (sizeof(double*)+sizeof(double)*nr)*power + //alloc K
                      (sizeof(double*)+sizeof(double)*ni)*ni+(sizeof(double*)+sizeof(double)*nj)*nj))) //alloc U and V
        goto LABEL_MEMORY_ERROR;
      diag=vp, sdiag=diag+nr, p=sdiag+nr, _t=p+nk;
      K=(void*)_t+sizeof(double)*nk;                        for (K[0]=(void*)K+sizeof(double*)*power, _i=1;_i<power;_i++) K[_i]=K[_i-1]+nr;   
      U=(void*)K+(sizeof(double*)+sizeof(double)*nr)*power; for (U[0]=(void*)U+sizeof(double*)*ni, _i=1;_i<ni;_i++) U[_i]=U[_i-1]+ni;
      V=(void*)U+(sizeof(double*)+sizeof(double)*ni)*ni;    for (V[0]=(void*)V+sizeof(double*)*nj, _j=1;_j<nj;_j++) V[_j]=V[_j-1]+nj;
      //svd
      memset(U[0],FALSE,sizeof(double)*ni*ni); _i=ni; while (_i--) U[_i][_i]=1.;
      memset(V[0],FALSE,sizeof(double)*nj*nj); _j=nj; while (_j--) V[_j][_j]=1.;
      _d=Hausholder_bidiagonalization(diag,sdiag,ni,nj,nr,A,U,V);
      diagonalize(ni,nj,nr,sdiag,diag,V,U,_d*SMALL2*SMALL2*SMALL2);
      //Krylov
      multiple_transp_matrix_origin_vector(ni,ni,U,y,p);
      yy=sqrt(calc_vect_norm(ni,y));
      _j=nr; while(_j--) K[0][_j]=diag[_j]*p[_j]/yy;
      for (_i=1;_i<power;_i++) { _j=nr; while (_j--) K[_i][_j]=diag[_j]*K[_i-1][_j]; }
      gram_schmidt_ortonormalization(power,nr,K,_t);
      //Final evaluation
      multiple_transp_matrix_origin_vector(power,nr,K,p,_t);
      multiple_transp_matrix_origin_vector(power,nr,K,_t,p);
      _j=nr; while (_j--) p[_j]/=sqrt(diag[_j]);
      multiple_origin_matrix_origin_vector(nj,nr,V,p,x);

      free(vp); vp=0x0;
      }
    }
  }
ylib_errno=YERR_OK;
return TRUE;
}

