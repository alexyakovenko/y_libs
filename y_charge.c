//This file contains charge calculations plugin
#include "y_charge.h"

#define MAX_DENSE_SIZE 0xF //it is not wisely to use dense matrix method for molecules bigger than 16 atoms

extern unsigned int ylib_errno;


/********************************************* I O N I Z E D    C D S ******************************************/
//These functions compute exact solution for trully uncharged and single-virtual-edged molecules; otherwise they compute a numerical solution using the mono-valent counter ion hardnesses as the initial guess that is further fitted with nsteps (or till tol convergence is reached) of Gauss-Newton algorithm. 
//
// Trivial part of cds is q=([H^1]-I)x and q+x=[H^-1]x. 
// To define ionized H on one virtual edge (i-j) of a counterion for given x we need denote [H^-1] as IH, then Serman-Morrison fomula stays:
// (H+d*u(X)u^T)^-1=IH-(IH.d*u.u^T.IH/(1+d*uT.IH.u)=IH-b*[c(X)c]/(1+b*c.u), where vector c==(uT.IH)^T==IH.u and (X) means cross-product
// For purposes of cds, u is constructed as [0...,-1,0...1] and b - is the unknown virtula hardness to resolve. 
// So q+x=IH0.x-b/(1+b*c.u)*[c(X)c].x; q[j]+x[j]=x[j]-b/(1+b*k)*Q`[j], where Q`[j] is i-th component of [c(X)c].x: Q`[j]=c[j]*(c^T.x) vector, where c[j]==1. and k=u^T.c=c[j]-c[i]=1-c[i]
// Resolving we have -q[i]/Q`[i]=K=b[ij]/(1+b[ij]*k) and b[ij]=K/(1-K*k); K=-q[i]/(c[j]*{c.x})=-q[i]/(1-x[n+1]+SUMM{c[1...n]*x[1...n]})
// P.S. Technically c=backsubstitute(H,u) where u:{u[k]=(k==j)?-1:0} thus K=-q[i]/Q`[i]=-q[i]/(c[i]*(c^T.x))=-q[i]/(1*(c{i-1}^T.x+1*x[i]))=-q[i]/(x[i]+c{i-1}^T.x)
// k=1-c[i] and finally b[ij]=K/(1-K*k). 
//
// To have the b's ultimately resolved we need amount of constraints (more or) equal to the amount of variables. 
// Thus we are adding more constraints to handle the typical situation of more virtual edges than virtual atoms
// In particular all b[i] of b=1/(b[i]+b[j]) at the same virtual atom are equal.   
// For the reason above we form only one edge of the virtual atom to the central atom in the case of NH3(+), CN2H4(+), COO(-), SO3(-), PO3(2-), SIO2(-) etc groups.
// The solution for b's is found with conjugate gradient first derivative method
// The derivative is calculated as: dq/db[j]=(dIH/db).x=-IH.(dH(b)/db).IH.x=-IH.u(X)u^T.IH.x=-c(X)c^T.x=-_d*c, where c=IH.u and _d={c^T.x}
// and Err=Summ[i](q[i]-q0[i])^2 so dErr/db[j]=SUMM[i](2*(q[i]-q0[i])*dq[i]/db[j])




/************************************************   D E N S E   C D S   ***********************************************/ 

//This function construct Lower triangle dense matrix for ionized cds
inline void _build_ionized_CDS_tdm(unsigned int natoms,unsigned int nvatoms,unsigned int nedges,register t_edge *edges,unsigned int nvedges,unsigned int (*vedges)[2],register double *h,register double **dL)
{
register unsigned int _i, _j;
register union{
              double *dp;
              double  d;
              } _d;
_i=natoms+nvatoms; while (_i--) { _d.dp=dL[_i]; _j=_i; while (_j--) *_d.dp++=0.; *_d.dp=1.; }
_i=nedges;
while (_i--)
  {
  _d.d=*h;
  if ( edges->vertice[0]>edges->vertice[1] ) dL[edges->vertice[0]][edges->vertice[1]]=-_d.d; else dL[edges->vertice[1]][edges->vertice[0]]=-_d.d;
  dL[edges->vertice[0]][edges->vertice[0]]+=_d.d, dL[edges->vertice[1]][edges->vertice[1]]+=_d.d;
  h++, edges++;
  }
_i=nvedges;
while (_i--)
  {
  _d.d=*h;
  dL[natoms+(*vedges)[1]][(*vedges)[0]]=-_d.d;
  dL[(*vedges)[0]][(*vedges)[0]]+=_d.d, dL[natoms+(*vedges)[1]][natoms+(*vedges)[1]]+=_d.d;
  h++, vedges++;
  }
}
//This function calculates ionized charges
//NB! It expects x and dL to be pre-build/pre-calculated before the call
inline void _calc_ionized_CDS_charges_tdm(register double *q,unsigned int natoms,register double *x,double **dL) 
{
register unsigned int _i; 
memcpy(q,x,sizeof(double)*natoms);
bsubstitute_tdLDLT(natoms,dL,q);  
_i=natoms, q+=_i, x+=_i; while (_i--) { --q, --x, *q-=*x; }
}
//This function calculates beta parameter of CDS using Sherman-Morrison formula.
//NB! It expects x and dL to be pre-build/pre-calculated before the call
inline double _calc_ionized_CDS_beta_tdm(unsigned int beta_i,double beta_Q,double beta_x,register double *u,unsigned int natoms,register double *x,double **dL) 
{
register unsigned int _i;
register double _d; 
memset(u,0,sizeof(double)*natoms); u[beta_i]=-1.;                     //  u:: [0,0,..., -1,0,...,0,+1]
bsubstitute_tdLDLT(natoms,dL,u);                                      //  c <- L^-1.u
_d=0., _i=natoms, u+=_i, x+=_i; while (_i--) { --u, --x, _d+=*u**x; } // _d <- {c.x}
_d=-beta_Q/(beta_x+_d);                                               // K=-q[i]/(c[j]*({c.x}+c[n+1]*x[n+1])), c[j]::c[n+1]==1,
_d/=(1.-_d*(1.-u[beta_i]));                                           // 1/bij=K/(1-K*k), k=u^T.c=c[j]-c[i]=1-ci ; c[j]==1
return _d;                                                            // bij=1/(1/bij)
}

//This function calculates errors of virtual charges
double _calc_ionized_CDS_tdm(unsigned int nbs,double *bs,double *q,unsigned int natoms,unsigned int nvatoms,int *vatoms,double *x,unsigned int nedges,t_edge *edges,unsigned int nvedges,unsigned int (*vedges)[2],double *h,double **dL)
{
register unsigned int _i;
register double err;
//Stage I. Build operator
_i=nbs; while (_i--) h[nedges+_i]=bs[_i];
_build_ionized_CDS_tdm(natoms,nvatoms,nedges,edges,nvedges,vedges,h,dL);
cholesky_decomposition_tdmatrix('L',natoms+nvatoms,dL);
//Stage II. Calc the error
_calc_ionized_CDS_charges_tdm(q,natoms+nvatoms,x,dL); 
_i=nvatoms, err=0.; while (_i--) err+=sqrd(q[natoms+_i]-(double)vatoms[_i]);                           //Calc the error err=SUMM((q[i]-q0[i])^2)
return err;
}
//This function is a minimizer-compatible wrapper for ionized charges error calculation routine
char _calc_ionized_CDS_wrapper_tdm(double *err,unsigned int nbs,double *bs,double *dbs,double **G,va_list stack)
{
unsigned int natoms, nvatoms, nedges, nvedges, (*vedges)[2];
int *vatoms;
double *q, *x, *h, **dL;
t_edge *edges;
va_list _stack;
//Stage 0. Unwrap stack
va_copy(_stack,stack);
      q=va_arg(_stack,double*); 
 natoms=va_arg(_stack,unsigned int);
nvatoms=va_arg(_stack,unsigned int);
 vatoms=va_arg(_stack,int*);
      x=va_arg(_stack,double*); 
 nedges=va_arg(_stack,unsigned int);
  edges=va_arg(_stack,t_edge*);
nvedges=va_arg(_stack,unsigned int); 
 vedges=va_arg(_stack,unsigned int(*)[2]);
      h=va_arg(_stack,double*);
     dL=va_arg(_stack,double**);
//Stage I. Calculate error
*err=_calc_ionized_CDS_tdm(nbs,bs,q,natoms,nvatoms,vatoms,x,nedges,edges,nvedges,vedges,h,dL);
//Stage II. End with the stack and exit
va_end(_stack);
return TRUE;
}
 
//This function calculates errors and derivatives of virtual charges
double _calc_ionized_dCDS_tdm(unsigned int nbs,double *bs,double *dbs,double *q,unsigned int natoms,unsigned int nvatoms,int *vatoms,double *x,unsigned int nedges,t_edge *edges,unsigned int nvedges,unsigned int (*vedges)[2],double *h,double **dL,double *tmp_dq)
{
register unsigned int _i, _j;
register double _d, err;
//Stage I. Build operator
_j=nbs; while (_j--) h[nedges+_j]=bs[_j];
_build_ionized_CDS_tdm(natoms,nvatoms,nedges,edges,nvedges,vedges,h,dL);
cholesky_decomposition_tdmatrix('L',natoms+nvatoms,dL);
//Stage II. Calc the error
_calc_ionized_CDS_charges_tdm(q,natoms+nvatoms,x,dL); 
_i=nvatoms, err=0.; while (_i--) { tmp_dq[_i]=q[natoms+_i]-(double)vatoms[_i], err+=sqrd(tmp_dq[_i]); }  //Calc the error Err=SUMM((q[i]-q0[i])^2)
//Stage III. Calc derivatives
_j=nbs;
while (_j--)
  {
  memset(q,0,sizeof(double)*(natoms+nvatoms)), q[vedges[_j][0]]=-1., q[natoms+vedges[_j][1]]=+1.; //Build 'u' vector
  bsubstitute_tdLDLT(natoms+nvatoms,dL,q);                                                        //Calc c vector c=IL.u
  _i=natoms+nvatoms, q+=_i, x+=_i, _d=0.; while (_i--) { --q, --x, _d+=*q**x; }                   //Calc _d scalar _d={c.x}
  _i=nvatoms, dbs[_j]=0.; while (_i--) dbs[_j]-=2.*tmp_dq[_i]*_d*q[natoms+_i];                    //Calc dErr/db[j]=SUMM(2*(q[i]-q0[i])*dq[i]/db[j]); dq[i]/db[j]=-_d*c
  }
return err;
}
//This function is a minimizer-compatible wrapper for ionized charges errors and derivatives calculation routine
char _calc_ionized_dCDS_wrapper_tdm(double *err,unsigned int nbs,double *bs,double *dbs,double **G,va_list stack)
{
unsigned int natoms, nvatoms, nedges, nvedges, (*vedges)[2];
int *vatoms;
double *q, *x, *h, **dL, *tmp_dq;
t_edge *edges;
va_list _stack;
//Stage 0. Unwrap stack
va_copy(_stack,stack);
      q=va_arg(_stack,double*); 
 natoms=va_arg(_stack,unsigned int);
nvatoms=va_arg(_stack,unsigned int);
 vatoms=va_arg(_stack,int*);
      x=va_arg(_stack,double*); 
 nedges=va_arg(_stack,unsigned int);
  edges=va_arg(_stack,t_edge*);
nvedges=va_arg(_stack,unsigned int); 
 vedges=va_arg(_stack,unsigned int(*)[2]);
      h=va_arg(_stack,double*);
     dL=va_arg(_stack,double**);
 tmp_dq=va_arg(_stack,double*);
//Stage I. Calculate gradients
*err=_calc_ionized_dCDS_tdm(nbs,bs,dbs,q,natoms,nvatoms,vatoms,x,nedges,edges,nvedges,vedges,h,dL,tmp_dq);
//Stage II. End with the stack and exit
va_end(_stack);
return TRUE;
} 

//This function calculates molecules charge. The inverted operator of internal electronic strucutre of a molecule is stored as a real dense matrix.
//NB! The mol need to possess correct natoms, nvatoms, vatoms[], engs[], ^charges[], nedges, edges[], nvedges, vedges[], ^hrds[] massives
char calc_ionized_CDS_cholesky_tdm(unsigned int nsteps,double tol,t_mol *mol)
{
register unsigned int _j;
double *(bs[2]), *(dbs[2]), *p, *tmp_dq, err; //the temp is a storage for dErr/dq data
t_dmatrix *dL;
//Stage I. Prepare memory
if ( (mol->charges)) { free(mol->charges), mol->charges=0x0; }
if (!(mol->charges=(double*)malloc(sizeof(double)*(mol->natoms+mol->nvatoms)))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }  
if ( (mol->cmtype))
  {
       if (mol->cmtype==-1) { mol->cmtype=FALSE; if ( (mol->C.dL)) { free(mol->C.dL),         mol->C.dL=0x0; } }
  else if (mol->cmtype==+1) { mol->cmtype=FALSE; if ( (mol->C.sL)) { free_smatrix(mol->C.sL), mol->C.sL=0x0; } }
  else { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
  }
if (!(dL=(t_dmatrix*)alloc_tdmatrix('L',mol->natoms))) goto LABEL_MEMORY_ERROR;
else dL->ni=dL->nj=mol->natoms;
//Stage II. Build kernel dL matrix (lower triangle only)
_build_ionized_CDS_tdm(mol->natoms,0,mol->nedges,mol->edges,0,0x0,mol->hrds,dL->d);
cholesky_decomposition_tdmatrix('L',mol->natoms,dL->d);
//Stage III. Switch strategy
//Stage III.1. Skip not ionized
if ( (mol->nvedges))
  {
  //Stage III.2. Solve for beta(s)
  //Stage III.2.a. Gather charged statistics
  tmp_dq=&mol->charges[mol->natoms]; memset(tmp_dq,0,sizeof(double)*mol->nvatoms);
  _j=mol->nvedges; while (_j--) tmp_dq[mol->vedges[_j][1]]+=1.;
  _j=mol->nvatoms; while (_j--) if (!(tmp_dq[_j])) { free(dL), dL=0x0; ylib_errno=YERROR_INTERNAL_CODE; return FALSE; } //A unconnected counterion detected!
  //Stage III.2.b. Resolve beta(s)
  _j=mol->nvedges;
  while (_j--)
    mol->hrds[mol->nedges+_j]=_calc_ionized_CDS_beta_tdm(mol->vedges[_j][0],(double)mol->vatoms[mol->vedges[_j][1]]/tmp_dq[mol->vedges[_j][1]],mol->engs[mol->natoms+mol->vedges[_j][1]],mol->charges,mol->natoms,mol->engs,dL->d);
  //Stage III.3. Allocate memory for kernel dL matrix (lower triangle only)
  free(dL), dL=0x0;
  if (!(dL=(t_dmatrix*)alloc_tdmatrix('L',mol->natoms+mol->nvatoms))) goto LABEL_MEMORY_ERROR;
  else dL->ni=dL->nj=mol->natoms+mol->nvatoms;
  //Stage III.4. Skip analytic case
  if (mol->nvedges!=1)
    {//Stage III.5. Solve standard case
    if (!(p=(double*)malloc(sizeof(double)*(5*mol->nvedges+mol->nvatoms)))) goto LABEL_MEMORY_ERROR;
    else
      {
      bs[0]=p+mol->nvedges, bs[1]=bs[0]+mol->nvedges, dbs[0]=bs[1]+mol->nvedges, dbs[1]=dbs[0]+mol->nvedges, tmp_dq=dbs[1]+mol->nvedges;
      memcpy(bs[0],&mol->hrds[mol->nedges],sizeof(double)*mol->nvedges);
      }
    //Run Polak-Ribiere optimization routine
    if ( (!(polak_ribiere_(&err,nsteps,tol,tol*SMALL2,1.,mol->nvedges,bs,dbs,0x0,p,
                           _calc_ionized_CDS_wrapper_tdm,
                           _calc_ionized_dCDS_wrapper_tdm,
                           line_search_square_fapproximation_,0x0,
                           mol->charges,mol->natoms,mol->nvatoms,mol->vatoms,mol->engs,mol->nedges,mol->edges,mol->nvedges,mol->vedges,mol->hrds,dL->d,tmp_dq)))&&
         (ylib_errno!=YERROR_NCONVERGED)&&(ylib_errno!=YERROR_SUSPICIOUS) )
      { free(p); p=0x0; return FALSE; }
    _j=mol->nvedges; while (_j--) mol->hrds[mol->nedges+_j]=bs[0][_j]; free(p); p=0x0; //Copy the solution (before freeing the memory)
    }
  //Finalize the run
  _build_ionized_CDS_tdm(mol->natoms,mol->nvatoms,mol->nedges,mol->edges,mol->nvedges,mol->vedges,mol->hrds,dL->d);
  cholesky_decomposition_tdmatrix('L',mol->natoms+mol->nvatoms,dL->d);
  }
_calc_ionized_CDS_charges_tdm(mol->charges,mol->natoms+mol->nvatoms,mol->engs,dL->d);
mol->cmtype=+1, mol->C.dL=dL;
return TRUE;
}

/**************************************   S P A R S E    C D S   **************************************************************/


//This function construct lower triangle sparse matrix for NON-ionized cds
inline void _build_non_ionized_CDS_spm(unsigned int natoms,unsigned int nedges,register t_edge *edges,register double *h,t_smatrix *sLA)
{
register unsigned int _i, _j;
register double _d;
//Stage I. Map spm memory.
//NB! Using the diagonal as size massive
memset(sLA->i,0,sizeof(unsigned int)*(natoms+1));
_j=nedges; while (_j--) if (edges[_j].vertice[0]<edges[_j].vertice[1]) sLA->i[edges[_j].vertice[1]+1]++; else sLA->i[edges[_j].vertice[0]+1]++;
for (_i=1; _i<natoms+1; _i++) { sLA->i[_i]+=1+sLA->i[_i-1], sLA->j[sLA->i[_i]-1]=sLA->i[_i-1], sLA->d[sLA->i[_i]-1]=1.; } 
//Stage II. Fill the matrix
_j=nedges; 
while (_j--) 
  {
  _d=h[_j];
  if (edges[_j].vertice[0]>edges[_j].vertice[1])
    { sLA->d[sLA->j[sLA->i[edges[_j].vertice[0]+1]-1]]=-_d, sLA->j[sLA->j[sLA->i[edges[_j].vertice[0]+1]-1]]=edges[_j].vertice[1], sLA->j[sLA->i[edges[_j].vertice[0]+1]-1]++; }
  else
    { sLA->d[sLA->j[sLA->i[edges[_j].vertice[1]+1]-1]]=-_d, sLA->j[sLA->j[sLA->i[edges[_j].vertice[1]+1]-1]]=edges[_j].vertice[0], sLA->j[sLA->i[edges[_j].vertice[1]+1]-1]++; }
  sLA->d[sLA->i[edges[_j].vertice[0]+1]-1]+=_d, sLA->d[sLA->i[edges[_j].vertice[1]+1]-1]+=_d;
  }
//Stage III. Sort the neighbors and restore the diagonal
_i=natoms; 
while (_i--) 
  {
  ud_isort(sLA->i[_i+1]-sLA->i[_i]-1,&sLA->j[sLA->i[_i]],&sLA->d[sLA->i[_i]]); //The amount of neighbors is expected to be small: 1-4 items
  sLA->j[sLA->i[_i+1]-1]=_i; 
  }
}
//This function is used for qsort
int _compare_vedges(const void *va,const void *vb)
{
if ( (*((unsigned int(*)[2])va))[1]==(*((unsigned int(*)[2])vb))[1] ) return (int)(*((unsigned int(*)[2])va))[0]-(int)(*((unsigned int(*)[2])vb))[0];
else                                                                  return (int)(*((unsigned int(*)[2])va))[1]-(int)(*((unsigned int(*)[2])vb))[1];
}
//This function construct lower triangle sparse matrix for ionized cds
inline void _build_ionized_CDS_spm(unsigned int natoms,unsigned int nvatoms,unsigned int nedges,register t_edge *edges,unsigned int nvedges,unsigned int (*vedges)[2],register double *h,t_smatrix *sLA)
{
register unsigned int _i, _j;
register double _d;
//Stage 0. Sort virtual edges to facilitate sLA (re-)generation
qsort((void*)vedges,(size_t)nvedges,sizeof(unsigned int)*2,_compare_vedges); 
//Stage I. Map spm memory.
//NB! Using the diagonal as size massive
memset(sLA->i,0,sizeof(unsigned int)*(natoms+nvatoms+1));
_j=nedges; while (_j--) if (edges[_j].vertice[0]<edges[_j].vertice[1]) sLA->i[edges[_j].vertice[1]+1]++; else sLA->i[edges[_j].vertice[0]+1]++;
_j=nvedges; while (_j--) sLA->i[natoms+vedges[_j][1]+1]++;
for (_i=1; _i<natoms+1; _i++) { sLA->i[_i]+=1+sLA->i[_i-1], sLA->j[sLA->i[_i]-1]=sLA->i[_i-1], sLA->d[sLA->i[_i]-1]=1.; } 
for (_i=natoms+1; _i<natoms+nvatoms+1; _i++) { sLA->i[_i]+=1+sLA->i[_i-1], sLA->d[sLA->i[_i]-1]=1.; } //vedges are sorted - no need to keep temporary counters for virtuals
//Stage II. Fill the real part of the matrix
_j=nedges; 
while (_j--) 
  {
  _d=h[_j];
  if (edges[_j].vertice[0]>edges[_j].vertice[1])
    { sLA->d[sLA->j[sLA->i[edges[_j].vertice[0]+1]-1]]=-_d, sLA->j[sLA->j[sLA->i[edges[_j].vertice[0]+1]-1]]=edges[_j].vertice[1], sLA->j[sLA->i[edges[_j].vertice[0]+1]-1]++; }
  else
    { sLA->d[sLA->j[sLA->i[edges[_j].vertice[1]+1]-1]]=-_d, sLA->j[sLA->j[sLA->i[edges[_j].vertice[1]+1]-1]]=edges[_j].vertice[0], sLA->j[sLA->i[edges[_j].vertice[1]+1]-1]++; }
  sLA->d[sLA->i[edges[_j].vertice[0]+1]-1]+=_d, sLA->d[sLA->i[edges[_j].vertice[1]+1]-1]+=_d;
  }
//Stage III. Sort the real neighbors and restore the diagonal 
_i=natoms; 
while (_i--) 
  {
  ud_isort(sLA->i[_i+1]-sLA->i[_i]-1,&sLA->j[sLA->i[_i]],&sLA->d[sLA->i[_i]]); //The amount of neighbors is expected to be small: 1-4 items
  sLA->j[sLA->i[_i+1]-1]=_i; 
  }
//Stage IV. Fill the virtual part of the matrix
for (_j=0; _j<nvedges; _j++) //Advantage of sorted vedges
  {
  _d=h[nedges+_j];
  sLA->d[sLA->i[natoms]+_j+vedges[_j][1]]=-_d, sLA->j[sLA->i[natoms]+_j+vedges[_j][1]]=vedges[_j][0];
  sLA->d[sLA->i[vedges[_j][0]+1]-1]+=_d, sLA->d[sLA->i[natoms+vedges[_j][1]+1]-1]+=_d;
  }
//The vedges are sorted, just restore the diag
_i=nvatoms; while (_i--) sLA->j[sLA->i[natoms+_i+1]-1]=natoms+_i;
}
//This function updates ionized sparse operator with new bs. Only diagonal and virtual edges should be updated.
//NB!! vedges should be pre-sorted.
inline void _update_ionized_CDS_spm(unsigned int natoms,unsigned int nvatoms,unsigned int nedges,t_edge *edges,unsigned int nvedges,unsigned int (*vedges)[2],double *h,t_smatrix *sLA)
{
register unsigned int _i, _j;
register double _d;
//Stage I. Unit diagonal
_i=natoms+nvatoms+1; while (--_i) sLA->d[sLA->i[_i]-1]=1.;
//Stage II. Fill the matrix diagonal
_j=nedges; 
while (_j--) 
  {
  _d=h[_j];
  sLA->d[sLA->i[edges[_j].vertice[0]+1]-1]+=_d, sLA->d[sLA->i[edges[_j].vertice[1]+1]-1]+=_d;
  }
for (_j=0; _j<nvedges; _j++) //Advantage of sorted vedges
  {
  _d=h[nedges+_j];
  sLA->d[sLA->i[natoms]+_j+vedges[_j][1]]=-_d, sLA->d[sLA->i[vedges[_j][0]+1]-1]+=_d, sLA->d[sLA->i[natoms+vedges[_j][1]+1]-1]+=_d;
  }
}

//This function calculates ionized charges
//NB! It expects x and dL to be pre-build/pre-calculated before the call
inline void _calc_ionized_CDS_charges_spm(register double *q,unsigned int natoms,register double *x,t_smatrix *sL) 
{
register unsigned int _i; 
memcpy(q,x,sizeof(double)*natoms);
bsubstitute_sLDLT(sL,q);
_i=natoms, q+=_i, x+=_i; while (_i--) { --q, --x, *q-=*x; }
}
//This function calculates beta parameter of CDS using Sherman-Morrison formula.
//NB! It expects x and dL to be pre-build/pre-calculated before the call
inline double _calc_ionized_CDS_beta_spm(unsigned int beta_i,double beta_Q,double beta_x,register double *u,unsigned int natoms,register double *x,t_smatrix *sL) 
{
register unsigned int _i;
register double _d; 
memset(u,0,sizeof(double)*natoms); u[beta_i]=-1.;                     //  u:: [0,0,..., -1,0,...,0,+1]
bsubstitute_sLDLT(sL,u);                                              //  c <- L^-1.u
_d=0., _i=natoms, u+=_i, x+=_i; while (_i--) { --u, --x, _d+=*u**x; } // _d <- {c.x}
_d=-beta_Q/(beta_x+_d);                                               // K=-q[i]/(c[j]*({c.x}+c[n+1]*x[n+1])), c[j]::c[n+1]==1,
_d/=(1.-_d*(1.-u[beta_i]));                                           // 1/bij=K/(1-K*k), k=u^T.c=c[j]-c[i]=1-ci ; c[j]==1
return _d;                                                            // bij=1/(1/bij)
}

//This function calculates errors and derivatives of virtual charges
double _calc_ionized_CDS_spm(unsigned int nbs,double *bs,double *q,unsigned int natoms,unsigned int nvatoms,int *vatoms,double *x,unsigned int nedges,t_edge *edges,unsigned int nvedges,unsigned int (*vedges)[2],double *h,t_smatrix *sLA,t_smatrix *sL)
{
register unsigned int _i;
register double err;
//Stage I. Build operator
_i=nbs; while (_i--) h[nedges+_i]=bs[_i];
_update_ionized_CDS_spm(natoms,nvatoms,nedges,edges,nvedges,vedges,h,sLA);
calc_sparse_cholesky_LD(sLA,sL);
//Stage II. Calc the error
_calc_ionized_CDS_charges_spm(q,natoms+nvatoms,x,sL); 
_i=nvatoms, err=0.; while (_i--) err+=sqrd(q[natoms+_i]-(double)vatoms[_i]);  //Calc the error Err=SUMM((q[i]-q0[i])^2)
return err;
}
//This function is a minimizer-compatible wrapper for ionized charges errors and derivatives calculation routine
char _calc_ionized_CDS_spm_wrapper(double *err,unsigned int nbs,double *bs,double *dbs,double **G,va_list stack)
{
unsigned int natoms, nvatoms, nedges, nvedges, (*vedges)[2];
int *vatoms;
double *q, *x, *h;
t_smatrix *sLA, *sL;
t_edge *edges;
va_list _stack;
//Stage 0. Unwrap stack
va_copy(_stack,stack);
      q=va_arg(_stack,double*); 
 natoms=va_arg(_stack,unsigned int);
nvatoms=va_arg(_stack,unsigned int);
 vatoms=va_arg(_stack,int*);
      x=va_arg(_stack,double*); 
 nedges=va_arg(_stack,unsigned int);
  edges=va_arg(_stack,t_edge*);
nvedges=va_arg(_stack,unsigned int); 
 vedges=va_arg(_stack,unsigned int(*)[2]);
      h=va_arg(_stack,double*);
    sLA=va_arg(_stack,t_smatrix*);
     sL=va_arg(_stack,t_smatrix*);
//Stage I. Calculate error
*err=_calc_ionized_CDS_spm(nbs,bs,q,natoms,nvatoms,vatoms,x,nedges,edges,nvedges,vedges,h,sLA,sL);
//Stage II. End with the stack and exit
va_end(_stack);
return TRUE;
} 
//This function calculates errors and derivatives of virtual charges
double _calc_ionized_dCDS_spm(unsigned int nbs,double *bs,double *dbs,double *q,unsigned int natoms,unsigned int nvatoms,int *vatoms,double *x,unsigned int nedges,t_edge *edges,unsigned int nvedges,unsigned int (*vedges)[2],double *h,t_smatrix *sLA,t_smatrix *sL,double *tmp_dq)
{
register unsigned int _i, _j;
register double _d, err;
//Stage I. Build operator
_j=nbs; while (_j--) h[nedges+_j]=bs[_j];
_update_ionized_CDS_spm(natoms,nvatoms,nedges,edges,nvedges,vedges,h,sLA);
calc_sparse_cholesky_LD(sLA,sL);
//Stage II. Calc the error
_calc_ionized_CDS_charges_spm(q,natoms+nvatoms,x,sL); 
_i=nvatoms, err=0.; while (_i--) { tmp_dq[_i]=q[natoms+_i]-(double)vatoms[_i], err+=sqrd(tmp_dq[_i]); }  //Calc the error Err=SUMM((q[i]-q0[i])^2)
//Stage III. Calc derivatives
_j=nvedges;
while (_j--)
  {
  memset(q,0,sizeof(double)*(natoms+nvatoms)), q[vedges[_j][0]]=-1., q[natoms+vedges[_j][1]]=+1.; //Build 'u' vector
  bsubstitute_sLDLT(sL,q);                                                                        //Calc c vector c=IL.u
  _i=natoms+nvatoms, q+=_i, x+=_i, _d=0.; while (_i--) { --q, --x, _d+=*q**x; }                   //Calc _d scalar _d={c.x}
  _i=nvatoms, dbs[_j]=0.; while (_i--) dbs[_j]-=2.*tmp_dq[_i]*_d*q[natoms+_i];                    //Calc dErr/db[j]=SUMM(2*(q[i]-q0[i])*dq[i]/db[j]); dq[i]/db[j]=-_d*c
  }
return err;
}
//This function is a minimizer-compatible wrapper for ionized charges errors and derivatives calculation routine
char _calc_ionized_dCDS_spm_wrapper(double *err,unsigned int nbs,double *bs,double *dbs,double **G,va_list stack)
{
unsigned int natoms, nvatoms, nedges, nvedges, (*vedges)[2];
int *vatoms;
double *q, *x, *h, *tmp_dq;
t_smatrix *sLA, *sL;
t_edge *edges;
va_list _stack;
//Stage 0. Unwrap stack
va_copy(_stack,stack);
      q=va_arg(_stack,double*); 
 natoms=va_arg(_stack,unsigned int);
nvatoms=va_arg(_stack,unsigned int);
 vatoms=va_arg(_stack,int*);
      x=va_arg(_stack,double*); 
 nedges=va_arg(_stack,unsigned int);
  edges=va_arg(_stack,t_edge*);
nvedges=va_arg(_stack,unsigned int); 
 vedges=va_arg(_stack,unsigned int(*)[2]);
      h=va_arg(_stack,double*);
    sLA=va_arg(_stack,t_smatrix*);
     sL=va_arg(_stack,t_smatrix*);
 tmp_dq=va_arg(_stack,double*);
//Stage I. Calculate gradients
*err=_calc_ionized_dCDS_spm(nbs,bs,dbs,q,natoms,nvatoms,vatoms,x,nedges,edges,nvedges,vedges,h,sLA,sL,tmp_dq);
//Stage II. End with the stack and exit
va_end(_stack);
return TRUE;
} 

//This function calculates molecules charge. The inverted operator of internal electronic strucutre of a molecule is stored as a real dense matrix.
//NB! The mol need to possess correct natoms, nvatoms, vatoms[], engs[], ^charges[], nedges, edges[], nvedges, vedges[], ^hrds[] massives
char calc_ionized_CDS_cholesky_spm(unsigned int nsteps,double tol,t_mol *mol)
{
register unsigned int _j;
double *(bs[2]), *(dbs[2]), *p, *tmp_dq, err; //the temp is a storage for dErr/dq data
t_smatrix *sLA, *sL;
//Stage I. Prepare memory
if ( (mol->charges)) { free(mol->charges), mol->charges=0x0; }
if (!(mol->charges=(double*)malloc(sizeof(double)*(mol->natoms+mol->nvatoms)))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }  
if ( (mol->cmtype))
  {
       if (mol->cmtype==-1) { mol->cmtype=FALSE; if ( (mol->C.dL)) { free(mol->C.dL),         mol->C.dL=0x0; } }
  else if (mol->cmtype==+1) { mol->cmtype=FALSE; if ( (mol->C.sL)) { free_smatrix(mol->C.sL), mol->C.sL=0x0; } }
  else { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
  }
if (!(sLA=alloc_smatrix(TRUE,mol->natoms,mol->natoms,mol->nedges+mol->natoms))) goto LABEL_MEMORY_ERROR;
else sLA->ni=sLA->nj=mol->natoms;
//Stage II.1. Build kernel sL matrix (lower triangle only)
_build_non_ionized_CDS_spm(mol->natoms,mol->nedges,mol->edges,mol->hrds,sLA);
//Stage II.2. Build elimination tree of sL and solve the default matrix
if (!(sL=build_cholesky_portain(sLA))) { free_smatrix(sLA); sLA=0x0; return FALSE; }
sL->ni=sL->nj=mol->natoms, sLA->nnz=sLA->i[mol->natoms], sL->nnz=sL->i[mol->natoms];
calc_sparse_cholesky_LD(sLA,sL);
free_smatrix(sLA), sLA=0x0;
//Stage III. Switch strategy
//Stage III.1. Skip not ionized
if ( (mol->nvedges))
  {
  //Stage III.2. Solve for beta(s)
  //Stage III.2.a. Gather statistics
  tmp_dq=&mol->charges[mol->natoms]; memset(tmp_dq,0,sizeof(double)*mol->nvatoms);
  _j=mol->nvedges; while (_j--) tmp_dq[mol->vedges[_j][1]]+=1.;
  _j=mol->nvatoms; while (_j--) if (!(tmp_dq[_j])) { free_smatrix(sL), sL=0x0; ylib_errno=YERROR_INTERNAL_CODE; return FALSE; } //A unconnected counterion detected!
  //Stage III.2. Resolve beta(s)
  _j=mol->nvedges;
  while (_j--)
    mol->hrds[mol->nedges+_j]=_calc_ionized_CDS_beta_spm(mol->vedges[_j][0],(double)mol->vatoms[mol->vedges[_j][1]]/tmp_dq[mol->vedges[_j][1]],mol->engs[mol->natoms+mol->vedges[_j][1]],mol->charges,mol->natoms,mol->engs,sL);
  //Stage III.3. Allocate memory for kernel sL matrix (lower triangle only)
  free_smatrix(sL), sL=0x0;
  if (!(sLA=(t_smatrix*)alloc_smatrix(TRUE,mol->natoms+mol->nvatoms,mol->natoms+mol->nvatoms,mol->nedges+mol->nvedges+mol->natoms+mol->nvatoms))) goto LABEL_MEMORY_ERROR;
  else sLA->ni=sLA->nj=mol->natoms+mol->nvatoms;
  _build_ionized_CDS_spm(mol->natoms,mol->nvatoms,mol->nedges,mol->edges,mol->nvedges,mol->vedges,mol->hrds,sLA);
  if (!(sL=build_cholesky_portain(sLA))) { free_smatrix(sLA); sLA=0x0; return FALSE; }
  if (mol->nvedges!=1)
    {//Stage III.5. Solve standard case
    if (!(p=(double*)malloc(sizeof(double)*(5*mol->nvedges+mol->nvatoms)))) goto LABEL_MEMORY_ERROR;
    else 
      {
      bs[0]=p+mol->nvedges, bs[1]=bs[0]+mol->nvedges, dbs[0]=bs[1]+mol->nvedges, dbs[1]=dbs[0]+mol->nvedges, tmp_dq=dbs[1]+mol->nvedges;
      memcpy(bs[0],&mol->hrds[mol->nedges],sizeof(double)*mol->nvedges);
      }
    //Run Polak-Ribiere optimization routine
    if ( (!(polak_ribiere_(&err,nsteps,tol,tol*SMALL2,1.,mol->nvedges,bs,dbs,0x0,p,
                           _calc_ionized_CDS_spm_wrapper,
                           _calc_ionized_dCDS_spm_wrapper,
                           line_search_square_fapproximation_,0x0,
                           mol->charges,mol->natoms,mol->nvatoms,mol->vatoms,mol->engs,mol->nedges,mol->edges,mol->nvedges,mol->vedges,mol->hrds,sLA,sL,tmp_dq)))&&
         (ylib_errno!=YERROR_NCONVERGED)&&(ylib_errno!=YERROR_SUSPICIOUS) )
      { free(p); p=0x0; return FALSE; }
    _j=mol->nvedges; while (_j--) mol->hrds[mol->nedges+_j]=bs[0][_j]; free(p); p=0x0; //Copy the solution (before freeing the memory)
    }
  //Finalize the run
  _update_ionized_CDS_spm(mol->natoms,mol->nvatoms,mol->nedges,mol->edges,mol->nvedges,mol->vedges,mol->hrds,sLA);
  calc_sparse_cholesky_LD(sLA,sL);
  free_smatrix(sLA), sLA=0x0;
  }
_calc_ionized_CDS_charges_spm(mol->charges,mol->natoms+mol->nvatoms,mol->engs,sL);
mol->cmtype=-1, mol->C.sL=sL;
return TRUE;
}

//printf("\n\nsLA matrix:\n");
//for (_i=0;_i<sLA->ni;_i++)
//  {
//  printf("%1d: ",_i);
//  for (_j=sLA->i[_i];_j<sLA->i[_i+1];_j++)
//    printf("%1d(%f) ",sLA->j[_j],sLA->d[_j]);
//  printf("\n");
//  }
//printf("\n\nsL matrix:\n");
//for (_i=0;_i<sL->ni;_i++)
//  {
//  printf("%1d: ",_i);
//  for (_j=sL->i[_i];_j<sL->i[_i+1];_j++)
//    printf("%1d(%f) ",sL->j[_j],sL->d[_j]);
//  printf("\n");
//  }
//for (_i=0; _i<sL->ni; _i++) printf("%1d -> %f\n",_i,mol->engs[_i]); printf("\n");
//for (_i=0; _i<sL->ni; _i++) printf("\n%1d -> %f, %f",_i,ksi[_i],mol->charges[_i]); printf("\n");
//for (_i=0; _i<sL->ni; _i++) printf("\n%1d -> %f",_i,mol->charges[_i]); printf("\n");  //Stage III.4. Skip analytic case


/*******************************   A   D E F A U L T   C D S   W R A P P E R   F U N C T I I O N   ********************************************/

//This function calculates molecules charge. The inverted operator of internal electronic strucutre of a molecule is stored as a real dense matrix.
char calc_Oliferenko_ionized_CDS_cholesky_default(unsigned int nsteps,double tol,t_mol *mol,t_top *top)
{
register  unsigned int _i, *up;
//Stage I. Preallocate memory
if ( (mol->engs)) { free(mol->engs); mol->engs=0x0; }
if (!(mol->engs=(double*)malloc(sizeof(double)*(mol->natoms+mol->nvatoms))))
  { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
if ( (mol->hrds)) { free(mol->hrds); mol->hrds=0x0; }
if (!(mol->hrds=(double*)malloc(sizeof(double)*(mol->nedges+mol->nvedges)))) goto LABEL_MEMORY_ERROR;
//Stage II. Setup defaults
if (!(up=(unsigned int*)calloc(mol->nvatoms,sizeof(unsigned int)))) goto LABEL_MEMORY_ERROR;
_i=mol->nvedges; while (_i--) up[mol->vedges[_i][1]]++;
_i=mol->nvatoms; while (_i--) if (!(up[_i])) { free(up); up=0x0; ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; } else mol->engs[mol->natoms+_i]=0.;
_i=mol->nvedges; while (_i--) mol->engs[mol->natoms+mol->vedges[_i][1]]+=top->ff_a[mol->ytypes[mol->vedges[_i][0]]].Eng;
_i=mol->nvatoms; while (_i--) { mol->engs[mol->natoms+_i]/=(double)up[_i]; if (mol->vatoms[_i]>0) mol->engs[mol->natoms+_i]/=top->ff_a[0].Eng; else mol->engs[mol->natoms+_i]*=top->ff_a[0].Eng; }
free(up); up=0x0;
_i=mol->natoms; while (_i--) mol->engs[_i]=top->ff_a[mol->ytypes[_i]].Eng;
_i=mol->nedges; while (_i--) mol->hrds[_i]=1./(top->ff_a[mol->ytypes[mol->edges[_i].vertice[0]]].Hrd+top->ff_a[mol->ytypes[mol->edges[_i].vertice[1]]].Hrd);
//Stage III. Calculate charges
if (mol->natoms+mol->nvatoms < COMPLEXITY_LIMIT_CDS_TDM) return (!(calc_ionized_CDS_cholesky_tdm(nsteps,tol,mol))) ? (unsigned int)-1 : mol->nvatoms;
else                                                     return (!(calc_ionized_CDS_cholesky_spm(nsteps,tol,mol))) ? (unsigned int)-1 : mol->nvatoms;
}



//printf("\n\nsLA matrix:\n");
//for (_i=0;_i<sLA->ni;_i++)
//  {
//  printf("%1d: ",_i);
//  for (_j=sLA->i[_i];_j<sLA->i[_i+1];_j++)
//    printf("%1d(%f) ",sLA->j[_j],sLA->d[_j]);
//  printf("\n");
//  }
//printf("\n\nsL matrix:\n");
//for (_i=0;_i<sL->ni;_i++)
//  {
//  printf("%1d: ",_i);
//  for (_j=sL->i[_i];_j<sL->i[_i+1];_j++)
//    printf("%1d(%f) ",sL->j[_j],sL->d[_j]);
//  printf("\n");
//  }
//for (_i=0; _i<sL->ni; _i++) printf("\n%1d -> %f, %f",_i,ksi[_i],mol->charges[_i]); printf("\n");
//for (_i=0; _i<sL->ni; _i++) printf("\n%1d -> %f",_i,mol->charges[_i]); printf("\n");


/*******************************************************   C H A R G E    G R O U P S   ************************************************************************/

//This function, IDEALLY, should edit charges and define charged groups tp have a molecules' electric field be represented by a set of small zero-charged fragments.
//Note. I am doubt if this is possible at all and particularly have no ideas of how to implement it.
//The current plan is to keep the true best charges on atoms but split the molecule so to reduce the 'artificial' dipole errors. 
//Stage I. Move all in the anchor into the charged group
//Stage II. Join H-X-C into one groups
//That is it so far...
unsigned int *split_charges(t_clist *neighbors,unsigned int *anchor_id,t_mol *mol,t_top *top)
{
register unsigned int _i, _j, _l;
unsigned int *chrg_grps=0x0;

//Stage 0. Allocate the memory
if (!(chrg_grps=(unsigned int*)malloc(sizeof(unsigned int)*mol->natoms))) { ylib_errno=YERROR_MEMORY; return FALSE; }

//Stage I.a. Move all multivalent into a separate charged group
_j=0, _i=mol->natoms; while (_i--) if (neighbors->list[_i].size!=1) chrg_grps[_i]=_j++;
//Stage I.b. Move the onevalent into corresponding multivalent charged group
_i=mol->natoms; while (_i--) if (neighbors->list[_i].size==1) { if (neighbors->list[_l=*neighbors->list[_i].list].size!=1) chrg_grps[_i]=chrg_grps[_l]; else if (_i<_l) chrg_grps[_i]=chrg_grps[_l]=_j; } 

//Stage II. Join H-O-X and H-S-X (where X!=H && valence[X]>2) into one group
_i=mol->natoms; 
while (_i--) 
  if ( (neighbors->list[_i].size==2)&&( (mol->a[_i]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[_i]==CHEM_ATOM_TYPE_SULFUR) )&&
     ( ( (mol->a[neighbors->list[_i].list[0]]==CHEM_ATOM_TYPE_HYDROGEN)&&(neighbors->list[neighbors->list[_i].list[1]].size>2)&&(mol->a[neighbors->list[_i].list[1]]!=CHEM_ATOM_TYPE_HYDROGEN) )||
       ( (mol->a[neighbors->list[_i].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(neighbors->list[neighbors->list[_i].list[0]].size>2)&&(mol->a[neighbors->list[_i].list[0]]!=CHEM_ATOM_TYPE_HYDROGEN) ) ) )
    {
    if (mol->a[neighbors->list[_i].list[0]]==CHEM_ATOM_TYPE_HYDROGEN) { _l=neighbors->list[_i].list[0], neighbors->list[_i].list[0]=neighbors->list[_i].list[1], neighbors->list[_i].list[1]=_l; }
    if (chrg_grps[_i]!=--_j) { _l=mol->natoms; while (_l--) if (chrg_grps[_l]==_j) chrg_grps[_l]=chrg_grps[_i]; }
    chrg_grps[_i]=chrg_grps[neighbors->list[_i].list[1]]=chrg_grps[neighbors->list[_i].list[0]];
    }

//Exit
return chrg_grps;
}



/*******************************************************   I M P L IC I T     S O L V E N T   *******************************************************************************/


//This function calculates the solvatation energy of the system (accordingly with Robert S. Rizzo, Tiba Aynechi, David A. Case, Irwin D. Kunts "Estimate of Absolute Free Energies of Gidratation Using Continium Methods: Accuracy of Partial Charge Models and Optimization of Nonpolar Contribution", JCTC 2006 2 128-139  )
double calc_esolv(double *sasa,double *sava,t_list *it,t_mol *mol,t_top *top)
{
register unsigned int _i,_j,_k;
double *a,*dij2,aa,qq,gg;
//Get square distance matrix
dij2=(double*)malloc(sizeof(double)*it->size*(it->size-1)/2);
a=(double*)malloc(sizeof(double)*it->size);
for (_k=0,_i=0;_i<it->size;_i++)
  for (_j=_i+1;_j<it->size;_j++,_k++)
    dij2[_k]=calc_distance(&mol->r[it->list[_i]],&mol->r[it->list[_j]]);
//Calc Borns radii
calc_Born_radii(a,dij2,sava,it,mol,top);
//Evaluate qq
for (gg=0.00,qq=0.00,_i=0;_i<it->size;_i++)
  {
  gg+=sasa[_i]*top->ff_a[mol->ytypes[it->list[_i]]].sgamma;
  qq+=sqrd(mol->charges[it->list[_i]])/a[_i];
  _k=_i-1;
  for (_j=0;_j<_i;_j++,_k+=it->size-_j-1)    //_k=_j*it->size-_j*(_j+3)/2+_i-1;
    {
    aa=a[_i]*a[_j];
    qq+=mol->charges[it->list[_i]]*mol->charges[it->list[_j]]/sqrt(dij2[_k]+aa*exp(-0.25*dij2[_k]/aa));
    }
  _k=_i*it->size-_i*(_i+1)/2;                                     //Set iterator to nonredundant part of distance matrix
  for (_j=_i+1;_j<it->size;_j++,_k++)
    {
    aa=a[_i]*a[_j];
    qq+=mol->charges[it->list[_i]]*mol->charges[it->list[_j]]/sqrt(dij2[_k]+aa*exp(-0.25*dij2[_k]/aa));
    }
  }
free(dij2);
free(a);

//printf("\n");
//for (_i=0;_i<it->size;_i++)
//  printf("%d type=%2.d i=%8.4f j=%8.4f k=%8.4f q=%8.4f a=%8.4f surf=%8.4f vol=%8.4f\n",it->list[_i],mol->ytypes[it->list[_i]],mol->rvecs[it->list[_i]].i,mol->rvecs[it->list[_i]].j,mol->rvecs[it->list[_i]].k,mol->charges[it->list[_i]],a[_i],sasa[it->list[_i]],sava[it->list[_i]]);
//printf("E1=%8.4f E2=%8.4f ",-166.00*(1.00-1.00/ESOLVENT)*qq,gg);
//             Polar                    Apolar
return -166.00*(1.00-1.00/ESOLVENT)*qq;// + gg;
}

//This function calculates the Born radiuses (accordingly with Di Qiu,Peter S. Shenkin,Frank P. Hollinger, and W. Clark Still "The GB/SA Continuum Model for Solvation. A Fast Analitical Method for the Calculation of Aproximate Born Radii" J.Phys.Chem. 1997, 101,3005-3014)
//Note. Vj is now replaced with Solvent Acessible Volume Area that calculated from delaunay triangulation.
void calc_Born_radii(double *rBorn,double *dij2,double *sava,t_list *it,t_mol *mol,t_top *top)
{
unsigned int _i,_j,_k;
double vrij4;

for (_i=0;_i<it->size;_i++)
  {
  vrij4=0.00;   //The summ of Vj/rij^4
  _k=_i-1;
  for (_j=0;_j<_i;_j++,_k+=it->size-_j-1) //  vrij4+=sava[_j]/sqrd(dij[_j*it->size-_j*(_j+3)/2+_i-1]);
    vrij4+=sava[_j]/sqrd(dij2[_k]);                                      //Scan columns of distance matrix. Note it is a squares matrix!!!
  _k=_i*it->size-_i*(_i+1)/2;                                          //Set iterator to nonredundant part of distance matrix
  for (_j++;_j<it->size;_j++,_k++)
    vrij4+=sava[_j]/sqrd(dij2[_k]);                                       //Scan rows of distance matrix. Note it is a squares matrix!!!
  rBorn[_i]=-166.0/(-166.0/(top->ff_a[mol->ytypes[it->list[_i]]].rvdw+RSOLVENT)+vrij4);
  }
}



