#include "y_cluster.h"


//This function performs clustering of set of points with k-mean cluster algorithm (ISOCLUS)
//Note _t is of size; clust, D, N and Vmax are of nclus*2+1, Dij is [nclus*2+1][nclus*2+1], V is [nclus*2+1][ndim].
unsigned int cluster_kmean(unsigned int ndim,unsigned int size,double *r,unsigned int niter,unsigned int nclus,double *clus,
                           unsigned int samprm,double stdv,double lump,unsigned int maxpair,unsigned int *_t,unsigned int *N,double *D,double **Dij,double *V,unsigned int *Vmax)
{
unsigned int _i, _j, _k, _l, _p, _q, size_c;
double _d, __d, Do;
//Step 1. Setup trial centers 
memcpy(&clus[0],   r,sizeof(double)*ndim); //Min coords
memcpy(&clus[ndim],r,sizeof(double)*ndim); //Max coords
_i=size; while (--_i) { _l=ndim; while (_l--) if (clus[_l]>r[_i*ndim+_l]) clus[_l]=r[_i*ndim+_l]; else if (clus[(nclus-1)*ndim+_l]<r[_i*ndim+_l]) clus[(nclus-1)*ndim+_l]=r[_i*ndim+_l]; }
size_c=nclus, _i=size_c-1; while (--_i) { _l=ndim; while (_l--) clus[_i*ndim+_l]=clus[_l]+(double)_i*(clus[(nclus-1)*ndim+_l]-clus[_l])/(double)(nclus-1); }

while (niter--)
  {
  //Step 2. Distribute samples among cluster centers
  _j=size_c; while (_j--) N[_j]=0;
  _i=size; 
  while (_i--) 
    {
    _t[_i]=0, _d=0., _l=ndim; while (_l--) _d+=sqrd(clus[_l]-r[_i*ndim+_l]);
    _j=size_c; 
    while (--_j) 
      {
      __d=0., _l=ndim; while (_l--) __d+=sqrd(clus[_j*ndim+_l]-r[_i*ndim+_l]);
      if (__d<_d) { _d=__d, _t[_i]=_j; }
      }
    N[_t[_i]]++;
    }
  //Step 3. Discard samples with amount of points less then samprm
  _j=size_c; 
  while ((_j--)&&(size_c)) 
    if (N[_j]<samprm)
      {
      size_c--,  N[_j]=N[size_c], _l=ndim; while (_l--) clus[_j*ndim+_l]=clus[size_c*ndim+_l]; 
      _i=size; while (_i--) if (_t[_i]==_j) _t[_i]=(unsigned int)-1; else if (_t[_i]==size_c) _t[_i]=_j;
      } 
  //Step 4. Update cluster centers
  _j=size_c; while (_j--) { _l=ndim; while (_l--) clus[_j*ndim+_l]=0.; }
  _i=size; while (_i--) if (_t[_i]!=(unsigned int)-1) { _l=ndim; while (_l--) clus[_t[_i]*ndim+_l]+=r[_i*ndim+_l]; }
  _j=size_c; while (_j--) { _l=ndim; while (_l--) clus[_j*ndim+_l]/=(double)N[_j]; }
  //Stage 5. Compute average distance from center of cluster to its members
  _i=size; while (_i--) if (_t[_i]!=(unsigned int)-1) { _d=0., _l=ndim; while (_l--) _d+=sqrd(clus[_t[_i]*ndim+_l]-r[_i*ndim+_l]); D[_t[_i]]+=sqrt(_d); }
  _j=size_c; while (_j--) D[_j]/=(double)N[_j]; 
  //Stage 6. Compurte overall average distance from samples to ceneter of cluster
  Do=0., _i=0, _j=size_c; while (_j--) { Do+=D[_j]*(double)N[_j], _i+=N[_j]; } Do/=(double)_i;
  //Stage 7. Switching
  if (!niter) { lump=0; goto STEP_11; } 
  if ( ( (!(niter%2))&&(nclus/2<size_c) )||(3*nclus/2<=size_c) ) goto STEP_11;
  //Stage 8. Find standard deviation for each cluster 
  _j=size_c; while (_j--) { _l=ndim; while (_l--) V[_j*ndim+_l]=0.; }
  _i=size; while (_i--) if (_t[_i]!=(unsigned int)-1) { _l=ndim; while (_l--) V[_t[_i]*ndim+_l]+=sqrd(clus[_t[_i]*ndim+_l]-r[_i*ndim+_l]); }
  _j=size_c; while (_j--) { _l=ndim; while (_l--) V[_j*ndim+_l]=sqrt(V[_j*ndim+_l]/(double)N[_j]); }
  //Stage 9. Find maximal components of deviation vectors
  _j=size_c; while (_j--) { Vmax[_j]=0, _l=ndim; while (--_l) if (V[_j*ndim+Vmax[_j]]<V[_j*ndim+_l]) Vmax[_j]=_l; }
  //Stage 10. Split clusters
  _d=FALSE, _j=size_c; 
  while ( (_j--)&&(size_c<nclus) )
    if ( (V[_j*ndim+Vmax[_j]]>stdv)&&(D[_j]>Do)&&(N[_j]>2*samprm+1) )
      {
      _l=ndim; while (_l--) clus[size_c*ndim+_l]=clus[_j*ndim+_l];
      clus[size_c*ndim+Vmax[_j]]+=.5*V[_j*ndim+Vmax[_j]], V[_j*ndim+Vmax[_j]]-=.5*V[_j*ndim+Vmax[_j]];
      size_c++; 
      _d=TRUE;
      } 
  if ((_d)) continue;
  //Step 11. Compute distances between all cluster centers
  STEP_11: 
  _j=size_c;
  while (_j--)
    {
    _i=_j;
    while (_i--)
      {
      Dij[_i][_j]=0., _l=ndim; while (_l--) Dij[_i][_j]+=sqrd(clus[_j*ndim+_l]-clus[_i*ndim+_l]);
      Dij[_i][_j]=Dij[_j][_i]=sqrt(Dij[_i][_j]);
      }
    }
  //Step 12+13. Lump clusters 
  _k=0, _i=size_c;
  while ((_i--)&&((!_k)||(_k<nclus/10)))
    {
    _j=_i;
    while (_j--)
      if (Dij[_i][_j]<lump) 
        {//merging pair avialable
        //Find closest i,j pair
        _q=_j; while (_q--) if (Dij[_i][_q]<Dij[_i][_j]) _j=_q;
        _p=_i; while (_p--) { _q=_p; while (_q--) if (Dij[_p][_q]<Dij[_i][_j]) { _i=_p, _j=_q; } }
        //Split centers (i+j)/2->j
        _l=ndim; while (_l--) clus[_j*ndim+_l]=((double)N[_i]*clus[_i*ndim+_l]+(double)N[_j]*clus[_j*ndim+_l])/(double)(N[_i]+N[_j]); N[_j]+=N[_i];
        //Do i <- size_c-_k-1
        N[_i]=N[size_c-_k-1], _l=ndim; while (_l--) clus[_i*ndim+_l]=clus[(size_c-_k-1)*ndim+_l]; 
        _l=size_c-_k-1; while (_l--) { Dij[_i][_l]=Dij[size_c-_k-1][_l], Dij[_l][_i]=Dij[_l][size_c-_k-1]; } 
        //Do size-_k <- --size_c
        size_c--;
        N[size_c-_k]=N[size_c], _l=ndim; while (_l--) clus[(size_c-_k)*ndim+_l]=clus[size_c*ndim+_l]; 
        //Do j <-> size_c-_k-2
        _k++;
        _p=N[_j], N[_j]=N[size_c-_k], N[size_c-_k]=_p, _l=ndim; while (_l--) { _d=clus[_i*ndim+_l], clus[_i*ndim+_l]=clus[(size_c-_k)*ndim+_l], clus[(size_c-_k)*ndim+_l]=_d; }
        //do j <- size_c-_k-2
        _l=size_c-_k; while (_l--) { Dij[_j][_l]=Dij[size_c-_k][_l], Dij[_l][_j]=Dij[_l][size_c-_k]; }
        //Restart
        _i=size_c-_k;
        break; 
        } 
    }
  }
return size_c;  
}


//This function clusters the grid
t_cluster3D *cluster_kmean_dgrid(unsigned int nsteps,unsigned int size_cl,unsigned int lcutoff,double stdv,double lump,t_len *len,double ***x)
{
t_cluster3D *cl=0x0;
unsigned int _i, _j, _k, _l, _t, *num=0x0;
t_vec *v=0x0, *v2=0x0;
double *d=0x0, _d, __d, ___d;

//Step 0. Alloc memory for clusters
if (!size_cl) { LABEL_USER_ERROR: ylib_errno=YERROR_USER; return FALSE; }
if (!(cl=(t_cluster3D*)malloc(sizeof(t_cluster3D)+sizeof(t_lvec)*(size_cl*2)))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
else cl->size=size_cl, cl->_size=size_cl*2, cl->lvec=(void*)cl+sizeof(t_cluster3D);
if (!(d=(double*)malloc(sizeof(double)*cl->_size)))                                           { free(cl); goto LABEL_MEMORY_ERROR; }
if (!(v=(t_vec*)malloc(sizeof(t_vec)*cl->_size)))                                    { free(d); free(cl); goto LABEL_MEMORY_ERROR; } //+1 required to proper qsort
if (!(v2=(t_vec*)malloc(sizeof(t_vec)*cl->_size)))                          { free(v); free(d); free(cl); goto LABEL_MEMORY_ERROR; }
if (!(num=(unsigned int*)malloc(sizeof(unsigned int)*cl->_size))) { free(v2); free(v); free(d); free(cl); goto LABEL_MEMORY_ERROR; }

//Step 1. Assign random size_cl clusters over the set. The centers are assigned into equal boxes
_i=(len->j*len->k), _l=len->i*_i/size_cl, cl->lvec[0].i=cl->lvec[0].j=cl->lvec[0].k=0.;
_t=size_cl; while (--_t) { _k=(_t*_l), cl->lvec[_t].i=(double)(_k/_i), _k=_k%_i, cl->lvec[_t].j=(double)(_k/len->k), cl->lvec[_t].k=(double)(_k%len->k); }

//Do main cycle
lump*=lump;
while (nsteps--)
  {
  //Step 1. Distribute samples among clusters
  _l=cl->size; while (_l--) { num[_l]=0, cl->lvec[_l].w=v[_l].i=v[_l].j=v[_l].k=0.; }
  for (_i=0;_i<len->i;_i++)
    for (_j=0;_j<len->j;_j++)
      for (_k=0;_k<len->k;_k++)
        if ( (!(isnan(x[_i][_j][_k])))&&( (x[_i][_j][_k])) ) 
          {
          _t=0, _d=sqrd(cl->lvec[0].i-(double)_i)+sqrd(cl->lvec[0].j-(double)_j)+sqrd(cl->lvec[0].k-(double)_k);  
          _l=cl->size; while (--_l) if ((__d=(sqrd(cl->lvec[_l].i-(double)_i)+sqrd(cl->lvec[_l].j-(double)_j)+sqrd(cl->lvec[_l].k-(double)_k)))<_d) { _t=_l, _d=__d; }
          cl->lvec[_t].w+=x[_i][_j][_k], v[_t].i+=x[_i][_j][_k]*(double)_i, v[_t].j+=x[_i][_j][_k]*(double)_j, v[_t].k+=x[_i][_j][_k]*(double)_k;  
          num[_t]++;
          }
  //Step 2. Discard small clusters
  _t=_l=cl->size; 
  while (_l--) 
    if (num[_l]<lcutoff)
      {
      while (--cl->size!=_l) if (num[cl->size]>lcutoff) break;
      cl->lvec[_l].i=cl->lvec[cl->size].i,  cl->lvec[_l].j=cl->lvec[cl->size].j, cl->lvec[_l].k=cl->lvec[cl->size].k, cl->lvec[_l].w=cl->lvec[cl->size].w;
      v[_l].i=v[cl->size].i, v[_l].j=v[cl->size].j, v[_l].k=v[cl->size].k, num[_l]=num[cl->size];
      }
  if (!cl->size) { free(cl); free(v); free(d); free(num); goto LABEL_USER_ERROR; }
  //Step 3. Update centers
 _l=cl->size; while (_l--) { v[_l].i/=cl->lvec[_l].w, v[_l].j/=cl->lvec[_l].w, v[_l].k/=cl->lvec[_l].w; }
  if (cl->size==_t)
    {
    //Step 6. Switch startegy
    if (!nsteps) continue;
    if ( (nsteps>0xF)&&( (cl->size<size_cl/2)||(!(nsteps%2)) ) )
      {
      //Step 4. Compute average distances
      _l=cl->size; while (_l--) { d[_l]=v2[_l].i=v2[_l].j=v2[_l].k=0.; }
      for (_i=0;_i<len->i;_i++)
        for (_j=0;_j<len->j;_j++)
          for (_k=0;_k<len->k;_k++)
            if ( (!(isnan(x[_i][_j][_k])))&&( (x[_i][_j][_k])) ) 
              {
              _t=0, _d=sqrd(cl->lvec[0].i-(double)_i)+sqrd(cl->lvec[0].j-(double)_j)+sqrd(cl->lvec[0].k-(double)_k);  
              _l=cl->size; while (--_l) if ((__d=(sqrd(cl->lvec[_l].i-(double)_i)+sqrd(cl->lvec[_l].j-(double)_j)+sqrd(cl->lvec[_l].k-(double)_k)))<_d) { _t=_l, _d=__d; }
              d[_t]+=x[_i][_j][_k]*sqrt(sqrd(v[_t].i-(double)_i)+sqrd(v[_t].j-(double)_j)+sqrd(v[_t].k-(double)_k));
              v2[_t].i+=x[_i][_j][_k]*sqrd(v[_t].i-(double)_i), v2[_t].j+=x[_i][_j][_k]*sqrd(v[_t].j-(double)_j), v2[_t].k=x[_i][_j][_k]*sqrd(v[_t].k-(double)_k);
              }
      //Step 5. Calculate new centers of clusters 
      _d=0., _l=cl->size; while (_l--) { _d+=d[_l], d[_l]/=cl->lvec[_l].w, cl->lvec[_l].i=v[_l].i, cl->lvec[_l].j=v[_l].j, cl->lvec[_l].k=v[_l].k; } 
      _d/=(double)cl->size;

      //Step 7. Find principal components
      _l=cl->size; while (_l--) { v2[_l].i/=cl->lvec[_l].w, v2[_l].j/=cl->lvec[_l].w, v2[_l].k/=cl->lvec[_l].w; }
      //Step 8. Separete clusters
      _l=cl->size;
      while (_l--) 
        if ( ( (v2[_l].i>stdv)||(v2[_l].j>stdv)||(v2[_l].k>stdv) )&&( ( (d[_l]>_d)&&(num[_l]>2*lcutoff+1) )||(cl->size<=size_cl/2) ) )  
          {
          //Step 9. Split the cluster along maximum deviation direction
             if (v2[_l].i>v2[_l].j)
               {
               if (v2[_l].i>v2[_l].k) { cl->lvec[cl->size].i=cl->lvec[_l].i+d[_l], cl->lvec[_l].i-=d[_l], cl->lvec[cl->size].j=cl->lvec[_l].j, cl->lvec[cl->size].k=cl->lvec[_l].k; }
               else                   { cl->lvec[cl->size].i=cl->lvec[_l].i, cl->lvec[cl->size].j=cl->lvec[_l].j, cl->lvec[cl->size].k=cl->lvec[_l].k+d[_l], cl->lvec[_l].k-=d[_l]; } 
               }
          else if (v2[_l].j>v2[_l].k) { cl->lvec[cl->size].i=cl->lvec[_l].i, cl->lvec[cl->size].j=cl->lvec[_l].j+d[_l], cl->lvec[_l].j-=d[_l], cl->lvec[cl->size].k=cl->lvec[_l].k; } 
               else                   { cl->lvec[cl->size].i=cl->lvec[_l].i, cl->lvec[cl->size].j=cl->lvec[_l].j, cl->lvec[cl->size].k=cl->lvec[_l].k+d[_l], cl->lvec[_l].k-=d[_l]; } 
          cl->size++;
          }
      } 
    else if ( (nsteps>0xF)&&(cl->size>1)&&( (cl->size>=3*size_cl/2)||( (nsteps%2)) ) ) 
           {
           //Step 5. Calculate new centers of clusters 
           _d=0., _l=cl->size; while (_l--) { cl->lvec[_l].i=v[_l].i, cl->lvec[_l].j=v[_l].j, cl->lvec[_l].k=v[_l].k; }   

           //Step 10, 11, 12. Compute pair-wise distances between all clusters and join them  
           _l=cl->size; while (_l--) d[_l]=0;
           __d=(double)NAN, *num=size_cl/10+1;
           while ((*num)--)
             {
             for (_d=(double)NAN,_t=_i=0;_i<cl->size;_i++)
               for (_j=_i+1;_j<cl->size;_j++)
                 if ( ( ( (isnan( _d)))||( _d<(___d=sqrd(cl->lvec[_i].i-cl->lvec[_j].i)+sqrd(cl->lvec[_i].j-cl->lvec[_j].j)+sqrd(cl->lvec[_i].k-cl->lvec[_j].k))) )&&
                      ( ( (isnan(__d)))||(__d<___d)) ) { _k=_i, _l=_j, _d=___d; }
             if (_d>lump) break;
             if ( (!(isnan(d[_k])))&&(!(isnan(d[_l]))) )
               {
               if (_k>_l) { _i=_k, _k=_l, _l=_i; }
               ___d=1./cl->lvec[_k].w+cl->lvec[_l].w;
               cl->lvec[_k].i=__d*(cl->lvec[_k].w*cl->lvec[_k].i+cl->lvec[_l].w*cl->lvec[_l].i), cl->lvec[_k].j=__d*(cl->lvec[_k].w*cl->lvec[_k].j+cl->lvec[_l].w*cl->lvec[_l].j), cl->lvec[_k].k=__d*(cl->lvec[_k].w*cl->lvec[_k].k+cl->lvec[_l].w*cl->lvec[_l].k);
               if (--cl->size!=_l) cl->lvec[_l].i=cl->lvec[cl->size].i, cl->lvec[_l].j=cl->lvec[cl->size].j, cl->lvec[_l].k=cl->lvec[cl->size].k, cl->lvec[_l].w=cl->lvec[cl->size].w, d[_l]=d[cl->size];
               d[_k]=(double)NAN;
               }
             else __d=_d;
             }
           }
    }
  else { _l=cl->size; while (_l--) { cl->lvec[_l].i=v[_l].i, cl->lvec[_l].j=v[_l].j, cl->lvec[_l].k=v[_l].k; } }
  }
//Free some memory and exit
free(num); free(d); free(v); free(v2); 
return cl;
}


//Utility function for sorting of clusters distnces
int compare_cluster_3Dpair(const void *va,const void *vb)
{
if ( ((t_vec*)va)->k<((t_vec*)vb)->k) return -1;
else return (int)((t_vec*)va)->k!=((t_vec*)vb)->k;
}

//This function do clustering with ISOCLUS algorithm (variaton of unweighted k-means )
//Note. t is a working massive of fitted ints.
t_cluster3D *cluster_kmean_grid(unsigned int nsteps,unsigned int size_cl,unsigned int lcutoff,double stdv,double lump,t_len *len,char ***x)
{
t_cluster3D *cl=0x0;
unsigned int _i, _j, _k, _l, _t;
t_vec r, *v=0x0, *v2=0x0;
double *d=0x0, _d, _rr;
size_t n;
void *vp;

//Step 0. Alloc memory for clusters
n=(size_t)len->i*(size_t)len->j*(size_t)len->k;
if (!size_cl) { LABEL_USER_ERROR: ylib_errno=YERROR_USER; return FALSE; }
if (!(cl=(t_cluster3D*)malloc(sizeof(t_cluster3D)+sizeof(t_lvec)*(_i=(size_cl<0x8F) ? 0xFF : size_cl*2+1)))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
else cl->size=0x0, cl->_size=_i, cl->lvec=(void*)cl+sizeof(t_cluster3D);
if (!(d=(double*)malloc(sizeof(double)*cl->_size)))                  { free(cl); goto LABEL_MEMORY_ERROR; }
if (!(v=(t_vec*)malloc(sizeof(t_vec)*(cl->_size+1)))) {       free(d); free(cl); goto LABEL_MEMORY_ERROR; } //+1 required to proper qsort
if (!(v2=(t_vec*)malloc(sizeof(t_vec)*cl->_size))) { free(v); free(d); free(cl); goto LABEL_MEMORY_ERROR; }
//Step 1. Assign random size_cl clusters over the set. The centers are assigned into equal boxes
_l=(unsigned int)cube_root((double)size_cl);
r.i=(double)len->i/(double)_l, r.j=(double)len->j/(double)_l, r.k=(double)len->k/(double)_l;
for (_i=0;_i<_l;_i++)
  for (_j=0;_j<_l;_j++)
    for (_k=0;_k<_l;_k++)
      {
      cl->lvec[cl->size].i=r.i*(0.50+(double)_i);
      cl->lvec[cl->size].j=r.j*(0.50+(double)_j);
      cl->lvec[cl->size].k=r.k*(0.50+(double)_k);
      cl->lvec[cl->size].w=0.00;
      cl->size++;
      }
//Calculate number of items in system
n=0;
_i=len->i;
while (_i--)
  {
  _j=len->j;
  while (_j--)
    {
    _k=len->k;
    while (_k--)
      if ((x[_i][_j][_k])) n++;
    }
  }

//Do executional hypercycle
lump*=lump;
while (--nsteps)
  {
  //Step 2.  Assign points to nearest clusters XinSj if  ||X-Zj|| < ||X<Zi||   i=1...Nc, j!=i
  STEP_2:
  _j=cl->size; while (_j--) { cl->lvec[_j].w=v[_j].i=v[_j].j=v[_j].k=v2[_j].i=v2[_j].j=v2[_j].k=0.00; }
  //Scan grid using hash table of neares neigbors along each axis
  _i=len->i;
  while (_i--)
    {
    _j=len->j;
    while (_j--)
      {
      _k=len->k;
      while (_k--)
        {
        if (!x[_i][_j][_k]) continue;
        _t=0;
        _d=sqrd(cl->lvec[0].i-(double)_i)+sqrd(cl->lvec[0].j-(double)_j)+sqrd(cl->lvec[0].k-(double)_k);
        _l=cl->size;
        while (--_l)
          if ((_rr=sqrd(cl->lvec[_l].i-(double)_i)+sqrd(cl->lvec[_l].j-(double)_j)+sqrd(cl->lvec[_l].k-(double)_k))<_d) { _d=_rr; _t=_l; }
        cl->lvec[_t].w+=1.00;
        v[_t].i+=(double)_i, v[_t].j+=(double)_j, v[_t].k+=(double)_k;    //Summ members
        v2[_t].i+=sqrd(cl->lvec[_t].i-(double)_i), v2[_t].j+=sqrd(cl->lvec[_t].j-(double)_j), v2[_t].k+=sqrd(cl->lvec[_t].k-(double)_k);  //Summ members squares
        }
      }
    }
  //Step 3. Discard small enought clusters   if (Nj<lcutoff) then  discard Sj and Nc=Nc-1
  _t=cl->size;
  for (_j=0;_j<cl->size;_j++)
    if ((unsigned int)cl->lvec[_j].w<lcutoff)
      {
      REMOVE_SMALL_CLUSTERS: if (--cl->size==_j) break; //Goto used to jump out of two cycles with one 'break'
                             if (cl->lvec[cl->size].w<lcutoff) goto REMOVE_SMALL_CLUSTERS;
      //delete cluster _j
      cl->lvec[_j].i=cl->lvec[cl->size].i;
      cl->lvec[_j].j=cl->lvec[cl->size].j;
      cl->lvec[_j].k=cl->lvec[cl->size].k;
      cl->lvec[_j].w=cl->lvec[cl->size].w;
      }
  if (cl->size!=_t)
    {
    if (!cl->size) { free(cl); free(d); goto LABEL_USER_ERROR; }
    else goto STEP_2; //Reassign samples
    }
  //Step 4. Compute new clusters centers Zj=1/Nj*SUMM[XinSj](X)   j=1...Nc
  _j=cl->size;
  while (_j--)
    {
    v[_j].i/=cl->lvec[_j].w, v[_j].j/=cl->lvec[_j].w, v[_j].k/=cl->lvec[_j].w;
    _d=v[_j].i, v[_j].i=cl->lvec[_j].i, cl->lvec[_j].i=_d;
    _d=v[_j].j, v[_j].j=cl->lvec[_j].j, cl->lvec[_j].j=_d;
    _d=v[_j].k, v[_j].k=cl->lvec[_j].k, cl->lvec[_j].k=_d;
    }
  //Step 5. Compute the average distance in cluster domain from their center  Dj=1/Nj*SUMM[XinSj](||X-Zj||)   j=1...Nc
  //Step 6. Compute average distance of the samples from their centers Do=1/N*SUMM[1...Nc](Nj*Dj)    j=1...Nc
  _d=0.00;
  _j=cl->size;
  while (_j--)
    {
     d[_j]=(v2[_j].i+v2[_j].j+v2[_j].k)/cl->lvec[_j].w-sqrd(cl->lvec[_j].i-v[_j].i)-sqrd(cl->lvec[_j].j-v[_j].j)-sqrd(cl->lvec[_j].k-v[_j].k);
    _d+=d[_j]*cl->lvec[_j].w;
    }
  _d/=3.00*(double)n;
  //Step 7. If this is the last iteration set LUMP=0 and goto step 11. Otherwise if Nc<=size_cl/2 goto step 8. Otherwise if even-numbered itration  or if Nc>=size_cl*2 goto step 11. Oterwise continue;
       if ((nsteps)&&(cl->size<=size_cl/2)) goto STEP_8;
  else if ((nsteps)&&((!nsteps%2)||(cl->size>=size_cl*2))) goto STEP_11;
  else continue;
  //Step 8. Find standart deviation vector for each samples subset    Vij=sqrt( 1/Nj*SUMM((Xik-Zij)**2)  )
  STEP_8:
  _j=cl->size; while (_j--) { v[_j].i=v2[_j].i/cl->lvec[_j].w-sqrd(cl->lvec[_j].i-v[_j].i), v[_j].j=v2[_j].j/cl->lvec[_j].w-sqrd(cl->lvec[_j].j-v[_j].j), v[_j].k=v2[_j].k/cl->lvec[_j].w-sqrd(cl->lvec[_j].k-v[_j].k); }
  //Step 9.   Find the maximal deviation (joined with step 10 for efficiency)
  //Step 10. Splite cluster  If vmax[_j]>average_STDV &&  ( ( (Dj>Do)&&(Nj>2*lcutoff+1) ) || (Nc<=size_cl/2) ) then split Zj into Zj+ and Zj-
  _l=_j=cl->size;
  while (_j--)
    if (v[_j].i>v[_j].j)
      {
      if (v[_j].i>v[_j].k)
        {//Split along X
        if ( (v[_j].i>stdv) && (v[_j].i>_d) && ( (cl->size<=size_cl/2) || ( (d[_j]>_d) && ((unsigned int)cl->lvec[_j].w>2*lcutoff+1) ) ) )
          {//Split
          if (cl->size==cl->_size)
            {
            if ( (vp=(void*)realloc(cl,sizeof(t_cluster3D)+sizeof(t_lvec)*(cl->_size+=0xF))))
              { LABEL_REALLOC_FAIL: free(cl); free(d); free(v); free(v2); goto LABEL_MEMORY_ERROR; }
            else cl=(t_cluster3D*)vp;
            if ( (vp=(void*)realloc(d,sizeof(double)*cl->_size))) goto LABEL_REALLOC_FAIL;
            else d=(double*)vp;
            if (!(vp=(void*)realloc(v,sizeof(t_vec)*(cl->_size+1)))) goto LABEL_REALLOC_FAIL; //+1 required to proper qsort
            else v=(t_vec*)vp;
            if (!(vp=(void*)realloc(v2,sizeof(t_vec)*(cl->_size)))) goto LABEL_REALLOC_FAIL;
            else v2=(t_vec*)vp;
            }
          v[_j].i=0.50*sqrt(v[_j].i);
          cl->lvec[cl->size].i=cl->lvec[_j].i+v[_j].i, cl->lvec[cl->size].j=cl->lvec[_j].j, cl->lvec[cl->size].k=cl->lvec[_j].k;
          cl->lvec[_j].i-=v[_j].i;
          cl->size++;
          }
        }
      else goto SPLIT_Z;
      }
    else if (v[_j].j>v[_j].k)
           {//Split along Y
           if ( (v[_j].j>stdv) && (v[_j].j>_d) && ( (cl->size<=size_cl/2) || ( (d[_j]>_d) && ((unsigned int)cl->lvec[_j].w>2*lcutoff+1) ) ) )
             {//Split
             if (cl->size==cl->_size)
               {
               if ( (vp=(void*)realloc(cl,sizeof(t_cluster3D)+sizeof(t_lvec)*(cl->_size+=0xF)))) goto LABEL_REALLOC_FAIL;
               else cl=(t_cluster3D*)vp;
               if ( (vp=(void*)realloc(d,sizeof(double)*cl->_size))) goto LABEL_REALLOC_FAIL;
               else d=(double*)vp;
               if (!(vp=(void*)realloc(v,sizeof(t_vec)*(cl->_size+1)))) goto LABEL_REALLOC_FAIL; //+1 required to proper qsort
               else v=(t_vec*)vp;
               if (!(vp=(void*)realloc(v2,sizeof(t_vec)*(cl->_size)))) goto LABEL_REALLOC_FAIL;
               else v2=(t_vec*)vp;
               }
             v[_j].j=0.50*sqrt(v[_j].j);
             cl->lvec[cl->size].i=cl->lvec[_j].i, cl->lvec[cl->size].j=cl->lvec[_j].j+v[_j].j, cl->lvec[cl->size].k=cl->lvec[_j].k;
             cl->lvec[_j].j-=v[_j].j;
             cl->size++;
             }
           }
         else
           {//Split along Z
           SPLIT_Z:
           if ( (v[_j].k>stdv) && (v[_j].k>_d) && ( (cl->size<=size_cl/2) || ( (d[_j]>_d) && ((unsigned int)cl->lvec[_j].w>2*lcutoff+1) ) ) )
             {//Split
             if (cl->size==cl->_size)
               {
               if ( (vp=(void*)realloc(cl,sizeof(t_cluster3D)+sizeof(t_lvec)*(cl->_size+=0xF)))) goto LABEL_REALLOC_FAIL;
               else cl=(t_cluster3D*)vp;
               if ( (vp=(void*)realloc(d,sizeof(double)*cl->_size))) goto LABEL_REALLOC_FAIL;
               else d=(double*)vp;
               if (!(vp=(void*)realloc(v,sizeof(t_vec)*(cl->_size+1)))) goto LABEL_REALLOC_FAIL; //+1 required to proper qsort
               else v=(t_vec*)vp;
               if (!(vp=(void*)realloc(v2,sizeof(t_vec)*(cl->_size)))) goto LABEL_REALLOC_FAIL;
               else v2=(t_vec*)vp;
               }
             v[_j].k=0.50*sqrt(v[_j].k);
             cl->lvec[cl->size].i=cl->lvec[_j].i, cl->lvec[cl->size].j=cl->lvec[_j].j, cl->lvec[cl->size].k=cl->lvec[_j].k+v[_j].k;
             cl->lvec[_j].k-=v[_j].k;
             cl->size++;
             }
           }
  if (_l!=cl->size) continue; //Splits occured
  //Step 11. Compute pairwise distances between cluster centers
  //Step 12. Compare distances with lump and arrange them in ascending order. Note we gather no more than cl->_size closest clusters
  STEP_11:
  _l=0;
  _j=cl->size;
  while (_j--)
    {
    _i=_j;
    while (_i--)
      {
      if ((_d=sqrd(cl->lvec[_j].i-cl->lvec[_i].i)+sqrd(cl->lvec[_j].j-cl->lvec[_i].j)+sqrd(cl->lvec[_j].k-cl->lvec[_i].k))<lump)
        {
        if (_l<cl->_size)
          { v[_l].i=(double)_i, v[_l].j=(double)_j, v[_l].k=_d, _l++;	} //Just add new member
        else //replace biggest
          {
          _t=0;
          _k=cl->_size;
          while (--_k)
            if (v[_t].k<v[_k].k) _t=_k;
          }
        if (v[_t].k>_d) { v[_t].i=(double)_i, v[_t].j=(double)_j, v[_t].k=_d; } //replace it
        }
      }
    }
  qsort(v,(size_t)_l,sizeof(t_vec),compare_cluster_3Dpair);
  //Step 13. Merge clusters
  _t=_l;
  for (_k=0;_k<_l;_k++)
    {
    //Merge
    _i=(unsigned int)v[_k].i;
    _j=(unsigned int)v[_k].j;
    _d=(cl->lvec[_i].w+cl->lvec[_j].w);
    cl->lvec[_i].i=(cl->lvec[_i].w*cl->lvec[_i].i+cl->lvec[_j].w*cl->lvec[_j].i)/_d;
    cl->lvec[_i].j=(cl->lvec[_i].w*cl->lvec[_i].j+cl->lvec[_j].w*cl->lvec[_j].j)/_d;
    cl->lvec[_i].k=(cl->lvec[_i].w*cl->lvec[_i].k+cl->lvec[_j].w*cl->lvec[_j].k)/_d;
    cl->lvec[_i].w=_d;
    //Del _j-th cluster
    if (_j!=--cl->size)
      {
      cl->lvec[_j].i=cl->lvec[cl->size].i, cl->lvec[_j].j=cl->lvec[cl->size].j, cl->lvec[_j].k=cl->lvec[cl->size].k, cl->lvec[_j].w=cl->lvec[cl->size].w;
      for (_i=_k+1;_i<_l;_i--)
        if ( (v[_k].i==v[_i].i)||(v[_k].i==v[_i].j)||(v[_k].j==v[_i].i)||(v[_k].j==v[_i].j) )
          {
          if (--_l!=_i) { v[_i].i=v[_l].i, v[_i].j=v[_l].k, v[_i].k=v[_l].k, _i--; }
          }
        else if ((unsigned int)v[_i].i==cl->size) v[_i].i=v[_k].j;
        else if ((unsigned int)v[_i].j==cl->size) v[_i].j=v[_k].j;
      }
    }
  if (_t!=_l) goto STEP_2;
  //Step 14. if last terminal other-wise continue
  }

//Exit
free(d);
free(v);
return cl;
}


//This function do clustering with ISOCLUS algorithm (variaton of unweighted k-means )
//Note. t is a working massive of fitted ints.
t_cluster3D *cluster_kmean_3D(unsigned int nsteps,unsigned int size_cl,unsigned int lcutoff,double stdv,double lump,unsigned int n,t_vec *x)
{
t_cluster3D *cl=0x0;
unsigned int _i, _j, _k, _l, _t;
t_vec r_min, r_max, *v=0x0, *v2=0x0;
double *d=0x0, _d, _rr;
void *vp;

//Step 0. Alloc memory for clusters
if (!size_cl) { LABEL_USER_ERROR: ylib_errno=YERROR_USER; return FALSE; }
if (!(cl=(t_cluster3D*)malloc(sizeof(t_cluster3D)+sizeof(t_lvec)*(_i=(size_cl<0x8F) ? 0xFF : size_cl*2+1)))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
else cl->size=0x0, cl->_size=_i, cl->lvec=(void*)cl+sizeof(t_cluster3D);
if (!(d=(double*)malloc(sizeof(double)*cl->size)))                     { free(cl); goto LABEL_MEMORY_ERROR; }
if (!(v=(t_vec*)malloc(sizeof(t_vec)*(cl->size+1)))) {          free(d); free(cl); goto LABEL_MEMORY_ERROR; } //+1 required to proper qsort
if (!(v2=(t_vec*)malloc(sizeof(t_vec)*cl->size)))    { free(v); free(d); free(cl); goto LABEL_MEMORY_ERROR; } //+1 required to proper qsort
//Step 1. Assign random size_cl clusters over the set. The centers are assigned into equal boxes
r_min.i=r_max.i=x[0].i;
r_min.j=r_max.j=x[0].j;
r_min.k=r_max.k=x[0].k;
_i=n;
while(_i--)
  {
       if (r_min.i>x[_i].i) r_min.i=x[_i].i;
  else if (r_max.i<x[_i].i) r_max.i=x[_i].i;
       if (r_min.j>x[_i].j) r_min.j=x[_i].j;
  else if (r_max.j<x[_i].j) r_max.j=x[_i].j;
       if (r_min.k>x[_i].k) r_min.k=x[_i].k;
  else if (r_max.k<x[_i].k) r_max.k=x[_i].k;
  }
_l=(unsigned int)cube_root((double)size_cl);
r_max.i=(r_min.i-r_max.i)/(double)_l, r_max.j=(r_min.j-r_max.j)/(double)_l, r_max.k=(r_min.k-r_max.k)/(double)_l;
for (_i=0;_i<_l;_i++)
  for (_j=0;_j<_l;_j++)
    for (_k=0;_k<_l;_k++)
      {
      cl->lvec[cl->size].i=r_max.i*(0.50+(double)_i)+r_min.i;
      cl->lvec[cl->size].j=r_max.j*(0.50+(double)_j)+r_min.j;
      cl->lvec[cl->size].k=r_max.k*(0.50+(double)_k)+r_min.k;
      cl->lvec[cl->size].w=0.00;
      cl->size++;
      }
//Do executional hypercycle
lump*=lump;
while (--nsteps)
  {
  //Step 2.  Assign points to nearest clusters XinSj if  ||X-Zj|| < ||X<Zi||   i=1...Nc, j!=i
  STEP_2:
  _j=cl->size; while (_j--) { cl->lvec[_j].w=v[_j].i=v[_j].j=v[_j].k=v2[_j].i=v2[_j].j=v2[_j].k=0.00; }
  //Scan grid using hash table of neares neigbors along each axis
  _i=n;
  while (_i--)
    {
    _t=0;
    _d=sqrd(cl->lvec[0].i-x[_i].i)+sqrd(cl->lvec[0].j-x[_i].j)+sqrd(cl->lvec[0].k-x[_i].k);
    _l=cl->size;
    while (--_l)
      if ((_rr=sqrd(cl->lvec[_l].i-x[_i].i)+sqrd(cl->lvec[_l].j-x[_i].j)+sqrd(cl->lvec[_l].k-x[_i].k))<_d) { _d=_rr; _t=_l; }
    cl->lvec[_t].w+=1.00;
    v[_t].i+=x[_i].i, v[_t].j+=x[_i].j, v[_t].k+=x[_i].k;    //Summ members
    v2[_t].i+=sqrd(cl->lvec[_t].i-x[_i].i), v2[_t].j+=sqrd(cl->lvec[_t].j-x[_i].j), v2[_t].k+=sqrd(cl->lvec[_t].k-x[_i].k);  //Summ members squares
    }
  //Step 3. Discard small enought clusters   if (Nj<lcutoff) then  discard Sj and Nc=Nc-1
  _t=cl->size;
  for (_j=0;_j<cl->size;_j++)
    if ((unsigned int)cl->lvec[_j].w<lcutoff)
      {
      REMOVE_SMALL_CLUSTERS: if (--cl->size==_j) break; //Goto used to jump out of two cycles with one 'break'
                             if (cl->lvec[cl->size].w<lcutoff) goto REMOVE_SMALL_CLUSTERS;
      //delete cluster _j
      cl->lvec[_j].i=cl->lvec[cl->size].i;
      cl->lvec[_j].j=cl->lvec[cl->size].j;
      cl->lvec[_j].k=cl->lvec[cl->size].k;
      _j--;
      }
  if (cl->size!=_t)
    {
    if (!cl->size) { free(cl); free(d); goto LABEL_USER_ERROR; }
    else goto STEP_2; //Reassign samples
    }
  //Step 4. Compute new clusters centers Zj=1/Nj*SUMM[XinSj](X)   j=1...Nc
  _j=cl->size;
  while (_j--)
    {
    v[_j].i/=cl->lvec[_j].w, v[_j].j/=cl->lvec[_j].w, v[_j].k/=cl->lvec[_j].w;
    _d=v[_j].i, v[_j].i=cl->lvec[_j].i, cl->lvec[_j].i=_d;
    _d=v[_j].j, v[_j].j=cl->lvec[_j].j, cl->lvec[_j].j=_d;
    _d=v[_j].k, v[_j].k=cl->lvec[_j].k, cl->lvec[_j].k=_d;
    }
  //Step 5. Compute the average distance in cluster domain from their center  Dj=1/Nj*SUMM[XinSj](||X-Zj||)   j=1...Nc
  //Step 6. Compute average distance of the samples from their centers Do=1/N*SUMM[1...Nc](Nj*Dj)    j=1...Nc
  _d=0.00; _j=cl->size; while (_j--) { d[_j]=(_rr=v2[_t].i+v2[_t].j+v2[_t].k)/cl->lvec[_j].w; _d+=_rr; }
  _d/=(double)n;
  //Step 7. If this is the last iteration set LUMP=0 and goto step 11. Otherwise if Nc<=size_cl/2 goto step 8. Otherwise if even-numbered itration  or if Nc>=size_cl*2 goto step 11. Oterwise continue;
       if ((nsteps)&&(cl->size<=size_cl/2)) goto STEP_8;
  else if ((nsteps)&&( (!(nsteps%2))||(cl->size>=size_cl*2) )) goto STEP_11;
  else continue;
  //Step 8. Find standart deviation vector for each samples subset    Vij=sqrt( 1/Nj*SUMM((Xik-Zij)**2)  )
  STEP_8:
  _j=cl->size; while (_j--) { v[_j].i=sqrt(v2[_j].i/cl->lvec[_j].w), v[_j].j=sqrt(v2[_j].j/cl->lvec[_j].w), v[_j].k=sqrt(v2[_j].k/cl->lvec[_j].w); }
  //Step 9.   Find the maximal deviation (joined with step 10 for efficiency)
  //Step 10. Splite cluster  If vmax[_j]>STDV &&  ( ( (Dj>Do)&&(Nj>2*lcutoff+1) ) || (Nc<=size_cl/2) ) then split Zj into Zj+ and Zj-
  _l=_j=cl->size;
  while (_j--)
    if (v[_j].i>v[_j].j)
      {
      if (v[_j].i>v[_j].k)
        {//Split along X
        if ( (v[_j].i>stdv) && ( (cl->size<=size_cl/2) || ( (d[_j]>_d) && ((unsigned int)cl->lvec[_j].w>2*lcutoff+1) ) ) )
          {//Split
          if (cl->size==cl->_size)
            {
            if ( (vp=(void*)realloc(cl,sizeof(t_cluster3D)+sizeof(t_lvec)*(cl->_size+=0xF))))
              { LABEL_REALLOC_FAIL: free(cl); free(d); free(v); free(v2); goto LABEL_MEMORY_ERROR; }
            else cl=(t_cluster3D*)vp;
            if ( (vp=(void*)realloc(d,sizeof(double)*cl->_size))) goto LABEL_REALLOC_FAIL;
            else d=(double*)vp;
            if (!(vp=(void*)realloc(v,sizeof(t_vec)*(cl->_size+1)))) goto LABEL_REALLOC_FAIL; //+1 required to proper qsort
            else v=(t_vec*)vp;
            if (!(vp=(void*)realloc(v2,sizeof(t_vec)*(cl->_size)))) goto LABEL_REALLOC_FAIL;
            else v2=(t_vec*)vp;
            }
          v[_j].i*=0.50;
          cl->lvec[cl->size].i=cl->lvec[_j].i+v[_j].i, cl->lvec[cl->size].j=cl->lvec[_j].j, cl->lvec[cl->size].k=cl->lvec[_j].k;
          cl->lvec[_j].i-=v[_j].i;
          cl->size++;
          }
        }
      else goto SPLIT_Z;
      }
    else if (v[_j].j>v[_j].k)
           {//Split along Y
           if ( (v[_j].j>stdv) && ( (cl->size<=size_cl/2) || ( (d[_j]>_d) && ((unsigned int)cl->lvec[_j].w>2*lcutoff+1) ) ) )
             {//Split
             if (cl->size==cl->_size)
               {
               if ( (vp=(void*)realloc(cl,sizeof(t_cluster3D)+sizeof(t_lvec)*(cl->_size+=0xF)))) goto LABEL_REALLOC_FAIL;
               else cl=(t_cluster3D*)vp;
               if ( (vp=(void*)realloc(d,sizeof(double)*cl->_size))) goto LABEL_REALLOC_FAIL;
               else d=(double*)vp;
               if (!(vp=(void*)realloc(v,sizeof(t_vec)*(cl->_size+1)))) goto LABEL_REALLOC_FAIL; //+1 required to proper qsort
               else v=(t_vec*)vp;
               if (!(vp=(void*)realloc(v2,sizeof(t_vec)*(cl->_size)))) goto LABEL_REALLOC_FAIL;
               else v2=(t_vec*)vp;
               }
             v[_j].j*=0.50;
             cl->lvec[cl->size].i=cl->lvec[_j].i, cl->lvec[cl->size].j=cl->lvec[_j].j+v[_j].j, cl->lvec[cl->size].k=cl->lvec[_j].k;
             cl->lvec[_j].j-=v[_j].j;
             cl->size++;
             }
           }
         else
           {//Split along Z
           SPLIT_Z:
           if ( (v[_j].k>stdv) && ( (cl->size<=size_cl/2) || ( (d[_j]>_d) && ((unsigned int)cl->lvec[_j].w>2*lcutoff+1) ) ) )
             {//Split
             if (cl->size==cl->_size)
               {
               if ( (vp=(void*)realloc(cl,sizeof(t_cluster3D)+sizeof(t_lvec)*(cl->_size+=0xF)))) goto LABEL_REALLOC_FAIL;
               else cl=(t_cluster3D*)vp;
               if ( (vp=(void*)realloc(d,sizeof(double)*cl->_size))) goto LABEL_REALLOC_FAIL;
               else d=(double*)vp;
               if (!(vp=(void*)realloc(v,sizeof(t_vec)*(cl->_size+1)))) goto LABEL_REALLOC_FAIL; //+1 required to proper qsort
               else v=(t_vec*)vp;
               if (!(vp=(void*)realloc(v2,sizeof(t_vec)*(cl->_size)))) goto LABEL_REALLOC_FAIL;
               else v2=(t_vec*)vp;
               }
             v[_j].k*=0.50;
             cl->lvec[cl->size].i=cl->lvec[_j].i, cl->lvec[cl->size].j=cl->lvec[_j].j, cl->lvec[cl->size].k=cl->lvec[_j].k+v[_j].k;
             cl->lvec[_j].k-=v[_j].k;
             cl->size++;
             }
           }
  if (_l!=cl->size) continue; //Splits occured
  //Step 11. Compute pairwise distances between cluster centers
  //Step 12. Compare distances with lump and arrange them in ascending order. Note we gather no more than cl->_size closest clusters
  STEP_11:
  _l=0;
  _j=cl->size;
  while (_j--)
    {
    _i=_j;
    while (_i--)
      {
      if ((_d=sqrd(cl->lvec[_j].i-cl->lvec[_i].i)+sqrd(cl->lvec[_j].j-cl->lvec[_i].j)+sqrd(cl->lvec[_j].k-cl->lvec[_i].k))<lump)
        {
        if (_l<cl->_size)
          { v[_l].i=(double)_i, v[_l].j=(double)_j, v[_l].k=_d, _l++;	} //Just add new member
        else //replace biggest
          {
          _t=0;
          _k=cl->_size;
          while (--_k)
            if (v[_t].k<v[_k].k) _t=_k;
          }
        if (v[_t].k>_d) { v[_t].i=(double)_i, v[_t].j=(double)_j, v[_t].k=_d; } //replace it
        }
      }
    }
  qsort(v,(size_t)_l,sizeof(t_vec),compare_cluster_3Dpair);
  //Step 13. Merge clusters
  _t=_l;
  for (_k=0;_k<_l;_k++)
    {
    //Merge
    _i=(unsigned int)v[_k].i;
    _j=(unsigned int)v[_k].j;
    _d=(cl->lvec[_i].w+cl->lvec[_j].w);
    cl->lvec[_i].i=(cl->lvec[_i].w*cl->lvec[_i].i+cl->lvec[_j].w*cl->lvec[_j].i)/_d;
    cl->lvec[_i].j=(cl->lvec[_i].w*cl->lvec[_i].j+cl->lvec[_j].w*cl->lvec[_j].j)/_d;
    cl->lvec[_i].k=(cl->lvec[_i].w*cl->lvec[_i].k+cl->lvec[_j].w*cl->lvec[_j].k)/_d;
    cl->lvec[_i].w=_d;
    //Del _j-th cluster
    if (_j!=--cl->size)
      {
      cl->lvec[_j].i=cl->lvec[cl->size].i, cl->lvec[_j].j=cl->lvec[cl->size].j, cl->lvec[_j].k=cl->lvec[cl->size].k, cl->lvec[_j].w=cl->lvec[cl->size].w;
      for (_i=_k+1;_i<_l;_i--)
        if ( (v[_k].i==v[_i].i)||(v[_k].i==v[_i].j)||(v[_k].j==v[_i].i)||(v[_k].j==v[_i].j) )
          {
          if (--_l!=_i) { v[_i].i=v[_l].i, v[_i].j=v[_l].k, v[_i].k=v[_l].k, _i--; }
          }
        else if ((unsigned int)v[_i].i==cl->size) v[_i].i=v[_k].j;
        else if ((unsigned int)v[_i].j==cl->size) v[_i].j=v[_k].j;
      }
    }
  if (_t!=_l) goto STEP_2;
  //Step 14. if last terminal other-wise continue
  }

//Exit
free(d);
free(v);
return cl;
}




int compare_edgef(const void *e1,const void *e2)
{
if (((t_edgef*)e1)->d<((t_edgef*)e2)->d) return  -1;
else return (int)(((t_edgef*)e1)->d!=((t_edgef*)e2)->d);
}
//Warning. It make inverse order of items
int inv_compare_edgef(const void *e1,const void *e2)
{
if (((t_edgef*)e1)->d>((t_edgef*)e2)->d) return  -1;
else return (int)(((t_edgef*)e1)->d!=((t_edgef*)e2)->d);
}

double calc_disimilarity_vec(const void *e1,const void *e2)
{
return sqrd(((t_vec*)e1)->i-((t_vec*)e2)->i)+sqrd(((t_vec*)e1)->j-((t_vec*)e2)->j)+sqrd(((t_vec*)e1)->k-((t_vec*)e2)->k);
}

double calc_disimilarity_lvec(const void *e1,const void *e2)
{
return sqrd(((t_lvec*)e1)->i-((t_lvec*)e2)->i)+sqrd(((t_lvec*)e1)->j-((t_lvec*)e2)->j)+sqrd(((t_lvec*)e1)->k-((t_lvec*)e2)->k);
}

//This function performs 3D clustering accordingly to complete-link clustering algorithm
//Note. size_cl isan optinal pointer
//Note2. sort flag makes it longer by 2*log(m) times [m - is number of clusters] but output clusters are sorted by their size
unsigned int cluster_complete_link(char sort,unsigned int *id_cl,unsigned int *sizes_cl,double cutoff,unsigned int n,size_t m,void *x,double (*calc_disimilarity)(const void *e1,const void *e2))
{
t_edgef *edgef;
unsigned int *_sizes_cl,size_cl;
unsigned int _i,_j,_k,_l;
double _d;

//Alloc memory
if (!(edgef=(t_edgef*)malloc(sizeof(t_edgef)*((size_t)n*((size_t)n-1)/2+1)))) {           MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
if   ((sizes_cl)) _sizes_cl=sizes_cl;
else if (!(_sizes_cl=(unsigned int*)malloc(sizeof(unsigned int)*(size_t)n))) { free(edgef); goto MEMORY_ERROR; }

//Calculate and sort edges
for (_d=0,_l=0,_i=0;_i<n;_i++)
  for (_j=_i+1;_j<n;_j++)
    {
    edgef[_l].vertice[0]=_i;
    edgef[_l].vertice[1]=_j;
    edgef[_l].d=calc_disimilarity(x+m*(size_t)_i,x+m*(size_t)_j);
    _d+=edgef[_l].d;
    _l++;
    }
_d/=(double)_l;
qsort(edgef,(size_t)_l,sizeof(t_edgef),inv_compare_edgef); //Algorithm expecting sorting from bigger to lower disimilarities

//Gather clusters
size_cl=n; _i=n; while(_i--) { id_cl[_i]=_i; _sizes_cl[_i]=1; }
_l=n*(n-1)/2; cutoff*=_d;
while ( (_l--)&&(edgef[_l].d<=cutoff) )
  {
  //try to merge clusters edges[_l].vertice[0] and cluster[_l].vertice[1]
  _k=_l;
  while (_k--)
    if ( ( (id_cl[edgef[_k].vertice[0]]==id_cl[edgef[_l].vertice[0]])&&(id_cl[edgef[_k].vertice[1]]==id_cl[edgef[_l].vertice[1]]) )||
         ( (id_cl[edgef[_k].vertice[0]]==id_cl[edgef[_l].vertice[1]])&&(id_cl[edgef[_k].vertice[1]]==id_cl[edgef[_l].vertice[0]]) ) )
      goto NEXT_L; //Bigger pair exists
  //Merge two clusters
  size_cl--;
  if (id_cl[edgef[_l].vertice[0]]>id_cl[edgef[_l].vertice[1]]) { _i=id_cl[edgef[_l].vertice[1]]; _j=id_cl[edgef[_l].vertice[0]]; }
  else                                                         { _i=id_cl[edgef[_l].vertice[0]]; _j=id_cl[edgef[_l].vertice[1]]; }
  _sizes_cl[_i]+=_sizes_cl[_j];
  if (_j==size_cl) { _k=n; while (_k--) if (id_cl[_k]==_j) id_cl[_k]=_i; }
  else             { _k=n; while (_k--) if (id_cl[_k]==_j) id_cl[_k]=_i; else if (id_cl[_k]==size_cl) id_cl[_k]=_j; _sizes_cl[_j]=_sizes_cl[size_cl];  }
  _sizes_cl[size_cl]=0;
  NEXT_L: ;
  }

//Sort clusters if asked
if (sort)
  {
  _k=size_cl; while (_k--) { edgef[_k].vertice[0]=_k; edgef[_k].d=(double)_sizes_cl[_k]; }
  qsort(edgef,(size_t)size_cl,sizeof(t_edgef),inv_compare_edgef);
  _k=size_cl; while (_k--) { edgef[_k].vertice[1]=_k; _sizes_cl[_k]=edgef[_k].d; edgef[_k].d=(double)edgef[_k].vertice[0]; }
  qsort(edgef,(size_t)size_cl,sizeof(t_edgef),compare_edgef);
  _k=n; while (_k--) id_cl[_k]=edgef[id_cl[_k]].vertice[1];
  }
  
//Free memory and exit
free(edgef);
if (!sizes_cl) free(_sizes_cl);
return size_cl;
}




