#define Y_CLUSTER 0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif

#ifndef Y_FFSYS
#include "y_ffsys.h"
#endif

typedef struct{
              unsigned int  size;  //Number of centers
              unsigned int _size;  //Number of reserved centers
              t_lvec       *lvec;  //Clusters centers and their 'mass'
              }t_cluster3D;

typedef struct{
              unsigned int  size;  //Number of centers
              unsigned int _size;  //Number of reserved centers
              unsigned int     m;  //Dimmensions of space
              double         **x;  //Clusters centers and their 'mass'
              }t_cluster;

typedef struct{
              unsigned int vertice[2];
              double d;
              }t_edgef;

#define CLUSTER_HPH  0
#define CLUSTER_HBA -1
#define CLUSTER_HBD +1
#define CLUSTER_NONE 128

//This function performs clustering of set of points with k-mean cluster algorithm (ISOCLUS)
//Note _t is of size; clust, D, N and Vmax are of nclus*2+1, Dij is [nclus*2+1][nclus*2+1], V is [nclus*2+1][ndim].
unsigned int cluster_kmean(unsigned int ndim,unsigned int size,double *r,unsigned int niter,unsigned int nclus,double *clus,
                           unsigned int samprm,double stdv,double lump,unsigned int maxpair,unsigned int *_t,unsigned int *N,double *D,double **Dij,double *V,unsigned int *Vmax);

//This function clusters the grid
t_cluster3D *cluster_kmean_dgrid(unsigned int nsteps,unsigned int size_cl,unsigned int lcutoff,double stdv,double lump,t_len *len,double ***x);

//This function do clustering with ISOCLUS algorithm (variaton of unweighted k-means )
//Note. t is a working massive of fitted ints.
t_cluster3D *cluster_kmean_grid(unsigned int nsteps,unsigned int size_cl,unsigned int lcutoff,double stdv,double lump,t_len *len,char ***x);

//This function do clustering with ISOCLUS algorithm (variaton of unweighted k-means )
//Note. t is a working massive of fitted ints.
t_cluster3D *cluster_kmean_3D(unsigned int nsteps,unsigned int size_cl,unsigned int lcutoff,double stdv,double lump,unsigned int n,t_vec *x);

double calc_disimilarity_vec(const void *e1,const void *e2);
double calc_disimilarity_lvec(const void *e1,const void *e2);

//This function performs 3D clustering accordingly to complete-link clustering algorithm
//Note. size_cl isan optinal pointer
//Note2. sort flag makes it longer by 2*log(m) times [m - is number of clusters] but output clusters are sorted by their size
unsigned int cluster_complete_link(char sort,unsigned int *id_cl,unsigned int *sizes_cl,double cutoff,unsigned int n,size_t m,void *x,double (*calc_disimilarity)(const void *e1,const void *e2));

