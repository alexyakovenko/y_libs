
#include "y_solvent.h"

//-------------------------------------------------- D T R    F I L L E R S --------------------------------------------------------//

//This function allocate memory for t_wat_db
//NOTE. t_wat_db iw simply destroyed with free(*it)
t_wat_db *alloc_wat_db(unsigned int n_cfgs,unsigned int _n_E,unsigned int _n_r)
{
t_wat_db *w_db;
if (!(w_db=(t_wat_db*)malloc(sizeof(t_wat_db)+sizeof(t_wat_cfg)*n_cfgs+sizeof(double)*_n_E+sizeof(t_vec)*(_n_E+_n_r)))) 
  { ylib_errno=YERROR_MEMORY; return FALSE; }
else { w_db->n_cfgs=n_cfgs, w_db->_n_E=_n_E, w_db->_n_r=_n_r; }
w_db->cfgs=(void*)w_db+sizeof(t_wat_db);
w_db->_e=(void*)w_db+sizeof(t_wat_db)+sizeof(t_wat_cfg)*n_cfgs;
w_db->_E=(void*)w_db+sizeof(t_wat_db)+sizeof(t_wat_cfg)*n_cfgs+_n_E* sizeof(double);
w_db->_r=(void*)w_db+sizeof(t_wat_db)+sizeof(t_wat_cfg)*n_cfgs+_n_E*(sizeof(double)+sizeof(t_vec));

return w_db;
}

//This function reads wat_db from file
//DB file structure
// [ [ n_cfgs ], [ _n_E ], [ _n_r ] ]
// [ [ n_E ], [ n_r ], [ th ], [ d2 ] ] x n_cfgs
// [ _E ] x _n_E, [ _r ] x _n_r
t_wat_db *read_wat_db(FILE *in)
{
unsigned int n_cfgs, _n_E, _n_r;
double qO, qH;
t_wat_db *w_db;

//Stage I. Allocate structure
if ( (fread(&n_cfgs,sizeof(unsigned int),0x1,in)!=0x1)||(n_cfgs!=Y_MAGIC)||(fread(&qO,sizeof(double),0x1,in)!=0x1)||(fread(&qH,sizeof(double),0x1,in)!=0x1)||
     (fread(&n_cfgs,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&_n_E,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&_n_r,sizeof(unsigned int),0x1,in)!=0x1)    ) 
  { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if (!(w_db=alloc_wat_db(n_cfgs,_n_E,_n_r))) return FALSE;
else { w_db->qO=qO, w_db->qH=qH; }
//Stage II. Read data into wat_db
for (w_db->n_cfgs=w_db->_n_E=w_db->_n_r=0; w_db->n_cfgs<n_cfgs; w_db->n_cfgs++)
  {
  w_db->cfgs[w_db->n_cfgs].e=&w_db->_e[w_db->_n_E];
  w_db->cfgs[w_db->n_cfgs].E=&w_db->_E[w_db->_n_E];
  w_db->cfgs[w_db->n_cfgs].r=&w_db->_r[w_db->_n_r];
  if ( (fread(&w_db->cfgs[w_db->n_cfgs].n_E,sizeof(unsigned int),0x1,in)!=0x1)||((w_db->_n_E+=w_db->cfgs[w_db->n_cfgs].n_E)>_n_E)||
       (fread(&w_db->cfgs[w_db->n_cfgs].n_r,sizeof(unsigned int),0x1,in)!=0x1)||((w_db->_n_r+=w_db->cfgs[w_db->n_cfgs].n_E*w_db->cfgs[w_db->n_cfgs].n_r)>_n_r)||
       (fread(w_db->cfgs[w_db->n_cfgs].th,sizeof(t_vec),0x4,in)!=0x4) )
    { LABEL_LOAD_FAILURE: free(w_db); goto LABEL_IO_ERROR; }
  }
if ( (fread(w_db->_e,sizeof(double),w_db->_n_E,in)!=_n_E)||(fread(w_db->_E,sizeof(t_vec),w_db->_n_E,in)!=_n_E)||(fread(w_db->_r,sizeof(t_vec),w_db->_n_r,in)!=_n_r) )
  goto LABEL_LOAD_FAILURE; 
//Everything is fine: fill d2 and exit
while (n_cfgs--)
  {
  w_db->cfgs[n_cfgs].d2[0][1]=w_db->cfgs[n_cfgs].d2[1][0]=sqrt(calc_distance(&w_db->cfgs[n_cfgs].th[0],&w_db->cfgs[n_cfgs].th[1]));
  w_db->cfgs[n_cfgs].d2[0][2]=w_db->cfgs[n_cfgs].d2[2][0]=sqrt(calc_distance(&w_db->cfgs[n_cfgs].th[0],&w_db->cfgs[n_cfgs].th[2]));
  w_db->cfgs[n_cfgs].d2[0][3]=w_db->cfgs[n_cfgs].d2[3][0]=sqrt(calc_distance(&w_db->cfgs[n_cfgs].th[0],&w_db->cfgs[n_cfgs].th[3]));
  w_db->cfgs[n_cfgs].d2[1][2]=w_db->cfgs[n_cfgs].d2[2][1]=sqrt(calc_distance(&w_db->cfgs[n_cfgs].th[1],&w_db->cfgs[n_cfgs].th[2]));
  w_db->cfgs[n_cfgs].d2[1][3]=w_db->cfgs[n_cfgs].d2[3][1]=sqrt(calc_distance(&w_db->cfgs[n_cfgs].th[1],&w_db->cfgs[n_cfgs].th[3]));
  w_db->cfgs[n_cfgs].d2[2][3]=w_db->cfgs[n_cfgs].d2[3][2]=sqrt(calc_distance(&w_db->cfgs[n_cfgs].th[2],&w_db->cfgs[n_cfgs].th[3]));
  }
//Check if there are empty or overcrouded configurations in db
n_cfgs=w_db->n_cfgs;
while (n_cfgs--)
  if ( (w_db->cfgs[n_cfgs].n_r<3)||(w_db->cfgs[n_cfgs].n_r>3*YSOLVENT_MAX_WAT_MOLS_IN_DB)||( (w_db->cfgs[n_cfgs].n_r%3)) ) { free(w_db); w_db=0x0; return FALSE; }
//Anything is cool, exit
return w_db;
}
//This function write wat_db to HDD 
char write_wat_db(FILE *out,t_wat_db *w_db)
{
unsigned int i=Y_MAGIC;

if ( (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&w_db->qO,sizeof(double),0x1,out)!=0x1)||(fwrite(&w_db->qH,sizeof(double),0x1,out)!=0x1)||
     (fwrite(&w_db->n_cfgs,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&w_db->_n_E,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&w_db->_n_r,sizeof(unsigned int),0x1,out)!=0x1) )
  { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
for (i=0; i<w_db->n_cfgs; i++)
  if ( (fwrite(&w_db->cfgs[i].n_E,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&w_db->cfgs[i].n_r,sizeof(unsigned int),0x1,out)!=0x1)||
       (fwrite(w_db->cfgs[i].th,sizeof(t_vec),0x4,out)!=0x4) ) goto LABEL_IO_ERROR;
if ( (fwrite(w_db->_e,sizeof(double),w_db->_n_E,out)!=w_db->_n_E)||
     (fwrite(w_db->_E,sizeof(t_vec),w_db->_n_E,out)!=w_db->_n_E) ||
     (fwrite(w_db->_r,sizeof(t_vec),w_db->_n_r,out)!=w_db->_n_r)  ) goto LABEL_IO_ERROR; 
return TRUE;
}

//This function imports y_wat.db from it's standard location at ../top/
t_wat_db *import_wat_db() 
{
t_wat_db *w_db;
FILE *in;
if (!(in=fopen("../top/y_wat.db","r"))) { ylib_errno=YERROR_USER; return FALSE; }
if (!(w_db=read_wat_db(in))) { fclose(in); return FALSE;}
else                         { fclose(in); return w_db; }
}

//------------------------------------ F I L L E R S ----------------------------------------------------------//


//This function fill a union of ths with a single water
unsigned int fill_th_union(t_vec *rvecs,t_vec *th_cm,t_vec *E,t_wat_db *w_db)
{
register unsigned int _i, _j, _k;
register double _d, __d;

//Stage 0. Confirm if the recond with only a water molecule is zero
for (_i=0;_i<w_db->n_cfgs;_i++) if (w_db->cfgs[_i].n_r==3) break; 
if (_i==w_db->n_cfgs) { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; } //solvation without a records is not implemented
//Stage I. Find the closest E vector
_d=calc_distance(E,&w_db->cfgs[_i].E[0]), _j=0, _k=w_db->cfgs[_i].n_E; 
while (--_k)
  if ((__d=calc_distance(E,&w_db->cfgs[_i].E[_k]))<_d) { _d=__d, _j=_k; }  
//Stage II. Copy (with direct transformation) water molecules coordinates into a structure 
rvecs[0].i=th_cm->i+w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+0].i,
rvecs[0].j=th_cm->j+w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+0].j,
rvecs[0].k=th_cm->k+w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+0].k;
rvecs[1].i=th_cm->i+w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+1].i,
rvecs[1].j=th_cm->j+w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+1].j,
rvecs[1].k=th_cm->k+w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+1].k;
rvecs[2].i=th_cm->i+w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+2].i,
rvecs[2].j=th_cm->j+w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+2].j,
rvecs[2].k=th_cm->k+w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+2].k;
//Stage IV. Exit.
return 1;
}

//This function fills a thetrahedron with water using empirical DB
//NOTE. th define rvecs, r stores results (should be of MAX_WATS reserved) and returns amount of added molecules.
//NB! thetrahedrons should have prechecked volume, here we just find appropriate fillings.
unsigned int fill_th(t_vec *rvecs,unsigned int th[4],t_lvec *lvecs,t_vec *E,t_wat_db *w_db)
{
unsigned int _i, _j, _k, _l;
unsigned int _th[4];
double _d, __d, ___d, d2[4][4];
t_vec tr, mx, my;
t_tensor T, U, V;

//Stage I. Align tetrahedron.
//Order vertices of tetrahedron so the longest edge is _th[0]..._th[1]
_th[0]=th[0], _th[1]=th[1], _th[2]=th[2], _th[3]=th[3];
     d2[1][0]=d2[0][1]=sqrt(fabs(calc_distance((t_vec*)&lvecs[_th[0]],(t_vec*)&lvecs[_th[1]]))),               _i=0, _j=1;
if ((d2[2][0]=d2[0][2]=sqrt(fabs(calc_distance((t_vec*)&lvecs[_th[0]],(t_vec*)&lvecs[_th[2]]))))>d2[_i][_j]) { _i=0, _j=2; }
if ((d2[3][0]=d2[0][3]=sqrt(fabs(calc_distance((t_vec*)&lvecs[_th[0]],(t_vec*)&lvecs[_th[3]]))))>d2[_i][_j]) { _i=0, _j=3; }
if ((d2[2][1]=d2[1][2]=sqrt(fabs(calc_distance((t_vec*)&lvecs[_th[1]],(t_vec*)&lvecs[_th[2]]))))>d2[_i][_j]) { _i=3, _j=1; }
if ((d2[3][1]=d2[1][3]=sqrt(fabs(calc_distance((t_vec*)&lvecs[_th[1]],(t_vec*)&lvecs[_th[3]]))))>d2[_i][_j]) { _i=3, _j=1; }
if ((d2[3][2]=d2[2][3]=sqrt(fabs(calc_distance((t_vec*)&lvecs[_th[2]],(t_vec*)&lvecs[_th[3]]))))>d2[_i][_j]) { _i=2, _j=3; }
if (_i!=0) //_i!=1 by construction
  {
  if (_i==2)
    {
    _l=_th[0], _th[0]=_th[2], _th[2]=_l;
    _d=d2[1][0], d2[0][1]=d2[1][0]=d2[1][2], d2[2][1]=d2[1][2]=_d;
    _d=d2[3][0], d2[0][3]=d2[3][0]=d2[3][2], d2[2][3]=d2[3][2]=_d;
    }
  else
    {
    _l=_th[0], _th[0]=_th[3], _th[3]=_l;
    _d=d2[1][0], d2[0][1]=d2[1][0]=d2[1][3], d2[3][1]=d2[1][3]=_d;
    _d=d2[2][0], d2[0][2]=d2[2][0]=d2[2][3], d2[3][2]=d2[2][3]=_d;
    }
  }
if (_j!=1) //_j!=0 by construction
  {
  if (_j==2)
    {
    _l=_th[1], _th[1]=_th[2], _th[2]=_l;
    _d=d2[0][1], d2[1][0]=d2[0][1]=d2[0][2], d2[2][0]=d2[0][2]=_d;
    _d=d2[3][1], d2[1][3]=d2[3][1]=d2[3][2], d2[2][3]=d2[3][2]=_d;
    }
  else
    {
    _l=_th[1], _th[1]=_th[3], _th[3]=_l;
    _d=d2[0][1], d2[1][0]=d2[0][1]=d2[0][3], d2[3][0]=d2[0][3]=_d;
    _d=d2[2][1], d2[1][2]=d2[2][1]=d2[2][3], d2[3][2]=d2[2][3]=_d;
    }
  }
//Find the most appropriate thetrahedron
_i=0, _j=0, _d=sqrd(w_db->cfgs[0].d2[0][1]-d2[0][1])+sqrd(w_db->cfgs[0].d2[0][2]-d2[0][2])+sqrd(w_db->cfgs[0].d2[0][3]-d2[0][3])+
               sqrd(w_db->cfgs[0].d2[1][2]-d2[1][2])+sqrd(w_db->cfgs[0].d2[1][3]-d2[1][3])+sqrd(w_db->cfgs[0].d2[2][3]-d2[2][3]);
_k=w_db->n_cfgs;
while (_k--)
  {
  ___d=sqrd(w_db->cfgs[_k].d2[0][1]-d2[0][1])+sqrd(w_db->cfgs[_k].d2[2][3]-d2[2][3]);
  if ((__d=sqrd(w_db->cfgs[_k].d2[0][2]-d2[0][2])+sqrd(w_db->cfgs[_k].d2[0][3]-d2[0][3])+
           sqrd(w_db->cfgs[_k].d2[1][2]-d2[1][2])+sqrd(w_db->cfgs[_k].d2[1][3]-d2[1][3])+___d)<_d) { _d=__d, _i=_k, _j=0; }
  if ((__d=sqrd(w_db->cfgs[_k].d2[0][2]-d2[0][3])+sqrd(w_db->cfgs[_k].d2[0][3]-d2[0][2])+
           sqrd(w_db->cfgs[_k].d2[1][2]-d2[1][3])+sqrd(w_db->cfgs[_k].d2[1][3]-d2[1][2])+___d)<_d) { _d=__d, _i=_k, _j=1; } 
  if ((__d=sqrd(w_db->cfgs[_k].d2[0][2]-d2[1][2])+sqrd(w_db->cfgs[_k].d2[0][3]-d2[1][3])+
           sqrd(w_db->cfgs[_k].d2[1][2]-d2[0][2])+sqrd(w_db->cfgs[_k].d2[1][3]-d2[0][3])+___d)<_d) { _d=__d, _i=_k, _j=2; } 
  if ((__d=sqrd(w_db->cfgs[_k].d2[0][2]-d2[1][3])+sqrd(w_db->cfgs[_k].d2[0][3]-d2[1][2])+
           sqrd(w_db->cfgs[_k].d2[1][2]-d2[0][3])+sqrd(w_db->cfgs[_k].d2[1][3]-d2[0][2])+___d)<_d) { _d=__d, _i=_k, _j=3; } 
  }
if ( (_j)) //Align _th, d2 matrix is needed no more
  { 
       if (_j==1) { _l=_th[2], _th[2]=_th[3], _th[3]=_l; } 
  else if (_j==2) { _l=_th[0], _th[0]=_th[1], _th[1]=_l; }
  else            { _l=_th[0], _th[0]=_th[1], _th[1]=_l; _l=_th[2], _th[2]=_th[3], _th[3]=_l; }
  }

//Stage II. Find optimal RT transformation
//Th aligned, find optimal RT transformation
mx.i=w_db->cfgs[_i].th[0].i+w_db->cfgs[_i].th[1].i+w_db->cfgs[_i].th[2].i+w_db->cfgs[_i].th[3].i,
mx.j=w_db->cfgs[_i].th[0].j+w_db->cfgs[_i].th[1].j+w_db->cfgs[_i].th[2].j+w_db->cfgs[_i].th[3].j,
mx.k=w_db->cfgs[_i].th[0].k+w_db->cfgs[_i].th[1].k+w_db->cfgs[_i].th[2].k+w_db->cfgs[_i].th[3].k;
mx.i/=4., mx.j/=4., mx.k/=4.;
my.i=lvecs[_th[0]].i+lvecs[_th[1]].i+lvecs[_th[2]].i+lvecs[_th[3]].i,
my.j=lvecs[_th[0]].j+lvecs[_th[1]].j+lvecs[_th[2]].j+lvecs[_th[3]].j,
my.k=lvecs[_th[0]].k+lvecs[_th[1]].k+lvecs[_th[2]].k+lvecs[_th[3]].k;
my.i/=4., my.j/=4., my.k/=4.;
T[0][0]=T[0][1]=T[0][2]=T[1][0]=T[1][1]=T[1][2]=T[2][0]=T[2][1]=T[2][2]=0.;

T[0][0]+=(lvecs[_th[0]].i-my.i)*(w_db->cfgs[_i].th[0].i-mx.i),
T[0][1]+=(lvecs[_th[0]].i-my.i)*(w_db->cfgs[_i].th[0].j-mx.j),
T[0][2]+=(lvecs[_th[0]].i-my.i)*(w_db->cfgs[_i].th[0].k-mx.k),
T[1][0]+=(lvecs[_th[0]].j-my.j)*(w_db->cfgs[_i].th[0].i-mx.i),
T[1][1]+=(lvecs[_th[0]].j-my.j)*(w_db->cfgs[_i].th[0].j-mx.j),
T[1][2]+=(lvecs[_th[0]].j-my.j)*(w_db->cfgs[_i].th[0].k-mx.k),
T[2][0]+=(lvecs[_th[0]].k-my.k)*(w_db->cfgs[_i].th[0].i-mx.i),
T[2][1]+=(lvecs[_th[0]].k-my.k)*(w_db->cfgs[_i].th[0].j-mx.j),
T[2][2]+=(lvecs[_th[0]].k-my.k)*(w_db->cfgs[_i].th[0].k-mx.k);

T[0][0]+=(lvecs[_th[1]].i-my.i)*(w_db->cfgs[_i].th[1].i-mx.i),
T[0][1]+=(lvecs[_th[1]].i-my.i)*(w_db->cfgs[_i].th[1].j-mx.j),
T[0][2]+=(lvecs[_th[1]].i-my.i)*(w_db->cfgs[_i].th[1].k-mx.k),
T[1][0]+=(lvecs[_th[1]].j-my.j)*(w_db->cfgs[_i].th[1].i-mx.i),
T[1][1]+=(lvecs[_th[1]].j-my.j)*(w_db->cfgs[_i].th[1].j-mx.j),
T[1][2]+=(lvecs[_th[1]].j-my.j)*(w_db->cfgs[_i].th[1].k-mx.k),
T[2][0]+=(lvecs[_th[1]].k-my.k)*(w_db->cfgs[_i].th[1].i-mx.i),
T[2][1]+=(lvecs[_th[1]].k-my.k)*(w_db->cfgs[_i].th[1].j-mx.j),
T[2][2]+=(lvecs[_th[1]].k-my.k)*(w_db->cfgs[_i].th[1].k-mx.k);

T[0][0]+=(lvecs[_th[2]].i-my.i)*(w_db->cfgs[_i].th[2].i-mx.i),
T[0][1]+=(lvecs[_th[2]].i-my.i)*(w_db->cfgs[_i].th[2].j-mx.j),
T[0][2]+=(lvecs[_th[2]].i-my.i)*(w_db->cfgs[_i].th[2].k-mx.k),
T[1][0]+=(lvecs[_th[2]].j-my.j)*(w_db->cfgs[_i].th[2].i-mx.i),
T[1][1]+=(lvecs[_th[2]].j-my.j)*(w_db->cfgs[_i].th[2].j-mx.j),
T[1][2]+=(lvecs[_th[2]].j-my.j)*(w_db->cfgs[_i].th[2].k-mx.k),
T[2][0]+=(lvecs[_th[2]].k-my.k)*(w_db->cfgs[_i].th[2].i-mx.i),
T[2][1]+=(lvecs[_th[2]].k-my.k)*(w_db->cfgs[_i].th[2].j-mx.j),
T[2][2]+=(lvecs[_th[2]].k-my.k)*(w_db->cfgs[_i].th[2].k-mx.k);

T[0][0]+=(lvecs[_th[3]].i-my.i)*(w_db->cfgs[_i].th[3].i-mx.i),
T[0][1]+=(lvecs[_th[3]].i-my.i)*(w_db->cfgs[_i].th[3].j-mx.j),
T[0][2]+=(lvecs[_th[3]].i-my.i)*(w_db->cfgs[_i].th[3].k-mx.k),
T[1][0]+=(lvecs[_th[3]].j-my.j)*(w_db->cfgs[_i].th[3].i-mx.i),
T[1][1]+=(lvecs[_th[3]].j-my.j)*(w_db->cfgs[_i].th[3].j-mx.j),
T[1][2]+=(lvecs[_th[3]].j-my.j)*(w_db->cfgs[_i].th[3].k-mx.k),
T[2][0]+=(lvecs[_th[3]].k-my.k)*(w_db->cfgs[_i].th[3].i-mx.i),
T[2][1]+=(lvecs[_th[3]].k-my.k)*(w_db->cfgs[_i].th[3].j-mx.j),
T[2][2]+=(lvecs[_th[3]].k-my.k)*(w_db->cfgs[_i].th[3].k-mx.k);

T[0][0]/=4., T[0][1]/=4., T[0][2]/=4., T[1][0]/=4., T[1][1]/=4., T[1][2]/=4., T[2][0]/=4., T[2][1]/=4., T[2][2]/=4.;
 
//Calc svd
tensor_svd(TRUE,TRUE,TRUE,SMALL2,&T,&U,&tr,&V);
if (fabs(1.-fabs(TENSOR_DET(U)))>SMALL2) round_udet_tensor(&U,SMALL2);
if (fabs(1.-fabs(TENSOR_DET(V)))>SMALL2) round_udet_tensor(&V,SMALL2);
if (TENSOR_DET(U)*TENSOR_DET(V)<0.) { U[0][2]=-U[0][2], U[1][2]=-U[1][2], U[2][2]=-U[2][2]; }
multiple_origin_tensor_origin_tensor(&T,&U,&V); //Calc T=U.S.V^T
//Calc b=my-R.mx
multiple_origin_tensor_origin_vec(&tr,&T,&mx);
tr.i=my.i-tr.i, tr.j=my.j-tr.j, tr.k=my.k-tr.k;

//Stage III. Compute desired configuration by external el. field
//Find E with inverse transformation
multiple_transp_tensor_origin_vec(&my,&T,E);
//Find the closest E vector
_d=calc_distance(&my,&w_db->cfgs[_i].E[0]), _j=0, _k=w_db->cfgs[_i].n_E; while (--_k) if ((__d=calc_distance(&my,&w_db->cfgs[_i].E[_k]))<_d) { _d=__d, _j=_k; }  

//Stage IV. Copy (with direct transformation) water molecules coordinates into a structure 
_k=w_db->cfgs[_i].n_r;
while (_k--) 
  {
  multiple_transp_tensor_origin_vec(&rvecs[_k],&T,&w_db->cfgs[_i].r[w_db->cfgs[_i].n_r*_j+_k]); 
  rvecs[_k].i+=tr.i, rvecs[_k].j+=tr.j, rvecs[_k].k+=tr.k;
  }

//Stage V. Exit.
return w_db->cfgs[_i].n_r/3;
}

//This function solvates complex explicitely
//It solvates thetragedrons sorted by their E. theragedrons are filled using db or by gouping witha common facet or common edge of sufficient lenght.
//NOTE. This function does max_nwats+/-biggest ensable from db
char solvate_explicit(unsigned int *nwats,unsigned int max_nwats,double min_th_ir,double min_th_vol,double max_th_vol,t_list *set,unsigned int n,t_vec *rvecs,double *q,t_delaunay *dtr,double *radii,t_wat_db *w_db)
{
unsigned int *th_id, nwat; 
t_vec *th_cm, *E, _r, p;
double *ee, *vol;
char *th_st;

unsigned int _i, _j, _k, _l;
double _ir, _d;

//Stage 0. Init thx
*nwats=0;
if (!(th_id=(unsigned int*)malloc((sizeof(unsigned int)+(sizeof(t_vec)+sizeof(double))*2+sizeof(char))*dtr->size_th))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { E=(t_vec*)((void*)th_id+sizeof(unsigned int)*dtr->size_th); th_cm=E+dtr->size_th; ee=(double*)((void*)th_cm+sizeof(t_vec)*dtr->size_th); vol=ee+dtr->size_th, th_st=(char*)((void*)vol+sizeof(double)*dtr->size_th); }

//Stage I.1 Define borders and the SET thetrahedrons
_ir=2.*min_th_ir, _l=0, memset(th_st,0x0,sizeof(char)*dtr->size); 
_j=0, _i=dtr->size_th;
while (_i--)
  if ( (dtr->th[_i][0]<4)||(dtr->th[_i][1]<4)||(dtr->th[_i][2]<4)||(dtr->th[_i][3]<4) ) th_st[_i]=(unsigned char)-126;
  else if ( (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][1]])<sqrd(_ir+radii[dtr->th[_i][0]]+radii[dtr->th[_i][1]]))&&
            (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][2]])<sqrd(_ir+radii[dtr->th[_i][0]]+radii[dtr->th[_i][2]]))&&
            (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][3]])<sqrd(_ir+radii[dtr->th[_i][0]]+radii[dtr->th[_i][3]]))&&
            (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][2]],(t_vec*)&dtr->lvec[dtr->th[_i][3]])<sqrd(_ir+radii[dtr->th[_i][2]]+radii[dtr->th[_i][3]]))&&
            (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][2]])<sqrd(_ir+radii[dtr->th[_i][1]]+radii[dtr->th[_i][2]]))&&
            (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][3]])<sqrd(_ir+radii[dtr->th[_i][1]]+radii[dtr->th[_i][3]])) )  
         th_st[_i]=(unsigned char)-127;
  else if ( ( (find_in_list(FALSE,dtr->th[_i][0],set)))||( (find_in_list(FALSE,dtr->th[_i][1],set)))||
            ( (find_in_list(FALSE,dtr->th[_i][2],set)))||( (find_in_list(FALSE,dtr->th[_i][3],set))) ) { th_id[_l++]=_i, th_st[_i]=(unsigned char)-1; }

//Stage I.2. Define 'local' thetrahedrons
for (_j=0; _j<_l; _j++)
  {
  _i=th_id[_j];
  if ( (!th_st[dtr->nth[_i][0]])&&
     ( (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][2]])>sqrd(_ir+radii[dtr->th[_i][1]]+radii[dtr->th[_i][2]]))&&
       (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][3]])>sqrd(_ir+radii[dtr->th[_i][1]]+radii[dtr->th[_i][3]]))&&
       (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][2]],(t_vec*)&dtr->lvec[dtr->th[_i][3]])>sqrd(_ir+radii[dtr->th[_i][2]]+radii[dtr->th[_i][3]])) ) )
    {
    _k=dtr->nth[_i][0];
    if ( ((unsigned char)th_st[dtr->nth[_k][0]]!=(unsigned char)-126)&&((unsigned char)th_st[dtr->nth[_k][1]]!=(unsigned char)-126)&&
         ((unsigned char)th_st[dtr->nth[_k][2]]!=(unsigned char)-126)&&((unsigned char)th_st[dtr->nth[_k][3]]!=(unsigned char)-126) )
      { th_id[_l++]=_k, th_st[_k]=(unsigned char)-1; }
    }
  if ( (!th_st[dtr->nth[_i][1]])&&
     ( (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][2]])>sqrd(_ir+radii[dtr->th[_i][0]]+radii[dtr->th[_i][2]]))&&
       (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][3]])>sqrd(_ir+radii[dtr->th[_i][0]]+radii[dtr->th[_i][3]]))&&
       (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][2]],(t_vec*)&dtr->lvec[dtr->th[_i][3]])>sqrd(_ir+radii[dtr->th[_i][2]]+radii[dtr->th[_i][3]])) ) )
    {
    _k=dtr->nth[_i][1];
    if ( ((unsigned char)th_st[dtr->nth[_k][0]]!=(unsigned char)-126)&&((unsigned char)th_st[dtr->nth[_k][1]]!=(unsigned char)-126)&&
         ((unsigned char)th_st[dtr->nth[_k][2]]!=(unsigned char)-126)&&((unsigned char)th_st[dtr->nth[_k][3]]!=(unsigned char)-126) )
      { th_id[_l++]=_k, th_st[_k]=(unsigned char)-1; }
    }
  if ( (!th_st[dtr->nth[_i][2]])&&
     ( (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][1]])>sqrd(_ir+radii[dtr->th[_i][0]]+radii[dtr->th[_i][1]]))&&
       (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][3]])>sqrd(_ir+radii[dtr->th[_i][0]]+radii[dtr->th[_i][3]]))&&
       (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][3]])>sqrd(_ir+radii[dtr->th[_i][1]]+radii[dtr->th[_i][3]])) ) )
    {
    _k=dtr->nth[_i][2];
    if ( ((unsigned char)th_st[dtr->nth[_k][0]]!=(unsigned char)-126)&&((unsigned char)th_st[dtr->nth[_k][1]]!=(unsigned char)-126)&&
         ((unsigned char)th_st[dtr->nth[_k][2]]!=(unsigned char)-126)&&((unsigned char)th_st[dtr->nth[_k][3]]!=(unsigned char)-126) )
      { th_id[_l++]=_k, th_st[_k]=(unsigned char)-1; }
    }
  if ( (!th_st[dtr->nth[_i][3]])&&
     ( (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][1]])>sqrd(_ir+radii[dtr->th[_i][0]]+radii[dtr->th[_i][1]]))&&
       (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][2]])>sqrd(_ir+radii[dtr->th[_i][0]]+radii[dtr->th[_i][2]]))&&
       (calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][2]])>sqrd(_ir+radii[dtr->th[_i][1]]+radii[dtr->th[_i][2]])) ) )
    {
    _k=dtr->nth[_i][3];
    if ( ((unsigned char)th_st[dtr->nth[_k][0]]!=(unsigned char)-126)&&((unsigned char)th_st[dtr->nth[_k][1]]!=(unsigned char)-126)&&
         ((unsigned char)th_st[dtr->nth[_k][2]]!=(unsigned char)-126)&&((unsigned char)th_st[dtr->nth[_k][3]]!=(unsigned char)-126) )
      { th_id[_l++]=_k, th_st[_k]=(unsigned char)-1; }
    }
  }

//Stage I.3. Reset marks of local thetrahedrons
_i=dtr->size_th; 
while (_i--) 
  if ((unsigned char)th_st[_i]==(unsigned char)-126) th_st[_i]=(unsigned char)-127;
  else if ((unsigned char)th_st[_i]!=(unsigned char)-127) 
         {
         if (!th_st[_i]) th_st[_i]=(unsigned char)-127;
         else th_st[_i]=(unsigned char)-1;
         }

//Stage II.1. Precompute volumes and cms thetrahedrons 
_i=dtr->size_th;
while(_i--)
  {
  th_id[_i]=_i;
  if ((unsigned char)th_st[_i]!=(unsigned char)-127)
    {
    vol[_i]=fabs(calc_th_volume((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][2]],(t_vec*)&dtr->lvec[dtr->th[_i][3]]));
    if (vol[_i]>max_th_vol) th_st[_i]=(unsigned char)-127;
    else 
      {
      calc_th_cm(&th_cm[_i],(t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][2]],(t_vec*)&dtr->lvec[dtr->th[_i][3]]);
      if (vol[_i]>min_th_vol) th_st[_i]=1;
      }
    }
  }

//Stage II.2. Screen unions with common facet
_i=dtr->size_th;
while(--_i)
  if ((unsigned char)th_st[_i]==(unsigned char)-1)
    {//Chose pair with the biggest r and suitable vol 
    _j=(unsigned int)-1, _ir=0.;
    if ( (dtr->nth[_i][0]<_i)&&(th_st[dtr->nth[_i][0]]==-1)&&(vol[_i]+vol[dtr->nth[_i][0]]>min_th_vol) )
      {
      _ir=calc_triangle_iradii((t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][2]],(t_vec*)&dtr->lvec[dtr->th[_i][3]]);
      _j=0; 
      }
    if ( (dtr->nth[_i][1]<_i)&&(th_st[dtr->nth[_i][1]]==-1)&&(vol[_i]+vol[dtr->nth[_i][1]]>min_th_vol) )
      {
      _d=calc_triangle_iradii((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][2]],(t_vec*)&dtr->lvec[dtr->th[_i][3]]);
      if (_d>_ir) { _j=1, _ir=_d; }
      }
    if ( (dtr->nth[_i][2]<_i)&&(th_st[dtr->nth[_i][2]]==-1)&&(vol[_i]+vol[dtr->nth[_i][2]]>min_th_vol) )
      {
      _d=calc_triangle_iradii((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][3]]);
      if (_d>_ir) { _j=2, _ir=_d; }
      }
    if ( (dtr->nth[_i][3]<_i)&&(th_st[dtr->nth[_i][3]]==-1)&&(vol[_i]+vol[dtr->nth[_i][3]]>min_th_vol) )
      {
      _d=calc_triangle_iradii((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][2]]);
      if (_d>_ir) { _j=3, _ir=_d; }
      }
    if (_ir>min_th_ir) 
      {
      th_st[_i]=0, th_st[dtr->nth[_i][_j]]=(unsigned char)-127; 
      th_cm[_i].i+=th_cm[dtr->nth[_i][_j]].i, th_cm[_i].j+=th_cm[dtr->nth[_i][_j]].j, th_cm[_i].k+=th_cm[dtr->nth[_i][_j]].k; 
      th_cm[_i].i/=2., th_cm[_i].j/=2., th_cm[_i].k/=2.;
      }
    }

//Stage II.3. Screen unions with common edge
_i=dtr->size_th;
while(_i--)
  if ((unsigned char)th_st[_i]==(unsigned char)-1)
    {
    //Chose union with the biggest r and suitable vol 
        _ir=calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][1]]),        _j=0;
    if ((_d=calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][2]]))>_ir) { _j=1, _ir=_d; }
    if ((_d=calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][0]],(t_vec*)&dtr->lvec[dtr->th[_i][3]]))>_ir) { _j=2, _ir=_d; }
    if ((_d=calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][2]]))>_ir) { _j=3, _ir=_d; }
    if ((_d=calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][1]],(t_vec*)&dtr->lvec[dtr->th[_i][3]]))>_ir) { _j=4, _ir=_d; }
    if ((_d=calc_distance((t_vec*)&dtr->lvec[dtr->th[_i][2]],(t_vec*)&dtr->lvec[dtr->th[_i][3]]))>_ir) { _j=5, _ir=_d; }
    switch (_j)
      {
      case  0 : {//No replacements 
                break; }
      case  1 : {//b <---> c
                _k=dtr->th[_i][1],  dtr->th[_i][1]=dtr->th[_i][2],    dtr->th[_i][2]=_k; 
                _k=dtr->nth[_i][1], dtr->nth[_i][1]=dtr->nth[_i][2],  dtr->nth[_i][2]=_k; 
                break; }
      case  2 : {//b <---> d
                _k=dtr->th[_i][1],  dtr->th[_i][1]=dtr->th[_i][3],    dtr->th[_i][3]=_k; 
                _k=dtr->nth[_i][1], dtr->nth[_i][1]=dtr->nth[_i][3],  dtr->nth[_i][3]=_k; 
                break; }
      case  3 : {//a <---> c
                _k=dtr->th[_i][0],  dtr->th[_i][0]=dtr->th[_i][2],    dtr->th[_i][2]=_k; 
                _k=dtr->nth[_i][0], dtr->nth[_i][0]=dtr->nth[_i][2],  dtr->nth[_i][2]=_k; 
                break; }
      case  4 : {//a <---> d
                _k=dtr->th[_i][0],  dtr->th[_i][0]=dtr->th[_i][3],    dtr->th[_i][3]=_k; 
                _k=dtr->nth[_i][0], dtr->nth[_i][0]=dtr->nth[_i][3],  dtr->nth[_i][3]=_k; 
                break; }
      case  5 : {//a <---> c && b <---> d
                _k=dtr->th[_i][0],  dtr->th[_i][0]=dtr->th[_i][2],    dtr->th[_i][2]=_k; 
                _k=dtr->nth[_i][0], dtr->nth[_i][0]=dtr->nth[_i][2],  dtr->nth[_i][2]=_k; 
                _k=dtr->th[_i][1],  dtr->th[_i][1]=dtr->th[_i][3],    dtr->th[_i][3]=_k; 
                _k=dtr->nth[_i][1], dtr->nth[_i][1]=dtr->nth[_i][3],  dtr->nth[_i][3]=_k; 
                break; }
      }
    //Join all tetrahedrons
    _j=_i, _k=3, _l=dtr->th[_i][2], _d=vol[_i];
    while (dtr->nth[_j][_k]!=_i)
      {
      if ((unsigned char)th_st[dtr->nth[_j][_k]]!=(unsigned char)-1) goto SKIP_EDGE_UNION;
      _j=dtr->nth[_j][_k];
           if (dtr->th[_j][0]==_l) _k=0;
      else if (dtr->th[_j][1]==_l) _k=1;
      else if (dtr->th[_j][2]==_l) _k=2;
      else                         _k=3;   
           if ( (_k!=0)&&(dtr->th[_j][0]!=dtr->th[_i][0])&&(dtr->th[_j][0]!=dtr->th[_i][1]) ) _l=dtr->th[_j][0];
      else if ( (_k!=1)&&(dtr->th[_j][1]!=dtr->th[_i][0])&&(dtr->th[_j][1]!=dtr->th[_i][1]) ) _l=dtr->th[_j][1];
      else if ( (_k!=2)&&(dtr->th[_j][2]!=dtr->th[_i][0])&&(dtr->th[_j][2]!=dtr->th[_i][1]) ) _l=dtr->th[_j][2];
      else                                                                                    _l=dtr->th[_j][3];
      _d+=vol[_j];
      }
    if (_d>min_th_vol)
      {
      //Mark ths 
      _ir=1., _j=_i, _k=3, _l=dtr->th[_i][2], th_st[_i]=0;
      while (dtr->nth[_j][_k]!=_i)
        {
        _j=dtr->nth[_j][_k];
             if (dtr->th[_j][0]==_l) _k=0;
        else if (dtr->th[_j][1]==_l) _k=1;
        else if (dtr->th[_j][2]==_l) _k=2;
        else                         _k=3;   
             if ( (_k!=0)&&(dtr->th[_j][0]!=dtr->th[_i][0])&&(dtr->th[_j][0]!=dtr->th[_i][1]) ) _l=dtr->th[_j][0];
        else if ( (_k!=1)&&(dtr->th[_j][1]!=dtr->th[_i][0])&&(dtr->th[_j][1]!=dtr->th[_i][1]) ) _l=dtr->th[_j][1];
        else if ( (_k!=2)&&(dtr->th[_j][2]!=dtr->th[_i][0])&&(dtr->th[_j][2]!=dtr->th[_i][1]) ) _l=dtr->th[_j][2];
        else                                                                                    _l=dtr->th[_j][3];
        //Update cm
        th_cm[_i].i+=th_cm[_j].i, th_cm[_i].j+=th_cm[_j].j, th_cm[_i].k+=th_cm[_j].k, _ir+=1.;
        th_st[_j]=(unsigned char)-127;
        }
      th_cm[_i].i/=_ir, th_cm[_i].j/=_ir, th_cm[_i].k/=_ir;
      }
    else SKIP_EDGE_UNION: ; //no circular th sequence available or total volume is insufficient
    }

//Stage III. Calc Es.
_i=dtr->size_th;
while (_i--)
  {
  if ( ((unsigned char)th_st[_i]!=(unsigned char)-127)&&((unsigned char)th_st[_i]!=(unsigned char)-1) )
    {
    _l=n, E[_i].i=E[_i].j=E[_i].k=0.;
    while (_l--)
      {
      _r.i=rvecs[_l].i-th_cm[_i].i, _r.j=rvecs[_l].j-th_cm[_i].j, _r.k=rvecs[_l].k-th_cm[_i].k; 
      _d=1./(calc_vec_norm(&_r)+SMALL2), _d*=sqrt(_d);
      E[_i].i+=_d*q[_l]*_r.i, E[_i].j+=_d*q[_l]*_r.j, E[_i].k+=_d*q[_l]*_r.k;
      }
    multiple_vec_scalar(&E[_i],&E[_i],COULOMB_K), ee[_i]=calc_vec_norm(&E[_i]);
    }
  else ee[_i]=-1.;
  }

//Stage IV. Fill th and keep them sorted
//Stage IV.1. Presort ths to reduce job on nonfillable ths
if (!(di_qsort(dtr->size_th,ee,(int*)th_id))) { ERROR_EXIT: free(th_id); return FALSE; }
_i=dtr->size_th; while (ee[_i-1]<0.) _i--; 
//Stage IV.2. Fill the rest
while ( (_i)&&(*nwats<max_nwats) )
  {
  //Solvating th_id[0] - the biggest E vector lenght
  if (!th_st[*th_id])
    {//Paste a single wat
    if (!(nwat=fill_th_union(&rvecs[n+*nwats*3],&th_cm[*th_id],&E[*th_id],w_db))) goto ERROR_EXIT;
    }
  else 
    {//Solvate accordingly to db
    if (!(nwat=fill_th(&rvecs[n+*nwats*3],dtr->th[*th_id],dtr->lvec,&E[*th_id],w_db))) goto ERROR_EXIT;
    }
  //Update E, ee and sort thetrahedrons
  //Calculate dipole
  p.i=p.j=p.k=0., _k=nwat;
  while (_k--) 
    {
    p.i+=w_db->qO*rvecs[n+(*nwats+_k)*3+0].i+w_db->qH*rvecs[n+(*nwats+_k)*3+1].i+w_db->qH*rvecs[n+(*nwats+_k)*3+2].i;
    p.j+=w_db->qO*rvecs[n+(*nwats+_k)*3+0].j+w_db->qH*rvecs[n+(*nwats+_k)*3+1].j+w_db->qH*rvecs[n+(*nwats+_k)*3+2].j;
    p.k+=w_db->qO*rvecs[n+(*nwats+_k)*3+0].k+w_db->qH*rvecs[n+(*nwats+_k)*3+1].k+w_db->qH*rvecs[n+(*nwats+_k)*3+2].k;
    }
  //Calculate dipole's electric field 
  _j=_i;
  while (--_j)
    {
    _r.i=th_cm[*th_id].i-th_cm[th_id[_j]].i, _r.j=th_cm[*th_id].j-th_cm[th_id[_j]].j, _r.k=th_cm[*th_id].k-th_cm[th_id[_j]].k;
    _d=1./calc_vec_norm(&_r), _ir=sqrt(_d), multiple_vec_scalar(&_r,&_r,_ir); 
    _d*=COULOMB_K*_ir, multiple_vec_scalar(&_r,&_r,3.*calc_vec_vec_scalar_product(&_r,&p));
    E[th_id[_j]].i+=_d*(_r.i-p.i), E[th_id[_j]].j+=_d*(_r.j-p.j), E[th_id[_j]].k+=_d*(_r.k-p.k);
    }
//  _j=_i;
//  while (--_j)
//    {
//    _k=nwat;
//    while (_k--)
//      {
//      _r.i=rvecs[n+(*nwats+_k)*3+0].i-th_cm[th_id[_j]].i, _r.j=rvecs[n+(*nwats+_k)*3+0].j-th_cm[th_id[_j]].j, _r.k=rvecs[n+(*nwats+_k)*3+0].k-th_cm[th_id[_j]].k;
//      _d=1./(calc_vec_norm(&_r)+SMALL2), _d*=COULOMB_K*sqrt(_d);
//      E[th_id[_j]].i+=_d*w_db->qO*_r.i, E[th_id[_j]].j+=_d*w_db->qO*_r.j, E[th_id[_j]].k+=_d*w_db->qO*_r.k;
//      _r.i=rvecs[n+(*nwats+_k)*3+1].i-th_cm[th_id[_j]].i, _r.j=rvecs[n+(*nwats+_k)*3+1].j-th_cm[th_id[_j]].j, _r.k=rvecs[n+(*nwats+_k)*3+1].k-th_cm[th_id[_j]].k;
//      _d=1./(calc_vec_norm(&_r)+SMALL2), _d*=COULOMB_K*sqrt(_d);
//      E[th_id[_j]].i+=_d*w_db->qH*_r.i, E[th_id[_j]].j+=_d*w_db->qH*_r.j, E[th_id[_j]].k+=_d*w_db->qH*_r.k;
//      _r.i=rvecs[n+(*nwats+_k)*3+2].i-th_cm[th_id[_j]].i, _r.j=rvecs[n+(*nwats+_k)*3+2].j-th_cm[th_id[_j]].j, _r.k=rvecs[n+(*nwats+_k)*3+2].k-th_cm[th_id[_j]].k;
//      _d=1./(calc_vec_norm(&_r)+SMALL2), _d*=COULOMB_K*sqrt(_d);
//      E[th_id[_j]].i+=_d*w_db->qH*_r.i, E[th_id[_j]].j+=_d*w_db->qH*_r.j, E[th_id[_j]].k+=_d*w_db->qH*_r.k;
//      }
//    ee[_j]=calc_vec_norm(&E[th_id[_j]]);
//    }
  *nwats+=nwat;

  //Bubble top th into the bottom
  for (_k=*th_id, _d=*ee, _j=1; _j<_i; _j++) { th_id[_j-1]=th_id[_j], ee[_j-1]=ee[_j]; }
  th_id[--_i]=_k, ee[_i]=_d;
//  if (!(di_qsort(_i,ee,(int*)th_id))) goto ERROR_EXIT;
  }

//Stage V. Filling done, free && exit
free(th_id);
ylib_errno=YERR_OK;
return TRUE;
}

//-------------------------------------- M O L E C U L E S     O P E R A T I O N    R O U T I N E S -------------------------------------//

