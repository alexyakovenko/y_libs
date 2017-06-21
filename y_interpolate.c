
//This module keeps various interpolation routines

#include "y_interpolate.h"

extern unsigned int ylib_errno;

//Inverse quadratic interpolation i.e. extreme x of parabola spanned over points a < b < c
inline char inverse_quadratic_interpolation(register double *x,register double xa,register double xb,register double xc,register double fa,register double fb,register double fc)
{
register double _d;
if (fabs(_d=(xb-xa)*(fb-fc)-(xb-xc)*(fb-fa))<TINY) return FALSE; //colinear case
*x=xb-((xb-xa)*(xb-xa)*(fb-fc)-(xb-xc)*(xb-xc)*(fb-fa))/(2.*_d);
return TRUE; 
}

//------------------------- M O N O T O N E   C U B I C   I N T E R P O L A T I O N --------------------------------

//The interpolation polynome derivative model on the left
inline double _calc_mc_interpolation_derivative_left(register double x[4],register double f[4])  
{
return ((2.*(x[1]-x[0])+x[3]-x[2])*(f[1]-f[0])/(x[1]-x[0])-(x[1]-x[0])*(f[3]-f[2])/(x[3]-x[2]))/(x[1]-x[0]+x[3]-x[2]);
}
//The interpolation polynome derivative model
inline double _calc_mc_interpolation_derivative(register double s0,register double s1,register double dx0,register double dx1) 
{
return (3.*s1*s0/(s1*s1+s1*s0+s0*s0))*(dx0*s1+dx1*s0)/(dx0+dx1);
}
//The interpolation polynome derivative model on the right
inline double _calc_mc_interpolation_derivative_right(register double x[4],register double f[4]) 
{
return ((2.*(x[3]-x[2])+(x[1]-x[0]))*(f[3]-f[2])/(x[3]-x[2])-(x[3]-x[2])*(f[1]-f[0])/(x[1]-x[0]))/(x[3]-x[2]+x[1]-x[0]);
}
//The interpolation polynome
inline void _calc_mc_interpolation_polynome(double c[4],register double f0,register double g0,register double g1,register double s,register double dx0)
{
c[0]=f0, c[1]=g0, c[2]=((s-g1)+2.*(s-g0))/dx0, c[3]=-((s-g1)+(s-g0))/dx0/dx0;
}
//The interpolation itself
inline double _calc_mc_interpolation(register double c[4],register double x0,register double x)
{
register double dx=(x-x0);
return ((c[3]*dx+c[2])*dx+c[1])*dx+c[0];
}			  
			  
//The calculator. It generates aligned massive of values for the common t0+n*dt checkpoints using arbitrary t,v sets of datapoints.
//To get missing values, monotonic three-cubic interpolation is used (gl is the 'left monotonic' gradient).
//NOTE. The uncovered area in w is marked as NAN;
char cubic_align_dataset_in_time(register double t0,register double dt,register unsigned int n,register double *w,register unsigned int size,register double *t,register double *v,double *g)
{
register unsigned int _i, _j;
register double dt0, dt1, s0, s1, g1, g0;
double c[4];

//Stage 1.1. Check points availability
if ( (size<4)||(t0>t[size-1])||(t[0]+(double)n*dt<t[0]) ) { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
else //Do cubic interpolation
    {
    //Stage 1. INITIALIZATION
    _i=0; while (t0+(double)_i*dt<t[0]) w[_i++]=(double)NAN;
    if ( (t0>=t[0])&&(t0<t[1]) ) 
      {
      _j=0;
      dt0=t[1]-t[0], s0=(v[1]-v[0])/dt0;
      if ( (isnan(*g))) g0=_calc_mc_interpolation_derivative_left(&t[0],&v[0]); else g0=*g; 
      }
    else 
      {
      _j=0; while ( (t0>=t[_j])&&(t0<t[_j+1])) _j++;
      dt0=t[_j+1]-t[_j], s0=(v[_j+1]-v[_j])/dt0;
      g0=_calc_mc_interpolation_derivative((v[_j]-v[_j-1])/(t[_j]-t[_j-1]),s0,(t[_j]-t[_j-1]),dt0);
      }
    //Stage 2. INTERPOLATION
    for (; _j<size-1; _j++, dt0=dt1, s0=s1, g0=g1)
      {
      dt1=t[_j+2]-t[_j+1];
      s1=(v[_j+2]-v[_j+1])/dt1;
      g1=_calc_mc_interpolation_derivative(s0,s1,dt0,dt1);
      if (t0+(double)_i*dt<t[_j+1]) 
        {
        _calc_mc_interpolation_polynome(c,v[_j],g0,g1,s0,dt0);
        do { w[_i]=_calc_mc_interpolation(c,t[_j],t0+(double)_i*dt); if (++_i==n) return TRUE; } while (t0+(double)_i*dt<t[_j+1]);
        }
      }
    //Stage 3. FINALIZATION
    *g=g1=_calc_mc_interpolation_derivative_right(&t[size-4],&v[size-4]);
    if (t0+(double)_i*dt<=t[size-1]) 
      {
      _calc_mc_interpolation_polynome(c,v[_j],g0,g1,s0,dt0);
      do { w[_i]=_calc_mc_interpolation(c,t[_j],t0+(double)_i*dt); if (++_i==n) return TRUE; } while (t0+(double)_i*dt<=t[size-1]);
      }
    while (_i!=n) w[_i++]=(double)NAN;    
    }
//Exiting
return TRUE;
}

//--------------------- I N T E R P O L A T I O N     O N     G R I D -----------------------------------------

//This function creates grid for 3D char grid
t_cgrid* alloc_cgrid(unsigned int ni,unsigned int nj,unsigned int nk)
{
unsigned int _i,_j;
t_cgrid *cgrid;
if (!(cgrid=malloc(sizeof(t_cgrid)+(sizeof(char**)+(sizeof(char*)+sizeof(char)*(size_t)nk)*(size_t)nj)*(size_t)ni))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { cgrid->c=(void*)cgrid+sizeof(t_cgrid), cgrid->len.i=ni, cgrid->len.j=nj, cgrid->len.k=nk; }
for (_i=0;_i<ni;_i++)
  {
  cgrid->c[_i]=(void*)cgrid->c+sizeof(char**)*(size_t)ni+sizeof(char*)*(size_t)nj*(size_t)_i;
  for (_j=0;_j<nj;_j++)
    cgrid->c[_i][_j]=(void*)cgrid->c+sizeof(char**)*(size_t)ni+sizeof(char*)*(size_t)ni*(size_t)nj+(sizeof(char)*(size_t)nk)*((size_t)nj*(size_t)_i+(size_t)_j);
  }
return cgrid;
}
//This function creates grid for 3D char grid
t_igrid* alloc_igrid(unsigned int ni,unsigned int nj,unsigned int nk)
{
unsigned int _i,_j;
t_igrid *igrid;
if (!(igrid=malloc(sizeof(t_igrid)+(sizeof(int**)+(sizeof(int*)+sizeof(int)*(size_t)nk)*(size_t)nj)*(size_t)ni))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { igrid->i=(void*)igrid+sizeof(t_igrid), igrid->len.i=ni, igrid->len.j=nj, igrid->len.k=nk; }
for (_i=0;_i<ni;_i++)
  {
  igrid->i[_i]=(void*)igrid->i+sizeof(int**)*(size_t)ni+sizeof(int*)*(size_t)nj*(size_t)_i;
  for (_j=0;_j<nj;_j++)
    igrid->i[_i][_j]=(void*)igrid->i+sizeof(int**)*(size_t)ni+sizeof(int*)*(size_t)ni*(size_t)nj+(sizeof(int)*(size_t)nk)*((size_t)nj*(size_t)_i+(size_t)_j);
  }
return igrid;
}
//This function creates grid for 3D double interpolation
t_dgrid* alloc_dgrid(unsigned int ni,unsigned int nj,unsigned int nk)
{
unsigned int _i,_j;
t_dgrid *dgrid;
if (!(dgrid=malloc(sizeof(t_dgrid)+(sizeof(double**)+(sizeof(double*)+sizeof(double)*(size_t)nk)*(size_t)nj)*(size_t)ni))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { dgrid->d=(void*)dgrid+sizeof(t_dgrid), dgrid->len.i=ni, dgrid->len.j=nj, dgrid->len.k=nk; }
for (_i=0;_i<ni;_i++)
  {
  dgrid->d[_i]=(void*)dgrid->d+sizeof(double**)*(size_t)ni+sizeof(double*)*(size_t)nj*(size_t)_i;
  for (_j=0;_j<nj;_j++)
    dgrid->d[_i][_j]=(void*)dgrid->d+sizeof(double**)*(size_t)ni+sizeof(double*)*(size_t)ni*(size_t)nj+(sizeof(double)*(size_t)nk)*((size_t)nj*(size_t)_i+(size_t)_j);
  }
return dgrid;
}
//This function creates grid fot tricubic interpolation
t_tcgrid* alloc_tcgrid(unsigned int ni,unsigned int nj,unsigned int nk)
{
unsigned int _i,_j;
t_tcgrid *tcgrid;
if (!(tcgrid=malloc(sizeof(t_tcgrid)+(sizeof(double**)+(sizeof(double*)+sizeof(double)*(size_t)nk*64)*(size_t)nj)*(size_t)ni))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { tcgrid->d=(void*)tcgrid+sizeof(t_tcgrid), tcgrid->len.i=ni, tcgrid->len.j=nj, tcgrid->len.k=nk; }
for (_i=0;_i<ni;_i++)
  {
  tcgrid->d[_i]=(void*)tcgrid->d+sizeof(double**)*(size_t)ni+sizeof(double*)*(size_t)nj*(size_t)_i;
  for (_j=0;_j<nj;_j++)
    tcgrid->d[_i][_j]=(void*)tcgrid->d+sizeof(double**)*(size_t)ni+sizeof(double*)*(size_t)ni*(size_t)nj+(sizeof(double)*64*(size_t)nk)*((size_t)nj*(size_t)_i+(size_t)_j);
  }
return tcgrid;
}
//This function read 3D double grid file
t_cgrid *read_cgrid(FILE *in)
{
t_len len;
t_vec ori;
t_cgrid *cgrid;
double sp;
if ( (fread(&len.i,sizeof(int),0x1,in)!=0x1)||(len.i!=Y_MAGIC)||(fread(&len,sizeof(t_len),0x1,in)!=0x1)||
     (fread(&ori,sizeof(t_vec),0x1,in)!=0x1)||(fread(&sp,sizeof(double),0x1,in)!=0x1)                   )
  { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if ( ((int)len.i<=0)||((int)len.j<=0)||((int)len.k<=0) ) { ylib_errno=YERROR_USER; return FALSE; }
if (!(cgrid=(t_cgrid*)alloc_cgrid(len.i,len.j,len.k))) return FALSE;
else { cgrid->ori.i=ori.i, cgrid->ori.j=ori.j, cgrid->ori.k=ori.k, cgrid->sp=sp; }
if (fread(cgrid->c[0][0],sizeof(char),(size_t)len.i*(size_t)len.j*(size_t)len.k,in)!=(size_t)len.i*(size_t)len.j*(size_t)len.k) { free(cgrid->c); free(cgrid); goto LABEL_IO_ERROR; }
return cgrid;
}
//This function read 3D double grid file
t_igrid *read_igrid(FILE *in)
{
t_len len;
t_vec ori;
t_igrid *igrid;
double sp;
if ( (fread(&len.i,sizeof(int),0x1,in)!=0x1)||(len.i!=Y_MAGIC)||(fread(&len,sizeof(t_len),0x1,in)!=0x1)||
     (fread(&ori,sizeof(t_vec),0x1,in)!=0x1)||(fread(&sp,sizeof(double),0x1,in)!=0x1)                   )
  { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if ( ((int)len.i<=0)||((int)len.j<=0)||((int)len.k<=0) ) { ylib_errno=YERROR_USER; return FALSE; }
if (!(igrid=(t_igrid*)alloc_igrid(len.i,len.j,len.k))) return FALSE;
else { igrid->ori.i=ori.i, igrid->ori.j=ori.j, igrid->ori.k=ori.k, igrid->sp=sp; }
if (fread(igrid->i[0][0],sizeof(int),(size_t)len.i*(size_t)len.j*(size_t)len.k,in)!=(size_t)len.i*(size_t)len.j*(size_t)len.k) { free(igrid->i); free(igrid); goto LABEL_IO_ERROR; }
return igrid;
}
//This function read 3D double grid file
t_dgrid *read_dgrid(FILE *in)
{
t_len len;
t_vec ori;
t_dgrid *dgrid;
double sp;
if ( (fread(&len.i,sizeof(int),0x1,in)!=0x1)||(len.i!=Y_MAGIC)||(fread(&len,sizeof(t_len),0x1,in)!=0x1)||
     (fread(&ori,sizeof(t_vec),0x1,in)!=0x1)||(fread(&sp,sizeof(double),0x1,in)!=0x1)                   )
  { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if ( ((int)len.i<=0)||((int)len.j<=0)||((int)len.k<=0) ) { ylib_errno=YERROR_USER; return FALSE; }
if (!(dgrid=(t_dgrid*)alloc_dgrid(len.i,len.j,len.k))) return FALSE;
else { dgrid->ori.i=ori.i, dgrid->ori.j=ori.j, dgrid->ori.k=ori.k, dgrid->sp=sp; }
if (fread(dgrid->d[0][0],sizeof(double),(size_t)len.i*(size_t)len.j*(size_t)len.k,in)!=(size_t)len.i*(size_t)len.j*(size_t)len.k) { free(dgrid->d); free(dgrid); goto LABEL_IO_ERROR; }
return dgrid;
}
//This function read grid file
t_tcgrid *read_tcgrid(FILE *in)
{
t_len len;
t_vec ori;
t_tcgrid *tcgrid;
double sp;
if ( (fread(&len.i,sizeof(int),0x1,in)!=0x1)||(len.i!=Y_MAGIC)||(fread(&len,sizeof(t_len),0x1,in)!=0x1)||
     (fread(&ori,sizeof(t_vec),0x1,in)!=0x1)||(fread(&sp,sizeof(double),0x1,in)!=0x1)                   )
  { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if ( ((int)len.i<=0)||((int)len.j<=0)||((int)len.k<=0) ) { ylib_errno=YERROR_USER; return FALSE; }
if (!(tcgrid=(t_tcgrid*)alloc_tcgrid(len.i,len.j,len.k))) return FALSE;
else { tcgrid->ori.i=ori.i, tcgrid->ori.j=ori.j, tcgrid->ori.k=ori.k, tcgrid->sp=sp; }
if (fread(tcgrid->d[0][0],sizeof(double)*64,(size_t)len.i*(size_t)len.j*(size_t)len.k,in)!=(size_t)len.i*(size_t)len.j*(size_t)len.k) { free(tcgrid->d); free(tcgrid); goto LABEL_IO_ERROR; }
return tcgrid;
}

//This function writes 3D char grid to file
char write_cgrid(FILE *out,t_cgrid *cgrid)
{
unsigned int i;
i=Y_MAGIC;
if ( (!out)||(fwrite(&i,sizeof(int),0x1,out)!=0x1)||(fwrite(&cgrid->len,sizeof(t_len),0x1,out)!=0x1)  ||
     (fwrite(&cgrid->ori,sizeof(t_vec),0x1,out)!=0x1)||(fwrite(&cgrid->sp,sizeof(double),0x1,out)!=0x1)||
     (fwrite(cgrid->c[0][0],sizeof(char),(size_t)cgrid->len.i*(size_t)cgrid->len.j*(size_t)cgrid->len.k,out)!=(size_t)cgrid->len.i*(size_t)cgrid->len.j*(size_t)cgrid->len.k) )
  { ylib_errno=YERROR_IO; return FALSE; }
return TRUE;
}
//This function writes 3D int grid to file
char write_igrid(FILE *out,t_igrid *igrid)
{
unsigned int i;
i=Y_MAGIC;
if ( (!out)||(fwrite(&i,sizeof(int),0x1,out)!=0x1)||(fwrite(&igrid->len,sizeof(t_len),0x1,out)!=0x1)  ||
     (fwrite(&igrid->ori,sizeof(t_vec),0x1,out)!=0x1)||(fwrite(&igrid->sp,sizeof(double),0x1,out)!=0x1)||
     (fwrite(igrid->i[0][0],sizeof(int),(size_t)igrid->len.i*(size_t)igrid->len.j*(size_t)igrid->len.k,out)!=(size_t)igrid->len.i*(size_t)igrid->len.j*(size_t)igrid->len.k) )
  { ylib_errno=YERROR_IO; return FALSE; }
return TRUE;
}
//This function writes 3D double grid to file
char write_dgrid(FILE *out,t_dgrid *dgrid)
{
unsigned int i;
i=Y_MAGIC;
if ( (!out)||(fwrite(&i,sizeof(int),0x1,out)!=0x1)||(fwrite(&dgrid->len,sizeof(t_len),0x1,out)!=0x1)  ||
     (fwrite(&dgrid->ori,sizeof(t_vec),0x1,out)!=0x1)||(fwrite(&dgrid->sp,sizeof(double),0x1,out)!=0x1)||
     (fwrite(dgrid->d[0][0],sizeof(double),(size_t)dgrid->len.i*(size_t)dgrid->len.j*(size_t)dgrid->len.k,out)!=(size_t)dgrid->len.i*(size_t)dgrid->len.j*(size_t)dgrid->len.k) )
  { ylib_errno=YERROR_IO; return FALSE; }
return TRUE;
}
//This function writes tricubic grid to file
char write_tcgrid(FILE *out,t_tcgrid *tcgrid)
{
unsigned int i;
i=Y_MAGIC;
if ( (!out)||(fwrite(&i,sizeof(int),0x1,out)!=0x1)||(fwrite(&tcgrid->len,sizeof(t_len),0x1,out)!=0x1)    ||
     (fwrite(&tcgrid->ori,sizeof(t_vec),0x1,out)!=0x1)||(fwrite(&tcgrid->sp,sizeof(double),0x1,out)!=0x1)||
     (fwrite(tcgrid->d[0][0],sizeof(double)*64,(size_t)tcgrid->len.i*(size_t)tcgrid->len.j*(size_t)tcgrid->len.k,out)!=(size_t)tcgrid->len.i*(size_t)tcgrid->len.j*(size_t)tcgrid->len.k) )
  { ylib_errno=YERROR_IO; return FALSE; }
return TRUE;
}


//--------------------- T R I L I N E A R     P A R T ------------------------------//

//This function calculates values and its derivaives over grid
inline double calc_trilinear_interpolation(double i, double j, double k,
                                           double f000, double f100, double f010, double f001,
                                           double f110, double f101, double f011, double f111)
{
return ((f000*(1.00-k)+f001*k)*(1.00-j)+(f010*(1.00-k)+f110*k)*j)*(1.00-i)+((f100*(1.00-k)+f101*k)*(1.00-j)+(f011*(1.00-k)+f111*k)*j)*i;
}
//This function calculates value of trilinear interpolation on grid
char calc_linear_interpolation_value_dfindif(double *f,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map)
{
t_vec _r;
t_len _n;
if ( ((_r.i=r->i-ori->i)<sp)||((_n.i=(int)(_r.i/sp))>ni-2)||((_r.j=r->j-ori->j)<sp)||((_n.j=(int)(_r.j/sp))>nj-2)||((_r.k=r->k-ori->k)<sp)||((_n.k=(int)(_r.k/sp))>nk-2) ) return FALSE;
else { _r.i=(_r.i-_n.i*sp)/sp, _r.j=(_r.j-_n.j*sp)/sp, _r.k=(_r.k-_n.k*sp)/sp; } // r is {0...1}
*f=calc_trilinear_interpolation(_r.i,_r.j,_r.k,a_map[_n.i][_n.j][_n.k],a_map[_n.i+1][_n.j][_n.k],a_map[_n.i][_n.j+1][_n.k],a_map[_n.i][_n.j][_n.k+1],
                                a_map[_n.i+1][_n.j+1][_n.k],a_map[_n.i+1][_n.j][_n.k+1],a_map[_n.i][_n.j+1][_n.k+1],a_map[_n.i+1][_n.j+1][_n.k+1]); 
return TRUE;
}
//This function calculates values and its derivaives over grid
inline double calc_dtrilinear_interpolation(double i, double j, double k,
                                            double f000, double f100, double f010, double f001,
                                            double f110, double f101, double f011, double f111,t_vec *df)
{
df->i=((f100-f000)*(1.00-k)+(f101-f001)*k)*(1.00-j)+((f110-f010)*(1.00-k)+(f111-f011)*k)*j;
df->j=((f010-f000)*(1.00-k)+(f011-f001)*k)*(1.00-i)+((f110-f100)*(1.00-k)+(f111-f101)*k)*i;
df->k=((f001-f000)*(1.00-j)+(f011-f010)*j)*(1.00-i)+((f101-f100)*(1.00-j)+(f111-f110)*j)*i;
return ((f000*(1.00-k)+f001*k)*(1.00-j)+(f010*(1.00-k)+f011*k)*j)*(1.00-i)+((f100*(1.00-k)+f101*k)*(1.00-j)+(f110*(1.00-k)+f111*k)*j)*i;
}
//This function calculates gradients of trilinear interpolation on grid
char calc_linear_interpolation_grad_dfindif(double *f,t_vec *g,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map)
{
t_vec _r;
t_len _n;
if ( ((_r.i=r->i-ori->i)<sp)||((_n.i=(int)(_r.i/sp))>ni-2)||((_r.j=r->j-ori->j)<sp)||((_n.j=(int)(_r.j/sp))>nj-2)||((_r.k=r->k-ori->k)<sp)||((_n.k=(int)(_r.k/sp))>nk-2) ) return FALSE;
else { _r.i=(_r.i-_n.i*sp)/sp, _r.j=(_r.j-_n.j*sp)/sp, _r.k=(_r.k-_n.k*sp)/sp; } // r is {0...1}
*f=calc_dtrilinear_interpolation(_r.i,_r.j,_r.k,a_map[_n.i][_n.j][_n.k],a_map[_n.i+1][_n.j][_n.k],a_map[_n.i][_n.j+1][_n.k],a_map[_n.i][_n.j][_n.k+1],
                                 a_map[_n.i+1][_n.j+1][_n.k],a_map[_n.i+1][_n.j][_n.k+1],a_map[_n.i][_n.j+1][_n.k+1],a_map[_n.i+1][_n.j+1][_n.k+1],g); 
return TRUE;
}



//----------------------   T H R E E - C U B I C     I N T E R P O L A T I O N     P A R T ------------------------------------
// Lekien F. and Marsden J. Tricubic interpolation in three dimensions Int.J.Numer.Meth.Engng 2005, 63:455-471		  


//This function calculates threecubic proto-map
char calc_threecubic_protomap(unsigned int ni,unsigned int nj,unsigned int nk,t_vec *ori,double sp,double (***p_map)[8],char (*funct)(double *f,double *dfdx,double *dfdy,double *dfdz,double *d2fdxdy,double *d2fdxdz,double *d2fdydz,double *d3fdxdydz,double rr,t_vec *r,char label,va_list stack),char label, ... )
{
va_list stack;
unsigned int _i,_j,_k;
t_vec r;
va_start(stack,label);

for (r.i=ori->i, _i=0;_i<ni;_i++, r.i+=sp)
  for (r.j=ori->j, _j=0;_j<nj;_j++, r.j+=sp)
    for (r.k=ori->k, _k=0;_k<nk;_k++, r.k+=sp)
      {
      if (!(funct(&p_map[_i][_j][_k][0],&p_map[_i][_j][_k][1],&p_map[_i][_j][_k][2],&p_map[_i][_j][_k][3],&p_map[_i][_j][_k][4],&p_map[_i][_j][_k][5],&p_map[_i][_j][_k][6],&p_map[_i][_j][_k][7],2.*sp,&r,label,stack))) return FALSE;
      p_map[_i][_j][_k][1]*=sp;
      p_map[_i][_j][_k][2]*=sp;
      p_map[_i][_j][_k][3]*=sp;
      p_map[_i][_j][_k][4]*=sp*sp;
      p_map[_i][_j][_k][5]*=sp*sp;
      p_map[_i][_j][_k][6]*=sp*sp;
      p_map[_i][_j][_k][7]*=sp*sp*sp;
      }
va_end(stack);
return TRUE;
}


//
//    z   p011-----p111        z   p7-------p8 
//    ^  /|       /|           ^  /|       /|
//    | / |   *y / |           | / |   *y / |
//    |/  |  /  /  |           |/  |  /  /  |
//    p001+----p101|           p5--+----p6  |
//    |   |/   |   |           |   |/   |   |
//    |   p010-+---p110        |   p3---+---p4
//    |  /     |  /            |  /     |  /
//    | /      | /             | /      | /
//    |/       |/              |/       |/
//    p000-----p100--> x       p1-------p2-----> x
//
//This function transforms interpolation proto-map [nx][ny][nz][8] of {fx,dx,dy,dz,dxdy,dydz,dzdx,dxdydz} into threecubic cooficients [nx-1][ny-1][nz-1][64]    
void transform_threecubic_protomap(unsigned int ni,unsigned int nj,unsigned int nk,double (***a_map)[64],double (***p_map)[8])
{
char IA[64][64] = {
{  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{ -3,  3,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  2, -2,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{ -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  9, -9, -9,  9,  0,  0,  0,  0,  6,  3, -6, -3,  0,  0,  0,  0,  6, -6,  3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  2,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{ -6,  6,  6, -6,  0,  0,  0,  0, -3, -3,  3,  3,  0,  0,  0,  0, -4,  4, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2, -2, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{ -6,  6,  6, -6,  0,  0,  0,  0, -4, -2,  4,  2,  0,  0,  0,  0, -3,  3, -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2, -1, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  4, -4, -4,  4,  0,  0,  0,  0,  2,  2, -2, -2,  0,  0,  0,  0,  2, -2,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  9, -9, -9,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  3, -6, -3,  0,  0,  0,  0,  6, -6,  3, -3,  0,  0,  0,  0,  4,  2,  2,  1,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3, -3,  3,  3,  0,  0,  0,  0, -4,  4, -2,  2,  0,  0,  0,  0, -2, -2, -1, -1,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -4, -2,  4,  2,  0,  0,  0,  0, -3,  3, -3,  3,  0,  0,  0,  0, -2, -1, -2, -1,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2, -2, -2,  0,  0,  0,  0,  2, -2,  2, -2,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0}, 
{ -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  9, -9,  0,  0, -9,  9,  0,  0,  6,  3,  0,  0, -6, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  2,  0,  0,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{ -6,  6,  0,  0,  6, -6,  0,  0, -3, -3,  0,  0,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -4,  4,  0,  0, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2, -2,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  9, -9,  0,  0, -9,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  3,  0,  0, -6, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  3, -3,  0,  0,  4,  2,  0,  0,  2,  1,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3, -3,  0,  0,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -4,  4,  0,  0, -2,  2,  0,  0, -2, -2,  0,  0, -1, -1,  0,  0}, 
{  9,  0, -9,  0, -9,  0,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  0,  3,  0, -6,  0, -3,  0,  6,  0, -6,  0,  3,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  2,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  9,  0, -9,  0, -9,  0,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  0,  3,  0, -6,  0, -3,  0,  6,  0, -6,  0,  3,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  2,  0,  2,  0,  1,  0}, 
{-27, 27, 27,-27, 27,-27,-27, 27,-18, -9, 18,  9, 18,  9,-18, -9,-18, 18, -9,  9, 18,-18,  9, -9,-18, 18, 18,-18, -9,  9,  9, -9,-12, -6, -6, -3, 12,  6,  6,  3,-12, -6, 12,  6, -6, -3,  6,  3,-12, 12, -6,  6, -6,  6, -3,  3, -8, -4, -4, -2, -4, -2, -2, -1}, 
{ 18,-18,-18, 18,-18, 18, 18,-18,  9,  9, -9, -9, -9, -9,  9,  9, 12,-12,  6, -6,-12, 12, -6,  6, 12,-12,-12, 12,  6, -6, -6,  6,  6,  6,  3,  3, -6, -6, -3, -3,  6,  6, -6, -6,  3,  3, -3, -3,  8, -8,  4, -4,  4, -4,  2, -2,  4,  4,  2,  2,  2,  2,  1,  1}, 
{ -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0, -3,  0,  3,  0,  3,  0, -4,  0,  4,  0, -2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -2,  0, -1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0, -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0, -3,  0,  3,  0,  3,  0, -4,  0,  4,  0, -2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -2,  0, -1,  0, -1,  0}, 
{ 18,-18,-18, 18,-18, 18, 18,-18, 12,  6,-12, -6,-12, -6, 12,  6,  9, -9,  9, -9, -9,  9, -9,  9, 12,-12,-12, 12,  6, -6, -6,  6,  6,  3,  6,  3, -6, -3, -6, -3,  8,  4, -8, -4,  4,  2, -4, -2,  6, -6,  6, -6,  3, -3,  3, -3,  4,  2,  4,  2,  2,  1,  2,  1}, 
{-12, 12, 12,-12, 12,-12,-12, 12, -6, -6,  6,  6,  6,  6, -6, -6, -6,  6, -6,  6,  6, -6,  6, -6, -8,  8,  8, -8, -4,  4,  4, -4, -3, -3, -3, -3,  3,  3,  3,  3, -4, -4,  4,  4, -2, -2,  2,  2, -4,  4, -4,  4, -2,  2, -2,  2, -2, -2, -2, -2, -1, -1, -1, -1}, 
{  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{ -6,  6,  0,  0,  6, -6,  0,  0, -4, -2,  0,  0,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  4, -4,  0,  0, -4,  4,  0,  0,  2,  2,  0,  0, -2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -4, -2,  0,  0,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0, -2, -1,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0, -2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0}, 
{ -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0, -4,  0, -2,  0,  4,  0,  2,  0, -3,  0,  3,  0, -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0, -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -4,  0, -2,  0,  4,  0,  2,  0, -3,  0,  3,  0, -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0, -2,  0, -1,  0}, 
{ 18,-18,-18, 18,-18, 18, 18,-18, 12,  6,-12, -6,-12, -6, 12,  6, 12,-12,  6, -6,-12, 12, -6,  6,  9, -9, -9,  9,  9, -9, -9,  9,  8,  4,  4,  2, -8, -4, -4, -2,  6,  3, -6, -3,  6,  3, -6, -3,  6, -6,  3, -3,  6, -6,  3, -3,  4,  2,  2,  1,  4,  2,  2,  1}, 
{-12, 12, 12,-12, 12,-12,-12, 12, -6, -6,  6,  6,  6,  6, -6, -6, -8,  8, -4,  4,  8, -8,  4, -4, -6,  6,  6, -6, -6,  6,  6, -6, -4, -4, -2, -2,  4,  4,  2,  2, -3, -3,  3,  3, -3, -3,  3,  3, -4,  4, -2,  2, -4,  4, -2,  2, -2, -2, -1, -1, -2, -2, -1, -1}, 
{  4,  0, -4,  0, -4,  0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0, -2,  0, -2,  0,  2,  0, -2,  0,  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
{  0,  0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0, -4,  0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0, -2,  0, -2,  0,  2,  0, -2,  0,  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0}, 
{-12, 12, 12,-12, 12,-12,-12, 12, -8, -4,  8,  4,  8,  4, -8, -4, -6,  6, -6,  6,  6, -6,  6, -6, -6,  6,  6, -6, -6,  6,  6, -6, -4, -2, -4, -2,  4,  2,  4,  2, -4, -2,  4,  2, -4, -2,  4,  2, -3,  3, -3,  3, -3,  3, -3,  3, -2, -1, -2, -1, -2, -1, -2, -1}, 
{  8, -8, -8,  8, -8,  8,  8, -8,  4,  4, -4, -4, -4, -4,  4,  4,  4, -4,  4, -4, -4,  4, -4,  4,  4, -4, -4,  4,  4, -4, -4,  4,  2,  2,  2,  2, -2, -2, -2, -2,  2,  2, -2, -2,  2,  2, -2, -2,  2, -2,  2, -2,  2, -2,  2, -2,  1,  1,  1,  1,  1,  1,  1,  1}};
double _t[64];
unsigned int _i,_j,_k,_I,_J,_K,_p,_q;

//Calculate three-cubic map
for (_i=0;_i<ni-1;_i++)
  for (_I=_i+1, _j=0;_j<nj-1;_j++)
    for (_J=_j+1, _k=0;_k<nk-1;_k++)
      {
      _K=_k+1;
      _t[ 0]=p_map[_i][_j][_k][0],  //f
      _t[ 1]=p_map[_I][_j][_k][0],  //f
      _t[ 2]=p_map[_i][_J][_k][0],  //f
      _t[ 3]=p_map[_I][_J][_k][0],  //f
      _t[ 4]=p_map[_i][_j][_K][0],  //f
      _t[ 5]=p_map[_I][_j][_K][0],  //f
      _t[ 6]=p_map[_i][_J][_K][0],  //f
      _t[ 7]=p_map[_I][_J][_K][0],  //f
      _t[ 8]=p_map[_i][_j][_k][1],  //df/dx
      _t[ 9]=p_map[_I][_j][_k][1],  //df/dx
      _t[10]=p_map[_i][_J][_k][1],  //df/dx
      _t[11]=p_map[_I][_J][_k][1],  //df/dx
      _t[12]=p_map[_i][_j][_K][1],  //df/dx
      _t[13]=p_map[_I][_j][_K][1],  //df/dx
      _t[14]=p_map[_i][_J][_K][1],  //df/dx
      _t[15]=p_map[_I][_J][_K][1],  //df/dx
      _t[16]=p_map[_i][_j][_k][2],  //df/dy
      _t[17]=p_map[_I][_j][_k][2],  //df/dy
      _t[18]=p_map[_i][_J][_k][2],  //df/dy
      _t[19]=p_map[_I][_J][_k][2],  //df/dy
      _t[20]=p_map[_i][_j][_K][2],  //df/dy
      _t[21]=p_map[_I][_j][_K][2],  //df/dy
      _t[22]=p_map[_i][_J][_K][2],  //df/dy
      _t[23]=p_map[_I][_J][_K][2],  //df/dy
      _t[24]=p_map[_i][_j][_k][3],  //df/dz
      _t[25]=p_map[_I][_j][_k][3],  //df/dz
      _t[26]=p_map[_i][_J][_k][3],  //df/dz
      _t[27]=p_map[_I][_J][_k][3],  //df/dz
      _t[28]=p_map[_i][_j][_K][3],  //df/dz
      _t[29]=p_map[_I][_j][_K][3],  //df/dz
      _t[30]=p_map[_i][_J][_K][3],  //df/dz
      _t[31]=p_map[_I][_J][_K][3],  //df/dz
      _t[32]=p_map[_i][_j][_k][4],  //d2f/dx/dy
      _t[33]=p_map[_I][_j][_k][4],  //d2f/dx/dy
      _t[34]=p_map[_i][_J][_k][4],  //d2f/dx/dy
      _t[35]=p_map[_I][_J][_k][4],  //d2f/dx/dy
      _t[36]=p_map[_i][_j][_K][4],  //d2f/dx/dy
      _t[37]=p_map[_I][_j][_K][4],  //d2f/dx/dy
      _t[38]=p_map[_i][_J][_K][4],  //d2f/dx/dy
      _t[39]=p_map[_I][_J][_K][4],  //d2f/dx/dy
      _t[40]=p_map[_i][_j][_k][5],  //d2f/dx/dz
      _t[41]=p_map[_I][_j][_k][5],  //d2f/dx/dz
      _t[42]=p_map[_i][_J][_k][5],  //d2f/dx/dz
      _t[43]=p_map[_I][_J][_k][5],  //d2f/dx/dz
      _t[44]=p_map[_i][_j][_K][5],  //d2f/dx/dz
      _t[45]=p_map[_I][_j][_K][5],  //d2f/dx/dz
      _t[46]=p_map[_i][_J][_K][5],  //d2f/dx/dz
      _t[47]=p_map[_I][_J][_K][5],  //d2f/dx/dz
      _t[48]=p_map[_i][_j][_k][6],  //d2f/dy/dz
      _t[49]=p_map[_I][_j][_k][6],  //d2f/dy/dz
      _t[50]=p_map[_i][_J][_k][6],  //d2f/dy/dz
      _t[51]=p_map[_I][_J][_k][6],  //d2f/dy/dz
      _t[52]=p_map[_i][_j][_K][6],  //d2f/dy/dz
      _t[53]=p_map[_I][_j][_K][6],  //d2f/dy/dz
      _t[54]=p_map[_i][_J][_K][6],  //d2f/dy/dz
      _t[55]=p_map[_I][_J][_K][6],  //d2f/dy/dz
      _t[56]=p_map[_i][_j][_k][7],  //d3f/dx/dy/dz
      _t[57]=p_map[_I][_j][_k][7],  //d3f/dx/dy/dz
      _t[58]=p_map[_i][_J][_k][7],  //d3f/dx/dy/dz
      _t[59]=p_map[_I][_J][_k][7],  //d3f/dx/dy/dz
      _t[60]=p_map[_i][_j][_K][7],  //d3f/dx/dy/dz
      _t[61]=p_map[_I][_j][_K][7],  //d3f/dx/dy/dz
      _t[62]=p_map[_i][_J][_K][7],  //d3f/dx/dy/dz
      _t[63]=p_map[_I][_J][_K][7];  //d3f/dx/dy/dz
      //Do kernel transformation
      for (_p=0;_p<64;_p++)
        for (a_map[_i][_j][_k][_p]=0., _q=0;_q<64;_q++)
          a_map[_i][_j][_k][_p]+=(double)(IA[_p][_q])*_t[_q]; 
      }
}

//This function calculates threecubic iterpolated value of function
char calc_threecubic_interpolation_value(double *f,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double (***a_map)[64])
{
double i, j, k, ii, jj, kk, ij, ik, jk, iij, ijj, jjk, jkk, iii, jjj, kkk, *_t;
unsigned int _i,_j,_k;

if ( ((i=r->i-ori->i)<0.)||((j=r->j-ori->j)<0.)||((k=r->k-ori->k)<0.) ) { DATA_CONSISTMENT_ERROR: ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
if ( ((_i=(int)(i/sp))>ni)||((_j=(int)(j/sp))>nj)||((_k=(int)(k/sp))>nk) ) goto DATA_CONSISTMENT_ERROR;
else 
  {
  i-=(double)(_i*sp), j-=(double)(_j*sp), k-=(double)(_k*sp);
  i/=sp, j/=sp, k/=sp;
  }
_t=&a_map[_i][_j][_k][0];
ii=i*i, jj=j*j, kk=k*k, ij=i*j, ik=i*k, jk=j*k, iij=ii*j, ijj=i*jj, jjk=jj*k, jkk=j*kk, iii=ii*i, jjj=jj*j, kkk=kk*k; // 15 flops
*f =_t[ 0];
*f+=_t[ 1]*i;
*f+=_t[ 2]*ii;
*f+=_t[ 3]*iii;
*f+=_t[ 4]*j;
*f+=_t[ 5]*ij;
*f+=_t[ 6]*iij;
*f+=_t[ 7]*iii*j;
*f+=_t[ 8]*jj;
*f+=_t[ 9]*ijj;
*f+=_t[10]*ii*jj;
*f+=_t[11]*iii*jj;
*f+=_t[12]*jjj;
*f+=_t[13]*i*jjj;
*f+=_t[14]*ii*jjj;
*f+=_t[15]*iii*jjj;
*f+=_t[16]*k;
*f+=_t[17]*ik;
*f+=_t[18]*ii*k;
*f+=_t[19]*iii*k;
*f+=_t[20]*jk;
*f+=_t[21]*ij*k;
*f+=_t[22]*ii*jk;
*f+=_t[23]*iii*jk;
*f+=_t[24]*jjk;
*f+=_t[25]*i*jjk;
*f+=_t[26]*ii*jjk;
*f+=_t[27]*iii*jjk;
*f+=_t[28]*jjj*k;
*f+=_t[29]*ijj*jk;
*f+=_t[30]*iij*jjk;
*f+=_t[31]*iii*jjj*k;
*f+=_t[32]*kk;
*f+=_t[33]*i*kk;
*f+=_t[34]*ii*kk;
*f+=_t[35]*iii*kk;
*f+=_t[36]*jkk;
*f+=_t[37]*i*jkk;
*f+=_t[38]*ii*jkk;
*f+=_t[39]*iii*jkk;
*f+=_t[40]*jj*kk;
*f+=_t[41]*ijj*kk;
*f+=_t[42]*iij*jkk;
*f+=_t[43]*iii*jj*kk;
*f+=_t[44]*jjj*kk;
*f+=_t[45]*ijj*jkk;
*f+=_t[46]*ii*jjj*kk;
*f+=_t[47]*iii*jjj*kk;
*f+=_t[48]*kkk;
*f+=_t[49]*i*kkk;
*f+=_t[50]*ii*kkk;
*f+=_t[51]*iii*kkk;
*f+=_t[52]*j*kkk;
*f+=_t[53]*ij*kkk;
*f+=_t[54]*iij*kkk;
*f+=_t[55]*iii*j*kkk;
*f+=_t[56]*jj*kkk;
*f+=_t[57]*ijj*kkk;
*f+=_t[58]*ii*jj*kkk;
*f+=_t[59]*iii*jj*kkk;
*f+=_t[60]*jjj*kkk;
*f+=_t[61]*i*jjj*kkk;
*f+=_t[62]*ii*jjj*kkk;
*f+=_t[63]*iii*jjj*kkk;
return TRUE; //total 195 flops (63*2+54+15)
}
//This function calculates threecubic iterpolated value of function and its first derivative
char calc_threecubic_interpolation_derivative(double *f,t_vec *g,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double (***a_map)[64])
{
double i, j, k, ii, jj, kk, ij, ik, jk, ijk, *_t;
unsigned int _i,_j,_k;

if ( ((i=r->i-ori->i)<0.)||((j=r->j-ori->j)<0.)||((k=r->k-ori->k)<0.) ) { DATA_CONSISTMENT_ERROR: ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
if ( ((_i=(int)(i/sp))>ni)||((_j=(int)(j/sp))>nj)||((_k=(int)(k/sp))>nk) ) goto DATA_CONSISTMENT_ERROR;
else 
  {
  i-=(double)(_i*sp), j-=(double)(_j*sp), k-=(double)(_k*sp);
  i/=sp, j/=sp, k/=sp;
  }
_t=&a_map[_i][_j][_k][0];
ii=i*i, jj=j*j, kk=k*k, ij=i*j, ik=i*k, jk=j*k, ijk=i*jk;
*f =_t[ 0];
*f+=_t[ 1]*i,            g->i =_t[ 1];
*f+=_t[ 2]*ii,           g->i+=_t[ 2]*2.*i;
*f+=_t[ 3]*ii*i,         g->i+=_t[ 3]*3.*ii;    
*f+=_t[ 4]*j,                                        g->j =_t[ 4];
*f+=_t[ 5]*ij,           g->i+=_t[ 5]*j,             g->j+=_t[ 5]*i;    
*f+=_t[ 6]*ij*i,         g->i+=_t[ 6]*2.*ij,         g->j+=_t[ 6]*ii;
*f+=_t[ 7]*ij*ii,        g->i+=_t[ 7]*3.*ii*j,       g->j+=_t[ 7]*ii*i;
*f+=_t[ 8]*jj,                                       g->j+=_t[ 8]*2.*j;
*f+=_t[ 9]*jj*i,         g->i+=_t[ 9]*jj,            g->j+=_t[ 9]*2.*ij;
*f+=_t[10]*ij*ij,        g->i+=_t[10]*2.*jj*i,       g->j+=_t[10]*2.*ii*j;
*f+=_t[11]*ij*ij*i,      g->i+=_t[11]*3.*ij*ij,      g->j+=_t[11]*2.*ij*ii;
*f+=_t[12]*jj*j,                                     g->j+=_t[12]*3.*jj;
*f+=_t[13]*jj*ij,        g->i+=_t[13]*jj*j,          g->j+=_t[13]*3.*jj*i;
*f+=_t[14]*ij*ij*j,      g->i+=_t[14]*2.*ij*jj,      g->j+=_t[14]*3.*ij*ij;
*f+=_t[15]*ij*ij*ij,     g->i+=_t[15]*3.*ij*ij*j,    g->j+=_t[15]*3.*ij*ij*i;
*f+=_t[16]*k,                                                                      g->k =_t[16];
*f+=_t[17]*ik,           g->i+=_t[17]*k,                                           g->k+=_t[17]*i;
*f+=_t[18]*ii*k,         g->i+=_t[18]*2.*ik,                                       g->k+=_t[18]*ii;
*f+=_t[19]*ii*ik,        g->i+=_t[19]*3.*ii*k,                                     g->k+=_t[19]*ii*i;
*f+=_t[20]*jk,                                       g->j+=_t[20]*k,               g->k+=_t[20]*j;
*f+=_t[21]*ijk,          g->i+=_t[21]*jk,            g->j+=_t[21]*ik,              g->k+=_t[21]*ij;
*f+=_t[22]*ij*ik,        g->i+=_t[22]*2.*ijk,        g->j+=_t[22]*ii*k,            g->k+=_t[22]*ij*i;
*f+=_t[23]*ijk*ii,       g->i+=_t[23]*3.*ijk*i,      g->j+=_t[23]*ik*ii,           g->k+=_t[23]*ij*ii;
*f+=_t[24]*jj*k,                                     g->j+=_t[24]*2.*jk,           g->k+=_t[24]*jj;
*f+=_t[25]*jj*ik,        g->i+=_t[25]*jj*k,          g->j+=_t[25]*2.*ijk,          g->k+=_t[25]*jj*i;
*f+=_t[26]*ij*ijk,       g->i+=_t[26]*2.*ijk*j,      g->j+=_t[26]*2.*ijk*i,        g->k+=_t[26]*ij*ij;
*f+=_t[27]*jj*ii*ik,     g->i+=_t[27]*3.*ij*ijk,     g->j+=_t[27]*2.*ii*ijk,       g->k+=_t[27]*jj*ii*i;
*f+=_t[28]*jj*jk,                                    g->j+=_t[28]*3.*jj*k,         g->k+=_t[28]*jj*j;
*f+=_t[29]*jj*ijk,       g->i+=_t[29]*jj*jk,         g->j+=_t[29]*3.*ijk*j,        g->k+=_t[29]*jj*ij;
*f+=_t[30]*ij*ij*jk,     g->i+=_t[30]*2.*jj*ijk,     g->j+=_t[30]*3.*ij*ijk,       g->k+=_t[30]*ij*ij*j;
*f+=_t[31]*ij*ij*ijk,    g->i+=_t[31]*3.*ij*ij*ik,   g->j+=_t[31]*3.*ij*ii*jk,     g->k+=_t[31]*ij*ij*ij;
*f+=_t[32]*kk,                                                                     g->k+=_t[32]*2.*k;
*f+=_t[33]*kk*i,         g->i+=_t[33]*kk,                                          g->k+=_t[33]*2.*ik;
*f+=_t[34]*ii*kk,        g->i+=_t[34]*2.*kk*i,                                     g->k+=_t[34]*2.*ik;
*f+=_t[35]*ik*ik*i,      g->i+=_t[35]*3.*ik*ik,                                    g->k+=_t[35]*2.*ik*ii;
*f+=_t[36]*jk*k,                                     g->j+=_t[36]*kk,              g->k+=_t[36]*2.*jk;
*f+=_t[37]*ij*kk,        g->i+=_t[37]*kk*j,          g->j+=_t[37]*kk*i,            g->k+=_t[37]*2.*ijk;
*f+=_t[38]*ik*ijk,       g->i+=_t[38]*2.*ijk*k,      g->j+=_t[38]*ik*ik,           g->k+=_t[38]*2.*ijk*i;
*f+=_t[39]*ijk*ii*k,     g->i+=_t[39]*3.*ijk*ik,     g->j+=_t[39]*ii*kk*i,         g->k+=_t[39]*2.*ijk*ii;
*f+=_t[40]*jj*kk,                                    g->j+=_t[40]*2.*kk*j,         g->k+=_t[40]*2.*jj*k;
*f+=_t[41]*jk*ijk,       g->i+=_t[41]*jk*jk,         g->j+=_t[41]*2.*ijk*k,        g->k+=_t[41]*2.*ijk*j;
*f+=_t[42]*ijk*ijk,      g->i+=_t[42]*2.*ijk*jk,     g->j+=_t[42]*2.*ijk*ik,       g->k+=_t[42]*2.*ijk*ij;
*f+=_t[43]*ijk*ijk*i,    g->i+=_t[43]*3.*ijk*ijk,    g->j+=_t[43]*2.*ijk*ik*i,     g->k+=_t[44]*2.*ijk*ij*i;
*f+=_t[44]*jk*jk*j,                                  g->j+=_t[44]*3.*jk*jk,        g->k+=_t[44]*2.*jk*jj;
*f+=_t[45]*jk*jk*ij,     g->i+=_t[45]*jk*jk*j,       g->j+=_t[45]*3.*jk*ijk,       g->k+=_t[45]*2.*jj*ijk;
*f+=_t[46]*jj*ijk*ik,    g->i+=_t[46]*2.*ij*jk*jk,   g->j+=_t[46]*3.*ijk*ijk,      g->k+=_t[46]*2.*ijk*jj*i;
*f+=_t[47]*ijk*ijk*ij,   g->i+=_t[47]*3.*ijk*ijk*j,  g->j+=_t[47]*3.*ijk*ijk*i,    g->k+=_t[47]*2.*ijk*ij*ij;
*f+=_t[48]*kk*k,                                                                   g->k+=_t[48]*3.*kk;
*f+=_t[49]*kk*ik,        g->i+=_t[49]*kk*k,                                        g->k+=_t[49]*3.*kk*i; 
*f+=_t[50]*kk*ik*i,      g->i+=_t[50]*2.*kk*ik,                                    g->k+=_t[50]*3.*kk*ii;
*f+=_t[51]*kk*ik*ii,     g->i+=_t[51]*3.*kk*ik*i,                                  g->k+=_t[51]*3.*ik*ik*i;
*f+=_t[52]*kk*jk,                                     g->j+=_t[52]*kk*k,           g->k+=_t[52]*3.*kk*j;   
*f+=_t[53]*kk*ijk,       g->i+=_t[53]*kk*jk,          g->j+=_t[53]*kk*ik,          g->k+=_t[53]*3.*kk*ij;
*f+=_t[54]*kk*ik*ij,     g->i+=_t[54]*2.*kk*ijk,      g->j+=_t[54]*kk*ii*k,        g->k+=_t[54]*3.*ik*ijk;
*f+=_t[55]*kk*ii*ijk,    g->i+=_t[55]*3.*kk*ii*jk,    g->j+=_t[55]*kk*ii*ik,       g->k+=_t[55]*3.*kk*ii*ij;    
*f+=_t[56]*kk*jj*k,                                   g->j+=_t[56]*2.*kk*jk,       g->k+=_t[56]*3.*kk*jj;
*f+=_t[57]*kk*jj*ik,     g->i+=_t[57]*kk*jk*j,        g->j+=_t[57]*2.*kk*ijk,      g->k+=_t[57]*3.*ijk*jk;
*f+=_t[58]*kk*ij*ijk,    g->i+=_t[58]*2.*kk*ijk*j,    g->j+=_t[58]*2.*ik*ik*jk,    g->k+=_t[58]*3.*jk*ik*ij;
*f+=_t[59]*ik*ijk*ijk,   g->i+=_t[59]*3.*ijk*ijk*k,   g->j+=_t[59]*2.*ik*ik*ijk,   g->k+=_t[59]*3.*ijk*ik*ij;
*f+=_t[60]*kk*jj*jk,                                  g->j+=_t[60]*3.*kk*jj*k,     g->k+=_t[60]*3.*kk*jj*j;
*f+=_t[61]*kk*jj*ijk,    g->i+=_t[61]*kk*jj*jk,       g->j+=_t[61]*3.*kk*j*ijk,    g->k+=_t[61]*3.*k*jj*ijk;
*f+=_t[62]*ijk*ijk*jk,   g->i+=_t[62]*2.*ijk*jk*jk,   g->j+=_t[62]*3.*ijk*ijk*k,   g->k+=_t[62]*3.*ijk*ijk*j;
*f+=_t[63]*ijk*ijk*ijk,  g->i+=_t[63]*3.*ijk*ijk*jk,  g->j+=_t[63]*3.*ijk*ijk*ik,  g->k+=_t[63]*3.*ijk*ijk*ij;
//Correct derivatives
g->i/=sp, g->j/=sp, g->k/=sp;
return TRUE;  // about 512 flops (4*4*64)
}
//This function calculates threecubic iterpolated value of function, its first and second derivatives
char calc_threecubic_interpolation_derivatives(double *f,t_vec *g,t_tensor *G,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double (***a_map)[64])
{
ylib_errno=YERROR_NIMPLEMENTED;
return FALSE;
}


//---------------------------------- H E R M E T I A N       S P L I N E ---------------------------------------------/

//This function calculates b-spline constants using finite differnces approach
inline void calc_cubic_spline_dfindif(register t_quaternion *dq,register double x)
{
dq->i=(4.-3.*x)*x-1, dq->j=x*(9.*x-10.), dq->k=(8.-9.*x)*x+1., dq->w=x*(3.*x-2.); // 14 flops
}
//This function calculates b-spline derivative constants using finite differnces approach
inline void calc_cubic_spline_findif(register t_quaternion *q,register double x)
{
q->i=x*((2.-x)*x-1.), q->j=x*x*(3.*x-5.)+1., q->k=x*((4.-3.*x)*x+1.), q->w=x*x*(x-1.); // 21 flops
}
//This function calculates tricubic interpolation value using finite differences instead of derivatives
char calc_tricubic_interpolation_value_dfindif(double *f,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map)
{
t_quaternion q;
t_vec _r;
t_len _n;
double b_map[4][4];

if ( ((_r.i=r->i-ori->i)<sp)||((_n.i=(int)(_r.i/sp))>ni-2)||((_r.j=r->j-ori->j)<sp)||((_n.j=(int)(_r.j/sp))>nj-2)||((_r.k=r->k-ori->k)<sp)||((_n.k=(int)(_r.k/sp))>nk-2) ) return FALSE;
else { _r.i=(_r.i-_n.i*sp)/sp, _r.j=(_r.j-_n.j*sp)/sp, _r.k=(_r.k-_n.k*sp)/sp; } // r is {0...1}
calc_cubic_spline_findif(&q,_r.k);
b_map[0][0]=q.i*a_map[_n.i-1][_n.j-1][_n.k-1]+q.j*a_map[_n.i-1][_n.j-1][_n.k-0]+q.k*a_map[_n.i-1][_n.j-1][_n.k+1]+q.w*a_map[_n.i-1][_n.j-1][_n.k+2];
b_map[0][1]=q.i*a_map[_n.i-1][_n.j-0][_n.k-1]+q.j*a_map[_n.i-1][_n.j-0][_n.k-0]+q.k*a_map[_n.i-1][_n.j-0][_n.k+1]+q.w*a_map[_n.i-1][_n.j-0][_n.k+2];
b_map[0][2]=q.i*a_map[_n.i-1][_n.j+1][_n.k-1]+q.j*a_map[_n.i-1][_n.j+1][_n.k-0]+q.k*a_map[_n.i-1][_n.j+1][_n.k+1]+q.w*a_map[_n.i-1][_n.j+1][_n.k+2];
b_map[0][3]=q.i*a_map[_n.i-1][_n.j+2][_n.k-1]+q.j*a_map[_n.i-1][_n.j+2][_n.k-0]+q.k*a_map[_n.i-1][_n.j+2][_n.k+1]+q.w*a_map[_n.i-1][_n.j+2][_n.k+2];
b_map[1][0]=q.i*a_map[_n.i-0][_n.j-1][_n.k-1]+q.j*a_map[_n.i-0][_n.j-1][_n.k-0]+q.k*a_map[_n.i-0][_n.j-1][_n.k+1]+q.w*a_map[_n.i-0][_n.j-1][_n.k+2];
b_map[1][1]=q.i*a_map[_n.i-0][_n.j-0][_n.k-1]+q.j*a_map[_n.i-0][_n.j-0][_n.k-0]+q.k*a_map[_n.i-0][_n.j-0][_n.k+1]+q.w*a_map[_n.i-0][_n.j-0][_n.k+2];
b_map[1][2]=q.i*a_map[_n.i-0][_n.j+1][_n.k-1]+q.j*a_map[_n.i-0][_n.j+1][_n.k-0]+q.k*a_map[_n.i-0][_n.j+1][_n.k+1]+q.w*a_map[_n.i-0][_n.j+1][_n.k+2];
b_map[1][3]=q.i*a_map[_n.i-0][_n.j+2][_n.k-1]+q.j*a_map[_n.i-0][_n.j+2][_n.k-0]+q.k*a_map[_n.i-0][_n.j+2][_n.k+1]+q.w*a_map[_n.i-0][_n.j+2][_n.k+2];
b_map[2][0]=q.i*a_map[_n.i+1][_n.j-1][_n.k-1]+q.j*a_map[_n.i+1][_n.j-1][_n.k-0]+q.k*a_map[_n.i+1][_n.j-1][_n.k+1]+q.w*a_map[_n.i+1][_n.j-1][_n.k+2];
b_map[2][1]=q.i*a_map[_n.i+1][_n.j-0][_n.k-1]+q.j*a_map[_n.i+1][_n.j-0][_n.k-0]+q.k*a_map[_n.i+1][_n.j-0][_n.k+1]+q.w*a_map[_n.i+1][_n.j-0][_n.k+2];
b_map[2][2]=q.i*a_map[_n.i+1][_n.j+1][_n.k-1]+q.j*a_map[_n.i+1][_n.j+1][_n.k-0]+q.k*a_map[_n.i+1][_n.j+1][_n.k+1]+q.w*a_map[_n.i+1][_n.j+1][_n.k+2];
b_map[2][3]=q.i*a_map[_n.i+1][_n.j+2][_n.k-1]+q.j*a_map[_n.i+1][_n.j+2][_n.k-0]+q.k*a_map[_n.i+1][_n.j+2][_n.k+1]+q.w*a_map[_n.i+1][_n.j+2][_n.k+2];
b_map[3][0]=q.i*a_map[_n.i+2][_n.j-1][_n.k-1]+q.j*a_map[_n.i+2][_n.j-1][_n.k-0]+q.k*a_map[_n.i+2][_n.j-1][_n.k+1]+q.w*a_map[_n.i+2][_n.j-1][_n.k+2];
b_map[3][1]=q.i*a_map[_n.i+2][_n.j-0][_n.k-1]+q.j*a_map[_n.i+2][_n.j-0][_n.k-0]+q.k*a_map[_n.i+2][_n.j-0][_n.k+1]+q.w*a_map[_n.i+2][_n.j-0][_n.k+2];
b_map[3][2]=q.i*a_map[_n.i+2][_n.j+1][_n.k-1]+q.j*a_map[_n.i+2][_n.j+1][_n.k-0]+q.k*a_map[_n.i+2][_n.j+1][_n.k+1]+q.w*a_map[_n.i+2][_n.j+1][_n.k+2];
b_map[3][3]=q.i*a_map[_n.i+2][_n.j+2][_n.k-1]+q.j*a_map[_n.i+2][_n.j+2][_n.k-0]+q.k*a_map[_n.i+2][_n.j+2][_n.k+1]+q.w*a_map[_n.i+2][_n.j+2][_n.k+2]; //112 flops
calc_cubic_spline_findif(&q,_r.j);
b_map[0][0]=q.i*b_map[0][0]+q.j*b_map[0][1]+q.k*b_map[0][2]+q.w*b_map[0][3];
b_map[0][1]=q.i*b_map[1][0]+q.j*b_map[1][1]+q.k*b_map[1][2]+q.w*b_map[1][3];
b_map[0][2]=q.i*b_map[2][0]+q.j*b_map[2][1]+q.k*b_map[2][2]+q.w*b_map[2][3];
b_map[0][3]=q.i*b_map[3][0]+q.j*b_map[3][1]+q.k*b_map[3][2]+q.w*b_map[3][3]; // 28 flops
calc_cubic_spline_findif(&q,_r.i);
*f=.125*(q.i*b_map[0][0]+q.j*b_map[0][1]+q.k*b_map[0][2]+q.w*b_map[0][3]); //7 flops
return TRUE; // total 189 flops
}
//This function calculates tricubic interpolation and its derivative using finite differences instead of derivatives
char calc_tricubic_interpolation_derivative_dfindif(double *f,t_vec *g,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map)
{
t_quaternion qx, qy, qz;
t_vec _r;
t_len _n;
double b_map[4][4], c_map[4];

if ( ((_r.i=r->i-ori->i)<sp)||((_n.i=(int)(_r.i/sp))>ni-2)||((_r.j=r->j-ori->j)<sp)||((_n.j=(int)(_r.j/sp))>nj-2)||((_r.k=r->k-ori->k)<sp)||((_n.k=(int)(_r.k/sp))>nk-2) ) return FALSE;
else { _r.i=(_r.i-_n.i*sp)/sp, _r.j=(_r.j-_n.j*sp)/sp, _r.k=(_r.k-_n.k*sp)/sp; } // r is {0...1}
calc_cubic_spline_findif(&qz,_r.k);
b_map[0][0]=qz.i*a_map[_n.i-1][_n.j-1][_n.k-1]+qz.j*a_map[_n.i-1][_n.j-1][_n.k-0]+qz.k*a_map[_n.i-1][_n.j-1][_n.k+1]+qz.w*a_map[_n.i-1][_n.j-1][_n.k+2];
b_map[0][1]=qz.i*a_map[_n.i-1][_n.j-0][_n.k-1]+qz.j*a_map[_n.i-1][_n.j-0][_n.k-0]+qz.k*a_map[_n.i-1][_n.j-0][_n.k+1]+qz.w*a_map[_n.i-1][_n.j-0][_n.k+2];
b_map[0][2]=qz.i*a_map[_n.i-1][_n.j+1][_n.k-1]+qz.j*a_map[_n.i-1][_n.j+1][_n.k-0]+qz.k*a_map[_n.i-1][_n.j+1][_n.k+1]+qz.w*a_map[_n.i-1][_n.j+1][_n.k+2];
b_map[0][3]=qz.i*a_map[_n.i-1][_n.j+2][_n.k-1]+qz.j*a_map[_n.i-1][_n.j+2][_n.k-0]+qz.k*a_map[_n.i-1][_n.j+2][_n.k+1]+qz.w*a_map[_n.i-1][_n.j+2][_n.k+2];
b_map[1][0]=qz.i*a_map[_n.i-0][_n.j-1][_n.k-1]+qz.j*a_map[_n.i-0][_n.j-1][_n.k-0]+qz.k*a_map[_n.i-0][_n.j-1][_n.k+1]+qz.w*a_map[_n.i-0][_n.j-1][_n.k+2];
b_map[1][1]=qz.i*a_map[_n.i-0][_n.j-0][_n.k-1]+qz.j*a_map[_n.i-0][_n.j-0][_n.k-0]+qz.k*a_map[_n.i-0][_n.j-0][_n.k+1]+qz.w*a_map[_n.i-0][_n.j-0][_n.k+2];
b_map[1][2]=qz.i*a_map[_n.i-0][_n.j+1][_n.k-1]+qz.j*a_map[_n.i-0][_n.j+1][_n.k-0]+qz.k*a_map[_n.i-0][_n.j+1][_n.k+1]+qz.w*a_map[_n.i-0][_n.j+1][_n.k+2];
b_map[1][3]=qz.i*a_map[_n.i-0][_n.j+2][_n.k-1]+qz.j*a_map[_n.i-0][_n.j+2][_n.k-0]+qz.k*a_map[_n.i-0][_n.j+2][_n.k+1]+qz.w*a_map[_n.i-0][_n.j+2][_n.k+2];
b_map[2][0]=qz.i*a_map[_n.i+1][_n.j-1][_n.k-1]+qz.j*a_map[_n.i+1][_n.j-1][_n.k-0]+qz.k*a_map[_n.i+1][_n.j-1][_n.k+1]+qz.w*a_map[_n.i+1][_n.j-1][_n.k+2];
b_map[2][1]=qz.i*a_map[_n.i+1][_n.j-0][_n.k-1]+qz.j*a_map[_n.i+1][_n.j-0][_n.k-0]+qz.k*a_map[_n.i+1][_n.j-0][_n.k+1]+qz.w*a_map[_n.i+1][_n.j-0][_n.k+2];
b_map[2][2]=qz.i*a_map[_n.i+1][_n.j+1][_n.k-1]+qz.j*a_map[_n.i+1][_n.j+1][_n.k-0]+qz.k*a_map[_n.i+1][_n.j+1][_n.k+1]+qz.w*a_map[_n.i+1][_n.j+1][_n.k+2];
b_map[2][3]=qz.i*a_map[_n.i+1][_n.j+2][_n.k-1]+qz.j*a_map[_n.i+1][_n.j+2][_n.k-0]+qz.k*a_map[_n.i+1][_n.j+2][_n.k+1]+qz.w*a_map[_n.i+1][_n.j+2][_n.k+2];
b_map[3][0]=qz.i*a_map[_n.i+2][_n.j-1][_n.k-1]+qz.j*a_map[_n.i+2][_n.j-1][_n.k-0]+qz.k*a_map[_n.i+2][_n.j-1][_n.k+1]+qz.w*a_map[_n.i+2][_n.j-1][_n.k+2];
b_map[3][1]=qz.i*a_map[_n.i+2][_n.j-0][_n.k-1]+qz.j*a_map[_n.i+2][_n.j-0][_n.k-0]+qz.k*a_map[_n.i+2][_n.j-0][_n.k+1]+qz.w*a_map[_n.i+2][_n.j-0][_n.k+2];
b_map[3][2]=qz.i*a_map[_n.i+2][_n.j+1][_n.k-1]+qz.j*a_map[_n.i+2][_n.j+1][_n.k-0]+qz.k*a_map[_n.i+2][_n.j+1][_n.k+1]+qz.w*a_map[_n.i+2][_n.j+1][_n.k+2];
b_map[3][3]=qz.i*a_map[_n.i+2][_n.j+2][_n.k-1]+qz.j*a_map[_n.i+2][_n.j+2][_n.k-0]+qz.k*a_map[_n.i+2][_n.j+2][_n.k+1]+qz.w*a_map[_n.i+2][_n.j+2][_n.k+2]; //112 flops
calc_cubic_spline_findif(&qy,_r.j);
c_map[0]=qy.i*b_map[0][0]+qy.j*b_map[0][1]+qy.k*b_map[0][2]+qy.w*b_map[0][3];
c_map[1]=qy.i*b_map[1][0]+qy.j*b_map[1][1]+qy.k*b_map[1][2]+qy.w*b_map[1][3];
c_map[2]=qy.i*b_map[2][0]+qy.j*b_map[2][1]+qy.k*b_map[2][2]+qy.w*b_map[2][3];
c_map[3]=qy.i*b_map[3][0]+qy.j*b_map[3][1]+qy.k*b_map[3][2]+qy.w*b_map[3][3]; // 28 flops
calc_cubic_spline_findif(&qx,_r.i);
*f=.125*(qx.i*c_map[0]+qx.j*c_map[1]+qx.k*c_map[2]+qx.w*c_map[3]); //7 flops
calc_cubic_spline_dfindif(&qz,_r.i); //Here after qz is used as derivative
g->i=.125*(qz.i*c_map[0]+qz.j*c_map[1]+qz.k*c_map[2]+qz.w*c_map[3]); //7 flops
calc_cubic_spline_dfindif(&qz,_r.j);
c_map[0]=qz.i*b_map[0][0]+qz.j*b_map[0][1]+qz.k*b_map[0][2]+qz.w*b_map[0][3];
c_map[1]=qz.i*b_map[1][0]+qz.j*b_map[1][1]+qz.k*b_map[1][2]+qz.w*b_map[1][3];
c_map[2]=qz.i*b_map[2][0]+qz.j*b_map[2][1]+qz.k*b_map[2][2]+qz.w*b_map[2][3];
c_map[3]=qz.i*b_map[3][0]+qz.j*b_map[3][1]+qz.k*b_map[3][2]+qz.w*b_map[3][3]; //28 flops
g->j=.125*(qx.i*c_map[0]+qx.j*c_map[1]+qx.k*c_map[2]+qx.w*c_map[3]); //7 flops
calc_cubic_spline_dfindif(&qz,_r.k);
b_map[0][0]=qz.i*a_map[_n.i-1][_n.j-1][_n.k-1]+qz.j*a_map[_n.i-1][_n.j-1][_n.k-0]+qz.k*a_map[_n.i-1][_n.j-1][_n.k+1]+qz.w*a_map[_n.i-1][_n.j-1][_n.k+2];
b_map[0][1]=qz.i*a_map[_n.i-1][_n.j-0][_n.k-1]+qz.j*a_map[_n.i-1][_n.j-0][_n.k-0]+qz.k*a_map[_n.i-1][_n.j-0][_n.k+1]+qz.w*a_map[_n.i-1][_n.j-0][_n.k+2];
b_map[0][2]=qz.i*a_map[_n.i-1][_n.j+1][_n.k-1]+qz.j*a_map[_n.i-1][_n.j+1][_n.k-0]+qz.k*a_map[_n.i-1][_n.j+1][_n.k+1]+qz.w*a_map[_n.i-1][_n.j+1][_n.k+2];
b_map[0][3]=qz.i*a_map[_n.i-1][_n.j+2][_n.k-1]+qz.j*a_map[_n.i-1][_n.j+2][_n.k-0]+qz.k*a_map[_n.i-1][_n.j+2][_n.k+1]+qz.w*a_map[_n.i-1][_n.j+2][_n.k+2];
b_map[1][0]=qz.i*a_map[_n.i-0][_n.j-1][_n.k-1]+qz.j*a_map[_n.i-0][_n.j-1][_n.k-0]+qz.k*a_map[_n.i-0][_n.j-1][_n.k+1]+qz.w*a_map[_n.i-0][_n.j-1][_n.k+2];
b_map[1][1]=qz.i*a_map[_n.i-0][_n.j-0][_n.k-1]+qz.j*a_map[_n.i-0][_n.j-0][_n.k-0]+qz.k*a_map[_n.i-0][_n.j-0][_n.k+1]+qz.w*a_map[_n.i-0][_n.j-0][_n.k+2];
b_map[1][2]=qz.i*a_map[_n.i-0][_n.j+1][_n.k-1]+qz.j*a_map[_n.i-0][_n.j+1][_n.k-0]+qz.k*a_map[_n.i-0][_n.j+1][_n.k+1]+qz.w*a_map[_n.i-0][_n.j+1][_n.k+2];
b_map[1][3]=qz.i*a_map[_n.i-0][_n.j+2][_n.k-1]+qz.j*a_map[_n.i-0][_n.j+2][_n.k-0]+qz.k*a_map[_n.i-0][_n.j+2][_n.k+1]+qz.w*a_map[_n.i-0][_n.j+2][_n.k+2];
b_map[2][0]=qz.i*a_map[_n.i+1][_n.j-1][_n.k-1]+qz.j*a_map[_n.i+1][_n.j-1][_n.k-0]+qz.k*a_map[_n.i+1][_n.j-1][_n.k+1]+qz.w*a_map[_n.i+1][_n.j-1][_n.k+2];
b_map[2][1]=qz.i*a_map[_n.i+1][_n.j-0][_n.k-1]+qz.j*a_map[_n.i+1][_n.j-0][_n.k-0]+qz.k*a_map[_n.i+1][_n.j-0][_n.k+1]+qz.w*a_map[_n.i+1][_n.j-0][_n.k+2];
b_map[2][2]=qz.i*a_map[_n.i+1][_n.j+1][_n.k-1]+qz.j*a_map[_n.i+1][_n.j+1][_n.k-0]+qz.k*a_map[_n.i+1][_n.j+1][_n.k+1]+qz.w*a_map[_n.i+1][_n.j+1][_n.k+2];
b_map[2][3]=qz.i*a_map[_n.i+1][_n.j+2][_n.k-1]+qz.j*a_map[_n.i+1][_n.j+2][_n.k-0]+qz.k*a_map[_n.i+1][_n.j+2][_n.k+1]+qz.w*a_map[_n.i+1][_n.j+2][_n.k+2];
b_map[3][0]=qz.i*a_map[_n.i+2][_n.j-1][_n.k-1]+qz.j*a_map[_n.i+2][_n.j-1][_n.k-0]+qz.k*a_map[_n.i+2][_n.j-1][_n.k+1]+qz.w*a_map[_n.i+2][_n.j-1][_n.k+2];
b_map[3][1]=qz.i*a_map[_n.i+2][_n.j-0][_n.k-1]+qz.j*a_map[_n.i+2][_n.j-0][_n.k-0]+qz.k*a_map[_n.i+2][_n.j-0][_n.k+1]+qz.w*a_map[_n.i+2][_n.j-0][_n.k+2];
b_map[3][2]=qz.i*a_map[_n.i+2][_n.j+1][_n.k-1]+qz.j*a_map[_n.i+2][_n.j+1][_n.k-0]+qz.k*a_map[_n.i+2][_n.j+1][_n.k+1]+qz.w*a_map[_n.i+2][_n.j+1][_n.k+2];
b_map[3][3]=qz.i*a_map[_n.i+2][_n.j+2][_n.k-1]+qz.j*a_map[_n.i+2][_n.j+2][_n.k-0]+qz.k*a_map[_n.i+2][_n.j+2][_n.k+1]+qz.w*a_map[_n.i+2][_n.j+2][_n.k+2]; //112 flops
c_map[0]=qy.i*b_map[0][0]+qy.j*b_map[0][1]+qy.k*b_map[0][2]+qy.w*b_map[0][3];
c_map[1]=qy.i*b_map[1][0]+qy.j*b_map[1][1]+qy.k*b_map[1][2]+qy.w*b_map[1][3];
c_map[2]=qy.i*b_map[2][0]+qy.j*b_map[2][1]+qy.k*b_map[2][2]+qy.w*b_map[2][3];
c_map[3]=qy.i*b_map[3][0]+qy.j*b_map[3][1]+qy.k*b_map[3][2]+qy.w*b_map[3][3]; //28 flops
g->k=.125*(qx.i*c_map[0]+qx.j*c_map[1]+qx.k*c_map[2]+qx.w*c_map[3]); //7 flops
return TRUE; // total 460 flops
}

//--------- M O N O T O N I C Y - P R E S E R V I N G    P A R T ------------------//

//Hyman J.M. Accurate monotonicity preserving cubic interpolation SIAM J.Sci.Stat.Comput., 4, 645-654, 1983
//Fritsch F.N. and Carlson R.E. Monotonicity preserving bicubic interpolation: A progress report Comp.Aid.Geom.Des., 2, 117-121, 1985 
//Carlson R.E. and Fritsch F.N. Monotone piecewise bicubic interpolation SIAM J.Numer.Anal., 22, 386-400, 1985 

//           
//           x^    ^z
//            |   /
//            |  /                   |<gx1>|
//         f2 * * f3              |<gx0>|  |
//      f7 f4 |/ f1               |  |  |  |
//    --*--*--*--*----> x       --*--*--*--*---->
//           /|f0                    |__|f0|       
//       f6 * * f5                    s0|__|
//         /  |                          s1 
//     f9 *   * f8
//       /    |  
//

//This function calculate finite difference derivatives
inline double calc_fdif_derivative_preserving_monotonocity(double f0,double f1,double f2)
{
register double s0, s1;
s0=f1-f0, s1=f2-f1;
if ((s0==0.)&&(s1==0.)) return 0.;
else return (3.*s1*s0/(s1*s1+s0*s0+s1*s0))*(f2-f0)/2.;
}
inline void calc_derivative_fdif_derivative_preserving_monotonocy(double dg[3],double f0,double f1,double f2)
{
register double _denominator;
_denominator=2.*sqrd(sqrd(f2-f1)+sqrd(f1-f0)+(f2-f1)*(f1-f0))/3.;
if (_denominator==0.) { dg[0]=dg[1]=dg[2]=0.; }
else
  {
  dg[0]=+(2.*f0-f1-f2)*cubed(f2-f1)/_denominator,
  dg[1]=+(f0-2.*f1+f2)*cubed(f2-f0)/_denominator,
  dg[2]=-(f0+f1-2.*f2)*cubed(f1-f0)/_denominator;
  }
}
//This function calculates polynome constants
inline void calc_polinomial_constants(double *c0,double *c1,double *c2,double *c3,double f0,double f1,double f2,double f3)
{
register double s, df0, df1;
s=f2-f1, df0=calc_fdif_derivative_preserving_monotonocity(f0,f1,f2), df1=calc_fdif_derivative_preserving_monotonocity(f1,f2,f3);
*c0=f1, *c1=df0, *c2=3.*s-df1-2.*df0, *c3=df1+df0-2.*s;
}
//This function calculates derivatives of polinomial constants over f values
inline void calc_polinomial_constants_derivative(double dc0[4],double dc1[4],double dc2[4],double dc3[4],double f0,double f1,double f2,double f3)
{
double dg1[3], dg2[3];
//Calc derivatives of derivatives
calc_derivative_fdif_derivative_preserving_monotonocy(dg1,f0,f1,f2);
calc_derivative_fdif_derivative_preserving_monotonocy(dg2,f1,f2,f3);
//Calculate dc0/df
dc0[0]=0, dc0[1]=1., dc0[2]=0., dc0[3]=0.; 
//Calculate dc1/df
dc1[0]=dg1[0], dc1[1]=dg1[1], dc1[2]=dg1[2], dc1[3]=0.;
//Calculate dc2/df
dc2[0]=-2.*dg1[0], dc2[1]=-3.-2.*dg1[1]-dg2[0], dc2[2]=3.-2.*dg1[2]-dg2[1], dc2[3]=-dg2[2];
//Calculate dc3/df
dc3[0]=dg1[0], dc3[1]=2.+dg1[1]+dg2[0], dc3[2]=-2.+dg1[2]+dg2[1], dc3[3]=dg2[2];
}
//This function interpolates value of a unit cubic polynome (x is 0...1)
inline double interpolate_cubic_polynome(double x,double c0,double c1,double c2,double c3)
{
return c0+x*(c1+x*(c2+c3*x));
}
//This function calculates derivatives of polynome constants as a function of x
inline double calc_derivative_of_cubic_polynome(double x,double c0,double c1,double c2,double c3)
{
return c1+x*(2.*c2+3.*c3*x);
}

//This function calculates tricubic polinome values with satisfing monotonicity and sign-preserving inequalities
char calc_tricubic_interpolation_wp_monotonicity(double *f,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map)
{
t_vec _r;
t_len _n;
double c0, c1, c2, c3, b_map[4][4], c_map[4];

if ( ((_r.i=r->i-ori->i)<sp)||((_n.i=(int)(_r.i/sp))>ni-3)||((_r.j=r->j-ori->j)<sp)||((_n.j=(int)(_r.j/sp))>nj-3)||((_r.k=r->k-ori->k)<sp)||((_n.k=(int)(_r.k/sp))>nk-3) ) return FALSE;
else { _r.i=(_r.i-_n.i*sp)/sp, _r.j=(_r.j-_n.j*sp)/sp, _r.k=(_r.k-_n.k*sp)/sp; } // r is {0...1}
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-1][_n.j-1][_n.k-1],a_map[_n.i-1][_n.j-1][_n.k-0],a_map[_n.i-1][_n.j-1][_n.k+1],a_map[_n.i-1][_n.j-1][_n.k+2]), b_map[0][0]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-1][_n.j-0][_n.k-1],a_map[_n.i-1][_n.j-0][_n.k-0],a_map[_n.i-1][_n.j-0][_n.k+1],a_map[_n.i-1][_n.j-0][_n.k+2]), b_map[0][1]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-1][_n.j+1][_n.k-1],a_map[_n.i-1][_n.j+1][_n.k-0],a_map[_n.i-1][_n.j+1][_n.k+1],a_map[_n.i-1][_n.j+1][_n.k+2]), b_map[0][2]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-1][_n.j+2][_n.k-1],a_map[_n.i-1][_n.j+2][_n.k-0],a_map[_n.i-1][_n.j+2][_n.k+1],a_map[_n.i-1][_n.j+2][_n.k+2]), b_map[0][3]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-0][_n.j-1][_n.k-1],a_map[_n.i-0][_n.j-1][_n.k-0],a_map[_n.i-0][_n.j-1][_n.k+1],a_map[_n.i-0][_n.j-1][_n.k+2]), b_map[1][0]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-0][_n.j-0][_n.k-1],a_map[_n.i-0][_n.j-0][_n.k-0],a_map[_n.i-0][_n.j-0][_n.k+1],a_map[_n.i-0][_n.j-0][_n.k+2]), b_map[1][1]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-0][_n.j+1][_n.k-1],a_map[_n.i-0][_n.j+1][_n.k-0],a_map[_n.i-0][_n.j+1][_n.k+1],a_map[_n.i-0][_n.j+1][_n.k+2]), b_map[1][2]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-0][_n.j+2][_n.k-1],a_map[_n.i-0][_n.j+2][_n.k-0],a_map[_n.i-0][_n.j+2][_n.k+1],a_map[_n.i-0][_n.j+2][_n.k+2]), b_map[1][3]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+1][_n.j-1][_n.k-1],a_map[_n.i+1][_n.j-1][_n.k-0],a_map[_n.i+1][_n.j-1][_n.k+1],a_map[_n.i+1][_n.j-1][_n.k+2]), b_map[2][0]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+1][_n.j-0][_n.k-1],a_map[_n.i+1][_n.j-0][_n.k-0],a_map[_n.i+1][_n.j-0][_n.k+1],a_map[_n.i+1][_n.j-0][_n.k+2]), b_map[2][1]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+1][_n.j+1][_n.k-1],a_map[_n.i+1][_n.j+1][_n.k-0],a_map[_n.i+1][_n.j+1][_n.k+1],a_map[_n.i+1][_n.j+1][_n.k+2]), b_map[2][2]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+1][_n.j+2][_n.k-1],a_map[_n.i+1][_n.j+2][_n.k-0],a_map[_n.i+1][_n.j+2][_n.k+1],a_map[_n.i+1][_n.j+2][_n.k+2]), b_map[2][3]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+2][_n.j-1][_n.k-1],a_map[_n.i+2][_n.j-1][_n.k-0],a_map[_n.i+2][_n.j-1][_n.k+1],a_map[_n.i+2][_n.j-1][_n.k+2]), b_map[3][0]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+2][_n.j-0][_n.k-1],a_map[_n.i+2][_n.j-0][_n.k-0],a_map[_n.i+2][_n.j-0][_n.k+1],a_map[_n.i+2][_n.j-0][_n.k+2]), b_map[3][1]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+2][_n.j+1][_n.k-1],a_map[_n.i+2][_n.j+1][_n.k-0],a_map[_n.i+2][_n.j+1][_n.k+1],a_map[_n.i+2][_n.j+1][_n.k+2]), b_map[3][2]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+2][_n.j+2][_n.k-1],a_map[_n.i+2][_n.j+2][_n.k-0],a_map[_n.i+2][_n.j+2][_n.k+1],a_map[_n.i+2][_n.j+2][_n.k+2]), b_map[3][3]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,b_map[0][0],b_map[0][1],b_map[0][2],b_map[0][3]), c_map[0]=interpolate_cubic_polynome(_r.j,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,b_map[1][0],b_map[1][1],b_map[1][2],b_map[1][3]), c_map[1]=interpolate_cubic_polynome(_r.j,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,b_map[2][0],b_map[2][1],b_map[2][2],b_map[2][3]), c_map[2]=interpolate_cubic_polynome(_r.j,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,b_map[3][0],b_map[3][1],b_map[3][2],b_map[3][3]), c_map[3]=interpolate_cubic_polynome(_r.j,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,c_map[0],c_map[1],c_map[2],c_map[3]), *f=interpolate_cubic_polynome(_r.i,c0,c1,c2,c3);
return TRUE;
}
//This function calculates tricubic polinome values with satisfing monotonicity and sign-preserving inequalities
char calc_tricubic_interpolation_derivative_wp_monotonicity(double *f,t_vec *g,t_vec *ori,double sp,unsigned int ni,unsigned int nj,unsigned int nk,t_vec *r,double ***a_map)
{
t_vec _r;
t_len _n;
double c0, c1, c2, c3, dc0[4], dc1[4], dc2[4], dc3[4], b_map[4][4], db_map[4][4], c_map[4], dc_map[4], dd_map[4];

if ( ((_r.i=r->i-ori->i)<sp)||((_n.i=(int)(_r.i/sp))>ni-3)||((_r.j=r->j-ori->j)<sp)||((_n.j=(int)(_r.j/sp))>nj-3)||((_r.k=r->k-ori->k)<sp)||((_n.k=(int)(_r.k/sp))>nk-3) ) return FALSE;
else { _r.i=(_r.i-_n.i*sp)/sp, _r.j=(_r.j-_n.j*sp)/sp, _r.k=(_r.k-_n.k*sp)/sp; } // r is {0...1}
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-1][_n.j-1][_n.k-1],a_map[_n.i-1][_n.j-1][_n.k-0],a_map[_n.i-1][_n.j-1][_n.k+1],a_map[_n.i-1][_n.j-1][_n.k+2]), b_map[0][0]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[0][0]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-1][_n.j-0][_n.k-1],a_map[_n.i-1][_n.j-0][_n.k-0],a_map[_n.i-1][_n.j-0][_n.k+1],a_map[_n.i-1][_n.j-0][_n.k+2]), b_map[0][1]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[0][1]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-1][_n.j+1][_n.k-1],a_map[_n.i-1][_n.j+1][_n.k-0],a_map[_n.i-1][_n.j+1][_n.k+1],a_map[_n.i-1][_n.j+1][_n.k+2]), b_map[0][2]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[0][2]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-1][_n.j+2][_n.k-1],a_map[_n.i-1][_n.j+2][_n.k-0],a_map[_n.i-1][_n.j+2][_n.k+1],a_map[_n.i-1][_n.j+2][_n.k+2]), b_map[0][3]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[0][3]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-0][_n.j-1][_n.k-1],a_map[_n.i-0][_n.j-1][_n.k-0],a_map[_n.i-0][_n.j-1][_n.k+1],a_map[_n.i-0][_n.j-1][_n.k+2]), b_map[1][0]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[1][0]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-0][_n.j-0][_n.k-1],a_map[_n.i-0][_n.j-0][_n.k-0],a_map[_n.i-0][_n.j-0][_n.k+1],a_map[_n.i-0][_n.j-0][_n.k+2]), b_map[1][1]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[1][1]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-0][_n.j+1][_n.k-1],a_map[_n.i-0][_n.j+1][_n.k-0],a_map[_n.i-0][_n.j+1][_n.k+1],a_map[_n.i-0][_n.j+1][_n.k+2]), b_map[1][2]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[1][2]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i-0][_n.j+2][_n.k-1],a_map[_n.i-0][_n.j+2][_n.k-0],a_map[_n.i-0][_n.j+2][_n.k+1],a_map[_n.i-0][_n.j+2][_n.k+2]), b_map[1][3]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[1][3]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+1][_n.j-1][_n.k-1],a_map[_n.i+1][_n.j-1][_n.k-0],a_map[_n.i+1][_n.j-1][_n.k+1],a_map[_n.i+1][_n.j-1][_n.k+2]), b_map[2][0]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[2][0]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+1][_n.j-0][_n.k-1],a_map[_n.i+1][_n.j-0][_n.k-0],a_map[_n.i+1][_n.j-0][_n.k+1],a_map[_n.i+1][_n.j-0][_n.k+2]), b_map[2][1]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[2][1]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+1][_n.j+1][_n.k-1],a_map[_n.i+1][_n.j+1][_n.k-0],a_map[_n.i+1][_n.j+1][_n.k+1],a_map[_n.i+1][_n.j+1][_n.k+2]), b_map[2][2]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[2][2]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+1][_n.j+2][_n.k-1],a_map[_n.i+1][_n.j+2][_n.k-0],a_map[_n.i+1][_n.j+2][_n.k+1],a_map[_n.i+1][_n.j+2][_n.k+2]), b_map[2][3]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[2][3]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+2][_n.j-1][_n.k-1],a_map[_n.i+2][_n.j-1][_n.k-0],a_map[_n.i+2][_n.j-1][_n.k+1],a_map[_n.i+2][_n.j-1][_n.k+2]), b_map[3][0]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[3][0]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+2][_n.j-0][_n.k-1],a_map[_n.i+2][_n.j-0][_n.k-0],a_map[_n.i+2][_n.j-0][_n.k+1],a_map[_n.i+2][_n.j-0][_n.k+2]), b_map[3][1]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[3][1]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+2][_n.j+1][_n.k-1],a_map[_n.i+2][_n.j+1][_n.k-0],a_map[_n.i+2][_n.j+1][_n.k+1],a_map[_n.i+2][_n.j+1][_n.k+2]), b_map[3][2]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[3][2]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
calc_polinomial_constants(&c0,&c1,&c2,&c3,a_map[_n.i+2][_n.j+2][_n.k-1],a_map[_n.i+2][_n.j+2][_n.k-0],a_map[_n.i+2][_n.j+2][_n.k+1],a_map[_n.i+2][_n.j+2][_n.k+2]), b_map[3][3]=interpolate_cubic_polynome(_r.k,c0,c1,c2,c3), db_map[3][3]=calc_derivative_of_cubic_polynome(_r.k,c0,c1,c2,c3);
//Calc dfx/dz
calc_polinomial_constants(&c0,&c1,&c2,&c3,b_map[0][0],b_map[0][1],b_map[0][2],b_map[0][3]), c_map[0]=interpolate_cubic_polynome(_r.j,c0,c1,c2,c3), dc_map[0]=calc_derivative_of_cubic_polynome(_r.j,c0,c1,c2,c3);
calc_polinomial_constants_derivative(dc0,dc1,dc2,dc3,b_map[0][0],b_map[0][1],b_map[0][2],b_map[0][3]);
dd_map[0]=interpolate_cubic_polynome(_r.j,dc0[0]*db_map[0][0]+dc0[1]*db_map[0][1]+dc0[2]*db_map[0][2]+dc0[3]*db_map[0][3],dc1[0]*db_map[0][0]+dc1[1]*db_map[0][1]+dc1[2]*db_map[0][2]+dc1[3]*db_map[0][3],dc2[0]*db_map[0][0]+dc2[1]*db_map[0][1]+dc2[2]*db_map[0][2]+dc2[3]*db_map[0][3],dc3[0]*db_map[0][0]+dc3[1]*db_map[0][1]+dc3[2]*db_map[0][2]+dc3[3]*db_map[0][3]);
calc_polinomial_constants(&c0,&c1,&c2,&c3,b_map[1][0],b_map[1][1],b_map[1][2],b_map[1][3]), c_map[1]=interpolate_cubic_polynome(_r.j,c0,c1,c2,c3), dc_map[1]=calc_derivative_of_cubic_polynome(_r.j,c0,c1,c2,c3);
calc_polinomial_constants_derivative(dc0,dc1,dc2,dc3,b_map[1][0],b_map[1][1],b_map[1][2],b_map[1][3]);
dd_map[1]=interpolate_cubic_polynome(_r.j,dc0[0]*db_map[1][0]+dc0[1]*db_map[1][1]+dc0[2]*db_map[1][2]+dc0[3]*db_map[1][3],dc1[0]*db_map[1][0]+dc1[1]*db_map[1][1]+dc1[2]*db_map[1][2]+dc1[3]*db_map[1][3],dc2[0]*db_map[1][0]+dc2[1]*db_map[1][1]+dc2[2]*db_map[1][2]+dc2[3]*db_map[1][3],dc3[0]*db_map[1][0]+dc3[1]*db_map[1][1]+dc3[2]*db_map[1][2]+dc3[3]*db_map[1][3]);
calc_polinomial_constants(&c0,&c1,&c2,&c3,b_map[2][0],b_map[2][1],b_map[2][2],b_map[2][3]), c_map[2]=interpolate_cubic_polynome(_r.j,c0,c1,c2,c3), dc_map[2]=calc_derivative_of_cubic_polynome(_r.j,c0,c1,c2,c3);
calc_polinomial_constants_derivative(dc0,dc1,dc2,dc3,b_map[2][0],b_map[2][1],b_map[2][2],b_map[2][3]);
dd_map[2]=interpolate_cubic_polynome(_r.j,dc0[0]*db_map[2][0]+dc0[1]*db_map[2][1]+dc0[2]*db_map[2][2]+dc0[3]*db_map[2][3],dc1[0]*db_map[2][0]+dc1[1]*db_map[2][1]+dc1[2]*db_map[2][2]+dc1[3]*db_map[2][3],dc2[0]*db_map[2][0]+dc2[1]*db_map[2][1]+dc2[2]*db_map[2][2]+dc2[3]*db_map[2][3],dc3[0]*db_map[2][0]+dc3[1]*db_map[2][1]+dc3[2]*db_map[2][2]+dc3[3]*db_map[2][3]);
calc_polinomial_constants(&c0,&c1,&c2,&c3,b_map[3][0],b_map[3][1],b_map[3][2],b_map[3][3]), c_map[3]=interpolate_cubic_polynome(_r.j,c0,c1,c2,c3), dc_map[3]=calc_derivative_of_cubic_polynome(_r.j,c0,c1,c2,c3);
calc_polinomial_constants_derivative(dc0,dc1,dc2,dc3,b_map[3][0],b_map[3][1],b_map[3][2],b_map[3][3]);
dd_map[3]=interpolate_cubic_polynome(_r.j,dc0[0]*db_map[3][0]+dc0[1]*db_map[3][1]+dc0[2]*db_map[3][2]+dc0[3]*db_map[3][3],dc1[0]*db_map[3][0]+dc1[1]*db_map[3][1]+dc1[2]*db_map[3][2]+dc1[3]*db_map[3][3],dc2[0]*db_map[3][0]+dc2[1]*db_map[3][1]+dc2[2]*db_map[3][2]+dc2[3]*db_map[3][3],dc3[0]*db_map[3][0]+dc3[1]*db_map[3][1]+dc3[2]*db_map[3][2]+dc3[3]*db_map[3][3]);
//Calc p
calc_polinomial_constants(&c0,&c1,&c2,&c3,c_map[0],c_map[1],c_map[2],c_map[3]), *f=interpolate_cubic_polynome(_r.i,c0,c1,c2,c3);
//Calc dcx/dfy
calc_polinomial_constants_derivative(dc0,dc1,dc2,dc3,c_map[0],c_map[1],c_map[2],c_map[3]);
//Calc derivatives
g->i=calc_derivative_of_cubic_polynome(_r.i,c0,c1,c2,c3);
g->j=interpolate_cubic_polynome(_r.i,dc0[0]*dc_map[0]+dc0[1]*dc_map[1]+dc0[2]*dc_map[2]+dc0[3]*dc_map[3],dc1[0]*dc_map[0]+dc1[1]*dc_map[1]+dc1[2]*dc_map[2]+dc1[3]*dc_map[3],dc2[0]*dc_map[0]+dc2[1]*dc_map[1]+dc2[2]*dc_map[2]+dc2[3]*dc_map[3],dc3[0]*dc_map[0]+dc3[1]*dc_map[1]+dc3[2]*dc_map[2]+dc3[3]*dc_map[3]);
g->k=interpolate_cubic_polynome(_r.i,dc0[0]*dd_map[0]+dc0[1]*dd_map[1]+dc0[2]*dd_map[2]+dc0[3]*dd_map[3],dc1[0]*dd_map[0]+dc1[1]*dd_map[1]+dc1[2]*dd_map[2]+dc1[3]*dd_map[3],dc2[0]*dd_map[0]+dc2[1]*dd_map[1]+dc2[2]*dd_map[2]+dc2[3]*dd_map[3],dc3[0]*dd_map[0]+dc3[1]*dd_map[1]+dc3[2]*dd_map[2]+dc3[3]*dd_map[3]);
g->i/=sp, g->j/=sp, g->k/=sp; //Calculate derivative in external units
return TRUE;
}


//------------------------------------------------------ O B J E C T S     R O U N D I N G -------------------------------------------------//

//This function numerically perturbates a tensor so its determinant become |Det(T)^2-1.| < tol
unsigned int _round_udet_tensor(double *s,unsigned int n,double *x,double *g,double **G,va_list stack)
{
register double _d;
t_tensor *xR=(t_tensor*)x, *gR=(t_tensor*)g;

*s=calc_det3x3((*xR)[0][0],(*xR)[0][1],(*xR)[0][2],(*xR)[1][0],(*xR)[1][1],(*xR)[1][2],(*xR)[2][0],(*xR)[2][1],(*xR)[2][2]);
_d=4.**s*(*s**s-1.);
(*gR)[0][0]=_d*calc_det2x2((*xR)[1][1],(*xR)[1][2],(*xR)[2][1],(*xR)[2][2]);
(*gR)[0][1]=_d*calc_det2x2((*xR)[1][2],(*xR)[1][0],(*xR)[2][2],(*xR)[2][0]);
(*gR)[0][2]=_d*calc_det2x2((*xR)[1][0],(*xR)[1][1],(*xR)[2][0],(*xR)[2][1]);
(*gR)[1][0]=_d*calc_det2x2((*xR)[2][1],(*xR)[2][2],(*xR)[0][1],(*xR)[0][2]);
(*gR)[1][1]=_d*calc_det2x2((*xR)[2][2],(*xR)[2][0],(*xR)[0][2],(*xR)[0][0]);
(*gR)[1][2]=_d*calc_det2x2((*xR)[2][0],(*xR)[2][1],(*xR)[0][0],(*xR)[0][1]);
(*gR)[2][0]=_d*calc_det2x2((*xR)[0][1],(*xR)[0][2],(*xR)[1][1],(*xR)[1][2]);
(*gR)[2][1]=_d*calc_det2x2((*xR)[0][2],(*xR)[0][0],(*xR)[1][2],(*xR)[1][0]);
(*gR)[2][2]=_d*calc_det2x2((*xR)[0][0],(*xR)[0][1],(*xR)[1][0],(*xR)[1][1]);

*s*=*s, *s-=1., *s*=*s;
return YERROR_OK;
}
void round_udet_tensor(t_tensor *R,double tol)
{
double s, *x[2], *g[2], p[9*5];

x[0]=&p[9], x[1]=&p[18], g[0]=&p[27], g[1]=&p[36];
x[0][0]=(*R)[0][0], x[0][1]=(*R)[0][1], x[0][2]=(*R)[0][2],
x[0][3]=(*R)[1][0], x[0][4]=(*R)[1][1], x[0][5]=(*R)[1][2],
x[0][6]=(*R)[2][0], x[0][7]=(*R)[2][1], x[0][8]=(*R)[2][2];
_round_udet_tensor(&s,9,x[0],g[0],0x0,0x0);
if (s>tol*tol)
  {
  if (calc_vect_norm(9,g[0])<SMALL2)
    {//No derivatives - turn into unit tensor
    (*R)[0][0]=1., (*R)[0][1]=0., (*R)[0][2]=0.,
    (*R)[1][0]=0., (*R)[1][1]=1., (*R)[1][2]=0.,
    (*R)[2][0]=0., (*R)[2][1]=0., (*R)[2][2]=1.;
    }
  else
    {//Minimize given tensor into unit one
    polak_ribiere(&s,(unsigned int)-1,tol,0.1*SMALL2,4.*sqrt(s),9,x,g,0x0,p,_round_udet_tensor,line_search_square_fapproximation,0x0);
    (*R)[0][0]=x[0][0], (*R)[0][1]=x[0][1], (*R)[0][2]=x[0][2],
    (*R)[1][0]=x[0][3], (*R)[1][1]=x[0][4], (*R)[1][2]=x[0][5],
    (*R)[2][0]=x[0][6], (*R)[2][1]=x[0][7], (*R)[2][2]=x[0][8];
    }
  }
}




