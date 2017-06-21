#include "y_mbuild.h"

extern unsigned int ylib_errno;

//RESOLVERS

//1-resolver
void _resolve_one_monovalent_at_one_coordinated_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
r[group[1]].i=r[*group].i, r[group[1]].j=r[*group].j, r[group[1]].k=r[*group].k+v[0][1];  
}
//2-resolver
void _resolve_two_monovalent_at_two_coordinated_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
r[group[1]].i=r[*group].i,                      r[group[1]].j=r[*group].j, r[group[1]].k=r[*group].k+v[0][1];  
r[group[2]].i=r[*group].i+v[0][2]*sin(v[1][2]), r[group[2]].j=r[*group].j, r[group[2]].k=r[*group].k+v[0][2]*cos(v[1][2]);  
}
//This (help) function takes an orthogonal vector with smallest projection onto x, y or z (global) coordinate axis
void _get_orthogonal_vector(t_vec *p,t_vec *u)
{
//Chose the leading axis: min(u.i*i_+u.j*j_+u.k*k_)
if (fabs(u->i)<fabs(u->j))
  {
  if (fabs(u->i)<fabs(u->k))
    {//x is the leader
    p->i=0., p->j=-u->k, p->k=+u->j; //p = u x x = ( 0.; -u.k; +u.j )
    }
  else
    {//z is the leader
    LABEL_Z_LEADER: p->i=-u->j, p->j=+u->i, p->k=0.; //p = u x z = ( -u.j; +u.i; 0. )
    }
  }
else
  {
  if (fabs(u->j)<fabs(u->k))
    {//y is the leader
    p->i=+u->k, p->j=0., p->k=-u->i; //p = u x y = ( u.k; 0.; -u.i )
    }
  else goto LABEL_Z_LEADER;
  }
}
void _resolve_one_monovalent_at_two_coordinated_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
t_vec u, p;
u.i=r[*group].i-r[group[1]].i, u.j=r[*group].j-r[group[1]].j, u.k=r[*group].k-r[group[1]].k; multiple_vec_scalar(&u,&u,1./sqrt(calc_vec_norm(&u)));
_get_orthogonal_vector(&p,&u); multiple_vec_scalar(&p,&p,1./sqrt(calc_vec_norm(&p)));
r[group[2]].i=r[*group].i+u.i*v[0][2]*cos(v[1][2])+p.i*v[0][2]*sin(v[1][2]), r[group[2]].j=r[*group].j+u.j*v[0][2]*cos(v[1][2])+p.j*v[0][2]*sin(v[1][2]), r[group[2]].k=r[*group].k+u.k*v[0][2]*cos(v[1][2])+p.k*v[0][2]*sin(v[1][2]); //r=r0+u/|u|*d*cosT+p/|p|*d*sinT
}
//3-resolvers (planar)
void _resolve_three_monovalent_at_three_coordinated_planar_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
r[group[1]].i=r[*group].i,                       r[group[1]].j=r[*group].j, r[group[1]].k=r[*group].k+v[0][1];  
r[group[2]].i=r[*group].i+v[0][2]*sin(+v[1][2]), r[group[2]].j=r[*group].j, r[group[2]].k=r[*group].k+v[0][2]*cos(+v[1][2]);  
r[group[3]].i=r[*group].i+v[0][3]*sin(-v[1][3]), r[group[3]].j=r[*group].j, r[group[3]].k=r[*group].k+v[0][3]*cos(-v[1][3]);  
}
void _resolve_two_monovalent_at_three_coordinated_planar_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
t_vec u, p;
u.i=r[*group].i-r[group[1]].i, u.j=r[*group].j-r[group[1]].j, u.k=r[*group].k-r[group[1]].k; multiple_vec_scalar(&u,&u,1./sqrt(calc_vec_norm(&u)));
_get_orthogonal_vector(&p,&u); multiple_vec_scalar(&p,&p,1./sqrt(calc_vec_norm(&p)));
r[group[2]].i=r[*group].i+u.i*v[0][2]*cos(+v[1][2])+p.i*v[0][2]*sin(+v[1][2]), r[group[2]].j=r[*group].j+u.j*v[0][2]*cos(v[1][2])+p.j*v[0][2]*sin(v[1][2]), r[group[2]].k=r[*group].k+u.k*v[0][2]*cos(+v[1][2])+p.k*v[0][2]*sin(+v[1][2]);
r[group[3]].i=r[*group].i+u.i*v[0][3]*cos(-v[1][3])+p.i*v[0][3]*sin(-v[1][3]), r[group[3]].j=r[*group].j+u.j*v[0][3]*cos(v[1][3])+p.j*v[0][3]*sin(v[1][3]), r[group[3]].k=r[*group].k+u.k*v[0][3]*cos(-v[1][3])+p.k*v[0][3]*sin(-v[1][3]); //r=r0+u/|u|*d*cosT+p/|p|*d*sinT
}
void _resolve_one_monovalent_at_three_coordinated_planar_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
t_vec u;
u.i=r[group[1]].i+r[group[2]].i-2.*r[*group].i, u.j=r[group[1]].j+r[group[2]].j-2.*r[*group].j, u.k=r[group[1]].k+r[group[2]].k-2.*r[*group].k;
multiple_vec_scalar(&u,&u,v[0][3]/sqrt(calc_vec_norm(&u)));
r[group[3]].i=r[*group].i-u.i, r[group[3]].j=r[*group].j-u.j, r[group[3]].k=r[*group].k-u.k;
}
//3-resolvers (nonplanar)
void _resolve_three_monovalent_at_three_coordinated_nonplanar_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
_resolve_three_monovalent_at_three_coordinated_planar_center(n,m,v,group,r);
}
void _resolve_two_monovalent_at_three_coordinated_nonplanar_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
_resolve_two_monovalent_at_three_coordinated_planar_center(n,m,v,group,r);
}
void _resolve_one_monovalent_at_three_coordinated_nonplanar_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
_resolve_one_monovalent_at_three_coordinated_planar_center(n,m,v,group,r);
}
//4-resolvers
void _resolve_four_monovalent_at_four_coordinated_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
t_vec u, p, q;
u.i=0., u.j=0., u.k=1.; 
p.i=0., p.j=1., p.k=0.;
q.i=1., q.j=0., q.k=0.;
r[group[1]].i=r[*group].i+0., r[group[1]].j=r[*group].j+0., r[group[1]].k=r[*group].k-v[0][1]; 
r[group[2]].i=r[*group].i+u.i*v[0][2]*cos(v[1][2])+q.i*v[0][2]*sin(v[1][2]),                                              r[group[2]].j=r[*group].j+u.j*v[0][2]*cos(v[1][2])+q.j*v[0][2]*sin(v[1][2]),                                             r[group[2]].k=r[*group].k+u.k*v[0][2]*cos(v[1][2])+q.k*v[0][2]*sin(v[1][2]);
r[group[3]].i=r[*group].i+u.i*v[0][3]*cos(v[1][3])-q.i*(1./2.)*v[0][3]*sin(v[1][3])+p.i*(SQRT_3/2.)*v[0][3]*sin(v[1][3]), r[group[3]].j=r[*group].j+u.j*v[0][3]*cos(v[1][3])-q.j*(1./2.)*v[0][3]*sin(v[1][3])+p.j*(SQRT_3/2.)*v[0][3]*sin(v[1][3]), r[group[3]].k=r[*group].k+u.k*v[0][3]*cos(v[1][3])-q.k*(1./2.)*v[0][3]*sin(v[1][3])+p.k*(SQRT_3/2.)*v[0][3]*sin(v[1][3]);
r[group[4]].i=r[*group].i+u.i*v[0][4]*cos(v[1][4])-q.i*(1./2.)*v[0][4]*sin(v[1][4])-p.i*(SQRT_3/2.)*v[0][4]*sin(v[1][4]), r[group[4]].j=r[*group].j+u.j*v[0][4]*cos(v[1][4])-q.j*(1./2.)*v[0][4]*sin(v[1][4])-p.j*(SQRT_3/2.)*v[0][4]*sin(v[1][4]), r[group[4]].k=r[*group].k+u.k*v[0][4]*cos(v[1][4])-q.k*(1./2.)*v[0][4]*sin(v[1][4])-p.k*(SQRT_3/2.)*v[0][4]*sin(v[1][4]);
}
void _resolve_three_monovalent_at_four_coordinated_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)  
{
t_vec u, p, q;
u.i=r[*group].i-r[group[1]].i, u.j=r[*group].j-r[group[1]].j, u.k=r[*group].k-r[group[1]].k; multiple_vec_scalar(&u,&u,1./sqrt(calc_vec_norm(&u)));
_get_orthogonal_vector(&p,&u); multiple_vec_scalar(&p,&p,1./sqrt(calc_vec_norm(&p)));
vec_vec_vmult(&q,&u,&p); multiple_vec_scalar(&q,&q,1./sqrt(calc_vec_norm(&q)));
r[group[2]].i=r[*group].i+u.i*v[0][2]*cos(v[1][2])+q.i*v[0][2]*sin(v[1][2]),                                              r[group[2]].j=r[*group].j+u.j*v[0][2]*cos(v[1][2])+q.j*v[0][2]*sin(v[1][2]),                                             r[group[2]].k=r[*group].k+u.k*v[0][2]*cos(v[1][2])+q.k*v[0][2]*sin(v[1][2]);
r[group[3]].i=r[*group].i+u.i*v[0][3]*cos(v[1][3])-q.i*(1./2.)*v[0][3]*sin(v[1][3])+p.i*(SQRT_3/2.)*v[0][3]*sin(v[1][3]), r[group[3]].j=r[*group].j+u.j*v[0][3]*cos(v[1][3])-q.j*(1./2.)*v[0][3]*sin(v[1][3])+p.j*(SQRT_3/2.)*v[0][3]*sin(v[1][3]), r[group[3]].k=r[*group].k+u.k*v[0][3]*cos(v[1][3])-q.k*(1./2.)*v[0][3]*sin(v[1][3])+p.k*(SQRT_3/2.)*v[0][3]*sin(v[1][3]);
r[group[4]].i=r[*group].i+u.i*v[0][4]*cos(v[1][4])-q.i*(1./2.)*v[0][4]*sin(v[1][4])-p.i*(SQRT_3/2.)*v[0][4]*sin(v[1][4]), r[group[4]].j=r[*group].j+u.j*v[0][4]*cos(v[1][4])-q.j*(1./2.)*v[0][4]*sin(v[1][4])-p.j*(SQRT_3/2.)*v[0][4]*sin(v[1][4]), r[group[4]].k=r[*group].k+u.k*v[0][4]*cos(v[1][4])-q.k*(1./2.)*v[0][4]*sin(v[1][4])-p.k*(SQRT_3/2.)*v[0][4]*sin(v[1][4]);
}
void _resolve_two_monovalent_at_four_coordinated_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
double csA;
t_vec u, p, q;
u.i=r[*group].i-r[group[1]].i, u.j=r[*group].j-r[group[1]].j, u.k=r[*group].k-r[group[1]].k; multiple_vec_scalar(&u,&u,1./sqrt(calc_vec_norm(&u)));
q.i=r[group[2]].i-r[*group].i, q.j=r[group[2]].j-r[*group].j, q.k=r[group[2]].k-r[*group].k; csA=calc_vec_vec_scalar_product(&u,&q);
p.i=q.i-u.i*csA, p.j=q.j-u.j*csA, p.k=q.k-u.k*csA; multiple_vec_scalar(&p,&p,1./sqrt(calc_vec_norm(&p)));
vec_vec_vmult(&q,&u,&p); multiple_vec_scalar(&q,&q,1./sqrt(calc_vec_norm(&q)));
r[group[3]].i=r[*group].i+u.i*v[0][3]*cos(v[1][3])-q.i*(1./2.)*v[0][3]*sin(v[1][3])+p.i*(SQRT_3/2.)*v[0][3]*sin(v[1][3]), r[group[3]].j=r[*group].j+u.j*v[0][3]*cos(v[1][3])-q.j*(1./2.)*v[0][3]*sin(v[1][3])+p.j*(SQRT_3/2.)*v[0][3]*sin(v[1][3]), r[group[3]].k=r[*group].k+u.k*v[0][3]*cos(v[1][3])-q.k*(1./2.)*v[0][3]*sin(v[1][3])+p.k*(SQRT_3/2.)*v[0][3]*sin(v[1][3]);
r[group[4]].i=r[*group].i+u.i*v[0][4]*cos(v[1][4])-q.i*(1./2.)*v[0][4]*sin(v[1][4])-p.i*(SQRT_3/2.)*v[0][4]*sin(v[1][4]), r[group[4]].j=r[*group].j+u.j*v[0][4]*cos(v[1][4])-q.j*(1./2.)*v[0][4]*sin(v[1][4])-p.j*(SQRT_3/2.)*v[0][4]*sin(v[1][4]), r[group[4]].k=r[*group].k+u.k*v[0][4]*cos(v[1][4])-q.k*(1./2.)*v[0][4]*sin(v[1][4])-p.k*(SQRT_3/2.)*v[0][4]*sin(v[1][4]);
}
void _resolve_one_monovalent_at_four_coordinated_center(unsigned int n,unsigned int m,double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1],unsigned int group[MAX_NEIGHBORS+1],t_vec *r)
{
t_vec u;
u.i=r[group[1]].i+r[group[2]].i+r[group[3]].i-3.*r[*group].i, u.j=r[group[1]].j+r[group[2]].j+r[group[3]].j-3.*r[*group].j, u.k=r[group[1]].k+r[group[2]].k+r[group[3]].k-3.*r[*group].k; multiple_vec_scalar(&u,&u,1./sqrt(calc_vec_norm(&u)));
r[group[4]].i=r[*group].i-u.i*v[0][4], r[group[4]].j=r[*group].j-u.j*v[0][4], r[group[4]].k=r[*group].k-u.k*v[0][4];
}

//This function resolves coordinates of onevalent atoms (mostrly hydrogens)
//Note. It does NOT perform conformation search
unsigned int resolve_onevalent(t_clist *neighbors,t_vec *r,t_mol *mol)
{
register unsigned int _i, _j, _t;
unsigned int group[MAX_NEIGHBORS+1], xid, n, m, item0, item2, count=0;
double v[MAX_NEIGHBORS+1][MAX_NEIGHBORS+1];
 
_t=mol->natoms;
while (_t--)
  if ( ( (isnan(r[_t].i)))||( (isnan(r[_t].j)))||( (isnan(r[_t].k))) )
    {
    if ( (neighbors->list[_t].size!=1)||(neighbors->list[xid=*neighbors->list[_t].list].size==1)||(neighbors->list[xid].size>MAX_NEIGHBORS) )
      { LABEL_ERROR_NIMPLMENTED: ylib_errno=YERROR_NIMPLEMENTED; return (unsigned int)-1; }
    //Fill in group
    group[0]=xid, n=0, m=0;
    _i=neighbors->list[xid].size;
    while (_i--)
      if ( ( (isnan(r[neighbors->list[xid].list[_i]].i)))||( (isnan(r[neighbors->list[xid].list[_i]].j)))||( (isnan(r[neighbors->list[xid].list[_i]].k))) )
        {//The group: | xid | #m | #n |
        if (neighbors->list[neighbors->list[xid].list[_i]].size!=1) goto LABEL_ERROR_NIMPLMENTED;
        group[1+m+n++]=neighbors->list[xid].list[_i]; 
        }
      else 
        {
        _j=n; while (_j--) group[1+m+_j+1]=group[1+m+_j];
        group[1+m++]=neighbors->list[xid].list[_i];
        } 
    //Fill in bonds
    memset(v,0,sizeof(double)*(neighbors->list[xid].size+1)*(MAX_NEIGHBORS+1));
    _i=mol->size_b;
    while (_i--) 
           if (mol->ff_b[_i].atom[0]==xid)
             {
             if ( ((item0=find_in_row(mol->ff_b[_i].atom[1],n+m,&group[1]))==(unsigned int)-1)||( (v[++item0][0]))||( (v[0][item0])) ) goto LABEL_ERROR_IMPOSSIBLE;
             v[0][item0]=v[item0][0]=mol->ff_b[_i].v; 
             }
      else if (mol->ff_b[_i].atom[1]==xid)
             {
             if ( ((item2=find_in_row(mol->ff_b[_i].atom[0],n+m,&group[1]))==(unsigned int)-1)||( (v[++item2][0]))||( (v[0][item2])) ) goto LABEL_ERROR_IMPOSSIBLE;
             v[0][item2]=v[item2][0]=mol->ff_b[_i].v; 
             }
    _i=neighbors->list[xid].size+1; while (--_i) if ( (!(v[0][_i]))||(!(v[_i][0])) ) { LABEL_ERROR_IMPOSSIBLE: ylib_errno=YERROR_IMPOSSIBLE; return (unsigned int)-1; }
    //Fill in angles
    _i=mol->size_g;
    while (_i--) 
      if (mol->ff_g[_i].atom[1]==xid)
        {
        if ( ((item0=find_in_row(mol->ff_g[_i].atom[0],n+m,&group[1]))==(unsigned int)-1)||((item2=find_in_row(mol->ff_g[_i].atom[2],n+m,&group[1]))==(unsigned int)-1)||(v[++item2][++item0])||(v[item0][item2]) ) goto LABEL_ERROR_IMPOSSIBLE;
        v[item0][item2]=v[item2][item0]=mol->ff_g[_i].v;
        }
    _i=neighbors->list[xid].size+1; while (--_i) { _j=_i; while (--_j) if ( (!(v[_i][_j]))||(!(v[_j][_i])) ) goto LABEL_ERROR_IMPOSSIBLE; }
    //Resolve coords
    switch (neighbors->list[xid].size) // M+N
      {
      case  1 : {//M==0
        _resolve_one_monovalent_at_one_coordinated_center(n,m,v,group,r);
        break;  }
      case  2 : {//M==0 || M==1
        if (!m) _resolve_two_monovalent_at_two_coordinated_center(n,m,v,group,r);
        else    _resolve_one_monovalent_at_two_coordinated_center(n,m,v,group,r);
        break;  }
      case  3 : {//M==0 || M==1 || M==2
        //Check planarity
        _j=(unsigned int)-1, _i=mol->size_i;
        while (_i--)
          if ( (!mol->ff_i[_i].type)&&(mol->ff_i[_i].atom[0]==xid) )
            {
            if ( (_j!=(unsigned int)-1)||(find_in_row(mol->ff_i[_i].atom[1],n+m,&group[1])==(unsigned int)-1)||(find_in_row(mol->ff_i[_i].atom[2],n+m,&group[1])==(unsigned int)-1)||(find_in_row(mol->ff_i[_i].atom[3],n+m,&group[1])==(unsigned int)-1) ) goto LABEL_ERROR_IMPOSSIBLE;
            v[0][0]=fmod(mol->ff_i[_i].v,PI); if ( ( (v[0][0]<0.+SMALL)&&(v[0][0]>0-SMALL) )||(mol->ff_i[_i].v>+PI-SMALL)||(mol->ff_i[_i].v<-PI+SMALL) ) _j=0; else _j=1; 
            }
        if (!_j) 
          {
               if (!m)   _resolve_three_monovalent_at_three_coordinated_planar_center(n,m,v,group,r);
          else if (m==1) _resolve_two_monovalent_at_three_coordinated_planar_center(n,m,v,group,r);
          else           _resolve_one_monovalent_at_three_coordinated_planar_center(n,m,v,group,r);
          }
        else
          { 
               if (!m)   _resolve_three_monovalent_at_three_coordinated_nonplanar_center(n,m,v,group,r);
          else if (m==1) _resolve_two_monovalent_at_three_coordinated_nonplanar_center(n,m,v,group,r);
          else           _resolve_one_monovalent_at_three_coordinated_nonplanar_center(n,m,v,group,r);
          }
        break;  }
      case  4 : {//M==0 || M==1 || M==2 || M==3
             if (!m)   _resolve_four_monovalent_at_four_coordinated_center(n,m,v,group,r);
        else if (m==1) _resolve_three_monovalent_at_four_coordinated_center(n,m,v,group,r);  
        else if (m==2) _resolve_two_monovalent_at_four_coordinated_center(n,m,v,group,r);
        else           _resolve_one_monovalent_at_four_coordinated_center(n,m,v,group,r);
        break;  }
      default : { ylib_errno=YERROR_NIMPLEMENTED; return (unsigned int)-1; }
      }
    }
return count;
}




