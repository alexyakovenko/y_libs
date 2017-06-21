
#include "../y_libs/y_ff.h"


//------------------------------------    B O N D E D   P R I M I T I V E S   ------------------------------------------

//This function calculates bond energy
inline double _calc_benergy_yff1(double K,double V,t_vec *ri,t_vec *rj)
{
double _R;

_R=sqrt(calc_distance(ri,rj));
return K*sqrd(_R-V);
}
//This function calculates bond energy and its gradient
inline double _calc_bgrad_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *gi,t_vec *gj)
{
register double _d;
double R, dr[6];

calc_bond_derivative(&R,ri,rj,dr);
_d=2.*K*(R-V);
gi->i+=dr[0]*_d, gi->j+=dr[1]*_d, gi->k+=dr[2]*_d;
gj->i+=dr[3]*_d, gj->j+=dr[4]*_d, gj->k+=dr[5]*_d;
return K*sqrd(R-V);
}
//This function calculates angle energy
inline double _calc_genergy_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk)
{
register double _A;

_A=acos(calc_cos(ri,rj,rk));
return K*sqrd(_A-V);
}
//This function calculates angle energy and its gradient
inline double _calc_ggrad_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *gi,t_vec *gj,t_vec *gk)
{
register double _A, _d;
double csA, snA, dr[9];

calc_angle_derivative(&csA,&snA,ri,rj,rk,dr);
_A=acos(csA);
_d=2.*K*(_A-V);
gi->i+=dr[0]*_d, gi->j+=dr[1]*_d, gi->k+=dr[2]*_d;
gj->i+=dr[3]*_d, gj->j+=dr[4]*_d, gj->k+=dr[5]*_d;
gk->i+=dr[6]*_d, gk->j+=dr[7]*_d, gk->k+=dr[8]*_d;
return K*sqrd(_A-V);
}
//This function calculates improper energy
inline double _calc_ienergy_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl)
{
register double _F;
double csF, snF;

calc_dih_angle(&csF,&snF,ri,rj,rk,rl);
_F=calc_trig_angle(csF,snF);
return K*sqrd(_F-V);
}
//This function calculates improper energy and its gradient
inline double _calc_igrad_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl,t_vec *gi,t_vec *gj,t_vec *gk,t_vec *gl)
{
register double _d, _F;
double csF, snF, dr[0xC];

calc_dih_derivative(&csF,&snF,ri,rj,rk,rl,dr);
_F=calc_trig_angle(csF,snF);
_d=2.*K*(_F-V);
gi->i+=dr[0x0]*_d, gi->j+=dr[0x1]*_d, gi->k+=dr[0x2]*_d;
gj->i+=dr[0x3]*_d, gj->j+=dr[0x4]*_d, gj->k+=dr[0x5]*_d;
gk->i+=dr[0x6]*_d, gk->j+=dr[0x7]*_d, gk->k+=dr[0x8]*_d;
gl->i+=dr[0x9]*_d, gl->j+=dr[0xA]*_d, gl->k+=dr[0xB]*_d;
return K*sqrd(_F-V);
}
//This function calculates improper energy with respect to circular rotation
inline double _calc_ienergy_PI_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl)
{
register double _F;
double csF, snF;

calc_dih_angle(&csF,&snF,ri,rj,rk,rl);
_F=calc_trig_angle(csF,snF);
if (_F<0.) return K*sqrd(_F+fabs(V)); //approach -|V|
else       return K*sqrd(_F-fabs(V)); //approach +|V|
}
//This function calculates improper energy and its gradient
inline double _calc_igrad_PI_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl,t_vec *gi,t_vec *gj,t_vec *gk,t_vec *gl)
{
register double _d, _F, _V;
double csF, snF, dr[0xC];

calc_dih_derivative(&csF,&snF,ri,rj,rk,rl,dr);
_F=calc_trig_angle(csF,snF);
if (_F<0.) _V=-fabs(V); //approach -|V|
else       _V=+fabs(V); //approach +|V|
_d=2.*K*(_F-_V);
gi->i+=dr[0x0]*_d, gi->j+=dr[0x1]*_d, gi->k+=dr[0x2]*_d;
gj->i+=dr[0x3]*_d, gj->j+=dr[0x4]*_d, gj->k+=dr[0x5]*_d;
gk->i+=dr[0x6]*_d, gk->j+=dr[0x7]*_d, gk->k+=dr[0x8]*_d;
gl->i+=dr[0x9]*_d, gl->j+=dr[0xA]*_d, gl->k+=dr[0xB]*_d;
return K*sqrd(_F-_V);
}
//This function calculates torsion energy
inline double _calc_denergy_yff1(double K[SIZE_DIH],double V[SIZE_DIH],double N[SIZE_DIH],t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl)
{
register unsigned int _i;
register double _F, _e=0.;
double csF, snF;

calc_dih_angle(&csF,&snF,ri,rj,rk,rl);
_F=calc_trig_angle(csF,snF);
_i=SIZE_DIH; while (_i--) _e+=K[_i]*(1.+cos(N[_i]*_F-V[_i]));
return _e;
}
//This function calculates torsion energy and its gradient
inline double _calc_dgrad_yff1(double K[SIZE_DIH],double V[SIZE_DIH],double N[SIZE_DIH],t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl,t_vec *gi,t_vec *gj,t_vec *gk,t_vec *gl)
{
register unsigned int _i;
register double _d, _F, _e=0.;
double csF, snF, dr[0xC];

calc_dih_derivative(&csF,&snF,ri,rj,rk,rl,dr);
_F=calc_trig_angle(csF,snF);
_i=SIZE_DIH;
while (_i--)
  {
  _e+=K[_i]*(1.  +cos(N[_i]*_F-V[_i]));
  _d=-K[_i]*N[_i]*sin(N[_i]*_F-V[_i]);
  gi->i+=dr[0x0]*_d, gi->j+=dr[0x1]*_d, gi->k+=dr[0x2]*_d;
  gj->i+=dr[0x3]*_d, gj->j+=dr[0x4]*_d, gj->k+=dr[0x5]*_d;
  gk->i+=dr[0x6]*_d, gk->j+=dr[0x7]*_d, gk->k+=dr[0x8]*_d;
  gl->i+=dr[0x9]*_d, gl->j+=dr[0xA]*_d, gl->k+=dr[0xB]*_d;
  }
return _e;
}

//This function summs energy of bonded interations of all atoms of the molecule
inline double summ_bnmol_energy_yff1(t_mol *mol,register t_vec *r)
{
register unsigned int _i;
register union {
               t_ff_b *_b;
               t_ff_g *_g;
               t_ff_i *_i;
               t_ff_d *_d;
               }ff;
register double _e=0.;

//Stage 1. Summ dihs (the weakest term)
_i=mol->size_d, ff._d=mol->ff_d; while (_i--) { _e+=_calc_denergy_yff1(ff._d->k,ff._d->v,ff._d->n,&r[ff._d->atom[0]],&r[ff._d->atom[1]],&r[ff._d->atom[2]],&r[ff._d->atom[3]]), ff._d++; }
//Stage 2. Summ imprs
_i=mol->size_i, ff._i=mol->ff_i;
while (_i--) 
  {//Handle the circular rotations
  if ( ((ff._i->v>-PI-SMALL)&&(ff._i->v<-0.9*PI))||((ff._i->v>+0.9*PI)&&(ff._i->v<+PI+SMALL)) )
    _e+=_calc_ienergy_PI_yff1(ff._i->k,ff._i->v,&r[ff._i->atom[0]],&r[ff._i->atom[1]],&r[ff._i->atom[2]],&r[ff._i->atom[3]]);
  else
    _e+=_calc_ienergy_yff1(ff._i->k,ff._i->v,&r[ff._i->atom[0]],&r[ff._i->atom[1]],&r[ff._i->atom[2]],&r[ff._i->atom[3]]);
  ff._i++;
  }
//Stage 3. Summ constraints
_i=mol->size_c, ff._b=mol->ff_c; while (_i--) { _e+=_calc_benergy_yff1(ff._b->k,ff._b->v,&r[ff._b->atom[0]],&r[ff._b->atom[1]]), ff._b++; }
//Stage 4. Summ angles
_i=mol->size_g, ff._g=mol->ff_g; while (_i--) { _e+=_calc_genergy_yff1(ff._g->k,ff._g->v,&r[ff._g->atom[0]],&r[ff._g->atom[1]],&r[ff._g->atom[2]]), ff._g++; }
//Stage 5. Summ bonds (the hardest term)
_i=mol->size_b, ff._b=mol->ff_b; while (_i--) { _e+=_calc_benergy_yff1(ff._b->k,ff._b->v,&r[ff._b->atom[0]],&r[ff._b->atom[1]]), ff._b++; }
return _e;
}
//This function summs energies and gradients of bonded interations of all atoms of the molecule
inline double summ_bnmol_grad_yff1(register t_vec *g,t_mol *mol,register t_vec *r)
{
register unsigned int _i, _ai, _aj, _ak, _al;
register union {
               t_ff_b *_b;
               t_ff_g *_g;
               t_ff_i *_i;
               t_ff_d *_d;
               }ff;
register double _e=0.;

//Stage 1. Summ dihs (the weakest term)
_i=mol->size_d, ff._d=mol->ff_d; while (_i--) { _ai=ff._d->atom[0], _aj=ff._d->atom[1], _ak=ff._d->atom[2], _al=ff._d->atom[3], _e+=_calc_dgrad_yff1(ff._d->k,ff._d->v,ff._d->n,&r[_ai],&r[_aj],&r[_ak],&r[_al],&g[_ai],&g[_aj],&g[_ak],&g[_al]), ff._d++; }
//Stage 2. Summ imprs
_i=mol->size_i, ff._i=mol->ff_i; 
while (_i--)
  {//Handle the circular rotations
  _ai=ff._i->atom[0], _aj=ff._i->atom[1], _ak=ff._i->atom[2], _al=ff._i->atom[3];
  if ( ((ff._i->v>-PI-SMALL)&&(ff._i->v<-0.9*PI))||((ff._i->v>+0.9*PI)&&(ff._i->v<+PI+SMALL)) )
    _e+=_calc_igrad_PI_yff1(ff._i->k,ff._i->v,&r[_ai],&r[_aj],&r[_ak],&r[_al],&g[_ai],&g[_aj],&g[_ak],&g[_al]);
  else 
    _e+=_calc_igrad_yff1(ff._i->k,ff._i->v,&r[_ai],&r[_aj],&r[_ak],&r[_al],&g[_ai],&g[_aj],&g[_ak],&g[_al]);
  ff._i++;
  }
//Stage 3. Summ constraints
_i=mol->size_c, ff._b=mol->ff_c; while (_i--) { _ai=ff._b->atom[0], _aj=ff._b->atom[1]; _e+=_calc_bgrad_yff1(ff._b->k,ff._b->v,&r[_ai],&r[_aj],&g[_ai],&g[_aj]), ff._b++; }
//Stage 4. Summ angles
_i=mol->size_g, ff._g=mol->ff_g; while (_i--) { _ai=ff._g->atom[0], _aj=ff._g->atom[1], _ak=ff._g->atom[2]; _e+=_calc_ggrad_yff1(ff._g->k,ff._g->v,&r[_ai],&r[_aj],&r[_ak],&g[_ai],&g[_aj],&g[_ak]), ff._g++; }
//Stage 5. Summ bonds (the hardest term)
_i=mol->size_b, ff._b=mol->ff_b; while (_i--) { _ai=ff._b->atom[0], _aj=ff._b->atom[1]; _e+=_calc_bgrad_yff1(ff._b->k,ff._b->v,&r[_ai],&r[_aj],&g[_ai],&g[_aj]), ff._b++; }
return _e;
}


//------------------------------------    N O N B O N D E D   P R I M I T I V E S   ------------------------------------------


//This function calculates energy of classic nonbonded interatons in YFF1. 
inline double _calc_atom__nbenergy_yff1(register double Aij,register double Bij,double qi,double qj,t_vec *ri,t_vec *rj)
{
register double _r, _rr, _rrrrrr;
if ((_rr=sqrd(ri->i-rj->i)+sqrd(ri->j-rj->j)+sqrd(ri->k-rj->k))>SMALL2*SMALL2) { _r=sqrt(_rr), _rrrrrr=_rr*_rr*_rr; }
else { _r=SMALL2, _rr=SMALL2*SMALL2, _rrrrrr=SMALL2*SMALL2*SMALL2*SMALL2*SMALL2*SMALL2; } //Singularity!!!
return (Aij/_rrrrrr-Bij)/_rrrrrr+COULOMB_K*qi*qj/_r;
}
//This function calculates energy and force of classic nonbonded interatons in YFF1. 
inline double _calc_atom__nbgrad_yff1(register double Aij,register double Bij,double qi,double qj,t_vec *ri,t_vec *rj,t_vec *gi,t_vec *gj)
{
register double _r, _rr, _rrrrrr, Qij, _d, _vi, _vj, _vk;
_vi=ri->i-rj->i, _vj=ri->j-rj->j, _vk=ri->k-rj->k;
if ((_rr=_vi*_vi+_vj*_vj+_vk*_vk)>SMALL2*SMALL2) { _r=sqrt(_rr), _rrrrrr=_rr*_rr*_rr; }
else { _r=SMALL2, _rr=SMALL2*SMALL2, _rrrrrr=SMALL2*SMALL2*SMALL2*SMALL2*SMALL2*SMALL2; } //Singularity!!!
Qij=COULOMB_K*qi*qj, _d=((-12.*Aij/_rrrrrr+6.*Bij)/_rrrrrr-Qij/_r)/_rr;
_vi*=_d, _vj*=_d, _vk*=_d, gi->i+=_vi, gi->j+=_vj, gi->k+=_vk, gj->i-=_vi, gj->j-=_vj, gj->k-=_vk;
return (Aij/_rrrrrr-Bij)/_rrrrrr+Qij/_r;
}

//This function summs nonbonded interations energies of all anchors of the same molecule
inline double summ_nbmol_energy_yff1(t_mol *mol,t_vec *r,double **A,double **B)
{
register unsigned int _i, _j, _k, *tj, *ep;
register double _e=0., *qi, *qj, *Ai, *Bi;
register t_vec *ri, *rj;

_i=mol->natoms; 
while (_i--)
  {
  Ai=A[mol->fftypes[_i]], Bi=B[mol->fftypes[_i]], _j=_i, tj=&mol->fftypes[_i], qj=qi=&mol->charges[_i], rj=ri=&r[_i];
  //Duff's device on amount of exclusions to cover up to 4x4 cases
  _k=mol->excl[_i].size, ep=&mol->excl[_i].list[_k];
  switch (_k)
    {
    case 16 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case 15 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case 14 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case 13 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case 12 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case 11 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case 10 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  9 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  8 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  7 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  6 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  5 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  4 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  3 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  2 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  1 : { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; }
    case  0 : { ZERO_EXCL: while (_j--) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } break; }
    default : { do { --ep; while (--_j!=*ep) { --tj, --qj, --rj, _e+=_calc_atom__nbenergy_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj); } tj--, qj--, rj--; } while (--_k); goto ZERO_EXCL; } //Unexpected general case (mol->excl[_i].size>16)
    }
  }
return _e;
}
//This function summs nonbonded interations energies and gradients of all anchors of the same molecule
inline double summ_nbmol_grad_yff1(t_vec *g, t_mol *mol, t_vec *r, double **A, double**B)
{
register unsigned int _i, _j, _k, *tj, *ep;
register double _e=0., *qi, *qj, *Ai, *Bi;
register t_vec *ri, *rj, *gi, *gj;

_i=mol->natoms; 
while (_i--)
  {
  Ai=A[mol->fftypes[_i]], Bi=B[mol->fftypes[_i]], _j=_i, tj=&mol->fftypes[_i], qj=qi=&mol->charges[_i], gj=gi=&g[_i], rj=ri=&r[_i];
  //Duff's device on amount of exclusions to cover ep to 4x4 cases
  _k=mol->excl[_i].size, ep=&mol->excl[_i].list[_k];
  switch (_k)
    {
    case 16 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case 15 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case 14 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case 13 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case 12 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case 11 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case 10 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  9 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  8 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  7 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  6 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  5 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  4 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  3 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  2 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  1 : { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; }
    case  0 : { ZERO_EXCL: while (_j--) { 
--tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); 
} break; }
    default : { do { --ep; while (--_j!=*ep) { --tj, --qj, --gj, --rj, _e+=_calc_atom__nbgrad_yff1(Ai[*tj],Bi[*tj],*qi,*qj,ri,rj,gi,gj); } tj--, qj--, gj--, rj--; } while (--_k); goto ZERO_EXCL; } //Unexpected general case (mol->excl[_i].size>16)
    }
  }
return _e;
}


//------------------------------------    U T I L I T I E S   ------------------------------------------


//This function calculates energy and its gradient for a stand alone mol
double calc_mol_energy_yff1(t_mol *mol,t_vec *rvecs,double **A,double **B)
{
register double _e=0.;
//Calculate energies
_e+=summ_nbmol_energy_yff1(mol,rvecs,A,B);
_e+=summ_bnmol_energy_yff1(mol,rvecs);
return _e;
}
//This function calculates energy and its gradient for a stand alone mol
char calc_mol_grad_yff1(double *e,unsigned int n,double *x,double *g,double **G,va_list stack)
{
t_vec *rvecs, *gvecs;
t_mol *mol;
double **A, **B;

va_list _stack;

//Stage 0. Unwrap stack
va_copy(_stack,stack);
mol=va_arg(_stack,t_mol*);
A=va_arg(_stack,double**);
B=va_arg(_stack,double**);

//Stage I. Calculate gradients
rvecs=(t_vec*)x, gvecs=(t_vec*)g; 
(*e)=0., memset(gvecs,FALSE,sizeof(t_vec)*mol->natoms);
(*e)+=summ_nbmol_grad_yff1(gvecs,mol,rvecs,A,B);
(*e)+=summ_bnmol_grad_yff1(gvecs,mol,rvecs);

//End with the stack and exit
va_end(_stack);
return TRUE;
}
//This function calculates energy and its gradient for a stand alone mol
char calc_mol_grad_yff1_restraints(double *e,unsigned int n,double *x,double *g,double **G,va_list stack)
{
register unsigned int _i, _j;
t_vec *rvecs, *gvecs;
t_mol *mol;
double **A, **B, *restraints_k, *restraints_x;
t_list *restraints;

va_list _stack;

//Stage 0. Unwrap stack
va_copy(_stack,stack);
mol=va_arg(_stack,t_mol*);
A=va_arg(_stack,double**);
B=va_arg(_stack,double**);
restraints=va_arg(_stack,t_list*);
restraints_k=va_arg(_stack,double*);
restraints_x=va_arg(_stack,double*);

//Stage I. Calculate gradients
rvecs=(t_vec*)x, gvecs=(t_vec*)g; 
(*e)=0., memset(gvecs,FALSE,sizeof(t_vec)*mol->natoms);
(*e)+=summ_nbmol_grad_yff1(gvecs,mol,rvecs,A,B);
(*e)+=summ_bnmol_grad_yff1(gvecs,mol,rvecs);
//Stage II. Add restraints
_i=restraints->size;
while (_i--)
  {
  _j=restraints->list[_i];
  (*e) +=restraints_k[_i]*(x[_j]-restraints_x[_i])*(x[_j]-restraints_x[_i]);
  g[_j]+=restraints_k[_i]*(x[_j]-restraints_x[_i])*2.;
  }
//End with the stack and exit
va_end(_stack);
return TRUE;
}

//This function perform energy minimization of stand alone molecule
char optimize_mol(unsigned int nsteps,double tol,t_vec *rvecs,char *label,t_mol *mol,double **A,double **B)  
{
double e, *x[3], *g[2], *p;
//Allocate memory
if (!(p=(double*)malloc(sizeof(t_vec)*mol->natoms*(1+3+2)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { x[0]=p+mol->natoms*3, x[1]=x[0]+mol->natoms*3, x[2]=x[1]+mol->natoms*3, g[0]=x[2]+mol->natoms*3, g[1]=g[0]+mol->natoms*3; memcpy(x[0],rvecs,sizeof(t_vec)*mol->natoms); }
//Run the minimization
//if (!(polak_ribiere(&e,nsteps,tol,.1*SMALL,.1,mol->natoms*3,x,g,0x0,p,calc_mol_grad_yff1,line_search_square_besier,label,mol,A,B)))
if (!(polak_ribiere_points(&e,nsteps,tol,SMALL2,.1,mol->natoms*3,x,g,0x0,p,calc_mol_grad_yff1,line_search_square_fapproximation,label,mol,A,B)))
  {
  if ( (ylib_errno==YERROR_NCONVERGED)||(ylib_errno==YERROR_SUSPICIOUS) )
    { 
    memcpy(rvecs,x[0],sizeof(t_vec)*mol->natoms); 
    free(p);
    return NTNF;
    }
  else return FALSE;
  }
//Save results, free memory and exit routine
memcpy(rvecs,x[0],sizeof(t_vec)*mol->natoms);
free(p);
return TRUE;
}
//The same as previous but can restraint subset of atomic coordinates
char optimize_mol_restraints(unsigned int nsteps,double tol,t_vec *rvecs,char *label,t_mol *mol,double **A,double **B,t_list *restraints,double *restraints_k)
{
register unsigned int _i;
double e, *x[3], *g[2], *p, *restraints_x;
//Allocate memory
if (!(p=(double*)malloc(sizeof(t_vec)*mol->natoms*(1+3+2)+sizeof(double)*restraints->size))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else 
  {
  x[0]=p+mol->natoms*3, x[1]=x[0]+mol->natoms*3, x[2]=x[1]+mol->natoms*3, g[0]=x[2]+mol->natoms*3, g[1]=g[0]+mol->natoms*3, restraints_x=g[1]+mol->natoms*3;
  memcpy(x[0],rvecs,sizeof(t_vec)*mol->natoms); _i=restraints->size; while (_i--) restraints_x[_i]=x[0][restraints->list[_i]];
  }
//Run the minimization
//if (!(polak_ribiere(&e,nsteps,tol,.1*SMALL,.1,mol->natoms*3,x,g,0x0,p,calc_mol_grad_yff1_restraints,line_search_square_besier,label,mol,A,B,restraints,restraints_k,restraints_x)))
if (!(polak_ribiere_points(&e,nsteps,tol,SMALL2,.1,mol->natoms*3,x,g,0x0,p,calc_mol_grad_yff1_restraints,line_search_square_fapproximation,label,mol,A,B,restraints,restraints_k,restraints_x)))
  {
  if ( (ylib_errno==YERROR_NCONVERGED)||(ylib_errno==YERROR_SUSPICIOUS) )
    { 
    memcpy(rvecs,x[0],sizeof(t_vec)*mol->natoms); 
    free(p);
    return NTNF;
    }
  else return FALSE;
  }
//Save results, free memory and exit routine
memcpy(rvecs,x[0],sizeof(t_vec)*mol->natoms);
free(p);
return TRUE;
}


//This function perturb cys peptide bond
void _perturb_cis_peptide_bond(t_vec *r1,t_vec *o,t_vec *c,t_vec *n,t_vec *h,t_vec *r2)
{
t_tensor R;
t_vec u, t;
double rr;
//Stage I. Rotate
u.i=c->i-n->i, u.j=c->j-n->j, u.k=c->k-n->k;
if ((rr=calc_vec_norm(&u))>SMALL2)
  {
  rr=1./sqrt(rr); u.i*=rr, u.j*=rr, u.k*=rr;
//  rotate_around_uvector(&R,&u,0.,+1); //rotate +PI/2.
  rotate_around_uvector(&R,&u,-0.500,+0.866); //rotate +PI/2+PI/6.
  o->i-=c->i, o->j-=c->j, o->k-=c->k;
  multiple_origin_tensor_origin_vec(&t,&R,o);
  o->i=t.i+c->i, o->j=t.j+c->j, o->k=t.k+c->k;
//  rotate_around_uvector(&R,&u,0.,-1.); //rotate -PI/2.
  rotate_around_uvector(&R,&u,-0.500,-0.866); //rotate -PI/2-PI/6.
  h->i-=n->i, h->j-=n->j, h->k-=n->k;
  multiple_origin_tensor_origin_vec(&t,&R,h);
  h->i=t.i+n->i, h->j=t.j+n->j, h->k=t.k+n->k;
  }
//Stage II. Shift
//t.i=((r1->i+r2->i)-(c->i+n->i))/2., t.j=((r1->j+r2->j)-(c->j+n->j))/2., t.k=((r1->k+r2->k)-(c->k+n->k))/2.;
//o->i+=t.i, o->j+=t.j, o->k+=t.k;
//c->i+=t.i, c->j+=t.j, c->k+=t.k;
//n->i+=t.i, n->j+=t.j, n->k+=t.k;
//h->i+=t.i, h->j+=t.j, h->k+=t.k;
}
//This function check molecular structure for evident geometrical errors. It return amount of detected errors or (unsigned int)-1 on failure.
//So far it checks:
// 1. If any amide bonds in molecule A is in cys conformation. If so it swaps H and X atoms around N~CR=O atom.
// 2. If any <=6 cycle in molecule A is pierced by a bond of molecule B (A can be equal to B). If so the bond forming atoms are translated to the closes edge +1 A away.
unsigned int check_awful_geometry_errors(unsigned int nsteps,char (*ordera)[4],t_clist *neighborsa,t_vec *ra,t_mol *mola,char (*orderb)[4],t_clist *neighborsb,t_vec *rb,t_mol *molb,double **A,double **B) 
{
register unsigned int _i, _j, _k, _l;
t_list *restraints;
double *restraints_k;
unsigned int count=0;
void *vp;

//Check 1. If any amide bonds in molecule A is in cys conformation. If so it swaps H and X atoms around N~CR=O atom.
if (!(restraints=alloc_list(0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return (unsigned int)-1; } else restraints->size=0;
_i=mola->natoms;
while (_i--)
  if ( (mola->a[_i]==CHEM_ATOM_TYPE_CARBON)&&(neighborsa->list[_i].size==3)&&(!(ordera[_i][3]))&&(ordera[_i][2]==1)&&(neighborsa->list[neighborsa->list[_i].list[0]].size==1)&&
     ( (mola->a[neighborsa->list[_i].list[0]]==CHEM_ATOM_TYPE_OXYGEN)||(mola->a[neighborsa->list[_i].list[0]]==CHEM_ATOM_TYPE_SULFUR) ) )
    {   //Check neighbor #1
    _k=mola->cycles->size; while (_k--) if ( (mola->cycles->list[_k].size<8)&&(find_in_list(_i,&mola->cycles->list[_k])!=(unsigned int)-1) ) goto NEXT_I;
    _k=2, _l=1;
    TRY_NEXT_NEIGHBOUR: _j=neighborsa->list[_i].list[_k];
    if ( (mola->a[_j]==CHEM_ATOM_TYPE_NITROGEN)&&(neighborsa->list[_j].size==3)&&(!(ordera[_j][3]))&&(!(ordera[_j][2])) )
      {
           if ( (mola->a[neighborsa->list[_j].list[0]]!=CHEM_ATOM_TYPE_HYDROGEN)&&(mola->a[neighborsa->list[_j].list[1]]!=CHEM_ATOM_TYPE_HYDROGEN)&&(mola->a[neighborsa->list[_j].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) )
             {//H==[2]
             if (fabs(calc_dih_angle_value(&ra[neighborsa->list[_i].list[0]],&ra[_i],&ra[_j],&ra[neighborsa->list[_j].list[2]]))<PI/2.)
               {//Swap
               if (neighborsa->list[_j].list[0]==_i) 
                 _perturb_cis_peptide_bond(&ra[neighborsa->list[_i].list[_l]],&ra[neighborsa->list[_i].list[0]],&ra[_i],&ra[_j],&ra[neighborsa->list[_j].list[2]],&ra[neighborsa->list[_j].list[1]]);
               else
                 _perturb_cis_peptide_bond(&ra[neighborsa->list[_i].list[_l]],&ra[neighborsa->list[_i].list[0]],&ra[_i],&ra[_j],&ra[neighborsa->list[_j].list[2]],&ra[neighborsa->list[_j].list[0]]);
               if ((restraints->size%0xFF)<0xC)
                 { if (!(vp=realloc_list(restraints,restraints->size+(restraints->size%0xFF)+0xFF))) { LABEL_MEMORY_ERROR: free(restraints), restraints=0x0; goto LABEL_MEMORY_ERROR_0; } else restraints=(t_list*)vp; }
               restraints->list[restraints->size++]=neighborsa->list[_i].list[0]*3+0, restraints->list[restraints->size++]=neighborsa->list[_i].list[0]*3+1; restraints->list[restraints->size++]=neighborsa->list[_i].list[0]*3+2; //O coords
               restraints->list[restraints->size++]=neighborsa->list[_j].list[2]*3+0, restraints->list[restraints->size++]=neighborsa->list[_j].list[2]*3+1, restraints->list[restraints->size++]=neighborsa->list[_j].list[2]*3+2; //H coords
               restraints->list[restraints->size++]=_i*3+0, restraints->list[restraints->size++]=_i*3+1, restraints->list[restraints->size++]=_i*3+2; //C coords
               restraints->list[restraints->size++]=_j*3+0, restraints->list[restraints->size++]=_j*3+1, restraints->list[restraints->size++]=_j*3+2; //N coords
               count++; 
               }
             }
      else if ( (mola->a[neighborsa->list[_j].list[0]]!=CHEM_ATOM_TYPE_HYDROGEN)&&(mola->a[neighborsa->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mola->a[neighborsa->list[_j].list[2]]!=CHEM_ATOM_TYPE_HYDROGEN) )
             {//H==[1]
             if (fabs(calc_dih_angle_value(&ra[neighborsa->list[_i].list[0]],&ra[_i],&ra[_j],&ra[neighborsa->list[_j].list[1]]))<PI/2.)
               {//Swap
               if (neighborsa->list[_j].list[0]==_i) 
                 _perturb_cis_peptide_bond(&ra[neighborsa->list[_i].list[_l]],&ra[neighborsa->list[_i].list[0]],&ra[_i],&ra[_j],&ra[neighborsa->list[_j].list[1]],&ra[neighborsa->list[_j].list[2]]);
               else
                 _perturb_cis_peptide_bond(&ra[neighborsa->list[_i].list[_l]],&ra[neighborsa->list[_i].list[0]],&ra[_i],&ra[_j],&ra[neighborsa->list[_j].list[1]],&ra[neighborsa->list[_j].list[0]]);
               if ((restraints->size%0xFF)<0xC)
                 { if (!(vp=realloc_list(restraints,restraints->size+(restraints->size%0xFF)+0xFF))) goto LABEL_MEMORY_ERROR; else restraints=(t_list*)vp; }
               restraints->list[restraints->size++]=neighborsa->list[_i].list[0]*3+0, restraints->list[restraints->size++]=neighborsa->list[_i].list[0]*3+1; restraints->list[restraints->size++]=neighborsa->list[_i].list[0]*3+2; //O coords
               restraints->list[restraints->size++]=neighborsa->list[_j].list[1]*3+0, restraints->list[restraints->size++]=neighborsa->list[_j].list[1]*3+1, restraints->list[restraints->size++]=neighborsa->list[_j].list[1]*3+2; //H coords
               restraints->list[restraints->size++]=_i*3+0, restraints->list[restraints->size++]=_i*3+1, restraints->list[restraints->size++]=_i*3+2; //C coords
               restraints->list[restraints->size++]=_j*3+0, restraints->list[restraints->size++]=_j*3+1, restraints->list[restraints->size++]=_j*3+2; //N coords
               count++;
               } 
             }
      else if ( (mola->a[neighborsa->list[_j].list[0]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mola->a[neighborsa->list[_j].list[1]]!=CHEM_ATOM_TYPE_HYDROGEN)&&(mola->a[neighborsa->list[_j].list[2]]!=CHEM_ATOM_TYPE_HYDROGEN) )
             {//H==[0]
             if (fabs(calc_dih_angle_value(&ra[neighborsa->list[_i].list[0]],&ra[_i],&ra[_j],&ra[neighborsa->list[_j].list[0]]))<PI/2.)
               {//Swap
               if (neighborsa->list[_j].list[1]==_i) 
                 _perturb_cis_peptide_bond(&ra[neighborsa->list[_i].list[_l]],&ra[neighborsa->list[_i].list[0]],&ra[_i],&ra[_j],&ra[neighborsa->list[_j].list[0]],&ra[neighborsa->list[_j].list[2]]);
               else
                 _perturb_cis_peptide_bond(&ra[neighborsa->list[_i].list[_l]],&ra[neighborsa->list[_i].list[0]],&ra[_i],&ra[_j],&ra[neighborsa->list[_j].list[0]],&ra[neighborsa->list[_j].list[1]]);
               if ((restraints->size%0xFF)<0xC)
                 { if (!(vp=realloc_list(restraints,restraints->size+(restraints->size%0xFF)+0xFF))) goto LABEL_MEMORY_ERROR; else restraints=(t_list*)vp; }
               restraints->list[restraints->size++]=neighborsa->list[_i].list[0]*3+0, restraints->list[restraints->size++]=neighborsa->list[_i].list[0]*3+1; restraints->list[restraints->size++]=neighborsa->list[_i].list[0]*3+2; //O coords
               restraints->list[restraints->size++]=neighborsa->list[_j].list[0]*3+0, restraints->list[restraints->size++]=neighborsa->list[_j].list[0]*3+1, restraints->list[restraints->size++]=neighborsa->list[_j].list[0]*3+2; //H coords
               restraints->list[restraints->size++]=_i*3+0, restraints->list[restraints->size++]=_i*3+1, restraints->list[restraints->size++]=_i*3+2; //C coords
               restraints->list[restraints->size++]=_j*3+0, restraints->list[restraints->size++]=_j*3+1, restraints->list[restraints->size++]=_j*3+2; //N coords
               count++;
               }
             }
      }
    if ( (--_k)) { _l++; goto TRY_NEXT_NEIGHBOUR; }
    NEXT_I: ;
    }
if ( (count))
  {
  if (!(restraints_k=(double*)malloc(sizeof(double)*restraints->size))) goto LABEL_MEMORY_ERROR;
  else { _i=restraints->size; while (_i--) restraints_k[_i]=AMIDE_BONDS_RESTRAINTS_K; }
  //Preoptimize restrained coordinates
  if (!(optimize_mol_restraints(nsteps,SMALL,ra,"restrained optimization",mola,A,B,restraints,restraints_k)))
    {
    free(restraints), restraints=0x0; free(restraints_k), restraints_k=0x0;
    return (unsigned int)-1;
    }
  else { free(restraints), restraints=0x0; free(restraints_k); restraints_k=0x0; }
  }
else
  { free(restraints), restraints=0x0; }
//Exit
return count;
}


