#define Y_FFSYS 0x1

#ifndef Y_MOL
#include "y_mol.h"
#endif
#ifndef Y_CHARGE
#include "y_charge.h"
#endif
#ifndef Y_INTERPOLATE
#include "y_interpolate.h"
#endif
#ifndef Y_FF
#include "y_ff.h"
#endif

enum { FFTYPE_NONE, FFTYPE_YFF1, FFTYPE_YFF2 };          //Enumerate Force Fields

#define N_ESLICES 0x9   //Amount of slices for accuracy of summation, should be odd!

//Active layer description
typedef struct{
              t_list *active_a;             //Active anchors list
              t_list *active_b;             //Active bonds list
              t_list *active_g;             //Active angles list
              t_list *active_i;             //Active impropers list
              t_list *active_d;             //Active dihedrals list
              t_list *active_p;             //Active pairs list             | <-  size_p -> |
              unsigned int size_p;          //Amount of Active-Active pairs | Active-Active | Active-Passive |
              }t_activem;                   //Active molecules components   | <-      active_p->size      -> |

//NOTE. linear grids are not suiteble for gradient calculations
typedef struct{
              unsigned int type;            //Grid types: FFTYPE_NONE, FFTYPE_YFF1, FFTYPE_YFF2
              t_vec ori;                    //Origin of grid {x,y,z} (left frontal bottom corner) |/_ 
              t_len len;                    //Lenght of grid {x,y,z} 
              double sp;                    //Grid spacering
              double ***A;                  //                                                        ,- repulive exponent grid
              double ***B;                  // YFF1 tricubic hermite spline w numerical derivatives: -|- attractive exponent grid
              double ***Q;                  //                                                        `- charge grid
              }t_ffdgrid;
typedef struct{
              unsigned int type;            //Grid types: FFTYPE_NONE, FFTYPE_YFF1, FFTYPE_YFF2
              t_vec ori;                    //Origin of grid {x,y,z} (left frontal bottom corner) |/_ 
              t_len len;                    //Lenght of grid {x,y,z} 
              double sp;                    //Grid spacering
              double (***A)[64];            //                                                        ,- repulive exponent grid
              double (***B)[64];            // YFF1 tricubic hermite spline w numerical derivatives: -|- attractive exponent grid
              double (***Q)[64];            //                                                        `- charge grid
              }t_fftcgrid;

//General ffsys structure
typedef struct{
              unsigned int     nmols;       //Number of molecules in mechanical system
              t_mol           **mols;       //Pointers to mols topologies
              unsigned int    natoms;       //Number of atoms in mechanical system
              unsigned int   naatoms;       //Number of active atoms in ffsys
              t_activem     *activem;       //Global massive of molecular active elements 
              unsigned int       *nr;       //Global massive of molecular starting atoms ids
              t_vec               *r;       //Global massive of atoms radius-vectors
              t_vec  *g, (*_g)[N_ESLICES];   //Global massive of atoms forces and its [N_ESLICES] submassive: g==_g[0]
              double              *q;       //Global massive of atoms charges
              char                *a;       //Global massive of atoms types
              t_rtree       **rtrees;       //Rotational trees of molecules
              }t_ffsys;                     //Molecular mechanical force field structure



#define SMALL_DISTORTION SMALL

//YFF0
#define YFF0_COULOMB_K 0.06   //coulomb scaling factor 
#define YFF0_Rc 14.0          //coulomb cutoff
#define YFF0_Rb  5.6          //vdW zone 
#define YFF0_sD  0.1          //smooth factor
#define YFF0_Rc2 (YFF0_Rc*YFF0_Rc)
#define YFF0_Rd  (YFF0_Rc-YFF0_Rb)
#define YFF0_C3  ((3.-YFF0_Rd/YFF0_Rb)/YFF0_Rb)
#define YFF0_C4  ((2.-YFF0_Rd/YFF0_Rb)/YFF0_Rb)

//YFF1
#define YFF1_14_SCALE 0.25
#define YFF1_Rc 14.0
#define YFF1_Rb 14.0
#define YFF1_Rb2 (YFF1_Rb*YFF1_Rb)
#define YFF1_Rc2 (YFF1_Rc*YFF1_Rc)
//Shifting factors of YFF1 - interpolation with a polinome ax^4+bx^3+c*x^2+d*x+e to maintain V(Rb), dV(Rb), d2V(Rb), V(Rc)=0., dV(Rc)=0.;
#define YFF1_Rd  (YFF1_Rc-YFF1_Rb)
#define YFF1_Rd2 (YFF1_Rd*YFF1_Rd)
#define YFF1_AA (-(15.*YFF1_Rc-12.*YFF1_Rb)/(YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc*YFF1_Rd2)) 
#define YFF1_AB (+(14.*YFF1_Rc-12.*YFF1_Rb)/(YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc*YFF1_Rd2*YFF1_Rd)) 
#define YFF1_AF(r,rr,rrrrrr) (-12./rrrrrr/rrrrrr/r+YFF1_AA*(r-YFF1_Rb)*2.         +YFF1_AB*(r-YFF1_Rb)*(r-YFF1_Rb)*3.         )
#define YFF1_AV(r,rr,rrrrrr) (+ 1./rrrrrr/rrrrrr  +YFF1_AA*(r-YFF1_Rb)*(r-YFF1_Rb)+YFF1_AB*(r-YFF1_Rb)*(r-YFF1_Rb)*(r-YFF1_Rb))
#define YFF1_BA (-(9.*YFF1_Rc-6.*YFF1_Rb)/(YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc*YFF1_Rd2)) 
#define YFF1_BB (+(8.*YFF1_Rc-6.*YFF1_Rb)/(YFF1_Rc2*YFF1_Rc2*YFF1_Rc2*YFF1_Rc*YFF1_Rd2*YFF1_Rd)) 
#define YFF1_BF(r,rr,rrrrrr) (-6./rrrrrr/r+YFF1_BA*(r-YFF1_Rb)*2.         +YFF1_BB*(r-YFF1_Rb)*(r-YFF1_Rb)*3.         )
#define YFF1_BV(r,rr,rrrrrr) (+1./rrrrrr  +YFF1_BA*(r-YFF1_Rb)*(r-YFF1_Rb)+YFF1_BB*(r-YFF1_Rb)*(r-YFF1_Rb)*(r-YFF1_Rb))
#define YFF1_QA (-(4.*YFF1_Rc-YFF1_Rb)/(YFF1_Rc2*YFF1_Rd2)) 
#define YFF1_QB (+(3.*YFF1_Rc-YFF1_Rb)/(YFF1_Rc2*YFF1_Rd2*YFF1_Rd)) 
#define YFF1_QF(r,rr,rrrrrr) (-1./rr+YFF1_QA*(r-YFF1_Rb)*2.         +YFF1_QB*(r-YFF1_Rb)*(r-YFF1_Rb)*3.         )
#define YFF1_QV(r,rr,rrrrrr) (+1./r +YFF1_QA*(r-YFF1_Rb)*(r-YFF1_Rb)+YFF1_QB*(r-YFF1_Rb)*(r-YFF1_Rb)*(r-YFF1_Rb))


//Usage of shifting factors
//_x=YFF1_Rx(_r), E=YFF1_pE(_A,_B,_Q), D=YFF1_pD(_A,_B,_Q), C=YFF1_pC(_A,_B,_Q), B=YFF1_pB(C,D,E), A=YFF1_pA(C,D,E), _e=YFF1_pV(_x,A,B,C,D,E), _f=YFF1_pF(_x,A,B,C,D)


//This function compile activem block
//Note. It DO NOT copy active_a list, just a pointer to it.
char create_activem(t_activem *activem,t_list *active_a,t_mol *mol);

//This function removes global root elements from activem
void remove_global_roots_activem(char bonded,char nonbonded,t_activem *activem,t_rtree *rtree,t_mol *mol);

//This function deletes activem lists
void free_activem(t_activem *activem);

//This function reads activem from hdd
t_activem *read_activem(FILE *in,unsigned int *nmols);

//This function writes activem to hdd
char write_activem(FILE *out,unsigned int nmols,t_activem *activem);

//     F F S Y S       C O M P I L A T O R     P A R T

//This function empty ffsys data structure
void free_ffsys(t_ffsys *ffsys);
//This function reads ffsys from hdd
t_ffsys *read_ffsys(FILE *in);
//This function writes ffsys to hdd
char write_ffsys(FILE *out,t_ffsys *ffsys);

//This function compile individual mol into ffsys after allocing memory for it
//It DO NOT copy rvecs
char ffsys_add_mol(unsigned int nmols,t_mol *mol,t_activem *activem,t_rtree *rtree,t_ffsys *ffsys);

//This function add water molucules template to ffsys (atom's coords are uninitialized)
char ffsys_add_sol(unsigned int nsols,t_top *top,t_ffsys *ffsys);

//-----------------------------   Y F F 0    P A R T   ---------------------------------------------

// Inspirated by 
// Gehlhaar D, Verkhivker G., Rejto P., Sherman C., Fogel D., Fogel L., Freer S., Chem Biol 1995, 2, 317.
// and 
// Fuhrmaann J., Rurainski A., Lenhof H. and Neumann D., J.Comput.Chem., 2009, 30, 1371-1378.

//Setup yff0.

//   ^ vdW                                    ^ Q
// F +.                                       |
//   ||                                       |       sQ              tQ     YFF0_Rc
//   | .                                      |       |                |       |
//   | |                                    F |       v                |       |
//   |  .                                  +- +-------.                |       |
//   |  |                                   2 |       |`-._            v       | 
//   |   .                                    |       |    ```-----....        |
//   |  A|   B       C    D                   |       |           C  D ``--..  v
// --+----.--.---------.-.._______---->     --+----.--+------------.---.-----``____---->        
// 0 |    ^\            / ^                  0|   A  B|                 _..--``       
// E +    | `__________`  |                   |       |      ____...--``
//   |    | ^          ^  |                 F |       |.-````   
//   |    | |          |  |               - - +-------`
//   |    A  B         C  D  } +/-Delta {   2 |  
//
//This function setups YFF0
inline void setup_yff0(t_top *top);

//------------------------------------    N O N B O N D E D   P R I M I T I V E S   ------------------------------------------

//This function calculates energy and force of nonbonded interatons in YFF0. The flag defines if energy(forces) adds of subtracts
//Note g_vecs and r_vecs migh be ffsys->g and ffsys->r correspondingly
inline double calc_atom__nb_yff0(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top);
inline double calc_atom__nb_enrg_yff0(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top);

//This function calculates nonbonded vdW energy in yff0. Do not call it for hydrogens!
inline double calc_atom__nb_vdW_enrg_yff0(register double _r,t_top *top);
//This function calculates nonbonded coulomb energy in yff0
inline double calc_atom__nb_coulomb_enrg_yff0(register double _Q,register double _r,t_top *top);

//------------------------------------    Y F F 0    G L O B A L     H A M I L T O N I A N   ------------------------------------------------

//This function calculates energy in the active layer and marks atoms to exclude from grid summing
inline void calc_nbenergy_yff0(double *e,t_ffsys *ffsys,t_top *top);
//Note it scales repulsion as E*=E*(3*x^2-2*x^3) in range 0 @ px==0. ... E @ px=1. dE*/dr=(3*x^2-2*x^3)*dE/dr
//This function calculates energy in the active layer and marks atoms to exclude from grid summing
inline void calc_nbgrad_yff0(double *e,t_ffsys *ffsys,t_top *top);

//!!!!!!!!!!!!! This function can be improved by excluding the exact anchor in explicite summation of ongrid failure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculates yff1 energy on dgrid
//It uses marks made by calc_yff1_grad(t_ffsys *ffsys,t_top *top) to skip frozen atoms
inline void calc_energy_on_grid_yff0(double *e,t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***Q);
//!!!!!!!!!!!!! This function can be improved by excluding the exact anchor in explicite summation of ongrid failure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculates yff1 energy on dgrid
//It uses marks made by calc_yff1_grad(t_ffsys *ffsys,t_top *top) to skip frozen atoms
inline void calc_grad_on_grid_yff0(double *e,t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***Q);

//This function calculates yff0 energy on grid
inline double calc_ffsys_nbenergy_on_dgrid_yff0(t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *Q);
//This function calculates yff0 energy gradient on grid
inline double calc_ffsys_nbgrad_on_dgrid_yff0(t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *Q);

//-----------------------------   Y F F 1    P A R T   ---------------------------------------------

//-----------------------------   G R I D    P A R T   ---------------------------------------------

//This function calculates values of yff1 A exponent
double inline calc_ffgrid_yff1_A(t_vec *r,t_ffsys *ffsys,t_top *top);
//This function calculates protomap element for yff1 A exponent
char calc_ffgrid_yff1_protoA(double *f,double *dfdx,double *dfdy,double *dfdz,double *d2fdxdy,double *d2fdxdz,double *d2fdydz,double *d3fdxdydz,double rr,t_vec *r,char label,va_list stack);
//This function calculates values of yff1 B exponent
double inline calc_ffgrid_yff1_B(t_vec *r,t_ffsys *ffsys,t_top *top);
//This function calculates protomap element for yff1 B exponent
char calc_ffgrid_yff1_protoB(double *f,double *dfdx,double *dfdy,double *dfdz,double *d2fdxdy,double *d2fdxdz,double *d2fdydz,double *d3fdxdydz,double rr,t_vec *r,char label,va_list stack);
//This function calculates values of yff1 electrostatic energy
double inline calc_ffgrid_yff1_Q(t_vec *r,t_ffsys *ffsys,t_top *top);
//This function calculates protomap element for yff1 electrostatic
char calc_ffgrid_yff1_protoQ(double *f,double *dfdx,double *dfdy,double *dfdz,double *d2fdxdy,double *d2fdxdz,double *d2fdydz,double *d3fdxdydz,double rr,t_vec *r,char label,va_list stack);

//This function calculates energy and gradient for atom on tricubic grid and return TRUE otherwise it returns FALSE 
char calc_grad_yff1_on_threecubic_tcgrid(unsigned int a_i,double q_i,t_vec *r_i,t_vec *g_i,double *e,t_vec *ori,t_len *len,double sp,double (***A)[64],double (***B)[64],double (***Q)[64],t_top *top);
//This function calculates energy and gradient for atom on dgrid and return TRUE otherwise it returns FALSE 
char calc_energy_yff1_on_tricubic_dgrid(unsigned int a_i,double q_i,t_vec *r_i,double *e,t_vec *ori,t_len *len,double sp,double ***A,double ***B,double ***Q,t_top *top);
//This function calculates energy and gradient for atom on dgrid and return TRUE otherwise it returns FALSE 
char calc_grad_yff1_on_tricubic_dgrid(unsigned int a_i,double q_i,t_vec *r_i,t_vec *g_i,double *e,t_vec *ori,t_len *len,double sp,double ***A,double ***B,double ***Q,t_top *top);

//-----------------------------   B O N D E D   ---------------------------------------------------   

//This function calculates energy and forces of bond
inline double _calc_gbond_gaff(double r,double *de,t_ff_b *ff_b);
//This function calculates energies and forces of angle
inline double _calc_gangle_gaff(double csA,double *de,t_ff_g *ff_g);
//This function calculates energies and forces of improper angle
inline double _calc_gimpr_gaff(double csF,double snF,double *de,t_ff_i *ff_i);
//This function calculates energies and forces of torsion
inline double _calc_gdihs_gaff(double csF,double snF,double *de,t_ff_d *ff_d);
inline double _calc_gdihs_enrg_gaff(double csF,double snF,t_ff_d *ff_d);

//This function calculated bonded interactions for whole molecule
inline double calc__bgrad_gaff(t_vec *g,t_vec *r,t_mol *mol,t_activem *activem);
inline void _calc__bgrad_gaff(double e[N_ESLICES],t_vec (*g)[N_ESLICES],t_vec *r,t_mol *mol,t_activem *activem);

//This function calculates torsional energy for internal coordinates gradients
inline double calc_mol_tbenergy_gaff(t_vec *r,t_mol *mol,t_activem *activem);
inline void _calc_mol_tbenergy_gaff(double e[N_ESLICES],t_vec *r,t_mol *mol,t_activem *activem);
inline double calc_mol_tbgrad_gaff(t_vec *g,t_vec *r,t_mol *mol,t_activem *activem);
inline void _calc_mol_tbgrad_gaff(double e[N_ESLICES],t_vec (*g)[N_ESLICES],t_vec *r,t_mol *mol,t_activem *activem);

//------------------------------------    N O N B O N D E D   ------------------------------------------

//This function calculates energy and force of nonbonded interatons in YFF1. The flag defines if energy(forces) adds of subtracts
//Note g_vecs and r_vecs migh be ffsys->g and ffsys->r correspondingly
inline double calc_atom__nb_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top);
inline double calc_atom__nb_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top);
//This function calculates zero-level yff approximation
inline double calc_atom__znb_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top);
inline double calc_atom__znb_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top);
//This function calculates 1-4 coulomb interactions only
inline double calc_atom__14_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top);
inline double calc_atom__14_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top);
//This function calculates negative difference between 1-4 and normal coulomb interactions only
inline double calc_atom_n14_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top);
inline double calc_atom_n14_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top);
//This function calculates potential on grid
inline double calc_atom__znb_yff1_on_grid(double a,double b,double q,t_vec *da,t_vec *db,t_vec *dq,unsigned int a_i,double q_i,t_vec *g_i,t_top *top);
inline double calc_atom__znb_enrg_yff1_on_grid(double a,double b,double q,unsigned int a_i,double q_i,t_top *top);

//------------------------------------   Y F F 1     G L O B A L     H A M I L T O N I A N    ------------------------------------------------

//This function calculates torsions energy
inline void calc_torss_energy_yff1(double e[N_ESLICES],t_ffsys *ffsys);
//This function calculates torsions energy and its derivative
inline void calc_torss_grad_yff1(double e[N_ESLICES],t_ffsys *ffsys);

//!!!!!!!!!!!!! This function can be improved by excluding the exact anchor in explicite summation of ongrid failure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculates yff1 energy on dgrid
//It uses marks made by calc_yff1_grad(t_ffsys *ffsys,t_top *top) to skip frozen atoms
inline void calc_energy_on_grid_yff1(double e[N_ESLICES],t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***B,double ***Q);
//This function calculates yff1 energy on dgrid
//It uses marks made by calc_yff1_grad(t_ffsys *ffsys,t_top *top) to skip frozen atoms
inline void calc_grad_on_grid_yff1(double e[N_ESLICES],t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***B,double ***Q);

//This function calculates yff1 energy on grid
inline double calc_ffsys_nbenergy_on_dgrid_yff1(t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *B,t_dgrid *Q);
//This function calculates yff1 energy on grid
inline double calc_ffsys_nbgrad_on_dgrid_yff1(t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *B,t_dgrid *Q);


//This function update rotational coordinates to be always withing radii<=PI sphere
char sync_ic_coords(double *_e,unsigned int n,double *x,va_list stack);

//--------------------------     L I N E A R     C O O R D I N A T E S    P A R T   -----------------------------------------------------


//This function gather statistics for IC system
void get_ic_complexity(unsigned int *n,unsigned int *m,unsigned int *natoms,t_ffsys *ffsys);

//This function update rotational coordinates to be always withing radii<=PI sphere
char sync_ic_coords(double *_e,unsigned int n,double *x,va_list stack);

//This function initializes ic massives for use in functions below
char initialize_ic_minimizer(unsigned int method,unsigned int nwats,unsigned int *n,unsigned int *m,unsigned int *naatoms,
                             t_tensor **R,t_vec **cm,t_vec **rvec,double *x[2],double *g[2],double ***G,double **p,t_ffsys *ffsys);

//This function minimizes system in internal coordinates on the grid
char optimize_ic_ffsys_on_grid_yff0(double *e,unsigned int method,unsigned int nsteps,double tol,unsigned int n,double **x,double **g,double **G,double *p,unsigned int m,t_ffsys *ffsys,t_dgrid *A,t_dgrid *Q,t_vec *rvecs,t_vec *cm,t_tensor *R,t_top *top);
//This function minimizes system in internal coordinates on the grid
char optimize_ic_ffsys_on_grid_yff1(double *e,unsigned int method,unsigned int nsteps,double tol,unsigned int n,double **x,double **g,double **G,double *p,unsigned int m,t_ffsys *ffsys,t_dgrid *A,t_dgrid *B,t_dgrid *Q,t_vec *rvecs,t_vec *cm,t_tensor *R,t_top *top);

