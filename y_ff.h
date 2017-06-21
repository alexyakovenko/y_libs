#define Y_FF 0x1

#ifndef Y_MOL
#include "y_mol.h"
#endif
#ifndef Y_CHARGE
#include "y_charge.h"
#endif
#ifndef Y_INTERPOLATE
#include "y_minimize.h"
#endif

//------------------------------------    B O N D E D   P R I M I T I V E S   ------------------------------------------

//This function calculates bond energy
inline double _calc_benergy_yff1(double K,double V,t_vec *ri,t_vec *rj);
//This function calculates bond energy and its gradient
inline double _calc_bgrad_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *gi,t_vec *gj);
//This function calculates angle energy
inline double _calc_genergy_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk);
//This function calculates angle energy and its gradient
inline double _calc_ggrad_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *gi,t_vec *gj,t_vec *gk);
//This function calculates improper energy
inline double _calc_ienergy_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl);
//This function calculates improper energy and its gradient
inline double _calc_igrad_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl,t_vec *gi,t_vec *gj,t_vec *gk,t_vec *gl);
//This function calculates improper energy with respect to circular rotation
inline double _calc_ienergy_PI_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl);
//This function calculates improper energy and its gradient
inline double _calc_igrad_PI_yff1(double K,double V,t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl,t_vec *gi,t_vec *gj,t_vec *gk,t_vec *gl);
//This function calculates torsion energy
inline double _calc_denergy_yff1(double K[SIZE_DIH],double V[SIZE_DIH],double N[SIZE_DIH],t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl);
//This function calculates torsion energy and its gradient
inline double _calc_dgrad_yff1(double K[SIZE_DIH],double V[SIZE_DIH],double N[SIZE_DIH],t_vec *ri,t_vec *rj,t_vec *rk,t_vec *rl,t_vec *gi,t_vec *gj,t_vec *gk,t_vec *gl);

//This function summs energy of bonded interations of all atoms of the molecule
inline double summ_bnmol_energy_yff1(t_mol *mol,t_vec *r);

//This function summs energies and gradients of bonded interations of all atoms of the molecule
inline double summ_bnmol_grad_yff1(t_vec *g,t_mol *mol,t_vec *r);


//------------------------------------    N O N B O N D E D   P R I M I T I V E S   ------------------------------------------


//This function calculates energy of classic nonbonded interatons in YFF1. 
inline double _calc_atom__nbenergy_yff1(register double Aij,register double Bij,double qi,double qj,t_vec *ri,t_vec *rj);

//This function summs nonbonded interations energies of all anchors of the same molecule
inline double summ_nbatom_mol_energy_yff1(t_mol *mol,t_vec *r,double **A,double **B);
//This function summs nonbonded interations energies and gradients of all anchors of the same molecule
inline double summ_nbmol_grad_yff1(t_vec *g, t_mol *mol, t_vec *r, double **A, double**B);


//------------------------------------    U T I L I T I E S   ------------------------------------------


//This function calculates energy and its gradient for a stand alone mol
double calc_mol_energy_yff1(t_mol *mol,t_vec *rvecs,double **A,double **B);
//This function calculates energy and its gradient for a stand alone mol
char calc_mol_grad_yff1(double *e,unsigned int n,double *x,double *g,double **G,va_list stack);

//This function perform energy minimization of stand alone molecule
char optimize_mol(unsigned int nsteps,double tol,t_vec *rvecs,char *label,t_mol *mol,double **A,double **B);
//The same as previous but can retraint subset of atomic coordinates
char optimize_mol_restraints(unsigned int nsteps,double tol,t_vec *rvecs,char *label,t_mol *mol,double **A,double **B,t_list *restraints,double *restraints_k);

#define AMIDE_BONDS_RESTRAINTS_K 4000.

//This function check molecular structure for evident geometrical errors
//So far it checks:
// 1. If any amide bonds in molecule A is in cys conformation. If so it swaps H and X atoms around N~CR=O atom.
// 2. If any <=6 cycle in molecule A is pierced by a bond of molecule B (A can be equal to B). If so the bond forming atoms are translated to the closes edge +1 A away.
unsigned int check_awful_geometry_errors(unsigned int nsteps,char (*ordera)[4],t_clist *neighborsa,t_vec *ra,t_mol *mola,char (*orderb)[4],t_clist *neighborsb,t_vec *rb,t_mol *molb,double **A,double **B);



