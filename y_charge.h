//This file contains charge calculations plugin
#define Y_CHARGE 0x1
#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_MATRIX
#include "y_matrix.h"
#endif
#ifndef Y_MINIIZE
#include "y_minimize.h"
#endif
#ifndef Y_MOL
#include "y_mol.h"
#endif


//Coulomb constant
#define COULOMB_K  1389.35486               //kJ*mol^-1*A*e^2
//Dielectric constant for water
#define ESOLVENT 80.00
//Water radii (offset)
#define RSOLVENT 1.4

//Size of problem to switch from TDM to SPM
#define COMPLEXITY_LIMIT_CDS_TDM 0xF
//#define COMPLEXITY_LIMIT_CDS_TDM 0x4 //Debuging


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
// The derivative is calculated as: dq/db[j]=(dIH/db[j]).x=-IH.(dH(b)/db).IH.x=-IH.u(X)u^T.IH.x=-c(X)c^T.x=-_d*c, where c=IH.u and _d={c^T.x}
// and Err=Summ[i](q[i]-q0[i])^2 so dErr/db[j]=SUMM[i](2*(q[i]-q0[i])*dq[i]/db[j])


/************************************************   D E N S E   C D S   ***********************************************/ 

//This function calculates molecules charge. The inverted operator of internal electronic strucutre of a molecule is stored as a real dense matrix.
//NB! The mol need to possess correct natoms, nvatoms, vatoms[], engs[], ^charges[], nedges, edges[], nvedges, vedges[], ^hrds[] massives
char calc_ionized_CDS_cholesky_tdm(unsigned int nsteps,double tol,t_mol *mol);

/**************************************   S P A R S E    C D S   **************************************************************/

//This function calculates molecules charge. The inverted operator of internal electronic strucutre of a molecule is stored as a real dense matrix.
//NB! The mol need to possess correct natoms, nvatoms, vatoms[], engs[], ^charges[], nedges, edges[], nvedges, vedges[], ^hrds[] massives
char calc_ionized_CDS_cholesky_spm(unsigned int nsteps,double tol,t_mol *mol);

/*******************************   A   D E F A U L T   C D S   W R A P P E R   F U N C T I I O N   ********************************************/

//This function calculates molecules charge. The inverted operator of internal electronic strucutre of a molecule is stored as a real dense matrix.
char calc_Oliferenko_ionized_CDS_cholesky_default(unsigned int nsteps,double tol,t_mol *mol,t_top *top);

/*******************************************************   C H A R G E    G R O U P S   ************************************************************************/

//This function, IDEALLY, should edit charges and define charged groups tp have a molecules' electric field be represented by a set of small zero-charged fragments.
//Note. I am doubt if this is possible at all and particularly have no ideas of how to implement it.
//The current plan is to keep the true best charges on atoms but split the molecule so to reduce the 'artificial' dipole errors. 
//Stage I. Move all in the anchor into the charged group
//Stage II. Join H-X-C into one groups
//That is it so far...
unsigned int *split_charges(t_clist *neighbors,unsigned int *anchor_id,t_mol *mol,t_top *top);

/******************************************************* D E S O L V A T A T I O N *******************************************************************************/

//This function calculates the solvatation energy of the system (accordingly with Robert S. Rizzo, Tiba Aynechi, David A. Case, Irwin D. Kunts "Estimate of Absolute Free Energies of Gidratation Using Continium Methods: Accuracy of Partial Charge Models and Optimization of Nonpolar Contribution", JCTC 2006 2 128-139  )
double calc_esolv(double *sasa,double *sava,t_list *it,t_mol *mol,t_top *top);

//This function calculates the Born radiuses (accordingly with Di Qiu,Peter S. Shenkin,Frank P. Hollinger, and W. Clark Still "The GB/SA Continuum Model for Solvation. A Fast Analitical Method for the Calculation of Aproximate Born Radii" J.Phys.Chem. 1997, 101,3005-3014)
//Note. Vj is now replaced with Solvent Acessible Volume Area that calculated from delaunay triangulation.
void calc_Born_radii(double *rBorn,double *dij2,double *sava,t_list *it,t_mol *mol,t_top *top);


