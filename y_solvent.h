//This module take care about solvation of molecules with water

#define Y_SOLVENT 0x1


#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_LIST
#include "y_list.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif
#ifndef Y_VECTOR
#include "y_vector.h"
#endif
#ifndef Y_GEOMETRY
#include "y_geometry.h"
#endif
#ifndef Y_DELAUNAY
#include "y_delaunay.h"
#endif
#ifndef Y_MOL
#include "y_mol.h" 
#endif
#ifndef Y_CHARGE
#include "y_charge.h"
#endif

#define YSOLVENT_MAX_WAT_MOLS_IN_DB   5 //Maximal amount of water molecues in a configuration
#define MAX_NWATS  YSOLVENT_MAX_WAT_MOLS_IN_DB

//Water atoms charges
#define YSOLVENT_FF_qH    +0.360368
#define YSOLVENT_FF_qO    -0.720736

//This module construct explicit solvent shell using precomputed DB of waters configurations

//A water configuration
typedef struct{
              unsigned int       n_E;                   //Amount of external E vectors per configuration
              t_vec               *E;                   //Values of external vector 
              double              *e;                   //Potential energy of the configuration
              unsigned int       n_r;                   //Amount of water molecules coordinates in configuration
              t_vec               *r;                   //Water's coordinates
              t_vec            th[4];                   //The tetrahedron that inscribes a configuration (vertices th[0] and th[1] has the longest edge)
              double        d2[4][4];                   //Matrix of distances between thetrahedron vertices  
              }t_wat_cfg;

//Water configurations db
typedef struct{
              double          qO, qH;                   //Charges of water's oxigen and hydrogen 
// --------------------------- M E M O R Y    M A N A G E M E N T -----------------------------------------------//          
              unsigned int      _n_E;                   //Total amount of E vectors
              t_vec              *_E;                   //Massive of E vectors
              double             *_e;                   //Massive of potential energies
              unsigned int      _n_r;                   //Total amount of water molecule's coordinates
              t_vec              *_r;                   //Massive of water coordinates   
// ------------------------------------ D B    I T S E L F -----------------------------------------------------//          
              unsigned int    n_cfgs;                   //Amount of configurations in DB
              t_wat_cfg        *cfgs;                   //Water configurations                  
              }t_wat_db;

//-------------------------------------------------- D T R    F I L L E R S --------------------------------------------------------//

//This function allocate memory for t_wat_db
//NOTE. t_wat_db iw simply destroyed with free(*it)
t_wat_db *alloc_wat_db(unsigned int n_cfgs,unsigned int _n_E,unsigned int _n_r);

//This function reads wat_db from file
//DB file structure
// [ [ n_cfgs ], [ _n_E ], [ _n_r ] ]
// [ [ n_E ], [ n_r ], [ th ], [ d2 ] ] x n_cfgs
// [ _E ] x _n_E, [ _r ] x _n_r
t_wat_db *read_wat_db(FILE *in);

//This function write wat_db to HDD 
char write_wat_db(FILE *out,t_wat_db *w_db);

//This function imports y_wat.db from it's standard location at ../top/
t_wat_db *import_wat_db();

#define free_wat_db(w_db)   free(w_db) 

//---------------------------------------------------- F I L L E R S ----------------------------------------------------------------//


//This function fill a union of ths with a single water
unsigned int fill_th_union(t_vec *rvecs,t_vec *th_cm,t_vec *E,t_wat_db *w_db);

//This function fills a thetrahedron with water using empirical DB
//NOTE. th define rvecs, r stores results (should be of MAX_WATS reserved) and returns amount of added molecules.
//NB! thetrahedrons should have prechecked volume, here we just find appropriate fillings.
unsigned int fill_th(t_vec *rvecs,unsigned int th[4],t_lvec *lvecs,t_vec *E,t_wat_db *w_db);

//This function solvates complex explicitely
//It solvates thetragedrons sorted by their E. theragedrons are filled using db or by gouping witha common facet or common edge of sufficient lenght.
//NOTE. This function does max_nwats+/-biggest ensable from db
char solvate_explicit(unsigned int *nwats,unsigned int max_nwats,double min_th_ir,double min_th_vol,double max_th_vol,t_list *set,unsigned int n,t_vec *rvec,double *q,t_delaunay *dtr,double *radii,t_wat_db *w_db);


//-------------------------------------- M O L E C U L E S     O P E R A T I O N    R O U T I N E S -------------------------------------//


