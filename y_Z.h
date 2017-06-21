//This file contains charge calculations plugin
#define Y_Z 0x1
#ifndef Y_MOL
#include "y_mol.h"
#endif
#ifndef Y_DELAUNAY
#include "y_delaunay.h"
#endif
#ifndef Y_GEOMETRY
#include "y_geometry.h"
#endif
#ifndef Y_MATRIX
#include "y_matrix.h"
#endif
#ifndef  Y_FFSYS
#include "y_ffsys.h"
#endif



/*
typedef struct{
              unsigned int       size;         //Number of rotation units
              unsigned int    *aedges;         //Anchors joining edge
              t_vec        (*u0)[0x2];         //The vector of rotation axis base
              t_vec              *ru0;         //The vector to rotation axes base (joining middle point)
              t_list *(*anchors)[0x2];         //The anchors list that belong to current rotation unit
              }t_runits;

//This function builds internal coordinates rotation tree for single fragment based on neighbouring and temp massive
unsigned int build_rotation_tree_with_neighbors(t_list *neighbors,unsigned int *step,unsigned int *stack);

//This function free memory from runits structure
void free_runits(t_runits *runits);

//This function defines rotation units for given molecule.
t_runits *define_runits(t_mol *mol);

//THis function defines rotation units for given list of given molecule.
//t_runits *define_runits_with_list(t_list *list,t_mol *mol);

//This function define vector into runits
inline void calc_runits(register t_vec *r,register t_runits *runits,register t_mol *mol);

//This function calculates internal distribution functions
double calc_internal_distribution(t_runits **runits,t_ffsys *ffsys);
*/

//This function calculates mass of given anchor fragment
inline void calc_anchors_cm(t_vec *cm,t_vec *r,unsigned int anchor_id,t_mol *mol,t_top *top);

//This function move system to the molecule cm
inline void calc_mols_cm(t_vec *cm,t_vec *r_vecs,t_mol *mol,t_top *top);

//This function calculates external distribution functions for 6D dimensions
double calc_external_distribution(double *e,double Z[0x6][0x6],unsigned int mol_id,t_ffsys *ffsys,t_top *top);

