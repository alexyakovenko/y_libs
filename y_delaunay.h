#define Y_DELAUNAY 0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_VECTOR
#include "y_vector.h"
#endif
#ifndef Y_LIST
#include "y_list.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif
#ifndef Y_GEOMETRY
#include "y_geometry.h"
#endif


//This structure describes Delaunay triangulation
typedef struct{
              unsigned int   size;                   //The number of vertices in superstructure
              unsigned int  _size;                   //The amount of allocated memory for vertices in superstructure
              t_lvec        *lvec;                   //The nodes lifter coordinates
              unsigned int   size_th;                //The number of tetrahedrons in current triangulation
              unsigned int  _size_th;                //The amount of allocated memory for tetrahedrons
              unsigned int  (*th)[0x4];              //The tetrahedrons of triangulation
              unsigned int (*nth)[0x4];              //The neighbours associated with tetrahedron
              }t_delaunay;


#ifndef Y_MOL
#include "y_mol.h"
#endif

//This function uploads Delaunay
t_delaunay *read_delaunay(FILE *in);

//This function saves Delaunay
char write_delaunay(FILE *out,t_delaunay *d);

//This function looks for the tetrahedron that contain current point in and return its id (direction algorithm is applied)
char ffind_tetrahedron(t_lvec *p,unsigned int *th,t_delaunay *d);

//This function complete Delaunay triangulation. It require the regular delaunay of certain fragment and number of points to be triangulated.
char triangulate_3d_delaunay(unsigned int nsize,t_delaunay *d);

//This function construct Delaunay triangulation structure
t_delaunay *alloc_delaunay(unsigned int size,unsigned int size_th);

//This function clone delaunay triangulation extended with size centers and size_th tetrahedrons
t_delaunay *clone_delaunay(unsigned int size,unsigned int size_th,t_delaunay *d);

//This ia a destuctor for Delaunay triangulation
void free_delaunay(t_delaunay *d);

//This function show tetrahedrons list of triangulation
inline void show_delaunay(t_delaunay *d);

//This function builds 3d triangulation on th set of given points
char triangulate_3d_set(t_delaunay *d);

//This function performs delaunay triangulation via mol interface
t_delaunay *triangulate_3d_mol(t_list *it,t_vec *r,t_mol *mol,t_top *top);

//This function performs delaunay triangulation via mol interface
char add_mol_to_triangulation(t_list *it,t_vec *r,t_mol *mol,t_top *top,t_delaunay *d);

//---------------------------------- T E T R A HE D R O N S   P A R T --------------------------------

//This function calculates tetrahedron volume
//Note. These function returns 'signed' volume
double calc_th_volume(t_vec *a,t_vec *b,t_vec *c,t_vec *d); 
double calc_thetrahedron_volume(t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d);

//This function calculates the center of mass of solid gomogenic tetrahedron as 1/3 of the line that connect a mediane intersection of triangle to opposite vertice
void calc_th_cm(t_vec *cm,t_vec *a,t_vec *b,t_vec *c,t_vec *d);

//This function calculates intrinsic radii of sphere between four other spheres and returns signed volume of tetrahedron
//E.N.Sickafus and Neil A. Mackie "Interstitial space in hard-shere clusters" Acta Cryst. A30, 850-851 1974
//Note if tetrahedrons volume less than cutoff calculations are not done and 0.0000 is returned
double calc_thetrahedron_insphere(double cutoff,t_lvec *insphere,t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d);

//This empirical function for evaluation of tetrahedron inscribed sphere
//Note. This function inscribes cyrcle in tetrahedron bisector plane triangle that is opposite to longest edge of tetrahedron
double embed_sphere_in_tetrahedron(t_vec *o,t_vec *a,t_vec *b,t_vec *c,t_vec *d);


