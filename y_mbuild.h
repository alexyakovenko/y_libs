#define Y_BUILD 1

//This file keeps routines for 3D molecules binding using parameterized topology

#ifndef Y_MOL
#include "y_mol.h"
#endif

#define MAX_NEIGHBORS 4

//This function resolves coordinates of onevalent atoms (mostrly hydrogens)
//Note. It does NOT perform conformation search
unsigned int resolve_onevalent(t_clist *neighbors,t_vec *r,t_mol *mol);

