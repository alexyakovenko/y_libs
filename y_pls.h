//This module contain PLS routine operatins
//It is mainly based on Sijmen de Jong, Barry M. Wise and N.Lawrence Ricker "Canonical partial least squares and continuum power regression" J.Chemometrics, 15, 85-100, 2001.

//Note. this module deals with 3D operation. If you are interested in somthing more complex check y_vec.c/y_vec.h
#define Y_PLS 0x1

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
#ifndef Y_MATRIX
#include "y_matrix.h"
#endif
#ifndef Y_PLIBS
#include "y_plibs.h"
#endif
#ifndef Y_PMATRIX
#include "y_pmatrix.h"
#endif

//This function calculates PLS model serially. flag_T determines A or A^T is submitted.
char calc_serial_pls_model(unsigned int power,unsigned int ni,unsigned int nj,double **A,char flag_T,double *x,double *y);

