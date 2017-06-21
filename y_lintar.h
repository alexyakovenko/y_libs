// Header contains routines for arithmetical operations with long fixed-point
// numbers represented as int arrays, as well as some auxiliary functions.

#define Y_LINTAR 0x1

#ifndef Y_SYS 
#include "y_system.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif

#define SIGN(A) (((A) > 0) - ((A) < 0))

//===============================-GENERAL-======================================

// Aux constants for double representations in x86
// Quantity of significant binary digits in double mantissa including the
// hidden bit.
#define DOUBLE_SIGNIF_DIGITS 53
// Exponent bias.
#define DOUBLE_EX_BIAS 1023
// Quantity of bits that store exponent.
#define DOUBLE_EX_DIGITS 11

// It is assumed, that A[0] of array int *A contains length of a number (each
// int stores up to 30 binary digits), signed with number sign.
#define LINTAR_DIG_CONST 28

// Position of decimal point (if relevant) is passed into functions explicitly.
// To convert data to another precision format, use convert_precision().

//========================-MEMORY AND CONVERSIONS-==============================

// Converts double into long number.
// Call allocates memory for l, old value is NOT freed and is lost.
char double2lint(int **l, double a, int precision);

// Converts long number into double.
double lint2double(int *l, int precision);

// Converts long number in another precision format, saving value or rounding
// if needed. If conversion is to larger array size, memory is allocated anew.
// NB! In that case, free is called on pointer l!
// NB! Pointer of pointer in first argument.
char convert_precision(int **l, int former_p, int target_p);

//================================-COMPARISON-==================================

// In comparison functions precision is assumed to be the same for both numbers.

// Compares absolute values of two numbers. Returns 1 if |a| > |b|, -1 if
//|a| < |b| and 0 if |a| = |b|.
char us_lint_com(int *a, int *b);

// Compares two numbers. Returns 1 if a > b, -1 if a < b and 0 if a = b.
char lint_com(int *a, int *b);

//==========================-ADDITIVE ARITHMETICS-==============================

// In additive functions precision is assumed to be the same for both arguments
// and result.

// Unsigned plus r = |a| + |b|.
// Call allocates memory for r, old value is NOT freed and is lost.
// a and b may point to same memory.
char us_lint_add(int **r, int *a, int *b);

// Unsigned minus r = ||a| - |b||.
// Call allocates memory for r, old value is NOT freed and is lost.
// a and b may point to same memory.
char us_lint_sub(int **r, int *a, int *b);

// Adds two numbers r = a + b.
// Call allocates memory for r, old value is NOT freed and is lost.
// a and b may point to same memory.
char lint_add(int **r, int *a, int *b);

// Subtracts two numbers r = a - b.
// Call allocates memory for r, old value is NOT freed and is lost.
// a and b may point to same memory.
char lint_sub(int **r, int *a, int *b);

//=========================-MULTIPLICATIVE ARITHMETICS-=========================

// Multiplies numbers r = a * b using Karatsuba's recursive algorithm.
// Numbers consisting of 4 or less ints are not divided further, but multiplied
// straightly.
// Function allocates and frees necessary memory only once.
// Precision is relevant, so it is passed in the function. It is applied both
// to arguments and to result.
// Call allocates memory for r, old value is NOT freed and is lost.
// a and b may point to same memory.
char lint_fast_mult(int **r, int *a, int *b, int precision);

//===========================-DETERMINANT SIGNS-================================

// In determinant functions precision is relevant, so it is passed in.
// It is applied both to arguments and to result.

// Evaluates determinant 2x2.
// Call allocates memory for r, old value is NOT freed and is lost.
// aij may point to same memory.
void lint_det2(int **r, int *a00, int *a01,
                        int *a10, int *a11, int precision);

// Checks if sign of 2-determinant is more or less than zero.
// Uses SoS (always returns unambiguous and determined result).
char lint_SoS_udet2(double a00,
                    double a10, int precision);

// Checks if sign of 4-determinant is more or less than zero.
// Uses SoS (always returns unambiguous and determined result).
char lint_SoS_det2(double a00, double a01,
                   double a10, double a11, int precision);

// Evaluates determinant 3x3.
// Call allocates memory for r, old value is NOT freed and is lost.
// aij may point to same memory.
void lint_det3(int **r, int *a00, int *a01, int *a02,
                        int *a10, int *a11, int *a12,
                        int *a20, int *a21, int *a22, int precision);

// Checks if sign of 6-determinant is more or less than zero.
// Uses SoS (always returns unambiguous and determined result).
char lint_SoS_udet3(double a00, double a01,
                    double a10, double a11,
                    double a20, double a21, int precision);

// Checks if sign of 9-determinant is more or less than zero.
// Uses SoS (always returns unambiguous and determined result).
char lint_SoS_det3(double a00, double a01, double a02,
                   double a10, double a11, double a12, 
                   double a20, double a21, double a22, int precision);

// Evaluates determinant 4x4.
// Call allocates memory for r, old value is NOT freed and is lost.
// aij may point to same memory.
void lint_det4(int **r, int *a00, int *a01, int *a02, int *a03,
                        int *a10, int *a11, int *a12, int *a13,
                        int *a20, int *a21, int *a22, int *a23,
                        int *a30, int *a31, int *a32, int *a33, int precision);

// Checks if sign of 12-determinant is more or less than zero.
// Uses SoS (always returns unambiguous and determined result).
char lint_SoS_udet4(double a00, double a01, double a02,
                    double a10, double a11, double a12,
                    double a20, double a21, double a22,
                    double a30, double a31, double a32, int precision);

// Checks if sign of 16-determinant is more or less than zero.
// Uses SoS (always returns unambiguous and determined result).
char lint_SoS_det4(double a00, double a01, double a02, double a03,
                   double a10, double a11, double a12, double a13,
                   double a20, double a21, double a22, double a23,
                   double a30, double a31, double a32, double a33, int precision);

// Evaluates determinant 5x5.
// Call allocates memory for r, old value is NOT freed and is lost.
// aij may point to same memory.
void lint_det5(int **r, int *a00, int *a01, int *a02, int *a03, int *a04,
                        int *a10, int *a11, int *a12, int *a13, int *a14,
                        int *a20, int *a21, int *a22, int *a23, int *a24,
                        int *a30, int *a31, int *a32, int *a33, int *a34,
                        int *a40, int *a41, int *a42, int *a43, int *a44, int precision);
                        
// Checks if sign of 20-determinant is more or less than zero.
// Uses SoS (always returns unambiguous and determined result).
char lint_SoS_udet5(double a00, double a01, double a02, double a03,
                    double a10, double a11, double a12, double a13,
                    double a20, double a21, double a22, double a23,
                    double a30, double a31, double a32, double a33,
                    double a40, double a41, double a42, double a43, int precision);

// Checks if sign of 25-determinant is more or less than zero.
// Uses SoS (always returns unambiguous and determined result).
char lint_SoS_det5(double a00, double a01, double a02, double a03, double a04,
                   double a10, double a11, double a12, double a13, double a14,
                   double a20, double a21, double a22, double a23, double a24,
                   double a30, double a31, double a32, double a33, double a34,
                   double a40, double a41, double a42, double a43, double a44, int precision);
