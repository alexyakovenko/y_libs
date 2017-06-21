#include "y_lintar.h"

#define YERROR_NO_ERRORS YERROR_OK

char double2lint(int **l, double a, int precision)
{
double t;
int _i, sign, ndigs, ex, target, size;
unsigned long long int man;

if (!(isfinite(a)) || (precision < 0))
  return YERROR_DATA_FORMAT;                   // invalid input

sign = (a > 0 ? 1 : -1);                       // extracting sign

a = fabs(a);
if (a < SMALL2)
  {
  *l = malloc(sizeof(int) * 2);
  (*l)[0] = 0;                                 // zero or small input
  return YERROR_NO_ERRORS;
  }

t = frexp(a, &ex);
man = fabs(t) * pow(2, DOUBLE_SIGNIF_DIGITS);
ex--;
if (ex < -precision * LINTAR_DIG_CONST)
  {
  *l = malloc(sizeof(int) * 2);
  (*l)[0] = 0;                                 // small input
  }
size = (ex / LINTAR_DIG_CONST + precision + (ex < 0 ? 0 : 1));
*l = malloc(sizeof(int) * (size + 2));
(*l)[0] = size * sign;
for (_i = 1; _i <= abs((*l)[0]); ++_i)
  (*l)[_i] = 0;
_i = DOUBLE_SIGNIF_DIGITS; 
while (_i > 0)
  {
  if (_i >= LINTAR_DIG_CONST)
    {
    ndigs = ex % LINTAR_DIG_CONST + 1;
    target = ex / LINTAR_DIG_CONST + precision + 1;
    if (ndigs <= 0)
      {
      ndigs += LINTAR_DIG_CONST;
      target -= 1;
      }
    if (target > 0)
      (*l)[target] = (man >> (_i - ndigs));
    man -= (man >> (_i - ndigs)) << (_i - ndigs);
    ex -= ndigs;
    _i -= ndigs;
    }
  else
    {
    target = ex / LINTAR_DIG_CONST + precision + (ex < 0 ? 0 : 1);
    if(target > 0)
      (*l)[target] = (man << (LINTAR_DIG_CONST - _i));
    _i = 0;
    }
  }
return YERROR_NO_ERRORS;
}

double lint2double(int *l, int precision)
{
int ndigs, f, _i, ex;
double res;
unsigned long long int man;

if (l[0] == 0)
  return 0;                                    // zero length returns zero

if (precision < 0)
  return YERROR_DATA_FORMAT;

_i = abs(l[0]);                                // extracting significand and exponent
f = DOUBLE_SIGNIF_DIGITS;
ndigs = 0;
while (l[_i] >= (1LL << ndigs))
  ++ndigs;
ex = (_i - precision - 1) * LINTAR_DIG_CONST + ndigs;
man = l[_i];
f -= ndigs;
--_i;
while ((f > 0) && (_i > 0))
  {
  if (f >= LINTAR_DIG_CONST)
    {
    man = (man << LINTAR_DIG_CONST) + l[_i];
    f -= LINTAR_DIG_CONST;
    }
  else
    {
    man = (man << f) + (l[_i] >> (LINTAR_DIG_CONST - f));
    f = 0;
    }
  --_i;
  }
man <<= f;

res = pow(2, ex - DOUBLE_SIGNIF_DIGITS);
res *= man;

if (l[0] < 0)                                  // applying sign
  res *= -1;

return res;
}

char convert_precision(int **l, int former_p, int target_p)
{

int _i, _step, half;
char round_up, flaw;

if ((former_p < 0) || (target_p < 0))
  return YERROR_DATA_FORMAT;

if (former_p > target_p)
  {                                            // Rounding conversion.
  _step = former_p - target_p;
  half = (1 << (LINTAR_DIG_CONST - 1));
  if ((*l)[_step] < half)                      // Cheking first digit to be deleted.
    round_up = FALSE;
  else
    round_up = TRUE;
  for (_i = 1; _i <= abs((*l)[0]) - _step; ++_i)
    (*l)[_i] = (*l)[_i + _step];
  (*l)[0] = SIGN((*l)[0]) * (abs((*l)[0]) - _step + 1);
  (*l)[abs((*l[0]) - 1) + 1] = 0;
  if (round_up)
    ++(*l)[1];
  flaw = ((*l)[1] == (half << 1));
  _i = 1;
  while (flaw)
    {                                          // Overflow not possible due to
    (*l)[_i] = 0;                              // the fact at least one int freed during rounding.
    ++(*l)[_i + 1];
    flaw = ((*l)[_i + 1] == (half << 1)); ++_i;
    }
  if ((*l)[abs((*l[0])) + 1] == 0)
    (*l)[0] = SIGN((*l)[0]) * (abs((*l)[0]) - 1);
  }

if (target_p > former_p)
  {                                            // New allocating conversion.
  _step = target_p - former_p;
  int *newmem, alloc_size;
  char is_zero = FALSE;
  
  if ((*l)[0] == 0)
    is_zero = TRUE;
  
  if (abs((*l)[0]) + _step + 1 > target_p + 2)
    alloc_size = abs((*l)[0]) + _step + 1;
  else
    alloc_size = target_p + 2;
  newmem = malloc(sizeof(int) * alloc_size);
  for (_i = 1; _i < alloc_size; ++_i)
    newmem[_i] = 0;
  for (_i = 1 + _step; _i <= abs((*l)[0]) + _step; ++_i)
    newmem[_i] = (*l)[_i - _step];
  newmem[0] = SIGN(l[0]) * (abs((*l)[0]) + _step);
  if (is_zero)
    newmem[0] = 0;
  free(*l);
  *l = newmem;
  }

// No conversion needed.
return YERROR_NO_ERRORS;
}

char us_lint_com(int *a, int *b)
{
if (abs(a[0]) != abs(b[0]))
  if (abs(a[0]) > abs(b[0]))
    return 1;
  else
    return -1;
else
  {
  int _i;
  for (_i = abs(a[0]); _i > 0; --_i)
    if (a[_i] != b[_i])
    {                                       // comparing digits
      if (a[_i] > b[_i])
        return 1;
      else
        return -1;
    }
  }
return 0;
}

char lint_com(int *a, int *b)
{
if (a[0] * b[0] > 0)
  return SIGN(a[0]) * us_lint_com(a, b);    // same signs
if (a[0] * b[0] < 0)
  return SIGN(a[0]);                        // different signs
if (b[0] == 0)
  return SIGN(a[0]);                        // comparing with zero second operand
return -SIGN(b[0]);                         // comparing with zero first operand
}

char us_lint_add(int **r, int *a, int *b)
{
int maxl, minl, _i, *max, andconst;
char flaw;
andconst = (1 << LINTAR_DIG_CONST) - 1;

max = (abs(a[0]) > abs(b[0]) ? a : b);
maxl = (abs(a[0]) > abs(b[0]) ? abs(a[0]) : abs(b[0]));
minl = (abs(a[0]) < abs(b[0]) ? abs(a[0]) : abs(b[0]));
if (maxl == 0)
  {                                              // adding two zeros
  *r = malloc(sizeof(int) * 2);
  (*r)[0] = 0;
  return YERROR_NO_ERRORS;
  }
*r = malloc(sizeof(int) * (maxl + 3));

(*r)[1] = 0;                                     // adding common number parts
for (_i = 1; _i < minl; ++_i)
  {
  (*r)[_i] += (a[_i] + b[_i]) & andconst;
  (*r)[_i + 1] = (a[_i] + b[_i]) >> LINTAR_DIG_CONST;
  }

(*r)[minl] += (a[minl] + b[minl]) & andconst;
flaw = (((a[minl] + b[minl]) >> LINTAR_DIG_CONST) > 0);
for (_i = minl + 1; _i <= maxl; ++_i)
  {                                              // pasting longer number part
  if (flaw)
    (*r)[_i] = 1;
  else
    (*r)[_i] = 0;
  flaw = ((((*r)[_i] + max[_i]) >> LINTAR_DIG_CONST) > 0);
  (*r)[_i] = ((*r)[_i] + max[_i]) & andconst;
  }

if (flaw)
  {
  (*r)[0] = maxl + 1;                            // result length exceeds operand lengths
  (*r)[maxl + 1] = 1;
  }
else
  (*r)[0] = maxl;                                // result length equals longer operand length

return YERROR_NO_ERRORS;
}

char us_lint_sub(int **r, int *a, int *b)
{
int maxl, minl, _i, *max, *min, cmpres, modconst, *tmpmem;
char flaw;
modconst = 1 << LINTAR_DIG_CONST;

if (abs(a[0]) > abs(b[0]))
  {
  maxl = abs(a[0]);
  minl = abs(b[0]);
  }
else
  {
  maxl = abs(b[0]);
  minl = abs(a[0]);
  }
cmpres = us_lint_com(a, b);
if (cmpres == 0)
  {
  *r = malloc(sizeof(int) * 2);            // subtracting two equal numbers
  (*r)[0] = 0;
  return YERROR_NO_ERRORS;
  }
tmpmem = malloc(sizeof(int) * (maxl + 1));
if (cmpres == 1)
  {
  max = a;
  min = b;
  }
else
  {
  max = b;
  min = a;
  }

tmpmem[1] = 0;                              // subtracting common number parts
for (_i = 1; _i < minl; ++_i)
  {
  tmpmem[_i] += max[_i] - min[_i];
  if (tmpmem[_i] < 0)
    {
    tmpmem[_i + 1] = -1;
    tmpmem[_i] += modconst;
    }
  else
    tmpmem[_i + 1] = 0;
  }
tmpmem[minl] += max[minl] - min[minl];

flaw = (tmpmem[minl] < 0);                  // pasting longer number part
for (_i = minl + 1; _i <= maxl; ++_i)
  {
  tmpmem[_i] = max[_i];
  if (flaw)
    {
    --tmpmem[_i];
    tmpmem[_i - 1] += 1 << LINTAR_DIG_CONST;
    flaw = (tmpmem[_i] < 0);
    }
  }
_i = maxl;
while ((_i > 0) && (tmpmem[_i] == 0))
  --_i;
tmpmem[0] = _i;

*r = malloc(sizeof(int) * (abs(tmpmem[0]) + 2));
for (_i = 0; _i <= abs(tmpmem[0]); ++_i)
  (*r)[_i] = tmpmem[_i];

free(tmpmem);
return YERROR_NO_ERRORS;
}

char lint_add(int **r, int *a, int *b)
{
int cmpres, callres, *max, *min, _i;

if (a[0] == 0)
  {
  *r = malloc(sizeof(int) * (abs(b[0]) + 2));
  for (_i = 0; _i <= abs(b[0]); ++_i)
    (*r)[_i] = b[_i];
  return YERROR_NO_ERRORS;
  }
if (b[0] == 0)
  {
  *r = malloc(sizeof(int) * (abs(a[0]) + 2));
  for (_i = 0; _i <= abs(a[0]); ++_i)
    (*r)[_i] = a[_i];
  return YERROR_NO_ERRORS;
  }

cmpres = us_lint_com(a, b);
if (cmpres == 1)
  {
  max = a;
  min = b;
  }
else
  {
  max = b;
  min = a;
  }

if (a[0] * b[0] > 0)
  {
  callres = us_lint_add(r, a, b);            // same signs result in adding
  (*r)[0] *= SIGN(a[0]);
  return callres;
  }
else
  {
  callres = us_lint_sub(r, max, min);        // different signs result in subtracting
  (*r)[0] *= SIGN(max[0]);
  return callres;
  }
}

char lint_sub(int **r, int *a, int *b)
{
int callres, _i, sign;

if (a[0] == 0)
  {
  *r = malloc(sizeof(int) * (abs(b[0]) + 2));
  for (_i = 0; _i <= abs(b[0]); ++_i)
    (*r)[_i] = b[_i];
  (*r)[0] *= -1;
  return YERROR_NO_ERRORS;
  }
if (b[0] == 0)
  {
  *r = malloc(sizeof(int) * (abs(a[0]) + 2));
  for (_i = 0; _i <= abs(a[0]); ++_i)
    (*r)[_i] = a[_i];
  return YERROR_NO_ERRORS;
  }


if (a[0] * b[0] >= 0)
  {
  if (us_lint_com(a, b) == 1)
    sign = SIGN(a[0]);
  else
    sign = -SIGN(a[0]);
  callres = us_lint_sub(r, a, b);            // same signs result in subtracting
  (*r)[0] *= sign;
  return callres;
  }
else
  {
  callres = us_lint_add(r, a, b);            // different signs result in adding
  (*r)[0] *= SIGN(a[0]);;
  return callres;
  }
}

// Calculates total necessary memory needed for lin_karatsuba_mult to multiply
// two numberes of length 2^pow2.
int nec_mem(int pow2)
{
int _i, res;
res = 0;

for (_i = 3; _i <= pow2; ++_i)
  res = 4 * res + (1 << (_i + 2));

return res;
}

// Multiples two numbers containing 4 ints each. r = a * b.
// Neither a[0], b[0] nor r[0] does not contain signed length, but digtis are
// stored in all 4 ints (8 ints in case of r).
void straight4mult(int *r, int *a, int *b)
{
int _i, _j, modconst, andconst;
long long int tmpmem[8];
modconst = 1 << LINTAR_DIG_CONST;
andconst = modconst - 1;

for (_i = 0; _i < 8; ++_i)
  tmpmem[_i] = 0;

for (_i = 0; _i < 4; ++_i)
for (_j = 0; _j < 4; ++_j)
  {
  tmpmem[_i + _j] += (long long int)(a[_i]) * b[_j];
  if (tmpmem[_i + _j] > modconst)
    {
    tmpmem[_i + _j + 1] += tmpmem[_i + _j] >> LINTAR_DIG_CONST;
    tmpmem[_i + _j] &= andconst;
    }
  }
for (_i = 0; _i < 8; ++_i)
  r[_i] = tmpmem[_i];
}

// Recursive part of multiplication. Does not allocate or free memory.
// Here number format is different than in rest of code. Number lengths are not
// stored in x[0], but log2 of inputs length is stored in independent variable
// pow2.
// Instead, digits are stored starting with x[0], not x[1].
// sslab indicates beginning of allocated dynamic memory used for intermediate
// computations.
void lint_karatsuba_mult(int *r, int *a, int *b, int pow2, int *sslab)
{
int lslab, _i, rlen, dlen, modconst, andconst;
modconst = 1 << LINTAR_DIG_CONST;
andconst = modconst - 1;

lslab = nec_mem(pow2);
dlen = 1 << pow2;
rlen = 2 * dlen;

if (pow2 == 2)
  {                                         // recursion bottom
  straight4mult(r, a, b);
  return;
  }

                                            // evaluating intermediate products
lint_karatsuba_mult(sslab + 0 * lslab / 4,
                   a, b, pow2 - 1,
                   sslab + 0 * lslab / 4 + dlen);
lint_karatsuba_mult(sslab + 1 * lslab / 4,
                   a + dlen/2, b, pow2 - 1,
                   sslab + 1 * lslab / 4 + dlen);
lint_karatsuba_mult(sslab + 2 * lslab / 4,
                   a, b + dlen/2, pow2 - 1,
                   sslab + 2 * lslab / 4 + dlen);
lint_karatsuba_mult(sslab + 3 * lslab / 4,
                   a + dlen/2, b + dlen/2, pow2 - 1,
                   sslab + 3 * lslab / 4 + dlen);

for (_i = 0; _i < rlen; ++_i)               // initializing result
  r[_i] = 0;

for (_i = 0; _i < rlen; ++_i)
  {                                         // adding parts together
  if (_i < rlen / 2)
    r[_i] += sslab[_i];
                                            // left*right + right*left part
  if ((_i < 3 * rlen / 4) && (_i >= 1 * rlen /4))
    r[_i] += sslab[_i - rlen /4 + 1 * lslab / 4] + sslab[_i - rlen /4 + 2 * lslab / 4];
                                            // left*left part
  if (_i >= rlen / 2)
    r[_i] += sslab[_i - rlen /2 + 3 * lslab / 4];
  
  if (r[_i] > modconst)
    {
    r[_i + 1] += r[_i] >> LINTAR_DIG_CONST;
    r[_i] &= andconst;
    }
  }
}

char lint_fast_mult(int **r, int *a, int *b, int precision)
{
int maxl, pow2, *sslab, _i, dlen, sign, len;
maxl = (abs(a[0]) > abs(b[0]) ? abs(a[0]) : abs(b[0]));

pow2 = 2;                                   // finding the least pow2 such than maxl <= 2^pow2.
while ((1 << pow2) < maxl)
  ++pow2;
dlen = 1 << pow2;
                                            // allocating memory for storing formated inputs and outputs
                                            // and for multiplication
sslab = malloc(sizeof(int) * (4 * dlen + nec_mem(pow2)));
if (sslab == NULL)
  return YERROR_MEMORY;

for (_i = 0; _i < dlen; ++_i)
  {                                         // storing formated inputs in allocated memory
  if (_i + 1 <= abs(a[0]))
    sslab[2 * dlen + _i] = a[_i + 1];
  else
    sslab[2 * dlen + _i] = 0;
  if (_i + 1 <= abs(b[0]))
    sslab[3 * dlen + _i] = b[_i + 1];
  else
    sslab[3 * dlen + _i] = 0;
  }
                                             // recursive function call
lint_karatsuba_mult(sslab, sslab + 2 * dlen, sslab + 3 * dlen, pow2, sslab + 4 * dlen);
sign =  SIGN(a[0]) * SIGN(b[0]);
len = abs(a[0]) + abs(b[0]) - precision;
if (len < precision + 2)
  len = precision + 2;

*r = malloc(sizeof(int) * (len + 2));

(*r)[0] = 0;                                // getting formated outputs out of allocated memory
for (_i = 0; _i < 2 * dlen - precision; ++_i)
  {
  if(_i + 1 <= len)
    {
    (*r)[_i + 1] = sslab[_i + precision];
    if (((*r)[_i + 1] != 0) && ((*r)[0] < _i + 1))
      (*r)[0] = _i + 1;                     // updating result length
    }
  }
  
(*r)[0] *= sign;                            // applying result sign

free(sslab);
return YERROR_NO_ERRORS;
}

void lint_det2(int **r, int *a00, int *a01,
                        int *a10, int *a11, int precision)
{
int *t1, *t2;

lint_fast_mult(&t1, a00, a11, precision);
lint_fast_mult(&t2, a01, a10, precision);

lint_sub(r, t1, t2);

free(t1);
free(t2);
}

void lint_det3(int **r, int *a00, int *a01, int *a02,
                        int *a10, int *a11, int *a12,
                        int *a20, int *a21, int *a22, int precision)
{
int *t0, *t0m, *t1, *t1m, *t2, *t2m;

lint_det2(&t0, a10, a11, a20, a21, precision);
lint_fast_mult(&t0m, t0, a02, precision);
free(t0);

lint_det2(&t1, a00, a01, a20, a21, precision);
t1[0] *= -1;
lint_fast_mult(&t1m, t1, a12, precision);
free(t1);

lint_det2(&t2, a00, a01, a10, a11, precision);
lint_fast_mult(&t2m, t2, a22, precision);
free(t2);

lint_add(&t0, t0m, t1m);
lint_add(r, t0, t2m);
free(t0);
}

void lint_det4(int **r, int *a00, int *a01, int *a02, int *a03,
                        int *a10, int *a11, int *a12, int *a13,
                        int *a20, int *a21, int *a22, int *a23,
                        int *a30, int *a31, int *a32, int *a33, int precision)
{
int *t0, *t0m, *t1, *t1m, *t2, *t2m, *t3, *t3m;

lint_det3(&t0, a10, a11, a12, a20, a21, a22, a30, a31, a32, precision);
t0[0] *= -1;
lint_fast_mult(&t0m, t0, a03, precision);
free(t0);

lint_det3(&t1, a00, a01, a02, a20, a21, a22, a30, a31, a32, precision);
lint_fast_mult(&t1m, t1, a13, precision);
free(t1);

lint_det3(&t2, a00, a01, a02, a10, a11, a12, a30, a31, a32, precision);
t2[0] *= -1;
lint_fast_mult(&t2m, t2, a23, precision);
free(t2);

lint_det3(&t3, a00, a01, a02, a10, a11, a12, a20, a21, a22, precision);
lint_fast_mult(&t3m, t3, a33, precision);
free(t3);

lint_add(&t0, t0m, t1m);
lint_add(&t2, t2m, t3m);
lint_add(r, t0, t2);
free(t0);
free(t2);
}

void lint_det5(int **r, int *a00, int *a01, int *a02, int *a03, int *a04,
                        int *a10, int *a11, int *a12, int *a13, int *a14,
                        int *a20, int *a21, int *a22, int *a23, int *a24,
                        int *a30, int *a31, int *a32, int *a33, int *a34,
                        int *a40, int *a41, int *a42, int *a43, int *a44, int precision)
{
int *t0, *t0m, *t1, *t1m, *t2, *t2m, *t3, *t3m, *t4, *t4m;

lint_det4(&t0, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33, a40, a41, a42, a43, precision);
lint_fast_mult(&t0m, t0, a04, precision);
free(t0);

lint_det4(&t1, a00, a01, a02, a03, a20, a21, a22, a23, a30, a31, a32, a33, a40, a41, a42, a43, precision);
t1[0] *= -1;
lint_fast_mult(&t1m, t1, a14, precision);
free(t1);

lint_det4(&t2, a00, a01, a02, a03, a10, a11, a12, a13, a30, a31, a32, a33, a40, a41, a42, a43, precision);
lint_fast_mult(&t2m, t2, a24, precision);
free(t2);

lint_det4(&t3, a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a40, a41, a42, a43, precision);
t3[0] *= -1;
lint_fast_mult(&t3m, t3, a34, precision);
free(t3);

lint_det4(&t4, a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33, precision);
lint_fast_mult(&t4m, t4, a44, precision);
free(t4);

lint_add(&t0, t0m, t1m);
lint_add(&t1, t2m, t3m);
lint_add(&t2, t0, t1);
lint_add(r, t2, t4m);

free(t0);
free(t1);
free(t2);
}

// Increments v to identify next term of significance in SoS search.
// Returns 1 if increment is impossible (last term).
// Returns -1 if input is invalid
// Returns 0 otherwise.
char inc_v(int v_size, int *v)
{
int _i, _incs = 0;

if ((v[v_size - 1] != v_size - 1) && (v[v_size - 1] != v_size - 2))
  return -1;                                // invalid input

for (_i = 0; _i < v_size - 1; ++_i)
  {
  if (v[_i] < 0)
    return -1;
  if (v[_i] < v[_i + 1])
    ++_incs;                                // counting rows to delete
  else
    if (v[_i] > v[_i + 1])
      return -1;                            // invalid input
  }
if ((v[v_size - 1] == v_size - 2) && (_incs == v_size - 2))
  return 1;                                 // last term
else
  if (_incs == v_size - 1)
    return 1;                               // last term

_i = 0;                                     // searching for first non-zero
while (v[_i] == 0)
  ++_i;
--v[_i];  --_i;                             // decrementing component
while (_i >= 0)                             // updating zeros
  {
  v[_i] = v[_i + 1];
  --_i;
  }

return 0;
}

// Returns sign of significant term in SoS search, which is identified by vector v.
// m_dim_size denotes dimension of the matrix (2 for 2x2, 3 for 3x3 and so on).
// precision denotes format of lintars in matrix.
char term_sign(int m_dim_size, int **matrix, int *v, int precision)
{
int **tmp_matrix, _i, _j, *res, res_dim_size, probe, sign;
char *active_rows, *active_cols;
active_rows = malloc(sizeof(char) * m_dim_size);
active_cols = malloc(sizeof(char) * m_dim_size);

res_dim_size = m_dim_size;
sign = 1;
for (_i = 0; _i < m_dim_size; ++_i)         //scanning which rows and columns
  {                                         //to be crossed out
  active_rows[_i] = TRUE;
  active_cols[_i] = TRUE;
  }
for (_i = 0; _i < m_dim_size; ++_i)
  if (v[_i] < v[_i + 1])
    {
    active_rows[_i]    = FALSE;
    active_cols[v[_i]] = FALSE;
    --res_dim_size;
    sign *= (((_i + v[_i]) % 2 == 0) ? 1 : -1);
    }

if (res_dim_size <= 0)
  {
  free(active_rows);
  free(active_cols);
  return 1;                                 //trivial case
  }

if (res_dim_size > 5)                       //minor too large
  {
  free(active_rows);
  free(active_cols);
  error_exit("\nERROR! Such determinant size not implemented!\n");
  return 0;
  }

tmp_matrix = malloc(sizeof(int*) * res_dim_size * res_dim_size);
probe = 0;
for (_i = 0; _i < m_dim_size; ++_i)         //filling minor
for (_j = 0; _j < m_dim_size; ++_j)
  if (active_rows[_i] && active_cols[_j])
    tmp_matrix[probe++] = matrix[_i * m_dim_size + _j];

free(active_rows);
free(active_cols);

if (res_dim_size == 5)
  lint_det5(&res, tmp_matrix[0],  tmp_matrix[1],  tmp_matrix[2],  tmp_matrix[3],  tmp_matrix[4],
                  tmp_matrix[5],  tmp_matrix[6],  tmp_matrix[7],  tmp_matrix[8],  tmp_matrix[9],
                  tmp_matrix[10], tmp_matrix[11], tmp_matrix[12], tmp_matrix[13], tmp_matrix[14],
                  tmp_matrix[15], tmp_matrix[16], tmp_matrix[17], tmp_matrix[18], tmp_matrix[19],
                  tmp_matrix[20], tmp_matrix[21], tmp_matrix[22], tmp_matrix[23], tmp_matrix[24],
                  precision);

if (res_dim_size == 4)
  lint_det4(&res, tmp_matrix[0],  tmp_matrix[1],  tmp_matrix[2],  tmp_matrix[3],
                  tmp_matrix[4],  tmp_matrix[5],  tmp_matrix[6],  tmp_matrix[7],
                  tmp_matrix[8],  tmp_matrix[9],  tmp_matrix[10], tmp_matrix[11],
                  tmp_matrix[12], tmp_matrix[13], tmp_matrix[14], tmp_matrix[15],
                  precision);
if (res_dim_size == 3)
  lint_det3(&res, tmp_matrix[0],  tmp_matrix[1],  tmp_matrix[2],
                  tmp_matrix[3],  tmp_matrix[4],  tmp_matrix[5],
                  tmp_matrix[6],  tmp_matrix[7],  tmp_matrix[8],
                  precision);
if (res_dim_size == 2)
  lint_det2(&res, tmp_matrix[0],  tmp_matrix[1],
                  tmp_matrix[2],  tmp_matrix[3],
                  precision);
if (res_dim_size == 1)
  {
  res = malloc(sizeof(int));
  res[0] = tmp_matrix[0][0];
  }

sign *= SIGN(res[0]);
free(res);
free(tmp_matrix);
return sign;
}

char lint_SoS_udet2(double a00,
                    double a10, int precision)
{
char res = 0;
int *matrix[4], _i, v[3];
precision *= 2;                             //now we guarantee absolute presicion

double2lint(&(matrix[0]), a00, precision);  //init data
double2lint(&(matrix[1]), 1.,  precision);
double2lint(&(matrix[2]), a10, precision);
double2lint(&(matrix[3]), 1.,  precision);
v[0] = v[1] = v[2] = 1;

while (res == 0)                            //SoS search
  {
  res = term_sign(2, matrix, v, precision);
  inc_v(3, v);
  }

for (_i = 0; _i < 4; ++_i)
  free(matrix[_i]);

return res;
}

char lint_SoS_det2(double a00, double a01,
                   double a10, double a11, int precision)
{
char res = 0;
int *matrix[4], _i, v[3];
precision *= 2;                             //now we guarantee absolute presicion

double2lint(&(matrix[0]), a00, precision);  //init data
double2lint(&(matrix[1]), a01, precision);
double2lint(&(matrix[2]), a10, precision);
double2lint(&(matrix[3]), a11, precision);
v[0] = v[1] = v[2] = 2;

while (res == 0)                            //SoS search
  {
  res = term_sign(2, matrix, v, precision);
  inc_v(3, v);
  }

for (_i = 0; _i < 4; ++_i)
  free(matrix[_i]);

return res;
}

char lint_SoS_udet3(double a00, double a01,
                    double a10, double a11,
                    double a20, double a21, int precision)
{
char res = 0;
int *matrix[9], _i, v[4];
precision *= 3;                             //now we guarantee absolute presicion

double2lint(&(matrix[0]), a00, precision);  //init data
double2lint(&(matrix[1]), a01, precision);
double2lint(&(matrix[2]), 1.,  precision);
double2lint(&(matrix[3]), a10, precision);
double2lint(&(matrix[4]), a11, precision);
double2lint(&(matrix[5]), 1.,  precision);
double2lint(&(matrix[6]), a20, precision);
double2lint(&(matrix[7]), a21, precision);
double2lint(&(matrix[8]), 1.,  precision);
v[0] = v[1] = v[2] = v[3] = 2;

while (res == 0)                            //SoS search
  {
  res = term_sign(3, matrix, v, precision);
  inc_v(4, v);
  }

for (_i = 0; _i < 9; ++_i)
  free(matrix[_i]);

return res;
}

char lint_SoS_det3(double a00, double a01, double a02,
                   double a10, double a11, double a12, 
                   double a20, double a21, double a22, int precision)
{
char res = 0;
int *matrix[9], _i, v[4];
precision *= 3;                             //now we guarantee absolute presicion

double2lint(&(matrix[0]), a00, precision);  //init data
double2lint(&(matrix[1]), a01, precision);
double2lint(&(matrix[2]), a02, precision);
double2lint(&(matrix[3]), a10, precision);
double2lint(&(matrix[4]), a11, precision);
double2lint(&(matrix[5]), a12, precision);
double2lint(&(matrix[6]), a20, precision);
double2lint(&(matrix[7]), a21, precision);
double2lint(&(matrix[8]), a22, precision);
v[0] = v[1] = v[2] = v[3] = 3;

while (res == 0)                            //SoS search
  {
  res = term_sign(3, matrix, v, precision);
  inc_v(4, v);
  }

for (_i = 0; _i < 9; ++_i)
  free(matrix[_i]);

return res;
}

char lint_SoS_udet4(double a00, double a01, double a02,
                    double a10, double a11, double a12,
                    double a20, double a21, double a22,
                    double a30, double a31, double a32, int precision)
{
char res = 0;
int *matrix[16], _i, v[5];
precision *= 4;                             //now we guarantee absolute presicion

double2lint(&(matrix[0]),  a00, precision); //init data
double2lint(&(matrix[1]),  a01, precision);
double2lint(&(matrix[2]),  a02, precision);
double2lint(&(matrix[3]),  1.,  precision);
double2lint(&(matrix[4]),  a10, precision);
double2lint(&(matrix[5]),  a11, precision);
double2lint(&(matrix[6]),  a12, precision);
double2lint(&(matrix[7]),  1.,  precision);
double2lint(&(matrix[8]),  a20, precision);
double2lint(&(matrix[9]),  a21, precision);
double2lint(&(matrix[10]), a22, precision);
double2lint(&(matrix[11]), 1.,  precision);
double2lint(&(matrix[12]), a30, precision);
double2lint(&(matrix[13]), a31, precision);
double2lint(&(matrix[14]), a32, precision);
double2lint(&(matrix[15]), 1.,  precision);
v[0] = v[1] = v[2] = v[3] = v[4] = 3;

while (res == 0)                            //SoS search
  {
  res = term_sign(4, matrix, v, precision);
  inc_v(5, v);
  }

for (_i = 0; _i < 16; ++_i)
  free(matrix[_i]);

return res;
}

char lint_SoS_det4(double a00, double a01, double a02, double a03,
                   double a10, double a11, double a12, double a13,
                   double a20, double a21, double a22, double a23,
                   double a30, double a31, double a32, double a33, int precision)
{
char res = 0;
int *matrix[16], _i, v[5];
precision *= 4;                             //now we guarantee absolute presicion

double2lint(&(matrix[0]),  a00, precision); //init data
double2lint(&(matrix[1]),  a01, precision);
double2lint(&(matrix[2]),  a02, precision);
double2lint(&(matrix[3]),  a03, precision);
double2lint(&(matrix[4]),  a10, precision);
double2lint(&(matrix[5]),  a11, precision);
double2lint(&(matrix[6]),  a12, precision);
double2lint(&(matrix[7]),  a13, precision);
double2lint(&(matrix[8]),  a20, precision);
double2lint(&(matrix[9]),  a21, precision);
double2lint(&(matrix[10]), a22, precision);
double2lint(&(matrix[11]), a23, precision);
double2lint(&(matrix[12]), a30, precision);
double2lint(&(matrix[13]), a31, precision);
double2lint(&(matrix[14]), a32, precision);
double2lint(&(matrix[15]), a33, precision);
v[0] = v[1] = v[2] = v[3] = v[4] = 4;

while (res == 0)                            //SoS search
  {
  res = term_sign(4, matrix, v, precision);
  inc_v(5, v);
  }

for (_i = 0; _i < 16; ++_i)
  free(matrix[_i]);

return res;
}

char lint_SoS_udet5(double a00, double a01, double a02, double a03,
                    double a10, double a11, double a12, double a13,
                    double a20, double a21, double a22, double a23,
                    double a30, double a31, double a32, double a33,
                    double a40, double a41, double a42, double a43, int precision)
{
char res = 0;
int *matrix[25], _i, v[6];
precision *= 5;                             //now we guarantee absolute presicion

double2lint(&(matrix[0]),  a00, precision); //init data
double2lint(&(matrix[1]),  a01, precision);
double2lint(&(matrix[2]),  a02, precision);
double2lint(&(matrix[3]),  a03, precision);
double2lint(&(matrix[4]),  1.,  precision);
double2lint(&(matrix[5]),  a10, precision);
double2lint(&(matrix[6]),  a11, precision);
double2lint(&(matrix[7]),  a12, precision);
double2lint(&(matrix[8]),  a13, precision);
double2lint(&(matrix[9]),  1.,  precision);
double2lint(&(matrix[10]), a20, precision);
double2lint(&(matrix[11]), a21, precision);
double2lint(&(matrix[12]), a22, precision);
double2lint(&(matrix[13]), a23, precision);
double2lint(&(matrix[14]), 1.,  precision);
double2lint(&(matrix[15]), a30, precision);
double2lint(&(matrix[16]), a31, precision);
double2lint(&(matrix[17]), a32, precision);
double2lint(&(matrix[18]), a33, precision);
double2lint(&(matrix[19]), 1.,  precision);
double2lint(&(matrix[20]), a40, precision);
double2lint(&(matrix[21]), a41, precision);
double2lint(&(matrix[22]), a42, precision);
double2lint(&(matrix[23]), a43, precision);
double2lint(&(matrix[24]), 1.,  precision);
v[0] = v[1] = v[2] = v[3] = v[4] = v[5] = 4;

while (res == 0)                            //SoS search
  {
  res = term_sign(5, matrix, v, precision);
  inc_v(6, v);
  }

for (_i = 0; _i < 25; ++_i)
  free(matrix[_i]);

return res;
}

char lint_SoS_det5(double a00, double a01, double a02, double a03, double a04,
                   double a10, double a11, double a12, double a13, double a14,
                   double a20, double a21, double a22, double a23, double a24,
                   double a30, double a31, double a32, double a33, double a34,
                   double a40, double a41, double a42, double a43, double a44, int precision)
{
char res = 0;
int *matrix[25], _i, v[6];
precision *= 5;                             //now we guarantee absolute presicion

double2lint(&(matrix[0]),  a00, precision); //init data
double2lint(&(matrix[1]),  a01, precision);
double2lint(&(matrix[2]),  a02, precision);
double2lint(&(matrix[3]),  a03, precision);
double2lint(&(matrix[4]),  a04, precision);
double2lint(&(matrix[5]),  a10, precision);
double2lint(&(matrix[6]),  a11, precision);
double2lint(&(matrix[7]),  a12, precision);
double2lint(&(matrix[8]),  a13, precision);
double2lint(&(matrix[9]),  a14, precision);
double2lint(&(matrix[10]), a20, precision);
double2lint(&(matrix[11]), a21, precision);
double2lint(&(matrix[12]), a22, precision);
double2lint(&(matrix[13]), a23, precision);
double2lint(&(matrix[14]), a24, precision);
double2lint(&(matrix[15]), a30, precision);
double2lint(&(matrix[16]), a31, precision);
double2lint(&(matrix[17]), a32, precision);
double2lint(&(matrix[18]), a33, precision);
double2lint(&(matrix[19]), a34, precision);
double2lint(&(matrix[20]), a40, precision);
double2lint(&(matrix[21]), a41, precision);
double2lint(&(matrix[22]), a42, precision);
double2lint(&(matrix[23]), a43, precision);
double2lint(&(matrix[24]), a44, precision);
v[0] = v[1] = v[2] = v[3] = v[4] = v[5] = 5;

while (res == 0)                            //SoS search
  {
  res = term_sign(5, matrix, v, precision);
  inc_v(6, v);
  }

for (_i = 0; _i < 25; ++_i)
  free(matrix[_i]);

return res;
}
