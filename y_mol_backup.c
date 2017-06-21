//This file contain universal molecule distribution in y_system and connected structure manipulation routines
#include "y_mol.h"

extern unsigned int ylib_errno;

#define N_ATOMS 63
#define PDB_STR_LEN            0x38
#define TOP_STR_LEN            0xFE
#define MOL2_STR_LEN           0x50


#define UNSET_FLAG 9999.9999


//This script permit to atchive end of the buffer
#define GET_FILE_STRING_END(buffer,in) { if (!(check_string_end(buffer))) fstrskip(in); }

//********************************     T O P O L O G Y     P A R T     ******************************************/

//This function free top structure
void  free_top(t_top *top)
{
if (top->ff_a) free(top->ff_a);
if (top->A) free(top->A);
if (top->B) free(top->B);
if (top->ff_n) free(top->ff_n);
if (top->ff_b) free(top->ff_b);
if (top->ff_g) free(top->ff_g);
if (top->ff_i) free(top->ff_i);
if (top->ff_d) free(top->ff_d);
if (top->_atoms) free(top->_atoms);
if (top->_ctypes) free(top->_ctypes);
if (top->_edges) free(top->_edges);
if (top->ress)   free(top->ress);
if (top->res)    free(top->res);
if (top) free(top);
}

//This functions reorder members of the topological units so they gives the lowest possible value in compare_xxxx functions
//This function order atom types in the topology primitives 
inline void order_bond(unsigned int *x0,unsigned int *x1,unsigned int a0,unsigned int a1)
{
register unsigned int _t;
if (a0>a1) { _t=*x0, *x0=*x1, *x1=_t; }
}
inline void order_angle(unsigned int *x0,unsigned int *x1,unsigned int *x2,unsigned int a0,unsigned int a1,unsigned int a2)
{
register unsigned int _t;
if (a0>a2) { _t=*x0, *x0=*x2, *x2=_t; }
}
inline void order_impr(unsigned int type,unsigned int *x0,unsigned int *x1,unsigned int *x2,unsigned int *x3,unsigned int a0,unsigned int a1,unsigned int a2,unsigned int a3)
{
register unsigned int _t;
if (!(type))
  {// (A0)->A1-A2-A3
  if (a1>a2)
    {
    if (a1>a3)
      {
      _t=*x1, *x1=*x3, *x3=_t;  //a1 is the biggest
      if (a2>a3)                          ;   //a3<a2<a1 
      else       { _t=*x1, *x1=*x2, *x2=_t; } //a2<a3<a1
      }
    else         { _t=*x1, *x1=*x2, *x2=_t; } //a2<a1<a3
    }
  else
    {
    if (a2>a3) 
      {
      _t=*x2, *x2=*x3, *x3=_t; //a2 is the biggest
      if (a1>a3) { _t=*x1, *x1=*x2, *x2=_t; }  //a3<a1<a2 
      else                                ;    //a1<a3<a2
      }
    else                                  ;    //a1<a2<a3
    }
  }
else
  {//A0<-A1-A2->A3
       if (a1>a2)               { _t=*x0, *x0=*x3, *x3=_t; _t=*x1, *x1=*x2, *x2=_t;  }
  else if ( (a1==a2)&&(a0>a3) ) { _t=*x0, *x0=*x3, *x3=_t; _t=*x1, *x1=*x2, *x2=_t;  }
  }
}
inline void order_dih(unsigned int *x0,unsigned int *x1,unsigned int *x2,unsigned int *x3,unsigned int a0,unsigned int a1,unsigned int a2,unsigned int a3)
{
order_impr(TRUE,x0,x1,x2,x3,a0,a1,a2,a3);
}

//This function compare two angles topolies to decode whos is lower(greater)
int compare_bonds(const void *b1,const void *b2)
{
     if (((t_ff_b*)b1)->atom[0]==((t_ff_b*)b2)->atom[0])
       {//The first atom do not solve comparison problem
            if (((t_ff_b*)b1)->atom[1]==((t_ff_b*)b2)->atom[1]) return  0; //The second atom solve comparison problem anyway
       else if (((t_ff_b*)b1)->atom[1]> ((t_ff_b*)b2)->atom[1]) return +1;
       else                                                     return -1;
       }
else if (((t_ff_b*)b1)->atom[0]>((t_ff_b*)b2)->atom[0])         return +1;
else                                                            return -1;
}
//This function compare two angles topolies to decode whos is lower(greater)
int compare_angles(const void *a1,const void *a2)
{
     if (((t_ff_g*)a1)->atom[1]==((t_ff_g*)a2)->atom[1])
       {//The second atom do not solve comparison problem
            if (((t_ff_g*)a1)->atom[0]==((t_ff_g*)a2)->atom[0])
              {//The first atom do not solve comparison problem
                   if (((t_ff_g*)a1)->atom[2]==((t_ff_g*)a2)->atom[2]) return  0; //The third atom solve comparison problem anyway
              else if (((t_ff_g*)a1)->atom[2]> ((t_ff_g*)a2)->atom[2]) return +1;
              else                                                     return -1;
              }
       else if (((t_ff_g*)a1)->atom[0]>((t_ff_g*)a2)->atom[0])         return +1;
       else                                                            return -1;
       }
else if (((t_ff_g*)a1)->atom[1]>((t_ff_g*)a2)->atom[1])                return +1;
else                                                                   return -1;
}

//This function compare two impropers topolies to decode who is lower(greater)
int compare_imprs(const void *i1,const void *i2)
{
     if (((t_ff_i*)i1)->atom[0]==((t_ff_i*)i2)->atom[0])
       { //First item do not solve comparison problem
            if (((t_ff_i*)i1)->atom[1]==((t_ff_i*)i2)->atom[1])
              { //Second item do not solve comparison problem
                   if (((t_ff_i*)i1)->atom[2]==((t_ff_i*)i2)->atom[2])
                     { //Third item do not solve comparison problem
                          if (((t_ff_i*)i1)->atom[3]==((t_ff_i*)i2)->atom[3]) return  0;
                     else if (((t_ff_i*)i1)->atom[3]>((t_ff_i*)i2)->atom[3])  return +1;
                     else                                                     return -1;
                     }
              else if (((t_ff_i*)i1)->atom[2]>((t_ff_i*)i2)->atom[2])         return +1;
              else                                                            return -1;
              } 
       else if (((t_ff_i*)i1)->atom[1]>((t_ff_i*)i2)->atom[1])                return +1;
       else                                                                   return -1;
       }
else if (((t_ff_i*)i1)->atom[0]>((t_ff_i*)i2)->atom[0])                       return +1;
else                                                                          return -1;
}

//This function compare two torsions topolies to decode who is lower(greater); topology version using the above function
int compare_dihs(const void *d1,const void *d2)
{
     if (((t_ff_d*)d1)->atom[1]==((t_ff_d*)d2)->atom[1])
       { //First item do not solve comparison problem
            if (((t_ff_d*)d1)->atom[2]==((t_ff_d*)d2)->atom[2])
              { //Second item do not solve comparison problem
                   if (((t_ff_d*)d1)->atom[0]==((t_ff_d*)d2)->atom[0])
                     { //Third item do not solve comparison problem
                          if (((t_ff_d*)d1)->atom[3]==((t_ff_d*)d2)->atom[3]) return  0;
                     else if (((t_ff_d*)d1)->atom[3]>((t_ff_d*)d2)->atom[3])  return +1;
                     else                                                     return -1;
                     }
              else if (((t_ff_d*)d1)->atom[0]>((t_ff_d*)d2)->atom[0])         return +1;
              else                                                            return -1;
              }
       else if (((t_ff_d*)d1)->atom[2]>((t_ff_d*)d2)->atom[2])                return +1;
       else                                                                   return -1;
       }
else if (((t_ff_d*)d1)->atom[1]>((t_ff_d*)d2)->atom[1])                       return +1;
else                                                                          return -1;
}
//This function compare two torsions topolies to decode who is lower(greater)
int compare_torss(const void *t1,const void *t2)
{
     if (((t_ff_t*)t1)->atom[1]==((t_ff_t*)t2)->atom[1])
       { //First item do not solve comparison problem
            if (((t_ff_t*)t1)->atom[2]==((t_ff_t*)t2)->atom[2])
              { //Second item do not solve comparison problem
                   if (((t_ff_t*)t1)->atom[0]==((t_ff_t*)t2)->atom[0])
                     { //Third item do not solve comparison problem
                          if (((t_ff_t*)t1)->atom[3]==((t_ff_t*)t2)->atom[3]) return  0; //Forth solve problem anyway
                     else if (((t_ff_t*)t1)->atom[3]> ((t_ff_t*)t2)->atom[3]) return +1;
                     else                                                     return -1;
                     }
              else if (((t_ff_t*)t1)->atom[0]>((t_ff_t*)t2)->atom[0])         return +1;
              else                                                            return -1;
              }
       else if (((t_ff_t*)t1)->atom[2]>((t_ff_t*)t2)->atom[2])                return +1;
       else                                                                   return -1;
       }
else if (((t_ff_t*)t1)->atom[1]>((t_ff_t*)t2)->atom[1])                       return +1;
else                                                                          return -1;
}

//This function upload topoly from files
t_top *read_top(FILE *in_atoms_yff1,FILE *in_pairs_yff1,FILE *in_atoms_gaff,FILE *in_bonds_gaff,FILE *in_angles_gaff,FILE *in_imprs_gaff,FILE *in_dihs_gaff,FILE *in_ress)
{	
t_top *top;
extern char buffer[];
char *lexem, _c, *file_type;
register unsigned int _i, _j, count;
register void *vp, *_vp;

double d;
unsigned int id, ff_n;
t_ff_b ff_b;
t_ff_g ff_g;
t_ff_i ff_i;
t_ff_d ff_d;


//Stage 0. Allocate memory for topology
if (!(top=calloc(sizeof(t_top),0x1))) 
  {
  lexem="top structure";
  ERROR_MEMORY_EXIT:
  ylib_errno=YERROR_MEMORY; 
  return FALSE;
  }

//Stage 1. Upload ff atoms with 0xFF buffer
top->size_a=count=0, file_type="yff1 atom types";
if (!(top->ff_a=(t_ff_a*)malloc(sizeof(t_ff_a)*0xFF))) 
  {
  lexem="ff_a";
  LABEL_ERROR_MEMORY: 
  ylib_errno=YERROR_MEMORY;
  free_top(top); top=0x0;
  goto ERROR_MEMORY_EXIT;
  }
while (fgets(buffer,TOP_STR_LEN,in_atoms_yff1))
  {//Check format
  count++; 
  CUT_COMMENTS(_i,buffer);
  if (!(check_lexem(0x1,buffer))) continue;
  if ( (!(check_lexem(0x12,buffer)))||( (check_lexem(0x13,buffer))) ) 
    { 
    lexem="wrong line structure"; 
    LABEL_ERROR_DATA_FORMAT: 
    yprintf(YPRINTF_ERROR,"Wrong format of data in %s file, line %d (%s)\n",file_type, count, lexem);
    ylib_errno=YERROR_DATA_FORMAT;
    free_top(top); top=0x0;
    return FALSE;
    }
  //Check serial number
  if ( (!(lexem=get_lex(0x1,buffer)))||(!(check_int_type(lexem)))||((unsigned int)atoi(lexem)!=top->size_a) ) 
    { lexem="wrong id"; goto LABEL_ERROR_DATA_FORMAT; }
  //Upload chem id
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||(abs((int)(_i=atoi(lexem)))>0xFE) ) 
    { lexem="wrong chem_id"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].chem_id=(char)((int)_i);
  //Upload atom gaff type
  if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) goto LABEL_ERROR_DATA_FORMAT;
  else { 
       _j=0; do ((char*)&top->ff_a[top->size_a].gtype)[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
       while (_j<sizeof(unsigned int)) ((char*)&top->ff_a[top->size_a].gtype)[_j++]=' ';
       if (((char*)&top->ff_a[top->size_a].gtype)[3]!=' ')
         { lexem="The 4-th character of gtype is reserved for compilator's of chemical envinronment suffixes"; goto LABEL_ERROR_DATA_FORMAT; }
       }
  //Upload atom yff1 generic type
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem))>sizeof(unsigned int)) ) goto LABEL_ERROR_DATA_FORMAT;
  else top->ff_a[top->size_a].generic=atoi(lexem);
  //Upload atom mass
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem)))||((top->ff_a[top->size_a].mass=atof(lexem))<0.) )
    { lexem="wrong mass"; goto LABEL_ERROR_DATA_FORMAT; }
  //Upload topology parameters
  //Upload pi electrons number
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||((int)(top->ff_a[top->size_a].ppn=atoi(lexem))<0) )
    { lexem="wrong ppn";  goto LABEL_ERROR_DATA_FORMAT; }
  //Upload number of nonsigma edges
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem))) ) 
    { lexem="wrong pi-edges number"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].pnn=atoi(lexem);
  //Upload number of valent electrons
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||((int)(_i=atoi(lexem))>0x13)||((int)_i<0) ) 
    { lexem="wrong amount of valent electrons"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].vn=(char)_i;
  //Upload number of coordinations
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||((int)(_i=atoi(lexem))>0x0B)||((int)_i<0) )
    { lexem="wrong neighbors number"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].nn=(char)_i;
  //Upload involvement in resonance
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem))) ) 
    { lexem="wrong resonance number"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].rn=(char)(atoi(lexem));
  //Upload cycle size
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||(abs((int)(_i=atoi(lexem)))>0x06)||( (_i)&&(abs((int)_i)<0x03) ) )
    { lexem="wrong cycle number"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].cycle=(char)_i;
  //Upload atoms' hydrogen bonding donors amount (+1 for hydrogens up to -3 on fluorine)
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||((int)abs(_i=atoi(lexem))>0x03) )
    { lexem="wrong donors number"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].nD=(char)_i;
  //Upload nonnbonded FF parameters
  //Upload vdw A (+12. power constant)
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong A nonbonded constant"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].A=atof(lexem)*1.e+12;
  //Upload vdw A (+6.  power constant)
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong B nonbonded constant"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].B=atof(lexem)*1.e+6;
  //Upload Oliferenko's electronegativity
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong electronegativity"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].Eng=atof(lexem);
  //Upload Oliferenko's hardness
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong hardness"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].Hrd=atof(lexem);
  //Upload atoms meta-data
  //Upload vdW radii
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong vdW radii"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].rvdw=atof(lexem);
  //Upload Born gamma parameter for explicite solvent models
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) )
    { lexem="wrong Born gamma"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->ff_a[top->size_a].sgamma=atof(lexem);
  //manage the memory
  if (!((++top->size_a+0x1)%0xFF))
    {
    if (!(vp=realloc(top->ff_a,sizeof(t_ff_a)*(top->size_a+0x1+0xFF)))) 
      { lexem="yff1 types"; goto LABEL_ERROR_MEMORY; }
    else top->ff_a=(t_ff_a*)vp;
    }
  }
//Check generics
_i=top->size_a;
while (_i--)
  if ( ((int)top->ff_a[_i].generic<0)||(top->ff_a[_i].generic>=top->size_a) )
    {
    yprintf(YPRINTF_ERROR,"Wrong generic if of atom %1d\n",_i);
    ylib_errno=YERROR_DATA_FORMAT;
    free_top(top); top=0x0;
    return FALSE; 
    }
//Sync names
if (!(vp=realloc(top->ff_a,sizeof(t_ff_a)*top->size_a))) 
  { lexem="yff1 atoms"; goto LABEL_ERROR_MEMORY; }
else top->ff_a=(t_ff_a*)vp;

//Alloc memory for vdw parameters of yff1
if ( (!(top->A=(double**)malloc((sizeof(double*)+sizeof(double)*top->size_a)*top->size_a)))||
     (!(top->B=(double**)malloc((sizeof(double*)+sizeof(double)*top->size_a)*top->size_a))) )
  { lexem="A and B nonbonded parameters"; goto LABEL_ERROR_MEMORY; }
top->A[0x0]=(void*)top->A+sizeof(double*)*top->size_a;
top->B[0x0]=(void*)top->B+sizeof(double*)*top->size_a;
for (_i=0x1;_i<top->size_a;_i++)
  {
  top->A[_i]=(double*)((void*)top->A[_i-1]+sizeof(double)*top->size_a);
  top->B[_i]=(double*)((void*)top->B[_i-1]+sizeof(double)*top->size_a);
  }
//Fill zero row especially
_i=top->size_a;
while (--_i)
  {
  top->A[_i][0]=top->A[0][_i]=sqrt(top->ff_a[_i].A);
  top->B[_i][0]=top->B[0][_i]=sqrt(top->ff_a[_i].B);
  }
top->A[0][0]=top->B[0][0]=0.; //init zero-atom
//Fill the matrix as a prodct
_i=top->size_a;
while (--_i)
  {
  top->A[_i][_i]=top->ff_a[_i].A;
  top->B[_i][_i]=top->ff_a[_i].B;
  _j=_i;
  while (--_j)
    {
    top->A[_i][_j]=top->A[_j][_i]=top->A[0][_i]*top->A[0][_j];
    top->B[_i][_j]=top->B[_j][_i]=top->B[0][_i]*top->B[0][_j];
    }
  }

//Stage 2. Upload YFF corrections
top->nnb_corrections=count=0, file_type="yff1 nonbonded corrections";
while (fgets(buffer,TOP_STR_LEN,in_pairs_yff1))
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(check_lexem(0x1,buffer))) continue;
  if ( (!(check_lexem(0x4,buffer)))||( (check_lexem(0x5,buffer))) ) 
    { lexem="wrong line format"; goto LABEL_ERROR_DATA_FORMAT; }
  //Check type of correction
  if ( (!(lexem=get_lex(0x1,buffer)))||(lexlen(lexem)!=1) ) 
    { lexem="wrong type"; goto LABEL_ERROR_DATA_FORMAT; }
  else _c=*lexem;
  //Upload atom id #0
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||((_i=(int)atoi(lexem))<=0)||(_i>=top->size_a) )
    { lexem="wrong atom#0 type"; goto LABEL_ERROR_DATA_FORMAT; }
  //Upload atom id #1
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||((_j=(int)atoi(lexem))<=0)||(_j>=top->size_a) ) 
    { lexem="wrong atom#1 type"; goto LABEL_ERROR_DATA_FORMAT; }
  //Upload atom's pair interactions correction
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong correcion"; goto LABEL_ERROR_DATA_FORMAT; }
  //Apply the correction
       if ( (_c=='B')||(_c=='b') ) top->A[_i][_j]=top->A[_j][_i]=atof(lexem);
  else if ( (_c=='A')||(_c=='a') ) top->B[_i][_j]=top->B[_j][_i]=atof(lexem); 
  else { lexem="wrong type"; goto LABEL_ERROR_DATA_FORMAT; }
  top->nnb_corrections++;
  }
//Add safe-repulsive corrections
_i=top->size_a;
while (--_i)
  {
  if ((d=64.*cubed(sqrd(top->ff_a[_i].rvdw)))<2.) d=2;
  if (top->A[_i][_i]<d*top->B[_i][_i]) 
    { top->A[_i][_i]=d*top->B[_i][_i]; }                  //Apply a repulsive correction: A/r^6-B = 0 @ 2*rvdw_i, A <- r^6*B or A-B=B @ 1A
  _j=_i;
  while (--_j)
    {
    if ((d=cubed(sqrd(top->ff_a[_i].rvdw+top->ff_a[_j].rvdw)))<2.) d=2.;
    if (top->A[_i][_j]<d*top->B[_i][_j])
      { top->A[_i][_j]=top->A[_j][_i]=d*top->B[_i][_j]; } //Apply a repulsive correction: A/r^6-B = 0 @ rvdw_i+rvdw_j, A <- r^6*B or A-B=B @ 1A
    }
  }

//Stage 3. Upload gaff names space with buffer of 0xFF
if (!(top->ff_n=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) 
  { lexem="gaff types"; goto LABEL_ERROR_MEMORY; }
top->size_n=count=0, file_type="gaff atom types";
while (fgets(buffer,TOP_STR_LEN,in_atoms_gaff))
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(check_lexem(0x1,buffer))) continue;
  if ( (!(check_lexem(0x3,buffer)))||( (check_lexem(0x4,buffer))) ) 
    { lexem="wrong line format"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (!(lexem=get_lex(0x1,buffer)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wrong type name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {
       _j=0; do ((char*)&ff_n)[_j]=tolower((int)lexem[_j]); while (++_j!=_i); 
       while (_j<sizeof(unsigned int)) ((char*)&ff_n)[_j++]=' ';
       } 
  if ( (!(lexem=get_lex(0x2,lexem))) ||(!(check_double_type(lexem))) ) 
    { lexem="wrong mass"; goto LABEL_ERROR_DATA_FORMAT; }
  //else mass=atof(lexem);
  if ( (!(lexem=get_lex(0x2,lexem))) ||(!(check_double_type(lexem))) ) 
    { lexem="wrong electronegativity"; goto LABEL_ERROR_DATA_FORMAT; }
  //else electronegativity=atof(lexem);
  //keep the massive sorted for rapid search for dublicates
  if ( (find_in_sorted_lth_urow(&id,ff_n,top->size_n,top->ff_n)))
    { lexem="dublicated gaff name"; goto LABEL_ERROR_DATA_FORMAT; }
  else
    {
    _j=top->size_n;
    while (id!=_j--)
      {
      top->ff_n[_j+1]=top->ff_n[_j];
      }
    top->ff_n[id]=ff_n;
    }
  //manage the memory
  if (!((++top->size_n+0x1)%0xFF))
    {
    if (!(vp=realloc(top->ff_n,sizeof(unsigned int)*(top->size_n+0x1+0xFF)))) 
      { lexem="gaff types"; goto LABEL_ERROR_MEMORY; }
    else top->ff_n=(unsigned int*)vp;
    }
  }
//Sync names
if (!(vp=realloc(top->ff_n,sizeof(unsigned int)*top->size_n))) 
  { lexem="gaff atoms"; goto LABEL_ERROR_MEMORY; }
else top->ff_n=(unsigned int*)vp;

//Stage 4. Upload ff bonds with buffer of 0xFF
top->size_b=count=0, file_type="gaff bonds";
if (!(top->ff_b=(t_ff_b*)malloc(sizeof(t_ff_b)*0xFF))) 
  { lexem="gaff bonds"; goto LABEL_ERROR_MEMORY; }
while (fgets(buffer,TOP_STR_LEN,in_bonds_gaff))
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(check_lexem(0x1,buffer))) continue;
  if ( (!(check_lexem(0x4,buffer)))||( (check_lexem(0x5,buffer))) ) 
    { lexem="wrong line format"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (!(lexem=get_lex(0x1,buffer)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wrong atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {
       _j=0; do ((char*)&ff_b.atom[0])[_j]=tolower((int)lexem[_j]); while (++_j!=_i); 
       while (_j<sizeof(unsigned int)) ((char*)&ff_b.atom[0])[_j++]=' ';
       } 
  if ( (!(lexem=get_lex(0x2,lexem))) ||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wroing atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {
       _j=0; do ((char*)&ff_b.atom[1])[_j]=tolower((int)lexem[_j]); while (++_j!=_i); 
       while (_j<sizeof(unsigned int)) ((char*)&ff_b.atom[1])[_j++]=' ';
       } 
  if ( (!(lexem=get_lex(0x2,lexem))) ||(!(check_double_type(lexem))) ) 
    { lexem="wrong K"; goto LABEL_ERROR_DATA_FORMAT; }
  else ff_b.k=atof(lexem);
  if ( (!(lexem=get_lex(0x2,lexem))) ||(!(check_double_type(lexem))) )
    { lexem="wrong v0"; goto LABEL_ERROR_DATA_FORMAT; }
  else ff_b.v=atof(lexem);
  //check atom types
  if ( (ff_b.atom[0]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_b.atom[0],top->size_n,top->ff_n))) )
    { lexem="unknown atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (ff_b.atom[1]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_b.atom[1],top->size_n,top->ff_n))) ) 
    { lexem="unknown atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
  //reorder bond
  order_bond(&ff_b.atom[0],&ff_b.atom[1],ff_b.atom[0],ff_b.atom[1]);
  //keep the massive sorted for rapid search for dublicates
  if ( (find_in_sorted_lth_objects(&id,&ff_b,top->size_b,(void*)top->ff_b,sizeof(t_ff_b),compare_bonds)))
    { lexem="dublicated gaff bond"; goto LABEL_ERROR_DATA_FORMAT; }
  else
    {
    _j=top->size_b;
    while (id!=_j--)
      {
      top->ff_b[_j+1].atom[0]=top->ff_b[_j].atom[0],
      top->ff_b[_j+1].atom[1]=top->ff_b[_j].atom[1],
      top->ff_b[_j+1].k=top->ff_b[_j].k,
      top->ff_b[_j+1].v=top->ff_b[_j].v;
      }
    top->ff_b[id].atom[0]=ff_b.atom[0],
    top->ff_b[id].atom[1]=ff_b.atom[1],
    top->ff_b[id].k=ff_b.k,
    top->ff_b[id].v=ff_b.v;
    }
  //manage the memory
  if (!((++top->size_b+0x1)%0xFF))
    {
    if (!(vp=realloc(top->ff_b,sizeof(t_ff_b)*(top->size_b+0x1+0xFF)))) 
      { lexem="bonds"; goto LABEL_ERROR_MEMORY; }
    else top->ff_b=(t_ff_b*)vp;
    }
  }
//Sync bonds
if (!(vp=realloc(top->ff_b,sizeof(t_ff_b)*(top->size_b)))) 
  { lexem="bonds"; goto LABEL_ERROR_MEMORY; }
else top->ff_b=(t_ff_b*)vp;

//Stage 5. Upload ff angles with buffer of 0xFF
top->size_g=count=0, file_type="gaff angles";
if (!(top->ff_g=(t_ff_g*)malloc(sizeof(t_ff_g)*0xFF))) 
  { lexem="angles"; goto LABEL_ERROR_MEMORY; }
while (fgets(buffer,TOP_STR_LEN,in_angles_gaff))
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(check_lexem(0x1,buffer))) continue;
  if ( (!(check_lexem(0x5,buffer)))||( (check_lexem(0x6,buffer))) ) 
    { lexem="wrong line format"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (!(lexem=get_lex(0x1,buffer)))||((_i=lexlen(lexem))>sizeof(unsigned int)) )
    { lexem="wrong atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {
       _j=0; do ((char*)&ff_g.atom[0])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
       while (_j<sizeof(unsigned int)) ((char*)&ff_g.atom[0])[_j++]=' ';
       } 
  if ( (!(lexem=get_lex(0x2,lexem))) ||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wrong atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {
       _j=0; do ((char*)&ff_g.atom[1])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
       while (_j<sizeof(unsigned int)) ((char*)&ff_g.atom[1])[_j++]=' ';
       } 
  if ( (!(lexem=get_lex(0x2,lexem))) ||((_i=lexlen(lexem))>sizeof(unsigned int)) )
    { lexem="wrong atom#2 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {
       _j=0; do ((char*)&ff_g.atom[2])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
       while (_j<sizeof(unsigned int)) ((char*)&ff_g.atom[2])[_j++]=' ';
       }
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong K"; goto LABEL_ERROR_DATA_FORMAT; }
  else ff_g.k=atof(lexem);
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong v0"; goto LABEL_ERROR_DATA_FORMAT; }
  else ff_g.v=atof(lexem);
  //check atom types
  if ( (ff_g.atom[0]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_g.atom[0],top->size_n,top->ff_n))) )
    { lexem="unknown atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (ff_g.atom[1]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_g.atom[1],top->size_n,top->ff_n))) ) 
    { lexem="unknown atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (ff_g.atom[2]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_g.atom[2],top->size_n,top->ff_n))) ) 
    { lexem="unknown atom#2 name"; goto LABEL_ERROR_DATA_FORMAT; }
  //reorder angle
  order_angle(&ff_g.atom[0],&ff_g.atom[1],&ff_g.atom[2],ff_g.atom[0],ff_g.atom[1],ff_g.atom[2]);
  //keep the massive sorted for rapid search for dublicates
  if ( (find_in_sorted_lth_objects(&id,&ff_g,top->size_g,(void*)top->ff_g,sizeof(t_ff_g),compare_angles)))
    { lexem="dublicated gaff angles"; goto LABEL_ERROR_DATA_FORMAT; }
  else
    {
    _j=top->size_g;
    while (id!=_j--)
      {
      top->ff_g[_j+1].atom[0]=top->ff_g[_j].atom[0],
      top->ff_g[_j+1].atom[1]=top->ff_g[_j].atom[1],
      top->ff_g[_j+1].atom[2]=top->ff_g[_j].atom[2],
      top->ff_g[_j+1].k=top->ff_g[_j].k,
      top->ff_g[_j+1].v=top->ff_g[_j].v;
      }
    top->ff_g[id].atom[0]=ff_g.atom[0],
    top->ff_g[id].atom[1]=ff_g.atom[1],
    top->ff_g[id].atom[2]=ff_g.atom[2],
    top->ff_g[id].k=ff_g.k,
    top->ff_g[id].v=ff_g.v;
    }
  //manage the memory
  if (!((++top->size_g+0x1)%0xFF))
    {
    if (!(vp=realloc(top->ff_g,sizeof(t_ff_g)*(top->size_g+0x1+0xFF))))
      { lexem="angles"; goto LABEL_ERROR_MEMORY; }
    else top->ff_g=(t_ff_g*)vp;
    }
  }
//Sync angles
if (!(vp=realloc(top->ff_g,sizeof(t_ff_g)*(top->size_g))))
  { lexem="angles"; goto LABEL_ERROR_MEMORY; }
else top->ff_g=(t_ff_g*)vp;

//Stage 6. Upload ff impropers with 0xFF buffer
top->size_i=count=0, file_type="gaff impropers";
if (!(top->ff_i=(t_ff_i*)malloc(sizeof(t_ff_i)*0xFF))) 
  { lexem="impropers"; goto LABEL_ERROR_MEMORY; }
while (fgets(buffer,TOP_STR_LEN,in_imprs_gaff))
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(check_lexem(0x1,buffer))) continue;
  if ( (!(check_lexem(0x7,buffer)))||( (check_lexem(0x8,buffer))) ) 
    { lexem="wrong line format"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (!(lexem=get_lex(0x1,buffer)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wrong atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {//Correct for the atom order, for unknown reason gaff has the root atom third in a card
       _j=0; do ((char*)&ff_i.atom[2])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
       while (_j<sizeof(unsigned int)) ((char*)&ff_i.atom[2])[_j++]=' ';
       } 
  if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) )
    { lexem="wrong atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {
       _j=0; do ((char*)&ff_i.atom[1])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
       while (_j<sizeof(unsigned int)) ((char*)&ff_i.atom[1])[_j++]=' ';
       } 
  if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wrong atom#2 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {//Correct for the atom order, for unknown reason gaff has the root atom third in a card
       _j=0; do ((char*)&ff_i.atom[0])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
       while (_j<sizeof(unsigned int)) ((char*)&ff_i.atom[0])[_j++]=' ';
       } 
  if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wrong atom#3 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else {
       _j=0; do ((char*)&ff_i.atom[3])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
       while (_j<sizeof(unsigned int)) ((char*)&ff_i.atom[3])[_j++]=' ';
       } 
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong K"; goto LABEL_ERROR_DATA_FORMAT; }
  else ff_i.k=atof(lexem);
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong v0"; goto LABEL_ERROR_DATA_FORMAT; }
  else ff_i.v=atof(lexem);
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong n"; goto LABEL_ERROR_DATA_FORMAT; }
  //else ff_i.n=atof(lexem);
  //check atom types
  if ( (ff_i.atom[0]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_i.atom[0],top->size_n,top->ff_n))) )
    { lexem="unknown atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (ff_i.atom[1]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_i.atom[1],top->size_n,top->ff_n))) )
    { lexem="unknown atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (ff_i.atom[2]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_i.atom[2],top->size_n,top->ff_n))) )
    { lexem="unknown atom#2 name"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (ff_i.atom[3]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_i.atom[3],top->size_n,top->ff_n))) )
    { lexem="unknown atom#3 name"; goto LABEL_ERROR_DATA_FORMAT; }
  ff_i.type=1;
  //arrange impropers
  order_impr(ff_i.type,&ff_i.atom[0],&ff_i.atom[1],&ff_i.atom[2],&ff_i.atom[3],ff_i.atom[0],ff_i.atom[1],ff_i.atom[2],ff_i.atom[3]);
  //keep the massive sorted for rapid search for dublicates
  if ( (find_in_sorted_lth_objects(&id,&ff_i,top->size_i,(void*)top->ff_i,sizeof(t_ff_i),compare_imprs)))
    { lexem="dublicated gaff impropers"; goto LABEL_ERROR_DATA_FORMAT; }
  else
    {
    _j=top->size_i;
    while (id!=_j--)
      {
      top->ff_i[_j+1].atom[0]=top->ff_i[_j].atom[0];
      top->ff_i[_j+1].atom[1]=top->ff_i[_j].atom[1];
      top->ff_i[_j+1].atom[2]=top->ff_i[_j].atom[2];
      top->ff_i[_j+1].atom[3]=top->ff_i[_j].atom[3];
      top->ff_i[_j+1].k=top->ff_i[_j].k;
      top->ff_i[_j+1].v=top->ff_i[_j].v;
      top->ff_i[_j+1].type=top->ff_i[_j].type;
      }
    top->ff_i[id].atom[0]=ff_i.atom[0];
    top->ff_i[id].atom[1]=ff_i.atom[1];
    top->ff_i[id].atom[2]=ff_i.atom[2];
    top->ff_i[id].atom[3]=ff_i.atom[3];
    top->ff_i[id].k=ff_i.k;
    top->ff_i[id].v=ff_i.v;
    top->ff_i[id].type=ff_i.type;
    }
  //manage the memory
  if (!((++top->size_i+0x1)%0xFF))
    {
    if (!(vp=realloc(top->ff_i,sizeof(t_ff_i)*(top->size_i+0x1+0xFF))))
      { lexem="impropers"; goto LABEL_ERROR_MEMORY; }
    else top->ff_i=(t_ff_i*)vp;
    }
  }
//Sync imps
if (!(vp=realloc(top->ff_i,sizeof(t_ff_i)*(top->size_i))))
  { lexem="impropers"; goto LABEL_ERROR_MEMORY; }
else top->ff_i=(t_ff_i*)vp;

//Stage 7. Upload gaff dihedralss with 0xFF buffer using c as a dihedrals counter in complex dihedral
_c=0, top->size_d=count=0, file_type="gaff dihedrals";
if (!(top->ff_d=(t_ff_d*)malloc(sizeof(t_ff_d)*0xFF)))
  { lexem="dihedrals"; goto LABEL_ERROR_MEMORY; }
while (fgets(buffer,TOP_STR_LEN,in_dihs_gaff))
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(check_lexem(0x1,buffer))) continue;
  if ( (!(check_lexem(0x8,buffer)))||( (check_lexem(0x9,buffer))) ) 
    { lexem="wrong line format"; goto LABEL_ERROR_DATA_FORMAT; }
  //Process atom names
  if (!_c)
    {//Begin upload of a (potentially) complex torsion
    if ( (!(lexem=get_lex(0x1,buffer)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
      { lexem="wrong atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
    else {          
         _j=0; do ((char*)&ff_d.atom[0])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
         while (_j<sizeof(unsigned int)) ((char*)&ff_d.atom[0])[_j++]=' ';
         } 
    if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) )
      { lexem="wrong atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
    else {
         _j=0; do ((char*)&ff_d.atom[1])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
         while (_j<sizeof(unsigned int)) ((char*)&ff_d.atom[1])[_j++]=' ';
         } 
    if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
      { lexem="wrong atom#2 name"; goto LABEL_ERROR_DATA_FORMAT; }
    else {
         _j=0; do ((char*)&ff_d.atom[2])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
         while (_j<sizeof(unsigned int)) ((char*)&ff_d.atom[2])[_j++]=' ';
         } 
    if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
      { lexem="wrong atom#3 name"; goto LABEL_ERROR_DATA_FORMAT; }
    else {
         _j=0; do ((char*)&ff_d.atom[3])[_j]=tolower((int)lexem[_j]); while (++_j!=_i);
         while (_j<sizeof(unsigned int)) ((char*)&ff_d.atom[3])[_j++]=' ';
         }
    //check atom types
    if ( (ff_d.atom[0]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_d.atom[0],top->size_n,top->ff_n))) )
      { lexem="unknown atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
    if ( (ff_d.atom[1]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_d.atom[1],top->size_n,top->ff_n))) )
      { lexem="unknown atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
    if ( (ff_d.atom[2]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_d.atom[2],top->size_n,top->ff_n))) )
      { lexem="unknown atom#2 name"; goto LABEL_ERROR_DATA_FORMAT; }
    if ( (ff_d.atom[3]!=*((unsigned int*)&("x   ")))&&(!(find_in_sorted_lth_urow(FALSE,ff_d.atom[3],top->size_n,top->ff_n))) )
      { lexem="unknown atom#3 name"; goto LABEL_ERROR_DATA_FORMAT; }
    } 
  else 
    {//Continue upload of a complex torsion
    if ( (!(lexem=get_lex(0x1,buffer)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
      { lexem="wrong atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
    else {          
         _j=0; do if (((char*)&ff_d.atom[0])[_j]!=tolower((int)lexem[_j])) goto LABEL_COMPLEX_TORSION_NAME_FAILURE_0; while (++_j!=_i);
         while (_j<sizeof(unsigned int)) 
           if (((char*)&ff_d.atom[0])[_j++]!=' ') { LABEL_COMPLEX_TORSION_NAME_FAILURE_0: lexem="complex torsion atom#0 mismatch"; goto LABEL_ERROR_DATA_FORMAT; }
         } 
    if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) )
      { lexem="wrong atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
    else {
         _j=0; do if (((char*)&ff_d.atom[1])[_j]!=tolower((int)lexem[_j])) goto LABEL_COMPLEX_TORSION_NAME_FAILURE_1; while (++_j!=_i);
         while (_j<sizeof(unsigned int)) 
           if (((char*)&ff_d.atom[1])[_j++]!=' ') { LABEL_COMPLEX_TORSION_NAME_FAILURE_1: lexem="complex torsion atom#1 mismatch"; goto LABEL_ERROR_DATA_FORMAT; }
         } 
    if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
      { lexem="wrong atom#2 name"; goto LABEL_ERROR_DATA_FORMAT; }
    else {
         _j=0; do if (((char*)&ff_d.atom[2])[_j]!=tolower((int)lexem[_j])) goto LABEL_COMPLEX_TORSION_NAME_FAILURE_2; while (++_j!=_i);
         while (_j<sizeof(unsigned int)) 
           if (((char*)&ff_d.atom[2])[_j++]!=' ') { LABEL_COMPLEX_TORSION_NAME_FAILURE_2: lexem="complex torsion atom#2 mismatch"; goto LABEL_ERROR_DATA_FORMAT; }
         } 
    if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
      { lexem="wrong atom#3 name"; goto LABEL_ERROR_DATA_FORMAT; }
    else {
         _j=0; do if (((char*)&ff_d.atom[3])[_j]!=tolower((int)lexem[_j])) goto LABEL_COMPLEX_TORSION_NAME_FAILURE_3; while (++_j!=_i);
         while (_j<sizeof(unsigned int)) 
           if (((char*)&ff_d.atom[3])[_j++]!=' ') { LABEL_COMPLEX_TORSION_NAME_FAILURE_3: lexem="complex torsion atom#3 mismatch"; goto LABEL_ERROR_DATA_FORMAT; }
         }
    }
  //Process parameters 
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong D"; goto LABEL_ERROR_DATA_FORMAT; }
  else if (!(ff_d.n[(int)_c]=atof(lexem))) { lexem="zero D"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong K"; goto LABEL_ERROR_DATA_FORMAT; }
  else { ff_d.k[(int)_c]=atof(lexem); } //ff_d.k[(int)_c]/=ff_d.n[(int)_c]; }
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong v0"; goto LABEL_ERROR_DATA_FORMAT; }
  else ff_d.v[(int)_c]=atof(lexem);
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) 
    { lexem="wrong n"; goto LABEL_ERROR_DATA_FORMAT; }
  else 
    {
    ff_d.n[(int)_c]=atof(lexem);
    if (ff_d.n[(int)_c]!=(long long)ff_d.n[(int)_c]) { lexem="wrong non-integer N"; goto LABEL_ERROR_DATA_FORMAT; }
    }
  //Check for multiple torsions
  if (ff_d.n[(int)_c]<0.)
    { ff_d.n[(int)_c]=fabs(ff_d.n[(int)_c]); if (++_c==SIZE_DIH) { lexem="too many complex torsion replicas"; goto LABEL_ERROR_DATA_FORMAT; } }
  else
    {//Arrange uploaded Dihedral
    order_dih(&ff_d.atom[0],&ff_d.atom[1],&ff_d.atom[2],&ff_d.atom[3],ff_d.atom[0],ff_d.atom[1],ff_d.atom[2],ff_d.atom[3]);
    //Keep the massive sorted for rapid search for dublicates
    if ( (find_in_sorted_lth_objects(&id,&ff_d,top->size_d,(void*)top->ff_d,sizeof(t_ff_d),compare_dihs)))
      { lexem="dublicated gaff dihedrals"; goto LABEL_ERROR_DATA_FORMAT; }
    else
      {
      _j=top->size_d;
      while (id!=_j--)
        {
        top->ff_d[_j+1].atom[0]=top->ff_d[_j].atom[0];
        top->ff_d[_j+1].atom[1]=top->ff_d[_j].atom[1];
        top->ff_d[_j+1].atom[2]=top->ff_d[_j].atom[2];
        top->ff_d[_j+1].atom[3]=top->ff_d[_j].atom[3];
        _i=SIZE_DIH; 
        while (_i--) 
          {
          top->ff_d[_j+1].k[_i]=top->ff_d[_j].k[_i];
          top->ff_d[_j+1].v[_i]=top->ff_d[_j].v[_i];
          top->ff_d[_j+1].n[_i]=top->ff_d[_j].n[_i];
          }
        }
      top->ff_d[id].atom[0]=ff_d.atom[0];
      top->ff_d[id].atom[1]=ff_d.atom[1];
      top->ff_d[id].atom[2]=ff_d.atom[2];
      top->ff_d[id].atom[3]=ff_d.atom[3];
      _i=SIZE_DIH; while (--_i!=(int)_c) top->ff_d[id].k[_i]=top->ff_d[id].v[_i]=top->ff_d[id].n[_i]=0.;
      do { top->ff_d[id].k[(int)_c]=ff_d.k[(int)_c];
           top->ff_d[id].v[(int)_c]=ff_d.v[(int)_c];
           top->ff_d[id].n[(int)_c]=ff_d.n[(int)_c]; } while (_c--); _c=0;
      }    
    //Manage the memory
    if (!((++top->size_d+0x1)%0xFF))
      {
      if (!(vp=realloc(top->ff_d,sizeof(t_ff_d)*(top->size_d+0x1+0xFF)))) 
        { lexem="dihedrals"; goto LABEL_ERROR_MEMORY; }
      else top->ff_d=(t_ff_d*)vp;
      }
    }
  }
if ( (_c)) { lexem="incomplete complex dihedral"; goto LABEL_ERROR_DATA_FORMAT; } 
//Sync dihs 
if (!(vp=realloc(top->ff_d,sizeof(t_ff_d)*(top->size_d)))) 
  { lexem="dihedrals"; goto LABEL_ERROR_MEMORY; }
else top->ff_d=(t_ff_d*)vp;

//Stage 8. Upload residues topologies with 0xFF buffer
if (!(top->ress=alloc_list(0xFF)))
  { lexem="ressidues list"; goto LABEL_ERROR_MEMORY; }
if (!(top->res=(t_res*)malloc(sizeof(t_res)*0xFF)))
  { lexem="ressidues structure"; goto LABEL_ERROR_MEMORY; }
if (!(top->_ctypes=(char*)malloc(sizeof(char)*0xFF)))
  { lexem="residues chemtypes"; goto LABEL_ERROR_MEMORY; }
if (!(top->_atoms=(char(*)[sizeof(int)])malloc(sizeof(int)*0xFF)))
 { lexem="residues atoms"; goto LABEL_ERROR_MEMORY; }
if (!(top->_edges=(t_edge*)malloc(sizeof(t_edge)*0xFF)))
 { lexem="residues edges"; goto LABEL_ERROR_MEMORY; }
top->ress->size=top->_natoms=top->_nedges=count=0; file_type="yff1 residues";
while ((fgets(buffer,TOP_STR_LEN,in_ress)))
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(lexem=get_lex(0x1,buffer))) continue;
  if ( (lexem[0]=='[')&&(lexem[1]==' ')&&(lexem[2])&&(lexem[3])&&(lexem[4])&&(lexem[5])&&(lexem[6]==' ')&&(lexem[7]==']') ) goto UPLOAD_RESIDUE;
  else { lexem="loading the first residue card"; goto LABEL_ERROR_DATA_FORMAT; }
  }
{ lexem="absent the first residue card"; goto LABEL_ERROR_DATA_FORMAT; }

UPLOAD_RESIDUE:
lexem[2]=tolower((int)lexem[2]), lexem[3]=tolower((int)lexem[3]), lexem[4]=tolower((int)lexem[4]), lexem[5]=tolower((int)lexem[5]);
if (find_in_list(*((unsigned int*)&buffer[2]),top->ress)!=(unsigned int)-1) 
  { lexem="dublicated residue name"; goto LABEL_ERROR_DATA_FORMAT; }
//Prepare memory for a new residue
top->ress->list[top->ress->size]=*((unsigned int*)&lexem[2]);
top->res[top->ress->size].atoms.size=top->res[top->ress->size].nedges=0x0;
top->res[top->ress->size].atoms.list=(unsigned int*)&top->_atoms[top->_natoms];
top->res[top->ress->size].ctypes=&top->_ctypes[top->_natoms];
top->res[top->ress->size].edges=&top->_edges[top->_nedges];
//Upload residues atoms
while ( (fgets(buffer,TOP_STR_LEN,in_ress))) 
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(lexem=get_lex(0x1,buffer))) continue;
  if ( (lexem[0]=='[')&&(lexem[1]==' ')&&((lexem[2]=='A')||(lexem[2]=='a'))&&((lexem[3]=='T')||(lexem[3]=='t'))&&((lexem[4]=='O')||(lexem[4]=='o'))&&((lexem[5]=='M')||(lexem[5]=='m'))&&((lexem[6]=='S')||(lexem[6]=='s'))&&(lexem[7]==' ')&&(lexem[8]==']') )
    goto UPLOAD_RESIDUE_ATOMS;
  else { lexem="missing [ atom ] card declaration"; goto LABEL_ERROR_DATA_FORMAT; }
  }
{ lexem="incomplete residue card"; goto LABEL_ERROR_DATA_FORMAT; }
UPLOAD_RESIDUE_ATOMS:
while ( (fgets(buffer,TOP_STR_LEN,in_ress))) 
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(lexem=get_lex(0x1,buffer))) continue;
  if ( (lexem[0]=='[')&&(lexem[1]==' ')&&((lexem[2]=='B')||(lexem[2]=='b'))&&((lexem[3]=='O')||(lexem[3]=='o'))&&((lexem[4]=='N')||(lexem[4]=='n'))&&((lexem[5]=='D')||(lexem[5]=='d'))&&((lexem[6]=='S')||(lexem[6]=='s'))&&(lexem[7]==' ')&&(lexem[8]==']') )
    break; //Bonds section is comming  
  if (!(check_lexem(0x3,lexem))) 
    { lexem="wrong atom line format"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (!(check_int_type(lexem)))||((unsigned int)atoi(lexem)!=top->res[top->ress->size].atoms.size+1) ) 
    { lexem="wrong atom id"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wrong atom name"; goto LABEL_ERROR_DATA_FORMAT; }
  else 
    {
    _j=0; do { ((char*)&top->_atoms[top->_natoms])[_j]=tolower((int)lexem[_j]); } while (_i!=_j++);
    while (_j<sizeof(unsigned int)) ((char*)&top->_atoms[top->_natoms])[_j++]=' ';
    }
  if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem))) ) 
    { lexem="wrong atom chemtype"; goto LABEL_ERROR_DATA_FORMAT; }
  else top->_ctypes[top->_natoms]=(char)atoi(lexem);
    top->res[top->ress->size].atoms.size++;
  //Sync _res_natoms
  if (!(++top->_natoms%0xFF))
    {
    if (!(vp=realloc(top->_atoms,sizeof(unsigned int)*(top->_natoms+0xFF))))
      { lexem="residues atoms"; goto LABEL_ERROR_MEMORY; } 
    if (vp!=(void*)top->_atoms)
      {
      vp-=(size_t)top->_atoms;
      _i=top->ress->size;
      do { _vp=top->res[_i].atoms.list, _vp+=(size_t)vp, top->res[_i].atoms.list=_vp; } while (_i--);
      top->_atoms=(char(*)[sizeof(int)])top->res[0x0].atoms.list;
      }
    if (!(vp=realloc(top->_ctypes,sizeof(char)*(top->_natoms+0xFF))))
      { lexem="residues chemtypes"; goto LABEL_ERROR_MEMORY; }
    if (vp!=top->_ctypes)
      {
      vp-=(size_t)top->_ctypes;
      _i=top->ress->size;
      do { _vp=top->res[_i].ctypes, _vp+=(size_t)vp, top->res[_i].ctypes=_vp; } while (_i--);
      top->_ctypes=top->res[0x0].ctypes;
      } 
    }
  }
if (!top->res[top->ress->size].atoms.size) 
  { lexem="empty atoms list"; goto LABEL_ERROR_DATA_FORMAT; } //Empty atoms list of a residues is not permitted!
while ( (fgets(buffer,TOP_STR_LEN,in_ress))) //Empty bonds lists are permited in contradiction to banned empty atoms list
  {
  count++;
  CUT_COMMENTS(_i,buffer);
  if (!(lexem=get_lex(0x1,buffer))) continue; //Skip the empty strings
  else
    if ( (lexem[0]=='[')&&(lexem[1]==' ')&&(lexem[6]==' ')&&(lexem[7]==']') ) 
      {//Probably a new residue to upload
       //Sync ress
      if (!(++top->ress->size%0xFF))
        {
        if (!(vp=realloc_list(top->ress,top->ress->size+0xFF))) 
          { lexem="ressidues list"; goto LABEL_ERROR_MEMORY; }
        else top->ress=(t_list*)vp;
        }
      goto UPLOAD_RESIDUE;
      } 
  if (!(check_lexem(0x4,lexem))) 
    { lexem="wrong bonds line format"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (!(check_int_type(lexem)))||((unsigned int)atoi(lexem)!=top->res[top->ress->size].nedges+1) )
    { lexem="wrong bond id"; goto LABEL_ERROR_DATA_FORMAT; }
  if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wrong atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else 
    {
    _j=0; do ((char*)&top->_edges[top->_nedges].vertice[0x0])[_j]=tolower((int)lexem[_j]); while (_i!=_j++);
    while (_j<sizeof(unsigned int)) ((char*)&top->_edges[top->_nedges].vertice[0x0])[_j++]=' ';
    }  
  if ( (*lexem!='+')&&(*lexem!='-')&&(*lexem!='*') )
    {
    if ((top->_edges[top->_nedges].vertice[0]=find_in_list(top->_edges[top->_nedges].vertice[0],&top->res[top->ress->size].atoms))==(unsigned int)-1) 
      { lexem="absent atom#0 name"; goto LABEL_ERROR_DATA_FORMAT; }
    }
  else 
    {//Check for nonresidue linear bonds 
    if (*lexem!='*') 
      for (_i=0;_i<top->res[top->ress->size].nedges;_i++) 
        if (*lexem==*((char*)&top->res[top->ress->size].edges[_i].vertice[0]))  
          { lexem="only one \'+\' or \'-\' sign are allowed for linear polymerization"; goto LABEL_ERROR_DATA_FORMAT; } 
    }
  if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) 
    { lexem="wrong atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
  else 
    {
    _j=0; do ((char*)&top->_edges[top->_nedges].vertice[0x1])[_j]=tolower((int)lexem[_j]); while (_i!=_j++);
    while (_j<sizeof(unsigned int)) ((char*)&top->_edges[top->_nedges].vertice[0x1])[_j++]=' ';
    }
  if ( (*lexem!='+')&&(*lexem!='-')&&(*lexem!='*') )
    {
    if ((top->_edges[top->_nedges].vertice[1]=find_in_list(top->_edges[top->_nedges].vertice[1],&top->res[top->ress->size].atoms))==(unsigned int)-1) 
      { lexem="absent atom#1 name"; goto LABEL_ERROR_DATA_FORMAT; }
    }
  else 
    {
    if ( (*((char*)&top->res[top->ress->size].edges[_i].vertice[0])=='+')||
         (*((char*)&top->res[top->ress->size].edges[_i].vertice[0])=='-')||
         (*((char*)&top->res[top->ress->size].edges[_i].vertice[0])=='*') ) 
      { lexem="only transbonding is allowed"; goto LABEL_ERROR_DATA_FORMAT; } //No trans-bonding  permitted in residues
    if (*lexem!='*') //Check if the only one entrance of '+' and '-' 
      for (_i=0;_i<top->res[top->ress->size].nedges;_i++) 
        if (*lexem==*((char*)&top->res[top->ress->size].edges[_i].vertice[0])) 
          { lexem="only one \'+\' or \'-\' sign are allowed for linear polymerization"; goto LABEL_ERROR_DATA_FORMAT; }
    _i=top->_edges[top->_nedges].vertice[0], top->_edges[top->_nedges].vertice[0]=top->_edges[top->_nedges].vertice[1], top->_edges[top->_nedges].vertice[1]=_i;
    }
  if ( (!(lexem=get_lex(0x2,lexem)))||(lexlen(lexem)!=sizeof(char)) )
    { lexem="wrong bond type"; goto LABEL_ERROR_DATA_FORMAT; }
  if (_j) top->_edges[top->_nedges].type=-(*lexem);
  else    top->_edges[top->_nedges].type=+(*lexem);
  top->res[top->ress->size].nedges++;
  //Sync _res_nedges
  if (!(++top->_nedges%0xFF))
    {
    if (!(vp=realloc(top->_edges,sizeof(t_edge)*(top->_nedges+0xFF)))) 
      { lexem="residues edges"; goto LABEL_ERROR_MEMORY; }
    if (vp!=top->_edges)
      {
      vp-=(size_t)top->_edges;
      _i=top->ress->size;
      do { _vp=top->res[_i].edges, _vp+=(size_t)vp, top->res[_i].edges=_vp; } while (_i--);
      top->_edges=top->res[0x0].edges;
      }
    }
  }
//Sync ress
top->ress->size++;
if (!(vp=realloc_list(top->ress,top->ress->size))) 
  { lexem="ressidues list"; goto LABEL_ERROR_MEMORY; }
else top->ress=(t_list*)vp;
if (!(vp=realloc(top->res,sizeof(t_res)*top->ress->size))) 
  { lexem="ressidues structure"; goto LABEL_ERROR_MEMORY; }
else top->res=(t_res*)vp;
//Sync atoms, types and bonds
if (!(vp=realloc(top->_ctypes,sizeof(char)*top->_natoms))) 
  { lexem="residues chemtypes"; goto LABEL_ERROR_MEMORY; }
if (vp!=top->_ctypes)
  {
  vp-=(size_t)top->_ctypes;
  _i=top->ress->size;
  while (_i--) { _vp=top->res[_i].ctypes, _vp+=(size_t)vp, top->res[_i].ctypes=_vp; }
  top->_ctypes=top->res[0x0].ctypes;
  }
if (!(vp=realloc(top->_atoms,sizeof(unsigned int)*top->_natoms)))
  { lexem="residues atoms"; goto LABEL_ERROR_MEMORY; }
if (vp!=top->_atoms)
  {
  vp-=(size_t)top->_atoms;
  _i=top->ress->size;
  while (_i--) { _vp=top->res[_i].atoms.list, _vp+=(size_t)vp, top->res[_i].atoms.list=_vp; }
  top->_atoms=(char(*)[sizeof(int)])top->res[0x0].atoms.list;
  }
if (!(vp=realloc(top->_edges,sizeof(t_edge)*top->_nedges)))
 { lexem="residues edges"; goto LABEL_ERROR_MEMORY; }
if (vp!=top->_edges)
  {
  vp-=(size_t)top->_edges;
  _i=top->ress->size;
  while (_i--) { _vp=top->res[_i].edges, _vp+=(size_t)vp, top->res[_i].edges=_vp; }
  top->_edges=top->res[0x0].edges;
  }

//Stage 9. Store various constants

//Stage 9.1. Store GAFF empirical parameters
//NOTE. The atom types
// H  0
// C  1
// N  2
// O  3
// F  4
// Cl 5
// Br 6
// I  7
// P  8
// S  9
// Si 8 (*)
//Bonds table
//Stage 9.1.1. GAFF empirical parameters for bonds
//The table (Table 3) from Wang et al "Development and testing of a general amber force field" JComputChem 2004
//  #  A1 A2  rij    lnKij
//  1  H  H   0.738  4.661
//  2  C  C   1.526  7.643
//  3  N  N   1.441  7.634
//  4  O  O   1.460  7.561
//  5  F  F   1.406  7.358
//  6  Cl Cl  2.031  8.648
//  7  Br Br  2.337  9.012
//  8  I  I   2.836  9.511
//  9  P  P   2.324  8.805
// 10  S  S   2.038  8.316
// 11  H  C   1.090  6.217
// 12  H  N   1.010  6.057
// 13  H  O   0.960  5.794
// 14  H  F   0.920  5.600
// 15  H  Cl  1.280  6.937
// 16  H  Br  1.410  7.301
// 17  H  I   1.600  7.802
// 18  H  P   1.410  7.257
// 19  H  S   1.340  7.018
// 20  C  N   1.470  7.504
// 21  C  O   1.440  7.347
// 22  C  F   1.370  7.227
// 23  C  Cl  1.800  8.241
// 24  C  Br  1.940  8.478
// 25  C  I   2.160  8.859
// 26  C  P   1.830  8.237
// 27  C  S   1.820  8.117
// 28  N  O   1.420  7.526
// 29  N  F   1.420  7.475
// 30  N  Cl  1.750  8.266
// 31  N  Br  1.930  8.593
// 32  N  I   2.120  8.963
// 33  N  P   1.720  8.212
// 34  N  S   1.690  8.073
// 35  O  F   1.410  7.375
// 36  O  Cl  1.700  8.097
// 37  O  Br  1.790  8.276
// 38  O  I   2.110  8.854
// 39  O  P   1.640  7.957
// 40  O  S   1.650  7.922
// 41  F  Cl  1.648  7.974
// 42  F  Br  1.872  8.185
// 43  F  I   2.121  8.435
// 44  F  P   1.500  7.592
// 45  F  S   1.580  7.733
// 46  Cl Br  2.184  8.830
// 47  Cl I   2.550  9.309
// 48  Cl P   2.040  8.656
// 49  Cl S   2.030  8.619
// 50  Br I   2.671  9.380
// 51  Br P   2.240  8.729
// 52  Br S   2.210  8.728
// 53  I  P   2.490  9.058
// 54  I  S   2.560  9.161
// 55  P  S   2.120  8.465
top->br[0][0]=0.738, top->bK[0][0]=exp(4.661);
top->br[0][1]=top->br[1][0]=1.090, top->bK[0][1]=top->bK[1][0]=exp(6.217);
top->br[0][2]=top->br[2][0]=1.010, top->bK[0][2]=top->bK[2][0]=exp(6.057);
top->br[0][3]=top->br[3][0]=0.960, top->bK[0][3]=top->bK[3][0]=exp(5.794);
top->br[0][4]=top->br[4][0]=0.920, top->bK[0][4]=top->bK[4][0]=exp(5.600);
top->br[0][5]=top->br[5][0]=1.280, top->bK[0][5]=top->bK[5][0]=exp(6.937);
top->br[0][6]=top->br[6][0]=1.410, top->bK[0][6]=top->bK[6][0]=exp(7.301);
top->br[0][7]=top->br[7][0]=1.600, top->bK[0][7]=top->bK[7][0]=exp(7.802);
top->br[0][8]=top->br[8][0]=1.410, top->bK[0][8]=top->bK[8][0]=exp(7.257);
top->br[0][9]=top->br[9][0]=1.340, top->bK[0][9]=top->bK[9][0]=exp(7.018);
top->br[1][1]=1.526, top->bK[1][1]=exp(7.643);
top->br[1][2]=top->br[2][1]=1.470, top->bK[1][2]=top->bK[2][1]=exp(7.504);
top->br[1][3]=top->br[3][1]=1.440, top->bK[1][3]=top->bK[3][1]=exp(7.347);
top->br[1][4]=top->br[4][1]=1.370, top->bK[1][4]=top->bK[4][1]=exp(7.227);
top->br[1][5]=top->br[5][1]=1.800, top->bK[1][5]=top->bK[5][1]=exp(8.241);
top->br[1][6]=top->br[6][1]=1.940, top->bK[1][6]=top->bK[6][1]=exp(8.478);
top->br[1][7]=top->br[7][1]=2.160, top->bK[1][7]=top->bK[7][1]=exp(8.859);
top->br[1][8]=top->br[8][1]=1.830, top->bK[1][8]=top->bK[8][1]=exp(8.237);
top->br[1][9]=top->br[9][1]=1.820, top->bK[1][9]=top->bK[9][1]=exp(8.117);
top->br[2][2]=1.441, top->bK[2][2]=exp(7.634);
top->br[2][3]=top->br[3][2]=1.420, top->bK[2][3]=top->bK[3][2]=exp(7.526);
top->br[2][4]=top->br[4][2]=1.420, top->bK[2][4]=top->bK[4][2]=exp(7.475);
top->br[2][5]=top->br[5][2]=1.750, top->bK[2][5]=top->bK[5][2]=exp(8.266);
top->br[2][6]=top->br[6][2]=1.930, top->bK[2][6]=top->bK[6][2]=exp(8.593);
top->br[2][7]=top->br[7][2]=2.120, top->bK[2][7]=top->bK[7][2]=exp(8.963);
top->br[2][8]=top->br[8][2]=1.720, top->bK[2][8]=top->bK[8][2]=exp(8.212);
top->br[2][9]=top->br[9][2]=1.690, top->bK[2][9]=top->bK[9][2]=exp(8.073);
top->br[3][3]=1.460, top->bK[3][3]=exp(7.561);
top->br[3][4]=top->br[4][3]=1.410, top->bK[3][4]=top->bK[4][3]=exp(7.375);
top->br[3][5]=top->br[5][3]=1.700, top->bK[3][5]=top->bK[5][3]=exp(8.097);
top->br[3][6]=top->br[6][3]=1.790, top->bK[3][6]=top->bK[6][3]=exp(8.276);
top->br[3][7]=top->br[7][3]=2.110, top->bK[3][7]=top->bK[7][3]=exp(8.854);
top->br[3][8]=top->br[8][3]=1.640, top->bK[3][8]=top->bK[8][3]=exp(7.957);
top->br[3][9]=top->br[9][3]=1.650, top->bK[3][9]=top->bK[9][3]=exp(7.922);
top->br[4][4]=1.406, top->bK[4][4]=exp(7.358);
top->br[4][5]=top->br[5][4]=1.648, top->bK[4][5]=top->bK[5][4]=exp(7.974);
top->br[4][6]=top->br[6][4]=1.872, top->bK[4][6]=top->bK[6][4]=exp(8.185);
top->br[4][7]=top->br[7][4]=2.121, top->bK[4][7]=top->bK[7][4]=exp(8.435);
top->br[4][8]=top->br[8][4]=1.500, top->bK[4][8]=top->bK[8][4]=exp(7.592);
top->br[4][9]=top->br[9][4]=1.580, top->bK[4][9]=top->bK[9][4]=exp(7.733);
top->br[5][5]=2.031, top->bK[5][5]=exp(8.648);
top->br[5][6]=top->br[6][5]=2.184, top->bK[5][6]=top->bK[6][5]=exp(8.830);
top->br[5][7]=top->br[7][5]=2.550, top->bK[5][7]=top->bK[7][5]=exp(9.309);
top->br[5][8]=top->br[8][5]=2.040, top->bK[5][8]=top->bK[8][5]=exp(8.656);
top->br[5][9]=top->br[9][5]=2.030, top->bK[5][9]=top->bK[9][5]=exp(8.619);
top->br[6][6]=2.337, top->bK[6][6]=exp(9.012);
top->br[6][7]=top->br[7][6]=2.671, top->bK[6][7]=top->bK[7][6]=exp(9.380);
top->br[6][8]=top->br[8][6]=2.240, top->bK[6][8]=top->bK[8][6]=exp(8.729);
top->br[6][9]=top->br[9][6]=2.210, top->bK[6][9]=top->bK[9][6]=exp(8.728);
top->br[7][7]=2.836, top->bK[7][7]=exp(9.511);
top->br[7][8]=top->br[8][7]=2.490, top->bK[7][8]=top->bK[8][7]=exp(9.058);
top->br[7][9]=top->br[9][7]=2.560, top->bK[7][9]=top->bK[9][7]=exp(9.161);
top->br[8][8]=2.324, top->bK[8][8]=exp(8.805);
top->br[8][9]=top->br[9][8]=2.120, top->bK[8][9]=top->bK[9][8]=exp(8.465);
top->br[9][9]=2.038, top->bK[9][9]=exp(8.316);
//Stage 9.1.2. GAFF empirical parameters for angles
//The table (Table 4) from Wang et al "Development and testing of a general amber force field" JComputChem 2004
// a[ABC]=(a[ABA]+a[CBC])/2, K[ABC]=143.9*Z[A]*C[B]*Z[C]*exp<-2*(r[AB]-r[BC])^2/(r[AB]+r[BC])^2>/((r[AB]+r[BC])*(a[ABC])^2)
// #  A    C      Z
// 1  H    -     0.784
// 2  C   1.339  1.183
// 3  N   1.300  1.212
// 4  O   1.249  1.219
// 5  F   1.*    1.166
// 6  Cl  1.*    1.272
// 7  Br  1.*    1.378
// 8  I   1.*    1.398
// 9  P   0.906  1.620
//10  S   1.448  1.280
top->gC[0]=(double)NAN, top->gZ[0]=0.784;
top->gC[1]=1.339,       top->gZ[1]=1.183;
top->gC[2]=1.300,       top->gZ[2]=1.212;
top->gC[3]=1.249,       top->gZ[3]=1.219;
top->gC[4]=1.,          top->gZ[4]=1.166;
top->gC[5]=1,           top->gZ[5]=1.272;
top->gC[6]=1.,          top->gZ[6]=1.378;
top->gC[7]=1.,          top->gZ[7]=1.398;
top->gC[8]=0.906,       top->gZ[8]=1.620;
top->gC[9]=1.448,       top->gZ[9]=1.280;

//Stage 9.2. Convert GAFF constants 
_i=top->size_b; while (_i--) { top->ff_b[_i].k*=4.184; }                                    // kcal/mol/A^2  -> kJ/mol/A^2
_i=top->size_g; while (_i--) { top->ff_g[_i].k*=4.184, top->ff_g[_i].v*=PI/180.; }          // kcal/mol/rad^2 -> kJ/mol/Deg^2
_i=top->size_i; while (_i--) { top->ff_i[_i].k*=4.184*180./PI, top->ff_i[_i].v*=PI/180.; }  // kcal/mol/rad^2 -> kJ/mol/Deg^2, Deg -> Rad ??        
_i=top->size_d; while (_i--) { _j=SIZE_DIH; while (_j--) { top->ff_d[_i].k[_j]*=4.184, top->ff_d[_i].v[_j]*=PI/180.; } } // kcal/mol -> kJ/mol, Deg -> Rad         

//Job finished
return top;
}


//This function performs default topology import from ../top/*
t_top *import_top()
{
FILE *in[0x8];
t_top *top=0x0;

if (!(in[0x0]=fopen("/home/ayakovenko/y_project/top/yff1.a","r"))) { yprintf(YPRINTF_ERROR,"Can't open file yff1.a for reading");
ylib_errno=YERROR_USER; return FALSE; }
if (!(in[0x1]=fopen("/home/ayakovenko/y_project/top/yff1.p","r"))) { yprintf(YPRINTF_ERROR,"Can't open file yff1.p for reading"); ylib_errno=YERROR_USER; goto LABEL_ERROR_0; }
if (!(in[0x2]=fopen("/home/ayakovenko/y_project/top/gaff.a","r"))) { yprintf(YPRINTF_ERROR,"Can't open file gaff.a for reading"); ylib_errno=YERROR_USER; goto LABEL_ERROR_1; }
if (!(in[0x3]=fopen("/home/ayakovenko/y_project/top/gaff.b","r"))) { yprintf(YPRINTF_ERROR,"Can't open file gaff.b for reading"); ylib_errno=YERROR_USER; goto LABEL_ERROR_2; }
if (!(in[0x4]=fopen("/home/ayakovenko/y_project/top/gaff.g","r"))) { yprintf(YPRINTF_ERROR,"Can't open file gaff.g for reading"); ylib_errno=YERROR_USER; goto LABEL_ERROR_3; }
if (!(in[0x5]=fopen("/home/ayakovenko/y_project/top/gaff.i","r"))) { yprintf(YPRINTF_ERROR,"Can't open file gaff.i for reading"); ylib_errno=YERROR_USER; goto LABEL_ERROR_4; }
if (!(in[0x6]=fopen("/home/ayakovenko/y_project/top/gaff.d","r"))) { yprintf(YPRINTF_ERROR,"Can't open file gaff.d for reading"); ylib_errno=YERROR_USER; goto LABEL_ERROR_5; }
if (!(in[0x7]=fopen("/home/ayakovenko/y_project/top/yff1.r","r"))) { yprintf(YPRINTF_ERROR,"Can't open file yff1.r for reading"); ylib_errno=YERROR_USER; goto LABEL_ERROR_6; }
top=read_top(in[0x0],in[0x1],in[0x2],in[0x3],in[0x4],in[0x5],in[0x6],in[0x7]);
               fclose(in[0x7]);
LABEL_ERROR_6: fclose(in[0x6]);
LABEL_ERROR_5: fclose(in[0x5]);
LABEL_ERROR_4: fclose(in[0x4]);
LABEL_ERROR_3: fclose(in[0x3]);
LABEL_ERROR_2: fclose(in[0x2]);
LABEL_ERROR_1: fclose(in[0x1]);
LABEL_ERROR_0: fclose(in[0x0]);
return top;
}



//************************************    I / O     P R O C E S S O R S       *********************************/

//This function free memory of str
//Note it doesn't delete t_str itself
inline void free_str(register t_str *str)
{
if ( (str))
  {
  if ( (str->name))   free(str->name);
  if ( (str->ress)) 
    {
    free(str->ress);
    if ( (str->rsize))  free(str->rsize);
    }
  else str->rsize=0x0;
  if ( (str->a))      free(str->a); 
  if ( (str->anames)) free(str->anames); 
  if ( (str->r))      free(str->r); 
  if ( (str->edges))  free(str->edges);
  }
}

//This function free mol structure
void free_mol(t_mol *mol)
{
if (mol)
  {
  //Stage 1. free YSTR data
  if (mol->name)      free(mol->name),    mol->name=0x0;
  if (mol->ress)      free(mol->ress),    mol->ress=0x0;
  if (mol->rsize)     free(mol->rsize),   mol->rsize=0x0;
  if (mol->anames)    free(mol->anames),  mol->anames=0x0;
  if (mol->a)         free(mol->a),       mol->a=0x0;
  if (mol->r)         free(mol->r),       mol->r=0x0;
  if (mol->edges)     free(mol->edges),   mol->edges=0x0;
  //Stage 2. free Tolology data
  if (mol->anchors)   free(mol->anchors), mol->anchors=0x0;
  if (mol->aedges)    free(mol->aedges),  mol->aedges=0x0;
  if (mol->cycles)    free(mol->cycles),  mol->cycles=0x0;
  if ( (mol->fftypes)&&(mol->fftypes!=mol->ytypes) ) free(mol->fftypes), mol->fftypes=0x0;
  if (mol->ytypes)    free(mol->ytypes),  mol->ytypes=0x0;
  //Stage 3. free CDS data
       if (mol->cmtype==+1) { if (mol->C.dL)         free(mol->C.dL), mol->C.dL=0x0; }
  else if (mol->cmtype==-1) { if (mol->C.sL) free_smatrix(mol->C.sL), mol->C.sL=0x0; }
  if (mol->vatoms)    free(mol->vatoms),  mol->vatoms=0x0;
  if (mol->vedges)    free(mol->vedges),  mol->vedges=0x0;
  if (mol->engs)      free(mol->engs),    mol->engs=0x0;
  if (mol->hrds)      free(mol->hrds),    mol->hrds=0x0;
  if (mol->charges)   free(mol->charges), mol->charges=0x0;
  //Stage 4. free FF data
  if (mol->ff_b)      free(mol->ff_b),    mol->ff_b=0x0;
  if (mol->ff_g)      free(mol->ff_g),    mol->ff_g=0x0;
  if (mol->ff_i)      free(mol->ff_i),    mol->ff_i=0x0;
  if (mol->ff_d)      free(mol->ff_d),    mol->ff_d=0x0;
  if (mol->excl)      free(mol->excl),    mol->excl=0x0;
  if (mol->ff_c)      free(mol->ff_c),    mol->ff_c=0x0;
  free(mol), mol=0x0;
  }
}

//This function allocate memory for str
inline t_str *alloc_str(t_str *str,register unsigned int size_r,register unsigned int natoms,register unsigned int nedges)
{
if (!natoms) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
if (!(str))
  {
  if (!(str=(t_str*)calloc(1,sizeof(t_str))))                            { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return FALSE; }
  if (!(str->a=(char*)malloc(sizeof(char)*natoms)))                      { LABEL_MEMORY_ERROR_1a: free(str);    str=0x0;              goto LABEL_MEMORY_ERROR_0;  } 
  if (!(str->anames=(char (*)[sizeof(int)])malloc(sizeof(int)*natoms)))  { LABEL_MEMORY_ERROR_2a: free(str->a); str->a=0x0;           goto LABEL_MEMORY_ERROR_1a; }
  if (!(str->r=(t_vec*)malloc(sizeof(t_vec)*natoms)))                    { LABEL_MEMORY_ERROR_3a: free(str->anames); str->anames=0x0; goto LABEL_MEMORY_ERROR_2a; }
  if ( (nedges))
    {
    if (!(str->edges=(t_edge*)malloc(sizeof(t_edge)*nedges)))            { LABEL_MEMORY_ERROR_4a: free(str->r); str->r=0x0;           goto LABEL_MEMORY_ERROR_3a; }
    if ( (!size_r)||(size_r==1) )
      {
      str->ress=0x0;
      str->rsize=0x0;
      }
    else 
      {
      if (!(str->rsize=(unsigned int*)malloc(sizeof(unsigned int)*(size_r+0x1)))) { LABEL_MEMORY_ERROR_4d: free(str->edges); str->edges=0x0; goto LABEL_MEMORY_ERROR_4a; }
      if (!(str->ress=(t_list*)alloc_list(size_r))) { free(str->rsize); str->rsize=0x0; goto LABEL_MEMORY_ERROR_4d; }
      }
    }
  else
    {  
    if ( (!size_r)||(size_r==1) )
      {
      str->ress=0x0;
      str->rsize=0x0;
      }
    else 
      {
      if (!(str->rsize=(unsigned int*)malloc(sizeof(unsigned int)*(size_r+0x1))))       goto LABEL_MEMORY_ERROR_4a;
      if (!(str->ress=(t_list*)alloc_list(size_r))) { free(str->rsize); str->rsize=0x0; goto LABEL_MEMORY_ERROR_4a; }
      }
    }
  }
else
  {
  if (!(str->a=(char*)malloc(sizeof(char)*natoms)))                      { LABEL_MEMORY_ERROR_1b: goto LABEL_MEMORY_ERROR_0; } 
  if ( (str->anames)) { free(str->anames); str->anames=0x0; } 
  if (!(str->anames=(char (*)[sizeof(int)])malloc(sizeof(int)*natoms)))  { LABEL_MEMORY_ERROR_2b: free(str->a); str->a=0x0;           goto LABEL_MEMORY_ERROR_1b; }
  if ( (str->r)) { free(str->r); str->r=0x0; } 
  if (!(str->r=(t_vec*)malloc(sizeof(t_vec)*natoms)))                    { LABEL_MEMORY_ERROR_3b: free(str->anames); str->anames=0x0; goto LABEL_MEMORY_ERROR_2b; }
  if ( (nedges))
    {
    if ( (str->edges)) { free(str->edges); str->edges=0x0; } 
    if (!(str->edges=(t_edge*)malloc(sizeof(t_edge)*nedges)))            { LABEL_MEMORY_ERROR_4b: free(str->r); str->r=0x0;           goto LABEL_MEMORY_ERROR_3b; }
    if ( (str->rsize)) { free(str->rsize); str->rsize=0x0; } 
    if ( (str->ress))  { free(str->ress);  str->ress=0x0;  } 
    if ( (!size_r)||(size_r==1) )
      {
      str->ress=0x0;
      str->rsize=0x0;
      }
    else 
      {
      if (!(str->rsize=(unsigned int*)malloc(sizeof(unsigned int)*(size_r+0x1)))) { LABEL_MEMORY_ERROR_4e: free(str->edges); str->edges=0x0; goto LABEL_MEMORY_ERROR_4b; }
      if (!(str->ress=(t_list*)alloc_list(size_r))) { free(str->rsize); str->rsize=0x0; goto LABEL_MEMORY_ERROR_4e; }
      }
    }
  else
    {  
    if ( (str->rsize)) { free(str->rsize); str->rsize=0x0; } 
    if ( (str->ress))  { free(str->ress);  str->ress=0x0;  } 
    if ( (!size_r)||(size_r==1) )
      {
      str->ress=0x0;
      str->rsize=0x0;
      }
    else 
      {
      if (!(str->rsize=(unsigned int*)malloc(sizeof(unsigned int)*(size_r+0x1))))       goto LABEL_MEMORY_ERROR_4b;
      if (!(str->ress=(t_list*)alloc_list(size_r))) { free(str->rsize); str->rsize=0x0; goto LABEL_MEMORY_ERROR_4b; }
      }
    }
  }
return str;
}
//This function allocate solid str
inline t_str *alloc_solid_str(register unsigned int len,register unsigned int size_r,register unsigned int natoms,register unsigned int nedges)
{
register t_str *str;
register size_t align1, align2;
//Allocate solid structure in memory
if ( ((len+1)%sizeof(unsigned int))) align1=sizeof(unsigned int)-((len+1)%sizeof(unsigned int)); else align1=0x0;
if ( (natoms%sizeof(unsigned int))) align2=sizeof(unsigned int)-(natoms%sizeof(unsigned int)); else align2=0x0; 
if (!(str=(t_str*)malloc(sizeof(t_str)+sizeof(char)*(len+1)+align1+((size_r==1) ? 0x0 : sizeof(t_list)+sizeof(unsigned int)*(size_r+size_r+1))+(sizeof(t_vec)+sizeof(int))*natoms+sizeof(char)*natoms+align2+sizeof(t_edge)*nedges))) { ylib_errno=YERROR_MEMORY; return FALSE; } 
//Map memory
str->name=(char*)((void*)str+sizeof(t_str));
if (size_r==1)
  {
  str->ress=0x0, str->rsize=0x0;
  str->anames=(char(*)[4])((void*)str->name+sizeof(char)*(len+1)+align1);
  }
else
  {
  str->ress=(t_list*)((void*)str->name+sizeof(char)*(len+1)+align1);
  str->ress->size=size_r, str->ress->list=(unsigned int*)((void*)str->ress+sizeof(t_list));
  str->rsize=(unsigned int*)((void*)str->ress->list+sizeof(unsigned int)*size_r);
  str->anames=(char(*)[4])((void*)str->rsize+sizeof(unsigned int)*(size_r+1));
  }
str->a=(char*)((void*)str->anames+sizeof(int)*natoms);
str->r=(t_vec*)((void*)str->a+sizeof(char)*natoms+align2);
if (nedges) str->edges=(t_edge*)((void*)str->r+sizeof(t_vec)*natoms);
else str->edges=0x0;
return str;
}

//This function define chemical type from the name
unsigned char name_to_chemid(register char *lexem)
{
//Get atoms type
while (*lexem)
  if ( ((*lexem>='a')&&(*lexem<='z'))||((*lexem>='A')&&(*lexem<='Z')) )
    break;
  else
    lexem++;
switch(*lexem)
  {
  case 'A' : 
  case 'a' : {
             switch (*(++lexem))
               {
               case 'G' :
               case 'g' : return 47; //Ag
               case 'L' :
               case 'l' : return 13; //Al
               case 'R' :
               case 'r' : return 18; //Ar
               case 'S' :
               case 's' : return 33; //As
               case 'T' :
               case 't' : return 85; //At
               case 'U' :
               case 'u' : return 79; //Au
               default  : return FALSE;
               }
             }
  case 'B' : 
  case 'b' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return  5; //B
               case 'A' :
               case 'a' : return 56; //Ba
               case 'E' :
               case 'e' : return  4; //Be
               case 'I' :
               case 'i' : return 83; //Bi
               case 'R' :
               case 'r' : return 35; //Br
               default  : return FALSE;
               }
             }
  case 'C' : 
  case 'c' : {
             switch (*(++lexem))
               {
               case'\t' :
               case ' ' :
               case '.' : return  6; //C
               case 'A' :
               case 'a' : return 20; //Ca
               case 'D' :
               case 'd' : return 48; //Cd
               case 'E' :
               case 'e' : return 58; //Ce
               case 'L' :
               case 'l' : return 17; //Cl
               case 'O' :
               case 'o' : return 27; //Co
               case 'R' :
               case 'r' : return 24; //Cr
               case 'S' :
               case 's' : return 55; //Cs
               case 'U' :
               case 'u' : return 29; //Cu
               default  : return FALSE;
               }
             }
  case 'D' :
  case 'd' : {
             switch (*(++lexem))
               {
               case 'Y' :
               case 'y' : return 66; //Dy
               default  : return FALSE;
               }
             }
  case 'E' :
  case 'e' : {
             switch (*(++lexem))
               {
               case 'R' :
               case 'r' : return 68; //Er
               case 'U' :
               case 'u' : return 63; //Eu
               default  : return FALSE;
               }
             }
  case 'F' :
  case 'f' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return  9; //F
               case 'E' :
               case 'e' : return 26; //Fe
               default  : return FALSE;
               }
             }
  case 'G' :
  case 'g' : {
             switch (*(++lexem))
               {
               case 'A' :
               case 'a' : return 31; //Ga
               case 'D' :
               case 'd' : return 64; //Gd
               case 'E' :
               case 'e' : return 32; //Ge
               default  : return FALSE;
               }
             }
  case 'H' :
  case 'h' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return  1; //H
               case 'E' :
               case 'e' : return  2; //He
               case 'F' :
               case 'f' : return 72; //Hf
               case 'G' :
               case 'g' : return 80; //Hg
               case 'O' :
               case 'o' : return 67; //Ho
               default  : return FALSE;
               }
             }
  case 'I' :
  case 'i' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return 53; //I
               case 'N' :
               case 'n' : return 49; //In
               case 'R' :
               case 'r' : return 77; //Ir
               default  : return FALSE;
               }
             }
  case 'K' :
  case 'k' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return 19; //K
               case 'R' :
               case 'r' : return 36; //Kr
               default  : return FALSE;
               }
             }
  case 'L' :
  case 'l' : {
             switch (*(++lexem))
               {
               case 'A' :
               case 'a' : return  57; //La
               case 'I' :
               case 'i' : return  3; //Li
               case 'U' :
               case 'u' : return 71; //Lu
               default  : return FALSE;
               }
             }
  case 'M' :
  case 'm' : {
             switch (*(++lexem))
               {
               case 'G' :
               case 'g' : return 12; //Mg
               case 'N' :
               case 'n' : return 25; //Mm
               case 'O' :
               case 'o' : return 42; //Mo
               default  : return FALSE;
               }
             }
  case 'N' : 
  case 'n' : {
             switch (*(++lexem))
               {
               case'\t' :
               case ' ' :
               case '.' : return  7; //N
               case 'A' :
               case 'a' : return 11; //Na
               case 'B' :
               case 'b' : return 41; //Nb
               case 'D' :
               case 'd' : return 60; //Nd
               case 'E' :
               case 'e' : return 10; //Ne
               case 'I' :
               case 'i' : return 28; //Ni
               default  : return FALSE;
               }
             }
  case 'O' :
  case 'o' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return  8; //O
               case 'S' :
               case 's' : return 76; //Os
               default  : return FALSE;
               }
             }
  case 'P' :
  case 'p' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return 15; //P
               case 'B' :
               case 'b' : return 82; //Pb
               case 'D' :
               case 'd' : return 46; //Pd
               case 'M' :
               case 'm' : return 61; //Pm
               case 'O' :
               case 'o' : return 84; //Po
               case 'R' :
               case 'r' : return 59; //Pr
               case 'T' :
               case 't' : return 78; //Pt
               default  : return FALSE;
               }
             }
  case 'R' :
  case 'r' : {
             switch (*(++lexem))
               {
               case 'B' :
               case 'b' : return 37; //Rb
               case 'E' :
               case 'e' : return 75; //Re
               case 'H' :
               case 'h' : return 45; //Rh
               case 'N' :
               case 'n' : return 86; //Rn
               case 'U' :
               case 'u' : return 44; //Ru
               default  : return FALSE;
               }
             }
  case 'S' : 
  case 's' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return 16; //S
               case 'B' :
               case 'b' : return 51; //Sb
               case 'C' :
               case 'c' : return 21; //Sc
               case 'E' :
               case 'e' : return 34; //Se
               case 'I' :
               case 'i' : return 14; //Si
               case 'M' :
               case 'm' : return 62; //Sm
               case 'N' :
               case 'n' : return 50; //Sn
               case 'R' :
               case 'r' : return 38; //Sr
               default  : return FALSE;
               }
             }
  case 'T' :
  case 't' : {
             switch (*(++lexem))
               {
               case 'A' :
               case 'a' : return 73; //Ta
               case 'B' :
               case 'b' : return 65; //Tb
               case 'C' :
               case 'c' : return 43; //Tc
               case 'E' :
               case 'e' : return 52; //Te
               case 'I' :
               case 'i' : return 22; //Ti
               case 'L' :
               case 'l' : return 81; //Tl
               case 'M' :
               case 'm' : return 69; //Tm
               default  : return FALSE;
               }
             }
  case 'V' :
  case 'v' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' : 
               case '.' : return 23; //V
               default  : return FALSE;
               }
             }
  case 'W' :
  case 'w' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return 74; //W
               default  : return FALSE;
               }
             }
  case 'X' :
  case 'x' : {
             switch (*(++lexem))
               {
               case 'E' :
               case 'e' : return 54; //Xe
               default  : return FALSE;
               }
             }
  case 'Y' :
  case 'y' : {
             switch (*(++lexem))
               {
               case '\t' :
               case ' ' :
               case '.' : return 39; //Y
               case 'B' :
               case 'b' : return 70; //Yb
               default  : return FALSE;
               }
             }
  case 'Z' :
  case 'z' : {
             switch (*(++lexem))
               {
               case 'N' :
               case 'n' : return 30; //Zn
               case 'R' :
               case 'r' : return 40; //Zr
               default  : return FALSE;
               }
             }
  default  : return FALSE;                                 //ERROR exit
  }
}
//This function define chemical name from the type
int chemid_to_name(register char id)
{
//Get atoms type
switch (id)
  {
  case  -1 :
  case   0 : return *((int*)"e  ");
  case   1 : return *((int*)"H  ");
  case   2 : return *((int*)"He ");
  case   3 : return *((int*)"Li ");
  case   4 : return *((int*)"Be ");
  case   5 : return *((int*)"B  ");
  case   6 : return *((int*)"C  ");
  case   7 : return *((int*)"N  ");
  case   8 : return *((int*)"O  ");
  case   9 : return *((int*)"F  ");
  case  10 : return *((int*)"Ne ");
  case  11 : return *((int*)"Na ");
  case  12 : return *((int*)"Mg ");
  case  13 : return *((int*)"Al ");
  case  14 : return *((int*)"Si ");
  case  15 : return *((int*)"P  ");
  case  16 : return *((int*)"S  ");
  case  17 : return *((int*)"Cl ");
  case  18 : return *((int*)"Ar ");
  case  19 : return *((int*)"K  ");
  case  20 : return *((int*)"Ca ");
  case  21 : return *((int*)"Sc ");
  case  22 : return *((int*)"Ti ");
  case  23 : return *((int*)"V  ");
  case  24 : return *((int*)"Cr ");
  case  25 : return *((int*)"Mn ");
  case  26 : return *((int*)"Fe ");
  case  27 : return *((int*)"Co ");
  case  28 : return *((int*)"Ni ");
  case  29 : return *((int*)"Cu ");
  case  30 : return *((int*)"Zn ");
  case  31 : return *((int*)"Ga ");
  case  32 : return *((int*)"Ge ");
  case  33 : return *((int*)"As ");
  case  34 : return *((int*)"Se ");
  case  35 : return *((int*)"Br ");
  case  36 : return *((int*)"Kr ");
  case  50 : return *((int*)"Sn ");
  case  51 : return *((int*)"Sb ");
  case  52 : return *((int*)"Te ");
  case  53 : return *((int*)"I  ");
  case  54 : return *((int*)"Xe ");
  default  : return FALSE;                                 //ERROR exit
  }
}

//This function assign unique namme for atoms in str. The names should be formatted as lexems.
char get_unique_atom_name(unsigned int names_size,char (*anames)[sizeof(int)])
{
const char * const letters="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
register unsigned int _i, _j, _k, _l, _n, _m;
unsigned int combinatorial_amount[sizeof(int)];

for (combinatorial_amount[0]=strlen(letters), _i=1; _i<sizeof(int); _i++) combinatorial_amount[_i]=combinatorial_amount[_i-1]*combinatorial_amount[0];
for (_i=1; _i<names_size; _i++)
  {
  _j=0, _k=1; while ( (_k!=sizeof(int))&&(anames[_i][_k]!=' ')&&(anames[_i][_k]!='\t') ) _k++;
  while (find_in_row(*((unsigned int*)(anames[_i])),_i,(unsigned int*)anames)!=(unsigned int)-1)
    {
    EDIT_NAME: if (_j==combinatorial_amount[sizeof(int)-_k]) return FALSE; 
    _n=0; while (_j>=combinatorial_amount[_n]) _n++; //Get mantissa lenght
    _m=_j, _l=0; while (_l!=_n) { anames[_i][_k+_l]=letters[_m/combinatorial_amount[_l]], _m%=combinatorial_amount[_l], _l++; } anames[_i][_k+_l]=letters[_m];
    //Skip names reserved for chemical element
    _j++;
    if ( (!(_n))&&(_k==1)&&( (name_to_chemid(anames[_i]))) ) goto EDIT_NAME;
    }
  }
return TRUE;
}

//ONE DAY TO IMPLEMENT READ_MOL2 WITH yfgets() function
//This function load molecule from mol2 file
t_str *read_mol2(FILE *in)
{
register unsigned int _i, _j, count;
register char *lexem;
unsigned int natoms, nedges, size_r;
t_str *str=0x0;
char *name=0x0;

//Catch structure start
while ( (fgets(buffer,0xFE,in)))
//Upload molecule 
  if (!(strncmp(buffer,"@<TRIPOS>MOLECULE",0x11)))
    {//Upload molecule
    count++;
    if ( (!(fgets(buffer,0xFE,in)))||(!(check_lexem(0x1,buffer))) )
      { 
      lexem="wrong molecule name"; 
      LABEL_DATA_FORMAT_ERROR_0: ylib_errno=YERROR_DATA_FORMAT;
      yprintf(YPRINTF_ERROR,"Wrong format of data (%s) in line %d of the input mol2 file\n",lexem,count);
      return FALSE;
      }
    else _i=lexlen(buffer);
    if (!(name=(char*)malloc(sizeof(char)*(_i+1))))
      {
      LABEL_ERROR_MEMORY: ylib_errno=YERROR_MEMORY;
      yprintf(YPRINTF_ERROR,"There is a memory allocation problem of read_mol() routine.\n",lexem,count);
      return FALSE;
      }
    else { name[_i]='\0', strncpy(name,buffer,_i); }
    count++;
    if ( (!(fgets(buffer,0xFE,in)))||
         (!(lexem=get_lex(0x1,buffer)))||(!(check_int_type(lexem)))||((int)(natoms=atoi(lexem))<=0)||
         (!(lexem=get_lex(0x2,lexem))) ||(!(check_int_type(lexem)))||((int)(nedges=atoi(lexem))< 0)||
         (!(lexem=get_lex(0x2,lexem))) ||(!(check_int_type(lexem)))||((int)(size_r=atoi(lexem))< 0)||
         (!(check_lexem(0x3,lexem)))||( (check_lexem(0x4,lexem))) )
      { free(name); name=0x0; lexem="wrong molecule properties description"; goto LABEL_DATA_FORMAT_ERROR_0; }
    if (!size_r) size_r=1;
    if (!(str=(t_str*)alloc_str(FALSE,size_r,natoms,nedges)))
      { free(name); name=0x0; goto LABEL_ERROR_MEMORY; }
    else str->name=name;
//Upload atoms
    if ( (str->ress)) str->rsize[0]=0;
    while ( (fgets(buffer,0xFE,in)))
      if (!(strncmp(buffer,"@<TRIPOS>ATOM",0xD)))
        {//Upload atoms
        GET_FILE_STRING_END(buffer,in); count++;
        //Upload the first atom 
        if (!(fgets(buffer,0xFE,in)))
          {
          LABEL_IO_ERROR: ylib_errno=YERROR_IO;
          yprintf(YPRINTF_ERROR,"Unexpeced input file truncation at line %1d\n",count);
          free_str(str); free(str); str=0x0;
          return FALSE;
          }
        else count++; 
        if ( (!(lexem=get_lex(0x1,buffer)))||(!(check_int_type(lexem)))||(atoi(lexem)!=1) )   
          { 
          lexem="wrong atom id";
          LABEL_DATA_FORMAT_ERROR_1: free_str(str); free(str); str=0x0;
          goto LABEL_DATA_FORMAT_ERROR_0;
          }
        if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) { lexem="wrong atom name"; goto LABEL_DATA_FORMAT_ERROR_1; }
        { _j=sizeof(int); while (_i!=_j--) ((char*)&str->anames[0])[_j]=' '; do { ((char*)&str->anames[0])[_j]=lexem[_j]; } while (_j--) ; } 
        if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) { lexem="wrong atom x coordinate"; goto LABEL_DATA_FORMAT_ERROR_1; }
        str->r[0].i=atof(lexem);
        if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) { lexem="wrong atom y coordinate"; goto LABEL_DATA_FORMAT_ERROR_1; }
        str->r[0].j=atof(lexem);
        if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) { lexem="wrong atom z coordinate"; goto LABEL_DATA_FORMAT_ERROR_1; }
        str->r[0].k=atof(lexem);
        if ( (!(lexem=get_lex(0x2,lexem)))||(!(str->a[0]=(unsigned char)name_to_chemid(lexem))) ) { lexem="wrong chemical type of atom"; goto LABEL_DATA_FORMAT_ERROR_1; }
        if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||((_i=atoi(lexem))<=0) )   { lexem="wrong residue id"; goto LABEL_DATA_FORMAT_ERROR_1; }
        str->start_rid=_i;
        if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) { lexem="wrong residue name"; goto LABEL_DATA_FORMAT_ERROR_1; }
          {
          if (!str->ress) { _j=sizeof(unsigned int); while (_i!=_j--) ((char*)&str->rsize)[_j]=' '; do { ((char*)&str->rsize)[_j]=lexem[_j]; } while (_j--) ; }
          else { _j=sizeof(unsigned int); while (_i!=_j--) ((char*)&str->ress->list[0])[_j]=' '; do { ((char*)&str->ress->list[0])[_j]=lexem[_j]; } while (_j--) ; } 
          }
        if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) { lexem="wrong charge"; goto LABEL_DATA_FORMAT_ERROR_1; } //just skip 'their' charges :D
        GET_FILE_STRING_END(buffer,in);
        str->natoms=1; 
        //Upload the rest of the atoms
        while (natoms!=str->natoms)
          {
          if (!(fgets(buffer,0xFE,in))) goto LABEL_IO_ERROR;
          count++;
          if ( (!(lexem=get_lex(0x1,buffer)))||(!(check_int_type(lexem)))||(atoi(lexem)!=str->natoms+1) ) { lexem="gap in atom ids"; goto LABEL_DATA_FORMAT_ERROR_1; }
          if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) )  { lexem="wrong atom name"; goto LABEL_DATA_FORMAT_ERROR_1; }
          { _j=sizeof(int); while (_i!=_j--) ((char*)&str->anames[str->natoms])[_j]=' '; do { ((char*)&str->anames[str->natoms])[_j]=lexem[_j]; } while (_j--) ; } 
          if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) { lexem="wrong atom x coordinate"; goto LABEL_DATA_FORMAT_ERROR_1; }
          str->r[str->natoms].i=atof(lexem);
          if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) { lexem="wrong atom y coordinate"; goto LABEL_DATA_FORMAT_ERROR_1; }
          str->r[str->natoms].j=atof(lexem);
          if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) { lexem="wrong atom z coordinate"; goto LABEL_DATA_FORMAT_ERROR_1; }
          str->r[str->natoms].k=atof(lexem);
          if ( (!(lexem=get_lex(0x2,lexem)))||(!(str->a[str->natoms]=(unsigned char)name_to_chemid(lexem))) ) { lexem="wrong chemical type of atom"; goto LABEL_DATA_FORMAT_ERROR_1; }
          if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_int_type(lexem)))||((_i=atoi(lexem))<=0) )    { lexem="wrong residue id"; goto LABEL_DATA_FORMAT_ERROR_1; }
          if (!str->ress)
            {
            if (_i!=str->start_rid) { lexem="residue number is bigger that it was declared in the header"; goto LABEL_DATA_FORMAT_ERROR_1; }
            else
              {// Keep loading the current residue, confirm the residue name
              if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) { lexem="wrong residue name"; goto LABEL_DATA_FORMAT_ERROR_1; } 
              _j=sizeof(unsigned int);
              while (_i!=_j--) if (((char*)&str->rsize)[_j]!=' ') { lexem="residue name doesn't match"; goto LABEL_DATA_FORMAT_ERROR_1; }
              do { if (((char*)&str->rsize)[_j]!=lexem[_j]) { lexem="residue name doesn't match"; goto LABEL_DATA_FORMAT_ERROR_1; } } while (_j--) ;
              }
            }
          else
            {
            if (_i!=str->start_rid+str->ress->size)
              {// A new residue, store the new name
              if (_i!=str->start_rid+str->ress->size+1) { lexem="sequence gaps are detected"; goto LABEL_DATA_FORMAT_ERROR_1; }
              if (_i>size_r) { lexem="residue number is bigger that it was declared in the header"; goto LABEL_DATA_FORMAT_ERROR_1; }
              str->rsize[str->ress->size++]=str->natoms;  
              if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) { lexem="wrong residue name"; goto LABEL_DATA_FORMAT_ERROR_1; }
              { _j=sizeof(unsigned int); while (_i!=_j--) ((char*)&str->ress->list[str->ress->size])[_j]=' '; do { ((char*)&str->ress->list[str->ress->size])[_j]=lexem[_j]; } while (_j--) ; }
              }
            else
              {// Keep loading the current residue, confirm the residue name
              if ( (!(lexem=get_lex(0x2,lexem)))||((_i=lexlen(lexem))>sizeof(unsigned int)) ) { lexem="wrong residue name"; goto LABEL_DATA_FORMAT_ERROR_1; } 
              _j=sizeof(unsigned int);
              while (_i!=_j--) if (((char*)&str->ress->list[str->ress->size])[_j]!=' ') { lexem="residue name doesn't match"; goto LABEL_DATA_FORMAT_ERROR_1; }
              do { if (((char*)&str->ress->list[str->ress->size])[_j]!=lexem[_j]) { lexem="residue name doesn't match"; goto LABEL_DATA_FORMAT_ERROR_1; } } while (_j--) ;
              }
            }
          if ( (!(lexem=get_lex(0x2,lexem)))||(!(check_double_type(lexem))) ) { lexem="wrong charge"; goto LABEL_DATA_FORMAT_ERROR_1; } //just skip 'their' charges :D
          GET_FILE_STRING_END(buffer,in);
          str->natoms++; 
          }
        if ( (str->ress)) str->rsize[str->ress->size++]=str->natoms;
//Upload bonds
        if (!nedges) return str; //molecule may consist of a single atom
        while ((fgets(buffer,0xFE,in)))
          if (!(strncmp(buffer,"@<TRIPOS>BOND",0xD)))
            {//Upload bonds
            GET_FILE_STRING_END(buffer,in); count++;
            while (nedges!=str->nedges)
              { 
              if (!(fgets(buffer,0xFE,in))) goto LABEL_IO_ERROR;
              count++;
              if ( (!(lexem=get_lex(0x1,buffer)))||(!(check_int_type(lexem)))||(atoi(lexem)!=str->nedges+1) ) { lexem="wrong bond id"; goto LABEL_DATA_FORMAT_ERROR_1; }
              if ( (!(lexem=get_lex(0x2,lexem )))||(!(check_int_type(lexem)))||((str->edges[str->nedges].vertice[0]=atoi(get_lex(0x2,buffer)))>str->natoms) ) { lexem="wrong atom #1 index of bond description"; goto LABEL_DATA_FORMAT_ERROR_1; }
              if ( (!(lexem=get_lex(0x2,lexem )))||(!(check_int_type(lexem)))||((str->edges[str->nedges].vertice[1]=atoi(get_lex(0x3,buffer)))>str->natoms) ) { lexem="wrong atom #2 index of bond description"; goto LABEL_DATA_FORMAT_ERROR_1; }
              if (!(lexem=get_lex(0x2,lexem )))  { lexem="wrong bond type description"; goto LABEL_DATA_FORMAT_ERROR_1; }
              switch (lexem[0x0]) //define bond type from lexem type
                {
                case '1' :
                case '2' :
                case '3' : {
                           if (lexlen(lexem)!=1) { lexem="wrong bond type description"; goto LABEL_DATA_FORMAT_ERROR_1; }
                           str->edges[str->nedges].type=(int)lexem[0x0];
                           break;
                           }
                case 'a' : {
                           if (lexlen(lexem)!=2) { lexem="wrong bond type description"; goto LABEL_DATA_FORMAT_ERROR_1; }
                           switch (lexem[0x1])
                             {
                             case 'm' : {//Amide bond 
                                        str->edges[str->nedges].type=(int)('m');
                                        break;
                                        }
                             case 'r' : {//Aromatic bond
                                        str->edges[str->nedges].type=(int)('a');
                                        break;
                                        }
                             default  : { lexem="wrong bond type description"; goto LABEL_DATA_FORMAT_ERROR_1; }
                             }
                           break;
                           }
                default  : { lexem="wrong bond type description"; goto LABEL_DATA_FORMAT_ERROR_1; }
                }
              //Update zero-start enumeration
              str->edges[str->nedges].vertice[0]--;
              str->edges[str->nedges].vertice[1]--;
              str->nedges++;
              GET_FILE_STRING_END(buffer,in);
              }
            //Everything was successful! 
            return str;
            }//End of upload bonds
          else { GET_FILE_STRING_END(buffer,in); count++; }
        }//End of upload atoms
      else { GET_FILE_STRING_END(buffer,in); count++; }
    free_str(str); free(str); str=0x0;
    }//End of upload molecule description
  else { GET_FILE_STRING_END(buffer,in); count++; }

return FALSE;  
}

//This function writes pdb to hdd. It can be applied to dei_onized mol only!
//This function exports mol2 files with charges
char write_mol2(FILE *out,char (*order)[4],t_mol *mol)
{
unsigned int _i, _j, _n;
char ch_type[0x4],mol2_type[0x4];

fprintf(out,"\n@<TRIPOS>MOLECULE\nYLIG\n %4.d %4.d    1    0     0\nSMALL\nYCDS_CHARGES\n\n\n@<TRIPOS>ATOM\n",mol->natoms,mol->nedges);
for (_j=0,_i=0;_i<mol->natoms;_i++)
  {
  _n=order[_i][0]+order[_i][1]+order[_i][2]+order[_i][3];
  switch (mol->a[_i])
    {
    case  1 : { memcpy(mol2_type,"H  ",sizeof(char)*0x4); memcpy(ch_type,"H  ",sizeof(char)*0x4); break; }
    case  6 : {
      switch (_n)
        {
        case 4 : { memcpy(mol2_type,"C.3",sizeof(char)*0x4); break; }
        case 3 : { memcpy(mol2_type,"C.2",sizeof(char)*0x4); break; }
        case 2 :
        case 1 : { memcpy(mol2_type,"C.1",sizeof(char)*0x4); break; }
        default: goto UNKNOWN_TYPE;
        }
      memcpy(ch_type,"C  ",sizeof(char)*0x4); break;
      }
    case  7 : {
      switch (_n)
        {
        case 4 : { memcpy(mol2_type,"N.4",sizeof(char)*0x4); break; }
        case 3 : { memcpy(mol2_type,"N.3",sizeof(char)*0x4); break; }
        case 2 : { memcpy(mol2_type,"N.2",sizeof(char)*0x4); break; }
        case 1 : { memcpy(mol2_type,"N.1",sizeof(char)*0x4); break; }
        default: goto UNKNOWN_TYPE;
        }
      memcpy(ch_type,"N  ",sizeof(char)*0x4); break;
      }
    case  8 : {
      switch (_n)
        {
        case 2 : { memcpy(mol2_type,"O.3",sizeof(char)*0x4); break; }
        case 1 : { memcpy(mol2_type,"O.2",sizeof(char)*0x4); break; }
        default: goto UNKNOWN_TYPE;
        }
      memcpy(ch_type,"O  ",sizeof(char)*0x4); break;
      }
    case  9 : { memcpy(mol2_type,"F  ",sizeof(char)*0x4); memcpy(ch_type,"F  ",sizeof(char)*0x4); break; }
    case 14 : { memcpy(mol2_type,"Si ",sizeof(char)*0x4); memcpy(ch_type,"Si ",sizeof(char)*0x4); break; }
    case 15 : {
      switch (_n)
        {
        case 4 : { memcpy(mol2_type,"P.4",sizeof(char)*0x4); break; }
        case 3 : { memcpy(mol2_type,"P.3",sizeof(char)*0x4); break; }
        case 2 : { memcpy(mol2_type,"P.2",sizeof(char)*0x4); break; }
        case 1 : { memcpy(mol2_type,"P.1",sizeof(char)*0x4); break; }
        default: goto UNKNOWN_TYPE;
        }
      memcpy(ch_type,"P  ",sizeof(char)*0x4); break;
      }
    case 16 : {
      switch (_n)
        {
        case 4 :
        case 3 : { memcpy(mol2_type,"S  ",sizeof(char)*0x4); break; }
        case 2 : { memcpy(mol2_type,"S.3",sizeof(char)*0x4); break; }
        case 1 : { memcpy(mol2_type,"S.2",sizeof(char)*0x4); break; }
        default: goto UNKNOWN_TYPE;
        }
      memcpy(ch_type,"S  ",sizeof(char)*0x4); break;
      }
    case 17 : { memcpy(ch_type,"Cl ",sizeof(char)*0x4); memcpy(mol2_type,"Cl ",sizeof(char)*0x4); break; }
    case 35 : { memcpy(ch_type,"Br ",sizeof(char)*0x4); memcpy(mol2_type,"Br ",sizeof(char)*0x4); break; }
    case 53 : { memcpy(ch_type,"I  ",sizeof(char)*0x4); memcpy(mol2_type,"I  ",sizeof(char)*0x4); break; }
    default : {
        switch (mol->a[_i])
          {
          case  -1 :
          case  -6 :
          case  -7 :
          case  -8 :
          case  -9 :
          case -14 :
          case -15 :
          case -16 :
          case -17 :
          case -35 :
          case -53 : { memcpy(ch_type,"H  ",sizeof(char)*0x4); memcpy(mol2_type,"H  ",sizeof(char)*0x4); break; }
          default  : { UNKNOWN_TYPE: memcpy(ch_type,"X  ",sizeof(char)*0x4); memcpy(mol2_type,"??? ",sizeof(char)*0x4); } //Unknown atom type
          }
      }
    }
  if (_i==mol->rsize[_j+1]) _j++;
  fprintf(out," %5.d %3.3s        %8.4f  %8.4f  %8.4f %3.3s    %4.4s %1d       %8.4f \n",_i+1,(char*)&mol->anames[_i],mol->r[_i].i,mol->r[_i].j,mol->r[_i].k,mol2_type,(char*)&mol->ress->list[_j],mol->start_rid+_j-1,(mol->charges) ? mol->charges[_i] : 0.000 );
  }
fprintf(out,"@<TRIPOS>BOND\n");
for (_i=0;_i<mol->nedges;_i++)
       if (mol->edges[_i].type==(int)'a') fprintf(out," %5.d %5.d %5.d    ar\n",_i+1,mol->edges[_i].vertice[0]+1,mol->edges[_i].vertice[1]+1);
  else if (mol->edges[_i].type==(int)'m') fprintf(out," %5.d %5.d %5.d    am\n",_i+1,mol->edges[_i].vertice[0]+1,mol->edges[_i].vertice[1]+1);
  else   fprintf(out," %5.d %5.d %5.d %c\n",_i+1,mol->edges[_i].vertice[0]+1,mol->edges[_i].vertice[1]+1,mol->edges[_i].type);
fprintf(out,"@<TRIPOS>SUBSTRUCTURE\n   1 XXX           1 TEMP              0 ****  ****      0 ROOT\n");
fprintf(out,"\n");
return TRUE;
}

//This function reads first structure from sdf file. It is not safe theoretically relaying onto fitting of every line of molstructure into 0xFE size char buffer!
t_str *read_sdf(unsigned int *count,FILE *in)
{
register unsigned int _i, _j, _k;
char *string, *name;
t_str *str=0x0;
extern char buffer[];

//Stage I. Read the header (3 lines)
if (!(string=yfgets(in))) { LABEL_IO_ERROR_0: ylib_errno=YERROR_IO; return FALSE; } else (*count)++;
if ((_i=strlen(string))==1) name=0x0; //No mol name 
else
  {
  string[_i]='\0'; 
  if (string!=buffer) name=string;
  else { if (!(name=(char*)malloc(sizeof(char)*(_i+1)))) { ylib_errno=YERROR_MEMORY; return FALSE; } else memcpy(name,string,_i); }
  }
fstrskip(in); fstrskip(in); (*count)+=2;
//Stage II. Read record data string
if ( (!(string=yfgets(in)))||(strlen(string)<6) ) { if ( ( (string))&&(string!=buffer) ) { free(string), string=0x0; } goto LABEL_IO_ERROR; } else (*count)++;
string[6]='\0', _j=atoi(&string[3]); string[3]='\0', _i=atoi(&string[0]); if (string!=buffer) { free(string), string=0x0; }
if (!(_i)) { LABEL_ERROR_DATA_FORMAT_0: yprintf(YPRINTF_ERROR,"Wrong sdf data format, line %d\n",*count); ylib_errno=YERROR_DATA_FORMAT; return FALSE; } 
if (!(str=alloc_str(0x0,1,_i,_j))) { if (name) { free(name), name=0; } return FALSE; }
else { str->name=name, str->ress=0x0, *((unsigned int*)&str->rsize)=*((unsigned int*)"DRG "), str->start_rid=1, str->natoms=_i, str->nedges=_j; } 
//Stage III. Read atoms
for (_k=0;_k<str->natoms;_k++)
  if ( (!(string=yfgets(in)))||(strlen(string)<34) )
    { LABEL_IO_ERROR: if ( ( (string))&&(string!=buffer) ) { free(string), string=0x0; } free_str(str), str=0x0; goto LABEL_IO_ERROR_0; } 
  else
    {
    (*count)++, string[10]=string[20]=string[30]=string[34]='\0';
    if ( (!(check_double_type(&string[01])))||(!(check_double_type(&string[11])))||(!(check_double_type(&string[21])))||(!(str->a[_k]=name_to_chemid(&string[31]))) )
      { LABEL_ERROR_DATA_FORMAT: if ( ( (string))&&(string!=buffer) ) { free(string), string=0x0; } free_str(str), str=0x0; goto LABEL_ERROR_DATA_FORMAT_0; } 
    MAKE_INT_NAME(_i,_j,((char*)str->anames[_k]),((char*)&string[31]));
    str->r[_k].i=atof(&string[01]), str->r[_k].j=atof(&string[11]), str->r[_k].k=atof(&string[21]); if (string!=buffer) { free(string), string=0x0; }
    }
//Stage IV. Read bonds
for (_k=0;_k<str->nedges;_k++)
  if ( (!(string=yfgets(in)))||(strlen(string)<9) ) goto LABEL_IO_ERROR;
  else
    { 
    (*count)++;
    string[9]='\0';
    if (!(name=get_lex(0x1,&string[6])))  goto LABEL_ERROR_DATA_FORMAT;
    switch (name[0]) //define bond type from lexem type
      {
      case '0' : {//Skip the edge
        _k--, str->nedges--;
        continue;
        } 
      case '1' :
      case '2' :
      case '3' : {
                 if (lexlen(name)!=1) goto LABEL_ERROR_DATA_FORMAT;
                 str->edges[_k].type=(int)name[0];
                 break;
                 }
      case 'a' : {
                 if (lexlen(name)!=2) goto LABEL_ERROR_DATA_FORMAT;
                 switch (name[1])
                   {
                   case 'm' : {//Amide bond 
                              str->edges[_k].type=(int)('m');
                              break;
                              }
                   case 'r' : {//Aromatic bond
                              str->edges[_k].type=(int)('a');
                              break;
                              }
                    default  : goto LABEL_ERROR_DATA_FORMAT;
                    }
                 }
      default : goto LABEL_ERROR_DATA_FORMAT;
      }
    string[6]='\0';
    if ( (!(check_uint_type(&string[3])))||((str->edges[_k].vertice[1]=atoi(&string[3]))>str->natoms) ) goto LABEL_ERROR_DATA_FORMAT; else str->edges[_k].vertice[1]--;
    string[3]='\0';
    if ( (!(check_uint_type(&string[0])))||((str->edges[_k].vertice[0]=atoi(&string[0]))>str->natoms) ) goto LABEL_ERROR_DATA_FORMAT; else str->edges[_k].vertice[0]--;
    if (string!=buffer) { free(string); string=0x0; }
    } 
//Stage V. Find end of record. NB! Use name as  temporary vatriable here.
while ( (string=yfgets(in)))
  {
  (*count)++;
  if ( ( (name=get_lex(0x1,string)))&&(lexlen(name)==4)&&( (name[0]=='$')||(name[1]=='$')||(name[2]=='$')||(name[3]=='$') )&&(!(get_lex(0x2,name))) )
    {
    if (string!=buffer) { free(string), string=0x0; }
    break;
    }
  if (string!=buffer) { free(string), string=0x0; }
  }
//Exit
return str;
}

/*
//This function replaces residues in mol before edges compilation. Source rsidues at position r_id is replaced with r_name from db. 
char replace_residue_from_db(unsigned int r_id,char r_name[4],t_mol *mol,t_top *top)
{
register unsigned int _i, _j;
register int _k;
unsigned int item;
void *vp;

//Find target residue
if ((item=find_in_list(*((unsigned int*)r_name),top->ress))==(unsigned int)-1) { ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
else if (item==mol->ress->list[r_id]) return TRUE; 
//Edit atoms
if ((_k=top->res[item].atoms.size-top->res[mol->ress->list[r_id]].atoms.size)>0)
  {//The new residue is bigger
  if (!(vp=(t_list*)realloc(mol->atoms,sizeof(t_list)+sizeof(unsigned int)*(mol->atoms->size+_k)))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
  else mol->atoms=(t_list*)vp;
  if (!(vp=(t_vec*)realloc(mol->r,sizeof(t_vec)*(mol->atoms->size+_k)))) goto LABEL_MEMORY_ERROR;
  else mol->r=(t_vec*)vp;
  //Move coords
  _i=mol->atoms->size;
  while (mol->rsize[r_id+1]!=_i--) 
    { mol->atoms->list[_i+_k]=mol->atoms->list[_i], *((int*)&mol->anames[_i+_k])=*((int*)&mol->anames[_i]), mol->r[_i+_k].i=mol->r[_i].i, mol->r[_i+_k].j=mol->r[_i].j, mol->r[_i+_k].k=mol->r[_i].k; }
  mol->atoms->size+=_k;
  //Sort out atoms and coords
  if (!(vp=(t_vec*)malloc(sizeof(t_vec)*top->res[item].atoms.size))) goto LABEL_MEMORY_ERROR;
  _i=top->res[item].atoms.size;
  while (_i--)
    {//Copy into a buffer
    _j=top->res[mol->ress->list[r_id]].atoms.size; 
    while (_j--) 
      if (top->res[item].atoms.list[_i]==top->res[mol->ress->list[r_id]].atoms.list[_j])
        {
        ((t_vec*)vp)[_i].i=mol->r[_j+mol->rsize[r_id]].i, ((t_vec*)vp)[_i].j=mol->r[_j+mol->rsize[r_id]].j, ((t_vec*)vp)[_i].k=mol->r[_j+mol->rsize[r_id]].k;
        goto NEXT_I_BIGGER;
        }
    ((t_vec*)vp)[_i].i=((t_vec*)vp)[_i].j=((t_vec*)vp)[_i].k=(double)NAN;
    NEXT_I_BIGGER: ;
    }
  //Copy from the buffer
  _i=top->res[item].atoms.size; while (_i--) { mol->atoms->list[_i+mol->rsize[r_id]]=top->res[item].ctypes[_i], *((int*)&mol->anames[_i+mol->rsize[r_id]])=*((int*)&top->res[item].atoms.list[_i]), mol->r[_i+mol->rsize[r_id]].i=((t_vec*)vp)[_i].i, mol->r[_i+mol->rsize[r_id]].j=((t_vec*)vp)[_i].j, mol->r[_i+mol->rsize[r_id]].k=((t_vec*)vp)[_i].k; } 
  free(vp);
  //Moves rsize
  _i=r_id; while (mol->ress->size!=_i) mol->rsize[++_i]+=_k; mol->rsize[++_i]+=_k;
  }
else
  {//The new residue is smaller
  //Sort out atoms and coords
  if (!(vp=(t_vec*)malloc(sizeof(t_vec)*top->res[item].atoms.size))) goto LABEL_MEMORY_ERROR;
  _i=top->res[item].atoms.size;
  while (_i--)
    {//Copy into a buffer
    _j=top->res[mol->ress->list[r_id]].atoms.size;
    while (_j--) 
      if (top->res[item].atoms.list[_i]==top->res[mol->ress->list[r_id]].atoms.list[_j])
        {
        ((t_vec*)vp)[_i].i=mol->r[_j+mol->rsize[r_id]].i, ((t_vec*)vp)[_i].j=mol->r[_j+mol->rsize[r_id]].j, ((t_vec*)vp)[_i].k=mol->r[_j+mol->rsize[r_id]].k;
        goto NEXT_I_SMALLER;
        }
    ((t_vec*)vp)[_i].i=((t_vec*)vp)[_i].j=((t_vec*)vp)[_i].k=(double)NAN;
    NEXT_I_SMALLER: ;
    }
  //Copy from the buffer
  _i=top->res[item].atoms.size; while (_i--) { mol->atoms->list[_i+mol->rsize[r_id]]=top->res[item].ctypes[_i], *((int*)&mol->anames[_i+mol->rsize[r_id]])=*((int*)&top->res[item].atoms.list[_i]), mol->r[_i+mol->rsize[r_id]].i=((t_vec*)vp)[_i].i, mol->r[_i+mol->rsize[r_id]].j=((t_vec*)vp)[_i].j, mol->r[_i+mol->rsize[r_id]].k=((t_vec*)vp)[_i].k; } 
  free(vp);
  //Move coords
  _i=mol->rsize[r_id+1]-1;  
  while (++_i!=mol->atoms->size)
    {
    mol->atoms->list[_i+_k]=mol->atoms->list[_i], *((int*)&mol->anames[_i+_k])=*((int*)&mol->anames[_i]), mol->r[_i+_k].i=mol->r[_i].i, mol->r[_i+_k].j=mol->r[_i].j, mol->r[_i+_k].k=mol->r[_i].k;
    } 
  mol->atoms->size+=_k;
  //Moves rsize
  _i=r_id; while (mol->ress->size!=_i) mol->rsize[++_i]+=_k; mol->rsize[++_i]+=_k;   
  }
//anyhing is cool
mol->ress->list[r_id]=item;
return TRUE;
}


//This function upload structure from pdb file
//NOTE. The rules is that each residue has to have '-' and '+' bonds for linear polimerization (if any) and only one '*' and '^' pair in case of nonlinear polimerization (if any). 
t_mol *read_pdb(unsigned int start_type,unsigned int end_type,FILE *in, t_top *top)
{
unsigned int *nline_res=0x0, nline_size=0, item;
int HE2, HD1, __t;
extern char buffer[];
register unsigned int _i, _j, _k, _l, _p, __p, _q, _t;
void *vp;
t_mol *mol;
register double _r; 
double _rmin;
t_edge *nline_edges=0x0;
t_vec u, _u;

if (!(mol=(t_mol*)calloc(sizeof(t_mol),0x1))) { free(mol); ylib_errno=YERROR_MEMORY; return FALSE; }
if (!(mol->ress=(t_list*)alloc_list(0xFF))) { LABEL_MEMORY_ERROR: free_mol(mol); ylib_errno=YERROR_MEMORY; return FALSE; }
else mol->ress->size=(unsigned int)-1;
if ( (!(mol->atoms=(t_list*)alloc_list(0xFF)))||(!(mol->r=(t_vec*)malloc(sizeof(t_vec)*0xFF)))||
     (!(mol->rsize=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) ) goto LABEL_MEMORY_ERROR;
else mol->atoms->size=mol->rsize[0]=0;
if (!(mol->name=(char*)malloc(sizeof(char)*(strlen("YMOL")+1)))) goto LABEL_MEMORY_ERROR;
else memcpy(mol->name,"YMOL",sizeof(char)*(strlen("YMOL")+1));
//Do main cycle
while (fgets(buffer,0xFE,in))
  {
  CUT_COMMENTS(_i,buffer);
  if ( ( (!(strncmp(buffer,"ATOM  ",0x6)))||(!(strncmp(buffer,"HETATM",0x6))) )&&(strlen(buffer)>=PDB_STR_LEN-1) )
    {//Load new atom
    if ( (!check_int_type(&buffer[22]))||(!check_double_type(&buffer[30]))||(!check_double_type(&buffer[38]))||(!check_double_type(&buffer[46])) ) //Process normal atoms record
      { LABEL_DATA_FORMAT_ERROR: free_mol(mol); ylib_errno=YERROR_DATA_FORMAT; return FALSE; }
    if ((item=find_in_list(*((unsigned int*)&buffer[17]),top->ress))==(unsigned int)-1) goto LABEL_DATA_FORMAT_ERROR; 
    if ((mol->ress->size!=(unsigned int)-1))
      {//load new residue
      if ((_i=(unsigned int)atoi(&buffer[22]))!=nline_size)
        {//New residue detected
        if (++nline_size!=_i) //Is new residue consequtive?
          { 
          fseek(in,(long)-strlen(buffer),SEEK_CUR); //step back and finish uploading
          break; 
          }
        if (!((++mol->ress->size+1)%0xFF))
          {//Realloc residues if necessary
          if (!(vp=(void*)realloc_list(mol->ress,mol->ress->size+0xFF))) goto LABEL_MEMORY_ERROR;
          else mol->ress=(t_list*)vp;
          if (!(vp=(void*)realloc(mol->rsize,sizeof(unsigned int)*(mol->ress->size+0xFF+0x1)))) goto LABEL_MEMORY_ERROR;
          else mol->rsize=(unsigned int*)vp;
          }
        if ((mol->ress->list[mol->ress->size]=find_in_list(*((unsigned int*)&buffer[17]),top->ress))==(unsigned int)-1) goto LABEL_DATA_FORMAT_ERROR; 
        mol->rsize[mol->ress->size+1]=mol->rsize[mol->ress->size]+top->res[mol->ress->list[mol->ress->size]].atoms.size;
        //Realloc atoms if necessary
        if (((mol->atoms->size%0xFF)+top->res[mol->ress->list[mol->ress->size]].atoms.size)>=0xFF)
          {
          _j=((mol->atoms->size+top->res[mol->ress->list[mol->ress->size]].atoms.size)/0xFF+1)*0xFF;
          if (!(vp=(void*)realloc_list(mol->atoms,_j))) goto LABEL_MEMORY_ERROR;
          else mol->atoms=(t_list*)vp;
          if (!(vp=(void*)realloc(mol->r,sizeof(t_vec)*_j))) goto LABEL_MEMORY_ERROR;
          else mol->r=(t_vec*)vp;
          }
        //Set new residues atoms coords to NAN and copy their atoms to its type
        for (_j=0;_j<top->res[mol->ress->list[mol->ress->size]].atoms.size;_j++, mol->atoms->size++)
          { 
          mol->atoms->list[mol->atoms->size]=top->res[mol->ress->list[mol->ress->size]].ctypes[_j];
          mol->r[mol->atoms->size].i=mol->r[mol->atoms->size].j=mol->r[mol->atoms->size].k=(double)NAN;
          } 
        }
      }
    else  
      { //Handle start
      mol->start_rid=nline_size=(int)atoi(&buffer[22]);
      mol->ress->size=0, *mol->ress->list=item;
      mol->rsize[0]=0, mol->rsize[1]=top->res[item].atoms.size;
      _j=top->res[item].nedges;
      while (_j--)
        if (((char*)&top->res[*mol->ress->list].edges[_j].vertice[0])[0]=='-')
          {
          if ( (!start_type)||((item=find_in_list(start_type,top->ress))==(unsigned int)-1) )
            { //Lock it with protons
            START_LOCK_WITH_PROTONS:
            if ( (top->res[*mol->ress->list].edges[_j].type!=1)&&(top->res[*mol->ress->list].edges[_j].type!=2)&&(top->res[*mol->ress->list].edges[_j].type!=3) )
              goto LABEL_DATA_CONSISTMENT_ERROR;
            mol->ress->list[1]=mol->ress->list[0], mol->ress->list[0]=(unsigned int)-1;
            mol->rsize[2]=mol->rsize[1]+top->res[*mol->ress->list].edges[_j].type, mol->rsize[1]=top->res[*mol->ress->list].edges[_j].type; 
            mol->atoms->size=mol->rsize[2];
            if (mol->atoms->size>=0xFF)
              {
              _l=(mol->atoms->size/0xFF+1)*0xFF;
              if (!(vp=(void*)realloc_list(mol->atoms,_l))) goto LABEL_MEMORY_ERROR;
              else mol->atoms=(t_list*)vp;
              if (!(vp=(void*)realloc(mol->r,sizeof(t_vec)*_l))) goto LABEL_MEMORY_ERROR;
              else mol->r=(t_vec*)vp;
              }
            _l=top->res[*mol->ress->list].edges[_j].type; while(_l--) { mol->r[_l].i=mol->r[_l].j=mol->r[_l].k=(double)NAN; mol->atoms->list[_j]=1; }
            }
          else 
            { //Lock it with terminal residue
            _l=top->res[item].atoms.size;
            while (_l--) if (!(strncmp(&((char*)&top->res[*mol->ress->list].edges[_j].vertice[0])[1],(char*)&top->res[item].atoms.list[_l],3))) goto PROCESS_FIRST_RESIDUE;
            goto START_LOCK_WITH_PROTONS;
            PROCESS_FIRST_RESIDUE:
            mol->ress->list[1]=mol->ress->list[0], mol->ress->list[0]=item;
            mol->rsize[2]=mol->rsize[1]+top->res[item].atoms.size, mol->rsize[1]=top->res[item].atoms.size; 
            mol->atoms->size=mol->rsize[2]; 
            if (mol->atoms->size>=0xFF)
              {
              _l=(mol->atoms->size/0xFF+1)*0xFF;
              if (!(vp=(void*)realloc_list(mol->atoms,_l))) goto LABEL_MEMORY_ERROR;
              else mol->atoms=(t_list*)vp;
              if (!(vp=(void*)realloc(mol->r,sizeof(t_vec)*_l))) goto LABEL_MEMORY_ERROR;
              else mol->r=(t_vec*)vp;
              }
            _l=top->res[item].atoms.size; while(_l--) { mol->atoms->list[_l]=top->res[item].ctypes[_l]; mol->r[_l].i=mol->r[_l].j=mol->r[_l].k=(double)NAN; }
            }
          _j=top->res[mol->ress->list[1]].atoms.size;
          while (_j--)
            {
            mol->atoms->list[_j+top->res[item].atoms.size]=top->res[mol->ress->list[1]].ctypes[_j];
            mol->r[_j+top->res[item].atoms.size].i=mol->r[_j+top->res[item].atoms.size].j=mol->r[_j+top->res[item].atoms.size].k=(double)NAN;
            }
          mol->ress->size=1;
          goto END_FIRST_PROCESSING;
          }
      //No linear polimerization required
      _j=top->res[*mol->ress->list].atoms.size;
      while (_j--)
        {
        mol->atoms->list[_j]=top->res[mol->ress->list[1]].ctypes[_j];
        mol->r[_j].i=mol->r[_j].j=mol->r[_j].k=(double)NAN;
        }
      END_FIRST_PROCESSING: ;
      }
    //Prepare name
    _j=0;
    while ( (_j<sizeof(unsigned int))&&(buffer[12+_j]==' ') ) _j++;
    if (_j)
      {
      if (_j==sizeof(unsigned int)) continue; //No name given at all but hopely file damaged only slightly
      for (_i=0;_i<sizeof(unsigned int)-_j;_i++) buffer[12+_i]=buffer[12+_i+_j];
      while (_i<sizeof(unsigned int)) buffer[12+_i++]=' ';
      }
    if ((item=find_in_list(*((unsigned int*)&buffer[12]),&top->res[mol->ress->list[mol->ress->size]].atoms))==(unsigned int)-1)  goto LABEL_DATA_CONSISTMENT_ERROR; //Store atoms data into its place
    mol->r[mol->rsize[mol->ress->size]+item].i=atof(&buffer[30]);
    mol->r[mol->rsize[mol->ress->size]+item].j=atof(&buffer[38]);
    mol->r[mol->rsize[mol->ress->size]+item].k=atof(&buffer[46]);
    }
  else if ( (!(strncmp(buffer,"TER",sizeof(char)*0x3)))||(!(strncmp(buffer,"END",sizeof(char)*0x3))) ) break;
  }

//Process terminal residue
if  (mol->ress->size==(unsigned int)-1) { LABEL_DATA_CONSISTMENT_ERROR: ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
_j=top->res[mol->ress->list[mol->ress->size]].nedges;
while (_j--)
  if (((char*)&top->res[mol->ress->list[mol->ress->size]].edges[_j].vertice[0])[0]=='+')
    {//Synchronize/realloc residues
    mol->ress->size++;
    if (!(vp=(void*)realloc_list(mol->ress,mol->ress->size+1))) goto LABEL_MEMORY_ERROR;
    else mol->ress=(t_list*)vp;
    if (!(vp=(void*)realloc(mol->rsize,sizeof(unsigned int)*(mol->ress->size+2)))) goto LABEL_MEMORY_ERROR;
    else mol->rsize=(unsigned int*)vp;
    if ( (end_type)&&((item=find_in_list(end_type,top->ress))!=(unsigned int)-1) )
      {//Set terminal item
      _j=top->res[item].nedges; 
      while (_j--)
        if (((char*)&top->res[item].edges[_j].vertice[0])[0]=='-')
          {
          _l=top->res[mol->ress->list[mol->ress->size-1]].atoms.size;
          while (_l--)
            if (!(strncmp((char*)&top->res[mol->ress->list[mol->ress->size-1]].atoms.list[_l],&((char*)&top->res[item].edges[_j].vertice[0])[1],3))) goto PROCESS_END_RESIDUE;
          }
      goto END_LOCK_WITH_PROTONS;
      PROCESS_END_RESIDUE:
      if (((mol->atoms->size%0xFF)+top->res[item].atoms.size)>=0xFF)
        {//Realloc atoms is necessary
        _j=((mol->atoms->size+top->res[item].atoms.size)/0xFF+1)*0xFF;
        if (!(vp=(void*)realloc_list(mol->atoms,_j))) goto LABEL_MEMORY_ERROR;
        else mol->atoms=(t_list*)vp;
        if (!(vp=(void*)realloc(mol->r,sizeof(t_vec)*_j))) goto LABEL_MEMORY_ERROR;
        else mol->r=(t_vec*)vp;
        }
      mol->ress->list[mol->ress->size]=item;
      mol->rsize[mol->ress->size+1]=mol->atoms->size+=top->res[item].atoms.size;
      _j=top->res[item].atoms.size;
      while (_j--) { _k=mol->rsize[mol->ress->size]+_j; mol->atoms->list[_k]=top->res[item].ctypes[_j]; mol->r[_k].i=mol->r[_k].j=mol->r[_k].k=(double)NAN; }
      }
    else
      {
      END_LOCK_WITH_PROTONS: 
      if ( (top->res[mol->ress->size-1].edges[_j].type!=1)&&(top->res[mol->ress->size-1].edges[_j].type!=2)&&(top->res[mol->ress->size-1].edges[_j].type!=3) ) goto LABEL_DATA_CONSISTMENT_ERROR;
      if (((mol->atoms->size%0xFF)+top->res[*mol->ress->list].edges[_j].type)>=0xFF)
        {
        _i=mol->atoms->size+top->res[*mol->ress->list].edges[_i].type;
        if (!(vp=(void*)realloc_list(mol->atoms,_i))) goto LABEL_MEMORY_ERROR;
        else mol->atoms=(t_list*)vp;
        if (!(vp=(void*)realloc(mol->r,sizeof(t_vec)*_i))) goto LABEL_MEMORY_ERROR;
        else mol->r=(t_vec*)vp;
        }
      mol->ress->list[mol->ress->size]=(unsigned int)-1;
      mol->rsize[mol->ress->size+1]=mol->atoms->size+=top->res[*mol->ress->list].edges[_j].type;
      for (_j=mol->rsize[mol->ress->size];_j<mol->rsize[mol->ress->size+1];_j++) { mol->atoms->list[_j]=1; mol->r[_j].i=mol->r[_j].j=mol->r[_j].k=(double)NAN; }
      }
    mol->ress->size++;
    goto END_PROCESSED;
    }
//Synchronize residues
if (!(vp=(void*)realloc_list(mol->ress,mol->ress->size))) goto LABEL_MEMORY_ERROR;
else mol->ress=(t_list*)vp;
if (!(vp=(void*)realloc(mol->rsize,sizeof(unsigned int)*(mol->ress->size+1)))) goto LABEL_MEMORY_ERROR;
else mol->rsize=(unsigned int*)vp;
END_PROCESSED: ;

//Fill atom names
if (!(mol->anames=(char (*)[sizeof(int)])malloc(sizeof(int)*mol->atoms->size))) goto LABEL_MEMORY_ERROR;
_i=mol->ress->size; 
while (_i--) { _k=top->res[mol->ress->list[_i]].atoms.size; while (_k--) { *((int*)&mol->anames[_k+mol->rsize[_i]])=*((int*)&top->res[mol->ress->list[_i]].atoms.list[_k]); } }

//Determine HIS protonation considering nearest h-bonding.
_i=mol->ress->size;
while (_i--)
  if ( (!(strncmp("HID ",(char*)&top->ress->list[mol->ress->list[_i]],strlen("HID "))))||
       (!(strncmp("HIE ",(char*)&top->ress->list[mol->ress->list[_i]],strlen("HIE "))))||
       (!(strncmp("HIP ",(char*)&top->ress->list[mol->ress->list[_i]],strlen("HIP "))))||
       (!(strncmp("HIS ",(char*)&top->ress->list[mol->ress->list[_i]],strlen("HIS ")))) )
    {
    HD1=HE2=0;
    //Check if NE2 require protonation
    for (_j=(unsigned int)-1, _l=mol->rsize[_i];_l<mol->rsize[_i+1];_l++) //find NE2
      if (!strncmp((char*)&top->res[mol->ress->list[_i]].atoms.list[_l-mol->rsize[_i]],"NE2 ",strlen("NE2 "))) { _j=_l; break; }
    if ( (_j==(unsigned int)-1)||( (isnan(mol->r[_j].i)))||( (isnan(mol->r[_j].j)))||( (isnan(mol->r[_j].k))) ) goto END_OF_HE2;
    else { u.i=mol->r[_j].i, u.j=mol->r[_j].j, u.k=mol->r[_j].k; }
    for (_k=(unsigned int)-1, _l=mol->rsize[_i];_l<mol->rsize[_i+1];_l++) //find CE1
      if (!strncmp((char*)&top->res[mol->ress->list[_i]].atoms.list[_l-mol->rsize[_i]],"CE1 ",strlen("CE1 "))) { _k=_l; break; }
    if ( (_k==(unsigned int)-1)||( (isnan(mol->r[_k].i)))||( (isnan(mol->r[_k].j)))||( (isnan(mol->r[_k].k))) ) goto END_OF_HE2;
    else { u.i-=mol->r[_k].i/2., u.j-=mol->r[_k].j/2., u.k-=mol->r[_k].k/2.; }
    for (_k=(unsigned int)-1, _l=mol->rsize[_i];_l<mol->rsize[_i+1];_l++) //find CE1
      if (!strncmp((char*)&top->res[mol->ress->list[_i]].atoms.list[_l-mol->rsize[_i]],"CD2 ",strlen("CD2 "))) { _k=_l; break; }
    if ( (_k==(unsigned int)-1)||( (isnan(mol->r[_k].i)))||( (isnan(mol->r[_k].j)))||( (isnan(mol->r[_k].k))) ) goto END_OF_HE2;
    else { u.i-=mol->r[_k].i/2., u.j-=mol->r[_k].j/2., u.k-=mol->r[_k].k/2.; multiple_vec_scalar(&u,&u,1./sqrt(calc_vec_norm(&u))); } 
    _t=0, _l=mol->atoms->size;
    while (_l--)
      if ( ( (mol->atoms->list[_l]==7)||(mol->atoms->list[_l]==8)||(mol->atoms->list[_l]==9)||(mol->atoms->list[_l]==16) ) &&
           (!(isnan(mol->r[_l].i)))&&(!(isnan(mol->r[_l].j)))&&(!(isnan(mol->r[_l].k)) )&&((_r=calc_distance(&mol->r[_j],&mol->r[_l]))<MAX_HB_LENGHT2)&&
           ( (u.i*(mol->r[_l].i-mol->r[_j].i)+u.j*(mol->r[_l].j-mol->r[_j].j)+u.k*(mol->r[_l].k-mol->r[_j].k))/sqrt(_r) > MIN_HB_COS ) ) 
        {//Scann nearby heavy atoms
        //Check H-X bond rotability of the atom: it should be in form of (H)n-X-Y, where X and Y are heavy atoms (the search is done inside its residue only and polymerizing atoms are excluded) 
        _p=0; while (mol->rsize[_p+1]<_l) _p++;
        _k=(unsigned int)-1, _q=top->res[mol->ress->list[_p]].nedges;
        while (_q--) 
               if (top->res[mol->ress->list[_p]].edges[_q].vertice[0]==_l-mol->rsize[_p])
                 {
                 if ( (*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[1])=='+')||(*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[1])=='-')||
                      (*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[1])=='*')||(*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[1])=='/')||
                      ( (_k!=(unsigned int)-1)&&(top->res[mol->ress->list[_p]].ctypes[top->res[mol->ress->list[_p]].edges[_q].vertice[1]]!=1) ) ) _k=(unsigned int)-2; 
                 else
                   {
                   if (top->res[mol->ress->list[_p]].ctypes[top->res[mol->ress->list[_p]].edges[_q].vertice[1]]==1) _t=TRUE;
                   else _k=top->res[mol->ress->list[_p]].edges[_q].vertice[1]+mol->rsize[_p];
                   }
                 }
          else if (top->res[mol->ress->list[_p]].edges[_q].vertice[1]==_l-mol->rsize[_p])
                 {
                 if ( (*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[0])=='+')||(*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[0])=='-')||
                      (*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[0])=='*')||(*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[0])=='/')||
                      ( (_k!=(unsigned int)-1)&&(top->res[mol->ress->list[_p]].ctypes[top->res[mol->ress->list[_p]].edges[_q].vertice[0]]!=1) ) ) _k=(unsigned int)-2;
                 else
                   { 
                   if (top->res[mol->ress->list[_p]].ctypes[top->res[mol->ress->list[_p]].edges[_q].vertice[0]]==1) _t=TRUE;
                   else _k=top->res[mol->ress->list[_p]].edges[_q].vertice[1]+mol->rsize[_p];
                   }
                 }
        if ( (_k!=(unsigned int)-1)&&(_k!=(unsigned int)-2)&&(_t) )
          {//The bond is rotable Y-X-H(n) configuration, scann its surrounding
          _u.i=mol->r[_l].i-mol->r[_k].i, _u.j=mol->r[_l].j-mol->r[_k].j, _u.k=mol->r[_l].k-mol->r[_k].k;
          multiple_vec_scalar(&_u,&_u,1./sqrt(calc_vec_norm(&_u)));
          __t=0, _k=mol->atoms->size;
          while (_k--)
            if ( ( (mol->atoms->list[_k]==7)||(mol->atoms->list[_k]==8)||(mol->atoms->list[_k]==9)||(mol->atoms->list[_k]==16) )&&
                 (!(isnan(mol->r[_k].i)))&&(!(isnan(mol->r[_k].j)))&&(!(isnan(mol->r[_k].k)) )&&((_r=calc_distance(&mol->r[_l],&mol->r[_k]))<MAX_HB_LENGHT2)&&
                 ( (u.i*(mol->r[_k].i-mol->r[_l].i)+u.j*(mol->r[_k].j-mol->r[_l].j)+u.k*(mol->r[_k].k-mol->r[_l].k))/sqrt(_r) > MIN_HB_COS ) )
              {
              _t=FALSE, __p=0; while (mol->rsize[__p+1]<_k) __p++;
              _q=top->res[mol->ress->list[__p]].nedges;
              while (_q--) 
                      if ( (*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='+')&&(*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='-')&&
                           (*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='*')&&(*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='/')&&
                           (top->res[mol->ress->list[__p]].edges[_q].vertice[0]==_k-mol->rsize[__p])&&(top->res[mol->ress->list[__p]].ctypes[top->res[mol->ress->list[__p]].edges[_q].vertice[0]]==1) ) _t=TRUE;
                 else if ( (*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[1])!='+')&&(*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[1])!='-')&&
                           (*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[1])!='*')&&(*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[1])!='/')&&
                           (top->res[mol->ress->list[__p]].edges[_q].vertice[1]==_k-mol->rsize[__p])&&(top->res[mol->ress->list[__p]].ctypes[top->res[mol->ress->list[__p]].edges[_q].vertice[1]]==1) ) _t=TRUE;
              if (_t) __t++; else __t--;  
              }
          if (__t>0) HE2++; else HE2--;
          }
        else
          { //The bond is not rotable, summ it directly
          if (_t) HE2++; //Potential acceptor
          else    HE2--; //Potential donor
          }
        }
    END_OF_HE2: ;
    //Check if HD1 requires protonation
    for (_j=(unsigned int)-1, _l=mol->rsize[_i];_l<mol->rsize[_i+1];_l++) //find ND1
      if (!strncmp((char*)&top->res[mol->ress->list[_i]].atoms.list[_l-mol->rsize[_i]],"ND1 ",strlen("ND1 "))) { _j=_l; break; }
    if ( (_j==(unsigned int)-1)||( (isnan(mol->r[_j].i)))||( (isnan(mol->r[_j].j)))||( (isnan(mol->r[_j].k))) ) goto END_OF_HD1;
    else { u.i=mol->r[_j].i, u.j=mol->r[_j].j, u.k=mol->r[_j].k; }
    for (_k=(unsigned int)-1, _l=mol->rsize[_i];_l<mol->rsize[_i+1];_l++) //find CE1
      if (!strncmp((char*)&top->res[mol->ress->list[_i]].atoms.list[_l-mol->rsize[_i]],"CE1 ",strlen("CE1 "))) { _k=_l; break; }
    if ( (_k==(unsigned int)-1)||( (isnan(mol->r[_k].i)))||( (isnan(mol->r[_k].j)))||( (isnan(mol->r[_k].k))) ) goto END_OF_HD1;
    else { u.i-=mol->r[_k].i/2., u.j-=mol->r[_k].j/2., u.k-=mol->r[_k].k/2.; }
    for (_k=(unsigned int)-1, _l=mol->rsize[_i];_l<mol->rsize[_i+1];_l++) //find CG
      if (!strncmp((char*)&top->res[mol->ress->list[_i]].atoms.list[_l-mol->rsize[_i]],"CG  ",strlen("CG  "))) { _k=_l; break; }
    if ( (_k==(unsigned int)-1)||( (isnan(mol->r[_k].i)))||( (isnan(mol->r[_k].j)))||( (isnan(mol->r[_k].k))) ) goto END_OF_HD1;
    else { u.i-=mol->r[_k].i/2., u.j-=mol->r[_k].j/2., u.k-=mol->r[_k].k/2.; multiple_vec_scalar(&u,&u,1./sqrt(calc_vec_norm(&u))); } 
    _l=mol->atoms->size;
    while (_l--)
      if ( ( (mol->atoms->list[_l]==7)||(mol->atoms->list[_l]==8)||(mol->atoms->list[_l]==9)||(mol->atoms->list[_l]==16) ) &&
           (!(isnan(mol->r[_l].i)))&&(!(isnan(mol->r[_l].j)))&&(!(isnan(mol->r[_l].k)) )&&((_r=calc_distance(&mol->r[_j],&mol->r[_l]))<MAX_HB_LENGHT2)&&
           ( (u.i*(mol->r[_l].i-mol->r[_j].i)+u.j*(mol->r[_l].j-mol->r[_j].j)+u.k*(mol->r[_l].k-mol->r[_j].k))/sqrt(_r) > MIN_HB_COS ) ) 
        {//Scann nearby heavy atoms
        //Check H-X bond rotability of the atom: it should be in form of (H)n-X-Y, where X and Y are heavy atoms (the search is done inside its residue only and polymerizing atoms are excluded) 
        _t=FALSE, _p=0; while (mol->rsize[_p+1]<_l) _p++;
        _k=(unsigned int)-1, _q=top->res[mol->ress->list[_p]].nedges;
        while (_q--) 
               if (top->res[mol->ress->list[_p]].edges[_q].vertice[0]==_l-mol->rsize[_p])
                 {
                 if ( (*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[1])=='+')||(*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[1])=='-')||
                      (*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[1])=='*')||(*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[1])=='/')||
                      ( (_k!=(unsigned int)-1)&&(top->res[mol->ress->list[_p]].ctypes[top->res[mol->ress->list[_p]].edges[_q].vertice[1]]!=1) ) ) _k=(unsigned int)-2; 
                 else 
                   {
                   if (top->res[mol->ress->list[_p]].ctypes[top->res[mol->ress->list[_p]].edges[_q].vertice[1]]==1) _t=TRUE;
                   else _k=top->res[mol->ress->list[_p]].edges[_q].vertice[1]+mol->rsize[_p];
                   }
                 }
          else if (top->res[mol->ress->list[_p]].edges[_q].vertice[1]==_l-mol->rsize[_p])
                 {
                 if ( (*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[0])=='+')||(*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[0])=='-')||
                      (*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[0])=='*')||(*((char*)&top->res[mol->ress->list[_p]].edges[_q].vertice[0])=='/')||
                      ( (_k!=(unsigned int)-1)&&(top->res[mol->ress->list[_p]].ctypes[top->res[mol->ress->list[_p]].edges[_q].vertice[0]]!=1) ) ) _k=(unsigned int)-2; 
                 else 
                   {
                   if (top->res[mol->ress->list[_p]].ctypes[top->res[mol->ress->list[_p]].edges[_q].vertice[0]]==1) _t=TRUE;
                   _k=top->res[mol->ress->list[_p]].edges[_q].vertice[1]+mol->rsize[_p];
                   }
                 }
        if ( (_k!=(unsigned int)-1)&&(_k!=(unsigned int)-2)&&(_t) )
          {//The bond is rotable X-A-H(n) configuration, scann its surrounding
          _u.i=mol->r[_l].i-mol->r[_k].i, _u.j=mol->r[_l].j-mol->r[_k].j, _u.k=mol->r[_l].k-mol->r[_k].k;
          multiple_vec_scalar(&_u,&_u,1./sqrt(calc_vec_norm(&_u)));
          __t=0, _k=mol->atoms->size;
          while (_k--)
            if ( ( (mol->atoms->list[_k]==7)||(mol->atoms->list[_k]==8)||(mol->atoms->list[_k]==9)||(mol->atoms->list[_k]==16) ) &&
                 (!(isnan(mol->r[_k].i)))&&(!(isnan(mol->r[_k].j)))&&(!(isnan(mol->r[_k].k)) )&&((_r=calc_distance(&mol->r[_l],&mol->r[_k]))<MAX_HB_LENGHT2) &&
                 ( (u.i*(mol->r[_k].i-mol->r[_l].i)+u.j*(mol->r[_k].j-mol->r[_l].j)+u.k*(mol->r[_k].k-mol->r[_l].k))/sqrt(_r) > MIN_HB_COS ) )
              {
              _t=FALSE, __p=0; while (mol->rsize[__p+1]<_k) __p++;
              _q=top->res[mol->ress->list[__p]].nedges;
              while (_q--) 
                      if ( (*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='+')&&(*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='-')&&
                           (*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='*')&&(*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='/')&&
                           (top->res[mol->ress->list[__p]].edges[_q].vertice[0]==_k-mol->rsize[__p])&&(top->res[mol->ress->list[__p]].ctypes[top->res[mol->ress->list[__p]].edges[_q].vertice[0]]==1) ) _t=TRUE;
                 else if ( (*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[1])!='+')&&(*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='-')&&
                           (*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[1])!='*')&&(*((char*)&top->res[mol->ress->list[__p]].edges[_q].vertice[0])!='/')&&
                           (top->res[mol->ress->list[__p]].edges[_q].vertice[1]==_k-mol->rsize[__p])&&(top->res[mol->ress->list[__p]].ctypes[top->res[mol->ress->list[__p]].edges[_q].vertice[1]]==1) ) _t=TRUE;
              if (_t) __t++; else __t--;  
              }
          if (__t>0) HD1++; else HD1--;
          }
        else
          { //The bond is not rotable, summ it directly
          if (_t) HD1++; //Potential acceptor
          else    HD1--; //Potential donor
          }
        }
    END_OF_HD1: ;
    //Select and replace the residue
    if (HD1<0)
      {
      if (HE2<0)
        { if (!(replace_residue_from_db(_i,"HIP ",mol,top))) return FALSE; }
      else
        { if (!(replace_residue_from_db(_i,"HID ",mol,top))) return FALSE; }
      }
    else
      {
      if (HE2<0)
        { if (!(replace_residue_from_db(_i,"HIE ",mol,top))) return FALSE; }
      else
        {//unprotonated 'HIS' is not allowed!
        if (HD1<HE2) { if (!(replace_residue_from_db(_i,"HID ",mol,top))) return FALSE; }
        else         { if (!(replace_residue_from_db(_i,"HIE ",mol,top))) return FALSE; }
        }
      }
    }

//Get bonds number and nonlinear polimerizatioin matrix
nline_size=_k=0;
for (_j=0;_j<mol->ress->size;_j++)
  if (mol->ress->list[_j]!=(unsigned int)-1)
    {//Setup bonds
    _i=top->res[mol->ress->list[_j]].nedges;
    while (_i--)
      if (top->res[mol->ress->list[_j]].edges[_i].type<0)
        { //Nonlinear polimerization detected
             if (((char*)&top->res[mol->ress->list[_j]].edges[_i].vertice[0])[0]=='*')
               {//Only geomety is used to restore nonlinear bonding
               if ( (mol->r[mol->rsize[_j]+top->res[mol->ress->list[_j]].edges[_i].vertice[1]].i==(double)NAN)||
                    (mol->r[mol->rsize[_j]+top->res[mol->ress->list[_j]].edges[_i].vertice[1]].j==(double)NAN)||
                    (mol->r[mol->rsize[_j]+top->res[mol->ress->list[_j]].edges[_i].vertice[1]].k==(double)NAN) ) goto LABEL_DATA_CONSISTMENT_ERROR;
               if (!(nline_size%0xFF))
                 {
                 if (!(vp=realloc(nline_edges,sizeof(t_edge)*(nline_size+0xFF))))     { free(nline_res); free(nline_edges); goto LABEL_MEMORY_ERROR; }
                 else nline_edges=(t_edge*)vp;
                 if (!(vp=realloc(nline_res,sizeof(unsigned int)*(nline_size+0xFF)))) { free(nline_res); free(nline_edges); goto LABEL_MEMORY_ERROR; }
                 else nline_res=(unsigned int*)vp;
                 }
               nline_res[nline_size]=_j;
               nline_edges[nline_size].vertice[0]=top->res[mol->ress->list[_j]].edges[_i].vertice[0];
               nline_edges[nline_size].vertice[1]=top->res[mol->ress->list[_j]].edges[_i].vertice[1];
               nline_edges[nline_size].type=-top->res[mol->ress->list[_j]].edges[_i].type;
               ((char*)&nline_edges[nline_size].vertice[0])[0]=((char*)&nline_edges[nline_size].vertice[0])[1];
               ((char*)&nline_edges[nline_size].vertice[0])[1]=((char*)&nline_edges[nline_size].vertice[0])[2];
               ((char*)&nline_edges[nline_size].vertice[0])[2]=((char*)&nline_edges[nline_size].vertice[0])[3];
               ((char*)&nline_edges[nline_size].vertice[0])[3]=' ';
               _k+=nline_edges[nline_size].type;
               nline_size++;
               }
        else if (((char*)&top->res[mol->ress->list[_j]].edges[_i].vertice[0])[0]=='+')
               {
               if (mol->ress->list[_j+1]==(unsigned int)-1) _k+=mol->rsize[_j+1]-mol->rsize[_j];
               }
        else   {
               if (mol->ress->list[_j-1]==(unsigned int)-1) _k+=mol->rsize[_j-1]-mol->rsize[_j];
               else _k++;
               }
        }
      else _k++;
    }
if (!(mol->edges=(t_edge*)malloc(sizeof(t_edge)*_k))) { free(nline_edges); goto LABEL_MEMORY_ERROR; } // _k - is upper limit of bonds in the system

//Compile bonds
for (_j=0;_j<mol->ress->size;_j++)
  if (mol->ress->list[_j]!=(unsigned int)-1)
    {
    for (_i=0;_i<top->res[mol->ress->list[_j]].nedges;_i++)
      if (top->res[mol->ress->list[_j]].edges[_i].type<0)
        { //Polymerize nonlineary
             if (((char*)&top->res[mol->ress->list[_j]].edges[_i].vertice[0])[0]=='*') continue; //Nonlinear is already stored separately
        //Polymerize lineary
        else if (((char*)&top->res[mol->ress->list[_j]].edges[_i].vertice[0])[0]=='+')
               {//if residue is nonterminal with proton locks than check if it has required atoms
               if (_j==mol->ress->size) goto LABEL_DATA_CONSISTMENT_ERROR;
               if (mol->ress->list[_j+1]==(unsigned int)-1) 
                 {//Add all bonds 
                 for (_l=mol->rsize[_j+1];_l<mol->rsize[_j+2];_l++, mol->nedges++)
                   {
                   mol->edges[mol->nedges].vertice[0]=_l;
                   mol->edges[mol->nedges].vertice[1]=mol->rsize[_j]+top->res[mol->ress->list[_j]].edges[_i].vertice[1];
                   mol->edges[mol->nedges].type='1';
                   }  
                 }
               else
                 {//Check if its '-' correspons to given '+'
                 for (_l=0;_l<top->res[mol->ress->list[_j+1]].nedges;_l++) if (((char*)&top->res[mol->ress->list[_j+1]].edges[_l].vertice[0])[0]=='-') goto CHECK_PLUS_MINUS_BOND;
                 goto LABEL_DATA_CONSISTMENT_ERROR;
                 CHECK_PLUS_MINUS_BOND: //Check names and type correspondence...
                 if ( (((char*)&top->res[mol->ress->list[_j]].edges[_i].vertice[0])[1]==((char*)&top->res[mol->ress->list[_j+1]].atoms.list[top->res[mol->ress->list[_j+1]].edges[_l].vertice[1]])[0])&&
                      (((char*)&top->res[mol->ress->list[_j]].edges[_i].vertice[0])[2]==((char*)&top->res[mol->ress->list[_j+1]].atoms.list[top->res[mol->ress->list[_j+1]].edges[_l].vertice[1]])[1])&&
                      (((char*)&top->res[mol->ress->list[_j]].edges[_i].vertice[0])[3]==((char*)&top->res[mol->ress->list[_j+1]].atoms.list[top->res[mol->ress->list[_j+1]].edges[_l].vertice[1]])[2])&&
                      (((char*)&top->res[mol->ress->list[_j+1]].atoms.list[top->res[mol->ress->list[_j+1]].edges[_l].vertice[1]])[2]==' ')&&
                      (top->res[mol->ress->list[_j]].edges[_i].type==top->res[mol->ress->list[_j+1]].edges[_l].type)                      &&
                      (((char*)&top->res[mol->ress->list[_j+1]].edges[_l].vertice[0])[1]==((char*)&top->res[mol->ress->list[_j]].atoms.list[top->res[mol->ress->list[_j]].edges[_i].vertice[1]])[0])&&
                      (((char*)&top->res[mol->ress->list[_j+1]].edges[_l].vertice[0])[2]==((char*)&top->res[mol->ress->list[_j]].atoms.list[top->res[mol->ress->list[_j]].edges[_i].vertice[1]])[1])&&
                      (((char*)&top->res[mol->ress->list[_j+1]].edges[_l].vertice[0])[3]==((char*)&top->res[mol->ress->list[_j]].atoms.list[top->res[mol->ress->list[_j]].edges[_i].vertice[1]])[2])&&
                      (((char*)&top->res[mol->ress->list[_j]].atoms.list[top->res[mol->ress->list[_j]].edges[_i].vertice[1]])[2]==' ') )
                   {
                   mol->edges[mol->nedges].vertice[0]=mol->rsize[_j]+top->res[mol->ress->list[_j]].edges[_i].vertice[1];
                   mol->edges[mol->nedges].vertice[1]=mol->rsize[_j+1]+top->res[mol->ress->list[_j+1]].edges[_l].vertice[1];
                   mol->edges[mol->nedges++].type=-top->res[mol->ress->list[_j]].edges[_i].type;
                   }
                 else goto LABEL_DATA_CONSISTMENT_ERROR; // '+' '-' mismatch
                 }
               }
        else if (((char*)&top->res[mol->ress->list[_j]].edges[_i].vertice[0])[0]=='-')
               {
               if (!_j) goto LABEL_DATA_CONSISTMENT_ERROR;
               if (mol->ress->list[_j-1]==(unsigned int)-1) 
                 {//Add all bonds 
                 for (_l=mol->rsize[_j-1];_l<mol->rsize[_j];_l++, mol->nedges++)
                   {
                   mol->edges[mol->nedges].vertice[0]=_l;
                   mol->edges[mol->nedges].vertice[1]=mol->rsize[_j]+top->res[mol->ress->list[_j]].edges[_i].vertice[1];
                   mol->edges[mol->nedges].type='1';
                   }  
                 }
               }  
        }
      else
        { //Polimerize in residue
        mol->edges[mol->nedges].vertice[0]=top->res[mol->ress->list[_j]].edges[_i].vertice[0]+mol->rsize[_j];
        mol->edges[mol->nedges].vertice[1]=top->res[mol->ress->list[_j]].edges[_i].vertice[1]+mol->rsize[_j];
        mol->edges[mol->nedges++].type=top->res[mol->ress->list[_j]].edges[_i].type;
        }
    }

//Polymerize nonlineary
while (nline_size)
  {//Scan matrix
  _rmin=(double)NAN;
  _j=nline_size;
  while (--_j)
    if ( (nline_edges[_j].type==nline_edges[ 0].type)                                                                 && 
         (top->res[mol->ress->list[nline_res[_j]]].atoms.list[nline_edges[_j].vertice[1]]==nline_edges[ 0].vertice[0])&&
         (top->res[mol->ress->list[nline_res[ 0]]].atoms.list[nline_edges[ 0].vertice[1]]==nline_edges[_j].vertice[0]) )
      { //Topological correspondence atchieved, check geometrical
      if ((_r=calc_distance(&mol->r[mol->rsize[nline_res[ 0]]+nline_edges[ 0].vertice[1]],
                            &mol->r[mol->rsize[nline_res[_j]]+nline_edges[_j].vertice[1]]))<MAX_CHEM_BOND_LENGHT*MAX_CHEM_BOND_LENGHT) //Cutoff for bond formation
        if ( (isnan(_rmin))||(_r<_rmin) ) { _rmin=_r, item=_j; }
      }
  if (!isnan(_rmin)) 
    { 
    //Plomerize nonlinearly
    mol->edges[mol->nedges].vertice[0]=mol->rsize[nline_res[   0]]+nline_edges[   0].vertice[1];
    mol->edges[mol->nedges].vertice[1]=mol->rsize[nline_res[item]]+nline_edges[item].vertice[1];
    mol->edges[mol->nedges].type=nline_edges[   0].type;
    mol->nedges++;
    //Move j-th bond
    if (--nline_size!=item)
      {
      nline_res[item]=nline_res[nline_size];
      nline_edges[item].vertice[0]=nline_edges[nline_size].vertice[0];
      nline_edges[item].vertice[1]=nline_edges[nline_size].vertice[1];
      nline_edges[item].type=nline_edges[nline_size].type;
      }
    }
  else
    {//Lock it with a protons
    if (((mol->atoms->size%0xFF)+nline_edges[   0].type)>=0xFF)
      {
      _i=mol->atoms->size+nline_edges[   0].type;
      if (!(vp=(void*)realloc_list(mol->atoms,_i))) goto LABEL_MEMORY_ERROR;
      else mol->atoms=(t_list*)vp;
      if (!(vp=(void*)realloc(mol->anames,sizeof(int)*_i))) goto LABEL_MEMORY_ERROR;
      else mol->anames=(char (*)[sizeof(int)])vp; 
      if (!(vp=(void*)realloc(mol->r,sizeof(t_vec)*_i))) goto LABEL_MEMORY_ERROR;
      else mol->r=(t_vec*)vp;
      }
    //Move atoms
    _k=nline_edges[   0].type-'0', _l=mol->atoms->size; 
    while (_l!=mol->rsize[nline_res[0]+1]) 
      {
      --_l;
      mol->atoms->list[_l+_k]=mol->atoms->list[_l], mol->r[_l+_k].i=mol->r[_l].i, mol->r[_l+_k].j=mol->r[_l].j, mol->r[_l+_k].k=mol->r[_l].k;
      *((int*)&mol->anames[_l+_k])=*((int*)&mol->anames[_l]); 
      }
    //Move edges
    _l=mol->nedges;
    while (_l--)
      {
      if (mol->edges[_l].vertice[0]>=mol->rsize[nline_res[0]+1]) mol->edges[_l].vertice[0]+=_k;
      if (mol->edges[_l].vertice[1]>=mol->rsize[nline_res[0]+1]) mol->edges[_l].vertice[1]+=_k;
      }
    //Insert atoms and edges
    for (_l=0; _l<_k; _l++, mol->atoms->size++, mol->nedges++)
      {
      mol->edges[mol->nedges].vertice[0]=mol->rsize[nline_res[   0]]+nline_edges[   0].vertice[1];
      mol->edges[mol->nedges].vertice[1]=mol->rsize[nline_res[   0]+1]+_l;
      mol->edges[mol->nedges].type='1';
      mol->atoms->list[mol->rsize[nline_res[   0]+1]+_l]=1;
      mol->anames[mol->rsize[nline_res[   0]+1]+_l][0]=mol->anames[mol->rsize[nline_res[   0]+1]+_l][1]='_', mol->anames[mol->rsize[nline_res[   0]+1]+_l][2]='H', mol->anames[mol->rsize[nline_res[   0]+1]+_l][3]=(char)(_l+'1');
      mol->r[mol->rsize[nline_res[   0]+1]+_l].i=mol->r[mol->rsize[nline_res[   0]+1]+_l].j=mol->r[mol->rsize[nline_res[   0]+1]+_l].k=(double)NAN;
      }
    _l=mol->ress->size; while (_l!=nline_res[0]) { mol->rsize[_l]+=_k, _l--; }
    }
  nline_size--;
  nline_res[0]=nline_res[nline_size];
  nline_edges[0].vertice[0]=nline_edges[nline_size].vertice[0];
  nline_edges[0].vertice[1]=nline_edges[nline_size].vertice[1];
  nline_edges[0].type=nline_edges[nline_size].type;
  }


//Syncronize all
if (!(vp=(void*)realloc_list(mol->atoms,mol->atoms->size))) goto LABEL_MEMORY_ERROR;
else mol->atoms=(t_list*)vp;
if (!(vp=(void*)realloc(mol->r,sizeof(t_vec)*mol->atoms->size))) goto LABEL_MEMORY_ERROR;
else mol->r=(t_vec*)vp;
if (!(vp=(void*)realloc(mol->edges,sizeof(t_edge)*mol->nedges))) goto LABEL_MEMORY_ERROR;
else mol->edges=(t_edge*)vp;

if (nline_edges) free(nline_edges);
if (nline_res)   free(nline_res);

return mol;
}

*/

//This function read structure from pdb/gro file. The str is ready for conversion into a mol.
/*t_str *read_pdb(FILE *in)
{
t_str *str;
extern char buffer[0xFF];
char *line;

//Stage I. Read straingforward str
while (!(line=yfgets(in)))
  if ( ( (line[0]=='A')&&(str[1]=='T')&&(line[2]=='O')&&(line[3]=='M')&&(line[4]==' ')&&(line[5]==' ')&&(line[6]==' ') )||
       ( (line[0]=='H')&&(str[1]=='E')&&(line[2]=='T')&&(line[3]=='A')&&(line[4]=='T')&&(line[5]=='M')&&(line[6]==' ') ) )
    {
    *((int*)str->anames[str->size_a])=*((int*)&line[]);
    }
}
*/
//This function import coordinates from PDBQT file (id-th structure; it return ylib_errno==YERROR_LEGAL if there is no id-th structure in the input file).
//It doesn't produce a bond matrix, just residues, atom names and corresponding coordinates
t_str *read_pdbqt(unsigned int id,FILE *in)
{
register unsigned int _i, _j, len;
t_str *str;
extern char buffer[];
char *string;
unsigned int resid;
void *vp;

//Skip till id-th structure 
while (id)
  { 
  if (!(string=yfgets(in))) { ylib_errno=YERROR_LEGAL; return FALSE; }
  if ( (strlen(string)>=strlen("ENDMDL"))&&(!(strncmp(string,"ENDMDL",strlen("ENDMDL")))) ) id--;
  if (string!=buffer) { free(string), string=0x0; }
  }
//Preallocate memory
if (!(str=(t_str*)malloc(sizeof(t_str)))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return FALSE; } 
else { str->natoms=0, str->name=0x0, str->nedges=0, str->edges=0x0; }
if (!(str->ress=alloc_list(0xFF))) { LABEL_MEMORY_ERROR_1: free(str); goto LABEL_MEMORY_ERROR_0; }
else str->ress->size=0;
if (!(str->rsize=(unsigned int*)malloc(sizeof(unsigned int)*(0xFF+0x1)))) { LABEL_MEMORY_ERROR_2: free(str->ress); goto LABEL_MEMORY_ERROR_1; }
if (!(str->a=(char*)malloc(sizeof(char)*0xFF))) { LABEL_MEMORY_ERROR_3: free(str->rsize); goto LABEL_MEMORY_ERROR_2; }
if (!(str->anames=(char(*)[sizeof(int)])malloc(sizeof(int)*0xFF))) { LABEL_MEMORY_ERROR_4: free(str->a); goto LABEL_MEMORY_ERROR_3; }
if (!(str->r=(t_vec*)malloc(sizeof(t_vec)*0xFF))) { LABEL_MEMORY_ERROR_5: free(str->anames); goto LABEL_MEMORY_ERROR_4; }
//Upload id-th structure
while ( (string=yfgets(in)))
  {
  len=strlen(string);
  if ( (len>=strlen("ENDMDL"))&&(!(strncmp(string,"ENDMDL",strlen("ENDMDL")))) ) 
    { if (string!=buffer) { free(string), string=0x0; } break; }
  if ( (len==80)&&( (!(strncmp(string,"ATOM  ",strlen("ATOM  "))))||(!(strncmp(string,"HETATM",strlen("HETATM")))) ) )
    {
    //Handle residue data
    string[26]='\0', resid=atoi(&string[20]);
    if (!(str->ress->size)) { str->start_rid=resid; str->ress->list[0]=*((unsigned int*)&string[16]), str->ress->size=1; }
    else
      {
      if ( (str->start_rid+str->ress->size!=1+resid)||(str->ress->list[str->ress->size-1]!=*((unsigned int*)&string[16])) )
        {
        if ((str->start_rid+str->ress->size!=2+resid)) { LABEL_DATA_CONSISTMENT_FAILURE: free_str(str); ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; } //No gaps within residues sequence are permitted 
        str->rsize[str->ress->size]=str->natoms, str->ress->list[str->ress->size]=*((unsigned int*)&string[16]);
        if (!(++str->ress->size%0xFF))
          {
          if (!(vp=(t_list*)realloc_list(str->ress,str->ress->size+0xFF))) { LABEL_MEMORY_ERROR_6: free(str->r); goto LABEL_MEMORY_ERROR_5; }
          else str->ress=(t_list*)vp;
          if (!(vp=realloc(str->rsize,sizeof(unsigned int)*(str->ress->size+0xFF+0x1)))) goto LABEL_MEMORY_ERROR_6;
          else str->rsize=(unsigned int*)vp;
          }
        }
      }
    //Handle atom data 
    switch(string[77])
      {//Atom type
      case 'A' : { if (string[78]==' ') str->a[str->natoms]=CHEM_ATOM_TYPE_CARBON;   else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'B' : { if (string[78]=='r') str->a[str->natoms]=CHEM_ATOM_TYPE_BROMINE;  else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'C' : { if (string[78]==' ') str->a[str->natoms]=CHEM_ATOM_TYPE_CARBON;   else if (string[78]=='l') str->a[str->natoms]=CHEM_ATOM_TYPE_CHLORINE; else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'F' : { if (string[78]==' ') str->a[str->natoms]=CHEM_ATOM_TYPE_FLUORINE; else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'G' : { if ( (string[78]==' ')||(string[78]=='A') ) str->a[str->natoms]=CHEM_ATOM_TYPE_CARBON; else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'H' : { if ( (string[78]==' ')||(string[78]=='S')||(string[78]=='D') ) str->a[str->natoms]=CHEM_ATOM_TYPE_HYDROGEN; else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'I' : { if (string[78]==' ') str->a[str->natoms]=CHEM_ATOM_TYPE_IODINE;   else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'J' : { if (string[78]==' ') str->a[str->natoms]=CHEM_ATOM_TYPE_CARBON;   else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'N' : { if ( (string[78]==' ')||(string[78]=='A')||(string[78]=='S') ) str->a[str->natoms]=CHEM_ATOM_TYPE_NITROGEN; else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'O' : { if ( (string[78]==' ')||(string[78]=='S')||(string[78]=='A') ) str->a[str->natoms]=CHEM_ATOM_TYPE_OXYGEN;   else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'P' : { if (string[78]==' ') str->a[str->natoms]=CHEM_ATOM_TYPE_PHOSPHOR; else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'Q' : { if (string[78]==' ') str->a[str->natoms]=CHEM_ATOM_TYPE_CARBON;   else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      case 'S' : { if ( (string[78]==' ')||(string[78]=='A') ) str->a[str->natoms]=CHEM_ATOM_TYPE_SULFUR; else if (string[78]=='i') str->a[str->natoms]=CHEM_ATOM_TYPE_SILICON; else goto LABEL_DATA_CONSISTMENT_FAILURE; break; }
      default  : goto LABEL_DATA_CONSISTMENT_FAILURE;
      }
    //Atom name
    MAKE_INT_NAME(_i,_j,str->anames[str->natoms],(&string[12]));
    //Atom coords
    string[8]='\0';
    memcpy(string,&string[30],sizeof(char)*8); if ( (check_double_type(string))) str->r[str->natoms].i=atof(string); else goto LABEL_DATA_CONSISTMENT_FAILURE; 
    memcpy(string,&string[38],sizeof(char)*8); if ( (check_double_type(string))) str->r[str->natoms].j=atof(string); else goto LABEL_DATA_CONSISTMENT_FAILURE; 
    memcpy(string,&string[46],sizeof(char)*8); if ( (check_double_type(string))) str->r[str->natoms].k=atof(string); else goto LABEL_DATA_CONSISTMENT_FAILURE; 
    if (!(++str->natoms%0xFF))
      {//Memory management
      if (!(vp=realloc(str->a,sizeof(char)*(str->natoms+0xFF)))) goto LABEL_MEMORY_ERROR_6;
      else str->a=(char*)vp;
      if (!(vp=realloc(str->anames,sizeof(int)*(str->natoms+0xFF)))) goto LABEL_MEMORY_ERROR_6;
      else str->anames=(char(*)[sizeof(int)])vp;
      if (!(vp=realloc(str->r,sizeof(t_vec)*(str->natoms+0xFF)))) goto LABEL_MEMORY_ERROR_6;
      else str->r=(t_vec*)vp;
      }
    }
  if (string!=buffer) { free(string), string=0x0; } //Free memory 
  }
if (!(str->natoms)) { ylib_errno=YERROR_LEGAL; free_str(str); return FALSE; } else str->rsize[str->ress->size]=str->natoms;
//Sync memory & Exit
if (str->ress->size==1)
  {
  free(str->rsize);
  *((unsigned int*)&str->rsize)=*str->ress->list;
  free(str->ress); str->ress=0x0;
  }
else
  {
  if (!(vp=realloc(str->rsize,sizeof(unsigned int)*(str->ress->size+1)))) goto LABEL_MEMORY_ERROR_6;
  else str->rsize=(unsigned int*)vp;
  if (!(vp=(void*)realloc_list(str->ress,str->ress->size))) goto LABEL_MEMORY_ERROR_6;
  else str->ress=(t_list*)vp;
  }
if (!(vp=realloc(str->a,sizeof(char)*str->natoms))) goto LABEL_MEMORY_ERROR_6;
else str->a=(char*)vp;
if (!(vp=realloc(str->anames,sizeof(int)*str->natoms))) goto LABEL_MEMORY_ERROR_6;
else str->anames=(char(*)[sizeof(int)])vp;
if (!(vp=realloc(str->r,sizeof(t_vec)*str->natoms))) goto LABEL_MEMORY_ERROR_6;
else str->r=(t_vec*)vp;
return str;
}


//This function write pdbqt file
//Note. This function requires superanchor and superrtree to save correctly 
inline char _write_pdbqt_atom_record(FILE *out,unsigned int _id,unsigned int *atom_id,char (*order)[4],t_clist *neighbors,t_mol *mol)
{
register unsigned int _i, _j, _l;
char atd[2];

//Stage I. Define the autodock atom type
switch (mol->a[_id])
  {
  case CHEM_ATOM_TYPE_HYDROGEN : { //Hydrgen types "H ", "HD"
    atd[0]='H';
    if (neighbors->list[_id].size!=1) goto LABEL_DATA_CONSISTMENT_FAILURE;
    atd[1]=( (mol->a[_i=*neighbors->list[_id].list]==CHEM_ATOM_TYPE_NITROGEN)||(mol->a[_i]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[_i]==CHEM_ATOM_TYPE_SULFUR) ) ? 'D' : ' '; 
    break; }
  case CHEM_ATOM_TYPE_CARBON   : { //Carbon types "C " or "A "
    atd[0]= (order[_id][0]+order[_id][1]==neighbors->list[_id].size) ? 'C' : 'A', atd[1]=' ';
    break; }
  case CHEM_ATOM_TYPE_NITROGEN : { //Nitrogen types "N ", "NA" or "NS"
    atd[0]='N', atd[1]=' ';
    if (order[_id][0]==neighbors->list[_id].size)
//      { if (neighbors->list[_id].size==3) { _j=0, _l=neighbors->list[_id].size; while (_l--) if ( ( (order[neighbors->list[_id].list[_l]][2]))||( (order[neighbors->list[_id].list[_l]][3])) ) _j++; if (!_j) atd[1]='S'; } } //NS is not valid for vina
      { if (neighbors->list[_id].size==3) { _j=0, _l=neighbors->list[_id].size; while (_l--) if ( ( (order[neighbors->list[_id].list[_l]][2]))||( (order[neighbors->list[_id].list[_l]][3])) ) _j++; if (!_j) atd[1]=' '; } }
    else { if ( ( (order[_id][2]==1)&&(neighbors->list[_id].size==2) )||( (order[_id][3]==1)&&(neighbors->list[_id].size==1) ) ) atd[1]='A'; }
    break;
    }
  case CHEM_ATOM_TYPE_OXYGEN   : { //Oxygen types "OA" or "O "
    atd[0]='O', atd[1]=(neighbors->list[_id].size==1) ? 'A' : ' ';
    break; }
  case CHEM_ATOM_TYPE_FLUORINE : { atd[0]='F', atd[1]=' '; break; } //Fluorine, only "F "
  case CHEM_ATOM_TYPE_SILICON  : { atd[0]='P', atd[1]=' '; break; } //Silicone, only "P "
  case CHEM_ATOM_TYPE_PHOSPHOR : { atd[0]='P', atd[1]=' '; break; } //Phosphor, only "P "
  case CHEM_ATOM_TYPE_SULFUR   : { //Sulfur "SA" or "S "
    atd[0]='S', atd[1]=(neighbors->list[_id].size==1) ? 'A' : ' ';
    break; }
  case CHEM_ATOM_TYPE_CHLORINE : { atd[0]='C', atd[1]='l'; break; } //Chlorine, only "Cl"
  case CHEM_ATOM_TYPE_BROMINE  : { atd[0]='B', atd[1]='r'; break; } //Bromine,  only "Br"
  case CHEM_ATOM_TYPE_IODINE   : { atd[0]='I', atd[1]=' '; break; } //Iodane,   only "I "
  default                      : { LABEL_DATA_CONSISTMENT_FAILURE: ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
  }
//Stage II. Get resid and resname
_i=1; while (mol->rsize[_i]<=_id) _i++; //get _i+1==resid
//Stage III. Export the record and exit
return (fprintf(out,"ATOM  %5.1d %4.4s %4.4s%5d    %8.3f%8.3f%8.3f  1.00  0.00    %6.3f %2.2s\n",atom_id[_id],mol->anames[_id],(char*)&mol->ress->list[_i-1],_i,mol->r[_id].i,mol->r[_id].j,mol->r[_id].k,mol->charges[_id],atd)>0);
}
char write_pdbqt(FILE *out,char (*order)[4],t_clist *neighbors,t_clist *sanchors,t_rtree *srtree,t_mol *mol)
{
register unsigned int _i, _l, atom_count=1;
unsigned int *atom_id; 

//Stage I. Prepare memory
if (!(atom_id=(unsigned int*)calloc(mol->natoms,sizeof(unsigned int)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
//Stage II. Export atoms
//Stage II.1. Export root atoms
_l=sanchors->list[_i=srtree->root].size; while (_l--) atom_id[sanchors->list[_i].list[_l]]=atom_count++;
if (fwrite((const char*)"ROOT \n",sizeof(char)*strlen("ROOT \n"),1,out)!=1) { LABEL_IO_ERROR: ylib_errno=YERROR_IO; free(atom_id); return FALSE; }
_l=sanchors->list[_i].size; while (_l--) if (!(_write_pdbqt_atom_record(out,sanchors->list[_i].list[_l],atom_id,order,neighbors,mol))) { LABEL_ERROR: free(atom_id); return FALSE; }
if (fwrite((const char*)"ENDROOT \n",sizeof(char)*strlen("ENDROOT \n"),1,out)!=1) goto LABEL_IO_ERROR;
srtree->rbranch[_i].edge.type=srtree->rbranch[_i].nrbranch;
while (srtree->rbranch[_i].edge.type--)
  {
  //Stage II.2. Export root BRANCHES
  _i=srtree->rbranch[_i].rbranch[srtree->rbranch[_i].edge.type];
  srtree->rbranch[_i].edge.type=srtree->rbranch[_i].nrbranch;
  _l=sanchors->list[_i].size; while (_l--) atom_id[sanchors->list[_i].list[_l]]=atom_count++;
  if (fprintf(out,"BRANCH %5.1d %5.1d \n",atom_id[srtree->rbranch[_i].edge.vertice[1]],atom_id[srtree->rbranch[_i].edge.vertice[0]])<0) goto LABEL_IO_ERROR;
  _l=sanchors->list[_i].size; while (_l--) if (!(_write_pdbqt_atom_record(out,sanchors->list[_i].list[_l],atom_id,order,neighbors,mol))) goto LABEL_ERROR;
  while (_i!=srtree->root)
    {//Stage II.3. Export non-root BRANCHES
    while(--srtree->rbranch[_i].edge.type)
      {
      _i=srtree->rbranch[_i].rbranch[srtree->rbranch[_i].edge.type];
      srtree->rbranch[_i].edge.type=srtree->rbranch[_i].nrbranch;
      _l=sanchors->list[_i].size; while (_l--) atom_id[sanchors->list[_i].list[_l]]=atom_count++;
      if (fprintf(out,"BRANCH %5.1d %5.1d \n",atom_id[srtree->rbranch[_i].edge.vertice[1]],atom_id[srtree->rbranch[_i].edge.vertice[0]])<0) goto LABEL_IO_ERROR;
      _l=sanchors->list[_i].size; while (_l--) if (!(_write_pdbqt_atom_record(out,sanchors->list[_i].list[_l],atom_id,order,neighbors,mol))) goto LABEL_ERROR;
      }
    //Stage II.4. Export ENDBRANCH record
    if (fprintf(out,"ENDBRANCH %5.1d %5.1d \n",atom_id[srtree->rbranch[_i].edge.vertice[1]],atom_id[srtree->rbranch[_i].edge.vertice[0]])<0) goto LABEL_IO_ERROR;
    _i=*srtree->rbranch[_i].rbranch;
    }
  }
//Stage II.5 Export TORSDOF
if (fprintf(out,"TORSDOF %5.1d \n\n",srtree->nidofs)<0) goto LABEL_IO_ERROR;
//Free some memory and exit
free(atom_id);
return TRUE;
}

//This function define amber99sb-ildn - compatible atom type
int _define_amber99sb_ildn_atom_type(unsigned int _id,char (*order)[4],t_clist *neighbors,t_mol *mol)
{
register unsigned int _i, _j;
switch (mol->a[_id])
  {
  case CHEM_ATOM_TYPE_HYDROGEN : { // H, HC, HA, HO, HS
    if (neighbors->list[_id].size==1)
      switch (mol->a[*neighbors->list[_id].list])
        {
        case CHEM_ATOM_TYPE_CARBON   : { if ( ( (order[*neighbors->list[_id].list][2]))||( (order[*neighbors->list[_id].list][3])) ) return *((int*)"HA "); else return *((int*)"HC "); }
        case CHEM_ATOM_TYPE_NITROGEN : { return *((int*)"H  "); } 
        case CHEM_ATOM_TYPE_OXYGEN   : { return *((int*)"HO "); }
        case CHEM_ATOM_TYPE_SULFUR   : { return *((int*)"HS "); }
        default                      : { return *((int*)"HC "); }
        }
    else return *((int*)"HC ");
    }
  case CHEM_ATOM_TYPE_CARBON   : { // C, CA, CT
    if ( (order[_id][2])) { _i=mol->size_ar; while (_i--) { _j=mol->cycles->list[_i].size; while (_j--) if (mol->cycles->list[_i].list[_j]==_id) return *((int*)"CA "); } } 
    if ( ( (order[_id][2]))||( (order[_id][3])) ) return *((int*)"C  ");
    else                                          return *((int*)"CT ");
    }     
  case CHEM_ATOM_TYPE_NITROGEN : { // N, NA, N3, NC, N2
    _i=mol->size_ar; while (_i--) { _j=mol->cycles->list[_i].size; while (_j--) if (mol->cycles->list[_i].list[_j]==_id) { if (neighbors->list[_id].size==2) return *((int*)"NC "); else return *((int*)"NA "); } } 
    _i=mol->nvedges; while (_i--) if ( (mol->vedges[_i][0]==_id)&&(mol->vatoms[mol->vedges[_i][1]]>0) ) return *((int*)"N3 ");
    if (neighbors->list[_id].size==2) return *((int*)"N2 "); else return *((int*)"N  ");
    } 
  case CHEM_ATOM_TYPE_OXYGEN   : { // O, OS, OH
         if (neighbors->list[_id].size==1) return *((int*)"O  "); 
    else if ( (neighbors->list[_id].size==2)&&
            ( (mol->a[neighbors->list[_id].list[0]]==CHEM_ATOM_TYPE_HYDROGEN)||(mol->a[neighbors->list[_id].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) ) ) return *((int*)"OH ");
    else return *((int*)"OS ");
    }
  case CHEM_ATOM_TYPE_FLUORINE : { return *((int*)"F  "); } 
  case CHEM_ATOM_TYPE_SILICON  : { return *((int*)"P  "); } //There is no Si atom type, using P
  case CHEM_ATOM_TYPE_PHOSPHOR : { return *((int*)"P  "); }
  case CHEM_ATOM_TYPE_SULFUR   : { // S, SH
    if ( (neighbors->list[_id].size==2)&&
       ( (mol->a[neighbors->list[_id].list[0]]==CHEM_ATOM_TYPE_HYDROGEN)||(mol->a[neighbors->list[_id].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) ) ) return *((int*)"SH ");
    else return *((int*)"S  ");
    } 
  case CHEM_ATOM_TYPE_CHLORINE : { return *((int*)"Cl "); }
  case CHEM_ATOM_TYPE_BROMINE  : { return *((int*)"Br "); }
  case CHEM_ATOM_TYPE_IODINE   : { return *((int*)"I  "); }
  default                      : {
    switch (neighbors->list[_id].size) 
      {
      case 4  : return *((int*)"CT ");
      case 3  : return *((int*)"C  ");
      case 2  : return *((int*)"C  ");
      case 1  : return *((int*)"Cl ");
      default : return *((int*)"Cl "); 
      }
    }
  }
}
//This function writes gromacs itp from parameterized mol
char export_itp(FILE *out,char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_mol *mol,t_top *top)
{
register unsigned int _i, _j;
unsigned int amber_atom_type, *chrg_grps;
double a, b, c;
char flag;
//Define (simple) charged groups
if (!(chrg_grps=split_charges(neighbors,anchor_id,mol,top))) return FALSE;

//Write header 
if (fprintf(out,";\n;The topology is generated by %s for amber99sb_ildn all-atoms FF\n;\n",Y_APPLICATION)<0) { LABEL_IO_ERROR_0: free(chrg_grps); goto LABEL_IO_ERROR; }
//Write moleculetype 
if (fprintf(out,"\n[ moleculetype ]\n")<0) goto LABEL_IO_ERROR_0;
if (fprintf(out,"; Name      nrexcl\n")<0) goto LABEL_IO_ERROR_0;
if (fprintf(out,"DRG           3 \n")<0)   goto LABEL_IO_ERROR_0;

//Write atoms
if (fprintf(out,"\n[ atoms ]\n")<0) goto LABEL_IO_ERROR_0;
if (fprintf(out,";   nr   type  resnr  resid  atom    cgnr   charge      mass\n")<0) goto LABEL_IO_ERROR_0;
for (_i=_j=0;_i<mol->natoms;_i++)
  {
  amber_atom_type=_define_amber99sb_ildn_atom_type(_i,order,neighbors,mol);
  if (fprintf(out," %5d  %4.4s  %5d  %4.4s  %4.4s  %5d  %8.6f  %8.4f\n",_i+1,(char*)&amber_atom_type,mol->start_rid+_j,( (mol->ress->list[_j])) ? (char*)&mol->ress->list[_j] : "DRG",mol->anames[_i],chrg_grps[_i]+1,mol->charges[_i],top->ff_a[mol->ytypes[_i]].mass)<0) goto LABEL_IO_ERROR_0;
  if (_i==mol->rsize[_j+1]) _j++; 
  }
free(chrg_grps), chrg_grps=0x0; 
//Write bonds
if (fprintf(out,"\n[ bonds ]\n")<0) { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if (fprintf(out,";  ai      aj    fu    c0,     c1, ...\n")<0) goto LABEL_IO_ERROR;
for (_i=0; _i<mol->size_b; _i++)
  if (fprintf(out,"%5d  %5d  2  %12.8e  %12.8e  %12.8e  %12.8e  ;  %4.4s - %4.4s\n",mol->ff_b[_i].atom[0]+1,mol->ff_b[_i].atom[1]+1,0.1*mol->ff_b[_i].v,10000.0*mol->ff_b[_i].k,0.1*mol->ff_b[_i].v,10000.*mol->ff_b[_i].k,mol->anames[mol->ff_b[_i].atom[0]],mol->anames[mol->ff_b[_i].atom[1]])<0) goto LABEL_IO_ERROR; //Convert it into kJ/nm and nm
//Write X-[O,S]-H constraints if necessary
flag=FALSE; _i=mol->natoms;
while (_i--)
  if ( (neighbors->list[_i].size==2)&&( (mol->a[_i]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[_i]==CHEM_ATOM_TYPE_SULFUR) )&&
     ( (mol->a[neighbors->list[_i].list[0]]==CHEM_ATOM_TYPE_HYDROGEN)||(mol->a[neighbors->list[_i].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
    {
    if (!(flag)) { flag=TRUE; if (fprintf(out,"\n[ constraints ]\n")<0) goto LABEL_IO_ERROR; }
    _j=mol->size_b;
    while (--_j)
      if ( ( (mol->ff_b[_j].atom[0]==_i)&&(mol->ff_b[_j].atom[1]==neighbors->list[_i].list[0]) )||
           ( (mol->ff_b[_j].atom[0]==neighbors->list[_i].list[0])&&(mol->ff_b[_j].atom[1]==_i) ) ) { a=mol->ff_b[_j].v; break; }
    _j=mol->size_b;
    while (--_j)
      if ( ( (mol->ff_b[_j].atom[0]==_i)&&(mol->ff_b[_j].atom[1]==neighbors->list[_i].list[1]) )||
           ( (mol->ff_b[_j].atom[0]==neighbors->list[_i].list[1])&&(mol->ff_b[_j].atom[1]==_i) ) ) { b=mol->ff_b[_j].v; break; }
    _j=mol->size_g;
    while (--_j)
      if ( (mol->ff_g[_j].atom[1]==_i)&&( ( (mol->ff_g[_j].atom[0]==neighbors->list[_i].list[0])&&(mol->ff_g[_j].atom[2]==neighbors->list[_i].list[1]) )||
                                          ( (mol->ff_g[_j].atom[2]==neighbors->list[_i].list[0])&&(mol->ff_g[_j].atom[0]==neighbors->list[_i].list[1]) ) ) ) { c=cos(mol->ff_g[_j].v); break; }
    c=0.1*sqrt(a*a+b*b-2.*a*b*cos(mol->ff_g[_j].v));
    if (fprintf(out,"%5d  %5d  1  %12.8e  %12.8e\n",neighbors->list[_i].list[0]+1,neighbors->list[_i].list[1]+1,c,c)<0) goto LABEL_IO_ERROR;
    }

//Write angles
if (fprintf(out,"\n[ angles ]\n")<0) goto LABEL_IO_ERROR;
if (fprintf(out,";  ai      aj      ak    fu     c0,         c1, ...\n")<0) goto LABEL_IO_ERROR;
for (_i=0; _i<mol->size_g; _i++)
  if (fprintf(out,"%5d  %5d  %5d  2  %5.1f  %12.8e  %5.1f  %12.8e ;  %4.4s - %4.4s - %4.4s\n",mol->ff_g[_i].atom[0]+1,mol->ff_g[_i].atom[1]+1,mol->ff_g[_i].atom[2]+1,mol->ff_g[_i].v*180.0/PI,2.0*mol->ff_g[_i].k,mol->ff_g[_i].v*180.0/PI,2.0*mol->ff_g[_i].k,mol->anames[mol->ff_g[_i].atom[0]],mol->anames[mol->ff_g[_i].atom[1]],mol->anames[mol->ff_g[_i].atom[2]])<0) goto LABEL_IO_ERROR;

//Write dihedrals
if (fprintf(out,"\n[ dihedrals ]\n")<0) goto LABEL_IO_ERROR;
//Write impropers
if (fprintf(out,"; impropers\n")<0) goto LABEL_IO_ERROR;
if (fprintf(out,";  ai      aj      ak      al    fu     c0,    c1, ...\n")<0) goto LABEL_IO_ERROR;
//Save type0 imprs
for (_i=0; _i<mol->size_i; _i++)
  if (!(mol->ff_i[_i].type))
    if (fprintf(out,"%5d  %5d  %5d  %5d  2  %12.8e %12.8e  %12.8e  %12.8e  ;  optic imp0  %4.4s -> ( %4.4s %4.4s %4.4s )\n",mol->ff_i[_i].atom[0]+1,mol->ff_i[_i].atom[1]+1,mol->ff_i[_i].atom[2]+1,mol->ff_i[_i].atom[3]+1,mol->ff_i[_i].v*180./PI,mol->ff_i[_i].k,mol->ff_i[_i].v*180./PI,mol->ff_i[_i].k,mol->anames[mol->ff_i[_i].atom[0]],mol->anames[mol->ff_i[_i].atom[1]],mol->anames[mol->ff_i[_i].atom[2]],mol->anames[mol->ff_i[_i].atom[3]])<0) goto LABEL_IO_ERROR; 
//Save type1 imprs
for (_i=0; _i<mol->size_i; _i++)
  if ( (mol->ff_i[_i].type))
    if (fprintf(out,"%5d  %5d  %5d  %5d  2  %12.8e %12.8e  %12.8e  %12.8e  ;  imp1  %4.4s - %4.4s - %4.4s - %4.4s \n",mol->ff_i[_i].atom[0]+1,mol->ff_i[_i].atom[1]+1,mol->ff_i[_i].atom[2]+1,mol->ff_i[_i].atom[3]+1,mol->ff_i[_i].v*180./PI,mol->ff_i[_i].k,mol->ff_i[_i].v*180./PI,mol->ff_i[_i].k,mol->anames[mol->ff_i[_i].atom[0]],mol->anames[mol->ff_i[_i].atom[1]],mol->anames[mol->ff_i[_i].atom[2]],mol->anames[mol->ff_i[_i].atom[3]])<0) goto LABEL_IO_ERROR; 
//Write torsions
if (fprintf(out,"; torsions\n")<0) goto LABEL_IO_ERROR;
if (fprintf(out,";  ai      aj      ak      al    fu     c0,    c1,     m ...\n")<0) goto LABEL_IO_ERROR;
for (_i=0; _i<mol->size_d; _i++)
  if (!(mol->ff_d[_i].d))
    {
    for (_j=0; _j<SIZE_DIH; _j++)
      if (mol->ff_d[_i].k[_j]!=0.)
        {
        if (fprintf(out,"%5d  %5d  %5d  %5d  4  %12.8e %12.8e  %5.2f  %12.8e  %12.8e  %5.2f  ;  quasy-imp0  %4.4s -> ( %4.4s %4.4s %4.4s )\n",mol->ff_d[_i].atom[0]+1,mol->ff_d[_i].atom[1]+1,mol->ff_d[_i].atom[2]+1,mol->ff_d[_i].atom[3]+1,mol->ff_d[_i].v[_j]*180./PI,mol->ff_d[_i].k[_j],mol->ff_d[_i].n[_j],mol->ff_d[_i].v[_j]*180./PI,mol->ff_d[_i].k[_j],mol->ff_d[_i].n[_j],mol->anames[mol->ff_d[_i].atom[0]],mol->anames[mol->ff_d[_i].atom[1]],mol->anames[mol->ff_d[_i].atom[2]],mol->anames[mol->ff_d[_i].atom[3]])<0) goto LABEL_IO_ERROR; 
        }
    }
  else
    {
    for (_j=0; _j<SIZE_DIH; _j++)
      if (mol->ff_d[_i].k[_j]!=0.)
        {
        if (fprintf(out,"%5d  %5d  %5d  %5d  9  %12.8e %12.8e  %5.2f  %12.8e  %12.8e  %5.2f;  dih  %4.4s - %4.4s - %4.4s - %4.4s \n",mol->ff_d[_i].atom[0]+1,mol->ff_d[_i].atom[1]+1,mol->ff_d[_i].atom[2]+1,mol->ff_d[_i].atom[3]+1,mol->ff_d[_i].v[_j]*180./PI,mol->ff_d[_i].k[_j],mol->ff_d[_i].n[_j],mol->ff_d[_i].v[_j]*180./PI,mol->ff_d[_i].k[_j],mol->ff_d[_i].n[_j],mol->anames[mol->ff_d[_i].atom[0]],mol->anames[mol->ff_d[_i].atom[1]],mol->anames[mol->ff_d[_i].atom[2]],mol->anames[mol->ff_d[_i].atom[3]])<0) goto LABEL_IO_ERROR; 
        }
    }
//Save posres. Cant be done with file - "include" mechanism in gromacs
//if (fprintf(out,"\n\n; Include position restraints\n")<0) goto LABEL_IO_ERROR; 
//if (fprintf(out,"#ifdef POSRES\n")<0) goto LABEL_IO_ERROR; 
//if (fprintf(out,"\n[ position_restraints ]\n")<0) goto LABEL_IO_ERROR;
//if (fprintf(out,";  ai    type    fx      fy      fz\n")<0) goto LABEL_IO_ERROR; 
//for (_i=0; _i<mol->natoms; _i++)
//  if (mol->a[_i]!=CHEM_ATOM_TYPE_HYDROGEN)
//    if (fprintf(out,"%5d   1   1000.   1000.   1000.\n",_i+1)<0) goto LABEL_IO_ERROR; 
//if (fprintf(out,"#endif\n")<0) goto LABEL_IO_ERROR; 
if (fprintf(out,"\n")<0) goto LABEL_IO_ERROR; 
//Exit
return TRUE;
}

//This function writes gro file for topology generated with a function above
char write_gro(FILE *out,t_mol *mol)
{
register unsigned int _i, _j;

if (fprintf(out,"GROMACS COORDINATES FROM YMOL\n")<0) { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if (fprintf(out,"%6.1d\n",mol->natoms)<0) goto LABEL_IO_ERROR;
if ( (mol->ress->size==1)&&(!(*mol->ress->list)) )
  {
  for (_i=0; _i<mol->natoms; _i++)
    if (fprintf(out,"%5dDRG    %4.4s%4.1d %7.3f %7.3f %7.3f\n",mol->start_rid,mol->anames[_i],mol->start_rid,0.1*mol->r[_i].i,0.1*mol->r[_i].j,0.1*mol->r[_i].k)<0) goto LABEL_IO_ERROR;
  }
else
  {
  for (_i=0, _j=0; _i<mol->natoms; _i++)
    {
    if (fprintf(out,"%5d%3.3s    %4.4s%4.1d %7.3f %7.3f %7.3f\n",mol->start_rid+_j,(char*)&mol->ress->list[_j],mol->anames[_i],_i+1,0.1*mol->r[_i].i,0.1*mol->r[_i].j,0.1*mol->r[_i].k)<0) goto LABEL_IO_ERROR;
    if (_i==mol->rsize[_j+1]) _j++;
    } 
  }
if (fprintf(out,"   0.00000   0.00000   0.00000\n")<0) goto LABEL_IO_ERROR;
return TRUE;
}


//This function reads str
//NOTE. input str==0 means create new str
t_str *read_str(t_str *str,FILE *in)
{
register char flag;
unsigned int i;

if (fread(&i,sizeof(unsigned int),0x1,in)!=0x1) { IO_ERROR_0: ylib_errno=YERROR_IO; return FALSE; }
if (i!=Y_MAGIC) { ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }  //Check YMAGIC
if ( (str)) 
  { //Init str
  flag=FALSE;
  memset(str,0x0,sizeof(t_str));
  }
else
  { //Create str 
  flag=TRUE;
  if (!(str=(t_str*)calloc(0x1,sizeof(t_str)))) { MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return FALSE; }
  }
if ( (fread(&i,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&str->start_rid,sizeof(unsigned int),0x1,in)!=0x1)||(!(str->start_rid))||
     (fread(&str->natoms,sizeof(unsigned int),0x1,in)!=0x1)||(!(str->natoms))|(fread(&str->nedges,sizeof(unsigned int),0x1,in)!=0x1) )
  { IO_ERROR: free_str(str); if (flag) free(str); goto IO_ERROR_0; }
//Upload str name
if (i)
  {
  if (!(str->name=(char*)malloc(sizeof(char)*(i+1)))) { MEMORY_ERROR: free_str(str); if (flag) free(str); goto MEMORY_ERROR_0; }
  else str->name[i]='\0';
  if (fread(str->name,sizeof(char),i,in)!=i) goto IO_ERROR;
  }
//Upload residues
if (str->start_rid==1)
  {//Upload compressed residues
  if (fread(&str->rsize,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR; else str->ress=0x0;
  }
else
  {//Upload rest of the list
  if ( (!(str->ress=alloc_list(str->start_rid)))||(!(str->rsize=(unsigned int*)malloc(sizeof(unsigned int)*(str->ress->size+1)))) ) goto MEMORY_ERROR; 
  else
    {
    if ( (fread(str->ress->list,sizeof(unsigned int),str->ress->size,  in)!=str->ress->size)||
         (fread(str->rsize,     sizeof(unsigned int),str->ress->size+1,in)!=str->ress->size+1) ) goto IO_ERROR;
    }
  }
if (fread(&str->start_rid,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
//Upload atoms
if (!(str->anames=(char (*)[sizeof(int)])malloc(sizeof(int)*str->natoms))) goto MEMORY_ERROR;
if (!(str->a=(char*)malloc(sizeof(char)*str->natoms))) goto MEMORY_ERROR;
if (!(str->r=(t_vec*)malloc(sizeof(t_vec)*str->natoms))) goto MEMORY_ERROR;
if (fread(str->anames,sizeof(int),str->natoms,in)!=str->natoms) goto IO_ERROR;
if (fread(str->a,sizeof(char),str->natoms,in)!=str->natoms) goto IO_ERROR;
if (fread(str->r,sizeof(t_vec),str->natoms,in)!=str->natoms) goto IO_ERROR;
//Upload edges
if (str->nedges)
  {
  if (!(str->edges=(t_edge*)malloc(sizeof(t_edge)*str->nedges))) goto MEMORY_ERROR;
  if (fread(str->edges,sizeof(t_edge),str->nedges,in)!=str->nedges) goto IO_ERROR;
  }
return str;
}
//This function reads str in compact solid form
t_str *read_solid_str(FILE *in)
{
t_str *str;
unsigned int i, name_len, size_r, natoms, nedges;

//Read constants
if ( (fread(&i,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&name_len,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&size_r,sizeof(unsigned int),0x1,in)!=0x1)||
     (fread(&natoms,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&nedges,sizeof(unsigned int),0x1,in)!=0x1) ) { LABEL_IO_ERROR_0: ylib_errno=YERROR_IO; return FALSE; }
if ( (i!=Y_MAGIC)||(!size_r)||(!natoms) ) { ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
//Allocate solid structure in memory
if (!(str=(t_str*)alloc_solid_str(name_len,size_r,natoms,nedges))) return FALSE; 
//Upload name
if (fread(str->name,sizeof(char),name_len,in)!=name_len) { LABEL_IO_ERROR_1: free(str); str=0x0; goto LABEL_IO_ERROR_0; }
else str->name[name_len]='\0';
//Upload residues
if (size_r==1) { if (fread(&str->rsize,sizeof(unsigned int),0x1,in)!=0x1) goto LABEL_IO_ERROR_1; else str->ress=0x0; }
else           { if ( (fread(str->ress->list,sizeof(unsigned int),size_r,in)!=size_r)||(fread(str->rsize,sizeof(unsigned int),size_r,in)!=size_r) ) goto LABEL_IO_ERROR_1; }
if (fread(&str->start_rid,sizeof(unsigned int),0x1,in)!=0x1) goto LABEL_IO_ERROR_1;
//Upload atoms
str->natoms=natoms;
if (fread(str->anames,sizeof(int),str->natoms,in)!=str->natoms) goto LABEL_IO_ERROR_1;
if (fread(str->a,sizeof(char),str->natoms,in)!=str->natoms) goto LABEL_IO_ERROR_1;
if (fread(str->r,sizeof(t_vec),str->natoms,in)!=str->natoms) goto LABEL_IO_ERROR_1;
//Upload edges
str->nedges=nedges;
if (str->nedges) { if (fread(str->edges,sizeof(t_edge),str->nedges,in)!=str->nedges) goto LABEL_IO_ERROR_1; }
else str->edges=0x0;
return str;
}

//This function writes str
char write_str(FILE *out,t_str *str)
{
unsigned int i, j;

i=Y_MAGIC; if (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1) { IO_ERROR: ylib_errno=YERROR_IO; return FALSE; } //Write YMAGIC
i=0; if ( (str->name)&&(*str->name) ) { while ( (str->name[i])) i++; }
if (!(str->ress))
  {//Compressed residues
  j=1;
  if ( (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1)||
       (fwrite(&str->natoms,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&str->nedges,sizeof(unsigned int),0x1,out)!=0x1) ) goto IO_ERROR;
  if ( (fwrite(str->name,sizeof(char),i,out)!=i)||(fwrite(&str->rsize,sizeof(unsigned int),0x1,out)!=0x1) ) goto IO_ERROR;
  }
else
  {//Full residues 
  if ( (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&str->ress->size,sizeof(unsigned int),0x1,out)!=0x1)||
       (fwrite(&str->natoms,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&str->nedges,sizeof(unsigned int),0x1,out)!=0x1) ) goto IO_ERROR;
  if ( (fwrite(str->name,sizeof(char),i,out)!=i)||
       (fwrite(str->ress->list,sizeof(unsigned int),str->ress->size,out)!=str->ress->size)||
       (fwrite(str->rsize,sizeof(unsigned int),str->ress->size+1,out)!=str->ress->size+1)   ) goto IO_ERROR; //Write residues
  }
if (fwrite(&str->start_rid,sizeof(unsigned int),0x1,out)!=0x1)      goto IO_ERROR; //Save starting residue id of the sequence
if (fwrite(str->anames,sizeof(int),str->natoms,out)!=str->natoms)   goto IO_ERROR; //Write atoms names
if (fwrite(str->a,sizeof(char),str->natoms,out)!=str->natoms)       goto IO_ERROR; //Write atoms names
if (fwrite(str->r,sizeof(t_vec),str->natoms,out)!=str->natoms)      goto IO_ERROR; //Write coords
if (fwrite(str->edges,sizeof(t_edge),str->nedges,out)!=str->nedges) goto IO_ERROR; //Write edges
return TRUE;
}

//This function reads str from a solid memory block
t_str *read_str_from_memory(void *vp)
{
t_str *str;
register unsigned int name_len, size_r, natoms, nedges;
size_r=*((unsigned int*)vp),                       vp+=sizeof(unsigned int);
natoms=*((unsigned int*)vp),                       vp+=sizeof(unsigned int);
nedges=*((unsigned int*)vp),                       vp+=sizeof(unsigned int);
name_len=*((unsigned int*)vp),                     vp+=sizeof(unsigned int);
if (!(str=(t_str*)alloc_solid_str(name_len,size_r,natoms,nedges))) return FALSE; 
memcpy(str->name,vp,sizeof(char)*name_len),        vp+=sizeof(char)*name_len,   str->name[name_len]='\0';
if (!(size_r)) { *((unsigned int*)&str->rsize)=*((unsigned int*)vp), vp+=sizeof(unsigned int); }
else
  {
  memcpy(str->ress->list,vp,sizeof(unsigned int)*str->ress->size), vp+=sizeof(unsigned int)*str->ress->size;
  memcpy(str->rsize,vp,sizeof(unsigned int)*(str->ress->size+1)), vp+=sizeof(unsigned int)*(str->ress->size+1);
  }
str->start_rid=*((unsigned int*)vp),               vp+=sizeof(unsigned int);
memcpy(str->anames,vp,sizeof(int)*str->natoms),    vp+=sizeof(int)*str->natoms;
memcpy(str->a,vp,sizeof(char)*str->natoms),        vp+=sizeof(char)*str->natoms;
memcpy(str->edges,vp,sizeof(t_edge)*str->nedges),  vp+=sizeof(t_edge)*str->nedges;
memcpy(str->r,vp,sizeof(t_vec)*str->natoms),       vp+=sizeof(t_vec)*str->natoms;
return str;
}
//This function writes str into solid memory block
void *write_str_to_memory(size_t *size,t_str *str)
{
void *vp;
register unsigned int _i;
_i=0; if ( (str->name)&&(*str->name) ) { while ( (str->name[_i])) _i++; }
if (!(str->ress))
  {
  *size=sizeof(unsigned int)*0x4+sizeof(char)*_i+sizeof(unsigned int)+sizeof(unsigned int)+(sizeof(char)+sizeof(int)+sizeof(t_vec))*str->natoms+sizeof(t_edge)*str->nedges;
  if (!(vp=(void*)malloc(*size)))
    { ylib_errno=YERROR_MEMORY; return FALSE; }
  *((unsigned int*)vp)=1,                             vp+=sizeof(unsigned int);
  *((unsigned int*)vp)=str->natoms,                   vp+=sizeof(unsigned int);
  *((unsigned int*)vp)=str->nedges,                   vp+=sizeof(unsigned int);
  *((unsigned int*)vp)=_i,                            vp+=sizeof(unsigned int);
  memcpy(vp,str->name,sizeof(char)*_i),               vp+=sizeof(char)*_i;
  *((unsigned int*)vp)=*((unsigned int*)&str->rsize), vp+=sizeof(unsigned int);
  }
else
  {
  *size=sizeof(unsigned int)*0x4+sizeof(char)*_i+sizeof(unsigned int)*(0x2*str->ress->size+0x1)+sizeof(unsigned int)+(sizeof(char)+sizeof(int)+sizeof(t_vec))*str->natoms+sizeof(t_edge)*str->nedges;
  if (!(vp=(void*)malloc(*size)))
    { ylib_errno=YERROR_MEMORY; return FALSE; }
  *((unsigned int*)vp)=str->ress->size,               vp+=sizeof(unsigned int);
  *((unsigned int*)vp)=str->natoms,                   vp+=sizeof(unsigned int);
  *((unsigned int*)vp)=str->nedges,                   vp+=sizeof(unsigned int);
  *((unsigned int*)vp)=_i,                            vp+=sizeof(unsigned int);
  memcpy(vp,str->name,sizeof(char)*_i),               vp+=sizeof(char)*_i;
  memcpy(vp,str->ress->list,sizeof(unsigned int)*str->ress->size), vp+=sizeof(unsigned int)*str->ress->size;
  memcpy(vp,str->rsize,sizeof(unsigned int)*(str->ress->size+1)), vp+=sizeof(unsigned int)*(str->ress->size+1);
  }
*((unsigned int*)vp)=str->start_rid,                  vp+=sizeof(unsigned int);
memcpy(vp,str->anames,sizeof(int)*str->natoms),       vp+=sizeof(int)*str->natoms;
memcpy(vp,str->a,sizeof(char)*str->natoms),           vp+=sizeof(char)*str->natoms;
memcpy(vp,str->edges,sizeof(t_edge)*str->nedges),     vp+=sizeof(t_edge)*str->nedges;
memcpy(vp,str->r,sizeof(t_vec)*str->natoms),          vp+=sizeof(t_vec)*str->natoms;
return vp-*size;
}

//This function reads ymol
t_mol *read_ymol(FILE *in)
{
unsigned int i;
t_mol *mol;
if (!(mol=(t_mol*)calloc(sizeof(t_mol),0x1))) { ylib_errno=YERROR_MEMORY; return FALSE; }
if (fread(&i,sizeof(unsigned int),0x1,in)!=0x1) { IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if (i!=Y_MAGIC) { DATA_CONSISTMENT_ERROR: free_mol(mol); ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }  //Check YMAGIC
if (fread(&i,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
//                YSTR part
if (!(i)) mol->name=0x0;
else
  { 
  if (!(mol->name=(char*)malloc(sizeof(char)*(i+1)))) { MEMORY_ERROR: free_mol(mol); ylib_errno=YERROR_MEMORY; return FALSE; }
  if (fread(mol->name,sizeof(char),i,in)!=i) goto IO_ERROR;
  else mol->name[i]='\0';
  }
if (!(mol->ress=read_list(in))) { EXIT: free_mol(mol); return FALSE; }
if (!mol->ress->size) goto DATA_CONSISTMENT_ERROR;
if (!(mol->rsize=(unsigned int*)malloc(sizeof(unsigned int)*(mol->ress->size+1)))) goto MEMORY_ERROR;
if (fread(mol->rsize,sizeof(unsigned int),mol->ress->size+1,in)!=mol->ress->size+1) goto IO_ERROR;
if (fread(&mol->start_rid,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if (fread(&mol->natoms,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if (!(mol->natoms)) goto DATA_CONSISTMENT_ERROR;
if (!(mol->anames=(char (*)[sizeof(int)])malloc(sizeof(int)*mol->natoms))) goto MEMORY_ERROR;
if (!(mol->a=(char*)malloc(sizeof(char)*mol->natoms))) goto MEMORY_ERROR;
if (!(mol->r=(t_vec*)malloc(sizeof(t_vec)*mol->natoms))) goto MEMORY_ERROR;
if (fread(mol->anames,sizeof(int),mol->natoms,in)!=mol->natoms) goto IO_ERROR;
if (fread(mol->a,sizeof(char),mol->natoms,in)!=mol->natoms) goto IO_ERROR;
if (fread(mol->r,sizeof(t_vec),mol->natoms,in)!=mol->natoms) goto IO_ERROR;
if (fread(&mol->nedges,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (mol->nedges))
  {
  if (!(mol->edges=(t_edge*)malloc(sizeof(t_edge)*mol->nedges))) goto MEMORY_ERROR;
  if (fread(mol->edges,sizeof(t_edge),mol->nedges,in)!=mol->nedges) goto IO_ERROR;
  }
//                Topology part
if (fread(&mol->size_ar,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if (!(mol->cycles=read_clist(in))) goto EXIT;
if (!(mol->anchors=read_clist(in))) goto EXIT;
if (!(mol->anchors->size)) goto DATA_CONSISTMENT_ERROR;
if (fread(&mol->naedges,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (mol->naedges))
  {
  if (!(mol->aedges=(t_edge*)malloc(sizeof(t_edge)*mol->naedges))) goto MEMORY_ERROR;
  if (fread(mol->aedges,sizeof(t_edge),mol->naedges,in)!=mol->naedges) goto IO_ERROR;
  }
if (!(mol->ytypes=(unsigned int*)malloc(sizeof(unsigned int)*mol->natoms))) goto MEMORY_ERROR;
if (fread(mol->ytypes,sizeof(unsigned int),mol->natoms,in)!=mol->natoms) goto IO_ERROR;
//                Molecular Mechanics FF part
//                Electron distribution functions 
if (fread(&mol->nvatoms,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (mol->nvatoms))
  {
  if (!(mol->vatoms=(int*)malloc(sizeof(int)*mol->nvatoms))) goto MEMORY_ERROR;
  if (fread(mol->vatoms,sizeof(int),mol->nvatoms,in)!=mol->nvatoms) goto IO_ERROR;
  }
if (!(mol->engs=(double*)malloc(sizeof(double)*(mol->natoms+mol->nvatoms)))) goto MEMORY_ERROR;
if (fread(mol->engs,sizeof(double),mol->natoms+mol->nvatoms,in)!=mol->natoms+mol->nvatoms) goto IO_ERROR;
if (fread(&mol->nvedges,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (mol->nvedges))
  {
  if (!(mol->vedges=(unsigned int(*)[2])malloc(sizeof(unsigned int)*2*mol->nvedges))) goto MEMORY_ERROR;
  if (fread(mol->vedges,sizeof(unsigned int)*2,mol->nvedges,in)!=mol->nvedges) goto IO_ERROR;
  }
if ( (mol->nedges+mol->nvedges))
  {
  if (!(mol->hrds=(double*)malloc(sizeof(double)*(mol->nedges+mol->nvedges)))) goto MEMORY_ERROR;
  if (fread(mol->hrds,sizeof(double),mol->nedges+mol->nvedges,in)!=mol->nedges+mol->nvedges) goto IO_ERROR; 
  }
if (fread(&mol->cmtype,sizeof(char),0x1,in)!=0x1) goto IO_ERROR;
     if (mol->cmtype==+1) {      if (!(mol->C.dL=(t_dmatrix*)read_tdmatrix(in,(char*)&i))) goto IO_ERROR; 
                            else if ( (*((char*)&i)!='L')&&(*((char*)&i)!='l') ) goto DATA_CONSISTMENT_ERROR; }
else if (mol->cmtype==-1) { if (!(mol->C.sL=(t_smatrix*)read_smatrix(in)))  goto IO_ERROR; }
else mol->C.dL=0x0;
if (!(mol->charges=(double*)malloc(sizeof(double)*(mol->natoms+mol->nvatoms)))) goto MEMORY_ERROR;
if (fread(mol->charges,sizeof(double),mol->natoms+mol->nvatoms,in)!=mol->natoms+mol->nvatoms) goto IO_ERROR;
//                FF Parameters
if (fread(&i,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if (i==(unsigned int)-1) mol->fftypes=mol->ytypes; //NB! Default assignement
else
  {
  if (!(mol->fftypes=(unsigned int*)malloc(sizeof(unsigned int)*mol->natoms))) goto MEMORY_ERROR;
  else mol->fftypes[0]=i;
  if (mol->natoms>1)
    {
    if (fread(&mol->fftypes[1],sizeof(unsigned int),mol->natoms-1,in)!=mol->natoms-1) goto IO_ERROR;
    } 
  }
if (fread(&mol->size_b,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (mol->size_b))
  {
  if (!(mol->ff_b=(t_ff_b*)malloc(sizeof(t_ff_b)*mol->size_b))) goto MEMORY_ERROR;
  if (fread(mol->ff_b,sizeof(t_ff_b),mol->size_b,in)!=mol->size_b) goto IO_ERROR;
  }
if (fread(&mol->size_g,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (mol->size_g))
  {
  if (!(mol->ff_g=(t_ff_g*)malloc(sizeof(t_ff_g)*mol->size_g))) goto MEMORY_ERROR;
  if (fread(mol->ff_g,sizeof(t_ff_g),mol->size_g,in)!=mol->size_g) goto IO_ERROR;
  }
if (fread(&mol->size_i,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (mol->size_i))
  {
  if (!(mol->ff_i=(t_ff_i*)malloc(sizeof(t_ff_i)*mol->size_i))) goto MEMORY_ERROR;
  if (fread(mol->ff_i,sizeof(t_ff_i),mol->size_i,in)!=mol->size_i) goto IO_ERROR;
  }
if (fread(&mol->size_d,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (mol->size_d))
  {
  if (!(mol->ff_d=(t_ff_d*)malloc(sizeof(t_ff_d)*mol->size_d))) goto MEMORY_ERROR;
  if (fread(mol->ff_d,sizeof(t_ff_d),mol->size_d,in)!=mol->size_d) goto IO_ERROR;
  }
if (fread(&i,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (i))
  {
  if (!(mol->excl=(t_list*)malloc(sizeof(t_list)*mol->natoms+sizeof(unsigned int)*i))) goto MEMORY_ERROR;
  else
    {
    mol->excl[0].list=(unsigned int*)((void*)mol->excl+sizeof(t_list)*mol->natoms); 
    if (fread(mol->excl[0].list,sizeof(unsigned int),i,in)!=i) goto IO_ERROR;
    }
  for (i=0; i<mol->natoms; i++)
    if (fread(&mol->excl[i].size,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR; 
  for (i=1; i<mol->natoms; i++) 
    mol->excl[i].list=mol->excl[i-1].list+mol->excl[i-1].size;
  }
if (fread(&mol->size_c,sizeof(unsigned int),0x1,in)!=0x1) goto IO_ERROR;
if ( (mol->size_c))
  {
  if (!(mol->ff_c=(t_ff_b*)malloc(sizeof(t_ff_b)*mol->size_c))) goto MEMORY_ERROR;
  if (fread(mol->ff_c,sizeof(t_ff_b),mol->size_c,in)!=mol->size_c) goto IO_ERROR;
  }
return mol;
}

//This function writtes ymol
char write_ymol(FILE *out,t_mol *mol)
{
unsigned int i=Y_MAGIC;
//                YSTR part
if (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1) { IO_ERROR: ylib_errno=YERROR_IO; return FALSE; } //write YMAGIC
i=0;
if (!(mol->name)) { if (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR; } //write zero- name_len
else
  {
  while ( (mol->name[i])) i++;
  if ( (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(mol->name,sizeof(char),i,out)!=i) ) goto IO_ERROR; // write name_len and name 
  }
if ( (!(write_list(out,mol->ress)))||(fwrite(mol->rsize,sizeof(unsigned int),mol->ress->size+1,out)!=mol->ress->size+1) ) goto IO_ERROR; //write residues
if (fwrite(&mol->start_rid,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR; //Save starting residue id of the sequence
if (fwrite(&mol->natoms,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR;
if (fwrite(mol->anames,sizeof(int),mol->natoms,out)!=mol->natoms) goto IO_ERROR; //write atoms names
if (fwrite(mol->a,sizeof(char),mol->natoms,out)!=mol->natoms) goto IO_ERROR; //write atoms names
if (fwrite(mol->r,sizeof(t_vec),mol->natoms,out)!=mol->natoms) goto  IO_ERROR; //write coords
if (fwrite(&mol->nedges,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR;
if ( (mol->nedges)) { if (fwrite(mol->edges,sizeof(t_edge),mol->nedges,out)!=mol->nedges) goto IO_ERROR; } //write edges
//                Topology part
if (fwrite(&mol->size_ar,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR; //write amount of aromatic cycles
if (!(write_clist(out,mol->cycles))) goto IO_ERROR; //write cycles
if (!(write_clist(out,mol->anchors))) goto IO_ERROR; //write anchors 
if (fwrite(&mol->naedges,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR;
if ( (mol->naedges)) { if (fwrite(mol->aedges,sizeof(t_edge),mol->naedges,out)!=mol->naedges) goto IO_ERROR; } //write aedges
if (fwrite(mol->ytypes,sizeof(unsigned int),mol->natoms,out)!=mol->natoms) goto IO_ERROR; //write atoms names
//                Molecular Mechanics FF part
//                Electron distribution functions 
if (fwrite(&mol->nvatoms,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR;
if ( (mol->nvatoms)) { if (fwrite(mol->vatoms,sizeof(int),mol->nvatoms,out)!=mol->nvatoms) goto IO_ERROR; }
if (fwrite(mol->engs,sizeof(double),mol->natoms+mol->nvatoms,out)!=mol->natoms+mol->nvatoms) goto IO_ERROR; //write (v)atoms
if (fwrite(&mol->nvedges,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR;
if ( (mol->nvedges)) { if (fwrite(mol->vedges,sizeof(unsigned int)*2,mol->nvedges,out)!=mol->nvedges) goto IO_ERROR; }
if ( (mol->nedges+mol->nvedges)) { if (fwrite(mol->hrds,sizeof(double),mol->nedges+mol->nvedges,out)!=mol->nedges+mol->nvedges) goto IO_ERROR; } //write (v)bonds
if (fwrite(&mol->cmtype,sizeof(char),0x1,out)!=0x1) goto IO_ERROR;
     if (mol->cmtype==+1) { if (!(write_tdmatrix(out,'L',mol->C.dL))) goto IO_ERROR; }
else if (mol->cmtype==-1) { if (!(write_smatrix(out,mol->C.sL)))  goto IO_ERROR; }
if (fwrite(mol->charges,sizeof(double),mol->natoms+mol->nvatoms,out)!=mol->natoms+mol->nvatoms) goto IO_ERROR; //write cmantrix and charges
//                FF Parameters
if (mol->fftypes==mol->ytypes) //NB! Default assignement
  { i=(unsigned int)-1; if (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR; }
else { if (fwrite(&mol->fftypes[1],sizeof(unsigned int),mol->natoms,out)!=mol->natoms) goto IO_ERROR; }
if ( (fwrite(&mol->size_b,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(mol->ff_b,sizeof(t_ff_b),mol->size_b,out)!=mol->size_b) ) goto IO_ERROR; //write ff_b
if ( (fwrite(&mol->size_g,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(mol->ff_g,sizeof(t_ff_g),mol->size_g,out)!=mol->size_g) ) goto IO_ERROR; //write ff_g
if ( (fwrite(&mol->size_i,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(mol->ff_i,sizeof(t_ff_i),mol->size_i,out)!=mol->size_i) ) goto IO_ERROR; //write ff_i
if ( (fwrite(&mol->size_d,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(mol->ff_d,sizeof(t_ff_d),mol->size_d,out)!=mol->size_d) ) goto IO_ERROR; //write ff_d
i=(unsigned int)(mol->excl[mol->natoms-1].list-mol->excl[0].list)+mol->excl[mol->natoms-1].size;
if (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR;
if ( (i)) //Write exclusions
  {
  if (fwrite(mol->excl[0].list,sizeof(unsigned int),i,out)!=i) goto IO_ERROR;
  for (i=0; i<mol->natoms; i++)
    if (fwrite(&mol->excl[i].size,sizeof(unsigned int),0x1,out)!=0x1) goto IO_ERROR; 
  }
if ( (fwrite(&mol->size_c,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(mol->ff_c,sizeof(t_ff_b),mol->size_c,out)!=mol->size_c) ) goto IO_ERROR; //write ff_c
return TRUE;
}

//START OF DBSTR PART

//Tips: dbstr is solid - free(dbstr) is OK
//This function allocates solid dbstr structure
t_dbstr *alloc_solid_dbstr(unsigned int size_r,unsigned int natoms,unsigned int nedges)
{
t_dbstr *dbstr;
if (size_r==1)
  {
  if (!(dbstr=(t_dbstr*)malloc(sizeof(t_dbstr)+sizeof(int)*(natoms*sizeof(char)/sizeof(int)+(( (natoms*sizeof(char)%sizeof(int))) ? 1 : 0))+sizeof(t_edge)*nedges))) { ylib_errno=YERROR_MEMORY; return FALSE; } 
  else { dbstr->a=(char*)((void*)dbstr+sizeof(t_dbstr)), dbstr->edges=(t_edge*)((void*)dbstr->a+sizeof(int)*(natoms*sizeof(char)/sizeof(int)+(( (natoms*sizeof(char)%sizeof(int))) ? 1 : 0))); }
  dbstr->ress=0x0, dbstr->natoms=natoms, dbstr->nedges=nedges;
  }
else
  {
  if (!(dbstr=(t_dbstr*)malloc(sizeof(t_dbstr)+sizeof(int)*(natoms*sizeof(char)/sizeof(int)+(( (natoms*sizeof(char)%sizeof(int))) ? 1 : 0))+sizeof(t_edge)*nedges+sizeof(t_list)+sizeof(unsigned int)*(2*size_r+1)))) { ylib_errno=YERROR_MEMORY; return FALSE; } 
  else { dbstr->a=(char*)((void*)dbstr+sizeof(t_dbstr)), dbstr->edges=(t_edge*)((void*)dbstr->a+sizeof(int)*(natoms*sizeof(char)/sizeof(int)+(( (natoms*sizeof(char)%sizeof(int))) ? 1 : 0))), dbstr->ress=(t_list*)((void*)dbstr->edges+sizeof(t_edge)*dbstr->nedges), dbstr->ress->list=(unsigned int*)((void*)dbstr->ress+sizeof(t_list)), dbstr->rsize=dbstr->ress->list+size_r; }
  dbstr->ress->size=size_r, dbstr->natoms=natoms, dbstr->nedges=nedges;
  }
return dbstr;
}
//This function reads dbstr structure from file [edited from y_mol.c]
t_dbstr *read_dbstr(FILE *in)
{
unsigned int i, size_r, natoms, nedges;
t_dbstr *dbstr;
//Read constants
if ( (fread(&i,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&size_r,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&natoms,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&nedges,sizeof(unsigned int),0x1,in)!=0x1) ) { LABEL_IO_ERROR_0: ylib_errno=YERROR_IO; return FALSE; }
if ( (i!=Y_MAGIC)||(!size_r)||(!natoms) ) { ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
//Allocate solid structure in memory and upload variable ress + rsize
if (!(dbstr=alloc_solid_dbstr(size_r,natoms,nedges))) return FALSE;
if (size_r==1)
  {
  if (fread(&dbstr->rsize,sizeof(unsigned int),0x1,in)!=0x1) { LABEL_IO_ERROR_1: free(dbstr); goto LABEL_IO_ERROR_0; }
  }
else
  {
  if ( (fread(dbstr->ress->list,sizeof(unsigned int),size_r,in)!=size_r)||(fread(dbstr->rsize,sizeof(unsigned int),size_r,in)!=size_r) ) goto LABEL_IO_ERROR_1;
  }
//Upload atoms
if (fread(dbstr->a,sizeof(char),natoms,in)!=natoms) goto LABEL_IO_ERROR_1;
//Upload edges
if (dbstr->nedges) { if (fread(dbstr->edges,sizeof(t_edge),dbstr->nedges,in)!=dbstr->nedges) goto LABEL_IO_ERROR_1; } else dbstr->edges=0x0;
return dbstr;
}
//This function write dbstr to hdd (including the header) [edited from y_mol.c]
char write_dbstr(t_dbstr *dbstr,FILE *out)
{
unsigned int i;
i=Y_MAGIC; if (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1) { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; } //Write YMAGIC
if ( ( (dbstr->ress))&&( (dbstr->ress->size))&&(dbstr->ress->size!=1) )
  {//Full residues 
  if ( (fwrite(&dbstr->ress->size,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&dbstr->natoms,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&dbstr->nedges,sizeof(unsigned int),0x1,out)!=0x1)||
       (fwrite(dbstr->ress->list,sizeof(unsigned int),dbstr->ress->size,out)!=dbstr->ress->size)||(fwrite(dbstr->rsize,sizeof(unsigned int),dbstr->ress->size+1,out)!=dbstr->ress->size+1) ) goto LABEL_IO_ERROR;
  }
else
  {//Compressed residues
  i=1;
  if ( (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&dbstr->natoms,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&dbstr->nedges,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&dbstr->rsize,sizeof(unsigned int),0x1,out)!=0x1) ) goto LABEL_IO_ERROR;
  }
//Write atoms
if (fwrite(dbstr->a,sizeof(char),dbstr->natoms,out)!=dbstr->natoms) goto LABEL_IO_ERROR; 
//Write edges
if (fwrite(dbstr->edges,sizeof(t_edge),dbstr->nedges,out)!=dbstr->nedges) goto LABEL_IO_ERROR; 
return TRUE;
}

//This function does str_to_dbstr conversion
t_dbstr *str_to_dbstr(t_str *str)
{
t_dbstr *dbstr;
if (!str->ress)
  {
  if (!(dbstr=alloc_solid_dbstr(1,str->natoms,str->nedges))) return FALSE; 
  else { dbstr->ress=0x0, *((unsigned int*)&dbstr->rsize)=*((unsigned int*)&str->rsize); }
  }
else
  {
  if (!(dbstr=alloc_solid_dbstr(str->ress->size,str->natoms,str->nedges))) return FALSE;
  else { dbstr->ress->size=str->ress->size; }
  memcpy(dbstr->ress->list,str->ress->list,sizeof(unsigned int)*dbstr->ress->size);
  memcpy(dbstr->rsize,str->rsize,sizeof(unsigned int)*(dbstr->ress->size+1));
  }
dbstr->natoms=str->natoms, dbstr->nedges=str->nedges;
memcpy(dbstr->a,str->a,sizeof(char)*dbstr->natoms);
memcpy(dbstr->edges,str->edges,sizeof(t_edge)*dbstr->nedges);
return dbstr;
}
//This function does dbstr_to_str conversion (since dbstr is smaller not all data is recovered)
t_str *dbstr_to_str(t_vec  *r,t_dbstr *dbstr)
{
register unsigned int _l;
t_str *str;
if (!(dbstr->ress)) 
  {
  if (!(str=alloc_solid_str(sizeof(unsigned int)+1,0x1,dbstr->natoms,dbstr->nedges))) return FALSE; 
  else { str->ress=0x0, *((unsigned int*)&str->rsize)=*((unsigned int*)&dbstr->rsize); }
  }
else 
  {
  if (!(str=alloc_solid_str(sizeof(unsigned int)+1,dbstr->ress->size,dbstr->natoms,dbstr->nedges))) return FALSE; 
  str->ress->size=dbstr->ress->size;
  memcpy(str->ress->list,dbstr->ress->list,sizeof(unsigned int)*dbstr->ress->size); 
  memcpy(str->rsize,dbstr->rsize,sizeof(unsigned int)*(dbstr->ress->size+1)); 
  }
*((int*)str->name)=dbstr->id, str->start_rid=1;
str->natoms=dbstr->natoms;
memcpy(str->a,dbstr->a,sizeof(char)*dbstr->natoms);
memcpy(str->r,r,sizeof(t_vec)*dbstr->natoms);
_l=dbstr->natoms; while (_l--) *((int*)str->anames[_l])=chemid_to_name(dbstr->a[_l]);
str->nedges=dbstr->nedges;
memcpy(str->edges,dbstr->edges,sizeof(t_edge)*dbstr->nedges);
return str;
}

//END OF DBSTR PART


/************************************      S T R      C O M P I L A T I O N     P A R T        *****************************************************/





//We DON'T use neighbors here becuse the molecular graph is rebuilt upon the callings 

//This function is a service module for the protonate_str
//It returns the id of the atom to be inserted
inline unsigned int str_add_hydrogen(register unsigned int atom_id,register t_str *str)
{
register unsigned int _i, _j;
void *vp;

//Stage 1. Prepare memory
if (!(str->natoms%0xFF))
  {
  if (!(vp=realloc(str->a,(_j=str->natoms+0xFF)*sizeof(char))))
    {
    LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; 
    return (unsigned int)-1;
    }
  else str->a=(char*)vp;
  if (!(vp=(void*)realloc(str->anames,_j*sizeof(int)))) goto LABEL_MEMORY_ERROR;
  else str->anames=(char (*)[sizeof(int)])vp; 
  if (!(vp=realloc(str->r,_j*sizeof(t_vec)))) goto LABEL_MEMORY_ERROR;
  else str->r=(t_vec*)vp; 
  }
if (!(str->nedges%0xFF))
  {
  if (!(vp=realloc(str->edges,(str->nedges+0xFF)*sizeof(t_edge)))) goto LABEL_MEMORY_ERROR;
  else str->edges=(t_edge*)vp; 
  }

//Stage 2. Locate the position to paste the new hydrogen
if ( (str->ress))
  {//Insert hydrogen into the body of molecule at the end of parent residue - a lot of work to be done
  _j=0; while (str->rsize[_j]<=atom_id) _j++;
  //Shift residues 
  _i=str->rsize[_j]; do { str->rsize[_j]++; } while (str->ress->size!=_j++) ;
  //Shift atoms
  _j=str->natoms; 
  do{ 
    str->a[_j]=str->a[_j-1], str->r[_j].i=str->r[_j-1].i, str->r[_j].j=str->r[_j-1].j, str->r[_j].k=str->r[_j-1].k;
    str->anames[_j][0]=str->anames[_j-1][0], str->anames[_j][1]=str->anames[_j-1][1], str->anames[_j][2]=str->anames[_j-1][2], str->anames[_j][3]=str->anames[_j-1][3];
    } while (_i!=_j--);
  //Update edges
  _j=str->nedges;
  while (_j--) 
    {
    if (str->edges[_j].vertice[0]>=_i) str->edges[_j].vertice[0]++;
    if (str->edges[_j].vertice[1]>=_i) str->edges[_j].vertice[1]++;
    }
  }
else _i=str->natoms; //Insert hydrogen in the end of the molecule - no much job to do 

//Stage 3. Paste the new hydrogen
str->a[_i]=CHEM_ATOM_TYPE_HYDROGEN;
str->anames[_i][0]='H', str->anames[_i][1]=str->anames[_i][2]=str->anames[_i][3]=' ';
str->r[_i].i=str->r[_i].j=str->r[_i].k=(double)NAN;
str->edges[str->nedges].vertice[0]=atom_id, str->edges[str->nedges].vertice[1]=_i, str->edges[str->nedges].type='1';
str->natoms++, str->nedges++;

return _i;
}

//This function deletes a hydrohen from str.
//It returns the atom_id that was deleted or (unsigned int)-1 on failure 
inline unsigned int str_del_hydrogen(register unsigned int atom_id,register t_str *str)
{
register unsigned int _i;
if (atom_id>=str->natoms) { ylib_errno=YERROR_IMPOSSIBLE; return (unsigned int)-1; }
//shift residues
if ( (str->ress)) { _i=0; while (str->rsize[_i]<=atom_id) _i++; do { str->rsize[_i]--; } while (str->ress->size!=_i++); }
str->natoms--;
//shift atoms
_i=atom_id; 
while (_i!=str->natoms)
  {
  str->a[_i]=str->a[_i+1], str->r[_i].i=str->r[_i+1].i, str->r[_i].j=str->r[_i+1].j, str->r[_i].k=str->r[_i+1].k;
  str->anames[_i][0]=str->anames[_i+1][0], str->anames[_i][1]=str->anames[_i+1][1], str->anames[_i][2]=str->anames[_i+1][2], str->anames[_i][3]=str->anames[_i+1][3], _i++;
  } 
//update edges
_i=str->nedges; 
while (_i--) 
  if ( (str->edges[_i].vertice[0]==atom_id)||(str->edges[_i].vertice[1]==atom_id) )
    { str->nedges--, str->edges[_i].vertice[0]=str->edges[str->nedges].vertice[0], str->edges[_i].vertice[1]=str->edges[str->nedges].vertice[1], str->edges[_i].type=str->edges[str->nedges].type; }
  else 
    {
    if (str->edges[_i].vertice[0]>atom_id) str->edges[_i].vertice[0]--;    
    if (str->edges[_i].vertice[1]>atom_id) str->edges[_i].vertice[1]--;    
    }
return atom_id;
}


//This is a primitive function to full-fill missing hydogens
unsigned int protonate_str(register t_str *str)
{
register unsigned int _i, count;
char *nv;

//Stage 1. Prepare memory and gather bonding statistics
if (!(nv=(char*)calloc(str->natoms,sizeof(char)))) { ylib_errno=YERROR_MEMORY; return (unsigned int)-1; }
else 
  {
  _i=str->nedges;
  while (_i--)
    switch (str->edges[_i].type)
      {
      case '1' : { nv[str->edges[_i].vertice[0]]++;  nv[str->edges[_i].vertice[1]]++;  break; }
      case '2' : { nv[str->edges[_i].vertice[0]]+=2; nv[str->edges[_i].vertice[1]]+=2; break; }
      case '3' : { nv[str->edges[_i].vertice[0]]+=3, nv[str->edges[_i].vertice[1]]+=3; break; }
      case 'm' : { nv[str->edges[_i].vertice[0]]++,  nv[str->edges[_i].vertice[1]]++;  break; }
      case 'a' : { nv[str->edges[_i].vertice[0]]++,  nv[str->edges[_i].vertice[1]]++;  break; } 
      default  : { free(nv); nv=0x0; ylib_errno=YERROR_IMPOSSIBLE;  return (unsigned int)-1; } 
      }
  }

//Stage 2. Add the formal hydrogens
count=0, _i=str->natoms;
while (_i--)
  {
  switch (str->a[_i])
    {
    case CHEM_ATOM_TYPE_HYDROGEN : {
      while (nv[_i]<1) if (str_add_hydrogen(_i,str)==(unsigned int)-1) { LABEL_ERROR_PROTON_ADDING: free(nv); nv=0x0; return FALSE; } else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_CARBON   : {
      while (nv[_i]<1) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_NITROGEN : {
      while (nv[_i]<3) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_OXYGEN   : {
      while (nv[_i]<2) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_FLUORINE : {
      while (nv[_i]<1) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_SILICON  : {
      while (nv[_i]<4) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_PHOSPHOR : {
      while (nv[_i]<3) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_SULFUR   : {
      while (nv[_i]<2) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_CHLORINE : {
      while (nv[_i]<1) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_BROMINE  : {
      while (nv[_i]<1) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    case CHEM_ATOM_TYPE_IODINE   : {
      while (nv[_i]<1) if (str_add_hydrogen(_i,str)==(unsigned int)-1) goto LABEL_ERROR_PROTON_ADDING; else { nv[_i]++, count++; }
      break;
      }
    default                      : ; //Do nothing by default
    }
  }

//Exit
free(nv); nv=0x0;
return count;
}

//This function edits R-(N=O)-OH -> R-N(=O)2
unsigned int update_NO2(char (*order)[4],t_str *str)
{
register unsigned int _i, _j, _k, _t, _n1, count=0;

//Search for NO2 group
_i=str->nedges;
while(_i--) 
  if (str->edges[_i].type=='2')
    {//case (N=X)-Y: NO2 if ( X=='O' and Y == OH )
    if ( (str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_NITROGEN)&&(str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_OXYGEN) )
      { _t=str->edges[_i].vertice[0], str->edges[_i].vertice[0]=str->edges[_i].vertice[1], str->edges[_i].vertice[1]=_t; goto ANALYZE_NO; }
    if ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_NITROGEN)&&(str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN) )
      {
      ANALYZE_NO: ; 
      if ( ( (order[str->edges[_i].vertice[0]][3]))||( (order[str->edges[_i].vertice[1]][3]))||
             (order[str->edges[_i].vertice[1]][2]+order[str->edges[_i].vertice[1]][1]+order[str->edges[_i].vertice[1]][0]!=3)||
             (order[str->edges[_i].vertice[1]][2]+order[str->edges[_i].vertice[1]][1]+order[str->edges[_i].vertice[1]][0]!=1) )
        goto NEXT_I;
      //Switch O=N-O or O=N=O 
           if (order[str->edges[_i].vertice[0]][2]==1)
             {
             //Check for N-O 
             _n1=2, _j=str->nedges;
             while (_j--) 
               {
               if ( (str->edges[_j].vertice[1]==str->edges[_i].vertice[0])&&(_i!=_j) )
                 { _t=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=str->edges[_j].vertice[1], str->edges[_j].vertice[1]=_t; goto ANALYZE_ONOH; }
               if ( (str->edges[_j].vertice[0]==str->edges[_i].vertice[0])&&(_i!=_j) )
                 {
                 ANALYZE_ONOH: ;
                 if ( (!(order[str->edges[_j].vertice[1]][3]))&&(!(order[str->edges[_j].vertice[1]][2]))&&(order[str->edges[_j].vertice[1]][1]+order[str->edges[_j].vertice[1]][0]==2)&&(str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN) )
                   {
                   //Check for O-H
                   _k=str->nedges;
                   while (_k--)
                     { 
                     if ( (str->edges[_k].vertice[1]==str->edges[_j].vertice[1])&&(_i!=_j) )
                       { _t=str->edges[_k].vertice[0], str->edges[_k].vertice[0]=str->edges[_k].vertice[1], str->edges[_k].vertice[1]=_t; goto ANALYZE_OH; }
                     if ( (str->edges[_k].vertice[0]==str->edges[_j].vertice[1])&&(_k!=_j) )
                       {
                       ANALYZE_OH: ; 
                       if ( (!(order[str->edges[_k].vertice[1]][3]))&&(!(order[str->edges[_k].vertice[1]][2]))&&(order[str->edges[_k].vertice[1]][1]+order[str->edges[_k].vertice[1]][0]==1)&&(str->a[str->edges[_k].vertice[1]]==CHEM_ATOM_TYPE_HYDROGEN) ) 
                         {
                         if (str->edges[_j].type=='1') { order[str->edges[_j].vertice[0]][0]--, order[str->edges[_j].vertice[1]][0]--; }
                         else                          { order[str->edges[_j].vertice[0]][1]--, order[str->edges[_j].vertice[1]][1]--; }
                         order[str->edges[_j].vertice[0]][2]++, order[str->edges[_j].vertice[1]][2]++, str->edges[_j].type='2';                         
                         if (str->edges[_k].type!='1') order[str->edges[_j].vertice[1]][1]--; else order[str->edges[_j].vertice[1]][0]--;
                         if ((_k=str_del_hydrogen(str->edges[_k].vertice[1],str))==(unsigned int)-1) return (unsigned int)-1;
                         else { while (_k!=str->natoms) { order[_k][0]=order[_k+1][0], order[_k][1]=order[_k+1][1], order[_k][2]=order[_k+1][2], order[_k][3]=order[_k+1][3], _k++; } }
                         count++;
                         }
                       break; 
                       }
                     }
                   }
                 if (!(--_n1)) break;
                 }
               }   
             }
      else if (order[str->edges[_i].vertice[0]][2]==2)
             {
             //Check & Skip O=N=O
             _j=str->nedges;
             while (_j--) 
               if ( (str->edges[_j].type=='2')&&(_i!=_j) )
                 {
                 if ( (str->edges[_j].vertice[1]==str->edges[_i].vertice[0])&&(str->a[str->edges[_j].vertice[0]]==CHEM_ATOM_TYPE_OXYGEN) )
                   { _t=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=str->edges[_j].vertice[1], str->edges[_j].vertice[1]=_t; goto ANALYZE_ONO; } 
                 if ( (str->edges[_j].vertice[0]==str->edges[_i].vertice[0])&&(str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN) )
                   {
                   ANALYZE_ONO: ;
                   if ( (!(order[str->edges[_j].vertice[1]][3]))&&(order[str->edges[_j].vertice[1]][2]==1)&&(!(order[str->edges[_j].vertice[1]][1]))&&(!(order[str->edges[_j].vertice[1]][0])) )
                     { 
                     order[str->edges[_j].vertice[1]][2]++, order[str->edges[_j].vertice[0]][1]++;
                     order[str->edges[_j].vertice[1]][2]=0, order[str->edges[_j].vertice[1]][1]=1;
                     str->edges[_j].type='a', count++;
                     }
                   break;
                   }
                 }
             }
      }
    NEXT_I: ;
    }
return count;
}

//This function edits amide bond: R1-(C-X-H)=N-R2 -> R1-(C=X)-(N-H)-R2 where X={O,S}
unsigned int update_amides(char (*order)[4],t_str *str)
{
register unsigned int _i, _j, _k, _t, _n1, _n2, count=0;

//The first run over the bonds: search for N=C-O-H
_i=str->nedges;
while(_i--) 
  if (str->edges[_i].type=='2')
    {//case R2-(N=(C-Y))-R1: amide if ( Y == OH or Y == SH )
    if ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_NITROGEN)&&(str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_CARBON) )
      { _t=str->edges[_i].vertice[0], str->edges[_i].vertice[0]=str->edges[_i].vertice[1], str->edges[_i].vertice[1]=_t; goto ANALYZE_CN; }
    if ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_CARBON)&&(str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_NITROGEN) )
      {   
      ANALYZE_CN: ;
      if ( (!(order[str->edges[_i].vertice[0]][3]))&&(order[str->edges[_i].vertice[0]][2]==1)&&(order[str->edges[_i].vertice[0]][1]+order[str->edges[_i].vertice[0]][0]==2) )
        {
        //Check for C-O 
        _n1=2, _j=str->nedges;
        while (_j--) 
          {
          if ( (str->edges[_j].vertice[1]==str->edges[_i].vertice[0])&&(_j!=_i) )
            { _t=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=str->edges[_j].vertice[1], str->edges[_j].vertice[1]=_t; goto ANALYZE_CO; }
          if ( (str->edges[_j].vertice[0]==str->edges[_i].vertice[0])&&(_j!=_i) )
            {
            ANALYZE_CO: ; 
            if ( (!(order[str->edges[_j].vertice[1]][3]))&&(!(order[str->edges[_j].vertice[1]][2]))&&(order[str->edges[_j].vertice[1]][1]+order[str->edges[_j].vertice[1]][0]==2)&&( (str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
              {
              //Check for O-H
              _n2=1, _k=str->nedges;
              while (_k--) 
                {
                if ( (str->edges[_k].vertice[1]==str->edges[_j].vertice[1])&&(_k!=_j)&&(_k!=_i) )
                  { _t=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=str->edges[_j].vertice[1], str->edges[_j].vertice[1]=_t; goto ANALYZE_OH; }
                if ( (str->edges[_k].vertice[0]==str->edges[_j].vertice[1])&&(_k!=_j)&&(_k!=_i) )
                  {
                  ANALYZE_OH: ;
                  if ( (!(order[str->edges[_k].vertice[1]][3]))&&(!(order[str->edges[_k].vertice[1]][2]))&&(order[str->edges[_k].vertice[1]][1]+order[str->edges[_k].vertice[1]][0]==1)&&(str->a[str->edges[_k].vertice[1]]==CHEM_ATOM_TYPE_HYDROGEN) )
                    {
                    str->r[str->edges[_k].vertice[1]].i=str->r[str->edges[_k].vertice[1]].j=str->r[str->edges[_k].vertice[1]].k=(float)NAN;
                    if (str->edges[_k].type!='1')
                      {
                      order[str->edges[_k].vertice[0]][1]--; 
                      order[str->edges[_k].vertice[1]][1]--, order[str->edges[_k].vertice[1]][0]++; 
                      str->edges[_k].type='1'; 
                      }
                    else order[str->edges[_k].vertice[0]][0]--;
                    order[str->edges[_k].vertice[0]][2]++;
                    order[str->edges[_i].vertice[1]][2]--, order[str->edges[_i].vertice[1]][1]++; 
                    str->edges[_k].vertice[0]=str->edges[_i].vertice[1]; 
                    str->edges[_j].type='2';
                    if (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_CARBON) str->edges[_i].type=(int)'m';
                    goto NEXT_I;
                    }
                  if (!(--_n2)) break;
                  }
                }
              }
            if (!(--_n1)) break;
            }      
          }
        }
      }
    NEXT_I: ;
    }
//The second run over the bonds
_i=str->nedges;
while(_i--) 
  if (str->edges[_i].type=='2')
    {//case R1-((C=Y)-NH)-R2: amide if ( Y == O or Y == S )
    if ( (str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_CARBON)&&
       ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_SULFUR) ) )
      { _t=str->edges[_i].vertice[0], str->edges[_i].vertice[0]=str->edges[_i].vertice[1], str->edges[_i].vertice[1]=_t; goto ANALYZE__CO; }
    if ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_CARBON)&&
       ( (str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
      {   
      ANALYZE__CO: ;
      if ( (!(order[str->edges[_i].vertice[0]][3]))&&(order[str->edges[_i].vertice[0]][2]==1)&&(order[str->edges[_i].vertice[0]][1]+order[str->edges[_i].vertice[0]][0]==2)&&
           (!(order[str->edges[_i].vertice[1]][3]))&&(order[str->edges[_i].vertice[1]][2]==1)&&(!(order[str->edges[_i].vertice[1]][1]))&&(!(order[str->edges[_i].vertice[1]][0])) )
        {//Check for C-N
        _n1=2;
        _j=str->nedges;
        while (_j--) 
          {
          if ( (str->edges[_j].vertice[1]==str->edges[_i].vertice[0])&&(_i!=_j) )
            { _t=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=str->edges[_j].vertice[1], str->edges[_j].vertice[1]=_t; goto ANALYZE__CN; }
          if ( (str->edges[_j].vertice[0]==str->edges[_i].vertice[0])&&(_i!=_j) )
            {
            ANALYZE__CN: ;
            if ( (!(order[str->edges[_j].vertice[1]][3]))&&(!(order[str->edges[_j].vertice[1]][2]))&&(order[str->edges[_j].vertice[1]][1]+order[str->edges[_j].vertice[1]][0]==3)&&(str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_NITROGEN) )
              {
              //Check for N-H
              _n2=2, _k=str->nedges;
              while (_k--) 
                {
                if ( (str->edges[_k].vertice[1]==str->edges[_j].vertice[1])&&(_k!=_j) )
                  { _t=str->edges[_k].vertice[0], str->edges[_k].vertice[0]=str->edges[_k].vertice[1], str->edges[_k].vertice[1]=_t; goto ANALYZE__NH; }
                if ( (str->edges[_k].vertice[0]==str->edges[_j].vertice[1])&&(_k!=_j) )
                  {
                  ANALYZE__NH: ;
                  if ( (!(order[str->edges[_k].vertice[1]][3]))&&(!(order[str->edges[_k].vertice[1]][2]))&&(order[str->edges[_k].vertice[1]][1]+order[str->edges[_k].vertice[1]][0]==1)&&(str->a[str->edges[_k].vertice[1]]==CHEM_ATOM_TYPE_HYDROGEN) )
                    {
                    if (str->edges[_k].type!='1')
                      {
                      order[str->edges[_k].vertice[0]][1]--, order[str->edges[_k].vertice[0]][0]++; 
                      order[str->edges[_k].vertice[1]][1]--, order[str->edges[_k].vertice[1]][0]++; 
                      str->edges[_k].type='1';
                      }
                    if (str->edges[_j].type=='1') 
                      {
                      order[str->edges[_j].vertice[0]][1]++, order[str->edges[_j].vertice[0]][0]--; 
                      order[str->edges[_j].vertice[1]][1]++, order[str->edges[_j].vertice[1]][0]--; 
                      }
                    str->edges[_j].type=(int)'m';
                    count++;
                    break;
                    }
                  if (!(--_n2)) break;  
                  }
                }
              }
            if (!(--_n1)) break;
            }    
          }
        }
      }
    }
return count;
}

//Rule 0. H(?)N=(C-NH2)-NH-R -> H2N=(+C~NH2)~NH-R 
//Rule 1. H(?)N=(C~NH2)-R -> H2N=(+C~NH2)-R
//Rule 2. H(?)N=(C-R)-NH-R2 -> H2N-(C-R)=N-R2
char _update_guanidines_rule012(unsigned int _k,unsigned int cid, unsigned int nid,t_clist *_neighbors,char (**order)[4],t_str *str)
{
register unsigned int _i, _j, _l;
void *vp;
     if (_neighbors->list[cid].list[0]==nid) { _i=_neighbors->list[cid].list[1], _j=_neighbors->list[cid].list[2]; }
else if (_neighbors->list[cid].list[1]==nid) { _i=_neighbors->list[cid].list[0], _j=_neighbors->list[cid].list[2]; }
else                                         { _i=_neighbors->list[cid].list[0], _j=_neighbors->list[cid].list[1]; } 
if ( (str->a[_i]!=CHEM_ATOM_TYPE_NITROGEN)||(_neighbors->list[_i].size!=1)||( ((*order)[_i][3]))||( ((*order)[_i][2]))||(_neighbors->list[_i].size==(*order)[_i][1]+(*order)[_i][0]) )
  {
  if ( (str->a[_j]!=CHEM_ATOM_TYPE_NITROGEN)||(_neighbors->list[_j].size!=1)||( ((*order)[_j][3]))||( ((*order)[_j][2]))||(_neighbors->list[_i].size==(*order)[_j][1]+(*order)[_j][0]) ) 
    {//There is no C-NH2
    //Rule 2. H(?)N=(C-R)-NH-R2 -> H2N-(C-R)=N-R2
    if ( (str->a[_i]!=CHEM_ATOM_TYPE_NITROGEN)||(_neighbors->list[_i].size!=2)||( ((*order)[_i][3]))||( ((*order)[_i][2]))||((*order)[_i][1]+(*order)[_i][0]!=3) )
      {
      if ( (str->a[_j]!=CHEM_ATOM_TYPE_NITROGEN)||(_neighbors->list[_j].size!=2)||( ((*order)[_j][3]))||( ((*order)[_j][2]))||((*order)[_j][1]+(*order)[_j][0]!=3) ) return NTNF;
      else _i=_j;
      }
    //Edit  _i and nid
    _l=str->nedges;
    while (_l--)
           if ( (str->a[str->edges[_l].vertice[0]]==CHEM_ATOM_TYPE_HYDROGEN)&&(str->edges[_l].vertice[1]==_i) ) { str->edges[_l].vertice[1]=nid; break; }
      else if ( (str->edges[_l].vertice[0]==_i)&&(str->a[str->edges[_l].vertice[1]]==CHEM_ATOM_TYPE_HYDROGEN) ) { str->edges[_l].vertice[0]=nid; break; }
    _l=str->nedges;
    while (_l--)
      if ( ( (str->edges[_l].vertice[0]==cid)&&(str->edges[_l].vertice[1]==_i) )||( (str->edges[_l].vertice[0]==_i)&&(str->edges[_l].vertice[1]==cid) ) )
        { 
        if ( (str->edges[_l].type==(int)'m')||(str->edges[_l].type==(int)'a') )
          { (*order)[cid][1]--, (*order)[cid][0]++, (*order)[_i][1]--, (*order)[_i][0]++, str->edges[_l].type='1'; }
        str->edges[_l].vertice[0]=cid, str->edges[_l].vertice[1]=nid;
        str->edges[_k].vertice[0]=cid, str->edges[_k].vertice[1]=_i;
        (*order)[_i][0]--, (*order)[_i][2]++, (*order)[nid][2]--, (*order)[nid][0]++;
        break;
        }
    return NTNF;
    }
  else { _l=_i, _i=_j, _j=_l; } //C-NH2 is detected: Rule 1. N=(C~NH2)-R -> H2N=(+C~NH2)-R || Rule 2. N=(C-R)-NH-R2 -> H2N-(C-R)=N-R2
  }    
if ( (str->a[_j]==CHEM_ATOM_TYPE_NITROGEN)&&(!((*order)[_j][3]))&&(!((*order)[_j][2]))&&(_neighbors->list[_i].size!=(*order)[_j][1]+(*order)[_j][0]) )
  {//Rule 0. HN=(C-NH2)-NH-R -> H2N=(+C~NH2)~NH-R
  _l=str->nedges;
  while (_l--) //Update cid-_j edge for double case
    if ( ( (str->edges[_l].vertice[0]==cid)&&(str->edges[_l].vertice[1]==_j) )||( (str->edges[_l].vertice[0]==_j)&&(str->edges[_l].vertice[1]==cid) ) ) 
      {
      if (str->edges[_l].type=='1')
        { (*order)[cid][0]--, (*order)[cid][1]++, (*order)[_j][0]--, (*order)[_j][1]++; }
      str->edges[_l].type='a';
      break;
      }
  }
_l=str->nedges;
while (_l--) //Update cid-_i edge in any case
  if ( ( (str->edges[_l].vertice[0]==cid)&&(str->edges[_l].vertice[1]==_i) )||( (str->edges[_l].vertice[0]==_i)&&(str->edges[_l].vertice[1]==cid) ) )
    {
    if (str->edges[_l].type=='1')
      { (*order)[cid][0]--, (*order)[cid][1]++, (*order)[_i][0]--, (*order)[_i][1]++; }
    str->edges[_l].type='a';
    break;
    }
//Edit HN=C -> H2N=C if needed
if ((*order)[nid][0]+(*order)[nid][1]==1)
  {
  (*order)[nid][0]++;
  if ((_k=str_add_hydrogen(nid,str))==(unsigned int)-1) return FALSE;
  if (!(vp=(void*)realloc((*order),str->natoms*sizeof(unsigned int)))) { ylib_errno=YERROR_MEMORY; return FALSE; } 
  else 
    {//Add (+) proton
    (*order)=(char(*)[4])vp;
    _l=str->natoms; while (--_l!=_k) *((unsigned int*)&(*order)[_l])=*((unsigned int*)&(*order)[_l-1]);
    (*order)[_l][3]=(*order)[_l][2]=(*order)[_l][1]=0, (*order)[_l][0]=1;
    }
  }
return TRUE;
}
//Rule 3. Update if it is R-N=C(NH2)2 -> R-NH~(+C~NH2)=NH2
char _update_guanidines_rule3(unsigned int _k,unsigned int cid, unsigned int nid,t_clist *_neighbors,char (**order)[4],t_str *str)
{
register unsigned int _i, _j, _l;
void *vp;
     if (_neighbors->list[cid].list[0]==nid) { _i=_neighbors->list[cid].list[1], _j=_neighbors->list[cid].list[2]; }
else if (_neighbors->list[cid].list[1]==nid) { _i=_neighbors->list[cid].list[0], _j=_neighbors->list[cid].list[2]; }
else                                                { _i=_neighbors->list[cid].list[0], _j=_neighbors->list[cid].list[1]; } 
if ( ( (str->a[_i]!=CHEM_ATOM_TYPE_NITROGEN)||(_neighbors->list[_i].size!=1)||( ((*order)[_i][3]))||( ((*order)[_i][2]))||((*order)[_i][1]+(*order)[_i][0]!=3) )||
     ( (str->a[_j]!=CHEM_ATOM_TYPE_NITROGEN)||(_neighbors->list[_j].size!=1)||( ((*order)[_j][3]))||( ((*order)[_j][2]))||((*order)[_j][1]+(*order)[_j][0]!=3) ) )
  return NTNF;
//Edit _i and nid 
_l=str->nedges;
while (_l--)
  if ( ( (str->edges[_l].vertice[0]==cid)&&(str->edges[_l].vertice[1]==_i) )||( (str->edges[_l].vertice[0]==_i)&&(str->edges[_l].vertice[1]==cid) ) )
    {
    if (str->edges[_l].type=='1')
      { (*order)[cid][0]--, (*order)[cid][1]++, (*order)[_i][0]--, (*order)[_i][1]++; }
    str->edges[_l].vertice[0]=cid, str->edges[_l].vertice[1]=nid, str->edges[_l].type='a';
    str->edges[_k].vertice[0]=cid, str->edges[_k].vertice[1]=_i;
    (*order)[_i][1]--, (*order)[_i][2]++, (*order)[nid][2]--, (*order)[nid][1]++;
    break;
    }
//Edit _j
_l=str->nedges;
while (_l--)
  if ( ( (str->edges[_l].vertice[0]==cid)&&(str->edges[_l].vertice[1]==_j) )||( (str->edges[_l].vertice[0]==_j)&&(str->edges[_l].vertice[1]==cid) ) )
    {
    if (str->edges[_l].type=='1')
      {  (*order)[cid][0]--, (*order)[cid][1]++, (*order)[_j][0]--, (*order)[_j][1]++; }
    str->edges[_l].type='a';
    break;
    }
//Add proton on _i
(*order)[_i][0]++;
if ((_k=str_add_hydrogen(nid,str))==(unsigned int)-1) return FALSE;
if (!(vp=(void*)realloc((*order),str->natoms*sizeof(unsigned int)))) { ylib_errno=YERROR_MEMORY; return FALSE; } 
else 
  {//Add (+) proton
  (*order)=(char(*)[4])vp;
  _l=str->natoms; while (--_l!=_k) *((unsigned int*)&(*order)[_l])=*((unsigned int*)&(*order)[_l-1]);
  (*order)[_l][3]=(*order)[_l][2]=(*order)[_l][1]=0, (*order)[_l][0]=1;
  }
return TRUE;
}

//This function edits R-N=C-(NH2)2 -> R-NH-(C=NH)-NH2
//Note. It should be called after amide bond definition rutine
//Rule 0. H(?)N=(C-NH2)-NH-R -> H2N=(+C~NH2)~NH-R 
//Rule 1. H(?)N=(C~NH2)-R -> H2N=(+C~NH2)-R
//Rule 2. H(?)N=(C-R)-NH-R2 -> H2N-(C-R)=N-R2
//Rule 3. Update if it is R-N=C(NH2)2 -> R-NH~(+C~NH2)=NH2
unsigned int update_guanidines(char (**order)[4],t_clist *_neighbors,t_str *str)
{
register unsigned int _i, _j, _k, _t, count=0;

_k=str->nedges;
while (_k--)
  if (str->edges[_k].type=='2')
    {
    _i=str->edges[_k].vertice[0], _j=str->edges[_k].vertice[1];
    if ( (str->a[_i]==CHEM_ATOM_TYPE_NITROGEN)&&(str->a[_j]==CHEM_ATOM_TYPE_CARBON) )
      { _t=_i, _i=_j, _j=_t; goto ANALYZE_CN; }
    if ( (str->a[_i]==CHEM_ATOM_TYPE_CARBON)&&(str->a[_j]==CHEM_ATOM_TYPE_NITROGEN) )
      {
      ANALYZE_CN: ;
      if ( ( ( ((*order)[_i][3]))||( ((*order)[_i][2]!=1))||(_neighbors->list[_i].size!=3)||((*order)[_i][0]+(*order)[_i][1]!=2) )||
           ( ( ((*order)[_j][3]))||( ((*order)[_j][2]!=1))||(_neighbors->list[_j].size>2)||( ((*order)[_j][0]+(*order)[_j][1]!=1)&&((*order)[_j][0]+(*order)[_j][1]!=2) ) ) )
        continue; //Confirm it's [1-3]N=C[3]
      if (_neighbors->list[_j].size==1)
        {//It's H[1-2]N=C[3]
        if (!(_t=(unsigned int)_update_guanidines_rule012(_k,_i,_j,_neighbors,order,str))) return (unsigned int)-1;
        else if (_t==TRUE) count++; //Rule 0,1,2.
        }
      else 
        {//It's R-N=C[3]
        if (!(_t=(unsigned int)_update_guanidines_rule3  (_k,_i,_j,_neighbors,order,str))) return (unsigned int)-1;
        else if (_t==TRUE) count++; //Rule 3.
        }
      }
    }
return count;
}

//This function edits R(0-3)-N -> R(0-3)-NH(+)
unsigned int update_amines(char (**order)[4],t_clist *_neighbors,t_str *str)
{
register unsigned int _i, _j, _k, count=0;
void *vp; 

_i=str->natoms;
while (_i--)
  if ( (str->a[_i]==CHEM_ATOM_TYPE_NITROGEN)&&(!((*order)[_i][3]))&&(!((*order)[_i][2]))&&(!((*order)[_i][1])) )
    {
    _j=_neighbors->list[_i].size;
    while (_j--)
      if ( (str->a[_k=_neighbors->list[_i].list[_j]]!=CHEM_ATOM_TYPE_CARBON)||( ((*order)[_k][3]))||( ((*order)[_k][2])) )
        goto NEXT_I;
         if ((*order)[_i][0]==3)
           {
           (*order)[_i][0]++;
           if ((_j=str_add_hydrogen(_i,str))==(unsigned int)-1) return (unsigned int)-1;
           if (!(vp=(void*)realloc((*order),str->natoms*sizeof(unsigned int)))) { ylib_errno=YERROR_MEMORY; return (unsigned int)-1; } 
           else 
             {//Add (+) proton
             (*order)=(char(*)[4])vp;
             _k=str->natoms;
             while (--_k!=_j) *((unsigned int*)&(*order)[_k])=*((unsigned int*)&(*order)[_k-1]);
             (*order)[_k][3]=(*order)[_k][2]=(*order)[_k][1]=0, (*order)[_k][0]=1;
             }
           count++;
           }
    else if ((*order)[_i][0]==4) count++;
    NEXT_I: ;
    }
return count;
}

//This function edits R-CO-OH -> R-CO~O(-) and R-CO-O(-) -> R-CO~O(-) 
unsigned int update_carboxyles(char (*order)[4],t_str *str)
{
register unsigned int _i, _j, _k, _t, _n1, count=0;

_i=str->nedges;
while (_i--)
  if (str->edges[_i].type=='2')
    {
    if ( ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_SULFUR) )&&(str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_CARBON) )
      { _t=str->edges[_i].vertice[0], str->edges[_i].vertice[0]=str->edges[_i].vertice[1], str->edges[_i].vertice[1]=_t; goto ANALYZE_CO; }
    if ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_CARBON)&&( (str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
      {
      ANALYZE_CO: ;
      if ( (!(order[str->edges[_i].vertice[0]][3]))&&(order[str->edges[_i].vertice[0]][2]==1)&&(order[str->edges[_i].vertice[0]][1]+order[str->edges[_i].vertice[0]][0]==2)&&
           (!(order[str->edges[_i].vertice[1]][3]))&&(order[str->edges[_i].vertice[1]][2]==1)&&(!(order[str->edges[_i].vertice[1]][1]))&&(!(order[str->edges[_i].vertice[1]][0])) )
        {
        _n1=2, _j=str->nedges;
        while (_j--)
          {
          if ( (str->edges[_j].vertice[1]==str->edges[_i].vertice[0])&&(_j!=_i) )
            { _t=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=str->edges[_j].vertice[1], str->edges[_j].vertice[1]=_t; goto ANALYZE_O; }
          if ( (str->edges[_j].vertice[0]==str->edges[_i].vertice[0])&&(_j!=_i) )
            {
            ANALYZE_O: ;
            if ( (!(order[str->edges[_j].vertice[1]][3]))&&(!(order[str->edges[_j].vertice[1]][2]))&&( (order[str->edges[_j].vertice[1]][1]+order[str->edges[_j].vertice[1]][0]==2)||(order[str->edges[_j].vertice[1]][1]+order[str->edges[_j].vertice[1]][0]==1) )&&( (str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
              {
              if (order[str->edges[_j].vertice[1]][1]+order[str->edges[_j].vertice[1]][0]==2)
                {
                _k=str->nedges;
                while (_k--)
                  {
                  if ( (str->edges[_k].vertice[1]==str->edges[_j].vertice[1])&&(_k!=_j) )
                    { _t=str->edges[_k].vertice[0], str->edges[_k].vertice[0]=str->edges[_k].vertice[1], str->edges[_k].vertice[1]=_t; goto ANALYZE_OH; }
                  if ( (str->edges[_k].vertice[0]==str->edges[_j].vertice[1])&&(_k!=_j) )
                    {
                    ANALYZE_OH: ;
                    if ( (!(order[str->edges[_k].vertice[1]][3]))&&(!(order[str->edges[_k].vertice[1]][2]))&&(order[str->edges[_k].vertice[1]][1]+order[str->edges[_k].vertice[1]][0]==1)&&(str->a[str->edges[_k].vertice[1]]==CHEM_ATOM_TYPE_HYDROGEN) )
                      {//Edit O=C-O-H -> O=C~O(-)
                      if (str->edges[_j].type=='1') 
                        { 
                        order[str->edges[_j].vertice[0]][0]--, order[str->edges[_j].vertice[0]][1]++;
                        order[str->edges[_j].vertice[1]][0]--, order[str->edges[_j].vertice[1]][1]++;
                        str->edges[_j].type='a';
                        }
                      if (str->edges[_k].type=='1') { order[str->edges[_k].vertice[0]][0]--; } else { order[str->edges[_k].vertice[0]][1]--; }
                      if ((_k=str_del_hydrogen(str->edges[_k].vertice[1],str))==(unsigned int)-1) return (unsigned int)-1;
                      else { while (_k!=str->natoms) { *((unsigned int*)&order[_k])=*((unsigned int*)&order[_k+1]), _k++; } }
                      count++;
                      //We could ocassionaly move the bonds when deleted the hydrogen so restart the search to be on a safe side
                      if (_i<str->nedges) _i++; else _i=str->nedges;
                      goto NEXT_I;
                      }
                    break; 
                    }
                  }
                }
              else if (str->edges[_j].type=='1') str->edges[_j].type='a';
              }
            if (!(--_n1)) break; 
            }
          }
        }
      }
    NEXT_I: ;
    }
return count;
}

//This function edits R-(Si=O)-OH -> R->(Si=O)-O(-)
unsigned int update_silicates(char (*order)[4],t_str *str)
{
register unsigned int _i, _j, _k, _t, _n1, count=0;

_i=str->nedges;
while (_i--)
  if (str->edges[_i].type=='2')
    {
    if ( ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_SULFUR) )&&(str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_SILICON) )
      { _t=str->edges[_i].vertice[0], str->edges[_i].vertice[0]=str->edges[_i].vertice[1], str->edges[_i].vertice[1]=_t; goto ANALYZE_SiO; }
    if ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_SILICON)&&( (str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
      {
      ANALYZE_SiO: ;
      if ( (!(order[str->edges[_i].vertice[0]][3]))&&(order[str->edges[_i].vertice[0]][2]==1)&&(order[str->edges[_i].vertice[0]][1]+order[str->edges[_i].vertice[0]][0]==2)&&
           (!(order[str->edges[_i].vertice[1]][3]))&&(order[str->edges[_i].vertice[1]][2]==1)&&(!(order[str->edges[_i].vertice[1]][1]))&&(!(order[str->edges[_i].vertice[1]][0])) )
        {
        _n1=2, _j=str->nedges;
        while (_j--)
          {
          if ( (str->edges[_j].vertice[1]==str->edges[_i].vertice[0])&&(_j!=_i) )
            { _t=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=str->edges[_j].vertice[1], str->edges[_j].vertice[1]=_t; goto ANALYZE_O; }
          if ( (str->edges[_j].vertice[0]==str->edges[_i].vertice[0])&&(_j!=_i) )
            {
            ANALYZE_O: ;
            if ( (!(order[str->edges[_j].vertice[1]][3]))&&(!(order[str->edges[_j].vertice[1]][2]))&&(order[str->edges[_j].vertice[1]][1]+order[str->edges[_j].vertice[1]][0]==2)&&( (str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
              {
              _k=str->nedges;
              while (_k--)
                {
                if ( (str->edges[_k].vertice[1]==str->edges[_j].vertice[1])&&(_k!=_j) )
                  { _t=str->edges[_k].vertice[0], str->edges[_k].vertice[0]=str->edges[_k].vertice[1], str->edges[_k].vertice[1]=_t; goto ANALYZE_OH; }
                if ( (str->edges[_k].vertice[0]==str->edges[_j].vertice[1])&&(_k!=_j) )
                  {
                  ANALYZE_OH: ;
                  if ( (!(order[str->edges[_k].vertice[1]][3]))&&(!(order[str->edges[_k].vertice[1]][2]))&&(order[str->edges[_k].vertice[0]][1]+order[str->edges[_k].vertice[0]][0]==1)&&(str->a[str->edges[_k].vertice[1]]==CHEM_ATOM_TYPE_HYDROGEN) )
                    {//Edit O=C-O-H -> O=CaO(-)
                    if (str->edges[_j].type=='1') 
                      { 
                      order[str->edges[_j].vertice[0]][0]--, order[str->edges[_j].vertice[0]][1]++;
                      order[str->edges[_j].vertice[1]][0]--, order[str->edges[_j].vertice[1]][1]++;
                      str->edges[_j].type='a';
                      }
                    if (str->edges[_k].type=='1') { order[str->edges[_k].vertice[0]][0]--; } else { order[str->edges[_k].vertice[0]][1]--; }
                    if ((_k=str_del_hydrogen(str->edges[_k].vertice[1],str))==(unsigned int)-1) return (unsigned int)-1;
                    else { while (_k!=str->natoms) { *((unsigned int*)&order[_k])=*((unsigned int*)&order[_k+1]), _k++; } }
                    count++;
                    }
                  break; 
                  }
                }
              }
            if (!(--_n1)) break; 
            }
          }
        }
      }
    }
return count;
}

//This function edits R-(P=O)-(OH)2 -> R->(P=O)-O2(2-)
unsigned int update_phosphates(char (*order)[4],t_str *str)
{
register unsigned int _i, _j, _k, _t, _n1, count=0;

_i=str->nedges;
while (_i--)
  if (str->edges[_i].type=='2')
    {
    if ( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_PHOSPHOR)&&( (str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
      { _t=str->edges[_i].vertice[0], str->edges[_i].vertice[0]=str->edges[_i].vertice[1], str->edges[_i].vertice[1]=_t; goto ANALYZE_PO; }
    if ( (str->a[str->edges[_i].vertice[1]]==CHEM_ATOM_TYPE_PHOSPHOR)&&( (str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_i].vertice[0]]==CHEM_ATOM_TYPE_SULFUR) ) )
      {//P to attached P=O 
      ANALYZE_PO:
      if ( (!(order[str->edges[_i].vertice[0]][3]))&&(2*order[str->edges[_i].vertice[0]][2]+order[str->edges[_i].vertice[0]][1]+order[str->edges[_i].vertice[0]][0]==4)&&
           (!(order[str->edges[_i].vertice[1]][3]))&&(order[str->edges[_i].vertice[1]][2]==1)&&(!(order[str->edges[_i].vertice[1]][1]))&&(!(order[str->edges[_i].vertice[1]][0])) )
        {
        _n1=3, _j=str->nedges;
        while (_j--)
          {
          if ( (str->edges[_j].vertice[1]==str->edges[_i].vertice[0])&&(_i!=_j) )
            { _t=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=str->edges[_j].vertice[1], str->edges[_j].vertice[1]=_t; goto ANALYZE_P5; }
          if ( (str->edges[_j].vertice[0]==str->edges[_i].vertice[0])&&(_i!=_j) )
            {//O=P(+5) look for attached OH
            ANALYZE_P5: ;
            if ( (!(order[str->edges[_j].vertice[1]][3]))&&(!(order[str->edges[_j].vertice[1]][2]))&&(order[str->edges[_j].vertice[1]][1]+order[str->edges[_j].vertice[1]][0]==2)&&( (str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
              {
              _k=str->nedges;
              while (_k--)
                {
                if ( (str->edges[_k].vertice[1]==str->edges[_j].vertice[1])&&(_j!=_k) )
                  { _t=str->edges[_k].vertice[0], str->edges[_k].vertice[0]=str->edges[_k].vertice[1], str->edges[_k].vertice[1]=_t; goto ANALYZE_OH; }
                if ( (str->edges[_k].vertice[0]==str->edges[_j].vertice[1])&&(_j!=_k) ) 
                  {
                  ANALYZE_OH: ;
                  if ( (!(order[str->edges[_k].vertice[1]][3]))&&(!(order[str->edges[_k].vertice[1]][2]))&&(order[str->edges[_k].vertice[1]][1]+order[str->edges[_k].vertice[1]][0]==1)&&(str->a[str->edges[_k].vertice[1]]==CHEM_ATOM_TYPE_HYDROGEN) )
                    {//O=P[+5]-O-H -> O=P[+5]-O(-)
                    if (str->edges[_k].type=='1') order[str->edges[_k].vertice[0]][0]--; else order[str->edges[_k].vertice[0]][1]--;
                    if (str->edges[_j].type=='1') 
                      {
                      order[str->edges[_j].vertice[0]][0]--, order[str->edges[_j].vertice[0]][1]++;  
                      order[str->edges[_j].vertice[1]][0]--, order[str->edges[_j].vertice[1]][1]++;
                      str->edges[_j].type='a';  
                      }
                    if ((_k=str_del_hydrogen(str->edges[_k].vertice[1],str))==(unsigned int)-1) return (unsigned int)-1;
                    else { while (_k!=str->natoms) { *((unsigned int*)&order[_k])=*((unsigned int*)&order[_k+1]), _k++; } }
                    count++;
                    }
                  break;
                  }
                }
              }
            }
          if (!(--_n1)) break;
          }
        }
      }
    }
return count;
} 

//This function edits R-SO2-OH -> R-SO2-O(-) and R-SO-OH -> R-SO-O(-)
unsigned int update_sulfaites(char (*order)[4],t_str *str)
{
register unsigned int _i, _j, _k, _t, _n1, count=0;

_i=str->natoms;
while (_i--)
  if ( (str->a[_i]==CHEM_ATOM_TYPE_SULFUR)&&(!(order[_i][3]))&&( (order[_i][2]))&&
     ( ( (2*order[_i][2]+order[_i][1]+order[_i][0]==6)&&(order[_i][2]+order[_i][1]+order[_i][0]>=4) )||
       ( (2*order[_i][2]+order[_i][1]+order[_i][0]==4)&&(order[_i][2]+order[_i][1]+order[_i][0]==3) ) ) ) 
    {//S(+4/+6) look for attached S=O 
    _n1=order[_i][3]+order[_i][2]+order[_i][1]+order[_i][0], _j=str->nedges;
    while (_j--)
      {
      if (str->edges[_j].vertice[1]==_i)
        { str->edges[_j].vertice[1]=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=_i; goto ANALYZE_S; }
      if (str->edges[_j].vertice[0]==_i)
        {
        ANALYZE_S: ;
        if ( (str->edges[_j].type=='2')&&(!(order[str->edges[_j].vertice[1]][3]))&&(order[str->edges[_j].vertice[1]][2]==1)&&(!(order[str->edges[_j].vertice[1]][1]))&&(!(order[str->edges[_j].vertice[1]][0]))&&( (str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
          {
          _n1=order[_i][3]+order[_i][2]+order[_i][1]+order[_i][0], _j=str->nedges;
          while (_j--)
            {
            if (str->edges[_j].vertice[1]==_i)
              { str->edges[_j].vertice[1]=str->edges[_j].vertice[0], str->edges[_j].vertice[0]=_i; goto ANALYZE_SO; }
            if (str->edges[_j].vertice[0]==_i)
              {
              ANALYZE_SO: ;
              if ( (!(order[str->edges[_j].vertice[1]][3]))&&(!(order[str->edges[_j].vertice[1]][2]))&&(order[str->edges[_j].vertice[1]][1]+order[str->edges[_j].vertice[1]][0]==2)&&( (str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[str->edges[_j].vertice[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
                {
                _k=str->nedges;
                while (_k--)
                  {
                  if ( (str->edges[_k].vertice[1]==str->edges[_j].vertice[1])&&(_j!=_k) )
                    { _t=str->edges[_k].vertice[0], str->edges[_k].vertice[0]=str->edges[_k].vertice[1], str->edges[_k].vertice[1]=_t; goto ANALYZE_OH; }
                  if ( (str->edges[_k].vertice[0]==str->edges[_j].vertice[1])&&(_j!=_k) ) 
                    {
                    ANALYZE_OH: ;
                    if ( (!(order[str->edges[_k].vertice[1]][3]))&&(!(order[str->edges[_k].vertice[1]][2]))&&(order[str->edges[_k].vertice[1]][1]+order[str->edges[_k].vertice[1]][0]==1)&&(str->a[str->edges[_k].vertice[1]]==CHEM_ATOM_TYPE_HYDROGEN) )
                      {// S[+6]-O-H -> S[+6]-O(-)
                      if (str->edges[_k].type=='1') order[str->edges[_k].vertice[0]][0]--; else order[str->edges[_k].vertice[0]][1]--;
                      if (str->edges[_j].type=='1') 
                        {
                        order[str->edges[_j].vertice[0]][0]--, order[str->edges[_j].vertice[0]][1]++;  
                        order[str->edges[_j].vertice[1]][0]--, order[str->edges[_j].vertice[1]][1]++;
                        str->edges[_j].type='a';  
                        }
                      if ((_k=str_del_hydrogen(str->edges[_k].vertice[1],str))==(unsigned int)-1) return (unsigned int)-1;
                      else { while (_k!=str->natoms) { *((unsigned int*)&order[_k])=*((unsigned int*)&order[_k+1]), _k++; } }
                      count++;
                      }
                    break;
                    }
                  }
                }
              if (!(--_n1)) break;
              }
            }
          }
        else _n1--;
        if (!(_n1)) break;
        }
      }
    }
return count;
} 


//Since here we activating neighbors because the algorithms are more complex

//The order[_i] consists of 4 chars -> | amount_of_3_bonds | amount_of_2_bonds | amount of amide/quasy conjugated bonds | amount_of_1_bonds | and OPTIONALLY the neighbours in list are sorted correspondingly

//This funcion calculates chemical 'order' of atoms
char *calc_order(char (**order)[4],unsigned int natoms,unsigned int nedges,t_edge *edges)
{
char (*_order)[4];
if ( ( (order))&&( (*order)) ) { _order=(*order); memset(_order,0,sizeof(char)*4*natoms); }
else if (!(_order=(char(*)[4])calloc(natoms,sizeof(char)*0x4))) { ylib_errno=YERROR_MEMORY; return FALSE; }
while (nedges--) 
  if ( (edges[nedges].vertice[0]<natoms)||(edges[nedges].vertice[1]<natoms) )
    switch (edges[nedges].type)
      {
      case '1' : { _order[edges[nedges].vertice[0]][0]++, _order[edges[nedges].vertice[1]][0]++; break; }
      case 'a' :
      case 'm' : { _order[edges[nedges].vertice[0]][1]++, _order[edges[nedges].vertice[1]][1]++; break; }
      case '2' : { _order[edges[nedges].vertice[0]][2]++, _order[edges[nedges].vertice[1]][2]++; break; }
      case '3' : { _order[edges[nedges].vertice[0]][3]++, _order[edges[nedges].vertice[1]][3]++; break; }
      default  : { LABEL_ERROR_DATA_CONSISTMENT: if ( (!(order))||(!(*order)) ) free(_order); ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
      }
  else goto LABEL_ERROR_DATA_CONSISTMENT;
if ( ( (order))&&(!(*order))) *order=(char (*)[4])_order;
return (char*)_order;
}

//This function create ordered neighbors list
void arrange_neighbors(char (*order)[4],t_clist *neighbors,unsigned int nvertices,unsigned int nedges,t_edge *edges)
{
register unsigned int _i, _j, _p, _id0, _id1;

//Stage 0. Init orders
_i=nvertices; while (_i--) { order[_i][0]=order[_i][1]=order[_i][2]=order[_i][3]=0; }
//Stage 1. Shift triple bonds
_i=nedges;
while (_i--)
  if (edges[_i].type=='3')
    {
    _id0=edges[_i].vertice[0], _id1=edges[_i].vertice[1];
    if (neighbors->list[_id0].list[_p=order[_id0][3]]!=_id1)
      {
      _j=neighbors->list[_id0].size; while (--_j!=_p) if (neighbors->list[_id0].list[_j]==_id1) { neighbors->list[_id0].list[_j]=neighbors->list[_id0].list[_p]; break; }
      neighbors->list[_id0].list[_p]=_id1;
      }
    if (neighbors->list[_id1].list[_p=order[_id1][3]]!=_id0)
      {
      _j=neighbors->list[_id1].size; while (--_j!=_p) if (neighbors->list[_id1].list[_j]==_id0) { neighbors->list[_id1].list[_j]=neighbors->list[_id1].list[_p]; break; }
      neighbors->list[_id1].list[_p]=_id0;
      }
    order[_id0][3]++, order[_id1][3]++;
    }
//Stage 2. Shift double bonds
_i=nedges;
while (_i--)
  if (edges[_i].type=='2')
    {
    _id0=edges[_i].vertice[0], _id1=edges[_i].vertice[1];
    if (neighbors->list[_id0].list[_p=order[_id0][3]+order[_id0][2]]!=_id1)
      {
      _j=neighbors->list[_id0].size; while (--_j!=_p) if (neighbors->list[_id0].list[_j]==_id1) { neighbors->list[_id0].list[_j]=neighbors->list[_id0].list[_p]; break; }
      neighbors->list[_id0].list[_p]=_id1;
      }
    if (neighbors->list[_id1].list[_p=order[_id1][3]+order[_id1][2]]!=_id0)
      {
      _j=neighbors->list[_id1].size; while (--_j!=_p) if (neighbors->list[_id1].list[_j]==_id0) { neighbors->list[_id1].list[_j]=neighbors->list[_id1].list[_p]; break; }
      neighbors->list[_id1].list[_p]=_id0;
      }
    order[_id0][2]++, order[_id1][2]++;
    }
//Stage 2. Shift resonance bonds
_i=nedges;
while (_i--)
  if ( (edges[_i].type==(int)'a')||(edges[_i].type==(int)'m') )
    {
    _id0=edges[_i].vertice[0], _id1=edges[_i].vertice[1];
    if (neighbors->list[_id0].list[_p=order[_id0][3]+order[_id0][2]+order[_id0][1]]!=_id1)
      {
      _j=neighbors->list[_id0].size; while (--_j!=_p) if (neighbors->list[_id0].list[_j]==_id1) { neighbors->list[_id0].list[_j]=neighbors->list[_id0].list[_p]; break; }
      neighbors->list[_id0].list[_p]=_id1;
      }
    if (neighbors->list[_id1].list[_p=order[_id1][3]+order[_id1][2]+order[_id1][1]]!=_id0)
      {
      _j=neighbors->list[_id1].size; while (--_j!=_p) if (neighbors->list[_id1].list[_j]==_id0) { neighbors->list[_id1].list[_j]=neighbors->list[_id1].list[_p]; break; }
      neighbors->list[_id1].list[_p]=_id0;
      }
    order[_id0][1]++, order[_id1][1]++;
    }
//Stage 3. Set single bonds amount
_i=nvertices; while (_i--) order[_i][0]=(char)neighbors->list[_i].size-order[_i][3]-order[_i][2]-order[_i][1];
}

//This function orders neighbours list accordingly to bonds' order 
char* define_neighbors_order(t_clist *neighbors,unsigned int nvertices,unsigned int nedges,t_edge *edges)
{
char (*order)[4];
register unsigned int _i, _id0, _id1;

if (!(order=(char (*)[4])malloc(sizeof(char)*0x4*nvertices))) { ylib_errno=YERROR_MEMORY; return FALSE; }
if ( (neighbors))
  {
  arrange_neighbors(order,neighbors,nvertices,nedges,edges);
  return (char*)order;
  }  
else
  {
  _i=nvertices; while (_i--) { order[_i][0]=order[_i][1]=order[_i][2]=order[_i][3]=0; }
  while (nedges--)
    {
    if ( ((_id0=edges[nedges].vertice[0])>nvertices)||((_id1=edges[nedges].vertice[1])>nvertices) )
      { LABEL_ERROR_IMPOSSIBLE: ylib_errno=YERROR_IMPOSSIBLE; free(order); order=0x0; return FALSE; }
    switch (edges[nedges].type)
      {
      case '1' : { order[_id0][0]++, order[_id1][0]++; break; }
      case 'a' :
      case 'm' : { order[_id0][1]++, order[_id1][1]++; break; }
      case '2' : { order[_id0][2]++, order[_id1][2]++; break; }
      case '3' : { order[_id0][3]++, order[_id1][3]++; break; }
      default  : { goto LABEL_ERROR_IMPOSSIBLE; }
      }
    }
  }
return (char*)order;
}

//
// This function need to be completed. It should resolve double bonds topology.
//
//This function correct topology of pi-electron subgraph
char update_double_bonds(char (*order)[4],t_clist *neighbors,t_str *str,t_top *top)
{
return FALSE;
}


//This function do compilation of initial str (apply all the above heuristics)
char compile_str(double pH,char verbose,char (**order)[4],t_clist **neighbors,t_str *str,t_top *top)
{
register unsigned int _i, _j;
t_clist *_neighbors;  //The list of heavy atom neighbors
void *vp;

//Init memory managemnent
if ( (str->natoms%0xFF))
  {
  _i=str->natoms+(0xFF-str->natoms%0xFF);
  if (!(vp=realloc(str->a,_i*sizeof(char)))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
  else str->a=(char*)vp;
  if (!(vp=realloc(str->anames,_i*sizeof(int)))) goto LABEL_MEMORY_ERROR;
  else str->anames=(char (*)[sizeof(int)])vp; 
  if (!(vp=realloc(str->r,_i*sizeof(t_vec)))) goto LABEL_MEMORY_ERROR;
  else str->r=(t_vec*)vp; 
  }
if ( (str->nedges%0xFF))
  {
  _i=str->nedges+(0xFF-str->nedges%0xFF);  
  if (!(vp=realloc(str->edges,_i*sizeof(t_edge)))) goto LABEL_MEMORY_ERROR;
  else str->edges=(t_edge*)vp; 
  }
//Compile heavy-atoms neighbours list (protons are about to rearrange)
_i=0, _j=str->nedges; while (_j--) if ( (str->a[str->edges[_j].vertice[0]]!=CHEM_ATOM_TYPE_HYDROGEN)&&(str->a[str->edges[_j].vertice[1]]!=CHEM_ATOM_TYPE_HYDROGEN) ) _i++;
if (!(_neighbors=alloc_clist(str->natoms,_i*2))) goto LABEL_MEMORY_ERROR;
_i=str->natoms; while (_i--) _neighbors->list[_i].size=0;
_j=str->nedges;
while (_j--)
  if ( (str->a[str->edges[_j].vertice[0]]!=CHEM_ATOM_TYPE_HYDROGEN)&&(str->a[str->edges[_j].vertice[1]]!=CHEM_ATOM_TYPE_HYDROGEN) )
    { _neighbors->list[str->edges[_j].vertice[0]].size++, _neighbors->list[str->edges[_j].vertice[1]].size++; }
for (_neighbors->list[0].list=_neighbors->_items, _i=1; _i<str->natoms; _i++) _neighbors->list[_i].list=_neighbors->list[_i-1].list+_neighbors->list[_i-1].size;
_i=str->natoms; while (_i--) _neighbors->list[_i].size=0;
_j=str->nedges;
while (_j--)
  if ( (str->a[str->edges[_j].vertice[0]]!=CHEM_ATOM_TYPE_HYDROGEN)&&(str->a[str->edges[_j].vertice[1]]!=CHEM_ATOM_TYPE_HYDROGEN) )
    {
    _neighbors->list[str->edges[_j].vertice[0]].list[_neighbors->list[str->edges[_j].vertice[0]].size++]=str->edges[_j].vertice[1],
    _neighbors->list[str->edges[_j].vertice[1]].list[_neighbors->list[str->edges[_j].vertice[1]].size++]=str->edges[_j].vertice[0];
    }
 
//Process STR
if ( (verbose))
  {
  if ((_i=protonate_str(str))==(unsigned int)-1) 
    { yprintf(YPRINTF_ERROR,"error in protonate_str(), err_code=%s.\n",get_yerrno(ylib_errno)); free(_neighbors), _neighbors=0x0; return FALSE; }
  else yprintf(YPRINTF_INFO,"protonate_str() has successfully added %d formal protons.\n",_i);
  if (!((*order)=(char(*)[4])define_neighbors_order(0x0,str->natoms,str->nedges,str->edges)))
    { yprintf(YPRINTF_ERROR,"error in define_neighbors_order(), err_code=%s.\n",get_yerrno(ylib_errno)); free(_neighbors), _neighbors=0x0; return FALSE; }
  if ((_i=update_NO2((*order),str))==(unsigned int)-1) 
    {  yprintf(YPRINTF_ERROR,"error in update_NO2(), err_code=%s.\n",get_yerrno(ylib_errno)); LABEL_ERROR_1: free(*order); (*order)=0x0; free(_neighbors), _neighbors=0x0; return FALSE; }
  else yprintf(YPRINTF_INFO,"update_NO2() has successfully edited %d NO2 groups.\n",_i);
  if ((_i=update_amides((*order),str))==(unsigned int)-1) 
    {  yprintf(YPRINTF_ERROR,"error in update_amides(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_1; }
  else yprintf(YPRINTF_INFO,"update_amides() has successfully edited %d amide bonds.\n",_i);
  if (pH<=13.6)
    {
    if ((_i=update_guanidines(order,_neighbors,str))==(unsigned int)-1) 
      {  yprintf(YPRINTF_ERROR,"error in update_guanidines(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_1; }
    else yprintf(YPRINTF_INFO,"requested pH is %.3f which is below 13.6 so update_guanidines() has been called and it successfully protonated %d guanidine groups.\n",pH,_i);
    }
  if (pH<=12.0)
    {
    if ((_i=update_amines(order,_neighbors,str))==(unsigned int)-1) 
      {  yprintf(YPRINTF_ERROR,"error in update_amines(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_1; }
    else yprintf(YPRINTF_INFO,"requested pH is %.3f, which is below 12.0 so update_amines() has been called and it successfully protonated %d amine groups.\n",pH,_i);
    }
  if (pH>=3.0)
    {
    if ((_i=update_carboxyles((*order),str))==(unsigned int)-1) 
      {  yprintf(YPRINTF_ERROR,"error in update_carboxyles(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_1; }
    else yprintf(YPRINTF_INFO,"requested pH is %.3f, which is above 3.0 so update_carboxyles() has been called and it successfully deprotonated %d carboxyl groups.\n",pH,_i);
    }
  if (pH>=4.0) 
    {
    if ((_i=update_silicates((*order),str))==(unsigned int)-1) 
      {  yprintf(YPRINTF_ERROR,"error in update_silicates(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_1; }
    else yprintf(YPRINTF_INFO,"requested pH is %.3f, which is above 4.0 so update_silicates() has been called and it successfully deprotonated %d silicate groups.\n",pH,_i);
    }  
  if (pH>=2.0)
    {
    if ((_i=update_phosphates((*order),str))==(unsigned int)-1) 
      {  yprintf(YPRINTF_ERROR,"error in update_phosphates(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_1; }
    else yprintf(YPRINTF_INFO,"requested pH is %.3f, which is above 2.0 so update_phosphates() has been called and it successfully deprotonated %d phosphate groups.\n",pH,_i);
    } 
  if (pH>=2.0)
    {
    if ((_i=update_sulfaites((*order),str))==(unsigned int)-1) 
      {  yprintf(YPRINTF_ERROR,"error in update_sulfaites(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_1; }
    else yprintf(YPRINTF_INFO,"requested pH is %.3f, which is above 2.0 so update_sulfaites() has been called and it successfully deprotonated %d sulf(a/i)te groups.\n",pH,_i);
    }
  if (!((*neighbors)=define_neighbors(FALSE,str->natoms,str->nedges,str->edges))) 
    {  yprintf(YPRINTF_ERROR,"error in define_neighbors(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_1; }
  else arrange_neighbors((*order),(*neighbors),str->natoms,str->nedges,str->edges);
  if ((_i=update_double_bonds((*order),(*neighbors),str,top))==(unsigned int)-1)
    { 
    yprintf(YPRINTF_ERROR,"error in update_double_bonds(), err_code=%s.\n",get_yerrno(ylib_errno));
    LABEL_ERROR_2: free(neighbors); neighbors=0x0; goto LABEL_ERROR_1;
    }
  else yprintf(YPRINTF_INFO,"[NOTE. It's yet NOT IMPLEMENTED] update_double_bonds() successfully edited %d edges in the molecular graph.\n",pH,_i);
  }
else
  {
  if ((_i=protonate_str(str))==(unsigned int)-1) { free(_neighbors), _neighbors=0x0; return FALSE; }
  if (!((*order)=(char(*)[4])define_neighbors_order(0x0,str->natoms,str->nedges,str->edges))) { free(_neighbors), _neighbors=0x0; return FALSE; }
  if ((_i=update_NO2((*order),str))==(unsigned int)-1) goto LABEL_ERROR_1;
  if ((_i=update_amides((*order),str))==(unsigned int)-1) goto LABEL_ERROR_1;
  if (pH<=13.6)
    {
    if ((_i=update_guanidines(order,_neighbors,str))==(unsigned int)-1) goto LABEL_ERROR_1;
    }
  if (pH<=12.0)
    {
    if ((_i=update_amines(order,_neighbors,str))==(unsigned int)-1) goto LABEL_ERROR_1;
    }
  if (pH>=3.5)
    {
    if ((_i=update_carboxyles((*order),str))==(unsigned int)-1) goto LABEL_ERROR_1;
    }
  if (pH>=2.0)
    {
    if ((_i=update_phosphates((*order),str))==(unsigned int)-1) goto LABEL_ERROR_1;
    } 
  if (pH>=4.0) 
    {
    if ((_i=update_silicates((*order),str))==(unsigned int)-1) goto LABEL_ERROR_1;
    }  
  if (pH>=2.0)
    {
    if ((_i=update_sulfaites((*order),str))==(unsigned int)-1) goto LABEL_ERROR_1;
    }
  if (!((*neighbors)=define_neighbors(FALSE,str->natoms,str->nedges,str->edges))) goto LABEL_ERROR_1;
  else arrange_neighbors((*order),(*neighbors),str->natoms,str->nedges,str->edges);
  if ((_i=update_double_bonds((*order),(*neighbors),str,top))==(unsigned int)-1) goto LABEL_ERROR_2;
  }

//Sync memory
if (!(vp=realloc(str->a,str->natoms*sizeof(char)))) goto LABEL_MEMORY_ERROR;
else str->a=(char*)vp;
if (!(vp=realloc(str->anames,str->natoms*sizeof(int)))) goto LABEL_MEMORY_ERROR;
else str->anames=(char (*)[sizeof(int)])vp; 
if (!(vp=realloc(str->r,str->natoms*sizeof(t_vec)))) goto LABEL_MEMORY_ERROR;
else str->r=(t_vec*)vp; 
if (!(vp=realloc(str->edges,str->nedges*sizeof(t_edge)))) goto LABEL_MEMORY_ERROR;
else str->edges=(t_edge*)vp; 
free(_neighbors), _neighbors=0x0;

return TRUE;
}

//This routine calculates brutto formula of a compound (fragment) in form of a text string
unsigned int calc_brutto_f(char **brt_f,unsigned int natoms,char *a)
{
register unsigned int _i, _j;
unsigned int item, ch[0x7F]={0};
size_t len;

while (natoms--) if (a[natoms]>0) ch[(unsigned int)a[natoms]]++; else ch[0]++;
natoms=0, _i=0x7F; while (_i--) if ( (ch[_i])) natoms++;
for (len=0x0, _i=_j=0; _j!=natoms && _i<0x7F; _i++) if ( (ch[_i])) { item=chemid_to_name(_i); len+=snprintf(0x0,0,"%s%d",(const char*)&item,ch[_i]); _j++; }
if (!((*brt_f)=(char*)malloc(sizeof(char)*(len+1)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
for (len=0x0, _i=_j=0; _j!=natoms && _i<0x7F; _i++) if ( (ch[_i])) { item=chemid_to_name(_i); len+= sprintf(&(*brt_f)[len],"%s%d",(const char*)&item,ch[_i]); _j++; }
return len+1;
}





/************************************      M O L      C O M P I L A T I O N     P A R T        *****************************************************/





//It is just a routine for the following function 
inline void _shift_not_aromatic_cycle(register unsigned int _i,register unsigned int *size_al,register t_clist *cycles)
{
register unsigned int *_ip, _t;
if ((*size_al)--!=_i)
  {
   _t=cycles->list[_i].size, cycles->list[_i].size=cycles->list[*size_al].size, cycles->list[*size_al].size=_t;
  _ip=cycles->list[_i].list, cycles->list[_i].list=cycles->list[*size_al].list, cycles->list[*size_al].list=_ip;
  }
}
//This function determines aromacity of cycles.
unsigned int classify_cycles(char (*order)[4],t_clist *neighbors,t_clist *cycles,t_str *str)
{
register unsigned int _i, _j, _k, _l, size_ar;
unsigned int *cg, *ca, ne, np, size_al;

//Stage 1. The first run over cycles.
// Eliminate definitely non-aromatic cycles and pop those are definitely aromatic
// A cycle is *definitely* non-aromatic if the cycle:
// a) a cycle of size >8 (due to  steric strain and angular strain) or contain atom with triple bond or atom with valence count >4 or neighbors >3
// b) contain d-elements
// c) contain an atom in the cycle that is either non-resonance (for carbon, silicon or phosphor) or non-resonance (nitrogen, oxigen or sulfur) which is not perfectly conjugated from both sides
// d) if one cycle enclose another one completely then remove a bigger cycle from considerations
//A cycle that is not *definitely* non-aromatic and has electron count 4n+2 is then aromatic
_i=0; size_al=cycles->size;
while (_i!=size_al)
  if ( ((_j=cycles->list[_i].size)>8)||(_j<5) )
    LABEL_SHIFT_NOT_AROMATIC_CYCLE: _shift_not_aromatic_cycle(_i,&size_al,cycles); //Move non-aromatic cycle to the end
  else
    {
    while (_j--)
      {
      _l=cycles->list[_i].list[_j];
      //Rule a) and b)
           if ( ( (order[_l][3]))||(order[_l][0]+order[_l][1]+2*order[_l][2]>4)||(order[_l][2]>1)||(neighbors->list[_l].size>3) ) goto LABEL_SHIFT_NOT_AROMATIC_CYCLE;
        //Rule c)
      else if (!order[_l][2])
             {
                  if ( (str->a[_l]==CHEM_ATOM_TYPE_CARBON)||(str->a[_l]==CHEM_ATOM_TYPE_SILICON)||(str->a[_l]==CHEM_ATOM_TYPE_PHOSPHOR) ) goto LABEL_SHIFT_NOT_AROMATIC_CYCLE;
             else if ( (str->a[_l]==CHEM_ATOM_TYPE_NITROGEN)||(str->a[_l]==CHEM_ATOM_TYPE_OXYGEN)||(str->a[_l]==CHEM_ATOM_TYPE_SULFUR) )
                    {
                    _k=neighbors->list[_l].size;
                    while (_k--)
                      if ( ((find_in_list(neighbors->list[_l].list[_k],&cycles->list[_i])!=(unsigned int)-1))&&(!(order[neighbors->list[_l].list[_k]][2])) ) goto LABEL_SHIFT_NOT_AROMATIC_CYCLE;
                    }
             }
      }
    _i++;
    }

//Stage 2. Run over conjugates
//Stage 2.0. Prepare some memory and init cycles->size
if (!(cg=(unsigned int*)malloc(sizeof(unsigned int)*size_al))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return (unsigned int)-1; }
_l=size_al, _k=0; while (_l--) { _k+=cycles->list[_l].size, cycles->list[_l].size=-(int)cycles->list[_l].size; }
if (_k>str->natoms) _k=str->natoms;
if (!(ca=(unsigned int*)malloc(sizeof(unsigned int)*_k))) { free(cg), cg=0x0; goto LABEL_MEMORY_ERROR; }

//Stage 2.1. Generate pi-conjugates 
//Note. Previously determined aromatic cycles are marked with sign bit of size
for (_i=0; _i<size_al; _i++)
  if ((int)(cycles->list[_i].size)<0)
    {
    //Stage 2.1.1. Construct a conjugate.
    _k=0, size_ar=1, *cg=_i, cycles->list[_i].size=-(int)cycles->list[_i].size;
    do{
      _j=size_al;
      while (--_j!=_i)
        if ((int)(cycles->list[_j].size)<0)
          {//Scan smaller-in-bigger to save on 'enclosure' test
          if (cycles->list[cg[_k]].size>-(int)cycles->list[_j].size)
            {
            _l=-(int)cycles->list[_j].size;
            while (--_l) //We are looking for 2-atoms intersection
              if (find_in_row(cycles->list[_j].list[_l],cycles->list[cg[_k]].size,cycles->list[cg[_k]].list)!=(unsigned int)-1)
                {
                CONJUGATE_JTH_CYCLE: cg[size_ar++]=_j, cycles->list[_j].size=-(int)cycles->list[_j].size;     //Reset conjugate member
                break;
                }
            }
          else
            {
            _l=cycles->list[cg[_k]].size;
            while (--_l) //We are looking for 2-atoms intersection
              if (find_in_row(cycles->list[cg[_k]].list[_l],-(int)cycles->list[_j].size,cycles->list[_j].list)!=(unsigned int)-1)
                goto CONJUGATE_JTH_CYCLE;
            }
          }
      }while (++_k!=size_ar);

    //Stage 2.1.2. Analyse the i-th conjugate
    if (size_ar!=1)   
      {//Stage 2.1.2.1. Create conjugates' unique atom list
      _j=_l=cycles->list[*cg].size; while (_l--) ca[_l]=cycles->list[*cg].list[_l];
      while (--_k) 
        { 
        _l=cycles->list[cg[_k]].size; 
        while (_l--)
          if (find_in_row(cycles->list[cg[_k]].list[_l],_j,ca)==(unsigned int)-1)
            ca[_j++]=cycles->list[cg[_k]].list[_l];
        }
      //Stage 2.1.2.2. Summ conjugates' electrons
      np=0; _l=ne=_j;
      while (_l--) 
        if (!(order[ca[_l]][2])) ne++; //if there is a perfectly conjugated sp3 element or external double bond then add extra electron
        else if (find_in_row(*neighbors->list[ca[_l]].list,_j,ca)==(unsigned int)-1) np++;
      switch (ne%4)
        {
        case 0 : {//Donate/drain 2e from/into external pi-system 
          if (np>=2) goto NEXT_CYCLE_I;
          break; }
        case 1 : {//Donate 1e into external pi-system
          if (np>=1) goto NEXT_CYCLE_I;
          break; }
        case 2 : {//Just aromaic conjugate itself
          goto NEXT_CYCLE_I; 
          break; }
        case 3 : {//Drain 1e from external pi-system 
          if (np>=1) goto NEXT_CYCLE_I;
          break; }
        }
      }

    //The whole conjugate is not aromatic, parse all the cycles
    //Stage 2.1.2.3. Consider elementary cycles which the conjugate consists of to save them if it's possible
    _k=size_ar;
    while (_k--)
      {
      np=0, _j=ne=cycles->list[cg[_k]].size;
      while (_j--)
        {
        _l=cycles->list[cg[_k]].list[_j];   
        if (!(order[_l][2])) ne++;
        else if (find_in_row(*neighbors->list[_l].list,cycles->list[cg[_k]].size,cycles->list[cg[_k]].list)==(unsigned int)-1) np++;
        }
      switch (ne%4)
        {
        case 0 : {//Donate/drain 2e from/into external pi-system 
          if (np<2) 
            {//Check if size_al is in the conjugate and update the id is it is so
            _shift_not_aromatic_cycle(cg[_k],&size_al,cycles);
            _l=_k; while (_l--) if (cg[_l]==size_al) { cg[_l]=cg[_k]; break; }
            }  
          break; }
        case 1 : {//Donate 1e into external pi-system
          if (!np)
            {
            _shift_not_aromatic_cycle(cg[_k],&size_al,cycles);
            _l=_k; while (_l--) if (cg[_l]==size_al) { cg[_l]=cg[_k]; break; }
            }  
          break; }
        case 2 : {//Just aromaic itself
          break; }
        case 3 : {//Drain 1e from external pi-system 
          if (!np)
            {
            _shift_not_aromatic_cycle(cg[_k],&size_al,cycles);
            _l=_k; while (_l--) if (cg[_l]==size_al) { cg[_l]=cg[_k]; break; }
            }
          break; }
        }
      } 
    NEXT_CYCLE_I: ;
    }

//Stage 2.2. Free the memory
free(cg), cg=0x0;
free(ca), ca=0x0;

//Exit
return size_al;
}


//*************   T H E     F I R S T      O R D E R       C O M P I L E R   ************/


//This function define YFF1 atom types of mol from chemical description
t_mol *compile_mol_YFF1(char (*order)[4],t_clist *neighbors,unsigned int size_ar,t_clist *cycles,t_str *str,t_top *top)
{
t_mol *mol;
register unsigned int _i, _j, _k;
unsigned int v_num, n_num, r_num;
int *c_num;

//Stage 0. Create empty mol
if (!(mol=(t_mol*)calloc(0x1,sizeof(t_mol))))                               { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return FALSE; }
else mol->size_ar=size_ar;
if (!(mol->ytypes=(unsigned int*)calloc(str->natoms,sizeof(unsigned int)))) { LABEL_MEMORY_ERROR_1: free(mol);         mol=0x0;         goto LABEL_MEMORY_ERROR_0; }
if (!(c_num=(int*)calloc(str->natoms,sizeof(int))))                         { LABEL_MEMORY_ERROR_2: free(mol->ytypes); mol->ytypes=0x0; goto LABEL_MEMORY_ERROR_1; }

//Stage I. Define atoms global c_num 
_i=cycles->size; 
while (size_ar!=_i) 
  {
  _j=cycles->list[--_i].size;
       if (_j==3) while (_j--) { { _k=cycles->list[_i].list[_j]; if ( (!c_num[_k])||(c_num[_k]==4) ) c_num[_k]=3; } }
  else if (_j==4) while (_j--) { { _k=cycles->list[_i].list[_j]; if   (!c_num[_k])                   c_num[_k]=4; } }
  } 
while (_i)
  {
  _j=cycles->list[--_i].size;
  if (_j==5) { while (_j--) { c_num[cycles->list[_i].list[_j]]=-5; } }
  else       { while (_j--) {    _k=cycles->list[_i].list[_j]; if (c_num[_k]!=-5) c_num[_k]=-6; } }
  }

//Stage II. Define heavy atoms types first
_i=str->natoms;
while (_i--)
  if (str->a[_i]!=CHEM_ATOM_TYPE_HYDROGEN) //hydrogens will be compiled later on the base of the root heavy atom
    {
    //Calculate atoms' descriptors
    v_num=3*order[_i][3]+2*order[_i][2]+order[_i][1]+order[_i][0];
    n_num=neighbors->list[_i].size;
    if ( (!(order[_i][2]))&&(!(order[_i][3]))&&(!(order[_i][1])) )
      { r_num=0, _k=n_num; while (_k--) { _j=neighbors->list[_i].list[_k]; if ( ( (order[_j][2]))||( (order[_j][3]))||( (order[_j][1])) )  { r_num-= 1; break; } } }
    else r_num=+1;
    //Find appropriate atom accordingly to calculated descriptors. Do two searches, one with the actual cycles descriptor and the second with the assumption c_num==0
    _j=top->size_a;
    while (_j--)
      if ( (top->ff_a[_j].chem_id==str->a[_i])&&(n_num==top->ff_a[_j].nn)&&(v_num==top->ff_a[_j].vn)&&(r_num==top->ff_a[_j].rn)&&(c_num[_i]==top->ff_a[_j].cycle) )
        {
        mol->ytypes[_i]=_j;
        goto NEXT_I;
        }
    _j=top->size_a;
    while (_j--)
      if ( (top->ff_a[_j].chem_id==str->a[_i])&&(n_num==top->ff_a[_j].nn)&&(v_num==top->ff_a[_j].vn)&&(r_num==top->ff_a[_j].rn)&&(!top->ff_a[_j].cycle) )
        {
        mol->ytypes[_i]=_j;
        goto NEXT_I;
        }
    LABEL_DATA_CONSISTMENT: 
    yprintf(YPRINTF_ERROR,"ERROR. No atom type defined for atom %d! (ch_type=%2d n_num=%1d v_num=%1d r_num=%1d c_num=%1d)\n",_i,str->a[_i],n_num,v_num,r_num,c_num[_i]);
    ylib_errno=YERROR_DATA_CONSISTMENT; free(mol->ytypes); free(mol); free(c_num); return FALSE;
    NEXT_I: ;
    }

//Stage III. Define hydrogens types then
_i=str->natoms;
while (_i--)
  if (str->a[_i]==CHEM_ATOM_TYPE_HYDROGEN) //hydrogens are compiled depending on the root heavy atom
    {
         if (neighbors->list[_i].size==1)
           {  
           _j=*neighbors->list[_i].list;
           r_num=top->ff_a[mol->ytypes[_j]].rn;
           //Find appropriate atom accordingly to calculated descriptors
           _k=top->size_a;
           while (_k--)
             if ( (top->ff_a[_k].chem_id==-(int)str->a[_j])&&(top->ff_a[_k].rn==r_num) )  
               {
               mol->ytypes[_i]=_k; 
               goto NEXT_HI;
               }
           yprintf(YPRINTF_ERROR,"ERROR. No atom type defined for hydrogen %d! (neighbour ch_type=%2d, r_num=%1d)\n",_i,-(int)str->a[_k],r_num); 
           goto LABEL_DATA_CONSISTMENT;
           }
    else if (!neighbors->list[_i].size)
           {
           if (str->natoms==1) 
             {//It's a hydrogen ion
             r_num=0, _k=top->size_a;
             while (_k--)
               if ( (top->ff_a[_k].chem_id==1)&&(top->ff_a[_k].rn==r_num) )  
                 {
                 mol->ytypes[_i]=_k;
                 goto NEXT_HI;
                 }
             yprintf(YPRINTF_ERROR,"ERROR. No atom type defined for hydrogen %d! (neighbour ch_type=1, r_num=%1d)\n",_i,r_num); 
             goto LABEL_DATA_CONSISTMENT;
             }
           else
             {
             yprintf(YPRINTF_ERROR,"ERROR. A hydrogen (%d), that is not attached to molecular graph, is detected.\n",_i);
             goto LABEL_DATA_CONSISTMENT;
             }
           }
    else   {
           yprintf(YPRINTF_ERROR,"ERROR. YFFx - is a proton-reactive FF family, i.e. there is no atom type is reserved for multivalent hydrogen (like %d), sorry.\n",_i);
           goto LABEL_DATA_CONSISTMENT;
           }
    NEXT_HI: ;
    }
    
//Stage IV. Move data from str into mol
mol->name=str->name;
mol->start_rid=str->start_rid;
if (!(str->ress))
  {//Compressed form of str
  if (!(mol->ress=(t_list*)alloc_list(0x1))) { LABEL_MEMORY_ERROR_3: free(c_num); c_num=0x0; goto LABEL_MEMORY_ERROR_2; }
  else { mol->ress->size=1, mol->ress->list[0]=*((unsigned int*)&str->rsize); }
  if (!(mol->rsize=(unsigned int*)malloc(sizeof(unsigned int)*2))) { free(mol->ress);   mol->ress=0x0;   goto LABEL_MEMORY_ERROR_3; }
  else { mol->rsize[0]=0, mol->rsize[1]=str->natoms; }
  }
else
  {//Full form of str
  mol->ress=str->ress;
  mol->rsize=str->rsize;
  }
mol->natoms=str->natoms;
mol->a=str->a;
mol->anames=str->anames;
mol->r=str->r;
mol->nedges=str->nedges;
mol->edges=str->edges;

//Done. Sved cycles and exit.
free(c_num); c_num=0x0;
mol->cycles=cycles;
return mol;
}

//This function implements topological proton transfering on the base of empirical rules.
//Please note that actual molecular graph should be modified previously at structure level (see numeous functions above).
//Rule 1. R-[C==[NH2]2(+)] -e is bound to each nitrogen.
//Rule 2. R3-N(sp3)-H(+) -e is bound to nitrogen.
//Rule 3. R-[O=C-O(-)] +e is bound to each oxigen.
//Rule 4. R-[O=S(+6)-O(-)] +e is bound to each oxigen.
//Rule 5. R-[O=S(+4)-O(-)] +e is bound to each oxigen
//Rule 6. R-[O=P(+5)-O(-)n] +ne are bound to each oxigen.
//Rule 7. R-[O=Si(+4)-O(-)n] +ne are bound to each oxigen.
char ionize_mol_empirically(char (*order)[4],t_clist *neighbors,t_mol *mol,t_top *top)
{
register unsigned int _i, _j, _k, _l;
unsigned char O1_CRG, N4_CRG, N3_CRG;
void *vp;

//Stage 0. Allocate some memory and define key atoms types
if ( (mol->vedges)) { free(mol->vedges); mol->vedges=0x0; }
if (!(mol->vedges=(unsigned int(*)[2])malloc(sizeof(unsigned int)*2*0xFF))) { LABEL_MEMORY_ERROR_1: ylib_errno=YERROR_MEMORY; return FALSE; }
if ( (mol->vatoms)) { free(mol->vatoms); mol->vatoms=0x0; }
if (!(mol->vatoms=(int*)malloc(sizeof(int)*0xFF))) { LABEL_MEMORY_ERROR_2: free(mol->vedges); mol->vedges=0x0; goto LABEL_MEMORY_ERROR_1; }

O1_CRG=N4_CRG=N3_CRG=0, _i=top->size_a;
while (_i--)
  {
       if (top->ff_a[_i].chem_id==CHEM_ATOM_TYPE_NITROGEN)
         {//N
              if ( (top->ff_a[_i].nn==4)&&(top->ff_a[_i].rn==+0)&&(top->ff_a[_i].vn==4)&&(!(top->ff_a[_i].cycle)) ) N4_CRG=_i;
         else if ( (top->ff_a[_i].nn==3)&&(top->ff_a[_i].rn==+1)&&(top->ff_a[_i].vn==4)&&(!(top->ff_a[_i].cycle)) ) N3_CRG=_i;
         }
  else if (top->ff_a[_i].chem_id==CHEM_ATOM_TYPE_OXYGEN)
         {//O
              if ( (top->ff_a[_i].nn==1)&&(top->ff_a[_i].rn==+1)&&(top->ff_a[_i].vn==1)&&(!(top->ff_a[_i].cycle)) ) O1_CRG=_i;
         }
  }
if ( (!O1_CRG)||(!N4_CRG)||(!N3_CRG) ) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }

//Stage I. Ionizing. Note sign bit is used here to avoid atoms reconsideing
//Rule 1. R-[C==[NH2]2(+)] -e is bound to each nitrogen. 
_i=mol->natoms;
while (_i--)
  if (mol->ytypes[_i]==N3_CRG)
    {
    //Charged atom detected - add one counterion anyway
    if (!(mol->nvedges%0xFF))
      {
      if (!(vp=(void*)realloc(mol->vedges,sizeof(unsigned int)*2*(mol->nvedges+0xFF)))) { LABEL_MEMORY_ERROR: free(mol->vatoms); mol->vatoms=0x0; goto LABEL_MEMORY_ERROR_2; }
      else mol->vedges=(unsigned int(*)[2])vp;
      }
    if (!(mol->nvatoms%0xFF))
      {
      if (!(vp=(void*)realloc(mol->vatoms,sizeof(int)*(mol->nvatoms+0xFF)))) goto LABEL_MEMORY_ERROR;
      else mol->vatoms=(int*)vp;
      }
    mol->vedges[mol->nvedges][0]=_i;
    mol->vedges[mol->nvedges][1]=mol->nvatoms;
    mol->nvedges++;
    mol->vatoms[mol->nvatoms]=-1;
    mol->nvatoms++;
    //In addition...
    if ( (!(order[_i][3]))&&(order[_i][2]==1)&&(!(order[_i][1]))&&(order[_i][0]==2)&&
         (mol->a[neighbors->list[_i].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mol->a[neighbors->list[_i].list[2]]==CHEM_ATOM_TYPE_HYDROGEN)&&
         (!(order[_j=*neighbors->list[_i].list][3]))&&(order[_j][2]==1)&&( (order[_j][1]))&&(order[_j][1]+order[_j][0]==2) )
      {
      _k=neighbors->list[_j].list[1], _l=neighbors->list[_j].list[2];
      if ( (mol->a[_k]==CHEM_ATOM_TYPE_NITROGEN)&&(!(order[_k][3]))&&(!(order[_k][2]))&&(order[_k][1]==1)&&(order[_k][0]==2)&&
           (mol->a[neighbors->list[_k].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mol->a[neighbors->list[_k].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) )
        {
        if ( (order[_j][1]==2)&&(mol->a[_l]==CHEM_ATOM_TYPE_NITROGEN)&&(!(order[_l][3]))&&(!(order[_l][2]))&&(order[_l][1]==1)&&(order[_l][0]==2)&&
             (mol->a[neighbors->list[_l].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mol->a[neighbors->list[_l].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) )
          {//H2N-a-[H2N=C]-a-NH2
          mol->ytypes[_k]=-(int)N3_CRG, mol->ytypes[_l]=-(int)N3_CRG;
          mol->vedges[mol->nvedges-1][0]=_j;
          }
        else
          {//R-[H2N=C]-a-NH2
          EDIT_AMIDINE: mol->ytypes[_k]=-(int)N3_CRG;
          mol->vedges[mol->nvedges-1][0]=_j; //Recomutate counterion to the central atom of the group
          }
        }
      else if ( (order[_j][1]==2)&&(mol->a[_l]==CHEM_ATOM_TYPE_NITROGEN)&&(!(order[_l][3]))&&(!(order[_l][2]))&&(order[_l][1]==1)&&(order[_l][0]==2)&&
                (mol->a[neighbors->list[_l].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mol->a[neighbors->list[_l].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) )
             { _k=_l; goto EDIT_AMIDINE; }
      }
    }
_i=mol->natoms; while (_i--) if ((int)mol->ytypes[_i]==-(int)N3_CRG) mol->ytypes[_i]=N3_CRG; 
//Rule 2. R-[X-N-Y] -> R-[X-N-Y]-H(+) +e is bound to nitrogen.
_i=mol->natoms;
while (_i--)
  if (mol->ytypes[_i]==N4_CRG)
    {//Charged atom detected - add counterion anyway
    if (!(mol->nvedges%0xFF))
      {
      if (!(vp=(void*)realloc(mol->vedges,sizeof(unsigned int)*2*(mol->nvedges+0xFF))))  goto LABEL_MEMORY_ERROR;
      else mol->vedges=(unsigned int(*)[2])vp;
      }
    if (!(mol->nvatoms%0xFF))
      {
      if (!(vp=(void*)realloc(mol->vatoms,sizeof(int)*(mol->nvatoms+0xFF))))    goto LABEL_MEMORY_ERROR;
      else mol->vatoms=(int*)vp;
      }
    mol->vedges[mol->nvedges][0]=_i;
    mol->vedges[mol->nvedges][1]=mol->nvatoms;
    mol->nvedges++;
    mol->vatoms[mol->nvatoms]=-1;
    mol->nvatoms++;
    }
//Rule 3-7. R-[O=X=O(-)n] +ne are bound to each oxigen.
_i=mol->natoms;
while (_i--)
  if ( (mol->ytypes[_i]==O1_CRG)&&(neighbors->list[_i].size==1) )
    {//Join gegnion to the central group atom
    _j=*neighbors->list[_i].list;
    //Check if there is the multivalent counterions (like PO3(2-))
    _k=mol->nvedges;
    while (_k--) 
      if ( (mol->vedges[_k][0]==_j)&&(mol->vatoms[mol->vedges[_k][1]]<0.) )
        {
        mol->vatoms[mol->vedges[_k][1]]--;
        goto NEXT_ACID_I;
        }
    //Charged atom detected - add counterion anyway
    if (!(mol->nvedges%0xFF))
      {
      if (!(vp=(void*)realloc(mol->vedges,sizeof(unsigned int)*2*(mol->nvedges+0xFF))))  goto LABEL_MEMORY_ERROR;
      else mol->vedges=(unsigned int(*)[2])vp;
      }
    if (!(mol->nvatoms%0xFF))
      {
      if (!(vp=(void*)realloc(mol->vatoms,sizeof(int)*(mol->nvatoms+0xFF))))    goto LABEL_MEMORY_ERROR;
      else mol->vatoms=(int*)vp;
      }
    mol->vedges[mol->nvedges][0]=_j;
    mol->vedges[mol->nvedges][1]=mol->nvatoms;
    mol->nvedges++;
    mol->vatoms[mol->nvatoms]=+1;
    mol->nvatoms++;
    //In addition symmeterize oxigens around the group
    _k=(unsigned int)order[_j][2]; 
    while (_k--)
      if ( (neighbors->list[neighbors->list[_j].list[_k+order[_j][3]]].size==1)&&
           (mol->a[neighbors->list[_j].list[_k+order[_j][3]]]==CHEM_ATOM_TYPE_OXYGEN) )
        mol->ytypes[neighbors->list[_j].list[_k+order[_j][3]]]=-(int)O1_CRG; //The minus to avoid multiple entrances on this atom
    NEXT_ACID_I: ;
    }
_i=mol->natoms; while (_i--) if ((int)mol->ytypes[_i]==-(int)O1_CRG) mol->ytypes[_i]=O1_CRG; 

//Stage III. Sync memory and calc charges at last
if ( (mol->nvatoms))
  {
  if (!(vp=(void*)realloc(mol->vatoms,sizeof(int)*mol->nvatoms)))    goto LABEL_MEMORY_ERROR;
  else mol->vatoms=(int*)vp;
  if (!(vp=(void*)realloc(mol->vedges,sizeof(unsigned int)*2*mol->nvedges)))  goto LABEL_MEMORY_ERROR;
  else mol->vedges=(unsigned int(*)[2])vp;
  }
else { free(mol->vatoms), mol->vatoms=0x0; free(mol->vedges), mol->vedges=0x0; }
if ( (mol->cmtype))
  {
       if (mol->cmtype==-1) { if (mol->C.dL) { free(mol->C.dL);         mol->C.dL=0x0; } }
  else if (mol->cmtype==+1) { if (mol->C.sL) { free_smatrix(mol->C.sL); mol->C.sL=0x0; } }
  else { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
  mol->cmtype=0;
  }

//Done
return calc_Oliferenko_ionized_CDS_cholesky_default(0xFF,SMALL2*sqrt((double)mol->nvatoms),mol,top);
}


//*************   T H E     S E C O N D      O R D E R       C O M P I L E R S   ************/


//This function shows created anchors
void show_anchors(register t_mol *mol)
{
register unsigned int _j,_i;
_i=mol->anchors->size;
while(_i--)
  {
  printf("Anchor %1d:  ",_i);
  _j=mol->anchors->list[_i].size;
  while (_j--)
    printf("%1d ",mol->anchors->list[_i].list[_j]+1);
  printf("\n");
  }
}

//This function defines topological part of mol - anchors and their graph. Neighbors hyper-list is used.
char disassemble_mol(unsigned int **_anchor_id,t_clist *neighbors,t_mol *mol)
{
register unsigned int _i, _j, _k, _l, _t;
unsigned int *anchor_id, nanchors=0;

//Stage 0. Prepare memory
if (!(anchor_id=(unsigned int*)malloc(sizeof(unsigned int)*mol->natoms)))  { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return FALSE; }
_l=mol->natoms; while (_l--) anchor_id[_l]=(unsigned int)-1;

//Stage 1. Move cycles into separate anchors (joining them if necessary)
for (_i=0; _i<mol->cycles->size; _i++)
  if ( ((_k=mol->cycles->list[_i].size)<5)||(_i<mol->size_ar) )
    {//Move cycle into separate anchor
    _j=nanchors;
    while (_k--)
      {
      if ( ((_t=anchor_id[mol->cycles->list[_i].list[_k]])!=(unsigned int)-1)&&(_t!=_j) )
        {//Do joining
        if (_t<_j) { _l=_t, _t=_j, _j=_l; } //Set _j to assignance and _t to replacement 
        _l=mol->natoms; while (_l--) if (anchor_id[_l]==_t) anchor_id[_l]=_j; //Set all already known to (new) _j
        if (_t!=nanchors--)
          { _l=mol->natoms; while (_l--) if (anchor_id[_l]==nanchors) anchor_id[_l]=_t; } //Shift the most recent anchor
        }
      anchor_id[mol->cycles->list[_i].list[_k]]=_j;
      }
    nanchors++;
    }
//Stage 2. Move all resonance bonds into separate anchors (joining them if necessary)
_k=mol->nedges;
while (_k--)
  if (mol->edges[_k].type!='1')
    {//Move resonance bond into separate anchor
         if ((_i=anchor_id[mol->edges[_k].vertice[0]])==(unsigned int)-1)
           {
           if ((_j=anchor_id[mol->edges[_k].vertice[1]])==(unsigned int)-1)
             { anchor_id[mol->edges[_k].vertice[0]]=anchor_id[mol->edges[_k].vertice[1]]=nanchors++; }
           else anchor_id[mol->edges[_k].vertice[0]]=_j;
           }
    else if ((_j=anchor_id[mol->edges[_k].vertice[1]])==(unsigned int)-1)
           { anchor_id[mol->edges[_k].vertice[1]]=_i; }
    else if (_i!=_j)
           {//Do joining 
           if (_i<_j) { _l=_i, _i=_j, _j=_l; } ///Set _j to assignance and _i to replacement 
           _l=mol->natoms; while (_l--) if (anchor_id[_l]==_i) anchor_id[_l]=_j; //Set all already known to (new) _j
           if (_i!=--nanchors) { _l=mol->natoms; while (_l--) if (anchor_id[_l]==nanchors) anchor_id[_l]=_i; } //Shift the most recent anchor
           }
    }

//Stage 3. Move multivalent to separate anchors
_i=mol->natoms; while (_i--) if ( (anchor_id[_i]==(unsigned int)-1)&&(neighbors->list[_i].size!=1) ) anchor_id[_i]=nanchors++;

//Stage 4. Move unsigned onevalent to corresponding anchors
_i=mol->natoms; while (_i--) if ( (anchor_id[_i]==(unsigned int)-1)&&(neighbors->list[_i].size==1) ) anchor_id[_i]=anchor_id[*neighbors->list[_i].list];

//Stage 5. Find the cases of two anchors connected with more than 1 bond and merge these anchors (like in CH2=(NH-C=O)2)
_t=mol->nedges; 
while (--_t)
  if ((_i=anchor_id[mol->edges[_t].vertice[0]])!=(_j=anchor_id[mol->edges[_t].vertice[1]]))
    {
    _k=_t;
    while (_k--) 
      if ( ( (_i==anchor_id[mol->edges[_k].vertice[0]])&&(_j==anchor_id[mol->edges[_k].vertice[1]]) )||
           ( (_j==anchor_id[mol->edges[_k].vertice[0]])&&(_i==anchor_id[mol->edges[_k].vertice[1]]) ) )
        {//Merge the anchors
        if (_i>_j) { _k=_i, _i=_j, _j=_k; }
        _l=mol->natoms; while (_l--) { if (anchor_id[_l]==_j) anchor_id[_l]=_i; else if (anchor_id[_l]>_j) anchor_id[_l]--; }
        nanchors--, _t=mol->nedges; break;
        }
    }

//Stage 6. Convert anchors into mol->anchors hyperlist and mol->aedges hypergraph
if ( (mol->anchors)) { free(mol->anchors); mol->anchors=0x0; } //Remove the old hyperlist
if (!(mol->anchors=alloc_clist(nanchors,mol->natoms))) { LABEL_MEMORY_ERROR_1: free(anchor_id); anchor_id=0x0; goto LABEL_MEMORY_ERROR_0; }
else
  {//Set it up!
  _i=mol->anchors->size=nanchors; while (_i--) mol->anchors->list[_i].size=0;
  _i=mol->natoms; while (_i--) mol->anchors->list[anchor_id[_i]].size++;
  for (mol->anchors->list[0].list=mol->anchors->_items,_i=1;_i<nanchors;_i++) mol->anchors->list[_i].list=mol->anchors->list[_i-1].list+mol->anchors->list[_i-1].size;
  _i=nanchors; while (_i--) mol->anchors->list[_i].size=0;
  _i=mol->natoms; while (_i--) mol->anchors->list[anchor_id[_i]].list[mol->anchors->list[anchor_id[_i]].size++]=_i;
  }
if ( (mol->aedges)) { free(mol->aedges); mol->aedges=0x0; } //Remove the old aedges
_i=mol->nedges, mol->naedges=0; while (_i--) if (anchor_id[mol->edges[_i].vertice[0]]!=anchor_id[mol->edges[_i].vertice[1]]) mol->naedges++;
if (!(mol->aedges=(t_edge*)malloc(sizeof(t_edge)*mol->naedges))) goto LABEL_MEMORY_ERROR_1;
else
  {//Setup aedges
  _i=mol->nedges, _j=0; 
  while (_i--)
    if (anchor_id[mol->edges[_i].vertice[0]]!=anchor_id[mol->edges[_i].vertice[1]]) 
      { mol->aedges[_j].vertice[0]=anchor_id[mol->edges[_i].vertice[0]], mol->aedges[_j].vertice[1]=anchor_id[mol->edges[_i].vertice[1]], mol->aedges[_j++].type=_i; }
  }

//Free memory and exit
if (!(*_anchor_id)) { free(anchor_id); anchor_id=0x0; } else (*_anchor_id)=anchor_id;
return TRUE;
}

//This is a service function to define the most suitable atom for torsion (one per anchor pair)
//One torsion per rotable bonds: atoms priority C.2 > C.1 > C.3 > N.2 > N.1 > N.3 > O.1 > O.2 > S.6 > S.4 > S.1 > S.2 > P.5 > P.2 > P.3 > Si.2 > Si.4 > I > Br > Cl > F > H > X
unsigned int _define_appropriate_torsion_atom(t_clist *neighbors,unsigned int atom1,unsigned int atom2,char *a)
{
register unsigned int _i, _j, _id=(unsigned int)-1;
_i=neighbors->list[atom1].size;
while (_i--)            
  if ((_j=neighbors->list[atom1].list[_i])!=atom2)
    {
    if (_id==(unsigned int)-1) _id=_j;
    else
      switch (a[_j])
        {
        case CHEM_ATOM_TYPE_CARBON   : {
            if (a[_id]==CHEM_ATOM_TYPE_CARBON)
              {// C.2 > C.1 > C.3 
                   if (neighbors->list[_j].size==2) { if (neighbors->list[_id].size!=2) _id=_j; }
              else if (neighbors->list[_j].size==1) { if (neighbors->list[_id].size!=1) _id=_j; }
              else if (neighbors->list[_j].size==3) { if (neighbors->list[_id].size!=3) _id=_j; }
              }
            else _id=_j;
          break; }
        case CHEM_ATOM_TYPE_NITROGEN : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)                                                                     )
            {
            if (a[_id]==CHEM_ATOM_TYPE_NITROGEN)
              {// N.2 > N.1 > N.3 
                   if (neighbors->list[_j].size==2) { if (neighbors->list[_id].size!=2) _id=_j; }
              else if (neighbors->list[_j].size==1) { if (neighbors->list[_id].size!=1) _id=_j; }
              else if (neighbors->list[_j].size==3) { if (neighbors->list[_id].size!=3) _id=_j; }
              } 
            else _id=_j;
            }
          break; }
        case CHEM_ATOM_TYPE_OXYGEN   : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)&&(a[_id]!=CHEM_ATOM_TYPE_NITROGEN)                                  )
            {
            if (a[_id]==CHEM_ATOM_TYPE_OXYGEN)
              {// O.1 > O.2 
                   if (neighbors->list[_j].size==1) { if (neighbors->list[_id].size!=1) _id=_j; }
              else if (neighbors->list[_j].size==2) { if (neighbors->list[_id].size!=2) _id=_j; }
              }
            else _id=_j;
            }
          break; }
        case CHEM_ATOM_TYPE_SULFUR   : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)&&(a[_id]!=CHEM_ATOM_TYPE_NITROGEN)&&(a[_id]!=CHEM_ATOM_TYPE_OXYGEN) )
            {
            if (a[_id]==CHEM_ATOM_TYPE_SULFUR)
              {// S.6 > S.4 > S.1 > S.2
                   if (neighbors->list[_j].size==4) { if (neighbors->list[_id].size!=4) _id=_j; }
              else if (neighbors->list[_j].size==3) { if (neighbors->list[_id].size!=3) _id=_j; }
              else if (neighbors->list[_j].size==1) { if (neighbors->list[_id].size!=1) _id=_j; }
              else if (neighbors->list[_j].size==2) { if (neighbors->list[_id].size!=2) _id=_j; }
              }
            else _id=_j;
            }
          break; }
        case CHEM_ATOM_TYPE_PHOSPHOR : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)&&(a[_id]!=CHEM_ATOM_TYPE_NITROGEN)&&(a[_id]!=CHEM_ATOM_TYPE_OXYGEN)&&
               (a[_id]!=CHEM_ATOM_TYPE_SULFUR)                                                                     )
            {
            if (a[_id]==CHEM_ATOM_TYPE_PHOSPHOR) 
              {// P.5 > P.2 > P.1 > P.3
                   if (neighbors->list[_j].size==4) { if (neighbors->list[_id].size!=4) _id=_j; }
              else if (neighbors->list[_j].size==2) { if (neighbors->list[_id].size!=2) _id=_j; }
              else if (neighbors->list[_j].size==1) { if (neighbors->list[_id].size!=1) _id=_j; }
              else if (neighbors->list[_j].size==3) { if (neighbors->list[_id].size!=3) _id=_j; }
              }
            else _id=_j;
            }
          break; }
        case CHEM_ATOM_TYPE_SILICON  : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)&&(a[_id]!=CHEM_ATOM_TYPE_NITROGEN)&&(a[_id]!=CHEM_ATOM_TYPE_OXYGEN)&&
               (a[_id]!=CHEM_ATOM_TYPE_SULFUR)&&(a[_id]!=CHEM_ATOM_TYPE_PHOSPHOR)                                  )
            {
            if (a[_id]!=CHEM_ATOM_TYPE_SILICON)
              {// Si.2 > Si.4
                   if (neighbors->list[_j].size==3) { if (neighbors->list[_id].size!=3) _id=_j; }
              else if (neighbors->list[_j].size==4) { if (neighbors->list[_id].size!=4) _id=_j; }
              }
            else _id=_j;
            }
          break; }
        case CHEM_ATOM_TYPE_IODINE   : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)&&(a[_id]!=CHEM_ATOM_TYPE_NITROGEN)&&(a[_id]!=CHEM_ATOM_TYPE_OXYGEN) &&
               (a[_id]!=CHEM_ATOM_TYPE_SULFUR)&&(a[_id]!=CHEM_ATOM_TYPE_PHOSPHOR)&&(a[_id]!=CHEM_ATOM_TYPE_SILICON) )
            {// Just I
            if (a[_id]!=CHEM_ATOM_TYPE_IODINE) 
              _id=_j;
            }
          break; }
        case CHEM_ATOM_TYPE_BROMINE  : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)&&(a[_id]!=CHEM_ATOM_TYPE_NITROGEN)&&(a[_id]!=CHEM_ATOM_TYPE_OXYGEN) &&
               (a[_id]!=CHEM_ATOM_TYPE_SULFUR)&&(a[_id]!=CHEM_ATOM_TYPE_PHOSPHOR)&&(a[_id]!=CHEM_ATOM_TYPE_SILICON)&&
               (a[_id]!=CHEM_ATOM_TYPE_IODINE)                                                                      )
            {// Just Br
            if (a[_id]!=CHEM_ATOM_TYPE_BROMINE)
              _id=_j;
            }
          break; }
        case CHEM_ATOM_TYPE_CHLORINE : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)&&(a[_id]!=CHEM_ATOM_TYPE_NITROGEN)&&(a[_id]!=CHEM_ATOM_TYPE_OXYGEN) &&
               (a[_id]!=CHEM_ATOM_TYPE_SULFUR)&&(a[_id]!=CHEM_ATOM_TYPE_PHOSPHOR)&&(a[_id]!=CHEM_ATOM_TYPE_SILICON)&&
               (a[_id]!=CHEM_ATOM_TYPE_IODINE)&&(a[_id]!=CHEM_ATOM_TYPE_BROMINE)                                    )
            {// Just Cl
            if (a[_id]!=CHEM_ATOM_TYPE_CHLORINE)
              _id=_j;
            }
          break; }
        case CHEM_ATOM_TYPE_FLUORINE : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)&&(a[_id]!=CHEM_ATOM_TYPE_NITROGEN)&&(a[_id]!=CHEM_ATOM_TYPE_OXYGEN)  &&
               (a[_id]!=CHEM_ATOM_TYPE_SULFUR)&&(a[_id]!=CHEM_ATOM_TYPE_PHOSPHOR)&&(a[_id]!=CHEM_ATOM_TYPE_SILICON) &&
               (a[_id]!=CHEM_ATOM_TYPE_IODINE)&&(a[_id]!=CHEM_ATOM_TYPE_BROMINE) &&(a[_id]!=CHEM_ATOM_TYPE_CHLORINE) )
            {// Just F
            if (a[_id]!=CHEM_ATOM_TYPE_FLUORINE)
              _id=_j;
            }
          break; }
        case CHEM_ATOM_TYPE_HYDROGEN : {
          if ( (a[_id]!=CHEM_ATOM_TYPE_CARBON)&&(a[_id]!=CHEM_ATOM_TYPE_NITROGEN)&&(a[_id]!=CHEM_ATOM_TYPE_OXYGEN)  &&
               (a[_id]!=CHEM_ATOM_TYPE_SULFUR)&&(a[_id]!=CHEM_ATOM_TYPE_PHOSPHOR)&&(a[_id]!=CHEM_ATOM_TYPE_SILICON) &&
               (a[_id]!=CHEM_ATOM_TYPE_IODINE)&&(a[_id]!=CHEM_ATOM_TYPE_BROMINE) &&(a[_id]!=CHEM_ATOM_TYPE_CHLORINE)&&(a[_id]!=CHEM_ATOM_TYPE_FLUORINE) )
            {// Just H
            if (a[_id]!=CHEM_ATOM_TYPE_HYDROGEN)
              _id=_j;
            }
          break; }
        }
    }
return _id;
}
//This function construct molecular mechanic force field of molecule. Neighbors and anchors hyper-lists are used and in_aromatic_cycle complex list is created.
//NB! We initialize impropers and dihedrals with dividing constant and aromatic cycles impropers with values - be aware not to overwrite them 
char compose_mol(char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_adjacency *adjacency,t_clist **_inacycles,t_mol *mol,t_top *top)
{
register unsigned int _i, _j, _k, _l, _n, _p, _q, _t; 
t_clist *inacycles;
void *vp;
//Stage 0. General tests
if ( (!(mol->natoms))||( (mol->size_b))||( (mol->size_g))||( (mol->size_i))||( (mol->size_d)) ) { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
if (mol->natoms==1) return TRUE; //No bonded topology for single-atom ion
else
  {//Build inacycles clist
  _i=mol->size_ar, _j=0; while (_i--) _j+=mol->cycles->list[_i].size;
  if (!(inacycles=(t_clist*)alloc_clist(mol->natoms,_j))) { ylib_errno=YERROR_MEMORY; return FALSE; }
  _i=mol->natoms; while (_i--) inacycles->list[_i].size=0;
  _i=mol->size_ar; while (_i--) { _j=mol->cycles->list[_i].size; while (_j--) inacycles->list[mol->cycles->list[_i].list[_j]].size++; }
  for (inacycles->list->list=inacycles->_items, _i=1; _i<mol->natoms; _i++) inacycles->list[_i].list=inacycles->list[_i-1].list+inacycles->list[_i-1].size;
  _i=mol->natoms; while (_i--) inacycles->list[_i].size=0;
  _i=mol->size_ar; while (_i--) { _j=mol->cycles->list[_i].size; while (_j--) { _k=mol->cycles->list[_i].list[_j], inacycles->list[_k].list[inacycles->list[_k].size++]=_i; } }  
  }
//Stage 1. Get statistics for bonded part
//Run over atoms
_i=mol->natoms;
while(_i--)
  if ((_j=neighbors->list[_i].size)!=1)
    {
    //Define angles
    mol->size_g+=_j*(_j-1)/2;
    //Define qualy-impropers
    if (_j==3)
      {
      if ( ( ( (order[_i][2]))||( (order[_i][3])) )||( (mol->a[_i]==CHEM_ATOM_TYPE_NITROGEN)&&( ( (inacycles->list[_i].size))||( ( (order[_i][1]))&&
         ( ( (order[neighbors->list[_i].list[0]][2]))||( (order[neighbors->list[_i].list[1]][2]))||( (order[neighbors->list[_i].list[2]][2])) ) ) ) ) )
        mol->size_d++; //_j*(_j-1)*(_j-2)/6;
      }
    //Define optical hydrogen-based impropers
    else
      if (_j==4)
        {
        _n=0;
        if (mol->a[neighbors->list[_i].list[0]]==CHEM_ATOM_TYPE_HYDROGEN) { _l=0, _n++; }
        if (mol->a[neighbors->list[_i].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) { _l=1, _n++; }
        if (mol->a[neighbors->list[_i].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) { _l=2, _n++; }
        if (mol->a[neighbors->list[_i].list[3]]==CHEM_ATOM_TYPE_HYDROGEN) { _l=3, _n++; }
        if (_n==1) mol->size_i++;
        }
    }
//Run over bonds
_t=mol->nedges;
while (_t--)
  if ( (neighbors->list[_i=mol->edges[_t].vertice[0]].size!=1)&&(neighbors->list[_j=mol->edges[_t].vertice[1]].size!=1) )
    {
         if (anchor_id[_i]!=anchor_id[_j]) //Rule 2.1. Cross-anchor
           {//Skip potential A-B-C-A dihedrals
           mol->size_d++;
           }
    else if ( (mol->edges[_t].type!='1')&&( (!(inacycles->list[_i].size))||(!(inacycles->list[_j].size)) ) ) //Rule 2.1. In anchor
           {
           //Rule 2.0. Amide bond with (O,S)=C-NH2
           if ( (mol->edges[_t].type==(int)'m')&&(neighbors->list[_i=mol->edges[_t].vertice[0]].size==3)&&(neighbors->list[_j=mol->edges[_t].vertice[1]].size==3) ) 
             { 
             if (mol->a[_i]==CHEM_ATOM_TYPE_NITROGEN) { _i=mol->edges[_t].vertice[1], _j=mol->edges[_t].vertice[0]; }
             if ( (mol->a[_i]==CHEM_ATOM_TYPE_CARBON)&&(!(order[_j][3]))&&(!(order[_j][2]))&&(neighbors->list[*neighbors->list[_i].list].size==1)&&
                ( (mol->a[*neighbors->list[_i].list]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[*neighbors->list[_i].list]==CHEM_ATOM_TYPE_SULFUR) )       &&
                ( (mol->a[_j]==CHEM_ATOM_TYPE_NITROGEN)&&(!(order[_j][3]))&&(!(order[_j][2]))&&( (order[_j][1])) )                               &&
                ( (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mol->a[neighbors->list[_j].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
               { mol->size_d+=2; continue; } //Symmetric case - use dihedral quasy-improper instead
             }
           //Skip potential A-B-C-A impropers
           if (neighbors->list[_i].size<neighbors->list[_j].size)
             { 
             _k=_l=neighbors->list[_i].size, _k--; 
             while (_l--) if ( (neighbors->list[_i].list[_l]!=_j)&&(find_in_list(neighbors->list[_i].list[_l],&neighbors->list[_j])!=(unsigned int)-1) ) _k--;
             _k*=(neighbors->list[_j].size-1);
             }
           else
             { 
             _k=_l=neighbors->list[_j].size, _k--; 
             while (_l--) if ( (neighbors->list[_j].list[_l]!=_i)&&(find_in_list(neighbors->list[_j].list[_l],&neighbors->list[_i])!=(unsigned int)-1) ) _k--;
             _k*=(neighbors->list[_i].size-1);
             }
           mol->size_i+=_k;
           }
    else if ( ( (inacycles->list[_i].size))&&( (inacycles->list[_j].size)) )
           {//Run over inacycles: screen all overlaps in aromatic cycles 
           _p=inacycles->list[_i].size;
           while (_p--)
             {
             _n=inacycles->list[_i].list[_p], _q=inacycles->list[_j].size;
             while (_q--)
               if (_n==inacycles->list[_j].list[_q])
                 {
                 _k=neighbors->list[_i].size;
                 while (_k--) 
                   if ( (_j!=neighbors->list[_i].list[_k])&&( (inacycles->list[neighbors->list[_i].list[_k]].size))&&(find_in_list(_n,&inacycles->list[neighbors->list[_i].list[_k]])!=(unsigned int)-1) )
                     {
                     _l=neighbors->list[_j].size;
                     while (_l--) 
                       if ( (_i!=neighbors->list[_j].list[_l])&&( (inacycles->list[neighbors->list[_j].list[_l]].size))&&(find_in_list(_n,&inacycles->list[neighbors->list[_j].list[_l]])!=(unsigned int)-1) )
                         mol->size_i++;
                     }
                 }
             }
           }
    }
//Stage 2. Alloc memory
if ( (mol->nedges)&&(!(mol->ff_b=(t_ff_b*)malloc(sizeof(t_ff_b)*mol->nedges))) ) { LABEL_MEMORY_ERROR_1: free(inacycles); ylib_errno=YERROR_MEMORY; return FALSE; }
else mol->size_b=mol->nedges;
if ( (mol->size_g)&&(!(mol->ff_g=(t_ff_g*)malloc(sizeof(t_ff_g)*mol->size_g))) ) { LABEL_MEMORY_ERROR_2: if (mol->size_b)  free(mol->ff_b); mol->ff_b=0x0; goto LABEL_MEMORY_ERROR_1; }
else mol->size_g=0x0;
if ( (mol->size_i)&&(!(mol->ff_i=(t_ff_i*)malloc(sizeof(t_ff_i)*mol->size_i))) ) { LABEL_MEMORY_ERROR_3: if (mol->size_g) free(mol->ff_g); mol->ff_g=0x0; goto LABEL_MEMORY_ERROR_2; }
else mol->size_i=0x0;
if ( (mol->size_d)&&(!(mol->ff_d=(t_ff_d*)malloc(sizeof(t_ff_d)*mol->size_d))) ) { LABEL_MEMORY_ERROR_4: if (mol->size_i) free(mol->ff_i); mol->ff_i=0x0; goto LABEL_MEMORY_ERROR_3; }
else mol->size_d=0x0;
//Stage 3. Generate list of bonded terms 
//Run over atoms
_i=mol->natoms;
while(_i--)
  if ((_j=neighbors->list[_i].size)!=1) 
    {
    //Define angles
    _k=_j;
    while (_k--)
      {
      _l=_k;
      while (_l--)
        {
        mol->ff_g[mol->size_g].atom[1]=_i, mol->ff_g[mol->size_g].atom[0]=neighbors->list[_i].list[_k], mol->ff_g[mol->size_g].atom[2]=neighbors->list[_i].list[_l];
        mol->ff_g[mol->size_g].k=mol->ff_g[mol->size_g].v=(double)NAN; 
        mol->size_g++;
        }
      }
    //Define quasy-impropers
    if (_j==3)
      {
      if ( ( ( (order[_i][2]))||( (order[_i][3])) )||( (mol->a[_i]==CHEM_ATOM_TYPE_NITROGEN)&&( ( (inacycles->list[_i].size))||( ( (order[_i][1]))&&
         ( ( (order[neighbors->list[_i].list[0]][2]))||( (order[neighbors->list[_i].list[1]][2]))||( (order[neighbors->list[_i].list[2]][2])) ) ) ) ) )
        {//Planar 3-center quasy-improper
        mol->ff_d[mol->size_d].atom[0]=_i, mol->ff_d[mol->size_d].atom[1]=neighbors->list[_i].list[0], mol->ff_d[mol->size_d].atom[2]=neighbors->list[_i].list[1], mol->ff_d[mol->size_d].atom[3]=neighbors->list[_i].list[2];
        mol->ff_d[mol->size_d].d=0., mol->ff_d[mol->size_d].k[0]=(double)NAN, mol->ff_d[mol->size_d].v[0]=PI, mol->ff_d[mol->size_d].n[0]=2.; 
        mol->size_d++;
        }
      }
    //Define optical hydrogen-based impropers
    else
      if (_j==4)
        {
        _n=0;
        if (mol->a[neighbors->list[_i].list[0]]==CHEM_ATOM_TYPE_HYDROGEN) { _l=0, _n++; }
        if (mol->a[neighbors->list[_i].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) { _l=1, _n++; }
        if (mol->a[neighbors->list[_i].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) { _l=2, _n++; }
        if (mol->a[neighbors->list[_i].list[3]]==CHEM_ATOM_TYPE_HYDROGEN) { _l=3, _n++; }
        if (_n==1) 
          {
          _n=4, mol->ff_i[mol->size_i].atom[0]=_i;
          while (_n--)
                 if (_n>_l) mol->ff_i[mol->size_i].atom[_n  ]=neighbors->list[_i].list[_n];
            else if (_n<_l) mol->ff_i[mol->size_i].atom[_n+1]=neighbors->list[_i].list[_n];
          mol->ff_i[mol->size_i].type=0, mol->ff_i[mol->size_i].k=GAFF_GENERIC_IMPROPER_K*4.184*180./PI, mol->ff_i[mol->size_i].v=35.3*PI/180.; //Pyramidal
          mol->size_i++;
          }
        }
    }
//Run over bonds
_t=mol->nedges;
while (_t--)
  {
  //Define bonds
  mol->ff_b[_t].atom[0]=_i=mol->edges[_t].vertice[0], mol->ff_b[_t].atom[1]=_j=mol->edges[_t].vertice[1];
  mol->ff_b[_t].k=mol->ff_b[_t].v=(double)NAN; 
  //Define dihedrals 
  if ( (neighbors->list[_i].size!=1)&&(neighbors->list[_j].size!=1) )
    {
         if (anchor_id[_i]!=anchor_id[_j])
           {
           if ( ((mol->ff_d[mol->size_d].atom[0]=_define_appropriate_torsion_atom(neighbors,_i,_j,mol->a))==(unsigned int)-1)||
                ((mol->ff_d[mol->size_d].atom[3]=_define_appropriate_torsion_atom(neighbors,_j,_i,mol->a))==(unsigned int)-1) )
             { 
             ylib_errno=YERROR_DATA_CONSISTMENT;
             free(mol->ff_b), mol->ff_b=0x0;
             free(mol->ff_g), mol->ff_g=0x0;
             free(mol->ff_i), mol->ff_i=0x0;
             free(mol->ff_d), mol->ff_d=0x0;
             return FALSE;
             }
           mol->ff_d[mol->size_d].atom[1]=_i, mol->ff_d[mol->size_d].atom[2]=_j, mol->ff_d[mol->size_d].d=neighbors->list[_i].size*neighbors->list[_j].size;
           _p=SIZE_DIH; while (_p--) { mol->ff_d[mol->size_d].k[_p]=mol->ff_d[mol->size_d].v[_p]=mol->ff_d[mol->size_d].n[_p]=(double)NAN; }
           mol->size_d++;
           }
    else if ( ( (inacycles->list[_i].size))&&( (inacycles->list[_j].size))&&( (overlap_list(FALSE,&inacycles->list[_i],FALSE,&inacycles->list[_j]))) )
           {//Type 1 torsions in cycles. Run over inacycles: screen all overlaps in aromatic cycles 
           _p=inacycles->list[_i].size;
           while (_p--)
             {
             _n=inacycles->list[_i].list[_p], _q=inacycles->list[_j].size;
             while (_q--)
               if (_n==inacycles->list[_j].list[_q])
                 {
                 _k=neighbors->list[_i].size;
                 while (_k--) 
                   if ( (_j!=neighbors->list[_i].list[_k])&&( (inacycles->list[neighbors->list[_i].list[_k]].size))&&(find_in_list(_n,&inacycles->list[neighbors->list[_i].list[_k]])!=(unsigned int)-1) )
                     {
                     _l=neighbors->list[_j].size;
                     while (_l--) 
                       if ( (_i!=neighbors->list[_j].list[_l])&&( (inacycles->list[neighbors->list[_j].list[_l]].size))&&(find_in_list(_n,&inacycles->list[neighbors->list[_j].list[_l]])!=(unsigned int)-1) )
                         {
                         mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[_k], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[_l];
                         mol->ff_i[mol->size_i].type=1, mol->ff_i[mol->size_i].k=(double)NAN, mol->ff_i[mol->size_i].v=0.; 
                         mol->size_i++;
                         }
                     }
                 }
             }
           }
    else if ( (mol->edges[_t].type==(int)'m')&&(neighbors->list[_i=mol->edges[_t].vertice[0]].size==3)&&(neighbors->list[_j=mol->edges[_t].vertice[1]].size==3) ) 
           { //Rule 2.0. Amide bond
           if (mol->a[_i]==CHEM_ATOM_TYPE_NITROGEN) { _i=mol->edges[_t].vertice[1], _j=mol->edges[_t].vertice[0]; }
           if ( (mol->a[_i]==CHEM_ATOM_TYPE_CARBON)&&(!(order[_j][3]))&&(!(order[_j][2]))&&(neighbors->list[*neighbors->list[_i].list].size==1)&&
              ( (mol->a[*neighbors->list[_i].list]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[*neighbors->list[_i].list]==CHEM_ATOM_TYPE_SULFUR) )       &&
              ( (mol->a[_j]==CHEM_ATOM_TYPE_NITROGEN)&&(!(order[_j][3]))&&(!(order[_j][2]))&&( (order[_j][1])) )                               &&
              ( (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)||(mol->a[neighbors->list[_j].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
             {
             //Rule 2.0.0. (O, S)=C-NH2
             if ( (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mol->a[neighbors->list[_j].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) )
               {//Quasy-dihedrals
               mol->ff_d[mol->size_d].atom[0]=*neighbors->list[_i].list, mol->ff_d[mol->size_d].atom[1]=_i, mol->ff_d[mol->size_d].atom[2]=_j, mol->ff_d[mol->size_d].atom[3]=neighbors->list[_j].list[1];
               mol->ff_d[mol->size_d].d=0., mol->ff_d[mol->size_d].k[0]=(double)NAN, mol->ff_d[mol->size_d].v[0]=PI, mol->ff_d[mol->size_d].n[0]=2.; 
               mol->size_d++;
               mol->ff_d[mol->size_d].atom[0]=*neighbors->list[_i].list, mol->ff_d[mol->size_d].atom[1]=_i, mol->ff_d[mol->size_d].atom[2]=_j, mol->ff_d[mol->size_d].atom[3]=neighbors->list[_j].list[2];
               mol->ff_d[mol->size_d].d=0., mol->ff_d[mol->size_d].k[0]=(double)NAN, mol->ff_d[mol->size_d].v[0]=PI, mol->ff_d[mol->size_d].n[0]=2.; 
               mol->size_d++;
               }
             else
               {
               //Rule 2.0.1. Inside a cycle
               _l=mol->cycles->size;
               while (mol->size_ar!=_l--)
                 if ( (mol->cycles->list[_l].size<8)&&(find_in_list(_i,&mol->cycles->list[_l])!=(unsigned int)-1)&&(find_in_list(_j,&mol->cycles->list[_l])!=(unsigned int)-1) )           
                   {//Rule 2.0.2. In a cycle
                   if (neighbors->list[_i].list[1]==_j)
                     {
                     if (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)
                       {//neighbors->list[_i].list[2]==R1, neighbors->list[_j].list[2]==R2
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N -H
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[2], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N -H
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N-R2
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[2], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N-R2
                       }
                     else
                       {//neighbors->list[_i].list[2]==R1, neighbors->list[_j].list[1]==R2
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N -H
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[2], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N -H
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N-R2
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[2], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N-R2
                       }
                     }
                   else
                     {
                     if (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)
                       {//neighbors->list[_i].list[1]==R1, neighbors->list[_j].list[2]==R2
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N -H
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[1], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N -H
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N-R2
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[1], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N-R2
                       }
                     else
                       {//neighbors->list[_i].list[1]==R1, neighbors->list[_j].list[1]==R2
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N -H
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[1], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N -H
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N-R2
                       mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[1], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N-R2
                       }
                     }
                   goto NEXT_BOND;
                   }
               //Rule 2.0.2. Out of a cycle
               if (neighbors->list[_i].list[1]==_j)
                 {
                 if (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)
                   {//neighbors->list[_i].list[2]==R1, neighbors->list[_j].list[2]==R2
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N -H
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[2], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N -H
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N-R2
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[2], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N-R2
                   }
                 else
                   {//neighbors->list[_i].list[2]==R1, neighbors->list[_j].list[1]==R2
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N -H
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[2], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N -H
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N-R2
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[2], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N-R2
                   }
                 }
               else
                 {
                 if (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)
                   {//neighbors->list[_i].list[1]==R1, neighbors->list[_j].list[2]==R2
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N -H
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[1], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N -H
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N-R2
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[1], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N-R2
                   }
                 else
                   {//neighbors->list[_i].list[1]==R1, neighbors->list[_j].list[1]==R2
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N -H
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[1], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[2], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N -H
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[0], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=0., mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //O -C-N-R2
                   mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[1], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[1], mol->ff_i[mol->size_i].type=4, mol->ff_i[mol->size_i].v=PI, mol->ff_i[mol->size_i].k=(double)NAN, mol->size_i++; //R1-C-N-R2
                   }
                 goto NEXT_BOND;
                 }
               }
             }
           else goto DEFINE_GENERIC_IMPROPER;
           }
    else if (mol->edges[_t].type!='1') //Rule 2.1. In anchor
           {//Improper
           DEFINE_GENERIC_IMPROPER: ; 
           if (neighbors->list[_i].size<neighbors->list[_j].size)
             { 
             //Get impropers count
             _n=0, _l=neighbors->list[_i].size; 
             while (_l--) 
               if ( (neighbors->list[_i].list[_l]!=_j)&&(find_in_list(neighbors->list[_i].list[_l],&neighbors->list[_j])==(unsigned int)-1) )
                 {
                 _k=neighbors->list[_j].size;
                 while (_k--)
                   if (neighbors->list[_j].list[_k]!=_i)
                     _n++;
                 }
             //Fill angles 
             _l=neighbors->list[_i].size; 
             while (_l--) 
               if ( (neighbors->list[_i].list[_l]!=_j)&&(find_in_list(neighbors->list[_i].list[_l],&neighbors->list[_j])==(unsigned int)-1) )
                 {
                 _k=neighbors->list[_j].size;
                 while (_k--)
                   if (neighbors->list[_j].list[_k]!=_i)
                     {
                     mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[_l], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[_k];
                     mol->ff_i[mol->size_i].type=_n, mol->ff_i[mol->size_i].k=mol->ff_i[mol->size_i].v=(double)NAN; 
                     mol->size_i++;
                     }
                 }
             }
           else
             { 
             //Get impropers count
             _n=0, _l=neighbors->list[_j].size; 
             while (_l--) 
               if ( (neighbors->list[_j].list[_l]!=_i)&&(find_in_list(neighbors->list[_j].list[_l],&neighbors->list[_i])==(unsigned int)-1) )
                 {
                 _k=neighbors->list[_i].size;
                 while (_k--)
                   if (neighbors->list[_i].list[_k]!=_j)
                     _n++;
                 }
             //Fill angles
             _l=neighbors->list[_j].size; 
             while (_l--) 
               if ( (neighbors->list[_j].list[_l]!=_i)&&(find_in_list(neighbors->list[_j].list[_l],&neighbors->list[_i])==(unsigned int)-1) )
                 {
                 _k=neighbors->list[_i].size;
                 while (_k--)
                   if (neighbors->list[_i].list[_k]!=_j)
                     {
                     mol->ff_i[mol->size_i].atom[0]=neighbors->list[_i].list[_k], mol->ff_i[mol->size_i].atom[1]=_i, mol->ff_i[mol->size_i].atom[2]=_j, mol->ff_i[mol->size_i].atom[3]=neighbors->list[_j].list[_l];
                     mol->ff_i[mol->size_i].type=_n, mol->ff_i[mol->size_i].k=mol->ff_i[mol->size_i].v=(double)NAN; 
                     mol->size_i++;
                     }
                 }
             }
           }
    NEXT_BOND: ;
    }
  }
//Stage 4. Compose excludes
if (!(mol->excl=(t_list*)malloc(sizeof(t_list)*mol->natoms+sizeof(unsigned int)*(mol->size_b+mol->size_g)))) { LABEL_MEMORY_ERROR_5: if (mol->size_d) free(mol->ff_d); mol->ff_d=0x0; goto LABEL_MEMORY_ERROR_4; }
else { _i=mol->natoms; while (_i--) mol->excl[_i].size=0; }
_i=mol->size_b; while (_i--) if (mol->ff_b[_i].atom[0]>mol->ff_b[_i].atom[1]) mol->excl[mol->ff_b[_i].atom[0]].size++; else mol->excl[mol->ff_b[_i].atom[1]].size++;
_i=mol->size_g; while (_i--) if (mol->ff_g[_i].atom[0]>mol->ff_g[_i].atom[2]) mol->excl[mol->ff_g[_i].atom[0]].size++; else mol->excl[mol->ff_g[_i].atom[2]].size++;
for (mol->excl[0].list=(unsigned int*)((void*)mol->excl+sizeof(t_list)*mol->natoms), _k=1; _k<mol->natoms; _k++) mol->excl[_k].list=mol->excl[_k-1].list+mol->excl[_k-1].size;
_t=FALSE; _k=mol->natoms; while (_k--) mol->excl[_k].size=0;
_k=mol->size_b;
while (_k--)
  {
  _i=mol->ff_b[_k].atom[0], _j=mol->ff_b[_k].atom[1];
  if (_i<_j) { if (find_in_list(_i,&mol->excl[_j])==(unsigned int)-1) mol->excl[_j].list[mol->excl[_j].size++]=_i; else _t=TRUE; }
  else       { if (find_in_list(_j,&mol->excl[_i])==(unsigned int)-1) mol->excl[_i].list[mol->excl[_i].size++]=_j; else _t=TRUE; }
  }
_k=mol->size_g;
while (_k--)
  {
  _i=mol->ff_g[_k].atom[0], _j=mol->ff_g[_k].atom[2];
  if (_i<_j) { if (find_in_list(_i,&mol->excl[_j])==(unsigned int)-1) mol->excl[_j].list[mol->excl[_j].size++]=_i; else _t=TRUE; }
  else       { if (find_in_list(_j,&mol->excl[_i])==(unsigned int)-1) mol->excl[_i].list[mol->excl[_i].size++]=_j; else _t=TRUE; }
  }
if ( (_t))
  {//Sync excludes if needed
  for (_i=_j=0; _i<mol->natoms; _j+=mol->excl[_i++].size)
    if (mol->excl[0].list+_j!=mol->excl[_i].list) 
      for (_k=0; _k<mol->excl[_i].size; _k++) //One by one copy because the memory may overlap
        *(mol->excl[0].list+_j+_k)=mol->excl[_i].list[_k];
  if (!(vp=realloc(mol->excl,sizeof(t_list)*mol->natoms+sizeof(unsigned int)*_j))) {                       free(mol->excl); mol->excl=0x0; goto LABEL_MEMORY_ERROR_5; }
  else for (mol->excl=(t_list*)vp, mol->excl[0].list=(unsigned int*)(vp+sizeof(t_list)*mol->natoms), _i=1; _i<mol->natoms; _i++) mol->excl[_i].list=mol->excl[_i-1].list+mol->excl[_i-1].size;
  }
_i=mol->natoms; while (_i--) if (!(u_qsort(mol->excl[_i].size,mol->excl[_i].list))) u_isort(mol->excl[_i].size,mol->excl[_i].list); //Sort the exclusions lists

//Exiting
if ( (_inacycles)) (*_inacycles)=inacycles; else free(inacycles);
return TRUE;
}


//GAFF parameterizators

//This service function compose gaff atom type from its chemical id and given enviroment suffix
unsigned int _compose_gaff_type(char gsuffix,unsigned int ytype,t_top *top)
{
unsigned int _i;
char gtype[sizeof(unsigned int)];

if (!(gsuffix))
  {//No enviromental index
  *(unsigned int*)gtype=top->ff_a[ytype].gtype; //Generic type
  _i=1; while (gtype[_i]!=' ') _i++;
  while(++_i!=sizeof(unsigned int)) gtype[_i]=' '; 
  }
else 
  {
  _i=1; while(++_i!=sizeof(unsigned int)) gtype[_i]=' '; //Special type
  switch(top->ff_a[ytype].chem_id)
    {
    case CHEM_ATOM_TYPE_CARBON   : { gtype[0]='c', gtype[1]=gsuffix; break; }
    case CHEM_ATOM_TYPE_NITROGEN : { gtype[0]='n'; if (gsuffix!='n') gtype[1]=gsuffix; break; }
    case CHEM_ATOM_TYPE_OXYGEN   : { gtype[0]='o', gtype[1]=gsuffix; break; }
    case CHEM_ATOM_TYPE_FLUORINE : { gtype[0]='f', gtype[1]=gsuffix; break; }
    case CHEM_ATOM_TYPE_SILICON  : { gtype[0]='p', gtype[1]=gsuffix; break; } //!!!!!!!!!!!!!!!!!!!!
    case CHEM_ATOM_TYPE_PHOSPHOR : { gtype[0]='p', gtype[1]=gsuffix; break; }
    case CHEM_ATOM_TYPE_SULFUR   : { gtype[0]='s', gtype[1]=gsuffix; break; } 
    case CHEM_ATOM_TYPE_CHLORINE : { gtype[0]='c', gtype[1]='l', gtype[2]=gsuffix; break; }
    case CHEM_ATOM_TYPE_BROMINE  : { gtype[0]='b', gtype[1]='r', gtype[2]=gsuffix; break; }
    case CHEM_ATOM_TYPE_IODINE   : { gtype[0]='i', gtype[1]=gsuffix; break; }
    default                      : { if (top->ff_a[ytype].chem_id<0) { gtype[0]='h'; gtype[1]=gsuffix; } //hydrogen
                                     else { ylib_errno=YERROR_INTERNAL_CODE; return (unsigned int)-1; } }
    }
  }
return *((unsigned int*)gtype); 
}

//This function defines parameters indexes from chemical type of the atom
inline unsigned int _get_gaff_top_index(unsigned char chem_type)
{
switch (chem_type)
  {
  case CHEM_ATOM_TYPE_HYDROGEN : return 0;
  case CHEM_ATOM_TYPE_CARBON   : return 1;
  case CHEM_ATOM_TYPE_NITROGEN : return 2;
  case CHEM_ATOM_TYPE_OXYGEN   : return 3;
  case CHEM_ATOM_TYPE_FLUORINE : return 4;
  case CHEM_ATOM_TYPE_SILICON  : return 8;
  case CHEM_ATOM_TYPE_PHOSPHOR : return 8;
  case CHEM_ATOM_TYPE_SULFUR   : return 9;
  case CHEM_ATOM_TYPE_CHLORINE : return 5;
  case CHEM_ATOM_TYPE_BROMINE  : return 6;
  case CHEM_ATOM_TYPE_IODINE   : return 7;
  default : { ylib_errno=YERROR_EXTERNAL_CODE; return (unsigned int)-1; }
  }
}
//This function does bonds compilation
//NOTE. ff_b atoms should contain ytypes of the bond
inline char parameterize_bonds_gaff(char generic,unsigned int type,t_ff_b *ff_b,char gtypes0,char gtypes1,t_mol *mol,t_top *top)
{
register unsigned int _i, _j;
unsigned int id;
t_ff_b _ff_b; 
//Try searching approach first
if ((_ff_b.atom[0]=_compose_gaff_type(gtypes0,mol->ytypes[ff_b->atom[0]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
if ((_ff_b.atom[1]=_compose_gaff_type(gtypes1,mol->ytypes[ff_b->atom[1]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
order_bond(&_ff_b.atom[0],&_ff_b.atom[1],_ff_b.atom[0],_ff_b.atom[1]);
if ( (find_in_sorted_lth_objects(&id,&_ff_b,top->size_b,(void*)top->ff_b,sizeof(t_ff_b),compare_bonds)))
  { ff_b->k=top->ff_b[id].k, ff_b->v=top->ff_b[id].v; return TRUE; }
else
  {//Try generic types
  if (generic)
    {
    if ((_ff_b.atom[0]=_compose_gaff_type(' ',top->ff_a[mol->ytypes[ff_b->atom[0]]].generic,top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
    if ((_ff_b.atom[1]=_compose_gaff_type(' ',top->ff_a[mol->ytypes[ff_b->atom[1]]].generic,top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
    order_bond(&_ff_b.atom[0],&_ff_b.atom[1],_ff_b.atom[0],_ff_b.atom[1]);
    if ( (find_in_sorted_lth_objects(&id,&_ff_b,top->size_b,(void*)top->ff_b,sizeof(t_ff_b),compare_bonds)))
      { ff_b->k=top->ff_b[id].k, ff_b->v=top->ff_b[id].v; }
    else //Try missing bonds
      {
      if (type=='1')
        {
        if ((_i=_get_gaff_top_index(mol->a[ff_b->atom[0]]))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
        if ((_j=_get_gaff_top_index(mol->a[ff_b->atom[1]]))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
        ff_b->k=top->bK[_i][_j], ff_b->v=top->br[_i][_j]; 
        }
      else 
        return FALSE;
      }
    }      
  return NTNF;
  }
LABEL_PARAMETERIZATION_FAILURE: ylib_errno=YERROR_EXTERNAL_CODE;
ff_b->k=0., ff_b->v=0.;
return FALSE;
}

//This function does angles compilation
inline char parameterize_angles_gaff(char generic,t_ff_g *ff_g,char gtypes0,char gtypes1,char gtypes2,t_mol *mol,t_top *top)
{
register unsigned int _l;
unsigned int id;
t_ff_g ff_gi, ff_gk;

//Try searching approach first
if ((ff_gi.atom[0]=_compose_gaff_type(gtypes0,mol->ytypes[ff_g->atom[0]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
if ((ff_gi.atom[1]=_compose_gaff_type(gtypes1,mol->ytypes[ff_g->atom[1]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
if ((ff_gi.atom[2]=_compose_gaff_type(gtypes2,mol->ytypes[ff_g->atom[2]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
order_angle(&ff_gi.atom[0],&ff_gi.atom[1],&ff_gi.atom[2],ff_gi.atom[0],ff_gi.atom[1],ff_gi.atom[2]);
if ( (find_in_sorted_lth_objects(&id,&ff_gi,top->size_g,(void*)top->ff_g,sizeof(t_ff_g),compare_angles)))
  { ff_g->v=top->ff_g[id].v, ff_g->k=top->ff_g[id].k; return TRUE; }
else 
  {
  if ( (generic))
    {//Try extrapolate with empirical rules 
    //Find ABA and CBC angles
    ff_gk.atom[0]=ff_gk.atom[2]=ff_gi.atom[2], ff_gi.atom[2]=ff_gi.atom[0], ff_gk.atom[1]=ff_gi.atom[1]; 
    order_angle(&ff_gi.atom[0],&ff_gi.atom[1],&ff_gi.atom[2],ff_gi.atom[0],ff_gi.atom[1],ff_gi.atom[2]);
    if ( (find_in_sorted_lth_objects(&id,&ff_gi,top->size_g,(void*)top->ff_g,sizeof(t_ff_g),compare_angles)))
      {
      ff_gi.v=top->ff_g[id].v, ff_gi.k=top->ff_g[id].k;
      }
    else
      {//Switch to generic type on failure
      ff_gi.v=0., ff_gi.k=0., _l=top->size_g; while (_l--) if (ff_gi.atom[1]==top->ff_g[_l].atom[1]) { ff_gi.k+=1., ff_gi.v+=(top->ff_g[_l].v-ff_gi.v)/ff_gi.k; }
      if (!ff_gi.k) goto LABEL_PARAMETERIZATION_FAILURE;
      }
    order_angle(&ff_gk.atom[0],&ff_gk.atom[1],&ff_gk.atom[2],ff_gk.atom[0],ff_gk.atom[1],ff_gk.atom[2]);
    if ( (find_in_sorted_lth_objects(&id,&ff_gk,top->size_g,(void*)top->ff_g,sizeof(t_ff_g),compare_angles)))
      {
      ff_gk.v=top->ff_g[id].v, ff_gi.k=top->ff_g[id].k;
      }
    else
      {//Switch to generic type
      ff_gk.v=0., ff_gk.k=0., _l=top->size_g; while (_l--) if (ff_gk.atom[1]==top->ff_g[_l].atom[1]) { ff_gk.k+=1., ff_gk.v+=(top->ff_g[_l].v-ff_gk.v)/ff_gk.k; }
      if (!ff_gk.k) goto LABEL_PARAMETERIZATION_FAILURE;
      }
    //Equlibrate
    ff_g->v=.5*(ff_gi.v+ff_gk.v), ff_g->k=.5*(ff_gi.k+ff_gk.k);
    }
  return NTNF;
  }
LABEL_PARAMETERIZATION_FAILURE: ylib_errno=YERROR_EXTERNAL_CODE;
ff_g->k=0., ff_g->v=0.;
return FALSE;
}

//This function does impropers compilation
//NOTE. ff_i atoms should contain ytypes of the angle
inline char parameterize_impropers_gaff(char generic,t_ff_i *ff_i,char gtypes0,char gtypes1,char gtypes2,char gtypes3,t_mol *mol,t_top *top)
{
unsigned int _i, _j, _k, _l, _id, _n, _s, _S=0;
t_ff_i _ff_i;

if ((_ff_i.atom[0]=_compose_gaff_type(gtypes0,mol->ytypes[ff_i->atom[0]],top))==(unsigned int)-1) { LABEL_PARAMETERIZATION_FAILURE: ylib_errno=YERROR_EXTERNAL_CODE; ff_i->k=0., ff_i->v=0.; return FALSE; }
if ((_ff_i.atom[1]=_compose_gaff_type(gtypes1,mol->ytypes[ff_i->atom[1]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
if ((_ff_i.atom[2]=_compose_gaff_type(gtypes2,mol->ytypes[ff_i->atom[2]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
if ((_ff_i.atom[3]=_compose_gaff_type(gtypes3,mol->ytypes[ff_i->atom[3]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
order_impr(ff_i->type,&_ff_i.atom[0],&_ff_i.atom[1],&_ff_i.atom[2],&_ff_i.atom[3],_ff_i.atom[0],_ff_i.atom[1],_ff_i.atom[2],_ff_i.atom[3]);
if (!(ff_i->type))
  {//type A1->(A0-A2-A3)
  _l=top->size_i;
  while (_l--)
    if ( (top->ff_i[_l].type))
      {
           if (top->ff_i[_l].atom[0]==_ff_i.atom[0])
             {
             _s=5;
             SCANN_IMPROPER_I: ;
             _n=0;
             _i=4; while ( (--_i)) if   (top->ff_i[_l].atom[_i]==_ff_i.atom[1])                       { _n++, _s++; break; }
             _j=4; while ( (--_j)) if ( (top->ff_i[_l].atom[_j]==_ff_i.atom[2])&&(_j!=_i) )           { _n++, _s++; break; }
             _k=4; while ( (--_k)) if ( (top->ff_i[_l].atom[_k]==_ff_i.atom[3])&&(_k!=_i)&&(_k!=_j) ) { _n++, _s++; break; }
             if (_n!=3) { _i=4; while ( (--_i)) if (top->ff_i[_l].atom[_i]==*((unsigned int*)"x   ")) _n++; } 
             if ( (_n==3)&&(_s>_S) ) { _id=_l; if ((_S=_s)==8) break; }
             }
      else if ( (_S<4)&&(top->ff_i[_l].atom[0]==*((unsigned int*)"x   ")) )
             { _s=1; goto SCANN_IMPROPER_I; }
      }
  } 
//Use generic impropers for chains - they are sorted to be type==0
//else
//  {//type A0-A1-A2-A3. at parameterization stage it is determined up to +/-PI
//  _l=top->size_i;
//  while (_l--)
//    {
//        if (top->ff_i[_l].atom[1]==_ff_i.atom[1])
//          {
//                if (top->ff_i[_l].atom[2]==_ff_i.atom[2])
//                  { 
//                  _s=7;
//                  SCANN_IMPROPER_0: ;
//                  if   (top->ff_i[_l].atom[0]==_ff_i.atom[0]) _s++; 
//                  else { if (top->ff_i[_l].atom[0]!=*((unsigned int*)"x   ")) continue; }  
//                  if   (top->ff_i[_l].atom[3]==_ff_i.atom[3]) _s++;
//                  else { if (top->ff_i[_l].atom[3]!=*((unsigned int*)"x   ")) continue; }  
//                  if (_s>_S) { _id=_l; if ((_S=_s)==9) break; }
//                  }
//           else if ((top->ff_i[_l].atom[2]==*((unsigned int*)"x   "))&&(_S<6))
//                  { _s=4; goto SCANN_IMPROPER_0; }
//           }
//    else if (top->ff_i[_l].atom[1]==*((unsigned int*)"x   "))
//           {
//                if ((top->ff_i[_l].atom[2]==_ff_i.atom[2])&&(_S<6))
//                  { _s=4; goto SCANN_IMPROPER_I; }
//           else if ((top->ff_i[_l].atom[2]==*((unsigned int*)"x   "))&&(_S<3))
//                  { _s=1; goto SCANN_IMPROPER_0; }
//           }
//    }
//  }
//Assign the parameters including previous partial setup in compile_mol()
if ( (_S)) 
  {
  ff_i->k=(!(ff_i->type)) ? top->ff_i[_id].k : top->ff_i[_id].k/(double)ff_i->type; 
  ff_i->v=(!(isnan(ff_i->v))) ? ff_i->v : top->ff_i[_id].v;
  return TRUE;
  }
else
  {
  if (generic)
    {
    ff_i->k=(!(ff_i->type)) ? 4.184*GAFF_GENERIC_IMPROPER_K*180./PI : 4.184*GAFF_GENERIC_IMPROPER_K*180./PI/2.; 
    ff_i->v=(!(isnan(ff_i->v))) ? ff_i->v : 0.;
    }
  return NTNF;
  }
}

//This function does torsions compilation
//NOTE. ff_d atoms should contain ytypes of the angle
inline char parameterize_dihedrals_gaff(char generic,t_ff_d *ff_d,char gtypes0,char gtypes1,char gtypes2,char gtypes3,t_mol *mol,t_top *top)
{
unsigned int _l, _id, _s, _S=0;
t_ff_d _ff_d;
t_ff_i _ff_i;
//Parameterize quasy-dihedrals
if (!(ff_d->d))
  {
  _ff_i.atom[0]=ff_d->atom[0], _ff_i.atom[1]=ff_d->atom[1], _ff_i.atom[2]=ff_d->atom[2], _ff_i.atom[3]=ff_d->atom[3], _ff_i.type=0, _ff_i.k=ff_d->k[0], _ff_i.v=ff_d->v[0];
  if (!(parameterize_impropers_gaff(generic,&_ff_i,gtypes0,gtypes1,gtypes2,gtypes3,mol,top))) return FALSE;
  ff_d->k[0]=_ff_i.k, ff_d->v[0]=_ff_i.v, ff_d->n[0]=2., _l=SIZE_DIH; while (--_l) { ff_d->k[_l]=ff_d->v[_l]=0., ff_d->n[_l]=2.; }
  return TRUE;
  }
//Parameterize "normal" dihedrals
if ((_ff_d.atom[0]=_compose_gaff_type(gtypes0,mol->ytypes[ff_d->atom[0]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
if ((_ff_d.atom[1]=_compose_gaff_type(gtypes1,mol->ytypes[ff_d->atom[1]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
if ((_ff_d.atom[2]=_compose_gaff_type(gtypes2,mol->ytypes[ff_d->atom[2]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
if ((_ff_d.atom[3]=_compose_gaff_type(gtypes3,mol->ytypes[ff_d->atom[3]],top))==(unsigned int)-1) goto LABEL_PARAMETERIZATION_FAILURE;
order_dih(&_ff_d.atom[0],&_ff_d.atom[1],&_ff_d.atom[2],&_ff_d.atom[3],_ff_d.atom[0],_ff_d.atom[1],_ff_d.atom[2],_ff_d.atom[3]);
_l=top->size_d;
while (_l--)
  {
       if (top->ff_d[_l].atom[1]==_ff_d.atom[1])
         {
              if (top->ff_d[_l].atom[2]==_ff_d.atom[2])
                { 
                _s=7;
                SCANN_DIHEDRAL: ;
                if   (top->ff_d[_l].atom[0]==_ff_d.atom[0]) _s++; 
                else { if (top->ff_d[_l].atom[0]!=*((unsigned int*)"x   ")) continue; }  
                if   (top->ff_d[_l].atom[3]==_ff_d.atom[3]) _s++;
                else { if (top->ff_d[_l].atom[3]!=*((unsigned int*)"x   ")) continue; }  
                if (_s>_S) { _id=_l; if ((_S=_s)==9) break; }
                }
         else if ((top->ff_d[_l].atom[2]==*((unsigned int*)"x   "))&&(_S<6))
                { _s=4; goto SCANN_DIHEDRAL; }
         }
  else if (top->ff_d[_l].atom[1]==*((unsigned int*)"x   "))
         {
              if ((top->ff_d[_l].atom[2]==_ff_d.atom[2])&&(_S<6))
                { _s=4; goto SCANN_DIHEDRAL; }
         else if ((top->ff_d[_l].atom[2]==*((unsigned int*)"x   "))&&(_S<3))
                { _s=1; goto SCANN_DIHEDRAL; }
         }
  }
//Assign the parameters including the angles amount from compose_mol()
if ( (_S))
  {
  _l=SIZE_DIH; while (_l--) { ff_d->k[_l]=top->ff_d[_id].k[_l], ff_d->v[_l]=top->ff_d[_id].v[_l], ff_d->n[_l]=top->ff_d[_id].n[_l]; }
  return TRUE;
  }
else 
  {  
  if (generic) { _l=SIZE_DIH; while (_l--) { ff_d->k[_l]=ff_d->v[_l]=0., ff_d->n[_l]=2.; } }
  return NTNF;
  }
LABEL_PARAMETERIZATION_FAILURE: ylib_errno=YERROR_EXTERNAL_CODE;
_l=SIZE_DIH; while (_l--) { ff_d->k[_l]=ff_d->v[_l]=0., ff_d->n[_l]=2.; }
return FALSE;
}


//This function compile gaff part of the FF model
//Note. It demands that aromatic cycles can't contain X_=_Y or Y=X=Z atoms
//Note. It attempts to avoid 5-ring paradox so is compiling the molecule in fragments: cycles -> intercycles -> resonance chanins -> rest
unsigned int parameterize_mol_GAFF(char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_clist *inacycles,t_mol *mol,t_top *top)
{
register unsigned int _i, _j, _k, _l, _s, _t;
char _c, _c0, _c1, _c2, _c3;
unsigned int *temp, count=0;
char *gtypes;

//Stage 0. Preparations.
//Setup GAFF types
if (!(gtypes=(char*)calloc(mol->natoms,sizeof(char)))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return (unsigned int)-1; }
if (!(temp=(unsigned int*)malloc(sizeof(unsigned int)*mol->natoms))) { free(gtypes); gtypes=0x0; goto LABEL_MEMORY_ERROR; }
//Stage I. Enviromental assignance #1. Onevalent (tio-)ketones and hydrogens.
_i=mol->natoms;
while (_i--)
  //Stage I.1. Insert (tio-)ketones
       if (mol->a[_i]==CHEM_ATOM_TYPE_CARBON)
         {
         if ( (neighbors->list[_i].size==3)&&(order[_i][2]==1)&&(neighbors->list[_j=*neighbors->list[_i].list].size==1) )
           {
                if (mol->a[_j]==CHEM_ATOM_TYPE_OXYGEN)
                  {
                  gtypes[_i]=' ', gtypes[_j]=' '; // 'c   ' && 'o   '
                  _k=neighbors->list[_i].list[1], _l=neighbors->list[_i].list[2];
                  if (neighbors->list[_k].size==1) 
                    { 
                         if (mol->a[_k]==CHEM_ATOM_TYPE_OXYGEN) gtypes[_k]=' ';
                    else if (mol->a[_k]==CHEM_ATOM_TYPE_SULFUR) gtypes[_k]='2';
                    } 
                  if (neighbors->list[_l].size==1) 
                    { 
                         if (mol->a[_l]==CHEM_ATOM_TYPE_OXYGEN) gtypes[_l]=' ';
                    else if (mol->a[_l]==CHEM_ATOM_TYPE_SULFUR) gtypes[_l]='2';
                    } 
                  }
           else if (mol->a[_j]==CHEM_ATOM_TYPE_SULFUR) 
                  {
                  gtypes[_i]=' ', gtypes[_j]='2';  // 'c   ' && 's2  '
                  _k=neighbors->list[_i].list[1], _l=neighbors->list[_i].list[2];
                  if (neighbors->list[_k].size==1) 
                    { 
                         if (mol->a[_k]==CHEM_ATOM_TYPE_OXYGEN) gtypes[_k]=' ';
                    else if (mol->a[_k]==CHEM_ATOM_TYPE_SULFUR) gtypes[_k]='2';
                    } 
                  if (neighbors->list[_l].size==1) 
                    { 
                         if (mol->a[_l]==CHEM_ATOM_TYPE_OXYGEN) gtypes[_l]=' ';
                    else if (mol->a[_l]==CHEM_ATOM_TYPE_SULFUR) gtypes[_l]='2';
                    } 
                  }
           }
         }
  //Stage I.2. Insert NO2
  else if (mol->a[_i]==CHEM_ATOM_TYPE_NITROGEN)
         {
         if ( (neighbors->list[_i].size==3)&&(order[_i][2]==2)&&(mol->a[_j=neighbors->list[_i].list[0]]==CHEM_ATOM_TYPE_OXYGEN)&&(mol->a[_k=neighbors->list[_i].list[1]]==CHEM_ATOM_TYPE_OXYGEN)&&(neighbors->list[_j].size==1)&&(neighbors->list[_k].size==1) )
           { gtypes[_i]='o', gtypes[_j]=' ', gtypes[_k]=' '; }
         }
  //Stage I.3. Define hydrogens
  else if (mol->a[_i]==CHEM_ATOM_TYPE_HYDROGEN)
         {
         if (neighbors->list[_i].size==1)  
           switch (mol->a[*neighbors->list[_i].list])
             {
             case CHEM_ATOM_TYPE_CARBON   : { if ( ( (order[*neighbors->list[_i].list][3]))||( (order[*neighbors->list[_i].list][2])) ) gtypes[_i]='a'; break; }
             case CHEM_ATOM_TYPE_NITROGEN : { gtypes[_i]='n'; break; }
             case CHEM_ATOM_TYPE_OXYGEN   : { gtypes[_i]='o'; break; }
             case CHEM_ATOM_TYPE_FLUORINE : { gtypes[_i]='o'; break; }
             case CHEM_ATOM_TYPE_SULFUR   : { gtypes[_i]='s'; break; }
             case CHEM_ATOM_TYPE_CHLORINE : { gtypes[_i]='s'; break; }
             case CHEM_ATOM_TYPE_PHOSPHOR : { gtypes[_i]='p'; break; }
             //NOTE. Don't define hc because it interfere with cycles 'c'/'d' flipping (do it later at the second enviromental assistance).
             }
         }

//Stage II.1. Compile cycles cycle-by-cycle to avoid theoretical 'c'/'d' flipping loop.
_k=mol->size_ar; 
while (_k--)
  {//Work on _k-th cycle
  if ((_l=mol->cycles->list[_k].size)==6) 
    {
    while (_l--)
      { //If all double bonds are inside the cycle then mark it as 'a'
      _i=mol->cycles->list[_k].list[_l], _j=*neighbors->list[_i].list;
      if ( ( (order[_i][3]))||(order[_i][2]!=1)||(find_in_list(_j,&mol->cycles->list[_k])==(unsigned int)-1) ) goto MARK_KTH_CYCLE_C;
      switch (mol->a[_i])
        {
        case CHEM_ATOM_TYPE_CARBON   : { if (neighbors->list[_j].size==3) break; else goto MARK_KTH_CYCLE_C; } 
        case CHEM_ATOM_TYPE_NITROGEN : 
        case CHEM_ATOM_TYPE_PHOSPHOR : { if (neighbors->list[_j].size==2) break; else goto MARK_KTH_CYCLE_C; }
        default                      : goto MARK_KTH_CYCLE_C;
        }
      }
    _l=6; while (_l--) gtypes[mol->cycles->list[_k].list[_l]]='a'; //It's aromatic
    }
  else
    { 
    MARK_KTH_CYCLE_C:     _l=mol->cycles->list[_k].size;
    while (_l--) 
      {
      _i=mol->cycles->list[_k].list[_l];
      if ( (!(gtypes[_i]))&&( (mol->a[_i]==CHEM_ATOM_TYPE_CARBON)||(mol->a[_i]==CHEM_ATOM_TYPE_NITROGEN)||(mol->a[_i]==CHEM_ATOM_TYPE_PHOSPHOR) ) )
        gtypes[_i]='c';
      }
    }
  }
//Stage II.2. Compile in-cycles (switching c/d on a fly)
_t=mol->size_b;
while (_t--)
  if ( ( (isnan(mol->ff_b[_t].k)))&&( (gtypes[_i=mol->ff_b[_t].atom[0]]))&&( (gtypes[_j=mol->ff_b[_t].atom[1]])) )
    {
    if ( (gtypes[_i]=='c')&&(gtypes[_j]=='c') )
      { 
      _s=0; while ( ((mol->edges[_s].vertice[0]!=_i)||(mol->edges[_s].vertice[1]!=_j))&&((mol->edges[_s].vertice[1]!=_i)||(mol->edges[_s].vertice[0]!=_j)) ) _s++;
      _c0='c', _c1=(mol->edges[_s].type!='2') ? 'c' : 'd';  
      }
    else { _c0=gtypes[_i], _c1=gtypes[_j]; }
    if (!(_c=parameterize_bonds_gaff(TRUE,mol->edges[_t].type,&mol->ff_b[_t],_c0,_c1,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE; 
    else if (_c==NTNF) count++;
    }
_t=mol->size_g;
while (_t--)
  if ( ( (isnan(mol->ff_g[_t].k)))&&( (gtypes[_i=mol->ff_g[_t].atom[0]]))&&( (gtypes[_j=mol->ff_g[_t].atom[1]]))&&( (gtypes[_k=mol->ff_g[_t].atom[2]])) )
    {
    _c1=gtypes[_j];
    if (gtypes[_j]=='c')
      {
      if (gtypes[_i]=='c')
        {
        _s=0; while ( ((mol->edges[_s].vertice[0]!=_i)||(mol->edges[_s].vertice[1]!=_j))&&((mol->edges[_s].vertice[1]!=_i)||(mol->edges[_s].vertice[0]!=_j)) ) _s++;
        if (mol->edges[_s].type!='2') _c0='c'; else _c0='d';
        }
      else _c0=gtypes[_i];
      if (gtypes[_k]=='c')
        {
        _s=0; while ( ((mol->edges[_s].vertice[0]!=_k)||(mol->edges[_s].vertice[1]!=_j))&&((mol->edges[_s].vertice[1]!=_k)||(mol->edges[_s].vertice[0]!=_j)) ) _s++;
        if (mol->edges[_s].type!='2') _c2='c'; else _c2='d';
        }
      else _c2=gtypes[_k];
      }
    else
      { _c0=gtypes[_i], _c2=gtypes[_k]; }
    if (!(_c=parameterize_angles_gaff(TRUE,&mol->ff_g[_t],_c0,_c1,_c2,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
    else if (_c==NTNF) count++; 
    }
_t=mol->size_i;
while (_t--)
  if ( ( (isnan(mol->ff_i[_t].k)))&&( (gtypes[_i=mol->ff_i[_t].atom[0]]))&&( (gtypes[_j=mol->ff_i[_t].atom[1]]))&&( (gtypes[_k=mol->ff_i[_t].atom[2]]))&&( (gtypes[_l=mol->ff_i[_t].atom[3]])) )
    {
    if (!(mol->ff_i[_t].type))
      {// A->(B, C, D)
      _c0=(gtypes[_i]=='c') ? 'a' : gtypes[_i], _c1=(gtypes[_j]=='c') ? 'a' : gtypes[_j], _c2=(gtypes[_k]=='c') ? 'a' : gtypes[_k], _c3=(gtypes[_l]=='c') ? 'a' : gtypes[_l];  
      if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_t],_c0,_c1,_c2,_c3,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
      else if (_c==NTNF) count++; 
      }
    else
      {// A-B-C-D anchor[B]==anchor[C], otherwise impropers are pointless
      if (gtypes[_j]=='c')
        {
        if (gtypes[_k]=='c') 
          {//?~'c'~'c'~?
          _s=0; while ( ((mol->edges[_s].vertice[0]!=_k)||(mol->edges[_s].vertice[1]!=_j))&&((mol->edges[_s].vertice[1]!=_k)||(mol->edges[_s].vertice[0]!=_j)) ) _s++;
          if (mol->edges[_s].type=='2')
            {//(?) 'c'~'c'='d'~'d' (?)
            _c0=(gtypes[_i]=='c') ? 'c' : gtypes[_i], _c1='c', _c2='d', _c3=(gtypes[_l]=='c') ? 'd' : gtypes[_l];
            }
          else
            {// ?~'c'-'c'~?
            _c1='c', _c2='c';
            if (gtypes[_i]=='c')
              {//(?) 'd'='c'~?~?
              _s=0; while ( ((mol->edges[_s].vertice[0]!=_i)||(mol->edges[_s].vertice[1]!=_j))&&((mol->edges[_s].vertice[1]!=_i)||(mol->edges[_s].vertice[0]!=_j)) ) _s++;
              _c0=(mol->edges[_s].type=='2') ? 'd': 'c';
              }
            else _c0=gtypes[_i];
            if (gtypes[_l]=='c')
              {
              _s=0; while ( ((mol->edges[_s].vertice[0]!=_l)||(mol->edges[_s].vertice[1]!=_k))&&((mol->edges[_s].vertice[1]!=_l)||(mol->edges[_s].vertice[0]!=_k)) ) _s++;
              _c3=(mol->edges[_s].type=='2') ? 'd' : 'c';
              }
            else _c3=gtypes[_l];
            }
          }
        else
          {//?~'c'~?~?
          _c1='c', _c2=gtypes[_k], _c3=gtypes[_l];
          if (gtypes[_i]=='c')
            {//(?) 'd'='c'~?~? 
            _s=0; while ( ((mol->edges[_s].vertice[0]!=_i)||(mol->edges[_s].vertice[1]!=_j))&&((mol->edges[_s].vertice[1]!=_i)||(mol->edges[_s].vertice[0]!=_j)) ) _s++;
            _c1=(mol->edges[_s].type=='2') ? 'd' : 'c';
            }
          else _c0=gtypes[_i];
          }
        }
      else
        {
        if (gtypes[_k]=='c')
          {//?~?~'c'~?
          _c0=gtypes[_i], _c1=gtypes[_j], _c2='c';
          if (gtypes[_l]=='c')
            {//?~?~'c'='d' (?)
            _s=0; while ( ((mol->edges[_s].vertice[0]!=_l)||(mol->edges[_s].vertice[1]!=_k))&&((mol->edges[_s].vertice[1]!=_l)||(mol->edges[_s].vertice[0]!=_k)) ) _s++;
            _c3=(mol->edges[_s].type=='2') ? 'd' : 'c';
            }
          else _c3=gtypes[_l];
          }
        else { _c0=gtypes[_i], _c1=gtypes[_j], _c2=gtypes[_k], _c3=gtypes[_l]; }
        }
      if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_t],_c0,_c1,_c2,_c3,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
      else if (_c==NTNF) count++; 
      }
   }
//Stage II.3. Set all into 'a' for "external" compilation
_i=mol->natoms; while (_i--) if (gtypes[_i]=='c') gtypes[_i]='a';

//Stage II.3. Re-compile cross-cycles parameters
_t=mol->size_b;
while (_t--) 
  if ( (mol->edges[_t].type=='1')&&(!(isnan(mol->ff_b[_t].k)))&&( (inacycles->list[_i=mol->ff_b[_t].atom[0]].size))&&( (inacycles->list[_j=mol->ff_b[_t].atom[1]].size)) )
    {
    if (anchor_id[_i]!=anchor_id[_j])
      {
      if (!(_c=parameterize_bonds_gaff(FALSE,mol->edges[_t].type,&mol->ff_b[_t],'p','p',mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
      else if (_c==NTNF) count++; 
      }
    }
_t=mol->size_g;
while (_t--)
  if ( (!(isnan(mol->ff_g[_t].k)))&&( (inacycles->list[_i=mol->ff_g[_t].atom[0]].size))&&( (inacycles->list[_j=mol->ff_g[_t].atom[1]].size))&&( (inacycles->list[_k=mol->ff_g[_t].atom[2]].size)) )
    {
         if (anchor_id[_i]!=anchor_id[_j]) 
           {
           if ( (gtypes[_k]))
             {
             if (!(_c=parameterize_angles_gaff(FALSE,&mol->ff_g[_t],'p','p',gtypes[_k],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
             else if (_c==NTNF) count++;
             }
           }
    else if (anchor_id[_j]!=anchor_id[_k])
           {
           if ( (gtypes[_i]))
             {
             if (!(_c=parameterize_angles_gaff(FALSE,&mol->ff_g[_t],gtypes[_i],'p','p',mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
             else if (_c==NTNF) count++;
             }
           }
    else   {
           if ( ( (gtypes[_i]))&&( (gtypes[_k])) )
             {
             _l=neighbors->list[_j].size;
             while (_l--)
               if ( ( (inacycles->list[neighbors->list[_j].list[_l]].size))&&(anchor_id[_j]!=anchor_id[neighbors->list[_j].list[_l]]) )
                 {
                 if (!(_c=parameterize_angles_gaff(FALSE,&mol->ff_g[_t],gtypes[_i],'p',gtypes[_k],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                 else if (_c==NTNF) count++;
                 break;
                 }
             }
           }
    }
_t=mol->size_i;
while (_t--)
  {
  _i=mol->ff_i[_t].atom[0], _j=mol->ff_i[_t].atom[1], _k=mol->ff_i[_t].atom[2], _l=mol->ff_i[_t].atom[3];
  if (!(mol->ff_i[_t].type))
    {//A->(B-C-D) angle
         if ( (anchor_id[_i]!=anchor_id[_j])&&(!(isnan(mol->ff_i[_k].k)))&&( (inacycles->list[_i].size))&&( (inacycles->list[_j].size))&&( (gtypes[_k]))&&( (gtypes[_l])) )
           {
           if (!(_c=parameterize_impropers_gaff(FALSE,&mol->ff_i[_t],'p','p',gtypes[_k],gtypes[_l],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
           else if (_c==NTNF) count++;
           }
    else if ( (anchor_id[_i]!=anchor_id[_k])&&(!(isnan(mol->ff_i[_t].k)))&&( (inacycles->list[_i].size))&&( (gtypes[_j]))&&( (inacycles->list[_k].size))&&( (gtypes[_l])) )
           {
           if (!(_c=parameterize_impropers_gaff(FALSE,&mol->ff_i[_t],'p',gtypes[_j],'p',gtypes[_l],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
           else if (_c==NTNF) count++;  
           }
    else if ( (anchor_id[_i]!=anchor_id[_l])&&(!(isnan(mol->ff_i[_t].k)))&&( (inacycles->list[_i].size))&&( (gtypes[_j]))&&( (gtypes[_k]))&&( (inacycles->list[_l].size)) )
           {
           if (!(_c=parameterize_impropers_gaff(FALSE,&mol->ff_i[_l],'p',gtypes[_j],gtypes[_k],'p',mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
           else if (_c==NTNF) count++;  
           }
    }
  else
    { //A-B-C-D angle
         if ( (anchor_id[_i]!=anchor_id[_j])&&(!(isnan(mol->ff_i[_t].k)))&&( (inacycles->list[_i].size))&&( (inacycles->list[_j].size))&&( (gtypes[_k]))&&( (gtypes[_l])) )
           {
           if ( (anchor_id[_k]!=anchor_id[_l])&&(!(isnan(mol->ff_i[_t].k)))&&( (inacycles->list[_k].size))&&( (inacycles->list[_l].size)) )
             {
             if (!(_c=parameterize_impropers_gaff(FALSE,&mol->ff_i[_t],'p','p','q','q',mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
             else if (_c==NTNF) count++;
             }
           else
             {
             if (!(_c=parameterize_impropers_gaff(FALSE,&mol->ff_i[_t],'p','p',gtypes[_k],gtypes[_l],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
             else if (_c==NTNF) count++;
             }
           }
    else if ( (anchor_id[_k]!=anchor_id[_l])&&(!(isnan(mol->ff_i[_t].k)))&&( (gtypes[_i]))&&( (gtypes[_j]))&&( (inacycles->list[_k].size))&&( (inacycles->list[_l].size)) )
           {
           if (!(_c=parameterize_impropers_gaff(FALSE,&mol->ff_i[_t],gtypes[_i],gtypes[_j],'p','p',mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
           else if (_c==NTNF) count++;
           }
    }
  }
//Parameterize inter-cycle dihedrals
_t=mol->size_d;
while (_t--)
  {
  _i=mol->ff_d[_t].atom[0], _j=mol->ff_d[_t].atom[1], _k=mol->ff_d[_t].atom[2], _l=mol->ff_d[_t].atom[3];
  if ( ( (gtypes[_i]))&&( (inacycles->list[_j].size))&&( (inacycles->list[_k].size))&&( (gtypes[_l])) )
    {
         if ( (inacycles->list[_i].size))
           {
           if (gtypes[_i]=='s') _c0='s'; else _c0='a';
           if ( (inacycles->list[_l].size))
             {
             if (gtypes[_l]=='s') _c3='s'; else _c3='a';
             if (!(_c=parameterize_dihedrals_gaff(TRUE,&mol->ff_d[_t],_c0,'p','p',_c3,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE; 
             else if (_c==NTNF) count++; 
             }
           else 
             {
             if (!(_c=parameterize_dihedrals_gaff(TRUE,&mol->ff_d[_t],_c0,'p','p',gtypes[_l],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE; 
             else if (_c==NTNF) count++; 
             }
           }  
    else if ( (inacycles->list[_l].size))
           {
           if (gtypes[_l]=='s') _c3='s'; else _c3='a';
           if (!(_c=parameterize_dihedrals_gaff(TRUE,&mol->ff_d[_t],gtypes[_i],'p','p',_c3,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE; 
           else if (_c==NTNF) count++; 
           }
    }
  }

//Stage III. Enviromental assignance #2.
//Stage III.1. Insert (sulf-)amides N 
_i=mol->natoms;
while (_i--)
  if ( (!gtypes[_i])&&(mol->a[_i]==CHEM_ATOM_TYPE_NITROGEN)&&(!order[_i][3])&&(!order[_i][2])&&(neighbors->list[_i].size==3) )
    {
         if ( (!order[_j=neighbors->list[_i].list[0]][3])&&(order[_j][2])&&
            ( (mol->a[_j]==CHEM_ATOM_TYPE_CARBON)||(mol->a[_j]==CHEM_ATOM_TYPE_SULFUR)||(mol->a[_j]==CHEM_ATOM_TYPE_PHOSPHOR) )                  &&
            (neighbors->list[_k=*neighbors->list[_j].list].size==1)&&( (mol->a[_k]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[_k]==CHEM_ATOM_TYPE_SULFUR) ) )
           gtypes[_i]=' '; //amide N 
    else if ( (!order[_j=neighbors->list[_i].list[1]][3])&&(order[_j][2])&&
            ( (mol->a[_j]==CHEM_ATOM_TYPE_CARBON)||(mol->a[_j]==CHEM_ATOM_TYPE_SULFUR)||(mol->a[_j]==CHEM_ATOM_TYPE_PHOSPHOR) )                  &&
            (neighbors->list[_k=*neighbors->list[_j].list].size==1)&&( (mol->a[_k]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[_k]==CHEM_ATOM_TYPE_SULFUR) ) )
           gtypes[_i]=' '; //amide N
    else if ( (!order[_j=neighbors->list[_i].list[2]][3])&&(order[_j][2])&&
            ( (mol->a[_j]==CHEM_ATOM_TYPE_CARBON)||(mol->a[_j]==CHEM_ATOM_TYPE_SULFUR)||(mol->a[_j]==CHEM_ATOM_TYPE_PHOSPHOR) )                  &&
            (neighbors->list[_k=*neighbors->list[_j].list].size==1)&&( (mol->a[_k]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[_k]==CHEM_ATOM_TYPE_SULFUR) ) )
           gtypes[_i]=' '; //amide N 

    }
//Stage III.2. N-H on aromatic
_i=mol->natoms;
while (_i--)
  if ( (!gtypes[_i])&&(mol->a[_i]==CHEM_ATOM_TYPE_NITROGEN)&&(!order[_i][3])&&(!order[_i][2]) )
    {
    if ( (order[_i][1])) gtypes[_i]=' ';
    else 
      {
      if ((_j=neighbors->list[_i].size)==3) gtypes[_i]='3'; else if (_j==4) gtypes[_i]='4'; 
      while (_j--)
        if (mol->a[neighbors->list[_i].list[_j]]==CHEM_ATOM_TYPE_HYDROGEN)
          {
          _j=neighbors->list[_i].size; while (_j--) if ( (inacycles->list[neighbors->list[_i].list[_j]].size)) { gtypes[_i]='h'; break; }
          break;
          }
      }
    }
//Stage III.3. H-bonded to aliphatic carbons
_i=mol->natoms;
while (_i--)
  if ( (mol->a[_i]==CHEM_ATOM_TYPE_HYDROGEN)&&(!gtypes[_i])&&(neighbors->list[_i].size==1)&&(mol->a[*neighbors->list[_i].list]==CHEM_ATOM_TYPE_CARBON)&&(!(inacycles->list[*neighbors->list[_i].list].size)) )
    gtypes[_i]='c'; //hc can't be define earlier because it would interfere with cycles 'c'/'d' flipping
//Stage III.4. H-bonded atoms OH and SH
_i=mol->natoms;
while (_i--)
  if ( (mol->a[_i]==CHEM_ATOM_TYPE_HYDROGEN)&&(neighbors->list[_i].size==1)&&(!(gtypes[_j=*neighbors->list[_i].list]))&&( (mol->a[_j]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[_j]==CHEM_ATOM_TYPE_SULFUR) ) ) 
    gtypes[_j]='h';
//Stage III.5. Compile alifatic cycles
for (_i=mol->size_ar; _i<mol->cycles->size; _i++) 
  {
  _j=mol->cycles->list[_i].size;
       if (_j==3)
         {
         while (_j--)
           if ( (!(gtypes[_l=mol->cycles->list[_i].list[_j]]))||(gtypes[_l]=='u')||(gtypes[_l]=='x')||(gtypes[_l]=='v')||(gtypes[_l]=='y') )
             {
             if ( ( (order[_l][1]))||( (order[_l][2]))||( (order[_l][3])) ) gtypes[_l]='u';
             else                                                           gtypes[_l]='x';
             }
         }
  else if (_j==4)
         {
         while (_j--)
           if ( (!(gtypes[_l=mol->cycles->list[_i].list[_j]]))||(gtypes[_l]=='v')||(gtypes[_l]=='y') )
             {
             if ( ( (order[_l][1]))||( (order[_l][2]))||( (order[_l][3])) ) gtypes[_l]='v';
             else                                                           gtypes[_l]='y';
             }
          }
  }
//Stage III.6. Edit hypervalent P and S near resonance
_i=mol->natoms;
while (_i--)
       if (mol->a[_i]==CHEM_ATOM_TYPE_PHOSPHOR)
         {
              if (neighbors->list[_i].size==3)
                {
                if ( (!(order[_i][3]))&&(order[_i][2]==1)&&( (!(gtypes[_i]))||(gtypes[_i]=='4') ) )
                  {//Look for P(4)-X=Y  
                  _j=neighbors->list[_i].size;
                  while (_j--)
                    if ( (*neighbors->list[_k=neighbors->list[_i].list[_j]].list!=_i)&&( ((order[_k][3]==1)&&(!(order[_k][2])))||((!(order[_k][3]))&&(order[_k][2]==1)) ) )
                      { gtypes[_i]='x'; break; }
                  }
                }
         else if (neighbors->list[_i].size==4)
                {
                if ( (!(order[_i][3]))&&( (order[_i][2]))&&( (!(gtypes[_i]))||(gtypes[_i]=='5') ) )
                  {//Look for P(5)-X=Y 
                  _j=neighbors->list[_i].size;
                  while (_j--)
                    if ( (*neighbors->list[_k=neighbors->list[_i].list[_j]].list!=_i)&&( ((order[_k][3]==1)&&(!(order[_k][2])))||((!(order[_k][3]))&&(order[_k][2]==1)) ) )
                      { gtypes[_i]='y'; break; }
                  }
                }
         } 
  else if (mol->a[_i]==CHEM_ATOM_TYPE_SULFUR)
         {
              if (neighbors->list[_i].size==3)
                {
                if ( (!(order[_i][3]))&&(order[_i][2]==1)&&( (!(gtypes[_i]))||(gtypes[_i]=='4') ) ) 
                  {//Look for S(4)-X=Y  
                  _j=neighbors->list[_i].size;
                  while (_j--)
                    if ( (*neighbors->list[_k=neighbors->list[_i].list[_j]].list!=_i)&&( ((order[_k][3]==1)&&(!(order[_k][2])))||((!(order[_k][3]))&&(order[_k][2]==1)) ) )
                      { gtypes[_i]='x'; break; }
                  }
                }
         else if (neighbors->list[_i].size==4)
                {
                if ( (!(order[_i][3]))&&( (order[_i][2]))&&( (!(gtypes[_i]))||(gtypes[_i]=='6') ) )
                  {//Look for S(6)-X=Y  
                  _j=neighbors->list[_i].size;
                  while (_j--)
                    if ( (*neighbors->list[_k=neighbors->list[_i].list[_j]].list!=_i)&&( ((order[_k][3]==1)&&(!(order[_k][2])))||((!(order[_k][3]))&&(order[_k][2]==1)) ) )
                      { gtypes[_i]='y'; break; }
                  }                   
                }
         }

//Stage IV. Compile of resonance chains using local A0(_)=X-Y=(_)A1 pattern
_t=mol->nedges;
while (_t--)
  if ( (mol->edges[_t].type=='1')&&( (isnan(mol->ff_b[_t].k)))&&(anchor_id[_i=mol->edges[_t].vertice[0]]!=anchor_id[_j=mol->edges[_t].vertice[1]])&&(order[_i][0]+order[_i][1]!=neighbors->list[_i].size)&&(order[_j][0]+order[_j][1]!=neighbors->list[_j].size) )
    {
    //Compose i-th GAFF suffix
    _c0=gtypes[_i];
    switch (mol->a[_i])
      {
      case CHEM_ATOM_TYPE_CARBON   : {
             if ( (order[_i][3]==1)&&(!order[_i][2]) ) 
               { if ( (!(inacycles->list[_i].size))&&( (!(gtypes[_i]))||(gtypes[_i]=='1') ) ) _c0='g'; } 
        else if ( (!(order[_i][3]))&&(order[_i][2]==1) )
               { if ( (!(inacycles->list[_i].size))&&( (!(gtypes[_i]))||(gtypes[_i]=='2') ) ) _c0='e'; } 
        break; } 
      case CHEM_ATOM_TYPE_NITROGEN : {
        if ( (neighbors->list[_i].size==2)&&(!(order[_i][3]))&&(order[_i][2]==1)&&(!(inacycles->list[_i].size))&&( (!(gtypes[_i]))||(gtypes[_i]=='2') ) ) _c0='e';
        break; }
      case CHEM_ATOM_TYPE_SILICON  : {
            if ( (order[_i][3]==1)&&(!(order[_i][2])) )
              { if ( (!(inacycles->list[_i].size))&&( (!(gtypes[_i]))||(gtypes[_i]=='1') ) ) _c0='g'; } 
        else if ( (!(order[_i][3]))&&(order[_i][2]==1) )
               { if ( (!(inacycles->list[_i].size))&&( (!(gtypes[_i]))||(gtypes[_i]=='2') ) ) _c0='e'; } 
        break; }
      case CHEM_ATOM_TYPE_PHOSPHOR : {
        if ( (neighbors->list[_i].size==2)&&(!(order[_i][3]))&&(order[_i][2]==1)&&(!(inacycles->list[_i].size))&&( (!(gtypes[_i]))||(gtypes[_i]=='2') ) ) _c0='e';
        break; } 
      }
    //Compose j-th GAFF suffix
    _c1=gtypes[_j];
    switch (mol->a[_j])
      {
      case CHEM_ATOM_TYPE_CARBON   : {
             if ( (order[_j][3]==1)&&(!order[_j][2]) ) 
               { if ( (!(inacycles->list[_j].size))&&( (!(gtypes[_j]))||(gtypes[_j]=='1') ) ) _c1='h'; } 
        else if ( (!(order[_j][3]))&&(order[_j][2]==1) )
               { if ( (!(inacycles->list[_j].size))&&( (!(gtypes[_j]))||(gtypes[_j]=='2') ) ) _c1='f'; } 
        break; } 
      case CHEM_ATOM_TYPE_NITROGEN : {
        if ( (neighbors->list[_j].size==2)&&(!(order[_j][3]))&&(order[_j][2]==1)&&(!(inacycles->list[_j].size))&&( (!(gtypes[_j]))||(gtypes[_j]=='2') ) ) _c1='f';
        break; }
      case CHEM_ATOM_TYPE_SILICON  : {
             if ( (order[_j][3]==1)&&(!(order[_j][2])) )
               { if ( (!(inacycles->list[_j].size))&&( (!(gtypes[_j]))||(gtypes[_j]=='1') ) ) _c1='h'; } 
        else if ( (!(order[_j][3]))&&(order[_j][2]==1) )
               { if ( (!(inacycles->list[_j].size))&&( (!(gtypes[_j]))||(gtypes[_j]=='2') ) ) _c1='f'; } 
        break; }
      case CHEM_ATOM_TYPE_PHOSPHOR : { 
        if ( (neighbors->list[_j].size==2)&&(!(order[_j][3]))&&(order[_j][2]==1)&&(!(inacycles->list[_j].size))&&( (!(gtypes[_j]))||(gtypes[_j]=='2') ) ) _c1='f';
        break; } 
      }
    //parameterize key bond
    if (!(_c=parameterize_bonds_gaff(TRUE,mol->edges[_k].type,&mol->ff_b[_t],_c0,_c1,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
    else if (_c==NTNF) count++; 
    //parameterize angles
    _l=mol->size_g;
    while (_l--)
           if (mol->ff_g[_l].atom[1]==_i)  
             {
             if ( (isnan(mol->ff_g[_l].k)))
               {
                    if (mol->ff_g[_l].atom[0]==_j)
                      {
                      if (!(_c=parameterize_angles_gaff(TRUE,&mol->ff_g[_l],_c1,_c0,gtypes[mol->ff_g[_l].atom[2]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                      else if (_c==NTNF) count++; 
                      }
               else if (mol->ff_g[_l].atom[2]==_j)   
                      {
                      if (!(_c=parameterize_angles_gaff(TRUE,&mol->ff_g[_l],gtypes[mol->ff_g[_l].atom[0]],_c0,_c1,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                      else if (_c==NTNF) count++;
                      }
               }
             }
      else if (mol->ff_g[_l].atom[1]==_j)
             {
             if ( (isnan(mol->ff_g[_l].k)))
               {
                    if (mol->ff_g[_l].atom[0]==_i)
                      {
                      if (!(_c=parameterize_angles_gaff(TRUE,&mol->ff_g[_l],_c0,_c1,gtypes[mol->ff_g[_l].atom[2]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                      else if (_c==NTNF) count++;
                      }
               else if (mol->ff_g[_l].atom[2]==_i)
                      {
                      if (!(_c=parameterize_angles_gaff(TRUE,&mol->ff_g[_l],gtypes[mol->ff_g[_l].atom[0]],_c1,_c0,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                      else if (_c==NTNF) count++;
                      }
               }
             }
    //parameterize impropers, works but it is pointless IMHO
    _l=mol->size_i;
    while (_l--)
      if (!(mol->ff_i[_l].type))
        {//A<-(B-C-D) angle
             if (mol->ff_i[_l].atom[0]==_i)
               {
               if ( (isnan(mol->ff_i[_l].k))) 
                 {
                      if (mol->ff_i[_l].atom[1]==_j)
                        {
                        if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_l],_c0,_c1,gtypes[mol->ff_i[_l].atom[2]],gtypes[mol->ff_i[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                        else if (_c==NTNF) count++; 
                        }
                 else if (mol->ff_i[_l].atom[2]==_j)
                        {
                        if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_l],_c0,gtypes[mol->ff_i[_l].atom[1]],_c1,gtypes[mol->ff_i[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                        else if (_c==NTNF) count++;
                        }
                 else if (mol->ff_i[_l].atom[3]==_j)
                        {
                        if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_l],_c0,gtypes[mol->ff_i[_l].atom[1]],gtypes[mol->ff_i[_l].atom[2]],_c1,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                        else if (_c==NTNF) count++;
                        }
                 }
               }
        else if (mol->ff_i[_l].atom[0]==_j) 
               {
               if ( (isnan(mol->ff_i[_l].k)))
                 {
                      if (mol->ff_i[_l].atom[1]==_i)
                        {
                        if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_l],_c1,_c0,gtypes[mol->ff_i[_l].atom[2]],gtypes[mol->ff_i[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                        else if (_c==NTNF) count++;
                        }
                 else if (mol->ff_i[_l].atom[2]==_i)
                        {
                        if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_l],_c1,gtypes[mol->ff_i[_l].atom[1]],_c0,gtypes[mol->ff_i[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                        else if (_c==NTNF) count++;
                        }
                 else if (mol->ff_i[_l].atom[3]==_i)
                        {
                        if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_l],_c1,gtypes[mol->ff_i[_l].atom[1]],gtypes[mol->ff_i[_l].atom[2]],_c0,mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                        else if (_c==NTNF) count++;
                        }
                 }
               }
        }
      else
        {//A-B-C-D angle
             if ( (mol->ff_i[_l].atom[1]==_i)&&(mol->ff_i[_l].atom[2]==_j) )
               {
               if ( (isnan(mol->ff_i[_l].k)))
                 {
                 if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_l],gtypes[mol->ff_i[_l].atom[0]],_c0,_c1,gtypes[mol->ff_i[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                 else if (_c==NTNF) count++; 
                 }
               }
        else if ( (mol->ff_i[_l].atom[1]==_j)&&(mol->ff_i[_l].atom[2]==_i) )
               {
               if ( (isnan(mol->ff_i[_l].k)))
                 {
                 if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_l],gtypes[mol->ff_i[_l].atom[0]],_c1,_c0,gtypes[mol->ff_i[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
                 else if (_c==NTNF) count++; 
                 }
               }
        }
    //parameterize dihedrals          
    _l=mol->size_d;
    while (_l--)
      {
           if ( (mol->ff_d[_l].atom[1]==_i)&&(mol->ff_d[_l].atom[2]==_j) )
             {
             if ( (isnan(mol->ff_d[_l].k[0])))
               {
               if (!(_c=parameterize_dihedrals_gaff(TRUE,&mol->ff_d[_l],gtypes[mol->ff_d[_l].atom[0]],_c0,_c1,gtypes[mol->ff_d[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
               else if (_c==NTNF) count++; 
               }
             }
      else if ( (mol->ff_d[_l].atom[1]==_j)&&(mol->ff_d[_l].atom[2]==_i) )  
             {
             if ( (isnan(mol->ff_d[_l].k[0])))
               {
               if ((_c=parameterize_dihedrals_gaff(TRUE,&mol->ff_d[_l],gtypes[mol->ff_d[_l].atom[0]],_c1,_c0,gtypes[mol->ff_d[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
               else if (_c==NTNF) count++; 
               }
             }
      }
    }

//Stage V. Classic parameterization run
//Stage V.1. Parameterize classic bonds
_l=mol->size_b;
while (_l--)
  if ( (isnan(mol->ff_b[_l].k)))
    {
    if (!(_c=parameterize_bonds_gaff(TRUE,mol->edges[_l].type,&mol->ff_b[_l],gtypes[mol->ff_b[_l].atom[0]],gtypes[mol->ff_b[_l].atom[1]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
    else if (_c==NTNF) count++; 
    }
//Stage V.2. Parameterize classic angles
_l=mol->size_g;
while (_l--)
  if ( (isnan(mol->ff_g[_l].k)))
    {
    if (!(_c=parameterize_angles_gaff(TRUE,&mol->ff_g[_l],gtypes[mol->ff_g[_l].atom[0]],gtypes[mol->ff_g[_l].atom[1]],gtypes[mol->ff_g[_l].atom[2]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
    else if (_c==NTNF) count++; 
    }
//Stage V.3. Parameterize classic impropers
_l=mol->size_i;
while (_l--)
  if ( (isnan(mol->ff_i[_l].k)))
    {
    if (!(_c=parameterize_impropers_gaff(TRUE,&mol->ff_i[_l],gtypes[mol->ff_i[_l].atom[0]],gtypes[mol->ff_i[_l].atom[1]],gtypes[mol->ff_i[_l].atom[2]],gtypes[mol->ff_i[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE;
    else if (_c==NTNF) count++;
    }
//Stage V.4. Parameterize classic dihedrals
_l=mol->size_d;
while (_l--)
  if ( (isnan(mol->ff_d[_l].k[0])))
    {
    if (!(_c=parameterize_dihedrals_gaff(TRUE,&mol->ff_d[_l],gtypes[mol->ff_d[_l].atom[0]],gtypes[mol->ff_d[_l].atom[1]],gtypes[mol->ff_d[_l].atom[2]],gtypes[mol->ff_d[_l].atom[3]],mol,top))) goto LABEL_PARAMETERIZATION_FAILURE; 
    else if (_c==NTNF) count++;
    }

//Everything is done. Free memory and exit.
free(gtypes); gtypes=0x0;
free(temp); temp=0x0;
return count;
LABEL_PARAMETERIZATION_FAILURE: ylib_errno=YERROR_EXTERNAL_CODE; 
free(gtypes); gtypes=0x0;
free(temp); temp=0x0;
return (unsigned int)-1;
}

//This function adjust parameters values to general logic, runtime problems and (optional) match given geometry
unsigned int adjust_mol_parameters(char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_clist *inacycles,t_mol *mol,t_vec *r,t_top *top)
{
register unsigned int _i, _j, _k, _l, count=0;
union { double alpha; char flag; } a;

//Stage I. Remove non-existing dihedrals (i.e. at open angles) X0 <- C_=_C -> X1
//Stage I.1. Remove impropers
_i=mol->size_i;
while (_i--)
  if ( (mol->ff_i[_i].type))
    {
    //Check X0-A-B angle
    _k=mol->ff_i[_i].atom[1], _j=mol->size_g;
    while (_j--)
      if ( (mol->ff_g[_j].atom[1]==_k)&&( ( (mol->ff_g[_j].atom[0]==mol->ff_i[_i].atom[0])&&(mol->ff_g[_j].atom[2]==mol->ff_i[_i].atom[2]) )||
                                          ( (mol->ff_g[_j].atom[2]==mol->ff_i[_i].atom[0])&&(mol->ff_g[_j].atom[0]==mol->ff_i[_i].atom[2]) ) ) )
        {
        if ( (mol->ff_g[_j].v<PI*005./180.)||(mol->ff_g[_j].v>PI*175./180.) ) // value { 0, PI } +/-5deg
          { REMOVE_OPEN_IDIH: if (!(--mol->size_i)) { free(mol->ff_i), mol->ff_i=0x0; } else memcpy(&mol->ff_i[_i],&mol->ff_i[_i+1],sizeof(t_ff_i)*(mol->size_i-_i)); goto NEXT_IDIH; }
        else break;
        }
    //Check A-B-X1 angle
    _k=mol->ff_i[_i].atom[2], _j=mol->size_g;
    while (_j--)
      if ( (mol->ff_g[_j].atom[1]==_k)&&( ( (mol->ff_g[_j].atom[0]==mol->ff_i[_i].atom[1])&&(mol->ff_g[_j].atom[2]==mol->ff_i[_i].atom[3]) )||
                                          ( (mol->ff_g[_j].atom[2]==mol->ff_i[_i].atom[1])&&(mol->ff_g[_j].atom[0]==mol->ff_i[_i].atom[3]) ) ) )
        {
        if ( (mol->ff_g[_j].v<PI*005./180.)||(mol->ff_g[_j].v>PI*175./180.) ) goto REMOVE_OPEN_IDIH; // value { 0, PI } +/-5deg
        else break;
        }
    NEXT_IDIH : ; 
    }
//Stage I.2. Remove torsions
_i=mol->size_d;
while (_i--)
  if ( (mol->ff_d[_i].d))
    {
    //Check X0-A-B angle
    _k=mol->ff_d[_i].atom[1], _j=mol->size_g;
    while (_j--)
      if ( (mol->ff_g[_j].atom[1]==_k)&&( ( (mol->ff_g[_j].atom[0]==mol->ff_d[_i].atom[0])&&(mol->ff_g[_j].atom[2]==mol->ff_d[_i].atom[2]) )||
                                          ( (mol->ff_g[_j].atom[2]==mol->ff_d[_i].atom[0])&&(mol->ff_g[_j].atom[0]==mol->ff_d[_i].atom[2]) ) ) )
        {
        if ( (mol->ff_g[_j].v<PI*005./180.)||(mol->ff_g[_j].v>PI*175./180.) ) // value { 0, PI } +/-5deg
          { REMOVE_OPEN_DDIH: if (!(--mol->size_d)) { free(mol->ff_d), mol->ff_d=0x0; } else memcpy(&mol->ff_d[_i],&mol->ff_d[_i+1],sizeof(t_ff_d)*(mol->size_d-_i)); goto NEXT_DDIH; }
        else break;
        }
    //Check A-B-X1 angle
    _k=mol->ff_d[_i].atom[2], _j=mol->size_g;
    while (_j--)
      if ( (mol->ff_g[_j].atom[1]==_k)&&( ( (mol->ff_g[_j].atom[0]==mol->ff_d[_i].atom[1])&&(mol->ff_g[_j].atom[2]==mol->ff_d[_i].atom[3]) )||
                                          ( (mol->ff_g[_j].atom[2]==mol->ff_d[_i].atom[1])&&(mol->ff_g[_j].atom[0]==mol->ff_d[_i].atom[3]) ) ) )
        {
        if ( (mol->ff_g[_j].v<PI*005./180.)||(mol->ff_g[_j].v>PI*175./180.) ) goto REMOVE_OPEN_DDIH; // value { 0, PI } +/-5deg
        else break;
        }
    NEXT_DDIH : ; 
    }

//Stage II. Set all out-of-cycles to their initial isomer values
if ( (r))
  {
  _i=mol->size_i;
  while (_i--)
    if (!(mol->ff_i[_i].type))
      {//Adjust optical isomers to the geometry
      if (calc_dih_angle_value(&r[mol->ff_i[_i].atom[0]],&r[mol->ff_i[_i].atom[1]],&r[mol->ff_i[_i].atom[2]],&r[mol->ff_i[_i].atom[3]])<0.) 
        { _l=mol->ff_i[_i].atom[2], mol->ff_i[_i].atom[2]=mol->ff_i[_i].atom[3], mol->ff_i[_i].atom[3]=_l; }
      }
    else
      {//A-B-C-D case
      if ( (!(inacycles->list[mol->ff_i[_i].atom[1]].size))||(!(inacycles->list[mol->ff_i[_i].atom[2]].size))||(!(overlap_list(FALSE,&inacycles->list[mol->ff_i[_i].atom[1]],FALSE,&inacycles->list[mol->ff_i[_i].atom[2]]))) )
        {
        //Exclude amide bonds which conformation is a priory determined 
        if ( ( (mol->a[mol->ff_i[_i].atom[1]]==CHEM_ATOM_TYPE_CARBON)&&(mol->a[mol->ff_i[_i].atom[2]]==CHEM_ATOM_TYPE_NITROGEN) )||
             ( (mol->a[mol->ff_i[_i].atom[2]]==CHEM_ATOM_TYPE_CARBON)&&(mol->a[mol->ff_i[_i].atom[1]]==CHEM_ATOM_TYPE_NITROGEN) ) ) 
          {
          _j=mol->nedges; do { _j--; } while( ( (mol->edges[_j].vertice[0]!=mol->ff_i[_i].atom[1])||((mol->edges[_j].vertice[1]!=mol->ff_i[_i].atom[2])) )&&( (mol->edges[_j].vertice[0]!=mol->ff_i[_i].atom[2])||((mol->edges[_j].vertice[1]!=mol->ff_i[_i].atom[1])) ) );
          if (mol->edges[_j].type==(int)'m') continue;
          }
        a.alpha=calc_dih_angle_value(&r[mol->ff_i[_i].atom[0]],&r[mol->ff_i[_i].atom[1]],&r[mol->ff_i[_i].atom[2]],&r[mol->ff_i[_i].atom[3]]);
             if ( (mol->ff_i[_i].v>-0.1*PI)&&(mol->ff_i[_i].v<+0.1*PI) )             { if (fabs(a.alpha)>PI/2.) { mol->ff_i[_i].v=PI, count++; } } //Tune isomer value
        else if ( (fabs(mol->ff_i[_i].v)>0.9*PI)&&(fabs(mol->ff_i[_i].v)<PI+SMALL) ) { if (fabs(a.alpha)<PI/2.) { mol->ff_i[_i].v=0., count++; } } //Tune isomer value
        }
      }
  }

//Exit 
return count;
}



//This function compiles a mechanical model of a molecule from the scratch (i.e. str)
t_mol *compile_mol(char verbose,double ph,char (*order)[4],t_clist *neighbors,unsigned int **_anchor_id,t_adjacency **_adjacency,t_vec *r,t_str *str,t_top *top)
{
register unsigned int _i;
unsigned int *anchor_id;
t_clist *cycles, *_neighbors, *inacycles;
t_mol *mol;
t_adjacency *adjacency;

//Now compile ymol
if ( (verbose))
  {
  if (!(_neighbors=copy_clist(neighbors))) 
    { 
    yprintf(YPRINTF_ERROR,"error in copy_clist(), err_code=%s.\n",get_yerrno(ylib_errno));
    LABEL_ERROR_1: return FALSE; }
  //Then compile ymols
  if (!(cycles=define_cycles_with_neighbors(_neighbors)))
    { 
    yprintf(YPRINTF_ERROR,"error in define_cycles_with_neighbors(), err_code=%s.\n",get_yerrno(ylib_errno));
    free(_neighbors); _neighbors=0x0; goto LABEL_ERROR_1; }
  else 
    {
    yprintf(YPRINTF_INFO,"define_cycles_with_neighbors() found %d cycles in the molecular graph.\n",cycles->size);
    free(_neighbors); _neighbors=0x0;
    }
  if ((_i=classify_cycles(order,neighbors,cycles,str))==(unsigned int)-1) 
    { 
    yprintf(YPRINTF_ERROR,"error in classify_cycles(), err_code=%s.\n",get_yerrno(ylib_errno));
    LABEL_ERROR_2: free(cycles); cycles=0x0; goto LABEL_ERROR_1; }
  else yprintf(YPRINTF_INFO,"classify_cycles() has detected %d aromatic rings.\n",_i);
  //First order compiler
  if (!(mol=compile_mol_YFF1(order,neighbors,_i,cycles,str,top))) 
    { yprintf(YPRINTF_ERROR,"error in compile_mol_YFF1(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_2; }
  else { mol->cycles=cycles; yprintf(YPRINTF_INFO,"compile_mol_YFF1() has successfully constructed model of nonboned interaction of the molecule.\n",_i); }
  if ((_i=ionize_mol_empirically(order,neighbors,mol,top))==(unsigned int)-1) 
    { 
    yprintf(YPRINTF_ERROR,"error in ionize_mol_empirically(), err_code=%s.\n",get_yerrno(ylib_errno));
    LABEL_ERROR_3: free_mol(mol); mol=0x0; goto LABEL_ERROR_2; }
  else yprintf(YPRINTF_INFO,"ionize_mol_empirically() has added %d virtual counterions into molecular graph.\n",_i);
  //Second order compiler
  if (!(disassemble_mol(&anchor_id,neighbors,mol))) 
    {  yprintf(YPRINTF_ERROR,"error in disassemble_mol(), err_code=%s.\n",get_yerrno(ylib_errno)); goto LABEL_ERROR_3; }
  else yprintf(YPRINTF_INFO,"disassemble_mol() has successfully disassemblem kinestatic model of the molecule into a set of independent anchors.\n");
  if (!(adjacency=construct_adjacency(anchor_id,mol->anchors,mol->natoms,mol->nedges,mol->edges)))
    {
    yprintf(YPRINTF_ERROR,"error in construct_adjacency(), err_code=%s.\n",get_yerrno(ylib_errno)); 
    LABEL_ERROR_4: free(anchor_id); anchor_id=0x0; goto LABEL_ERROR_3; }
  else yprintf(YPRINTF_INFO,"construct_adjacencyl() has successfully constructed anchors connectivity net.\n");
  if (!(compose_mol(order,neighbors,anchor_id,adjacency,&inacycles,mol,top)))
    { 
    yprintf(YPRINTF_ERROR,"error in compose_mol(), err_code=%s.\n",get_yerrno(ylib_errno));
    LABEL_ERROR_5: free(adjacency); adjacency=0x0; goto LABEL_ERROR_4; }
  else yprintf(YPRINTF_INFO,"compose_mol() has successfully constructed molecular mechanic model of the bonded molecular interactions.\n");
  if ((_i=parameterize_mol_GAFF(order,neighbors,anchor_id,inacycles,mol,top))==(unsigned int)-1)
    { yprintf(YPRINTF_ERROR,"error in parameterize_mol_GAFF(), err_code=%s.\n",get_yerrno(ylib_errno)); free(inacycles); goto LABEL_ERROR_5; }
  else yprintf(YPRINTF_INFO,"parameterize_mol_GAFF() has successfully constructed rough parameterization of the bonded molecur interactions (using %1.d generic parameters out of %1.d total).\n",_i,mol->size_b+mol->size_g+mol->size_i+mol->size_d); 
  //Check if molecular geometry is 'solvable'
  _i=mol->natoms; while (_i--) if ( (neighbors->list[_i].size>1)&&( ( (isnan(mol->r[_i].i)))||( (isnan(mol->r[_i].j)))||( (isnan(mol->r[_i].k))) ) ) { yprintf(YPRINTF_ERROR,"there are multivalent atoms with undefined coordinates, but denovo construction is not (yet) implemented.\n"); ylib_errno=YERROR_NIMPLEMENTED; free(inacycles); goto LABEL_ERROR_5; } //We can resolve only onevalent nans
  if ((_i=resolve_onevalent(neighbors,mol->r,mol))==(unsigned int)-1) 
    { yprintf(YPRINTF_ERROR,"error in resolve_onevalent(), err_code=%s.\n",get_yerrno(ylib_errno)); free(inacycles); goto LABEL_ERROR_5; }
  //Adjust molecular geometry
  if ((_i=adjust_mol_parameters(order,neighbors,anchor_id,inacycles,mol,r,top))==(unsigned int)-1)
    { yprintf(YPRINTF_ERROR,"error in adjust_mol_impropers(), err_code=%s.\n",get_yerrno(ylib_errno)); free(inacycles); goto LABEL_ERROR_5; }
  else { yprintf(YPRINTF_INFO,"adjust_mol_impropers() has successfully edited %1d improper values (out of %1d total).\n",_i,mol->size_i); }
  }
else 
  {
  //Process STR
  if (!(_neighbors=copy_clist(neighbors))) goto LABEL_ERROR_1;
  if (!(cycles=define_cycles_with_neighbors(_neighbors))) { free(_neighbors); _neighbors=0x0; goto LABEL_ERROR_1; }
  else { free(_neighbors); _neighbors=0x0; }  
  if ((_i=classify_cycles(order,neighbors,cycles,str))==(unsigned int)-1) { free(cycles); cycles=0x0; goto LABEL_ERROR_2; } 
  //First order compiler
  if (!(mol=compile_mol_YFF1(order,neighbors,_i,cycles,str,top))) { free(cycles); cycles=0x0; goto LABEL_ERROR_2; }
  if ((_i=ionize_mol_empirically(order,neighbors,mol,top))==(unsigned int)-1) {free(cycles); cycles=0x0; goto LABEL_ERROR_2; }
  //Second order compiler
  if (!(disassemble_mol(&anchor_id,neighbors,mol))) goto LABEL_ERROR_3;
  if (!(adjacency=construct_adjacency(anchor_id,mol->anchors,mol->natoms,mol->nedges,mol->edges))) goto LABEL_ERROR_3;
  if (!(compose_mol(order,neighbors,anchor_id,adjacency,&inacycles,mol,top))) goto LABEL_ERROR_4; 
  if ((_i=parameterize_mol_GAFF(order,neighbors,anchor_id,inacycles,mol,top))==(unsigned int)-1) { free(inacycles); goto LABEL_ERROR_5; }
  //Check if molecular geometry is 'solvable'
  _i=mol->natoms; while (_i--) if ( (neighbors->list[_i].size>1)&&( ( (isnan(mol->r[_i].i)))||( (isnan(mol->r[_i].j)))||( (isnan(mol->r[_i].k))) ) ) { ylib_errno=YERROR_NIMPLEMENTED; free(inacycles); goto LABEL_ERROR_5; } //We can resolve only onevalent nans
  if ((_i=resolve_onevalent(neighbors,mol->r,mol))==(unsigned int)-1) { free(inacycles); goto LABEL_ERROR_5; }
  //Adjust molecular geometry
  if ((_i=adjust_mol_parameters(order,neighbors,anchor_id,inacycles,mol,r,top))==(unsigned int)-1) { free(inacycles); goto LABEL_ERROR_5; } 
  }

//Save meta-data if needded
//if ( (_inacycles)) (*_inacycles)=inacycles; else { free(inacycles); inacycles=0x0; }
free(inacycles);
if ( (_anchor_id)) (*_anchor_id)=anchor_id; else { free(anchor_id); anchor_id=0x0; }
if ( (_adjacency)) (*_adjacency)=adjacency; else { free(adjacency); adjacency=0x0; }
//Everything is cool, exiting
return mol;
}





// ***********************   M O L E C U L E S    A N D    S U B S T R U C T U R E S     C O N S T R U C T I O N    P A R T   ************************************  
//NB!! Only single bond can separate two anchors
//NB!! Each gag starts with root atom, the edges are [acceptor_id],[root_id]:[1]



//This function construct water molecule t_mol structure
t_str *construct_wat(t_str *str,t_top *top)
{
if (!(str=alloc_str(str,0x1,0x3,0x2))) return FALSE;
str->rsize=((unsigned int*)YMOL_WAT_NAME);
str->a[0]=CHEM_ATOM_TYPE_OXYGEN, str->a[1]=CHEM_ATOM_TYPE_HYDROGEN, str->a[2]=CHEM_ATOM_TYPE_HYDROGEN; 
str->edges[0].vertice[0]=0, str->edges[0].vertice[1]=1, str->edges[0].type='1'; 
str->edges[1].vertice[0]=0, str->edges[1].vertice[1]=2, str->edges[1].type='1'; 
return str;
}


//Gags trigger works as follows:
// | geometry bit | topology bit | output bit* |
// |       1/0    |       1/0    |     1/0     |
// *Means 0 - output atoms number, 1 - output bonds number

//This is service function to define approximate single bond lenght between two atom with defined chemical types
//The supported types H,C,N,O,F,Si,P,S,Cl,Br,I; all other pairs are very approximately detemined on their row in periodic table
double get_approximate_single_bond_lenght(char chem_type_1,char chem_type_2)
{
register char _c;
register double _d;
if (chem_type_1>chem_type_2) { _c=chem_type_1, chem_type_1=chem_type_2, chem_type_2=_c; }
switch (chem_type_1)
  {
  case CHEM_ATOM_TYPE_HYDROGEN : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_HYDROGEN : return 0.74; //H-H  is 0.74A
      case CHEM_ATOM_TYPE_CARBON   : return 1.09; //H-C  is 1.09A
      case CHEM_ATOM_TYPE_NITROGEN : return 1.01; //H-N  is 1.01A
      case CHEM_ATOM_TYPE_OXYGEN   : return 0.96; //H-O  is 0.96A
      case CHEM_ATOM_TYPE_FLUORINE : return 0.92; //H-F  is 0.92A
      case CHEM_ATOM_TYPE_SILICON  : return 1.43; //H-Si is 1.43A
      case CHEM_ATOM_TYPE_PHOSPHOR : return 1.38; //H-P  is 1.38A
      case CHEM_ATOM_TYPE_SULFUR   : return 1.34; //H-S  is 1.34A
      case CHEM_ATOM_TYPE_CHLORINE : return 1.31; //H-Cl is 1.31A
      case CHEM_ATOM_TYPE_BROMINE  : return 1.46; //H-Br is 1.46A
      case CHEM_ATOM_TYPE_IODINE   : return 1.65; //H-I  is 1.62A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_CARBON : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_CARBON   : return 1.54; //C-C  is 1.54A
      case CHEM_ATOM_TYPE_NITROGEN : return 1.47; //C-N  is 1.47A
      case CHEM_ATOM_TYPE_OXYGEN   : return 1.43; //C-O  is 1.43A
      case CHEM_ATOM_TYPE_FLUORINE : return 1.36; //C-F  is 1.36A
      case CHEM_ATOM_TYPE_SILICON  : return 1.88; //C-Si is 1.88A
      case CHEM_ATOM_TYPE_PHOSPHOR : return 1.83; //C-P  is 1.83A
      case CHEM_ATOM_TYPE_SULFUR   : return 1.79; //C-S  is 1.79A
      case CHEM_ATOM_TYPE_CHLORINE : return 1.76; //C-Cl is 1.76A
      case CHEM_ATOM_TYPE_BROMINE  : return 1.91; //C-Br is 1.91A
      case CHEM_ATOM_TYPE_IODINE   : return 2.10; //C-I  is 2.10A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_NITROGEN : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_NITROGEN : return 1.45; //N-N  is 1.45A
      case CHEM_ATOM_TYPE_OXYGEN   : return 1.36; //N-O  is 1.36A
      case CHEM_ATOM_TYPE_FLUORINE : return 1.34; //N-F  is 1.34A
      case CHEM_ATOM_TYPE_SILICON  : return 1.86; //N-Si is 1.86A
      case CHEM_ATOM_TYPE_PHOSPHOR : return 1.81; //N-P  is 1.81A
      case CHEM_ATOM_TYPE_SULFUR   : return 1.77; //N-S  is 1.77A
      case CHEM_ATOM_TYPE_CHLORINE : return 1.74; //N-Cl is 1.74A
      case CHEM_ATOM_TYPE_BROMINE  : return 1.89; //N-Br is 1.89A
      case CHEM_ATOM_TYPE_IODINE   : return 2.08; //N-I  is 2.08A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_OXYGEN : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_OXYGEN   : return 1.48; //O-O  is 1.48A
      case CHEM_ATOM_TYPE_FLUORINE : return 1.45; //O-F  is 1.45A
      case CHEM_ATOM_TYPE_SILICON  : return 1.84; //O-Si is 1.84A
      case CHEM_ATOM_TYPE_PHOSPHOR : return 1.79; //O-P  is 1.79A
      case CHEM_ATOM_TYPE_SULFUR   : return 1.75; //O-S  is 1.75A
      case CHEM_ATOM_TYPE_CHLORINE : return 1.72; //O-Cl is 1.72A
      case CHEM_ATOM_TYPE_BROMINE  : return 1.87; //O-Br is 1.87A
      case CHEM_ATOM_TYPE_IODINE   : return 2.06; //O-I  is 2.06A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_FLUORINE : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_FLUORINE : return 1.42; //F-F  is 1.42A
      case CHEM_ATOM_TYPE_SILICON  : return 1.83; //F-Si is 1.83A
      case CHEM_ATOM_TYPE_PHOSPHOR : return 1.78; //F-P  is 1.78A
      case CHEM_ATOM_TYPE_SULFUR   : return 1.74; //F-S  is 1.74A
      case CHEM_ATOM_TYPE_CHLORINE : return 1.71; //F-Cl is 1.71A
      case CHEM_ATOM_TYPE_BROMINE  : return 1.86; //F-Br is 1.86A
      case CHEM_ATOM_TYPE_IODINE   : return 2.05; //F-I  is 2.05A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_SILICON : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_SILICON  : return 2.22; //Si-Si is 2.22A
      case CHEM_ATOM_TYPE_PHOSPHOR : return 2.17; //Si-P  is 2.17A
      case CHEM_ATOM_TYPE_SULFUR   : return 2.13; //Si-S  is 2.13A
      case CHEM_ATOM_TYPE_CHLORINE : return 2.10; //Si-Cl is 2.10A
      case CHEM_ATOM_TYPE_BROMINE  : return 2.25; //Si-Br is 2.25A
      case CHEM_ATOM_TYPE_IODINE   : return 2.44; //Si-I  is 2.44A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_PHOSPHOR : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_PHOSPHOR : return 2.12; //P-P  is 2.12A
      case CHEM_ATOM_TYPE_SULFUR   : return 2.08; //P-S  is 2.08A
      case CHEM_ATOM_TYPE_CHLORINE : return 2.05; //P-Cl is 2.05A
      case CHEM_ATOM_TYPE_BROMINE  : return 2.20; //P-Br is 2.20A
      case CHEM_ATOM_TYPE_IODINE   : return 2.39; //P-I  is 2.39A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_SULFUR : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_SULFUR   : return 2.04; //S-S  is 2.04A
      case CHEM_ATOM_TYPE_CHLORINE : return 2.01; //S-Cl is 2.01A
      case CHEM_ATOM_TYPE_BROMINE  : return 2.16; //S-Br is 2.16A
      case CHEM_ATOM_TYPE_IODINE   : return 2.35; //S-I  is 2.35A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_CHLORINE : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_CHLORINE : return 1.98; //Cl-Cl is 1.98A
      case CHEM_ATOM_TYPE_BROMINE  : return 2.13; //Cl-Br is 2.13A
      case CHEM_ATOM_TYPE_IODINE   : return 2.32; //Cl-I  is 2.32A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_BROMINE : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_BROMINE  : return 2.28; //Br-Br is 2.28A
      case CHEM_ATOM_TYPE_IODINE   : return 2.47; //Br-I  is 2.47A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  case CHEM_ATOM_TYPE_IODINE : {
    switch (chem_type_2)
      {
      case CHEM_ATOM_TYPE_IODINE   : return 2.66; //I-I  is 2.66A
      default                      : { goto TABLE_ROW_ESTIMATE; }  //Select using row of periodic table - very approximate
      }
    }
  }
TABLE_ROW_ESTIMATE: ; //Select using row of periodic table - very approximate
     if (chem_type_1<3)  _d=1.0;
else if (chem_type_1<11) _d=1.5;
else if (chem_type_1<19) _d=2.0;
else if (chem_type_1<37) _d=2.3;
else if (chem_type_1<55) _d=2.7;
else                     _d=2.9;
     if (chem_type_2<3)  _d+=1.0;
else if (chem_type_2<11) _d+=1.5;
else if (chem_type_2<19) _d+=2.0;
else if (chem_type_2<37) _d+=2.3;
else if (chem_type_2<55) _d+=2.7;
else                     _d+=2.9; 
return _d/2.;
}

//This is service function for gag construction which is used to align default - OZ oriented - gag onto bonding vector
//Note n - is amount of atoms to transform, r - coordinates to transform, r0 - reference aceptor, v - accepting bond vector, bond_lenght - desired bond lenght 
void _construct_gag_geometry(unsigned int n,t_vec *r,t_vec *r0,t_vec *v,double bond_lenght)
{
register unsigned int _i;
register double _d;
double csA, snA;
t_tensor R;
t_vec u;

//Check the direction vector 
if ((_d=sqrt(v->i*v->i+v->j*v->j+v->k*v->k))<SMALL)
  { _i=n; while (_i--) { r[_i].i+=r0->i, r[_i].j+=r0->j, r[_i].k+=r0->k; } } //Vector v doesn't exists
else
  {//Vector v does exist
  v->i/=_d, v->j/=_d, v->k/=_d; //Normalize the vector for origin shift below
  //Construct rotation tensor (don't forget that the vectors directions should be opposite!!!)
  csA=v->k; 
       if (csA>+1.-SMALL2)
         { R[0][0]=R[1][1]=+1., R[2][2]=-1., R[0][1]=R[0][2]=R[1][0]=R[1][2]=R[2][0]=R[2][1]=0.; } //Co-directed v and OZ
  else if (csA<-1.+SMALL2)
         { R[0][0]=R[1][1]=+1., R[2][2]=+1., R[0][1]=R[0][2]=R[1][0]=R[1][2]=R[2][0]=R[2][1]=0.; } //Anti-directed v and OZ
  else 
    {
    //Calculate u vector as u = [v X [OZ]]/|[v X [OZ]]|
    u.i=-v->j, u.j=+v->i, u.k=0.;
    _d=sqrt(u.i*u.i+u.j*u.j);
    u.i/=_d, u.j/=_d;
    //Note the sign of the rotation angle is determined using Determinant of [v, [OZ], u] matrix
    snA=sqrt(1.-csA*csA);
    _d=v->j*u.i-v->i*u.j;
    if (_d>0.) 
      { rotate_around_uvector(&R,&u,csA,+snA); } //Do -Sin rotation
    else 
      { rotate_around_uvector(&R,&u,csA,-snA); } //Do +Sin rotation
    }
  //Apply the rotation
  _i=n;  while (_i--) { multiple_origin_tensor_origin_vec(&u,&R,&r[_i]), r[_i].i=u.i, r[_i].j=u.j, r[_i].k=u.k; }
  //Shift the origin
  multiple_vec_scalar(&u,v,bond_lenght-1.), u.i+=r0->i, u.j+=r0->j, u.k+=r0->k;
  _i=n; while (_i--) { r[_i].i+=u.i, r[_i].j+=u.j, r[_i].k+=u.k; }
  }
}

// C A R B O N    G A G S
//This function makes carboxyl gag insted of an atom0->atom1 :: atom0->COO(-)
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(C)-C = 1.55342A, C-O = 1.23445A, X(C)-C-O = 115.3Deg, O-C-O = 129.4Deg
//X   0.0000  0.0000  0.0000
//C   0.0000  0.0000  1.0000 (1.5268)
//O1  0.0000  1.4044  1.6639
//O2  0.0000 -1.4044  1.6639 
void make_carboxyl_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if ( (mode&2))
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"cxyg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='C';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='O';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='O';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id),   str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id),   str->edges[(*edge_id)+1].type='2';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id),   str->edges[(*edge_id)+2].type='a'; 
  }
//Stage II. Resolve geometry if needed
if ( (mode&4))
  {
  //Setup default methyl coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.0000, str->r[(*atom_id)+0].k=+1.0000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+1.4044, str->r[(*atom_id)+1].k=+1.6639;
  str->r[(*atom_id)+2].i=+0.0000, str->r[(*atom_id)+2].j=-1.4044, str->r[(*atom_id)+2].k=+1.6639;
  //Construct gag geometry
  _construct_gag_geometry(3,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_CARBON));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"cxyg");
  if ( (atom_id)) (*atom_id)=3;
  if ( (edge_id)) (*edge_id)=3;
  }
}//Exit
//This function makes cyane gag insted of an atom0->atom1 :: atom0->[NH2-C(+)=NH2]
//It adds 7 atoms and 7 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-C = 1.49999A, C-N = 1.3085A, X(C)-C-N = 120.0Deg, N-H = 0.99803A, C-N-H = 121.6Deg, all in one plane
//X    0.0000  0.00000  0.00000
//C    0.0000  0.00000  1.00000 (1.49999)
//N1   0.0000  1.13319  1.65425
//H11  0.0000  1.16106  2.65189
//H12  0.0000  2.01111  2.12894
//N2   0.0000 -1.13319  1.65425
//H21  0.0000 -1.16106  2.65189
//H22  0.0000 -2.01111  2.12894
void make_cguanidine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if ( (mode&2))
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"cgug");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='C';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='N';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='N';
  str->a[(*atom_id)+5]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+5][0]=str->anames[(*atom_id)+5][1]=str->anames[(*atom_id)+5][2]=' ', str->anames[(*atom_id)+5][3]='H';
  str->a[(*atom_id)+6]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+6][0]=str->anames[(*atom_id)+6][1]=str->anames[(*atom_id)+6][2]=' ', str->anames[(*atom_id)+6][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0,   str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0,   str->edges[(*edge_id)+1].type='2';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+1,   str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1,   str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+0,   str->edges[(*edge_id)+1].type='a';
  str->edges[(*edge_id)+5].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+4,   str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+6].vertice[0]=(*atom_id)+6, str->edges[(*edge_id)+6].vertice[1]=(*atom_id)+4,   str->edges[(*edge_id)+1].type='1';
  }
//Stage II. Resolve geometry if needed
if ( (mode&4))
  {
  //Setup default methyl coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+1.13319, str->r[(*atom_id)+1].k=+1.65425;
  str->r[(*atom_id)+2].i=+0.0000, str->r[(*atom_id)+2].j=+1.16106, str->r[(*atom_id)+2].k=+2.65189;
  str->r[(*atom_id)+3].i=+0.0000, str->r[(*atom_id)+3].j=+2.01111, str->r[(*atom_id)+3].k=+2.12894;
  str->r[(*atom_id)+4].i=+0.0000, str->r[(*atom_id)+4].j=-1.13319, str->r[(*atom_id)+4].k=+1.65425;
  str->r[(*atom_id)+5].i=+0.0000, str->r[(*atom_id)+5].j=-1.16106, str->r[(*atom_id)+5].k=+2.65189;
  str->r[(*atom_id)+6].i=+0.0000, str->r[(*atom_id)+6].j=-2.01111, str->r[(*atom_id)+6].k=+2.12894;
  //Construct gag geometry
  _construct_gag_geometry(7,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_CARBON));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"cgug");
  if ( (atom_id)) (*atom_id)=7;
  if ( (edge_id)) (*edge_id)=7;
  }
}//Exit
//This function makes cyane gag insted of an atom0->atom1 :: atom0->C_=_N
//It adds 2 atoms and 2 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all on OZ line, X(C)-C = 1.467A, C-N = 1.13478A, X(C)-C-N = 180.0Deg
//X   0.0000  0.0000  0.0000
//C   0.0000  0.0000  1.0000 (1.5268)
//N   0.0000  0.0000  2.1348
void make_cyane_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if ( (mode&2))
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"cyng");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='C';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='N';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id),   str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id),   str->edges[(*edge_id)+1].type='1';
  }
//Stage II. Resolve geometry if needed
if ( (mode&4))
  {
  //Setup default methyl coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.0000, str->r[(*atom_id)+0].k=+1.0000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+0.0000, str->r[(*atom_id)+1].k=+2.1348;
  //Construct gag geometry
  _construct_gag_geometry(2,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_CARBON));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"cyng");
  if ( (atom_id)) (*atom_id)=2;
  if ( (edge_id)) (*edge_id)=2;
  }
}//Exit
//This function makes amide gag insted of an atom0->atom1 :: atom0->CO-NH2
//It adds 5 atoms and 5 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz-plane, X(C)-C = 1.51314A, C-O = 1.19775A, X(C)-C-O = 122.8Deg, C-N = 1.35592A, X(C)-C-N = 115.0Deg, O-C-N = 122.2Deg, N-H = 0.9920A, C-N-H = 120.4Deg, H-N-H = 119.2Deg 
//X   0.0000  0.0000  0.0000
//C   0.0000  0.0000  1.0000 (1.51314)
//O   0.0000  1.0068  1.6488 
//N   0.0000 -1.2105  1.5645
//H1  0.0000 -1.3039  2.5521 
//H2  0.0000 -2.0271  1.0012
//Confirmed visually
void make_camide_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if ( (mode&2))
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"cmdg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='C';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='O';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='N';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id),   str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id),   str->edges[(*edge_id)+1].type='2';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id),   str->edges[(*edge_id)+2].type='1'; 
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+4].type='1';
  }
//Stage II. Resolve geometry if needed
if ( (mode&4))
  {
  //Setup default methyl coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.0000, str->r[(*atom_id)+0].k=+1.0000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+1.0068, str->r[(*atom_id)+1].k=+1.6488;
  str->r[(*atom_id)+2].i=+0.0000, str->r[(*atom_id)+2].j=-1.2105, str->r[(*atom_id)+2].k=+1.5645;
  str->r[(*atom_id)+3].i=+0.0000, str->r[(*atom_id)+3].j=-1.3039, str->r[(*atom_id)+3].k=+2.5521;
  str->r[(*atom_id)+4].i=+0.0000, str->r[(*atom_id)+4].j=-2.0271, str->r[(*atom_id)+4].k=+1.0012;
  //Construct gag geometry
  _construct_gag_geometry(5,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_CARBON));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"cmdg");
  if ( (atom_id)) (*atom_id)=5;
  if ( (edge_id)) (*edge_id)=5;
  }
}//Exit
//This function makes methyl gag insted of an atom0->atom1 :: atom0->CH3
//It adds 4 atoms and 4 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-C = 1.52679A, C-H = 1.08579A, X(C)-C-H = 111.217Deg, z-dih = PI*2/3, H-C-H = 107.671Deg
//X   0.0000  0.0000  0.0000
//C   0.0000  0.0000  1.0000 (1.5268)
//H1  0.0000  1.0122  1.3930
//H2 -0.8766 -0.5061  1.3930 
//H3  0.8766 -0.5061  1.3930
void make_methyl_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if ( (mode&2))
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"metg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='C';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id), str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id), str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id), str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id), str->edges[(*edge_id)+3].type='1';
  }
//Stage II. Resolve geometry if needed
if ( (mode&4))
  {
  //Setup default methyl coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.0000, str->r[(*atom_id)+0].k=+1.0000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+1.0122, str->r[(*atom_id)+1].k=+1.3930;
  str->r[(*atom_id)+2].i=-0.8766, str->r[(*atom_id)+2].j=-0.5061, str->r[(*atom_id)+2].k=+1.3930;
  str->r[(*atom_id)+3].i=+0.8766, str->r[(*atom_id)+3].j=-0.5061, str->r[(*atom_id)+3].k=+1.3930;
  //Construct gag geometry
  _construct_gag_geometry(4,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_CARBON));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"metg");
  if ( (atom_id)) (*atom_id)=4;
  if ( (edge_id)) (*edge_id)=4;
  }
}//Exit
//This function makes ethene gag insted of an atom0->atom1 :: R->CH=CH2
//It adds 5 atoms and 5 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-C = 1.51213A, C-H = 1.08782A, X(C)-C-H = 116.6Deg, X(C)-C-C = 124.8Deg, C-C = 1.31853A, C-H = 1.0775A, C-C-H = 121.6Deg, all in the same plane
//X    0.0000  0.00000  0.00000
//C1   0.0000  0.00000  1.00000 (1.51213)
//H1   0.0000 -0.97268  1.48708
//C2   0.0000  1.08271  1.75250
//H21  0.0000  1.02256  2.82832 
//H22  0.0000  2.07009  1.32112
//Confirmed visually
void make_ethene_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if ( (mode&2))
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"ethg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='C';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='C';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='2';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+4].type='1';
  }
//Stage II. Resolve geometry if needed
if ( (mode&4))
  {
  //Setup default methyl coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=-0.97268, str->r[(*atom_id)+1].k=+1.48708;
  str->r[(*atom_id)+2].i=+0.0000, str->r[(*atom_id)+2].j=+1.08271, str->r[(*atom_id)+2].k=+1.75250;
  str->r[(*atom_id)+3].i=+0.0000, str->r[(*atom_id)+3].j=+1.02256, str->r[(*atom_id)+3].k=+2.82832;
  str->r[(*atom_id)+4].i=+0.0000, str->r[(*atom_id)+4].j=+2.07009, str->r[(*atom_id)+4].k=+1.32112;
  //Construct gag geometry
  _construct_gag_geometry(5,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_CARBON));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"ethg");
  if ( (atom_id)) (*atom_id)=5;
  if ( (edge_id)) (*edge_id)=5;
  }
}//Exit

// N I T R O G E N    G A G S
//This function makes charged amine gag insted of an atom0->atom1 :: atom0->NH3(+)
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: CH3-NH3 = 1.50585A, N-H = 1.01021A, X(CH3)-N-H = 111.472Deg, z-dih = PI*2/3, H-N-H = 107.398Deg
//X   0.0000   0.0000  0.0000
//N   0.0000   0.0000  1.0000 (1.50585)
//H1  0.0000   0.9401  1.3698
//H2  0.81415 -0.47005 1.3698 
//H3 -0.81415 -0.47005 1.3698 
void make_camine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"camg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='N';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+3].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  {
  //Setup default amine coordinates
  str->r[(*atom_id)+0].i=+0.0000,  str->r[(*atom_id)+0].j=+0.0000,  str->r[(*atom_id)+0].k=+1.0000;
  str->r[(*atom_id)+1].i=+0.0000,  str->r[(*atom_id)+1].j=+0.9401,  str->r[(*atom_id)+1].k=+1.3698;
  str->r[(*atom_id)+2].i=-0.81415, str->r[(*atom_id)+2].j=-0.47005, str->r[(*atom_id)+2].k=+1.3698;
  str->r[(*atom_id)+3].i=+0.81415, str->r[(*atom_id)+3].j=-0.47005, str->r[(*atom_id)+3].k=+1.3698;
  //Construct gag geometry
  _construct_gag_geometry(4,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_NITROGEN));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"camg");
  if ( (atom_id)) (*atom_id)=4;
  if ( (edge_id)) (*edge_id)=4;
  }
}//Exit
//This function makes charged nguanidine gag insted of an atom0->atom1 :: atom0->NH~(HN=C(+)~NH2)
//It adds 9 atoms and 9 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-NH = 1.46284A, N-H = 0.996A, X(CH3)-N-H = 117.40Deg, N-C = 1.31841A, C(X)-N-C = 125.17Deg, C-N = 1.32676A, N-C-N = 120.0Deg, C-N-H = 121.60Deg, all in one plane.
//X    0.00000  0.00000  0.00000
//N0   0.00000  0.00000  1.00000 (1.46284)
//H0   0.00000 -0.88426  1.45836
//C    0.00000  1.07773  1.75941 
//N1   0.00000  0.95817  3.08077 
//H11  0.00000  0.06628  3.52410 
//H12  0.00000  1.75601  2.48456 
//N2   0.00000  2.28184  2.31655
//H21  0.00000  3.11172  2.86730 
//H22  0.00000  2.39925  3.30561 
void make_nguanidine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"ngug");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='N';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='C';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='N';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->a[(*atom_id)+5]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+5][0]=str->anames[(*atom_id)+5][1]=str->anames[(*atom_id)+5][2]=' ', str->anames[(*atom_id)+5][3]='H';
  str->a[(*atom_id)+6]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+6][0]=str->anames[(*atom_id)+6][1]=str->anames[(*atom_id)+6][2]=' ', str->anames[(*atom_id)+6][3]='N';
  str->a[(*atom_id)+7]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+7][0]=str->anames[(*atom_id)+7][1]=str->anames[(*atom_id)+7][2]=' ', str->anames[(*atom_id)+7][3]='H';
  str->a[(*atom_id)+8]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+8][0]=str->anames[(*atom_id)+8][1]=str->anames[(*atom_id)+8][2]=' ', str->anames[(*atom_id)+8][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='a';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+3].type='2';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+3, str->edges[(*edge_id)+4].type='1';
  str->edges[(*edge_id)+5].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+3, str->edges[(*edge_id)+5].type='1';
  str->edges[(*edge_id)+6].vertice[0]=(*atom_id)+6, str->edges[(*edge_id)+6].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+6].type='a';
  str->edges[(*edge_id)+7].vertice[0]=(*atom_id)+7, str->edges[(*edge_id)+7].vertice[1]=(*atom_id)+6, str->edges[(*edge_id)+7].type='1';
  str->edges[(*edge_id)+8].vertice[0]=(*atom_id)+8, str->edges[(*edge_id)+8].vertice[1]=(*atom_id)+6, str->edges[(*edge_id)+8].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  {
  //Setup default amine coordinates
  str->r[(*atom_id)+0].i=+0.0000,  str->r[(*atom_id)+0].j=+0.00000,  str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.0000,  str->r[(*atom_id)+1].j=-0.88426,  str->r[(*atom_id)+1].k=+1.45836;
  str->r[(*atom_id)+2].i=+0.0000,  str->r[(*atom_id)+2].j=+1.07773,  str->r[(*atom_id)+2].k=+1.75941;
  str->r[(*atom_id)+3].i=+0.0000,  str->r[(*atom_id)+3].j=+0.95817,  str->r[(*atom_id)+3].k=+3.08077;
  str->r[(*atom_id)+4].i=+0.0000,  str->r[(*atom_id)+4].j=+0.06628,  str->r[(*atom_id)+4].k=+3.52410;
  str->r[(*atom_id)+5].i=+0.0000,  str->r[(*atom_id)+5].j=+1.75601,  str->r[(*atom_id)+5].k=+2.48456;
  str->r[(*atom_id)+6].i=+0.0000,  str->r[(*atom_id)+6].j=+2.28184,  str->r[(*atom_id)+6].k=+2.31655;
  str->r[(*atom_id)+7].i=+0.0000,  str->r[(*atom_id)+7].j=+3.11172,  str->r[(*atom_id)+7].k=+2.86730;
  str->r[(*atom_id)+8].i=+0.0000,  str->r[(*atom_id)+8].j=+2.39925,  str->r[(*atom_id)+8].k=+3.30561;
  //Construct gag geometry
  _construct_gag_geometry(9,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_NITROGEN));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"ngug");
  if ( (atom_id)) (*atom_id)=9;
  if ( (edge_id)) (*edge_id)=9;
  }
}//Exit
//This function makes amide gag insted of an atom0->atom1 :: atom0->NH2
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(CO)-NH2 = 1.34798A, N-H = 0.993921A, X(CO)-N-H = 119.053Deg, H-N-H = 119.336Deg
//X   0.0000  0.0000  0.0000
//N   0.0000  0.0000  1.0000 (1.3480)
//H1  0.0000  0.8689  1.4827
//H2  0.0000 -0.8689  1.4827 
void make_pamine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"pamg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='N';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='1';
  }
//Stage II. Resolve geometry if needed 
if (mode&4)
  {
  //Setup default amide coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.0000, str->r[(*atom_id)+0].k=+1.0000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+0.8689, str->r[(*atom_id)+1].k=+1.4827;
  str->r[(*atom_id)+2].i=+0.0000, str->r[(*atom_id)+2].j=-0.8689, str->r[(*atom_id)+2].k=+1.4827;
  //Construct gag geometry
  _construct_gag_geometry(3,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_NITROGEN));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"pamg");
  if ( (atom_id)) (*atom_id)=3;
  if ( (edge_id)) (*edge_id)=3;
  }
}//Exit
//This function makes amide gag insted of an atom0->atom1 :: atom0->N=CH2
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-NH2 = 1.45565A, N-C = 1.2809A, X(C)-N-C = 122.47Deg, C-H = 1.0875A, N-C-H = 121.65Deg, all in one plane
//X   0.0000  0.00000  0.00000
//N   0.0000  0.00000  1.00000 (1.45565)
//C   0.0000  1.08030  1.68823
//H1  0.0000  1.06417  2.77561
//H2  0.0000  2.05898  2.16240 
void make_nimine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"pimg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='N';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='2';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+3].type='1';
  }
//Stage II. Resolve geometry if needed 
if (mode&4)
  {
  //Setup default amide coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+1.08030, str->r[(*atom_id)+1].k=+1.68823;
  str->r[(*atom_id)+2].i=+0.0000, str->r[(*atom_id)+2].j=+1.06417, str->r[(*atom_id)+2].k=+2.77561;
  str->r[(*atom_id)+3].i=+0.0000, str->r[(*atom_id)+3].j=+2.05898, str->r[(*atom_id)+3].k=+2.16240;
  //Construct gag geometry
  _construct_gag_geometry(4,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_NITROGEN));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"pimg");
  if ( (atom_id)) (*atom_id)=4;
  if ( (edge_id)) (*edge_id)=4;
  }
}//Exit
//This function makes namide gag insted of an atom0->atom1 :: atom0->NH-(C=O)H
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(C)-N = 1.44772A, N-H = 0.9931A, X(C)-N-H = 119.56Deg, X(C)-N-C = 121.76Deg, H-N-C = 118.68Deg, N-C = 1.34474A, C=O = 1.19569A, C-H = 1.09128A, N-C-O = 124.76Deg, N-C-H = 113.0Deg, O-C-H = 122.24Deg, H-N-C-O = 180.0Deg
//X   0.0000  0.00000  0.0000
//N   0.0000  0.00000  1.0000  (1.44772)
//H   0.0000  0.86384  1.48993
//C   0.0000 -1.14338  1.70782
//O   0.0000 -2.12570  2.38953
//H   0.0000 -1.13885  2.13422 
void make_namide_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"nadg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='N';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='C';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='O';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type=(int)'1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type=(int)'1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type=(int)'m';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+3].type=(int)'2';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+4].type=(int)'1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  {
  //Setup default keton coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+0.86384, str->r[(*atom_id)+1].k=+1.48993;
  str->r[(*atom_id)+2].i=+0.00000, str->r[(*atom_id)+2].j=-1.14338, str->r[(*atom_id)+2].k=+1.70782;
  str->r[(*atom_id)+3].i=+0.00000, str->r[(*atom_id)+3].j=-2.12570, str->r[(*atom_id)+3].k=+2.38953;
  str->r[(*atom_id)+4].i=+0.00000, str->r[(*atom_id)+4].j=-1.13885, str->r[(*atom_id)+4].k=+2.13422;
  //Construct gag geometry
  _construct_gag_geometry(5,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_NITROGEN));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"nadg");
  if ( (atom_id)) (*atom_id)=5;
  if ( (edge_id)) (*edge_id)=5;
  }
}//Exit
//This function makes nitro gag insted of an atom0->atom1 :: atom0-> NO2
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, all symmetrical on OZ, X(C)-N = 1.47864A, N-O = 1.19220A, X(C)-N-O = 117.1Deg, O-N-O = 125.8Deg
//X   0.0000  0.00000  0.0000
//N   0.0000  0.00000  1.0000  (1.47864)
//O   0.0000  1.06131  1.5431
//O   0.0000 -1.06131  1.5431
void make_nitro_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"no2g");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_NITROGEN;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='N';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='O';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='O';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='a';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='2';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  {
  //Setup default keton coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.0000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.06131, str->r[(*atom_id)+1].k=+1.5431;
  str->r[(*atom_id)+2].i=+0.00000, str->r[(*atom_id)+2].j=-1.06131, str->r[(*atom_id)+2].k=+1.5431;
  //Construct gag geometry
  _construct_gag_geometry(3,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_NITROGEN));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"no2g");
  if ( (atom_id)) (*atom_id)=3;
  if ( (edge_id)) (*edge_id)=3;
  }
}//Exit

// O X Y G E N    G A G S 
//This function makes hydroxyl gag insted of an atom0->atom1 :: atom0->OH
//It adds 2 atoms and 2 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(C)-O = 1.39965A, O-H = 0.94630A, X(C)-O-H = 109.45Deg
//X   0.0000  0.0000  0.0000
//O   0.0000  0.0000  1.0000 (1.39965)
//H   0.0000  0.8923  1.3151
void make_hydroxyl_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"hxyg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='O';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.0000, str->r[(*atom_id)+0].k=+1.0000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+0.8923, str->r[(*atom_id)+1].k=+1.3151;
  //Construct gag geometry
  _construct_gag_geometry(2,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_OXYGEN));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"hxyg");
  if ( (atom_id)) (*atom_id)=2;
  if ( (edge_id)) (*edge_id)=2;
  }
}//Exit
//This function makes methoxy gag insted of an atom0->atom1 :: atom0->O-CH3
//It adds 2 atoms and 2 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all but H2&H3 in yz plane, X(C)-O = 1.3922A, O-C = 1.3922A, X(C)-O-C = 115.70Deg, C-H = 1.0849A, O-C-H = 109.7Deg, O->C__L_H=120Deg.
//X   0.00000  0.00000  0.00000
//O   0.00000  0.00000  1.00000 (1.3922)
//C   0.00000  1.25448  1.60374
//H1  0.00000  2.02696  0.84197
//H2  0.20185  2.17857  1.91914
//H3 -0.20185  2.17857  1.91914
void make_methoxy_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"omeg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='O';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+4].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.25448, str->r[(*atom_id)+1].k=+1.60374;
  str->r[(*atom_id)+2].i=+0.00000, str->r[(*atom_id)+2].j=+2.02696, str->r[(*atom_id)+2].k=+0.84197;
  str->r[(*atom_id)+3].i=+0.20185, str->r[(*atom_id)+3].j=+2.17857, str->r[(*atom_id)+3].k=+1.91914;
  str->r[(*atom_id)+4].i=-0.20185, str->r[(*atom_id)+4].j=+2.17857, str->r[(*atom_id)+4].k=+1.91914;
  //Construct gag geometry
  _construct_gag_geometry(5,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_OXYGEN));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"omeg");
  if ( (atom_id)) (*atom_id)=5;
  if ( (edge_id)) (*edge_id)=5;
  }
}//Exit

// S I L I C O N     G A G S 
//This function makes silane gag insted of an atom0->atom1 :: atom0->Si-H3
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-Si = 1.88831, Si-H = 1.08615, X(C)-Si-H = 111.12Deg, H-Si-H = 107.78Deg, H...X(C)-Si..H = 2*PI/3
//X    0.00000  0.00000  0.00000
//Si   0.00000  0.00000  1.00000 (1.88831)
//H1   0.00000  1.01319  1.39136
//H2  -0.87745 -0.50660  1.39136
//H3   0.87745 -0.50660  1.39136
void make_silane_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"slag");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_SILICON;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=' ', str->anames[(*atom_id)+0][2]='S', str->anames[(*atom_id)+0][3]='i';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+3].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.01319, str->r[(*atom_id)+1].k=+1.39136;
  str->r[(*atom_id)+2].i=-0.87745, str->r[(*atom_id)+2].j=-0.50660, str->r[(*atom_id)+2].k=+1.39136;
  str->r[(*atom_id)+3].i=+0.87745, str->r[(*atom_id)+3].j=-0.50660, str->r[(*atom_id)+3].k=+1.39136;
  //Construct gag geometry
  _construct_gag_geometry(4,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_SILICON));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"slag");
  if ( (atom_id)) (*atom_id)=4;
  if ( (edge_id)) (*edge_id)=4;
  }
}//Exit
//This function makes silyl gag insted of an atom0->atom1 :: atom0->Si-(CH3)3
//It adds 13 atoms and 13 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-Si = 1.89364, Si-C = 1.89364, X(C)-Si-C = 109.47Deg, C-H = 1.08748A, Si-C-H = 111.46Deg, H-C-H = 107.41Deg, C...X(C)-Si..C = 2*PI/3, H...Si-C...H = 2*PI/3
//X     0.00000  0.00000  0.00000
//Si    0.00000  0.00000  1.00000 (1.89364)
//C1    0.00000  1.78535  1.63118
//C2   -1.54616 -0.89268  1.63118
//C3    1.54616 -0.89268  1.63118
//H11   0.00000  2.87232  1.59778
//H12  -0.28649  2.70693  2.13242
//H13   0.28649 -2.70693  2.13242
//H21  -2.48750 -1.43616  1.59778
//H22  -2.20102 -1.60157  2.13242
//H23  -2.48750 -1.10535  2.13242
//H31   2.48750 -1.43616  1.59778
//H32   2.48750 -1.10535  2.13242
//H33   2.20102 -1.60157  2.13242
void make_silyl_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"slyg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_SILICON;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=' ', str->anames[(*atom_id)+0][2]='S', str->anames[(*atom_id)+0][3]='i';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='C';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='C';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->a[(*atom_id)+5]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+5][0]=str->anames[(*atom_id)+5][1]=str->anames[(*atom_id)+5][2]=' ', str->anames[(*atom_id)+5][3]='H';
  str->a[(*atom_id)+6]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+6][0]=str->anames[(*atom_id)+6][1]=str->anames[(*atom_id)+6][2]=' ', str->anames[(*atom_id)+6][3]='H';
  str->a[(*atom_id)+7]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+7][0]=str->anames[(*atom_id)+7][1]=str->anames[(*atom_id)+7][2]=' ', str->anames[(*atom_id)+7][3]='H';
  str->a[(*atom_id)+8]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+8][0]=str->anames[(*atom_id)+8][1]=str->anames[(*atom_id)+8][2]=' ', str->anames[(*atom_id)+8][3]='H';
  str->a[(*atom_id)+9]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+9][0]=str->anames[(*atom_id)+9][1]=str->anames[(*atom_id)+9][2]=' ', str->anames[(*atom_id)+9][3]='H';
  str->a[(*atom_id)+10]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+10][0]=str->anames[(*atom_id)+10][1]=str->anames[(*atom_id)+10][2]=' ', str->anames[(*atom_id)+10][3]='H';
  str->a[(*atom_id)+11]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+11][0]=str->anames[(*atom_id)+11][1]=str->anames[(*atom_id)+11][2]=' ', str->anames[(*atom_id)+11][3]='H';
  str->a[(*atom_id)+12]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+12][0]=str->anames[(*atom_id)+12][1]=str->anames[(*atom_id)+12][2]=' ', str->anames[(*atom_id)+12][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+4].type='1';
  str->edges[(*edge_id)+5].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+5].type='1';
  str->edges[(*edge_id)+6].vertice[0]=(*atom_id)+6, str->edges[(*edge_id)+6].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+6].type='1';
  str->edges[(*edge_id)+7].vertice[0]=(*atom_id)+7, str->edges[(*edge_id)+7].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+7].type='1';
  str->edges[(*edge_id)+8].vertice[0]=(*atom_id)+8, str->edges[(*edge_id)+8].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+8].type='1';
  str->edges[(*edge_id)+9].vertice[0]=(*atom_id)+9, str->edges[(*edge_id)+9].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+9].type='1';
  str->edges[(*edge_id)+10].vertice[0]=(*atom_id)+10, str->edges[(*edge_id)+10].vertice[1]=(*atom_id)+3, str->edges[(*edge_id)+10].type='1';
  str->edges[(*edge_id)+11].vertice[0]=(*atom_id)+11, str->edges[(*edge_id)+11].vertice[1]=(*atom_id)+3, str->edges[(*edge_id)+11].type='1';
  str->edges[(*edge_id)+12].vertice[0]=(*atom_id)+12, str->edges[(*edge_id)+12].vertice[1]=(*atom_id)+3, str->edges[(*edge_id)+12].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.78535, str->r[(*atom_id)+1].k=+1.63118;
  str->r[(*atom_id)+2].i=-1.54616, str->r[(*atom_id)+2].j=-0.89268, str->r[(*atom_id)+2].k=+1.63118;
  str->r[(*atom_id)+3].i=+1.54616, str->r[(*atom_id)+3].j=-0.89268, str->r[(*atom_id)+3].k=+1.63118;
  str->r[(*atom_id)+4].i=+0.00000, str->r[(*atom_id)+4].j=+2.87232, str->r[(*atom_id)+4].k=+1.59778;
  str->r[(*atom_id)+5].i=-0.28649, str->r[(*atom_id)+5].j=+2.70693, str->r[(*atom_id)+5].k=+2.13242;
  str->r[(*atom_id)+6].i=+0.28649, str->r[(*atom_id)+6].j=-2.70693, str->r[(*atom_id)+6].k=+2.13242;
  str->r[(*atom_id)+7].i=-2.48750, str->r[(*atom_id)+7].j=-1.43616, str->r[(*atom_id)+7].k=+1.59778;
  str->r[(*atom_id)+8].i=-2.20102, str->r[(*atom_id)+8].j=-1.60157, str->r[(*atom_id)+8].k=+2.13242;
  str->r[(*atom_id)+9].i=-2.48750, str->r[(*atom_id)+9].j=-1.10535, str->r[(*atom_id)+9].k=+2.13242;
  str->r[(*atom_id)+10].i=+2.48750, str->r[(*atom_id)+10].j=-1.43616, str->r[(*atom_id)+10].k=+1.59778;
  str->r[(*atom_id)+11].i=+2.48750, str->r[(*atom_id)+11].j=-1.10535, str->r[(*atom_id)+11].k=+2.13242;
  str->r[(*atom_id)+12].i=+2.20102, str->r[(*atom_id)+12].j=-1.60157, str->r[(*atom_id)+12].k=+2.13242;
  //Construct gag geometry
  _construct_gag_geometry(13,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_SILICON));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"slyg");
  if ( (atom_id)) (*atom_id)=13;
  if ( (edge_id)) (*edge_id)=13;
  }
}
// P H O S P H O R     G A G S 
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->P-H2
//It adds 3 atoms and 3 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.86035, P-H = 1.40416, X(C)-P-H = 98.65Deg, H-P-H = 95.09Deg, H...X(C)-P..H = 96.5Deg
//X    0.00000  0.00000  0.00000
//P    0.00000  0.00000  1.00000 (1.86035)
//H1   0.00000  1.38819  1.21118
//H2  -1.37926 -0.15715  1.21118
void make_hphosphine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"hphg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_PHOSPHOR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='P';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.38819, str->r[(*atom_id)+1].k=+1.21118;
  str->r[(*atom_id)+2].i=-1.37926, str->r[(*atom_id)+2].j=-0.15715, str->r[(*atom_id)+2].k=+1.21118;
  //Construct gag geometry
  _construct_gag_geometry(3,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_PHOSPHOR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"hphg");
  if ( (atom_id)) (*atom_id)=3;
  if ( (edge_id)) (*edge_id)=3;
  }
}//Exit
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->P-(CH3)2
//It adds 9 atoms and 9 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.85344, P-C = 1.85344, X(C)-P-C = 100.03Deg, C-P-C = 100.03Deg, C...X(C)-P..C = 102.2Deg, C-H = 1.08534A, P-C-H = 109.8Deg, H-C-H = 107.67Deg, H...P->C-H = 2PI/3
//X    0.00000  0.00000  0.00000
//P    0.00000  0.00000  1.00000 (1.85344)
//C1   0.00000  1.82606  1.32297
//C2  -1.78482 -0.38589  1.32297  
//H11  0.00000  2.36594  0.38143
//H12 -1.63996  3.42353  2.12402
//H13  1.63996  3.42353  2.12402
//H21 -2.31250 -0.49998  0.38143
//H22 -2.99964 -2.32640  2.12402
//H23 -3.69277  0.87945  2.12402
void make_mphosphine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"mphg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_PHOSPHOR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='P';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='C';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->a[(*atom_id)+5]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+5][0]=str->anames[(*atom_id)+5][1]=str->anames[(*atom_id)+5][2]=' ', str->anames[(*atom_id)+5][3]='H';
  str->a[(*atom_id)+6]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+6][0]=str->anames[(*atom_id)+6][1]=str->anames[(*atom_id)+6][2]=' ', str->anames[(*atom_id)+6][3]='H';
  str->a[(*atom_id)+7]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+7][0]=str->anames[(*atom_id)+7][1]=str->anames[(*atom_id)+7][2]=' ', str->anames[(*atom_id)+7][3]='H';
  str->a[(*atom_id)+8]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+8][0]=str->anames[(*atom_id)+8][1]=str->anames[(*atom_id)+8][2]=' ', str->anames[(*atom_id)+8][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+4].type='1';
  str->edges[(*edge_id)+5].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+5].type='1';
  str->edges[(*edge_id)+6].vertice[0]=(*atom_id)+6, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+6].type='1';
  str->edges[(*edge_id)+7].vertice[0]=(*atom_id)+7, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+7].type='1';
  str->edges[(*edge_id)+8].vertice[0]=(*atom_id)+8, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+8].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.82606, str->r[(*atom_id)+1].k=+1.32297;
  str->r[(*atom_id)+2].i=-1.78482, str->r[(*atom_id)+2].j=-0.38589, str->r[(*atom_id)+2].k=+1.32297;
  str->r[(*atom_id)+3].i=+0.00000, str->r[(*atom_id)+3].j=+2.36594, str->r[(*atom_id)+3].k=+0.38143;
  str->r[(*atom_id)+4].i=-1.63996, str->r[(*atom_id)+4].j=+3.42353, str->r[(*atom_id)+4].k=+2.12402;
  str->r[(*atom_id)+5].i=+1.63996, str->r[(*atom_id)+5].j=+3.42353, str->r[(*atom_id)+5].k=+2.12402;
  str->r[(*atom_id)+6].i=-2.31250, str->r[(*atom_id)+6].j=-0.49998, str->r[(*atom_id)+6].k=+0.38143;
  str->r[(*atom_id)+7].i=-2.99964, str->r[(*atom_id)+7].j=-2.32640, str->r[(*atom_id)+7].k=+2.12402;
  str->r[(*atom_id)+8].i=-3.69277, str->r[(*atom_id)+8].j=+0.87945, str->r[(*atom_id)+8].k=+2.12402;
  //Construct gag geometry
  _construct_gag_geometry(9,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_PHOSPHOR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"mphg");
  if ( (atom_id)) (*atom_id)=9;
  if ( (edge_id)) (*edge_id)=9;
  }
}//Exit
//This function makes phosphane gag insted of an atom0->atom1 :: atom0->P=CH2
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.8624A, P-C = 1.64612A, X(C)-P-C = 103.825Deg, C-H = 1.07666A, P-C-H = 122.17Deg, H-C-H = 115.66Deg, H..(P-C)..H = PI.
//X   0.00000  0.00000  0.00000
//P   0.00000  0.00000  1.00000 (1.89757)
//C   0.00000  1.59843  1.39335
//H1  0.00000  1.93730  2.41529
//H2  0.00000  2.37285  2.14133
void make_phosphane_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"ppng");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_PHOSPHOR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='P';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='2';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+3].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+1.59843, str->r[(*atom_id)+1].k=+1.39335;
  str->r[(*atom_id)+2].i=+0.0000, str->r[(*atom_id)+2].j=+1.93730, str->r[(*atom_id)+2].k=+2.41529;
  str->r[(*atom_id)+3].i=+0.0000, str->r[(*atom_id)+3].j=+2.37285, str->r[(*atom_id)+3].k=+2.14133;
  //Construct gag geometry
  _construct_gag_geometry(4,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_PHOSPHOR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"ppng");
  if ( (atom_id)) (*atom_id)=4;
  if ( (edge_id)) (*edge_id)=4;
  }
}//Exit
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->(P=O)-O(-)2
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.89757A, P-O = 1.5190A, X(C)-P-O = 103.15Deg, O...X(C)-P..O = 2PI/3
//X   0.00000  0.00000  0.00000
//P   0.00000  0.00000  1.00000 (1.89757)
//O1  0.00000  1.47917  1.34557
//O2 -1.28100 -0.73958  1.34557
//O3  1.28100 -0.73958  1.34557
void make_phosphate2_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"pp2g");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_PHOSPHOR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='P';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='O';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='O';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='O';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='2';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='a';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+3].type='a';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.0000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.0000, str->r[(*atom_id)+1].j=+1.47917, str->r[(*atom_id)+1].k=+1.34557;
  str->r[(*atom_id)+2].i=-1.2810, str->r[(*atom_id)+2].j=-0.73958, str->r[(*atom_id)+2].k=+1.34557;
  str->r[(*atom_id)+3].i=+1.2810, str->r[(*atom_id)+3].j=-0.73958, str->r[(*atom_id)+3].k=+1.34557;
  //Construct gag geometry
  _construct_gag_geometry(4,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_PHOSPHOR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"pp2g");
  if ( (atom_id)) (*atom_id)=4;
  if ( (edge_id)) (*edge_id)=4;
  }
}//Exit
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->(O=P~O(-))-CH3
//It adds 7 atoms and 7 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.84644A, X(C)-P-C = 101.32Deg, C-H = 1.087324A, P-C-H = 111.60Deg, P-O = 1.48887A, X(C)-P-O = 107.88Deg, C/O...X(C)-P..O = 2PI/3
//X    0.00000  0.00000  0.00000
//P    0.00000  0.00000  1.00000 (1.81915)
//C1   0.00000  1.81052  1.36243
//H11  0.00000  2.00456  2.43230
//H12  0.87552  2.30222  0.94535 
//H13 -0.87552  2.30222  0.94535
//O1   1.22712 -0.70848  1.45707 
//O2  -1.22712 -0.70848  1.45707
void make_phosphate1_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"pp1g");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_PHOSPHOR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='P';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->a[(*atom_id)+5]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+5][0]=str->anames[(*atom_id)+5][1]=str->anames[(*atom_id)+5][2]=' ', str->anames[(*atom_id)+5][3]='O';
  str->a[(*atom_id)+6]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+6][0]=str->anames[(*atom_id)+6][1]=str->anames[(*atom_id)+6][2]=' ', str->anames[(*atom_id)+6][3]='O';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+4].type='1';
  str->edges[(*edge_id)+5].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+5].type='2';
  str->edges[(*edge_id)+6].vertice[0]=(*atom_id)+6, str->edges[(*edge_id)+6].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+6].type='a';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.81052, str->r[(*atom_id)+1].k=+1.36243;
  str->r[(*atom_id)+2].i=+0.00000, str->r[(*atom_id)+2].j=+2.00456, str->r[(*atom_id)+2].k=+2.43230;
  str->r[(*atom_id)+3].i=+0.87552, str->r[(*atom_id)+3].j=+2.30222, str->r[(*atom_id)+3].k=+0.94535;
  str->r[(*atom_id)+4].i=-0.87552, str->r[(*atom_id)+4].j=+2.30222, str->r[(*atom_id)+4].k=+0.94535;
  str->r[(*atom_id)+5].i=+1.22712, str->r[(*atom_id)+5].j=-0.70848, str->r[(*atom_id)+5].k=+1.45707;
  str->r[(*atom_id)+6].i=-1.22712, str->r[(*atom_id)+6].j=-0.70848, str->r[(*atom_id)+6].k=+1.45707;
  //Construct gag geometry
  _construct_gag_geometry(7,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_PHOSPHOR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"pp1g");
  if ( (atom_id)) (*atom_id)=7;
  if ( (edge_id)) (*edge_id)=7;
  }
}//Exit
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->(P=O)-(CH3)2
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.81915A, X(C)-P-C = 104.80Deg, C-H = 1.08524A, P-C-H = 109.44Deg, H..P-C..H = 2PI/3, P-O = 1.47371A, X(C)-P-O = 113.80Deg, C...X(C)-P..O = 2PI/3
//X    0.00000  0.00000  0.00000
//P    0.00000  0.00000  1.00000 (1.81915)
//C1   0.00000  1.75880  1.48513
//H11  0.00000  1.92393  2.56681
//H12  0.88627  2.31606  1.08268 
//H13 -0.88627  2.31606  1.08268
//C2   1.52316 -0.87940  1.48513 
//H21  1.66617  0.96197  2.56681
//H22  1.56263 -1.92556  1.08268
//H23  2.44890  0.39050  1.08268 
//O   -1.16774 -0.67419  1.59471
void make_phosphate0_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"pp0g");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_PHOSPHOR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='P';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->a[(*atom_id)+5]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+5][0]=str->anames[(*atom_id)+5][1]=str->anames[(*atom_id)+5][2]=' ', str->anames[(*atom_id)+5][3]='C';
  str->a[(*atom_id)+6]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+6][0]=str->anames[(*atom_id)+6][1]=str->anames[(*atom_id)+6][2]=' ', str->anames[(*atom_id)+6][3]='H';
  str->a[(*atom_id)+7]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+7][0]=str->anames[(*atom_id)+7][1]=str->anames[(*atom_id)+7][2]=' ', str->anames[(*atom_id)+7][3]='H';
  str->a[(*atom_id)+8]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+8][0]=str->anames[(*atom_id)+8][1]=str->anames[(*atom_id)+8][2]=' ', str->anames[(*atom_id)+8][3]='H';
  str->a[(*atom_id)+9]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+9][0]=str->anames[(*atom_id)+9][1]=str->anames[(*atom_id)+9][2]=' ', str->anames[(*atom_id)+9][3]='O';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+4].type='1';
  str->edges[(*edge_id)+5].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+5].type='1';
  str->edges[(*edge_id)+6].vertice[0]=(*atom_id)+6, str->edges[(*edge_id)+6].vertice[1]=(*atom_id)+5, str->edges[(*edge_id)+6].type='1';
  str->edges[(*edge_id)+7].vertice[0]=(*atom_id)+7, str->edges[(*edge_id)+7].vertice[1]=(*atom_id)+5, str->edges[(*edge_id)+7].type='1';
  str->edges[(*edge_id)+8].vertice[0]=(*atom_id)+8, str->edges[(*edge_id)+8].vertice[1]=(*atom_id)+5, str->edges[(*edge_id)+8].type='1';
  str->edges[(*edge_id)+9].vertice[0]=(*atom_id)+9, str->edges[(*edge_id)+9].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+9].type='2';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.75880, str->r[(*atom_id)+1].k=+1.48513;
  str->r[(*atom_id)+2].i=+0.00000, str->r[(*atom_id)+2].j=+1.92393, str->r[(*atom_id)+2].k=+2.56681;
  str->r[(*atom_id)+3].i=+0.88627, str->r[(*atom_id)+3].j=+2.31606, str->r[(*atom_id)+3].k=+1.08268;
  str->r[(*atom_id)+4].i=-0.88627, str->r[(*atom_id)+4].j=+2.31606, str->r[(*atom_id)+4].k=+1.08268;
  str->r[(*atom_id)+5].i=+1.52316, str->r[(*atom_id)+5].j=-0.87940, str->r[(*atom_id)+5].k=+1.48513;
  str->r[(*atom_id)+6].i=+1.66617, str->r[(*atom_id)+6].j=+0.96197, str->r[(*atom_id)+6].k=+2.56681;
  str->r[(*atom_id)+7].i=+1.56263, str->r[(*atom_id)+7].j=-1.92556, str->r[(*atom_id)+7].k=+1.08268;
  str->r[(*atom_id)+8].i=+2.44890, str->r[(*atom_id)+8].j=+0.39050, str->r[(*atom_id)+8].k=+1.08268;
  str->r[(*atom_id)+9].i=-1.16774, str->r[(*atom_id)+9].j=-0.67419, str->r[(*atom_id)+9].k=+1.59471;
  //Construct gag geometry
  _construct_gag_geometry(10,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_PHOSPHOR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"pp0g");
  if ( (atom_id)) (*atom_id)=10;
  if ( (edge_id)) (*edge_id)=10;
  }
}//Exit


// S U L F U R     G A G S 
//This function makes sulfide gag insted of an atom0->atom1 :: atom0->S-CH3
//It adds 5 atoms and 5 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.8087A, S-C = 1.8087A, X(C)-S-C = 100.0Deg, C-H = 1.08304A, S-C-H = 111.16Deg, H-C-H = 109.5Deg, H...S-C..H = 2PI/3
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.8087)
//C   0.00000  1.78122  1.31408
//H1  0.00000  2.31454  0.37145
//H2 -1.59915  3.26166  2.09345
//H3  1.59915  3.26166  2.09345
void make_sulfide_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"sidg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_SULFUR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='S';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+4].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.78122, str->r[(*atom_id)+1].k=+1.31408;
  str->r[(*atom_id)+2].i=+0.00000, str->r[(*atom_id)+2].j=+2.31454, str->r[(*atom_id)+2].k=+0.37145;
  str->r[(*atom_id)+3].i=-1.59915, str->r[(*atom_id)+3].j=+3.26166, str->r[(*atom_id)+3].k=+2.09345;
  str->r[(*atom_id)+4].i= 1.59915, str->r[(*atom_id)+4].j=+3.26166, str->r[(*atom_id)+4].k=+2.09345;
  //Construct gag geometry
  _construct_gag_geometry(5,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_SULFUR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"sidg");
  if ( (atom_id)) (*atom_id)=5;
  if ( (edge_id)) (*edge_id)=5;
 }
}//Exit
//This function makes thiol gag insted of an atom0->atom1 :: atom0->S-H
//It adds 2 atoms and 2 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(C)-S = 1.81768A, S-H = 1.32663A, X(C)-S-H = 97.9Deg
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.81768)
//H   0.00000  1.31404  1.18234
void make_thiol_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"thig");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_SULFUR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='S';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.31404, str->r[(*atom_id)+1].k=+1.18234;
  //Construct gag geometry
  _construct_gag_geometry(2,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_SULFUR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"thig");
  if ( (atom_id)) (*atom_id)=2;
  if ( (edge_id)) (*edge_id)=2;
  }
}//Exit
//This function makes sulfite gag insted of an atom0->atom1 :: atom0->(S=O)-O(-)
//It adds 3 atoms and 3 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.82014A, S-O = 1.4944A, X(C)-S-O = 101.0Deg, O-S-O = 113.82Deg, O...X(C)-S..O = 117.1
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.82014)
//O1  0.00000  1.46694  1.28514
//O2 -1.30589 -0.66826  1.28514
void make_isulfite_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"isig");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_SULFUR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='S';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='O';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='O';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='2';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='2';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.46694, str->r[(*atom_id)+1].k=+1.28514;
  str->r[(*atom_id)+2].i=-1.30589, str->r[(*atom_id)+2].j=-0.66826, str->r[(*atom_id)+2].k=+1.28514;
  //Construct gag geometry
  _construct_gag_geometry(3,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_SULFUR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"isig");
  if ( (atom_id)) (*atom_id)=3;
  if ( (edge_id)) (*edge_id)=3;
  }
}//Exit
//This function makes methyl sulfite gag insted of an atom0->atom1 :: atom0->(S=O)-CH3
//It adds 6 atoms and 6 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.79647A, C-H = 1.08300A, S-C-H = 109.68Deg, X(C)-S-C = 97.66Deg, S-O = 1.48524, X(C)-S-O = 106.72Deg, C...X(C)-S..O = 110.9Deg
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.82014)
//C   0.00000  1.78044  1.23946
//H1  0.00000  2.00598  2.29871
//H2  0.88312  2.20987  0.78275
//H3 -0.88312  2.20987  0.78275
//O2 -1.32886 -0.50744  1.42730
void make_msulfite_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"msig");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_SULFUR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='S';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->a[(*atom_id)+5]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+5][0]=str->anames[(*atom_id)+5][1]=str->anames[(*atom_id)+5][2]=' ', str->anames[(*atom_id)+5][3]='O';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+4].type='1';
  str->edges[(*edge_id)+5].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+5].type='2';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.78044, str->r[(*atom_id)+1].k=+1.23946;
  str->r[(*atom_id)+2].i=+0.00000, str->r[(*atom_id)+2].j=+2.00598, str->r[(*atom_id)+2].k=+2.29871;
  str->r[(*atom_id)+3].i=+0.88312, str->r[(*atom_id)+3].j=+2.20987, str->r[(*atom_id)+3].k=+0.78275;
  str->r[(*atom_id)+4].i=-0.88312, str->r[(*atom_id)+4].j=+2.20987, str->r[(*atom_id)+4].k=+0.78275;
  str->r[(*atom_id)+5].i=-1.32886, str->r[(*atom_id)+5].j=-0.50744, str->r[(*atom_id)+5].k=+1.42730;
  //Construct gag geometry
  _construct_gag_geometry(6,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_SULFUR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"msig");
  if ( (atom_id)) (*atom_id)=6;
  if ( (edge_id)) (*edge_id)=6;
  }
}//Exit
//This function makes sulfate gag insted of an atom0->atom1 :: atom0->(O=S=O)-O(-)
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.78722A, S-O = 1.45473A, X(C)-S-O = 104.46Deg, O...X(C)-S..O = 2PI/3
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.78722)
//O1  0.00000  1.40865  1.36325
//O2 -1.21992 -0.70432  1.36325
//O3  1.21992 -0.70432  1.36325
void make_isulfate_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"isag");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_SULFUR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='S';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='O';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='O';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='O';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='2';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+2].type='2';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+3].type='a';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.40865, str->r[(*atom_id)+1].k=+1.36325;
  str->r[(*atom_id)+2].i=-1.21992, str->r[(*atom_id)+2].j=-0.70432, str->r[(*atom_id)+2].k=+1.36325;
  str->r[(*atom_id)+3].i=+1.21992, str->r[(*atom_id)+3].j=-0.70432, str->r[(*atom_id)+3].k=+1.36325;
  //Construct gag geometry
  _construct_gag_geometry(4,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_SULFUR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"isag");
  if ( (atom_id)) (*atom_id)=4;
  if ( (edge_id)) (*edge_id)=4;
  }
}//Exit
//This function makes methyl-sulfate gag insted of an atom0->atom1 :: atom0->(O=S=O)-CH3
//It adds 7 atoms and 7 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.78722A, X(C)-S-C = 104.32Deg, C-H = 1.08124A, S-C-H = 109.76Deg, H-(S->C)-H = 2PI/3, S-O = 1.45473A, X(C)-S-O = 104.46Deg, O...X(C)-S..O = 2PI/3
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.78722)
//C1  0.00000  1.73060  1.44628
//H1  0.00000  1.83048  2.52289
//H2  0.88124  2.21162  1.04489
//H3 -0.88124  2.21162  1.04489
//O2 -1.21992 -0.70432  1.36325
//O3  1.21992 -0.70432  1.36325
void make_msulfate_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"msag");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_SULFUR;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='S';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='H';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->a[(*atom_id)+5]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+5][0]=str->anames[(*atom_id)+5][1]=str->anames[(*atom_id)+5][2]=' ', str->anames[(*atom_id)+5][3]='O';
  str->a[(*atom_id)+6]=CHEM_ATOM_TYPE_OXYGEN;
  str->anames[(*atom_id)+6][0]=str->anames[(*atom_id)+6][1]=str->anames[(*atom_id)+6][2]=' ', str->anames[(*atom_id)+6][3]='O';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='1';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+2, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+3].type='1';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+4, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+4].type='1';
  str->edges[(*edge_id)+5].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+5].type='2';
  str->edges[(*edge_id)+6].vertice[0]=(*atom_id)+6, str->edges[(*edge_id)+6].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+6].type='2';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.00000, str->r[(*atom_id)+0].j=+0.00000, str->r[(*atom_id)+0].k=+1.00000;
  str->r[(*atom_id)+1].i=+0.00000, str->r[(*atom_id)+1].j=+1.73060, str->r[(*atom_id)+1].k=+1.44628;
  str->r[(*atom_id)+2].i=+0.00000, str->r[(*atom_id)+2].j=+1.83048, str->r[(*atom_id)+2].k=+2.52289;
  str->r[(*atom_id)+3].i=+0.88124, str->r[(*atom_id)+3].j=+2.21162, str->r[(*atom_id)+3].k=+1.04489;
  str->r[(*atom_id)+4].i=-0.88124, str->r[(*atom_id)+4].j=+2.21162, str->r[(*atom_id)+4].k=+1.04489;
  str->r[(*atom_id)+5].i=-1.21992, str->r[(*atom_id)+5].j=-0.70432, str->r[(*atom_id)+5].k=+1.36325;
  str->r[(*atom_id)+6].i=+1.21992, str->r[(*atom_id)+6].j=-0.70432, str->r[(*atom_id)+6].k=+1.36325;
  //Construct gag geometry
  _construct_gag_geometry(7,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_SULFUR));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"msag");
  if ( (atom_id)) (*atom_id)=7;
  if ( (edge_id)) (*edge_id)=7;
  }
}//Exit


// M I S C E L L A N E O U S     G A G S 
//This function makes phenyl gag insted of an atom0->atom1 :: atom0->C6H5
//It adds 11 atoms and 12 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-C = 1.56A, X(C)-C-C(H) = 120.0Deg, C-C = 1.34A, C-H = 1.08A all in the XZ plane
//X   0.000000  0.000000  0.000000
//C1  0.000000  0.000000  1.000000 (1.56)
//C2  1.160474  0.000000  1.670000
//H2  2.095781  0.000000  1.130000
//C3  1.160474  0.000000  3.010000
//H3  2.095781  0.000000  3.550000
//C4  0.000000  0.000000  3.680000
//H4  0.000000  0.000000  4.760000
//C5 -1.160474  0.000000  3.010000
//H5 -2.095781  0.000000  3.550000
//C6 -1.160474  0.000000  1.670000
//H6 -2.095781  0.000000  1.130000
//Confirmed visualy
void make_phenyl_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str)
{
//Stage I. Resolve topology if needed 
if (mode&2)
  {
  str->ress->list[(*ress_id)]=*((unsigned int*)&"pheg");
  str->a[(*atom_id)+0]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+0][0]=str->anames[(*atom_id)+0][1]=str->anames[(*atom_id)+0][2]=' ', str->anames[(*atom_id)+0][3]='C';
  str->a[(*atom_id)+1]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+1][0]=str->anames[(*atom_id)+1][1]=str->anames[(*atom_id)+1][2]=' ', str->anames[(*atom_id)+1][3]='C';
  str->a[(*atom_id)+2]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+2][0]=str->anames[(*atom_id)+2][1]=str->anames[(*atom_id)+2][2]=' ', str->anames[(*atom_id)+2][3]='H';
  str->a[(*atom_id)+3]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+3][0]=str->anames[(*atom_id)+3][1]=str->anames[(*atom_id)+3][2]=' ', str->anames[(*atom_id)+3][3]='C';
  str->a[(*atom_id)+4]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+4][0]=str->anames[(*atom_id)+4][1]=str->anames[(*atom_id)+4][2]=' ', str->anames[(*atom_id)+4][3]='H';
  str->a[(*atom_id)+5]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+5][0]=str->anames[(*atom_id)+5][1]=str->anames[(*atom_id)+5][2]=' ', str->anames[(*atom_id)+5][3]='C';
  str->a[(*atom_id)+6]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+6][0]=str->anames[(*atom_id)+6][1]=str->anames[(*atom_id)+6][2]=' ', str->anames[(*atom_id)+6][3]='H';
  str->a[(*atom_id)+7]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+7][0]=str->anames[(*atom_id)+7][1]=str->anames[(*atom_id)+7][2]=' ', str->anames[(*atom_id)+7][3]='C';
  str->a[(*atom_id)+8]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+8][0]=str->anames[(*atom_id)+8][1]=str->anames[(*atom_id)+8][2]=' ', str->anames[(*atom_id)+8][3]='H';
  str->a[(*atom_id)+9]=CHEM_ATOM_TYPE_CARBON;
  str->anames[(*atom_id)+9][0]=str->anames[(*atom_id)+9][1]=str->anames[(*atom_id)+9][2]=' ', str->anames[(*atom_id)+9][3]='C';
  str->a[(*atom_id)+10]=CHEM_ATOM_TYPE_HYDROGEN;
  str->anames[(*atom_id)+10][0]=str->anames[(*atom_id)+10][1]=str->anames[(*atom_id)+10][2]=' ', str->anames[(*atom_id)+10][3]='H';
  str->edges[(*edge_id)+0].vertice[0]=acceptor_id,  str->edges[(*edge_id)+0].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+0].type='1';
  str->edges[(*edge_id)+1].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+1].vertice[1]=(*atom_id)+0, str->edges[(*edge_id)+1].type='2';
  str->edges[(*edge_id)+2].vertice[0]=(*atom_id)+1, str->edges[(*edge_id)+2].vertice[1]=(*atom_id)+2, str->edges[(*edge_id)+2].type='1';
  str->edges[(*edge_id)+3].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+3].vertice[1]=(*atom_id)+1, str->edges[(*edge_id)+3].type='a';
  str->edges[(*edge_id)+4].vertice[0]=(*atom_id)+3, str->edges[(*edge_id)+4].vertice[1]=(*atom_id)+4, str->edges[(*edge_id)+4].type='1';
  str->edges[(*edge_id)+5].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+5].vertice[1]=(*atom_id)+3, str->edges[(*edge_id)+5].type='2';
  str->edges[(*edge_id)+6].vertice[0]=(*atom_id)+5, str->edges[(*edge_id)+6].vertice[1]=(*atom_id)+6, str->edges[(*edge_id)+6].type='1';
  str->edges[(*edge_id)+7].vertice[0]=(*atom_id)+7, str->edges[(*edge_id)+7].vertice[1]=(*atom_id)+5, str->edges[(*edge_id)+7].type='a';
  str->edges[(*edge_id)+8].vertice[0]=(*atom_id)+7, str->edges[(*edge_id)+8].vertice[1]=(*atom_id)+8, str->edges[(*edge_id)+8].type='1';
  str->edges[(*edge_id)+9].vertice[0]=(*atom_id)+9, str->edges[(*edge_id)+9].vertice[1]=(*atom_id)+7, str->edges[(*edge_id)+9].type='2';
  str->edges[(*edge_id)+10].vertice[0]=(*atom_id)+9, str->edges[(*edge_id)+10].vertice[1]=(*atom_id)+10, str->edges[(*edge_id)+10].type='1';
  str->edges[(*edge_id)+11].vertice[0]=(*atom_id)+0, str->edges[(*edge_id)+11].vertice[1]=(*atom_id)+9, str->edges[(*edge_id)+11].type='a';
  }
//Stage II. Resolve geometry if needed
if (mode&4)
  { 
  //Setup default hydroxyl coordinates
  str->r[(*atom_id)+0].i=+0.000000, str->r[(*atom_id)+0].j=+0.000000, str->r[(*atom_id)+0].k=+1.000000;
  str->r[(*atom_id)+1].i=+1.160474, str->r[(*atom_id)+1].j=+0.000000, str->r[(*atom_id)+1].k=+1.670000;
  str->r[(*atom_id)+2].i=+2.095781, str->r[(*atom_id)+2].j=+0.000000, str->r[(*atom_id)+2].k=+1.130000;
  str->r[(*atom_id)+3].i=+1.160474, str->r[(*atom_id)+3].j=+0.000000, str->r[(*atom_id)+3].k=+3.010000;
  str->r[(*atom_id)+4].i=+2.095781, str->r[(*atom_id)+4].j=+0.000000, str->r[(*atom_id)+4].k=+3.550000;
  str->r[(*atom_id)+5].i=+0.000000, str->r[(*atom_id)+5].j=+0.000000, str->r[(*atom_id)+5].k=+3.680000;
  str->r[(*atom_id)+6].i=+0.000000, str->r[(*atom_id)+6].j=+0.000000, str->r[(*atom_id)+6].k=+4.760000;
  str->r[(*atom_id)+7].i=-1.160474, str->r[(*atom_id)+7].j=+0.000000, str->r[(*atom_id)+7].k=+3.010000;
  str->r[(*atom_id)+8].i=-2.095781, str->r[(*atom_id)+8].j=+0.000000, str->r[(*atom_id)+8].k=+3.550000;
  str->r[(*atom_id)+9].i=-1.160474, str->r[(*atom_id)+9].j=+0.000000, str->r[(*atom_id)+9].k=+1.670000;
  str->r[(*atom_id)+10].i=-2.095781, str->r[(*atom_id)+10].j=+0.000000, str->r[(*atom_id)+10].k=+1.130000;
  //Construct gag geometry
  _construct_gag_geometry(11,&str->r[(*atom_id)],&str->r[acceptor_id],r,get_approximate_single_bond_lenght(str->a[acceptor_id],CHEM_ATOM_TYPE_CARBON));
  }
if ( (mode&1))
  {
  if ( (ress_id)) (*ress_id)=*((unsigned int*)&"pheg");
  if ( (atom_id)) (*atom_id)=11;
  if ( (edge_id)) (*edge_id)=12;
  }
}//Exit


//The gags are named as separate residues  
//The gags are one of the list:
//  function_name           chemistry             gag_name   description
//  make_carboxy_gag()    : R<-[O=C~O(-)]          [cxyg]  (carboxylic group is kept)
//  make_cguanidine_gag() : R<-[NH2~C(+)=NH2]      [cgug]  (guanidine group is kept)
//  make_cyane_gag()      : R<-C_=_N               [cyng]  (cyane group is kept)
//  make_camide_gag()     : R<-(C=O)~NH2           [cmdg]  (C near/in resonance)
//  make_methyl_gag()     : R<-CH3                 [metg]  (generic aliphatic C)
//  make_ethene_gag()     : R<-CH=CH2              [ethg]  (alkene carbon) 
//  make_camine_gag()     : R<-NH3+                [camg]  (aliphatic charged N)
//  make_nguanidine_gag() : R<-XN~[NH2~C(+)=NH2]   [ngug]  (guanidine group is kept)
//  make_pamine_gag()     : R<-NH2                 [pamg]  (N near aromatic if X is sp2 atom)
//  make_nimine_gag()     : R<-N=CH2               [pimg]  (generic sp2 N)
//  make_namide_gag()     : R<-NH~(C=O)H           [nadg]  (N near resonance if X is sp3 atom)
//  make_nitro_gag()      : R<-NO2                 [no2g]  (NO2 group is kept)
//  make_hydroxy_gag()    : R<-OH                  [hxyg]  (hydroxy gag) 
//  make_methoxy_gag()    : R<-O-CH3               [omeg]  (generic oxygen gag)
//  make_silane_gag()     : R<-SiH3                [slag]  (silane is kept)
//  make_silyl_gag()      : R<-Si(CH3)3            [slyg]  (silyls are generic silicone gag)
//  make_hphosphine_gag() : R<-PH2                 [hphg]  (3-coordinated phosphor bound to hydrogen)
//  make_mphosphine_gag() : R<-P-(CH3)2            [mphg]  (3-coordinated phosphor bound to methyls)
//  make_phosphane_gag()  : R<-P=CH2               [ppng]  (2-coordinated P near double bond)
//  make_phosphate2_gag() : R<-(P=O)~[O(-)]2       [pp2g]  (5-coordinated phosphor ion -2) 
//  make_phosphate1_gag() : R<-(H3C-P=O)~O(-)      [pp1g]  (5-coordinated phosphor ion -1) 
//  make_phosphate0_gag() : R<-(O=P)-(CH3)2        [pp0g]  (5-coordinated phosphor) 
//  make_sulfide_gag()    : R<-S-CH3               [sidg]  (2-coordinated sulfur)
//  make_thiol_gag()      : R<-SH                  [thig]  (2-coordinated sulfur bound to hydrogen)
//  make_isulfite_gag()   : R<-(S=O)~O(-)          [isig]  (4-coordinated sulfur ion)
//  make_msulfite_gag()   : R<-(S=O)-CH3           [msig]  (4-coordinated sulfur)
//  make_isulfate_gag()   : R<-(O=S=O)~O(-)        [isag]  (6-coordinated sulfur ion)
//  make_msulfate_gag()   : R<-(O=S=O)-CH3         [msag]  (6-coordinated sulfur)
//The gags bits modes (4+2+1 is a valid mode) 
// 1 - info
// 2 - topology
// 4 - geometry
//
//NB!! Only single (nonresonance, nonaromatic, nonamide etc) bond is allowed to separate anchors
//This function calls specific gag depending on it's name and pass the parameters through
char use_gags(unsigned int gag_name,char mode,unsigned int *res_name,unsigned int *atom_num,unsigned int *edge_num,t_vec *r,unsigned int acceptor_id,t_str *str)
{
     if (gag_name==*(unsigned int*)"cxyg") make_carboxyl_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);   // R<-[O=C~O(-)]        (carboxylic group is kept)
else if (gag_name==*(unsigned int*)"cgug") make_cguanidine_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str); // R<-[NH2~C(+)=NH2]    (guanidine group is kept)
else if (gag_name==*(unsigned int*)"cyng") make_cyane_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);      // R<-C_=_N             (cyane group is kept)
else if (gag_name==*(unsigned int*)"cmdg") make_camide_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);     // R<-(C=O)~NH2         (C near/in resonance)
else if (gag_name==*(unsigned int*)"metg") make_methyl_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);     // R<-CH3               (generic aliphatic C)
else if (gag_name==*(unsigned int*)"ethg") make_ethene_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);     // R<-CH=CH2            (alkene carbon)
else if (gag_name==*(unsigned int*)"camg") make_camine_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);     // R<-NH3+              (aliphatic charged N)
else if (gag_name==*(unsigned int*)"ngug") make_nguanidine_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str); // R<-XN~[NH2~C(+)=NH2] (guanidine group is kept)
else if (gag_name==*(unsigned int*)"pamg") make_pamine_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);     // R<-NH2               (N near aromatic if X is sp2 atom)
else if (gag_name==*(unsigned int*)"pimg") make_nimine_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);     // R<-N=CH2             (generic sp2 N)
else if (gag_name==*(unsigned int*)"nadg") make_namide_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);     // R<-NH~(C=O)H         (N near r-nce if X is sp3 atom)
else if (gag_name==*(unsigned int*)"no2g") make_nitro_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);      // R<-NO2               (NO2 group is kept)
else if (gag_name==*(unsigned int*)"hxyg") make_hydroxyl_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);   // R<-OH                (hydroxy gag) 
else if (gag_name==*(unsigned int*)"omeg") make_methoxy_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);    // R<-O-CH3             (generic oxygen gag)
else if (gag_name==*(unsigned int*)"slag") make_silane_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);     // R<-SiH3              (silane is kept)
else if (gag_name==*(unsigned int*)"slyg") make_silyl_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);      // R<-Si(CH3)3          (silyls are generic silicone gag)
else if (gag_name==*(unsigned int*)"hphg") make_hphosphine_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str); // R<-PH2               (3-coordinated P bound to H)
else if (gag_name==*(unsigned int*)"mphg") make_mphosphine_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str); // R<-P-(CH3)2          (3-coordinated P bound to Me)
else if (gag_name==*(unsigned int*)"ppng") make_phosphane_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);  // R<-P=CH2             (3-coordinated P near double bond)
else if (gag_name==*(unsigned int*)"pp2g") make_phosphate2_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str); // R<-(P=O)~[O(-)]2     (5-coordinated -2 phosphor) 
else if (gag_name==*(unsigned int*)"pp1g") make_phosphate1_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str); // R<-(H3C-P=O)~O(-)    (5-coordinated -1 phosphor) 
else if (gag_name==*(unsigned int*)"pp0g") make_phosphate0_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str); // R<-(O=P)-(CH3)2      (5-coordinated  0 phosphor) 
else if (gag_name==*(unsigned int*)"sidg") make_sulfide_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);    // R<-S-CH3             (2-coordinated sulfur)
else if (gag_name==*(unsigned int*)"thig") make_thiol_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);      // R<-SH                (2-coordinated sulfur bound to H)
else if (gag_name==*(unsigned int*)"isig") make_isulfite_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);   // R<-(S=O)~O(-)        (4-coordinated sulfur ion)
else if (gag_name==*(unsigned int*)"msig") make_msulfite_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);   // R<-(S=O)-CH3         (4-coordinated sulfur)
else if (gag_name==*(unsigned int*)"isag") make_isulfate_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);   // R<-(O=S=O)~O(-)      (6-coordinated sulfur ion)
else if (gag_name==*(unsigned int*)"msag") make_msulfate_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);   // R<-(O=S=O)-CH3       (6-coordinated sulfur)
else if (gag_name==*(unsigned int*)"pheg") make_phenyl_gag(mode,res_name,atom_num,edge_num,r,acceptor_id,str);     // R<-C6H5              (phenyl gag)
else { ylib_errno=YERROR_IMPOSSIBLE; return FALSE; } //No such gag exists - we shouldn't ever get here
return TRUE;
}

//This function identifies the most appropriate gag for _i -> _j bond
//Note. It returns TRUE on success and FALSE it he generic gag was used
char _identify_gag(unsigned int *gag_name,unsigned int *atom_num,unsigned int *edge_num, register unsigned int _id, register unsigned int _i, register unsigned int _j,int *charged,char (*order)[4],t_clist *neighbors,t_mol *mol)
{

//Stage I.1.2. Update counts
switch (mol->a[_j])
  {
  case CHEM_ATOM_TYPE_CARBON   : { //Carbon types
         if (order[_j][3]==1)
           {
           if (neighbors->list[_j].size==2) make_cyane_gag(0x1,gag_name,atom_num,edge_num,0x0,_i,0x0);
           else return FALSE;
           }
    else if ( (order[_j][2])) 
           {
           if (neighbors->list[_j].size==3)
             {
//R<-CO2(-): neighbors.size==3, order[2]==1, order[1]==1, neighbors[0&1].size==1, neighbors[0&1].chem_type=OXYGEN
                  if ( (charged[_j]==+1)&&(order[_j][2]==1)&&(order[_j][1]==1)&&
                     ( (neighbors->list[neighbors->list[_j].list[0]].size==1)&&( (mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_SULFUR) ) )&& 
                     ( (neighbors->list[neighbors->list[_j].list[1]].size==1)&&( (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_SULFUR) ) ) ) 
                    make_carboxyl_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
//R<-C(NH2)2(+): neighbors.size==3, order[2]==1, order[1]==1, neighbors[0&1].size==3, neighbors[0&1].chem_type==NITROGEN, each neighbors[2&3]_of_neighbors[0&1]==HYDROGEN
             else if ( (charged[_j]==-1)&&(order[_j][2]==1)&&(order[_j][1]==1)&&
                     ( (neighbors->list[neighbors->list[_j].list[0]].size==3)&&(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_NITROGEN)&&(mol->a[neighbors->list[neighbors->list[_j].list[0]].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mol->a[neighbors->list[neighbors->list[_j].list[0]].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) )&& 
                     ( (neighbors->list[neighbors->list[_j].list[1]].size==3)&&(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_NITROGEN)&&(mol->a[neighbors->list[neighbors->list[_j].list[1]].list[1]]==CHEM_ATOM_TYPE_HYDROGEN)&&(mol->a[neighbors->list[neighbors->list[_j].list[1]].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) ) ) 
                    make_cguanidine_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             else if ( (mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_CARBON)&&
                     ( (neighbors->list[_j].list[1]==_i)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_CARBON)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) )&&
                     ( (neighbors->list[_j].list[2]==_i)||(mol->a[neighbors->list[_j].list[2]]==CHEM_ATOM_TYPE_CARBON)||(mol->a[neighbors->list[_j].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
                    {
                    if ( (neighbors->list[_i].size==2)&&( (mol->a[_i]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[_i]==CHEM_ATOM_TYPE_SULFUR) )&&
                       ( (mol->a[neighbors->list[_i].list[0]]==CHEM_ATOM_TYPE_HYDROGEN)||(mol->a[neighbors->list[_i].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
                      make_phenyl_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                    else make_ethene_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                    }
             else   {
                    if ( (neighbors->list[_i].size==2)&&( (mol->a[_i]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[_i]==CHEM_ATOM_TYPE_SULFUR) )&&
                       ( (mol->a[neighbors->list[_i].list[0]]==CHEM_ATOM_TYPE_HYDROGEN)||(mol->a[neighbors->list[_i].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
                      make_phenyl_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                    else make_camide_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                    }
             }
           else return FALSE;
           }
    else if ( (order[_j][1])) return FALSE;
    else   {
           if (neighbors->list[_j].size==4) make_methyl_gag(0x1,gag_name,atom_num,edge_num,0x0,_i,0x0);
           else return FALSE;
           }
    return TRUE; }
  case CHEM_ATOM_TYPE_NITROGEN : { //Nitrogen types
         if ( (order[_j][3])) return FALSE; //There is no R<-N_=_X (yet)
    else if ( (order[_j][2])) //R<-N=X
           {
           if ( (charged[_j]))
             { //R<-NO2: neighbors.size==3, order[2]==1, order[1]==1, neighbors[0&1].size==1, neighbors[0&1].chem_type=OXYGEN
                  if ( (neighbors->list[_j].size==3)&&(order[_j][2]==1)&&(order[_j][1]==1)&&
                     ( (neighbors->list[neighbors->list[_j].list[0]].size==1)&&(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_OXYGEN) )&& 
                     ( (neighbors->list[neighbors->list[_j].list[1]].size==1)&&(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_OXYGEN) ) ) 
                    make_nitro_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             else if (neighbors->list[_j].size==2) make_nimine_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0); 
             else return FALSE;
             }
           else return FALSE;
           }
    else if ( (order[_j][1])) //R<-N...X
           {//R<-XN~[NH2~C(+)=NH2]
           if (neighbors->list[_j].size==3)
             {
             if ( (order[_j][1]==1)&&(order[_j][0]==2)&&
                (charged[*neighbors->list[_j].list]==+1)&&(mol->a[*neighbors->list[_j].list]==CHEM_ATOM_TYPE_CARBON)&&(neighbors->list[*neighbors->list[_j].list].size==3)&&(order[*neighbors->list[_j].list][2]==1)&&(order[*neighbors->list[_j].list][1]==2)&&
                (mol->a[neighbors->list[*neighbors->list[_j].list].list[0]]==CHEM_ATOM_TYPE_NITROGEN)&&(mol->a[neighbors->list[*neighbors->list[_j].list].list[1]]==CHEM_ATOM_TYPE_NITROGEN)&&(mol->a[neighbors->list[*neighbors->list[_j].list].list[2]]==CHEM_ATOM_TYPE_NITROGEN) )
               make_nguanidine_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             else make_namide_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             }
           else return FALSE;
           }
    else   {// R<-N(-X)2
           if (neighbors->list[_j].size==4)
             {
             if ( (charged[_j]==+1)) make_camine_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             else return FALSE;
             }
           else if (neighbors->list[_j].size==3) 
             {
             if ( (order[_i][0]!=neighbors->list[_i].size)) make_pamine_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             else make_namide_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             }
           else return FALSE;
           }
    return TRUE; }
  case CHEM_ATOM_TYPE_OXYGEN   : { //Oxygen types
    if ( (order[_j][0]==neighbors->list[_j].size))
      {
      if (neighbors->list[_j].size==2)
        {
        if ( ( (neighbors->list[_j].list[0]==_i)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_HYDROGEN) )&&
           ( (neighbors->list[_j].list[1]==_i)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
          make_hydroxyl_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
        else make_methoxy_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
        }
      else return FALSE;
      }
    else return FALSE;
    return TRUE; }
  case CHEM_ATOM_TYPE_SILICON  : { //Silicon types
    if (neighbors->list[_j].size==4)
      {
      if ( ( (neighbors->list[_j].list[0]==_i)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_HYDROGEN) )&&
           ( (neighbors->list[_j].list[1]==_i)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) )&&
           ( (neighbors->list[_j].list[2]==_i)||(mol->a[neighbors->list[_j].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) )&&
           ( (neighbors->list[_j].list[3]==_i)||(mol->a[neighbors->list[_j].list[3]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
        make_silane_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
      else make_silyl_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
      }
    else return FALSE;
    return TRUE; }
  case CHEM_ATOM_TYPE_PHOSPHOR : { //Phosphor types
         if ( (order[_j][3])) return FALSE;
    else if ( (order[_j][2]))
           {
                if ( (neighbors->list[_j].size==4)&&(order[_j][2]==1)&&
                   ( (mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_OXYGEN)  ||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_SULFUR)  ||
                     (mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_NITROGEN)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_PHOSPHOR) ) )
                  {
                       if ( (charged[_j]==-2)&&(order[_j][1]==2) ) make_phosphate2_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                  else if ( (charged[_j]==-1)&&( (order[_j][1])) ) make_phosphate1_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                  else if (!(charged[_j]))                         make_phosphate0_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                  else return FALSE;
                  }
           else if (neighbors->list[_j].size==2) make_phosphane_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
           else return FALSE; 
           }
    else if ( (order[_j][1])) return FALSE;
    else   {
           if (neighbors->list[_j].size==3)
             {
             if ( ( (neighbors->list[_j].list[0]==_i)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_HYDROGEN) )&&
                  ( (neighbors->list[_j].list[1]==_i)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) )&&
                  ( (neighbors->list[_j].list[2]==_i)||(mol->a[neighbors->list[_j].list[2]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
               make_hphosphine_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             else make_mphosphine_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             }
           else return FALSE;
           }
    return TRUE; }
  case CHEM_ATOM_TYPE_SULFUR   : { //Sulfur types
         if ( (order[_j][3])) return FALSE;
    else if ( (order[_j][2]))
           {
           if (neighbors->list[_j].size==3)
                  {
                       if ( (charged[_j]==-1)&(order[_j][2]==1)&&(order[_j][1]==1)&&
                          ( (mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_SULFUR) )&&
                          ( (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_SULFUR) ) )
                         make_isulfite_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                  else if ( (!(charged[_j]))&(order[_j][2]==1)&&
                          ( (mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_SULFUR)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_NITROGEN) ) )
                         make_msulfite_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                  else return FALSE; 
                  } 
           else if (neighbors->list[_j].size==4)
                  {
                       if ( (charged[_j]==-1)&(order[_j][2]==2)&&(order[_j][1]==1)&&
                          ( (mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_SULFUR) )&&
                          ( (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_SULFUR) )&&
                          ( (mol->a[neighbors->list[_j].list[2]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[2]]==CHEM_ATOM_TYPE_SULFUR) ) )
                         make_isulfate_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                  else if ( (!(charged[_j]))&(order[_j][2]==2)&&
                          ( (mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_SULFUR)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_NITROGEN) )&&
                          ( (mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_OXYGEN)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_SULFUR)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_NITROGEN) ) )
                         make_msulfate_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
                  else return FALSE; 
                  }
           else return FALSE;
           }
    else if ( (order[_j][1])) return FALSE;
    else   {
           if (neighbors->list[_j].size==2)
             {
             if ( ( (neighbors->list[_j].list[0]==_i)||(mol->a[neighbors->list[_j].list[0]]==CHEM_ATOM_TYPE_HYDROGEN) )&&
                ( (neighbors->list[_j].list[1]==_i)||(mol->a[neighbors->list[_j].list[1]]==CHEM_ATOM_TYPE_HYDROGEN) ) )
               make_thiol_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             else make_sulfide_gag(0x1,gag_name,atom_num,edge_num,0x0,0x0,0x0);
             }
           else return FALSE;
           } 
    return TRUE; }
  default : return FALSE;
  }
}

//This function generates set of mrr_strs: one per anchor, each consists of the anchor itself and all functional gags added into the end of the list.
//Return amount of manchors of -1 on error
//NB! it uses sstr->ress->list[0] to store id of the anchor in molecule and start_rid to store amount of in-anchor bonds
//NB! it uses anames field of each gag atom to store actual gagged neighboring atom from original mol
//The errors are YERROR_NIMPLEMENTED, YERROR_IMPOSSIBLE, YERROR_LEGAL, YERROR_MEMORY
unsigned int _splitt_anchors_mrr_str(t_str ***mrr_str,unsigned int *err_count,char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_adjacency *adjacency,t_mol *mol,unsigned int default_gag_name)
{
register unsigned int _i, _j, _k, _l, _t;
unsigned int res_name, atom_num, edge_num;
unsigned int *vertice, *neighbor_id;
int *charged;
char *brt_f;
t_vec *v;
t_str sstr;

//Stage 0. Prepare memory
if ( (!(*mrr_str))&&(!((*mrr_str)=(t_str**)calloc(adjacency->size_v,sizeof(t_str*)))) ) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return (unsigned int)-1; }
_i=0, _j=adjacency->size_v; while (--_j) if (adjacency->nn[_j]>adjacency->nn[_i]) _i=_j;
if (!(charged=(int*)malloc(mol->natoms*sizeof(int)+adjacency->nn[_i]*(sizeof(t_vec)+sizeof(unsigned int)+sizeof(unsigned int))))) 
  { LABEL_MEMORY_ERROR_1: free(*mrr_str); *mrr_str=0x0; goto LABEL_MEMORY_ERROR_0; }
else
  {
  vertice=(unsigned int*)((void*)charged+sizeof(int)*mol->natoms);
  neighbor_id=vertice+adjacency->nn[_i];
  v=(t_vec*)((void*)neighbor_id+sizeof(unsigned int)*adjacency->nn[_i]);
  //Mark charged atoms
  memset(charged,0,mol->natoms), _l=mol->nvedges; while (_l--) charged[mol->vedges[_l][0]]=mol->vatoms[mol->vedges[_l][1]]; 
  }
//Stage I. Create mono-strs
for (_i=0; _i<adjacency->size_v; _i++)
  {
  //Stage I.1. Define gags
  //Stage I.1.1. Init sstr
  memset(&sstr,0,sizeof(t_str));
  if (!(sstr.ress=(t_list*)alloc_list(adjacency->nn[_i]+1)))
    { 
    LABEL_MEMORY_ERROR_2: free(charged);
    while (_i--) free((*mrr_str)[_i]); goto LABEL_MEMORY_ERROR_1; 
    }
  sstr.ress->list[0]=_i, sstr.ress->size=1;
  if (!(sstr.rsize=(unsigned int*)malloc(sizeof(unsigned int)*(adjacency->nn[_i]+1+1))))
    { if (sstr.ress) free(sstr.ress), sstr.ress=0x0; goto LABEL_MEMORY_ERROR_2; }
  else { sstr.rsize[0]=0, sstr.rsize[1]=sstr.natoms=mol->anchors->list[_i].size; }
  sstr.nedges=sstr.start_rid=0; 
  //Stage I.1.2. Process anchor neighbors
  _j=mol->anchors->list[_i].size;
  while (_j--)
    {
    _k=mol->anchors->list[_i].list[_j], _l=neighbors->list[_k].size;
    while (_l--)
      {
      if (anchor_id[_t=neighbors->list[_k].list[_l]]==_i) { if (_k<_t) sstr.start_rid++; }
      else
        {//Ideentify the most appropriate gag
        vertice[sstr.ress->size-1]=_j;
        v[sstr.ress->size-1].i=mol->r[_t].i-mol->r[_k].i, v[sstr.ress->size-1].j=mol->r[_t].j-mol->r[_k].j, v[sstr.ress->size-1].k=mol->r[_t].k-mol->r[_k].k;
        neighbor_id[sstr.ress->size-1]=_t;
        if (calc_vec_norm(&v[sstr.ress->size-1])<SMALL2)
          {//Geometry ERROR
          ylib_errno=YERROR_LEGAL; free(sstr.ress);  sstr.ress=0x0; free(sstr.rsize); sstr.rsize=0x0;
          while (_i--) free((*mrr_str)[_i]); free(*mrr_str), *mrr_str=0x0; free(charged); 
          return (unsigned int)-1;
          }   
        if (!(_identify_gag(&res_name,&atom_num,&edge_num,_i,_k,_t,charged,order,neighbors,mol))) 
          {
          (*err_count)++; //Have no idea of what the proper gag is, using generic gag 
          if (!(use_gags(default_gag_name,0x1,&res_name,&atom_num,&edge_num,0x0,_k,&sstr))) 
            {
            LABEL_YERROR_IMPOSSIBLE: ylib_errno=YERROR_IMPOSSIBLE; free(sstr.ress);  sstr.ress=0x0; free(sstr.rsize); sstr.rsize=0x0;
            while (_i--) free((*mrr_str)[_i]); free(*mrr_str), *mrr_str=0x0; free(charged); 
            return (unsigned int)-1;
            }  
          }
        sstr.ress->list[sstr.ress->size]=res_name, sstr.rsize[++sstr.ress->size]=_j, sstr.natoms+=atom_num, sstr.nedges+=edge_num;
        }
      }
    } 
  //Stage I.2. Compose sparse mrr_str
  if (!(sstr.anames=(char (*)[sizeof(int)])malloc(sizeof(int)*sstr.natoms)))    
    { LABEL_MEMORY_ERROR_3:      if (sstr.ress) { free(sstr.ress), sstr.ress=0x0; free(sstr.rsize),   sstr.rsize=0x0; } goto LABEL_MEMORY_ERROR_2; }
  if (!(sstr.a=(char*)malloc(sizeof(char)*sstr.natoms)))   { LABEL_MEMORY_ERROR_4: free(sstr.anames), sstr.anames=0x0;  goto LABEL_MEMORY_ERROR_3; }
  if (!(sstr.r=(t_vec*)malloc(sizeof(t_vec)*sstr.natoms))) { LABEL_MEMORY_ERROR_5: free(sstr.a),      sstr.a=0x0;       goto LABEL_MEMORY_ERROR_4; }
  if (!(_j=sstr.nedges+sstr.start_rid)) sstr.edges=0x0;
  else if (!(sstr.edges=(t_edge*)malloc(sizeof(t_edge)*_j)))                    { LABEL_MEMORY_ERROR_6: free(sstr.r),      sstr.r=0x0;      goto LABEL_MEMORY_ERROR_5; }
  //Stage I.3. Fill str
  //Stage I.3.1. Copy core anchor
  sstr.natoms=mol->anchors->list[_i].size, sstr.nedges=0; 
  for (_j=0;_j<mol->anchors->list[_i].size;_j++)
    {
    _k=mol->anchors->list[_i].list[_j];
    sstr.a[_j]=mol->a[_k];
    *((int*)sstr.anames[_j])=*((int*)mol->anames[_k]);
    sstr.r[_j].i=mol->r[_k].i, sstr.r[_j].j=mol->r[_k].j, sstr.r[_j].k=mol->r[_k].k;
    for (_l=0;_l<neighbors->list[_k].size;_l++)
      {//This logic is applicable here as we are not interested in gags
      _t=neighbors->list[_k].list[_l];
      if ( (_t>_k)&&(anchor_id[_t]==_i)&&((atom_num=find_in_list(_t,&mol->anchors->list[_i]))!=(unsigned int)-1) )
        {
        sstr.edges[sstr.nedges].vertice[0]=_j, sstr.edges[sstr.nedges].vertice[1]=atom_num;
             if (_l<order[_k][3])              sstr.edges[sstr.nedges].type=(int)'3';
        else if (_l<order[_k][3]+order[_k][2]) sstr.edges[sstr.nedges].type=(int)'2';
        else                                   sstr.edges[sstr.nedges].type=(int)'1'; 
        sstr.nedges++;
        }
      } 
    }
  if (sstr.start_rid!=sstr.nedges) goto LABEL_YERROR_IMPOSSIBLE; //We shouldn't ever get here
  //Stage I.3.2. Fill gags
  for (_j=1;_j<sstr.ress->size;_j++)
    {
    res_name=_j, atom_num=sstr.natoms, edge_num=sstr.nedges;
    if (!(use_gags(sstr.ress->list[_j],0x7,&res_name,&atom_num,&edge_num,&v[_j-1],vertice[_j-1],&sstr))) goto LABEL_YERROR_IMPOSSIBLE;
    else 
      {//Update sstr
       //Save neighboring id into anames field
      *((unsigned int*)&sstr.anames[sstr.natoms])=neighbor_id[_j-1];
      sstr.natoms+=atom_num, sstr.nedges+=edge_num, sstr.rsize[_j+1]=sstr.natoms;
      } 
    }
  //Stage I.4. Copy solid memory mrr_str
  if (!(_j=calc_brutto_f(&brt_f,sstr.natoms,sstr.a))) { LABEL_MEMORY_ERROR_7: if ( (sstr.nedges)) free(sstr.edges), sstr.edges=0x0; goto LABEL_MEMORY_ERROR_6; }
  if (!((*mrr_str)[_i]=(t_str*)alloc_solid_str(sizeof(char)*_j,sstr.ress->size,sstr.natoms,sstr.nedges))) { free(brt_f); goto LABEL_MEMORY_ERROR_7; }
  (*mrr_str)[_i]->natoms=sstr.natoms, (*mrr_str)[_i]->nedges=sstr.nedges; 
  memcpy((*mrr_str)[_i]->name,brt_f,sizeof(char)*_j); free(brt_f);
  if (sstr.ress->size==1) 
    { (*mrr_str)[_i]->ress=0x0, *((unsigned int*)&(*mrr_str)[_i]->rsize)=*sstr.ress->list; }
  else
    {
    (*mrr_str)[_i]->ress->size=sstr.ress->size;
    memcpy((*mrr_str)[_i]->ress->list,sstr.ress->list,sizeof(unsigned int)*sstr.ress->size);
    memcpy((*mrr_str)[_i]->rsize,sstr.rsize,sizeof(unsigned int)*(sstr.ress->size+1));
    }
  memcpy((*mrr_str)[_i]->anames,sstr.anames,sizeof(char)*0x4*sstr.natoms); 
  memcpy((*mrr_str)[_i]->a,sstr.a,sizeof(char)*sstr.natoms); 
  memcpy((*mrr_str)[_i]->r,sstr.r,sizeof(t_vec)*sstr.natoms); 
  if ( (sstr.nedges)) memcpy((*mrr_str)[_i]->edges,sstr.edges,sizeof(t_edge)*sstr.nedges);
  else (*mrr_str)[_i]->edges=0x0;
  (*mrr_str)[_i]->start_rid=sstr.start_rid;
  free_str(&sstr);
  }
//Copy memory and exit
free(charged); charged=0x0;
return _i;
}
//This function generates set of brr_strs: one per anchor-anchor bond, each consists of two anchors itself and all functional gags added into the end of the list.
//Return amount of brr_strs of -1 on error
//NB! it uses brr_str->ress->list[0] and brr_str->ress->list[1] to store id of anchors in the molecule, and start_rid to store amount of in-anchor bonds. The brr_str->edges[0] is the edge connecting two anchors.
//The errors are YERROR_NIMPLEMENTED, YERROR_IMPOSSIBLE, YERROR_MEMORY
unsigned int _splitt_anchors_brr_str(t_str ***brr_str,unsigned int *err_count,char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_adjacency *adjacency,t_mol *mol,unsigned int default_gag_name)
{
register unsigned int _i, _j, _k, _l, _p, _q, _t;
unsigned int res_name, atom_num, edge_num;
unsigned int *vertice, a1, a2, brr_size;
int *charged;
char  *brt_f;
t_vec *v;
t_str cstr;

//Stage 0. Prepare memory
if ( (!(*brr_str))&&(!((*brr_str)=(t_str**)calloc(adjacency->size_e,sizeof(t_str*)))) ) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return (unsigned int)-1; }
_i=0, _j=adjacency->size_v; while (--_j) if (adjacency->nn[_j]>adjacency->nn[_i]) _i=_j;
if (!(charged=(int*)malloc((mol->natoms*sizeof(int)+2*adjacency->nn[_i]*(sizeof(t_vec)+sizeof(unsigned int)))))) { LABEL_MEMORY_ERROR_1: free(*brr_str); *brr_str=0x0; goto LABEL_MEMORY_ERROR_0; }
else
  {
  vertice=(unsigned int*)((void*)charged+sizeof(int)*mol->natoms);
  v=(t_vec*)((void*)vertice+sizeof(unsigned int)*2*adjacency->nn[_i]);
  //Mark charged atoms
  memset(charged,0,mol->natoms), _l=mol->nvedges; while (_l--) charged[mol->vedges[_l][0]]=mol->vatoms[mol->vedges[_l][1]]; 
  }
//Stage I. Generate couples
for ( brr_size=0, _i=0; _i<adjacency->size_v; _i++)
  for ( _k=0; _k<adjacency->nn[_i]; _k++)
    if ((_j=adjacency->nv[_i][_k])>_i)
      {
      //Stage I.1. Default edge
      //Stage I.1.1. Init cstr
      memset(&cstr,0,sizeof(t_str));
      cstr.natoms=mol->anchors->list[_i].size+mol->anchors->list[_j].size;
      //Stage I.1.1.a.Compose ress
      if (!(cstr.ress=(t_list*)alloc_list(adjacency->nn[_i]+adjacency->nn[_j]+1)))
        { LABEL_MEMORY_ERROR_2: free(charged); while (brr_size--) free((*brr_str)[brr_size]); goto LABEL_MEMORY_ERROR_1; }
      else { cstr.ress->list[0]=_i, cstr.ress->list[1]=_j, cstr.ress->size=2; }
      //Stage I.1.1.a.Compose rsizes
      if (!(cstr.rsize=(unsigned int*)malloc(sizeof(unsigned int)*(adjacency->nn[_i]+adjacency->nn[_j]+2))))
        { LABEL_MEMORY_ERROR_3: free(cstr.ress), cstr.ress=0x0; goto LABEL_MEMORY_ERROR_2; } 
      else { cstr.rsize[0]=0, cstr.rsize[1]=mol->anchors->list[_i].size, cstr.rsize[2]=mol->anchors->list[_i].size+mol->anchors->list[_j].size; }
      //Stage I.1.2. Process anchor#1 neighbors
      _l=mol->anchors->list[_i].size;
      while (_l--)
        {
        _p=mol->anchors->list[_i].list[_l], _q=neighbors->list[_p].size;
        while (_q--)
          {
          _t=neighbors->list[_p].list[_q];
               if (anchor_id[_t]==_i) { if (_p<_t) cstr.start_rid++; }
          else if (anchor_id[_t]==_j) { a1=_p, a2=_t; }
          else
            {//Ideentify the most appropriate gag
            vertice[cstr.ress->size-2]=_l;
            v[cstr.ress->size-2].i=mol->r[_t].i-mol->r[_p].i, v[cstr.ress->size-2].j=mol->r[_t].j-mol->r[_p].j, v[cstr.ress->size-2].k=mol->r[_t].k-mol->r[_p].k;
            if (calc_vec_norm(&v[cstr.ress->size-2])<SMALL2)
              {//Geometry ERROR
              ylib_errno=YERROR_LEGAL; free(cstr.ress);  cstr.ress=0x0; free(cstr.rsize); cstr.rsize=0x0;
              while (brr_size--) free((*brr_str)[brr_size]); free(*brr_str), *brr_str=0x0; free(charged); 
              return (unsigned int)-1;
              }   
            if (!(_identify_gag(&res_name,&atom_num,&edge_num,_i,_p,_t,charged,order,neighbors,mol))) 
              {
              (*err_count)++; //Have no idea of what the proper gag is, using generic gag 
              if (!(use_gags(default_gag_name,0x1,&res_name,&atom_num,&edge_num,0x0,_p,&cstr))) 
                {
                LABEL_YERROR_IMPOSSIBLE: ylib_errno=YERROR_IMPOSSIBLE; free(cstr.ress);  cstr.ress=0x0; free(cstr.rsize); cstr.rsize=0x0;
                while (brr_size--) free((*brr_str)[brr_size]); free(*brr_str), *brr_str=0x0; free(charged); 
                return (unsigned int)-1;
                }  
              }
            cstr.ress->list[cstr.ress->size]=res_name, cstr.rsize[++cstr.ress->size]=_j, cstr.natoms+=atom_num, cstr.nedges+=edge_num;
            }
          }
        }
      //Stage I.1.3. Process anchor#2 neighbors
      _l=mol->anchors->list[_j].size;
      while (_l--)
        {
        _p=mol->anchors->list[_j].list[_l], _q=neighbors->list[_p].size;
        while (_q--)
          {
          _t=neighbors->list[_p].list[_q];
               if (anchor_id[_t]==_j) { if (_p<_t) cstr.start_rid++; }
          else if (anchor_id[_t]==_i) { if ( (a2!=_p)||(a1!=_t) ) goto LABEL_YERROR_IMPOSSIBLE; }
          else
            {//Ideentify the most appropriate gag
            vertice[cstr.ress->size-2]=mol->anchors->list[_i].size+_l;
            v[cstr.ress->size-2].i=mol->r[_t].i-mol->r[_p].i, v[cstr.ress->size-2].j=mol->r[_t].j-mol->r[_p].j, v[cstr.ress->size-2].k=mol->r[_t].k-mol->r[_p].k;
            if (calc_vec_norm(&v[cstr.ress->size-2])<SMALL2)
              {//Geometry ERROR
              ylib_errno=YERROR_LEGAL; free(cstr.ress);  cstr.ress=0x0; 
              while (brr_size--) free((*brr_str)[brr_size]); free(*brr_str), *brr_str=0x0; free(charged); 
              return (unsigned int)-1;
              }   
            if (!(_identify_gag(&res_name,&atom_num,&edge_num,_j,_p,_t,charged,order,neighbors,mol))) 
              {
              (*err_count)++; //Have no idea of what the proper gag is, using generic gag 
              if (!(use_gags(default_gag_name,0x1,&res_name,&atom_num,&edge_num,0x0,_p,&cstr))) goto LABEL_YERROR_IMPOSSIBLE;
              }
            cstr.ress->list[cstr.ress->size]=res_name, cstr.rsize[++cstr.ress->size]=_j, cstr.natoms+=atom_num, cstr.nedges+=edge_num;
            }
          }
        }
      //Stage II. Compose sparse brr_str | anchor_i | anchor_j | gags_i | gags_j | 
      if (!(cstr.anames=(char (*)[sizeof(int)])malloc(sizeof(int)*cstr.natoms))) { LABEL_MEMORY_ERROR_4: free(cstr.rsize),  cstr.rsize=0x0;  goto LABEL_MEMORY_ERROR_3; }
      if (!(cstr.a=(char*)malloc(sizeof(char)*cstr.natoms)))                     { LABEL_MEMORY_ERROR_5: free(cstr.anames), cstr.anames=0x0; goto LABEL_MEMORY_ERROR_4; }
      if (!(cstr.r=(t_vec*)malloc(sizeof(t_vec)*cstr.natoms)))                   { LABEL_MEMORY_ERROR_6: free(cstr.a),      cstr.a=0x0;      goto LABEL_MEMORY_ERROR_5; }
      if (!(cstr.edges=(t_edge*)malloc(sizeof(t_edge)*(cstr.nedges+cstr.start_rid+1)))) { LABEL_MEMORY_ERROR_7: free(cstr.r), cstr.r=0x0;    goto LABEL_MEMORY_ERROR_6; }
      //Stage II.1. Fill str
      cstr.natoms=mol->anchors->list[_i].size+mol->anchors->list[_j].size, cstr.nedges=1; //Reserve 0-bond for the inter-anchor bond 
      if ( (((*cstr.edges).vertice[0]=find_in_list(a1,&mol->anchors->list[_i]))==(unsigned int)-1)||(((*cstr.edges).vertice[1]=find_in_list(a2,&mol->anchors->list[_j]))==(unsigned int)-1) ) goto LABEL_YERROR_IMPOSSIBLE;
      else { (*cstr.edges).vertice[1]+=mol->anchors->list[_i].size, (*cstr.edges).type=(int)'1'; }
      //Stage II.1.1. Copy core anchor#1
      for (_l=0; _l<mol->anchors->list[_i].size; _l++)
        {
        _p=mol->anchors->list[_i].list[_l];
        cstr.a[_l]=mol->a[_p];
        *((int*)cstr.anames[_l])=*((int*)mol->anames[_p]);
        cstr.r[_l].i=mol->r[_p].i, cstr.r[_l].j=mol->r[_p].j, cstr.r[_l].k=mol->r[_p].k;
        for (_q=0;_q<neighbors->list[_p].size;_q++)
          {//This logic is applicable here as we are not interested in gags
          _t=neighbors->list[_p].list[_q];
          if ( (_t>_p)&&(anchor_id[_t]==_i)&&((atom_num=find_in_list(_t,&mol->anchors->list[_i]))!=(unsigned int)-1) )
            {
            cstr.edges[cstr.nedges].vertice[0]=_l, cstr.edges[cstr.nedges].vertice[1]=atom_num;
                 if (_q<order[_p][3])              cstr.edges[cstr.nedges].type=(int)'3';
            else if (_q<order[_p][3]+order[_p][2]) cstr.edges[cstr.nedges].type=(int)'2';
            else                                   cstr.edges[cstr.nedges].type=(int)'1'; 
            cstr.nedges++;
            }
          } 
        }
      //Stage II.1.2. Copy core anchor#2
      for (_l=0; _l<mol->anchors->list[_j].size; _l++)
        {
        _p=mol->anchors->list[_j].list[_l];
        cstr.a[mol->anchors->list[_i].size+_l]=mol->a[_p];
        *((int*)cstr.anames[mol->anchors->list[_i].size+_l])=*((int*)mol->anames[_p]);
        cstr.r[mol->anchors->list[_i].size+_l].i=mol->r[_p].i, cstr.r[mol->anchors->list[_i].size+_l].j=mol->r[_p].j, cstr.r[mol->anchors->list[_i].size+_l].k=mol->r[_p].k;
        for (_q=0;_q<neighbors->list[_p].size;_q++)
          {//This logic is applicable here as we are not interested in gags
          _t=neighbors->list[_p].list[_q];
          if ( (_t>_p)&&(anchor_id[_t]==_j)&&((atom_num=find_in_list(_t,&mol->anchors->list[_j]))!=(unsigned int)-1) )
            {
            cstr.edges[cstr.nedges].vertice[0]=mol->anchors->list[_i].size+_l, cstr.edges[cstr.nedges].vertice[1]=mol->anchors->list[_i].size+atom_num;
                 if (_q<order[_p][3])                                          cstr.edges[cstr.nedges].type=(int)'3';
            else if (_q<order[_p][3]+order[_p][2])                             cstr.edges[cstr.nedges].type=(int)'2';
            else                                                               cstr.edges[cstr.nedges].type=(int)'1'; 
            cstr.nedges++;
            }
          } 
        }
      if (1+cstr.start_rid!=cstr.nedges) goto LABEL_YERROR_IMPOSSIBLE; //We shouldn't ever get here
      //Stage I.3.2. Fill gags
      for (_l=2;_l<cstr.ress->size;_l++)
        {
        res_name=_l, atom_num=cstr.natoms, edge_num=cstr.nedges;
        if (!(use_gags(cstr.ress->list[_l],0x7,&res_name,&atom_num,&edge_num,&v[_l-2],vertice[_l-2],&cstr))) goto LABEL_YERROR_IMPOSSIBLE;
        else { cstr.natoms+=atom_num, cstr.nedges+=edge_num, cstr.rsize[_l+1]=cstr.natoms; } //Update cstr 
        }
      //Stage I.3. Copy solid memory brr_str
      if (!(_j=calc_brutto_f(&brt_f,cstr.natoms,cstr.a)))
        { LABEL_MEMORY_ERROR_8: free(cstr.edges), cstr.edges=0x0; goto LABEL_MEMORY_ERROR_7; }
      if (!((*brr_str)[brr_size]=(t_str*)alloc_solid_str(sizeof(char)*_j,cstr.ress->size,cstr.natoms,cstr.nedges)))
        {                       free(brt_f),      brt_f=0x0;      goto LABEL_MEMORY_ERROR_8; }
      (*brr_str)[brr_size]->natoms=cstr.natoms, (*brr_str)[brr_size]->nedges= cstr.nedges; 
      memcpy((*brr_str)[brr_size]->name,brt_f,sizeof(char)*_j); free(brt_f), brt_f=0x0;
      (*brr_str)[brr_size]->ress->size=cstr.ress->size, memcpy((*brr_str)[brr_size]->ress->list,cstr.ress->list,sizeof(unsigned int)*cstr.ress->size);
      memcpy((*brr_str)[brr_size]->rsize,cstr.rsize,sizeof(unsigned int)*(cstr.ress->size+1));
      memcpy((*brr_str)[brr_size]->anames,cstr.anames,sizeof(int)*cstr.natoms); 
      memcpy((*brr_str)[brr_size]->a,cstr.a,sizeof(char)*cstr.natoms); 
      memcpy((*brr_str)[brr_size]->r,cstr.r,sizeof(t_vec)*cstr.natoms); 
      memcpy((*brr_str)[brr_size]->edges,cstr.edges,sizeof(t_edge)*cstr.nedges); 
      free_str(&cstr);
      brr_size++;
      }
//All structures are done
free(charged); charged=0x0;
return brr_size;
}
//This function split molecule into a collection of mono- and bi-root strs
//It is a wrapper to call the previous two
char splitt_anchors(unsigned int *mrr_size,t_str ***mrr_str,unsigned int *brr_size,t_str ***brr_str,unsigned int maxerrcount,char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_adjacency *adjacency,t_mol *mol,unsigned int default_gag_name)
{
unsigned int err_count=0;
//Do mono-root split first (its results are required for the second stage)
if (((*mrr_size)=_splitt_anchors_mrr_str(mrr_str,&err_count,order,neighbors,anchor_id,adjacency,mol,default_gag_name))==(unsigned int)-1) return FALSE;
if (err_count>=maxerrcount) { ylib_errno=YERROR_LEGAL; goto LABEL_ERROR_EXIT; }
if (((*brr_size)=_splitt_anchors_brr_str(brr_str,&err_count,order,neighbors,anchor_id,adjacency,mol,default_gag_name))!=(unsigned int)-1) return TRUE;
else { LABEL_ERROR_EXIT: while ((*mrr_size)--) free((*mrr_str)[(*mrr_size)]); free(*mrr_str); *mrr_str=0x0; return FALSE; }
}




/*************************************   R T R E E     P A R T    **********************************/





//This function alloc rotation tree data structure
t_rtree *alloc_rtree(unsigned int nrbranch,unsigned int nnrbranch)
{
t_rtree *rtree;
if (!(rtree=(t_rtree*)malloc(sizeof(t_rtree)+sizeof(unsigned int)*nnrbranch*2+(nrbranch+1)*sizeof(t_rbranch)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
rtree->rbranch=(void*)rtree+sizeof(t_rtree);
rtree->_rbranch=(void*)rtree+sizeof(t_rtree)+sizeof(t_rbranch)*(nrbranch+1);
rtree->_nrbranch=nnrbranch*2;
rtree->nrbranch=0;
rtree->rbranch[nrbranch].nrbranch=0; //init superroot
return rtree;
}

//This function writes rtree to hdd 
char write_rtree(FILE *out,t_rtree *rtree)
{
unsigned int i;
i=Y_MAGIC;
if ( (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&rtree->nrbranch,sizeof(unsigned int),0x1,out)!=0x1)           ||   
     (fwrite(&rtree->_nrbranch,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&rtree->root,sizeof(unsigned int),0x1,out)!=0x1)||
     (fwrite(&rtree->nidofs,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(rtree->_rbranch,sizeof(unsigned int),rtree->_nrbranch,out)!=rtree->_nrbranch) )
  { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
for (i=0;i<rtree->nrbranch;i++)  
  if ( (fwrite(&rtree->rbranch[i].edge,sizeof(t_edge),0x1,out)!=0x1)||(fwrite(&rtree->rbranch[i].nrbranch,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&rtree->rbranch[i].cross,sizeof(unsigned int),0x1,out)!=0x1) )
    goto LABEL_IO_ERROR;  
return TRUE;
}

//This function reads rtree from hdd
t_rtree *read_rtree(FILE *in)
{
unsigned int _i, nrbranch, _nrbranch;
t_rtree *rtree;
if (fread(&_nrbranch,sizeof(unsigned int),0x1,in)!=0x1) { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if (_nrbranch!=Y_MAGIC) { ylib_errno=YERROR_USER; return FALSE; }
if ( (fread(&nrbranch,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&_nrbranch,sizeof(unsigned int),0x1,in)!=0x1) ) goto LABEL_IO_ERROR;
if (!(rtree=alloc_rtree(nrbranch,_nrbranch))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { rtree->nrbranch=nrbranch; rtree->_nrbranch=_nrbranch; }
if ( (fread(&rtree->root,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&rtree->nidofs,sizeof(unsigned int),0x1,in)!=0x1)||(fread(rtree->_rbranch,sizeof(unsigned int),rtree->_nrbranch,in)!=rtree->_nrbranch) )
  { free(rtree); goto LABEL_IO_ERROR; }
for (_i=nrbranch=0; _i<rtree->nrbranch; nrbranch+=rtree->rbranch[_i].nrbranch+rtree->rbranch[_i].cross, _i++)  
  if ( (fread(&rtree->rbranch[_i].edge,sizeof(t_edge),0x1,in)!=0x1)||(fread(&rtree->rbranch[_i].nrbranch,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&rtree->rbranch[_i].cross,sizeof(unsigned int),0x1,in)!=0x1) )
    { free(rtree); goto LABEL_IO_ERROR; }  
  else rtree->rbranch[_i].rbranch=&rtree->_rbranch[nrbranch]; 
rtree->rbranch[rtree->nrbranch].nrbranch=0, rtree->rbranch[rtree->nrbranch].rbranch=0x0;
return rtree;
}

//This function build rtree for given anchors complex list and massive of atomic edges
t_rtree *build_rtree(unsigned int root,unsigned int *anchor_id,t_clist *anchors,unsigned int nedges,t_edge *edges)
{
register unsigned int _i, _j, _k, _l;
unsigned int *data;
t_rtree *rtree=0x0;

//Stage I. Alloc memory
if (root>anchors->size) { ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
else { _k=0, _l=nedges; while (_l--) if (anchor_id[edges[_l].vertice[0]]!=anchor_id[edges[_l].vertice[1]]) _k++; } //Amount of inter-anchors edges
if (!(rtree=alloc_rtree(anchors->size,_k))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
else 
  {
  rtree->root=root, rtree->nrbranch=anchors->size; 
  _i=anchors->size+1; while (_i--) { rtree->rbranch[_i].edge.type=rtree->rbranch[_i].nrbranch=rtree->rbranch[_i].cross=0; }
  }
//Stage II. Build basic rtree
rtree->nidofs=0;
if (anchors->size==1) { rtree->_rbranch=0; return rtree; }
_l=nedges; while (_l--) if ((_i=anchor_id[edges[_l].vertice[0]])!=(_j=anchor_id[edges[_l].vertice[1]])) { rtree->rbranch[_i].nrbranch++, rtree->rbranch[_j].nrbranch++; }
for (_l=_k=0;_l<anchors->size;_k+=rtree->rbranch[_l++].nrbranch) { rtree->rbranch[_l].rbranch=&rtree->_rbranch[_k], rtree->rbranch[_l].edge.type=0; }
rtree->rbranch[rtree->nrbranch].rbranch=&rtree->_rbranch[_k], rtree->rbranch[rtree->nrbranch].edge.type=0;
_l=nedges; while (_l--) if ((_i=anchor_id[edges[_l].vertice[0]])!=(_j=anchor_id[edges[_l].vertice[1]])) { rtree->rbranch[_i].rbranch[rtree->rbranch[_i].edge.type++]=_j, rtree->rbranch[_j].rbranch[rtree->rbranch[_j].edge.type++]=_i; }
//Stage III. Order branches and map crosses
//Stage III.1. Order basic rtree 
if (!(data=(unsigned int*)calloc(anchors->size,sizeof(unsigned int)))) { free(rtree); goto LABEL_MEMORY_ERROR; } 
data[root]=1, rtree->rbranch[root].edge.type=rtree->rbranch[root].nrbranch; 
while (rtree->rbranch[root].edge.type--)
  {
  _i=rtree->rbranch[root].rbranch[rtree->rbranch[root].edge.type];
  data[_i]=data[root]+1;
  if (*rtree->rbranch[_i].rbranch!=root)  
    {//Find root and move it to the head
    if ((_l=find_in_row(root,rtree->rbranch[_i].nrbranch,rtree->rbranch[_i].rbranch))!=(unsigned int)-1) { rtree->rbranch[_i].rbranch[_l]=*rtree->rbranch[_i].rbranch, *rtree->rbranch[_i].rbranch=root; }
    else { LABEL_INTERNAL_CODE_ERROR: ylib_errno=YERROR_INTERNAL_CODE; free(data); free(rtree); return FALSE; }
    }
  rtree->rbranch[_i].edge.type=rtree->rbranch[_i].nrbranch;
  do{
    while(--rtree->rbranch[_i].edge.type)
      {
      _j=rtree->rbranch[_i].rbranch[rtree->rbranch[_i].edge.type];
      if ( (data[_j]))
        {//Remove j-th connect
        _k=--rtree->rbranch[_i].nrbranch, rtree->rbranch[_i].cross++, _l=rtree->rbranch[_i].edge.type;
        if (_k!=_l) { rtree->rbranch[_i].rbranch[_l]=rtree->rbranch[_i].rbranch[_k], rtree->rbranch[_i].rbranch[_k]=_j; }
        if ( (!(_l=xfind_in_row(_i,rtree->rbranch[_j].nrbranch,rtree->rbranch[_j].rbranch)))&&(_j!=root) ) goto LABEL_INTERNAL_CODE_ERROR;
        _k=--rtree->rbranch[_j].nrbranch, rtree->rbranch[_j].cross++; if (_l<=rtree->rbranch[_j].edge.type) rtree->rbranch[_j].edge.type--;
        while (_k!=_l++) rtree->rbranch[_j].rbranch[_l-1]=rtree->rbranch[_j].rbranch[_l]; rtree->rbranch[_j].rbranch[_k]=_i;        
        }
      else
        {
        if (*rtree->rbranch[_j].rbranch!=_i)  
          {//Find root and move it to the head
          if (!(_l=xfind_in_row(_i,rtree->rbranch[_j].nrbranch,rtree->rbranch[_j].rbranch))) goto LABEL_INTERNAL_CODE_ERROR;
          rtree->rbranch[_j].rbranch[_l]=*rtree->rbranch[_j].rbranch; *rtree->rbranch[_j].rbranch=_i;
          }
        data[_j]=data[_i]+1, _i=_j, rtree->rbranch[_i].edge.type=rtree->rbranch[_i].nrbranch;
        }
      }
    }while ((_i=*rtree->rbranch[_i].rbranch)!=root); 
  }
//Stage III.2. Check for islands
_l=anchors->size; while (_l--) if (!(data[_l])) { free(data); free(rtree); ylib_errno=YERROR_LEGAL; return FALSE; } //an island detected 
//Stage IV. Setup edges 
rtree->nidofs=0; 
_l=nedges;
while (_l--) //Run to order edges 
  if ((_i=anchor_id[edges[_l].vertice[0]])!=(_j=anchor_id[edges[_l].vertice[1]]))
    {
    if ( (data[_i]>data[_j])&&(*rtree->rbranch[_i].rbranch==_j) )
      {
      rtree->rbranch[_i].edge.vertice[0]=edges[_l].vertice[0], rtree->rbranch[_i].edge.vertice[1]=edges[_l].vertice[1];
      rtree->nidofs++;
      }
    if ( (data[_j]>data[_i])&&(*rtree->rbranch[_j].rbranch==_i) )
      {
      rtree->rbranch[_j].edge.vertice[0]=edges[_l].vertice[1], rtree->rbranch[_j].edge.vertice[1]=edges[_l].vertice[0];
      rtree->nidofs++;
      }
    }
//That is it, free some memory and exit
free(data);
return rtree;
}

//This function changes the root of a rtree
char change_root(unsigned int new_root,t_rtree *rtree)
{
unsigned int _i, anchor_id, next_anchor_id, prev_anchor_id;
t_edge next_edge, edge;

if (new_root>=rtree->nrbranch) { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
anchor_id=new_root, next_anchor_id=*rtree->rbranch[new_root].rbranch;
edge.vertice[0]=rtree->rbranch[anchor_id].edge.vertice[0], edge.vertice[1]=rtree->rbranch[anchor_id].edge.vertice[1], next_edge.vertice[0]=rtree->rbranch[next_anchor_id].edge.vertice[0], next_edge.vertice[1]=rtree->rbranch[next_anchor_id].edge.vertice[1];
while (anchor_id!=rtree->root)
  {
  //Swap the master <-> slave anchor pair
  prev_anchor_id=anchor_id, anchor_id=next_anchor_id, next_anchor_id=*rtree->rbranch[anchor_id].rbranch;
  rtree->rbranch[anchor_id].edge.vertice[0]=edge.vertice[1], rtree->rbranch[anchor_id].edge.vertice[1]=edge.vertice[0]; 
  edge.vertice[0]=next_edge.vertice[0], edge.vertice[1]=next_edge.vertice[1];
  next_edge.vertice[0]=rtree->rbranch[next_anchor_id].edge.vertice[0], next_edge.vertice[1]=rtree->rbranch[next_anchor_id].edge.vertice[1];
  if (*rtree->rbranch[anchor_id].rbranch!=prev_anchor_id)
    {
    _i=rtree->rbranch[anchor_id].nrbranch; do { if (!(_i--)) { ylib_errno=YERROR_INTERNAL_CODE; free(rtree); return FALSE; } } while (rtree->rbranch[anchor_id].rbranch[_i]!=prev_anchor_id);
    rtree->rbranch[anchor_id].rbranch[_i]=*rtree->rbranch[anchor_id].rbranch, *rtree->rbranch[anchor_id].rbranch=prev_anchor_id;
    }
  }
rtree->root=new_root;
return TRUE;
}

//This function builds super-anchors stucture from given rtree and anchors list
t_clist *build_sanchors(unsigned int natoms,unsigned int *anchor_id,t_clist *anchors,t_rtree *rtree)
{
register unsigned int _i, _j, _k, _l, _n, _t, count;
unsigned int *data, *sanchor_id;
t_clist *sanchors;
//Stage I. Allocate memory for data and number anchors
if (!(data=(unsigned int*)malloc(sizeof(unsigned int)*anchors->size*2))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { _l=anchors->size, sanchor_id=data+_l; while (_l--) sanchor_id[_l]=(unsigned int)-1; } //Init sanchors
data[rtree->root]=1;                                              //Enumerate data
rtree->rbranch[rtree->root].edge.type=rtree->rbranch[rtree->root].nrbranch;
while (rtree->rbranch[rtree->root].edge.type--)
  {
  _i=rtree->rbranch[rtree->root].rbranch[rtree->rbranch[rtree->root].edge.type];
  data[_i]=2, rtree->rbranch[_i].edge.type=rtree->rbranch[_i].nrbranch;
  do{
    while(--rtree->rbranch[_i].edge.type)
      {
      _j=rtree->rbranch[_i].rbranch[rtree->rbranch[_i].edge.type];
      rtree->rbranch[_j].edge.type=rtree->rbranch[_j].nrbranch;
      data[_j]=data[_i]+1, _i=_j;
      }
    }while ((_i=*rtree->rbranch[_i].rbranch)!=rtree->root);
  }
//Stage II. Aggregare anchors
count=0, _t=anchors->size;
while (_t--)
  {
  _n=rtree->rbranch[_t].cross;
  while (_n--)
    if (data[_t]<data[_j=rtree->rbranch[_t].rbranch[rtree->rbranch[_t].nrbranch+_n]]) 
      {//Skip double calls at the same cross
      _i=_t;
      if (sanchor_id[_i]==(unsigned int)-1) sanchor_id[_i]=count++; 
      do{
        if (data[_j]<data[_i]) { _l=_i, _i=_j, _j=_l; } //swap i and j
        if (sanchor_id[_j]!=(unsigned int)-1)
          {//sanchor_id[_j]!=(unsigned int)-1
               if (sanchor_id[_i]<sanchor_id[_j]) 
                 {
                 count--, _k=sanchor_id[_j], _l=anchors->size; while (_l--) { if (sanchor_id[_l]==_k) sanchor_id[_l]=sanchor_id[_i]; else if (sanchor_id[_l]==count) sanchor_id[_l]=_k; } 
                 }
          else if (sanchor_id[_i]>sanchor_id[_j])
                 {
                 count--, _k=sanchor_id[_i], _l=anchors->size; while (_l--) { if (sanchor_id[_l]==_k) sanchor_id[_l]=sanchor_id[_j]; else if (sanchor_id[_l]==count) sanchor_id[_l]=_k; } 
                 }
          }
        else 
          sanchor_id[_j]=sanchor_id[_i]; //sanchor_id[_j]==(unsigned int)-1 
        }while ((_j=*rtree->rbranch[_j].rbranch)!=_i);
      }
  }
//Stage III. Form atomic aggregates
_l=anchors->size; while (_l--) if (sanchor_id[_l]==(unsigned int)-1) sanchor_id[_l]=count++;
if (!(sanchors=alloc_clist(count,natoms))) { free(data); return FALSE; } else { sanchors->size=_l=count; while (_l--) sanchors->list[_l].size=0; }
_l=anchors->size; while (_l--) sanchors->list[sanchor_id[_l]].size+=anchors->list[_l].size;
for (sanchors->list[0].list=sanchors->_items, _l=1; _l<count; _l++) sanchors->list[_l].list=sanchors->list[_l-1].list+sanchors->list[_l-1].size;
_l=count; while (_l--) sanchors->list[_l].size=0;
_l=anchors->size; while (_l--) { memcpy(&sanchors->list[sanchor_id[_l]].list[sanchors->list[sanchor_id[_l]].size],anchors->list[_l].list,sizeof(unsigned int)*anchors->list[_l].size), sanchors->list[sanchor_id[_l]].size+=anchors->list[_l].size; }
//Free some memory and exit
free(data);
return sanchors;
}

//This function computes supertree (i.e. joins anchors so there is no crosses)
t_rtree *build_srtree(unsigned int sroot,unsigned int natoms,t_clist *sanchors,unsigned int nedges,t_edge *edges,t_rtree *rtree)
{
register unsigned int _i, _j;
unsigned int *sanchor_id;
t_rtree *srtree;

//Stage I. Create sanchor_id;
if (!(sanchor_id=(unsigned int*)malloc(sizeof(unsigned int)*natoms))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { _i=sanchors->size; while (_i--) { _j=sanchors->list[_i].size; while (_j--) sanchor_id[sanchors->list[_i].list[_j]]=_i; } }
//Stage II. Built the supertree
srtree=build_rtree(sroot,sanchor_id,sanchors,nedges,edges);
//Exit
free(sanchor_id);
return srtree;
}

/*
//This function create superroot for given rtree and active_a list.
//It cleans the active list from cycles and crosses.
char compile_superroot(unsigned int *nrbranch,unsigned int **rbranch,t_list *active_a,t_rtree *rtree) 
{
unsigned int _i, _j, _k, _l, _n;

//Stage 1. Alloc memory
if ( (!(active_a))||(!(active_a->size))||(active_a->size==rtree->nrbranch) ) { *nrbranch=0x0, *rbranch=0x0; return TRUE; }
//Mark probably active anchors as '-'
_i=active_a->size;
while (_i--) 
       if ( (rtree->rbranch[active_a->list[_i]].cross)) active_a->list[_i]=active_a->list[--active_a->size]; //Remove crosses
  else if (((int)rtree->rbranch[active_a->list[_i]].aedge.type)<0) active_a->list[_i]=active_a->list[--active_a->size]; //Remove duplicates from active list
  else   rtree->rbranch[active_a->list[_i]].aedge.type=-(int)rtree->rbranch[active_a->list[_i]].aedge.type; //Mark active anchor 
if (active_a->size==rtree->nrbranch) { free(*rbranch); return NTNF; } //Entire rtree is active

//Stage 2. Check if the rtree->root is in static fragment of a tree and change the root if it is not
if ((int)rtree->rbranch[rtree->root].aedge.type<0)
  {//A root must be moved to the static layer
  _i=0; while ( ((int)rtree->rbranch[_i].aedge.type<0)) _i++; //A new root choosen
  _j=_i, _i=*rtree->rbranch[_i].rbranch;
  while (_j!=rtree->root)
    {//Invert path to previous root
    _l=rtree->rbranch[_i].nrbranch; do { _l--; } while(rtree->rbranch[_i].rbranch[_l]!=_j);
    rtree->rbranch[_i].rbranch[_l]=*rtree->rbranch[_i].rbranch, *rtree->rbranch[_i].rbranch=_j;
    _j=_i, _i=rtree->rbranch[_i].rbranch[_l];
    }
  }

//Stage 3. Setup superroot. Our strategy uses benefits of static root position so flexible fragments are downstream of it in the rtree. So we chosing promising vertice and try to scann a subtree downstream of it.
for (_k=0; _k<active_a->size; _k++) 
  {//Scann rtree
  _l=rtree->rbranch[active_a->list[_k]].edge.type=rtree->rbranch[active_a->list[_k]].nrbranch;
  while (--_l)  
    if ((int)rtree->rbranch[rtree->rbranch[active_a->list[_k]].rbranch[_l]].aedge.type>0) //neighbor is static, delete the candidate
      {
      DEL_K:;
      rtree->rbranch[active_a->list[_k]].aedge.type=-(int)rtree->rbranch[active_a->list[_k]].aedge.type;
      active_a->list[_k--]=active_a->list[--active_a->size]; 
      goto NEXT_K;
      }
  while (--rtree->rbranch[active_a->list[_k]].edge.type)
    if (find_in_row(rtree->rbranch[active_a->list[_k]].rbranch[rtree->rbranch[active_a->list[_k]].edge.type],_k,active_a->list)==(unsigned int)-1)
      {
      _i=rtree->rbranch[active_a->list[_k]].rbranch[rtree->rbranch[active_a->list[_k]].edge.type];
      _l=rtree->rbranch[_i].edge.type=rtree->rbranch[_i].nrbranch;
      while (_l--)
        if ((int)rtree->rbranch[rtree->rbranch[_i].rbranch[_l]].aedge.type>0)
           {//neighbor is static, current and previous 
           rtree->rbranch[_i].aedge.type=-(int)rtree->rbranch[_i].aedge.type;
           _j=active_a->size; do { _j--; } while (active_a->list[_j]!=_i); active_a->list[_j]=active_a->list[--active_a->size]; 
           goto DEL_K;
           }
      //Continue scann (if necessary)
      do{
        while (--rtree->rbranch[_i].edge.type)
          if (find_in_row(rtree->rbranch[_i].rbranch[rtree->rbranch[_i].edge.type],_k,active_a->list)==(unsigned int)-1)
            {
            _i=rtree->rbranch[_i].rbranch[rtree->rbranch[_i].edge.type];
            rtree->rbranch[_i].edge.type=rtree->rbranch[_i].nrbranch;
            //Check downstream
            _l=rtree->rbranch[_i].nrbranch;
            while (--_l)  
              if ((int)rtree->rbranch[rtree->rbranch[_i].rbranch[_l]].aedge.type>0) 
                {//neighbor is static, delete entire route 
                do{
                  rtree->rbranch[_i].aedge.type=-(int)rtree->rbranch[_i].aedge.type;
                  _j=active_a->size; do { _j--; } while (active_a->list[_j]!=_i); active_a->list[_j]=active_a->list[--active_a->size]; 
                  _i=*rtree->rbranch[_i].rbranch;
                  } while (_i!=active_a->list[_k]);
                goto DEL_K;
                }
            }
        _i=*rtree->rbranch[_i].rbranch;
        }while (_i!=active_a->list[_k]);
      }
   
  //The subtree is active, move it into next to _k zone
  _n=1,  rtree->rbranch[active_a->list[_k]].edge.type=rtree->rbranch[active_a->list[_k]].nrbranch;
  while (--rtree->rbranch[active_a->list[_k]].edge.type)
    if (find_in_row(rtree->rbranch[active_a->list[_k]].rbranch[rtree->rbranch[active_a->list[_k]].edge.type],_k,active_a->list)==(unsigned int)-1)
      {
      _i=rtree->rbranch[active_a->list[_k]].rbranch[rtree->rbranch[active_a->list[_k]].edge.type];
      rtree->rbranch[_i].edge.type=rtree->rbranch[_i].nrbranch;
      _j=active_a->size; do { _j--; } while (active_a->list[_j]!=_i); 
      if (_j!=_k+_n) { _l=active_a->list[_j], active_a->list[_j]=active_a->list[_k+_n], active_a->list[_k+_n]=_l; } 
      _n++;
      do{
        while (--rtree->rbranch[_i].edge.type)
          if (find_in_row(rtree->rbranch[_i].rbranch[rtree->rbranch[_i].edge.type],_k,active_a->list)==(unsigned int)-1)
            {
            _i=rtree->rbranch[_i].rbranch[rtree->rbranch[_i].edge.type];
            rtree->rbranch[_i].edge.type=rtree->rbranch[_i].nrbranch;
            _j=active_a->size; do { _j--; } while (active_a->list[_j]!=_i); 
            if (_j!=_k+_n) { _l=active_a->list[_j], active_a->list[_j]=active_a->list[_k+_n], active_a->list[_k+_n]=_l; } 
            _n++;
            }
        _i=*rtree->rbranch[_i].rbranch;
        }while (_i!=active_a->list[_k]);
      }
  _k+=_n-1;
  NEXT_K: ;
  } 

//Stage 4. Sync memory and exit
if (!active_a->size) { *nrbranch=0, *rbranch=0x0; return NTNF; } //Could happens, no active layer at all...
for (_j=_i=0; _i<active_a->size; _i++)
  if (rtree->rbranch[*rtree->rbranch[active_a->list[_i]].rbranch].edge.type>0) 
    {
    if (_j!=_i) { _k=active_a->list[_j], active_a->list[_j]=active_a->list[_i], active_a->list[_i]=_k; }
    _j++;
    }
if (!((*rbranch)=(unsigned int*)malloc(sizeof(unsigned int)*_j))) { ylib_errno=YERROR_MEMORY; return FALSE; }
*nrbranch=_j, _k=active_a->size; 
while (_j!=_k--) { _i=active_a->list[_k], rtree->rbranch[_i].aedge.type=-(int)rtree->rbranch[_i].aedge.type; } //Unmark active layer
do {               _i=active_a->list[_k], rtree->rbranch[_i].aedge.type=-(int)rtree->rbranch[_i].aedge.type, (*rbranch)[_k]=active_a->list[_k]; } while (_k--); //Copy & unmark superroot

return TRUE;
}
*/

/*************************************   M I S S C E L L A N E O U S    P A R T    **********************************/


//This function gather an achors of given residues list
t_list *get_ress_anchors(t_list *ress,t_mol *mol)
{
register unsigned int _i, _j, _k;
unsigned int *anchors, nanchors;
t_list *active_a;
void *vp;

//Stage 1. Gather all the anchors from list (with step 0xFF)
if (!(anchors=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
nanchors=0;
_k=ress->size;
while (_k--)
  {
  if (ress->list[_k]>mol->ress->size) { ylib_errno=YERROR_USER; return FALSE; }
  _i=mol->anchors->size;
  while (_i--)
    {
    _j=mol->anchors->list[_i].size;
    while (_j--)
      if ( (mol->anchors->list[_i].list[_j]>=mol->rsize[ress->list[_k]])&&(mol->anchors->list[_i].list[_j]<mol->rsize[ress->list[_k]+1])&&(find_in_row(_i,nanchors,anchors)==(unsigned int)-1) )
        {
        anchors[nanchors]=_i; //Note. Dublicates are not allowed as no overlaping anchors are accepted
        if (++nanchors%0xFF)
          {
          if (!(vp=realloc(anchors,sizeof(unsigned int)*(nanchors+0xFF)))) { free(anchors); goto LABEL_MEMORY_ERROR; }
          else anchors=(unsigned int*)vp;
          }
        }
    }
  }
//Stage 2. gnerate list
if (!(active_a=alloc_list(nanchors))) { free(anchors); goto LABEL_MEMORY_ERROR; }
memcpy(active_a->list,anchors,sizeof(unsigned int)*nanchors);
free (anchors);
return active_a;
}


