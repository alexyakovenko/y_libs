#include "y_ffsys.h"

extern int ylib_errno;

const double LG_SERIES[]={ //9 items
                         1./EXP_G/EXP_G/EXP_G/EXP_G/EXP_G/EXP_G/EXP_G/EXP_G,
                         1./EXP_G/EXP_G/EXP_G/EXP_G/EXP_G/EXP_G,
                         1./EXP_G/EXP_G/EXP_G/EXP_G,
                         1./EXP_G/EXP_G,
                         1.,
                            EXP_G*EXP_G,
                            EXP_G*EXP_G*EXP_G*EXP_G,
                            EXP_G*EXP_G*EXP_G*EXP_G*EXP_G*EXP_G,
                            EXP_G*EXP_G*EXP_G*EXP_G*EXP_G*EXP_G*EXP_G*EXP_G,
                         }; //10^-8...10^+8
#define E_SERIES LG_SERIES


//This function reads activem from hdd
char read_activem(FILE *in,t_activem *activem)
{ 
return (char) ( ( (activem->active_a=read_list(in)))&&( (activem->active_b=read_list(in)))&&( (activem->active_g=read_list(in)))&&
                ( (activem->active_i=read_list(in)))&&( (activem->active_t=read_list(in)))&&( (activem->active_p=read_list(in)))&&
                (fread(&activem->size_p,sizeof(unsigned int),0x1,in)==0x1) );
}

//This function writes activem to hdd
char write_activem(FILE *out,unsigned int n,t_activem *activem)
{
unsigned int _i;
unsigned int  j;
_i=0;
while (n--)
  {
  if ( (!activem[_i].active_a)||(!activem[_i].active_a->size) ) 
    {
    j=0x0;
    if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) return FALSE;
    }	
  else
    {
    j=Y_MAGIC; 
    if ( (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1)||	
         (!write_list(out,activem[_i].active_a)) ||(!write_list(out,activem[_i].active_b))||(!write_list(out,activem[_i].active_g))||
         (!write_list(out,activem[_i].active_i)) ||(!write_list(out,activem[_i].active_t))||(!write_list(out,activem[_i].active_p))||
         (fwrite(&activem[_i].size_p,sizeof(unsigned int),0x1,out)!=0x1) ) 
      return FALSE;
    }
  _i++;
  }
return TRUE;
}

//This function compile activem block
//Note. It DO NOT copy active_a list, just a pointer to it.
char create_activem(t_activem *activem,t_list *active_a,t_mol *mol)
{
unsigned int _k, _l;
//Init activem
activem->active_a=activem->active_b=activem->active_g=activem->active_i=activem->active_t=activem->active_p=0x0;
if ( ( (active_a))&&( (active_a->size)) )
  {
  //Copy active_a list
  activem->active_a=active_a;
  //Mark active atoms
  _k=active_a->size; while (_k--) { _l=mol->anchors->list[active_a->list[_k]].size; while (_l--) { mol->atoms->list[mol->anchors->list[active_a->list[_k]].list[_l]]=(unsigned int)-mol->atoms->list[mol->anchors->list[active_a->list[_k]].list[_l]]; } }
  //Create active_b, active_g, active_i and active_t lists
  _l=0;_k=mol->size_b; 
  while (_k--) 
    if ( ((int)mol->atoms->list[mol->ff_b[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_b[_k].atom[1]]<0) ) _l++;
  if (_l)
    {
    if (!(activem->active_b=alloc_list(_l))) { MEMORY_ERROR_0:  activem->active_a=0x0, ylib_errno=YERROR_MEMORY, _l=mol->atoms->size; while (_l--) if ((int)mol->atoms->list[_l]<0) mol->atoms->list[_l]=abs(mol->atoms->list[_k]); return FALSE; }
    _l=0, _k=mol->size_b; 
    while (_k--)
      if ( ((int)mol->atoms->list[mol->ff_b[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_b[_k].atom[1]]<0) ) 
        activem->active_b->list[_l++]=_k;
    _l=0, _k=mol->size_g; 
    while (_k--)
      if ( ((int)mol->atoms->list[mol->ff_g[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_g[_k].atom[1]]<0)||((int)mol->atoms->list[mol->ff_g[_k].atom[2]]<0) ) _l++;
    if (_l)
      {
      if (!(activem->active_g=alloc_list(_l))) { MEMORY_ERROR_1: free(activem->active_b); activem->active_b=0x0; goto MEMORY_ERROR_0; }
      _l=0, _k=mol->size_g; 
      while (_k--)
        if ( ((int)mol->atoms->list[mol->ff_g[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_g[_k].atom[1]]<0)||((int)mol->atoms->list[mol->ff_g[_k].atom[2]]<0) ) 
          activem->active_g->list[_l++]=_k;
      _l=0, _k=mol->size_i; 
      while (_k--)
        if ( ((int)mol->atoms->list[mol->ff_i[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_i[_k].atom[1]]<0)||((int)mol->atoms->list[mol->ff_i[_k].atom[2]]<0)||((int)mol->atoms->list[mol->ff_i[_k].atom[3]]<0) ) _l++;
      if (_l)  
        {
        if (!(activem->active_i=alloc_list(_l))) { MEMORY_ERROR_2: free(activem->active_g); activem->active_g=0x0; goto MEMORY_ERROR_1; }
        _l=0, _k=mol->size_i; 
        while (_k--)
          if ( ((int)mol->atoms->list[mol->ff_i[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_i[_k].atom[1]]<0)||((int)mol->atoms->list[mol->ff_i[_k].atom[2]]<0)||((int)mol->atoms->list[mol->ff_i[_k].atom[3]]<0) )
            activem->active_i->list[_l++]=_k;
        }
      _l=0, _k=mol->size_t; 
      while (_k--)
        if ( ((int)mol->atoms->list[mol->ff_t[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_t[_k].atom[1]]<0)||((int)mol->atoms->list[mol->ff_t[_k].atom[2]]<0)||((int)mol->atoms->list[mol->ff_t[_k].atom[3]]<0) ) _l++;
      if (_l)
        {
        if (!(activem->active_t=alloc_list(_l))) { MEMORY_ERROR_3: free(activem->active_i); activem->active_i=0x0; goto MEMORY_ERROR_2; }
        _l=0, _k=mol->size_t; 
        while (_k--)
          if ( ((int)mol->atoms->list[mol->ff_t[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_t[_k].atom[1]]<0)||((int)mol->atoms->list[mol->ff_t[_k].atom[2]]<0)||((int)mol->atoms->list[mol->ff_t[_k].atom[3]]<0) )
            activem->active_t->list[_l++]=_k;
        }
      }
    }
  _l=0;_k=mol->size_p; 
  while(_k--)
    if ( ((int)mol->atoms->list[mol->ff_p[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_p[_k].atom[1]]<0) ) _l++;
  if (_l)
    {
    if (!(activem->active_p=alloc_list(_l))) { free(activem->active_t); activem->active_t=0x0; goto MEMORY_ERROR_3; } 	  
    else activem->size_p=0;        
    _l=0, _k=mol->size_p;
    while (_k--)
      if ((int)mol->atoms->list[mol->ff_p[_k].atom[0]]<0)
        {
        if ((int)mol->atoms->list[mol->ff_p[_k].atom[1]]<0)
          {//Active-Active pair
          activem->active_p->list[_l++]=activem->active_p->list[activem->size_p];
          activem->active_p->list[activem->size_p++]=_k;
          }  
        else //Active-Passive pair
          { AP_PAIR: activem->active_p->list[_l++]=_k; } 
        }
      else if ((int)mol->atoms->list[mol->ff_p[_k].atom[1]]<0) goto AP_PAIR;  
    }   
  //Unmark all
  _k=mol->atoms->size; while (_k--) if ((int)mol->atoms->list[_k]<0) mol->atoms->list[_k]=abs(mol->atoms->list[_k]);
  }
return TRUE;
}


//This function removes global root elements from activem
void remove_global_roots_activem(char bonded,char nonbonded,t_activem *activem,t_rtree *rtree,t_mol *mol)
{
unsigned int _id, _k, _l;
if ( (activem)&&(activem->active_a)&&(mol)&&(rtree)&&(rtree->rbranch[rtree->nrbranch].nrbranch) )
  {
  //Mark all
  _k=activem->active_a->size; while (_k--) { _l=mol->anchors->list[activem->active_a->list[_k]].size; while (_l--) { _id=mol->anchors->list[activem->active_a->list[_k]].list[_l]; mol->atoms->list[_id]=(unsigned int)-mol->atoms->list[_id]; } }

  //Pass nonbonded
  if (nonbonded)
    {
    _k=rtree->rbranch[rtree->nrbranch].nrbranch;
    while (_k--)
      {
      _id=rtree->rbranch[rtree->rbranch[rtree->nrbranch].rbranch[_k]].edge.vertice[1];
      _l=activem->active_p->size;
      while (_l--)
        if ( ( (mol->ff_p[activem->active_p->list[_l]].atom[0]==_id)&&((int)mol->atoms->list[mol->ff_p[activem->active_p->list[_l]].atom[1]]>0) )||
             ( (mol->ff_p[activem->active_p->list[_l]].atom[1]==_id)&&((int)mol->atoms->list[mol->ff_p[activem->active_p->list[_l]].atom[0]]>0) ) )
          activem->active_p->list[_l]=activem->active_p->list[--activem->active_p->size];
      }
    }
  //Pass bonded
  if (bonded)
    {
    _k=rtree->rbranch[rtree->nrbranch].nrbranch;
    while (_k--)
      {
      _id=rtree->rbranch[rtree->rbranch[rtree->nrbranch].rbranch[_k]].edge.vertice[1];
      _l=activem->active_b->size;
      while (_l--)
        if ( ( (mol->ff_b[activem->active_b->list[_l]].atom[0]==_id)&&((int)mol->atoms->list[mol->ff_b[activem->active_b->list[_l]].atom[1]]>0) )||
             ( (mol->ff_b[activem->active_b->list[_l]].atom[1]==_id)&&((int)mol->atoms->list[mol->ff_b[activem->active_b->list[_l]].atom[0]]>0) ) )
          activem->active_b->list[_l]=activem->active_b->list[--activem->active_b->size];
      _l=activem->active_g->size;
      while (_l--)
        if ( ( (mol->ff_g[activem->active_g->list[_l]].atom[0]==_id)&&((int)mol->atoms->list[mol->ff_g[activem->active_g->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_g[activem->active_g->list[_l]].atom[2]]>0) ) ||
             ( (mol->ff_g[activem->active_g->list[_l]].atom[1]==_id)&&((int)mol->atoms->list[mol->ff_g[activem->active_g->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_g[activem->active_g->list[_l]].atom[2]]>0) ) ||
             ( (mol->ff_g[activem->active_g->list[_l]].atom[2]==_id)&&((int)mol->atoms->list[mol->ff_g[activem->active_g->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_g[activem->active_g->list[_l]].atom[1]]>0) )  )
          activem->active_g->list[_l]=activem->active_g->list[--activem->active_g->size];
      _l=activem->active_i->size;
      while (_l--)
        if ( ( (mol->ff_i[activem->active_i->list[_l]].atom[0]==_id)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[2]]>0)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[3]]>0) ) ||
             ( (mol->ff_i[activem->active_i->list[_l]].atom[1]==_id)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[2]]>0)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[3]]>0) ) ||
             ( (mol->ff_i[activem->active_i->list[_l]].atom[2]==_id)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[3]]>0) ) ||
             ( (mol->ff_i[activem->active_i->list[_l]].atom[3]==_id)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_i[activem->active_i->list[_l]].atom[2]]>0) )  )
          activem->active_i->list[_l]=activem->active_i->list[--activem->active_i->size];
      _l=activem->active_t->size;
      while (_l--)
        if ( ( (mol->ff_t[activem->active_t->list[_l]].atom[0]==_id)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[2]]>0)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[3]]>0) ) ||
             ( (mol->ff_t[activem->active_t->list[_l]].atom[1]==_id)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[2]]>0)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[3]]>0) ) ||
             ( (mol->ff_t[activem->active_t->list[_l]].atom[2]==_id)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[3]]>0) ) ||
             ( (mol->ff_t[activem->active_t->list[_l]].atom[3]==_id)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_t[activem->active_t->list[_l]].atom[2]]>0) )  )
          activem->active_t->list[_l]=activem->active_t->list[--activem->active_t->size];
      }
    }
  //Unmark all
  _k=activem->active_a->size; while (_k--) { _l=mol->anchors->list[activem->active_a->list[_k]].size; while (_l--) { _id=mol->anchors->list[activem->active_a->list[_k]].list[_l]; mol->atoms->list[_id]=(unsigned int)-mol->atoms->list[_id]; } }
  }
}

//This function deletes activem lists
void free_activem(t_activem *activem)
{
if (activem->active_a) { free(activem->active_a); activem->active_a=0x0; }
if (activem->active_b) { free(activem->active_b); activem->active_b=0x0; }
if (activem->active_g) { free(activem->active_g); activem->active_g=0x0; }
if (activem->active_i) { free(activem->active_i); activem->active_i=0x0; }
if (activem->active_t) { free(activem->active_t); activem->active_t=0x0; }
if (activem->active_p) { free(activem->active_p); activem->active_p=0x0; }
}
 
//     F F S Y S       C O M P I L A T O R     P A R T

//This function empty ffsys data structure
void free_ffsys(t_ffsys *ffsys)
{
register unsigned int _i;
if (ffsys)
  {
  _i=ffsys->nmols;
  while (_i--)
    {
    if ( (!_i)||(ffsys->mols[_i]!=ffsys->mols[_i-1]) ) { free_mol(ffsys->mols[_i]); ffsys->mols[_i]=0x0; }
    if ( (!_i)||(ffsys->activem[_i].active_a!=ffsys->activem[_i-1].active_a) )
      {
      free_activem(&ffsys->activem[_i]);
      if (ffsys->rtrees[_i]) { free(ffsys->rtrees[_i]); ffsys->rtrees[_i]=0x0; }
      }
   }
  if (ffsys->mols)    { free(ffsys->mols);    ffsys->mols=0x0;    }
  if (ffsys->nr)      { free(ffsys->nr);      ffsys->nr=0x0;      }
  if (ffsys->activem) { free(ffsys->activem); ffsys->activem=0x0; }
  if (ffsys->rtrees)  { free(ffsys->rtrees);  ffsys->rtrees=0x0;  }
  if (ffsys->a)       { free(ffsys->a);       ffsys->a=0x0;       }
  if (ffsys->q)       { free(ffsys->q);       ffsys->q=0x0;       }
  if (ffsys->g)       { free(ffsys->g);       ffsys->g=0x0;       }
  if (ffsys->_g)      { free(ffsys->_g);      ffsys->_g=0x0;      }
  if (ffsys->r)       { free(ffsys->r);       ffsys->r=0x0;       }
  free(ffsys);
  }
}

//This function reads ffsys from hdd
t_ffsys *read_ffsys(FILE *in)
{
t_ffsys *ffsys=0x0;
unsigned int i, j;
//Read general data
if ( (fread(&i,sizeof(unsigned int),0x1,in)!=0x1)||(i!=Y_MAGIC) ) { _LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
if (!(ffsys=(t_ffsys*)calloc(sizeof(t_ffsys),0x1))) { _LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
if ( (fread(&ffsys->nmols,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&ffsys->natoms,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&ffsys->naatoms,sizeof(unsigned int),0x1,in)!=0x1) ) { LABEL_IO_ERROR: free_ffsys(ffsys); goto _LABEL_IO_ERROR; }
if (!(ffsys->mols=(t_mol**)malloc(sizeof(t_mol*)*ffsys->nmols))) goto _LABEL_MEMORY_ERROR;
//Read mols data
for (i=0;i<ffsys->nmols;i++)
  {
  if ( (!(ffsys->mols[i]=read_ymol(in)))||(fread(&j,sizeof(unsigned int),0x1,in)!=0x1) ) goto LABEL_IO_ERROR;
  if (i+j>ffsys->nmols) goto LABEL_IO_ERROR;
  while (--j) { ffsys->mols[i+1]=ffsys->mols[i], i++; }
  }
if ( (!(ffsys->activem=(t_activem*)calloc(sizeof(t_activem),ffsys->nmols)))||(!(ffsys->rtrees=(t_rtree**)calloc(sizeof(t_rtree*),ffsys->nmols)))||
     (!(ffsys->nr=(unsigned int*)malloc(sizeof(unsigned int)*ffsys->nmols))) ) { LABEL_MEMORY_ERROR: free_ffsys(ffsys); goto _LABEL_MEMORY_ERROR; }
for (i=0;i<ffsys->nmols;i++)
  {
  if (fread(&j,sizeof(unsigned int),0x1,in)!=0x1) goto LABEL_IO_ERROR;
  if (j)
    {
    if (j==Y_MAGIC)
      {
      if ( (!(ffsys->activem[i].active_a=read_list(in)))||(!(ffsys->activem[i].active_b=read_list(in)))||(!(ffsys->activem[i].active_g=read_list(in)))||
           (!(ffsys->activem[i].active_i=read_list(in)))||(!(ffsys->activem[i].active_t=read_list(in)))||(!(ffsys->activem[i].active_p=read_list(in))) ) goto LABEL_IO_ERROR;
      }
    else goto LABEL_IO_ERROR;
    }
  if (fread(&j,sizeof(unsigned int),0x1,in)!=0x1) goto LABEL_IO_ERROR;
  if (j)
    {
    if (j==Y_MAGIC)
      {
      if (!(ffsys->rtrees[i]=read_rtree(in))) return FALSE;
      }
    else goto LABEL_IO_ERROR;
    }
  }
if (fread(ffsys->nr,sizeof(unsigned int),ffsys->nmols,in)!=ffsys->nmols) goto LABEL_IO_ERROR;
//Read global massives
if ( (!(ffsys->a=(char*)malloc(sizeof(char)*ffsys->natoms)))||(!(ffsys->q=(double*)malloc(sizeof(double)*ffsys->natoms)))||
     (!(ffsys->g=(t_vec*)malloc(sizeof(t_vec)*ffsys->natoms)))||(!(ffsys->_g=(t_vec (*)[N_ESLICES])malloc(N_ESLICES*sizeof(t_vec)*ffsys->natoms)))||
     (!(ffsys->r=(t_vec*)malloc(sizeof(t_vec)*ffsys->natoms))) ) goto LABEL_MEMORY_ERROR;
if ( (fread(ffsys->a,sizeof(char),ffsys->natoms,in)!=ffsys->natoms)||(fread(ffsys->q,sizeof(double),ffsys->natoms,in)!=ffsys->natoms)||(fread(ffsys->r,sizeof(t_vec),ffsys->natoms,in)!=ffsys->natoms) ) goto LABEL_IO_ERROR;
return ffsys;
}


//This function writes ffsys to hdd
char write_ffsys(FILE *out,t_ffsys *ffsys)
{
unsigned int i,j;
//write general data
i=Y_MAGIC;
if ( (fwrite(&i,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&ffsys->nmols,sizeof(unsigned int),0x1,out)!=0x1)||
     (fwrite(&ffsys->natoms,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&ffsys->naatoms,sizeof(unsigned int),0x1,out)!=0x1) ) { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
//write mol data
for (j=1, i=0; i<ffsys->nmols; i++, j++)
  {
  if ( (i==ffsys->nmols)||(ffsys->mols[i]!=ffsys->mols[i+1]) )
    {
    if (!(write_ymol(out,ffsys->mols[i]))) return FALSE;
    if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;
    else j=0;
    }	
  }
for (i=0;i<ffsys->nmols;i++)
  {
  if (ffsys->activem[i].active_a)
    {
    j=Y_MAGIC;
    if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;	
    if ( (!write_list(out,ffsys->activem[i].active_a))||(!write_list(out,ffsys->activem[i].active_b))||(!write_list(out,ffsys->activem[i].active_g))||
         (!write_list(out,ffsys->activem[i].active_i))||(!write_list(out,ffsys->activem[i].active_t))||(!write_list(out,ffsys->activem[i].active_p)) ) 
      goto LABEL_IO_ERROR;
    }
  else
    {  
    j=0;
    if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;  
    }
  if (ffsys->rtrees[i]) 
    {
    j=Y_MAGIC;
    if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;	
    if (!write_rtree(out,ffsys->rtrees[i])) goto LABEL_IO_ERROR;
    }
  else	
    {  
    j=0;
    if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;  
    }
  }
if (fwrite(ffsys->nr,sizeof(unsigned int),ffsys->nmols,out)!=ffsys->nmols) goto LABEL_IO_ERROR;
//write global massives  
if ( (fwrite(ffsys->a,sizeof(char),ffsys->natoms,out)!=ffsys->natoms)||(fwrite(ffsys->q,sizeof(double),ffsys->natoms,out)!=ffsys->natoms)||(fwrite(ffsys->r,sizeof(t_vec),ffsys->natoms,out)!=ffsys->natoms) ) goto LABEL_IO_ERROR;
return TRUE;
}


//This function compile individual mol into ffsys after allocing memory for it
//It DO NOT copy rvecs
char ffsys_add_mol(unsigned int nmols,t_mol *mol,t_activem *activem,t_rtree *rtree,t_ffsys *ffsys)
{
unsigned int _j, _k;
//Book memory
if (ffsys->g)  { free(ffsys->g);  ffsys->g=0x0;  }
if (ffsys->_g) { free(ffsys->_g); ffsys->_g=0x0; }
if ( (!(ffsys->mols=(t_mol**)realloc(ffsys->mols,sizeof(t_mol)*(ffsys->nmols+nmols))))                          ||
     (!(ffsys->nr=(unsigned int*)realloc(ffsys->nr,sizeof(unsigned int)*(ffsys->nmols+nmols))))                 || 
     (!(ffsys->rtrees=(t_rtree**)realloc(ffsys->rtrees,sizeof(t_rtree*)*(ffsys->nmols+nmols))))                 ||
     (!(ffsys->activem=(t_activem*)realloc(ffsys->activem,sizeof(t_activem)*(ffsys->nmols+nmols))))             ||
     (!(ffsys->r=(t_vec*)realloc(ffsys->r,sizeof(t_vec)*(ffsys->natoms+nmols*mol->atoms->size))))               ||
     (!(ffsys->g=(t_vec*)malloc((sizeof(t_vec)*(ffsys->natoms+nmols*mol->atoms->size)))))                       ||
     (!(ffsys->_g=(t_vec (*)[N_ESLICES])malloc(N_ESLICES*sizeof(t_vec)*(ffsys->natoms+nmols*mol->atoms->size))))||
     (!(ffsys->q=(double*)realloc(ffsys->q,sizeof(double)*(ffsys->natoms+nmols*mol->atoms->size))))             || 
     (!(ffsys->a=(char*)realloc(ffsys->a,sizeof(char)*(ffsys->natoms+nmols*mol->atoms->size))))                  ) { ylib_errno=YERROR_MEMORY; return FALSE; }

//Mark active components and create active_a lists
for (_j=0;_j<nmols;_j++,ffsys->nmols++)
  {
  ffsys->mols[ffsys->nmols]=mol, ffsys->nr[ffsys->nmols]=ffsys->natoms;
  for (_k=0;_k<mol->atoms->size;_k++) 
    {
    ffsys->a[ffsys->natoms]=(unsigned char)mol->atoms->list[_k];
    ffsys->q[ffsys->natoms]=mol->charges[_k];
    ffsys->natoms++;
    }	
  if ( ( (activem))&&( (activem->active_a)) )
    {	
    _k=activem->active_a->size; while (_k--) ffsys->naatoms+=mol->anchors->list[activem->active_a->list[_k]].size;
    ffsys->activem[ffsys->nmols].active_a=activem->active_a;
    ffsys->activem[ffsys->nmols].active_b=activem->active_b;
    ffsys->activem[ffsys->nmols].active_g=activem->active_g;
    ffsys->activem[ffsys->nmols].active_i=activem->active_i;
    ffsys->activem[ffsys->nmols].active_t=activem->active_t;
    ffsys->activem[ffsys->nmols].active_p=activem->active_p;
    _k=activem->active_a->size; while (_k--) ffsys->naatoms+=mol->anchors->list[activem->active_a->list[_k]].size; 
    }
  else
    {//Whole molecule is active
    memset(&ffsys->activem[ffsys->nmols],0x0,sizeof(t_activem));
    ffsys->naatoms+=mol->atoms->size;
    }
  //Setup rtree
  ffsys->rtrees[ffsys->nmols]=rtree; 
  }
return TRUE;
}

//This function add water molucules template to ffsys (atom's coords are uninitialized)
char ffsys_add_sol(unsigned int nsols,t_top *top,t_ffsys *ffsys)
{
unsigned int _j, _k;
t_mol *sol=0x0;
t_rtree *sol_rtree=0x0;
t_clist *neighbors=0x0;

//Stage I. Create sol template 
//Stage I.1. Create chemical structure
if ( (!(sol=calloc(sizeof(t_mol),1)))||(!(sol->name=(char*)malloc(sizeof(char)*4)))                      ||
     (!(sol->ress=(t_list*)alloc_list(1)))||(!(sol->rsize=(unsigned int*)malloc(sizeof(unsigned int)*2)))||
     (!(sol->atoms=(t_list*)alloc_list(3)))||(!(sol->rvecs=(t_vec*)malloc(sizeof(t_vec)*3)))||(!(sol->charges=(double*)malloc(sizeof(double)*3)))||
     (!(sol->edges=(t_edge*)malloc(sizeof(t_edge)*(sol->nedges=2)))) ) 
  {
  ylib_errno=YERROR_MEMORY;
  FREE_SOL:
  if ( (sol->ress))  { free(sol->ress);  sol->ress=0x0; }
  if ( (sol->rsize)) { free(sol->rsize); sol->rsize=0x0; }
  if ( (sol->atoms)) { free(sol->atoms); sol->atoms=0x0; }
  if ( (sol->rvecs)) { free(sol->rvecs); sol->rvecs=0x0; }
  if ( (sol->charges)) { free(sol->charges); sol->charges=0x0; }
  if ( (sol->edges)) { free(sol->edges); sol->edges=0x0; }
  if ( (sol)) { free(sol); sol=0x0; }
  return FALSE;
  }
sol->name[0]='W', sol->name[1]='A', sol->name[2]='T', sol->name[3]='\0';
memcpy(sol->ress->list,sol->name,sizeof(unsigned int));
sol->rsize[0]=0, sol->rsize[1]=sol->atoms->size;
sol->atoms->list[0]=8, sol->atoms->list[1]=-8, sol->atoms->list[2]=-8;
sol->edges[0].vertice[0]=0, sol->edges[0].vertice[1]=1, sol->edges[0].type=1, sol->edges[1].vertice[0]=0, sol->edges[1].vertice[1]=2, sol->edges[1].type=1;
//Stage I.2. Compile sol template
if ( (!(neighbors=define_neighbors(FALSE,3,2,sol->edges)))                                  ||
     (!(compile_mol(neighbors,sol,top)))||(!(ionize_mol_empirically(7.,&neighbors,sol,top)))||
     (!(compose_mol(TRUE,neighbors,sol,top)))||(!(parameterize_mol(sol,top)))                )
  {
  if (!(neighbors)) { ylib_errno=YERROR_MEMORY; goto FREE_SOL; }
  else        { free(neighbors); neighbors=0x0; goto FREE_SOL; }
  }
free(neighbors); neighbors=0x0;
//Stage I.3. Create rtree
if (!(sol_rtree=build_rtree(0x0,sol))) goto FREE_SOL;

//Stage II. Integrate waters into ffsys
if (ffsys->g)  { free(ffsys->g);  ffsys->g=0x0;  }
if (ffsys->_g) { free(ffsys->_g); ffsys->_g=0x0; }
if ( (!(ffsys->mols=(t_mol**)realloc(ffsys->mols,sizeof(t_mol)*(ffsys->nmols+nsols))))                          ||
     (!(ffsys->nr=(unsigned int*)realloc(ffsys->nr,sizeof(unsigned int)*(ffsys->nmols+nsols))))                 || 
     (!(ffsys->rtrees=(t_rtree**)realloc(ffsys->rtrees,sizeof(t_rtree*)*(ffsys->nmols+nsols))))                 ||
     (!(ffsys->activem=(t_activem*)realloc(ffsys->activem,sizeof(t_activem)*(ffsys->nmols+nsols))))             ||
     (!(ffsys->r=(t_vec*)realloc(ffsys->r,sizeof(t_vec)*(ffsys->natoms+nsols*sol->atoms->size))))               ||
     (!(ffsys->g=(t_vec*)malloc(sizeof(t_vec)*(ffsys->natoms+nsols*sol->atoms->size))))                         ||
     (!(ffsys->_g=(t_vec (*)[N_ESLICES])malloc(N_ESLICES*sizeof(t_vec)*(ffsys->natoms+nsols*sol->atoms->size))))||
     (!(ffsys->q=(double*)realloc(ffsys->q,sizeof(double)*(ffsys->natoms+nsols*sol->atoms->size))))             || 
     (!(ffsys->a=(char*)realloc(ffsys->a,sizeof(char)*(ffsys->natoms+nsols*sol->atoms->size))))                  ) { ylib_errno=YERROR_MEMORY; free(sol_rtree); sol_rtree=0x0; goto FREE_SOL; }
//Mark active components and create active_a lists
for (_j=0;_j<nsols;_j++, ffsys->naatoms+=sol->atoms->size, ffsys->nmols++)
  {
  ffsys->mols[ffsys->nmols]=sol, ffsys->nr[ffsys->nmols]=ffsys->natoms;
  for (_k=0;_k<sol->atoms->size;_k++, ffsys->natoms++) 
    {
    ffsys->a[ffsys->natoms]=(unsigned char)sol->atoms->list[_k];
    ffsys->q[ffsys->natoms]=sol->charges[_k];
    }	
  //Whole solvent molecules are active plus handle rtrees
  memset(&ffsys->activem[ffsys->nmols],0x0,sizeof(t_activem));
  ffsys->rtrees[ffsys->nmols]=sol_rtree; 
  }
return TRUE;
}

//-----------------------------   Y F F 1    P A R T   ---------------------------------------------

//-----------------------------   G R I D    P A R T   ---------------------------------------------

//This function calculates values of yff1 A exponent
double inline calc_ffgrid_yff1_A(t_vec *r,t_ffsys *ffsys,t_top *top)
{
register unsigned int _i;
register double _a=0., _rr;
_i=ffsys->natoms;
while (_i--)
  if (ffsys->a[_i]>0) //Avoid marked active atoms
    {
    _rr=sqrd(r->i-ffsys->r[_i].i)+sqrd(r->j-ffsys->r[_i].j)+sqrd(r->k-ffsys->r[_i].k)+SMALL2;
    _a+=top->A[(unsigned int)ffsys->a[_i]][0]/sqrd(_rr*_rr*_rr);
    }
return _a;
}
//This function calculates protomap element for yff1 A exponent
char calc_ffgrid_yff1_protoA(double *f,double *dfdx,double *dfdy,double *dfdz,double *d2fdxdy,double *d2fdxdz,double *d2fdydz,double *d3fdxdydz,double rr,t_vec *r,char label,va_list stack)
{
unsigned int _i;
double _rr,_d;
t_vec _r;
t_top *top;
t_ffsys *ffsys;
va_list _stack;

va_copy(_stack,stack);
ffsys=va_arg(_stack,t_ffsys*);
top=va_arg(_stack,t_top*);
*f=*dfdx=*dfdy=*dfdz=*d2fdxdy=*d2fdxdz=*d2fdydz=*d3fdxdydz=0.;	
_i=ffsys->natoms;
while (_i--)
  if (ffsys->a[_i]>0) //Avoid marked active atoms
    {
    _r.i=r->i-ffsys->r[_i].i, _r.j=r->j-ffsys->r[_i].j, _r.k=r->k-ffsys->r[_i].k;
    if ((_rr=_r.i*_r.i+_r.j*_r.j+_r.k*_r.k)<rr) _rr=0.5*(rr+_rr);
    _rr+=SMALL2;
    _d=top->A[(unsigned int)ffsys->a[_i]][0]/sqrd(_rr*_rr*_rr);
    *f+=_d, _d/=_rr;
    *dfdx-=12.*_d*_r.i,
    *dfdy-=12.*_d*_r.j,
    *dfdz-=12.*_d*_r.k, _d/=_rr;
    *d2fdxdy+=168.*_d*_r.i*_r.j,
    *d2fdxdz+=168.*_d*_r.i*_r.k,
    *d2fdydz+=168.*_d*_r.j*_r.k, _d/=_rr;
    *d3fdxdydz-=2688.*_d*_r.i*_r.j*_r.k;	
    }
va_end(stack);
return TRUE;
}

//This function calculates values of yff1 B exponent
double inline calc_ffgrid_yff1_B(t_vec *r,t_ffsys *ffsys,t_top *top)
{
register unsigned int _i;
register double _b=0., _rr;
_i=ffsys->natoms;
while (_i--)
  if (ffsys->a[_i]>0) //Avoid marked active atoms
    {
    _rr=sqrd(r->i-ffsys->r[_i].i)+sqrd(r->j-ffsys->r[_i].j)+sqrd(r->k-ffsys->r[_i].k)+SMALL2;
    _b+=top->B[(unsigned int)ffsys->a[_i]][0]/_rr/_rr/_rr;
    }
return _b;
}
//This function calculates protomap element for yff1 B exponent
char calc_ffgrid_yff1_protoB(double *f,double *dfdx,double *dfdy,double *dfdz,double *d2fdxdy,double *d2fdxdz,double *d2fdydz,double *d3fdxdydz,double rr,t_vec *r,char label,va_list stack)
{
unsigned int _i;
double _rr,_d;
t_vec _r;
t_top *top;
t_ffsys *ffsys;
va_list _stack;

va_copy(_stack,stack);
ffsys=va_arg(_stack,t_ffsys*);
top=va_arg(_stack,t_top*);
*f=*dfdx=*dfdy=*dfdz=*d2fdxdy=*d2fdxdz=*d2fdydz=*d3fdxdydz=0.;	
_i=ffsys->natoms;
while (_i--)
  {
  if (ffsys->a[_i]>0) //Avoid marked active atoms
    {
    _r.i=r->i-ffsys->r[_i].i, _r.j=r->j-ffsys->r[_i].j, _r.k=r->k-ffsys->r[_i].k;
    if ((_rr=_r.i*_r.i+_r.j*_r.j+_r.k*_r.k)<rr) _rr=0.5*(rr+_rr);
    _rr+=SMALL2;
    _d=top->B[(unsigned int)ffsys->a[_i]][0]/_rr/_rr/_rr;
    *f+=_d, _d/=_rr,
    *dfdx-=6.*_d*_r.i,
    *dfdy-=6.*_d*_r.j,
    *dfdz-=6.*_d*_r.k, _d/=_rr,
    *d2fdxdy+=48.*_d*_r.i*_r.j,
    *d2fdxdz+=48.*_d*_r.i*_r.k,
    *d2fdydz+=48.*_d*_r.j*_r.k, _d/=_rr;
    *d3fdxdydz-=480.*_d*_r.i*_r.j*_r.k;	
    }
  }
va_end(stack);
return TRUE;
}

//This function calculates values of yff1 electrostatic energy
double inline calc_ffgrid_yff1_Q(t_vec *r,t_ffsys *ffsys,t_top *top)
{
register unsigned int _i;
register double _e=0.;
_i=ffsys->natoms;
while (_i--)
  if (ffsys->a[_i]>0) //Avoid marked active atoms
    _e+=ffsys->q[_i]/sqrt(sqrd(r->i-ffsys->r[_i].i)+sqrd(r->j-ffsys->r[_i].j)+sqrd(r->k-ffsys->r[_i].k)+SMALL2);
return _e;
}
//This function calculates protomap element for yff1 electrostatic
char calc_ffgrid_yff1_protoQ(double *f,double *dfdx,double *dfdy,double *dfdz,double *d2fdxdy,double *d2fdxdz,double *d2fdydz,double *d3fdxdydz,double rr,t_vec *r,char label,va_list stack)
{
unsigned int _i;
double _rr,_d;
t_vec _r;
t_top *top;
t_ffsys *ffsys;
va_list _stack;

va_copy(_stack,stack);
ffsys=va_arg(_stack,t_ffsys*);
top=va_arg(_stack,t_top*);
*f=*dfdx=*dfdy=*dfdz=*d2fdxdy=*d2fdxdz=*d2fdydz=*d3fdxdydz=0.;	
_i=ffsys->natoms;
while (_i--)
  {
  if (ffsys->a[_i]>0) //Avoid marked active atoms
    {
    _r.i=r->i-ffsys->r[_i].i, _r.j=r->j-ffsys->r[_i].j, _r.k=r->k-ffsys->r[_i].k;
    if ((_rr=_r.i*_r.i+_r.j*_r.j+_r.k*_r.k)<rr) _rr=0.5*(rr+_rr);
    _rr+=SMALL2;
    _d=COULOMB_K*ffsys->q[_i]/sqrt(_rr);
    *f+=_d, _d/=_rr,
    *dfdx-=_d*_r.i,
    *dfdy-=_d*_r.j,
    *dfdz-=_d*_r.k, _d/=_rr,
    *d2fdxdy+=3.*_d*_r.i*_r.j,
    *d2fdxdz+=3.*_d*_r.i*_r.k,
    *d2fdydz+=3.*_d*_r.j*_r.k, _d/=_rr;
    *d3fdxdydz-=15.*_d*_r.i*_r.j*_r.k;	
    }
  }
va_end(stack);
return TRUE;
}

//This function calculates energy and gradient for atom on tricubic grid and return TRUE otherwise it returns FALSE 
char calc_grad_yff1_on_threecubic_tcgrid(unsigned int a_i,double q_i,t_vec *r_i,t_vec *g_i,double *e,t_vec *ori,t_len *len,double sp,double (***A)[64],double (***B)[64],double (***Q)[64],t_top *top)
{
double a, b, q;
t_vec da, db, dq;
//Stage I. Calculate position on grid
if (!(calc_threecubic_interpolation_derivative(&a,&da,ori,sp,len->i,len->j,len->k,r_i,A))) return FALSE;
if (!(calc_threecubic_interpolation_derivative(&b,&db,ori,sp,len->i,len->j,len->k,r_i,B))) return FALSE;
if (!(calc_threecubic_interpolation_derivative(&q,&dq,ori,sp,len->i,len->j,len->k,r_i,Q))) return FALSE;
//Stage II. Calculate energies from interpolated components *f+=A*top->A[a_i][0]-B*top->B[a_i][0]+q_i*Q;
*e+=a*top->A[a_i][0]-b*top->B[a_i][0]+q*q_i; 
g_i->i+=da.i*top->A[a_i][0]-db.i*top->B[a_i][0]+dq.i*q_i, g_i->j+=da.j*top->A[a_i][0]-db.j*top->B[a_i][0]+dq.j*q_i, g_i->k+=da.k*top->A[a_i][0]-db.k*top->B[a_i][0]+dq.k*q_i;
return TRUE;
}

//This function calculates energy and gradient for atom on dgrid and return TRUE otherwise it returns FALSE 
char calc_energy_yff1_on_tricubic_dgrid(unsigned int a_i,double q_i,t_vec *r_i,double *e,t_vec *ori,t_len *len,double sp,double ***A,double ***B,double ***Q,t_top *top)
{
double a, b, q;
//Stage I. Calculate position on grid
if (!(calc_tricubic_interpolation_wp_monotonicity(&a,ori,sp,len->i,len->j,len->k,r_i,A))) return FALSE;
if (!(calc_tricubic_interpolation_wp_monotonicity(&b,ori,sp,len->i,len->j,len->k,r_i,B))) return FALSE;
if (!(calc_tricubic_interpolation_wp_monotonicity(&q,ori,sp,len->i,len->j,len->k,r_i,Q))) return FALSE;
//Stage II. Calculate energies from interpolated components
*e=a*top->A[a_i][0]-b*top->B[a_i][0]+q*q_i;
return TRUE;
}
//This function calculates energy and gradient for atom on dgrid and return TRUE otherwise it returns FALSE 
char calc_grad_yff1_on_tricubic_dgrid(unsigned int a_i,double q_i,t_vec *r_i,t_vec *g_i,double *e,t_vec *ori,t_len *len,double sp,double ***A,double ***B,double ***Q,t_top *top)
{
double a, b, q;
t_vec da, db, dq;
//Stage I. Calculate position on grid
if (!(calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,len->i,len->j,len->k,r_i,A))) return FALSE;
if (!(calc_tricubic_interpolation_derivative_wp_monotonicity(&b,&db,ori,sp,len->i,len->j,len->k,r_i,B))) return FALSE;
if (!(calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,len->i,len->j,len->k,r_i,Q))) return FALSE;
//Stage II. Calculate energies from interpolated components
*e+=a*top->A[a_i][0]-b*top->B[a_i][0]+q*q_i;
g_i->i+=da.i*top->A[a_i][0]-db.i*top->B[a_i][0]+dq.i*q_i, g_i->j+=da.j*top->A[a_i][0]-db.j*top->B[a_i][0]+dq.j*q_i, g_i->k+=da.k*top->A[a_i][0]-db.k*top->B[a_i][0]+dq.k*q_i;
return TRUE;
}

//This function required for debugging of threecubic interpolation
char calc_tricubic_tcgrid_validate(double *f,double *dfdx,double *dfdy,double *dfdz,double *d2fdxdy,double *d2fdxdz,double *d2fdydz,double *d3fdxdydz,t_vec *r,char label,va_list stack)
{
double c0=.5, c1=12., c2=3., c3=-75., c4=-21., c5=44., c6=15.5, c7=22.3, c8=-6.6667, c9=17.;
double x=r->i, y=r->j, z=r->k;
*f=c0*x*x*x+c1*x*x*y+c2*x*x*z+c3*x*y*y+c4*y*y*z+c5*y*y*y+c6*x*z*z+c7*y*z*z+c8*z*z*z+c9*x*y*z;
*dfdx=3.*c0*x*x+2.*c1*x*y+2.*c2*x*z+c3*y*y+c6*z*z+c9*y*z;
*dfdy=c1*x*x+2.*c3*x*y+2.*c4*y*z+3.*c5*y*y+c7*z*z+c9*x*z;
*dfdz=c2*x*x+c4*y*y+2.*c6*x*z+2.*c7*y*z+3.*c8*z*z+c9*x*y;
*d2fdxdy=2.*c1*x+2.*c3*y+c9*z;
*d2fdxdz=2.*c2*x+2.*c6*z+c9*y;
*d2fdydz=2.*c4*y+2.*c7*z+c9*x;
*d3fdxdydz=c9;
return TRUE;
}

//-----------------------------   B O N D E D   ---------------------------------------------------   

//This function calculates energy and forces of bond
inline double _calc_gbond_yff1(double r,double *de,t_ff_b *ff_b)
{
r-=ff_b->v;
*de=2.*ff_b->k*r;
return ff_b->k*r*r;
}
//This function calculates energies and forces of angle
inline double _calc_gangle_yff1(double csA,double snA,double *de,t_ff_g *ff_g)
{
double alpha;
alpha=asin(snA*cos(ff_g->v*PI/180.)-csA*sin(ff_g->v*PI/180.));
*de=2.*180.*180./PI/PI*ff_g->k*alpha;
return 180.*180./PI/PI*ff_g->k*alpha*alpha;
}
//This function calculates energies and forces of improper angle
inline double _calc_gimpr_yff1(double csF,double snF,double *de,t_ff_i *ff_i)
{
double phi;
if (fabs(csF)<SMALL2)
  {
  if (snF>0) phi=+PI/2.;
  else       phi=-PI/2.;
  }
else phi=atan(snF/csF);
phi=phi*180./PI-ff_i->v;
*de=2.*180./PI*ff_i->k*phi;
return ff_i->k*phi*phi;
}
//This function calculates energies and forces of torsion
inline double _calc_gtors_yff1(double csF,double snF,double *de,t_ff_t *ff_t)
{
*de=-(ff_t->k[0]-3.*ff_t->k[2]+4.*(ff_t->k[1]+3.*ff_t->k[2]*csF)*csF)*snF;
return (ff_t->k[0]-3.*ff_t->k[2]+2.*(ff_t->k[1]+2.*ff_t->k[2]*csF)*csF)*csF-ff_t->k[1];
}
inline double _calc_gtors_enrg_yff1(double csF,double snF,t_ff_t *ff_t)
{
return (ff_t->k[0]-3.*ff_t->k[2]+2.*(ff_t->k[1]+2.*ff_t->k[2]*csF)*csF)*csF-ff_t->k[1];
}

//This function calculated bonded interactions for whole molecule
inline double calc__bgrad_yff1(t_vec *g,t_vec *r,t_mol *mol,t_activem *activem)
{
unsigned int _i, _j;
double e=0., de, csA, snA, _g[12];

//Do summation
if (!(activem->active_a))
  {
  _i=mol->size_b;
  while (_i--)
    {
    calc_bond_derivative(&csA,&r[mol->ff_b[_i].atom[0]],&r[mol->ff_b[_i].atom[1]],_g);
    e+=_calc_gbond_yff1(csA,&de,&mol->ff_b[_i]);
    g[mol->ff_b[_i].atom[0]].i+=de*_g[0], g[mol->ff_b[_i].atom[0]].j+=de*_g[1], g[mol->ff_b[_i].atom[0]].k+=de*_g[2];
    g[mol->ff_b[_i].atom[1]].i+=de*_g[3], g[mol->ff_b[_i].atom[1]].j+=de*_g[4], g[mol->ff_b[_i].atom[1]].k+=de*_g[5];
    }
  _i=mol->size_g;
  while (_i--)
    {
    calc_angle_derivative(&csA,&snA,&r[mol->ff_g[_i].atom[0]],&r[mol->ff_g[_i].atom[1]],&r[mol->ff_g[_i].atom[2]],_g);
    e+=_calc_gangle_yff1(csA,snA,&de,&mol->ff_g[_i]);
    g[mol->ff_g[_i].atom[0]].i+=de*_g[0], g[mol->ff_g[_i].atom[0]].j+=de*_g[1], g[mol->ff_g[_i].atom[0]].k+=de*_g[2];
    g[mol->ff_g[_i].atom[1]].i+=de*_g[3], g[mol->ff_g[_i].atom[1]].j+=de*_g[4], g[mol->ff_g[_i].atom[1]].k+=de*_g[5];
    g[mol->ff_g[_i].atom[2]].i+=de*_g[6], g[mol->ff_g[_i].atom[2]].j+=de*_g[7], g[mol->ff_g[_i].atom[2]].k+=de*_g[8];
    }
  _i=mol->size_i;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_i[_i].atom[0]],&r[mol->ff_i[_i].atom[1]],&r[mol->ff_i[_i].atom[2]],&r[mol->ff_i[_i].atom[3]],_g);
    e+=_calc_gimpr_yff1(csA,snA,&de,&mol->ff_i[_i]);
    g[mol->ff_i[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_i[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_i[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_i[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_i[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_i[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_i[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_i[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_i[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_i[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_i[_i].atom[3]].j+=de*_g[10], g[mol->ff_i[_i].atom[3]].k+=de*_g[11];
    }
  _i=mol->size_t;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_t[_i].atom[0]],&r[mol->ff_t[_i].atom[1]],&r[mol->ff_t[_i].atom[2]],&r[mol->ff_t[_i].atom[3]],_g);
    e+=_calc_gtors_yff1(csA,snA,&de,&mol->ff_t[_i]);
    g[mol->ff_t[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_t[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_t[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_t[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_t[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_t[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_t[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_t[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_t[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_t[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_t[_i].atom[3]].j+=de*_g[10], g[mol->ff_t[_i].atom[3]].k+=de*_g[11];
    }
  }
else
  {
  _j=activem->active_b->size;
  while (_j--)
    {
    _i=activem->active_b->list[_j];
    calc_bond_derivative(&csA,&r[mol->ff_b[_i].atom[0]],&r[mol->ff_b[_i].atom[1]],_g);
    e+=_calc_gbond_yff1(csA,&de,&mol->ff_b[_i]);
    g[mol->ff_b[_i].atom[0]].i+=de*_g[0], g[mol->ff_b[_i].atom[0]].j+=de*_g[1], g[mol->ff_b[_i].atom[0]].k+=de*_g[2];
    g[mol->ff_b[_i].atom[1]].i+=de*_g[3], g[mol->ff_b[_i].atom[1]].j+=de*_g[4], g[mol->ff_b[_i].atom[1]].k+=de*_g[5];
    }
  _j=activem->active_g->size;
  while (_j--)
    {
    _i=activem->active_g->list[_j];
    calc_angle_derivative(&csA,&snA,&r[mol->ff_g[_i].atom[0]],&r[mol->ff_g[_i].atom[1]],&r[mol->ff_g[_i].atom[2]],_g);
    e+=_calc_gangle_yff1(csA,snA,&de,&mol->ff_g[_i]);
    g[mol->ff_g[_i].atom[0]].i+=de*_g[0], g[mol->ff_g[_i].atom[0]].j+=de*_g[1], g[mol->ff_g[_i].atom[0]].k+=de*_g[2];
    g[mol->ff_g[_i].atom[1]].i+=de*_g[3], g[mol->ff_g[_i].atom[1]].j+=de*_g[4], g[mol->ff_g[_i].atom[1]].k+=de*_g[5];
    g[mol->ff_g[_i].atom[2]].i+=de*_g[6], g[mol->ff_g[_i].atom[2]].j+=de*_g[7], g[mol->ff_g[_i].atom[2]].k+=de*_g[8];
    }
  _j=activem->active_i->size;
  while (_j--)
    {
    _i=activem->active_i->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_i[_i].atom[0]],&r[mol->ff_i[_i].atom[1]],&r[mol->ff_i[_i].atom[2]],&r[mol->ff_i[_i].atom[3]],_g);
    e+=_calc_gimpr_yff1(csA,snA,&de,&mol->ff_i[_i]);
    g[mol->ff_i[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_i[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_i[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_i[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_i[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_i[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_i[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_i[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_i[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_i[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_i[_i].atom[3]].j+=de*_g[10], g[mol->ff_i[_i].atom[3]].k+=de*_g[11];
    }
  _j=activem->active_t->size;
  while (_j--)
    {
    _i=activem->active_t->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_t[_i].atom[0]],&r[mol->ff_t[_i].atom[1]],&r[mol->ff_t[_i].atom[2]],&r[mol->ff_t[_i].atom[3]],_g);
    e+=_calc_gtors_yff1(csA,snA,&de,&mol->ff_t[_i]);
    g[mol->ff_t[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_t[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_t[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_t[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_t[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_t[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_t[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_t[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_t[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_t[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_t[_i].atom[3]].j+=de*_g[10], g[mol->ff_t[_i].atom[3]].k+=de*_g[11];
    }
  }
return e;
}
inline void _calc__bgrad_yff1(double e[N_ESLICES],t_vec (*g)[N_ESLICES],t_vec *r,t_mol *mol,t_activem *activem)
{
unsigned int _i, _j, _L;
double _e, de, csA, snA, _g[12];

//Do summation
if (!(activem->active_a))
  {
  _i=mol->size_b;
  while (_i--)
    {
    calc_bond_derivative(&csA,&r[mol->ff_b[_i].atom[0]],&r[mol->ff_b[_i].atom[1]],_g);
    _e=_calc_gbond_yff1(csA,&de,&mol->ff_b[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_b[_i].atom[0]][_L].i+=de*_g[0], g[mol->ff_b[_i].atom[0]][_L].j+=de*_g[1], g[mol->ff_b[_i].atom[0]][_L].k+=de*_g[2],
    g[mol->ff_b[_i].atom[1]][_L].i+=de*_g[3], g[mol->ff_b[_i].atom[1]][_L].j+=de*_g[4], g[mol->ff_b[_i].atom[1]][_L].k+=de*_g[5];
    }
  _i=mol->size_g;
  while (_i--)
    {
    calc_angle_derivative(&csA,&snA,&r[mol->ff_g[_i].atom[0]],&r[mol->ff_g[_i].atom[1]],&r[mol->ff_g[_i].atom[2]],_g);
    _e=_calc_gangle_yff1(csA,snA,&de,&mol->ff_g[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_g[_i].atom[0]][_L].i+=de*_g[0], g[mol->ff_g[_i].atom[0]][_L].j+=de*_g[1], g[mol->ff_g[_i].atom[0]][_L].k+=de*_g[2],
    g[mol->ff_g[_i].atom[1]][_L].i+=de*_g[3], g[mol->ff_g[_i].atom[1]][_L].j+=de*_g[4], g[mol->ff_g[_i].atom[1]][_L].k+=de*_g[5],
    g[mol->ff_g[_i].atom[2]][_L].i+=de*_g[6], g[mol->ff_g[_i].atom[2]][_L].j+=de*_g[7], g[mol->ff_g[_i].atom[2]][_L].k+=de*_g[8];
    }
  _i=mol->size_i;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_i[_i].atom[0]],&r[mol->ff_i[_i].atom[1]],&r[mol->ff_i[_i].atom[2]],&r[mol->ff_i[_i].atom[3]],_g);
    _e=_calc_gimpr_yff1(csA,snA,&de,&mol->ff_i[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_i[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_i[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_i[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_i[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_i[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_i[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_i[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_i[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_i[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_i[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_i[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_i[_i].atom[3]][_L].k+=de*_g[11];
    }
  _i=mol->size_t;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_t[_i].atom[0]],&r[mol->ff_t[_i].atom[1]],&r[mol->ff_t[_i].atom[2]],&r[mol->ff_t[_i].atom[3]],_g);
    _e=_calc_gtors_yff1(csA,snA,&de,&mol->ff_t[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_t[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_t[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_t[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_t[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_t[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_t[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_t[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_t[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_t[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_t[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_t[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_t[_i].atom[3]][_L].k+=de*_g[11];
    }
  }
else
  {
  _j=activem->active_b->size;
  while (_j--)
    {
    _i=activem->active_b->list[_j];
    calc_bond_derivative(&csA,&r[mol->ff_b[_i].atom[0]],&r[mol->ff_b[_i].atom[1]],_g);
    _e=_calc_gbond_yff1(csA,&de,&mol->ff_b[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_b[_i].atom[0]][_L].i+=de*_g[0], g[mol->ff_b[_i].atom[0]][_L].j+=de*_g[1], g[mol->ff_b[_i].atom[0]][_L].k+=de*_g[2],
    g[mol->ff_b[_i].atom[1]][_L].i+=de*_g[3], g[mol->ff_b[_i].atom[1]][_L].j+=de*_g[4], g[mol->ff_b[_i].atom[1]][_L].k+=de*_g[5];
    }
  _j=activem->active_g->size;
  while (_j--)
    {
    _i=activem->active_g->list[_j];
    calc_angle_derivative(&csA,&snA,&r[mol->ff_g[_i].atom[0]],&r[mol->ff_g[_i].atom[1]],&r[mol->ff_g[_i].atom[2]],_g);
    _e=_calc_gangle_yff1(csA,snA,&de,&mol->ff_g[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_g[_i].atom[0]][_L].i+=de*_g[0], g[mol->ff_g[_i].atom[0]][_L].j+=de*_g[1], g[mol->ff_g[_i].atom[0]][_L].k+=de*_g[2],
    g[mol->ff_g[_i].atom[1]][_L].i+=de*_g[3], g[mol->ff_g[_i].atom[1]][_L].j+=de*_g[4], g[mol->ff_g[_i].atom[1]][_L].k+=de*_g[5],
    g[mol->ff_g[_i].atom[2]][_L].i+=de*_g[6], g[mol->ff_g[_i].atom[2]][_L].j+=de*_g[7], g[mol->ff_g[_i].atom[2]][_L].k+=de*_g[8];
    }
  _j=activem->active_i->size;
  while (_j--)
    {
    _i=activem->active_i->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_i[_i].atom[0]],&r[mol->ff_i[_i].atom[1]],&r[mol->ff_i[_i].atom[2]],&r[mol->ff_i[_i].atom[3]],_g);
    _e=_calc_gimpr_yff1(csA,snA,&de,&mol->ff_i[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_i[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_i[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_i[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_i[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_i[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_i[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_i[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_i[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_i[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_i[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_i[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_i[_i].atom[3]][_L].k+=de*_g[11];
    }
  _j=activem->active_t->size;
  while (_j--)
    {
    _i=activem->active_t->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_t[_i].atom[0]],&r[mol->ff_t[_i].atom[1]],&r[mol->ff_t[_i].atom[2]],&r[mol->ff_t[_i].atom[3]],_g);
    _e=_calc_gtors_yff1(csA,snA,&de,&mol->ff_t[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_t[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_t[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_t[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_t[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_t[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_t[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_t[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_t[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_t[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_t[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_t[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_t[_i].atom[3]][_L].k+=de*_g[11];
    }
  }
}

//This function calculates torsional energy for internal coordinates gradients
inline double calc_mol_tbgrad_yff1(t_vec *g,t_vec *r,t_mol *mol,t_activem *activem)
{
unsigned int _i,_j;
double e=0., de, csA, snA, _g[12];

if (!(activem->active_a))
  {
  _i=mol->size_t;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_t[_i].atom[0]],&r[mol->ff_t[_i].atom[1]],&r[mol->ff_t[_i].atom[2]],&r[mol->ff_t[_i].atom[3]],_g);
    e+=_calc_gtors_yff1(csA,snA,&de,&mol->ff_t[_i]);
    g[mol->ff_t[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_t[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_t[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_t[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_t[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_t[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_t[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_t[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_t[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_t[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_t[_i].atom[3]].j+=de*_g[10], g[mol->ff_t[_i].atom[3]].k+=de*_g[11];
    }
  }
else
  {
  _j=activem->active_t->size;
  while (_j--)
    {
    _i=activem->active_t->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_t[_i].atom[0]],&r[mol->ff_t[_i].atom[1]],&r[mol->ff_t[_i].atom[2]],&r[mol->ff_t[_i].atom[3]],_g);
    e+=_calc_gtors_yff1(csA,snA,&de,&mol->ff_t[_i]);
    g[mol->ff_t[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_t[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_t[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_t[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_t[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_t[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_t[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_t[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_t[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_t[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_t[_i].atom[3]].j+=de*_g[10], g[mol->ff_t[_i].atom[3]].k+=de*_g[11];
    }
  }
return e;
}
inline void _calc_mol_tbgrad_yff1(double e[N_ESLICES],t_vec (*g)[N_ESLICES],t_vec *r,t_mol *mol,t_activem *activem)
{
unsigned int _i, _j, _L;
double _e, de, csA, snA, _g[12];

if (!(activem->active_a))
  {
  _i=mol->size_t;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_t[_i].atom[0]],&r[mol->ff_t[_i].atom[1]],&r[mol->ff_t[_i].atom[2]],&r[mol->ff_t[_i].atom[3]],_g);
    _e=_calc_gtors_yff1(csA,snA,&de,&mol->ff_t[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_t[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_t[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_t[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_t[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_t[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_t[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_t[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_t[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_t[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_t[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_t[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_t[_i].atom[3]][_L].k+=de*_g[11];
    }
  }
else
  {
  _j=activem->active_t->size;
  while (_j--)
    {
    _i=activem->active_t->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_t[_i].atom[0]],&r[mol->ff_t[_i].atom[1]],&r[mol->ff_t[_i].atom[2]],&r[mol->ff_t[_i].atom[3]],_g);
    _e=_calc_gtors_yff1(csA,snA,&de,&mol->ff_t[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_t[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_t[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_t[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_t[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_t[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_t[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_t[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_t[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_t[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_t[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_t[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_t[_i].atom[3]][_L].k+=de*_g[11];
    }
  }
}

//------------------------------------    N O N B O N D E D   ------------------------------------------

//This function calculates energy and force of nonbonded interatons in YFF1. The flag defines if energy(forces) adds of subtracts
//Note g_vecs and r_vecs migh be ffsys->g and ffsys->r correspondingly
inline double calc_atom__nb_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top)
{
double _r, _x, _rr, _rrrrrr, _a, _b, _q, _d, A, B, C, D, E;
t_vec _u;
_u.i=r_i->i-r_j->i, _u.j=r_i->j-r_j->j, _u.k=r_i->k-r_j->k, _rr=_u.i*_u.i+_u.j*_u.j+_u.k*_u.k+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _a=top->A[a_i][a_j], _b=top->B[a_i][a_j], _q=COULOMB_K*q_i*q_j, _rrrrrr=_rr*_rr*_rr, _r=sqrt(_rr);
  if (_r<YFF1_Rb)
    {//Do classic interaction
    _d=(-12.*_a/_rrrrrr+6.*_b)/_rrrrrr/_r-_q/_rr;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return (_a/_rrrrrr-_b)/_rrrrrr+_q/_r;
    }
  else 
    {//Shift interaction function
    _x=YFF1_Rx(_r), E=YFF1_pE(_a,_b,_q), D=YFF1_pD(_a,_b,_q), C=YFF1_pC(_a,_b,_q), B=YFF1_pB(C,D,E), A=YFF1_pA(C,D,E), _d=YFF1_pF(_x,A,B,C,D);
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return YFF1_pV(_x,A,B,C,D,E);
    }
  }
else return 0.;
}
inline double calc_atom__nb_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top)
{
double _r, _rr, _rrrrrr, _a, _b, _q, C, D, E;
_rr=sqrd(r_i->i-r_j->i)+sqrd(r_i->j-r_j->j)+sqrd(r_i->k-r_j->k)+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _a=top->A[a_i][a_j], _b=top->B[a_i][a_j], _q=COULOMB_K*q_i*q_j, _rrrrrr=_rr*_rr*_rr, _r=sqrt(_rr);
  if (_r<YFF1_Rb) return (_a/_rrrrrr-_b)/_rrrrrr+_q/_r;                                                                                          //Do classic interaction
  else { E=YFF1_pE(_a,_b,_q), D=YFF1_pD(_a,_b,_q), C=YFF1_pC(_a,_b,_q); return YFF1_pV(YFF1_Rx(_r),YFF1_pA(C,D,E),YFF1_pB(C,D,E),C,D,E); } //Do shift interaction
  }
else return 0.;
}
inline double calc_atom__znb_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top)
{
double _r, _x, _rr, _rrrrrr, _a, _b, _q, _d, A, B, C, D, E;
t_vec _u;
_u.i=r_i->i-r_j->i, _u.j=r_i->j-r_j->j, _u.k=r_i->k-r_j->k, _rr=_u.i*_u.i+_u.j*_u.j+_u.k*_u.k+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _a=top->A[a_i][0]*top->A[a_j][0], _b=top->B[a_i][0]*top->B[a_j][0], _q=COULOMB_K*q_i*q_j, _rrrrrr=_rr*_rr*_rr, _r=sqrt(_rr);
  if (_r<YFF1_Rb)
    {//Do classic interaction
    _d=(-12.*_a/_rrrrrr+6.*_b)/_rrrrrr/_r-_q/_rr;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return (_a/_rrrrrr-_b)/_rrrrrr+_q/_r;
    }
  else 
    {//Shift interaction function
    _x=YFF1_Rx(_r), E=YFF1_pE(_a,_b,_q), D=YFF1_pD(_a,_b,_q), C=YFF1_pC(_a,_b,_q), B=YFF1_pB(C,D,E), A=YFF1_pA(C,D,E), _d=YFF1_pF(_x,A,B,C,D);
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return YFF1_pV(_x,A,B,C,D,E);
    }
  }
else return 022.;
}
inline double calc_atom__znb_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top)
{
double _r, _rr, _rrrrrr, _a, _b, _q, C, D, E;
_rr=sqrd(r_i->i-r_j->i)+sqrd(r_i->j-r_j->j)+sqrd(r_i->k-r_j->k)+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _a=top->A[a_i][0]*top->A[a_j][0], _b=top->B[a_i][0]*top->B[a_j][0], _q=COULOMB_K*q_i*q_j, _rrrrrr=_rr*_rr*_rr, _r=sqrt(_rr);
  if (_r<YFF1_Rb) return (_a/_rrrrrr-_b)/_rrrrrr+_q/_r;                                                                                          //Do classic interaction
  else { E=YFF1_pE(_a,_b,_q), D=YFF1_pD(_a,_b,_q), C=YFF1_pC(_a,_b,_q); return YFF1_pV(YFF1_Rx(_r),YFF1_pA(C,D,E),YFF1_pB(C,D,E),C,D,E); } //Do shift interaction
  }
else return 0.;
}
//This function calculates 1-4 coulomb interactions only
inline double calc_atom__14_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top)
{
double _x, _r, _rr, _q, _d, A, B, C, D, E;
t_vec _u;
_u.i=r_i->i-r_j->i, _u.j=r_i->j-r_j->j, _u.k=r_i->k-r_j->k, _rr=_u.i*_u.i+_u.j*_u.j+_u.k*_u.k+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _q=COULOMB_K*q_i*q_j, _r=sqrt(_rr);
  if (_r<YFF1_Rb)
    {//Do classic interaction
    _d=-YFF1_14_SCALE*_q/_rr;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return YFF1_14_SCALE*_q/_r;
    }
  else 
    {//Shift interaction function
    _x=YFF1_Rx(_r), E=YFF1_14_SCALE*_q/YFF1_Rb, D=-YFF1_14_SCALE*_q/YFF1_Rb2, C=YFF1_14_SCALE*_q/YFF1_Rb2/YFF1_Rb, B=YFF1_pB(C,D,E), A=YFF1_pA(C,D,E), _d=YFF1_pF(_x,A,B,C,D);
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return YFF1_pV(_x,A,B,C,D,E);
    }
  }
else return 0.;
}
inline double calc_atom__14_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top)
{
double _r, _rr, _q, C, D, E;
_rr=sqrd(r_i->i-r_j->i)+sqrd(r_i->j-r_j->j)+sqrd(r_i->k-r_j->k)+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _q=COULOMB_K*q_i*q_j, _r=sqrt(_rr);
  if (_r<YFF1_Rb) return YFF1_14_SCALE*_q/_rr;                           //Do classic interaction
  else { E=YFF1_14_SCALE*_q/YFF1_Rb, D=-YFF1_14_SCALE*_q/YFF1_Rb2, C=YFF1_14_SCALE*_q/YFF1_Rb2/YFF1_Rb; return YFF1_pV(YFF1_Rx(_r),YFF1_pA(C,D,E),YFF1_pB(C,D,E),C,D,E); } //Shift interaction function
  }
else return 0.;
}
//This function calculates negative difference between 1-4 and normal coulomb interactions only
inline double calc_atom_n14_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top)
{
double _x, _r, _rr, _q, _d, A, B, C, D, E;
t_vec _u;
_u.i=r_i->i-r_j->i, _u.j=r_i->j-r_j->j, _u.k=r_i->k-r_j->k, _rr=_u.i*_u.i+_u.j*_u.j+_u.k*_u.k+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _q=COULOMB_K*q_i*q_j, _r=sqrt(_rr);
  if (_r<YFF1_Rb)
    {//Do classic interaction
    _d=-(1.*YFF1_14_SCALE)*_q/_rr;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return (1.-YFF1_14_SCALE)*_q/_r;
    }
  else 
    {//Shift interaction function
    _x=YFF1_Rx(_r), E=(1.-YFF1_14_SCALE)*_q/YFF1_Rb, D=-(1.-YFF1_14_SCALE)*_q/YFF1_Rb2, C=(1.-YFF1_14_SCALE)*_q/YFF1_Rb2/YFF1_Rb, B=YFF1_pB(C,D,E), A=YFF1_pA(C,D,E), _d=YFF1_pF(_x,A,B,C,D);
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return YFF1_pV(_x,A,B,C,D,E);
    }
  }
else return 0.;
}
inline double calc_atom_n14_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top)
{
double _r, _rr, _q, C, D, E;
_rr=sqrd(r_i->i-r_j->i)+sqrd(r_i->j-r_j->j)+sqrd(r_i->k-r_j->k)+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _q=COULOMB_K*q_i*q_j, _r=sqrt(_rr);
  if (_r<YFF1_Rb) return (1.-YFF1_14_SCALE)*_q/_rr;                           //Do classic interaction
  else { E=(1.-YFF1_14_SCALE)*_q/YFF1_Rb, D=-(1.-YFF1_14_SCALE)*_q/YFF1_Rb2, C=(1.-YFF1_14_SCALE)*_q/YFF1_Rb2/YFF1_Rb; return YFF1_pV(YFF1_Rx(_r),YFF1_pA(C,D,E),YFF1_pB(C,D,E),C,D,E); } //Shift interaction function
  }
else return 0.;
}
//This function calculates potential on grid
inline double calc_atom__znb_yff1_on_grid(double a,double b,double q,t_vec *da,t_vec *db,t_vec *dq,unsigned int _a_i,double q_i,t_vec *g_i,t_top *top)
{
double a_i=top->A[_a_i][0], b_i=top->B[_a_i][0];
g_i->i=a_i*da->i-b_i*db->i+q_i*dq->i, g_i->j=a_i*da->j-b_i*db->j+q_i*dq->j, g_i->k=a_i*da->k-b_i*db->k+q_i*dq->k;
return a_i*a-b_i*b+q_i*q; 
}
inline double calc_atom_enrg__znb_yff1_on_grid(double a,double b,double q,unsigned int a_i,double q_i,t_top *top)
{
return top->A[a_i][0]*a-top->B[a_i][0]*b+q_i*q; 
}

//Note it scales repulsion as E*=E*(3*x^2-2*x^3) in range 0 @ px==0. ... E @ px=1. dE*/dr=(3*x^2-2*x^3)*dE/dr
//This function calculates energy in the active layer and marks atoms to exclude from grid summing
inline void calc_yff1_grad(double e[N_ESLICES],t_ffsys *ffsys,t_top *top)
{
unsigned int _i, _j, _k, _l, mol_i, mol_j, active_i, active_j, a_i, a_j, _L;
t_vec g_i, g_j;
double _e;

mol_i=ffsys->nmols;
while (mol_i--)
  if (!(ffsys->activem[mol_i].active_a))
    {
    active_i=ffsys->mols[mol_i]->anchors->size;
    while (active_i--)
      {//Summ inside molecule
      //Do all but bonded
      _k=active_i;
      LOOP_A: _l=0, _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch; while (_i--) { _j=ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]; if ((_j<_k)&&(_j>=_l)) _l=_j+1; }
      if (!(_l))
        {//All neighbors are scanned, summ till zero
        active_j=_k;
        while (active_j--)
          {
          _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
          while (_i--)
            {
            a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]; 
            _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
            while (_j--)
              {
              a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
              if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                } 
              }
            }
          }
        }
      else
        {//Scann till neighbour
        _l--, active_j=_k;
        while (--active_j!=_l)
          {
          _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
          while (_i--)
            {
            a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]; 
            _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
            while (_j--)
              {
              a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
              if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                }
              }
            }
          }
        _k=_l;
        goto LOOP_A; 
        }
      //Do bonded anchors
      if (active_i!=ffsys->rtrees[mol_i]->root)
        {
        active_j=*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch;
        _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (--_i)
          {
          a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]; 
          _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
          while (_j--)
            {
            a_j=ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
            if (a_j!=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0])
              {
              a_j+=ffsys->nr[mol_i]; 
              if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                }
              }
            }
          }
        }
      }
    //Edit 1-4
    active_i=ffsys->mols[mol_i]->size_t;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_t[active_i].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_t[active_i].atom[3]; 
      if ( (_e=calc_atom_n14_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
        {
        if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
        else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
        e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
        }
      }
    //Remove pairs
    active_i=ffsys->mols[mol_i]->size_p;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_p[active_i].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_p[active_i].atom[1]; 
      if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
        {
        if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
        else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
        e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
        }
      }
    //Summ cross molecules
    mol_j=mol_i; 
    while (mol_j--)
      {
      if (!(ffsys->activem[mol_j].active_a))
        {
        _i=ffsys->mols[mol_i]->atoms->size;
        while (_i--)
          {
          a_i=ffsys->nr[mol_i]+_i;
          _j=ffsys->mols[mol_j]->atoms->size;
          while (_j--)
            {
            a_j=ffsys->nr[mol_j]+_j;
            if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
              {
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
              }
            }
          }
        }
      else
        {
        _i=ffsys->mols[mol_i]->atoms->size;
        while (_i--)
          {
          a_i=ffsys->nr[mol_i]+_i;
          _l=ffsys->activem[mol_j].active_a->size;
          while (_l--)
            {
            active_j=ffsys->activem[mol_j].active_a->list[_l];
            _j=ffsys->mols[mol_j]->anchors->list[active_j].size; 
            while (_j--)
              {
              a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j];
              if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                }
              }
            }
          }
        }
      }
    }
  else
    {//Summ active mol to all
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {//Do inside mol
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      _l=_k;
      while (_l--)
        {
        active_j=ffsys->activem[mol_i].active_a->list[_l];
             if (*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch==active_j)
               {
               _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
               while (--_i)
                 {
                 a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
                 _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
                 while (_j--)
                   {
                   a_j=ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
                   if (a_j!=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0])
                     {
                     a_j+=ffsys->nr[mol_i];
                     if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                       {
                       if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                       else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                       e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                       }
                     }
                   }
                 }
               //Mark as processed neighbour
               *ffsys->rtrees[mol_i]->rbranch[active_i].rbranch=-(int)(*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch); 
               }
        else if (*ffsys->rtrees[mol_i]->rbranch[active_j].rbranch==active_i)
               {
               _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
               while (_i--)
                 {
                 a_i=ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
                 if (a_i!=ffsys->rtrees[mol_i]->rbranch[active_j].edge.vertice[0])
                   {
                   a_i+=ffsys->nr[mol_i], _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
                   while (--_j)
                     {
                     a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
                     if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                       {
                       if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                       else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                       e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                       }
                     }
                   }
                 }
               //Mark as processed neighbour
               *ffsys->rtrees[mol_i]->rbranch[active_j].rbranch=-(int)(*ffsys->rtrees[mol_i]->rbranch[active_j].rbranch); 
               }
        else   {//General case: anchors are not connected
               _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
               while (_i--)
                 {
                 a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
                 _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
                 while (_j--)
                   {
                   a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
                   if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                     {
                     if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                     else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                     e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                     } 
                   }
                 }
               }
        }
      //Remove cross-anchors !root
//      _l=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch; 
//      while (--_l)
//        if ((int)*ffsys->rtrees[mol_i]->rbranch[active_j=ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_l]].rbranch>0) //Unmarked root, remove it explicitly
//          {
//          a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_j].edge.vertice[0];
//          _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
//          while (--_j)
//            {
//            a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
//            if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
//              {
//              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
//              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
//              e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
//              }
//            }
//          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_j].edge.vertice[1];
//          _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
//          while (_i--)
//            {
//            a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
//            if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
//              {
//              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
//              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
//              e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
//              }
//            }
//          }
      }
    //Do inside mol: remove cross-anchor == root // This function is partially supressed for accuracy purposes: no rescore + no grid on directly bonded atom
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      if ((int)(*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch)>0)
        {
        active_j=*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch;
        a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0];
        _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (--_i)
          {
          a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
          if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
            }
          }
    //    a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[1];
    //    _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
    //    while (_j--)
    //      {
    //      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
    //      if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
    //        {
    //        if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    //        else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    //        e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
    //        }
    //      }
        }
      }
    //Edit torsions
    active_i=ffsys->activem[mol_i].active_t->size;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_t[ffsys->activem[mol_i].active_t->list[active_i]].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_t[ffsys->activem[mol_i].active_t->list[active_i]].atom[3]; 
      if ( (_e=calc_atom_n14_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
        {
        if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
        else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
        e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
        }
      }
    //Remove pairs
    active_i=ffsys->activem[mol_i].active_p->size;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_p[ffsys->activem[mol_i].active_p->list[active_i]].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_p[ffsys->activem[mol_i].active_p->list[active_i]].atom[1]; 
      if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
        {
        if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
        else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
        e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
        } 
      }
    //Summ cross molecules                    !!!!! Not debugged yet !!!!!!!!
    mol_j=mol_i;
    while (mol_j--)
      {
      if (!(ffsys->activem[mol_j].active_a))
        {
        _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (_i--)
          {
          a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
          _j=ffsys->mols[mol_j]->atoms->size;
          while (_j--)
            {
            a_j=ffsys->nr[mol_j]+_j; 
            if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
              {
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
              }
            }
          }
        }
      else
        {
        _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (_i--)
          {
          a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
          _l=ffsys->activem[mol_j].active_a->size;
          while (_l--)
            {
            active_j=ffsys->activem[mol_j].active_a->list[_l];
            _j=ffsys->mols[mol_j]->anchors->list[active_j].size;
            while (_j--)
              {
              a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j];
              if (!(_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                }   
              } 
            }
          }
        }  
      }
    }
}
//The same as previous but with exclusions
inline void calc_yff1_grad_wexclusions(double e[N_ESLICES],unsigned char *exclude,t_ffsys *ffsys,t_top *top)
{
unsigned int _i, _j, _k, _l, mol_i, mol_j, active_i, active_j, a_i, a_j, _L;
t_vec g_i, g_j;
double _e;

mol_i=ffsys->nmols;
while (mol_i--)
  if (!(ffsys->activem[mol_i].active_a))
    {
    active_i=ffsys->mols[mol_i]->anchors->size;
    while (active_i--)
      {//Summ inside molecule
      //Do all but bonded
      _k=active_i;
      LOOP_A: _l=0, _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch; while (_i--) { _j=ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]; if ((_j<_k)&&(_j>=_l)) _l=_j+1; }
      if (!(_l))
        {//Scann till zero
        active_j=_k;
        while (active_j--)
          {
          _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
          while (_i--)
            {
            if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]])
              { 
              _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
                if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                  } 
                }
              }
            else 
              { //Skip exlcude-exclude
              _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
              while (_j--)
                if (!exclude[a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j]])
                  {
                  if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                    } 
                  }
              }
            }
          }
        }
      else
        {//Scann till neighbour
        _l--, active_j=_k;
        while (--active_j!=_l)
          {
          _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
          while (_i--)
            if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]]) 
              {
              _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
                if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                  }
                }
              }
            else
              {
              _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
              while (_j--)
                if (!exclude[a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j]])
                  { //Skip exlcude-exclude
                  if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                    }
                  }
              }
          }
        _k=_l;
        goto LOOP_A; 
        }
      //Do bonded anchors
      if (active_i!=ffsys->rtrees[mol_i]->root)
        {
        active_j=*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch;
        _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (--_i)
          {
          if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]]) 
            {
            _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
            while (_j--)
              {
              a_j=ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
              if (a_j!=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0])
                {
                a_j+=ffsys->nr[mol_i]; 
                if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                  }
                }
              }
            }
          else
            {
            _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
            while (_j--)
              {
              a_j=ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
              if (a_j!=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0])
                if (!exclude[a_j+=ffsys->nr[mol_i]]) //skip excl-excl 
                  if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                    }
              }
            }
          }
        }
      }
    //Edit 1-4
    active_i=ffsys->mols[mol_i]->size_t;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_t[active_i].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_t[active_i].atom[3];
      if ( (!exclude[a_i])||(!exclude[a_j]) ) //skip excl-excl 
        if ( (_e=calc_atom_n14_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
          {
          if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
          else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
          e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
          }
      }
    //Remove pairs
    active_i=ffsys->mols[mol_i]->size_p;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_p[active_i].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_p[active_i].atom[1]; 
      if ( (!exclude[a_i])||(!exclude[a_j]) ) //skip excl-excl 
        if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
          {
          if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
          else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
          e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
          }
      }
    //Summ cross molecules
    mol_j=mol_i; 
    while (mol_j--)
      {
      if (!(ffsys->activem[mol_j].active_a))
        {
        _i=ffsys->mols[mol_i]->atoms->size;
        while (_i--)
          if (!exclude[a_i=ffsys->nr[mol_i]+_i]) 
            {
            _j=ffsys->mols[mol_j]->atoms->size;
            while (_j--)
              {
              a_j=ffsys->nr[mol_j]+_j;
              if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                }
              }
            }
          else
            {
            _j=ffsys->mols[mol_j]->atoms->size;
            while (_j--)
              if (!exclude[a_j=ffsys->nr[mol_j]+_j]) 
                {//skip excl-excl
                if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                  }
                }
            } 
        }
      else
        {
        _i=ffsys->mols[mol_i]->atoms->size;
        while (_i--)
          if (!exclude[a_i=ffsys->nr[mol_i]+_i]) 
            {
            _l=ffsys->activem[mol_j].active_a->size;
            while (_l--)
              {
              active_j=ffsys->activem[mol_j].active_a->list[_l];
              _j=ffsys->mols[mol_j]->anchors->list[active_j].size; 
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j];
                if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                  }
                }
              } 
            }
          else 
            {
            _l=ffsys->activem[mol_j].active_a->size;
            while (_l--)
              {
              active_j=ffsys->activem[mol_j].active_a->list[_l];
              _j=ffsys->mols[mol_j]->anchors->list[active_j].size; 
              while (_j--)
                if (!exclude[a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j]])
                  { //skip excl-excl
                  if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                    }
                  }
              } 
            }
        }
      }
    }
  else
    {//Summ active mol to all
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {//Do inside mol
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      _l=_k;
      while (_l--)
        {
        active_j=ffsys->activem[mol_i].active_a->list[_l];
             if (*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch==active_j)
               {
               _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
               while (--_i)
                 if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]])
                   {
                   _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
                   while (_j--)
                     {
                     a_j=ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
                     if (a_j!=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0])
                       {
                       a_j+=ffsys->nr[mol_i];
                       if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                         {
                         if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                         else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                         e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                         }
                       }
                     }
                   }
                 else
                   {
                   _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
                   while (--_j)
                     {
                     a_j=ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
                     if (a_j!=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0])
                       if (!exclude[a_j+=ffsys->nr[mol_i]]) // skip excl-excl
                         if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                           {
                           if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                           else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                           e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                           }
                     }
                   }
               //Mark as processed neighbour
               *ffsys->rtrees[mol_i]->rbranch[active_i].rbranch=-(int)(*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch); 
               }
        else if (*ffsys->rtrees[mol_i]->rbranch[active_j].rbranch==active_i)
               {
               _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
               while (_i--)
                 {
                 a_i=ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
                 if (a_i!=ffsys->rtrees[mol_i]->rbranch[active_j].edge.vertice[0])
                   {
                   if (!exclude[a_i+=ffsys->nr[mol_i]])
                     {
                     _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
                     while (--_j)
                       {
                       a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
                       if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                         {
                         if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                         else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                         e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                         }
                       }
                     }
                   else
                     {
                     _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
                     while (--_j)
                       if (!exclude[a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j]]) //skip excl-excl
                         if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                           {
                           if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                           else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                           e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                           }
                     }
                   }
                 }
               //Mark as processed neighbour
               *ffsys->rtrees[mol_i]->rbranch[active_j].rbranch=-(int)(*ffsys->rtrees[mol_i]->rbranch[active_j].rbranch); 
               }
        else   {
               _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
               while (_i--)
                 if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]])
                   {
                   _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
                   while (_j--)
                     {
                     a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
                     if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                       {
                       if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                       else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                       e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                       } 
                     }
                   }
                 else
                   {
                   _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
                   while (_j--)
                     if (!exclude[a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j]]) //skip excl-excl
                       if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                         {
                         if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                         else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                         e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                         } 
                   }
               }
        }
      //Remove cross-anchors !root
//      _l=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch; 
//      while (--_l)
//        if ((int)*ffsys->rtrees[mol_i]->rbranch[active_j=ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_l]].rbranch>0) //Unmarked root, remove it explicitly
//          {
//          if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_j].edge.vertice[0]])
//            {
//            _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
//            while (--_j)
//              {
//              a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
//              if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
//                {
//                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
//                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
//                e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
//                }
//              }
//            }
//          else
//            {
//            _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
//            while (--_j)
//              if (!exclude[a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j]]) //skip exlc-excl
//                if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
//                  {
//                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
//                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
//                  e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
//                  }
//            }
//          if (!exclude[a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_j].edge.vertice[1]])
//            { 
//            _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
//            while (_i--)
//              {
//              a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
//              if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
//                {
//                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
//                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
//                e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
//                }
//              }
//            }
//          else
//            { 
//            _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
//            while (_i--)
//              if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]]) //skip excl-excl
//                if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
//                  {
//                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
//                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
//                  e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
//                  }
//            } 
//          }
      }
    //Do inside mol: remove cross-anchor == root // This function is partially supressed for accuracy purposes: no rescore + no grid on directly bonded atom
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      if ((int)(*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch)>0)
        {
        active_j=*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch;
        if (!exclude[a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0]])
          {
          _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
          while (--_i)
            {
            a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
            if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
              {
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
              }
            }
          }
        else
          {
          _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
          while (--_i)
            if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]]) //skip excl-excl
              if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
                }
          }   
    //    a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[1];
    //    _j=ffsys->mols[mol_i]->anchors->list[active_j].size;
    //    while (_j--)
    //      {
    //      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_j].list[_j];
    //      if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
    //        {
    //        if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    //        else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    //        e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
    //        }
    //      }
        }
      }
    //Edit torsions
    active_i=ffsys->activem[mol_i].active_t->size;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_t[ffsys->activem[mol_i].active_t->list[active_i]].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_t[ffsys->activem[mol_i].active_t->list[active_i]].atom[3]; 
      if ( (!exclude[a_i])||(!exclude[a_j]) ) //skip excl-excl
        if ( (_e=calc_atom_n14_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
          {
          if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
          else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
          e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
          }
      }
    //Remove pairs
    active_i=ffsys->activem[mol_i].active_p->size;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_p[ffsys->activem[mol_i].active_p->list[active_i]].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_p[ffsys->activem[mol_i].active_p->list[active_i]].atom[1]; 
      if ( (!exclude[a_i])||(!exclude[a_j]) ) //skip excl-excl
        if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
          {
          if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
          else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
          e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
          } 
      }
    //Summ cross molecules                    !!!!! Not debugged yet !!!!!!!!
    mol_j=mol_i;
    while (mol_j--)
      {
      if (!(ffsys->activem[mol_j].active_a))
        {
        _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (_i--)
          if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]])
            {
            _j=ffsys->mols[mol_j]->atoms->size;
            while (_j--)
              {
              a_j=ffsys->nr[mol_j]+_j; 
              if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                }
              }
            }
          else
            {
            _j=ffsys->mols[mol_j]->atoms->size;
            while (_j--)
              if (!exclude[a_j=ffsys->nr[mol_j]+_j]) //skip excl-excl
                if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                  }
            }
        }
      else
        {
        _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (_i--)
          if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]])
            {
            _l=ffsys->activem[mol_j].active_a->size;
            while (_l--)
              {
              active_j=ffsys->activem[mol_j].active_a->list[_l];
              _j=ffsys->mols[mol_j]->anchors->list[active_j].size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j];
                if (!(_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                  }   
                } 
              }
            }
          else
            {
            _l=ffsys->activem[mol_j].active_a->size;
            while (_l--)
              {
              active_j=ffsys->activem[mol_j].active_a->list[_l];
              _j=ffsys->mols[mol_j]->anchors->list[active_j].size;
              while (_j--)
                if (!exclude[a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j]]) //skip excl-excl
                  if (!(_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k, ffsys->_g[a_j][_L].i+=g_j.i, ffsys->_g[a_j][_L].j+=g_j.j, ffsys->_g[a_j][_L].k+=g_j.k;
                    }   
                 
              }
            }
        }  
      }
    }
}

//This function calculates torsions energy and its derivative
inline void calc_yff1_torss_grad(double e[N_ESLICES],t_ffsys *ffsys)
{
register unsigned int mol_id;
mol_id=ffsys->nmols; while (mol_id--) _calc_mol_tbgrad_yff1(e,&ffsys->_g[ffsys->nr[mol_id]],&ffsys->r[ffsys->nr[mol_id]],ffsys->mols[mol_id],&ffsys->activem[mol_id]);
}

//!!!!!!!!!!!!! This function can be improved by excluding the exact anchor in explicite summation of ongrid failure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculates yff1 energy on dgrid
//It uses marks made by calc_yff1_grad(t_ffsys *ffsys,t_top *top) to skip frozen atoms
inline void calc_yff1_grad_on_grid(double e[N_ESLICES],t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***B,double ***Q)
{
unsigned int _i, _j, _k, _l, _L, mol_i, mol_j, active_i, active_j, a_i, a_j;
double _e, a, b, q;
t_vec g_i, g_j, da, db, dq;

mol_i=ffsys->nmols;
while (mol_i--)
  if (!(ffsys->activem[mol_i].active_a))
    {
    _i=ffsys->mols[mol_i]->atoms->size;
    while (_i--)
      {
      a_i=ffsys->nr[mol_i]+_i;
      if ( ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,ni,nj,nk,&ffsys->r[a_i],A)))&&
           ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&b,&db,ori,sp,ni,nj,nk,&ffsys->r[a_i],B)))&&
           ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) )
        {
        //Summ on grid
        if ( (_e=top->A[(unsigned int)ffsys->a[a_i]][0]*a-top->B[(unsigned int)ffsys->a[a_i]][0]*b+ffsys->q[a_i]*q))
          { 
          multiple_vec_scalar(&da,&da,top->A[(unsigned int)ffsys->a[a_i]][0]);
          multiple_vec_scalar(&db,&db,top->B[(unsigned int)ffsys->a[a_i]][0]);
          multiple_vec_scalar(&dq,&dq,ffsys->q[a_i]);
          if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
          else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
          e[_L]+=_e, ffsys->_g[a_i][_L].i+=da.i-db.i+dq.i, ffsys->_g[a_i][_L].j+=da.j-db.j+dq.j, ffsys->_g[a_i][_L].k+=da.k-db.k+dq.k;
          }
        }
      else
        {
        //Cross summ
        for (mol_j=0; mol_j<ffsys->nmols;mol_j++)
          if ( (ffsys->activem[mol_j].active_a)) 
            {
            //Minus active
            active_j=ffsys->activem[mol_j].active_a->size;
            while (active_j--)
              {
              _j=ffsys->mols[mol_j]->anchors->list[ffsys->activem[mol_j].active_a->list[active_j]].size;
              while (_j--) 
                {
                a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[ffsys->activem[mol_j].active_a->list[active_j]].list[_j];
                _e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top);
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k;
                }
              }
            //Summ all 
            _j=ffsys->mols[mol_j]->atoms->size;
            while (_j--)
              {
              a_j=ffsys->nr[mol_j]+_j;
              _e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top);
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k;
              }
            }
        }
      }
    }
  else
    {
    //Special case - check for mark from calc_yff1_grad(t_ffsys *ffsys,t_top *top) and skip if required (it removes marks)
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch<0) //skip frozen atoms otherwise
        {
        *ffsys->rtrees[mol_i]->rbranch[active_i].rbranch=-(int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch;
        a_i=ffsys->nr[mol_i]+*ffsys->mols[mol_i]->anchors->list[active_i].list;
        if ( ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,ni,nj,nk,&ffsys->r[a_i],A)))&&
             ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&b,&db,ori,sp,ni,nj,nk,&ffsys->r[a_i],B)))&&
             ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) )
          {
          //Summ on grid
          if ( (_e=top->A[(unsigned int)ffsys->a[a_i]][0]*a-top->B[(unsigned int)ffsys->a[a_i]][0]*b+ffsys->q[a_i]*q))
            {
            multiple_vec_scalar(&da,&da,top->A[(unsigned int)ffsys->a[a_i]][0]);
            multiple_vec_scalar(&db,&db,top->B[(unsigned int)ffsys->a[a_i]][0]);
            multiple_vec_scalar(&dq,&dq,ffsys->q[a_i]);
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]+=_e, ffsys->_g[a_i][_L].i+=da.i-db.i+dq.i, ffsys->_g[a_i][_L].j+=da.j-db.j+dq.j, ffsys->_g[a_i][_L].k+=da.k-db.k+dq.k;
            }
          }
        else
          {
          //Cross summ for those forming the grid
          mol_j=ffsys->nmols;
          while (mol_j--)
            if ( (ffsys->activem[mol_j].active_a)) 
              {
              //Minus active
              _l=ffsys->activem[mol_j].active_a->size;
              while (_l--)
                {
                active_j=ffsys->activem[mol_j].active_a->list[_l];
                _j=ffsys->mols[mol_j]->anchors->list[active_j].size;
                while (_j--) 
                  {
                  a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j];
                  if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k;
                    }
                  }
                }
              //Summ all 
              _j=ffsys->mols[mol_j]->atoms->size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+_j;
                if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k;
                  }
                }
              }
          }
        }
      }
    //Summ rest as usually (copy of previous block)
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
      while (--_i)
        {
        a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
        if ( ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,ni,nj,nk,&ffsys->r[a_i],A)))&&
             ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&b,&db,ori,sp,ni,nj,nk,&ffsys->r[a_i],B)))&&
             ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) )
          {
          //Summ on grid
          if ( (_e=top->A[(unsigned int)ffsys->a[a_i]][0]*a-top->B[(unsigned int)ffsys->a[a_i]][0]*b+ffsys->q[a_i]*q))
            {
            multiple_vec_scalar(&da,&da,top->A[(unsigned int)ffsys->a[a_i]][0]);
            multiple_vec_scalar(&db,&db,top->B[(unsigned int)ffsys->a[a_i]][0]);
            multiple_vec_scalar(&dq,&dq,ffsys->q[a_i]);
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]+=_e, ffsys->_g[a_i][_L].i+=da.i-db.i+dq.i, ffsys->_g[a_i][_L].j+=da.j-db.j+dq.j, ffsys->_g[a_i][_L].k+=da.k-db.k+dq.k;
            }
          }
        else
          {
          //Cross summ for those forming the grid
          mol_j=ffsys->nmols;
          while (mol_j--)
            if ( (ffsys->activem[mol_j].active_a)) 
              {
              //Minus active
              _l=ffsys->activem[mol_j].active_a->size;
              while (_l--)
                {
                active_j=ffsys->activem[mol_j].active_a->list[_l];
                _j=ffsys->mols[mol_j]->anchors->list[active_j].size;
                while (_j--) 
                  {
                  a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j];
                  if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k;
                    }
                  }
                }
              //Summ all 
              _j=ffsys->mols[mol_j]->atoms->size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+_j;
                if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k;
                  }
                }
              }
          }
        }
      }
    }
}
//The same as previous but with skipping list in atoms
inline void calc_yff1_grad_on_grid_wexclusions(double e[N_ESLICES],unsigned char *exclude,t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***B,double ***Q)
{
unsigned int _i, _j, _k, _l, _L, mol_i, mol_j, active_i, active_j, a_i, a_j;
double _e, a, b, q;
t_vec g_i, g_j, da, db, dq;

mol_i=ffsys->nmols;
while (mol_i--)
  if (!(ffsys->activem[mol_i].active_a))
    {
    _i=ffsys->mols[mol_i]->atoms->size;
    while (_i--)
      if (!exclude[a_i=ffsys->nr[mol_i]+_i])
        {
        if ( ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,ni,nj,nk,&ffsys->r[a_i],A)))&&
             ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&b,&db,ori,sp,ni,nj,nk,&ffsys->r[a_i],B)))&&
             ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) )
          {
          //Summ on grid
          if ( (_e=top->A[(unsigned int)ffsys->a[a_i]][0]*a-top->B[(unsigned int)ffsys->a[a_i]][0]*b+ffsys->q[a_i]*q))
            { 
            multiple_vec_scalar(&da,&da,top->A[(unsigned int)ffsys->a[a_i]][0]);
            multiple_vec_scalar(&db,&db,top->B[(unsigned int)ffsys->a[a_i]][0]);
            multiple_vec_scalar(&dq,&dq,ffsys->q[a_i]);
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]+=_e, ffsys->_g[a_i][_L].i+=da.i-db.i+dq.i, ffsys->_g[a_i][_L].j+=da.j-db.j+dq.j, ffsys->_g[a_i][_L].k+=da.k-db.k+dq.k;
            }
          }
         else
          {
          //Cross summ
          for (mol_j=0; mol_j<ffsys->nmols;mol_j++)
            if ( (ffsys->activem[mol_j].active_a)) 
              {
              //Minus active
              active_j=ffsys->activem[mol_j].active_a->size;
              while (active_j--)
                {
                _j=ffsys->mols[mol_j]->anchors->list[ffsys->activem[mol_j].active_a->list[active_j]].size;
                while (_j--) 
                  {
                  a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[ffsys->activem[mol_j].active_a->list[active_j]].list[_j];
                  _e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top);
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k;
                  }
                }
              //Summ all 
              _j=ffsys->mols[mol_j]->atoms->size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+_j;
                _e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top);
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k;
                }
              }
          }
        }
    }
  else
    {
    //Special case - check for mark from calc_yff1_grad(t_ffsys *ffsys,t_top *top) and skip if required (it removes marks)
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch<0) //skip frozen atoms otherwise
        {
        *ffsys->rtrees[mol_i]->rbranch[active_i].rbranch=-(int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch;
        if (!exclude[a_i=ffsys->nr[mol_i]+*ffsys->mols[mol_i]->anchors->list[active_i].list])
          {
          if ( ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,ni,nj,nk,&ffsys->r[a_i],A)))&&
               ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&b,&db,ori,sp,ni,nj,nk,&ffsys->r[a_i],B)))&&
               ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) )
            {
            //Summ on grid
            if ( (_e=top->A[(unsigned int)ffsys->a[a_i]][0]*a-top->B[(unsigned int)ffsys->a[a_i]][0]*b+ffsys->q[a_i]*q))
              {
              multiple_vec_scalar(&da,&da,top->A[(unsigned int)ffsys->a[a_i]][0]);
              multiple_vec_scalar(&db,&db,top->B[(unsigned int)ffsys->a[a_i]][0]);
              multiple_vec_scalar(&dq,&dq,ffsys->q[a_i]);
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]+=_e, ffsys->_g[a_i][_L].i+=da.i-db.i+dq.i, ffsys->_g[a_i][_L].j+=da.j-db.j+dq.j, ffsys->_g[a_i][_L].k+=da.k-db.k+dq.k;
              }
            }
          else
            {
            //Cross summ for those forming the grid
            mol_j=ffsys->nmols;
            while (mol_j--)
              if ( (ffsys->activem[mol_j].active_a)) 
                {
                //Minus active
                _l=ffsys->activem[mol_j].active_a->size;
                while (_l--)
                  {
                  active_j=ffsys->activem[mol_j].active_a->list[_l];
                  _j=ffsys->mols[mol_j]->anchors->list[active_j].size;
                  while (_j--) 
                    {
                    a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j];
                    if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                      {
                      if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                      else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                      e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k;
                      }
                    }
                  }
                //Summ all 
                _j=ffsys->mols[mol_j]->atoms->size;
                while (_j--)
                  {
                  a_j=ffsys->nr[mol_j]+_j;
                  if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k;
                    }
                  }
                }
            }
          }
        }
      }
    //Summ rest as usually (copy of previous block)
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
      while (--_i)
        if (!exclude[a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i]])
          {
          if ( ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,ni,nj,nk,&ffsys->r[a_i],A)))&&
               ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&b,&db,ori,sp,ni,nj,nk,&ffsys->r[a_i],B)))&&
               ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) )
            {
            //Summ on grid
            if ( (_e=top->A[(unsigned int)ffsys->a[a_i]][0]*a-top->B[(unsigned int)ffsys->a[a_i]][0]*b+ffsys->q[a_i]*q))
              {
              multiple_vec_scalar(&da,&da,top->A[(unsigned int)ffsys->a[a_i]][0]);
              multiple_vec_scalar(&db,&db,top->B[(unsigned int)ffsys->a[a_i]][0]);
              multiple_vec_scalar(&dq,&dq,ffsys->q[a_i]);
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]+=_e, ffsys->_g[a_i][_L].i+=da.i-db.i+dq.i, ffsys->_g[a_i][_L].j+=da.j-db.j+dq.j, ffsys->_g[a_i][_L].k+=da.k-db.k+dq.k;
              }
            }
          else
            {
            //Cross summ for those forming the grid
            mol_j=ffsys->nmols;
            while (mol_j--)
              if ( (ffsys->activem[mol_j].active_a)) 
                {
                //Minus active
                _l=ffsys->activem[mol_j].active_a->size;
                while (_l--)
                  {
                  active_j=ffsys->activem[mol_j].active_a->list[_l];
                  _j=ffsys->mols[mol_j]->anchors->list[active_j].size;
                  while (_j--) 
                    {
                    a_j=ffsys->nr[mol_j]+ffsys->mols[mol_j]->anchors->list[active_j].list[_j];
                    if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                      {
                      if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                      else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                      e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k;
                      }
                    }
                  }
                //Summ all 
                _j=ffsys->mols[mol_j]->atoms->size;
                while (_j--)
                  {
                  a_j=ffsys->nr[mol_j]+_j;
                  if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]+=_e, ffsys->_g[a_i][_L].i+=g_i.i, ffsys->_g[a_i][_L].j+=g_i.j, ffsys->_g[a_i][_L].k+=g_i.k;
                    }
                  }
                }
            }
          }
      }
    }
}

//Note. Exclusions concept - atoms present and mobile (not in precalculated grid layer) but at the moment they are frozen.

//This function calculates yff1 energy on grid
inline double calc_ffsys_nbgrad_yff1_on_dgrid(t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *B,t_dgrid *Q)
{
unsigned int _L;
//Init energies and forces
double  e[N_ESLICES]={[0 ... N_ESLICES-1]=0.};
memset(ffsys->_g,0x0,N_ESLICES*ffsys->natoms*sizeof(t_vec));
//Summ active layers
calc_yff1_grad(e,ffsys,top);
calc_yff1_torss_grad(e,ffsys);
calc_yff1_grad_on_grid(e,ffsys,top,&A->ori,A->len.i,A->len.j,A->len.k,A->sp,A->d,B->d,Q->d);
//Summ energy
for (_L=1; _L<N_ESLICES; _L++) e[_L]+=e[_L-1];
return e[_L-1];
}
//The smae as previous but with exclusions
inline double calc_ffsys_nbgrad_yff1_on_dgrid_wexclusions(unsigned char *exclude,t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *B,t_dgrid *Q)
{
unsigned int _L;
//Init energies and forces
double  e[N_ESLICES]={[0 ... N_ESLICES-1]=0.};
memset(ffsys->_g,0x0,N_ESLICES*ffsys->natoms*sizeof(t_vec));
//Summ active layers
calc_yff1_grad_wexclusions(e,exclude,ffsys,top);
calc_yff1_torss_grad(e,ffsys);
calc_yff1_grad_on_grid_wexclusions(e,exclude,ffsys,top,&A->ori,A->len.i,A->len.j,A->len.k,A->sp,A->d,B->d,Q->d);
//Summ energy
for (_L=1; _L<N_ESLICES; _L++) e[_L]+=e[_L-1];
return e[N_ESLICES-1];
}


//-----------------------------------   I N T E R N A L   C O O R D I N A T E S    P A R T     ----------------------------------------------


//This function calculates anchor position in tree
inline void orient_anchor(t_vec *r_vec,t_tensor *R,t_tensor *R_i,double alpha,unsigned int anchor_id,t_vec *r,t_mol *mol,t_rtree *rtree)
{
unsigned int _j;
t_tensor _R;
t_vec n, tr;

tr.i=r_vec[rtree->rbranch[anchor_id].edge.vertice[1]].i;
tr.j=r_vec[rtree->rbranch[anchor_id].edge.vertice[1]].j;
tr.k=r_vec[rtree->rbranch[anchor_id].edge.vertice[1]].k;
n.i=tr.i-r_vec[rtree->rbranch[anchor_id].edge.vertice[0]].i;
n.j=tr.j-r_vec[rtree->rbranch[anchor_id].edge.vertice[0]].j;
n.k=tr.k-r_vec[rtree->rbranch[anchor_id].edge.vertice[0]].k;
multiple_vec_scalar(&n,&n,1./sqrt(calc_vec_norm(&n)));
rotate_around_uvector(&_R,&n,cos(alpha),sin(alpha));
multiple_origin_tensor_origin_tensor(R,&_R,R_i);

//Copy root itself
_j=mol->anchors->list[anchor_id].size;
while (_j--)
  if (mol->anchors->list[anchor_id].list[_j]!=rtree->rbranch[anchor_id].edge.vertice[1])
    {
    n.i=r[mol->anchors->list[anchor_id].list[_j]].i-r[rtree->rbranch[anchor_id].edge.vertice[1]].i;
    n.j=r[mol->anchors->list[anchor_id].list[_j]].j-r[rtree->rbranch[anchor_id].edge.vertice[1]].j;
    n.k=r[mol->anchors->list[anchor_id].list[_j]].k-r[rtree->rbranch[anchor_id].edge.vertice[1]].k;
    multiple_origin_tensor_origin_vec(&r_vec[mol->anchors->list[anchor_id].list[_j]],R,&n);
    r_vec[mol->anchors->list[anchor_id].list[_j]].i+=tr.i;
    r_vec[mol->anchors->list[anchor_id].list[_j]].j+=tr.j;
    r_vec[mol->anchors->list[anchor_id].list[_j]].k+=tr.k;
    }
//copy associated edges
_j=rtree->rbranch[anchor_id].nrbranch;
while (--_j)
  {
  n.i=r[rtree->rbranch[rtree->rbranch[anchor_id].rbranch[_j]].edge.vertice[1]].i-r[rtree->rbranch[anchor_id].edge.vertice[1]].i;
  n.j=r[rtree->rbranch[rtree->rbranch[anchor_id].rbranch[_j]].edge.vertice[1]].j-r[rtree->rbranch[anchor_id].edge.vertice[1]].j;
  n.k=r[rtree->rbranch[rtree->rbranch[anchor_id].rbranch[_j]].edge.vertice[1]].k-r[rtree->rbranch[anchor_id].edge.vertice[1]].k;
  multiple_origin_tensor_origin_vec(&r_vec[rtree->rbranch[rtree->rbranch[anchor_id].rbranch[_j]].edge.vertice[1]],R,&n);
  r_vec[rtree->rbranch[rtree->rbranch[anchor_id].rbranch[_j]].edge.vertice[1]].i+=tr.i;
  r_vec[rtree->rbranch[rtree->rbranch[anchor_id].rbranch[_j]].edge.vertice[1]].j+=tr.j;
  r_vec[rtree->rbranch[rtree->rbranch[anchor_id].rbranch[_j]].edge.vertice[1]].k+=tr.k;
  }
}

//This function calculates new orientation from given internal coordinates vector
unsigned int construct_internal_coords_mols(double *x,t_vec *cm,t_vec *r_vecs,t_vec *r,t_rtree *rtree,t_mol *mol,t_tensor *R)
{
unsigned int _j, _l, _id, x_id, root, anchor_id;
t_quaternion q;
t_vec n, tr;

//Stage 1. Reorient whole mol
if (!rtree->rbranch[rtree->nrbranch].nrbranch)
  {
  root=rtree->root;
  //Get derivatives for unit quaternion
  map_exp_uquaternion(&q,(t_vec*)x); //Get derivatives for unit quaternion
  calc_R_from_unit_quaternion(&R[0],&q);
  tr.i=cm->i+x[3], tr.j=cm->j+x[4], tr.k=cm->k+x[5];
  //Move in-root
  _j=mol->anchors->list[root].size;
  while (_j--)
    {
    r_vecs[mol->anchors->list[root].list[_j]].i=r[mol->anchors->list[root].list[_j]].i-cm->i;
    r_vecs[mol->anchors->list[root].list[_j]].j=r[mol->anchors->list[root].list[_j]].j-cm->j;
    r_vecs[mol->anchors->list[root].list[_j]].k=r[mol->anchors->list[root].list[_j]].k-cm->k;
    multiple_origin_tensor_origin_vec(&n,&R[0],&r_vecs[mol->anchors->list[root].list[_j]]);
    r_vecs[mol->anchors->list[root].list[_j]].i=n.i+tr.i;
    r_vecs[mol->anchors->list[root].list[_j]].j=n.j+tr.j;
    r_vecs[mol->anchors->list[root].list[_j]].k=n.k+tr.k;
    }
  //Move cross-root
  _j=rtree->rbranch[root].nrbranch;
  while (_j--)
    {
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].i=r[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].i-cm->i;
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].j=r[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].j-cm->j;
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].k=r[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].k-cm->k;
    multiple_origin_tensor_origin_vec(&n,&R[0],&r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]]);
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].i=n.i+tr.i;
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].j=n.j+tr.j;
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].k=n.k+tr.k;
    }
  //Stage c. Update _id
  x_id=6; //Store rotational quaternion
  }
else
  {
  root=rtree->nrbranch;
  //Move cross-root
  _j=rtree->rbranch[root].nrbranch;
  while (_j--)
    {
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[0]].i=r[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[0]].i,
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[0]].j=r[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[0]].j,
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[0]].k=r[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[0]].k;
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].i=r[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].i,
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].j=r[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].j,
    r_vecs[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].k=r[rtree->rbranch[rtree->rbranch[root].rbranch[_j]].edge.vertice[1]].k;
    }
  R[0][0][0]=R[0][1][1]=R[0][2][2]=1., R[0][0][1]=R[0][1][0]=R[0][0][2]=R[0][2][0]=R[0][1][2]=R[0][2][1]=0.;
  x_id=0;
  }
//Construct coordinates
_l=rtree->rbranch[root].nrbranch;
while (_l--)
  {
  anchor_id=rtree->rbranch[root].rbranch[_l];
  _id=0;
  LOOP_A:
  rtree->rbranch[anchor_id].edge.type=rtree->rbranch[anchor_id].nrbranch;
  orient_anchor(r_vecs,&R[_id+1],&R[_id],x[x_id],anchor_id,r,mol,rtree);
  _id++, x_id++;
  LOOP_B:
  while (--rtree->rbranch[anchor_id].edge.type)
    {
    anchor_id=rtree->rbranch[anchor_id].rbranch[rtree->rbranch[anchor_id].edge.type];
    goto LOOP_A;
    }
  if ( (--_id))
    {
    anchor_id=*rtree->rbranch[anchor_id].rbranch;
    goto LOOP_B;
    }
  }
return x_id;
}

//This function calculates gradient over internal coords
//NOTE. Torss in active list should be sorted previously!
void convert_ic_grad(double *g,t_vec *r_vecs,t_vec *g_vecs,t_rtree *rtree,t_mol *mol,t_tensor *R)
{
unsigned int _i, _j, _l, atom_i, root, _id, g_id=0, anchor_id;
typedef struct{
              t_vec n, o;
              unsigned int g_id;
              }t_rnode;
t_rnode *rnode;
t_vec p, dp;
rnode=(t_rnode*)R;

//Convert nonbonded grads
root= (!rtree->rbranch[rtree->nrbranch].nrbranch) ? rtree->root : rtree->nrbranch;
_l=rtree->rbranch[root].nrbranch;
while (_l--)
  {
  anchor_id=rtree->rbranch[root].rbranch[_l];
  _id=0;
  LOOP_A:
  rtree->rbranch[anchor_id].edge.type=rtree->rbranch[anchor_id].nrbranch;
  //Setup rotations
  g[g_id]=0.; 
  rnode[_id].g_id=g_id;
  rnode[_id].o.i=r_vecs[rtree->rbranch[anchor_id].edge.vertice[0]].i, 
  rnode[_id].o.j=r_vecs[rtree->rbranch[anchor_id].edge.vertice[0]].j, 
  rnode[_id].o.k=r_vecs[rtree->rbranch[anchor_id].edge.vertice[0]].k;
  rnode[_id].n.i=r_vecs[rtree->rbranch[anchor_id].edge.vertice[1]].i-rnode[_id].o.i, 
  rnode[_id].n.j=r_vecs[rtree->rbranch[anchor_id].edge.vertice[1]].j-rnode[_id].o.j, 
  rnode[_id].n.k=r_vecs[rtree->rbranch[anchor_id].edge.vertice[1]].k-rnode[_id].o.k;
  multiple_vec_scalar(&rnode[_id].n,&rnode[_id].n,1./sqrt(calc_vec_norm(&rnode[_id].n))); 
  //Summ gradients
  _i=mol->anchors->list[anchor_id].size;
  while (_i--)
    {
    atom_i=mol->anchors->list[anchor_id].list[_i];
    _j=_id;
    do{
      p.i=r_vecs[atom_i].i-rnode[_j].o.i, p.j=r_vecs[atom_i].j-rnode[_j].o.j, p.k=r_vecs[atom_i].k-rnode[_j].o.k;
      vec_vec_vmult(&dp,&rnode[_j].n,&p);
      g[rnode[_j].g_id]-=calc_vec_vec_scalar_product(&dp,&g_vecs[atom_i]);
      }while (_j--);
    }
  _id++, g_id++;
  //Continue rotations
  LOOP_B:
  while ( (--rtree->rbranch[anchor_id].edge.type))
    {
    anchor_id=rtree->rbranch[anchor_id].rbranch[rtree->rbranch[anchor_id].edge.type];
    goto LOOP_A;
    }
  if ( (--_id))
    {
    anchor_id=*rtree->rbranch[anchor_id].rbranch;
    goto LOOP_B;    
    }
  }
}

//This function calculates derivatives over rotational/translational DOF
//Note t should be of anchors.size x 2
void convert_RT_coords_grad(double *g,t_vec *v,t_vec *tr,t_vec *cm,t_vec *r_vecs,t_vec *g_vecs,t_mol *mol)
{
unsigned int _i;
double _d,__d;
t_vec r, _r;
t_quaternion qo;
t_qtensor dq;
t_tensor T, dRi, dRj, dRk, IR;

//Move to RT gradients (using quaternions)
if ((_d=v->i*v->i+v->j*v->j+v->k*v->k)>PI*PI)
  {
  _d=sqrt(_d);
  __d=1.-2.*PI*(double)((int)(_d/(2.*PI)))/_d; //Drop continius value to 0...2*PI
  if (_d*__d>PI) __d*=1.-2.*PI/(_d*__d);
  multiple_vec_scalar(&r,v,__d); //Avoid singularity
  //calc inverse matrix
  map_exp_uquaternion(&qo,&r);   //Get derivatives for unit quaternion
  dimap_exp_quaternion(&dq,&r);
  }
else
  {//Get derivatives for unit quaternion
  map_exp_uquaternion(&qo,v); 
  dimap_exp_quaternion(&dq,v);
  }
calc_dR_dquaternion(&dRi,&dRj,&dRk,&dq,&qo);
//Calculate inverse rotation matrix and vector: IR=R[-q], v2=IR.(v1-cm-tr)
qo.i=-qo.i, qo.j=-qo.j, qo.k=-qo.k;
calc_R_from_unit_quaternion(&IR,&qo);
qo.i=tr->i+cm->i, qo.j=tr->j+cm->j, qo.k=tr->k+cm->k;
//Get rotation quaternion
_i=mol->atoms->size;
while(_i--)
  {
  _r.i=r_vecs[_i].i-qo.i, _r.j=r_vecs[_i].j-qo.j, _r.k=r_vecs[_i].k-qo.k;  //  3 flops
  multiple_origin_tensor_origin_vec(&r,&IR,&_r);            // 15 flops
  multiple_origin_tensor_origin_vec((t_vec*)T[0],&dRi,&r);  // 15 flops
  multiple_origin_tensor_origin_vec((t_vec*)T[1],&dRj,&r);  // 15 flops
  multiple_origin_tensor_origin_vec((t_vec*)T[2],&dRk,&r);  // 15 flops
  multiple_origin_tensor_origin_vec(&r,&T,&g_vecs[_i]);     // 15 flops
  g[0]+=r.i, g[1]+=r.j, g[2]+=r.k, g[3]+=g_vecs[_i].i, g[4]+=g_vecs[_i].j, g[5]+=g_vecs[_i].k; // 6 flops, total 84 flops
  }
}

//This function calculates derivatives over internal + RT DOFs
//Note. g here organized as: | INTERNAL | RT |
//NB! if molecule has active list than it do NOT has RT DOFs but has several built-in rtrees
char calc_ic_grad_yff1_on_dgrid(double *_e,unsigned int n,double *x,double *g,double **G,va_list stack)
{
unsigned int _id, mol_id;
t_vec *rvecs, *cm;
t_ffsys *ffsys;
t_dgrid *A, *B, *Q;
t_top *top;
t_tensor *R;
void *rp;
va_list _stack;

va_copy(_stack,stack);

//Stage 0. Unwrap stack
ffsys=va_arg(_stack,t_ffsys*);
A=va_arg(_stack,t_dgrid*);
B=va_arg(_stack,t_dgrid*);
Q=va_arg(_stack,t_dgrid*);
rvecs=va_arg(_stack,t_vec*);
cm=va_arg(_stack,t_vec*);
R=va_arg(_stack,t_tensor*);
top=va_arg(_stack,t_top*);

//Stage 1. Construct get molecular coords
for (_id=0, mol_id=0; mol_id<ffsys->nmols; mol_id++)
  _id+=construct_internal_coords_mols(&x[_id],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->r[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);

//Stage 2. Do nonbonded linear gradient & energy calculation
memset(ffsys->_g,0x0,N_ESLICES*sizeof(t_vec)*ffsys->natoms);
rp=ffsys->r, ffsys->r=rvecs, rvecs=rp; 
*_e=calc_ffsys_nbgrad_yff1_on_dgrid(ffsys,top,A,B,Q); while(n--) g[n]=0.;
rp=rvecs, rvecs=ffsys->r, ffsys->r=rp;

//Stage 3. Convert linear into internal-coords gradient
for (_id=mol_id=0; mol_id<ffsys->nmols; mol_id++)
  {
  if (!ffsys->activem[mol_id].active_a)
    { //Stage a. Calculate RT gradients (using quaternions)
    convert_RT_coords_grad(&g[_id],(t_vec*)&x[_id+0],(t_vec*)&x[_id+3],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->g[ffsys->nr[mol_id]],ffsys->mols[mol_id]);
    _id+=6; //Store rotational quaternion
    }
  //Stage b. Calculate torsional gradients and update with soft-bonded energy
  convert_ic_grad(&g[_id],&rvecs[ffsys->nr[mol_id]],&ffsys->g[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);
  _id+=ffsys->rtrees[mol_id]->nidofs;
  }

va_end(_stack);
return TRUE;
}
//The same as previous but with exclusions
char calc_ic_grad_yff1_on_dgrid_wexclusions(double *_e,unsigned int n,double *x,double *g,double **G,va_list stack)
{
unsigned int _id, _L, mol_id, _i, _j, _k, _l, count;
unsigned char *exclude;
t_vec *rvecs, *cm, _r;
t_ffsys *ffsys;
t_dgrid *A, *B, *Q;
t_cgrid *bump;
t_top *top;
t_tensor *R;
void *rp;
va_list _stack;

va_copy(_stack,stack);

//Stage 0. Unwrap stack
exclude=va_arg(_stack,unsigned char*);
ffsys=va_arg(_stack,t_ffsys*);
A=va_arg(_stack,t_dgrid*);
B=va_arg(_stack,t_dgrid*);
Q=va_arg(_stack,t_dgrid*);
bump=va_arg(_stack,t_cgrid*);
rvecs=va_arg(_stack,t_vec*);
cm=va_arg(_stack,t_vec*);
R=va_arg(_stack,t_tensor*);
top=va_arg(_stack,t_top*);

//Stage 1. Construct get molecular coords
for (_id=0, mol_id=0; mol_id<ffsys->nmols; mol_id++)
  _id+=construct_internal_coords_mols(&x[_id],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->r[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);

//Stage 2. Calculate amount of excludes (terminate on zero)
while(n--) g[n]=0.; 
if ( (isnan(*_e)))
  {
  count=0, mol_id=ffsys->nmols;
  while (mol_id--)
    {
    //Mark dynamic layer, static layer is suggested to be '-1'
    if ( (!ffsys->activem[mol_id].active_a)||(ffsys->activem[mol_id].active_a->size==ffsys->mols[mol_id]->anchors->size) ) 
      {
      _l=0, _L=ffsys->mols[mol_id]->atoms->size; 
      while (_L--) 
        {
        _id=ffsys->nr[mol_id]+_L;
        _r.i=rvecs[_id].i-bump->ori.i, _r.j=rvecs[_id].j-bump->ori.j, _r.k=rvecs[_id].k-bump->ori.k; 
        if ( (_r.i>0.)&&(_r.j>0.)&&(_r.k>0.)&&((_i=(unsigned int)(_r.i/bump->sp))<bump->len.i-2)&&((_j=(unsigned int)(_r.j/bump->sp))<bump->len.j-2)&&((_k=(unsigned int)(_r.k/bump->sp))<bump->len.k-2) )
          {//Check bumps on grid
          if ((_r.i-=(double)_i*bump->sp)<.5*bump->sp)
            { if ( (bump->c[_i  ][_j  ][_k  ])||(bump->c[_i  ][_j+1][_k  ])||(bump->c[_i  ][_j  ][_k+1])||(bump->c[_i  ][_j+1][_k+1]) ) { exclude[_id]=1, _l++; goto NEXT_ID_0; } }
          else
            { if ( (bump->c[_i+1][_j  ][_k  ])||(bump->c[_i+1][_j+1][_k  ])||(bump->c[_i+1][_j  ][_k+1])||(bump->c[_i+1][_j+1][_k+1]) ) { exclude[_id]=1, _l++; goto NEXT_ID_0; } }
          if ((_r.j-=(double)_j*bump->sp)<.5*bump->sp)
            { if ( (bump->c[_i  ][_j  ][_k  ])||(bump->c[_i+1][_j  ][_k  ])||(bump->c[_i  ][_j  ][_k+1])||(bump->c[_i+1][_j  ][_k+1]) ) { exclude[_id]=1, _l++; goto NEXT_ID_0; } }
          else
            { if ( (bump->c[_i  ][_j+1][_k  ])||(bump->c[_i+1][_j+1][_k  ])||(bump->c[_i  ][_j+1][_k+1])||(bump->c[_i+1][_j+1][_k+1]) ) { exclude[_id]=1, _l++; goto NEXT_ID_0; } }
          if ((_r.k-=(double)_k*bump->sp)<.5*bump->sp)
            { if ( (bump->c[_i  ][_j  ][_k  ])||(bump->c[_i+1][_j  ][_k  ])||(bump->c[_i  ][_j+1][_k  ])||(bump->c[_i+1][_j+1][_k  ]) ) { exclude[_id]=1, _l++; goto NEXT_ID_0; } }
          else
            { if ( (bump->c[_i  ][_j  ][_k+1])||(bump->c[_i+1][_j  ][_k+1])||(bump->c[_i  ][_j+1][_k+1])||(bump->c[_i+1][_j+1][_k+1]) ) { exclude[_id]=1, _l++; goto NEXT_ID_0; } }
          } 
        else
          {//Check bumps manually
          _i=ffsys->natoms; 
          while (--_i!=_id) 
            if ( (exclude[_i]==(unsigned char)-1)&&(fabs(rvecs[_i].i-_r.i)<2.*bump->sp)&&(fabs(rvecs[_i].j-_r.j)<2.*bump->sp)&&(fabs(rvecs[_i].k-_r.k)<2.*bump->sp) ) { exclude[_id]=+1, _l++; goto NEXT_ID_0; }
          while (_i--)
            if ( (exclude[_i]==(unsigned char)-1)&&(fabs(rvecs[_i].i-_r.i)<2.*bump->sp)&&(fabs(rvecs[_i].j-_r.j)<2.*bump->sp)&&(fabs(rvecs[_i].k-_r.k)<2.*bump->sp) ) { exclude[_id]=+1, _l++; goto NEXT_ID_0; }
          }
        exclude[_id]=0;
        NEXT_ID_0: ; 
        } 
      if (_l==ffsys->mols[mol_id]->atoms->size) { va_end(_stack); return FALSE; }
      else count+=_l;
      } 
    else
      {
      _l=ffsys->activem[mol_id].active_a->size;
      while (_l--)
        {
        _L=ffsys->mols[mol_id]->anchors->list[ffsys->activem[mol_id].active_a->list[_l]].size;
        while (_L--)
          {
          _id=ffsys->nr[mol_id]+ffsys->mols[mol_id]->anchors->list[ffsys->activem[mol_id].active_a->list[_l]].list[_L];
          _r.i=rvecs[_id].i-bump->ori.i, _r.j=rvecs[_id].j-bump->ori.j, _r.k=rvecs[_id].k-bump->ori.k; 
          if ( (_r.i>0.)&&(_r.j>0.)&&(_r.k>0.)&&((_i=(unsigned int)(_r.i/bump->sp))<bump->len.i-2)&&((_j=(unsigned int)(_r.j/bump->sp))<bump->len.j-2)&&((_k=(unsigned int)(_r.k/bump->sp))<bump->len.k-2) )
            {//Check bumps on grid
            if ((_r.i-=(double)_i*bump->sp)<.5*bump->sp)
              { if ( (bump->c[_i  ][_j  ][_k  ])||(bump->c[_i  ][_j+1][_k  ])||(bump->c[_i  ][_j  ][_k+1])||(bump->c[_i  ][_j+1][_k+1]) ) { exclude[_id]=1, count++; goto NEXT_ID_1; } }
            else
              { if ( (bump->c[_i+1][_j  ][_k  ])||(bump->c[_i+1][_j+1][_k  ])||(bump->c[_i+1][_j  ][_k+1])||(bump->c[_i+1][_j+1][_k+1]) ) { exclude[_id]=1, count++; goto NEXT_ID_1; } }
            if ((_r.j-=(double)_j*bump->sp)<.5*bump->sp)
              { if ( (bump->c[_i  ][_j  ][_k  ])||(bump->c[_i+1][_j  ][_k  ])||(bump->c[_i  ][_j  ][_k+1])||(bump->c[_i+1][_j  ][_k+1]) ) { exclude[_id]=1, count++; goto NEXT_ID_1; } }
            else
              { if ( (bump->c[_i  ][_j+1][_k  ])||(bump->c[_i+1][_j+1][_k  ])||(bump->c[_i  ][_j+1][_k+1])||(bump->c[_i+1][_j+1][_k+1]) ) { exclude[_id]=1, count++; goto NEXT_ID_1; } }
            if ((_r.k-=(double)_k*bump->sp)<.5*bump->sp)
              { if ( (bump->c[_i  ][_j  ][_k  ])||(bump->c[_i+1][_j  ][_k  ])||(bump->c[_i  ][_j+1][_k  ])||(bump->c[_i+1][_j+1][_k  ]) ) { exclude[_id]=1, count++; goto NEXT_ID_1; } }
            else
              { if ( (bump->c[_i  ][_j  ][_k+1])||(bump->c[_i+1][_j  ][_k+1])||(bump->c[_i  ][_j+1][_k+1])||(bump->c[_i+1][_j+1][_k+1]) ) { exclude[_id]=1, count++; goto NEXT_ID_1; } }
            } 
          else
            {//Check bumps manually
            _i=ffsys->natoms; 
            while (--_i!=_id) 
              if ( (exclude[_i]==(unsigned char)-1)&&(fabs(rvecs[_i].i-_r.i)<2.*bump->sp)&&(fabs(rvecs[_i].j-_r.j)<2.*bump->sp)&&(fabs(rvecs[_i].k-_r.k)<2.*bump->sp) ) { exclude[_id]=+1, count++; goto NEXT_ID_1; }
            while (_i--)
              if ( (exclude[_i]==(unsigned char)-1)&&(fabs(rvecs[_i].i-_r.i)<2.*bump->sp)&&(fabs(rvecs[_i].j-_r.j)<2.*bump->sp)&&(fabs(rvecs[_i].k-_r.k)<2.*bump->sp) ) { exclude[_id]=+1, count++; goto NEXT_ID_1; }
            }
          exclude[_id]=0;
          NEXT_ID_1: ; 
          }
        }
      }
    }
  if (!count) { goto END; } //Exit on zero excludes (gradients vector is zero now!)
  }
  
//Stage 3. Do nonbonded linear gradient & energy calculation
rp=ffsys->r, ffsys->r=rvecs, rvecs=rp; 
memset(ffsys->_g,0x0,sizeof(t_vec)*N_ESLICES*ffsys->natoms);
*_e=calc_ffsys_nbgrad_yff1_on_dgrid_wexclusions(exclude,ffsys,top,A,B,Q);
rp=rvecs, rvecs=ffsys->r, ffsys->r=rp;

//Stage 4. Convert linear into internal-coords gradient
_id=ffsys->natoms; 
while (_id--) 
  {
  ffsys->g[_id].i=ffsys->g[_id].j=ffsys->g[_id].k=0.;
  if (!exclude[_id]) 
    for (_L=0;_L<N_ESLICES;_L++) { ffsys->g[_id].i+=ffsys->_g[_id][_L].i, ffsys->g[_id].j+=ffsys->_g[_id][_L].j, ffsys->g[_id].k+=ffsys->_g[_id][_L].k; }
  }
for (_id=mol_id=0; mol_id<ffsys->nmols; mol_id++)
  {
  if (!ffsys->activem[mol_id].active_a)
    { //Stage a. Calculate RT gradients (using quaternions)
    convert_RT_coords_grad(&g[_id],(t_vec*)&x[_id+0],(t_vec*)&x[_id+3],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->g[ffsys->nr[mol_id]],ffsys->mols[mol_id]);
    _id+=6; //Store rotational quaternion
    }
  //Stage b. Calculate torsional gradients and update with soft-bonded energy
  convert_ic_grad(&g[_id],&rvecs[ffsys->nr[mol_id]],&ffsys->g[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);
  _id+=ffsys->rtrees[mol_id]->nidofs;
  }

END: va_end(_stack);
return TRUE;
}


//--------------------------     L I N E A R     C O O R D I N A T E S    P A R T   -----------------------------------------------------

//    G E N E R A L - P U R P O S E     M I N I M I Z A T I O N     R O U T I N E

//This function initializes internal coordinates minimizer
char initialize_ic_minimizer(unsigned int method,unsigned int max_nwats,t_tensor **R,t_vec **cm,t_vec **rdocked_m,unsigned char **exclude,double ***x,double ***g,double ***G,double **p,
                             unsigned int *naatoms_om,unsigned int **id_om,unsigned int *nidofs_om,unsigned int *max_nidofs_om,t_vec **rdocked_x,t_rtree **mrtrees,t_ffsys *ffsys)
{
unsigned int _i, _j, _k, _root, _id;

//Gather statistics
*nidofs_om=*max_nidofs_om=0, *naatoms_om=0, _i=ffsys->nmols; 
while (_i--)
  {
  //Scan mrtrees[_i]
  if ( (mrtrees[_i]))
    {
    if ( (mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].edge.type=mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].nrbranch))
      {
      *nidofs_om+=mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].edge.type;
      while (mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].edge.type--)
        {
        _k=1, _root=_id=mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].rbranch[mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].edge.type], (*naatoms_om)+=ffsys->mols[_i]->anchors->list[_id].size;
        mrtrees[_i]->rbranch[_id].edge.type=mrtrees[_i]->rbranch[_id].nrbranch; 
        MRTREE_LOOP_0:
        while ( (--mrtrees[_i]->rbranch[_id].edge.type))
          {
          (*nidofs_om)++;
          _id=mrtrees[_i]->rbranch[_id].rbranch[mrtrees[_i]->rbranch[_id].edge.type], mrtrees[_i]->rbranch[_id].edge.type=mrtrees[_i]->rbranch[_id].nrbranch;
          (*naatoms_om)+=ffsys->mols[_i]->anchors->list[_id].size, _k++; //for static molecules only
          }
        if (_k>(*max_nidofs_om)) (*max_nidofs_om)=_k;
        _k--;
        if (_id!=_root) { _id=*mrtrees[_i]->rbranch[_id].rbranch; goto MRTREE_LOOP_0; }
        }
      }
    else *nidofs_om=mrtrees[_i]->nidofs+6; //All rtree + RT freedom
    }
  }
(*nidofs_om)+=max_nwats*6; //Reserve space for waters

//Alloc memory
if (!((*id_om)=(unsigned int*)malloc(sizeof(unsigned int)*(*naatoms_om))))     { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return FALSE; }
if (!((*rdocked_m)=(t_vec*)malloc(sizeof(t_vec)*(*naatoms_om))))               { LABEL_MEMORY_ERROR_1: free(*id_om); *id_om=0x0; goto LABEL_MEMORY_ERROR_0; } 
if (!((*rdocked_x)=(t_vec*)malloc(sizeof(t_vec)*ffsys->natoms+max_nwats*3)))   { LABEL_MEMORY_ERROR_2: free(*rdocked_m); *rdocked_m=0x0; goto LABEL_MEMORY_ERROR_1; } 
if (!((*cm)=(t_vec*)malloc(sizeof(t_vec)*ffsys->nmols)))                       { LABEL_MEMORY_ERROR_3: free(*rdocked_x); *rdocked_x=0x0; goto LABEL_MEMORY_ERROR_2; } 
if (!((*R)=(t_tensor*)malloc(sizeof(t_tensor)*(*max_nidofs_om+1))))            { LABEL_MEMORY_ERROR_4: free(*cm); *cm=0x0; goto LABEL_MEMORY_ERROR_3; } 
if (!((*exclude)=(unsigned char*)malloc(sizeof(unsigned char)*ffsys->natoms))) { LABEL_MEMORY_ERROR_5: free(*R); *R=0x0; goto LABEL_MEMORY_ERROR_4; } 

//Fill ids and rrs
_k=0, _i=ffsys->nmols;
while (_i--)
  //Scann rtrees
  if ( (mrtrees[_i]))
    {
    //Scan mrtrees[_i]
    mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].edge.type=mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].nrbranch;
    while (mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].edge.type--)
      {
      _root=_id=mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].rbranch[mrtrees[_i]->rbranch[mrtrees[_i]->nrbranch].edge.type];
      mrtrees[_i]->rbranch[_id].edge.type=mrtrees[_i]->rbranch[_id].nrbranch;
      for (_j=0; _j<ffsys->mols[_i]->anchors->list[_id].size; _j++, _k++) { (*id_om)[_k]=ffsys->nr[_i]+ffsys->mols[_i]->anchors->list[_id].list[_j]; }
      RTREE_LOOP_1:
      while ( (--mrtrees[_i]->rbranch[_id].edge.type))
        {
        _id=mrtrees[_i]->rbranch[_id].rbranch[mrtrees[_i]->rbranch[_id].edge.type], mrtrees[_i]->rbranch[_id].edge.type=mrtrees[_i]->rbranch[_id].nrbranch;
        for (_j=0; _j<ffsys->mols[_i]->anchors->list[_id].size; _j++, _k++) { (*id_om)[_k]=ffsys->nr[_i]+ffsys->mols[_i]->anchors->list[_id].list[_j]; }
        }
      if (_id!=_root) { _id=*mrtrees[_i]->rbranch[_id].rbranch; goto RTREE_LOOP_1; }
      }
    }

//Alloc memory
_i=(*nidofs_om)+max_nwats*6;
switch (method)
  {
  case MINIMIZE_STEEPEST :
  case MINIMIZE_POLAKRB  :
  case MINIMIZE_LBFGS    : {
                           if (!((*p)=(double*)malloc(sizeof(double)*_i*5+sizeof(double*)*4))) { free(*exclude), *exclude=0x0; goto LABEL_MEMORY_ERROR_5; }
                           *x=(void*)(*p)+sizeof(double)*_i, (*x)[0]=(void*)*x+sizeof(double*)*2, (*x)[1]=(*x)[0]+_i, *g=(void*)(*x)[1]+sizeof(double)*_i, (*g)[0]=(void*)(*g)+sizeof(double*)*2, (*g)[1]=(*g)[0]+_i;
                           *G=0x0;
                           break;
                           }
  case MINIMIZE_NEWTON   : {
                           if (!((*p)=(double*)malloc((sizeof(double)*_i+sizeof(double*))*_i+sizeof(double)*_i*5+sizeof(double*)*4))) { free(*exclude), *exclude=0x0; goto LABEL_MEMORY_ERROR_5; }
                           *x=(void*)(*p)+sizeof(double)*_i, (*x)[0]=(void*)(*x)+sizeof(double*)*2, (*x)[1]=(*x)[0]+_i, *g=(void*)(*x)[1]+sizeof(double)*_i, (*g)[0]=(void*)(*g)+sizeof(double*)*2, (*g)[1]=(*g)[0]+_i;
                           for (*G=(void*)(*g)[1]+sizeof(double)*_i, (*G)[0]=(void*)(*G)+sizeof(double*)*_i, _j=1;_j<_i;_j++) (*G)[_j]=(*G)[_i-1]+_i;
                           break;
                           }
  default                : { ylib_errno=YERROR_NIMPLEMENTED; free(*id_om); *id_om=0x0; free(*rdocked_m); *rdocked_m=0x0; free(*rdocked_x); *rdocked_x=0x0; free(*cm); *cm=0x0; free(*R); *R=0x0; return FALSE; }
  }
return TRUE;
}


//This function setsup exclisions for minimizing
void setup_ic_minimizer_wexclusions(unsigned int *n_excl,unsigned char *exclude,t_ffsys *ffsys)
{
unsigned int mol_id, _id, _i, _j;

//Stage 0. Init
*n_excl=0, memset(exclude,(unsigned char)-1,sizeof(unsigned char)*ffsys->natoms);

//Stage 1. Scan mols
mol_id=ffsys->nmols;
while (mol_id--)
  if (!(ffsys->activem[mol_id].active_a))
    {
    (*n_excl)+=ffsys->rtrees[mol_id]->nidofs+6;
    memset(&exclude[ffsys->nr[mol_id]],0x0,sizeof(unsigned char)*ffsys->mols[mol_id]->atoms->size);
    }
  else if ( (_i=ffsys->activem[mol_id].active_a->size)) //Scann rtrees
         {
         while (_i--)
           {
           _id=ffsys->activem[mol_id].active_a->list[_i], _j=ffsys->mols[mol_id]->anchors->list[_id].size;
           while (_j--) exclude[ffsys->mols[mol_id]->anchors->list[_id].list[_j]]=0;
           }
         _i=ffsys->rtrees[mol_id]->rbranch[ffsys->rtrees[mol_id]->nrbranch].nrbranch, (*n_excl)+=_i;
         while (_i--)
           {
           _id=ffsys->rtrees[mol_id]->rbranch[ffsys->rtrees[mol_id]->nrbranch].rbranch[_i];
           ffsys->rtrees[mol_id]->rbranch[_id].edge.type=ffsys->rtrees[mol_id]->rbranch[_id].nrbranch;
           SCANN_RTREE: while (--ffsys->rtrees[mol_id]->rbranch[_id].edge.type)
                          {
                          _id=ffsys->rtrees[mol_id]->rbranch[_id].rbranch[ffsys->rtrees[mol_id]->rbranch[_id].edge.type];
                          ffsys->rtrees[mol_id]->rbranch[_id].edge.type=ffsys->rtrees[mol_id]->rbranch[_id].nrbranch;
                          (*n_excl)++;                
                          }
           if (_id!=ffsys->rtrees[mol_id]->rbranch[ffsys->rtrees[mol_id]->nrbranch].rbranch[_i]) { _id=*ffsys->rtrees[mol_id]->rbranch[_id].rbranch; goto SCANN_RTREE; }
           }
         }
}

//This function minimizes system in internal coordinates on the grid
char optimize_ic_ffsys_on_grid(double *e,unsigned int method,unsigned int nsteps,double tol,unsigned int n,double **x,double **g,double **G,double *p,unsigned int lbfgs_m,t_ffsys *ffsys,t_dgrid *A,t_dgrid *B,t_dgrid *Q,t_vec *rvecs,t_vec *cm,t_tensor *R,t_top *top)
{
unsigned int mol_i, _i, _j;
double m, _m;
//Calculate cms
mol_i=ffsys->nmols;
while (mol_i--) 
  {
  cm[mol_i].i=cm[mol_i].j=cm[mol_i].k=0.; 
  if (!(ffsys->activem[mol_i].active_a))
    {
    m=0., _i=ffsys->mols[mol_i]->anchors->list[ffsys->rtrees[mol_i]->root].size;
    while (_i--) 
      { 
      _j=ffsys->mols[mol_i]->anchors->list[ffsys->rtrees[mol_i]->root].list[_i], 
      _m=top->ff_a[(unsigned int)ffsys->a[_j]].mass, m+=_m, cm[mol_i].i+=_m*ffsys->r[_j].i, cm[mol_i].j+=_m*ffsys->r[_j].j, cm[mol_i].k+=_m*ffsys->r[_j].k;
      }
    cm[mol_i].i/=m, cm[mol_i].j/=m, cm[mol_i].k/=m;
    }
  }
//Init x
memset(x[0],0x0,sizeof(double)*n);
//Do minimization
switch (method)
  {
  case MINIMIZE_STEEPEST : { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
  case MINIMIZE_POLAKRB  : { 
    return polak_ribiere(e,nsteps,tol,PI/2.,n,x,g,G,p,calc_ic_grad_yff1_on_dgrid,line_search_square_fapproximation,0x0,ffsys,A,B,Q,rvecs,cm,R,top);
    }
  case MINIMIZE_LBFGS    : {
    return lbfgs(e,nsteps,tol,PI/2.,n,x,g,G,p,lbfgs_m,calc_ic_grad_yff1_on_dgrid,line_search_square_fapproximation,0x0,ffsys,A,B,Q,rvecs,cm,R,top);
    }
  case MINIMIZE_NEWTON   : { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
  default                : { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
  }
}

//This function minimizes system in internal coordinates on the grid
char optimize_ic_ffsys_on_grid_wexclusions(double *e,unsigned int nsteps,double tol,unsigned int n,double **x,double **g,double **G,double *p,unsigned char *exclude,t_ffsys *ffsys,t_dgrid *A,t_dgrid *B,t_dgrid *Q,t_cgrid *bump,t_vec *rvecs,t_vec *cm,t_tensor *R,t_top *top)
{
unsigned int mol_i, _i, _j;
double m, _m;
char status;
//Calculate cms
mol_i=ffsys->nmols;
while (mol_i--) 
  {
  cm[mol_i].i=cm[mol_i].j=cm[mol_i].k=0.; 
  if (!(ffsys->activem[mol_i].active_a))
    {
    m=0., _i=ffsys->mols[mol_i]->anchors->list[ffsys->rtrees[mol_i]->root].size;
    while (_i--) 
      { 
      _j=ffsys->mols[mol_i]->anchors->list[ffsys->rtrees[mol_i]->root].list[_i], 
      _m=top->ff_a[(unsigned int)ffsys->a[_j]].mass, m+=_m, cm[mol_i].i+=_m*ffsys->r[_j].i, cm[mol_i].j+=_m*ffsys->r[_j].j, cm[mol_i].k+=_m*ffsys->r[_j].k;
      }
    cm[mol_i].i/=m, cm[mol_i].j/=m, cm[mol_i].k/=m;
    }
  }
//Init x
memset(x[0],0x0,sizeof(double)*n);
//Do iterative minimization
return steepest_iterative(e,nsteps,tol,PI/2.,n,x,g,G,p,calc_ic_grad_yff1_on_dgrid_wexclusions,line_search_square_fapproximation,0x0,exclude,ffsys,A,B,Q,bump,rvecs,cm,R,top);
}

