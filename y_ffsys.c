#include "y_ffsys.h"


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



//This function compile activem block
//Note. It DO NOT copy active_a list, just a pointer to it.
char create_activem(t_activem *activem,t_list *active_a,t_mol *mol)
{
unsigned int _k, _l;
//Init activem
activem->active_a=activem->active_b=activem->active_g=activem->active_i=activem->active_d=activem->active_p=0x0;
if ( ( (active_a))&&( (active_a->size)) )
  {
  //Copy active_a list
  activem->active_a=active_a;
  //Mark active atoms
  _k=active_a->size; while (_k--) { _l=mol->anchors->list[active_a->list[_k]].size; while (_l--) { mol->atoms->list[mol->anchors->list[active_a->list[_k]].list[_l]]=(unsigned int)-mol->atoms->list[mol->anchors->list[active_a->list[_k]].list[_l]]; } }
  //Create active_b, active_g, active_i and active_d lists
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
      else { if (!(activem->active_i=alloc_list(0x0))) goto MEMORY_ERROR_2; else activem->active_i->size=0; }
      _l=0, _k=mol->size_d; 
      while (_k--)
        if ( ((int)mol->atoms->list[mol->ff_d[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_d[_k].atom[1]]<0)||((int)mol->atoms->list[mol->ff_d[_k].atom[2]]<0)||((int)mol->atoms->list[mol->ff_d[_k].atom[3]]<0) ) _l++;
      if (_l)
        {
        if (!(activem->active_d=alloc_list(_l))) { MEMORY_ERROR_3: free(activem->active_i); activem->active_i=0x0; goto MEMORY_ERROR_2; }
        _l=0, _k=mol->size_d; 
        while (_k--)
          if ( ((int)mol->atoms->list[mol->ff_d[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_d[_k].atom[1]]<0)||((int)mol->atoms->list[mol->ff_d[_k].atom[2]]<0)||((int)mol->atoms->list[mol->ff_d[_k].atom[3]]<0) )
            activem->active_d->list[_l++]=_k;
        }
      else { if (!(activem->active_d=alloc_list(0x0))) goto MEMORY_ERROR_3; else activem->active_d->size=0; } 
      }
    else
      {
      if (!(activem->active_g=alloc_list(0x0))) goto MEMORY_ERROR_1; else activem->active_g->size=0;
      if (!(activem->active_i=alloc_list(0x0))) goto MEMORY_ERROR_2; else activem->active_i->size=0;
      if (!(activem->active_d=alloc_list(0x0))) goto MEMORY_ERROR_3; else activem->active_d->size=0;
      }
    }
  else
    {
    if (!(activem->active_b=alloc_list(0x0))) goto MEMORY_ERROR_0; else activem->active_b->size=0; 
    if (!(activem->active_g=alloc_list(0x0))) goto MEMORY_ERROR_1; else activem->active_g->size=0;
    if (!(activem->active_i=alloc_list(0x0))) goto MEMORY_ERROR_2; else activem->active_i->size=0;
    if (!(activem->active_d=alloc_list(0x0))) goto MEMORY_ERROR_3; else activem->active_d->size=0;
    }
  _l=0;_k=mol->size_p; 
  while(_k--)
    if ( ((int)mol->atoms->list[mol->ff_p[_k].atom[0]]<0)||((int)mol->atoms->list[mol->ff_p[_k].atom[1]]<0) ) _l++;
  if (_l)
    {
    if (!(activem->active_p=alloc_list(_l))) { free(activem->active_d); activem->active_d=0x0; goto MEMORY_ERROR_3; } 	  
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
      _l=activem->active_d->size;
      while (_l--)
        if ( ( (mol->ff_d[activem->active_d->list[_l]].atom[0]==_id)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[2]]>0)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[3]]>0) ) ||
             ( (mol->ff_d[activem->active_d->list[_l]].atom[1]==_id)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[2]]>0)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[3]]>0) ) ||
             ( (mol->ff_d[activem->active_d->list[_l]].atom[2]==_id)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[3]]>0) ) ||
             ( (mol->ff_d[activem->active_d->list[_l]].atom[3]==_id)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[0]]>0)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[1]]>0)&&((int)mol->atoms->list[mol->ff_d[activem->active_d->list[_l]].atom[2]]>0) )  )
          activem->active_d->list[_l]=activem->active_d->list[--activem->active_d->size];
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
if (activem->active_d) { free(activem->active_d); activem->active_d=0x0; }
if (activem->active_p) { free(activem->active_p); activem->active_p=0x0; }
}
 
//This function reads stack of activems from hdd
t_activem *read_activem(FILE *in,unsigned int *nmols)
{
t_activem *activem;
register unsigned int _i;
unsigned int j;

if ( (fread(&j,sizeof(unsigned int),0x1,in)!=0x1)||(j!=Y_MAGIC) ) { LABEL_IO_ERROR_0: ylib_errno=YERROR_IO; return FALSE; }
if (fread(&j,sizeof(unsigned int),0x1,in)!=0x1) goto LABEL_IO_ERROR_0;
(*nmols)=j;
if (!(activem=(t_activem*)malloc(sizeof(t_activem)*(*nmols)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
for (_i=0;_i<(*nmols);_i++)
  {
  if (fread(&j,sizeof(unsigned int),0x1,in)!=0x1) { LABEL_IO_ERROR_1: while (_i--) free_activem(&activem[_i]); free(activem); goto LABEL_IO_ERROR_0; }
  if (j)
    {
    if (j==Y_MAGIC)
      {
      if ( (!(activem[_i].active_a=read_list(in)))||(!(activem[_i].active_b=read_list(in)))||(!(activem[_i].active_g=read_list(in)))||
           (!(activem[_i].active_i=read_list(in)))||(!(activem[_i].active_d=read_list(in)))||(!(activem[_i].active_p=read_list(in))) ) goto LABEL_IO_ERROR_1;
      }
    else goto LABEL_IO_ERROR_1;
    }
  else memset(&activem[_i],0x0,sizeof(t_activem));
  }
return activem;
}

//This function writes stack of activems to hdd
char write_activem(FILE *out,unsigned int nmols,t_activem *activem)
{
register unsigned int _i;
unsigned int j;
j=Y_MAGIC;
if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) { LABEL_IO_ERROR: ylib_errno=YERROR_IO; return FALSE; }
j=nmols;
if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;
for (_i=0;_i<nmols;_i++)
  if (activem[_i].active_a)
    {
    j=Y_MAGIC;
    if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;	
    if ( (!write_list(out,activem[_i].active_a))||(!write_list(out,activem[_i].active_b))||(!write_list(out,activem[_i].active_g))||
         (!write_list(out,activem[_i].active_i))||(!write_list(out,activem[_i].active_d))||(!write_list(out,activem[_i].active_p)) ) 
      return FALSE;
    }
  else
    {  
    j=0;
    if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;  
    }
return TRUE;
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
      if (ffsys->activem[_i].active_a) { free_activem(&ffsys->activem[_i]); ffsys->activem[_i].active_a=0x0; }
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
      if (!(ffsys->rtrees[i]=read_rtree(in))) return FALSE;
      }
    else goto LABEL_IO_ERROR;
    }
  else ffsys->rtrees[i]=0x0;
  }
//Read global massives
if ( (fread(&i,sizeof(unsigned int),0x1,in)!=0x1)||(i!=Y_MAGIC) ) goto LABEL_IO_ERROR;
if (fread(ffsys->nr,sizeof(unsigned int),ffsys->nmols,in)!=ffsys->nmols) goto LABEL_IO_ERROR;
if ( (!(ffsys->a=(char*)malloc(sizeof(char)*ffsys->natoms)))||(!(ffsys->q=(double*)malloc(sizeof(double)*ffsys->natoms)))||
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
  if ( (i==ffsys->nmols-1)||(ffsys->mols[i]!=ffsys->mols[i+1]) )
    {
    if (!(write_ymol(out,ffsys->mols[i]))) return FALSE;
    if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;
    else j=0;
    }	
  }
for (i=0;i<ffsys->nmols;i++)
  {
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
//write global massives  
j=Y_MAGIC;
if (fwrite(&j,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_IO_ERROR;
if (fwrite(ffsys->nr,sizeof(unsigned int),ffsys->nmols,out)!=ffsys->nmols) goto LABEL_IO_ERROR;
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
    ffsys->activem[ffsys->nmols].active_d=activem->active_d;
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
t_mol *sol=0x0;
t_rtree *sol_rtree=0x0;

//Stage I. Create sol template 
if (!(sol=construct_wat(top))) return FALSE;

//Stage I.3. Create rtree
if (!(sol_rtree=build_rtree(0x0,sol))) { LABEL_FAILURE_0: free_mol(sol); sol=0x0; return FALSE; }

//Stage II. Integrate waters into ffsys
if (!(ffsys_add_mol(nsols,sol,0x0,sol_rtree,ffsys))) { free(sol_rtree); sol_rtree=0x0; goto LABEL_FAILURE_0; }

return TRUE;
}

//-----------------------------   Y F F 0    P A R T   ---------------------------------------------

// Inspirated by 
// Gehlhaar D, Verkhivker G., Rejto P., Sherman C., Fogel D., Fogel L., Freer S., Chem Biol 1995, 2, 317.
// and 
// Fuhrmaann J., Rurainski A., Lenhof H. and Neumann D., J.Comput.Chem., 2009, 30, 1371-1378.

//Setup yff0.

//   ^ vdW                                    ^ Q
// F +.                                       |
//   ||                                       |       sQ              tQ     YFF0_Rc
//   | .                                      |       |                |       |
//   | |                                    F |       v                |       |
//   |  .                                  +- +-------.                |       |
//   |  |                                   2 |       |`-._            v       | 
//   |   .                                    |       |    ```-----....        |
//   |  A|   B       C    D                   |       |           C  D ``--..  v
// --+----.--.---------.-.._______---->     --+----.--+------------.---.-----``____---->        
// 0 |    ^\            / ^                  0|   A  B|                 _..--``       
// E +    | `__________`  |                   |       |      ____...--``
//   |    | ^          ^  |                 F |       |.-````   
//   |    | |          |  |               - - +-------`
//   |    A  B         C  D  } +/-Delta {   2 |  
//
//This function setups YFF0
inline void setup_yff0(t_top *top)
{
top->yff0_sA[0]=3.2;
top->yff0_sB[0]=3.6;
top->yff0_sC[0]=4.5;
top->yff0_sD[0]=YFF0_Rb-YFF0_sD;
top->yff0_sE=-0.4;
top->yff0_sF=4.*YFF0_COULOMB_K*COULOMB_K/(top->yff0_sB[0]-YFF0_sD); //~96.0
// sQ[0]=sB[0], { sQ[1]=1/(B+d), sQ[2]=0, sQ[3]=3*(1/(B+d)-1/(B-d))+(B-A)/(B+d)^2, sQ[4]=-2*(1/(B+d)-1/(B-d))-(B-A)/(B+d)^2 } x Q
top->yff0_sQ[0]=+1./(top->yff0_sB[0]-YFF0_sD),
top->yff0_sQ[1]=0., 
top->yff0_sQ[2]=+3.*(1./(top->yff0_sB[0]+YFF0_sD)-1./(top->yff0_sB[0]-YFF0_sD))+2.*YFF0_sD/sqrd(top->yff0_sB[0]+YFF0_sD),
top->yff0_sQ[3]=-2.*(1./(top->yff0_sB[0]+YFF0_sD)-1./(top->yff0_sB[0]-YFF0_sD))-2.*YFF0_sD/sqrd(top->yff0_sB[0]+YFF0_sD);
// sD[0]=D, sD[1]=E/(D-C), sD[2]=-2*d*E/(D-C), sD[3]=d*E/(D-C)
top->yff0_sD[1]=+YFF0_sD*top->yff0_sE/(top->yff0_sD[0]-top->yff0_sC[0]),
top->yff0_sD[2]=-2.*top->yff0_sD[1],
top->yff0_sD[3]=top->yff0_sD[1];
// sC[0]=C, sC[1]=E, sC[2]=0, sC[3]=-E/(D-C)
top->yff0_sC[1]=top->yff0_sE,
top->yff0_sC[2]=0.,
top->yff0_sC[3]=-top->yff0_sD[1];
// sB[0]=B, sB[1]=E-d*E/(B-A), sB[2]=2*d*E/(B-A), sB[3]=-d*E/(B-A)
top->yff0_sB[3]=YFF0_sD*top->yff0_sE/(top->yff0_sA[0]-top->yff0_sB[0]),
top->yff0_sB[2]=-2.*top->yff0_sB[3],
top->yff0_sB[1]=top->yff0_sE+top->yff0_sB[3];
// sA[0]=A, sA[1]=d*F/A, sA[2]=-2*d*F/A, sA[3]=d*E/(B-A)+d*F/A
top->yff0_sA[1]=+YFF0_sD*top->yff0_sF/top->yff0_sA[0],
top->yff0_sA[2]=-2.*top->yff0_sA[1],
top->yff0_sA[3]=top->yff0_sA[1]-top->yff0_sB[3];
}

//------------------------------------    N O N B O N D E D   P R I M I T I V E S   ------------------------------------------

//This function calculates energy and force of nonbonded interatons in YFF0. The flag defines if energy(forces) adds of subtracts
//Note g_vecs and r_vecs migh be ffsys->g and ffsys->r correspondingly
inline double calc_atom__nb_yff0(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top)
{
register double _d, _rr, _r, _x, _Q;
t_vec _u;

_u.i=r_i->i-r_j->i, _u.j=r_i->j-r_j->j, _u.k=r_i->k-r_j->k, _rr=_u.i*_u.i+_u.j*_u.j+_u.k*_u.k+SMALL2; //To keep divisibility

if (_rr<YFF0_Rc2)
  {  
  _r=sqrt(_rr), _u.i/=_r, _u.j/=_r, _u.k/=_r, _Q=YFF0_COULOMB_K*COULOMB_K*q_i*q_j;
  if (_r>top->yff0_sD[0]+YFF0_sD)
    {//Coulomb cutoff spline
    _x=(_r-YFF0_Rc)/YFF0_Rd;
    _d=_Q*(2.*YFF0_C3+3.*YFF0_C4*_x)*_x/YFF0_Rd;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return _Q*(YFF0_C3+YFF0_C4*_x)*_x*_x;
    }
  //Skip protons in vdW interactions
  if ( (top->ff_a[a_i].chem_id<2)||(top->ff_a[a_j].chem_id<2) ) 
    {//Do coulomb interactions only
         if (_r>top->yff0_sB[0]+YFF0_sD)
           {//Line B->C Derivative is zero here
           _d=-_Q/_rr;
           _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return _Q/_r;
           }
    else if (_r>top->yff0_sB[0]-YFF0_sD)
           {//Spline B
           _x=.5*(_r-top->yff0_sB[0]+YFF0_sD)/YFF0_sD;
           _d=_Q*(top->yff0_sQ[2]+1.5*top->yff0_sQ[3]*_x)*_x/YFF0_sD;
           _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return _Q*(top->yff0_sQ[0]+(top->yff0_sQ[2]+top->yff0_sQ[3]*_x)*_x*_x);
           }
    else   {//Line 0...A
           _d=g_i->i=g_i->j=g_i->k=g_j->i=g_j->j=g_j->k=0.;
           return _Q*top->yff0_sQ[0];
           }
    }
  else
    {//Do general stuff
         if (_r>top->yff0_sD[0]-YFF0_sD)
           {//Spline D
           _x=.5*(_r-top->yff0_sD[0]+YFF0_sD)/YFF0_sD;
           _d=(top->yff0_sD[2]/2.+top->yff0_sD[3]*_x)/YFF0_sD-_Q/_rr;
           _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return top->yff0_sD[1]+(top->yff0_sD[2]+top->yff0_sD[3]*_x)*_x+_Q/_r;
           }
    else if (_r>top->yff0_sC[0]+YFF0_sD)
           {//Line C->D
           _d=-top->yff0_sD[1]/YFF0_sD-_Q/_rr;
           _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return (top->yff0_sD[0]-_r)*top->yff0_sD[1]/YFF0_sD+_Q/_r;
           }
    else if (_r>top->yff0_sC[0]-YFF0_sD)
           {//Splane C
           _x=.5*(_r-top->yff0_sC[0]+YFF0_sD)/YFF0_sD;
           _d=top->yff0_sC[3]*_x/YFF0_sD-_Q/_rr;
           _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return top->yff0_sC[1]+top->yff0_sC[3]*_x*_x+_Q/_r;
           }
    else if (_r>top->yff0_sB[0]+YFF0_sD)
           {//Line B->C Derivative is zero here
           _d=0.-_Q/_rr;
           _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return top->yff0_sE+_Q/_r;
           }
    else if (_r>top->yff0_sB[0]-YFF0_sD)
           {//Spline B
           _x=.5*(_r-top->yff0_sB[0]+YFF0_sD)/YFF0_sD;
           _d=(.5*top->yff0_sB[2]+top->yff0_sB[3]*_x+_Q*(top->yff0_sQ[2]+1.5*top->yff0_sQ[3]*_x)*_x)/YFF0_sD;
           _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return top->yff0_sB[1]+(top->yff0_sB[2]+top->yff0_sB[3]*_x)*_x+_Q*(top->yff0_sQ[0]+(top->yff0_sQ[2]+top->yff0_sQ[3]*_x)*_x*_x);
           }
    else if (_r>top->yff0_sA[0]+YFF0_sD)
           {//Line A->B
           _d=.5*top->yff0_sB[2]/YFF0_sD;
           _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return _d*(_r-top->yff0_sA[0])+_Q*top->yff0_sQ[0];
           }
    else if (_r>top->yff0_sA[0]-YFF0_sD)
           {//Spline A
           _x=.5*(_r-top->yff0_sA[0]+YFF0_sD)/YFF0_sD;
           _d=.5*(top->yff0_sA[2]+2.*top->yff0_sA[3]*_x)/YFF0_sD;
           _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return top->yff0_sA[1]+(top->yff0_sA[2]+top->yff0_sA[3]*_x)*_x+_Q*top->yff0_sQ[0];
           }
    else   {//Line 0...A
           _d=-top->yff0_sA[1]/YFF0_sD, _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
           return top->yff0_sF+_d*_r+_Q*top->yff0_sQ[0];
           }
    }
  }
else return 0.;
}
inline double calc_atom__nb_enrg_yff0(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top)
{
register double _rr, _r, _x, _Q;

_rr=calc_distance(r_i,r_j)+SMALL2; //To keep divisibility

if (_rr<YFF0_Rc2)
  {  
  _r=sqrt(_rr), _Q=YFF0_COULOMB_K*COULOMB_K*q_i*q_j;
  if (_r>top->yff0_sD[0]+YFF0_sD)
    {//Coulomb cutoff spline
    _x=(_r-YFF0_Rc)/YFF0_Rd;
    return _Q*(YFF0_C3+YFF0_C4*_x)*_x*_x;
    }
  //Skip protons in vdW interactions
  if ( (top->ff_a[a_i].chem_id<2)||(top->ff_a[a_j].chem_id<2) ) 
    {//Do coulomb interactions only
         if (_r>top->yff0_sB[0]+YFF0_sD)
           {//Line B->C Derivative is zero here
           return _Q/_r;
           }
    else if (_r>top->yff0_sB[0]-YFF0_sD)
           {//Spline B
           _x=.5*(_r-top->yff0_sB[0]+YFF0_sD)/YFF0_sD;
           return _Q*(top->yff0_sQ[0]+(top->yff0_sQ[2]+top->yff0_sQ[3]*_x)*_x*_x);
           }
    else   {//Line 0...A
           return _Q*top->yff0_sQ[0];
           }
    }
  else
    {//Do general stuff
         if (_r>top->yff0_sD[0]-YFF0_sD)
           {//Spline D
           _x=.5*(_r-top->yff0_sD[0]+YFF0_sD)/YFF0_sD;
           return top->yff0_sD[1]+(top->yff0_sD[2]+top->yff0_sD[3]*_x)*_x+_Q/_r;
           }
    else if (_r>top->yff0_sC[0]+YFF0_sD)
           {//Line C->D
           return (top->yff0_sD[0]-_r)*top->yff0_sD[1]/YFF0_sD+_Q/_r;
           }
    else if (_r>top->yff0_sC[0]-YFF0_sD)
           {//Splane C
           _x=.5*(_r-top->yff0_sC[0]+YFF0_sD)/YFF0_sD;
           return top->yff0_sC[1]+top->yff0_sC[3]*_x*_x+_Q/_r;
           }
    else if (_r>top->yff0_sB[0]+YFF0_sD)
           {//Line B->C Derivative is zero here
           return top->yff0_sE+_Q/_r;
           }
    else if (_r>top->yff0_sB[0]-YFF0_sD)
           {//Spline B
           _x=.5*(_r-top->yff0_sB[0]+YFF0_sD)/YFF0_sD;
           return top->yff0_sB[1]+(top->yff0_sB[2]+top->yff0_sB[3]*_x)*_x+_Q*(top->yff0_sQ[0]+(top->yff0_sQ[2]+top->yff0_sQ[3]*_x)*_x*_x);
           }
    else if (_r>top->yff0_sA[0]+YFF0_sD)
           {//Line A->B
           return .5*top->yff0_sB[2]*(_r-top->yff0_sA[0])/YFF0_sD+_Q*top->yff0_sQ[0];
           }
    else if (_r>top->yff0_sA[0]-YFF0_sD)
           {//Spline A
           _x=.5*(_r-top->yff0_sA[0]+YFF0_sD)/YFF0_sD;
           return top->yff0_sA[1]+(top->yff0_sA[2]+top->yff0_sA[3]*_x)*_x+_Q*top->yff0_sQ[0];
           }
    else   {//Line 0...A
           return top->yff0_sF-top->yff0_sA[1]*_r/YFF0_sD+_Q*top->yff0_sQ[0];
           }
    }
  }
else return 0.;
}

//This function calculates nonbonded vdW energy in yff0. Do not call it for hydrogens!
inline double calc_atom__nb_vdW_enrg_yff0(register double _r,t_top *top)
{
register double _x;

if (_r<YFF0_Rb)
  {//Do general stuff  
       if (_r>top->yff0_sD[0]-YFF0_sD)
         {//Spline D
         _x=.5*(_r-top->yff0_sD[0]+YFF0_sD)/YFF0_sD;
         return top->yff0_sD[1]+(top->yff0_sD[2]+top->yff0_sD[3]*_x)*_x;
         }
  else if (_r>top->yff0_sC[0]+YFF0_sD)
         {//Line C->D
         return (top->yff0_sD[0]-_r)*top->yff0_sD[1]/YFF0_sD;
         }
  else if (_r>top->yff0_sC[0]-YFF0_sD)
         {//Splane C
         _x=.5*(_r-top->yff0_sC[0]+YFF0_sD)/YFF0_sD;
         return top->yff0_sC[1]+top->yff0_sC[3]*_x*_x;
         }
  else if (_r>top->yff0_sB[0]+YFF0_sD)
         {//Line B->C Derivative is zero here
         return top->yff0_sE;
         }
  else if (_r>top->yff0_sB[0]-YFF0_sD)
         {//Spline B
         _x=.5*(_r-top->yff0_sB[0]+YFF0_sD)/YFF0_sD;
         return top->yff0_sB[1]+(top->yff0_sB[2]+top->yff0_sB[3]*_x)*_x;
         }
  else if (_r>top->yff0_sA[0]+YFF0_sD)
         {//Line A->B
         return .5*top->yff0_sB[2]*(_r-top->yff0_sA[0])/YFF0_sD;
         }
  else if (_r>top->yff0_sA[0]-YFF0_sD)
         {//Spline A
         _x=.5*(_r-top->yff0_sA[0]+YFF0_sD)/YFF0_sD;
         return top->yff0_sA[1]+(top->yff0_sA[2]+top->yff0_sA[3]*_x)*_x;
         }
  else   {//Line 0...A
         return top->yff0_sF-top->yff0_sA[1]*_r/YFF0_sD;
         }
  }
return 0.;
}
//This function calculates nonbonded coulomb energy in yff0
inline double calc_atom__nb_coulomb_enrg_yff0(register double _Q,register double _r,t_top *top)
{
register double _x;

if (_r<YFF0_Rc)
  {//Do general stuff
       if (_r>top->yff0_sD[0]+YFF0_sD)
         {//Coulomb cutoff spline
         _x=(_r-YFF0_Rc)/YFF0_Rd;
         return _Q*(YFF0_C3+YFF0_C4*_x)*_x*_x;
         }
  else if (_r>top->yff0_sB[0]+YFF0_sD)
         {//curve B->Rb
         return _Q/_r;
         }
  else if (_r>top->yff0_sB[0]-YFF0_sD)
         {//Spline B
         _x=.5*(_r-top->yff0_sB[0]+YFF0_sD)/YFF0_sD;
         return _Q*(top->yff0_sQ[0]+(top->yff0_sQ[2]+top->yff0_sQ[3]*_x)*_x*_x);
         }
  else   {//Line 0...B
         return _Q*top->yff0_sQ[0];
         }
  }
return 0.;
}

//------------------------------------  E N E R G Y   F U N C T I O N S   Y F F 0  ------------------------------------------

//This function calculates energy in the active layer and marks atoms to exclude from grid summing
inline void calc_nbenergy_yff0(double *e,t_ffsys *ffsys,t_top *top)
{
unsigned int _i, _j, _k, _l, mol_i, mol_j, active_i, active_j, a_i, a_j;
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
              if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                *e+=_e;
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
              if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                *e+=_e;
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
              if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                *e+=_e;
                }
              }
            }
          }
        }
      }
    //Do inside mol: remove anchor-[anchor == root]-anchor and remove pairs
    if ((active_i=ffsys->mols[mol_i]->anchors->size)!=1)
      while (active_i--)
        {
        //Do root-rest
        _j=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
        if (active_i!=ffsys->rtrees[mol_i]->root) a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0];
        else a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch].edge.vertice[1];                                     
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            *e-=_e;
            }
          }
        //Do rest-rest
        _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
        while (--_i)
          {
          a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]].edge.vertice[1];
          _j=_i;
          while (--_j)  
            {
            a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
            if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
              {
              *e-=_e;
              }
            }
          }
        }
    //Edit 1-4: No such interactions in yff0
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
            if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
              {
              *e+=_e;
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
              if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                *e+=_e;
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
                     if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                       {
                       *e+=_e;
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
                     if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                       {
                       *e+=_e;
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
                   if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                     {
                     *e+=_e;
                     } 
                   }
                 }
               }
        }
      }
    //Do inside mol: remove anchor-[anchor == root]-anchor and remove pairs
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      //Do root-rest
      _j=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
      a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0];
      if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch>0)
        {
        _l=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (--_l)
          {
          a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_l];
          if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            *e-=_e;
            }
          }
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            *e-=_e;
            }
          }
        }
      else
        while (--_j)
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];  
          if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            *e-=_e;
            }
          }
      //Do rest-rest
      _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
      while (--_i)
        {
        a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]].edge.vertice[1];
        _j=_i;
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            *e-=_e;
            }
          }
        }
      }
    //Edit 1-4 : NO SUCH INTERACTIONS IN YFF0
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
            if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
              {
              *e+=_e;
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
              if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                *e+=_e;
                }   
              } 
            }
          }
        }  
      }
    }
}
//Note it scales repulsion as E*=E*(3*x^2-2*x^3) in range 0 @ px==0. ... E @ px=1. dE*/dr=(3*x^2-2*x^3)*dE/dr
//This function calculates energy in the active layer and marks atoms to exclude from grid summing
inline void calc_nbgrad_yff0(double *e,t_ffsys *ffsys,t_top *top)
{
unsigned int _i, _j, _k, _l, mol_i, mol_j, active_i, active_j, a_i, a_j;
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
              if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
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
              if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
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
              if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
                }
              }
            }
          }
        }
      }
    //Do inside mol: remove anchor-[anchor == root]-anchor and remove pairs
    if ((active_i=ffsys->mols[mol_i]->anchors->size)!=1)
      while (active_i--)
        {
        //Do root-rest
        _j=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
        if (active_i!=ffsys->rtrees[mol_i]->root) a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0];
        else a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch].edge.vertice[1];                                     
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            *e-=_e, ffsys->g[a_i].i-=g_i.i, ffsys->g[a_i].j-=g_i.j, ffsys->g[a_i].k-=g_i.k, ffsys->g[a_j].i-=g_j.i, ffsys->g[a_j].j-=g_j.j, ffsys->g[a_j].k-=g_j.k;
            }
          }
        //Do rest-rest
        _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
        while (--_i)
          {
          a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]].edge.vertice[1];
          _j=_i;
          while (--_j)  
            {
            a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
            if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
              {
              *e-=_e, ffsys->g[a_i].i-=g_i.i, ffsys->g[a_i].j-=g_i.j, ffsys->g[a_i].k-=g_i.k, ffsys->g[a_j].i-=g_j.i, ffsys->g[a_j].j-=g_j.j, ffsys->g[a_j].k-=g_j.k;
              }
            }
          }
        }
    //Edit 1-4, No such interactions in yff0
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
            if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
              {
              *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
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
              if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
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
                     if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                       {
                       *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
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
                     if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                       {
                       *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
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
                   if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                     {
                     *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
                     } 
                   }
                 }
               }
        }
      }
    //Do inside mol: remove anchor-[anchor == root]-anchor and remove pairs
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      //Do root-rest
      _j=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
      a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0];
      if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch>0)
        {
        _l=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (--_l)
          {
          a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_l];
          if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            *e-=_e, ffsys->g[a_i].i-=g_i.i, ffsys->g[a_i].j-=g_i.j, ffsys->g[a_i].k-=g_i.k, ffsys->g[a_j].i-=g_j.i, ffsys->g[a_j].j-=g_j.j, ffsys->g[a_j].k-=g_j.k;
            }
          }
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            *e-=_e, ffsys->g[a_i].i-=g_i.i, ffsys->g[a_i].j-=g_i.j, ffsys->g[a_i].k-=g_i.k, ffsys->g[a_j].i-=g_j.i, ffsys->g[a_j].j-=g_j.j, ffsys->g[a_j].k-=g_j.k;
            }
          }
        }
      else
        while (--_j)
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];  
          if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            *e-=_e, ffsys->g[a_i].i-=g_i.i, ffsys->g[a_i].j-=g_i.j, ffsys->g[a_i].k-=g_i.k, ffsys->g[a_j].i-=g_j.i, ffsys->g[a_j].j-=g_j.j, ffsys->g[a_j].k-=g_j.k;
            }
          }
      //Do rest-rest
      _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
      while (--_i)
        {
        a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]].edge.vertice[1];
        _j=_i;
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            *e-=_e, ffsys->g[a_i].i-=g_i.i, ffsys->g[a_i].j-=g_i.j, ffsys->g[a_i].k-=g_i.k, ffsys->g[a_j].i-=g_j.i, ffsys->g[a_j].j-=g_j.j, ffsys->g[a_j].k-=g_j.k;
            }
          }
        }
      }
    //Edit 1-4. No such interactions in yff0
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
            if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
              {
              *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
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
              if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k, ffsys->g[a_j].i+=g_j.i, ffsys->g[a_j].j+=g_j.j, ffsys->g[a_j].k+=g_j.k;
                }   
              } 
            }
          }
        }  
      }
    }
}


//!!!!!!!!!!!!! This function can be improved by excluding the exact anchor in explicite summation of ongrid failure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculates yff0 energy on dgrid
//It uses marks made by calc_yff0_grad(t_ffsys *ffsys,t_top *top) to skip frozen atoms
inline void calc_energy_on_grid_yff0(double *e,t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***Q)
{
unsigned int _i, _j, _k, _l, mol_i, mol_j, active_i, active_j, a_i, a_j;
double _e, a, q;

mol_i=ffsys->nmols;
while (mol_i--)
  if (!(ffsys->activem[mol_i].active_a))
    {//Score all
    _i=ffsys->mols[mol_i]->atoms->size;
    while (_i--)
      {
      a_i=ffsys->nr[mol_i]+_i, a=0.;
      if (                                                       ( (calc_tricubic_interpolation_wp_monotonicity(&q,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) &&
           ( (top->ff_a[(unsigned int)ffsys->a[a_i]].chem_id<2)||( (calc_tricubic_interpolation_wp_monotonicity(&a,ori,sp,ni,nj,nk,&ffsys->r[a_i],A))) ) )
        {
        //Summ on grid
        if ( (_e=a+ffsys->q[a_i]*q))
          { 
          *e+=_e;
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
                if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                  {
                  *e-=_e;
                  }
                }
              }
            //Summ all 
            _j=ffsys->mols[mol_j]->atoms->size;
            while (_j--)
              {
              a_j=ffsys->nr[mol_j]+_j;
              if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                *e+=_e;
                }
              }
            }
        }
      }
    }
  else
    {
    //Special case - check for mark from calc_grad_yff0(t_ffsys *ffsys,t_top *top) and skip if required (it removes marks)
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch<0) //skip frozen atoms otherwise
        {
        *ffsys->rtrees[mol_i]->rbranch[active_i].rbranch=-(int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch;
        a_i=ffsys->nr[mol_i]+*ffsys->mols[mol_i]->anchors->list[active_i].list, a=0.;
        if (                                                       ( (calc_tricubic_interpolation_wp_monotonicity(&q,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) &&
             ( (top->ff_a[(unsigned int)ffsys->a[a_i]].chem_id<2)||( (calc_tricubic_interpolation_wp_monotonicity(&a,ori,sp,ni,nj,nk,&ffsys->r[a_i],A))) ) )
          {
          //Summ on grid
          if ( (_e=a+ffsys->q[a_i]*q))
            {
            *e+=_e;
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
                  if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                    {
                    *e-=_e;
                    }
                  }
                }
              //Summ all
              _j=ffsys->mols[mol_j]->atoms->size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+_j;
                if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                  {
                  *e+=_e;
                  }
                }
              }
          }
        }
     //Summ the rest 
     _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
      while (--_i)
        {
        a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i], a=0.;
        if (                                                       ( (calc_tricubic_interpolation_wp_monotonicity(&q,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) &&
             ( (top->ff_a[(unsigned int)ffsys->a[a_i]].chem_id<2)||( (calc_tricubic_interpolation_wp_monotonicity(&a,ori,sp,ni,nj,nk,&ffsys->r[a_i],A))) ) )
          {
          //Summ on grid
          if ( (_e=a+ffsys->q[a_i]*q))
            {
            *e+=_e;
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
                  if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                    {
                    *e-=_e;
                    }
                  }
                }
              //Summ all 
              _j=ffsys->mols[mol_j]->atoms->size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+_j;
                if ( (_e=calc_atom__nb_enrg_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                  {
                  *e+=_e;
                  }
                }
              }
          }
        }
      }
    }
}
//!!!!!!!!!!!!! This function can be improved by excluding the exact anchor in explicite summation of ongrid failure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculates yff0 energy on dgrid
//It uses marks made by calc_yff0_grad(t_ffsys *ffsys,t_top *top) to skip frozen atoms
inline void calc_grad_on_grid_yff0(double *e,t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***Q)
{
unsigned int _i, _j, _k, _l, mol_i, mol_j, active_i, active_j, a_i, a_j;
double _e, a, q;
t_vec g_i, g_j, da, dq;

mol_i=ffsys->nmols;
while (mol_i--)
  if (!(ffsys->activem[mol_i].active_a))
    {
    _i=ffsys->mols[mol_i]->atoms->size;
    while (_i--)
      {
      a_i=ffsys->nr[mol_i]+_i, a=da.i=da.j=da.k=0.;
      if (                                                       ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) &&
           ( (top->ff_a[(unsigned int)ffsys->a[a_i]].chem_id<2)||( (calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,ni,nj,nk,&ffsys->r[a_i],A))) ) )
        {
        //Summ on grid
        if ( (_e=a+ffsys->q[a_i]*q))
          { 
          multiple_vec_scalar(&dq,&dq,ffsys->q[a_i]);
          *e+=_e, ffsys->g[a_i].i+=da.i+dq.i, ffsys->g[a_i].j+=da.j+dq.j, ffsys->g[a_i].k+=da.k+dq.k;
          }
        }
      else
        {
        //Cross summ
        for (mol_j=0; mol_j<ffsys->nmols; mol_j++)
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
                if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  *e-=_e, ffsys->g[a_i].i-=g_i.i, ffsys->g[a_i].j-=g_i.j, ffsys->g[a_i].k-=g_i.k;
                  }
                }
              }
            //Summ all 
            _j=ffsys->mols[mol_j]->atoms->size;
            while (_j--)
              {
              a_j=ffsys->nr[mol_j]+_j;
              if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                {
                *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k;
                }
              }
            }
        }
      }
    }
  else
    {
    //Special case - check for mark from calc_yff0_grad(t_ffsys *ffsys,t_top *top) and skip if required (it removes marks)
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch<0) //skip frozen atoms otherwise
        {
        *ffsys->rtrees[mol_i]->rbranch[active_i].rbranch=-(int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch;
        a_i=ffsys->nr[mol_i]+*ffsys->mols[mol_i]->anchors->list[active_i].list, a=da.i=da.j=da.k=0.;
        if (                                                       ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) &&
             ( (top->ff_a[(unsigned int)ffsys->a[a_i]].chem_id<2)||( (calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,ni,nj,nk,&ffsys->r[a_i],A))) ) )
          {
          //Summ on grid
          if ( (_e=a+ffsys->q[a_i]*q))
            {
            multiple_vec_scalar(&dq,&dq,ffsys->q[a_i]);
            *e+=_e, ffsys->g[a_i].i+=da.i+dq.i, ffsys->g[a_i].j+=da.j+dq.j, ffsys->g[a_i].k+=da.k+dq.k;
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
                  if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    *e-=_e, ffsys->g[a_i].i-=g_i.i, ffsys->g[a_i].j-=g_i.j, ffsys->g[a_i].k-=g_i.k;
                    }
                  }
                }
              //Summ all 
              _j=ffsys->mols[mol_j]->atoms->size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+_j;
                if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k;
                  }
                }
              }
          }
        }
     //Summ the rest 
     _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
      while (--_i)
        {
        a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i], a=da.i=da.j=da.k=0.;
        if (                                                       ( (calc_tricubic_interpolation_derivative_wp_monotonicity(&q,&dq,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) &&
             ( (top->ff_a[(unsigned int)ffsys->a[a_i]].chem_id<2)||( (calc_tricubic_interpolation_derivative_wp_monotonicity(&a,&da,ori,sp,ni,nj,nk,&ffsys->r[a_i],A))) ) )
          {
          //Summ on grid
          if ( (_e=a+ffsys->q[a_i]*q))
            {
            multiple_vec_scalar(&dq,&dq,ffsys->q[a_i]);
            *e+=_e, ffsys->g[a_i].i+=da.i+dq.i, ffsys->g[a_i].j+=da.j+dq.j, ffsys->g[a_i].k+=da.k+dq.k;
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
                  if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                    {
                    *e-=_e, ffsys->g[a_i].i-=g_i.i, ffsys->g[a_i].j-=g_i.j, ffsys->g[a_i].k-=g_i.k;
                    }
                  }
                }
              //Summ all 
              _j=ffsys->mols[mol_j]->atoms->size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+_j;
                if ( (_e=calc_atom__nb_yff0(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
                  {
                  *e+=_e, ffsys->g[a_i].i+=g_i.i, ffsys->g[a_i].j+=g_i.j, ffsys->g[a_i].k+=g_i.k;
                  }
                }
              }
          }
        }
      }
    }
}

//This function calculates torsions energy using gaff primitives
inline void calc_dihs_energy_gaff(double *e,t_ffsys *ffsys)
{
register unsigned int mol_id;
mol_id=ffsys->nmols; while (mol_id--) *e+=calc_mol_tbenergy_gaff(&ffsys->r[ffsys->nr[mol_id]],ffsys->mols[mol_id],&ffsys->activem[mol_id]);
}
//This function calculates torsions energy and its derivative using gaff primitives
inline void calc_dihs_grad_gaff(double *e,t_ffsys *ffsys)
{
register unsigned int mol_id;
mol_id=ffsys->nmols; while (mol_id--) *e+=calc_mol_tbgrad_gaff(&ffsys->g[ffsys->nr[mol_id]],&ffsys->r[ffsys->nr[mol_id]],ffsys->mols[mol_id],&ffsys->activem[mol_id]);
}

//This function calculates yff0 energy on grid
inline double calc_ffsys_nbenergy_on_dgrid_yff0(t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *Q)
{
//Init energies and forces
double  e=0.;
//Summ active layers
calc_nbenergy_yff0(&e,ffsys,top);
calc_dihs_energy_gaff(&e,ffsys);
calc_energy_on_grid_yff0(&e,ffsys,top,&A->ori,A->len.i,A->len.j,A->len.k,A->sp,A->d,Q->d);

return e;
}
//This function calculates yff0 energy gradient on grid
inline double calc_ffsys_nbgrad_on_dgrid_yff0(t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *Q)
{
//Init energies and forces
double  e=0.;
memset(ffsys->g,0x0,ffsys->natoms*sizeof(t_vec));
//Summ active layers
calc_nbgrad_yff0(&e,ffsys,top);
calc_dihs_grad_gaff(&e,ffsys);
calc_grad_on_grid_yff0(&e,ffsys,top,&A->ori,A->len.i,A->len.j,A->len.k,A->sp,A->d,Q->d);

return e;
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
inline double _calc_gbond_gaff(double r,double *de,t_ff_b *ff_b)
{
r-=ff_b->v;
*de=2.*ff_b->k*r;
return ff_b->k*r*r;
}
//This function calculates energies and forces of angle
inline double _calc_gangle_gaff(double csA,double *de,t_ff_g *ff_g)
{
double alpha;
alpha=acos(csA)*180./PI-ff_g->v;
*de=2.*ff_g->k*alpha;
return ff_g->k*alpha*alpha;
}
//This function calculates energies and forces of improper angle
inline double _calc_gimpr_gaff(double csF,double snF,double *de,t_ff_i *ff_i)
{
double phi;
//Switch the most popular cases
     if (ff_i->n==1.)
       {//(x-y)
       *de=ff_i->k*ff_i->n*(snF*ff_i->csF0-csF*ff_i->snF0);
       return ff_i->k*(1.+(csF*ff_i->csF0+snF*ff_i->snF0));
       }
else if (ff_i->n==2.) 
       {//(2x-y)
       *de=ff_i->k*ff_i->n*(2.*snF*csF*ff_i->csF0-(2.*csF*csF-1.)*ff_i->snF0);
       return ff_i->k*(1.+((2.*csF*csF-1.)*ff_i->csF0+2.*snF*csF*ff_i->snF0));
       }
else if (ff_i->n==3.)
       {//(3x-y)
       *de=ff_i->k*ff_i->n*((3.-4.*snF)*snF*ff_i->csF0-(4.*csF-3.)*csF*ff_i->snF0);
       return ff_i->k*(1.+((4.*csF-3.)*csF*ff_i->csF0+(3.-4.*snF)*snF*ff_i->snF0));
       }
else   {//An arbitrary case
       if (fabs(csF)<SMALL2)
         {
         if (snF>0) phi=+ff_i->n*PI/2.-ff_i->v;
         else       phi=-ff_i->n*PI/2.-ff_i->v;
         }
       else phi=ff_i->n*atan(snF/csF)-ff_i->v;
       *de=-ff_i->k*ff_i->n*sin(phi);
       return ff_i->k*(1.+cos(phi));
       }
}
//This function calculates energies and forces of torsion
inline double _calc_gdihs_gaff(double csF,double snF,double *de,t_ff_d *ff_d)
{
double phi;
//Switch the most popular cases
     if (ff_d->n==1.)
       {//(x-y)
       *de=ff_d->k*ff_d->n*(snF*ff_d->csF0-csF*ff_d->snF0);
       return ff_d->k*(1.+(csF*ff_d->csF0+snF*ff_d->snF0));
       }
else if (ff_d->n==2.) 
       {//(2x-y)
       *de=ff_d->k*ff_d->n*(2.*snF*csF*ff_d->csF0-(2.*csF*csF-1.)*ff_d->snF0);
       return ff_d->k*(1.+((2.*csF*csF-1.)*ff_d->csF0+2.*snF*csF*ff_d->snF0));
       }
else if (ff_d->n==3.)
       {//(3x-y)
       *de=ff_d->k*ff_d->n*((3.-4.*snF)*snF*ff_d->csF0-(4.*csF-3.)*csF*ff_d->snF0);
       return ff_d->k*(1.+((4.*csF-3.)*csF*ff_d->csF0+(3.-4.*snF)*snF*ff_d->snF0));
       }
else   {//An arbitrary case
       if (fabs(csF)<SMALL2)
         {
         if (snF>0) phi=+ff_d->n*PI/2.-ff_d->v;
         else       phi=-ff_d->n*PI/2.-ff_d->v;
         }
       else phi=ff_d->n*atan(snF/csF)-ff_d->v;
       *de=-ff_d->k*ff_d->n*sin(phi);
       return ff_d->k*(1.+cos(phi));
       }
}
inline double _calc_gdihs_enrg_gaff(double csF,double snF,t_ff_d *ff_d)
{
double phi;
//Switch the most popular cases
     if (ff_d->n==1.)
       {//(x-y)
       return ff_d->k*(1.+(csF*ff_d->csF0+snF*ff_d->snF0));
       }
else if (ff_d->n==2.) 
       {//(2x-y)
       return ff_d->k*(1.+((2.*csF*csF-1.)*ff_d->csF0+2.*snF*csF*ff_d->snF0));
       }
else if (ff_d->n==3.)
       {//(3x-y)
       return ff_d->k*(1.+((4.*csF-3.)*csF*ff_d->csF0+(3.-4.*snF)*snF*ff_d->snF0));
       }
else   {//An arbitrary case
       if (fabs(csF)<SMALL2)
         {
         if (snF>0) phi=+ff_d->n*PI/2.-ff_d->v;
         else       phi=-ff_d->n*PI/2.-ff_d->v;
         }
       else phi=ff_d->n*atan(snF/csF)-ff_d->v;
       return ff_d->k*(1.+cos(phi));
       }
}

//This function calculated bonded interactions for whole molecule
inline double calc__bgrad_gaff(t_vec *g,t_vec *r,t_mol *mol,t_activem *activem)
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
    e+=_calc_gbond_gaff(csA,&de,&mol->ff_b[_i]);
    g[mol->ff_b[_i].atom[0]].i+=de*_g[0], g[mol->ff_b[_i].atom[0]].j+=de*_g[1], g[mol->ff_b[_i].atom[0]].k+=de*_g[2];
    g[mol->ff_b[_i].atom[1]].i+=de*_g[3], g[mol->ff_b[_i].atom[1]].j+=de*_g[4], g[mol->ff_b[_i].atom[1]].k+=de*_g[5];
    }
  _i=mol->size_g;
  while (_i--)
    {
    calc_angle_derivative(&csA,&snA,&r[mol->ff_g[_i].atom[0]],&r[mol->ff_g[_i].atom[1]],&r[mol->ff_g[_i].atom[2]],_g);
    e+=_calc_gangle_gaff(csA,&de,&mol->ff_g[_i]);
    g[mol->ff_g[_i].atom[0]].i+=de*_g[0], g[mol->ff_g[_i].atom[0]].j+=de*_g[1], g[mol->ff_g[_i].atom[0]].k+=de*_g[2];
    g[mol->ff_g[_i].atom[1]].i+=de*_g[3], g[mol->ff_g[_i].atom[1]].j+=de*_g[4], g[mol->ff_g[_i].atom[1]].k+=de*_g[5];
    g[mol->ff_g[_i].atom[2]].i+=de*_g[6], g[mol->ff_g[_i].atom[2]].j+=de*_g[7], g[mol->ff_g[_i].atom[2]].k+=de*_g[8];
    }
  _i=mol->size_i;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_i[_i].atom[0]],&r[mol->ff_i[_i].atom[1]],&r[mol->ff_i[_i].atom[2]],&r[mol->ff_i[_i].atom[3]],_g);
    e+=_calc_gimpr_gaff(csA,snA,&de,&mol->ff_i[_i]);
    g[mol->ff_i[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_i[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_i[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_i[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_i[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_i[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_i[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_i[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_i[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_i[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_i[_i].atom[3]].j+=de*_g[10], g[mol->ff_i[_i].atom[3]].k+=de*_g[11];
    }
  _i=mol->size_d;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]],_g);
    e+=_calc_gdihs_gaff(csA,snA,&de,&mol->ff_d[_i]);
    g[mol->ff_d[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_d[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_d[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_d[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_d[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_d[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_d[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_d[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_d[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_d[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_d[_i].atom[3]].j+=de*_g[10], g[mol->ff_d[_i].atom[3]].k+=de*_g[11];
    }
  }
else
  {
  _j=activem->active_b->size;
  while (_j--)
    {
    _i=activem->active_b->list[_j];
    calc_bond_derivative(&csA,&r[mol->ff_b[_i].atom[0]],&r[mol->ff_b[_i].atom[1]],_g);
    g[mol->ff_b[_i].atom[0]].i+=de*_g[0], g[mol->ff_b[_i].atom[0]].j+=de*_g[1], g[mol->ff_b[_i].atom[0]].k+=de*_g[2];
    g[mol->ff_b[_i].atom[1]].i+=de*_g[3], g[mol->ff_b[_i].atom[1]].j+=de*_g[4], g[mol->ff_b[_i].atom[1]].k+=de*_g[5];
    }
  _j=activem->active_g->size;
  while (_j--)
    {
    _i=activem->active_g->list[_j];
    calc_angle_derivative(&csA,&snA,&r[mol->ff_g[_i].atom[0]],&r[mol->ff_g[_i].atom[1]],&r[mol->ff_g[_i].atom[2]],_g);
    e+=_calc_gangle_gaff(csA,&de,&mol->ff_g[_i]);
    g[mol->ff_g[_i].atom[0]].i+=de*_g[0], g[mol->ff_g[_i].atom[0]].j+=de*_g[1], g[mol->ff_g[_i].atom[0]].k+=de*_g[2];
    g[mol->ff_g[_i].atom[1]].i+=de*_g[3], g[mol->ff_g[_i].atom[1]].j+=de*_g[4], g[mol->ff_g[_i].atom[1]].k+=de*_g[5];
    g[mol->ff_g[_i].atom[2]].i+=de*_g[6], g[mol->ff_g[_i].atom[2]].j+=de*_g[7], g[mol->ff_g[_i].atom[2]].k+=de*_g[8];
    }
  _j=activem->active_i->size;
  while (_j--)
    {
    _i=activem->active_i->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_i[_i].atom[0]],&r[mol->ff_i[_i].atom[1]],&r[mol->ff_i[_i].atom[2]],&r[mol->ff_i[_i].atom[3]],_g);
    e+=_calc_gimpr_gaff(csA,snA,&de,&mol->ff_i[_i]);
    g[mol->ff_i[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_i[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_i[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_i[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_i[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_i[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_i[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_i[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_i[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_i[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_i[_i].atom[3]].j+=de*_g[10], g[mol->ff_i[_i].atom[3]].k+=de*_g[11];
    }
  _j=activem->active_d->size;
  while (_j--)
    {
    _i=activem->active_d->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]],_g);
    e+=_calc_gdihs_gaff(csA,snA,&de,&mol->ff_d[_i]);
    g[mol->ff_d[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_d[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_d[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_d[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_d[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_d[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_d[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_d[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_d[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_d[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_d[_i].atom[3]].j+=de*_g[10], g[mol->ff_d[_i].atom[3]].k+=de*_g[11];
    }
  }
return e;
}
inline void _calc__bgrad_gaff(double e[N_ESLICES],t_vec (*g)[N_ESLICES],t_vec *r,t_mol *mol,t_activem *activem)
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
    _e=_calc_gbond_gaff(csA,&de,&mol->ff_b[_i]);
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
    _e=_calc_gangle_gaff(csA,&de,&mol->ff_g[_i]);
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
    _e=_calc_gimpr_gaff(csA,snA,&de,&mol->ff_i[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_i[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_i[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_i[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_i[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_i[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_i[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_i[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_i[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_i[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_i[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_i[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_i[_i].atom[3]][_L].k+=de*_g[11];
    }
  _i=mol->size_d;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]],_g);
    _e=_calc_gdihs_gaff(csA,snA,&de,&mol->ff_d[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_d[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_d[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_d[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_d[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_d[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_d[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_d[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_d[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_d[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_d[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_d[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_d[_i].atom[3]][_L].k+=de*_g[11];
    }
  }
else
  {
  _j=activem->active_b->size;
  while (_j--)
    {
    _i=activem->active_b->list[_j];
    calc_bond_derivative(&csA,&r[mol->ff_b[_i].atom[0]],&r[mol->ff_b[_i].atom[1]],_g);
    _e=_calc_gbond_gaff(csA,&de,&mol->ff_b[_i]);
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
    _e=_calc_gangle_gaff(csA,&de,&mol->ff_g[_i]);
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
    _e=_calc_gimpr_gaff(csA,snA,&de,&mol->ff_i[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_i[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_i[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_i[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_i[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_i[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_i[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_i[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_i[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_i[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_i[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_i[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_i[_i].atom[3]][_L].k+=de*_g[11];
    }
  _j=activem->active_d->size;
  while (_j--)
    {
    _i=activem->active_d->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]],_g);
    _e=_calc_gdihs_gaff(csA,snA,&de,&mol->ff_d[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_d[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_d[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_d[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_d[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_d[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_d[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_d[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_d[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_d[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_d[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_d[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_d[_i].atom[3]][_L].k+=de*_g[11];
    }
  }
}

//These functions calculate torsional energy
inline double calc_mol_tbenergy_gaff(t_vec *r,t_mol *mol,t_activem *activem)
{
unsigned int _i,_j;
double e=0., csA, snA;

if (!(activem->active_a))
  {
  _i=mol->size_d;
  while (_i--)
    {
    calc_dih_angle(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]]);
    e+=_calc_gdihs_enrg_gaff(csA,snA,&mol->ff_d[_i]);
    }
  }
else
  {
  _j=activem->active_d->size;
  while (_j--)
    {
    _i=activem->active_d->list[_j];
    calc_dih_angle(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]]);
    e+=_calc_gdihs_enrg_gaff(csA,snA,&mol->ff_d[_i]);
    }
  }
return e;
}
inline void _calc_mol_tbenergy_gaff(double e[N_ESLICES],t_vec *r,t_mol *mol,t_activem *activem)
{
unsigned int _i, _j, _L;
double _e, csA, snA;

if (!(activem->active_a))
  {
  _i=mol->size_d;
  while (_i--)
    {
    calc_dih_angle(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]]);
    _e=_calc_gdihs_enrg_gaff(csA,snA,&mol->ff_d[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e;
    }
  }
else
  {
  _j=activem->active_d->size;
  while (_j--)
    {
    _i=activem->active_d->list[_j];
    calc_dih_angle(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]]);
    _e=_calc_gdihs_enrg_gaff(csA,snA,&mol->ff_d[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e;
    }
  }
}
//These functions calculate torsional energy for internal coordinates gradients
inline double calc_mol_tbgrad_gaff(t_vec *g,t_vec *r,t_mol *mol,t_activem *activem)
{
unsigned int _i,_j;
double e=0., de, csA, snA, _g[12];

if (!(activem->active_a))
  {
  _i=mol->size_d;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]],_g);
    e+=_calc_gdihs_gaff(csA,snA,&de,&mol->ff_d[_i]);
    g[mol->ff_d[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_d[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_d[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_d[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_d[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_d[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_d[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_d[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_d[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_d[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_d[_i].atom[3]].j+=de*_g[10], g[mol->ff_d[_i].atom[3]].k+=de*_g[11];
    }
  }
else
  {
  _j=activem->active_d->size;
  while (_j--)
    {
    _i=activem->active_d->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]],_g);
    e+=_calc_gdihs_gaff(csA,snA,&de,&mol->ff_d[_i]);
    g[mol->ff_d[_i].atom[0]].i+=de*_g[ 0], g[mol->ff_d[_i].atom[0]].j+=de*_g[ 1], g[mol->ff_d[_i].atom[0]].k+=de*_g[ 2];
    g[mol->ff_d[_i].atom[1]].i+=de*_g[ 3], g[mol->ff_d[_i].atom[1]].j+=de*_g[ 4], g[mol->ff_d[_i].atom[1]].k+=de*_g[ 5];
    g[mol->ff_d[_i].atom[2]].i+=de*_g[ 6], g[mol->ff_d[_i].atom[2]].j+=de*_g[ 7], g[mol->ff_d[_i].atom[2]].k+=de*_g[ 8];
    g[mol->ff_d[_i].atom[3]].i+=de*_g[ 9], g[mol->ff_d[_i].atom[3]].j+=de*_g[10], g[mol->ff_d[_i].atom[3]].k+=de*_g[11];
    }
  }
return e;
}
inline void _calc_mol_tbgrad_gaff(double e[N_ESLICES],t_vec (*g)[N_ESLICES],t_vec *r,t_mol *mol,t_activem *activem)
{
unsigned int _i, _j, _L;
double _e, de, csA, snA, _g[12];

if (!(activem->active_a))
  {
  _i=mol->size_d;
  while (_i--)
    {
    calc_dih_derivative(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]],_g);
    _e=_calc_gdihs_gaff(csA,snA,&de,&mol->ff_d[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_d[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_d[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_d[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_d[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_d[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_d[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_d[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_d[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_d[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_d[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_d[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_d[_i].atom[3]][_L].k+=de*_g[11];
    }
  }
else
  {
  _j=activem->active_d->size;
  while (_j--)
    {
    _i=activem->active_d->list[_j];
    calc_dih_derivative(&csA,&snA,&r[mol->ff_d[_i].atom[0]],&r[mol->ff_d[_i].atom[1]],&r[mol->ff_d[_i].atom[2]],&r[mol->ff_d[_i].atom[3]],_g);
    _e=_calc_gdihs_gaff(csA,snA,&de,&mol->ff_d[_i]);
    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
    e[_L]+=_e,
    g[mol->ff_d[_i].atom[0]][_L].i+=de*_g[ 0], g[mol->ff_d[_i].atom[0]][_L].j+=de*_g[ 1], g[mol->ff_d[_i].atom[0]][_L].k+=de*_g[ 2],
    g[mol->ff_d[_i].atom[1]][_L].i+=de*_g[ 3], g[mol->ff_d[_i].atom[1]][_L].j+=de*_g[ 4], g[mol->ff_d[_i].atom[1]][_L].k+=de*_g[ 5],
    g[mol->ff_d[_i].atom[2]][_L].i+=de*_g[ 6], g[mol->ff_d[_i].atom[2]][_L].j+=de*_g[ 7], g[mol->ff_d[_i].atom[2]][_L].k+=de*_g[ 8],
    g[mol->ff_d[_i].atom[3]][_L].i+=de*_g[ 9], g[mol->ff_d[_i].atom[3]][_L].j+=de*_g[10], g[mol->ff_d[_i].atom[3]][_L].k+=de*_g[11];
    }
  }
}

//------------------------------------    N O N B O N D E D   ------------------------------------------

//------------------------------------    N O N B O N D E D   P R I M I T I V E S   ------------------------------------------

//This function calculates energy and force of nonbonded interatons in YFF1. The flag defines if energy(forces) adds of subtracts
//Note g_vecs and r_vecs migh be ffsys->g and ffsys->r correspondingly
inline double calc_atom__nb_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top)
{
double _r, _rr, _rrrrrr, _a, _b, _q, _d;
t_vec _u;
_u.i=r_i->i-r_j->i, _u.j=r_i->j-r_j->j, _u.k=r_i->k-r_j->k, _rr=_u.i*_u.i+_u.j*_u.j+_u.k*_u.k+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _a=top->A[a_i][a_j], _b=top->B[a_i][a_j], _q=COULOMB_K*q_i*q_j, _rrrrrr=_rr*_rr*_rr, _r=sqrt(_rr);
  if (_r<YFF1_Rb)
    {//Do classic interaction
    _d=(-12.*_a/_rrrrrr+6.*_b)/_rrrrrr/_rr-_q/_rr/_r;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return (_a/_rrrrrr-_b)/_rrrrrr+_q/_r;
    }
  else 
    {//Shift interaction function
    _d=(_a*YFF1_AF(_r,_rr,_rrrrrr)-_b*YFF1_BF(_r,_rr,_rrrrrr)+_q*YFF1_QF(_r,_rr,_rrrrrr))/_r;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return _a*YFF1_AV(_r,_rr,_rrrrrr)-_b*YFF1_BV(_r,_rr,_rrrrrr)+_q*YFF1_QV(_r,_rr,_rrrrrr);
    }
  }
else return 0.;
}
inline double calc_atom__nb_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top)
{
double _r, _rr, _rrrrrr, _a, _b, _q;
_rr=sqrd(r_i->i-r_j->i)+sqrd(r_i->j-r_j->j)+sqrd(r_i->k-r_j->k)+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _a=top->A[a_i][a_j], _b=top->B[a_i][a_j], _q=COULOMB_K*q_i*q_j, _rrrrrr=_rr*_rr*_rr, _r=sqrt(_rr);
  if (_r<YFF1_Rb) return (_a/_rrrrrr-_b)/_rrrrrr+_q/_r;                                                    //Do classic interaction
  else            return _a*YFF1_AV(_r,_rr,_rrrrrr)-_b*YFF1_BV(_r,_rr,_rrrrrr)+_q*YFF1_QV(_r,_rr,_rrrrrr); //Do shifted function
  }
else return 0.;
}
inline double calc_atom__znb_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top)
{
double _r, _rr, _rrrrrr, _a, _b, _q, _d;
t_vec _u;
_u.i=r_i->i-r_j->i, _u.j=r_i->j-r_j->j, _u.k=r_i->k-r_j->k, _rr=_u.i*_u.i+_u.j*_u.j+_u.k*_u.k+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _a=top->A[a_i][0]*top->A[a_j][0], _b=top->B[a_i][0]*top->B[a_j][0], _q=COULOMB_K*q_i*q_j, _rrrrrr=_rr*_rr*_rr, _r=sqrt(_rr);
  if (_r<YFF1_Rb)
    {//Do classic interaction
    _d=(-12.*_a/_rrrrrr+6.*_b)/_rrrrrr/_rr-_q/_rr/_r;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return (_a/_rrrrrr-_b)/_rrrrrr+_q/_r;
    }
  else 
    {//Shift interaction function
    _d=(_a*YFF1_AF(_r,_rr,_rrrrrr)-_b*YFF1_BF(_r,_rr,_rrrrrr)+_q*YFF1_QF(_r,_rr,_rrrrrr))/_r;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return _a*YFF1_AV(_r,_rr,_rrrrrr)-_b*YFF1_BV(_r,_rr,_rrrrrr)+_q*YFF1_QV(_r,_rr,_rrrrrr);
    }
  }
else return 0.;
}
inline double calc_atom__znb_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top)
{
double _r, _rr, _rrrrrr, _a, _b, _q;
_rr=sqrd(r_i->i-r_j->i)+sqrd(r_i->j-r_j->j)+sqrd(r_i->k-r_j->k)+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _a=top->A[a_i][0]*top->A[a_j][0], _b=top->B[a_i][0]*top->B[a_j][0], _q=COULOMB_K*q_i*q_j, _rrrrrr=_rr*_rr*_rr, _r=sqrt(_rr);
  if (_r<YFF1_Rb) return (_a/_rrrrrr-_b)/_rrrrrr+_q/_r;                                                     //Do classic interaction
  else            return _a*YFF1_AV(_r,_rr,_rrrrrr)-_b*YFF1_BV(_r,_rr,_rrrrrr)+_q*YFF1_QV(_r,_rr,_rrrrrr);  //Do shifted function
  }
else return 0.;
}
//This function calculates 1-4 coulomb interactions only
inline double calc_atom__14_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top)
{
double _r, _rr, _q, _d;
t_vec _u;
_u.i=r_i->i-r_j->i, _u.j=r_i->j-r_j->j, _u.k=r_i->k-r_j->k, _rr=_u.i*_u.i+_u.j*_u.j+_u.k*_u.k+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _q=COULOMB_K*q_i*q_j, _r=sqrt(_rr);
  if (_r<YFF1_Rb)
    {//Do classic interaction
    _d=-YFF1_14_SCALE*_q/_rr/_r;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return YFF1_14_SCALE*_q/_r;
    }
  else 
    {//Shift interaction function
    _d=YFF1_14_SCALE*_q*YFF1_QF(_r,_rr,0.)/_r;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return YFF1_14_SCALE*_q*YFF1_QV(_r,_rr,0.);
    }
  }
else return 0.;
}
inline double calc_atom__14_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top)
{
double _r, _rr, _q;
_rr=sqrd(r_i->i-r_j->i)+sqrd(r_i->j-r_j->j)+sqrd(r_i->k-r_j->k)+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _q=COULOMB_K*q_i*q_j, _r=sqrt(_rr);
  if (_r<YFF1_Rb) return YFF1_14_SCALE*_q/_r;                       //Do classic interaction
  else            return YFF1_14_SCALE*_q*YFF1_QV(_r,_rr,_rrrrrr);  //Do shifted function
  }
else return 0.;
}
//This function calculates negative difference between 1-4 and normal coulomb interactions only
inline double calc_atom_n14_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_vec *g_i,t_vec *g_j,t_top *top)
{
double _r, _rr, _q, _d;
t_vec _u;
_u.i=r_i->i-r_j->i, _u.j=r_i->j-r_j->j, _u.k=r_i->k-r_j->k, _rr=_u.i*_u.i+_u.j*_u.j+_u.k*_u.k+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _q=COULOMB_K*q_i*q_j, _r=sqrt(_rr);
  if (_r<YFF1_Rb)
    {//Do classic interaction
    _d=-(1.-YFF1_14_SCALE)*_q/_rr/_r;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return (1.-YFF1_14_SCALE)*_q/_r;
    }
  else 
    {//Shift interaction function
    _d=(1.-YFF1_14_SCALE)*_q*YFF1_QF(_r,_rr,0.)/_r;
    _u.i*=_d, _u.j*=_d, _u.k*=_d, g_i->i=+_u.i, g_i->j=+_u.j, g_i->k=+_u.k, g_j->i=-_u.i, g_j->j=-_u.j, g_j->k=-_u.k;
    return (1.-YFF1_14_SCALE)*_q*YFF1_QV(_r,_rr,0.);
    }
  }
else return 0.;
}
inline double calc_atom_n14_enrg_yff1(unsigned int a_i,unsigned int a_j,double q_i,double q_j,t_vec *r_i,t_vec *r_j,t_top *top)
{
double _r, _rr, _q;
_rr=sqrd(r_i->i-r_j->i)+sqrd(r_i->j-r_j->j)+sqrd(r_i->k-r_j->k)+SMALL2; //To keep divisibility
if (_rr<YFF1_Rc2)
  {
  _q=COULOMB_K*q_i*q_j, _r=sqrt(_rr);
  if (_r<YFF1_Rb) return (1.-YFF1_14_SCALE)*_q/_r;                  //Do classic interaction
  else            return (1.-YFF1_14_SCALE)*_q*YFF1_QV(_r,_rr,0.);  //Do shifted function
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
inline double calc_atom__znb_enrg_yff1_on_grid(double a,double b,double q,unsigned int a_i,double q_i,t_top *top)
{
return top->A[a_i][0]*a-top->B[a_i][0]*b+q_i*q; 
}

//------------------------------------  E N E R G Y   F U N C T I O N S   Y F F 1  ------------------------------------------

//This function calculates energy in the active layer and marks atoms to exclude from grid summing
inline void calc_nbenergy_yff1(double e[N_ESLICES],t_ffsys *ffsys,t_top *top)
{
unsigned int _i, _j, _k, _l, mol_i, mol_j, active_i, active_j, a_i, a_j, _L;
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
              if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e;
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
              if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e;
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
              if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e;
                }
              }
            }
          }
        }
      }
    //Do inside mol: remove anchor-[anchor == root]-anchor and remove pairs
    if ((active_i=ffsys->mols[mol_i]->anchors->size)!=1)
      while (active_i--)
        {
        //Do root-rest
        _j=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
        if (active_i!=ffsys->rtrees[mol_i]->root) a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0];
        else a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch].edge.vertice[1];                                     
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e;
            }
          }
        //Do rest-rest
        _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
        while (--_i)
          {
          a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]].edge.vertice[1];
          _j=_i;
          while (--_j)  
            {
            a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
            if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
              {
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]-=_e;
              }
            }
          }
        }
    //Edit 1-4
    active_i=ffsys->mols[mol_i]->size_d;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_d[active_i].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_d[active_i].atom[3]; 
      if ( (_e=calc_atom_n14_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
        {
        if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
        else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
        e[_L]-=_e;
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
            if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
              {
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]+=_e;
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
              if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e;
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
                     if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                       {
                       if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                       else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                       e[_L]+=_e;
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
                     if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                       {
                       if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                       else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                       e[_L]+=_e;
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
                   if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                     {
                     if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                     else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                     e[_L]+=_e;
                     } 
                   }
                 }
               }
        }
      }
    //Do inside mol: remove anchor-[anchor == root]-anchor and remove pairs
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      //Do root-rest
      _j=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
      a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0];
      if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch>0)
        {
        _l=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (--_l)
          {
          a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_l];
          if ( (_e=calc_atom__znb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e;
            }
          }
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__znb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e;
            }
          }
        }
      else
        while (--_j)
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];  
          if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e;
            }
          }
      //Do rest-rest
      _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
      while (--_i)
        {
        a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]].edge.vertice[1];
        _j=_i;
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e;
            }
          }
        }
      }
    //Edit 1-4
    _i=ffsys->activem[mol_i].active_a->size; while (_i--) { active_i=ffsys->activem[mol_i].active_a->list[_i]; if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch>0) { a_i=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0], ffsys->mols[mol_i]->atoms->list[a_i]=(unsigned int)-ffsys->mols[mol_i]->atoms->list[a_i]; } } //Mark roots
    active_i=ffsys->activem[mol_i].active_d->size;
    while (active_i--)
      if ( ((int)ffsys->mols[mol_i]->atoms->list[a_i=ffsys->mols[mol_i]->ff_d[ffsys->activem[mol_i].active_d->list[active_i]].atom[0]]>0)&&((int)ffsys->mols[mol_i]->atoms->list[a_j=ffsys->mols[mol_i]->ff_d[ffsys->activem[mol_i].active_d->list[active_i]].atom[3]]>0)&&
           ((int)ffsys->mols[mol_i]->atoms->list[    ffsys->mols[mol_i]->ff_d[ffsys->activem[mol_i].active_d->list[active_i]].atom[1]]>0)&&((int)ffsys->mols[mol_i]->atoms->list[    ffsys->mols[mol_i]->ff_d[ffsys->activem[mol_i].active_d->list[active_i]].atom[2]]>0) ) 
        {//Correct only torsions in pure active layer
        a_i+=ffsys->nr[mol_i], a_j+=ffsys->nr[mol_i];
        if ( (_e=calc_atom_n14_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
          {
          if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
          else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
          e[_L]-=_e;
          }
        }
    _i=ffsys->activem[mol_i].active_a->size; while (_i--) { active_i=ffsys->activem[mol_i].active_a->list[_i]; if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch>0) { a_i=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0], ffsys->mols[mol_i]->atoms->list[a_i]=(unsigned int)-ffsys->mols[mol_i]->atoms->list[a_i]; } } //Unmark roots
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
            if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
              {
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]+=_e;
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
              if ( (_e=calc_atom__nb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e;
                }   
              } 
            }
          }
        }  
      }
    }
}
//Note it scales repulsion as E*=E*(3*x^2-2*x^3) in range 0 @ px==0. ... E @ px=1. dE*/dr=(3*x^2-2*x^3)*dE/dr
//This function calculates energy in the active layer and marks atoms to exclude from grid summing
inline void calc_nbgrad_yff1(double e[N_ESLICES],t_ffsys *ffsys,t_top *top)
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
    //Do inside mol: remove anchor-[anchor == root]-anchor and remove pairs
    if ((active_i=ffsys->mols[mol_i]->anchors->size)!=1)
      while (active_i--)
        {
        //Do root-rest
        _j=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
        if (active_i!=ffsys->rtrees[mol_i]->root) a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0];
        else a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch].edge.vertice[1];                                     
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
            }
          }
        //Do rest-rest
        _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
        while (--_i)
          {
          a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]].edge.vertice[1];
          _j=_i;
          while (--_j)  
            {
            a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
            if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
              {
              if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
              else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
              e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
              }
            }
          }
        }
    //Edit 1-4
    active_i=ffsys->mols[mol_i]->size_d;
    while (active_i--)
      {
      a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_d[active_i].atom[0],
      a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->ff_d[active_i].atom[3]; 
      if ( (_e=calc_atom_n14_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
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
      }
    //Do inside mol: remove anchor-[anchor == root]-anchor and remove pairs
    _k=ffsys->activem[mol_i].active_a->size;
    while (_k--)
      {
      active_i=ffsys->activem[mol_i].active_a->list[_k];
      //Do root-rest
      _j=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
      a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0];
      if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch>0)
        {
        _l=ffsys->mols[mol_i]->anchors->list[active_i].size;
        while (--_l)
          {
          a_j=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_l];
          if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
            }
          }
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__znb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
            }
          }
        }
      else
        while (--_j)
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];  
          if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
            }
          }
      //Do rest-rest
      _i=ffsys->rtrees[mol_i]->rbranch[active_i].nrbranch;
      while (--_i)
        {
        a_i=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_i]].edge.vertice[1];
        _j=_i;
        while (--_j)  
          {
          a_j=ffsys->nr[mol_i]+ffsys->rtrees[mol_i]->rbranch[ffsys->rtrees[mol_i]->rbranch[active_i].rbranch[_j]].edge.vertice[1];
          if ( (_e=calc_atom__nb_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
            }
          }
        }
      }
    //Edit 1-4
    _i=ffsys->activem[mol_i].active_a->size; while (_i--) { active_i=ffsys->activem[mol_i].active_a->list[_i]; if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch>0) { a_i=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0], ffsys->mols[mol_i]->atoms->list[a_i]=(unsigned int)-ffsys->mols[mol_i]->atoms->list[a_i]; } } //Mark roots
    active_i=ffsys->activem[mol_i].active_d->size;
    while (active_i--)
      if ( ((int)ffsys->mols[mol_i]->atoms->list[a_i=ffsys->mols[mol_i]->ff_d[ffsys->activem[mol_i].active_d->list[active_i]].atom[0]]>0)&&((int)ffsys->mols[mol_i]->atoms->list[a_j=ffsys->mols[mol_i]->ff_d[ffsys->activem[mol_i].active_d->list[active_i]].atom[3]]>0)&&
           ((int)ffsys->mols[mol_i]->atoms->list[    ffsys->mols[mol_i]->ff_d[ffsys->activem[mol_i].active_d->list[active_i]].atom[1]]>0)&&((int)ffsys->mols[mol_i]->atoms->list[    ffsys->mols[mol_i]->ff_d[ffsys->activem[mol_i].active_d->list[active_i]].atom[2]]>0) ) 
        {//Correct only torsions in pure active layer
        a_i+=ffsys->nr[mol_i], a_j+=ffsys->nr[mol_i];
        if ( (_e=calc_atom_n14_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],&g_i,&g_j,top)))
          {
          if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
          else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
          e[_L]-=_e, ffsys->_g[a_i][_L].i-=g_i.i, ffsys->_g[a_i][_L].j-=g_i.j, ffsys->_g[a_i][_L].k-=g_i.k, ffsys->_g[a_j][_L].i-=g_j.i, ffsys->_g[a_j][_L].j-=g_j.j, ffsys->_g[a_j][_L].k-=g_j.k;
          }
        }
    _i=ffsys->activem[mol_i].active_a->size; while (_i--) { active_i=ffsys->activem[mol_i].active_a->list[_i]; if ((int)*ffsys->rtrees[mol_i]->rbranch[active_i].rbranch>0) { a_i=ffsys->rtrees[mol_i]->rbranch[active_i].edge.vertice[0], ffsys->mols[mol_i]->atoms->list[a_i]=(unsigned int)-ffsys->mols[mol_i]->atoms->list[a_i]; } } //Unmark roots
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
}

//!!!!!!!!!!!!! This function can be improved by excluding the exact anchor in explicite summation of ongrid failure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculates yff1 energy on dgrid
//It uses marks made by calc_yff1_grad(t_ffsys *ffsys,t_top *top) to skip frozen atoms
inline void calc_energy_on_grid_yff1(double e[N_ESLICES],t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***B,double ***Q)
{
unsigned int _i, _j, _k, _l, _L, mol_i, mol_j, active_i, active_j, a_i, a_j;
double _e, a, b, q;

mol_i=ffsys->nmols;
while (mol_i--)
  if (!(ffsys->activem[mol_i].active_a))
    {
    _i=ffsys->mols[mol_i]->atoms->size;
    while (_i--)
      {
      a_i=ffsys->nr[mol_i]+_i;
      if ( ( (calc_tricubic_interpolation_wp_monotonicity(&a,ori,sp,ni,nj,nk,&ffsys->r[a_i],A)))&&
           ( (calc_tricubic_interpolation_wp_monotonicity(&b,ori,sp,ni,nj,nk,&ffsys->r[a_i],B)))&&
           ( (calc_tricubic_interpolation_wp_monotonicity(&q,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) )
        {
        //Summ on grid
        if ( (_e=top->A[(unsigned int)ffsys->a[a_i]][0]*a-top->B[(unsigned int)ffsys->a[a_i]][0]*b+ffsys->q[a_i]*q))
          { 
          if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
          else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
          e[_L]+=_e;
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
                if ( (_e=calc_atom__znb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]-=_e;
                  }
                }
              }
            //Summ all 
            _j=ffsys->mols[mol_j]->atoms->size;
            while (_j--)
              {
              a_j=ffsys->nr[mol_j]+_j;
              if ( (_e=calc_atom__znb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                {
                if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                e[_L]+=_e;
                }
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
        if ( ( (calc_tricubic_interpolation_wp_monotonicity(&a,ori,sp,ni,nj,nk,&ffsys->r[a_i],A)))&&
             ( (calc_tricubic_interpolation_wp_monotonicity(&b,ori,sp,ni,nj,nk,&ffsys->r[a_i],B)))&&
             ( (calc_tricubic_interpolation_wp_monotonicity(&q,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) )
          {
          //Summ on grid
          if ( (_e=top->A[(unsigned int)ffsys->a[a_i]][0]*a-top->B[(unsigned int)ffsys->a[a_i]][0]*b+ffsys->q[a_i]*q))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]+=_e;
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
                  if ( (_e=calc_atom__znb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]-=_e;
                    }
                  }
                }
              //Summ all
              _j=ffsys->mols[mol_j]->atoms->size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+_j;
                if ( (_e=calc_atom__znb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e;
                  }
                }
              }
          }
        }
     //Summ the rest 
     _i=ffsys->mols[mol_i]->anchors->list[active_i].size;
      while (--_i)
        {
        a_i=ffsys->nr[mol_i]+ffsys->mols[mol_i]->anchors->list[active_i].list[_i];
        if ( ( (calc_tricubic_interpolation_wp_monotonicity(&a,ori,sp,ni,nj,nk,&ffsys->r[a_i],A)))&&
             ( (calc_tricubic_interpolation_wp_monotonicity(&b,ori,sp,ni,nj,nk,&ffsys->r[a_i],B)))&&
             ( (calc_tricubic_interpolation_wp_monotonicity(&q,ori,sp,ni,nj,nk,&ffsys->r[a_i],Q))) )
          {
          //Summ on grid
          if ( (_e=top->A[(unsigned int)ffsys->a[a_i]][0]*a-top->B[(unsigned int)ffsys->a[a_i]][0]*b+ffsys->q[a_i]*q))
            {
            if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
            else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
            e[_L]+=_e;
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
                  if ( (_e=calc_atom__znb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                    {
                    if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                    else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                    e[_L]-=_e;
                    }
                  }
                }
              //Summ all 
              _j=ffsys->mols[mol_j]->atoms->size;
              while (_j--)
                {
                a_j=ffsys->nr[mol_j]+_j;
                if ( (_e=calc_atom__znb_enrg_yff1(ffsys->a[a_i],ffsys->a[a_j],ffsys->q[a_i],ffsys->q[a_j],&ffsys->r[a_i],&ffsys->r[a_j],top)))
                  {
                  if (fabs(_e)>E_SERIES[_L=N_ESLICES/2]) { while (++_L<N_ESLICES) if (fabs(_e)<E_SERIES[_L]) break; _L--; }
                  else                                   { while (--_L)           if (fabs(_e)>E_SERIES[_L]) break;       }
                  e[_L]+=_e;
                  }
                }
              }
          }
        }
      }
    }
}
//!!!!!!!!!!!!! This function can be improved by excluding the exact anchor in explicite summation of ongrid failure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculates yff1 energy on dgrid
//It uses marks made by calc_yff1_grad(t_ffsys *ffsys,t_top *top) to skip frozen atoms
inline void calc_grad_on_grid_yff1(double e[N_ESLICES],t_ffsys *ffsys,t_top *top,t_vec *ori,unsigned int ni,unsigned int nj,unsigned int nk,double sp,double ***A,double ***B,double ***Q)
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
     //Summ the rest 
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

//This function calculates yff1 energy on grid
inline double calc_ffsys_nbenergy_on_dgrid_yff1(t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *B,t_dgrid *Q)
{
unsigned int _L;
//Init energies and forces
double  e[N_ESLICES]={[0 ... N_ESLICES-1]=0.};
//Summ active layers
calc_nbenergy_yff1(e,ffsys,top);
calc_dihs_energy_gaff(e,ffsys);
//_L=ffsys->rtrees[0]->nrbranch; while (_L--) if ((int)*ffsys->rtrees[0]->rbranch[_L].rbranch<0) *ffsys->rtrees[0]->rbranch[_L].rbranch=(unsigned int)(-(int)*ffsys->rtrees[0]->rbranch[_L].rbranch);
calc_energy_on_grid_yff1(e,ffsys,top,&A->ori,A->len.i,A->len.j,A->len.k,A->sp,A->d,B->d,Q->d);
//Summ energy and gradient
for (_L=1; _L<N_ESLICES; _L++) e[_L]+=e[_L-1];
return e[N_ESLICES-1];
}
//This function calculates yff1 energy gradient on grid
inline double calc_ffsys_nbgrad_on_dgrid_yff1(t_ffsys *ffsys,t_top *top,t_dgrid *A,t_dgrid *B,t_dgrid *Q)
{
unsigned int _i, _L;
//Init energies and forces
double  e[N_ESLICES]={[0 ... N_ESLICES-1]=0.};
memset(ffsys->_g,0x0,N_ESLICES*ffsys->natoms*sizeof(t_vec));
//Summ active layers
calc_nbgrad_yff1(e,ffsys,top);
calc_dihs_grad_gaff(e,ffsys);
//_L=ffsys->rtrees[0]->nrbranch; while (_L--) if ((int)*ffsys->rtrees[0]->rbranch[_L].rbranch<0) *ffsys->rtrees[0]->rbranch[_L].rbranch=(unsigned int)(-(int)*ffsys->rtrees[0]->rbranch[_L].rbranch);
calc_grad_on_grid_yff1(e,ffsys,top,&A->ori,A->len.i,A->len.j,A->len.k,A->sp,A->d,B->d,Q->d);
//Summ energy and gradient
for (_L=1; _L<N_ESLICES; _L++) e[_L]+=e[_L-1];
_i=ffsys->natoms; 
while (_i--) 
  {
  for (_L=1; _L<N_ESLICES; _L++) { ffsys->_g[_i][_L].i+=ffsys->_g[_i][_L-1].i, ffsys->_g[_i][_L].j+=ffsys->_g[_i][_L-1].j, ffsys->_g[_i][_L].k+=ffsys->_g[_i][_L-1].k; }
  ffsys->g[_i].i=ffsys->_g[_i][_L-1].i, ffsys->g[_i].j=ffsys->_g[_i][_L-1].j, ffsys->g[_i].k=ffsys->_g[_i][_L-1].k;
  }
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
  //Copy all coords to keep static layer
  _j=mol->atoms->size, r_vecs+=_j, r+=_j; while (_j--) { r_vecs--, r--, (*r_vecs).i=(*r).i, (*r_vecs).j=(*r).j, (*r_vecs).k=(*r).k; }
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
//NOTE. Torss in active list should be sorted previously.
//NB! This function is only valid for calculation around origin!
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

//----------------------------------   M I N I M I Z A T I O N     R O U T I N E S   ----------------------------------------//

//This function calculates energies in internal + RT DOFs
//Note. x here organized as: | INTERNAL | RT |
//NB! if molecule has active list than it do NOT has RT DOFs but has several built-in rtrees
char calc_ic_energy_on_dgrid_yff0(double *_e,unsigned int n,double *x,va_list stack)
{
unsigned int _id, mol_id;
t_vec *rvecs, *cm;
t_ffsys *ffsys;
t_dgrid *A, *Q;
t_top *top;
t_tensor *R;
void *rp;
va_list _stack;

va_copy(_stack,stack);

//Stage 0. Unwrap stack
ffsys=va_arg(_stack,t_ffsys*);
A=va_arg(_stack,t_dgrid*);
Q=va_arg(_stack,t_dgrid*);
rvecs=va_arg(_stack,t_vec*);
cm=va_arg(_stack,t_vec*);
R=va_arg(_stack,t_tensor*);
top=va_arg(_stack,t_top*);

//Stage 1. Construct get molecular coords
for (_id=0, mol_id=0; mol_id<ffsys->nmols; mol_id++)
  _id+=construct_internal_coords_mols(&x[_id],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->r[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);

//Stage 2. Do nonbonded linear gradient & energy calculation
rp=ffsys->r, ffsys->r=rvecs, rvecs=rp; 
*_e=calc_ffsys_nbenergy_on_dgrid_yff0(ffsys,top,A,Q); 
rp=rvecs, rvecs=ffsys->r, ffsys->r=rp;

va_end(_stack);
return TRUE;
}
//This function calculates derivatives over internal + RT DOFs
//Note. x and g here organized as: | INTERNAL | RT |
//NB! if molecule has active list than it do NOT has RT DOFs but has several built-in rtrees
char calc_ic_grad_on_dgrid_yff0(double *_e,unsigned int n,double *x,double *g,double **G,va_list stack)
{
unsigned int _id, mol_id;
t_vec *rvecs, *cm;
t_ffsys *ffsys;
t_dgrid *A, *Q;
t_top *top;
t_tensor *R;
void *rp;
va_list _stack;

va_copy(_stack,stack);

//Stage 0. Unwrap stack
ffsys=va_arg(_stack,t_ffsys*);
A=va_arg(_stack,t_dgrid*);
Q=va_arg(_stack,t_dgrid*);
rvecs=va_arg(_stack,t_vec*);
cm=va_arg(_stack,t_vec*);
R=va_arg(_stack,t_tensor*);
top=va_arg(_stack,t_top*);

//Stage 1. Construct get molecular coords
for (_id=0, mol_id=0; mol_id<ffsys->nmols; mol_id++)
  _id+=construct_internal_coords_mols(&x[_id],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->r[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);

//Stage 2. Do nonbonded linear gradient & energy calculation
rp=ffsys->r, ffsys->r=rvecs, rvecs=rp; 
*_e=calc_ffsys_nbgrad_on_dgrid_yff0(ffsys,top,A,Q); 
rp=rvecs, rvecs=ffsys->r, ffsys->r=rp;

//Stage 3. Convert linear into internal-coords gradient
while(n--) g[n]=0.;
for (_id=mol_id=0; mol_id<ffsys->nmols; mol_id++)
  {
  if (!ffsys->activem[mol_id].active_a)
    { //Stage a. Calculate RT gradients (using quaternions)
    convert_RT_coords_grad(&g[_id],(t_vec*)&x[_id+0],(t_vec*)&x[_id+3],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->g[ffsys->nr[mol_id]],ffsys->mols[mol_id]);
    //for DEBUGGING g[_id+0]=g[_id+1]=g[_id+2]=g[_id+3]=g[_id+4]=g[_id+5]=0.;
    _id+=6; //Store rotational quaternion
    }
  //Stage b. Calculate torsional gradients and update with soft-bonded energy
  convert_ic_grad(&g[_id],&rvecs[ffsys->nr[mol_id]],&ffsys->g[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);
  _id+=ffsys->rtrees[mol_id]->nidofs;
  }

va_end(_stack);
return TRUE;
}


//This function calculates energies in internal + RT DOFs
//Note. x here organized as: | INTERNAL | RT |
//NB! if molecule has active list than it do NOT has RT DOFs but has several built-in rtrees
char calc_ic_energy_on_dgrid_yff1(double *_e,unsigned int n,double *x,va_list stack)
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
rp=ffsys->r, ffsys->r=rvecs, rvecs=rp; 
*_e=calc_ffsys_nbenergy_on_dgrid_yff1(ffsys,top,A,B,Q); 
rp=rvecs, rvecs=ffsys->r, ffsys->r=rp;

va_end(_stack);
return TRUE;
}
//This function calculates derivatives over internal + RT DOFs
//Note. x and g here organized as: | INTERNAL | RT |
//NB! if molecule has active list than it do NOT has RT DOFs but has several built-in rtrees
char calc_ic_grad_on_dgrid_yff1(double *_e,unsigned int n,double *x,double *g,double **G,va_list stack)
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
rp=ffsys->r, ffsys->r=rvecs, rvecs=rp; 
*_e=calc_ffsys_nbgrad_on_dgrid_yff1(ffsys,top,A,B,Q); 
rp=rvecs, rvecs=ffsys->r, ffsys->r=rp;

//Stage 3. Convert linear into internal-coords gradient
while(n--) g[n]=0.;
for (_id=mol_id=0; mol_id<ffsys->nmols; mol_id++)
  {
  if (!ffsys->activem[mol_id].active_a)
    { //Stage a. Calculate RT gradients (using quaternions)
    convert_RT_coords_grad(&g[_id],(t_vec*)&x[_id+0],(t_vec*)&x[_id+3],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->g[ffsys->nr[mol_id]],ffsys->mols[mol_id]);
    //for DEBUGGING g[_id+0]=g[_id+1]=g[_id+2]=g[_id+3]=g[_id+4]=g[_id+5]=0.;
    _id+=6; //Store rotational quaternion
    }
  //Stage b. Calculate torsional gradients and update with soft-bonded energy
  convert_ic_grad(&g[_id],&rvecs[ffsys->nr[mol_id]],&ffsys->g[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);
  _id+=ffsys->rtrees[mol_id]->nidofs;
  }

va_end(_stack);
return TRUE;
}

//This function update rotational coordinates to be always withing radii<=PI sphere
char sync_ic_coords(double *_e,unsigned int n,double *x,va_list stack)
{
unsigned int _id, mol_id;
double _d;
t_ffsys *ffsys;
va_list _stack;

va_copy(_stack,stack);

//Stage 0. Unwrap stack
ffsys=va_arg(_stack,t_ffsys*);

//Stage 1. Sync rotational vector
for (_id=mol_id=0; mol_id<ffsys->nmols; mol_id++)
  {
  if (!ffsys->activem[mol_id].active_a)
    {
    if ((_d=calc_vec_norm((t_vec*)&x[_id]))>=9.*PI*PI/4.) 
      multiple_vec_scalar((t_vec*)&x[_id],(t_vec*)&x[_id],.5-PI/sqrt(_d));
    _id+=6;
    }
  _id+=ffsys->rtrees[mol_id]->nidofs;
  }

va_end(_stack);
return TRUE;
} 

//--------------------------     L I N E A R     C O O R D I N A T E S    P A R T   -----------------------------------------------------

//------------------------------- G E N E R A L - P U R P O S E     M I N I M I Z A T I O N     R O U T I N E S  ----------------------------------------------

//This function gather statistics for IC system
void get_ic_complexity(unsigned int *n,unsigned int *m,unsigned int *naatoms,t_ffsys *ffsys)
{
unsigned int _i, _j, _k;
//Get *n (amount of idofs) and *m (longest path)
*naatoms=*n=*m=0; _j=0, _i=ffsys->nmols; 
while (_i--)
       if (!ffsys->activem[_i].active_a)
         {
         (*n)+=6+ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->root].nrbranch; if (!*m) *m=1; 
         ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->root].edge.type=ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->root].nrbranch;
         *naatoms=ffsys->mols[_i]->anchors->list[ffsys->rtrees[_i]->root].size;
         while (ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->root].edge.type--)
           {
           _k=1, _j=ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->root].rbranch[ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->root].edge.type];
           *naatoms+=ffsys->mols[_i]->anchors->list[_j].size;
           RESET_SCANN_0:
           ffsys->rtrees[_i]->rbranch[_j].edge.type=ffsys->rtrees[_i]->rbranch[_j].nrbranch;
           do{
             while (--ffsys->rtrees[_i]->rbranch[_j].edge.type)     
               {
               (*n)++; if (++_k>*m) *m=_k;
               _j=ffsys->rtrees[_i]->rbranch[_j].rbranch[ffsys->rtrees[_i]->rbranch[_j].edge.type];
               *naatoms+=ffsys->mols[_i]->anchors->list[_j].size;
               goto RESET_SCANN_0;
               }
             _k--, _j=*ffsys->rtrees[_i]->rbranch[_j].rbranch;
             } while (_j!=ffsys->rtrees[_i]->root) ;
           }
         }
  else if (ffsys->activem[_i].active_a->size)
         {
         (*n)+=ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->nrbranch].nrbranch; if (*m<2) *m=2;
         ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->nrbranch].edge.type=ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->nrbranch].nrbranch;
         while (ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->nrbranch].edge.type--)
           {
           _k=2, _j=ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->nrbranch].rbranch[ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->nrbranch].edge.type];
           *naatoms+=ffsys->mols[_i]->anchors->list[_j].size;
           RESET_SCANN_1:
           ffsys->rtrees[_i]->rbranch[_j].edge.type=ffsys->rtrees[_i]->rbranch[_j].nrbranch;
           do{
             while (--ffsys->rtrees[_i]->rbranch[_j].edge.type)     
               {
               (*n)++; if (++_k>*m) *m=_k;
               _j=ffsys->rtrees[_i]->rbranch[_j].rbranch[ffsys->rtrees[_i]->rbranch[_j].edge.type];
               *naatoms+=ffsys->mols[_i]->anchors->list[_j].size;
               goto RESET_SCANN_1;
               }
             _k--, _j=*ffsys->rtrees[_i]->rbranch[_j].rbranch;
             } while (_j!=*ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->nrbranch].rbranch[ffsys->rtrees[_i]->rbranch[ffsys->rtrees[_i]->nrbranch].edge.type]].rbranch) ;
           }
         }
}

//This function initializes ic massives for use in functions below
char initialize_ic_minimizer(unsigned int method,unsigned int nwats,unsigned int *n,unsigned int *m,unsigned int *naatoms,
                             t_tensor **R,t_vec **cm,t_vec **rvec,double *x[2],double *g[2],double ***G,double **p,t_ffsys *ffsys) 
{

//Stage I. Gather statistics
get_ic_complexity(n,m,naatoms,ffsys);
*n+=nwats*6, *naatoms+=nwats*3;

//Stage II. Allocate memory
switch (method)
  {
  case MINIMIZE_STEEPEST : { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
  case MINIMIZE_POLAKRB  : { 
    if (!((*p)=(double*)malloc(sizeof(t_tensor)**m+
                               sizeof(t_vec)*(ffsys->nmols+nwats)+
                               sizeof(t_vec)*(ffsys->natoms+nwats*3)+
                               sizeof(double)**n*5))) { ylib_errno=YERROR_MEMORY; return FALSE; }
    x[0]=*p+*n, x[1]=x[0]+*n, g[0]=x[1]+*n, g[1]=g[0]+*n; 
    *R=(void*)*p+sizeof(double)**n*5; 
    *cm=(void*)*R+sizeof(t_tensor)**m;
    *rvec=*cm+(ffsys->nmols+nwats);
    return TRUE;
    }
  case MINIMIZE_LBFGS    : 
  case MINIMIZE_NEWTON   : 
  default                : { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
  }

}

//This function minimizes system in internal coordinates on the grid using YFF0
char optimize_ic_ffsys_on_grid_yff0(double *e,unsigned int method,unsigned int nsteps,double tol,unsigned int n,double **x,double **g,double **G,double *p,unsigned int lbfgs_m,t_ffsys *ffsys,t_dgrid *A,t_dgrid *Q,t_vec *rvecs,t_vec *cm,t_tensor *R,t_top *top)
{
unsigned int mol_id, _i, _j;
double m, _m;
//Calculate cms
mol_id=ffsys->nmols;
while (mol_id--) 
  {
  cm[mol_id].i=cm[mol_id].j=cm[mol_id].k=0.; 
  if (!(ffsys->activem[mol_id].active_a))
    {
    m=0., _i=ffsys->mols[mol_id]->anchors->list[ffsys->rtrees[mol_id]->root].size;
    while (_i--) 
      { 
      _j=ffsys->nr[mol_id]+ffsys->mols[mol_id]->anchors->list[ffsys->rtrees[mol_id]->root].list[_i], 
      _m=top->ff_a[(unsigned int)ffsys->a[_j]].mass, m+=_m, cm[mol_id].i+=_m*ffsys->r[_j].i, cm[mol_id].j+=_m*ffsys->r[_j].j, cm[mol_id].k+=_m*ffsys->r[_j].k;
      }
    cm[mol_id].i/=m, cm[mol_id].j/=m, cm[mol_id].k/=m;
    }
  }
//Init x
memset(x[0],0x0,sizeof(double)*n);
//Do minimization
switch (method)
  {
  case MINIMIZE_STEEPEST : { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
  case MINIMIZE_POLAKRB  : { 
    if ( (!(polak_ribiere_sync(e,nsteps,tol,PI/6.,n,x,g,G,p,calc_ic_energy_on_dgrid_yff0,calc_ic_grad_on_dgrid_yff0,line_search_square_fapproximation,sync_ic_coords," ",ffsys,A,Q,rvecs,cm,R,top)))&&(ylib_errno!=YERROR_NCONVERGED) ) return FALSE;
    //Some improovement is gained, construct new coordinates
    for (_i=0, mol_id=0; mol_id<ffsys->nmols; mol_id++)
      _i+=construct_internal_coords_mols(&x[0x0][_i],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->r[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);
    _i=ffsys->natoms; while (_i--) { ffsys->r[_i].i=rvecs[_i].i, ffsys->r[_i].j=rvecs[_i].j, ffsys->r[_i].k=rvecs[_i].k; }
    return TRUE;
    }
  case MINIMIZE_LBFGS    : 
  case MINIMIZE_NEWTON   : 
  default                : { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
  }
}
//This function minimizes system in internal coordinates on the grid using YFF1 
char optimize_ic_ffsys_on_grid_yff1(double *e,unsigned int method,unsigned int nsteps,double tol,unsigned int n,double **x,double **g,double **G,double *p,unsigned int lbfgs_m,t_ffsys *ffsys,t_dgrid *A,t_dgrid *B,t_dgrid *Q,t_vec *rvecs,t_vec *cm,t_tensor *R,t_top *top)
{
unsigned int mol_id, _i, _j;
double m, _m;
//Calculate cms
mol_id=ffsys->nmols;
while (mol_id--) 
  {
  cm[mol_id].i=cm[mol_id].j=cm[mol_id].k=0.; 
  if (!(ffsys->activem[mol_id].active_a))
    {
    m=0., _i=ffsys->mols[mol_id]->anchors->list[ffsys->rtrees[mol_id]->root].size;
    while (_i--) 
      { 
      _j=ffsys->nr[mol_id]+ffsys->mols[mol_id]->anchors->list[ffsys->rtrees[mol_id]->root].list[_i], 
      _m=top->ff_a[(unsigned int)ffsys->a[_j]].mass, m+=_m, cm[mol_id].i+=_m*ffsys->r[_j].i, cm[mol_id].j+=_m*ffsys->r[_j].j, cm[mol_id].k+=_m*ffsys->r[_j].k;
      }
    cm[mol_id].i/=m, cm[mol_id].j/=m, cm[mol_id].k/=m;
    }
  }
//Init x
memset(x[0],0x0,sizeof(double)*n);
//Do minimization
switch (method)
  {
  case MINIMIZE_STEEPEST : { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
  case MINIMIZE_POLAKRB  : { 
    if ( (!(polak_ribiere_sync(e,nsteps,tol,PI/6.,n,x,g,G,p,calc_ic_energy_on_dgrid_yff1,calc_ic_grad_on_dgrid_yff1,line_search_square_fapproximation,sync_ic_coords," ",ffsys,A,B,Q,rvecs,cm,R,top)))&&(ylib_errno!=YERROR_NCONVERGED) ) return FALSE;
    //Some improovement is gained, construct new coordinates
    for (_i=0, mol_id=0; mol_id<ffsys->nmols; mol_id++)
      _i+=construct_internal_coords_mols(&x[0x0][_i],&cm[mol_id],&rvecs[ffsys->nr[mol_id]],&ffsys->r[ffsys->nr[mol_id]],ffsys->rtrees[mol_id],ffsys->mols[mol_id],R);
    _i=ffsys->natoms; while (_i--) { ffsys->r[_i].i=rvecs[_i].i, ffsys->r[_i].j=rvecs[_i].j, ffsys->r[_i].k=rvecs[_i].k; }
    return TRUE;
    }
  case MINIMIZE_LBFGS    : 
  case MINIMIZE_NEWTON   : 
  default                : { ylib_errno=YERROR_NIMPLEMENTED; return FALSE; }
  }
}





