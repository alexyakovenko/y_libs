  //This file contain routines for distribution function calculations
#include "y_Z.h"
#include <string.h>
#include <stdlib.h>

#define ROTABLE_BOND_LIMIT 1e-2          // 0.01A
#define COLINEAR_LIMIT 0.995562          // 5.4 deg == 0.03 PI


/*
//This function builds internal coordinates rotation tree for single fragment based on neighbouring and temp massive
unsigned int build_rotation_tree_with_neighbors(t_list *neighbors,unsigned int *step,unsigned int *stack)
{
register unsigned int _i=0x0,_j,_k;
register bool cycle;
unsigned int base;
//Take arbitrary lovest-valent point
_j=0x0;
while (neighbors->size>_j++)
  {
  _i=neighbors->size;
  while (_i--)
    if (((t_list*)neighbors->list[_i])->size==_j) goto EXIT_LOWVALENT;
  }
EXIT_LOWVALENT: ;
//Get first the longest vertice;
if ((_i=enumerate_vertices_with_neighbors(_i,(unsigned int)-1,neighbors,step,stack))==(unsigned int)-1) return FALSE; //Fragmewnt detected
else _i=enumerate_vertices_with_neighbors(stack[neighbors->size-0x1],(unsigned int)-1,neighbors,step,stack);  //Get the longest edge
//Get the middle edge of the graph
if (_i%2) _j=--_i/2; //not odd
else      _j=  _i/2; //yes odd
_i=neighbors->size;
while (--_i)
  if (step[_i]==_j) break;
//Reorder neighbors into tree
_j=enumerate_vertices_with_neighbors(base=_i,(unsigned int)-1,neighbors,step,stack);
//invert vectors
_i=neighbors->size;
while (_i--)
  {
  cycle=FALSE;
  _j=((t_list*)neighbors->list[_i])->size;
  while (--_j)
    if (step[_i]>=step[((t_list*)neighbors->list[_i])->list[_j]])
      {//Upsteream detected
      if (cycle)
        return (unsigned int)-1; //Cycles are not implemented jet
      else
        {//Set root as zero neighbor
        cycle=TRUE;
        _k=*((t_list*)neighbors->list[_i])->list;
        *((t_list*)neighbors->list[_i])->list=((t_list*)neighbors->list[_i])->list[_j];
        ((t_list*)neighbors->list[_i])->list[_j]=_k;
        }
      }
  }
//inform about base
return base;
}




//This function shows rotation units data
void show_runits(t_runits *runits,t_mol *mol)
{
unsigned int _i,_j;
for (_i=0;_i<runits->size;_i++)
  {
  printf("Runit %1d. Joinig anchors %d -- %d with middle (%f %f %f)\n",_i,mol->aedges[runits->aedges[_i]].vertice[0x0],mol->aedges[runits->aedges[_i]].vertice[0x1],runits->ru0[_i].i,runits->ru0[_i].j,runits->ru0[_i].k);
  for (printf("RU[0x0] is:"),_j=0;_j<runits->anchors[_i][0x0]->size;_j++)
    printf("%d ",runits->anchors[_i][0x0]->list[_j]);
  printf("\n");
  for (printf("RU[0x1] is:"),_j=0;_j<runits->anchors[_i][0x1]->size;_j++)
    printf("%d ",runits->anchors[_i][0x1]->list[_j]);
  printf("\n");
  }
}

//This function free memory from runits structure
void free_runits(t_runits *runits)
{
if (runits)
  {
  if (runits->anchors)
    {
    runits->anchors+=runits->size;
    while(runits->size--)
      {
      runits->anchors--;
      if ( (*runits->anchors)[0x0])
        free((*runits->anchors)[0x0]);
      if ( (*runits->anchors)[0x1])
        free((*runits->anchors)[0x1]);
      }
    free(runits->anchors);
    }
  if (runits->ru0) free(runits->ru0);
  if (runits->u0)  free(runits->u0);
  if (runits->aedges) free(runits->aedges);
  free(runits);
  }
}

//This function defines rotation units for given molecule. It returns the list of rotation units.
t_runits *define_runits(t_mol *mol)
{
unsigned int _i,_j,_k;
t_list *neighbors=0x0;
t_runits *runits=0x0;

//Stage 1. Init memory and build anchors graph
if (!(runits=(t_runits*)calloc(sizeof(t_runits),0x1))) goto ERROR;
if (!(runits->size=mol->naedges)) return runits; //nothing ot do - return empty structure
if ( (!(runits->aedges=(unsigned int*)malloc(sizeof(unsigned int)*runits->size)))||(!(runits->ru0=(t_vec*)malloc(sizeof(t_vec)*runits->size)))             ||
     (!(runits->u0=(t_vec(*)[0x2])malloc(sizeof(t_vec)*runits->size*0x2)))||(!(runits->anchors=(t_list *(*)[0x2])calloc(sizeof(t_list*),runits->size*0x2)))||
     (!(neighbors=(t_list*)define_neighbors(FALSE,mol->nanchors,mol->naedges,mol->aedges))) )
  goto ERROR;

//Stage 2. Generate rotable bond weighted anchors
_k=mol->naedges;
while(_k--)
  {
  runits->aedges[_k]=_k;
  runits->anchors[_k][0x0]=(t_list*)alloc_list(((t_list*)neighbors->list[mol->aedges[_k].vertice[0x0]])->size);
  //Init 0x0 entrance
  _j=0x1;
  _i=((t_list*)neighbors->list[mol->aedges[_k].vertice[0x0]])->size;
  while(_i--)
    if (((t_list*)neighbors->list[mol->aedges[_k].vertice[0x0]])->list[_i]==mol->aedges[_k].vertice[0x1])
      {
      while (_i--) //copy the rest
        runits->anchors[_k][0x0]->list[_j++]=((t_list*)neighbors->list[mol->aedges[_k].vertice[0x0]])->list[_i];
      *runits->anchors[_k][0x0]->list=mol->aedges[_k].vertice[0x0];
      break;
      }
    else runits->anchors[_k][0x0]->list[_j++]=((t_list*)neighbors->list[mol->aedges[_k].vertice[0x0]])->list[_i]; //copy unique item
  //Grow 0x0 entrance
  _i=0x0;
  while(++_i!=runits->anchors[_k][0x0]->size)
    {
    _j=((t_list*)neighbors->list[runits->anchors[_k][0x0]->list[_i]])->size;
    while (_j--)
      if (((t_list*)neighbors->list[runits->anchors[_k][0x0]->list[_i]])->list[_j]==mol->aedges[_k].vertice[0x1])
        goto ERROR; //Cycle detected!!!
      else if (!(find_in_list(FALSE,((t_list*)neighbors->list[runits->anchors[_k][0x0]->list[_i]])->list[_j],runits->anchors[_k][0x0])))
             if (!(list_add(((t_list*)neighbors->list[runits->anchors[_k][0x0]->list[_i]])->list[_j],&runits->anchors[_k][0x0]))) goto ERROR;
    }

  //Init 0x1 entrance
  runits->anchors[_k][0x1]=(t_list*)alloc_list(((t_list*)neighbors->list[mol->aedges[_k].vertice[0x1]])->size);
  _j=0x1;
  _i=((t_list*)neighbors->list[mol->aedges[_k].vertice[0x1]])->size;
  while(_i--)
    if (((t_list*)neighbors->list[mol->aedges[_k].vertice[0x1]])->list[_i]==mol->aedges[_k].vertice[0x0])
      {
      while (_i--) //copy the rest
        runits->anchors[_k][0x1]->list[_j++]=((t_list*)neighbors->list[mol->aedges[_k].vertice[0x1]])->list[_i];
      *runits->anchors[_k][0x1]->list=mol->aedges[_k].vertice[0x1];
      break;
      }
    else runits->anchors[_k][0x1]->list[_j++]=((t_list*)neighbors->list[mol->aedges[_k].vertice[0x1]])->list[_i]; //copy unique item
  _i=0x0;
  while(++_i!=runits->anchors[_k][0x1]->size)
    {
    _j=((t_list*)neighbors->list[runits->anchors[_k][0x1]->list[_i]])->size;
    while (_j--)
      if (((t_list*)neighbors->list[runits->anchors[_k][0x1]->list[_i]])->list[_j]==mol->aedges[_k].vertice[0x0])
        goto ERROR; //Cycle detected!!!
      else if (!(find_in_list(FALSE,((t_list*)neighbors->list[runits->anchors[_k][0x1]->list[_i]])->list[_j],runits->anchors[_k][0x1])))
             if (!(list_add(((t_list*)neighbors->list[runits->anchors[_k][0x1]->list[_i]])->list[_j],&runits->anchors[_k][0x1]))) goto ERROR;
    }
  }
//Return runits
free(neighbors);
return runits; //show_runits(runits,mol)
ERROR: free_runits(runits);
if (neighbors) free(neighbors);
return FALSE;
}

//This function define vector into runits
inline void calc_runits(register t_vec *r,register t_runits *runits,register t_mol *mol)
{
register unsigned int _k;
_k=runits->size;
while(_k--)
  {
  runits->u0[_k][0x0].i=r[mol->edges[mol->aedges[runits->aedges[_k]].type].vertice[0x0]].i;
  runits->u0[_k][0x0].j=r[mol->edges[mol->aedges[runits->aedges[_k]].type].vertice[0x0]].j;
  runits->u0[_k][0x0].k=r[mol->edges[mol->aedges[runits->aedges[_k]].type].vertice[0x0]].k;
  runits->ru0[_k].i=0.5*(runits->u0[_k][0x0].i+r[mol->edges[mol->aedges[runits->aedges[_k]].type].vertice[0x1]].i);
  runits->ru0[_k].j=0.5*(runits->u0[_k][0x0].j+r[mol->edges[mol->aedges[runits->aedges[_k]].type].vertice[0x1]].j);
  runits->ru0[_k].k=0.5*(runits->u0[_k][0x0].k+r[mol->edges[mol->aedges[runits->aedges[_k]].type].vertice[0x1]].k);
  runits->u0[_k][0x0].i-=runits->ru0[_k].i;
  runits->u0[_k][0x0].j-=runits->ru0[_k].j;
  runits->u0[_k][0x0].k-=runits->ru0[_k].k;
  multiple_vec_scalar(&runits->u0[_k][0x0],&runits->u0[_k][0x0],1.00/sqrt(calc_vec_norm(&runits->u0[_k][0x0])));
  runits->u0[_k][0x1].i=-runits->u0[_k][0x0].i;
  runits->u0[_k][0x1].j=-runits->u0[_k][0x0].j;
  runits->u0[_k][0x1].k=-runits->u0[_k][0x0].k;
  }
}


//This function calculates simply distribution over simple angle
double calc_diagonal_internal_dde(t_vec *r_i,t_vec *u_i,t_vec *ru_i,t_list *anchors_i,t_list *anchors_j,t_mol *mol_i,t_ffsys *ffsys)
{
double rra,rrb,r,csA,e,de,dde,dd=0.00;
t_vec ra,rb,rc,**r_j;
unsigned int _i,_j,_p,_q,nmols;
t_mol **mol_j;


//Process it list
e=0.00;
for (_i=0;_i<anchors_i->size;_i++)
  for (_p=0;_p<mol_i->anchors[anchors_i->list[_i]]->size;_p++)
    {
    //Calculate ra vector
    SUBT_VEC(((t_vec*)&rc),((t_vec*)&r_i[mol_i->anchors[anchors_i->list[_i]]->list[_p]]),ru_i); //rc temporary is ra0 vector
    multiple_vec_scalar(&ra,u_i,calc_vec_vec_scalar_product(u_i,&rc));
    SUBT_VEC(((t_vec*)&ra),((t_vec*)&rc),((t_vec*)&ra));
    rra=calc_vec_vec_scalar_product(&ra,&ra);
    if (rra>ROTABLE_BOND_LIMIT*ROTABLE_BOND_LIMIT) //Process iff the rotable bond really exists
      {
      r_j=ffsys->r;
      mol_j=ffsys->mols;
      nmols=ffsys->nmols;
      while(nmols--)
        {
        if ((*mol_j)==mol_i)
          {//summ over j-th anchor
          for (_j=0;_j<anchors_j->size;_j++)
            for (_q=0;_q<mol_i->anchors[anchors_j->list[_j]]->size;_q++)
              {
              //Calculate rb vector
              SUBT_VEC(((t_vec*)&rc),((t_vec*)&r_i[mol_i->anchors[anchors_j->list[_j]]->list[_q]]),ru_i);   //rc temporary is rb0 vector
              multiple_vec_scalar(&rb,u_i,calc_vec_vec_scalar_product(u_i,&rc));
              SUBT_VEC(((t_vec*)&rb),((t_vec*)&rc),((t_vec*)&rb));
              rrb=calc_vec_vec_scalar_product(&rb,&rb);
              if (rrb>ROTABLE_BOND_LIMIT*ROTABLE_BOND_LIMIT) //Process iff the rotable bond really exists
                {
                //Calculate centers distance and energy derivatives
                SUBT_VEC(((t_vec*)&rc),((t_vec*)&r_i[mol_i->anchors[anchors_i->list[_i]]->list[_p]]),((t_vec*)&r_i[mol_i->anchors[anchors_j->list[_j]]->list[_q]]));
                r=calc_vec_vec_scalar_product(&rc,&rc);
                if (r<FF_CUTOFF)
                  {
                  r=sqrt(r);
                  e+=ffsys->fnbonded(r,&de,&dde,ffsys->e[mol_i->atoms->list[mol_i->anchors[anchors_i->list[_i]]->list[_p]]][mol_i->atoms->list[mol_i->anchors[anchors_j->list[_j]]->list[_q]]],
                                                ffsys->R[mol_i->atoms->list[mol_i->anchors[anchors_i->list[_i]]->list[_p]]][mol_i->atoms->list[mol_i->anchors[anchors_j->list[_j]]->list[_q]]],
                                     mol_i->charges[mol_i->atoms->list[mol_i->anchors[anchors_i->list[_i]]->list[_p]]],mol_i->charges[mol_i->atoms->list[mol_i->anchors[anchors_j->list[_j]]->list[_q]]]);
                  de/=r;
                  dde=(dde-de)/r/r;
                  //Calculate relative angle
                  r=sqrt(rra*rrb);
                  csA=calc_vec_vec_scalar_product(&ra,&rb)/r;
                  //Calculate reference angle
                  dd+=dde*rra*rrb*(1.00-csA*csA)+de*r*csA;
                  }
                }
              }
          }
        else
          {//summ over the rest of molecules
          for (_j=0;_j<(*mol_j)->nanchors;_j++)
            for (_q=0;_q<(*mol_j)->anchors[_j]->size;_q++)
              {
              //Calculate rb vector
              SUBT_VEC(((t_vec*)&rc),((t_vec*)&(*r_j)[(*mol_j)->anchors[_j]->list[_q]]),ru_i);   //rc temporary is rb0 vector
              multiple_vec_scalar(&rb,u_i,calc_vec_vec_scalar_product(u_i,&rc));
              SUBT_VEC(((t_vec*)&rb),((t_vec*)&rc),((t_vec*)&rb));
              rrb=calc_vec_vec_scalar_product(&rb,&rb);
              if (rrb>ROTABLE_BOND_LIMIT*ROTABLE_BOND_LIMIT) //Process iff the rotable bond really exists
                {
                //Calculate centers distance and energy derivatives
                SUBT_VEC(((t_vec*)&rc),((t_vec*)&r_i[mol_i->anchors[anchors_i->list[_i]]->list[_p]]),((t_vec*)&(*r_j)[(*mol_j)->anchors[_j]->list[_q]]));
                r=calc_vec_vec_scalar_product(&rc,&rc);
                if (r<FF_CUTOFF)
                  {
                  r=sqrt(r);
                  e+=ffsys->fnbonded(r,&de,&dde,ffsys->e[mol_i->atoms->list[mol_i->anchors[anchors_i->list[_i]]->list[_p]]][(*mol_j)->atoms->list[(*mol_j)->anchors[_j]->list[_q]]],
                                                ffsys->R[mol_i->atoms->list[mol_i->anchors[anchors_i->list[_i]]->list[_p]]][(*mol_j)->atoms->list[(*mol_j)->anchors[_j]->list[_q]]],
                                                mol_i->charges[mol_i->anchors[anchors_i->list[_i]]->list[_p]],(*mol_j)->charges[(*mol_j)->anchors[_j]->list[_q]]);
                  de/=r;
                  dde=(dde-de)/r/r;
                  //Calculate relative angle
                  r=sqrt(rra*rrb);
                  csA=calc_vec_vec_scalar_product(&ra,&rb)/r;
                  //Calculate reference angle
                  dd+=dde*rra*rrb*(1.00-csA*csA)+de*r*csA;
                  }
                }
              }
          }
        r_j++;
        mol_j++;
        }
      }
    }
return dd;
}


//This function calculates independent distribution functions
double calc_independent_internal_dde(t_vec *r_i,t_vec *r_j,t_mol *mol_i,t_mol *mol_j,t_vec *u_i,t_vec *u_j,t_vec *ru_i,t_vec *ru_j,t_list *anchors_i,t_list *anchors_j,t_ffsys *ffsys)
{
double rra,rrb,r,csP,snP,csA,snA,csB,snB,e,de,dde,dd=0.00;
unsigned int _i,_j,_p,_q;
t_vec ra,rb,rc,rx,ry,h;

//Calc abscise axis, csP and snP
csP=calc_vec_vec_scalar_product(u_i,u_j);
//Real two angles interaction
snP=sqrt(1.00-csP*csP);
vec_vec_vmult(&rx,u_i,u_j);
vec_vec_vmult(&ry,u_i,&rx);
multiple_vec_scalar(&rx,&rx,1.00/sqrt(calc_vec_vec_scalar_product(&rx,&rx)));
multiple_vec_scalar(&ry,&ry,1.00/sqrt(calc_vec_vec_scalar_product(&ry,&ry)));
for (_i=0;_i<anchors_i->size;_i++)
  for (_p=0;_p<mol_i->anchors[anchors_i->list[_i]]->size;_p++)
    {
    //Calculate ra vector
    SUBT_VEC(((t_vec*)&rc),((t_vec*)&r_i[mol_i->anchors[anchors_i->list[_i]]->list[_p]]),ru_i);  //rc temporary is ra0 vector
    multiple_vec_scalar(&ra,u_i,calc_vec_vec_scalar_product(u_i,&rc));
    SUBT_VEC(((t_vec*)&ra),((t_vec*)&rc),((t_vec*)&ra));
    rra=calc_vec_vec_scalar_product(&ra,&ra);
    if (rra>ROTABLE_BOND_LIMIT*ROTABLE_BOND_LIMIT) //Do evaluation iff the rotation around alfa is really observed
      {
      //Calculate csA and snA
      rra=sqrt(rra);
      csA=calc_vec_vec_scalar_product(&ra,&rx)/rra;
      snA=calc_vec_vec_scalar_product(&ra,&ry)/rra;
      for (_j=0;_j<anchors_j->size;_j++)
        for (_q=0;_q<mol_j->anchors[anchors_j->list[_j]]->size;_q++)
          {
          //Calculate rb vector
          SUBT_VEC(((t_vec*)&rc),((t_vec*)&r_j[mol_j->anchors[anchors_j->list[_j]]->list[_q]]),ru_j);  //rc temporary is rb0 vector
          multiple_vec_scalar(&rb,u_j,calc_vec_vec_scalar_product(u_j,&rc));
          SUBT_VEC(((t_vec*)&rb),((t_vec*)&rc),((t_vec*)&rb));
          rrb=calc_vec_vec_scalar_product(&rb,&rb);
          if (rrb>ROTABLE_BOND_LIMIT*ROTABLE_BOND_LIMIT) //Do evaluation iff the rotation around beta is really observed
            {
            //Calculate csB and snB
            rrb=sqrt(rrb);
            csB=calc_vec_vec_scalar_product(&rb,&rx)/rrb;
            snB=sqrt(1.00-csB*csB);
            //Calculate rc vector, A-B distance and energy derivatives
            SUBT_VEC(((t_vec*)&rc),((t_vec*)&r_i[mol_i->anchors[anchors_i->list[_i]]->list[_p]]),((t_vec*)&r_j[mol_j->anchors[anchors_j->list[_j]]->list[_q]]));
            r=calc_vec_vec_scalar_product(&rc,&rc);
            if (r<FF_CUTOFF)
              {
              r=sqrt(r);
              e+=ffsys->fnbonded(r,&de,&dde,ffsys->e[mol_i->atoms->list[mol_i->anchors[anchors_i->list[_i]]->list[_p]]][mol_j->atoms->list[mol_j->anchors[anchors_j->list[_j]]->list[_q]]],
                                            ffsys->R[mol_i->atoms->list[mol_i->anchors[anchors_i->list[_i]]->list[_p]]][mol_j->atoms->list[mol_j->anchors[anchors_j->list[_j]]->list[_q]]],
                                            mol_i->charges[mol_i->anchors[anchors_i->list[_i]]->list[_p]],mol_j->charges[mol_j->anchors[anchors_j->list[_j]]->list[_q]]);
              de/=r;
              dde=(dde-de)/r/r;
              //Calculate vector of rotation centers
              SUMM_VEC(((t_vec*)&rc),((t_vec*)&rc),((t_vec*)&rb));
              SUBT_VEC(((t_vec*)&rc),((t_vec*)&rc),((t_vec*)&ra));
              //Calculate projections
              h.i=calc_vec_vec_scalar_product(&rc,&rx); //calc_vec_vec_scalar_product(&rc,&rc);  calc_vec_vec_scalar_product(&rb,&runits->u0[a_id]);
              h.j=calc_vec_vec_scalar_product(&rc,&ry);
              h.k=calc_vec_vec_scalar_product(&rc,u_i);
              //Calculate derivatives   sqrd(rra*csA-rrb*csB+h.i)+sqrd(rra*snA+rrb*snB*csP+h.j)+sqrd(-rrb*snB*snP+h.k)-r*r         h.i*h.i+h.j*h.j+h.k*h.k
              if (calc_vec_vec_scalar_product(&rb,u_i)<0.00)
                {
                if (fabs(sqrd(rra*csA-rrb*csB+h.i)+sqrd(rra*snA+rrb*snB*csP+h.j)+sqrd( rrb*snB*snP+h.k)-r*r)>0.00001)
                  printf("ERROR!!!\n");
                dd+=dde*rra*rrb*((rrb*csB-h.i)*snA+csA*(h.j+rrb*snB*csP))*((rra*csA+h.i)*snB+csB*((rra*snA+h.j)*csP+h.k*snP))+de*rra*rrb*(csA*csB*csP-snA*snB);
                }
              else
                {
                if (fabs(sqrd(rra*csA-rrb*csB+h.i)+sqrd(rra*snA-rrb*snB*csP+h.j)+sqrd(-rrb*snB*snP+h.k)-r*r)>0.00001)
                  printf("ERROR!!!\n");
                dd+=dde*rra*rrb*((rrb*csB-h.i)*snA+csA*(h.j-rrb*snB*csP))*((rra*csA+h.i)*snB-csB*((rra*snA+h.j)*csP+h.k*snP))-de*rra*rrb*(csA*csB*csP+snA*snB);
                }
              }
            }
          }
      }
    }
  
//Return derivative value
return dd;
}


//This function calculates dependent distribution functions for _i rotatoin against _j rotation.
double calc_dependent_internal_dde(t_vec **r_i,t_mol **mol_i,t_vec *u_i,t_vec *u_j,t_vec *ru_i,t_vec *ru_j,t_list *anchors_i,t_list *anchors_j,t_ffsys *ffsys)
{
double rra,rrb,r,csP,snP,csA,snA,csB,snB,e,de,dde,dd=0.00;
unsigned int _i,_j,_p,_q,nmols;
t_vec ra,rb,rc,rx,ry,h,**r_j;
t_mol **mol_j;

//Calc abscise axis, csP and snP
csP=calc_vec_vec_scalar_product(u_i,u_j);
//Real two angles interaction
snP=sqrt(1.00-csP*csP);
vec_vec_vmult(&rx,u_i,u_j);
vec_vec_vmult(&ry,u_i,&rx);
multiple_vec_scalar(&rx,&rx,1.00/sqrt(calc_vec_vec_scalar_product(&rx,&rx)));
multiple_vec_scalar(&ry,&ry,1.00/sqrt(calc_vec_vec_scalar_product(&ry,&ry)));
for (_i=0;_i<anchors_i->size;_i++)
  for (_p=0;_p<(*mol_i)->anchors[anchors_i->list[_i]]->size;_p++)
    {
    //Calculate ra vector
    SUBT_VEC(((t_vec*)&rc),((t_vec*)&(*r_i)[(*mol_i)->anchors[anchors_i->list[_i]]->list[_p]]),ru_i);  //rc temporary is ra0 vector
    multiple_vec_scalar(&ra,u_i,calc_vec_vec_scalar_product(u_i,&rc));
    SUBT_VEC(((t_vec*)&ra),((t_vec*)&rc),((t_vec*)&ra));
    rra=calc_vec_vec_scalar_product(&ra,&ra);
    if (rra>ROTABLE_BOND_LIMIT*ROTABLE_BOND_LIMIT) //Do evaluation iff the rotation around alfa is really observed
      {
      //Calculate csA and snA
      rra=sqrt(rra);
      csA=calc_vec_vec_scalar_product(&ra,&rx)/rra;
      snA=calc_vec_vec_scalar_product(&ra,&ry)/rra;
      r_j=ffsys->r;
      mol_j=ffsys->mols;
      nmols=ffsys->nmols;
      while (nmols--)
        {
        if (mol_j==mol_i)
          //Summ over current molecule
          for (_j=0;_j<anchors_j->size;_j++)
            for (_q=0;_q<(*mol_i)->anchors[anchors_j->list[_j]]->size;_q++)
              {
              //Calculate rb vector
              SUBT_VEC(((t_vec*)&rc),((t_vec*)&(*r_i)[(*mol_i)->anchors[anchors_j->list[_j]]->list[_q]]),ru_j);  //rc temporary is rb0 vector
              multiple_vec_scalar(&rb,u_j,calc_vec_vec_scalar_product(u_j,&rc));
              SUBT_VEC(((t_vec*)&rb),((t_vec*)&rc),((t_vec*)&rb));
              rrb=calc_vec_vec_scalar_product(&rb,&rb);
              if (rrb>ROTABLE_BOND_LIMIT*ROTABLE_BOND_LIMIT) //Do evaluation iff the rotation around beta is really observed
                {
                //Calculate csB and snB
                rrb=sqrt(rrb);
                csB=calc_vec_vec_scalar_product(&rb,&rx)/rrb;
                snB=sqrt(1.00-csB*csB);
                //Calculate rc vector, A-B distance and energy derivatives
                SUBT_VEC(((t_vec*)&rc),((t_vec*)&(*r_i)[(*mol_i)->anchors[anchors_i->list[_i]]->list[_p]]),((t_vec*)&(*r_i)[(*mol_i)->anchors[anchors_j->list[_j]]->list[_q]]));
                r=calc_vec_vec_scalar_product(&rc,&rc);
                if (r<FF_CUTOFF)
                  {
                  r=sqrt(r);
                  e+=ffsys->fnbonded(r,&de,&dde,ffsys->e[(*mol_i)->atoms->list[(*mol_i)->anchors[anchors_i->list[_i]]->list[_p]]][(*mol_i)->atoms->list[(*mol_i)->anchors[anchors_j->list[_j]]->list[_q]]],
                                                ffsys->R[(*mol_i)->atoms->list[(*mol_i)->anchors[anchors_i->list[_i]]->list[_p]]][(*mol_i)->atoms->list[(*mol_i)->anchors[anchors_j->list[_j]]->list[_q]]],
                                                (*mol_i)->charges[(*mol_i)->anchors[anchors_i->list[_i]]->list[_p]],(*mol_i)->charges[(*mol_i)->anchors[anchors_j->list[_j]]->list[_q]]);
                  de/=r;
                  dde=(dde-de)/r/r;
                  //Calculate vector of rotation centers
                  SUMM_VEC(((t_vec*)&rc),((t_vec*)&rc),((t_vec*)&rb));
                  SUBT_VEC(((t_vec*)&rc),((t_vec*)&rc),((t_vec*)&ra));
                  //Calculate projections
                  h.i=calc_vec_vec_scalar_product(&rc,&rx); //calc_vec_vec_scalar_product(&rc,&rc);  calc_vec_vec_scalar_product(&rb,&runits->u0[a_id]);
                  h.j=calc_vec_vec_scalar_product(&rc,&ry);
                  h.k=calc_vec_vec_scalar_product(&rc,u_i);
                  //Calculate derivatives   sqrd(rra*csA-rrb*csB+h.i)+sqrd(rra*snA+rrb*snB*csP+h.j)+sqrd(-rrb*snB*snP+h.k)-r*r         h.i*h.i+h.j*h.j+h.k*h.k
                  if (calc_vec_vec_scalar_product(&rb,u_i)<0.00)
                    {
                    if (fabs(sqrd(rra*csA-rrb*csB+h.i)+sqrd(rra*snA+rrb*snB*csP+h.j)+sqrd( rrb*snB*snP+h.k)-r*r)>0.00001)
                      printf("ERROR!!!\n");
                    dd+=dde*rra*rrb*((rrb*csB-h.i)*snA+csA*(h.j+rrb*snB*csP))*((rra*csA+h.i)*snB+csB*((rra*snA+h.j)*csP+h.k*snP))+de*rra*rrb*(csA*csB*csP-snA*snB);
                    }
                  else
                    {
                    if (fabs(sqrd(rra*csA-rrb*csB+h.i)+sqrd(rra*snA-rrb*snB*csP+h.j)+sqrd(-rrb*snB*snP+h.k)-r*r)>0.00001)
                      printf("ERROR!!!\n");
                    dd+=dde*rra*rrb*((rrb*csB-h.i)*snA+csA*(h.j-rrb*snB*csP))*((rra*csA+h.i)*snB-csB*((rra*snA+h.j)*csP+h.k*snP))-de*rra*rrb*(csA*csB*csP+snA*snB);
                    }
                  }
                }
              }
        else  //Summ over all other molecules
          for (_j=0;_j<(*mol_j)->nanchors;_j++)
            for (_q=0;_q<(*mol_j)->anchors[_j]->size;_q++)
              {
              //Calculate rb vector
              SUBT_VEC(((t_vec*)&rc),((t_vec*)&(*r_j)[(*mol_j)->anchors[_j]->list[_q]]),ru_j);  //rc temporary is rb0 vector
              multiple_vec_scalar(&rb,u_j,calc_vec_vec_scalar_product(u_j,&rc));
              SUBT_VEC(((t_vec*)&rb),((t_vec*)&rc),((t_vec*)&rb));
              rrb=calc_vec_vec_scalar_product(&rb,&rb);
              if (rrb>ROTABLE_BOND_LIMIT*ROTABLE_BOND_LIMIT) //Do evaluation iff the rotation around beta is really observed
                {
                //Calculate csB and snB
                rrb=sqrt(rrb);
                csB=calc_vec_vec_scalar_product(&rb,&rx)/rrb;
                snB=sqrt(1.00-csB*csB);
                //Calculate rc vector, A-B distance and energy derivatives
                SUBT_VEC(((t_vec*)&rc),((t_vec*)&(*r_i)[(*mol_i)->anchors[anchors_i->list[_i]]->list[_p]]),((t_vec*)&(*r_j)[(*mol_j)->anchors[_j]->list[_q]]));
                r=calc_vec_vec_scalar_product(&rc,&rc);
                if (r<FF_CUTOFF)
                  {
                  r=sqrt(r);
                  e+=ffsys->fnbonded(r,&de,&dde,ffsys->e[(*mol_i)->atoms->list[(*mol_i)->anchors[anchors_i->list[_i]]->list[_p]]][(*mol_j)->atoms->list[(*mol_j)->anchors[_j]->list[_q]]],
                                                ffsys->R[(*mol_i)->atoms->list[(*mol_i)->anchors[anchors_i->list[_i]]->list[_p]]][(*mol_j)->atoms->list[(*mol_j)->anchors[_j]->list[_q]]],
                                                (*mol_i)->charges[(*mol_i)->anchors[anchors_i->list[_i]]->list[_p]],(*mol_j)->charges[(*mol_j)->anchors[_j]->list[_q]]);
                  de/=r;
                  dde=(dde-de)/r/r;
                  //Calculate vector of rotation centers
                  SUMM_VEC(((t_vec*)&rc),((t_vec*)&rc),((t_vec*)&rb));
                  SUBT_VEC(((t_vec*)&rc),((t_vec*)&rc),((t_vec*)&ra));
                  //Calculate projections
                  h.i=calc_vec_vec_scalar_product(&rc,&rx); //calc_vec_vec_scalar_product(&rc,&rc);  calc_vec_vec_scalar_product(&rb,&runits->u0[a_id]);
                  h.j=calc_vec_vec_scalar_product(&rc,&ry);
                  h.k=calc_vec_vec_scalar_product(&rc,u_i);
                  //Calculate derivatives   sqrd(rra*csA-rrb*csB+h.i)+sqrd(rra*snA+rrb*snB*csP+h.j)+sqrd(-rrb*snB*snP+h.k)-r*r         h.i*h.i+h.j*h.j+h.k*h.k
                  if (calc_vec_vec_scalar_product(&rb,u_i)<0.00)
                    {
                    if (fabs(sqrd(rra*csA-rrb*csB+h.i)+sqrd(rra*snA+rrb*snB*csP+h.j)+sqrd( rrb*snB*snP+h.k)-r*r)>0.00001)
                      printf("ERROR!!!\n");
                    dd+=dde*rra*rrb*((rrb*csB-h.i)*snA+csA*(h.j+rrb*snB*csP))*((rra*csA+h.i)*snB+csB*((rra*snA+h.j)*csP+h.k*snP))+de*rra*rrb*(csA*csB*csP-snA*snB);
                    }
                  else
                    {
                    if (fabs(sqrd(rra*csA-rrb*csB+h.i)+sqrd(rra*snA-rrb*snB*csP+h.j)+sqrd(-rrb*snB*snP+h.k)-r*r)>0.00001)
                      printf("ERROR!!!\n");
                    dd+=dde*rra*rrb*((rrb*csB-h.i)*snA+csA*(h.j-rrb*snB*csP))*((rra*csA+h.i)*snB-csB*((rra*snA+h.j)*csP+h.k*snP))-de*rra*rrb*(csA*csB*csP+snA*snB);
                    }
                  }
                }
              }
        r_j++;
        mol_j++;
        }
      }
  }
//Return derivative value
return dd;
}


//This function calculates internal distribution functions
double calc_internal_distribution(t_runits **runits,t_ffsys *ffsys)
{
t_matrix *z=0x0;
unsigned int _i,_j,size_i,size_j,nmols_i,nmols_j;
t_mol **mol_i,**mol_j;
t_runits **runit_i,**runit_j;
t_vec **r_i,**r_j;
double _det=1.00;

//Stage 1. Gather data about system complexity
size_i=0;
runit_i=runits;
nmols_i=ffsys->nmols;
while(nmols_i--)
  {
  size_i+=(*runit_i)->size;
  runit_i++;
  }
if (size_i)
  {
  //Stage 2. Alloc memory for matrix
  if (!(z=(t_matrix*)alloc_matrix(size_i<<1,size_i<<1))) return FALSE;
  //Stage 3. Fill z
  set_matrix(z,FALSE);
  nmols_i=ffsys->nmols;
  r_i=&ffsys->r[ffsys->nmols];
  mol_i=&ffsys->mols[ffsys->nmols];
  while(nmols_i--)
    {
    r_i--;
    mol_i--;
    runit_i--;
    size_i-=(*runit_i)->size;
    _i=(*runit_i)->size;
//    show_runits(*runit_i,*mol_i);
    while(_i--)
      {
      //Pass inter molecular distributions
      runit_j=runit_i;
      z->matrix[size_i+(_i<<1)+0x0][size_i+(_i<<1)+0x0]=calc_diagonal_internal_dde(*r_i,&(*runit_i)->u0[_i][0x0],&(*runit_i)->ru0[_i],(*runit_i)->anchors[_i][0x0],(*runit_i)->anchors[_i][0x1],*mol_i,ffsys);
      z->matrix[size_i+(_i<<1)+0x1][size_i+(_i<<1)+0x1]=calc_diagonal_internal_dde(*r_i,&(*runit_i)->u0[_i][0x1],&(*runit_i)->ru0[_i],(*runit_i)->anchors[_i][0x1],(*runit_i)->anchors[_i][0x0],*mol_i,ffsys);
      _j=_i;
      while(_j--)
        if ((find_in_list(FALSE,*(*runit_i)->anchors[_j][0x0]->list,(*runit_i)->anchors[_i][0x0])))
          {
          if ((find_in_list(FALSE,*(*runit_i)->anchors[_i][0x1]->list,(*runit_i)->anchors[_j][0x0])))
            {
            //        |     i1           |        i0     |
            //        |                      j0  |  j1   |
            z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x0]=calc_independent_internal_dde(*r_i,*r_i,*mol_i,*mol_i,&(*runit_i)->u0[_i][0x1],&(*runit_i)->u0[_j][0x1],&(*runit_i)->ru0[_i],&(*runit_i)->ru0[_j],(*runit_i)->anchors[_i][0x1],(*runit_i)->anchors[_j][0x1],ffsys); //SYMETRICAL
            z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x1]=  calc_dependent_internal_dde(r_i,mol_i,&(*runit_i)->u0[_j][0x1],&(*runit_i)->u0[_i][0x0],&(*runit_i)->ru0[_j],&(*runit_i)->ru0[_i],(*runit_i)->anchors[_j][0x1],(*runit_i)->anchors[_i][0x1],ffsys); //ANTI-
            z->matrix[((size_i+_i)<<1)+0x1][((size_i+_j)<<1)+0x0]=  calc_dependent_internal_dde(r_i,mol_i,&(*runit_i)->u0[_i][0x1],&(*runit_i)->u0[_j][0x0],&(*runit_i)->ru0[_i],&(*runit_i)->ru0[_j],(*runit_i)->anchors[_i][0x1],(*runit_i)->anchors[_j][0x1],ffsys); //   SYMETRICAL
            z->matrix[((size_i+_i)<<1)+0x1][((size_i+_j)<<1)+0x1]=z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x0];
            }
          else
            {
            //        |     i1           |        i0     |
            //        |                      j1  |  j0   |
            z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x1]=calc_independent_internal_dde(*r_i,*r_i,*mol_i,*mol_i,&(*runit_i)->u0[_i][0x1],&(*runit_i)->u0[_j][0x0],&(*runit_i)->ru0[_i],&(*runit_i)->ru0[_j],(*runit_i)->anchors[_i][0x1],(*runit_i)->anchors[_j][0x0],ffsys); //SYMETRICAL
            z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x0]=  calc_dependent_internal_dde(r_i,mol_i,&(*runit_i)->u0[_j][0x0],&(*runit_i)->u0[_i][0x0],&(*runit_i)->ru0[_j],&(*runit_i)->ru0[_i],(*runit_i)->anchors[_j][0x0],(*runit_i)->anchors[_i][0x1],ffsys); //ANTI-
            z->matrix[((size_i+_i)<<1)+0x1][((size_i+_j)<<1)+0x1]=  calc_dependent_internal_dde(r_i,mol_i,&(*runit_i)->u0[_i][0x1],&(*runit_i)->u0[_j][0x1],&(*runit_i)->ru0[_i],&(*runit_i)->ru0[_j],(*runit_i)->anchors[_i][0x1],(*runit_i)->anchors[_j][0x0],ffsys); //   SYMETRICAL
            z->matrix[((size_i+_i)<<1)+0x1][((size_i+_j)<<1)+0x0]=z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x1];
            }
          }
        else if ((find_in_list(FALSE,*(*runit_i)->anchors[_i][0x0]->list,(*runit_i)->anchors[_j][0x0])))
               {
               //        |     i0           |        i1     |
               //        |                      j0  |  j1   |
               z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x1]=calc_independent_internal_dde(*r_i,*r_i,*mol_i,*mol_i,&(*runit_i)->u0[_i][0x0],&(*runit_i)->u0[_j][0x1],&(*runit_i)->ru0[_i],&(*runit_i)->ru0[_j],(*runit_i)->anchors[_i][0x0],(*runit_i)->anchors[_j][0x1],ffsys); //SYMETRICAL
               z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x0]=  calc_dependent_internal_dde(r_i,mol_i,&(*runit_i)->u0[_i][0x0],&(*runit_i)->u0[_j][0x0],&(*runit_i)->ru0[_i],&(*runit_i)->ru0[_j],(*runit_i)->anchors[_i][0x0],(*runit_i)->anchors[_j][0x1],ffsys); //ANTI-
               z->matrix[((size_i+_i)<<1)+0x1][((size_i+_j)<<1)+0x1]=  calc_dependent_internal_dde(r_i,mol_i,&(*runit_i)->u0[_j][0x1],&(*runit_i)->u0[_i][0x1],&(*runit_i)->ru0[_j],&(*runit_i)->ru0[_i],(*runit_i)->anchors[_j][0x1],(*runit_i)->anchors[_i][0x0],ffsys); //   SYMETRICAL
               z->matrix[((size_i+_i)<<1)+0x1][((size_i+_j)<<1)+0x0]=z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x1];
               }
             else
               {
               //        |     i0           |        i1     |
               //        |                      j1  |  j0   |
               z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x0]=calc_independent_internal_dde(*r_i,*r_i,*mol_i,*mol_i,&(*runit_i)->u0[_i][0x0],&(*runit_i)->u0[_j][0x0],&(*runit_i)->ru0[_i],&(*runit_i)->ru0[_j],(*runit_i)->anchors[_i][0x0],(*runit_i)->anchors[_j][0x0],ffsys); //SYMETRICAL
               z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x1]=  calc_dependent_internal_dde(r_i,mol_i,&(*runit_i)->u0[_i][0x0],&(*runit_i)->u0[_j][0x1],&(*runit_i)->ru0[_i],&(*runit_i)->ru0[_j],(*runit_i)->anchors[_i][0x0],(*runit_i)->anchors[_j][0x0],ffsys); //ANTI-
               z->matrix[((size_i+_i)<<1)+0x1][((size_i+_j)<<1)+0x0]=  calc_dependent_internal_dde(r_i,mol_i,&(*runit_i)->u0[_j][0x0],&(*runit_i)->u0[_i][0x1],&(*runit_i)->ru0[_j],&(*runit_i)->ru0[_i],(*runit_i)->anchors[_j][0x0],(*runit_i)->anchors[_i][0x0],ffsys); //   SYMETRICAL
               z->matrix[((size_i+_i)<<1)+0x1][((size_i+_j)<<1)+0x1]=z->matrix[((size_i+_i)<<1)+0x0][((size_i+_j)<<1)+0x0];
               }
      //Stage 4. Current mol to rest molecules runits
      nmols_j=nmols_i;
      r_j=r_i;
      mol_j=mol_i;
      size_j=size_i;
      runit_j=runit_i;
      while(nmols_j--)
        {
        r_j--;
        mol_j--;
        runit_j--;
        size_j-=(*runit_j)->size;
        _j=(*runit_j)->size;
        while(_j--)
          {
          z->matrix[((size_i+_i)<<1)+0x0][((size_j+_j)<<1)+0x0]=calc_independent_internal_dde(*r_i,*r_j,*mol_i,*mol_j,&(*runit_i)->u0[_i][0x0],&(*runit_j)->u0[_j][0x0],&(*runit_i)->ru0[_i],&(*runit_j)->ru0[_j],(*runit_i)->anchors[_i][0x0],(*runit_j)->anchors[_j][0x1],ffsys);
          z->matrix[((size_i+_i)<<1)+0x0][((size_j+_j)<<1)+0x1]=calc_independent_internal_dde(*r_i,*r_j,*mol_i,*mol_j,&(*runit_i)->u0[_i][0x0],&(*runit_j)->u0[_j][0x1],&(*runit_i)->ru0[_i],&(*runit_j)->ru0[_j],(*runit_i)->anchors[_i][0x0],(*runit_j)->anchors[_j][0x0],ffsys);
          z->matrix[((size_i+_i)<<1)+0x1][((size_j+_j)<<1)+0x0]=calc_independent_internal_dde(*r_i,*r_j,*mol_i,*mol_j,&(*runit_i)->u0[_i][0x1],&(*runit_j)->u0[_j][0x0],&(*runit_i)->ru0[_i],&(*runit_j)->ru0[_j],(*runit_i)->anchors[_i][0x1],(*runit_j)->anchors[_j][0x1],ffsys);
          z->matrix[((size_i+_i)<<1)+0x1][((size_j+_j)<<1)+0x1]=calc_independent_internal_dde(*r_i,*r_j,*mol_i,*mol_j,&(*runit_i)->u0[_i][0x1],&(*runit_j)->u0[_j][0x1],&(*runit_i)->ru0[_i],&(*runit_j)->ru0[_j],(*runit_i)->anchors[_i][0x1],(*runit_j)->anchors[_j][0x0],ffsys);
          }
        }
      }
    }
  //Stage 5. Calc P(B)
  _i=z->x;
  while(_i--)
    {
    _j=_i;
    while(_j--)
      z->matrix[_j][_i]=z->matrix[_i][_j]; //Symmeterize matrix
    }
  //show_matrix(z);
  diagonalization(0xFF,z);
  _i=z->x;
  while(_i--)
    if (z->matrix[_i][_i]>0.00)
      _det*=z->matrix[_i][_i];
    else  if (z->matrix[_i][_i]<0.00)
           _det/=-z->matrix[_i][_i];
  //show_matrix(z);
  //LU_decomposition (z->matrix,z->x,p_order,twoness,temp);
  //_f=LU_determinant(z->matrix,z->x,p_order,twoness[0x0]);
  //Stage 5. Free memory;
  free(z);
  }

return 0.5*log(_det);
}

*/

//This function calculates mass of given anchor fragment
inline void calc_anchors_cm(t_vec *cm,t_vec *r,unsigned int anchor_id,t_mol *mol,t_top *top)
{
register unsigned int _i;
register double _m;
_m=cm->i=cm->j=cm->k=0.;
_i=mol->anchors->list[anchor_id].size;
while (_i--)
  {
  cm->i+=top->ff_a[mol->atoms->list[mol->anchors->list[anchor_id].list[_i]]].mass*r[mol->anchors->list[anchor_id].list[_i]].i;
  cm->j+=top->ff_a[mol->atoms->list[mol->anchors->list[anchor_id].list[_i]]].mass*r[mol->anchors->list[anchor_id].list[_i]].j;
  cm->k+=top->ff_a[mol->atoms->list[mol->anchors->list[anchor_id].list[_i]]].mass*r[mol->anchors->list[anchor_id].list[_i]].k;
  _m+=top->ff_a[mol->atoms->list[mol->anchors->list[anchor_id].list[_i]]].mass;
  }
cm->i/=_m, cm->j/=_m, cm->k/=_m;
}

//This function move system to the molecule cm
inline void calc_mols_cm(t_vec *cm,t_vec *r_vecs,t_mol *mol,t_top *top)
{
register unsigned int _i;
register double m;

m=cm->i=cm->j=cm->k=0.00;
_i=mol->atoms->size;
while(_i--)
  {
  cm->i+=top->ff_a[mol->atoms->list[_i]].mass*r_vecs[_i].i;
  cm->j+=top->ff_a[mol->atoms->list[_i]].mass*r_vecs[_i].j;
  cm->k+=top->ff_a[mol->atoms->list[_i]].mass*r_vecs[_i].k;
  m+=top->ff_a[mol->atoms->list[_i]].mass;
  }
cm->i/=m;
cm->j/=m;
cm->k/=m;
}

//This function calculates second derivaties matrix of distance square over six freedom dimensions
//Note the system have to be located in the center of mass of it
double calc_3D_RT_Hessian_yff1(double Z[0x6][0x6],unsigned int mol_id,t_ffsys *ffsys,t_top *top)
{
unsigned int _i,_j,_k,_l,mols,*list_i,*list_j;
double e,de,dde,r,rr,rrrrrr,A,B,Q,d[0x6],*q_i,*q_j;
t_vec *r_i,*r_j,cm;

//Init matrix
e=0.00;
Z[0][0]=Z[0][1]=Z[0][2]=Z[0][3]=Z[0][4]=Z[0][5]=0.00;
Z[1][0]=Z[1][1]=Z[1][2]=Z[1][3]=Z[1][4]=Z[1][5]=0.00;
Z[2][0]=Z[2][1]=Z[2][2]=Z[2][3]=Z[2][4]=Z[2][5]=0.00;
Z[3][0]=Z[3][1]=Z[3][2]=Z[3][3]=Z[3][4]=Z[3][5]=0.00;
Z[4][0]=Z[4][1]=Z[4][2]=Z[4][3]=Z[4][4]=Z[4][5]=0.00;
Z[5][0]=Z[5][1]=Z[5][2]=Z[5][3]=Z[5][4]=Z[5][5]=0.00;
calc_mols_cm(&cm,&ffsys->r[ffsys->nr[mol_id]],ffsys->mols[mol_id],top);
_i=ffsys->natoms; while (_i--) { ffsys->r[_i].i-=cm.i, ffsys->r[_i].j-=cm.j, ffsys->r[_i].k-=cm.k; }

//Loop over all atoms of mol_id
_i=ffsys->mols[mol_id]->atoms->size;
r_i=ffsys->r+ffsys->nr[mol_id]+_i;
q_i=ffsys->mols[mol_id]->charges+_i;
list_i=ffsys->mols[mol_id]->atoms->list+_i;
while (_i--)
  {
  r_i--;
  q_i--;
  list_i--;
  mols=ffsys->nmols; //loop over molecules in system
  while (mols--)
  //while (--mols)
    if (mols!=mol_id)  //skip self eenrgy
      { //loop over all atoms of mols
      _j=ffsys->mols[mols]->atoms->size;
      r_j=ffsys->r+ffsys->nr[mols]+_j;
      q_j=ffsys->mols[mols]->charges+_j;
      list_j=ffsys->mols[mols]->atoms->list+_j;
      while (_j--)
        {
        r_j--;
        q_j--;
        list_j--;
        rr=sqrd(r_i->i-r_j->i)+sqrd(r_i->j-r_j->j)+sqrd(r_i->k-r_j->k)+SMALL_DISTORTION;
     //   if (rr<FF_CUTOFF2)
          {
          //Calculate second energy derivative
          r=sqrt(rr);
          rrrrrr=rr*rr*rr;
          A=top->A[*list_i][*list_j]/rrrrrr/rrrrrr;
          B=top->B[*list_i][*list_j]/rrrrrr;
          Q=COULOMB_K**q_i**q_j/r;
          e+=A-B+Q; 
          de=(-12.00*A+6.00*B-Q)/r;
          dde=(156.00*A-42.00*B+2.00*Q)/rr;
          //Prepare derivatives
          de/=r;
          dde=(dde-de)/rr;
          //Stage 3. Update hessian matrix wiht k*D
          //row0
          Z[0x0][0x0]+=de;
          Z[0x0][0x4]-=de*r_j->k;
          Z[0x0][0x5]+=de*r_j->j;
          //row1
          Z[0x1][0x1]+=de;
          Z[0x1][0x3]+=de*r_j->k;
          Z[0x1][0x5]-=de*r_j->i;
          //row2
          Z[0x2][0x2]+=de;
          Z[0x2][0x3]-=de*r_j->j;
          Z[0x2][0x4]+=de*r_j->i;
          //row 3
          Z[0x3][0x1]+=de*r_j->k;
          Z[0x3][0x2]-=de*r_j->j;
          Z[0x3][0x3]+=de*(r_i->j*r_j->j+r_i->k*r_j->k);
          Z[0x3][0x4]-=de*r_i->i*r_j->j;
          Z[0x3][0x5]-=de*r_i->i*r_j->k;
          //row 4
          Z[0x4][0x0]-=de*r_j->k;
          Z[0x4][0x2]+=de*r_j->i;
          Z[0x4][0x3]-=de*r_i->j*r_j->i;
          Z[0x4][0x4]+=de*(r_i->i*r_j->i+r_i->k*r_j->k);
          Z[0x4][0x5]-=de*r_i->j*r_j->k;
          //row 5
          Z[0x5][0x0]+=de*r_j->j;
          Z[0x5][0x1]-=de*r_j->i;
          Z[0x5][0x3]-=de*r_i->k*r_j->i;
          Z[0x5][0x4]-=de*r_i->k*r_j->j;
          Z[0x5][0x5]+=de*(r_i->i*r_j->i+r_i->j*r_j->j);
          //Stage 5. Calculate derivative vector
          d[0x0]=r_i->i-r_j->i;
          d[0x1]=r_i->j-r_j->j;
          d[0x2]=r_i->k-r_j->k;
          d[0x3]=r_i->j*r_j->k-r_i->k*r_j->j;
          d[0x4]=r_i->k*r_j->i-r_i->i*r_j->k;
          d[0x5]=r_i->i*r_j->j-r_i->j*r_j->i;
          //Stage 6. Update hessian matrix with K*dT.d
          for (_k=0x6;_k--;)
            for (_l=0x6;_l--;)
              Z[_k][_l]+=dde*d[_l]*d[_k];
          }
        }
      }
  }

//Show Z
//Z[0][3]=Z[0][4]=Z[0][5]=Z[1][3]=Z[1][4]=Z[1][5]=Z[2][3]=Z[2][4]=Z[2][5]=0.;
//Z[3][0]=Z[4][0]=Z[5][0]=Z[3][1]=Z[4][1]=Z[5][1]=Z[3][2]=Z[4][2]=Z[5][2]=0.;
//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",Z[0][0],Z[0][1],Z[0][2],Z[0][3],Z[0][4],Z[0][5]);
//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",Z[1][0],Z[1][1],Z[1][2],Z[1][3],Z[1][4],Z[1][5]);
//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",Z[2][0],Z[2][1],Z[2][2],Z[2][3],Z[2][4],Z[2][5]);
//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",Z[3][0],Z[3][1],Z[3][2],Z[3][3],Z[3][4],Z[3][5]);
//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",Z[4][0],Z[4][1],Z[4][2],Z[4][3],Z[4][4],Z[4][5]);
//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",Z[5][0],Z[5][1],Z[5][2],Z[5][3],Z[5][4],Z[5][5]);
//Move system back
_i=ffsys->natoms; while (_i--) { ffsys->r[_i].i+=cm.i, ffsys->r[_i].j+=cm.j, ffsys->r[_i].k+=cm.k; }
return e;
}


//This function calculates external distribution functions for 6D dimensions
double calc_external_distribution(double *e,double Z[0x6][0x6],unsigned int mol_id,t_ffsys *ffsys,t_top *top)
{
double _det;
//Stage 1. Calculate Hessian matrix
*e=calc_3D_RT_Hessian_yff1(Z,mol_id,ffsys,top); 
//Stage 2. Calculate determinant of Hessian matrix
_det=calc_det6x6(Z[0][0],Z[0][1],Z[0][2],Z[0][3],Z[0][4],Z[0][5],
                 Z[1][0],Z[1][1],Z[1][2],Z[1][3],Z[1][4],Z[1][5],
                 Z[2][0],Z[2][1],Z[2][2],Z[2][3],Z[2][4],Z[2][5],
                 Z[3][0],Z[3][1],Z[3][2],Z[3][3],Z[3][4],Z[3][5],
                 Z[4][0],Z[4][1],Z[4][2],Z[4][3],Z[4][4],Z[4][5],
                 Z[5][0],Z[5][1],Z[5][2],Z[5][3],Z[5][4],Z[5][5] );
return log(fabs(_det));
//return log(cube(THERMODYNAMIC_R*THERMODYNAMIC_T)/MOLAR_BOX_VOLUME)-0.25*log(_det);
}




