//This module contain various routines connected with graphs
#include <stdlib.h>
#include <string.h>
#include "y_system.h"
#include "y_txt.h"
#include "y_list.h"
#include "y_graph.h"


//This function try to find an edge of two vertices suggesting it can be either A-B or B-A. It returns (unsigned int)-1 if there is no edge
inline unsigned int find_edge_unordered(register unsigned int _a,register unsigned int _b,register unsigned int _n,register t_edge *_e)
{
while (_n--) if ( ( (_e[_n].vertice[0]==_a)&&(_e[_n].vertice[1]==_b) )||( (_e[_n].vertice[0]==_b)&&(_e[_n].vertice[1]==_a) ) ) break; 
return _n;
} 

//This function defines neighbors
//Note. it doesn't check preallocated memory [should be sizeof(t_clist)+sizeof(t_list)*nvertice+sizeof(unsigned int)*0x2*nedges]
t_clist *define_neighbors(t_clist *neighbors,unsigned int nvertices,unsigned int nedges,t_edge *edges)
{
register unsigned int _i;
t_clist *_neighbors;
if (!nvertices) { LABEL_DATA_CONSISTMENT_ERROR: ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
//Allocate memory if needed
if ( (neighbors)) _neighbors=neighbors;
else if (!(_neighbors=(t_clist*)malloc(sizeof(t_clist)+sizeof(t_list)*nvertices+sizeof(unsigned int)*0x2*nedges))) { ylib_errno=YERROR_MEMORY; return FALSE; }
//Map memory
_neighbors->size=nvertices;
_neighbors->_size=nedges*2;
_neighbors->list=(void*)_neighbors+sizeof(t_clist);
_neighbors->_items=(void*)_neighbors+sizeof(t_clist)+sizeof(t_list)*nvertices;
_i=nvertices; while(_i--) _neighbors->list[_i].size=0;
_i=nedges;    
while(_i--) 
  if ( (edges[_i].vertice[0]<nvertices)&&(edges[_i].vertice[1]<nvertices) )
    { _neighbors->list[edges[_i].vertice[0]].size++, _neighbors->list[edges[_i].vertice[1]].size++; }
  else { if (!(neighbors)) free(_neighbors); goto LABEL_DATA_CONSISTMENT_ERROR; }
for (_neighbors->list[0].list=_neighbors->_items, _i=1; _i<nvertices; _i++) 
  { _neighbors->list[_i].list=_neighbors->list[_i-1].list+_neighbors->list[_i-1].size, _neighbors->list[_i-1].size=0; }
_neighbors->list[_i-1].size=0;
_i=nedges;
while (_i--)
  {
  _neighbors->list[edges[_i].vertice[0]].list[_neighbors->list[edges[_i].vertice[0]].size++]=edges[_i].vertice[1];
  _neighbors->list[edges[_i].vertice[1]].list[_neighbors->list[edges[_i].vertice[1]].size++]=edges[_i].vertice[0];
  }
return _neighbors;
}

//This function enumerate vertices from start till end with neighbors. It returns amount of steps used (or (unsigned int)-1 to denote infinity/no path situation).
//NOTE. If you wona total enumeration use end==(unsigned int)-1.
unsigned int enumerate_vertices_with_neighbors(unsigned int start,unsigned int end,t_clist *neighbors,unsigned int *step,unsigned int *stack)
{
register unsigned int _i,_j,_k;
unsigned int curr_step, curr_size, next_size;

_i=neighbors->size;
while (_i--) step[_i]=(unsigned int)-1;
///Init massives
_i=step[stack[0x0]=start]=0x0;
next_size=curr_step=curr_size=0x1;
if (!(_j=neighbors->list[start].size)) return FALSE; //Single item submitted
while (_j--)
  step[stack[next_size++]=neighbors->list[start].list[_j]]=curr_step;
//Choose search strategy
if (end!=(unsigned int)-1)
  {
  while (++_i!=next_size)
    {
    if (_i==curr_size)
      {
      curr_size=next_size;
      curr_step++;
      }
    _j=neighbors->list[stack[_i]].size;
    while (_j--)
      if (step[_k=neighbors->list[stack[_i]].list[_j]]==(unsigned int)-1)
        {
        step[_k]=curr_step;
        stack[next_size++]=_k;
        if (_k==end) return curr_step; //difference in starategies
        }
    }
  return (unsigned int)-1; //No joining betweenb start and end so edge lenght = inf.
  }
else
  {
  do{
    if (++_i==curr_size)
      {
      curr_size=next_size;
      curr_step++;
      }
    _j=neighbors->list[stack[_i]].size;
    while (_j--)
      if (step[_k=neighbors->list[stack[_i]].list[_j]]==(unsigned int)-1)
        {
        step[_k]=curr_step;
        stack[next_size++]=_k;
        }
    }while ( (next_size!=neighbors->size)&&(next_size!=_i) );
  if (next_size==neighbors->size) return curr_step;
  else                            return (unsigned int)-1; //Fragment detected
  }
}

//This routine calculates intra-graph distances matrix for subset dist_v of vertices in graph size_v, neighbors
char calculate_intragraph_distance_matrix(unsigned int **d,unsigned int size_vd,unsigned int *vd,t_clist *neighbors)
{
register unsigned int _i, _j; 
unsigned int *step;
if (!(step=(unsigned int*)malloc(sizeof(unsigned int)*neighbors->size*0x2))) { ylib_errno=YERROR_MEMORY; return FALSE; } 
_i=size_vd;
while (--_i)
  {
  enumerate_vertices_with_neighbors(vd[_i],(unsigned int)-1,neighbors,step,&step[neighbors->size]);
  d[_i][_i]=0, _j=_i; while (_j--) d[_i][_j]=d[_j][_i]=step[vd[_j]];
  }
d[0][0]=0;
free(step);
return TRUE;
}

//This function do enumeration of atoms with neighbors, builds quasy-rtree structure of neighbors and maps cycles edges stack[...]==size
//Note. it return -1 if graph is not solid
unsigned int enumerate_cycles_with_neighbors (unsigned int root,unsigned int size,unsigned int *enumerate,t_clist *neighbors,unsigned int *cycles,unsigned int *stack)
{
register unsigned int _i, _j, _k,_l;
_i=size; while (_i--) { enumerate[_i]=(unsigned int)-1, cycles[_i]=0; }
*stack=root, enumerate[root]=0, _i=0, _j=1, _k=neighbors->list[root].size; while (_k--) { stack[_j]=neighbors->list[root].list[_k], enumerate[stack[_j]]=1, _j++; }
while (++_i!=_j)
  {
  _k=neighbors->list[stack[_i]].size;
  while (--_k) 
    if (enumerate[neighbors->list[stack[_i]].list[_k]]==(unsigned int)-1) { ENUMERATE: stack[_j]=neighbors->list[stack[_i]].list[_k], enumerate[stack[_j]]=enumerate[stack[_i]]+1, _j++; }
    else 
      {
      if (enumerate[*neighbors->list[stack[_i]].list]==(unsigned int)-1) //root detected
        {
        _l=*neighbors->list[stack[_i]].list, *neighbors->list[stack[_i]].list=neighbors->list[stack[_i]].list[_k], neighbors->list[stack[_i]].list[_k]=_l; //Swap
        goto ENUMERATE;
        }
      else //Cycle detected
        {
        //Choose the smallest as the root
        if (enumerate[neighbors->list[stack[_i]].list[_k]]<enumerate[*neighbors->list[stack[_i]].list])
          { _l=neighbors->list[stack[_i]].list[_k], neighbors->list[stack[_i]].list[_k]=*neighbors->list[stack[_i]].list, *neighbors->list[stack[_i]].list=_l; } //Swap
        //Find element in upper massive and swap
        cycles[neighbors->list[stack[_i]].list[_k]]++, _l=--neighbors->list[neighbors->list[stack[_i]].list[_k]].size;
        while (_l--)
          if (neighbors->list[neighbors->list[stack[_i]].list[_k]].list[_l]==stack[_i])
            {
            neighbors->list[neighbors->list[stack[_i]].list[_k]].list[_l]=neighbors->list[neighbors->list[stack[_i]].list[_k]].list[neighbors->list[neighbors->list[stack[_i]].list[_k]].size], 
            neighbors->list[neighbors->list[stack[_i]].list[_k]].list[neighbors->list[neighbors->list[stack[_i]].list[_k]].size]=stack[_i]; 
            break;
            }   
        //Move inpoint position
        cycles[stack[_i]]++;
        if (--neighbors->list[stack[_i]].size!=_k) 
          { //Swap
          _l=neighbors->list[stack[_i]].list[neighbors->list[stack[_i]].size], 
          neighbors->list[stack[_i]].list[neighbors->list[stack[_i]].size]=neighbors->list[stack[_i]].list[_k],
          neighbors->list[stack[_i]].list[_k]=_l;
          } 
        } 
      }
  if (enumerate[*neighbors->list[stack[_i]].list]==(unsigned int)-1) return (unsigned int)-1; //Graph is not solid
  }
if (_i!=size) return (unsigned int)-1; //Check if graph is solid
return enumerate[stack[_i-1]];
}

//This function cuts onevalent vertices
unsigned int cut_onevalent_with_neighbors(t_list *ids,t_clist *neighbors)
{
register unsigned int monovalent, _i,_j,_k;
unsigned int temp;
_i=ids->size;
while (_i--)
       if (!neighbors->list[ids->list[_i]].size) ids->list[_i]=ids->list[--ids->size];
  else if ( neighbors->list[ids->list[_i]].size==1)
         {
         monovalent=ids->list[_i];
         ids->list[_i]=ids->list[--ids->size];
         do{//remove monovalents chain
           _k=ids->size;
           while (--_k)
             if (*neighbors->list[monovalent].list==ids->list[_k]) break;
           if ((temp=find_in_list(monovalent,&neighbors->list[ids->list[_k]]))==(unsigned int)-1) { ylib_errno=YERROR_INTERNAL_CODE; return (unsigned int)-1; }
           _j=neighbors->list[ids->list[_k]].list[temp];
           neighbors->list[ids->list[_k]].list[temp]=neighbors->list[ids->list[_k]].list[--neighbors->list[ids->list[_k]].size];
           neighbors->list[ids->list[_k]].list[neighbors->list[ids->list[_k]].size]=_j;
           neighbors->list[monovalent].size=0;
           if (!neighbors->list[ids->list[_k]].size) ids->list[_k]=ids->list[--ids->size];
           else if (neighbors->list[ids->list[_k]].size==1) { monovalent=ids->list[_k]; ids->list[_k]=ids->list[--ids->size]; }
                else break;
           }while((ids->size));
         _i=(_i>ids->size) ? ids->size : _i;
         }
return ids->size;
}

//This function find the closest path in graph between point A and point B setted as 0 and 1 element of path list
unsigned int find_closest_path_with_neighbors(unsigned int start,unsigned int end,t_clist *neighbors,unsigned int *route)
{
register unsigned int _i,_j,_k,_l;

if ((_l=enumerate_vertices_with_neighbors(start,end,neighbors,&route[neighbors->size],route))==(unsigned int)-1) //show_list(vertices)
  return FALSE;
//go backward
route[0x0]=start;
_j=route[_i=_l]=end;
while (--_i)
  {
  _k=neighbors->list[_j].size;
  while (--_k)
    if (route[neighbors->size+neighbors->list[_j].list[_k]]<route[neighbors->size+_j]) break;
  _j=route[_i]=neighbors->list[_j].list[_k];
  }
return _l+1;
}

//This function found all cycles inside the graph with aid of neighbors massive
t_clist *define_cycles_with_neighbors(t_clist *neighbors)
{
register unsigned int _i,_k;
unsigned int *sizes=0x0,order,end;
t_clist *cycles=0x0;
t_list *cycle=0x0, *ids=0x0;

//Prepare memory firstly
_i=neighbors->size;
if (!(cycles=(t_clist*)alloc_clist(FALSE,FALSE)))             return FALSE;
if (!(ids=(t_list*)alloc_list(_i)))                           goto LABEL_ERROR_0;
if (!(cycle=(t_list*)alloc_list(_i<<1)))                      goto LABEL_ERROR_1;
if (!(sizes=(unsigned int*)malloc(sizeof(unsigned int)*_i)))  { ylib_errno=YERROR_MEMORY; goto LABEL_ERROR_2; }
//Save sizes to rescue neighbours than
ids->size=_i;
while (_i--)
  {
  sizes[_i]=neighbors->list[_i].size;
  ids->list[_i]=_i;
  }
//Search for cycles
while ((cut_onevalent_with_neighbors(ids,neighbors)))
  {
  order=0x1; //init lowest vertices order
  //Stage B. Find cycles
  while (++order<ids->size)
    {//try i-th atom
    _i=ids->size;
    while (_i--)
      if (neighbors->list[ids->list[_i]].size==order)
        {//try j-th growing direction
        //Cut trial edge and init cycle growing
        end=ids->list[neighbors->list[ids->list[_i]].list[--neighbors->list[ids->list[_i]].size]];
        //hide edge
        _k=--neighbors->list[end].size;
        while (_k--)
          if (neighbors->list[end].list[_k]==ids->list[_i])
            {
            neighbors->list[end].list[_k]=neighbors->list[end].list[neighbors->list[end].size];
            neighbors->list[end].list[neighbors->list[end].size]=ids->list[_i];
            break;
            }
        //Growing cycle
        if ( (cycle->size=find_closest_path_with_neighbors(ids->list[_i],end,neighbors,cycle->list)))
          if (!(clist_add(cycle,&cycles))) goto LABEL_ERROR_3; //Define cycle
        goto NEXT_SEARCH;
        }
    }
  NEXT_SEARCH: ;
  }

//Restore neighbours
_i=neighbors->size;
while(_i--)
  neighbors->list[_i].size=sizes[_i];

free(cycle);
free(sizes);
free(ids);
return cycles;
LABEL_ERROR_3:
_i=neighbors->size;
while(_i--)
  neighbors->list[_i].size=sizes[_i];
                      free(sizes);
LABEL_ERROR_2: free(cycle);
LABEL_ERROR_1: free(ids);
LABEL_ERROR_0: free(cycles);
return FALSE;
}

//This function found all cycles inside the graph
t_clist *define_cycles(unsigned int nvertices,unsigned int nedges,t_edge *edges,t_clist **neighbors)
{
register t_clist *cycles;
if (!((*neighbors)=define_neighbors(FALSE,nvertices,nedges,edges))) return FALSE;
else cycles=define_cycles_with_neighbors(*neighbors);
return cycles;
}

//This function builds list of vertices that lays on the vertice side of edge
t_list *get_fragment_with_neighbors(t_clist *neighbors,unsigned int root_vertice,t_edge *edge)
{
t_list *fragment;
char *t;
register unsigned int _i,_j;
register t_list *vp;

if (!(fragment=(t_list*)alloc_list(neighbors->size))) return FALSE;
if (!(t=(char*)calloc(neighbors->size,sizeof(char)))) { LABEL_ERROR_MEMORY: free(fragment); ylib_errno=YERROR_MEMORY; return FALSE; }
//Create init list excluding opposite edge's vertice
fragment->size=1;
  *fragment->list=edge->vertice[root_vertice];
t[*fragment->list]=TRUE;
_i=neighbors->list[*fragment->list].size;
while (_i--)
  if (neighbors->list[*fragment->list].list[_i]!=edge->vertice[(unsigned int)(!root_vertice)])
    {
    fragment->list[fragment->size]=neighbors->list[*fragment->list].list[_i]; 
    t[fragment->list[fragment->size++]]=TRUE;
    }
//Collect neighbors
_i=0;
while(++_i!=fragment->size)
  {
  _j=neighbors->list[fragment->list[_i]].size;
  while (_j--)
    if (!t[neighbors->list[fragment->list[_i]].list[_j]])
      {
      fragment->list[fragment->size++]=neighbors->list[fragment->list[_i]].list[_j];
      t[neighbors->list[fragment->list[_i]].list[_j]]=TRUE;
      }
  }
//Free some memory
free(t);
if (!(vp=(t_list*)realloc(fragment,sizeof(t_list)+sizeof(unsigned int)*fragment->size))) goto LABEL_ERROR_MEMORY;
else { fragment=(t_list*)vp; fragment->list=(void*)fragment+sizeof(t_list); }
return fragment;
}


//This function colours all subgraphs into a concequtive colours set (starting #1). The amount of used colours is returned.
//NB! It stores amount of painted by each colour vertices in the ${count} first cells of the buff massive
inline unsigned int paint_subgraph_with_neighbors(t_clist *neighbors,unsigned int *buff,unsigned int *color)
{
register unsigned int _i, _j, _k, _n, size;
unsigned int count=0;
memset(color,0x0,sizeof(unsigned int)*neighbors->size);
_n=0;
while (_n!=neighbors->size)
  {
  ++count, _i=0; while ( (color[_i])) _i++; _n++, buff[count-1]=_i, color[buff[count-1]]=count, _i=count-1, size=count; 
  do { _j=neighbors->list[buff[_i]].size; while (_j--) if (!(color[_k=neighbors->list[buff[_i]].list[_j]])) { _n++, buff[size++]=_k, color[_k]=count; } } while (++_i!=size) ;
  buff[count-1]=size; //Save the amount of painted vertices
  }
return count;
}

//This function cuts all but the biggest connected vertices fragment
//It returns FALSE on error, NTNF if not editions to the input graph were made and TRUE if there were successful corrections
char cut_the_biggest_subgraph(unsigned int *nvertices,int *vertices,unsigned int *nedges,t_edge *edges) 
{
register unsigned int _i, _j; 
unsigned int *color, *id;
t_clist *neighbors;
if (!(neighbors=define_neighbors(NULL,*nvertices,*nedges,edges))) return FALSE;
if (!(id=(unsigned int*)malloc(2*sizeof(unsigned int)*(*nvertices)))) { ylib_errno=YERROR_MEMORY; free(neighbors); return FALSE; } else color=&id[*nvertices];
if ((_i=paint_subgraph_with_neighbors(neighbors,id,color))!=1)
  {//Many fragments
  _j=0; while ( (--_i)) if (id[_j]<id[_i]) _j=_i; _j++; //Find the 'largest' color in subgraph
  _i=0, *nvertices=0; do { --_i; if (color[_i]==_j) color[_i]=(*nvertices)++; else color[_i]=(unsigned int)-1;  } while (*nvertices!=id[_j]) ; //Generate new ids for vertices
  _i=_j=0; do { if (color[_i]!=(unsigned int)-1) vertices[_j++]=vertices[_i]; _i++; } while ( (_i!=*nvertices)) ;  //Shift vertices
  _i=_j=0; while (_i!=*nedges)  //Shift edges
             {
             if ( (edges[_i].vertice[0]<(*nvertices))&&(color[edges[_i].vertice[0]]!=(unsigned int)-1)&&
                  (edges[_i].vertice[1]<(*nvertices))&&(color[edges[_i].vertice[1]]!=(unsigned int)-1) )  
               { edges[_j].vertice[0]=color[edges[_i].vertice[0]], edges[_j].vertice[1]=color[edges[_i].vertice[1]], edges[_j++].type=edges[_i].type; }
             _i++;
             }
  *nedges=_j; free(id), free(neighbors); return TRUE;
  }
else { free(id), free(neighbors); return NTNF; }
}





/*********************    ( S U B - ) S T R U C T U R E S      S E A R C H I N G      P A R T    ********************************************/




//Heuristics #1. Convert the graph into a sorted bi-int string of virtual colors and compare two virtual colors lines.

//This function returns negative, 0 or positive integer if vcolor0 is less, equal or bigger than vcolor1 correspondingly
inline int _compare_vcolor_descriptors(const void *vcolor0,const void *vcolor1)
{
if (((int*)vcolor0)[0]==((int*)vcolor1)[0]) return (((int*)vcolor0)[1]-((int*)vcolor1)[1]);
else                                        return (((int*)vcolor0)[0]-((int*)vcolor1)[0]);
}
//This function converts graph into a sorted array of "virtually collored" vertices
//The synthetic 'color' scheme is | (vertice_type, n_neighbors) | (neighbor_type_#1, edge_type_#1) | (neighbor_type_#2, edge_type_#2 ) | ... 
//The array is sorted by desceding inner sorts
inline unsigned int **paint_virtual_colors_with_neighbors(register t_clist *neighbors,register int *vertices,register unsigned int nedges,register t_edge *edges)
{
register unsigned int _i;
register unsigned int **vcolors;
//Stage I. Allocate memory
if (!(vcolors=(unsigned int**)malloc(sizeof(int*)*neighbors->size+0x2*sizeof(int)*(neighbors->size+nedges*2)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else for (vcolors[0]=(unsigned int*)((void*)vcolors+sizeof(int*)*neighbors->size), _i=1; _i<neighbors->size; _i++) vcolors[_i]=vcolors[_i-1]+2*(1+neighbors->list[_i-1].size);
//Stage II. Gather colors
_i=neighbors->size; while (_i--) { vcolors[_i][0]=vertices[_i], vcolors[_i][1]=0; }
_i=nedges;
while (_i--)
  {
  vcolors[edges[_i].vertice[0]][2*vcolors[edges[_i].vertice[0]][1]+0]=edges[_i].vertice[1], vcolors[edges[_i].vertice[0]][2*vcolors[edges[_i].vertice[0]][1]+1]=edges[_i].type, vcolors[edges[_i].vertice[0]][1]++;
  vcolors[edges[_i].vertice[1]][2*vcolors[edges[_i].vertice[1]][1]+0]=edges[_i].vertice[0], vcolors[edges[_i].vertice[1]][2*vcolors[edges[_i].vertice[1]][1]+1]=edges[_i].type, vcolors[edges[_i].vertice[1]][1]++;
  }
//Stage III. Sort in-line color descriptors
_i=neighbors->size; while (_i--) if (!(v_qsort(vcolors[_i][1],(void*)&vcolors[_i][2],_compare_vcolor_descriptors))) { free(vcolors), vcolors=0x0; return FALSE; }
return vcolors;  
}

//This function returns negative, 0 or positive integer if vcolor0 is less, equal or bigger than vcolor1 correspondingly
inline int _compare_vcolor_vectors(const void *vcolor0,const void *vcolor1)
{
register unsigned int _i; 
if (*((int*)vcolor0)==*((int*)vcolor1))
  {
  vcolor0+=sizeof(int), vcolor1+=sizeof(int);
  if (*((int*)vcolor0)==*((int*)vcolor1))
    {
    _i=*((int*)vcolor0)*2; while (_i--) { vcolor0+=sizeof(int), vcolor1+=sizeof(int); if (*((int*)vcolor0)!=*((int*)vcolor1)) return (*((int*)vcolor0)-*((int*)vcolor1)); }
    return FALSE;
    }
  } 
return (*((int*)vcolor0)-*((int*)vcolor1));
}
//This function creates arranged vector (useful for in-graph comparing/search)
//inline unsigned int *arrange_vcolors(register unsigned int size,register unsigned int **vcolors)
//{
//register unsigned int _i, *_id;
//if (!(_id=(unsigned int*)malloc(sizeof(unsigned int)*size))) { ylib_errno=YERROR_MEMORY; return FALSE; } else { _i=size, _id+=_i; while (_i--) *(--_id)=_i; }
//if (!(object_id_qsort(size,_id,(void**)vcolors,_compare_vcolor_vectors))) { free(_id), _id=0x0; return FALSE; }
//return _id;
//}

//This function compares two lines of vcolors
inline int compare_vcolors(unsigned int size,unsigned int *id0,unsigned int **vcolors0,unsigned int *id1,unsigned int **vcolors1)
{
register unsigned int _i, _j, _k, _l;
_k=size;
while (_k--) 
  {
  _i=id0[_k], _j=id1[_k];
  if (vcolors0[_i][0]==vcolors1[_j][0])
     {
     if (vcolors0[_i][1]==vcolors1[_j][1])
       { _l=vcolors0[_i][0]; while (_l--) if (vcolors0[_i][_l+2]!=vcolors1[_j][_l+2]) return (vcolors0[_i][_l+2]-vcolors1[_j][_l+2]); }
     else return (vcolors0[_i][1]-vcolors1[_j][1]);
     }
  else return (vcolors0[_i][0]-vcolors1[_j][0]);
  }
return 0; 
}


//This function efficiently compares two strutures. It returns -1 on failure, TRUE if the structures match and FALSE otherwise.
//It eventually comes to configurational scanning but there are heuristics to facilitate the search in practice.
//Heuristics #1. Define multiplicity of identical elements in searched object. If there is at least an element with insufficient matching then the blocks are not equivalent.
//Heuristics #2. Figure out the optimal search tree is before actual searcheding by trying all vertices (starting from rarer one)
//Heuristics #3. Guild search with tree-like structure to minimize disturbance of already perfectly fitted brnches.
//Eventually do optimal tree-guided deep first search 
//Note. b?[_i][_j] is a massive of bonds that join given atom and it's corresponding neighbors from n?->list[_i]->list[_j] 
char compare_two_chemical_structures(unsigned int *accordance,unsigned int nv,int *vi,t_clist *ni,unsigned char **bi,int *vj,t_clist *nj,unsigned char **bj)
{
register unsigned int _i, _j, _k, _l, _t;
register float _f, max_f;
t_clist *correspondence;
unsigned int *step, *stack, *size, *mantissa, hroot, stack_size;
float *score;
char *flag;
void *vp;

//Heuristic #1. Find amount of suitable correspondences
if (!(correspondence=(t_clist*)malloc(sizeof(t_clist)+sizeof(t_list)*nv+sizeof(unsigned int)*nv*nv))) 
  { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return (char)-1; }
else 
  {
  correspondence->size=nv, correspondence->_size=0;
  correspondence->list=(void*)correspondence+sizeof(t_clist);
  correspondence->_items=(void*)correspondence+sizeof(t_clist)+sizeof(t_list)*correspondence->size;
  }

for (_i=0; _i<nv; correspondence->_size+=correspondence->list[_i++].size) 
  {
  correspondence->list[_i].size=0, correspondence->list[_i].list=&correspondence->_items[correspondence->_size];
   _j=nv;
  while (_j--)
    if (ni->list[_i].size==nj->list[_j].size)
      {
      _k=ni->list[_i].size;
      while (_k--)
        {
        _l=nj->list[_j].size;
        while (_l--) 
          if ( (bj[_j][nj->list[_j].list[_l]]>0)&&(bi[_i][ni->list[_i].list[_k]]==bj[_j][nj->list[_j].list[_l]])&&(vi[ni->list[_i].list[_k]]==vj[nj->list[_j].list[_l]]) )
            { bj[_j][nj->list[_j].list[_l]]=-bj[_j][nj->list[_j].list[_l]]; goto NEXT_K; }
        _l=nj->list[_j].size;
        while (_l--)
          if ((char)bj[_j][nj->list[_j].list[_l]]<0)
            { bj[_j][nj->list[_j].list[_l]]=-bj[_j][nj->list[_j].list[_l]]; }
        goto NEXT_I;
        NEXT_K: ;
        }
      _l=nj->list[_j].size;
      while (_l--)
        if ((char)bj[_j][nj->list[_j].list[_l]]<0)
          { bj[_j][nj->list[_j].list[_l]]=-bj[_j][nj->list[_j].list[_l]]; }
      correspondence->list[_i].list[correspondence->list[_i].size++]=_j;
      }
  if (!(correspondence->list[_i].size)) { free(correspondence); correspondence=0x0; return FALSE; } //Graphs are not compatible
  NEXT_I: ;
  }
if (!(vp=(t_clist*)realloc(correspondence,sizeof(t_clist)+sizeof(t_list)*correspondence->size+sizeof(unsigned int)*correspondence->_size))) 
  { LABEL_MEMORY_ERROR_1: free(correspondence); correspondence=0x0; goto LABEL_MEMORY_ERROR_0; }
else if ((t_clist*)vp!=correspondence)
       {
       correspondence=(t_clist*)vp;
       correspondence->list=(void*)correspondence+sizeof(t_clist);
       correspondence->_items=(void*)correspondence+sizeof(t_clist)+sizeof(t_list)*correspondence->size;
       for (correspondence->list[0].list=correspondence->_items, _i=1; _i<correspondence->size; _i++)
         correspondence->list[_i].list=correspondence->list[_i-1].list+correspondence->list[_i-1].size;
       }

//Heuristic #2. Find the most efficient hierarchical tree score=П(i=1...N)[exp(step[i]*num[i])] but summ with implicit exponent order
if (sizeof(float)>sizeof(unsigned int)) { if (!(stack=(unsigned int*)malloc(sizeof(float)*nv)))        goto LABEL_MEMORY_ERROR_1; }
else                                    { if (!(stack=(unsigned int*)malloc(sizeof(unsigned int)*nv))) goto LABEL_MEMORY_ERROR_1; }
if (!(step=(unsigned int*)malloc(sizeof(unsigned int)*nv))) { LABEL_MEMORY_ERROR_2: free(stack); stack=0x0; goto LABEL_MEMORY_ERROR_1; }

if (enumerate_vertices_with_neighbors(0,(unsigned int)-1,ni,step,stack)==(unsigned int)-1) 
  { LABEL_MEMORY_ERROR_3: free(step); step=0x0; goto LABEL_MEMORY_ERROR_2; }
else { max_f=0., _j=nv; while (_j--) max_f+=(float)step[_j]*(float)correspondence->list[_j].size; }
_i=nv, hroot=0;
while (--_i)
  if (enumerate_vertices_with_neighbors(_i,(unsigned int)-1,ni,step,stack)==(unsigned int)-1) goto LABEL_MEMORY_ERROR_3;
  else { _f=0., _j=nv; while (_j--) _f+=(float)step[_j]*(float)correspondence->list[_j].size; if (_f>max_f) { max_f=_f, hroot=_i; } }
if (enumerate_vertices_with_neighbors(hroot,(unsigned int)-1,ni,step,stack)==(unsigned int)-1) goto LABEL_MEMORY_ERROR_3;

//Heuristic #3. Constuct the best hierarchical tree
if (!(score=(float*)calloc(nv,sizeof(float)))) goto LABEL_MEMORY_ERROR_3;
if (!(size=(unsigned int*)malloc(sizeof(unsigned int)*nv))) 
  { LABEL_MEMORY_ERROR_4: free(score); score=0x0; goto LABEL_MEMORY_ERROR_3; }
else { _i=nv; while (_i--) size[_i]=ni->list[_i].size; }
  
//Stage H.3.1. Construct default rtree
_i=nv; 
while (_i--) 
  if (_i!=hroot)
    {
    _j=size[_i];
    while (--_j) 
      {
           if (step[_i]<step[ni->list[_i].list[_j]])
             {
             if (step[*ni->list[_i].list]<step[_j]) goto RM_CYCLE_IN_HIERARCHY_TREE;
             else { _t=ni->list[_i].list[_j], ni->list[_i].list[_j]=ni->list[_i].list[size[_i]], ni->list[_i].list[size[_i]]=_t; }
             }
      else if (step[_i]==step[ni->list[_i].list[_j]])
             { RM_CYCLE_IN_HIERARCHY_TREE: _t=ni->list[_i].list[_j], ni->list[_i].list[_j]=ni->list[_i].list[size[_i]], ni->list[_i].list[size[_i]]=_t; size[_i]--; }
      }
    }
//Stage H.3.2. Score h-tree
for (_k=0; _k<size[hroot]; _k++)
  {
  _i=ni->list[hroot].list[_k];
  step[_i]=size[_i];
  SCANN_HTREE_H32:
  while (--stack[_i])
    {//forward
    _i=ni->list[_i].list[stack[_i]];
    step[_i]=size[_i];
    goto SCANN_HTREE_H32;
    }
  //Reverse
  _j=_i, _i=*ni->list[_j].list;
  score[_j]+=logf((float)correspondence->list[_j].size), score[_i]+=score[_j];
  if ((_i)!=hroot) goto SCANN_HTREE_H32;
  }
//Stage H.3.3. Reconstruct optimal h-tree using stakc and step as temporary storages
_j=size[hroot]; while (_j--) { step[_j]=ni->list[hroot].list[_j], ((float*)stack)[_j]=score[ni->list[hroot].list[_j]]; }
if (!(fi_qsort(size[hroot],((float*)stack),(int*)step))) goto LABEL_MEMORY_ERROR_4;
_j=size[hroot]; while (_j--) stack[_j]=(unsigned int)bi[hroot][step[_j]];
_j=size[hroot]; while (_j--) { bi[hroot][_j]=(unsigned char)stack[_j], ni->list[hroot].list[_j]=step[_j]; }
_i=nv;
while (_i--)
  if (_i!=hroot)
    {
    _j=size[_i]; while (--_j) { step[_j]=ni->list[_i].list[_j], ((float*)stack)[_j]=score[ni->list[_i].list[_j]]; }
    if (!(fi_qsort(size[_i]-1,&((float*)stack)[1],(int*)&step[1]))) goto LABEL_MEMORY_ERROR_4; 
    _j=size[_i]; while (--_j) stack[_j]=(unsigned int)bi[_i][step[_j]];
    _j=size[_i]; while (--_j) { bi[_i][_j]=(unsigned char)stack[_j]; ni->list[_i].list[_j]=step[_j]; }
    }
//Free some memory
free(score); score=0x0;

//Do deep first search 
if (!(mantissa=(unsigned int*)malloc(sizeof(unsigned int)*nv))) 
  { LABEL_MEMORY_ERROR_5: free(size);       size=0x0;       goto LABEL_MEMORY_ERROR_3; }
if (!(flag=(char*)calloc(nv,sizeof(char)))) 
  {                       free(mantissa);   mantissa=0x0;   goto LABEL_MEMORY_ERROR_5; }

mantissa[hroot]=correspondence->list[hroot].size;
NEXT_MANTISSA_HROOT: stack[0]=hroot;
while (mantissa[hroot]--)
  {
  for (flag[accordance[hroot]=correspondence->list[hroot].list[mantissa[hroot]]]=TRUE, stack_size=_k=0; _k<size[hroot]; _k++)
    {
    _i=ni->list[hroot].list[_k], stack[++stack_size]=_i, mantissa[_i]=correspondence->list[_i].size;
    while (mantissa[_i]--)
      if ( (!(flag[_j=correspondence->list[_i].list[mantissa[_i]]]))&&(find_in_row(accordance[hroot],nj->list[_j].size,nj->list[_j].list)!=(unsigned int)-1) )
        { //Success
        flag[accordance[_i]=correspondence->list[_i].list[mantissa[_i]]]=TRUE, step[_i]=size[_i];
        goto SCANN_HTREE_DFS;
        }
    goto SCANN_FAILURE_DFS; //Failure   
    SCANN_HTREE_DFS: ;
    while (--step[_i])
      {//forward
      _i=ni->list[_i].list[step[_i]], stack[++stack_size]=_i, mantissa[_i]=correspondence->list[_i].size;
      do{
        while (mantissa[_i]--)
          if ( (!(flag[_j=correspondence->list[_i].list[mantissa[_i]]]))&&(find_in_row(accordance[*ni->list[_i].list],nj->list[_j].size,nj->list[_j].list)!=(unsigned int)-1) )
            { flag[accordance[_i]=correspondence->list[_i].list[mantissa[_i]]]=TRUE, step[_i]=size[_i]; goto SCANN_HTREE_DFS; } //Success
        SCANN_FAILURE_DFS: ; //Failure  
        _i=stack[--stack_size];
        flag[accordance[_i]]=FALSE;
        } while ( (stack_size)) ;
      goto NEXT_MANTISSA_HROOT;
      }
    //Reverse
    _j=_i, _i=*ni->list[_j].list;
    if ((_i)!=hroot) goto SCANN_HTREE_DFS;
    }
  //The exact matching of the to given graphs is found!
  free(correspondence); correspondence=0x0;
  free(step); step=0x0;
  free(stack); stack=0x0;
  free(size); size=0x0;
  free(mantissa); mantissa=0x0;
  free(flag); flag=0x0;
  return TRUE;
  }  
//There is no exact matching between two given graphs :(
free(correspondence); correspondence=0x0;
free(step); step=0x0;
free(stack); stack=0x0;
free(size); size=0x0;
free(mantissa); mantissa=0x0;
free(flag); flag=0x0;
return FALSE;
}





//*****************************************               A D J A C E N C Y      P A R T              *****************************************/





//This function construct an adjacency structure for the clustered graph
// | adjacency_itself | &nv | &ne | size_nv | &&nv | &&ne |   
//Note. ajacency is solid massive which can be deleted with simple free(adjacency);
t_adjacency *construct_adjacency(unsigned int *_cvertices,t_clist *clusters,unsigned int nvertices,unsigned int nedges,t_edge *edges)
{
unsigned int _i, _j, _ne;
unsigned int *cvertices, *csize;
t_adjacency *adjacency=0x0;

//Stage I.1. Prepare memory.
if (!(csize=(unsigned int*)calloc(clusters->size,sizeof(unsigned int)))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY; return FALSE; }
if (!(_cvertices))
  {
  if (!(cvertices=(unsigned int*)malloc(sizeof(unsigned int)*nvertices))) { LABEL_MEMORY_ERROR_1: free(csize); csize=0x0; goto LABEL_MEMORY_ERROR_0; }
  //Stage I.2. Map clusters on vertices.
  _i=clusters->size; while (_i--) { _j=clusters->list[_i].size; while (_j--) cvertices[clusters->list[_i].list[_j]]=_i; }
  }
else cvertices=_cvertices;
//Stage I.3. Gather statistics.
_ne=0, _i=nedges;
while (_i--)
  if (cvertices[edges[_i].vertice[0]]!=cvertices[edges[_i].vertice[1]])
    { csize[cvertices[edges[_i].vertice[0]]]++, csize[cvertices[edges[_i].vertice[1]]]++, _ne++; }
//Stage I.4. Alloc memory
if (!(adjacency=(t_adjacency*)malloc(sizeof(t_adjacency)+clusters->size*(sizeof(unsigned int)+sizeof(unsigned int*)+sizeof(t_edge*))+0x2*_ne*(sizeof(unsigned int)+sizeof(t_edge)))))
  { if (!(_cvertices)) { free(cvertices); cvertices=0x0; } goto LABEL_MEMORY_ERROR_1; }
else
  {
  adjacency->size_v=clusters->size, adjacency->size_e=_ne;
  adjacency->nn=(void*)adjacency+sizeof(t_adjacency);
  adjacency->nv=(void*)adjacency+sizeof(t_adjacency)+adjacency->size_v*(sizeof(unsigned int));
  adjacency->ne=(void*)adjacency+sizeof(t_adjacency)+adjacency->size_v*(sizeof(unsigned int)+sizeof(unsigned int*));
  adjacency->nv[0]=(void*)adjacency+sizeof(t_adjacency)+adjacency->size_v*(sizeof(unsigned int)+sizeof(unsigned int*)+sizeof(t_edge*));
  adjacency->ne[0]=(void*)adjacency+sizeof(t_adjacency)+adjacency->size_v*(sizeof(unsigned int)+sizeof(unsigned int*)+sizeof(t_edge*))+0x2*adjacency->size_e*sizeof(unsigned int);
  for (adjacency->nn[0]=0, _i=1; _i<clusters->size; _i++)
    { adjacency->nn[_i]=0, adjacency->nv[_i]=adjacency->nv[_i-1]+csize[_i-1], adjacency->ne[_i]=adjacency->ne[_i-1]+csize[_i-1]; }
  }
//Stage II. Fill adjacency
_i=nedges;
while (_i--)
  if (cvertices[edges[_i].vertice[0]]!=cvertices[edges[_i].vertice[1]])
    {
    adjacency->nv[cvertices[edges[_i].vertice[0]]][adjacency->nn[cvertices[edges[_i].vertice[0]]]]=cvertices[edges[_i].vertice[1]];
    adjacency->ne[cvertices[edges[_i].vertice[0]]][adjacency->nn[cvertices[edges[_i].vertice[0]]]].vertice[0]=edges[_i].vertice[0];
    adjacency->ne[cvertices[edges[_i].vertice[0]]][adjacency->nn[cvertices[edges[_i].vertice[0]]]].vertice[1]=edges[_i].vertice[1];
    adjacency->ne[cvertices[edges[_i].vertice[0]]][adjacency->nn[cvertices[edges[_i].vertice[0]]]].type=_i;
    adjacency->nn[cvertices[edges[_i].vertice[0]]]++;
    adjacency->nv[cvertices[edges[_i].vertice[1]]][adjacency->nn[cvertices[edges[_i].vertice[1]]]]=cvertices[edges[_i].vertice[0]];
    adjacency->ne[cvertices[edges[_i].vertice[1]]][adjacency->nn[cvertices[edges[_i].vertice[1]]]].vertice[0]=edges[_i].vertice[1];
    adjacency->ne[cvertices[edges[_i].vertice[1]]][adjacency->nn[cvertices[edges[_i].vertice[1]]]].vertice[1]=edges[_i].vertice[0];
    adjacency->ne[cvertices[edges[_i].vertice[1]]][adjacency->nn[cvertices[edges[_i].vertice[1]]]].type=_i;
    adjacency->nn[cvertices[edges[_i].vertice[1]]]++;
    }
//Free memory and exit.
free(csize); csize=0x0;
if (!(_cvertices)) { free(cvertices); cvertices=0x0; }
return adjacency;
}


/*
//This function performs optimized memory allocation neighbours search in two passes.
//NOTE. If you use fragments=list it will reorder edges to improve its performance
t_clist *define_neighbors(t_list *vertex,unsigned int nvertices,unsigned int nedges,t_edge *edges)
{
register unsigned int _i;
register t_clist *neighbors;
unsigned int j,k;

if (!nvertices) { LABEL_ERROR_DATA_CONSISTMENT: ylib_errno=YERROR_DATA_CONSISTMENT; return FALSE; }
if (!vertex)
  {
  if (!(neighbors=(t_clist*)alloc_clist(nvertices,nedges*2))) { LABEL_ERROR_MEMORY: ylib_errno=YERROR_MEMORY; return FALSE; }
  _i=nvertices;
  while(_i--) neighbors->list[_i].size=0x0;
  _i=nedges;
  while (_i--) { neighbors->list[edges[_i].vertice[0x0]].size++; neighbors->list[edges[_i].vertice[0x1]].size++; }
  neighbors->list[0].list=neighbors->_items;
  for (_i=1;_i<nvertices;_i++)
    {
    neighbors->list[_i].list=neighbors->list[_i-1].list+neighbors->list[_i-1].size;
    neighbors->list[_i-1].size=0;
    }
  neighbors->list[_i-1].size=0;
  _i=nedges;
  while (_i--)
    {
    j=edges[_i].vertice[0x0];
    k=edges[_i].vertice[0x1];
    neighbors->list[k].list[neighbors->list[k].size++]=j;
    neighbors->list[j].list[neighbors->list[j].size++]=k;
    }
  }
else
  {
  if (!(neighbors=(t_clist*)alloc_clist(vertex->size,nedges*2))) goto LABEL_ERROR_MEMORY;
  _i=vertex->size;
  while (_i--) neighbors->list[_i].size=0x0;
  for (_i=0;_i<nedges;)
    if ( ( (find_in_list(&j,edges[_i].vertice[0x0],vertex)))&&( (find_in_list(&k,edges[_i].vertice[0x1],vertex))) )
      {
      neighbors->list[j].size++;
      neighbors->list[k].size++;
      _i++;
      }
    else
      {//update ne to minimize searches during second pass
      j=edges[_i].vertice[0x0], edges[_i].vertice[0x0]=edges[nedges].vertice[0x0], edges[nedges].vertice[0x0]=j;
      j=edges[_i].vertice[0x1], edges[_i].vertice[0x1]=edges[nedges].vertice[0x1], edges[nedges].vertice[0x1]=j;
      j=edges[_i].type,         edges[_i].type=edges[nedges].type,                 edges[nedges].type=j;
      nedges--;
      }
  if (!nedges) { free(neighbors); goto LABEL_ERROR_DATA_CONSISTMENT; }
  for (neighbors->list[_i].list=neighbors->_items, _i=1;_i<nvertices;_i++)
    {
    neighbors->list[_i].list=neighbors->list[_i-1].list+neighbors->list[_i-1].size;
    neighbors->list[_i-1].size=0;
    }
  neighbors->list[_i-1].size=0;
  _i=nedges;
  while (_i--)
    {
    if ( (!(find_in_list(&j,edges[_i].vertice[0x0],vertex)))||(!(find_in_list(&k,edges[_i].vertice[0x1],vertex))) ) { free(neighbors); goto LABEL_ERROR_DATA_CONSISTMENT; }
    neighbors->list[k].list[neighbors->list[k].size++]=j;
    neighbors->list[j].list[neighbors->list[j].size++]=k;
    }
  }
return neighbors;
}
*/


