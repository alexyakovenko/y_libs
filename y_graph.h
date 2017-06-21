//This module contain various routines connected with graphs

#define Y_GRAPH 0x1

#ifndef Y_LIST
#include "y_list.h"
#endif


//One edge structure
typedef struct{
              unsigned int  vertice[2];      //two points that connected with edge
              int  type;                     //edge type or multiplexity
              }t_edge;                       //universal edge structure

//This function try to find an edge of two vertices suggesting it can be either A-B or B-A. It returns (unsigned int)-1 if there is no edge
inline unsigned int find_edge_unordered(register unsigned int _a,register unsigned int _b,register unsigned int _n,register t_edge *_e);

//This function performs optimized memory allocation neighbours search in two passes.
//NOTE. If you use fragments=list it will reorder edges to improve its performance
t_clist *define_neighbors(t_clist *neighbors,unsigned int nvertices,unsigned int nedges,t_edge *edges);

//This function enumerate vertices from start till end with neighbors. It returns (unsigned int)-1
//NOTE. If you wona total enumeration use end==(unsigned int)-1
unsigned int enumerate_vertices_with_neighbors(unsigned int start,unsigned int end,t_clist *neighbors,unsigned int *step,unsigned int *stack);

//This routine calculates intra-graph distances matrix for subset dist_v of vertices in graph size_v, neighbors
char calculate_intragraph_distance_matrix(unsigned int **d,unsigned int size_vd,unsigned int *vd,t_clist *neighbors);

//This function do enumeration of atoms with neighbors, builds quasy-rtree structure of neighbors and maps cycles edges temp[...]==size
unsigned int enumerate_cycles_with_neighbors (unsigned int root,unsigned int size,unsigned int *enumerate,t_clist *neighbors,unsigned int *cycles,unsigned int *stack);

//This function cuts onevalent vertices
unsigned int cut_onevalent_with_neighbors(t_list *ids,t_clist *neighbors);

//This function find the closest path in graph between point A and point B setted as 0 and 1 element of path list
unsigned int find_closest_path_with_neighbors(unsigned int start,unsigned int end,t_clist *neighbors,unsigned int *route);

//This function found all cycles inside the graph with aid of neighbors massive
t_clist *define_cycles_with_neighbors(t_clist *neighbors);

//This function found all cycles inside the graph
t_clist *define_cycles(unsigned int nvertices,unsigned int nedges,t_edge *edges,t_clist **neighbors);

//This function builds list of vertices that lays on the vertice side of edge
t_list *get_fragment_with_neighbors(t_clist *neighbors,unsigned int root_vertice,t_edge *edge);

//This function colours all subgraphs into a concequtive colours set (starting #1). The amount of used colours is returned.
//NB! It stores amount of painted by each colour vertices in the ${count} first cells of the buff massive
inline unsigned int paint_subgraph_with_neighbors(t_clist *neighbors,unsigned int *buff,unsigned int *colors);

//This function cuts all but the biggest connected vertices fragment
//It returns FALSE on error, NTNF if not editions to the input graph were made and TRUE if there were successful corrections
char cut_the_biggest_subgraph(unsigned int *nvertices,int *vertices,unsigned int *nedges,t_edge *edges);





/*********************    ( S U B - ) S T R U C T U R E S      S E A R C H I N G      P A R T    ********************************************/





//This function converts graph into a sorted array of "virtually collored" vertices
//The synthetic 'color' scheme is | (vertice_type, n_neighbors) | (neighbor_type_#1, edge_type_#1) | (neighbor_type_#2, edge_type_#2 ) | ... 
//The array is sorted by desceding inner sorts
inline unsigned int **paint_virtual_colors_with_neighbors(register t_clist *neighbors,register int *vertices,register unsigned int nedges,register t_edge *edges);
//This function creates arranged vector (useful for in-graph comparing/search)
inline unsigned int *arrange_vcolors(register unsigned int size,register unsigned int **vcolors);
//This function compares two lines of vcolors
inline int compare_vcolors(unsigned int size,unsigned int *id0,unsigned int **vcolors0,unsigned int *id1,unsigned int **vcolors1);

//This function efficiently compares two strutures. It returns -1 on failure, TRUE if the structures match and FALSE otherwise.
//It eventually comes to configurational scanning but there are heuristics to facilitate the search in practice.
//Heuristics #1. Define multiplicity of identical elements in searched object. If there is at least an element with insufficient matching then the blocks are not equivalent.
//Heuristics #2. Figure out the optimal search tree is before actual searcheding by trying all vertices (starting from rarer one)
//Heuristics #3. Guild search with tree-like structure to minimize disturbance of already perfectly fitted brnches.
//Eventually do optimal tree-guided deep first search 
//Note. b?[_i][_j] is a massive of bonds that join given atom and it's corresponding neighbors from n?->list[_i].list[_j] 
char compare_two_chemical_structures(unsigned int *accordance,unsigned int nv,int *vi,t_clist *ni,unsigned char **bi,int *vj,t_clist *nj,unsigned char **bj);





//*****************************************               A D J A C E N C Y      P A R T              *****************************************/





typedef struct{
              unsigned int size_v, size_e;     //amount of vertices and edges in the adjacency
              unsigned int  *nn;               //neighbors number of a 'cluster'
              unsigned int **nv;               //neighbors of a 'cluster'
              t_edge       **ne;               //edges between vertices of neighboring clusters: vertice[0] is own, vertixe[1] is neighbour 
              }t_adjacency;                    //adjacency of a clustered graph

//This function construct an adjacency structure for the clustered graph
// | adjacency_itself | &nv | &ne | size_nv | &&nv | &&ne |   
//Note. ajacency is solid massive which can be deleted with simple free(adjacency);
t_adjacency *construct_adjacency(unsigned int *_cvertices,t_clist *clusters,unsigned int nvertices,unsigned int nedges,t_edge *edges);


