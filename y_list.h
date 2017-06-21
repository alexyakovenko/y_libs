//This file contain various t_list operation routines

#define Y_LIST   0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_TXT
#include "y_txt.h"
#endif
#ifndef Y_FILE
#include "y_file.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif

typedef struct{
              unsigned int size;           //Number of elements in list
              unsigned int *list;          //elements of list
              } t_list;                    //odinary list structure

typedef struct{
              unsigned int size;           //number of lists in complex list
              unsigned int *_items;        //members of complex list
              unsigned int  _size;         //number of members of complex list
              t_list       *list;          //massive of lists in complex list
              }t_clist;                    //Complex list strcture

typedef struct{
              unsigned int it[2];          //elements of pair   
              }t_pair;                     //pair of int elements

typedef struct{
              unsigned int size;           //Size of objects
              void      *object;           //Object itself
              }t_object;

//This structure defines the branch structutre
typedef struct{
              void         *it;                    //The pointer to object that belongs to current node of tree
              unsigned int  root;                  //The root branch
              unsigned int  size;                  //The number of children branches
              unsigned int *children;              //The children branches id in container
              }t_branch;
//This structure defines the tree
typedef struct{
              unsigned int  size;                   //The number of branches on the tree
              t_branch     *branches;               //The branches array
              }t_tree;



//t_list  EMPTY_LIST={ 0, 0x0 };
//t_clist EMPTY_CLIST={ 0, 0x0, 0, 0x0 }; 
//t_clist EMPTY_UCLIST={ 1, 0x0, 0x0, &EMPTY_LIST };

#define ISORT_ALGORITHM_CUTOFF 0xF

//Empty list
#define NULL_LIST (t_list*)("\0\0\0\0\0\0\0\0") // {list->size=0; list->list=0;}

//Redefine operator macross
#define LIST_OPR(A,OPERATOR,B) (      IF (OPERATOR=='+')      or_list(A,B) \
                                 ELSE IF (OPERATOR=='-')     not_list(A,B) \
                                 ELSE IF (OPERATOR=='&')     and_list(A,B) \
                                 ELSE IF (OPERATOR=='|') overlap_list(A,B) )


/*
 * There are three types of lists:
 * 1. LIST - just a list of unsigned ints
 * 2. COMPLEX_LIST - list of pointers to objects (of the variable lenght)
 * 3. OBJECT_LIST - list of pointers into included in list objects (of the same lenght)
 */

//This is debuging function
void show_list(t_list *list);

//This is a lazy function for list importing
t_list *import_list(char *file_name);

//This function reads list from file
t_list *read_list(register FILE *in);

//This function reads clist from file
t_clist *read_clist(register FILE *in);

//This function write list to binary file
inline char write_list(register FILE *out,register t_list *list);

//This function reads clist from file
char write_clist(register FILE *out,t_clist *clist);

//This function exports list to txt
char export_list(FILE *out,t_list *list);

//This function alloc memory for t_list
inline t_list* alloc_list(register unsigned int size);

//This function alloc memory  for t_clist. Memory is not initialized!
inline t_clist *alloc_clist(register unsigned int size,register unsigned int _size);

//This function makes copy of complex list
t_clist *copy_clist(t_clist *src);

//THis function alloc objects
t_object *alloc_object(unsigned int size);

//This function realloc t_list.
inline t_list* realloc_list(register t_list* list,register unsigned int size);


//Search functions

//This function find object in row of values (any sizeof(int) datatypes)
inline unsigned int  find_in_row(register unsigned int item,register unsigned int size,register unsigned int *row);
//This function find object in row of values (any sizeof(int) datatypes) but eXcludes zero element from consideraton
inline unsigned int xfind_in_row(register unsigned int item,register unsigned int size,register unsigned int *row);

//Fast finders. They require lists to be pre-sorted. They returns TRUE and the position if match is found and FALSE and the position to paste a new value if there is no matching.

//This function finds the value in a row of _sorted_ unsigned ints. It expects the sort order is from low to high
inline char find_in_sorted_lth_urow(unsigned int *_id,register unsigned int item,register unsigned int size,register unsigned int *urow);
//This function finds the value in a row of _sorted_ unsigned ints. It expects the sort order is from low to high
inline char find_in_sorted_lth_irow(unsigned int *_id,register int item,register unsigned int size,register unsigned int *irow);

//This function finds an object in the stack of objects
inline char find_in_sorted_lth_objects(unsigned int *_id,register void *item,unsigned int size,register void *objects,register size_t object_size,int (*compare_objects)(const void* object0,const void* object1));

//The list wrappers

//This function find object in list massive
inline unsigned int  find_in_list(register unsigned int item,register t_list *list);
//This function find object in list massive but eXcldes zero element from consideration
inline unsigned int xfind_in_list(register unsigned int item,register t_list *list);

//This function finds object in sorted list
inline char find_in_sorted_lth_list(unsigned int *_id,register unsigned int item,register t_list *list);

//This function find object in the row of lists
inline char find_in_clist(register unsigned int *_id,register unsigned int item,register t_clist *clist);


//This funtion check is two rows has any overlap elements
//NOTE it return TRUE if there is an overlap and two (optional) unsigned ints with id of overlaping element a_row[i]==b_row[j]
inline char overlap_row(unsigned int *i, unsigned int a_size,unsigned int *a_row,unsigned int *j,unsigned int b_size,unsigned int *b_row);

//This functions finds min and max values
//This function find the lowest value in a row of doubles
unsigned int find_min_d(register unsigned int size,register double *_f);

//This function find the lowest value by module in a row of doubles
unsigned int find_min_absd(register unsigned int size,register double *_f);

//This function find next min in row of doubles
char find_next_min_d(register unsigned int *id,unsigned int size,register double *_f);

//This function find next min by module in row of doubles
char find_next_min_absd(register unsigned int *id,unsigned int size,register double *_f);

//This function find the lowest value in a row of ints
unsigned int find_min_i(register unsigned int size,register int *_i);

//This function find next min in row of ints
char find_next_min_i(register unsigned int *id,unsigned int size,register int *_i);

//This function find the largest value in a row of doubles
unsigned int find_max_d(register unsigned int size,register double *_f);

//This function find the largest by module value in a row of doubles
unsigned int find_max_absd(register unsigned int size,register double *_f);

//This function find next max in row of doubles
char find_next_max_d(register unsigned int *id,unsigned int size,register double *_f);

//This function find next max by module in row of doubles
char find_next_max_absd(register unsigned int *id,unsigned int size,register double *_f);

//This function find the largest value in a row of ints
unsigned int find_max_i(register unsigned int size,register int *_i);

//This function find next max in row of ints
char find_next_max_i(register unsigned int *id,unsigned int size,register int *_i);

//This function add an element to list massive
inline char list_add(unsigned int item,t_list **list);

//This function add an list to lists massive
char clist_add(t_list *list,t_clist **clist);

//This function delete a list from lists massive
void clist_del(unsigned int id,t_clist *clist);

//This function finds summ of double row
inline double summ_in_rowf(register unsigned int num,register double *_f);

//This function finds summ of int row
inline int summ_in_rowi(unsigned int num,int *_i);



// ------------------- L O G I C   P A R T ---------------------



//This fuinction performs logic operaion or  at two t_lists (Union) c=(a+b)
t_list *or_list(register t_list *a,register t_list *b);

//This function performs logic operation not at two t_list (Exclude) c=(a-b)
t_list *not_list(register t_list *a,register t_list *b);

//This function performs logic operation and at two t_list (Cross) c=(a-(a-b))
t_list *and_list(register t_list *a,register t_list *b);

//This funtion check is two lists has any overlap elements
//NOTE it return TRUE if there is an overlap and two (optional) unsigned ints with id of overlaping elements a->list[i]==b->list[j]
inline char overlap_list(register unsigned int *i,register t_list *a,register unsigned int *j,register t_list *b);



//--------------------- ( Q U I C K )   S O R T I N G   A L G O R I T H M S   ---------------------



//This function sorts 3 ints with efficient if algorithm
inline void sort_iif3(int *i0,int *i1,int *i2);
//This function sorts 3 doubles with efficient if algorithm
inline void sort_dif3(double *d0,double *d1,double *d2);
//This function sorts 3 doubles-int pairs with efficient if algorithm
inline void sort_diif3(double *d0,double *d1,double *d2,int *i0,int *i1, int *i2);
//This function sorts 3 doubles-double pairs with efficient if algorithm
inline void sort_ddif3(double *d0,double *d1,double *d2,double *dd0,double *dd1,double *dd2);



//This function performs sorting of int array. 
//Sorting is performed with qsort algorithm from http://algolist.manual.ru/sort/quick_sort.php

//Sorting int array
inline void i_isort(register unsigned int size,register int *i);
inline char i_qsort(register unsigned int size,register int *i);

//Sorting unsigned int array. 
inline void u_isort(register unsigned int size,register unsigned int *u);
inline char u_qsort(register unsigned int size,register unsigned int *u);

//This function sorts array of doubles and moves the associated int array accordingly
inline void di_isort(register unsigned int size,register double *d,register int *ii);
inline char di_qsort(register unsigned int size,register double *d,register int *ii);

//This function sorts array of floats and moves the associated int array accordingly
inline void fi_isort(register unsigned int size,register float *f,register int *ii);
inline char fi_qsort(register unsigned int size,register float *f,register int *ii);

//This function sorts massives of ints and moves the associated double massive accordingly
inline void id_isort(register unsigned int size,register int *i,register double *dd);
inline char id_qsort(register unsigned int size,register int *i,register double *dd);

//This function sorts massives of ints and moves the associated int massive accordingly
inline void ii_isort(register unsigned int size,register int *i,register int *ii);
inline char ii_qsort(register unsigned int size,register int *i,register int *ii);

//This function performs sorting of unsigned int array and moves associated double massive. 
inline void ud_isort(register unsigned int size,register unsigned int *u,register double *dd);
inline char ud_qsort(register unsigned int size,register unsigned int *u,register double *dd);

//This function sorts massives of unsigned ints and moves the associated (unsigned) int massive accordingly
inline void ui_isort(register unsigned int size,register unsigned int *u,register int *ii);
inline char ui_qsort(register unsigned int size,register unsigned int *u,register int *ii);

//This function sorts massives of objects and moves the associated (unsigned) int massive accordingly
inline void object_id_isort(unsigned int size,size_t obj_size,void *obj,unsigned int *id,int (*compare_obj)(const void *obj0,const void *obj1));
//Note for lbstack[size] and ubstack[size] size=2+(unsigned int)(ln((float)amount)/LN_2), amount >= 0.
inline char object_id_qsort(unsigned int size,size_t obj_size,void *obj,unsigned int *id,int (*compare_obj)(const void *obj0,const void *obj1));

//Sorting of objects array moved by pointers
inline void v_isort(register unsigned int size,void **objects,int (*compare_objects)(const void* object0,const void* object1));
inline char v_qsort(register unsigned int size,void **objects,int (*compare_objects)(const void* object0,const void* object1));

//Sorting of objects array and then move physical items
inline void w_isort(register unsigned int size,size_t obj_size,void *objects,int (*compare_objects)(const void* object0,const void* object1),void (*swap_objects)(const void* object0,const void* object1));
inline char w_qsort(register unsigned int size,size_t obj_size,void *objects,int (*compare_objects)(const void* object0,const void* object1),void (*swap_objects)(const void* object0,const void* object1));



