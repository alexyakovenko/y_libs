//This file contain various t_list operation routines
#include <stdlib.h>
#include <string.h>
#include "y_list.h"


//This is debuging function
void show_list(t_list *list)
{
register unsigned int _i;
printf("\nList size = %d :\n",list->size);
for (_i=0;_i<list->size;_i++)
  printf("%d ",list->list[_i]);
printf("\n");
}

//This is a lazy function for list importing
t_list *import_list(char *file_name)
{
register t_list *list=0x0;
FILE *in=0x0;
register unsigned int _i;
char buffer[0xFF];
register char *lexem;

if ( ( (access(file_name,R_OK)))||(in=fopen(file_name,"r")) )
  {
  if ( ((fgets(buffer,0xFE,in)))&&((lexem=get_lex(0x1,buffer)))&&((check_int_type(lexem)))&&((list=(t_list*)alloc_list((unsigned int)atoi(lexem)))) )
    {
    for (_i=0;_i<list->size;_i++)
      if ( (!(fgets(buffer,0xFE,in)))||(!(lexem=get_lex(0x1,buffer)))||(!(check_int_type(lexem))) )
        {
        free(list);
        list=0x0;
        break;
        }
      else
        list->list[_i]=(unsigned int)atoi(lexem);
    }
  fclose(in);
  }
return list;
}

//This function reads list from file
t_list *read_list(register FILE *in)
{
unsigned int size;
t_list *list;
if ( (fread(&size,sizeof(unsigned int),0x1,in)==0x1)&&( (list=(t_list*)alloc_list(size))) )
  {
  if (fread(list->list,sizeof(unsigned int),list->size,in)==list->size)
    return list;
  else
    free(list);
  }
return FALSE;
}

//This function reads clist from file
t_clist *read_clist(register FILE *in)
{
unsigned int size,_size;
register t_clist *clist;
if ( (fread(&size,sizeof(unsigned int),0x1,in)!=0x1)||(fread(&_size,sizeof(unsigned int),0x1,in)!=0x1) )
  { LABEL_ERROR_IO: ylib_errno=YERROR_IO; return FALSE; }
if (!(clist=alloc_clist(size,_size))) return FALSE;
else clist->_size=_size;
if (fread(clist->_items,sizeof(unsigned int),clist->_size,in)!=clist->_size) { ERROR: free(clist); goto LABEL_ERROR_IO; }
if (clist->size)
  {
  if (fread(&clist->list[0x0].size,sizeof(unsigned int),0x1,in)!=0x1) goto ERROR;
  else clist->list[0x0].list=clist->_items;
  for (_size=1;_size<size;_size++)
    {
    if (fread(&clist->list[_size].size,sizeof(unsigned int),0x1,in)!=0x1) goto ERROR;
    else clist->list[_size].list=clist->list[_size-1].list+clist->list[_size-1].size;
    }
  }
return clist;
}

//This function write list to binary file
inline char write_list(register FILE *out,register t_list *list)
{
return ( (fwrite(&list->size,sizeof(unsigned int),0x1,out)==0x1)              &&
         (fwrite( list->list,sizeof(unsigned int),list->size,out)==list->size) );
}

//This function reads clist from file
char write_clist(register FILE *out,t_clist *clist)
{
register unsigned int _i;
if ( (fwrite(&clist->size,sizeof(unsigned int),0x1,out)!=0x1)||(fwrite(&clist->_size,sizeof(unsigned int),0x1,out)!=0x1) )
  { LABEL_ERROR_IO: ylib_errno=YERROR_IO; return FALSE; }
if (fwrite(clist->_items,sizeof(unsigned int),clist->_size,out)!=clist->_size) goto LABEL_ERROR_IO;
for (_i=0;_i<clist->size;_i++)
  if (fwrite(&clist->list[_i].size,sizeof(unsigned int),0x1,out)!=0x1) goto LABEL_ERROR_IO;
return TRUE;
}

//This function exports list to txt
char export_list(register FILE *out,t_list *list) //it shopuld be FILE
{
unsigned int _i;
if (!(fprintf(out,"%d\n",list->size)))
  return FALSE;
for (_i=0;_i<list->size;_i++)
  if (!(fprintf(out,"%d\n",list->list[_i])))
    return FALSE;
return TRUE;
}

//This function alloc memory for t_list
inline t_list* alloc_list(register unsigned int size)
{
register t_list* list;
extern unsigned int ylib_errno;
if (!(list=(t_list*)malloc(sizeof(t_list)+sizeof(unsigned int)*size))) { ylib_errno=YERROR_MEMORY; return FALSE; } 
list->list=(unsigned int*)((void*)list+sizeof(t_list));
list->size=size;
return list;
}

//This function alloc memory  for t_clist. Memory is not initialized!
inline t_clist* alloc_clist(register unsigned int size,register unsigned int _size)
{
register t_clist *clist;
if (!(clist=(t_clist*)malloc(sizeof(t_clist)+sizeof(t_list)*size+sizeof(unsigned int)*_size))) { ylib_errno=YERROR_MEMORY; return FALSE; }
clist->size=size;
clist->_size=_size;
clist->_items=(void*)clist+sizeof(t_clist);
clist->list=(void*)clist+sizeof(t_clist)+sizeof(unsigned int)*_size; //Move it in the end to avoid ambigious copiing in realloc
return clist;
}

//This function makes copy of complex list
t_clist *copy_clist(t_clist *srcc)
{
register unsigned int _i, _j;
t_clist *ccopy;

if (!(ccopy=alloc_clist(srcc->size,srcc->_size))) return FALSE;
else
  {
  memcpy(ccopy->_items,srcc->_items,sizeof(unsigned int)*srcc->_size);
  for (_i=_j=0; _i<ccopy->size; _i++)
    {
    ccopy->list[_i].list=&ccopy->_items[_j];
    _j+=ccopy->list[_i].size=srcc->list[_i].size;
    }
  }
return ccopy;
}


//THis function alloc objects
t_object *alloc_object(unsigned int size)
{
t_object *object;
object=(t_object*)malloc(size+sizeof(t_object));
object->object=(void*)object+sizeof(t_object);
object->size=size;
return object;
}

//This function realloc t_list.
inline t_list* realloc_list(register t_list* list,register unsigned int size)
{
extern unsigned int ylib_errno;
if (!(list=(t_list*)realloc(list,sizeof(t_list)+sizeof(unsigned int)*size))) { ylib_errno=YERROR_MEMORY; return FALSE; }
list->list=(unsigned int*)((void*)list+sizeof(t_list));
return list;
}

//This function realloc t_object
t_object* realloc_object(t_object *object,unsigned int size)
{
object=(t_object*)realloc(object,sizeof(object)+size);
object->object=(void*)object+sizeof(object);
return object;
}

//Search functions

//This function find object in row of values (any sizeof(int) datatypes). It returns (unsigned int)-1 on the absence.
inline unsigned int  find_in_row(register unsigned int item,register unsigned int size,register unsigned int *row)
{
row+=size; while (size--) if (*(--row)==item) break; return size;
}
//This function find object in row of values (any sizeof(int) datatypes) but eXcludes zero element from consideraton. It returns zero on the absence.
inline unsigned int xfind_in_row(register unsigned int item,register unsigned int size,register unsigned int *row)
{
row+=size; while (--size) if (*(--row)==item) break; return size;
}

//Fast finders. They require lists to be pre-sorted. They returns TRUE and the position if match is found and FALSE and the position to paste a new value if there is no matching.

//This function finds the value in a row of _sorted_ unsigned ints. It expects the sort order is from low to high
inline char find_in_sorted_lth_urow(unsigned int *_id,register unsigned int item,register unsigned int size,register unsigned int *urow)
{
register unsigned int start, mid, end;

if ( (!size--)||(!(urow)) ) { if (_id) *_id=0; return FALSE; } //Exception
start=0, end=size;
while (start!=end)
  {
  mid=start+(end-start)/2;
       if (urow[mid] < item) start=mid+1;
  else if (urow[mid] > item) end  =mid;
  else { if (_id) *_id=mid; return TRUE; } //Item is found in the urow
  }
//Item is not found in the urow
     if (urow[size] < item) 
       {
       if (_id)
         {
         while ( (end!=size)&&(urow[end]==urow[end+1]) ) end++; //In the case of non-unique objects in massive
         *_id=size+1;
         }
       return FALSE;
       }
else if (urow[size] > item)
       {
       if (_id)
         {
         while ( ( (end))&&(urow[end]==urow[end-1]) ) end--; //In the case of non-unique objects in massive
         *_id=end;
         }
       return FALSE;
       }
else   { if (_id) *_id=size;   return TRUE;  }
}
//This function finds the value in a row of _sorted_ unsigned ints. It expects the sort order is from low to high
inline char find_in_sorted_lth_irow(unsigned int *_id,register int item,register unsigned int size,register unsigned int *irow)
{
register unsigned int start, mid, end;
if ( (!size--)||(!(irow)) ) { if (_id) *_id=0; return FALSE; } //Exception
start=0, end=size;
while (start!=end)
  {
  mid=start+(end-start)/2;
       if (irow[mid] < item) start=mid+1;
  else if (irow[mid] > item) end  =mid;
  else { if (_id) *_id=mid; return TRUE; } //Item is found in the irow
  }
//Item is not found in the irow
     if (irow[size] < item) 
       {
       if (_id)
         {
         while ( (end!=size)&&(irow[end]==irow[end+1]) ) end++; //In the case of non-unique objects in massive
         *_id=size+1;
         }
       return FALSE;
       }
else if (irow[size] > item)
       {
       if (_id)
         {
         while ( ( (end))&&(irow[end]==irow[end-1]) ) end--; //In the case of non-unique objects in massive
         *_id=end;
         }
       return FALSE;
       }
else   { if (_id) *_id=size;   return TRUE;  }
}

//This function finds an object in the stack of objects
inline char find_in_sorted_lth_objects(unsigned int *_id,register void *item,unsigned int size,register void *objects,register size_t object_size,int (*compare_objects)(const void* object0,const void* object1))
{
register unsigned int start, mid, end;
register int _t;
if ( (!size--)||(!(objects)) ) { if (_id) *_id=0; return FALSE; } //Exception
start=0, end=size;
while (start!=end)
  {
  mid=start+(end-start)/2;
  _t=compare_objects(objects+(size_t)mid*object_size,item);
       if ( _t < 0 ) start=mid+1;
  else if ( _t > 0 ) end  =mid;
  else { if (_id) *_id=mid; return TRUE; } //Item is found in the objects
  }
//Item is not found in the objects
_t=compare_objects(objects+(size_t)end*object_size,item);
     if (_t < 0) 
       {
       if (_id) 
         {
         while ( (end!=size)&&(!(compare_objects(objects+(size_t)end*object_size,objects+(size_t)(end+1)*object_size))) ) end++; //In the case of non-unique objects in massive
         *_id=end+1;
         }
       return FALSE;
       }
else if (_t > 0) 
       {
       if (_id)
         {
         while ( ( (end))&&(!(compare_objects(objects+(size_t)end*object_size,objects+(size_t)(end-1)*object_size))) ) end--; //In the case of non-unique objects in massive
         *_id=end;
         }
       return FALSE;
       }
else   { if (_id) *_id=size;   return TRUE;  }
}

//The list wrappers

//This function find object in list massive
inline unsigned int  find_in_list(register unsigned int item,register t_list *list)
{
return find_in_row(item,list->size,list->list);
}
//This function find object in list massive but eXcldes zero element from consideration
inline unsigned int xfind_in_list(register unsigned int item,register t_list *list)
{
return xfind_in_row(item,list->size,list->list);
}

//This function finds object in sorted list
inline char find_in_sorted_lth_list(unsigned int *_id,register unsigned int item,register t_list *list)
{
return find_in_sorted_lth_urow(_id,item,list->size,list->list);
}


//This function find object in the row of lists
inline char find_in_clist(register unsigned int *_id,register unsigned int item,register t_clist *clist)
{
register unsigned int size;
size=clist->size;
while (size--)
  if ((_id[1]=find_in_row(item,clist->_size,clist->_items))!=(unsigned int)-1)
    {
    _id[0]=size;
    return TRUE;
    }
return FALSE;
}


//This funtion check is two rows has any overlap elements
//NOTE it return TRUE if there is an overlap and write two (optional) unsigned ints with id of overlaping element, a_row[i]==b_row[j]
inline char overlap_row(unsigned int *i, unsigned int a_size,unsigned int *a_row,unsigned int *j,unsigned int b_size,unsigned int *b_row)
{
register unsigned int *row, size;

if (a_size<b_size)
  {
  if ( (a_size))
    {
    a_row+=a_size;
    while (a_size--)
      {
      a_row--;
      row=b_row;
      size=b_size, row+=size;
      while (size--)
        {
        row--;
        if (*a_row==*row)
          {
          if ( (i)) *i=a_size;
          if ( (j)) *j=size;
          return TRUE;
          }
        }
      }
    }
  }
else
  {
  if ( (b_size))
    {
    b_row+=b_size;
    while (b_size--)
      {
      b_row--;
      row=a_row;
      size=a_size, row+=size;
      while (size--)
        {
        row--;
        if (*b_row==*row)
          {
          if ( (i)) *i=size;
          if ( (j)) *j=b_size;
          return TRUE;
          }
        }
      }
    }
  }
return FALSE; //no overlap
}


//This functions finds min and max values
//This function find the lowest value in a row of doubles
unsigned int find_min_d(register unsigned int size,register double *_f)
{
register unsigned int id;
if (!size)
  return FALSE;
id=0;
while (--size)
  if (_f[size]<_f[id])
    id=size;
return id;
}

//This function find the lowest value by module in a row of doubles
unsigned int find_min_absd(register unsigned int size,register double *_f)
{
register unsigned int id;
id=0;
while (--size)
  if (fabs(_f[size])<fabs(_f[id]))
    id=size;
return id;
}

//This function find next min in row of doubles
char find_next_min_d(register unsigned int *id,unsigned int size,register double *_f)
{
register unsigned int _j;
register double _min;

//Get first suitable
for (_min=_f[(*id)],(*id)++;(*id)<size;(*id)++)
  if (_f[(*id)]==_min)
    return TRUE;
for ((*id)=0;(*id)<size;(*id)++)  //Sign first suitable
  if (_f[(*id)]>_min)
    goto SEEK;
return FALSE;   //Nothing sutable found
SEEK: for (_j=(*id)+1;_j<size;_j++)
        if ((_f[_j]<_f[(*id)])&&(_f[_j]>_min))
          (*id)=_j;
return TRUE;
}

//This function find next min by module in row of doubles
char find_next_min_absd(register unsigned int *id,unsigned int size,register double *_f)
{
register unsigned int _j;
register double _min;

//Get first suitable
for (_min=fabs(_f[(*id)]),(*id)++;(*id)<size;(*id)++)
  if (fabs(_f[(*id)])==_min)
    return TRUE;
for ((*id)=0;(*id)<size;(*id)++)  //Sign first suitable
  if (fabs(_f[(*id)])>_min)
    goto SEEK;
return FALSE;   //Nothing sutable found
SEEK: for (_j=(*id)+1;_j<size;_j++)
        if ((fabs(_f[_j])<fabs(_f[(*id)]))&&(fabs(_f[_j])>_min))
          (*id)=_j;
return TRUE;
}

//This function find the lowest value in a row of ints
unsigned int find_min_i(register unsigned int size,register int *_i)
{
register unsigned int id;
if (!size)
  return FALSE;
id=0;
while (--size)
  if (_i[size]<_i[id])
    id=size;
return id;
}

//This function find next min in row of ints
char find_next_min_i(register unsigned int *id,unsigned int size,register int *_i)
{
register unsigned int _j;
register int _min;

//Get first suitable
for (_min=_i[(*id)],(*id)++;(*id)<size;(*id)++)
  if (_i[(*id)]==_min)
    return TRUE;
for ((*id)=0;(*id)<size;(*id)++)  //Sign first suitable
  if (_i[(*id)]>_min)
    goto SEEK;
return FALSE;   //Nothing sutable found
SEEK: for (_j=(*id)+1;_j<size;_j++)
        if ((_i[_j]<_i[(*id)])&&(_i[_j]>_min))
          (*id)=_j;
return TRUE;
}


//This function find the largest value in a row of doubles
unsigned int find_max_d(register unsigned int size,register double *_f)
{
register unsigned int id;
id=0;
while (--size)
  if (_f[size]>_f[id])
    id=size;
return id;
}

//This function find the largest by module value in a row of doubles
unsigned int find_max_absd(register unsigned int size,register double *_f)
{
register unsigned int id;
id=0;
while (--size)
  if (fabs(_f[size])>fabs(_f[id]))
    id=size;
return id;
}

//This function find next max in row of doubles
char find_next_max_d(register unsigned int *id,unsigned int size,register double *_f)
{
register unsigned int _j;
register double _max;

//Get first suitable
for (_max=_f[(*id)],(*id)++;(*id)<size;(*id)++)
  if (_f[(*id)]==_max)
    return TRUE;
for ((*id)=0;(*id)<size;(*id)++)  //Sign first suitable
  if (_f[(*id)]<_max)
    goto SEEK;
return FALSE;   //Nothing sutable found
SEEK: for (_j=(*id)+1;_j<size;_j++)
        if ((_f[_j]>_f[(*id)])&&(_f[_j]<_max))
          (*id)=_j;
return TRUE;
}


//This function find next max by module in row of doubles
char find_next_max_absd(register unsigned int *id,unsigned int size,register double *_f)
{
register unsigned int _j;
register double _max;

//Get first suitable
for (_max=fabs(_f[(*id)]),(*id)++;(*id)<size;(*id)++)
  if (fabs(_f[(*id)])==_max)
    return TRUE;
for ((*id)=0;(*id)<size;(*id)++)  //Sign first suitable
  if (fabs(_f[(*id)])<_max)
    goto SEEK;
return FALSE;   //Nothing sutable found
SEEK: for (_j=(*id)+1;_j<size;_j++)
        if ((fabs(_f[_j])>fabs(_f[(*id)]))&&(fabs(_f[_j])<_max))
          (*id)=_j;
return TRUE;
}

//This function find the largest value in a row of ints
unsigned int find_max_i(register unsigned int size,register int *_i)
{
register unsigned int id;
if (!size)
  return FALSE;
id=0;
while (--size)
  if (_i[size]>_i[id])
    id=size;
return id;
}

//This function find next max in row of ints
char find_next_max_i(register unsigned int *id,unsigned int size,register int *_i)
{
register unsigned int _j;
register int _max;

//Get first suitable
for (_max=_i[(*id)],(*id)++;(*id)<size;(*id)++)
  if (_i[(*id)]==_max)
    return TRUE;
for ((*id)=0;(*id)<size;(*id)++)  //Sign first suitable
  if (_i[(*id)]<_max)
    goto SEEK;
return FALSE;   //Nothing sutable found
SEEK: for (_j=(*id)+1;_j<size;_j++)
        if ((_i[_j]>_i[(*id)])&&(_i[_j]<_max))
          (*id)=_j;
return TRUE;
}

//This function add an element to list massive
inline char list_add(unsigned int item,t_list **list)
{
register t_list *vp;
if (!(vp=(t_list*)realloc_list(*list,(*list)->size+1))) return FALSE;
else *list=vp;
vp->list[vp->size++]=item;
return TRUE;
}

//This function add an list to lists massive
char clist_add(t_list *list,t_clist **clist)
{
register unsigned int _i, _j, _k;
void *vp;
_k=abs((int)list->size);
if (!(vp=realloc(*clist,sizeof(t_clist)+sizeof(t_list)*((*clist)->size+1)+sizeof(unsigned int)*((*clist)->_size+_k)))) { ylib_errno=YERROR_MEMORY; return 0x0; }
else { (*clist)=(t_clist*)vp; (*clist)->_items=vp+sizeof(t_clist); (*clist)->list=vp+sizeof(t_clist)+sizeof(unsigned int)*((*clist)->_size+_k); }
vp+=sizeof(t_clist)+sizeof(unsigned int)*(*clist)->_size;
_i=(*clist)->size;
_j=(*clist)->_size;
while(_i--) (*clist)->list[_i].list=&(*clist)->_items[_j-=((*clist)->list[_i].size=((t_list*)vp)[_i].size)]; //Copy lists
//Add new list
(*clist)->list[(*clist)->size].size=list->size;
(*clist)->list[(*clist)->size].list=&(*clist)->_items[(*clist)->_size];
(*clist)->_size+=_k; //Sync clist
while (_k--) (*clist)->list[(*clist)->size].list[_k]=list->list[_k]; //Copy new members
(*clist)->size++;  //Sync clist
return TRUE;
}


//This function delete an list to lists massive
inline void clist_del(unsigned int id,t_clist *clist)
{
register unsigned int _i, _j, _k;
register unsigned int *_pi;

clist->_size-=abs((int)clist->list[id].size), clist->size--;
for (_i=id, _pi=clist->list[_i].list; _i!=clist->size; _i++)
  {
  clist->list[_i].list=_pi, clist->list[_i].size=clist->list[_i+1].size, _k=abs((int)clist->list[_i].size), _pi+=_k;
  for (_j=0; _j<_k; _j++) clist->list[_i].list[_j]=clist->list[_i+1].list[_j];
  }
}

//This function finds summ of double row
inline double summ_in_rowf(register unsigned int num,register double *_f)
{
register double summ=0;
while (num) summ+=_f[--num];
return summ;
}

//This function finds summ of int row
inline int summ_in_rowi(unsigned int num,int *_i)
{
register int _summ=0;
while (num--) _summ+=_i[num];
return _summ;
}

// ------------------- L O G I C   P A R T ---------------------

//This fuinction performs logic operaion or  at two t_lists (Union) c=(a+b)
t_list *or_list(register t_list *a,register t_list *b)
{
register unsigned int _i, _j, _k;
register t_list *list, *_list;
extern unsigned int ylib_errno;

if (!(list=(t_list*)alloc_list(a->size+b->size))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
list->size=a->size;
memcpy(list->list,a->list,a->size*sizeof(unsigned int)); 
_j=b->size; while (_j--) { _k=b->list[_j], _i=a->size;    while (_i--) { if (a->list[_i]==_k)     goto NEXT_J;  } list->list[list->size++]=_k; NEXT_J:  ; } 
//Sync memory
if (!(_list=(t_list*)realloc_list(list,list->size))) { free(list); list=0x0; goto LABEL_MEMORY_ERROR; }
if (_list!=list) 
  {
  _list->size=list->size, _list->list=(unsigned int*)((void*)list+sizeof(t_list)); 
  memcpy(_list->list,list->list,sizeof(unsigned int)*list->size); 
  free(list); list=0x0;
  }
return _list;
}
//This function performs logic operation not at two t_list (Exclude) c=(a-b)
t_list *not_list(register t_list *a,register t_list *b)
{
register unsigned int _i, _j, _k;
register t_list *list, *_list;

if (!(list=(t_list*)alloc_list(a->size))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
list->size=0, _i=a->size; while (_i--) { _k=a->list[_i]; _j=b->size; while (_j--) { if (b->list[_j]==_k) goto NEXT_I; } list->list[list->size++]=_k; NEXT_I: ; }
//Sync memory
if (!(_list=(t_list*)realloc_list(list,list->size))) { free(list); list=0x0; goto LABEL_MEMORY_ERROR; }
if (_list!=list) 
  {
  _list->size=list->size, _list->list=(unsigned int*)((void*)list+sizeof(t_list)); 
  memcpy(_list->list,list->list,sizeof(unsigned int)*list->size); 
  free(list); list=0x0;
  }
return _list;
}
//This function performs logic operation and at two t_list (Cross) c=(a&b))
t_list *and_list(register t_list *a,register t_list *b)
{
register unsigned int _i, _j, _k;
register t_list *list, *_list;

if (!(list=(t_list*)alloc_list( (a->size<b->size) ? a->size : b->size ))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
list->size=0, _i=a->size; while (_i--) { _k=a->list[_i]; _j=b->size; while (_j--) { if (b->list[_j]==_k) { list->list[list->size++]=_k; break; } } }
//Sync memory
if (!(_list=(t_list*)realloc_list(list,list->size))) { free(list); list=0x0; goto LABEL_MEMORY_ERROR; }
if (_list!=list) 
  {
  _list->size=list->size, _list->list=(unsigned int*)((void*)list+sizeof(t_list)); 
  memcpy(_list->list,list->list,sizeof(unsigned int)*list->size); 
  free(list); list=0x0;
  }
return _list;
}

//This funtion check is two lists has any overlap elements
//NOTE it return TRUE if there is and overlap and two unsigned ints with id of overlaping element a->list[i]==b->list[j]
inline char overlap_list(register unsigned int *i,register t_list *a,register unsigned int *j,register t_list *b)
{
return overlap_row(i,a->size,a->list,j,b->size,b->list);
}






//Sorting algorithms

//This function sorts 3 ints with efficient if algorithm
inline void sort_iif3(int *i0,int *i1,int *i2)
{
register int _i;
if (*i0>*i1) 
  {
  if (*i1>*i2) ; //0,1,2
  else if (*i0>*i2)  { _i=*i1, *i1=*i2, *i2=_i; } //0,2,1
       else { _i=*i0, *i0=*i2, *i2=*i1, *i1=_i; } //2,0,1
  }
else 
  {
  if (*i0>*i2) { _i=*i0, *i0=*i1, *i1=_i; } //1,0,2
  else if (*i1>*i2) { _i=*i0, *i0=*i1, *i1=*i2, *i2=_i; } //1,2,0
       else                  { _i=*i0, *i0=*i2, *i2=_i; } //2,1,0 
  }
}
//This function sorts 3 doubles with efficient if algorithm
inline void sort_dif3(double *d0,double *d1,double *d2)
{
register double _d;
if (*d0>*d1) 
  {
  if (*d1>*d2) ; //0,1,2
  else if (*d0>*d2)  { _d=*d1, *d1=*d2, *d2=_d; } //0,2,1
       else { _d=*d0, *d0=*d2, *d2=*d1, *d1=_d; } //2,0,1
  }
else 
  {
  if (*d0>*d2) { _d=*d0, *d0=*d1, *d1=_d; } //1,0,2
  else if (*d1>*d2) { _d=*d0, *d0=*d1, *d1=*d2, *d2=_d; } //1,2,0
       else                  { _d=*d0, *d0=*d2, *d2=_d; } //2,1,0 
  }
}
//This function sorts 3 doubles-int pairs with efficient if algorithm
inline void sort_diif3(double *d0,double *d1,double *d2,int *i0,int *i1, int *i2)
{
register double _d;
register int _i;
if (*d0>*d1) 
  {
  if (*d1>*d2) ; //0,1,2
  else if (*d0>*d2)  { _d=*d1, *d1=*d2, *d2=_d, _i=*i1, *i1=*i2, *i2=_i; } //0,2,1
       else { _d=*d0, *d0=*d2, *d2=*d1, *d1=_d, _i=*i0, *i0=*i2, *i2=*i1, *i1=_i; } //2,0,1
  }
else 
  {
  if (*d0>*d2) { _d=*d0, *d0=*d1, *d1=_d, _i=*i0, *i0=*i1, *i1=_i; } //1,0,2
  else if (*d1>*d2) { _d=*d0, *d0=*d1, *d1=*d2, *d2=_d, _i=*i0, *i0=*i1, *i1=*i2, *i2=_i; } //1,2,0
       else                  { _d=*d0, *d0=*d2, *d2=_d, _i=*i0, *i0=*i2, *i2=_i; } //2,1,0 
  }
}
//This function sorts 3 doubles-double pairs with efficient if algorithm
inline void sort_ddif3(double *d0,double *d1,double *d2,double *dd0,double *dd1,double *dd2)
{
register double _d;
if (*d0>*d1) 
  {
  if (*d1>*d2) ; //0,1,2
  else if (*d0>*d2)  { _d=*d1, *d1=*d2, *d2=_d, _d=*dd1, *dd1=*dd2, *dd2=_d; } //0,2,1
       else { _d=*d0, *d0=*d2, *d2=*d1, *d1=_d, _d=*dd0, *dd0=*dd2, *dd2=*dd1, *dd1=_d; } //2,0,1
  }
else 
  {
  if (*d0>*d2) { _d=*d0, *d0=*d1, *d1=_d, _d=*dd0, *dd0=*dd1, *dd1=_d; } //1,0,2
  else if (*d1>*d2) { _d=*d0, *d0=*d1, *d1=*d2, *d2=_d, _d=*dd0, *dd0=*dd1, *dd1=*dd2, *dd2=_d; } //1,2,0
       else                  { _d=*d0, *d0=*d2, *d2=_d, _d=*dd0, *dd0=*dd2, *dd2=_d; } //2,1,0 
  }
}


//This function performs sorting of int array. 
//Sorting is performed with qsort algorithm from http://algolist.manual.ru/sort/quick_sort.php

//Sorting is performed with insertion sort algorithm
inline void i_isort(register unsigned int size,register int *i)
{
register unsigned int _i, _j, _t;

for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ((_j)&&(i[_j-1]>i[_j]))
      {
      _t=i[_j-1], i[_j-1]=i[_j], i[_j]=_t;
      _j--;
      }
  }
}
inline char i_qsort(register unsigned int size,register int *i)
{
register unsigned int _i, _j;  //stack indexes
register int pivot, _t;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack;
//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  i_isort(size,i);
else
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        i_isort(_i+1,&i[lb]);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (i[_i]>i[ppos])
             {
             _t=i[_i], i[_i]=i[ppos], i[ppos]=_t;
             if (i[ppos]>i[_j])
               {
               _t=i[_j], i[_j]=i[ppos], i[ppos]=_t;
               if (i[_i]>i[ppos])
                 {//3<2<1
                 _t=i[_i], i[_i]=i[ppos], i[ppos]=_t;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (i[ppos]>i[_j])
             {
             _t=i[_j], i[_j]=i[ppos], i[ppos]=_t;
             if (i[_i]>i[ppos])
               {//3<1<2
               _t=i[_i], i[_i]=i[ppos], i[ppos]=_t;
               }
             //else 1<3<2
             }
           //else 1<2<3

      pivot=i[ppos];
      //Do inner loop
      do{
        do _i++; while (i[_i]<pivot);
        do _j--; while (i[_j]>pivot);
        if (_i<=_j)
          {
          _t=i[_i], i[_i]=i[_j], i[_j]=_t;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            }  
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }
//Checker
_i=size; while (--_i) if (i[_i-1]>i[_i]) yprintf(YPRINTF_WARNING,"ERROR in i_qsort()!!\n");
return TRUE;
}

//This function performs sorting of unsigned int array. 
inline void u_isort(register unsigned int size,register unsigned int *u)
{
register unsigned int _i, _j, _t;

for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ((_j)&&(u[_j-1]>u[_j]))
      {
      _t=u[_j-1], u[_j-1]=u[_j], u[_j]=_t;
      _j--;  
      }
  }
}
//Note for lbstack[size] and ubstack[size] size=2+(unsigned int)(ln((float)amount)/LN_2), amount >= 0.
inline char u_qsort(register unsigned int size,register unsigned int *u)
{
register unsigned int _i, _j;  //stack indexes
register unsigned int pivot, _t;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack;
//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  u_isort(size,u);
else
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        u_isort(_i+1,&u[lb]);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (u[_i]>u[ppos])
             {
             _t=u[_i], u[_i]=u[ppos], u[ppos]=_t;
             if (u[ppos]>u[_j])
               {
               _t=u[_j], u[_j]=u[ppos], u[ppos]=_t;
               if (u[_i]>u[ppos])
                 {//3<2<1
                 _t=u[_i], u[_i]=u[ppos], u[ppos]=_t;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (u[ppos]>u[_j])
             {
             _t=u[_j], u[_j]=u[ppos], u[ppos]=_t;
             if (u[_i]>u[ppos])
               {//3<1<2
               _t=u[_i], u[_i]=u[ppos], u[ppos]=_t;
               }
             //else 1<3<2
             }
           //else 1<2<3

      pivot=u[ppos];
      //Do inner loop
      do{
        do _i++; while (u[_i]<pivot);
        do _j--; while (u[_j]>pivot);
        if (_i<=_j)
          {
          _t=u[_i], u[_i]=u[_j], u[_j]=_t;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            }  
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }
//Checker
_i=size; while (--_i) if (u[_i-1]>u[_i]) yprintf(YPRINTF_WARNING,"ERROR in u_qsort()!!\n");
return TRUE;
}

//This function sorts massives of doubles and moves the associated int massive accordingly
inline void di_isort(register unsigned int size,register double *d,register int *ii)
{
register unsigned int _i, _j;
register int _t;
register double _d;
for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ((_j)&&(d[_j-1]>d[_j]))
      {
       _d=d[_j-1],  d[_j-1]=d[_j],   d[_j]=_d;
      _t=ii[_j-1], ii[_j-1]=ii[_j], ii[_j]=_t;
      _j--;
      }
  }
}
//Note for lbstack[size] and ubstack[size] size=2+(unsigned int)(ln((float)amount)/LN_2), amount >= 0.
inline char di_qsort(register unsigned int size,register double *d,register int *ii)
{
register unsigned int _i, _j;  //stack indexes
register int _t;
register double _d, pivot;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack;
//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  di_isort(size,d,ii);
else
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        di_isort(_i+1,&d[lb],&ii[lb]);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (d[_i]>d[ppos])
             {
             _d=d[_i],   d[_i]=d[ppos],   d[ppos]=_d;
             _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
             if (d[ppos]>d[_j])
               {
               _d=d[_j],   d[_j]=d[ppos],   d[ppos]=_d;
               _t=ii[_j], ii[_j]=ii[ppos], ii[ppos]=_t;
               if (d[_i]>d[ppos])
                 {//3<2<1
                 _d=d[_i],   d[_i]=d[ppos],   d[ppos]=_d;
                 _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (d[ppos]>d[_j])
             {
             _d=d[_j],   d[_j]=d[ppos],   d[ppos]=_d;
             _t=ii[_j], ii[_j]=ii[ppos], ii[ppos]=_t;
             if (d[_i]>d[ppos])
               {//3<1<2
               _d=d[_i],   d[_i]=d[ppos],   d[ppos]=_d;
               _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
               }
             //else 1<3<2
             }
           //else 1<2<3 

      pivot=d[ppos];
      //Do inner loop
      do{
        do _i++; while (d[_i]<pivot);
        do _j--; while (d[_j]>pivot);
        if (_i<=_j)
          {
          _d=d[_i],   d[_i]=d[_j],   d[_j]=_d;
          _t=ii[_i], ii[_i]=ii[_j], ii[_j]=_t;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            }  
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }
//Checker
_i=size; while (--_i) if (d[_i-1]>d[_i]) yprintf(YPRINTF_WARNING,"ERROR in di_qsort()!!\n");
return TRUE;
}

//This function sorts massives of floats and moves the associated int massive accordingly
inline void fi_isort(register unsigned int size,register float *f,register int *ii)
{
register unsigned int _i, _j;
register int _t;
register float _f;
for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ((_j)&&(f[_j-1]>f[_j]))
      {
       _f=f[_j-1],  f[_j-1]=f[_j],   f[_j]=_f;
      _t=ii[_j-1], ii[_j-1]=ii[_j], ii[_j]=_t;
      _j--;
      }
  }
}
//Note for lbstack[size] and ubstack[size] size=2+(unsigned int)(ln((float)amount)/LN_2), amount >= 0.
inline char fi_qsort(register unsigned int size,register float *f,register int *ii)
{
register unsigned int _i, _j;  //stack indexes
register int _t;
register float _f, pivot;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack; //stacks

//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  fi_isort(size,f,ii);
else
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        fi_isort(_i+1,&f[lb],&ii[lb]);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (f[_i]>f[ppos])
             {
             _f=f[_i],   f[_i]=f[ppos],   f[ppos]=_f;
             _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
             if (f[ppos]>f[_j])
               {
               _f=f[_j],   f[_j]=f[ppos],   f[ppos]=_f;
               _t=ii[_j], ii[_j]=ii[ppos], ii[ppos]=_t;
               if (f[_i]>f[ppos])
                 {//3<2<1
                 _f=f[_i],   f[_i]=f[ppos],   f[ppos]=_f;
                 _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (f[ppos]>f[_j])
             {
             _f=f[_j],   f[_j]=f[ppos],   f[ppos]=_f;
             _t=ii[_j], ii[_j]=ii[ppos], ii[ppos]=_t;
             if (f[_i]>f[ppos])
               {//3<1<2
               _f=f[_i],   f[_i]=f[ppos],   f[ppos]=_f;
               _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
               }
             //else 1<3<2
             }
           //else 1<2<3

      pivot=f[ppos];
      //Do inner loop
      do{
        do _i++; while (f[_i]<pivot);
        do _j--; while (f[_j]>pivot);
        if (_i<=_j)
          {
          _f=f[_i],   f[_i]=f[_j],   f[_j]=_f;
          _t=ii[_i], ii[_i]=ii[_j], ii[_j]=_t;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            }  
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }

//Checker
_i=size; while (--_i) if (f[_i-1]>f[_i]) yprintf(YPRINTF_WARNING,"ERROR in fi_qsort()!!\n");
return TRUE;
}

//This function sorts massives of ints and moves the associated double massive accordingly
inline void id_isort(register unsigned int size,register int *i,register double *dd)
{
register unsigned int _i, _j;
register int _t;
register double _d;
for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ((_j)&&(i[_j-1]>i[_j]))
      {
      _t=i[_j-1], i[_j-1]=i[_j], i[_j]=_t;
      _d=dd[_j-1], dd[_j-1]=dd[_j], dd[_j]=_d;
      _j--;
      }
  }
}
//Note for lbstack[size] and ubstack[size] size=2+(unsigned int)(ln((float)amount)/LN_2), amount >= 0.
inline char id_qsort(register unsigned int size,register int *i,register double *dd)
{
register unsigned int _i, _j;  //stack indexes
register double _d;
register int pivot, _t;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack;

//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  id_isort(size,i,dd); 
else
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        id_isort(_i+1,&i[lb],&dd[lb]);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (i[_i]>i[ppos])
             {
             _t=i[_i],  i[_i]=i[ppos], i[ppos]=_t;
             _d=dd[_i], dd[_i]=dd[ppos], dd[ppos]=_d;
             if (i[ppos]>i[_j])
               {
               _t=i[_j], i[_j]=i[ppos], i[ppos]=_t;
               _d=dd[_j], dd[_j]=dd[ppos], dd[ppos]=_d;
               if (i[_i]>i[ppos])
                 {//3<2<1
                 _t=i[_i], i[_i]=i[ppos], i[ppos]=_t;
                 _d=dd[_i], dd[_i]=dd[ppos], dd[ppos]=_d;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (i[ppos]>i[_j])
             {
             _t=i[_j], i[_j]=i[ppos], i[ppos]=_t;
             _d=dd[_j], dd[_j]=dd[ppos], dd[ppos]=_d;
             if (i[_i]>i[ppos])
               {//3<1<2
               _t=i[_i], i[_i]=i[ppos], i[ppos]=_t;
               _d=dd[_i], dd[_i]=dd[ppos], dd[ppos]=_d;
               }
             //else 1<3<2
             }
           //else 1<2<3

      pivot=i[ppos];
      //Do inner loop
      do{
        do _i++; while (i[_i]<pivot);
        do _j--; while (i[_j]>pivot);
        if (_i<=_j)
          {
          _t=i[_i], i[_i]=i[_j], i[_j]=_t;
          _d=dd[_i], dd[_i]=dd[_j], dd[_j]=_d;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            }  
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }

//Checker
_i=size; while (--_i) if (i[_i-1]>i[_i]) yprintf(YPRINTF_WARNING,"ERROR in id_qsort()!!\n");
return TRUE;
}

//This function sorts massives of ints and moves the associated int massive accordingly
inline void ii_isort(register unsigned int size,register int *i,register int *ii)
{
register unsigned int _i, _j;
register int _t;
for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ((_j)&&(i[_j-1]>i[_j]))
      {
      _t=i[_j-1],   i[_j-1]=i[_j],   i[_j]=_t;
      _t=ii[_j-1], ii[_j-1]=ii[_j], ii[_j]=_t;
      _j--;
      }
  }
}
//Note for lbstack[size] and ubstack[size] size=2+(unsigned int)(ln((float)amount)/LN_2), amount >= 0.
inline char ii_qsort(register unsigned int size,register int *i,register int *ii)
{
register unsigned int _i, _j;  //stack indexes
register int pivot, _t;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack;

//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  ii_isort(size,i,ii); 
else
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        ii_isort(_i+1,&i[lb],&ii[lb]);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (i[_i]>i[ppos])
             {
             _t=i[_i],   i[_i]=i[ppos],   i[ppos]=_t;
             _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
             if (i[ppos]>i[_j])
               {
               _t=i[_j],   i[_j]=i[ppos],   i[ppos]=_t;
               _t=ii[_j], ii[_j]=ii[ppos], ii[ppos]=_t;
               if (i[_i]>i[ppos])
                 {//3<2<1
                 _t=i[_i],   i[_i]=i[ppos],   i[ppos]=_t;
                 _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (i[ppos]>i[_j])
             {
             _t=i[_j],   i[_j]=i[ppos],   i[ppos]=_t;
             _t=ii[_j], ii[_j]=ii[ppos], ii[ppos]=_t;
             if (i[_i]>i[ppos])
               {//3<1<2
               _t=i[_i],   i[_i]=i[ppos],   i[ppos]=_t;
               _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
               }
             //else 1<3<2
             }
           //else 1<2<3
      pivot=i[ppos];
      //Do inner loop
      do{
        do _i++; while (i[_i]<pivot);
        do _j--; while (i[_j]>pivot);
        if (_i<=_j)
          {
          _t=i[_i],   i[_i]=i[_j],   i[_j]=_t;
          _t=ii[_i], ii[_i]=ii[_j], ii[_j]=_t;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }

//Checker
_i=size; while (--_i) if (i[_i-1]>i[_i]) yprintf(YPRINTF_WARNING,"ERROR in ii_qsort()!!\n");
return TRUE;
}

//This function performs sorting of unsigned int array and moves associated double massive. 
inline void ud_isort(register unsigned int size,register unsigned int *u,register double *dd)
{
register unsigned int _i, _j, _t;
register double _d;
for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ( (_j)&&(u[_j-1]>u[_j]) )
      {
      _t=u[_j-1],   u[_j-1]=u[_j],   u[_j]=_t;
      _d=dd[_j-1], dd[_j-1]=dd[_j], dd[_j]=_d;
      _j--;
      }
  }
}
//Note for lbstack[size] and ubstack[size] size=2+(unsigned int)(ln((float)amount)/LN_2), amount >= 0.
inline char ud_qsort(register unsigned int size,register unsigned int *u,register double *dd)
{
register unsigned int _i, _j, _t;  //stack indexes
register double _d;
unsigned int pivot;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack;
//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  ud_isort(size,u,dd);
else 
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        ud_isort(_i+1,&u[lb],&dd[lb]);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (u[_i]>u[ppos])
             {
             _t=u[_i],   u[_i]=u[ppos],   u[ppos]=_t;
             _d=dd[_i], dd[_i]=dd[ppos], dd[ppos]=_d;
             if (u[ppos]>u[_j])
               {
               _t=u[_j],   u[_j]=u[ppos],   u[ppos]=_t;
               _d=dd[_j], dd[_j]=dd[ppos], dd[ppos]=_d;
               if (u[_i]>u[ppos])
                 {//3<2<1
                 _t=u[_i],   u[_i]=u[ppos],   u[ppos]=_t;
                 _d=dd[_i], dd[_i]=dd[ppos], dd[ppos]=_d;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (u[ppos]>u[_j])
             {
             _t=u[_j],   u[_j]=u[ppos],   u[ppos]=_t;
             _d=dd[_j], dd[_j]=dd[ppos], dd[ppos]=_d;
             if (u[_i]>u[ppos])
               {//3<1<2
               _t=u[_i],   u[_i]=u[ppos],   u[ppos]=_t;
               _d=dd[_i], dd[_i]=dd[ppos], dd[ppos]=_d;
               }
             //else 1<3<2
             }
           //else 1<2<3

      pivot=u[ppos];
      //Do inner loop
      do{
        do _i++; while (u[_i]<pivot);
        do _j--; while (u[_j]>pivot);
        if (_i<=_j)
          {
          _t=u[_i],   u[_i]=u[_j],   u[_j]=_t;
          _d=dd[_i], dd[_i]=dd[_j], dd[_j]=_d;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }

//Checker
_i=size; while (--_i) if (u[_i-1]>u[_i]) yprintf(YPRINTF_WARNING,"ERROR in ud_qsort()!!\n");
return TRUE;
}

//This function sorts massives of unsigned ints and moves the associated (unsigned) int massive accordingly
inline void ui_isort(register unsigned int size,register unsigned int *u,register int *ii)
{
register unsigned int _i, _j, _t;
for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ((_j)&&(u[_j-1]>u[_j]))
      {
      _t=u[_j-1],   u[_j-1]=u[_j],   u[_j]=_t;
      _t=ii[_j-1], ii[_j-1]=ii[_j], ii[_j]=_t;
      _j--;
      }
  }
}
//Note for lbstack[size] and ubstack[size] size=2+(unsigned int)(ln((float)amount)/LN_2), amount >= 0.
inline char ui_qsort(register unsigned int size,register unsigned int *u,register int *ii)
{
register unsigned int _i, _j, _t;  //stack indexes
int pivot;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack;
//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  ui_isort(size,u,ii); 
else
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        ui_isort(_i+1,&u[lb],&ii[lb]);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (u[_i]>u[ppos])
             {
             _t=u[_i],   u[_i]=u[ppos],   u[ppos]=_t;
             _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
             if (u[ppos]>u[_j])
               {
               _t=u[_j],   u[_j]=u[ppos],   u[ppos]=_t;
               _t=ii[_j], ii[_j]=ii[ppos], ii[ppos]=_t;
               if (u[_i]>u[ppos])
                 {//3<2<1
                 _t=u[_i],   u[_i]=u[ppos],   u[ppos]=_t;
                 _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (u[ppos]>u[_j])
             {
             _t=u[_j],   u[_j]=u[ppos],   u[ppos]=_t;
             _t=ii[_j], ii[_j]=ii[ppos], ii[ppos]=_t;
             if (u[_i]>u[ppos])
               {//3<1<2
               _t=u[_i],   u[_i]=u[ppos],   u[ppos]=_t;
               _t=ii[_i], ii[_i]=ii[ppos], ii[ppos]=_t;
               }
             //else 1<3<2
             }
           //else 1<2<3

      pivot=u[ppos];
      //Do inner loop
      do{
        do _i++; while (u[_i]<pivot);
        do _j--; while (u[_j]>pivot);
        if (_i<=_j)
          {
          _t=u[_i],   u[_i]=u[_j],   u[_j]=_t;
          _t=ii[_i], ii[_i]=ii[_j], ii[_j]=_t;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }

//Checker
_i=size; while (--_i) if (u[_i-1]>u[_i]) yprintf(YPRINTF_WARNING,"ERROR in ui_qsort()!!\n");
return TRUE;
}

//This function sorts massives of unsigned ints and moves the associated (unsigned) int massive accordingly
inline void object_id_isort(unsigned int size,size_t obj_size,void *obj,unsigned int *id,int (*compare_obj)(const void *obj0,const void *obj1))
{
register unsigned int _i, _j, _t;
for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ( (_j)&&(compare_obj(obj+obj_size*(size_t)id[_j-1],obj+obj_size*(size_t)id[_j])>0) )
    {
    _t=id[_j-1],   id[_j-1]=id[_j],   id[_j]=_t;
    _j--;
    }
  }
}
//Note for lbstack[size] and ubstack[size] size=2+(unsigned int)(ln((float)amount)/LN_2), amount >= 0.
inline char object_id_qsort(unsigned int size,size_t obj_size,void *obj,unsigned int *id,int (*compare_obj)(const void *obj0,const void *obj1))
{
register unsigned int _i, _j, _t;  //stack indexes
int pivot;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack;
//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  object_id_isort(size,obj_size,obj,id,compare_obj); 
else
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        object_id_isort(_i+1,obj_size,obj,&id[lb],compare_obj);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (compare_obj(obj+obj_size*(size_t)id[_i],obj+obj_size*(size_t)id[ppos])>0)
             {
             _t=id[_i],   id[_i]=id[ppos],   id[ppos]=_t;
             if (compare_obj(obj+obj_size*(size_t)id[ppos],obj+obj_size*(size_t)id[_j])>0)
               {
               _t=id[_j],   id[_j]=id[ppos],   id[ppos]=_t;
               if (compare_obj(obj+obj_size*(size_t)id[_i],obj+obj_size*(size_t)id[ppos])>0)
                 {//3<2<1
                 _t=id[_i],   id[_i]=id[ppos],   id[ppos]=_t;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (compare_obj(obj+obj_size*(size_t)id[ppos],obj+obj_size*(size_t)id[_j])>0)
             {
             _t=id[_j],   id[_j]=id[ppos],   id[ppos]=_t;
             if (compare_obj(obj+obj_size*(size_t)id[_i],obj+obj_size*(size_t)id[ppos])>0)
               {//3<1<2
               _t=id[_i],   id[_i]=id[ppos],   id[ppos]=_t;
               }
             //else 1<3<2
             }
           //else 1<2<3

      pivot=id[ppos];
      //Do inner loop
      do{
        do _i++; while (compare_obj(obj+obj_size*(size_t)id[_i],obj+obj_size*(size_t)id[pivot])<0);
        do _j--; while (compare_obj(obj+obj_size*(size_t)id[_j],obj+obj_size*(size_t)id[pivot])>0);
        if (_i<=_j)
          {
          _t=id[_i],   id[_i]=id[_j],   id[_j]=_t;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }

//Checker
_i=size; while (--_i) if (compare_obj(obj+obj_size*(size_t)id[_i-1],obj+obj_size*(size_t)id[_i])>0) yprintf(YPRINTF_WARNING,"ERROR in object_id_qsort()!!\n");
return TRUE;
}

//Sorting of objects array moved by pointers
inline void v_isort(register unsigned int size,void **objects,int (*compare_objects)(const void* object0,const void* object1))
{
register unsigned int _i, _j;
register void *object;

for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ( (_j)&&(compare_objects(objects[_j-1],objects[_j])>0) )
    {
    object=objects[_j-1], objects[_j-1]=objects[_j], objects[_j]=object;
    _j--;
    }
  }
}
inline char v_qsort(register unsigned int size,void **objects,int (*compare_objects)(const void* object0,const void* object1))
{
register unsigned int _i, _j;  //stack indexes
register void *object, *pivot;
unsigned int lb, ub, stackpos, ppos;  // stack positions
unsigned int *up, *lbstack, *ubstack;
//Make it numerically stable
if (size<=ISORT_ALGORITHM_CUTOFF)
  v_isort(size,objects,compare_objects);
else
  {//Do qsort
  if (!(lbstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_0: ylib_errno=YERROR_MEMORY;   return FALSE; }
  if (!(ubstack=(unsigned int*)malloc(sizeof(unsigned int)*0xFF))) { LABEL_MEMORY_ERROR_1: free(lbstack); goto LABEL_MEMORY_ERROR_0; }
  stackpos=0x1;
  lbstack[stackpos]=0x0;
  ubstack[stackpos]=size-0x1;
  do{
    // Get lb and ub limits of currennt massive from stack
    lb=lbstack[stackpos];
    ub=ubstack[stackpos];
    stackpos--;
    while (lb<ub)
      {
      if ((_i=ub-lb)<=ISORT_ALGORITHM_CUTOFF)
        {//Apply isort when the massive is small enough
        v_isort(_i+1,&objects[lb],compare_objects);
        break;
        }
      _i=lb, _j=ub;
      // Stage 1. Devide accordingly to pivot
      ppos=(lb+ub)>>0x1;
      //sort borders
           if (compare_objects(objects[_i],objects[ppos])>0)
             {
             object=objects[_i], objects[_i]=objects[ppos], objects[ppos]=object;
             if (compare_objects(objects[ppos],objects[_j])>0)
               {
               object=objects[_j], objects[_j]=objects[ppos], objects[ppos]=object;
               if (compare_objects(objects[_i],objects[ppos])>0)
                 {//3<2<1
                 object=objects[_i], objects[_i]=objects[ppos], objects[ppos]=object;
                 }
               //else 2<3<1
               }
             //else 2<1<3
             }
      else if (compare_objects(objects[ppos],objects[_j])>0)
             {
             object=objects[_j], objects[_j]=objects[ppos], objects[ppos]=object;
             if (objects[_i]>objects[ppos])
               {//3<1<2
               object=objects[_i], objects[_i]=objects[ppos], objects[ppos]=object;
               }
             //else 1<3<2
             }
           //else 1<2<3

      pivot=objects[ppos];
      //Do inner loop
      do{
        do _i++; while (compare_objects(objects[_i],pivot)<0);
        do _j--; while (compare_objects(objects[_j],pivot)>0);
        if (_i<=_j)
          {
          object=objects[_i], objects[_i]=objects[_j], objects[_j]=object;
          }
        }while (_i<_j);

      //Now i point to the begin of right submassive
      // j - end of the left submassive, lb > j > i > ub.
      // It is possible thar i and j are out of massive

      //Stage 2 and 3.Move larger part into stack and move lb, ub as well
      if (_i<ppos)
        {     // right part is larger
        if (_i<ub)
          {     //  if it hase more than 1 element it is needs to
          if (!(++stackpos%0xFF))   //  sort, move request into stack
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) { LABEL_MEMORY_ERROR_2: free(ubstack); goto LABEL_MEMORY_ERROR_1; } else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            }  
          lbstack[stackpos]=_i;
          ubstack[stackpos]=ub;
          }
        ub=_j; // Next deviding iteration will deals with left part
        }
      else
        {      // left part larger
        if (_j>lb)
          {
          if (!(++stackpos%0xFF))
            {
            if (!(up=(unsigned int*)realloc(lbstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else lbstack=up;
            if (!(up=(unsigned int*)realloc(ubstack,sizeof(unsigned int)*(stackpos+0xFF)))) goto LABEL_MEMORY_ERROR_2; else ubstack=up;
            } 
          lbstack[stackpos]=lb;
          ubstack[stackpos]=_j;
          }
        lb=_i;
        }
      }
    }while ((stackpos));     // while stack request exists
  free(lbstack);
  free(ubstack);
  }
//Checker
_i=size; while (--_i) if (compare_objects(objects[_i-1],objects[_i])>0) yprintf(YPRINTF_WARNING,"ERROR in v_qsort()!!\n");
return TRUE;
}

//Sorting of objects array and then move physical items
inline void w_isort(unsigned int size,size_t obj_size,void *obj,int (*compare_obj)(const void* obj0,const void* obj1),void (*swap_obj)(const void* obj0,const void* obj1))
{
register unsigned int _i, _j;
for (_i=1; _i<size; _i++)
  {
  _j=_i;
  while ( (_j)&&(compare_obj(obj+obj_size*(size_t)(_j-1),obj+obj_size*(size_t)(_j))>0) )
    {
    swap_obj(obj+obj_size*(size_t)(_j-1),obj+obj_size*(size_t)(_j));
    _j--;
    }
  }
}
//We are using id-associated sort to optimize amount of real objects translocations
inline char w_qsort(unsigned int size,size_t obj_size,void *obj,int (*compare_obj)(const void* obj0,const void* obj1),void (*swap_obj)(const void* obj0,const void* obj1))
{
register unsigned int _i;
unsigned int *id, *items, *poses;
if (!(id=(unsigned int*)malloc(sizeof(unsigned int)*size*3))) { ylib_errno=YERROR_MEMORY; return FALSE; } else { items=id+size, poses=items+size; }
id+=size, items+=size, poses+=size, _i=size; while (_i--) { --id, *id=_i, --items, *items=_i, --poses, *poses=_i; } //Enumerate translocation massives
if (!(object_id_qsort(size,obj_size,obj,id,compare_obj))) { free(id); return FALSE; }
while (--size) //item=id[size], pose=size
  {
  if (id[size]!=items[size])
    {
    swap_obj(obj+obj_size*(size_t)size,obj+obj_size*(size_t)poses[id[size]]);
    items[poses[id[size]]]=items[size], poses[items[size]]=poses[id[size]];
//    items[size]=id[size], poses[id[size]]=size;
    }
  }
free(id); return TRUE;
}



