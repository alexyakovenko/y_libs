#include "y_delaunay.h"
#include <stdlib.h>
#include <string.h>

unsigned char MASK[0x8]={1,2,4,8,16,32,64,128};
//This script update the source neighbours
#define REPAIR_NEIGHBOUR_LIST(_id,n_id,vx_id,th_id) {                                                                                          \
                                                    if (d->nth[th_id][vx_id]!=0xFFFFFFFF)                                                      \
                                                      {                                                                                        \
                                                      if ((_id=find_in_row(n_id,0x4,d->nth[d->nth[th_id][vx_id]]))==(unsigned int)-1)          \
                                                        {                                                                                      \
                                                        printf("ERROR cant find %d neighbour in %d tetrahedron\n",th_id,d->nth[th_id][vx_id]); \
                                                        return FALSE;                                                                          \
                                                        }                                                                                      \
                                                      else                                                                                     \
                                                        d->nth[d->nth[th_id][vx_id]][_id]=th_id;                                               \
                                                      }                                                                                        \
                                                    }

//Triangulation part

typedef struct{
              unsigned int   th[0x2];   // The ids of tetrahedrons that forms common facet
              unsigned char sth[0x2];   // The status of tetrahedrond (0 - deleted, 1 - active)
              }t_fstack;                // Facets stack structure



//This function uploads Delaunay
t_delaunay *read_delaunay(FILE *in)
{
unsigned int size,size_th;
t_delaunay *d;
if ( (fread(&size,sizeof(unsigned int),0x1,in)!=0x1)||(size!=Y_MAGIC)||  //Check YMAGIC
     (fread(&size,sizeof(unsigned int),0x1,in)!=0x1)                 ||  //Read points number
     (fread(&size_th,sizeof(unsigned int),0x1,in)!=0x1)              ||  //Read tetrahedrons number
     (!(d=alloc_delaunay(size+0x4,size_th)))                          ) return FALSE; //Alloc memory
else { d->size=size, d->size_th=size_th; }
if ( (fread(d->lvec,sizeof(t_lvec),d->size,in)!=d->size)               ||   //Read point`s lvecs
     (fread( d->th,sizeof(unsigned int)*0x4,d->size_th,in)!=d->size_th)||   //Read tetrahedrons
     (fread(d->nth,sizeof(unsigned int)*0x4,d->size_th,in)!=d->size_th) )   //Read tetrahedrons neighbors
  {
  free_delaunay(d);
  return FALSE;
  }
return d;
}


//This function saves Delaunay
char write_delaunay(FILE *out,t_delaunay *d)
{
unsigned int i=Y_MAGIC;
return ( (fwrite(&i,sizeof(unsigned int),0x1,out)==0x1)                       &&   //write YMAGIC
         (fwrite(&d->size,sizeof(unsigned int),0x1,out)==0x1)                 &&   //Write points number
         (fwrite(&d->size_th,sizeof(unsigned int),0x1,out)==0x1)              &&   //Wtrite tetrahedrons number
         (fwrite( d->lvec,sizeof(t_lvec),d->size,out)==d->size)               &&   //Write lifted points
         (fwrite( d->th,sizeof(unsigned int)*0x4,d->size_th,out)==d->size_th) &&   //Write tetrahedrons
         (fwrite( d->nth,sizeof(unsigned int)*0x4,d->size_th,out)==d->size_th) );  //Write tetrahedrons neighbors
}

// ------------------------------- T R A D E    F U N C T I O N S ---------------------------------


//Perform trade from one to four tetrahedrons
char trade_one_four (unsigned int p,unsigned int aid,unsigned int *fsize,t_fstack **fstack,t_delaunay *d)
{
register unsigned int _id;

//Prepare the memory
d->th=(unsigned int (*)[0x4])realloc(d->th,sizeof(unsigned int)*4*(d->size_th+=4));
d->nth=(unsigned int (*)[0x4])realloc(d->nth,sizeof(unsigned int)*4*(d->size_th));
(*fstack)=(t_fstack*)realloc((*fstack),sizeof(t_fstack)*((*fsize)+4)); //four link facets will be added

//Trade 1->4
//Create new tetrahedrons
//abdp
d->th[d->size_th-0x4][0x0]=d->th[aid][0x0];   // a
d->th[d->size_th-0x4][0x1]=d->th[aid][0x1];   // b
d->th[d->size_th-0x4][0x2]=d->th[aid][0x3];   // d
d->th[d->size_th-0x4][0x3]=p;                 // p
d->nth[d->size_th-0x4][0x0]=d->size_th-0x2;   // Facet bdp : bcdp
d->nth[d->size_th-0x4][0x1]=d->size_th-0x3;   // Facet adp : acdp
d->nth[d->size_th-0x4][0x2]=d->size_th-0x1;   // Facet abp : abcp
d->nth[d->size_th-0x4][0x3]=d->nth[aid][0x2]; // Facet abd : opposite source 'c'
REPAIR_NEIGHBOUR_LIST(_id,aid,0x3,d->size_th-0x4);
//acdp
d->th[d->size_th-0x3][0x0]=d->th[aid][0x0];   // a
d->th[d->size_th-0x3][0x1]=d->th[aid][0x2];   // c
d->th[d->size_th-0x3][0x2]=d->th[aid][0x3];   // d
d->th[d->size_th-0x3][0x3]=p;                 // p
d->nth[d->size_th-0x3][0x0]=d->size_th-0x2;   // Facet cdp : bcdp
d->nth[d->size_th-0x3][0x1]=d->size_th-0x4;   // Facet adp : abdp
d->nth[d->size_th-0x3][0x2]=d->size_th-0x1;   // Facet acp : abcp
d->nth[d->size_th-0x3][0x3]=d->nth[aid][0x1]; // Facet acd
REPAIR_NEIGHBOUR_LIST(_id,aid,0x3,d->size_th-0x3);
//bcdp
d->th[d->size_th-0x2][0x0]=d->th[aid][0x1];   // b
d->th[d->size_th-0x2][0x1]=d->th[aid][0x2];   // c
d->th[d->size_th-0x2][0x2]=d->th[aid][0x3];   // d
d->th[d->size_th-0x2][0x3]=p;                 // p
d->nth[d->size_th-0x2][0x0]=d->size_th-0x3;   // Facet cdp : acdp
d->nth[d->size_th-0x2][0x1]=d->size_th-0x4;   // Facet bdp : abdp
d->nth[d->size_th-0x2][0x2]=d->size_th-0x1;   // Facet bcp : abcp
d->nth[d->size_th-0x2][0x3]=d->nth[aid][0x0]; // Facet bcd
REPAIR_NEIGHBOUR_LIST(_id,aid,0x3,d->size_th-0x2);
//abcp
d->th[d->size_th-0x1][0x0]=d->th[aid][0x0];   // a
d->th[d->size_th-0x1][0x1]=d->th[aid][0x1];   // b
d->th[d->size_th-0x1][0x2]=d->th[aid][0x2];   // c
d->th[d->size_th-0x1][0x3]=p;                 // p
d->nth[d->size_th-0x1][0x0]=d->size_th-0x2;   // Facet bcp : bcdp
d->nth[d->size_th-0x1][0x1]=d->size_th-0x3;   // Facet acp : acdp
d->nth[d->size_th-0x1][0x2]=d->size_th-0x4;   // Facet abp : abdp
d->nth[d->size_th-0x1][0x3]=d->nth[aid][0x3]; // Facet abc
REPAIR_NEIGHBOUR_LIST(_id,aid,0x3,d->size_th-0x1);
//Remove tetrahedron from triangulation
memset(d->th[aid],0xFF,sizeof(unsigned int)*0x4);
//Update fstack. Add facets abc, abd, acd and bcd
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x4;                     //abdp
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x4][0x3];        //Facet abd
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x3;                     //acdp
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x3][0x3];        //Facet acd
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x2;                     //abcp
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x2][0x3];        //Facet abc
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x1;                     //bcdp
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x1][0x3];        //Facet bcd
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;

return TRUE;
}

//Perform trade from four to one tetrahedron. Note, tetrahedrons aid and bid sugested to be aligned previously. The common point is the first one. The cid = abop and did = acop.
char trade_four_one (unsigned int aid,unsigned int bid,unsigned int cid,unsigned int did,unsigned int *fsize,t_fstack **fstack,t_delaunay *d)
{
register unsigned int _id;

//Prepare memory
d->th=(unsigned int (*)[0x4])realloc(d->th, sizeof(unsigned int)*4*(d->size_th+=1));
d->nth=(unsigned int (*)[0x4])realloc(d->nth, sizeof(unsigned int)*4*d->size_th);
(*fstack)=(t_fstack*)realloc((*fstack),sizeof(t_fstack)*((*fsize)+1));//only one new link facets will be added
//Create tetrahedrons
d->th[d->size_th-0x1][0x0]=d->th[aid][0x1]; // b
d->th[d->size_th-0x1][0x1]=d->th[aid][0x2]; // c
d->th[d->size_th-0x1][0x2]=d->th[bid][0x3]; // o
d->th[d->size_th-0x1][0x3]=d->th[aid][0x3]; // p
_id=find_in_row(d->th[bid][0x0],0x4,d->th[did]);
d->nth[d->size_th-0x1][0x0]=d->nth[did][_id]; // Facet cop
_id=find_in_row(d->th[aid][0x0],0x4,d->th[cid]);
d->nth[d->size_th-0x1][0x1]=d->nth[cid][_id]; // Facet bop
d->nth[d->size_th-0x1][0x2]=d->nth[aid][0x0]; // Facet bcp
d->nth[d->size_th-0x1][0x3]=d->nth[bid][0x0]; // Facet bco
REPAIR_NEIGHBOUR_LIST(_id,did,0x0,d->size_th-0x1);
REPAIR_NEIGHBOUR_LIST(_id,cid,0x1,d->size_th-0x1);
REPAIR_NEIGHBOUR_LIST(_id,aid,0x2,d->size_th-0x1);
REPAIR_NEIGHBOUR_LIST(_id,bid,0x3,d->size_th-0x1);
//Remove tetrahedrons from triangulation
memset(d->th[aid],0xFF,sizeof(unsigned int)*0x4);
memset(d->th[bid],0xFF,sizeof(unsigned int)*0x4);
memset(d->th[cid],0xFF,sizeof(unsigned int)*0x4);
memset(d->th[did],0xFF,sizeof(unsigned int)*0x4);
//Mark facets in the stack
for (_id=0;_id<(*fsize);_id++)
  {
  if ( ((*fstack)[_id].th[0x0]==aid)||((*fstack)[_id].th[0x0]==bid)||((*fstack)[_id].th[0x0]==cid)||((*fstack)[_id].th[0x0]==did) )
    (*fstack)[_id].sth[0x0]=FALSE;
  if ( ((*fstack)[_id].th[0x1]==aid)||((*fstack)[_id].th[0x1]==bid)||((*fstack)[_id].th[0x1]==cid)||((*fstack)[_id].th[0x1]==did) )
    (*fstack)[_id].sth[0x1]=FALSE;
  }
//Update the stack elements (as we have just one facet added we can save the memory allocation)
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x1;                         //bcop
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x1][0x3];            //Facet bco
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;
return TRUE;
}

//Trade from two to three tetrahedrons. Note, tetrahedrons aid and bid sugested to be aligned previously.
char trade_two_three(unsigned int aid, unsigned int bid,unsigned int *fsize,t_fstack **fstack,t_delaunay *d)
{
register unsigned int _id;
//Prepare the memory
d->th=(unsigned int (*)[0x4])realloc(d->th, sizeof(unsigned int)*4*(d->size_th+=3));
d->nth=(unsigned int (*)[0x4])realloc(d->nth, sizeof(unsigned int)*4*d->size_th);
(*fstack)=(t_fstack*)realloc((*fstack),sizeof(t_fstack)*(*fsize+3));  //three link facets will be added

//We do not need to order aid and bid tetrahedrons because of they was already aligned in tradebility and regularity
//Trade 2->3
//Build one tetrahedron and relpace the  two
//abpo
d->th[d->size_th-0x3][0x0]=d->th[aid][0x0];     // a
d->th[d->size_th-0x3][0x1]=d->th[aid][0x1];     // b
d->th[d->size_th-0x3][0x2]=d->th[aid][0x3];     // p
d->th[d->size_th-0x3][0x3]=d->th[bid][0x3];     // o
d->nth[d->size_th-0x3][0x0]=d->size_th-0x2;     // Facet bpo : bcpo
d->nth[d->size_th-0x3][0x1]=d->size_th-0x1;     // Facet apo : acpo
d->nth[d->size_th-0x3][0x2]=d->nth[bid][0x2];   // Facet abo
d->nth[d->size_th-0x3][0x3]=d->nth[aid][0x2];   // Facet abp
REPAIR_NEIGHBOUR_LIST(_id,bid,0x2,d->size_th-0x3);
REPAIR_NEIGHBOUR_LIST(_id,aid,0x3,d->size_th-0x3);
//bcpo
d->th[d->size_th-0x2][0x0]=d->th[aid][0x1];     // b
d->th[d->size_th-0x2][0x1]=d->th[aid][0x2];     // c
d->th[d->size_th-0x2][0x2]=d->th[aid][0x3];     // p
d->th[d->size_th-0x2][0x3]=d->th[bid][0x3];     // o
d->nth[d->size_th-0x2][0x0]=d->size_th-0x1;     // Facet cpo : acpo
d->nth[d->size_th-0x2][0x1]=d->size_th-0x3;     // Facet bpo : abpo
d->nth[d->size_th-0x2][0x2]=d->nth[bid][0x0];   // Facet bco
d->nth[d->size_th-0x2][0x3]=d->nth[aid][0x0];   // Facet bcp
REPAIR_NEIGHBOUR_LIST(_id,bid,0x2,d->size_th-0x2);
REPAIR_NEIGHBOUR_LIST(_id,aid,0x3,d->size_th-0x2);
//acpo
d->th[d->size_th-0x1][0x0]=d->th[aid][0x0];     // a
d->th[d->size_th-0x1][0x1]=d->th[aid][0x2];     // c
d->th[d->size_th-0x1][0x2]=d->th[aid][0x3];     // p
d->th[d->size_th-0x1][0x3]=d->th[bid][0x3];     // o
d->nth[d->size_th-0x1][0x0]=d->size_th-0x2;     // Facet cpo : bcpo
d->nth[d->size_th-0x1][0x1]=d->size_th-0x3;     // Facet apo : abpo
d->nth[d->size_th-0x1][0x2]=d->nth[bid][0x1];   // Facet aco
d->nth[d->size_th-0x1][0x3]=d->nth[aid][0x1];   // Facet acp
REPAIR_NEIGHBOUR_LIST(_id,bid,0x2,d->size_th-0x1);
REPAIR_NEIGHBOUR_LIST(_id,aid,0x3,d->size_th-0x1);
//Remove tetrahedrons from triangulation
memset(d->th[aid],0xFF,sizeof(unsigned int)*0x4);
memset(d->th[bid],0xFF,sizeof(unsigned int)*0x4);
//Mark facets in the stack
for (_id=0;_id<(*fsize);_id++)
  {
  if ( ((*fstack)[_id].th[0x0]==aid)||((*fstack)[_id].th[0x0]==bid) )
    (*fstack)[_id].sth[0x0]=FALSE;
  if ( ((*fstack)[_id].th[0x1]==aid)||((*fstack)[_id].th[0x1]==bid) )
    (*fstack)[_id].sth[0x1]=FALSE;
  }
//Update fstack
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x3;                     //abpo
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x3][0x2];        //Facet abo
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x2;                         //bcpo
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x2][0x2];            //Facet bco
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x1;                         //acpo
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x1][0x2];            //Facet aco
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;
return TRUE;
}

//Perform trade from three to two tetrahedrons. Note, tetrahedrons aid and bid sugested to be aligned previously. The common two points are the first ones.
char trade_three_two (unsigned int aid,unsigned int bid,unsigned int cid,unsigned int *fsize,t_fstack **fstack,t_delaunay *d)
{
register unsigned int _id;

//Prepare the memory. Trade 3->2.
d->th=(unsigned int (*)[0x4])realloc(d->th,sizeof(unsigned int)*0x4*(d->size_th+=2));
d->nth=(unsigned int (*)[0x4])realloc(d->nth,sizeof(unsigned int)*0x4*d->size_th);
(*fstack)=(t_fstack*)realloc((*fstack),sizeof(t_fstack)*(*fsize+2)); //two new link facets will be added

//Create a new tetrahedron. Note abc is a common facet of abcp and abco
//abop
d->th[d->size_th-0x2][0x0]=d->th[aid][0x0];   // a
d->th[d->size_th-0x2][0x1]=d->th[aid][0x2];   // c
d->th[d->size_th-0x2][0x2]=d->th[bid][0x3];   // o
d->th[d->size_th-0x2][0x3]=d->th[aid][0x3];   // p
d->nth[d->size_th-0x2][0x0]=d->size_th-0x1;   // Facet cop : bcpo
_id=find_in_row(d->th[aid][0x1],0x4,d->th[cid]);
d->nth[d->size_th-0x2][0x1]=d->nth[cid][_id];   // Facet aop
d->nth[d->size_th-0x2][0x2]=d->nth[aid][0x1];   // Facet acp
d->nth[d->size_th-0x2][0x3]=d->nth[bid][0x1];   // Facet aco
REPAIR_NEIGHBOUR_LIST(_id,cid,0x1,d->size_th-0x2);
REPAIR_NEIGHBOUR_LIST(_id,aid,0x2,d->size_th-0x2);
REPAIR_NEIGHBOUR_LIST(_id,bid,0x3,d->size_th-0x2);
//bcop
d->th[d->size_th-0x1][0x0]=d->th[bid][0x1];   // b
d->th[d->size_th-0x1][0x1]=d->th[bid][0x2];   // c
d->th[d->size_th-0x1][0x2]=d->th[bid][0x3];   // o
d->th[d->size_th-0x1][0x3]=d->th[aid][0x3];   // p
d->nth[d->size_th-0x1][0x0]=d->size_th-0x2;   // Facet cop : acpo
_id=find_in_row(d->th[bid][0x0],0x4,d->th[cid]);
d->nth[d->size_th-0x1][0x1]=d->nth[cid][_id];   // Facet bop
d->nth[d->size_th-0x1][0x2]=d->nth[aid][0x0];   // Facet bcp
d->nth[d->size_th-0x1][0x3]=d->nth[bid][0x0];   // Facet bco
REPAIR_NEIGHBOUR_LIST(_id,cid,0x1,d->size_th-0x1);
REPAIR_NEIGHBOUR_LIST(_id,aid,0x2,d->size_th-0x1);
REPAIR_NEIGHBOUR_LIST(_id,bid,0x3,d->size_th-0x1);
//Remove tetrahedrons from triangulation
memset(d->th[aid],0xFF,sizeof(unsigned int)*0x4);
memset(d->th[bid],0xFF,sizeof(unsigned int)*0x4);
memset(d->th[cid],0xFF,sizeof(unsigned int)*0x4);
//Mark facets in the stack
for (_id=0;_id<(*fsize);_id++)
  {
  if ( ((*fstack)[_id].th[0x0]==aid)||((*fstack)[_id].th[0x0]==bid)||((*fstack)[_id].th[0x0]==cid) )
    (*fstack)[_id].sth[0x0]=FALSE;
  if ( ((*fstack)[_id].th[0x1]==aid)||((*fstack)[_id].th[0x1]==bid)||((*fstack)[_id].th[0x1]==cid) )
    (*fstack)[_id].sth[0x1]=FALSE;
  }
//Update fstack
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x2;                     //acop
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x2][0x3];        //Facet aco
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;
(*fstack)[(*fsize)].th[0x0]=d->size_th-0x1;                     //bcop
(*fstack)[(*fsize)].th[0x1]=d->nth[d->size_th-0x1][0x3];        //Facet bco
(*fstack)[(*fsize)].sth[0x0]=(*fstack)[(*fsize)].sth[0x1]=TRUE; //Setup active facet
(*fsize)++;
return TRUE;
}


//------------------------------------ T R I A N G U L A T I O N    S E A R C H    P A R T ------------------------------------


//This function looks for the tetrahedron that contain current point in and return its id (direction algorithm is applied)
char ffind_tetrahedron(t_lvec *p,unsigned int *th,t_delaunay *d)
{
double *_p[5];
char d4,f[4];
unsigned int _i;
//Climb throught thetrahedrons
f[0]=f[1]=f[2]=f[3]=TRUE;
*th=0;
do{
  //Calculate original tetrahedron chirality
  if ((bjsort4(&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][1]],&d->lvec[d->th[*th][2]],&d->lvec[d->th[*th][3]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)) //Check abcd
    d4=!sign_ldet4(_p);
  else
    d4= sign_ldet4(_p);
  //Check all facets keeping clockwise rule
  if ( ( f[0]=(((bjsort4(p,&d->lvec[d->th[*th][1]],&d->lvec[d->th[*th][3]],&d->lvec[d->th[*th][2]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)!=sign_ldet4(_p))!=d4) ) && //Check pbdc
       ( f[1]=(((bjsort4(p,&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][2]],&d->lvec[d->th[*th][3]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)!=sign_ldet4(_p))!=d4) ) && //Check pacd
       ( f[2]=(((bjsort4(p,&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][3]],&d->lvec[d->th[*th][1]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)!=sign_ldet4(_p))!=d4) ) && //Check padb
       ( f[3]=(((bjsort4(p,&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][1]],&d->lvec[d->th[*th][2]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)!=sign_ldet4(_p))!=d4) )  ) //Check pabc
    //Check redundance of point (is it has any uncovered space inside tetrahedron)
    return ((bjsort5(p,&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][1]],&d->lvec[d->th[*th][2]],&d->lvec[d->th[*th][3]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3],(void**)&_p[4])%2)!=sign_ldet5(_p))==d4;
  //get direction
  for (_i=0;_i<3;_i++) if (!f[_i]) { f[_i]=TRUE; break; }
  *th=d->nth[*th][_i];
  }while (*th!=0xFFFFFFFF);
error_exit("ERROR. Superstructure inscription failed in ffind_tetrahedron(). Check.code.\n");
return FALSE;
}



//This function looks for the tetrahedron that contain current point in and return its id
char find_tetrahedron(unsigned int p,unsigned int *th,t_delaunay *d)
{
char d4;
double *_p[0x5];
//Climb throught thetrahedrons
for ((*th)=0;(*th)<d->size_th;(*th)++)
  {
  //Calculate original tetrahedron chirality
  if ((bjsort4(&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][1]],&d->lvec[d->th[*th][2]],&d->lvec[d->th[*th][3]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)) //Check abcd
    d4=!sign_ldet4(_p);
  else
    d4= sign_ldet4(_p);
  //Check all facets keeping clockwise rule
  if ( ( ((bjsort4(&d->lvec[p],&d->lvec[d->th[*th][1]],&d->lvec[d->th[*th][3]],&d->lvec[d->th[*th][2]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)!=sign_ldet4(_p))!=d4 ) &&     //Check pbdc
       ( ((bjsort4(&d->lvec[p],&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][2]],&d->lvec[d->th[*th][3]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)!=sign_ldet4(_p))!=d4 ) &&     //Check pacd
       ( ((bjsort4(&d->lvec[p],&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][3]],&d->lvec[d->th[*th][1]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)!=sign_ldet4(_p))!=d4 ) &&     //Check padb
       ( ((bjsort4(&d->lvec[p],&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][1]],&d->lvec[d->th[*th][2]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3])%2)!=sign_ldet4(_p))!=d4 )  )     //Check pabc
    //Check redundance of point (is it has any uncovered space inside tetrahedron)
    return ((bjsort5(&d->lvec[p],&d->lvec[d->th[*th][0]],&d->lvec[d->th[*th][1]],&d->lvec[d->th[*th][2]],&d->lvec[d->th[*th][3]],(void**)&_p[0],(void**)&_p[1],(void**)&_p[2],(void**)&_p[3],(void**)&_p[4])%2)!=sign_ldet5(_p))==d4;
  }
printf("ERROR. Superstructure inscription failed. Check.code.\n");
exit(0);
}


//This function extract two tetrahedrons from the facet stack description
//It return FALSE if facet is a convec hull
char get_tetrahedrons(t_fstack *fstack,t_delaunay *d)
{
unsigned int _i,_j;

//Align two tetrahedrons first
for (_i=0;_i<0x3;_i++)
  if (d->th[fstack->th[0x0]][_i]!=d->th[fstack->th[0x1]][_i])
    {
    for (_j=_i+1;_j<0x4;_j++)
      if (d->th[fstack->th[0x0]][_i]==d->th[fstack->th[0x1]][_j])
        {
        memexchange(&d->th [fstack->th[0x1]][_i],&d->th [fstack->th[0x1]][_j],sizeof(unsigned int));  //Copy vertex
        memexchange(&d->nth[fstack->th[0x1]][_i],&d->nth[fstack->th[0x1]][_j],sizeof(unsigned int));  //Copy neighbours
        goto NEXT_I;
        }
    memexchange(&d->th [fstack->th[0x0]][_i],&d->th [fstack->th[0x0]][0x3],sizeof(unsigned int));  //Copy vertex
    memexchange(&d->nth[fstack->th[0x0]][_i],&d->nth[fstack->th[0x0]][0x3],sizeof(unsigned int));  //Copy neighbours
    _i--;
    NEXT_I:;
    }
//Should we do anything?
return ( ((fstack->sth[0x0])||((fstack->th[0x0]=d->nth[fstack->th[0x1]][0x3])!=0xFFFFFFFF))&&((fstack->sth[0x1])||((fstack->th[0x1]=d->nth[fstack->th[0x0]][0x3])!=0xFFFFFFFF)) );
}


//This function chceck the delaunay criteria. Note, tetrahedrons aid and bid suggested to be already aligned!
char check_regularity(unsigned int aid,unsigned int bid,t_delaunay *d)
{
char d4,d5;
double *p[5];
//Calculate det for 4 matrix
if (bjsort4(&d->lvec[d->th[aid][0]],&d->lvec[d->th[aid][1]],&d->lvec[d->th[aid][2]],&d->lvec[d->th[aid][3]],(void**)&p[0],(void**)&p[1],(void**)&p[2],(void**)&p[3])%2)
  d4=!sign_ldet4(p);
else
  d4= sign_ldet4(p);
if (bjsort5(&d->lvec[d->th[aid][0]],&d->lvec[d->th[aid][1]],&d->lvec[d->th[aid][2]],&d->lvec[d->th[aid][3]],&d->lvec[d->th[bid][3]],(void**)&p[0],(void**)&p[1],(void**)&p[2],(void**)&p[3],(void**)&p[4])%2)
  d5=!sign_ldet5(p);
else
  d5= sign_ldet5(p);
return d4!=d5;
}

//This function check two tetrahedrons for their fliability. Note this two tetrahedrons should have aligned sequences already!!!!
//It is equile to check if the same sign has all dets: det[a,b,p,o], det[b,c,p,o], det[c,a,p,o] and the det[a,b,c,p] has an opposite sign
char check_tradeability(unsigned char *orientation,unsigned int th1,unsigned int th2,t_delaunay *d)
{
double *p[4];
// Test abpo
if (bjsort4(&d->lvec[d->th[th1][0]],&d->lvec[d->th[th1][1]],&d->lvec[d->th[th1][3]],&d->lvec[d->th[th2][3]],(void**)&p[0],(void**)&p[1],(void**)&p[2],(void**)&p[3])%2)
  orientation[0]=(char)(!sign_ldet4(p));
else
  orientation[0]=(char)( sign_ldet4(p));
// Test bcpo
if (bjsort4(&d->lvec[d->th[th1][1]],&d->lvec[d->th[th1][2]],&d->lvec[d->th[th1][3]],&d->lvec[d->th[th2][3]],(void**)&p[0],(void**)&p[1],(void**)&p[2],(void**)&p[3])%2)
  orientation[1]=(char)(!sign_ldet4(p));
else
  orientation[1]=(char)( sign_ldet4(p));
// Test capo
if (bjsort4(&d->lvec[d->th[th1][2]],&d->lvec[d->th[th1][0]],&d->lvec[d->th[th1][3]],&d->lvec[d->th[th2][3]],(void**)&p[0],(void**)&p[1],(void**)&p[2],(void**)&p[3])%2)
  orientation[2]=(char)(!sign_ldet4(p));
else
  orientation[2]=(char)( sign_ldet4(p));
// Test pabc
if (bjsort4(&d->lvec[d->th[th1][0]],&d->lvec[d->th[th1][1]],&d->lvec[d->th[th1][2]],&d->lvec[d->th[th1][3]],(void**)&p[0],(void**)&p[1],(void**)&p[2],(void**)&p[3])%2)
  orientation[3]=(char)(!sign_ldet4(p));
else
  orientation[3]=(char)( sign_ldet4(p));
//Check tradeability
return ( ( !orientation[3] && orientation[0] && orientation[1] && orientation[2] )||( orientation[3] && !orientation[0] && !orientation[1] && !orientation[2] ) );
}

//This function check the tradeability of edges in reflexes. Note the tetranedrons shopuld be aligned already!
char check_edges_tradeability(unsigned int *tradable_th,unsigned char *orientation,unsigned int th1,unsigned int th2,t_delaunay *d)
{
unsigned char _i,_j;
//Check ab edge
if (orientation[0x0]==orientation[0x3]) // abpo == pabc
  {// The vertice 'a' is a reflex
  for (tradable_th[0x0]=0;tradable_th[0x0]<d->size_th;tradable_th[0x0]++) //abpo exists?
    if ( (d->nth[tradable_th[0x0]][0x0]!=0xFFFFFFFF)||(d->nth[tradable_th[0x0]][0x1]!=0xFFFFFFFF)||(d->nth[tradable_th[0x0]][0x2]!=0xFFFFFFFF)||(d->nth[tradable_th[0x0]][0x3]!=0xFFFFFFFF) )
      if ( (find_in_row(d->th[th1][0x0],0x4,d->th[tradable_th[0x0]])!=(unsigned int)-1)&&(find_in_row(d->th[th1][0x1],0x4,d->th[tradable_th[0x0]])!=(unsigned int)-1)&&
           (find_in_row(d->th[th1][0x3],0x4,d->th[tradable_th[0x0]])!=(unsigned int)-1)&&(find_in_row(d->th[th2][0x3],0x4,d->th[tradable_th[0x0]])!=(unsigned int)-1) ) //the reflex orientation and abpo exists
        goto A;
  return FALSE;
  A: orientation[0x0]=TRUE;
  }
else
  orientation[0x0]=FALSE;
//Check bc edge
if (orientation[0x1]==orientation[0x3]) // bcpo == pabc
  {// The vertice 'b' is a reflex
  for (tradable_th[0x1]=0;tradable_th[0x1]<d->size_th;tradable_th[0x1]++) // bcpo exists?
    if ( (d->nth[tradable_th[0x1]][0x0]!=0xFFFFFFFF)||(d->nth[tradable_th[0x1]][0x1]!=0xFFFFFFFF)||(d->nth[tradable_th[0x1]][0x2]!=0xFFFFFFFF)||(d->nth[tradable_th[0x1]][0x3]!=0xFFFFFFFF) )
      if ( (find_in_row(d->th[th1][0x1],0x4,d->th[tradable_th[0x1]])!=(unsigned int)-1)&&(find_in_row(d->th[th1][0x2],0x4,d->th[tradable_th[0x1]])!=(unsigned int)-1)&&
           (find_in_row(d->th[th1][0x3],0x4,d->th[tradable_th[0x1]])!=(unsigned int)-1)&&(find_in_row(d->th[th2][0x3],0x4,d->th[tradable_th[0x1]])!=(unsigned int)-1) ) //the reflex orientation and abpo exists
        goto B;
  return FALSE;
  B: orientation[0x1]=TRUE;
  }
else
  orientation[0x1]=FALSE;
//Check ca edge
if (orientation[0x2]==orientation[0x3]) // acpo == pabc
  {
  for (tradable_th[0x2]=0;tradable_th[0x2]<d->size_th;tradable_th[0x2]++) // acpo exists?
    if ( (d->nth[tradable_th[0x2]][0x0]!=0xFFFFFFFF)||(d->nth[tradable_th[0x2]][0x1]!=0xFFFFFFFF)||(d->nth[tradable_th[0x2]][0x2]!=0xFFFFFFFF)||(d->nth[tradable_th[0x2]][0x3]!=0xFFFFFFFF) )
      if ( (find_in_row(d->th[th1][0x0],0x4,d->th[tradable_th[0x2]])!=(unsigned int)-1)&&(find_in_row(d->th[th1][0x2],0x4,d->th[tradable_th[0x2]])!=(unsigned int)-1)&&
           (find_in_row(d->th[th1][0x3],0x4,d->th[tradable_th[0x2]])!=(unsigned int)-1)&&(find_in_row(d->th[th2][0x3],0x4,d->th[tradable_th[0x2]])!=(unsigned int)-1) ) //the reflex orientation and abpo exists
        goto C;
  return FALSE;
  C: orientation[0x2]=TRUE;
  }
else
  orientation[0x2]=FALSE;

//Order items
for (_i=0,_j=0;_i<0x3;_i++)
  if (orientation[_i])
    orientation[_j++]=_i;
return _j;
}


//this function print the tetraherons in triangulation
inline void show_delaunay(t_delaunay *d)
{
unsigned int _th,_i,_j;

for (_th=0;_th<d->size_th;_th++)
  {
  printf("\t%d:\t%d\t%d\t%d\t%d\n",_th,d->th[_th][0x0],d->th[_th][0x1],d->th[_th][0x2],d->th[_th][0x3]);
  for (_i=0;_i<0x4;_i++)
    for (_j=_i+1;_j<0x4;_j++)
      if (d->th[_th][_i]==d->th[_th][_j])
        printf("ERROR! Tetrahedron %d corrupted!!!!!\n",_th);
  }
}


//This function complete Delaunay triangulation. It require the regular delaunay of certain fragment and number of points to be triangulated.
char triangulate_3d_delaunay(unsigned int nsize,t_delaunay *d)
{
unsigned int _t[0x4],fsize,f_id;
unsigned char   _tc[0x4];
t_fstack   *fstack;

//Prepare memory
fstack=(t_fstack*)malloc(sizeof(t_fstack)*0xFF);

//Perform triangulation cycle
while (d->size!=nsize)
  //Find tetrahedron that keeps new point inside and check if there any space left for this point inside tetrahedron
//  if ((find_tetrahedron(d->size,&_t[0x0],d)))
  if ((ffind_tetrahedron(&d->lvec[d->size],&_t[0x0],d)))
    {
    //Init stack
    fsize=0;
    f_id=0;
    //Perform one four trade
    if (!(trade_one_four(d->size,_t[0x0],&fsize,&fstack,d))) goto ERROR;
    //Perform series of trades
    do{
      if ((fstack[f_id].th[0x0]!=0xFFFFFFFF)&&(fstack[f_id].th[0x1]!=0xFFFFFFFF)&&(d->th[fstack[f_id].th[0x0]][0x0]!=0xFFFFFFFF)&&(d->th[fstack[f_id].th[0x1]][0x0]!=0xFFFFFFFF)&&((get_tetrahedrons(&fstack[f_id],d))) ) //Note. This function aligns them as well.
        {
        //Chcek facet local nonregularity and its tradebility
        if (!(check_regularity(fstack[f_id].th[0x0],fstack[f_id].th[0x1],d)))
          {
          if ( (check_tradeability(_tc,fstack[f_id].th[0x0],fstack[f_id].th[0x1],d)) ) //Is the system is convex. Note, we use id as massive of four chars.
            {//Perform two to three trade
            if (!(trade_two_three(fstack[f_id].th[0x0],fstack[f_id].th[0x1],&fsize,&fstack,d))) goto ERROR;
            }
          else //Get reflexes flipaility
            switch (check_edges_tradeability(_t,_tc,fstack[f_id].th[0x0],fstack[f_id].th[0x1],d)) //Note, we use id as a massive of four chars
              {
              case 0x0 : {//Nothing to do with this case.
                         break;
                         }
              case 0x1 : { //The redundant link detected. Eliminate it.
                         switch (_tc[0x0])
                           {
                           case 0x0 : {                                                                                                            break;} //ab-link
                           case 0x1 : { _t[0x0]=_t[0x1];
                                        memexchange(&d->th [fstack[f_id].th[0x0]][0x0], &d->th[fstack[f_id].th[0x0]][0x2],sizeof(unsigned int));
                                        memexchange(&d->th [fstack[f_id].th[0x1]][0x0], &d->th[fstack[f_id].th[0x1]][0x2],sizeof(unsigned int));
                                        memexchange(&d->nth[fstack[f_id].th[0x0]][0x0],&d->nth[fstack[f_id].th[0x0]][0x2],sizeof(unsigned int));
                                        memexchange(&d->nth[fstack[f_id].th[0x1]][0x0],&d->nth[fstack[f_id].th[0x1]][0x2],sizeof(unsigned int)); break;} //bc-link
                           case 0x2 : { _t[0x0]=_t[0x2];
                                        memexchange( &d->th[fstack[f_id].th[0x0]][0x1], &d->th[fstack[f_id].th[0x0]][0x2],sizeof(unsigned int));
                                        memexchange( &d->th[fstack[f_id].th[0x1]][0x1], &d->th[fstack[f_id].th[0x1]][0x2],sizeof(unsigned int));
                                        memexchange(&d->nth[fstack[f_id].th[0x0]][0x1],&d->nth[fstack[f_id].th[0x0]][0x2],sizeof(unsigned int));
                                        memexchange(&d->nth[fstack[f_id].th[0x1]][0x1],&d->nth[fstack[f_id].th[0x1]][0x2],sizeof(unsigned int)); break;} //ca-link
                           default  : {printf("ERROR. 3D Delaunay has detecter strange vertice to 3->2 flip!!!!\n"); goto ERROR; }
                           }
                         if (!(trade_three_two(fstack[f_id].th[0x0],fstack[f_id].th[0x1],_t[0x0],&fsize,&fstack,d))) goto ERROR;
                         break;
                         }
              case 0x2 : {//The redundant point at a common link is detected. Eliminate thetrahedron.
                         switch(_tc[0x0]+_tc[0x1])
                           {
                           case 0x1 : { _t[0x0]=_t[0x1];
                                        memexchange( &d->th[fstack[f_id].th[0x0]][0x0], &d->th[fstack[f_id].th[0x0]][0x1],sizeof(unsigned int));
                                        memexchange( &d->th[fstack[f_id].th[0x1]][0x0], &d->th[fstack[f_id].th[0x1]][0x1],sizeof(unsigned int));
                                        memexchange(&d->nth[fstack[f_id].th[0x0]][0x0],&d->nth[fstack[f_id].th[0x0]][0x1],sizeof(unsigned int));
                                        memexchange(&d->nth[fstack[f_id].th[0x1]][0x0],&d->nth[fstack[f_id].th[0x1]][0x1],sizeof(unsigned int)); break; } // b point
                           case 0x2 : { _t[0x1]=_t[0x2];                                                                                         break; } // a point
                           case 0x3 : { _t[0x0]=_t[0x1];
                                        _t[0x1]=_t[0x2];
                                        memexchange( &d->th[fstack[f_id].th[0x0]][0x0], &d->th[fstack[f_id].th[0x0]][0x2],sizeof(unsigned int));
                                        memexchange( &d->th[fstack[f_id].th[0x1]][0x0], &d->th[fstack[f_id].th[0x1]][0x2],sizeof(unsigned int));
                                        memexchange(&d->nth[fstack[f_id].th[0x0]][0x0],&d->nth[fstack[f_id].th[0x0]][0x2],sizeof(unsigned int));
                                        memexchange(&d->nth[fstack[f_id].th[0x1]][0x0],&d->nth[fstack[f_id].th[0x1]][0x2],sizeof(unsigned int)); break; } // c point
                           default  : {printf("ERROR. 3D Delaunay has detecter strange vertice to 4->1 flip!!!!\n"); goto ERROR; }
                           }
                         if (!(trade_four_one(fstack[f_id].th[0x0],fstack[f_id].th[0x1],_t[0x0],_t[0x1],&fsize,&fstack,d))) goto ERROR;
                         break;
                         }
              default  : {//This case should never occur. Raport about error.
                         printf("ERROR. 3D delaunay has detected three edges flip!\n"); goto ERROR;
                         }
              }
          }
        }
      }while(++f_id!=fsize);
    //Rebuild triangulation (discard redundant tetraiders)
    _t[0x0]=d->size_th;
    while (_t[0x0]--)
      if (d->th[_t[0x0]][0x0]==0xFFFFFFFF) //Redundance marker
        {
        //Copy tetrahedrons
        memcpy(d->th[_t[0x0]],d->th[--d->size_th],sizeof(unsigned int)*0x4);
        memcpy(d->nth[_t[0x0]],d->nth[d->size_th],sizeof(unsigned int)*0x4);
        //Update tetrahedrons neighbour list (four entries)
        for (_t[0x1]=0;_t[0x1]<0x4;_t[0x1]++)
          if (d->nth[_t[0x0]][_t[0x1]]!=0xFFFFFFFF)
            for (_t[0x2]=0;_t[0x2]<0x4;_t[0x2]++)
              if (d->nth[d->nth[_t[0x0]][_t[0x1]]][_t[0x2]]==d->size_th)
                {
                d->nth[d->nth[_t[0x0]][_t[0x1]]][_t[0x2]]=_t[0x0];
                break;
                }
        }

    //Check if all neighbours are localy regular
    for (_t[0x0]=0;_t[0x0]<d->size_th;_t[0x0]++)
      for (_t[0x1]=0;_t[0x1]<0x4;_t[0x1]++)
        if ( (d->nth[_t[0x0]][_t[0x1]]>_t[0x0])&&(d->nth[_t[0x0]][_t[0x1]]!=0xFFFFFFFF) )
          for (_t[0x2]=0;_t[0x2]<0x4;_t[0x2]++)
            if (d->nth[d->nth[_t[0x0]][_t[0x1]]][_t[0x2]]==_t[0x0])
              {
              memexchange(&d->th[_t[0x0]][_t[0x1]],&d->th[_t[0x0]][0x3],sizeof(unsigned int));
              memexchange(&d->nth[_t[0x0]][_t[0x1]],&d->nth[_t[0x0]][0x3],sizeof(unsigned int));
              memexchange(&d->th[d->nth[_t[0x0]][0x3]][_t[0x2]],&d->th[d->nth[_t[0x0]][0x3]][0x3],sizeof(unsigned int));
              memexchange(&d->nth[d->nth[_t[0x0]][0x3]][_t[0x2]],&d->nth[d->nth[_t[0x0]][0x3]][0x3],sizeof(unsigned int));
              if (!(check_regularity(_t[0x0],d->nth[_t[0x0]][0x3],d)))
                { printf("ERROR. Nonregularity detected in Delaunay triangulation! NP=%d\n",d->size); goto ERROR;}
              memexchange(&d->nth[d->nth[_t[0x0]][0x3]][_t[0x2]],&d->nth[d->nth[_t[0x0]][0x3]][0x3],sizeof(unsigned int));
              memexchange(&d->th[d->nth[_t[0x0]][0x3]][_t[0x2]],&d->th[d->nth[_t[0x0]][0x3]][0x3],sizeof(unsigned int));
              memexchange(&d->nth[_t[0x0]][_t[0x1]],&d->nth[_t[0x0]][0x3],sizeof(unsigned int));
              memexchange(&d->th[_t[0x0]][_t[0x1]],&d->th[_t[0x0]][0x3],sizeof(unsigned int));
              break;
              }

    //Add new point in triangulation   show_delaunay(d);
    d->size++;
    if (!(d->size%0xFF)) printf("%d atoms successfuly added to triangulation\n",d->size);
    }
  else d->size++;
//Free some memory
free(fstack);
return TRUE;
ERROR: free(fstack);
return FALSE;
}

//This function construct Delaunay triangulation structure
t_delaunay *alloc_delaunay(unsigned int size,unsigned int size_th)
{
t_delaunay *d;
if ( ((d=(t_delaunay*)calloc(sizeof(t_delaunay),0x1)))                     &&
     ((d->lvec=(t_lvec*)malloc(sizeof(t_lvec)*(size+4))))                  &&
     ((d->th =(unsigned int (*)[4])malloc(sizeof(unsigned int)*size_th*4)))&&
     ((d->nth=(unsigned int (*)[4])malloc(sizeof(unsigned int)*size_th*4))) )
  {
  d->_size=d->size;
  d->_size_th=size_th;
  return d;
  }
free_delaunay(d);
return FALSE;
}


//This function clone delaunay triangulations
t_delaunay *clone_delaunay(unsigned int size,unsigned int size_th,t_delaunay *d_o)
{
t_delaunay *d;
if (!(d=alloc_delaunay(size,size_th))) return FALSE;
memcpy(d->lvec,d_o->lvec,sizeof(t_lvec)*d_o->size);
d->size=d_o->size;
memcpy(d->th,d_o->th,sizeof(unsigned int)*4*d_o->size_th);
memcpy(d->nth,d_o->nth,sizeof(unsigned int)*4*d_o->size_th);
d->size_th=d_o->size_th;
return d;
}


//This ia a destuctor for Delaunay triangulation
void free_delaunay(t_delaunay *d)
{
if (d)
  {
  if (d->lvec) free(d->lvec);
  if (d->th)   free(d->th);
  if (d->nth)  free(d->nth);
  free(d);
  }
}

//This function builds 3d triangulation on th set of given points
//NOTE. The first 4 points has to be reserved for superstructure
char triangulate_3d_set(t_delaunay *d)
{
unsigned int size;
//Fill the superstructure that is a reclinear tetrahedron
inscribe_tetrahedron(1.e3,&d->lvec[0],&d->lvec[1],&d->lvec[2],&d->lvec[3],d->size-4,&d->lvec[4]);
d->lvec[0].w=calc_vec_norm((t_vec*)&d->lvec[0]), d->lvec[1].w=calc_vec_norm((t_vec*)&d->lvec[1]), d->lvec[2].w=calc_vec_norm((t_vec*)&d->lvec[2]), d->lvec[3].w=calc_vec_norm((t_vec*)&d->lvec[3]);
//Mount initial tetrahedron
d->th[0][0]=0, d->th[0][1]=1, d->th[0][2]=2, d->th[0][3]=3;
memset(d->nth[0],0xFF,sizeof(unsigned int)*0x4);
d->size_th=1;
size=d->size;
d->size=4;
//Finish the triangulation
if ((triangulate_3d_delaunay(size,d))) return TRUE;
else                                                             return FALSE;
}


//This function performs delaunay triangulation via mol interface
t_delaunay *triangulate_3d_mol(t_list *it,t_vec *r,t_mol *mol,t_top *top)
{
unsigned int _i,_j;
t_delaunay *d;
//Prepare memory

//Copy and lift coords into lvec
if (it)
  {
  _j=0;
  _i=it->size; while (_i--) _j+=mol->anchors->list[it->list[_i]].size;
  if (!(d=alloc_delaunay(_j+0x4,0x8))) return FALSE;
  d->size=4;
  _i=it->size;
  while (_i--)
    {
    _j=mol->anchors->list[it->list[_i]].size;
    while (_j--)
      {
      d->lvec[d->size].i=r[mol->anchors->list[it->list[_i]].list[_j]].i;
      d->lvec[d->size].j=r[mol->anchors->list[it->list[_i]].list[_j]].j;
      d->lvec[d->size].k=r[mol->anchors->list[it->list[_i]].list[_j]].k;
      d->lvec[d->size].w=calc_vec_norm((t_vec*)&d->lvec[d->size])-sqrd(top->ff_a[mol->ytypes[mol->anchors->list[it->list[_i]].list[_j]]].rvdw);
      d->size++;
      }
    }
  _i=d->size;
  }
else
  {
  if (!(d=alloc_delaunay(mol->natoms+0x4,0x8))) return FALSE;
  for (d->size=0;d->size<mol->natoms;d->size++)
    {
    d->lvec[d->size+4].i=r[d->size].i;
    d->lvec[d->size+4].j=r[d->size].j;
    d->lvec[d->size+4].k=r[d->size].k;
    d->lvec[d->size+4].w=calc_vec_norm(&r[d->size])-sqrd(top->ff_a[mol->ytypes[d->size]].rvdw);
    }
  }
//Triangulate selection
d->size=4;
if ( (triangulate_3d_set(d))) return d;
else  { free_delaunay(d); return FALSE; }
}



//This function performs delaunay triangulation via mol interface
char add_mol_to_triangulation(t_list *it,t_vec *r,t_mol *mol,t_top *top,t_delaunay *d)
{
register unsigned int _i,_j,_shift;
void *vp;

if (it)
  {
  _shift=0; _i=it->size; while(_i--) _shift+=mol->anchors->list[it->list[_i]].size;
  if (!(vp=realloc(d->lvec,sizeof(t_lvec)*(_shift+d->size+0x4)))) return FALSE;
  else d->lvec=(t_lvec*)vp;
  //Setup new points lvecs
  _i=it->size;
  while (_i--)
    {
    _j=mol->anchors->list[it->list[_i]].size;
    while (_j--)
      {
      d->lvec[d->size].i=r[mol->anchors->list[it->list[_i]].list[_j]].i;
      d->lvec[d->size].j=r[mol->anchors->list[it->list[_i]].list[_j]].j;
      d->lvec[d->size].k=r[mol->anchors->list[it->list[_i]].list[_j]].k;
      d->lvec[d->size].w=calc_vec_norm((t_vec*)&d->lvec[d->size])-sqrd(top->ff_a[mol->ytypes[mol->anchors->list[it->list[_i]].list[_j]]].rvdw);
      d->size++;
      }
    }
  d->size-=_shift;
  }
else
  {
  _shift=mol->natoms;
  if (!(vp=realloc(d->lvec,sizeof(t_lvec)*(_shift+d->size)))) return FALSE;
  else d->lvec=(t_lvec*)vp;
  //Setup new points lvecs
  _i=mol->natoms;
  while (_i--)
    {
    d->lvec[d->size+_i].i=r[_i].i, d->lvec[d->size+_i].j=r[_i].j, d->lvec[d->size+_i].k=r[_i].k;
    d->lvec[d->size+_i].w=calc_vec_norm(&r[_i])-sqrd(top->ff_a[mol->ytypes[_i]].rvdw);
    }
  }
//update triangulation
return ((triangulate_3d_delaunay(d->size+_shift,d)));
}

//---------------------------------- T E T R A H E D R O N S   P A R T --------------------------------

//This function calculates tetrahedron volume
//Note. this function returns 'signed' volume
double calc_th_volume(t_vec *a,t_vec *b,t_vec *c,t_vec *d) 
{
t_tensor T;
T[0][0]=a->i-b->i, T[0][1]=a->j-b->j, T[0][2]=a->k-b->k;
T[1][0]=a->i-c->i, T[1][1]=a->j-c->j, T[1][2]=a->k-c->k;
T[2][0]=a->i-d->i, T[2][1]=a->j-d->j, T[2][2]=a->k-d->k;
return TENSOR_DET(T)/6.00;
}
double calc_thetrahedron_volume(t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d)
{
return calc_th_volume((t_vec*) a,(t_vec*) b,(t_vec*) c,(t_vec*) d);
}

//This function calculates the center of mass of solid gomogenic tetrahedron as 1/3 of the line that connect a mediane intersection of triangle to opposite vertice
void calc_th_cm(t_vec *cm,t_vec *a,t_vec *b,t_vec *c,t_vec *d)
{
cm->i=(a->i+b->i+c->i+d->i)/4., cm->j=(a->j+b->j+c->j+d->j)/4., cm->k=(a->k+b->k+c->k+d->k)/4.; //Funny, right?
}

//This function calculates intrinsic radii of sphere between four other spheres and returns signed volume of tetrahedron
//E.N.Sickafus and Neil A. Mackie "Interstitial space in hard-shere clusters" Acta Cryst. A30, 850-851 1974
//Note if tetrahedrons volume less than cutoff calculations are not done and 0.0000 is returned
double calc_thetrahedron_insphere(double cutoff,t_lvec *insphere,t_lvec *a,t_lvec *b,t_lvec *c,t_lvec *d)
{
double Ra,Rb,Rc,Rd,Det,E,F,G,H,P,Q,U,V,T;
t_tensor t;

//Calculate volume
t[0][0]=a->i-b->i, t[0][1]=a->i-c->i, t[0][2]=a->i-d->i;
t[1][0]=a->j-b->j, t[1][1]=a->j-c->j, t[1][2]=a->j-d->j;
t[2][0]=a->k-b->k, t[2][1]=a->k-c->k, t[2][2]=a->k-d->k;
Det=TENSOR_DET(t);
if (fabs(Det)<cutoff*6.00) return 0.000;
Ra=sqrt(a->i*a->i+a->j*a->j+a->k*a->k-a->w);
Rb=sqrt(b->i*b->i+b->j*b->j+b->k*b->k-b->w);
Rc=sqrt(c->i*c->i+c->j*c->j+c->k*c->k-c->w);
Rd=sqrt(d->i*d->i+d->j*d->j+d->k*d->k-d->w);
//Calculate sub-determinants
t[0][0]=a->w-b->w, t[0][1]=a->w-c->w, t[0][2]=a->w-d->w;
E=0.5*TENSOR_DET(t);
t[0][0]=Rb-Ra,     t[0][1]=Rc-Ra,     t[0][2]=Rd-Ra;
F=TENSOR_DET(t);
t[0][0]=a->i-b->i, t[0][1]=a->i-c->i, t[0][2]=a->i-d->i;

t[1][0]=a->w-b->w, t[1][1]=a->w-c->w, t[1][2]=a->w-d->w;
G=0.5*TENSOR_DET(t);
t[1][0]=Rb-Ra,     t[1][1]=Rc-Ra,     t[1][2]=Rd-Ra;
H=TENSOR_DET(t);
t[1][0]=a->j-b->j, t[1][1]=a->j-c->j, t[1][2]=a->j-d->j;

t[2][0]=a->w-b->w, t[2][1]=a->w-c->w, t[2][2]=a->w-d->w;
P=0.5*TENSOR_DET(t);
t[2][0]=Rb-Ra,     t[2][1]=Rc-Ra,     t[2][2]=Rd-Ra;
Q=TENSOR_DET(t);

T=sqrd(F/Det)+sqrd(H/Det)+sqrd(Q/Det)-1.00;
U=(a->i-E/Det)*F/Det+(a->j-G/Det)*H/Det+(a->k-P/Det)*Q/Det+Ra;
V=sqrd(a->i-E/Det)+sqrd(a->j-G/Det)+sqrd(a->k-P/Det)-Ra*Ra;

if (U/T<0.00) insphere->w=U/T*(1.00-sqrt(1.00-V*T/U/U));
else          insphere->w=U/T*(1.00+sqrt(1.00-V*T/U/U));
insphere->i=E/Det+F/Det*insphere->w;
insphere->j=G/Det+H/Det*insphere->w;
insphere->k=P/Det+Q/Det*insphere->w;

Ra=sqrd(insphere->i-a->i)+sqrd(insphere->j-a->j)+sqrd(insphere->k-a->k)-sqrd(Ra+insphere->w);
Rb=sqrd(insphere->i-b->i)+sqrd(insphere->j-b->j)+sqrd(insphere->k-b->k)-sqrd(Rb+insphere->w);
Rc=sqrd(insphere->i-c->i)+sqrd(insphere->j-c->j)+sqrd(insphere->k-c->k)-sqrd(Rc+insphere->w);
Rd=sqrd(insphere->i-d->i)+sqrd(insphere->j-d->j)+sqrd(insphere->k-d->k)-sqrd(Rd+insphere->w);
return Det/6.00;
}

//This empirical function for evaluation of tetrahedron inscribed sphere
double embed_sphere_in_tetrahedron(t_vec *o,t_vec *a,t_vec *b,t_vec *c,t_vec *d)
{
t_vec ba, ca, da, x, w, xw;
double ax, bx, ab, ab2, ri;
double _d;

//Calculate X coords
ba.i=b->i-a->i, ba.j=b->j-a->j, ba.k=b->k-a->k;                            //  3 flops
ca.i=c->i-a->i, ca.j=c->j-a->j, ca.k=c->k-a->k;                            //  3 flops
da.i=d->i-a->i, da.j=d->j-a->j, da.k=d->k-a->k;                            //  3 flops
ab2=calc_vec_norm(&ba);                                                    //  5 flops 
_d=calc_vec_vec_scalar_product(&ca,&ba)/ab2;                               //  6 flops
ax=sqrd(ca.i-ba.i*_d)+sqrd(ca.j-ba.j*_d)+sqrd(ca.k-ba.k*_d);               // 12 flops  
_d=calc_vec_vec_scalar_product(&da,&ba)/ab2;                               //   6 flops 
bx=sqrd(da.i-ba.i*_d)+sqrd(da.j-ba.j*_d)+sqrd(da.k-ba.k*_d);               // 12 flops  
_d=sqrt(bx/ax)+1.;                                                         //   2 flops + 1 root 
x.i=c->i+(d->i-c->i)/_d, x.j=c->j+(d->j-c->j)/_d, x.k=c->k+(d->k-c->k)/_d; //   9 flops
//Calculate W coords
ax=sqrt(sqrd(a->i-x.i)+sqrd(a->j-x.j)+sqrd(a->k-x.k));                     //  8 flops + 1 root
bx=sqrt(sqrd(b->i-x.i)+sqrd(b->j-x.j)+sqrd(b->k-x.k));                     //  8 flops + 1 root
_d=bx/ax+1.;                                                               //  2 flops
w.i=a->i+ba.i/_d, w.j=a->j+ba.j/_d, w.k=a->k+ba.k/_d;                      //  6 flops  
//Calculate ri
xw.i=x.i-w.i, xw.j=x.j-w.j, xw.k=x.k-w.k;                                  //  3 flops 
_d=calc_vec_norm(&xw);                                                     //  5 flops 
ab=sqrt(ab2);                                                              //                  1 root
ri=(ax+bx+ab)/2.;                                                          //   3 flops 
ri=sqrt((ri-ax)*(ri-bx)*(ri-ab)/ri);                                       //   6 flops + 1 root
//Calculate insphere center coordinates
_d*=(1.-sqrd(calc_vec_vec_scalar_product(&ba,&xw))/ab2/_d);                //  10 flops
multiple_vec_scalar(&xw,&xw,ri/sqrt(_d));                                  //   4 flops + 1 root
o->i=w.i+xw.i, o->j=w.j+xw.j, o->k=w.k+xw.k;                               //   3 flops
return ri;                                                                 // 119 flops + 6 roots  total  
}





