#include "y_pmatrix.h"

/*
The logic:
0. a) split in nthread pieces if possible
   b) if it too big split in bigger amount of pieces each of MAX_... size
   c) if it too small split it in smaller amount of pieces each of MIN_... size
1. a) if tail+1 is <2.MIN_... calculate <MIN_... piece there and do n of lenght
   b) otherwise split into two almost equal pieces and push them into the end of queue
*/

//This is worker function 
unsigned int _pcalc_vect_norm(char *pcard,void *vp)
{
t_ypsystem *ypsys=(t_ypsystem*)vp;
register unsigned int pjob_id, size;
register double *vect, _d;
pjob_id=((unsigned int*)pcard)[0];                  //Set job id
size=((unsigned int*)pcard)[1];                     //Set job lenght
vect=((double**)(pcard+sizeof(unsigned int)*2))[0]; //Set vectors fragment address   
_d=*vect**vect; while (--size) { vect++, _d+=*vect**vect; }         //Piece of parallel job to do
if ( (pthread_mutex_lock  (&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
*(((double**)(pcard+sizeof(unsigned int)*2))[1])+=_d; //Summ vector lenght address
size=++ypsys->pjobs[pjob_id].m_pcards;                          //Standard job incremental part
if ( (pthread_mutex_unlock(&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
return size;
}
//This function calculates vectors square norm in parallel
inline char pcalc_vect_norm(double *norm,unsigned int size,double *vect,unsigned int pjob_id,t_ypsystem *ypsys)
{
register unsigned int _l, len, num;
void* hpcard=ypsys->hpcards+ypsys->sizeof_hpcard*ypsys->_hpstore;

*norm=0.;
     if (size>=MAX_PQANT_pcalc_vect_norm*ypsys->n_threads) 
       {
       len=MAX_PQANT_pcalc_vect_norm, num=size/len;
       //Pass the last two fragments
       SUBMIT_JOB:
       if (!(size%len))
         {
         //Do main body
         MAIN_BODY_LOOP:
         if (!(push_pjob(pjob_id,num,_pcalc_vect_norm,ypsys))) return FALSE; 
         _l=num; 
         while (_l--) 
           {
           ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
           ((unsigned int*)hpcard)[1]=len;      //Set job lenght
           hpcard+=sizeof(unsigned int)*2;
           ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
           ((double**)hpcard)[1]=norm;          //Set vector lenght address  
           if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
           else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
           }
         }
       else  
         {
         _l=size-num*len;
         if (_l<2*MIN_PQANT_pcalc_vect_norm)
           {
           for (_l=len*num;_l!=size;_l++) (*norm)+=vect[_l]*vect[_l]; //Do the smallest fragment
           goto MAIN_BODY_LOOP; 
           }
         else
           {
           //Do main body first to store short cards in the end of stack
           if (!(push_pjob(pjob_id,--num+2,_pcalc_vect_norm,ypsys))) return FALSE; 
           _l=num; 
           while (_l--) 
             {
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=len;      //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             }
           _l=size-num*len;
           if (!(_l%2))
             {//Store two the same cards
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num+_l]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             }
           else
             {
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l+1;     //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num+_l]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             }
           } 
         }
       }
else if (size>=MIN_PQANT_pcalc_vect_norm*ypsys->n_threads) { len=size/ypsys->n_threads, num=ypsys->n_threads;  goto SUBMIT_JOB; }
else if (size>=MIN_PQANT_pcalc_vect_norm) { len=MIN_PQANT_pcalc_vect_norm, num=size/MIN_PQANT_pcalc_vect_norm; goto SUBMIT_JOB; }
else { while (size--) { (*norm)+=*vect**vect, vect++; } return NTNF; } //Do the small fragment serially
return TRUE;
}

//This is worker function 
unsigned int _pcalc_vect_scaled_norm(char *pcard,void *vp)
{
t_ypsystem *ypsys=(t_ypsystem*)vp;
register unsigned int pjob_id, size;
register double *vect, _d, scale;
pjob_id=((unsigned int*)pcard)[0];                  //Set job id
size=((unsigned int*)pcard)[1];                     //Set job lenght
vect=((double**)(pcard+sizeof(unsigned int)*2))[0]; //Set vectors fragment address
scale=*(double*)(pcard+(sizeof(unsigned int)+sizeof(double*))*2);   
*vect*=scale, _d=*vect**vect; while (--size) { vect++, *vect*=scale, _d+=*vect**vect; }     //Piece of parallel job to do
if ( (pthread_mutex_lock  (&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
*(((double**)(pcard+sizeof(unsigned int)*2))[1])+=_d; //Summ vector lenght address
size=++ypsys->pjobs[pjob_id].m_pcards;                          //Standard job incremental part
if ( (pthread_mutex_unlock(&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
return size;
}
//This function calculates scaled vector lenght in parallel
inline char pcalc_vect_scaled_norm(double *norm,double scale,unsigned int size,double *vect,unsigned int pjob_id,t_ypsystem *ypsys)
{
register unsigned int _l, len, num;
void* hpcard=ypsys->hpcards+ypsys->sizeof_hpcard*ypsys->_hpstore;

*norm=0.;
     if (size>=MAX_PQANT_pcalc_vect_norm*ypsys->n_threads) 
       {
       len=MAX_PQANT_pcalc_vect_norm, num=size/len;
       //Pass the last two fragments
       SUBMIT_JOB:
       if (!(size%len))
         {
         //Do main body
         MAIN_BODY_LOOP:
         if (!(push_pjob(pjob_id,num,_pcalc_vect_scaled_norm,ypsys))) return FALSE; 
         _l=num; 
         while (_l--) 
           {
           ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
           ((unsigned int*)hpcard)[1]=len;      //Set job lenght
           hpcard+=sizeof(unsigned int)*2;
           ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
           ((double**)hpcard)[1]=norm;          //Set vector lenght address  
           hpcard+=sizeof(double*)*2;
           *(double*)hpcard=scale;              //Set vectors scale  
           if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
           else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-(sizeof(unsigned int)+sizeof(double*))*2; }
           }
         }
       else  
         {
         _l=size-num*len;
         if (_l<2*MIN_PQANT_pcalc_vect_norm)
           {
           for (_l=len*num;_l!=size;_l++) { vect[_l]*=scale, (*norm)+=vect[_l]*vect[_l]; } //Do the smallest fragment
           goto MAIN_BODY_LOOP; 
           }
         else
           {
           //Do main body first to store short cards in the end of stack
           if (!(push_pjob(pjob_id,--num+2,_pcalc_vect_scaled_norm,ypsys))) return FALSE; 
           _l=num; 
           while (_l--) 
             {
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=len;      //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             hpcard+=sizeof(double*)*2;
             *(double*)hpcard=scale;              //Set vectors scale  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-(sizeof(unsigned int)+sizeof(double*))*2; }
             }
           _l=size-num*len;
           if (!(_l%2))
             {//Store two the same cards
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             hpcard+=sizeof(double*)*2;
             *(double*)hpcard=scale;              //Set vectors scale  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-(sizeof(unsigned int)+sizeof(double*))*2; }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num+_l]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             hpcard+=sizeof(double*)*2;
             *(double*)hpcard=scale;              //Set vectors scale  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-(sizeof(unsigned int)+sizeof(double*))*2; }
             }
           else
             {
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             hpcard+=sizeof(double*)*2;
             *(double*)hpcard=scale;              //Set vectors scale  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-(sizeof(unsigned int)+sizeof(double*))*2; }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l+1;     //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num+_l]; //Set vectors fragment address   
             ((double**)hpcard)[1]=norm;          //Set vector lenght address  
             hpcard+=sizeof(double*)*2;
             *(double*)hpcard=scale;              //Set vectors scale  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-(sizeof(unsigned int)+sizeof(double*))*2; }
             }
           } 
         }
       }
else if (size>=MIN_PQANT_pcalc_vect_norm*ypsys->n_threads) { len=size/ypsys->n_threads, num=ypsys->n_threads;  goto SUBMIT_JOB; }
else if (size>=MIN_PQANT_pcalc_vect_norm) { len=MIN_PQANT_pcalc_vect_norm, num=size/MIN_PQANT_pcalc_vect_norm; goto SUBMIT_JOB; }
else { while (size--) { *vect*=scale, (*norm)+=*vect**vect, vect++; } return NTNF; } //Do the small fragment serially
return TRUE;
}
void get_stack_pcalc_vect_scaled_norm(unsigned int size,unsigned int n_threads,unsigned int *n_hpcards,size_t *hpcards_size,unsigned int *n_pcards,size_t *pcards_size)
{
*n_hpcards=(n_threads>size/MAX_PQANT_pcalc_vect_norm) ? n_threads : size/MAX_PQANT_pcalc_vect_norm;
*n_pcards=0;
*hpcards_size=sizeof(unsigned int)*2+sizeof(double*)*2+sizeof(double);
*pcards_size=0x0;
(*n_hpcards)++;
}

//This is worker function 
unsigned int _pcalc_vect_unsigned_summ(char *pcard,void *vp)
{
t_ypsystem *ypsys=(t_ypsystem*)vp;
register unsigned int pjob_id, size;
register double *vect, _d;
pjob_id=((unsigned int*)pcard)[0];                  //Set job id
size=((unsigned int*)pcard)[1];                     //Set job lenght
vect=((double**)(pcard+sizeof(unsigned int)*2))[0]; //Set vectors fragment address   
_d=fabs(*vect); while (--size) { vect++; _d+=fabs(*vect); }         //Piece of parallel job to do
if ( (pthread_mutex_lock  (&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
*(((double**)(pcard+sizeof(unsigned int)*2))[1])+=_d; //Summ vector lenght address
size=++ypsys->pjobs[pjob_id].m_pcards;                          //Standard job incremental part
if ( (pthread_mutex_unlock(&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
return size;
}
//This function calculates summ of vector's unsigned elements in parallel
inline char pcalc_vect_unsigned_summ(double *summ,unsigned int size,double *vect,unsigned int pjob_id,t_ypsystem *ypsys)
{
register unsigned int _l, len, num;
void* hpcard=ypsys->hpcards+ypsys->sizeof_hpcard*ypsys->_hpstore;

*summ=0.;
     if (size>=MAX_PQANT_pcalc_vect_norm*ypsys->n_threads) 
       {
       len=MAX_PQANT_pcalc_vect_norm, num=size/len;
       //Pass the last two fragments
       SUBMIT_JOB:
       if (!(size%len))
         {
         //Do main body
         MAIN_BODY_LOOP:
         if (!(push_pjob(pjob_id,num,_pcalc_vect_unsigned_summ,ypsys))) return FALSE; 
         _l=num; 
         while (_l--) 
           {
           ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
           ((unsigned int*)hpcard)[1]=len;      //Set job lenght
           hpcard+=sizeof(unsigned int)*2;
           ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
           ((double**)hpcard)[1]=summ;          //Set vector lenght address  
           if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
           else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
           }
         }
       else  
         {
         _l=size-num*len;
         if (_l<2*MIN_PQANT_pcalc_vect_norm)
           {
           for (_l=len*num;_l!=size;_l++) (*summ)+=fabs(vect[_l]); //Do the smallest fragment
           goto MAIN_BODY_LOOP; 
           }
         else
           {
           //Do main body first to store short cards in the end of stack
           if (!(push_pjob(pjob_id,--num+2,_pcalc_vect_unsigned_summ,ypsys))) return FALSE; 
           _l=num; 
           while (_l--) 
             {
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=len;      //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
             ((double**)hpcard)[1]=summ;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             }
           _l=size-num*len;
           if (!(_l%2))
             {//Store two the same cards
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num]; //Set vectors fragment address   
             ((double**)hpcard)[1]=summ;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num+_l]; //Set vectors fragment address   
             ((double**)hpcard)[1]=summ;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             }
           else
             {
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num]; //Set vectors fragment address   
             ((double**)hpcard)[1]=summ;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l+1;     //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*num+_l]; //Set vectors fragment address   
             ((double**)hpcard)[1]=summ;          //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             }
           } 
         }
       }
else if (size>=MIN_PQANT_pcalc_vect_norm*ypsys->n_threads) { len=size/ypsys->n_threads, num=ypsys->n_threads;  goto SUBMIT_JOB; }
else if (size>=MIN_PQANT_pcalc_vect_norm) { len=MIN_PQANT_pcalc_vect_norm, num=size/MIN_PQANT_pcalc_vect_norm; goto SUBMIT_JOB; }
else { while (size--) { (*summ)+=fabs(*vect), vect++; } return NTNF; } //Do the small fragment serially
return TRUE;
}
void get_stack_pcalc_vect_unsigned_summ(unsigned int size,unsigned int n_threads,unsigned int *n_hpcards,size_t *hpcards_size,unsigned int *n_pcards,size_t *pcards_size)
{
*n_hpcards=(n_threads>size/MAX_PQANT_pcalc_vect_norm) ? n_threads : size/MAX_PQANT_pcalc_vect_norm;
*n_pcards=0;
*hpcards_size=sizeof(unsigned int)*2+sizeof(double*)*2;
*pcards_size=0x0;
(*n_hpcards)++;
}

//This is worker function 
unsigned int _pcalc_vect_scalar_product(char *pcard,void *vp)
{
t_ypsystem *ypsys=(t_ypsystem*)vp;
register unsigned int pjob_id, size;
register double *vect, scalar;
pjob_id=((unsigned int*)pcard)[0];                  //Set job id
size=((unsigned int*)pcard)[1];                     //Set job lenght
vect=((double**)(pcard+sizeof(unsigned int)*2))[0]; //Set vectors fragment address   
scalar=*((double*)(pcard+sizeof(unsigned int)*2+sizeof(double)));
*vect*=scalar; while (--size) { vect++; *vect*=scalar; }        //Piece of parallel job to do
if ( (pthread_mutex_lock  (&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
size=++ypsys->pjobs[pjob_id].m_pcards;                          //Standard job incremental part
if ( (pthread_mutex_unlock(&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
return size;
}
//This function multiples vector on scalar in parallel
inline char pcalc_vect_scalar_product(double scalar,unsigned int size,double *vect,unsigned int pjob_id,t_ypsystem *ypsys)
{
register unsigned int _l, len, num;
void* hpcard=ypsys->hpcards+ypsys->sizeof_hpcard*ypsys->_hpstore;

     if (size>=MAX_PQANT_pcalc_vect_norm*ypsys->n_threads) 
       {
       len=MAX_PQANT_pcalc_vect_norm, num=size/len;
       //Pass the last two fragments
       SUBMIT_JOB:
       if (!(size%len))
         {
         //Do main body
         MAIN_BODY_LOOP:
         if (!(push_pjob(pjob_id,num,_pcalc_vect_scalar_product,ypsys))) return FALSE; 
         _l=num; 
         while (_l--) 
           {
           ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
           ((unsigned int*)hpcard)[1]=len;      //Set job lenght
           hpcard+=sizeof(unsigned int)*2;
           ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
           hpcard+=sizeof(double*);
           *((double*)hpcard)=scalar;           //Set scalar  
           if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
           else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2-sizeof(double); }
           }
         }
       else  
         {
         _l=size-num*len;
         if (_l<2*MIN_PQANT_pcalc_vect_norm)
           {
           for (_l=len*num;_l!=size;_l++) vect[_l]*=scalar; //Do the smallest fragment
           goto MAIN_BODY_LOOP; 
           }
         else
           {
           //Do main body first to store short cards in the end of stack
           if (!(push_pjob(pjob_id,--num+2,_pcalc_vect_scalar_product,ypsys))) return FALSE; 
           _l=num; 
           while (_l--) 
             {
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=len;      //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
             hpcard+=sizeof(double*);
             *((double*)hpcard)=scalar;           //Set scalar  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2-sizeof(double); }
             }
           _l=size-num*len;
           if (!(_l%2))
             {//Store two the same cards
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=len;      //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
             hpcard+=sizeof(double*);
             *((double*)hpcard)=scalar;           //Set scalar  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2-sizeof(double); }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=len;      //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
             hpcard+=sizeof(double*);
             *((double*)hpcard)=scalar;           //Set scalar  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2-sizeof(double); }
             }
           else
             {
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=len;      //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
             hpcard+=sizeof(double*);
             *((double*)hpcard)=scalar;           //Set scalar  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2-sizeof(double); }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=len;      //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect[len*_l]; //Set vectors fragment address   
             hpcard+=sizeof(double*);
             *((double*)hpcard)=scalar;           //Set scalar  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2-sizeof(double); }
             }
           } 
         }
       }
else if (size>=MIN_PQANT_pcalc_vect_norm*ypsys->n_threads) { len=size/ypsys->n_threads, num=ypsys->n_threads;  goto SUBMIT_JOB; }
else if (size>=MIN_PQANT_pcalc_vect_norm) { len=MIN_PQANT_pcalc_vect_norm, num=size/MIN_PQANT_pcalc_vect_norm; goto SUBMIT_JOB; }
else { _l=size; while (_l--) vect[_l]*=scalar; return NTNF; } //Do the small fragment serially
return TRUE;
}

//This function transforms set of vectors: it calculates c[_i][_k]=A[_i].A[_i+_k] for all _k: _i+1...ni. 
//Amount of cards is ni*(nj/MAX_PQANT_pcalc_vect_norm+1)+1. Sizeof cards sizeof(unsigned int)*2+sizeof(double*)*3
//This is worker function 
unsigned int _pcalc_vect_vect_dot_product(char *pcard,void *vp)
{
t_ypsystem *ypsys=(t_ypsystem*)vp;
register unsigned int pjob_id, size;
register double *vect_a, *vect_b, _d;
pjob_id=((unsigned int*)pcard)[0];                  //Set job id
size=((unsigned int*)pcard)[1];                     //Set job lenght
vect_a=((double**)(pcard+sizeof(unsigned int)*2))[0]; //Set vectors fragment address   
vect_b=((double**)(pcard+sizeof(unsigned int)*2))[1]; //Set vectors fragment address   
_d=*vect_a**vect_b; while (--size) { vect_a++, vect_b++; _d+=*vect_a**vect_b; }        //Piece of parallel job to do
if ( (pthread_mutex_lock  (&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
*(((double**)(pcard+sizeof(unsigned int)*2))[2])+=_d; //Summ vector lenght address
size=++ypsys->pjobs[pjob_id].m_pcards;                          //Standard job incremental part
if ( (pthread_mutex_unlock(&ypsys->pjobs[pjob_id].j_access))) return (unsigned int)-1; 
return size;
}
//This function calculates vectors dot product in parallel
inline char pcalc_vect_vect_dot_product(double *product,unsigned int size,double *vect_a,double *vect_b,unsigned int pjob_id,t_ypsystem *ypsys)
{
register unsigned int _l, len, num;
void* hpcard=ypsys->hpcards+ypsys->sizeof_hpcard*ypsys->_hpstore;

*product=0.;
     if (size>=MAX_PQANT_pcalc_vect_norm*ypsys->n_threads) 
       {
       len=MAX_PQANT_pcalc_vect_norm, num=size/len;
       //Pass the last two fragments
       SUBMIT_JOB:
       if (!(size%len))
         {
         //Do main body
         MAIN_BODY_LOOP:
         if (!(push_pjob(pjob_id,num,_pcalc_vect_vect_dot_product,ypsys))) return FALSE; 
         _l=num; 
         while (_l--) 
           {
           ((unsigned int*)hpcard)[0]=pjob_id;          //Set job id
           ((unsigned int*)hpcard)[1]=len;              //Set job lenght
           hpcard+=sizeof(unsigned int)*2;
           ((double**)hpcard)[0]=&vect_a[len*_l];       //Set vectors1 fragment address   
           ((double**)hpcard)[1]=&vect_b[len*_l];       //Set vectors2 fragment address  
           ((double**)hpcard)[2]=product;               //Set vector lenght address  
           if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
           else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
           }
         }
       else  
         {
         _l=size-num*len;
         if (_l<2*MIN_PQANT_pcalc_vect_norm)
           {
           for (_l=len*num;_l!=size;_l++) (*product)+=vect_a[_l]*vect_b[_l]; //Do the smallest fragment
           goto MAIN_BODY_LOOP; 
           }
         else
           {
           //Do main body first to store short cards in the end of stack
           if (!(push_pjob(pjob_id,--num+2,_pcalc_vect_vect_dot_product,ypsys))) return FALSE; 
           _l=num; 
           while (_l--) 
             {
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=len;      //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect_a[len*_l];       //Set vectors1 fragment address   
             ((double**)hpcard)[1]=&vect_b[len*_l];       //Set vectors2 fragment address  
             ((double**)hpcard)[2]=product;               //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             }
           _l=size-num*len;
           if (!(_l%2))
             {//Store two the same cards
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect_a[len*num];       //Set vectors1 fragment address   
             ((double**)hpcard)[1]=&vect_b[len*num];       //Set vectors2 fragment address  
             ((double**)hpcard)[2]=product;               //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect_a[len*num+_l];       //Set vectors1 fragment address   
             ((double**)hpcard)[1]=&vect_b[len*num+_l];       //Set vectors2 fragment address  
             ((double**)hpcard)[2]=product;               //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             }
           else
             {
             _l/=2;
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l;       //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect_a[len*num];       //Set vectors1 fragment address   
             ((double**)hpcard)[1]=&vect_b[len*num];       //Set vectors2 fragment address  
             ((double**)hpcard)[2]=product;               //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             ((unsigned int*)hpcard)[0]=pjob_id;  //Set job id
             ((unsigned int*)hpcard)[1]=_l+1;     //Set job lenght
             hpcard+=sizeof(unsigned int)*2;
             ((double**)hpcard)[0]=&vect_a[len*num+_l];       //Set vectors1 fragment address   
             ((double**)hpcard)[1]=&vect_b[len*num+_l];       //Set vectors2 fragment address  
             ((double**)hpcard)[2]=product;               //Set vector lenght address  
             if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
             else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
             }
           } 
         }
       }
else if (size>=MIN_PQANT_pcalc_vect_norm*ypsys->n_threads) { len=size/ypsys->n_threads, num=ypsys->n_threads;  goto SUBMIT_JOB; }
else if (size>=MIN_PQANT_pcalc_vect_norm) { len=MIN_PQANT_pcalc_vect_norm, num=size/MIN_PQANT_pcalc_vect_norm; goto SUBMIT_JOB; }
else { _l=size; while (_l--) (*product)+=vect_a[_l]*vect_b[_l]; return NTNF; } //Do the small fragment serially
return TRUE;
}




/*
The matrix:
,-------------------. 
|                   |  A - a block
'-------------------'
,-------------------. 
|                   |  B - a block
'-------------------'
,-------------------.
|                   |
|                   |  C - a matrix
|                   |
'-------------------'
The logic:
....................................
 A[0]*             |       |     
 |                 |       |
 '-> B[0] + C[0] <-'       |
     |       '-------------+---------.
     '-> A[1]*             |         |
         |                 |         |
         '-> B[1] + C[1] <-'         |
              |      '---------------+---------.
              '-> A[2]*              |         |
                   |                 |         |
                   '-> B[2] + C[2] <-'         |
                        |      '---------------+--------.
                        '-> A[3]*              |        |
                        ..................................... so on
The algorithm;
HANDLE_BEGIN_OF_MATRIX;
_i=ni-2;
while (--_i>MIN_BLOCK_SIZE)
  {
  job_C-2=job_C-1, job_C-1=job_C0, job_B-1=job_B0;
                                  job_A0=calc_hpA(i);
  if (await_job_finished(job_A0)) job_B0=calc_hpB(i);
  if ( (await_job_finished(job_B0))&&( (job_C-2==-1)||(await_job_finished(job_C-2)) ) )
    job_C0=calc_pC(i);
  }
HANDLE_END_OF_MATRIX;

*/


//This function calculates product of submatrix [i...ni,j...nj] with a vector p=A.y in parallel
inline char pcalc_banded_dmatrix_origin_vect_product(unsigned int i,unsigned int j,unsigned int ni,unsigned int nj,double **A,double *y,double *p,unsigned int pjob_id,t_ypsystem *ypsys)
{
register unsigned int _l, len, num;
register double _d, *vect_a, *vect_b;
void* hpcard=ypsys->hpcards+ypsys->sizeof_hpcard*ypsys->_hpstore;

if (nj-j>=MAX_PQANT_pcalc_vect_norm) 
  {
  len=MAX_PQANT_pcalc_vect_norm, num=(nj-j)/len;
  //Pass the last two fragments
  if (!((nj-j)%MAX_PQANT_pcalc_vect_norm))
    {
    if (!(push_pjob(pjob_id,num*(ni-i),_pcalc_vect_vect_dot_product,ypsys))) return FALSE; 
    while (--ni!=i)
      {
      p[ni]=0.; 
      _l=num; 
      while (_l--) 
        {
        ((unsigned int*)hpcard)[0]=pjob_id;         //Set job id
        ((unsigned int*)hpcard)[1]=len;             //Set job lenght
        hpcard+=sizeof(unsigned int)*2;
        ((double**)hpcard)[0]=&A [i][_l*num+j];     //Set vectors1 fragment address   
        ((double**)hpcard)[1]=&A[ni][_l*num+j];     //Set vectors2 fragment address  
        ((double**)hpcard)[2]=&p[ni];               //Set vector lenght address  
        if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
        else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
        }
      }
    }
  else 
    {
    if (!(push_pjob(pjob_id,(--num+2)*(ni-i),_pcalc_vect_vect_dot_product,ypsys))) return FALSE; 
    while (--ni!=i)
      {
      p[ni]=0.; 
      _l=num; 
      while (_l--) 
        {
        ((unsigned int*)hpcard)[0]=pjob_id;         //Set job id
        ((unsigned int*)hpcard)[1]=len;             //Set job lenght
        hpcard+=sizeof(unsigned int)*2;
        ((double**)hpcard)[0]=&A [i][_l*num+j];     //Set vectors1 fragment address   
        ((double**)hpcard)[1]=&A[ni][_l*num+j];     //Set vectors2 fragment address  
        ((double**)hpcard)[2]=&p[ni];               //Set vector lenght address  
        if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
        else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
        }
      }
    if (!((nj-j-num*len)%2))
      {
      _l=(nj-j-num*len)/2;
      while (--ni!=i)
        {
        //Store two identical terminal cards
        ((unsigned int*)hpcard)[0]=pjob_id;         //Set job id
        ((unsigned int*)hpcard)[1]=_l;              //Set job lenght
        hpcard+=sizeof(unsigned int)*2;
        ((double**)hpcard)[0]=&A [i][j+len*num];    //Set vectors1 fragment address   
        ((double**)hpcard)[1]=&A[ni][j+len*num];    //Set vectors2 fragment address  
        ((double**)hpcard)[2]=&p[ni];               //Set vector lenght address  
        if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
        else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
        ((unsigned int*)hpcard)[0]=pjob_id;         //Set job id
        ((unsigned int*)hpcard)[1]=_l;              //Set job lenght
        hpcard+=sizeof(unsigned int)*2;
        ((double**)hpcard)[0]=&A [i][j+_l+len*num]; //Set vectors1 fragment address   
        ((double**)hpcard)[1]=&A[ni][j+_l+len*num]; //Set vectors2 fragment address  
        ((double**)hpcard)[2]=&p[ni];               //Set vector lenght address  
        if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
        else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
        }
      }
    else
      {
      _l=(nj-j-num*len)/2;
      while (--ni!=i)
        {
        //Store two similar terminal cards
        ((unsigned int*)hpcard)[0]=pjob_id;         //Set job id
        ((unsigned int*)hpcard)[1]=_l;              //Set job lenght
        hpcard+=sizeof(unsigned int)*2;
        ((double**)hpcard)[0]=&A [i][j+len*num];    //Set vectors1 fragment address   
        ((double**)hpcard)[1]=&A[ni][j+len*num];    //Set vectors2 fragment address  
        ((double**)hpcard)[2]=&p[ni];               //Set vector lenght address  
        if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
        else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
        ((unsigned int*)hpcard)[0]=pjob_id;         //Set job id
        ((unsigned int*)hpcard)[1]=_l+1;            //Set job lenght
        hpcard+=sizeof(unsigned int)*2;
        ((double**)hpcard)[0]=&A [i][j+_l+len*num]; //Set vectors1 fragment address   
        ((double**)hpcard)[1]=&A[ni][j+_l+len*num]; //Set vectors2 fragment address  
        ((double**)hpcard)[2]=&p[ni];               //Set vector lenght address  
        if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
        else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
        }
      }
    } 
  } 
else
  {
  if ( (nj-j>=MIN_PQANT_pcalc_vect_scalar_product))
    {
    len=nj-j;
    if (!(push_pjob(pjob_id,ni-i,_pcalc_vect_vect_dot_product,ypsys))) return FALSE; 
    while (--ni!=i)
      {
      p[ni]=0.; 
      ((unsigned int*)hpcard)[0]=pjob_id;         //Set job id
      ((unsigned int*)hpcard)[1]=len;             //Set job lenght
      hpcard+=sizeof(unsigned int)*2;
      ((double**)hpcard)[0]=&A [i][j];            //Set vectors1 fragment address   
      ((double**)hpcard)[1]=&A[ni][j];            //Set vectors2 fragment address  
      ((double**)hpcard)[2]=&p[ni];               //Set vector lenght address  
      if (ypsys->_hpstore==ypsys->n_hpcards-1) { ypsys->_hpstore=0, hpcard=ypsys->hpcards; }
      else { ypsys->_hpstore++, hpcard+=ypsys->sizeof_hpcard-sizeof(unsigned int)*2; }
      }
    }
  else 
    {//Fragments are too small, do them serially
    while (--ni!=i)
      {
      _d=0., vect_a=&A[i][j], vect_b=&A[ni][j], _l=nj-j;
      while (nj!=_l)
        {
        _d+=*vect_a**vect_b;
        vect_a++, vect_b++, _l++;
        }
      p[ni]=_d; 
      }
    return NTNF;
    }
  }
return TRUE;
}


/*
//This is the leader function of parallel QL factorization: A=L.Q and Q=I-Z.T.Z^T, where A - nixnj (ni<=nj), L - lower triangle njxnj, T - lower triangle njxnj and Z is upper-trapezoidal stored in zeroed entries of A.
char pLQ_factorization(unsigned int ni,unsigned int nj,unsigned int nr,double **A,double **T,t_ypsystem *ypsystem)
{
unsigned int mi, mj;
double u, c[BLOCK_QUANT_pQL_factorization][];
//Stage I. Setup algorithm

//mi=(!((ni-1)%MAX_PQUANT)) ? (ni-1)/MAX_PQUANT_pQL_factorization : (ni-1)/MAX_PQUANT_pQL_factorization+1;
//mj=(!((nj-1)%MAX_PQUANT)) ? (nj-1)/MAX_PQUANT_pQL_factorization : (nj-1)/MAX_PQUANT_pQL_factorization+1;
//if (!setup_ypsys(ypsys,4,2*mi,2*mi*mj,sizeof(unsigned int)*2+sizeof(*double)*2,sizeof(unsigned int)*8+sizeof(double**))) return FALSE;

//Stage II. Do palgorithm
//Stage II. 1. Scann blocks
for (i=0;i<nj-
  {
  //Prepare next block
  block_size=0;
  while ( (BLOCK_QUANT_pQL_factorization!=nr--)&&(BLOCK_QUANT_pQL_factorization!=block_size) ) 
    {//Gather BLOCK transformations vector by vector
    if ( (pcalc_vect_norm(&u,nj-i,&A[i][i],0,ypsys)!=NTNF)&&(wait_pjob_done(0,ypsys)) ) goto LABEL_ERROR; //Calc vect norm
    if (u<TINY) continue; 
    else { u=sqrt(u)*scale; A[nr][nr]*=scale; }
    if (A[i][i]>0.) { diag[i]=A[i][i]-=u); u=1./sqrt(2*u*(u-A[i][i])); }
    else            { diag[i]=A[i][i]+=u); u=1./sqrt(2*u*(u+A[i][i])); }
    if ( (pcalc_vect_scalar_product(u,nj-nr,&A[i][i],0,ypsys)!=NTNF)&&(wait_pjob_done(0,ypsys)) )     goto LABEL_ERROR; //Normalize vector
    //Calc contants
    if ( (pcalc_banded_dmatrix_origin_vect_product(i+1,i,ni,nj,A,A[i],&c[i][i+1],0,ypsys)!=NTNF)&&(wait_pjob_done(0,ypsys)) ) goto LABEL_ERROR; //Transform all vectors in block
    }
  //Store C matrix 
  }
*/

//hpcard_id=0;
//job_CP=job_DP=job_DPP=(unsigned int)-1; 
//Stage II. 2. Do main body
//_i=ni;
//while (--_i)
//  {//Scan row-by-row, but skip the bottomest vector
//
//  //Submit job_A to stack 
//
//  //Submit job_B to stack
//  }
//
//Stage II. 3. Do post-calculations
//
//Stage III. Exit sucessfully
//return TRUE;
//}





