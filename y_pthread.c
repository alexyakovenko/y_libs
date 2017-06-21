#include "y_pthread.h"

//The BROADCAST means unblocking of threads sleeping at conditional variable:
#define ypsem_BROADCAST_NO    FALSE //Nobody
#define ypsem_BROADCAST_ONE   TRUE  //Only one thread
#define ypsem_BROADCAST_ALL   NTNF  //All threads

//This function initialize ypsem
inline char ypsem_init(t_ypsem *ypsem)
{
ypsem->v=0;
if ( ( (pthread_mutex_init(&ypsem->m,NULL)))||( (pthread_cond_init(&ypsem->c,NULL)) ) ) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
return TRUE;
}
//This function delete ypsem
inline char ypsem_destroy(t_ypsem *ypsem)
{
if ( ( (pthread_mutex_destroy(&ypsem->m)))||( (pthread_cond_destroy(&ypsem->c))) ) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
return TRUE;
}

//This function increase condition variable by value
inline char ypsem_incr(register t_ypsem *ypsem,register unsigned int value,register int broadcast,register unsigned int trigger)
{
if ( (pthread_mutex_lock(&ypsem->m))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
if (((unsigned int)-1)-ypsem->v>=value) ypsem->v+=value; else { ylib_errno=YERROR_IMPOSSIBLE; return FALSE; }
if ( ( (broadcast))&&(ypsem->v>=trigger) )
  {
       if (broadcast==ypsem_BROADCAST_ONE) { if ( (pthread_cond_signal(&ypsem->c)))    { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; } }
  else if (broadcast==ypsem_BROADCAST_ALL) { if ( (pthread_cond_broadcast(&ypsem->c))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; } }
  else                                                                                 { ylib_errno=YERROR_NIMPLEMENTED;  return FALSE; } 
  }
if ( (pthread_mutex_unlock(&ypsem->m))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
return TRUE;
}
//This function decrease condition variable by value
inline char ypsem_decr(register t_ypsem *ypsem,register unsigned int value,register int  broadcast,register unsigned int trigger)
{
if ( (pthread_mutex_lock(&ypsem->m)))   { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
if (ypsem->v>=value) ypsem->v-=value; else { ylib_errno=YERROR_IMPOSSIBLE; return FALSE; }
if ( ( (broadcast))&&(ypsem->v<=trigger) )
  {
       if (broadcast==ypsem_BROADCAST_ONE) { if ( (pthread_cond_signal(&ypsem->c)))    { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; } }
  else if (broadcast==ypsem_BROADCAST_ALL) { if ( (pthread_cond_broadcast(&ypsem->c))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; } }
  else                                                                                 { ylib_errno=YERROR_NIMPLEMENTED;  return FALSE; } 
  }
if ( (pthread_mutex_unlock(&ypsem->m))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
return TRUE;
}

//These are waiting functions
inline char ypsem_wait_value(register t_ypsem *ypsem,register unsigned int trigger)
{
if ( (pthread_mutex_lock(&ypsem->m))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
while (ypsem->v<trigger) if ( (pthread_cond_wait(&ypsem->c,&ypsem->m))) { pthread_mutex_unlock(&ypsem->m); ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
if ( (pthread_mutex_unlock(&ypsem->m))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
return TRUE;
}
//The same as above but with trigger == zero 
inline char ypsem_wait__zero(register t_ypsem *ypsem)
{
if ( (pthread_mutex_lock(&ypsem->m))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
while ( (ypsem->v)) if ( (pthread_cond_wait(&ypsem->c,&ypsem->m))) { pthread_mutex_unlock(&ypsem->m); ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
if ( (pthread_mutex_unlock(&ypsem->m))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
return TRUE;
}

//This function activates a group of sleepeing threads and waits while they do their job
char master_ypsem_activate_and_wait(unsigned int nthreads,t_ypsem *ypsem_s,t_ypsem *ypsem_c)
{
return ( ( (ypsem_incr(ypsem_s,nthreads,ypsem_BROADCAST_ALL,nthreads)))&&( (ypsem_wait__zero(ypsem_s)))&&
         ( (ypsem_incr(ypsem_c,nthreads,ypsem_BROADCAST_ALL,nthreads)))&&( (ypsem_wait__zero(ypsem_c))) );
}
//These two functions are the split of previous one (for nonblocking runs)
char master_ypsem_activate(unsigned int nthreads,t_ypsem *ypsem_s)
{
return ( (ypsem_incr(ypsem_s,nthreads,ypsem_BROADCAST_ALL,nthreads)));
} 
char master_ypsem_wait(unsigned int nthreads,t_ypsem *ypsem_s,t_ypsem *ypsem_c)
{
return ( ( (ypsem_wait__zero(ypsem_s)))&&( (ypsem_incr(ypsem_c,nthreads,ypsem_BROADCAST_ALL,nthreads)))&&( (ypsem_wait__zero(ypsem_c))) );
}



//--------------------------------------------------   E N D   O F   P Y S E M   F U N C T I O N S   --------------------------------------------------//


//This function does nothing, it is just a flag to keep a thread
char ypthread_function_nope(unsigned int worker_id,t_ypthread *yphread)
{
return TRUE;
}
//This function does nothing just a flag to quit a thread
char ypthread_function_quit(unsigned int worker_id,t_ypthread *ypthread)
{
extern unsigned int ylib_errno;
ylib_errno=YERROR_IMPOSSIBLE; return FALSE; //This function should never be executed!
}

//This function free ypthreads
void free_ypthread(t_ypthread *ypthread)
{
if ( (ypthread))
  {
  //Stop workers
  if (ypthread->threads) finish_ypthreads(ypthread); 
  //Synchronization stuff
  ypsem_destroy(&ypthread->ypsem_s);
  ypsem_destroy(&ypthread->ypsem_c);
  pthread_mutex_destroy(&ypthread->mutex);
  pthread_spin_destroy(&ypthread->spinlock); 
  free(ypthread);
  }
}

//This function initilaize pthreads
t_ypthread *init_ypthread(unsigned int nthreads)
{
t_ypthread *ypthread;
//Allocate memory
if (!(ypthread=(t_ypthread*)malloc(sizeof(t_ypthread)+sizeof(pthread_t)*nthreads))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { ypthread->threads=(pthread_t*)((void*)ypthread+sizeof(t_ypthread)); ypthread->nthreads=nthreads, ypthread->funct=ypthread_function_nope; }
//Init sub-structures
if (!(ypsem_init(&ypthread->ypsem_s))) { LABEL_ERROR_0: free(ypthread); ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
if (!(ypsem_init(&ypthread->ypsem_c))) { LABEL_ERROR_1: ypsem_destroy(&ypthread->ypsem_s); goto LABEL_ERROR_0; }
if ( (pthread_mutex_init(&ypthread->mutex,NULL))) { LABEL_ERROR_2: ypsem_destroy(&ypthread->ypsem_c); goto LABEL_ERROR_1; }
if ( (pthread_spin_init(&ypthread->spinlock,PTHREAD_PROCESS_PRIVATE))) { pthread_mutex_destroy(&ypthread->mutex); goto LABEL_ERROR_2; }
return ypthread;
}

//This is the main function of ypthreads
void *ypthread_main(void *arg)
{
unsigned int worker_id;
pthread_t self_pthread;
t_ypthread *ypthread=(t_ypthread*)arg;

//Stage I. Get own worker_id
self_pthread=pthread_self();
worker_id=0; while (!(pthread_equal(self_pthread,ypthread->threads[worker_id]))) worker_id++;
//Stage II. Do job & sleep while asked
while ( ( (ypsem_wait_value(&ypthread->ypsem_s,1)))&&(ypthread->funct!=ypthread_function_quit) )
  if ( (!(ypsem_decr(&ypthread->ypsem_s,1,ypsem_BROADCAST_ONE,0)))||(!(ypthread->funct(worker_id,ypthread)))||
       (!(ypsem_wait_value(&ypthread->ypsem_c,1)))||(!(ypsem_decr(&ypthread->ypsem_c,1,ypsem_BROADCAST_ONE,0))) ) 
    exit(EXIT_FAILURE); //VERY unsuccessful thread termination
return FALSE; //Exit OK
}

//This function runs numbrer of threads and put them on hold
char launch_ypthreads(t_ypthread *ypthread)
{
register unsigned int _i, _j;
ypthread->funct=ypthread_function_nope;
for (_i=_j=0; _i<ypthread->nthreads; _i++)
  if (!(pthread_create(&ypthread->threads[_j],NULL,ypthread_main,(void*)ypthread))) 
    _j++;
if (!(_j))
  { yprintf(YPRINTF_ERROR,"Failed to create even a single worker thread, forced to quit.\n"); ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
else
  if (_j!=ypthread->nthreads)
    {//Run in not-full-thread mode
     yprintf(YPRINTF_WARNING,"Failed to create more than %d threads, working with whatever I have.\n",_j);
     ypthread->nthreads=_j;
     }
return TRUE;
}

//This function execute (*funct) with pthreads in parallel (blocking the calling thread!) 
inline char exec_client_ypthreads(t_ypthread *ypthread,char (*funct)(unsigned int ,struct ypthread_t*))
{
ypthread->funct=funct; return ( (master_ypsem_activate_and_wait(ypthread->nthreads,&ypthread->ypsem_s,&ypthread->ypsem_c)));
}
//These two functions are a nonblocking split of the previous
inline char exec_client_ypthreads_nonblock(t_ypthread *ypthread,char (*funct)(unsigned int ,struct ypthread_t*))
{
ypthread->funct=funct; return master_ypsem_activate(ypthread->nthreads,&ypthread->ypsem_s);
}
inline char wait_exec_client_ypthreads(t_ypthread *ypthread)
{
return master_ypsem_wait(ypthread->nthreads,&ypthread->ypsem_s,&ypthread->ypsem_c);
}


//This function terminates launched threads
char finish_ypthreads(t_ypthread *ypthread)
{
register unsigned int _i;
ypthread->funct=ypthread_function_quit;
if (!(ypsem_incr(&ypthread->ypsem_s,ypthread->nthreads,ypsem_BROADCAST_ALL,ypthread->nthreads))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
_i=ypthread->nthreads; while (_i--) if ( (pthread_join(ypthread->threads[_i],NULL))) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
return TRUE;
}


