#include "y_plibs.h"
#include <errno.h>
#include <sys/types.h>
extern unsigned int ylib_errno;

struct sembuf sem_wait={  0,  0,  0 };
struct sembuf sem_decr={  0, -1, IPC_NOWAIT };
struct sembuf sem_incr={  0, +1, IPC_NOWAIT };
struct sembuf sem_decr_wait={  0, -1, 0 };

union {int val; struct semid_ds *buff; unsigned short *array; struct seminfo *_buf;} semun;

//This function makes threads pool and inits structures
t_ypsystem *init_ypsystem(unsigned int nthreads)
{
unsigned int _i;
t_ypsystem *ypsys;

//Set defaults
if (!(nthreads)) { ylib_errno=YERROR_USER; return FALSE; }
//Alloc memory
if (!(ypsys=malloc(sizeof(t_ypsystem)+sizeof(pthread_t)*nthreads))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else 
  {//Init memory
  ypsys->n_threads=nthreads;
  ypsys->t_ids=(void*)ypsys+sizeof(t_ypsystem);
  ypsys->n_pjobs=0;
  ypsys->_n_pjobs=0;
  ypsys->n_hpcards=0;
  ypsys->sizeof_hpcard=0;
  ypsys->_sizeof_hpcards=0;
  ypsys->hpcards=0x0;
  ypsys->pjobs=0x0;
  ypsys->n_pcards=0;
  ypsys->sizeof_pcard=0;
  ypsys->_sizeof_pcards=0;
  ypsys->pcards=0x0;
  }
//Setup semaphore and mutexes
semun.val=0; //Init threads off
if ((ypsys->sem_threads=semget(IPC_PRIVATE,1,IPC_CREAT|IPC_EXCL|0660))==-1) { LABEL_ERROR_1: ylib_errno=YERROR_EXTERNAL_CODE; free(ypsys); return FALSE; }
else if (semctl(ypsys->sem_threads,0,SETVAL,semun)==-1) { LABEL_ERROR_2: semctl(ypsys->sem_threads,0,IPC_RMID); goto LABEL_ERROR_1; } //Init semaphore
 ypsys->sem_pjobs=(unsigned int)-1;
if ( (pthread_mutex_init(&ypsys-> c_access,NULL))) goto LABEL_ERROR_2;			  
if ( (pthread_mutex_init(&ypsys->hc_access,NULL))) { pthread_mutex_destroy(&ypsys-> c_access); goto LABEL_ERROR_2; }			  
//Start threads
_i=nthreads;
while (_i--)
  if ( (pthread_create(&ypsys->t_ids[_i],0x0,pthread_worker,(void*)ypsys)))
    {
    while (++_i<nthreads) pthread_cancel(ypsys->t_ids[_i]);
    pthread_mutex_destroy(&ypsys-> c_access);
    pthread_mutex_destroy(&ypsys->hc_access);
    goto LABEL_ERROR_2;
    } 
//Exit successfuly
return ypsys;
}

//This function sets up the parallel data structure for particular algorithm
char setup_ypsys(t_ypsystem *ypsys,unsigned int n_pjobs,unsigned int n_hpcards,unsigned int n_pcards,size_t sizeof_hpcard,size_t sizeof_pcard)
{
void *vp;
if ( ((n_hpcards)&&(sizeof_hpcard<sizeof(unsigned int)))||((n_pcards)&&(sizeof_pcard<sizeof(unsigned int)))||(n_pjobs>n_hpcards+n_pcards) )
  { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; } //At least job id is required (if changed affect space for SETALL job semphore 
//Update pjobs
if (ypsys->_n_pjobs<n_pjobs) 
  {
  if (!(vp=realloc(ypsys->pjobs,sizeof(t_ypjob)*n_pjobs))) { LABEL_MEMORY_ERROR: ylib_errno=YERROR_MEMORY; return FALSE; }
  else ypsys->pjobs=(t_ypjob*)vp;
  if (ypsys->sem_pjobs!=(unsigned int)-1)
    {
    if (semctl(ypsys->sem_pjobs,0,IPC_RMID)==1) goto LABEL_EXTERNAL_ERROR;
    else ypsys->sem_pjobs=(unsigned int)-1;
    }
  if ((ypsys->sem_pjobs=semget(IPC_PRIVATE,n_pjobs,IPC_CREAT|IPC_EXCL|0660))==-1) { LABEL_EXTERNAL_ERROR: ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
  do{
    ypsys->pjobs[ypsys->_n_pjobs].n_pcards=ypsys->pjobs[ypsys->_n_pjobs].m_pcards=0;
    if ( (pthread_mutex_init(&ypsys->pjobs[ypsys->_n_pjobs].j_access,NULL))) goto LABEL_EXTERNAL_ERROR;  
    }while (++ypsys->_n_pjobs!=n_pjobs); 
  }
else ypsys->n_pjobs=n_pjobs;
//Update hpcards
if ( (!ypsys->hpcards)||(ypsys->_sizeof_hpcards<n_hpcards*sizeof_hpcard) ) 
  {
  if (ypsys->hpcards) { free(ypsys->hpcards); ypsys->hpcards=0x0; }
  if (!(ypsys->hpcards=(char*)malloc(sizeof_hpcard*n_hpcards))) goto LABEL_MEMORY_ERROR;
  else ypsys->_sizeof_hpcards=n_hpcards*sizeof_hpcard;
  }
ypsys->n_hpcards=n_hpcards;
ypsys->sizeof_hpcard=sizeof_hpcard;
//Update pcards
if ( (!ypsys->pcards)||(ypsys->_sizeof_pcards<n_pcards*sizeof_pcard) ) 
  {
  if (ypsys->pcards) { free(ypsys->pcards); ypsys->pcards=0x0; }
  if (!(ypsys->pcards=(char*)malloc(sizeof_pcard*n_pcards))) goto LABEL_MEMORY_ERROR;
  else ypsys->_sizeof_pcards=n_pcards*sizeof_pcard;
  }
ypsys->n_pcards=n_pcards;
ypsys->sizeof_pcard=sizeof_pcard;
//Update the semaphore
     if (n_hpcards>=n_pjobs/2) semun.array=ypsys->hpcards;
else if ( n_pcards>=n_pjobs/2) semun.array=ypsys->pcards;
else
  {
  if (!(semun.array=calloc(n_pjobs,sizeof(short)))) goto LABEL_MEMORY_ERROR;
  if (semctl(ypsys->sem_pjobs,n_pjobs,SETALL,semun)==-1) { semctl(ypsys->sem_pjobs,0,IPC_RMID); ypsys->sem_pjobs=-1; goto LABEL_EXTERNAL_ERROR; }
  else free(semun.array);
  return TRUE;  
  }
while(n_pjobs--) semun.array[n_pjobs]=0;
if (semctl(ypsys->sem_pjobs,ypsys->n_pjobs,SETALL,semun)==-1) { semctl(ypsys->sem_pjobs,0,IPC_RMID); ypsys->sem_pjobs=-1; goto LABEL_EXTERNAL_ERROR; }
return TRUE;
}

//This function do work in ypsys
void *pthread_worker(void *vp)
{
t_ypsystem *ypsys=(t_ypsystem*)vp;
int pcard_id, pjob_id, n_pcards, m_pcards, _i, _j;
struct sembuf sembuf_pjob_free_leader={  0,  -1,  IPC_NOWAIT };
struct sembuf sembuf_worker_wait[2]={ {  0,  -1,  0 }, {  0, +1,  IPC_NOWAIT } };
//Infinite wait for new jobs
INFINITE_LOOP:
while (semop(ypsys->sem_threads,sembuf_worker_wait,2)==-1) if ( (errno!=EAGAIN)&&(errno!=EINTR) ) ERROR_LABEL: return 0x0;
//There is some job for me, do it
//Check high-priority jobs first
//printf("I am awaked %d worker\n",pthread_self());
FINITE_LOOP:
if (ypsys->_hpexec!=ypsys->_hpstore)
  {
  if ( (pthread_mutex_lock(&ypsys->hc_access))) goto ERROR_LABEL;
  if (ypsys->_hpexec!=ypsys->_hpstore)
    {//High-priority job is probably there - execute it
    pcard_id=ypsys->_hpexec;
    if (++ypsys->_hpexec==ypsys->n_hpcards) ypsys->_hpexec=0;
    if ( (pthread_mutex_unlock(&ypsys->hc_access))) goto ERROR_LABEL; //Release mutex as fast as posible
    pjob_id=*(unsigned int*)((void*)ypsys->hpcards+ypsys->sizeof_hpcard*pcard_id); 
    n_pcards=ypsys->pjobs[pjob_id].n_pcards;
    if ((m_pcards=ypsys->pjobs[pjob_id].funct(ypsys->hpcards+ypsys->sizeof_hpcard*pcard_id,ypsys))==-1) goto ERROR_LABEL; //Worker code has crashed 
    PROCESS_PJOB:
//    printf("I am a %d worker, did %d job card\n",pthread_self(),m_pcards);
    if (m_pcards==n_pcards)
      {//I have to free release the job slot
      //Relize semaphore to leader and remove the job from scheduller
      sembuf_pjob_free_leader.sem_num=pjob_id;
      while (semop(ypsys->sem_threads,&sem_decr,1)==-1)              if ( (errno!=EAGAIN)&&(errno!=EINTR) ) goto ERROR_LABEL;  
      while (semop(ypsys->sem_pjobs,&sembuf_pjob_free_leader,1)==-1) if ( (errno!=EAGAIN)&&(errno!=EINTR) ) goto ERROR_LABEL;
      }
    goto FINITE_LOOP;
    }
  else 
    {
    if ( (pthread_mutex_unlock(&ypsys->hc_access))) goto ERROR_LABEL; //Release mutex as fast as posible
    else goto CHECK_LOW_PRIORITY_JOBS;
    }
  }
else
  {//Try to exec low-priority job
  CHECK_LOW_PRIORITY_JOBS:
  if (ypsys->_pexec!=ypsys->_pstore)
    {
    if ( (pthread_mutex_lock(&ypsys->c_access))) goto ERROR_LABEL; 
    if (ypsys->_pexec!=ypsys->_pstore)
      {//High-priority job is probably there - execute it
      pcard_id=ypsys->_pexec;
      if (++ypsys->_pexec==ypsys->n_pcards) ypsys->_pexec=0;
      if ( (pthread_mutex_unlock(&ypsys->c_access))) goto ERROR_LABEL; //Release mutex as fast as posible
      pjob_id=*(unsigned int*)((void*)ypsys->pcards+ypsys->sizeof_pcard*pcard_id); 
      n_pcards=ypsys->pjobs[pjob_id].n_pcards;
      if ((m_pcards=ypsys->pjobs[pjob_id].funct(ypsys->pcards+ypsys->sizeof_pcard*pcard_id,ypsys))==-1) goto ERROR_LABEL; //Worker code has crashed 
      goto PROCESS_PJOB; 
      }
    else 
      if ( (pthread_mutex_unlock(&ypsys->c_access))) goto ERROR_LABEL; //Release mutex as fast as posible
    }
  }
goto INFINITE_LOOP; //The thread must be removed from pull with a signal
return 0x0;
}

//This function push new job into the stack
inline char push_pjob(unsigned int pjob_id,unsigned int n_cards,unsigned int (*funct)(char *card,void *vp),t_ypsystem *ypsys)
{
struct sembuf sem_incr={ 0, +1, IPC_NOWAIT };
ypsys->pjobs[pjob_id].n_pcards=n_cards, ypsys->pjobs[pjob_id].m_pcards=0, ypsys->pjobs[pjob_id].funct=funct;
semun.val=+1;
while (semctl(ypsys->sem_pjobs,pjob_id,SETVAL,semun)==-1) if ( (errno!=EAGAIN)&&(errno!=EINTR) ) { ERROR_EXIT: ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
while (semop(ypsys->sem_threads,&sem_incr,1)==-1) if ( (errno!=EAGAIN)&&(errno!=EINTR) ) goto ERROR_EXIT;
return TRUE;
}

//This function awaits the submitted job is done
inline char wait_pjob_done(unsigned int pjob_id,t_ypsystem *ypsys)
{
int _i,_j,_k,_l;
struct sembuf sembuf_pjob_wait;
sembuf_pjob_wait.sem_num=(unsigned short)pjob_id, sembuf_pjob_wait.sem_op=0, sembuf_pjob_wait.sem_flg=0;
if (pjob_id!=(unsigned int)-1) {
while (semop(ypsys->sem_pjobs,&sembuf_pjob_wait,1)==-1) if ( (errno!=EAGAIN)&&(errno!=EINTR) ) return FALSE;
return TRUE;
              }
_i=semctl(ypsys->sem_pjobs,pjob_id,GETVAL);
_j=semctl(ypsys->sem_threads,0,GETVAL);
_k=pthread_mutex_lock(&ypsys->hc_access);
return FALSE;
}


