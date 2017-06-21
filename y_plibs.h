#define Y_PLIBS 0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif

//Semaphores stuff
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>

//Pthreads stuff
#include <pthread.h>


//This module do parallel implemenation of ylibs

//One job description
typedef struct{
              pthread_mutex_t                 j_access;     //Mutex to control access to jobs resources
              unsigned int                    n_pcards;     //Total amount of asked cards in jobs infrastructure
              unsigned int                    m_pcards;     //Total amount of processed cards in jobs infrastrucute
              unsigned int (*funct)(char *pcard,void *ypsys);  //Function to provide calculations
              }t_ypjob;


//general parallel structure
typedef struct{
              //Threads
              int                         sem_threads;     //Semaphore to start jobs (up to 32768-1 is allowed)
              unsigned int                  n_threads;     //Amount of threads
              pthread_t                        *t_ids;     //Thread ids
              //Regular jobs
              int                           sem_pjobs;     //Job acessment semaphore  
              unsigned int                    n_pjobs;     //Amount of active pjobs <=_n_pjobs
              unsigned int                   _n_pjobs;     //Amount of allowed pjobs <=_n_hpcards+n_pcards
              t_ypjob                          *pjobs;     //Pointer to parallel jobs storage
              //Cards
              pthread_mutex_t               hc_access;     //Mutex to control access to jobs cards
              unsigned int                  n_hpcards;     //Amount of cards to store the job <= _n_pcards 
              void                           *hpcards;     //Cards for parallel jobs       
              pthread_mutex_t                c_access;     //Mutex to control access to high-priority jobs cards                    
              unsigned int                   n_pcards;     //Amount of high-priority cards to store the job <= _n_hpcards 
              void                            *pcards;     //High-priority cards for parallel 
              size_t      sizeof_pcard, sizeof_hpcard;     //Lenght of a card 
              size_t  _sizeof_pcards, _sizeof_hpcards;     //amount of memory alloced for cards        
              unsigned int  _hpexec, _hpstore, _pexec, _pstore;  //iterators to store and execute new cards
              }t_ypsystem;

//Global structure pointer
t_ypsystem *ypsys;

//This function makes threads pool and inits structures
t_ypsystem *init_ypsystem(unsigned int nthreads);

//This function sets up the parallel data structure for particular algorithm
char setup_ypsys(t_ypsystem *ypsys,unsigned int n_pjobs,unsigned int n_hpcards,unsigned int n_pcards,size_t sizeof_hpcard,size_t sizeof_pcard);

//This function do work in ypsys
void *pthread_worker(void *vp);

//This function push new job into the stack
inline char push_pjob(unsigned int pjob_id,unsigned int n_cards,unsigned int (*funct)(char *card,void *vp),t_ypsystem *ypsys);

//This function awaits the submitted job is done
inline char wait_pjob_done(unsigned int pjob_id,t_ypsystem *ypsys);

