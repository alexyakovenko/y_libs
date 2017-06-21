
#include "y_lib.h"

//This module do parallel implemenation of ylibs

//Card description
#define YPCARD_LEN 0xFF

//One job description
typedef struct{
              unsigned int                    n_pcards;    //Total amount of asked cards in jobs infrastructure
              unsigned int                    m_pcards;    //Total amount of processed cards in jobs infrastrucute
              char (*funct)(char (*pcard)[YPCARD_LEN]);   //Function to provide calculations
              }t_ypjob;

//Thread description
typedef struct{
              pthread_t                           t_id;    //Thread id
              unsigned int                      job_id;    //Job id thread is executing at the moment (-1 means sleeping)
              }t_ypthread;

//general parallel structure
typedef struct{
              int                              sem_id;     //Semaphore to start jobs (up to 32768-1 is allowed)
              //Threads
              unsigned int                 n_pthreads;     //Amount of threads
              t_ypthread                    *pthreads;     //Pointer to threads storage
              //Regular jobs
              pthread_mutex_t                j_access;     //Job acessment   
              unsigned int                    n_pjobs;     //Amount of active pjobs <=_n_pjobs
              unsigned int                   _n_pjobs;     //Amount of allowed pjobs <=_n_hpcards+n_pcards        .default= 4*n_pthreads 
              t_ypjob                          *pjobs;     //Pointer to parallel jobs storage
              //Cards
              pthread_mutex_t               hc_access;     //Mutex to control access to jobs cards
              unsigned int                   n_pcards;     //Amount of cards to store the job <= _n_pcards 
              unsigned int                  _n_pcards;     //Maximal amount of job cards in system                                 -.   .default= 32*n_pthreads
              char              (*pcards)[YPCARD_LEN];     //Cards for parallel jobs                                                |
              pthread_mutex_t              hpc_access;     //Mutex to control access to high-priority jobs cards                    |-> together <=32768-1
              unsigned int                  n_hpcards;     //Amount of high-priority cards to store the job <= _n_hpcards           |
              unsigned int                 _n_hpcards;     //Maximal amount of high-priority job cards in system  jobs             -'   .default=  2*n_pthreads
              char             (*hpcards)[YPCARD_LEN];     //High-priority cards for parallel 
              }t_ypsystem;

/*
The logic:
0. Worker sleeps                                                                     | 0. Commander setup the job (determine m_pcards, set_mpcards=0 and set funct to work on)
1. Worker wakes up and pickup the last card                                          | 1. Commander generates cards and updates sem_id value by one after each card is supplied 
2. Worker process the job                                                            |
3. Worker updates m_pcards with ++ and if it == n_pcards delete the job.             |
   Worker check sem_id and if it zeroes it sleeps                                    |
                           otherwise it grabs new card                               |
                                                                                     | 4. Wait while sem_id==0, jobs are done
 
*/

//This function makes threads pool and inits structures
t_ypsystem *init_ypsystem(unsigned int nthreads,unsigned int n_pcards,unsigned int n_pcards,unsigned int n_hpcards,unsigned int n_pjobs);

//This function sets up a job 
//Note it waits if quenes are full unless O_NONBLOCK is set (FALSE is returned then)
char push_job(t_ypsystem *ypsys,char (*funct)(char (*card)[YPCARD_LEN]) ,char (*card)[YPCARD_LEN],int flag);

//This function supply card to the workers
//Note it waits if quenes are full unless O_NONBLOCK is set (FALSE is returned then)
char push_pcard(t_ypsystem *ypsys,unsigned int job_id,char (*card)[YPCARD_LEN],int flag);


