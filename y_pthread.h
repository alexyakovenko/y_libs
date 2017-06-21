
#define Y_PTHREAD 0x1

#ifndef Y_SYSTEM
#include "../y_libs/y_system.h"
#endif
#include <pthread.h>


//The ypthread parallel model is:
// 0. Create a pool of workers and put them to sleep
// 1. While needed, wake up the pool to execute funct() in parallel and then put the pool back to sleep 
// 2. Finish all threads in the parallel pool

//y_POSIX_semaphore structure
typedef struct {
               pthread_mutex_t  m;    //Synchronizetion mutex
               pthread_cond_t   c;    //The condition
               unsigned int     v;    //Value of the semaphore 
               }t_ypsem;

//ypthread parallel infrastructure
typedef struct ypthread_t{
               t_ypsem                         ypsem_s, ypsem_c;   //y pthread semaphores for server and clients
               unsigned int                            nthreads;   //amount of working threads
               pthread_t                               *threads;   //threads ids (adjusted to main structure)
               pthread_mutex_t                            mutex;   //y_pthread synchronization mutex used for (potentially) slow synchronization operations
               pthread_spinlock_t                      spinlock;   //y_pthread synchronization spinlock used for extremelly quick (one-variable) synchronization operations
               volatile int                                  id;   //problem id to work on
               char  (*funct)(unsigned int ,struct ypthread_t*);   //the function to execute in parallel
               union{ 
                    unsigned int                          reg_u;   //coupled general purposes unsigned integer synchronization register
                             int                          reg_i;   //coupled general purposes signed   integer synchronization register   
                    };
               union{
                    double                                reg_d;   //coupled general purposes double real synchronization register
                    float                                 reg_f;   //coupled general purposes float  real synchronization register
                    };                                    
               void                                       *data;   //general purposes data storage pointer
               }t_ypthread; //it has to be allocated via malloc to accomodate threads massive 


//--------------------------------------------------   E N D   O F   P Y S E M   F U N C T I O N S   --------------------------------------------------//

//This function increase condition variable by value
inline char ypsem_incr(register t_ypsem *ypsem,register unsigned int value,register int broadcast,register unsigned int trigger);
//This function decrease condition variable by value
inline char ypsem_decr(register t_ypsem *ypsem,register unsigned int value,register int  broadcast,register unsigned int trigger);

//These are waiting functions
inline char ypsem_wait_value(register t_ypsem *ypsem,register unsigned int trigger);
//The same as above but with trigger == zero 
inline char ypsem_wait__zero(register t_ypsem *ypsem); 

//This function activates a group of sleepeing threads and waits while they do their job
char master_ypsem_activate_and_wait(unsigned int nthreads,t_ypsem *ypsem_s,t_ypsem *ypsem_c);
//These two functions are the split of previous one (for nonblocking runs)
char master_ypsem_activate(unsigned int nthreads,t_ypsem *ypsem_s);
char master_ypsem_wait(unsigned int nthreads,t_ypsem *ypsem_s,t_ypsem *ypsem_c);

//--------------------------------------------------   E N D   O F   P Y S E M   F U N C T I O N S   --------------------------------------------------//

//This function does nothing, it is just a flag to keep a thread
char ypthread_function_nope(unsigned int worker_id,t_ypthread *yphread);
//This function does nothing just a flag to quit a thread
char ypthread_function_quit(unsigned int worker_id,t_ypthread *ypthread);

//This function free ypthreads
void free_ypthread(t_ypthread *ypthread);
//This function initilaize pthreads
t_ypthread *init_ypthread(unsigned int nthreads);

//This function runs ypthread->nthreads of threads and put them on hold
char launch_ypthreads(t_ypthread *ypthread);

//This function execute (*funct) with pthreads in parallel (blocking the calling thread!) 
inline char exec_client_ypthreads(t_ypthread *ypthread,char (*funct)(unsigned int ,struct ypthread_t*));
//These two functions are a nonblocking split of the previous
inline char exec_client_ypthreads_nonblock(t_ypthread *ypthread,char (*funct)(unsigned int ,struct ypthread_t*));
inline char wait_exec_client_ypthreads(t_ypthread *ypthread);

//This function terminates launched ypthread->nthreads threads
char finish_ypthreads(t_ypthread *ypthread);




