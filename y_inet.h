//This is a header of SMP interprocess communications in y_system
#define Y_INET 0x1

#include <fcntl.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>

#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_COMMAND
#include "y_command.h"
#endif


#define YSYSTEM_INET_NAME "YMAGIC_2004_"
#define MAX_INET_CONNECTION_ATTEMPTS 16
#define MAX_YSHMEM_SIZE 32768  // 32k

#define SVSEM_MODE 0660
#define SVSHM_MODE 0660

//The internal server structure
typedef struct{
              int                yappl;    //Server type 
              t_yaddr            yaddr;    //Servers' address
              unsigned int      status;    //Servers' status
              int                semid;    //System V bi-semaphore id
              int                shmid;    //Shared memory id
              void              *shmem;    //Shared memory pointer
              size_t        shmem_size;    //Shared memory size
              size_t             _size;    //Service field - amount of data transfered over the chanel so far
              unsigned int       _ycmd;    //Service field - id of a command associated with the iserver (recv/send)
              }t_yinet;                    //Structure to work with internal connections

int ykernel_semid;

//This function deletes internal server
char remove_master_islot(t_yinet *yinet);
//This function create internal communication block for clients to inquire connection requests
char open_master_igate(char flag,int i_id,t_yinet *yinet);
//This function connect client to master
char master_connect(unsigned int *i_id,unsigned int *size_i,t_yinet **yinet);

//This function create internal communication block for clients to inquire connection requests
char create_client_bislot(unsigned int yappl,t_yinet *yinet);
//This function connects to y_kernel via inet interface.
//It returns TRUE on success; on failure it returns NTNF if reconnection attempt is feasiabble (later) and FALSE otherwise.
char client_connect(unsigned int yappl,t_yinet *yinet);

//This function process kernel servers. 
char master_process_iserver(char block,unsigned int *i_id,unsigned int *size_i,t_yinet **yinet,int *size_c,int *_size_c,t_ycommand **ycommands);
//This function process kernel servers.
char client_process_iserver(char block,t_yinet *yinet,int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int *size_c,int *_size_c,t_ycommand **ycommands);

//------------------------------------------ S T A R T   O F   Y - L O O P   P A R T --------------------------------------

//This is the main loop of clients in Y
char loop_ycmd(t_yinet *yinet,int *size_c,int *_size_c,t_ycommand **cmd, int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),
                                          char (*YAPPL_EXEC_ACTIVE_CMD) (int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int ,int* ,int* ,t_ycommand** ,char ,va_list ),
                                          char (*YAPPL_EXEC_DONE_CMD)   (int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int ,int* ,int* ,t_ycommand** ,char ,va_list ),
                                          char (*YAPPL_EXEC_RESPAWN_CMD)(int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int ,int* ,int* ,t_ycommand** ,char ,va_list ),char verbose, ... );

//------------------------------------------ E N D   O F   Y - L O O P   P A R T --------------------------------------

