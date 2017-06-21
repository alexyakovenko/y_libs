//This is a module of SMP interprocess communications in y_system
#include "y_inet.h"

//Logic:
// [ 0 ][ 0 ] - bsemaphore free
// [ 1 ][ 1 ] = bsemaphore grabbed
// [ 1 ][ 0 ] - shared for master 
// [ 0 ][ 1 ] - shared for client  
//Grab  shared  bin semaphore
static struct sembuf  bsem_grab[4]        = { 
                                              { 0,  0, IPC_NOWAIT },
                                              { 1,  0, IPC_NOWAIT },
                                              { 0, +1, IPC_NOWAIT },
                                              { 1, +1, IPC_NOWAIT }
                                            }; //[0/0] -> [1,1]
//static struct sembuf  bsem_wgrab[4]       = {
//                                              { 0,  0, 0 },
//                                              { 1,  0, 0 },
//                                              { 0, +1, 0 },
//                                              { 1, +1, 0 }
//                                            }; //[0/0] -> [1,1]
//Get grabbed bin semaphore
static struct sembuf  bsem_master_get[4]  = { 
                                              { 0, -1, IPC_NOWAIT },
                                              { 1,  0, IPC_NOWAIT },
                                              { 0, +1, IPC_NOWAIT },
                                              { 1, +1, IPC_NOWAIT }
                                            }; //[1/0] -> [1/1]
static struct sembuf  bsem_client_get[4]  = {
                                              { 0,  0, IPC_NOWAIT },
                                              { 1, -1, IPC_NOWAIT },
                                              { 0, +1, IPC_NOWAIT },
                                              { 1, +1, IPC_NOWAIT }
                                            }; //[0/1] -> [1,1]
static struct sembuf  bsem_client_wget[4] = {
                                              { 0,  0, 0 },
                                              { 1, -1, 0 },
                                              { 0, +1, 0 },
                                              { 1, +1, 0 }
                                            }; //[0/1] -> [1,1]
//Put grabbed bin semaphore
static struct sembuf  bsem_master_put[4]  = { 
                                              { 0, -1, IPC_NOWAIT },
                                              { 1, -1, IPC_NOWAIT },
                                              { 0,  0, IPC_NOWAIT },
                                              { 1, +1, IPC_NOWAIT }
                                            }; //[1/1] -> [0,1]
static struct sembuf  bsem_client_put[4]  = {
                                              { 0, -1, IPC_NOWAIT },
                                              { 1, -1, IPC_NOWAIT },
                                              { 0, +1, IPC_NOWAIT },
                                              { 1,  0, IPC_NOWAIT }
                                            }; //[1/1] -> [1,0]
//Share grabbed bin semaphore
static struct sembuf  bsem_share[4]       = { 
                                              { 0, -1, IPC_NOWAIT },
                                              { 1, -1, IPC_NOWAIT },
                                              { 0,  0, IPC_NOWAIT },
                                              { 1,  0, IPC_NOWAIT }
                                            }; //[1/1]->[0,0]
//Kernel wait till inet request
static struct sembuf  tsem_wait[1]        = { 
                                              { 2, -1, 0 }
                                            }; //[?/?/?]->[?/?/?-1]
//Info kernel about an inet request
static struct sembuf  tsem_info[1]        = { 
                                              { 2, +1, IPC_NOWAIT }
                                            }; //[?/?/?]->[?/?/?+1]

union semum{
           int val;
           struct semid_ds  *buf;
           unsigned short *array;
           struct seminfo *__buf;  
           };

//The internal connections structure
typedef struct{
              t_yaddr                   yaddr;
              int         yappl, semid, shmid;
              size_t               shmem_size;        //This field used in connection response as DENY_OF_ISERVICE flag
              }t_iconnect;

unsigned short bsem[2], tsem[3]; 
struct semid_ds sem_info;

//This function deletes internal server
char remove_master_islot(t_yinet *yinet)
{
return (char)( (semctl(yinet->semid,0,IPC_RMID)!=-1)&&(shmdt(yinet->shmem)!=-1)&&(shmctl(yinet->shmid,IPC_RMID,0x0)!=-1) );
}

//This function create internal communication block for clients to inquire connection requests
char open_master_igate(char flag,int i_id,t_yinet *yinet)
{
char sflag;
key_t key;
//Stage I. Create igate
if (!yinet) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
else yinet->shmem_size=sizeof(t_iconnect);
if ((key=ftok("y_kernel",Y_MAGIC))==(key_t)-1) { LABEL_EXTERNAL_CODE_ERROR_0: ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
sflag=flag;
LABEL_SEMSET:
if ((yinet->semid=semget(key,0x3,IPC_CREAT|IPC_EXCL|SVSEM_MODE))==-1)
  {
  if ( (errno==EEXIST)&&(sflag)&&((yinet->semid=semget(key,0x3,SVSEM_MODE))!=-1) )
    { 
    semctl(yinet->semid,0,IPC_RMID);
    if ((yinet->shmid=shmget(key,yinet->shmem_size,SVSHM_MODE))!=-1) shmctl(yinet->shmid,IPC_RMID,0x0);
    sflag=FALSE;
    goto LABEL_SEMSET;
    }
  else return FALSE;
  }
sflag=flag;
LABEL_SHMSET:
if ((yinet->shmid=shmget(key,yinet->shmem_size,IPC_CREAT|SVSHM_MODE))==-1)
  {
  if ( (errno==EEXIST)&&(sflag)&&((yinet->semid=shmget(key,yinet->shmem_size,IPC_CREAT|SVSHM_MODE))!=-1) )
    {
    semctl(yinet->semid,0,IPC_RMID);
    sflag=FALSE;
    goto LABEL_SHMSET;    
    }
  LABEL_EXTERNAL_CODE_ERROR_1: semctl(yinet->semid,0,IPC_RMID); goto LABEL_EXTERNAL_CODE_ERROR_0;
  }
if ((yinet->shmem=shmat(yinet->shmid,NULL,0x0))==(void*)-1) { LABEL_EXTERNAL_CODE_ERROR_2: shmctl(yinet->shmid,0x0,IPC_RMID); goto LABEL_EXTERNAL_CODE_ERROR_1; }
//Stage II. Initialize the igate
tsem[0]=tsem[1]=tsem[2]=0;
((t_iconnect*)yinet->shmem)->yappl=YAPPLICATION_KERNEL, ((t_iconnect*)yinet->shmem)->yaddr.e_id=my_yaddr.e_id, ((t_iconnect*)yinet->shmem)->yaddr.i_id=i_id;
if ( (semctl(yinet->semid,0x0,SETALL,tsem)==-1)||(semop(yinet->semid,bsem_grab,0x4)==-1)||(semop(yinet->semid,bsem_master_put,0x4)==-1) ) { shmdt(yinet->shmem); goto LABEL_EXTERNAL_CODE_ERROR_2; }
else yinet->status=YSTATUS_AWAIT;
return TRUE;
}

//This function connect client to master
char master_connect(unsigned int *i_id,unsigned int *size_i,t_yinet **yinet)
{
register unsigned int _i;
void *vp;
//Check for the incoming connection request
if (semop((*yinet)->semid,bsem_master_get,0x4)==-1)
  {//Nop
  if (errno==EAGAIN) { ylib_errno=YERROR_OK; return NTNF; } //It's OK - just nobody willing to connect at the moment
  else { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; } //Something wrong with semaphore...
  }
if ( (((t_iconnect*)(*yinet)->shmem)->yappl==YAPPLICATION_KERNEL)||(((t_iconnect*)(*yinet)->shmem)->shmem_size<sizeof(t_ycommand)) ) 
  {//Denay of service
  LABEL_DENAY_OF_ISERVICE:;
  if ( (semctl(((t_iconnect*)(*yinet)->shmem)->semid,0,IPC_RMID)!=-1)||(shmctl(((t_iconnect*)(*yinet)->shmem)->shmid,IPC_RMID,0x0)!=-1) ) return NTNF;
  ((t_iconnect*)(*yinet)->shmem)->yappl=YAPPLICATION_KERNEL, ((t_iconnect*)(*yinet)->shmem)->yaddr.e_id=my_yaddr.e_id, ((t_iconnect*)(*yinet)->shmem)->yaddr.i_id=*i_id;
  if (semop((*yinet)->semid,bsem_master_put,0x4)==-1) return NTNF;
  else return FALSE;
  }
//We have a client knocking...!
//Realloc iservers
if (!(vp=(void*)realloc((*yinet),sizeof(t_yinet)*(*size_i+1))))
  { ylib_errno=YERROR_MEMORY; goto LABEL_DENAY_OF_ISERVICE; }
else 
  {
  (*yinet)=(t_yinet*)vp;
  (*yinet)[(*size_i)].yappl=((t_iconnect*)(*yinet)->shmem)->yappl;
  (*yinet)[(*size_i)].yaddr.e_id=((t_iconnect*)(*yinet)->shmem)->yaddr.e_id, (*yinet)[(*size_i)].yaddr.i_id=((t_iconnect*)(*yinet)->shmem)->yaddr.i_id; 
  (*yinet)[(*size_i)].semid=(unsigned int)((t_iconnect*)(*yinet)->shmem)->semid;
  (*yinet)[(*size_i)].shmid=(unsigned int)((t_iconnect*)(*yinet)->shmem)->shmid;
  (*yinet)[(*size_i)].shmem_size=(unsigned int)((t_iconnect*)(*yinet)->shmem)->shmem_size;
  if (((*yinet)[(*size_i)].shmem=shmat((*yinet)[(*size_i)].shmid,NULL,0))==(void*)-1) goto LABEL_DENAY_OF_ISERVICE;
  else (*yinet)[(*size_i)].status=YSTATUS_ACTIVE; 
  (*size_i)++;
  LABEL_NEXT_IID: (*i_id)++, _i=(*size_i-1); while (_i--) { if ((*yinet)[_i].yaddr.i_id==(*i_id)) goto LABEL_NEXT_IID; }
  ((t_iconnect*)(*yinet)->shmem)->yappl=YAPPLICATION_KERNEL, ((t_iconnect*)(*yinet)->shmem)->yaddr.e_id=my_yaddr.e_id, ((t_iconnect*)(*yinet)->shmem)->yaddr.i_id=*i_id;
  return (semop((*yinet)->semid,bsem_master_put,0x4)!=-1);
  }
}

//This function create internal communication block for clients to inquire connection requests
char create_client_bislot(unsigned int yappl,t_yinet *yinet)
{
extern unsigned int ylib_errno;
//Stage I. Create private bislot
if (!yinet) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
if ((yinet->semid=semget(IPC_PRIVATE,0x2,IPC_CREAT|IPC_EXCL|SVSEM_MODE))==-1) { LABEL_EXTERNAL_CODE_ERROR_0: ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
if ((yinet->shmid=shmget(IPC_PRIVATE,yinet->shmem_size,IPC_CREAT|SVSHM_MODE))==-1) { LABEL_EXTERNAL_CODE_ERROR_1: semctl(yinet->semid,0,IPC_RMID); goto LABEL_EXTERNAL_CODE_ERROR_0; }
if ((yinet->shmem=shmat(yinet->shmid,NULL,0x0))==(void*)-1) { LABEL_EXTERNAL_CODE_ERROR_2: shmctl(yinet->shmid,0x0,IPC_RMID); goto LABEL_EXTERNAL_CODE_ERROR_1; }
//Stage II. Initialize the semaphores
bsem[0]=bsem[1]=0;
if (semctl(yinet->semid,0x0,SETALL,bsem)==-1) { shmdt(yinet->shmem); goto LABEL_EXTERNAL_CODE_ERROR_2; }
else { yinet->yappl=yappl, yinet->status=YSTATUS_ACTIVE, yinet->yaddr.e_id=my_yaddr.e_id, yinet->yaddr.i_id=0; }
return TRUE;
}

//This function connects to y_kernel via inet interface.
//It returns TRUE on success; on failure it returns NTNF if reconnection attempt is feasiabble (later) and FALSE otherwise.
char client_connect(unsigned int yappl,t_yinet *yinet)
{
register int _i;
key_t key;
int semid, shmid;
void *shmem;
//Stage I. Briefly check data integrity
if ( (!yinet)||(yappl==YAPPLICATION_KERNEL)||(yinet->shmem_size<sizeof(t_ycommand))||(yinet->shmem_size>MAX_YSHMEM_SIZE) ) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
//Stage II. Connect to the gate
if ((key=ftok("y_kernel",Y_MAGIC))==(key_t)-1) goto LABEL_EXTERNAL_CODE_ERROR;
if ((semid=semget(key,0x3,SVSEM_MODE))==-1) goto LABEL_EXTERNAL_CODE_ERROR; //No such slot (yet)
else ykernel_semid=semid;
for (_i=0; _i<MAX_INET_CONNECTION_ATTEMPTS; _i++)
  {
  if (semctl(semid,0x0,IPC_STAT,&sem_info)==-1) goto LABEL_EXTERNAL_CODE_ERROR; 
  if (!sem_info.sem_otime) sleep(1); //wait till initialization
  else goto LABEL_INIT_CLIENT;
  }
LABEL_EXTERNAL_CODE_ERROR: ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; //Too many tries
LABEL_INIT_CLIENT: ;
if ((shmid=shmget(key,sizeof(t_iconnect),SVSHM_MODE))==-1) goto LABEL_EXTERNAL_CODE_ERROR;
if ((shmem=shmat(shmid,NULL,0))==(void*)-1) goto LABEL_EXTERNAL_CODE_ERROR;
//Stage III. Negotiate the request
CLIENT_WGET: if (semop(semid,bsem_client_wget,0x4)==-1)
               {
                    if   (errno==EINTR) goto 	CLIENT_WGET; // signal received
               else if ( (errno==EIDRM)&&(shmdt(shmem)!=-1) ) { ylib_errno=YERROR_EXTERNAL_CODE; return NTNF; } //semaphore is removed by kernel
               else return FALSE;
               }
if (((t_iconnect*)shmem)->yappl!=YAPPLICATION_KERNEL) { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
else { my_yaddr.e_id=((t_iconnect*)shmem)->yaddr.e_id, my_yaddr.i_id=((t_iconnect*)shmem)->yaddr.i_id; }
//Stage IV. Create private islot
if (!(create_client_bislot(yappl,yinet))) 
  {
  LABEL_ERROR_EXTERNAL_CODE_0: ylib_errno=YERROR_EXTERNAL_CODE;
  if ( (shmdt(shmem)!=-1)||(semop(semid,bsem_master_put,0x4)!=-1) ) return NTNF; else return FALSE;
  }
((t_iconnect*)shmem)->yappl=yinet->yappl, ((t_iconnect*)shmem)->semid=yinet->semid, ((t_iconnect*)shmem)->shmid=yinet->shmid, ((t_iconnect*)shmem)->shmem_size=yinet->shmem_size;
if ( (semop(semid,bsem_client_put,0x4)==-1)||(semop(semid,tsem_info,0x1)==-1)||(shmdt(shmem)==-1) )
  {
  ((t_iconnect*)shmem)->yappl=YAPPLICATION_KERNEL;
  semctl(yinet->semid,0,IPC_RMID);
  shmdt(yinet->shmem);
  shmctl(yinet->shmid,IPC_RMID,0x0);
  goto LABEL_ERROR_EXTERNAL_CODE_0;
  }
return TRUE;
}

//This function process kernel servers.
char master_process_iserver(char block,unsigned int *i_id,unsigned int *size_i,t_yinet **yinet,int *size_c,int *_size_c,t_ycommand **ycommands)
{
register int _i, _j;
register char _c, flag=NTNF;
struct sembuf tsem;

//if (*size_i>1) printf("[%1d:%1d] block=%1d, cmd[2].status=%1d\n",semctl((*yinet)[1].semid,0,GETVAL),semctl((*yinet)[1].semid,1,GETVAL),block,(*ycommands)[2].status);

//Stage I. Block on core semaphore
if ((_j=semctl((*yinet)->semid,0x2,GETVAL))==-1) { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
if (_j) 
  { 
  tsem.sem_num=0x2, tsem.sem_op=-_j,  tsem.sem_flg=IPC_NOWAIT; 
  if (semop((*yinet)->semid,&tsem,0x1)==-1) { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
//  printf("_j=%1d\n",_j);
  }
else if ( (block))
       {//There is nothing to do now - block till something happens
       if (semop((*yinet)->semid,tsem_wait,0x1)==-1)
         { if ( (errno!=EAGAIN)&&(errno!=EINTR) ) { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; } } 
       } 
//Stage II. Check the rest of iservers
_i=*size_i;
while (--_i)
       if ((*yinet)[_i].status==YSTATUS_ACTIVE)
         {
         //Stage II.1. Try to send something
         _j=*size_c; 
         while (_j--)
           if ( ((*ycommands)[_j].status==YSTATUS_ACTIVE)&&((*yinet)[_i].yaddr.e_id==(*ycommands)[_j].t_addr.e_id)&&((*yinet)[_i].yaddr.i_id==(*ycommands)[_j].t_addr.i_id) )
             {//Initiating send
             flag=TRUE;
             if (semop((*yinet)[_i].semid,bsem_grab,0x4)==-1)
               {
               if (errno!=EAGAIN)
                 { //Problem with semaphore - kill iserver 
                 LABEL_KILL_ISERVER: flag=TRUE;
                 if ( ((*yinet)[_i].status==YSTATUS_SEND)||((*yinet)[_i].status==YSTATUS_RECV) )
                   { (*ycommands)[(*yinet)[_i]._ycmd].status=YSTATUS_DELETE; delete_ycommand(size_c,(*yinet)[_i]._ycmd,(*ycommands)); }
                 if (!(remove_master_islot(&(*yinet)[_i]))) { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
                 //Delete iserver
                 if ((_j=get_free_ycommand_safe(size_c,_size_c,ycommands))==-1) return FALSE;
                 form_ycommand(YSTATUS_ACTIVE,my_yaddr.e_id,my_yaddr.i_id,(*yinet)[_i].yaddr.e_id,(*yinet)[_i].yaddr.i_id,YCOMMAND_EVENT_DEL_ISERVER,0x0,0x0,0x0,0x0,&(*ycommands)[_j]);
                 if (_i!=--(*size_i)) memcpy(&(*yinet)[_i],&(*yinet)[*size_i],sizeof(t_yinet));
                 }  
               else break; //Hm... hold on - the client overrun me
               }
             else
               {//Send new command
               (*yinet)[_i].status=YSTATUS_SEND, (*yinet)[_i]._ycmd=_j, (*ycommands)[_j].status=YSTATUS_SEND;
               memcpy((*yinet)[_i].shmem,&(*ycommands)[_j],sizeof(t_ycommand));
               if ( ( ((*ycommands)[_j].r_command))&&((*ycommands)[_j].command!=YCOMMAND_REPORT) ) //Update r_command for AWAIT scripts and simple execution mode
                 { ((t_ycommand*)(*yinet)[_i].shmem)->r_command=_j, (*ycommands)[_j].r_command=-(*ycommands)[_j].r_command; } 
               if (!((*ycommands)[_j].size)) (*yinet)[_i]._size=sizeof(t_ycommand);
               else
                 {
                 if (sizeof(t_ycommand)+(*ycommands)[_j].size<=(*yinet)[_i].shmem_size)
                   {//Send entire attachment
                   memcpy((*yinet)[_i].shmem+sizeof(t_ycommand),(*ycommands)[_j].attachment,(*ycommands)[_j].size);
                   (*yinet)[_i]._size=(*ycommands)[_j].size;
                   }
                 else
                   {//send part of attachment
                   memcpy((*yinet)[_i].shmem+sizeof(t_ycommand),(*ycommands)[_j].attachment,(*yinet)[_i].shmem_size-sizeof(t_ycommand));
                   (*yinet)[_i]._size=(*yinet)[_i].shmem_size-sizeof(t_ycommand);
                   }
                 }
               if (semop((*yinet)[_i].semid,bsem_master_put,0x4)==-1) goto LABEL_KILL_ISERVER; 
               }
             goto NEXT_ISERVER;
             }
         //Stage II.2. Try to receive something
//         _j=semctl((*yinet)[_i].semid,0,GETVAL);
//         _j=semctl((*yinet)[_i].semid,1,GETVAL);
         if (semop((*yinet)[_i].semid,bsem_master_get,0x4)==-1)
           { if (errno!=EAGAIN) goto LABEL_KILL_ISERVER; }
         else
           {//Incoming message, begin command receiving
           flag=TRUE;
           if ((_j=get_free_ycommand_safe(size_c,_size_c,ycommands))==-1) return FALSE;
           memcpy(&(*ycommands)[_j],(*yinet)[_i].shmem,sizeof(t_ycommand));
           if (!((*ycommands)[_j].size)) 
             {
             (*ycommands)[_j].attachment=0x0;
             (*ycommands)[_j].status=YSTATUS_ACTIVE;
             *((size_t*)(*yinet)[_i].shmem)=sizeof(t_ycommand);
             }
           else
             {
             if (!((*ycommands)[_j].attachment=malloc((*ycommands)[_j].size))) { ylib_errno=YERROR_MEMORY; return FALSE; }
             if (sizeof(t_ycommand)+(*ycommands)[_j].size<=(*yinet)[_i].shmem_size)
               {//Download entire attachment at once 
               memcpy((*ycommands)[_j].attachment,(*yinet)[_i].shmem+sizeof(t_ycommand),(*ycommands)[_j].size);
               (*ycommands)[_j].status=YSTATUS_ACTIVE;
               *((size_t*)(*yinet)[_i].shmem)=(*ycommands)[_j].size;
               }
             else
               {//Download a part of attachment
               (*yinet)[_i].status=YSTATUS_RECV, (*yinet)[_i]._ycmd=_j; (*yinet)[_i]._size=(*yinet)[_i].shmem_size-sizeof(t_ycommand);
               (*ycommands)[_j].status=YSTATUS_RECV;
               memcpy((*ycommands)[_j].attachment,(*yinet)[_i].shmem+sizeof(t_ycommand),(*yinet)[_i]._size);
               *((size_t*)(*yinet)[_i].shmem)=(*yinet)[_i]._size;
//               printf(">%s\n",(char*)(*ycommands)[_j].attachment);
               }
             }
           if (semop((*yinet)[_i].semid,bsem_master_put,0x4)==-1) goto LABEL_KILL_ISERVER; 
           }  
         NEXT_ISERVER: ;
         }
  else if ((*yinet)[_i].status==YSTATUS_SEND)
         {//Continue/finalize command sending
         if (semop((*yinet)[_i].semid,bsem_master_get,0x4)==-1)
           { if (errno!=EAGAIN) goto LABEL_KILL_ISERVER; }
         else
           {
           flag=TRUE; 
           if (*((size_t*)(*yinet)[_i].shmem)!=(*yinet)[_i]._size) goto LABEL_KILL_ISERVER; //Transaction missmatch
           else _j=(*yinet)[_i]._ycmd;
           if ( ((*ycommands)[_j].status!=YSTATUS_SEND)&&((*ycommands)[_j].status!=YSTATUS_KILL) ) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
           if ( (!(*ycommands)[_j].size)||((*yinet)[_i]._size==(*ycommands)[_j].size) )
             {//Command is sent
             if ((*ycommands)[_j].status==YSTATUS_KILL) (*ycommands)[_j].status=YSTATUS_DELETE;  
             else (*ycommands)[_j].status=YSTATUS_SENT;
             (*yinet)[_i].status=YSTATUS_ACTIVE;
             if (semop((*yinet)[_i].semid,bsem_share,0x4)==-1) goto LABEL_KILL_ISERVER;
             }
           else
             {//Coninue trasmission
             if ((*ycommands)[_j].size<=(*yinet)[_i]._size+(*yinet)[_i].shmem_size)
               {//Send the last piece
               memcpy((*yinet)[_i].shmem,(*ycommands)[_j].attachment+(*yinet)[_i]._size,(*ycommands)[_j].size-(*yinet)[_i]._size);
               (*yinet)[_i]._size=(*ycommands)[_j].size;
               }
             else
               {//Send yet another piece
               memcpy((*yinet)[_i].shmem,(*ycommands)[_j].attachment+(*yinet)[_i]._size,(*yinet)[_i].shmem_size);
               (*yinet)[_i]._size+=(*yinet)[_i].shmem_size;
               }
             if (semop((*yinet)[_i].semid,bsem_master_put,0x4)==-1) goto LABEL_KILL_ISERVER;
             }
           }
         }
  else if ((*yinet)[_i].status==YSTATUS_RECV)
         {//Continue/finalize command receiving
         if (semop((*yinet)[_i].semid,bsem_master_get,0x4)==-1)
           { if (errno!=EAGAIN) goto LABEL_KILL_ISERVER; }
         else 
           {
           flag=TRUE;
           _j=(*yinet)[_i]._ycmd;
           if ( ((*ycommands)[_j].status!=YSTATUS_RECV)&&(((*ycommands)[_j].status!=YSTATUS_KILL)) ) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
           //Coninue trasmission
           if ((*ycommands)[_j].size<=(*yinet)[_i]._size+(*yinet)[_i].shmem_size)
             {//Receive the last piece
             memcpy((*ycommands)[_j].attachment+(*yinet)[_i]._size,(*yinet)[_i].shmem,(*ycommands)[_j].size-(*yinet)[_i]._size);
             *((size_t*)(*yinet)[_i].shmem)=(*ycommands)[_j].size;
             if ((*ycommands)[_j].status==YSTATUS_KILL) (*ycommands)[_j].status=YSTATUS_DELETE;  
             else (*ycommands)[_j].status=YSTATUS_ACTIVE;
             (*yinet)[_i].status=YSTATUS_ACTIVE;
             }
           else
             {//Receive yet another piece
             memcpy((*ycommands)[_j].attachment+(*yinet)[_i]._size,(*yinet)[_i].shmem,(*yinet)[_i].shmem_size);
             (*yinet)[_i]._size+=(*yinet)[_i].shmem_size;
             *((size_t*)(*yinet)[_i].shmem)=(*yinet)[_i]._size;
             }
//           printf(">%s\n",(char*)(*ycommands)[_j].attachment);
           if (semop((*yinet)[_i].semid,bsem_master_put,0x4)==-1) { delete_ycommand(size_c,(*yinet)[_i]._ycmd,(*ycommands)); goto LABEL_KILL_ISERVER; }
           }
         }
  else goto LABEL_KILL_ISERVER;
//Stage III. Check igate and connect a new server
     if (!(_c=master_connect(i_id,size_i,yinet)))
       {//For some reason i-gate is not working. Try to restart it.
       if ( (!(remove_master_islot(*yinet)))||(!(open_master_igate(TRUE,*i_id,*yinet))) ) return FALSE;
       }
else if (_c==TRUE)
       {//New server is connected
       flag=TRUE;
       if ((_j=get_free_ycommand_safe(size_c,_size_c,ycommands))==-1) return FALSE;
       else form_ycommand(YSTATUS_ACTIVE,(*yinet)[(*size_i)-1].yaddr.e_id,(*yinet)[(*size_i)-1].yaddr.i_id,my_yaddr.e_id,my_yaddr.i_id,YCOMMAND_EVENT_NEW_ISERVER,0x0,(*yinet)[(*size_i)-1].yappl,0x0,0x0,&(*ycommands)[_j]);
       }
//Exit 
return flag;
}

//This function process kernel servers. It DOESN'T block on semaphore.
char client_process_iserver(char block,t_yinet *yinet,int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int *size_c,int *_size_c,t_ycommand **ycommands)
{
register int _j;
register char flag=NTNF;
//Stage I. Check the iserver

//printf("[%1d:%1d] block=%1d, cmd[1].status=%1d\n",semctl(yinet->semid,0,GETVAL),semctl(yinet->semid,1,GETVAL),block,(*ycommands)[1].status);

     if (yinet->status==YSTATUS_ACTIVE)
       {
       //Stage I.1. Try to send something
       _j=*size_c; 
       while (_j--)
         if ( ((*ycommands)[_j].status==YSTATUS_ACTIVE)&&( (my_yaddr.e_id!=(*ycommands)[_j].t_addr.e_id)||(my_yaddr.i_id!=(*ycommands)[_j].t_addr.i_id) ) )
           {//Initiating send
           flag=TRUE;
           if (semop(yinet->semid,bsem_grab,0x4)==-1)
             {
             if (errno!=EAGAIN) { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#1\n"); return FALSE; } //Problem with semaphore - self-kill 
             else break; //Hm... hold on - the server overrun me
             }
           else
             {//Send new command
             yinet->status=YSTATUS_SEND, yinet->_ycmd=_j, (*ycommands)[_j].status=YSTATUS_SEND;
             memcpy(yinet->shmem,&(*ycommands)[_j],sizeof(t_ycommand));
             if ( ( ((*ycommands)[_j].r_command))&&((*ycommands)[_j].command!=YCOMMAND_REPORT) ) //Update r_command for AWAIT scripts or simple execution mode
               { ((t_ycommand*)yinet->shmem)->r_command=_j, (*ycommands)[_j].r_command=-(*ycommands)[_j].r_command; } 
             if (!((*ycommands)[_j].size)) yinet->_size=sizeof(t_ycommand);
             else
               {
               if (sizeof(t_ycommand)+(*ycommands)[_j].size<=yinet->shmem_size)
                 {//Send entire attachment
                 memcpy(yinet->shmem+sizeof(t_ycommand),(*ycommands)[_j].attachment,(*ycommands)[_j].size);
                 yinet->_size=(*ycommands)[_j].size;
                 }
               else
                 {//send part of attachment
                 memcpy(yinet->shmem+sizeof(t_ycommand),(*ycommands)[_j].attachment,yinet->shmem_size-sizeof(t_ycommand));
                 yinet->_size=yinet->shmem_size-sizeof(t_ycommand);
                 }
               }
//             printf(">%s\n",(char*)(*ycommands)[_j].attachment);
             if ( (semop(yinet->semid,bsem_client_put,0x4)==-1)||(semop(ykernel_semid,tsem_info,0x1)==-1) ) { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#2\n"); return FALSE; }
             }
           goto NEXT_ISERVER;
           }
       //Stage I.2. Try to receive something 
       if (semop(yinet->semid,( (block)) ? bsem_client_wget : bsem_client_get, 0x4)==-1)
         { if ( (errno!=EAGAIN)&&(errno!=EINTR) )  { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#3\n"); return FALSE; } else return flag; } //!!!To avoid races
       else
         {//Incoming message, begin command receiving
         flag=TRUE;
         if ((_j=GET_FREE_YCOMMAND(size_c,_size_c,ycommands))==-1) return FALSE;
         memcpy(&(*ycommands)[_j],yinet->shmem,sizeof(t_ycommand));
         if ( ((*ycommands)[_j].t_addr.e_id!=my_yaddr.e_id)||((*ycommands)[_j].t_addr.i_id!=my_yaddr.i_id) )
           { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#4\n"); return FALSE; } //Wrong delivery address
         if (!((*ycommands)[_j].size)) 
           {
           (*ycommands)[_j].attachment=0x0;
           (*ycommands)[_j].status=YSTATUS_ACTIVE;
           *((size_t*)yinet->shmem)=sizeof(t_ycommand);
           }
         else
           {
           if (!((*ycommands)[_j].attachment=malloc((*ycommands)[_j].size))) { ylib_errno=YERROR_MEMORY; printf("EERE#5\n"); return FALSE; }
           if (sizeof(t_ycommand)+(*ycommands)[_j].size<=yinet->shmem_size)
             {//Download entire attachment at once 
             memcpy((*ycommands)[_j].attachment,yinet->shmem+sizeof(t_ycommand),(*ycommands)[_j].size);
             (*ycommands)[_j].status=YSTATUS_ACTIVE;
             *((size_t*)yinet->shmem)=(*ycommands)[_j].size;
             }
           else
             {//Download a part of attachment
             yinet->status=YSTATUS_RECV, yinet->_ycmd=_j; yinet->_size=yinet->shmem_size-sizeof(t_ycommand);
             memcpy((*ycommands)[_j].attachment,yinet->shmem+sizeof(t_ycommand),yinet->_size);
             (*ycommands)[_j].status=YSTATUS_RECV;
             *((size_t*)yinet->shmem)=yinet->_size;
             }
           }
         if ( (semop(yinet->semid,bsem_client_put,0x4)==-1)||(semop(ykernel_semid,tsem_info,0x1)==-1) ) { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#6\n"); return FALSE; } 
         }  
       NEXT_ISERVER: ;
       }
else if (yinet->status==YSTATUS_SEND)
       {//Continue/finalize command sending
       if (semop(yinet->semid,( (block)) ? bsem_client_wget : bsem_client_get, 0x4)==-1)
         { if ( (errno!=EAGAIN)&&(errno!=EINTR) )  { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#7\n"); return FALSE; } else return NTNF; }
       else
         {
         flag=TRUE;
         if (*((size_t*)yinet->shmem)!=yinet->_size) { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#8\n"); return FALSE; } //Transaction missmatch
         else _j=yinet->_ycmd;
         if ( ((*ycommands)[_j].status!=YSTATUS_SEND)&&((*ycommands)[_j].status!=YSTATUS_KILL) ) { ylib_errno=YERROR_INTERNAL_CODE; printf("EERE#9\n"); return FALSE; } //Command mismatch
         if ( (!(*ycommands)[_j].size)||(yinet->_size==(*ycommands)[_j].size) )
           {//Command sent
           if ((*ycommands)[_j].status==YSTATUS_KILL) (*ycommands)[_j].status=YSTATUS_DELETE;
           else (*ycommands)[_j].status=YSTATUS_SENT;
           yinet->status=YSTATUS_ACTIVE;
           if (semop(yinet->semid,bsem_share,0x4)==-1) { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#10\n"); return FALSE; }
           }
         else
           {//Coninue trasmission
           if ((*ycommands)[_j].size<=yinet->_size+yinet->shmem_size)
             {//Send the last piece
             memcpy(yinet->shmem,(*ycommands)[_j].attachment+yinet->_size,(*ycommands)[_j].size-yinet->_size);
             yinet->_size=(*ycommands)[_j].size;
             }
           else
             {//Send yet another piece
             memcpy(yinet->shmem,(*ycommands)[_j].attachment+yinet->_size,yinet->shmem_size);
             yinet->_size+=yinet->shmem_size;
             }
//           printf(">%s\n",(char*)(*ycommands)[_j].attachment);
           if ( (semop(yinet->semid,bsem_client_put,0x4)==-1)||(semop(ykernel_semid,tsem_info,0x1)==-1) ) { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#11\n"); return FALSE; }
           }
         }
       }
else if (yinet->status==YSTATUS_RECV)
       {//Continue/finalize command receiving
       if (semop(yinet->semid,( (block)) ? bsem_client_wget : bsem_client_get, 0x4)==-1)
         { if ( (errno!=EAGAIN)&&(errno!=EINTR) )  { ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#12\n"); return FALSE; } else return NTNF; }
       else 
         {
         flag=TRUE;
         _j=yinet->_ycmd;
         if ( ((*ycommands)[_j].status!=YSTATUS_RECV)&&((*ycommands)[_j].status!=YSTATUS_KILL) ) { ylib_errno=YERROR_INTERNAL_CODE; printf("EERE#13\n"); return FALSE; }
         //Coninue trasmission
         if ((*ycommands)[_j].size<=yinet->_size+yinet->shmem_size)
           {//Receive the last piece
           memcpy((*ycommands)[_j].attachment+yinet->_size,yinet->shmem,(*ycommands)[_j].size-yinet->_size);
           *((size_t*)yinet->shmem)=(*ycommands)[_j].size;
           if ((*ycommands)[_j].status==YSTATUS_KILL) (*ycommands)[_j].status=YSTATUS_DELETE;
           else (*ycommands)[_j].status=YSTATUS_ACTIVE;
           yinet->status=YSTATUS_ACTIVE;
           }
         else
           {//Receive yet another piece
           memcpy((*ycommands)[_j].attachment+yinet->_size,yinet->shmem,yinet->shmem_size);
           yinet->_size+=yinet->shmem_size;
           *((size_t*)yinet->shmem)=yinet->_size;
           }
         if ( (semop(yinet->semid,bsem_client_put,0x4)==-1)||(semop(ykernel_semid,tsem_info,0x1)==-1) )
           { delete_ycommand(size_c,yinet->_ycmd,(*ycommands)); ylib_errno=YERROR_EXTERNAL_CODE; printf("EERE#14\n"); return FALSE; }
         }
       }
else   { ylib_errno=YERROR_INTERNAL_CODE; printf("EERE#15\n"); return FALSE; }
return flag;
}

//------------------------------------------ S T A R T   O F   Y - L O O P   P A R T --------------------------------------

//This is the main loop of clients in Y (except blocking ones like yconsole)
char loop_ycmd(t_yinet *yinet,int *size_c,int *_size_c,t_ycommand **cmd, int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),
                                          char (*YAPPL_EXEC_ACTIVE_CMD) (int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int ,int* ,int* ,t_ycommand** ,char ,va_list ),
                                          char (*YAPPL_EXEC_DONE_CMD)   (int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int ,int* ,int* ,t_ycommand** ,char ,va_list ),
                                          char (*YAPPL_EXEC_RESPAWN_CMD)(int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int ,int* ,int* ,t_ycommand** ,char ,va_list ),char verbose, ... )
{
register int _i, _j;
register char cflag, sflag ,_c;
va_list stack;

va_start(stack,verbose);
cflag=TRUE;
do{
  //Process commands
  sflag=FALSE; 
  _i=*size_c; 
  while (_i--)
    switch ((*cmd)[_i].status)
      {
      case YSTATUS_FREE    : continue; //The command slot is empty at he moment
      case YSTATUS_ACTIVE  : { //The command is quened to be executed
        if ( ((*cmd)[_i].t_addr.e_id==my_yaddr.e_id)&&((*cmd)[_i].t_addr.i_id==my_yaddr.i_id) ) 
          {//Execute ACTIVE command
          switch ((*cmd)[_i].command)
            {
            case YCOMMAND_QUIT   : { va_end(stack); return TRUE; } //Quit command
            case YCOMMAND_REPORT : { //Report command
              if ( ( (_j=(*cmd)[_i].r_command))&&(_j<*size_c)&&((*cmd)[_j].t_addr.e_id==(*cmd)[_i].s_addr.e_id)&&((*cmd)[_j].t_addr.i_id==(*cmd)[_i].s_addr.i_id) )
                {
                if ((*cmd)[_j].status==YSTATUS_BUSY)
                  {
                  (*cmd)[_j].data=(*cmd)[_i].data;
                  if ((*cmd)[_j].size) free((*cmd)[_j].attachment);
                  (*cmd)[_j].size=(*cmd)[_i].size,             (*cmd)[_i].size=0x0;
                  (*cmd)[_j].attachment=(*cmd)[_i].attachment, (*cmd)[_i].attachment=0x0;
                  if ((*cmd)[_j].r_command<0) (*cmd)[_j].status=YSTATUS_DONE;
                  else //Handle AWAIT command
                    {
                    (*cmd)[_j].status=YSTATUS_RECD;
                    if ( ((*cmd)[(*cmd)[_j].r_command].status==YSTATUS_AWAIT)&&((*cmd)[(*cmd)[_j].r_command].data>0) )
                      { if (!--(*cmd)[(*cmd)[_j].r_command].data) (*cmd)[(*cmd)[_j].r_command].status=YSTATUS_DONE; } //update the branching command
                    else delete_ycommand(size_c,_j,(*cmd)); //a failure in scripting instructions 
                    }
                  }
                else if ((*cmd)[_j].status==YSTATUS_DELETE) { (*cmd)[_j].status=YSTATUS_DONE; delete_ycommand(size_c,_j,(*cmd)); }
                }
              delete_ycommand(size_c,_i,(*cmd));
              break; }
            default              : {//Exec the command
                                   if (!(YAPPL_EXEC_ACTIVE_CMD(GET_FREE_YCOMMAND,_i,size_c,_size_c,cmd,verbose,stack))) 
                                     { if (verbose) printf("ERROR occured during execution of ACTIVE command (ylib_errno=%s).\n",get_yerrno(ylib_errno)); va_end(stack); return FALSE; }
                                   } //Call the handler
            }
          sflag=TRUE;
          } 
        break; }      
      case YSTATUS_AWAIT   : continue; //The command is awaiting for the rest to arrive  
      case YSTATUS_BUSY    : continue; //The command is awaiting responses from other Y-components         
      case YSTATUS_DELETE  : { //Delete the command
        delete_ycommand(size_c,_i,(*cmd));
        break; } 
      case YSTATUS_DONE    : { //Done means the command that receivd everything it needs from outside and now ready to run (or fail due to failure reports)
        if (!(YAPPL_EXEC_DONE_CMD(GET_FREE_YCOMMAND,_i,size_c,_size_c,cmd,verbose,stack))) 
          { if (verbose) printf("ERROR occured during execution of DONE command (ylib_errno=%s).\n",get_yerrno(ylib_errno)); va_end(stack); return FALSE; }
        sflag=TRUE;    
        break; }
      case YSTATUS_SEND    : { sflag=TRUE; continue; } //The command is being sent
      case YSTATUS_SENT    : { //The sending is completed
        if ( ((*cmd)[_i].s_addr.e_id==my_yaddr.e_id)&&((*cmd)[_i].s_addr.i_id==my_yaddr.i_id)&&((*cmd)[_i].command!=YCOMMAND_REPORT)&&
             ( ((*cmd)[_i].r_command))&&((*cmd)[_i].r_command<*size_c)&&((*cmd)[(*cmd)[_i].r_command].status=YSTATUS_AWAIT) ) (*cmd)[_i].status=YSTATUS_BUSY;
        else delete_ycommand(size_c,_i,(*cmd));
        break; }
      case YSTATUS_RECV    : { sflag=TRUE; continue; } //The command is being received
      case YSTATUS_RECD    : continue; //Execution of this command is finished
      case YSTATUS_SLEEP   : continue; //The command is received and waiting for processing
      case YSTATUS_RESPAWN : { //The command is permanently active through the application life-time
        if ( ((*cmd)[_i].t_addr.e_id==my_yaddr.e_id)&&((*cmd)[_i].t_addr.i_id==my_yaddr.i_id) ) 
          {
          if (!(_c=YAPPL_EXEC_RESPAWN_CMD(GET_FREE_YCOMMAND,_i,size_c,_size_c,cmd,verbose,stack))) 
            { if (verbose) printf("ERROR occured during execution of RESPAWN command (ylib_errno=%s).\n",get_yerrno(ylib_errno)); va_end(stack); return FALSE; }
          else if (_c!=NTNF) sflag=TRUE; 
          }
        else { if (verbose) fprintf(stderr,"The RESPAWN command to some OTHER server is encountered.\n"); delete_ycommand(size_c,_i,(*cmd)); } //Should never get here
        break; }                
      default              : { //The unknown status of a command
        fprintf(stderr,"The impossible situation occured (unknown command status). Quiting, sorry.\n");
        va_end(stack); return FALSE; }
      }
  //Process i-server (either bloking or not depending on flag)
       if (!(_c=client_process_iserver(!cflag,yinet,GET_FREE_YCOMMAND,size_c,_size_c,cmd)))
         {//Abnormal termination due to IPC failure
         if (verbose) fprintf(stderr,"\n%s is completed abnormally due to IPC faulire\n",Y_APPLICATION);
         va_end(stack); return FALSE;
         }
  else if (_c==TRUE) sflag=TRUE;
  cflag=sflag;
  }while ((*cmd)->status==YSTATUS_RESPAWN);

va_end(stack); //free stack
return TRUE;
}

//------------------------------------------ E N D   O F   Y - L O O P   P A R T --------------------------------------

