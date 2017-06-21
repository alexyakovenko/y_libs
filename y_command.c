#include "y_command.h"

extern unsigned int ylib_errno;


//This function finds the type of iserver in constant able
const char *get_yserver_type(unsigned int yserver_type)
{
switch (yserver_type)
  {
  case YAPPLICATION_KERNEL    : return "y_kernel";     //00 y_kernel application identifier in Y  --- LOCALY UNIOQUE in Y ---
  case YAPPLICATION_CONSOLE   : return "y_console";    //01 y_console application identifier in Y
  case YAPPLICATION_TBDB      : return "y_tbdb";       //02 y_db topology database identifier in Y  --- GLOBALLY UNIQUE in Y ---
  case YAPPLICATION_SPLITTER  : return "y_splitter";   //04 y_splitter application identifier in Y
  default                     : return "unknown"; 
  }
}

//------------------------------------------------ C O M M A N D    P A R T ---------------------------------------------

//This rather script checks if current command is of this yapplication
inline char is_me(register unsigned int e_id,register unsigned int i_id)
{
return ( (e_id==my_yaddr.e_id)||(i_id==my_yaddr.i_id) );
}

//This function setup command structure and zero QUIT command into waitining mode
char init_ycommands(int size_c)
{
if (!size_c) { ylib_errno=YERROR_INTERNAL_CODE; return FALSE; }
if (!(ycommands=(t_ycommand*)malloc(sizeof(t_ycommand)*size_c))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { size_ycommands=1, _size_ycommands=size_c; }
while (--size_c) { ycommands[size_c].status=YSTATUS_FREE; }
ycommands->status=YSTATUS_RESPAWN; 
ycommands->command=YCOMMAND_QUIT;
ycommands->s_addr.e_id=my_yaddr.e_id, ycommands->s_addr.i_id=my_yaddr.i_id;
ycommands->t_addr.e_id=my_yaddr.e_id, ycommands->t_addr.i_id=my_yaddr.i_id; //Always addressed to itself
ycommands->r_command=0;
ycommands->data=FALSE;
ycommands->size=0x0;
ycommands->attachment=0x0;
return TRUE;
}

//This  function get free command ( a new one if required )
int get_free_ycommand(int *size_c,int *_size_c,t_ycommand **cmd)
{
register int _i;
void *vp;
for (_i=0; _i<(*_size_c); _i++) if ((*cmd)[_i].status==YSTATUS_FREE) { if (_i>=(*size_c)) (*size_c)=_i+1; return _i; }
if (!(vp=(void*)realloc(*cmd,sizeof(t_ycommand)*((*_size_c)+1)))) { ylib_errno=YERROR_MEMORY; return -1; } //FALSE; }
else { (*cmd)=(t_ycommand*)vp, (*cmd)[(*_size_c)].status=YSTATUS_FREE, (*size_c)=++(*_size_c); }
return (*size_c)-1;
}
//This function gets free command in a safe way (with temporary blocked signals)
int get_free_ycommand_safe(int *size_c,int *_size_c,t_ycommand **cmd)
{
register int _i;
sigset_t old_set, new_set;
void *vp;
//Try it in a simple way first
for (_i=0; _i<(*_size_c); _i++) if ((*cmd)[_i].status==YSTATUS_FREE) { if (_i>=(*size_c)) (*size_c)=_i+1; return _i; }
//OK we need more events, allocate them safely
if ( (sigfillset(&new_set)==-1)||(sigprocmask(SIG_BLOCK,&new_set,&old_set)==-1) ) { LABEL_ERROR: ylib_errno=YERROR_EXTERNAL_CODE; return -1; } //FALSE; }
if (!(vp=realloc((*cmd),sizeof(t_ycommand)*((*_size_c)+1)))) 
  {
  if (sigprocmask(SIG_SETMASK,&old_set,NULL)==-1) goto LABEL_ERROR; //Restore signals mask
  else { ylib_errno=YERROR_MEMORY; return -1; }
  }
else { (*cmd)=(t_ycommand*)vp, (*cmd)[(*_size_c)].status=YSTATUS_FREE, (*size_c)=++(*_size_c); }
if (sigprocmask(SIG_SETMASK,&old_set,NULL)==-1) goto LABEL_ERROR; //Restore signals mask
return (*size_c)-1;
}

//This function forms commands in Y
void form_ycommand(unsigned char status,unsigned int s_e_id,unsigned int s_i_id,unsigned int t_e_id,unsigned int t_i_id,unsigned int command,int r_command,int data,size_t size,void *attachment,t_ycommand *cmd)
{
cmd->status=status;
cmd->t_addr.e_id=t_e_id, cmd->t_addr.i_id=t_i_id;
cmd->s_addr.e_id=s_e_id, cmd->s_addr.i_id=s_i_id;
cmd->command=command;
cmd->r_command=r_command;
cmd->data=data;
cmd->size=size;
cmd->attachment=attachment;
}

//This function free command. If the command is for/from transition, the attachment destructor is free() otherwise the destructor is cmd->afree(). 
inline void _delete_ycommand(int *size_c,int _i,t_ycommand *cmd)
{
if ( (cmd[_i].attachment))
  {
  if ( (cmd[_i].t_addr.e_id==cmd[_i].s_addr.e_id)&&(cmd[_i].t_addr.i_id==cmd[_i].s_addr.i_id) ) cmd[_i].afree(cmd[_i].attachment); else free(cmd[_i].attachment);
  cmd[_i].attachment=0x0;
  }
cmd[_i].status=YSTATUS_FREE;
while (cmd[(*size_c)-1].status==YSTATUS_FREE) (*size_c)--;
}
inline void delete_ycommand(int *size_c,int id,t_ycommand *cmd)
{
register int _i;
switch (cmd[id].status)
  {
  case YSTATUS_ACTIVE  : { _delete_ycommand(size_c,id,cmd); break; }
  case YSTATUS_AWAIT   : { //Remove all dependent commands
    LABEL_DELETE_AWAIT_COMMAND: _i=(*size_c); while (--_i) if (cmd[_i].r_command==id) _delete_ycommand(size_c,_i,cmd);  _delete_ycommand(size_c,id,cmd); break; }
  case YSTATUS_BUSY    : { _i=cmd[id].r_command; if (cmd[_i].status==YSTATUS_AWAIT) { id=_i; goto LABEL_DELETE_AWAIT_COMMAND; } else _delete_ycommand(size_c,id,cmd); break; }
  case YSTATUS_DELETE  : { _delete_ycommand(size_c,id,cmd); break; }
  case YSTATUS_DONE    : { _delete_ycommand(size_c,id,cmd); break; }
  case YSTATUS_FREE    : break; //do nothing
  case YSTATUS_KILL    : break; //do nothing
  case YSTATUS_SEND    : { cmd[id].status=YSTATUS_KILL; break; } 
  case YSTATUS_SENT    : { _delete_ycommand(size_c,id,cmd); break; } 
  case YSTATUS_SLEEP   : { _delete_ycommand(size_c,id,cmd); break; } 
  case YSTATUS_RECV    : { cmd[id].status=YSTATUS_KILL; break; } 
  case YSTATUS_RECD    : { _i=cmd[id].r_command; if (cmd[_i].status==YSTATUS_AWAIT) { id=_i; goto LABEL_DELETE_AWAIT_COMMAND; } else _delete_ycommand(size_c,id,cmd); break; }
  case YSTATUS_RESPAWN : { _delete_ycommand(size_c,id,cmd); break; }
  default              : { _delete_ycommand(size_c,id,cmd); break; } //Unspecified elimination
  }
}

//This is built-in NOP gag of Y
char YAPPL_NOP(int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int i,int *size_c,int *_size_c,t_ycommand **cmd,va_list stack)
{
unsigned int j;

if ( ((*cmd)[i].command!=YCOMMAND_EVENT_PRINT)&&((*cmd)[i].command!=YCOMMAND_EVENT_NEW_ISERVER)&&((*cmd)[i].command!=YCOMMAND_EVENT_DEL_ISERVER)&&( ((*cmd)[i].r_command)) ) 
  {
  if (!(j=GET_FREE_YCOMMAND(size_c,_size_c,cmd))) { delete_ycommand(size_c,i,(*cmd)); return FALSE; }
  form_ycommand(YSTATUS_ACTIVE,(*cmd)[i].s_addr.e_id,(*cmd)[i].s_addr.i_id,my_yaddr.e_id,my_yaddr.i_id,YCOMMAND_REPORT,(*cmd)[i].r_command,YERROR_NIMPLEMENTED,0x0,0x0,&(*cmd)[j]);
  }
delete_ycommand(size_c,i,(*cmd));
return TRUE;
}


//------------------------------------------ E N D   O F   C O M M A N D   P A R T --------------------------------------



//------------------------------------------------ T I M E R    P A R T ---------------------------------------------

//SIG_ALRM handler
void sig_alrm(int sig)
{
register int _i;
time_t curr_time, next_time=0;

curr_time=time(NULL);
//Exec events
_i=size_yevents;
while (--_i)
  if (yevents[_i].status==YSTATUS_ACTIVE)
    {
    if (yevents[_i].time<=curr_time)
      {
      if (ycommands[yevents[_i].command_id].status==YSTATUS_SLEEP) ycommands[yevents[_i].command_id].status=YSTATUS_ACTIVE;
      yevents[_i].status=YSTATUS_FREE;
      }
    else if ( (!next_time)||(next_time>yevents[_i].time) ) next_time=yevents[_i].time;
    }
//Restore the time
if (yevents->status==YSTATUS_RESPAWN) alarm(1);
else { if ( (next_time)) alarm(next_time-curr_time); }
}
			  
//This function initializes 
char init_yevents(unsigned int size_e)
{
alarm(0);
if (!(yevents=(t_yevent*)malloc(sizeof(t_yevent)*size_e))) { ylib_errno=YERROR_MEMORY; return FALSE; }
else { size_yevents=0, _size_yevents=size_e; }
while (--size_e) yevents[size_e].status=YSTATUS_FREE;
yevents->status=YSTATUS_SLEEP;
yevents->time=1;
yevents->command_id=-1;
return TRUE;
}
			  
//This function gets free event in a safe way (with temporary blocked signals)
int get_new_event_safe(int *size_e,int *_size_e,t_yevent **events)
{
register int _i;
sigset_t old_set, new_set;
void *vp;

//Try it in a simple way first
for (_i=1; _i<(*_size_e); _i++) if ((*events)[_i].status==YSTATUS_FREE) { if (_i>(*size_e)) (*size_e)=_i; return _i; } 
//OK we need more events, allocate them safely
if ( (sigfillset(&new_set)==-1)||(sigprocmask(SIG_BLOCK,&new_set,&old_set)==-1) ) { LABEL_ERROR: ylib_errno=YERROR_EXTERNAL_CODE; return -1; } // FALSE; }
if (!(vp=realloc((*events),sizeof(t_yevent)*((*_size_e)+1)))) 
  {
  if (sigprocmask(SIG_SETMASK,&old_set,NULL)==-1) goto LABEL_ERROR; //Restore signals mask
  else { ylib_errno=YERROR_MEMORY; return -1; } // NTNF; }
  }
else { (*events)=(t_yevent*)vp, (*events)[(*_size_e)].status=YSTATUS_FREE, (*size_e)=++(*_size_e); }
if (sigprocmask(SIG_SETMASK,&old_set,NULL)==-1) goto LABEL_ERROR; //Restore signals mask
return (*size_e)-1;
}


//------------------------------------------ E N D   O F   T I M E R    P A R T --------------------------------------


