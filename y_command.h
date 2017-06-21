//ver 0.0.4
#define Y_COMMAND 0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif

#include <signal.h>
#include <time.h>


/*
Instructions in brief.
0. Zero command can't be changed or 'awaked' ever 
1. Each application has zero command YCOMMAND_QUIT working in a YSTATUS_RESPAWN mode throught entire application life
2. The respawn commands has FALSE on ERRORS, NTNF on SKIP and something other (probably TRUE) otherwise
3. Each non-local command is removed with aid of free() while every local command is removed with aid of cmd->afree()
4. The negative r_command means that r_command will be set to current command only at destination server but the transferred command itself will keep positive r_command at home side and undergo YSTATUS_RECD on the return (used for AWAITs). The positive r_command means it is set to negative value at home side after the transfer and undergo YSTATUS_DONE on the return (simple execution).
*/

//The yserver parameters
typedef struct{ unsigned int e_id, i_id; }t_yaddr; //external/internal address of application in y_system (assigned by y_kernel)

t_yaddr my_yaddr;

//                 ,->(X)    ,->(X)    ,->(X)      ,->(X)    ,->(X)
//Sending process: ACTIVE -> SEND   -> SENT     -> BUSY   -> RECD
//                           `-> KILL -> DELETE -> (X) 

//                               ,->(X) 
//Receiving process: RECV     -> ACTIVE
//                   `-> KILL -> DELETE -> (X)        
  
enum {
     YSTATUS_ACTIVE,                //00 command is pending for execution
     YSTATUS_AWAIT,                 //01 command awaiting for BUSY chain
     YSTATUS_BUSY,                  //02 command is being executing now (await report)
     YSTATUS_DELETE,                //03 the command is marked for normal deletion
     YSTATUS_DONE,                  //04 command execution is finished
     YSTATUS_FREE,                  //05 this command slot is empty
     YSTATUS_KILL,                  //06 this commands marks a sending command for deletion  
     YSTATUS_RECV,                  //07 receiving command
     YSTATUS_RECD,                  //08 satisfied BUSY command
     YSTATUS_RESPAWN,               //09 the respawning command/even    
     YSTATUS_SEND,                  //10 sending command
     YSTATUS_SENT,                  //11 the command have been sent
     YSTATUS_SLEEP                  //12 the command is waitining to be waked up
     };

//Normally all commands sends report to the command source
enum {
     // U N I V E R S A L
     YCOMMAND_REPORT,               //00 this command is a UNIVERSAL report command in Y (error status, including YERR_OK, if given in .data field)
     YCOMMAND_EVENT_PRINT,          //01 this event is printing message if verbose is ON at destination server
     YCOMMAND_EVENT_NEW_ISERVER,    //02 this event inform Y about new internal server; field DATA carry the appl type of the new server. 
     YCOMMAND_EVENT_DEL_ISERVER,    //03 this event is about deleted iserver in Y
     // Y _ K E R N E L
     YKERNEL_FIND_SERVER_GLOABL,    //16 this command is a Y_KERNEL (GLOBAL) REQUEST for a list of y_addrs of this kind of servers (server type is given in .data field)
     YKERNEL_FIND_SERVER_LOCAL,     //17 this command is a Y_KERNEL (LOCAL) REQUEST for a list of y_addrs of this kind of servers (server type is given in .data field)
     // Y _ C O N S O L E
     // Y _ T B D B
     YTBDB_ADD_PROBE_STR,           //33 this command adds new probe str from ystr file (the file name is in attachment)
//     YTBDB_UPDATE_PROBE_MRR,        //34 this command adds probe Mono-root record from *.out FFly file (the filename is in attachemnt, probe id is in data)
     YTBDB_PARAMETERIZE_MRR,        //35 this command adds new mrr_str from attachment or return its parameters if the structure already exists
     YTBDB_PARAMETERIZE_BRR,        //36 this command adds new brr_str from attachment or return its parameters if the structure already exists
     // Y _ E X E C    
     // Y _ S P L I T T E R
     YSPLITTER_PARAMETERIZE,        // this command is used to parameterize molecules in Y
     // G E N E R I C    T E R M I N A L    C O M M A N D
     YCOMMAND_QUIT                  //the very last; this command is a UNIVERSAL command that (un-)quiting application in Y. It's typical respawning zero command.
     };

//Types of YAPPLICATIONS
enum {
     YAPPLICATION_KERNEL,           //00 y_kernel application identifier in Y  --- LOCALY UNIOQUE in Y ---
     YAPPLICATION_CONSOLE,          //01 y_console application identifier in Y
     YAPPLICATION_TBDB,             //02 y_db topology database identifier in Y  --- GLOBALLY UNIQUE in Y ---
     YAPPLICATION_SPLITTER          //04 y_splitter application identifier in Y
     };

//The command structure
typedef struct{
              unsigned char         status;        //Command status
              t_yaddr       s_addr, t_addr;        //Source and target yserver address
              unsigned int         command;        //Command itself
              int                r_command;        //Command to be reported [r_command<0 means use this id during the first submition]
              int                     data;        //Supplementary data of the command (ylib_errno code for reports)
              union{
                   size_t             size;        //Lenght of attached solid data (default meaning)
                   void    (*afree)(void*);        //If the command is for/from transition, the attachment destructor is free() otherwise the destructor is cmd->afree().
                   };
              void             *attachment;        //Attached data
              }t_ycommand;


//Global variables for access from signal handlers 
int size_ycommands, _size_ycommands;
t_ycommand *ycommands; //NOTE zero command is reserved


//This function finds the type of iserver in constant table
const char *get_yserver_type(unsigned int yserver_type);

//------------------------------------------------ C O M M A N D    P A R T ---------------------------------------------

//This function checks if current address is of this yapplication
inline char is_me(register unsigned int e_id,register unsigned int i_id);

//This function setup command structure and zero QUIT command into waitining mode
char init_ycommands(int size_c);

//This  function get free command ( a new one if required )
int get_free_ycommand(int *size_c,int *_size_c,t_ycommand **ycommands);

//This  function get free command ( a new one if required ) with blocking of signal handlers
int get_free_ycommand_safe(int *size_c,int *_size_c,t_ycommand **ycommands);

//This function forms commands in Y
void form_ycommand(unsigned char status,unsigned int s_e_id,unsigned int s_i_id,unsigned int t_e_id,unsigned int t_i_id,unsigned int command,int r_command,int data,size_t size,void *attachment,t_ycommand *cmd);

//This function free command. If the command is for/from transition, the attachment destructor is free() otherwise the destructor is cmd->afree().
void delete_ycommand(int *size_c,int i,t_ycommand *ycommand);

//This is the NOP gag of Y
char YAPPL_NOP(int (*GET_FREE_YCOMMAND)(int* ,int* ,t_ycommand** ),int i,int *size_c,int *_size_c,t_ycommand **cmd,va_list stack);

//------------------------------------------ E N D   O F   C O M M A N D   P A R T --------------------------------------


//------------------------------------------------ T I M E R    P A R T ---------------------------------------------


//The timer structure
typedef struct{
              char              status;
              time_t              time;
              unsigned int  command_id;
              }t_yevent;


//Global variables for access from signal handlers 
int size_yevents, _size_yevents;
t_yevent *yevents; //NOTE zero event is reserved


//SIG_ALRM handler
void sig_alrm(int sig);
			  
//This function initializes 
char init_yevents(unsigned int size_e);
			  
//This function gets free event in a safe way (with temporary blocked signals)
int get_new_event_safe(int *size_e,int *_size_e,t_yevent **events);

//------------------------------------------ E N D   O F   T I M E R    P A R T --------------------------------------



