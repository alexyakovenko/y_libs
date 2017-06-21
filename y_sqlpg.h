//This library  contain data for interacting with PostgreSQL

#define Y_SQLPG 0x1


//Include PostgreSQL interface
#include <libpq-fe.h> 
#ifndef Y_MOL
  #include "y_mol.h"
#endif
//Our model of interaction with DB assumes one conection per application; these are the pointers
extern char            PG_buff[0x1000];   //The statically preallocated 4K+1bit memory buffer for (SQL) Postgres commands call without malloc  
extern char           *PG_line;           //The line pointer for (SQL) Postgres commands
extern PGconn         *PG_conn;           //The connection structure for Postgres library
extern PGresult       *PG_result;         //The result structure for Postgres library functions
extern PGcancel       *PG_cancel;         //The 'service' cancel structure
extern int             PQ_socket;         //The connection file descriptor
extern ExecStatusType  PG_ExecStatusType; //The status of query result


//This function establish connection to Postgre database server
char connect_postgre_db_server(const char *name_db,const char *name_appl);
void disconnect_postgree_db_server(void);

//This function sets (non)blocking connection mode with postgre server
char set____blocking_postgre_connection(void);
char set_nonblocking_postgre_connection(void);

//This function converts solid binary object into a string of 2 hexadecimal digits per byte
inline void encode_obj_sql_hexadecimal(register size_t size,register char *str,register void *obj);
inline char decode_obj_sql_hexadecimal(register size_t size,register char *str,register void *obj);


//This function converts str structure into a string of 2 hexadecimal digits per byte
//Note. it doesn't encode molecule name as DB has a separate field for it. Additionally edges are compressed, the atomn id filed is choosen either char if (size_a<=0xFF), short (if size_a<=0xFFF) and as int olny otherwise
inline size_t mono_str_line_len(register t_str *str);
inline char encode_mono_str_sql_hexadecimal(register char *line,register t_str *str); //NB! str here is "structure". The encode: | name_len | size_a | nedges | name | a | edges |
inline size_t line_mono_str_len(register char *line);
inline char decode_mono_str_sql_hexadecimal(register t_str *str,register char *line); //NB! str here is "string"; str should has size_a assigned before the call



