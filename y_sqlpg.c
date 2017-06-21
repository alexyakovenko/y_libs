//This library  contain data for interacting with PostgreSQL

#include "y_sqlpg.h"

extern unsigned int ylib_errno; 


char            PG_buff[0x1000];   //The statically preallocated 4K+1bit memory buffer for (SQL) Postgres commands call without malloc  
char           *PG_line=0x0;       //The line pointer for (SQL) Postgres commands
PGconn         *PG_conn=0x0;       //The connection structure for Postgres library
PGresult       *PG_result=0x0;     //The result structure for Postgres library functions
PGcancel       *PG_cancel=0x0;     //The 'service' cancel structure
int             PQ_socket=-1;      //The connection file descriptor
ExecStatusType  PG_ExecStatusType; //The status of query result


//This function establish connection to Postgre database server
char connect_postgre_db_server(const char *name_db,const char *name_appl)
{
const char const * DB_NAME_LEXEM="dbname = ";
const char const * APPL_NAME_LEXEM="application_name = ";
register unsigned int len_db, len_appl;
len_db=strlen(name_db), len_appl=strlen(name_appl);
if (1+len_db+strlen(DB_NAME_LEXEM)+len_appl+strlen(APPL_NAME_LEXEM)<0xFFF)
  {//Static memory buffer
  memcpy((void*) PG_buff,(const void *)DB_NAME_LEXEM,sizeof(char)*strlen(DB_NAME_LEXEM)), memcpy((void*)&PG_buff[strlen(DB_NAME_LEXEM)],(const void *)name_db,sizeof(char)*len_db), PG_buff[strlen(DB_NAME_LEXEM)+len_db]=' ';
  memcpy(&PG_buff[1+strlen(DB_NAME_LEXEM)+len_db],APPL_NAME_LEXEM,sizeof(char)*strlen(APPL_NAME_LEXEM)), memcpy(&PG_buff[1+strlen(DB_NAME_LEXEM)+len_db+strlen(APPL_NAME_LEXEM)],name_appl,sizeof(char)*len_appl), PG_buff[1+strlen(DB_NAME_LEXEM)+len_db+strlen(APPL_NAME_LEXEM)+len_appl]='\0';
  PG_conn=PQconnectdb(PG_buff); //Make a connection to the database
  }
else
  {//Memory allocation routine
  if (!(PG_line=(char*)malloc(sizeof(char)*(1+strlen(DB_NAME_LEXEM)+len_db+strlen(APPL_NAME_LEXEM)+len_appl+1)))) { ylib_errno=YERROR_MEMORY; return FALSE; }
  memcpy( PG_line,DB_NAME_LEXEM,sizeof(char)*strlen(DB_NAME_LEXEM)), memcpy(&PG_line[strlen(DB_NAME_LEXEM)],name_db,sizeof(char)*len_db), PG_line[strlen(DB_NAME_LEXEM)+len_db]=' ';
  memcpy(&PG_line[1+strlen(DB_NAME_LEXEM)+len_db],APPL_NAME_LEXEM,sizeof(char)*strlen(APPL_NAME_LEXEM)), memcpy(&PG_line[1+strlen(DB_NAME_LEXEM)+len_db+strlen(APPL_NAME_LEXEM)],name_appl,sizeof(char)*len_appl), PG_line[1+strlen(DB_NAME_LEXEM)+len_db+strlen(APPL_NAME_LEXEM)+len_appl]='\0';
  PG_conn=PQconnectdb(PG_line); //Make a connection to the database
  free(PG_line), PG_line=0x0;
  }
//Check to see that the backend connection was successfully made
if (!(PG_conn)) { ylib_errno=YERROR_MEMORY; return FALSE; }
if ( (PQstatus(PG_conn)!=CONNECTION_OK)||((PQ_socket=PQsocket(PG_conn))==-1) )
  {
  if ( (PG_conn)) PQfinish(PG_conn);
  PQ_socket=-1, ylib_errno=YERROR_POSTGRES; return FALSE;
  }
else return TRUE;
}
//This function disconnects from Postgre databse server
void disconnect_postgree_db_server(void)
{
if ( (PG_conn)) PQfinish(PG_conn); 
PQ_socket=-1;
}


//This function sets (non)blocking connection mode with postgre server
char set____blocking_postgre_connection(void)
{
return (!(PQsetnonblocking(PG_conn, 0x0)));
}
char set_nonblocking_postgre_connection(void)
{
return (!(PQsetnonblocking(PG_conn, 0x1)));
}



//This function converts solid binary object into a string of 2 hexadecimal digits per byte
inline void encode_obj_sql_hexadecimal(register size_t size,register char *str,register void *obj) 
{
const char HEXADECIMAL[0x10] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' }; 
while (size--)
  {
  *str=(char)HEXADECIMAL[(unsigned int)(*(unsigned char*)obj)/0x10], str++; 
  *str=(char)HEXADECIMAL[(unsigned int)(*(unsigned char*)obj)%0x10], str++;
  obj++;
  }
}
inline char decode_obj_sql_hexadecimal(register size_t size,register char *str,register void *obj) //str should has assigned str[sizeof(obj)*2]='\0'
{
register char byte;
while (size--)
  {
  switch (*str)
    {
    case (char)'0' : { byte=0x0*0x10; break; }
    case (char)'1' : { byte=0x1*0x10; break; }
    case (char)'2' : { byte=0x2*0x10; break; }
    case (char)'3' : { byte=0x3*0x10; break; }
    case (char)'4' : { byte=0x4*0x10; break; }
    case (char)'5' : { byte=0x5*0x10; break; }
    case (char)'6' : { byte=0x6*0x10; break; }
    case (char)'7' : { byte=0x7*0x10; break; }
    case (char)'8' : { byte=0x8*0x10; break; }
    case (char)'9' : { byte=0x9*0x10; break; }
    case (char)'A' : { byte=0xA*0x10; break; }
    case (char)'B' : { byte=0xB*0x10; break; }
    case (char)'C' : { byte=0xC*0x10; break; }
    case (char)'D' : { byte=0xD*0x10; break; }
    case (char)'E' : { byte=0xE*0x10; break; }
    case (char)'F' : { byte=0xF*0x10; break; }
    case (char)'a' : { byte=0xA*0x10; break; }
    case (char)'b' : { byte=0xB*0x10; break; }
    case (char)'c' : { byte=0xC*0x10; break; }
    case (char)'d' : { byte=0xD*0x10; break; }
    case (char)'e' : { byte=0xE*0x10; break; }
    case (char)'f' : { byte=0xF*0x10; break; }
    default  : { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
    }
  str++;
  switch (*str)
    {
    case (char)'0' : { byte+=0x0; break; }
    case (char)'1' : { byte+=0x1; break; }
    case (char)'2' : { byte+=0x2; break; }
    case (char)'3' : { byte+=0x3; break; }
    case (char)'4' : { byte+=0x4; break; }
    case (char)'5' : { byte+=0x5; break; }
    case (char)'6' : { byte+=0x6; break; }
    case (char)'7' : { byte+=0x7; break; }
    case (char)'8' : { byte+=0x8; break; }
    case (char)'9' : { byte+=0x9; break; }
    case (char)'A' : { byte+=0xA; break; }
    case (char)'B' : { byte+=0xB; break; }
    case (char)'C' : { byte+=0xC; break; }
    case (char)'D' : { byte+=0xD; break; }
    case (char)'E' : { byte+=0xE; break; }
    case (char)'F' : { byte+=0xF; break; }
    case (char)'a' : { byte+=0xA; break; }
    case (char)'b' : { byte+=0xB; break; }
    case (char)'c' : { byte+=0xC; break; }
    case (char)'d' : { byte+=0xD; break; }
    case (char)'e' : { byte+=0xE; break; }
    case (char)'f' : { byte+=0xF; break; }
    default  : { ylib_errno=YERROR_EXTERNAL_CODE; return FALSE; }
    }
  str++;
  *(char*)obj=byte, obj++;
  }
return TRUE;
}


//This function converts str structure into a string of 2 hexadecimal digits per byte
//Note. it doesn't encode molecule name as DB has a separate field for it. Additionally edges are compressed, the atomn id filed is choosen either char if (natoms<=0xFF), short (if natoms<=0xFFF) and as int olny otherwise
typedef struct { uint8_t  v0, v1; char type; } t_cedge;
typedef struct { uint16_t v0, v1; char type; } t_sedge;
typedef struct { uint32_t v0, v1; char type; } t_iedge;
inline size_t mono_str_line_len(register t_str *str) 
{
     if (str->natoms<=0xFF  ) return 0x2*(sizeof(uint32_t)*0x2+str->natoms*sizeof(char)+str->nedges*sizeof(t_cedge))+0x1;
else if (str->natoms<=0xFFFF) return 0x2*(sizeof(uint32_t)*0x2+str->natoms*sizeof(char)+str->nedges*sizeof(t_sedge))+0x1;
else                          return 0x2*(sizeof(uint32_t)*0x2+str->natoms*sizeof(char)+str->nedges*sizeof(t_iedge))+0x1;
}
inline char encode_mono_str_sql_hexadecimal(register char *line,register t_str *str) //NB! str here is "structure". The encode: | name_len | natoms | nedges | name | a | edges |
{
register unsigned int _i;
t_cedge cedge;
t_sedge sedge;
t_iedge iedge;
if ( ( (str->ress))||(!(str->natoms)) ) { ylib_errno=YERROR_LEGAL; return FALSE; } //This function is for decoding of 1-residue compounds only
encode_obj_sql_hexadecimal(sizeof(uint32_t), line,                      &str->natoms);
encode_obj_sql_hexadecimal(sizeof(uint32_t),&line[sizeof(uint32_t)*0x2],&str->nedges);
encode_obj_sql_hexadecimal(str->natoms,     &line[sizeof(uint32_t)*0x4],str->a);
     if (str->natoms<=0xFF  ) 
       {
       for (_i=0; _i!=str->nedges; _i++)
         { 
         cedge.v0=(uint8_t)str->edges[_i].vertice[0], cedge.v1=(uint8_t)str->edges[_i].vertice[1],  cedge.type=(uint8_t)str->edges[_i].type;
         encode_obj_sql_hexadecimal(sizeof(t_cedge),&line[0x2*(sizeof(uint32_t)*0x2+str->natoms+_i*sizeof(t_cedge))],(void*)&cedge);
         }
       line[0x2*(sizeof(uint32_t)*0x2+str->natoms+str->nedges*sizeof(t_cedge))]='\0';
       }
else if (str->natoms<=0xFFFF) 
       {
       for (_i=0; _i!=str->nedges; _i++) 
         {
         sedge.v0=(uint16_t)str->edges[_i].vertice[0], sedge.v1=(uint16_t)str->edges[_i].vertice[1], sedge.type=(char)str->edges[_i].type;
         encode_obj_sql_hexadecimal(sizeof(t_sedge),&line[0x2*(sizeof(uint32_t)*0x2+str->natoms+_i*sizeof(t_sedge))],(void*)&sedge);
         }
       line[0x2*(sizeof(uint32_t)*0x2+str->natoms+str->nedges*sizeof(t_sedge))]='\0';
       }
else   {
       for (_i=0; _i!=str->nedges; _i++)
         { 
         iedge.v0=(uint32_t)str->edges[_i].vertice[0], iedge.v1=(uint32_t)str->edges[_i].vertice[1], iedge.type=(char)str->edges[_i].type;
         encode_obj_sql_hexadecimal(sizeof(t_iedge),&line[0x2*(sizeof(uint32_t)*0x2+str->natoms+_i*sizeof(t_iedge))],(void*)&iedge);
         }
       line[0x2*(sizeof(uint32_t)*0x2+str->natoms+str->nedges*sizeof(t_iedge))]='\0';
       } 
return TRUE; 
}
inline size_t line_mono_str_len(register char *line) 
{
uint32_t natoms, nedges;
decode_obj_sql_hexadecimal(sizeof(uint32_t), line,                      &natoms);
decode_obj_sql_hexadecimal(sizeof(uint32_t),&line[sizeof(uint32_t)*0x2],&nedges);
return sizeof(t_str)+natoms*sizeof(unsigned int)+((!(natoms%sizeof(int))) ? natoms : sizeof(unsigned int)*(1+natoms/sizeof(int)))+nedges*sizeof(t_edge);
}
inline char decode_mono_str_sql_hexadecimal(register t_str *str,register char *line) //NB! str here is "string"; str should has natoms assigned before the call
{
uint32_t i, j;
t_cedge cedge;
t_sedge sedge;
t_iedge iedge;
decode_obj_sql_hexadecimal(sizeof(uint32_t), line,                      &i);
decode_obj_sql_hexadecimal(sizeof(uint32_t),&line[sizeof(uint32_t)*0x2],&j);
str->natoms=(unsigned int)i, str->nedges=(unsigned int)j;
if (!(str->natoms)) { ylib_errno=YERROR_LEGAL; return FALSE; } else { str->name=0x0, str->ress=0x0, str->rsize=0x0, str->r=0x0, str->start_rid=0; }
str->a=(char*)((void*)str+sizeof(t_str));
str->anames=(char(*)[sizeof(int)])((void*)str->a+((!(str->natoms%sizeof(int))) ? str->natoms : sizeof(int)*(1+str->natoms/sizeof(int))));
str->edges=(t_edge*)((void*)str->anames+sizeof(int)*str->natoms);
decode_obj_sql_hexadecimal(str->natoms,     &line[sizeof(uint32_t)*0x4],str->a);
     if (str->natoms<=0xFF  )
       for (i=0; i<str->nedges; i++)   
         {
         decode_obj_sql_hexadecimal(sizeof(t_cedge),&line[sizeof(uint32_t)*0x4+str->natoms*0x2+sizeof(t_cedge)*i*0x2],(void*)&cedge);
         str->edges[i].vertice[0]=(unsigned int)cedge.v0, str->edges[i].vertice[1]=(unsigned int)cedge.v1, str->edges[i].type=(int)cedge.type;
         }
else if (str->natoms<=0xFFFF)
       for (i=0; i<str->nedges; i++)   
         {
         decode_obj_sql_hexadecimal(sizeof(t_sedge),&line[sizeof(uint32_t)*0x4+str->natoms*0x2+sizeof(t_sedge)*i*0x2],(void*)&sedge);
         str->edges[i].vertice[0]=(unsigned int)sedge.v0, str->edges[i].vertice[1]=(unsigned int)sedge.v1, str->edges[i].type=(int)sedge.type;
         }
else   for (i=0; i<str->nedges; i++)   
         {
         decode_obj_sql_hexadecimal(sizeof(t_iedge),&line[sizeof(uint32_t)*0x4+str->natoms*0x2+sizeof(t_iedge)*i*0x2],(void*)&iedge);
         str->edges[i].vertice[0]=(unsigned int)iedge.v0, str->edges[i].vertice[1]=(unsigned int)iedge.v1, str->edges[i].type=(int)iedge.type;
         }
return TRUE;
}






