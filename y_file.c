#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include "y_file.h"

//This module perform read/write file operatioons with various file format of y_so programm


//----------------- U N I V E R S A L    F I L E     O P E R A T I O N S ------------------------

//This function reads arbitrary length string into memory.
//It returns a pointer to buffer[] (from y_systems), if the string fits into it (thus avoiding allocation and being almost as quick as fgets()), or store the string in (re-)allocated memory to ensure that it was read completely.
//NB! It allocates memory in 0xFF quants and may require trimming if the string are stores
inline char *yfgets(FILE *in)
{
register unsigned int _l;
register size_t size;
register char *cp, *string;
extern char buffer[0xFF];
//Stage I. Quickly try to upload a line into the buffer and return it
if (!(fgets(buffer,0xFF,in))) return 0x0;
if ( (strlen(buffer)!=0xFF-0x1)||(buffer[0xFF-0x2]=='\n') ) return buffer;
//Stage II. Slow memory (re-allocation) process with extern char buffer[0xFF] as a reading buffer
if (!(string=(char*)malloc(sizeof(char)*0xFF))) { ylib_errno=YERROR_MEMORY; return 0x0; }
else { memcpy(string,buffer,sizeof(char)*0xFF); size=0xFF-0x1; }
while ( (fgets(buffer,0xFF,in)))
  {
  if ( ((_l=strlen(buffer))!=0xFF-0x1)||(buffer[0xFF-0x2]=='\n') ) 
    {
    if (!(cp=(char*)realloc(string,sizeof(char)*(++_l+size)))) { free(string); ylib_errno=YERROR_MEMORY; return 0x0; }
    else { string=cp; memcpy(&string[size],buffer,_l); return string; }
    }
  //Even more room is needed
  if (!(cp=(char*)realloc(string,sizeof(char)*(size+0xFF)))) { free(string); ylib_errno=YERROR_MEMORY; return 0x0; } 
  else { string=cp; memcpy(&string[size],buffer,0xFF-0x1); size+=0xFF-0x1; }
  }
return string;
}

//This function read file while it possiable and '\n' not readed
inline char fstrskip(FILE *in)
{
extern char buffer[0xFF];
while ( (fgets(buffer,0xFE,in))) if (buffer[strlen(buffer)-1]=='\n') return TRUE;
return FALSE;
}

//This function copy string from one file into other;
char fstrcpy(FILE *in,FILE *out)
{
extern char buffer[];
register unsigned int _l;
while ( (fgets(buffer,0xFE,in)))
  {
  _l=strlen(buffer);  
  if (fwrite(buffer,sizeof(char),_l,out)!=_l) { ylib_errno=YERROR_IO; return FALSE; }
  if (buffer[_l-1]=='\n') break;
  }
return TRUE;
}


