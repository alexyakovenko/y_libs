//This file contain data about file operations.
#define Y_FILE 0x1

#ifndef Y_SYSTEM
#include "y_system.h"
#endif


//----------------- U N I V E R S A L    F I L E     O P E R A T I O N S ------------------------

//This function reads arbitrary length string into memory.
//It returns a pointer to buffer[] (from y_systems) if the string fits in and thus avoid allocation being almost as quick as fgets() or store the string in (re-)allocated memory to ensure that it was read completely
inline char *yfgets(FILE *in);

//This function read file while it possible and '\n' not read
inline char fstrskip(FILE *in);

//This function copy string from one file into other;
char fstrcpy(FILE *in,FILE *out);

