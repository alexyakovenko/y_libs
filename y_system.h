#define Y_SYSTEM 0x1

#define Y_MAGIC 2004

#include <limits.h> 
#include <float.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>   //include stack
#include <stdint.h>   //include support for (u)int8_t, (u)int16_t, (u)int32_t etc
#include <stdio.h>    //include files
#include <stdlib.h>   //include memory
#include <string.h>   //include char massives
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <inttypes.h>

unsigned int ylib_errno;  //Code of error occured during calling of ylib functions

//Computer constants

#define TRUE   1
#define FALSE  0
#define NTNF  -1

#define ERROR_CODE -1

#define MAX_UINT 0xFFFFFFFF
#define MAX_INT  2147483647
#define MAX_DOUBLE +1.e72
#define MAX_CHAR 0xFF
#define MAX_INT_STR_LEN 0xA

#define MIN_DOUBLE -MAX_DOUBLE
#define PRECISION TINY

//Physical constants

//The Boltzmann constant, measured as J/K
#define THERMODYNAMIC_k     1.380648813e-23 
//The Avogadro constant, measured as 1/mol
#define THERMODYNAMIC_N     6.0221417930e+23 
//A temperature of human body, measured as K (310K=36.85C)
#define THERMODYNAMIC_T   310.00
#define THERMODYNAMIC_TT   0.310
//Gas R constant measured as  J/K/mol
#define THERMODYNAMIC_R     8.314459648
//Volume per molecule at concentration of 1M, measured as A^3
#define MOLAR_BOX_VOLUME 1660.00

//ylibs_errno values description:
// System level errors
#define YERROR_OK                  0  // #0 all are correct 
#define YERROR_MEMORY              1  // #1 memory allocation/reallocation problem 
#define YERROR_IO                  2  // #2 input/output problem
#define YERROR_USER               14  // #14 error specefied by external (relatively to source code) reasons
#define YERROR_NIMPLEMENTED       15  // #15 requested function is not implemented (yet)
// Algorithms-related errors
#define YERROR_LEGAL              16  // #16 an algorithm run into its natural limitations - all is correct, but some other method/algorithm need to be used
#define YERROR_INTERNAL_CODE      17  // #17 an algorithm works as it shouldn`t be
#define YERROR_EXTERNAL_CODE      18  // #18 an algorithm fail due to external (relatively to particular algorithm) reasons
#define YERROR_NCONVERGED         19  // #19 an algorithm worked as it should be but requested exit conditions was not meet (yet)
#define YERROR_COMPLEXITY         20  // #20 an algorithm run into complexity limitations: the calculations aren't gona be done before the end of the days... 
#define YERROR_SUSPICIOUS         21  // #21 an algorithm worked as it should be but results looks like there is an error occured
#define YERROR_IMPOSSIBLE         22  // #22 an algorithm run into conditions that considered as imposible during development 
// Data-related errors
#define YERROR_DATA_CONSISTMENT   64  // #64 error occured in ambundant data check
#define YERROR_DATA_DUPLICATE     65  // #65 error occured in (extended) data analysis due to detected duplication of data records
#define YERROR_DATA_FORMAT        66  // #66 error occured in (input) data formatting
// Communication errors 
#define YERROR_POSITIVE_OUTCOME  128  // #128 everything is fine, the request is successfully satisfied   
#define YERROR_NEGATIVE_OUTCOME  129  // #129 everything is fine, but the request can't be satisfied (say there is no such record is found in db)
#define YERROR_NO_SERVER         130  // #130 there is no required server
#define YERROR_NO_YSERVER        131  // #131 there is no y-server      
// External libraries errors
#define YERROR_POSTGRES          256  // #256 an error occured inside of Postgres library 

#define YERR_OK   YERROR_OK 


char buffer[0xFF];
extern const char const * Y_APPLICATION;


//This function performs error exit
inline void error_exit(register char *templatee, ... );

//yprintf statuses
#define YPRINTF_NOTHING               0 //just print the message
#define YPRINTF_INFO                  1
#define YPRINTF_WARNING               2
#define YPRINTF_ERROR                 3
#define YPRINTF_NOTE                  4

//This function simply pass the data to printf.
//It is assumed that it will pass the message to the server in future.
inline void yprintf(register unsigned int yprintf_status,register char *templatee, ... );

//This is debugging function
char * const get_errno(unsigned int _errno);

//This is debugging function for ylib errors
char * const get_yerrno(unsigned int _yerrno);



