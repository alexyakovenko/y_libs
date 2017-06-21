#include "y_system.h"


//This function performs error exit
inline void error_exit(char *templatee, ... )
{
va_list stack;
va_start(stack,templatee);
printf("ERROR. ");
vprintf(templatee,stack);
printf("\n");
va_end(stack);
exit(FALSE);
}

//This function simply pass the data to printf.
//It is assumed that it will pass the message to the server in future.
inline void yprintf(unsigned int yprintf_status,char *templatee, ... )
{
va_list stack;
va_start(stack,templatee);
switch (yprintf_status)
  {
  case YPRINTF_NOTHING     : {                                        break; } 
  case YPRINTF_INFO        : { printf("%s INFO. ",    Y_APPLICATION); break; }
  case YPRINTF_WARNING     : { printf("%s WARNING. ", Y_APPLICATION); break; }
  case YPRINTF_ERROR       : { printf("%s ERROR. ",   Y_APPLICATION); break; }
  case YPRINTF_NOTE        : { printf("%s NOTE. ",    Y_APPLICATION); break; }
  default                  : ; //show nothing extra
  }
vprintf(templatee,stack);
va_end(stack);
}

//This function shows error code
char * const get_errno(unsigned int _errno)
{
switch (_errno)
  {
  case EACCES       : return "EACCES";
  case EAGAIN       : return "EAGAIN";
  case EBADF        : return "EBADF";
  case ECONNABORTED : return "ECONNABORTED";
  case EFAULT       : return "EFAULT";
  case EIDRM        : return "EIDRM";
  case EINTR        : return "EINTR";
  case EINVAL       : return "EINVAL";
  case EMFILE       : return "EMFILE";
  case ENFILE       : return "ENFILE";
  case ENOTSOCK     : return "ENOTSOCK";
  case EOPNOTSUPP   : return "EOPNOTSUPP";
  case EPERM        : return "EPERM";
  case ERANGE       : return "ERANGE";
  default           : return "Unknown error code received, sorry.";
  }
}

//This function shows yerror code
char * const get_yerrno(unsigned int _yerrno)
{
switch (_yerrno)
  {
// System level errors
  case YERROR_OK               : return "YERROR_OK";                     // #0 all are correct 
  case YERROR_MEMORY           : return "YERROR_MEMORY";                 // #1 memory allocation/reallocation problem 
  case YERROR_IO               : return "YERROR_IO";                     // #2 input/output problem
  case YERROR_USER             : return "YERROR_USER";                   // #14 error specefied by external (relatively to source code) reasons
  case YERROR_NIMPLEMENTED     : return "YERROR_NIMPLEMENTED";           // #15 requested function is not implemented (yet)
// Algorithms-related errors
  case YERROR_LEGAL            : return "YERROR_LEGAL";                  // #16 an algorithm run into its natural limitations - all is correct, but some other method/algorithm need to be used
  case YERROR_INTERNAL_CODE    : return "YERROR_INTERNAL_CODE";          // #17 an algorithm works as it shouldn`t be
  case YERROR_EXTERNAL_CODE    : return "YERROR_EXTERNAL_CODE";          // #18 an algorithm fail due to external (relatively to particular algorithm) reasons
  case YERROR_NCONVERGED       : return "YERROR_NCONVERGED";             // #19 an algorithm worked as it should be but requested exit conditions was not meet (yet)
  case YERROR_COMPLEXITY       : return "YERROR_COMPLEXITY";             // #20 an algorithm run into complexity limitations: the calculations aren't gona be done before the end of the days... 
  case YERROR_SUSPICIOUS       : return "YERROR_SUSPICIOUS";             // #21 an algorithm worked as it should be but results looks like there is an error occured
  case YERROR_IMPOSSIBLE       : return "YERROR_IMPOSSIBLE";             // #22 an algorithm run into conditions that considered as imposible during development 
// Data-related errors
  case YERROR_DATA_CONSISTMENT : return "YERROR_DATA_CONSISTMENT";       // #64 error occured in ambundant data check
  case YERROR_DATA_DUPLICATE   : return "YERROR_DATA_DUPLICATE";         // #65 error occured in (extended) data analysis due to detected duplication of data records
  case YERROR_DATA_FORMAT      : return "YERROR_DATA_FORMAT";            // #66 error occured in (input) data formatting
// Communication errors
  case YERROR_POSITIVE_OUTCOME : return "YERROR_POSITIVE_OUTCOME";       // #128 everything is fine, the request is successfully satisfied   
  case YERROR_NEGATIVE_OUTCOME : return "YERROR_NEGATIVE_OUTCOME";       // #129 everything is fine, but the request can't be satisfied (say there is no such record is found in db)
  case YERROR_NO_SERVER        : return "YERROR_NO_ SERVER";             // #130 there is no required server
  case YERROR_NO_YSERVER       : return "YERROR_NO_YSERVER";             // #131 there is no y-server      
// External libraries errors
  case YERROR_POSTGRES         : return "YERROR_POSTGRES";               // #256 an error occured inside of Postgres library 
  default                      : return "An unknown error code received, sorry.";
  }
}

#define YERR_OK   YERROR_OK 



