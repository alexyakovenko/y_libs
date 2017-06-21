 //Text operations routines

#ifndef Y_SYSTEM
#include "y_system.h"
#endif

#define YNAME_LEN 0x1F

//Type data checking
//This function define is given char a number
#define CHECK_NUMBER_TYPE(c) ( (c=='0')||(c=='1')||(c=='2')||(c=='3')||(c=='4')||(c=='5')||(c=='6')||(c=='7')||(c=='8')||(c=='9') )

//Redefine operator macross
#define STR_OPS(A,OPERATOR,B) (      IF (OPERATOR=='+') AddString(A,B)   \
                                ELSE IF (OPERATOR=='-') SubString(A,B)   \
                                ELSE IF (OPERATOR=='?') FindPattern(A,B) )

//This script delete commants from char buffer pointer using temp_ui temporary unsigned int variable
#define CUT_COMMENTS(nlen,buffer) do{                                                                            \
                                    nlen=0x0;                                                                    \
                                    while (buffer[nlen]!='\0')                                                   \
                                      {                                                                          \
                                      if ( (buffer[nlen]=='#')||((buffer[nlen]=='/')&&(buffer[nlen+1]=='/')) )   \
                                        {                                                                        \
                                        buffer[nlen]='\0';                                                       \
                                        break;                                                                   \
                                        }                                                                        \
                                      else nlen++;                                                               \
                                      }                                                                          \
                                    }while(0)

//This script creates name string from the prefix and the suffix
#define MAKE_NAME(_i,_j,name,prefix_name,suffix_name) do{                                            \
                                                        _j=0;                                        \
                                                        while (prefix_name[_j])                      \
                                                          name[_j]=prefix_name[_j++];                \
                                                        _i=0;                                        \
                                                        do                                           \
                                                          name[_j++]=suffix_name[_i];                \
                                                        while (suffix_name[_i++]);                   \
                                                        }while(0)

//This script cut spaces from 'int' names
#define MAKE_INT_NAME(_i,_j,int_name,name) do{ _i=0, _j=0; while ( (_i!=sizeof(int))&&( (name[_i]==' ')||(name[_i]=='\t') ) ) _i++;                               \
                                             while (_i!=sizeof(int)&&(name[_i])&&(name[_i]!=' ')&&(name[_i]!='\t')&&(name[_i]!='\n') ) int_name[_j++]=name[_i++]; \
                                             while (_j!=sizeof(int)) int_name[_j++]=' '; }while(0)
                                             

// -------------------- S I M P L E    S T R I N G --------------------
//This function search for new string '\n' in the buffer till '\0'
inline char check_string_end(register char *buffer);

//This function check the lexem presence
char check_lexem(register unsigned int lexem_number,register char *string);

//This function gets the number of lexems in string
inline unsigned int get_lex_num(register char *string);

//This function finds begin of the lexem of an asked number in the sting or return 0 otherwise
char *get_lex (register unsigned int lexem_number,register char *string);

//This function calculates lexem lenght
unsigned int lexlen(char *lexem);

//This function compates lexem alike str_cmp
inline char lex_cmp(register char *lexem1,register char *lexem2);

//This function compares couple strings
inline char str_cmp(register char *string1,register char *string2);

//This function compares part of couple strings
inline char strn_cmp(register char *string1,register char *string2,register unsigned int n);

//This function remove odd symbols from start and end of string
//Note string and lexem might be the same
void cut_odd_symbols(char *lexem,char *string);

//This function check is given symbol is a number symbol or not
char check_number_type(register char c);

//This function finds the begin of the first lexem in the string or return 0x0 if there in no the first lexem
inline char *get_first_lex(char *string);

//This function check about accordance of string data to double (in accorance with postgre rules)
//Rules (where digits are '0', '1', '2', '3', '4', '5', '6', '7', '8' and '9'):
// 1. [+/-]digits['.'[digits]['e'[+/-]digits]]
// 2. [+/-].digits['e'[+/-]digits]
inline char check_real_type(register char *lexem);
#define check_double_type check_real_type 

//This gunction check is it a int
/* Rules:
1. All simbols numbers or sign if firs position
*/
char check_int_type(register char *lexem);
//This gunction check is it a unsigned int
char check_uint_type(register char *lexem);

//These function sets ylib_errno=YERROR_SUSPICIOUS or ylib_errno=YERROR_OK to indicate success of the conversion.
//This function convert decimal string into int
inline int strntoi(register unsigned int len,register char *str);
//This function convert decimal string into double 
inline double strntod(register unsigned int len,register char *str);


//All to upper case
void all_to_upper(register char *string);
void lex_to_upper(register char *lexem);

//All to lower case
void all_to_lower(register char *string);
void lex_to_lower(register char *lexem);

// ----------------------- D U A L    S T R I N G -----------------------------

//This function add the lexem to end of string. It return the lenght of combined string.
unsigned int add_string(register char *A,register char *B);

//This function subtrace string B from string A, it return the lenght of string A after operation
unsigned int sub_string(register char *A,register char *B);

//This function looking for pattern in string
//It return the letter number of first pattern character in the sting or (unsigned int)-1 if the pattern was not found
unsigned int find_pattern(register char *pattern,register char *string);

//This function extend (from lexem) or compress (into lexem) line to given lenght
void equalize_string(register char* lexem,unsigned int lenght,char *string);

//This function exchange two zones of memory. It's dangerous: due to speedup reasons the function doesn't check if the exchanged areas overlap.
inline void memexchange(register void *a, register void *b, register size_t size);

// ------------------- T R I P L E    S T R I N G -----------------------

//This function exchange all patterns of pattern1 to pattern2
//It return the number of maded changes
unsigned int replace_pattern(char *pattern1,char *pattern2,char *string);

// ------------------------- Q U A D R O    S T I N G ---------------------------

//This function take value stored between first found pattern1 and first found pattern2
unsigned int find_lexem_inside_patterns(unsigned int *start,char *pattern1,char *pattern2,char *string);

