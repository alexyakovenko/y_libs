#include "y_txt.h"
#include <limits.h>
#include <math.h>


// -------------------- S I M P L E    S T R I N G --------------------

//This function search for new string in the buffer or upload it from file
inline char check_string_end(register char *buff)
{
do { if (*buff=='\n') return TRUE; } while ( (*buff++)); return FALSE;
}

//This function check the lexem presence
inline char check_lexem (register unsigned int lexem_number,register char *string)
{
register unsigned int n=0;

while ( (string[n])&&(string[n]!='\n') )
  if ( (string[n]==' ')||(string[n]=='\t')||(string[n]=='\r') )
    n++;
  else
    {
    if (!(--lexem_number))
      return TRUE;
    else
      do n++; while ( (string[n])&&(string[n]!='\n')&&(string[n]!=' ')&&(string[n]!='\t')&&(string[n]!='\r') );
    }
return FALSE;
}

//This function gets the number of lexems in string
inline unsigned int get_lex_num(register char *string)
{
register unsigned int n=0, lexem_number=0;

while ( (string[n])&&(string[n]!='\n') )
  if ( (string[n]==' ')||(string[n]=='\t')||(string[n]=='\r') )
    n++;
  else
    {
    while ( (string[n])&&(string[n]!='\n')&&(string[n]!=' ')&&(string[n]!='\t')&&(string[n]!='\r') ) 
      n++;
    lexem_number++;
    }
return lexem_number;
}

//This function finds begin of the lexem of an asked number in the sting or return 0 otherwise
inline char *get_lex(register unsigned int lexem_number,register char *string)
{
register unsigned int n=0;

while ( (string[n])&&(string[n]!='\n'))
  if ( (string[n]==' ')||(string[n]=='\t')||(string[n]=='\r') )
    n++;
  else
    {   
    if (!(--lexem_number))
      return &string[n];
    else
      while ( (string[n])&&(string[n]!='\n')&&(string[n]!=' ')&&(string[n]!='\t')&&(string[n]!='\r') )
        n++;
    }
return FALSE;
}

//This function calculates lexem lenght
inline unsigned int lexlen(register char *lexem)
{
register unsigned int _n=0;
while ( ( (*lexem))&&(*lexem!='\n')&&(*lexem!=' ')&&(*lexem!='\t')&&(*lexem!='\r') ) { lexem++, _n++; }
return _n;
}

//This function finds the begin of the first lexem in the string or return FALSE if there in no the first lexem
inline char *get_first_lex(register char *lexem)
{
if (!(lexem)) return FALSE;
while ( ( (*lexem))&&(*lexem!='\n') ) if ( (*lexem!=' ')&&(*lexem!='\t')&&(*lexem!='\r') ) return lexem; else lexem++;
return FALSE;
}

//This function compates lexem alike str_cmp
inline char lex_cmp(register char *lexem1,register char *lexem2)
{
while ((*lexem1!=' ')&&(*lexem1!='\0')&&(*lexem1!='\t')&&(*lexem1!='\r')&&(*lexem1!='\n'))
  if (*lexem1!=*lexem2)
    return FALSE;
  else
    {
    lexem1++;
    lexem2++;
    }
return ((*lexem2==' ')||(*lexem2=='\0')||(*lexem2=='\t')||(*lexem2=='\r')||(*lexem2=='\n'));
}

//This function compares couple strings
inline char str_cmp(register char *string1,register char *string2)
{
while (*string1)
  if (*string1!=*string2)
    return FALSE;
  else
    {
    string1++;
    string2++;
    }
return (!*string2);
}

//This function compares part of couple strings
inline char strn_cmp(register char *string1,register char *string2,register unsigned int n)
{
while ((*string1)&&(n))
  if (*string1!=*string2)
    return FALSE;
  else 
    {
    string1++;
    string2++;
    n--;
    }
return (*string2==*string1);
}

//This function remove odd symbols from start and end of string
//Note string and lexem might be the same
void cut_odd_symbols(register char *lexem,register char *string)
{
register unsigned int _i;
_i=strlen(string);
while ((--_i))
  if ( ((string[_i]!='\r')&&(string[_i]!='\t')&&(string[_i]!='\n')) )
    break;
if (_i)
  {
  memcpy(lexem,string,_i+1);
  lexem[_i+1]='\0';
  }
else
  lexem[0]='\0';
}



//This function check about accordance of string data to double (in accorance with postgre rules)
/* Rules (where digits are '0', '1', '2', '3', '4', '5', '6', '7', '8' and '9'):
1. [+/-]digits['.'[digits]['e'[+/-]digits]]
2. [+/-].digits['e'[+/-]digits]
*/
inline char check_real_type(register char *lexem)
{
register unsigned int _n;
register char flag;
if ( (!(lexem=get_first_lex(lexem)))||(!(_n=lexlen(lexem))) ) return FALSE; else flag=FALSE;
if ( (*lexem=='-')||(*lexem=='+') ) lexem++, _n--;
while ( ( (_n))&&( (CHECK_NUMBER_TYPE(*lexem))) ) { lexem++, _n--, flag=TRUE; }
if (!(_n)) return flag;
if (*lexem=='.') { lexem++, _n--; while ( ( (_n))&&( (CHECK_NUMBER_TYPE(*lexem))) ) { lexem++, _n--, flag=TRUE; } }
if (!(_n)) return flag;
if ( (!(flag))||( (*lexem!='e')&&(*lexem!='E') ) ) return FALSE;
lexem++, _n--, flag=FALSE;
if ( (*lexem=='-')||(*lexem=='+') ) lexem++, _n--;
while ( ( (_n))&&( (CHECK_NUMBER_TYPE(*lexem))) ) { lexem++, _n--, flag=TRUE; }
return (!(_n)) ? flag : FALSE; 
}


//This gunction check is it a int
/* Rules:
1. All simbols numbers or a sign in the first position
*/
char check_int_type(register char *lexem)
{
register unsigned int n=0;
while ((lexem[n]==' ')||(lexem[n]=='\t')||(lexem[n]=='\r'))
  n++;
if ((!(CHECK_NUMBER_TYPE(lexem[n])))&&(lexem[n]!='+')&&(lexem[n]!='-')) return FALSE;
while (lexem[++n])
  if (!(CHECK_NUMBER_TYPE(lexem[n])))
    {
    if ((!n)||(lexem[n-1]=='+')||(lexem[n-1]=='-')||((lexem[n]!=' ')&&(lexem[n]!='\t')&&(lexem[n]!='\r')&&(lexem[n]!='\n')))
      return FALSE;
    else
      break;
    }
return TRUE;
}

//This gunction check is it a unsigned int
char check_uint_type(register char *lexem)
{
register unsigned int n=0;
while ((lexem[n]==' ')||(lexem[n]=='\t')||(lexem[n]=='\r'))
  n++;
if ((!(CHECK_NUMBER_TYPE(lexem[n])))&&(lexem[n]!='+')) return FALSE;
while (lexem[++n])
  if (!(CHECK_NUMBER_TYPE(lexem[n])))
    {
    if ((!n)||(lexem[n-1]=='+')||((lexem[n]!=' ')&&(lexem[n]!='\t')&&(lexem[n]!='\r')&&(lexem[n]!='\n')))
      return FALSE;
    else
      break;
    }
return TRUE;
}

//This function convert decimal string into int
inline int strntoi(register unsigned int len,register char *str)
{
extern unsigned int ylib_errno;
register unsigned int _i, _s, _n;
if (!(len)) return 0;
while ( (isspace(*str)))  {  str++; if (!(--len)) { ylib_errno=YERROR_SUSPICIOUS; return 0; } }
_s=0; if (*str=='+') {       str++; if (!(--len)) { ylib_errno=YERROR_SUSPICIOUS; return 0; } }
else  if (*str=='-') { _s=1, str++; if (!(--len)) { ylib_errno=YERROR_SUSPICIOUS; return 0; } }
TRIM_ZEROES: ; //Cut leading zeroes 
switch (*str)
  {
  case '0' : { str++; if (!(--len)) { ylib_errno=YERROR_SUSPICIOUS; return 0; } else goto TRIM_ZEROES; }
  case '1' : { _i=1; break; }
  case '2' : { _i=2; break; }
  case '3' : { _i=3; break; }
  case '4' : { _i=4; break; }
  case '5' : { _i=5; break; }
  case '6' : { _i=6; break; }
  case '7' : { _i=7; break; }
  case '8' : { _i=8; break; }
  case '9' : { _i=9; break; }
  default  : { ylib_errno=YERROR_SUSPICIOUS; return 0; } //Not a digit
  }
while (len--)
  {
  switch (*++str)
    {
    case '0' : { _n=0; break; }
    case '1' : { _n=1; break; }
    case '2' : { _n=2; break; }
    case '3' : { _n=3; break; }
    case '4' : { _n=4; break; }
    case '5' : { _n=5; break; }
    case '6' : { _n=6; break; }
    case '7' : { _n=7; break; }
    case '8' : { _n=8; break; }
    case '9' : { _n=9; break; }
    default  : goto LABEL_RETURN_VALUE; //Not a digit
    }
  if (INT_MAX/_i>=10) _i*=10; else { ylib_errno=YERROR_SUSPICIOUS; return (!(_s)) ? INT_MAX : INT_MIN; } //int overflow
  if (INT_MAX-_i>=_n) _i+=_n; else { ylib_errno=YERROR_SUSPICIOUS; return (!(_s)) ? INT_MAX : INT_MIN; } //int overflow
  }
LABEL_RETURN_VALUE: ylib_errno=YERROR_OK; return (!(_s)) ? _i : -_i;
}

//This function convert decimal string into double 
inline double strntod(register unsigned int len,register char *str)
{
extern unsigned int ylib_errno;
register double _d;
register unsigned int _n;
register int _es, _s, _p;
if (!(len)) return 0.;
while ( (isspace(*str)))  {  str++; if (!(--len)) { ylib_errno=YERROR_SUSPICIOUS; return 0.; } }
_s=0; if (*str=='+') {       str++; if (!(--len)) { ylib_errno=YERROR_SUSPICIOUS; return 0.; } }
else  if (*str=='-') { _s=1, str++; if (!(--len)) { ylib_errno=YERROR_SUSPICIOUS; return 0.; } }
_p=0;
TRIM_ZEROES_VAL: ; //Cut leading zeroes 
switch (*str)
  {
  case '.' : { if (!(_p)) { ylib_errno=YERROR_SUSPICIOUS; return 0.; } else _p=1; }
  case '0' : { if ( (_p)) { if (_p==INT_MIN) { ylib_errno=YERROR_OK; return 0.; } else _p--; } str++; if (!(--len)) { ylib_errno=YERROR_SUSPICIOUS; return 0.; } else goto TRIM_ZEROES_VAL; }
  case '1' : { _d=1.; break; }
  case '2' : { _d=2.; break; }
  case '3' : { _d=3.; break; }
  case '4' : { _d=4.; break; }
  case '5' : { _d=5.; break; }
  case '6' : { _d=6.; break; }
  case '7' : { _d=7.; break; }
  case '8' : { _d=8.; break; }
  case '9' : { _d=9.; break; }
  default  : { ylib_errno=YERROR_SUSPICIOUS; return 0.; } //Not a digit
  }
while (len--)
  {
  str++;
  if ( (*str=='e')||(*str=='E') )
    {
    _es=0; if (*str=='+') {        str++; if (!(--len)) goto LABEL_RETURN_VALUE; }
    else   if (*str=='-') { _es=1, str++; if (!(--len)) goto LABEL_RETURN_VALUE; }
    if (len--)
      {
      TRIM_ZEROES_EXP: ; //Cut leading exp zeroes 
      switch (*str++)
        {
        case '0' : { if (!(--len)) goto LABEL_RETURN_VALUE; else goto TRIM_ZEROES_EXP; }
        case '1' : { _n=1; break; }
        case '2' : { _n=2; break; }
        case '3' : { _n=3; break; }
        case '4' : { _n=4; break; }
        case '5' : { _n=5; break; }
        case '6' : { _n=6; break; }
        case '7' : { _n=7; break; }
        case '8' : { _n=8; break; }
        case '9' : { _n=9; break; }
        default  : goto LABEL_RETURN_VALUE; //Not a number
        }
      if (!(_es)) { if (INT_MAX-_p<_n) { ylib_errno=YERROR_OK; return (!(_s)) ? INFINITY : -INFINITY; } else _p+=_n; }
      else        { if (_p-INT_MIN<_n) { ylib_errno=YERROR_OK; return 0.;                             } else _p-=_n; }
      while (len--)
        {
        switch (*str++)
          {
          case '0' : { _n=0; break; }
          case '1' : { _n=1; break; }
          case '2' : { _n=2; break; }
          case '3' : { _n=3; break; }
          case '4' : { _n=4; break; }
          case '5' : { _n=5; break; }
          case '6' : { _n=6; break; }
          case '7' : { _n=7; break; }
          case '8' : { _n=8; break; }
          case '9' : { _n=9; break; }
          default  : goto LABEL_RETURN_VALUE; //Not a numkber
          }
        if (!(_es)) 
          {
          if (INT_MAX/_p<10) { ylib_errno=YERROR_OK; return (!(_s)) ? INFINITY : -INFINITY; } else _p*=10;
          if (INT_MAX-_p<_n) { ylib_errno=YERROR_OK; return (!(_s)) ? INFINITY : -INFINITY; } else _p+=_n;
          }
        else 
          {
          if (INT_MIN/_p<10) { ylib_errno=YERROR_OK; return 0.;                             } else _p*=10;
          if (_p-INT_MIN<_n) { ylib_errno=YERROR_OK; return 0.;                             } else _p-=_n; 
          }
        }  
      }
    break;
    }
  else
    {
    if ( (_p)) { if (_p==INT_MIN) { ylib_errno=YERROR_OK; return 0.; } else _p--; } 
    switch (*str)
      {
      case '0' : { _n=0; break; }
      case '1' : { _n=1; break; }
      case '2' : { _n=2; break; }
      case '3' : { _n=3; break; }
      case '4' : { _n=4; break; }
      case '5' : { _n=5; break; }
      case '6' : { _n=6; break; }
      case '7' : { _n=7; break; }
      case '8' : { _n=8; break; }
      case '9' : { _n=9; break; }
      case '.' : { if (!(_p)) { _p=1; continue; } } 
      default  : { goto LABEL_RETURN_VALUE; } //Not a digit
      }            
    _d*=10., _d+=(double)_n;
    }
  }
LABEL_RETURN_VALUE: if (abs(_p)>8) _d*=powl(10.0,(double)_p); else { if (_p>0) { while (_p--) _d*=10.; } else { while (_p++) _d*=0.1; } }
ylib_errno=YERROR_OK; return (!(_s)) ? _d : -_d;
}


//All to upper case
void all_to_upper(register char *string)
{
register unsigned int _i=0;
while (string[_i]!='0')
  {
  if (islower((int)(string[_i])))
    string[_i]=(char)(toupper(string[_i]));
  _i++;
  }
}
void lex_to_upper(register char *lexem)
{
register unsigned int _i=0;
while ( (lexem[_i]!=' ')&&(lexem[_i]!='0')&&(lexem[_i]!='\n')&&(lexem[_i]!='\t') )
  {
  if (islower((int)(lexem[_i])))
    lexem[_i]=(char)(toupper(lexem[_i]));
  _i++;
  }
}

//All to lower case
void all_to_lower(register char *string)
{
register unsigned int _i=0;
while (string[_i]!='0')
  {
  if (isupper((int)(string[_i])))
    string[_i]=(char)(tolower(string[_i]));
  _i++;  
  }
}
void lex_to_lower(register char *lexem)
{
register unsigned int _i=0;
while ( (lexem[_i]!=' ')&&(lexem[_i]!='0')&&(lexem[_i]!='\n')&&(lexem[_i]!='\t') )
  {
  if (isupper((int)(lexem[_i])))
    lexem[_i]=(char)(tolower(lexem[_i]));
  _i++;
  }
}

// ----------------------- D U A L    S T R I N G -----------------------------

//This function add the lexem to end of string. It return the lenght of combined string.
unsigned int add_string(register char *A,register char *B)
{
register unsigned int _i,_j=0;

A=(char*)realloc(A,sizeof(char*)*((_i=strlen(A))+strlen(B)+1));
while(B[_j]!='\0')
  A[_i++]=B[_j++];
A[_i++]='\0';
return _i;
}

//This function subtrace string B from string A, it return the lenght of string A after operation
unsigned int sub_string(register char *A,register char *B)
{
register unsigned int _i,_j;

if ( ((_i=find_pattern(A,B))!=(unsigned int)-1)&&(_j=strlen(B)))
  {
  _i--;
  memcpy(&A[_i],&A[_i+_j],_j);
  return TRUE;
  }
return FALSE;
}

//This function looking for pattern in string
//It return the letter number of first pattern character in the sting or (unsigned int)-1 if the pattern was not found
unsigned int find_pattern(register char *pattern,register char *string)
{
register unsigned int _i,_j=0;

if (pattern[0]=='\0') return TRUE;
while (string[_j]!='\0')
  if (string[_j]==pattern[0])
    {
    _i=1;
    while (pattern[_i]!='\0')
      {
      if (string[_j+_i]!=pattern[_i]) goto NEXT_J;
      _i++;
      }
    return _j;
    }
  else
    NEXT_J: _j++;
return (unsigned int)-1;
}

//This function extend (from lexem) or compress (into lexem) line to given lenght
void equalize_string(register char* lexem,unsigned int lenght,register char *string)
{
register unsigned int _j,_i;

if ((_j=strlen(string))<lenght)
  do{
    _i=strlen(lexem);
    memcpy(&string[_j],lexem,(_i>lenght-_j) ? lenght-_j : _i);
    }while(_j+=_i<lenght);
else
  {
  if (lexem)
    memcpy(lexem,string,lenght-_j);
  lexem[_j]='\0';
  }
string[lenght]='\0';
}


//This function exchange two zones of memory. It's dangerous: due to speedup reasons the function doesn't check if the exchanged areas overlap.
inline void memexchange(register void *a, register void *b, register size_t size)
{
register int buff;

while (size>sizeof(int))
  {
  buff=*(int*)a, *(int*)a=*(int*)b, *(int*)b=buff;
  a+=sizeof(int), b+=sizeof(int), size-=sizeof(int);
  }
while (size)
  {
  buff=*(char*)a, *(char*)a=*(char*)b, *(char*)b=buff;
  a+=sizeof(char), b+=sizeof(char), size-=sizeof(char);
  }
}


// ------------------- T R I P L E    S T R I N G -----------------------

//This function exchange all patterns of pattern1 to pattern2
//It return the number of made changes (or (unsigned int)-1 on failure
unsigned int replace_pattern(char *pattern1,char *pattern2,char *string)
{
unsigned int _i,len0,len1,len2,count=0;
char *cp;

len0=strlen(string);
len0++;
len1=strlen(pattern1);
len2=strlen(pattern2);
while ((_i=find_pattern(pattern1,string))!=(unsigned int)-1)
  {
  if ( (cp=(char*)realloc(string,sizeof(char)*(len0+len2-len1)))) { ylib_errno=YERROR_MEMORY; return (unsigned int)-1; }
  else string=cp;
  //Shift string to patter2
  memcpy(&string[_i],&string[_i+len2-len1],len0-_i);
  //Insert pattern2
  memcpy(&string[_i-len1],pattern2,len2);
  count++;
  }
return count;
}

// ------------------------- Q U A D R O    S T I N G ---------------------------

//This function take value stored between first found pattern1 and first found pattern2
unsigned int find_lexem_inside_patterns(unsigned int *start,char *pattern1,char *pattern2,char *string)
{
register unsigned int end;

while ((*start=find_pattern(pattern1,string))!=(unsigned int)-1)
  while ((end=find_pattern(pattern2,string))!=(unsigned int)-1)
    if (*start<end)
      return *start-end-strlen(pattern2);
return FALSE;
}












