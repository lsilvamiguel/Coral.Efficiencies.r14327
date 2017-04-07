#include <wordexp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern "C" int expnm_( char* name, char* expname, long name_len, long
expname_len ) {

  int i;
  for( i=0; name[i] != ' ' && i<name_len; i++ ); name[i] = '\0';

  int status;
  wordexp_t* exp = (wordexp_t*) malloc( 512 );
  if( wordexp( name, exp, 0 ) == 0 ) {
    strcpy( expname, exp->we_wordv[0] );
    status = 0;
  }
  else {
    status = 1;
 }
  
  wordfree( exp );
  return status;

}
