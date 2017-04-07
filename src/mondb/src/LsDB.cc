#include "Parameters.h"
#include <string>
#include <stdlib.h>
#include <unistd.h>

MYSQL *connection;
MYSQL mysql;

int main(int argc,char **args) {
  
  // Connection to the MySQL DataBase
 
  if(!DBOpen(connection, &mysql,READ)) return 1;

  // find database location 
  string dirname = GetLocation();

  //cd DB directory
  char curdir[100];
  getcwd(curdir,100);
  if(chdir(dirname.c_str()))
    return 1;

  // list directory
  string ls = "ls ";
  if(argc == 2) {
    ls += string(args[1]);
  }
  if(system(ls.c_str())) {
    return 1;
  } 

  // come back
  if(chdir(curdir))
    return 1;  
  
  return 0;
}




