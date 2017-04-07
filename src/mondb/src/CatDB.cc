#include "Parameters.h"
#include <string>
#include <fstream>

MYSQL *connection;
MYSQL mysql;

int main(int argc,char **args) {

  if(argc!=2) {
    cout<<"usage : CatDB <filename>"<<endl;
    return 1;
  }
  
  if(!DBOpen(connection, &mysql,READ)) return 1;
   
  string dirname = GetLocation();
   
  string filename = dirname;
  filename += "/";
  filename += args[1];
  ifstream in(filename.c_str());
  if(!in) 
    cerr<<dirname<<" doesn't exist in "<<dirname<<endl;
  
  while(in) {
    char s[200];
    in.getline(s,200);
    cout<<s<<endl;
  } 
  return 0;
}




