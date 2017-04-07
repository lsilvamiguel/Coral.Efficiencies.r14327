#include "Parameters.h"
#include <string>
#include <ctype.h>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

MYSQL *connection;
MYSQL mysql;


int main(int argc,char **args) {
 
  if(argc!=5) {
    cout<<"usage : Add2DB <filename> <TBName> <start time> <stop time>"<<endl;
    cout<<"time format accepted :  yyyy-mm-dd-hh:mm:ss"<<endl;
    return 1;
  }


  string filename = args[1];
  string tbname = args[2];  
  string starttime = MakeSQLTime(string(args[3]));
  string stoptime = MakeSQLTime(string(args[4]));

  // checks that this file exists, and that it is ASCII
  ifstream in(filename.c_str());
  if(!in) {
    cerr<<filename<<" doesn't exist in this directory"<<endl;
    return 1;
  }
  else {
    unsigned char c;
    in.get(c);
    if(!isascii(c)) {
      cerr<<"data in "<<filename<<" is not ASCII"<<endl;
      return 1;
    }
  }
  in.close();

  // connect to a MySQL database /////////////////////////////////////////////

  if(!DBOpen(connection, &mysql, WRITE)) return 1;
  
  // find database location
  string dirname = GetLocation();
  if(dirname.empty()) {
    cerr<<"FATAL : Could not get location of physical calibration files."<<endl;
    exit(1);
  }
 
  string newname = MakeFileName(tbname, starttime, 
				stoptime);  
  if(newname.empty()) {
    cerr<<"FATAL : cannot generate filename"<<endl;
    exit(1);
  }

  // enters filename in database 
  int authorid = DBEntryAuthor();
  if(authorid) {
    if(!DBEntry(starttime, stoptime, string(tbname,0,2),tbname, newname, authorid)) {
      cerr<<"FATAL : Cannot make entry"<<endl;
      return 1;
    }
  }
  else {
    cerr<<"FATAL : Author not registered"<<endl;
    return 1;
  }

  // copy this file to the DB physical directory
  string cp = "cp ";
  cp += filename;
  cp += " ";
  cp += dirname +"/";
  cp += newname;

  cout<<cp.c_str()<<endl;
  if(system(cp.c_str())) {
    cerr<<"cannot write in "<<dirname<<endl;
    return 1;
  }  

  mysql_close(connection);
  return 0;
}



