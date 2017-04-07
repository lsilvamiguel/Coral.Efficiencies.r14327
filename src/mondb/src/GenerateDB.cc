#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include "errno.h"
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "Parameters.h"

MYSQL *connection;
MYSQL mysql;


extern bool DBEntry(const string& starttime, const string& stoptime, 
		    const string& dettype, const string& detname, 
		    const string& name);

bool ParseFileName(char* filename) {
  
  string name = filename;

  string dettype; dettype.assign(name,0,2);
  //cout<<"DetType = "<<dettype<<" ";

  string detname; detname.assign(name,0,8);
  //cout<<"DetName = "<<detname<<" ";
  
  int pos = name.find("~~start-",0,8); pos += 8;
  string starttime; starttime.assign(name,pos,19);  
  pos=0;
  while((pos=starttime.find("-",pos,1)) != -1) 
    starttime.erase(pos,1);
  pos=0;
  while((pos=starttime.find(":",pos,1)) != -1) 
    starttime.erase(pos,1);
  //cout<<"StartTime = "<<starttime<<" ";

  pos = name.find("~~finish-",0,9); pos += 9;
  string stoptime; stoptime.assign(name,pos,19);  
  pos=0;
  while((pos=stoptime.find("-",pos,1)) != -1) 
    stoptime.erase(pos,1);
  pos=0;
  while((pos=stoptime.find(":",pos,1)) != -1) 
    stoptime.erase(pos,1);
  //cout<<"StopTime = "<<stoptime<<" ";  
  //cout<<endl; 
 
  return DBEntry(starttime, stoptime, dettype, detname, name);
}

void PrintFile(char* filename) {
  

  FILE *file=fopen(filename,"r");
  
  if(file!=NULL) {
    float value=0;
    cout<<endl;
    int n=0;
    while(fscanf(file,"%f",&value) != EOF) {
      if((n%8)==0) cout<<endl;
      cout.form("%5.3f  ",value);
      n++;
    }
    fclose(file);
    cout<<endl<<endl;
  }
  else {
    perror("erreur ");
    cout<<filename<<endl;
  }
}


int main(int argc,char **args) {

  MYSQL_ROW row;

  int k, state, Count, Match[5], Rows;
  //char cmd[120];
  char * Min, * Max;

  // connect to a MySQL database /////////////////////////////////////////////
  cout<<"user "<<USER<<" connects to database "<<DBNAME<<" @ "<<HOST<<" ... "; 
  mysql_init(&mysql);
  connection = mysql_real_connect(&mysql,HOST,USER , "daqna58",
				  DBNAME, 0, 0, 0);

  if (connection == NULL) {
    printf(mysql_error(&mysql));
    return 1;
  }
  else cout<<"OK"<<endl;

  // user ID /////////////////////////////////////////////////////////////////  
  

  // current working directory ///////////////////////////////////////////////
  char cwd[200];
  if (getcwd(cwd,200)== NULL) {
    printf("cannot get cwd name. Increase buffer.");
  }
  else cout<<"Current working directory ... "<<cwd<<endl;

  // scan cwd 
  struct stat info;
  struct dirent *filedir;
  string filedirname;
  DIR* dirfd;
  cout<< "Scanning calibration files in "<<cwd<<endl;
  dirfd = opendir(cwd);
  while ((filedir = readdir(dirfd))) {
    if(filedir != NULL) {
      if( 0!=stat(filedir->d_name,&info) ) {
	cerr << "Can't stat file "<<filedir->d_name<<endl;
	return 1;
      }
      if (S_ISREG(info.st_mode)) {
	cout << "reading "<<filedir->d_name<<" ... "<<endl;
	if(ParseFileName(filedir->d_name)) {
	  cout<<"                                   ... OK"<<endl;
	}
      }
    }
  }
  closedir(dirfd);

  mysql_close(connection);
}







