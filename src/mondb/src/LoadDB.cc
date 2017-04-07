#include <mysql.h>
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

bool ParseFileName(char* filename, string& name, string &dettype, 
		   string& detname, string& starttime, string& stoptime) {
  
  name = filename;

  dettype.assign(name,0,2);
  cout<<"DetType = "<<dettype<<" ";

  detname.assign(name,0,8);
  cout<<"DetName = "<<detname<<" ";
  
  int pos = name.find("~~start-",0,8); pos += 8;
  starttime.assign(name,pos,19);  
  pos=0;
  while((pos=starttime.find("-",pos,1)) != -1) 
    starttime.erase(pos,1);
  pos=0;
  while((pos=starttime.find(":",pos,1)) != -1) 
    starttime.erase(pos,1);
  cout<<"StartTime = "<<starttime<<" ";

  pos = name.find("~~finish-",0,9); pos += 9;
  stoptime.assign(name,pos,19);  
  pos=0;
  while((pos=stoptime.find("-",pos,1)) != -1) 
    stoptime.erase(pos,1);
  pos=0;
  while((pos=stoptime.find(":",pos,1)) != -1) 
    stoptime.erase(pos,1);
  cout<<"StopTime = "<<stoptime<<" ";  
  cout<<endl;  
  return true;
}

void PrintFile(char* filename) {
  
//    ifstream infile(filename);
 
//    if(infile) {
//      string s;
//      while(infile>>s) {cout<<s<<endl;}
  
//      infile.close();
//    }
//    else cerr<<"Cannot open file : "<<filename<<endl;

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


  MYSQL_RES *result;
  MYSQL_ROW row;
  MYSQL *connection, mysql;
  int k, state, Count, Match[5], Rows;
  //char cmd[120];
  char * Min, * Max;

  // connect to a MySQL database 
  cout<<"user "<<USER<<" connects to database "<<DBNAME<<" @ "<<HOST<<" ... "; 
  mysql_init(&mysql);
  connection = mysql_real_connect(&mysql, "localhost", "truc", "",
				  "CDB", 0, 0, 0);

  /* check for connection error */
  if (connection == NULL) {
    printf(mysql_error(&mysql));
    return 1;
  }
  else cout<<"OK"<<endl;

  // current working directory
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


  string cmd="select StartTime,DetType,DetName,FileName from data where DetType = \"FI\"";
  state = mysql_query(connection,cmd.c_str());
  if (state != 0) {
    printf(mysql_error(connection));
    return 1;
  }   
  result = mysql_store_result(connection);
  Rows = mysql_num_rows(result);
  
  row = mysql_fetch_row(result);
  while(row!=NULL) {
    printf("%s %s %s %s\n",row[0],row[1],row[2],row[3]);
    //PrintFile(row[3]);
    
    row = mysql_fetch_row(result);
  }

  mysql_free_result(result);
  mysql_close(connection);
}







