#include <string>
#include <set>
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

bool ParseFileName(char* filename) {
  
  string name = filename;

  string dettype; dettype.assign(name,0,2);
  //cout<<"DetType = "<<dettype<<" ";

  int startpoint = name.find("~~start-",0,8); 
  string detname; detname.assign(name,0,startpoint);  
  startpoint += 8;
  int stoppoint = name.find("~~finish-",0,9); 
  string starttime; starttime.assign(name,startpoint,stoppoint-startpoint);  
  starttime = MakeTime(starttime);
  string stoptime; stoptime.assign(name,stoppoint+9);  
  stoptime = MakeTime(stoptime);
  
  if(starttime.empty()) {
    cerr<<"ParseFileName() bad start time - entry skipped"<<endl;
    return false;
  }
  if(stoptime.empty()) {
    cerr<<"ParseFileName() bad stop time - entry skipped"<<endl;
    return false;
  }

  if(int authorid = DBEntryAuthor()) { 
    if(DBEntry(starttime, stoptime, dettype, detname, name, authorid)) {
      //cerr<<"adding "<<filename<<" to database"<<endl;
      return true;  
    }
  }
  else {
    cerr<<"WARNING : Author unregistered"<<endl;
  }
  return false;
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


int main(int argc,char **argv) {

  bool doit = true;


  char cdirname[200]; getcwd(cdirname,200);
  string dirname(cdirname);

  int pargc = 1;
  while(pargc < argc) {
    char *pargv, *pargvn;
    pargv =  argv[pargc];
    pargvn = 0;
    
    if ((pargc + 1) < argc) pargvn =  argv[pargc+1];

    if (pargv[0] == '-') {   // it is an option
      if (strcmp(pargv, "-n") == 0) {  
	cerr<<"TEST MODE ACTIVATED. No modifications of physical file planned"<<endl;
	doit=false;
	continue;
      }      
      if (strcmp(pargv, "-dir") == 0) {  // specify directory
	if (pargvn && (pargvn[0] != '-')) dirname = pargvn;
	else cerr<<"badly formed -dir option"<<endl;
	pargc+=2;
	continue;
      }
      if (strcmp(pargv, "--help") == 0) {  
	cerr<<"Usage SyncDB [-n]"<<endl;
	exit(0);
      }      
    }
  }

  cout<<dirname<<endl;
  exit(0);

  // connect to a MySQL database ////////////////////////////////////////////

  if(!DBOpen(connection, &mysql,ADMIN)) {
    cerr<<"FATAL : cannot open database"<<endl;
    return 1;
  }

  // find database location 
  // string dirname = GetLocation();

  //cd DB directory
  char curdir[100];
  getcwd(curdir,1000);  
  if(chdir(dirname.c_str())) {
    cerr<<"FATAL : cannot go to directory "<<dirname<<endl;    
    return 1;  
  }
  
  //scan directory
  set<string> files;
  struct stat info;
  struct dirent *filedir;
  DIR* dirfd;
  cout<< "Scanning calibration files in "<<dirname<<endl;
  dirfd = opendir(dirname.c_str());
  while ((filedir = readdir(dirfd))) {
    if(filedir != NULL) {
      if( 0!=stat(filedir->d_name,&info) ) {
	cerr << "Can't stat file "<<filedir->d_name<<endl;
	return 1;
      }
      if (S_ISREG(info.st_mode)) {
	//cout << "reading "<<filedir->d_name<<" ... "<<endl;
	files.insert(string(filedir->d_name));
	ParseFileName(filedir->d_name);
      }
    }
  }
  closedir(dirfd);
  
  if(!doit) { 
    // come back
    if(chdir(curdir))
      return 1;  
    mysql_close(connection);
    return 0;
  }

  // scan DB and look for unlinked files or updates /////////////////////////
  cout<< "Scanning DB files for unlinked files and updates"<<endl;
  dirname+="/";

  string cmd = "select ID,StartTime,EndTime,DetType,DetName,FileName from tb_calibration where 1"; 
  int state = mysql_query(connection,cmd.c_str());
  if (state != 0) {
    cout<<mysql_error(connection)<<endl;
    return 1;
  }  
  MYSQL_RES *result = mysql_store_result(connection);
  MYSQL_ROW row=mysql_fetch_row(result);
  while(row != NULL) {
    if(files.find(row[5]) == files.end()) { 
      // the file has been removed -> the entry must disappear

      cout<<row[5]<<" UNLINKED ! ";

      string del = "delete from tb_calibration where ID=";
      del += row[0];
      bool state = mysql_query(connection,del.c_str());
      if (state != 0) {
	cout<<mysql_error(connection)<<endl;
	cout<<" cannot delete -> diy !"<<endl;
      } else {
	cout<<" -> entry deleted"<<endl;
      } 
    }
    else {
      string newname=MakeFileName(string(row[4]), string(row[1]), string(row[2]));
      if(newname != row[5]) {
	// something has been changed in the DB, changing the filename.
	cout<<"changing filename "<<row[5]<<" to "<<newname<<endl;
	string mv = "mv ";
	mv += dirname + row[5];
	mv += " ";
	mv += dirname + newname;
	if(system(mv.c_str())) {
	  cerr<<"cannot write in "<<dirname<<endl;
	  return 1;
	}	
	
	string up = "update tb_calibration set ID='"; up += row[0]; up+="',";
	up+="StartTime='"; up += row[1]; up+="',";
	up+="EndTime='"; up += row[2]; up+="',";
	up+="DetType='"; up += row[3]; up+="',";
	up+="DetName='"; up += row[4]; up+="',";
	up+="FileName='"; up += newname; up+="',";
	up+="authorID = '1' WHERE  ID = '"; up += row[0]; up+="'";
	bool state = mysql_query(connection,up.c_str());
	if (state != 0) {
	  cout<<mysql_error(connection)<<endl;
	  cout<<" cannot update -> diy !"<<endl;
	} 
      }
    }
    row=mysql_fetch_row(result);
  }
  mysql_free_result(result);  

  // come back
  if(chdir(curdir))
    return 1;  
  mysql_close(connection);

  return 0;
}








