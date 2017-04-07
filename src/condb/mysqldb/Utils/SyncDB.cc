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
#include "MyInterface.h"

MYSQL *connection;
MYSQL mysql;

bool ParseFileName(char* filename) {
  
  std::string name = filename;

  std::string dettype; dettype.assign(name,0,2);
  //std::cout<<"DetType = "<<dettype<<" ";

  int startpoint = name.find("~~start-",0,8); 
  std::string detname; detname.assign(name,0,startpoint);  
  startpoint += 8;
  int stoppoint = name.find("~~finish-",0,9); 
  std::string starttime; starttime.assign(name,startpoint,stoppoint-startpoint);  
  starttime = MakeTime(starttime);
  std::string stoptime; stoptime.assign(name,stoppoint+9,std::string::npos);  
  stoptime = MakeTime(stoptime);
  
  if(starttime.empty()) {
    std::cerr<<"ParseFileName() bad start time - entry skipped"<<std::endl;
    return false;
  }
  if(stoptime.empty()) {
    std::cerr<<"ParseFileName() bad stop time - entry skipped"<<std::endl;
    return false;
  }

  if(int authorid = DBEntryAuthor()) { 
    if(DBEntry(starttime, stoptime, dettype, detname, name, authorid)) {
      //std::cerr<<"adding "<<filename<<" to database"<<std::endl;
      return true;  
    }
  }
  else {
    std::cerr<<"WARNING : Author unregistered"<<std::endl;
  }
  return false;
}

void PrintFile(char* filename) {
  

  FILE *file=fopen(filename,"r");
  
  if(file!=NULL) {
    float value=0;
    std::cout<<std::endl;
    int n=0;
    while(fscanf(file,"%f",&value) != EOF) {
      if((n%8)==0) std::cout<<std::endl;
      std::cout<<value<<"  ";
      n++;
    }
    fclose(file);
    std::cout<<std::endl<<std::endl;
  }
  else {
    perror("erreur ");
    std::cout<<filename<<std::endl;
  }
}


int main(int argc,char **argv) {

  bool doit = true;

  char cdirname[200]; getcwd(cdirname,200);
  std::string dirname(cdirname);

  int pargc = 1;
  while(pargc < argc) {
    char *pargv, *pargvn;
    pargv =  argv[pargc];
    pargvn = 0;
    
    if ((pargc + 1) < argc) pargvn =  argv[pargc+1];

    if (pargv[0] == '-') {   // it is an option
      if (strcmp(pargv, "-n") == 0) {  
	std::cerr<<pargc<<" "<<argc<<"TEST MODE ACTIVATED. No modifications of physical file planned"<<std::endl;
	doit=false;
	pargc++;
	continue;
      }      
      if (strcmp(pargv, "-dir") == 0) {  // specify directory
	if (pargvn && (pargvn[0] != '-')) dirname = pargvn;
	else std::cerr<<"badly formed -dir option"<<std::endl;
	pargc+=2;
	continue;
      }
      if (strcmp(pargv, "--help") == 0) {  
	std::cerr<<"Usage SyncDB [-n]"<<std::endl;
	exit(0);
      }      
    }
  }


  // connect to a MySQL database ////////////////////////////////////////////
  
  MyInterface *db;
  try {
    db = new MyInterface(HOST,WRITER,PASSWORD,DBNAME);
  }
  catch(const char* errmsg) {
    std::cout<<"exception "<<errmsg<<std::endl;
    exit(1);
  }
  
  if(!db->connect()) exit(1);
  
  if(!db->dirWritable(dirname)) return 1;

  //cd and scan DB directory
  char curdir[100];
  getcwd(curdir,1000);  
  if(chdir(dirname.c_str())) {
    std::cerr<<"FATAL : cannot go to directory "<<dirname<<std::endl;    
    return 1;  
  }

  std::set<std::string> files;
  struct stat info;
  struct dirent *filedir;
  DIR* dirfd;
  std::cout<< "Scanning calibration files in "<<dirname<<std::endl;
  dirfd = opendir(dirname.c_str());
  while ((filedir = readdir(dirfd))) {
    if(filedir != NULL) {
      if( 0!=stat(filedir->d_name,&info) ) {
	std::cerr << "Can't stat file "<<filedir->d_name<<std::endl;
	return 1;
      }
      if (S_ISREG(info.st_mode)) {
	std::cout << "reading "<<filedir->d_name<<" ... "<<std::endl;
	files.insert(std::string(filedir->d_name));
	std::string errmsg = db->cdbUpload(filedir->d_name,dirname,"5");
	std::cout<<errmsg<<std::endl;
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

  exit(0);

  // scan DB and look for unlinked files or updates /////////////////////////
  std::cout<< "Scanning DB files for unlinked files and updates"<<std::endl;
  dirname+="/";

  std::string cmd = "select ID,StartTime,EndTime,DetType,DetName,FileName from tb_calibration where 1"; 
  int state = mysql_query(connection,cmd.c_str());
  if (state != 0) {
    std::cout<<mysql_error(connection)<<std::endl;
    return 1;
  }  
  MYSQL_RES *result = mysql_store_result(connection);
  MYSQL_ROW row=mysql_fetch_row(result);
  while(row != NULL) {
    if(files.find(row[5]) == files.end()) { 
      // the file has been removed -> the entry must disappear

      std::cout<<row[5]<<" UNLINKED ! ";

      std::string del = "delete from tb_calibration where ID=";
      del += row[0];
      bool state = mysql_query(connection,del.c_str());
      if (state != 0) {
	std::cout<<mysql_error(connection)<<std::endl;
	std::cout<<" cannot delete -> diy !"<<std::endl;
      } else {
	std::cout<<" -> entry deleted"<<std::endl;
      } 
    }
    else {
      std::string newname=MakeFileName(std::string(row[4]), std::string(row[1]), std::string(row[2]));
      if(newname != row[5]) {
	// something has been changed in the DB, changing the filename.
	std::cout<<"changing filename "<<row[5]<<" to "<<newname<<std::endl;
	std::string mv = "mv ";
	mv += dirname + row[5];
	mv += " ";
	mv += dirname + newname;
	if(system(mv.c_str())) {
	  std::cerr<<"cannot write in "<<dirname<<std::endl;
	  return 1;
	}	
	
	std::string up = "update tb_calibration set ID='"; up += row[0]; up+="',";
	up+="StartTime='"; up += row[1]; up+="',";
	up+="EndTime='"; up += row[2]; up+="',";
	up+="DetType='"; up += row[3]; up+="',";
	up+="DetName='"; up += row[4]; up+="',";
	up+="FileName='"; up += newname; up+="',";
	up+="authorID = '1' WHERE  ID = '"; up += row[0]; up+="'";
	bool state = mysql_query(connection,up.c_str());
	if (state != 0) {
	  std::cout<<mysql_error(connection)<<std::endl;
	  std::cout<<" cannot update -> diy !"<<std::endl;
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








