#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <exception>
#include <ctype.h>

#include <cstdlib>

#include "MyInterface.h"

using namespace std;

MyInterface::MyInterface (const char* server, const char* username,
			  const char* password, const char* dbname)
  : MySQLInterface (server, username, password, dbname)
  { }

MyInterface::~MyInterface() { }


//*************** some convenient functions *************************

bool MyInterface::dirWritable(const string& dirname) const {
  // not too elegant... testing if dirname is writable by
  // touching and removing a dummy file there

  string test = dirname; test += "/t";
  string touch = "touch "; touch += test;
  if(system(touch.c_str())) {
    cout<<dirname<<" not writable"<<endl;
    return false;
  }
  else {
    string rm = "rm "; rm += test;
    system(rm.c_str());
    return true;
  }
}


string MyInterface::cdbUpload(const string& start, const string& stop,
			      const string& tbname,
			      const string& file,
			      const string& dir,
			      const string& authorid) {

  // if dir is an empty string, current directory is used.

  string dbfile = fileName(tbname,start, stop);
  if(dbfile.empty()) return "Bad tbname or time format";
  if(!dirWritable(dir)) return "DB directory not writable";
  if(!checkTBName(tbname)) return "Bad tbname";
  if(!checkLocalFile(file)) return "Incorrect local file";
  if(!cdbCheckAuthorID(authorid)) return "author not registered";

//   string qry = "INSERT INTO runlb.tb_calibration (ID, StartTime, EndTime, DetType, DetName, FileName, dirID, authorID) VALUES ('','";
  string qry = "INSERT INTO tb_calibration (ID, StartTime, EndTime, DetType, DetName, FileName, dirID, authorID) VALUES ('','";
  qry += start + "','";
  qry += stop + "','";
  string dettype(tbname);
  qry += string(dettype,0,2);
  qry += "','";
  qry += tbname + "','";
  qry += dbfile + "','";
  qry += cdbDirID(dir);
  qry += "','";
  qry += authorid; qry += "')";

  if(query(qry)) {

    string destfile = dir; destfile += "/"; destfile += dbfile;
    if(file != destfile) {
      string cp = "cp ";
      cp += file;
      cp += " ";
      cp += destfile;

      if(system(cp.c_str())) {
	return "DISCREPANCY !!!";
      }
    }

    endQuery();
    return "";
  }
  else
    return "Bad query";
}


string MyInterface::cdbUpload(const string& filename,
			      const string& dir,
			      const string& authorid) {

  string localdir = filename;
  int sfn = localdir.rfind("/");

  string name;
  name = filename.substr(sfn+1);
  if(name.size()<25) return name + " too short : \nminimal name is xxxxxxxx~~start-~~finish-";

  string dettype; dettype.assign(name,0,2);
  int startpoint = name.find("~~start-",0,8);
  string detname; detname.assign(name,0,startpoint);
  startpoint += 8;
  int stoppoint = name.find("~~finish-",0,9);
  string starttime; starttime.assign(name,startpoint,stoppoint-startpoint);
  string stoptime;
  stoptime = name.substr(stoppoint+9);

  return cdbUpload(starttime, stoptime, detname, filename, dir, authorid);
}


string MyInterface::cdbDirID(const string& dir) {
//   string qry = "select ID from runlb.tb_directories where runlb.tb_directories.Directory = '"; qry += dir; qry += "'";
  string qry = "select ID from tb_directories where tb_directories.Directory = '"; qry += dir; qry += "'";
  if(query(qry)) {
    if(getNRows())
      return getCol(0);
    endQuery();
  }

  cerr<<dir<<" not found"<<endl;
  return "";
}


string MyInterface::fileName(const string& tbname,
			     const string& start,
			     const string& stop) {
  string filename;
  if(!tbname.size()==8 ||
     !checkFormat(start) ||
     !checkFormat(stop))
    return filename;

  filename = tbname; filename += "~~start-";
  filename += start ; filename += "~~finish-";
  filename += stop ;
  return filename;
}


bool MyInterface::checkFormat(const string& time) {

  if(time.size() != MyInterface::fCdbTimeFormat.size()) {
    cerr<<time<<" has a bad format. use "<<MyInterface::fCdbTimeFormat<<endl;
    return false;
  }
  string field;
  for(unsigned i = 0; i<time.size(); i++) {
    char cref = fCdbTimeFormat[i];
    char c = time[i];
    if(!isalnum(cref)) {
      // separator. must be the same as in timeformat
      if(isalnum(c)) {
	cerr<<time<<" has a bad format. use "<<MyInterface::fCdbTimeFormat<<endl;
	return false;
      }
      field = "";
    }
    else field += c;
  }
  return true;
}


bool MyInterface::checkTBName(const string& tbname) {
  if(tbname.size() == 8) return true;
  cerr<<tbname<<" is a bad TBName. must have 8 letters"<<endl;
  return false;
}


bool MyInterface::checkLocalFile(const string& filename) {
  // checks that this file exists. some more checks to be added
  ifstream in(filename.c_str());
  if(!in) {
    cerr<<filename<<" doesn't exist"<<endl;
    return false;
  }
  return true;
}


string MyInterface::authorID(const string& login,
			     const string& firstname,
			     const string& lastname,
			     const string& email) {

//   string qry = "select ID from runlb.tb_authors where Login = '";
  string qry = "select ID from tb_authors where Login = '";
  qry += login; qry += "' and FirstName = '";
  qry += firstname; qry += "' and LastName = '";
  qry += lastname; qry += "' and email = '";
  qry += email; qry += "'";

  if(query(qry)) {
    if(getNRows()) {
      endQuery();
      return getCol(0);
    }
    else {

//       qry = "insert into runlb.tb_authors (ID,Login,FirstName,Lastname,email) values ('','";
      qry = "insert into tb_authors (ID,Login,FirstName,Lastname,email) values ('','";
      qry += login; qry += "','";
      qry += firstname; qry += "','";
      qry += lastname; qry += "','";
      qry += email; qry += "')";

      query(qry);
      endQuery();

      char id[10];
      sprintf(id,"%d",static_cast<unsigned int>(mysql_insert_id(&fMysql)));
      return id;
    }
  }
  else return "";
}


bool MyInterface::cdbCheckAuthorID(const string& id) {
//   string qry = "select * from runlb.tb_authors where ID = '";
  string qry = "select * from tb_authors where ID = '";
  qry += id; qry += "'";

  if(query(qry)) {
    if(getNRows()) {
      endQuery();
      return true;
    }
    else
      return false;
  }
  cerr<<"Bad query : "<<qry<<endl;
  return false;
}


string MyInterface::fileOwner(const string& start, const string& stop,
			    const string& file, const string& authorid) {

//   string qry = "select authorID, LastName from runlb.tb_calibration,runlb.tb_authors where StartTime = '";
  string qry = "select authorID, LastName from tb_calibration,tb_authors where StartTime = '";
  qry += start; qry += "' and EndTime = '";
  qry += stop; qry += "' and FileName = '";
  qry += file; qry += "' and tb_authors.ID = tb_calibration.authorID";

  string owner = "";
  if(query(qry)) {
    if(getNRows()) {
      if(authorid == getCol(0))
	owner = "myself";
      else {
	owner = getCol(1);
      }
    }
    else
      return "";
    endQuery();
    return owner;
  }
  cerr<<"Bad query : "<<qry<<endl;
  return "";
}
