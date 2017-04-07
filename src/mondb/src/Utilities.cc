#include <unistd.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include "Parameters.h"


bool DBOpen(MYSQL *conx, MYSQL* sql, int mode) {
  // Connection to the MySQL DataBase
 
  string user;
  string pwd;
  switch(mode) {
  case READ:
    user = READER;
    pwd="";
    cout<<"user "<<user<<" connects to database "<<DBNAME<<" @ "<<HOST<<" ... "; 
    break;
  case WRITE:
    user = WRITER;
    cout<<"user "<<user<<" connects to database "<<DBNAME<<" @ "<<HOST<<" ... "; 
    cout<<"password : ";
    cin>>pwd;
    break;
  case ADMIN:
    user = ADMINISTRATOR;
    cout<<"user "<<user<<" connects to database "<<DBNAME<<" @ "<<HOST<<" ... "; 
    cout<<"password : ";
    cin>>pwd;
    break;
  }
  mysql_init(sql);
  connection = mysql_real_connect(sql,HOST,user.c_str(),pwd.c_str(),
				  DBNAME, 0, 0, 0);
  if (connection == NULL) {
    cout<<mysql_error(sql)<<endl;
    return false;
  }
  else {
    cout<<"OK"<<endl;
    return true;
  }
}

bool DBEntry(const string& starttime, const string& stoptime, 
	     const string& dettype, const string& detname, const string& name,
	     int authorid) {
//    cout<<"ADDING "<<name<<endl;
//    cout<<DBNAME<<" "
//        <<HOST<<" "
//        <<USER<<" "
//        <<starttime<<" "
//        <<stoptime<<" "
//        <<dettype<<" "
//        <<detname<<endl;

  char aid[10];
  sprintf(aid,"%d",authorid);

  string cmd = "INSERT INTO tb_calibration (ID, StartTime, EndTime, DetType, DetName, FileName, authorID) VALUES ('','";
  cmd += starttime + "','";
  cmd += stoptime + "','";
  cmd += dettype + "','";
  cmd += detname + "','";
  cmd += name + "','";
  cmd += aid; cmd += "')";

  if(starttime.empty() || stoptime.empty()) {
    cout<<cmd<<" NOT EXECUTED !!!"<<endl;
    return false;
  }

  bool state = mysql_query(connection,cmd.c_str());
  if (state != 0) {
    //cout<<mysql_error(connection)<<endl;
    return false;
  }  

  return true;
}


int DBEntryAuthor() {

  string login = getlogin();

  string firstname;
  string lastname;
  string xwho = "phone ";
  xwho += login;
  FILE* scratch;    
  if( !(scratch = popen ( xwho.c_str(), "r" )) ){
    cerr<<"Cannot phone user."<<endl;
  }
  else {
    char ln[50];
    char fn[50];
    fscanf(scratch,"%s%s",ln,fn);
    firstname = fn;
    lastname = ln;
  }
  pclose( scratch );


  string email = firstname +"."+ lastname + "@cern.ch";
  
  // Looks in the DB if this author is registered : 

  string cmd = "select ID from tb_authors where Login='";
  cmd += login + "'";
  int state = mysql_query(connection, cmd.c_str());
  if (state != 0) {
    cout<<mysql_error(connection)<<endl;
    return 0;
  }  
  
  MYSQL_RES *result = mysql_store_result(connection);
  MYSQL_ROW row=mysql_fetch_row(result);
  if(row != NULL) { // this entry exists 
    return atoi(row[0]);
  } else { // adding this author to the list
    cmd = "INSERT INTO tb_authors (ID,login,LastName,FirstName,email) VALUES ('','";
    cmd += login + "','";
    cmd += lastname + "','";
    cmd += firstname + "','";
    cmd += email + "')";
    bool state = mysql_query(connection,cmd.c_str());
    if (state != 0) {
      cout<<mysql_error(connection)<<endl;
      return 0;
    }  
    mysql_free_result(result);
    return mysql_insert_id(&mysql);
  }
  return  0;
}



bool DecodeTime(vector<string>& results, const string& in) {
  switch(in.size()) {
  case 19:
    // Sasha or mysql DateTime
    results.push_back(string(in,0,4));
    results.push_back(string(in,5,2));
    results.push_back(string(in,8,2));
    results.push_back(string(in,11,2));
    results.push_back(string(in,14,2));
    results.push_back(string(in,17,2));
    break;
  default:
    cerr<<"DecodeTime : format not recognized"<<endl;
    cerr<<"use           yyyy-mm-dd-hh:mm:ss"<<endl;
    cerr<<in<<endl;
    return false;
  }
  return true;
}

string MakeTime(const string& in) {


  string out;
  string year; string month; string day; string hour; string min; string sec;

  int pos = in.find("-",1);
  if(pos != -1) { // Sasha's format ?
    vector<string> results;
    if(DecodeTime(results,in)) {
      
      year=results[0];
      month=results[1];
      day=results[2];
      hour=results[3];
      min=results[4];
      sec=results[5];
      
      if(atoi(hour.c_str())==24) {
	hour = "23";
	min = "59";
      }
	

      out += year; out += "-";
      out += month; out += "-";
      out += day; out += "-";
      out += hour; out += ":";
      out += min; out += ":";
      out += sec; 
    }
  }
  return out;
}

string MakeSQLTime(const string& in) {


  string out;
  string year; string month; string day; string hour; string min; string sec;

  vector<string> results;
  if(DecodeTime(results,in)) {
    
    year=results[0];
    month=results[1];
    day=results[2];
    hour=results[3];
    min=results[4];
    sec=results[5];
      
    out += year; out += "-";
    out += month; out += "-";
    out += day; out += " ";
    out += hour; out += ":";
    out += min; out += ":";
    out += sec; 
  }
  assert(!out.empty());
  return out;
}

const char* MakeFileName(const string& detname, const string& starttime, 
			 const string& stoptime) {

  string filename = detname;
  string start = MakeTime(starttime);
  string stop = MakeTime(stoptime);

  if(start.empty()||stop.empty())
    return "";

  filename += "~~start-" + start; 
  filename += "~~finish-" + stop; 
  
  return filename.c_str();
}

string GetLocation() {
  string dirname;
  string cmd = "select Directory from tb_directories where tb_directories.key = 'calibration'"; 
  int state = mysql_query(connection, cmd.c_str());
  if (state != 0) {
    cout<<mysql_error(connection)<<endl;
    return dirname;
  }  
  
  MYSQL_RES *result = mysql_store_result(connection);
  MYSQL_ROW row=mysql_fetch_row(result);
  dirname=row[0];

  mysql_free_result(result);
  return dirname;
}





