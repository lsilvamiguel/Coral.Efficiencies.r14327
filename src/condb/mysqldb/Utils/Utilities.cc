#include <unistd.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include "Parameters.h"


bool DBOpen(MYSQL *conx, MYSQL* sql, int mode) {
  // Connection to the MySQL DataBase
 
  std::string user;
  std::string pwd;
  switch(mode) {
  case READ:
    user = READER;
    pwd="";
    std::cout<<"user "<<user<<" connects to database "<<DBNAME<<" @ "<<HOST<<" ... "; 
    break;
  case WRITE:
    user = WRITER;
    std::cout<<"user "<<user<<" connects to database "<<DBNAME<<" @ "<<HOST<<" ... "; 
    std::cout<<"password : ";
    std::cin>>pwd;
    break;
  case ADMIN:
    user = ADMINISTRATOR;
    std::cout<<"user "<<user<<" connects to database "<<DBNAME<<" @ "<<HOST<<" ... "; 
    std::cout<<"password : ";
    std::cin>>pwd;
    break;
  }
  mysql_init(sql);
  connection = mysql_real_connect(sql,HOST,user.c_str(),pwd.c_str(),
				  DBNAME, 0, 0, 0);
  if (connection == NULL) {
    std::cout<<mysql_error(sql)<<std::endl;
    return false;
  }
  else {
    std::cout<<"OK"<<std::endl;
    return true;
  }
}

bool DBEntry(const std::string& starttime, const std::string& stoptime, 
	     const std::string& dettype, const std::string& detname, const std::string& name,
	     int authorid) {
//    std::cout<<"ADDING "<<name<<std::endl;
//    std::cout<<DBNAME<<" "
//        <<HOST<<" "
//        <<USER<<" "
//        <<starttime<<" "
//        <<stoptime<<" "
//        <<dettype<<" "
//        <<detname<<std::endl;

  char aid[10];
  sprintf(aid,"%d",authorid);

  std::string cmd = "INSERT INTO tb_calibration (ID, StartTime, EndTime, DetType, DetName, FileName, authorID) VALUES ('','";
  cmd += starttime + "','";
  cmd += stoptime + "','";
  cmd += dettype + "','";
  cmd += detname + "','";
  cmd += name + "','";
  cmd += aid; cmd += "')";

  if(starttime.empty() || stoptime.empty()) {
    std::cout<<cmd<<" NOT EXECUTED !!!"<<std::endl;
    return false;
  }

  bool state = mysql_query(connection,cmd.c_str());
  if (state != 0) {
    //std::cout<<mysql_error(connection)<<std::endl;
    return false;
  }  

  return true;
}


int DBEntryAuthor() {

  std::string login = getlogin();

  std::string firstname;
  std::string lastname;
  std::string xwho = "phone ";
  xwho += login;
  FILE* scratch;    
  if( !(scratch = popen ( xwho.c_str(), "r" )) ){
    std::cerr<<"Cannot phone user."<<std::endl;
  }
  else {
    char ln[50];
    char fn[50];
    fscanf(scratch,"%s%s",ln,fn);
    firstname = fn;
    lastname = ln;
  }
  pclose( scratch );


  std::string email = firstname +"."+ lastname + "@cern.ch";
  
  // Looks in the DB if this author is registered : 

  std::string cmd = "select ID from tb_authors where Login='";
  cmd += login + "'";
  int state = mysql_query(connection, cmd.c_str());
  if (state != 0) {
    std::cout<<mysql_error(connection)<<std::endl;
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
      std::cout<<mysql_error(connection)<<std::endl;
      return 0;
    }  
    mysql_free_result(result);
    return mysql_insert_id(&mysql);
  }
  return  0;
}



bool DecodeTime(std::vector<std::string>& results, const std::string& in) {
  switch(in.size()) {
  case 19:
    // Sasha or mysql DateTime
    results.push_back(std::string(in,0,4));
    results.push_back(std::string(in,5,2));
    results.push_back(std::string(in,8,2));
    results.push_back(std::string(in,11,2));
    results.push_back(std::string(in,14,2));
    results.push_back(std::string(in,17,2));
    break;
  default:
    std::cerr<<"DecodeTime : format not recognized"<<std::endl;
    std::cerr<<"use           yyyy-mm-dd-hh:mm:ss"<<std::endl;
    std::cerr<<in<<std::endl;
    return false;
  }
  return true;
}

std::string MakeTime(const std::string& in) {


  std::string out;
  std::string year; std::string month; std::string day; std::string hour; std::string min; std::string sec;

  int pos = in.find("-",1);
  if(pos != -1) { // Sasha's format ?
    std::vector<std::string> results;
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

std::string MakeSQLTime(const std::string& in) {


  std::string out;
  std::string year; std::string month; std::string day; std::string hour; std::string min; std::string sec;

  std::vector<std::string> results;
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

const char* MakeFileName(const std::string& detname, const std::string& starttime, 
			 const std::string& stoptime) {

  std::string filename = detname;
  std::string start = MakeTime(starttime);
  std::string stop = MakeTime(stoptime);

  if(start.empty()||stop.empty())
    return "";

  filename += "~~start-" + start; 
  filename += "~~finish-" + stop; 
  
  return filename.c_str();
}

std::string GetLocation() {
  std::string dirname;
  std::string cmd = "select Directory from tb_directories where tb_directories.key = 'calibration'"; 
  int state = mysql_query(connection, cmd.c_str());
  if (state != 0) {
    std::cout<<mysql_error(connection)<<std::endl;
    return dirname;
  }  
  
  MYSQL_RES *result = mysql_store_result(connection);
  MYSQL_ROW row=mysql_fetch_row(result);
  dirname=row[0];

  mysql_free_result(result);
  return dirname;
}





