#ifndef __MyInterface__
#define __MyInterface__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <typeinfo>
#include <mysql.h>


class MyInterface {

 private:
  
  string    fServer;
  string    fUserName;
  string    fPassword;
  string    fDBName;
  MYSQL           fMysql;
  MYSQL          *fConnection;
  bool            fConnected;
  MYSQL_RES      *fResult;
  bool            fResInMem;
  MYSQL_ROW       fRow;
  int             fNRows;
  
  static const string fCdbTimeFormat;

 public:
  
  /// Destructor
  virtual  ~MyInterface (void);

  MyInterface (const char* server, const char* username,
	       const char* passwd, const char* dbname);
 

  // interface to mysql 
  bool   connect(void);
  void   disconnect(void);
  bool   isConnected(void)   { return fConnected; }
  bool   query(string& request);
  char*  getCol(int i)  { if (fRow) return fRow[i]; else return 0; }
  void   endQuery(void);
  bool   getNextRow(void);
  int    getNRows() {return fNRows;}

  // some tools
  bool   dirWritable(const string& dirname) const; 
  
  string fileName(const string& tbname,const string& start, 
		  const string& stop);
  
  // convenient functions for cdb
  
  string   cdbUpload(const string& start, const string& stop, 
		     const string& tbname,
		     const string& file, const string& dir, 
		     const string& authorid);

  string   cdbUpload(const string& file,
		     const string& authorid);

  string cdbDirID(const string& dir);

  string authorID(const string& login,
		  const string& firstname,
		  const string& lastname,
		  const string& email);
  
  bool checkFormat(const string& time);

  bool checkTBName(const string& tbname);

  bool checkLocalFile(const string& filename);

  bool cdbCheckAuthorID(const string& id);

  
};

#endif 


