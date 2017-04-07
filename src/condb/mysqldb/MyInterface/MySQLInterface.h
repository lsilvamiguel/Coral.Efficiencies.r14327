#ifndef __MySQLInterface__
#define __MySQLInterface__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <typeinfo>
#include <mysql.h>
#include <time.h>


class MySQLInterface {

 protected:
  
  // connection parameters
  std::string    fServer;
  std::string    fUserName;
  std::string    fPassword;
  std::string    fDBName;
  int            fNumPort;

  // mysql variables
  MYSQL           fMysql;
  MYSQL          *fConnection;

  // true if connected
  bool            fConnected;

  // mysql results
  MYSQL_RES      *fResult;

  // true if some results are in memory
  bool            fResInMem;
  
  // mysql row
  MYSQL_ROW       fRow;

  // number of rows returned by the query
  int             fNRows;
  
  // time format used in filenames
  static const std::string fCdbTimeFormat;
  static const char* fCdbTimeFormatStr;
  static const char* fMySQLTimeFormatStr;

 public:
  
  /// Destructor
  virtual  ~MySQLInterface (void);

  MySQLInterface (const char* server, const char* username,
	          const char* passwd, const char* dbname);
 

  // interface to mysql -------------------------------------------
  bool   connect(void);
  void   disconnect(void);
  bool   isConnected(void)   { return fConnected; }
  int    selectDB (void)
           { return mysql_select_db(&fMysql, fDBName.c_str()); }
  bool   query(const std::string& request);
  char*  getCol(int i)  { if (fRow) return fRow[i]; else return 0; }
  void   endQuery(void);
  bool   getNextRow(void);
  int    getNRows() {return fNRows;}
  long   getLastInsertID() { return mysql_insert_id(&fMysql); }
  void   setNumPort(int nport)  { fNumPort = nport; }
  int    getNumPort(void)  { return fNumPort; }

 public:
  static std::string toMySQLtime (const char* cdbtime);
  static std::string toMySQLtime (tm* stm);
  static std::string toMySQLtime (time_t *timep);
  static std::string toCDBtime   (const char* mysqltime);


};

#endif 


