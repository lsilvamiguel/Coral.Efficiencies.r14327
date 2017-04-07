#ifndef __CS_MIGR_ORA_SESSION_H__
#define __CS_MIGR_ORA_SESSION_H__

#include <vector>
#include <string>
#include <occi.h>

class CsOraSession
{

  static CsOraSession* _instance;
 public:
  static CsOraSession* instance();
  
  bool init(int runNumber, bool isDstProduction);
  bool initInfoSession();
  bool terminate();
  bool commit();
  bool rollback();

  oracle::occi::Statement*    stmt()          { return _runStatus?_stmt:NULL; }
  oracle::occi::Connection*   conn()          { return _runStatus?_conn:NULL; }
  const std::string&	getLastSqlStr() { return _sqlStr; }

  // execute a query on the default statement
  oracle::occi::ResultSet*	executeQuery(const std::string& sql);
  oracle::occi::ResultSet*	executeQuery(const std::string& sql, const int val);
  bool                      executeUpdate(const std::string& sql);

  virtual 	~CsOraSession();

  bool          stopSession();
  bool          startSession();

  std::string decodePassword(std::string encodedPassword);

private:
  CsOraSession() : _env(NULL), _conn(NULL), _stmt(NULL) {};
  
protected:  

  oracle::occi::Environment*		_env;
  oracle::occi::Connection*		_conn;
  oracle::occi::Statement*		_stmt;
  short			_initStatus;
  short                 _runStatus;

  std::string		_sqlStr;

private:
   static std::string _svcServiceName;
	std::string	_svcUser;
	std::string _svcPwd;
};

#endif
