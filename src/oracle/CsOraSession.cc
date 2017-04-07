#include <iostream>
#include <vector>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>


#include "CsOraSession.h"
#include "RunInfo.h"

#include "CsOpt.h"

using namespace oracle::occi; 
using namespace std; 

CsOraSession* CsOraSession::_instance = NULL;
string CsOraSession::_svcServiceName("");

CsOraSession* CsOraSession::instance()
{
 
	if(!_instance) {
		_instance = new CsOraSession;
		_svcServiceName = ""; 
	}
 return _instance;
}

bool CsOraSession::initInfoSession()
{
  _initStatus = 0;
  _runStatus  = 0;

  if ( _svcServiceName == "" ) { 
	  bool bDBParametersGotFromOptFile(false);
	  CsOpt *DBsvcOpts(CsOpt::Instance());
	  if ( DBsvcOpts->getOpt("RUNDB", "user", _svcUser) ) { // first try with included opt file
		  if ( DBsvcOpts->getOpt("RUNDB", "name", _svcServiceName) ) {
			  if ( DBsvcOpts->getOpt("RUNDB", "key", _svcPwd) ) {
				  _svcPwd=decodePassword(_svcPwd);
				  bDBParametersGotFromOptFile = true;
			  }
		  }
	  }
	  if ( ! bDBParametersGotFromOptFile ) { // try default or environment-driven settings
		  cout << "RUN-DB access. No relevant coral-options cards found. Trying to get settings from environment variables or defaults" << endl;
		  _svcUser        = "compass_all";
		  if ( getenv("ORACLE_HOOKUPUSER") ) _svcUser = getenv("ORACLE_HOOKUPUSER"); // VD061218
		  if ( getenv("RUNDB_USER") )  _svcUser = getenv("RUNDB_USER");
		  _svcPwd = decodePassword("9dLfVdJ$+dKeM4k]3!");
		  if ( getenv("ORACLE_HOOKUPPWD") ) _svcPwd = getenv("ORACLE_HOOKUPPWD");    // VD061218
		  
		  if ( getenv("RUNDB_KEY") ) {
			  _svcPwd = getenv("RUNDB_KEY");	
			  _svcPwd=decodePassword(_svcPwd);
		  }
		  _svcServiceName    = "comp_analysis";
		  if ( getenv("ORACLE_HOOKUPDB") ) _svcServiceName = getenv("ORACLE_HOOKUPDB"); // VD03xxxx
		  if ( getenv("RUNDB_NAME") ) _svcServiceName = getenv("RUNDB_NAME");
	  }
	  cout << "RUNDB name: " << _svcServiceName << endl;
	  cout << "RUNDB user: " << _svcUser << endl;
	  cout << "RUNDB pwd: " << /* _svcPwd */ "*******************" << endl;
	  cout << endl;
  }
	  
  return startSession();
}

bool CsOraSession::init(int runNumber, bool isDstProduction)
{
  _initStatus = 0;
  _runStatus  = 0;
  try
    {
      RunInfo myRun(runNumber,false,isDstProduction);
      _svcUser          = myRun.username;
      _svcPwd           = myRun.password;
      _svcServiceName   = myRun.database;
      unsigned year = myRun.year;
      string period = myRun.period;

      cout << "Year:                       " << year     << endl;
      cout << "Period:                     " << period   << endl;
      cout << "DB username:                " << _svcUser << endl;
      cout << "Database:                   " << _svcServiceName << endl;
    }
  catch(...)
    {
      return false;
    }
  return startSession();
}

bool CsOraSession::terminate()
{
  _initStatus = -1;
  try
    {
      //  cout << "terminate statement" << endl;
      if(_conn && _stmt)
	_conn->terminateStatement(_stmt); 
      _stmt = NULL;
      //  cout << "terminate connection" << endl;
      if(_conn && _env)
	_env->terminateConnection(_conn); 
      _conn = NULL;
      //  cout << "terminate environment" << endl;
      if(_env)
	Environment::terminateEnvironment(_env);
      _env = NULL;
    }
  catch(SQLException& e)
    {
      cerr << e.getMessage() << endl;
      return false;
    }
  catch(...)
    {
      cerr << "CsOraSession::terminate: Unknown exception." << endl;
      return false;
    }
  _initStatus = 0;
  return true;
}

bool CsOraSession::commit()
{
  if(!_runStatus) return false;
  try
    {
      _conn->commit();
    }
  catch(SQLException& e)
    {
      cerr << e.getMessage() << endl;
      return false;
    }
  catch(...)
    {
      cerr << "CsOraSession::commit: Unknown exception." << endl;
      return false;
    }
  return true;
}

bool CsOraSession::rollback()
{
  if(!_runStatus) return false;
  try
    {
      _conn->rollback();
      return true;
    }
  catch(SQLException& e)
    {
      cerr << e.getMessage() << endl;
      return false;
    }
  catch(...)
    {
      cerr << "CsOraSession::commit: Unknown exception." << endl;
      return false;
    }
  return false;
}

ResultSet* CsOraSession::executeQuery(const string& sql)
{
  if(!_runStatus) return NULL;
  _sqlStr = sql;
  ResultSet* rs = NULL;
  try
    {
      rs = _stmt->executeQuery(sql);
    }
  catch(SQLException& e)
    {
      cerr << e.getMessage() << endl;
      return NULL;
    }
  catch(...)
    {
      cerr << "CsOraSession::commit: Unknown exception." << endl;
      return NULL;
    }
  return rs;
}

ResultSet* CsOraSession::executeQuery(const string& sql, const int val)
{
  if(!_runStatus) return NULL;
  _sqlStr = sql;
  ResultSet* rs = NULL;
  try
    {
      _stmt->setSQL(sql);
      _stmt->setInt( 1, val );
      rs = _stmt->executeQuery();
    }
  catch(SQLException& e)
    {
      cerr << e.getMessage() << endl;
      return NULL;
    }
  catch(...)
    {
      cerr << "CsOraSession::commit: Unknown exception." << endl;
      return NULL;
    }
  return rs;
}

bool CsOraSession::executeUpdate(const string& sql)
{
  if(!_runStatus) return false;
  _sqlStr = sql;
  try
    {
      _stmt->executeUpdate(sql);
    }
  catch(SQLException& e)
    {
      cerr << e.getMessage() << endl;
      return false;
    }
  catch(...)
    {
      cerr << "CsOraSession::commit: Unknown exception." << endl;
      return false;
    }
  return true;
}

CsOraSession::~CsOraSession()
{
  if(_initStatus == 1 && !uncaught_exception()) 
    {
      if(!stopSession())
	cerr << "Cannot terminate Oracle session." << endl;
    }
  _instance = NULL;
}

bool CsOraSession::stopSession()
{
  if(_runStatus != 1) return true;
  _runStatus = 0;
  return terminate();
}

bool CsOraSession::startSession()
{
  if(_runStatus != 0) return true;
  _initStatus = -1;
  try
    {
      _env = Environment::createEnvironment(Environment::OBJECT);
      if(!_env)
	{
	  cerr << "Cannot create Oracle environment." << endl;
	  return false;
	}
      _conn = _env->createConnection( _svcUser, _svcPwd, _svcServiceName );
      if(!_conn)
	{
	  cerr << "Cannot create Oracle connection." << endl;
	  return false;
	}
      _stmt = _conn->createStatement();
      if(!_stmt)
	{
	  cerr << "Cannot create Oracle statement." << endl;
	  return false;
	}
      _stmt->setAutoCommit(false);
      _initStatus = 1;
      _runStatus = 1;
    }
  catch(SQLException& e)
    {
      cerr << e.getMessage() << endl;
      return false;
    }
  catch(...)
    {
      cerr << "CsOraSession::startSession: Unknown exception." << endl;
      return false;
    }
  return true;
}

string CsOraSession::decodePassword(string encodedPassword) {
	string decodedPwd(encodedPassword);
	for (uint cps(0); cps < encodedPassword.length(); ++cps) {
		decodedPwd[cps] = encodedPassword[cps] ^ ((char)(cps) % 255);
	}
	return decodedPwd;
}
