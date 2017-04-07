#include <iostream>
#include <cstdio>
#include "RunInfo.h"
#include "CsOraSession.h"
#include "OracleErrorCode.h"
#include "CsStoreMisc.h"

using namespace std;
using namespace oracle::occi;

int RunInfo::GetFirstRunOfPeriod(const string& period)
{
  CsOraSession* session = CsOraSession::instance();
  bool retVal = session->initInfoSession();
  if(!retVal)
    {
      cerr << "RunInfo::GetFirstRunOfPeriod: Oracle DB session error." << endl;
      return 0;
    }
  string query = string("select min(run_number) from RUNS_ALL_MV where PERIOD like '%") + period
    + "%'";
  ResultSet* result = session->executeQuery(query);

  int run_number = 0;

  if(result && result->next())
    {
      run_number = result->getInt(1);
    }
  if(!session->stopSession())
    {
      cerr << "RunInfo::GetFirstRunOfPeriod: Cannot close Oracle DB." << endl;
      return 0;
    }
  return run_number;
}

RunInfo::RunInfo(int _runNumber, bool isTestDB, bool isDstProduction)
{
  runNumber = _runNumber;
  bool isFound = false;
  CsOraSession* session = CsOraSession::instance();
  bool retVal = session->initInfoSession();
  if(!retVal)
    {
      cerr << "RunInfo::RunInfo: Oracle DB session error." << endl;
      throw ORACLE_ERROR;
    }

  string query = string("select DB_NAME, PERIOD from RUNS_ALL_MV where RUN_NUMBER=:1");
  ResultSet* result = session->executeQuery(query, runNumber);
  if(result && result->next())
    {
      database = result->getString(1);
      period   = result->getString(2);
      if(isDstProduction)
	{
	  username = "COMPDST_";
	  password = "dstwriter";
	}
      else
	{
	  username = "COMPANL_";
	  password = "compass";
	}
		// avoid hard-codings: read "decoded" pwds from special table RUNDB_SCHEMAS
	   query = string("select SCHEMA_KEY from RUNDB_SCHEMAS where SCHEMA_NAME='"+ username +"'");
		result = session->executeQuery(query);
		if(result && result->next())
			{
				isFound = true;
				// cout << "===================== RUN-DB. Pwd read from DB (pattern name)" << endl;
				password = session->decodePassword(result->getString(1));
			} // silently skip this error in case no prefix-only username cannot be found

		username += period; 

		if ( ! isFound ) // if prefix-only username's pwd was not found try with full username
			{
				query = string("select SCHEMA_KEY from RUNDB_SCHEMAS where SCHEMA_NAME='"+ username +"'");
				result = session->executeQuery(query);
				if(result && result->next())
					{
						isFound = true;
						//cout << "===================== RUN-DB. Pwd read from DB (exact name)" << endl;
						password = session->decodePassword(result->getString(1));
					} 
				else
					{
						cerr << "RunInfo::RunInfo: Oracle error. Cannot connect as " << username << endl;
						throw ORACLE_ERROR;
					}	
			}
      char p[16];
      if(sscanf(period.c_str(),"%02d%s",&year,p) != 2)
	{
	  cerr << "RunInfo::RunInfo: Incorrect period name: " << period << endl;
	  throw ORACLE_ERROR;
	}
      year += 2000;
      period = string(p);
      if(!session->stopSession())
	{
	  cerr << "RunInfo::RunInfo: Cannot close Oracle DB." << endl;
	  throw ORACLE_ERROR;
	}
    }
  else
    {
      cerr << "RunInfo::RunInfo: Oracle DB error." << endl;
      throw ORACLE_ERROR;
    }
  if(!isFound)
    {
      cerr << "RunInfo::RunInfo: Run number " << runNumber << " not found in database." << endl;
      throw RUN_NOT_FOUND_IN_DB;
    }
}
