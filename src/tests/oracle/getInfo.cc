#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "CsOraSession.h"
#include "CsTime.h"
#include "RunInfo.h"

using namespace std;
using namespace oracle::occi;

static CsOraSession* session = NULL;

string itostr(int i) {
   char	buff[128];
   sprintf(buff, "%d", i);
   return buff;
}


struct EventInfo
{
  unsigned event_number;
  unsigned event_in_burst;
  unsigned trigger_mask;
  unsigned time_sec;
  unsigned time_usec;
  unsigned error_code;
  unsigned event_size;
  unsigned rawev_file_id;
  unsigned rawev_file_offset;
};

int main(int argc, char* argv[])
{
  if(argc == 1)
    {
      cout << "Usage getInfo run_num burst_num ev_in_burst" << endl;
      return 0;
    }
  string runN = argv[1];
  string burstN = argv[2];
  string evInBurst = argv[3];

  int runNumber = 0;
  int burstNumber = 0;
  int eventInBurst = 0;

  sscanf(runN.c_str(),"%d",&runNumber);
  sscanf(burstN.c_str(),"%d",&burstNumber);
  sscanf(evInBurst.c_str(),"%d",&eventInBurst);

  cout << "Wait, please" << endl;
  try
    {
      RunInfo run(runNumber,false);
      string userName = run.username;
      string connectString = run.database;
      session = CsOraSession::instance();
      string password = run.password;
      session->init(runNumber,false);

      string query = string("select file_name,file_dir from file_maps ")
	+ "where file_id in (select rawev_file_id from event_headers where"
	+ " rawev_file_id in (select file_id from file_maps where run_number="
	+ runN + ") and burst=" + burstN + " and event_in_burst=" + evInBurst + ")";
//      cout << query << endl;
      ResultSet* rs = session->executeQuery(query);

      if(!rs || !rs->next())
	{
	  cerr << "No any information about event " << evInBurst << " in burst "
	       << burstN << " of run number " << runN << endl;
	  return -1;
	}
      string fileName = rs->getString(1);
      string fileDir = rs->getString(2);
      cout << endl;
      cout << "Run number     = " << runN << endl
	   << "database       = " << run.database << endl
	   << "account        = " << run.username << endl
 	   << "Burst number   = " << burstN << endl
	   << "Event in burst = " << evInBurst << endl
	   << "RAW file name  = " << fileName << endl
	   << "Directory name = " << fileDir << endl;
    }
  catch (SQLException& e)
    {
      cerr << "ORACLE session error: " << e.getMessage() << endl;
      return -2;
    }
  catch(...)
    {
      cerr << "Unexpected error." << endl;
    }
  return 0;
}
