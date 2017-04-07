#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <ctime>

#include "CsOraSession.h"
#include "CsStoreMisc.h"
#include "RunInfo.h"

using namespace std;
using namespace oracle::occi;

static CsOraSession* session = NULL;

struct DstInfo
{
  double variable[3];
};

int main(int argc, char* argv[])
{
  if(argc == 1)
    {
      cout << "Usage: runInfo run_number" << endl;
      return 1;
    }
  int runNumber = atoi(argv[1]);
  try
    {
      session = CsOraSession::instance();
      session->init(runNumber,false);      
      cout << "Get information about run " << runNumber << endl;
      cout << "."; cout.flush();
      string query = string("select t_min, t_max from runs where run_number=")
	+ itostr(runNumber);
      ResultSet* rs = session->executeQuery(query);
      time_t start = 0;
      time_t stop = 0;
      if(rs && rs->next())
	{
	  start = rs->getUInt(1);
	  stop = rs->getUInt(2);
	}
      
      query = string("select count(*) from event_headers where rawev_file_id in "
		     "(select file_id from file_maps where run_number=")+itostr(runNumber)
	+ " and file_type='RAW')";
      cout << "."; cout.flush();
      rs = session->executeQuery(query);
      int numOfRawEvents = 0;
      if(rs && rs->next())
	{
	  numOfRawEvents = rs->getInt(1);
	}
      query = string("select count(*) from FILE_MAPS where "
		     "FILE_TYPE='RAW' and RUN_NUMBER=") + itostr(runNumber);
      cout << "."; cout.flush();
      rs = session->executeQuery(query);
      int numOfRawChunks = 0;
      if(rs && rs->next())
	{
	  numOfRawChunks = rs->getInt(1);
	}
      query = string("select count(*) from DST_FILES where ")
	+ "RUN_NUMBER=" + itostr(runNumber);
      cout << "."; cout.flush();
      // cout << "Query: " << query << endl;
      rs = session->executeQuery(query);
      int numOfDstChunks = 0;
      int rawFileId = 0;
      if(rs && rs->next())
	{
	  numOfDstChunks = rs->getInt(1);
	}
      query = string("select dst_version, dst_type_id, value1_is, value2_is, value3_is")
	+ " from dst_files where run_number=" + itostr(runNumber);
      int slotNumber = 0;
      int dstVersion = 0;
      string varName1;
      string varName2;
      string varName3;
      cout << "."; cout.flush();
      // cout << "Query: " << query << endl;
      rs = session->executeQuery(query);
      if(rs && rs->next())
	{
	  slotNumber = rs->getInt(1);
	  dstVersion = rs->getInt(2);
	  varName1 = rs->getString(3);
	  varName2 = rs->getString(4);
	  varName3 = rs->getString(5);
	}
      vector<DstInfo> dstInfo;
      query = string("select value1, value2, value3 from dst where file_id in ")
	+ "(select file_id from dst_files where run_number=" + itostr(runNumber) + ")";
      Statement* stmt = session->conn()->createStatement();
      stmt->setSQL(query);
      stmt->setPrefetchRowCount(1000);
      rs = stmt->executeQuery(query);
      if(rs)
	{
	  while(rs->next())
	    {
	      DstInfo di;
	      di.variable[0] = rs->getDouble(1);
	      di.variable[1] = rs->getDouble(2);
	      di.variable[2] = rs->getDouble(3);
	      dstInfo.push_back(di);
	    }
	}
      stmt->closeResultSet(rs);
      session->conn()->terminateStatement(stmt);
      cout << "."; cout.flush();
      int nEventsWith1 = 0;
      int nEventsWith2 = 0;
      int nEventsWith3 = 0;
      int nDstEvents = (int)dstInfo.size();
      for(int i = 0; i < nDstEvents; i++)
	{
	  if(dstInfo[i].variable[0] > 0.)
	    nEventsWith1++;
	  if(dstInfo[i].variable[1] > 0.)
	    nEventsWith2++;
	  if(dstInfo[i].variable[2] > 0.)
	    nEventsWith3++;
	}
      cout << endl;
      cout << "Run number:                         " << runNumber << endl;
      cout << "Started at:                         " << ctime(&start);
      cout << "Finished at:                        " << ctime(&stop);
      cout << "Number of chunks:                   " << numOfRawChunks << endl;
      cout << "Number of raw events:               " << numOfRawEvents << endl;
      cout << "Number of dst chunks:               " << numOfDstChunks << endl;
      cout << "Number of dst events:               " << nDstEvents << endl;
      cout << "DST version:                        " << dstVersion << endl;
      cout << "DST slot number:                    " << slotNumber << endl;
      cout << "External variable #1 is:            " << varName1 << endl;
      cout << "External variable #2 is:            " << varName2 << endl;
      cout << "External variable #3 is:            " << varName3 << endl;
      cout << "Number of dst events with var1>0:   " << nEventsWith1 << endl;
      cout << "Number of dst events with var2>0:   " << nEventsWith2 << endl;
      cout << "Number of dst events with var3>0:   " << nEventsWith3 << endl;
      cout << "For additional information use tora." << endl;
    }
  catch (SQLException& e)
    {
      cerr << "ORACLE session error: " << e.getMessage() << endl;
      return 2;
    }
  return 0;
}
