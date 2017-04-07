#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <string>
#include <vector>
#include <iostream>
#define linux 1
#define _LARGEFILE64_SOURCE 1
#define __USE_LARGEFILE64 1
#include <sys/types.h>
#include <rfio.h>
#include <serrno.h>
#include <dirent.h>

#include "CsOraSession.h"
#include "CsStoreMisc.h"
#include "RunInfo.h"

using namespace std;
using namespace oracle::occi;

struct DstInfo
{
  int      runNumber;
  string   file_dir;
  string   file_name;
  unsigned file_size;
};

static CsOraSession* session = NULL;

int main(int argc, char* argv[])
{
  if(argc == 1)
    {
      cout << "Usage: oradbTest period [-p|-k]" << endl;
      cout << "\t-p  does not ask user before deleting the chunk from DB" << endl;
      cout << "\t-k  keep all chunks in DB" << endl;
      return 1;
    }
  bool prompt = true;
  bool keep = false;
  bool DB_changed = false;
  if(argc == 3 && strcmp(argv[2],"-p") == 0)
    prompt = false;
  if(argc == 3 && strcmp(argv[2],"-k") == 0)
    keep = true;
  string period = argv[1];
  int first = RunInfo::GetFirstRunOfPeriod(period);
  if(!first)
    {
      cerr << "No runs are found in DB." << endl;
      return 2;
    }
  RunInfo run(first,false,true); // get info for openning session in  READ-WRITE mode
  try
    {
      session = CsOraSession::instance();
      session->init(first,true);
    }
  catch (SQLException& e)
    {
      cerr << "ORACLE session error: " << e.getMessage() << endl;
      return 4;
    }  
 
  string query = string("select run_number from runs order by run_number");
  ResultSet* rs = session->executeQuery(query);
  if(!rs || !rs->next())
    {
      cerr << "No any runs in period" << argv[1] << " in DB." << endl;
      return 5;
    }
  vector<int> runs;
  do
    {
      runs.push_back(rs->getInt(1));
    } while(rs->next());
  vector<DstInfo> problems;
  for(unsigned i = 0; i < runs.size(); i++)
    {
      cout << "Run: " << runs[i] << endl;
      query = string("select file_dir, file_name, file_size from dst_files where run_number=")
	+ itostr(runs[i])+" and dst_version=1"; // only slot=1 for the moment
      rs = session->executeQuery(query);
      if(!rs || !rs->next())
	{
	  //	  cout << "No DST for the run number " << runs[i] << endl;
	  continue;
	}
      vector<DstInfo> dstInfo;
      do
	{
	  DstInfo info;
	  info.runNumber = runs[i];
	  info.file_dir = rs->getString(1);
	  info.file_name = rs->getString(2);
	  info.file_size = rs->getUInt(3);
	  dstInfo.push_back(info);
	} while(rs->next());
      if(dstInfo.size() == 0)
	{
	  //	  cout << "No DST for the run number " << runs[i] << endl;
	  continue;
	}
      for(unsigned j = 0; j < dstInfo.size(); j++)
	{
	  struct stat st;
	  string fileName = dstInfo[j].file_dir + dstInfo[j].file_name;
	  int retVal = rfio_stat(const_cast<char*>(fileName.c_str()),&st);
	  if(retVal != 0)
	    {
	      cerr << "File: " << fileName << " does not exist. serrno=" << serrno << endl;
	      problems.push_back(dstInfo[j]);
	      continue;
	    }
	  if(st.st_size != (int)dstInfo[j].file_size)
	    {
	      cerr << "Size of file: " << fileName << " on castor is differ with DB size: " <<
		st.st_size << " : " << dstInfo[j].file_size << endl;
	      if(!keep)
		{
		  bool remove = true;
		  bool quit = false;
		  if(prompt)
		    {
		      cout << "Do you want remove this chunk from database? (Yes|No|Quit):  " << endl;
		      int c = getchar();
		      if(c == 'y' || c == 'Y')
			remove = true;
		      else
			remove = false;
		      if(c == 'q' || c == 'Q')
			quit = true;
		    }
		  if(quit) exit(0);
		  if(remove)
		    {
		      query = string("delete from dst_files where file_name='") + dstInfo[j].file_name
			+ "' and RUN_NUMBER=" + itostr(dstInfo[j].runNumber)
			+ " and FILE_DIR='" + dstInfo[j].file_dir + "'";
		      session->executeUpdate(query);
		      DB_changed = true;
		    }
		}
	      problems.push_back(dstInfo[j]);
	    }
	}
    }
  if(problems.size())
    {
      FILE* output = fopen("dst_error.dat","w");
      for(unsigned i = 0; i < problems.size(); i++)
	{
	  fprintf(output,"%s%s\n",problems[i].file_dir.c_str(),problems[i].file_name.c_str());
	}
      fclose(output);
    }
  if(DB_changed)
    session->commit();
  return 0;
}
