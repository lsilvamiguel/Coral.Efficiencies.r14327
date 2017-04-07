#include <string>
#include <stdio.h>
#include "RunInfo.h"
#include "CsStoreMisc.h"
#include "CsOraSession.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc != 2)
    {
      printf("USAGE: removeTBnames <run number>\n");
      return 1;
    }
  int runNumber = 0;
  int retVal = sscanf(argv[1],"%d",&runNumber);
  if(retVal != 1 || runNumber < 11961 || runNumber > 99999)
    {
      printf("Incorrect run number: %s\n",argv[1]);
      return 2;
    }
  printf("Program will erase TBnames string for run number %d from Oracle DB.\n Are you sure? (y/n)\n", runNumber);
  int answer = getchar();
  if(answer == 'y' || answer == 'Y');
  else
    {
      printf("Erasing aborted.\n");
      return 3;
    }
  try 
    {
      RunInfo run(runNumber);
      CsOraSession* session = CsOraSession::instance();
      if(!session->init(run.runNumber,true))
	{
	  return 4;
	}
      string query = string("update runs set log_info=NULL where run_number=") + itostr(run.runNumber);
      if(!session->executeUpdate(query))
	{
	  return 5;
	}
      session->commit();
    }
  catch( int err)
    {
      return 3;
    }
  printf("Done.\n");
  return 0;
}

