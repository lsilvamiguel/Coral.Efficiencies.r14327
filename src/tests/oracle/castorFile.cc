#include <string>
#include "RunInfo.h"
#include "CsStoreMisc.h"
#include "CsOraSession.h"
#undef __STRICT_ANSI__
#define __USE_MISC
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <signal.h>
#include <stage_api.h>
#include <rtcp_server.h>
#include <serrno.h>
#define VERBOSE 0

extern "C" void stcplog(int level, char* msg)
{
#if VERBOSE
  fprintf(stdout, "%s", msg);
#endif
}

bool IsOnDisk(const char* fileName)
{
  int nstcp_output;
  struct stgcat_entry *stcp_output = NULL;
  int rc, i, isstaged = 0;
  char* fname = const_cast<char*>(fileName);
  static bool first_entry = true;
  if(first_entry)
    {
      stage_setlog(&stcplog);
      first_entry = false;
    }
  rc = stage_qry_Hsm(STAGE_FILENAME|STAGE_NORETRY,/* Flags */
		     NULL,                        /* Hostname */
		     NULL,                        /* Poolname */
		     fname,                       /* HSM filename */
		     &nstcp_output,               /* Nb stcp output */
		     &stcp_output,                /* Stcp output */
		     NULL,                        /* Nb stpp output */
		     NULL);                       /* Stpp output */
  if (rc != 0) 
    {
      //printf("### stage_qry_Hsm error (%s)", sstrerror(serrno));
      return false;
    } 
  else 
    {
      for (i = 0; i < nstcp_output; i++) 
	{
	  if (ISSTAGED((&stcp_output[i]))) isstaged++;
	}
    }
  if (stcp_output != NULL) free (stcp_output);   /* User responsability ! */

  return (isstaged!=0);
}

using namespace std;

int main(int argc, char* argv[])
{
  if(argc != 2)
    {
      printf("USAGE: castorFile run_number\n");
      return 1;
    }
  int runNumber = 0;
  int retVal = sscanf(argv[1],"%d",&runNumber);
  if(retVal != 1 || runNumber < 11961 || runNumber > 99999)
    {
      printf("Incorrect run number: %s\n",argv[1]);
      return 2;
    }
  try 
    {
      RunInfo run(runNumber);
      printf("Looking for CASTOR location of:\n"
	     "Run number:                 %d\n",
	     run.runNumber);
      CsOraSession* session = CsOraSession::instance();
      if(!session->init(run.runNumber,false))
	{
	  return 4;
	}
      string query = string("select t_min,t_max from runs where run_number=") + itostr(run.runNumber);
      ResultSet* result = session->executeQuery(query);
      if(!result || !result->next())
	{
	  printf("Run number %d has not been found in Oracle DB\n",run.runNumber);
	  return 5;
	}
      else
	{
	  time_t t_min = result->getUInt(1);
	  time_t t_max = result->getUInt(2);
	  printf("Start time:                 %s",ctime(&t_min));
	  printf("Stop  time:                 %s",ctime(&t_max));
	}
      query = string("select file_dir,file_name,file_size from file_maps where run_number=") + itostr(run.runNumber)
	+ " order by file_name";
      result = session->executeQuery(query);
      if(!result || !result->next())
	{
	  printf("No any chunks found for the run number %d in DB\n", run.runNumber);
	  return 6;
	}
      printf("----------------------------------------------------------------------------------------------------\n");
      printf("   CASTOR file name                                                                 Size  Location  \n");
      printf("----------------------------------------------------------------------------------------------------\n");
      int nFiles = 0;
      unsigned long long total_size = 0;
      do
	{
	  string file_dir    = result->getString(1);
	  string file_name   = result->getString(2);
	  unsigned file_size = result->getUInt(3);
	  nFiles++;
	  total_size += file_size;
	  string onDisk = "UNKNOWN";
	  if(IsOnDisk((file_dir+file_name).c_str()))
	    onDisk = "onDISK";
	  else
	    onDisk = "onTAPE";
	  printf("%-74s   %11u  %s\n",(file_dir+file_name).c_str(),file_size,onDisk.c_str());
	} while(result->next());
      printf("----------------------------------------------------------------------------------------------------\n");
      printf("Total number of files: %d\n",nFiles);
      printf("Total size:            %.3f GB (1GB = %u bytes)\n",float(total_size/1024/1024)/1024.,1024*1024*1024U); 
      printf("----------------------------------------------------------------------------------------------------\n");
    }
  catch( int err)
    {
      return 3;
    }
  return 0;
}
