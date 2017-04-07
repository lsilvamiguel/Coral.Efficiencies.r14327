#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "CsOraSession.h"
#include "RunInfo.h"

using namespace std;
using namespace oracle::occi;

static CsOraSession* session = NULL;

string itostr(int i) {
   char	buff[128];
   sprintf(buff, "%d", i);
   return buff;
}

string ftostr(float i) {
   char	buff[128];
   if((i-int(i)) != 0.)
     sprintf(buff, "%.2f", i);
   else
     sprintf(buff, "%d", int(i));     
   return buff;
}

void Usage()
{
  cout << "Usage: mDstInfo dst_info_file_name" << endl;
}

const char* scratch_dir = "/afs/cern.ch/compass/scratch/mdst";

int main(int argc, char* argv[])
{
  if(argc != 2)
    {
      Usage();
      cerr << "mDstInfo parameters are wrong." << endl;
      return -1;
    }
  int runNumber = 0;
  string rn = getenv("DST_RUN_NUMBER");
  sscanf(rn.c_str(),"%d",&runNumber);
  RunInfo run(runNumber,false);
  map<int,string> dst_chunks;
  map<int,int> mdst_chunks;
  try
    {
      session = CsOraSession::instance();
      session->init(runNumber,false);
    }
  catch (SQLException& e)
    {
      cerr << "ORACLE session error: " << e.getMessage() << endl;
      return -2;
    }  
  string query = string("select count(*) from FILE_MAPS where "
			"FILE_TYPE='RAW' and RUN_NUMBER=") + rn;
  ResultSet* rs = session->executeQuery(query);
  if(!rs || !rs->next())
    {
      cerr << "No any chunks for the run:" << runNumber << " in DB." << endl;
      return -3;
    }
  int nChunks = rs->getInt(1);
  string slotNumber = getenv("DST_SLOT_NUMBER");
  query = string("select file_name from dst_files where run_number=") + rn
    + " and DST_VERSION=" + slotNumber;
  rs = session->executeQuery(query);
  while(rs && rs->next())
    {
      string chunk_name = rs->getString(1);
      int run_num = 0;
      int chunk_num = 0;
      int slot_num = 0;
      if(sscanf(chunk_name.c_str(),"cdr%05d-%05d.dst%d", &chunk_num, &run_num,
		&slot_num) == 3 && run_num == runNumber 
	 && itostr(slot_num) == slotNumber)
	{
	  dst_chunks[chunk_num] = chunk_name;
	}
    }
  string output_dir = string(scratch_dir) + "1/";
  FILE* dstInfo = fopen(argv[1], "r");
  if(!dstInfo)
    {
      cerr << "Cannot find dst info file: " << argv[1] << endl;
      return -4;
    }
  char str[1024];
  int nline = 0;
  vector<string> mDstVec;
  while(fgets(str,1023,dstInfo))
    {
      nline++;
      if(nline == 1 && strstr(str,"Total number of chunks:")) // for compatibility with previous version
	{
	  continue;
	}
      if(!strstr(str,"mDST"))
	continue;
      mDstVec.push_back(string(str));
    }
  fclose(dstInfo);
  vector<string> goodChunks;
  int nGoodChunks = 0;
  for(unsigned i = 0; i < mDstVec.size(); i++)
    {
      string delimiter = " ";
      std::string::size_type idx = mDstVec[i].find(delimiter);
      if(idx != std::string::npos)
	{
	  int value = -1;
	  int retVal = sscanf(&(mDstVec[i].c_str()[idx]),"%d",&value);
	  if(retVal == 1 && (value == 0 || value == 1))
	    {
	      if(value)
		{
		  string name = mDstVec[i].substr(0,idx);
		  goodChunks.push_back(name);
		  string d = string("/mDST-") + rn;
		  std::string::size_type _idx = name.find(d);
		  if(_idx != std::string::npos)
		    {
		      int slot = 0;
		      int r = 0;
		      int ch = 0;
		      float version = 0.;
		      if(sscanf(name.substr(_idx).c_str(),"/mDST-%5d-%5d-%d-%f", 
				&r, &ch, &slot, &version) != 4)
			{
			  cerr << "Invalid mDST file name: " << name << endl;
			}
		      else
			{
			  mdst_chunks[ch]++;
			  if(mdst_chunks[ch] == 1)
			    nGoodChunks++;
			}
		    }
		}
	    }
	  else
	    {
	      cerr << "mDstInfo: Bad line: " << mDstVec[i] << endl;
	    }
	}
      else
	{
	  cerr << "mDstInfo: Bad line: " << mDstVec[i] << endl;
	}
    }
  cout << "mDstInfo: Total number of chunks what have been finished: " << goodChunks.size() << endl;
  if(dst_chunks.size() == (unsigned)nChunks)
    {
      if(nGoodChunks == nChunks && goodChunks.size() == (unsigned)nChunks)
	{
	  string arg = string(argv[1]);
	  string run = arg.substr(arg.find("Run_"),9);
	  string out_fn = output_dir + "__" + run + ".info";
	  FILE* output = fopen(out_fn.c_str(),"a");
	  for(unsigned i = 0; i < goodChunks.size(); i++)
	    fprintf(output,"%s\n",goodChunks[i].c_str());
	  fclose(output);
	  return 2;
	}
      else
	{
	  string arg = string(argv[1]);
	  string run = arg.substr(arg.find("Run_"),9);
	  string err_fn = output_dir+"__"+run+".error";
	  FILE* error_file = fopen(err_fn.c_str(),"w");
	  map<int,string>::iterator it;
	  for(it = dst_chunks.begin(); it != dst_chunks.end(); ++it)
	    {
	      if(mdst_chunks[it->first] == 0)
		{
		  fprintf(error_file,"mDST for the chunk %s could not be found\n", it->second.c_str());
		}
	      else if(mdst_chunks[it->first] > 1)
		{
		  fprintf(error_file,"%d mDST chunks for chunk %s in list.\n", mdst_chunks[it->first],
			  it->second.c_str());
		}
	    }
	  fclose(error_file);
	}
    }
  return 0;
}
