#undef __STRICT_ANSI__
#include "CsOraSession.h"
#include "RunInfo.h"
#include <string>
#include <cstdio>
#include <iostream>

using namespace std;
using namespace oracle::occi;

string itostr(int i) {
   const int buflen = 100;
   char	buff[buflen];
   snprintf(buff, buflen, "%d", i);
   return buff;
}

int main(int argc, char* argv[])
{
  CsOraSession* session = CsOraSession::instance();
  string runNumber;
  string chunkNumber;
  string slotNumber;
  string password;
  if(argc == 5)
    {
      runNumber = argv[1];
      chunkNumber = argv[2];
      slotNumber = argv[3];
      password = argv[4];
    }
  else if(argc == 4)
    {
      runNumber = argv[1];
      chunkNumber = argv[2];
      slotNumber = argv[3];
      cout << "Enter ORACLE password, please: ";
      cin >> password;
   }
  else if(argc == 3)
    {
      runNumber = argv[1];
      chunkNumber = argv[2];
      cout << "Enter slot number: ";
      cin >> slotNumber;
      cout << "Enter ORACLE password, please: ";
      cin >> password;
    }
  else if(argc == 2)
    {
      runNumber = argv[1];
      cout << "Enter chunk number: ";
      cin >> chunkNumber;
      cout << "Enter slot number: ";
      cin >> slotNumber;
      cout << "Enter ORACLE password, please: ";
      cin >> password;
    }
  else
    {
      cout << "Enter run number: ";
      cin >> runNumber;
      cout << "Enter chunk number: ";
      cin >> chunkNumber;
      cout << "Enter slot number: ";
      cin >> slotNumber;
      cout << "Enter ORACLE password, please: ";
      cin >> password;
    }
  int rn = 0;
  sscanf(runNumber.c_str(), "%5d", &rn);
  RunInfo run(rn,false,true); // get info for openning session in READ-WRITE mode
  cout << "Remove " << chunkNumber << " chunk of run number " << runNumber <<
    " slot number " << slotNumber << ".\nConnection: " << run.username << "@" << run.database <<
    " ORACLE password: " << password << endl;
  if(password != run.password)
    {
      cerr << "Incorrect password: " << password << endl;
      return -1;
    }
  try
    {
      string userName = run.username;
      string connectString = run.database;
      
      cout << "Try to create connection to the ORACLE database "
	   << connectString << " as user: " << userName << "...." << endl;

      session->init(rn,true);

      string fileName;

      if(chunkNumber == "all")
	{
	  chunkNumber = "ALL chunks";
	  fileName = string("cdr")+"*"+"-"+runNumber+".dst"+slotNumber;
	}
      else
	{
	  fileName = string("cdr")+chunkNumber+"-"+runNumber+".dst"+slotNumber;
	  chunkNumber = string("chunk ") + chunkNumber;
	}
      cout << "WARNING: All DST data for " << chunkNumber << ", slot " << slotNumber
	   << " of the run number " << runNumber
	   << " will be destoyed!!!" << endl;
      cout << "Are you sure that you would like to remove DST file: " << fileName << " (yes/no) ";
      string answer;
      cin >> answer;
      if(answer == "yes")
	{
	  cout << "As you wish." << endl;
	}
      else
	{
	  cout << "Good bye!" << endl;
	  return 0;
	}

      string query;

      if(chunkNumber == "ALL chunks")
	{
	  query = string("select file_id, file_name from dst_files where RUN_NUMBER=")
	    + runNumber + " and DST_VERSION=" + slotNumber;
	  ResultSet* rs = session->executeQuery(query);
	  vector<int> fileid;
	  vector<string> filenames;
	  if(!rs || !rs->next())
	    {
	      cerr << "No any chunks for the run " << runNumber << " in DB." << endl;
	    }
	  else
	    {
	      do
		{
		  fileid.push_back(rs->getInt(1));
		  filenames.push_back(rs->getString(2));
		} while(rs->next());
	    }
	  for(unsigned i = 0; i < fileid.size(); i++)
	    {
	      query = string("delete from DST_FILES where FILE_ID=")
		+ itostr(fileid[i]);
	      cout << "Execute query:\n" << query << endl;
	      session->executeUpdate(query);
	      session->commit();
	      cout << "File " << filenames[i] << " has been removed from DB." << endl;
	    }
	}
      else
	{
	  query = string("delete from DST_FILES where FILE_NAME='")
	    + fileName + "' and RUN_NUMBER=" + runNumber +
	    " and DST_VERSION=" + slotNumber;
	  cout << "Execute query:\n" << query << endl;
	  session->executeUpdate(query);
	  cout << "done!" << endl;
	}
      session->commit();
    }
  catch (SQLException& e)
    {
      cerr << "ORACLE session: " << e.getMessage() << endl;
      return -1;
    }
  return 0;
}
