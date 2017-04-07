// Basic example of data storage.

#ifndef NO_OBJY  
#include "oo.h"
#include "ConditionsDB/ICondDBMgr.h"
#include "ConditionsDB/CondDBObjyDBMgrFactory.h"
#include "ConditionsDB/CondDBObjFactory.h"
#endif

#include <string>
#include <iostream>
#include <strstream>
#include <fstream>
#include <unistd.h>
#include <sys/param.h>


int main ( int argc, char* argv[] )
{

#ifndef NO_OBJY  
  try {

    string user = getenv("USER");
    
#ifdef STANDALONE
    ooNoLock();
#endif

    ICondDBMgr* condDBmgr = CondDBObjyDBMgrFactory::createCondDBMgr();
    condDBmgr->init(argv[1]);
    ICondDBDataAccess* condDataAccess = condDBmgr->getCondDBDataAccess();
    ICondDBFolderMgr*  condFolderMgr  = condDBmgr->getCondDBFolderMgr(); 

    condDBmgr->startUpdate();    
    condDBmgr->createCondDB();
 
    condDBmgr->openDatabase();

    ifstream ifs;
    ifs.open(argv[2]);
    string folder = "";
    string atribute;
    string dbFile = "";
    string description = "";
    ifs>>folder;
    if(ifs.is_open()){
      do{
	dbFile+="."+user;
	//	cout << folder << " " << dbFile << " " << description << endl;
      	condFolderMgr->createCondDBFolder(folder,atribute,description);
	ifs>>folder;
      }while(!ifs.eof());
      ifs.close();
    }

    condDBmgr->commit();

  return 0;
  }
  catch (CondDBException &e)
    {
      cerr << "*** ConditionsDB exception caught: " << e.getMessage() << "\n"
           << "***   error code: " << e.getErrorCode() << endl;
      
      return 1;  // return failure
    }

  catch(char *pErr){
    cout << "Error:" << pErr << endl; 
  }
  catch(string str){
    cout << "Error:" << str << endl; 
  }
#else
    cout << " You need to compile with Objectivity/DB!" << endl;
#endif
}



