#if USE_FileStore
#include "CsFileStore.h"
#include "CsFileRun.h"
#include "CsOpt.h"
#include "CsRegistry.h"
#include "CsGeom.h"

#include <sys/types.h>
#include <dirent.h>
#include "rfio_api.h"

using namespace std;

static bool TestTBNameString(const string& tbNames, const string& curTBNames, int dstSlot)
{
  const string startTBNStr = ":Start of TBNames string:";
  const string endTBNStr   = ":End of TBNames string:";
  if(tbNames == "")
    {
      cerr << "CsOraStore::scan: No TBNames string into RUNS.log_info "
	   << "for the slot " << dstSlot << "\n"
	   << "You have to write it before DST production. " << endl;
      return false;
    }
  string tmp = tbNames;
  std::string::size_type idx = 0;
  for(int i = 0; i < 3; i++)
    {
      std::string::size_type idx_start = tmp.find( startTBNStr , idx );
      if(idx_start == std::string::npos) break;
      std::string::size_type idx_end = tmp.find(endTBNStr, idx_start);
      if(idx_end == std::string::npos)
	{
	  cerr << "TestTBNameString(): error during test of TBNames string: end-tag not found." << endl;
	  return false;
	}
      idx_end += endTBNStr.length();
      string TBNames = tmp.substr(idx_start, idx_end-idx_start);
      if(i == (dstSlot-1))
	{
	  if(TBNames == curTBNames) 
	    return true;
	  cerr << "The TBNames string into RUNS.log_info "
	       << "is not equil TBNames for the current slot number.\n"
	       << "You have to change detectors.dat file for the current "
	       << "DST production or update RUNS.log_info." << endl;
	  cerr << "TBNames into RUNS.log_info:\n"
	       << tbNames << "\n\n";
	  cerr << "TBNames for current detectors.dat: \n"
	       << curTBNames << endl;
	  return false;
	}
      idx = idx_end;
    }
  cerr << "TestTBNameString(): TBNames string for slot " << dstSlot << " not found." << endl;
  return false;
}

CsFileStore* CsFileStore::Instance()
{
  static  CsFileStore* _instance = new CsFileStore;
  return _instance;
}

CsFileStore::CsFileStore()
{
  isInit = isScanned = false;
  runNumber  = 0;
  dstSlot    = 0;
  dstVersion = 0;
  fileRun    = NULL;
}

CsFileStore::~CsFileStore()
{
  delete fileRun;
  isInit = isScanned = false;
  runNumber  = 0;
  dstSlot    = 0;
  dstVersion = 0;
  fileRun    = NULL;
}

bool CsFileStore::init()
{
  if(isInit) 
    {
      cerr << "CsFileStore::init init already performed " << endl;
      return(true);
    }

  CsOpt*  myOptions = CsOpt::Instance();
  
  string dType = "";

  if(!myOptions->getOpt("Data", "type", dType)) 
    {
      cerr << "CsFileStore::init: Data type not selected..." << endl;
      return(false);
    }
  else
    {
      for(unsigned p = 0; p < dType.size(); p++) dType[p] = tolower(dType[p]);
      if(dType != "dst")
	{
	  cerr << "CsFileStore::init: incorrect data type selected: "
	       << dType << endl;
	  return false;
	}
    }

  string directory;
  
  if(!myOptions->getOpt("Data", "directory", directory ))
    {
      cerr << "CsFileStore::init: Data directory parameter is not defined." << endl;
      return false;
    }

  string container;
  if(myOptions->getOpt("Data", "container", container )) 
    {
      chunkName = container.substr(3);
      cout << "Data chunk " << chunkName << " has been selected.\n" 
	   << "Selected run number will be ignored." << endl;
      cout << "Container name:             " << container << ".raw" << endl;
      cout << "Chunk number:               " << chunkName << endl;
    }

  if(chunkName == "" && myOptions->getOpt("Data run", "select", runNumber))
    {
      cout << "All chunks for run " << runNumber << " have been selected." << endl;
    }
  else if(chunkName == "")
    {
      cerr << "CsFileStore::init: No any chunk selected." << endl;
      return false;
    }
  if(chunkName != "") // chunk has been selected
    {
      int chunkNum = 0;
      if(sscanf(chunkName.c_str(),"%05d-%d",&chunkNum, &runNumber) != 2)
	{
	  cerr << "CsFileStore::init: Incorrect chunk name: " <<
	    chunkName << endl;
	  return false;
	}
    }
  if(CsOpt::Instance()->getOpt( "", "store TBNames on DST only" ) )
    {
      cerr << "CsFileStore::init(): Incorrect option: \"store TBNames on DST only\"" << endl;
      return false;
    }

  if(myOptions->getOpt("make", "DST", dstVersion))
    {
      cerr << "CsFileStore::init(): incorrect option \"make DST " << dstVersion << "\"" << endl;
      return false;
    }

  if(myOptions->getOpt("Dst select", "slot", dstSlot))
    {
      if(dstSlot < 1 || dstSlot > 3)
	{
	  cerr << "CsFileStore::init: Incorrect DST slot " << dstSlot
	       << " has been selected." << endl;
	  return false;
	}
      cout << "DST slot:                   " << dstSlot << endl;
    }
  else
    {
      cerr << "CsFileStore::init: You have to select DST slot number." << endl;
      return false;
    }

  if(!ScanDirectory(directory,runNumber,chunkName))
    {
      cerr << "Cannot find data file(s) in directory: " << directory << endl;
      return false;
    }

  // These are the two key strings to find TBNames string 
  const string startTBNStr = ":Start of TBNames string:";
  const string endTBNStr   = ":End of TBNames string:";
  
  // I've to do three things: 
  // - build a string of detector TBnames to be stored on Run header
  // - build a map of detectors
  // - store on Run HEader DST version
  list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
  list<CsDetector*>::iterator id;
  _allTBNames = startTBNStr;
  for( id=dets.begin(); id!=dets.end(); id++ ) 
    {
      _allTBNames += (*id)->GetTBName();
    }
  _allTBNames += endTBNStr;

  // Please TM is needed
  CsRegistry reg;
  reg.EOJRegistration(this);
  
  isInit = true;
  return(true);   
}

bool CsFileStore::scan()
{   
  if(!isInit) 
    {
      cerr << "CsFileStore::scan: CsFileStore should be init first" << endl;
      return false;
    }

  if(isScanned) 
    {
      cerr << "CsFileStore::scan: Multiple scan ignored" << endl;
      return(true);
    }


  fileRun = new CsFileRun(runNumber);

  bool initSuccess = fileRun->Init(file_names, chunkName);

  if(!initSuccess)
    {
      cerr << "CsFileStore::scan: cannot initialize CsFileRun" << endl;
      return false;
    }

  string tbNames = fileRun->GetRunLog();
  if(!TestTBNameString(tbNames,_allTBNames,1))
    {
      cerr << "CsFileStore::scan: false test of the TBNames string."
	   << endl;
      return false;
    }
 
  isScanned = true;
  return true;
}

bool CsFileStore::next()
{
  // next() positions the iterator
  // on the next event (the first
  // time also...) returning true; 
  // if it returns false, the last 
  // event was the previous one.
  if(!isScanned) 
    {
      cerr << "CsFileStore::next no scan performed" << endl;
      cerr << "CsFileStore::next abort and exit..." << endl;
      return(false); // return or exit?
    }

  while(!fileRun->NextEvent())
    {
      if(!fileRun->SelectNextChunk(1,chunkName))
	return false;
    }
  return true;
}

bool CsFileStore::downloadDST(CsRecoEvent* event, int version)
{
  if(!dst())
    {
      cerr << "CsFileStore::downloadDST: no DST data for the current event." << endl;
      return false;
    }
  CsDstEvent* dstEvent = fileRun->GetDstEvent();
  if(!dstEvent)
    {
      cerr << "CsOraStore::downloadDST: DST downloading error." << endl;
      return false;
    }
  dstEvent->GetEvent(*event);
  return true;
}

bool CsFileStore::downloadTBNamesString(string &log)
{
  if(!fileRun) return false;
  log = fileRun->GetRunLog();
  return true;
}

bool CsFileStore::ScanDirectory(std::string& dir, int runNum, std::string& chunk)
{
  DIR* pdir = rfio_opendir(const_cast<char*>(dir.c_str()));
  if(!pdir)
    {
      cerr << "CsFileStore::ScanDirectory:  Cannot find directory: " << dir << endl;
      return false;
    }
  struct dirent* direntry = NULL;
  string mask = string("-")+itostr(runNumber)+".dst"+itostr(dstSlot);
  //  cout << "Mask: " << mask << endl;
  //  cout << "Dir: " << dir << endl;
  while((direntry = rfio_readdir(pdir))!= NULL)
    {
      string fileName = dir;
      fileName += direntry->d_name;
      // cout << "File name: " << fileName << endl;
      if(chunk == "")
	{
	  if(strstr(direntry->d_name,"cdr") && strstr(direntry->d_name,mask.c_str()))
	    {
	      file_names.push_back(fileName);
	      cout << fileName << endl;
	    }
	}
      else if(strstr(direntry->d_name,chunk.c_str()) && strstr(direntry->d_name,mask.c_str()))
	{
	  file_names.push_back(fileName);
	  cout << fileName << endl;
	  break;
	}
    }
  rfio_closedir(pdir);
  if(file_names.size() == 0) return false;
  return true;
}

#endif //#if USE_FileStore
