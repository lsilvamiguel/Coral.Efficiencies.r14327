#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include <cctype>

#include "CsOraStore.h"
#include "CsOpt.h"
#include "CsRegistry.h"
#include "CsGeom.h"

using namespace std;

vector<string>  CsOraStore::_extVarNames(3);
vector<float32> CsOraStore::_extVars(3);

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

bool CsOraStore::GetDstSelectionCriteria(CsOpt* myOptions)
{
  dstSelectionCriteria = "";
  string key = "DST event";
  string tag = "selection criteria";
  string criteria;
  while(myOptions->getOpt(key, tag, criteria))
    {
      if(criteria.size() == 0) break;
      if(dstSelectionCriteria != "") dstSelectionCriteria += string(" ");
      tag += string(" ") + criteria;
      dstSelectionCriteria += criteria;
      criteria = "";
    }
  if(dstSelectionCriteria != "" && isDstProduction)
    {
      cerr << "CsOraStore::GetDstSelectionCriteria: "
	   << "we should not have selection criteria for DST production."
	   << endl;
      return false;
    }
  else if(dstSelectionCriteria != "")
    {
      cout << "DST event selection criteria: \n"
	   << dstSelectionCriteria << "\n"
	   << "has been selected." << endl;
    }
  return true;
}

bool CsOraStore::GetRawSelectionCriteria(CsOpt* myOptions)
{
  rawSelectionCriteria = "";
  string key = "Raw event";
  string tag = "selection criteria";
  string criteria;
  while(myOptions->getOpt(key, tag, criteria))
    {
      if(criteria.size() == 0) break;
      if(rawSelectionCriteria != "") rawSelectionCriteria += string(" ");
      tag += string(" ") + criteria;
      rawSelectionCriteria += criteria;
      criteria = "";
    }
  if(rawSelectionCriteria != "" && isDstProduction)
    {
      cerr << "CsOraStore::GetRawSelectionCriteria: "
	   << "we should not have selection criteria for DST production."
	   << endl;
      return false;
    }
  else if(rawSelectionCriteria != "")
    {
      cout << "Raw event selection criteria: \n"
	   << rawSelectionCriteria << "\n"
	   << "has been selected." << endl;
    }
  return true;
}

bool CsOraStore::ParseDstCriteriaString()
{
  if(dstSelectionCriteria == "") return true;
  string criteria = dstSelectionCriteria;
  const string delims("!@#$%^&*()+`-=\t{}|[]\\;':\",./<>? \n");
  string tmp = "";
  string::size_type endIdx = 0;
  while(1)
    {
      string::size_type begIdx = criteria.find("EV.",endIdx);
      tmp += criteria.substr(endIdx,begIdx-endIdx);
      if(begIdx == string::npos) 
	{
	  break;
	}
      begIdx += 3; // skip "EV." substring
      endIdx = criteria.find_first_of(delims,begIdx);
      string variableName = criteria.substr(begIdx,endIdx-begIdx);
      unsigned n = 0;
      for(n = 0; n < _extVarNames.size(); n++)
	{
	  if(variableName == _extVarNames[n]) break;
	}
      if(n == _extVarNames.size())
	{
	  cerr << "DST selection criteria fault: the external variable: "
	       << variableName << " has not been found." << endl;
	  cerr << "Correct variable names are:\n";
	  for(unsigned i = 0; i < _extVarNames.size(); i++)
	    {
	      cerr << "\t" << _extVarNames[i] << "\n";
	    }
	  cerr << endl;
	  return false;
	}
      tmp += string("VALUE") + utostr(n+1);
    }
  cout << "DST selection criteria after parsing: " << tmp << endl;
  dstSelectionCriteria = tmp;
  return true;
}


CsOraStore* CsOraStore::Instance()
{
  static  CsOraStore* _instance = new CsOraStore;
  return _instance;
}

CsOraStore::CsOraStore() 
{
  dataType = DATA_TYPE_ERROR;
  dstSelectionCriteria = "";
  rawSelectionCriteria = "";
  chunkName = "";
  runNumber = 0;
  dstSlot = 0;
  isDstProduction = false;
  oraRun = NULL;
  isInit = false;
  isScanned = false;

  _preuploadFunc = NULL;
}

bool CsOraStore::init()
{
  if(isInit) 
    {
      cerr << "CsOraStore::init init already performed " << endl;
      return(true);
    }

  CsOpt*  myOptions = CsOpt::Instance();
  
  string dType = "";

  if(!myOptions->getOpt("Data", "type", dType )) 
    {
      cerr << "CsOraStore::init: Data type not selected..." << endl;
      return(false);
    }
  else
    {
      for(unsigned p = 0; p < dType.size(); p++) dType[p] = tolower(dType[p]);
      if (dType == "raw")
	{
	  dataType = RAW_DATA_TYPE;
	}
      else if(dType == "dst")
	{
	  dataType = DST_DATA_TYPE;
	}
      else
	{
	  cerr << "CsOraStore::init: incorrect data type selected: "
	       << dType << endl;
	  return false;
	}
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
      cerr << "CsOraStore::init: No any chunk selected." << endl;
      return false;
    }
  else if(chunkName != "") // chunk has been selected
    {
      int chunkNum = 0;
      if(sscanf(chunkName.c_str(),"%05d-%d",&chunkNum, &runNumber) != 2)
	{
	  cerr << "CsOraStore::init: Incorrect chunk name: " <<
	    chunkName << endl;
	  return false;
	}
    }
  dstVersion = 0;

  bool saveTBNamesOnly = false;

  if(CsOpt::Instance()->getOpt( "", "store TBNames on DST only" ) )
    {
      isDstProduction = true;
      saveTBNamesOnly = true;
    }

  if(!saveTBNamesOnly && myOptions->getOpt("make", "DST", dstVersion) && dataType == RAW_DATA_TYPE)
    {
      if(chunkName == "")
	{
	  cerr << "CsOraStore::init: chunk name for DST production is not selected."
	       << endl;
	  return false;
	}
      if(dstVersion == 0 || dstVersion > LAST_DST_VERSION)
	{
	  cerr << "CsOraStore::init: Incorrect DST version " << dstVersion 
	       << " has been selected." << endl;
	  return false;
	}
      cout << "DST" << dstVersion << " production in Oracle DB has been selected."
	   << endl;
      isDstProduction = true;
    }
  else if(dstVersion && dataType != RAW_DATA_TYPE)
    {
      cerr << "CsOraStore::init: For DST production, data type must be RAW!" << endl;
      return false;
    }

  if(myOptions->getOpt("Dst select", "slot", dstSlot))
    {
      if(dstSlot < 1 || dstSlot > 3)
	{
	  cerr << "CsOraStore::init: Incorrect DST slot " << dstSlot
	       << " has been selected." << endl;
	  return false;
	}
      cout << "DST slot:                   " << dstSlot << endl;
    }
  else if(dataType == DST_DATA_TYPE || isDstProduction)
    {
      cerr << "CsOraStore::init: You have to select DST slot number." << endl;
      return false;
    }
  else if(dataType == RAW_DATA_TYPE && !isDstProduction)
    {
      dstSlot = 1; // default slot number in case of reading RAW
    }
  // Read preselection criteria from coral.option file:
  if(!GetRawSelectionCriteria(myOptions))
    return false;
  if(!GetDstSelectionCriteria(myOptions))
    return false;
  
  if(dataType == DST_DATA_TYPE || isDstProduction)
    {
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
    }

  // Please TM is needed
  CsRegistry reg;
  reg.EOJRegistration(this);
  
  isInit = true;
  return(true);   
}

bool CsOraStore::scan()
{   
  if(!isInit) 
    {
      cerr << "CsOraStore::scan: CsOraStore should be init first" << endl;
      return false;
    }

  if(isScanned) 
    {
      cerr << "CsOraStore::scan: Multiple scan ignored" << endl;
      return(true);
    }


  oraRun = new CsOraRun(runNumber);

  bool saveTBNamesOnly = false;

  if(CsOpt::Instance()->getOpt( "", "store TBNames on DST only" ) )
    {
      saveTBNamesOnly = true;
    }

  bool initSuccess = false;
  if(isDstProduction)
    {
      if(saveTBNamesOnly)
	initSuccess = oraRun->Init(IS_RUN_LOG_UPDATE);
      else if(dstVersion == 4)
	initSuccess = oraRun->Init(IS_DST_PRODUCTION_04);
    }
  else if(dataType == DST_DATA_TYPE)
    {
      initSuccess = oraRun->Init(IS_DST_READING);
    }
  else if(dataType == RAW_DATA_TYPE)
    {
      initSuccess = oraRun->Init(IS_RAW_READING);
    }
  else
    {
      cerr << "CsOraStore::scan: unknown data type." << endl;
      return false;
    }
  if(!initSuccess)
    {
      cerr << "CsOraStore::scan: cannot initialize CsOraRun" << endl;
      return false;
    }
  else if(saveTBNamesOnly)
    return true;
  if(dataType == DST_DATA_TYPE || isDstProduction)
    {
      string tbNames = oraRun->GetRunLog();
      if(!TestTBNameString(tbNames,_allTBNames,dstSlot))
	{
	  cerr << "CsOraStore::scan: false test of the TBNames string."
	       << endl;
	  return false;
	}
    }
  else
    {
      oraRun->dstVersion = 1;
    }
#ifdef CS_DUMP

  oraRun->Print(cout);

#endif
 
  isScanned = true;
  return oraRun->SelectNextChunk(dstSlot-1,chunkName);
}

bool CsOraStore::next() 
{
  // next() positions the iterator
  // on the next event (the first
  // time also...) returning true; 
  // if it returns false, the last 
  // event was the previous one.
  if(!isScanned) 
    {
      cerr << "CsOraStore::next no scan performed" << endl;
      cerr << "CsOraStore::next abort and exit..." << endl;
      abortAndExit();
      return(false); // return or exit?
    }

  for(unsigned i = 0; i < _extVars.size(); i++)
    _extVars[i] = 0.;

  while(!oraRun->NextEvent())
    {
      if(!oraRun->SelectNextChunk(dstSlot-1,chunkName))
	return false;
    }
  return true;
}

bool CsOraStore::downloadDST(CsRecoEvent* event, int version)
{
  if(!dst())
    {
      cerr << "CsOraStore::downloadDST: no DST data for the current event." << endl;
      return false;
    }
  CsDstEvent* dstEvent = oraRun->GetDstEvent();
  if(!dstEvent)
    {
      cerr << "CsOraStore::downloadDST: DST downloading error." << endl;
      return false;
    }
  dstEvent->GetEvent(*event);
  for(int i = 0; i < 3; i++)
    {
      _extVars[i] = dstEvent->GetExternalVariable(i);
    }
  delete dstEvent;
  return true;
}

bool CsOraStore::uploadDST(CsRecoEvent* event, int fatness, int version)
{
  if(!isDstProduction)
    {
      cerr << "CsOraStore::uploadDST: incorrect usage of the method." << endl;
      return false;
    }

  CsDstEvent* dstEvent = NULL;

  switch(version)
    {
    case 4:
      dstEvent = new CsDstEvent_04(fatness);
      break;
    default:
      cerr << "CsOraStore::uploadDST: incorrect version number: " << version << endl;
      return false;
    }

  dstEvent->Init(*event);

  if(_preuploadFunc)
    {
      _preuploadFunc();
      for(unsigned i = 0; i < _extVars.size(); i++)
	{
	  dstEvent->SetExternalVariable(i,_extVars[i]);
	}
    }

  if(!oraRun->SaveDstEvent(dstEvent))
    {
      delete dstEvent;
      return false;
    }
  delete dstEvent;
  return true;
}

bool CsOraStore::downloadTBNamesString(string &log)
{
  if(!oraRun) return false;
  log = oraRun->GetRunLog();
  return true;
}

bool CsOraStore::uploadTBNamesString(string& log)
{
  if(log.find(":Start of DST version string:") != std::string::npos)
    {
      cout << "DST version string is not necessary anymore. Skip update." << endl;
      return true;
    }
  if(isDstProduction)
    {
      return oraRun->UpdateRunLog(log,dstSlot-1);
    }
  return true;
}

void CsOraStore::saveAndExit()
{
  if(!isDstProduction) return;
  if(!oraRun->CloseDst())
    {
      cerr << "CsOraStore::saveAndExit: FATAL ERROR." << endl;
      throw;
    }
}

bool CsOraStore::end()
{
  try
    {
      saveAndExit();
    }
  catch(...)
    {
      cerr << "CsOraStore::end: saveAndExit ERROR." << endl;
      abortAndExit();
      return false;
    }
  return true;
}

void CsOraStore::abortAndExit()
{
  if(!isDstProduction) return;
  oraRun->Abort();
}

void CsOraStore::setExtVarName(int n, const string& name)
{
  if(n < 0 || n >=3 ) // 3 ext. vars only for DST4
    {
      cerr << "CsOraStore::setExtVarName: maximum number of external variable is 2"
	   << endl;
      return;
    }
  _extVarNames[n] = name;
  if( _extVarNames[n].size() >= 12 )
    {
      _extVarNames[n].erase(11); // Maximum name length is 12 characters
      _extVarNames[n][11] = 0;
    }
}

CsOraStore::~CsOraStore()
{
  delete oraRun;
}

