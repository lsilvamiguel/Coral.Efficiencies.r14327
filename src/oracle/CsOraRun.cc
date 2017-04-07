#include <CsTime.h>
#include <vector>
#include <algorithm>
#include "FlatFile.h"
#include "CsOraRun.h"
#include "RunInfo.h"
#include "CsOraStore.h"
#include "OracleErrorCode.h"
#include "CsOpt.h"
#include "runflag.h"

using namespace std;
using namespace oracle::occi;

static bool UpdateTBNamesString(string& runLog, const string& lg, int dstSlot) //dstSlot = 0..2
{
  const string startTBNStr = ":Start of TBNames string:";
  const string endTBNStr   = ":End of TBNames string:";
  string tmp = runLog;
  runLog = "";
  std::string::size_type idx = 0;
  for(int i = 0; i < 3; i++)
    {
      if(i ==  dstSlot)
	runLog += lg;
      std::string::size_type idx_start = tmp.find( startTBNStr , idx );
      if(idx_start == std::string::npos)
	{
	  if(i >= dstSlot)
	    break;
	  else if(i < dstSlot)
	    {
	      runLog += startTBNStr+endTBNStr;
	      continue;
	    }
	}
      std::string::size_type idx_end = tmp.find(endTBNStr, idx_start);
      if(idx_end == std::string::npos)
	{
	  cerr << "UpdateTBNamesString(): TBNames string error: end-tag not found." << endl;
	  return false;
	}
      idx_end += endTBNStr.length();
      if(i !=  dstSlot)
	runLog += tmp.substr(idx_start, idx_end-idx_start);
      idx = idx_end;
    }
  return true;
}

CsOraRun::CsOraRun(int runnum)
{
  session = NULL;
  isDstProduction = isDstReading = isRawReading = isRunLogUpdate = false;
  runInfo = new OraRunInfo;
  runNumber = runnum;
  if( runNumber > 32832 ) { runGT32832 = true; } else { runGT32832 = false; }
  runGT46692 = (runNumber > 46692); //used to limit # of rich1probs to be read from DSTs
  try
    {
      RunInfo run(runNumber);
      runYear = run.GetRunYear();
      runPeriod = run.GetRunPeriod();
    }
  catch(int nError)
    {
      cerr << "CsOraRun::CsOraRun: run " << runNumber << " is not found in DB." << endl;
      throw ORACLE_ERROR;
    }
  curChunkInfo = NULL;
  runInfo->sorFileId = 0;
}

CsOraRun::~CsOraRun()
{
  for(unsigned i = 0; i < runInfo->chunkInfo.size(); i++)
    {
      for(int j = 0; j < 3; j++)
	{
	  delete runInfo->chunkInfo[i].dstChunkInfo[j];
	}
    }
  delete runInfo;
}

bool CsOraRun::Init(int todo)
{
  switch(todo)
    {
    case IS_DST_PRODUCTION_04:
      dstVersion = 4;
      isDstProduction = true;
      break;
    case IS_RAW_READING:
      isRawReading = true;
      break;
    case IS_DST_READING:
      isDstReading = true;
      break;
    case IS_RUN_LOG_UPDATE:
      isRunLogUpdate = true;
      break;
    default:
      cerr << "CsOraRun::Init: incorrect parameter: " << todo << endl;
      return false;
    }
  session = CsOraSession::instance();
  if(!session->init(runNumber,isDstProduction||isRunLogUpdate))
    return false;
  
  
  string query = string("select T_MIN, T_MAX, LOG_INFO, SOR_FILE_ID "
			"from RUNS where RUN_NUMBER=:1");
  ResultSet* result = session->executeQuery(query,runNumber);
  if(!result || !result->next())
    {
      cerr << "CsOraRun::Init: Run number: " << runNumber
	   << " not found in database." << endl;
      return false;
    }
  start = result->getUInt(1);
  stop  = result->getUInt(2);
  runInfo->sorFileId = result->getInt(4);
  { // to call ~Clob before the session will be stopped
    Clob clob = result->getClob(3);
    if(clob.isNull() || clob.length()<1)
      runLog = "";
    else
      {
	clob.open(OCCI_LOB_READONLY);
	vector<unsigned char>	buffer;
	int	len = clob.length();
	buffer.resize( len+1 );
	clob.read(len, &buffer[0], len);
	buffer[len] = 0;
	runLog = (const char *)(&buffer[0]);
	clob.close();
      }
  }
  if(isRunLogUpdate) return true;
  query = string("select file_id, file_dir, file_name, file_size, file_type from file_maps where run_number=:1");
  result = session->executeQuery(query,runNumber);
  if(!result || !result->next())
    {
      cerr << "CsOraRun::Init: No any chunks found for run number " << runNumber << endl;
      return false;
    }
  do
    {
      OraRawChunkInfo rawChunk;
      rawChunk.fileId = result->getUInt(1);
      rawChunk.dirName = result->getString(2);
      rawChunk.fileName = result->getString(3);
      rawChunk.fileSize = result->getUInt(4);
      string fileType = result->getString(5);
      if(fileType == "RAW")
	{
	  rawChunk.dstChunkInfo[0] = NULL;
	  rawChunk.dstChunkInfo[1] = NULL;
	  rawChunk.dstChunkInfo[2] = NULL;
	  runInfo->chunkInfo.push_back(rawChunk);
	}
      else if(fileType == "SOR")
	{
	  cout << "SOR file: " << rawChunk.dirName << rawChunk.fileName << endl;
	}
    } while(result->next());
  cout << "Total number of chunks:     " << runInfo->chunkInfo.size() << endl;
  if(isDstProduction || isDstReading)
    {
      for(unsigned i = 0; i < runInfo->chunkInfo.size(); i++)
	{
	  OraRawChunkInfo& rawChunk = runInfo->chunkInfo[i];
	  query = string("select file_id, file_dir, file_name, file_size, dst_version, dst_type_id,"
			 " value1_is, value2_is, value3_is from dst_files where rawev_file_id=:1");
	  result = session->executeQuery(query,rawChunk.fileId);
	  if(!result || !result->next()) continue;
	  do
	    {
	      unsigned file_id     = result->getUInt(1);
	      string   file_dir    = result->getString(2);
	      string   file_name   = result->getString(3);
	      unsigned file_size   = result->getUInt(4);
	      unsigned dst_slot    = result->getUInt(5)-1;
	      unsigned dst_version = result->getUInt(6);
	      string   extVar0     = result->getString(7);
	      string   extVar1     = result->getString(8);
	      string   extVar2     = result->getString(9);
	      rawChunk.dstChunkInfo[dst_slot] = new OraDstChunkInfo;
	      OraDstChunkInfo& dstInfo = *rawChunk.dstChunkInfo[dst_slot];
	      dstInfo.fileId     = file_id;
	      dstInfo.dirName    = file_dir;
	      dstInfo.fileName   = file_name;
	      dstInfo.fileSize   = file_size;
	      dstInfo.dstVersion = dst_version;
	      dstInfo.extVarName.push_back(extVar0);
	      dstInfo.extVarName.push_back(extVar1);
	      dstInfo.extVarName.push_back(extVar2);
	    } while(result->next());
	}
    }
  return session->stopSession();
}

bool CsOraRun::SelectNextChunk(int slotNum, const std::string& chunkName)
{
  slotNumber = slotNum;
  if(isRunLogUpdate) return true;
  OraRawChunkInfo* chunkInfo = NULL;
  if(rawFile.isOpen())
    {
      rawFile.Close();
    }
  if(dstFile.isOpen())
    {
      dstFile.Close();
    }
  if(chunkName == "") // full run reading
	  {
		  if(!curChunkInfo) // select first chunk
			  {
				  if(runInfo->chunkInfo.size() == 0)
					  {
						  cerr << "CsOraRun::SelectNextChunk: no any chunk for run " << runNumber << endl;
						  return false;
					  }
				  chunkInfo = &runInfo->chunkInfo[0];
			  }
		  else // select next chunk
			  {
	  for(unsigned i = 0; i < runInfo->chunkInfo.size()-1; i++)
	    {
	      chunkInfo = &runInfo->chunkInfo[i];
	      if(curChunkInfo == chunkInfo)
		{
		  // Free previous chunk memory:
		  chunkInfo->eventInfo.clear();
		  delete chunkInfo->dstChunkInfo[slotNumber];
		  chunkInfo->dstChunkInfo[slotNumber] = NULL;

		  // Open new chunk:
		  chunkInfo = &runInfo->chunkInfo[i+1];

		  // Non produced DST chunks are skipped  >> SILENTLY <<
		  if(isDstReading && !chunkInfo->dstChunkInfo[slotNumber]) {
			  curChunkInfo = chunkInfo;
			  return SelectNextChunk(slotNum /* chunkName="" */);
		  }
		  break;
		}
	      chunkInfo = NULL;
	    }
	}
    }
  else // select chunk chunkName
    {
      if(curChunkInfo)
	return false;
      for(unsigned i = 0; i < runInfo->chunkInfo.size(); i++)
	{
	  chunkInfo = &runInfo->chunkInfo[i];
	  if(chunkInfo->fileName.find(chunkName) != std::string::npos)
	    {
	      break;
	    }
	  chunkInfo = NULL;
	}
      if(!chunkInfo)
	{
	  cerr << "CsOraRun::SelectNextChunk: Chunk " << chunkName << " is not found in DB." << endl;
	  return false;
	}
    }
  if(!chunkInfo)
    {
      cout << "CsOraRun::SelectNextChunk: No more chunks for run " << runNumber << endl;
      return false;
    }
  cout << "Selected chunk:             " << chunkInfo->fileName << endl;

  // chunkName == "" means reading a whole chain; this error should not occur there
  if(isDstReading && !chunkInfo->dstChunkInfo[slotNumber] && chunkName != "" )
    {
      cerr << "CsOraRun::SelectNextChunk: DST for chunk " << chunkName << " slot number " << slotNumber+1 
	   << " is not found in DB." << endl;
      return false;
    }
  if(isDstProduction && chunkInfo->dstChunkInfo[slotNumber])
    {
      cerr << "CsOraRun::SelectNextChunk: DST for chunk " << chunkName << " slot number " << slotNumber+1
	   << " is exist in DB already. \nYou have to delete it before start new production." << endl;
      return false;
    }
  if(!session->startSession())
    return false;
  string query = string("select event_number,burst,event_in_burst,trigger_mask,time_sec,time_usec,"
			"error_code,event_size,rawev_file_offset from event_headers where rawev_file_id=:1");

  if (CsOraStore::Instance()->RawSelectionCriteria() != "" ) {
    query += " and " + CsOraStore::Instance()->RawSelectionCriteria();
  }

  Statement* stmt;
  ResultSet* result;
  try {
	  stmt = session->conn()->createStatement();
	  stmt->setSQL(query);
	  stmt->setInt(1,chunkInfo->fileId);
	  stmt->setPrefetchRowCount(1000);
	  result = stmt->executeQuery();
  } catch(SQLException& e) {
	  cerr << "CsOraRun::SelectNextChunk: the Raw Selection Criteria \"" 
			 << CsOraStore::Instance()->RawSelectionCriteria()
			 << "\" turned out to be uncorrect - ERROR "
			 << e.getMessage();
	  return false;
  } 

   if(!result) {
     cerr << "CsOraRun::SelectNextChunk: No any RAW event found for the chunk " << chunkName << endl;
     return false;
   }
   if ( !result->next()) {
     if ( CsOraStore::Instance()->RawSelectionCriteria() != "" ) {
       cout << "CsOraRun::SelectNextChunk: No any RAW event matching criteria" << endl;
       curChunkInfo = chunkInfo;
       stmt->closeResultSet(result);
       session->conn()->terminateStatement(stmt);
       return SelectNextChunk(slotNum); //try with the next chunk
     } else {
       cerr << "CsOraRun::SelectNextChunk: No any RAW event found for the chunk " << chunkName << endl;
       return false;
     }
   }

  do
    {
      OraRawEventInfo eventInfo;
      eventInfo.eventNumber  = result->getUInt(1);
      eventInfo.burstNumber  = result->getUInt(2);
      eventInfo.eventInBurst = result->getUInt(3);
      eventInfo.triggerMask  = result->getUInt(4);
      eventInfo.timeInSec    = result->getUInt(5);
      eventInfo.timeInUsec   = result->getUInt(6);
      eventInfo.errorCode    = result->getUInt(7);
      eventInfo.eventSize    = result->getUInt(8);
      eventInfo.fileOffset   = result->getUInt(9);
      eventInfo.eventToRead  = true;
      chunkInfo->eventInfo.push_back(eventInfo);
    } while(result->next());
  stmt->closeResultSet(result);
  session->conn()->terminateStatement(stmt);
  cout << "Total number of RAW events: " << chunkInfo->eventInfo.size() << endl;
  if(isDstReading)
    {
      query = string("select event_number, dst_size, file_offset, trigger_mask, value1, value2, value3"
		     " from dst where file_id=:1");
      Statement* stmt = session->conn()->createStatement();
      stmt->setSQL(query);
      stmt->setInt(1,chunkInfo->dstChunkInfo[slotNumber]->fileId);
      stmt->setPrefetchRowCount(1000);
      ResultSet* result = stmt->executeQuery();
      if(!result || !result->next())
	{
	  cerr << "CsOraRun::SelectNextChunk: No any DST event found for the chunk " << chunkName << endl;
	  return false;
	}
      do
	{
	  OraDstEventInfo eventInfo;
	  eventInfo.eventNumber = result->getUInt(1);
	  eventInfo.eventSize   = result->getUInt(2);
	  eventInfo.eventOffset = result->getUInt(3);
	  eventInfo.triggerMask = result->getUInt(4);
	  eventInfo.extVar.push_back(result->getDouble(5));
	  eventInfo.extVar.push_back(result->getDouble(6));
	  eventInfo.extVar.push_back(result->getDouble(7));
	  chunkInfo->dstChunkInfo[slotNumber]->eventInfo.push_back(eventInfo);
	} while(result->next());
      stmt->closeResultSet(result);
      session->conn()->terminateStatement(stmt);
      cout << "Total number of DST events: " << chunkInfo->dstChunkInfo[slotNumber]->eventInfo.size() << endl;
    }
  string SkippedEventFileName;
  CsOpt*  myOptions = CsOpt::Instance();
  if(myOptions && myOptions->getOpt("Events to", "read", SkippedEventFileName))
    {
      FILE* file = fopen(SkippedEventFileName.c_str(),"r");
      vector<unsigned long> eventsToRead;
      if(file)
	{
	  char line[1024];
	  while(fgets(line,1023,file))
	    {
	      unsigned long ev_num = 0;
	      if(sscanf(line,"%lu",&ev_num) != 1)
		continue;
	      eventsToRead.push_back(ev_num);
	    }
	  fclose(file);
	  if(eventsToRead.size() != 0)
	    {
	      for(unsigned i = 0; i < chunkInfo->eventInfo.size(); i++)
		{
		  vector<unsigned long>::iterator pos = 
		    find(eventsToRead.begin(),eventsToRead.end(),chunkInfo->eventInfo[i].eventNumber);
		  if(pos != eventsToRead.end())
		    {
		      chunkInfo->eventInfo[i].eventToRead = true;
		      eventsToRead.erase(pos);
		    }
		  else
		    chunkInfo->eventInfo[i].eventToRead = false;
		}
	    }
	}
      else
	{
	  cerr << "******************************************************************************\n";
	  cerr << "CsOraRun::SelectNextChunk: Cannot open file: " << SkippedEventFileName << "\nNo skipped events.\n";
	  cerr << "******************************************************************************\n" << endl;
	}
    }
  delete CsDstRun::dstChunk;
  CsDstRun::dstChunk = NULL;
  session->stopSession();
  if(isDstReading)
    {
      dstVersion = chunkInfo->dstChunkInfo[slotNumber]->dstVersion;
      string curChunkName = chunkInfo->fileName.substr(3,11);
      switch(dstVersion)
	{
	case 4:
	  CsDstRun::dstChunk = new CsDstChunk_04;
	  CsDstRun::dstChunk->Init(curChunkName.c_str());
	  break;
	default:
	  cerr << "CsOraRun::SelectNextChunk: Invalid DST version: " << dstVersion << endl;
	  return false;
	}
      dstFile.Open(chunkInfo->dstChunkInfo[slotNumber]->dirName,chunkInfo->dstChunkInfo[slotNumber]->fileName);
      if(!dstFile.isOpen())
	{
		if ( chunkName == "" ) {
			curChunkInfo = chunkInfo;
			cout << "CsOraRun::SelectNextChunk. Skipping " 
				  << chunkInfo->dstChunkInfo[slotNumber]->fileName << endl;
			return SelectNextChunk(slotNumber /* chunkName="" */);
		}
	  return false;
	}
      if(!CsDstRun::dstChunk->Load(dstFile))
	{
		if ( chunkName == "" ) {
			curChunkInfo = chunkInfo;
			cout << "CsOraRun::SelectNextChunk. Skipping "
				  << chunkInfo->dstChunkInfo[slotNumber]->fileName 
				  << " for CsDstRun::dstChunk->Load failed" << endl;
			return SelectNextChunk(slotNumber /* chunkName="" */);
		}
	  return false;
	}
    }
  else if(isDstProduction)
    {
      switch(dstVersion)
	{
	case 4:
	  CsDstRun::dstChunk = new CsDstChunk_04;
	  CsDstRun::dstChunk->Init(chunkName.c_str());
	  break;
	default:
	  cerr << "CsOraRun::SelectNextChunk: Invalid DST version: " << dstVersion << endl;
	  return false;
	}
      char XXX[10];
      sprintf(XXX,"%02d",int(((runNumber%100)/10)*10));
      chunkInfo->dstChunkInfo[slotNumber] = new OraDstChunkInfo;
      chunkInfo->dstChunkInfo[slotNumber]->dstVersion = dstVersion;
      chunkInfo->dstChunkInfo[slotNumber]->fileSize   = 0;
      chunkInfo->dstChunkInfo[slotNumber]->fileId     = 0;
      chunkInfo->dstChunkInfo[slotNumber]->fileName   = string("cdr") + chunkName + ".dst" + itostr(slotNumber+1);
      chunkInfo->dstChunkInfo[slotNumber]->dirName    = string("/castor/cern.ch/compass/data/") + itostr(runYear)
	+ "/oracle_dst/" + runPeriod + "/" + XXX + "/slot" + itostr(slotNumber+1) + "/";
      chunkInfo->dstChunkInfo[slotNumber]->extVarName.resize(3);
      if(!dstFile.Create(chunkInfo->dstChunkInfo[slotNumber]->dirName,chunkInfo->dstChunkInfo[slotNumber]->fileName))
	{
	  return false;
	}
      for(int i = 0; i < 3; i++)
	{
	  chunkInfo->dstChunkInfo[slotNumber]->extVarName[i] = CsOraStore::getExtVarName(i);
	  CsDstRun::dstChunk->SetExtVarName(i,CsOraStore::getExtVarName(i));
	}
      CsDstRun::dstChunk->SetTBNames(CsOraStore::Instance()->getTBNames());
      vector<uint8> buffer;
      if(!CsDstRun::dstChunk->Save(dstFile))
	{
	  return false;
	}
    }
  eventMap.clear();
  for(unsigned i = 0; i < chunkInfo->eventInfo.size(); i++)
    {
      eventMap[chunkInfo->eventInfo[i].eventNumber].rawEventInfo = &chunkInfo->eventInfo[i];
    }
  if(isDstReading)
    {
      for(unsigned i = 0; i < chunkInfo->dstChunkInfo[slotNumber]->eventInfo.size(); i++)
	{
	  uint32 eventNumber = chunkInfo->dstChunkInfo[slotNumber]->eventInfo[i].eventNumber;
	  if(eventMap[eventNumber].rawEventInfo == NULL)
	    {
	      cerr << "CsOraRun::SelectNextChunk: ORACLE DB ERROR: No RAW event number " << eventNumber << " in database." << endl;
	      return false;
	    }
	  eventMap[eventNumber].dstEventInfo = &chunkInfo->dstChunkInfo[slotNumber]->eventInfo[i];
	}
    }
  isNewChunk = true;
  curChunkInfo = chunkInfo;
  return true;
}

bool CsOraRun::NextEvent()
{
  if(isNewChunk)
    {
      isNewChunk = false;
      eventIt = eventMap.begin();
    }
  else
    {
      eventIt++;
    }
  while(eventIt != eventMap.end())
    {
      OraRawEventInfo* eventInfo = eventIt->second.rawEventInfo;
      if(eventInfo->eventToRead)
	break;
      eventIt++;      
    }
  if(eventIt == eventMap.end())
    return false;
  return true;
}

uint8* CsOraRun::GetRawEvent()
{
  if(eventIt == eventMap.end())
    return NULL;
  OraRawEventInfo* eventInfo = eventIt->second.rawEventInfo;
  if(!eventInfo) // Debug check, should be never.
    { // should be never
      cerr << "CsOraRun::GetRawEvent: no DB information about the event." << endl;
      throw;
    }
  if(rawBuffer.size() < eventInfo->eventSize) rawBuffer.resize(eventInfo->eventSize);
  if(!rawFile.isOpen())
    {
      string rawFileDir  = curChunkInfo->dirName;
      string rawFileName = curChunkInfo->fileName;
      if(rawFileDir[rawFileDir.length()-1] != '/')
	rawFileDir += "/";
      string fname = rawFileDir + rawFileName;
      cout << "Open read-only raw file: " << fname << endl;
      rawFile.Open(rawFileDir, rawFileName);
      cout << "done." << endl;
    }
  if(!rawFile.Read(&rawBuffer[0], eventInfo->eventSize, eventInfo->fileOffset))
    {
      throw;
    }
  return &rawBuffer[0];
}

CsDstEvent* CsOraRun::GetDstEvent()
{
  if(!dstFile.isOpen()) return NULL;
  CsDstEvent* event = NULL;
  switch(dstVersion)
    {
    case 4:
      event = new CsDstEvent_04(0); // fatness is not important here
      break;
    default:
      cerr << "CsOraRun::GetDstEvent: Incorrect DST version: " << dstVersion << endl;
      throw;
    }
  if(!event)
    {
      cerr << "CsOraRun::GetDstEvent: Memory allocation error. " << endl;
      throw;
    }
  dstFile.Seek(eventIt->second.dstEventInfo->eventOffset);
  if(!event->Load(dstFile))
    {
      cerr << "CsOraRun::GetDstEvent: Event loading error." << endl;
      throw;
    }
  event->SetRunNumber(runNumber);
  return event;
}

bool CsOraRun::SaveDstEvent(const CsDstEvent* event)
{
  if(!event)
    {
      cerr << "CsOraRun::SaveDstEvent: NULL pointer." << endl;
      return false;
    }

  vector<uint8> buffer;
  buffer.reserve(10000);
  event->GetBuffer(buffer);
  OraDstEventInfo eventInfo;
  eventInfo.eventNumber = event->GetEventNumber();
  eventInfo.eventSize   = buffer.size();
  eventInfo.eventOffset = dstFile.getFileOffset();
  eventInfo.triggerMask = event->GetTriggerMask();
  for(int i = 0; i < 3; i++)
    {
      eventInfo.extVar.push_back(event->GetExternalVariable(i));
    }

  if(!dstFile.Append(&buffer[0], buffer.size()))
    return false;

  curChunkInfo->dstChunkInfo[slotNumber]->eventInfo.push_back(eventInfo);
  curChunkInfo->dstChunkInfo[slotNumber]->fileSize += eventInfo.eventSize;

  return true;
}

bool CsOraRun::UpdateRunLog(const string& tbnames, int slot)
{
  slotNumber = slot;
  if(!isRunLogUpdate)
    {
      cerr << "CsOraRun::UpdateRunLog: Incorrect usage of method." << endl;
      return false;
    }
  if(tbnames.find(":Start of DST version string:") != std::string::npos)
    {
      cout << "DST version string is not necessary anymore. Skip update." << endl;
      return true;
    }
  if(!session->startSession())
    return false;
  try
    {
      string query = string("select log_info from runs where run_number=:1");
      session->stmt()->setSQL(query);
      session->stmt()->setUInt(1,runNumber);
      ResultSet* result = session->stmt()->executeQuery();
      if(!result || !result->next())
	return false;
      Clob clob = result->getClob(1);
      runLog = "";
      if(clob.isNull())
	{
	  Clob newClob(session->conn());
	  newClob.setEmpty();
	  query = string("update RUNS set LOG_INFO=:1 where RUN_NUMBER=:2");
	  session->stmt()->setSQL(query);
	  session->stmt()->setClob(1,newClob);
	  session->stmt()->setUInt(2,runNumber);
	  session->stmt()->executeUpdate();
	}
      else if(clob.length() > 0)
	{
	  clob.open(OCCI_LOB_READONLY);
	  vector<unsigned char>	buffer;
	  int	len = clob.length();
	  buffer.resize( len+1 );
	  clob.read(len, &buffer[0], len);
	  buffer[len] = 0;
	  runLog = (const char *)(&buffer[0]);
	  clob.close();
	}
      if(!UpdateTBNamesString(runLog,tbnames,slotNumber))
	return false;
      query = string("select log_info from runs where run_number=:1 for update");
	  session->stmt()->setSQL(query);
	  session->stmt()->setUInt(1,runNumber);
      result = session->stmt()->executeQuery();
      result->next();
      clob=result->getClob(1);
      clob.write( runLog.size(),(unsigned char *)runLog.c_str(), runLog.size());
    }
  catch(SQLException& e)
    {
      cerr << e.getMessage() << endl;
      return false;
    }
  cout << "##################################" << endl;
  cout << "String:\n";
  cout << tbnames << "\n";
  cout << "has been appended to runs.log_info" << endl;
  cout << "##################################" << endl;
  
  return true;
}


bool CsOraRun::CloseDst()
{
  if(isRunLogUpdate)
    {
      if(!session->commit())
	{
	  return false;
	}
      return session->stopSession();
    }
  if(!isDstProduction)
    return false;
try {
  if(!session->startSession()) return false;

  string query = string("select file_maps_seq.nextval from dual");
  ResultSet* result = session->executeQuery(query);

  if(!result || !result->next()) return false;

  Statement* stmt = session->conn()->createStatement();
  
  curChunkInfo->dstChunkInfo[slotNumber]->fileId = result->getNumber(1);
  query = string("insert into DST_FILES(FILE_ID,RAWEV_FILE_ID,FILE_DIR,FILE_NAME,FILE_SIZE,DST_VERSION,DST_TYPE_ID,VALUE1_IS,VALUE2_IS,VALUE3_IS,RUN_NUMBER) values(:1,:2,:3,:4,:5,:6,:7,:8,:9,:10,:11)");

  stmt->setSQL(query);
  stmt->setUInt(1,curChunkInfo->dstChunkInfo[slotNumber]->fileId);
  stmt->setUInt(2,curChunkInfo->fileId);
  stmt->setString(3,curChunkInfo->dstChunkInfo[slotNumber]->dirName);
  stmt->setString(4,curChunkInfo->dstChunkInfo[slotNumber]->fileName);
  stmt->setUInt(5,curChunkInfo->dstChunkInfo[slotNumber]->fileSize);
  stmt->setUInt(6,slotNumber+1);
  stmt->setUInt(7,curChunkInfo->dstChunkInfo[slotNumber]->dstVersion);
  stmt->setString(8,curChunkInfo->dstChunkInfo[slotNumber]->extVarName[0]);
  stmt->setString(9,curChunkInfo->dstChunkInfo[slotNumber]->extVarName[1]);
  stmt->setString(10,curChunkInfo->dstChunkInfo[slotNumber]->extVarName[2]);
  stmt->setUInt(11,runNumber);

  if( ! stmt->executeUpdate() ) return false;

  stmt->setSQL("insert into DST(file_id, event_number, dst_size, file_offset, trigger_mask, value1, value2, value3) values(:1,:2,:3,:4,:5,:6,:7,:8)");
  const int MAX_ITERATION_NUM = 1000;
  stmt->setMaxIterations(MAX_ITERATION_NUM+10);
  stmt->setAutoCommit(false);
  int iterCount = 0;
  for(unsigned i = 0; i < curChunkInfo->dstChunkInfo[slotNumber]->eventInfo.size(); i++) {
    stmt->setUInt(1,curChunkInfo->dstChunkInfo[slotNumber]->fileId);
    stmt->setUInt(2,curChunkInfo->dstChunkInfo[slotNumber]->eventInfo[i].eventNumber);
    stmt->setUInt(3,curChunkInfo->dstChunkInfo[slotNumber]->eventInfo[i].eventSize);
    stmt->setUInt(4,curChunkInfo->dstChunkInfo[slotNumber]->eventInfo[i].eventOffset);
    stmt->setUInt(5,curChunkInfo->dstChunkInfo[slotNumber]->eventInfo[i].triggerMask);
    for(int j = 0; j < 3; j++) {
      stmt->setDouble(6+j,curChunkInfo->dstChunkInfo[slotNumber]->eventInfo[i].extVar[j]);
    }
    stmt->addIteration();
    iterCount++;
    if(iterCount >= MAX_ITERATION_NUM) {
      if(!stmt->executeArrayUpdate(iterCount))
	return false;
      iterCount = 0;
    }
  }
  if(iterCount != 0) {
    if(!stmt->executeArrayUpdate(iterCount))
      return false;
  }
  if(!dstFile.Close()) return false;
  if(!session->commit()) return false;
  }
catch(SQLException& e)
  {
    cerr << e.getMessage() << endl;
    return false;
  }
  cout << curChunkInfo->dstChunkInfo[slotNumber]->eventInfo.size() << " DST events have been written in DB." << endl;
  return session->stopSession();
}

void CsOraRun::Abort()
{
  if(!isDstProduction && !isRunLogUpdate)
    return;
  session->rollback();
  dstFile.Discard();
}

void CsOraRun::Print(ostream& msgout) const
{
  msgout << "Run number:      " << runNumber << endl;
  msgout << "DST slot number: " << slotNumber << endl;
  msgout << "Started at:      " << CsTime(start,0) << endl;
  msgout << "Finished at:     " << CsTime(stop,0) << endl;
  msgout << "Run logs:        \n\"" << runLog << "\"" << endl;
}
