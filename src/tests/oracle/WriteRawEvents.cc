#include <iostream>
#include <cstdio>
#include <string>
#include <vector>

#include "CsOraSession.h"
#include "CsStoreMisc.h"
#include "RunInfo.h"
#include "FlatFile.h"

using namespace std;
using namespace oracle::occi;

union EventNumber
{
  unsigned long event_number;
  struct
  {
    unsigned event_in_spill : 20;
    unsigned spill_number   : 12;
  } spill_event;
};

class EventInfo
{
  friend class EventRawFile;
  static FlatFile* inputFile;
protected:
  EventNumber   eventNumber;
  unsigned      runNumber;
  string        directory;
  string        fileName;
  unsigned long fileOffset;
  unsigned long eventSize;
public:
  EventInfo(unsigned run_number, unsigned spill_number, unsigned event_in_spill)
  {
    runNumber = run_number;
    eventNumber.spill_event.spill_number   = spill_number;
    eventNumber.spill_event.event_in_spill = event_in_spill;
    fileOffset = eventSize = 0;
  }
  EventInfo(unsigned run_number, unsigned unique_event_number)
  {
    runNumber = run_number;
    eventNumber.event_number = unique_event_number;
    fileOffset = eventSize = 0;
  }
  ~EventInfo() {};
  bool Init();
  bool ReadEvent(vector<unsigned char>& buffer);
};

FlatFile* EventInfo::inputFile = NULL;

bool EventInfo::Init()
{
  CsOraSession* session = CsOraSession::instance();
  cout << "Run number: " << runNumber;
  cout << " Spill: " << eventNumber.spill_event.spill_number;
  cout << " Event in spill: " << eventNumber.spill_event.event_in_spill << endl;
  string sql_query = string("select RAWEV_FILE_ID,EVENT_SIZE,RAWEV_FILE_OFFSET from EVENT_HEADERS ")
    + "where RUN_NUMBER=" + itostr(runNumber) + " and EVENT_NUMBER=" + utostr(eventNumber.event_number);
  ResultSet* rs = session->executeQuery(sql_query);
  if(!rs || !rs->next())
    {
      cerr << "Cannot find event in database:\n";
      cerr << "Run number:     " << runNumber << "\n";
      cerr << "Spill number:   " << eventNumber.spill_event.spill_number << "\n";
      cerr << "Event in spill: " << eventNumber.spill_event.event_in_spill << endl;
      return false;
    }
  unsigned file_id = rs->getUInt(1);
  eventSize        = rs->getUInt(2);
  fileOffset       = rs->getUInt(3);

  sql_query = string("select FILE_DIR,FILE_NAME from FILE_MAPS where FILE_ID=") + utostr(file_id);
  rs = session->executeQuery(sql_query);
  if(!rs || !rs->next())
    {
      cerr << "DB error: incorrect file_id=" << file_id << " for the event:\n";
      cerr << "Run number:     " << runNumber << "\n";
      cerr << "Spill number:   " << eventNumber.spill_event.spill_number << "\n";
      cerr << "Event in spill: " << eventNumber.spill_event.event_in_spill << endl;
      return false;
    }
  directory = rs->getString(1);
  fileName  = rs->getString(2);
  if(!inputFile || inputFile->getDir() != directory || inputFile->getFileName() != fileName)
    {
      delete inputFile;
      inputFile = new FlatFile;
      cout << "Read event from file: " << directory << fileName << endl;
      if(!inputFile->Open(directory,fileName))
	{
	  cerr << "Cannot open castor file: " << directory << fileName << endl;
	  return false;
	}
    }
  return true;
}

bool EventInfo::ReadEvent(vector<unsigned char>& buffer)
{
  if(!inputFile || !inputFile->Seek(fileOffset))
    {
      cerr << "CASTOR error, cannot positioning file: " << directory << fileName << endl;
      return false;
    }
  buffer.resize(eventSize);
  if(!inputFile->Read((void*)(&buffer[0]),eventSize))
    {
      cerr << "CASTOR error, cannot read event from the file: " << directory << fileName << endl;
      return false;
    }
  return true;
}

class EventRawFile
{
  vector<EventInfo*> events;
  CsOraSession*      session;
  FlatFile*          outputFile;
public:
  EventRawFile(int run_number, const string& list_file_name, const string& castor_dir);
  ~EventRawFile();
};

EventRawFile::EventRawFile(int run_number, const string& list_file_name, const string& castor_dir)
{
  try
    {
      session = CsOraSession::instance();
      session->init(run_number,false);
    }
  catch(...)
    {
      throw 5;
    }
  FILE* list_file = fopen(list_file_name.c_str(),"r");
  if(!list_file)
    {
      cerr << "Cannot open event list file: " << list_file_name << endl;
      throw 3;
    }
  char line[1024];
  int nline = 0;
  while(fgets(line,1023,list_file))
    {
      nline++;
      unsigned rn = 0;
      unsigned sn = 0;
      unsigned en = 0;
      int retVal = sscanf(line,"%u %u %u", &rn, &sn, &en);
      if(retVal == 3 && rn != (unsigned)run_number) continue;
      if(retVal == 3)
	{
	  EventInfo* event = new EventInfo(rn,sn,en);
	  events.push_back(event);
	}
      else if(retVal == 2) // spill_number<space>event_in_spill format
	{
	  EventInfo* event = new EventInfo(run_number,rn,sn);
	  events.push_back(event);
	} 
      else if(retVal == 1) // (spill_number<<20)|event_in_spill format
	{
	  EventInfo* event = new EventInfo(run_number,rn);
	  events.push_back(event);
	}
      else
	{
	  cerr << "Error in line " << nline << " of " << list_file_name << endl;
	  throw 4;
	}
    }
  fclose(list_file);
  try
    {
      outputFile = new FlatFile;
      string output_file_name = string("Run_")+itostr(run_number)+".raw";
      if(!outputFile || !outputFile->Create(castor_dir,output_file_name))
	{
	  cerr << "Cannot create castor file: " << castor_dir << output_file_name << endl;
	  throw 8;
	}
      if(events.size() != 0)
	{
	  for(unsigned i = 0; i < events.size(); i++)
	    {
	      if(!events[i]->Init())
		{
		  throw 7;
		}
	      vector<unsigned char> buffer;
	      if(!events[i]->ReadEvent(buffer))
		{
		  cerr << "ERROR: Cannot read event:\n";
		  cerr << "Run number:     " << run_number << "\n";
		  cerr << "Spill number:   " << events[i]->eventNumber.spill_event.spill_number << "\n";
		  cerr << "Event in spill: " << events[i]->eventNumber.spill_event.event_in_spill << endl;
		  throw 9;
		}
	      if(buffer.size() == 0)
		{
		  cerr << "Zerro event size ERROR:\n";
		  cerr << "Run number:     " << run_number << "\n";
		  cerr << "Spill number:   " << events[i]->eventNumber.spill_event.spill_number << "\n";
		  cerr << "Event in spill: " << events[i]->eventNumber.spill_event.event_in_spill << endl;
		  continue;
		}
	      if(!outputFile->Append((void*)(&buffer[0]),buffer.size()))
		{
		  outputFile->Discard();
		  throw 10;
		}
	    }
	  delete EventInfo::inputFile;
	  EventInfo::inputFile = NULL;
	}
      else
	{
	  outputFile->Discard();
	}
      session->stopSession();
    }
  catch(int& nError)
    {
      throw nError;
    }
  catch(...)
    {
      throw 6;
    }
}

EventRawFile::~EventRawFile()
{
  for(unsigned i = 0; i < events.size(); i++)
    delete events[i];
  try
    {
      delete outputFile;
    }
  catch(...)
    {
      throw 11;
    }
}

int main(int argc, char* argv[])
{
  if(argc != 4)
    {
      cout << "Usage: WriteRawEvents <Run number> <Event list file> <castor directory>" << endl;
      return 1;
    }
  int run_number = 0;
  if(sscanf(argv[1],"%d", &run_number) != 1 || run_number == 0)
    {
      cerr << "Incorrect run number: " << argv[1] << endl;
      return 2;
    }
  string list_file_name = argv[2];
  string castor_dir     = argv[3];
  try
    {
      EventRawFile* raw_file = new EventRawFile(run_number, list_file_name, castor_dir);
      delete raw_file;
    }
  catch (int& error_number)
    {
      return error_number;
    }
  return 0;
}
