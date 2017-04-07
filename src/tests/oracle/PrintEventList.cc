#include <dirent.h>
#include <stdio.h>

#include "DaqDataDecoding/DaqEventsManager.h"

int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      printf("Usage: PrintEventList <raw file name> <event list file>\n");
      return 1;
    }
  std::string raw_file_name = argv[1];
  std::string event_list_name = argv[2];
  FILE* event_list_file = fopen(event_list_name.c_str(),"a");
  if(!event_list_file)
    {
      printf("Cannot open file: %s\n",event_list_name.c_str());
      return 2;
    }
  CS::DaqEventsManager* eventManager = new CS::DaqEventsManager();
  if(eventManager)
    {
      eventManager->AddDataSource(raw_file_name.c_str());
      while(eventManager->ReadEvent())
	{
	  const CS::DaqEvent& event = eventManager->GetEvent();
	  const CS::DaqEvent::Header& header = event.GetHeader();
	  unsigned run_number     = header.GetRunNumber();
	  unsigned spill_number   = header.GetBurstNumber();
	  unsigned event_in_spill = header.GetEventNumberInBurst();
	  fprintf(event_list_file,"%u %u %u\n",run_number, spill_number, event_in_spill);
	}
    }
  delete eventManager;
  fclose(event_list_file);
  return 0;
}
