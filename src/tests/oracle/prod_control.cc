#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

void fflush(FILE*& logFile, char* fileName)
{
  fclose(logFile);
  logFile = fopen(fileName,"a");
}

int main(int argc, char* argv[])
{
  int retCode = 0;
  char slot_number[32];
  char run_number[32];
  char chunk_number[32];
  char* env = getenv("DST_SLOT_NUMBER");
  if(!env) 
    {
      printf("$DST_SLOT_NUMBER is not defined.\n");
      return 1;
    }
  strcpy(slot_number,env);
  env = getenv("DST_RUN_NUMBER");
  if(!env) 
    {
      printf("$DST_RUN_NUMBER is not defined.\n");
      return 2;
    }
  strcpy(run_number,env);
  env = getenv("DST_CHUNK_NUMBER");
  if(!env) 
    {
      printf("$DST_CHUNK_NUMBER is not defined.\n");
      return 3;
    }
  strcpy(chunk_number,env);
  if(argc != 2) 
    {
      printf("Usage: prodcontrol <log directory name>\n");
      return 4;
    }
  char log_file_dir[512];
  sprintf(log_file_dir,"%s/%s",argv[1],run_number);
  DIR* directory = opendir(log_file_dir);
  if(!directory)
    {
      mode_t mode = S_IWGRP|S_IRGRP|S_IRUSR|S_IWUSR|S_IXUSR|S_IXGRP|S_IROTH|S_IXOTH;
      int ret_val = mkdir(log_file_dir,mode);
      if(ret_val != 0) 
	{
	  printf("Cannot create log directory: %s\n",log_file_dir);
	  return 6;
	}
    }
  closedir(directory);
  char log_file_name[512];
  sprintf(log_file_name,"%s/%s-%s-%s.log",log_file_dir,run_number,chunk_number,slot_number);
  FILE* logFile = fopen(log_file_name,"w");
  if(!logFile)
    {
      printf("Cannot create log file: %s\n", log_file_name);
      return 5;
    }
  time_t job_start_time      = time(NULL);
  time_t event_start_time    = 0;
  time_t first_event_start   = 0;
  unsigned last_event_number = 0;
  bool is_summary            = false;
  char input[1024];
  while(fgets(input,1023,stdin))
    {
      printf(input); fflush(stdout);
      if(strstr(input,"CORAL RETURN CODE:"))
	{
	  sscanf(input,"CORAL RETURN CODE: %d",&retCode);
	  time_t curTime = time(NULL);
	  fprintf(logFile, "CORAL END:     %u\t\tERROR CODE:  %d\n\n", curTime,retCode);
	  fprintf(logFile, "RUN TIME:      %u sec.\n",curTime-job_start_time);
	  fflush(logFile,log_file_name);
	  
	  is_summary = false;
	}
      else if(is_summary)
	{
	  fprintf(logFile,input);
	  fflush(logFile,log_file_name);
	}
      else if(strstr(input,"Total initialization time"))
	{
	  time_t curTime = time(NULL);
	  fprintf(logFile, "CORAL START:   %u\n\n",curTime);
	  fflush(logFile,log_file_name);
	  job_start_time = curTime;
	}
      else if(strstr(input,"DATABASE INIT started"))
	{
	  time_t curTime = time(NULL);
	  fprintf(logFile, "DB_INIT START: %u\n\n",curTime);
	  fflush(logFile,log_file_name);
	}
      else if(strstr(input,"DATABASE INIT ended"))
	{
	  time_t curTime = time(NULL);
	  fprintf(logFile, "DB_INIT END:   %u\n\n",curTime);
	  fflush(logFile,log_file_name);
	}
      else if(strstr(input,"DATABASE SCAN started"))
	{
	  time_t curTime = time(NULL);
	  fprintf(logFile, "DB_SCAN START: %u\n\n",curTime);
	  fflush(logFile,log_file_name);
	}
      else if(strstr(input,"DATABASE SCAN ended"))
	{
	  time_t curTime = time(NULL);
	  fprintf(logFile, "DB_SCAN END:   %u\n\n",curTime);
	  fflush(logFile,log_file_name);
	}
      else if(strstr(input,"1ST EVENT UPLOAD started") || strstr(input,"1ST EVENT DOWNLOAD started"))
	{
	  time_t curTime = time(NULL);
	  fprintf(logFile, "FIRST_EVENT:   %u\n",curTime);
	  fflush(logFile,log_file_name);
	  event_start_time  = curTime;
	  first_event_start = curTime;
	  last_event_number = 1;
	}
      else if(strstr(input,"--  Run #") && strstr(input,"  Event #"))
	{
	  int i1,i2,i3,i4,i5,i6,i7;
	  unsigned u1;
	  unsigned event_number = 0;
	  if(sscanf(input,"[%d:%d:%d] %d.%d.%d  --  Run # %d  Event # %u  (%d)",
		    &i1,&i2,&i3,&i4,&i5,&i6,&i7,&u1,&event_number) == 9)
	    {
	      time_t curTime = time(NULL);
	      float time_per_event = float(curTime-event_start_time)/(event_number-last_event_number);
	      float total_tpe = float(curTime-first_event_start)/event_number;
	      event_start_time  = curTime;
	      last_event_number = event_number;
	      fprintf(logFile, "EVENT %7d: %u\t\t%.3f\t\t%.3f\n",event_number,curTime,time_per_event,total_tpe);
	      fflush(logFile,log_file_name);
	    }
	}
      else if(strstr(input,"Total time in CORAL"))
	{
	  is_summary = true;
	}
    }
  fclose(logFile);
  return retCode;
}
