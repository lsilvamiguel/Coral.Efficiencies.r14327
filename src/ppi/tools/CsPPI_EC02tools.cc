#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
//#include "CsPPI_EC02time.h"
#include "data_struct.h"
#include "CsPPI_EC02tools.h"
#include <iostream>
#include<vector>
using namespace std;


int ReadFile(const char* fname,const bool first,int &global_run_nmb,unsigned int data_array[200][6][8][8][5],unsigned int global_checksum[200]) {
  int file =open(fname,O_RDONLY);
  if(file==-1) {
    cerr<<"cannot open file: "<<fname<<" aborting!\n";
    exit(-101);
  }


  unsigned int l_data_array[200][6][8][8][5];
  unsigned int l_checksum[200];
  unsigned int bin[5];
  ec02_tcorr_header_t header1;
  ec02_tcorr_header2_t checksum_header;
  
//   int time;
  int run;
  ssize_t bytes_read=0;
  bytes_read = pread(file, &header1, 20,0);
  tm * ptm=gmtime((time_t*)&(header1.time));
  time_t ttt=*(time_t*)&(header1.time);
//   time(&ttt);
  cout<<"Run Number:  "<<header1.run<<"   creation time: "<<   ctime((time_t*)(&ttt))<<endl;
//check if you are merging files from only one run
// ----------------------------------------
  if(first)
    global_run_nmb=header1.run;
  else {
    if(global_run_nmb!=header1.run) {
      cerr<<"You try to merge files from different runs\n run1: "<< global_run_nmb<<" run2: "<<header1.run<<endl;
      exit(-102);
    }
  }
// ----------------------------------------


  bytes_read += pread(file, &checksum_header,800,bytes_read);


//get data for 200 spills
// ----------------------------------------
  for(int spill=0;spill<200;spill++) {
    l_checksum[spill]=0;
    for(int srcid=0;srcid<6;srcid++) {
      for(int port=0;port<8;port++) {
        for(int ch_grp=0;ch_grp<8;ch_grp++) {
          bytes_read+= pread(file, l_data_array[spill][srcid][port][ch_grp], 20,bytes_read);
          for(int i=0; i<5; i++) {
            l_checksum[spill]+=l_data_array[spill][srcid][port][ch_grp][i];
          }
        }
      }
    }
  }
// ----------------------------------------



//check if checksum in the current file is correct
// ----------------------------------------
  for(int i=0;i<200;i++) {
    if(l_checksum[i]!=checksum_header.checksum[i]) {
      cerr << "checksum for file "<< fname<<" does not match, aborting!\n" ;
      exit(-100);
    }
  }
// ----------------------------------------
  close(file);


//add data from this chunk to all data
// ----------------------------------------
  for(int spill=0;spill<200;spill++)  {
    global_checksum[spill]+=l_checksum[spill];
    for(int srcid=0;srcid<6;srcid++) 
      for(int port=0;port<8;port++) 
        for(int ch_grp=0;ch_grp<8;ch_grp++) 
          for(int bin=0; bin<5; bin++) 
            data_array[spill][srcid][port][ch_grp][bin]+=l_data_array[spill][srcid][port][ch_grp][bin];
  }
// ----------------------------------------



}

int WriteFile(const char* fname,int &global_run_nmb,unsigned int data_array[200][6][8][8][5],unsigned int global_checksum[200]) {
  int result;
  ec02_tcorr_data_t* ec02_tcorr;
  int nentries(200*6*8*8);
  int filesize((nentries+41)*sizeof(ec02_tcorr_data_t));
  int output_file;

  output_file = open(fname, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0644);
  if (output_file == -1) {
    cerr<<"Error opening output file "<<  __FILE__<< __LINE__<<endl;
    exit(-200);
  }
  result = lseek(output_file, filesize-1, SEEK_SET);
  if (result == -1) {
    close(output_file);
    cerr<<"Stretching output file to requested size "<<  __FILE__<< __LINE__<<endl;
    exit(-201);
  }
  result = write(output_file, "", 1);
  if (result != 1) {
    close(output_file);
    cerr<<"Error writing to output fil "<<  __FILE__<< __LINE__<<endl;
    exit(-202);
  }

  time_t ttt;//=*(time_t*)&(header1.time);
  time(&ttt);
  ec02_tcorr_header_t header;
  header.version=1;
  header.time=ttt;
  header.run=global_run_nmb;
  header.nu_1=0;
  header.nu_2=0;
  /* mmapp.
   *      */
  ec02_tcorr = (ec02_tcorr_data_t*)mmap(0, filesize, PROT_READ | PROT_WRITE, MAP_SHARED, output_file, 0);
  if (ec02_tcorr == MAP_FAILED) {
    close(output_file);
    cerr<<" Error mmapping the output file "<<  __FILE__<< __LINE__<<endl;
           exit(-203);
  }
  ec02_tcorr[0]=*((ec02_tcorr_data_t*)&header);
  for(int spill=0;spill<200;spill++){
    ((unsigned int*)ec02_tcorr)[5+spill]=global_checksum[spill];
  }
  
  for(int spill=0;spill<200;spill++) 
    for(int srcid=0;srcid<6;srcid++) 
      for(int port=0;port<8;port++) 
        for(int ch_grp=0;ch_grp<8;ch_grp++) 
          ec02_tcorr[41+(spill*384+64*(srcid)+8*port+ch_grp)]=*(ec02_tcorr_data_t*)data_array[spill][srcid][port][ch_grp];

  if (munmap(ec02_tcorr, filesize) == -1) {
    cerr<< "Error un-mmapping the output file"<<  __FILE__<< __LINE__<<endl;
  }
  close(output_file);



  return 0;

}

void
init_d_array(unsigned int data_array[200][6][8][8][5],unsigned int global_checksum[200]) {
  for(int spill=0;spill<200;spill++) {
    global_checksum[spill]=0; 
    for(int srcid=0;srcid<6;srcid++) 
      for(int port=0;port<8;port++) 
        for(int ch_grp=0;ch_grp<8;ch_grp++) 
          for(int bin=0; bin<5; bin++) 
            data_array[spill][srcid][port][ch_grp][bin]=0;
  }
}
