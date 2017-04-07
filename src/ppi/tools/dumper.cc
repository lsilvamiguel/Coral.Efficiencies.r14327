#include <stdio.h>
#include <time.h>
#include <getopt.h>
#include <stdlib.h>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include "CsPPI_EC02tools.h"
//#include "CsPPI_EC02time.h"
#include "data_struct.h"
#include "merger.h"
#include <iostream>
#include <vector>
using namespace std;




int main(int argc, char* argv[]) { 
  bool first=true;
  int global_run_nmb;
  unsigned int global_checksum[200];
  unsigned int data_array[200][6][8][8][5];

  float qcut=0;
  unsigned int stcut=0;
  bool cut=false;


  init_d_array(data_array,global_checksum);

  vector<string> InputFiles;
  if( argc <2) {
    cerr<<"to few arguments \n";
    PrintUsage(argv[0]);
    exit(-105);
  }
  char c;
  while( (c=getopt(argc,argv, "s:q:o:h")) != -1) {
    switch(c) {
      case 'h':
        PrintUsage(argv[0]);
        exit(-108);
        break;
      case 'q':
        qcut = atof(optarg);
        cut=true;
        break;
      case 's':
        stcut = atoi(optarg);
        cut=true;
        break;
      default:
        break;

    }
  }
  int nargs=argc-optind;
  for( int i=0; i<nargs;i++) {
    string tfname(argv[optind+i]);
    bool stored= false;
    for(unsigned int k=0; k < InputFiles.size(); k++){ 
      if(tfname== InputFiles.at(k)) {
        stored=true;
        break;
      }
    }
    if(!stored)
      InputFiles.push_back(tfname);
  }

  if(InputFiles.size()==0) {
    cerr<<"No Input files specified\n";
    exit(-108);
  }
  cout<<qcut<<endl;
  //loop over input files
  //---------------------------------------
  for(int i=0; i<InputFiles.size();i++)
  {
    ReadFile(InputFiles.at(i).c_str(),first,global_run_nmb,data_array,global_checksum);
    if(first) {
      first=false;
    }
    int shift;
    unsigned int stat;
    
    for(int spill=0;spill<200;spill++){
      for(int srcid=0;srcid<6;srcid++) {
        for(int port=0;port<8;port++) {
          for(int ch_grp=0;ch_grp<8;ch_grp++) {
            stat=0;
            shift=0;
            float qual;
            for(int bin=0;bin<5;bin++) {
              stat+= data_array[spill][srcid][port][ch_grp][bin];
              if(data_array[spill][srcid][port][ch_grp][bin]>data_array[spill][srcid][port][ch_grp][shift])
                shift=bin;
            }
            if(stat!=0)
              qual=(float)data_array[spill][srcid][port][ch_grp][shift]/stat;
            else
              qual=0;
            shift-=2;

            if( (qcut>=qual  || stcut>=stat)|| !cut){
            printf("spill: %3u  srcid: %u  port: %u chgrp: %u ||  ",spill+1,srcid+616,port,ch_grp);
            printf("%7u | %7u | %7u | %7u | %7u ||\t",data_array[spill][srcid][port][ch_grp][0],
                data_array[spill][srcid][port][ch_grp][1],data_array[spill][srcid][port][ch_grp][2],
                data_array[spill][srcid][port][ch_grp][3],data_array[spill][srcid][port][ch_grp][4]);
            printf("shift: %2i    statistics %7u   quality %1.3f\n",shift,stat,qual);
            }
          }
        }
      }
    }
 
  }
  //---------------------------------------
}

void PrintUsage(char* exe) {
  cout<<endl;
  cout<<"Usage:"<<endl;
  cout<<exe<<" [ -option(s)] [in_files] "<<endl;
  cout<<"Arguments:"<<endl;
  cout<<"  in_files      input binary files"<<endl;
  cout<<"Options:"<<endl;
  cout<<endl;

}
