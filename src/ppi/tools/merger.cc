#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "data_struct.h"
// #include "CsPPI_EC02time.h"
#include "CsPPI_EC02tools.h"
#include "merger.h"
#include <iostream>
#include<vector>
using namespace std;



int main(int argc, char* argv[]) { 
  bool first=true;
  int global_run_nmb;
  unsigned int global_checksum[200];
  unsigned int data_array[200][6][8][8][5];



  init_d_array(data_array,global_checksum);

  char* ofname=NULL;
  char* odir=NULL;
  vector<string> InputFiles;
  if( argc <2) {
    cout<<"to few arguments \n";
    PrintUsage(argv[0]);
    exit(-105);
  }
  char c;
  while( (c=getopt(argc,argv, "o:h")) != -1) {
    switch(c) {
      case 'h':
        PrintUsage(argv[0]);
        exit(-108);
        break;
      case 'o':
        odir=optarg;
        break;
      default:
        PrintUsage(argv[0]);
        exit(-108);
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
  if(odir==NULL) {
    cerr<<"No Output directory specified\n";
    exit(-108);
  }


  //loop over input files
  //---------------------------------------
  for(int i=0; i<(int)InputFiles.size();i++)
  {
    cout<<InputFiles.at(i)<<endl;
    ReadFile(InputFiles.at(i).c_str(),first,global_run_nmb,data_array,global_checksum);
    if(first) {
      first=false;
    }
  }
  //---------------------------------------
  asprintf(&ofname,"%s/EC02_timecorr_run%u_merged",odir,global_run_nmb);

  //check if file exists
  //----------------------------------------
  char* pre_fname;
  int fname_counter=0;
  asprintf(&pre_fname,"%s",ofname);
  {
    bool flag = false;
    do {
      fstream fin;
      fin.open(pre_fname,ios::in);
      if( fin.is_open() )
      {
        fname_counter++;
        flag=true;
        free(pre_fname);
        asprintf(&pre_fname,"%s_%u",ofname,fname_counter);
      }
      else {
        flag=false;
      }
      fin.close();
    } while(flag==true);
  }
  if(fname_counter!=0)
    ofname=pre_fname;

  //----------------------------------------
  //

  WriteFile(ofname,global_run_nmb,data_array,global_checksum);
  free(ofname);
}









void PrintUsage(char* exe) {
  cout<<endl;
  cout<<"Usage:"<<endl;
  cout<<exe<<" [ -option(s)] [in_files] "<<endl;
  cout<<"Arguments:"<<endl;
  cout<<"  in_files      input binary files"<<endl;
  cout<<"Options:"<<endl;
  cout<<"-o output file"<<endl;
  cout<<endl;



}
