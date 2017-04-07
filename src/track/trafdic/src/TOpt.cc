// $Id: TOpt.cc 14069 2015-09-17 20:44:46Z lsilva $

/*! 
  Takes TraFFiC related information from CORAL options file
  and stores it in TOpt static data members
*/

#include <vector>
#include <iterator>
#include <typeinfo>
#include <stdlib.h>
#include <cassert>
#include <cstring>

#include "CsOpt.h"
#define   INIT_STATIC_HERE
#include "TOpt.h"

using namespace std;

bool TOpt::getOptions()
{

  TOpt::Print[0] = -1; // default (printing ON)
  
  cout<<endl<<"TOpt::getOptions ==> TRAFFIC job control options :"<<endl<<endl;;

  CsOpt* opt = CsOpt::Instance();

  int size;
  string key;

  { // int arrays
    int* arr; int s = sizeof(int);
    key = "Print";  arr = TOpt::Print;  size=sizeof(TOpt::Print)/s;
    getOptArray(key,arr,size,'n');
    key = "Graph";  arr = TOpt::Graph;  size=sizeof(TOpt::Graph)/s;
    getOptArray(key,arr,size,'n');
    key = "DetOff"; arr = TOpt::DetOff; size=sizeof(TOpt::DetOff)/s; 
    getOptArray(key,arr,size,'n');
    key = "Hist";   arr = TOpt::Hist;   size=sizeof(TOpt::Hist)/s;
    getOptArray(key,arr,size,'n');
    key = "iCut";   arr = TOpt::iCut;   size=sizeof(TOpt::iCut)/s;
    getOptArray(key,arr,size);
    key = "ReMode"; arr = TOpt::ReMode; size=sizeof(TOpt::ReMode)/s;
    getOptArray(key,arr,size);
    key = "iPRpar"; arr = TOpt::iPRpar; size=sizeof(TOpt::iPRpar)/s;
    getOptArray(key,arr,size);
    key = "CAOpt"; arr = TOpt::CAOpt; size=sizeof(TOpt::CAOpt)/s;
    getOptArray(key,arr,size,'n');
  }

  { // double arrays
    double* arr; int s = sizeof(double);
    key = "DefView";  arr = TOpt::DefView;  size=sizeof(TOpt::DefView)/s;
    getOptArray(key,arr,size);
    key = "dCut";     arr = TOpt::dCut;     size=sizeof(TOpt::dCut)/s;
    getOptArray(key,arr,size);
    key = "dPRpar";   arr = TOpt::dPRpar;   size=sizeof(TOpt::dPRpar)/s;
    getOptArray(key,arr,size);
    key = "Target";   arr = TOpt::Target;   size=sizeof(TOpt::Target)/s;
    if (getOptArray(key,arr,size,'n')) // Mandatory or not depends upon COMGeant being >= v7.3...
      // ...but we want to keep the TOpt class as much independent as possible
      // from the rest of coral. Therefore we postpone checking the availability
      // of the Target option against COMGeant version until "TSetup::Init"...
      TOpt::targetOpt = true;  // ...in the mean time, let's keep track of it.
    key = "MuWall";   arr = TOpt::MuonWall; size=sizeof(TOpt::MuonWall)/s;
    getOptArray(key,arr,size);
    key = "Calo";     arr = TOpt::Calo;     size=sizeof(TOpt::Calo)/s;
    getOptArray(key,arr,size);
    key = "RIPipe";   arr = TOpt::RICHPipe; size=sizeof(TOpt::RICHPipe)/s;
    getOptArray(key,arr,size,'n');
    key = "SmoothDet"; if(opt->getOpt("TraF", key, TOpt::SmoothDet)){
      if(TOpt::Print[0] != 0) 
	cout<<"--- SmoothDet : "<<TOpt::SmoothDet<<endl;
    }
    else {
      // "SmoothDet" takes precedence over "SmoothPos"
      key = "SmoothPos";arr = TOpt::SmoothPos;size=sizeof(TOpt::SmoothPos)/s;
      getOptArray(key,arr,size,'n');
    }
    key = "CAOptD"; arr = TOpt::CAOptD; size=sizeof(TOpt::CAOptD)/s;
    getOptArray(key,arr,size,'n');
  }
  { // string lists
    string* arr; int s = sizeof(string);
    key = "DetNameOff"; arr =TOpt::DetNameOff; size=sizeof(TOpt::DetNameOff)/s;
    getOptStrings(key,arr,size,'n');
    key = "Det2Go2Fit"; arr =TOpt::Det2Go2Fit; size=sizeof(TOpt::Det2Go2Fit)/s;
    getOptStrings(key,arr,size,'n');
    key = "DZisActive"; arr =TOpt::DZisActive; size=sizeof(TOpt::DZisActive)/s;
    getOptStrings(key,arr,size,'n');
    key = "Det2Ignore"; arr =TOpt::Det2Ignore; size=sizeof(TOpt::Det2Ignore)/s;
    getOptStrings(key,arr,size,'n');
  }
  // NB: argument 'n' in above getOptArray() calls means "not mandatory option"

  // strings
  key = "PSdir"; 
  if(opt->getOpt("TraF", key, TOpt::PSdir)){
    if(TOpt::Print[0] != 0) 
      cout<<"--- PSdir : "<<TOpt::PSdir<<endl;
  }
  

  // other stuff
  key = "Dicofit";
  if (TOpt::ReMode[14]&0x2 && !opt->getOpt("TraF", key, TOpt::Dicofit))
    goto err;
  opt->getOpt("TraF","GEMSpacers",TOpt::GEMSpacers);
  TOpt::ELossStraggling = false;
  opt->getOpt( "TraF", "ELossStraggling", TOpt::ELossStraggling );
  DY_InAcceptance = 1;
  opt->getOpt( "TraF", "DY_InAcceptance", TOpt::DY_InAcceptance);
  DY_VD_Chi2 = 0.3;
  opt->getOpt( "TraF", "DY_VD_Chi2", TOpt::DY_VD_Chi2);
  DY_VD_time = 5;
  opt->getOpt( "TraF", "DY_VD_time", TOpt::DY_VD_time);

  if(TOpt::Print[0] != 0) 
    cout<<endl<<"------ End of TRAFFIC job control options  -------"<<endl<<endl;

  // ********** SOME ACTIONS WITH FLAGS **********

  if((getenv("LSB_QUEUE") ||                     // ***** LXBATCH...
      getenv("PBS_ENVIRONMENT") ||               // ***** BATCH @ gridKa...
      getenv("ENVIRONMENT") &&                   // ***** BATCH @ Lyon...
      // "ENVIRONMENT" env. var. can be either of "BATCH", in the BQS system, or
      // "SEQUENTIAL_BATCH | PARALLEL_BATCH | INTERACTIVE_BATCH" in the GE one.
      strstr(getenv("ENVIRONMENT"),"BATCH") &&
      strcmp(getenv("ENVIRONMENT"),"INTERACTIVE_BATCH")) &&
      Graph[0]) {      // ***** ...SWITCH OFF GRAPHICS IN BATCH MODE *****
    if (TOpt::Print[0]) {
      if(getenv("LSB_QUEUE"))
	cout<<"--- env. \"LSB_QUEUE\" = "<<getenv("LSB_QUEUE")<<endl;
      else if (getenv("PBS_ENVIRONMENT"))
	cout<<"--- env. \"PBS_ENVIRONMENT\" = "<<getenv("PBS_ENVIRONMENT")<<endl;
      else if (getenv("ENVIRONMENT") &&
	       strstr(getenv("ENVIRONMENT"),"BATCH"))
	cout<<"--- env. \"ENVIRONMENT\" = "<<getenv("ENVIRONMENT")<<endl;
      cout<<"--- Seems it's a batch job ==> Event Display is set to 'OFF' !"
	<<endl<<endl;
    }
    Graph[0] = 0;
  }


  if(TOpt::Hist[0] == 0) { // all histograms are OFF
    size=sizeof(TOpt::Hist)/sizeof(int);
    for(int i = 1; i < size; i++) TOpt::Hist[i] = 0;
  }

  if(TOpt::Print[0] == 0) { // all prints are OFF
    size=sizeof(TOpt::Print)/sizeof(int);
    for(int i = 1; i < size; i++) TOpt::Print[i] = 0;
  }

  if(TOpt::Graph[6] > 0) {
    TOpt::ReMode[1] = 2; 
    if(TOpt::Print[0] != 0) cout<<"TOpt::getOptions ==> Warning: ReMode[1] is set to 2"<<endl;
  }
  return(true);

 err:
  cout<<endl<<endl
      <<"TOpt::getOptions ==> "<<endl
      <<"Option <"<<key<<">, needed for TRAFFIC package, "<<endl
      <<"is missing in your options file. "<<endl
      <<"Please check lines with 'TraF' tag"<<endl<<endl<<endl;
  assert(false); exit(1); return(false);

}


//-----------------------------------------------------------------

template <class T>
bool TOpt::getOptArray(string key, T* a, int size, char mandatory)
  //
  // Store values of the found "key" to array a[size] of type T
  //
{
  CsOpt* opt = CsOpt::Instance();
  string tag("TraF");
  int k;
  
  // START: bg 2006/03/15
  typedef vector<T>   myv;
  myv v;
  typename myv::iterator iv;
  //vector<T>            v;
  //vector<T>::iterator iv;
  // END: bg 2006/03/15  

  if( !opt->getOpt( tag, key, v)){ // not found ?
    if(mandatory == 'm') {
      cout<<endl<<endl
	  <<"TOpt::getOptArray ==> "<<endl
	  <<"Option's array <"<<key<<">, needed for TRAFFIC package, "<<endl
	  <<"is missing in your options file or not correctly specified. "<<endl
	  <<"Please check lines with 'TraF' tag"<<endl<<endl<<endl;
      assert(false); exit(1);
    } else { //just print warning
      if(TOpt::Print[0] != 0) 
	cout<<"TOpt::getOptArray ==> Warning: error in option's array <"
	    <<key<<"> of TRAFFIC package."<<endl;
    }
    return(false);
  } else if (int(v.size()) >  size){   // wrong size ?
    cout<<"TOpt::getOptArray ==> tag = "<<tag<<"  key = "<<key<<endl;
    cout<<"Maximum number of elements = "<< size <<". Was found "<<v.size()<<endl;
    cout<<endl;
    return(false);
  } else { // OK
    if(TOpt::Print[0] != 0) cout<<"--- "<<key<<" : ";
    for(k=0,iv=v.begin(); iv != v.end(); iv++,k++){
      a[k]=(*iv); //strore to array
      if(TOpt::Print[0] != 0) {
	if(k%10 == 0) cout<<endl;
	cout<<a[k]<<"\t ";
      }
    }
    cout<<endl;
 
  }
  
  return(true);
}

//-----------------------------------------------------------------

template <class T>
bool TOpt::getOptStrings(string key, T* a, int size, char mandatory)
  //
  // Store values of the found "key" to array a[size] of type T
  //
{
  CsOpt* opt = CsOpt::Instance();
  string tag("TraF");
  int k;
  
  // START: bg 2006/03/15
  typedef list<T>   myl;
  myl v;
  typename myl::iterator iv;
  //list<T> v;
  //list<T>::iterator iv;
  // END: bg 2006/03/15  
  
  if( !opt->getOpt( tag, key, v)){ // not found ?
    if(mandatory == 'm') {
      cout<<endl<<endl
	  <<"TOpt::getOptStrings ==> "<<endl
	  <<"Option's array <"<<key<<">, needed for TRAFFIC package, "<<endl
	  <<"is missing in your options file. "<<endl
	  <<"Please check lines with 'TraF' tag"<<endl<<endl<<endl;
      assert(false);
    }
    return(false);
  } else if (int(v.size()) >  size){   // wrong size ?
    cout<<"TOpt::getOptStringArray ==> tag = "<<tag<<"  key = "<<key<<endl;
    cout<<"Maximum number of elements = "<< size <<". Was found "<<v.size()<<endl;
    cout<<endl;
    return(false);
  } else { // OK
    if(TOpt::Print[0] != 0) cout<<"--- "<<key<<" : ";
    for(k=0,iv=v.begin(); iv != v.end(); iv++,k++){
      a[k]=(*iv); //store to array
      if(TOpt::Print[0] != 0) {
	if(k%10 == 0) cout<<endl;
	cout<<a[k]<<"\t ";
      }
    }
    cout<<endl;
  }
  
  return(true);
}



