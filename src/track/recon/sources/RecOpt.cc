#include "Coral.h"
#define   INI_STATIC_HERE
#include "RecOpt.h"
#include <iostream>
#include "CsOpt.h"

using namespace std;

       bool     RecOpt::getRecOptions()
{
      
cout<<endl<<"RecOpt::getRecOptions ==> Recon job control options :"<<endl;

  CsOpt* opt = CsOpt::Instance();

  //  int    n;
  //vector<int> vi;
  //vector<int>::iterator ivi;

  string tag("ReC");
  string key;
  int size;


  //  if(!opt->getOpt( tag, "ReHits",vi)) goto err;
  //  cout << tag  << " " << key  <<" : " << endl;
  
  //   n = 0;
  //  for( ivi=vi.begin(); ivi!=vi.end(); ivi++ ) {
  //    n++;
  //   cout << *ivi << " "; 
  //    if( n%10 == 0 ) cout << endl;
  //  }

  // int number 
  if(!opt->getOpt( tag, "Switch", RecOpt::Switch)) goto err;

  //int array
  size=sizeof(RecOpt::McPar)/sizeof(int);
  if(!RecOpt::getOptArray(tag,"McPar",RecOpt::McPar   ,   size)) goto err;

 //int array
  size=sizeof(RecOpt::inPar)/sizeof(int);
  if(!RecOpt::getOptArray(tag,"inPar",RecOpt::inPar   ,   size)) goto err;
 
 //int array
   size=sizeof(RecOpt::nbPlanes)/sizeof(int);
   if(!RecOpt::getOptArray(tag,"nbPlanes",RecOpt::nbPlanes  ,   size)) goto err;
  

  //float array
    size=sizeof(RecOpt::activeRegion)/sizeof(float);
  if(!RecOpt::getOptArray(tag,"activeRegion",RecOpt::activeRegion  ,   size)) goto err;
  
  //float array
    size=sizeof(RecOpt::roadPar)/sizeof(float);
  if(!RecOpt::getOptArray(tag,"roadPar",RecOpt::roadPar  ,   size)) goto err;
  
    size=sizeof(RecOpt::ReProj)/sizeof(float);
  if(!RecOpt::getOptArray(tag,"ReProj",RecOpt::ReProj   ,   size)) goto err;
   
  size=sizeof(RecOpt::roadLengh)/sizeof(float);
  if(!RecOpt::getOptArray(tag,"roadLengh",RecOpt::roadLengh   ,   size)) goto err;
  
  size=sizeof(RecOpt::chiPar)/sizeof(float);
  if(!RecOpt::getOptArray(tag,"chiPar",RecOpt::chiPar   ,   size)) goto err;
  
  size=sizeof(RecOpt::ReTarget)/sizeof(float);
  if(!RecOpt::getOptArray(tag,"ReTarget",RecOpt::ReTarget   ,   size)) goto err;

  //int array
  size=sizeof(RecOpt::blacklist)/sizeof(int);
  if(!RecOpt::getOptArray(tag,"blacklist",RecOpt::blacklist  ,   size)) goto err;

  //int array
   size=sizeof(RecOpt::deadPlanes1)/sizeof(int);
   if(!RecOpt::getOptArray(tag,"deadPlanes1",RecOpt::deadPlanes1  ,   size)) goto err;
 
  //int array
   size=sizeof(RecOpt::deadPlanes2)/sizeof(int);
   if(!RecOpt::getOptArray(tag,"deadPlanes2",RecOpt::deadPlanes2  ,   size)) goto err;

  cout<<endl<<"-----------------------------------------"<<endl<<endl;

  return(true);



 err:
  cout<<endl<<endl
      <<"RecOpt::getRecOptions ==> "<<endl
      <<"One or more job control flags or parameters, needed for Recon package, "<<endl
      <<"are missing in your options file. "<<endl
      <<"Please check lines with 'ReC' tag"<<endl<<endl<<endl;
  exit(1);
  return(false);


};

// from S.G 

template <class Typ>

bool RecOpt::getOptArray(string tag, string key, Typ* a, int size)
  //
  // Store values of the found "key" to int array a[size]
  //
{
  CsOpt* opt = CsOpt::Instance();
  int k;
  vector<Typ> v;
  typename vector<Typ>::iterator iv;
  bool rc;

  if( !opt->getOpt( tag, key, v)){ // not found ?

    return(false);
  } else if (int(v.size()) >  size){   // wrong size ?
    cout<<"RecOpt::getOptArray ==> tag = "<<tag<<"  key = "<<key<<endl;
    cout<<"Maximum number of elements = "<< size <<". Was found "<<v.size()<<endl;
    cout<<endl;
    return(false);
  } else { // OK
    //print
    cout<<"--- "<<key<<" : ";
    for(k=0,iv=v.begin(); iv != v.end(); iv++,k++){
      a[k]=(*iv); //store to array
      if(k%10 == 0) cout<<endl;
      cout<<a[k]<<"\t ";
    }
    cout<<endl;
  }
  return(true);
}
















