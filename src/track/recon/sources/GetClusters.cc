#include <iostream>
#include "Coral.h"
#include "CsHelix.h"
#include "Recon.h"
#include "RecCall.h"
#include "CsTrafficPrepattern.h"
#include <list>
#include <algorithm>

using namespace std;

 void  RecCall::GetClusters(list<CsCluster>& listClus) 
{

  list<CsCluster*> UnusedClus; //unused clusters  after Traffic prepattern
  listClus.clear();   // all clusters from Traffic (segments without mom. and unused)
  UnusedClus.clear();
  	
   
  // clusters used by Traffic in segments without momentum before SM2
  CsEvent* event = CsEvent::Instance(); 
  list<CsTrack*> SegmTra=event->getTracks();   list<CsTrack*>::iterator  It;
  
  for( It=SegmTra.begin(); It!=SegmTra.end(); It++ ) {
    vector< CsHelix> lhelix=(*It)->getHelices();
    if ( fabs((lhelix[0]).getCop()) > 0.) continue;  // take only tracks  without P <==> 1/mom==0.0
    if ( (lhelix[1]).getZ() > Recon::ref().magType[1] || (lhelix[1]).getZ() < 0.0 ) continue; // and before SM2
    // printf("PED  %f  %f %f \n", (lhelix[0]).getCop(),fabs(1/(lhelix[0]).getCop()),(lhelix[1]).getZ());
    list<CsCluster*> clus=(*It)->getClusters();  list<CsCluster*>::iterator Ic;
    
    for( Ic=clus.begin(); Ic!=clus.end(); Ic++ ) {
      listClus.push_back(*(*Ic) ); 
    } // clusters loop
  } //   tracks loop
  
  // clusters unused by  Traffic before SM2      
  CsTrafficPrepattern   UnCl;
  UnusedClus=UnCl.CsTrafficPrepattern::getUnusedClusters(); list<CsCluster*>::const_iterator Icu;
  for( Icu=UnusedClus.begin(); Icu!=UnusedClus.end(); Icu++ ) {
    if ( (*Icu)->getW()  > Recon::ref().magType[1] || (*Icu)->getW() < 0.0 ) continue;
    listClus.push_back(*(*Icu));  }
 

  listClus.sort();	

 

//      printf( "GetClusters: All   %d \n",listClus.size());
   
};












