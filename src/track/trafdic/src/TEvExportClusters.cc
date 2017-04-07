#include "CsCluster.h"
#include "TEv.h"
#include "TOpt.h"

std::list<CsCluster*> TEv::ExportClusters(std::string opt)
{
  
  if(opt != "all" && opt != "unused") {
    std::cout<<"TEv::ExportClusters ==> unknown argumet 'opt' value : "<<opt<<std::endl;
    assert(false);
  }

  std::list<CsCluster*> lCsCl;
  if(TOpt::ReMode[30] > 0) return(lCsCl); // export nothing

  for(int ih = 0; ih < int(vecHit.size()); ih++){ // loop over all hits
    if( opt == "unused" && !vecHit[ih].sTrackID().empty() ) continue; // hit is associated with some track
    CsCluster* pCsCl =  vecHit[ih].PtrClus(); 
    if(pCsCl == NULL) {
      std::cout<<"TEv::ExportClusters ==> Hit # "<<ih
	  <<" has NULL pointer to CsCluster"<<std::endl;
      assert(false);
    }
    lCsCl.push_back(pCsCl);
    
    // mirror cluster for drift detector
    CsCluster* pCsClMirr =  vecHit[ih].PtrClusMirr();
    if( pCsClMirr != NULL) lCsCl.push_back(pCsClMirr);

  } // end of loop over hits
  
  return(lCsCl);
}









