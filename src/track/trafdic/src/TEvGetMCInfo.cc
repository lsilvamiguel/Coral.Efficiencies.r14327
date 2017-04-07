/*!
  This function get MC information, important for 
  debugging and tests of tracking algorithms
  from the CORAL and store it in  corresponding containers of the
  TEv class object

  \warning  Do not put any tracks/hits selection in this function. 
  It can destroy MCTrack - MCHits - MCVetrex cross-references 

 */

#include "CsEvent.h"
#include "TEv.h"
#include "TOpt.h"
#include "partab.h"

using namespace std;

void TEv::GetMCInfo()
{

  CsEvent*  ev = CsEvent::Instance();
  const TSetup& setup = TSetup::Ref();
  //
  // fill vecKine
  // 
  list<CsMCTrack*> lt = ev->getMCTracks();
  list<CsMCTrack*>::iterator it;
  for(it=lt.begin(); it != lt.end(); it++){ // loop over CsMCTracks

    // cross-check MC track -> MC vertices -> incoming MC track 
    list<CsMCVertex*> lv = (*it)->getOutVertices();
    list<CsMCVertex*>::iterator iv;
    for(iv=lv.begin(); iv != lv.end(); iv++){    // loop over daughter CsMCVertices
      const CsMCTrack* intr = (*iv)->getInTrack();     // pointer to mother track
      if( intr != (*it) ) { // the same track?
	cout<<"TEv::GetMCInfo ==> Track -> daughter vertices -> Track mismatch: "<<endl;
	cout<<intr<<"  "<<(*it)<<endl;
	assert(false);
      }
    }
    
    TKine  k;
    k.ptrTrk  = (*it);
    k.e       = (*it)->getE();
    k.p[0]    = (*it)->getPZ();  
    k.p[1]    = (*it)->getPX();  
    k.p[2]    = (*it)->getPY();  

    // k.iVtxOrig    - will be filled in vertex block
    // k.vecVtxMCRef - will be filled in vertex block   
    // k.vecHitRef;  - will be filled later TEv::ImportClusters()

    // k.setTrackIDs - will be filled after track finding in TTrack::FindKine()

    
    // some checks

    if( k.ptrTrk->getInVertex() == NULL){
      cout<<"TEv::GetMCInfo ==> CsMCTrack: NULL pointer to CsMCVertex object "<<endl;
      assert(false);
    }
    if (k.ptrTrk->getParticle() == NULL) {
      cout<<"TEv::GetMCInfo ==> CsMCTrack: NULL pointer to CsMCParticle object "<<endl;
    }
    if (k.ptrTrk->getParticle()->getNumber() == 0) {
      int ignum = k.ptrTrk->getParticle()->getGeantNumber();
      cout<<"TEv::GetMCInfo ==>  particle with Geant number "<<ignum<<" not in PDG table. Mass (sqrt(E^2 -P^2)) = "<<k.MCmass()<<endl;
    }
    int ipart = (*it)->getParticle()->getGeantNumber(); // GEANT particle type
    if( Id2q[ipart] != (*it)->getParticle()->getCharge()){
      /*
	cout<<"TEv::GetMCInfo ==> PDG charge = "
	<< (*it)->getParticle()->getCharge()
	<<" PDG   name: "
	<<(*it)->getParticle()->getName()
	<<" ( PDG   Number = "
	<<(*it)->getParticle()->getNumber()<<" )"<<endl
	<<"             BUT Geant charge = "
	<<Id2q[ipart]
	<<" Geant name: "
	<<Id2nam[ipart]
	<<" ( Geant Number = "<<ipart<<" )"<<endl;
      */
    }
    if( fabs(k.MCmass()-k.Mass()) > 1.E-4 ){
      cout<<"TEv::GetMCInfo ==> PDG mass - GEANT mass > 0.1 Mev : "<<endl;
      cout<<k.Name()<<"\t E = "<<k.E<<"\t MCmass = "<<k.MCmass()<<"\t PDG mass = "<<k.Mass()<<endl;
    }


    //
    // fill vecHitMC
    // 
    list<CsMCHit*> lh = (*it)->getMCHits();
    list<CsMCHit*>::iterator ih;
    for(ih=lh.begin(); ih != lh.end(); ih++){  // loop over CsMCHit* of the MC track

      if( (*it) != (*ih)->getMCTrack() ) { 
	cout<<"TEv::GetMCInfo ==> MCTrack -> MCHit -> MCTrack mismatch "<<endl;
	assert(false);
      }

      THitMC h;
      h.ptrHit = (*ih);
      h.iOrig  = (*ih)->getOrigin();
      h.deltaT = (*ih)->getDTime();
      CsDetector* csdet = dynamic_cast<CsDetector*>((*ih)->getDet());
      if( csdet != NULL ) {
	h.iDet   = csdet->GetID();
	// check if detector to be ignored
	if(setup.sIgnoreIDs().find(h.iDet) !=  setup.sIgnoreIDs().end()) {
	  continue;
	}
      } else { // can't cast to CsDetector*. Not a tracking detector's hit (e.g. RICH) Skip it.
	continue;
      }
      h.xyz[0] = (*ih)->getZ()/10.;
      h.xyz[1] = (*ih)->getX()/10.;
      h.xyz[2] = (*ih)->getY()/10.;
      // h.setHitRef will be filled in TEv::ImportClusters

      TPlane& pl = const_cast<TPlane&>(setup.Id2Plane(h.iDet));

#ifdef DW_FIX_DZ_BUG
      if (h.IDet>=0) {
	const TDetect &d = setup.vDetect(pl.IDetRef);
	if (d.IType==16 /* i.e. DW */ && !d.InActive(h.xyz[1],h.xyz[2]))
	  continue;
      }
#endif

      int ih = vecHitMC.size();
      h.iKine = vecKine.size();                 // store THitMC -> TKine  reference  
      k.vecHitMCRef.push_back(ih);              // store TKine  -> THitMC references
      vecHitMC.push_back(h); // store MC hit 
      mapCsMCHit2HitMC[h.ptrHit]=ih;            // fill CsMCHit* -> THitMC index map


      pl.addHitMCRef(ih);                       // store TPlane -> THitMC objects references
   
   } // end of loop over CsMCHit*
    
    vecKine.push_back(k); // store it
  } // end of loop over CsMCTracks
  
  //
  // fill vecVtxMC
  //
  int nprim=0;
  list<CsMCVertex*> lv = ev->getMCVertices(); 
  list<CsMCVertex*>::iterator iv;
  for(iv=lv.begin(); iv != lv.end(); iv++){ // loop over CsMCVertices

    // cross-checks MC vertex -> MC tracks -> MC vertex

    const CsMCTrack* intr = (*iv)->getInTrack(); // pointer to mother track
    if( intr != NULL ) { // vtx _with_ incoming track: NOT a primary vtx , NOT a pileup track mother vtx
      // check OutVertices()
      list<CsMCVertex*> lov = intr->getOutVertices();
      list<CsMCVertex*>::iterator iov;
      if(lov.empty()){
	cout<<"TEv::GetMCInfo ==> Vertex -> mother Track -> OutVertices gives empty list"<<endl;
	assert(false);
      }
      int nsame=0;
      for(iov = lov.begin(); iov != lov.end(); iov++){
	if( (*iv) == (*iov) ) nsame++;
      }
      if(nsame != 1){
	cout<<"TEv::GetMCInfo ==> Vertex -> mother Track -> OutVertices -> Vertex mismatch. "
	    <<nsame<<endl;
	assert(false);
      }
      // check  OutTracks()
      list<CsMCTrack*> ltt = intr->getOutTracks();     // list of outgoing tracks
      list<CsMCTrack*>::iterator itt;
      for(itt = ltt.begin(); itt != ltt.end(); itt++){ // loop over outgoing tracks 
	if( (*itt) == NULL ){
	  cout<<"TEv::GetMCInfo ==> NULL pointer in the list, returned by CsMCTrack::getOutTracks()"<<endl;
	  assert(false);
	}
      }
    } else { // vertex _without_ incoming track : primary vtx or pileup track mother vtx

      
      if((*iv)->getOutTracks().size() == 0){
	cout<<"TEv::GetMCInfo ==> primary or 'pileup mother' vertex (#"<<(*iv)->getGnum()
	    <<", X = "<<(*iv)->getZ()/10.<<") has NO outgoing tracks."<<endl;
      }
      if((*iv)->getOutTracks().size() > 1) { // > 1 tracks in the vertex. Looks like primary. 
	nprim++;
      }
    } 

    TVtxMC v;
    v.ptrVtx = (*iv);
    v.tof    = (*iv)->getT();
    v.v[0]   = (*iv)->getZ()/10.;
    v.v[1]   = (*iv)->getX()/10.;
    v.v[2]   = (*iv)->getY()/10.;

    v.iOrig = -1;    int nn=0;
    for(int i = 0; i < int(vecKine.size()); i++){   // loop over vecKine
      if( (*iv) == vecKine[i].ptrTrk->getInVertex() ) { 
	v.vecTrkRef.push_back(i);                   // store TVtxMC -> daughter tracks references
	vecKine[i].iVtxOrig = vecVtxMC.size();      // store TKine track -> TVtxMC-of-origin reference 
      }
      if( (*iv)->getInTrack() == vecKine[i].ptrTrk ){
	v.iOrig = i;                                // store TVtxMC -> mother TKine track reference
	nn++;
      }
    } //end of loop over vecKine

    if( nn > 1 ){
      cout<<"TEv::GetMCInfo ==> Number of found mother tracks for the vertex is "<<nn<<endl;
      assert(false);
    }

    vecVtxMC.push_back(v);
  } // end of loop over CsMCVertices

  if( nprim > 1 ){
    cout<<"TEv::GetMCInfo ==> "<<nprim<<" vertices with NO incomming track and >= 2 outgoing tracks"<<endl
	<<"It looks like more then 1 primary vertex in the event"<<endl;
  }



  // Yet another pass through vecKine and vecVtxMC :-(
  for(int itk = 0; itk < int(vecKine.size()); itk++){ // loop over vecKine 
    list<CsMCVertex*> lv = vecKine[itk].ptrTrk->getOutVertices();
    list<CsMCVertex*>::iterator iv;
    for(iv=lv.begin(); iv != lv.end(); iv++){    // loop over daughter CsMCVertices
      for(int ivt = 0; ivt < int(vecVtxMC.size()); ivt++){ // find index of current CsMCVertex in vecVtxMC
	if( vecVtxMC[ivt].ptrVtx == (*iv) ) {
	  vecKine[itk].vecVtxMCRef.push_back(ivt);         // store TKine track  -> daughters TVtxMC references
	  break; 
	}
      }
    }// end of loop over daughter CsMCVertices
  } // end of loop over vecKine

  // Printouts

  if(TOpt::Print[1] > 0) this->PrintMC(0);

}















