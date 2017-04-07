// $Id: TEv.cc 13556 2012-08-19 22:56:21Z ybedfer $

#include "CsInit.h"
#include "CsEvent.h"
#include "TEv.h"
#include "TSetup.h"
#include "TOpt.h"
#include "TDisplay.h"

using namespace std;

/*!
  \file
  Constructor(s), destructor and trivial methods
  of the class TEv
*/

/*
  (Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TEv.cc":
  i) Do not read MC data upon "ReMode[1]==3".)
  ii) Init "eventTime", "eventTRef".
  iii) "TrigMaskString" lists master trigger first.
*/

TEv* TEv::address = 0;                   // init static pointer

unsigned int  TTrack::TrackCounter = 0;  // init reconstructed track counter

//! Constructor
TEv::TEv():
  isMC(false)
{
  if (address==0) {                    // Protection against multiple instances
    address = this;
    TTrack::TrackCounter = 0;          // Reset reconstructed track counter
    if (CsEvent::Instance()->isAMonteCarloEvent() && // ***** If MC... *****
	TOpt::ReMode[1]!=3) {// ...and depending upon TRaFFiC's secret trap
      isMC = true; this->GetMCInfo();  // Read in MC data
      BMSSmearing = 0; // Reset BMS smearing, cf. "TEv::TracksFit2"
    }
    hadronJob = CsInit::Instance()->IsAHadronJob();
    if (TOpt::Graph[6]>0)
      TDisplay::Ref().Draw(1);         // Draw setup
  }

  // Get run#, event #, trigger mask from Coral
  run        = unsigned(CsEvent::Instance()->getRunNumber());
  event      = unsigned(CsEvent::Instance()->getEventNumberInRun());
  ev_in_burst= unsigned(CsEvent::Instance()->getEventNumberInBurst());
  trig_mask  = unsigned(CsEvent::Instance()->getTriggerMask());

  eventTime = 0; eventTRef = 0; // Default event time and time reference = 0
  reTracking = false;           // Default status: not currently in re-tracking.
}

//! Destructor
TEv::~TEv()
{
  // reset TPlane -> THitMC, THit references in the setup.vecPlane
  TSetup::Ptr()->Reset();

  address = 0;
}

//! Accessor to the object (by pointer)
TEv* TEv::Ptr()
{
  // No check on existence of the object. Normaly, TEv::Ref() has to be used.
  return(address);
}

//! Accessor to the object (by reference)
TEv& TEv::Ref()
{
  if (address != 0) {
    return(*address);
  }
  cout<<"TEv::Ref() ==> The object of TEv class is not yet created or already destructed"<<endl;
  assert(false);
  return(*address); // just to get rid of compiler warnings
}

CsEvent* TEv::ptrEvt()
{
  return(CsEvent::Instance());
}


unsigned int TEv::TrigMask() const
{
  return trig_mask;
}

#include "DaqDataDecoding/DaqEvent.h"
#include "DaqDataDecoding/TriggerTime.h"

string TEv::TrigMaskString() const
{
  unsigned int tr = TrigMask();
  string s; char c[3];

  if (!IsMC()) {
    //                                                  ***** GET MASTER TRIGGER
    // (I.e. trigger used to set master time by the decoding library on view
    // of the trigger pattern TDC.)
    const CsEvent *csEv = const_cast<TEv*>(this)->ptrEvt();
    const CS::Trigger *masterTrig = csEv->getDaqEvent().GetTT().GetTriggerMT();
    int masterTrigBit = masterTrig ? (int)masterTrig->GetBit() : -1;
    if (!(1<<masterTrigBit&tr)) {
      CsErrLog::msg(elError,__FILE__,__LINE__,"Master trigger (=%d) not in trigger pattern (=0x%x)",masterTrigBit,tr);
    }
    else {
      sprintf(c,"%2u",masterTrigBit); s += string(c)+" ";
      tr &= ~(1<<masterTrigBit);
    }
  }

  //                                               ***** REST of TRIGGER PATTERN
  for (int i=0; i<32; i++) {
    if ( (tr>>i)&1 ) { sprintf(c,"%2u",i); s += string(c)+" "; }
  }
  return s;
}


//! Print MC information block
void TEv::PrintMC(int mode)
{
  cout<<endl;
  cout<<"----------------------------------------------------"<<endl;
  if(mode == 0) {
    cout<<endl<<"\t\t\t Track block "<<endl<<endl;
    for(int i = 0; i < int(vecKine.size()); i++){
      TKine& k = vecKine[i];
      cout<<"Track # "<<i<<" ("<<k.Name()<<")  \t mother vertex "<<k.iVtxOrig<<"  \t daughter vertices : ";
      for(int j = 0; j < int(k.vecVtxMCRef.size()); j++) cout<<k.vecVtxMCRef[j]<<" ";
      cout<<endl;
    }
    cout<<endl<<"\t\t\t Vertex block "<<endl<<endl;
    for(int i = 0; i < int(vecVtxMC.size()); i++){
      TVtxMC& v = vecVtxMC[i];
      cout<<"Vertex # "<<i<<" X = "<<v.V(0)<<"  \t mother track "<<v.iOrig<<"  \t daughter tracks :";
      for(int j = 0; j < int(v.vecTrkRef.size()); j++) cout<<v.vecTrkRef[j]<<" ";
      cout<<endl;
    }
  }

  if(mode == 1) {
    cout<<endl<<"\t\t\t MC Track <--> Reconstructed tracks "<<endl<<endl;
    for(int i = 0; i < int(vecKine.size()); i++){
      TKine& k = vecKine[i];
      cout<<"Track # "<<i<<" ("<<k.Name()<<") \t  P = "<<setw(10)<<1./k.Pinv();
      if(!k.sTrackID().empty()){
	cout<<" <--> ";
	set<unsigned int, less<unsigned int> >::iterator iID;
	for(iID = k.sTrackID().begin(); iID != k.sTrackID().end(); iID++){
	  cout<<(*iID)<<" ";
	}
      }
      cout<<endl;
    }
  }

  cout<<"----------------------------------------------------"<<endl;
}

//! Print reconstructed tracks related information
void TEv::PrintRecTracks(int mode)
{
  cout<<endl<<"\t\t\t Reconstructed tracks "<<endl<<endl;
  list<TTrack>::iterator it;
  for(it = listTrack.begin(); it != listTrack.end(); it++){
    TTrack& t = (*it);
    cout<<"Track ID "<<t.Id<<" NDFs = "<<t.NDFs<<" Chi2 = "<<setw(8)<<t.Chi2tot<<" IKine = "<<t.IKine<<endl;
  }

}


/*
  MM cluster position correction for magnetic field effects and geometry of the detector
*/
void TEv::CorrectMMClusters(TTrack& track)
{
  const TSetup &setup = TSetup::Ref();

  if( setup.vIplFirst().size() < 2 ) return; // no zones!

  // first and last plane in zone 0 (target -- SM1)
  const int& ipFirstZone = setup.vIplFirst()[0];
  const int& ipLastZone  = setup.vIplLast()[0];

  // find first and last hit in group 0
  THit *firstHit=0, *lastHit=0;
  THlx firstHlx, lastHlx;

  //--- loop on hits to have first and last hit in the zone
  vector<THit*> vMMHits;
  list<int>::iterator ih = track.lHitPat.begin();
  list<int>::iterator ip = track.lPlnRef.begin();
  for( ; ih != track.lHitPat.end(); ++ih,++ip){
    if( *ip > ipLastZone ) break;
    if( *ih < 0 || *ip < ipFirstZone ) continue;
    if( !firstHit ){
      firstHit = &vecHit[*ih];
    }
    lastHit = &vecHit[*ih];
    string TBname = vecHit[*ih].DetRef().PtrDet()->GetTBName().substr(0,2);
    if (TBname == "MP"){ // add to the if :" (TBname == "MM") || " for consider magnetic fiel correction
      vMMHits.push_back(&vecHit[*ih]);
    }
  }

  // exit if no MM hit found
  if (vMMHits.size() == 0 ){
    //cout<<"TTrack::CorrectMMClusters: No MM hits, exiting."<<endl;
    return;
  }

  // calculate helices for first and last hit
  if(firstHit && lastHit){
    firstHlx(0) = firstHit->DetRef().X(0);  // take position of first hit along beam axis
    track.Hfirst.Extrapolate(firstHlx);
    //printf("    first hit in %s (%e,%e,%e)\n",firstHit->DetRef().Name.data(),firstHlx(0),firstHlx(1),firstHlx(2));
    lastHlx(0) = lastHit->DetRef().X(0);  // take position of last hit along beam axis
    track.Hfirst.Extrapolate(lastHlx);
    //printf("    last hit in %s (%e,%e,%e)\n",lastHit->DetRef().Name.data(),lastHlx(0),lastHlx(1),lastHlx(2));
  }else{
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "TTrack::CorrectMMClusters: first hit or last hit in zone not found.");
    return;
  }

  double x0[] = {firstHlx(0),firstHlx(1),firstHlx(2)};
  double x1[] = { lastHlx(0), lastHlx(1), lastHlx(2)};
  double dX = x1[0] - x0[0];
  double xh[3], w0=1.0,w1=0.0;
  bool error;

  //printf("    x0=(%e,%e,%e)  x1=(%e,%e,%e)  dz=%e\n",x0[0],x0[1],x0[2],x1[0],x1[1],x1[2],dX);
  //cout<<"TTrack::CorrectMMClusters: entering second loop "<<endl;

  //--- loop on MM hits for correcting cluster positions
  vector<THit*>::iterator it;
  for( it=vMMHits.begin(); it!=vMMHits.end(); ++it){

    // find hit position interpolating linearly between first and last helix in zone
    const TDetect &d = (*it)->DetRef();
    xh[0] = d.X(0);
    if( dX > 0 ){
      w0 = (x1[0] - xh[0])/dX;
      w1 = 1. - w0;
    }else if ( dX == 0){ // if there is only one hit in the zone
      w0 = 1.;
      w1 = 0.;
    }else{ // negative dX should not be possible
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "TTrack::CorrectMMClusters: dX<0 contradicts plane ordering along beam axis");
    }
    for(int i=1;i<3;++i) xh[i] = w0*x0[i] + w1*x1[i];

    // correct cluster position
    CsDetector *det = d.PtrDet();
    CsCluster *c = (*it)->PtrClus();
    double u = c->getU();
    double up = det->getCorrU(c,xh[1]*10.,xh[2]*10.,0.0 /*tT0*/,error);

    //printf(" ++++ Found hit on %s u = %e up = %e  hit.u = %e\n",det->GetTBName().data(),u,up,(*it)->U);

    // update hit and cluster
    (*it)->u = up/10.;
    c->setU(up);
  }

}
