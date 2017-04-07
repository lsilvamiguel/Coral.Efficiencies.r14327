// $Id: Traffic.cc 13509 2012-07-10 17:22:19Z ybedfer $

#include <csignal>
#include <cassert>
#include <iostream>

#include "CsInit.h"
#include "CsEvent.h"
#include "TConstants.h"
#include "Traffic.h"
#include "TDisplay.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "TTrack.h"

using namespace std;

int  UserSignal(0);
void SignalHandler(int isig)
{
  cout<<endl<<"=====>  User signal  # "<<isig<<" had been received"<<endl;
  UserSignal=isig;
};


/*!
  \file
  Constructor(s), destructor and trivial methods
  of the class Traffic
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/Traffic":
  i) Handle case where "Remode[1]==3", viz. see too that MC be read in at last,
  ii) More TWatches printed.
*/

Traffic* Traffic::address = 0; // init static pointer
int      Traffic::Nevt    = 0;

/*! 
  external HBOOK/ZEBRA initialization from CORAL
*/
extern "C" int inithbook_( int& );
   
//! Constructor (Traffic package initialisation)
Traffic::Traffic()
{
  if(address == 0){ // if Traffic object not yet created
    address = this;

    signal(SIGUSR1, &SignalHandler);
    signal(SIGUSR2, &SignalHandler);

    // reset counters
    nevs=0;  
    ntracks=0; ntracks_with_P=0;
    ntracks_beam=0; ntracks_beam_with_P=0;
    n_selected_MC=0; n_reconstructed_MC=0;

    cout<<endl<<endl
	<<"------>  TRAFFIC version "<<TConstants_TraFFiC_Version<<"  <------" 
	<<endl<<endl;
    
    // Get options for TraFFiC reconstruction 
    TOpt::getOptions();
    if(TOpt::Print[0] == 0) cout<<"         (Traffic printing is switched off)";
    cout<<endl;

    if(TOpt::ReMode[0] == 1) {
      cout<<endl<<endl
	  <<"------>  TRAFFIC track finding is OFF by TraF ReMode[0] option "<<endl<<endl;
    }

    if(TOpt::ReMode[0] >= 2) {
      cout<<endl<<endl
	  <<"------>  TRAFFIC is switched OFF by TraF ReMode[0] option "<<endl<<endl;
      return;
    }

    // Register for "end of job" call
    CsRegistry reg;
    reg.EOJRegistration( this );
    // Register for "end of event" call
    reg.EOERegistration( this );

    // Create Tracking setup object
    new TSetup;
    
    // create TDisplay object.
    // (but if the event display is OFF this is dummy object)
    if(TOpt::Graph[0] > 0){
      if(! CsInit::Instance()->initZebraDone()){
	int unit(0); (void) inithbook_(unit);  // this is CORAL Fortran function 
	cout<<"------>  ZEBRA had been initialized for HIGZ graphics"<<endl;
      }
    }
    new TDisplay;
  } // end if
 
  // do nothing if Traffic object already created

}


//! Destructor
Traffic::~Traffic()
{
  if(TDisplay::Ptr() != NULL) delete TDisplay::Ptr(); // delete TDisplay 
  if(TSetup::  Ptr() != NULL) delete TSetup  ::Ptr(); // delete TSetup
  address = 0;
}



// "End Of Event" method
bool Traffic::eoe() {

  if(TOpt::ReMode[0] >= 2) return(true); // Traffic is switched OFF completely


  CsEvent* csev = CsEvent::Instance();

  bool print = false;

  if(print){ // EOF statistics
    cout<<"--------> N CsMCHits   = "<<csev->getMCHits().  size()<<endl;
    cout<<"--------> N CsDigits   = "<<csev->getDigits().  size()<<endl;
    cout<<"--------> N CsClusters = "<<csev->getClusters().size()<<endl;
    list<CsTrack*> lCsTr = csev->getTracks();
    cout<<"--------> N CsTracks   = "<<lCsTr.size()<<endl;
    list<CsTrack*>::iterator it;
    for( it = lCsTr.begin(); it != lCsTr.end(); it++){
      cout<<(*it)->getId()<<"   ";
    }
    cout<<" <-------- IDs of CsTracks seen from Traffic EOE "<< endl;
  }


  // Get information for End-Of-Event actions

  // Special case where MC info has been kept hidden (at least partially) so far
  int hidden_MC = TOpt::ReMode[1];
  if (hidden_MC==3 || hidden_MC==4)
    TOpt::ReMode[1] = 0;  // => Reset "ReMode[1]" so as to retrieve it now

  if (TOpt::ReMode[1]==0 && TEv::Ptr()!=NULL) delete TEv::Ptr();
  if (TEv::Ptr()==NULL) {
    if (TOpt::ReMode[1]==2)
      // No TEv object existing while we have Traffic modularity = 2, i.e.
      // TEv being kept alive for all of the reconstruction. => Tracking has not
      // been activated for the current CsEvent. This can be because we are
      // skipping the BoS (I don't know of any other reason anyway). => There's
      // nothing but an empty event to show, nor is there any statistics to be
      // updated. => Exit.
      return true;
    new TEv; // Create TEv object instance (if not yet exists)
  }

  TEv &ev = TEv::Ref();
  const TSetup &setup = TSetup::Ref();

  if (TOpt::ReMode[1]==2) {
    if (ev.IsMC())
      CsErrLog::mes(elWarning,"ReMode[1]==2 => MC Info is not reliable");
  }
  else {
    list<CsTrack*> lTr = csev->getTracks();
    ev.ImportTracks(lTr, "cp");
    list<CsCluster*> lClust = csev->getClusters();
    TOpt::ReMode[1] = hidden_MC;
    ev.ImportClusters(lClust);
  }
  
  Traffic::Ref().Stopwatch.Start(4);

  // ********** STATISTICS **********
  // #Tracks: Overall and for this->event
  static int nTAll = 0, nT0x1 = 0, nT0x3 = 0, nT0x6 = 0, nT0xe = 0, nT0xf = 0;
  int        ntAll = 0, nt0x1 = 0, nt0x3 = 0, nt0x6 = 0, nt0xe = 0, nt0xf = 0;
  // Idem for tracks associated to a TKine
  static int nKAll = 0, nK0x1 = 0, nK0x3 = 0, nK0x6 = 0, nK0xe = 0, nK0xf = 0;
  int        nkAll = 0, nk0x1 = 0, nk0x3 = 0, nk0x6 = 0, nk0xe = 0, nk0xf = 0;
  nevs++;
  list<TTrack>::const_iterator it = ev.lTrack().begin();
  for (; it != ev.lTrack().end(); it++) {
    const TTrack &t = *it;

    if (t.H('d')(0)>setup.TargetCenter[0]) { // Spectrometer tracks
      ntracks++;
      if (t.Pinv()!=0 &&  // Require P
	  // ...and exclude 0xc (which are granted P = beam's P from options
	  // so their extension through the mu-absorber can be fitted w/ MS ON).
	  t.GetType()!=0xc) {
	if (t.NGroups()>1) ntracks_with_P++;  
	//#define DEBUG_TraFDic
#ifdef DEBUG_TraFDic
	if (TOpt::Print[8]&0x4 && ev.IsMC()) {
	  printf(" #%d(%d,%x),",t.GetId(),t.GetIKine(),t.GetType());
	}
#endif
	if (TOpt::Print[8]&0x2) {  // ***** UPON OPTION...
	  // ***** PRINT #TTrack's PER TYPE for CURRENT EVENT, CUMULATED *****
	  ntAll++; bool associated = t.GetIKine()>=0; if (associated) nkAll++;
	  if (t.GetType()==0x1)       { nt0x1++; if (associated) nk0x1++; }
	  if ((t.GetType()&0x3)==0x3) { nt0x3++; if (associated) nk0x3++; }
	  if ((t.GetType()&0x6)==0x6) {
	    nt0x6++; if (associated) nk0x6++;
	    if (t.GetType()&0x8) { nt0xe++; if (associated) nk0xe++; }
	  }
	  if (t.GetType()==0xf) { nt0xf++; if (associated) nk0xf++; }
	}
      }
    }
    if (t.H('d')(0)<setup.TargetCenter[0]) { // Beam tracks
      ntracks_beam++;
      if (t.Pinv()!=0) ntracks_beam_with_P++;  
    }  
  }
  // ***** UPON OPTION: PRINT *****
  if (TOpt::Print[8]&0x1) 
    cout<<"TTracks/Evt: "<<(double)ntracks_with_P/nevs<<endl;
  if (TOpt::Print[8]&0x2) {
    nTAll += ntAll; nT0x1 += nt0x1;
    nT0x3 += nt0x3; nT0x6 += nt0x6; nT0xe += nt0xe; nT0xf += nt0xf;
    printf("#TTracks: %2d  %2d %2d %2d %2d %2d    %5d  %5d %5d %5d %5d %5d\n",
	   ntAll,nt0x1,nt0x3,nt0x6,nt0xe,nt0xf,
	   nTAll,nT0x1,nT0x3,nT0x6,nT0xe,nT0xf);
    if (TOpt::Print[8]&0x4) {
      nKAll += nkAll; nK0x1 += nk0x1;
      nK0x3 += nk0x3; nK0x6 += nk0x6; nK0xe += nk0xe; nK0xf += nk0xf;
      printf("#TTracks: %2d  %2d %2d %2d %2d %2d    %5d  %5d %5d %5d %5d %5d\n",
	     nkAll,nk0x1,nk0x3,nk0x6,nk0xe,nk0xf,
	     nKAll,nK0x1,nK0x3,nK0x6,nK0xe,nK0xf);
    }
  }


  // Internal consistency checks
  if(TOpt::ReMode[27] > 0) ev.Checks();
  // simple V0 search (like a final tracking crosscheck)
  //if(TOpt::ReMode[28] > 0) ev.SearchV0();

  // monitoring histograms
  ev.Monitor();

  // Event Display
  // (following drawing methods do nothing if the event display is off)
 draw_again:
  TDisplay& display = TDisplay::Ref();
  int iret(0);
  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
  unsigned int evTrig = ev.TrigMask()&allTrigs; // Cut away trailing bits
  if (TOpt::Graph[0] &&
      (TOpt::iCut[0]==0 ||
       (evTrig&TOpt::iCut[0]) && // Require trig to be included in selection...
       (evTrig&(~TOpt::iCut[0]))==0))  // ...AND strictly included
    iret=display.Draw(0);    // Draw everything, wait menu click
  if(TOpt::Graph[0] != 0 && UserSignal == SIGUSR2) { // if in "movie" mode
    UserSignal = 0 ;
    TOpt::Graph[7] *= -1;
    goto draw_again;
  }

  Traffic::Ref().Stopwatch.Stop(4);

  // final EOE print
  int n = TOpt::Print[9];
  if(n != 0) {
    if(TOpt::Graph[0] > 0) n=1;
    if( ++Nevt%n == 0 ) {
      cout<<"------- Processed by Traffic:  Run # " <<  csev->getRunNumber() 
	  <<"  Event # " <<  csev->getEventNumberInRun()<<"  ("<<Nevt<<")";
      if(csev->isADataEvent()) cout<<" taken "<<csev->getEventTime();
      cout<<endl;
    }
  }

  // Cleaning
  if(TEv::Ptr() != NULL) delete TEv::Ptr(); // delete TEv object at the and-of-event 

  if (iret==1)              return(false); // if "exit" was requested in the Event Display
  if(UserSignal == SIGUSR1) return(false); // terminate by user signal
  return( true );
}



// "End of Job" method 
bool Traffic::end() {

  cout.setf(ios::fixed,ios::scientific);
  cout.unsetf(ios::showpos);
  cout<<setprecision(7);
  
  if(TOpt::ReMode[0] == 0 && nevs > 0) { // if Traffic is ON
    cout<<"-----------> Traffic End-Of-Job statistics"<<endl;
    cout<<"Time in Pre-Pattern   = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(1)/nevs<<"  sec/ev"<<endl;
    if (TOpt::ReMode[11]==2){   // Alternative PrePatter #2
      cout<<" ... in Proj          = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(5)/nevs<<"  sec/ev"<<endl;
      cout<<" ... in Space         = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(6)/nevs<<"  sec/ev"<<endl;
      cout<<" ... Cleaning         = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(7)/nevs<<"  sec/ev"<<endl;
      cout<<" ... Paraxial         = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(16)/nevs<<"  sec/ev"<<endl;
      for (int zone = 0; zone<3; zone++)
	cout<<" ... in Zone "<<zone<<"        = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(17+zone)/nevs<<"  sec/ev"<<endl;
      cout<<" ... in Zones 3|4     = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(20)/nevs<<"  sec/ev"<<endl;
    }
    cout<<"Time in Bridging      = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(2)/nevs<<"  sec/ev"<<endl;
    if (TOpt::ReMode[12]==2) {  // Alternative Bridging #2
      if (CsInit::Instance()->IsAMonteCarloJob() &&
	  TOpt::ReMode[1]!=3) {  // Otherwise no MC info to tell good from wrong
	cout<<" ... SM1              = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(9)/nevs<<"  sec/ev"<<endl;
	cout<<" ... SM2              = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(10)/nevs<<"  sec/ev"<<endl;
	cout<<" ... SM1 wrong        = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(11)/nevs<<"  sec/ev"<<endl;
	cout<<" ... SM2 wrong        = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(12)/nevs<<"  sec/ev"<<endl;
      }
      else {
	cout<<" ... SM1              = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(11)/nevs<<"  sec/ev"<<endl;
	cout<<" ... SM2              = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(12)/nevs<<"  sec/ev"<<endl;
      }
      cout<<" ... Fitting          = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(13)/nevs<<"  sec/ev"<<endl;
    }
    cout<<"Time in Track Fit     = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(3)/nevs<<"  sec/ev"<<endl;
    cout<<"Time in End-of-Event  = "<<setprecision(4)<<setw(7)<<Stopwatch.SumTime(4)/nevs<<"  sec/ev"<<endl;
    float tsum = 0;
    tsum += Stopwatch.SumTime(1);
    tsum += Stopwatch.SumTime(2);
    tsum += Stopwatch.SumTime(3);
    tsum += Stopwatch.SumTime(4);
    cout<<"Total time in TRAFFIC = "<<setprecision(4)<<setw(7)<<tsum/nevs<<"  sec/ev"<<endl;
    cout<<"Total number of TRAFFIC beam  tracks / ev          = "<<double(ntracks_beam)/nevs<<endl;
    cout<<"Total number of TRAFFIC event tracks / ev          = "<<double(ntracks)     /nevs<<endl;
    cout<<"Number of TRAFFIC beam  tracks with momentum / ev  = "<<double(ntracks_beam_with_P)/nevs<<endl;
    cout<<"Number of TRAFFIC event tracks with momentum / ev  = "<<double(ntracks_with_P)     /nevs<<endl;
    if(n_selected_MC !=0){
      cout<<"Overal track finding efficiency in 2-50 GeV range = "
	  <<100.0*n_reconstructed_MC/n_selected_MC<<" %"<<endl; 
    }
      
  }

  if(TEv::Ptr() != NULL) delete TEv::Ptr(); // delete TEv at the end-of-job
  if(address    != NULL) delete this; // if Traffic still alive, commit suicide :-( 
  return( true );
}
