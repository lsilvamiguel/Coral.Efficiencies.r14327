// $Id: Main.traf.cc,v 1.18 2010/09/07 18:27:00 tnagel Exp $
/*!
   \file    Main.traf.cc
   \brief   Runs traffic/trafdic on modified option file to generate alignment tree
   \author  Hugo Pereira
   \version $Revision: 1.18 $
   \date    $Date: 2010/09/07 18:27:00 $
*/


#include "DaqDataDecoding/Exception.h"
#include "Coral.h"
#include "CsInit.h"
#include "CsOpt.h"
#include "CsGeom.h"
#include "CsField.h"
#include "CsEvent.h"
#include "CsRegistrySing.h"
#include "CsErrLog.h"
#include "CsZone.h"

#include "CsDetector.h"
#include "CsDriftChamberDetector.h"
#include "CsDriftTubeDetector.h"
#include "CsStrawTubesDetector.h"
#include "CsTrack.h"
#include "CsCluster.h"
#include "CsDigit.h"
#include "CsHelix.h"
#include "CsHistograms.h"
#include "Defs.h"

#include <signal.h>

using namespace std;
using namespace CLHEP;

/*! \fn void InitAlignmentTree( void );
  \brief Init alignment Tree
*/
  void InitAlignmentTree( void );

/*! \fn bool FillAlignmentTree( unsigned int evt=0 )
  \brief fill the tree according to what's in CsEvent::Instance()
  \param evt the event index. It is stored in the tree
*/
  bool FillAlignmentTree( unsigned int evt=0 );

  void RequestInteruption(int dummy); //!< to catch interruptions, try waiting to the end of event
  void ForceInteruption(int dummy);   //!< to catch interruptions, try waiting to the end of event
  bool InteruptionRequested;          //!< true when interruption is requested

  bool do_alignment = true;   //!< if false, tree is not booked nor filled
  bool require_cop = false;   //!< if true, only tracks with cop are stored to the tree
  bool magnets_on = false;    //!< control the structure of the tree
  double chi2_cut = -1;       //!< Cut on chi2

/*! \fn int main( int argc, char *argv[] )
  \brief main file. 
  \param argc number of arguments
  \param argv arguments must at least contain option file
*/

int main( int argc, char *argv[] ) {


  try
  {                   
    InteruptionRequested = false;

  	// redirect interuptions 
  	signal(SIGTERM,RequestInteruption);
  	signal(SIGQUIT,RequestInteruption);
  	signal(SIGABRT,RequestInteruption);
  	signal(SIGALRM,RequestInteruption);

    // Package Initialization ===
    Coral::init( argc, argv );    
    
    CsEvent* event = CsEvent::Instance();
    int nevt=0;
    int refresh; 
    if( ! CsOpt::Instance()->getOpt("events", "refresh rate", refresh) ) refresh = 1;
   
    // Check for alignment options
    do_alignment = CsOpt::Instance()->getOpt("main", "do alignment");
    magnets_on  = CsOpt::Instance()->getOpt("main", "magnets on");
    require_cop = CsOpt::Instance()->getOpt("main", "require cop");
    if( !CsOpt::Instance()->getOpt("main", "chi2 cut", chi2_cut) ) chi2_cut = -1;
 
    // Dump options
    cout << endl;
    cout << "main - INFO: ================================\n";
    cout << "main - INFO: = Magnets are asumed to be " << ((magnets_on)?"on ":"off") << endl;
    if( !magnets_on ){
      cout << "main - INFO: = RequireCop forced to false " << endl;
      require_cop = false;
    } else 
    cout << "main - INFO: = RequireCop is " << ((require_cop)?"on":"off")<< endl;
    cout << "main - INFO: = Chi2Cut: ";
      if( chi2_cut > 0 ) cout << chi2_cut << endl;
      else cout << "none.\n";
    cout << "main - INFO: ================================\n";
    cout << endl;

    int ngood = 0;

    // Loop on events ===   
    while( (!InteruptionRequested) && event->getNextEvent() ) { 
      bool accept = true;  
      
      if( do_alignment ) {     
        if( nevt == 0 )  InitAlignmentTree( );           
        accept = FillAlignmentTree( nevt );
      } 
      
      if( accept ) ngood++;
      if( !(CsRegistrySing::Instance()->callEoeMethods()) ) break;
      if((++nevt)%refresh == 0 ) cout << "Event: " << nevt << " (" << ngood << ")" << endl;
    }

    // End session ===
    CsRegistrySing::Instance()->callEndMethods();
    CsErrLog::Instance()->dump( elDebugging );
    	               
  }  
  catch(std::exception &e ){ cerr << "Exception:\n" << e.what() << endl; }
  catch(const char *e) { std::cerr << "Exception:\n" << e << "\n"; }
  catch( ... ){ cerr << "Unknown exception!" << endl; }

  return 0;
}

#include "TTree.h"
#include "TMath.h"
#include "TObjString.h"
TTree* T_all_;            //!< Alignment tree for all detectors used by traffic

//===============
// Alignment tree
//===============
static unsigned int T_evt;     //!< Event number
static unsigned int T_trigMsk; //!< Event trigger mask
static unsigned int T_zone;    //!< Event trigger mask
static unsigned int T_cmlt;    //!< track multiplicity
static double T_chi2;          //!< track chisquare
static int    T_ndf;           //!< track number of degree of freedom
static double T_prob;          //!< track chisquare propability
static double T_cop;           //!< charge over momentum
static double T_meanT;         //!< track time

static unsigned int T_nDets;           //!< number of fired detectors
static unsigned int T_detVect[NPLAN];  //!< detector ID
static double T_uVect[NPLAN];  //!< u_cluster, coordinate perp to the wire
static double T_vVect[NPLAN];  //!< v_cluster, coordinate along the wires (2nd proj for pixel detectors)
static double T_duVect[NPLAN]; //!< u_cluster - u_track
static double T_dvVect[NPLAN];    // v_cluster - v_track (for pixel detector)

static double T_rVect[NPLAN];  //!< u_cluster - u_wire     (drift like detectors only)
static double T_tVect[NPLAN];  //!< t_cluster cluster time (drift like detectors only) 

//!<  branches specific to magnets_on_ = false
static double T_xLx;                //!< track X
static double T_yLx;                //!< track Y
static double T_zLx;                //!< track Z
static double T_txLx;               //!< track dX/dZ
static double T_tyLx;               //!< track dY/dZ

//!<  branches specific to magnets_on_ = true
static double T_xVect[NPLAN];  //!< track X         
static double T_yVect[NPLAN];  //!< track Y         
static double T_zVect[NPLAN];  //!< track Z         
static double T_txVect[NPLAN]; //!< track dX/dZ     
static double T_tyVect[NPLAN]; //!< track dY/dZ     
static double T_dzVect[NPLAN]; //!< Extrapolation length 
static double T_BVect[NPLAN];  //!< Magnetic field along the wires

//____________________________________________________
void InitAlignmentTree( void ) 
{
  static bool firstCall = true;
  if( !firstCall ) return;
  
  CsHistograms::SetCurrentPath("/CsAlignment");
  cout << "InitAlignmentTree ... \n"; 
  string TreeName( "T_align_" );
  T_all_ = new TTree( TreeName.c_str(), TreeName.c_str() );
  T_all_->Branch("T_evt",     &T_evt,       "T_evt/i",     BUFFER_SIZE);
  T_all_->Branch("T_trigMsk", &T_trigMsk,   "T_trigMsk/i", BUFFER_SIZE);
  T_all_->Branch("T_zone",    &T_zone,      "T_zone/i",    BUFFER_SIZE);
  T_all_->Branch("T_cmlt",    &T_cmlt,      "T_cmlt/i",    BUFFER_SIZE);
  T_all_->Branch("T_chi2",    &T_chi2,      "T_chi2/D",    BUFFER_SIZE);
  T_all_->Branch("T_ndf",     &T_ndf,       "T_ndf/I",     BUFFER_SIZE);
  T_all_->Branch("T_prob",    &T_prob,      "T_prob/D",    BUFFER_SIZE);
  T_all_->Branch("T_cop",     &T_cop,       "T_cop/D",     BUFFER_SIZE);
  T_all_->Branch("T_meanT",   &T_meanT,     "T_meanT/D",   BUFFER_SIZE);

  T_all_->Branch("T_nDets",    &T_nDets,       "T_nDets/i", BUFFER_SIZE);
  T_all_->Branch("T_detVect",  &T_detVect[0],  "T_detVect[T_nDets]/I", BUFFER_SIZE);

  if( magnets_on ) {
    //!<  branches specific to magnets_on_ = true
    //!< Track parameters extrapolated on each detector
    T_all_->Branch("T_xVect",    &T_xVect[0],    "T_xVect[T_nDets]/D",   BUFFER_SIZE);  
    T_all_->Branch("T_yVect",    &T_yVect[0],    "T_yVect[T_nDets]/D",   BUFFER_SIZE);  
    T_all_->Branch("T_zVect",    &T_zVect[0],    "T_zVect[T_nDets]/D",   BUFFER_SIZE);  
    T_all_->Branch("T_txVect",   &T_txVect[0],   "T_txVect[T_nDets]/D",  BUFFER_SIZE);  
    T_all_->Branch("T_tyVect",   &T_tyVect[0],   "T_tyVect[T_nDets]/D",  BUFFER_SIZE);  
    T_all_->Branch("T_dzVect",   &T_dzVect[0],   "T_dzVect[T_nDets]/D",  BUFFER_SIZE);  
    T_all_->Branch("T_BVect",    &T_BVect[0],    "T_BVect[T_nDets]/D",   BUFFER_SIZE);  
  } else {
    //!<  branches specific to magnets_on_ = false
    //!<  straight track parameters
    T_all_->Branch("T_xLx",     &T_xLx,       "T_xLx/D", BUFFER_SIZE);
    T_all_->Branch("T_yLx",     &T_yLx,       "T_yLx/D", BUFFER_SIZE);
    T_all_->Branch("T_zLx",     &T_zLx,       "T_zLx/D", BUFFER_SIZE);
    T_all_->Branch("T_txLx",    &T_txLx,      "T_txLx/D", BUFFER_SIZE);
    T_all_->Branch("T_tyLx",    &T_tyLx,      "T_tyLx/D", BUFFER_SIZE);
  }
    
  // Hits in detector
  T_all_->Branch("T_uVect",    &T_uVect[0],    "T_uVect[T_nDets]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_vVect",    &T_vVect[0],    "T_vVect[T_nDets]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_duVect",   &T_duVect[0],   "T_duVect[T_nDets]/D",  BUFFER_SIZE);  
  T_all_->Branch("T_dvVect",   &T_dvVect[0],   "T_dvVect[T_nDets]/D",  BUFFER_SIZE);  

  // This branches are filled with non 0 values for drift-like detectors only
  T_all_->Branch("T_rVect",    &T_rVect[0],    "T_rVect[T_nDets]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_tVect",    &T_tVect[0],    "T_tVect[T_nDets]/D",   BUFFER_SIZE);  

  T_nDets = 0;
  firstCall = false;
  
  cout << "InitAlignmentTree. done.\n"; 
  return;
}
  

//____________________________________________________________________
#include <functional>
//! To sort helices according to their z position
struct sortLx : public binary_function< CsHelix, CsHelix, bool > 
{ 
  bool operator() ( CsHelix lx0, CsHelix lx1 ) 
  { return ( lx0.getZ() < lx1.getZ() ); } 
};

//____________________________________________________
#include <list>
bool FillAlignmentTree( unsigned int evt )
{
  T_evt = evt; 
  T_trigMsk = (unsigned int) CsEvent::Instance()->getTriggerMask();

  // Get Tracks - Check number of tracks
  list< CsTrack* > tr = CsEvent::Instance()->getTracks();
  if( !tr.size() ) return false;

  for( list< CsTrack* >::iterator It = tr.begin(); It != tr.end(); It++ ) {
    if( !(*It)->getHelices().size() ) continue;
    
    // Put helices into a list, and sort along Z
    list< CsHelix > lx;
    for( unsigned int i=0; i<(*It)->getHelices().size(); i++ ) {
      lx.push_back( (*It)->getHelices()[i] );
      lx.sort( sortLx() );
    }
    
    // Get track parameters  
    T_cop = lx.front().getCop();                   // helix Charge/Momentum
    if( require_cop && T_cop == 0 ) continue;      // check T_cop && cop required
    
    T_chi2 = (*It)->getChi2();                     // track chi2
    T_cmlt = (*It)->getClusters().size();          // track cluster size

    if( fabs(T_cop) < 0.0001 ) T_ndf = T_cmlt - 4; // not bridged track
    else T_ndf = T_cmlt - 5;                       // bridged track    
    if( chi2_cut > 0 && T_chi2/T_ndf > chi2_cut ) continue;  // Check chi2/ndf

    T_prob = TMath::Prob( T_chi2, T_ndf );         // track chi2 probability
    T_meanT = ( (*It)->hasMeanTime() ) ? (*It)->getMeanTime():0; // track mean time
    
    if( !magnets_on ) {
      T_xLx = lx.front().getX();
      T_yLx = lx.front().getY();
      T_zLx = lx.front().getZ();
      T_txLx = lx.front().getDXDZ();
      T_tyLx = lx.front().getDYDZ();
    }
    
    // Get List of zones associated to track
    list< CsZone* > zones = (*It)->getZones();
    T_zone = 0;
    for( list< CsZone* >::iterator Iz = zones.begin(); Iz != zones.end(); Iz++ ) T_zone += (1<<(*Iz)->getId());    
      
    // loop over detectors in the tracks
    T_nDets = 0;
    list< CsCluster* > c = (*It)->getClusters();
    list< CsCluster* >::iterator Ic;
    for( Ic = c.begin(); Ic != c.end(); Ic++ ) {
      
      // Get detector associated to cluster
      CsDetector *d = (*Ic)->getDetsList().front();
      HepMatrix iR = d->getRotWRSInv();
      double zD = d->getZcm();
      double u;
      double v;
      
      if( magnets_on ) {
        
        // Get helix closest to detector
        double dzMin = 0;
        CsHelix* lxMin = 0;
        for( list< CsHelix >::iterator Ilx = lx.begin(); Ilx != lx.end(); Ilx ++ )
        if( (!lxMin) || fabs( (*Ilx).getZ()-zD ) < fabs(dzMin) ) {
          dzMin = (*Ilx).getZ()-zD;
          lxMin = &(*Ilx );
        } else break;
        
        // Extrapolate lx to detector
        CsHelix *lxDet = new CsHelix();
        lxMin->Extrapolate( zD, *lxDet ); 
        
        // Fill detectorWise branches
        
        // Track position
        double x  = T_xVect[T_nDets]  = lxDet->getX();
        double y  = T_yVect[T_nDets]  = lxDet->getY();
        double z  = T_zVect[T_nDets]  = lxDet->getZ();
        double tx = T_txVect[T_nDets] = lxDet->getDXDZ();
        double ty = T_tyVect[T_nDets] = lxDet->getDYDZ();
        u = iR(1,1)*( x + tx*( zD - z) ) + iR(1,2)*( y + ty*( zD - z ) );
	v = iR(2,1)*( x + tx*( zD - z) ) + iR(2,2)*( y + ty*( zD - z ) );
        T_dzVect[T_nDets] = dzMin;
      
        // Magnetic field
        float bx, by, bz;
        CsGeom::Instance()->getCsField()->getField( x, y, z, bx, by, bz ); 
        T_BVect[T_nDets]= iR(2,1)*double(bx)+iR(2,2)*double(by);        
      
        delete lxDet;
        
      } else {
        
        // Straight trach extrapolation to detector
        u = iR(1,1)*( T_xLx + T_txLx*( zD-T_zLx) ) + iR(1,2)*( T_yLx + T_tyLx*( zD-T_zLx) );
	v = iR(2,1)*( T_xLx + T_txLx*( zD-T_zLx) ) + iR(2,2)*( T_yLx + T_tyLx*( zD-T_zLx) );
      
      }
      
      //==== detector
      T_detVect[T_nDets]=d->GetID();      //!< Detector id
      T_uVect[T_nDets]  = (*Ic)->getU();  //!< coordinate perp to the wire
      T_duVect[T_nDets] = (*Ic)->getU() - u; //!< residual perp to the wire

      T_vVect[T_nDets]  = (*Ic)->getV();  //!< coordinate along the wire
      T_dvVect[T_nDets] = (*Ic)->getV() - v; //!< residual along wire
            
      // T_rVect and T_tVect are filled only for drift like detectors
      CsDriftChamberDetector* dc = dynamic_cast<CsDriftChamberDetector*>(d);
      CsDriftTubeDetector*    dt = dynamic_cast<CsDriftTubeDetector*>(d);
      CsStrawTubesDetector*   st = dynamic_cast<CsStrawTubesDetector*>(d);
      if( dc || dt || st || d->hasDrift() ) {
        int wire  = (*Ic)->getDigitsList().front()->getAddress();
        double uw = d->getWirD() 
          + iR(1,1) * d->getDeltaXCorrection()
          + iR(1,2) * d->getDeltaYCorrection() 
          + wire * d->getWirP();
        T_rVect[T_nDets] = (*Ic)->getU()-uw;
        T_tVect[T_nDets] = (*Ic)->getDigitsList().front()->getDatum();
      } else {
        T_rVect[T_nDets] = 0;
        T_tVect[T_nDets] = 0;
      }
        
      T_nDets++;

      if( T_nDets >= NPLAN ) {
        cout << "FillAlignmentTree - ERROR: too many detectors in track. EoT skipped.\n";
        break;
      } // test on number of clusters
    }   // loop over clusters

    T_all_->Fill();

  } // loop over tracks
  return true;
}  

//____________________________________________________________________
void ForceInteruption(int dummy)
{
	printf("\nForceInteruption - INFO: Got signal (%i).\n",dummy);
	printf("ForceInteruption - INFO: Try to finish safely.\n");
  try {                   
    CsRegistrySing::Instance()->callEoeMethods();  // try finish the event
    CsRegistrySing::Instance()->callEndMethods();  // try finish the run
	  printf("\n\n\aterminated.\n");
  }  
  catch(std::exception &e ){ cerr << "Exception:\n" << e.what() << endl; }
  catch( ... ){ cerr << "Unknown exception!" << endl; }
  exit(0);

	return;
}

//____________________________________________________________________
//! to catch interruptions, try waiting to the end of event
void RequestInteruption(int dummy)
{
	printf("\nRequestInteruption - INFO: Got signal (%i).\n",dummy);
	printf("RequestInteruption - INFO: Try waiting for end of event to finish.\n");
  InteruptionRequested = true;
	return;
}
