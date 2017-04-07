// $Id: Main.traf.cc,v 1.13 2010/09/07 18:27:02 tnagel Exp $

#include "DaqDataDecoding/Exception.h"
#include "Coral.h"
#include "CsInit.h"
#include "CsOpt.h"
#include "CsGeom.h"
#include "CsEvent.h"
#include "CsRegistrySing.h"
#include "CsErrLog.h"

#include "CsDetector.h"
#include "CsDriftChamberDetector.h"
#include "CsDriftTubeDetector.h"
#include "CsStrawTubesDetector.h"
#include "CsTrack.h"
#include "CsCluster.h"
#include "CsHelix.h"
#include "CsHistograms.h"
#include "Defs.h"

#include <signal.h>

#include <TObjString.h>

using namespace std;
using namespace CLHEP;

// Calibration Tree is for Residuals, RTRelation, efficiency, ...
void InitCalibrationTree( void );               
bool FillCalibrationTree( unsigned int evt=0 );

void RequestInteruption(int dummy);
void ForceInteruption(int dummy);
bool InteruptionRequested;
bool require_cop = false;   //!< if true, only tracks with cop are stored to the tree
int require_vertex = 0;   //!< if > 0, only events with primary vertex are stored to the tree. Vertex has to have at least require_vertex tracks
bool do_calibration = true;   //!< if false, tree is not booked nor filled
double chi2_cut = -1;       //!< Cut on chi2

//______________________________________________________________________
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

    // Check for calibration options
    int nevt=0;
    int ngood = 0;
    int refresh; 
    if( ! CsOpt::Instance()->getOpt("events", "refresh rate", refresh) ) refresh = 1;
    do_calibration = CsOpt::Instance()->getOpt("main", "do calibration");
    require_cop = CsOpt::Instance()->getOpt("main", "require cop");
    if( !CsOpt::Instance()->getOpt("main", "require vertex", require_vertex) ) require_vertex = 0;
    if( !CsOpt::Instance()->getOpt("main", "chi2 cut", chi2_cut) ) chi2_cut = -1;
 
    // Dump options
    cout << endl;
    cout << "main - INFO: ================================\n";
    cout << "main - INFO: = RequireCop is " << ((require_cop)?"on":"off")<< endl;
    cout << "main - INFO: = RequireVertex is " << ((require_vertex)?"on":"off")<< endl;
    if(require_vertex)
	cout << "main - INFO: = MinimumNuberOfTracksInVertex: "<< require_vertex<< endl;
    cout << "main - INFO: = Chi2Cut: ";
      if( chi2_cut > 0 ) cout << chi2_cut << endl;
      else cout << "none.\n";
    cout << "main - INFO: ================================\n";
    cout << endl;

    // Loop on events ===   
    while( (!InteruptionRequested) && event->getNextEvent() ) { 
      
      bool accept = true;       
      
      if( do_calibration ){ 
        if( nevt == 0 ) InitCalibrationTree(); 
        accept = FillCalibrationTree( nevt );
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
  catch( ... ){ cerr << "Unknown exception!" << endl; }

  return 0;
}

//=================
// Calibration tree
//=================
 
vector< CsDetector* > prDets_;   // list of detectors to be probed
vector< double > t0_;            // detector t0 (for drift like detectors) 0 otherwise
vector< bool > isOff_;           // detectors off for traffic

#include <TTree.h>
#include <TMath.h>
TTree* T_eff_ = 0;            // Alignment tree for all detectors used by traffic

//============
// Calibration tree
//============
static unsigned int T_evt;     // Event number
static unsigned int T_trigMsk; // Event trigger mask
static unsigned int T_zone;    // Event zone mask
static unsigned int T_cmlt;    // track multiplicity
static double T_chi2;          // track chisquare
static int    T_ndf;           // track number of degree of freedom
static double T_prob;          // track chisquare propability
static double T_cop;           // charge over momentum
static double T_meanT;         // track time

static unsigned int T_nDets;           // number of fired detectors
static unsigned int T_detVect[NPLAN];  // detector ID
static unsigned int T_inActive[NPLAN];  // 1 if track is in detector active area, 0 otherwise
static unsigned int T_fnd[NPLAN];     // 1 if a close enough cluster was found in the detector, is in detector active area, 0 otherwise
static double T_duVect[NPLAN];    // u_cluster - u_track

static double T_dvVect[NPLAN];    // v_cluster - v_track (for pixel detector)

static double D_U[NPLAN];    // pixel coordinates
static double D_V[NPLAN];    // 

static double T_rVect[NPLAN];     // u_cluster - u_wire     
static double T_tVect[NPLAN];     // t_cluster cluster time 

static double T_xVect[NPLAN];  // track X         
static double T_yVect[NPLAN];  // track Y         
static double T_zVect[NPLAN];  // track Z         
static double T_txVect[NPLAN]; // track dX/dZ     
static double T_tyVect[NPLAN]; // track dY/dZ     

static double T_uLx[NPLAN];  // track U         
static double T_vLx[NPLAN];  // track V         
static double T_tuLx[NPLAN]; // track dU/dZ         
static double T_tvLx[NPLAN]; // track dV/dZ         

static double T_dzVect[NPLAN]; // Extrapolation length (helix - detector) 
static double T_BVect[NPLAN];  // Magnetic field along the wires

//____________________________________________________
bool Match( string s1, string s2 )
{
  if( s2.size() > s1.size() ) return false;
  bool accept = true;
  for( unsigned int i=0; i < s2.size(); i++ ) 
  if( s2[i] != '*' && s2[i] != s1[i] ) {
    accept = false;
    break;
  }
  return accept;
}

//____________________________________________________
void InitCalibrationTree() 
{
  static bool firstCall = true;
  if( !firstCall ) return;
  
  // init detector vectors ===
  list<CsDetector*> dets = CsGeom::Instance()->getDetectors(); 
  prDets_.clear();
  t0_.clear();
  
  isOff_.clear();
   
  // Load Traffic desactivated detectors
  CsOpt* opt = CsOpt::Instance();
  list< string > detOff; detOff.clear();
  opt->getOpt( "TraF","DetNameOff", detOff );

  // init string of desactivated detectors TBNames
  string prDetsNames = "";
  
  // fill list of  probed detectors ===
  cout << endl << "=== Calibration tree === " << endl;
  list< string > TBN;
  opt->getOpt( "main", "do calibration", TBN );
  for( list< string >::iterator I = TBN.begin(); I != TBN.end(); I++ ) {
    string TBName = *I;
    bool found = false;
    list< CsDetector* >::iterator Id;
    for( Id = dets.begin(); Id != dets.end() && found == false; Id++ ) 
    if( Match( (*Id)->GetTBName(), TBName ) ) {
            
      // only for traffic desactivated detectors ===
      bool isOff = false;
      for( list< string >::iterator IDO = detOff.begin(); IDO != detOff.end(); IDO++ ) {
        string TB( (*Id)->GetTBName(), 0, (*IDO).size() );
        if( (*IDO ) != TBName && (*IDO) != TB ) continue;
        isOff = true;
      }
      
      if( !isOff ) continue;
      
      cout << "InitCalibrationTree - INFO: add Det " << (*Id)->GetTBName()
        << " to prDets_. "
        << "<off for traffic>.";

      // fill detector Info vectors ===
      CsDriftChamberDetector* dc = dynamic_cast<CsDriftChamberDetector*>(*Id);
      CsDriftTubeDetector*    dt = dynamic_cast<CsDriftTubeDetector*>(*Id);
      CsStrawTubesDetector*   st = dynamic_cast<CsStrawTubesDetector*>(*Id);
      prDets_.push_back( *Id );
      prDetsNames += (*Id)->GetTBName()+" ";
      if( dc ) cout << "<DC> ";
      if( dt ) cout << "<DT> ";
      if( st ) cout << "<ST> ";
      
      if( dc || st || dt || (*Id)->hasDrift() ) t0_.push_back( (*Id)->getT0() );
      else t0_.push_back( 0 );
      cout << " Done." << endl;
    }
  }
    
  // book the tree
  CsHistograms::SetCurrentPath( ROOT_PATH );
  cout << "InitCalibrationTree ... \n"; 
  string TreeName( TREE_NAME );
  T_eff_ = new TTree( TreeName.c_str(), TreeName.c_str() );
  T_eff_->Branch("T_evt",     &T_evt,       "T_evt/i",     BUFFER_SIZE);
  T_eff_->Branch("T_trigMsk", &T_trigMsk,   "T_trigMsk/i", BUFFER_SIZE);
  T_eff_->Branch("T_zone",    &T_zone,      "T_zone/i",    BUFFER_SIZE);
  T_eff_->Branch("T_cmlt",    &T_cmlt,      "T_cmlt/i",    BUFFER_SIZE);
  T_eff_->Branch("T_chi2",    &T_chi2,      "T_chi2/D",    BUFFER_SIZE);
  T_eff_->Branch("T_ndf",     &T_ndf,       "T_ndf/I",     BUFFER_SIZE);
  T_eff_->Branch("T_prob",    &T_prob,      "T_prob/D",    BUFFER_SIZE);
  T_eff_->Branch("T_cop",     &T_cop,       "T_cop/D",     BUFFER_SIZE);
  T_eff_->Branch("T_meanT",   &T_meanT,     "T_meanT/D",   BUFFER_SIZE);

  T_eff_->Branch("T_nDets",    &T_nDets,       "T_nDets/i", BUFFER_SIZE);
  T_eff_->Branch("T_detVect",  &T_detVect[0],  "T_detVect[T_nDets]/i",  BUFFER_SIZE);
  T_eff_->Branch("T_inActive", &T_inActive[0], "T_inActive[T_nDets]/i", BUFFER_SIZE);  
  T_eff_->Branch("T_fnd",      &T_fnd[0],      "T_fnd[T_nDets]/i",      BUFFER_SIZE);  

  // Track parameters extrapolated on each detector
  T_eff_->Branch("T_xVect",    &T_xVect[0],    "T_xVect[T_nDets]/D",   BUFFER_SIZE);  
  T_eff_->Branch("T_yVect",    &T_yVect[0],    "T_yVect[T_nDets]/D",   BUFFER_SIZE);  
  T_eff_->Branch("T_zVect",    &T_zVect[0],    "T_zVect[T_nDets]/D",   BUFFER_SIZE);  
  T_eff_->Branch("T_txVect",   &T_txVect[0],   "T_txVect[T_nDets]/D",  BUFFER_SIZE);  
  T_eff_->Branch("T_tyVect",   &T_tyVect[0],   "T_tyVect[T_nDets]/D",  BUFFER_SIZE);  
  T_eff_->Branch("T_dzVect",   &T_dzVect[0],   "T_dzVect[T_nDets]/D",  BUFFER_SIZE);  
  T_eff_->Branch("T_BVect",    &T_BVect[0],    "T_BVect[T_nDets]/D",   BUFFER_SIZE);  

  T_eff_->Branch("T_uLx",    &T_uLx[0],    "T_uLx[T_nDets]/D",   BUFFER_SIZE);  
  T_eff_->Branch("T_vLx",    &T_vLx[0],    "T_vLx[T_nDets]/D",   BUFFER_SIZE);  
  T_eff_->Branch("T_tuLx",   &T_tuLx[0],   "T_tuLx[T_nDets]/D",  BUFFER_SIZE);  
  T_eff_->Branch("T_tvLx",   &T_tvLx[0],   "T_tvLx[T_nDets]/D",  BUFFER_SIZE);  
    
  // Hits in detector
  T_eff_->Branch("T_duVect",   &T_duVect[0],   "T_duVect[T_nDets]/D",   BUFFER_SIZE);  
  T_eff_->Branch("T_dvVect",   &T_dvVect[0],   "T_dvVect[T_nDets]/D",   BUFFER_SIZE);  
  T_eff_->Branch("T_rVect",    &T_rVect[0],    "T_rVect[T_nDets]/D",   BUFFER_SIZE);  
  T_eff_->Branch("T_tVect",    &T_tVect[0],    "T_tVect[T_nDets]/D",   BUFFER_SIZE);  

  T_eff_->Branch("D_U",   &D_U[0],   "D_U[T_nDets]/D",   BUFFER_SIZE);  
  T_eff_->Branch("D_V",   &D_V[0],   "D_V[T_nDets]/D",   BUFFER_SIZE);  

  T_nDets = 0;
  firstCall = false;
  
  // write prDetsNames to the rootfile
  TObjString *T_ostr = new TObjString( prDetsNames.c_str() );
  T_ostr->Write( LIST_NAME );
  
  CsHistograms::SetCurrentPath( "/" );
  cout << "InitCalibrationTree. done.\n"; 
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

//_________________________________________________________________
bool FillCalibrationTree( unsigned int evt )
{
  
  T_evt = evt; 
  T_trigMsk = (unsigned int) CsEvent::Instance()->getTriggerMask();

  // Get Tracks - Check number of tracks
  list< CsTrack* > tr = CsEvent::Instance()->getTracks();
  if( !tr.size() ) return false;

  // Store only events with primary vertices
  if( require_vertex > 0 ) {
      list< CsVertex* > vt = CsEvent::Instance()->getVertices();
      if( !vt.size() ) return false;
      bool primary_found = false;
      for( list< CsVertex* >::iterator It = vt.begin(); It != vt.end(); It++ ) {
	  if((*It)->isPrimary())
	      if((*It)->getNTracks() >= require_vertex) {
		  primary_found = true;
		  break;
	      }
      }
      if(!primary_found) return false;
  }

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
  
    // Get List of zones associated to track
    list< CsZone* > zones = (*It)->getZones();
    T_zone = 0;
    for( list< CsZone* >::iterator Iz = zones.begin(); Iz != zones.end(); Iz++ ) T_zone += (1<<(*Iz)->getId());    
    
    // loop over desactivated detectors
    T_nDets = 0;
    for( unsigned int id = 0; id < prDets_.size(); id++ ) { 
      CsDetector *d = prDets_[id];
      HepMatrix iR = d->getRotWRSInv();
      double zD = d->getZcm();
      double wirP = d->getWirP();
    
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
            
      // Track position
      double x  = T_xVect[T_nDets]  = lxDet->getX();
      double y  = T_yVect[T_nDets]  = lxDet->getY();
      double z  = T_zVect[T_nDets]  = lxDet->getZ();
      double tx = T_txVect[T_nDets] = lxDet->getDXDZ();
      double ty = T_tyVect[T_nDets] = lxDet->getDYDZ();
      double uLx = iR(1,1)*( x + tx*( zD - z) ) + iR(1,2)*( y + ty*( zD - z ) );
      double vLx = iR(2,1)*( x + tx*( zD - z) ) + iR(2,2)*( y + ty*( zD - z ) );
      T_dzVect[T_nDets] = dzMin;

      T_uLx[T_nDets]  = iR(1,1)*T_xVect[T_nDets]  + iR(1,2)*T_yVect[T_nDets];
      T_vLx[T_nDets]  = iR(2,1)*T_xVect[T_nDets]  + iR(2,2)*T_yVect[T_nDets];
      T_tuLx[T_nDets] = iR(1,1)*T_txVect[T_nDets] + iR(1,2)*T_tyVect[T_nDets];
      T_tvLx[T_nDets] = iR(2,1)*T_txVect[T_nDets] + iR(2,2)*T_tyVect[T_nDets];
    
      // Magnetic field
      float bx, by, bz;
      CsGeom::Instance()->getCsField()->getField( x, y, z, bx, by, bz ); 
      T_BVect[T_nDets]= iR(2,1)*double(bx)+iR(2,2)*double(by);        
      
      // see of helix is in detector inActive Area
      T_inActive[T_nDets] = prDets_[id]->inActiveArea( T_xVect[T_nDets], T_yVect[T_nDets] );
    
      delete lxDet;
                     
      // Get Detectors List of clusters
      list< CsCluster* > cDet = d->getMyClusters();
      T_detVect[T_nDets] = d->GetID();
      T_fnd[T_nDets]   = 0;
      T_duVect[T_nDets] = 1e6;
      T_dvVect[T_nDets] = 1e6;
      T_rVect[T_nDets]  = 1e6;
      
      CsCluster* clFnd = 0;      
      bool firstCl = true;
      for( list<CsCluster*>::iterator Ic = cDet.begin(); Ic != cDet.end(); Ic++ ) {
        
        double uCl = (*Ic)->getU(); 
        double vCl = (*Ic)->getV(); 
        if( fabs( uCl - uLx ) < wirP ) T_fnd[T_nDets] = 1;

        // keep closest cluster ===
        if( firstCl || fabs(uCl - uLx) < fabs(T_duVect[T_nDets]) ){

	  D_U[T_nDets] = uCl;
	  D_V[T_nDets] = vCl;          
          T_duVect[T_nDets] = uCl - uLx;
          T_dvVect[T_nDets] = vCl - vLx;
          (*Ic)->getTime( T_tVect[T_nDets] );
          
          // Add detector t0 to tMin to have 'Absolute Times'
          T_tVect[T_nDets] += t0_[id];
          
          //=== Look for distance to wire
          int wire  = (*Ic)->getDigitsList().front()->getAddress();
          double uw = d->getWirD() 
            + iR(1,1) * d->getDeltaXCorrection()
            + iR(1,2) * d->getDeltaYCorrection() 
            + wire * wirP;
          T_rVect[T_nDets] = uCl-uw;
          
          firstCl = false;
          clFnd = (*Ic);
        } // test on cluster in track
      }   // loop over clusters in detector
      T_nDets++;
      if( T_nDets >= NPLAN ) {
        cout << "FillCalibrationTree - ERROR: too many detectors required. skipped.\n";
        break;
      }   // test on number of clusters
    }     // loop over detectors
    
    T_eff_->Fill();
    
  }       // loop over tracks

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
void RequestInteruption(int dummy)
{
	printf("\nRequestInteruption - INFO: Got signal (%i).\n",dummy);
	printf("RequestInteruption - INFO: Try waiting for end of event to finish.\n");
  InteruptionRequested = true;
	return;
}
