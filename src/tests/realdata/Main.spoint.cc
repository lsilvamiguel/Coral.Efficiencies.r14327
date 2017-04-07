// $Id: Main.spoint.cc,v 1.6 2010/09/07 18:27:02 tnagel Exp $

#include "TMath.h"
#include "Coral.h"
#include "CsOpt.h"
#include "CsGeom.h"
#include "CsStopwatch.h"
#include "CsRegistrySing.h"

#include "CsSPMaker.h"
#include "CsSPUtils.h"
#include "CsDetFamily.h"
#include "CsResOpt.h"
#include "CsSpacePoint.h"
#include "CsCalSpacePoint.h"
#include "CsDriftChamberDetector.h"
#include "CsDriftTubeDetector.h"
#include "CsStrawTubesDetector.h"

#include "Defs.h"
#include <signal.h>

using namespace std;
using namespace CLHEP;

//=== prototypes
void FillCalibrationTree( vector< CsCalSpacePoint* > spl, unsigned int evt = 0 );
void RequestInteruption(int dummy);
void ForceInteruption(int dummy);
bool InteruptionRequested;

bool do_calibration = true;   //!< if false, tree is not booked nor filled
double chi2_cut = -1;       //!< Cut on chi2

//______________________________________________________________________
int main( int argc, char *argv[] )
{
  try
  {
    InteruptionRequested = false;

  	// redirect interuptions 
  	signal(SIGINT ,ForceInteruption);
  	signal(SIGTERM,ForceInteruption);
  	signal(SIGABRT,RequestInteruption);
  	signal(SIGALRM,RequestInteruption);
        
    //=== Package Initialization === 
    Coral* coral        = Coral::init( argc, argv );
    list<CsZone*> zones = CsGeom::Instance()->getZones();         //get zones
    list<CsDetector*> dets = CsGeom::Instance()->getDetectors();  //get Detectors    
    vector<CsDetFamily*> dfV = CsSPMaker::Instance()->getDetFamilies();
 
    int nevt = 0;
    int refresh; 
    bool do_calibration = ( CsOpt::Instance()->getOpt("main", "do calibration") );
    if( !CsOpt::Instance()->getOpt("main", "chi2 cut", chi2_cut) ) chi2_cut = -1;
    if( !CsOpt::Instance()->getOpt("events", "refresh rate", refresh) ) refresh = 1;

    while( (!InteruptionRequested) && coral->getNextEvent() ) {
      nevt++;
      for( unsigned int i = 0; i < dfV.size(); i++ )
      if( do_calibration ) FillCalibrationTree( dfV[i]->getCalibrationSpacePoints( ) , nevt);
      if( nevt%refresh == 0 ) cout << "Event: " << nevt << endl;
      
    }		    // loop over events
    
    CsRegistrySing::Instance()->callEndMethods();
    CsErrLog::Instance()->dump( elDebugging );
  }
	
  catch( std::exception &e ) { cerr << "Exception (std::exception) was caught:\n" << e.what() << "\n";  }
  catch(const char *s) { cerr << "Exception (const char *) was caught:\n" << s << "\n"; }
  catch(...) { cerr << "Unknown exception was caught.\n"; }
  return 0;
}


#include <TTree.h>
#include <TObjString.h>
TTree* T_eff_ = 0;          // Efficicency tree for all detector families
vector< TTree* > T_align_;  // Alignment tree for each detector family

//============
// Calibration tree
//============

static unsigned int T_evt;          // Event number
static unsigned int T_trigMsk;      // Event trigger mask
static unsigned int T_cmlt;    // track multiplicity
static double T_chi2;               // track chisquare
static int T_ndf;                   // track number of degrees of freedom
static double T_prob;               // track chisquare propability
static unsigned int T_full;         // 1 if full minimisation was used

static unsigned int T_nDets;           // number of fired detectors
static unsigned int T_detVect[NPLAN];  // detector ID
static unsigned int T_inActive[NPLAN]; // 1 if track is in detector active area, 0 otherwise
static unsigned int T_fnd[NPLAN];      // 1 if a close enough cluster was found in the detector, is in detector active area, 0 otherwise

static double T_duVect[NPLAN]; // u_cluster - u_track
static double T_rVect[NPLAN];  // u_cluster - u_wire     
static double T_tVect[NPLAN];  // t_cluster cluster time 
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

//__________________________________________________
TTree *AddCalibrationTree( void )
{
  CsHistograms::SetCurrentPath( ROOT_PATH );
  cout << "AddCalibrationTree ...\n"; 
  
  string TreeName( TREE_NAME );
  TTree *tree = new TTree( TreeName.c_str(), TreeName.c_str() );
  tree->Branch("T_evt",     &T_evt,     "T_evt/i",     BUFFER_SIZE);
  tree->Branch("T_trigMsk", &T_trigMsk, "T_trigMsk/i", BUFFER_SIZE);
  tree->Branch("T_cmlt",    &T_cmlt,    "T_cmlt/i",    BUFFER_SIZE);
  tree->Branch("T_chi2",    &T_chi2,    "T_chi2/D",    BUFFER_SIZE);
  tree->Branch("T_ndf",     &T_ndf,     "T_ndf/I",     BUFFER_SIZE);
  tree->Branch("T_prob",    &T_prob,    "T_prob/D",    BUFFER_SIZE);
  tree->Branch("T_full",    &T_full,    "T_full/i",     BUFFER_SIZE);
  
  tree->Branch("T_nDets",    &T_nDets,       "T_nDets/i", BUFFER_SIZE);
  tree->Branch("T_detVect",  &T_detVect[0],  "T_detVect[T_nDets]/i",  BUFFER_SIZE);
  tree->Branch("T_inActive", &T_inActive[0], "T_inActive[T_nDets]/i", BUFFER_SIZE); 
  tree->Branch("T_fnd",      &T_fnd[0],      "T_fnd[T_nDets]/i",      BUFFER_SIZE); 


  // Track parameters extrapolated on each detector
  tree->Branch("T_xVect",    &T_xVect[0],    "T_xVect[T_nDets]/D",   BUFFER_SIZE);  
  tree->Branch("T_yVect",    &T_yVect[0],    "T_yVect[T_nDets]/D",   BUFFER_SIZE);  
  tree->Branch("T_zVect",    &T_zVect[0],    "T_zVect[T_nDets]/D",   BUFFER_SIZE);  
  tree->Branch("T_txVect",   &T_txVect[0],   "T_txVect[T_nDets]/D",  BUFFER_SIZE);  
  tree->Branch("T_tyVect",   &T_tyVect[0],   "T_tyVect[T_nDets]/D",  BUFFER_SIZE);  
  tree->Branch("T_dzVect",   &T_dzVect[0],   "T_dzVect[T_nDets]/D",  BUFFER_SIZE);  

  tree->Branch("T_uLx",    &T_uLx[0],    "T_uLx[T_nDets]/D",   BUFFER_SIZE);  
  tree->Branch("T_vLx",    &T_vLx[0],    "T_vLx[T_nDets]/D",   BUFFER_SIZE);  
  tree->Branch("T_tuLx",   &T_tuLx[0],   "T_tuLx[T_nDets]/D",  BUFFER_SIZE);  
  tree->Branch("T_tvLx",   &T_tvLx[0],   "T_tvLx[T_nDets]/D",  BUFFER_SIZE);  
    
  // Hits in detector
  tree->Branch("T_duVect",   &T_duVect[0],   "T_duVect[T_nDets]/D",   BUFFER_SIZE);  
  tree->Branch("T_rVect",    &T_rVect[0],    "T_rVect[T_nDets]/D",   BUFFER_SIZE);  
  tree->Branch("T_tVect",    &T_tVect[0],    "T_tVect[T_nDets]/D",   BUFFER_SIZE);  
  T_nDets = 0;
  
  T_eff_ = tree ;
  
  // write prDetsNames to the rootfile
  string prDetsNames = "";
  vector<CsDetFamily*> dfV = CsSPMaker::Instance()->getDetFamilies();
  for( unsigned int i=0; i < dfV.size(); i++ ) {
    vector< CsDetector* > dets = dfV[i]->getDetectors();
    for( unsigned int j=0; j < dets.size(); j++ ) 
    prDetsNames += dets[j]->GetTBName()+" ";
  }
  
  TObjString *T_ostr = new TObjString( prDetsNames.c_str() );
  T_ostr->Write( LIST_NAME );
  
  CsHistograms::SetCurrentPath( "/" );
  cout << "AddCalibrationTree. done.\n"; 
  
  return tree;
}

//__________________________________________________
TTree * GetCalibrationTree( CsDetector *d )
{ return ( T_eff_ )? T_eff_:AddCalibrationTree(); }

//__________________________________________________
void FillCalibrationTree( vector< CsCalSpacePoint* > spl, unsigned int evt )
{

  T_evt = evt; 
  T_trigMsk = (unsigned int) CsEvent::Instance()->getTriggerMask();
  for( unsigned int i=0; i < spl.size(); i++ ) {
    
    // Get off Detector
    CsDetector *d = spl[i]->detOff_;
    if( !d ){
      cout << "FillCalibrationTree - ERROR: calSpacepoint has no detOff.\n";
      continue;
    }
    
    // check minimisation method
    
    T_nDets = 0;
    T_detVect[T_nDets] = d->GetID();
    if( !( (T_full = spl[i]->minimised() ) || spl[i]->minimised_Fast() ) ) {
      cout << "FillCalibrationTree - ERROR: calSpacepoint is neither fast nor full minimized.\n";
      continue;
    }
    

    // retrieve detector usefull informations
    HepMatrix iR_ = d->getRotWRSInv();
    double wirP_ = d->getWirP(); 
    double wirD_cor_ = d->getWirD() 
       + iR_(1,1) * d->getDeltaXCorrection()
       + iR_(1,2) * d->getDeltaYCorrection();

    CsDriftChamberDetector* dc = dynamic_cast<CsDriftChamberDetector*>(d);
    CsDriftTubeDetector*    dt = dynamic_cast<CsDriftTubeDetector*>(d);
    CsStrawTubesDetector*   st = dynamic_cast<CsStrawTubesDetector*>(d);
    double t0_ = ( dc || st || dt || (d)->hasDrift() ) ? d->getT0():0;
   
    T_cmlt = spl[i]->cSize();
    spl[i]->getZ( T_zVect[T_nDets] );
    T_dzVect[T_nDets] = T_zVect[T_nDets]-d->getZcm();
    
    if( T_full ) {
      spl[i]->getChi2( T_chi2 ); 
      T_ndf = T_cmlt-4;
      
      // update chi2 as CsSpacePoints objects has chi2/ndf
      T_chi2*= T_ndf ;    
      
      spl[i]->getX(  T_xVect[T_nDets] );
      spl[i]->getY(  T_yVect[T_nDets] );
      spl[i]->getTx( T_txVect[T_nDets] );
      spl[i]->getTy( T_tyVect[T_nDets] );
      
    } else {
      spl[i]->getChi2_Fast( T_chi2 ); 
      T_ndf = T_cmlt-2;
      
      // update chi2 as CsSpacePoints objects has chi2/ndf
      T_chi2*= T_ndf ;
      
      spl[i]->getX_Fast( T_xVect[T_nDets] );
      spl[i]->getY_Fast( T_yVect[T_nDets] );
      T_txVect[T_nDets] = ( spl[i]->getMode() ) ? 0:T_xVect[T_nDets]/T_zVect[T_nDets];
      T_tyVect[T_nDets] = ( spl[i]->getMode() ) ? 0:T_yVect[T_nDets]/T_zVect[T_nDets];
    }

    // Check chi2/ndf
    if( chi2_cut > 0 && T_chi2/T_ndf > chi2_cut ) continue;
    
    T_prob = TMath::Prob( T_chi2, T_ndf );
    
    T_uLx[T_nDets]  = iR_(1,1)*T_xVect[T_nDets] + iR_(1,2)*T_yVect[T_nDets];
    T_vLx[T_nDets]  = iR_(2,1)*T_xVect[T_nDets] + iR_(2,2)*T_yVect[T_nDets];
    
    if( T_full ) {
      T_tuLx[T_nDets]  = iR_(1,1)*T_txVect[T_nDets] + iR_(1,2)*T_tyVect[T_nDets];
      T_tvLx[T_nDets]  = iR_(2,1)*T_txVect[T_nDets] + iR_(2,2)*T_tyVect[T_nDets];
    } else {
      T_tuLx[T_nDets] = ( spl[i]->getMode() ) ? 0:T_uLx[T_nDets]/T_zVect[T_nDets];
      T_tvLx[T_nDets] = ( spl[i]->getMode() ) ? 0:T_vLx[T_nDets]/T_zVect[T_nDets];
    }
    
    int wire = ( (T_uLx[T_nDets] - wirD_cor_)/wirP_ < 0 ) ? 
      int( (T_uLx[T_nDets] - wirD_cor_)/wirP_ -0.5 ):
      int( (T_uLx[T_nDets] - wirD_cor_)/wirP_ +0.5 );
    
    T_fnd[T_nDets]      = spl[i]->found_;
    T_inActive[T_nDets] = d->inActiveArea( T_xVect[T_nDets], T_yVect[T_nDets] );
    T_duVect[T_nDets] = 1e6;
    T_tVect[T_nDets]  = 1e6;
    T_rVect[T_nDets]  = 1e6;
    
    if( T_fnd[T_nDets] ) {
      double uCl = spl[i]->clFound_->getU(); 
      spl[i]->clFound_->getTime( T_tVect[T_nDets] );
      T_tVect[T_nDets] += t0_;
      T_duVect[T_nDets] = uCl - T_uLx[T_nDets];
    }
    
    // increment number of stored detectors (always 1!)
    T_nDets++;
    
    TTree *tree = 0;
    if( ( tree = GetCalibrationTree( d ) ) ) tree->Fill();
  }   
       
  return;
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
