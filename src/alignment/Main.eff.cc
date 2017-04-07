// $Id: Main.eff.cc,v 1.14 2010/09/07 18:27:00 tnagel Exp $



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

#include "DetFileManager.h"

#include <signal.h>

using namespace std;
using namespace CLHEP;

// Calibration Tree is for Residuals, RTRelation, efficiency, ...
void InitCalibrationTree( void );               
bool FillCalibrationTree( unsigned int evt=0 );

void RequestInteruption(int dummy);
void ForceInteruption(int dummy);
bool InteruptionRequested;

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
    int nevt=0;
    int refresh; 
    if( ! CsOpt::Instance()->getOpt("events", "refresh rate", refresh) ) refresh = 1;
    bool doCalibration = CsOpt::Instance()->getOpt("main", "do calibration");

    // check for number of requested _good_ events (i.e. events with tracks)
    int maxgood;
    int ngood = 0;
    if( ! CsOpt::Instance()->getOpt("main", "good events", maxgood) ) maxgood = 0;
    if( maxgood ) cout << "main - INFO: Requested number of matching events: " << maxgood << ".\n";

    // Loop on events ===   
    while( (!InteruptionRequested) && event->getNextEvent() ) { 
      
      bool accept = true;       
      
      if( doCalibration ){ 
        if( nevt == 0 ) InitCalibrationTree(); 
        accept = FillCalibrationTree( nevt );
      }

      if( accept ) ngood++;
      if( !(CsRegistrySing::Instance()->callEoeMethods()) ) break;
      if((++nevt)%refresh == 0 ) cout << "Event: " << nevt << " (" << ngood << ")" << endl;

      if( maxgood && ngood > maxgood ) break;
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

#include "TTree.h"
#include "TMath.h"
#include "TObjString.h"
vector< TTree* > T_eff_;            // Efficicency tree for each off traffic detector
TTree* T_all_;            // Alignment tree for all detectors used by traffic

int nomomok = 0;

//============
// Calibration tree
//============
static unsigned int T_evt;     // Event number
static unsigned int T_trigMsk; // Event trigger mask
static unsigned int T_zone;    // Event zone mask
// static unsigned int T_cmlt;    // track multiplicity
static double T_chi2;          // track chisquare
static int    T_ndf;           // track number of degree of freedom
// static double T_prob;          // track chisquare propability
static double T_cop;           // charge over momentum
static double T_meanT;         // track time

static unsigned int T_nDets;           // number of detectors in the tree
static unsigned int T_detVect[NPLAN];  // detector ID
static unsigned int T_inActive[NPLAN]; // 1 if track is in detector active area, 0 otherwise
static unsigned int T_fnd[NPLAN];      // 1 if a close enough cluster was found in the detector, is in detector active area, 0 otherwise

static unsigned int T_nClust;      // total number of clusters

static double T_duVect[NCLUST*NPLAN];    // u_cluster - u_track
static double T_rVect[NCLUST*NPLAN];     // u_cluster - u_wire     
static double T_tVect[NCLUST*NPLAN];     // t_cluster 
static double T_idVect[NCLUST*NPLAN];    // plane id of cluster
static double T_sVect[NCLUST*NPLAN];     // cluster size

#define TREE_CLUSTER_A1
#ifdef TREE_CLUSTER_A1
// add some more cluster info
static double T_aVect[NCLUST*NPLAN];     // cluster analog value
#endif

//#define TREE_DIGITS
#ifdef TREE_DIGITS
static unsigned T_didVect[NCLUST*NPLAN]; // index of first digit from cluster in digit array

static unsigned int T_nDigit;                    // total number of digits
static double T_idDigit[NDIGIT*NCLUST*NPLAN];    // position of my cluster
static double T_tDigit[NDIGIT*NCLUST*NPLAN];     // digit time
static double T_chDigit[NDIGIT*NCLUST*NPLAN];    // digit channel
static double T_aDigit[NDIGIT*NCLUST*NPLAN];     // digit amplitude
#endif

//  branches specific to magnets_on_ = true
static double T_xVect[NPLAN];  // track X         
static double T_yVect[NPLAN];  // track Y         
static double T_zVect[NPLAN];  // track Z         
static double T_txVect[NPLAN]; // track dX/dZ     
static double T_tyVect[NPLAN]; // track dY/dZ     
static double T_siguVect[NPLAN]; // track sigma u at det position     

static double T_uLx[NPLAN];  // track U         
static double T_vLx[NPLAN];  // track V         
static double T_tuLx[NPLAN]; // track dU/dZ         
static double T_tvLx[NPLAN]; // track dV/dZ         

static double T_dzVect[NPLAN]; // Extrapolation length 
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
  if(gDirectory) {
    TDirectory *curdir = gDirectory;
    gDirectory->cd("/");
    DetFileManager df( CsInit::Instance()->getDetectorTable()->c_str() );
    df.Write();
    curdir->cd();
  }
  else {
    cerr<<"InitCalibration NO ROOT FILE OPENED !"<<endl;
    exit(1);
  }
  

  static bool firstCall = true;
  if( !firstCall ) return;
  
  // init detector vectors ===
  list<CsDetector*> dets = CsGeom::Instance()->getDetectors(); 
  prDets_.clear();
  t0_.clear();
  
  isOff_.clear();
  T_eff_.clear();

  CsOpt* opt = CsOpt::Instance();
  // momentum not requested ?
  opt->getOpt( "EFF","No_momentum_ok", nomomok);
    
  // Load Traffic desactivated detectors
  list< string > detOff; detOff.clear();
  opt->getOpt( "TraF","DetNameOff", detOff );

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

      if( dc ) cout << "<DC> ";
      if( dt ) cout << "<DT> ";
      if( st ) cout << "<ST> ";
      
      if( dc || st || dt || (*Id)->hasDrift() ) t0_.push_back( (*Id)->getT0() );
      else t0_.push_back( 0 );
      cout << " Done." << endl;
    }
  }
    
  // book the tree
  CsHistograms::SetCurrentPath("/CsEfficiency");
  cout << "InitCalibrationTree ... \n"; 
  string TreeName( "T_eff_" );
  T_all_ = new TTree( TreeName.c_str(), TreeName.c_str() );
  T_all_->Branch("T_evt",     &T_evt,       "T_evt/i",     BUFFER_SIZE);
  T_all_->Branch("T_trigMsk", &T_trigMsk,   "T_trigMsk/i", BUFFER_SIZE);
  T_all_->Branch("T_zone",    &T_zone,      "T_zone/i",    BUFFER_SIZE);
  //  T_all_->Branch("T_cmlt",    &T_cmlt,      "T_cmlt/i",    BUFFER_SIZE);
  T_all_->Branch("T_chi2",    &T_chi2,      "T_chi2/D",    BUFFER_SIZE);
  T_all_->Branch("T_ndf",     &T_ndf,       "T_ndf/I",     BUFFER_SIZE);
  //  T_all_->Branch("T_prob",    &T_prob,      "T_prob/D",    BUFFER_SIZE);
  T_all_->Branch("T_cop",     &T_cop,       "T_cop/D",     BUFFER_SIZE);
  T_all_->Branch("T_meanT",   &T_meanT,     "T_meanT/D",   BUFFER_SIZE);

  T_all_->Branch("T_nDets",    &T_nDets,       "T_nDets/i", BUFFER_SIZE);
  T_all_->Branch("T_detVect",  &T_detVect[0],  "T_detVect[T_nDets]/i",  BUFFER_SIZE);
  T_all_->Branch("T_inActive", &T_inActive[0], "T_inActive[T_nDets]/i", BUFFER_SIZE);  
  T_all_->Branch("T_fnd",      &T_fnd[0],      "T_fnd[T_nDets]/i",      BUFFER_SIZE);  

  // Track parameters extrapolated on each detector
  T_all_->Branch("T_xVect",    &T_xVect[0],    "T_xVect[T_nDets]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_yVect",    &T_yVect[0],    "T_yVect[T_nDets]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_zVect",    &T_zVect[0],    "T_zVect[T_nDets]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_txVect",   &T_txVect[0],   "T_txVect[T_nDets]/D",  BUFFER_SIZE);  
  T_all_->Branch("T_tyVect",   &T_tyVect[0],   "T_tyVect[T_nDets]/D",  BUFFER_SIZE);  
  T_all_->Branch("T_siguVect",   &T_siguVect[0],   "T_siguVect[T_nDets]/D",  BUFFER_SIZE);  


  T_all_->Branch("T_dzVect",   &T_dzVect[0],   "T_dzVect[T_nDets]/D",  BUFFER_SIZE);  
  T_all_->Branch("T_BVect",    &T_BVect[0],    "T_BVect[T_nDets]/D",   BUFFER_SIZE);  

  T_all_->Branch("T_uLx",    &T_uLx[0],    "T_uLx[T_nDets]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_vLx",    &T_vLx[0],    "T_vLx[T_nDets]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_tuLx",   &T_tuLx[0],   "T_tuLx[T_nDets]/D",  BUFFER_SIZE);  
  T_all_->Branch("T_tvLx",   &T_tvLx[0],   "T_tvLx[T_nDets]/D",  BUFFER_SIZE);  
  // number of clusters in detector
  T_all_->Branch("T_nClust", &T_nClust,       "T_nClust/i", BUFFER_SIZE);
  
  T_all_->Branch("T_duVect", &T_duVect[0],   "T_duVect[T_nClust]/D",   BUFFER_SIZE);  
  
  // This branches are filled with non 0 values for drift-like detectors only
  T_all_->Branch("T_rVect", &T_rVect[0],    "T_rVect[T_nClust]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_tVect", &T_tVect[0],    "T_tVect[T_nClust]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_idVect",    &T_idVect[0],    "T_tVect[T_nClust]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_sVect", &T_sVect[0],    "T_sVect[T_nClust]/D",   BUFFER_SIZE);  

#ifdef TREE_CLUSTER_A1
  T_all_->Branch("T_aVect", &T_aVect[0],    "T_aVect[T_nClust]/D",   BUFFER_SIZE);  
#endif
#ifdef TREE_DIGITS
  T_all_->Branch("T_didVect", &T_didVect[0],    "T_didVect[T_nClust]/i",   BUFFER_SIZE);  

  T_all_->Branch("T_nDigit", &T_nDigit,       "T_nDigit/i", BUFFER_SIZE);
  T_all_->Branch("T_chDigit", &T_chDigit[0],    "T_chDigit[T_nDigit]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_idDigit", &T_idDigit[0],    "T_idDigit[T_nDigit]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_tDigit", &T_tDigit[0],    "T_tDigit[T_nDigit]/D",   BUFFER_SIZE);  
  T_all_->Branch("T_aDigit",    &T_aDigit[0],    "T_aDigit[T_nDigit]/D",   BUFFER_SIZE); 
#endif
  T_nDets = 0;
  firstCall = false;
  
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
    T_chi2 = (*It)->getChi2();                     // track chi2
    
    int cmlt = (*It)->getClusters().size();          // track cluster size

    if( ! nomomok && fabs(T_cop) < 0.0001 ) {
      T_ndf = cmlt - 4; // not bridged track
      continue; // we keep only tracks with mom
    }
    else T_ndf = cmlt - 5;                       // bridged track    

    //T_prob = TMath::Prob( T_chi2, T_ndf );         // track chi2 probability
    T_meanT = ( (*It)->hasMeanTime() ) ? (*It)->getMeanTime():0; // track mean time
  
    // Get List of zones associated to track
    list< CsZone* > zones = (*It)->getZones();
    T_zone = 0;
    for( list< CsZone* >::iterator Iz = zones.begin(); Iz != zones.end(); Iz++ ) T_zone += (1<<(*Iz)->getId());    
    
    // loop over desactivated detectors
    T_nDets = 0;
    T_nClust = 0;

#ifdef TREE_DIGITS
    T_nDigit = 0;
#endif

    // flag set to true if the track is crossing the active zone 
    // of at least one of the dets being processed.
    // if this is not the case, the track will not enter the tree
    bool inactiveatleastone = false; 
    
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
      
      T_dzVect[T_nDets] = dzMin;
      T_uLx[T_nDets]  = iR(1,1)*T_xVect[T_nDets]  + iR(1,2)*T_yVect[T_nDets];
      T_vLx[T_nDets]  = iR(2,1)*T_xVect[T_nDets]  + iR(2,2)*T_yVect[T_nDets];
      T_tuLx[T_nDets] = iR(1,1)*T_txVect[T_nDets] + iR(1,2)*T_tyVect[T_nDets];
      T_tvLx[T_nDets] = iR(2,1)*T_txVect[T_nDets] + iR(2,2)*T_tyVect[T_nDets];
    
      
      // determination of sigma u : error on the helix U at the z position of
      // the detector :

      const double* cov = lxDet->getCov();
      
      // u  =  cost  sint * x
      // v    -sint  cost   y

      T_siguVect[T_nDets] = 
	cov[0] * iR(1,1) * iR(1,1) +
	+ 2*cov[1] * iR(1,1) * iR(1,2) +
	cov[2] * iR(1,2) * iR(1,2);   

      // Magnetic field
      float bx, by, bz;
      CsGeom::Instance()->getCsField()->getField( x, y, z, bx, by, bz ); 
      T_BVect[T_nDets]= iR(2,1)*double(bx)+iR(2,2)*double(by);        
      
      // see of helix is in detector inActive Area
      T_inActive[T_nDets] = prDets_[id]->inActiveArea( T_xVect[T_nDets], T_yVect[T_nDets] );
    
      if(!inactiveatleastone && T_inActive[T_nDets]) 
	inactiveatleastone = true;

      delete lxDet;
                     
      // Get Detectors List of clusters
      list< CsCluster* > cDet = d->getMyClusters();
      T_detVect[T_nDets] = d->GetID();
      T_fnd[T_nDets]   = 0;
      
      bool firstCl = true;
      // int closestpos = -1; // unused variable (jj)

      for( list<CsCluster*>::iterator Ic = cDet.begin(); Ic != cDet.end(); Ic++ ) {
        
        double uCl = (*Ic)->getU(); 
	
	// only fairly close clusters are kept
        if( fabs( uCl - uLx ) > wirP*DURANGE) continue;
	
        if( fabs( uCl - uLx ) < wirP ) T_fnd[T_nDets] = 1;
	
	
	//cluster time 
	T_duVect[T_nClust] = uCl - uLx;
	(*Ic)->getTime( T_tVect[T_nClust] );
	// Add detector t0 to tMin to have 'Absolute Times' ONLY Drift dets
	T_tVect[T_nClust] += t0_[id];

	//=== Look for distance to wire
	int wire  = (*Ic)->getDigitsList().front()->getAddress();
	double uw = d->getWirD() 
	  + iR(1,1) * d->getDeltaXCorrection()
	  + iR(1,2) * d->getDeltaYCorrection() 
	  + wire * wirP;
	T_rVect[T_nClust] = uCl-uw;
	
	// store detid
	T_idVect[T_nClust] = T_detVect[T_nDets];

	// cluster size 
	T_sVect[T_nClust] = (*Ic)->getDigitsList().size();
	
#ifdef TREE_CLUSTER_A1
	// cluster 1st analog value
	(*Ic)->getAnalog( T_aVect[T_nClust] );
#endif
        // looking for closest cluster 
        if( firstCl || 
	    fabs(uCl - uLx) < fabs(T_duVect[T_nClust]) ) {  
          firstCl = false;
	} 

#ifdef TREE_DIGITS
	T_didVect[T_nClust] = T_nDigit;
	int ndigit = 0;
	const list<CsDigit*>& digits = (*Ic)->getDigitsList();
	for(list<CsDigit*>::const_iterator id = digits.begin();
	    id != digits.end(); id++) {
	  
	  int dsiz = (*id)->getDataSize();
	  double* data = (*id)->getData();
	  switch (dsiz) {
	  case 2:
	    T_aDigit[ T_nDigit] = data[1];  // no break ...
	  case 1:
	    T_tDigit[ T_nDigit] = data[0];
	    T_chDigit[ T_nDigit] = (*id)->getAddress();	    
	    break;
	  case 0:
	    cerr<<"FillCalibrationTree - ERROR: empty digit."<<endl;
	    T_aDigit[ T_nDigit] = 0;
	    T_tDigit[ T_nDigit] = 0;
	    T_chDigit[ T_nDigit] = -1;	    
	    break;
	  default:
	    break;
	  }
	  T_idDigit[ T_nDigit] = T_nClust;
	  T_nDigit++;

	  if(T_nDigit >= NDIGIT*NCLUST*NPLAN) {
	    cout << "FillCalibrationTree - ERROR: too many digits required. skipped.\n";
	    break;
	  }
	}

#endif	
	T_nClust++;

	if( T_nClust >= NCLUST*NPLAN ) {
	  cout << "FillCalibrationTree - ERROR: too many clusters required. skipped.\n";
	  break;
	}   // test on number of clusters
      }   // loop over clusters in detector
      T_nDets++;
      
      
      if( T_nDets >= NPLAN ) {
        cout << "FillCalibrationTree - ERROR: too many detectors required. skipped.\n";
        break;
      }   // test on number of detectors
    }     // loop over detectors

    if(inactiveatleastone)
      T_all_->Fill();
    
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

