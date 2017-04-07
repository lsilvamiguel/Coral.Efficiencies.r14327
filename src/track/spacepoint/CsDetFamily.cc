
// $Id: CsDetFamily.cc,v 1.43 2003/04/24 07:23:24 benigno Exp $

/*!
   \file    CsDetFamily.cc
   \brief   Plane list corresponding to one space point + associated cuts.
   \author  Hugo Pereira
   \version $Revision: 1.43 $
   \date    $Date: 2003/04/24 07:23:24 $
*/

#include "CsDetFamily.h"
#include "CsErrLog.h"
#include "CsGeom.h"
#include "CsDetector.h"
#include "CsCluster.h"
#include "CsEventUtils.h"
#include "CsSpacePoint.h"
#include "CsCalSpacePoint.h"
#include "CsOpt.h"
#include "CsDriftChamberDetector.h"
#include "CsStrawTubesDetector.h"
#include "CsRTRelation.h"

using namespace std;
using namespace CLHEP;

//_PUBLIC METHODS______________________________________________________________
CsDetFamily::CsDetFamily( const int id  ) : 
  id_( id ), name_( "NoName" ),
  z_( 0 ),
  nClCut_( -1 ), 
  chi2cut_Fast_( -1 ),
  chi2cut_( -1 ),
  mode_( STRAIGHT ),
  bigMask_( 0 ),
  geomOK_( false ),
  nCmb_( 0 ),
  hLevel_( None ),
  histoBooked_( false ),
  calHistoBooked_( false ),
  nTBins_( 50 )
{

  // clear detector information vectors
  d_.clear(); 
  dc_.clear();
  st_.clear();
  d_A_.clear();
  iR_.clear();
  wirP_.clear();
  wirD_cor_.clear();
  maxMlt_.clear();
  
  // clear cluster information vectors
  c_.clear();
  cd_.clear();
  uc_.clear();
  wc_.clear();
  rc2_.clear();

  // clear minimisation matrices
  A_Fast_.clear(); 
  B_Fast_.clear();   
  A_.clear(); 
  B_.clear();
  
  // clear associated cluster information vectors
  ic_A_.clear();
  c_A_.clear();
  cd_A_.clear();
  uc_A_.clear();
  wc_A_.clear();
  rc2_A_.clear();
  A_Fast_A_.clear();
  B_Fast_A_.clear();
  A_A_.clear();
  B_A_.clear();
  
  // clear indexes
  nc_.clear();
  fc_.clear();
  ic_.clear();

  // clear _bookCalHisto vectors 
  H_du_.clear();         // residual     (deltaU = u-uRef)
  H_r_vs_t_.clear();    // rt relation  (rRef vs t)
  H_ra_vs_t_.clear();   // rt relation  (fabs(rRef) vs t)
  H_du_vs_u_.clear();   // correlations (deltaU vs uRef)
  H_du_vs_v_.clear();   // correlations (deltaU vs vRef)
  H_du_vs_r_.clear();   // correlations (deltaU vs rRef)
  H_dt_vs_r_.clear();   // correlations (DeltaT vs rRef)
  H_wir_.clear();       // efficiency   (wiRef)
  H_wid_.clear();       // efficiency   (wiDet)
  H_chi2_cal_.clear();  // ref point chi2 
  H_spMlt_cal_.clear(); // spacepoint multiplicity
   
  // clear spacepoint lists
  sp_.clear();
  sp_cal_.clear();
  
  // Get settings
  useLRFilter_   = CsOpt::Instance()->getOpt("SP", "useLRFilter");  
  removeADet_    = CsOpt::Instance()->getOpt("SP", "removeADet");  
  useFastMinim_  = CsOpt::Instance()->getOpt("SP", "useFastMinimisation");  
  useFullMinim_  = CsOpt::Instance()->getOpt("SP", "useFullMinimisation"); 
  useClusterAsc_ = CsOpt::Instance()->getOpt("SP", "useClusterAssociation"); 
  tgCenter_ = CsGeom::Instance()->getTargetCenter();  

  // Get Histogram level
  string hLevStr("none");
  if( CsOpt::Instance()->getOpt("SP", "hist level", hLevStr ) ) {
    
    if( hLevStr == "none" )        hLevel_ = None;
    else if( hLevStr == "normal" ) hLevel_ = Normal;
    else if( hLevStr == "high" )   hLevel_ = High;
    else CsErrLog::msg(elError, __FILE__, __LINE__, 
      "Unrecognised histogram level: \"%s\". Default is \"none\".",
      hLevStr.c_str() );
      
  } else hLevel_ = None;

  double arg;
  if( CsOpt::Instance()->getOpt("RT","nTBins", arg ) ) nTBins_ = (unsigned int) arg;
    
  // consistancy checks
  if( !( useFastMinim_ || useFullMinim_ ) ) 
  CsErrLog::Instance()->mes(elFatal,"CsDetFamily::CsDetFamily. At least \"SP useFastMinimisation\" or \"SP useFastMinimisation\" must be set.");
  
  // Warning. This check is ABSOLUTELY necessary otherwise residuals have bias leading to artificialy good resolution.
  if( useClusterAsc_ && (!removeADet_) ) {
    CsErrLog::Instance()->mes(elError,"CsDetFamily::CsDetFamily. \"SP useClusterAssociation\" selected. \"removeADet\" forced to true.");
    removeADet_ = true;
  }
  
  return;  
}

//_____________________________________________________________________________
CsDetFamily& CsDetFamily::operator=( const CsDetFamily & df)
{
  if( this != &df ) *this = df;
  return *this;
}

//_____________________________________________________________________________
bool CsDetFamily::operator==( const CsDetFamily & df) const 
{
  return ( 
    id_ == df.id_ &&
    name_ == df.name_ ) ? true : false;
}

//_____________________________________________________________________________
void CsDetFamily::addDetector( CsDetector &det, int maxMlt )
{
  CsDriftChamberDetector* dc = dynamic_cast<CsDriftChamberDetector*> (&det);
  CsStrawTubesDetector*   st = dynamic_cast<CsStrawTubesDetector*> (&det);
  
  // add detector and related informations
  d_.push_back( &det );
  dc_.push_back( dc );
  st_.push_back( st );
  d_A_.push_back(-1);
  maxMlt_.push_back( maxMlt );
  iR_.push_back( det.getRotWRSInv() );
  wirP_.push_back( det.getWirP() );
  wirD_cor_.push_back( det.getWirD()          
    + det.getRotWRSInv()(1,1) * det.getDeltaXCorrection()
    + det.getRotWRSInv()(1,2) * det.getDeltaYCorrection() 
  );
  
  if( ! (removeADet_ || useClusterAsc_) ) return;
  
  // fill d_assoc
  unsigned int id = d_.size() - 1;
  if( dc_[id] != NULL && dc_[id]->hasAssociateDet() ) {
    string aName = dc_[id]->getAssociateDet()->GetTBName();
    for( unsigned int i=0; i< d_.size()-1; i++ ) 
    if( d_[i]->GetTBName() == aName ) {
      d_A_[i] = id;
      d_A_[id] = i;
      if( useClusterAsc_ ) bigMask_ |= (1 << id);
    }
  } else if( st_[id] != NULL && st_[id]->hasAssociateDet() ) {
    string aName = st_[id]->getAssociateDet()->GetTBName();
    for( unsigned int i=0; i< d_.size()-1; i++ ) 
    if( d_[i]->GetTBName() == aName ) {
      d_A_[i] = id;
      d_A_[id] = i;
      if( useClusterAsc_ ) bigMask_ |= (1 << id);
    }
  }
  
  return;
}
  
//_____________________________________________________________________________
vector <CsSpacePoint*> CsDetFamily::getSpacePoints( void )
{  
  if( hLevel_ >= Normal && !histoBooked_ ) _bookHistograms();
  
  if (! _importClusters() ) {
    CsErrLog::Instance()->mes(elWarning,"CsDetFamily::buildSpacePoints. Too few clusters in event. Skipped.");
    return sp_;
  }
  
  if(! _checkMlt(0) ) {
    CsErrLog::Instance()->mes(elWarning,"CsDetFamily::buildSpacePoints. Too large multiplicities. Event Skipped.");
    return sp_;
  }
  
  double chi2_min = -1;
  
  while( hasClusterCmb_ ) {
        
    double x_Fast, y_Fast, chi2_Fast;
    if( useFastMinim_ ) {
      
      if( !_minimise_Fast( x_Fast, y_Fast, chi2_Fast, 0 ) ) {
 	  	  CsErrLog::Instance()->mes(elWarning,"CsDetFamily::buildSpacePoints. Troubles with _minimise_Fast().");
         _getNextClusterCmb( 0 );
         continue;
      }
     
      // fill chi2_Fast distribution
      if( histoBooked_ ) H_chi2_Fast_->Fill( chi2_Fast );
     
      // reject too large chi2_fast
      if( histoBooked_ && chi2cut_Fast_ > 0 && chi2_Fast > chi2cut_Fast_ ) H_status_->Fill( 2 );
      if( chi2_Fast < 0 || ( chi2cut_Fast_ > 0 && chi2_Fast > chi2cut_Fast_ ) ) {
        _getNextClusterCmb( 0 );
        continue;
      }
      
      if( histoBooked_ )            H_status_->Fill( 0 );
      if( histoBooked_ && geomOK_ ) H_profile_Fast_->Fill( x_Fast, y_Fast );         
          
    }      // succesfull Fast minimisation
        
    double x, y, tx, ty, chi2;
    if( useFullMinim_ ) {
      if( !_minimise( x, y, tx, ty, chi2, 0 ) ) {
	  	  CsErrLog::Instance()->mes(elWarning,"CsDetFamily::buildSpacePoints. Troubles with _minimise().");
        _getNextClusterCmb( 0 );
        continue;
      }
            
      if( histoBooked_ ) H_chi2_->Fill( chi2 );
      
      if( chi2_min < 0 || ( chi2 > 0 && chi2 < chi2_min ) ) chi2_min = chi2;

      // reject too large chi2_fast
      if( histoBooked_ && chi2cut_ > 0 && chi2 > chi2cut_ ) H_status_->Fill( 3 );
      if( chi2 < 0 || ( chi2cut_ > 0 && chi2 > chi2cut_ ) ) {
        _getNextClusterCmb( 0 );
        continue;
      }
    
      if( histoBooked_ ) H_status_->Fill( 1 );
      if( histoBooked_ && geomOK_ ){
        H_profile_->Fill( x, y );
        H_angles_->Fill( tx, ty );
        H_tpx_->Fill( x, tx );
        H_tpy_->Fill( y, ty );
      }
      
    } // succesfull Full minimisation
  
    // build new space point
    list< CsCluster* > c; c.clear();
    
    for( unsigned int i=0; i<d_.size(); i++ )
    if( ic_[i] >= fc_[i] && 
      ic_[i] < fc_[i]+nc_[i] && 
      !(bigMask_ & (1 << i)) ) { 
      
      register unsigned int ic = ic_[i];
      c.push_back( c_[ic] );
      
      if( useClusterAsc_ && ic_A_[ic] >= 0 ) 
      c.push_back( c_A_[ (unsigned int) ic_A_[ic] ] );
 
    }
    
    CsSpacePoint* sp = new CsSpacePoint( *this, c, z_, (unsigned int) mode_ );
    if( useFastMinim_ ) sp->fill_Fast( x_Fast, y_Fast, chi2_Fast );
    if( useFullMinim_ ) sp->fill( x, y, tx, ty, chi2 );
    sp_.push_back( sp );
    
    if( histoBooked_ ) H_spSize_->Fill( c.size() );
    
    _getNextClusterCmb( 0 );
  }

  if( histoBooked_ ) H_spMlt_->Fill( sp_.size() );
  if( histoBooked_ && useFullMinim_ ) H_chi2_min_->Fill( chi2_min );
   
  return sp_;
}   
   
//_____________________________________________________________________________
vector< CsCalSpacePoint* > CsDetFamily::getCalibrationSpacePoints( void )
{
  static bool firstCall = true;
  if( firstCall ) { firstCall = false; cout << "CsDetFamily::getCalibrationSpacePoints - WARNING: removeADet is forced ON.\n"; }
  if( hLevel_ >= High && !calHistoBooked_ ) _bookCalHistograms();
  
  if (! _importClusters() ){
    CsErrLog::Instance()->mes(elWarning,"CsDetFamily::monitorEff. Too few clusters in event. Skipped.");
    return sp_cal_;
  }
  
  //=== Looop over detectors in the family ===  
  for( unsigned int i = 0; i<d_.size(); i++ ) {
    
    //=== Desactivate detector and (if needed) it's associate ===
    unsigned int mask = 0 | ( 1 << i );      
//    if( removeADet_ && d_A_[i] >= 0 ) mask |= ( 1 << d_A_[i] ); 
    if( d_A_[i] >= 0 ) mask |= ( 1 << d_A_[i] ); 
    
    //=== check multiplicities for active detectors ===
    if( ! _checkMlt( mask ) ) continue; 
    
    nCmb_ = 0;
    list< CsCluster* > cl = d_[i]->getMyClusters();
    list< CsCluster* >::iterator Icl;
    unsigned int calMlt = 0;
    
    //=== Loop over combinations ===
    while( hasClusterCmb_ ) {
      double x, y, tx, ty, chi2=-1;
      
      //=== perform minimisation
      bool result = ( useFullMinim_ ) ? _minimise( x, y, tx, ty, chi2, mask ) : _minimise_Fast( x, y, chi2, mask );
      
      //=== test the spacepoint
      if( !result ){
        _getNextClusterCmb( mask );
        continue;
      }
      
      //=== Check (and keep) chi2 (reduced)
      if( calHistoBooked_ ) H_chi2_cal_[i]->Fill( chi2 );
      if( chi2 < 0  || 
        ( (!useFullMinim_)  && chi2cut_Fast_ >= 0 && chi2 > chi2cut_Fast_ ) ||
        ( chi2cut_ >=0 && chi2 > chi2cut_ ) ) {
        _getNextClusterCmb( mask );
        continue;
      }
      
      //=== correct from angles ===
      //=== full minimization
      if( useFullMinim_ ){
        x+=(d_[i]->getZcm()-z_)*tx;
        y+=(d_[i]->getZcm()-z_)*ty;
     
      //=== fast minimization, target pointing mode
      } else if( mode_ == TPOINT ) {
        x *= (d_[i]->getZcm()-tgCenter_)/(z_-tgCenter_);
        y *= (d_[i]->getZcm()-tgCenter_)/(z_-tgCenter_);
      }
      
      //=== Get Space point (U,V) coordinates ===
      double uRef = iR_[i](1,1)*x + iR_[i](1,2)*y;
      double vRef = iR_[i](2,1)*x + iR_[i](2,2)*y;
      int wiRef = (int) ( (uRef-wirD_cor_[i])/wirP_[i] );

      // check if uRef is consistant with detector size
      if( wiRef < 0 || wiRef >= (int) d_[i]->getNWir() ) {
        _getNextClusterCmb( mask );
        continue;
      } else if( calHistoBooked_ ) H_wir_[i]->Fill( wiRef );
      
      //=== book new calibration spacepoint ===
      list< CsCluster* > clist; clist.clear();
      
      for( unsigned int icl=0; icl<d_.size(); icl++ )
      if( ic_[icl] >= fc_[icl] && 
        ic_[icl] < fc_[icl]+nc_[icl] && 
//        !(bigMask_ & (1 << icl)) ) { 
        !( (mask | bigMask_) & (1 << icl)) ) { 
        
        register unsigned int ic = ic_[icl];
        clist.push_back( c_[ic] );
        
        if( useClusterAsc_ && ic_A_[ic] >= 0 ) 
        clist.push_back( c_A_[ (unsigned int) ic_A_[ic] ] );
   
      }
      
      CsCalSpacePoint* sp = new CsCalSpacePoint( *this, clist, d_[i]->getZcm(), (unsigned int) mode_ );
      if( useFastMinim_ ) sp->fill_Fast( x, y, chi2 );
      if( useFullMinim_ ) sp->fill( x, y, tx, ty, chi2 );
      sp->detOff_ = d_[i]; 
      sp_cal_.push_back( sp );
      
      //=== scan detector i clusters ===
      CsCluster* c = NULL;
      Icl = cl.begin();
      double duMin = -1;
      while( Icl != cl.end() && ( duMin < 0 || fabs( uRef - (*Icl)->getU() ) < duMin ) ) {
        c = (*Icl);
        duMin = fabs( uRef - (*Icl)->getU() );
        Icl++;
      }
      
      if( c != NULL ) {
        
        calMlt++;
        double u = c->getU();
        double t = c->getDigitsList().front()->getDatum();
        int wire = (int) c->getDigitsList().front()->getAddress();
        double uW = wirD_cor_[i] + wirP_[i]*wire ;
        
        // fill histograms
        if( calHistoBooked_ ) {
          H_du_[i]->Fill( u-uRef );
          H_du_vs_r_[i]->Fill( uRef-uW, u-uRef );
        
          //=== rtRelations histograms for drift like detectors
          if( d_[i]->hasDrift() ) {
            H_r_vs_t_[i]->Fill(    t+d_[i]->getT0(), uRef-uW );
            H_ra_vs_t_[i]->Fill(   t+d_[i]->getT0(), fabs(uRef-uW) );
            H_r_vs_t_uniq_->Fill(  t+d_[i]->getT0(), uRef-uW );
            H_ra_vs_t_uniq_->Fill( t+d_[i]->getT0(), fabs(uRef-uW) );
          }
            
          if( geomOK_ ){
            H_du_vs_u_[i]->Fill( uRef, u-uRef ); 
            H_du_vs_v_[i]->Fill( vRef, u-uRef ); 
          } 
        }
        
        //=== Check distance between cluster and spacepoint,
        //=== forward to calibration space point, fill histo if any
        if( duMin > 0 && duMin < wirP_[i] ){
          sp->found_   = true;
          sp->clFound_ = c;
          if( calHistoBooked_ ) H_wid_[i]->Fill( wiRef );   
        }
       
        //=== For DC detectors, fill dt (estimated drift time difference) 
        if( calHistoBooked_ && dc_[i]!= NULL && dc_[i]->hasRTRelation() ) {
           CsRTRelation* rt = dc_[i]->getRTRelation();
           if( rt->hasGrid() ) { 
             bool error;         
             H_dt_vs_r_[i]->Fill( uRef-uW,  t-rt->getTfromR( fabs(uRef-uW), error ) ); 
           }
        }
        
        //=== For ST detectors, fill dt (estimated drift time difference) 
        if( calHistoBooked_ && st_[i]!= NULL && st_[i]->hasRTRelation() ) {
           CsRTRelation* rt = st_[i]->getRTRelation();
           if( rt->hasGrid() ) { 
             bool error;         
             H_dt_vs_r_[i]->Fill( uRef-uW,  t-rt->getTfromR( fabs(uRef-uW), error ) ); 
           }
        } 
        
        
      }   // test on best cluster found
      _getNextClusterCmb( mask );
    }    // loop over clusters combinations
    
    //=== Keep cluster multiplicity ===
    if( calHistoBooked_ ) H_spMlt_cal_[i]->Fill( calMlt );
    
    //=== reset cluster combinations ===
    for( unsigned int j = 0; j < d_.size(); j++ ) ic_[j] = fc_[j];
    hasClusterCmb_ = true;

  }   // loop over detectors
  return sp_cal_;
    
}      

//_____________________________________________________________________________
void CsDetFamily::cleanEvent( void )
{
  
  // clean spacepoints
  if( !sp_.empty() )
  for( unsigned int i = 0; i < sp_.size(); i++ ) delete sp_[i];
  sp_.clear();
  
  // clean calibration spacepoints
  if( !sp_cal_.empty() )
  for( unsigned int i = 0; i < sp_.size(); i++ ) delete sp_cal_[i];
  sp_cal_.clear();

  nCmb_ = 0;
  
  // clear cluster informations
  c_.clear();
  cd_.clear();
  uc_.clear();
  wc_.clear();
  rc2_.clear();
  A_Fast_.clear();
  B_Fast_.clear();
  A_.clear();
  B_.clear();
  
  if( useClusterAsc_ ) {
    ic_A_.clear();
    c_A_.clear();
    cd_A_.clear();
    uc_A_.clear();
    wc_A_.clear();
    rc2_A_.clear();
    A_Fast_A_.clear();
    B_Fast_A_.clear();
    A_A_.clear();
    B_A_.clear();
  }
  
  // clear cluster combination indexes
  fc_.clear();
  ic_.clear();
  nc_.clear();
  hasClusterCmb_ = true;
  
  return;
}  

//___DEBUGING__________________________________________________________________
void CsDetFamily::dumpConfig( void ) 
{

  cout << "Family#"<< id_ << " (" << name_.c_str() << "): " << d_.size() << " detectors" << endl;  
  
  // detector informations
  for( unsigned int i = 0; i < d_.size(); i++ ) {
    if( useClusterAsc_ && (bigMask_ & ( 1 << i )) ) continue;
    cout << " " << d_[i]->GetTBName().c_str() 
      << " (MxMlt=" << maxMlt_[i] << ")"; 
    if( d_A_[i] >= 0 ) 
    cout << " <-> " << d_[d_A_[i]]->GetTBName().c_str();
    cout << endl;
  }
  
  cout << " zrec                  :"<< z_            << endl;
  cout << " chi2Cut_Fast          :"<< chi2cut_Fast_ << endl;
  cout << " chi2Cut               :"<< chi2cut_      << endl;
  cout << " NClCut                :"<< nClCut_       << endl;
  cout << " mode                  :"<< (( mode_== STRAIGHT ) ? "STRAIGHT":"TPOINT") << endl;
  cout << " target center         :"<< tgCenter_ << "mm\n";
  
  cout << " Geometry              :";
  if( geomOK_ ) dumpGeometry();
  else cout << "none" << endl;
  
  cout << " useLRFilter           :" <<( ( useLRFilter_ )   ? "yes":"no") << endl;
  cout << " removeADet            :" <<( ( removeADet_ )    ? "yes":"no") << endl;
  cout << " useClusterAssociation :" <<( ( useClusterAsc_ ) ? "yes":"no") << endl;
  cout << " useFastMinimisation   :" <<( ( useFastMinim_ )  ? "yes":"no") << endl;
  cout << " useFullMinimisation   :" <<( ( useFullMinim_ )  ? "yes":"no") << endl;
  cout << " Histogram Level       :" <<  hLevel_ << endl;
  cout << endl;  
  
  return;
}

//_____________________________________________________________________________
void CsDetFamily::dumpGeometry( void ) 
{
  if( geomOK_ ) 
  cout << "(" << xMin_ << "," << xMax_ << ")x(" 
    << yMin_ << "," << yMax_ << ")"
    << endl;
  return;
}

//_____________________________________________________________________________
void CsDetFamily::dumpAssociationTable( void ) 
{
  cout << "CsDetFamily::dumpAssociationTable." << endl;
  for( unsigned int i=0; i < d_A_.size(); i++ )
  cout << "  " << i << "<->" << d_A_[i] << endl;

  return;
}

//_____________________________________________________________________________
void CsDetFamily::dumpClusters( void )
{
  cout << "CsDetFamily::dumpClusters. Family#" << name_.c_str() << endl;
  for( unsigned int i = 0; i < c_.size(); i++ )  
  cout << "  " << c_[i]->getDetsList().front()->GetTBName() 
    << ": u=" << c_[i]->getU() << "(" << uc_[i] << ")"
    << "  w=" << c_[i]->getW() << "(" << wc_[i] << ")"    
    << "  p=" << c_[i]->getLRProb() 
    << "  cov=" << rc2_[i]
    << endl;

  return;
}

//_____________________________________________________________________________
void CsDetFamily::dumpClusterSizes( void )
{
  cout << "Family#" << id_ << "::dumpClusterSizes" << endl;
  for( unsigned int i = 0; i < d_.size(); i++ )  
  cout << "  " << d_[i]->GetTBName().c_str() << " has " << nc_[i]
    << " cls." << endl;   
}    
      
//______________________________________________________________________________
void CsDetFamily::dumpCmbIndexes( void )
{
  if( hasClusterCmb_ ) {
    cout << "cmb : ";
    for( unsigned int i = 0; i < d_.size(); i++ )
    cout << "[" << ic_[i]-fc_[i] << "," << ic_[i] << "] ";
    cout << endl;
  }
  else return;
}
      
//______________________________________________________________________________
void CsDetFamily::dumpNCmb( void )
{
  if( hasClusterCmb_ ) return;
  cout << "CsDetFamily::dumpNCmb. " << nCmb_ << endl;
  return;
}

//___PRIVATE METHODS___________________________________________________________

void CsDetFamily::_bookHistograms( void )
{
  cout << "CsDetFamily::bookHistograms. Booking Space Points histograms." << endl;
  CsHistograms::SetCurrentPath("/CsSPMonitor");
  char name[100];
  
  sprintf( name, "%s_spMlt_", getName().c_str() );
  H_spMlt_ = new CsHist1S( name, name, 10, 0, 10 );
  
  sprintf( name, "%s_spSize_", getName().c_str() );
  H_spSize_ = new CsHist1S( name, name, d_.size()+1, 0, d_.size()+1 );
  
  sprintf( name, "%s_status_", getName().c_str() );
  H_status_ = new CsHist1S( name, name, 10, 0, 10 );
   
  if( useFastMinim_ ) {
    sprintf( name, "%s_chi2_Fast_", getName().c_str() );
    H_chi2_Fast_ = new CsHist1S( name, name, 60, 0, 30 );
  
    if( geomOK_ ) {
      sprintf( name, "%s_profile_Fast_", getName().c_str() );
      H_profile_Fast_ = new CsHist2S( 
        name, name, 
        (int) fabs( xMax_ - xMin_ ), xMin_, xMax_, 
        (int) fabs( yMax_ - yMin_ ), yMin_, yMax_          
      );
    } // geom histograns
  }   // fast minimisation histograms
  
  if( useFullMinim_ ) {  
    sprintf( name, "%s_chi2_", getName().c_str() );
    H_chi2_ = new CsHist1S( name, name, 100, 0, 50 );
  
    sprintf( name, "%s_chi2_min_", getName().c_str() );
    H_chi2_min_ = new CsHist1S( name, name, 100, 0, 50 );
  
    if( geomOK_ ) {
  
      sprintf( name, "%s_profile_", getName().c_str() );
      H_profile_ = new CsHist2S( 
        name, name, 
        (int) fabs( xMax_ - xMin_ )/10, xMin_, xMax_, 
        (int) fabs( yMax_ - yMin_ )/10, yMin_, yMax_          
      );
    
      sprintf( name, "%s_angles_", getName().c_str() );
      H_angles_ = new CsHist2S( 
        name, name, 
        (int) fabs( xMax_ - xMin_ )/10, 2*xMin_/z_, 2*xMax_/z_, 
        (int) fabs( yMax_ - yMin_ )/10, 2*yMin_/z_, 2*yMax_/z_          
      );
     
      sprintf( name, "%s_tpx_", getName().c_str() );
      H_tpx_ = new CsHist2S( 
        name, name, 
        (int) fabs( xMax_ - xMin_ )/10, xMin_, xMax_, 
        (int) fabs( xMax_ - xMin_ )/10, 2*xMin_/z_, 2*xMax_/z_          
      );
     
      sprintf( name, "%s_tpy_", getName().c_str() );
      H_tpy_ = new CsHist2S( 
        name, name, 
        (int) fabs( yMax_ - yMin_ )/10, yMin_, yMax_, 
        (int) fabs( yMax_ - yMin_ )/10, 2*yMin_/z_, 2*yMax_/z_          
      );
    } // geom histograms
  }   // full minimisation histograms
  
  histoBooked_ = true;
}  

//_____________________________________________________________________________
void CsDetFamily::_bookCalHistograms( void )
{
  if( calHistoBooked_ ) return;
  cout << "CsDetFamily::_bookCalHistograms. Booking RTcalibration histograms." << endl;
  
  // check if only one Grid is used all drift chambers
  double tGateMax = 0;
  double wirPMax = 0;
  
  double t0Min = 0;
  double t0Max = 0;
  cout << "CsDetFamily::_bookCalHistograms - INFO: Absolute drift times are used.\n";
     
  for( unsigned int i=0; i < d_.size(); i++ ) {
    char name[100];
    
    //=== RT Relations ===
    //=== Filled only for detectors with drift time informations ===
    CsHistograms::SetCurrentPath("/CsRTRelations");    
    
    if( d_[i]->hasDrift() ) {
      // r_vs_t (RT relation - algebric)
      sprintf( name, "%s_r_vs_t_", d_[i]->GetTBName().c_str() );
      CsHist2S* h2s = new CsHist2S( 
        name, name, 
        nTBins_, d_[i]->getT0(), d_[i]->getTGate()+d_[i]->getT0(),
        100, -0.6*wirP_[i], 0.6*wirP_[i]
      );
      H_r_vs_t_.push_back( h2s );
     
      // ra_vs_t_ (RT relation - positive)
      sprintf( name, "%s_ra_vs_t_", d_[i]->GetTBName().c_str() );
      h2s = new CsHist2S( 
        name, name, 
        nTBins_, d_[i]->getT0(), d_[i]->getTGate()+d_[i]->getT0(),
        50, 0, 0.6*wirP_[i]
      );    
      H_ra_vs_t_.push_back( h2s );
    } else {
      H_r_vs_t_.push_back( 0 );
      H_ra_vs_t_.push_back( 0 );
    } 
    //=== Residuals ===
    CsHistograms::SetCurrentPath("/CsResiduals");    
    
    // du (residuals)
    sprintf( name, "%s_du_", d_[i]->GetTBName().c_str() );
    CsHist1S* hs = new CsHist1S( name, name, 70, -0.2*d_[i]->getWirP(), 0.2*d_[i]->getWirP() );  
    H_du_.push_back( hs );
  
    // du_vs_r (residual vs dist to wire)
    sprintf( name, "%s_du_vs_r_", d_[i]->GetTBName().c_str() );
    CsHist2S* h2s = new CsHist2S( 
      name, name, 
      200, -0.5*d_[i]->getWirP(), 0.5*d_[i]->getWirP(),
      200, -0.5*d_[i]->getWirP(), 0.5*d_[i]->getWirP()         
    );
    H_du_vs_r_.push_back( h2s );

    if( geomOK_ ) { 
      // du_vs_u_ (residual vs uref)
      sprintf( name, "%s_du_vs_u_", d_[i]->GetTBName().c_str() );
      h2s = new CsHist2S( 
        name, name, 
        int( fabs( xMax_ - xMin_ )/d_[i]->getWirP() ), xMin_, xMax_, 
        200, -0.5*d_[i]->getWirP(), 0.5*d_[i]->getWirP()         
      );
      H_du_vs_u_.push_back( h2s );

      // du_vs_v_ (residual vs vref, along the wire)
      sprintf( name, "%s_du_vs_v_", d_[i]->GetTBName().c_str() );
      h2s = new CsHist2S( 
        name, name, 
        int( fabs( xMax_ - xMin_ )/d_[i]->getWirP() ), xMin_, xMax_, 
        200, -0.5*d_[i]->getWirP(), 0.5*d_[i]->getWirP()         
      );
      H_du_vs_v_.push_back( h2s );
    }
  
    // dt_vs_r (time difference vs dist to wire)
    sprintf( name, "%s_dt_vs_r_", d_[i]->GetTBName().c_str() );
    h2s = new CsHist2S( 
      name, name, 
      200, -0.5*d_[i]->getWirP(), 0.5*d_[i]->getWirP(),
      int(d_[i]->getTGate()), -0.5*d_[i]->getTGate(), 0.5*d_[i]->getTGate()         
    );
    H_dt_vs_r_.push_back( h2s );
    
    //=== efficiency histograms ===
    CsHistograms::SetCurrentPath("/CsEfficiency");
    
    // wiRef (wire profile from space point)
    sprintf( name, "%s_wiRef_", d_[i]->GetTBName().c_str() );
    CsHist1D* hd = new CsHist1D( name, name, d_[i]->getNWir(), 0, d_[i]->getNWir() );  
    H_wir_.push_back( hd );
    
    // wiDet (wire profile from detector)
    sprintf( name, "%s_wiDet_", d_[i]->GetTBName().c_str() );
    hd = new CsHist1D( name, name, d_[i]->getNWir(), 0, d_[i]->getNWir() );  
    H_wid_.push_back( hd );    
    
    //=== other space points related histograms ===
    CsHistograms::SetCurrentPath("/CsSPCalibration");
  
    // spMlt (space point multiplicity)
    sprintf( name, "%s_%s_spCalMlt_", getName().c_str(), d_[i]->GetTBName().c_str() );
    hs = new CsHist1S( name, name, 10, 0, 10 );
    H_spMlt_cal_.push_back( hs );
       
    // chi2 (chi2 distribution)
    sprintf( name, "%s_%s_chi2_cal_", getName().c_str(), d_[i]->GetTBName().c_str() );
    hs = new CsHist1S( name, name, 100, 0, 25 );
    H_chi2_cal_.push_back( hs );

    if( i == 0 || tGateMax < d_[i]->getTGate() ) tGateMax = d_[i]->getTGate();
    if( i == 0 || wirPMax < wirP_[i] )                            wirPMax  = wirP_[i];    
    if( i == 0 || ( d_[i]->hasDrift() && t0Min > d_[i]->getT0() ) ) t0Min = d_[i]->getT0();
    if( i == 0 || ( d_[i]->hasDrift() && t0Max < d_[i]->getT0() ) ) t0Max = d_[i]->getT0();
    
  }
  
  //=== unique RT Relations  
  CsHistograms::SetCurrentPath("/CsRTRelations");    
  
  char name[100];
  sprintf( name, "%s_r_vs_t_Unique_", getName().c_str() );
  
  H_r_vs_t_uniq_ = new CsHist2S( 
    name, name, 
    nTBins_, t0Min, tGateMax+t0Max,
    100, -0.6*wirPMax, 0.6*wirPMax
  );

  sprintf( name, "%s_ra_vs_t_Unique_", getName().c_str() );
  H_ra_vs_t_uniq_ = new CsHist2S( 
    name, name, 
    nTBins_, t0Min, tGateMax+t0Max,
    50, 0, 0.6*wirPMax
  );
    
  calHistoBooked_ = true;
  
}

//_____________________________________________________________________________
bool CsDetFamily::_importClusters( void )
{

  // clean previous event
  cleanEvent();
  
  // import clusters then fill ic_, nc_ and fc_.
  for( unsigned int i = 0; i < d_.size(); i++ ) {
    if( useClusterAsc_ && (bigMask_ & ( 1<<i )) ) {
      nc_.push_back( 0 );
      fc_.push_back( c_.size() );
      ic_.push_back( c_.size() );
      continue;
    }
    
    list<CsCluster*> c = d_[i]->getMyClusters();    
    if( useLRFilter_ ) CsEventUtils::useLRFilter( c );
    fc_.push_back( c_.size() );
    ic_.push_back( c_.size() );
    nc_.push_back( c.size() );
     
    list<CsCluster*>::iterator Ic;
    for( Ic = c.begin(); Ic != c.end(); Ic++ ) {
      c_.push_back( *Ic );
      cd_.push_back( i );
      uc_.push_back( (*Ic)->getU() );
      wc_.push_back( (*Ic)->getW() );
      HepMatrix cov = (*Ic)->getCov(); 
	 	  rc2_.push_back( cov(1,1) );
      
      // book associated cluster informations in separate lists
      if( !useClusterAsc_ ) continue;
      if( ! ( (*Ic)->hasAssociateClusters() && d_A_[i] >= 0 ) ) ic_A_.push_back( -1 );
      else {
        CsCluster* ac = (*Ic)->getAssociateClusters().front();
        c_A_.push_back( ac );
        cd_A_.push_back( (unsigned int) d_A_[i] );
        uc_A_.push_back( ac->getU() );
        wc_A_.push_back( ac->getW() );
        cov = ac->getCov();
        rc2_A_.push_back( cov(1,1) );
        ic_A_.push_back( c_A_.size()-1 );
      }
            
    }
  }
    
  // check number of clusters
  if( nClCut_ > 0 && c_.size() < (unsigned int) nClCut_ ) return false;
  if( useFastMinim_ ) _fillMatrix_Fast();
  if( useFullMinim_ ) _fillMatrix(); 
  return true;
  
}
      
//_____________________________________________________________________________
void CsDetFamily::_getNextClusterCmb( unsigned int mask )
{
  nCmb_++;
  if( !hasClusterCmb_ ) return;
  
  for( unsigned int i = 0; i <= d_.size(); i++ ) 
  if( i==d_.size() ) hasClusterCmb_ = false;  // All combinations scanned 
  else { 
   if( nc_[i] == 0 || (( mask | bigMask_ ) & ( 1 << i )) ) continue;
   ic_[i]++;
    if( ic_[i] < nc_[i]+fc_[i] ) break;
    else ic_[i] = fc_[i];
  }
  return;
  
}  

//_____________________________________________________________________________
unsigned int CsDetFamily::_getCmbSize( unsigned int mask )
{
  if( !hasClusterCmb_ ) return 0;
	
  unsigned int cs = 0;
  for( unsigned int i = 0; i<d_.size(); i++ ) {
    register unsigned int ic = ic_[i];
    if( ic < fc_[i] || ic >= fc_[i]+nc_[i] ||
      (( mask | bigMask_ ) & ( 1 << i )) ) continue;
    cs++;    
    if( useClusterAsc_ && ic_A_[ic] >= 0 ) cs++;
  }
  
  return cs;
}

//_____________________________________________________________________________
bool CsDetFamily::_checkMlt( unsigned int mask )
{
  bool accept = true;
  for( unsigned int i = 0; i < d_.size() && accept; i++ ) {
    if( ( mask | bigMask_ ) & ( 1 << i ) ) continue;
    if( maxMlt_[i] > 0 && nc_[i] > (unsigned int) maxMlt_[i] )
    accept = false;
  }
  
  return accept;
}

//_______________________________________________________________________________
//___________________FAST minimisation METHODS (position only)___________________
//_______________________________________________________________________________

bool CsDetFamily::_minimise_Fast( double &x, double &y, double &chi2, unsigned int mask )
{
	// initialize minimisation matrices
	HepMatrix A(2,2,0);
	HepMatrix B(2,1,0);	
  
  int nc = 0; 
  if( !hasClusterCmb_ ) return false;
  
	// loop over clusters
	for( unsigned int i = 0; i<d_.size(); i++ ) {
    if( ic_[i] < fc_[i] || ic_[i] >= fc_[i]+nc_[i] ||
      (( mask | bigMask_ ) & ( 1 << i )) ) continue;
    register unsigned int ic = ic_[i];
    A += A_Fast_[ic];
    B += B_Fast_[ic];
    nc++;
    
    if( useClusterAsc_ && ic_A_[ic] >= 0 ) {
      register unsigned int ica = (unsigned int) ic_A_[ic];
      A += A_Fast_A_[ica]; 
      B += B_Fast_A_[ica]; 
      nc++;
    }
    
  }
  
  // check number of used clusters
  if( nc < 2 || ( nClCut_ > 0 && nc < nClCut_ ) ) {
		if( histoBooked_ ) H_status_->Fill( 4 );
    CsErrLog::Instance()->mes(elWarning,"Not enough clusters for fast minim.");
    return false;
  }
  
  // fill positions
  x = (A(2,2)*B(1,1) - A(1,2)*B(2,1))/(A(1,1)*A(2,2)-A(1,2)*A(1,2));
  y = (A(1,1)*B(2,1) - A(1,2)*B(1,1))/(A(1,1)*A(2,2)-A(1,2)*A(1,2));
  
  // fill Chi2
  chi2 = 0;
	for( unsigned int i = 0; i<d_.size(); i++ ) {
    if( ic_[i] < fc_[i] || ic_[i] >= fc_[i]+nc_[i] || 
      (( mask | bigMask_ ) & ( 1 << i )) ) continue;
    register unsigned int ic = ic_[i];
    double a = ( mode_ == STRAIGHT ) ? 1 : ( wc_[ic]-tgCenter_ )/( z_-tgCenter_ );
    chi2 += pow( iR_[i](1,1)*x*a + 
      iR_[i](1,2)*y*a +
      iR_[i](1,3)*wc_[ic] - uc_[ic] ,2 ) / rc2_[i];

    if( useClusterAsc_ && ic_A_[ic] >= 0 ) {
      register unsigned int ica = (unsigned int) ic_A_[ic];
      double aa = ( mode_ == STRAIGHT ) ? 1 : ( wc_A_[ica]-tgCenter_ ) / ( z_-tgCenter_ );
      chi2 += pow(iR_[ cd_A_[ica] ](1,1)*x*aa + 
        iR_[ cd_A_[ica] ](1,2)*y*aa +
        iR_[ cd_A_[ica] ](1,3)*wc_A_[ica] - uc_A_[ica] ,2 ) / rc2_A_[ica];      
    } 

  }
  
  int ndf = nc-2;
  if( ndf > 0 ) chi2 /= ndf;
  else {
		if( histoBooked_ ) H_status_->Fill( 6 );
    chi2 = -1;
  }
  
  return true;    
}

//_____________________________________________________________________________
void CsDetFamily::_fillMatrix_Fast( void )
{
  for( unsigned int i = 0; i < c_.size(); i++ ) {
    register unsigned int j = cd_[i];
    A_Fast_.push_back( _getA_Fast( iR_[j], uc_[i], wc_[i], rc2_[i]) );
    B_Fast_.push_back( _getB_Fast( iR_[j], uc_[i], wc_[i], rc2_[i]) );
  }
  
  // fill matrix for associated clusters
  if( useClusterAsc_ )
  for( unsigned int i = 0; i< c_A_.size(); i++ ) {
    register unsigned int j = cd_A_[i];
    A_Fast_A_.push_back( _getA_Fast( iR_[j], uc_A_[i], wc_A_[i], rc2_A_[i]) );
    B_Fast_A_.push_back( _getB_Fast( iR_[j], uc_A_[i], wc_A_[i], rc2_A_[i]) );
  }
      
} 

//_______________________________________________________________________________
HepMatrix CsDetFamily::_getA_Fast( HepMatrix iR, double uc, double wc, double rc2 )
{
  HepMatrix A(2,2,0);
  double a = ( mode_ == STRAIGHT ) ? 1 : ( wc-tgCenter_ )/( z_-tgCenter_ );

  A(1,1) = pow( iR(1,1)*a,2 ) / rc2;
  A(2,2) = pow( iR(1,2)*a,2 ) / rc2;
  A(1,2) = iR(1,1)*iR(1,2)*a*a / rc2;
  return A;
}
  
//_________________________________________________________________________________
HepMatrix CsDetFamily::_getB_Fast( HepMatrix iR, double uc, double wc, double rc2 )
{
  HepMatrix B(2,1,0);	
  double a = ( mode_ == STRAIGHT ) ? 1:( wc-tgCenter_ )/( z_-tgCenter_ );

  B(1,1) = iR(1,1)*(uc-wc*iR(1,3))*a/rc2;
	B(2,1) = iR(1,2)*(uc-wc*iR(1,3))*a/rc2;
  return B;
}  


//_______________________________________________________________________________
//___________________FULL minimisation METHOD (position + angles)________________
//_______________________________________________________________________________

bool CsDetFamily::_minimise( double &x, double &y, double &tx, double &ty, double &chi2, unsigned int mask )
{
	// initialize minimisation matrices
	HepMatrix A(4,4,0);
	HepMatrix B(4,1,0);	
  
  int nc = 0; 
  if( !hasClusterCmb_ ) return false;
  
	// loop over clusters
	for( unsigned int i = 0; i<d_.size(); i++ ) {
    if( ic_[i] < fc_[i] || ic_[i] >= fc_[i]+nc_[i] || 
      (( mask | bigMask_ ) & ( 1 << i )) ) continue;
    register unsigned int ic = ic_[i];
    A += A_[ic];
    B += B_[ic];
    nc++;
    
    if( useClusterAsc_ && ic_A_[ic] >= 0 ) {
      register unsigned int ica = (unsigned int) ic_A_[ic];
      A += A_A_[ica]; 
      B += B_A_[ica]; 
      nc++;
    }

  }
  
  if( nc < 4 || ( nClCut_ > 0 && nc < nClCut_ ) ) {
    if( histoBooked_ )  H_status_->Fill( 5 );
		CsErrLog::Instance()->mes(elWarning,"Not enough clusters for full minim.");
    return false;
  }
  
  // inverse A
  int err;
	A.invert( err );
	if( err ) {
    if( histoBooked_ )  H_status_->Fill( 8 );
		CsErrLog::Instance()->mes(elError,"CsDetFamily::_minimise. Singular Matrix.");
		return false;
	}

	// calculate x, y, tx, ty
	x = 0; 	
	y = 0; 
	tx = 0;	
	ty = 0;
	for( unsigned int i = 1; i<= 4; i++ ){
		x  += A(1,i)*B(i,1);
		y  += A(2,i)*B(i,1);
		tx += A(3,i)*B(i,1);
		ty += A(4,i)*B(i,1);
	}
  
  // calculate chi2
  chi2 = 0;
	for( unsigned int i = 0; i<d_.size(); i++ ) {
    if( ic_[i] < fc_[i] || ic_[i] >= fc_[i]+nc_[i]  || 
      (( mask | bigMask_ ) & ( 1 << i )) ) continue;
    
    register unsigned int ic = ic_[i];
    chi2 += pow( 
      (x+(wc_[ic]-z_)*tx) * iR_[i](1,1) +     
      (y+(wc_[ic]-z_)*ty) * iR_[i](1,2) +  
      wc_[ic]*iR_[i](1,3)-uc_[ic] , 2 ) / rc2_[ic];  
    
    if( useClusterAsc_ && ic_A_[ic] >= 0 ) {
      register unsigned int ica = (unsigned int) ic_A_[ic];
      chi2 += pow( 
        (x+(wc_A_[ica]-z_)*tx) * iR_[ cd_A_[ica] ](1,1) +     
        (y+(wc_A_[ica]-z_)*ty) * iR_[ cd_A_[ica] ](1,2) +  
        wc_A_[ica]*iR_[ cd_A_[ica] ](1,3)-uc_A_[ica] , 2 ) / rc2_A_[ica];      
    } 
    
  }
  
  // correct chi2 from number of degrees of freedom
  int ndf = nc - 4;
  if( ndf > 0 ) chi2 /= (double) ndf;
  else {
    if( histoBooked_ ) H_status_->Fill( 7 );
    chi2 = -1;
  }
    
  return true;
}

//_____________________________________________________________________________
void CsDetFamily::_fillMatrix( void )
{
  for( unsigned int i = 0; i < c_.size(); i++ ) {
    register unsigned int j = cd_[i];   
    A_.push_back( _getA( iR_[j], uc_[i], wc_[i], rc2_[i] ) );
    B_.push_back( _getB( iR_[j], uc_[i], wc_[i], rc2_[i] ) );
  }    

  // fill matrix for associated clusters
  if( useClusterAsc_ )
  for( unsigned int i = 0; i< c_A_.size(); i++ ) {
    register unsigned int j = cd_A_[i];
    A_A_.push_back( _getA( iR_[j], uc_A_[i], wc_A_[i], rc2_A_[i]) );
    B_A_.push_back( _getB( iR_[j], uc_A_[i], wc_A_[i], rc2_A_[i]) );
  }

  

  return; 
}  

//_____________________________________________________________________________
HepMatrix CsDetFamily::_getA( HepMatrix iR, double uc, double wc, double rc2 )
{
  HepMatrix A(4,4,0);

	// first line (dChi2/dx)
  A(1,1) += iR(1,1)*iR(1,1)/rc2;           // x
  A(1,2) += iR(1,1)*iR(1,2)/rc2;           // y
  A(1,3) += iR(1,1)*iR(1,1)*(wc-z_)/rc2;   // tx
  A(1,4) += iR(1,1)*iR(1,2)*(wc-z_)/rc2;   // tx

	// second line (dChi2/dy)
  A(2,2) += iR(1,2)*iR(1,2)/rc2;           // y
  A(2,4) += iR(1,2)*iR(1,2)*(wc-z_)/rc2;   // ty

  // third line (dChi2/dtx)
  A(3,3) += iR(1,1)*(wc-z_)*iR(1,1)*(wc-z_)/rc2; // tx
  A(3,4) += iR(1,1)*(wc-z_)*iR(1,2)*(wc-z_)/rc2; // ty

  // fourth line (dChi2/dty)
  A(4,4) += iR(1,2)*(wc-z_)*iR(1,2)*(wc-z_)/rc2; // ty

  // complete A, using simetry properties
  A(2,1) = A(1,2);
  A(2,3) = A(1,4);
  A(3,1) = A(1,3);
  A(3,2) = A(2,3);
  A(4,1) = A(1,4); 
  A(4,2) = A(2,4); 
  A(4,3) = A(3,4);
  
  return A;
}
    
//_____________________________________________________________________________
HepMatrix CsDetFamily::_getB( HepMatrix iR, double uc, double wc, double rc2 )
{
 	HepMatrix B(4,1,0);	
  B(1,1) += iR(1,1)*(uc-wc*iR(1,3))/rc2;
  B(2,1) += iR(1,2)*(uc-wc*iR(1,3))/rc2;
  B(3,1) += iR(1,1)*(wc-z_)*(uc-wc*iR(1,3))/rc2;
  B(4,1) += iR(1,2)*(wc-z_)*(uc-wc*iR(1,3))/rc2;

  return B;
}


      
