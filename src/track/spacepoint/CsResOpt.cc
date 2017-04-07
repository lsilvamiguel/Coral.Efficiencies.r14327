// $Id: CsResOpt.cc,v 1.14 2010/06/18 10:44:22 tnagel Exp $

/*!
   \file    CsResOpt.cc
   \brief   Compass Resolution optimiser from Space point Residuals
   \author  Hugo Pereira
   \version $Revision: 1.14 $
   \date    $Date: 2010/06/18 10:44:22 $
*/

#include "CsResOpt.h"
#include "Coral.h"
#include "CsInit.h"
#include "CsOpt.h"
#include "CsErrLog.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CsHistograms.h"
#include "CsSPUtils.h"
#include "CsDetFamily.h"
#include "CsDetector.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMinuit.h"

using namespace std;
using namespace CLHEP;

static CsResOpt* resOpt;

//_____________________________________________________________________________
CsResOpt::CsResOpt( const CsDetFamily& df ): 
  df_ (&df),
  residOK_( false ), 
  resIOK_( false ), 
  resFOK_( false )
{
  ;
  resid_.clear();

  // test root Implementation
  if( !(CsHistograms::GetImplementation()==CsHistograms::ROOT) ) 
  CsErrLog::Instance()->mes(elFatal,"CsResOpt::CsResOpt. Root package needed.");

  // init Fit resolutions
  resFAve_ = 0;
  resF_.clear();
  
  // try getting initial resolutions
  resIAve_ = 0;
  resI_.clear();
  vector< CsDetector* > d = df_->getDetectors(); 
  if( d.size() == 0 ) 
  CsErrLog::Instance()->mes(elFatal,"CsResOpt::CsResOpt. Empty det family.");
  
  // try getting starting point resolutions
  for( unsigned int i=0; i< d.size(); i++ ){
    double res = d[i]->getVel()*d[i]->getSpSli();
    resI_.push_back( res );
    resIAve_+=res;
  }
  resIAve_ /= d.size();
  resIOK_ = true;
  resOpt = NULL;
}

//_____________________________________________________________________________
bool CsResOpt::setResidsFromHists( void )
{
  // get histograms
  if( !df_->calHistoBooked() ) {
    CsErrLog::Instance()->mes(elError,"CsResOpt::SetResidsFromHists. detFamily calibration histograms not booked.");
    return false;
  }
  
  vector< CsHist1S* > h=df_->getResHists( );
  if( h.size() != df_->detSize() ){ 
    CsErrLog::Instance()->mes(elError,"CsResOpt::SetResidsFromHists. Histograms and Detector sizes do not match.");
    return false;
  }
  
  residOK_ = false;
  resid_.clear();
  
  cout << endl << "CsResOpt::SetResidsFromHists. " << df_->getName().c_str() << endl;
  for( unsigned int i=0; i< h.size(); i++ ){

    TH1S* hRoot = reinterpret_cast<TH1S*>(h[i]->ThePointer());
    
    if( hRoot == NULL ) {
      CsErrLog::Instance()->mes(elError,"CsResOpt::SetResidsFromHists. Troubles to get Root histogram.");
      resid_.clear();
      return false;
    }
    
    hRoot->Fit("gaus","0Q");
    double resid =  hRoot->GetFunction("gaus")->GetParameter( 2 );
    cout << "Residual for " << df_->getDetectors()[i]->GetTBName().c_str() << ": " << resid << endl;
    resid_.push_back( resid );
  }
  residOK_ = true;
  return true;
}

//_____________________________________________________________________________
bool CsResOpt::setResidsFrom2DHists( void )
{
  // get histograms
  if( !df_->calHistoBooked() ) {
    CsErrLog::Instance()->mes(elError,"CsResOpt::setResidsFrom2DHists. detFamily calibration histograms not booked.");
    return false;
  }
  
  vector< CsHist2S* > h=df_->getRes2DHists( );
  if( h.size() != df_->detSize() ){ 
    CsErrLog::Instance()->mes(elError,"CsResOpt::setResidsFrom2DHists. Histograms and Detector sizes do not match.");
    return false;
  }
  
  residOK_ = false;
  resid_.clear();
  
  cout << endl << "CsResOpt::setResidsFrom2DHists. " << df_->getName().c_str() << endl;
  for( unsigned int i=0; i< h.size(); i++ ){

    TH2S* hRoot = reinterpret_cast<TH2S*>(h[i]->ThePointer());
    
    if( hRoot == NULL ) {
      CsErrLog::Instance()->mes(elError,"CsResOpt::setResidsFrom2DHists. Troubles to get Root histogram.");
      resid_.clear();
      return false;
    }
    
    unsigned int n = hRoot->GetNbinsX();
    unsigned int n1 = n/2;
    unsigned int n2 = (n%2) ? ( n/2 + 2 ): ( n/2 + 1 );
    char rName[128];
    char lName[128];
    sprintf( lName,"%s_l",hRoot->GetName() );
    sprintf( rName,"%s_r",hRoot->GetName() );
    TH1D* hl = hRoot->ProjectionY( lName, 1, n1 );
    TH1D* hr = hRoot->ProjectionY( rName, n2, n );
    hl->Fit("gaus","0Q");  hr->Fit("gaus","0Q");
    double sigl = hl->GetFunction("gaus")->GetParameter( 2 );
    double sigr = hr->GetFunction("gaus")->GetParameter( 2 );
    double resid =  0.5*(sigl+sigr);
    cout << "Residual for " << df_->getDetectors()[i]->GetTBName().c_str() 
      << ". L:"  << sigl
      << " R:"   << sigr
      << " ave:" << resid
      << endl;
      
    resid_.push_back( resid );
  }
  residOK_ = true;
  return true;
}

//_____________________________________________________________________________
bool CsResOpt::setResidsByHand( vector< double > resid )
{
  if( resid.size() != df_->detSize() ) {
    CsErrLog::Instance()->mes(elError,"CsResOpt::SetResidsByHand. resid and det sizes do not match.");
    return false;
  }

  resid_.clear();
  for( unsigned int i = 0; i < resid.size(); i++ ) resid_.push_back( resid[i] );
  residOK_ = true;
  return true;
}    

//_____________________________________________________________________________
bool CsResOpt::setResolutionsByHand( vector< double > res )
{
  if( res.size() != df_->detSize() ) {
    CsErrLog::Instance()->mes(elError,"CsResOpt::SetResolutionsByHand. res and det sizes do not match.");
    return false;
  }

  resIAve_ = 0;
  resI_.clear();
  for( unsigned int i = 0; i < res.size(); i++ ) {
    resIAve_+= res[i];
    resI_.push_back( res[i] );
  }
  
  resIAve_ /= res.size();
  resIOK_ = true;
  return true;
}    

//_____________________________________________________________________________
void CsResOpt::testResCalc( void )
{

  // check resolutions
  if( ! resIOK_ ) {
    CsErrLog::Instance()->mes(elError,"CsResOpt::testResCalc. Resolutions not set.");
    return;
  }
  
  cout << endl << "CsResOpt::testResCalc " << df_->getName().c_str() << endl;
  cout << "resAve=" << resIAve_ << "mm." << endl;
  
  vector< CsDetector* > d = df_->getDetectors();
  for( unsigned int i=0; i<d.size(); i++ )
  cout << d[i]->GetTBName().c_str() 
    << ": res=" << resI_[i] 
    << "mm, sglresid=" << resCalc( i, resIAve_ ) 
    << "mm, mltresid=" << resCalc( i, resI_ )
    << endl;
  return;
}

//___________________________________________________________________________
bool CsResOpt::fitResAve( double& resAve, double& resAveErr  )
{
  
  // check resolutions
  if( ! resIOK_ ) {
    CsErrLog::Instance()->mes(elError,"CsResOpt::fitResAve. Resolutions not set.");
    return false;
  }

  // check residuals
  if( ! residOK_ ) {
    CsErrLog::Instance()->mes(elError,"CsResOpt::fitResAve. Residuals not set.");
    return false;
  }
 
  resAve = resIAve_;
  resAveErr = 0;
  resOpt = this;
  
  TMinuit *gMinuit = new TMinuit(1); 
  gMinuit->mninit(5,6,7);     // Logical units
  gMinuit->SetFCN(fcnAve);

  cout << "fitResAve. TMinuit initialization performed." << endl;
  
  double dSig = .00001;
  int error;
  gMinuit->mnparm(0,"ResAve",resAve,dSig, 0,0,error);
  if (error) CsErrLog::Instance()->mes(elFatal,"CsResOpt::fitResAve. Troubles defining parameters.");

  double flag0 = 0, flag1 = 1, flag3 = 3;
  gMinuit->mnexcm("CALL FCN",  &flag1 ,1,error);
  gMinuit->mnexcm("SET PRINT", &flag0 ,1,error);
  gMinuit->mnexcm("SET NOGradient", &flag0 ,1,error);
  gMinuit->mnexcm("MIGRAD",   &flag0 ,0,error);
  gMinuit->mnexcm("MINOS",    &flag0 ,0,error);
  gMinuit->mnexcm("CALL FCN", &flag3 ,1,error);
  
  cout << endl << "--------------------------------------------" << endl;
  cout << "fitResAve. TMinuit minimisation performed." << endl;
  
  gMinuit->GetParameter(0, resAve, resAveErr);
  resFAve_ = resAve;

  cout << "fitResAve. aveRes: " << resAve << " (" 
    << resAveErr << ") mm." 
    << endl;
  cout << "fitResAve. Residual <begin> <end>" << endl;
  for( unsigned int i=0; i<df_->detSize(); i++ )
  cout << df_->getDetectors()[i]->GetTBName().c_str() 
    << " resid= " << resid_[i]
    << " vs " << sqrt( pow( resCalc( i, resAve ), 2 ) + pow( resAve, 2 ) )
    << endl;
  
  return true;
}

//___________________________________________________________________________
bool CsResOpt::fitRes( FILE* resOut )
{
  
  double resAve, resAveErr;
  
  #define USE_RES_AVE
  #ifdef USE_RES_AVE
  if( !fitResAve( resAve, resAveErr ) ) { 
    CsErrLog::Instance()->mes(elError,"CsResOpt::fitRes. Troubles with fitResAve()" );
    return false;
  }
  #else
  resAveErr = -1.000;
  resAve = resIAve_;
  resOpt = this;
  #endif

  vector<CsDetector* > d = df_->getDetectors();
  TMinuit *gMinuit = new TMinuit((int)d.size()); 
  gMinuit->mninit(5,6,7);     // Logical units
  gMinuit->SetFCN(fcn);

  cout << "fitRes. TMinuit initialization performed." << endl;
  
  double dSig = .00001;
  int error;
    
  for( unsigned int i=0; i< d.size(); i++ ) {
    char pName[15];
    sprintf(pName,"res%s",d[i]->GetTBName().c_str() );
    gMinuit->mnparm(i,pName,resAve,dSig, 0,0,error);
    if (error) CsErrLog::Instance()->mes(elFatal,"fitRes. Troubles defining parameters.");
  }
  
  
  double flag0 = 0, flag1 = 1, flag3 = 3;
  gMinuit->mnexcm("CALL FCN",  &flag1 ,1,error);
  gMinuit->mnexcm("SET PRINT", &flag0 ,1,error);
  gMinuit->mnexcm("SET NOGradient", &flag0 ,1,error);
  gMinuit->mnexcm("MIGRAD",   &flag0 ,0,error);
  gMinuit->mnexcm("MINOS",    &flag0 ,0,error);
  gMinuit->mnexcm("CALL FCN", &flag3 ,1,error);
  
  cout << endl << "--------------------------------------------" << endl;
  cout << "fitRes. TMinuit minimisation performed." << endl;
  
  vector< double > resV;    resV.clear();
  vector< double > resErrV; resErrV.clear();
  resF_.clear();
  
  for( unsigned int i=0; i<d.size(); i++ ) {
    double res, resErr;
    gMinuit->GetParameter(i, res, resErr);
    resV.push_back( res );
    resErrV.push_back( resErr );
    resF_.push_back( res );
  }
  
  if( resOut == NULL ) resOut = stdout;
  fprintf(resOut,"\n//%s %s - CsResOpt::fitRes. DF:%s Run:%i maxEvt:%i \n", 
    CsSPUtils::Instance()->getDate().c_str(), 
    CsSPUtils::Instance()->getTime().c_str(),
    df_->getName().c_str(),
    Coral::Instance()->getRunNumber(),
    CsInit::Instance()->getMaxEvents() );
  fprintf(resOut,"%8s %10.5f %10.5f %10.5f %10.5f %10.5f\n",
    "AVERAGE", 
    resAve,  resAveErr, 
    -1.0,
    -1.0,  
    -1.0 );
  for( unsigned int i=0; i<d.size(); i++ ) 
  fprintf(resOut,"%8s %10.5f %10.5f %10.5f %10.5f %10.5f\n",
    d[i]->GetTBName().c_str(), 
    resV[i],  resErrV[i], 
    resCalc( i, resV ),
    sqrt( pow( resCalc( i, resV ), 2 ) + pow( resV[i], 2 ) ),
    resid_[i]  
  );
  
  fprintf(resOut, "%8s %10s %10s %10s %10s %10s\n\n", "TBName", "res", "error", "tr_res", "residF" ,"residI");
  
  return true;
}

//_____________________________________________________________________________
double CsResOpt::resCalc( unsigned int it, double res )
{
  vector< CsDetector* > d = df_->getDetectors();
  vector< double > resV; resV.clear();
  for( unsigned int i =  0; i< d.size(); i++ ) resV.push_back( res );
  return resCalc( it, resV );
}    

//_____________________________________________________________________________
double CsResOpt::resCalc( unsigned int it, vector< double > res )
{
  vector< CsDetector* > d = df_->getDetectors();
  
  if( res.size() != d.size() ) 
  CsErrLog::Instance()->mes(elFatal,"CsResOpt::resCalc. resolution and detector sizes does not match. ");
  
  if( it >= d.size() ) 
  CsErrLog::Instance()->mes(elFatal,"CsResOpt::resCalc. Wrong detector index. ");
  
  HepMatrix iV_Fast(2,2,0);
  HepMatrix iV(4,4,0);
  
  double tz = d[it]->getZcm();
  bool mode = df_->useFullMinimisation();
  if( ! mode && !df_->useFastMinimisation() ) 
  CsErrLog::Instance()->mes(elFatal,"CsResOpt::resCalc. unrecognized minimisation procedure for det family. ");
    
  for( unsigned int i = 0; i < d.size(); i++ ) {
    if( it == i || ( df_->removeADet() &&  df_->getADetIndexes()[it] == (int) i ) ) continue;
    HepMatrix iR = d[i]->getRotWRSInv();
    double r2 = pow( res[i], 2 );
    if( mode ) {      // Resolution corresponding to full minimisation  
      double rz = d[i]->getZcm();
      iV(1,1) += iR(1,1)*iR(1,1)/r2;
      iV(1,2) += iR(1,1)*iR(1,2)/r2;
      iV(1,3) += iR(1,1)*iR(1,1)*(rz-tz)/r2;
      iV(1,4) += iR(1,1)*iR(1,2)*(rz-tz)/r2;

      iV(2,2) += iR(1,2)*iR(1,2)/r2;   
      iV(2,4) += iR(1,2)*iR(1,2)*(rz-tz)/r2;

      iV(3,3) += iR(1,1)*(rz-tz)*iR(1,1)*(rz-tz)/r2;   
      iV(3,4) += iR(1,1)*(rz-tz)*iR(1,2)*(rz-tz)/r2;   
 
      iV(4,4) += iR(1,2)*(rz-tz)*iR(1,2)*(rz-tz)/r2;   
    } else {         // Resolution corresponding to fast minimisation
      iV_Fast(1,1) += iR(1,1)*iR(1,1)/r2;
      iV_Fast(2,2) += iR(1,2)*iR(1,2)/r2;
      iV_Fast(1,2) += iR(1,1)*iR(1,2)/r2;
      iV_Fast(2,1) += iR(1,1)*iR(1,2)/r2;
    }
  }
  
  if( mode ){   // complete iV using symetries in case of full minimisation resolution    
    iV(2,1) = iV(1,2);
    iV(2,3) = iV(1,4);
    iV(3,1) = iV(1,3);
    iV(3,2) = iV(2,3);
    iV(4,1) = iV(1,4);
    iV(4,2) = iV(2,4);
    iV(4,3) = iV(3,4);
  }
  
  int err;
  if( mode ) {   //Full minimisation 
    
    HepMatrix V = iV.inverse( err );
    if( err ) {
		  CsErrLog::Instance()->mes(elError,"CsResOpt::resCalc. Singular Matrix.");
      return -1;
    } else {
      HepMatrix iR = d[it]->getRotWRSInv();
      return sqrt( 
        V(1,1)*iR(1,1)*iR(1,1)
        + V(2,2)*iR(1,2)*iR(1,2)
        + 2*V(1,2)*iR(1,1)*iR(1,2)
      );
    }
    
  } else {       //Fast minimisation
    
    HepMatrix V_Fast = iV_Fast.inverse( err );
    if( err ) {
		  CsErrLog::Instance()->mes(elError,"CsResOpt::resCalc. Singular Matrix.");
      return -1;
    } else {
      HepMatrix iR = d[it]->getRotWRSInv();
      return sqrt( 
        V_Fast(1,1)*iR(1,1)*iR(1,1)
        + V_Fast(2,2)*iR(1,2)*iR(1,2)
        + 2*V_Fast(1,2)*iR(1,1)*iR(1,2)
      );
    }

  }
}    

//_____________________________________________________________________________
double CsResOpt::getZoptX( double res ) {
  vector< double > resV; resV.clear();
  for( unsigned int i = 0; i <= df_->getDetectors().size(); i++ ) resV.push_back( res );
  return getZoptX( resV );
}

//_____________________________________________________________________________
double CsResOpt::getZoptX( vector< double > res )
{
  vector< CsDetector* > d = df_->getDetectors();
  if( d.size() < 5 ) 
  CsErrLog::Instance()->mes( elFatal, "CsResOpt::getZoptX. Zopt is meaningfull only for full minimisation. At least 5 detectors required.\n");
  
  if( res.size() != d.size() ) 
  CsErrLog::Instance()->mes(elFatal,"CsResOpt::getZoptX. resolution and detector sizes does not match. ");
  
  vector< double > resXTrs; resXTrs.clear();
  
  // get resolutions projected on X
  unsigned int nPl = 0;
  for( unsigned int i = 0; i < d.size(); i++ ) {
    HepMatrix iR = d[i]->getRotWRSInv();
    if( fabs(iR(1,1))>0.0001 ) {
      resXTrs.push_back( res[i]/fabs(iR(1,1)) ); 
      nPl ++;
    } else resXTrs.push_back( -1 );
  }
  cout << "CsResOpt::getZoptX. " << nPl << " detectors OK." << endl;
  
  double A = 0, B = 0; // analytic terms
  
  for( unsigned int i = 0; i < d.size(); i++ ) {
    if( resXTrs[i] < 0 ) continue;
    double SZijj=0, SZij = 0;
    double zi = d[i]->getZcm();
    for( unsigned int j = 0; j < d.size(); j++ ) {
      if( resXTrs[j] < 0 ) continue;
      double zj = d[j]->getZcm();
      SZij  +=  ( zi - zj )/pow( resXTrs[j], 2 );
      SZijj += zj*( zi - zj )/pow( resXTrs[j], 2 );
    }   
    A += pow( SZij/resXTrs[i], 2 );
    B += SZij*SZijj/pow( resXTrs[i], 2 );
    
  }
  
  return B/A;
}

//_____________________________________________________________________________
double CsResOpt::getZoptY( double res ) {
  vector< double > resV; resV.clear();
  for( unsigned int i = 0; i <= df_->getDetectors().size(); i++ ) resV.push_back( res );
  return getZoptY( resV );
}
  
//_____________________________________________________________________________
double CsResOpt::getZoptY( vector< double > res )
{
  vector< CsDetector* > d = df_->getDetectors();
  if( d.size() < 5 ) 
  CsErrLog::Instance()->mes( elFatal, "CsResOpt::getZoptX. Zopt is meaningfull only for full minimisation. At least 5 detectors required.\n");
  
  if( res.size() != d.size() ) 
  CsErrLog::Instance()->mes(elFatal,"CsResOpt::getZoptX. resolution and detector sizes does not match. ");
  
  vector< double > resYTrs; resYTrs.clear();
  
  // get resolutions projected on Y
  unsigned int nPl = 0; 
  for( unsigned int i = 0; i < d.size(); i++ ) {
    HepMatrix iR = d[i]->getRotWRSInv();
    if( fabs(iR(1,2))>0.0001 ) {
      resYTrs.push_back( res[i]/fabs(iR(1,2)) ); 
      nPl ++;
    } else resYTrs.push_back( -1 );
  }
  cout << "CsResOpt::getZoptY. " << nPl << " detectors OK." << endl;
  
  double A = 0, B = 0; // analytic terms
  
  for( unsigned int i = 0; i < d.size(); i++ ) {
    if( resYTrs[i] < 0 ) continue;
    double SZijj=0, SZij = 0;
    double zi = d[i]->getZcm();
    for( unsigned int j = 0; j < d.size(); j++ ) {
      if( resYTrs[j] < 0 ) continue;
      double zj = d[j]->getZcm();
      SZij  +=  ( zi - zj )/pow( resYTrs[j], 2 );
      SZijj += zj*( zi - zj )/pow( resYTrs[j], 2 );
    }   
    A += pow( SZij/resYTrs[i], 2 );
    B += SZij*SZijj/pow( resYTrs[i], 2 );
    
  }
  
  return B/A;
}

//_____________________________________________________________________________
bool CsResOpt::dumpResAtZ( double z, double res )
{
  vector< double > resV; resV.clear();
  for( unsigned int i = 0; i < df_->getDetectors().size(); i++ ) resV.push_back( res );
  return dumpResAtZ( z, resV );

}

//_____________________________________________________________________________
bool CsResOpt::dumpResAtZ( double z, vector< double > res )
{
  vector< CsDetector* > d = df_->getDetectors();
  if( d.size() < 5 ) 
  CsErrLog::Instance()->mes( elFatal, "CsResOpt::dumpResAtZ is meaningfull only for full minimisation. At least 5 detectors required.\n");
  
  if( res.size() != d.size() ) 
  CsErrLog::Instance()->mes(elFatal,"CsResOpt::dumpResAtZres. resolution and detector sizes does not match. ");
  HepMatrix iV(4,4,0);
  for( unsigned int i = 0; i < d.size(); i++ ) {
    HepMatrix iR = d[i]->getRotWRSInv();
    double r2 = pow( res[i], 2 );
    double rz = d[i]->getZcm();
    iV(1,1) += iR(1,1)*iR(1,1)/r2;
    iV(1,2) += iR(1,1)*iR(1,2)/r2;
    iV(1,3) += iR(1,1)*iR(1,1)*(rz-z)/r2;
    iV(1,4) += iR(1,1)*iR(1,2)*(rz-z)/r2;
    iV(2,2) += iR(1,2)*iR(1,2)/r2;   
    iV(2,4) += iR(1,2)*iR(1,2)*(rz-z)/r2;
    iV(3,3) += iR(1,1)*(rz-z)*iR(1,1)*(rz-z)/r2;   
    iV(3,4) += iR(1,1)*(rz-z)*iR(1,2)*(rz-z)/r2;   
    iV(4,4) += iR(1,2)*(rz-z)*iR(1,2)*(rz-z)/r2;   
  }
  iV(2,1) = iV(1,2);
  iV(2,3) = iV(1,4);
  iV(3,1) = iV(1,3);
  iV(3,2) = iV(2,3);
  iV(4,1) = iV(1,4);
  iV(4,2) = iV(2,4);
  iV(4,3) = iV(3,4);
  
  int err;
  HepMatrix V = iV.inverse( err );
  if( err ) {
		CsErrLog::Instance()->mes(elError,"resCalc. Singular Matrix.");
    return false;
  } 
  
  fprintf(stdout,"\n//%s %s - CsResOpt::dumpResAtZ. DF:%s Run:%i maxEvt:%i\n", 
    CsSPUtils::Instance()->getDate().c_str(), 
    CsSPUtils::Instance()->getTime().c_str(),
    df_->getName().c_str(),
    Coral::Instance()->getRunNumber(),
    CsInit::Instance()->getMaxEvents() 
  );
  
  fprintf(stdout,"z = %f\n",z);
  fprintf(stdout,"resX  = %f\n",sqrt(V(1,1)) );
  fprintf(stdout,"resY  = %f\n",sqrt(V(2,2)));
  fprintf(stdout,"resTx = %f\n",sqrt(V(3,3)));
  fprintf(stdout,"resTY = %f\n",sqrt(V(4,4)));
  
  return true;
}
 
//_____________________________________________________________________________
void fcnAve( int &npar, double *gin, double &f, double *x, int flag )
{

  if( resOpt == NULL ) 
  CsErrLog::Instance()->mes(elFatal,"fcnAve. resOpt is NULL.");

  f = 0;
  for( unsigned int i=0; i<resOpt->getDetectors().size(); i++ )
  f += pow( resOpt->getResids()[i] - sqrt( pow( resOpt->resCalc( i, x[0] ),2 )+ pow( x[0],2 ) ), 2 ); 
  f *= 100;
  return;
}

//_____________________________________________________________________________
void fcn( int &npar, double *gin, double &f, double *x, int flag )
{
  vector< double > res; res.clear();
  for( unsigned int i=0; i<resOpt->getDetectors().size(); i++ ) res.push_back( x[i] );

  if( resOpt == NULL ) 
  CsErrLog::Instance()->mes(elFatal,"fcn. resOpt is NULL.");
  
  f = 0;
  for( unsigned int i=0; i<resOpt->getDetectors().size(); i++ )
  f += pow( resOpt->getResids()[i] - sqrt( pow( resOpt->resCalc( i, res ),2 )+ pow( res[i],2 ) ), 2 ); 
  f *= 100;
  return;
}
