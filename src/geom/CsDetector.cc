// $Id: CsDetector.cc,v 1.61 2010/07/12 13:01:19 tnagel Exp $

/*!
   \file    CsDetector.cc
   \brief   Compass Generic Tracking Detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.61 $
   \date    $Date: 2010/07/12 13:01:19 $
*/

#include <math.h>
#include <string.h>
#include <strings.h>

#include <cstdlib>
#include <functional>

#include "CsDetector.h"
#include "CsZone.h"
#include "CsDigit.h"
#include "CsOpt.h"
#include "CsEvent.h"
#include "CsDigit.h"
#include "CsMCDigit.h"
#include "CsComgNtCommons.h"
#include "CsGeom.h"
#include "CsMCTrkHit.h"
#include "CsGeant3.h"
#include "CsMCTrack.h"

using namespace std;
using namespace CLHEP;
using CS::DetID;

extern QhitType Qhit;

CsDetector::CsDetector( const int    row,
		        const int    id,    const char* name,    const char *TBname,
			const int    unit,  const int    type,
			const double rdLen, const double xsiz,  
			const double ysiz,  const double zsiz,
			const double xcm,   const double ycm,   
			const double zcm,   const HepMatrix rotDRS,
			const HepMatrix rotWRS,
			const double wirD,  const double ang,   
			const int    nWir,  const double wirP, 
			const double eff,   const double bkg,
			const double tGate ) :
  CsDet( DetID(name,id), TBname ),

  // histogram switches initialisation
  hLevel_( None ),

  row_(row), unit_(unit), type_(type), rdLen_(rdLen), 
  xsiz_(xsiz), ysiz_(ysiz), zsiz_(zsiz), xcm_(xcm), ycm_(ycm), zcm_(zcm),
  wirD_(wirD), ang_(ang), nWir_(nWir), wirP_(wirP), hasVarP_(false), eff_(eff),
  bkg_(bkg), tGate_(tGate),
  
  // Dead zone initialization: point-like rectangular...
  dZtype_(0x11),dZdim_(0),dZydim_(0),
  dZxdrs_(0), dZydrs_(0),      // ... w/ same position
  rotD2DZ_(HepMatrix(3,3,1))   // ...and same orientation as detector's
{

  ids_.insert(id);
  xcms_[id]=xcm;
  ycms_[id]=ycm;
  zcms_[id]=zcm;

  _xshift = 0;
  _yshift = 0;
  _deltax = 0;
  _deltay = 0;
  sinAng_ = sin( ang_ * (M_PI) / 180. );
  cosAng_ = cos( ang_ * (M_PI) / 180. );
  rotDRS_ = rotDRS; 
  rotWRS_ = rotWRS;    
  int err; 
  rotWRS_inv_ = rotWRS.inverse(err);    
  myZones_.clear();
  myClusters_.clear();

  wirPs_[0] = wirP_;
 
  HepMatrix m(3,3);
  m(1,1) =  getCosAng(); m(1,2) = -getSinAng(); m(1,3) = 0;
  m(2,1) =  getSinAng(); m(2,2) = getCosAng(); m(2,3) = 0;
  m(3,1) =            0; m(3,2) =           0; m(3,3) = 1;
  rotDRS2WRS_ = m;
  
  decodingDone_ = false;
  decode_ = false;
  clusteringDone_ = false;
  
}

//___________________________________________________________________________
void CsDetector::ReadHistLevel( void )
{
  
  //=== Look for histogram keys:
  string tbn = GetTBName();
  string tb( tbn, 0, 2 );
  
  string hLevStr("none");
  if( CsOpt::Instance()->getOpt(tbn.c_str(), "hist level", hLevStr ) ||
      CsOpt::Instance()->getOpt( tb.c_str(), "hist level", hLevStr ) ) {
    
    if( hLevStr == "none" )        hLevel_ = None;
    else if( hLevStr == "normal" ) hLevel_ = Normal;
    else if( hLevStr == "high" )   hLevel_ = High;
    else CsErrLog::msg(elError, __FILE__, __LINE__, 
      "\"%s\" - Unrecognised histogram level: \"%s\". Default is \"none\".",
      GetTBName().c_str(),
      hLevStr.c_str() );
      
  } else hLevel_ = None;

  if( hLevel_ > None ) 
  cout << "CsDetector::CsDetector - INFO: Histogram level for \"" << tbn
       << "\" is " << hLevStr 
       << " (" << hLevel_ << ").\n";

}  

//___________________________________________________________________________
bool CsDetector::operator==( const CsDetector& det ) const
{
  return GetID() == det.GetID();
}

bool CsDetector::operator<( const CsDetector& det ) const
{
  if( zcm_ < det.zcm_ )
    return( true );
  else
    if( zcm_ == det.zcm_ )
      return GetID().GetNumber() < det.GetID().GetNumber();
    else
      return false;
}


void CsDetector::AddSubDetector( const int    row,
				 const int    id,    const char* name,    
				 const char *TBname,
				 const int    unit,  const int    type,
				 const double rdLen, const double xsiz,  
				 const double ysiz,  const double zsiz,
				 const double xcm,   const double ycm,   
				 const double zcm,   const HepMatrix rotDRS,
				 const HepMatrix rotWRS,
				 const double wirD,  const double ang,   
				 const int    nWir,  const double wirP, 
				 const double eff,   const double bkg,
				 const double tGate ) {

  typedef map<int,double>::iterator IP;


  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"***********   ADDING SUBDETECTOR TO %s",GetTBName().c_str());
  
  // these parameters must be the same in all subdetectors.

  if(strcmp(TBname,GetTBName().c_str())!=0)
    throw CS::Exception("CsDetector::addSubDetector : bad TBName"); 
  if(type!=getType())
    throw CS::Exception("CsDetector::addSubDetector : bad type");
  if(rdLen-rdLen_)
    throw CS::Exception("CsDetector::addSubDetector : bad rdLen");
  if(!(rotDRS==rotDRS_))
    throw CS::Exception("CsDetector::addSubDetector : bad rotDRS_");
  if(!(rotWRS==rotWRS_))
    throw CS::Exception("CsDetector::addSubDetector : bad rotWRS_");
  if(ang!=ang_)
    throw CS::Exception("CsDetector::addSubDetector : bad wire angle");
  if(eff!=eff_)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "%s: Efficiency: Newly entered %f != Previous %f",TBname,eff,eff_);
  if(bkg!=bkg_)
    throw CS::Exception("CsDetector::addSubDetector : bad bkg");
  if(tGate!=tGate_)
    throw CS::Exception("CsDetector::addSubDetector : bad tGate");


  //          ***** SIZE of THE OVERALL DETECTOR *****
  // (The size in the dimension along the measured coordinate axis (and hence
  // along the axis of the shift from one sub-detector to the next) has to be
  // updated. It's either "xsiz_" or "ysiz_", depending upon whether the DRS is
  // (or, rather, tends to be) parallel (case of all MMs and 'X' hodos) or
  // perpendicular (case of the 'Y' hodos) to the WRS. In order not to restrict
  // the code to the above-listed cases and to make it more general (cases where
  // the wires are tilted (by design that is, as opposed to small tilts due to
  // misalignment) w.r.t. the DRS axes are still not provided for), the angle
  // between the shift axes in DRS and WRS is computed.)
  double x[3], xdrs[3];
  rotateVectMRS2WRS(xcm-xcm_,ycm-ycm_,zcm-zcm_, x[0],   x[1],   x[2]);
  rotateVectMRS2DRS(xcm-xcm_,ycm-ycm_,zcm-zcm_, xdrs[0],xdrs[1],xdrs[2]);
  if (fabs(x[0]*xdrs[0]+x[1]*xdrs[1]+x[2]*xdrs[2])/(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])>.04) { // Cutting on the angle being 90+/-2.3 to allow for some misalignment
    xsiz_ += xsiz;
#define CsDetector_XCHECK_SIZE_UPDATE
#ifdef CsDetector_XCHECK_SIZE_UPDATE
    // Since failing to properly handle this size update leads to a error that
    // does not manifests itself conspicuously (and the case already happened
    // once, w/ an older, cruder version of the code), let's X-check it, on the
    // MM and H cases. If, at some point, some other variable pitch detectors
    // are introduced, one would possibly have to undefine the present macro.
    if (TBname[0]=='H' && TBname[4]=='Y')
      CsErrLog::msg(elFatal,__FILE__,__LINE__,"%s: Size update goes wrong", TBname);
#endif
  }
  else {
    ysiz_ += ysiz;
#ifdef CsDetector_XCHECK_SIZE_UPDATE
    if ( TBname[0]=='M' || (TBname[0]=='H' && TBname[4]=='X') )
      CsErrLog::msg(elFatal,__FILE__,__LINE__,"%s: Size update goes wrong", TBname);
#endif
  }

  //               ***** CONTINUITY TEST (in WRS) ***** 
  double tolerance = wirP/10.;

  double halfxs = 0;  // Old half size of the active zone (WRS)
  int lastw = 0; double lastp = 0;
  for (IP i = wirPs_.begin(); i!=wirPs_.end(); i++) {
    halfxs += (i->first-lastw)*lastp;
    lastw = i->first; lastp = i->second;
  }
  halfxs += (nWir_ - lastw)*lastp; halfxs = halfxs/2.;

  double dist = sqrt((fabs(x[0])-(halfxs+nWir*wirP/2.))
		   *(fabs(x[0])-(halfxs+nWir*wirP/2.))
		   + x[1]*x[1] + x[2]*x[2]);  
  if (fabs(dist)>tolerance)     // ***** TEST FAILS => ERROR (non fatal) MESSAGE
    CsErrLog::msg(elError,__FILE__,__LINE__,
  "Continuity test failed for %s, check \"detectors.dat\"",GetTBName().c_str());

  //      ***** STORING COORDINATES OF THE SUBDETECTOR CENTER *****
  xcms_[id]=xcm; ycms_[id]=ycm; zcms_[id]=zcm;

  if(x[0]>0) {
    wirPs_[nWir_] = wirP;
  }
  else {
    wirD_ -= ((nWir-0.5)*wirP+wirPs_[0]/2.);

    map<int, double> newpitch;
    newpitch[0] = wirP;
    
    for(IP i=wirPs_.begin(); i!=wirPs_.end(); i++) {
      newpitch[i->first + nWir] = i->second;
    }
    
    wirPs_ = newpitch;
  }

  // Update "wirP_" (understood as the average pitch of a variable pitch det)
  wirP_ = (nWir_*wirP_+nWir*wirP)/(nWir_+nWir);

  nWir_ += nWir;  // Update # of wires

  double cmWRS[3];
  rotateVectMRS2WRS(xcm_,ycm_,zcm_,cmWRS[0],cmWRS[1],cmWRS[2]);
  
  // new half size of the active zone
  halfxs += nWir*wirP/2.;

  double deltacm[3];
  rotateVectWRS2MRS(halfxs + wirD_ - wirPs_[0]/2. - cmWRS[0], 0, 0,
		    deltacm[0], deltacm[1], deltacm[2]);
  
  xcm_+=deltacm[0];
  ycm_+=deltacm[1];
  zcm_+=deltacm[2];

  // Output to CsErrLog
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"%s New center position %f %f %f\n         New xyz sizes %f %f %f",
		GetTBName().c_str(),xcm_,ycm_,zcm_,xsiz_,ysiz_,zsiz_);

  // storing this piece's detid
  ids_.insert(id);

  hasVarP_ = true;
  
}


void CsDetector::SetVarPitch( const std::map<int,double> *varP ) {

  hasVarP_ = true;

  std::map<int,double> varPP = *varP;

  int    wire1 = varPP.begin()->first;
  double corr1 = varPP.begin()->second;
  int    wire2 = 0;
  double corr2 = 0.0;

  wirD_ = wirD_ - corr1 * wirP_;    // put the first wire position really to zero!

  for (std::map<int,double>::iterator ipitch = varPP.begin(); ipitch != varPP.end(); ipitch++ ) {
    wire2 = ipitch->first;
    corr2 = ipitch->second;
    if ( wire1 == wire2 && corr1 == corr2 ) continue;
//     wirPs_[wire1] = wirP_ * ( 1.0 - corr2 + corr1 );

    double diff = corr1 - corr2;
//     if ( diff < -1.0 ) {
//       cout << "Wire = " << wire1 << "     corr1 = " << corr1 << "     corr2 = " << corr2 << "     diff = " << diff << endl;
// //       diff = -0.99;
//     }
    wirPs_[wire1] = wirP_ * ( 1.0 + diff );

    wire1 = wire2;
    corr1 = corr2;
  }

  //#define CsDetect_DEBUG_VarPitch
#ifdef CsDetect_DEBUG_VarPitch
  // print out ABSOLUTE values of variable-sized pitches
  cout << "CsDetector::VarPitch" ;
  for (std::map<int,double>::iterator ipitch = wirPs_.begin(); ipitch != wirPs_.end(); ipitch++ ) {
    cout << "  " << ipitch->second;
  }
  cout << endl;

  // print out RELATIVE values of variable-sized pitches
  cout << "CsDetector::SetVarPitch" ;
  for (std::map<int,double>::iterator ipitch = *varP.begin(); ipitch !=* varP.end(); ipitch++ ) {
    cout << "  " << ipitch->second;
  }
  cout << endl;
#endif
}


double CsDetector::Wire2Pos(double wire) const {
  // Converts distance to first wire from wire units to mm.
  // (In case of a uniform pitch <P>, it returns: wire*P.)
  map<int,double>::const_iterator iPitch = wirPs_.begin();
  int lastCh = iPitch->first;
  double lastP = iPitch->second, firstP = lastP, dist;
  for (iPitch++, dist = 0; iPitch!=wirPs_.end(); iPitch++) {
    int ch = iPitch->first; if (wire<ch) break;
    dist += (ch-lastCh)*lastP;
    lastP = iPitch->second; lastCh = ch;
  }

  //  There are 2 different ways to envision a sequence of pitcthes:
  //  w       0<=w<ch1       ch1<=w<ch2
  //  d (I)     p0*w   p0*(ch1-1) +p0/2+p1/2 +p1*(w-ch1)
  //  d (II)    p0*w   p0*ch1+p1*(w-ch1)
  // And, propably, and unfortunately, both are at play in coral.      
  //  (I) has been implemented in "Wire2Pos" from its inception (by Colin).
  //     It's been hence used for the alignment of MMs ever since. Therefore,
  //    whether it corresponds to the actual MM scheme or not, we cannot easily
  //    change it.
  //    (Note still that I(Y.B.) slightly modified this ab origo scheme to
  //    correct the case of argument "wire" < 0, where initially "-firstP/2"
  //    was returned.)
  dist += (wire-lastCh+0.5)*lastP - firstP*0.5;
  // (II) seems to be used by Albert for the VarPitch option, cf.
  //     CsDetector.cc,v1.50,  where the above was replaced by:
  //   dist += (wire-lastCh)*lastP;

  return dist;
}

int CsDetector::Pos2Wire(double x) {
  // Converts distance to the first wire (WRS) to wire number
  // (This is the ``inverse'' of "Wire2Pos". It accepts as argument a
  // position w.r.t. 1st wire, not a U coordinate. Therefore it does NOT take
  // into account the distance "wirD_" of the detector's center to 1st wire.)
  double pos = x;

  map<int,double>::const_iterator iPitch = wirPs_.begin();
  int lastCh = iPitch->first, wire;
  double lastP = iPitch->second;
  x += 0.5 * lastP; //x is now referenced to the edge of the active zone
  double lastX = x;

  for (iPitch++, wire = 0; iPitch!=wirPs_.end(); iPitch++) {
    int deltaCh = iPitch->first - lastCh;
    x -= deltaCh*lastP; if (x<0) break;
    wire += deltaCh;
    lastP = iPitch->second; lastCh = iPitch->first;
    lastX = x;
  }
  
  if (lastX<0.) { 
    CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
		  "%s: Hit outside range %f",GetTBName().c_str(),pos);
    wire = -1;
  }
  else {
    wire += int(lastX/lastP);
    if (wire>=nWir_) { 
      CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
		    "%s: Hit outside range %f",GetTBName().c_str(),pos);
    }
  }
  return wire;
}

double CsDetector::Pitch(double wire) {
  double pitch = 0; map<int,double>::iterator ipitch;
  for (ipitch = wirPs_.begin(), pitch = 0; ipitch!=wirPs_.end();
       ipitch++ ) {
    if (wire>=ipitch->first-ipitch->second/2) pitch = ipitch->second;
  }
  if (wire>=nWir_+pitch/2) pitch = 0;
  return pitch;
}

bool CsDetector::IsMyHit(int hitid) {
  
  set<int>::iterator subid = ids_.find(hitid);
  if( subid != ids_.end()) return true;
  else return false;
}

struct CsDetector::sortZones_ : 
  public binary_function<CsZone*, CsZone*, bool> {
  bool operator() ( CsZone* z1, CsZone* z2 ) {
    if( (*z1) < (*z2) ) {
      return( true );
    }
    else {
      return( false );
      }
  }
};

void CsDetector::addZone( CsZone& zone ) {
  myZones_.push_back( &zone );
  myZones_.sort( CsDetector::sortZones_() );
}

void CsDetector::addMCHit( CsMCHit& hit ) {
  myMCHits_.push_back( &hit );
}

void CsDetector::clearMCHitList() {
  myMCHits_.clear();
  myDigits_.clear();
  decodingDone_ = false;
  decode_ = false;
}

void CsDetector::rotateVectDRS2MRS
        (const double u, const double v, const double w, 
        double & x, double & y, double & z ) const {
    HepMatrix rotM(3,3); 
    rotM = getRotDRS();        
    x = rotM(1,1)*u+rotM(1,2)*v+rotM(1,3)*w;
    y = rotM(2,1)*u+rotM(2,2)*v+rotM(2,3)*w;
    z = rotM(3,1)*u+rotM(3,2)*v+rotM(3,3)*w;
}

void CsDetector::rotateVectMRS2DRS
        (const double x, const double y, const double z, 
        double & u, double & v, double & w )  const {
//  HepMatrix rotM(3,3) = getRotDRS();        
  HepMatrix rotM(3,3);
  rotM = getRotDRS();        
  int err;
  HepMatrix irotM(3,3);
  irotM = rotM.inverse( err );
  u = irotM(1,1)*x+irotM(1,2)*y+irotM(1,3)*z;
  v = irotM(2,1)*x+irotM(2,2)*y+irotM(2,3)*z;
  w = irotM(3,1)*x+irotM(3,2)*y+irotM(3,3)*z;
}


void CsDetector::rotatePointMRS2DRSOppanCOMGEANTStyle(const double x, const double y, const double z, double& u, double& v, double& w,int hitID)
{  
	// Find center of this hit's subdetector
	double xcm=xcm_, ycm=ycm_, zcm=zcm_;
	map<int, double>::iterator i;
	if((i=xcms_.find(hitID))!=xcms_.end())
	  xcm = i->second;
	if((i=ycms_.find(hitID))!=ycms_.end())
	  ycm = i->second;
	if((i=zcms_.find(hitID))!=zcms_.end())
	  zcm = i->second;
	double xx = x - xcm;
	double yy = y - ycm;
	double zz = z - zcm;
    rotateVectMRS2DRS(xx,yy,zz,u,v,w);
}


void CsDetector::rotatePointDRS2MRS
        (const double u, const double v, const double w, 
        double & x, double & y, double & z ) const {
  rotateVectDRS2MRS(u,v,w,x,y,z);
  x=x+getXcm(); // perform translation
  y=y+getYcm();
  z=z+getZcm();
}

void CsDetector::rotatePointMRS2DRS
        (const double x, const double y, const double z, 
        double & u, double & v, double & w ) const {
  double xx=x-getXcm(); // perform translation
  double yy=y-getYcm();
  double zz=z-getZcm();
  rotateVectMRS2DRS(xx,yy,zz,u,v,w);
}

void CsDetector::rotateVectWRS2MRS
        (const double u, const double v, const double w, 
        double & x, double & y, double & z ) const {
    x = rotWRS_(1,1)*u+rotWRS_(1,2)*v+rotWRS_(1,3)*w; 
    y = rotWRS_(2,1)*u+rotWRS_(2,2)*v+rotWRS_(2,3)*w;
    z = rotWRS_(3,1)*u+rotWRS_(3,2)*v+rotWRS_(3,3)*w;
}

void CsDetector::rotateVectMRS2WRS
        (const double x, const double y, const double z, 
        double & u, double & v, double & w )  const {
  u = rotWRS_inv_(1,1)*x+rotWRS_inv_(1,2)*y+rotWRS_inv_(1,3)*z; 
  v = rotWRS_inv_(2,1)*x+rotWRS_inv_(2,2)*y+rotWRS_inv_(2,3)*z;
  w = rotWRS_inv_(3,1)*x+rotWRS_inv_(3,2)*y+rotWRS_inv_(3,3)*z;
}

void CsDetector::rotatePointWRS2MRS
        (const double u, const double v, const double w, 
        double & x, double & y, double & z ) const {
  rotateVectWRS2MRS(u,v,w,x,y,z);
  x=x+getXcm(); // perform translation
  y=y+getYcm();
  z=z+getZcm();
}

void CsDetector::rotatePointMRS2WRS
        (const double x, const double y, const double z, 
        double & u, double & v, double & w ) const {
  double xx=x-getXcm(); // perform translation
  double yy=y-getYcm();
  double zz=z-getZcm();
  rotateVectMRS2WRS(xx,yy,zz,u,v,w);
}

void CsDetector::rotateVectWRS2DRS(const double x, const double y, const
				   double z, double & u, double & v, double & w ) const {
  const HepMatrix& rotM=getRotDRS2WRS();
  u = rotM(1,1)*x+rotM(1,2)*y+rotM(1,3)*z; 
  v = rotM(2,1)*x+rotM(2,2)*y+rotM(2,3)*z;
  w = rotM(3,1)*x+rotM(3,2)*y+rotM(3,3)*z;  
}

struct CsDetector::sortClusters_ :
	public binary_function<CsCluster*, CsCluster*, bool> {
	bool operator() ( CsCluster* c1, CsCluster* c2 ) {
		if( c1->getU() < c2->getU() ) return true;
		else return false;
	}
};

void CsDetector::sortClusters( void ) 
{
	myClusters_.sort( sortClusters_() );
}

bool CsDetector::inActiveArea(double x, double y) const {

  // Returns true if argument MRS(x,y) is in active aread, i.e. w/in frame AND
  // outside dead zone. (Nota bene: this routine is slow! But fairly general.
  // A faster calculation is provided by TraFFiC's TDetect::InFrame, which works
  // for the restricted cases where:
  // - dead zone center coincides w/ detector's,
  // - either circular shape or, rectangular aligned with detector's sides.

  double xdrs, ydrs, z;
  // rotate argunent (x,y) to frame system (don't use "rotatePointMRS2DRS"
  // 'cause we don't want to bother about z)
  x -= getXcm(); y -= getYcm(); rotateVectMRS2DRS(x,y,0,xdrs,ydrs,z);
 

  if (fabs(xdrs)<getXsiz()/2 && fabs(ydrs)<getYsiz()/2) {
    // W/in frame => check dead zone...
    if (outDeadZone(xdrs,ydrs)) return true;
    else                        return false;
  }
  else                          return false;
}

bool CsDetector::outDeadZone(double xdrs, double ydrs) const {

  // This routine is made standalone so that TraFFiC's TDetect::InFrame
  // defaults to it when confronted to a case where none of its
  // simplications applies.

  double xdzrs, ydzrs, zdzrs;
  // translate/rotate to DZRS
  xdrs -= getDZXdrs(); ydrs -= getDZYdrs();
  HepMatrix rotM(3,3); rotM = getRotD2DZ();
  xdzrs = rotM(1,1)*xdrs+rotM(1,2)*ydrs;

  if      ((getDZType()&0xf)==0x1) {                 // ...rectangular dead zone
    ydzrs = rotM(2,1)*xdrs+rotM(2,2)*ydrs;
    if (fabs(xdzrs)>getDZDim() || fabs(ydzrs)>getDZYdim()) return true;
    else                                                   return false;
  }
  else if ((getDZType()&0xf)==0x5) {                 // ...circular dead zone
    zdzrs = rotM(3,1)*xdrs+rotM(3,2)*ydrs;
    if (xdzrs*xdzrs+zdzrs*zdzrs>getDZDim())                return true;
    else                                                   return false;
  }
  else {
    // This is expected not to ever happen (cf. assignments @ initialization
    // and when reading from file). But as long as this is ``slow'' routine,
    // make it safe and check.
    CsErrLog::mes(elFatal,"Unexpected Dead Zone shape!\n");
    return false;
  }
}


struct CsDetector::sortMyDigits_ : 
  public binary_function<CsDigit*, CsDigit*, bool> {
  bool operator() (CsDigit* d1, CsDigit* d2) { 

    int addr1 = d1->getAddress();
    int addr2 = d2->getAddress();

    if( addr1 < addr2 ) {
      return( true );
    }
    else {
      return( false );
    }
  }
};

void CsDetector::sortMyDigits( void ) { 
  if( myDigits_.size() > 1 ) {
    myDigits_.sort( sortMyDigits_() ); 
  }
}

void CsDetector::restoreU(CsCluster *left, CsCluster *right) {
  // The CsClusters of drift detectors may be modified by the tracking software
  // (propagation time correction, coalescence). This method allows to restore
  // the original U. It acts on the 2 mirrors simultaneoulsy. Left (i.e.
  // Salève or bottom) is expected to supplied first and will be assigned
  // a <0 drift, so that the ordering of CsClusters w/ increasing U is
  // preserved.
  if (left->getMirrorCluster()!=right)
    CsErrLog::mes(elFatal,"Arg CsClusters are not mirrors of each other!");
  const std::list<CsDigit*> &digits = left->getDigitsList();
  if (digits.size()!=1)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "Arg CsCluster w/ digits list of size %d != 1",digits.size());
  int wire = digits.front()->getAddress();
  int err; HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
  double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 
  double u = wireDCorr + Wire2Pos(wire);
  if (!left->hasTime()) CsErrLog::mes(elFatal,"Arg CsCluster w/o time info");
  double time; left->getTime(time);
  bool error; double drift = getDistToWire(time,error);
  left->setU(u-drift); right->setU(u+drift);
}
