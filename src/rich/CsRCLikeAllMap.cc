/*!
   \File    CsRCLikeAllMap.cc
   \----------------------------
   \brief   CsRCLikeAllMap class implementation.
   \author  Stefano Panebianco
   \version 1.0
   \date    19 March 2004
*/


//- WARNING: OBSOLETE!


  #include <iostream>
  #include <ostream>
  #include <cstdio>

  #include "CsOpt.h"
  #include "CsHist.h"
  #include "CsHistograms.h"

//------------------------------
  #include "CsRCLikeAll.h"
  #include "CsRCLikeAllMap.h"

  #include "CsRCEventPartPhotons.h"
  #include "CsRCPartPhotons.h"
  #include "CsRCPhoton.h"
  #include "CsRCParticle.h"
  #include "CsRCCluster.h"

  #include "CsRCHistos.h"

  #include "CsRCRecConst.h"
  #include "CsRCExeKeys.h"
  #include "CsRCMirrors.h"

  #include "TLorentzVector.h"
  #include <fstream>

//------------------------------

  using namespace std;

  CsRCLikeAllMap::CsRCLikeAllMap( const CsRCPartPhotons* papho ):
    CsRCLikeAll() {

    pPartPhot_ = papho;

  }

  CsRCLikeAllMap::~CsRCLikeAllMap() {}


  extern "C" {
    float prob_( const float&, const int& );
    float erf_( const float& );
  }



//===========================================================================
  double CsRCLikeAllMap::normSignal( const double theIpo , double N0) {
//---------------------------------------------------------


//- compute Signal normalization for Likelihood
//  -------------------------------------------
//- 

    CsRCRecConst* cons = CsRCRecConst::Instance();

    N0=2.2;	// when beta=1, (n_ph/eff)=N0 L sin^2(55 mrad)=0.90 N0=20
    //const float L = cons->Lenght();        // paolo 050204
    const float L = pPartPhot_->pPart()->pMirPart()->RR()/2.;
    const float eff = 0.73;
    
    double normS = eff*N0*L*sin(theIpo/1000.)*sin(theIpo/1000.);

//	- N0 should depend on the detector
// 	- this gives n_j but we need also s_j<=n_j  which takes into account that 
//	may be only part of the ring is in the detector
//	(in particular when a ring is split between up and down detectors).

    if( !(theIpo > 0.) ) normS = 0.;
    
    return  normS;

  }


//===========================================================================
  double CsRCLikeAllMap::normBackgr(std::list<CsRCPhoton*, std::allocator<CsRCPhoton*> >& lPhotons) {
//---------------------------------------------------------


    CsRCRecConst* cons = CsRCRecConst::Instance();
    double normB = 0.;
    list<CsRCPhoton*>::iterator ih;
    
    
    
    for(ih=lPhotons.begin(); ih!=lPhotons.end(); ih++){
      
      double xPhot = (*ih)->ptToClu()->xc();
      double yPhot = (*ih)->ptToClu()->yc();
           
      normB += likeBackgr(xPhot, yPhot);
   	     
    }
    
//----------------------------------------------------------------------
    if( normB < 0. ) normB = 0.;

    return  normB;

  }


//===========================================================================
  double CsRCLikeAllMap::likeSignal( const double thePhot,
//-------------------------------------------------------
				    const double theIpo,
				    const double sigPhot ) {


    CsRCRecConst* cons = CsRCRecConst::Instance();
    double TwoPI = cons->TwoPI();
    double sq2pi = sqrt( TwoPI );
    //float F = cons->Focal();        // paolo 050204
    float F = pPartPhot_->pPart()->pathLen();
    
    double qSw = (thePhot - theIpo) / sigPhot;
    double signal = exp( - 0.5* qSw*qSw ) / (sigPhot*F*sq2pi/1000);
    if( !(theIpo > 0.) ) signal = 0.;

    return  signal;

  }


//===========================================================================
  double CsRCLikeAllMap::likeBackgr( const double xPhot,
//-------------------------------------------------------
				    const double yPhot ) {


    CsRCRecConst* cons = CsRCRecConst::Instance();
    CsRCEventPartPhotons* evpphot = CsRCEventPartPhotons::Instance();
    
    int i = int((xPhot+1600.)/8 +1) ;
    int j = int((yPhot+1600.)/8 +1) ;
    //double backgr = (evpphot->backgrmapval(i,j)/(8*8))/10.;
    double backgr = (cons->backgrmapval(i,j)/(8*8))/10.;       //   ???
    

    return  backgr;

  }

//===========================================================================
double CsRCLikeAllMap::likeSignalCorr( double x0,  double y0,  double tx0,  double ty0, 
 double theta,  double phi, double xph, double yph, const int top) {

// For a given track and given photon angles (theta, phi)
// returns the fraction of the photons which will go
// in the top or botton detector.
// Inputs are: distances in mm, theta in mrad, phi in degrees
// top = 1 ( in upper detector ),  top = 0 ( in lower detector ), top = -1 (no up-down correction, for normalis.)
  
  CsRCRecConst* cons = CsRCRecConst::Instance();
  
  double R=50.;
  //const float L = cons->Lenght();        // paolo 050204
  const float L = pPartPhot_->pPart()->pMirPart()->RR()/2.;
  double TwoPI = cons->TwoPI();
  double PI = TwoPI/2.;
  phi*=PI/180.;

  theta/=1000.;
  
//  reference frame with u along the track and v horizontal pointing to +   
  TVector3 u,v,w;
  u.SetXYZ(tx0,ty0,1.);
  u.SetMag(1.);
  v.SetXYZ(1.,0.,-tx0);
  v.SetMag(1.);
  w.SetXYZ(-tx0*ty0,1.+tx0*tx0,-ty0);
  w.SetMag(1.);

// 	k = unit vector in the direction of the emitted photon 
//	(phi=0 corresponds to v)
  if (cos(theta)==0) {
    printf( "photon emitted at 90 degrees !!!!!");
    return -1;
  }
  TVector3 k=u+ tan(theta) * (cos(phi)*v+ sin(phi)*w);
  k.SetMag(1.);

  double sum=0.;
  double ff=0.;
  int i, N=50;
  for (i=0;i<N;i++){
  
    double z=L*(i+0.5)/N;
// 	photon emitted in M(x0+tx0*z, y0+ty0*z, z)
// 	parametrize the photon line as MP= lambda k
// 	P(x0+tx0*z+lambda*k.x, y0+ty0*z+lambda*k.y, z+lambda*k.z)
// 	d^2(P,ox)= (x0+tx0*z+lambda*k.x)^2 + (y0+ty0*z+lambda*k.y)^2
// 	and find lambda which minimizes d^2
// 	partial d^2/partial lambda = 
//	    2*(x0+tx0*z+lambda*k.x)*k.x + 2*(y0+ty0*z+lambda*k.y)*k.y=0
    double lambda_m = -(x0+tx0*z)*k.X()-(y0+ty0*z)*k.Y();
    if (k.X()*k.X() +k.Y()*k.Y()==0) lambda_m=0;	
  		// k//z d independent from lambda
    else lambda_m /= k.X()*k.X() +k.Y()*k.Y();

// 	now, 0<M.z<L so 0<lambda< (L-z)/k.z, 
// 	and we define lambda_r which minimizes d^2 in this range for lambda
    double lambda_r=lambda_m;
    if (lambda_m <0) lambda_r=0;
    if (k.Z()==0) {
      printf( "photon direction perpendicular to z !!!!!");
      return -1;
    }
    if (lambda_m >(L-z)/k.Z()) lambda_r=(L-z)/k.Z();

// 	if d(lambda_r)<R, the photon crosses the pipe and is killed
    double d2=pow(x0+tx0*z+lambda_r*k.X(),2) + pow(y0+ty0*z+lambda_r*k.y(),2);
    double f;
    if (d2 < R*R) f=0;
    else f=1;   
  
// 	check top or bottom detector, corresponds to z+lambda k.z=L
    double y_mirror=y0+ty0*z +(L-z)*k.y()/k.z();
    if ((y_mirror<0)&&(top==1)) f=0;
    if ((y_mirror>0)&&(top==0)) f=0;
    
    if(i==0){

   // 	check if inside lunette
   //    double x_mirror=x0+tx0*z +(L-z)*k.x()/k.z();
       double ThetaPhi2XY = thetaphitoxy( x0,  y0,  tx0,  ty0,  z,  theta,  phi,  xph,  yph);
       double x=fabs(xph);		//check
       double y=fabs(yph)-1535;   	// check 1535 and for lower detector 
       if (x>201) ff=1;
       if (y>179.132-0.0945265*x+4.67622e-04*x*x-1.69670e-05*x*x*x) ff=1;
       if ((x>180)&&(y>12+14.2568*sqrt(201-x))) ff=1;
       if (y<12-11.2031*sqrt(201-x)+2.67675*x-1.33172e-02*x*x) ff=1;
       if ((x>18)&&(y<110.234+7.45164e-03*x-2.56130e-03*x*x-2.93800e-06*x*x*x)) ff=1;
       if ((x<18)&&(y<45)) ff=1.;  

       // check for the strips of dead zones.
   // horizontal strips (2) from up to down
       if (y>848.220 && y < 871.872 ) ff = 0.;
       if (y>-871.082 && y<-840.429) ff = 0.;
   // vertical strips (3 long ones + 1 short one) from left to right
       if (x>-663.758 && x<-608.014) ff = 0.;
       if (x>-31.468 && x<31.748) ff = 0.;
       if (x>608.288 && x<663.832) ff = 0.;
       if (x>711.938 && x<759.491 && y>848.220) ff = 0.;

    }
    
    f=f*ff;

//  cout <<"u= "<< u.x() <<" "<< u.y() <<" "<< u.z() << "\n";
//  cout <<"v= "<< v.x() <<" "<< v.y() <<" "<< v.z() << "\n";
//  cout <<"w= "<< w.x() <<" "<< w.y() <<" "<< w.z() << "\n";
//  cout <<"k= "<< k.x() <<" "<< k.y() <<" "<< k.z() << "\n";
//  cout <<"lambda_m= "<< lambda_m << "\n";
//  cout <<"d= "<< sqrt(d2) << "\n";
//  cout <<"y_mirror= "<< y_mirror <<" "<<y0<<" "<<k.y()/k.z()<<" "<<z<< "\n";
//  if (top) cout << "in the top detector \n";
//  else cout << "in the bottom detector \n";
//  cout <<"f= "<< f << "\n";
    sum+=f;
  }
  sum=sum/N;
  //cout <<"sum= "<< sum << "\n";
  return sum;

}
   
//===========================================================================
double CsRCLikeAllMap::normSignalCorr(const double x0, const double y0, const double tx0, const double ty0, 
const double theta, const double xph, const double yph) {

CsRCRecConst* cons = CsRCRecConst::Instance();
const double twopi = cons->TwoPI();
const int top = -1;
const double bin = 3.6;
const double nbins = 360./bin;

double S = 0;
	
for(double i=0; i<360.; i+=bin){
	
	S += likeSignalCorr(x0, y0, tx0, ty0, theta, i, xph, yph, top);

}

double f = S/nbins;

  return f;
}
   
//===========================================================================
  double CsRCLikeAllMap::thetaphitoxy(const double x0, const double y0, const double tx0, const double ty0, 
  double z, double theta, double phi, double &xph, double &yph) {
  
//-----------------------------------------
// input: tracks defined by x0, y0, tx0, ty0 at RICH entrance (z0)
//		which emits in z a photon defined by theta, phi
// output: x and y of the reflected photon on detector plane

//  double xC=0.,yC=1600.,zC=2385.;	// coordinate of mirror center 

 
  CsRCRecConst* cons = CsRCRecConst::Instance();
    
  //const float L = cons->Lenght();        // paolo 050204
  const float L = pPartPhot_->pPart()->pMirPart()->RR()/2.;
  double TwoPI = cons->TwoPI();
  
  
  TVector3 C(0.,1600.,2385.);	//  mirror center
  double R=6600;	// mirror radius
  double z0=6000; 	// RICH entrance at which x0 and y0 are defined
  double zm=C.z()+sqrt(R*R-C.y()*C.y());	// z of mirror on z axis
//  double L=2798.;	//8798-z0
//  double yD=1931, zD=5941;	// detector center
  TVector3 D(0.,1931.,5941.);	//  detector center
  double sinThetaD=0.14943;	// detector angle
  double cosThetaD=sqrt(1-sinThetaD*sinThetaD);		// 0.97767
  TVector3 n(0,-sinThetaD,cosThetaD);	// normal to detector plane

//----------	compute k = unit vector in the direction of the emitted photon
//  reference frame with u along the track and v horizontal pointing to +   
  TVector3 u(tx0,ty0,1.);
  u.SetMag(1.);
  TVector3 v(1.,0.,-tx0);
  v.SetMag(1.);
  TVector3 w(-tx0*ty0,1.+tx0*tx0,-ty0);
  w.SetMag(1.);
  if (cos(theta)==0) {
    cout << "photon emitted at 90 degrees !!!!!";
    return -1;
  }
//	(phi=0 corresponds to v)
  TVector3 k= u +tan(theta) * (cos(phi)*v +sin(phi)*w);
  k.SetMag(1.);
  //cout <<"k: "<<k.x()<<" "<<k.y()<<" "<<k.z()<<endl;
  
//---------- check top or bottom mirror 
  TVector3 P(x0+tx0*(z-z0), y0+ty0*(z-z0), z);	//  point of photon emission
//  if (P.y()+(k.y()/k.z())*(zm-z)<0){	// commented for test
  if (yph<0) {				// added for test
    C.SetY(-C.y());	// use bottom mirror
    D.SetY(-D.y());	// use bottom detector
    n.SetY(-n.y());	// normal to det plane
    //cout << "bottom detector" << endl;
  }

//---------- compute intersection M with mirror sphere
//	photon emitted in P, let's find intersection M of photon line with mirror
//	PM= lambda k => CM= CP +lambda k 
//	CM^2=R^2 => lambda^2 +2k.CP lambda +CP^2-R^2=0
//	delta'=(k.CP)^2-CP^2+R^2
//	=> lambda = -k.CP + sqrt(delta')	(positive solution)
  TVector3 CP=P-C;
  double deltap=pow(k.Dot(CP),2) -CP.Mag2() +R*R;
  double lambda = -k.Dot(CP) + sqrt(deltap);	
  //cout <<"lambda "<<lambda<<" "<<zm-z<<endl;	// check lambda \approx zm-z
  TVector3 M=P+lambda*k;
//---------- compute unit vector l for reflected photon 
// 	k_paral = (k.CM/CM^2) CM
//	l=k_perp-k_paral=k-2k_paral
  TVector3 CM=M-C;
  TVector3 kpar=((k.Dot(CM))/CM.Mag2())*CM;
  TVector3 l=k-2*kpar;	
  //cout <<"l: "<<l.x()<<" "<<l.y()<<" "<<l.z()<<endl;
//---------- compute intersection I with detector plane 
//	DI.n=0 	=> MI.n=-DM.n
//	and since MI=mu*l	=> mu=-DM.n/l.n
  double mu=-n.Dot(M-D)/n.Dot(l);
  //cout <<"mu "<<mu<<endl;
  //cout <<"M: "<<M.x()<<" "<<M.y()<<" "<<M.z()<<endl;
  TVector3 I=M+mu*l;
  xph=I.x();
  yph=(I.y()-D.y())*cosThetaD;  // coordinates wrt detector center
  if (I.y()>0) yph+=860;	// shift to get Paolo's coordinate
  if (I.y()<0) yph-=860;	// 860 has to be checked
//  if (I.y()>0) yph=I.y()-1070;
//  if (I.y()<0) yph=I.y()+1070;
//  cout<<"xph="<<xph<<"  yph="<<yph<<endl;
//	
//	I=M+mu*l   in detector plane: 
//	sin(theta)*(M.y+mu*l.y-yD) -cos(theta)(M.z+mu*l.z-zD)=0
//	mu = sin(theta)*yD-cos(theta)*zD
  return 1;
  
}


