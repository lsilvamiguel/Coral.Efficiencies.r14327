// $Id: CsHelix.cc,v 1.14 2010/10/15 12:02:38 schluter Exp $

/*!
   \file CsHelix.h
   \brief Compass Helix ParametersClass.
   \author  Benigno Gobbo
   \version $Revision: 1.14 $
   \date    $Date: 2010/10/15 12:02:38 $
*/

#include "CsHelix.h"
#include "CsErrLog.h"

using namespace std;
using namespace CLHEP;

CsHelix::CsHelix() : x_(0), y_(0), z_(0), dxdz_(0), dydz_(0), cop_(0)
{
  for(int i=0; i<15; i++ ) cov_[i] = 0;
}

CsHelix::CsHelix( double x, double y, double z, double dxdz, double dydz,
	 double cop, const double* cov ) :
x_(x), y_(y), z_(z), dxdz_(dxdz), dydz_(dydz), cop_(cop)
{
  for(int i=0; i<15; i++ ) cov_[i] = cov[i];
}

CsHelix::CsHelix( const CsHelix& helix ) :
  x_(helix.x_), y_(helix.y_), z_(helix.z_), dxdz_(helix.dxdz_),
  dydz_(helix.dydz_), cop_(helix.cop_)
{
  for(int i=0; i<15; i++ ) cov_[i] = helix.cov_[i];
}

CsHelix& CsHelix::operator=( const CsHelix& helix ) {
  if( this != &helix ) {
    x_    = helix.x_;
    y_    = helix.y_;
    z_    = helix.z_;
    dxdz_ = helix.dxdz_;
    dydz_ = helix.dydz_;
    cop_  = helix.cop_;
    for(int i=0; i<15; i++ ) cov_[i] = helix.cov_[i];
  }
  return( *this );
}

bool CsHelix::operator==( const CsHelix& helix ) const {
  if( x_    == helix.x_    &&
      y_    == helix.y_    &&
      z_    == helix.z_    &&
      dxdz_ == helix.dxdz_ &&
      dydz_ == helix.dydz_ &&
      cop_  == helix.cop_  )
    return( true );
  else
    return( false );
}

HepMatrix CsHelix::getCovMat() const {

  HepMatrix mat(5,5);
  for( int i=0; i<25; i++ ) {
    mat( i/5+1, i%5+1 ) = 
      cov_[ ((i%5)>(i/5))*((i%5)*((i%5)+1)/2 + i/5)+
        ((i%5)<=(i/5))*((i%5)+(i/5)*((i/5)+1)/2) ];
    if( (i%5+1) != (i/5+1) ) mat( i%5+1, i/5+1 ) = mat( i/5+1, i%5+1 );
  }
  return( mat );
}


double& CsHelix::operator () (const int i, const int j)
{
  double *tmp(0);

  if(j==0){
    if((i < 0) || (i > 5))
      cout<<"CsHelix::double& operator()  ==> Index out of range : "
        <<i<<endl;
  }
  else{
    if( (i < 1) || (i > 5) || (j < 1) || (j > 5) )
      cout<<"CsHelix::double& operator()  ==> Indecies out of range : "
        <<i<<" "<<j<<endl;
  }


  if(j==0){
    if( i == 0 ) tmp = &z_; 
    else if( i == 1 ) tmp = &x_;
    else if( i == 2 ) tmp = &y_;
    else if( i == 3 ) tmp = &dxdz_;
    else if( i == 4 ) tmp = &dydz_;
    else if( i == 5 ) tmp = &cop_;
  }
  else {
    tmp = ( j > i ? &cov_[(i-1)+((j-1)*j)/2] : &cov_[(j-1)+((i-1)*i)/2] );
  }

  return *tmp;
}


void CsHelix::Print(const char* str) const
{

  cout<<str<<endl;
  cout.fill(' ');
  cout.setf(ios::showpos);
  cout<<" Z = "<<setprecision(7)<<setw(10)<<
    z_<<" cm"<<endl;
  cout.setf(ios::scientific);
  cout<<setprecision(3)<<setw(10)<<
    x_<<" | "<<cov_[0]<<endl;
  cout<<setprecision(3)<<setw(10)<<
    y_<<" | "<<cov_[1]<<" "<<cov_[2]<<endl;
  cout<<setprecision(3)<<setw(10)<<
    dxdz_<<" | "<<cov_[3]<<" "<<cov_[4]<<" "<<cov_[5]<<endl;
  cout<<setprecision(3)<<setw(10)<<
    dydz_<<" | "<<cov_[6]<<" "<<cov_[7]<<" "<<cov_[8]<<" "<<cov_[9]<<endl;
  cout<<setprecision(3)<<setw(10)<<
    cop_<<" | "<<cov_[10]<<" "<<cov_[11]<<" "<<cov_[12]<<" "<<cov_[13]<<" "<<cov_[14]<<endl;

  //Reset some format flags 
  cout.setf(ios::fixed,ios::scientific);
  cout.unsetf(ios::showpos);
  cout<<setprecision(7);

  if( cop_ != 0 ) cout<<"P = "<<1/cop_<<" +- "<< sqrt(cov_[14]) / (cop_*cop_) <<endl;

}



#include "Traffic.h"
#include "THlx.h"

bool CsHelix::Extrapolate(double Z, CsHelix& hout, bool mmap) const
{

  if(Traffic::Ptr() == NULL) {
    cout<<endl<<"TRAFFIC package will be initialized, as CsHelix::Extrapolate() use it"<<endl;
    new Traffic;
  }
  bool err;
  THlx Hin,Hout; Hout(0)=Z/10.;
  Hin.ImportHelix(*this);
  if(Hin.empty_cov()) CsErrLog::Instance()->mes( elFatal, "Empty covariance matrix");
  err  = Hin. Extrapolate(Hout,mmap); // do the extrapolation
  hout = Hout.ExportHelix();
  return(err);
}


bool CsHelix::FindCDA( CsHelix& H2 )
{
  if(Traffic::Ptr() == NULL) {
    cout<<endl<<"TRAFFIC package will be initialized, as CsHelix::FindCDA() use it"<<endl;
    new Traffic;
  }
  bool err;
  THlx TH1,TH2;
  TH1.ImportHelix( *this );
  TH2.ImportHelix( H2 );
  err  = TH1.FindCDA( TH2 );
  (*this) = TH1.ExportHelix();
  H2      = TH2.ExportHelix();
  return(err);
  
}


