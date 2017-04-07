
//File :CsVTrack.cc

#include "coral_config.h"
#include "CsVTrack.h"
#include "TMtx.h"
#include <iostream>

using namespace std;
using namespace CLHEP;

// Constructor
CsVTrack::CsVTrack():
  TrkRef_(0),
  status_(1),
  Chi2Tr_(-1.),
  RadLength_(0),
  ELoss_(0),
  Pk(5),  
  Gk(5,5),   
  
  Ak(5,3),
  AkT(3,5),
  Bk(5,3),
  BkT(3,5),
  Cke(5),
  GkB(5,5),
  Wk(3,3),
  
  Qnk(3),
  Pnk(5),

  Dnk(3,3),
  Enk(3,3)
{
  for(int i=0;i<5;i++) Res_[i]=0;
};


// Copy constructor
CsVTrack::CsVTrack(const CsVTrack& t):

  TrkRef_  (t.TrkRef_),
  status_    (t.status_),
  Chi2Tr_  (t.Chi2Tr_),
  RadLength_(t.RadLength_),
  ELoss_(t.ELoss_),
  
  Pk      (t.Pk),  
  Gk      (t.Gk),   

  Ak     (t.Ak),
  AkT    (t.AkT),
  Bk     (t.Bk),
  BkT    (t.BkT),
  Cke    (t.Cke),
  GkB     (t.GkB),
  Wk      (t.Wk),
  
  Qnk     (t.Qnk),
  Pnk     (t.Pnk),

  Dnk     (t.Dnk),
  Enk     (t.Enk)

{
  for(int i=0;i<5;i++) Res_[i]=t.Res_[i];
}


CsVTrack& CsVTrack::operator = (const CsVTrack& t)
{

  TrkRef_    = t.TrkRef_;
  Chi2Tr_    = t.Chi2Tr_;
  ELoss_     = t.ELoss_;
  RadLength_ = t.RadLength_;
  status_    = t.status_;
  //for(int i=0;i<5;i++) Res_[i]=t.Res_[i];
  Res_[0]=t.Res_[0];
  Res_[1]=t.Res_[1];
  Res_[2]=t.Res_[2];
  Res_[3]=t.Res_[3];
  Res_[4]=t.Res_[4];

  Pk       = t.Pk;  
  Gk       = t.Gk;   
  
  Ak       = t.Ak;
  AkT      = t.AkT;
  Bk       = t.Bk;
  BkT      = t.BkT;
  Cke      = t.Cke;
  GkB      = t.GkB;
  Wk       = t.Wk;
  
  Qnk      = t.Qnk;
  Pnk      = t.Pnk;
  
  Enk      = t.Enk;
  Dnk      = t.Dnk;

  return (*this);
}


double CsVTrack::getMom(){
  if( Qnk(3) != 0 ) return fabs(1/Qnk(3));
  return 9999999;
}

Hep3Vector CsVTrack::Vec () const {
  double cop = fabs(Qnk(3)), p = cop ? 1/cop : 9999999;
  double tgx = Qnk(1), tgy = Qnk(2);
  double pz = p/sqrt(1+tgx*tgx+tgy*tgy);
  return Hep3Vector(tgx*pz,tgy*pz,pz);
}

HepLorentzVector CsVTrack::LzVec (const double mass) const {
  double cop = fabs(Qnk(3)), p = cop ? 1/cop : 9999999;
  double tgx = Qnk(1), tgy = Qnk(2);
  double pz = p/sqrt(1+tgx*tgx+tgy*tgy);
  double e = sqrt(p*p+mass*mass);
  return HepLorentzVector(tgx*pz,tgy*pz,pz,e);
}
