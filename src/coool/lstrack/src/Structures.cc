#include "Structures.h"

#include <iostream>


ClassImp(Cluster);
ClassImp(Point);
ClassImp(PointN);
ClassImp(Track);
ClassImp(BridgedTrack);
ClassImp(Digit);

PointN::PointN(Int_t dim, Float_t* values) {
  
  fDim=dim;
  fValues=new Float_t[fDim];
  for(Int_t i=0;i<fDim;i++) fValues[i]=values[i];
}

PointN::~PointN() {
  delete fValues;
}

Float_t PointN::operator[] (Int_t i) {

  if(i<fDim) return fValues[i];
  else std::cout<<"PointN::operator[] overflow. dim : "<<fDim;
  return 0;
}

void PointN::Dump() {

  for(Int_t i=0;i<fDim;i++) {
    std::cout<<fValues[i]<<"\t";
  }
}

/////////////////////////      Track      //////////////////////////////////


Track::Track(Float_t x0,Float_t y0,Float_t z0,Float_t tany,Float_t tanz, Float_t p) {
  fX=x0;
  fY=y0;
  fZ=z0;
  fTany=tany;
  fTanz=tanz;
  fP=p;
  
  fLine=new TPolyLine3D(2);
}

Track::Track(Float_t* p1, Float_t* p2, Float_t p) {
  fX=p1[0];
  fY=p1[1];
  fZ=p1[2];  
  
  Float_t dx=(p2[0]-p1[0]);
  fTany=(p2[1]-p1[1])/dx;
  fTanz=(p2[2]-p1[2])/dx;

  fP=p;
  fLine=new TPolyLine3D(2);
}

Track::Track(float x1, float y1, float z1, float x2, float y2, float z2, float p) :
  fX(x1), fY(y1), fZ(z1), fP(p) {
  
  float dx=x2-x1;
  fTany=(y2-y1)/dx;
  fTanz=(z2-z1)/dx;
  fLine=new TPolyLine3D(2);  
}

Track::Track(Track& track) {

  fX=track.fX;
  fY=track.fY;
  fZ=track.fZ;
  fTany=track.fTany;
  fTanz=track.fTanz;
  fP=track.fP;
  
  fLine=new TPolyLine3D(2);
}

Track::~Track() {
  delete fLine;
}

void Track::Print() {

  std::cout<<"track : "<<fX<<" "<<fY<<" "<<fZ<<" "<<fTany<<" "<<fTanz<<" "<<fP<<std::endl;
}

void Track::Draw(Float_t x1,Float_t x2, Float_t xoffset) {

  this->Move(x1);
  fLine->SetPoint(0,fX - xoffset,fY,fZ);

  this->Move(x2);
  fLine->SetPoint(1,fX - xoffset,fY,fZ);

  fLine->SetLineColor(5);
  fLine->Draw();
}


BridgedTrack::BridgedTrack(Track* track, Float_t beta0,Float_t p0, Float_t xm) 
  : TObject() {
  fTr1 = new Track(*track);
  
  Float_t zp = fTr1->GetZp();
  Float_t p  = fTr1->GetP();
  Float_t yp = fTr1->GetYp() + beta0 * p0/p;
  fXm  = xm;
  Float_t z  = fTr1->GetZ() + fTr1->GetZp()*(fXm-fTr1->GetX());
  Float_t y  = fTr1->GetY() + fTr1->GetYp()*(fXm-fTr1->GetX());  

  fTr2 = new Track(fXm,y,z,yp,zp,p);
}

BridgedTrack::~BridgedTrack() {

  delete fTr1;
  delete fTr2;
}

void BridgedTrack::Draw(Float_t x1,Float_t x2, Float_t xoffset){

  if(x1 > fXm) fTr2->Draw(x1,x2,xoffset);
  else 
    if(x2 < fXm) fTr1->Draw(x1,x2,xoffset);
    else {
      fTr1->Draw(x1,fXm,xoffset);
      fTr2->Draw(fXm,x2,xoffset);      
    }
}





