// A set of highly reusable classes inheriting from TObject
// -> can be stored in TCollections.
// CBernet

#ifndef __Structures__
#define __Structures__

#include "TObject.h"
#include "TPolyLine3D.h"
#include <vector>

class Tracker;

class Point3 {
 public:
  float fX;
  float fY;
  float fZ;
  Point3(float x, float y, float z) : fX(x),fY(y),fZ(z) {}
};

class Cluster : public TObject { 
  
 private:
  Tracker* fTracker;
  
 public:
  float fPos;
  float fT;
  float fTot;
  float fSize;
  float fRes;
  //int fZone;         // 0 : inner pitch, 1 : outer pitch
  float fT2;           // second time information (ampl. ratio,APV only)
  float fT_APV_0;      // cluster time from r0
  float fT_APV_1;      // cluster time from r1
  float fT_APV_c;      // cluster time (weighted mean, APV only)
  float fResidual;

  Cluster(float x, float t, float tot, float size, float res,  Tracker* p) 
    : fPos(x),fT(t),fTot(tot),fSize(size),fRes(res)
    {SetTracker(p);}
  
  void SetTracker(Tracker* fTracker_) {fTracker=fTracker_; }
  Tracker* GetTracker() { return fTracker; }
  
  ClassDef(Cluster,0)
    
};

class PointN : public TObject {
  
 private:
  Int_t fDim;
  Float_t *fValues;
  
 public:
  PointN() {fDim=0;}
  PointN(Int_t dim, Float_t* values);
  ~PointN();
  Float_t* GetValues() {return fValues;}
  Float_t operator[] (Int_t i);
  void Dump();
  
  ClassDef(PointN,0)
};
    
    
class Point : public TObject {
  
 private:
  Int_t point[2];
  
 public:
  Point() : TObject() {}
  Point(Int_t x, Int_t y) {
    point[0]=x; point[1]=y;
  }
  
  ~Point(){}
  
  Int_t GetX() {return point[0];}
  Int_t GetY() {return point[1];}
  
  ClassDef(Point,0)
 };
    
	    
 class Digit : public TObject {
   
 private:
   int fChannel;
   float fData;
   
 public:
   Digit(){}
   Digit(int channel, float data):TObject() {
     fChannel=channel;
     fData=data;}
   ~Digit(){}
   int GetChannel() {return fChannel;}
   float GetData() {return fData;}  
   
   ClassDef(Digit,0)
     
};

class Track : public TObject {
  
 private:
  Float_t fX,fY,fZ;     // current reference point
  Float_t fTany,fTanz;  // tan(thetax) and tan(thetay)  
  Float_t fP;           // impulsion (GeV)
  
  TPolyLine3D *fLine;
  
  std::vector<Cluster*> fClusters;
  Float_t chi2;
  Float_t chi2p;

 public:
  
  Track(Float_t x0,Float_t y0,Float_t z0,Float_t tany,Float_t tanz, Float_t p);
  Track(float x1, float y1, float z1, float x2, float y2, float z2, float p);

  Track(Float_t* p1, Float_t* p2, Float_t p);
  Track(Track& track);
  
  void SetClusters (const std::vector<Cluster*>& v) {fClusters = v;}
  std::vector<Cluster*>& GetClusters() {return fClusters;}

  void SetChi2(float chi2_) { chi2 = chi2_; } 
  void SetChi2Prob(float chi2p_) { chi2p = chi2p_; } 
  float GetChi2() {return chi2;}
  float GetChi2Prob() {return chi2p;}
  
  ~Track();
  
  void Move(Float_t x) {
    Float_t dx=x-fX;
    fY+=dx*fTany;
    fZ+=dx*fTanz;
    fX=x;
  }

  Float_t GetX() {return fX;}
  Float_t GetY() {return fY;}
  Float_t GetZ() {return fZ;}
  
  Float_t GetYp() {return fTany;}
  Float_t GetZp() {return fTanz;}
  Float_t GetP() {return fP;}
  
  void Draw(Float_t x1,Float_t x2, Float_t xoffset);
  void Print();

  ClassDef(Track,0)
};


class BridgedTrack : public TObject{

 private:
  Track *fTr1; // track segment before pivotal point
  Track *fTr2; // track segment after pivotal point
  Float_t fXm; // pivotal point
  
 public:
  BridgedTrack(Track* track, Float_t beta0,Float_t p0, Float_t xm);
  virtual ~BridgedTrack();

  Track* GetTrack(Int_t i) 
    {
      if(i==1) return fTr1;
      else if(i==2) return fTr2;
      else return NULL;
    }

  void Draw(Float_t x1,Float_t x2, Float_t xoffset);
  ClassDef(BridgedTrack,0)

};

#endif











