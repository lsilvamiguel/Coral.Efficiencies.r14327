#ifndef __TRACK_RECONSTRUCTOR_H__
#define __TRACK_RECONSTRUCTOR_H__

#include <stdio.h>
#include <math.h>
#include <vector>
#include <stdlib.h>

const int X_PLANE = 0;
const int Y_PLANE = 1;
const int U_PLANE = 2;
const int V_PLANE = 3;

class CIntersection;
class CIntersectionArray;
class CPlane;
class CTrPointArray;

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

inline float MAX(float a, float b)
{
  return (a>b)?a:b;
}

inline float MAX(float a, float b, float c)
{
  return MAX(MAX(a,b),c);
}

inline float POW2(float a)
{
  return a*a;
}

struct PlaneGeometry
{
  int   planeType;
  int   nWires;
  float pitch;
  float slope;
  float width;
  float height;
};

extern PlaneGeometry XPlane;
extern PlaneGeometry UPlane;
extern PlaneGeometry VPlane;
extern PlaneGeometry YPlane;


const int MAX_POINTS_NUM = 20;

struct MwpcEvent
{
  int  fNhits[4];
  int* fHits[4];
  int  fNpoints;
  float fX[MAX_POINTS_NUM];
  float fY[MAX_POINTS_NUM];
  float fErr[MAX_POINTS_NUM];
};

class CObject
{
 public:
  CObject() {};
  virtual ~CObject() {};
};

class CAbstractArray
{
 protected:
  int size;
//  CObject** array;
  std::vector<CObject*> array;
  bool isEmpty;
 public:
//   CAbstractArray() { array = NULL; size = 0; isEmpty = true; }
  CAbstractArray() { array.clear(); size = 0; isEmpty = true; }
  virtual ~CAbstractArray() { Destroy(); }
  CObject* operator[](int n) { return GetObject(n); }
  void  AddObject(CObject* object);
  CObject* GetObject(int idx) { return ((idx >= size)?NULL:array[idx]); }
  void  Destroy();
  void  ClearArray();
  int   GetSize() { return size; }
  void  RemoveObject(int n);
  bool  IsEmpty() const { return isEmpty; }
};

class CWireArray;

class CWire : public CObject
{
 public:
  static float   xcenter;
  static float   ycenter;
  static float   ucenter;
  static float   vcenter;
  static float   uvclaim;
 protected:
  float          number;
  float          nTrig;
  float          nXtalk;
  float          nNoise;
  float          nEff;
  float          nIneff;
 public:
  int            plane;
  CWire() {};
  virtual ~CWire() {};
  void SetWire(int planeType, int _number);
  CIntersection* GetIntersection(CWire* wire);
  float GetNumber() const { return number; }
  int   GetPlane() const { return plane; }
  void  ClearStatistic();
  bool  HaveIntersection(CWire* wire);
  virtual void TriggeredEvent()    { nTrig++; }
  virtual void EffectiveEvent(float weight)    { nEff += weight; }
  virtual void IneffectiveEvent(float weight)  { nIneff += weight; }
  virtual void NoiseEvent(float weight)        { nNoise += weight; }
  virtual void XtalkEvent(float weight)        { nXtalk += weight; }
  float GetNTrig()    const  { return nTrig; }
  float GetNEff()     const  { return nEff; }
  float GetNIneff()   const  { return nIneff; }
  float GetNNoise()   const  { return nNoise; }
  bool  IsTriggered() const  { return (nTrig != 0.); }
  bool  IsMultipleIntersected(CWireArray* array, int nIntersections = 2);
};

class CTrPoint : public CObject
{
 public:
  float x;
  float y;
  float error;
  float test;
 public:
  CTrPoint() : x(0.), y(0.), error(0.0) {};
  CTrPoint(CTrPoint& p) : x(p.x), y(p.y) {};
  CTrPoint(float _x, float _y);
  virtual ~CTrPoint() {}
  void  SetX(float _x) { x = _x; }
  void  SetY(float _y) { y = _y; }
  float GetDistance(CTrPoint* p) 
    { return sqrt(POW2(x-p->x)+POW2(y-p->y)); }
  float GetDistance(CWire* w);
  virtual float GetX() { return x; }
  virtual float GetY() { return y; }
  float GetError() const { return error; }
  void Print(FILE* outstream) 
    { fprintf(outstream,"X = %.1f; Y = %.1f, err = %.1f\n",x,y,error); }
};

class CTrPointArray : public CAbstractArray
{
 public:
  CTrPointArray() {};
  ~CTrPointArray() { Destroy(); }
  CTrPoint* operator[](int idx) { return (CTrPoint*)GetObject(idx); }
  void Print(FILE* outstream);
};

class CIntersection : public CTrPoint
{
 public:
  CWire* wire[2];
 public:
  CIntersection() { wire[0] = wire[1] = NULL; }
  virtual ~CIntersection() {};
  bool IsEqual(CIntersection* obj)
    {
      return ((obj)&&((obj->wire[0] == wire[0] && obj->wire[1] == wire[1]) ||
	(obj->wire[0] == wire[1] && obj->wire[1] == wire[0])));
    }
  virtual float GetX();
  virtual float GetY();
  virtual float GetU();
  virtual float GetV();
};
  
class CMwpcCluster : public CWire
{
  CWire** cluster;
  int    size;
 public:
  CMwpcCluster(int clusterSize, CWire* hits);
//   ~CMwpcCluster() { delete cluster; }
  ~CMwpcCluster() { free(cluster); }
  float  GetSize() const { return size; }
  CWire* operator[](int n) { return cluster[n]; }
  void  TriggeredEvent();
  void  EffectiveEvent(float weight);
  void  NonEffectiveEvent(float weight);
  void  NoiseEvent(float weight);
};

class CIntersectionArray : public CAbstractArray
{
 public:
  CIntersectionArray() {};
  ~CIntersectionArray() { Destroy(); }
  CIntersection* operator[](int idx) { return (CIntersection*)GetObject(idx); }
  void operator+=(CIntersectionArray* array);
};

class CWireArray : public CAbstractArray
{
 public:
  CWireArray() {};
  virtual ~CWireArray() { ClearArray(); }
  CWire* operator[](int n) { return (CWire*)GetObject(n); }
  CIntersectionArray* GetIntersections(CWireArray* );
  void Print(FILE* outstream);
};

class CMwpcClusterArray : public CWireArray
{
 public:
  CMwpcClusterArray() {};
  virtual ~CMwpcClusterArray() { Destroy(); }
  CMwpcCluster* operator[](int n) { return (CMwpcCluster*)GetObject(n); }
};

class CPlane
{
  friend class CChamber;
  int planeType;
  int nWires;
  int nHits;
  int minNum;
  int maxNum;
  int nReconstrEvents;
 public:
  int*            hits;
  int*            hitTime;
  CWire*          wires;
 public:
  CPlane(PlaneGeometry& geometry);
  virtual ~CPlane();
 public:
  void PutHit(int nWire, int time);
  void ClearHits();
  void ClearStatistic();
  bool UneffectiveCluster(float n);
  CMwpcClusterArray* GetClusters();
  float GetEfficiency(); // Returns mean of summary efficiency for all wires
  float GetNoise();      // Returns sum of noising events for all wires
  int   GetNRecEvents() const { return nReconstrEvents; }
  int   GetHits(int* buffer);
};

class CChamber
{
  MwpcEvent* fReconData;
 public:
  int nReconstrEvents;
  CPlane* planes[4];
  float   efficiency[4];
 public:
  CChamber(PlaneGeometry* config[4]);
  ~CChamber();
  void PutHit(int nPlane, int nWire, int time) 
    { planes[nPlane]->PutHit(nWire,time); }
  CTrPointArray* GetTrackPoints();
  void ClearHits();
  void ClearStatistic();
  bool FindTrackPoint(CMwpcCluster* u, CMwpcCluster* v, CMwpcCluster* x, 
		      CTrPointArray* array);
  bool IsTriplet(CMwpcCluster* u, CMwpcCluster* v, CMwpcCluster* x);
  MwpcEvent* GetReconPoints();
  MwpcEvent* GetHits();
  bool CalculateEfficiency();
  float GetEfficiency(int plane) { return efficiency[plane]; }
};


#endif
