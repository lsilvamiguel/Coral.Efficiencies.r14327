#ifndef __Tracker__
#define __Tracker__

#include<string>
#include<vector>
#include<map>

#include "TObject.h"
#include "TMatrix.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TNode.h"
#include "TFile.h"

#include "Structures.h"

class Tracker : public TObject {

 private:
  
  static const float sqrt12;

  bool fClustering;   //true if clustering if AddCluster is not called

  int fNeff;          //number of good candidates 
  int fNtot;          //total number of track hits
  int fNbackg;        //number of background clusters

  bool fIsActive;     //state
  int fNsigmat;

  float fCheckWidth;  //width of the road for efficiency calculation
  float fCheckMin;  
  float fCheckMax; 

  unsigned int fId;   //id
  std::string fName;       //name

  float fX,fY,fZ;   //chamber center coordinates in cm
                      // x : up->downstream
                      // y : Saleve->Jura 
                      // z : down->up

  float fYbck,fZbck; 

  float fDx, fDy, fDz; // chamber dimensions (DRS)

  float fAngle;     //rotation angle of the wires along x with respect to z axis

  float fHalfsize;
  
  float   fR2Min;      // "software" beam killer

  float fSin;       // sin(fAngle)
  float fCos;       // cos(fAngle)
  float   fXc;        // fX*fCos
  float   fC2;        // fCos*fCos
  float   fXc2;       // fX*fCos*fCos
  float   fX2c2;      // fX*fX*fCos*fCos
  float   fXs;        // fX*fSin
  float   fS2;        // fSin*fSin
  float   fXs2;       // fX*fSin*fSin
  float   fX2s2;      // fX*fX*fSin*fSin
  float   fSc;        // fSin*fCos
  float   fXsc;       // fX*fSin*fCos
  float   fX2sc;      // fX*fX*fSin*fCos
  float   fA;         // fY*fSin - fX*fCos
  float   fAc;        // fA*fCos
  float   fXac;       // fX*fA*fCos
  float   fAs;        // fA*fSin
  float   fXas;       // fX*fA*fSin
  
  int fNwires;        //number of wires
  float fPitch;     //spacing of the wires 
  
  float fIp;          //inner pitch
  float fOp;          //outer pitch
  std::map<int, double> fPitchs;  // all pitchs 

  int   fNwireshit;   //number of hits
  int *fHitPat;   
  int *fHitToTs;
  int *fHitTimes;
  int fPrevChan;

  std::vector<Cluster*> fClusters;

  float fCTMin;
  float fCTMax;
  float fCTotMin;

  TRotMatrix *fRotMatrix; 

 public:
  TH1F *fHres, *fHresl, *fHresz;
  TH2F *fHrespos;
  TProfile *fHresVSpos;
  TH1F *fHctime;   
  TH1F *fHcpos;
  TH1F *fHdcpos;
  TH1F *fHctot;
  TH1F *fHchits;

  TH2F *fHeff;
  TH2F *fHtot;
  TH2F *fHbackg;
  TH2F *fHeffcor;

  TH1F *fH1deff;
  TH1F *fH1dtot;
  TH1F *fH1dbackg;
  TH1F *fH1deffcor;
  TH1F *fHtres;
  
  std::vector<TH1*> fHists;
  
  //length in cm, angle in degrees.
  Tracker(int id, const char *name, Int_t nwires, float x, 
	  float y, float z, float dx, float dy, float dz,
	  float angle, float inpitch, float outpitch);

  // copy constructor
  Tracker(const Tracker& other);
 
  void Init();
  
  // calculates tracking parameters from geometry
  void TrackingParameters();

  ~Tracker();

  // adds another detector piece
  void AddSubDetector(int id, const char *name, Int_t nwires, float x, 
		      float y, float z, float dx, float dy, float dz,
		      float angle, float inpitch, float outpitch); 
  
  

  // finds time cuts
  void TimeCuts(int nsigma, TFile *file);

  // activates the plane
  void Activate(bool active); 

  // returns state of the plane
  bool IsActive() {return fIsActive;}

  // add one digit
  void AddDigit(int chan, int time, int tot);    // wire ranges from 1 to fNwires  

  // add one cluster (pos is given in cm, origin at the detector center)
  void AddClusterWRS(float pos, int size, float res);   
  
  // add one cluster (pos is given in cm, origin at the detector center)
  void AddClusterWRS(float pos, float time, int size, float res);   
  
  // add one cluster (pos is given in channel units)
  void AddCluster(float pos, float time, float tot, float size);   

  // clustering 
  Int_t Clusterize();

  // checks if there is a cluster within checkwidth from (y,z)
  // returns 0 for not found, (i+1) if i-th cluster is in limits
  Int_t CheckCandidate(float y,float z, float checkmin, float checkmax,
		       bool fillHistos=true);
  Int_t CheckCandidate(Track* track, float checkmin, float checkmax,
		       bool fillHistos=true);
  Int_t CheckCandidate(Track* track, bool fillHistos=true);
  
  // checks if the track hits the active zone
  bool InActiveZone(Track* track);
  bool InActiveZone(float y, float z);

  // updates efficiency
  // return value see CheckCandidate
  int Efficiency(Track* track, bool fillHistos=true);
  
  // cluster residual
  float GetResidualOfCluster(float pos, float y, float z);
  
  // residuals for a given track
  void Residuals(float y,float z);
  
  // sets checkrange
  void SetCheckRange(float min, float max) {fCheckMin=min; fCheckMax=max;}

  // set time cuts
  void SetCTimeRange(int min, int max) {fCTMin=min; fCTMax=max;}

  // set tot cuts
  void SetCTotMin(int ctotmin) {fCTotMin = ctotmin;}

  // called before each event
  void ResetHits();
  
  // resets all histograms
  void ResetHistos();

  // print information
  void Dump();

  // draws the plane (frame + clusters)
  void Draw(TNode *worldnode);

  // returns numebr of clusters
  unsigned int GetNclusters() const {return fClusters.size();}

  // returns all clusters
  const std::vector<Cluster*>& GetClusters() const {return fClusters;}

  // returns plane id
  unsigned int     GetId() const {return fId;}

  // returns name
  const char*   GetName() const {return fName.c_str();}

  // accessors to geometrical parameters
  float GetAngle() const {return fAngle;}  //rad
  float GetSin() const  {return fSin;}
  float GetCos() const  {return fCos;}
  float GetX() const  {return fX;}
  float GetY() const  {return fY;}
  float GetZ() const  {return fZ;}
  float GetDX() const  {return fDx;}
  float GetDY() const  {return fDy;}
  float GetDZ() const  {return fDz;}
  
  void SetX(float x);
  void SetY(float y);
  void SetZ(float z);

  // move the detector in the direction perp to wires.
  void MovePerp(float offset);

  // Cancels last move
  void MoveBack();  

  float GetXc() const {return fXc;}
  float GetC2() const {return fC2;}
  float GetXc2() const {return fXc2;}
  float GetX2c2() const {return fX2c2;}
  float GetXs() const {return fXs;}
  float GetS2() const {return fS2;}
  float GetXs2() const {return fXs2;}
  float GetX2s2() const {return fX2s2;}
  float GetSc() const {return fSc;}
  float GetXsc() const {return fXsc;}
  float GetX2sc() const {return fX2sc;}
  float GetA() const {return fA;}
  float GetAc() const {return fAc;}
  float GetXac() const {return fXac;}
  float GetAs() const {return fAs;}
  float GetXas() const {return fXas;}

  float   GetIp() const {return fIp;}
  float   GetOp() const {return fOp;}
  float   GetPitch() const {return (fIp+fOp)/2.;}

  float GetTMin() const {return fCTMin;}
  float GetTMax() const {return fCTMax;}
  float GetTotMin() const {return fCTotMin;}

  float GetHalfSize() const {return fHalfsize;}

  float GetOffset() const {return fY*fCos + fZ*fSin;}

  // returns total number of tracks in active zone 
  int GetNTot() const {return fNtot;}

  // returns seen tracks  
  int GetNEff() const {return fNeff;}

  // returns background  
  int GetNBackg() const {return fNbackg;}

  // returns efficiency (na : -1)
  float GetEfficiency() const;

  // returns error on efficiency calculation
  float GetDEfficiency() const;

  // returns background probability
  float GetBGEfficiency() const; 

  // returns rotation matrix
  TRotMatrix* GetRotMatrix() {return fRotMatrix;}

  // returns all histograms
  std::vector<TH1*>& GetHistograms() {return fHists;}

  // sets BK radius
  void SetRMin(float rmin) {fR2Min = rmin*rmin;}

  // plots efficiency profile
  void PlotEffProfile();

  // calculates efficiency profile
  TH2F* EffProfile();

  // writes efficiency value on current pad
  void WriteEfficiency();

  // plots residuals 
  void PlotRes();

  // plots residuals cuts
  void DrawResCuts();

  // plots cluster time with applied cuts
  void PlotCTime();

  // plots cluster tot with applied cuts
  void PlotCTot();

  // plots cluster pos with applied cuts
  void PlotCPos();

  //float Wire2cm(int wire) { return (wire - fNwires/2. - 0.5)*fPitch;}
  float C2x(float c); // Channel number to Position X
 
  int GetNwires() const {return fNwires;}

  int GetNsigmat() const {return fNsigmat;}

  float GetCheckMin() const {return fCheckMin;}

  float GetCheckMax() const {return fCheckMax;}

  ClassDef(Tracker,0)
};

#endif








