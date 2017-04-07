#ifndef __PlaneMuonWallA__
#define __PlaneMuonWallA__
#include <stdio.h>
#include "Plane.h"
#include "PlanePanel.h"
#include "TThread.h"

#include "TTree.h"

class PlaneMuonWallAPanel;

// Muon Wall A, S.Trusov

class PlaneMuonWallA : public  Plane {
  static TFile* refFile;
 protected:
#ifndef __CINT__
  std::map<int,CDigit1> fDigits;
  /// class needed for calibration purposes, here for any chance TSV
     class ChannelCalib {
     public:
       int ch;
       float t0;
       ChannelCalib() : ch(0),t0(0) {}
       ChannelCalib(const char *s) {
	 if(2 != sscanf(s,"%d%f",&ch,&t0)) {
	   throw CS::Exception("PlaneMuonWallA::ChannelCalib : bad line \"%s\"",s);
	   std::cerr<<"bad line, exception not caught !"<<std::endl;
	 }
       }
     };
     std::vector<ChannelCalib> calib_data;  // calibration for every channel
  
     friend istream& operator>>(istream& in, PlaneMuonWallA::ChannelCalib &c) {
       in>>c.ch;
       in>>c.t0;
       return in;
     }
     /// end calibrations
#endif  // __CINT__
  static const int   fMAX_MULT;     // maximum multiplicity per tdc channel
  static const float fF1_TICK;      // F1 TDC tick in ms
  static const int   fRATE_UPDATE;  // Period in event numbers of rate histo update
  static const float fHIT_TMIN, fHIT_TMAX; // Hit time window
  static const float fMT_T0;               // Coinc time
  static const float fCL_TMIN, fCL_TMAX;   // Cluster time window

  int          fNchan;        // number of channels
//   int       *fChannel;        // channel for each digit
//   int          *fTime;        // data for each digit
//   int          *fKey;         // signature 0/1/-1 for each digit
//   bool fIsHR;  
  int fCurChan;          // current channel looked at


  std::set<int> fNoisyCh;          // Noisy channels
  int fNCl;                   // number of clusters passing the cuts

  TH1F *fHhits;              // Multiplicity histogram
  //  TH1F *fHch;                 // Hit profile
  TH1F *fHt;                  // Time distribution

  TH1F *fHch_up;
  TH1F *fHch_down;

  TH1F *fHdet_up;
  TH1F *fHdet_down;
  
  Variable *fVch, *fVt;       
  Variable *fVch_up, *fVch_down;
  Variable *fVdet_up, *fVdet_down;
/*    Variable *fVcch, *fVct;   */ // clusters

  std::set<int> fMissingChans;          // missing channels
  std::set<int> fNoisyChans;            // noisy channels
 public:
  PlaneMuonWallA(const char *detname,int nchan, int center, int width);
  virtual ~PlaneMuonWallA();

  void StoreDigit(int channel, int data);  /// Passes a digit to the plane
  void StoreDigit(int channel, const std::vector<double>& data);
#ifndef __CINT__
  void StoreDigit(CS::Chip::Digit* digit);

  /*! Applies the cuts, fills histograms and tree. 
   * This function must be called at the end of the event
   */  
  void EndEvent(const CS::DaqEvent &event);  
#endif

  void print_hits() {
    for(std::map<int,CDigit1>::const_iterator i= fDigits.begin();
	i!=fDigits.end(); i++) {
      if(i->second.dt[0]) {
  	std::cout<<"TSV: "<<GetName()
  	    <<" "<<i->second.ch
  	    <<" "<<i->second.dt[0]
  	    <<" "<<i->second.dt[1]
  	    <<'\n';
      }
    }
  }

  static bool OpenRefFile(const char* fileName);
  static void CloseRefFile();

  virtual void Clusterize();                    // do the clustering 
  virtual void Init(TTree* tree =0);            // book histograms and branch the tree
  int GetNchannels() {return fNchan;}  // return number of channels
  void ResetHistograms();

  virtual void ControlPanel(const TGWindow *p, const TGWindow *main);
  void TextOutput(ostream& out); // ??

  float* GetTimeValues() {return fVt->GetValues();}
  float* GetChanValues() {return fVch->GetValues();}
  void SetChannel(int channel) {fCurChan=channel;}
//   void SetHRTime(bool ishr) {fIsHR=ishr;}

  ClassDef(PlaneMuonWallA, 1)
};

#endif




