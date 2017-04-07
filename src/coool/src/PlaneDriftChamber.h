#ifndef __PlaneDriftChamber__
#define __PlaneDriftChamber__

#include "Plane1V.h"

#include "TTree.h"
#include "TThread.h"

/// Plane for drift chambers
class PlaneDriftChamber : public Plane1V {
 
 private:
  
  //  bool fIsFirstEvent; 
   
/*    // F1 TDC tick in ms */
/*    static const float fF1_TICK; */

/*    // indexes in histo and variable vectors  */
/*    Variable *fOnTrigTVarID; */
/*    Variable *fOffTrigTVarID; */

/*    TH1F *fHonTT; //on-Trigger time distribution */
/*    TH1F *fHonTP; //on-Trigger profile */
/*    TH1F *fHoffTT; //off-Trigger time distribution */
/*    TH1F *fHoffTP; //off-Trigger profile */
/*    TH1F *fHoncorrTP; //corrected on-trigger profile */
/*    TH1F *fHrates; //rates */
  
/*    int fNhitsKeptOnTrig;   // number of hits kept in the on-trigger time window */
/*    int fNhitsKeptOffTrig;   // number of hits kept in the off-trigger time window */
 

/*    // Period in event numbers of rate histo update */
/*    static const int fRATE_UPDATE; */

public:
  
  PlaneDriftChamber( const char *detname,int nchan, int center, int width );

  void NewRun(int runnumber);
  virtual void ControlPanel(const TGWindow* p, const TGWindow* main); 

  //  void Init( TTree* tree =0 );
  //void Reset( void );
  //void ResetHistograms( void );
  
  //  #ifndef __CINT__
  //void EndEvent( int reftime, const CS::DaqEvent &event );
  //#endif
  
  //~PlaneDriftChamber() {}

  ClassDef(PlaneDriftChamber,1)
};

#endif




