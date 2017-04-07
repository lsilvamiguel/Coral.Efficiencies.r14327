#include "PlaneDriftChamber.h"
#include "PlaneDriftChamberPanel.h"
#include "ChipF1.h"

ClassImp(PlaneDriftChamber);

// implement your member functions here
//const float PlaneDriftChamber::fF1_TICK = 128.e-9; // in ms
//const int PlaneDriftChamber::fRATE_UPDATE = 100;

//________________________________________________________________________________________________
PlaneDriftChamber::PlaneDriftChamber(const char *detname,int nchan, int center, int width)
  :Plane1V(detname,nchan,center,width)
{
  
  // new variable for on and off trigger time
  // fOnTrigTVarID = AddVariable("_t_on_trigger",100,center-width,center+width, fNchan*fMAX_MULT);
  // fOffTrigTVarID = AddVariable("_t_off_trigger",100,center-width,center+width, fNchan*fMAX_MULT);
}

//________________________________________________________________________________________________
void PlaneDriftChamber::NewRun(int runnumber) {
  std::vector<TH1*> histlist = GetHistoList();
  for( unsigned int i = 0; i < histlist.size(); i++ ) {
    char title[128];
    sprintf(title,"%s_Run%i",histlist[i]->GetName(),runnumber );
    histlist[i]->SetTitle( title );
  }
}

//________________________________________________________________________________________________
void PlaneDriftChamber::ControlPanel(const TGWindow* p, const TGWindow* main) 
{
  if (!fControlPanel) fControlPanel = new PlaneDriftChamberPanel(p, main, 100, 100, this);
}

//________________________________________________________________________________________________
//  void PlaneDriftChamber::Init(TTree* tree ) {

//  Plane1V::Init(tree);
//   fRateCounter = 0;
  
//    // On-Trigger Time Distribution
//    string otname = fName + "_t_on_trig";
//    fHonTT=new TH1F(otname.c_str(),otname.c_str(),   //onTT is "off Trigger Time"
//  		  fOnTrigTVarID->GetNbins(),
//  		  fOnTrigTVarID->GetMin(),
//  		  fOnTrigTVarID->GetMax());
//    AddHistogram(fHonTT);

//    // Off-Trigger Time Distribution
//    string tname = fName + "_t_off_trig";
//    fHoffTT=new TH1F(tname.c_str(),tname.c_str(),   //offTP is "off Trigger Time"
//  		   fOffTrigTVarID->GetNbins(),
//  		   fOffTrigTVarID->GetMin(),
//  		   fOffTrigTVarID->GetMax());
//    AddHistogram(fHoffTT);

  
//    // On-Trigger Profile
//    string cname = fName + "_ch_on_trig";
//    fHonTP=new TH1F(cname.c_str(),cname.c_str(),   //onTP is "on Trigger Profile"
//  		  fVch->GetNbins(),
//  		  fVch->GetMin(),
//  		  fVch->GetMax());
//    AddHistogram(fHonTP);
  
//    // Off-Trigger Profile
//    string pname = fName + "_ch_off_trig";
//    fHoffTP=new TH1F(pname.c_str(),pname.c_str(),   //offTP is "off Trigger Profile"
//  		   fVch->GetNbins(),
//  		   fVch->GetMin(),
//  		   fVch->GetMax());
//    AddHistogram(fHoffTP);
  

//    // Corrected On-Trigger Profile
//    string corrname = fName + "_ch_on_corr_trig";
//    fHoncorrTP=new TH1F(corrname.c_str(),corrname.c_str(),   //oncorrTP is "on Trigger Time Profile (corrected)"
//  		      fVch->GetNbins(),
//  		      fVch->GetMin(),
//  		      fVch->GetMax());
//    AddHistogram(fHoncorrTP);

//    // Rates HISTOGRAMS
//    string rname = fName + "_rates";
//    fHrates=new TH1F(rname.c_str(),rname.c_str(),
//  		   fVch->GetNbins(),
//  		   fVch->GetMin(),
//  		   fVch->GetMax());
//    fHrates->SetYTitle("rates per channel (kHz)");
//    fHrates->SetTitleOffset(1.2,"Y");
//    AddHistogram(fHrates);
  
//  }

//  //________________________________________________________________________________________________
//  void PlaneDriftChamber::Reset( void )
//  {
//    Plane::Reset();
//    fNhitsKeptOnTrig = 0;
//    fNhitsKeptOffTrig = 0;
//    return;
//  }
  

//  //________________________________________________________________________________________________
//  void PlaneDriftChamber::ResetHistograms() {
//    Plane::ResetHistograms();
//    fIsFirstEvent = true;
//  }
  
//  #ifndef __CINT__
//  //________________________________________________________________________________________________
//  void PlaneDriftChamber::EndEvent( int reftime, const CS::DaqEvent &event) {

//    // get Run number at first event
//    if( fIsFirstEvent ) {
//      fIsFirstEvent = false;
    
//      // change all histogram titles
//      vector<TH1*> histlist = GetHistoList();
//      for( unsigned int i = 0; i < histlist.size(); i++ ) {
//        char title[128];
//        sprintf(title,"%s_Run%i",histlist[i]->GetName(),event.GetRunNumber() );
//        histlist[i]->SetTitle( title );
//      }
//    }

//    if (thr_flag) TThread::Lock();

//    for (int i=0; i<fNhits; i++) {

//      int time = CS::ChipF1::TimeDifference(fTime[i],reftime);
//      int channel=fChannel[i];


//      if(fVch->Test(channel) ) {
//        if( fVt->Test(time) ) {
//          fVt->Store(time);
//          fHt->Fill(time);          // time distrib
//          fHtvsch->Fill(channel,time);  // channel vs time
//          fVch->Store(channel);
//          fHch->Fill(channel);       // profile
//          fNhitsKept++;
//        }
      
//        if( fOnTrigTVarID->Test(time) )  {
//          fOnTrigTVarID->Store(time);
//  	fHonTT->Fill(time);        // on trigger time distribution
//          fHonTP->Fill(channel);  // on trigger profile
//          fNhitsKeptOnTrig++;
//        }
      
//        if( fOffTrigTVarID->Test(time) )  {
//          fOffTrigTVarID->Store(time);
//          fHoffTT->Fill(time);        // off trigger time distribution
//          fHoffTP->Fill(channel);  // off trigger profile
//          fNhitsKeptOffTrig++;
//        }
//      }
//    }
//    fHhits->Fill(fNhitsKept);         // multiplicities

//    if((fRateCounter%PlaneDriftChamber::fRATE_UPDATE)==0) {

    
//      // corrected on-trigger profile
//      if( fOffTrigTVarID->GetMax() != fOffTrigTVarID->GetMin() ) {
//        float scale = (float) ( fOnTrigTVarID->GetMax() - fOnTrigTVarID->GetMin() )/( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() );
//        for(register int i=1; i<=fHonTP->GetNbinsX(); i++)
//        fHoncorrTP->SetBinContent(i, (float) fHonTP->GetBinContent(i) - scale*fHoffTP->GetBinContent(i) );
//      }  
    
//      // rates histogram
//      float timewin = ( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() )*fF1_TICK*fRateCounter;
//      if( timewin > 0 )
//      for(register int i=1; i<=fHonTP->GetNbinsX(); i++) {
//  	    fHrates->SetBinContent(i,(float) fHoffTP->GetBinContent(i)/timewin);
//  	    fHrates->SetBinError(i,(float) fHoffTP->GetBinError(i)/timewin);
//      }
    
//    }
    
//    if (thr_flag) TThread::UnLock();
//  }
//  #endif
