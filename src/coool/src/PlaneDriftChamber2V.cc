#include "PlaneDriftChamber2V.h"
#include "PlaneDriftChamberPanel.h"
#include "ChipSinica.h"

ClassImp(PlaneDriftChamber2V);

const float PlaneDriftChamber2V::fSINICA_TICK = 1071.e-9; // in ms

//________________________________________________________________________________________________
PlaneDriftChamber2V::PlaneDriftChamber2V(const char *detname,int nchan, int center, int width)
  :Plane1V(detname,nchan,center,width)
{
  fVtemperature = AddVariable("_temperature",100,center-width,center+width, fNchan*fMAX_MULT);

  fVport = AddVariable("_port",100,center-width,center+width, fNchan*fMAX_MULT);
}

//________________________________________________________________________________________________
void PlaneDriftChamber2V ::NewRun(int runnumber) {
  std::vector<TH1*> histlist = GetHistoList();
  for( unsigned int i = 0; i < histlist.size(); i++ ) {
    char title[128];
    sprintf(title,"%s_Run%i",histlist[i]->GetName(),runnumber );
    histlist[i]->SetTitle( title );
  }
}

//________________________________________________________________________________________________
void PlaneDriftChamber2V::ControlPanel(const TGWindow* p, const TGWindow* main) 
{
  if (!fControlPanel) fControlPanel = new PlaneDriftChamberPanel(p, main, 100, 100, this);
}

//________________________________________________________________________________________________
void PlaneDriftChamber2V::Init(TTree* tree ) {

  Plane1V::Init(tree);
  
  std::string hitsname = fName + "_hits";

  //branch : temperature
  std::string temperaturename = fName + fVtemperature->GetName();
  fHtemperature=new TH2F(temperaturename.c_str(),temperaturename.c_str(),
			     fVport->GetNbins(),
			     fVport->GetMin(),
			     fVport->GetMax(),
			     fVtemperature->GetNbins(),
			     fVtemperature->GetMin(),
			     fVtemperature->GetMax());
  fHtemperature->SetXTitle("Port# (FEM#)");
  fHtemperature->SetYTitle("Temperature [C]");
  fHtemperature->SetOption("colz");
  AddHistogram(fHtemperature);
}



//________________________________________________________________________________________________
void PlaneDriftChamber2V::StoreDigit(CS::Chip::Digit* digit) {
  lDigits.push_back(digit);
  fNhits++;
  CS::ChipSinica::Digit *ds = dynamic_cast<CS::ChipSinica::Digit*>(digit);

  if (!ds) {
    std::cerr<<"Plane1V::StoreDigit ("<<GetName()<<"): a digit is not a Sinica one, strange...\n";
    return;
  }

  if (fFillfDigits) {
	  register int time,channel;
	   if (ds) {
		  time = ds->GetHitTime();
		  channel = ds->GetWire();
	  }
    if (fNhits==0) fDigits.clear();
    std::vector<double> data;
    data.push_back(time);
    fDigits.insert(make_pair(channel,CDigit1(channel,data)));
  }
}

//________________________________________________________________________________________________
void PlaneDriftChamber2V::EndEvent(const CS::DaqEvent &event) {

  int reftemp = -1;
  
  if (thr_flag) TThread::Lock();

  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {

	  CS::ChipSinica::Digit *ds = dynamic_cast<CS::ChipSinica::Digit*>(*ii);

	  if (!ds) {
	    std::cerr<<"Plane1V::StoreDigit ("<<GetName()<<"): a digit is not a Sinica one, strange...\n";
	    continue;
	  }

	  register double time,channel,port,temperature;
    if (ds){
      CS::ChipSinica::DataID id = ds->GetDataID();
      port = id.u.s.port;

      /* 
       * when dividing by time unit "time" is in digits.
       * leave as
       * time = (double) (ds->GetTimeDecoded()
       * to have time in nanoseconds
      */
      time = (double) (ds->GetTimeDecoded() / ds->GetTimeUnit()); 
      channel = ds->GetWire();
      temperature = ds->GetTemperature();
    }

    if( fVch->Test(channel) ) {

      if( fVt->Test(time) ) {
	fVch->Store(channel);
	fVt->Store(time);
	fHch->Fill(channel);
	fHt->Fill(time);
	fHtvsch->Fill(channel,time);
	fNhitsKept++;
      }

      if( fOnTrigTVarID->Test(time) )  {
        fOnTrigTVarID->Store(time);
        fHonTT->Fill(time);        // on trigger time distribution
        fHonTP->Fill(channel);     // on trigger profile
        fNhitsKeptOnTrig++;
      }

      if( fOffTrigTVarID->Test(time) )  {
        fOffTrigTVarID->Store(time);
        fHoffTT->Fill(time);        // off trigger time distribution
        fHoffTP->Fill(channel);     // off trigger profile
        fNhitsKeptOffTrig++;
      }
    }

    if((fRateCounter%Plane1V::fRATE_UPDATE)==0) {  // to limit update rate
      // Since every port/FEM has same temperature for all digits,
      // make sure we dont keep repeating entries. 
      // only adds entries when temp changes (will never happen on same FEM)
      if(reftemp != temperature){     
	fVtemperature->Store(temperature);
	fHtemperature->Fill(port, temperature);
	reftemp = temperature;
      }
    }
  }
  fHhits->Fill(fNhitsKept);

  if((fRateCounter%Plane1V::fRATE_UPDATE)==0) {

    // corrected on-trigger profile
    if( fOffTrigTVarID->GetMax() != fOffTrigTVarID->GetMin() ) {
      float scale = (float) ( fOnTrigTVarID->GetMax() - fOnTrigTVarID->GetMin() )/( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() );
      for(register int i=1; i<=fHonTP->GetNbinsX(); i++)
	fHoncorrTP->SetBinContent(i, (float) fHonTP->GetBinContent(i) - scale*fHoffTP->GetBinContent(i) );
    }

    // rates histogram
    float timewin = ( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() )*fSINICA_TICK*fRateCounter;
    if( timewin > 0 )
      for(register int i=1; i<=fHonTP->GetNbinsX(); i++) {
	fHrates->SetBinContent(i,(float) fHoffTP->GetBinContent(i)/timewin);
	fHrates->SetBinError(i,(float) fHoffTP->GetBinError(i)/timewin);
      }
  }
  if (thr_flag) TThread::UnLock();
}


