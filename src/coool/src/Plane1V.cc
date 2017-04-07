#include "Plane1V.h"
#include "Plane1VPanel.h"
#include "ChipF1.h"
#include "ChipGandalf.h"
#include "TriggerTime.h"
#include "TFile.h"
#include "TWebFile.h"
#include "Reference.h"

#include "TF1.h"

ClassImp(Plane1V);

const int Plane1V::fMAX_MULT = 8;
const float Plane1V::fF1_TICK = 128.e-9; // in ms
const int Plane1V::fRATE_UPDATE = 100;


Plane1V::Plane1V(const char *detname,int nchan, int center, int width)
  :
  Plane(detname), fFillfDigits(true), fNchan(nchan),
  fCenter(center), fWidth(width) {


  fVch = AddVariable("_ch",fNchan,0,fNchan,fNchan*fMAX_MULT);
  fVt = AddVariable("_t",100,center-width,center+width,fNchan*fMAX_MULT);
  fOnTrigTVarID = AddVariable("_t_on_trigger",100,center-width,center+width, fNchan*fMAX_MULT);
  fOffTrigTVarID = AddVariable("_t_off_trigger",100,center-width,center+width, fNchan*fMAX_MULT);

//   fChannel=new int[fNchan*fMAX_MULT];
//   fTime=new int[fNchan*fMAX_MULT];
}

Plane1V::~Plane1V() {
//   delete fTime; delete fChannel;
}

void Plane1V::Init(TTree *tree) {

  fRateCounter = 0;

  // hits --------------------------------------------------------

  //multiplicity within the time cut on the whole detector
  std::string hitsname = fName + "_hits";
  fHhits=new TH1F_Ref(hitsname.c_str(),hitsname.c_str(),10,0,10,fRateCounter);
  AddHistogram(fHhits);
  std::string hitsleavlist = hitsname + "/I";

  //branch : hit map
  std::string chname = fName + fVch->GetName();
  fHch=new TH1F_Ref(chname.c_str(),chname.c_str(),
		    fVch->GetNbins(),
		    fVch->GetMin(),
		    fVch->GetMax(),
		    fRateCounter);
  AddHistogram(fHch);
  std::string chleavlist = chname + "[" + hitsname +"]/F";

  //branch : time
  std::string tname = fName + fVt->GetName();
  fHt=new TH1F_Ref(tname.c_str(),tname.c_str(),
		   fVt->GetNbins(),
		   fVt->GetMin(),
		   fVt->GetMax(),
		   fRateCounter);
  AddHistogram(fHt);
  std::string tleavlist = tname + "[" + hitsname +"]/F";

  // time vs channel histo
  std::string tcname = fName + "_tVSch";
  fHtvsch=new TH2F(tcname.c_str(),tcname.c_str(),
		   fVch->GetNbins(),
		   fVch->GetMin(),
		   fVch->GetMax(),
		   fVt->GetNbins(),
		   fVt->GetMin(),
		   fVt->GetMax());
  fHtvsch->SetOption("col");
  AddHistogram(fHtvsch);

  // On-Trigger Time Distribution
  std::string onttname = fName + "_t_on_trig";
  fHonTT=new TH1F_Ref(onttname.c_str(),onttname.c_str(),   //onTT is "on Trigger Time"
		      fOnTrigTVarID->GetNbins(),
		      fOnTrigTVarID->GetMin(),
		      fOnTrigTVarID->GetMax(),
		      fRateCounter);
  AddHistogram(fHonTT);

  // Off-Trigger Time Distribution
  std::string offttname = fName + "_t_off_trig";
  fHoffTT=new TH1F(offttname.c_str(),offttname.c_str(),   //offTP is "off Trigger Time"
		   fOffTrigTVarID->GetNbins(),
		   fOffTrigTVarID->GetMin(),
		   fOffTrigTVarID->GetMax());
  AddHistogram(fHoffTT);


  // On-Trigger Profile
  std::string offtcname = fName + "_ch_on_trig";
  fHonTP=new TH1F_Ref(offtcname.c_str(),offtcname.c_str(),   //onTP is "on Trigger Profile"
		  fVch->GetNbins(),
		  fVch->GetMin(),
		  fVch->GetMax(),fRateCounter);
  AddHistogram(fHonTP);

  // Off-Trigger Profile
  std::string offtpname = fName + "_ch_off_trig";
  fHoffTP=new TH1F(offtpname.c_str(),offtpname.c_str(),   //offTP is "off Trigger Profile"
		   fVch->GetNbins(),
		   fVch->GetMin(),
		   fVch->GetMax());
  AddHistogram(fHoffTP);

  // Corrected On-Trigger Profile
  std::string ontcorrname = fName + "_ch_on_corr_trig";
  fHoncorrTP=new TH1F(ontcorrname.c_str(),ontcorrname.c_str(),
		      fVch->GetNbins(),
		      fVch->GetMin(),
		      fVch->GetMax());
  AddHistogram(fHoncorrTP);

  // Rates HISTOGRAMS
  std::string rname = fName + "_rates";
  fHrates=new TH1F_Ref(rname.c_str(),rname.c_str(),
		   fVch->GetNbins(),
		   fVch->GetMin(),
		   fVch->GetMax(),
                                   fRateCounter,true);
  fHrates->SetYTitle("rates per channel (kHz)");
  fHrates->SetTitleOffset(1.2,"Y");
  fHrates->SetOption("hist");
  AddHistogram(fHrates);

  if (tree) {
    fIsInTree = true;
    tree->Branch(hitsname.c_str(),&fNhitsKept,
		 hitsleavlist.c_str(),32000);
    tree->Branch(chname.c_str(),fVch->GetValues(),
		 chleavlist.c_str(),32000);
    tree->Branch(tname.c_str(),fVt->GetValues(),
		 tleavlist.c_str(),32000);
  }
  if (fReferenceDirectory) {
    ((TH1F_Ref*)fHhits)->SetReference(fReferenceDirectory);
    ((TH1F_Ref*)fHch)->SetReference(fReferenceDirectory);
    ((TH1F_Ref*)fHt)->SetReference(fReferenceDirectory);
    ((TH1F_Ref*)fHonTT)->SetReference(fReferenceDirectory);
    ((TH1F_Ref*)fHonTP)->SetReference(fReferenceDirectory);
    ((TH1F_Ref*)fHrates)->SetReference(fReferenceDirectory);
  }
}

void Plane1V::Reset( void ) {
  Plane::Reset();
  fNhitsKeptOnTrig = 0;
  fNhitsKeptOffTrig = 0;
  return;
}

void Plane1V::ResetHistograms() {
  Plane::ResetHistograms();
}


void Plane1V::StoreDigit(int channel, int time) {
  std::cerr<<"Plane1V::StoreDigit(int, int) ("<<GetName()<<"): Error ! this method is obsolete !!!\n";
//   if (fNhits==0) fDigits.clear();
//   std::vector<double> data;
//   data.push_back(time);
//   fDigits.insert(make_pair(channel,CDigit1(channel,data)));
// 
//   if (fNhits < fNchan*fMAX_MULT) {
//     fTime[fNhits]=time;
//     fChannel[fNhits]=channel;
//     fNhits++;
//   }
}


void Plane1V::StoreDigit(CS::Chip::Digit* digit) {
//   std::vector<float> data=digit->GetNtupleData();
//   this->StoreDigit((int)data[0],(int)data[1]);
  lDigits.push_back(digit);
  fNhits++;
  CS::ChipF1::Digit *df1 = dynamic_cast<CS::ChipF1::Digit*>(digit);
  CS::ChipGandalf::DigitGADC *dg = dynamic_cast<CS::ChipGandalf::DigitGADC*>(digit);
  if (!df1 && !dg) {
    std::cerr<<"Plane1V::StoreDigit ("<<GetName()<<"): a digit is neither a F1 nor a Gandalf one, strange...\n";
    return;
  }

  if (fFillfDigits) {
	  register int time=-1, channel=-1;
	  if (df1) {
		  time = df1->GetTime();
		  channel = df1->GetChannel();
	  }
	  if (dg) {
		  time = dg->getTime();
		  channel = dg->getChannel();
	  }
    if (fNhits==0) fDigits.clear();
    std::vector<double> data;
    data.push_back(time);
    fDigits.insert(make_pair(channel,CDigit1(channel,data)));
  }
}


void Plane1V::EndEvent(const CS::DaqEvent &event) {

//   const float reftime = fIsHR ? event.GetTT().GetTimeHigh() :
//                                 event.GetTT().GetTimeNorm();

  if (thr_flag) TThread::Lock();

//   for (int i=0; i<fNhits; i++) {
//     double time = CS::ChipF1::TimeDifference(fTime[i],reftime);
//     int channel=fChannel[i];
  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {

	  CS::ChipF1::Digit *df1 = dynamic_cast<CS::ChipF1::Digit*>(*ii);
	  CS::ChipGandalf::DigitGADC *dg = dynamic_cast<CS::ChipGandalf::DigitGADC*>(*ii);
	  if (!df1 && !dg) {
	    std::cerr<<"Plane1V::StoreDigit ("<<GetName()<<"): a digit is neither a F1 nor a Gandalf one, strange...\n";
	    continue;
	  }

    register double time=-1, channel=-1;
    if (df1){
    	time = (double) (df1->GetTimeDecoded() / df1->GetTimeUnit());
    	channel = df1->GetChannel();
    }
    if (dg) {
    	time = (double) (dg->GetTimeDecoded() / dg->GetTimeUnit());
    	channel = dg->getChannel();
    }

//     register int timeraw = iii->GetTime();
//     register double timeunit = iii->GetTimeUnit();

// std::cerr<<"Plane1V::StoreDigit: GetName "<<GetName()<<" time "<<time<<" channel "<<channel<<" timeraw "<<timeraw<<" timeunit "<<timeunit<<std::endl;

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
        fHonTP->Fill(channel);  // on trigger profile
        fNhitsKeptOnTrig++;
      }

      if( fOffTrigTVarID->Test(time) )  {
        fOffTrigTVarID->Store(time);
        fHoffTT->Fill(time);        // off trigger time distribution
        fHoffTP->Fill(channel);  // off trigger profile
        fNhitsKeptOffTrig++;
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
    float timewin = ( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() )*fF1_TICK*fRateCounter;
    if( timewin > 0 )
      for(register int i=1; i<=fHonTP->GetNbinsX(); i++) {
	fHrates->SetBinContent(i,(float) fHoffTP->GetBinContent(i)/timewin);
	fHrates->SetBinError(i,(float) fHoffTP->GetBinError(i)/timewin);
      }
  }
  if (thr_flag) TThread::UnLock();
}

void Plane1V::ControlPanel(const TGWindow* p, const TGWindow* main) {

  if (!fControlPanel) fControlPanel = new Plane1VPanel(p, main, 100, 100, this);
}

void Plane1V::TimeSpectrum() {

  std::string hdname = fName + "_1ch";
  if(fCurChan<fNchan) {
    fHtvsch->ProjectionY(hdname.c_str(),fCurChan+1,fCurChan+1,"")->Draw();
  }
}

void Plane1V::TextOutput(ostream& out) {

  if(ChannelProblems()) {
    out<<"\tmissing";
//     out<<fMissingChans<<std::endl;
    for (std::set<int>::iterator it=fMissingChans.begin(); it!=fMissingChans.end(); it++)
      out<<" "<<(*it);
    out<<std::endl;
    out<<"\tnoisy";
//     out<<fNoisyChans<<std::endl;
    for (std::set<int>::iterator it=fNoisyChans.begin(); it!=fNoisyChans.end(); it++)
      out<<" "<<(*it);
    out<<std::endl;
  } else {
    out<<"\tmissing ***"<<std::endl;
    out<<"\tnoisy ***"<<std::endl;
  }
}


bool Plane1V::ChannelProblems() {
  // fHch must have at least one entry, otherwise fit will fail
  if (fHch->GetEntries()==0)
    return false;

  const int minaverage = 30;

  // missing channels are found in hit profile histogram, if we have enough stat
  TF1 pol0("pol0","pol0");
  fHch->Fit(&pol0,"WNQ");
  float average = pol0.GetParameter(0);
  if(average<minaverage) return false;

  fMissingChans.clear();
  fNoisyChans.clear();

  float min = average/20;
  float max = average*5;

  for(int i=1; i<=fHch->GetNbinsX(); i++) {
    if(fHch->GetBinContent(i)<min) fMissingChans.insert(i-1);
    if(fHch->GetBinContent(i)>max) fNoisyChans.insert(i-1);
  }

  return true;
}






