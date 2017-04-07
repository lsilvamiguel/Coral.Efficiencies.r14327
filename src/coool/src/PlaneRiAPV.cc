#include <cmath>

#include "PlaneRiAPV.h"
#include "PlaneRiAPVPanel.h"
#include "ChipAPVRICH.h"
#include <TF1.h>
#include <algorithm>
#include <stdexcept>
#include <utility>

#define RiAPV_A2_CUT 12
#define RiAPV_A2_FIT_MIN 8.
#define RiAPV_A2_FIT_MAX 200.

ClassImp(PlaneRiAPV);


//-----------------------------------------------------------------------------


PlaneRiAPV::PlaneRiAPV(const char *detname,int nchan, int center, int width)
    : PlaneAPV(detname,nchan,center,width), fMinChan(0),
      fHampRatio(0), fHampSubRatio(0), fHtimeRatio(0),
      fHa1a2Ratio(0),
      fHchvsa2(0), fHchvsa2ped(0), fHchvsa2CM(0),
      fHoccup(0), fHampSubRap(0),
      fHa0CM(0), fHa1CM(0), fHa2CM(0),
      fHa0s(0), fHa1s(0), fHa2s(0),
      fHa0CMs(0), fHa1CMs(0), fHa2CMs(0),
      fHch2D(0), fHchamp2D(0),
      fHavgamp(0), fHsigamp(0), fHavgampCM(0), fHsigampCM(0),
      fHavgamp2D(0), fHsigamp2D(0), fHavgampCM2D(0), fHsigampCM2D(0),
      fHCMa2(0),
      fHCMa2apv(0),
      fHevtamp2D(0), fHevtampcut2D(0),
      fHtvsTCSph(0),
      fHapvvsa2mean(0),
      fHa2tmp(0),
      fCurChan(0),
      fOldAPVBoard(false), fCardRotated(false)
{
  fVxpx = AddVariable("_xpx",72,0,72,fNchan*fMAX_MULT);
  fVypx = AddVariable("_ypx",72,0,72,fNchan*fMAX_MULT);
  fVa0->SetRange(-100, 900);
  fVa1->SetRange(-100, 900);
  fVa2->SetRange(-100, 900);
//  fVt->SetRange(-500, 500);
  fVt->SetRange(0, 4);
  fVt->SetNbins(200);
//  fVxpx->SetNbins(100);
//  fVxpx->SetRange(-500, 500); 
#if USE_DATABASE==1
  calib_arr = 0;
#endif
}

//-----------------------------------------------------------------------------


void PlaneRiAPV::Init(TTree* tree) {

//  fNchan *= 128; fNchan /= 108;     // to take into account the 20 channels/APV which are not connected to the RICH
  if (fNchan > MAX_STAT_RiAPV) {
    std::cerr << "\n ERROR in PlaneRiAPV::Init: too many channels !!!!\n";
    std::cerr << "  nb channels = "<<fNchan<<" maximum is "<<MAX_STAT_RiAPV<<std::endl<<std::endl;
    throw "ERROR in PlaneRiAPV::Init: too many channels !!!!";
  }
  if ((fNchan%MAX_APV_NBCH) != 0) {
    std::cerr << "\n ERROR in PlaneRiAPV::Init: nb of channels not a multiple of 128 !!!!\n";
    std::cerr << "  nb channels = "<<fNchan<<std::endl<<std::endl;
    throw "ERROR in PlaneRiAPV::Init: nb of channels not a multiple of 128 !!!!";
  }

  PlaneAPV::Init(tree);
//  PlaneRiAPV::initconndata();

  for (register int iapv = 0; iapv < MAX_NBAPV_RiAPV; iapv++) {
    fCMapv[iapv][0] = 0.;
    fCMapv[iapv][1] = 0.;
    fCMapv[iapv][2] = 0.;
  }

  for (register int ii =0; ii < MAX_STAT_RiAPV; ii++) {
    fStat[ii].sum = 0.; fStat[ii].sum2 = 0.; fStat[ii].nb = 0;
    fStatCM[ii].sum = 0.; fStatCM[ii].sum2 = 0.; fStatCM[ii].nb = 0;
  }
//   fVa2CM = AddVariable("_a2CM",600,-100,500,fNchan*fMAX_MULT);

  // Time from a1/a2 vs tcsphase
  std::string name = fName + "_a1overa2_vs_TCSphase";
  fHtvsTCSph=new TH2F(name.c_str(),name.c_str(), 40, -10, 30,
		      100, 0. ,2.); 
  fHtvsTCSph->SetOption("colz");
  AddHistogram(fHtvsTCSph);

  // a1/a2 amplitude ratio
  name = fName + "_a1overa2";
  fHa1a2Ratio = new TH1F_Ref(name.c_str(),name.c_str(),
			100, 0., 2.,fRateCounter);
  ((TH1F_Ref*)fHa1a2Ratio)->SetReference(fReferenceDirectory);
  AddHistogram(fHa1a2Ratio);
  fHa1a2Ratio->GetXaxis()->SetTitle("a1/a2");

  if (fExpertHistos) {

    // Amplitude 1st sample
    std::string a0CMname = fName + fVa0->GetName() + "_CMcorr";
    fHa0CM=new TH1F(a0CMname.c_str(),a0CMname.c_str(),
		     fVa0->GetNbins(),
		     fVa0->GetMin(),
		     fVa0->GetMax()); 
    AddHistogram(fHa0CM);

    // Amplitude 2nd sample
    std::string a1CMname = fName + fVa1->GetName() + "_CMcorr";
    fHa1CM=new TH1F(a1CMname.c_str(),a1CMname.c_str(),
		    fVa1->GetNbins(),
		    fVa1->GetMin(),
		    fVa1->GetMax()); 
    AddHistogram(fHa1CM);

    // Amplitude 3rd sample
    std::string a2CMname = fName + fVa2->GetName() + "_CMcorr";
    fHa2CM=new TH1F(a2CMname.c_str(),a2CMname.c_str(),
		     fVa2->GetNbins(),
		     fVa2->GetMin(),
		     fVa2->GetMax()); 
    AddHistogram(fHa2CM);

    // Amplitude 1st sample
    std::string a0name = fName + fVa0->GetName() + "_vs_spill";
    fHa0s=new TH2F(a0name.c_str(),a0name.c_str(), 200, 0., 200.,
		     fVa0->GetNbins(),
		     fVa0->GetMin(),
		     fVa0->GetMax());
    AddHistogram(fHa0s);

    // Amplitude 2nd sample
    std::string a1name = fName + fVa1->GetName() + "_vs_spill";
    fHa1s=new TH2F(a1name.c_str(),a1name.c_str(), 200, 0., 200.,
		     fVa1->GetNbins(),
		     fVa1->GetMin(),
		     fVa1->GetMax());
    AddHistogram(fHa1s);

    // Amplitude 3rd sample
    std::string a2name = fName + fVa2->GetName() + "_vs_spill";
    fHa2s=new TH2F(a2name.c_str(),a2name.c_str(), 200, 0., 200.,
		     fVa2->GetNbins(),
		     fVa2->GetMin(),
		     fVa2->GetMax());
    AddHistogram(fHa2s);

    // Amplitude 1st sample CM
    name = fName + fVa0->GetName() + "CM_vs_spill";
    fHa0CMs=new TH2F(name.c_str(),name.c_str(), 200, 0., 200.,
		     fVa0->GetNbins(),
		     fVa0->GetMin(),
		     fVa0->GetMax());
    AddHistogram(fHa0CMs);

    // Amplitude 2nd sample CM
    name = fName + fVa1->GetName() + "CM_vs_spill";
    fHa1CMs=new TH2F(name.c_str(),name.c_str(), 200, 0., 200.,
		     fVa1->GetNbins(),
		     fVa1->GetMin(),
		     fVa1->GetMax());
    AddHistogram(fHa1CMs);

    // Amplitude 3rd sample CM
    name = fName + fVa2->GetName() + "CM_vs_spill";
    fHa2CMs=new TH2F(name.c_str(),name.c_str(), 200, 0., 200.,
		     fVa2->GetNbins(),
		     fVa2->GetMin(),
		     fVa2->GetMax());
    AddHistogram(fHa2CMs);

    // Channel amplitude ratio ("banana")
    name = fName + "_a0/a2_vs_a1/a2";
    fHampRatio = new TH2F(name.c_str(),name.c_str(),
			  100, -0.01, 1.99, 100, -0.01, 1.99);
    AddHistogram(fHampRatio);
    fHampRatio->GetXaxis()->SetTitle("a1/a2");
    fHampRatio->GetYaxis()->SetTitle("a0/a2");
    fHampRatio->SetOption("colz");

    // Channel amplitudes substraction correlation
    name = fName + "_a2-a1_vs_a1-a0";
    fHampSubRatio = new TH2F(name.c_str(),name.c_str(),
			  200, -100., 100., 200, -100., 100.);
    AddHistogram(fHampSubRatio);
    fHampSubRatio->GetXaxis()->SetTitle("a1-a0");
    fHampSubRatio->GetYaxis()->SetTitle("a2-a1");
    fHampSubRatio->SetOption("colz");

    // Channel amplitudes substraction ratio
    name = fName + "_a2-a1_dividedby_a1-a0";
    fHampSubRap = new TH1F(name.c_str(),name.c_str(), 400, -3., 5.);
    AddHistogram(fHampSubRap);
    fHampSubRap->GetXaxis()->SetTitle("(a2-a1)/(a1-a0)");
    fHampSubRap->SetOption("colz");

    // Times from channel amplitudes ratio
    name = fName + "_a0/a2_vs_a1/a2_a2_gt_20";
    fHtimeRatio = new TH2F(name.c_str(),name.c_str(),
			  150, 0., 3., 150, 0., 3.);
    AddHistogram(fHtimeRatio);
    fHtimeRatio->GetXaxis()->SetTitle("a1/a2");
    fHtimeRatio->GetYaxis()->SetTitle("a0/a2");
    fHtimeRatio->SetOption("colz");
  }

  // Channel vs strip amplitude
  name = fName + "_ch_vs_a2";
  fHchvsa2 = new TH2F(name.c_str(),name.c_str(),
		    fNchan, -0.5, fNchan-0.5,
		    400, -5, 995);
  AddHistogram(fHchvsa2);

  if (fExpertHistos) {

    // Channel vs strip amplitude
    name = fName + "_ch_vs_a2ped";
    fHchvsa2ped = new TH2F(name.c_str(),name.c_str(),
		      fNchan, -0.5, fNchan-0.5,
		      400, -400, 400);
    AddHistogram(fHchvsa2ped);

    // Channel vs strip amplitude
    name = fName + "_ch_vs_a2CM";
    fHchvsa2CM = new TH2F(name.c_str(),name.c_str(),
		      fNchan, -0.5, fNchan-0.5,
		      400, -400, 400);
    AddHistogram(fHchvsa2CM);

    // temporary histogram used for a2 spectrum fit for each APV
    fHa2tmpname = fName + "_a2_tmp_for_fit";
    fHa2tmp = new TH1D(fHa2tmpname.c_str(), fHa2tmpname.c_str(), 400, -5, 995);

    // 2D plot where bin value is the a2 mean for each apv (shown in x,y)
    name = fName + "_a2_mean_vs_apv";
    fHapvvsa2mean = new TH2F(name.c_str(),name.c_str(),
                             12, 0., 12., 4, 0., 4.);
    fHapvvsa2mean->SetOption("colz");
    AddHistogram(fHapvvsa2mean);
  }

  // Single channel occupancies
  name = fName + "_occupancies";
  fHoccup = new TH1F_Ref(name.c_str(),name.c_str(),
			   fVch->GetNbins(),
			   fVch->GetMin(), fVch->GetMax(),fRateCounter,true);
  ((TH1F_Ref*)fHoccup)->SetReference(fReferenceDirectory);
  AddHistogram(fHoccup);
  fHoccup->SetStats(false);

  // 2D Hit profile
  name = fName + "_2D_profile";
  fHch2D=new TH2F(name.c_str(),name.c_str(),72,0.,72.,
                72, 0., 72.);
  fHch2D->SetOption("colz");
  fHch2D->GetXaxis()->SetNdivisions(12,kFALSE);
  fHch2D->GetYaxis()->SetNdivisions(9,kFALSE);
  AddHistogram(fHch2D);

  // 2D Hit profile
  name = fName + "_2D_amp";
  fHchamp2D=new TH2F(name.c_str(),name.c_str(),72,0.,72.,
                72, 0., 72.);
  fHchamp2D->SetOption("colz");
  fHchamp2D->GetXaxis()->SetNdivisions(12,kFALSE);
  fHchamp2D->GetYaxis()->SetNdivisions(9,kFALSE);
  AddHistogram(fHchamp2D);

  // a2 average profile
  name = fName + "_a2_avg";
  fHavgamp=new TH1D(name.c_str(),name.c_str(),
			   fVch->GetNbins(),
			   fVch->GetMin(), fVch->GetMax());
  AddHistogram(fHavgamp);

  // a2 sigma profile
  name = fName + "_a2_sigma";
  fHsigamp=new TH1D(name.c_str(),name.c_str(),
			   fVch->GetNbins(),
			   fVch->GetMin(), fVch->GetMax());
  AddHistogram(fHsigamp);

  // CM a2 histo
  name = fName + "_CMa2";
  fHCMa2 = new TH1D(name.c_str(),name.c_str(),
		    fVa2->GetNbins(),
		    fVa2->GetMin(),
		    fVa2->GetMax()); 
  AddHistogram(fHCMa2);

  if (fExpertHistos) {

    // CM a2 histo
    name = fName + "_CMa2apv";
    fHCMa2apv = new TH2D(name.c_str(),name.c_str(), 48, 0., 48.,
		         fVa2->GetNbins(),
		         fVa2->GetMin(),
		         fVa2->GetMax()); 
    fHCMa2apv->SetOption("colz");
    AddHistogram(fHCMa2apv);

    // 2D a2 average profile
    name = fName + "_2D_a2_avg";
    fHavgamp2D=new TH2D(name.c_str(),name.c_str(),72,0.,72.,
                  72, 0., 72.);
    fHavgamp2D->SetOption("colz");
    AddHistogram(fHavgamp2D);

    // 2D a2 sigma profile
    name = fName + "_2D_a2_sigma";
    fHsigamp2D=new TH2D(name.c_str(),name.c_str(),72,0.,72.,
                  72, 0., 72.);
    fHsigamp2D->SetOption("colz");
    AddHistogram(fHsigamp2D);

    // a2CM average profile
    name = fName + "_a2CM_avg";
    fHavgampCM=new TH1D(name.c_str(),name.c_str(),
			     fVch->GetNbins(),
			     fVch->GetMin(), fVch->GetMax());
    AddHistogram(fHavgampCM);

    // a2CM sigma profile
    name = fName + "_a2CM_sigma";
    fHsigampCM=new TH1D(name.c_str(),name.c_str(),
			     fVch->GetNbins(),
			     fVch->GetMin(), fVch->GetMax());
    AddHistogram(fHsigampCM);

    // 2D a2CM average profile
    name = fName + "_2D_a2CM_avg";
    fHavgampCM2D=new TH2D(name.c_str(),name.c_str(),72,0.,72.,
                  72, 0., 72.);
    fHavgampCM2D->SetOption("col");
    AddHistogram(fHavgampCM2D);

    // 2D a2CM sigma profile
    name = fName + "_2D_a2CM_sigma";
    fHsigampCM2D=new TH2D(name.c_str(),name.c_str(),72,0.,72.,
                  72, 0., 72.);
    fHsigampCM2D->SetOption("col");
    AddHistogram(fHsigampCM2D);
  }

  // 2D event with a2 CM amplitude
  name = fName + "_2D_a2_event";
  fHevtamp2D=new TH2D(name.c_str(),name.c_str(),72,0.,72.,
                72, 0., 72.);
  fHevtamp2D->SetOption("col");
  AddHistogram(fHevtamp2D);

  // 2D event with a2 CM amplitude with cut
  name = fName + "_2D_a2_cut_event";
  fHevtampcut2D=new TH2D(name.c_str(),name.c_str(),72,0.,72.,
                72, 0., 72.);
  fHevtampcut2D->SetOption("col");
  AddHistogram(fHevtampcut2D);
  
  if(tree) {
    fIsInTree = true;
    std::string hitsname = fName + "_hitMultiplicity";
    std::string xpxname = fName + fVxpx->GetName();
    std::string xpxleavlist = xpxname + "[" + hitsname +"]/F";
    tree->Branch(xpxname.c_str(),fVxpx->GetValues(),
		 xpxleavlist.c_str(),32000);
    std::string ypxname = fName + fVypx->GetName();
    std::string ypxleavlist = ypxname + "[" + hitsname +"]/F";
    tree->Branch(ypxname.c_str(),fVypx->GetValues(),
		 ypxleavlist.c_str(),32000);
  }
    
}

//-----------------------------------------------------------------------------


void PlaneRiAPV::ResetHistograms() {
  PlaneAPV::ResetHistograms();

  for (register int ii =0; ii < MAX_STAT_RiAPV; ii++) {
    fStat[ii].sum = 0.; fStat[ii].sum2 = 0.; fStat[ii].nb = 0;
    fStatCM[ii].sum = 0.; fStatCM[ii].sum2 = 0.; fStatCM[ii].nb = 0;
  }
}

//-----------------------------------------------------------------------------


void PlaneRiAPV::ControlPanel(const TGWindow* p, const TGWindow* main) {

  if (!fControlPanel) fControlPanel = new PlaneRiAPVPanel(p, main, 100, 100, this);
}

//-----------------------------------------------------------------------------


void PlaneRiAPV::StoreDigit(CS::Chip::Digit* digit) {

  lDigits.push_back(digit);
  fNhits++;
}



//-----------------------------------------------------------------------------


#define APV_TIME 25.
#define DELAY_SAMPLES 5

void PlaneRiAPV::EndEvent(const CS::DaqEvent &event) {
  bool sparsefg = true, ampTotCleaned = false;

  if (thr_flag) TThread::Lock();
  register bool rateflag = false;
  if((fRateCounter&0x3f)==0) {
    rateflag = true;
    if (fHevtamp2D) fHevtamp2D->Reset("ICE");
    if (fHevtampcut2D) fHevtampcut2D->Reset("ICE");
  }

  // Get TCS phase
  try {
    tcsphase = event.GetTT().GetPhaseTCS() - 40.; // minus 40. to stay in agreement
                                                  // with previous definition of TCS
                                                  // phase in COOOL
  } catch (std::logic_error& e) {
    std::cerr << "PlaneRiAPV::EndEvent: TCS phase not available: " << e.what() << std::endl;
    tcsphase = 0;
  }

  if ( ! sparsefg) {
    memset(&ampTot, 0, sizeof(ampTot));
  }

  fHhits->Fill(fNhits);

  // Get digits from decoding lib
  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
    bool connectfg = true;
    CS::ChipAPVRICH::Digit* idgt = dynamic_cast<CS::ChipAPVRICH::Digit*> (*ii);
    if (!idgt) {
      std::cerr<<"PlaneRiAPV::EndEvent: a digit is not a APVRICH one, strange...\n";
      (*ii)->Print();
      continue;
    }

    if ( (idgt->GetPixelX() == -1) && (idgt->GetPixelY() == -1) ) { connectfg = false; }


//     int channel=digit.ch;
//     int a0=(int)digit.dt[0];
//     int a1=(int)digit.dt[1];
//     int a2=(int)digit.dt[2];
//    int channel=idgt->GetChannel();
    fOldAPVBoard = idgt->GetOldAPVBoard();
    if (idgt->GetCardRotation() == 's' || idgt->GetCardRotation() == 'S') {
      fCardRotated = true;
    }
    const CS::ChipAPV::DataID dID(idgt->GetDataID());
//    int channel = dID.u.s.chan + 0x80*(dID.u.s.chip_id&3) + 0x200*(idgt->GetPDx());
    int channel = idgt->GetChipChannel() + 0x80*(dID.u.s.chip_id&3) + 0x200*(idgt->GetPDx());
    int a0 = idgt->GetAmplitude()[0];
    int a1 = idgt->GetAmplitude()[1];
    int a2 = idgt->GetAmplitude()[2];
    sparsefg = idgt->IsSparsifed();

// fprintf(stderr, "%s:  ch 0x%x   a0 0x%x   a1 0x%x   a2 0x%x\n", GetName(), channel, a0, a1, a2);
// fprintf(stderr, "       xpx %d   ypx %d    xpd %d   ypd %d \n", idgt->GetPixelX(), idgt->GetPixelY(), idgt->GetPDx(), idgt->GetPDy());
// fprintf(stderr, "       chipch %d 0x%x    chip %d         addr %d 0x%x \n", idgt->GetChipChannel(), idgt->GetChipChannel(), idgt->GetChip(), idgt->GetAddress(), idgt->GetAddress());
//     int t =(int)digit.dt[3];
    float a12 = (a2==0) ? 0 : float(a1)/float(a2);
    float a02 = (a2==0) ? 0 : float(a0)/float(a2);
//    register double t01 = (a1==a0) ? 0 : DELAY_SAMPLES*APV_TIME * (1 - (a1/(double(a1 - a0))));
//    register double t12 = (a2==a1) ? 0 : DELAY_SAMPLES*APV_TIME * (1 - (a1/(double(a2 - a1))));
//     register double t = (t01 + t12)/2. - tcsphase;
//    register double t = (t01 + t12)/2.;
    register double t = a12 - 0.0064*tcsphase;


// cerr<<"PlaneRiAPV::EndEvent: channel "<<channel<<" a0 "<<a0<<" a1 "<<a1<<" a2 "<<a2<<endl;
    if(fVch->Test(channel) &&
       fVa0->Test(a0) &&
       fVa1->Test(a1) &&
       fVa2->Test(a2)) {


      register int xpx = idgt->GetPixelX();
      register int ypx = idgt->GetPixelY();

      fVch->Store(channel);
      fVa0->Store(a0);
      fVa1->Store(a1);
      fVa2->Store(a2);
      fVt->Store(t);
      fVxpx->Store(xpx);
      fVypx->Store(ypx);


      fHch->Fill(channel);
      fHa0->Fill(a0);
      fHa1->Fill(a1);
      fHa2->Fill(a2);
      if (fHa0s) fHa0s->Fill(event.GetBurstNumber(),a0);
      if (fHa1s) fHa1s->Fill(event.GetBurstNumber(),a1);
      if (fHa2s) fHa2s->Fill(event.GetBurstNumber(),a2);
      fHt->Fill(t);
      fHchvsa2->Fill(channel, a2);
      if (a0<5 && a1>2 && a2>20) {
	fHa1a2Ratio->Fill(a12);
	fHtvsTCSph->Fill(tcsphase, a12);
      }
      if (fHampRatio) fHampRatio->Fill(a12,a02);
      if (fHampSubRatio) fHampSubRatio->Fill(a1-a0,a2-a1);
      if (fHampSubRap) if (a1-a0) fHampSubRap->Fill((a2-a1)*1./(a1-a0));
      if (fHtimeRatio) if (a2>20) fHtimeRatio->Fill(a12,a02);
      fNhitsKept++;
      if (rateflag) {
	fHoccup->Reset("ICE");
	fHoccup->Add(fHch, 1/double(fRateCounter));
      }
//      register int xpd = idgt->GetPDx();
//      register int ypd = idgt->GetPDy();
      if ( connectfg ) {
	if (a2 > RiAPV_A2_CUT) fHch2D->Fill(xpx, ypx);
	fHchamp2D->Fill(xpx, ypx, a2);
      }

      fStat[channel].sum += a2;
      fStat[channel].sum2 += a2*a2;
      fStat[channel].nb++;
      fStat[channel].xpx = xpx;
      fStat[channel].ypx = ypx;

      if (rateflag && xpx != -1) {
	fHevtamp2D->Fill(xpx, ypx, a2);
	if (a2 > RiAPV_A2_CUT) fHevtampcut2D->Fill(xpx, ypx, a2);
      }

      if ( ! sparsefg) {
        if (!ampTotCleaned) {
	  memset(&ampTot, 0, sizeof(ampTot));
	  ampTotCleaned = true;
	}
	register int apvnb = (channel & ~0x7f) >> 7;
	register int chnb = channel & 0x7f;
	ampTot[apvnb].ampSort[0].amp[chnb] = a0;
	ampTot[apvnb].ampSort[1].amp[chnb] = a1;
	ampTot[apvnb].ampSort[2].amp[chnb] = a2;
	ampTot[apvnb].ampFlag[chnb] = 1;
      }
    }
  }
  fHhits->Fill(fNhitsKept);


//if (fUseCalib) {
//     register int apvnb, chnb;
//     for (register int ihit = 0; ihit < fVch->GetNvalues(); ihit++) {
//       apvnb = (int) fVch->GetValues()[ihit]; apvnb &= ~0x7f; apvnb >>= 7;
//       chnb = ((int) fVch->GetValues()[ihit]) & 0x7f;
//       ampTot[apvnb].ampSort[0].amp[chnb] = fVa0->GetValues()[ihit];
//       ampTot[apvnb].ampSort[1].amp[chnb] = fVa1->GetValues()[ihit];
//       ampTot[apvnb].ampSort[2].amp[chnb] = fVa2->GetValues()[ihit];
//       ampTot[apvnb].ampFlag[chnb] = 1;
//     }

  if ( ! sparsefg) {
    calcCommonModeCorrection(MAX_AMP_DEVIATION);

    for (register int ch = 0; ch < fNchan; ch++) {
      register int apvnb = (ch & ~0x7f) >> 7;
      register int chnb = ch & 0x7f;
      register double ped;

      if (ampTot[apvnb].ampFlag[chnb] != 1) continue;

      if (fUseCalib) {
#if USE_DATABASE==1
        ped = calib_arr[ch].pedestal;
#else
        ped = 0;
#endif
      } else {
        ped = 0;
      }

//       register double a0CM = -(fVa0->GetValues()[ch] - fCMapv[apvnb][0] - ped);
//       register double a1CM = -(fVa1->GetValues()[ch] - fCMapv[apvnb][1] - ped);
//       register double a2CM = -(fVa2->GetValues()[ch] - fCMapv[apvnb][2] - ped);
      register double a0CM = -(ampTot[apvnb].ampSort[0].amp[chnb] - fCMapv[apvnb][0] - ped);
      register double a1CM = -(ampTot[apvnb].ampSort[1].amp[chnb] - fCMapv[apvnb][1] - ped);
      register double a2CM = -(ampTot[apvnb].ampSort[2].amp[chnb] - fCMapv[apvnb][2] - ped);

      if (fHa0CM) fHa0CM->Fill(a0CM);
      if (fHa1CM) fHa1CM->Fill(a1CM);
      if (fHa2CM) fHa2CM->Fill(a2CM);
//	fHa0CM->Fill(ampTot[apvnb].ampSort[2].amp[chnb]);
//	fHa1CM->Fill(ampTot[apvnb].ampSort[2].amp[chnb] - ped);
      if (fHa0CMs) fHa0CMs->Fill(event.GetBurstNumber(),a0CM);
      if (fHa1CMs) fHa1CMs->Fill(event.GetBurstNumber(),a1CM);
      if (fHa2CMs) fHa2CMs->Fill(event.GetBurstNumber(),a2CM);
      if (fHchvsa2ped) fHchvsa2ped->Fill(ch, ampTot[apvnb].ampSort[2].amp[chnb]-ped);
      if (fHchvsa2CM) fHchvsa2CM->Fill(ch, a2CM);

      if (fHCMa2apv) fHCMa2apv->Fill(apvnb, fCMapv[apvnb][2]);
      fHCMa2->Fill(fCMapv[apvnb][2]);
      fStatCM[ch].sum += a2CM;
      fStatCM[ch].sum2 += a2CM*a2CM;
      fStatCM[ch].nb++;
    }
  }
//}


  if (rateflag) {
    for (register int ii=0; ii<MAX_STAT_RiAPV; ii++) {
      register int nb = fStat[ii].nb;
      if (nb == 0) continue;
      register int nbCM = fStatCM[ii].nb;
//      register int apvnb = (ii & ~0x7f) >> 7;
      register short xpx = fStat[ii].xpx;
      register short ypx = fStat[ii].ypx;
      register double avg = fStat[ii].sum / nb;
      register double sigma = sqrt(fStat[ii].sum2/nb - fStat[ii].sum*fStat[ii].sum/nb/nb);
      fHavgamp->SetBinContent(fHavgamp->GetXaxis()->FindBin(ii),avg);
      fHsigamp->SetBinContent(fHsigamp->GetXaxis()->FindBin(ii),sigma);
      if( xpx != -1 ) {
	if (fHavgamp2D) fHavgamp2D->SetBinContent(fHavgamp2D->GetXaxis()->FindBin(xpx),
                        	  fHavgamp2D->GetYaxis()->FindBin(ypx),avg);
	if (fHsigamp2D) fHsigamp2D->SetBinContent(fHsigamp2D->GetXaxis()->FindBin(xpx),
                        	  fHsigamp2D->GetYaxis()->FindBin(ypx),sigma);
      }
      if (nbCM == 0) continue;
      register double avgCM = fStatCM[ii].sum / nbCM;
      register double sigmaCM = sqrt(fStatCM[ii].sum2/nbCM - fStatCM[ii].sum*fStatCM[ii].sum/nbCM/nbCM);
      if (fHavgampCM) fHavgampCM->SetBinContent(fHavgampCM->GetXaxis()->FindBin(ii),avgCM);
      if (fHsigampCM) fHsigampCM->SetBinContent(fHsigampCM->GetXaxis()->FindBin(ii),sigmaCM);
      if (fHavgampCM2D) fHavgampCM2D->SetBinContent(fHavgampCM2D->GetXaxis()->FindBin(xpx),
                                fHavgampCM2D->GetYaxis()->FindBin(ypx),avgCM);
      if (fHsigampCM2D) fHsigampCM2D->SetBinContent(fHsigampCM2D->GetXaxis()->FindBin(xpx),
                                fHsigampCM2D->GetYaxis()->FindBin(ypx),sigmaCM);
    }

    if (fExpertHistos) {
      for (register int iapv = 0; iapv < fNchan/128; iapv++) {
	TH1D* htmp = fHchvsa2->ProjectionY(fHa2tmpname.c_str(), iapv*128 + 10 + 1, iapv*128 + 118, "");
	htmp->Fit("expo", "Q0", "", RiAPV_A2_FIT_MIN, RiAPV_A2_FIT_MAX);
	TF1* fit1 = htmp->GetFunction("expo");
	register double pente = fit1->GetParameter(1);
	if (pente < 0) {
	  register int apvx = (iapv>>2) + 1;
	  register int apvy = (fCardRotated) ? 4 - (iapv&3) : (iapv&3) + 1;
	  if (fHapvvsa2mean) fHapvvsa2mean->SetBinContent(apvx, apvy, -1./pente);
	}
      }
    }
  }



  if (thr_flag) TThread::UnLock();
}

//-----------------------------------------------------------------------------




void PlaneRiAPV::ChannelA2Spectrum() {
  std::string hdname = fName + "_a2ch";
  if(fCurChan<fNchan) {
    fHchvsa2->ProjectionY(hdname.c_str(),fCurChan+1,fCurChan+1,"")->Draw();
  }
}

//-----------------------------------------------------------------------------
  
//-----------------------------------------------------------------------------
// Calculate common mode noise correction
// (from gemMonitor written by B. Ketzer)
//-----------------------------------------------------------------------------
void PlaneRiAPV::calcCommonModeCorrection(double cut) {

  // Local array for sorted amplitudes of active channels
  ampS3 ampAPV;

  for (register int iapv = 0; iapv < fNchan/128; iapv++) {
    // Loop over channels
    int nb_activech = 0;
    register int iactch = 0;
    memset(&ampAPV, 0, sizeof(ampAPV));
    for ( register int ich = 0; ich < MAX_APV_NBCH; ich++ ) {

      register int ch = ich + iapv*MAX_APV_NBCH;
      // Consider only active channels
#if USE_DATABASE == 1
      if ( fUseCalib ) {
	if ( calib_arr[ch].flag == 1 ) {

          // Copy active channels into new array for sorting
//          ampSort[0].amp[iactch] = fVa0->GetValues()[ch] - calib_arr[ch].pedestal;
//          ampSort[1].amp[iactch] = fVa1->GetValues()[ch] - calib_arr[ch].pedestal;
//          ampSort[2].amp[iactch] = fVa2->GetValues()[ch] - calib_arr[ch].pedestal;
          ampAPV.ampSort[0].amp[iactch] = ampTot[iapv].ampSort[0].amp[ich] - calib_arr[ch].pedestal;
          ampAPV.ampSort[1].amp[iactch] = ampTot[iapv].ampSort[1].amp[ich] - calib_arr[ch].pedestal;
          ampAPV.ampSort[2].amp[iactch] = ampTot[iapv].ampSort[2].amp[ich] - calib_arr[ch].pedestal;
          iactch++;
	}
      } else {

          // Copy active channels into new array for sorting
          ampAPV.ampSort[0].amp[iactch] = ampTot[iapv].ampSort[0].amp[ich];
          ampAPV.ampSort[1].amp[iactch] = ampTot[iapv].ampSort[1].amp[ich];
          ampAPV.ampSort[2].amp[iactch] = ampTot[iapv].ampSort[2].amp[ich];
          iactch++;
      }
#endif
    }
    nb_activech = iactch;

    for (register int ismp=0; ismp<3; ismp++) {

      // Sort in ascending order
      double *ampEnd = ampAPV.ampSort[ismp].amp + nb_activech;
      std::sort(ampAPV.ampSort[ismp].amp, ampEnd);
      //	copy(ampSort, ampEnd, ostream_iterator<double>(cout, " "));

      // Calculate median
      register double median;
      register int ctr;
      if ( nb_activech%2 == 1 ) {
        ctr = (nb_activech-1)/2;
        median = ampAPV.ampSort[ismp].amp[ctr];
      }
      else {
        ctr = nb_activech/2;
        median = (ampAPV.ampSort[ismp].amp[ctr-1]+ampAPV.ampSort[ismp].amp[ctr])/2.;
      }

      // Loop over active channels 
      register double ampMean = 0;
      register int n = 0;
      for ( register int i=0; i<nb_activech; i++ ) {

        // Sum up channels with amplitudes around median
        if ( ((median-ampAPV.ampSort[ismp].amp[i]) < cut) && 
	     ((ampAPV.ampSort[ismp].amp[i]-ampAPV.ampSort[ismp].amp[ctr+i-n]) < cut) ) {
          ampMean += ampAPV.ampSort[ismp].amp[i];
          n++;
        }
      }

      // Calculate mean of channels without signal for CM correction
      fCMapv[iapv][ismp] = ampMean / n;
    }
  }
  
}

//-----------------------------------------------------------------------------


#if USE_DATABASE == 1

void PlaneRiAPV::ReadCalib(const tm &t)
{
  //  Read calibration data for each strip
  try{
    ReadFromDataBase(calib_data, t);
    // std::cout<<"PlaneRiAPV::ReadCalib() ==> "
    // <<this->GetName()<<" calibrations are found !"<<std::endl;

// cerr<<"PlaneRiAPV::ReadCalib: "<<GetName()<<" calib_data.size "<<calib_data.size()<<" fNchan "<<fNchan<<endl;

    if(calib_data.size()*128 != (unsigned) fNchan) {
      std::cerr<<GetName()<<": Size of Calibration File is not correct ! Should be : "
	  <<fNchan<<" Is "<<calib_data.size()*128<<" \n date:"
	  <<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	  <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec<<std::endl;
    }
    else
      fUseCalib = true;
  }

  catch(CS::Exception& e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(const std::exception &e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(...) {
    std::cout<<"PlaneRiAPV::ReadCalib() ==> "<<GetName()
	<<" calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	<<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	<<", not found in DB"<<std::endl;
  }

#ifndef __CINT__
  if (fUseCalib) {

    calib_arr     = new calib_s[fNchan];

    // need apvRich.xml map to sort calibration data
    typedef CS::Chip::Maps::const_iterator m_it;
    typedef std::map<APVref, APVcalib>::const_iterator map_APVcal;

    std::set<uint16> detsrcIDs;
    fRunMaps->GetSrcIDs(GetName(), detsrcIDs);
    for (std::set<uint16>::const_iterator its = detsrcIDs.begin(); its != detsrcIDs.end(); its++) {
      register unsigned int srcId = *its;

      if (srcId) {
	m_it m_low = fRunMaps->lower_bound(uint64(srcId)<<48);
        // fRunMaps->Print(std::cerr, "");
	for( m_it cc=m_low; ((*cc).first>>48)==srcId; cc++ ) {
          const CS::ChipAPVRICH::Digit *m =
	    dynamic_cast<const CS::ChipAPVRICH::Digit*>((*cc).second);
          if( m==NULL ) {
	    std::cerr << "ChipAPVRICH wrong map.\n";
	    continue;
          }
          if (strcmp(m->GetDetID().GetName().c_str(), GetName()))
	    continue;


          CS::ChipAPVRICH::DataID dataid(m->GetDataID());
          const APVcalib& a_cal = calib_data[APVref(0, dataid.u.s.src_id, dataid.u.s.adc_id, dataid.u.s.chip_id)];
          if (m->GetChipChannel() == 0) for (register int ch=0; ch<MAX_APV_NBCH; ch++) {
	    register int strip = ch + 0x80*(dataid.u.s.chip_id&3) + 0x200*(m->GetPDx());
            if (&(a_cal.channel[ch])) {
              const APVchannel& a_ch = a_cal.channel[ch];
              calib_arr[strip].flag     = a_ch.flag;
              calib_arr[strip].pedestal = a_ch.ped/2.;    // values in calib files multiplied by 2
              calib_arr[strip].sigma    = a_ch.sigma/2.;  // values in calib files multiplied by 2
            } else {
              calib_arr[strip].flag     = 0;
              calib_arr[strip].pedestal = 0;
              calib_arr[strip].sigma    = 0;
            }
	  }
	}
      }
    }


    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_RA")!=0) {
      std::cout << "Calibration values for RichAPV "<<fName<<":\n";
      std::cout << "  ch  flag   ped   sigma\n";
      for (register int ii=0; ii<fNchan; ii++) {
        std::cout<<"  "<<ii<<"  "<<calib_arr[ii].flag<<"   "<<calib_arr[ii].pedestal<<"   "<<calib_arr[ii].sigma<<std::endl;
      }
    }
  }  // fUseCalib
#endif


//     // assume same source id for all apvs/channels from one RiAPV
//     unsigned int srcId = 0;
//     for (map_APVcal ical = calib_data.begin(); ical != calib_data.end(); ical++) {
//       const APVcalib& cal0 = (*ical).second;
//       if (&(cal0.channel[0])) {
//         //const APVchannel& ch0 = cal0.channel[0];
//         srcId = cal0.ref.srcId;
//         break;
//       }
//     }

}

//   // Read calibration data for conversion of amplitude->time for each plane
//   if (fDataBase==NULL)
//     throw CS::Exception("PlaneRiAPV::ReadCalib():  data base is not opened.");
//   try{
//     //    std::cout << "Reading timing calibrations for " << GetName() << std::endl;
//     struct tm tt(t);
//     CDB::Time tp(mktime(&tt),0);
//     std::string strdata("");
//     fDataBase->read(fName,strdata,tp,"timing");
//     if (strdata == "") {
//       std::cerr << "PlaneRiAPV::ReadCalib() "<<GetName()<<" : no timing calibration file found"<<std::endl;
//       return;
//     }
// # if __GNUC__ > 3 || __GNUC__ == 3
//     std::istringstream istrdata(strdata.c_str());
// #else
//     istrstream istrdata(strdata.c_str());
// #endif
//     //    std::cout << strdata << std::endl;
//     istrdata >> calib_time;
//   }
//   catch(CS::Exception& e) {
//     std::cout<<"rethrowing"<<std::endl;
//     throw;
//   }
// 
//   if (calib_time.size()!=0) {
// 
//     // Just debug printouts
//     if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_RA")!=0) {
//       const float& r1_min = this->calib_time[0];
//       const float& r1_max = this->calib_time[1];
//       const float& t0_1 = this->calib_time[2];
//       const float& sl_1 = this->calib_time[3];
//       const float& tc_1 = this->calib_time[4];
//       const float& r2_min = this->calib_time[5];
//       const float& r2_max = this->calib_time[6];
//       const float& t0_2 = this->calib_time[7];
//       const float& sl_2 = this->calib_time[8];
//       const float& tc_2 = this->calib_time[9];
//       std::cout<<std::endl<<"Timing calibrations read in for "<<GetName()<<std::endl;
//       std::cout<<"r1_min = "<<r1_min<<std::endl;
//       std::cout<<"r1_max = "<<r1_max<<std::endl;
//       std::cout<<"r2_min = "<<r2_min<<std::endl;
//       std::cout<<"r2_max = "<<r2_max<<std::endl;
//       std::cout<<"t0_1   = "<<t0_1<<std::endl;
//       std::cout<<"t0_2   = "<<t0_2<<std::endl;
//       std::cout<<"sl_1   = "<<sl_1<<std::endl;
//       std::cout<<"sl_2   = "<<sl_2<<std::endl;
//       std::cout<<"tc_1   = "<<tc_1<<std::endl;
//       std::cout<<"tc_2   = "<<tc_2<<std::endl;
//       std::cout<<std::endl;
//     }
//   } else {
//     std::cout<<"PlaneRiAPV::ReadCalib() ==> "<<GetName()
// 	<<" timing calibrations, valid for ";
//     std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
// 	<<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
// 	<<", not found in DB"<<std::endl;
//   }

#endif //USE_DATABASE




void PlaneRiAPV::TextOutput(ostream& out)
/// get of means in hits histograms by Damien, july03
{
  out << "\thitmean "<<fHhits->GetMean(1)<<std::endl;

  fHa2->Fit("expo", "Q0", "", RiAPV_A2_FIT_MIN, RiAPV_A2_FIT_MAX);
  TF1* fit1 = fHa2->GetFunction("expo");
  if (fit1 != NULL) {
    Double_t pente = fit1->GetParameter(1);
    out << "\tpente "<<pente<<std::endl;
  } else {
    out << "\tpente ***"<<std::endl;
  }
  if (fExpertHistos) {
    for (register int iapv = 0; iapv < fNchan/128; iapv++) {
      TH1D* htmp = fHchvsa2->ProjectionY(fHa2tmpname.c_str(), iapv*128 + 10 + 1, iapv*128 + 118, "");
      htmp->Fit("expo", "Q0", "", RiAPV_A2_FIT_MIN, RiAPV_A2_FIT_MAX);
      TF1* fit1 = htmp->GetFunction("expo");
      if (fit1 != NULL) {
        double penteapv = fit1->GetParameter(1);
        out << "\tpenteapv  "<<iapv<<" "<<penteapv<<std::endl;
      } else {
        out << "\tpenteapv  "<<iapv<<" ***"<<std::endl;
      }
    }
  }
}



//-----------------------------------------------------------------------------

#if USE_DATABASE == 1


istream& operator>>(istream& in, std::map<const PlaneRiAPV::APVref, PlaneRiAPV::APVcalib> &c) {
  char firstChar;
  char chipname[10];
  int ldcId, srcId, adcId, chipId;
  PlaneRiAPV::APVref* rf = 0;
  c.clear();
//  for(int i=0; i<12; i++ ) c.push_back(PlaneRiAPV::APVcalib());


  while (in) {
    in >> firstChar;
    if (firstChar == '#') {
      in >> chipname >> ldcId >> srcId >> adcId  >> chipId;
      ldcId = 0;  // to avoid difficulties to get this value in data
      //std::cout << "Reading in "<<chipname<<" "<<ldcId<<" "<<srcId
      //   <<" "<<adcId<<" "<<chipId<<std::endl;
      rf = new PlaneRiAPV::APVref(ldcId, srcId, adcId, chipId);
      c.insert(std::pair<const PlaneRiAPV::APVref, const PlaneRiAPV::APVcalib>(*rf, PlaneRiAPV::APVcalib(*rf)));
    }
    else {
      int f = atoi(&firstChar);
      float p, s, cp, cs;
      in >> p >> s >> cp >> cs;
      if (rf) {
        c[*rf].channel.push_back(PlaneRiAPV::APVchannel(f, p, s, cp, cs));
      } else {
        throw "PlaneRiAPV calib: bad calibration file\n";
      }
    }
  }

  return in;
}


#endif //USE_DATABASE
