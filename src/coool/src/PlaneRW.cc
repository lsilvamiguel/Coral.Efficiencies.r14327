#include <cmath>
#include <string.h>
#include <TTree.h>
#include <TThread.h>
#include "PlaneRW.h"
#include "ChipF1.h"
#include "Plane1VPanel.h"

ClassImp(PlaneRichWall);

const float f1_reso_RW = 0.123063;
const float TIME_SHIFT = 2500; // in ns
const int nPlanes = 8;

struct RW_geometryStruct
{
  char fName[16];
  int  fNChannels;
};

const RW_geometryStruct RW_geometry[nPlanes] =
  {
    {"DR01X1", 32},
    {"DR01X2", 32},
    {"DR02X1", 32},
    {"DR02X2", 32},
    {"DR01Y1", 32},
    {"DR01Y2", 32},
    {"DR02Y1", 32},
    {"DR02Y2", 32}
  };

PlaneRichWall::PlaneRichWall(const char *detname,int nchan, int center, int width)
  :Plane1V(detname,nchan,center,width)
{
  Nevents =0;
  fPlaneNumber = -1;
  for(int i = 0; i < nPlanes; i++)
    {
      if(strncmp(RW_geometry[i].fName,detname,strlen(RW_geometry[i].fName)) == 0)
	{
	  fPlaneNumber = i;
	}
    }
  if(fPlaneNumber == -1)
    {
      std::cerr << "PlaneRichWall::PlaneRichWall: Incorrect detector name: " << detname << std::endl;
      throw CS::Exception("PlaneRichWall::PlaneRichWall: Incorrect detector name");
    }
}

void PlaneRichWall::Init(TTree *tree) 
{

  Nevents =0;

  fRateCounter = 0;

  //  Plane1V::Init(tree);

  
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

  if(tree) {
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

  //---------------------------------------------------------------------------------
  // Rates HISTOGRAMS
  std::string rname2 = fName + "_rates_2";
  fHrates2=new TH1F_Ref(rname2.c_str(),rname2.c_str(),
		   fVch->GetNbins(),
		   fVch->GetMin(),
		   fVch->GetMax(), fRateCounter, true);
  fHrates2->SetYTitle("rates per channel (kHz) test");
  fHrates2->SetTitleOffset(1.2,"Y");
  fHrates2->SetOption("hist");
  AddHistogram(fHrates2);
  ((TH1F_Ref*)fHrates2)->SetReference(fReferenceDirectory);
  

  //branch : hit map
  std::string chname4 = fName + fVch->GetName() + "_2";
   fHch2=new TH1F_Ref(chname4.c_str(),chname4.c_str(),
 		    fVch->GetNbins(),
 		    fVch->GetMin(),
 		    fVch->GetMax(),
 		    fRateCounter);
//   fHch2=new TH1F(chname.c_str(),chname.c_str(),
// 		    fVch->GetNbins(),
// 		    fVch->GetMin(),
// 		    fVch->GetMax());
  ((TH1F_Ref*)fHch2)->SetReference(fReferenceDirectory);
  AddHistogram(fHch2);

  std::string chname2 = fName + fVch->GetName() + "_noise%";
  fHnoise_percent=new TH1F(chname2.c_str(),chname2.c_str(),
		    fVch->GetNbins(),
		    fVch->GetMin(),
		    fVch->GetMax());
  AddHistogram(fHnoise_percent);

  std::string chname3 = fName + fVch->GetName() + "_time2";
  fHtime_rev=new TH1F_Ref(chname3.c_str(),chname3.c_str(),
		    3001,
		    -3000,
			  0,fRateCounter);
  ((TH1F_Ref*)fHtime_rev)->SetReference(fReferenceDirectory);
  AddHistogram(fHtime_rev);

  
}

#ifndef __CINT__
void PlaneRichWall::EndEvent(const CS::DaqEvent &event)
{

  if (thr_flag) TThread::Lock();

  if( 7 == event.GetHeader().GetEventType() ){
    
    ++Nevents;
    
    
    typedef std::list<CS::Chip::Digit*>::iterator lDIter;
    for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
      CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
      if (!iii) {
	std::cerr<<"Plane1V::EndEvent: a digit is not a F1 one, strange...\n";
	continue;
      }
      double time = (double)iii->GetTimeDecoded();
      int channel = iii->GetChannel();
      
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
	
	//----------------------------------
	fHch2->Fill(channel);	
	fHtime_rev->Fill(time);
	
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
      //  float timewin = ( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() )*fF1_TICK*fRateCounter;
      float timewin = ( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() )*fRateCounter;
      if( timewin > 0 ){
	for(register int i=1; i<=fHonTP->GetNbinsX(); i++) {
	  fHrates->SetBinContent(i,((float) fHoffTP->GetBinContent(i)/timewin)*1000000);
	  fHrates->SetBinError(i,((float) fHoffTP->GetBinError(i)/timewin)*1000000);
	}
      }
      
    }
    //---------------------------------------------
    
    for(int i=1; i<=fHch2->GetNbinsX(); i++) {
      float rate = fHch2->GetBinContent(i) * 2857.143/Nevents;
      //     if(fPlaneNumber==1 && fHch2->GetBinContent(i) != 0)
      //       printf("Nevents=%5d, %d = %lf, %10.3f\n",Nevents,i,fHch2->GetBinContent(i),rate);
      fHrates2->SetBinContent(i, rate );
      fHrates2->SetBinError(i,  sqrt(rate) );
      
      float noise = fHch2->GetBinContent(i)*100./Nevents ;
      fHnoise_percent->SetBinContent(i, noise );

    }
    
  }
  
  if (thr_flag) TThread::UnLock();
  
  // Plane1V::EndEvent(event);
}

void PlaneRichWall::Reset( void ) {
  Plane1V::Reset();
}

void PlaneRichWall::ResetHistograms() {
  Nevents = 0;
  Plane1V::ResetHistograms();
}


#endif
