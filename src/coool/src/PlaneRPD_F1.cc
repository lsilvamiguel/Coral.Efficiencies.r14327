#include "PlaneRPD_F1.h"
#include "PlaneRPD_F1_Panel.h"
#include "ChipF1.h"
#include "TriggerTime.h"
#include "TFile.h"
#include "TWebFile.h"
#include "Reference.h"

#include "TF1.h"

  //new from here
typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  // to here

ClassImp(PlaneRPD_F1);

PlaneRPD_F1::PlaneRPD_F1(const char *detname,int nchan, int center, int width)
  : Plane1V(detname,nchan, -16000, 1000)
{
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_PLANE_RPD_F1; ++iChannel)
    {
      for(unsigned short int iMult=0; iMult<F1_MULTIPLICITY; ++iMult)
	{
	  fDataSent[iChannel][iMult]=false;
	  fTdcRpd[iChannel][iMult]=0;
	}
    }
}

PlaneRPD_F1::~PlaneRPD_F1()
{
//   h1NumberOfSentData->~TH1F();
//   h2TdcVsChannel->~TH2F();
//   h2MultiplicityVsChannel->~TH2F();
}

void PlaneRPD_F1::Init(TTree *tree)
{
  //new from here
  Plane1V::Init(tree);
  // to here
  
  std::string histname = fName + "_NbOfData";
  std::string histtitle = fName + " Number of channels which sent data";
  h1NumberOfSentData=new TH1F_Ref(histname.c_str(),histtitle.c_str(),N_CHANNEL_PLANE_RPD_F1,0,N_CHANNEL_PLANE_RPD_F1,fRateCounter);
  h1NumberOfSentData->GetXaxis()->SetTitle("Number of sent data");
  AddHistogram(h1NumberOfSentData); // this method is used to put the histo (in argument) in the online monitor window
  
  histname = fName + "_TdcVsChannel";
  histtitle = fName + " Tdc spectrum Vs channels";
  h2TdcVsChannel=new TH2F(histname.c_str(),histtitle.c_str(),N_CHANNEL_PLANE_RPD_F1,0,N_CHANNEL_PLANE_RPD_F1,100,-17500,-14500);
  h2TdcVsChannel->SetOption("col");
  h2TdcVsChannel->GetXaxis()->SetTitle("Channels");
  h2TdcVsChannel->GetYaxis()->SetTitle("Tdc");
  AddHistogram(h2TdcVsChannel);

  histname = fName + "_ MultVsChannel";
  histtitle = fName + " Multiplicity spectrum Vs channels";
  h2MultiplicityVsChannel=new TH2F(histname.c_str(),histtitle.c_str(),
				   N_CHANNEL_PLANE_RPD_F1,0,N_CHANNEL_PLANE_RPD_F1,F1_MULTIPLICITY+1,0,F1_MULTIPLICITY+1);
  h2MultiplicityVsChannel->SetOption("col");
  h2MultiplicityVsChannel->GetXaxis()->SetTitle("Channels");
  h2MultiplicityVsChannel->GetYaxis()->SetTitle("Tdc");
  AddHistogram(h2MultiplicityVsChannel);
  
  if(tree)
    {

      // Here describe what will be saved in the out ttree

      fIsInTree = true;
      TString planeName = (TString)fName;
      TString branchName = "fTdcRpd[";branchName+=N_CHANNEL_PLANE_RPD_F1;branchName+="][";branchName+=F1_MULTIPLICITY;branchName+="]/D";
      tree->Branch("rawTdc"+planeName,fTdcRpd,branchName);
      branchName = "fDataSent[";branchName+=N_CHANNEL_PLANE_RPD_F1;branchName+="][";branchName+=F1_MULTIPLICITY;branchName+="]/B";
      tree->Branch("rawDataSent"+planeName,fDataSent,branchName);
    }
}

void PlaneRPD_F1::StoreDigit(CS::Chip::Digit* digit)
{
  //  double time=-1.e6;

//   static int counterStore(0);
//   std::cout<<"store "<<counterStore<<"\n";
//   ++counterStore;

//   Plane1V::StoreDigit(digit);
  CS::ChipF1::Digit *f1dig = dynamic_cast<CS::ChipF1::Digit*>(digit); // Here try to convert from a Chip object to a ChipF1 object
  if(f1dig != NULL)  // if succeed, save the digit
    {lDigits.push_back(digit);}
  else // if not, gives an error
    {
      std::cerr<<"PlaneRPD_F1::StoreDigit ("<<GetName()<<"): a digit is not a F1 one, strange...\n";
      return;
    }
}


void PlaneRPD_F1::EndEvent(const CS::DaqEvent &event)
{ 
  Plane1V::EndEvent(event);

  fEventTimeSecond = (double)event.GetTime().first;
  fEventTimeMicroSecond = (double)event.GetTime().second;
  fEventRunNumber = (int)event.GetRunNumber();
  fEventNumberInRun = (int)event.GetEventNumberInRun();

  double tdcOfCurrentHit;  

  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) // loop over digits saves by PlaneRPD_F1::EndEvent(const CS::DaqEvent &event)
    {
      CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii); //convert the digit
      uint16 channel = iii->GetChannel();
      for(unsigned short int iMult=0; iMult<F1_MULTIPLICITY; ++iMult)
	{
	  tdcOfCurrentHit=(double) (iii->GetTimeDecoded() / iii->GetTimeUnit()); // calculate the time from the digit
	  if(fDataSent[channel][iMult]==false )
	    {
	      fDataSent[channel][iMult]=true;
	      fTdcRpd[channel][iMult]=tdcOfCurrentHit; // put the calculated time into the member
	      h2TdcVsChannel->Fill(channel, fTdcRpd[channel][iMult]);
	      break;
	    }
	}
    }
  unsigned int numberOfDataSent=0;
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_PLANE_RPD_F1; ++iChannel)
    {
      if(fDataSent[iChannel][0])++numberOfDataSent;
      unsigned short int multiplicity=0;
      for(unsigned short int iMult=0; iMult<F1_MULTIPLICITY; ++iMult)
	if(fDataSent[iChannel][iMult]==true )
	  {++multiplicity;}  // calculate the multiplicity for each channel
	else
	  break;
      h2MultiplicityVsChannel->Fill(iChannel,multiplicity);
    }
  h1NumberOfSentData->Fill(numberOfDataSent);

  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_PLANE_RPD_F1/2; ++iChannel)
    {
      if(fDataSent[iChannel][0] && fDataSent[iChannel+N_CHANNEL_PLANE_RPD_F1/2][0])
	fTimeRpd[iChannel][0] = fTdcRpd[iChannel][0] - fTdcRpd[iChannel+N_CHANNEL_PLANE_RPD_F1/2][0];
    }  
}

void PlaneRPD_F1::Reset( void )
{
  Plane1V::Reset();
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_PLANE_RPD_F1; ++iChannel)
    {
      for(unsigned short int iMult=0; iMult<F1_MULTIPLICITY; ++iMult)
	{
	  fDataSent[iChannel][iMult]=false;
	  fTdcRpd[iChannel][iMult]=0;
	}
    }
}

void PlaneRPD_F1::TdcSpectrum()
{
  std::string hdname = fName + "_1ch";
  if(fCurChan_tdc<N_CHANNEL_PLANE_RPD_F1) {
    h2TdcVsChannel->ProjectionY(hdname.c_str(),fCurChan_tdc+1,fCurChan_tdc+1,"")->Draw();
  }
}

// void PlaneRPD_F1::MultiplicitySpectrum()
// {
//   std::string hdname = fName + "_1ch";
//   if(fCurChan_multiplicity<N_CHANNEL_PLANE_RPD_F1) {
//     h2MultiplicityVsChannel->ProjectionY(hdname.c_str(),fCurChan_multiplicity+1,fCurChan_multiplicity+1,"")->Draw();
//   }
// }

void PlaneRPD_F1::ControlPanel(const TGWindow* p, const TGWindow* main)
{
  if (!fControlPanel) fControlPanel = new PlaneRPD_F1_Panel(p, main, 100, 100, this);
}

