#include <math.h>
#include <cstdio>
#include "TThread.h"
#include "GroupRPD.h"
#include "DaqEvent.h"
#include "TLine.h"

ClassImp(GroupRPD);

using namespace std;
using namespace CS;
int Modulo(int value, int modulo);

GroupRPD::GroupRPD(const char* name) :
    Group(name),pf1_A(0),pf1_Bl(0),pf1_Bh(0),psadcl(0),psadch(0),fOkPlanes(false)
{
}

void GroupRPD::Init(void)
{
  Group::Init();

// code to recognize RPD planes has been placed on Init instead of Endevent (DN 15/9/2008)
  int ipf1A = -1;
  int ipf1Bl = -1;
  int ipf1Bh = -1;
  int ipsadcl = -1;
  int ipsadch = -1;
  int i_it = -1;
  fOkPlanes = false;
  
  for( vector<const Plane*>::const_iterator pp=fPlanes.begin(); pp!=fPlanes.end(); pp++ )
    {
      ++i_it;
      
      //       cout<<(*pp)->GetName()<<"\n";
      const PlaneRPD_F1 *pf1=dynamic_cast<const PlaneRPD_F1 *>(*pp); // just to verify that it is of the proper type
      if( pf1==NULL )
        {
	  const PlaneRPD_SADC *psadc=dynamic_cast<const PlaneRPD_SADC *>(*pp); // just to verify that it is of the proper type
	  if( psadc==NULL )
	    {continue;}
	  else
	    if((TString)psadc->GetName()=="RP01Ql__") ipsadcl = i_it;
	    else if((TString)psadc->GetName()=="RP01Qh__") ipsadch = i_it;
	    else cout<<"Warning in GroupRPD::Init():: Casting to PlaneRPD_SADC works but name "<<(TString)psadc->GetName()<<" is not forseen\n";
	}
      else
	{
	  if((TString)pf1->GetName()=="RP01TA__") ipf1A = i_it;
	  else if((TString)pf1->GetName()=="RP01TBl_") ipf1Bl = i_it;
	  else if((TString)pf1->GetName()=="RP01TBh_") ipf1Bh = i_it;
	  else cout<<"Warning in GroupRPD::Init():: Casting to PlaneRPD_F1 works but name "<<(TString)pf1->GetName()<<" is not forseen\n";
	}       	
    }

  if(ipf1A==-1 || ipf1Bl==-1 || ipf1Bh==-1 || ipsadcl==-1|| ipsadch==-1)
    {
      cerr<<"GroupRPD::Init(): No PlaneRPD_SADC or PlaneRPD_F1 found!\n";
      cerr<<"These planes have been found : ";
      for( vector<const Plane*>::const_iterator pp=fPlanes.begin(); pp!=fPlanes.end(); pp++ )
	{
	  cerr<<(*pp)->GetName()<<"\t";
	}
      cerr<<"\n";
      return;
    }

  pf1_A=dynamic_cast<const PlaneRPD_F1 *>(fPlanes[ipf1A]);
  pf1_Bl=dynamic_cast<const PlaneRPD_F1 *>(fPlanes[ipf1Bl]);
//  pf1_Bh=dynamic_cast<const PlaneRPD_F1 *>(fPlanes[ipf1Bh]);
  psadcl=dynamic_cast<const PlaneRPD_SADC *>(fPlanes[ipsadcl]);
//  psadch=dynamic_cast<const PlaneRPD_SADC *>(fPlanes[ipsadch]);
  fOkPlanes = true;

  TString histname;
  TString histtitle;
  
  h1Za = new TH1F*[12];
  for(int iScintillatorA=0; iScintillatorA<12; ++iScintillatorA)
    {
      histname = fName + "_Z_A"; histname+=iScintillatorA;
      histtitle = fName + " Z spectrum for A "; histtitle+=iScintillatorA;
      h1Za[iScintillatorA]=new TH1F(histname,histtitle,100,-50,50);
      h1Za[iScintillatorA]->GetXaxis()->SetTitle("Za (cm)");
      AddHistogram(h1Za[iScintillatorA]);
    }

  h1Zb = new TH1F*[24];
  for(int iScintillatorB=0; iScintillatorB<24; ++iScintillatorB)
    {
      histname = fName + "_Z_B"; histname+=iScintillatorB;
      histtitle = fName + " Z spectrum for B "; histtitle+=iScintillatorB;
      h1Zb[iScintillatorB]=new TH1F(histname,histtitle,100,-70,120);
      h1Zb[iScintillatorB]->GetXaxis()->SetTitle("Zb (cm)");
      AddHistogram(h1Zb[iScintillatorB]);
    }

  if(fExpertHistos) {
    h1TupMinusTdownInA = new TH1F*[12];
    for(int iScintillatorA=0; iScintillatorA<12; ++iScintillatorA)
      {
	histname = fName + "_Tup_Tdown_A"; histname+=iScintillatorA;
	histtitle = fName + " (Tup - Tdown) spectrum for A "; histtitle+=iScintillatorA;
	h1TupMinusTdownInA[iScintillatorA]=new TH1F(histname,histtitle,400,-50,50);
	h1TupMinusTdownInA[iScintillatorA]->GetXaxis()->SetTitle("Tup - Tdown (ns)");
	AddHistogram(h1TupMinusTdownInA[iScintillatorA]);
      }

    h1TupMinusTdownInB = new TH1F*[24];
    for(int iScintillatorB=0; iScintillatorB<24; ++iScintillatorB)
      {
	histname = fName + "_Tup_Tdown_B"; histname+=iScintillatorB;
	histtitle = fName + " (Tup - Tdown) spectrum for B "; histtitle+=iScintillatorB;
	h1TupMinusTdownInB[iScintillatorB]=new TH1F(histname,histtitle,400,-50,50);
	h1TupMinusTdownInB[iScintillatorB]->GetXaxis()->SetTitle("Tup - Tdown (ns)");
	AddHistogram(h1TupMinusTdownInB[iScintillatorB]);
      }
  }

  h1Tof = new TH1F*[12*3];
  h2EnergyLostInBVsBeta = new TH2F*[12*3];
  if(fExpertHistos) h2EnergyLostInAVsBeta = new TH2F*[12*3];
  for(int iScintillatorA=0; iScintillatorA<12; ++iScintillatorA)
    {
      histname = fName + "_T_A";
      histname += iScintillatorA;
      histname += "_T_B";
      histname += Modulo((2*iScintillatorA-1),24);
      histtitle = fName + " Tof between  A ";
      histtitle += iScintillatorA;
      histtitle += " and B ";
      histtitle += Modulo((2*iScintillatorA-1),24);
      h1Tof[3*iScintillatorA]=new TH1F(histname,histtitle,200,-80,80);
      h1Tof[3*iScintillatorA]->GetXaxis()->SetTitle("Tof");
//       if(fExpertHistos)
	AddHistogram(h1Tof[3*iScintillatorA]);

      histname = fName + "_T_A";
      histname += iScintillatorA;
      histname += "_T_B";
      histname += Modulo((2*iScintillatorA),24);
      histtitle = fName + " Tof between  A ";
      histtitle += iScintillatorA;
      histtitle += " and B ";
      histtitle += Modulo((2*iScintillatorA),24);
      h1Tof[3*iScintillatorA+1]=new TH1F(histname,histtitle,200,-80,80);
      h1Tof[3*iScintillatorA+1]->GetXaxis()->SetTitle("Tof");
//       if(fExpertHistos)
	AddHistogram(h1Tof[3*iScintillatorA+1]);

      histname = fName + "_T_A";
      histname += iScintillatorA;
      histname += "_T_B";
      histname += Modulo((2*iScintillatorA+1),24);
      histtitle = fName + " Tof between  A ";
      histtitle += iScintillatorA;
      histtitle += " and B ";
      histtitle += Modulo((2*iScintillatorA+1),24);
      h1Tof[3*iScintillatorA+2]=new TH1F(histname,histtitle,200,-80,80);
      h1Tof[3*iScintillatorA+2]->GetXaxis()->SetTitle("Tof");
//       if(fExpertHistos)
	AddHistogram(h1Tof[3*iScintillatorA+2]);
      
      histname = fName + "_EnergyBVsBeta_A";
      histname+=iScintillatorA;
      histname+="B";
      histname+=Modulo((2*iScintillatorA)-1,24);
      histtitle = fName + " Energy lost in B vs beta for A ";
      histtitle += iScintillatorA;
      histtitle += " and B ";
      histtitle += Modulo((2*iScintillatorA)-1,24);
      h2EnergyLostInBVsBeta[3*iScintillatorA]=new TH2F(histname,histtitle, 200,0,1.2, 100,0,4000);
//      h2EnergyLostInBVsBeta[3*iScintillatorA]->SetOption("col");
      h2EnergyLostInBVsBeta[3*iScintillatorA]->SetOption("cont0");
      h2EnergyLostInBVsBeta[3*iScintillatorA]->GetXaxis()->SetTitle("Beta");
      h2EnergyLostInBVsBeta[3*iScintillatorA]->GetYaxis()->SetTitle("Energy");
//       if(fExpertHistos)
	AddHistogram(h2EnergyLostInBVsBeta[3*iScintillatorA]);

      histname = fName + "_EnergyBVsBeta_A";
      histname+=iScintillatorA;
      histname+="B";
      histname+=Modulo((2*iScintillatorA),24);
      histtitle = fName + " Energy lost in B vs beta for A ";
      histtitle += iScintillatorA;
      histtitle += " and B ";
      histtitle += Modulo((2*iScintillatorA),24);
      h2EnergyLostInBVsBeta[3*iScintillatorA+1]=new TH2F(histname,histtitle, 200,0,1.2, 100,0,4000);
      h2EnergyLostInBVsBeta[3*iScintillatorA+1]->SetOption("cont0");
      h2EnergyLostInBVsBeta[3*iScintillatorA+1]->GetXaxis()->SetTitle("Beta");
      h2EnergyLostInBVsBeta[3*iScintillatorA+1]->GetYaxis()->SetTitle("Energy");
//       if(fExpertHistos)
	AddHistogram(h2EnergyLostInBVsBeta[3*iScintillatorA+1]);

      histname = fName + "_EnergyBVsBeta_A";
      histname+=iScintillatorA;
      histname+="B";
      histname+=Modulo((2*iScintillatorA)+1,24);
      histtitle = fName + " Energy lost in B vs beta for A ";
      histtitle += iScintillatorA;
      histtitle += " and B ";
      histtitle += Modulo((2*iScintillatorA)+1,24);
      h2EnergyLostInBVsBeta[3*iScintillatorA+2]=new TH2F(histname,histtitle, 200,0,1.2, 100,0,4000);
      h2EnergyLostInBVsBeta[3*iScintillatorA+2]->SetOption("cont0");
      h2EnergyLostInBVsBeta[3*iScintillatorA+2]->GetXaxis()->SetTitle("Beta");
      h2EnergyLostInBVsBeta[3*iScintillatorA+2]->GetYaxis()->SetTitle("Energy");
//       if(fExpertHistos)
	AddHistogram(h2EnergyLostInBVsBeta[3*iScintillatorA+2]);


      if(fExpertHistos) {
	histname = fName + "_EnergyAVsBeta_A";
	histname+=iScintillatorA;
	histname+="B";
	histname+=Modulo((2*iScintillatorA)-1,24);
	histtitle = fName + " Energy lost in A vs beta for A ";
	histtitle += iScintillatorA;
	histtitle += " and B ";
	histtitle += Modulo((2*iScintillatorA)-1,24);
	h2EnergyLostInAVsBeta[3*iScintillatorA]=new TH2F(histname,histtitle,200,0,1.2, 100,0,4000);
	h2EnergyLostInAVsBeta[3*iScintillatorA]->SetOption("col");
	h2EnergyLostInAVsBeta[3*iScintillatorA]->GetXaxis()->SetTitle("Beta");
	h2EnergyLostInAVsBeta[3*iScintillatorA]->GetYaxis()->SetTitle("Energy");
  //       if(fExpertHistos)
	  AddHistogram(h2EnergyLostInAVsBeta[3*iScintillatorA]);

	histname = fName + "_EnergyAVsBeta_A";
	histname+=iScintillatorA;
	histname+="B";
	histname+=Modulo((2*iScintillatorA),24);
	histtitle = fName + " Energy lost energy vs beta in A for A ";
	histtitle += iScintillatorA;
	histtitle += " and B ";
	histtitle += Modulo((2*iScintillatorA),24);
	h2EnergyLostInAVsBeta[3*iScintillatorA+1]=new TH2F(histname,histtitle, 200,0,1.2, 100,0,4000);
	h2EnergyLostInAVsBeta[3*iScintillatorA+1]->SetOption("col");
	h2EnergyLostInAVsBeta[3*iScintillatorA+1]->GetXaxis()->SetTitle("Beta");
	h2EnergyLostInAVsBeta[3*iScintillatorA+1]->GetYaxis()->SetTitle("Energy");
  //       if(fExpertHistos)
	  AddHistogram(h2EnergyLostInAVsBeta[3*iScintillatorA+1]);

	histname = fName + "_EnergyAVsBeta_A";
	histname+=iScintillatorA;
	histname+="B";
	histname+=Modulo((2*iScintillatorA)+1,24);
	histtitle = fName + " Energy lost energy vs beta in A for A ";
	histtitle += iScintillatorA;
	histtitle += " and B ";
	histtitle += Modulo((2*iScintillatorA)+1,24);
	h2EnergyLostInAVsBeta[3*iScintillatorA+2]=new TH2F(histname,histtitle, 200,0,1.2, 100,0,4000);
	h2EnergyLostInAVsBeta[3*iScintillatorA+2]->SetOption("col");
	h2EnergyLostInAVsBeta[3*iScintillatorA+2]->GetXaxis()->SetTitle("Beta");
	h2EnergyLostInAVsBeta[3*iScintillatorA+2]->GetYaxis()->SetTitle("Energy");
  //       if(fExpertHistos)
	  AddHistogram(h2EnergyLostInAVsBeta[3*iScintillatorA+2]);
      }
    }

  if(fExpertHistos) {
    h2AdcUpVsAdcDown = new TH2F*[12+24];
    for(int iChannel=0; iChannel<12; ++iChannel)
      {
	histname = fName + "_AupVsAdown_A"; histname+=iChannel;
	histtitle = fName + " Adc up vs Adc down for A"; histtitle+=iChannel;
	h2AdcUpVsAdcDown[iChannel]=new TH2F(histname,histtitle,100,0,4000,100,0,4000);
	h2AdcUpVsAdcDown[iChannel]->SetOption("col");
	h2AdcUpVsAdcDown[iChannel]->GetXaxis()->SetTitle("Adc up");
	h2AdcUpVsAdcDown[iChannel]->GetYaxis()->SetTitle("Adc down");
	AddHistogram(h2AdcUpVsAdcDown[iChannel]);
      }

    for(int iChannel=12; iChannel<36; ++iChannel)
      {
	histname = fName + "_BupVsBdown_B"; histname+=iChannel-12;
	histtitle = fName + " Adc up vs Adc down for B"; histtitle+=iChannel-12;
	h2AdcUpVsAdcDown[iChannel]=new TH2F(histname,histtitle,100,0,4000,100,0,4000);
	h2AdcUpVsAdcDown[iChannel]->SetOption("col");
	h2AdcUpVsAdcDown[iChannel]->GetXaxis()->SetTitle("Adc up");
	h2AdcUpVsAdcDown[iChannel]->GetYaxis()->SetTitle("Adc down");
	AddHistogram(h2AdcUpVsAdcDown[iChannel]);
      }

    h2AdcVsPosition = new TH2F*[72];
    for(int iChannel=0; iChannel<12; ++iChannel)
      {
	histname = fName + "_AupVsZ_A"; histname+=iChannel;
	histtitle = fName + " Position Vs Adc Up for A"; histtitle+=iChannel;
	h2AdcVsPosition[iChannel]=new TH2F(histname,histtitle,100,-50,50,500,0,4000);
	h2AdcVsPosition[iChannel]->GetXaxis()->SetTitle("Position");
	h2AdcVsPosition[iChannel]->GetYaxis()->SetTitle("Adc up");
	h2AdcVsPosition[iChannel]->SetOption("col");
	AddHistogram(h2AdcVsPosition[iChannel]);
      }
    for(int iChannel=0; iChannel<12; ++iChannel)
      {
	histname = fName + "_AdownVsZ_A"; histname+=iChannel;
	histtitle = fName + " Position Vs Adc Down for A"; histtitle+=iChannel;
	h2AdcVsPosition[12+iChannel]=new TH2F(histname,histtitle,100,-50, 50,500,0,4000);
	h2AdcVsPosition[12+iChannel]->GetXaxis()->SetTitle("Position");
	h2AdcVsPosition[12+iChannel]->GetYaxis()->SetTitle("Adc down");
	h2AdcVsPosition[12+iChannel]->SetOption("col");
	AddHistogram(h2AdcVsPosition[12+iChannel]);
      }
    for(int iChannel=0; iChannel<24; ++iChannel)
      {
	histname = fName + "_BupVsZ_B"; histname+=iChannel;
	histtitle = fName + " Position Vs Adc Up for B"; histtitle+=iChannel;
	h2AdcVsPosition[24+iChannel]=new TH2F(histname,histtitle,100,-70,70,500,0,4000);
	h2AdcVsPosition[24+iChannel]->GetXaxis()->SetTitle("Position");
	h2AdcVsPosition[24+iChannel]->GetYaxis()->SetTitle("Adc up");
	h2AdcVsPosition[24+iChannel]->SetOption("col");
	AddHistogram(h2AdcVsPosition[24+iChannel]);
      }
    for(int iChannel=0; iChannel<24; ++iChannel)
      {
	histname = fName + "_BdownVsZ_B"; histname+=iChannel;
	histtitle = fName + " Position Vs Adc Down for B"; histtitle+=iChannel;
	h2AdcVsPosition[48+iChannel]=new TH2F(histname,histtitle,100,-70,70,500,0,4000);
	h2AdcVsPosition[48+iChannel]->GetXaxis()->SetTitle("Position");
	h2AdcVsPosition[48+iChannel]->GetYaxis()->SetTitle("Adc down");
	h2AdcVsPosition[48+iChannel]->SetOption("col");
	AddHistogram(h2AdcVsPosition[48+iChannel]);
      }
  }

  histname = fName + "_Multiplicity_trigger";
  histtitle = fName + "_Number of possible triggers";
  h1MultiplicityOfTriggerAB = new TH1I(histname,histtitle,37,0,37);
  AddHistogram(h1MultiplicityOfTriggerAB);

  histname = fName + "_Efficiency_cuts";
  histtitle = fName + "_list of cuts applied to events";
  h1EfficiencyAfterCuts = new TH1I(histname,histtitle,7,0,7);
  AddHistogram(h1EfficiencyAfterCuts);
  
  histname = fName + "_distribution_sector_AiBj";
  histtitle = fName + "_trigger sector by sector";
  h1IdentifiedSectorForTriggerABWithCuts = new TH1I(histname,histtitle,36,0,36);
  AddHistogram(h1IdentifiedSectorForTriggerABWithCuts);

  h2EnergyAVsEnergyB = new TH2F*[N_CHANNEL_SADC_A/2*3];
  for(int iChannel=0; iChannel<N_CHANNEL_SADC_A/2; ++iChannel)
    {
      histname = fName + "_Energy_AVsB_A";
      histname += iChannel;
      histname += "B";
      histname += Modulo((2*iChannel-1),N_CHANNEL_SADC_B/2);
      histtitle = fName + " Energy in A ";
      histtitle += iChannel;
      histtitle += " Vs energy in B ";
      histtitle += Modulo((2*iChannel-1),N_CHANNEL_SADC_B/2);
      h2EnergyAVsEnergyB[3*iChannel]=new TH2F(histname,histtitle,200,0,4000,200,0,4000);
      h2EnergyAVsEnergyB[3*iChannel]->SetOption("cont0");
      h2EnergyAVsEnergyB[3*iChannel]->GetXaxis()->SetTitle("Energy in B");
      h2EnergyAVsEnergyB[3*iChannel]->GetYaxis()->SetTitle("Energy in A");
      AddHistogram(h2EnergyAVsEnergyB[3*iChannel]);

      histname = fName + "_Energy_AVsB_A";
      histname += iChannel;
      histname += "B";
      histname += Modulo((2*iChannel),N_CHANNEL_SADC_B/2);
      histtitle = fName + " Energy in A ";
      histtitle += iChannel;
      histtitle += " Vs energy in B ";
      histtitle += Modulo((2*iChannel),N_CHANNEL_SADC_B/2);
      h2EnergyAVsEnergyB[3*iChannel+1]=new TH2F(histname,histtitle,200,0,4000,200,0,4000);
      h2EnergyAVsEnergyB[3*iChannel+1]->SetOption("cont0");
      h2EnergyAVsEnergyB[3*iChannel+1]->GetXaxis()->SetTitle("Energy in B");
      h2EnergyAVsEnergyB[3*iChannel+1]->GetYaxis()->SetTitle("Energy in A");
      AddHistogram(h2EnergyAVsEnergyB[3*iChannel+1]);

      histname = fName + "_Energy_AVsB_A";
      histname += iChannel;
      histname += "B";
      histname += Modulo((2*iChannel+1),N_CHANNEL_SADC_B/2);
      histtitle = fName + " Energy in A ";
      histtitle += iChannel;
      histtitle += " Vs energy in B ";
      histtitle += Modulo((2*iChannel+1),N_CHANNEL_SADC_B/2);
      h2EnergyAVsEnergyB[3*iChannel+2]=new TH2F(histname,histtitle,200,0,4000,200,0,4000);
      h2EnergyAVsEnergyB[3*iChannel+2]->SetOption("cont0");
      h2EnergyAVsEnergyB[3*iChannel+2]->GetXaxis()->SetTitle("Energy in B");
      h2EnergyAVsEnergyB[3*iChannel+2]->GetYaxis()->SetTitle("Energy in A");
      AddHistogram(h2EnergyAVsEnergyB[3*iChannel+2]);
    }
  histname = fName + "_zvertex";
  histtitle = fName + " Position in target extrapolated";
  h1ZVertex=new TH1F(histname,histtitle,200,-80,80);
  h1ZVertex->GetXaxis()->SetTitle("Z (cm)");
  AddHistogram(h1ZVertex);

  h1tdcA0Down = new TH1F(fName+TString("tdca0do"),"tdc spectrum for A0 down with general condition",100,-16650,-16250);
  AddHistogram(h1tdcA0Down);
}

void GroupRPD::EndEvent(const CS::DaqEvent &event)
{
  Group::EndEvent();
  if (!fOkPlanes) return;

  if(thr_flag)
    TThread::Lock();
  
  const PlaneRPD_F1 *pf1_Bused = pf1_Bl; // balise nicole

  //   float timeUnit=0.0589;

  int debug=0;
  if (debug) printf("********************************************************\n");

  // List of cuts part
  // It will fill the histo h1EfficiencyAfterCuts which show how cuts applied when filling histograms
  // reduce the sample.

  // Rq, j'ai pas ete tres inspiré en l'ecrivant.... c un peu long pour pas grand chose

  bool conditionIsFulfill=false;
  if(event.GetTrigger()&1)
    h1EfficiencyAfterCuts->Fill(0); // total number of events
  for(unsigned int iA=0; iA<12; ++iA)
    {
      if(pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][0] && pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][0] &&
	 psadcl->fDataSent[Pm2Channel('A',iA,'U',"SADC")]  && psadcl->fDataSent[Pm2Channel('A',iA,'D',"SADC")])
	{
	  unsigned int iScintillatorBPrevious=Modulo((2*iA-1),24);
	  if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"F1")][0] && 
	      psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"SADC")]
	      ) ||
	     (
	      pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"F1")][0] &&
	      psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"SADC")]
	      )
	     )
	    {conditionIsFulfill=true;break;}
	  unsigned int iScintillatorBFace=Modulo((2*iA),24);
	  if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"F1")][0] && 
	      psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"SADC")]
	      ) ||
	     (
	      pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"F1")][0] &&
	      psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"SADC")]
	      )
	     )
	    {conditionIsFulfill=true;break;}
	  unsigned int iScintillatorBNext=Modulo((2*iA+1),24);
	  if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"F1")][0] && 
	      psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"SADC")]
	      ) ||
	     (
	      pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"F1")][0] &&
	      psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"SADC")]
	      )
	     )
	    {conditionIsFulfill=true;break;}
	}
    }

  if(conditionIsFulfill)
    {
      if(event.GetTrigger()&1)
	h1EfficiencyAfterCuts->Fill(1);
      conditionIsFulfill=false;
      for(unsigned int iA=0; iA<12; ++iA)
	{
	  if(pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][0] && pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][0] &&
	     psadcl->fDataSent[Pm2Channel('A',iA,'U',"SADC")]  && psadcl->fDataSent[Pm2Channel('A',iA,'D',"SADC")])
	    {
	      unsigned int iScintillatorBPrevious=Modulo((2*iA-1),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"F1")][0] && 
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"SADC")]
		  ) &&
		 (
		  pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"SADC")]
		  )
		 )
		{conditionIsFulfill=true;break;}
	      unsigned int iScintillatorBFace=Modulo((2*iA),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"F1")][0] && 
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"SADC")]
		  ) &&
		 (
		  pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"SADC")]
		  )
		 )
		{conditionIsFulfill=true;break;}
	      unsigned int iScintillatorBNext=Modulo((2*iA+1),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"F1")][0] && 
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"SADC")]
		  ) &&
		 (
		  pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"SADC")]
		  )
		 )
		{conditionIsFulfill=true;break;}
	    }
	}
    }
  
  if(conditionIsFulfill)
    {
      if(event.GetTrigger()&1)
	h1EfficiencyAfterCuts->Fill(2);
      conditionIsFulfill=false;
      for(unsigned int iA=0; iA<12; ++iA)
	{
	  if(pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][0] && pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][0] &&
	     psadcl->fDataSent[Pm2Channel('A',iA,'U',"SADC")]  && psadcl->fDataSent[Pm2Channel('A',iA,'D',"SADC")] &&
	     !pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][1] && !pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][1])
	    {
	      unsigned int iScintillatorBPrevious=Modulo((2*iA-1),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"F1")][0] && 
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"SADC")]
		  ) &&
		 (
		  pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"SADC")]
		  )
		 )
		{conditionIsFulfill=true;break;}
	      unsigned int iScintillatorBFace=Modulo((2*iA),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"F1")][0] && 
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"SADC")]
		  ) &&
		 (
		  pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"SADC")]
		  )
		 )
		{conditionIsFulfill=true;break;}
	      unsigned int iScintillatorBNext=Modulo((2*iA+1),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"F1")][0] && 
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"SADC")]
		  ) &&
		 (
		  pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"SADC")]
		  )
		 )
		{conditionIsFulfill=true;break;}
	    }
	}      
    }
  
  if(conditionIsFulfill)
    {
      if(event.GetTrigger()&1)
	h1EfficiencyAfterCuts->Fill(3);
      conditionIsFulfill=false;
      for(unsigned int iA=0; iA<12; ++iA)
	{
	  if(pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][0] && pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][0] &&
	     psadcl->fDataSent[Pm2Channel('A',iA,'U',"SADC")]  && psadcl->fDataSent[Pm2Channel('A',iA,'D',"SADC")] &&
	     !pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][1] && !pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][1])
	    {
	      unsigned int iScintillatorBPrevious=Modulo((2*iA-1),24);
	      if(pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"F1")][0] && pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"F1")][0] &&
		 psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"SADC")] && psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"SADC")] &&
		 !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"F1")][1] && !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"F1")][1]
		 )
		{conditionIsFulfill=true; break;}
	      unsigned int iScintillatorBFace=Modulo((2*iA),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"F1")][0] && pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"SADC")] && psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"SADC")] &&
		  !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"F1")][1] && !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"F1")][1]
		  )
		 )
		{conditionIsFulfill=true; break;}
	      unsigned int iScintillatorBNext=Modulo((2*iA+1),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"F1")][0] && pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"SADC")] && psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"SADC")] &&
		  !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"F1")][1] && !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"F1")][1]
		  )
		 )
		{conditionIsFulfill=true; break;}
	    }
	}
    }
  
  if(conditionIsFulfill)
    {
      if(event.GetTrigger()&1)
	h1EfficiencyAfterCuts->Fill(4);
      conditionIsFulfill=false;
      for(unsigned int iA=0; iA<12; ++iA)
	{
	  if(pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][0] && pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][0] &&
	     psadcl->fDataSent[Pm2Channel('A',iA,'U',"SADC")]  && psadcl->fDataSent[Pm2Channel('A',iA,'D',"SADC")] &&
	     !pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][1] && !pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][1] &&
	     fabs(GetPosition('A',iA, pf1_A, psadcl))<28
	     )
	    {
	      unsigned int iScintillatorBPrevious=Modulo((2*iA-1),24);
	      if(pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"F1")][0] && pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"F1")][0] &&
		 psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"SADC")] && psadcl->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"SADC")] &&
		 !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"F1")][1] && !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"F1")][1] &&
		 fabs(GetPosition('B',iScintillatorBPrevious, pf1_Bused, psadcl))<59
		 )
		{conditionIsFulfill=true; break;}
	      unsigned int iScintillatorBFace=Modulo((2*iA),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"F1")][0] && pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"SADC")] && psadcl->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"SADC")] &&
		  !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"F1")][1] && !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"F1")][1] &&
		  fabs(GetPosition('B',iScintillatorBFace, pf1_Bused, psadcl))<59
		  )
		 )
		{conditionIsFulfill=true; break;}
	      unsigned int iScintillatorBNext=Modulo((2*iA+1),24);
	      if((pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"F1")][0] && pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"F1")][0] &&
		  psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"SADC")] && psadcl->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"SADC")] &&
		  !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"F1")][1] && !pf1_Bused->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"F1")][1] &&
		  fabs(GetPosition('B',iScintillatorBNext, pf1_Bused, psadcl))<59
		  )
		 )
		{conditionIsFulfill=true; break;}
	    }
	}
    }

  if(conditionIsFulfill)
    {
      if(event.GetTrigger()&1)
	h1EfficiencyAfterCuts->Fill(5);
      conditionIsFulfill=false;
      unsigned int numberOfPossibleTriggerTmp=0;
      for(unsigned int iA=0; iA<12; ++iA)
	{
	  if(CanAiTrig(iA, pf1_A, pf1_Bused, psadcl))
	    ++numberOfPossibleTriggerTmp;
	}
      if(numberOfPossibleTriggerTmp==1)
	{
	  conditionIsFulfill=true;
	}
    }
  
  if(conditionIsFulfill)
    {
      if(event.GetTrigger()&1)
	h1EfficiencyAfterCuts->Fill(6);
      for(unsigned int iA=0; iA<12; ++iA)
	{
	  if(CanAiTrig(iA, pf1_A, pf1_Bused, psadcl, Modulo(2*iA-1,24))){if(event.GetTrigger()&1)h1IdentifiedSectorForTriggerABWithCuts->Fill(3*iA);}
	  if(CanAiTrig(iA, pf1_A, pf1_Bused, psadcl, Modulo(2*iA,24))) {if(event.GetTrigger()&1)h1IdentifiedSectorForTriggerABWithCuts->Fill(3*iA+1);}
	  if(CanAiTrig(iA, pf1_A, pf1_Bused, psadcl, Modulo(2*iA+1,24))) {if(event.GetTrigger()&1)h1IdentifiedSectorForTriggerABWithCuts->Fill(3*iA+2);}
	}
    }

  // Trigger part. (gjegou 03/08/08)
  // This online trigger do not represent the real trigger since it is based on tdclow which have lower threshold than trigger.
  // It should overestimate the multiplicity of triggers
  unsigned int numberOfPossibleTrigger=0;
  unsigned int numberOfPossibleTriggerAB=0;
  bool *AiWhichCanTrig = new bool[12]; // the element i is true if Ai contributed to the trigger (not yet tested)
  for(unsigned int iA=0; iA<12; ++iA)
    {
      if(CanAiTrig(iA, pf1_A, pf1_Bused, psadcl))
	{
	  AiWhichCanTrig[iA]=true;
	  ++numberOfPossibleTrigger;
	}
      else
	AiWhichCanTrig[iA]=false;

      if(CanAiTrig(iA, pf1_A, pf1_Bused, psadcl, Modulo(2*iA-1,24))) {++numberOfPossibleTriggerAB;}
      if(CanAiTrig(iA, pf1_A, pf1_Bused, psadcl, Modulo(2*iA,24))) {++numberOfPossibleTriggerAB;}
      if(CanAiTrig(iA, pf1_A, pf1_Bused, psadcl, Modulo(2*iA+1,24))) {++numberOfPossibleTriggerAB;}
    }
  if(event.GetTrigger()&1)
    h1MultiplicityOfTriggerAB->Fill(numberOfPossibleTriggerAB);  
  // End of trigger part

  for(unsigned int iScintillatorA=0 ; iScintillatorA<12; ++iScintillatorA)
    {
      int iAuF1=Pm2Channel('A',iScintillatorA,'U',"F1");
      int iAdF1=Pm2Channel('A',iScintillatorA,'D',"F1");
      int iAuSadc=Pm2Channel('A',iScintillatorA,'U',"SADC");
      int iAdSadc=Pm2Channel('A',iScintillatorA,'D',"SADC");
      int multUp=0;
      int multDown=0;
      for (int i=0; i<5; i++)
	{
	  if (pf1_A->fDataSent[iAuF1][i])++multUp;
	  if (pf1_A->fDataSent[iAdF1][i])++multDown;
	}
      
      if (debug && multUp>0 && multDown>0)
	printf(" -----------------------------------------------\n In A %i  - mult Up: %i Down: %i\n",iScintillatorA,multUp,multDown);

      //       if(debug && pf1_A->fDataSent[iAuF1][0] && pf1_A->fDataSent[iAdF1][0])	  
      //       printf("TDC : hit in A : %i z=%4.0f \n",iScintillatorA,GetPosition('A',iScintillatorA, pf1_A, psadcl));
      // 
      //       if(debug && pf1_A->fDataSent[iAuF1][1] && pf1_A->fDataSent[iAdF1][1])
      //       printf("TDC : hit in A : %i z=%4.0f \n",iScintillatorA,GetPosition('A',iScintillatorA, pf1_A, psadcl));
      // 
      //       if(debug && pf1_A->fDataSent[iAuF1][2] && pf1_A->fDataSent[iAdF1][2])	  
      //       printf("TDC : hit in A : %i z=%4.0f \n",iScintillatorA,GetPosition('A',iScintillatorA, pf1_A, psadcl));
     
      if(debug && psadcl->fDataSent[iAuSadc] && psadcl->fDataSent[iAdSadc])	  
	printf("-----------------------------------------------\n ADC : hit in A : %i\n",iScintillatorA);
      
      if(1
   	 && numberOfPossibleTriggerAB==1 && AiWhichCanTrig[iScintillatorA]   // balise for nicole; numberOfPossibleTrigger>=1
	 )
	{
	  if(1
	     && pf1_A->fDataSent[iAuF1][0] && pf1_A->fDataSent[iAdF1][0]
	     && psadcl->fDataSent[iAuSadc] && psadcl->fDataSent[iAdSadc]
	     && multUp==1 && multDown==1
	     )
	    {
	      h1Za[iScintillatorA]->Fill(GetPosition('A', iScintillatorA, pf1_A, psadcl));
	      if(fabs(GetPosition('A',iScintillatorA, pf1_A, psadcl))<28)
		{
		  if(AiWhichCanTrig[0])
		    h1tdcA0Down->Fill(pf1_A->fTdcRpd[Pm2Channel('A',0,'D',"F1")][0]);
		  if(fExpertHistos) {
		    h2AdcVsPosition[iScintillatorA]->Fill(GetPosition('A', iScintillatorA, pf1_A, psadcl), psadcl->fAdc[iAuSadc]);
		    h2AdcVsPosition[12+iScintillatorA]->Fill(GetPosition('A', iScintillatorA, pf1_A, psadcl), psadcl->fAdc[iAdSadc]);
		    h1TupMinusTdownInA[iScintillatorA]->Fill(GetTimeDifference('A', iScintillatorA, pf1_A, psadcl));
		  }
		  float ElossA=sqrt(psadcl->fAdc[iAuSadc]*psadcl->fAdc[iAdSadc]);

		  if (debug) printf("TDC && ADC hit in A : %i z=%4.0f E=%4.0f \n",iScintillatorA,GetPosition('A',iScintillatorA, pf1_A, psadcl),ElossA);
	  
		  // filling histograms for the A(i) - B(2i-1) coincidences
	  
		  int iScintillatorBe=Modulo(2*iScintillatorA-1,24);
		  int iBuF1=Pm2Channel('B',iScintillatorBe,'U',"F1");
		  int iBdF1=Pm2Channel('B',iScintillatorBe,'D',"F1");
		  int iBuSadc=Pm2Channel('B',iScintillatorBe,'U',"SADC");
		  int iBdSadc=Pm2Channel('B',iScintillatorBe,'D',"SADC");
		  int multUpB=0;
		  int multDownB=0;
		  for (int i=0; i<F1_MULTIPLICITY; i++){   // is it really 5 ? -> gjegou : 5 is the size of the array, if the F1 card send 6 hits, the 6th is not taken into account.
		    if (pf1_Bused->fDataSent[iBuF1][i])multUpB++;
		    if (pf1_Bused->fDataSent[iBdF1][i])multDownB++;
		  }

		  if (debug) printf(" In B %i  - mult Up: %i Down: %i\n",iScintillatorBe,multUpB,multDownB);

		  if(1
		     && pf1_Bused->fDataSent[iBuF1][0] && pf1_Bused->fDataSent[iBdF1][0] 
		     && psadcl->fDataSent[iBuSadc] && psadcl->fDataSent[iBdSadc]
		     && (multUpB==1 && multDownB==1)
		     && fabs(GetPosition('B',iScintillatorBe, pf1_Bused, psadcl))<59
		     )
		    {
		      float tof    = GetTof(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl);
		      float beta   = GetBeta(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl);
//		      float polarAngle = TMath::ATan(TMath::Abs(65/(GetPosition('A',iScintillatorA, pf1_A, psadcl)-GetPosition('B',iScintillatorBe, pf1_Bused, psadcl)))); // this variable is the polar angle of the trajectory BETWEEN 0 AND PI/2
		      float polarAngle = atan(TMath::Abs(65/(GetPosition('A',iScintillatorA, pf1_A, psadcl)-GetPosition('B',iScintillatorBe, pf1_Bused, psadcl)))); // this variable is the polar angle of the trajectory BETWEEN 0 AND PI/2
		      float ElossB = sqrt(psadcl->fAdc[iBuSadc]*psadcl->fAdc[iBdSadc])/sin(polarAngle);
		      ElossA = ElossA/sin(polarAngle);
		      h1Tof[3*iScintillatorA]->Fill(tof);
		      h1ZVertex->Fill(GetZPositionInTarget(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl));
		      h2EnergyAVsEnergyB[3*iScintillatorA]->Fill(sqrt(psadcl->fAdc[iBuSadc]*psadcl->fAdc[iBdSadc]),sqrt(psadcl->fAdc[iAuSadc]*psadcl->fAdc[iAdSadc]));
		      if (debug)  printf("   hit in B : %i  z=%4.0f E=%4.0f \n",iScintillatorBe,GetZPositionInTarget(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl),ElossB);
		      if(1
			 // 	     && ElossA<(beta-.75)*(beta-.75)*5000+500
			 )h2EnergyLostInBVsBeta[3*iScintillatorA]->Fill(beta,ElossB);
		      if(fExpertHistos) if(1        
			 // 	     && ElossA<(beta-.75)*(beta-.75)*5000+500	     
			 //       && ElossB >500
			 )h2EnergyLostInAVsBeta[3*iScintillatorA]->Fill(beta,ElossA);	    
		    }      
	  
		  // filling histograms for the A(i) - B(2i) coincidences
	  
		  iScintillatorBe=2*iScintillatorA;
		  iBuF1=Pm2Channel('B',iScintillatorBe,'U',"F1");
		  iBdF1=Pm2Channel('B',iScintillatorBe,'D',"F1");
		  iBuSadc=Pm2Channel('B',iScintillatorBe,'U',"SADC");
		  iBdSadc=Pm2Channel('B',iScintillatorBe,'D',"SADC");
		  multUpB=0;
		  multDownB=0;
		  for (int i=0; i<5; i++){
		    if (pf1_Bused->fDataSent[iBuF1][i])multUpB++;
		    if (pf1_Bused->fDataSent[iBdF1][i])multDownB++;
		  }
		  if (debug) printf(" In B %i  - mult Up: %i Down: %i\n",iScintillatorBe,multUpB,multDownB);
		  if( 1
		      && pf1_Bused->fDataSent[iBuF1][0] && pf1_Bused->fDataSent[iBdF1][0] 
		      && psadcl->fDataSent[iBuSadc] && psadcl->fDataSent[iBdSadc]
		      && (multUpB==1 && multDownB==1)
		      && fabs(GetPosition('B',iScintillatorBe, pf1_Bused, psadcl))<59
		      )
		    {
		      float tof= GetTof(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl);
		      float beta = GetBeta(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl);
// 		      float polarAngle = TMath::ATan(TMath::Abs(65/(GetPosition('A',iScintillatorA, pf1_A, psadcl)-GetPosition('B',iScintillatorBe, pf1_Bused, psadcl)))); // this variable is the polar angle of the trajectory BETWEEN 0 AND PI/2
		      float polarAngle = atan(TMath::Abs(65/(GetPosition('A',iScintillatorA, pf1_A, psadcl)-GetPosition('B',iScintillatorBe, pf1_Bused, psadcl)))); // this variable is the polar angle of the trajectory BETWEEN 0 AND PI/2
		      float ElossB=sqrt(psadcl->fAdc[iBuSadc]*psadcl->fAdc[iBdSadc])/sin(polarAngle);
		      ElossA = ElossA/sin(polarAngle);
		      h1Tof[3*iScintillatorA+1]->Fill(tof);
		      h1ZVertex->Fill(GetZPositionInTarget(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl));
		      h2EnergyAVsEnergyB[3*iScintillatorA+1]->Fill(sqrt(psadcl->fAdc[iBuSadc]*psadcl->fAdc[iBdSadc]),sqrt(psadcl->fAdc[iAuSadc]*psadcl->fAdc[iAdSadc]));
		      if (debug) printf("   hit in B : %i  z=%4.0f E=%4.0f \n",iScintillatorBe,GetZPositionInTarget(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl),ElossB);
		      if(1
			 //&& psadcl->fAmplitudeMaximum[iBuSadc]>40 && 
			 //&& psadcl->fAmplitudeMaximum[iBdSadc]>40 
			 //&& ElossA<(beta-.75)*(beta-.75)*5000+500
			 )h2EnergyLostInBVsBeta[3*iScintillatorA+1]->Fill(beta,ElossB);	    
		      if(fExpertHistos) if(1
			 //&& ElossB >500
			 //&& psadcl->fAdc[iAuSadc]>200 && psadcl->fAdc[iAdSadc]>200
			 //&& ElossA<(beta-.75)*(beta-.75)*5000+500	     
			 )h2EnergyLostInAVsBeta[3*iScintillatorA+1]->Fill(beta,ElossA);	    
		    }      
	  
		  // filling histograms for the A(i) - B(2i+1) coincidences
	  
		  iScintillatorBe=Modulo((2*iScintillatorA+1),24);
		  iBuF1=Pm2Channel('B',iScintillatorBe,'U',"F1");
		  iBdF1=Pm2Channel('B',iScintillatorBe,'D',"F1");
		  iBuSadc=Pm2Channel('B',iScintillatorBe,'U',"SADC");
		  iBdSadc=Pm2Channel('B',iScintillatorBe,'D',"SADC");
		  multUpB=0;
		  multDownB=0;
		  for (int i=0; i<5; i++){
		    if (pf1_Bused->fDataSent[iBuF1][i])multUpB++;
		    if (pf1_Bused->fDataSent[iBdF1][i])multDownB++;
		  }
		  if (debug) printf(" In B %i  - mult Up: %i Down: %i\n",iScintillatorBe,multUpB,multDownB);
		  if( 1
		      && pf1_Bused->fDataSent[iBuF1][0] && pf1_Bused->fDataSent[iBdF1][0] 
		      && psadcl->fDataSent[iBuSadc] && psadcl->fDataSent[iBdSadc]
		      && (multUpB==1 && multDownB==1)
		      && fabs(GetPosition('B',iScintillatorBe, pf1_Bused, psadcl))<59
		      )
		    {
		      float tof= GetTof(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl);
		      float beta = GetBeta(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl);
// 		      float polarAngle = TMath::ATan(TMath::Abs(65/(GetPosition('A',iScintillatorA, pf1_A, psadcl)-GetPosition('B',iScintillatorBe, pf1_Bused, psadcl)))); // this variable is the polar angle of the trajectory BETWEEN 0 AND PI/2
		      float polarAngle = atan(TMath::Abs(65/(GetPosition('A',iScintillatorA, pf1_A, psadcl)-GetPosition('B',iScintillatorBe, pf1_Bused, psadcl)))); // this variable is the polar angle of the trajectory BETWEEN 0 AND PI/2
		      float ElossB=sqrt(psadcl->fAdc[iBuSadc]*psadcl->fAdc[iBdSadc])/sin(polarAngle);
		      ElossA = ElossA/sin(polarAngle);
		      h1Tof[3*iScintillatorA+2]->Fill(tof);
		      h1ZVertex->Fill(GetZPositionInTarget(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl));
		      h2EnergyAVsEnergyB[3*iScintillatorA+2]->Fill(sqrt(psadcl->fAdc[iBuSadc]*psadcl->fAdc[iBdSadc]),sqrt(psadcl->fAdc[iAuSadc]*psadcl->fAdc[iAdSadc]));
		      if (debug) printf("   hit in B : %i  z=%4.0f E=%4.0f \n",iScintillatorBe,GetZPositionInTarget(iScintillatorA, iScintillatorBe, pf1_A, pf1_Bused, psadcl),ElossB);
		      if(1
			 //&& ElossA<(beta-.75)*(beta-.75)*5000+500
			 )
			h2EnergyLostInBVsBeta[3*iScintillatorA+2]->Fill(beta,ElossB);
		      if(fExpertHistos) if(1
			 //&& ElossB >500
			 //&& ElossA<(beta-.75)*(beta-.75)*5000+500	     
			 )
			h2EnergyLostInAVsBeta[3*iScintillatorA+2]->Fill(beta,ElossA);	    
		    }      
      
		  //if(fabs(GetPosition('A',iScintillatorA, pf1_A, psadcl))<5.)
		  // 	      if( 1
		  // 		  && psadcl->fDataSent[iAuSadc] && psadcl->fDataSent[iAdSadc]
		  // 		  )
		  if(fExpertHistos) {
		    h2AdcUpVsAdcDown[iScintillatorA]->Fill(psadcl->fAdc[iAuSadc],psadcl->fAdc[iAdSadc]);
	          }
		} 
	    }
	}
    }

  for(unsigned int iScintillatorB=0 ; iScintillatorB<24; ++iScintillatorB)
    {
      int iBuF1=Pm2Channel('B',iScintillatorB,'U',"F1");
      int iBdF1=Pm2Channel('B',iScintillatorB,'D',"F1");
      if(pf1_Bused->fDataSent[iBuF1][0] && pf1_Bused->fDataSent[iBdF1][0] && numberOfPossibleTrigger==1)
	{
	  h1Zb[iScintillatorB]->Fill(GetPosition('B',iScintillatorB,pf1_Bused, psadcl));
	  if(fExpertHistos) h1TupMinusTdownInB[iScintillatorB]->Fill(GetTimeDifference('B',iScintillatorB,pf1_Bused,psadcl));
	  
	  int iBuSadc=Pm2Channel('B',iScintillatorB,'U',"SADC");
	  int iBdSadc=Pm2Channel('B',iScintillatorB,'D',"SADC");
	  if(fExpertHistos) {
	    if(psadcl->fDataSent[iBuSadc] && psadcl->fDataSent[iBdSadc])
	      {
	          h2AdcVsPosition[24+iScintillatorB]->Fill(GetPosition('B', iScintillatorB, pf1_Bused, psadcl), psadcl->fAdc[iBuSadc]);
	          h2AdcVsPosition[24+24+iScintillatorB]->Fill(GetPosition('B', iScintillatorB, pf1_Bused, psadcl), psadcl->fAdc[iBdSadc]);
	          // 		  if(psadcl->fDataSent[iBuSadc] && psadcl->fDataSent[iBdSadc] /*&& fabs(GetPosition('B',iScintillatorB, pf1_Bused, psadcl))<5.*/)
	          h2AdcUpVsAdcDown[12+iScintillatorB]->Fill(psadcl->fAdc[iBuSadc],psadcl->fAdc[iBdSadc]);
	      }
	  }
	}
    }
  
  if(thr_flag)
    TThread::UnLock();
}

bool GroupRPD::CanAiTrig(unsigned int iA, const PlaneRPD_F1 *pf1_A, const PlaneRPD_F1 *pf1_B, const PlaneRPD_SADC *psadc, const unsigned int iB)
{
  unsigned int iScintillatorBPrevious=Modulo((2*iA-1),24);
  unsigned int iScintillatorBFace=Modulo((2*iA),24);
  unsigned int iScintillatorBNext=Modulo((2*iA+1),24);
  if(
     pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][0] && pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][0] &&
     psadc->fDataSent[Pm2Channel('A',iA,'U',"SADC")]  && psadc->fDataSent[Pm2Channel('A',iA,'D',"SADC")] &&
     !pf1_A->fDataSent[Pm2Channel('A',iA,'U',"F1")][1] && !pf1_A->fDataSent[Pm2Channel('A',iA,'D',"F1")][1]// &&
//      fabs(GetPosition('A',iA, pf1_A, psadc))<28
     )
    {
      if((iScintillatorBPrevious==iB || iB==100)&&
	 (
	  pf1_B->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"F1")][0] && 
	  psadc->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"SADC")] &&
	  !pf1_B->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'U',"F1")][1] 
	  ) &&
	 (
	  pf1_B->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"F1")][0] &&
	  psadc->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"SADC")] &&
	  !pf1_B->fDataSent[Pm2Channel('B',iScintillatorBPrevious,'D',"F1")][1]
	  )// &&
// 	 fabs(GetPosition('B',iScintillatorBPrevious, pf1_B, psadc))<59
	 )
	{
	  return true;
	}
      else if((iScintillatorBFace==iB || iB==100)&&
	      (
	       pf1_B->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"F1")][0] &&
	       psadc->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"SADC")] &&
	       !pf1_B->fDataSent[Pm2Channel('B',iScintillatorBFace,'U',"F1")][1]
	       ) &&
	      (
	       pf1_B->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"F1")][0] &&
	       psadc->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"SADC")] &&
	       !pf1_B->fDataSent[Pm2Channel('B',iScintillatorBFace,'D',"F1")][1]
	       )// &&
// 	      fabs(GetPosition('B',iScintillatorBFace, pf1_B, psadc))<59
	      )
	{
	  return true;
	}
      else if((iScintillatorBNext==iB || iB==100)&&
	      (
	       pf1_B->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"F1")][0] && 
	       psadc->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"SADC")] &&
	       !pf1_B->fDataSent[Pm2Channel('B',iScintillatorBNext,'U',"F1")][1] 
	       ) &&
	      (
	       pf1_B->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"F1")][0] && 
	       psadc->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"SADC")] &&
	       !pf1_B->fDataSent[Pm2Channel('B',iScintillatorBNext,'D',"F1")][1] 
	       )// &&
// 	      fabs(GetPosition('B',iScintillatorBNext, pf1_B, psadc))<59
	      )
	{
	  return true;
	}
      else
	return false;
    }
  else 
    return false;
}

int GroupRPD::Pm2Channel(char scint, unsigned int n, char side, TString type)
{
// mapping of elements in the TDC or ADC buffers
// tells the channel number in TDC or ADC array corresponding to the element that is wanted
// for the moment only valid for "low" thresholds
  if(type=="F1")
    {
      if(scint=='A')
	{
	  if(side=='U') return n;
	  else if(side=='D') return 32+n;
	}
      else if(scint=='B')
	{
	  if(n<12)
	    {
	      if(side=='U') return n;
	      else if(side=='D') return 32+n ;
	    }
	  else
	    {
	      if(side=='U') return 16+n-12;
	      else if(side=='D') return 48+n-12 ;	      
	    }
	}
      else
	{
	  cerr<<"ERROR in int GroupRPD::Pm2Channel(char scint, unsigned int n, char side, TString type), scint!='A' and scint!='B'";
	  return -1;
	}
    }
  else if(type=="SADC")
    {
      if(scint=='A')
	{
	  if(side=='U') return n;
	  else if(side=='D') return 12+n;
	}
      else if(scint=='B')
	{
	  if(side=='U') return 24+n;
	  else if(side=='D') return 24+24+n;
	}
      else
	{
	  cerr<<"ERROR in int GroupRPD::Pm2Channel(char scint, unsigned int n, char side, TString type), scint!='A' and scint!='B'";
	  return -1;
	}
    }
  else
    {
      cerr<<"ERROR in int GroupRPD::Pm2Channel(char scint, unsigned int n, char side, TString type), type!=\"F1\" and type!=\"SADC\"";
      return -1;
    }
  return -1;
}

float GroupRPD::GetTimeDifferenceOffset(char scintillator, unsigned int n)
{

  // In this method, index of array correspond to hardware index of the RPD => offset for Bi = timeDifferenceOffsetForB[i]
  if(scintillator=='A')
    {
      float timeDifferenceOffsetForA[12];
      timeDifferenceOffsetForA[0]=-1.4;
      timeDifferenceOffsetForA[1]=-0.4;
      timeDifferenceOffsetForA[2]=4.6;
      timeDifferenceOffsetForA[3]=3.7;
      timeDifferenceOffsetForA[4]=2.7;
      timeDifferenceOffsetForA[5]=5.4;
      timeDifferenceOffsetForA[6]=2.1;
      timeDifferenceOffsetForA[7]=1.53;
      timeDifferenceOffsetForA[8]=1.7;
      timeDifferenceOffsetForA[9]=-1.7;
      timeDifferenceOffsetForA[10]=0.6;
      timeDifferenceOffsetForA[11]=-0.18; 
      return timeDifferenceOffsetForA[n];
    }
  else if(scintillator=='B')
    {
      float timeDifferenceOffsetForB[24];
      timeDifferenceOffsetForB[0]=1.2;
      timeDifferenceOffsetForB[1]=-1.;
      timeDifferenceOffsetForB[2]=0.5;
      timeDifferenceOffsetForB[3]=2.;
      timeDifferenceOffsetForB[4]=3.4;
      timeDifferenceOffsetForB[5]=2.4;
      timeDifferenceOffsetForB[6]=-0.8;
      timeDifferenceOffsetForB[7]=0.2;
      timeDifferenceOffsetForB[8]=3.5;
      timeDifferenceOffsetForB[9]=-1.;
      timeDifferenceOffsetForB[10]=0.4;
      timeDifferenceOffsetForB[11]=-2.1;
      timeDifferenceOffsetForB[12]=-0.7;//-5.7;  // some problems to select a peak
      timeDifferenceOffsetForB[13]=3.8;
      timeDifferenceOffsetForB[14]=-3;
      timeDifferenceOffsetForB[15]=0.3;
      timeDifferenceOffsetForB[16]=-4.4;
      timeDifferenceOffsetForB[17]=0.5;
      timeDifferenceOffsetForB[18]=-0.3;
      timeDifferenceOffsetForB[19]=-1.2;
      timeDifferenceOffsetForB[20]=-2.1;
      timeDifferenceOffsetForB[21]=-0.9;
      timeDifferenceOffsetForB[22]=-3.2;
      timeDifferenceOffsetForB[23]=-5.;//-1.8; changed 05/07/08 by gjegou
      return timeDifferenceOffsetForB[n];
    }
  else
    {
      cerr<<"ERROR in float GroupRPD::GetTimeDifferenceOffset(char scintillator, unsigned int n) : scintillator!='A' and scintillator!='B'\n";
      cerr<<"0 returned\n";
      return 0;
    }
}

float GroupRPD::GetLightSpeed(char scintillator, unsigned int n)
{
  if(scintillator=='A')
    {
      float lightSpeedInA[12];
      lightSpeedInA[0]=13;
      lightSpeedInA[1]=13;
      lightSpeedInA[2]=13;
      lightSpeedInA[3]=13;
      lightSpeedInA[4]=13;
      lightSpeedInA[5]=13;
      lightSpeedInA[6]=13;
      lightSpeedInA[7]=13;
      lightSpeedInA[8]=13;
      lightSpeedInA[9]=13;
      lightSpeedInA[10]=13;
      lightSpeedInA[11]=13;
      return lightSpeedInA[n];
    }
  else if(scintillator=='B')
    {
      float lightSpeedInB[24];
      lightSpeedInB[0]=13;
      lightSpeedInB[1]=13;
      lightSpeedInB[2]=13;
      lightSpeedInB[3]=13;
      lightSpeedInB[4]=13;
      lightSpeedInB[5]=13;
      lightSpeedInB[6]=13;
      lightSpeedInB[7]=13;
      lightSpeedInB[8]=13;
      lightSpeedInB[9]=13;
      lightSpeedInB[10]=13;
      lightSpeedInB[11]=13;
      lightSpeedInB[12]=13;
      lightSpeedInB[13]=13;
      lightSpeedInB[14]=13;
      lightSpeedInB[15]=13;
      lightSpeedInB[16]=13;
      lightSpeedInB[17]=13;
      lightSpeedInB[18]=13;
      lightSpeedInB[19]=13;
      lightSpeedInB[20]=13;
      lightSpeedInB[21]=13;
      lightSpeedInB[22]=13;
      lightSpeedInB[23]=13;
      return lightSpeedInB[n];
    }
  else
    {
      cerr<<"ERROR in float GroupRPD::GetLightSpeed(char scintillator, unsigned int n) : scintillator!='A' and scintillator!='B'\n";
      cerr<<"0 returned\n";
      return 0;
    }
    
}

float GroupRPD::GetTimeDifference(char scintillator, int iScintillator, const PlaneRPD_F1 *pf1, const PlaneRPD_SADC *psadc) // input iScintillator is the hardware number, not the channel
{
  int itu=Pm2Channel(scintillator,iScintillator,'U',"F1");
  int itd=Pm2Channel(scintillator,iScintillator,'D',"F1");
  float timeUnit=0.058;
//   for (int iu=0;iu<5;iu++){
//      for (int id=0;id<5;id++){
//         if (pf1->fDataSent[itu][iu]&&pf1->fDataSent[itd][id]) printf("%i %i dt=%f \n",iu,id,(pf1->fTdcRpd[itu][iu]-pf1->fTdcRpd[itd][id])*timeUnit);
//      }
//   }
  
  return (pf1->fTdcRpd[itu][0]-GetWalkCorrection(scintillator, iScintillator, 'U', psadc)-(pf1->fTdcRpd[itd][0]-GetWalkCorrection(scintillator, iScintillator, 'D', psadc)))*timeUnit + GetTimeDifferenceOffset(scintillator,iScintillator);

}

float GroupRPD::GetTimeSum(char scintillator, int iScintillator, const PlaneRPD_F1 *pf1, const PlaneRPD_SADC *psadc) // input iScintillator is the hardware number, not the channel
{
  
  int itu=Pm2Channel(scintillator,iScintillator,'U',"F1");
  int itd=Pm2Channel(scintillator,iScintillator,'D',"F1");
  float timeUnit=0.058;
  
  return (pf1->fTdcRpd[itu][0]-GetWalkCorrection(scintillator, iScintillator, 'U', psadc)+
	  pf1->fTdcRpd[itd][0]-GetWalkCorrection(scintillator, iScintillator, 'D', psadc))*timeUnit;
}

float GroupRPD::GetPosition(char scintillator, int iScintillator, const PlaneRPD_F1 *pf1, const PlaneRPD_SADC *psadc)
{
  float bOffset=0;
  if(scintillator=='B')bOffset=10;  // the position of the laser for B scintillators is 10 cm downstream compare to A scintillators
  return GetTimeDifference(scintillator, iScintillator, pf1, psadc)*GetLightSpeed(scintillator,iScintillator)/2 + bOffset;
}

float GroupRPD::GetTof(int iScintillatorA,int iScintillatorB, const PlaneRPD_F1 *pf1_A, const PlaneRPD_F1 *pf1_B, const PlaneRPD_SADC *psadc)
{
  return  GetTimeSum('B', iScintillatorB, pf1_B, psadc)/2. - GetTimeSum('A', iScintillatorA, pf1_A, psadc)/2. + GetTofOffset(iScintillatorA, iScintillatorB)/2. +5/2.+5;
  // +5 comes from the calibration with laser, we measured calib constants with the laser and the tof comes from the difference of lenght of fibers
  // in time units ->5ns
}

float GroupRPD::GetBeta(int iScintillatorA,int iScintillatorB, const PlaneRPD_F1 *pf1_A, const PlaneRPD_F1 *pf1_B, const PlaneRPD_SADC *psadc)
{
  return sqrt(pow(65.,2)+pow(GetPosition('B',iScintillatorB, pf1_B, psadc)-GetPosition('A',iScintillatorA, pf1_A, psadc),2))
    /(GetTof(iScintillatorA, iScintillatorB, pf1_A, pf1_B, psadc))/30.;  
}

float GroupRPD::GetZPositionInTarget(int iScintillatorA,int iScintillatorB, const PlaneRPD_F1 *pf1_A, const PlaneRPD_F1 *pf1_B, const PlaneRPD_SADC *psadc)
{
  float zA=GetPosition('A',iScintillatorA, pf1_A, psadc);
  return 12.0/65.*(GetPosition('B', iScintillatorB, pf1_B, psadc)-zA)+zA;
}

float GroupRPD::GetTofOffset(int iScintillatorA, int iScintillatorB)
{

  float tofOffset[12][24];

  tofOffset[0][23]=-9.4+7;
  tofOffset[0][0]=-10.9;
  tofOffset[0][1]=-14.8;

  tofOffset[1][1]=-17.7+2;
  tofOffset[1][2]=-15.3;
  tofOffset[1][3]=-15.1;

  tofOffset[2][3]=-10.3;
  tofOffset[2][4]=-11.8;
  tofOffset[2][5]=-9.9;

  tofOffset[3][5]=-22.9;
  tofOffset[3][6]=-23.5;
  tofOffset[3][7]=-22.8;

  tofOffset[4][7]=-17.6;
  tofOffset[4][8]=-17.5;
  tofOffset[4][9]=-10.0;

  tofOffset[5][9]=-10.8;
  tofOffset[5][10]=-13.1;
  tofOffset[5][11]=-9.8;

  tofOffset[6][11]=-9.7;
  tofOffset[6][12]=-23.4;//-28.6;  // diffictult to select the peak
  tofOffset[6][13]=-12.9;

  tofOffset[7][13]=-17.0;
  tofOffset[7][14]=-21.2;
  tofOffset[7][15]=-16.9;
  
  tofOffset[8][15]=-8.7+6.;
  tofOffset[8][16]=-15.5+6.;
  tofOffset[8][17]=-2.9+6.;

  tofOffset[9][17]=-9.5;
  tofOffset[9][18]=-19.3;
  tofOffset[9][19]=-18.1;

  tofOffset[10][19]=-16.1;
  tofOffset[10][20]=-14.6;
  tofOffset[10][21]=-10.4;

  tofOffset[11][21]=-15.2;
  tofOffset[11][22]=-24.1;
  tofOffset[11][23]=-21.3+7;//

  return tofOffset[iScintillatorA][iScintillatorB];

}

float GroupRPD::GetWalkCorrection(char scintillator, int iScintillator, char side, const PlaneRPD_SADC *psadc)
{
  unsigned int iSadc=Pm2Channel(scintillator,iScintillator,side,"SADC");

  float riseTime = 8;
  float f1Treshold = 30;
  float walkCorrection=0;

  if(psadc->fDataSent[iSadc])
    walkCorrection = riseTime * f1Treshold/(2.2*1.5/0.7*psadc->fAmplitudeMaximum[iSadc]);

  return 0;
//   return walkCorrection;

}
