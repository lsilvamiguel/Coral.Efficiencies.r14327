/*
 * PlaneCamera.cc
 *
 * Each ring is mapped to one plane; it is then discriminated between up and down;
 * For ring A there are 2*12 Gandalf channels to process for RPD, 2*24 for Camera
 * For ring B there are also 2*24 channels to process;
 *
 * This class generates plots for an overview over the plane and in addition,
 * if 'expertPlots' are enabled, it instantiates a 'PlaneGandalf' object for each
 * channel, where detailed informations of the channel can be seen.
 *
 *  Created on: Apr 18, 2012
 *      Author: Matthias Gorzellik
 */

#include "PlaneCamera.h"
#include "Reference.h"

ClassImp(PlaneCamera);

PlaneCamera::PlaneCamera(const char *detname) : Plane(detname)
{
	N_CHANNEL=24;
	calibTime = 0;
}
void PlaneCamera::Init(TTree* tree)
{
	//Plane::Init(tree);

	//N_CHANNEL=24;

//#ifdef PROTO
//	if (getType()==PlaneCamera::Aup || getType()==PlaneCamera::Adown){
//		N_CHANNEL=12;
//	}
//#endif


	TString histname = fName + "_NbOfData";
	TString histtitle = fName + " Number of channels which sent data";
	h1NumberOfSentData=new TH1F_Ref(histname,histtitle,N_CHANNEL,0,N_CHANNEL,fRateCounter);
	h1NumberOfSentData->GetXaxis()->SetTitle("Number of sent data");
	AddHistogram(h1NumberOfSentData);

	histname = fName + "_BaselineVsChannel";
	histtitle = fName + " Baseline spectrum vs channels";
	h2BaselineVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 100, 0, 500);
	h2BaselineVsChannel->SetOption("colz");
	h2BaselineVsChannel->GetXaxis()->SetTitle("Channel");
	h2BaselineVsChannel->GetYaxis()->SetTitle("Baseline");
	AddHistogram(h2BaselineVsChannel);

	histname = fName + "_BaselineCompVsChannel";
	histtitle = fName + " Baseline computed spectrum vs channels";
	h2BaselineCompVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 100, 0, 500);
	h2BaselineCompVsChannel->SetOption("colz");
	h2BaselineCompVsChannel->GetXaxis()->SetTitle("Channel");
	h2BaselineCompVsChannel->GetYaxis()->SetTitle("Baseline");
	AddHistogram(h2BaselineCompVsChannel);

	histname = fName + "_FrameTimeVsChannel";
	histtitle = fName + " FrameTime spectrum vs channels";
	h2FrameTimeVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 200, 0, 200);
	h2FrameTimeVsChannel->SetOption("colz");
	h2FrameTimeVsChannel->GetXaxis()->SetTitle("Channel");
	h2FrameTimeVsChannel->GetYaxis()->SetTitle("FrameTime");
	AddHistogram(h2FrameTimeVsChannel);

	histname = fName + "_MaxAmplVsChannel";
	histtitle = fName + " Maximum amplitude spectrum vs channel";
	h2AmplitudeMaximumVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 200, 0, 4096);
	h2AmplitudeMaximumVsChannel->SetOption("colz");
	h2AmplitudeMaximumVsChannel->GetXaxis()->SetTitle("Channel");
	h2AmplitudeMaximumVsChannel->GetYaxis()->SetTitle("Maximum of amplitude");
	AddHistogram(h2AmplitudeMaximumVsChannel);

	histname = fName + "_IntegralVsChannel";
	histtitle = fName + " Integral spectrum vs channel";
	h2IntegralVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 200, 0, 100000);
	h2IntegralVsChannel->SetOption("colz");
	h2IntegralVsChannel->GetXaxis()->SetTitle("Channel");
	h2IntegralVsChannel->GetYaxis()->SetTitle("Integral");
	AddHistogram(h2IntegralVsChannel);

	histname = fName + "_HighResVsChannel";
	histtitle = fName + " HighRes Time spectrum vs channel";
	h2HighResVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 100, 0, 1024);
	h2HighResVsChannel->SetOption("colz");
	h2HighResVsChannel->GetXaxis()->SetTitle("Channel");
	h2HighResVsChannel->GetYaxis()->SetTitle("Integral");
	AddHistogram(h2HighResVsChannel);

	histname = fName + "_TimeDecodedVsChannel";
	histtitle = fName + " TimeDecoded vs channel";
	h2TimeVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 100, -2500,-1500);
	h2TimeVsChannel->SetOption("colz");
	h2TimeVsChannel->GetXaxis()->SetTitle("Channel");
	h2TimeVsChannel->GetYaxis()->SetTitle("TimeDecoded in ns");
	AddHistogram(h2TimeVsChannel);

	histname = fName + "_MultiVsChannel";
	histtitle = fName + " multiplicity spectrum vs channel";
	h2MultiVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 16, 0, 16);
	h2MultiVsChannel->SetOption("colz");
	h2MultiVsChannel->GetXaxis()->SetTitle("Channel");
	h2MultiVsChannel->GetYaxis()->SetTitle("Hits");
	AddHistogram(h2MultiVsChannel);

#ifdef PROTO
	histname = fName + "_MultDiffGanF1VsChannel";
	histtitle = fName + " difference in multiplicity (Ngan - Nf1) vs channel";
	h2MultiDiffGanF1VsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 11, -5.5, 5.5);
	h2MultiDiffGanF1VsChannel->SetOption("text");
	h2MultiDiffGanF1VsChannel->GetXaxis()->SetTitle("Channel");
	h2MultiDiffGanF1VsChannel->GetYaxis()->SetTitle("Hit Difference");
	AddHistogram(h2MultiDiffGanF1VsChannel);

	histname = fName + "_TimeDiffGanF1VsChannel";
	histtitle = fName + " difference in Time (Tgan - Tf1) vs channel";
	h2TimeDiffGanF1VsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 400, 0, 1000);
	h2TimeDiffGanF1VsChannel->SetOption("colz");	
	h2TimeDiffGanF1VsChannel->GetXaxis()->SetTitle("Channel");
	h2TimeDiffGanF1VsChannel->GetYaxis()->SetTitle("Time Difference in ns");
	AddHistogram(h2TimeDiffGanF1VsChannel);

	histname = fName + "_SendDataGanF1VsChannel";
	histtitle = fName + " Send Data from Gan/F1  vs channel";
	h2DataTransmission=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 4, 0, 4);
	h2DataTransmission->SetOption("text");
	h2DataTransmission->GetXaxis()->SetTitle("Channel");
	h2DataTransmission->GetYaxis()->SetTitle("Hit Difference");
	h2DataTransmission->GetYaxis()->SetBinLabel(1, "g&f");
	h2DataTransmission->GetYaxis()->SetBinLabel(2, "g&!f");
	h2DataTransmission->GetYaxis()->SetBinLabel(3, "!g&f");
	h2DataTransmission->GetYaxis()->SetBinLabel(4, "!g&!f");
	AddHistogram(h2DataTransmission);

	histname = fName + "_QuotMaxAmpGanSADCVsChannel";
	histtitle = fName + " (SADC_maxAmp / Gandalf_maxAmp) vs channel";
	h2MaxAmpQuotGanSADCVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 100, 0, 1);
	h2MaxAmpQuotGanSADCVsChannel->SetOption("colz");
	h2MaxAmpQuotGanSADCVsChannel->GetXaxis()->SetTitle("Channel");
	h2MaxAmpQuotGanSADCVsChannel->GetYaxis()->SetTitle("Quotient");
	AddHistogram(h2MaxAmpQuotGanSADCVsChannel);

	for(int i = 0;i<N_CHANNEL;i++) {
		histname = fName + "_MaxAmpSADCVsMaxAmpGan_";
		histname +=i;
		histtitle = fName + " SADC_maxAmp vs Gandalf_maxAmp channel ";
		histtitle +=i;
		h2MaxAmpSADCVsGandalf[i]=new TH2F(histname,histtitle,200,0,3000, 200, 0, 1000);
		h2MaxAmpSADCVsGandalf[i]->SetOption("colz");
		h2MaxAmpSADCVsGandalf[i]->GetXaxis()->SetTitle("Gandalf Amplitude");
		h2MaxAmpSADCVsGandalf[i]->GetYaxis()->SetTitle("SADC Amplitude");
		AddHistogram(h2MaxAmpSADCVsGandalf[i]);

		histname = fName + "_IntegralSADCVsMaxAmpGan_";
		histname +=i;
		histtitle = fName + " SADC_Integral vs Gandalf_maxAmp channel ";
		histtitle +=i;
		h2IntegralSADCVsGandalf[i]=new TH2F(histname,histtitle,200,0,3000, 200, 0, 4000);
		h2IntegralSADCVsGandalf[i]->SetOption("colz");
		h2IntegralSADCVsGandalf[i]->GetXaxis()->SetTitle("Gandalf Amplitude");
		h2IntegralSADCVsGandalf[i]->GetYaxis()->SetTitle("SADC Integral");
		AddHistogram(h2IntegralSADCVsGandalf[i]);

		histname = fName + "_QuotMaxAmpGanSADCVsMaxAmpGan_";
		histname +=i;
		histtitle = fName + " (SADC_maxAmp / Gandalf_maxAmp) vs Gandalf_maxAmp channel ";
		histtitle +=i;
		h2MaxAmpDiffQuotGanSADCVsGandalf[i]=new TH2F(histname,histtitle,200,0,3000, 200, 0, 0.5);
		h2MaxAmpDiffQuotGanSADCVsGandalf[i]->SetOption("colz");
		h2MaxAmpDiffQuotGanSADCVsGandalf[i]->GetXaxis()->SetTitle("Gandalf Amplitude");
		h2MaxAmpDiffQuotGanSADCVsGandalf[i]->GetYaxis()->SetTitle("Quotient");
		AddHistogram(h2MaxAmpDiffQuotGanSADCVsGandalf[i]);

		histname = fName + "_QuotMaxAmpGanSADCVsMaxAmpSADC_";
		histname +=i;
		histtitle = fName + " (SADC_maxAmp / Gandalf_maxAmp) vs SADC_maxAmp channel ";
		histtitle +=i;
		h2MaxAmpDiffQuotGanSADCVsSADC[i]=new TH2F(histname,histtitle,200,0,1500, 200, 0, 0.5);
		h2MaxAmpDiffQuotGanSADCVsSADC[i]->SetOption("colz");
		h2MaxAmpDiffQuotGanSADCVsSADC[i]->GetXaxis()->SetTitle("SADC Amplitude");
		h2MaxAmpDiffQuotGanSADCVsSADC[i]->GetYaxis()->SetTitle("Quotient");
		AddHistogram(h2MaxAmpDiffQuotGanSADCVsSADC[i]);
	}

	histname = fName + "_IntgralDiffGanSADCVsChannel";
	histtitle = fName + " difference in Integral (I_g - I_sadc) vs channel";
	h2IntegralDiffGanSADCVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 200, -2000, 50000);
	h2IntegralDiffGanSADCVsChannel->SetOption("colz");
	h2IntegralDiffGanSADCVsChannel->GetXaxis()->SetTitle("Channel");
	h2IntegralDiffGanSADCVsChannel->GetYaxis()->SetTitle("Integral Difference");
	AddHistogram(h2IntegralDiffGanSADCVsChannel);

	histname = fName + "_MaxAmpDiffGanSADCVsChannel";
	histtitle = fName + " difference in maxAmp (MaxAmp_g - MaxAmp_sadc) vs channel";
	h2MaxAmpDiffGanSADCVsChannel=new TH2F(histname,histtitle,N_CHANNEL,0,N_CHANNEL, 200, -1000, 3000);
	h2MaxAmpDiffGanSADCVsChannel->SetOption("colz");
	h2MaxAmpDiffGanSADCVsChannel->GetXaxis()->SetTitle("Channel");
	h2MaxAmpDiffGanSADCVsChannel->GetYaxis()->SetTitle("MaxAmp Difference");
	AddHistogram(h2MaxAmpDiffGanSADCVsChannel);

#endif

	//////////////////////////////////////////////////////////////
	/////////////////// PlaneGandalf plots for experts

	if (FRAME_PLOTS && fExpertHistos) {
		for(int i=0;i<N_CHANNEL;i++) {
			stringstream name;
			name<< fName << "g_" << i;
			g_planes[i] = new PlaneGandalf(name.str().c_str());
			g_planes[i]->Init(tree);
			for( vector<TH1*>::iterator it = g_planes[i]->GetHistoList().begin(); it != g_planes[i]->GetHistoList().end();it++)
				AddHistogram((*it));
		}
	}

}


void PlaneCamera::Reset(void)
{
	Plane::Reset();
	chanDigits.clear();
	chanFrames.clear();
	chanf1Digits.clear();
	chansadcDigits.clear();

	if (0 && fExpertHistos) {
		for (map<int,PlaneGandalf*>::iterator plane=g_planes.begin();plane!=g_planes.end();++plane)
			(*plane).second->Reset();
	}
}


void PlaneCamera::StoreDigit(CS::Chip::Digit* digit)
{
	//test if we have the expected digit type
	CS::ChipGandalf::DigitGADC *digitG = dynamic_cast<CS::ChipGandalf::DigitGADC*>(digit);// check if digit is from gandalf and not somewhere else
	CS::ChipF1::Digit *digitF1= NULL;
	CS::ChipSADC::Digit *digitSADC= NULL;

#ifdef PROTO
	digitF1 = dynamic_cast<CS::ChipF1::Digit*>(digit);// check if digit is from gandalf and not somewhere else

	digitSADC = dynamic_cast<CS::ChipSADC::Digit*>(digit);// check if digit is from gandalf and not somewhere else

	if (!digitG && !digitF1 && !digitSADC)
#else
	if (!digitG)
#endif

	{
		std::cerr<<"PlaneCamera::StoreDigit (" << GetName()<< "): a digit is not a Gandalf one, strange...\n";
		return;
	}
	//TODO: sanity checks

	if (digitG) {
		//add digit to digits list
		lDigits.push_back(digit);

		//add digit to channel only if its not the frame digit
		if (digitG->getNumSamples()==0)
			chanDigits[digitG->getChannel()].push_back(digitG);
		else
			chanFrames[digitG->getChannel()].push_back(digitG);

		//increment hit count
		fNhits++;

		if (FRAME_PLOTS && fExpertHistos) {
			if (g_planes.find(digitG->getChannel())==g_planes.end())
				g_planes[digitG->getChannel()] = new PlaneGandalf("Unknown");
			g_planes[digitG->getChannel()]->StoreDigit(digit);
		}
	} else if ( digitF1 && (digitF1->GetTimeDecoded()/digitF1->GetTimeUnit() < -16500) && (digitF1->GetTimeDecoded()/digitF1->GetTimeUnit() > -16850) ) {
		chanf1Digits[digitF1->GetChannel()].push_back(digitF1);
	} else if ( digitSADC) {
		chansadcDigits[digitSADC->GetX()].push_back(digitSADC);
	}
}

void PlaneCamera::EndEvent(const CS::DaqEvent &event)
{

#ifdef PROTO
	// comprare f1 number of send data to gandalf
	for(int i = 0 ; i< N_CHANNEL; i++) {
		bool gan_send = chanDigits.find(i) != chanDigits.end();
		bool f1_send = chanf1Digits.find(i) != chanf1Digits.end();
		if (gan_send && f1_send)
			h2DataTransmission->Fill(i,0);
		if (gan_send && !f1_send)
			h2DataTransmission->Fill(i,1);
		if (!gan_send && f1_send)
			h2DataTransmission->Fill(i,2);
		if (!gan_send && !f1_send)
			h2DataTransmission->Fill(i,3);
	}
#endif

	for (map<int, vector<CS::ChipGandalf::DigitGADC*> >::const_iterator i=chanDigits.begin();i!=chanDigits.end();++i)
	{
		// fill total number of send data for each channel
		h1NumberOfSentData->Fill((*i).first);
		// fill multiplicity histogram
		h2MultiVsChannel->Fill((*i).first, (*i).second.size() );

#ifdef PROTO
		if (chanf1Digits.find((*i).first) != chanf1Digits.end()) {
			h2MultiDiffGanF1VsChannel->Fill((*i).first, (*i).second.size() - chanf1Digits[(*i).first].size() );
		}
#endif

		// fill hit informations in histograms
		for(vector<CS::ChipGandalf::DigitGADC*>::const_iterator hit=(*i).second.begin();hit != (*i).second.end();++hit ) {
			h2AmplitudeMaximumVsChannel->Fill((*i).first,(*hit)->getMaxAmplitude());
			h2IntegralVsChannel->Fill((*i).first,(*hit)->getIntegral());
			h2BaselineVsChannel->Fill((*i).first,(*hit)->getBaseLine());
			h2HighResVsChannel->Fill((*i).first,(*hit)->getHiResTime());
			h2TimeVsChannel->Fill((*i).first,(*hit)->GetTimeDecoded());
			if ( (*hit)->getOpMode()==CS::ChipGandalf::GADC_DEBUG
					|| (*hit)->getOpMode()==CS::ChipGandalf::GADC_DEBUG_IL)
				h2FrameTimeVsChannel->Fill((*i).first, (*hit)->getFrameTime() );

#ifdef PROTO
			// compare to f1 data
			if (chanf1Digits.find((*i).first) != chanf1Digits.end()) {
				for(digitsf1::const_iterator f1hit=chanf1Digits[(*i).first].begin();f1hit != chanf1Digits[(*i).first].end();++f1hit ) {
					h2TimeDiffGanF1VsChannel->Fill((*i).first, ((*hit)->GetTimeDecoded()-(*f1hit)->GetTimeDecoded()) - 1.2928415378e11);
					std::cout << (*hit)->GetTimeDecoded() << std::endl;
					//std::cin.get();
				}
			}
			// compare to sadc data
			if (chansadcDigits.find((*i).first) != chansadcDigits.end()) {
				for(digitsSADC::const_iterator sadchit=chansadcDigits[(*i).first].begin();sadchit != chansadcDigits[(*i).first].end();++sadchit ) {
					double sadc_integral = 0,baseline=0;
					double sadc_maxampl = 0;
					int count = 0;
					for (std::vector<uint16>::iterator  sample = (*sadchit)->GetSamples().begin();sample != (*sadchit)->GetSamples().end();++sample){
						//first 5 samples make baseline
						if (count < 5) baseline+=(*sample);
						if ((*sample) > sadc_maxampl) sadc_maxampl = (*sample);
						sadc_integral += (*sample);
						count++;
					}
					baseline /= 5;
					sadc_integral -= (*sadchit)->GetSamples().size()*baseline;
					sadc_maxampl -= baseline;

					h2IntegralDiffGanSADCVsChannel->Fill((*i).first,(*hit)->getIntegral()-sadc_integral);
					h2MaxAmpDiffGanSADCVsChannel->Fill((*i).first,(*hit)->getMaxAmplitude()-sadc_maxampl);
					h2MaxAmpQuotGanSADCVsChannel->Fill((*i).first,sadc_maxampl/(double)(*hit)->getMaxAmplitude());
					h2MaxAmpSADCVsGandalf[(*i).first]->Fill((double)(*hit)->getMaxAmplitude(),sadc_maxampl);
					h2IntegralSADCVsGandalf[(*i).first]->Fill((double)(*hit)->getMaxAmplitude(),sadc_integral);
					h2MaxAmpDiffQuotGanSADCVsGandalf[(*i).first]->Fill((double)(*hit)->getMaxAmplitude(),sadc_maxampl/(double)(*hit)->getMaxAmplitude());
					h2MaxAmpDiffQuotGanSADCVsSADC[(*i).first]->Fill(sadc_maxampl,sadc_maxampl/(double)(*hit)->getMaxAmplitude());
				}
			}
#endif
		}
	}

	for (map<int, vector<CS::ChipGandalf::DigitGADC*> >::const_iterator i=chanFrames.begin();i!=chanFrames.end();++i)
	{
		// fill frame information in histograms
		for(vector<CS::ChipGandalf::DigitGADC*>::const_iterator hit=(*i).second.begin();hit != (*i).second.end();++hit )
			h2BaselineCompVsChannel->Fill((*i).first,PlaneGandalf::getBaseline( (*hit),1,3,0,1 ) );
	}

	if (0 && fExpertHistos) {
		for (int i=0;i<N_CHANNEL;i++)
			g_planes[i]->EndEvent(event);
	}
}


///////////////////////////////////////////////////////////////////////////////
//////////////////////////// Helper functions

PlaneCamera::Type PlaneCamera::getType() const {
	if( strncmp(GetName(), "CA_Au", 5 ) == 0 ) return PlaneCamera::Aup;
	if( strncmp(GetName(), "CA_Bu", 5 ) == 0 ) return PlaneCamera::Bup;
	if( strncmp(GetName(), "CA_Ad", 5 ) == 0 ) return PlaneCamera::Adown;
	if( strncmp(GetName(), "CA_Bd", 5 ) == 0 ) return PlaneCamera::Bdown;
	return PlaneCamera::UNKNOWN;
}


///////////////////////////////////////////////////////////////////////////////
//////////////////////////// Calibration Database


#if USE_DATABASE == 1
void PlaneCamera::ReadCalib(const tm &t)
{
  calibTime = mktime((tm *)&t);

  std::cout<<"PlaneCamera::ReadCalib() ==> "<<this->GetName()<<" reading calibrations !"<<std::endl;
  // read-in corresponding calibration constants
  try{
    ReadFromDataBase(calib_data,t);

    if(calib_data.size() <  (unsigned) N_CHANNEL) {
          std::cerr<<"Size of Calibration File is not correct ! Should be : "
      	    <<N_CHANNEL<<" Is "<<calib_data.size()<<" "
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
    std::cout<<"PlaneCamera::ReadCalib() ==> "<<GetName()
	<<" calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	<<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	<<", not found in DB"<<std::endl;
  }

}

#endif //USE_DATABASE

