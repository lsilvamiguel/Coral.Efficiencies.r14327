#include "PlaneGandalf.h"
#include "PlanePanel.h"
#include "ChipGandalf.h"
#include "Reference.h" //needed for TH1F_Ref etc.

#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <stdio.h>
#include <math.h>

#include <stdlib.h>

ClassImp(PlaneGandalf);

double hitsPerEvent;

PlaneGandalf::PlaneGandalf(const char *detname)
	: Plane(detname)
{
	time_min=-1;time_max=-1;
	timeDec_min=-1;timeDec_max=-1;
	firstEvent=-1;
}


void PlaneGandalf::Init(TTree* tree)
{

	//////////////// constants for the CFD
	////////////////////////////////////////////////////////////

	// weather to use custom constants for cfd or the ones gandalf transmits
	USE_CUSTOM_CONSTANTS = false;

	FRACTFACT	= 0.6;
	DELAY		= 1;
	THRESHOLD	= 20;

	// sample points before CFD zero crossing
	CFD_OFFSET	= 25;
	// number of sambles per pulse
	PULSE_WIDTH	= 50;

	BASELINE 	= 200;


	//reset event counter
	fRateCounter = 0;

	std::string name;

	//add histos to histogram list

	// !Hits per Event
	name = fName + "_hits_per_event";
	fHpEvtHist = new TH1F_Ref(name.c_str(), name.c_str(), 20, 0.0, 20.0, fRateCounter);
	((TH1F_Ref*)fHpEvtHist)->SetReference(fReferenceDirectory);
	AddHistogram(fHpEvtHist);	
	
	name = fName + "_hits_per_event_computed";
	fHpEvtHist_computed = new TH1F_Ref(name.c_str(), name.c_str(), 20, 0.0, 20.0, fRateCounter);
	((TH1F_Ref*)fHpEvtHist_computed)->SetReference(fReferenceDirectory);
	AddHistogram(fHpEvtHist_computed);

	// Max Amplitude histogram
	name = fName + "_Max_Amplitude";
	fAmpHist = new TH1F_Ref(name.c_str(), name.c_str(), 4096, 0.0, 4096.0,fRateCounter);
	fAmpHist->GetXaxis()->SetTitle("ADC LSB");
	((TH1F_Ref*)fAmpHist)->SetReference(fReferenceDirectory);
	AddHistogram(fAmpHist);
	
	// Integral histogram
	name = fName + "_Integral";
	fIntHist = new TH1F_Ref(name.c_str(), name.c_str(), 10000, 0.0, 100000.0,fRateCounter);
	((TH1F_Ref*)fIntHist)->SetReference(fReferenceDirectory);
	AddHistogram(fIntHist);
	
	// Hi Res Time histogram
	name = fName + "_High_Res_Time";
	fHRTHist = new TH1F_Ref(name.c_str(), name.c_str(), 1024, 0.0, 1024.0, fRateCounter);
	fHRTHist->GetXaxis()->SetTitle("1/1024 CT LSB");
	((TH1F_Ref*)fHRTHist)->SetReference(fReferenceDirectory);
	AddHistogram(fHRTHist);

	// Hi Res Time histogram computed
	name = fName + "_High_Res_Time_computed";
	fHRTHist_computed = new TH1F_Ref(name.c_str(), name.c_str(), 1024, 0.0, 1024.0, fRateCounter);
	fHRTHist_computed->GetXaxis()->SetTitle("1/1024 CT LSB");
	((TH1F_Ref*)fHRTHist_computed)->SetReference(fReferenceDirectory);
	AddHistogram(fHRTHist_computed);

	// !Frametimes
	name = fName + "_Frame_Time";
	fFrameTHist = new TH1F_Ref(name.c_str(), name.c_str(), 100, 0.0, 200.0, fRateCounter);
	fFrameTHist->GetXaxis()->SetTitle("n'th sample");
	((TH1F_Ref*)fFrameTHist)->SetReference(fReferenceDirectory);
	AddHistogram(fFrameTHist);

	// !Frametimes
	name = fName + "_Frame_Time_computed";
	fFrameTHist_computed = new TH1F_Ref(name.c_str(), name.c_str(), 100, 0.0, 200.0, fRateCounter);
	fFrameTHist_computed->GetXaxis()->SetTitle("n'th sample");
	((TH1F_Ref*)fFrameTHist_computed)->SetReference(fReferenceDirectory);
	AddHistogram(fFrameTHist_computed);

	// !Decoded Times
	name = fName + "_TimeDecoded";
	fTimeDecoded = new TH1F_Ref(name.c_str(), name.c_str(), 1000, 0, 0, fRateCounter);
	fTimeDecoded->GetXaxis()->SetTitle("ns");
	((TH1F_Ref*)fTimeDecoded)->SetReference(fReferenceDirectory);
	AddHistogram(fTimeDecoded);

    // !Decoded Times with fixed range
    name = fName + "_TimeDecoded_cut";
    fTimeDecodedCut = new TH1F_Ref(name.c_str(), name.c_str(), 3000, 2500.0, 3500.0, fRateCounter);
    fTimeDecodedCut->GetXaxis()->SetTitle("ns");
    ((TH1F_Ref*)fTimeDecodedCut)->SetReference(fReferenceDirectory);
    AddHistogram(fTimeDecodedCut);


	// !BaseLine
	name = fName + "_BaseLine";
	fBaseLine = new TH1F_Ref(name.c_str(), name.c_str(), 200, 0.0, 2000.0, fRateCounter);
	fBaseLine->GetXaxis()->SetTitle("ADC LSB");
	((TH1F_Ref*)fBaseLine)->SetReference(fReferenceDirectory);
	AddHistogram(fBaseLine);

	// !BaseLine
	name = fName + "_BaseLine_computed";
	fBaseLine_computed = new TH1F_Ref(name.c_str(), name.c_str(), 200, 0.0, 2000.0, fRateCounter);
	fBaseLine_computed->GetXaxis()->SetTitle("ADC LSB");
	((TH1F_Ref*)fBaseLine_computed)->SetReference(fReferenceDirectory);
	AddHistogram(fBaseLine_computed);

	//Frame histogram
	name = fName + "_Frame_Histo";
	fFrameHist = new TH1F_Ref(name.c_str(), name.c_str(), 100, 0.0, 99,fRateCounter);
	fFrameTHist->GetXaxis()->SetTitle("n'th sample");
	((TH1F_Ref*)fFrameHist)->SetReference(fReferenceDirectory);
	AddHistogram(fFrameHist);
		
	//Plotting histogram last event
	name = fName + "_Frame_Plot";
	fFramePlotHist = new TH1F_Ref( name.c_str(), name.c_str(), 100, 0.0, 100.0, fRateCounter );
	fFramePlotHist->GetXaxis()->SetTitle("n'th sample");
	fFramePlotHist->GetXaxis()->SetTitle("ADC LSB");
	((TH1F_Ref*)fFramePlotHist)->SetReference(fReferenceDirectory);
	AddHistogram(fFramePlotHist);
	
	//Plotting CFD histogram last event
	name = fName + "_CFD_Plot";
	fCFDPlotHist = new TH1F_Ref( name.c_str(), name.c_str(), 100, 0.0, 100.0, fRateCounter );
	fCFDPlotHist->GetXaxis()->SetTitle("n'th sample");
	((TH1F_Ref*)fCFDPlotHist)->SetReference(fReferenceDirectory);
	AddHistogram(fCFDPlotHist);

	//Plotting histogram last event above thresh
	name = fName + "_Frame_Plot_above_thresh";
	fFramePlotHistThresh = new TH1F_Ref( name.c_str(), name.c_str(), 100, 0.0, 100.0, fRateCounter );
	fFramePlotHistThresh->GetXaxis()->SetTitle("n'th sample");
	((TH1F_Ref*)fFramePlotHistThresh)->SetReference(fReferenceDirectory);
	AddHistogram(fFramePlotHistThresh);

	//Plotting CFD histogram last event above thresh
	name = fName + "_CFD_Plot_above_thresh";
	fCFDPlotHistThresh = new TH1F_Ref( name.c_str(), name.c_str(), 100, 0.0, 100.0, fRateCounter );
	fCFDPlotHistThresh->GetXaxis()->SetTitle("n'th sample");
	((TH1F_Ref*)fCFDPlotHistThresh)->SetReference(fReferenceDirectory);
	AddHistogram(fCFDPlotHistThresh);


	//Plotting histogram last event with hit
	name = fName + "_Frame_Plot_hit";
	fFramePlotHistLast = new TH1F_Ref( name.c_str(), name.c_str(), 100, 0.0, 100.0, fRateCounter );
	fFramePlotHistLast->GetXaxis()->SetTitle("n'th sample");
	((TH1F_Ref*)fFramePlotHistLast)->SetReference(fReferenceDirectory);
	AddHistogram(fFramePlotHistLast);

	//Plotting CFD histogram last event with hit
	name = fName + "_CFD_Plot_hit";
	fCFDPlotHistLast = new TH1F_Ref( name.c_str(), name.c_str(), 100, 0.0, 100.0, fRateCounter );
	fCFDPlotHistLast->GetXaxis()->SetTitle("n'th sample");
	((TH1F_Ref*)fCFDPlotHistLast)->SetReference(fReferenceDirectory);
	AddHistogram(fCFDPlotHistLast);

	//integral_computed
	name = fName + "_Integral_computed";
	integral_computed = new TH1F_Ref(name.c_str(), name.c_str(), 10000, 0.0, 100000.0, fRateCounter );
	((TH1F_Ref*)integral_computed)->SetReference(fReferenceDirectory);
	AddHistogram(integral_computed);

	//amplitude_computed
	name = fName + "_Max_Amplitude_computed";
	fAmpHist_computed = new TH1F_Ref(name.c_str(), name.c_str(), 4096, 0.0, 4096.0, fRateCounter );
	((TH1F_Ref*)fAmpHist_computed)->SetReference(fReferenceDirectory);
	AddHistogram(fAmpHist_computed);

	//diff integral_computed to integral_gandalf
	name = fName + "_diff_integral_computed_to_integral_gandalf";
	fDiffIntegral = new TH1F_Ref(name.c_str(), name.c_str(), 2058, -1024.0, 1024.0, fRateCounter );
	((TH1F_Ref*)fDiffIntegral)->SetReference(fReferenceDirectory);
	AddHistogram(fDiffIntegral);

	name = fName + "_diff_amplitude_computed_to_amplitude_gandalf";
	fDiffAmp = new TH1F_Ref(name.c_str(), name.c_str(), 1024, -512.0, 512.0, fRateCounter );
	((TH1F_Ref*)fDiffAmp)->SetReference(fReferenceDirectory);
	AddHistogram(fDiffAmp);

	name = fName + "_diff_highres_computed_to_highres_gandalf";
	fDiffHighRes = new TH1F_Ref(name.c_str(), name.c_str(), 2048, -1024, 1024, fRateCounter );
	((TH1F_Ref*)fDiffHighRes)->SetReference(fReferenceDirectory);
	AddHistogram(fDiffHighRes);

	name = fName + "_diff_frametime_computed_to_frametime_gandalf";
	fDiffFrameTime = new TH1F_Ref(name.c_str(), name.c_str(), 400, -200, 200, fRateCounter );
	((TH1F_Ref*)fDiffFrameTime)->SetReference(fReferenceDirectory);
	AddHistogram(fDiffFrameTime);

	name = fName + "_diff_baseline_computed_to_baseline_gandalf";
	fDiffBaseline = new TH1F_Ref(name.c_str(), name.c_str(), 400, -200, 200, fRateCounter );
	((TH1F_Ref*)fDiffBaseline)->SetReference(fReferenceDirectory);
	AddHistogram(fDiffBaseline);


//	name = fName + "_time_vs_event";
//	fTimeVSevent = new TH2D(name.c_str(), name.c_str(), 1000, 0, 0, 1000, 0 , 0);
//	fTimeVSevent->SetOption("colz");
//	AddHistogram(fTimeVSevent);
//
//	name = fName + "_timeDec_vs_event";
//	fTimeDecVSevent = new TH2D(name.c_str(), name.c_str(), 1000, 0, 0, 1000, 0 , 0);
//	fTimeDecVSevent->SetOption("colz");
//	AddHistogram(fTimeDecVSevent);

	if (fExpertHistos) {
		//statistics of occurences of bits in frame words
		name = fName + "_bit_statistic";
		fBitStat = new TH1F_Ref(name.c_str(), name.c_str(), 14, 0, 14, fRateCounter );
		((TH1F_Ref*)fBitStat)->SetReference(fReferenceDirectory);
		AddHistogram(fBitStat);

		//statistics of occurences of bits in frame words
		name = fName + "_bit_statistic_01";
		fBitStat_01 = new TH1F_Ref(name.c_str(), name.c_str(), 28, 0, 14, fRateCounter );
		((TH1F_Ref*)fBitStat_01)->SetReference(fReferenceDirectory);
		AddHistogram(fBitStat_01);

		name = fName + "_overlayed_pulses_half_sample";
		fOverlayedPulsesHS = new TH2F(name.c_str(), name.c_str(), 500, 0, 50, 400, 0, 4000 );
		fOverlayedPulsesHS->SetOption("colz");
		AddHistogram(fOverlayedPulsesHS);

	} else {
		name = fName + "_overlayed_pulses";
		fOverlayedPulses = new TH2F(name.c_str(), name.c_str(), 500, 0, 50, 400, 0, 4000 );
		fOverlayedPulses->SetOption("colz");
		AddHistogram(fOverlayedPulses);
	}




}

void PlaneGandalf::ResetHistograms()
{
	Plane::ResetHistograms();
}


struct cfd_info {

	vector<double>  highres;
	vector<int>  frametime;
	int maxamp;

};

///////////////////////////////////////////////////////////////////



///////////////// other constants
int				HIT_THRESHOLD	= 50; //threshold for basic above thresh hit detection

// if more than xxx% of data is equal, then use this value as baseline
double			QUOT_TRIGGER	= 0.7;


void PlaneGandalf::StoreDigit(CS::Chip::Digit* digit)
{
	//test if we have the expected digit type
	CS::ChipGandalf::DigitGADC *digitG = dynamic_cast<CS::ChipGandalf::DigitGADC*>(digit);

	if (!digitG)
	{
		std::cerr<<"PlaneGandalf::StoreDigit (" << GetName()<< "): a digit is not a Gandalf one, strange...\n";
		return;
	}

	//add digit to digits list
	lDigits.push_back(digit);

	//increment hit count
	fNhits++;

	channel = digitG->getChannel();

	// overwrite CFD constants with the one gandalf transmits
	if (!USE_CUSTOM_CONSTANTS && (digitG->getOpMode()==CS::ChipGandalf::GADC_DEBUG
			|| digitG->getOpMode()==CS::ChipGandalf::GADC_DEBUG_IL) && digitG->getNumSamples()==0 ) {
		if (digitG->getFracfact()!=0)
			FRACTFACT  	= 1/(double)(pow(2,digitG->getFracfact()));
		else FRACTFACT  = 1;
		DELAY			= digitG->getDelay();
		THRESHOLD		= digitG->getThreshold();
	}

}

/**
 * Calculates CFD out of a Gandalf Frame Digit
 * start 		the first sample point
 * increment 	the increment for the samples loop (used to check normal mode while having interleaved data)
 * integral		the integral calculated for the frame, used to normalize the 'overlayed pulses' plot
 */
cfd_info* PlaneGandalf::calc_cfd(CS::ChipGandalf::DigitGADC* digitGADC, int start, int increment, float integral, float baseline) {

	cfd_info *info = new cfd_info();


	for (int i = start; i <= digitGADC->getNumSamples()-(increment); i+=increment)
	{
		// CFD Algorythm
		if(i>2+(DELAY*increment)) {
			double k = FRACTFACT *(digitGADC->getSample(i-(increment * 2)) - BASELINE) - (digitGADC->getSample( i - (DELAY*increment) - (increment * 2) ) - BASELINE);
			double m = FRACTFACT *(digitGADC->getSample(i-(increment * 1)) - BASELINE) - (digitGADC->getSample( i - (DELAY*increment) - (increment * 1) ) - BASELINE);
			double n = FRACTFACT *(digitGADC->getSample(i) - BASELINE) - (digitGADC->getSample( i- (DELAY*increment) ) - BASELINE);

			if (increment==1) fCFDPlotHist->Fill(i,m);

			// detected hit
			if (  m > 0 && n < 0 && (m > THRESHOLD || n < -THRESHOLD) ){
				// times 1024 to fit gandalf output
				info->highres.push_back(1024 * m/(m-n));
				info->frametime.push_back(i - increment);
				if (increment==1) {
					fHRTHist_computed->Fill(info->highres[info->highres.size()-1]);
					fFrameTHist_computed->Fill(info->frametime[info->frametime.size()-1]);
				}
			} else if (m==0  && k > 0 && n < 0 && (k > THRESHOLD || n < -THRESHOLD) ) {
				info->highres.push_back(0);
				info->frametime.push_back(i - increment);
				if (increment==1) {
					fHRTHist_computed->Fill(info->highres[info->highres.size()-1]);
					fFrameTHist_computed->Fill(info->frametime[info->frametime.size()-1]);
				}
			}
		}
	}

	// we have (at least) a hit
	if (info->frametime.size()>0) {
		for(unsigned int l=0;(l<info->frametime.size()&&l<info->highres.size());l++) {
			for(int j=1;j<PULSE_WIDTH*increment && (info->frametime[l]-(CFD_OFFSET*increment)+j < digitGADC->getNumSamples());j+=increment)
				if (fExpertHistos)
					fOverlayedPulsesHS->Fill( (j+1)/increment - info->highres[l]/(double)1024, 10000 * (digitGADC->getSample( info->frametime[l]-(CFD_OFFSET*increment)+j) - baseline) / (double)integral);
				else
					fOverlayedPulses->Fill( (j+1)/increment - info->highres[l]/(double)1024, 10000 * (digitGADC->getSample( info->frametime[l]-(CFD_OFFSET*increment)+j) - baseline) / (double)integral);

		}
	}

	return info;
}

/**
 * Baseline calculation ala Tobi Baumann
 *
 * uint32 	sigma		bins taken into account for middle value (default 1)
 * uint32 	deltaBin	width of one bin (default 3)
 *
 */
float PlaneGandalf::getBaseline (
		CS::ChipGandalf::DigitGADC* digitGADC,
		uint32 sigma=1,
		uint32 deltaBin=3,
		uint32 startSample=0,
		uint32 deltaSample=1 ) {

	float baseline = 0; //return value
	uint32 minBin = 4096;
	uint32 maxBin = 0;

	// if more than xxx% of data is equal, then use this value as baseline
	const float quotTrigger = QUOT_TRIGGER;

	//set minimum bin and maximum bin
	for(int32 i=startSample;i<digitGADC->getNumSamples();i+=deltaSample) {
		if (digitGADC->getSample(i)>maxBin) maxBin = digitGADC->getSample(i);
		if (digitGADC->getSample(i)<minBin) minBin = digitGADC->getSample(i);
	}

	//set the number of bins
	int32 numBins = maxBin + 2*sigma*deltaBin;

	//array with number of bin entries
	int32* binsEntries = new int32[numBins];

	//array with array of entries
	uint16** bins = new uint16*[numBins];

	for(CS::int32 i = 0; i < numBins; i++) {
		binsEntries[i] = 0;
		bins[i] = new uint16[digitGADC->getNumSamples()];
	}

	//round data and fill them in related bin
	for(int32 i = startSample; i < digitGADC->getNumSamples(); i += deltaSample) {

		//round data and increase number of bin entry
		uint16 val = (uint16)round(digitGADC->getSample(i),deltaBin);
		binsEntries[val]++;

		//position of actual data point
		int32 pos = binsEntries[val] - 1;

		//std::cout << "start: "<< startSample << " deltaSample: " << deltaSample << " i: " << i << std::endl;

		//fill bin
		bins[val][pos] = digitGADC->getSample(i);
	}

	//get position of maximum filled bin
	int32 maxPos=0;
	for(int i=0;i<numBins;i++)
		if (binsEntries[i]>maxPos) maxPos = i;

	float sum = 0;
	int32 numDivBins = 0;

	float quot = (float)binsEntries[maxPos] / (digitGADC->getNumSamples()/(double)deltaSample);

//	std::cout << "min: " << minBin << " max: " << maxBin << " numBins: " << numBins << std::endl;
//	std::cout << "maxPos: " << maxPos << " entries: " << binsEntries[maxPos] << " quot: " << quot  << std::endl;


	//check if data is constant, and use only that bin if so
	if (quot > quotTrigger) {

		for (int i=0;i<binsEntries[maxPos];i++)
			sum+=bins[maxPos][i];

		numDivBins += binsEntries[maxPos];

	} else {

		//check that summing range is valid
		int32 start = (maxPos - (int32)sigma >= 0) ? maxPos - (int32)sigma : 0;
		int32 end = (maxPos + (int32)sigma <= numBins) ? maxPos + (int32)sigma : numBins;

		for(int32 i = start; i <= end; i++) {

			if (binsEntries[i] > 0) {
				for (int j=0;j<binsEntries[i];j++)
					sum+=bins[i][j];

				numDivBins += binsEntries[i];
			}

		}
	}

	//delete all allocated heap data
	delete[] binsEntries;

	for(int32 i = 0; i < numBins; i++) {
		delete[] bins[i];
	}

	delete[] bins;

	// something has gone wrong...
	if (sum==0||numDivBins==0) return -1;

	baseline = sum / numDivBins;

	return baseline;
}

void PlaneGandalf::EndEvent(const CS::DaqEvent &event)
{


	unsigned int 	gintegral = 0;
	float			integral = 0;
	int				gmax_amp=-1;
	float 			max_amp=-1;
	int 			num_samples=-1;
	int 			gbaseline=-1;
	float 			baseline=-1;
	vector<double>  ghighres,gframetime;

	int 			GANDALF_FRAMETIME_OFFSET = 0;

	cfd_info* 		cfdinfo = NULL;
	isDebug					= false;

	if (thr_flag) TThread::Lock();

	//fill histograms
	typedef std::list<CS::Chip::Digit*>::iterator lDIter;
	

	//loop over digits list
	for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++)
	{
		//we expect a ChipGandalf::DigitGADC Digit
		CS::ChipGandalf::DigitGADC* digitGADC = dynamic_cast<CS::ChipGandalf::DigitGADC*> (*ii);
		if (!digitGADC)
		{
			std::cerr<<"PlaneGandalf::EndEvent: a digit is not a GandalfADC one, strange...\n";
			continue;
		}

		// we have a digit with frame data
		if (digitGADC->getNumSamples()>0) {

			isDebug=true;

			fFramePlotHist->Reset("M"); //reset minimum and maximum (and content, of course)
			fCFDPlotHist->Reset("M");

			stringstream plot_title;
			plot_title << fCFDPlotHist->GetName();
			plot_title << "___FF_"<<FRACTFACT << "___DELAY_" << DELAY << "___THRESH_" << THRESHOLD;
			fCFDPlotHist->SetTitle(plot_title.str().c_str());

			num_samples = digitGADC->getNumSamples();

			//set #bins according to #frames, also min and max
			fFrameHist->SetBins( num_samples, 0, num_samples-1);
			fFramePlotHist->SetBins( num_samples, 0, num_samples-1);
			fCFDPlotHist->SetBins( num_samples, 0, num_samples-1);
			fFrameTHist->SetBins( num_samples, 0, num_samples-1);
			fFrameTHist_computed->SetBins( num_samples, 0, num_samples-1);

			for (int i = 0; i < num_samples; i++)
			{
				if (max_amp==-1 || max_amp < digitGADC->getSample(i))
					max_amp=digitGADC->getSample(i);
				fFramePlotHist->Fill( i+1, digitGADC->getSample(i) );
				fFrameHist->Fill(  i+1, digitGADC->getSample(i) );
				integral += digitGADC->getSample(i);

				if (fExpertHistos){
					// fill bit 'j' in its bin; if bit = 0, fill 'j', if bit = 1, fill 'j+0.5'
					for(int j=0;j<=13;j++)
					{
						short val= ((digitGADC->getSample(i) >> j )&1);
						fBitStat_01->Fill( j+ 0.5*val );
						if (val==1) fBitStat->Fill( j );
					}
				}
			}
			// calculate the CFD with half sampling (used to test 'normal mode' with interleaved data)
//			if (digitGADC->getOpMode()==CS::ChipGandalf::GADC_DEBUG_IL) calc_cfd(digitGADC,0,2,integral);
//			cfdinfo = calc_cfd(digitGADC,0,1,integral);

			//calculate the Baseline for the event
			baseline = PlaneGandalf::getBaseline(digitGADC,1,3,0,1);

			// we need to subtract the baseline from the calculated values;
			// we have to trust gandalf here that it writes the correct baseline into the digit
			integral-=(num_samples*baseline);
			max_amp-=baseline;

			// calculate the CFD with half sampling (used to test 'normal mode' with interleaved data)
			//if (digitGADC->getOpMode()==CS::ChipGandalf::GADC_DEBUG_IL) calc_cfd(digitGADC,0,2,integral);
			cfdinfo = calc_cfd(digitGADC,0,1,integral,baseline);
		}

		// we have a normal digit containing processed data; if this is a dbg digit, we have the frame time in addition
		else
		{
			if (firstEvent==-1) firstEvent = event.GetEventNumberInRun();

			double time = digitGADC->getTime()*digitGADC->GetTimeUnit();
			if (time_min == -1 || time_min > time ) time_min = time;
			if (time_max == -1 || time_max < time ) time_max = time;
//			if (fExpertHistos) {
//				//fill time vs event
//				fTimeVSevent->Fill(event.GetEventNumberInRun(),time);
//				fTimeVSevent->GetYaxis()->SetLimits(time_min-20,time_max+20);
//				fTimeVSevent->GetXaxis()->SetLimits(firstEvent-10,event.GetEventNumberInRun()+100);
//			}

//			if ( (timeDec_min != -1 && timeDec_max == -1 &&
//					(digitGADC->GetTimeDecoded() < (timeDec_min - 500)
//							|| digitGADC->GetTimeDecoded() > (timeDec_max + 500))) || strncmp(fName.c_str(),"GT830",5) == 0) {
//				std::cout << "Time Decoded differs by more than 500ns; timeDec_min=" << timeDec_min
//						<<";  timeDec_max=" << timeDec_max<< ";  timeDec=" << digitGADC->GetTimeDecoded() << std::endl;
//				std::cout <<fName << "  " << "Spill= " << event.GetBurstNumber() << ";  Event= " << event.GetEventNumberInBurst() << "\n" << std::endl;
//				//cin.get();
//			} else {
			if (timeDec_min == -1 || timeDec_min > digitGADC->GetTimeDecoded() ) timeDec_min = digitGADC->GetTimeDecoded();
			if (timeDec_max == -1 || timeDec_max < digitGADC->GetTimeDecoded() ) timeDec_max = digitGADC->GetTimeDecoded();


			fTimeDecoded->Fill( digitGADC->GetTimeDecoded() );
			fTimeDecoded->GetXaxis()->SetLimits(timeDec_min-20,timeDec_max+20);

//			fTimeDecVSevent->Fill(event.GetEventNumberInRun(),digitGADC->GetTimeDecoded());
//			fTimeDecVSevent->GetYaxis()->SetLimits(timeDec_min-20,timeDec_max+20);
//			fTimeDecVSevent->GetXaxis()->SetLimits(firstEvent-10,event.GetEventNumberInRun()+100);

			fTimeDecodedCut->Fill( digitGADC->GetTimeDecoded() );



			//fill max Amplitude
			fAmpHist->Fill( digitGADC->getMaxAmplitude() );
			if (gmax_amp==-1 || gmax_amp < (int32)digitGADC->getMaxAmplitude())
				gmax_amp = digitGADC->getMaxAmplitude();

			//fill Integral value
			gintegral = digitGADC->getIntegral();
			fIntHist->Fill( gintegral );
		
			gbaseline  = digitGADC->getBaseLine();
			fBaseLine->Fill(gbaseline);

			//fill high res time value
			fHRTHist->Fill( digitGADC->getHiResTime());

			if (digitGADC->getOpMode()==CS::ChipGandalf::GADC_DEBUG ||digitGADC->getOpMode()==CS::ChipGandalf::GADC_DEBUG_IL)
			{
				fFrameTHist->Fill( digitGADC->getFrameTime() - GANDALF_FRAMETIME_OFFSET);
			}

			ghighres.push_back(digitGADC->getHiResTime());
			gframetime.push_back(digitGADC->getFrameTime() - GANDALF_FRAMETIME_OFFSET);

		}

	}

	fNhitsKept=isDebug?lDigits.size()-1:lDigits.size();
	// multiplicity
	fHpEvtHist->Fill(fNhitsKept);

	// only fill integral/amplitude and the differences to the digit, if we really had a hit
	if (isDebug&&fNhitsKept>0) {

		integral_computed->Fill(integral);
		fAmpHist_computed->Fill(max_amp);
		fBaseLine_computed->Fill(baseline);

		fDiffIntegral->Fill(gintegral - integral);
		fDiffAmp->Fill(gmax_amp - max_amp);
		fDiffBaseline->Fill(gbaseline - baseline);


		if (max_amp > HIT_THRESHOLD) {
			stringstream plot_title,plot_name;
			plot_title << fFramePlotHistThresh->GetTitle();
			plot_name << fFramePlotHistThresh->GetName();
			fFramePlotHist->Copy((TObject&)*fFramePlotHistThresh);
			fFramePlotHistThresh->SetTitle(plot_title.str().c_str());
			fFramePlotHistThresh->SetName(plot_name.str().c_str());

			plot_title.str("");plot_name.str("");
			plot_title << fCFDPlotHistThresh->GetName();
			plot_name << fCFDPlotHistThresh->GetName();
			fCFDPlotHist->Copy((TObject&)*fCFDPlotHistThresh);
			plot_title << "___FF_"<<FRACTFACT << "___DELAY_" << DELAY << "___THRESH_" << THRESHOLD;
			fCFDPlotHistThresh->SetTitle(plot_title.str().c_str());
			fCFDPlotHistThresh->SetName(plot_name.str().c_str());
		}

		if (cfdinfo != NULL) {
			fHpEvtHist_computed->Fill(cfdinfo->frametime.size());

			for(unsigned int i=0; (i<cfdinfo->frametime.size()&&i<gframetime.size()) ;i++ )
				fDiffFrameTime->Fill(cfdinfo->frametime[i]-gframetime[i]);

			for(unsigned int i=0; (i<cfdinfo->highres.size()&&i<ghighres.size()) ;i++ )
				fDiffHighRes->Fill(cfdinfo->highres[i]-ghighres[i]);

			// if cfd found a hit, create new plots showing the last detected hit
			if (cfdinfo->frametime.size()>0){
				stringstream plot_title,plot_name;
				plot_title << fCFDPlotHistLast->GetName();
				plot_name << fCFDPlotHistLast->GetName();
				fCFDPlotHist->Copy((TObject&)*fCFDPlotHistLast);
				plot_title << "___FF_"<<FRACTFACT << "___DELAY_" << DELAY << "___THRESH_" << THRESHOLD;
				fCFDPlotHistLast->SetTitle(plot_title.str().c_str());
				fCFDPlotHistLast->SetName(plot_name.str().c_str());

				plot_title.str("");plot_name.str("");
				plot_title << fFramePlotHistLast->GetTitle();
				plot_name << fFramePlotHistLast->GetName();
				fFramePlotHist->Copy((TObject&)*fFramePlotHistLast);
				fFramePlotHistLast->SetTitle(plot_title.str().c_str());
				fFramePlotHistLast->SetName(plot_name.str().c_str());
			}

		}

		// expert plots in case there is a difference in integral and amplitude calculated from samples vs gandalf
		if (fExpertHistos) {
			if (fabs(gintegral - integral)>1){
				stringstream name;
				name << fName << "_Frame_Plot_diff_integral_by_" << (gintegral - integral) << "_evt_" << event.GetEventNumberInRun();
				TH1F *tmp = new TH1F(*fFramePlotHist);
				tmp->SetTitle(name.str().c_str());
				AddHistogram(tmp);
			}
			if (fabs(gmax_amp - max_amp)>1){
				stringstream name;
				name << fName << "_Frame_Plot_diff_amp_by_"<< (gmax_amp - max_amp) << "_evt_" << event.GetEventNumberInRun();
				TH1F *tmp = new TH1F(*fFramePlotHist);
				tmp->SetTitle(name.str().c_str());
				AddHistogram(tmp);
			}
		}
	}

	if (thr_flag) TThread::UnLock();

}

void PlaneGandalf::ControlPanel(const TGWindow* p, const TGWindow* main) 
{
    if (fControlPanel) fControlPanel = new PlanePanel(p, main, 100, 100, this);
}


int32 PlaneGandalf::roundToInt(float x) {

	if (x < 0) {
		return (int32)floor(x - 0.5);
	} else {
		return (int32)floor(x + 0.5);
	}

}

int32 PlaneGandalf::floorToInt(float x) {
	return (int32)floor(x);
}

int32 PlaneGandalf::ceilToInt(float x) {
	return (int32)ceil(x);
}

float PlaneGandalf::round(float x, float step) {

	float tmp;
	int32 rounded;

	tmp = x / step;

	rounded = roundToInt(tmp);

	return rounded * step;
}


