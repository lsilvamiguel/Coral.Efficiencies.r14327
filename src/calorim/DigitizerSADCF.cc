#include <algorithm>

#include "TVirtualFFT.h"
#include "TDirectory.h"
#include "TFile.h"

#include "DigitizerSADCF.h"
#include "CsCalorimeter.h"
#include "CsOpt.h"
#include "CsEvent.h"
#include "CsHistograms.h"


///////////////////////////////////////////////////////////////////////////////
DigitizerSADCF::DigitizerSADCF ( CsCalorimeter* parent, Type type , int ic)
              : DigitizerSADCBase(parent,  type), newtonPrecission_(0.01),
                newtonNIter_(10), numberFourierComponents_(8),
                integrationWindow_(10), parabolaShape_(-1), sampleSize_(0),
                writeTree_(false) {
	mycell_ = ic;
	if( CsOpt::Instance()->getOpt( parent->GetName(), "SADC_SHAPE_FILTER_WRITE_TREE" ) ){
		writeTree_ = true;
	}
}

///////////////////////////////////////////////////////////////////////////////
DigitizerSADCF::~DigitizerSADCF() {}

///////////////////////////////////////////////////////////////////////////////
bool DigitizerSADCF::Check( void ) const {
	// Nothing to Init just check
	if( GetParent() == NULL   ) {
		std::cerr <<" FATAL! DigitizerSADCF::Check(): GetParent() == NULL which is fatal in present code implementation " << std::endl;
		exit(1);
	}
	if( mycell_ < 0 ) {
		std::cerr <<" FATAL! DigitizerSADCF::Check(): mycell_ =" <<  mycell_ <<" but we need for a while this field to be properly initialized " << std::endl;
		exit(1);
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
bool DigitizerSADCF::FitBase( const std::vector<uint16> &sample ) {
	if( debug_ )
		std::cout << " CsDigitizerSADCF::FitMaxBase debug " << std::endl;
	// Reset
	Clear();
	if( debug_ )
		std::cout << " CsDigitizerSADCF::FitMaxBase Clear() OK " << std::endl;
	if( debug_ )
		PrintSettings();

	if(parabolaShape_ < 0) {
		parabolaShape_ = parent_->GetFFTwidth(mycell_);
	}

	// Do basline subtraction and base shape analysis
	GetCcSample( sample );

	sampleSize_ = csample_.size();
	if(hDer_.size() != csample_.size() * 4) {
		hDer_.resize(csample_.size() * 4);
	}
	for(unsigned int i = 0; i < hDer_.size(); ++i)
		hDer_[i] = 0;

	////// FOURIER
	FourierTransform();
	int derMaxJ = 0;
	double derMax = 0.;
	for (unsigned int j = 0; j < hDer_.size(); j++) { // will return the max of the derivative
		if ((hDer_[j]) > derMax) {
			derMax = hDer_[j];
			derMaxJ = j;
		}
	}
	double nwtn_root = (derMaxJ + 1) / 4. - 0.25;
	GetRootNewton(nwtn_root);
	double timefourier = nwtn_root * sadc_clock_;
	double ampfourier = IntegrateDerivative(derMaxJ, derMax);

	// Detect overflows
	TreatOverflow(sample);

	result_.ampl = ampfourier;
	result_.time = timefourier;
	result_.ready = 1;
	if( writeTree_ )
		ManageTree( ampfourier, timefourier);

	return true;
}

///////////////////////////////////////////////////////////////////////////////
void DigitizerSADCF::FourierTransform(){
	int min = *std::min_element(csample_.begin(), csample_.end());
	if(hMag_.size() != csample_.size()) {
		hMag_.resize(csample_.size());
		hPhase_.resize(csample_.size());
	}

	//accessing TVirtualFFT like it is done in TH1
	TString opt("R2C");
	int ndim[3];
	ndim[0] = csample_.size();
	ndim[1] = 1;
	ndim[2] = 1;
	TVirtualFFT *fft = TVirtualFFT::FFT(1, ndim, opt.Data());
	for (int binx = 0; binx < ndim[0]; binx++)
		fft->SetPoint(binx, csample_[binx] - min);
	fft->Transform();
	double re, im, ph;
	for (int binx = 0; binx < ndim[0]; binx++) {
		fft->GetPointComplex(binx, re, im);
		hMag_[binx] = TMath::Sqrt(re*re + im*im);
		if (TMath::Abs(re) > 1e-13){
			ph = TMath::ATan(im/re);
			//find the correct quadrant
			if (re<0 && im<0)
				ph -= TMath::Pi();
			if (re<0 && im>=0)
				ph += TMath::Pi();
		} else {
			if (TMath::Abs(im) < 1e-13)
				ph = 0;
			else if (im>0)
				ph = TMath::Pi()*0.5;
			else
				ph = -TMath::Pi()*0.5;
		}
		hPhase_[binx] =ph;
	}
	//fix for memory leak in root_version <=5.34.23
	if(TVirtualFFT::GetCurrentTransform()) {
		delete TVirtualFFT::GetCurrentTransform();
		TVirtualFFT::SetTransform(NULL);
	}

	for (int k = 0; k < numberFourierComponents_; k++) { // derivative calculation
		double amp = hMag_[k] * 2./ sampleSize_;
		if (k == 0){
			amp /= 2.;
		}
		double phase = hPhase_[k];
		for (int i = 0; i < (sampleSize_ * 4); i++) {
			double factor = -k * TMath::TwoPi() / sampleSize_;
			hDer_[i] += factor * amp * std::sin(k * TMath::TwoPi() * i / (sampleSize_ * 4.) + phase);
		}
	} // end derivative filling
}

///////////////////////////////////////////////////////////////////////////////
void DigitizerSADCF::GetRootNewton(double& root) {
	double ref = 0.;
	for (unsigned int q = 0; q < newtonNIter_; q++){
		double der2 = 0.;
		double der3 = 0.;
		if (std::fabs(root - ref) >= newtonPrecission_ || q == 0){
			ref = root;
			for (int k = 1; k < numberFourierComponents_; k++) {
				double amp = hMag_[k] * 2. / sampleSize_;
				double phase = hPhase_[k];
				double factor = k * TMath::TwoPi() / sampleSize_;
				der2 += factor * factor * amp * std::cos(k * TMath::TwoPi() * root / sampleSize_ + phase);
				der3 += factor * factor * factor * amp * std::sin(k * TMath::TwoPi() * root / sampleSize_ + phase);
			}
			root += (der2 / der3);
		}
		if(q == 9 && (std::fabs(root - ref) >= newtonPrecission_ )) {
			std::cout<<"In CsDigitizer Newton Algorithm Does Not Converge"<<std::endl;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
double DigitizerSADCF::IntegrateDerivative(const int jmax, const double &dmax) {
	double ampFourier = 0.;
	for (int p = -integrationWindow_; p < integrationWindow_; p++) {
		double poly = -dmax / parabolaShape_ * (p * p - parabolaShape_) ;
		if((jmax + p >= (int)hDer_.size()) || jmax + p < 0)
			ampFourier += poly;
		else if (hDer_[jmax + p] <= poly) {
			ampFourier += hDer_[jmax + p];
		}
		else {
			ampFourier += 2 * poly - hDer_[jmax + p];
		}
	}
	ampFourier /= 4.;
	return ampFourier;
}

///////////////////////////////////////////////////////////////////////////////
void DigitizerSADCF::TreatOverflow(const std::vector<uint16> &sample) {
	//this actually doesn't do anything right now
	if( debug_ ) std::cout << " CsDigitizerSADC::FitMaxAdvanced overflow_amplitude_overflow_amplitude_overflow_amplitude_ " << overflow_amplitude_ << std::endl;

	unsigned int new_over = DetectOverflow(sample);

	// Treat overflows
	if( new_over  ) {
		if( debug_ ) {
			std::cout << " NOverflow =  " << overflow_samples_.size() << " Overflow? " << std::endl;
			PrintSample(sample);
			std::cout << std::endl << " Csample: ";
			for ( int is=0; is < (int)sample.size(); is++ ) {
				std::cout << " " << (int)(10.*csample_[is]);
			}
			std::cout << std::endl;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void DigitizerSADCF::ManageTree( const double &ampfourier, const double &timefourier ) {
	static TTree* tree;
	static std::vector<double> *sample_for_tree;
	static std::vector<double> *magnitudes;
	static std::vector<double> *phases;
	static int ECnumber;
	static int iicell;
	static int cx;
	static int cy;
	static unsigned int trm;
	static double tcsp;
	static int runnb;
	static int eventnb;
	static double timeF;
	static double ampF;
	static double oldTime;
	static double oldAmpl;

	static bool initialized = false;
	if(!initialized) {
		initialized = true;
		CsHistograms::SetCurrentPath("/Calorimeter");
		tree = new TTree("CCS", "CCS");
		tree->Branch("sample", &sample_for_tree, 16000, 1);
		tree->Branch("magnitudes", &magnitudes, 16000, 1);
		tree->Branch("phases", &phases, 16000, 1);
		tree->Branch("EC", &ECnumber, "ECnumber/I");
		tree->Branch("iic", &iicell, "iic/I");
		tree->Branch("cx", &cx, "cx/I");
		tree->Branch("cy", &cy, "cy/I");
		tree->Branch("trm", &trm, "trm/i");
		tree->Branch("tcsp", &tcsp, "tcsp/D");
		tree->Branch("runnb", &runnb, "runnb/I");
		tree->Branch("eventnb", &eventnb, "eventnb/I");
		tree->Branch("timefourier", &timeF, "timefourier/D");
		tree->Branch("ampfourier", &ampF, "ampfourier/D");
		tree->Branch("timeOld", &oldTime, "timeOld/D");
		tree->Branch("ampOld", &oldAmpl, "timeOldr/D");
		sample_for_tree->resize(csample_.size());
		magnitudes->resize(hMag_.size() / 2);
		phases->resize(hPhase_.size() / 2);
	}

	for(unsigned int i = 0; i < csample_.size(); ++i){
		(*sample_for_tree)[i] = csample_[i];
	}
	for(unsigned int i = 0; i < hMag_.size() / 2; ++i){
		(*magnitudes)[i] = hMag_[i];
		(*phases)[i] = hPhase_[i];
	}
	ECnumber = GetParentName() == "EC01P1__" ? 1 : (GetParentName() == "EC02P1__" ? 2 : 0);
	iicell = mycell_;
	cx = parent_->GetColumnOfCell(iicell);
	cy = parent_->GetRowOfCell(iicell);
	trm = CsEvent::Instance()->getTriggerMask() & 0xffff;
	tcsp = CsEvent::Instance()->getTCSPhaseTime() - TCS_T0; //Artificial shift as used in the other calo classes
	runnb = CsEvent::Instance()->getRunNumber();
	eventnb = CsEvent::Instance()->getEventNumberInRun();
	timeF = timefourier;
	ampF = ampfourier;
	oldAmpl = csample_[result_.max_pos];

	tree->Fill();
}
