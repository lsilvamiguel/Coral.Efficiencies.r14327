#ifndef __PlaneGandalf__
#define __PlaneGandalf__

#include "Plane.h"

#ifndef __CINT__
#include "ChipGandalf.h"
#endif
/// Plane for detectors with a readout based on an GANDALF chip.

class TH1F;
class TH2F;
struct cfd_info;
class PlaneGandalf : public  Plane 
{

 	public:

  	PlaneGandalf(const char* name);
  	virtual ~PlaneGandalf() {};

#ifndef __CINT__
	virtual void StoreDigit(CS::Chip::Digit* digit);

	virtual void EndEvent(const CS::DaqEvent &event);
	virtual cfd_info* calc_cfd(CS::ChipGandalf::DigitGADC* digitGADC,int start=0, int increment=1, float integral=1., float baseline = 0.);
	static float getBaseline(CS::ChipGandalf::DigitGADC* digitGADC,
			CS::uint32 sigma,
			CS::uint32 deltaBin,
			CS::uint32 startSample,
			CS::uint32 deltaSample);
#endif

	/// Resets all histograms (and associated counters)
	virtual void ResetHistograms();


	virtual void Init(TTree* tree =0);

	virtual void ControlPanel(const TGWindow* p, const TGWindow* main);

	int channel;

	bool isDebug;

  	/// histograms
  	
  	TH1F *fAmpHist; //!Amplitudes histo
  	TH1F *fIntHist; //!Integral histo
  	TH1F *fHRTHist; //!Hi Res Time histo
  	TH1F *fFrameTHist; //!Hi Res Time histo
  	TH1F *fFrameHist; //!Frame histogram 
  	TH1F *fFramePlotHist; //!'Plotting' histogram last event
  	TH1F *fCFDPlotHist; //!'Plotting' CFD histogram last event
  	TH1F *fFramePlotHistThresh; //!'Plotting' histogram last above thresh
  	TH1F *fCFDPlotHistThresh; //!'Plotting' CFD histogram last above thresh
  	TH1F *fFramePlotHistLast; //!'Plotting' histogram last event with hit
  	TH1F *fCFDPlotHistLast; //!'Plotting' CFD histogram last event with hit
  	TH2D *fTimeVSevent;
  	TH2D *fTimeDecVSevent;
  	TH2F *fOverlayedPulses;
  	TH2F *fOverlayedPulsesHS;
	TH1F *fHpEvtHist;
	TH1F *fHpEvtHist_computed;
	TH1F *integral_computed;
	TH1F *fAmpHist_computed;
	TH1F *fFrameTHist_computed;
	TH1F *fHRTHist_computed;

	TH1F *fDiffFrameTime;
	TH1F *fDiffHighRes;

	TH1F *fTimeDecoded; //!decoded time histogram
	TH1F *fTimeDecodedCut; //!decoded time histogram with fixed range
	TH1F *fBaseLine; //!diff computed vs gandalf Integral
	TH1F *fBaseLine_computed; //!diff computed vs gandalf Integral
	TH1F *fDiffIntegral; //!diff computed vs gandalf Integral
	TH1F *fDiffAmp; //!diff computed vs gandalf Amplitude
	TH1F *fDiffBaseline; //!diff computed vs gandalf Amplitude
	TH1F *fBitStat; //!statistics of occurences of bits in frame words
	TH1F *fBitStat_01; //!statistics of occurences of bits in frame words



	protected:


	//////////////// constants for the CFD
	////////////////////////////////////////////////////////////

	// weather to use custom constants for cfd or the ones gandalf transmits
	bool USE_CUSTOM_CONSTANTS;

	double 			FRACTFACT;
	int				DELAY;
	int				THRESHOLD;

	// sample points before CFD zero crossing
	int				CFD_OFFSET;
	// number of sambles per pulse
	int				PULSE_WIDTH;

	int 			BASELINE;

	double time_min,time_max;
	double timeDec_min,timeDec_max;
	int firstEvent;

#ifndef __CINT__
	/*! \brief round float to integer
	 * 	\author Tobias Baumann
	 * 	\ingroup RoundFunctions
	 *
	 * 	Rounds an float value and returns an integer value.
	 *
	 * 	example: 3.67823 -> 4
	 */
	static CS::int32 roundToInt(float x);

	/*! \brief floor float to integer
	 * 	\author Tobias Baumann
	 * 	\ingroup RoundFunctions
	 *
	 * 	Floors an float value and returns an integer value.
	 *
	 * 	example: 3.67823 -> 3
	 */
	static CS::int32 floorToInt(float x);

	/*! \brief ceil float to integer
	 * 	\author Tobias Baumann
	 * 	\ingroup RoundFunctions
	 *
	 * 	Ceils an float value and returns an integer value.
	 *
	 * 	example: 3.11 -> 4
	 */
	static CS::int32 ceilToInt(float x);

	/*! \brief round float to next step value
	 * 	\author Tobias Baumann
	 * 	\ingroup RoundFunctions
	 *
	 * 	Round a float value to an given step value.
	 *
	 * 	example:
	 * 	step = 3 => 4.67823 -> 6
	 */
	static float round(float x, float step);
#endif

	ClassDef(PlaneGandalf,1)
};

#endif











