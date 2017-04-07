/*
 * PlaneCamera.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gotzl
 */

#ifndef PLANECAMERA_H_
#define PLANECAMERA_H_




////////////////////////////////////////////////////////////////////////////////////////////
////////////////// remove this when using the final CAMERA detector instead of the prototype
//#define PROTO
//----------------------------------------------------------------------------------------//

#define FRAME_PLOTS	false


#include "Plane.h"
#include "PlanePanel.h"
#include "PlaneGandalf.h"

#ifndef __CINT__
#include "ChipGandalf.h"
#include "ChipSADC.h"
#endif

class PlaneCamera: public Plane {

public:
	PlaneCamera(const char *detname);//constructor
	virtual ~PlaneCamera(){}; //destructor

#ifndef __CINT__

  	enum Type //enumerate, begin with -1
  	{
  		UNKNOWN			=		-1,
  		Aup,
  		Adown,
  		Bup,
  		Bdown
  	};


  	virtual void Init (TTree* tree = 0);
  	virtual void StoreDigit (CS::Chip::Digit* digit);
  	virtual void EndEvent (const CS::DaqEvent &event);
  	virtual void Reset (void);
	Type getType() const;

	// TODO: when
	struct ChannelInfo {

	};

	// map channel number to the list of digits in the channel
	typedef vector<CS::ChipGandalf::DigitGADC*> digits;
	map<int, digits > chanDigits;
	map<int, digits > chanFrames;
	typedef vector<CS::ChipF1::Digit*> digitsf1;
	map<int, digitsf1 > chanf1Digits;
	typedef vector<CS::ChipSADC::Digit*> digitsSADC;
	map<int, digitsSADC > chansadcDigits;

#endif

  	int N_CHANNEL;

#ifndef __CINT__
	map<int, PlaneGandalf*> g_planes;
#endif

#if USE_DATABASE == 1
	/// method to read calibration data
	void ReadCalib(const tm& t);
#endif
	time_t calibTime;  // the timestamp for the calibrations to use

protected:

#if USE_DATABASE == 1
  // Camera calibration reading: 1 calibration variable per channel
  class ChannelCalib {
  public:
    int ch;
    float t0;
    int flag;
    ChannelCalib() : ch(0),t0(0),flag(0) {}
    ChannelCalib(const char *s) {
      if(3 != sscanf(s,"%d%f%d",&ch,&t0,&flag)) {
	throw CS::Exception("PlaneCamera::ChannelCalib : bad line \"%s\"",s);
	std::cerr<<"bad line, exception not caught !"<<std::endl;
      }
    }
  };
  
  /// calibration for every channel
  std::vector<ChannelCalib> calib_data;

  friend istream& operator>>(istream& in, PlaneCamera::ChannelCalib &c) {
    in>>c.ch;
    in>>c.t0;
    in>>c.flag;
    return in;
  }
#endif

	TH1F *h1NumberOfSentData;
	TH2F *h2BaselineVsChannel;
	TH2F *h2BaselineCompVsChannel;
	TH2F *h2FrameTimeVsChannel;
	TH2F *h2AmplitudeMaximumVsChannel;
	TH2F *h2IntegralVsChannel;
	TH2F *h2HighResVsChannel;
	TH2F *h2TimeVsChannel;
	TH2F *h2MultiVsChannel;
	TH2F *h2MultiDiffGanF1VsChannel;
	TH2F *h2TimeDiffGanF1VsChannel;
	TH2F *h2DataTransmission;
	TH2F *h2IntegralDiffGanSADCVsChannel;
	TH2F *h2MaxAmpDiffGanSADCVsChannel;
	TH2F *h2MaxAmpQuotGanSADCVsChannel;
	map<int, TH2F*> h2MaxAmpSADCVsGandalf,h2IntegralSADCVsGandalf,h2MaxAmpDiffQuotGanSADCVsGandalf,h2MaxAmpDiffQuotGanSADCVsSADC;

	ClassDef(PlaneCamera,1)
};

#endif /* PLANECAMERA_H_ */
