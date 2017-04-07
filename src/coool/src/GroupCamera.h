/*
 * GroupCamera.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gotzl
 */

#ifndef GROUPCAMERA_H_
#define GROUPCAMERA_H_

#include "Group.h"
#include "PlaneCamera.h"
#include "TMath.h"

// forward declaration of classes used to better organize data in this group
class Szinti;
class ASzinti;

class GroupCamera: public Group {

public:
	GroupCamera(const char* name);
	virtual ~GroupCamera() {};

    void                Init               (void);

#ifndef __CINT__
    void                EndEvent                (const CS::DaqEvent &event);
    // map the plane 'name' to its pointer
    map<int, const PlaneCamera*> planes;

    // maps for szinits; id==wire
    map<int,ASzinti*> aszintis;
    map<int,Szinti*> bszintis;

    const CS::DaqEvent *event;

    // add histo is protected in Group.h ...
    void addHistogram(TH1* hist) {
    	AddHistogram(hist);
    }

    const PlaneCamera* getPlane(int plane) {
    	return planes[plane];
    }

    bool isExpert() {
    	return fExpertHistos;
    }


  	// stores the important information coming from the PMTs
  	struct Hit {
  		Hit (int id,const CS::ChipGandalf::DigitGADC* hit_up,const CS::ChipGandalf::DigitGADC* hit_down, double dt_calib, double speedOfLight, double zOffset=0)
  			: id(id),hit_up(hit_up),hit_down(hit_down),dt_calib(dt_calib),speedOfLight(speedOfLight) {

			mean		= (hit_up->GetTimeDecoded() + hit_down->GetTimeDecoded())/2.;
		 	dt			= ((hit_up->GetTimeDecoded() - hit_down->GetTimeDecoded()) +dt_calib);
		 	zpos		= dt*(speedOfLight)/2. + zOffset;
		 	energy		= TMath::Sqrt(hit_up->getMaxAmplitude() * hit_down->getMaxAmplitude());
		 	energyInt	= TMath::Sqrt(hit_up->getIntegral() * hit_down->getIntegral());
  		}

  		int id;

  		const CS::ChipGandalf::DigitGADC *hit_up,*hit_down;

  		double dt_calib;
  		double speedOfLight;

  		// mean time (T_up + T_down)/2
  		double mean;
  		// time difference up/down
  		double dt;

  		double zpos;
  		// TODO: Energy loss sqrt (E_up * E_down)
  		double energy;
  		// Integral
  		double energyInt;

  	};


  	class ProtonEvent
  	{
  	public:
  		ProtonEvent (Hit *Ahit,Hit *Bhit,double tof,double beta) : Ahit(Ahit),Bhit(Bhit),tof(tof),beta(beta) {}
  		~ProtonEvent () {
  			//delete Ahit;
  			//delete Bhit;
  		}

  		Hit *Ahit,*Bhit;

  		double tof;
  		double beta;


  	};

	vector<ProtonEvent*> protonCandidates;

#endif

    int B_CHANNEL,A_CHANNEL;
    bool initialized;

    TH1F *chan_multi;
    TH2F *chan_ampl;
    TH2F *chan_int;
    TH2F *chan_time;
    TH2F *chan_ftime;
    TH2F *chan_hrtime;
    TH2F *chan_zpos;
    TH2F *chan_tmean;
    TH2F *chan_tof;
    TH1F *track_mult;


protected:

#if USE_DATABASE == 1
  /// method to read calibration data
  void ReadCalib(void);

  /// reads calibration data from the database
  template <class T>
  void ReadFromDataBase(T& v, const tm& t) {
    if( fDataBase==NULL )
      throw CS::Exception("GroupCamera::ReadFromDataBase():  data base is not opened.");
    try{
      struct tm tt(t);
      CDB::Time tp(mktime(&tt),0);
      std::string strdata("");
      fDataBase->read(fName,strdata,tp);
      if (strdata == "") {
        std::cerr << "GroupCamera::ReadFromDataBase() "<<GetName()<<" : no calibration file found"<<std::endl;
        return;
      }
      //std::cerr << strdata;
      std::istringstream istrdata(strdata.c_str());

      char temp[256];
      while ( istrdata.peek()=='#' )  // remove comment lines starting with #
        istrdata.getline (temp,256);

      istrdata >> v;
    }
    catch(CS::Exception& e) {
      std::cout<<"rethrowing"<<std::endl;
      throw;
    }
  }

  // Camera calibration reading: this class holds the Mean, Diff and TOF calibration constants
  class CameraCalib {
  public:
    int ch;
    float Amean, Bmean, Adiff, Bdiff, TOF_AiBi, TOF_AiBj;
    CameraCalib() : ch(0), Amean(0), Bmean(0), Adiff(0), Bdiff(0), TOF_AiBi(0), TOF_AiBj(0) {}
    CameraCalib(const char *s) {
      if(7 != sscanf(s,"%d%f%f%f%f%f%f", &ch, &Amean, &Bmean, &Adiff, &Bdiff, &TOF_AiBi, &TOF_AiBj)) {
	throw CS::Exception("GroupCamera::CameraCalib : bad line \"%s\"",s);
	std::cerr<<"bad line, exception not caught !"<<std::endl;
      }
    }
  };
  
  /// calibration constants for Mean, Diff and TOF
  std::vector<CameraCalib> calib_data;

  friend istream& operator>>(istream& in, GroupCamera::CameraCalib &c) {
    in>>c.ch;
    in>>c.Amean;
    in>>c.Bmean;
    in>>c.Adiff;
    in>>c.Bdiff;
    in>>c.TOF_AiBi;
    in>>c.TOF_AiBj;
    return in;
  }
#endif //USE_DATABASE

};

#endif /* GROUPCAMERA_H_ */
