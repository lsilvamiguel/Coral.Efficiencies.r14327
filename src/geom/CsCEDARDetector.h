//============================================================================
// Name        : CsCEDARDetector.h
// Author      : Prometeusz Jasinski Promme@web.de
// Version     : 2.0.0 created: 11/06/08 modified 11/11/08
// Copyright   : No copyrights
// Description : declaration of CsCEDARDetector decoding library class
// Comments    :
//    (11/06/08)
//    - Building from scratch based on CsMiscDetector, CsCEDAR
//      (thanks to Vladimir Kolosov for this support), and the base class
//		CsDet
//    (11/11/08)
//    - I decided to derive this class from the successor class CsDetector
//      since all containers and standard calls to lists of detectors in
//      coral need this class. So I can easily add it to the list of detectors
//    -	This is no tracking detector class, just PID class!
//
//============================================================================

#ifndef CsCEDARDetector_h
#define CsCEDARDetector_h

#include "CsDetector.h"
#include "coral_config.h"
#include "DaqDataDecoding/ChipF1.h"
#include "DaqDataDecoding/ChipGandalf.h"
#include "DaqDataDecoding/ChipSADC.h"
#include "CsHist.h"
#include "CsHistograms.h"
#include "CsHelix.h"
#include "CDB.h"

class CsDigit;
class TF1;
class TH1D;
class TH2D;

class CsCEDARDetector : public CsDet {
  public:
    CsCEDARDetector( const std::string &TBname, const std::string &geom_file);
    ~CsCEDARDetector() {}

    bool operator==( const CsCEDARDetector& det) const { return GetTBName()==det.GetTBName(); }
    static bool compById(CsCEDARDetector* lhs, CsCEDARDetector* rhs) { return lhs->getId() < rhs->getId(); }

    void makeMCDecoding(){}
    void clusterize(){}

    unsigned int getHitMask() const { return hitMask_; }
    int getId() const { return cedar_; }

    void DecodeChipDigits(const CS::Chip::Digits &digits);

    void readCalibration(time_t timePoint);

    const std::vector<double>& GetTransportMatrix() const { return transpMat_;}

    void Clear();

    void getDXDY(const CsHelix* Hi, double& dxCEDAR, double& dyCEDAR);
    void getLikelihoods(std::vector<float> &cedarInfo, double &dxCEDAR, double &dyCEDAR);


  private:
    // 1 or 2
    int cedar_;
    /*
     * Since we were running the SADC with a special detector ID
     * namely "CE0XP1_s" this has to be taken into account when
     * decoding the data
     */
    std::string sadc_name_;

    // raw time from the tdc
    std::vector < std::pair<int,double> > raw_F1_time_;

    /*
     * given by t0 calibration for each PM
     * Every pm gets it's own calibration constant since
     * they differ quite big
     */
    double t0_time_[8];
    double t0_rms_[8];
    double nSigma_;
    double xCent_;
    double yCent_;
    double offsetPhi_;
    std::vector<double> transpMat_;
    unsigned int hitMask_;
    int remap_[8];
    bool storeHistos_;
    TH1D* hDx_;
    TH1D* hDy_;
    TH2D* hTime_;
    TH2D* hPhiNum_;

    double pilike_[9][3][13];
    double kaonlike_[9][3][13];

    TF1 *CLIM_[8];

    void DecodeChipDigit(const CS::Chip::Digit &digit);
    void DecodeChipF1Digit(const CS::ChipF1::Digit &digit);
    void DecodeChipGandalfDigit(const CS::ChipGandalf::DigitGADC &digit);
    void DecodeChipSADCDigit(const CS::ChipSADC::Digit &digit);
};

#endif // CsCEDARDetector_h
