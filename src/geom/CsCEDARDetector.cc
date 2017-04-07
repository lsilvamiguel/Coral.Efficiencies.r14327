//============================================================================
// Name        : CsCEDARDetector.cc
// Author      : Prometeusz Jasinski Promme@web.de
// Version     : 2.0.3 created: 11/06/08 modified 09/23/09
// Copyright   : No copyrights
// Description : Implementation of CsCEDARDetector decoding library class
// Comments    :
//    (11/06/08)
//    - Building from scratch based on CsMiscDetector, CsCEDAR
//      (thanks to Vladimir Kolosov for this support), and the base class
//		CsDet
//	  (11/11/08)
//    - the base class is now the derived CsDetector class
//    - testing the embedded cvs environment with eclipse... seems to work
//      well if this comment will be submitted
//    (11/19/08)
//    - testing the embedded cvs environment by eclipse on windows platform
//      if you see this comment, the submit worked :)
//    (11/20/08)
//    - corrected sadc name (mistake by having 9 characters instead of 8)
//    (09/23/09)
//    - with the kind help of tobi, alex and sergei I was able to run this code
//      and solved some small bugs
//    - still not fully validated
//    - after discussion with sergei, information is not stored as a cluster
//      anymore
//    - SADC decoding is turned off. Caused sometimes exeptions.
//    (15/01/08(
//    -complete reimplemetation of the code in order to perform real cedar reconstruction
//
//============================================================================
#include <math.h>
#include <bitset>

#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"

#include "CsCEDARDetector.h"
#include "CsEvent.h"
#include "CsCluster.h"
#include "CsHistograms.h"
#include "CsOpt.h"

CsCEDARDetector::CsCEDARDetector( const std::string &TBname, const std::string &geom_file) :
                 CsDet(CS::DetID(TBname),TBname)
{
  if( strncmp( TBname.c_str(), "CE01", 4 ) == 0 ) {
    cedar_ = 1;
    sadc_name_ = "CE01P1_s";
  }
  else if( strncmp( TBname.c_str(), "CE02", 4 ) == 0 ) {
    cedar_ = 2;
    sadc_name_ = "CE02P1_s";
  }
  else {
    throw CS::Exception("CsCEDARDetector::CsCEDARDetector(): TBname does not match!");
  }
  for(unsigned int i = 0; i < 8; ++i) {
    t0_time_[i] = 0.;
    t0_rms_[i] = 0.;
    CLIM_[i] = NULL;
    remap_[i] = i;
  }
  nSigma_ = 3;
  hitMask_ = 0;
  offsetPhi_ = 9.5;

  storeHistos_ = CsOpt::Instance()->getOpt( "", "store cedar histogramms");

  if(storeHistos_) {
    CsHistograms::SetCurrentPath("/CEDAR");
    hDx_ = new TH1D(Form("h%sdx", TBname.c_str()), Form("h%sdx", TBname.c_str()), 500, -4e-4, 4e-4);
    hDy_ = new TH1D(Form("h%sdy", TBname.c_str()), Form("h%sdy", TBname.c_str()), 500, -4e-4, 4e-4);
    hTime_ = new TH2D(Form("h%sTime", TBname.c_str()), Form("h%sTime", TBname.c_str()), 8, -.5, 7.5, 200, -100, 100);
    hPhiNum_ = new TH2D(Form("h%sPhiNumber", TBname.c_str()), Form("h%sPhiNumber", TBname.c_str()), 600, 0, 2*TMath::TwoPi(),16, -8.5, 7.5);
  }

}

void CsCEDARDetector::DecodeChipDigit(const CS::Chip::Digit &digit){
  // by casting dynamically we determine weather we have an SADC
  // or F1 TDC entry to call
  const CS::ChipF1::Digit *_tdc_digit = dynamic_cast <const CS::ChipF1::Digit*> (&digit);
  if (_tdc_digit){
    DecodeChipF1Digit(*_tdc_digit);
    return;
  }

  const CS::ChipSADC::Digit *_sadc_digit = dynamic_cast <const CS::ChipSADC::Digit*> (&digit);
  if (_sadc_digit){
    DecodeChipSADCDigit(*_sadc_digit); // currently empty method
    return;
  }
  const CS::ChipGandalf::DigitGADC *_gandalf_digit = dynamic_cast <const CS::ChipGandalf::DigitGADC*> (&digit);
  if (_gandalf_digit){
    DecodeChipGandalfDigit(*_gandalf_digit); // currently empty method
    return;
  }

  // if none of them was correct, than it is an unexpected type of digit
  throw CS::Exception("CsCEDARDetector::DecodeChipDigit(): no TDC F1 nor SADC chip. decoding failed!");
}

void CsCEDARDetector::DecodeChipDigits(const CS::Chip::Digits &digits){
  typedef std::multimap<CS::DetID, CS::Chip::Digit*>::const_iterator T_iterator;
  // searching for all sadc digits
  // this part makes me think about 2008 where I originally wanted to name
  // the SADC and F1 data with the same detectorID
  // for 2008 it should be fine here
  std::pair<T_iterator,T_iterator> iterator_range = digits.equal_range(sadc_name_);
  for ( T_iterator digit_iterator = iterator_range.first; digit_iterator != iterator_range.second; digit_iterator++){
    DecodeChipDigit(*digit_iterator->second);
  }
  // this part will be called for F1 chips and probably in 2008 also for
  // SADC chips
  iterator_range = digits.equal_range(GetTBName());
  for ( T_iterator digit_iterator = iterator_range.first; digit_iterator != iterator_range.second; digit_iterator++){
    DecodeChipDigit(*digit_iterator->second);
  }
}

void CsCEDARDetector::DecodeChipF1Digit(const CS::ChipF1::Digit &digit){
  int pm   = digit.GetChannel();
  if ((pm < 0)||(pm > 7)){
    throw CS::Exception("CsCEDARDetector::DecodeChipF1Digit(): # pm is wrong. Please check mapping.");
    return;
  }
  // time in TDC channels
  raw_F1_time_.push_back(std::pair<int,double>(pm, digit.GetTimeDecoded()));
  //a PM is hit if the TDC time is arround T0
  if(storeHistos_)
    hTime_->Fill(pm, digit.GetTimeDecoded() - t0_time_[pm]);
  if(fabs(digit.GetTimeDecoded() - t0_time_[pm]) < nSigma_ * t0_rms_[pm])
    hitMask_ |= 1 << remap_[pm];
}

void CsCEDARDetector::DecodeChipGandalfDigit(const CS::ChipGandalf::DigitGADC &digit){
  return;
}

void CsCEDARDetector::DecodeChipSADCDigit(const CS::ChipSADC::Digit &digit) {
  return;
}

void CsCEDARDetector::getDXDY(const CsHelix* Hi, double& dxCEDAR, double& dyCEDAR) {
  double x30    = Hi->getX() / 10. + (transpMat_[0] - Hi->getZ() / 10.) * Hi->getDXDZ();
  double y30    = Hi->getY() / 10. + (transpMat_[0] - Hi->getZ() / 10.) * Hi->getDYDZ();
  double x_up   = transpMat_[53] * x30 + transpMat_[15] * 100 * Hi->getDXDZ();
  double y_up   = transpMat_[57] * y30 + transpMat_[19] * 100 * Hi->getDYDZ();
  double x_down = transpMat_[14] * x30 + transpMat_[15] * 100 * Hi->getDXDZ();
  double y_down = transpMat_[18] * y30 + transpMat_[19] * 100 * Hi->getDYDZ();
  double dz     = (-transpMat_[52] + transpMat_[13]);
  dxCEDAR = (x_up - x_down) / dz;
  dyCEDAR = (y_up - y_down) / dz;
}

void CsCEDARDetector::Clear() {
  raw_F1_time_.clear();
  hitMask_ = 0;
}

void CsCEDARDetector::readCalibration(time_t timePoint){
  CDB::Time tp(timePoint,0);
  tm *t = localtime(&tp.first);
  bool foundT0 = false;
  bool foundT0RMS = false;
  try {
    std::string strdata("");
    cdb_->read(GetTBName(), strdata, tp);
    if(strdata != "") {
      std::istringstream is(strdata);
      std::string str;
      while(std::getline(is, str)) {
        std::istringstream isLine(str);
        std::string type;
        isLine >> type;
        if(type == "timeT0") {
          for(unsigned int i = 0; i < 8; ++i)
            isLine >> t0_time_[i];
          foundT0 = true;
        }
        else if(type == "timeRMS") {
          for(unsigned int i = 0; i < 8; ++i)
            isLine >> t0_rms_[i];
          foundT0RMS = true;
        }
        else if(type == "xCent")
          isLine >> xCent_;
        else if(type == "yCent")
          isLine >> yCent_;
        else if(type == "offsetPhi")
          isLine >> offsetPhi_;
        else if(type == "transp") {
          double tmp;
          while(isLine >> tmp){
            transpMat_.push_back(tmp);
          }
        }
        else if(type == "remap") {
          for(unsigned int i = 0; i < 8; ++i)
            isLine >> remap_[i];
        }
        //read limit parametrization
        else if(type == "LIMIT") {
          std::string formula;
          isLine >> type >> formula;
          for(unsigned int i = 0; i < 8; ++i) {
            if(TString(type) == Form("C%dG%d%s%sLIM", cedar_,
                                                     (i < 4)     ? 0   : 1,
                                                     (i % 4 < 2) ? "P" : "K",
                                                     (i % 2)     ? "P" : "K")) {
              CLIM_[i] = new TF1(type.c_str(), formula.c_str());
            }
          }
        }
        //read likelihood parametrization
        else {
          for(unsigned int particle = 0; particle < 2; particle++) {
            for(unsigned int group = 0; group < 9; group++) {
              for(unsigned int hits = 0; hits < 3; hits++) {
                if( Form("C%dG%d%s%d", cedar_, group, particle == 0 ? "P" : "K" ,hits) == TString(type)) {
                  int nParams;
                  isLine >> nParams;
                  if(nParams != 13) {
                    throw CS::Exception("CsCEDARDetector wrong number of parameters");
                  }
                  for(int i = 0; i < 13; ++i) {
                    if(particle == 0)
                      isLine >> pilike_[group][hits][i];
                    else if(particle == 1)
                      isLine >> kaonlike_[group][hits][i];
                  }
                }
              }
            }
          }
        }
      //---------------------------------------
      }
    }
    if(!foundT0)
      throw CS::Exception("no T0 found");
    if(!foundT0RMS)
      throw CS::Exception("no T0 RMS found");
  }
  catch( const std::exception &e ) {
    std::cerr << "CsCEDARDetector::readCalibration() " << getName() << " error in reading for time point ";
    std::cerr << t << ": " << e.what() << std::endl;
  }
}


void CsCEDARDetector::getLikelihoods(std::vector<float> &cedarInfo, double &dxCEDAR, double &dyCEDAR) {
  double dxCmax = dxCEDAR - xCent_;
  double dyCmax = dyCEDAR - yCent_;
  if(storeHistos_ && hitMask_ == 0xFF) {
    hDx_->Fill(dxCmax);
    hDy_->Fill(dyCmax);
  }
  std::bitset <sizeof(hitMask_) * 8 > CEPMn(hitMask_);
  double CEDr = std::sqrt(dxCmax * dxCmax + dyCmax * dyCmax);
  double CEDphi = 0.5*TMath::Pi() + std::atan2(dxCmax,dyCmax);
  if (CEDphi < 0)
    CEDphi += 2 * TMath::Pi();
  if(storeHistos_ && CEDr < 1e-3) {
    for(int u = 0 ; u < 8; ++u){
      if(CEPMn[u]) {
        if((5 - 1.27 * CEDphi) < u)
        {
          hPhiNum_->Fill(CEDphi, u);
          hPhiNum_->Fill(CEDphi + TMath::TwoPi(), u - 8);
          if((14 - 1.27 * (CEDphi + TMath::TwoPi())) < u - 8)
            hPhiNum_->Fill(CEDphi, u - 8);
        }
        if((14 - 1.27 * (CEDphi+TMath::TwoPi())) > u)
        {
          hPhiNum_->Fill(CEDphi + TMath::TwoPi(), u);
        }
      }
    }
  }
  double CSlf = offsetPhi_ - 8. / TMath::TwoPi() * CEDphi; // result of linear fit
  if (CSlf > 7.75)
    CSlf -= 8.0;
  CSlf += 0.25;
  int CS = int(2 * CSlf);
  bool CisC = false;
  if (int(CS / 2.0) == CS / 2.0)
    CisC = true;
  int CG[9];
  for (int i=0; i < 9; i++)
    CG[i]=-1;
  int Ch = CS / 2;

  if (CisC) {
    CG[0] = CEPMn[Ch];
    CG[2] = CEPMn[Ch + 1 < 8 ? Ch + 1 : Ch - 7]
          + CEPMn[Ch - 1 < 0 ? Ch + 7 : Ch - 1];
    CG[4] = CEPMn[Ch + 2 < 8 ? Ch + 2 : Ch - 6]
          + CEPMn[Ch - 2 < 0 ? Ch + 6 : Ch - 2];
    CG[6] = CEPMn[Ch + 3 < 8 ? Ch + 3 : Ch - 5]
          + CEPMn[Ch - 3 < 0 ? Ch + 5 : Ch - 3];
    CG[8] = CEPMn[Ch + 4 < 8 ? Ch + 4 : Ch - 4];
  }
  else {
    CG[1] = CEPMn[Ch]
          + CEPMn[Ch + 1 < 8 ? Ch + 1 : Ch - 7];
    CG[3] = CEPMn[Ch - 1 < 0 ? Ch + 7 : Ch - 1]
          + CEPMn[Ch + 2 < 8 ? Ch + 2 : Ch - 6];
    CG[5] = CEPMn[Ch - 2 < 0 ? Ch + 6 : Ch - 2]
          + CEPMn[Ch + 3 < 8 ? Ch + 3 : Ch - 5];
    CG[7] = CEPMn[Ch - 3 < 0 ? Ch + 5 : Ch - 3]
          + CEPMn[Ch + 4 < 8 ? Ch + 4 : Ch - 4];
  }

  double CPLOG = 0;
  double CKLOG = 0;
  for (int ig = 0; ig < 9; ++ig) {
    if((CisC && ig % 2 == 1) || (!CisC && ig % 2 == 0))
      continue;
    int in = CG[ig];
    {
      double CGP = pilike_[ig][in][0];
      double CGK = kaonlike_[ig][in][0];
      for (int ia = 0; ia < 10; ia += 3) {
        CGP += pilike_[ig][in][1 + ia] *
          std::atan((CEDr - pilike_[ig][in][2 + ia]) / pilike_[ig][in][3 + ia]);
        CGK += kaonlike_[ig][in][1 + ia] *
          std::atan((CEDr - kaonlike_[ig][in][2+ia]) / kaonlike_[ig][in][3 + ia]);
      }
      double CGPL = CGP > 1e-8 ? log(CGP) : log(1e-8);
      double CGKL = CGK > 1e-8 ? log(CGK) : log(1e-8);
      CPLOG += CGPL;
      CKLOG += CGKL;
    }
  }

  double CPKLIM, CPPLIM, CKKLIM, CKPLIM;
  if(CisC) {
    CPKLIM = CLIM_[0]->Eval(CEDr);
    CPPLIM = CLIM_[1]->Eval(CEDr);
    CKKLIM = CLIM_[2]->Eval(CEDr);
    CKPLIM = CLIM_[3]->Eval(CEDr);
  }
  else {
    CPKLIM = CLIM_[4]->Eval(CEDr);
    CPPLIM = CLIM_[5]->Eval(CEDr);
    CKKLIM = CLIM_[6]->Eval(CEDr);
    CKPLIM = CLIM_[7]->Eval(CEDr);
  }
  cedarInfo.push_back(CPLOG - CPKLIM);
  cedarInfo.push_back(CPLOG - CPPLIM);
  cedarInfo.push_back(CKLOG - CKKLIM);
  cedarInfo.push_back(CKLOG - CKPLIM);
  double CPLOGID  = ((CPLOG < CPKLIM) - (CPLOG > CPPLIM));
  double CKLOGID  = ((CKLOG > CKKLIM) - (CKLOG < CKPLIM));
  double CPID = ((CKLOGID + CPLOGID) > 0) - ((CKLOGID + CPLOGID) < 0);
}
