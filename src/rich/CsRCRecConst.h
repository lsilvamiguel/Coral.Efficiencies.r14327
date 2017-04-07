
#ifndef  RECCONST_H
#define  RECCONST_H

/*!
   \file    CsRCRecConst.h
   \brief   Constants for reconstruction for CsRichOne.
   \author  Paolo Schiavon
   \version $Revision: 1.28 $
   \date    $Date: 2007/04/11 13:06:42 $
*/

  #include <string>


  class CsRCRecConst {


  public:


    static CsRCRecConst* Instance();

    inline double TwoPI() const { return TwoPI_; };
    inline double RadDeg() const { return RadDeg_; };

    inline int outBufferSize() const { return outBufferSize_; };

    inline int kOffHi() const { return kOffHi_; };
    inline float momMinAcc() const { return momMinAcc_; };
    inline float momMaxAcc() const { return momMaxAcc_; };

    inline float PHmipCut() const { return PHmipCut_; };

    inline float ddCluUsex() const { return ddCluUsex_; };
    inline float ddCluUsey() const { return ddCluUsey_; };
    inline float maxCluSize() const { return maxCluSize_; };

    double CFRefInd();
    inline double CFRefIndUV() { return CFRefInd(); };
    double CFRefIndVS();
    inline void setCFRefInd( double CFRefInd ) { CFRefInd_ = CFRefInd; };
    inline void setCFRefIndVS( double CFRefInd ) { CFRefIndVS_ = CFRefInd; };
    inline double MCCFRefInd() const { return MCCFRefInd_; };
    inline double MCCFRefIndUV() const { return MCCFRefIndUV_; };
    inline double MCCFRefIndVS() const { return MCCFRefIndVS_; };
    inline void setMCCFRefInd( double CFRefInd ) {
      MCCFRefInd_ = CFRefInd;
      MCCFRefIndSet_ = true;
      MCCFRefIndUV_ = CFRefInd;
      MCCFRefIndSetUV_ = true;
    };
    inline void setMCCFRefIndUV( double CFRefInd ) {
      MCCFRefIndUV_ = CFRefInd;
      MCCFRefIndSetUV_ = true;
    };
    inline void setMCCFRefIndVS( double CFRefInd ) {
      MCCFRefIndVS_ = CFRefInd;
      MCCFRefIndSetVS_ = true;
    };
    inline float qzRefInd() const { return qzRefInd_; };

    inline float partPathFr() const { return partPathFr_; };

    inline std::string ringDefMode() { return ringDefMode_; };
    inline std::string peakSrcMode() { return peakSrcMode_; };
    inline int *mcanWind() { return mcanWind_; };
    inline int nMore() const { return nMore_; };

    inline int mcaScan() const { return mcaScan_; };
    inline float xlScan() const { return xlScan_; };
    inline float xuScan() const { return xuScan_; };
    inline float binScan() const { return binScan_; };

    inline int nPhotMin() const { return nPhotMin_; };
    inline float sigBfCut() const { return sigBfCut_; };
    inline float *nSigmaBf() { return nSigmaBf_; };

    inline float nBinPhi() const { return nBinPhi_; };
    inline float phiBin() const { return phiBin_; };
    inline int binThr() const { return binThr_; };
    inline float weight() const { return weight_; };

    inline float sig3SCut() const { return sig3SCut_; };

    inline float theAliMir() const { return theAliMir_; };
    inline float rrAliMir() const { return rrAliMir_; };
    inline float nPhoAliMir() const { return nPhoAliMir_; };
    inline float nSigAliMir() const { return nSigAliMir_; };

    inline double *massPartv() { return massPartv_; };
    inline int nPhotMinRing() const { return nPhotMinRing_; };
    inline float momMinProc() const { return momMinProc_; };
    inline float momMaxProc() const { return momMaxProc_; };
    inline float sigCut() const { return sigCut_; };
    inline double theMaxRgQ() const { return theMaxRgQ_; };

    inline float *nPhotAtSat() { return nPhotAtSat_; };
    inline bool bPhotAtSat() { return bPhotAtSat_; };
    inline float nZero() const { return nZero_; };
    inline float nZeroUV() const { return nZeroUV_; };
    inline float nZeroVS() const { return nZeroVS_; };
    inline std::string thetaType() const { return thetaType_; };
    inline std::string likeType() const { return likeType_; };
    inline std::string backgrType() const { return backgrType_; };
    inline void setBackgrType( std::string str ) { backgrType_ = str; };
    inline bool likeRatio() const { return likeRatio_; };
    inline bool likeLog() const { return likeLog_; };
    inline double likeDefVa() const { return likeDefVa_; };
    //map OBSOLETE!
    inline double backgrmapval(int i, int j) { return  backgrmap_[i][j]; };
    inline float *cathN0() { return cathN0_; };

    inline int nEveDSkip() const { return nEveDSkip_; };
    inline int nEveDisplay() const { return nEveDisplay_; };
    inline std::string dispMode() const { return dispMode_; };

    //inline float *cathodeThre() { return cathodeThre_; };
    float* cathodeThre();
    inline float *cathPadThre() { return cathPadThre_; };

    inline bool printConsts() const { return printConsts_; };

    void print();


  protected:


    CsRCRecConst();

    ~CsRCRecConst();


  private:


    static CsRCRecConst* instance_;

    double TwoPI_;
    double RadDeg_;

    int outBufferSize_;

    int kOffHi_;

    float momMinAcc_;
    float momMaxAcc_;

    float PHmipCut_;
    float ddCluUsex_, ddCluUsey_;
    float maxCluSize_;

    double CFRefInd_;
    double CFRefIndUV_;
    double CFRefIndVS_;
    double MCCFRefInd_;
    double MCCFRefIndUV_;
    double MCCFRefIndVS_;
    bool MCCFRefIndSet_;
    bool MCCFRefIndSetUV_;
    bool MCCFRefIndSetVS_;
    float qzRefInd_;

    float partPathFr_;

    std::string ringDefMode_;
    std::string peakSrcMode_;
    int mcanWind_[2];
    int nMore_;

    int mcaScan_;
    float  xlScan_;
    float  xuScan_;
    float  binScan_;

    int nPhotMin_;
    float sigBfCut_;
    float nSigmaBf_[2];

    float nBinPhi_;
    float phiBin_;
    int binThr_;
    float weight_;

    float sig3SCut_;

    float theAliMir_;
    float rrAliMir_;
    float nPhoAliMir_;
    float nSigAliMir_;

    double massPartv_[31];
    int nPhotMinRing_;
    float momMinProc_;
    float momMaxProc_;
    float sigCut_;
    double theMaxRgQ_;

    float nPhotAtSat_[16];
    bool bPhotAtSat_;
    float nZero_;
    float nZeroUV_;
    float nZeroVS_;
    std::string thetaType_;
    std::string likeType_;
    std::string backgrType_;
    bool likeRatio_;
    bool likeLog_;
    double likeDefVa_;
    //map OBSOLETE!
    double backgrmap_[400][400];
    float cathN0_[16];

    int nEveDSkip_;
    int nEveDisplay_;
    std::string dispMode_;

    bool printConsts_;

    float cathPadThre_[20];
    float cathodeThre_[20];

  };

#endif
