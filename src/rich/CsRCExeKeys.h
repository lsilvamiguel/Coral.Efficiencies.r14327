#ifndef  EXECKEYS_H
#define  EXECKEYS_H

/*!
   \file    CsRCExeKeys.h
   \brief   Execution Keys for CsRichOne.
   \author  Paolo Schiavon
   \version 0.02,  rev. 20/6/00
   \date    June  2000
*/


  #include "CsOpt.h"


  class CsRCExeKeys {


  public:


  static CsRCExeKeys* Instance();

    void acknoMethod( std::string );

    inline bool MCarloEvent() const { return MCarloEvent_; }
    inline bool DataEvent() const { return DataEvent_; }

    inline bool readMyFile() const { return readMyFile_; }
    inline std::list<std::string> &RmyFileName() { return RmyFileName_; }
    void fillMyFileList( int, int );
    inline int myFileType() const { return myFileType_; }
    inline bool writeMyFile() const { return writeMyFile_; }
    inline std::string WmyFileName() const { return WmyFileName_; }
    inline bool endWriteMyFile() const { return endWriteMyFile_; }
    inline void setEndWrMyFile() { endWriteMyFile_ = true; }
    inline bool physSelection() const { return physSelection_; }
    inline int  nProMyEvs() const { return nProMyEvs_; }
    inline int  nSkipMyEvs() const { return nSkipMyEvs_; }

    inline bool PartFilt() const { return PartFilt_; }
    inline bool PartSele() const { return PartSele_; }
    inline bool ExitWind() const { return ExitWind_; }
    inline bool PartCorr() const { return PartCorr_; }
    inline bool PartReso() const { return PartReso_; }

    inline bool DoClu() const { return DoClu_; }

    inline bool UsePadThresh() const { return UsePadThresh_; }
    inline void setUsePadThresh( bool bb ) { UsePadThresh_ = bb; }
    inline bool KillHaloPads() const { return KillHaloPads_; }
    inline bool KillHaloClus() const { return KillHaloClus_; }

    inline bool selPMTonly() const { return selPMTonly_; };
    inline bool selAPVonly() const { return selAPVonly_; };

    inline bool padScrambled() const { return padScrambled_; }

    //inline bool NoRingSplitting() const { return NoRingSplitting_; }
    inline bool UseSplitRings() const { return UseSplitRings_; }
    inline bool CorrQzW() const { return CorrQzW_; }
    inline bool CorrPMTOpt() const { return CorrPMTOpt_; }

    inline bool likeONLY() const { return likeONLY_; };
    inline bool doCheckLike() const { return doCheckLike_; };
    inline bool doThetaLikeMax() const { return doThetaLikeMax_; };
    inline bool thetaLikeMaxTrue() const { return thetaLikeMaxTrue_; };

    inline bool SearchForRings() const { return SearchForRings_; }

    inline bool UseCluPH() const { return UseCluPH_; }
    inline bool RejWroPeak() const { return RejWroPeak_; }
    inline bool DoBkgFil() const { return DoBkgFil_; }
    inline bool DoWgtAv() const { return DoWgtAv_; }
    inline bool SigmaCut() const { return SigmaCut_; }

    inline bool CorrMirror() const { return CorrMirror_; }
    inline bool AliMirror() const { return AliMirror_; }

    inline bool EventDisp() const { return EventDisp_; }
    inline bool EventROOTDisp() const { return EventROOTDisp_; }

    inline bool EventAnasy() const { return EventAnasy_; }
    inline bool PartPhoSelect() const { return PartPhoSelect_; }
    inline bool RingSelect() const { return RingSelect_; }
    inline bool MCMoni() const { return MCMoni_; }
    inline bool DataMoni() const { return DataMoni_; }
    inline bool PartIdent() const { return PartIdent_; }

    inline int AcknoMethods() const { return AcknoMethods_; }
    inline int kPrintRichRec() const { return kPrintRichRec_; }
    inline int kPrintEventParts() const { return kPrintEventParts_; }
    inline int kPrintEventSimul() const { return kPrintEventSimul_; }
    inline int kPrintEventPads() const { return kPrintEventPads_; }
    inline int kPrintEventClus() const { return kPrintEventClus_; }
    inline int kPrintPartPhotons() const { return kPrintPartPhotons_; }
    inline int kPrintEventRings() const { return kPrintEventRings_; }
    inline int kPrintEventDisplay() const { return kPrintEventDisplay_; }
    inline int kPrintEventAnalysis() const { return kPrintEventAnalysis_; }
    inline int kPrintPartProbs() const { return kPrintPartProbs_; }
    inline int kPrintRejections() const { return kPrintRejections_; }
    inline int printFrequency() const { return printFrequency_; }

    inline void setPrintRichRec( int k ) { kPrintRichRec_ = k; }
    inline void setPrintEventParts( int k ) { kPrintEventParts_ = k; }
    inline void setPrintEventSimul( int k ) { kPrintEventSimul_ = k; }
    inline void setPrintEventPads( int k ) { kPrintEventPads_ = k; }
    inline void setPrintEventClus( int k ) { kPrintEventClus_ = k; }
    inline void setPrintPartPhotons( int k ) { kPrintPartPhotons_ = k; }
    inline void setPrintEventRings( int k ) { kPrintEventRings_ = k; }
    inline void setPrintEventDisplay( int k ) { kPrintEventDisplay_ = k; }
    inline void setPrintEventAnalysis( int k ) { kPrintEventAnalysis_ = k; }
    inline void setPrintPartProbs( int k ) { kPrintPartProbs_ = k; }
    inline void setPrintRejections( int k ) { kPrintRejections_ = k; }

    inline bool writeNtup() const { return writeNtup_; }

    inline bool readGf5() const { return readGf5_; }

    void print();


  protected:


    CsRCExeKeys();

    ~CsRCExeKeys();


  private:


    static CsRCExeKeys* instance_;


    bool MCarloEvent_;
    bool DataEvent_;

    bool readMyFile_;
    std::list<std::string> RmyFileName_;
    int myFileType_;
    bool writeMyFile_;
    std::string WmyFileName_;
    bool endWriteMyFile_;
    bool physSelection_;
    int  nProMyEvs_;
    int  nSkipMyEvs_;

    bool PartFilt_;
    bool PartSele_;
    bool ExitWind_;
    bool PartCorr_;
    bool PartReso_;

    bool DoClu_;

    bool UsePadThresh_;
    bool KillHaloPads_;
    bool KillHaloClus_;

    bool selPMTonly_;
    bool selAPVonly_;

    bool padScrambled_;

    //bool NoRingSplitting_;
    bool UseSplitRings_;
    bool CorrQzW_;
    bool CorrPMTOpt_;

    bool likeONLY_;
    bool doCheckLike_;
    bool doThetaLikeMax_;
    bool thetaLikeMaxTrue_;

    bool SearchForRings_;

    bool UseCluPH_;
    bool RejWroPeak_;
    bool DoBkgFil_;
    bool DoWgtAv_;
    bool SigmaCut_;

    bool CorrMirror_;
    bool AliMirror_;

    bool EventDisp_;
    bool EventROOTDisp_;

    bool EventAnasy_;
    bool PartPhoSelect_;
    bool RingSelect_;
    bool MCMoni_;
    bool DataMoni_;
    bool PartIdent_;

    int AcknoMethods_;
    int kPrintRichRec_;
    int kPrintEventParts_;
    int kPrintEventSimul_;
    int kPrintEventPads_;
    int kPrintEventClus_;
    int kPrintPartPhotons_;
    int kPrintEventRings_;
    int kPrintEventDisplay_;
    int kPrintEventAnalysis_;
    int kPrintPartProbs_;
    int kPrintRejections_;

    int printFrequency_;

    bool writeNtup_;

    bool readGf5_;

  };

#endif
