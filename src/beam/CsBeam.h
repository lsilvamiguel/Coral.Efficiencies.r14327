// $Id: CsBeam.h,v 1.7 2010/05/25 21:18:04 ybedfer Exp $

/*!
   \file CsBeam.h
   \brief Compass Beam Track Class.
   \author  G. Khaustov
   \version $Revision: 1.7 $
   \date    $Date: 2010/05/25 21:18:04 $
*/
#ifndef CsBeam_h
#define CsBeam_h
#include "CsTrack.h"

/*  class CsBeam

    Compass Beam Track Class.

    This class describes a beam track. It inherits one standard Compass CsTrack with only
    one Helix inside. This Helix contains the beam track parameters found using the beam 
    fiber hodoscopes and the BMS station + some extra parameters. The space track points 
    are measured in the hodoskope nealist to the target.
*/
class CsBeam : public CsTrack{

 public:
//
//   Constructors
//

  CsBeam();                              // Default Constructor
  virtual ~CsBeam();

  CsBeam(double bmtracktime, int nbmshits, int ntotbmshits, 
        double bmstime,  double bmstmchi2, double bmspchi2,
        int nscfbhits,   int  ntotscfbhits, double scfbtime, 
        double scfbtimechi2);
//
//
  CsBeam( const CsTrack& );              // Copy Constructor

  CsBeam( const CsBeam& );              // Copy Constructor

  CsBeam& operator=( const CsBeam& );   // Assign Operator

  virtual bool IsBeamTrack() const { return true; }
  
  bool operator==( const CsBeam& );     // "Equal to" Operator
//
// Set of accessors
//
  inline double getTrackTime(void)   const { return BmTrackTime; }
  inline bool   getTrigLabel(void)   const { return TrigLabel; }
  inline bool   getChi2CutFlag(void) const { return Chi2CutFlag; }
  inline double getBackPropLH()      const { return BackPropLH; }
//  inline int getNTrackPoints(void) const { return(NBMShits+NScFbhits); }  //do not valib anymore
//
  inline double getBMStime(void)   const { return(BMStime); }
  inline double getBMStmchi2(void) const { return(BMStmchi2); }
  inline double getBMSpchi2(void)  const { return(BMSpchi2); }
  inline int getntotBMShits(void)  const { return(NtotBMShits); }
  inline int getnBMShits(void)     const { return(NBMShits); }
//
  inline int getntotScFbHits(void)    const { return(NtotScFbHits); }
  inline int getnScFbHits(void)       const { return(NScFbhits); }
  inline double getScFbTime(void)     const { return(ScFbTime); }
  inline double getScFbTimeChi2(void) const { return(ScFbTimeChi2); }
//
// Set of mutators
//
  inline void setTrackTime(const double tracktime) { BmTrackTime = tracktime;};
  inline void setTrigLabel(const bool triglabel)   { TrigLabel   = triglabel;};
  inline void setChi2CutFlag(const bool flag)      { Chi2CutFlag = flag;};
  inline void setBackPropLH(const double lh)       { BackPropLH  = lh;};
//
  inline  void setBMStime(const double bmstime){BMStime=bmstime;};
  inline  void setBMStmchi2(const double bmstmchi2){BMStmchi2=bmstmchi2;};
  inline  void setBMSpchi2(const double bmspchi2){BMSpchi2=bmspchi2;};
  inline  void setNtotBMShits(const int ntotbmshits){NtotBMShits=ntotbmshits;};
  inline  void setNBMShits(const int nbmshits){NBMShits=nbmshits;};
//
  inline  void setNtotScFbHits(const int ntotscfbhits){NtotScFbHits=ntotscfbhits;};
  inline  void setScFbHits(const int nscfbhits){NScFbhits=nscfbhits;};
  inline  void setScFbTime(const double scfbtime){ScFbTime=scfbtime;};
  inline  void setScFbTimeChi2(const double scfbtimechi2){ScFbTimeChi2=scfbtimechi2;};
//
void setBmdata(const double bmtracktime, 
               const int nbmshits,     const int ntotbmshits,    
               const double bmstime,   const double bmstmchi2,
               const double bmspchi2, 
               const int nscfbhits,    const int ntotscfbhits,
               const double scfbtime,  const double scfbtimechi2);

 private:
//
  double           BmTrackTime;    //  Track time with respect to the trigger time
  bool             TrigLabel;      //  Is track in a trigger time window?
  bool             Chi2CutFlag;    //  Is tracks back propagation P(chi^2) < 0.005
  double           BackPropLH;     // Back propagation P(chi^2), i.e. the likelihood of the association of a BMS track to a beamTelescope track, based on transport of the latter back to BMS through the optical system of the beam line.

//
//     BMS:
//
  int              NBMShits;       // Number of BMS hits used for this track
  int              NtotBMShits;    // The total number of fired BMS channels (in 4 hods)
  double           BMStime;        // BMS track time
  double           BMStmchi2;      // BMS track time Chi2 
  double           BMSpchi2;       // BMS momentum fit Chi2  
//
//    SciFb:
//
  int              NScFbhits;      // The number of the additional SciFb  hits in 3 sigma time window
  int              NtotScFbHits;   // The total number of the SciFb  hits (in 4 planes)
  double           ScFbTime;       // Beam Fiber Hodoscope track time
  double           ScFbTimeChi2;   // Beam Fiber Hodoscope track time  chiq

};

#endif // CsBeam_h
