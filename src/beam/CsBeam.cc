// $Id: CsBeam.cc,v 1.6 2010/05/24 21:02:45 ybedfer Exp $

/*!
   \file CsBeam.cc
   \brief Compass Beam Track Class.
   \author  G. Khaustov
   \version $Revision: 1.6 $
   \date    $Date: 2010/05/24 21:02:45 $
*/
#include "CsBeam.h"
//#include "CsEvent.h"
//#ifdef COMPASS_USE_OSPACE_STD
//#  include <ospace/std/algorithm>
//#else
//# include <algorithm>
//#endif

CsBeam::~CsBeam()
{
}

CsBeam::CsBeam()  {
  BmTrackTime = 0.;
  NBMShits    = 0 ;      
  NtotBMShits = 0 ;       
  BMStime     = 0.;        
  BMStmchi2   = 0.;      
  BMSpchi2    = 0.;        
  NScFbhits   = 0 ;
  NtotScFbHits= 0 ;     
  ScFbTime    = 0.; 
  ScFbTimeChi2= 0.;   
  TrigLabel   =false;
  BackPropLH  = -1; // I.e. unphysical value
  Chi2CutFlag =false;
}


CsBeam::CsBeam(double bmtracktime, int nbmshits, int ntotbmshits, 
        double bmstime,  double bmstmchi2, double bmspchi2,
        int nscfbhits,   int  ntotscfbhits, double scfbtime, 
        double scfbtimechi2) :  BmTrackTime(bmtracktime),   
        NBMShits(nbmshits), NtotBMShits(ntotbmshits), BMStime(bmstime),     
        BMStmchi2(bmstmchi2),  BMSpchi2(bmspchi2), NScFbhits(nscfbhits),
        NtotScFbHits(ntotscfbhits), ScFbTime(scfbtime),
        ScFbTimeChi2(scfbtimechi2){}     

CsBeam::CsBeam( const CsTrack& track) : CsTrack(track),
	BmTrackTime(0),
	TrigLabel(0), Chi2CutFlag(0), BackPropLH(-1),
        NBMShits(0), NtotBMShits(0),
        BMStime(0),   BMStmchi2(0),  
        BMSpchi2(0),
        NScFbhits(0),NtotScFbHits(0),
        ScFbTime(0),  ScFbTimeChi2(0) {}

CsBeam::CsBeam( const CsBeam& track ) : CsTrack( track ),
        BmTrackTime(track.BmTrackTime),
        TrigLabel(track.TrigLabel),
	Chi2CutFlag(track.Chi2CutFlag), BackPropLH(track.BackPropLH),
        NBMShits(track.NBMShits), NtotBMShits(track.NtotBMShits),
        BMStime(track.BMStime),   BMStmchi2(track.BMStmchi2),  
        BMSpchi2(track.BMSpchi2),
        NScFbhits(track.NScFbhits),NtotScFbHits(track.NtotScFbHits),
        ScFbTime(track.ScFbTime),  ScFbTimeChi2(track.ScFbTimeChi2) {}

CsBeam& CsBeam::operator=(const CsBeam& track ) {
  if( this   != &track ) {
        BmTrackTime =  track.BmTrackTime;
        NBMShits    =  track.NBMShits;
        NtotBMShits =  track.NtotBMShits;
        BMStime     =  track.BMStime;     
        BMStmchi2   =  track.BMStmchi2;
        BMSpchi2    =  track.BMSpchi2;
        NScFbhits   =  track.NScFbhits;
        NtotScFbHits=  track.NtotScFbHits;
        ScFbTime    =  track.ScFbTime;
        ScFbTimeChi2=  track.ScFbTimeChi2;
        TrigLabel   =  track.TrigLabel;
	BackPropLH  =  track.BackPropLH;
	Chi2CutFlag =  track.Chi2CutFlag;
   (CsTrack)(*this) =  (CsTrack)track;
  }
  return( *this );
}

bool CsBeam::operator==( const CsBeam& track ) {
  if( 
      BmTrackTime  == track.BmTrackTime  &&
      NBMShits     == track.NBMShits     &&
      NtotBMShits  == track.NtotBMShits  &&
      BMStime      == track.BMStime      &&
      BMStmchi2    == track.BMStmchi2    &&
      BMSpchi2     == track.BMSpchi2     &&
      NScFbhits    == track.NScFbhits    &&
      NtotScFbHits == track.NtotScFbHits &&
      ScFbTime     == track.ScFbTime     &&
      ScFbTimeChi2 == track.ScFbTimeChi2 &&
      TrigLabel    == track.TrigLabel    &&
      BackPropLH   == track.BackPropLH   &&
      Chi2CutFlag  == track.Chi2CutFlag  &&
      (CsTrack)(*this) == (CsTrack)track )  return( true );
  else
      return( false );
}
void CsBeam:: setBmdata(const double bmtracktime, 
                        const int nbmshits,     const int ntotbmshits,    
                        const double bmstime,   const double bmstmchi2,
                        const double bmspchi2, 
                        const int nscfbhits,   const int ntotscfbhits,
                        const double scfbtime, const double scfbtimechi2)
{
        BmTrackTime =bmtracktime;
        NBMShits    =nbmshits;
        NtotBMShits =ntotbmshits;
        BMStime     =bmstime;
        BMStmchi2   =bmstmchi2;
        BMSpchi2    =bmspchi2;
        NScFbhits   =nscfbhits;
        NtotScFbHits=ntotscfbhits;
        ScFbTime    =scfbtime;
        ScFbTimeChi2=scfbtimechi2;
}

