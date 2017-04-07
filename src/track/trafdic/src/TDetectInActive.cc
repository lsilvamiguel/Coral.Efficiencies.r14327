// $Id: TDetectInActive.cc,v 1.9 2011/02/27 10:29:40 clamar Exp $

#include "TDetect.h"
// Include "TEv" to have "make TRAFDIC=1" overwrite traffic (this is useless
// now (as of v1.4), as a "TDetect.h" specific to lattice exists).
#include "TEv.h"

/*!
  Fastest possible check that argument (double y, double z) is w/in active area.
  (This is a clone of "InActive(float,float)" supra. It ensures compatibility
  w/ TraFFiC (as opposed to TraFDic) code.)
*/

bool TDetect::InActive(double y, double z) const 
{
  // ***** TRANSFORM y,z to DRS *****
  float y_shift = y - x[1];
  float z_shift = z - x[2];
  float ydrs = cframe*y_shift + sframe*z_shift;
  float zdrs = cframe*z_shift - sframe*y_shift;

  if (fabs(ydrs)>=hsizY || fabs(zdrs)>=hsizZ) return false;   // out of frame
#ifdef InActive_CHECK_URANGE
  // U range only matters for stereo planes in DC like detectors, where
  // sensitive area has a rhombus shape. (And then it still isn't required
  // to prevent a hit from being picked-up during track search, for hits
  // outside U range won't be available in any case (by definition of U range).
  // The check is still needed to determine the ``firing expectancy'',
  // of a given detector, for a given track. Which we do in 2 cases:
  // when determining the ``expected'' hit map and when determining the
  // ``firing efficiency'' of a given candidate track (quantity that is
  // checked in the track search to qualify the candidate).)
  //  Note: TDetect::range is = size in the case of "HI05", cf. "TSetup::Init".
  float uwrs = ca*y+sa*z-UEdge;            // "u" referred to edge of detector
  if (uwrs<0 || Range<uwrs) return false;  // Out of measurement range
#endif

  // ***** CHECK DEAD ZONE *****
  float ydrs_shift = ydrs - DZydrs;
  float zdrs_shift = zdrs - DZzdrs;

  switch( DZtype )
    {
    case RECTANGULAR :
      { 
	float ydzrs = DZCadrs*ydrs_shift + DZSadrs*zdrs_shift;
	if(fabs(ydzrs)>DZydim) return true;
	float zdzrs = DZCadrs*zdrs_shift - DZSadrs*ydrs_shift;
	return (fabs(zdzrs)>DZzdim);
      }
    case CIRCULAR    : 
      return (ydrs_shift*ydrs_shift + zdrs_shift*zdrs_shift > DZydim);
    case NO          :
    default :
      return true;      
    }  
  return true;
}


bool TDetect::InActive(double c[2], double cov[3]) const
{

  return(true);
}

bool TDetect::WinRange(float y, float z) const 
{
  // This method only check that (y,z) is w/in the range in measured coordinate
  // covered by the detector. This range can differ sizeably from the size for
  // a detector w/ wires at an angle w.r.t. the frame. 

  float uwrs = ca*y+sa*z-UEdge;            // "u" refered to edge of detector
  if (uwrs<0 || Range<uwrs) return false;  // Out of measurement range
  return true;
}

/*!
  Fastest possible check that argument (y,z) is w/in massive area, i.e. active area + central dead zone if the latter is massive (e.g. MM, as opposed to ST, DR. MA, MB or Hs). The method is used is TraFDic's version of "TTrack::FullKF", in place of "TTrack::InFrame". It allows in particular not to process 2006/7 tracks through DRs as if they were systematically traversing their thickness of Pb.
*/

bool TDetect::InMassive(double y, double z) const 
{
  // ***** TRANSFORM y,z to DRS *****
  float y_shift = y - x[1];
  float z_shift = z - x[2];
  float ydrs = cframe*y_shift + sframe*z_shift;
  float zdrs = cframe*z_shift - sframe*y_shift;

  if (fabs(ydrs)>=hsizY || fabs(zdrs)>=hsizZ) return false;   // out of frame

  // ***** CHECK EMPTY ZONE *****
  float ydrs_shift = ydrs - DZydrs;
  float zdrs_shift = zdrs - DZzdrs;

  switch( EZtype )
    {
    case RECTANGULAR :
      { 
	float ydzrs = DZCadrs*ydrs_shift + DZSadrs*zdrs_shift;
	if(fabs(ydzrs)>DZydim) return true;
	float zdzrs = DZCadrs*zdrs_shift - DZSadrs*ydrs_shift;
	return (fabs(zdzrs)>DZzdim);
      }
    case CIRCULAR    : 
      return (ydrs_shift*ydrs_shift + zdrs_shift*zdrs_shift > DZydim);
    case NO          :
    default :
      return true;      
    }  
  return true;
}
