// $Id: CsMiscDetector.cc,v 1.4 2003/06/08 22:39:46 zvyagin Exp $

/*!
   \file    CsMiscDetector.cc
   \brief   Compass Fiber Hodoscope like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.4 $
   \date    $Date: 2003/06/08 22:39:46 $
*/

#include "CsMiscDetector.h"
#include "CsDigit.h"
#include "CsGeom.h"
#include "CsOpt.h"
#include "CsInit.h"
#include "CsEvent.h"
#include "DaqDataDecoding/ChipF1.h"

using namespace std;
using CS::DetID;

CsMiscDetector::CsMiscDetector( const string &TBname ) 
  : CsDet(DetID("MISC",0), TBname) {

  _decodingDone = false;
  _decode       = false;
  _decodeCard   = false;

  string tag = ""; 
  string key = "make decoding";
  bool status = CsOpt::Instance()->getOpt( tag, key );
  if( status ) {
    list<string> options;
    list<string>::iterator Is;
    status = CsOpt::Instance()->getOpt( tag, key, options );
    if( status ) {
      for( Is=options.begin(); Is!=options.end(); Is++ ) {
	if( *Is == "Misc" || *Is == "all" ) {
	  _decodeCard = true;
	}
      }
    }      
    else {
      _decodeCard = true;
    }
  }
}

bool CsMiscDetector::operator==( const CsMiscDetector& det ) const {
  return GetTBName()==det.GetTBName();
}


void CsMiscDetector::DecodeChipDigit(const CS::Chip::Digit &digit) {

  const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  if( d==NULL ) {
    throw CS::Exception("CsMiscDetector::DecodeRawData(): Wrong digit type");
  }

  // address
  int    addr = d->GetChannel();
  // ***** TIME RELATIVE to TRIGGER *****
  // Since I (Y.B.) don't know a priori which precision have misc chips
  // ...try successively the 2 possibilities (probably a loss of time)
  double time = d->GetTimeDecoded();

 // Create CORAL digit
  myDigits_.push_back( new CsDigit(*this, addr, &time, 1) );
}



