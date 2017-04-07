// $Id: CsMiscDetector.h,v 1.2 2003/04/09 15:49:30 benigno Exp $

/*!
   \file    CsMiscDetector.h
   \brief   Class to contain digits not associated to any detector.
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 2003/04/09 15:49:30 $
*/

#ifndef CsMiscDetector_h
#define CsMiscDetector_h

#include "coral_config.h"
#include "CsDet.h"
#include <list>
#include <set>
#include <string>

class CsDigit;


/*! \class CsMiscDetector 
    \brief Miscellanea Detectors Class.
*/

class CsMiscDetector : public CsDet {

 public:

  //! Constructor
  CsMiscDetector( const std::string &TBname );

  bool operator==( const CsMiscDetector& det ) const;

  //! Clear the list of detector digits.
  inline void clearDigitsList() { _decodingDone = false; myDigits_.clear(); } 

  //! Dummy method 
  void makeMCDecoding() {}

  //! Dummy method
  void clusterize() {}

  //! \c true if this detector data was arleady decoded
  bool decoded() const { return( _decodingDone ); }
  
  //! set this detector as decoded
  void setDecodingDone() { _decodingDone = true; }

  //! set decoding on this detector
  void setDecode() { _decode = true; }

  //! Decode raw data
  void DecodeChipDigit(const CS::Chip::Digit &digit);

 private:
  
  bool            _decodingDone;   //!< \c true id Detector decoded
  bool            _decodeCard;     //!< Decoding set by cards
  bool            _decode;         //!< Decoding set at run time

};

#endif // CsMiscDetector_h
