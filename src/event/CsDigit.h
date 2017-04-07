// $Id: CsDigit.h,v 1.2 2002/10/22 14:31:04 zvyagin Exp $

/*!
   \file    CsDigit.h 
   \brief   Compass Digit Class.
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 2002/10/22 14:31:04 $
*/

#ifndef CsDigit_h
#define CsDigit_h

#include "coral_config.h"

#ifdef COMPASS_USE_OSPACE
# include <ospace/std/vector>
#else
# include <vector>
#endif

#if USE_ObjectsCounter
#include "ObjectsCounter.h"
#endif

class CsDet;

/*! \class CsDigit 
    \brief Compass Digit Class.
*/


class CsDigit {

 public:

  //! Default constructor
  CsDigit(); 

  //! Constructor
  CsDigit( CsDet& det, int address ); 

  //! Constructor
  CsDigit( CsDet& det, int address, double* data, int datasize = 1 ); 

  //! Destructor
  virtual ~CsDigit();

  //! Copy Constructor
  CsDigit( const CsDigit& );               
  
  //! Assignment Operator
  CsDigit& operator=( const CsDigit& );    

  //! "equal to" Operator
  bool operator==( const CsDigit& ) const; 

  //! "less than" Operator 
  bool operator<( const CsDigit& ) const;  

  //! Returns the digit address 
  int getAddress() const { return( address_ ); }

  //! Returns 1st digit datum
  double getDatum() const;

  //! Returns the digit data
  double* getData() const { return( data_ ); }

  //! Returns the size of digit data
  int getDataSize() const { return( datasize_ ); }

  //! Replace 1st datum
  void replaceDatum( double datum );

  //! Returns the pointer to the associated CsDet object.
        CsDet* getDet()       { return( det_); }
  const CsDet* getDet() const { return( det_); }

  //! Returns \c true if this digit belongs to a cluster
  inline bool isClusterized() const { return( clusterized_ ); }

  // set this digit as clusterized
  inline void setClusterized() { clusterized_ = true; }

 private:

  int address_;         //!< digit address 
  int datasize_;        //!< digit data size
  double *data_;        //!< array of digit data
  CsDet* det_;          //!< pointer to detector associated to the digit
  bool clusterized_;    //!< \c true if already clusterized

  #if USE_ObjectsCounter
  ObjectsCounter<CsDigit> objects_counter;
  #endif

};

#endif // CsDigit_h
