// $Id: CsStopwatch.h,v 1.3 2007/02/05 10:18:54 gobbo Exp $

/*!
  \file    CsStopwatch.h
  \brief   Compass class for elapsed time measurements
  \author  Benigno Gobbo
  \version $Revision: 1.3 $
  \date    $Date: 2007/02/05 10:18:54 $
*/  

#ifndef CsStopwatch_h 
#define CsStopwatch_h 

#include "CsSTD.h"

/*!
  \class   CsStopwatch
  \brief   Compass class for elapsed time measurements.
*/

class CsStopwatch {

 public:
  
  CsStopwatch( );     //!< Default constructor: 10 chronometers 

  /* \fn CsStopwatch( int maxtim ); 
     \brief Constructor
     \param Number of chronometers
   */
  CsStopwatch( int maxtim ); 

  ~CsStopwatch( );   //!< Destructor                               

  CsStopwatch( const CsStopwatch& sw );                //!< Copy Constructor
  CsStopwatch& operator=( const CsStopwatch& sw );     //!< Assigment operator

  /* \fn int start();
     \brief Starts a chronometer and returns its id number
   */
  int start();

  /* \fn double inter( int chrono );
     \brief Returns the elapsed time without stopping the chronometer
     \param chrono Chronometer id number
   */
  double inter( int chrono );
  
  /* \fn double stop( int chrono );
     \brief Returns the elapsed time and stops the chronometer
     \param chrono Chronometer id number
   */
  double stop( int chrono );

 private:

  int* _startS;    // Start time (seconds).
  int* _startU;    // Start time (microseconds).
  bool* _used;      //
  int   _maxtim;  // Max number of stopwatches.
  int   _ntim;    // Number of used stopwatches.

};

#endif // CsStopwatch_h 
