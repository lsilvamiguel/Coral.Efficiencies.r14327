// $Id: CsTime.h,v 1.11 2010/01/24 16:06:26 suhl Exp $

/*!
  \file    CsTime.h
  \brief   Compass Time class
  \author  Benigno Gobbo
  \version $Revision: 1.11 $
  \date    $Date: 2010/01/24 16:06:26 $
*/  

#ifndef CsTime_h 
#define CsTime_h 

#include "coral_config.h"
#include <ctime>
#include <iostream>


/*!
  \class   CsTime
  \brief   Compass Time class.

  Time Class. Created to make easy dates and times comparisons. Also the
  ostream \<\< operator was overloaded to type dates correctly. 
  Seconds from epoch (1 January 1970 00:00:00 UTC) and microseconds are 
  internally store.
*/

std::ostream& dmy(std::ostream&); //!< method used in \<\< 
std::ostream& sec(std::ostream&); //!< method used in \<\< 

class CsTime {

 public:
  
  /*! \fn CsTime();
    \brief Default constructor. 
    Time is set to the local machine time.
   */
  CsTime( );                                   

  /*! \fn CsTime( unsigned int secFrEpoch, int usec ); 
    \brief Constructor. 
    \param secFrEpoch Seconds from 1 Jan. 1970 00:00:00 UTC
    \param usec       micro Seconds
   */
  CsTime( time_t secFrEpoch, int usec ); 

  /*! \fn CsTime(int year, int month, int day, 
                 int hour=0, int min=0, int sec=0, int usec=0); 
    \brief Constructor. 
    \param year  year number  (four digits, e.g. 1962)
    \param month month number (1-12)
    \param day   day (1-31)
    \param hour  hour (0-23), defaults to 0
    \param min   minute (0-59), defaults to 0
    \param sec   second (0-59), defaults to 0
    \param usec  microseconds (0,99999), defaults to 0
   */
  CsTime( int year, int month, int day, 
	  int hour = 0, int min = 0, 
	  int sec = 0, int usec = 0 ); 

  CsTime( const CsTime & );                //!< Copy Constructor
  CsTime& operator=( const CsTime & );     //!< Assigment operator
  
  /*! \fn inline unsigned int secFrEpoch();
    \brief Returns the number of seconds from 1 Jan. 1970 00:00:00 UTC
  */
  inline time_t secFrEpoch() { return( _secFrEpoch ); }

  /*! \fn inline int          usec();
    \brief Returns the number of micro-seconds 
  */
  inline int          usec() { return( _usec ); }

  bool operator==( const CsTime & ) const; //!< "equal to" operator
  bool operator!=( const CsTime & ) const; //!< "not equal to" operator
  bool operator>( const CsTime & ) const;  //!< "greather than" comparison operator
  bool operator<( const CsTime & ) const;  //!< "less than" comparison operator
  bool operator>=( const CsTime & ) const; //!< "greather or equal to" operator
  bool operator<=( const CsTime & ) const; //!< "less or equal to" operator

  /*! \fn  std::ostream & operator<<( ostream &, const CsTime & ); 
    \brief std::ostream << operator overload.
    \par usage sample:
    cout << dmy << myTime << endl; // dumps day of week, day, month, year, time
    <br>
    cout << sec << myTime << endl; // dumps seconds 
   */  
  friend std::ostream & operator<<( std::ostream &, const CsTime & ); 

  friend std::ostream& dmy(std::ostream&); //!< method used in \<\< 
  friend std::ostream& sec(std::ostream&); //!< method used in \<\< 

 private:
  
  static const int _usecInSec;  //!< number of microseconds in a second
  time_t           _secFrEpoch; //!< seconds from 1 Jan. 1970 00:00:00 UTC
  int              _usec;       //!< micro-seconds
  static bool       _dmy;       //!< print flag for \<\<

  /*! \fn inline static void _setDmy( const bool dmy );
    \brief sets the _dmy private attribute
    \param dmy the value to be set for dmy
  */
  inline static void _setDmy( const bool dmy ) { _dmy = dmy; }
};

#endif // CsTime_h 
