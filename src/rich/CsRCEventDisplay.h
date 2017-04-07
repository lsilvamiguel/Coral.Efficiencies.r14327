  #ifndef  EVENTDISPLAY_H
  #define  EVENTDISPLAY_H

/*!
   \file    CsRCEventDisplay.h
   \--------------------------
   \brief   CsRCEventDisplay class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    23 November 1999,  rev. October 2000
*/



  class CsRCEventDisplay  {


      public:

  CsRCEventDisplay();

  CsRCEventDisplay( const CsRCEventDisplay& );

  void doEveDisplay();

  void sideDisplay( int );

  inline  bool flag() const { return flag_; };

  void print();

  ~CsRCEventDisplay();
  

      private:

  bool flag_;


  };

  #endif
