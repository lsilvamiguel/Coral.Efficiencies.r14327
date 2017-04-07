  #ifndef  HISTOINT_H
  #define  HISTOINT_H

/*!
   \file    HistoInt.h
   \------------------
   \brief   Interface to HBOOK.
   \author  Paolo Schiavon
   \version 1.0
   \date    29 November 1999
*/


  extern "C" {
  void histobo1_( const int&, const int&, const float&, const float&,
                  const float& );
  void histobo2_( const int&, const int&, const float&, const float&, 
                  const int&, const float&, const float&, const float& );
  void histof1_( const int&, const float&, const float& );
  void histof2_( const int&, const float&, const float&, const float& );
  void histofill_( const int&, const float&, const float&, const float& );
  void histodel_( const int& );
  }


  #endif
