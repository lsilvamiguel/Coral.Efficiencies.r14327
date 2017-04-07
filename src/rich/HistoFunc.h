  #ifndef  HISTOFUNC_H
  #define  HISTOFUNC_H

/*!
   \file    HistoFUNC.h
   \-------------------
   \brief   Interface to HBOOK.
   \author  Paolo Schiavon
   \version 1.0
   \date    December 1999
*/


  extern "C" {

  void hbook1_( const int&, const char*, const int&, const float&,
                const float&, const float& );
  void hbook2_( const int&, const char*, const int&, const float&, 
                const float&, const int&, const float&, const float&,
                const float& );
  void hf1_( const int&, const float&, const float& );
  void hf2_( const int&, const float&, const float&, const float& );
  void hfill_( const int&, const float&, const float&, const float& );

  float hsum_( const int& );

  void hropen_(const int&, const char*, const char*, const char*, const int&, 
               const int& );
  void hrput_( const int&, const char*, const char* );
  void hrout_( const int&, const int&, const char* );
  void hdelet_( const int& );

  }


  #endif
