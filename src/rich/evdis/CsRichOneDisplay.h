#ifndef  CsRichOneDisplay_h
#define  CsRichOneDisplay_h

/*!
   \file    CsRichOneDisplay.h
   \--------------------------
*/


#include <iostream>
#include <ostream>

class CsRCGraph;
class CsRCDisplay;      
class CsRCIFtoDisplay;

class CsRichOneDisplay {

  
 public:

  static CsRichOneDisplay* Instance();

  CsRichOneDisplay( const CsRichOneDisplay& );

 
  CsRCDisplay*     getDisplay() { return( display_ ); }
  inline CsRCIFtoDisplay* getIFtoDisplay() { return( IFdisplay_ ); }


  void doRichOneDisplay();

 protected:

  CsRichOneDisplay();

  ~CsRichOneDisplay();
  

 private:

  static CsRichOneDisplay* instance_;

  CsRCGraph*        graph_;  //!< Pointer to the CsRCGraph singleton
  CsRCDisplay*      display_;
  CsRCIFtoDisplay*  IFdisplay_;


  };

#endif
