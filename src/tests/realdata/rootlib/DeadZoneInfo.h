// $Id: DeadZoneInfo.h,v 1.1 2002/07/08 20:34:22 hpereira Exp $

#ifndef DeadZoneInfo_h
#define DeadZoneInfo_h
 
#include <string>
#include <vector.h>
#include <TROOT.h>
#include <TObject.h>
#include <TMatrix.h>

class DeadZoneInfo: public TObject {
  public:
  DeadZoneInfo( 
    string detFileName = "", 
    unsigned int lNumber = 0, 
    string detInfoLine = ""
  );
  string Dump( void );

   
  //=== line number in detector.dat
  unsigned int lineNumber_;
  string Comment_;
  int id_;
  string TBName_; 
  string name_;
  int unit_;
  int shape_;

  //=== dead zone size (mm)
  double zsiz_;
  double xsiz_;
  double ysiz_;
  
  //=== dead zone center position (mm)
  double xcm_;     //!< x center position from detector.dat file   (mm)
  double ycm_;     //!< y center position from detector.dat file   (mm)
  double zcm_;     //!< z center position from detector.dat file   (mm)

  int rotMNum_;    //!< Rotation MRS to DRS
  TMatrix rotM_;   //!< rotation matrix given by detector.dat file
	TMatrix irotM_;  //!< inverse rotation matrix given by detector.dat file
  
  private:
	vector< TMatrix > DRSrotM_;

  ClassDef(DeadZoneInfo,0)

};
#endif   
