// $Id: DetFileInfo.h,v 1.9 2002/09/16 17:23:30 hpereira Exp $

#ifndef DetFileInfo_h
#define DetfileInfo_h
 
#include <string>
#include <vector.h>
#include <TROOT.h>
#include <TObject.h>
#include <TMatrix.h>
#include <iostream.h>

class DetFileInfo: public TObject {
  public:
  DetFileInfo( 
    string detFileName = "", 
    unsigned int lNumber = 0, 
    string detInfoLine = ""
  );
  unsigned int Add( DetFileInfo* dfi );
  void DetFileInfo::UpdateRotMatrices( void );
  void DUtoCenter_CM( double du );
  DetFileInfo* GetMain( void );

  string Dump( void );
  
  //=== line number in detector.dat
  unsigned int lineNumber_;
  string Comment_;
  
  int id_;
  string TBName_; 
  string name_;
  int unit_;
  int type_;       //!< detector type
  double rdlen_;   //!< radiation length

  //=== detector size (mm)
  double zsiz_;
  double xsiz_;
  double ysiz_;
  
  //=== detector center position (mm)
  double xcm_;     //!< x center position from detector.dat file   (mm)
  double ycm_;     //!< y center position from detector.dat file   (mm)
  double zcm_;     //!< z center position from detector.dat file   (mm)

  int rotMNum_;    //!< Rotation MRS to DRS
  double wirD_;    //!< first wire position from detector.dat file (mm)
  double ang_;     //!< wire orientation in DRS (deg)

  int nWir_;       //!< number of wires
  double wirP_;    //!< pitch from detector.dat file               (mm)
  TMatrix rotM_;   //!< rotation matrix given by detector.dat file
	TMatrix irotM_;  //!< inverse rotation matrix given by detector.dat file

  double eff_;     //!< efficiency (for MC)
  double bkg_;     //!< background (for MC)
  
  double tGate_;   //!< time gate
  
  //=== more keys for drift like detectors
  double vel_;     //!< velocity mm/F1unit  
  double t0_;      //!< time offset t0 from detector.dat file      (t0)
  double thRes_;   //!< ???
  double spSli_;   //!< resolution (in F1 unit)
  double tiSli_;   //!< time slice (ns)
  
  double res_;     //!< resolution given by detector.dat file      (mm)
  bool isMWPC_;    //!< to handle detectors as MWPCs

  //=== SubDetectors contains at least _this_  
  vector< DetFileInfo* > sub_;
  
  private:
	vector< TMatrix > DRSrotM_;

  ClassDef(DetFileInfo,0)
  
};
#endif

