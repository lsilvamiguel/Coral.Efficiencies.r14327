// $Id: MagnetInfo.h,v 1.7 2006/06/16 15:21:47 conrad Exp $

/*!
   \file    MagnetInfo.h
   \brief   store detector.dat magnet informations
   \author  Hugo Pereira
   \version $Revision: 1.7 $
   \date    $Date: 2006/06/16 15:21:47 $
*/

#ifndef MagnetInfo_h
#define MagnetInfo_h

#include <TROOT.h>
#include <TObject.h>
#include <string>
class Obj3D;

/*! \class MagnetSize
   \brief to store magnets inner/outer half dimensions
*/
class MagnetSize {
  public:
  //! constructor
  MagnetSize( double xsizIn, double ysizIn, 
    double xsizOut, double ysizOut, double zsiz ): 
    xsizIn_( xsizIn ), 
    ysizIn_( ysizIn ), 
    xsizOut_( xsizOut ), 
    ysizOut_( ysizOut ), 
    zsiz_( zsiz ) {}
  double xsizIn_;  //!< inner magnet half width  [mm] (along x)
  double ysizIn_;  //!< inner magnet half height [mm] (along y)
  double xsizOut_; //!< outer magnet half width  [mm] (along x)
  double ysizOut_; //!< outer magnet half height [mm] (along y)
  double zsiz_; //!< magnet half depth  [mm] (along z)
};

/*!
   \class   MagnetInfo
   \brief   store detector.dat magnet informations
*/
class MagnetInfo: public TObject {
  public:
    
  /*! \fn   MagnetInfo( 
    string detFileName = "", 
    unsigned int lNumber = 0, 
    string detInfoLine = ""
    );
    \brief constructor
    \param detFileName the name of the detector file from which the Dead zone is read
    \param lNumber the corresponding line number in the file. It is kept for rewriting
    \param detInfoLine the line corresponding to the detector, used to parse the parameters
  */
  MagnetInfo( 
    std::string detFileName = "", 
    unsigned int lNumber = 0, 
    std::string detInfoLine = ""
  );
  
  #ifndef __CINT__
  //! return 3DShape object corresponding to active area
  Obj3D* GetArea3D( void );
  #endif
  int lineNumber_; //!< line number in det.dat file from which the magnets info have been read
  
  int id_;      //!< magnet id
  double xcm_;  //!< center along x [mm]
  double ycm_;  //!< center along y [mm]
  double zcm_;  //!< center along z [mm]
  
  int rot_;         //!< ?
  double scale_;    //!< magnetic field scale
  double flag1_;    //!< ?
  int flag2_;       //!< ?
  double current_;  //!< magnet current
  bool gotCurrent_; //!< true if magnet current was read from det file
 
  private:
  
  Obj3D* area_; //!< area box description

  static const MagnetSize SM1Size_;  //!< SM1 magnet dimensions
  static const MagnetSize SM2Size_;  //!< Sm2 magnet dimensions
  ClassDef(MagnetInfo,0)
  
};

#endif
  
  
