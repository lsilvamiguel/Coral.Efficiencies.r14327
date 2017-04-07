// $Id: DeadZoneInfo.h,v 1.2 2007/06/01 09:13:51 conrad Exp $

/*!
   \file    DeadZoneInfo.h
   \brief   store detector.dat dead zone informations
   \author  Hugo Pereira
   \version $Revision: 1.2 $
   \date    $Date: 2007/06/01 09:13:51 $
*/

#ifndef DeadZoneInfo_h
#define DeadZoneInfo_h
 
#include <string>
#include <vector>

#include <TROOT.h>
#include <TObject.h>
#include <TMatrix.h>


class Obj3D;
class DetectorInfo;

/*! \class DeadZoneInfo 
    \brief store detector.dat dead zone informations
*/
class DeadZoneInfo: public TObject {
  public:
  /*! \fn   DeadZoneInfo( 
    string detFileName = "", 
    unsigned int lNumber = 0, 
    string detInfoLine = ""
    );
    \brief constructor
    \param detFileName the name of the detector file from which the Dead zone is read
    \param lNumber the corresponding line number in the file. It is kept for rewriting
    \param detInfoLine the line corresponding to the dead zone, used to parse the parameters
  */
  DeadZoneInfo( 
    std::string detFileName = "", 
    unsigned int lNumber = 0, 
    std::string detInfoLine = ""
  );
  
  /*! \fn string Dump( void ); 
     \brief returns a string containing all deadZone informations 
     in a detector.dat like format
  */
  std::string Dump( void );
  
  #ifndef __CINT__
  //! return 3DShape object corresponding to active area
  Obj3D* GetArea3D( void );
  #endif

  unsigned int lineNumber_;  //!< line number in det.dat file
  
  //! All members corresponding to the det.dat columns
  int id_;          //!< DeadZone unique Id                 
  std::string TBName_;   //!< associated detector TBName
  std::string name_;     //!< associated detector name
  int unit_;        //!< associated detector unit
  int shape_;       //!< shape index

  //! dead zone size (mm)
  double zsiz_;    //!< z size from det.dat (mm)
  double xsiz_;    //!< x size from det.dat (mm)
  double ysiz_;    //!< y size from det.dat (mm)
  
  //! dead zone center position (mm)
  double xcm_;     //!< x center position from detector.dat file   (mm)
  double ycm_;     //!< y center position from detector.dat file   (mm)
  double zcm_;     //!< z center position from detector.dat file   (mm)

  int rotMNum_;    //!< Rotation MRS to DRS
  TMatrix rotM_;   //!< rotation matrix given by detector.dat file
	TMatrix irotM_;  //!< inverse rotation matrix given by detector.dat file
 
  DetectorInfo* det_; //!< associated DetectorInfo Object if any, 0 otherwise
 
  private:
  
	std::vector< TMatrix > DRSrotM_; //!< DRS to MRS rotation matrix corresponding to rotMNum_
  Obj3D* area_;               //!< area description
  
  ClassDef(DeadZoneInfo,0)

};
#endif   
