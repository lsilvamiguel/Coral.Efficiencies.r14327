// $Id: DetectorInfo.h,v 1.17 2009/07/21 13:30:19 aaust Exp $

/*!
  \file    DetectorInfo.h
  \brief   store detector.dat detector informations
  \author  Hugo Pereira
  \version $Revision: 1.17 $
  \date    $Date: 2009/07/21 13:30:19 $
*/

#ifndef DetectorInfo_h
#define DetectorInfo_h
#include <string>
#include <vector>
#include <TROOT.h>
#include <TObject.h>
#include <TMatrix.h>
#include <iostream>

class Opt;
class Obj3D;
class DeadZoneInfo;

/*! \class DetectorInfo 
  \brief store detector.dat detector informations
*/

class DetectorInfo: public TObject {
public:
  
  /*! \fn   DetectorInfo( 
    string detFileName = "", 
    unsigned int lNumber = 0, 
    string detInfoLine = ""
    );
    \brief constructor
    \param detFileName the name of the detector file from which the Dead zone is read
    \param lNumber the corresponding line number in the file. It is kept for rewriting
    \param detInfoLine the line corresponding to the detector, used to parse the parameters
  */
  DetectorInfo( 
	       std::string detFileName = "", 
	       unsigned int lNumber = 0, 
	       std::string detInfoLine = ""
	       );
  
  /*! \fn unsigned int Add( DetectorInfo* dfi )
    \brief Add a subdetector to the current detectorInfo object. Returns the number of added subdetectors
    it is used for Micromegas and for straws.
    \param dfi pointer to the DetectorInfo to be added
  */
  unsigned int Add( DetectorInfo* dfi );
  
  //! Update detector rotations matrices according to DRS to MRS matrix and WRS to DRS wire angle
  void UpdateRotMatrices( void );
  
  /*! \fn  void DUtoCenter( double du )
    \brief update detector center position from an offset perp to the wire, DV for pixel detector
    \param du the offset [mm]
  */
  void DUtoCenter( double du );
  void DVtoCenter( double du );
  
  /*! \fn  int UtoWire( double u )
    \brief returns wire (if any) corresponding to u; -1 otherwise
    \param u coordinate perp to the wires [mm]
  */
  int UtoWire( double u );
  
  /*! \fn DetectorInfo* GetMain( void );
    \brief returns main detector in case of subdetectors. 
    For micromegas, this is the central region
    For straws, this is the b section
  */ 
  DetectorInfo* GetMain( void );
  
  /*! \fn string Dump( void ); 
    \brief returns a string containing all detector informations 
    in a detector.dat like format
  */
  using TObject::Dump; //jj to avoid compilation warning
  std::string Dump( void );

  /*! \fn void Print(); 
    \brief prints all detector informations 
    in a detector.dat like format
  */
  using TObject::Print; //jj to avoid compilation warning
  void Print() {std::cout<<Dump()<<std::endl;}

  /*! \fn SetBias( const char* biasL );
    \brief set detector bias for alignment debugging
    \param biasL the string containing all bias
    in the following order: "<BiasU(mm)> <BiasZ(mm)> <BiasT(deg)> <BiasP(noUnit)> <BiasR(mm)> <BiasL(noUnit)>"
  */
  void SetBias( const char* biasL );
  
#ifndef __CINT__
  //! return 3DShape object corresponding to active area
  Obj3D* GetArea3D( void );

  //! return 3DShape object corresponding wire iw, if any, 0 otherwise
  Obj3D* GetWire3D( int iw );
#endif
  
  //! line number in detector.dat
  unsigned int lineNumber_;
  
  //! All members corresponding to det.dat columns
  int id_;              //!< detector id
  std::string TBName_;  //!< detector TBName
  std::string name_;    //!< detector name
  int unit_;            //!< detector unit
  int type_;            //!< detector type
  double rdlen_;        //!< radiation length

  //! detector size (mm)
  double zsiz_;    //!< z size from det.dat (mm)
  double xsiz_;    //!< x size from det.dat (mm)
  double ysiz_;    //!< y size from det.dat (mm)
  
  //! detector center position (mm)
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
  double cosTheta_; //!< cosine of the WRS to MRS angle
  double sinTheta_; //!< sine of the WRS to MRS angle

  double eff_;     //!< efficiency (for MC)
  double bkg_;     //!< background (for MC)
  
  double tGate_;   //!< time gate
  
  //! drift like detectors members
  double vel_;     //!< velocity mm/F1unit  
  double t0_;      //!< time offset t0 from detector.dat file      (t0)
  double thRes_;   //!< ???
  double spSli_;   //!< resolution (in F1 unit)
  double tiSli_;   //!< time slice (ns)
  
  double res_;     //!< resolution given by detector.dat file      (mm)
  bool isMWPC_;    //!< to handle detectors as MWPCs
  
  int projection_; //!< projection id
  
  DeadZoneInfo* dead_; //!< associated DeadZoneInfo Object, if any, 0 otherwise.
  
  //=== Detector bias
  double biasU_;  //!< offset perp to the wires
  double biasV_;  //!< offset along "wires" (2nd projection for pixel detectors)
  double biasZ_;  //!< offset along the beam
  double biasT_;  //!< wire orientation offset
  double biasP_;  //!< pitch offset (pitch *= 1+biasP )
  double biasR_;  //!< R0 offset    (prop to T0 for drift like detectors)
  double biasL_;  //!< Lorentz angle offset

  std::vector< DetectorInfo* > sub_;  //!< vector of pointer to subdetectors. it contains at least this
   
private:
  std::vector< TMatrix > DRSrotM_;    //!< DRS to MRS rotation matrix corresponding to rotMNum_
  
  //=== 3D objects 
  Obj3D* area_;                //!< active area box description
  std::vector<Obj3D*> wires_;  //!< vector of wire box description

  //! Used to dump informations needed for the alignment
  friend std::ostream &operator << (std::ostream &o,const DetectorInfo &p);
  ClassDef(DetectorInfo,0)
  
    };
#endif

