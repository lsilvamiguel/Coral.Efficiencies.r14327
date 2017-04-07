// $Id: CsField.h,v 1.25 2010/06/26 01:08:54 ybedfer Exp $

/*!
   \file CsField.h
   \brief Compass Magnetic Filed Class.
   \author  Benigno Gobbo, Alexandre Korzenev
   \version $Revision: 1.25 $
   \date    $Date: 2010/06/26 01:08:54 $
*/

#ifndef CsField_h
#define CsField_h

#include "CsSTD.h"
#include "CsErrLog.h"
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Matrix/Matrix.h>

class CsField;

/*! \struct CsMagInfo
 *  \brief Magnet information structure.
 *
 *  It is filled from the Comgeant detectors data file output. 
 *  Please, don't change anything in this file !                 
 *  The flags should be absolutely identical to those used in COMGeant. 
 */
struct CsMagInfo {
  int mag;        //!< Magnet number
  double xcm;     //!< Centre position X coordinate (MRS) (mm)
  double ycm;     //!< Centre position Y coordinate (MRS) (mm)
  double zcm;     //!< Centre position Z coordinate (MRS) (mm)
  int    rot;     //!< Rotation ??? (see Comgeant)
  double fsc;     //!< Field Scale Factor (see Comgeant)
  /*!
    Possible values for target magnet type:
    \arg \e 1 (or 0) COMPASS solenoid.
    \arg \e 2 SMC solenoid.
    \arg \e 3 SMC dipole.
    \arg \e 4 COMPASS dipole.

    Possible values for gap size of SM1: 172, 132 and 82 cm.
    \arg \e 172 only one to have been used as of 2010.
    \arg \e 132 has been evaluated on MC D0: improves on mass resolution, but cuts into the acceptance, w/ net effect not convincing, even w/ SMC target magnet.
    \arg \e 82 was originally foreseen for the hadron program (?).
  */
  double flg1;    //!< Gap size for SM1 in cm or magnet type for target.
  /*!
  Possible values for SM1:
  \arg \e 0 old (calculated) map is used.
  \arg \e 1 new (measured) map is used.
  
  Possible values for SM2:
  \arg \e 0 there is not field expansion to yoke.
  \arg \e 1 field is expanded to yoke.
  */
  double flg2;    //!< Unused for Solenoid. Switch of new/old map for SM1. The field expansion in the yoke for SM2.
  int    curr;    //!< Current (applies only for SM2) in Amperes
  
  std::map<unsigned int, std::pair<double,double> > Polar; //!< Polarization of upstream and downstream cells. For solenoid only.
  
  //! Read target info, to retrieve polarization and update the array of the solenoid mag field to have it match the sign of its current.
  void readPolarization(CsField *Field, int Run = 0);
  
  //! Gives polarization of up and down cells for run \a run. Retuns false if no information about this run.
  bool getPolarization(unsigned int run, std::pair<double,double> &polar);

  /*! \fn void scaleSM2withNMR(CsField* Field, int Run = 0);
    \brief Get NMR for argument run# and update field map.
  */
  void scaleSM2withNMR(CsField *Field, int Run = 0);
  
  /*! \fn void scaleSolenoidMap(CsField *Field, float current);
    \brief Scale solenoid field map using argument "current"
  */
  void scaleSolenoidMap(CsField *Field, float current);
};

/*! \class CsField
 *  \brief Compass Magnetic Field Class.
 *
 *  It reads the informations related to the magnets from the Comgeant 
 *  detector data file output. It also give the components and gradient 
 *  of magnetic field in a space.
 */
class CsField {

  friend struct CsMagInfo;

  friend class PaField; 
  friend class PaEvent; 

 public:

  static float cvsRevision;              //!< Revision # of the CsField class

  CsField();                             //!< Default Constructor

  ~CsField();                            //!< Destructor

  CsField( const CsField& );             //!< Copy Constructor

  CsField& operator=( const CsField& );  //!< Assignment Operator

  /*! \fn void addMagInfo( CsMagInfo magp );
    \brief Add information about new magnet.
  */
  void addMagInfo( CsMagInfo magp );

  /*! \fn inline CsMagInfo* getMagInfo();
    \brief Returns the list of information structure available for the magnets.
  */
  inline CsMagInfo* getMagInfo() { return( magp_ ); }

  /*! \fn inline int getNumOfMags();
    \brief Returns the number of magnets
   */
  inline int getNumOfMags() { return( nmagp_ ); }

  /*! \fn bool ReadMaps();
    \brief Read maps of Solenoid, SM1 and SM2 from files definited in option 
    file. 

    The tag corresponded to field maps in option file is CsField. 
    Keys are SOL_field, SM1m_field, SM1m_field_measured 
    and SM2_field.

    Solenoid field map contains one half of the field of solenoid. It is the
    interval 0 < \e z  < 300 cm and 0 < \e r < 100 cm. 
    Inside the map the linear interpolation is used. The function \e fint2 is 
    responsible for this procedure.
    To receive another part of solenoid field (\e z < 0) the first part is mirrored. 
    The transverse component with respect to axe \e z have to change its sign.

    SM1 \b calculated field map contains only one quarter of SM1's field. It is the interval
    0 < \e x < 200 cm , 0 < \e y < 300 cm and 0 < \e z < 800 cm (magnet center corresponds
    to \e z = 350cm). 
    The grid size is: 10x10x10cm.
    For mirroring and interpolation we use the same procedure as in solenoid case.

    SM1 \b measured field map contains one half of SM1's field. The symmetry left/right    
    is supposed. This map is created in SM1 system. The map region is 0 < \e x < 160 cm,
    -160 < \e y < 160 cm, -416 < \e z < 440 cm. The grid size is 4x8x8cm.
    The first number in the map is the gap size of SM1. This value is used for
    crosscheck with information from \e detectors.dat file.
  */
  bool ReadMaps();


  /*! \fn bool getField( float pos_x, float pos_y, float pos_z,
    float& field_x, float& field_y, float& field_z, HepMatrix& grad);
    \brief Returns the magnetic field and gradient of that field in 
    a particular point in space.

    To find the value of field between nodes of grid linear interpolation is used. 

    For calculation of derivative the following formula is used:
    \f[
      \frac{\partial B(x,y,z)}{\partial x}=\frac{B(x+h,y,z)-B(x,y,z)}{h}
    \f]
    \param pos The position in space.
    \param field Return value: the magnetic field vector in \a pos.
    \param grad Return value: the gradient of magnetic field in \a pos.
    \f[
      grad =
       \left(
       \begin{array}{ccc}
         \partial B_x / \partial x & \partial B_x / \partial y & \partial B_x / \partial z \\
	 \partial B_y / \partial x & \partial B_y / \partial y & \partial B_y / \partial z \\
	 \partial B_x / \partial x & \partial B_x / \partial y & \partial B_z / \partial z \\
       \end{array}
       \right)
    \f]
  */
  bool getField( float pos_x, float pos_y, float pos_z, 
		 float& field_x, float& field_y, float& field_z,
		 CLHEP::HepMatrix& grad);

  /*! \fn bool getField( float pos_x, float pos_y, float pos_z,
    float& field_x, float& field_y, float& field_z);
    \brief Returns the magnetic field in a particular point in space.

    To find the value of field between nodes of grid linear interpolation is used. 

    \param pos The position in space.
    \param field Return value: the magnetic field vector in \a pos.
  */
  bool getField( float pos_x, float pos_y, float pos_z, 
		 float& field_x, float& field_y, float& field_z);

  // export \a this object into root field object \a rf  
//  void Export( RField * rf );

 private:

  CsMagInfo *magp_;    //!< array of magnet information structures.
  int       nmagp_;      //!< number of magnet information structures.

  // *************** SOLENOIDs ***************

  // ********** COMPASS **********
#define NZ_SOL 62
#define NR_SOL 22
  float zsol_[NZ_SOL];           //!< Z coordinate of points for Solenoid field.
  float rsol_[NR_SOL];           //!< R coordinate of points (cylindrical system) for Solenoid field.
  float bzsol_[NZ_SOL][NR_SOL];  //!< Bz: current, scaled, map
  float brsol_[NZ_SOL][NR_SOL];  //!< Br: current, scaled, map
  // Original (i.e. as read from file) maps, to be used as a fixed reference
  // when following the evolution of the solenoid field in a rotation run.
  float *bzsol0_, *brsol0_;      //!< BZ,Br: original map from file
  // N.B.: "sol" and "smc" (cf. infra) maps are mutually exclusive => Hence
  // the idea to allocate them only when needed. This is only implemented for
  // the original maps, which were added later.

  // ********** SMC **********
#define NZ_SMC 305
#define NR_SMC 59
  float zsmc_[NZ_SMC];           //!< Z coordinate of points for SMC solenoid field.
  float rsmc_[NR_SMC];           //!< R coordinate of points for SMC solenoid field. 
  float bzsmc_[NR_SMC][NZ_SMC];  //!< Bz: current, scaled, map
  float brsmc_[NR_SMC][NZ_SMC];  //!< Br: current, scaled, map
  // Original maps, cf. comments supra.
  float *bzsmc0_, *brsmc0_;      //!< BZ,Br: original map from file

  // ********** SOLENOIDs: GLOBAL **********
  float solMapCur_;              //!< Current corresponding to file map.
  // Uncertainty on current. For SMC, I(Y.B.) determine it from the deviation
  // from nominal observed in DB = .08A. Which turns out to correspond to
  // negligible (<.025%) tracking effect. => Use it as a scale to  single out
  // any anomaly worth noticing or taking into account.
  float soldCur_;                //!< Uncertainty on current


  // *************** TARGET DIPOLE: SMC ***************
#define NX_SMCD 42 
#define NY_SMCD 25 
#define NZ_SMCD 12 

  float xsmcd_[NX_SMCD];          //!< X coordinate of points for SMC dipole field. 
  float ysmcd_[NY_SMCD];          //!< Y coordinate of points for SMC dipole field. 
  float zsmcd_[NZ_SMCD];          //!< Z coordinate of points for SMC dipole field. 
  float bxsmcd_[NX_SMCD][NY_SMCD][NZ_SMCD];
  float bysmcd_[NX_SMCD][NY_SMCD][NZ_SMCD];
  float bzsmcd_[NX_SMCD][NY_SMCD][NZ_SMCD];

  // *************** TARGET DIPOLE: Oxford ***************
#define NX_OXD 27 
#define NY_OXD 27 
#define NZ_OXD 62
  
  float xoxd_[NX_OXD];          //!< X coordinate for the field of Oxford's dipole
  float yoxd_[NY_OXD];          //!< Y coordinate for the field of Oxford's dipole
  float zoxd_[NZ_OXD];          //!< Z coordinate for the field of Oxford's dipole
  float bxoxd_[NX_OXD][NY_OXD][NZ_OXD];
  float byoxd_[NX_OXD][NY_OXD][NZ_OXD];
  float bzoxd_[NX_OXD][NY_OXD][NZ_OXD];
  
  // *************** TARGET DIPOLE: SM1 ***************
#define NX_SM1c 21
#define NY_SM1c 31
#define NZ_SM1c 81

  float xsm1_[NX_SM1c];       //!< X coordinate of points for field of SM1.
  float ysm1_[NY_SM1c];       //!< Y coordinate of points for field of SM1.
  float zsm1_[NZ_SM1c];       //!< Z coordinate of points for field of SM1.
  float bxsm1_[NX_SM1c][NY_SM1c][NZ_SM1c]; //!< X components of the SM1 magnetic field.
  float bysm1_[NX_SM1c][NY_SM1c][NZ_SM1c]; //!< Y components of the SM1 magnetic field.
  float bzsm1_[NX_SM1c][NY_SM1c][NZ_SM1c]; //!< Z components of the SM1 magnetic field.

#define NX_SM1m 41
#define NY_SM1m 41
#define NZ_SM1m 108

  float xSM1_[NX_SM1m];       //!< X coordinate of points where field of SM1 was measured.
  float ySM1_[NY_SM1m];       //!< Y coordinate of points where field of SM1 was measured.
  float zSM1_[NZ_SM1m];       //!< Z coordinate of points where field of SM1 was measured.
  float bxSM1_[NX_SM1m][NY_SM1m][NZ_SM1m]; //!< X components of the SM1 measured magnetic field.
  float bySM1_[NX_SM1m][NY_SM1m][NZ_SM1m]; //!< Y components of the SM1 measured magnetic field.
  float bzSM1_[NX_SM1m][NY_SM1m][NZ_SM1m]; //!< Z components of the SM1 measured magnetic field.

  /*! \brief This and following array contain polypomial coefficients 
    that are result of fitting procedure for field of SM2 magnet.
  */
  float FSMAX1[1770],FSMAY1[2360],FSMAZ1[3540],FSMBX1[472],
    FSMBY1[300],FSMCX1[84],FSMCY1[126],FSMCZ1[168],FSMDX1[312],
    FSMDY1[390],FSMDZ1[468],FSMA01,FSMA11;
	
  int IFSMC1;                    //!< Current in SM2 (that value is taken from SM2 map file).

  /*! \fn bool getFieldSol( float pos_x, float pos_y, float pos_z,
    float& field_x, float& field_y, float& field_z);
    \brief Returns the Magnetic Filed of the COMPASS solenoid in a particular point in space.
    \param pos The position in space.
    \param field Return value: the magnetic field vector in \c pos.

    Inside the field map the linear interpolation is used.
    The function \e fint2 is called for that perpuse.
    If a given point is outside of the map range ( abs(\a z) > 300 or \a r > 100 cm ) 
    the interpolation for this coordinate is replaced by extrapolation
    suppressed by exponent.
  */
  bool getFieldSol( float pos_x, float pos_y, float pos_z,
		    float& field_x, float& field_y, float& field_z);

  /*! \fn bool getFieldSMC( float pos_x, float pos_y, float pos_z,
    float& field_x, float& field_y, float& field_z);
    \brief Returns the Magnetic Filed of the SMC solenoid in a particular point in space.
    \param pos The position in space.
    \param field Return value: the magnetic field vector in \c pos.
  */  
  bool getFieldSMC( float pos_x, float pos_y, float pos_z,
		    float& field_x, float& field_y, float& field_z);

  /*! \fn bool getFieldSMCdip( float pos_x, float pos_y, float pos_z, 
    float& field_x, float& field_y, float& field_z); 
    \brief Returns the Magnetic Filed of the SMC dipole in a particular point in space.
    \param pos The position in space. 
    \param field Return value: the magnetic field vector in \c pos. 
  */
  bool getFieldSMCdip( float pos_x, float pos_y, float pos_z,
		       float& field_x, float& field_y, float& field_z);

  /*! \fn bool getFieldDipOx( float pos_x, float pos_y, float pos_z, 
    float& field_x, float& field_y, float& field_z); 
    \brief Returns the Magnetic Filed of the Oxford's dipole in a particular point in space.
    \param pos The position in space. 
    \param field Return value: the magnetic field vector in \c pos. 
  */
  bool getFieldDipOx( float pos_x, float pos_y, float pos_z,
		      float& field_x, float& field_y, float& field_z);
  
  /*! \fn bool getFieldSm1( float pos_x, float pos_y, float pos_z,
    float& field_x, float& field_y, float& field_z);
    \brief Returns the Magnetic Filed of SM1 in a particular point in space.
    The old (calculated) map is used.
    \param pos The position in space.
    \param field Return value: the magnetic field vector in \c pos.
  */
  bool getFieldSm1( float pos_x, float pos_y, float pos_z,
		    float& field_x, float& field_y, float& field_z);

  /*! \fn bool getFieldSM1( float pos_x, float pos_y, float pos_z,
    float& field_x, float& field_y, float& field_z);
    \brief Returns the Magnetic Filed of SM1 in a particular point in space.
    The old (measured) map is used.
    \param pos The position in space.
    \param field Return value: the magnetic field vector in \c pos.
  */
  bool getFieldSM1( float pos_x, float pos_y, float pos_z,
		    float& field_x, float& field_y, float& field_z);

  /*! \fn bool getFieldSm2( float pos_x, float pos_y, float pos_z,
    float& field_x, float& field_y, float& field_z);
    \brief Returns the Magnetic Filed of SM2 in a particular point in space.
    \param pos The position in space.
    \param field Return value: the magnetic field vector in \c pos.
  */
  bool getFieldSm2( float pos_x, float pos_y, float pos_z,
		    float& field_x, float& field_y, float& field_z);

  /*! \fn bool getFieldSm2( float pos_x, float pos_y, float pos_z,
    float& field_x, float& field_y, float& field_z, int& IRET);
    \brief The same as the privies one but it does not extrapolate field 
    into the yoke.
    \param pos The position in space.
    \param field Return value: the magnetic field vector in \c pos.
    \param IRET Return the information about the current region.
           \arg \c IRET = 1 : In main box.
	   \arg \c IRET = 2 : Hit the magnet.
	   \arg \c IRET = 3 : Outside of magnetic field region.
  */
  bool getFieldSm2( float pos_x, float pos_y, float pos_z,
		    float& field_x, float& field_y, float& field_z, int& IRET);

  /*! \fn float fint2(float X[2],float A[4],float F[2][2]);
     \brief That function uses linear interpolaition method
     to evaluate a function  f( \a z, \a r) of 2 variables which has been
     tabulated at nodes of an 2-dimensional rectangular grid.
     \param X   contains the coordinates \a z, \a r of the point at which 
                the interpolation is to be performed.
     \param A   \e A[0] and \e A[1] contain tabulated values of \a z , 
                 \e A[2] and \e A[3] contain tabulated values of \a r.
     \param F   contains values of the function f at the nodes of 
                 the rectangular grid.
   */
  float fint2(float X[2],float A[4],float F[2][2]);

  /*! \fn float fint3(float X[3],float A[6],float F[2][2][2]);
     \brief That function uses linear interpolaition method
     to evaluate a function f( \a x, \a y, \a z) of 3 variables which has been
     tabulated at nodes of an 3-dimensional rectangular grid.
     \param X   contains the coordinates \a x, \a y, \a z of the point at which 
                the interpolation is to be performed.
     \param A   \e A[0] and \e A[1] contain tabulated values of \a x , 
                 \e A[2] and \e A[3] contain tabulated values of \a y,
		 \e A[4] and \e A[5] contain tabulated values of \a z.
     \param F   contains values of the function f at the nodes of 
                 the rectangular grid.
  */
  float fint3(float X[3],float A[6],float F[2][2][2]);
  
  /*! \fn float sign(float x);
    \brief Returns sign of the value x.
   */
  float sign(float x)
    { return (x < 0) ? -1 : 1; };
  
};

#  ifdef CsF_INIT_STATIC
float CsField::cvsRevision;
#    undef CsF_INIT_STATIC
#  endif

#endif // CsField_h
