// Revision: Adam Mielech 2001/25/11 in CsField::getFieldSm2
// removed rewriting FSMBY1 into CY1 --> FSMBY1 is directly set VY
// $Id: CsField.cc 13282 2012-03-19 14:43:23Z suhl $

/*!
   \file CsField.cc
   \brief Compass Magnetic Field Class.
   \author  Benigno Gobbo, Alexandre Korzenev
   \version $Revision: 13282 $
   \date    $Date: 2012-03-19 15:43:23 +0100 (Mon, 19 Mar 2012) $
*/

#include <string.h>
#include <stdio.h>
#include <cstdlib>   // for atoi()

#include "coral_config.h"
#include "CsInit.h"
#include "CsOpt.h"
#define   CsF_INIT_STATIC
#include "CsField.h"
#if USE_MySQL
# include "MySQLDB.h"
#endif

using namespace std;
using namespace CLHEP;

void CsMagInfo::readPolarization(CsField *Field, int Run)
{

  // - Polarisations AND solenoid current are retrieved.
  // - Either from DB or from file.
  // - If current is found to be at variance w/ the field retrieved from
  //  "detectors.dat", given known (value is built-in in the present software
  //  code, cf. "soldCur") precision on the measurement of the current):
  //   - Update field map arrays. The retrieved current is compared to the
  //    value associated w/ the field map. Which is built-in (cf. infra
  //    "solMapCur_" in "CsField::ReadMaps") and algebraic (hence allowing to
  //    derive the polarity of the field from the current value).
  //   - If, in addition, a change in the polarity of the current is detected:
  //    edit an ERROR message. For it might be that the "detectors.dat"
  //    specified in options file (which sets initially the sign of the
  //    current) was intended for the other polarity. (N.B.: This is only valid
  //    if the full path to the "detectors.dat" (as opposed to the path to
  //    a directory) has been specified: otherwise the polarity of the
  //    "detectors.dat" has been selected by "CsInit::GetSolenoidFieldSign",
  //    method which basically duplicates what's done here.)

  if (mag!=3)
    CsErrLog::mes(elFatal,"Solenoid must be the 1st in list of magnets (#3)!");

  float fscNew = fsc /* Save current field scaling factor */;
  float solCur = 0, solCurStatus = -1; bool solCurError = false;

  if (CsOpt::Instance()->getOpt("CsMagInfo","MySQLDB")) {

    // ******************** MySQL ********************
#if USE_MySQL
    if (!CsInit::Instance()->useMySQLDB())
      CsErrLog::mes
        (elFatal,"Can't read MySQLDB offline polarization without MySQLDB on!");

    pair<double,double> *tgtCurrents = 0;
    MySQLDB *_mysqldb = (MySQLDB*) CsInit::Instance()->getDB();
    bool connectfg = _mysqldb->isConnected();
    if (!connectfg) _mysqldb->ConnectDB();
    map<string,double> offpol = _mysqldb->getTgtOfflinePolar(Run);
    if (Run>33000) {  // I.e. if later or equal 2004
      tgtCurrents = new pair<double,double>;
      *tgtCurrents = _mysqldb->getTgtCurrents(Run);
    }
    if (!connectfg) _mysqldb->DisconnectDB();

    pair<double,double> polp(offpol["uppol"], offpol["downpol"]);
    Polar[Run] = polp;
    CsErrLog::msg(elWarning,__FILE__,__LINE__,"run #%d  up polar %f  down polar %f  solenoid %f",
		  Run,polp.first,polp.second,solCur);

    //         ********** GET TARGET CURRENTS **********

    // On 09/12/2007, Damien Neyret wrote: In fact, there are 2  sources of
    // informations on the solenoid current coming from the target system: one
    // associated with the polar measurements which is active only when the
    // polarizations can be measured (solenoid at ±2.5T, and 1T in 2006), and
    // one working permanently giving solenoid and dipole currents. Both
    // informations are stored in database in 2 different tables (tb_target and
    // tb_tgtcurrent). The 2 different values given by the logbook correspond
    // to these 2 sources, it's why one of the measurements can be different
    // from the other, the first one simply corresponds to the last polar
    // measurement, check the measurement date.
    // 
    // Concerning Coral the situation is a bit complicated, as a 3rd source of
    // information comes from offline polarization measurements. These
    // measurements are given by the target group some months after the end of
    // the data taking, and they are also associated with solenoid current
    // measurements (table tb_offlinepolar). There are 2 methods in Coral
    // MySQLDB to get the solenoid current: getTgtOfflinePolar and
    // getTgtCurrents (getTgtCurrentsTime is derivated from the last one). The
    // first one gives the polar measurements and solenoid current
    // corresponding to a given run by averaging all measurements found during
    // the run period in the offline table. If no entry is found during the run
    // period (like table not filled up for this year), the tb_target table is
    // used. There is here a risk, if the last polar measurement is too old the
    // solenoid current returned by this method can be false. A warning is
    // written to stderr output if the time difference is more than one hour.
    // 
    // The getTgtCurrents method in an other hand returns only solenoid and
    // dipole currents, it uses the tb_tgtcurrent table which is safer in that
    // point of view (but does not work for years < 2004). I think it would be
    // better to use this method for all magnet current measurements, and keep
    // getTgtOfflinePolar for polar measurements. I will see if I can improve
    // the behaviour of getTgtOfflinePolar to avoid such problems (like
    // returning 0 if values got from DB are not safe enough).

    solCur = offpol["solenoid"];
    solCurStatus = offpol["status"];
    // Status returned by "getTgtOfflinePolar":
    // 0 no value, 1 offline value, 2 online value, 3 old online value
    solCurError = solCurStatus==-1 || solCurStatus==0 || solCurStatus==3;

    if (Run<33000) {// If earlier than 2004: "getTgtOfflinePolar" only available
      if (solCurError) CsErrLog::mes(elFatal,"Solenoid current returned by \"getTgtOfflinePolar\" == 0. Suspicious!   =>  Supply solenoid current via FileDB.");
    }
    else {          // If later or equal 2004: use "getTgtCurrents"...
      // Still: compare the 2 values. And warn upon disagreement > 2*precision.
      if (fabs(tgtCurrents->first-solCur)>2*Field->soldCur_)
	CsErrLog::msg(elError,__FILE__,__LINE__,
	  "Run %d: Solenoid current in MySQL's \"tb_offlinePolar\" (=%.2f) and "
	  "\"tb_tgtcurrents\" (=%.2f) differ. The latter is retained.",
		      Run,solCur,tgtCurrents->first);
      solCur = tgtCurrents->first;
      float dipCur = tgtCurrents->second; // dipole value
      // "getTgtCurrents" returns no status. But 0 for both currents is
      // suspicious. => For the time being, as long as no explicit status word
      // is returned, consider it to be an error.
      solCurError = solCur==0 && dipCur==0;
      delete tgtCurrents;
      if (solCurError) CsErrLog::mes(elError,"Solenoid and dipole currents current returned by \"getTgtCurrents\" == 0. Suspicious!   => Supply target currents via FileDB.");
    }

    if (flg1==3) // Info exists for SMC dipole field, but not yet used...
      return;    // ...we rely on the "detectors.dat" value.

    Field->magp_[0].curr=int(1000.*solCur); // Save solenoid current in mA
    if (Field->solMapCur_) fscNew = solCur/Field->solMapCur_;
#else
    CsErrLog::mes
      (elError,"Can't read MySQLDB offline polarization without MySQL\n");
    CsErrLog::mes
      (elFatal,"Please reconfigure Coral with --with-MySQL option");
#endif
  }
  else {

    //       ******************** FileDB ********************

    fstream in; string path;
    if (!CsOpt::Instance()->getOpt("CsMagInfo","File",path)) {
      CsErrLog::mes(elError,
		    "Option \"CsMagInfo\" not found =>"
		    "Target polar and polarity unknown."); return;
    }
    in.open(path.c_str(),ios::in);
    if (!in ) CsErrLog::msg(elFatal,__FILE__,__LINE__,
			    "Error opening CsMagInfo file \"%s\"",path.c_str());

    // File is expected to have been created by
    //    "$CORAL/scripts/magInfo.pl"
    // Check this is indeed the case.
    // Skip first line == header == runnb	uppol	downpol	solencur
    string ctmp, label; in>>ctmp>>label;
    if (ctmp!="/*" || strncmp(label.c_str(),"magInfo.pl",strlen("magInfo.pl")))
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "CsMagInfo file \"%s\" NOT a \"magInfo.pl\" output file!",
		    path.c_str());

    int run; float field, dummy; pair<double,double> pol;
    bool found = false, comment = true;
    while (in.good() && !in.eof()) {
      in>>ctmp;
      if (comment) {
	if (ctmp=="*/") comment = false;
      }
      else if (ctmp=="/*") comment = true;
      else {
	run = atoi(ctmp.c_str());
	// One gets successively:
	//   type  soleno dipole     up       +/-     down       +/-  target#  nmr   nmr#     sm1   sm2
	in>>dummy>>solCur>>dummy>>pol.first>>dummy>>pol.second>>dummy>>dummy>>dummy>>dummy>>dummy>>dummy;
	if      (Run==0)     // What's run==0? Don't know!? No documentation...
	  Polar[run] = pol;  // ...Anyway, we then get the last entry in file.
	else if (Run==run) {
	  Polar[run] = pol; found = true;

	  CsErrLog::msg(elWarning,__FILE__,__LINE__,
			"run #%d  up polar %f  down polar %f  solenoid %f",
			Run,pol.first,pol.second,solCur);

	  if( flg1 == 3 ) break;                 // No info for SMC dipole field
	  Field->magp_[0].curr=int(1000.*solCur);// Save solenoid current in mA
	  if (Field->solMapCur_) fscNew = solCur/Field->solMapCur_;

	  break; // only one entry per run
	}
      }
    }
    if (Run!=0 && !found)
      CsErrLog::msg(elError,__FILE__,__LINE__,
   "No Target info for run #%d in \"%s\" => Target polar and polarity unknown.",
		    Run,path.c_str());
    in.close();
  }

  if (fsc==0) {
    // "fsc" is expected at this point to be the scaling factor read from
    // the geometry file => Therefore, "fsc==0" means that the solenoid map
    // has not been read in. And "fsc" has to be kept ==0 to convey the info
    // to the rest of coral.
    if (fabs(solCur)>Field->soldCur_)
      CsErrLog::msg(elError,__FILE__,__LINE__,
 "Solenoid current from MySQLDB %f for run #%d !=0 while scaling factor == 0 in detector-table => The latter may not be appropriate...",solCur,Run);
    return;
  }

  if (!solCurError /* This condition never occurs: cf. Fatal errors supra */ &&
      Field->solMapCur_ /* Require it to have been assigned a finite value */ &&
      fabs(fscNew-fsc)>Field->soldCur_/fabs(Field->solMapCur_)) {

    //     ***** CHANGE in SCALING FACTOR => UPDATE FIELD MAP *****

    if (fscNew*fsc<0)     // Change in polarity
      // "fsc" is expected at this point to be the scaling factor read from
      // the geometry file => warning.
      CsErrLog::msg(elError,__FILE__,__LINE__,
 "Solenoid scaling factor from MySQLDB %f for run #%d != detector-table's %f => The latter may not be appropriate...",fscNew,Run,fsc);
    scaleSolenoidMap(Field,solCur);
  }
}

void CsMagInfo::scaleSolenoidMap(CsField *Field, float current)
{
  if (Field->solMapCur_==0)
    CsErrLog::mes(elFatal,"\"solMapCur_\"==0 => Scaling impossible!");
  fsc = current/Field->solMapCur_;
  if (fabs(fsc-1)>.001) {
    printf("scaleSolenoidMap ==> Rescaling solenoid map by %.5f\n",fsc);
    if (flg1==1 || flg1==0) {                    // ***** COMPASS SOLENOID *****
      for (int j = 0; j<NR_SOL; j++) for (int i = 0; i<NZ_SOL; i++) {
	Field->bzsol_[i][j] = Field->bzsol0_[NR_SOL*i+j]*fsc;
	Field->brsol_[i][j] = Field->brsol0_[NR_SOL*i+j]*fsc;
      }
    }
    else if (flg1==2) {                          // ***** SMC SOLENOID *****
      for (int j=0; j<NR_SMC; j++) for (int i = 0; i<NZ_SMC; i++) {
	Field->bzsmc_[j][i] = Field->bzsmc0_[NZ_SMC*j+i]*fsc;
	Field->brsmc_[j][i] = Field->brsmc0_[NZ_SMC*j+i]*fsc;
      }
    }
    else
      CsErrLog::mes(elFatal,"Rescale not available for target dipole!");
  }
  else fsc = 1;  // If too small a rescaling: do not do anything
}

bool CsMagInfo::getPolarization(unsigned int run, pair<double,double> &pol )
{
  if( Polar.size() == 0 ) return false;
  map<unsigned int, pair<double,double> >::iterator ip = Polar.find( run );
  if( ip == Polar.end() ) return false;
  pol = (*ip).second;
  return true;
}

void CsMagInfo::scaleSM2withNMR(CsField* Field, int Run)
{
  // - SM2 NMR is retrieved.
  // - Either from DB or from file.
  // - If field scaling factor is found to differ from previous one:
  //   - Update field map arrays.

  if (mag!=5)
    CsErrLog::mes(elFatal,"SM2 must be the 3rd in list of magnets (#3)!");

  double sm2nmr = 0, rescaling = 1;
  if (!CsOpt::Instance()->getOpt("CsMagInfo","SM2",rescaling) || rescaling==0) {
    CsErrLog::mes(elError,
      "CsMagInfo rescaling factor == 0 => No taking into account NMR log.");
    return;
  }

  // ***** GET NMR either from DB or from file.... *****

  if (CsOpt::Instance()->getOpt("CsMagInfo","MySQLDB")) { // ...from MySQL *****

#if USE_MySQL
    if (!CsInit::Instance()->useMySQLDB())
      CsErrLog::mes
	(elFatal,"Can't read MySQLDB SM2 NMR without MySQLDB on!\n");

    MySQLDB* _mysqldb = (MySQLDB*) CsInit::Instance()->getDB();
    bool connectfg = _mysqldb->isConnected();
    if (!connectfg) _mysqldb->ConnectDB();
    sm2nmr = _mysqldb->getSM2NMR(Run);
    if (!connectfg) _mysqldb->DisconnectDB();

    if      (sm2nmr==0)
      CsErrLog::mes(elError,
		    "SM2 NMR not found in MySQL! => No SM2 rescaling.");
    else if (sm2nmr<0)
      CsErrLog::mes(elFatal,"Error reading SM2 NMR from MySQL!");
#else
    CsErrLog::mes
      (elError,"Can't read MySQLDB offline polarization without MySQL\n");
    CsErrLog::mes
      (elFatal,"Please reconfigure Coral with --with-MySQL option");
#endif
  }
  else {                                                   // ...from FILE
    fstream in; string path;

    if (!CsOpt::Instance()->getOpt("CsMagInfo","File",path)) {
      CsErrLog::mes(elError,
	"Option \"CsMagInfo\" not found => No SM2 rescaling."); return;
    }
    in.open(path.c_str(),ios::in);
    if (!in ) CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Error opening CsMagInfo file \"%s\"",path.c_str());

    // File is expected to have been created by
    //    "$CORAL/scripts/magInfo.pl"
    // Check this is indeed the case.
    string ctmp, label; in>>ctmp>>label;
    if (ctmp!="/*" || strncmp(label.c_str(),"magInfo.pl",strlen("magInfo.pl")))
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"CsMagInfo file \"%s\" NOT a \"magInfo.pl\" output file!",path.c_str());

    int run; float nmr, dummy; bool found = false, comment = true;
    while (in.good() && !in.eof()) {
      in>>ctmp;
      if (comment) {
	if (ctmp=="*/") comment = false;
      }
      else if (ctmp=="/*") comment = true;
      else {
	run = atoi(ctmp.c_str());
	// One gets successively:
	//  type  soleno dipole   up     +/-   down    +/-  target# nmr   nmr#   sm1    sm2
	in>>dummy>>dummy>>dummy>>dummy>>dummy>>dummy>>dummy>>dummy>>nmr>>dummy>>dummy>>dummy;
	if      (Run==0)    // What's run==0? Don't know!? No documentation...
	  sm2nmr = nmr;     // ...=> I reproduce what's done in readPolar supra.
	else if (Run==run) {
	  sm2nmr = nmr; break;
	}
      }
    }
    in.close();
    if (sm2nmr==0) CsErrLog::msg(elFatal,__FILE__,__LINE__,
        "SM2 NMR not found for run #%d in CsMagInfo File \"%s\" "
	"=> No SM2 rescaling.\nEither cancel the \"CsMagInfo SM2\" option "
	"or update your CsMagInfo File",
				 Run,path.c_str());
  }

  if (sm2nmr>0) {
    if ( fabs(rescaling-1.) > 1e-14 ) sm2nmr *= rescaling;

    // ***** DETERMINE FIELD SCALING FACTOR *****

    // Eric Weise wrote :
    //> I have the coordinates of the NMR probe relative to the center of SM2:
    const float x = -675./10; // mm -> cm
    const float y = -450./10; // mm -> cm
    const float z =   15./10; // mm -> cm
    //> where x direction is towards Jura side, y is upwards and z is downstream
    //> To my knowledge, the fieldmap gives a value of
    //>  By = 1.57856 T
    // I find 1.57904 T (Y.B.)
    // (The above concerns the 4000 A map. The value is any case retrieved from
    // the map itsef below.)

    float field_x, field_y, field_z;
    Field->getFieldSm2( x, y, z, field_x, field_y, field_z );
    const double sm2map = sqrt( (double)field_x*field_x + (double)field_y*field_y + (double)field_z*field_z );

  
    // ***** UPDATE FIELD SCALING FACTOR *****
    if ( fabs(rescaling-1.) > 1e-14 ) 
      CsErrLog::msg(elBasicInfo,__FILE__,__LINE__,
		    "SM2: NMR(=%.5f) * Option-supplied scale(=%.5f) / map(=%.5f) =>  Rescaling=%.5f while detectors.dat scale =%.5f",
		    sm2nmr/rescaling, rescaling, sm2map/fsc, fsc*sm2nmr/sm2map, fsc);
    else 
      CsErrLog::msg(elBasicInfo,__FILE__,__LINE__,
		    "SM2: NMR(=%.5f) / map(=%.5f) => Rescaling=%.5f  while detectors.dat scale =%.5f",
		    sm2nmr, sm2map/fsc, fsc*sm2nmr/sm2map, fsc);
    
    fsc *= sm2nmr/sm2map;
    
    // What should have been done: operate on the array of field components.
    // And do this consistently throughout CsField. Instead of having to
    // rescale the field value each time field is retrieved
    /*
      int i;
      for( i=0; i<1770; i++ ) Field->FSMAX1[i] *= fsc;
      for( i=0; i<472 ; i++ ) Field->FSMBX1[i] *= fsc;
      for( i=0; i<2360; i++ ) Field->FSMAY1[i] *= fsc;
      for( i=0; i<300 ; i++ ) Field->FSMBY1[i] *= fsc;
      for( i=0; i<3540; i++ ) Field->FSMAZ1[i] *= fsc;
      for( i=0; i<84  ; i++ ) Field->FSMCX1[i] *= fsc;
      for( i=0; i<126 ; i++ ) Field->FSMCY1[i] *= fsc;
      for( i=0; i<168 ; i++ ) Field->FSMCZ1[i] *= fsc;
      for( i=0; i<312 ; i++ ) Field->FSMDX1[i] *= fsc;
      for( i=0; i<390 ; i++ ) Field->FSMDY1[i] *= fsc;
      for( i=0; i<468 ; i++ ) Field->FSMDZ1[i] *= fsc;
    */
  }
}

CsField::CsField()
{
  nmagp_ = 0;
  magp_ = new CsMagInfo[3];
  // Assign a revision # to the CsField class.
  // - This is used by PHAST to document the mDSTs it outputs and so doing
  //  be able to ensure backward compatibility.
  // - It is stored in static data member "cvsRevision" of type float.
  // - It relies on cvs keyword substitution feature.
  char cvsRevString[] = "$Revision: 13282 $"; 
  if (sscanf(cvsRevString+10,  // i.e. starting right after the colon ":"
	     "%f",&CsField::cvsRevision)!=1)
    CsErrLog::mes(elFatal,
      "Error retrieving revision# from cvs \"$Revision: 13282 $\" keyword string");
}

CsField::~CsField()
{
  delete[] magp_;
}

CsField::CsField( const CsField& field ) : magp_(field.magp_) {
}

CsField& CsField::operator=( const CsField& field ) {
  if( this != &field ) {
    magp_ = field.magp_;
  }
  return( *this );
}

void CsField::addMagInfo( CsMagInfo magp ) {
  magp_[nmagp_] = magp;
  nmagp_++;
}

bool CsField::getField( float pos_x, float pos_y, float pos_z,
			float& field_x, float& field_y, float& field_z) {
  float mag_field_x,mag_field_y,mag_field_z;
  float x,y,z;

  field_x = 0;
  field_y = 0;
  field_z = 0;

  // mm -> cm
  pos_x = pos_x/10.;
  pos_y = pos_y/10.;
  pos_z = pos_z/10.;

  if( fabs(pos_z-magp_[0].zcm/10) < 400 ) {               //--- SOL + SM1 ---

    if( magp_[0].fsc != 0 ) {                    //-- SOLENOD --

      x=pos_x-magp_[0].xcm/10;
      y=pos_y-magp_[0].ycm/10;
      z=pos_z-magp_[0].zcm/10;

      if( magp_[0].flg1==1 || magp_[0].flg1==0 ) // OD/COMPASS SOLENOID
	getFieldSol( x, y, z, field_x, field_y, field_z);
      else if( magp_[0].flg1==2 )                // SMC SOLENOID
	getFieldSMC( x, y, z, field_x, field_y, field_z);
      else if( magp_[0].flg1==3 )                // SMC DIPOLE
	getFieldSMCdip( x, y, z, field_x, field_y, field_z);
      else if( magp_[0].flg1==4 )                // OD/COMPASS DIPOLE
	getFieldDipOx( x, y, z, field_x, field_y, field_z);
      else
	CsErrLog::mes(elFatal,
		      "Not clear COMPASS or SMC solenoid/dipole map to use.");
    }

    if( magp_[1].fsc != 0 ) {                    //-- SM1 --

      if(magp_[1].flg2==1.) {                    // measured

	x = pos_x - magp_[1].xcm/10;
	y = pos_y - magp_[1].ycm/10;
	z = pos_z - magp_[1].zcm/10;

	getFieldSM1( x, y, z, mag_field_x, mag_field_y, mag_field_z);

      } else if(magp_[1].flg2==0.) {             // calculated
	
	x = pos_x - magp_[1].xcm/10;
	y = pos_y - magp_[1].ycm/10;
	z = pos_z + 350 - magp_[1].zcm/10;

	getFieldSm1( x, y, z, mag_field_x, mag_field_y, mag_field_z);

      } else
	CsErrLog::mes(elFatal,
		      "Not clear MEASURED or CALCULATED SM1 map to use.");

      field_x += mag_field_x;
      field_y += mag_field_y;
      field_z += mag_field_z;
    }

  }else if( fabs(pos_z - magp_[1].zcm/10) < 450 ){        //--- SM1 ---

    if( magp_[1].fsc != 0 ) {

      if(magp_[1].flg2==1.) {                    // measured
	
	x = pos_x - magp_[1].xcm/10;
	y = pos_y - magp_[1].ycm/10;
	z = pos_z - magp_[1].zcm/10;
	
	getFieldSM1( x, y, z, field_x, field_y, field_z);

      }else if(magp_[1].flg2==0.) {              // calculated
	
	x = pos_x - magp_[1].xcm/10;
	y = pos_y - magp_[1].ycm/10;
	z = pos_z + 350 - magp_[1].zcm/10;
	
	getFieldSm1( x, y, z, field_x, field_y, field_z);

      } else
	CsErrLog::mes(elFatal,
		      "Not clear MEASURED or CALCULATED SM1 map to use.");
    }

  }else if( fabs(pos_z - magp_[2].zcm/10) < 350 ){        //--- SM2 ---

    x = pos_x - magp_[2].xcm/10;
    y = pos_y - magp_[2].ycm/10;
    z = pos_z - magp_[2].zcm/10;

    if( magp_[2].fsc != 0 )
      getFieldSm2( x, y, z, field_x, field_y, field_z);

  }

  return (true);
}

bool CsField::getField( float pos_x, float pos_y, float pos_z,
	       float& field_x, float& field_y, float& field_z,
	       HepMatrix& grad){
  int j;
  float displ=1;
  float field_N_x,field_N_y,field_N_z;

  getField(pos_x,pos_y,pos_z,field_x,field_y,field_z);

  for(j=0;j<3;j++){

    if(j==0)getField(pos_x+displ,pos_y,pos_z,field_N_x,field_N_y,field_N_z);
    else if(j==1)getField(pos_x,pos_y+displ,pos_z,field_N_x,field_N_y,field_N_z);
    else if(j==2)getField(pos_x,pos_y,pos_z+displ,field_N_x,field_N_y,field_N_z);

    grad[0][j]=(field_N_x-field_x)/displ;
    grad[1][j]=(field_N_y-field_y)/displ;
    grad[2][j]=(field_N_z-field_z)/displ;
  }

  return (true);
}

bool CsField::getFieldSm2( float pos_x, float pos_y, float pos_z,
			   float& field_x, float& field_y, float& field_z ) {
  int IRET=0,IERR=2;
  float dist;
  float pos1_x=0, pos1_y=0, pos1_z=0, pos2_x=0, pos2_y=0, pos2_z=0;
  float field1_x,field1_y,field1_z,field2_x,field2_y,field2_z;

  // CORAL -> COMGeant

  pos1_x = pos_z;
  pos1_y = pos_x;
  pos1_z = pos_y;


  if( !getFieldSm2( pos1_x, pos1_y, pos1_z,
		    field1_x, field1_y, field1_z, IRET) ) return false;		
	
  //--- Extrapolation required?

  if( magp_[2].flg2 != 0 ) {

    //--- Is the point somewhere in the yoke (but close to the fiducial area)?

    if( IRET == 2 ) {

      //---        Extrapolate

      pos1_x = 200.;
      pos1_y = 102.;
      pos1_z = 52.;

      if( fabs(pos_z) < pos1_x ) {
	pos2_x = pos_z;
	while( IERR == 2 ){
	  pos1_y = pos1_y - 0.5;
	  pos1_z = pos1_z - 0.5;
	  pos2_y = ((fabs(pos_x)<pos1_y)?fabs(pos_x):pos1_y) * sign(pos_x) ;
	  pos2_z = ((fabs(pos_y)<pos1_z)?fabs(pos_y):pos1_z) * sign(pos_y) ;



	  if( !getFieldSm2( pos2_x, pos2_y, pos2_z,
			    field2_x, field2_y, field2_z, IERR) ) return false;

			    		
	}
	dist = sqrt((pos2_x - pos_z) * (pos2_x - pos_z)+
		    (pos2_y - pos_x) * (pos2_y - pos_x)+
		    (pos2_z - pos_y) * (pos2_z - pos_y));
	if( dist < 150. ) {
	  field1_x = field2_x * exp(-dist*0.03);
	  field1_y = field2_y * exp(-dist*0.03);
	  field1_z = field2_z * exp(-dist*0.03);
	}
      }
    }
  }
  field_x = field1_y * magp_[2].fsc;
  field_y = field1_z * magp_[2].fsc;
  field_z = field1_x * magp_[2].fsc;

  return true;
}

bool CsField::getFieldSm2( float pos_x, float pos_y, float pos_z,
			   float& field_x, float& field_y, float& field_z,
			   int& IRET) {
// *** ALL INPUT/OUTPUT VARIABLES IN THE MAGNET COORDINATE SYSTEM
// *** THE TRANSLATION TO EMC IS MADE IN THE CALLING SEQUENCE
// *** INPUT/OUTPUT VARIABLE UNITS ARE METRES,GAUSS
// *** INTERNAL UNITS ARE CMS,GAUSS
// *** FOR A CALL ....   CALL MAGB(X,Y,Z,BX,BY,BZ)
// *** NOTE THAT X (MAGNET) = Y (EMC)
// *** NOTE THAT Y (MAGNET) = X (EMC)
// *** NOTE THAT Z (MAGNET) = -Z (EMC)

  float BXX[3+1][3+1],BYY[3+1][3+1],BZZ[3+1][3+1],BXXX[3+1],BYYY[3+1],
    BZZZ[3+1];
  double Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10;

  float XX,YY,ZZ,ZP,FC,X,Y,Z,BX,BY,BZ,XS,YS,ZS,DELY,
    DX,DY,XR,YR,VX,VY,VZ,A,B;
  int IFLAG,IX,JX,LX,MX,JX1,MX1,NX1,NY1,IY,JY,LY,JZ;

  // ------------------------------------------------------------------

  IRET=0;

  //---   Move to the magnet reference frame

  XX = pos_y;
  YY = pos_x;
  ZP = pos_z;

  //---   Here the old code is plugged in

  // *** CHANGE Z FOR EMC/MAGNET DIFFERENCE

  ZZ=-ZP;

  // *** FIELD CURRENT SIGN

  FC=-1.0;
  X=fabs(XX);
  Y=fabs(YY);
  Z=fabs(ZZ);


  //     *** HIT THE MAGNET ??

  if( Y <= 200  &&
      ( X > 101.5   ||  ( Y < 120.0 && Z > 50.0 )   ||
	( Y >= 120.0 && Z > 52.0 ) )    ) {
    IRET=2;
    field_x = 0;
    field_y = 0;
    field_z = 0;
    return true;
  }

  //     *** IN MAGNETIC FIELD REGION ??

  if( Y > 340.0  ||  Z > 69.0  ||  X > 144.0 ){
    IRET=3;
    field_x = 0;
    field_y = 0;
    field_z = 0;
    return true;
   }


  // *** PROTECT AGAINST UNDERFLOW (GENUINE ONES !!)

  if(Z < 1.0E-5) Z=0.0;

  IRET=0;
  field_x = 0.;
  field_y = 0.;
  field_z = 0.;


  //     *** COMPONENT SIGNS

  XS=sign(XX);
  ZS=sign(ZZ);

  //     *** IN MAIN BOX ??

  if( (X <= 101.5 || Y <= 200.)  &&
      (Z <= 50.0 || Y <= 200.0)  &&
      (Y < 222.0 && Z <= 52.0) ){
    //     *** MAIN FIELD - FIRST WORK OUT Y SHIFT AS FUNCTION OF X

    //     *** SET EXTRAPOLATION FLAG IF NECESSARY

    if(X > 96.0) IFLAG=1;

    //     *** LIMIT X IF NECESSARY

    if(X > 96.0) X=96.0;

    //     *** LIMIT Z IF NECESSARY

    if(Y > 120.0 && Z > 50.0){
      IFLAG=1;
      Z=50.0;
    }else{
      if(Y <= 120.0 && Z > 48.0){
	IFLAG=1;
	Z=48.0;
      }
    }
    DELY = FSMA01+FSMA11*(16.0-X/6.0)*(16.0-X/6.0);
    Y = -fabs(YY+DELY);

    //     *** Y COMPONENT SIGN

    YS = -sign(YY+DELY);

    //     *** WORK OUT X RANGE AND ARRAY NUMBER

    if( X <= 36. ){
      NX1=-1;
      DX=18.;
      XR=X;
    }else{
      if(X <= 72.){
	XR=X-36.0;
	NX1=(int)(XR/12.0);
	if(NX1 < 1) NX1=1;
	if(NX1 > 2) NX1=2;
	DX=12.0;
	XR=XR-(NX1-1)*12.0;
      }else{
	XR=X-54.0;
	NX1=(int)(XR/6.0);
	if(NX1 < 4) NX1=4;
	if(NX1 > 6) NX1=6;
	DX=6.0;
	XR=XR-(NX1-1)*6.0;
      }
    }

    //     *** WORK OUT Y RANGE AND ARRAY NUMBER

    if( Y <= -212. ){
      YR=Y+224.0;
      NY1=(int)(YR/6.0);
      if(NY1 < -1) NY1=-1;
      if(NY1 > 1) NY1=1;
      DY=6.0;
      YR=YR-(NY1-1)*6.0;
    }else{
      if(Y <= -80.){
	YR=Y+218.0;
	NY1=(int)(YR/3.0);
	if(NY1 < 3) NY1=3;
	if(NY1 > 45) NY1=45;
	DY=3.0;
	YR=YR-(NY1-1)*3.0;
      }else{
	if(Y <= -50.){
	  YR=Y+356.0;
	  NY1=(int)(YR/6.0);
	  if(NY1 < 47) NY1=47;
	  if(NY1 > 50) NY1=50;
	  DY=6.0;
	  YR=YR-(NY1-1)*6.0;
	}else{
	  YR=Y+662.0;
	  NY1=(int)(YR/12.0);
	  if(NY1 < 52) NY1=52;
	  if(NY1 > 55) NY1=55;
	  DY=12.0;
	  YR=YR-(NY1-1)*12.0;
	}
      }
    }

    //     *** FORM POWERS OF Z FOR LOCAL Z FITS

    Z2=Z*Z;
    Z3=Z2*Z;
    Z4=Z2*Z2;
    Z5=Z4*Z;
    Z6=Z5*Z;
    Z7=Z6*Z;
    Z8=Z7*Z;
    Z9=Z8*Z;
    Z10=Z9*Z;
    MX=NY1*30+NX1*3;
    MX1=NY1*8+2*NX1-12;
    LY=NY1+1;
    for( IY=0+1; IY<3+1; IY++){
      LY=LY+1;
      MX=MX+30;
      MX1=MX1+8;
      JX=MX;
      JX1=MX1;
      LX=NX1+1;
      for( IX=0+1; IX<3+1; IX++){
	LX=LX+1;
	JX=JX+3;
	JX1=JX1+2;
	JY=JX+JX/3;
	JZ=JX+JX;
	VX=FSMAX1[JX]*Z+FSMAX1[JX+1]*Z3+FSMAX1[JX+2]*Z5;
	VY=FSMAY1[JY]*Z+FSMAY1[JY+1]*Z3+FSMAY1[JY+2]*Z5+
	  FSMAY1[JY+3]*Z7;
	BZZ[IY][IX]=FSMAZ1[JZ]+FSMAZ1[JZ+1]*Z2+FSMAZ1[JZ+2]*Z4+
	  FSMAZ1[JZ+3]*Z6+FSMAZ1[JZ+4]*Z8+FSMAZ1[JZ+5]*Z10;
	
	//     *** EXTRA X TERMS
	
	if(LX > 6) VX=VX+FSMBX1[JX1]*Z7+FSMBX1[JX1+1]*Z9;
	
	//     *** EXTRA Y TERMS
	
	if(LY > 15 && LY < 46){
	 VY=VY+FSMBY1[10*(LY-16)+LX-1]*Z9;
	}
	BXX[IY][IX]=VX;
	BYY[IY][IX]=VY;
      }
    }
  }else{
    //     *** FRINGE FIELD EVALUATION

    //     *** Y COMPONENT SIGN

    YS=-sign(YY);

    //     *** MAKE NECESSARY Y SHIFT

    if(YY > 0.0) Y=-YY;
    if(YY < 0.0) Y=YY+2.0*(FSMA01+FSMA11*(16.0-X/6.0)*(16.0-X/6.0));

    //     *** ADJUST Y TO NEAREST FRINGE FIELD MEASUREMENT IF NECESSARY

    if(Y >= -220.0) IRET=1;
    if(Y > -220.0) Y=-220.0;

    //     *** FORM LOCAL Z POLYNOMIALS

    Z2=Z*Z;
    Z3=Z2*Z;
    Z4=Z2*Z2;
    Z5=Z3*Z2;
    Z6=Z4*Z2;
    if(Y <= -250.0){
      NX1=(int)(X/24.0-2.0);
      if(NX1 < -1) NX1=-1;
      if(NX1 > 3) NX1=3;
      YR=Y+340.0;
      NY1=(int)(YR/18.0-2.0);
      if(NY1 < -1) NY1=-1;
      if(NY1 > 2) NY1=2;
      MX=14*NY1+2*NX1;
      for( IY=0+1; IY<3+1; IY++){
	MX=MX+14;
	JX=MX;
	for( IX=0+1; IX<3+1; IX++){
	  JX=JX+2;
	  JY=JX+JX/2;
	  JZ=JX+JX;
	  BXX[IY][IX]=FSMCX1[JX]*Z+FSMCX1[JX+1]*Z3;
	  BYY[IY][IX]=FSMCY1[JY]*Z+FSMCY1[JY+1]*Z3+FSMCY1[JY+2]*Z5;
	  BZZ[IY][IX]=FSMCZ1[JZ]+FSMCZ1[JZ+1]*Z2+
	    FSMCZ1[JZ+2]*Z4+FSMCZ1[JZ+3]*Z6;
	}
      }
      DX=24.0;
      DY=18.0;
      YR=YR-(NY1+1)*18.0;
      XR=X-(NX1+1)*24.0;
    }else{
      NX1=(int)(X/12.0-2.0);
      if(NX1 < -1) NX1=-1;
      if(NX1 > 9) NX1=9;
      DX=12.0;
      YR=Y+250.0;
      NY1=(int)(YR/6.0-2.0);
      if(NY1 < -1) NY1=-1;
      if(NY1 > 2) NY1=2;
      DY=6.0;
      Z7=Z5*Z2;
      Z8=Z4*Z4;
      Z9=Z7*Z2;
      Z10=Z9*Z;
      MX=NY1*52+NX1*4;
      for( IY=0+1; IY<3+1; IY++){
	MX=MX+52;
	JX=MX;
	for( IX=0+1; IX<3+1; IX++){
	  JX=JX+4;
	  JY=JX+JX/4;
	  JZ=JX+JX/2;
	  BXX[IY][IX]=FSMDX1[JX]*Z+FSMDX1[JX+1]*Z3+
	    FSMDX1[JX+2]*Z5+FSMDX1[JX+3]*Z7;
	  BYY[IY][IX]=FSMDY1[JY]*Z+FSMDY1[JY+1]*Z3+
	    FSMDY1[JY+2]*Z5+FSMDY1[ JY+3]*Z7+FSMDY1[JY+4]*Z9;
	  BZZ[IY][IX]=FSMDZ1[JZ]+FSMDZ1[JZ+1]*Z2+
	    FSMDZ1[JZ+2]*Z4+FSMDZ1[JZ+3]*Z6+FSMDZ1[JZ+4]*Z8+
	    FSMDZ1[JZ+5]*Z10;
	}
      }
      YR=YR-(NY1+1)*6.0;
      XR=X-(NX1+1)*12.0;

      //     *** DO Y INTERPOLATION

    }
  }

  YR=YR/DY;
  for(IX=0+1;IX<3+1;IX++){   //IX=1,3
    A=(BXX[2+1][IX]+BXX[0+1][IX]-2.0*BXX[1+1][IX])/2.0;
    B=BXX[1+1][IX]-BXX[0+1][IX]-A;
    BXXX[IX]=(A*YR+B)*YR+BXX[0+1][IX];
    A=(BYY[2+1][IX]+BYY[0+1][IX]-2.0*BYY[1+1][IX])/2.0;
    B=BYY[1+1][IX]-BYY[0+1][IX]-A;
    BYYY[IX]=(A*YR+B)*YR+BYY[0+1][IX];
    A=(BZZ[2+1][IX]+BZZ[0+1][IX]-2.0*BZZ[1+1][IX])/2.0;
    B=BZZ[1+1][IX]-BZZ[0+1][IX]-A;
    BZZZ[IX]=(A*YR+B)*YR+BZZ[0+1][IX];
  }

  // *** DO X INTERPOLATION

  XR=XR/DX;
  A=(BXXX[2+1]+BXXX[0+1]-2.0*BXXX[1+1])/2.0;
  B=BXXX[1+1]-BXXX[0+1]-A;
  VX=(A*XR+B)*XR+BXXX[0+1];
  A=(BYYY[2+1]+BYYY[0+1]-2.0*BYYY[1+1])/2.0;
  B=BYYY[1+1]-BYYY[0+1]-A;
  VY=(A*XR+B)*XR+BYYY[0+1];
  A=(BZZZ[2+1]+BZZZ[0+1]-2.0*BZZZ[1+1])/2.0;
  B=BZZZ[1+1]-BZZZ[0+1]-A;
  VZ=(A*XR+B)*XR+BZZZ[0+1];
  // ***
  // *** BX,BY,BZ IN EMC COORDINATE SYSTEM
  // ***
  BX=VX*XS*ZS*FC;
  BY=VY*YS*ZS*FC;
  BZ=-VZ*FC;

  //---     Move to the MRS (master reference frame)
  //---     Convert to kGs

  field_x = -BY*0.0001;
  field_y = -BX*0.0001;
  field_z = -BZ*0.0001;

  return true;
}

bool CsField::ReadMaps()
{

  // ****************** READ MAPS FROM FILE (or swicth to DB) ******************

  char title[54]; int i,j,k; float x,y,z,r,bx,by,bz,br,gap,tmp;
  ifstream in; string path,msg;
  CsOpt* opt = CsOpt::Instance();

  // Switch data source to database
  string from; opt->getOpt("CsField","From",from);
  if (from==string("database")) {
    return true;
  }

  if (magp_[0].fsc!=0) {        // *************** TARGET FIELDs ***************

    //                                           ********** SOLENOIDs **********
    // 2 copies of the field map are stored
    //   - Original values read as is from file: they will later serve as a
    //    reference when it comes to rescale the field, e.g. during a target
    //    rotation.
    //     Corresponding memmory is allocated in this->method in order to profit
    //    from the fact that "sol" and "smc" are mutually exclusive to gain some
    //    space. (The scaled maps, cf. infra, were implemented earlier, and
    //    are allocated in any case, in the constructor.)
    //   - Scaled values (for the current scale factor).
    if (!opt->getOpt( "CsField", "SOL_field", path))
      CsErrLog::mes(elFatal,
		    "Tag `CsField` and key `SOL_field` were not found.");
    in.open(path.c_str(), ios::in);
    if (!in) CsErrLog::msg(elFatal,__FILE__,__LINE__,
			   "Map file for Solenoid \"%s\" was not found",
			   path.c_str());

    if (magp_[0].flg1==1 || magp_[0].flg1==0 ) { // ***** COMPASS SOLENOID *****

      cout<<"+++ ReadMaps: Reading COMPASS solenoid map... "<<flush;
      char TITRE1[7]; in>>TITRE1;
      if (strcmp(TITRE1,"MAGNET"))
	CsErrLog::mes(elFatal,"\nWrong COMPASS solenoid map!\n");
      bzsol0_ = new float[NR_SOL*NZ_SOL]; brsol0_ = new float[NR_SOL*NZ_SOL];
      in.seekg(0x1a0,ios::beg);

      for (j = 0; j<NR_SOL-1; j++) for (i = 0; i<NZ_SOL-1; i++) {
	in>>z>>r>>bz>>br>>tmp>>tmp>>tmp;
	zsol_[i] = z; rsol_[j] = r;
	bzsol0_[NR_SOL*i+j] = bz;        brsol0_[NR_SOL*i+j] = br;
	bzsol_ [i][j] = bz*magp_[0].fsc; brsol_ [i][j] = br*magp_[0].fsc;
      }
      cout<<endl;

      solMapCur_ = 646.5;  // Cf. e-mail from Frabrice of 04/10/2007
      // Uncertainty on current. For SMC, I(Y.B.) determined it from the
      // deviation from nominal observed in DB = .08A. Assuming same precision
      // for the COMPASS/OD case and rounding up =>
      soldCur_ = .1;
    }
    else if (magp_[0].flg1==2) {                     // ***** SMC SOLENOID *****

      cout<<"+++ ReadMaps: Reading SMC solenoid map... "<<flush;
      char TITRE1[11],TITRE2[4]; in>>TITRE1>>TITRE2;
      if (strcmp(TITRE1,"*SOLENOIDE") || strcmp(TITRE2,"SMC"))
	CsErrLog::mes(elFatal,"\nWrong SMC solenoid map!\n");
      bzsmc0_ = new float[NR_SMC*NZ_SMC]; brsmc0_ = new float[NR_SMC*NZ_SMC];
      for (j = 0; j<NR_SMC; j++) for (i = 0; i<NZ_SMC; i++) {
	in>>br>>bz;
	bzsmc0_[NZ_SMC*j+i] = bz;        brsmc0_[NZ_SMC*j+i] = br;
	bzsmc_ [j][i] = bz*magp_[0].fsc; brsmc_ [j][i] = br*magp_[0].fsc;
      }
      cout<<endl;

      // SET UP COORDINATE SYSTEM FOR SMC SOLENOID FIELD
      for( i=0; i<NR_SMC; i++ ) {
	if( i <= 7 )            r = (float)(i)*20.;
	if( i > 7  && i <= 47 ) r = (float)(i-7)*1.+7.*20.;
	if( i > 47 && i <= 58 ) r = (float)(i-47)*20.+7.*20.+40.*1.;
	rsmc_[i]=r/10.;           /*  mm->cm  */
      }
      for( i=0; i<NZ_SMC; i++ ) {
	if( i <= 18 )             z = (float)(i)*50.;
	if( i > 18  && i <= 28  ) z = (float)(i-18)*10.+18.*50.;
	if( i > 28  && i <= 228 ) z = (float)(i-28)*1.+18.*50.+10.*10.;
	if( i > 228 && i <= 278 ) z = (float)(i-228)*10.+18.*50.+10.* 10.+200.*1.;
	if( i > 278 && i <= 304 ) z = (float)(i-278)*50.+18.*50.+10.* 10.+200.*1.+50.*10.;
	zsmc_[i]=z/10.;           /*  mm->cm  */
      }
      solMapCur_ = -417;  // Or so says Fabrice
      // Uncertainty on current. For SMC, I(Y.B.) determine it from the
      // deviation from nominal observed in DB = .08A. Rounding up =>
      soldCur_ = .1;
    }
    else if (magp_[0].flg1==3) {                       // ***** SMC DIPOLE *****
	
      char TITRE1[10],TITRE2[10];
      double XSTEP,YSTEP,ZSTEP;
      float XMINDP,XMAXDP, YMINDP,YMAXDP, ZMINDP,ZMAXDP;
      int NDIPOX,NDIPOY,NDIPOZ;
	
      cout<<"+++ ReadMaps: Reading SMC dipole map..."<<flush;

      in>>TITRE1>>TITRE2;
      if (strcmp(TITRE1,"*SMC") || strcmp(TITRE2,"DIPOLE"))
	CsErrLog::mes(elFatal,"\nWrong SMC dipole map!\n");

      in>>NDIPOY>>YMINDP>>YMAXDP;
      in>>NDIPOX>>XMINDP>>XMAXDP;
      in>>NDIPOZ>>ZMINDP>>ZMAXDP;
      if (NX_SMCD!=NDIPOX)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "\nWrong X for SMC dipole map. NX_SMCD %i",NDIPOX);
      if (NY_SMCD!=NDIPOY )
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "\nWrong Y for SMC dipole map. NY_SMCD %i",NDIPOY);
      if (NZ_SMCD!=NDIPOZ)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "\nWrong Z for SMC dipole map. NZ_SMCD %i",NDIPOZ);

      for (k = 0; k<NZ_SMCD; k++) for (i = 0; i<NX_SMCD; i++)
	for (j = 0; j<NY_SMCD; j++) {
	  in>>bysmcd_[i][j][k]>>bxsmcd_[i][j][k]>>bzsmcd_[i][j][k];
	  bxsmcd_[i][j][k] *= magp_[0].fsc;
	  bysmcd_[i][j][k] *= magp_[0].fsc;
	  bzsmcd_[i][j][k] *= magp_[0].fsc;
	}
      cout<<endl;

      // Build space information
      XSTEP=(XMAXDP-XMINDP)/NDIPOX;
      for (i = 0; i<NDIPOX; i++) {
	xsmcd_[i] = XMINDP + i*XSTEP + 0.5*XSTEP;
	xsmcd_[i] = xsmcd_[i]/1000.;
      }
      YSTEP=(YMAXDP-YMINDP)/NDIPOY;
      for (j = 0; j<NDIPOY; j++) {
	ysmcd_[j] = YMINDP + j*YSTEP + 0.5*YSTEP;
	ysmcd_[j] = ysmcd_[j]/1000.;
      }
      ZSTEP=(ZMAXDP-ZMINDP)/NDIPOZ;
      for (k = 0; k<NDIPOZ; k++) {
	zsmcd_[k] = ZMINDP + k*ZSTEP + 0.5*ZSTEP;
	zsmcd_[k] = zsmcd_[k]/1000.;
      }
      // Uncertainty on solenoid current: it's needed to be able to interpret
      // the value returned by the DB, be it =0. As to its value, cf. supra.
      soldCur_ = .1;
    }
    else if (magp_[0].flg1==4) {               // ***** OD/COMPASS dipole *****

      cout<<"+++ ReadMaps: Reading OD/COMPASS dipole field map... "<<endl;
      char TITRE1[40],TITRE2[4]; in>>TITRE1>>TITRE2;
      if (strcmp(TITRE1,"Map") || strcmp(TITRE2,"for"))
	CsErrLog::mes(elFatal,"\nWrong OD/COMPASS dipole map!\n");
      
      for (int i = 0; i<21; i++) in>>TITRE1;
      
      for (int i = 0; i<NX_OXD; i++) {
	for (int j = 0; j<NY_OXD; j++) {
	  for (int k = 0; k<NZ_OXD; k++) {
	    in>>xoxd_[i]>>yoxd_[j]>>zoxd_[k];
	    in>>bxoxd_[i][j][k]>>byoxd_[i][j][k]>>bzoxd_[i][j][k];
	    //printf("%5.1f %5.1f %5.1f  %8.5f %8.5f %8.5f\n",xoxd_[i],yoxd_[j],zoxd_[k], bxoxd_[i][j][k],byoxd_[i][j][k],bzoxd_[i][j][k]);
	    bxoxd_[i][j][k] *= magp_[0].fsc;
	    byoxd_[i][j][k] *= magp_[0].fsc;
	    bzoxd_[i][j][k] *= magp_[0].fsc;
	  }
	}
      }
      // Uncertainty on solenoid current: it's needed to be able to interpret
      // the value returned by the DB, be it =0. As to its value, cf. supra.
      soldCur_ = .1;
    }
    
    in.close();
  }

  if( magp_[1].fsc != 0 ){

    if(magp_[1].flg2==1.){

      //------- NEW SM1 ------------------
      if(magp_[1].flg1!=172 && magp_[1].flg1!=132 && magp_[1].flg1!=82){
	CsErrLog::mes(elFatal,
		      "Possible gap sizes for SM1 magnet: 172cm, 132cm, 82cm.");
	msg =  "The gap size read from detectors.dat is ";
	char GAP[10];
	sprintf(GAP, "%.1f",magp_[1].flg1);
	msg += GAP;
	msg += "cm.  You are used old version of COMGeant or detectors.dat was edit.";
	CsErrLog::mes( elFatal,msg);
      }

      if( !opt -> getOpt( "CsField", "SM1m_field_measured",path) )
	CsErrLog::mes(elFatal,
		      "Tag `CsField` and key `SM1m_field_measured` were not found. ");
      in.open( path.c_str(), ios::in);
      if(!in ){
	msg =  "Map file for SM1 ";
	msg += path;
	msg += " was not found.";
	CsErrLog::mes( elFatal, msg);
	return (false);
      }
      else{
	CsErrLog::mes( elInfo, "+++ ReadMaps: Reading new (measured) SM1 map... ");
	in>>gap;
	
	if(gap!=magp_[1].flg1){
	  CsErrLog::mes(elFatal,"The gap size of SM1 in detectors.dat does't correspond to gap size in field map file. Check path in option file: `CsField SM1m_field_measured`.");
	}

	for( i=0; i<41; i++ ){
	  for( j=0; j<41; j++ ){
	    for( k=0; k<108; k++ ){
	      in>>x>>y>>z>>bx>>by>>bz;
	      xSM1_[i] = x ; // cm
	      ySM1_[j] = y ; // cm
	      zSM1_[k] = z ; // cm note that map offset = 350
	      bxSM1_[i][j][k] = bx * magp_[1].fsc; // T
	      bySM1_[i][j][k] = by * magp_[1].fsc; // T
	      bzSM1_[i][j][k] = bz * magp_[1].fsc; // T
	      //cout<<x<<" "<<y<<" "<<z<<" "<<bx<<" "<<by<<" "<<bz<<endl;
	    }
	  }
	}
	
	in.close();
      }
    } else if(magp_[1].flg2==0.){

      //------- OLD SM1 ------------------
      if( !opt -> getOpt( "CsField", "SM1m_field",path) )
	CsErrLog::mes(elError,
		      "Tag `CsField` and key `SM1m_field` were not found. ");
      in.open( path.c_str(), ios::in);
      if(!in ){
	msg =  "Map file for SM1 ";
	msg += path;
	msg += " was not found.";
	CsErrLog::mes( elFatal, msg);
	return (false);
      }
      else{
	CsErrLog::mes( elInfo, "+++ ReadMaps: Reading old (calculated) SM1 map... ");
	in.seekg( 0x49, ios::beg);
	
	for( i=0; i<21; i++ ){
	  for( j=0; j<31; j++ ){
	    for( k=0; k<81; k++ ){
	      in>>x>>y>>z>>bx>>by>>bz>>tmp;
	      xsm1_[i] = x * 100; // m -> cm
	      ysm1_[j] = y * 100; // m -> cm
	      zsm1_[k] = z * 100; // m -> cm note that map offset = 350
	      bxsm1_[i][j][k] = bx * magp_[1].fsc; // T
	      bysm1_[i][j][k] = by * magp_[1].fsc; // T
	      bzsm1_[i][j][k] = bz * magp_[1].fsc; // T
	    }
	  }
	}
	in.close();
      }
    } else {
      CsErrLog::mes(elError,
		    "Not clear new or old field map must be taken for SM1.");
      CsErrLog::mes(elFatal,"Please, check file detectors.dat, string `mag   4`. `flag2` must be equal to 0(old map) or 1(new map).");
    }
  }

  if (magp_[2].fsc != 0) {

    //         ******************** SM2 ********************

    //                                      ***** GET MAP FILE NAME FROM OPTIONS
    if (!opt->getOpt("CsField","SM2_field",path))
      CsErrLog::mes(elFatal,
		    "Option \"CsField SM2_field\" not found in options file!");
    in.open( path.c_str(), ios::in );
    if (!in)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "File \"%s\", specified as magnetic map file for SM2, "
		    "does not exist!",path.c_str());
    //       ***** CHECK ELECTRIC CURRENT from MAP AGAINST CURRENT from GEOMETRY
    CsErrLog::mes(elInfo,"+++ ReadMaps: Reading SM2 map...\n");
    in>>IFSMC1;
    if (IFSMC1!=magp_[2].curr)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "Inconsistent electrical currents in SM2 in "
		    "detectors.dat (=%dA) and magnetic map (=%dA)\n   => Check "
		    "your detectors.dat against map file \"%s\"",
		    magp_[2].curr,IFSMC1,path.c_str());

    for (i = 0; i<54; i++) in>>title[i];
    in>>FSMA01;
    in>>FSMA11;
    for (i=0; i<1770; i++) in>>FSMAX1[i];
    for (i=0; i<472 ; i++) in>>FSMBX1[i];
    for (i=0; i<2360; i++) in>>FSMAY1[i];
    for (i=0; i<300 ; i++) in>>FSMBY1[i];
    for (i=0; i<3540; i++) in>>FSMAZ1[i];
    for (i=0; i<84  ; i++) in>>FSMCX1[i];
    for (i=0; i<126 ; i++) in>>FSMCY1[i];
    for (i=0; i<168 ; i++) in>>FSMCZ1[i];
    for (i=0; i<312 ; i++) in>>FSMDX1[i];
    for (i=0; i<390 ; i++) in>>FSMDY1[i];
    for (i=0; i<468 ; i++) in>>FSMDZ1[i];
    in.close();
  }
  return (true);
}

inline float CsField::fint2( float X[2], float A[4], float F[2][2]) {
  double t,u,g1,g2;
  double x,y;
  
  x = X[0];
  y = X[1];

  if (A[1]==A[0]) {
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "Z coordinates to extrapolate from are = (%f)!",A[0]);
  }

  t = (x-A[0])/(A[1]-A[0]);
  u = (y-A[2])/(A[3]-A[2]);

  g1 = (1-t)*F[0][0] + t*F[1][0];
  g2 = (1-t)*F[0][1] + t*F[1][1];

  return (1-u)*g1+u*g2;
}

inline float CsField::fint3( float X[3], float A[6], float F[2][2][2])
{
  double t,u,v,g1,g2,f1,f2,h1,h2;
  double x,y,z;

  x = X[0];
  y = X[1];
  z = X[2];

  if( A[1] == A[0] ){
    CsErrLog::mes(elFatal,"fint3: Coordinates for different points along axe X have to be different.");
  }
  if( A[3] == A[2] ){
    CsErrLog::mes(elFatal,"fint3: Coordinates for different points along axe Y have to be different.");
  }
  if( A[5] == A[4] ){
    CsErrLog::mes(elFatal,"fint3: Coordinates for different points along axe Z have to be different.");
  }

  t = (x-A[0])/(A[1]-A[0]);
  u = (y-A[2])/(A[3]-A[2]);
  v = (z-A[4])/(A[5]-A[4]);

  g1 = (1-t)*F[0][0][0]+t*F[1][0][0];
  g2 = (1-t)*F[0][1][0]+t*F[1][1][0];
  h1 = (1-u)*g1+u*g2;

  f1 = (1-t)*F[0][0][1]+t*F[1][0][1];
  f2 = (1-t)*F[0][1][1]+t*F[1][1][1];
  h2 = (1-u)*f1+u*f2;

  return (1-v)*h1+v*h2;
}



bool CsField::getFieldSMCdip( float pos_x, float pos_y, float pos_z,
			      float& field_x, float& field_y, float& field_z) {


  double XTARG1,YTARG1,ZTARG1,COSRG1,SINRG1;

  double XYZ[3],BXYZE[3]={0,0,0};
  double IRSIGX,IRSIGY,IRSIGZ;
  double X,Y,Z;
  double HX,HY,HZ;
  double HXA,HYA,HZA;
  double HXB,HYB,HZB;
  double HXC,HYC,HZC;
  double HXD,HYD,HZD;
  double HXE,HYE,HZE;
  double HXF,HYF,HZF;
  double HXG,HYG,HZG;
  double HXH,HYH,HZH;
  double HX1AB,HY1AB,HZ1AB;
  double HX1BC,HY1BC,HZ1BC;
  double HX1CD,HY1CD,HZ1CD;
  double HX1AD,HY1AD,HZ1AD;
  double HX1EH,HY1EH,HZ1EH;
  double HX1FG,HY1FG,HZ1FG;
  double XSSU,XSSL,XUL;
  double YSSU,YSSL,YUL;
  double ZSSU,ZSSL,ZUL;
  int IX,IY,IZ;
  int IUX,IUY,IUZ;
  int ILX,ILY,ILZ;

  IRSIGX = 1;
  IRSIGY = 1;
  IRSIGZ = 1;

  XYZ[0] = pos_z;   // Z -> X
  XYZ[1] = pos_x;   // X -> Y
  XYZ[2] = pos_y;   // Y -> Z

  field_x = 0;
  field_y = 0;
  field_z = 0;

  XTARG1=0;
  YTARG1=0;
  ZTARG1=0;
  COSRG1=1;
  SINRG1=0;

  X = XYZ[0]/100-XTARG1;
  Y = XYZ[1]/100-YTARG1;
  Z = XYZ[2]/100-ZTARG1;

  X=X*COSRG1-Y*SINRG1;
  Y=Y*COSRG1+X*SINRG1;
  //
  //     Remember sign
  //
  if( X < 0 ) IRSIGX = -1;
  if( Y < 0 ) IRSIGY = -1;
  if( Z < 0 ) IRSIGZ = -1;
  X=fabs(X);
  Y=fabs(Y);
  Z=fabs(Z);
  //
  //     Search for the equivalent X coordinates in field map
  //
  if( X >= xsmcd_[NX_SMCD-1] ) return false;

  for( IX=0; IX<NX_SMCD; IX++ )
    if( xsmcd_[IX] > X ) break;

  IUX=IX;
  ILX=IUX-1;
  //
  //     Search for the equivalent Y coordinates in field map
  //
  if( Y >= ysmcd_[NY_SMCD-1] ) return false;

  for( IY=0; IY<NY_SMCD; IY++ )
    if( ysmcd_[IY] > Y ) break;

  IUY=IY;
  ILY=IUY-1;
  //
  //     Search for the equivalent Z coordinates in field map
  //
  if( Z >= zsmcd_[NZ_SMCD-1] ) return false;

  for( IZ=0; IZ<NZ_SMCD; IZ++ )
    if( zsmcd_[IZ] > Z ) break;

  IUZ=IZ;
  ILZ=IUZ-1;

  //
  //     8 coordinate points given by ( 8 corner of a box= nodes)
  //        X       Y        Z
  //        ILX    ILY      ILZ
  //        ILX    ILY      IUZ
  //        ILX    IUY      ILZ
  //        ILX    IUY      IUZ
  //        IUX    ILY      ILZ
  //        IUX    ILY      IUZ
  //        IUX    IUY      ILZ
  //        IUX    IUY      IUZ
  //     are points around desired point   (X,Y,Z)
  //     where field values stored (from the field map).
  //
  //
  //     Distancy in X between the nodes
  //
  XSSU = xsmcd_[IUX]-X;
  XSSL = X-xsmcd_[ILX];
  XUL  = XSSU+XSSL;
  //
  //     Distancy in Y between the nodes
  //
  YSSU = ysmcd_[IUY]-Y;
  YSSL = Y-ysmcd_[ILY];
  YUL  = YSSU+YSSL;
  //
  //     Distancy in Z between the nodes
  //

  ZSSU = zsmcd_[IUZ]-Z;
  ZSSL = Z-zsmcd_[ILZ];
  ZUL  = ZSSU+ZSSL;
  //
  //     interpolate in X to find field components HA
  //     at (X,YDIPOL(ILY),ZDIPOL(ILZ))
  //
  HXA=(bxsmcd_[ILX][ILY][ILZ]*XSSU+bxsmcd_[IUX][ILY][ILZ]*XSSL)/XUL;
  HYA=(bysmcd_[ILX][ILY][ILZ]*XSSU+bysmcd_[IUX][ILY][ILZ]*XSSL)/XUL;
  HZA=(bzsmcd_[ILX][ILY][ILZ]*XSSU+bzsmcd_[IUX][ILY][ILZ]*XSSL)/XUL;

  //
  //     interpolate in X to find field components HB
  //     at (X,YDIPOL(ILY),ZDIPOL(IUZ))
  //

  HXB=(bxsmcd_[ILX][ILY][IUZ]*XSSU+bxsmcd_[IUX][ILY][IUZ]*XSSL)/XUL;
  HYB=(bysmcd_[ILX][ILY][IUZ]*XSSU+bysmcd_[IUX][ILY][IUZ]*XSSL)/XUL;
  HZB=(bzsmcd_[ILX][ILY][IUZ]*XSSU+bzsmcd_[IUX][ILY][IUZ]*XSSL)/XUL;
  //
  //     interpolate in X to find field components HC
  //     at (X,YDIPOL(IUY),ZDIPOL(IUZ))
  //

  HXC=(bxsmcd_[ILX][IUY][IUZ]*XSSU+bxsmcd_[IUX][IUY][IUZ]*XSSL)/XUL;
  HYC=(bysmcd_[ILX][IUY][IUZ]*XSSU+bysmcd_[IUX][IUY][IUZ]*XSSL)/XUL;
  HZC=(bzsmcd_[ILX][IUY][IUZ]*XSSU+bzsmcd_[IUX][IUY][IUZ]*XSSL)/XUL;
  //
  //     interpolate in X to find field components HD
  //     at (X,YDIPOL(IUY),ZDIPOL(ILZ))
  //

  HXD=(bxsmcd_[ILX][IUY][ILZ]*XSSU+bxsmcd_[IUX][IUY][ILZ]*XSSL)/XUL;
  HYD=(bysmcd_[ILX][IUY][ILZ]*XSSU+bysmcd_[IUX][IUY][ILZ]*XSSL)/XUL;
  HZD=(bzsmcd_[ILX][IUY][ILZ]*XSSU+bzsmcd_[IUX][IUY][ILZ]*XSSL)/XUL;
  //
  //     interpolate in Y to find field components HE
  //     at (XDIPOL(ILX),Y,ZDIPOL(ILZ))
  //

  HXE=(bxsmcd_[ILX][ILY][ILZ]*YSSU+bxsmcd_[ILX][IUY][ILZ]*YSSL)/YUL;
  HYE=(bysmcd_[ILX][ILY][ILZ]*YSSU+bysmcd_[ILX][IUY][ILZ]*YSSL)/YUL;
  HZE=(bzsmcd_[ILX][ILY][ILZ]*YSSU+bzsmcd_[ILX][IUY][ILZ]*YSSL)/YUL;
  //
  //     interpolate in Y to find field components HH
  //     at (XDIPOL(ILX),Y,ZDIPOL(IUZ))
  //

  HXH=(bxsmcd_[ILX][ILY][IUZ]*YSSU+bxsmcd_[ILX][IUY][IUZ]*YSSL)/YUL;
  HYH=(bysmcd_[ILX][ILY][IUZ]*YSSU+bysmcd_[ILX][IUY][IUZ]*YSSL)/YUL;
  HZH=(bzsmcd_[ILX][ILY][IUZ]*YSSU+bzsmcd_[ILX][IUY][IUZ]*YSSL)/YUL;
  //
  //     interpolate in Y to find field components HG
  //     at (XDIPOL(IUX),Y,ZDIPOL(IUZ))
  //

  HXG=(bxsmcd_[IUX][ILY][IUZ]*YSSU+bxsmcd_[IUX][IUY][IUZ]*YSSL)/YUL;
  HYG=(bysmcd_[IUX][ILY][IUZ]*YSSU+bysmcd_[IUX][IUY][IUZ]*YSSL)/YUL;
  HZG=(bzsmcd_[IUX][ILY][IUZ]*YSSU+bzsmcd_[IUX][IUY][IUZ]*YSSL)/YUL;
  //
  //     interpolate in Y to find field components HF
  //     at (XDIPOL(IUX),Y,ZDIPOL(ILZ))
  //

  HXF=(bxsmcd_[IUX][ILY][ILZ]*YSSU+bxsmcd_[IUX][IUY][ILZ]*YSSL)/YUL;
  HYF=(bysmcd_[IUX][ILY][ILZ]*YSSU+bysmcd_[IUX][IUY][ILZ]*YSSL)/YUL;
  HZF=(bzsmcd_[IUX][ILY][ILZ]*YSSU+bzsmcd_[IUX][IUY][ILZ]*YSSL)/YUL;
  //
  //     The 6 face centers are AB,BC,CD,DA,HE,FG
  //
  //
  //     now interpolate HA and HD in Y to find 1st estimate of field, H1AD
  //
  HX1AD=(HXA*YSSU+HXD*YSSL)/YUL;
  HY1AD=(HYA*YSSU+HYD*YSSL)/YUL;
  HZ1AD=(HZA*YSSU+HZD*YSSL)/YUL;
  //
  //     now interpolate HC and HD in Y to find 2st estimate of field, H1BC
  //
  HX1BC=(HXB*YSSU+HXC*YSSL)/YUL;
  HY1BC=(HYB*YSSU+HYC*YSSL)/YUL;
  HZ1BC=(HZB*YSSU+HZC*YSSL)/YUL;
  //
  //     now interpolate HA and HB in Z to find 3st estimate of field, H1AB
  //
  HX1AB=(HXA*ZSSU+HXB*ZSSL)/ZUL;
  HY1AB=(HYA*ZSSU+HYB*ZSSL)/ZUL;
  HZ1AB=(HZA*ZSSU+HZB*ZSSL)/ZUL;
  //
  //     now interpolate HC and HD in Z to find 4st estimate of field, H1CD
  //
  HX1CD=(HXD*ZSSU+HXC*ZSSL)/ZUL;
  HY1CD=(HYD*ZSSU+HYC*ZSSL)/ZUL;
  HZ1CD=(HZD*ZSSU+HZC*ZSSL)/ZUL;
  //
  //     now interpolate HE and HH in Z to find 5st estimate of field, H1EH
  //
  HX1EH=(HXE*ZSSU+HXH*ZSSL)/ZUL;
  HY1EH=(HYE*ZSSU+HYH*ZSSL)/ZUL;
  HZ1EH=(HZE*ZSSU+HZH*ZSSL)/ZUL;
  //
  //     now interpolate HF and HG in Z to find 6st estimate of field, H1FG
  //
  HX1FG=(HXF*ZSSU+HXG*ZSSL)/ZUL;
  HY1FG=(HYF*ZSSU+HYG*ZSSL)/ZUL;
  HZ1FG=(HZF*ZSSU+HZG*ZSSL)/ZUL;
  //
  //     best estimate is mean of HAB,HBC,HCD,HDA,HEH,HFG
  //
  HX=(HX1AB+HX1BC+HX1CD+HX1AD+HX1EH+HX1FG)/6.0;
  HY=(HY1AB+HY1BC+HY1CD+HY1AD+HY1EH+HY1FG)/6.0;
  HZ=(HZ1AB+HZ1BC+HZ1CD+HZ1AD+HZ1EH+HZ1FG)/6.0;
  //
  //     now fill BXYZE ( field values )
  //
  BXYZE[0]=HX*COSRG1-HY*SINRG1;
  BXYZE[1]=HY*COSRG1+HX*SINRG1;
  BXYZE[2]=HZ;
  //
  //     Take into account the "right" symmetry :
  //
  //     b_x(-z)=-b_x(z)
  //     b_x(-x)=-b_x(x)
  //
  //     b_y(-z)=-b_y(z)
  //     b_x(-y)=-b_x(y)
  //
  if( IRSIGZ*IRSIGX <= 0 ) BXYZE[0] *= -1;

  if( IRSIGZ*IRSIGY <= 0 ) BXYZE[1] *= -1;


  field_x=BXYZE[1];  // Y -> X
  field_y=BXYZE[2];  // Z -> Y
  field_z=BXYZE[0];  // X -> Z

  return (true);
}



bool CsField::getFieldSMC( float pos_x, float pos_y, float pos_z,
			   float& bmx, float& bmy, float& bmz) {
  int JOUT;
  int ii,jj,ip,jp;

  float XFINT[2];
  float AFINT[4];
  float GFINT[2][2];

  float r,bmr,eps;

  r = sqrt( pos_x*pos_x + pos_y*pos_y );

  //--- Interpolation on a square grid (4 points)

  if( r >= rsmc_[NR_SMC-1] )
    jp=NR_SMC-1;
  else {
    for( jp=0; jp<NR_SMC; jp++ )
      if( rsmc_[jp] > r ) break;
  }
  jj = jp-1;

  if( fabs(pos_z) >= zsmc_[NZ_SMC-1] )
    ip=NZ_SMC-1;
  else {
    for( ip=0; ip<NZ_SMC; ip++ )
      if( zsmc_[ip] > fabs(pos_z) ) break;
  }
  ii = ip-1;

  JOUT = 0;
  if( ii > NZ_SMC || ii < 0 ) JOUT = 1;
  if( jj > NR_SMC || jj < 0 ) JOUT = 1;
  if( ip > NZ_SMC || ip < 0 ) JOUT = 1;
  if( jp > NR_SMC || jp < 0 ) JOUT = 1;
  if( JOUT != 0 ) {
    cout<<endl<<ii<<", "<<ip<<", "<<jj<<", "<<jp<<endl;
    CsErrLog::mes(elFatal,"getFieldSMC error.");
  }

  XFINT[0] = fabs(pos_z);
  XFINT[1] = r;

  AFINT[0] = zsmc_[ii];
  AFINT[1] = zsmc_[ip];
  AFINT[2] = rsmc_[jj];
  AFINT[3] = rsmc_[jp];

  GFINT[0][0] = bzsmc_[jj][ii];
  GFINT[1][0] = bzsmc_[jj][ip];
  GFINT[0][1] = bzsmc_[jp][ii];
  GFINT[1][1] = bzsmc_[jp][ip];

  bmz = fint2( XFINT, AFINT, GFINT);

  GFINT[0][0] = brsmc_[jj][ii];
  GFINT[1][0] = brsmc_[jj][ip];
  GFINT[0][1] = brsmc_[jp][ii];
  GFINT[1][1] = brsmc_[jp][ip];

  bmr = fint2( XFINT, AFINT, GFINT);
  // ---- decreasing --------------

  if( fabs(pos_z) > zsmc_[NZ_SMC-1] ){
    eps = exp( -(fabs(pos_z)-zsmc_[NZ_SMC-1])/60 );
    bmr = bmr * eps;
    bmz = bmz * eps;
  }
  if( r > rsmc_[NR_SMC-1] ){
    eps = exp( -(r-rsmc_[NR_SMC-1])/20 );
    bmr = bmr * eps;
    bmz = bmz * eps;
  }
  // -------------

  if( r > 0 ) {
    bmx = pos_x/r * bmr * sign(pos_z);
    bmy = pos_y/r * bmr * sign(pos_z);
  }else{
    bmx = 0;
    bmy = 0;
  }

  return (true);
}



bool CsField::getFieldSol( float pos_x, float pos_y, float pos_z,
			   float& bmx, float& bmy, float& bmz) {
  int ii,jj,ip,jp;

  float XFINT[2];
  float AFINT[4];
  float GFINT[2][2];

  float r,bmr,eps;

  r = sqrt( pos_x*pos_x + pos_y*pos_y );
  
  if( fabs(pos_z) < 300. )
    ii = (int)floor(fabs(pos_z)/5.);
  else
    ii = 59;
    
  if( r < 100. )
    jj = (int)floor(r/5.);
  else
    jj = 19;
    
  ip = ii + 1;
  jp = jj + 1;
  
#ifdef DEBUG
  JOUT = 0;
  if( ii > 60 || ii < 0 ) JOUT = 1;
  if( jj > 20 || jj < 0 ) JOUT = 1;
  if( ip > 60 || ip < 0 ) JOUT = 1;
  if( jp > 20 || jp < 0 ) JOUT = 1;
  if( JOUT != 0 ) {
    cout<<"getFieldSol error."<<endl;
    exit(0);
  }
#endif
  

  XFINT[0] = fabs(pos_z);
  XFINT[1] = r;

  AFINT[0] = zsol_[ii];
  AFINT[1] = zsol_[ip];
  AFINT[2] = rsol_[jj];
  AFINT[3] = rsol_[jp];

  GFINT[0][0] = bzsol_[ii][jj];
  GFINT[1][0] = bzsol_[ip][jj];
  GFINT[0][1] = bzsol_[ii][jp];
  GFINT[1][1] = bzsol_[ip][jp];

  bmz = fint2( XFINT, AFINT, GFINT);
  GFINT[0][0] = brsol_[ii][jj];
  GFINT[1][0] = brsol_[ip][jj];
  GFINT[0][1] = brsol_[ii][jp];
  GFINT[1][1] = brsol_[ip][jp];

  bmr = fint2( XFINT, AFINT, GFINT);
  // ---- decreasing --------------

  if( fabs(pos_z) > 300. ){
    eps = exp( -(fabs(pos_z)-300.)/60 );
    bmr = bmr * eps;
    bmz = bmz * eps;
  }
  if( r > 100 ){
    eps = exp( -(r-100.)/30 );
    bmr = bmr * eps;
    bmz = bmz * eps;
  }
  // -------------

  if( r > 0 ) {
    bmx = pos_x/r * bmr * sign(pos_z);
    bmy = pos_y/r * bmr * sign(pos_z);
  }else{
    bmx = 0;
    bmy = 0;
  }

  return (true);
}

bool CsField::getFieldSm1( float pos_x, float pos_y, float pos_z,
			   float& bmx, float& bmy, float& bmz) {
  int ii,jj,kk,ip,jp,kp;
  float eps;

  float XFINT[3];
  float AFINT[6];
  float FFINT[2][2][2];

  if( fabs(pos_x) < 200 && fabs(pos_y) < 300 && pos_z < 800 ){
    ii = (int)floor(fabs(pos_x)/10);
    jj = (int)floor(fabs(pos_y)/10);
    if( pos_z < 0 ){
      kk = 0;
    } else {
      kk = (int)floor(fabs(pos_z)/10);
    }
    if( pos_z < 0 ) kk = 0;

    //--- Interpolation on a cube (8 points)
    ip = ii + 1;
    jp = jj + 1;
    kp = kk + 1;

    XFINT[0] = fabs(pos_x);
    XFINT[1] = fabs(pos_y);
    if( pos_z < 0 ) XFINT[2] = 0;
    else            XFINT[2] = pos_z;

    AFINT[0] = xsm1_[ii];
    AFINT[1] = xsm1_[ip];
    AFINT[2] = ysm1_[jj];
    AFINT[3] = ysm1_[jp];
    AFINT[4] = zsm1_[kk];
    AFINT[5] = zsm1_[kp];

    FFINT[0][0][0] = bxsm1_[ii][jj][kk];
    FFINT[1][0][0] = bxsm1_[ip][jj][kk];
    FFINT[0][1][0] = bxsm1_[ii][jp][kk];
    FFINT[1][1][0] = bxsm1_[ip][jp][kk];
    FFINT[0][0][1] = bxsm1_[ii][jj][kp];
    FFINT[1][0][1] = bxsm1_[ip][jj][kp];
    FFINT[0][1][1] = bxsm1_[ii][jp][kp];
    FFINT[1][1][1] = bxsm1_[ip][jp][kp];
    bmx = fint3( XFINT, AFINT, FFINT);
    FFINT[0][0][0] = bysm1_[ii][jj][kk];
    FFINT[1][0][0] = bysm1_[ip][jj][kk];
    FFINT[0][1][0] = bysm1_[ii][jp][kk];
    FFINT[1][1][0] = bysm1_[ip][jp][kk];
    FFINT[0][0][1] = bysm1_[ii][jj][kp];
    FFINT[1][0][1] = bysm1_[ip][jj][kp];
    FFINT[0][1][1] = bysm1_[ii][jp][kp];
    FFINT[1][1][1] = bysm1_[ip][jp][kp];
    bmy = fint3( XFINT, AFINT, FFINT);
    FFINT[0][0][0] = bzsm1_[ii][jj][kk];
    FFINT[1][0][0] = bzsm1_[ip][jj][kk];
    FFINT[0][1][0] = bzsm1_[ii][jp][kk];
    FFINT[1][1][0] = bzsm1_[ip][jp][kk];
    FFINT[0][0][1] = bzsm1_[ii][jj][kp];
    FFINT[1][0][1] = bzsm1_[ip][jj][kp];
    FFINT[0][1][1] = bzsm1_[ii][jp][kp];
    FFINT[1][1][1] = bzsm1_[ip][jp][kp];
    bmz = fint3( XFINT, AFINT, FFINT);
  } else {
    bmx = 0.;
    bmy = 0.;
    bmz = 0.;
  }
  //--  decreasing ----------

  if( pos_z < 0 ){
    eps = exp( -pos_z * pos_z/2500. );
    bmx = bmx * eps;
    bmy = bmy * eps;
    bmz = bmz * eps;
  }
  //---------------------------

  bmx = bmx * sign(pos_y) * sign(pos_x);
  bmz = bmz * sign(pos_y);

  return (true);
}



bool CsField::getFieldSM1( float pos_x, float pos_y, float pos_z,
			   float& bmx, float& bmy, float& bmz) {
  int ii,jj,kk,ip,jp,kp;
  float eps;

  float XFINT[3];
  float AFINT[6];
  float FFINT[2][2][2];

  if( fabs(pos_x) < 160 && fabs(pos_y) < 160 &&
      pos_z > -416 && pos_z < 440){

    ii = (int)floor(fabs(pos_x)/4.);
    jj = (int)floor(pos_y/8.+20);
    kk = (int)floor(pos_z/8.+52);

    //--- Interpolation on a cube (8 points)
    ip = ii + 1;
    jp = jj + 1;
    kp = kk + 1;

    XFINT[0] = fabs(pos_x);
    XFINT[1] = pos_y;
    XFINT[2] = pos_z;

    AFINT[0] = xSM1_[ii];
    AFINT[1] = xSM1_[ip];
    AFINT[2] = ySM1_[jj];
    AFINT[3] = ySM1_[jp];
    AFINT[4] = zSM1_[kk];
    AFINT[5] = zSM1_[kp];

    FFINT[0][0][0] = bxSM1_[ii][jj][kk];
    FFINT[1][0][0] = bxSM1_[ip][jj][kk];
    FFINT[0][1][0] = bxSM1_[ii][jp][kk];
    FFINT[1][1][0] = bxSM1_[ip][jp][kk];
    FFINT[0][0][1] = bxSM1_[ii][jj][kp];
    FFINT[1][0][1] = bxSM1_[ip][jj][kp];
    FFINT[0][1][1] = bxSM1_[ii][jp][kp];
    FFINT[1][1][1] = bxSM1_[ip][jp][kp];
    bmx = fint3( XFINT, AFINT, FFINT);
    FFINT[0][0][0] = bySM1_[ii][jj][kk];
    FFINT[1][0][0] = bySM1_[ip][jj][kk];
    FFINT[0][1][0] = bySM1_[ii][jp][kk];
    FFINT[1][1][0] = bySM1_[ip][jp][kk];
    FFINT[0][0][1] = bySM1_[ii][jj][kp];
    FFINT[1][0][1] = bySM1_[ip][jj][kp];
    FFINT[0][1][1] = bySM1_[ii][jp][kp];
    FFINT[1][1][1] = bySM1_[ip][jp][kp];
    bmy = fint3( XFINT, AFINT, FFINT);
    FFINT[0][0][0] = bzSM1_[ii][jj][kk];
    FFINT[1][0][0] = bzSM1_[ip][jj][kk];
    FFINT[0][1][0] = bzSM1_[ii][jp][kk];
    FFINT[1][1][0] = bzSM1_[ip][jp][kk];
    FFINT[0][0][1] = bzSM1_[ii][jj][kp];
    FFINT[1][0][1] = bzSM1_[ip][jj][kp];
    FFINT[0][1][1] = bzSM1_[ii][jp][kp];
    FFINT[1][1][1] = bzSM1_[ip][jp][kp];
    bmz = fint3( XFINT, AFINT, FFINT);
  } else {
    bmx = 0.;
    bmy = 0.;
    bmz = 0.;
  }
  //--  decreasing ----------
  /*
  if( pos_z < 0 ){
    eps = exp( -pos_z * pos_z/2500. );
    bmx = bmx * eps;
    bmy = bmy * eps;
    bmz = bmz * eps;
  }
  */
  //---------------------------

  //bmx = bmx * sign(pos_y) * sign(pos_x);
  //bmz = bmz * sign(pos_y);

  bmx = bmx  * sign(pos_x);

  return (true);
}


bool CsField::getFieldDipOx( float pos_x, float pos_y, float pos_z,
			     float& bmx,  float& bmy,  float& bmz ) {
  //printf("%f %f %f  %f %f %f\n", pos_x,pos_y,pos_z,  bmx,bmy,bmz );

  if( fabs(pos_x) < xoxd_[NX_OXD-1] &&
      fabs(pos_y) < yoxd_[NY_OXD-1] &&
      fabs(pos_z) < zoxd_[NZ_OXD-1] ){
    
    int ii=0;
    while ( fabs(pos_x)>=xoxd_[ii] && ii<NX_OXD ) ii++;
    ii=ii-1;
    
    int jj=0;
    while ( fabs(pos_y)>=yoxd_[jj] && jj<NY_OXD ) jj++;
    jj=jj-1;
    
    int kk=0;
    while ( fabs(pos_z)>=zoxd_[kk] && kk<NZ_OXD ) kk++;
    kk=kk-1;
    
    //printf("%i %i %i  %f %f %e\n", ii,jj,kk, pos_x, xoxd_[0], pos_x-xoxd_[0] );
    //printf("%i %i %i  %f %f %f\n", ii,jj,kk, pos_x, pos_y, pos_z );
    
    if(ii<0 || jj<0 || kk<0) {
      char cmess[100];
      sprintf(cmess,"CsField::getFieldDipOx==> %i %i %i  %f %f %f\n", ii,jj,kk, pos_x,pos_y,pos_z);
      CsErrLog::mes(elFatal,cmess);
    }
    
    //--- Interpolation on a cube (8 points)
    int ip = ii + 1;
    int jp = jj + 1;
    int kp = kk + 1;

    float XFINT[3];
    XFINT[0] = fabs(pos_x);
    XFINT[1] = fabs(pos_y);
    XFINT[2] = fabs(pos_z);
    //if( pos_z < 0 ) XFINT[2] = 0;
    //else            XFINT[2] = pos_z;
    
    float AFINT[6];
    AFINT[0] = xoxd_[ii];
    AFINT[1] = xoxd_[ip];
    AFINT[2] = yoxd_[jj];
    AFINT[3] = yoxd_[jp];
    AFINT[4] = zoxd_[kk];
    AFINT[5] = zoxd_[kp];
    
    float FFINT[2][2][2];
    FFINT[0][0][0] = bxoxd_[ii][jj][kk];
    FFINT[1][0][0] = bxoxd_[ip][jj][kk];
    FFINT[0][1][0] = bxoxd_[ii][jp][kk];
    FFINT[1][1][0] = bxoxd_[ip][jp][kk];
    FFINT[0][0][1] = bxoxd_[ii][jj][kp];
    FFINT[1][0][1] = bxoxd_[ip][jj][kp];
    FFINT[0][1][1] = bxoxd_[ii][jp][kp];
    FFINT[1][1][1] = bxoxd_[ip][jp][kp];
    bmx = fint3( XFINT, AFINT, FFINT);
    FFINT[0][0][0] = byoxd_[ii][jj][kk];
    FFINT[1][0][0] = byoxd_[ip][jj][kk];
    FFINT[0][1][0] = byoxd_[ii][jp][kk];
    FFINT[1][1][0] = byoxd_[ip][jp][kk];
    FFINT[0][0][1] = byoxd_[ii][jj][kp];
    FFINT[1][0][1] = byoxd_[ip][jj][kp];
    FFINT[0][1][1] = byoxd_[ii][jp][kp];
    FFINT[1][1][1] = byoxd_[ip][jp][kp];
    bmy = fint3( XFINT, AFINT, FFINT);
    FFINT[0][0][0] = bzoxd_[ii][jj][kk];
    FFINT[1][0][0] = bzoxd_[ip][jj][kk];
    FFINT[0][1][0] = bzoxd_[ii][jp][kk];
    FFINT[1][1][0] = bzoxd_[ip][jp][kk];
    FFINT[0][0][1] = bzoxd_[ii][jj][kp];
    FFINT[1][0][1] = bzoxd_[ip][jj][kp];
    FFINT[0][1][1] = bzoxd_[ii][jp][kp];
    FFINT[1][1][1] = bzoxd_[ip][jp][kp];
    bmz = fint3( XFINT, AFINT, FFINT);
  } else {
    bmx = 0.;
    bmy = 0.;
    bmz = 0.;
  }
  
  //---------------------------

  bmx *= ( sign(pos_y) * sign(pos_x) );
  bmz *= ( sign(pos_y) * sign(pos_z) );
  
  return (true);
}
