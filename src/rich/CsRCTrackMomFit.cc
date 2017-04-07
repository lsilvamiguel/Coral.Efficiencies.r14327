/*!
   \file    CsRCTrackMomFit.cc
   \------------------------
   \brief   CsRCTrackMomFit class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    June 2005
*/



#include "CsGeom.h"
#include "CsField.h"

//----------------------------
#include "CsRCTrackMomFit.h"

//----------------------------

using namespace CLHEP;

  CsRCTrackMomFit* CsRCTrackMomFit::instance_ = 0;


//==========================================================================
  CsRCTrackMomFit* CsRCTrackMomFit::Instance() {
//----------------------------------------------
    if( instance_ == 0 ) instance_ = new CsRCTrackMomFit();
    return instance_;
  }


//==========================================================================
  CsRCTrackMomFit::CsRCTrackMomFit() {}
//-------------------------------------


//==========================================================================
  bool CsRCTrackMomFit::test( Hep3Vector& vPoCoo, Hep3Vector& vDcCoo,
//-------------------------------------------------------------------
			      double chgOvMom, double zExit ) {

    //Hep3Vector vPoCoo;
    //vPoCoo.setX( 100. );
    //vPoCoo.setY( 0. );
    //vPoCoo.setZ( 1000. );
    //Hep3Vector vDcCoo;
    //vDcCoo.setX( 0.000 );
    //vDcCoo.setY( 0.000 );
    //vDcCoo.setZ( 1.000 );
    //double chgOvMom = 0.0003/10.;
    float stepIn = 0.;
    int nCall = -1;
    int nCoInt = 3;
    std::vector<Hep3Vector> vPoInt;
    std::vector<Hep3Vector> vDcInt;
    Hep3Vector vPo( 0.), vDc( 0.);
    if( nCoInt > 0 ) {
      for(int ki=0;ki<nCoInt; ki++ ) {
	vPoInt.push_back( vPo );
	vDcInt.push_back( vDc );
      }
      vPoInt[0].setZ( 3000.);
      vPoInt[1].setZ( 3500.);
      vPoInt[2].setZ( 4000.);
    }
    Hep3Vector vPoExit( 0.);
    vPoExit.setZ( zExit );
    Hep3Vector vDcExit( 0.);

    //if( CsRCTrackThruField( vPoCoo, vDcCoo, chgOvMom,
    //		    stepIn, nCall, nCoInt,
    //		    vPoInt, vDcInt, vPoExit, vDcExit ) ) {
    if( CsRCTrackThruField( vPoCoo, vDcCoo, chgOvMom,
	                    stepIn,
			    vPoExit, vDcExit ) ) {
      vPoCoo = vPoExit;
      vDcCoo = vDcExit;
      return  true;
    }
    else  return  false;

  }


//==========================================================================
  bool CsRCTrackMomFit::CsRCTrackThruField(
//---------------------------------------------------
		Hep3Vector vPoCoo, Hep3Vector vDcCoo,
		double chgOvMom,
		float& stepIn, int& nCall, int nCoInt,
		std::vector<Hep3Vector>& vPoInt,
		std::vector<Hep3Vector>& vDcInt,
		Hep3Vector& vPoExit, Hep3Vector& vDcExit ) {

/*!
   \file    CsRCTrackThruField.cc
   \-----------------------------
   \brief   CsRCTrackThruField class author.
   \implementation  Paolo Schiavon
   \version 0.01
   \date    June 2005
*/


//---------  version 1.01 (ex 3 a).
//---------  date  22/06/05 (ex 8/98).
//
//---------  from  DMGTRACK = MGTRACK2 and DMGTRAFD
//---------  coord.s in mm, momentum in GeV/c.
//---------  upgraded step control.
//---------  corrected (1/2) formula.
//---------  check for backward loop.
//---------  main motion in the z-axis direction.
//---------  BdL calculated along z-axis. (NOT YET)!
//---------  chambers inside the mg-field.
//---------  corrected allfd(.),ammfd(.) (3/11/97).
//---------  normalized direction cos.s after each pass (1/12/97).
//---------  no need of improve ssr (1/12/97).
//---------  backward tracking added - 8/98 :
//           use appropriate direction cos.s at in;
//           for backward tracking of forward going particle,
//           change sign to ccp.
//---------  direction cosines at in and out - 9/98.
//---------  improved double-step tracking added - 9/98.
//---------  stepIn at input = 0. : defaults
//           stepIn at input > 0. : stepx = abs(stepIn), step = stepx
//           stepIn at input < 0. : stepm = abs(stepIn), step = stepm
//---------  last revision on 14/9/98.

//- assume the main field component along y-axis
//  and the main motion along z-axis !
//------------------------------------

    static bool print = true;
    static bool tprint = false;
    static int nCallMx = 5000;
//@@--------------------------

    CsField* field = CsGeom::Instance()->getCsField();

    bool errFlag = true;
    static std::string mess = "CsRCTrackMomFit::CsRCTrackThruField : ";


//- set field strenght for step control. (prov???)
//--------------------------------------
    float bbXfrac = 0.01;
    float bbYfrac = 0.05;
    float bbZfrac = 0.01;
    float bbsx, bbsy, bbsz;
    float bbssx, bbssy, bbssz;
    double bss = 0.;
    double bxdl = 0.;
    double bydl = 0.;
    double bzdl = 0.;

//- input quantities.
//-------------------
    Hep3Vector vPo = vPoCoo;
    Hep3Vector vDc = vDcCoo;
    vDc = vDc.unit();
    double zzpin = vPoCoo.z();
    int sigXin = 1;
    if( vDc.x() < 0.) sigXin = -1;
    int sigYin = 1;
    if( vDc.y() < 0.) sigYin = -1;
    int sigZin = 1;
    if( vDc.z() < 0.) sigZin = -1;
    int sigZv = sigZin;
    int kVers = sigZv;
    double zzpout = vPoExit.z();
    double ccp = chgOvMom;
    if( ( kVers == 1 && zzpout < zzpin ) || ( kVers == -1 && zzpin < zzpout ) ) {
      std::string err = "error in input conditions (zzin><zzout)";
      mess.append( err );
      errFlag = false;
      return  errFlag;
    }

//- conditional printing (input).
//-------------------------------
    if( print ) {
      std::cout << std::endl;
      std::cout << "CsRCTrackThruField input : " << std::endl;
      std::cout << "xxin, yyin, zzin = "
		<< vPoCoo.x() << "  " << vPoCoo.y() << "  "
		<< vPoCoo.z() << std::endl;
      std::cout << "dxxin, dyyin, dzzin = "
		<< vDcCoo.x() << "  " << vDcCoo.y() << "  "
		<< vDcCoo.z() << std::endl;
      std::cout << "zzInt = ";
      for(int k=0; k<nCoInt; k++) std::cout << vPoInt[k].z() << "  ";
      std:: cout << std::endl;
      std:: cout << "ccp, ssl = " << ccp << "  " << stepIn << std::endl;
      std::cout << "zzpout = " << zzpout << std::endl;
      //print *,'bbymax, bbyfra = ',bbymax,bbyfra
    }

//- check intermediate chambers layout.
//-------------------------------------
    int kCoInt = 0;                                     //   ???
    for( int kC=0; kC<nCoInt; kC++ ) {
      int kCk = kC;
      if( kVers == -1 ) kCk = nCoInt-1 - kC;            //   ???
      float zzpl = vPoInt[kCk].z();
      vPoInt[kCk].setX( 0.);
      vPoInt[kCk].setY( 0.);
      //- WARNING - tangents!
      vDcInt[kCk].setX( 0.);
      vDcInt[kCk].setY( 0.);
      vDcInt[kCk].setZ( 1.);
      //float zzpin = vPoCoo.z();                       //#
      if( kVers * (zzpl-zzpin) < 0. ) kCoInt++;
    }

//- set variable step limits.
//---------------------------
    static const float stepMn =  4.;          //      minimun step value
    static const float stepMx = 16.;          //      maximun step value
//@@-------------------------------
    float stepm = stepMn;
    float stepx = stepMx;
    float step = stepx / 2.;
    if( step <= 0.) {
      std::string err = "error in input conditions (step=0)";
      mess.append( err );
      errFlag = false;
      return  errFlag;
    }
    if( stepIn < 0. ) {
      stepm = fabs( stepIn );
      if( stepx < stepm ) stepx = stepm;
      step = stepm;
    }
    if( stepIn > 0. ) {
      stepx = fabs( stepIn );
      //if( stepm > stepx ) stepm = stepx;      //   new!
      stepm = stepx;;                           //   new!
      step = stepx;
    }
    int nCallMn = 2* int( (zzpout-zzpin) / step );
    if( nCallMn > nCallMx ) nCallMx = 2* nCallMn;
//@@--------------------------------------------


//- start step-by-step tracking.
//------------------------------
    int kStep = 0;
    int kCall = 0;
    double trajLen = 0.;
    bool bDivide = false;

    double xxs = vPo.x();
    double yys = vPo.y();
    double zzs = vPo.z();
    double dxxs = vDc.x();
    double dyys = vDc.y();
    double dzzs = vDc.z();
    double xxss = 0., yyss = 0., zzss = 0.;
    double dxxss = 0., dyyss = 0., dzzss = 0.;


    while ( kCall < nCallMx ) {

      //bool bDivide = false;                         //   ???

//--- call to field at the beginning of the step.
//-----------------------------------------------
      if( !field->getField( float(xxs), float(yys), float(zzs),
//-------------------------------------------------------------
			    bbsx, bbsy, bbsz ) ) {
	std::string err = "error in getField call";
	mess.append( err );
	errFlag = false;
	return  errFlag;
      }
      kCall++;

//--- determine position after half step.
//---------------------------------------
      double step2 = step / 2.;
      double xxr = dxxs * step2;
      double yyr = dyys * step2;
      double zzr = dzzs * step2;
      double dx2ss = ccp * (bbsz*yyr - bbsy*zzr);
      double dy2ss = ccp * (bbsx*zzr - bbsz*xxr);
      double dz2ss = ccp * (bbsy*xxr - bbsx*yyr);
      xxss = xxs + xxr + 0.5 * dx2ss * step2;
      yyss = yys + yyr + 0.5 * dy2ss * step2;
      zzss = zzs + zzr + 0.5 * dz2ss * step2;

//--- call to field in the middle of the step.
//--------------------------------------------
      if( !field->getField( float(xxss), float(yyss), float(zzss),
//----------------------------------------------------------------
			    bbssx, bbssy, bbssz ) ) {
	std::string err = "error in getField call";
	mess.append( err );
	errFlag = false;
	return  errFlag;
      }
      kCall++;

//--- determine step versus field variation.
//------------------------------------------
      if( fabs( bbssx-bbsx ) > bbXfrac ) bDivide = true;
      if( fabs( bbssy-bbsy ) > bbYfrac ) bDivide = true;
      if( fabs( bbssz-bbsz ) > bbZfrac ) bDivide = true;
      if( bDivide ) {
	if( fabs( step ) > stepm ) {
	  step /= 2.;
	  if( fabs( step ) < stepm ) {
	    if( step > 0.) step = stepm;
	    else step = - stepm;
	  }
	  continue;                       //   ???
	}
      }

//--- determine position and direction after full step.
//-----------------------------------------------------
      bbsx = bbssx;
      bbsy = bbssy;
      bbsz = bbssz;
      xxr = dxxs * step;
      yyr = dyys * step;
      zzr = dzzs * step;
      dx2ss = ccp * (bbsz*yyr - bbsy*zzr);
      dy2ss = ccp * (bbsx*zzr - bbsz*xxr);
      dz2ss = ccp * (bbsy*xxr - bbsx*yyr);
      dxxss = dxxs + dx2ss;
      dyyss = dyys + dy2ss;
      dzzss = dzzs + dz2ss;
      xxss = xxs + xxr + 0.5 * dx2ss * step;
      yyss = yys + yyr + 0.5 * dy2ss * step;
      zzss = zzs + zzr + 0.5 * dz2ss * step;

//--- check tracking inversion (vs z-axis).
//-----------------------------------------
      int sigZ = 1;
      if( dzzss < 0. ) sigZ = -1;
//??      kVers = sigZv;
//$$      if(ksgz.eq.-ksgzin) go to 620          !   no back-going particle
//$$      if(kver*(zzss-zzpin).le.0) go to 620   !   no inversion,   cbk
      if( sigZ == - sigZv ) {
	vPoExit.setZ( vPoCoo.z() );
	errFlag = false;
	std::string err = "particle tracking inversion";
	mess.append( err );
	break;
      }

//--- determine position and direction on intermediate chambers.
//--------------------------------------------------------------
      bool bCoInt = true;
      if( kCoInt < nCoInt ) {                           //   ???
	int kCk = kCoInt;
	if( kVers == -1) kCk = nCoInt-1 - kCoInt;       //   ???
	float zzpl = vPoInt[kCk].z();
	if( kVers * (zzss-zzpl) >= 0 ) {
	  float ffr = (zzpl-zzs) / (zzss-zzs);
	  xxss = xxs + ffr*(xxss-xxs);
	  yyss = yys + ffr*(yyss-yys);
	  zzss = zzpl;
	  dxxss = dxxs + ffr*(dxxss-dxxs);
	  dyyss = dyys + ffr*(dyyss-dyys);
	  dzzss = dzzs + ffr*(dzzss-dzzs);
	  vPoInt[kCk].setX( xxss );
	  vPoInt[kCk].setY( yyss );
	  //- WARNING - tangents!
	  vDcInt[kCk].setX( dxxss/dzzss );
	  vDcInt[kCk].setY( dyyss/dzzss );
	  vDcInt[kCk].setZ( 1. );
          trajLen += fabs( ffr*step );      //      trajectory length
	  //bss = fabs(ffr*ss);              //      bdl along the trajectory
	  bss = fabs(zzpl-zzs);              //      bdl along z-axis
	  bxdl = bxdl + bbsx*bss;
	  bydl = bydl + bbsy*bss;
	  bzdl = bzdl + bbsz*bss;
	  kStep++;
//------- conditional printing (int pl).
	  if( print ) {
	    std::cout << std::endl;
	    std::cout << "CsRCTrackThruField int pl (" << zzpl
		      << ") : " << std::endl;
            std::cout << "kStep, kCall, bDivide, step, kVers, trajLen = "
	    	      << kStep << "  " << kCall << "  " << bDivide << "  "
		      << step << "  " << kVers << " " << trajLen << std::endl;
            std::cout << "xxss, yyss, zzss = " << xxss << "  " << yyss << "  "
	  	      << zzs << std::endl;
            std::cout << "dxxss, dyyss, dzzss = " << dxxs << "  " << dyys
		      << "  " << dzzss << std::endl;
            std::cout << "Bx, By, Bz = " << bbssx << "  " << bbssy << "  "
		      << bbssz << std::endl;
	  }
	  kCoInt++;
	  bCoInt = false;
	}
      }

//--- check end-of-traking and
//--- determine position and direction on exit.
//---------------------------------------------
      if( kVers * (zzss-zzpout) >= 0 ) {
	if( !bCoInt ) break;                    // ???
	double ffr = (zzpout-zzs) / (zzss-zzs);
	xxss = xxs + ffr*(xxss-xxs);
	yyss = yys + ffr*(yyss-yys);
	zzss = zzpout;
	dxxss = dxxs + ffr*(dxxss-dxxs);
	dyyss = dyys + ffr*(dyyss-dyys);
	dzzss = dzzs + ffr*(dzzss-dzzs);
	double sqdxyz = sqrt( dxxss*dxxss + dyyss*dyyss + dzzss*dzzss );
	dxxss = dxxss / sqdxyz;
	dyyss = dyyss / sqdxyz;
	dzzss = dzzss / sqdxyz;
        trajLen += fabs( ffr*step );      //      trajectory length
        //bss = fabs(ffr*ss);              //      bdl along the trajectory
        bss = fabs(zzpout-zzs);            //      bdl along z-axis
        bxdl = bxdl + bbsx*bss;
        bydl = bydl + bbsy*bss;
        bzdl = bzdl + bbsz*bss;
        kStep++;
	break;
      }

//--- prepare for the next step.
//------------------------------
      double zzsv = zzs;
      dxxs = dxxss;
      dyys = dyyss;
      dzzs = dzzss;
      double sqdxyz = sqrt( dxxs*dxxs + dyys*dyys + dzzs*dzzs );
      dxxs = dxxs / sqdxyz;
      dyys = dyys / sqdxyz;
      dzzs = dzzs / sqdxyz;
      xxs = xxss;
      yys = yyss;
      zzs = zzss;
      sigZv = sigZ;

      if( bCoInt ) {
        trajLen += fabs( step );      //      trajectory length
        //bss = fabs(step);            //      bdl along the trajectory
        bss = fabs(zzss-zzsv);         //      bdl along z-axis
        bxdl = bxdl + bbsx*bss;
	bydl = bydl + bbsy*bss;
        bzdl = bzdl + bbsz*bss;
        kStep++;
      }

//--- conditional printing (tracking).
//------------------------------------
      if( tprint ) {
	std::cout << std::endl;
	std::cout << "CsRCTrackThruField tracking : " << std::endl;
	if( !bCoInt ) std::cout << "------ on int. plane ------" << std::endl;
        std::cout << "kStep, kCall, bDivide, step, kVers, trajLen = "
	  	  << kStep << "  " << kCall << "  " << bDivide << "  "
		  << step << "  " << kVers << "  " << trajLen << std::endl;
        std::cout << "xxs, yys, zzs = " << xxs << "  " << yys << "  "
	  	  << zzs << std::endl;
        std::cout << "dxxs, dyys, dzzs = " << dxxs << "  " << dyys << "  "
	  	  << dzzs << std::endl;
        std::cout << "Bx, By, Bz = " << bbsx << "  " << bbsy << "  "
	          << bbsz << std::endl;
      }

//--- go to the next step.
//------------------------
      if( !bCoInt) continue;                  //   ???
      if( fabs( step ) >= stepx) continue;
      if( bDivide ) continue;
      step *= 2.0;
      if( fabs( step ) > stepx) {
	if( step > 0.) step = stepx;
	else step = - stepx;
      }

    }
    if( kCall >= nCallMx ) {
      std::string err = "too many calls to mag field";
      mess.append( err );
      errFlag = false;
    }

//- output quantities.
//--------------------
    vPoExit.setX( xxss );
    vPoExit.setY( yyss );
    vDcExit.setX( dxxss );
    vDcExit.setY( dyyss );
    vDcExit.setZ( dzzss );
    stepIn = trajLen;

    if( !errFlag ) std::cout << std::endl << mess << std::endl;
//------------------------------------------------------------

//- conditional printing (output).
//--------------------------------
    if( print ) {
      std::cout << std::endl;
      std::cout << "CsRCTrackThruField output : " << std::endl;
      std::cout << "kStep, nCall, trajLen = " << kStep << "  "
		<< kCall << "  " << trajLen << std::endl;
      std::cout << "xout, yout, zout = " << vPoExit.x() << "  "
		<< vPoExit.y() << "  " << vPoExit.z() << std::endl;
      std::cout << "dxout, dyout, dzout = " << vDcExit.x() << "  "
		<< vDcExit.y() << "  " << vDcExit.z() << std::endl;
      std::cout << "bxdl, bydl, bzdl (T/mm) = " << bxdl << "  "
		<< bydl << "  " << bzdl << std::endl; 
    }


    return  errFlag;

  }


//==========================================================================
  bool CsRCTrackMomFit::CsRCTrackThruField(
//---------------------------------------------------
		Hep3Vector vPoCoo, Hep3Vector vDcCoo,
		double chgOvMom, float& stepIn,
		Hep3Vector& vPoExit, Hep3Vector& vDcExit ) {

/*!
   \file    CsRCTrackThruField.cc
   \-----------------------------
   \brief   CsRCTrackThruField class author.
   \implementation  Paolo Schiavon
   \version 0.01
   \date    June 2005
*/


//---------  version 1.01 (ex 3 a).
//---------  date  22/06/05 (ex 8/98).
//
//---------  from  DMGTRACK = MGTRACK2 and DMGTRAFD
//---------  coord.s in mm, momentum in GeV/c.
//---------  upgraded step control.
//---------  corrected (1/2) formula.
//---------  check for backward loop.
//---------  main motion in the z-axis direction.
//---------  normalized direction cos.s after each pass (1/12/97).
//---------  no need of improve ssr (1/12/97).
//---------  backward tracking added - 8/98 :
//           use appropriate direction cos.s at in;
//           for backward tracking of forward going particle,
//           change sign to ccp.
//---------  direction cosines at in and out - 9/98.
//---------  improved double-step tracking added - 9/98.
//---------  stepIn at input = 0. : defaults
//           stepIn at input > 0. : stepx = abs(stepIn), step = stepx
//           stepIn at input < 0. : stepm = abs(stepIn), step = stepm
//---------  last revision on 14/9/98.

//- assume the main field component along y-axis
//  and the main motion along z-axis !
//------------------------------------

    static bool print = false;
    static bool tprint = false;

    static int nCallMx = 5000;
//@@--------------------------

    CsField* field = CsGeom::Instance()->getCsField();

    bool errFlag = true;
    static std::string mess = "CsRCTrackMomFit::CsRCTrackThruField : ";


//- set field strenght for step control. (prov???)
//--------------------------------------
    Hep3Vector BBvar( 0.01, 0.05, 0.01 );
//@@------------------------------------
    Hep3Vector BBs( 0., 0., 0. );
    Hep3Vector BBss( 0., 0., 0. );
    float bbsx, bbsy, bbsz;
    float bbssx, bbssy, bbssz;

//- input quantities.
//-------------------
    Hep3Vector vPos = vPoCoo;
    Hep3Vector vDcs = vDcCoo;
    vDcs = vDcs.unit();
    double zzpin = vPoCoo.z();
    int sigXin = 1;
    if( vDcs.x() < 0.) sigXin = -1;
    int sigYin = 1;
    if( vDcs.y() < 0.) sigYin = -1;
    int sigZin = 1;
    if( vDcs.z() < 0.) sigZin = -1;
    int sigZv = sigZin;
    int kVers = sigZv;
    double zzpout = vPoExit.z();
    double ccp = chgOvMom;
    if( ( kVers == 1 && zzpout < zzpin ) || ( kVers == -1 && zzpin < zzpout ) ) {
      std::string err = "error in input conditions (zzin><zzout)";
      mess.append( err );
      errFlag = false;
      return  errFlag;
    }

//- conditional printing (input).
//-------------------------------
    if( print ) {
      std::cout << std::endl;
      std::cout << "CsRCTrackThruField input : " << std::endl;
      std::cout << "xxin, yyin, zzin = "
		<< vPoCoo.x() << "  " << vPoCoo.y() << "  "
		<< vPoCoo.z() << std::endl;
      std::cout << "dxxin, dyyin, dzzin = "
		<< vDcCoo.x() << "  " << vDcCoo.y() << "  "
		<< vDcCoo.z() << std::endl;
      std::cout << "dxxin, dyyin, dzzin = "
		<< vDcs.x() << "  " << vDcs.y() << "  "
		<< vDcs.z() << std::endl;
      std::cout << "ccp, ssl = " << ccp << "  " << stepIn << std::endl;
      std::cout << "zzpout = " << zzpout << std::endl;
      std::cout << "BBvar = " << BBvar.x() << "  " << BBvar.y()
		<< "  " << BBvar.z() << std::endl;
    }

//- set variable step limits.
//---------------------------
    static const float stepMn =  4.;          //      minimun step value
    static const float stepMx = 16.;          //      maximun step value
//@@-------------------------------
    float stepm = stepMn;
    float stepx = stepMx;
    float step = stepx / 2.;
    if( step <= 0.) {
      std::string err = "error in input conditions (step=0)";
      mess.append( err );
      errFlag = false;
      return  errFlag;
    }
    if( stepIn < 0. ) {
      stepm = fabs( stepIn );
      if( stepx < stepm ) stepx = stepm;
      step = stepm;
    }
    if( stepIn > 0. ) {
      stepx = fabs( stepIn );
      //if( stepm > stepx ) stepm = stepx;      //   mod!
      stepm = stepx;;                           //   mod!
      step = stepx;
    }
    int nCallMn = 2* int( (zzpout-zzpin) / step );
    if( nCallMn > nCallMx ) nCallMx = 2* nCallMn;
//@@--------------------------------------------


//- start step-by-step tracking.
//------------------------------
    int kStep = 0;
    int kCall = 0;
    Hep3Vector vPoss( 0., 0., 0. );
    Hep3Vector vDcss( 0., 0., 0. );


    while ( kCall < nCallMx ) {

      bool bDivide = false;

//--- call to field at the beginning of the step.
//-----------------------------------------------
      if( !field->getField( vPos.x(), vPos.y(), vPos.z(),
//-------------------------------------------------------
			    bbsx, bbsy, bbsz ) ) {
	std::string err = "error in getField call";
	mess.append( err );
	errFlag = false;
	return  errFlag;
      }
      BBs.set( bbsx, bbsy, bbsz );
      kCall++;

//--- determine position after half step.
//---------------------------------------
      double step2 = step / 2.;
      Hep3Vector vDcr = step2 * vDcs;
      Hep3Vector dv2ss = ccp * vDcr.cross( BBs );
      vPoss = vPos + vDcr + 0.5 * step2 * dv2ss;

//--- call to field in the middle of the step.
//--------------------------------------------
      if( !field->getField( vPoss.x(), vPoss.y(), vPoss.z(),
//----------------------------------------------------------
			    bbssx, bbssy, bbssz ) ) {
	std::string err = "error in getField call";
	mess.append( err );
	errFlag = false;
	return  errFlag;
      }
      BBss.set( bbssx, bbssy, bbssz );
      kCall++;

//--- determine step versus field variation.
//------------------------------------------
      if( fabs( bbssx-bbsx ) > BBvar.x() ) bDivide = true;
      if( fabs( bbssy-bbsy ) > BBvar.y() ) bDivide = true;
      if( fabs( bbssz-bbsz ) > BBvar.z() ) bDivide = true;
      if( bDivide ) {
	if( fabs( step ) > stepm ) {
	  step /= 2.;
	  if( fabs( step ) < stepm ) {
	    if( step > 0.) step = stepm;
	    else step = - stepm;
	  }
	  continue;
	}
      }

//--- determine position and direction after full step.
//-----------------------------------------------------
      vDcr = step * vDcs;
      dv2ss = ccp * vDcr.cross( BBss );
      vDcss = vDcs + dv2ss;
      vPoss = vPos + vDcr + 0.5 * step * dv2ss;

//--- check tracking inversion (vs z-axis).
//-----------------------------------------
      int sigZ = 1;
      if( vDcss.z() < 0. ) sigZ = -1;
      if( sigZ == - sigZv ) {
	//vPoExit.setZ( vPoCoo.z() );
	errFlag = false;
	std::string err = "particle tracking inversion";
	mess.append( err );
	//std::cout << "mom " << 0.0003/chgOvMom << std::endl;
	break;
      }

//--- check end-of-traking and
//--- determine position and direction on exit.
//---------------------------------------------
      double zzs = vPos.z();
      double zzss = vPoss.z();
      if( kVers * (zzss-zzpout) >= 0 ) {
	double ffr = (zzpout-zzs) / (zzss-zzs);
	vPoss = vPos + ffr * (vPoss-vPos);
	vPoss.setZ( zzpout );
	vDcss = vDcs + ffr * (vDcss-vDcs);
	vDcss = vDcss.unit();
        kStep++;
	break;
      }

//--- prepare for the next step.
//------------------------------
      vDcs = vDcss;
      vDcs = vDcs.unit();
      vPos = vPoss;
      sigZv = sigZ;
      kStep++;

//--- conditional printing (tracking).
//------------------------------------
      if( tprint ) {
	std::cout << std::endl;
	std::cout << "CsRCTrackThruField tracking : " << std::endl;
        std::cout << "kStep, kCall, bDivide, step, kVers = "
	  	  << kStep << "  " << kCall << "  " << bDivide << "  "
		  << step << "  " << kVers << std::endl;
        std::cout << "xxs, yys, zzs = " << vPos.x() << "  " << vPos.y()
		  << "  " << vPos.z() << std::endl;
        std::cout << "dxxs, dyys, dzzs = " << vDcs.x() << "  " << vDcs.y()
		  << "  " << vDcs.z() << std::endl;
        std::cout << "Bx, By, Bz = " << bbsx << "  " << bbsy << "  "
	          << bbsz << std::endl;
      }

//--- go to the next step.
//------------------------
      if( fabs( step ) >= stepx) continue;
      if( bDivide ) continue;
      step *= 2.0;
      if( fabs( step ) > stepx) {
	if( step > 0.) step = stepx;
	else step = - stepx;
      }

    }
    if( kCall >= nCallMx ) {
      std::string err = "too many calls to mag field";
      mess.append( err );
      errFlag = false;
    }

//- output quantities.
//--------------------
    vPoExit = vPoss;
    vDcExit = vDcss.unit();

    if( !errFlag ) std::cout << std::endl << mess << std::endl;
//------------------------------------------------------------

//- conditional printing (output).
//--------------------------------
    if( print ) {
      std::cout << std::endl;
      std::cout << "CsRCTrackThruField output : " << std::endl;
      std::cout << "kStep, nCall = " << kStep << "  "
		<< kCall << std::endl;
      std::cout << "xout, yout, zout = " << vPoExit.x() << "  "
		<< vPoExit.y() << "  " << vPoExit.z() << std::endl;
      std::cout << "dxout, dyout, dzout = " << vDcExit.x() << "  "
		<< vDcExit.y() << "  " << vDcExit.z() << std::endl;
    }


    return  errFlag;

  }


