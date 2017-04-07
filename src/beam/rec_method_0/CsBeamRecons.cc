// $Id: CsBeamRecons.cc,v 2.71 2010/11/22 00:19:50 ybedfer Exp $
/*!
   \file CsBeamRecons.cc
   \brief Compass Beam reconstruction Class.
   \author  G. Khaustov
   \version $Revision: 2.71 $
   \date    $Date: 2010/11/22 00:19:50 $
*/
#include<iomanip>
#include <stdio.h>
#include "CsOpt.h"
#include "CsErrLog.h"
#include "CsBeamRecons.h"
#include "CsMCDigit.h"
#include "CsGeom.h"
#include "CsInit.h"
#include "CsEvent.h"

#include "TMath.h"

using namespace std;
using namespace CLHEP;

BackPropagated::BackPropagated()
{
    _X.reserve(6);
    _X_err.resize(4, std::vector<double>(3));
    _DX.reserve(2);
    _Y.reserve(6);
    _Y_err.resize(4, std::vector<double>(4));
    _DY.reserve(4);
}

void BackPropagated::GetX(double (&X)[6])
{
    for(int i=0; i<6; i++) {
	X[i] = _X[i];
    }
}

void BackPropagated::GetXerr(double (&X_err)[4][3])
{
    //only for old planes 
    for(int i=0; i<4; i++)
	for(int j=0; j<3; j++) {
	    X_err[i][j] = _X_err[i][j];
	}
}

void BackPropagated::GetDX(double (&DX)[2])
{
    DX[0] = _DX[0];
    DX[1] = _DX[1];
}

void BackPropagated::GetY(double (&Y)[6])
{
    for(int i=0; i<6; i++) {
	Y[i] = _Y[i];
    }
}

void BackPropagated::GetYerr(double (&Y_err)[4][4])
{
    //only for old planes
    for(int i=0; i<4; i++)
	for(int j=0; j<4; j++) {
	    Y_err[i][j] = _Y_err[i][j];
	}

}

void BackPropagated::GetDY(double (&DY)[2])
{
    DY[0] = _DY[0];
    DY[1] = _DY[1];
}

bool BackPropagated::Extrapolate(double BMS_p, double BmResol,
				 const CsHelix &Hin)
{
    //return false : sth went wrong
    //return true : ok
    //general equations:
    //Y = R33 * SF_Y +R34 * SF_DYDZ +R36 * dpp
    //DY= R43 * SF_Y +R44 * SF_DYDZ +R46 * dpp
    CsHelix exHel;
    double SF_Y, SF_DYDZ, SF_Y_err, SF_DYDZ_err;
    double SF_X, SF_DXDZ, SF_X_err, SF_DXDZ_err;

    //Tarnsport equations are defined for x and y known at 0 of Lau coordinate system
    //Lau coordinate system is shifted by 500 cm from compass coordinate system (calculated comparing detectors.dat and printout of transport equations parameters)
    //We take helix closest to z=0 abd extrapolate from it to z_lau=0
    //Then we do some conversions (Ask Martin why...)
    //if(!(Hin.Extrapolate(-5000,exHel))) return false;
    if(!(Hin.Extrapolate(-5220,exHel))) return false;
    SF_Y = exHel.getY();
    SF_Y_err = sqrt(exHel.getCov()[2]);
    SF_DYDZ = - exHel.getDYDZ()*1000;
    SF_DYDZ_err = sqrt(exHel.getCov()[9])*1000;

    SF_X = exHel.getX();
    SF_X_err = sqrt(exHel.getCov()[0]);
    SF_DXDZ = - exHel.getDXDZ()*1000;
    SF_DXDZ_err = sqrt(exHel.getCov()[5])*1000;

    // dp/p in units of percent. Note that "BMS_p" is the momentum scaled by
    // "160/BmMoment". => Therefore the reference momentum w.r.t. which one
    // must compute "dp/p" is 160.
    double BMS_dpp = (BMS_p-160)/1.6;
    double BMS_dpp_err = BmResol*BMS_p/1.6;
    double R11[4];
    double R12[4];
    double R21[4];
    double R22[4];
    double R33[4];
    double R34[4];
    double R36[4];
    double R43[4];
    double R44[4];
    double R46[4];
    //matrix elements:   R11    R12       R21     R22       R33    R34       R43      R44       R16      R26     R36     R46
    //BMS1           *   -0.741 -20.395   0.055   0.170 *    0.561 -43.498   0.044  -1.640 *    0.000   0.000   13.596   0.538
    //BMS2           *   -1.482 -22.671   0.055   0.170 *   -0.032 -21.484   0.044  -1.640 *    0.000   0.000    6.376   0.538
    //BMS3           *    0.068  50.581  -0.018   1.580 *   -0.894  16.304   0.032  -1.698 *    0.000   0.000    0.610  -0.037
    //BMS4           *    0.287  30.947  -0.018   1.580 *   -1.288  37.401   0.032  -1.698 *    0.000   0.000    1.076  -0.037
    R11[0]=-0.741;
    R11[1]=-1.482;
    R11[2]=0.068;
    R11[3]=0.287;

    R12[0]=-20.395;
    R12[1]=-22.671;
    R12[2]=50.581;
    R12[3]=30.947;

    R21[0]=0.055;
    R21[1]=0.055;
    R21[2]=-0.018;
    R21[3]=-0.018;

    R22[0]=0.170;
    R22[1]=0.170;
    R22[2]=1.580;
    R22[3]=1.580;

    R33[0]=0.561;
    R33[1]=-0.032;
    R33[2]=-0.894;
    R33[3]=-1.288;

    R34[0]=-43.498;
    R34[1]=-21.484;
    R34[2]=16.304;
    R34[3]=37.401 ;

    R36[0]=13.596;
    R36[1]=6.376;
    R36[2]=0.610;
    R36[3]=1.076;

    R43[0]=0.044;
    R43[1]=0.044;
    R43[2]=0.032;
    R43[3]=0.032;

    R44[0]=-1.640;
    R44[1]=-1.640;
    R44[2]=-1.698;
    R44[3]=-1.698;

    R46[0]=0.538;
    R46[1]=0.538;
    R46[2]=-0.037;
    R46[3]=-0.037;
    for(int i=0;i<4;i++){
	_Y[i] = R33[i] * SF_Y +R34[i] * SF_DYDZ +R36[i] * BMS_dpp;
	_Y_err[i][0]=fabs(R33[i] * SF_Y_err);
	_Y_err[i][1]=fabs(R34[i] * SF_DYDZ_err);
	_Y_err[i][2]=fabs(R36[i] * BMS_dpp_err);
	_Y_err[i][3]=sqrt( pow( R33[i] * SF_Y_err , 2)+ pow( R34[i] * SF_DYDZ_err ,2) + pow( R36[i] * BMS_dpp_err ,2)+ pow(5/sqrt(12.),2)); //why sqr(5/sqrt(12)) ??
	if(i==0)_DY[0]= R43[i] * SF_Y +R44[i] * SF_DYDZ +R46[i] * BMS_dpp;
	if(i==2)_DY[1]= R43[i] * SF_Y +R44[i] * SF_DYDZ +R46[i] * BMS_dpp;

	_X[i] = R11[i] * SF_X + R12[i] * SF_DXDZ;
	_X_err[i][0] = fabs(R11[i] * SF_X_err);
	_X_err[i][1] = fabs(R12[i] * SF_DXDZ_err);
	_X_err[i][2] = sqrt( pow( R11[i] * SF_X_err , 2)+ pow( R12[i] * SF_DXDZ_err ,2) + pow(5/sqrt(12.),2)); //why sqr(5/sqrt(12)) ??
	if(i==0)_DX[0]= R21[i] * SF_X + R22[i] * SF_DXDZ;
	if(i==2)_DX[1]= R21[i] * SF_X + R22[i] * SF_DXDZ;
    }
    _Y[4] = _Y[1] - ((_Y[1]-_Y[0])/13428.0) * 6794.0;
    _X[4] = _X[1] - ((_X[1]-_X[0])/13428.0) * 6794.0;
    _Y[5] = _Y[3] - ((_Y[3]-_Y[2])/12418.0) * 9520.0;
    _X[5] = _X[3] - ((_X[3]-_X[2])/12418.0) * 9520.0;

    
    return true;
}

//
CsBeamRecons* CsBeamRecons::instance_ = 0;
//

CsBeamRecons* CsBeamRecons::Instance() {
  if( instance_ == 0 ) {
    instance_ = new CsBeamRecons();
  }
  return( instance_ );
}
//
CsBeamRecons:: CsBeamRecons(void)
{
        string tag, key, str;
        int n;
        CsOpt* opt = CsOpt::Instance();
        clear();
//
        tag="BeamRecons";
        key="BMSrec";
        if( opt->getOpt( tag, key, n ) ) BMSrec=n; else BMSrec=1;
	if(BMSrec) {
                BMS = new CsBMSrecons;
                cout << "CsBeamRecons::  BMS reconstruction is ON" << endl;
        }
        else   cout << "CsBeamRecons::  BMS reconstruction is OFF!!" << endl;

	//          ***** BEAM CHARGE and MOMENTUM *****
	// - Defaults are +160 GeV, i.e. what the BMS parameterization is tuned
	//  on.
	// - Yet, this parameterization is a good approximation for any momentum
	//  P... provided:
	//  - it's re-scaled w/:                   P/160.
	//  - given, the characteristics Bi(I) of the beam line magnets i, the
	//   beam line currents are re-scaled w/:  Bi^-1(Bi(I_160)*P/160).
	// - The approximation being valid up to non-linear terms in the
	//  multipole expansion of the magnets' S Bdl (resp. S Gdl for quads). 
	// - Best would be to have this beam charge and momentum parameters
	//  specified once and for all of coral. In the mean time...
        if (opt->getOpt(tag,"BmMoment",BmMoment)) BmMoment = fabs(BmMoment);
        else BmMoment = 160;
        if (opt->getOpt(tag,"BmChrg",  BmChrg ));
	else BmChrg = +1;
//
        key="BmResol";
        if( opt->getOpt( tag, key, BmResol ) ) ; else BmResol=0.005;
//
        key="MxTmdiff";
        if( opt->getOpt( tag, key, MxTmdiff ) ) MxTmdiff=fabs(MxTmdiff);
        else MxTmdiff=2.1;
	key="MxTmdiff_center";
	if( opt->getOpt( tag, key, MxTmdiff_cntr ) ) cout<<"CsBeamRecons::MxTmdiff_cntr :"<<MxTmdiff_cntr<<endl;
        else MxTmdiff_cntr=0;
	MxTmdiff_min=MxTmdiff_cntr-MxTmdiff;
	MxTmdiff_max=MxTmdiff_cntr+MxTmdiff;
	//
	double goodBMStrack_cntr=0;
	key="goodBMStrack_center";
	if( opt->getOpt( tag, key, goodBMStrack_cntr ) ) cout<<"CsBeamRecons::goodBMStrack_cntr :"<<goodBMStrack_cntr<<endl;
        else goodBMStrack_cntr=0;
	goodBMStrack_min=goodBMStrack_cntr-MxTmdiff;
	goodBMStrack_max=goodBMStrack_cntr+MxTmdiff;
	cout<<"Cut on BMStrackTime: +/-"<<MxTmdiff<<" around "<<goodBMStrack_cntr<<" => ["<<goodBMStrack_min<<","<<goodBMStrack_max<<"]"<<endl;
	//
	key="MaxSFHitTime";
	if( opt->getOpt( tag, key, MaxSFHitTime ) ) MaxSFHitTime=fabs(MaxSFHitTime);
        else MaxSFHitTime=1.5;
//
        key="Printlev";
        if( opt->getOpt( tag, key, Printlev ) ) ; else Printlev=0;
//
        if( CsInit::Instance()->IsAMonteCarloJob()) {
	  MC=true;
	  tag = "";
	  key = "make decoding";
	  string mode;
	  MCExect= CsOpt::Instance()->getOpt( tag, key, mode ) && (mode == "MCExact" );
        } else {
	  MC=false;
	  MCExect=false;
	}
//
	tag="BMSreco";
	key="RecMethod";
        if( opt->getOpt( tag, key, RecMethod) ) ; else RecMethod=0;
	if(RecMethod==3){
	  MaxSFHitTime=0;
	  cout<<"CsBeamRecons: RecMethod 3 is used with MaxSFHitTime: "<<MaxSFHitTime<<endl;
	}
	tag="BMSreco";
	key="thetaXmax";
	if( opt->getOpt( tag, key, thetaXcut) ) ; else thetaXcut=0;
	tag="BMSreco";
	key="thetaYmax";
	if( opt->getOpt( tag, key, thetaYcut) ) ; else thetaYcut=0;

	tag="BeamRecons";
	key="useTRAFFIC";
	if( opt->getOpt( tag, key, useTRAFFIC)) ; else useTRAFFIC=0;
	if(useTRAFFIC!=0){ //make sure that tracking is also done
	  if( !(opt->getOpt( "", "make tracking" ))) CsErrLog::Instance()->mes(elFatal," If you want to use tracks from TRAFFIC in beam-reconstruction, you have to have make tracking in options!"); 
#ifndef BMS_IDEAL_MC
	  if(MC) CsErrLog::Instance()->mes(elFatal," If you want to process MC, please set useTRAFFIC 0 in your beam.opt"); 
#endif
	}
	tag="BeamRecons";
	key="ThetaCorrCut4";
	if( opt->getOpt( tag, key, ThetaCorrCut4) ) ; else ThetaCorrCut4=0;
	tag="BeamRecons";
	key="ThetaCorrCut3";
	if( opt->getOpt( tag, key, ThetaCorrCut3) ) ; else ThetaCorrCut3=0;
	cout<<"CsBeamRecons: ThetaCorrCut3(4): "<<ThetaCorrCut3<<" ("<<ThetaCorrCut4<<")"<<endl;

	//     ***** READ IN BACKTRACKING COEFFs *****
	// - These coeffs are used to evaluate the quality of the association
	//  of a BMS and a beamTelescope tracks by backtracking the latter to
	//  the abscissae of the former, cf. "FillChi2" infra.
	// - 2 types of inputs...
        //  i) "BackTrack Default": defaults to built-in's,
	// ii) "BackTrack ResBMS/Corr": specifies actual coeffs.
	for (int i = 0; i<bTNPLANES; i++) { bTResBMS[i]=bTCorr[i] = 0; }
	tag = "BeamRecons";
	key = "BackTrack Default";                    // ***** BUILT-IN DEFAULTS
	int bTInput = 0;
	if (opt->getOpt(tag,key)) {
	  bTInput = 0x3; bTDefault = true;
	}
	else bTDefault = false;
	key = "BackTrack ResBMS";         // ***** COEFFs READ from OPTIONS FILE
	list<string> coeffs; if (opt->getOpt(tag,key,coeffs)) {
	  if (bTInput&0x1)
	    CsErrLog::mes(elFatal,
 "\"BackTrack ResBMS\" option entered while \"BackTrack Default\" also booked");
	  list<string>::iterator is; int ns;
	  for (is = coeffs.begin(), ns = 0; is!=coeffs.end(); is++, ns++)
	    istringstream(*is)>>bTResBMS[ns];
	  if (ns<bTNPLANES)
	    CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "# of \"BackTrack ResBMS\" coefficients entered = %d => < %d",
			  ns,bTNPLANES);
	  bTInput |= 0x1;
	}
	key = "BackTrack Corr";
	if (opt->getOpt(tag,key,coeffs)) {
	  if (bTInput&0x2)
	    CsErrLog::mes(elFatal,
 "\"BackTrack Corr\" option entered while \"BackTrack Default\" also booked");
	  list<string>::iterator is; int ns;
	  for (is = coeffs.begin(), ns = 0; is!=coeffs.end(); is++, ns++)
	    istringstream(*is)>>bTCorr[ns];
	  if (ns<bTNPLANES)
	    CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "# of \"BackTrack Corr\" coefficients entered = %d => < %d",
			  ns,bTNPLANES);
	  bTInput |= 0x2;
	}
	if (bTInput!=0x3 && !CsInit::Instance()->IsAMonteCarloJob())
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Backtracking option coeffs not fully entered, required are \"BackTrack ResBMS/Corr\" w/ %d values each or opt for \"BackTrack Default\"",bTNPLANES);
	if (!bTDefault) {
	  printf("=========BeamRecons BackTrack\n");
	  for (int i = 0; i<bTNPLANES; i++)
	    printf("%d %5.2f %5.2f\n",i,bTResBMS[i],bTCorr[i]);
	  printf("=============================\n");
	}

	tag="BeamRecons";
	key="Debug";
	if( opt->getOpt( tag, key, DebugMode) ) ; else DebugMode=0;
	if(DebugMode==1){
	  cout<<"CsBeamRecons: Debug-Mode is switched on!"<<endl;
	}
	tag="BeamRecons";
	key="HistogramLevel";
	if( opt->getOpt( tag, key, histoLevel) ) ; else histoLevel=0;
	cout<<"CsBeamRecons: histoLevel: "<<histoLevel<<endl;

	tag="BeamRecons";
	key="F1_Res";
	if( opt->getOpt( tag, key, F1bin) ) ; else F1bin=1;
	cout<<"CsBeamRecons: F1_Res: "<<F1bin<<endl;
	
	tag="BeamRecons";
	key="TotMeanCut";
	if( opt->getOpt( tag, key, TotMeanCut) ) ; else TotMeanCut=2.2;
	cout<<"CsBeamRecons: TotMeanCut: "<<TotMeanCut<<endl;
	double TotMeanCut_cntr;
	key="TotMeanCut_center";
	if( opt->getOpt( tag, key, TotMeanCut_cntr ) ) cout<<"CsBeamRecons::TotMeanCut_cntr :"<<TotMeanCut_cntr<<endl;
        else TotMeanCut_cntr=0;
	TotMeanCut_min=TotMeanCut_cntr-TotMeanCut;
	TotMeanCut_max=TotMeanCut_cntr+TotMeanCut;
	cout<<"Cut on TotalMeanTime: +/-"<<TotMeanCut<<" around "<<TotMeanCut_cntr<<" => ["<<TotMeanCut_min<<","<<TotMeanCut_max<<"]"<<endl;

	tag="BeamRecons";
	key="year";
	if( opt->getOpt( tag, key, year)) cout<<"CsBeamRecons: year: "<<year<<endl;
	else year=2001;
	
	tag="BeamRecons";
	key="doRandomAna";
	if( opt->getOpt( tag, key,doRandom)){
	  if(doRandom!=0) cout<<"Special Analysis of random trigger will be performed!"<<endl;
	}
	else doRandom=0;
	
	tag="BeamRecons";
	key="killBoS";
	if( opt->getOpt( tag, key,killBoS)){
	  if(killBoS!=0) cout<<"Beginning of Spill will be ignored (events<"<<killBoS<<")!"<<endl;
	}
	else killBoS=0;

	tag="BeamRecons";
	key="SFtrigger_MaxSFplanes";
	if( opt->getOpt( tag, key,max_randomSFplanes)) {
	    if(max_randomSFplanes > 19) max_randomSFplanes = 19;
	    cout<<"BeamRecons: Max number of SF planes for SFtrigger: "<<max_randomSFplanes<<endl;
	}
	else max_randomSFplanes=13;

	tag="BeamRecons";
	key="randomSFnumber";
	if( opt->getOpt( tag, key,randomSFnumber)){
	  if(randomSFnumber>max_randomSFplanes)randomSFnumber=max_randomSFplanes;
	  cout<<"BeamRecons: randomSFnumber"<<randomSFnumber<<endl;
	}
	else randomSFnumber=11;

	tag="BeamRecons";
	key="randomSFnumber_SAS";
	if( opt->getOpt( tag, key,randomSFnumber_SAS)) {
	  if(randomSFnumber_SAS > max_randomSFplanes-11)
	      randomSFnumber_SAS = max_randomSFplanes-11;
	  cout<<"BeamRecons: Min number of SF planes for SFtrigger (zone after SM2): "
	      <<randomSFnumber_SAS<<endl;
	}
	else randomSFnumber_SAS=0;

	tag="BeamRecons";
	key="randomSpecMom";
	if( opt->getOpt( tag, key,GoodTrackMomCut)){
	  cout<<"BeamRecons: GoodTrackMomCut"<<GoodTrackMomCut<<endl;
	}
	else GoodTrackMomCut=0;
	tag="BeamRecons";
	key="RandomTargetCut";
	if( opt->getOpt( tag, key,RandomTargetCut)){
	  cout<<"BeamRecons: RandomTargetCut"<<RandomTargetCut<<endl;
	}
	else RandomTargetCut=0;
	tag="BeamRecons";
	key="EffTimeCutJ";
	if( opt->getOpt( tag, key,EffTimeCutJ)){
	  cout<<"BeamRecons: EffTimeCutJ"<<EffTimeCutJ<<endl;
	}
	else EffTimeCutJ=0.571;
	double EffTimeCutJ_cntr=0;
	key="EffTimeCutJ_center";
	if( opt->getOpt( tag, key,EffTimeCutJ_cntr)){
	  cout<<"BeamRecons: EffTimeCutJ_center"<<EffTimeCutJ_cntr<<endl;
	}
	EffTimeCutJ_min=EffTimeCutJ_cntr-EffTimeCutJ;
	EffTimeCutJ_max=EffTimeCutJ_cntr+EffTimeCutJ;
	cout<<"Cut for Efficiency of jap SF: +/-"<<EffTimeCutJ<<" around "<<EffTimeCutJ_cntr<<" => ["<<EffTimeCutJ_min<<","<<EffTimeCutJ_max<<"]"<<endl;
	

	tag="BeamRecons";
	key="EffTimeCut15";
	if( opt->getOpt( tag, key,EffTimeCut15)){
	  cout<<"BeamRecons: EffTimeCut15"<<EffTimeCut15<<endl;
	}
	else EffTimeCut15=0.438;
	//
	//get cuts for Veto
	tag="BeamRecons";
	key="VT01center";
	if( opt->getOpt( tag, key,TimeCut_VT_cntr[0])){}
	else TimeCut_VT_cntr[0]=0;
	key="VT02center";
	if( opt->getOpt( tag, key,TimeCut_VT_cntr[1])){}
	else TimeCut_VT_cntr[1]=0;
	key="VT01cut";
	if( opt->getOpt( tag, key,TimeCutVt[0])){}
	else TimeCutVt[0]=0;
	key="VT02cut";
	if( opt->getOpt( tag, key,TimeCutVt[1])){}
	else TimeCutVt[1]=0;
	cout<<"BeamRecons: Cut on VTI1: "<<TimeCut_VT_cntr[0]<<"+/-"<<TimeCutVt[0]<<endl;
	cout<<"BeamRecons: Cut on VTI2: "<<TimeCut_VT_cntr[1]<<"+/-"<<TimeCutVt[1]<<endl;
	TimeCut_VT_cntr[2]=0;
	TimeCutVt[2]=0;
	TimeCut_VT_cntr[3]=0;
	TimeCutVt[3]=0;
	//
	tag="BeamRecons";
        key="RandomTrackTimeCut";
        if( opt->getOpt( tag, key, RandomTrackTimeCut ) )cout<<"BeamRecons: RandomTrackTimeCut: "<<RandomTrackTimeCut<<endl;
	else RandomTrackTimeCut=9999999;
	//
	tag="BeamRecons";
        key="Printlev";
        if( opt->getOpt( tag, key, verb ) ) ; else verb=0;
//
        if(Printlev>0) {
	    cout << "CsBeamRecons::  Beam charge                               " 
                 << BmChrg << endl;
	    cout << "CsBeamRecons::  Beam momentum resolution                  " 
                 << BmResol << endl;
            if(!BMSrec)   cout << "CsBeamRecons::  Beam Moment                             " 
                               << BmMoment << endl;
	    else cout << "CsBeamRecons::  Max. track time difference (BMS-BmHod)    " 
                      << MxTmdiff << endl;
        }
//
        if (histoLevel==2 || useTRAFFIC==0) bookhist();
	else {
	  if (histoLevel==1) bookhist1();
	}
	tag="BeamRecons";
	key="doRescue";
	double doRescue_in;
	if( opt->getOpt( tag, key,doRescue_in));else doRescue_in=0;
	if(doRescue_in!=0) doRescue=true;
	else doRescue=false;
	tag="BeamRecons";
	key="useMultiples";
	double Mult_in;
	if( opt->getOpt( tag, key,Mult_in));else Mult_in=0;
	if(Mult_in!=0) {
	  cout<<"CsBeamRecons: useMultiples is activated!"<<endl;
	  useMultiples=true;
	}
	else useMultiples=false;

	//set up pointers to SF01/SF02-planes:
	detname[0]="FI01X1__";
	detname[1]="FI01Y1__";
	detname[2]="FI02X1__";
	detname[3]="FI02Y1__";
	detname[4]="BM01P1__";
	detname[5]="BM02P1__";
	detname[6]="BM03P1__";
	detname[7]="BM04P1__";
	detname[8]="FI03X1__";
	detname[9]="FI03Y1__";
	detname[10]="FI03U1__";
	detname[11]="FI04X1__";
	detname[12]="FI04Y1__";
	detname[13]="FI04U1__";
	detname[14]="FI05X1__";
	detname[15]="FI05Y1__";
	detname[16]="FI06X1__";
	detname[17]="FI06Y1__";
	detname[18]="FI06V1__";
	detname[19]="FI07X1__";
	detname[20]="FI07Y1__";
	detname[21]="GM10U1__";
	detname[22]="GM10V1__";
	detname[23]="GM10Y1__";
	detname[24]="GM10X1__";
	detname[25]="FI08X1__";
	detname[26]="FI08Y1__";
	bool detok;
	list <CsDetector*>  mydets = CsGeom::Instance()->getDetectors();
	list<CsDetector*>::iterator imydet;
	for(int i=0; i<27; i++) {
	  SFpnt[i]=0;
	  detok=false;
	  for( imydet = mydets.begin(); imydet != mydets.end(); imydet++ ) {
	    if((*imydet)->GetTBName()==detname[i] || (*imydet)->getName()==detname[i]) {
	      if(i<4 || i>=8){
		SFpnt[i]=(*imydet);
		detok=true;
	      }
	      if(i>3){
		BMSzoneDets.push_back((*imydet));
	      }
	    }
	  } //loop over detectorslist
	}// loop over detectorpointer
	BMS_zone=new CsZone(-13719,-6128,"BMSzone",BMSzoneDets);
	// set up pointers to detectors in SF/SI-telescope:
	detnameUpT[0]="FI01X1__";
	detnameUpT[1]="FI01Y1__";
	detnameUpT[2]="FI02X1__";
	detnameUpT[3]="FI02Y1__";
	detnameUpT[4]="FI15X1__";
	detnameUpT[5]="FI15Y1__";
	detnameUpT[6]="SI01X1__";
	detnameUpT[7]="SI01Y1__";
	detnameUpT[8]="SI01U1__";
	detnameUpT[9]="SI01V1__";
	detnameUpT[10]="SI02X1__";
	detnameUpT[11]="SI02Y1__";
	detnameUpT[12]="SI02U1__";
	detnameUpT[13]="SI02V1__";
	for(int i=0; i<=13; i++) {
	  UpTpnt[i]=0;
	  detok=false;
	  for( imydet = mydets.begin(); imydet != mydets.end(); imydet++ ) {
	    if((*imydet)->GetTBName()==detnameUpT[i] || (*imydet)->getName()==detnameUpT[i]) {
	      cout<<"CsBeamRecons: "<<detnameUpT[i]<<" was found!"<<endl;
	      UpTpnt[i]=(*imydet);
	      detok=true;
	    }
	  } //loop over detectorslist
	  if(!detok) cout<<"CsBeamRecons: "<<detnameUpT[i]<<" was not found!"<<endl;
	}// loop over detectorpointer
	
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons::clear() {
  list<CsBeam*>::iterator i;
  for( i=bmtracks.begin(); i!=bmtracks.end(); i++ ) {
    delete (*i);
  }
  bmtracks.clear();
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons:: bmreconsTR(list<CsTrack*> tracks)
{
  if(verb!=0)cout<<"Entering bmreconsTR..."<<endl;
  //check whether beginning of Spill should be ignored:
  if(killBoS!=0){
    if((int)CsEvent::Instance()->getEventNumberInBurst()<killBoS) return;
  }
  double dt=0,thetaX,thetaY,X,Y;
  double BMSmom, BMSspchq, BMStrtime, BMStmchiq, totalMeanTime;
  int BMSspec;
  double BMSThetaY;
  double goodRandomTrigSFtime=0;
  int goodBMSAngle=0; //0: not usable, 3,4: ok= number of part BMShits
  int goodBMSy; //-1: not usable, 1:Ok, both Station3 and 4 are available, 3,4: 0k=number of station which gave info [the complementary 
                //station is not available, 0: more than 1 "old" station missing
  int BMSnhits, BMStotnhits,trackmult,trackmultcut;
  int ncut=0; //Number of trackcomb which pass dt-cut
  int nfincut=0;
  int ncut_theta=0;
  int nBMStrackscut=0;  //Number of different BMS-tracks,which pass the dt-cut
  int nSFtrackscut=0;
  int BMStracknumber=99; //stores BMS-Tracknumber that has passed dt-cut 
  int SFtracknumber=999;
  int BMSntracks=0, BMSntracks_resc=0;
  int ngoodSFtrack,nSFfired=0,nGoodBMS,nSFcluster;
  int nBMSpassed,nSFpassed,BMSiOld,nBMSpassed_resc,nSFpassed_resc;
  int SFhasTime=0;
  int SFhasNoTime=0;
  int ntracksSFcl[8];
  int nCombDt;
  int goodTrBMSfired=0;
  int beamMult=0;
  double goodRandomTrigSpecMom=0;
  double BMSY[4]={0,0,0,0};
  for(int i=0;i<8;i++) ntracksSFcl[i]=0;
  unsigned int SFtrackIold;
  bool goodtrack=false; 
  bool trackort=true;
  bool zoneBT=false;
  bool skipThisEvent;
  bool dtPassed=false;
  bool oneTrigger=false;
  bool SFgoodTr=false;
  bool BMSgoodTr=false;
  bool SFlist=false;
  bool BMSlist=false;
  bool SFdtcut=false;
  bool BMSdtcut=false;
  bool goodTrBMSok=false;
  bool goodTrBMSok_resc=false;
  bool goodTrSFok=false;
  CsHelix startexhel;
  CsHelix  endhel;
  CsHelix  beginhel;
  CsHelix  goodRandomHelix;
  string tbn;
  list<CsTrack*>::iterator Ti;
  vector<int> partBMStracks;
  bool trig=false;
  bool goodTrigger=false;
  SFtrack thisCand;
  BeamComb thisComb;
  BackPropagated back_prop;
  SFtracks.clear();
  BeamCombV.clear();
  clear();
  goodSFhelixes.clear();
  // do analysis of random trigger:
  goodRandomTrigger=false;
  if(doRandom!=0) goodRandomTrigger=doRandomTrAna(goodRandomTrigSFtime,goodRandomTrigSpecMom,goodRandomHelix);
  SFEfficiency(4,goodRandomTrigSFtime);
  //goodRandomHelix.Print();
  
  //
  //      Preparations...
  //
  IT=MT=LT=OT=IMT=CT=false;
  oneTrigger=findTriggerMask(year);

  if (BMSrec) {               //  ********** BMS RECONSTRUCTION **********
    if (verb) cout<<"Triggering BMS-reconstruction...\n";
    BMS->bmsrec();
    BMSntracks = BMS->getBMSnmb();
    if (verb) cout<<"Done with BMS-rec.\n";
  }

  //   *************** CREATE VECTOR OF SF TRACKS FROM TRAFFIC ***************

  ngoodSFtrack = 0; trackmult = 0; trackmultcut = 0;
  double zTarget = CsGeom::Instance()->getTargetCenter();  // in mm
  for (Ti = tracks.begin(); Ti!=tracks.end(); Ti++) {
    //        ********** LOOP OVER TRACKS **********
    trig = false;
    vector<CsHelix> helix = (*Ti)->getHelices(); vector<CsHelix>::iterator Hi;
    zoneBT=false; for (Hi = helix.begin(); Hi!=helix.end(); Hi++) {
      if ((Hi)->getZ()<zTarget) zoneBT=true;
    }
    if (zoneBT) {        // ********** HAS HELIX UPSTREAM of TRAGET **********
      nSFcluster=0;
      goodtrack=checktracktiming( (*Ti),nSFcluster,ntracksSFcl);
      if(thetaXcut!=0 || thetaYcut!=0)trackort=checktrackorientation( (*Ti));
      if(goodtrack){
	trig=true;
	if(RecMethod==3) trig=false;
	ngoodSFtrack++;
      }
      thisCand.track=(*Ti);
      thisCand.Trig=trig;
      thisCand.inComb=false;
      if((*Ti)->hasMeanTime()){
	thisCand.time=(*Ti)->getMeanTime();
	if(mH1[46]!=NULL) mH1[46]->Fill(thisCand.time);
	SFhasTime++;
      }
      else {
	thisCand.time=9999;
	SFhasNoTime++;
	continue;
      }
      if(fabs(thisCand.time)<2.2) SFlist=true;
      if(Htrig2[10]!=NULL){
	if(IT && oneTrigger) Htrig2[10]->Fill(thisCand.time,0);
	if(MT && oneTrigger) Htrig2[10]->Fill(thisCand.time,1);
	if(LT && oneTrigger) Htrig2[10]->Fill(thisCand.time,2);
	if(OT && oneTrigger) Htrig2[10]->Fill(thisCand.time,3);
	if(IMT && oneTrigger) Htrig2[10]->Fill(thisCand.time,4);
      }
      if(RecMethod==3){ 
	if(MaxSFHitTime!=0){ //to be cleaned!
	  if(goodtrack) {
	    SFtracks.push_back(thisCand);
	    trackmultcut++;
	  }
	}
	else if(trackort){
	  SFtracks.push_back(thisCand);
	  trackmultcut++;
	}
      }// RecMethod 3
      else if(trackort && goodtrack){
	SFtracks.push_back(thisCand);
	trackmultcut++;
      }
      trackmult++;
    }
  }   // End loop over tracks
  if(mH1[106]!=NULL) mH1[106]->Fill(SFhasTime);
  //plot TrackMultiplicity:
  if(mH1[63]!=NULL) mH1[63]->Fill(trackmult);
  if(mH1[100]!=NULL) mH1[100]->Fill(trackmultcut);
  //plot multiplicity of tracks with a certain number of involved SF-stations
  for(int i=0; i<=7;i++)if(mH2[108]!=NULL) mH2[108]->Fill(i,ntracksSFcl[i]);
  goodTrigger=false;

  if (DebugMode==1) {
    nSFfired=SFEfficiency(0,0);
    //determine wether this is a special trigger:
    //either 4 out of 6 planes of SF3/4
    //or goodRandomTrigger
    if(doRandom){
      if(goodRandomTrigger) {
	goodTrigger=true;
      }
    }
    else if(SFTrigger()>=4) goodTrigger=true;
    if(goodTrigger){
      goodTrBMSfired=BMS->SFTrigger(goodRandomTrigSFtime,0,0);
      if(Ratrig2[2]!=NULL)Ratrig2[2]->Fill(goodRandomTrigSpecMom,goodTrBMSfired);
      SFEfficiency(1,0);
      //SFEfficiency(2,goodRandomTrigSFtime);
      if(Htrig2[9]!=NULL)Htrig2[9]->Fill(BMSntracks,ngoodSFtrack);
      if(oneTrigger) BMS->TriggerChecks(Trigger,1);
    }
  }
  
  if(mH1[172]!=NULL && goodTrigger) mH1[172]->Fill(SFhasTime);
  
  if(oneTrigger)BMS->TriggerChecks(Trigger,0);
  //
  if(trackmult>50000 && mH1[90]!=NULL)mH1[90]->Fill(CsEvent::Instance()->getEventNumberInBurst());
  if(mH2[91]!=NULL)mH2[91]->Fill(CsEvent::Instance()->getEventNumberInBurst(),trackmult);

  // PLOT MULTIPLICITIES FOR SF-PLANES:
  if(mH1[65]!=NULL) mH1[65]->Fill(SFMult[0]);
  if(mH1[66]!=NULL) mH1[66]->Fill(SFMult[1]);
  if(mH1[67]!=NULL) mH1[67]->Fill(SFMult[2]);
  if(mH1[68]!=NULL) mH1[68]->Fill(SFMult[3]);

  if(DebugMode==1) {     // MAKE PLOT SFHitMultvsWiremap
    makeWiremap();
    if(mH2[81]!=NULL)mH2[81]->Fill(SFMult[0],BMS->getBMSMult(3));
    if(mH2[82]!=NULL)mH2[82]->Fill(SFMult[1],BMS->getBMSMult(3));
  
    if(mH2[77]!=NULL)mH2[77]->Fill(SFMult[0],SFMult[1]);
    if(mH2[78]!=NULL)mH2[78]->Fill(SFMult[2],SFMult[3]);
    if(mH2[79]!=NULL)mH2[79]->Fill(SFMult[0],SFMult[2]);
    if(mH2[80]!=NULL)mH2[80]->Fill(SFMult[0],SFMult[3]);
    //plot Multiplicity vs tracktimes:
    //for this loop over SF-tracks:
    for(unsigned int SFtrackI=0; SFtrackI<SFtracks.size(); SFtrackI++) {
      if(mH2[62]!=NULL) mH2[62]->Fill(SFtracks[SFtrackI].time,trackmult);
      if(trackmult<150) if(mH1[64]!=NULL) mH1[64]->Fill(SFtracks[SFtrackI].time);
      if(mH2[69]!=NULL) mH2[69]->Fill(SFtracks[SFtrackI].time,SFMult[0]);
      if(mH2[70]!=NULL) mH2[70]->Fill(SFtracks[SFtrackI].time,SFMult[1]);
      if(mH2[71]!=NULL) mH2[71]->Fill(SFtracks[SFtrackI].time,SFMult[2]);
      if(mH2[72]!=NULL) mH2[72]->Fill(SFtracks[SFtrackI].time,SFMult[3]);
    }
  }// DebugMode
  if(SFlist) if(mH1[58]!=NULL) mH1[58]->Fill(0);

  //
  //      Get list of combinations, which pass dt-cut
  //
  int nCombTot=SFtracks.size()*BMSntracks;

  // *************** LOOP OVER BMS-TRACKS AND SF-LIST ***************

  ncut=0;
  nfincut=ncut_theta = 0; nGoodBMS = 0; nBMSpassed = 0; nSFpassed = 0;
  nCombDt = 0;
  BMSiOld=SFtrackIold = 9999;
  int BM05track = 0, BM06track = 0;
  if (verb) cout<<"Entering loop over SF-list and BMS-list..."<<endl;
  for (int BMSi = 0; BMSi<BMSntracks; BMSi++) {
    //     ********** LOOP OVER BMS-TRACKS **********
    BMS->getBMStrack(BMSi,BMSmom, BMSspchq, BMStrtime, BMStmchiq, BMSnhits, BMStotnhits,BMSspec); 
    BMStrtime=BMStrtime*F1bin;
    goodBMSAngle=BMS->getBMSThetaY(BMSi,BMSThetaY);
    if(mH1[50]!=NULL) mH1[50]->Fill(BMStrtime);
    if(Htrig2[11]!=NULL){
	if(IT && oneTrigger) Htrig2[11]->Fill(BMStrtime,0);
	if(MT && oneTrigger) Htrig2[11]->Fill(BMStrtime,1);
	if(LT && oneTrigger) Htrig2[11]->Fill(BMStrtime,2);
	if(OT && oneTrigger) Htrig2[11]->Fill(BMStrtime,3);
	if(IMT && oneTrigger) Htrig2[11]->Fill(BMStrtime,4);
    }
    if(fabs(BMStrtime)<2.2) BMSlist=true;
    //check how many good BMStracks there are b4 the combination with SF
    if(BMStrtime<2.2) nGoodBMS++;
    //
    for (unsigned int SFtrackI = 0; SFtrackI<SFtracks.size(); SFtrackI++) {
      //       ********** LOOP OVER SF LIST **********
      dt = BMStrtime-SFtracks[SFtrackI].time;
      if(mD1[0]!=NULL)mD1[0]->Fill(dt);
      if(mD1[1]!=NULL)mD1[1]->Fill(dt);
      if(SFtracks[SFtrackI].Trig) if(mD1[2]!=NULL)mD1[2]->Fill(dt);
      if(mH1[129]!=NULL)mH1[129]->Fill((BMStrtime+SFtracks[SFtrackI].time)/2);
      //           ********** TIME CUT **********
      if (dt>MxTmdiff_min && dt<MxTmdiff_max) {
	dtPassed = true; nCombDt++;
	if(mD1[6]!=NULL)mD1[6]->Fill(dt);
	if (SFtracks[SFtrackI].Trig) { 
	  ncut++; thisComb.Trig = true;
	}
	else      thisComb.Trig = false;
	thisComb.dt=dt;
	if(mD1[27]!=NULL)mD1[27]->Fill(BMStrtime);
	totalMeanTime=99999;
	if(SFtracks[SFtrackI].time!=9999){
	  if(mD1[29]!=NULL)mD1[29]->Fill(SFtracks[SFtrackI].time);
	  totalMeanTime=(BMStrtime+SFtracks[SFtrackI].time)/2;
	  if(mD1[30]!=NULL)mD1[30]->Fill(totalMeanTime);
	  if(doRandom!=0){ 
	    if(Ratrig1[6]!=NULL) Ratrig1[6]->Fill(goodRandomTrigSFtime-totalMeanTime);
	    totalMeanTime=totalMeanTime-goodRandomTrigSFtime;
	  }
	}
	thisComb.BMStrackNr = BMSi; thisComb.SFtrackNr = SFtrackI;
	thisComb.time = totalMeanTime;
	if (RecMethod==3) {
	  //     ********** GET HELIX WITH LARGEST Z **********
	  vector<CsHelix>trhel =SFtracks[SFtrackI].track->getHelices(); 
	  vector<CsHelix>::iterator Hi;
	  double endz=-99999999;
	  for(Hi = trhel.begin();Hi != trhel.end(); Hi++ ){
	    if( (Hi)->getZ()>endz ){
	      endz=(Hi)->getZ(); endhel=*Hi;
	    }
	  }
	  thetaY=atan(endhel.getDYDZ());
	  if(goodBMSAngle!=0){ //3 or 4-hit track
	    //Assuming the follwing correlation: thetaY=-0.5 (-BMSThetaY), "-"BMSThetaY because
	    //BMS-Ref-System is upside down!
	    if(goodBMSAngle==4){
	      if(thetaY>((0.5*BMSThetaY)+ThetaCorrCut4) || thetaY<((0.5*BMSThetaY)-ThetaCorrCut4)){
	      //if(thetaY>((0.5*BMSThetaY)+ThetaCorrCut4/0.894) || thetaY<((0.5*BMSThetaY)-ThetaCorrCut4/0.894)){
		if(ThetaCorrCut4!=0) continue;
	      }
	    }
	    if(goodBMSAngle==3){
	      if(thetaY>((0.5*BMSThetaY)+ThetaCorrCut3) || thetaY<((0.5*BMSThetaY)-ThetaCorrCut3)){
	      //if(thetaY>((0.5*BMSThetaY)+ThetaCorrCut3/0.894) || thetaY<((0.5*BMSThetaY)-ThetaCorrCut3/0.894)){
		if(ThetaCorrCut3!=0) continue;
	      }
	    }
	  }
	  //  ********** CUT ON ABSOLUTE TRACK-TIME **********
	  double extra = // If event's trigger jitter is larger than reference's
	    CsEvent::Instance()->getExtraTimeWidth();
	  if (totalMeanTime>TotMeanCut_min-extra &&
	      totalMeanTime<TotMeanCut_max+extra){
	    thisComb.Trig=true;
	    if (BMSspec==0) {// Only use BM05 tracks if there is no alternative
	      nfincut++;//used to keep track of tracks composed or "old" BMS stations 
	      ncut++;
	      if (BMSi!=BMSiOld) {
		nBMSpassed++; BMSiOld = BMSi;
		partBMStracks.push_back(BMSi);
	      }
	      thisComb.spec=BMSspec;
	    }
	    else if (BMSspec==1) {
	      BM05track++; thisComb.spec = BMSspec;
	    }
	    else {
	      BM06track++; thisComb.spec = BMSspec;
	    }
	    if(SFtrackI!=SFtrackIold){
	      SFtracks[SFtrackI].inComb=true;
	    }
	  }//trigger correlated beam candidates
	  else{
	    if(fabs(SFtracks[SFtrackI].time)<2.2) SFgoodTr=true;
	    if(fabs(BMStrtime)<2.2) BMSgoodTr=true;
	  }
	  //perform backpropagation from SF to BMS
	  back_prop.Extrapolate(BMSmom,BmResol,endhel);
	  back_prop.GetY(thisComb.Ex_bms_Y);
	  back_prop.GetYerr(thisComb.Ex_bms_Y_err);
	  back_prop.GetDY(thisComb.Ex_bms_DY);
	  
	  BMS->getBMSY(thisComb.BMStrackNr,BMSY);
	  FillChi2(BMSY,thisComb.Ex_bms_Y,thisComb.BackProp_Chi2,thisComb.BackProp_LH);
	  if (goodTrigger && BMSspec==0){ // Histo BackTracking, even if ! random trigger
	    if (mH1[179]) mH1[179]->Fill(thisComb.BackProp_Chi2);
	    if (mH1[183]) mH1[183]->Fill(thisComb.BackProp_LH);
	  }
	  thisComb.UseIt=0;
	}//RecM 3
	thisComb.BmResol = BmResol; thisComb.rescue = false;
	BeamCombV.push_back(thisComb);
	
	if(fabs(BMStrtime)<2.2) BMSdtcut=true;
	if(fabs(SFtracks[SFtrackI].time)<2.2) SFdtcut=true;
      } // End of time cut
    } // End of loop over SF list
  } // End of loop over BMS-tracks

  if (ncut==0 && doRescue) {
    //well, it seems as if the conventional reconstruction did not yield a satisfactory 
    //result. So, try if there is a combination of SF/SI and a 2hit BMStrack (rescue algo)

    //    ********** LOOP OVER BMS-RESCUED TRACKS **********
    nBMSpassed_resc=nSFpassed_resc=0;
    BMSntracks_resc=BMS->track_resc();
    BMSiOld=SFtrackIold=9999;
    for(int BMSi=0; BMSi<BMSntracks_resc; BMSi++){
      double BmResol_resc=BmResol;
      BMS->getBMStrack_resc(BMSi,BMSmom, BMSspchq, BMStrtime, BMStmchiq, BMSnhits,BMStotnhits,BmResol_resc); 
      BMStrtime=BMStrtime*F1bin;
      //loop over SF/SI tracks
      for(unsigned int SFtrackI=0; SFtrackI<SFtracks.size(); SFtrackI++) {
	dt=BMStrtime-SFtracks[SFtrackI].time;
	if(mH1[123]!=NULL)mH1[123]->Fill(dt);
	//do dt cut:
	if(dt>MxTmdiff_min && dt <MxTmdiff_max){
	  if(SFtracks[SFtrackI].Trig){ 
	    ncut++;
	    thisComb.Trig=true;
	  }
	  else thisComb.Trig=false;
	  thisComb.dt=dt;
	  totalMeanTime=99999;
	  if(SFtracks[SFtrackI].time!=9999){
	    totalMeanTime=(BMStrtime+SFtracks[SFtrackI].time)/2;
	    if(mH1[124]!=NULL)mH1[124]->Fill(totalMeanTime);
	    if(doRandom!=0){ 
	      if(Ratrig1[23]!=NULL) Ratrig1[23]->Fill(goodRandomTrigSFtime-totalMeanTime);
	      totalMeanTime=totalMeanTime-goodRandomTrigSFtime;
	    }
	  }
	  thisComb.BMStrackNr=BMSi;
	  thisComb.SFtrackNr=SFtrackI;
	  thisComb.time=totalMeanTime;
	  if(RecMethod==3){
	    //do cut on absolute track-time of beam-cand:
	    //
	    double extra = // If event's trigger jitter is larger than reference's
	      CsEvent::Instance()->getExtraTimeWidth();
	    //double extra=0;
	    //if(fabs(totalMeanTime)<TotMeanCut){ 
	    //if(totalMeanTime>TotMeanCut_min && totalMeanTime <TotMeanCut_max){
	    if(totalMeanTime>TotMeanCut_min-extra && totalMeanTime <TotMeanCut_max+extra){
	      thisComb.Trig=true;
	      nfincut++;
	      ncut++;
	      if(BMSi!=BMSiOld) {
		nBMSpassed_resc++;
		BMSiOld=BMSi;
		partBMStracks.push_back(BMSi);
	      }
	      if(SFtrackI!=SFtrackIold){
		SFtracks[SFtrackI].inComb=true;
	      }
	      if(mH1[170]!=NULL)mH1[170]->Fill(BMSmom);
	      thisComb.spec=0;
	    }//trigger correlated beam candidates
	  }//RecM 3
	  thisComb.BmResol=BmResol_resc;
	  thisComb.rescue=true;
	  thisComb.UseIt=0;
	  BeamCombV.push_back(thisComb);
	}//dt-cut
      }
    }// end of for-loops    
    if(mH1[125]!=NULL)mH1[125]->Fill(ncut);
    if(goodTrigger) if(mH1[128]!=NULL)mH1[128]->Fill(ncut);
    //count number of SFtracks in Combinations:
    nSFpassed_resc=0;
    for(unsigned int SFtrackI=0; SFtrackI<SFtracks.size(); SFtrackI++) {
      if(SFtracks[SFtrackI].inComb)nSFpassed_resc++;
    }
    if(mH2[127]!=NULL)mH2[127]->Fill(nBMSpassed_resc,nSFpassed_resc);
    if(ncut==0){
      //let's see how many 2-hit tracks there are with 
      //at least one hit in BMS01/02
      if(mH1[126]!=NULL)mH1[126]->Fill(BMSntracks_resc);
    }
    //let's see whether there is a 2hit track in the BMS inside a given interval
    for(int BMSi=0; BMSi<BMSntracks_resc; BMSi++){
      double BmResol_resc=BmResol;
      BMS->getBMStrack_resc(BMSi,BMSmom, BMSspchq, BMStrtime, BMStmchiq, BMSnhits, BMStotnhits,BmResol_resc); 
      BMStrtime=BMStrtime*F1bin;
      if(Ratrig1[24]!=NULL)Ratrig1[24]->Fill(BMStrtime-goodRandomTrigSFtime);
      if(BMStrtime-goodRandomTrigSFtime>goodBMStrack_min && BMStrtime-goodRandomTrigSFtime<goodBMStrack_max){
	goodTrBMSok_resc=true;
      }
    }
    if(Ratrig1[25]!=NULL){
      if(!goodTrBMSok_resc && !goodTrSFok)Ratrig1[25]->Fill(3); 
      if(goodTrBMSok_resc && goodTrSFok){Ratrig1[25]->Fill(2);}
      else{ 
	if(goodTrBMSok_resc)Ratrig1[25]->Fill(0);
	else{ 
	  if(goodTrSFok)Ratrig1[25]->Fill(1);  
	}
      }
    }
  } // End of rescue

  // Decide whether BM05/BM06 tracks have to be used:
  if (ncut==0) {
    ncut = BM05track;
  }
  if (ncut==0) {
    ncut=BM06track;
  }
  //
  if (verb) cout<<"Leaving loops over the lists...\n";

  //                        ********** FILL SOME HISTOS **********
  if(mH1[58]!=NULL){
    if(BMSlist)  mH1[58]->Fill(2);
    if(SFdtcut)  mH1[58]->Fill(1);
    if(BMSdtcut) mH1[58]->Fill(3);
  }
  if(ncut==0){
    if(mH1[52]!=NULL){
      if(!BMSgoodTr) mH1[52]->Fill(0);
      if(!SFgoodTr) mH1[52]->Fill(1);
      if(!BMSgoodTr && !SFgoodTr) mH1[52]->Fill(2);
    }
    //plot nCombDt:
    if(doRandom) {
      if(goodTrigger){
	if(mH1[109]!=NULL) mH1[109]->Fill(nCombDt);
      }
    }
    else if(mH1[109]!=NULL) mH1[109]->Fill(nCombDt);
    if(DebugMode==1){
      if(nCombDt==0){
	//loop over the SF & BMS lists again and check the tracktimes:
	for(int BMSi=0; BMSi<BMSntracks; BMSi++){
	  BMS->getBMStrack(BMSi,BMSmom, BMSspchq, BMStrtime, BMStmchiq, BMSnhits, BMStotnhits,BMSspec); 
	  BMStrtime=BMStrtime*F1bin;
	  for(unsigned int SFtrackI=0; SFtrackI<SFtracks.size(); SFtrackI++) {
	    if(mH2[110]!=NULL) mH2[110]->Fill(BMStrtime,SFtracks[SFtrackI].track->getMeanTime());
	  }
	}//loop over BMS-tracks
      }
      if(doRandom && goodTrigger){
	//loop over BMS & SF tracks individually and subtract the SFtrackTime from the spectrometer:
	//furthermore record which setup provided information inside a given time-cut:
	for(int BMSi=0; BMSi<BMSntracks; BMSi++){
	  BMS->getBMStrack(BMSi,BMSmom, BMSspchq, BMStrtime, BMStmchiq, BMSnhits, BMStotnhits,BMSspec); 
	  BMStrtime=BMStrtime*F1bin;
	  if(Ratrig1[7]!=NULL)Ratrig1[7]->Fill(BMStrtime-goodRandomTrigSFtime);
	  if(BMStrtime-goodRandomTrigSFtime>goodBMStrack_min && BMStrtime-goodRandomTrigSFtime<goodBMStrack_max){
	    goodTrBMSok=true;
	  }
	}
	for(unsigned int SFtrackI=0; SFtrackI<SFtracks.size(); SFtrackI++) {
	  double SFtime=SFtracks[SFtrackI].track->getMeanTime();
	  if(Ratrig1[8]!=NULL)Ratrig1[8]->Fill(SFtime-goodRandomTrigSFtime);
	  //if(SFtime-goodRandomTrigSFtime>MxTmdiff_min && SFtime-goodRandomTrigSFtime<MxTmdiff_max){
	  if(fabs(SFtime-goodRandomTrigSFtime)<MxTmdiff){
	    goodTrSFok=true;
	  }
	}
	if(Ratrig1[11]!=NULL){
	  if(!goodTrBMSok && !goodTrSFok)Ratrig1[11]->Fill(3); 
	  if(goodTrBMSok && goodTrSFok){Ratrig1[11]->Fill(2);}
	  else{ 
	    if(goodTrBMSok)Ratrig1[11]->Fill(0);
	    else{ 
	      if(goodTrSFok)Ratrig1[11]->Fill(1);  
	    }
	  }
	}
	//What's happend if !goodTrSFok
	if(Ratrig1[12]!=NULL){
	  //record SF/SI-tracks with time (any time!)
	  if(!goodTrSFok && SFhasTime<20)Ratrig1[12]->Fill(SFhasTime);
	  if(!goodTrSFok && goodTrBMSok)Ratrig1[12]->Fill(SFhasTime+20);
	}
	if(Ratrig2[0]!=NULL){
	  //record SF/SI-tracks with time vs those w/o time
	  if(!goodTrSFok)Ratrig2[0]->Fill(SFhasNoTime,SFhasTime);
	}
	if(Ratrig2[1]!=NULL){
	  //record SF/SI-tracks with time vs those w/o time for goodTrBMSok
	  if(!goodTrSFok && goodTrBMSok)Ratrig2[1]->Fill(SFhasNoTime,SFhasTime);
	}
	if(!goodTrBMSok)BMS->SFTrigger(goodRandomTrigSFtime,1,0);
      }//doRandom && goodTrigger
    }//DebugMode==1
  }//ncut==0
  if(!doRandom)CheckVeto(0,0);
  if(DebugMode==1){
    if(doRandom && goodTrigger){
      CheckVeto(0, goodRandomTrigSFtime);	 
      SFEfficiency(2,goodRandomTrigSFtime);
      if(Ratrig1[28]!=NULL)Ratrig1[28]->Fill(goodRandomHelix.getZ());
      if(Ratrig2[15]!=NULL)Ratrig2[15]->Fill(goodRandomHelix.getX(),goodRandomHelix.getY());
      if(Ratrig1[30]!=NULL)Ratrig1[30]->Fill(1000*atan(sqrt( pow(goodRandomHelix.getDXDZ(),2)+ pow(goodRandomHelix.getDYDZ(),2))));
      if(!goodTrSFok && ncut==0){
	int nSFplanes=SFEfficiency(3,goodRandomTrigSFtime);
	if(nSFplanes==0){
	  CheckVeto(2, goodRandomTrigSFtime);
	  if(Ratrig2[14]!=NULL)Ratrig2[14]->Fill(goodRandomHelix.getX(),goodRandomHelix.getY());
	  if(Ratrig1[29]!=NULL)Ratrig1[29]->Fill(1000*atan(sqrt( pow(goodRandomHelix.getDXDZ(),2)+ pow(goodRandomHelix.getDYDZ(),2))));
	  //extrapolate track-helix of spectrometer track to SF1/SF2
	  CsHelix exHel;
	  double zSF1X=-7598;
	  double zSF1Y=-7582;
	  double zSF2X=-2851;
	  double zSF2Y=-2867;
	  bool inAceptance[4];
	  for(int k=0;k<4;k++)inAceptance[k]=false;
	  goodRandomHelix.Extrapolate(zSF1X,exHel);
	  if(SFpnt[0]->inActiveArea (exHel.getX(), exHel.getY()))inAceptance[0]=true;
	  goodRandomHelix.Extrapolate(zSF1Y,exHel);
	  if(SFpnt[1]->inActiveArea (exHel.getX(), exHel.getY()))inAceptance[1]=true;
	  goodRandomHelix.Extrapolate(zSF2X,exHel);
	  if(SFpnt[2]->inActiveArea (exHel.getX(), exHel.getY()))inAceptance[2]=true;
	  goodRandomHelix.Extrapolate(zSF2Y,exHel);
	  if(SFpnt[3]->inActiveArea (exHel.getX(), exHel.getY()))inAceptance[3]=true;
	  //
	  int nSFAcept=0;
	  for(int k=0;k<4;k++)if(inAceptance[k])nSFAcept++;
	  if(Ratrig2[11]!=NULL){
	    for(int k=0;k<4;k++)if(inAceptance[k])Ratrig2[11]->Fill(k,nSFAcept);
	  }
	  //record spill number:
	  if(Ratrig1[27]!=NULL)Ratrig1[27]->Fill(CsEvent::Instance()->getBurstNumber());
	}//no SF planes have fired
	else CheckVeto(1, goodRandomTrigSFtime);
      }
      else {
	CheckVeto(1, goodRandomTrigSFtime);
      }
    }//goodTrigger
    //else{SFEfficiency(4,goodRandomTrigSFtime);}
    //determine how many different SF-tracks made it into a Combination:
    nSFpassed=0;
    for(unsigned int SFtrackI=0; SFtrackI<SFtracks.size(); SFtrackI++) {
      //plot x vs y, dx,dy @ -350mm:
      vector<CsHelix>helix= SFtracks[SFtrackI].track->getHelices();
      CsHelix exHTarget;
      if(SFtracks[SFtrackI].inComb){
	helix[1].Extrapolate(-350,exHTarget);
	if(mH2[113]!=NULL) mH2[113]->Fill(exHTarget.getX(),exHTarget.getY());
	if(mH1[114]!=NULL) mH1[114]->Fill(exHTarget.getDXDZ());
	if(mH1[115]!=NULL) mH1[115]->Fill(exHTarget.getDYDZ());
	nSFpassed++;
	if(goodTrigger){
	  if(mH2[116]!=NULL) mH2[116]->Fill(exHTarget.getX(),exHTarget.getY());
	  if(mH1[117]!=NULL) mH1[117]->Fill(exHTarget.getDXDZ());
	  if(mH1[118]!=NULL) mH1[118]->Fill(exHTarget.getDYDZ());
	}
	helix[1].Extrapolate(-1000,exHTarget);
	//position of all SFtracks in Comb
	if(mH2[59]!=NULL) mH2[59]->Fill(exHTarget.getX(),exHTarget.getY());
	//if(mH2[60]!=NULL) mH1[60]->Fill(X);
	//if(mH2[61]!=NULL) mH1[61]->Fill(Y);
	//position only for ncut==1:
	if(ncut==1){
	  if(mH2[122]!=NULL)mH2[122]->Fill(exHTarget.getX(),exHTarget.getY());
	}
      }//SFtrack in Comb
    }
  }//DebugMode==1
  if(mH2[53]!=NULL) mH2[53]->Fill(nBMSpassed,nSFpassed);
  if(doRandom && goodTrigger) if(Ratrig2[6]!=NULL)Ratrig2[6]->Fill(nBMSpassed,nSFpassed);
  if(nBMSpassed==1) if(mH1[87]!=NULL)mH1[87]->Fill(nSFpassed);
  if(mH1[31]!=NULL)mH1[31]->Fill(nfincut);
  if(mH1[176]!=NULL && goodTrigger)mH1[176]->Fill(nfincut);
  if(mH1[173]!=NULL)mH1[173]->Fill(BM05track);
  if(mH1[177]!=NULL && goodTrigger)mH1[177]->Fill(BM05track);
  if(mH1[174]!=NULL)mH1[174]->Fill(BM06track);
  if(mH1[178]!=NULL && goodTrigger)mH1[178]->Fill(BM06track);
  if(mH1[175]!=NULL && nfincut)mH1[175]->Fill(0);
  if(mH1[175]!=NULL && BM05track)mH1[175]->Fill(1);
  if(mH1[175]!=NULL && BM06track)mH1[175]->Fill(2);
  if(mH1[175]!=NULL && BM05track && !nfincut)mH1[175]->Fill(3);
  if(mH1[175]!=NULL && BM06track && !BM05track && !nfincut)mH1[175]->Fill(4);
  if(goodTrigger) {
      if(mH1[192]!=NULL && nfincut)mH1[192]->Fill(0);
      if(mH1[192]!=NULL && BM05track==1)mH1[192]->Fill(1);
      if(mH1[192]!=NULL && BM06track==1)mH1[192]->Fill(2);
      if(mH1[192]!=NULL && BM05track==1 && !nfincut)mH1[192]->Fill(3);
      if(mH1[192]!=NULL && BM06track==1 && !BM05track && !nfincut)mH1[192]->Fill(4);
  }

  if(!dtPassed) if(mH1[32]!=NULL)mH1[32]->Fill(0); 
  //
  //How many goodBMS b4 combination?
  if(mH1[51]!=NULL)mH1[51]->Fill(nGoodBMS);
  if(ncut==0 && BMSntracks!=0){
    if ( !dtPassed ) {
      if (mH1[32]!=NULL) mH1[32]->Fill(1);
    }
    else
      if (mH1[32]!=NULL) mH1[32]->Fill(2);
 
   //let's find out number of fired SF-planes:
    if(mH1[23]!=NULL)mH1[23]->Fill(nSFfired);
      
  }
  //Determine the number of part. BMS-planes in the Comb-sample:
  if(DebugMode==1 && ncut>1){
    if(nBMSpassed>1){
      int nThreeBMSPlanes=0;
      int nFourBMSPlanes=0;
      double smallest_Chi2=9999999;
      double scnd_smallest_Chi2=9999999;
      for(unsigned int i=0;i<partBMStracks.size();i++){
	BMS->getBMStrack(partBMStracks[i],BMSmom, BMSspchq, BMStrtime, BMStmchiq, BMSnhits, BMStotnhits,BMSspec); 
	if(BMSnhits==4) {
	  nFourBMSPlanes++;
	  if(doRandom && goodTrigger){
	    //look at BMS-space Chi2:
	    if(Ratrig1[34]!=NULL) Ratrig1[34]->Fill(BMSspchq);
	    if(BMSspchq<smallest_Chi2){
	      scnd_smallest_Chi2=smallest_Chi2;
	      smallest_Chi2=BMSspchq;
	    }
	  }
	}
	if(BMSnhits==3) nThreeBMSPlanes++;
      }
      BMS->CommonHits(partBMStracks);
      if(mH2[111]!=NULL) mH2[111]->Fill(nThreeBMSPlanes,nFourBMSPlanes);
      BMS->MultBMStracks(partBMStracks);
      if(doRandom && goodTrigger && nFourBMSPlanes>0){
	if(Ratrig1[36]!=NULL)Ratrig1[36]->Fill(nFourBMSPlanes);
	//space-Chi2 ratio
	if(Ratrig1[35]!=NULL)Ratrig1[35]->Fill(smallest_Chi2/scnd_smallest_Chi2);
      }
    }
  }
  //Find out how many ncut>1-events can be recuperated:
  if(ncut>1 && BMSntracks==1){
    //let's see how many SF-tracks there are:
    if(mD1[25]!=NULL)mD1[25]->Fill(ngoodSFtrack);
  }
  //if(mH1[3]!=NULL)mH1[3]->Fill(BeamCombV.size());
  if(mH1[3]!=NULL)mH1[3]->Fill(ncut);
  if(mHisto1[0]!=NULL)mHisto1[0]->Fill(ncut);
  if(goodTrigger){
    if(mH1[105]!=NULL)mH1[105]->Fill(ncut);
  }
  //sort this info by trigger:
  if(Htrig2[8]!=NULL)Htrig2[8]->Fill(BMSntracks,ngoodSFtrack);
  if(Htrig2[0]!=NULL && ncut==0)Htrig2[0]->Fill(BMSntracks,ngoodSFtrack);
  //compute Filter-info:
  int FiEntry=0;
  if(ncut==0)FiEntry=0+Filter_Flag;
  if(ncut==1)FiEntry=2+Filter_Flag;
  if(ncut>1)FiEntry=4+Filter_Flag;
  if(oneTrigger){ 
    if(Htrig1[0]!=NULL)Htrig1[0]->Fill(ncut);
    if(Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,0);
    if(Htrig2[1]!=NULL && ncut==0)Htrig2[1]->Fill(BMSntracks,ngoodSFtrack);
  }
  if(!oneTrigger && Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,10);
  if(IT){
    if(!oneTrigger && Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,11);
    if(oneTrigger){
      if(Htrig1[1]!=NULL)Htrig1[1]->Fill(ncut);
      if(Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,1);
      if(Htrig2[2]!=NULL && ncut==0)Htrig2[2]->Fill(BMSntracks,ngoodSFtrack);
    }
  }
  if(MT){
    if(!oneTrigger && Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,12);
    if(oneTrigger){
      if(Htrig1[2]!=NULL)Htrig1[2]->Fill(ncut);
      if(Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,2);
      if(Htrig2[3]!=NULL && ncut==0)Htrig2[3]->Fill(BMSntracks,ngoodSFtrack);
    }
  }
  if(LT){
    if(!oneTrigger && Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,13);
    if(oneTrigger){
      if(Htrig1[3]!=NULL)Htrig1[3]->Fill(ncut);
      if(Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,3);
      if(Htrig2[4]!=NULL && ncut==0)Htrig2[4]->Fill(BMSntracks,ngoodSFtrack);
    }
  }
  if(OT){
    if(!oneTrigger && Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,14);
    if(oneTrigger){
      if(Htrig1[4]!=NULL)Htrig1[4]->Fill(ncut);
      if(Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,4);
      if(Htrig2[5]!=NULL && ncut==0)Htrig2[5]->Fill(BMSntracks,ngoodSFtrack);
    }
  }
  if(IMT){
    if(!oneTrigger && Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,15);
    if(oneTrigger && IMT){
      if(Htrig1[5]!=NULL)Htrig1[5]->Fill(ncut);
      if(Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,5);
      if(Htrig2[6]!=NULL && ncut==0)Htrig2[6]->Fill(BMSntracks,ngoodSFtrack);
    }
  }
  if(CT){
    if(!oneTrigger && Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,16);
    if(oneTrigger && CT){
      if(Htrig1[15]!=NULL)Htrig1[15]->Fill(ncut);
      if(Htrig2[13]!=NULL)Htrig2[13]->Fill(FiEntry,6);
      if(Htrig2[7]!=NULL && ncut==0)Htrig2[7]->Fill(BMSntracks,ngoodSFtrack);
    }
  }
  if(ncut==0) if(mH1[8]!=NULL)mH1[8]->Fill(BMSntracks);
   
  if(DebugMode==1){
    int add_this=0;
    if(ngoodSFtrack==0) add_this=0;
    else add_this=1;
    if(Htrig1[7]!=NULL){
      if(oneTrigger && IT) Htrig1[7]->Fill(add_this); //IT
      if(oneTrigger && MT) Htrig1[7]->Fill(10+add_this); //MT
      if(oneTrigger && LT) Htrig1[7]->Fill(20+add_this); //LT
      if(oneTrigger && OT) Htrig1[7]->Fill(30+add_this); //OT
      if(oneTrigger && IMT) Htrig1[7]->Fill(40+add_this); //IMT
      if(oneTrigger && IMT) Htrig1[7]->Fill(50+add_this); //CT
    }
    if(oneTrigger) if(Htrig1[8]!=NULL) Htrig1[8]->Fill(BMSntracks);
    if(oneTrigger && IT) if(Htrig1[9]!=NULL) Htrig1[9]->Fill(BMSntracks);
    if(oneTrigger && MT) if(Htrig1[10]!=NULL) Htrig1[10]->Fill(BMSntracks);
    if(oneTrigger && LT) if(Htrig1[11]!=NULL) Htrig1[11]->Fill(BMSntracks);
    if(oneTrigger && OT) if(Htrig1[12]!=NULL) Htrig1[12]->Fill(BMSntracks);
    if(oneTrigger && IMT) if(Htrig1[13]!=NULL) Htrig1[13]->Fill(BMSntracks);
    if(oneTrigger && CT) if(Htrig1[14]!=NULL) Htrig1[14]->Fill(BMSntracks);
  }//DebugMode

  //
  //
  //    loop over Combinations and fill them in CsBeam-Obj in the order of 
  //    increasing time-difference
  //    for now a CsBeam-Obj is submitted only if there is only one BeamComb, 
  //    which has the Trig-flag
  //
  //

  //so, check how many of these are there:
  skipThisEvent=false;
  int nTrigComb=0;
  
  for(unsigned int BMCombI=0; BMCombI<BeamCombV.size(); BMCombI++) {
    //if(SFtracks[BeamCombV[BMCombI].SFtrackNr].Trig) nTrigComb++;
    if(BeamCombV[BMCombI].Trig && BeamCombV[BMCombI].spec==0)nTrigComb++;     
  }
  if(nTrigComb>1)HandleAmbig(nBMSpassed,nSFpassed,goodTrigger);
  //
  //the combination is to be given to CORAL:
  //    * in case there is only one
  //    * in case there is only one BMS-measurement:
  if(nTrigComb==1) skipThisEvent=false;
  vector<double> Chi2_ncut_gt1_nSF_eq_1;
  vector<double> Prob_ncut_gt1_nSF_eq_1;
  vector<double> Chi2_ncut_gt1_nBMS_eq_1;
  vector<double> Prob_ncut_gt1_nBMS_eq_1;
  //
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ////loop over Combinations:
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //
  if(verb!=0)cout<<"Entering loop over combinations..."<<endl;
  //
  for(unsigned int BMCombI=0; BMCombI<BeamCombV.size(); BMCombI++) {
    int category=0;
    //  save Combination to CsBeam:
    //get helix with largest z:
    vector<CsHelix>trhel =SFtracks[BeamCombV[BMCombI].SFtrackNr].track->getHelices(); 
    vector<CsHelix>::iterator Hi;
    double endz=-99999999;
    for(Hi = trhel.begin();Hi != trhel.end(); Hi++ ){
      if( (Hi)->getZ()>endz ){
	endz=(Hi)->getZ();
	endhel=*Hi;
      }
    }
    endz=0;
    for(Hi = trhel.begin();Hi != trhel.end(); Hi++ ){
      if( (Hi)->getZ()<endz ){
	endz=(Hi)->getZ();
	beginhel=*Hi;
      }
    }
    goodBMSAngle=0;
    goodBMSy=-1;
    BMSspec=0;
    if(!BeamCombV[BMCombI].rescue){
      BMS->getBMStrack(BeamCombV[BMCombI].BMStrackNr,BMSmom, BMSspchq, BMStrtime, BMStmchiq, BMSnhits, BMStotnhits,BMSspec);
      goodBMSAngle=BMS->getBMSThetaY(BeamCombV[BMCombI].BMStrackNr,BMSThetaY);
      goodBMSy=BMS->getBMSY(BeamCombV[BMCombI].BMStrackNr,BMSY);
    }
    else {
      double BmResol_resc=BmResol;
      BMS->getBMStrack_resc(BeamCombV[BMCombI].BMStrackNr,BMSmom, BMSspchq, BMStrtime, BMStmchiq, BMSnhits, BMStotnhits,BmResol_resc);
      BMS->moment_resc(BeamCombV[BMCombI].BMStrackNr,beginhel.getDYDZ(),&BMSmom, &BMSspchq);
      BMSspec=0;
      //cout<<"Rescue-algo used for this combination!"<<endl;
      category=1;
    }
    BMStrtime=BMStrtime*F1bin;
    //
    //Decide about usage of BM05 info:
    //
    bool tmp_trig = BeamCombV[BMCombI].Trig;
    if(nTrigComb>0){//not necessary to use BM05 or BM06 tracks
      if(BMSspec==1) {
	  continue;
      }
      if(BMSspec==2) {
	  continue;
      }
    }
    if(nTrigComb==0){//we need BM05-tracks, but only one:
      if(BM05track>1) {
	  continue;
      }
      category=3;
    }
    if(BM05track>0){//not necessary to use BM06 tracks
      if(BMSspec==2) {
	  continue;
      }
    }
    if(nTrigComb==0 && BM05track==0){//we need BM06-tracks, but only one:
      if(BM06track>1) {
	  continue;
      }
      category=4;
    }
    if (nTrigComb>1) {// ********** MULTIPLES (sorted out by Backprop) **********
      if (!CsInit::Instance()->IsAMonteCarloJob()) {
	if (useMultiples) {
	  if (BeamCombV[BMCombI].UseIt==0) continue;
	}
	else continue;
      }
      category = 2;
    }
    CsBeam* pntr = new CsBeam(999.,BMSnhits,BMStotnhits,BMStrtime,BMStmchiq,BMSspchq,99,9999,SFtracks[BeamCombV[BMCombI].SFtrackNr].track->getMeanTime(),99);
    
    CsBeam* beamPntr = new CsBeam(*SFtracks[BeamCombV[BMCombI].SFtrackNr].track);
    beamPntr->setBmdata(BeamCombV[BMCombI].time,BMSnhits,BMStotnhits,BMStrtime,BMStmchiq,BMSspchq,99,9999,SFtracks[BeamCombV[BMCombI].SFtrackNr].track->getMeanTime(),99);  
    //beamPntr->setBmdata(999,BMSnhits,BMStotnhits,BMStrtime,BMStmchiq,BMSspchq,99,9999,SFtracks[BeamCombV[BMCombI].SFtrackNr].track->getMeanTime(),99);
    
    double BmResol_act = BeamCombV[BMCombI].BmResol;
    vector< CsHelix > beamHel;
    //   ***** 1ST HELIX *****
    // Set err(mom) in cov matrix
    cov = beginhel.getCov(); cov[14] = pow((BmResol_act*BmChrg/BmMoment),2);
    double cop;
#ifdef BMS_IDEAL_MC
    // - This patch assigns to SFtracks the momentum of its parent CsTrack
    // which is expected to have been set by trafdic (= to the momentum
    // of the MC track associated to the CsTrack smeared).
    // - It DOES NOT affect the behaviour of CsBeamRecons on RD data.
    if (CsEvent::Instance()->isAMonteCarloEvent()) {
      //prepare Helix-list with which CsTrack is then fed...
      vector< CsHelix > beamHel_mc;
      CsHelix tmp0_mc( beginhel.getX(),beginhel.getY(),beginhel.getZ(),beginhel.getDXDZ(),beginhel.getDYDZ(),0, beginhel.getCov());
      beamHel_mc.push_back(tmp0_mc);
      CsHelix tmp1_mc( endhel.getX(),endhel.getY(),endhel.getZ(),endhel.getDXDZ(),endhel.getDYDZ(),0, endhel.getCov());
      beamHel_mc.push_back(tmp1_mc);
      //
      cop =
	SFtracks[BeamCombV[BMCombI].SFtrackNr].track->getHelices()[0].getCop();
      SFtracks[BeamCombV[BMCombI].SFtrackNr].track->setHelices(beamHel_mc);
      if (!cop) cop = 1.0/BMSmom;
    }
    else
#endif
      cop = 160/BmMoment*BmChrg/BMSmom;
    CsHelix tmp0( beginhel.getX(),beginhel.getY(),beginhel.getZ(),beginhel.getDXDZ(),beginhel.getDYDZ(),cop,cov);
    pntr->addHelix(tmp0);
    beamHel.push_back(tmp0);
    //   ***** 2ND HELIX *****
    cov = endhel.getCov(); cov[14] = pow((BmResol_act*BmChrg/BmMoment),2);
    CsHelix tmp1( endhel.getX(),endhel.getY(),endhel.getZ(),endhel.getDXDZ(),endhel.getDYDZ(),cop,cov);
    pntr->addHelix(tmp1);
    beamHel.push_back(tmp1);
    beamPntr->setHelices (beamHel);

    //check helices:
    pntr->setChi2(99);
    pntr->setTrigLabel(BeamCombV[BMCombI].Trig);
    beamPntr->setTrigLabel(BeamCombV[BMCombI].Trig);
    pntr->setClusters(SFtracks[BeamCombV[BMCombI].SFtrackNr].track->getClusters());
    pntr->setZones(SFtracks[BeamCombV[BMCombI].SFtrackNr].track->getZones());

    // ***** ADD "BMS"-ZONE *****
    pntr->addZone(*BMS_zone); beamPntr->addZone(*BMS_zone);

    // ***** BACK-PROPAGATION STUFF *****
    beamPntr->setBackPropLH(BeamCombV[BMCombI].BackProp_LH);
    beamPntr->setChi2CutFlag(BeamCombV[BMCombI].BackProp_LH<0.005);

    if (!CsInit::Instance()->IsAMonteCarloJob())
      // Set BMS bits in CsBeam's hit patterns.
      BMS->CsBMSrecons::setHitPattern(beamPntr, BeamCombV[BMCombI].BMStrackNr, BeamCombV[BMCombI].rescue);

    //get Angles
    thetaX=atan(endhel.getDXDZ());
    thetaY=atan(endhel.getDYDZ());
    //this Condition will be removed:
    if(!skipThisEvent && BeamCombV[BMCombI].Trig){  
      if(verb!=0)cout<<"Pushing back..."<<endl;
      bmtracks.push_back(beamPntr);
      //monitor beam-multiplicity:
      beamMult++;
      //monitoring histos:
      //*category vs trigger:
      bool FillIt=false;
      if(doRandom){
	if(goodTrigger)FillIt=true;
      }
      else FillIt=true;
      if(mHisto2[0]!=NULL && FillIt){
	if(oneTrigger){
	  if(IT)mHisto2[0]->Fill(category,1);
	  if(MT)mHisto2[0]->Fill(category,2);
	  if(LT)mHisto2[0]->Fill(category,3);
	  if(IMT)mHisto2[0]->Fill(category,4);
	  if(CT)mHisto2[0]->Fill(category,5);
	  
	}
	else mHisto2[0]->Fill(category,0);
      }
      //*momentum vs category:
      if(mHisto2[1]!=NULL)mHisto2[1]->Fill(BMSmom,category);
      //fill some histos for unambiguous case:
      if(ncut==1){
	if(ncut==1 && !goodTrBMSok && DebugMode==1) BMS->SFTrigger(goodRandomTrigSFtime,2,BeamCombV[BMCombI].BMStrackNr); 
	//record momentum spectrum for these candidates:
	//cout<<BMSmom<<endl;
	if(mH1[88]!=NULL)mH1[88]->Fill(BMSmom);
	if(mHisto1[1]!=NULL)mHisto1[1]->Fill(BMSmom);
	if(doRandom && goodTrigger) {
	    if(mH1[188]!=NULL && BeamCombV[BMCombI].BackProp_LH > 0.005)mH1[188]->Fill(BMSmom);
	    if(mH1[189]!=NULL && BeamCombV[BMCombI].BackProp_LH <= 0.005)mH1[189]->Fill(BMSmom);
	}
	if(goodTrigger && Ratrig2[8]!=NULL)Ratrig2[8]->Fill(goodRandomTrigSpecMom,BMSmom);
	if(goodTrigger && Ratrig1[22]!=NULL)Ratrig1[22]->Fill(BMSmom-goodRandomTrigSpecMom);
	//record number of BMS planes fired:
	if(mH1[89]!=NULL)mH1[89]->Fill(BMSnhits);

	if(mD1[9]!=NULL)mD1[9]->Fill(thetaX);
	if(mD1[10]!=NULL)mD1[10]->Fill(thetaY);
	if(mD2[11]!=NULL)mD2[11]->Fill(thetaX,thetaY);
	if(mD2[12]!=NULL)mD2[12]->Fill(BMSmom,thetaX);
	if(mD2[13]!=NULL)mD2[13]->Fill(BMSmom,thetaY);
	if(mD1[16]!=NULL)mD1[16]->Fill(endhel.getX());
	if(mH1[17]!=NULL)mH1[17]->Fill(endhel.getY());
	if(goodBMSAngle!=0 && mH2[95]!=NULL)mH2[95]->Fill(-BMSThetaY,thetaY);
	if(goodBMSAngle!=0 && mH2[96]!=NULL)mH2[96]->Fill(-BMSThetaY,thetaX);
	if(goodBMSy==3 && mH2[120]!=NULL)mH2[120]->Fill(BMSY[3],endhel.getY());
	if(goodBMSy==4 && mH2[120]!=NULL)mH2[120]->Fill(BMSY[2],endhel.getY());

	if(mH1[102]!=NULL)mH1[102]->Fill(BeamCombV[BMCombI].dt);
	if(DebugMode==1){
	  //let's see SF-Hitmaps for these tracks:
	  list<CsCluster*> BTcluster = SFtracks[BeamCombV[BMCombI].SFtrackNr].track->getClusters();
	  list<CsCluster*>::iterator Ci;
	  for(Ci = BTcluster.begin();Ci != BTcluster.end(); Ci++ ){
	    list<CsDigit*> ClDigits=(*Ci)->getDigitsList();
	    list<CsDigit*>::iterator ClDigitsI;
	    for(ClDigitsI = ClDigits.begin();ClDigitsI != ClDigits.end(); ClDigitsI++ ){
	      tbn=(*ClDigitsI)->getDet()->GetTBName();
	      if(tbn=="FI01X1__"){
		if(mD1[18]!=NULL)mD1[18]->Fill((*ClDigitsI)->getAddress());
		if(mD1[22]!=NULL)mD1[22]->Fill((*ClDigitsI)->getDatum());
	      }
	      //if(mD1[22]!=NULL)mD1[22]->Fill((*ClDigitsI)->getDatum());
	      if(tbn=="FI01Y1__" && mD1[19]!=NULL) mD1[19]->Fill((*ClDigitsI)->getAddress());
	      if(tbn=="FI02X1__" && mD1[20]!=NULL) mD1[20]->Fill((*ClDigitsI)->getAddress());
	      if(tbn=="FI02Y1__" && mD1[21]!=NULL) mD1[21]->Fill((*ClDigitsI)->getAddress()); 
	    }//loop over digits
	  }// loop over cluster
	  if(doRandom && goodTrigger){
	    if(Ratrig1[9]!=NULL)Ratrig1[9]->Fill(BMStrtime-goodRandomTrigSFtime);
	    if(Ratrig1[10]!=NULL)Ratrig1[10]->Fill(SFtracks[BeamCombV[BMCombI].SFtrackNr].track->getMeanTime()-goodRandomTrigSFtime);
	    if(BMSnhits==4 &&Ratrig1[15]!=NULL)Ratrig1[15]->Fill(goodRandomTrigSpecMom);
	  }
	}//DebugMode==1
      } //this was for trigger-correlated beam-candidates (events with one Cand)
      if(DebugMode==1) {
	  vector<double> raw_time;
	  vector<double> raw_y;
	  vector<int> track_hits;

	  if(BMS->getBMStrackHits(BeamCombV[BMCombI].BMStrackNr,track_hits))
	      for(int i=0; i<BMS->getBMSnhod(); i++) {
		  if(BMS->getBMSRawHits(i,raw_time,raw_y) && track_hits[i] < 9999)
		      if(HiBeamHitTimes!=NULL) HiBeamHitTimes->Fill(raw_time[track_hits[i]]);
		  raw_time.clear();
		  raw_y.clear();
	      }
	  track_hits.clear();
      }
    }//events with !ignore && trigger correlated
    else delete pntr;
    //  get some histos for this comb:(all comb)

    //for events with ncut>1 let's see
    if(BeamCombV[BMCombI].Trig){
      //plot angles vs totalmeantime:
      if(mH2[54]!=NULL)mH2[54]->Fill(BeamCombV[BMCombI].time,thetaX);
      if(mH2[55]!=NULL)mH2[55]->Fill(BeamCombV[BMCombI].time,thetaY);
      //plot angles vs ncut:
      if(mH2[56]!=NULL)mH2[56]->Fill(ncut,thetaX);
      if(mH2[57]!=NULL)mH2[57]->Fill(ncut,thetaY);

    }
    //check MeanTrackTimes sorted by trigger:
    if(IT && oneTrigger){
      if(mH1[33]!=NULL)mH1[33]->Fill(BeamCombV[BMCombI].time);
      if(mH1[36]!=NULL)mH1[36]->Fill(BMStrtime);
      if(mH1[39]!=NULL)mH1[39]->Fill(SFtracks[BeamCombV[BMCombI].SFtrackNr].time);
      if(mH2[42]!=NULL)mH2[42]->Fill(BMStrtime,SFtracks[BeamCombV[BMCombI].SFtrackNr].time);
      
    }
    if(MT && oneTrigger){
      if(mH1[34]!=NULL)mH1[34]->Fill(BeamCombV[BMCombI].time);
      if(mH1[37]!=NULL)mH1[37]->Fill(BMStrtime);
      if(mH1[40]!=NULL)mH1[40]->Fill(SFtracks[BeamCombV[BMCombI].SFtrackNr].time);
      if(mH2[43]!=NULL)mH2[43]->Fill(BMStrtime,SFtracks[BeamCombV[BMCombI].SFtrackNr].time);
    }
    if(LT && oneTrigger){
      if(mH1[35]!=NULL)mH1[35]->Fill(BeamCombV[BMCombI].time);
      if(mH1[38]!=NULL)mH1[38]->Fill(BMStrtime);
      if(mH1[41]!=NULL)mH1[41]->Fill(SFtracks[BeamCombV[BMCombI].SFtrackNr].time);
      if(mH2[44]!=NULL)mH2[44]->Fill(BMStrtime,SFtracks[BeamCombV[BMCombI].SFtrackNr].time);
    }
    
    //get Angles for these tracks:
    thetaX=atan(endhel.getDXDZ());
    thetaY=atan(endhel.getDYDZ());
    
    if(mD1[14]!=NULL)mD1[14]->Fill(thetaX);
    if(mD1[15]!=NULL)mD1[15]->Fill(thetaY);


    //         ***** BACKTRACKING BMS <- bT: HISTOGRAMS

    if (BeamCombV[BMCombI].Trig && // I.e. combi in-time w.r.t. trigger
	ncut==1) {
      bool OneOldPlaneMissing = goodBMSy>0 && goodBMSy<5; double delta;

      if (mH2[139]) {         //@ BMS1
	delta = BMSY[0]-BeamCombV[BMCombI].Ex_bms_Y[0]-bTCorr[0];
	if (goodBMSy==0) mH2[139]->Fill(delta,0); // All BMS-stations fired
	if (OneOldPlaneMissing && goodBMSy!=1)
	  mH2[139]->Fill(delta,1);      // One "old", !BMS1, BMS-station missing
      }
      if (mH2[140]) {         //@ BMS2
	delta = BMSY[1]-BeamCombV[BMCombI].Ex_bms_Y[1]-bTCorr[1];
	if (goodBMSy==0) mH2[140]->Fill(delta,0); // All BMS-stations fired
	if (OneOldPlaneMissing && goodBMSy!=2)
	  mH2[140]->Fill((delta),1);    // One "old", !BMS2, BMS-station missing
      }
      if (mH2[141]) {         //@ BMS3
	delta = BMSY[2]-BeamCombV[BMCombI].Ex_bms_Y[2]-bTCorr[2];
	if (goodBMSy==0) mH2[141]->Fill(delta,0); // All BMS-stations fired
	if (OneOldPlaneMissing && goodBMSy!=3)
	  mH2[141]->Fill(delta,1);      // One "old", !BMS3, BMS-station missing
      }  
      if (mH2[142]) {         //@ BMS4
	delta = BMSY[3]-BeamCombV[BMCombI].Ex_bms_Y[3]-bTCorr[3];
	if (goodBMSy==0)mH2[142]->Fill(delta,0);  // All BMS-stations fired
	if (OneOldPlaneMissing && goodBMSy!=4)
	  mH2[142]->Fill(delta,1);      // One "old", !BMS4, BMS-station missing
      }
      if (mH2[143] && goodBMSAngle!=0) //@ BMS4: residual in angle
	mH2[143]->Fill((BMSThetaY*1000-BeamCombV[BMCombI].Ex_bms_DY[0]),1);

      if (goodBMSy==0) { // All BMS-stations fired => Look at BMS5 and 6
	vector<double> raw_time, raw_y;
	double delta=0, y_BM05=0, y_BM06=0;
	double delta_BM05 = 9999999999.0;
	double delta_BM06 = 9999999999.0;
	bool ok_BM05 = false, ok_BM06 = true;

	if (BMS->getBMSRawHits(4,raw_time,raw_y) && !raw_time.empty()) {
	  ok_BM05 = true;
	  for (unsigned int i = 0; i<raw_time.size(); i++) {
	    delta = fabs(BeamCombV[BMCombI].time-raw_time[i]);
	    if(delta<delta_BM05) {
	      delta_BM05=delta; y_BM05=raw_y[i];
	    }
	  }
	}
	raw_time.clear(); raw_y.clear();
	if (BMS->getBMSRawHits(5,raw_time,raw_y) && !raw_time.empty()) {
	  ok_BM06 = true;
	  for (unsigned int i = 0; i<raw_time.size(); i++) {
	    delta=fabs(BeamCombV[BMCombI].time-raw_time[i]);
	    if(delta<delta_BM06) {
	      delta_BM06=delta; y_BM06=raw_y[i];
	    }
	  }
	}

	if (ok_BM05 && delta_BM05>1.5) ok_BM05 = false;
	if (ok_BM06 && delta_BM06>1.5) ok_BM06 = false;

	if (mH1[193] && ok_BM05)
	  mH1[193]->Fill(y_BM05-BeamCombV[BMCombI].Ex_bms_Y[4]);//@ BMS5
	if (mH1[194] && ok_BM06)
	  mH1[194]->Fill(y_BM06-BeamCombV[BMCombI].Ex_bms_Y[5]);//@ BMS6

	if (mH1[180]) mH1[180]->Fill(BeamCombV[BMCombI].BackProp_Chi2);
	if (mH1[182]) mH1[182]->Fill(BeamCombV[BMCombI].BackProp_LH);
      }
    }
    if (DebugMode==1) {

      //   ***** FILL a few MORE HISTOGRAMS useful for debugging *****

      // Check whether ncut>1 events with nBMSpassed=1 can be rescued via 
      // ThetaY(SF)-ThetaY(BMS) correlation
      if (ncut>1 && BeamCombV[BMCombI].Trig) {
	if (goodBMSAngle!=0 && mH2[97])  mH2[97] ->Fill(-BMSThetaY,thetaY);
	if (goodBMSy==3     && mH2[121]) mH2[121]->Fill(BMSY[3],endhel.getY());
	if (goodBMSy==4     && mH2[121]) mH2[121]->Fill(BMSY[2],endhel.getY());
      }
      if (BeamCombV[BMCombI].Trig) {
	if (!doRandom || goodTrigger) {
	  if (goodBMSAngle) {
	    if (thetaY<0.5*BMSThetaY-0.00025+0.000602 &&
		thetaY>0.5*BMSThetaY-0.00025-0.000602)
	      ncut_theta++;
	  }
	  else ncut_theta++;
	}
      }

      if (ncut==1 && BeamCombV[BMCombI].Trig) {
	if (goodBMSAngle && mH1[101]) mH1[101]->Fill((thetaY-0.5*BMSThetaY)*0.447);
	if (goodBMSAngle && mH2[99])  mH2[ 99]->Fill(-BMSThetaY,thetaY);
	
	// Test BMS-rescue algo:
	if (ncut==1 && BMSnhits==4) BMS->testRescueMethod(BeamCombV[0].BMStrackNr,thetaY);

	// X vs. Y plot for SF02 for CLEAN EVENTS:
	double cluY = 99, cluX = 99;
	//FI02X:
	const list<CsCluster*> &clustersX = SFpnt[2]->getMyClusters();
	if (clustersX.size()==1) cluX = clustersX.front()->getU();
	//FI02Y:
	const list<CsCluster*> &clustersY = SFpnt[3]->getMyClusters();
	if (clustersY.size()==1) cluY = clustersY.front()->getU();
	if (mH2[94] && cluX!=99 && cluY!=99) mH2[94]->Fill(cluX,cluY);
      } //ncut==1 && trigger-correlation

      // Let's have a look at BMS/SF yAngles for RandomTrigger
      if (BeamCombV[BMCombI].Trig){ // Should also be Trigger-correlated
	// Let's have a look at the events with 1 SF track and 1 BMS track
	if (nBMSpassed==1 && nSFpassed==1) {
	  // yy
	  if (goodBMSAngle && Ratrig2[16]) Ratrig2[16]->Fill(-BMSThetaY,thetaY);
	  if (goodBMSAngle && Ratrig1[32]) Ratrig1[32]->Fill((thetaY-0.5*BMSThetaY+0.00025)*0.894);
	}
	// Let's have a look at the events with 1 SF track and >1 BMS track
	if (nBMSpassed>1 && nSFpassed==1) {
	  if(goodBMSAngle && Ratrig2[17]) Ratrig2[17]->Fill(-BMSThetaY,thetaY);
	  if(goodBMSAngle && Ratrig1[33]) Ratrig1[33]->Fill((thetaY-0.5*BMSThetaY+0.00025)*0.894);
	}

	if (ncut==1) { // Look at Backpropagation-results:
	  if (goodBMSy==0) {  // All BMS-stations fired
	    if (mH2[144])
	      mH2[144]->Fill(BMSY[0],BeamCombV[BMCombI].Ex_bms_Y[0]);
	    if (mH2[145])
	      mH2[145]->Fill(BMSY[1],BeamCombV[BMCombI].Ex_bms_Y[1]);
	    if (mH2[146])
	      mH2[146]->Fill(BMSY[2],BeamCombV[BMCombI].Ex_bms_Y[2]);
	    if (mH2[147])
	      mH2[147]->Fill(BMSY[3],BeamCombV[BMCombI].Ex_bms_Y[3]);

	    // Plot calculated error (~components):
	    if(mH2[153]) {	    //BM01
	      mH2[153]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[0][0],0);
	      mH2[153]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[0][1],1);
	      mH2[153]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[0][2],2);
	      mH2[153]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[0][3],3);
	    }
	    if(mH2[154]) {	    //BM02
	      mH2[154]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[1][0],0);
	      mH2[154]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[1][1],1);
	      mH2[154]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[1][2],2);
	      mH2[154]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[1][3],3);
	    }
	    if(mH2[155]) {	    //BM03
	      mH2[155]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[2][0],0);
	      mH2[155]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[2][1],1);
	      mH2[155]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[2][2],2);
	      mH2[155]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[2][3],3);
	    }
	    if(mH2[156]) {	    //BM04
	      mH2[156]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[3][0],0);
	      mH2[156]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[3][1],1);
	      mH2[156]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[3][2],2);
	      mH2[156]->Fill(BeamCombV[BMCombI].Ex_bms_Y_err[3][3],3);
	    }
	  }
	  //plot Chi2
	  double Chi2=0; int nDof=0;
	  for (int i = 0; i<4; i++){
	    if (BMSY[i]<99999) {
	      Chi2 += pow(((BMSY[i]-BeamCombV[BMCombI].Ex_bms_Y[i]-bTCorr[i])/bTResBMS[i]),2);
	      nDof++;
	    }
	  }
	  if (nDof>1) Chi2=Chi2/(nDof-1);
	  else Chi2=-1;
	  if (mH2[148]!=NULL && nDof==4) mH2[148]->Fill(Chi2,0);
	  if (mH2[148]!=NULL && nDof==3) mH2[148]->Fill(Chi2,1);
	}
	if (ncut>1 && nSFpassed==1){
	  //record Chi2
	  double Chi2 = 0; int nDof = 0;
	  for (int i = 0; i<4; i++) {
	    if (BMSY[i]<99999) {
	      Chi2 += pow(((BMSY[i]-BeamCombV[BMCombI].Ex_bms_Y[i]-bTCorr[i])/bTResBMS[i] ),2);
	      nDof++;
	    }
	  }
	  double Time[4];
	  int Channel[4];
	  BMS->getBMS_HitInfo(BeamCombV[BMCombI].BMStrackNr, Time, Channel);
	  //
	  //plot Probabilities:
	  if(nDof==4){
	    if(mH2[157]!=NULL)mH2[157]->Fill(TMath::Prob(Chi2,nDof),0);
	    if(mH2[157]!=NULL)mH2[157]->Fill(TMath::Prob(Chi2,nDof-1),1);
	    if(mH2[157]!=NULL)mH2[157]->Fill(TMath::Prob(Chi2,nDof-2),2);
	    if(mH2[157]!=NULL)mH2[157]->Fill(TMath::Prob(Chi2,nDof-3),3);
	    if(mH2[157]!=NULL)mH2[157]->Fill(TMath::Prob(Chi2,nDof-4),4);
	  }
	  if(nDof==3){
	    if(mH2[157]!=NULL)mH2[157]->Fill(TMath::Prob(Chi2,nDof),10);
	    if(mH2[157]!=NULL)mH2[157]->Fill(TMath::Prob(Chi2,nDof-1),11);
	    if(mH2[157]!=NULL)mH2[157]->Fill(TMath::Prob(Chi2,nDof-2),12);
	    if(mH2[157]!=NULL)mH2[157]->Fill(TMath::Prob(Chi2,nDof-3),13);
	  }  
	  // Reduced Chi2
	  if (nDof>1)
	    Prob_ncut_gt1_nSF_eq_1.push_back(TMath::Prob(Chi2,nDof-1));
	  else
	    Prob_ncut_gt1_nSF_eq_1.push_back(-1);
	  if (nDof>1) Chi2=Chi2/(nDof-1);
	  else        Chi2=-1;
	  Chi2_ncut_gt1_nSF_eq_1.push_back(Chi2);
	}//ncut>1, nSFpassed==1
	if(ncut>1 && nBMSpassed==1) {
	  //record Chi2
	  double Chi2 = 0; int nDof = 0;
	  for (int i = 0; i<4; i++) {
	    if (BMSY[i]<99999){
	      Chi2 += pow(((BMSY[i]-BeamCombV[BMCombI].Ex_bms_Y[i]-bTCorr[i])/bTResBMS[i] ),2);
	      nDof++;
	    }
	  }
	  if (nDof>1)
	    Prob_ncut_gt1_nBMS_eq_1.push_back(TMath::Prob(Chi2,nDof-1));
	  else
	    Prob_ncut_gt1_nBMS_eq_1.push_back(-1);
	  if (nDof>1) Chi2=Chi2/(nDof-1);
	  else        Chi2=-1;
	  Chi2_ncut_gt1_nBMS_eq_1.push_back(Chi2);
	}//ncut>1 && nBMSpassed==1
      }
    }//DebugMode
#warning "TODO: removal of all combinations which contain these tracks!"
  }//loop over Comb-vector
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //HandleAmbig(nBMSpassed,nSFpassed,goodTrigger);
  //plot Chi2 of best comb and Chi2_best/Chi2_2ndbest
  if(DebugMode==1 && doRandom && goodTrigger){
    double bestChi2=1.0e8;
    double secBest=-1;
    //ncut>1 && nSF==1:
    for(unsigned int j=0;j<Chi2_ncut_gt1_nSF_eq_1.size();j++){
      if(Chi2_ncut_gt1_nSF_eq_1[j]>=0 && Prob_ncut_gt1_nSF_eq_1[j]>=0.005){ 
	if(Chi2_ncut_gt1_nSF_eq_1[j]<bestChi2){
	  secBest=bestChi2; //there have to be at least two entries since ncut>1!!
	  bestChi2=Chi2_ncut_gt1_nSF_eq_1[j];
	}
	else{
	  if(Chi2_ncut_gt1_nSF_eq_1[j]<secBest)secBest=Chi2_ncut_gt1_nSF_eq_1[j];
	}
      }
    }//loop over Chi2s
    //plot:
    if(Chi2_ncut_gt1_nSF_eq_1.size()>1){
      if(mH1[149]!=NULL)mH1[149]->Fill(bestChi2);
      if(mH1[150]!=NULL)mH1[150]->Fill(bestChi2/secBest);
      if(mH2[166]!=NULL)mH2[166]->Fill(bestChi2,bestChi2/secBest);
      //if(bestChi2/secBest<0.01)cout<<"****Chi2-ratio: "<<bestChi2/secBest<<" Best: "<<bestChi2<<" 2nd-Best: "<<secBest<<endl;
    }
    //ncut>1 && nBMS==1:
    for(unsigned int j=0;j<Chi2_ncut_gt1_nBMS_eq_1.size();j++){
      if(Chi2_ncut_gt1_nBMS_eq_1[j]>=0 && Prob_ncut_gt1_nBMS_eq_1[j]>=0.005){
	 if(Chi2_ncut_gt1_nBMS_eq_1[j]<bestChi2){
	   secBest=bestChi2; //there have to be at least two entries since ncut>1!!
	   bestChi2=Chi2_ncut_gt1_nBMS_eq_1[j];
	 }
	 else{
	   if(Chi2_ncut_gt1_nBMS_eq_1[j]<secBest)secBest=Chi2_ncut_gt1_nBMS_eq_1[j];
	 }
      }
    }//loop over Chi2s
    //plot:
    if(Chi2_ncut_gt1_nBMS_eq_1.size()>1){
      if(mH1[151]!=NULL)mH1[151]->Fill(bestChi2);
      if(mH1[152]!=NULL)mH1[152]->Fill(bestChi2/secBest);
      if(mH2[167]!=NULL)mH2[167]->Fill(bestChi2,bestChi2/secBest);
    }
    //look at Probabilities:
    int Prob_cnt=0;
    double BestProb=-1;
    double SecProb=-1;
    for(unsigned int j=0;j<Prob_ncut_gt1_nSF_eq_1.size();j++){
      if(Prob_ncut_gt1_nSF_eq_1[j]>=0.005){ 
	Prob_cnt++;
	
      }//Prob>0.05
      if(Prob_ncut_gt1_nSF_eq_1[j]>BestProb){
	SecProb=BestProb;
	BestProb=Prob_ncut_gt1_nSF_eq_1[j];
      }
      else{
	if(Prob_ncut_gt1_nSF_eq_1[j]>SecProb)SecProb=Prob_ncut_gt1_nSF_eq_1[j];
      }
      //if(mH1[159]!=NULL)mH1[159]->Fill(Prob_cnt);
    }//loop over Prob [nBMS>1 && nSF==1]
    //find out more about Prob_cnt==0:
    if(Prob_ncut_gt1_nSF_eq_1.size()>0 && Prob_cnt==0){
      //cout<<"*******"<<Prob_ncut_gt1_nSF_eq_1.size()<<endl;
      //cout<<"****Best Prob: "<<BestProb<<" Best: "<<bestChi2<<" 2nd-Best: "<<secBest<<endl;
      for(unsigned int j=0;j<Prob_ncut_gt1_nSF_eq_1.size();j++){
	cout<<"******Prob: "<<Prob_ncut_gt1_nSF_eq_1[j]<<" Chi2: "<<Chi2_ncut_gt1_nSF_eq_1[j]<<endl;
      }
    }
    //
    if(Prob_ncut_gt1_nSF_eq_1.size()>0){
      if(mH1[159]!=NULL)mH1[159]->Fill(Prob_cnt);
      if(mH2[158]!=NULL)mH2[158]->Fill(BestProb,SecProb/BestProb);
      if(mH2[161]!=NULL)mH2[161]->Fill(BestProb,SecProb);
    }
    for(unsigned int j=0;j<Prob_ncut_gt1_nBMS_eq_1.size();j++){
      if(Prob_ncut_gt1_nBMS_eq_1[j]>=0.005){ 
	Prob_cnt++;
	if(Prob_ncut_gt1_nBMS_eq_1[j]>BestProb){
	  SecProb=BestProb;
	  BestProb=Prob_ncut_gt1_nBMS_eq_1[j];
	}
	else{
	  if(Prob_ncut_gt1_nBMS_eq_1[j]>SecProb)SecProb=Prob_ncut_gt1_nBMS_eq_1[j];
	}
      }//Prob>0.05
    }//loop over Prob [nBMS>1 && nSF==1]
    if(Prob_ncut_gt1_nBMS_eq_1.size()>0){
      if(mH1[160]!=NULL)mH1[160]->Fill(Prob_cnt);
    }
  }//DebugMode
  //
  if(mH1[119]!=NULL) mH1[119]->Fill(beamMult);
  if(mHisto1[3]!=NULL)mHisto1[3]->Fill(beamMult);
  if(mH1[171]!=NULL && goodTrigger) mH1[171]->Fill(beamMult);


// Check BMS geometrical acceptance for events with 0 beam candidates but with good trigger
//  if(!beamMult && goodTrigger && doRandom) BmsGeometryCheck();
// Check BMS geometrical acceptance for events with good trigger
  if(goodTrigger && doRandom && DebugMode>0) BmsGeometryCheck();

  //cout<<"CsBeamRecons: beamMult: "<<beamMult<<endl;
  //cout<<"CsBeamRecons: CsBeamMult: "<<bmtracks.size()<<endl;
  if(mH1[98]!=NULL) mH1[98]->Fill(ncut_theta);
  if(mH2[112]!=NULL){
    if(doRandom!=0){
      if(goodTrigger)mH2[112]->Fill(ncut,ncut_theta);
    }
    else mH2[112]->Fill(ncut,ncut_theta);
  }
  if(mH1[45]!=NULL){
    if(IT && oneTrigger) mH1[45]->Fill(0);
    if(MT && oneTrigger) mH1[45]->Fill(1);
    if(LT && oneTrigger) mH1[45]->Fill(2);
  }
  if(verb!=0)cout<<"Leaving CsBeamRecons..."<<endl;


}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons:: bmrecons(void)
//
{
        double x, y, z, dxdz, dydz, Htime, Htmchiq;
        int Hnhit, Htot;
        bool Htrig;
        int i, nbms, nbmhod;
        int ntrack;
//
        if(Printlev>1)  cout << "++++++++++++++++++++++++ NEW BEAM EVENT +++++++++++++++++++++++" << endl;
        clear();
        ntrack=0;
        if(MCExect)  {
                                   reconsMcExect();
                                   if(Printlev>1) beamprint();
                                   return;
        }
	
//
//      BmFIHod reconstruction
//
        BmHod.bhodrecons();
        nbmhod=BmHod.getBmHodnmb();
        HiNScFbtrack->Fill((double)nbmhod);
        if(!nbmhod)  {
                HiNtrack->Fill((double)ntrack);
                if(Printlev>1) beamprint();
                return; 
        }
//
//      BMS reconstruction
//
        if(BMSrec){
           BMS->bmsrec();
           nbms=BMS->getBMSnmb();
           HiNBMStrack->Fill((double)nbms);
           if(!nbms) { 
               HiNtrack->Fill((double)ntrack); 
               if(Printlev>1) beamprint();
               return; 
           }
        } 
        else{                 // NO BMS (for MC and Hadrons) - save tracks
             for(i=0; i<nbmhod; i++){  
               BmHod.getBmFiHodtrack( i, x, y, z, dxdz, dydz, Hnhit, Htot, Htrig, Htime, Htmchiq);
               FillCovMx(BmResol*BmChrg/BmMoment);
               CsBeam* pntr = new CsBeam;
               CsHelix tmp( x, y, z, dxdz, dydz, BmChrg/BmMoment, cov );
               pntr->addHelix(tmp);
               //tmp.~CsHelix();
               pntr->setChi2(Htmchiq);
               pntr->setTrigLabel(Htrig);
               pntr->setBmdata(Htime, 0, 0, 0., 0., 0., Hnhit, Htot, Htime, Htmchiq);
               pntr->setClusters(BmHod.getClusters(i));
               pntr->setZones( BmHod.getZones()); 
               bmtracks.push_back(pntr);
//	
               HiBmom->Fill(BmMoment);
               HiXprof->Fill(x); 
               HiYprof->Fill(y); 
               HiXincl->Fill(dxdz); 
               HiYincl->Fill(dydz); 
               HiTime->Fill(Htime); 
               HiTMchiq->Fill(Htmchiq); 
             }
           HiNtrack->Fill((double)nbmhod); 
           if(Printlev>1) beamprint();
           return;
        } //for MC &Hadrons
//
//       find tracks with best time correlation between BMS and BmHod
//
        tmcorrel(nbmhod, nbms, ntrack);
        if(Printlev>1) beamprint();
//
        HiNtrack->Fill((double)ntrack); 
        HiNBNhod->Fill((double)nbmhod,(double)nbms);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons:: HandleAmbig(int nBMSpassed,int nSFpassed,bool goodTrigger){
//void CsBeamRecons:: HandleAmbig(){
  
  vector<CombCand> Comb; //[BMSi,SFi,BMCombI,Chi2]
  vector<CombCand> DeletedComb; //[BMSi,SFi,BMCombI,Chi2]
  for(unsigned int BMCombI=0; BMCombI<BeamCombV.size(); BMCombI++) {
    if(BeamCombV[BMCombI].Trig && !BeamCombV[BMCombI].rescue && BeamCombV[BMCombI].spec==0){
      //if(doRandom && goodTrigger
      //int nMultComb=0;
      int BMSi=BeamCombV[BMCombI].BMStrackNr;
      int SFi=BeamCombV[BMCombI].SFtrackNr;
      double Chi2=BeamCombV[BMCombI].BackProp_Chi2;
      bool commonIndex=false;
      bool commonSF=false;
      int hasOverwritten=-1;
      bool haskilled=false;
      
      if(Comb.size()>0){
	for(unsigned int MultI=0;MultI<Comb.size();MultI++){
	  if(BMSi==Comb[MultI].BMSi || SFi==Comb[MultI].SFi) commonIndex=true;
	  if(commonIndex){
	    double ratio;
	    double BestChi2;
	    double mom1=BMS->getBMSmom(BeamCombV[Comb[MultI].BMCombI].BMStrackNr);
	    if(Chi2>Comb[MultI].Chi2){//new track is worse
	      ratio=Comb[MultI].Chi2/Chi2;
	      BestChi2=Comb[MultI].Chi2;
	      mom1=BMS->getBMSmom(BeamCombV[Comb[MultI].BMCombI].BMStrackNr);
	      if(hasOverwritten>-1){//put Comb[MultI] to the postion of hasOverwritten
		CombCand thisComb;
		//cout<<"Overwriting: "<<hasOverwritten<<": "<<Comb[MultI].BMSi<<" "<<Comb[MultI].SFi<<" "<<Comb[MultI].BMCombI<<" "<<Comb[MultI].Chi2<<" "<<endl;
		DeletedComb.push_back(Comb[MultI]);
		Comb[hasOverwritten]=Comb[MultI];
		//mark MultI for deletion:
		Comb[MultI].BMSi=-1;
		haskilled=true;
	      }
	      else continue; 
	    }
	    else{//new comb is better
	      ratio=Chi2/Comb[MultI].Chi2;
	      BestChi2=Chi2;
	      if(hasOverwritten==-1){
		CombCand thisComb;
		//cout<<"Overwriting: "<<MultI<<": "<<BMSi<<" "<<SFi<<" "<<BMCombI<<" "<<Chi2<<" "<<endl;
		DeletedComb.push_back(Comb[MultI]);
		thisComb.BMSi=BMSi;thisComb.SFi=SFi;thisComb.BMCombI=BMCombI;thisComb.Chi2=Chi2;
		Comb[MultI]=thisComb;
		hasOverwritten=MultI;
	      }
	      else{//candidate has already overwritten -> delete the current one:
		DeletedComb.push_back(Comb[MultI]);
		//cout<<"Deleting: "<<MultI<<": "<<Comb[MultI].BMSi<<" "<<Comb[MultI].SFi<<" "<<Comb[MultI].BMCombI<<" "<<Comb[MultI].Chi2<<" "<<endl;
		Comb[MultI].BMSi=-1;
		haskilled=true;
	      }
	    }//new comb is better
	    //record ratio of Chi2s:
	    if(doRandom && goodTrigger){
	      if(nBMSpassed>1 && nSFpassed>1){ 
		if(mH1[165]!=NULL)mH1[165]->Fill(ratio);
		if(mH2[168]!=NULL)mH2[168]->Fill(BestChi2,ratio);
	      }
	      if(nBMSpassed>1){ 
		double mom2=BMS->getBMSmom(BeamCombV[BMCombI].BMStrackNr);
		if(mH2[169]!=NULL)mH2[169]->Fill(mom1,mom1-mom2);
	      }
	    }
	  }//common index
	}//loop over exiting combinations
	//remove the deleted tracks from the list:
	if(haskilled){
	  //loop over the list of tracks to find the candidates:
	  for(unsigned int MultI=0;MultI<Comb.size();MultI++){
	    bool removeThis=false;
	    if(Comb[MultI].BMSi==-1)removeThis=true;
	    if(removeThis){
	      vector<CombCand>::iterator itRemove = Comb.begin() + MultI;
	      Comb.erase(itRemove);
	    }
	  }//there is a track to be deleted!
	}
      }
      if(!commonIndex){
	CombCand thisComb;
	thisComb.BMSi=BMSi;thisComb.SFi=SFi;thisComb.BMCombI=BMCombI;thisComb.Chi2=Chi2;
	Comb.push_back(thisComb);
      }
    }//trigger correlation  
  }//loop over BeamCombV
  //try to put back deleted combinations:
  for(unsigned int delMultI=0;delMultI<DeletedComb.size();delMultI++){
    bool common=false;
    for(unsigned int MultI=0;MultI<Comb.size();MultI++){
      if(DeletedComb[delMultI].BMSi==Comb[MultI].BMSi || DeletedComb[delMultI].SFi==Comb[MultI].SFi){
	common=true;
	break;
      }
    }
    if (!common) Comb.push_back(DeletedComb[delMultI]);
  }

  //mark the chosen combinations in BeamCombV
  int Cnt_SF1=0;
  int Cnt_BMS1=0;
  int Cnt_bothGt1=0;
  int Cnt_passed=0;
  for(unsigned int BMCombI=0; BMCombI<BeamCombV.size(); BMCombI++) {
    for(unsigned int MultI=0;MultI<Comb.size();MultI++){
      if(Comb[MultI].BMCombI==BMCombI) {
	  if(goodTrigger){
	    if(mH1[181]!=NULL)mH1[181]->Fill(BeamCombV[BMCombI].BackProp_Chi2);
	    if(mH1[184]!=NULL)mH1[184]->Fill(BeamCombV[BMCombI].BackProp_LH);
	  }
	  if(BeamCombV[BMCombI].BackProp_LH>=0.005){//cut on P(Chi2)
	      if(mH1[190]!=NULL && doRandom && goodTrigger)mH1[190]->Fill(BMS->getBMSmom(BeamCombV[BMCombI].BMStrackNr));
	  }
	  else if(mH1[191]!=NULL && doRandom && goodTrigger)mH1[191]->Fill(BMS->getBMSmom(BeamCombV[BMCombI].BMStrackNr));

	  //Cut on P(chi^2) disabled...
	  BeamCombV[BMCombI].UseIt=1;
	  if(nBMSpassed>1 && nSFpassed==1)Cnt_SF1++;
	  if(nSFpassed>1 && nBMSpassed==1)Cnt_BMS1++;
	  if(nSFpassed>1 && nBMSpassed>1)Cnt_bothGt1++;
	  Cnt_passed++;
      }
    }
  }
  if(doRandom && goodTrigger){
    if(mH1[185]!=NULL)mH1[185]->Fill(Comb.size());
    if(nBMSpassed>1 && nSFpassed==1 && mH1[162]!=NULL)mH1[162]->Fill(Cnt_SF1);
    if(nSFpassed>1 && nBMSpassed==1 && mH1[163]!=NULL)mH1[163]->Fill(Cnt_BMS1);
    if(nSFpassed>1 && nBMSpassed>1  && mH1[164]!=NULL)mH1[164]->Fill(Cnt_bothGt1);
    if(mH2[186]!=NULL) mH2[186]->Fill(Comb.size(),Cnt_passed);
    if(Comb.size()==2 && Cnt_passed==0 && mH2[187]!=NULL) {
	double bestProb, secProb;
	if(BeamCombV[Comb[0].BMCombI].BackProp_LH < BeamCombV[Comb[1].BMCombI].BackProp_LH){
	    bestProb=BeamCombV[Comb[1].BMCombI].BackProp_LH;
	    secProb=BeamCombV[Comb[0].BMCombI].BackProp_LH;
	}
	else {
	    bestProb=BeamCombV[Comb[0].BMCombI].BackProp_LH;
	    secProb=BeamCombV[Comb[1].BMCombI].BackProp_LH;
	}
	mH2[187]->Fill(bestProb,secProb);
    }
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons:: FillChi2(double BMSY[],double Ex_BMSY[],double &Chi2,double &LH){
  static bool first = bTDefault; if (first) {
    first = false;
    // Defaults are only known for 2003 and 2004
    int runNum = CsEvent::Instance()->getRunNumber();
    if      (27000<runNum && runNum<33000) { // 2003
      bTResBMS[0]=4.54357;
      bTResBMS[1]=3.45162;
      bTResBMS[2]=1.98665;
      bTResBMS[3]=1.87227;
      bTCorr[0]= 7.55605;
      bTCorr[1]= 3.72391;
      bTCorr[2]=-2.72167;
      bTCorr[3]=-7.05849;
    }
    else if (33000<runNum && runNum<45000) { // 2004
      //numbers for 2004
      bTResBMS[0]=4.57;
      bTResBMS[1]=3.621;
      bTResBMS[2]=2.102;
      bTResBMS[3]=2.073;
      bTCorr[0]= 4.69;
      bTCorr[1]= 2.065;
      bTCorr[2]=-2.266;
      bTCorr[3]=-5.342;
    }
    else
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "\"BackTrack Default\" option booked while no default available for run #%d",
			  runNum);
    printf("=========BeamRecons BackTrack\n");
    for (int i = 0; i<bTNPLANES; i++)
      printf("%d %5.2f %5.2f\n",i,bTResBMS[i],bTCorr[i]);
    printf("=============================\n");
  }

  Chi2 = 0; int nDof = 0; for(int i = 0; i<4; i++) {
    if (BMSY[i]<99999) {
      Chi2=Chi2+pow(((BMSY[i]-Ex_BMSY[i]-bTCorr[i])/bTResBMS[i]),2); nDof++;
    }
  }
  // Compute LH:
  if (nDof>1) LH = TMath::Prob(Chi2,nDof-1);
  else        LH = -1;
  // Reduced Chi2:
  if (nDof>1) Chi2 = Chi2/(nDof-1);
  else        Chi2 = -1;
  //cout<<"Chi2: "<<Chi2<<endl;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CsBeamRecons:: doRandomTrAna(double &GoodTrackTime,double &GoodTrackMom,CsHelix &GoodTrackStartHelix){
  //get tracks from Spec and check the momentum:
  //
  int nplanes=max_randomSFplanes;//number of planes considered for SFTrigger
  int SFfired_ran[nplanes];
  list<CsTrack*>tracks = CsEvent::Instance()->getTracks();
  list<CsTrack*>::iterator Ti;
  if(Ratrig1[0]!=NULL) Ratrig1[0]->Fill(tracks.size());
  int ntrackMom=0;
  int nGoodTracks=0;
  double mom;
  bool hasGoodTrack=false;
  bool hasHaloCand=false;
  bool goodTr=false;
  double HaloTrackTime;
  map<CsDetector*,int> det2bit = CsEvent::Instance()->getDetectorMap();
  map<CsDetector*,int>::const_iterator id2b;
  int bit;
  int word;
  static double GoodTrackTgtX,GoodTrackTgtY,GoodTrackTgtDx,GoodTrackTgtDy,GoodTrackChi2;
  static double downX,downY,upX,upY;
  GoodTrackTime=9999;
      
  for(Ti = tracks.begin();Ti != tracks.end(); Ti++ ){
      goodTr = false;
    if(fabs((*Ti)->getMeanTime())>RandomTrackTimeCut) continue;
    mom=0;
    for(int i=0;i<nplanes;i++) SFfired_ran[i]=0;
    vector<CsHelix>helix= (*Ti)->getHelices();
    vector<CsHelix>::iterator Hi;
    CsHelix* starthel=NULL;
    double startz=9999999;
    for(Hi = helix.begin();Hi != helix.end(); Hi++ ){
      double cop=(Hi)->getCop();
      if(cop==0)continue;
      mom=fabs(1/cop);
      //look for helix with smallest z-value:
      if( (Hi)->getZ()<startz ){
	startz=(Hi)->getZ();
	starthel=&(*Hi);
      }
    }// loop over helices
    if(starthel == NULL) continue;
    if(mom!=0)ntrackMom++;
    if(Ratrig1[1]!=NULL)Ratrig1[1]->Fill(mom);
    //check FiredBitMask for FI3->Fi8 participation in track:
    const unsigned int* Firedmap=(*Ti)->getFiredDetsBitmap();
    //loop over interesting detectors
    int det_max=nplanes+8;
    for(int j=8;j<=det_max;j++){
      CsDetector* d=SFpnt[j];
      //cout<<"Checking "<<SFpnt[j]->GetTBName()<<endl;
      // find current detector's word/bit numbers 
      id2b = det2bit.find(d);
      word = int(((*id2b).second)/32);
      bit =      ((*id2b).second)%32 ;
      //cout<<"CsBeamRecons:: doRandomTrAna Detname: "<<SFpnt[j]->GetTBName()<<endl;
      //cout<<"CsBeamRecons:: doRandomTrAna word: "<<word<<endl;
      //cout<<"CsBeamRecons:: doRandomTrAna bit: "<<bit<<endl;
      if((Firedmap[word]&(1<<bit))!=0){
	//cout<<SFpnt[j]->GetTBName()<<" has fired!"<<endl;
	SFfired_ran[j-8]=1;	
      }
      //cout<<"Firedmap["<<word<<"]: "<<Firedmap[word]<<endl;
    } //loop over Effdet-vector
    int allSF=1;
    int allSFp=0;
    int SASplanes=0;
    for(int i=0;i<nplanes;i++){ 
      if(Ratrig1[2]!=NULL && SFfired_ran[i]!=0)Ratrig1[2]->Fill(i);
      allSF*=SFfired_ran[i];
      if(i<11 || !randomSFnumber_SAS) allSFp+=SFfired_ran[i];
      //if randomSFnumber_SAS != 0 planes after SM2 are treated separetly
      else SASplanes+=SFfired_ran[i];
    }
    if(Ratrig1[2]!=NULL && allSF!=0)Ratrig1[2]->Fill(19);
    if(Ratrig1[3]!=NULL)Ratrig1[3]->Fill(allSFp);
    /*if(allSFp==nplanes-1){
      cout<<"allSF: "<<allSF<<endl;
      allSF=1;
      for(int i=0;i<nplanes-1;i++){ 
	cout<<"SFfired_ran["<<i<<"]"<<SFfired_ran[i]<<endl;
	allSF*=SFfired_ran[i];
	cout<<"allSF: "<<allSF<<endl;
      }
    }*/
    //extrapolate to -350 mm if this track has momentum:
    CsHelix exHelTarget;
    //CsHelix exHelTarget1;
    //CsHelix exHelTarget2;
    if(mom!=0){
      //record positions at both ends of the target:
      //cout<<"Extrapolating..."<<endl;
      starthel->Extrapolate(300,exHelTarget);
      downX=exHelTarget.getX();
      downY=exHelTarget.getY();
      //cout<<"Extrapolating..."<<endl;
      starthel->Extrapolate(-1000,exHelTarget);
      upX=exHelTarget.getX();
      upY=exHelTarget.getY();
      //cout<<"Extrapolating..."<<endl;
      starthel->Extrapolate(-350,exHelTarget);
      //plot allSFp vs dx and dy:
      if(Ratrig2[4]!=NULL)Ratrig2[4]->Fill(exHelTarget.getDXDZ(),allSFp);
      if(Ratrig2[5]!=NULL)Ratrig2[5]->Fill(exHelTarget.getDYDZ(),allSFp);
      for(int i=0;i<nplanes;i++){ 
	if(SFfired_ran[i]!=0 && Ratrig2[7]!=NULL)Ratrig2[7]->Fill(i,allSFp);
      }
      //cout<<"Done Extrapolating..."<<endl;
    }    
    //check if enough planes fired before and after SM2
    if(allSFp>=randomSFnumber-randomSFnumber_SAS && SASplanes >= randomSFnumber_SAS){
      //cout<<"treating allSFp>=11"<<endl;
      //if(mom>GoodTrackMomCut && !SFfired_ran[11] && !SFfired_ran[12]){
      if(mom>GoodTrackMomCut){
	//cout<<"SF7 is not participating..."<<endl;
	//check whether there are time correlated hits in SF07 X&Y: (19 & 20)
	bool hasInTimeHit[2];
	for(int i=0;i<2;i++) hasInTimeHit[i]=false;
	int histCnt=0;
	for(int i=19;i<=20;i++){
	  double clTime;
	  if(SFpnt[i]!=0){
	    list<CsCluster*>clList = SFpnt[i]->getMyClusters();
	    list<CsCluster*>::iterator Iclu;
	    if(!clList.size()) continue;
	    for( Iclu = clList.begin(); Iclu != clList.end(); Iclu ++ ) {
	      if((*Iclu)->getTime(clTime)){
		if(Ratrig1[18+histCnt]!=NULL)Ratrig1[18+histCnt]->Fill(clTime-(*Ti)->getMeanTime());
		if(fabs(clTime-(*Ti)->getMeanTime())<2.0) hasInTimeHit[histCnt]=true;
	      }
	    }//loop over clusters
	  }
	  histCnt++;
	}//loop over SF07X/Y
	//
	//
	if(hasInTimeHit[0] && hasInTimeHit[1]){ hasGoodTrack=true; goodTr = true;}
	//hasGoodTrack=true;
	//cout<<"Applying cut..."<<endl;
	if(RandomTargetCut!=0){
	  if(sqrt(downX*downX+downY*downY)>RandomTargetCut || sqrt(upX*upX+upY*upY)>RandomTargetCut) {
	      hasGoodTrack=false;
	      goodTr = false;
	  }
	}
	//cout<<"Done Applying cut..."<<endl;
	nGoodTracks++;
	GoodTrackTime=(*Ti)->getMeanTime();
	GoodTrackChi2=(*Ti)->getChi2()/((*Ti)->getNumberOfAssociatedClusters()-5);
	GoodTrackMom=mom;
	GoodTrackTgtX=exHelTarget.getX();
	GoodTrackTgtY=exHelTarget.getY();
	GoodTrackTgtDx=exHelTarget.getDXDZ();
	GoodTrackTgtDy=exHelTarget.getDYDZ();
	GoodTrackStartHelix=*starthel;

	if(goodTr) {
	    //Save good tracks helix
	    goodSFhelixes.push_back(*starthel);
	}
	/*else if(beamAtBMS != NULL) {
	    delete beamAtBMS;
	    beamAtBMS = NULL;
	}*/
	//starthel->Print();
	for(int i=0;i<nplanes;i++){ 
	  if(SFfired_ran[i]!=0 && Ratrig2[9]!=NULL)Ratrig2[9]->Fill(i,allSFp);
	}
      }
      else{
	//cout<<"is this halo cand?"<<endl;
	hasHaloCand=true;
	HaloTrackTime=(*Ti)->getMeanTime();
	//cout<<"Done with halo cand"<<endl;
      }
    }
    if(allSFp==nplanes-2 && Ratrig1[21]!=NULL)Ratrig1[21]->Fill((*Ti)->getChi2()/((*Ti)->getNumberOfAssociatedClusters()-5));
  }//loop over tracks
  
  if(hasHaloCand){
    //cout<<"Calling SFTrigger in BMS..."<<endl;
    //BMS->SFTrigger(HaloTrackTime,3,0);  
    //cout<<"Done with SFTrigger..."<<endl;
  }
  //record number of tracks with momentum for "good" events:
  if(hasGoodTrack){
    if(Ratrig1[4]!=NULL) Ratrig1[4]->Fill(ntrackMom);
    if(Ratrig1[5]!=NULL) Ratrig1[5]->Fill(GoodTrackTime);
    if(Ratrig1[13]!=NULL) Ratrig1[13]->Fill(GoodTrackMom);
    if(Ratrig1[14]!=NULL) Ratrig1[14]->Fill(nGoodTracks);
    if(Ratrig1[16]!=NULL) Ratrig1[16]->Fill(GoodTrackTgtDx);
    if(Ratrig1[17]!=NULL) Ratrig1[17]->Fill(GoodTrackTgtDy);
    if(Ratrig1[20]!=NULL) Ratrig1[20]->Fill(GoodTrackChi2);
    if(Ratrig2[3]!=NULL) Ratrig2[3]->Fill(GoodTrackTgtX,GoodTrackTgtY);
  }
  //cout<<"returning"<<endl;
  return hasGoodTrack;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool CsBeamRecons:: findTriggerMask(int year){
  bool oneTrigger=false;
  IT=MT=LT=OT=IMT=VT=VTO=VTI=CT=false;
  Trigger[0]=Trigger[1]=Trigger[2]=Trigger[3]=Trigger[4]=Trigger[5]=false;
  unsigned trigger_mask;
  trigger_mask=CsEvent::Instance()->getTriggerMask();

  if(year==2001){
    //0 0000 0000 - 0 :Empty
    //0 0000 0001 - 1 :Empty
    //0 0000 0010 - 2 :Inner Trigger  
    //0 0000 0100 - 4 :Middle Trigger 
    //0 0000 1000 - 8 :Ladder Trigger 
    int triggerMult=0;
    if( trigger_mask&2 ){ 
      IT=true; 
      triggerMult++;
    }
    if( trigger_mask&4 ){ 
      MT=true; 
      triggerMult++;
    }
    if( trigger_mask&8 ){ 
      LT=true; 
      triggerMult++;
    }
    if( trigger_mask&16 ){ 
      VT=true; 
      //triggerMult++;
    }
    if( trigger_mask&32 ){ 
      CT=true; 
      //triggerMult++;
    }
    if( trigger_mask&64 ){ 
      //RT=true; 
      //triggerMult++;
    }
    if( trigger_mask&128 ){ 
      //MMT=true; 
      //triggerMult++;
    }
    if( trigger_mask&256 ){ 
      //BT=true; 
      //triggerMult++;
    }
    if(triggerMult==1){ 
      oneTrigger=true;
      Trigger[0]=IT;
      Trigger[1]=MT;
      Trigger[2]=LT;
      Trigger[3]=OT;
      Trigger[4]=CT;
      Trigger[5]=IMT;
    }
  }
  if(year==2002 || year==2003 || year==2004){
    //0000 0000 0000 -   0 :Empty
    //0000 0000 0001 -   1 :Inner Trigger (ITCV)
    //0000 0000 0010 -   2 :Middle Trigger (MTCV)   
    //0000 0000 0100 -   4 :Ladder Trigger (LTCV)
    //0000 0000 1000 -   8 :Outer Trigger (OTV)
    //0000 0001 0000 -  16 :Calorimeter
    //0000 0010 0000 -  32 :Veto Inner
    //0000 0100 0000 -  64 :Veto Outer
    //0000 1000 0000 - 128 :Beam Trigger
    //0001 0000 0000 - 256 :incl MT (MTV)
    //0010 0000 0000 - 512 :unused
    //0100 0000 0000 -1024 :unused
    //1000 0000 0000 -2048 :Random Trigger
    int triggerMult=0;
    //cout<<"Trigger-Mask: "<<trigger_mask<<endl;
    if( trigger_mask&1 ){ 
      IT=true; 
      triggerMult++;
      //cout << "Inner Trigger"<<endl;
    }
    if( trigger_mask&2 ){ 
      MT=true; 
      triggerMult++;
      //cout << "Middle Trigger"<<endl;
    }
    if( trigger_mask&4 ){ 
      LT=true; 
      triggerMult++;
      //cout << "Ladder Trigger"<<endl;
    }
    if( trigger_mask&8 ){ 
      OT=true; 
      triggerMult++;
      //cout << "Outer Trigger"<<endl;
    }
    if( trigger_mask&16 ){ 
      CT=true; 
      triggerMult++;
      //cout << "Veto Trigger"<<endl;
    }
    if( trigger_mask&32 ){ 
      VTI=true; 
      //triggerMult++;
      //cout << "Calorimeter Trigger"<<endl;
    }
    if( trigger_mask&64 ){ 
      VTO=true; 
      //triggerMult++;
      //cout << "Random Trigger"<<endl;
    }
    if( trigger_mask&128 ){ 
      BT=true; 
      //triggerMult++;
      //cout << "Beam Trigger"<<endl;
    }
    if( trigger_mask&256 ){ 
      IMT=true; 
      triggerMult++;
      //cout << "inclusive Middle Trigger"<<endl;
    }
    
    if( trigger_mask&2048 ){ 
      //RT=true; 
      //triggerMult++;
      //cout << "Random Trigger"<<endl;
    }
    //printf("Trigger Mask: %08X \n",trigger_mask);
    //if( trigger_mask&(unsigned)(pow(2,30)))Filter_Flag=1;
    if( trigger_mask&(1<<30))Filter_Flag=1;
    else Filter_Flag=0;
    if(triggerMult==1){
      oneTrigger=true;
    }
    Trigger[0]=IT;
    Trigger[1]=MT;
    Trigger[2]=LT;
    Trigger[3]=OT;
    Trigger[4]=CT;
    Trigger[5]=IMT;
  }
  return oneTrigger;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBeamRecons:: FiredSF()
{
  int nfiredSF=0;
  int nhit=0;
  for(int deti=0; deti<4; deti++) {
    list<CsDigit*>detdigit = SFpnt[deti]->getMyDigits();
    if(!detdigit.size()) continue;
    list<CsDigit*>::iterator Idetdig;
    nhit=0;
    for(Idetdig = detdigit.begin();Idetdig != detdigit.end(); Idetdig ++ )
      {
	if(fabs((*Idetdig)->getDatum())<1.5) nhit++;
	mD1[24]->Fill((*Idetdig)->getDatum());
	
      } //loop over digits     
    if(nhit!=0) nfiredSF++;
  }// loop over detectors
  return nfiredSF;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool CsBeamRecons:: checktrackorientation(CsTrack* track)
{
  CsHelix  endhel;
  double thetaX,thetaY;
  bool ortOK=true;
  vector<CsHelix>trhel =track->getHelices(); 
  vector<CsHelix>::iterator Hi;
  double endz=-99999999;
  for(Hi = trhel.begin();Hi != trhel.end(); Hi++ ){
    if( (Hi)->getZ()>endz ){
      endz=(Hi)->getZ();
      endhel=*Hi;
    }
  }
  //get Angles
  thetaX=atan(endhel.getDXDZ());
  thetaY=atan(endhel.getDYDZ());
  if(fabs(thetaX)>=thetaXcut || fabs(thetaY)>=thetaYcut) ortOK=false;
  return ortOK;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool CsBeamRecons:: checktracktiming(CsTrack* track,int &nSFCl,int ntracks[])
{
  bool isgoodtrack=true;
  bool isFI=false;
  double clustertime;
  double meantime=9999;
  vector<double> clTimes;
  clTimes.clear();
  nSFCl=0;
  if(track->hasMeanTime()) meantime=track->getMeanTime();
  list<CsCluster*> BTcluster = track->getClusters();
  list<CsCluster*>::iterator Ci;
  for(Ci = BTcluster.begin();Ci != BTcluster.end(); Ci++ ){
    isFI=false;
    //See wether it is a good track:
    if((*Ci)->getW()<0){
      (*Ci)->getTime(clustertime);
      if(fabs(clustertime)>MaxSFHitTime) isgoodtrack=false;
      if(meantime!=9999 && mD1[28]!=NULL) mD1[28]->Fill(meantime-clustertime);
    }

    //
    //check which detectors are involved in this track:
    //

    list<CsDetector*> det = (*Ci)->getDetsList();
    if(det.empty()) continue;
    //list<CsDetector*>::iterator idet=det.begin(); 
    list<CsDetector*>::iterator idet;
    for(idet = det.begin();idet != det.end(); idet++ ){
      //cout<<"TBName: "<<(*idet)->GetTBName()<<": "<<clustertime<<" "<<BTcluster.size()<<endl;
      if((*idet)->GetTBName()=="FI01X1__") {
	if(mH2[107]!=NULL) mH2[107]->Fill(0,BTcluster.size());
	nSFCl++;
	isFI=true;
      }
      if((*idet)->GetTBName()=="FI01Y1__") {
	if(mH2[107]!=NULL) mH2[107]->Fill(1,BTcluster.size());
	nSFCl++;
	isFI=true;
      }
      if((*idet)->GetTBName()=="FI02X1__") {
	if(mH2[107]!=NULL) mH2[107]->Fill(2,BTcluster.size());
	nSFCl++;
	isFI=true;
      }
      if((*idet)->GetTBName()=="FI02Y1__") {
	if(mH2[107]!=NULL) mH2[107]->Fill(3,BTcluster.size());
	nSFCl++;
	isFI=true;
      }
      if((*idet)->GetTBName()=="FI15X1__") {
	if(mH2[107]!=NULL) mH2[107]->Fill(4,BTcluster.size());
	nSFCl++;
	isFI=true;
      }
      if((*idet)->GetTBName()=="FI15Y1__") {
	if(mH2[107]!=NULL) mH2[107]->Fill(5,BTcluster.size());
	nSFCl++;
	isFI=true;
      }
      if((*idet)->GetTBName()=="SI01X1__" && mH2[107]!=NULL) mH2[107]->Fill(6,BTcluster.size());
      if((*idet)->GetTBName()=="SI01Y1__" && mH2[107]!=NULL) mH2[107]->Fill(7,BTcluster.size());
      if((*idet)->GetTBName()=="SI01U1__" && mH2[107]!=NULL) mH2[107]->Fill(8,BTcluster.size());
      if((*idet)->GetTBName()=="SI01V1__" && mH2[107]!=NULL) mH2[107]->Fill(9,BTcluster.size());
      if((*idet)->GetTBName()=="SI02X1__" && mH2[107]!=NULL) mH2[107]->Fill(10,BTcluster.size());
      if((*idet)->GetTBName()=="SI02Y1__" && mH2[107]!=NULL) mH2[107]->Fill(11,BTcluster.size());
      if((*idet)->GetTBName()=="SI02U1__" && mH2[107]!=NULL) mH2[107]->Fill(12,BTcluster.size());
      if((*idet)->GetTBName()=="SI02V1__" && mH2[107]!=NULL) mH2[107]->Fill(13,BTcluster.size());
    }
    //be careful since SI give clusters with time now!
    if(isFI && (*Ci)->getTime(clustertime)) clTimes.push_back(clustertime);
  } // loop over clusters
  //record nSFCl vs BTcluster.size():
  if(mH2[130]!=NULL) mH2[130]->Fill(nSFCl,BTcluster.size());

  if(nSFCl<=6) ntracks[nSFCl]++;
  if(track->hasMeanTime()) if(fabs(track->getMeanTime())<3.06)ntracks[7]++;

  for(unsigned int i=0;i<clTimes.size();i++){
    for(unsigned int j=i+1;j<clTimes.size();j++){
      if(mH1[103]!=NULL)mH1[103]->Fill(clTimes[i]-clTimes[j]);
    }
  }
  return isgoodtrack;
}
//==========================================================================
int CsBeamRecons:: CheckVeto(int mode, double RefTime)
{
  //mode==0: Check for cleaned trigger relative to reference time
  //mode==1: just like (0) for goodTrigger with at least one SF fired
  //mode==2: just like (0) but for events with !goodTrSFok
  //What kind of vetos do we have:
  //VT01X1_O: Veto Outer 1  
  //Oben
  //       18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 
  //        |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
  //Saleve  |  |  |  |  |  |  |  |        |  |  |  |  |  |  |  |  Jura
  //        |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
  //        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
  //                        Unten
  //Channel 8,9,26,27 not meantimed
  //VT01X1mO: meantimer of Veto Outer 1 [16 Channels !]
  //
  //VT01P1_I: Veto inner 1  
  //VT02P1_I: Veto inner 2 
  //VT01P1bl: Beamline-Veto  
  //VT01P1sf: Zu beginn des Runs die zu 16-ner Gruppen zusammengefassten
  //          Horizontalen (Y-)Fasern von FI01, spaeter die Vertikalen (X-) von FI02 
  //VT02P1sf: Die zu 16-er Gruppen zusammengefassten Horizontalen (Y-)-Fasern von FI02 
  //VTsum:    77  0 770  16 1 16   0 15 1  
  //VT01X1dO:  
  //VT01X1uO:  
  //
  //get VI1,VI2 and VO1 detectors:
  static bool first(true);
  static CsDetector*  VTpnt[4];
  static std::string VTname[4];
  double time;
  int channel;
  if(first){
    cout<<"CheckVeto: Setting things up..."<<endl;
    VTname[0]="VT01P1_I";
    VTname[1]="VT02P1_I";
    VTname[2]="VT01X1_O";
    VTname[3]="VT01X1mO";
    list <CsDetector*>  mydets = CsGeom::Instance()->getDetectors();
    list<CsDetector*>::iterator imydet;
    for(int i=0; i<4; i++) {
      VTpnt[i]=0;
      for( imydet = mydets.begin(); imydet != mydets.end(); imydet++ ) {
	//cout<<(*imydet)->GetTBName()<<endl;
	if((*imydet)->GetTBName()==VTname[i]) {
	  cout<<"CsBeamRecons: "<<VTname[i]<<" was found!"<<endl;
	  VTpnt[i]=(*imydet);
	}
	//cout<<"CsBeamRecons: "<<VTname[i]<<" was not found!"<<endl;
      }
    }
    first=false;
  }
  //loop over VI1 and VI2 and VO1 channels:
  bool VTfired[4];
  for(int k=0;k<4;k++){
    VTfired[k]=false;
    if(VTpnt[k]!=0){
      //cout<<"CheckVeto: "<<VTname[k]<<" was found!"<<endl;
      list<CsDigit*>digit = VTpnt[k]->getMyDigits();
      if(!digit.size()) continue;  
      list<CsDigit*>::iterator Idig;
      for( Idig = digit.begin(); Idig != digit.end(); Idig ++ ) {
	time  = (*Idig)->getDatum();
	channel=(*Idig)->getAddress();
	//cout<<VTname[k]<<" "<<channel<<" "<<time<<endl;
	if(mode==0){
	  //record timing of 4 Channels of VI1 and VI2 and 20 Channels [16 + 4] of Outer1
	  if(k==0){
	    if(mH2[131]!=NULL)mH2[131]->Fill(time,channel);
	    if(mH2[135]!=NULL)mH2[135]->Fill(time,channel);
	  }
	  if(k==1){
	    if(mH2[132]!=NULL)mH2[132]->Fill(time,channel);
	    if(mH2[136]!=NULL)mH2[136]->Fill(time,channel);
	  }
	  if(k==2){
	    if(mH2[133]!=NULL){
	      if(channel==8)mH2[133]->Fill(time,channel-8);
	      if(channel==9)mH2[133]->Fill(time,channel-8);
	      if(channel==26)mH2[133]->Fill(time,channel-24);
	      if(channel==27)mH2[133]->Fill(time,channel-24);
	    }
	    if(mH2[137]!=NULL){
	      if(channel==8)mH2[137]->Fill(time,channel-8);
	      if(channel==9)mH2[137]->Fill(time,channel-8);
	      if(channel==26)mH2[137]->Fill(time,channel-24);
	      if(channel==27)mH2[137]->Fill(time,channel-24);
	    }
	  }
	  if(k==3){
	    if(mH2[134]!=NULL)mH2[134]->Fill(time,channel);
	    if(mH2[138]!=NULL)mH2[138]->Fill(time,channel);
	  }
	}//mode==0
	double TimeCut_VT_min=TimeCut_VT_cntr[k]-TimeCutVt[k];
	double TimeCut_VT_max=TimeCut_VT_cntr[k]+TimeCutVt[k];
	if(time>TimeCut_VT_min && time<TimeCut_VT_max){
	  VTfired[k]=true;
	}
	if(RecMethod==3){//RecMethod 3: every hit counts
	  if(time>EffTimeCutJ_min && time<EffTimeCutJ_max){
	    //VTfired[i]=true;
	  }
	}
      }//loop over digits
    }
    
  }//loop over vetos
  int nfiredVT=0;
  for(int k=0;k<4;k++){
    if(VTfired[k])nfiredVT++;
  }
  if(mode==1){//at least one of the SF stations has fired:
    //let's see what the vetos have:
    for(int k=0;k<4;k++){
      if(VTfired[k] && Ratrig2[13]!=NULL)Ratrig2[13]->Fill(k,nfiredVT);
    }
  }
  if(mode==2){//none of the SF stations have fired:
    //let's see what the vetos have:
    for(int k=0;k<4;k++){
      if(VTfired[k] && Ratrig2[12]!=NULL)Ratrig2[12]->Fill(k,nfiredVT);
    }
  }
  return true;
}
//==========================================================================
int CsBeamRecons:: SFEfficiency(int mode, double RefTime)
{
  //mode==0: "normal" efficiencies
  //mode==1: efficiencies for cleaned trigger
  //mode==2: efficiencies for cleaned trigger relative to reference time
  //mode==3: just like (2) but for events with !goodTrSFok
  //mode==4: just like (2) but for events with !goodTrigger
  //let's see which of the SF has fired:
  double time;
  int nfiredSF=0;
  int nhit=0;
  SFfired[0]=SFfired[1]=SFfired[2]=SFfired[3]=false;
  SFMult[0]=SFMult[1]=SFMult[2]=SFMult[3]=0;
  bool SFfiredR[6];
  bool StationOk[4];//to find which stations have hits for ncut==0 && !goodTrSFok
  for(int i=0;i<6;i++){
    SFfiredR[i]=0;
    if(i<4)StationOk[i]=false;
  }
  int mult=0;
  if(mode<2){
    for(int i=0;i<4;i++){
      mult=0;
      if(SFpnt[i]!=0){
	list<CsDigit*>digit = SFpnt[i]->getMyDigits();
	if(!digit.size()) continue;  
	list<CsDigit*>::iterator Idig;
	nhit=0;
	for( Idig = digit.begin(); Idig != digit.end(); Idig ++ ) {
	  time  = (*Idig)->getDatum();
	  //if(fabs(time)<1.5) SFfired[i]=true;  
	  if(fabs(time)<MaxSFHitTime){ //RecMethod 0: only in-time hits count
	    SFfired[i]=true; 
	    mult++;
	  }
	  if(RecMethod==3){//RecMethod 3: every hit counts
	    if(time>EffTimeCutJ_min && time<EffTimeCutJ_max){
	      SFfired[i]=true;
	      mult++;
	    }
	  }
	}//loop over digits
      }
      if(SFfired[i]) nfiredSF++;
      SFMult[i]=mult;
    }//loop over SF
    //look at time-difference between all SF-clusters:
    double time1,time2;
    if(Htrig1[18]!=NULL){
      for(int i=0;i<6;i++){
	if(UpTpnt[i]!=0){
	  list<CsCluster*> clusterL=UpTpnt[i]->getMyClusters();
	  if(!clusterL.size()) continue;  
	  list<CsCluster*>::iterator Icl;
	  for( Icl = clusterL.begin(); Icl != clusterL.end(); Icl ++ ) { 
	    if((*Icl)->getTime(time1)){
	      for(int j=i+1;j<6;j++){
		if(UpTpnt[j]!=0){
		  list<CsCluster*> clusterL2=UpTpnt[j]->getMyClusters();
		  if(!clusterL2.size()) continue;  
		  list<CsCluster*>::iterator Icl2;
		  for( Icl2 = clusterL2.begin(); Icl2 != clusterL2.end(); Icl2 ++ ){ 
		    if((*Icl2)->getTime(time2)){
		      Htrig1[18]->Fill(time1-time2);
		    }
		  }//loop over cluster of SFplane 2
		}
	      }//loop over SFplanes 2
	    }
	  }//loop over cluster of SFplane 1
	}
      }//loop over SF-planes 1
    }//histo check
  }//mode<2
  else{//mode==2 || 3 || 4
    for(int i=0;i<6;i++){
      mult=0;
      if(UpTpnt[i]!=0){
	list<CsDigit*>digit = UpTpnt[i]->getMyDigits();
	if(!digit.size()) continue;  
	list<CsDigit*>::iterator Idig;
	for( Idig = digit.begin(); Idig != digit.end(); Idig ++ ) {
	  time  = (*Idig)->getDatum();
	  if(mode==2 && Htrig2[12]!=NULL)Htrig2[12]->Fill(time-RefTime,i);
	  if(mode==3 && Htrig2[14]!=NULL)Htrig2[14]->Fill(time-RefTime,i);
	  if(mode==4 && Htrig2[15]!=NULL)Htrig2[15]->Fill(time-RefTime,i);
	  if(i<4){
	    //if(fabs(time-RefTime)<EffTimeCutJ)SFfiredR[i]=true; 
	    if((time-RefTime)>EffTimeCutJ_min && (time-RefTime)<EffTimeCutJ_max)SFfiredR[i]=true; 
	    if(fabs(time-RefTime)<MxTmdiff)StationOk[i]=true;
	  }
	  else{
	   if(fabs(time-RefTime)<EffTimeCut15)SFfiredR[i]=true; 
	  }
	}//loop over digits
	//look
      }
    }//loop over SF
    //Fill histos for !goodTrSFok-statistics:
    int nfiredSpec=0;
    for(int i=0;i<4;i++)if(StationOk[i])nfiredSpec++;
    if(Ratrig1[26]!=NULL && mode==3)Ratrig1[26]->Fill(nfiredSpec);
    if(Ratrig2[10]!=NULL && mode==3){
      for(int i=0;i<4;i++)if(StationOk[i])Ratrig2[10]->Fill(i,nfiredSpec);
    }
    nfiredSF=nfiredSpec;
    if(mode==3 && nfiredSF==0){
      if(Ratrig1[31]!=NULL)Ratrig1[31]->Fill(RefTime);

    }
  }//mode>=2
  
  
  //Fill histos for Efficiency
  if(mode==0){
    //for all planes:
    if((SFfired[0]&&SFfired[1])&&(SFfired[2]&&SFfired[3])){
      if(mH1[26]!=NULL)mH1[26]->Fill(0);
    }
    if(mH1[26]!=NULL){
      //planes in front of the target:
      //for plane 0:
      if(SFfired[1]&&(SFfired[2]&&SFfired[3])) mH1[26]->Fill(1); 
      //for plane 1:
      if(SFfired[0]&&(SFfired[2]&&SFfired[3])) mH1[26]->Fill(2);
      //for plane 2:
      if(SFfired[0]&&(SFfired[1]&&SFfired[3])) mH1[26]->Fill(3);
      //for plane 3:
      if(SFfired[0]&&(SFfired[1]&&SFfired[2])) mH1[26]->Fill(4);
    }
  }
  if(mode==1){
    if(Htrig1[6]!=NULL){
      //for all planes:
      if((SFfired[0]&&SFfired[1])&&(SFfired[2]&&SFfired[3]))Htrig1[6]->Fill(0);
      //planes in front of the target:
      //for plane 0:
      if(SFfired[1]&&(SFfired[2]&&SFfired[3])) Htrig1[6]->Fill(1); 
      //for plane 1:
      if(SFfired[0]&&(SFfired[2]&&SFfired[3])) Htrig1[6]->Fill(2);
      //for plane 2:
      if(SFfired[0]&&(SFfired[1]&&SFfired[3])) Htrig1[6]->Fill(3);
      //for plane 3:
      if(SFfired[0]&&(SFfired[1]&&SFfired[2])) Htrig1[6]->Fill(4);
    }
  }
  if(mode==2){
    if(Htrig1[16]!=NULL){
      //for all planes:
      if((SFfiredR[0]&&SFfiredR[1])&&(SFfiredR[2]&&SFfiredR[3]))Htrig1[16]->Fill(0);
      //planes in front of the target:
      //for plane 0:
      if(SFfiredR[1]&&(SFfiredR[2]&&SFfiredR[3])) Htrig1[16]->Fill(1); 
      //for plane 1:
      if(SFfiredR[0]&&(SFfiredR[2]&&SFfiredR[3])) Htrig1[16]->Fill(2);
      //for plane 2:
      if(SFfiredR[0]&&(SFfiredR[1]&&SFfiredR[3])) Htrig1[16]->Fill(3);
      //for plane 3:
      if(SFfiredR[0]&&(SFfiredR[1]&&SFfiredR[2])) Htrig1[16]->Fill(4);
      //normalisation for FI15:
      if((SFfiredR[1]&&SFfiredR[3])&&(SFfiredR[4]&&SFfiredR[5]))Htrig1[16]->Fill(5);
      //FI15X
      if(SFfiredR[1]&&(SFfiredR[3]&&SFfiredR[5])) Htrig1[16]->Fill(6);
      //FI15Y
      if(SFfiredR[1]&&(SFfiredR[3]&&SFfiredR[4])) Htrig1[16]->Fill(7);
    }
    //requiring all 6 planes:
    if(Htrig1[17]!=NULL){
      //for all planes:
      if((SFfiredR[0]&&SFfiredR[1])&&(SFfiredR[2]&&SFfiredR[3])&&(SFfiredR[4]&&SFfiredR[5])) Htrig1[17]->Fill(0);
      //for FI01X1:
      if(SFfiredR[1]&&(SFfiredR[2]&&SFfiredR[3])&&(SFfiredR[4]&&SFfiredR[5])) Htrig1[17]->Fill(1);
      //for FI01Y1:
      if(SFfiredR[0]&&(SFfiredR[2]&&SFfiredR[3])&&(SFfiredR[4]&&SFfiredR[5])) Htrig1[17]->Fill(2);
      //for FI02X1:
      if((SFfiredR[0]&&SFfiredR[1])&&SFfiredR[3]&&(SFfiredR[4]&&SFfiredR[5])) Htrig1[17]->Fill(3);
      //for FI02Y1:
      if((SFfiredR[0]&&SFfiredR[1])&&SFfiredR[2]&&(SFfiredR[4]&&SFfiredR[5])) Htrig1[17]->Fill(4);
      //for FI15X1:
      if((SFfiredR[0]&&SFfiredR[1])&&(SFfiredR[2]&&SFfiredR[3])&&SFfiredR[5]) Htrig1[17]->Fill(5);
      //for FI15Y1:
      if((SFfiredR[0]&&SFfiredR[1])&&(SFfiredR[2]&&SFfiredR[3])&&SFfiredR[4]) Htrig1[17]->Fill(6);
    }
  }
  return nfiredSF;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBeamRecons:: SFTrigger(void)
{
  //designed to clean the trigger with the help of SF stations 03/04, 
  //by requiring a certain number of hits.
  //let's see which of the SF has fired:
  double time;
  int nfiredSF=0;
  SFfired[4]=SFfired[5]=SFfired[6]=SFfired[7]=SFfired[8]=SFfired[9]=false;
  for(int i=8;i<=13;i++){
    if(SFpnt[i]!=0){
      list<CsDigit*>digit = SFpnt[i]->getMyDigits();
      if(!digit.size()) continue;  
      list<CsDigit*>::iterator Idig;
      for( Idig = digit.begin(); Idig != digit.end(); Idig ++ ) {
	time  = (*Idig)->getDatum();
	//if(fabs(time)<1.5) SFfired[i]=true;  
	if(RecMethod==0){
	  if(fabs(time)<MaxSFHitTime){ //RecMethod 0: only in-time hits count
	    SFfired[i-4]=true; 
	  }
	}
	//if(RecMethod==3){
	else {
	  if(fabs(time)<1.57){ //RecMethod 3: use generic width of distribution
	    SFfired[i-4]=true;
	  }
	}
      }//loop over digits
    }
    if(SFfired[i-4]) nfiredSF++;
  }//loop over SF
  //Fill histos for Efficiency
  if( ((SFfired[4]&&SFfired[5])&&(SFfired[6]&&SFfired[7]))&&(SFfired[8]&&SFfired[9])){  
    if(mH1[104]!=NULL)mH1[104]->Fill(0);
  }
  if(mH1[104]!=NULL){
    //for plane 3x:
    if( ((SFfired[5]&&SFfired[6]) && (SFfired[7]&&SFfired[8])) && SFfired[9]) mH1[104]->Fill(1); 
    //for plane 3y:
    if( ((SFfired[4]&&SFfired[6]) && (SFfired[7]&&SFfired[8])) &&SFfired[9]) mH1[104]->Fill(2);
    //for plane 3u:
    if( ((SFfired[4]&&SFfired[5]) && (SFfired[7]&&SFfired[8])) &&SFfired[9]) mH1[104]->Fill(3);
    //for plane 4x:
    if( ((SFfired[4]&&SFfired[5]) && (SFfired[6]&&SFfired[8])) &&SFfired[9]) mH1[104]->Fill(4);
    //for plane 4y:
    if( ((SFfired[4]&&SFfired[5]) && (SFfired[6]&&SFfired[7])) &&SFfired[9]) mH1[104]->Fill(5);
    //for plane 4u:
    if( ((SFfired[4]&&SFfired[5]) && (SFfired[6]&&SFfired[7])) &&SFfired[8]) mH1[104]->Fill(6);
  }
  return nfiredSF;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons:: makeWiremap()
{
  //let's see which of the SF has fired:
  double wire;
  int ii=73;
  for(int i=0;i<4;i++){
    if(SFpnt[i]!=0){
      list<CsDigit*>digit = SFpnt[i]->getMyDigits();
      if(!digit.size()) continue;  
      list<CsDigit*>::iterator Idig;
      for( Idig = digit.begin(); Idig != digit.end(); Idig ++ ) {
	wire  = (*Idig)->getAddress();
	if(mH2[ii]!=NULL)mH2[ii]->Fill(wire,SFMult[i]);
      }//loop over digits
      if(i==2){
	list<CsCluster*>clusterL = SFpnt[i]->getMyClusters();
	list<CsCluster*>::iterator Iclu;
	for( Iclu = clusterL.begin(); Iclu != clusterL.end(); Iclu ++ ) {
	  wire  = (*Iclu)->getU();
	  if(mH1[92]!=NULL)mH1[92]->Fill(wire);
	}//loop over cluster
      }
      if(i==3){
	list<CsCluster*>clusterL = SFpnt[i]->getMyClusters();
	list<CsCluster*>::iterator Iclu;
	for( Iclu = clusterL.begin(); Iclu != clusterL.end(); Iclu ++ ) {
	  wire  = (*Iclu)->getU();
	  if(mH1[93]!=NULL)mH1[93]->Fill(wire);
	}//loop over cluster
      }
    }
    ii++;
  }//loop over SF
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons:: tmcorrel(int nbmhod,int nbms,int & ntrack)
{
        int i, j, ntot, len;
        double delt=0;
        double x, y, z, dxdz, dydz, Htime, Htmchiq;
        double Bp, Btime, Bpchiq, Btmchiq;
        double Bm0, Bm1, Hm0, Hm1;
        double TrackTime, TrackChi2;
        int Bnhit, Bnhittot, Hnhit, Htot;
        bool Htrig;
        int ibms, ihod;
//
        ntrack=0;
        len=nbmhod*nbms;
//    
       int (*npos)[2] = new int [len][2];
       double *dlttm  = new double [len]; 
       for(ntot=0, i=0; i<nbms; i++){
             for(j=0; j<nbmhod; j++)  {
                 delt=BMS->getBMStm(i)-BmHod.getBmFiHodtm(j);
                 HiTBThod->Fill(BmHod.getBmFiHodtm(j), BMS->getBMStm(i));
                 HiTmdiff->Fill(delt);
                 delt=fabs(delt);
                 if(delt<MxTmdiff){
                      dlttm[ntot]=delt;
                      npos[ntot][0]=i;
                      npos[ntot][1]=j;
                      ntot++;
                 }
             }
        }
        if((nbms==1)&&(nbmhod==1)) {
                 HiTmdif1->Fill(delt); 
                 HiTBTh1->Fill(BMS->getBMStm(0),BmHod.getBmFiHodtm(0));
        }
        if(!ntot) {
             delete [] npos;
             delete [] dlttm; 
             return;
        }
//
//       select min delttm
//
        double mindelt;
        int jmin=0;
        int nleft=ntot;
        while(nleft){
           for(mindelt=10.e5, i=0; i<ntot; i++) {
                 if((mindelt>dlttm[i])&&(npos[i][0]>=0)) 
                          {jmin=i; mindelt=dlttm[i];}
           }
//
//    save track in CsBeam
// 
      ibms=npos[jmin][0]; 
      int BMSspec;
      BMS->getBMStrack(ibms, Bp, Bpchiq, Btime, Btmchiq, Bnhit, Bnhittot,BMSspec);
      BMS->getBMSmoments(ibms, Bm0, Bm1);
//
      ihod=npos[jmin][1];
      BmHod.getBmFiHodtrack(ihod, x, y, z, dxdz, dydz, Hnhit, Htot, Htrig, Htime, Htmchiq);
      BmHod.getBmHodmoments(ihod, Hm0, Hm1);
//
      TrackTime = (Bm1+ Hm1)/(Bm0+Hm0);
      TrackChi2 = BmHod.getNewChi2(ihod, TrackTime) + BMS->getNewChi2(ibms, TrackTime);
      FillCovMx(BmResol*BmChrg/Bp);
//
          CsBeam* pntr = new CsBeam(TrackTime, Bnhit, Bnhittot, Btime, Btmchiq, Bpchiq,  
                                    Hnhit, Htot, Htime, Htmchiq);
          CsHelix tmp( x, y, z, dxdz, dydz,BmChrg/Bp, cov );
          pntr->addHelix(tmp);
          //tmp.~CsHelix();   
          pntr->setChi2(TrackChi2);
          pntr->setTrigLabel(Htrig);
          pntr->setClusters(BmHod.getClusters(ihod));
          pntr->setZones( BmHod.getZones()); 
          bmtracks.push_back(pntr);
          ntrack++; 
//
//      some hists
//
	HiBmom->Fill(Bp);
        HiXprof->Fill(x); 
        HiYprof->Fill(y); 
        HiXincl->Fill(dxdz); 
        HiYincl->Fill(dydz); 
        HiTime->Fill(TrackTime); 
        HiTMchiq->Fill(TrackChi2); 
        HiSPchiq->Fill(Bpchiq); 
//
//   wipe out the combinations containing either of this pair
//
           ibms= npos[jmin][0];
           ihod= npos[jmin][1]; 
           for(i=0; i<ntot; i++) {
              if(npos[i][0]<0) continue;
              if((npos[i][0]==ibms)||(npos[i][1]==ihod)){
                  npos[i][0]=-10;
                  --nleft;
              }
           }
      }
      delete [] npos;
      delete [] dlttm; 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons:: beamprint(void)
{
      int ntrack;
      CsBeam* pntr;
      vector<CsHelix> wk;
//
        ntrack=bmtracks.size();
        cout << "CsBeamRecons:  " << ntrack << " beam track(s) found" << endl;
        if(!ntrack) return;
//
        cout << " Id    P      X      Y        Z    dx/dz   dy/dz"
             <<  "   Time TmChiq BMSTm BMStmch BMSpch HodTm Hdchq BMSht Hodht" << endl;
	cout << resetiosflags(ios::scientific);
        cout << setiosflags(ios::fixed|ios::showpoint);
//
        list<CsBeam*>::iterator ib;
        for( ib = bmtracks.begin(); ib != bmtracks.end(); ib++ ) {
        pntr = *ib;
        double TrackTime = pntr->getTrackTime();
        double BMStime = pntr->getBMStime();
        double BMStmchiq = pntr->getBMStmchi2();
        double BMSpchiq = pntr->getBMSpchi2();
        double Hodtime = pntr->getScFbTime();
        double Hodtmchiq = pntr->getScFbTimeChi2();
//        int nBMStotht = pntr->getntotBMShits();
        int nBMStrpnt = pntr->getnBMShits();
//        int nHodtotht = pntr->getntotScFbHits();
        int nHodtrpnt = pntr->getnScFbHits();
//
        int Id = pntr->getId();
        double TmChiq = pntr->getChi2();
        wk = pntr->getHelices();
          cout << setprecision(2) << setfill(' ');
          cout << setw(3)<< Id <<setw(7)<< 1./wk[0].getCop() << setw(7) << wk[0].getX() 
               << setw(7) << wk[0].getY() << setw(9) << setprecision(2) << wk[0].getZ() 
               << setw(8) << setprecision(4) << wk[0].getDXDZ()
               << setw(8) << wk[0].getDYDZ() << setprecision(2)
               << setw(6) << TrackTime << setw(7) << TmChiq << setw(6) << BMStime
               << setw(8) << BMStmchiq << setw(7)<< BMSpchiq
               << setw(6) << Hodtime << setw(6) <<Hodtmchiq<< setw(5) << nBMStrpnt
               << setw(5) << nHodtrpnt;
          cout << endl;
        }

}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons::FillCovMx(const double Dp)
//
//      Fill diagonal elements of Cov. Matrix  
//
{
int i, k;
      BmHod.getBmErrors(cov[0], cov[2], cov[5], cov[9]);
      cov[14] = Dp;                            // dp=a*p;  d(Q/p)=Q/p**2*dp = Q*a/p; 
      for(k=0, i=1; i<6; i++) {
               cov[k]=cov[k]*cov[k];
               k=k+i+1;
      }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CsBeamRecons::bookhist1()
//
// Book histograms for histoLevel1
//
{
  static bool first = true;
  if (first) {
    first = false;
    string pathname =  "/BeamRecons/BeamReconsTR";
    CsHistograms::SetCurrentPath(pathname);
    mHisto1[0]=new CsHist1D("BeamReco_03","Number of CombCandidates",50,0,50);
    mHisto1[1]=new CsHist1D("BeamReco_88","Momentum (1Cand-events)",200,0,200);
    mHisto1[3]=new CsHist1D("BeamReco_119","Multiplicity of CsBeam-objects",20,0,20);
    mHisto2[0]=new CsHist2D("BeamReco2_00","Trigger vs Category",10,0,10,10,0,10);
    mHisto2[1]=new CsHist2D("BeamReco2_01","Category vs momentum",200,50,250,5,0,5);

    Htrig1[1]=new CsHist1D("BeamTrig1_01","Number of CombCandidates (IT)",50,0,50);
    Htrig1[2]=new CsHist1D("BeamTrig1_02","Number of CombCandidates (MT)",50,0,50);
    Htrig1[3]=new CsHist1D("BeamTrig1_03","Number of CombCandidates (LT)",50,0,50);
    Htrig1[4]=new CsHist1D("BeamTrig1_04","Number of CombCandidates (OT)",50,0,50);
    Htrig1[5]=new CsHist1D("BeamTrig1_05","Number of CombCandidates (IMT)",50,0,50);

    //         ***** HISTOS for BACKTRACKING... *****
    // ...of beamTelescope track to BMS planes, possibly corrected for offset
    char hBTName[] = "hBT_BMS1";
    char hBTTitle[]  = "BackTrack: BMS1 - default offset, Reco type 0/1   ";
    //char bTTitle[] = "BackTrack: BMS1 - offset=-99.99, Reco type 0/1   ";
    for (int iBMS = 1; iBMS<=4; iBMS++) {
      int iH = 138+iBMS; // Histo index = [139,142] for BMS[1,4] 
      sprintf(hBTName,"hBT_BMS%d",iBMS);
      if (bTDefault)
	sprintf(hBTTitle,"BackTrack: BMS%d - default offset, Reco type 0/1",
		iBMS);
      else
	sprintf(hBTTitle,"BackTrack: BMS%d - offset=%6.2f, Reco type 0/1",
		iBMS,bTCorr[iBMS-1]);
      mH2[iH] = new CsHist2D(hBTName,hBTTitle,120,-30,30,2,0,2);
    }
    for (int iBMS = 5; iBMS<=6; iBMS++) {
      int iH = 188+iBMS; // Histo index = [193,194] for BMS[5,6] 
      sprintf(hBTName,"hBT_BMS%d",iBMS);
      sprintf(hBTTitle,"BackTrack: BMS%d",iBMS);
      mH1[iH] = new CsHist1D(hBTName,hBTTitle,120,-30,30);
    }
    mH2[143]=new CsHist2D("BeamReco_143","BMS_DY - Ex_BMS_DY @ BMS4",100,-5,5,2,0,2);
    mH1[179]=new CsHist1D("BeamReco_179","BMS Extrapolation: Chi2", 200,0,100);
    mH1[180]=new CsHist1D("BeamReco_180","BMS Extrapolation: Chi2 (ncut==1)", 200,0,100);
    mH1[181]=new CsHist1D("BeamReco_181","BMS Extrapolation: Chi2 (after cleaning)",200,0,100);
    mH1[182]=new CsHist1D("BeamReco_182","BMS Extrapolation: Probability (ncut==1)", 1000,0,1);
    mH1[183]=new CsHist1D("BeamReco_183","BMS Extrapolation: Probability", 1000,0,1);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
void CsBeamRecons::bookhist()
//
// Book histograms
//
{
  static bool first = true;
  char name[40];
  int nh = 1; 
  if(first){
    first = false;
    if(useTRAFFIC==0){
      string pathname =  "/BeamRecons";
      CsHistograms::SetCurrentPath(pathname);
      //
      sprintf(name,"BmRec%4.4i",nh++);
      HiNBMStrack = new CsHist1D(name, "Number of tracks/event in BMS", 30, 0, 30);
      sprintf(name,"BmRec%4.4i",nh++);
      HiNScFbtrack = new CsHist1D(name, "Number of tracks/event ScFb", 30, 0, 30);
      sprintf(name,"BmRec%4.4i",nh++);
      HiNtrack = new CsHist1D(name, "Number of the reconstructed beam tracks/event ", 30, 0, 30);
      sprintf(name,"BmRec%4.4i",nh++);
      HiBmom  =  new CsHist1D(name, "Beam momemtum ", 100, 0, 250);
      sprintf(name,"BmRec%4.4i",nh++);
      HiXprof =  new CsHist1D(name, "Beam X-profile at ScFb2 hodoscope", 96, -48*0.41, 48*0.41);
      sprintf(name,"BmRec%4.4i",nh++);
      HiYprof =  new CsHist1D(name, "Beam Y-profile at ScFb2 hodoscope", 96, -48*0.41, 48*0.41);
      sprintf(name,"BmRec%4.4i",nh++);
      HiXincl =  new CsHist1D(name, "Beam X-Z plane inclination", 100, -0.01, 0.01);
      sprintf(name,"BmRec%4.4i",nh++);
      HiYincl =  new CsHist1D(name, "Beam Y-Z plane inclination", 100, -0.01, 0.01);
      sprintf(name,"BmRec%4.4i",nh++);
      HiTime =   new CsHist1D(name, "Beam Track time", 100, -10, 10);
      sprintf(name,"BmRec%4.4i",nh++);
      HiTMchiq = new CsHist1D(name, "Beam track time chiq/point", 100, 0, 100);
      sprintf(name,"BmRec%4.4i",nh++);
      HiSPchiq = new CsHist1D(name, "Beam Track BMS space chiq", 100, 0, 100);
      sprintf(name,"BmRec%4.4i",nh++);
      HiNBNhod = new CsHist2D(name, "Number of BMS tracks vs numb. of SciFb tracks",15, 0, 14, 15, 0, 14);
      sprintf(name,"BmRec%4.4i",nh++);
      HiTBThod = new CsHist2D(name, "Time of BMS tracks vs Time of SciFb tracks",20, -10, 10, 20, -10, 10);
      sprintf(name,"BmRec%4.4i",nh++);
      HiTBTh1 = new CsHist2D(name, "Time BMS vs Time SciFb, one track events",20, -10, 10, 20, -10, 10);
      sprintf(name,"BmRec%4.4i",nh++);
      HiTmdiff = new CsHist1D(name, "Time diff. (BMS track-SciFb track)", 100, -10, 10);
      sprintf(name,"BmRec%4.4i",nh++);
      HiTmdif1 = new CsHist1D(name, "Time diff. (BMS track-SciFc track), one track events",100, -10, 10);
    }
    //
    ////
    //
    if(useTRAFFIC!=0){
      int ii;
      string pathname =  "/BeamRecons/BeamReconsTR";
      CsHistograms::SetCurrentPath(pathname);
      ii=0;
      mD1[ii]=new CsHist1D("BeamReco_00","Tracktime: BMS-SF(raw)",300,-150,150);
      ii++;//1
      mD1[ii]=new CsHist1D("BeamReco_01","Tracktime(zoom): BMS-SF(raw)",100,-10,10);
      ii++; //2
      mD1[ii]=new CsHist1D("BeamReco_02","Tracktime(Hitcut,zoom): BMS-SF",40,-10,10);
      ii++; //3
      mH1[ii]=new CsHist1D("BeamReco_03","Number of CombCandidates",50,0,50);
      ii++; //4
      mH1[ii]=new CsHist1D("BeamReco_04","Number of involved BMS-Trx",10,0,10);
      ii++; //5
      mH1[ii]=new CsHist1D("BeamReco_05","Number of involved SF-Trx",20,0,20);
      ii++; //6
      mD1[ii]=new CsHist1D("BeamReco_06","Tracktime(Hit+dt-cut,zoom): BMS-SF",100,-10,10);
      ii++;//7
      mD1[ii]=new CsHist1D("BeamReco_07","Tracktime: BMS(raw)",100,-5,5);
      ii++;//8
      mH1[ii]=new CsHist1D("BeamReco_08","BMS: ntracks(ncut=0!)",10,0,10);
      ii++;//9
      mD1[ii]=new CsHist1D("BeamReco_09","1-Cand events: thetaX (x-z-plane)",100,-0.01,0.01);
      ii++;//10
      mD1[ii]=new CsHist1D("BeamReco_10","1-Cand events: thetaY (x-y-plane)",100,-0.01,0.01);
      ii++;//11
      mD2[ii]=new CsHist2D("BeamReco_11","1-Cand events: thetaY vs thetaX",40,-0.02,0.02,40,-0.02,0.02);
      ii++;//12
      mD2[ii]=new CsHist2D("BeamReco_12","1-Cand events: thetaX vs mom",200,0,200,100,-0.005,0.005);
      ii++;//13
      mD2[ii]=new CsHist2D("BeamReco_13","1-Cand events: thetaY vs mom",200,0,200,100,-0.005,0.005);
      ii++;//14
      mD1[ii]=new CsHist1D("BeamReco_14","All candidates: thetaX",100,-0.01,0.01);
      ii++;//15
      mD1[ii]=new CsHist1D("BeamReco_15","All candidates: thetaY",100,-0.01,0.01);
      ii++;//16
      mD1[ii]=new CsHist1D("BeamReco_16","1-Cand events: X",200,-10,10);
      ii++;//17
      mH1[ii]=new CsHist1D("BeamReco_17","1-Cand events: Y",200,-10,10);
      ii++;//18
      mD1[ii]=new CsHist1D("BeamReco_18","1-track events: FI1X-ch",96,0,96);
      ii++;//19
      mD1[ii]=new CsHist1D("BeamReco_19","1-track events: FI1Y-ch",96,0,96);
      ii++;//20
      mD1[ii]=new CsHist1D("BeamReco_20","1-track events: FI2X-ch",96,0,96);
      ii++;//21
      mD1[ii]=new CsHist1D("BeamReco_21","1-track events: FI2Y-ch",96,0,96);
      ii++;//22
      mD1[ii]=new CsHist1D("BeamReco_22","1-track events: FI1X-t",600,-30,30);
      ii++;//23
      mH1[ii]=new CsHist1D("BeamReco_23","Nmb SFplanes (ncut=0,nBMStrx!=0)",5,0,5);
      ii++;//24
      mD1[ii]=new CsHist1D("BeamReco_24","SF TimeDis (ncut=0,nBMStrx!=0,nocut)",100,-50,50);
      ii++;//25
      mD1[ii]=new CsHist1D("BeamReco_25","ngoodSFTRX (ncut>1,nBMStrx=1)",100,0,100);
      ii++;//26
      mH1[ii]=new CsHist1D("BeamReco_26","SF-Efficiency",5,0,5);
      ii++;//27
      mD1[ii]=new CsHist1D("BeamReco_27","BMS-meantime(dt-cut)",300,-30,30);
      ii++;//28
      mD1[ii]=new CsHist1D("BeamReco_28","SF: meanTrTime-hitTime",300,-30,30);
      ii++;//29
      mD1[ii]=new CsHist1D("BeamReco_29","SF: meanTime(dt-cut)",300,-30,30);
      ii++;//30
      mD1[ii]=new CsHist1D("BeamReco_30","BMS & SF: TotalMeanTime",800,-80,80);
      ii++;//31
      mH1[ii]=new CsHist1D("BeamReco_31","BMS & SF: BMSold tracks",1000,0,1000);
      ii++;//32
      mH1[ii]=new CsHist1D("BeamReco_32","dt-cut-effect",10,0,10);
      
      ii++;//33
      mH1[ii]=new CsHist1D("BeamReco_33","TotalTrackTime(IT)",300,-30,30);
      ii++;//34
      mH1[ii]=new CsHist1D("BeamReco_34","TotalTrackTime(MT)",300,-30,30);
      ii++;//35
      mH1[ii]=new CsHist1D("BeamReco_35","TotalTrackTime(LT)",300,-30,30);
      ii++;//36
      mH1[ii]=new CsHist1D("BeamReco_36","BMSTrackTime(IT)",600,-60,60);
      ii++;//37
      mH1[ii]=new CsHist1D("BeamReco_37","BMSTrackTime(MT)",600,-60,60);
      ii++;//38
      mH1[ii]=new CsHist1D("BeamReco_38","BMSTrackTime(LT)",600,-60,60);
      ii++;//39
      mH1[ii]=new CsHist1D("BeamReco_39","SFTrackTime(IT)",800,-80,80);
      ii++;//40
      mH1[ii]=new CsHist1D("BeamReco_40","SFTrackTime(MT)",800,-80,80);
      ii++;//41
      mH1[ii]=new CsHist1D("BeamReco_41","SFTrackTime(LT)",800,-80,80);
      ii++;//42
      mH2[ii]=new CsHist2D("BeamReco_42","SFTrackTime vs BMSTrackTime (IT)",600,-60,60,800,-80,80);
      ii++;//43
      mH2[ii]=new CsHist2D("BeamReco_43","SFTrackTime vs BMSTrackTime (MT)",600,-60,60,800,-80,80);
      ii++;//44
      mH2[ii]=new CsHist2D("BeamReco_44","SFTrackTime vs BMSTrackTime (LT)",600,-60,60,800,-80,80);
      ii++;//45
      mH1[ii]=new CsHist1D("BeamReco_45","TriggerDistribution",5,0,5);
      ii++;//46
      mH1[ii]=new CsHist1D("BeamReco_46","SFTrackTime: SF-list",800,-80,80);
      ii++;//47
      mH1[ii]=new CsHist1D("BeamReco_47","SFTrackTime: SF-list(IT)",800,-80,80);
      ii++;//48
      mH1[ii]=new CsHist1D("BeamReco_48","SFTrackTime: SF-list(MT)",800,-80,80);
      ii++;//49
      mH1[ii]=new CsHist1D("BeamReco_49","SFTrackTime: SF-list(LT)",800,-80,80);
      ii++;//50
      mH1[ii]=new CsHist1D("BeamReco_50","BMSTrackTime: BMS-list",600,-60,60);
      ii++;//51
      mH1[ii]=new CsHist1D("BeamReco_51","BMS Number of good tracks b4 comb",50,0,50);
      ii++;//52
      mH1[ii]=new CsHist1D("BeamReco_52","ToBeBlamed totalmeancut",10,0,10);
      ii++;//53
      mH2[ii]=new CsHist2D("BeamReco_53","nSFpassed vs nBMSpassed",50,0,50,50,0,50);
      ii++;//54
      mH2[ii]=new CsHist2D("BeamReco_54","thetaX vs totalmeantime",100,-10,10,100,-0.05,0.05);
      ii++;//55
      mH2[ii]=new CsHist2D("BeamReco_55","thetaY vs totalmeantime",100,-10,10,100,-0.05,0.05);
      ii++;//56
      mH2[ii]=new CsHist2D("BeamReco_56","thetaX vs ncut",50,0,50,100,-0.05,0.05);
      ii++;//57
      mH2[ii]=new CsHist2D("BeamReco_57","thetaY vs ncut",50,0,50,100,-0.05,0.05);
      ii++;//58
      mH1[ii]=new CsHist1D("BeamReco_58","Cut-Effects on good tracks",10,0,10);
      ii++;//59
      mH2[ii]=new CsHist2D("BeamReco_59","All SF/SI from Comb: y vs x (z=-1000)",100,-50,50,100,-50,50);
      ii++;//60
      mH1[ii]=new CsHist1D("BeamReco_60","",20,-100,100);
      ii++;//61
      mH1[ii]=new CsHist1D("BeamReco_61","",20,-10,10);
      ii++;//62
      mH2[ii]=new CsHist2D("BeamReco_62","SFtrackmult vs tracktime",800,-80,80,1000,0,1000);
      ii++;//63
      mH1[ii]=new CsHist1D("BeamReco_63","SFtrackmult",3000,0,3000);
      ii++;//64
      mH1[ii]=new CsHist1D("BeamReco_64","SFTrackTime: SF-list(mult-cut)",800,-80,80);
      ii++;//65
      mH1[ii]=new CsHist1D("BeamReco_65","SF1x: HitMult",100,0,100);
      ii++;//66
      mH1[ii]=new CsHist1D("BeamReco_66","S1y: HitMult",100,0,100);
      ii++;//67
      mH1[ii]=new CsHist1D("BeamReco_67","SF2x: HitMult",100,0,100);
      ii++;//68
      mH1[ii]=new CsHist1D("BeamReco_68","SF2y: HitMult",100,0,100);
      ii++;//69
      mH2[ii]=new CsHist2D("BeamReco_69","SF1x-HitMult vs tracktime",800,-80,80,100,0,100);
      ii++;//70
      mH2[ii]=new CsHist2D("BeamReco_70","SF1y-HitMult vs tracktime",800,-80,80,100,0,100);
      ii++;//71
      mH2[ii]=new CsHist2D("BeamReco_71","SF2x-HitMult vs tracktime",800,-80,80,100,0,100);
      ii++;//72
      mH2[ii]=new CsHist2D("BeamReco_72","SF2y-HitMult vs tracktime",800,-80,80,100,0,100);
      ii++;//73
      mH2[ii]=new CsHist2D("BeamReco_73","SF1x-HitMult vs WireNumber",96,0,96,100,0,100);
      ii++;//74
      mH2[ii]=new CsHist2D("BeamReco_74","SF1y-HitMult vs WireNumber",96,0,96,100,0,100);
      ii++;//75
      mH2[ii]=new CsHist2D("BeamReco_75","SF2x-HitMult vs WireNumber",96,0,96,100,0,100);
      ii++;//76
      mH2[ii]=new CsHist2D("BeamReco_76","SF2y-HitMult vs WireNumber",96,0,96,100,0,100);
      ii++;//77
      mH2[ii]=new CsHist2D("BeamReco_77","SF1y-HitMult vs SF1x-HitMult",50,0,50,50,0,50);
      ii++;//78
      mH2[ii]=new CsHist2D("BeamReco_78","SF2y-HitMult vs SF2x-HitMult",50,0,50,50,0,50);
      ii++;//79
      mH2[ii]=new CsHist2D("BeamReco_79","SF2x-HitMult vs SF1x-HitMult",50,0,50,50,0,50);
      ii++;//80
      mH2[ii]=new CsHist2D("BeamReco_80","SF2y-HitMult vs SF1x-HitMult",50,0,50,50,0,50);
      ii++;//81
      mH2[ii]=new CsHist2D("BeamReco_81","SF1x-HitMult vs BMS_P4-HitMult",50,0,50,50,0,50);
      ii++;//82
      mH2[ii]=new CsHist2D("BeamReco_82","SF1y-HitMult vs BMS_P4-HitMult",50,0,50,50,0,50);
      ii++;//83
      mH1[ii]=new CsHist1D("BeamReco_83","Number of CombCandidates (1 Trigger)",50,0,50);
      ii++;//84
      mH1[ii]=new CsHist1D("BeamReco_84","Number of CombCandidates (IT)",50,0,50);
      ii++;//85
      mH1[ii]=new CsHist1D("BeamReco_85","Number of CombCandidates (MT)",50,0,50);
      ii++;//86
      mH1[ii]=new CsHist1D("BeamReco_86","Number of CombCandidates (LT)",50,0,50);
      ii++;//87
      mH1[ii]=new CsHist1D("BeamReco_87","nSFpassed (nBMSpassed=1)",50,0,50);
      ii++;//88
      mH1[ii]=new CsHist1D("BeamReco_88","Momentum (1Cand-events)",200,0,200);
      ii++;//89
      mH1[ii]=new CsHist1D("BeamReco_89","nBMShits: (ncut=1)",5,0,5);
      ii++;//90
      mH1[ii]=new CsHist1D("BeamReco_90","EventNr: nSFtrack>60",2000,0,20000);
      ii++;//91
      mH2[ii]=new CsHist2D("BeamReco_91","SFtrackMult vs EventNr",2000,0,20000,1800,0,900000);
      ii++;//92
      mH1[ii]=new CsHist1D("BeamReco_92","SF02X: cluster-pos",100,-20,20);
      ii++;//93
      mH1[ii]=new CsHist1D("BeamReco_93","SF02Y: cluster-pos",100,-20,20);
      ii++;//94
      mH2[ii]=new CsHist2D("BeamReco_94","SF02: clusterY vs clusterX",120,-30,30,120,-30,30);      ii++;//95
      mH2[ii]=new CsHist2D("BeamReco_95","ThetaY(SF)vsThetaY(BMS):ncut=1",80,-0.020,0.020,100,-0.005,0.005);
      ii++;//96
      mH2[ii]=new CsHist2D("BeamReco_96","ThetaX(SF)vsThetaY(BMS):ncut=1",80,-0.020,0.020,100,-0.005,0.005);
      ii++;//97
      mH2[ii]=new CsHist2D("BeamReco_97","ThetaY(SF)vsThetaY(BMS):[ncut>1]",80,-0.020,0.020,100,-0.005,0.005);
      ii++;//98
      mH1[ii]=new CsHist1D("BeamReco_98","ncut_theta",50,0,50);
      ii++;//99
      mH2[ii]=new CsHist2D("BeamReco_99","ThetaY(SF)vsThetaY(BMS):ncut_theta-cut",80,-0.020,0.020,100,-0.005,0.005);
      ii++;//100
      mH1[ii]=new CsHist1D("BeamReco_100","ncut_theta",1000,0,1000);
      ii++;//101
      mH1[ii]=new CsHist1D("BeamReco_101","Distance from ThetaCorr",200,-0.002,0.002);
      ii++;//102
      mH1[ii]=new CsHist1D("BeamReco_102","Tracktime(1Cand): BMS-SF(raw)",100,-10,10);
      ii++;//103
      mH1[ii]=new CsHist1D("BeamReco_103","SciFi: ClusterTimeDiff",100,-10,10);
      ii++;//104
      mH1[ii]=new CsHist1D("BeamReco_104","SF-Efficiency (post target)",7,0,7);
      ii++;//105
      mH1[ii]=new CsHist1D("BeamReco_105","Number of CombCandidates (goodTrigger)",50,0,50);
      ii++;//106
      mH1[ii]=new CsHist1D("BeamReco_106","SF: number of tracks with time per event",100,0,100);
      ii++;//107
      mH2[ii]=new CsHist2D("BeamReco_107","nCluster vs invol Det(SF 0-5,Si 6-11 ",14,0,14,20,0,20);
      ii++;//108
      mH2[ii]=new CsHist2D("BeamReco_108","ntrack vs # invol SF Det",8,0,8,50,0,50);
      ii++;//109
      mH1[ii]=new CsHist1D("BeamReco_109","ncut=0: Nr of Comb after dt-cut",20,0,20);
      ii++;//110
      mH2[ii]=new CsHist2D("BeamReco_110","ncut=0 & nComb=0: SFtracktime vs BMStracktime",200,-20,20,50,-5,5);
      ii++;//111
      mH2[ii]=new CsHist2D("BeamReco_111","ncut>1 & nBMStracks>1: n4fired vs n3fired",20,0,20,20,0,20);
      ii++;//112
      mH2[ii]=new CsHist2D("BeamReco_112","ncut_theta vs ncut",30,0,30,30,0,30);
      ii++;//113
      mH2[ii]=new CsHist2D("BeamReco_113","all SF tracks in Comb: x vs y @-350mm",100,-50,50,100,-50,50);
      ii++;//114
      mH1[ii]=new CsHist1D("BeamReco_114","all SF tracks in Comb: dx @-350mm",200,-0.01,0.01);
      ii++;//115
      mH1[ii]=new CsHist1D("BeamReco_115","all SF tracks in Comb: dy @-350mm",200,-0.01,0.01);
      ii++;//116
      mH2[ii]=new CsHist2D("BeamReco_116","goodEvents: all SF tracks in Comb: x vs y @-350mm",100,-50,50,100,-50,50);
      ii++;//117
      mH1[ii]=new CsHist1D("BeamReco_117","goodEvents: all SF tracks in Comb: dx @-350mm",200,-0.01,0.01);
      ii++;//118
      mH1[ii]=new CsHist1D("BeamReco_118","goodEvents: all SF tracks in Comb: dy @-350mm",200,-0.01,0.01);
      ii++;//119
      mH1[ii]=new CsHist1D("BeamReco_119","Multiplicity of CsBeam-objects",20,0,20);
      ii++;//120
      mH2[ii]=new CsHist2D("BeamReco_120","Y(SF)vsY(BMS): ncut=1",80,-200,200,200,-100,100);
      ii++;//121
      mH2[ii]=new CsHist2D("BeamReco_121","Y(SF)vsY(BMS): ncut>1",80,-200,200,200,-100,100);
      ii++;//122
      mH2[ii]=new CsHist2D("BeamReco_122","SF/SI from Comb, ncut==1: y vs x (z=-1000)",100,-50,50,100,-50,50);
      ii++;//123
      mH1[ii]=new CsHist1D("BeamReco_123","Rescue: Tracktime(zoom): BMS-SF(raw)",100,-10,10);
      ii++;//124
      mH1[ii]=new CsHist1D("BeamReco_124","Rescue: BMS & SF: TotalMeanTime",800,-80,80);
      ii++;//125
      mH1[ii]=new CsHist1D("BeamReco_125","Rescue: ncut",100,0,100);
      ii++;//126
      mH1[ii]=new CsHist1D("BeamReco_126","Rescue: ncut==0,ntrack_resc",50,0,50);
      ii++;//127
      mH2[ii]=new CsHist2D("BeamReco_127","Rescue: nSFpassed vs nBMSpassed",50,0,50,50,0,50);
      ii++;//128
      mH1[ii]=new CsHist1D("BeamReco_128","Rescue: ncut(goodTrigger) ",100,0,100);
      ii++;//129
      mH1[ii]=new CsHist1D("BeamReco_129","BMS & SF: TotalMeanTime (before cuts)",800,-80,80);
      ii++;//130
      mH2[ii]=new CsHist2D("BeamReco_130","nCluster vs nSFcluster ",6,0,6,14,0,14);
      
      mH2[131]=new CsHist2D("BeamReco_131","VI1 timing - refTime",300,-1200,0,4,0,4);
      mH2[132]=new CsHist2D("BeamReco_132","VI2 timing - refTime",300,-1200,0,4,0,4);
      mH2[133]=new CsHist2D("BeamReco_133","VO1(center) timing - refTime",300,-1200,0,4,0,4);
      mH2[134]=new CsHist2D("BeamReco_134","VO1(meantimed part)timing - refTime",300,-1200,0,18,0,18);
      double histcenter=-920;
      mH2[135]=new CsHist2D("BeamReco_135","VI1 timing - refTime(cal)",300,histcenter-30,histcenter+30,4,0,4);
      mH2[136]=new CsHist2D("BeamReco_136","VI2 timing - refTime(cal)",300,histcenter-30,histcenter+30,4,0,4);
      mH2[137]=new CsHist2D("BeamReco_137","VO1(center) timing - refTime(cal)",300,histcenter-30,histcenter+30,4,0,4);
      mH2[138]=new CsHist2D("BeamReco_138","VO1(meantimed part)timing - refTime(cal)",300,histcenter-30,histcenter+30,18,0,18);

      //         ***** HISTOS for BACKTRACKING... *****
      // ...of beamTelescope track to BMS planes, possibly corrected for offset
      char hBTName[] = "hBT_BMS1";
      char hBTTitle[]  = "BackTrack: BMS1 - default offset, Reco type 0/1   ";
      //char bTTitle[] = "BackTrack: BMS1 - offset=-99.99, Reco type 0/1   ";
      for (int iBMS = 1; iBMS<=4; iBMS++) {
	int iH = 138+iBMS; // Histo index = [139,142] for BMS[1,4] 
	sprintf(hBTName,"hBT_BMS%d",iBMS);
	if (bTDefault)
	  sprintf(hBTTitle,"BackTrack: BMS%d - default offset, Reco type 0/1",
		  iBMS);
	else
	  sprintf(hBTTitle,"BackTrack: BMS%d - offset=%6.2f, Reco type 0/1",
		  iBMS,bTCorr[iBMS-1]);
	mH2[iH] = new CsHist2D(hBTName,hBTTitle,120,-30,30,2,0,2);
      }
      for (int iBMS = 5; iBMS<=6; iBMS++) {
	int iH = 188+iBMS; // Histo index = [193,194] for BMS[5,6] 
	sprintf(hBTName,"hBT_BMS%d",iBMS);
	sprintf(hBTTitle,"BackTrack: BMS%d",iBMS);
	mH1[iH] = new CsHist1D(hBTName,hBTTitle,120,-30,30);
      }
      mH2[143]=new CsHist2D("BeamReco_143","BMS_DY - Ex_BMS_DY @ BMS4",100,-5,5,2,0,2);
      mH2[144]=new CsHist2D("BeamReco_144","BMS_Y - Ex_BMS_Y @ BMS1",115,-115,115,115,-115,115);
      mH2[145]=new CsHist2D("BeamReco_145","BMS_Y - Ex_BMS_Y @ BMS2",115,-115,115,115,-115,115);
      mH2[146]=new CsHist2D("BeamReco_146","BMS_Y - Ex_BMS_Y @ BMS3",115,-115,115,115,-115,115);
      mH2[147]=new CsHist2D("BeamReco_147","BMS_Y - Ex_BMS_Y @ BMS4",115,-115,115,115,-115,115);
      mH2[148]=new CsHist2D("BeamReco_148","BMS Extrapolation: nHits(BMS) vs Chi2",50,0,100,2,0,2);      
      mH1[149]=new CsHist1D("BeamReco_149","BMS Extrapolation: Best Chi2 [ncut>1,nSFpassed==1]",40,0,40);
      mH1[150]=new CsHist1D("BeamReco_150","BMS Extrapolation: Ratio Best/2nd Chi2 [ncut>1,nSFpassed==1]",100,0,1);
      mH1[151]=new CsHist1D("BeamReco_151","BMS Extrapolation: Best Chi2 [ncut>1,nBMSpassed==1]",40,0,40);
      mH1[152]=new CsHist1D("BeamReco_152","BMS Extrapolation: Ratio Best/2nd Chi2 [ncut>1,nBMSpassed==1]",100,0,1);
      mH2[153]=new CsHist2D("BeamReco_153","Errors @ BMS1 [0: Y,1: DY,2: Dpp,3: All]",50,0,10,4,0,4);
      mH2[154]=new CsHist2D("BeamReco_154","Errors @ BMS2 [0: Y,1: DY,2: Dpp,3: All]",50,0,10,4,0,4);
      mH2[155]=new CsHist2D("BeamReco_155","Errors @ BMS3 [0: Y,1: DY,2: Dpp,3: All]",50,0,10,4,0,4);
      mH2[156]=new CsHist2D("BeamReco_156","Errors @ BMS4 [0: Y,1: DY,2: Dpp,3: All]",50,0,10,4,0,4);
      mH2[157]=new CsHist2D("BeamReco_157","Probabilities",200,0,1,15,0,15);
      mH2[158]=new CsHist2D("BeamReco_158","Best/SecBest Probabilities vs Best Prob",50,0,1,50,0,1);
      mH1[159]=new CsHist1D("BeamReco_159","BMS Extrapolation: Prob Multiplicties [ncut>1,nSFpassed==1]",10,0,10);
      mH1[160]=new CsHist1D("BeamReco_160","BMS Extrapolation: Prob Multiplicties [ncut>1,nBMSpassed==1]",10,0,10);
      mH2[161]=new CsHist2D("BeamReco_161","2ndBest Probabilities vs Best Prob",50,-3,0,50,-3,0);
      mH1[162]=new CsHist1D("BeamReco_162","BMS Extrapolation: Cand Mult(best Chi2) [ncut>1,nSFpassed==1]",10,0,10);
      mH1[163]=new CsHist1D("BeamReco_163","BMS Extrapolation: Cand Mult(best Chi2) [ncut>1,nBMSpassed==1]",10,0,10);
      mH1[164]=new CsHist1D("BeamReco_164","BMS Extrapolation: Cand Mult(best Chi2) [ncut>1,both>1]",10,0,10);
      mH1[165]=new CsHist1D("BeamReco_165","BMS Extrapolation: Ratio Best/2nd Chi2 [nBMS>1,nSFpassed>1]",100,0,1);
      mH2[166]=new CsHist2D("BeamReco_166","Best/SecBest Chi2 vs Best Chi [ncut>1,nSFpassed==1]",50,0,10,50,0,1);
      mH2[167]=new CsHist2D("BeamReco_167","Best/SecBest Chi2 vs Best Chi [ncut>1,nBMSpassed==1]",50,0,10,50,0,1);
      mH2[168]=new CsHist2D("BeamReco_168","Best/SecBest Chi2 vs Best Chi [both >1]",50,0,10,50,0,1);
      mH2[169]=new CsHist2D("BeamReco_169","Chi2 ratio > 0.7: mom1 vs mom2 [ncut>1,nBMSpassed==1]",30,140,200,30,-10,10);
      mH1[170]=new CsHist1D("BeamReco_170","Rescued tracks [Rolands method]: momentum",120,100,220);
      mH1[171]=new CsHist1D("BeamReco_171","GoodTrigger CsBeam-Multiplicity",10,0,10);
      mH1[172]=new CsHist1D("BeamReco_172","SF: number of tracks with time per event (GoodTrigger)",100,0,100);
      mH1[173]=new CsHist1D("BeamReco_173","BMS & SF: BMS05 tracks",1000,0,1000);
      mH1[174]=new CsHist1D("BeamReco_174","BMS & SF: BMS06 tracks",1000,0,1000);
      mH1[175]=new CsHist1D("BeamReco_175","Events: 0-BMSold track, 1-BMS05 track, 2-BMS06 track, 3-GoodBMS05, 4-GoodBMS06",5,0,5);
      mH1[176]=new CsHist1D("BeamReco_176","BMS & SF: BMSold tracks (good trigger)",1000,0,1000);
      mH1[177]=new CsHist1D("BeamReco_177","BMS & SF: BMS05 tracks (good trigger)",1000,0,1000);
      mH1[178]=new CsHist1D("BeamReco_178","BMS & SF: BMS06 tracks (good trigger)",1000,0,1000);
      mH1[179]=new CsHist1D("BeamReco_179","BMS Extrapolation: Chi2", 200,0,100);
      mH1[180]=new CsHist1D("BeamReco_180","BMS Extrapolation: Chi2 (ncut==1)", 200,0,100);
      mH1[181]=new CsHist1D("BeamReco_181","BMS Extrapolation: Chi2 (after cleaning)",200,0,100);
      mH1[182]=new CsHist1D("BeamReco_182","BMS Extrapolation: Probability (ncut==1)", 1000,0,1);
      mH1[183]=new CsHist1D("BeamReco_183","BMS Extrapolation: Probability", 1000,0,1);
      mH1[184]=new CsHist1D("BeamReco_184","BMS Extrapolation: Probability (after cleaning)",1000,0,1);
      mH1[185]=new CsHist1D("BeamReco_185","BMS Extrapolation: Candidates (before prob. cut)",50,0,50);
      mH2[186]=new CsHist2D("BeamReco_186","Cand Mult.(after Chi2 cut) vs Cand Mult.(before Chi2 cut)",6,0,6,6,0,6);
      mH2[187]=new CsHist2D("BeamReco_187","Best Prob vs Second Prob [2 unambigous candidates killed]",1000,0,0.005,1000,0,0.005);
      mH1[188]=new CsHist1D("BeamReco_188","BMS momentum: P(chi2) > 0.5% (1 Cand)", 100, 120, 220);
      mH1[189]=new CsHist1D("BeamReco_189","BMS momentum: P(chi2) < 0.5% (1 Cand)", 100, 120, 220);
      mH1[190]=new CsHist1D("BeamReco_190","BMS momentum: P(chi2) > 0.5% (>1 Cand)", 100, 120, 220);
      mH1[191]=new CsHist1D("BeamReco_191","BMS momentum: P(chi2) < 0.5% (>1 Cand)", 100, 120, 220);
      mH1[192]=new CsHist1D("BeamReco_192","Events: 0-BMSold track, 1-BMS05 track, 2-BMS06 track, 3-GoodBMS05, 4-GoodBMS06 (good trigger)",5,0,5);
      HiBeamHitTimes = new CsHist1D("BeamHitTimes","Hit times in BMS planes for BeamCand", 500, -100,100);
    }
    CsHistograms::SetCurrentPath("/BeamRecons/BeamReconsTR/TriggerStudy/");
    int ii=0;//0
    Htrig1[ii]=new CsHist1D("BeamTrig1_00","Number of CombCandidates (1 Trigger)",50,0,50);
    ii++;//1
    Htrig1[ii]=new CsHist1D("BeamTrig1_01","Number of CombCandidates (IT)",50,0,50);
    ii++;//2
    Htrig1[ii]=new CsHist1D("BeamTrig1_02","Number of CombCandidates (MT)",50,0,50);
    ii++;//3
    Htrig1[ii]=new CsHist1D("BeamTrig1_03","Number of CombCandidates (LT)",50,0,50);
    ii++;//4
    Htrig1[ii]=new CsHist1D("BeamTrig1_04","Number of CombCandidates (OT)",50,0,50);
    ii++;//5
    Htrig1[ii]=new CsHist1D("BeamTrig1_05","Number of CombCandidates (IMT)",50,0,50);
    ii++;//6
    Htrig1[ii]=new CsHist1D("BeamTrig1_06","SF-Efficiency(cleaned Trigger)",5,0,5);
    ii++;//7
    Htrig1[ii]=new CsHist1D("BeamTrig1_07", "RecPerformance vs trigger",100,0,100);
    ii++;//8
    Htrig1[ii]=new CsHist1D("BeamTrig1_08","BMS ntracks (one trigger)",30,0,30);
    ii++;//9
    Htrig1[ii]=new CsHist1D("BeamTrig1_09","BMS ntracks (IT)",30,0,30);
    ii++;//10
    Htrig1[ii]=new CsHist1D("BeamTrig1_10","BMS ntracks (MT)",30,0,30);
    ii++;//11
    Htrig1[ii]=new CsHist1D("BeamTrig1_11","BMS ntracks (LT)",30,0,30);
    ii++;//12
    Htrig1[ii]=new CsHist1D("BeamTrig1_12","BMS ntracks (OT)",30,0,30);
    ii++;//13
    Htrig1[ii]=new CsHist1D("BeamTrig1_13","BMS ntracks (IMT)",30,0,30);
    ii++;//14
    Htrig1[ii]=new CsHist1D("BeamTrig1_14","BMS ntracks (CT)",30,0,30);
    ii++;//15
    Htrig1[ii]=new CsHist1D("BeamTrig1_15","Number of CombCandidates (CT)",50,0,50);
    ii++;//16
    Htrig1[ii]=new CsHist1D("BeamTrig1_16","SF-Efficiency(cleaned Trigger,rel to RefTime, 4 planes)",10,0,10);
    ii++;//17
    Htrig1[ii]=new CsHist1D("BeamTrig1_17","SF-Efficiency(cleaned Trigger,rel to RefTime, 6 planes)",10,0,10);
    ii++;//18
    Htrig1[ii]=new CsHist1D("BeamTrig1_18","PreT SF: timeDiff between clusters",100,-10,10);
    //
    ii=0;
    Htrig2[ii]=new CsHist2D("BeamTrig2_00","nGoodSFtracks Vs nBMStracks (all Trigger)",30,0,30,30,0,30);
    ii++;//1
    Htrig2[ii]=new CsHist2D("BeamTrig2_01","nGoodSFtracks Vs nBMStracks (one Trigger)",30,0,30,30,0,30);
    ii++;//2
    Htrig2[ii]=new CsHist2D("BeamTrig2_02","nGoodSFtracks Vs nBMStracks (IT)",30,0,30,30,0,30);
    ii++;//3
    Htrig2[ii]=new CsHist2D("BeamTrig2_03","nGoodSFtracks Vs nBMStracks (MT)",30,0,30,30,0,30);
    ii++;//4
    Htrig2[ii]=new CsHist2D("BeamTrig2_04","nGoodSFtracks Vs nBMStracks (LT)",30,0,30,30,0,30);
    ii++;//5
    Htrig2[ii]=new CsHist2D("BeamTrig2_05","nGoodSFtracks Vs nBMStracks (OT)",30,0,30,30,0,30);
    ii++;//6
    Htrig2[ii]=new CsHist2D("BeamTrig2_06","nGoodSFtracks Vs nBMStracks (IMT)",30,0,30,30,0,30);
    ii++;//7
    Htrig2[ii]=new CsHist2D("BeamTrig2_07","nGoodSFtracks Vs nBMStracks (CT)",30,0,30,30,0,30);
    ii++;//8
    Htrig2[ii]=new CsHist2D("BeamTrig2_08","nGoodSFtracks Vs nBMStracks (all events)",30,0,30,30,0,30);
    ii++;//9
    Htrig2[ii]=new CsHist2D("BeamTrig2_09","nGoodSFtracks Vs nBMStracks (clean Trigger)",30,0,30,30,0,30);
    ii++;//10
    Htrig2[ii]=new CsHist2D("BeamTrig2_10","SFtrackTime vs trigger (SFlist)",200,-20,20,10,0,10);
    ii++;//11
    Htrig2[ii]=new CsHist2D("BeamTrig2_11","BMStrackTime vs trigger (BMSlist)",200,-20,20,10,0,10);
    ii++;//12
    Htrig2[ii]=new CsHist2D("BeamTrig2_12","SF station vs HitTime-RefTime [goodTrigger,!(ncut==0&&!GoodSFTrack]",150,-15, 15,7,0,7);
    ii++;//13
    Htrig2[ii]=new CsHist2D("BeamTrig2_13","Filter-Effect: (0/1:0Cand,2/3:1Cand,4/5:>1Cand)vs Trigger",6,0,6,20,0,20);
    Htrig2[14]=new CsHist2D("BeamTrig2_14","SF station vs HitTime-RefTime [ncut==0,!GoodSFtrack]",150,-15, 15,7,0,7);
    Htrig2[15]=new CsHist2D("BeamTrig2_15","SF station vs HitTime-RefTime [!goodTrigger]",150,-15, 15,7,0,7);
    CsHistograms::SetCurrentPath("/");
    if(doRandom){
      CsHistograms::SetCurrentPath("/BeamRecons/BeamReconsTR/RandomTrigger/");
      Ratrig1[0]=new CsHist1D("RaTrig1_00","number of tracks",30,0,30);
      Ratrig1[1]=new CsHist1D("RaTrig1_01","tracks in spec: momentum",200,0,200);
      Ratrig1[2]=new CsHist1D("RaTrig1_02","Participating SF-stations (19:all)",21,0,21);
      Ratrig1[3]=new CsHist1D("RaTrig1_03","Number of Participating SF-stations",20,0,20);
      Ratrig1[4]=new CsHist1D("RaTrig1_04","nSF>=13: number of tracks with momentum",30,0,30);
      Ratrig1[5]=new CsHist1D("RaTrig1_05","GoodRandomTrigSFtime",100,-10,10);
      Ratrig1[6]=new CsHist1D("RaTrig1_06","GoodRandomTrigSFtime-totalMeanTime",200,-20,20);
      Ratrig1[7]=new CsHist1D("RaTrig1_07","nComb=0: BMStrtime-GoodRandomTrigSFtime",400,-40,40);
      Ratrig1[8]=new CsHist1D("RaTrig1_08","nComb=0: SFtrtime-GoodRandomTrigSFtime",200,-20,20);
      Ratrig1[9]=new CsHist1D("RaTrig1_09","ncut=1: BMStrtime-GoodRandomTrigSFtime",200,-20,20);
      Ratrig1[10]=new CsHist1D("RaTrig1_10","ncut=1: SFtrtime-GoodRandomTrigSFtime",200,-20,20);
      Ratrig1[11]=new CsHist1D("RaTrig1_11","nComb=0: 0:BMSok,1:SFok,2:both ok,3:both !ok",10,0,10);
      Ratrig1[12]=new CsHist1D("RaTrig1_12","nComb=0,!SFok: 0->19:SFhastime,20->40:dito[BMSok]",40,0,40);
      Ratrig1[13]=new CsHist1D("RaTrig1_13","good tracks in spec: momentum",200,0,200);
      Ratrig1[14]=new CsHist1D("RaTrig1_14","number of good tracks per event",20,0,20);
      Ratrig1[15]=new CsHist1D("RaTrig1_15","goodEvent,ncut=1,nBMS=4: momentum in spec",200,0,200);
      Ratrig1[16]=new CsHist1D("RaTrig1_16","goodEvent: dx/dz",200,-0.01,0.01);
      Ratrig1[17]=new CsHist1D("RaTrig1_17","goodEvent: dy/dz",200,-0.01,0.01);
      Ratrig1[18]=new CsHist1D("RaTrig1_18","11 SF,noSF7: clusterTime-RefTime in SF07x",200,-20,20);
      Ratrig1[19]=new CsHist1D("RaTrig1_19","11 SF,noSF7: clusterTime-RefTime in SF07y",200,-20,20);
      Ratrig1[20]=new CsHist1D("RaTrig1_20","goodEvent: Chi2",100,0,50);
      Ratrig1[21]=new CsHist1D("RaTrig1_21","nSF==13: Chi2",100,0,50);
      Ratrig1[22]=new CsHist1D("RaTrig1_22","ncut=1: BMSmom-SFmom",800,-100,100);
      Ratrig1[23]=new CsHist1D("RaTrig1_23","Rescue: GoodRandomTrigSFtime-totalMeanTime",200,-20,20);
      Ratrig1[24]=new CsHist1D("RaTrig1_24","Rescue: nComb=0: BMStrtime-GoodRandomTrigSFtime",400,-40,40);
      Ratrig1[25]=new CsHist1D("RaTrig1_25","Rescue: nComb=0: 0:BMSok,1:SFok,2:both ok,3:both !ok",10,0,10);
      Ratrig1[26]=new CsHist1D("RaTrig1_26","nComb=0,!SFok: nfired SF", 10,0, 10);
      Ratrig1[27]=new CsHist1D("RaTrig1_27","nComb=0,!SFok: nfiredSF=0: SpillNumber", 100,0, 100);
      Ratrig1[28]=new CsHist1D("RaTrig1_28","nComb=0: z @ SF3 ", 200,100, 300);
      Ratrig1[29]=new CsHist1D("RaTrig1_29","nComb=0:nfired SF=0: theta @ SF3 ", 100,0, 10);
      Ratrig1[30]=new CsHist1D("RaTrig1_30","nComb=0:nfired SF>0: theta @ SF3 ", 100,0, 10);
      Ratrig1[31]=new CsHist1D("RaTrig1_31","nComb=0,nfired SF=0: GoodRandomTrigSFtime",100,-5,5);
      Ratrig1[32]=new CsHist1D("RaTrig1_32","Distance from ThetaCorr[ncut==1]",200,-0.002,0.002);
      Ratrig1[33]=new CsHist1D("RaTrig1_33","Distance from ThetaCorr[nSFpassed==1 && nBMSpassed>1]",200,-0.002,0.002);
      Ratrig1[34]=new CsHist1D("RaTrig1_34","ncut>1: 4BMS-planes: Space Chi2",100,0,50);
      Ratrig1[35]=new CsHist1D("RaTrig1_35","ncut>1: 4BMS-planes: Space Chi2/ratio",100,0,1);
      Ratrig1[36]=new CsHist1D("RaTrig1_36","ncut>1: 4BMS-planes: Multiplicity",10,0,10);

      Ratrig2[0]=new CsHist2D("RaTrig2_00","nComb=0,!SFok: SF-telescope: ntrack with time vs w/o time",20,0,20,20,0,20);
      Ratrig2[1]=new CsHist2D("RaTrig2_01","nComb=0,!SFok,BMSok: SF-telescope: ntrack with time vs w/o time",20,0,20,20,0,20);
      Ratrig2[2]=new CsHist2D("RaTrig2_02","good Event: nBMSfired vs spec-mom",200,0,200,5,0,5);
      Ratrig2[3]=new CsHist2D("RaTrig2_03","good Event: BeamSpot at -350 (yvs x)",300,-150,150,300,-150,150);
      Ratrig2[4]=new CsHist2D("RaTrig2_04","Number of part. SpecSF vs dx(from SpecTrack)",200,-0.01,0.01,20,0,20);
      Ratrig2[5]=new CsHist2D("RaTrig2_05","Number of part. SpecSF vs dy(from SpecTrack)",200,-0.01,0.01,20,0,20);
      Ratrig2[6]=new CsHist2D("RaTrig2_06","goodEvents: nSFpassed vs nBMSpassed",50,0,50,50,0,50);
      Ratrig2[7]=new CsHist2D("RaTrig2_07","Number of part. SpecSF vs part SFstationsnSF (all tacks)",20,0,20,20,0,20);
      Ratrig2[8]=new CsHist2D("RaTrig2_08","ncut==1:SpecMom vs BMSmom",200,0,200,200,0,200);
      Ratrig2[9]=new CsHist2D("RaTrig2_09","Number of part. SpecSF vs part SFstationsnSF (goodTracks)",20,0,20,20,0,20);
      Ratrig2[10]=new CsHist2D("RaTrig2_10","nComb=0,!SFok: nfired SF vs which SF fired",5,0,5,5,0,5);
      Ratrig2[11]=new CsHist2D("RaTrig2_11","nComb=0,no SF plane fired: nAcept vs SFplane !in acept",5,0,5,5,0,5);
      Ratrig2[12]=new CsHist2D("RaTrig2_12","nComb=0,no SF plane fired: nfiredVT vs which VT has fired",5,0,5,5,0,5);
      Ratrig2[13]=new CsHist2D("RaTrig2_13","nComb=0,>0 SF plane fired: nfiredVT vs which VT has fired",5,0,5,5,0,5);
      Ratrig2[14]=new CsHist2D("RaTrig2_14","nComb=0,no SF plane fired: y vs x @ SF3",80,-20,20,80,-20,20);
      Ratrig2[15]=new CsHist2D("RaTrig2_15","nComb=0,>0 SF plane fired: y vs x @ SF3",80,-20,20,80,-20,20);
      Ratrig2[16]=new CsHist2D("RaTrig2_16","ThetaY(SF)vsThetaY(BMS):[ncut==1]",80,-0.020,0.020,100,-0.005,0.005);
      Ratrig2[17]=new CsHist2D("RaTrig2_17","ThetaY(SF)vsThetaY(BMS):[nSFpassed=1 && nBMSpassed>1]",80,-0.020,0.020,100,-0.005,0.005);

      HiExBmsXY.resize(NBMSHODS,NULL);
      char title[100];
      for(int i=0; i<NBMSHODS; i++) {
	  sprintf(name,"ExBmsXY%4.4i",i+1);
	  sprintf(title,"Extrapolated Y(X) coordinates from SF: BMS%i",i+1);
	  HiExBmsXY[i] = new CsHist2D(name,title,400,-200,200,400,-200,200);
      }
      HiExBmsCorelX.resize(7,NULL);
      HiExBmsCorelY.resize(7,NULL);
      
      HiExBmsCorelX[0] = new CsHist2D("ExBmsCorelX01","Extrapolated BMS: X coordinate BMS01 vs BMS02",400,-200,200,400,-200,200);
      HiExBmsCorelX[1] = new CsHist2D("ExBmsCorelX02","Extrapolated BMS: X coordinate BMS02 vs BMS03",400,-200,200,400,-200,200);
      HiExBmsCorelX[2] = new CsHist2D("ExBmsCorelX03","Extrapolated BMS: X coordinate BMS03 vs BMS04",400,-200,200,400,-200,200);
      HiExBmsCorelX[3] = new CsHist2D("ExBmsCorelX04","Extrapolated BMS: X coordinate BMS01 vs BMS05",400,-200,200,400,-200,200);
      HiExBmsCorelX[4] = new CsHist2D("ExBmsCorelX05","Extrapolated BMS: X coordinate BMS02 vs BMS05",400,-200,200,400,-200,200);
      HiExBmsCorelX[5] = new CsHist2D("ExBmsCorelX06","Extrapolated BMS: X coordinate BMS03 vs BMS06",400,-200,200,400,-200,200);
      HiExBmsCorelX[6] = new CsHist2D("ExBmsCorelX07","Extrapolated BMS: X coordinate BMS04 vs BMS06",400,-200,200,400,-200,200);

      HiExBmsCorelY[0] = new CsHist2D("ExBmsCorelY01","Extrapolated BMS: Y coordinate BMS01 vs BMS02",400,-200,200,400,-200,200);
      HiExBmsCorelY[1] = new CsHist2D("ExBmsCorelY02","Extrapolated BMS: Y coordinate BMS02 vs BMS03",400,-200,200,400,-200,200);
      HiExBmsCorelY[2] = new CsHist2D("ExBmsCorelY03","Extrapolated BMS: Y coordinate BMS03 vs BMS04",400,-200,200,400,-200,200);
      HiExBmsCorelY[3] = new CsHist2D("ExBmsCorelY04","Extrapolated BMS: Y coordinate BMS01 vs BMS05",400,-200,200,400,-200,200);
      HiExBmsCorelY[4] = new CsHist2D("ExBmsCorelY05","Extrapolated BMS: Y coordinate BMS02 vs BMS05",400,-200,200,400,-200,200);
      HiExBmsCorelY[5] = new CsHist2D("ExBmsCorelY06","Extrapolated BMS: Y coordinate BMS03 vs BMS06",400,-200,200,400,-200,200);
      HiExBmsCorelY[6] = new CsHist2D("ExBmsCorelY07","Extrapolated BMS: Y coordinate BMS04 vs BMS06",400,-200,200,400,-200,200);
      
      
      HiExBmsFired = new CsHist1D("ExBmsFired","# of BMS planes \"fired\" by ex. tracks",7,0,7);
      HiExBmsOldFired = new CsHist1D("ExBmsOldFired","# of BMS planes \"fired\" by ex. tracks",5,0,5);
      HiExBmsStuff = new CsHist1D("ExBmsStuff","Extrapolated BMS: 0 = events with good tracks, 1 = events with good tracks (old BMS), 2 = BMS5 used, 3 = BMS6 used",4,0,4);
      HiExBmsEff = new CsHist1D("ExBmsEff","Extrapolated BMS: 0-5 = number of tracks that fired BM01-06",6,0,6);
    }  
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
CsBeamRecons::~CsBeamRecons()
{
        delete BMS;
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBeamRecons::reconsMcExect()
{  
static bool first=true; 
//
    const int nhod=4;
    static CsDetector* Id[nhod];
//
        if(first){
          int i, unit;
          string name;
          const string hodnames[nhod]= 
                          { "ScFbHod1X", "ScFbHod1Y", "ScFbHod2X", "ScFbHod2Y"};
          CsOpt* opt = CsOpt::Instance();
          list <CsDetector*>  det = CsGeom::Instance()->getDetectors();
          list<CsDetector*>::iterator idet;
          list<CsZone*> detzone; 
          for(i=0; i<nhod; i++){
             Id[i]=0; 
             if( opt->getOpt( hodnames[i], "Unit", unit) &&
                 opt->getOpt( hodnames[i], "Name", name)) {
                   for( idet = det.begin(); idet != det.end(); idet++ ) {
                   CsDetector* detpnt = *idet;
  	              if( detpnt ) {
                            if((detpnt->getName()==name)&&(detpnt->getUnit()==unit)) {
                               Id[i]=detpnt;
                            }
                      }
                   }  
              }
           }
           first=false;
         }
   if(!Id[0]||!Id[1]||!Id[2]||!Id[3]) return;
//
   CsMCHit* hit;
   vector<CsMCHit*> h0, h1, h2, h3;
   vector<CsCluster*> clh0, clh1, clh2, clh3;
   list<CsCluster*> clust = CsEvent::Instance()->getClusters();
   if(clust.empty()) return;
   list<CsCluster*>::iterator nextclust;
   h0.clear(); h1.clear(); h2.clear(); h3.clear();
   clh0.clear(); clh1.clear(); clh2.clear(); clh3.clear();  
   for( nextclust = clust.begin(); nextclust != clust.end(); nextclust++ ) {
      list<CsDetector*> det = (*nextclust)->getDetsList();
      if(det.empty()) continue;
      list<CsDetector*>::iterator idet=det.begin();         // It should be only 1 detector type here
          for(int i=0; i<nhod; i++){
              if((*idet)!=Id[i]) continue; 
              list<CsDigit*>::iterator dig;
              list<CsDigit*> digits = (*nextclust)->getDigitsList();
              if( digits.empty() ) continue;
              for( dig=digits.begin(); dig!=digits.end(); dig++ ) {
                  CsMCDigit* mcdig = dynamic_cast<CsMCDigit*>(*dig);
                  if( mcdig == 0 ) continue;               
                  hit=(mcdig->getHits()).front();         // It should be only 1 hit!!
                  if(hit->getOrigin()) continue;          // Take only hits from original track
                  if(i==0)      { h0.push_back(hit); clh0.push_back(*nextclust); }
                  else if(i==1) { h1.push_back(hit); clh1.push_back(*nextclust); }
                  else if(i==2) { h2.push_back(hit); clh2.push_back(*nextclust); }
                  else if(i==3) { h3.push_back(hit); clh3.push_back(*nextclust); }
              }
          }  
   }
//
   if(h0.empty()||h1.empty()||h2.empty()||h3.empty()) return;  
//
   unsigned int i, j, k, l;
   CsMCTrack* track;
   double x,y,z,x0,z0,x1,z1,x3,z3,dxdz,dydz, Pmom;
   double timebm;
   Hep3Vector p; 
   for( i = 0; i < h0.size(); i++ ) {
       track = (h0[i])->getMCTrack(); 
       x0=(clh0[i])->getU();     
       z0=(clh0[i])->getW();
       for( j = 0; j < h1.size(); j++ ) {
          if(track!=(h1[j])->getMCTrack()) continue;
          x1=(clh1[j])->getU();
          z1=(clh1[j])->getW();
          for( k = 0; k < h2.size(); k++ ) {
             if(track!=(h2[k])->getMCTrack()) continue;
             x =(clh2[k])->getU();
             z =(clh2[k])->getW();
             (clh2[k])->getTime(timebm);
             for( l = 0; l < h3.size(); l++ ) {
                 if(track!=(h3[l])->getMCTrack()) continue;
                 x3=(clh3[l])->getU();
                 z3=(clh3[l])->getW();
                 dxdz=(x-x0)/(z-z0);
                 dydz=(x3-x1)/(z3-z1);
                 y=x3 + dydz*(z-z3);
                 p=(h2[k])->getP();
                 Pmom=sqrt(p*p);
                 FillCovMx(BmResol*BmChrg/BmMoment);
                 CsBeam* pntr = new CsBeam;
                 CsHelix tmp( x, y, z, dxdz, dydz, BmChrg/Pmom, cov );
                 pntr->addHelix(tmp);
                 tmp.~CsHelix();
                 pntr->setChi2(0.);
                 bool Htrig=(-0.0001<timebm)&&(timebm<0.0001);     // set gate for Beam track
                 pntr->setTrigLabel(Htrig);
                 pntr->setBmdata(timebm, 0, 0, 0., 0., 0., 4, 4, timebm, 0.);
                 pntr->addCluster( *clh0[i] );
                 pntr->addCluster( *clh1[j] );
                 pntr->addCluster( *clh2[k] );
                 pntr->addCluster( *clh3[l] );
                 pntr->setZones( BmHod.getZones());
                 bmtracks.push_back(pntr);
//	
                 HiBmom->Fill(Pmom);
                 HiXprof->Fill(x); 
                 HiYprof->Fill(y); 
                 HiXincl->Fill(dxdz); 
                 HiYincl->Fill(dydz); 
                 HiTime->Fill(timebm); 
             }
          }
       }
   }
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CsBeamRecons::BmsGeometryCheck()
{
    BackPropagated back_prop;
    list<CsHelix>::iterator iter;
    int exp_fired = 0;
    int old_exp_fired = 0;
    bool track_ok = false;
    bool old_track_ok = false;
    bool BMS5_used = false;
    bool BMS6_used = false;
    vector<int> dimensionX(6);
    vector<int> dimensionY(6);
    vector<double> offsetY(6);
    vector<bool> exp_fired_planes(6);
    dimensionX[0] = 30;
    dimensionX[1] = 60;
    dimensionX[2] = 60;
    dimensionX[3] = 30;
    dimensionX[4] = 60;
    dimensionX[5] = 60;
    dimensionY[0] = 90;
    dimensionY[1] = 45;
    dimensionY[2] = 50;
    dimensionY[3] = 115;
    dimensionY[4] = 80;
    dimensionY[5] = 80;
    offsetY[0] = 4.69;
    offsetY[1] = 2.065;
    offsetY[2] = -2.266;
    offsetY[3] = -5.342;
    offsetY[4] = 4.524;
    offsetY[5] = -2.904;

    for(iter = goodSFhelixes.begin(); iter != goodSFhelixes.end(); iter++) {
	exp_fired = 0;
	old_exp_fired = 0;
	back_prop.Extrapolate(fabs(1/((*iter).getCop())), BmResol, *iter);//do backpropagation
	for(int i=0; i<6; i++) {//check if extrapolated track hits some bms planes
	    exp_fired_planes[i]=false;

	    if(HiExBmsXY[i] != NULL)
	    {
		HiExBmsXY[i]->Fill(back_prop.GetX(i),back_prop.GetY(i));
	    }
	    if(i<3) {
		if(HiExBmsCorelX[i] != NULL) {
		    HiExBmsCorelX[i]->Fill(back_prop.GetX(i),back_prop.GetX(i+1));
		}
		if(HiExBmsCorelY[i] != NULL) {
		    HiExBmsCorelY[i]->Fill(back_prop.GetY(i),back_prop.GetY(i+1));
		}
	    }
	    if(i<2) {
		if(HiExBmsCorelX[i+3] != NULL) {
		    HiExBmsCorelX[i+3]->Fill(back_prop.GetX(i),back_prop.GetX(4));
		}
		if(HiExBmsCorelY[i+3] != NULL) {
		    HiExBmsCorelY[i+3]->Fill(back_prop.GetY(i),back_prop.GetY(4));
		}
		if(HiExBmsCorelX[i+5] != NULL) {
		    HiExBmsCorelX[i+5]->Fill(back_prop.GetX(i+2),back_prop.GetX(5));
		}
		if(HiExBmsCorelY[i+5] != NULL) {
		    HiExBmsCorelY[i+5]->Fill(back_prop.GetY(i+2),back_prop.GetY(5));
		}
	    }


	    if (
		    fabs(back_prop.GetX(i)) < dimensionX[i] &&
		    fabs(back_prop.GetY(i)-offsetY[i]) < dimensionY[i]
	       )
	    {
		exp_fired++;
		if(i<4) old_exp_fired++;
		exp_fired_planes[i]=true;
		if(HiExBmsEff != NULL) HiExBmsEff->Fill(i);
	    }
	}

	//check if track could be reconstructed in BMS
	if(exp_fired > 1) {//we need at least 2 hits in bms
	    if(exp_fired > 3) {//when we have 4 hits it is posible to reconstruct a track
		track_ok = true;
	    }
	    else if(exp_fired == 2) {//only rescue here
		if (
			(exp_fired_planes[0] || exp_fired_planes[1]) &&
			(exp_fired_planes[2] || exp_fired_planes[3])
		   )
		{//two old planes on different sides of the magnet
		    track_ok = true;
		}
		else if(exp_fired_planes[0] && exp_fired_planes[1])
		{//two old planes before the magnet (is this true that this plus slope from SF/SI gives a track?)
		    track_ok = true;
		}
	    }
	    else if(old_exp_fired > 1)
	    {//here are 3 hit tracks - we need at least 2 old planes
		if (
			(exp_fired_planes[0] || exp_fired_planes[1]) &&
			(exp_fired_planes[2] || exp_fired_planes[3])
		   )
		{//two old plane on different sides of the magnet - it's not important where is the third
		    track_ok = true;
		}
		else if(exp_fired_planes[0] && exp_fired_planes[1])
		{//two old planes before the magnet - even when the third plane is BM05 we still can use rescue (is this true?)
		    track_ok = true;
		}
		else if(exp_fired_planes[2] && exp_fired_planes[3] && exp_fired_planes[4])
		{//BM03 + BM04 + BM05
		    track_ok = true;
		}
	    }

	    if(old_exp_fired  > 2) old_track_ok = true; //3 old planes

	    if (
		    old_exp_fired == 2 &&
		    (exp_fired_planes[0] + exp_fired_planes[1]) != 2 &&
		    exp_fired_planes[4]
	       )
	    {
		BMS5_used=true;
	    }

	    if (
		    old_exp_fired == 2 &&
		    (exp_fired_planes[2] + exp_fired_planes[3]) != 2 &&
		    exp_fired_planes[5] &&
		    !exp_fired_planes[4]
	       )
	    {
		BMS6_used=true;
	    }

	    if (
		    old_exp_fired == 2 &&
		    (exp_fired_planes[0] + exp_fired_planes[1]) == 2 &&
		    exp_fired_planes[5]
	       )
	    {
		BMS6_used=true;
	    }
	}

	//fill histos for tracks
	if(HiExBmsFired != NULL) HiExBmsFired->Fill(exp_fired);
	if(HiExBmsOldFired != NULL) HiExBmsOldFired->Fill(old_exp_fired);
    }
    if(HiExBmsStuff != NULL && track_ok) HiExBmsStuff->Fill(0);//event has a track that could be detected by BMS
    if(HiExBmsStuff != NULL && old_track_ok) HiExBmsStuff->Fill(1);//event has a track that could be detected by old BMS planes
    if(HiExBmsStuff != NULL && BMS5_used && !old_track_ok) HiExBmsStuff->Fill(2);//event has a track that could be detected by using BMS5
    if(HiExBmsStuff != NULL && BMS6_used && !old_track_ok && !BMS5_used) HiExBmsStuff->Fill(3);//event has a track that could be detected by using BMS6
}
