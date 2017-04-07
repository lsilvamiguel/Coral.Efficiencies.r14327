// $Id: CsBMSrecons.cc 13010 2011-11-11 14:56:33Z psznajde $

/*!
   \file CsBMSrecons.cc
   \brief Compass BMS reconstruction Class.
   \author  G. Khaustov
   \version $Revision: 13010 $
   \date    $Date: 2011-11-11 15:56:33 +0100 (Fri, 11 Nov 2011) $
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "CsTypes.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsBMSrecons.h"
#include "CsErrLog.h"
#include "CsOpt.h"
#include "DaqDataDecoding/DaqEvent.h"

#include <stdio.h>

using namespace std;

//
//
CsBMSrecons:: CsBMSrecons(void) 
{
        int ipl, ll, i, j, lbl;
        string tag, key;
        CsOpt* opt = CsOpt::Instance();
        tag="BeamRecons";
//
//  BMS hodoscope description
//
	//set nhod according to year-info
	//2001,2002: 4
	//2003: 5
	//2004: 6
	tag="BeamRecons";
	key="year";
	if( opt->getOpt( tag, key, year)) cout<<"CsBMSrecons: year: "<<year<<endl;
	else year=2001;
	if(year < 2003) nhod=4;
        else if(year < 2004) nhod=5;
	else nhod=NBMSHODS;
	cout<<"CsBMSrecons: Number of BMS-hodoscopes expected: "<<nhod<<endl;
	
	nhod_old=4;
	tag="BMSreco";
	key="useBM05";
	int optBM05;
	if( opt->getOpt( tag, key,optBM05)){
	  if(optBM05==0)useBM05=false;
	  else{cout<<"CsBMSrecons: BM05 will be used!"<<endl; useBM05=true;}
	}
	else useBM05=false;
	//should be: useBM05->nhod=5
	//
	tag="BMSreco";
	key="useBM06";
	int optBM06;
	if( opt->getOpt( tag, key,optBM06)){
	  if(optBM06==0)useBM06=false;
	  else{cout<<"CsBMSrecons: BM06 will be used!"<<endl; useBM06=true;}
	}
	else useBM06=false;
	//
	//
        hodnames[0]="BMSHod0";
        hodnames[1]="BMSHod1";
        hodnames[2]="BMSHod2";
        hodnames[3]="BMSHod3";
	hodnames[4]="BMSHod4";//in 2003 BM05
	hodnames[5]="BMSHod5";//in 2004 BM06
        for (i=0; i<NBMSHODS; i++) hodsize[i]=64;            // size of hods 
	hodsize[5]=128; // BM06 has 128 channels
        BMStbl();
	//checkBMStbl();
//
        tag="BMSreco";
        key="Printlev";
        if( opt->getOpt( tag, key, Printlev ) ) ; else Printlev=0;
//
//  BMS calibrations    
//
        itrg16=0;                        // Temporary. Should be replaced in a future with trigger type
        for(i=0; i<10; i++)  BMSMN1[i]=0.;
        TDCSLP=1.;
//
        nrunold=-1;  
//
//       Reconstruction Procedure parameters
//
        ntrackmx=NTRACKMX;
//
	tag="BeamRecons";
	key="useTRAFFIC";
	if( opt->getOpt( tag, key, useTRAFFIC)) ; else useTRAFFIC=0;
	cout<<"useTraffic: "<<useTRAFFIC<<endl;

        tag="BMSreco";
        key="TmWindow";
        if( opt->getOpt( tag, key, BMSpar[0]) ) ; else BMSpar[0]=10.;
//
        key="MxTmDiff";
        if( opt->getOpt( tag, key, BMSpar[1]) ) ; else BMSpar[1]=2.;
//
        key="MxTmChi4";
        if( opt->getOpt( tag, key, BMSpar[2]) ) ; else BMSpar[2]=30.;
//
        key="MxSpChiq";
        if( opt->getOpt( tag, key, BMSpar[3]) ) ; else BMSpar[3]=10.;
//
        key="MxTmChi3";
        if( opt->getOpt( tag, key, BMSpar[4]) ) ; else BMSpar[4]=30.;
//
        key="MxNmbHit";
        if( opt->getOpt( tag, key, mxnmbhit) ) ; else mxnmbhit=BMSBUFSIZE;
        if(mxnmbhit>BMSBUFSIZE) mxnmbhit=BMSBUFSIZE;
//
	//Sorry Yann, I like my version better:
	/*
	  if (!(opt->getOpt(tag,key,coeff_path)))
          CsErrLog::mes(elFatal,
	  "Path for BMS-Coeff-File is needed in options!");
	  opt->expand(coeff_path);
	  CsErrLog::msg(elInfo,__FILE__,__LINE__,
	  "The path for the CoeffFile is set to: %s",
	  coeff_path.c_str());
	*/
	//key="CoeffFile";
	//opt->getOpt(tag,key,coeff_path);
	//if(!(opt->getOpt(tag,key,coeff_path))) CsErrLog::Instance()->mes(elFatal," Path for BMS-Coeff-File is needed in options!");
	//cout<<"The path for the CoeffFile is set to: "<<coeff_path<<endl;
	char* tmp = getenv("CORAL");
	string coral = tmp;
	//char filepath[200];
	//sprintf(filepath,"%s/src/beam/M2_20.coeff",coral);
	//sprintf(coeff_path,"%s/src/beam/M2_20.coeff",coral);
	//sprintf(coeff_path,"%s/src/beam/M2_20_2004_rr.coef",coral);
	//sprintf(coeff_path1,"%s/src/beam/M2_BM05_1.coeff",coral);
	//sprintf(coeff_path2,"%s/src/beam/M2_BM05_2.coeff",coral);
	//cout<<"BMSrecons: The path for CoeffFile is set to "<<coeff_path<<endl;
	//cout<<"BMSrecons: The path for CoeffFile1 is set to "<<coeff_path1<<endl;
	//cout<<"BMSrecons: The path for CoeffFile2 is set to "<<coeff_path2<<endl;

	//Conrad 2005
	//Let's try one more time :)
	//Using new coeff files
	//rr  for 1,2  and 3,4
	//rc  for 1,2  and 3,6
	//rd  for 1,2  and 6,4
	//
	//ar  for 5,2  and 3,4
	//ac  for 5,2  and 3,6
	//ad  for 5,2  and 6,4
	//
	//br  for 1,5  and 3,4
	//bc  for 1,5  and 3,6
	//bd  for 1,5  and 6,4
	coeff_path[0]=coral + "/src/beam/rec_method_0/M2_20_2004_rr.coef";
	coeff_path[1]=coral + "/src/beam/rec_method_0/M2_20_2004_rc.coef";
	coeff_path[2]=coral + "/src/beam/rec_method_0/M2_20_2004_rd.coef";
	coeff_path[3]=coral + "/src/beam/rec_method_0/M2_20_2004_ar.coef";
	coeff_path[4]=coral + "/src/beam/rec_method_0/M2_20_2004_ac.coef";
	coeff_path[5]=coral + "/src/beam/rec_method_0/M2_20_2004_ad.coef";
	coeff_path[6]=coral + "/src/beam/rec_method_0/M2_20_2004_br.coef";
	coeff_path[7]=coral + "/src/beam/rec_method_0/M2_20_2004_bc.coef";
	coeff_path[8]=coral + "/src/beam/rec_method_0/M2_20_2004_bd.coef";


	
	tag="BMSreco";
	key="RecMethod";
        if( opt->getOpt( tag, key, method) ) ; else method=1;
	method3=false;
	if(method==3){
	  method3=true;
	  cout<<"Method3 is used!"<<endl;
	  //Not used any more. We need this cut to reduce combinatorics.
	  //It should be wide though.
	  //BMSpar[0]=9999999;  //open the time-window for digit selection
	}
	tag="BMSreco";
	key="doClustering";
        if( opt->getOpt( tag, key, clusterIt) )cout<<"CsBMSrecons: Clustering will be done!"<<endl; 
	else {
	    clusterIt=0;
	    //new planes not implemanted for no clustering code
	    nhod = 4;
	    useBM05 = false;
	    useBM06 = false;
	}
	
	tag="BMSreco";
	key="Cluster_ChnlLimit";
	if( opt->getOpt( tag, key, Cl_ChannelLimit) ) ; else Cl_ChannelLimit=3;
	tag="BMSreco";
	key="Cluster_TmLimit";
	if( opt->getOpt( tag, key, Cl_TimeLimit) ) ; else Cl_TimeLimit=3.5;
	if(clusterIt>0)cout<<"Clustering Parameters -- Time: "<<Cl_TimeLimit<<" y-Coordinate: "<<Cl_ChannelLimit<<endl;
	tag="BMSreco";
	key="CommonHits";
	if( opt->getOpt( tag, key, allowedCommon) ) ; else allowedCommon=1;
	cout<<"Number of common hits allowed on a track: "<<allowedCommon<<endl;
	tag="BeamRecons";
	key="Debug";
	if( opt->getOpt( tag, key, DebugMode) ) ; else DebugMode=0;
	if(DebugMode==1){
	  cout<<"CsBMSrecons: Debug-Mode is switched on!"<<endl;
	}
	tag="BeamRecons";
	key="HistogramLevel";
	if( opt->getOpt( tag, key, histoLevel) ) ; else histoLevel=0;
	//cout<<"CsBMSrecons: histoLevel: "<<histoLevel<<endl;
	tag="BMSreco";
	key="F1_Res";
        if( opt->getOpt( tag, key, TDCSLP) ) ; else TDCSLP=1;

	//get time-cuts for special-trigger analysis:
	tag="BeamRecons";
	key="MxTmdiff";
        if( opt->getOpt( tag, key, SpecTrTmdiff ) ) SpecTrTmdiff=fabs(SpecTrTmdiff);
        else SpecTrTmdiff=2.1;
	key="MxTmdiff_center";
	if( opt->getOpt( tag, key, SpecTrTmdiff_cntr ) ) cout<<"CsBMSrecons::SpecTrTmdiff_cntr :"<<SpecTrTmdiff_cntr<<endl;
        else SpecTrTmdiff_cntr=0;
	SpecTrTmdiff_min=SpecTrTmdiff_cntr-SpecTrTmdiff;
	SpecTrTmdiff_max=SpecTrTmdiff_cntr+SpecTrTmdiff;
	cout<<"BMSrecons special Trigger time-cut: ["<<SpecTrTmdiff_min<<","<<SpecTrTmdiff_max<<"]"<<endl;
	tag="BMSreco";
	key="EffTimeCut";
	double EffTimeCut =0;
	double EffTimeCut_cntr=0;
        if( opt->getOpt( tag, key, EffTimeCut ) ) EffTimeCut=fabs(EffTimeCut);
        else EffTimeCut=1.27;
	key="EffTimeCut_center";
	if( opt->getOpt( tag, key, EffTimeCut_cntr ) ) cout<<"CsBMSrecons::EffTimeCut_cntr :"<<EffTimeCut_cntr<<endl;
        else EffTimeCut_cntr=0;
	EffTimeCut_min=EffTimeCut_cntr-EffTimeCut;
	EffTimeCut_max=EffTimeCut_cntr+EffTimeCut;
	cout<<"BMSrecons Efficiency Time Cut: ["<<EffTimeCut_min<<","<<EffTimeCut_max<<"]"<<endl;
	//
	tag="BMSreco";
	key="BM05MxTmDiff_cntr";
	double BMS05MxTmDiff_cntr;
	if( opt->getOpt( tag, key, BMS05MxTmDiff_cntr ) ) cout<<"CsBMSrecons::BMS05MxTmDiff_cntr :"<<BMS05MxTmDiff_cntr<<endl;
	else BMS05MxTmDiff_cntr=0;
	key="BM05MxTmDiff";
	double BMS05MxTmDiff;
	if( opt->getOpt( tag, key, BMS05MxTmDiff ) ) cout<<"CsBMSrecons::BMS05MxTmDiff:"<<BMS05MxTmDiff<<endl;
	else BMS05MxTmDiff=0;
	if(BMS05MxTmDiff==0 && BMS05MxTmDiff_cntr==0){
	  BMS05TmDiff_min=-BMSpar[1];
	  BMS05TmDiff_max=BMSpar[1];
	}
	else{
	  BMS05TmDiff_min=BMS05MxTmDiff_cntr-BMS05MxTmDiff;
	  BMS05TmDiff_max=BMS05MxTmDiff_cntr+BMS05MxTmDiff;
	  cout<<"CsBMSrecons::BMS05TmDiff: ["<<BMS05TmDiff_min<<","<<BMS05TmDiff_max<<"]"<<endl;
	}
	//read BmResol for rescue-method:
	tag="BeamRecons";
	key="BmResol_1-2";
	if( opt->getOpt( tag, key, BmResol_1_2 ) ) ; else BmResol_1_2=0.013;
	key="BmResol_1-3";
	if( opt->getOpt( tag, key, BmResol_1_3 ) ) ; else BmResol_1_3=0.01;
	key="BmResol_1-4";
	if( opt->getOpt( tag, key, BmResol_1_4 ) ) ; else BmResol_1_4=0.008;
	key="BmResol_2-3";
	if( opt->getOpt( tag, key, BmResol_2_3 ) ) ; else BmResol_2_3=0.014;
	key="BmResol_2-4";
	if( opt->getOpt( tag, key, BmResol_2_4 ) ) ; else BmResol_2_4=0.015;
	if(Printlev>0) {
	  cout<<"Resolutions for Rescue-Algo: "<<endl;
	  cout<<"Plane 1&2: "<<BmResol_1_2<<endl;
	  cout<<"Plane 1&3: "<<BmResol_1_3<<endl;
	  cout<<"Plane 1&4: "<<BmResol_1_4<<endl;
	  cout<<"Plane 2&3: "<<BmResol_2_3<<endl;
	  cout<<"Plane 2&4: "<<BmResol_2_4<<endl;
	}
        if(Printlev>0) {
	cout << "CsBMSrecons:: Time window   " << BMSpar[0] << endl;
	cout << "CsBMSrecons:: Max hit time deviation   " << BMSpar[1] << endl;
	cout << "CsBMSrecons:: Max time Chiq for 4 hit tracks  " 
             <<  BMSpar[2] << endl;
	cout << "CsBMSrecons:: Max time Chiq for 3 hit tracks  " 
             <<  BMSpar[4] << endl;
	cout << "CsBMSrecons:: Max. Space Chiq   " 
             <<  BMSpar[3] << endl;
	cout << "CsBMSrecons:: Max. number of Hits/plane   " 
             <<  mxnmbhit << endl;
        }

	// required to avoid MC breakage, cf. emails on coral-weekly at 2011-10-25
	inithist();
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: getBMStrack(const int n, double &pmom, double &pchiq, 
                         double &mntime, double &tmchiq, int &nhit ,
                         int& nhttot, int &spec) const
{
          if((n<0)||(n>=ntrack)) return;
              pmom=TRKp[n];
              mntime=TRKtime[n];
              pchiq=SPchiq[n];
	      if(TRKhit[n]==5){//this means that BM05 was used!
		nhit=3;
		spec=1;//means BM05 has been used
	      }
	      else if(TRKhit[n]==6){//this means that BM06 was used!
		nhit=3;
		spec=2;//means BM06 has been used
	      }
	      else{
		nhit = TRKhit[n];
		spec=0; //means: nothing special
	      }
              tmchiq=TMchiq[n];
              nhttot=nhittotal;
              return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: getBMStrack_resc(const int n, double &pmom, double &pchiq, 
                         double &mntime, double &tmchiq, int &nhit ,
                         int& nhttot, double& BmResol_resc) const
{
          if((n<0)||(n>=ntrack_resc)) return;
              pmom=TRKp_resc[n];
              mntime=TRKtime_resc[n];
              pchiq=SPchiq_resc[n];
              nhit =  TRKhit_resc[n];
              tmchiq=TMchiq_resc[n];
              nhttot=nhittotal;
	      BmResol_resc=TRK_resol[n];
              return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBMSrecons:: getBMSThetaY(const int n, double &thetaY) const
{
  //returns thetaY calculated from Stations 3/4. If these stations are not available
  //(for some 3hit-tracks np is returned with 0.
  int nP=0;
          if((n<0)||(n>=ntrack)){
	    thetaY=99;
	    return nP;
	  }
	  if(TRKhit[n]==4) {
	    //thetaY=atan(-1*(TRKy[n][3]-TRKy[n][2])/12423);
	    thetaY=atan((TRKy[n][3]-TRKy[n][2])/12423);
	    nP=4;
	  }
	  else{
	    if(TRKy[n][3]<99999 && TRKy[n][2]<99999){
	      //thetaY=atan(-1*(TRKy[n][3]-TRKy[n][2])/12423);
	      thetaY=atan((TRKy[n][3]-TRKy[n][2])/12423);
	      nP=3;
	    }
	    else nP=0;
	  }
	  return nP;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBMSrecons:: getBMSY(const int n, double Y[]) const
{
  //returns Y-coordinate from Station1,2,3,4. 
  //nP: -1 -> wrong index
  //    0  -> all "old" planes have fired
  //   >0 && <5  -> missing "old" plane
  //    5  -> more than 1 missing "old" plane
  int nP=-1;
  if((n<0)||(n>=ntrack)){
    Y[0]=Y[1]=Y[2]=Y[3]=9999;
    return nP;
  }
  //loop over "old" planes:
  int nMissing=0;
  for(int i=0;i<4;i++){
    if(TRKy[n][i]<99999)Y[i]=-TRKy[n][i];
    else{
      Y[i]=99999;
      nP=i+1;
      nMissing++;
    }
  }
  if(nMissing==0)nP=0;
  if(nMissing>1)nP=5;
  return nP;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double CsBMSrecons::getBMSmom(const int n) const 
{
  double pmom;
  if((n<0)||(n>=ntrack)) return -1;
  pmom=TRKp[n];
  return pmom;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBMSrecons:: getBMS_HitInfo(const int n, double Time[], int Channel[]) const
{
  int nP=-1;
  if((n<0)||(n>=ntrack)){
    Time[0]=Time[1]=Time[2]=Time[3]=9999;
    Channel[0]=Channel[1]=Channel[2]=Channel[3]=9999;
    return nP;
  }
  //loop over "old" planes:
  int nMissing=0;
  for(int i=0;i<4;i++){
    int index=TRKindx[n][i];
    //cout<<"BMS-Hit: "<<i<<" Channel: "<<index<<endl;
    if(index<9999){
      Channel[i]=ibmz16[index][i];
      Time[i]=tbmz16[index][i];
    }
    else{
      Time[i]=99999;
      Channel[i]=99999;
      nP=i+1;
      nMissing++;
    }
  }
  return nP;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: getBMSmoments(const int n, double& m0, double& m1) const
{
int ipl, k;
          if((n<0)||(n>=ntrack)) return;
            for(m0=0., m1=0., ipl=0; ipl<nhod; ipl++){
              k = TRKindx[n][ipl];
              if(k<hodsize[ipl])  {
                  m1+=tbmz16[k][ipl]/tmresol[ipl];
                  m0+=1./tmresol[ipl];
              }
            }
            return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double CsBMSrecons:: getNewChi2(const int n, const double mntime) const
{
int ipl, k;
double chiq, temp;
          if((n<0)||(n>=ntrack)) return 0.;
            for(chiq=0., ipl=0; ipl<nhod; ipl++){
              k = TRKindx[n][ipl];
              if(k<hodsize[ipl])  {
                  temp=tbmz16[k][ipl] - mntime;
                  chiq+=temp*temp/disp[ipl];
              }
            }
            return chiq;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: bmsprnt(ostream &out) const
{
      int ipl,i;
      out << "CsBMSrecons: " << ntrack << "  BMS track(s) found" << endl;
      if(ntrack) {
        out <<"  #      P   TMchiq  SPchiq    Nhits   Time" << endl;
	out << resetiosflags(ios::scientific);
        out << setiosflags(ios::fixed|ios::showpoint);
        out << setprecision(1);
        for(i=0; i<ntrack; i++)
          out  << setw(3) << i << setw(8) << TRKp[i] << setw(8) << TMchiq[i]  
                << setw(8) << SPchiq[i] << setw(8) << TRKhit[i] << setw(8) 
                << TRKtime[i] << endl;
      }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBMSrecons:: getBMSMult(int plane)
{
  if(plane>=0 && plane<4) return khit16[plane];
  else return 0;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: testRescueMethod(int TrackId, double slopeTg)
{
  if((TrackId<0)||(TrackId>=ntrack)) return;
  if(TRKhit[TrackId]!=4) return;
  //test the rescue algorithm:
  //remove 2 hits from a track and compare resulting momentum:
  double z[4];
  double p_resc,chiq;
  double p=TRKp[TrackId];
  int result;
  for(int ipl=0; ipl<4; ipl++) z[ipl]=TRKy[TrackId][ipl];
  //remove plane 3,4 (keep 1,2):
  z[2]=z[3]=1000000;
  result=moment_resc(z, slopeTg, &p_resc, &chiq);
  if(RescueDeltaMom!=NULL && result!=0) RescueDeltaMom->Fill(p-p_resc,0);
  if(RescueMom2H!=NULL && result!=0)RescueMom2H->Fill(p_resc,0);
  //
  for(int ipl=0; ipl<4; ipl++) z[ipl]=TRKy[TrackId][ipl];
  //remove plane 2,4 (keep 1,3):
  z[1]=z[3]=1000000;
  result=moment_resc(z, slopeTg, &p_resc, &chiq);
  if(RescueDeltaMom!=NULL && result!=0) RescueDeltaMom->Fill(p-p_resc,1);
  if(RescueMom2H!=NULL && result!=0)RescueMom2H->Fill(p_resc,1);
  //
  for(int ipl=0; ipl<4; ipl++) z[ipl]=TRKy[TrackId][ipl];
  //remove plane 2,3 (keep 1,4):
  z[1]=z[2]=1000000;
  result=moment_resc(z, slopeTg, &p_resc, &chiq);
  if(RescueDeltaMom!=NULL && result!=0) RescueDeltaMom->Fill(p-p_resc,2);
  if(RescueMom2H!=NULL && result!=0)RescueMom2H->Fill(p_resc,2);
  //
  for(int ipl=0; ipl<4; ipl++) z[ipl]=TRKy[TrackId][ipl];
  //remove plane 1,4 (keep 2,3):
  z[0]=z[3]=1000000;
  result=moment_resc(z, slopeTg, &p_resc, &chiq);
  if(RescueDeltaMom!=NULL && result!=0) RescueDeltaMom->Fill(p-p_resc,3);
  if(RescueMom2H!=NULL && result!=0)RescueMom2H->Fill(p_resc,3);
  //
  for(int ipl=0; ipl<4; ipl++) z[ipl]=TRKy[TrackId][ipl];
  //remove plane 2,3 (keep 2,4):
  z[1]=z[2]=1000000;
  result=moment_resc(z, slopeTg, &p_resc, &chiq);
  if(RescueDeltaMom!=NULL && result!=0) RescueDeltaMom->Fill(p-p_resc,4);
  if(RescueMom2H!=NULL && result!=0)RescueMom2H->Fill(p_resc,4);

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: bmsrec(void)
{
  int n, ipl, k;

  //make sure that there's no garbage in TRKindx table
  for(n=0; n<NTRACKMX; n++)
      for(ipl=0; ipl<NBMSHODS; ipl++)
	  TRKindx[n][ipl]=9999;

#ifdef BMS_IDEAL_MC
  if (CsEvent::Instance()->isAMonteCarloEvent()) {
    ntrack = 1;
    for(ipl=0;ipl<nhod;ipl++) {
      TRKindx[0][ipl]=0;
      TRKy[0][ipl] = 0;
      khit16[ipl]=0;
    }
    TRKp[0] = 1; 
    TMchiq[0] = 0;
    SPchiq[0] = 0;
    TRKhit[0] = 4;
    TRKtime[0] = 0;
    return;
  }
#endif

  //cout<<"This is CsBMSrecons:: bmsrec "<<endl;
  double dt=-1;

  //
  // start of run and start of job Init
  //
  long int nrun = CsEvent::Instance()->getRunNumber();
  if ( nrun != nrunold ) {
    nrunold = nrun;
    
    inithist();
    if(histoLevel==2) bookhist();
    if(histoLevel==1) bookhist1();
    init();
    cout<<"Reading constants"<<endl;
    if(!readcalib()) CsErrLog::Instance()->mes(elError," Error to read calibrations"); 
  }

//
//    decode and select data
//
      ntrack=0;
      ntrack_resc=0;
      nfired=0;
//      test();
      if(Printlev>1) rawevprint(cout);
      if(clusterIt==1){
	if(!decodeCl()) return;
      }
      else{
	if(!decode()) return;
      }
      Efficiency(0,0);
      if(BMSPlanesfired!=NULL) BMSPlanesfired->Fill(nfired);

      //This one is for sure wrong
      //It has small impact though
      //Without it coral segfaults - wrong design somewhere :-(
      if(nfired<3) {
	if(NuTrpEv!=NULL) NuTrpEv->Fill(0);
	if(mH1[2]!=NULL) mH1[2]->Fill(0);
	if(BMSPlanesfiredNT!=NULL) BMSPlanesfiredNT->Fill(nfired);
	return;
      }
//
//       Fill some test hists
//
      double z[4], ptmp, chiq;
      if((nfired==4)&&
                (khit16[0]*khit16[1]*khit16[2]*khit16[3]==1)) {   // only for 1hit/hod events
	for(ipl=0; ipl<nhod_old; ipl++) {  
	  if(HiHtime1[ipl]!=NULL) HiHtime1[ipl]->Fill(tbmz16[0][ipl]);
	  z[ipl] = zhit16[0][ipl];
	}
	moment4(z, &ptmp, &chiq);                            //beam momentum, chiq
	if(chiq < BMSpar[3]) {    //compare with max spacechiq
	  for(ipl=0; ipl<nhod_old; ipl++) {  
	    if(HiHtime2[ipl]!=NULL) HiHtime2[ipl]->Fill(tbmz16[0][ipl]);
	    if(HiTmVsHt[ipl]!=NULL) HiTmVsHt[ipl]->Fill((double)ibmz16[0][ipl],tbmz16[0][ipl]);
	  } 
	  if(Hit0t1!=NULL)Hit0t1->Fill(tbmz16[0][0],tbmz16[0][1]);
	  if(Hit0t2!=NULL)Hit0t2->Fill(tbmz16[0][0],tbmz16[0][2]);
	  if(Hit0t3!=NULL)Hit0t3->Fill(tbmz16[0][0],tbmz16[0][3]);
	  if(Hit0t1dif!=NULL)Hit0t1dif->Fill(tbmz16[0][0]-tbmz16[0][1]);
	  if(Hit0t2dif!=NULL)Hit0t2dif->Fill(tbmz16[0][0]-tbmz16[0][2]);
	  if(Hit0t3dif!=NULL)Hit0t3dif->Fill(tbmz16[0][0]-tbmz16[0][3]);
	  
	}
      } //nfired == 4, 1 hit per hod.
      //
//    find BMS tracks
//

      //cout<<"CsBMSrecons: doing track searching..."<<endl;
      track();
      //cout<<"CsBMSrecons: done with track searching..."<<endl;


      if(Printlev>1) bmsprnt(cout);
//
//    determine BMS time resolution
//
      if(DebugMode!=0) TimeRes();
      if(P1Clus!=NULL && P4Clus!=NULL) checkClustering();
//
//       Fill output hists
//
      if(NuTrpEv!=NULL) NuTrpEv->Fill(ntrack);
      if(mH1[2]!=NULL) mH1[2]->Fill(ntrack);
      if(ntrack==0 && BMSPlanesfiredNT!=NULL) BMSPlanesfiredNT->Fill(nfired);
      if(!ntrack) return;
      if(ntrack>1 && DebugMode!=0){
	for(int trI=0; trI<ntrack; trI++){
	  //record y-pos of tracks in BMS01 vs BMS02 etc:
	  if(TRKhit[trI]==4 && HitBMS02VsHitBMS01y!=NULL && HitBMS03VsHitBMS02y!=NULL && HitBMS04VsHitBMS03y!=NULL) {
	    HitBMS02VsHitBMS01y->Fill(TRKy[trI][0],TRKy[trI][1]);
	    HitBMS03VsHitBMS02y->Fill(TRKy[trI][1],TRKy[trI][2]);
	    HitBMS04VsHitBMS03y->Fill(TRKy[trI][2],TRKy[trI][3]);
	  }
	  if(TRKhit[trI]==3 &&HitBMS02VsHitBMS01y3!=NULL && HitBMS03VsHitBMS02y3!=NULL && HitBMS04VsHitBMS03y3!=NULL){
	    if(TRKy[trI][0]<99999 && TRKy[trI][1]<99999)HitBMS02VsHitBMS01y3->Fill(TRKy[trI][0],TRKy[trI][1]);
	    if(TRKy[trI][1]<99999 && TRKy[trI][2]<99999)HitBMS03VsHitBMS02y3->Fill(TRKy[trI][1],TRKy[trI][2]);
	    if(TRKy[trI][2]<99999 && TRKy[trI][3]<99999)HitBMS04VsHitBMS03y3->Fill(TRKy[trI][2],TRKy[trI][3]);
	  }
	  double dt_smallest=9999;
	  for(int trJ=trI+1; trJ<ntrack; trJ++){
	    dt=TRKtime[trJ]-TRKtime[trI];
	    if(dt>0 && dt<dt_smallest) dt_smallest=dt;
	  }
	  if(dt_smallest!=9999 && BMSdtTRX!=NULL)BMSdtTRX->Fill(dt);
	}
      }
      int nhit, nhit3 = 0;
        double yhit; 
        for(n=0; n<ntrack; n++) { 
           nhit =  TRKhit[n];
           if(HiNhtev!=NULL) HiNhtev->Fill((double)nhit);
           if(nhit==4){                                             // 4 fired hod. hists
              if(HiBmom4!=NULL) HiBmom4->Fill(TRKp[n]);
	      if(mH1[0]!=NULL) mH1[0]->Fill(TRKp[n]);
              if(HiTMchi4!=NULL) HiTMchi4->Fill(TMchiq[n]);
              if(HiSPchi4!=NULL) HiSPchi4->Fill(SPchiq[n]);
              if(HiTime4!=NULL) HiTime4->Fill(TRKtime[n]);
              for(ipl=0; ipl<nhod_old; ipl++) {
                  k=TRKindx[n][ipl];
                  yhit=(double)ibmz16[k][ipl];
                  if(HiHpro4[ipl]!=NULL) HiHpro4[ipl]->Fill(yhit);
                  if(HiYpro4[ipl]!=NULL) HiYpro4[ipl]->Fill(TRKy[n][ipl]);
              }
              if(khit16[0]*khit16[1]*khit16[2]*khit16[3]==1) {       // 1hit/hod  hists
                   if(HiBmom1!=NULL) HiBmom1->Fill(TRKp[n]);
                   if(HiTMchi1!=NULL) HiTMchi1->Fill(TMchiq[n]);
                   if(HiSPchi1!=NULL) HiSPchi1->Fill(SPchiq[n]);
                   if(HiTime1!=NULL) HiTime1->Fill(TRKtime[n]);
                   for(ipl=0; ipl<nhod_old; ipl++) {
                       k=TRKindx[n][ipl];
                       yhit=(double)ibmz16[k][ipl];
                       if(HiHpro1[ipl]!=NULL) HiHpro1[ipl]->Fill(yhit);
                   }
              }
	      //test the rescue algorithm:
	      //remove 2 hits from a track and compare resulting momentum:
	      double z[4],p_resc;
	      for(ipl=0; ipl<nhod_old; ipl++) z[ipl]=TRKy[n][ipl];
	      //remove plane 3,4 (keep 1,2):
	      z[2]=z[3]=1000000;
	      
	      //remove plane 2,4 (keep 1,3):
	      z[1]=z[3]=1000000;
	      //remove plane 2,3 (keep 1,4):
	      z[1]=z[2]=1000000;
	      //remove plane 1,4 (keep 2,3):
	      z[0]=z[3]=1000000;
	      //remove plane 2,3 (keep 2,4):
	      z[1]=z[2]=1000000;
           }//4-hit-tracks
           else {
             if(HiBmom3!=NULL) HiBmom3->Fill(TRKp[n]);                       
	     if(mH1[1]!=NULL) mH1[1]->Fill(TRKp[n]);
	     //  3 fired hod. hists
	     if(HiTMchi3!=NULL) HiTMchi3->Fill(TMchiq[n]);
	     if(HiSPchi3!=NULL) HiSPchi3->Fill(SPchiq[n]);
	     if(HiTime3!=NULL) HiTime3->Fill(TRKtime[n]);
	     for(ipl=0; ipl<nhod; ipl++) {
	       k=TRKindx[n][ipl];
	       //if(k<hodsize[ipl] && HiHpro3[ipl]!=NULL) HiHpro3[ipl]->Fill((double)ibmz16[k][ipl]);
	       if(k<hodsize[ipl] && BeamH2D[15]!=NULL) BeamH2D[15]->Fill((double)ibmz16[k][ipl],nhod);
	     }
	     nhit3++;
	   }
	   if(nhit==5)if(BeamH1D[1]!=NULL) {
	       BeamH1D[1]->Fill(TRKp[n]);
//	       cout<<n<<':'
//		   <<TRKindx[n][0]<<','
//		   <<TRKindx[n][1]<<','
//		   <<TRKindx[n][2]<<','
//		   <<TRKindx[n][3]<<','
//		   <<TRKindx[n][4]<<','
//		   <<TRKindx[n][5]<<','
//		   <<endl;
	   }
	   if(nhit==6)if(BeamH1D[2]!=NULL) BeamH1D[2]->Fill(TRKp[n]);
	   //check y-correlation between BMS01/02:
	   if(TRKy[n][0]<99999 && TRKy[n][1]<99999) if(BeamH2D[5]!=NULL)BeamH2D[5]->Fill(TRKy[n][0],TRKy[n][1]);
	   //check y-correlation between BMS02/03:
	   if(TRKy[n][1]<99999 && TRKy[n][2]<99999) if(BeamH2D[6]!=NULL)BeamH2D[6]->Fill(TRKy[n][1],TRKy[n][2]);
	   //check y-correlation between BMS03/04:
	   if(TRKy[n][3]<99999 && TRKy[n][2]<99999) if(BeamH2D[7]!=NULL)BeamH2D[7]->Fill(TRKy[n][2],TRKy[n][3]);
	   //check y-correlation between BMS01/04:
	   if(TRKy[n][0]<99999 && TRKy[n][3]<99999) if(BeamH2D[8]!=NULL)BeamH2D[8]->Fill(TRKy[n][1],TRKy[n][3]);
	}
	if(HiN4h3h!=NULL) HiN4h3h->Fill((double)nhit3, (double)(ntrack-nhit3));
	//
	//cout<<"Leaving CsBMSrecons..."<<endl;
	//nhod checked up to here
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: track(void)
//
//      track reconstruction
//
//global variables set in this function:
// ntrack: number of tracks found 
{
  int badDec=0;
      int indx[NBMSHODS],kfl[NBMSHODS],i,j,k,l,ipl,ll,ilb;
      int ktrack, nn, ntot, nstart;
      double delt1, delt2, delt3, delt4, delt5, delt6, z[NBMSHODS], am;
      double p, chiq, tmean, tchiq;
      double t1,t2,t3,t4,t5,t6,delt;
      bool SkipThisComb=false;
      bool SaveThis=true;

      delt1=delt2=delt3=delt4=delt5=delt6=0;
//
//    put some dummies for 3 events with only 3 old planes fired
//    so that loop after hits will execute
//
      if(nfired==3) {
          for(ipl=0; ipl<nhod_old; ipl++) {
             if(!khit16[ipl]) {khit16[ipl]=1; iflc16[0][ipl]=0; }
          }
      }
//
//      four hit tracks
//
        if(nfired==4) {
        for(i=0; i<khit16[0]; i++) {
            indx[0]=i;
            for(j=0; j<khit16[1]; j++){
                indx[1]=j;
                for(k=0; k<khit16[2]; k++){
                    indx[2]=k;
                    for(l=0; l<khit16[3]; l++) {
		      t1=tbmz16[i][0];
		      t2=tbmz16[j][1];
		      t3=tbmz16[k][2];
		      t4=tbmz16[l][3];
		      tmean=(tbmz16[i][0]+tbmz16[j][1]+tbmz16[k][2]+tbmz16[l][3])/4.;
		      //do checks for method3:
		      //check time-differences between hits, if any is greater    
		      //than cut-> skip this combination:
		      SkipThisComb=false;
		      delt=t1-t2;
		      if(fabs(delt)>BMSpar[1]) SkipThisComb=true;
		      if(DeltaT!=NULL)DeltaT->Fill(delt);
		      delt=t1-t3;
		      if(fabs(delt)>BMSpar[1]) SkipThisComb=true;
		      if(DeltaT!=NULL)DeltaT->Fill(delt);
		      delt=t1-t4;
		      if(fabs(delt)>BMSpar[1]) SkipThisComb=true;
		      if(DeltaT!=NULL)DeltaT->Fill(delt);
		      delt=t2-t3;
		      if(fabs(delt)>BMSpar[1]) SkipThisComb=true;
		      if(DeltaT!=NULL) DeltaT->Fill(delt);
		      delt=t2-t4;
		      if(fabs(delt)>BMSpar[1]) SkipThisComb=true;
		      if(DeltaT!=NULL) DeltaT->Fill(delt);
		      delt=t3-t4;
		      if(fabs(delt)>BMSpar[1]) SkipThisComb=true;
		      if(DeltaT!=NULL) DeltaT->Fill(delt);

		      //do checks for method1 &2
		      delt1=tmean-tbmz16[i][0];
		      if(Delta!=NULL) Delta->Fill(delt1);
		      if(!SkipThisComb && DeltaCut!=NULL) DeltaCut->Fill(delt1);
		      delt1=fabs(delt1);
		      if(delt1>BMSpar[1] && !method3) {
			continue;         // timing check
		      }
		      delt2=tmean-tbmz16[j][1];
		      if(Delta!=NULL) Delta->Fill(delt2);
		      if(!SkipThisComb && DeltaCut!=NULL) DeltaCut->Fill(delt2);
		      delt2=fabs(delt2);
		      if(delt2>BMSpar[1] && !method3) {
			continue;
		      }
		      delt3=tmean-tbmz16[k][2];
		      if(Delta!=NULL) Delta->Fill(delt3);
		      if(!SkipThisComb &&  DeltaCut!=NULL) DeltaCut->Fill(delt3);
		      delt3=fabs(delt3);
		      if(delt3>BMSpar[1] && !method3) {
			continue;
		      }
		      delt4=tmean-tbmz16[l][3];
		      if(Delta!=NULL) Delta->Fill(delt4);
		      if(!SkipThisComb && DeltaCut!=NULL) DeltaCut->Fill(delt4);
		      delt4=fabs(delt4);
		      if(delt4>BMSpar[1] && !method3) {
			continue;
		      }
		      
		      if(SkipThisComb && method3) continue;

//
		      tchiq=delt1*delt1+delt2*delt2+delt3*delt3+delt4*delt4;
		      if(tchiq>BMSpar[2]) continue;        // time chiq cut
		      //
		      z[0]=zhit16[i][0];
		      z[1]=zhit16[j][1];
		      z[2]=zhit16[k][2];
		      z[3]=zhit16[l][3];
		      //                      
                       moment4(z,&p,&chiq);                 //beam momentum
                       if(chiq > BMSpar[3]) continue;                      // space chiq cut
//
                       indx[3]=l;
                       ktrack =  ntrack;
   
//
//ckeck for the same momentum/chiq already found - take the best time chiq
//and check for 2 common points - if find take track with better time chiq
		       //
		       if(TrackDec!=NULL)TrackDec->Fill(0);
		       SaveThis=true;
		       int hasOverwritten=-1;
		       bool haskilled=false;
                       if(ntrack>0) {
			 for(ll=0; ll<ntrack; ll++) {    
			   /*
                              if((fabs(p-TRKp[ll])<1.0e-6)
                                 &&(fabs(chiq-SPchiq[ll])<1.0e-6)) {
                                 if(tchiq>TMchiq[ll]) {
				   if(TrackDec!=NULL)TrackDec->Fill(2);
				   if(TchiqDiff!=NULL)TchiqDiff->Fill(fabs(tchiq-TMchiq[ll]));
				   if(fabs(tchiq-TMchiq[ll])<10) badDec++;
				   goto next;   
				 }
				 else{
				   if(TrackDec!=NULL)TrackDec->Fill(3);
				   if(TchiqDiff!=NULL)TchiqDiff->Fill(fabs(tchiq-TMchiq[ll]));
				   if(fabs(tchiq-TMchiq[ll])<10) badDec++;
				 }
                                 ktrack=ll;
                                 goto save;
			      }
			    */
			    for(nn=0, ipl=0;ipl<nhod_old;ipl++){
			      if(indx[ipl]==TRKindx[ll][ipl]){
				nn++;
			      }
			    }
			    if(nn>3){
			      if(TrackDec!=NULL)TrackDec->Fill(5);
			      goto next;			      
			    }
			    if(nn>allowedCommon) {
			      SaveThis=false;//there are common hits so this combination will not enter 
			                      //the list as additional new combination
			      if(tchiq>TMchiq[ll]){//new track is worse
				if(TrackDec!=NULL)TrackDec->Fill(7);
				if(TchiqDiff!=NULL)TchiqDiff->Fill(fabs(tchiq-TMchiq[ll]));
				if(fabs(tchiq-TMchiq[ll])<10) badDec++;
				//goto next; 
				
				if(hasOverwritten>-1){//move ll to the postion of hasOverwritten
				  ReplaceThisTrack(hasOverwritten,ll);
				  for(ipl=0;ipl<nhod;ipl++)TRKindx[ll][ipl]=-1;
				  haskilled=true;
				}
				else goto next; 
				
			      }
			      else{//new comb is better
				if(TrackDec!=NULL)TrackDec->Fill(8);
				if(TchiqDiff!=NULL)TchiqDiff->Fill(fabs(tchiq-TMchiq[ll]));
				if(fabs(tchiq-TMchiq[ll])<10) badDec++;
				//SaveThisTrack(4, ll, indx,kfl,tmean,tchiq,chiq,p,z);
				//goto next;
				
				if(hasOverwritten==-1){
				  SaveThisTrack(4, ll, indx,kfl,tmean,tchiq,chiq,p,z);
				  hasOverwritten=ll;
				}
				else{//delete the current one:
				  //cout<<"Killing "<<ll<<" Time("<<TRKtime[ll]<<") Index(";
				  //for(int u=0;u<4;u++)cout<<TRKindx[ll][u]<<",";
				  //cout<<")"<<endl;
				  for(ipl=0;ipl<nhod;ipl++)TRKindx[ll][ipl]=-1;
				  haskilled=true;
				}
				
			      }
			      //ktrack=ll; 
			      //goto save;
			    }
                          }//loop over existing tracks
			  //remove the deleted tracks from the list:
			  if(haskilled){
			    //loop over the list of tracks to find the candidates:
			    for(ll=0;ll<ntrack;ll++){
			      bool removeThis=true;
			      for(ipl=0;ipl<nhod;ipl++)if(TRKindx[ll][ipl]!=-1)removeThis=false;
			      if(removeThis){//move the last track in the list to this postion:
				if(ll==ntrack-1){//this is already the last track
				  ntrack=ntrack-1;
				}
				else{
				  bool exit=false;
				  for(int it=1;it<(ntrack-ll);it++){
				    if(!exit){
				      removeThis=true;
				      for(ipl=0;ipl<nhod;ipl++)if(TRKindx[ntrack-it][ipl]!=-1)removeThis=false;
				      if(!removeThis){//This track is ok and can replace ll:
					ReplaceThisTrack(ll,ntrack-it);
					ntrack=ntrack-1;
					exit=true;
				      }
				      else ntrack=ntrack-1;
				    }
				  }//looking for replacement track
				}
			      }
			    }//look for deleted tracks
			  }//there is a track to be deleted!
			  //check whether any of the deleted tracks can be re-introduced:
			  
                       }//ntrack>0
		       

//
//    save new track
//
		       if(SaveThis){//a new combination is to be saved
			 if(ntrack==NTRACKMX) return;
			 ntrack++;
			 SaveThisTrack(4, ktrack, indx,kfl,tmean,tchiq,chiq,p,z);
		       }
		       /*
		    save:
		    SaveThisTrack(4, ktrack, indx,kfl,tmean,tchiq,chiq,p,z);
		    
		    for(ipl=0;ipl<nhod_old;ipl++) {
		    TRKindx[ktrack][ipl]=indx[ipl];
		    TRKy[ktrack][ipl] = z[ipl];
		    }
                         TRKp[ktrack] = p; 
                         TMchiq[ktrack] = tchiq;
                         SPchiq[ktrack] = chiq;
                         TRKhit[ktrack] = 4;
                         TRKtime[ktrack] = tmean;
			 */
                next:    continue;  
		    }
                 }
             }
	} //outermost for-loop
	//
	
	if(BadDec4!=NULL)BadDec4->Fill(badDec);
	
	if( (khit16[0]*khit16[1]*khit16[2]*khit16[3]) ==1 ){
	  double deltas[4];
	  deltas[0]=delt1;
	  deltas[1]=delt2;
	  deltas[2]=delt3;
	  deltas[3]=delt4;
	  if(Delta1H1!=NULL) Delta1H1->Fill(delt1);
	  if(Delta1H2!=NULL) Delta1H2->Fill(delt2);
	  if(Delta1H3!=NULL) Delta1H3->Fill(delt3);
	  if(Delta1H4!=NULL) Delta1H4->Fill(delt4);
	  if(DeltaG!=NULL){
	    for(i=0;i<4;i++){
	      DeltaG->Fill(deltas[i]);
	    }
	  }
	}
////       set hit flags for track found, so 3 hit tracks procedure will not considere these hits 
//
         if(ntrack) {
             for(ipl=0; ipl<nhod_old; ipl++) {
                 for(i=0; i<ntrack; i++) {
                     ll=TRKindx[i][ipl];
                     iflc16[ll][ipl]=0;
                 }
                 for(kfl[ipl]=0, i=0; i<khit16[ipl]; i++){
                     if(iflc16[i][ipl]) kfl[ipl]=1;}
             }
	     //how many hits are left over? Only for !useBM05
             if(kfl[0]+kfl[1]+kfl[2]+kfl[3]<3) return;
         }
	}// treatment of 4-planes-cases
//
//       3 hit tracks
//
	badDec=0;
        nstart=ntrack;


        for(i=0; i<khit16[0]; i++) {
            for(j=0; j<khit16[1]; j++){
                for(k=0; k<khit16[2]; k++){
                    for(l=0; l<khit16[3]; l++) {
		      //check whether any of these hits has been used before
                       kfl[0]=iflc16[i][0];
                       kfl[1]=iflc16[j][1];
                       kfl[2]=iflc16[k][2];
                       kfl[3]=iflc16[l][3];
		       kfl[4]=0;
		       kfl[5]=0;
		       //
                       t1=tbmz16[i][0];
		       t2=tbmz16[j][1];
		       t3=tbmz16[k][2];
		       t4=tbmz16[l][3];
		       //
                       ntot=kfl[0]+kfl[1]+kfl[2]+kfl[3];
                       //
		       tmean=(tbmz16[i][0]*kfl[0]+tbmz16[j][1]*kfl[1]+
			      tbmz16[k][2]*kfl[2]+tbmz16[l][3]*kfl[3])/ntot;
		       delt1=fabs(tmean-tbmz16[i][0]);
		       delt2=fabs(tmean-tbmz16[j][1]);
		       delt3=fabs(tmean-tbmz16[k][2]);
		       delt4=fabs(tmean-tbmz16[l][3]);
		       delt5=0;
		       delt6=0;
		       //
		       ilb=0;
                       if(ntot==nhod_old){                    
			 //
			 // These are 4 hit tracks - must have been found before and failed
			 // hence remove the worst hit
			 //
                          am=delt1; ilb=0;                
                          if(delt2>am){am=delt2; ilb=1;}  
                          if(delt3>am){am=delt3; ilb=2;}
                          if(delt4>am) ilb=3;
                          kfl[ilb]=0;
                          ntot=3;
                          tmean=(tbmz16[i][0]*kfl[0]+tbmz16[j][1]*kfl[1]+
                                 tbmz16[k][2]*kfl[2]+tbmz16[l][3]*kfl[3])/ntot;
                          delt1=fabs(tmean-tbmz16[i][0]);
                          delt2=fabs(tmean-tbmz16[j][1]);
                          delt3=fabs(tmean-tbmz16[k][2]);
                          delt4=fabs(tmean-tbmz16[l][3]);
                       }
		       //
		       
		       int MinNrHits=3;
		       if(useBM05)MinNrHits=2;
		       if(ntot<MinNrHits)continue;
		       
		       //
		       
		       //for method3 choose hits by cutting on time difference
		       // between them
		       bool SkipThisComb=false;
		       if(method3){
			 if((kfl[0]!=0 && kfl[1]!=0) && fabs(t1-t2)>BMSpar[1]) SkipThisComb=true;
			 if((kfl[0]!=0 && kfl[2]!=0) && fabs(t1-t3)>BMSpar[1]) SkipThisComb=true;
			 if((kfl[0]!=0 && kfl[3]!=0) && fabs(t1-t4)>BMSpar[1]) SkipThisComb=true;
			 if((kfl[1]!=0 && kfl[2]!=0) && fabs(t2-t3)>BMSpar[1]) SkipThisComb=true;
			 if((kfl[1]!=0 && kfl[3]!=0) && fabs(t2-t4)>BMSpar[1]) SkipThisComb=true;
			 if((kfl[2]!=0 && kfl[3]!=0) && fabs(t3-t4)>BMSpar[1]) SkipThisComb=true;
		       }
		       else{
			 if(delt1>BMSpar[1]) kfl[0]=0;         // check timing
			 if(delt2>BMSpar[1]) kfl[1]=0;
			 if(delt3>BMSpar[1]) kfl[2]=0;
			 if(delt4>BMSpar[1]) kfl[3]=0;
			 //let's see how many hits are left
			 ntot=kfl[0]+kfl[1]+kfl[2]+kfl[3];
			 if(ntot<3) continue; //for now no recovery for this method
		       }		       
                       
		       int indxBM05=0;
		       int indxBM06=0;
		       if(SkipThisComb){
			   //this combination didn't corelate - maybe using BM05 or BM06 will help
			 if(!(useBM05 || useBM06))continue;
			 else{
			   //one of the usable combinations are too far apart in time
			   //->remove the one furthest from mean
			   am=0;
			   if(kfl[0]!=0 && delt1>am){ am=delt1; ilb=0;}                
			   if(kfl[1]!=0 && delt2>am){ am=delt2; ilb=1;}  
			   if(kfl[2]!=0 && delt3>am){ am=delt3; ilb=2;}
			   if(kfl[3]!=0 && delt4>am) ilb=3;
			   kfl[ilb]=0;
			   //now we should have two old stations left:
			   ntot=kfl[0]+kfl[1]+kfl[2]+kfl[3];
			   if(ntot>2){
			     CsErrLog::Instance()->mes(elFatal,"Problem in BM05-procedure: too many old stations left after removal!");
			   }
			   //check whether this pair is ok:
			   SkipThisComb=false;
			   if((kfl[0]!=0 && kfl[1]!=0) && fabs(t1-t2)>BMSpar[1]) SkipThisComb=true;
			   if((kfl[0]!=0 && kfl[2]!=0) && fabs(t1-t3)>BMSpar[1]) SkipThisComb=true;
			   if((kfl[0]!=0 && kfl[3]!=0) && fabs(t1-t4)>BMSpar[1]) SkipThisComb=true;
			   if((kfl[1]!=0 && kfl[2]!=0) && fabs(t2-t3)>BMSpar[1]) SkipThisComb=true;
			   if((kfl[1]!=0 && kfl[3]!=0) && fabs(t2-t4)>BMSpar[1]) SkipThisComb=true;
			   if((kfl[2]!=0 && kfl[3]!=0) && fabs(t3-t4)>BMSpar[1]) SkipThisComb=true;
			   //if also this pair is not satisfactory, forget it:
			   if(SkipThisComb)continue;
			   
			   //When we're using BM05 make sure that not both plane3 && plane4 are missing
			   if(useBM05 && (kfl[2]!=0 || kfl[3]!=0)){
			       //Now try to find a corresponding hit in BM05: loop over hits of BM05:
			       double TimeChiMin=-1;
			       int MinIndx=-1;
			       double MinDelt=-1;
			       double MinMean=-1;
			       for(indxBM05=0; indxBM05<khit16[4]; indxBM05++) {
				   //has this hit been used before?
				   //if(iflc16[indxBM05][4]==0) continue;
				   //ignore those too far away from the other 2 hits
				   t5=tbmz16[indxBM05][4];
				   SkipThisComb=false;
				   
				   if((kfl[0]!=0) && (t1-t5<BMS05TmDiff_min || t1-t5>BMS05TmDiff_max)) SkipThisComb=true;
				   if((kfl[1]!=0) && (t2-t5<BMS05TmDiff_min || t2-t5>BMS05TmDiff_max)) SkipThisComb=true;
				   if((kfl[2]!=0) && (t3-t5<BMS05TmDiff_min || t3-t5>BMS05TmDiff_max)) SkipThisComb=true;
				   if((kfl[3]!=0) && (t4-t5<BMS05TmDiff_min || t4-t5>BMS05TmDiff_max)) SkipThisComb=true;
				   if(SkipThisComb)continue;
				   //this BM05 hit fits to the old pair: 
				   kfl[4]=1;
				   //calculate mean and delti
				   ntot=kfl[0]+kfl[1]+kfl[2]+kfl[3]+kfl[4];
				   double res_old=1/(0.42*0.42);
				   double res_BM05=1/(0.99*0.99);
				   double weight=res_old*kfl[0]+res_old*kfl[1]+res_old*kfl[2]+res_old*kfl[3]+res_BM05*kfl[4];	
				   //tmean=(tbmz16[i][0]*kfl[0]+tbmz16[j][1]*kfl[1]+tbmz16[k][2]*kfl[2]+tbmz16[l][3]*kfl[3]+tbmz16[indxBM05][4]*kfl[4])/ntot;
				   tmean=(tbmz16[i][0]*res_old*kfl[0]+tbmz16[j][1]*res_old*kfl[1]+tbmz16[k][2]*res_old*kfl[2]+tbmz16[l][3]*res_old*kfl[3]+tbmz16[indxBM05][4]*res_BM05*kfl[4])/weight;
				   delt1=fabs(tmean-tbmz16[i][0]);
				   delt2=fabs(tmean-tbmz16[j][1]);
				   delt3=fabs(tmean-tbmz16[k][2]);
				   delt4=fabs(tmean-tbmz16[l][3]);
				   delt5=fabs(tmean-tbmz16[indxBM05][4]);
				   //save indxBM05 and time-chi2, if better than what we might have had before:
				   double TimeChi=(delt1*delt1)*kfl[0]+(delt2*delt2)*kfl[1]+(delt3*delt3)*kfl[2]+(delt4*delt4)*kfl[3]+(delt5*delt5)*kfl[4];
				   if(TimeChiMin==-1){
				       TimeChiMin=TimeChi;
				       MinIndx=indxBM05;
				       MinDelt=delt5;
				       MinMean=tmean;
				   }
				   else{
				       if(TimeChi<TimeChiMin){
					   TimeChiMin=TimeChi;
					   MinIndx=indxBM05;
					   MinDelt=delt5;
					   MinMean=tmean;
				       }
				       else continue;
				   }			     
			       }//loop over BM05 hits
			       //keep the one with the best time-chi2
			       //save best iterator of BM05-loop to indxBM05
			       if(MinIndx>=0){
				   indxBM05=MinIndx;
				   //compute tmean and compute delt5
				   tmean=MinMean;
				   delt5=MinDelt;
				   //cout<<"Found hit in BM05"<<endl;
			       }
			   }
			   //When we're using BM06 make sure that not both plane1 && plane2 are missing (and make sure that BM05 was not used for this combination).
			   if(useBM06  && (kfl[0]!=0 || kfl[1]!=0) && kfl[4]==0){
			       //Now try to find a corresponding hit in BM06: loop over hits of BM06:
			       double TimeChiMin=-1;
			       int MinIndx=-1;
			       double MinDelt=-1;
			       double MinMean=-1;
			       for(indxBM06=0; indxBM06<khit16[5]; indxBM06++) {
				   //has this hit been used before?
				   //if(iflc16[indxBM05][4]==0) continue;
				   //ignore those too far away from the other 2 hits
				   t6=tbmz16[indxBM06][5];
				   SkipThisComb=false;

#warning TODO: change BMS05TmDiff_ to BMS06TmDiff_
				   if((kfl[0]!=0) && (t1-t6<BMS05TmDiff_min || t1-t6>BMS05TmDiff_max)) SkipThisComb=true;
				   if((kfl[1]!=0) && (t2-t6<BMS05TmDiff_min || t2-t6>BMS05TmDiff_max)) SkipThisComb=true;
				   if((kfl[2]!=0) && (t3-t6<BMS05TmDiff_min || t3-t6>BMS05TmDiff_max)) SkipThisComb=true;
				   if((kfl[3]!=0) && (t4-t6<BMS05TmDiff_min || t4-t6>BMS05TmDiff_max)) SkipThisComb=true;
				   if(SkipThisComb)continue;
				   
				   
				   //this BM06 hit fits to the old pair: 
				   kfl[5]=1;
				   //calculate mean and delti
				   ntot=kfl[0]+kfl[1]+kfl[2]+kfl[3]+kfl[5];
				   double res_old=1.0/(0.42*0.42);
				   double res_BM06=1.0/(0.99*0.99);
				   double weight=res_old*kfl[0]+res_old*kfl[1]+res_old*kfl[2]+res_old*kfl[3]+res_BM06*kfl[5];
				   //tmean=(tbmz16[i][0]*kfl[0]+tbmz16[j][1]*kfl[1]+tbmz16[k][2]*kfl[2]+tbmz16[l][3]*kfl[3]+tbmz16[indxBM05][4]*kfl[4])/ntot;
				   tmean=(tbmz16[i][0]*res_old*kfl[0]+tbmz16[j][1]*res_old*kfl[1]+tbmz16[k][2]*res_old*kfl[2]+tbmz16[l][3]*res_old*kfl[3]+tbmz16[indxBM06][5]*res_BM06*kfl[5])/weight;
				   delt1=fabs(tmean-tbmz16[i][0]);
				   delt2=fabs(tmean-tbmz16[j][1]);
				   delt3=fabs(tmean-tbmz16[k][2]);
				   delt4=fabs(tmean-tbmz16[l][3]);
				   delt6=fabs(tmean-tbmz16[indxBM06][5]);
				   //save indxBM06 and time-chi2, if better than what we might have had before:
				   double TimeChi=(delt1*delt1)*kfl[0]+(delt2*delt2)*kfl[1]+(delt3*delt3)*kfl[2]+(delt4*delt4)*kfl[3]+(delt6*delt6)*kfl[5];
				   if(TimeChiMin==-1){
				       TimeChiMin=TimeChi;
				       MinIndx=indxBM06;
				       MinDelt=delt6;
				       MinMean=tmean;
				   }
				   else{
				       if(TimeChi<TimeChiMin){
					   TimeChiMin=TimeChi;
					   MinIndx=indxBM06;
					   MinDelt=delt6;
					   MinMean=tmean;
				       }
				       else continue;
				   }			     
			       }//loop over BM06 hits
			       //keep the one with the best time-chi2
			       //save best iterator of BM06-loop to indxBM06
			       if(MinIndx>=0){
				   indxBM06=MinIndx;
				   //compute tmean and compute delt5
				   tmean=MinMean;
				   delt6=MinDelt;
				   //cout<<"Found hit in BM06"<<endl;
			       }
			   }
			   //else continue;

			 }//use BM05 or BM06
		       }//3-hit comb does not fit together

		       //let's see how many hits we now have:
		       ntot=kfl[0]+kfl[1]+kfl[2]+kfl[3]+kfl[4]+kfl[5];
		       if(ntot<3) continue;
//
                       tchiq=(delt1*delt1)*kfl[0]+
                             (delt2*delt2)*kfl[1]+
                             (delt3*delt3)*kfl[2]+
                             (delt4*delt4)*kfl[3]+
			     (delt5*delt5)*kfl[4]+
		             (delt6*delt6)*kfl[5];
                       if(tchiq>BMSpar[4]) continue;        // time chiq cut
//
                       ktrack=ntrack;
                       indx[0]=(i+1)*kfl[0]-1;//in order to catch index=0 
                       indx[1]=(j+1)*kfl[1]-1;
                       indx[2]=(k+1)*kfl[2]-1; 
                       indx[3]=(l+1)*kfl[3]-1;
		       indx[4]=(indxBM05+1)*kfl[4]-1;
		       indx[5]=(indxBM06+1)*kfl[5]-1;
		       SaveThis=true;

                       if(ktrack!=nstart) {
//
// No common points allowed for 3 hit tracks -
// if find - take track with better time chiq
//			 
			 int hasOverwritten=-1;
			 bool haskilled=false;
			 //if(fabs(tmean)<1 && PrintIt) cout<<"Checking "<<ktrack<<" Time("<<tmean<<") Index("<<indx[0]<<","<<indx[1]<<","<<indx[2]<<","<<indx[3]<<")"<<endl;
			 //cout<<"Checking "<<ktrack<<" Time("<<tmean<<") Index("<<indx[0]<<","<<indx[1]<<","<<indx[2]<<","<<indx[3]<<")"<<endl;
			 int ntot_old;
			 for(ll=nstart;ll<ntrack;ll++) {
			   ntot_old=0;
			   ntot=0;
			   for(ipl=0;ipl<nhod;ipl++) {
			     if(kfl[ipl]&&(indx[ipl]==TRKindx[ll][ipl])){
			       ntot++;
			       if(ipl<4)ntot_old++;
 			       SaveThis=false;//there are common hits so this combination will not enter 
			                      //the list as additional new combination
			     }//there is common index
			   }
			   if(ntot==3) goto next3; //no need to check the other combinations
			   //if(fabs(tmean)<1 &&PrintIt)cout<<"Comparison with"<<ll<<": "<<ntot<<endl;
			   //cout<<"Comparison with"<<ll<<": "<<ntot<<endl;
			   if(ntot) {
			     if(ntot_old && nhod==5){
			       //if only one track has a BM05 or BM06 hit the other wins
			       if(kfl[4] && TRKindx[ll][4]==9999) goto next3;//new track has BM05 hit, the already existing doesn't
			       if(!kfl[4] && TRKindx[ll][4]!=9999){
				 if(hasOverwritten==-1){
				   SaveThisTrack(3, ll, indx,kfl, tmean,tchiq,0,0,z);//new track has no BM05 or BM06 hit, the already existin has
				   hasOverwritten=ll;
				 }
				 else{//delete the current one:
//				     cout<<"kill1"<<endl;
				   for(ipl=0;ipl<nhod;ipl++)TRKindx[ll][ipl]=-1;
				   haskilled=true;
				 }
				 continue;
			       }
			     }
			     else if(ntot_old && nhod==6){
			       //if only one track has a BM05 or BM06 hit the other wins
			       //if((kfl[4] && TRKindx[ll][4]==9999) || (kfl[5] && TRKindx[ll][5]==9999)) goto next3;//new track has BM05 hit, the already existing doesn't
			       if( ( (kfl[4] || kfl[5]) && (TRKindx[ll][4]==9999 && TRKindx[ll][5]==9999) )
                                   || (kfl[5] && TRKindx[ll][4])
				 ) goto next3;//new track has BM05 or BM06 hit, the already existing doesn't
			       //if((!kfl[4] && TRKindx[ll][4]!=9999) || (!kfl[5] && TRKindx[ll][5]!=9999)){
			       if( ( (!kfl[4] && !kfl[5]) && (TRKindx[ll][4]!=9999 || TRKindx[ll][5]!=9999) )
                                   || (kfl[4] && TRKindx[ll][5]!=9999)
				 ){
				 if(hasOverwritten==-1){
				   SaveThisTrack(3, ll, indx,kfl, tmean,tchiq,0,0,z);//new track has no BM05 or BM06 hit, the already existin has
				   hasOverwritten=ll;
				 }
				 else{//delete the current one:
				   for(ipl=0;ipl<nhod;ipl++)TRKindx[ll][ipl]=-1;
				   haskilled=true;
				 }
				 continue;
			       }
			     }
			     if(tchiq>TMchiq[ll]){//new combination is worse
			       if(TrackDec!=NULL)TrackDec->Fill(10);
			       if(TchiqDiff!=NULL)TchiqDiff3->Fill(fabs(tchiq-TMchiq[ll]));
			       if(fabs(tchiq-TMchiq[ll])<10) badDec++;
			       //goto next3;
			       if(hasOverwritten>-1){//move ll to the postion of hasOverwritten
				 ReplaceThisTrack(hasOverwritten,ll);
				 //...and mark ll for deletion:
				 //cout<<"Killing "<<ll<<" Time("<<TRKtime[ll]<<") Index(";
				 //for(int u=0;u<4;u++)cout<<TRKindx[ll][u]<<",";
				 //cout<<")"<<endl;
				 for(ipl=0;ipl<nhod;ipl++)TRKindx[ll][ipl]=-1;
				 haskilled=true;
				 
			       }
			       else goto next3;
			     }
			     else{//new combination is better
			       if(TrackDec!=NULL)TrackDec->Fill(11);
			       if(TchiqDiff!=NULL)TchiqDiff3->Fill(fabs(tchiq-TMchiq[ll]));
			       if(fabs(tchiq-TMchiq[ll])<10) badDec++;
			       if(hasOverwritten==-1){
				 SaveThisTrack(3, ll, indx,kfl, tmean,tchiq,0,0,z);
				 hasOverwritten=ll;
			       }
			       else{//delete the current one:
				 //cout<<"Killing "<<ll<<" Time("<<TRKtime[ll]<<") Index(";
				 //for(int u=0;u<4;u++)cout<<TRKindx[ll][u]<<",";
				 //cout<<")"<<endl;
				 for(ipl=0;ipl<nhod;ipl++)TRKindx[ll][ipl]=-1;
				 haskilled=true;
			       }
			     }
                           }//ntot!=0
			 }//loop over existing tracks
			 //remove the deleted tracks from the list:
			 if(haskilled){
			   //loop over the list of tracks to find the candidates:
			   for(ll=nstart;ll<ntrack;ll++){
			     bool removeThis=true;
			     for(ipl=0;ipl<nhod;ipl++)if(TRKindx[ll][ipl]!=-1)removeThis=false;
			     if(removeThis){//move the last track in the list to this postion:
			       if(ll==ntrack-1){//this is already the last track
				 ntrack=ntrack-1;
			       }
			       else{
				 bool exit=false;
				   for(int it=1;it<(ntrack-ll);it++){
				     if(!exit){
				       removeThis=true;
				       for(ipl=0;ipl<nhod;ipl++)if(TRKindx[ntrack-it][ipl]!=-1)removeThis=false;
				       if(!removeThis){//This track is ok and can replace ll:
					 ReplaceThisTrack(ll,ntrack-it);
					 ntrack=ntrack-1;
					 exit=true;
				       }
				       else ntrack=ntrack-1;
				     }
				   }//looking for replacement track
			       }
			     }
			   }//look for deleted tracks
			 }//there is a track to be deleted!
		       }//ktrack!=nstart
		       if(SaveThis && ntrack<NTRACKMX){//a new combination is to be saved
			 ntrack++;
			 SaveThisTrack(3, ktrack, indx,kfl, tmean,tchiq,0,0,z);
		       }
		       
		       /*
               save3: 
		       //SaveThisTrack(3, ktrack, indx,kfl, tmean,tchiq,z);
		       
                       for(ipl=0;ipl<nhod;ipl++) { 
                           if(kfl[ipl]) {
                                   TRKindx[ktrack][ipl]=indx[ipl]; 
                                   z[ipl]=zhit16[indx[ipl]][ipl];
                                   TRKy[ktrack][ipl] = z[ipl];
                           }
                           else  {
                                   z[ipl]=99999.;
                                   TRKindx[ktrack][ipl]=9999; 
                                   TRKy[ktrack][ipl] = z[ipl];
                           }
                       } 
                       TMchiq[ktrack] = tchiq;
                       TRKtime[ktrack] = tmean;  
                       moment3(z,&p,&chiq);
                       TRKp[ktrack] = p; 
                       SPchiq[ktrack] = fabs(p0-p);
                       TRKhit[ktrack] = 3;
		       if(kfl[4])TRKhit[ktrack] = 5;
		       //if(kfl[4])cout<<"Momentum: "<<p<<endl;
		       //if(fabs(tmean)<1)cout<<"Saving: "<<ktrack<<" Time("<<tmean<<") Index(";
		       cout<<"Saving: "<<ktrack<<" Time("<<tmean<<") Index(";
		       for(int u=0;u<4;u++)cout<<TRKindx[ktrack][u]<<",";
		       cout<<")"<<endl;
		       */
               next3:  continue;  
		    }
                }
            }
        } //outermost for-loop for 3-hit tracks
	if(BadDec3!=NULL)BadDec3->Fill(badDec);
	/*
	//loop over tracks to check for common hits:
	for(int I=0;I<ntrack;I++){
	  if(fabs(TRKtime[I])<1){
	    cout<<"BMS track: ID("<<I<<") Mom("<<TRKp[I]<<") SChi("<<SPchiq[I]<<") TChi("<<TMchiq[I]<<") Time("<<TRKtime[I]<<") Index(";
	    for(int P=0;P<4;P++)cout<<TRKindx[I][P]<<",";
	    cout<<")"<<endl;
	  }
	}
	*/
}
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//void CsBMSrecons:: StoreThisTrack(BMScomb &thisComb,int nplanes, int indx[],int kfl[],double tmean,double tchiq,double chiq,double p,double z[]){
void CsBMSrecons:: StoreThisTrack(BMScomb &thisComb,int TrIndex){
  if((TrIndex<0)||(TrIndex>=ntrack)) return;
  thisComb.tmean=TRKtime[TrIndex];
  thisComb.nplanes=TRKhit[TrIndex];
  for(int ipl=0;ipl<nhod;ipl++){
    thisComb.indx[ipl]=TRKindx[TrIndex][ipl];
    if(TRKindx[TrIndex][ipl]<9999)thisComb.kfl[ipl]=1;
    thisComb.z[ipl]=TRKy[TrIndex][ipl];
  }
  thisComb.tchiq=TMchiq[TrIndex];
  thisComb.chiq=SPchiq[TrIndex];
  thisComb.p=TRKp[TrIndex];
  //cout<<"Storing: "<<TRKtime[TrIndex]<<endl;
}
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: ReplaceThisTrack(int ToBeReplaced, int Replacing){
  //cout<<"Replacing: "<<ToBeReplaced<<" Time("<<TRKtime[ToBeReplaced]<<") Index(";
  //for(int u=0;u<4;u++)cout<<TRKindx[ToBeReplaced][u]<<",";
  //cout<<") with "<<Replacing<<" Time("<<TRKtime[Replacing]<<") Index(";
  //for(int u=0;u<4;u++)cout<<TRKindx[Replacing][u]<<",";
  //cout<<")"<<endl;
  if((Replacing<0)||(Replacing>=ntrack)) return;
  for(int ipl=0;ipl<nhod;ipl++) { 
    TRKindx[ToBeReplaced][ipl]=TRKindx[Replacing][ipl]; 
    TRKy[ToBeReplaced][ipl] =TRKy[Replacing][ipl];
  }
  TMchiq[ToBeReplaced] = TMchiq[Replacing];
  TRKtime[ToBeReplaced] = TRKtime[Replacing];  
  TRKp[ToBeReplaced] =TRKp[Replacing]; 
  SPchiq[ToBeReplaced] = SPchiq[Replacing];
  TRKhit[ToBeReplaced] = TRKhit[Replacing];
}
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: SaveThisTrack(int nplanes, int ktrack, int indx[],int kfl[],double tmean,double tchiq,double chiq,double p,double z[]){
  //kfl only needed for 3hit-tracks
  if(nplanes==4){
    for(int ipl=0;ipl<nhod_old;ipl++) {
      TRKindx[ktrack][ipl]=indx[ipl];
      TRKy[ktrack][ipl] = z[ipl];
    }
    TRKp[ktrack] = p; 
    TMchiq[ktrack] = tchiq;
    SPchiq[ktrack] = chiq;
    TRKhit[ktrack] = 4;
    TRKtime[ktrack] = tmean;
    //if(kfl[4])cout<<"Momentum: "<<p<<endl;
    //if(fabs(tmean)<1)cout<<"Saving: "<<ktrack<<" Time("<<tmean<<") Index(";
    //cout<<"Saving: "<<ktrack<<" Time("<<tmean<<") Index(";
    //for(int u=0;u<4;u++)cout<<TRKindx[ktrack][u]<<",";
    //cout<<")"<<endl;
    
  }//4planes involved

  if(nplanes==3){
    for(int ipl=0;ipl<nhod;ipl++) { 
      if(kfl[ipl]) {
	TRKindx[ktrack][ipl]=indx[ipl]; 
	z[ipl]=zhit16[indx[ipl]][ipl];
	TRKy[ktrack][ipl] = z[ipl];
      }
      else  {
	z[ipl]=99999.;
	TRKindx[ktrack][ipl]=9999; 
	TRKy[ktrack][ipl] = z[ipl];
      }
    } 
    TMchiq[ktrack] = tchiq;
    TRKtime[ktrack] = tmean;  
    moment3(z,&p,&chiq);
    TRKp[ktrack] = p; 
    SPchiq[ktrack] = fabs(p0-p);
    TRKhit[ktrack] = nplanes;
    if(kfl[4])TRKhit[ktrack] = 5;
    if(kfl[5])TRKhit[ktrack] = 6;
    //if(kfl[4]) {
//	cout<<"Momentum: "<<p<<endl;
	//if(fabs(tmean)<1)cout<<"Saving: "<<ktrack<<" Time("<<tmean<<") Index(";
//	cout<<"Saving: "<<ktrack<<" Time("<<tmean<<") Index(";
//	for(int u=0;u<5;u++)cout<<TRKindx[ktrack][u]<<",";
//	cout<<")"<<endl;
    //}
    
  }//3planes involved
} 
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBMSrecons:: track_resc(void){
  //
  int badDec=0;
  int indx[4],kfl[4],i,j,k,l,ipl,ll,ilb,ilb2;
  int ktrack, ntot;
  double delt1, delt2, delt3, delt4, z[4], am;
  double p, chiq;
  double tmean, tchiq;
  double t1,t2,t3,t4;
  //bool SkipThisComb=false;

  //
  //combine the remaining hits to 2-hit sets
  //which can be used by rescue algo:
  //
  for(i=0; i<khit16[0]; i++) {
    for(j=0; j<khit16[1]; j++){
      for(k=0; k<khit16[2]; k++){
	for(l=0; l<khit16[3]; l++) {
	  //check whether any of these hits has been used before
	  kfl[0]=iflc16[i][0];
	  kfl[1]=iflc16[j][1];
	  kfl[2]=iflc16[k][2];
	  kfl[3]=iflc16[l][3];
	  //
	  t1=tbmz16[i][0];
	  t2=tbmz16[j][1];
	  t3=tbmz16[k][2];
	  t4=tbmz16[l][3];
	  //
	  ntot=kfl[0]+kfl[1]+kfl[2]+kfl[3];
	  if(ntot<2) continue;               // select combination with >= 2 hits   
	  //
	  tmean=(tbmz16[i][0]*kfl[0]+tbmz16[j][1]*kfl[1]+
			      tbmz16[k][2]*kfl[2]+tbmz16[l][3]*kfl[3])/ntot;
	  delt1=fabs(tmean-tbmz16[i][0]);
	  delt2=fabs(tmean-tbmz16[j][1]);
	  delt3=fabs(tmean-tbmz16[k][2]);
	  delt4=fabs(tmean-tbmz16[l][3]);
	  //
	  if(ntot==nhod_old){                    
	    //
	    // These are 4 hit tracks - must have been found before and failed
	    // hence remove the two worst hits
	    //
	    am=delt1; ilb=0; ilb2=0;                
	    if(delt2>am){am=delt2; ilb=1;}  
	    if(delt3>am){am=delt3; ilb=2; ilb2=1;}
	    if(delt4>am){ilb=3; ilb2=2;}
	    kfl[ilb]=kfl[ilb2]=0;
	    ntot=2;
	    tmean=(tbmz16[i][0]*kfl[0]+tbmz16[j][1]*kfl[1]+
		   tbmz16[k][2]*kfl[2]+tbmz16[l][3]*kfl[3])/ntot;
	    delt1=fabs(tmean-tbmz16[i][0]);
	    delt2=fabs(tmean-tbmz16[j][1]);
	    delt3=fabs(tmean-tbmz16[k][2]);
	    delt4=fabs(tmean-tbmz16[l][3]);
	  }
	  if(ntot==3){                    
	    //
	    // These are 3 hit tracks - must have been found before and failed
	    // hence remove the worst hit
	    //
	    am=delt1; ilb=0;                
	    if(delt2>am){am=delt2; ilb=1;}  
	    if(delt3>am){am=delt3; ilb=2;}
	    if(delt4>am){ilb=3;}
	    kfl[ilb]=0;
	    ntot=2;
	    tmean=(tbmz16[i][0]*kfl[0]+tbmz16[j][1]*kfl[1]+
		   tbmz16[k][2]*kfl[2]+tbmz16[l][3]*kfl[3])/ntot;
	    delt1=fabs(tmean-tbmz16[i][0]);
	    delt2=fabs(tmean-tbmz16[j][1]);
	    delt3=fabs(tmean-tbmz16[k][2]);
	    delt4=fabs(tmean-tbmz16[l][3]);
	  }
	  //for method3 choose hits by cutting on time difference
	  // between them
	  if(method3){
	    if((kfl[0]!=0 && kfl[1]!=0) && fabs(t1-t2)>BMSpar[1]) continue;
	    if((kfl[0]!=0 && kfl[2]!=0) && fabs(t1-t3)>BMSpar[1]) continue;
	    if((kfl[0]!=0 && kfl[3]!=0) && fabs(t1-t4)>BMSpar[1]) continue;
	    if((kfl[1]!=0 && kfl[2]!=0) && fabs(t2-t3)>BMSpar[1]) continue;
	    if((kfl[1]!=0 && kfl[3]!=0) && fabs(t2-t4)>BMSpar[1]) continue;
	    if((kfl[2]!=0 && kfl[3]!=0) && fabs(t3-t4)>BMSpar[1]) continue;
	  }
	  else{
	    if(delt1>BMSpar[1]) kfl[0]=0;         // check timing
	    if(delt2>BMSpar[1]) kfl[1]=0;
	    if(delt3>BMSpar[1]) kfl[2]=0;
	    if(delt4>BMSpar[1]) kfl[3]=0;
	  }
	  ntot=kfl[0]+kfl[1]+kfl[2]+kfl[3];
	  if(ntot<2) continue;
	  if(kfl[0]==0 && kfl[1]==0) continue;
	  //
	  tchiq=(delt1*delt1)*kfl[0]+
	    (delt2*delt2)*kfl[1]+
	    (delt3*delt3)*kfl[2]+
	    (delt4*delt4)*kfl[3];                       
	  if(tchiq>BMSpar[4]) continue;        // time chiq cut
	  //
	  ktrack=ntrack_resc;
	  indx[0]=(i+1)*kfl[0]-1; 
	  indx[1]=(j+1)*kfl[1]-1;
	  indx[2]=(k+1)*kfl[2]-1; 
	  indx[3]=(l+1)*kfl[3]-1;
	  if(ktrack!=0) {
	    //
	    // No common points allowed for 2 hit tracks -
	    // if find - take track with better time chiq
	    //			       
	    for(ll=0;ll<ntrack_resc;ll++) {
	      for(ntot=0,ipl=0;ipl<nhod_old;ipl++) {
		if(kfl[ipl]&&(indx[ipl]==TRKindx_resc[ll][ipl])) ntot++;
	      }
	      if(ntot==3) goto next2;  
	      if(ntot) {
		if(tchiq>TMchiq_resc[ll]){
		  if(fabs(tchiq-TMchiq_resc[ll])<10) badDec++;
		  goto next2;
		}
		else{
		  if(fabs(tchiq-TMchiq_resc[ll])<10) badDec++;
		}
		ktrack=ll;
		goto save2;
	      }
	    }
	  }
	  if(ntrack_resc==NTRACKMX) return ntrack_resc;
	  ntrack_resc++;
	save2: 
	  bool partPlanes[4];
	  for(ipl=0;ipl<nhod_old;ipl++) { 
	    if(kfl[ipl]) {
	      TRKindx_resc[ktrack][ipl]=indx[ipl]; 
	      z[ipl]=zhit16[indx[ipl]][ipl];
	      TRKy_resc[ktrack][ipl] = z[ipl];
	      partPlanes[ipl]=true;
	    }
	    else  {
	      z[ipl]=1000000.;
	      TRKindx_resc[ktrack][ipl]=9999; 
	      TRKy_resc[ktrack][ipl] = z[ipl];
	      partPlanes[ipl]=false;
	    }
	  } 
          
          //the same for new planes
          for(ipl = 4; ipl < NBMSHODS; ipl++){
	      TRKindx_resc[ktrack][ipl]=9999; 
	      TRKy_resc[ktrack][ipl] = 1000000.;
          }
          
	  TMchiq_resc[ktrack] = tchiq;
	  TRKtime_resc[ktrack] = tmean;  
	  TRKp_resc[ktrack] = 0; 
	  SPchiq_resc[ktrack] = 999;
	  TRKhit_resc[ktrack] = 2;
	  if(partPlanes[0] && partPlanes[1])TRK_resol[ktrack]=BmResol_1_2;
	  if(partPlanes[0] && partPlanes[2])TRK_resol[ktrack]=BmResol_1_3;
	  if(partPlanes[0] && partPlanes[3])TRK_resol[ktrack]=BmResol_1_4;
	  if(partPlanes[1] && partPlanes[2])TRK_resol[ktrack]=BmResol_2_3;
	  if(partPlanes[1] && partPlanes[3])TRK_resol[ktrack]=BmResol_2_4;
	next2:  continue;  
	}
      }
    }
  } //outermost for-loop for 2-hit tracks
  //Fill some histos for these rescue-tracks:
  if(Resc1D[0]!=NULL)Resc1D[0]->Fill(ntrack_resc);
  for(int i=0;i<ntrack_resc;i++){
    if(Resc1D[1]!=NULL)Resc1D[1]->Fill(TMchiq_resc[i]);
    if(Resc1D[2]!=NULL)Resc1D[2]->Fill(TRKtime_resc[i]);
    
  }
  //for(int I=0;I<ntrack_resc;I++){
  //if(fabs(TRKtime_resc[I])<10){
  //cout<<"BMS track: ID("<<I<<") TChi("<<TMchiq_resc[I]<<") Time("<<TRKtime_resc[I]<<") Index(";
  //for(int P=0;P<4;P++)cout<<TRKindx_resc[I][P]<<",";
  //cout<<")"<<endl;
  //}
  //}
  return ntrack_resc;
//
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: moment_resc(int TrackId, double slopeTg, double* p, double* vchisq){
  double z[4];
  for(int ipl=0; ipl<nhod_old; ipl++) z[ipl]=TRKy_resc[TrackId][ipl];
  moment_resc(z,slopeTg,p,vchisq);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBMSrecons:: moment_resc(double z[], double slopeTg, double* p, double* vchisq){
  //*** BEAM MOMENTUM ANALYSIS for 2 hit tracks, using also the slope from 
  //the SF/SI telescope (slopeTg)
  //*------ Input variables:
  //       z(1-4) :  hits from bms hodos
  //                 2 of them are "real hits"
  //                 the other 2 are > 1.0E05  (indicating missing planes)
  //        =====> it is mandatory that at least one of the first 2 hits be REAL
  //       slopeTg  : slope of beam track in vertical plane, in front of target
  //                (as obtained from SciFi)
  //         =====> UNITS are mrad  <===========================
  //------- Output variables:
  //        p:    beam momentum (in GeV) calculated from 2 hits + slope
  //        vchisq: Chi2?
  double z2[4];
  int ihc = 0;
  int nout = 10;
  int imis2 = 0;
  for(int i=0;i<4;i++){
    if (z[i]<1.0E+05){
      z2[ihc] = z[i];
      ihc = ihc+1;
    }
    else{
      if(nout==10)nout=i;
      else{
        z2[ihc] = 0.;
        imis2 = ihc;
        ihc=ihc + 1;
      }
    }
  }
  //-----check one of the missing planes is downstream B6.....
  if (imis2<1){
    cout<<"CsBMSrecons: rescue warning----no plane in front of B6!!!"<<endl;
    //cout<<"z2:["<<z2[0]<<","<<z2[1]<<","<<z2[2]<<","<<z2[3]<<"], nout: "<<nout<<" ihc: "<<ihc<<"imis2: <<imis2<<endl;
    cout<<"P-beam set to nominal value:"<<p0<<endl;
    *p = p0;
    return 0;
  }
  //
  //
  int ncase=1+nout;
  slopeTg=slopeTg*1000; //conversion to mrad
  double zres = (
	  slopeTg*0.6
	  - coeff_coff16[0][0][0][1][ncase]
	  - coeff_coff16[0][0][1][1][ncase] * z2[0]
	  - coeff_coff16[0][0][2][1][ncase] * z2[1]
	  - coeff_coff16[0][0][3][1][ncase] * z2[2]
	  ) / coeff_coff16[0][0][imis2+1][1][ncase];
  z2[imis2] = zres;
  double pin=coeff_coff16[0][0][0][2][ncase];
  for(int i=0;i<3;i++) pin=pin+z2[i]*coeff_coff16[0][0][i+1][2][ncase];
  *p=1.0/pin;
  //cout<<"Rescue: momentum: "<<*p<<endl;
  return 1;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: moment3(double z[],double* p, double* vchisq)
//
// *** BEAM MOMENTUM ANALYSIS for 3 hit tracks
// *** USING PRINCIPLE COMPONENTS AND
// *** MULTI-DIMENSIONAL LEAST SQUARES FITTING.
//
//   input:  Z(5)   - hodoscope hits , if z[i]>10000 -> plane is not used for this track
//   output: p      - momentum
//           vchisq - chi squared
//   
{
	int nhit=0;
	int i, j, ncase=0;
	int missing_old=0;
	int first_missing=-1;
        double z2[nhod], new_z2[3], pb=0;
//	
	for(i=0; i<nhod; i++) {                //  HOW MANY HITS WE HAVE  ?
	    if(z[i]>10000.) {
	      if(i<4){ncase=i;missing_old++;}	 
	      if(first_missing==-1)first_missing=i;
	    }
	    else    z2[nhit++]=z[i];
	}
//
        if(nhit<3 || nhit>3) {                       // too few hits, return
	        *p=1.0e6;  *vchisq=1.0e6;
       	        return; 
        }
//
         *vchisq=0.0;  ncase++; 
//
	 if(missing_old<2){//only old stations used
	   pb = coeff_coff16[0][0][0][2][ncase];  
	   for(i=0; i<3; i++) {
	     pb = pb + z2[i]*coeff_coff16[0][0][i+1][2][ncase];
	   }
	   //correct momentum: [a la Roland]
	   pb=1./pb;
	   if(z[0]>10000.){//using planes 2,3,4
	     pb=pb - 0.0028*(pb-155)*(pb-171);
	     if(pb<150)pb=pb+0.1;
	     if(pb>172)pb=pb-0.1;
	   }
	   if(z[1]>10000.){//using planes 1,3,4
	     if(pb<153)pb=pb-0.12;
	     else pb=pb - 0.0018*(pb-156)*(pb-171);
	   }
	   if(z[2]>10000.){//using planes 1,2,4
	     if(pb<150)pb=pb-0.1;
	     else{
	       if(pb<171)pb=pb + (pb-158)*0.0071;
	       else pb=pb - (pb-173)*0.1;
	     }
	   }
	   if(z[3]>10000.){//using planes 1,2,3
	     pb=pb - 0.0014*(pb-155)*(pb-172);
	     if(pb<150)pb=pb-0.15;
	     if(pb>173)pb=pb-0.15;
	   }
	 }//only old stations used
	 else{//station 5 or 6 has to be used:
	   if(z[4]>10000.){//using BM06
	       if(ncase-1==2){//BM03 is missing, BM06 will be third one
		   //cout<<"Computing momentum for 0 case..."<<endl;
		   if(first_missing==(ncase-1)){
		       cout<<"Problem with BM06-procedure: only one old plane missing!"<<endl;
		       exit(1);
		   }
		   //reshuffle z2: now: z2[0]: (BM01 or BM02), z2[1]: BM04, z2[2]: BM06, BM06 should be at z2[1]
		   new_z2[0]=z2[0];
		   new_z2[1]=z2[2];
		   new_z2[2]=z2[1];
		   pb = coeff_coff16[0][2][0][2][first_missing+1];
		   for(i=0; i<3; i++) {
		       pb = pb + new_z2[i]*coeff_coff16[0][2][i+1][2][first_missing+1];
		   }
	       }
	       if(ncase-1==3){//BM04 is missing, BM06 will be the last one
		   //cout<<"Computing momentum for 1 case..."<<endl;
		   if(first_missing==(ncase-1)){
		       cout<<"Problem with BM06-procedure: only one old plane missing!"<<endl;
		       exit(1);
		   }

		   pb = coeff_coff16[0][1][0][2][first_missing+1];
		   for(i=0; i<3; i++) {
		       pb = pb + z2[i]*coeff_coff16[0][1][i+1][2][first_missing+1];
		   }
	       }//BM04 missing
	       pb=1./pb;
	   }
	   else{//using BM05
	       if(first_missing==0){//BM01 is missing, BM05 will be first one
		   //cout<<"Computing momentum for 0 case..."<<endl;
		   if(first_missing==(ncase-1)){
		       cout<<"Problem with BM05-procedure: only one old plane missing!"<<endl;
		       exit(1);
		   }
		   //reshuffle z2: now: z2[0]: (BM02 or BM03), z2[1]: (BM03 or BM04), z2[2]: BM05, BM05 should be at z2[0]
		   new_z2[0]=z2[2];
		   new_z2[1]=z2[0];
		   new_z2[2]=z2[1];
		   pb = coeff_coff16[1][0][0][2][ncase];
		   for(i=0; i<3; i++) {
		       pb = pb + new_z2[i]*coeff_coff16[1][0][i+1][2][ncase];
		   }
	       }//BM01 missing
	       if(first_missing==1){//BM02 is missing, BM05 will be second one
		   //cout<<"Computing momentum for 1 case..."<<endl;
		   if(first_missing==(ncase-1)){
		       cout<<"Problem with BM05-procedure: only one old plane missing!"<<endl;
		       exit(1);
		   }
		   //reshuffle z2: now: z2[0]: BM01, z2[1]: (BM03 or Bm04), z2[2]: BM05, BM05 should be at z2[1]
		   new_z2[0]=z2[0];
		   new_z2[1]=z2[2];
		   new_z2[2]=z2[1];
		   pb = coeff_coff16[2][0][0][2][ncase];
		   for(i=0; i<3; i++) {
		       pb = pb + new_z2[i]*coeff_coff16[2][0][i+1][2][ncase];
		   }
	       }//BM02 missing
	       pb=1./pb;
	   }
	 }
	 *p=pb;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: moment4(double z[], double* p, double* vchisq)
//
// *** BEAM MOMENTUM ANALYSIS for 4 hit tracks
// *** USING PRINCIPLE COMPONENTS AND
// *** MULTI-DIMENSIONAL LEAST SQUARES FITTING.
//
//   input:  Z(4)   - hodoscope hits
//   output: p      - momentum
//           vchisq - chi squared
//   
{
	int i, j;
	double z1[4];
        double z2[4], pb=0;
//
        for(i=0; i<4; i++) z1[i]=z[i]-coeff_avsg16[0][0][i];
        for(i=0; i<4; i++) {      
            for(j=0, z2[i]=0.; j<4; j++) 
                 z2[i]=z2[i]+z1[j]*coeff_evec16[0][0][j][i];
        }      
        *vchisq=z2[3]*z2[3];
//
        pb = coeff_coff16[0][0][0][2][0];  
        for(i=0; i<3; i++) {      
	    pb = pb + z2[i]*coeff_coff16[0][0][i+1][2][0];
        }
	//correct momentum: [a la Roland]
	pb=1./pb;
	if(pb>150 && pb<170){
	  pb=pb+(pb-160)*0.01;
	}
	else if(pb<150){
	  pb=pb-(pb-146)*0.025;
	}
	else{
	  pb=pb-(pb-173)*0.055;
	}
	//
        //*p=1./pb;
	*p=pb;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBMSrecons:: decode(void)
//
//   get, decode, calibrate and select data         
//   only old planes
//
{
  int ipl, i, nhit, nchnl;
  int rawhit[NBMSHODS];
  int RawHitCnt;
  double time, tm;
  int label = 1;
  Fired[0]=Fired[1]=Fired[2]=Fired[3]=Fired[4]=false;
//
  nfired=0;
  for(ipl=0; ipl<nhod; ipl++) { 
     khit16[ipl]=0; rawhit[ipl]=0; 
     //
     if(!Id[ipl]) continue;
     list<CsDigit*>digit = Id[ipl]->getMyDigits();
     if(!digit.size()) continue;  
     list<CsDigit*>::iterator Idig;
     for( Idig = digit.begin(); Idig != digit.end(); Idig ++ ) {
        CsDigit* digWB = *Idig;
//
	   nchnl = digWB->getAddress();
    	   time  = digWB->getDatum();
//       
          if((nchnl>hodsize[ipl])||(nchnl<0)) {
	         ostringstream message;
                 message << " CsBMScalib:: Wrong channel number:"<< nchnl
                         <<"   hodoscope:"<< ipl << endl;
                 CsErrLog::Instance()->mes(elInfo, message.str());                  
                 continue;
          }
//          if((time>65535.)||(time<0.)) {
//                 ostrstream message;
//                 message << " CsBMScalib:: Wrong time:"<< time<<" hodoscope:"<<ipl << endl;
//                 CsErrLog::Instance()->mes(elInfo, message.str());                  
//                 continue;
//          }
	  //tm=time-10*TMt0[nchnl][ipl]-BMSMN1[itrg16];
	                   //not needed since calib done before making digits
	  //tm=time;
	  time=time*TDCSLP;  //digits will come in ns!
	  //tm is allready in ns!
//
//        some hists for raw data
// 
             rawhit[ipl]++;        
             if(HiFired[ipl]!=NULL) HiFired[ipl]->Fill((double)nchnl);
             //if(HiProf[ipl]!=NULL) HiProf[ipl]->Fill(BMStable[nchnl][ipl]);	     
             if(HiHtime[ipl]!=NULL) HiHtime[ipl]->Fill(time);

	     if(BeamH2D[9]!=NULL) BeamH2D[9]->Fill((double)nchnl,ipl);
	     if(BeamH2D[10]!=NULL) BeamH2D[10]->Fill(BMStable[nchnl][ipl],ipl);
	     if(BeamH2D[11]!=NULL) BeamH2D[11]->Fill(time,ipl);
//
//   time cuts   
//
             if(fabs(time)<BMSpar[0]) {                       //   time window cut
               if(khit16[ipl]>=mxnmbhit) {label=0; continue;}
	       //cout<<"khit16= "<<khit16<<endl;
               nhit=khit16[ipl];
               tbmz16[nhit][ipl]=time;                 //save good hit
               ibmz16[nhit][ipl]=nchnl;
               zhit16[nhit][ipl]=-BMStable[nchnl][ipl];    
               iflc16[nhit][ipl]=1;
               khit16[ipl]++;
	       Fired[ipl]=true;
	       //if(HiHtimeTCut[ipl]!=NULL) HiHtimeTCut[ipl]->Fill(time);
	       if(HiCleanFired[ipl] != NULL) HiCleanFired[ipl]->Fill(nchnl);
	       if(HiProf[ipl]!=NULL) HiProf[ipl]->Fill(-BMStable[nchnl][ipl]);
	       if(BeamH2D[12]!=NULL)BeamH2D[12]->Fill(time,ipl);
             }
     }//loop over digits
        if(ipl<nhod_old)if(khit16[ipl]) nfired++;  
  }
  nhittotal=khit16[0]+khit16[1]+khit16[2]+khit16[3];
  if(nhod==5)nhittotal=khit16[0]+khit16[1]+khit16[2]+khit16[3]+khit16[4];
  for(ipl=0; ipl<nhod; ipl++)  {
    //if(HiNhit[ipl]!=NULL) HiNhit[ipl]->Fill((double)rawhit[ipl]);
    //if(HiNHitafter[ipl]!=NULL) HiNHitafter[ipl]->Fill((double)khit16[ipl]);
    if(BeamH2D[13]!=NULL)BeamH2D[13]->Fill(rawhit[ipl],ipl);
    if(BeamH2D[14]!=NULL)BeamH2D[14]->Fill(khit16[ipl],ipl);
  }
  if(NHitP2vsP1!=NULL) NHitP2vsP1->Fill((double)khit16[0],(double)khit16[1]);
  if(NHitP2vsP3!=NULL) NHitP2vsP3->Fill((double)khit16[2],(double)khit16[1]);
  if(NHitP4vsP3!=NULL) NHitP4vsP3->Fill((double)khit16[2],(double)khit16[3]);

  return label;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBMSrecons:: decodeCl(void)
//
//   get, decode, cluster and select data         
//   no clustering for BM05 and BM06!!
//
{
  //cout<<"Entering decodeCl()..."<<endl;
  int ipl, nhit,RawHitCnt,nchnl, nchnl1,nchnl2;
  int chnlNr,chnlNr1;
  int clCnt[NBMSHODS];
  int rawhit[NBMSHODS],clSize;
  double time, time1,time2, tm, clYPos,clYPos1,clTm,YPos;
  int label = 1;
  FiredHit[0]=FiredHit[1]=FiredHit[2]=FiredHit[3]=false;
//
  nfired=0;

  //Fill raw hits in appropriate arrays:
  for(ipl=0; ipl<nhod; ipl++) { 
    khit16[ipl]=0;
    rawhit[ipl]=0; 
    clCnt[ipl]=0;
    Fired[ipl]=false;

    //old planes
    if(ipl<nhod_old){
      //reset vector which contains cluster-flags (0 if hit has not been used in a cluster)
      for(int i=0; i<BMSBUFSIZE; i++){
	HitUsedCl[i][ipl]=0;
      }
      //reset number of acc hits on this plane, number of raw hits on this plane
      nRawH[ipl]=0;
      if(!Id[ipl]) continue;
      list<CsDigit*>digit = Id[ipl]->getMyDigits();
      if(!digit.size()) continue;  
      list<CsDigit*>::iterator Idig;
      for( Idig = digit.begin(); Idig != digit.end(); Idig ++ ) {
	CsDigit* digWB = *Idig;
	//
	nchnl1 = digWB->getAddress();
	time1  = digWB->getDatum();
	if((nchnl1>hodsize[ipl])||(nchnl1<0)) {
	  ostringstream message;
	  message << " CsBMScalib:: Wrong channel number:"<< nchnl1
		  <<"   hodoscope:"<< ipl << endl;
	  CsErrLog::Instance()->mes(elInfo, message.str());                  
	  continue;
	}
	//
	//        some hists for raw data
	// 
	//nRawH[ipl]++;    
	rawhit[ipl]++; 
        if(HiFired[ipl]!=NULL) HiFired[ipl]->Fill((double)nchnl1);
        if(HiHtime[ipl]!=NULL) HiHtime[ipl]->Fill((double)time1);
	
	if(BeamH2D[9]!=NULL) BeamH2D[9]->Fill((double)nchnl1,ipl);
	if(BeamH2D[10]!=NULL) BeamH2D[10]->Fill(BMStable[nchnl1][ipl],ipl);
	if(BeamH2D[11]!=NULL) BeamH2D[11]->Fill(time1,ipl);
	if(nRawH[ipl]>=mxnmbhit) {label=0; continue;}
	RawHitCnt=nRawH[ipl];
	rawHtime[RawHitCnt][ipl]=time1;
	rawHchnl[RawHitCnt][ipl]=nchnl1;
	nRawH[ipl]++;
	FiredHit[ipl]=true;
      }//loop over digits
      //}//loop over planes
      
      //loop over the arrays:
      for(int hitcnt1=0; hitcnt1<nRawH[ipl]; hitcnt1++){
	//check whether this hit has already been used in a cluster:
	if(HitUsedCl[hitcnt1][ipl]!=0)continue;
	clSize=1;
	chnlNr=rawHchnl[hitcnt1][ipl];
	chnlNr1=chnlNr;
	clYPos=-BMStable[chnlNr][ipl]; 
	clYPos1=clYPos;
	clTm=rawHtime[hitcnt1][ipl];
	HitUsedCl[hitcnt1][ipl]=1;
	for(int hitcnt2=0; hitcnt2<nRawH[ipl]; hitcnt2++){
	  chnlNr=rawHchnl[hitcnt2][ipl];
	  //if(fabs(rawHchnl[hitcnt1][ipl]-rawHchnl[hitcnt2][ipl])<=Cl_ChannelLimit){
	  if(fabs(clYPos1+BMStable[chnlNr][ipl])<=Cl_ChannelLimit){
	    if(HitUsedCl[hitcnt2][ipl]!=0)continue;
	    if(fabs(rawHtime[hitcnt1][ipl]-rawHtime[hitcnt2][ipl])<Cl_TimeLimit){
	      //record channel difference
	      if(Cl_chDiff!=NULL)Cl_chDiff->Fill(fabs(float(rawHchnl[hitcnt1][ipl])-rawHchnl[hitcnt2][ipl]));
	      //record time difference
	      if(Cl_tmDiff!=NULL)Cl_tmDiff->Fill(fabs(rawHtime[hitcnt1][ipl]-rawHtime[hitcnt2][ipl]));
	      //increment cluster size
	      clSize++;
	      //add y-position to sum
	      chnlNr=rawHchnl[hitcnt2][ipl];
	      clYPos=clYPos-BMStable[chnlNr][ipl];
	      //add time to sum
	      clTm=clTm+rawHtime[hitcnt2][ipl];
	      //flag hit as used
	      HitUsedCl[hitcnt2][ipl]=1;
	    }// if: time-difference condition
	  }// if: channel-difference condition
	  
	}//2nd loop over hits
	//ok call the result a cluster:
	clCnt[ipl]++;
	//calculate mean time
	tm=clTm/clSize;
	//calculate mean position
	YPos=clYPos/clSize;
	//histogram clustersize
	if(Cl_clSize!=NULL) Cl_clSize->Fill(clSize);
	
	//histogram position => done further down
	//histogram time => done further down
	
	//nchnl1 = digWB1->getAddress();
	nchnl=chnlNr1;
	
	//fill some histos for raw clusters:
	if(Cl_HiProf[ipl]!=NULL) Cl_HiProf[ipl]->Fill(YPos);
	if(Cl_HiHtime[ipl]!=NULL) Cl_HiHtime[ipl]->Fill(tm);
	//
	//   time cuts   
	//
	if(fabs(tm)<BMSpar[0]) {                       //   time window cut
	  if(khit16[ipl]>=mxnmbhit) {label=0; continue;}
	  //cout<<"khit16= "<<khit16<<endl;
	  nhit=khit16[ipl];
	  tbmz16[nhit][ipl]=tm;                 //save good hit
	  ibmz16[nhit][ipl]=nchnl;
	  //zhit16[nhit][ipl]=-BMStable[nchnl][ipl];    
	  zhit16[nhit][ipl]=YPos;
	  iflc16[nhit][ipl]=1; //this flag marks wheter hit/cluster has been used in a track
	  khit16[ipl]++;
	  Fired[ipl]=true;
	  if(HiCleanFired[ipl] != NULL) HiCleanFired[ipl]->Fill(nchnl);
	  if(HiProf[ipl]!=NULL) HiProf[ipl]->Fill(YPos);
	}
	
      }// 1st loop over hits
      if(khit16[ipl]) nfired++;  
    }//loop over old planes
    //new planes
    if(ipl >= nhod_old){
      if(!Id[ipl]) continue;
      list<CsDigit*>digit = Id[ipl]->getMyDigits();
      if(!digit.size()) continue;  
      list<CsDigit*>::iterator Idig;
      for( Idig = digit.begin(); Idig != digit.end(); Idig ++ ) {
        CsDigit* digWB = *Idig;
	//
	nchnl = digWB->getAddress();
	time  = digWB->getDatum();
	//       
	if((nchnl>hodsize[ipl])||(nchnl<0)) {
	  ostringstream message;
	  message << " CsBMScalib:: Wrong channel number:"<< nchnl
		  <<"   hodoscope:"<< ipl << endl;
	  CsErrLog::Instance()->mes(elInfo, message.str());
	  continue;
	}
	time=time*TDCSLP;  //digits will come in ns!
	//tm is allready in ns!
	//
	//        some hists for raw data
	// 
	rawhit[ipl]++;

        if(HiFired[ipl]!=NULL) HiFired[ipl]->Fill((double)nchnl);
        if(HiHtime[ipl]!=NULL) HiHtime[ipl]->Fill((double)time);
	//
	//   time cuts   
	//
	if(fabs(time)<BMSpar[0]) {                       //   time window cut
	  if(khit16[ipl]>=mxnmbhit) {label=0; continue;}
	  nhit=khit16[ipl];
	  tbmz16[nhit][ipl]=time;                 //save good hit
	  ibmz16[nhit][ipl]=nchnl;
	  zhit16[nhit][ipl]=-BMStable[nchnl][ipl];
	  iflc16[nhit][ipl]=1;
	  khit16[ipl]++;
	  Fired[ipl]=true;
	  if(HiProf[ipl]!=NULL) HiProf[ipl]->Fill(-BMStable[nchnl][ipl]);
	  if(HiCleanFired[ipl] != NULL) HiCleanFired[ipl]->Fill(nchnl);
	}
     }//loop over digits
    }//BM05 && BM06
  } // loop over planes
  nhittotal=khit16[0]+khit16[1]+khit16[2]+khit16[3];
  for(ipl=0; ipl<nhod; ipl++)  {
    if(HiNhit[ipl]!=NULL) HiNhit[ipl]->Fill((double)rawhit[ipl]);
    if(HiNHitafter[ipl]!=NULL) HiNHitafter[ipl]->Fill((double)khit16[ipl]);
  }
  for(ipl=0; ipl<nhod_old; ipl++)  {
    if(Cl_HiNhit[ipl]!=NULL) Cl_HiNhit[ipl]->Fill(clCnt[ipl]);
  }
  if(NHitP2vsP1!=NULL) NHitP2vsP1->Fill((double)khit16[0],(double)khit16[1]);
  if(NHitP2vsP3!=NULL) NHitP2vsP3->Fill((double)khit16[2],(double)khit16[1]);
  if(NHitP4vsP3!=NULL) NHitP4vsP3->Fill((double)khit16[2],(double)khit16[3]);
  
  //cout<<"Leaving decodeCl()..."<<endl;
  return label;
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: checkClustering(void)
{
  //does not consider BM05 & BM06!!
  bool cluster=false;
  double ch1,ch2,t1,t2,deltaCh,deltaT;
  int i,j;
  //plane 1:
  for(i=0; i<khit16[0]; i++) {
    for(j=i+1; j<khit16[0]; j++) {
      //ch1=ibmz16[i][0];
      //ch2=ibmz16[j][0];
      ch1=zhit16[i][0];
      ch2=zhit16[j][0];
      t1=tbmz16[i][0];
      t2=tbmz16[j][0];
      deltaCh=fabs(ch1-ch2);
      deltaT=t1-t2;
      //if((deltaCh==1||deltaChg==2)&&(fabs(deltaT)<4)) cluster=true;
      if(deltaCh<6 && fabs(deltaT)<4) cluster=true; 
      if(P1Clus!=NULL)P1Clus->Fill(deltaCh,deltaT);
    }
  }
  //plane 2:
  for(i=0; i<khit16[1]; i++) {
    for(j=i+1; j<khit16[1]; j++) {
      //ch1=ibmz16[i][1];
      //ch2=ibmz16[j][1];
      ch1=zhit16[i][1];
      ch2=zhit16[j][1];
      t1=tbmz16[i][1];
      t2=tbmz16[j][1];
      deltaCh=fabs(ch1-ch2);
      deltaT=t1-t2;
      //if((deltaCh==1)&&(fabs(deltaT)<4)) cluster=true;
      if(deltaCh<6 && fabs(deltaT)<4) cluster=true; 
      if(P2Clus!=NULL) P2Clus->Fill(deltaCh,deltaT);
    }
  }
  //plane 3:
  for(i=0; i<khit16[2]; i++) {
    for(j=i+1; j<khit16[2]; j++) {
      //ch1=ibmz16[i][2];
      //ch2=ibmz16[j][2];
      ch1=zhit16[i][2];
      ch2=zhit16[j][2];
      t1=tbmz16[i][2];
      t2=tbmz16[j][2];
      deltaCh=fabs(ch1-ch2);
      deltaT=t1-t2;
      //if((deltaCh==1)&&(fabs(deltaT)<4)) cluster=true;
      if(deltaCh<6 && fabs(deltaT)<4) cluster=true; 
      if(P3Clus!=NULL) P3Clus->Fill(deltaCh,deltaT);
    }
  }
  //plane 4:
  for(i=0; i<khit16[3]; i++) {
    for(j=i+1; j<khit16[3]; j++) {
      //ch1=ibmz16[i][3];
      //ch2=ibmz16[j][3];
      ch1=zhit16[i][3];
      ch2=zhit16[j][3];
      t1=tbmz16[i][3];
      t2=tbmz16[j][3];
      deltaCh=fabs(ch1-ch2);
      deltaT=t1-t2;
      //if((deltaCh==1||deltaChg==2)&&(fabs(deltaT)<4)) cluster=true;
      if(deltaCh<6 && fabs(deltaT)<4) cluster=true; 
      if(P4Clus!=NULL) P4Clus->Fill(deltaCh,deltaT);
      
    }
  }
  if(NeventsCl!=NULL){
    if(cluster)NeventsCl->Fill(1);
    else NeventsCl->Fill(0);
  }
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: CommonHits(vector<int> &vec){
  //is supposed to be called by CsBeamRecons::bmreconsTr in those cases, 
  //where ncut>1,nBMSpassed>1...
  //and is to determine whether these tracks share common hits!
  int add_this=0;
  vector<int> invTracks;
  for(unsigned int i=0;i<vec.size();i++)invTracks.push_back(0);
  for(int planeI=0;planeI<nhod;planeI++){
    int mult=1;
    int trackMultInHit=1;
    for(unsigned int TrackI=0;TrackI<vec.size();TrackI++){
      int trackMultInHit=1;
      for(unsigned int TrackJ=TrackI+1;TrackJ<vec.size();TrackJ++){
	if(TRKindx[TrackI][planeI]==TRKindx[TrackJ][planeI]){
	  mult++;
	  trackMultInHit++;
	  invTracks[TrackI]=1;
	  invTracks[TrackJ]=1;
	}//Comb has common hit
      }
      if(trackMultInHit>1 &&CommonHitsHist2!=NULL)CommonHitsHist2->Fill(trackMultInHit);
    }
    if(mult>1 && CommonHitsHist1!=NULL){
      CommonHitsHist1->Fill(mult+add_this);
    }
    add_this+=10;
  }//loop over planes
  //now count how many tracks were involved in common-hits:
  int nInvTracks=0;
  for(unsigned int i=0;i<invTracks.size();i++)if(invTracks[i]!=0)nInvTracks++;
  if(CommonHitsHist3!=NULL)CommonHitsHist3->Fill(vec.size(),nInvTracks);
  if(CommonHitsHist1!=NULL)CommonHitsHist1->Fill(0);
  /*
  //try to see what happens if constraint of 1 allowed common hit is tightened to 0
  //common hits:
  vector<int> accepted_tracks;
  for(unsigned int TrackI=0;TrackI<vec.size();TrackI++){
    //check wheter this track has a hit in common with the rest of the tracks:
    for(unsigned int acc=0;acc<accepted_tracks.size();acc++){
      for(nn=0, ipl=0;ipl<nhod;ipl++){
	if(TRKindx[TrackI][ipl]==TRKindx[acc][ipl]) nn++;
      }
      if(nn>0) {
	if(tchiq>TMchiq[ll]){
	  if(TrackDec!=NULL)TrackDec->Fill(7);
	  if(TchiqDiff!=NULL)TchiqDiff->Fill(fabs(tchiq-TMchiq[ll]));
	  if(fabs(tchiq-TMchiq[ll])<10) badDec++;
	  goto next; 
	}
      }
    }
  }
  */
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: MultBMStracks(vector<int> trackID)
{ //check properties of multiple tracks which take part in BeamCombs
  //plot BMS1y vs BMS2y and same for BMS3/4

  for(unsigned int i=0;i<trackID.size();i++){
    if(TRKy[trackID[i]][0]<99999 && TRKy[trackID[i]][1]<99999) if(BeamH2D[0]!=NULL)BeamH2D[0]->Fill(TRKy[trackID[i]][0],TRKy[trackID[i]][1]);
    if(TRKy[trackID[i]][3]<99999 && TRKy[trackID[i]][2]<99999) if(BeamH2D[1]!=NULL)BeamH2D[1]->Fill(TRKy[trackID[i]][2],TRKy[trackID[i]][3]);
    
    for(unsigned int j=i+1;j<trackID.size();j++){
      //look at individual planes plot dt vs dy:
      if(TRKy[trackID[i]][0]<99999 &&TRKy[trackID[j]][0]<99999 ){
	if(BeamH2D[2]!=NULL)BeamH2D[2]->Fill( fabs(TRKy[trackID[i]][0]-TRKy[trackID[j]][0]),tbmz16[TRKindx[trackID[i]][0]][0]-tbmz16[TRKindx[trackID[j]][0]][0]);
      }
      if(TRKy[trackID[i]][1]<99999 &&TRKy[trackID[j]][1]<99999 ){
	if(BeamH2D[3]!=NULL)BeamH2D[3]->Fill( fabs(TRKy[trackID[i]][1]-TRKy[trackID[j]][1]),tbmz16[TRKindx[trackID[i]][1]][1]-tbmz16[TRKindx[trackID[j]][1]][1]);
      }
      if(TRKy[trackID[i]][2]<99999 &&TRKy[trackID[j]][2]<99999 ){
	if(BeamH2D[4]!=NULL)BeamH2D[4]->Fill( fabs(TRKy[trackID[i]][2]-TRKy[trackID[j]][2]),tbmz16[TRKindx[trackID[i]][2]][2]-tbmz16[TRKindx[trackID[j]][2]][2]);
      }
      if(TRKy[trackID[i]][3]<99999 &&TRKy[trackID[j]][3]<99999 ){
	if(BeamH2D[5]!=NULL)BeamH2D[5]->Fill( fabs(TRKy[trackID[i]][3]-TRKy[trackID[j]][3]),tbmz16[TRKindx[trackID[i]][3]][3]-tbmz16[TRKindx[trackID[j]][3]][3]);
      }
      //look at momentum difference:
      if(BeamH1D[0]!=NULL)BeamH1D[0]->Fill(TRKp[trackID[i]]-TRKp[trackID[j]]);
    }
  }
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: TimeRes(void)
{
  double t1,t2,t3,t4,t5,t6,delt;
  int i,j,k,l;
  //Get time resolution from two dedicated channels:
  //Plot hitmap(BMS01) vs hitmap(BMS02):
  for(i=0; i<khit16[0]; i++) {
    for(j=0; j<khit16[1]; j++){
      if(fabs(tbmz16[i][0]-tbmz16[j][1])<15){
	if(HitBMS01VsHitBMS02!=NULL)HitBMS01VsHitBMS02->Fill(ibmz16[i][0],ibmz16[j][1]); 
      }
      if(ibmz16[i][0]==34 && ibmz16[j][1]==33){
	if(TimeRes1_0102!=NULL) TimeRes1_0102->Fill(tbmz16[i][0]-tbmz16[j][1]);
      }
      if(ibmz16[i][0]==27 && ibmz16[j][1]==25){
	if(TimeRes2_0102!=NULL) TimeRes2_0102->Fill(tbmz16[i][0]-tbmz16[j][1]);
      }
    }//loop over digits arrays
  }
  
  //Do the same for a pair for BMS03/BMS04
  for(i=0; i<khit16[2]; i++) {
    for(j=0; j<khit16[3]; j++){
      if(fabs(tbmz16[i][2]-tbmz16[j][3])<15){
	if(HitBMS03VsHitBMS04!=NULL)HitBMS03VsHitBMS04->Fill(ibmz16[i][2],ibmz16[j][3]); 
      }
      if(ibmz16[i][2]==36 && ibmz16[j][3]==33){
	if(TimeRes1_0304!=NULL) TimeRes1_0304->Fill(tbmz16[i][2]-tbmz16[j][3]);
      }
      if(ibmz16[i][2]==27 && ibmz16[j][3]==30){
	if(TimeRes2_0304!=NULL) TimeRes2_0304->Fill(tbmz16[i][2]-tbmz16[j][3]);
      }
    }//loop over digits arrays
  }

  //And the same for BMS01/BMS04
  for(i=0; i<khit16[0]; i++) {
    for(j=0; j<khit16[3]; j++){
      if(fabs(tbmz16[i][0]-tbmz16[j][3])<15){
	if(HitBMS01VsHitBMS04!=NULL)HitBMS01VsHitBMS04->Fill(ibmz16[i][0],ibmz16[j][3]); 
      }
      if(ibmz16[i][0]==34 && ibmz16[j][3]==33){
	if(TimeRes1_0104!=NULL) TimeRes1_0104->Fill(tbmz16[i][0]-tbmz16[j][3]);
      }
      if(ibmz16[i][0]==27 && ibmz16[j][3]==30){
	if(TimeRes2_0104!=NULL) TimeRes2_0104->Fill(tbmz16[i][0]-tbmz16[j][3]);
      }
    }//loop over digits arrays
  }

  //Take differences between hits on different planes (just once for each pair)
  
  //0<->1
  //0<->2
  //0<->3
  //0<->4 | in extra histo
  //0<->5 | in extra histo
  for(i=0; i<khit16[0]; i++) {
    for(j=0; j<khit16[1]; j++){
      t1=tbmz16[i][0];
      t2=tbmz16[j][1];
      delt=t1-t2;
      if(DeltaTclean!=NULL)DeltaTclean->Fill(delt);
    }
    for(k=0; k<khit16[2]; k++){
      t1=tbmz16[i][0];
      t3=tbmz16[k][2];
      delt=t1-t3;
      if(DeltaTclean!=NULL)DeltaTclean->Fill(delt);
    }
    for(l=0; l<khit16[3]; l++) {
      t1=tbmz16[i][0];
      t4=tbmz16[l][3];
      delt=t1-t4;
      if(DeltaTclean!=NULL)DeltaTclean->Fill(delt);
    }
    if(nhod>nhod_old){
      for(l=0; l<khit16[4]; l++) {
	t1=tbmz16[i][0];
	t5=tbmz16[l][4];
	delt=t1-t5;
	if(DeltaTcleanBMSnew!=NULL)DeltaTcleanBMSnew->Fill(delt);
      }
      if(nhod>5) {
	  for(l=0; l<khit16[5]; l++) {
	      t1=tbmz16[i][0];
	      t6=tbmz16[l][5];
	      delt=t1-t6;
	      if(DeltaTcleanBMSnew!=NULL)DeltaTcleanBMSnew->Fill(delt);
	  }
      }
    }
  }  
  //1<->2
  //1<->3
  //1<->4 | in extra histo
  //1<->5 | in extra histo
  for(j=0; j<khit16[1]; j++){
    for(k=0; k<khit16[2]; k++){
      t2=tbmz16[j][1];
      t3=tbmz16[k][2];
      delt=t2-t3;
      if(DeltaTclean!=NULL)DeltaTclean->Fill(delt);
    }
    for(l=0; l<khit16[3]; l++) {
      t2=tbmz16[j][1];
      t4=tbmz16[l][3];
      delt=t2-t4;
      if(DeltaTclean!=NULL)DeltaTclean->Fill(delt);
    }
    if(nhod>nhod_old){
      for(l=0; l<khit16[4]; l++) {
	t2=tbmz16[i][1];
	t5=tbmz16[l][4];
	delt=t2-t5;
	if(DeltaTcleanBMSnew!=NULL)DeltaTcleanBMSnew->Fill(delt);
      }
      if(nhod>5) {
	  for(l=0; l<khit16[5]; l++) {
	      t2=tbmz16[i][1];
	      t6=tbmz16[l][5];
	      delt=t2-t6;
	      if(DeltaTcleanBMSnew!=NULL)DeltaTcleanBMSnew->Fill(delt);
	  }
      }
    }
  }  
  //2<->3
  //2<->4 | in extra histo
  //2<->5 | in extra histo
  for(k=0; k<khit16[2]; k++){
    for(l=0; l<khit16[3]; l++) {
      t3=tbmz16[k][2];
      t4=tbmz16[l][3];
      delt=t3-t4;
      if(DeltaTclean!=NULL)DeltaTclean->Fill(delt);
    }
    if(nhod>nhod_old){
      for(l=0; l<khit16[4]; l++) {
	t3=tbmz16[i][2];
	t5=tbmz16[l][4];
	delt=t3-t5;
	if(DeltaTcleanBMSnew!=NULL)DeltaTcleanBMSnew->Fill(delt);
      }
      if(nhod>5) {
      for(l=0; l<khit16[5]; l++) {
	t3=tbmz16[i][2];
	t6=tbmz16[l][5];
	delt=t3-t6;
	if(DeltaTcleanBMSnew!=NULL)DeltaTcleanBMSnew->Fill(delt);
      }
      }
    }
  }
  //3<->4 | in extra histo
  //3<->5 | in extra histo
  if(nhod>nhod_old){
    for(k=0; k<khit16[3]; k++){
      for(l=0; l<khit16[4]; l++) {
	t4=tbmz16[i][3];
	t5=tbmz16[l][4];
	delt=t4-t5;
	if(DeltaTcleanBMSnew!=NULL)DeltaTcleanBMSnew->Fill(delt);
      }
      if(nhod>5) {
	  for(l=0; l<khit16[5]; l++) {
	      t4=tbmz16[i][3];
	      t6=tbmz16[l][5];
	      delt=t4-t6;
	      if(DeltaTcleanBMSnew!=NULL)DeltaTcleanBMSnew->Fill(delt);
	  }
      }
    }
  }
  //get correlation between thetaY and mom for 4-hit tracks:
  int nhit,n;
  double thetaY=0;
  for(n=0; n<ntrack; n++) { 
    nhit =  TRKhit[n];
    if(nhit==4){
      thetaY=atan((TRKy[n][3]-TRKy[n][2])/12423);
    }
    if(ThetaYVsMom!=NULL)ThetaYVsMom->Fill(TRKp[n],thetaY);
  }
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: Efficiency(int i, double RefTime)
{
  //i==0: "normal" efficiencies
  //i==1: efficiencies for cleaned trigger
  //i==2: efficiencies for special events with a special ref time
  if(i==0){
    if(BMSeff0!=NULL){
      int nHitsSt[nhod];
      for(int i=0;i<nhod;i++)nHitsSt[i]=0;
      //first check which planes have fired (cut relative to refTime)
      for(int planeI=0;planeI<nhod;planeI++){  
	for(int hitI=0;hitI<khit16[planeI];hitI++){
	  if(planeI<nhod_old) {
	      if(tbmz16[hitI][planeI]-RefTime>EffTimeCut_min && tbmz16[hitI][planeI]-RefTime<EffTimeCut_max){
		  nHitsSt[planeI]++;
	      }
	  }
#warning TODO: Check time cut for BM06
	  if(planeI>=4){//2003: Plane5: sigma 0.980225 center: -0.127013 ->[-3.07,2.81]
	    if(tbmz16[hitI][planeI]-RefTime>-3.0 && tbmz16[hitI][planeI]-RefTime<3.0){
	      nHitsSt[planeI]++;
	    }
	  }
	}
      }//loop over planes
      //for all planes:
	if((nHitsSt[0]!=0 &&nHitsSt[1]!=0)&&(nHitsSt[2]!=0 && nHitsSt[3]!=0)){
	  BMSeff0->Fill(0);
	}
	//for plane 0:
	if(nHitsSt[1]&&(nHitsSt[2]&&nHitsSt[3])) BMSeff0->Fill(1); 
	//for plane 1:
	if(nHitsSt[0]&&(nHitsSt[2]&&nHitsSt[3])) BMSeff0->Fill(2);
	//for plane 2:
	if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[3])) BMSeff0->Fill(3);
	//for plane 3:
	if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])) BMSeff0->Fill(4);
	//}
	//
	//including plane 4: (new BM05)
	//
	if(nhod>=5){
	  if((nHitsSt[0]!=0 &&nHitsSt[1]!=0)&&(nHitsSt[2]!=0 && nHitsSt[3]!=0)&& nHitsSt[4]){
	    BMSeff0->Fill(5);
	  }
	  //for plane 0:
	  if(nHitsSt[1]&&(nHitsSt[2]&&nHitsSt[3])&&nHitsSt[4]) BMSeff0->Fill(6); 
	  //for plane 1
	  if(nHitsSt[0]&&(nHitsSt[2]&&nHitsSt[3])&&nHitsSt[4]) BMSeff0->Fill(7);
	  //for plane 2:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[3])&&nHitsSt[4]) BMSeff0->Fill(8);
	  //for plane 3:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[4]) BMSeff0->Fill(9);
	  //for plane 4:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[3]) BMSeff0->Fill(10);
	}
	//
	//including plane 5: (new BM06)
	//
	if(nhod>=6){
	  if((nHitsSt[0]!=0 &&nHitsSt[1]!=0)&&(nHitsSt[2]!=0 && nHitsSt[3]!=0)&& nHitsSt[4] && nHitsSt[5]){
	    BMSeff0->Fill(11);
	  }
	  //for plane 0:
	  if(nHitsSt[1]&&(nHitsSt[2]&&nHitsSt[3])&&nHitsSt[4]&& nHitsSt[5]) BMSeff0->Fill(12); 
	  //for plane 1
	  if(nHitsSt[0]&&(nHitsSt[2]&&nHitsSt[3])&&nHitsSt[4]&& nHitsSt[5]) BMSeff0->Fill(13);
	  //for plane 2:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[3])&&nHitsSt[4]&& nHitsSt[5]) BMSeff0->Fill(14);
	  //for plane 3:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[4]&& nHitsSt[5]) BMSeff0->Fill(15);
	  //for plane 4:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[3]&& nHitsSt[5]) BMSeff0->Fill(16);
	  //for plane 5:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[3]&& nHitsSt[4]) BMSeff0->Fill(17);
	}
    }
  }
  if(i==1){
    if(BMSeff1!=NULL){
      //for all planes:
      if((Fired[0]&&Fired[1])&&(Fired[2]&&Fired[3])){
	BMSeff1->Fill(0);
      }
      //for plane 0:
      if(Fired[1]&&(Fired[2]&&Fired[3])) BMSeff1->Fill(1); 
      //for plane 1:
      if(Fired[0]&&(Fired[2]&&Fired[3])) BMSeff1->Fill(2);
      //for plane 2:
      if(Fired[0]&&(Fired[1]&&Fired[3])) BMSeff1->Fill(3);
      //for plane 3:
      if(Fired[0]&&(Fired[1]&&Fired[2])) BMSeff1->Fill(4);
    }
  }
  if(i==2){    
    if(BMSeff2!=NULL){
      int nHitsSt[nhod];
      for(int i=0;i<nhod;i++)nHitsSt[i]=0;
      //first check which planes have fired (cut relative to refTime)
      for(int planeI=0;planeI<nhod;planeI++){  
	for(int hitI=0;hitI<khit16[planeI];hitI++){
	  if(trig2H[3]!=NULL)trig2H[3]->Fill(tbmz16[hitI][planeI]-RefTime,planeI);
	  if(planeI<4){
	    if(tbmz16[hitI][planeI]-RefTime>EffTimeCut_min && tbmz16[hitI][planeI]-RefTime<EffTimeCut_max){
	      nHitsSt[planeI]++;
	    }
	  }
#warning TODO: Check time cut for BM06
	  if(planeI>=4){//2003: Plane5: sigma 0.980225 center: -0.127013 ->[-3.07,2.81]
	    if(tbmz16[hitI][planeI]-RefTime>-3.0 && tbmz16[hitI][planeI]-RefTime<3.0){
	      nHitsSt[planeI]++;
	    }
	  }
	}
      }//loop over planes
      //determine efficiency:
      //consider only events without multiple hits in a given plane:
      //if((nHitsSt[0]<=1 && nHitsSt[1]<=1)&&(nHitsSt[2]<=1 &&nHitsSt[3]<=1)){ 
      //for all planes:
      if((nHitsSt[0]!=0 &&nHitsSt[1]!=0)&&(nHitsSt[2]&&nHitsSt[3])){
	BMSeff2->Fill(0);
      }
      //for plane 0:
      if(nHitsSt[1]&&(nHitsSt[2]&&nHitsSt[3])) BMSeff2->Fill(1); 
      //for plane 1:
      if(nHitsSt[0]&&(nHitsSt[2]&&nHitsSt[3])) BMSeff2->Fill(2);
      //for plane 2:
      if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[3])) BMSeff2->Fill(3);
      //for plane 3:
      if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])) BMSeff2->Fill(4);
      //}
      //
      //including plane 4: (new BM05)
      //
      if(nhod>=5){
	if((nHitsSt[0]!=0 &&nHitsSt[1]!=0)&&(nHitsSt[2]!=0 && nHitsSt[3]!=0)&& nHitsSt[4]){
	  BMSeff2->Fill(5);
	}
	//for plane 0:
	if(nHitsSt[1]&&(nHitsSt[2]&&nHitsSt[3])&&nHitsSt[4]) BMSeff2->Fill(6); 
	//for plane 1
	if(nHitsSt[0]&&(nHitsSt[2]&&nHitsSt[3])&&nHitsSt[4]) BMSeff2->Fill(7);
	//for plane 2:
	if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[3])&&nHitsSt[4]) BMSeff2->Fill(8);
	//for plane 3:
	if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[4]) BMSeff2->Fill(9);
	//for plane 4:
	if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[3]) {
	  BMSeff2->Fill(10);
	}
      }
	//
	//including plane 5: (new BM06)
	//
	if(nhod>=6){
	  if((nHitsSt[0]!=0 &&nHitsSt[1]!=0)&&(nHitsSt[2]!=0 && nHitsSt[3]!=0)&& nHitsSt[4] && nHitsSt[5]){
	    BMSeff2->Fill(11);
	  }
	  //for plane 0:
	  if(nHitsSt[1]&&(nHitsSt[2]&&nHitsSt[3])&&nHitsSt[4]&& nHitsSt[5]) BMSeff2->Fill(12);
	  //for plane 1
	  if(nHitsSt[0]&&(nHitsSt[2]&&nHitsSt[3])&&nHitsSt[4]&& nHitsSt[5]) BMSeff2->Fill(13);
	  //for plane 2:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[3])&&nHitsSt[4]&& nHitsSt[5]) BMSeff2->Fill(14);
	  //for plane 3:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[4]&& nHitsSt[5]) BMSeff2->Fill(15);
	  //for plane 4:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[3]&& nHitsSt[5]) BMSeff2->Fill(16);
	  //for plane 5:
	  if(nHitsSt[0]&&(nHitsSt[1]&&nHitsSt[2])&&nHitsSt[3]&& nHitsSt[4]) BMSeff2->Fill(17);
	}
    }
  }//Mode 2
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CsBMSrecons:: SFTrigger(double RefTime,int Mode,int trackID)
{
  //called from CsBeamRecons in case of a special "trigger", eg
  // * a certain number of planes in SF03/SF04 have fired
  // * for random-trigger: a certain number of spec-SF have fired
  //RefTime: eg track-time of the SF-track in spectrometer
  //Mode: 0 called for all special events, returns number of planes 
  //        which have fired correlated to the RefTime
  //      1 called for special events with nComb=0 and !BMSok
  //      2 called for ncut==1 (1 beamparticle reconstructed)
  //      3 called for all good events with mom(spec)<momCut(~150GeV)
  //        performes similar things as 1!
  //record efficiency and number of tracks in this event:
  int returnFiredPlanes=0;
  bool StationOk[nhod];
  for(int i=0;i<nhod;i++)StationOk[i]=false;
  //
  if(Mode==0){
    Efficiency(2,RefTime);
    if(NuTrpEvClean!=NULL)NuTrpEvClean->Fill(ntrack);
    if(BMSPlanesfiredNTcl!=NULL && ntrack==0) BMSPlanesfiredNTcl->Fill(nfired);
    for(int planeI=0;planeI<nhod_old;planeI++){  
      for(int hitI=0;hitI<khit16[planeI];hitI++){
	if(tbmz16[hitI][planeI]-RefTime>SpecTrTmdiff_min && tbmz16[hitI][planeI]-RefTime<SpecTrTmdiff_max){
	  StationOk[planeI]=true;
	}
      }
    }
    //record how many planes have fired inside the given time interval:
    int nfiredSpec=0;
    for(int i=0;i<nhod;i++)if(StationOk[i])nfiredSpec++;

    //check planes channels
    if(StationOk[1] && StationOk[2] && StationOk[3])
       for(int i=0; i<khit16[0]; i++)
           if(HiPlanesCheck[0]!=NULL) HiPlanesCheck[0]->Fill((double)ibmz16[i][0]);
    if(StationOk[0] && StationOk[2] && StationOk[3])
       for(int i=0; i<khit16[1]; i++)
           if(HiPlanesCheck[1]!=NULL) HiPlanesCheck[1]->Fill((double)ibmz16[i][1]);
    if(StationOk[0] && StationOk[1] && StationOk[3])
       for(int i=0; i<khit16[2]; i++)
           if(HiPlanesCheck[2]!=NULL) HiPlanesCheck[2]->Fill((double)ibmz16[i][2]);
    if(StationOk[0] && StationOk[1] && StationOk[2])
       for(int i=0; i<khit16[3]; i++)
           if(HiPlanesCheck[3]!=NULL) HiPlanesCheck[3]->Fill((double)ibmz16[i][3]);

    returnFiredPlanes=nfiredSpec;
  }
  //
  //
  if(Mode==1){
    //record number of tracks:
    if(trigH[13]!=NULL)trigH[13]->Fill(ntrack);
    //record time-differences (w respect to RefTime) of all hits of different planes
    //for(int planeI=0;planeI<nhod;planeI++){
    for(int planeI=0;planeI<4;planeI++){  
      for(int hitI=0;hitI<khit16[planeI];hitI++){
	if(trigH[20]!=NULL)trigH[20]->Fill(tbmz16[hitI][planeI]-RefTime);
	if(tbmz16[hitI][planeI]-RefTime>SpecTrTmdiff_min && tbmz16[hitI][planeI]-RefTime<SpecTrTmdiff_max){
	  StationOk[planeI]=true;
	}
      }
    }
    //record how many planes have fired inside the given time interval:
    int nfiredSpec=0;
    for(int i=0;i<nhod;i++)if(StationOk[i])nfiredSpec++;
    if(trigH[18]!=NULL) trigH[18]->Fill(nfiredSpec);
    if(nfiredSpec<4 && trig2H[0]!=NULL){
      for(int i=0;i<nhod;i++){
	if(!StationOk[i])trig2H[0]->Fill(i,nfiredSpec);
      }
    }
    if(nfiredSpec==2 && trig2H[1]!=NULL){
      int notfired[2];
      int notfiredCnt=0;
      for(int i=0;i<nhod;i++){
	if(!StationOk[i]){
	  notfired[notfiredCnt]=i;
	  notfiredCnt++;
	}
      }
      trig2H[1]->Fill(notfired[0]+1,notfired[1]+1);
    }
  }//Mode1
  if(Mode==2){
    double hit_times[4];
    int indx = 0;
    for(int planeI=0;planeI<4;planeI++){  
      for(int hitI=0;hitI<khit16[planeI];hitI++){
	if(trigH[14+planeI]!=NULL)trigH[14+planeI]->Fill(tbmz16[hitI][planeI]-RefTime);
      }
	if(ntrack) {
		indx=TRKindx[trackID][planeI];
	}
	else if(ntrack_resc) {
		indx=TRKindx_resc[trackID][planeI];
	}
      
      if(trigH[19]!=NULL) trigH[19]->Fill(tbmz16[indx][planeI]-RefTime);
      hit_times[planeI]=tbmz16[indx][planeI];
    }
    for(int i=0;i<4;i++){
      for(int j=i+1;j<4;j++){
	if(trigH[21]!=NULL) trigH[21]->Fill(hit_times[i]-hit_times[j]);
      }
    }
  }//Mode 2
  //
  //
  //
  if(Mode==3){
    cout<<"Entering Mode 3"<<endl;
    cout<<"Reftime:"<<endl;
    cout<<RefTime<<endl;
    for(int planeI=0;planeI<4;planeI++){  
      cout<<"accessing plane:"<<planeI<<" "<<khit16[planeI]<< endl;
      for(int hitI=0;hitI<khit16[planeI];hitI++){
	//cout<<"accessing plane "<<planeI<<" hit:"<<hitI<<" ("<<khit16[planeI]<<endl;
	if(tbmz16[hitI][planeI]-RefTime>SpecTrTmdiff_min && tbmz16[hitI][planeI]-RefTime<SpecTrTmdiff_max){
	  //cout<<"Setting StationOk["<<planeI<<"] to true..."<<endl;
	  StationOk[planeI]=true;
	  //cout<<"Done setting StationOk..."<<endl;
	}
      }
    }
    cout<<"stop1"<<endl;
    //record how many planes have fired inside the given time interval:
    int nfiredSpec=0;
    for(int i=0;i<nhod;i++)if(StationOk[i])nfiredSpec++;
    if(trigH[22]!=NULL) trigH[22]->Fill(nfiredSpec);
    if(nfiredSpec==2 && trig2H[2]!=NULL){
      int notfired[2];
      int notfiredCnt=0;
      for(int i=0;i<nhod;i++){
	if(!StationOk[i]){
	  notfired[notfiredCnt]=i;
	  notfiredCnt++;
	}
      }
      trig2H[2]->Fill(notfired[0]+1,notfired[1]+1);
    }
    cout<<"Leaving Mode 3"<<endl;
  }//Mode3
  return returnFiredPlanes;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: TriggerChecks(bool Trigger[],int mode)
{
  if(mode==0){
    //called from CsBeamRecons to study eg the number of tracks per event
    //for the different trigger flavours
    //Only called if there was only one trigger in mask
    int add_this=0;
    if(ntrack==0) add_this=0;
    else add_this=1;
    if(RecPerfVStrig!=NULL){
      if(Trigger[0]) RecPerfVStrig->Fill(add_this); //IT
      if(Trigger[1]) RecPerfVStrig->Fill(10+add_this); //MT
      if(Trigger[2]) RecPerfVStrig->Fill(20+add_this); //LT
      if(Trigger[3]) RecPerfVStrig->Fill(30+add_this); //OT
      if(Trigger[4]) RecPerfVStrig->Fill(40+add_this); //CT
      if(Trigger[5]) RecPerfVStrig->Fill(50+add_this); //IMT
    }
    if(trigH[0]!=NULL)trigH[0]->Fill(ntrack);
    if(DebugMode==1){
      //record number of fired planes for six conditions:
      //oneTrigger,IT,MT,LT,OT,IMT
      int32 bSize;
      bSize=CsEvent::Instance()->getDaqEvent().GetLength();
      //cout<<bSize<<endl;
      int add_this=0;
      if(BMSPlanesfiredTr!=NULL){
	add_this=0;
	//if(((!Trigger[0] && !Trigger[1])&&(!Trigger[2] && !Trigger[3]))&&!Trigger[4]) BMSPlanesfiredTr->Fill(nfired+add_this);
	BMSPlanesfiredTr->Fill(nfired+add_this);
	if(nfired<3) EvSiVStri->Fill(bSize/1000,0);
	if(Trigger[0]) {
	  add_this=10;//IT
	  BMSPlanesfiredTr->Fill(nfired+add_this);
	  if(nfired<3) EvSiVStri->Fill(bSize/1000,1);
	  if(trigH[1]!=NULL)trigH[1]->Fill(ntrack);
	}
	if(Trigger[1]) {
	  add_this=20;//MT
	  BMSPlanesfiredTr->Fill(nfired+add_this);
	  if(nfired<3) EvSiVStri->Fill(bSize/1000,2);
	  if(trigH[2]!=NULL)trigH[2]->Fill(ntrack);
	  
	}
	if(Trigger[2]) {
	  add_this=30;//LT
	  BMSPlanesfiredTr->Fill(nfired+add_this);
	  if(nfired<3) EvSiVStri->Fill(bSize/1000,3);
	  if(trigH[3]!=NULL)trigH[3]->Fill(ntrack);
	}
	if(Trigger[3]) {
	  add_this=40;//OT
	  BMSPlanesfiredTr->Fill(nfired+add_this);
	  if(nfired<3) EvSiVStri->Fill(bSize/1000,4);
	  if(trigH[4]!=NULL)trigH[4]->Fill(ntrack);
	}
	if(Trigger[4]) {
	  add_this=50;//CT
	  BMSPlanesfiredTr->Fill(nfired+add_this);
	  if(nfired<3) EvSiVStri->Fill(bSize/1000,5);
	  if(trigH[5]!=NULL)trigH[5]->Fill(ntrack);
	}
	if(Trigger[5]) {
	  add_this=60;//IMT
	  BMSPlanesfiredTr->Fill(nfired+add_this);
	  if(nfired<3) EvSiVStri->Fill(bSize/1000,6);
	  if(trigH[6]!=NULL)trigH[6]->Fill(ntrack);
	}
      }
    }
  }
  if(mode==1){
    //if(trigH[0]!=NULL)trigH[0]->Fill(ntrack);
    if(Trigger[0]) { //IT
      if(trigH[7]!=NULL)trigH[7]->Fill(ntrack);
    }
    if(Trigger[1]) { //MT
      if(trigH[8]!=NULL)trigH[8]->Fill(ntrack);
    }
    if(Trigger[2]) { //LT
      if(trigH[9]!=NULL)trigH[9]->Fill(ntrack);
    }
    if(Trigger[3]) { //OT
      if(trigH[10]!=NULL)trigH[10]->Fill(ntrack);
    }
    if(Trigger[4]) { //CT
      if(trigH[11]!=NULL)trigH[11]->Fill(ntrack);
    }
    if(Trigger[5]) { //IMT
      if(trigH[12]!=NULL)trigH[12]->Fill(ntrack);
    }

  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: checkBMStbl(void){
  for(int i=0;i<nhod;i++){
    for(int ch=0;ch<64;ch++){
      cout<<"Plane "<<i<<",Channel "<<ch<<": "<<BMStable[ch][i]<<endl;
    
    }
  }

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: BMStbl(void)
//
//    decoding table for BMS hodoscopes
//
{
     int ipl, ll, chn;
         for(ipl=0; ipl < nhod; ipl++) {
             for(ll=1; ll <= hodsize[ipl]; ll++){
                chn = ll-1;
                switch(ipl){
                   case 0: if    (ll<=4)              BMStable[chn][ipl] = -87.5+5.*(ll-1); 
                           else if((ll>4)&&(ll<=60))  BMStable[chn][ipl] = -67.5+5.*((ll-5)/2);
                           else                       BMStable[chn][ipl] =  72.5+5.*(ll-61);
                           break;
                   case 1: if    (ll<=8)              BMStable[chn][ipl] = -42.5+5.*((ll-1)/2); 
                           else if((ll>8 )&&(ll<=20)) BMStable[chn][ipl] = -22.5+5.*((ll-9)/4);
                           else if((ll>20)&&(ll<=44)) BMStable[chn][ipl] =  -7.5+5.*((ll-21)/6);
                           else if((ll>44)&&(ll<=56)) BMStable[chn][ipl] =  12.5+5.*((ll-45)/4);
                           else                       BMStable[chn][ipl] =  27.5+5.*((ll-57)/2);
                           break;
                   case 2: if    (ll<=12)             BMStable[chn][ipl] = -47.5+5.*((ll-1)/2); 
                           else if((ll>12)&&(ll<=20)) BMStable[chn][ipl] = -17.5+5.*((ll-13)/4);
                           else if((ll>20)&&(ll<=44)) BMStable[chn][ipl] =  -7.5+5.*((ll-21)/6);
                           else if((ll>44)&&(ll<=54)) BMStable[chn][ipl] =  12.5+5.*((ll-45)/4);
                           else                       BMStable[chn][ipl] =  27.5+5.*((ll-55)/2);
                           break;
                   case 3: if    (ll<=14)             BMStable[chn][ipl] =-112.5+5.*(ll-1); 
		           else if((ll>14)&&(ll<=50)) BMStable[chn][ipl] = -42.5+5.*((ll-15)/2);
                           else                       BMStable[chn][ipl] =  47.5+5.*(ll-51);
                           break;
		   case 4: BMStable[chn][ipl] = (chn - 27.5) * 2.5 ; //According to Rainer Joosten
                           break;
	           case 5: BMStable[chn][ipl] = (chn - 62.5) * 1.25 ; //According to Rainer Joosten
                           break;
			   }
              }
         }
	 /*
	   BM05
	   Die Mitte des Detektors ist bei Kanal 29 (von 1 .. 64),
	   das spacing ist 2.5 mm.
	   D.h.:
	   (fuer Zaehlweise 0..63)
	   x = (ch - 28) * 2.5 mm
	   Natuerlich muss wieder fuer die Matrix die BMS Vorzeichen-Konvention 
	   beruecksichtigt werden.
	 */
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool  CsBMSrecons:: readcalib(void)
//
{
  int i,j;
  bool OK = true;
  //
  for(i=0; i<nhod; i++){
    tmresol[i]=1.; 
    // holds time resolution. Temporary. Sould be read from CsDetector
    disp[i]=tmresol[i]*tmresol[i];
  }
  for(i=0; i<3; i++)
      for(j=0; j<3; j++) {
	  readBMSconst(coeff_path[i*3+j],i,j);
      }
  if(Printlev>1) prntBMSconst();
//
//     Temporary. Should be replaced in a future by DB reading
//
  int ipl, ll;
  if(useTRAFFIC==0){         
    ifstream input_calib("bmst0.dat");
    if(input_calib.fail()){
      cerr << " Error during opening of input file "<< "bmst0.dat" << endl;
      OK=false;
      for(ipl=0; ipl<nhod; ipl++){
	for(ll=0; ll<hodsize[ipl]; ll++) TMt0[ll][ipl]=0;
      }
    }
    else {
      for(ipl=0; ipl<nhod; ipl++){
	for(ll=0; ll<hodsize[ipl]; ll++) input_calib >>TMt0[ll][ipl];
      }
      if(!input_calib.good()) OK=false;
    }
    input_calib.close();
    
    if(Printlev>1){         
      cout << endl << "BMS TDC t0s" << endl;
      prntclb(cout, TMt0);
    }
  }
//
  for(int i=0; i<nhod; i++) disp[i]=tmresol[i]*tmresol[i];
  return OK;
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: init(void)
{
  static bool first = true;
  if ( !first ) return;
  first = false;

  int i;
  string name;
  CsOpt* opt = CsOpt::Instance();
  const list <CsDetector*> &det = CsGeom::Instance()->getDetectors();
  for ( i=0; i<nhod; i++ ) {
    if ( !opt->getOpt( hodnames[i], "Name", name ) )
      CsErrLog::Instance()->msg( elFatal, __FILE__, __LINE__,
				 "Could not getOpt( hodnames[%i]=%s )!  Check CORAL option file!", i, hodnames[i].c_str() );

    Id[i] = 0;
    for ( list<CsDetector*>::const_iterator idet = det.begin(); idet != det.end(); idet++ ) {
      CsDetector* detpnt = *idet;
      if( detpnt ) {
	if(detpnt->getName()==name) {
	  //                          if(detpnt->getUnit()!=(i/2+1)) continue;  // Attention: for test only
	  Id[i] = detpnt;
	  //cout<<"CsBMSrecons: Found "<<name<<endl;
	  //                          hodsize[i] = detpnt->getNWir();  //Attention - test
	}
      } 
    }
    if ( !Id[i] )
      CsErrLog::Instance()->msg( elFatal, __FILE__, __LINE__,
				 "Detector plane %s not found!  Check CORAL option file and detectors.dat!", name.c_str() );
  }

  // BMStbl();
  return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons::inithist()
{
  for(int i=0; i<nhod_old; i++){
    HiHtime1[i] = NULL;
    HiHtime2[i] = HiHpro4[i] = HiYpro4[i] = HiHpro1[i] = NULL;
    HiTmVsHt[i]=NULL;
  }
  for(int i=0; i<nhod; i++){
    HiNhit[i] = HiHtime[i] = HiFired[i] = HiNHitafter[i] = NULL;
    HiPlanesCheck[i] = HiProf[i]  = NULL;
  }
  HiCleanFired.resize(nhod,NULL);
  //1D:
  Hit0t1dif =NULL;
  Hit0t2dif =NULL;
  Hit0t3dif =NULL;
  HiBmom4   =NULL;
  HiBmom3   =NULL;
  HiTMchi4  =NULL;
  HiTMchi3  =NULL;
  HiSPchi4  =NULL;
  HiSPchi3  =NULL;
  HiTime4   =NULL;
  HiTime3   =NULL;
  HiNhtev   =NULL;
  HiN4h3h   =NULL;
  HiBmom1   =NULL;
  HiTMchi1  =NULL;
  HiSPchi1  =NULL;
  HiTime1   =NULL;
  BMSPlanesfired  =NULL;
  NuTrpEv   =NULL;
  NuTrpEvClean =NULL;
  Delta1H1  =NULL;
  Delta1H2  =NULL;
  Delta1H3  =NULL;
  Delta1H4  =NULL;
  DeltaG    =NULL;
  BMSPlanesfiredNT  =NULL;
  BMSPlanesfiredNTcl  =NULL;
  BMSPlanesfiredTr =NULL;
  BMSdtTRX  =NULL;
  BMSeff0   =NULL;
  BMSeff1   =NULL;
  BMSeff2   =NULL;
  Delta     =NULL;
  DeltaT    =NULL;
  DeltaCut  =NULL;
  DeltaTclean =NULL;
  DeltaTcleanBMSnew =NULL;
  P1Clus    =NULL;
  P2Clus    =NULL;
  P3Clus    =NULL;
  P4Clus    =NULL;
  NeventsCl =NULL;
  RecPerfVStrig=NULL;
  TimeRes1_0102=NULL;
  TimeRes1_0304=NULL;
  TimeRes1_0104=NULL;
  TimeRes2_0102=NULL;
  TimeRes2_0304=NULL;
  TimeRes2_0104=NULL;
  HitBMS03VsHitBMS04=NULL;
  HitBMS01VsHitBMS02=NULL;
  HitBMS01VsHitBMS04=NULL;
  HitBMS04VsHitBMS03y=NULL;
  HitBMS02VsHitBMS01y=NULL;
  HitBMS03VsHitBMS02y=NULL;
  HitBMS04VsHitBMS03y3=NULL;
  HitBMS02VsHitBMS01y3=NULL;
  HitBMS03VsHitBMS02y3=NULL;
  EvSiVStri=NULL;
  CommonHitsHist1=NULL;
  CommonHitsHist2=NULL;
  CommonHitsHist3=NULL;

  //2D:
  Hit0t1 = Hit0t2 = Hit0t3 =HiN4h3h=NHitP2vsP1=NHitP2vsP3= NHitP4vsP3=NULL;
  ThetaYVsMom=NULL;
  TrackDec=NULL;
  TchiqDiff=NULL;
  TchiqDiff3=NULL;
  BadDec4=NULL;
  BadDec3=NULL;
  RescueDeltaMom=NULL;
  RescueMom2H=NULL;
  if(clusterIt!=0){
    Cl_chDiff=NULL;
    Cl_tmDiff=NULL;
    Cl_clSize=NULL;  
    for(int i=0;i<nhod_old;i++){
      Cl_HiProf[i]=NULL;
      Cl_HiHtime[i]=NULL;
      Cl_HiNhit[i]=NULL;
    }
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons::bookhist1()
{
// Book histograms for histolevel 1
//
  static bool first = true;
  if(first){
    first = false;
    int ii, nh = 1;
    char titl[80], name[40];
    string pathname =  "/BeamRecons/BMSrecond";
    CsHistograms::SetCurrentPath(pathname);

    sprintf(name,"BMS%4.4i",40);
    mH1[0]  = new CsHist1D(name, "BMS momemtum, 4hits tracks ", 100, 120, 220);
    sprintf(name,"BMS%4.4i",41);
    mH1[1]  = new CsHist1D(name, "BMS momemtum, 3hits trakcs ", 100, 120, 220);
    sprintf(name,"BMS%4.4i",101);
    mH1[2]  = new CsHist1D(name, "BMS Number of Tracks per event", 100, 0, 100);
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons::bookhist()
{
// Book histograms
//
  
  static bool first = true;
  if(first){
    first = false;
    int ii, nh = 1;
    char titl[80], name[40];
    string pathname =  "/BeamRecons/BMSrecond";
    CsHistograms::SetCurrentPath(pathname);
//

	for(int i=0; i<nhod; i++){
	  sprintf(titl,"Number of Hits, BMS%4.4i, after TimeCuts (not for method 3)",i+1);
	  sprintf(name,"BMShitsAfter%4.4i",i+1);
	  HiNHitafter[i] = new CsHist1D(name, titl, 65, 0, 64);
	  
          sprintf(name,"BMShits%4.4i",i+1);
          sprintf(titl,"Number of hits, BMS%i, raw events",i+1);
          HiNhit[i] = new CsHist1D(name, titl, 71, 0, 70);
	
          sprintf(name,"BMSchanels%4.4i",i+1);
          sprintf(titl,"Fired channels, BMS%i, raw events",i+1);
          HiFired[i] = new CsHist1D(name, titl, 131, 0, 130);

	  sprintf(name,"BMStime%4.4i",i+1);
	  sprintf(titl,"Hit time distr., BMS%i, raw events",i+1);
	  HiHtime[i] = new CsHist1D(name, titl, 200, -100, 100);

	  sprintf(name,"BMS_PlanesCheck%2.2i",i+1);
	  sprintf(titl,"Fired channels., BMS%i, events with beam",i+1);
	  HiPlanesCheck[i] = new CsHist1D(name, titl, 131, 0, 130);
	  
	  sprintf(name,"BMScleanChanels%4.4i",i+1);
	  sprintf(titl,"Fired channels, BMS%i, cleaned events",i+1);
	  HiCleanFired[i] = new CsHist1D(name, titl, 131, 0, 130);

	  sprintf(name,"BMScleanProfile%4.4i",i+1);
	  sprintf(titl,"Beam y-profile, BMS%i, cleaned events",i+1);
	  if(i < nhod_old) HiProf[i] = new CsHist1D(name, titl, 52, -130, 130);
	  else HiProf[i] = new CsHist1D(name, titl, 104, -130, 130);
	} 
	
      for(int i=0; i<nhod_old; i++){
        ii=i+1;

        sprintf(name,"BMS%4.4i",i+17);
        sprintf(titl,"Hit time distribution, BMS%i, 1hit/hod events",ii);
        HiHtime1[i] = new CsHist1D(name, titl, 60, -30, 30);
        sprintf(name,"BMS%4.4i",i+21);
        sprintf(titl,"Hittime distr.,BMS%i, 1hit/hod evnt, SPchiq cut",ii);
        HiHtime2[i] = new CsHist1D(name, titl, 60, -30, 30);
        //25->28
	sprintf(name,"BMS%4.4i",i+25);
        sprintf(titl,"Time vs channel number, BMS%i, 1hit/hod events",ii);
        HiTmVsHt[i] =  new CsHist2D(name, titl, 64, 0, 64, 40, -10, 10);
	//35->38
	//sprintf(name,"BMS%4.4i",i+35);
	//sprintf(titl,"Hit time distr., BMS%i, post time-cut",ii);
	//HiHtimeTCut[i] = new CsHist1D(name, titl, 700, -70, 70);
//
        sprintf(name,"BMS%4.4i",54+i);
        sprintf(titl,"Hit distrb. for 4hit tracks, BMS%i",ii);
        HiHpro4[i] =  new CsHist1D(name, titl, 64, 0, 64);
        //sprintf(name,"BMS%4.4i",58+i);
        //sprintf(titl,"Hit distrb. for 3hit tracks, BMS%i",ii);
        //HiHpro3[i] =  new CsHist1D(name, titl, 64, 0, 64);
        sprintf(name,"BMS%4.4i",62+i);
        sprintf(titl,"Track Y-profile, BMS%i",ii);
        HiYpro4[i] =  new CsHist1D(name, titl, 64, -32*5, 32*5);
//
        sprintf(name,"BMS%4.4i",75+i);
        sprintf(titl,"Hit distrb. for 1hit/hod tracks, BMS%i",ii);
        HiHpro1[i] =  new CsHist1D(name, titl, 64, 0, 64);
      }
//
        nh=29;
        sprintf(name,"BMS%4.4i",nh++);
        Hit0t1 = new CsHist2D(name, "BMS1 hittime vs BMS2 hittime, 1hit/hod events",
                       40, -20, 20, 40, -20, 20);
        sprintf(name,"BMS%4.4i",nh++);
        Hit0t2 = new CsHist2D(name, "BMS1 hittime vs BMS3 hittime, 1hit/hod events",
                       40, -20, 20, 40, -20, 20);
        sprintf(name,"BMS%4.4i",nh++);
        Hit0t3 = new CsHist2D(name, "BMS1 hittime vs BMS4 hittime, 1hit/hod events",
                       40, -20, 20, 40, -20, 20);
        sprintf(name,"BMS%4.4i",nh++);
        Hit0t1dif  = new CsHist1D(name,"BMS1 hittime - BMS2 hittime, 1hit/hod events", 50, -20, 20);
        sprintf(name,"BMS%4.4i",nh++);
        Hit0t2dif  = new CsHist1D(name,"BMS1 hittime - BMS3 hittime, 1hit/hod events", 50, -20, 20);
        sprintf(name,"BMS%4.4i",nh++);
        Hit0t3dif  = new CsHist1D(name,"BMS1 hittime - BMS4 hittime, 1hit/hod events", 50, -20, 20);
//
        nh=40;
        sprintf(name,"BMS%4.4i",nh++);
	HiBmom4  = new CsHist1D(name, "BMS momemtum, 4hits tracks ", 100, 120, 220);
        sprintf(name,"BMS%4.4i",nh++);
        HiBmom3  = new CsHist1D(name, "BMS momemtum, 3hits trakcs ", 100, 120, 220);
        sprintf(name,"BMS%4.4i",nh++);
        HiTMchi4 = new CsHist1D(name, "BMS track time chiq, 4 hit tracks", 100, 0, 100);
        sprintf(name,"BMS%4.4i",nh++);
        HiTMchi3 = new CsHist1D(name, "BMS track time chiq, 3 hit tracks", 100, 0, 100);
        sprintf(name,"BMS%4.4i",nh++);
        HiSPchi4 = new CsHist1D(name, "BMS Track space chiq, 4 hit tracks", 100, 0, 100);
        sprintf(name,"BMS%4.4i",nh++);
        HiSPchi3 = new CsHist1D(name, "BMS Track space chiq, 3 hit tracks", 100, 0, 100);
        sprintf(name,"BMS%4.4i",nh++);
        HiTime4  = new CsHist1D(name, "BMS track time, 4 hit tracks", 40, -20, 20);
        sprintf(name,"BMS%4.4i",nh++);
        HiTime3  = new CsHist1D(name, "BMS track time, 3 hit tracks", 40, -20, 20);
        sprintf(name,"BMS%4.4i",nh++);
        HiNhtev  = new CsHist1D(name, "Number of BMS hits/track", 5, 0, 5);
        sprintf(name,"BMS%4.4i",nh++);
        HiN4h3h  = new CsHist2D(name, "Number of 4hit tracks vs >numb. of 3hit tracks",
                       15, 0, 14, 15, 0, 14);
//
        nh=70;
        sprintf(name,"BMS%4.4i",nh++);//70
        HiBmom1  = new CsHist1D(name, "BMS momemtum, 1hit/hod events", 50, 120, 220);
        sprintf(name,"BMS%4.4i",nh++);//71
        HiTMchi1 = new CsHist1D(name, "BMS track time chiq, 1hit/hod events", 100, 0, 100);
        sprintf(name,"BMS%4.4i",nh++);//72
        HiSPchi1 = new CsHist1D(name, "BMS space chiq, 1hit/hod events", 100, 0, 100);
        sprintf(name,"BMS%4.4i",nh++);//73
        HiTime1  = new CsHist1D(name, "BMS track time, 1hit/hod events", 40, -20, 20);
	sprintf(name,"BMS%4.4i",94);
	EvSiVStri= new CsHist2D(name, "trigger vs eventsize (no recble)",80, 0, 80, 7, 0, 7);
	sprintf(name,"BMS%4.4i",95);
	BMSPlanesfiredTr=new CsHist1D(name, "Spectrum of fired planes VS trigger",100,0,100);
	sprintf(name,"BMS%4.4i",96);
	RecPerfVStrig=new CsHist1D(name, "RecPerformance vs trigger",100,0,100); 

	sprintf(name,"BMS%4.4i",97);
	BMSPlanesfiredNTcl  = new CsHist1D(name, "BMS # of fired planes(ntrack=0, clean Trigger)", 5, 0, 5);
	//98
	sprintf(name,"BMS%4.4i",98);
	NuTrpEvClean= new CsHist1D(name, "BMS Number of Tracks per event (goodTrigger)", 100, 0, 100);
	//99
	sprintf(name,"BMS%4.4i",99);
	BMSeff1  = new CsHist1D(name, "BMS Efficiency (goodTrigger)", 5, 0, 5);
	sprintf(name,"BMS%4.4i",100);
        BMSPlanesfired  = new CsHist1D(name, "BMS # of fired planes", 5, 0, 5);
	sprintf(name,"BMS%4.4i",101);
	NuTrpEv  = new CsHist1D(name, "BMS Number of Tracks per event", 100, 0, 100);	
       	//102->105
	//*** moved to another loop - conrad *** 
	//for(int m=0; m<nhod_old; m++){
	//  sprintf(titl,"Number of Hits, BMS%4.4i, after TimeCuts",m+1);
	//  sprintf(name,"BMS%4.4i",102+m);
	//  HiNHitafter[m] = new CsHist1D(name, titl, 65, 0, 64);
	//}
	sprintf(name,"BMS%4.4i",106);
        Delta1H1  = new CsHist1D(name, "BMS1 DeltaT(HitT<->MeanT)OneHitHod", 20, 0, 20);
	sprintf(name,"BMS%4.4i",107);
	Delta1H2  = new CsHist1D(name, "BMS2 DeltaT(HitT<->MeanT)OneHitHod", 20, 0, 20);
	sprintf(name,"BMS%4.4i",108);
	Delta1H3  = new CsHist1D(name, "BMS3 DeltaT(HitT<->MeanT)OneHitHod", 20, 0, 20);
	sprintf(name,"BMS%4.4i",109);
	Delta1H4  = new CsHist1D(name, "BMS4 DeltaT(HitT<->MeanT)OneHitHod", 20, 0, 20);
	sprintf(name,"BMS%4.4i",110);
	DeltaG  = new CsHist1D(name,"All_BMS DeltaT(HitT<->MeanT)OneHitHod", 20, 0, 20);
	sprintf(name,"BMS%4.4i",111);
        BMSPlanesfiredNT  = new CsHist1D(name, "BMS # of fired planes(ntrack=0)", 5, 0, 5);
	sprintf(name,"BMS%4.4i",112);
	BMSdtTRX  = new CsHist1D(name, "BMS timediff between tracks", 320,0, 80);
	sprintf(name,"BMS%4.4i",113);
	BMSeff0  = new CsHist1D(name, "BMS Efficiency", 18, 0,18);
	sprintf(name,"BMS%4.4i",114);
	Delta  = new CsHist1D(name, "BMS DeltaT(HitT<->MeanT)", 200, -100, 100);
	sprintf(name,"BMS%4.4i",115);
	DeltaT = new CsHist1D(name, "BMS DeltaT of all hit-combinations", 80, -5, 5);
	sprintf(name,"BMS%4.4i",116);
	DeltaCut  = new CsHist1D(name, "BMS DeltaT(HitT<->MeanT) after HitTimeDiffCut", 200, -100, 100);
	sprintf(name,"BMS%4.4i",117);
	NHitP2vsP1= new CsHist2D(name, "BMS_P2 Mult vs BMS_P1 Mult",50,0,50,50,0,50);
	sprintf(name,"BMS%4.4i",118);
	NHitP2vsP3= new CsHist2D(name, "BMS_P2 Mult vs BMS_P3 Mult",50,0,50,50,0,50);
	sprintf(name,"BMS%4.4i",119);
	NHitP4vsP3= new CsHist2D(name, "BMS_P4 Mult vs BMS_P3 Mult",50,0,50,50,0,50);
	sprintf(name,"BMS%4.4i",120);
	DeltaTclean = new CsHist1D(name, "BMS DeltaT of all hit-combinations (clean)", 128, -8, 8);
	sprintf(name,"BMS%4.4i",121);
	P1Clus = new CsHist2D(name, "Plane1: deltaT vs deltaCh",36, 0, 180, 300, -30, 30);
	sprintf(name,"BMS%4.4i",122);
	P2Clus = new CsHist2D(name, "Plane2: deltaT vs deltaCh",36, 0, 180, 300, -30, 30);
	sprintf(name,"BMS%4.4i",123);
	P3Clus = new CsHist2D(name, "Plane3: deltaT vs deltaCh",36, 0, 180, 300, -30, 30);
	sprintf(name,"BMS%4.4i",124);
	P4Clus = new CsHist2D(name, "Plane4: deltaT vs deltaCh",36, 0, 180, 300, -30, 30);
	sprintf(name,"BMS%4.4i",125);
	NeventsCl= new CsHist1D(name, "nEvents with clusters", 5, 0,5);
	sprintf(name,"BMS%4.4i",126);
	HitBMS01VsHitBMS02= new CsHist2D(name, "Hitmap: BMS02 vs BMS01 (ch)",64, 0, 64, 64, 0, 64);
	sprintf(name,"BMS%4.4i",127);
	HitBMS03VsHitBMS04= new CsHist2D(name, "Hitmap: BMS03 vs BMS04 (ch)",64, 0, 64, 64, 0, 64);
	sprintf(name,"BMS%4.4i",128);
	HitBMS01VsHitBMS04= new CsHist2D(name, "Hitmap: BMS01 vs BMS04 (ch)",64, 0, 64, 64, 0, 64);
	sprintf(name,"BMS%4.4i",129);
	TimeRes1_0102= new CsHist1D(name, "TimeRes: BMS01_34,BMS02_33", 160, -10,10);
	sprintf(name,"BMS%4.4i",130);
	TimeRes2_0102= new CsHist1D(name, "TimeRes: BMS01_27,BMS02_25", 160, -10,10);
	sprintf(name,"BMS%4.4i",131);
	TimeRes1_0304= new CsHist1D(name, "TimeRes: BMS03_36,BMS04_33", 160, -10,10);
	sprintf(name,"BMS%4.4i",132);
	TimeRes2_0304= new CsHist1D(name, "TimeRes: BMS03_27,BMS04_30", 160, -10,10);
	sprintf(name,"BMS%4.4i",133);
	TimeRes1_0104= new CsHist1D(name, "TimeRes: BMS01_34,BMS04_33", 160, -10,10);
	sprintf(name,"BMS%4.4i",134);
	TimeRes2_0104= new CsHist1D(name, "TimeRes: BMS01_27,BMS04_30", 160, -10,10);
	sprintf(name,"BMS%4.4i",135);
	ThetaYVsMom= new CsHist2D(name, "ThetaYVsMom: 4-hit Tracks", 200, 0,200,100,-0.020,0.020);
	sprintf(name,"BMS%4.4i",136);
	TrackDec= new CsHist1D(name, "Track Decisions", 12, 0,12);
	sprintf(name,"BMS%4.4i",137);
	TchiqDiff= new CsHist1D(name, "Tchiq differences: 4hit tracks", 100, 0,10);
	sprintf(name,"BMS%4.4i",138);
	TchiqDiff3= new CsHist1D(name, "Tchiq differences: 3hit tracks", 100, 0,10);
	sprintf(name,"BMS%4.4i",139);
	BadDec4= new CsHist1D(name, "Number of bad decisions: 4hit tracks", 50, 0,50);
	sprintf(name,"BMS%4.4i",140);
	BadDec3= new CsHist1D(name, "Number of bad decisions: 3hit tracks", 50, 0,50);
	sprintf(name,"BMS%4.4i",141);
	BadDec3= new CsHist1D(name, "Number of bad decisions: 3hit tracks", 50, 0,50);
	sprintf(name,"BMS%4.4i",142);
	CommonHitsHist1= new CsHist1D(name, "Number inv tracks per plane", 60, 0,60);
	sprintf(name,"BMS%4.4i",143);
	CommonHitsHist2= new CsHist1D(name, "Trackmult in common hits", 40, 0,40);
	sprintf(name,"BMS%4.4i",144);
	CommonHitsHist3= new CsHist2D(name, "inv Tracks(in common hits) vs nBMSpassed",20,0,20,20,0,20);
	sprintf(name,"BMS%4.4i",145);
	BMSeff2  = new CsHist1D(name, "BMS Efficiency (goodTrigger,RefTime)", 18, 0, 18);
	sprintf(name,"BMS%4.4i",146);
	HitBMS02VsHitBMS01y= new CsHist2D(name, "Hitmap: BMS02 vs BMS01 (y)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",147);
	HitBMS04VsHitBMS03y= new CsHist2D(name, "Hitmap: BMS04 vs BMS03 (y)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",148);
	HitBMS03VsHitBMS02y= new CsHist2D(name, "Hitmap: BMS03 vs BMS02 (y)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",149);
	HitBMS02VsHitBMS01y3= new CsHist2D(name, "Hitmap: BMS02 vs BMS01 (y,3-hit-track)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",150);
	HitBMS04VsHitBMS03y3= new CsHist2D(name, "Hitmap: BMS04 vs BMS03 (y,3-hit-track)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",151);
	HitBMS03VsHitBMS02y3= new CsHist2D(name, "Hitmap: BMS03 vs BMS02 (y,3-hit-track)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",152);
	RescueDeltaMom= new CsHist2D(name, "Rescue: mom(4hit)-mom(2hit)",200,-50,50,5,0,5);
	sprintf(name,"BMS%4.4i",153);
	RescueMom2H= new CsHist2D(name, "Rescue: mom(2hit)",200,0,200,5,0,5);

	sprintf(name,"BMS%4.4i",154);
	BeamH2D[0] = new CsHist2D(name, "MultBeam: BMS02 vs BMS01 (y)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",155);
	BeamH2D[1] = new CsHist2D(name, "MultBeam: BMS04 vs BMS03 (y)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",156);
	BeamH2D[2] = new CsHist2D(name, "MultBeam: deltaT vs deltaCh (BMS01)",36,0,180,300,-150,150);
	sprintf(name,"BMS%4.4i",157);
	BeamH2D[3] = new CsHist2D(name, "MultBeam: deltaT vs deltaCh (BMS02)",36,0,180,300,-150,150);
	sprintf(name,"BMS%4.4i",158);
	BeamH2D[4] = new CsHist2D(name, "MultBeam: deltaT vs deltaCh (BMS03)",36,0,180,300,-150,150);
	sprintf(name,"BMS%4.4i",159);
	BeamH2D[5] = new CsHist2D(name, "tracks: BMS2 vs BMS01 (y)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",160);
	BeamH2D[6] = new CsHist2D(name, "tracks: BMS3 vs BMS02 (y)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",161);
	BeamH2D[7] = new CsHist2D(name, "tracks: BMS4 vs BMS03 (y)",80,-200,200,80,-200,200);
	sprintf(name,"BMS%4.4i",162);
	BeamH2D[8] = new CsHist2D(name, "tracks: BMS4 vs BMS01 (y)",80,-200,200,80,-200,200);

	sprintf(name,"BMS%4.4i",163);
	DeltaTcleanBMSnew = new CsHist1D(name, "BMS DeltaT of all hit-combinations with BM05 and BM06 (clean)", 384, -24, 24);

	BeamH2D[9] = new CsHist2D("BeamH2D_09","BMS-plane vs fired channels, raw events", 65, 0, 64,5,0,5);
	BeamH2D[10] = new CsHist2D("BeamH2D_10","BMS planes vs Beam y-profile, raw events" , 64, -32*5, 32*5,5,0,5);
	BeamH2D[11] = new CsHist2D("BeamH2D_11","BMS planes vs Hit time distr., raw events",200, -100, 100,5,0,5);
	BeamH2D[12] = new CsHist2D("BeamH2D_12","BMS planes vs Hit time distr., post time-cut", 700, -70, 70,5,0,5);
	BeamH2D[13] = new CsHist2D("BeamH2D_13","BMS planes vs Number of hits, raw events", 70, 0, 70,5,0,5);
	BeamH2D[14] = new CsHist2D("BeamH2D_14","BMS planes vs Number of Hits, after TimeCuts", 65, 0, 65,5,0,5);
	BeamH2D[15] = new CsHist2D("BeamH2D_15","BMS planes vs Hit distrb. for 3hit tracks", 64, 0, 64,5,0,5);
	//
	BeamH1D[0] = new CsHist1D("BeamH1D_0", "MultBeam: momentum difference",60,-30,30);
	BeamH1D[1] = new CsHist1D("BeamH1D_1", "Track with BM05: momentum",100,120,220);
	BeamH1D[2] = new CsHist1D("BeamH1D_2", "Track with BM06: momentum",100,120,220);
	
	if(clusterIt!=0){
	  string pathname =  "/BeamRecons/BMSrecond/Clustering";
	  CsHistograms::SetCurrentPath(pathname);
	  sprintf(name,"BMS%4.4i",00);
	  Cl_chDiff= new CsHist1D(name, "Channel differences in clusters", 10, 0,10);
	  sprintf(name,"BMS%4.4i",01);
	  Cl_tmDiff= new CsHist1D(name, "Time differences in clusters", 10, 0,10);
	  sprintf(name,"BMS%4.4i",02);
	  Cl_clSize=new CsHist1D(name, "Number of hits in a cluster", 20, 0,20);
	  
	  for(int i=0; i<nhod_old; i++){
	    ii=i+1;
	    sprintf(name,"BMS%4.4i",i+3);
	    sprintf(titl,"Number of hits, BMS%i, raw events",ii);
	    Cl_HiNhit[i] = new CsHist1D(name, titl, 71, 0, 70);
	    sprintf(name,"BMS%4.4i",i+7);
	    sprintf(titl,"Beam y-profile, BMS%i, raw events",ii);
	    Cl_HiProf[i] = new CsHist1D(name, titl, 64, -32*5, 32*5);
	    sprintf(name,"BMS%4.4i",i+11);
	    sprintf(titl,"Hit time distr., BMS%i, raw events",ii);
	    Cl_HiHtime[i] = new CsHist1D(name, titl, 200, -100, 100);
	  }
	}
	pathname =  "/BeamRecons/BMSrecond/TriggerStudy";
	CsHistograms::SetCurrentPath(pathname);
	trigH[0]  = new CsHist1D("BMS_trig_00","BMS Number of Tracks per event (all unamb)", 100, 0, 100);
	trigH[1]  = new CsHist1D("BMS_trig_01","BMS Number of Tracks per event (ITCV)", 100, 0, 100);
	trigH[2]  = new CsHist1D("BMS_trig_02","BMS Number of Tracks per event (MTCV)", 100, 0, 100);
	trigH[3]  = new CsHist1D("BMS_trig_03","BMS Number of Tracks per event (LTCV)", 100, 0, 100);
	trigH[4]  = new CsHist1D("BMS_trig_04","BMS Number of Tracks per event (OT)", 100, 0, 100);
	trigH[5]  = new CsHist1D("BMS_trig_05","BMS Number of Tracks per event (CT)", 100, 0, 100);
	trigH[6]  = new CsHist1D("BMS_trig_06","BMS Number of Tracks per event (MTV)", 100, 0, 100);
	trigH[7]  = new CsHist1D("BMS_trig_07","BMS Number of Tracks per event (ITCV)", 100, 0, 100);
	trigH[8]  = new CsHist1D("BMS_trig_08","BMS Number of Tracks per event (MTCV)", 100, 0, 100);
	trigH[9]  = new CsHist1D("BMS_trig_09","BMS Number of Tracks per event (LTCV)", 100, 0, 100);
	trigH[10]  = new CsHist1D("BMS_trig_10","BMS Number of Tracks per event (OT)", 100, 0, 100);
	trigH[11]  = new CsHist1D("BMS_trig_11","BMS Number of Tracks per event (CT)", 100, 0, 100);
	trigH[12]  = new CsHist1D("BMS_trig_12","BMS Number of Tracks per event (MTV)", 100, 0, 100);
	trigH[13]  = new CsHist1D("BMS_trig_13","nComb=0:BMS number of tracks per event", 20, 0, 20);
	trigH[14]  = new CsHist1D("BMS_trig_14","ncut=1,BM01: HitTime-RefTime", 150,-15, 15);
	trigH[15]  = new CsHist1D("BMS_trig_15","ncut=1,BM02: HitTime-RefTime", 150,-15, 15);
	trigH[16]  = new CsHist1D("BMS_trig_16","ncut=1,BM03: HitTime-RefTime", 150,-15, 15);
	trigH[17]  = new CsHist1D("BMS_trig_17","ncut=1,BM04: HitTime-RefTime", 150,-15, 15);
	trigH[18]  = new CsHist1D("BMS_trig_18","nComb=0,!BMSok: nfired BMS", 10,0, 10);
	trigH[19]  = new CsHist1D("BMS_trig_19","ncut=1:for involved hits:HitTime-RefTime", 150,-15, 15);
	trigH[20]  = new CsHist1D("BMS_trig_20","nComb=0:HitTime-RefTime", 150,-15, 15);
	trigH[21]  = new CsHist1D("BMS_trig_21","ncut=1: differences between hittimes", 100,-10, 10);
	trigH[22]  = new CsHist1D("BMS_trig_22","SpecMom<cut,nComb=0,!BMSok: nfired BMS", 10,0, 10);
	
	trig2H[0]  = new CsHist2D("BMS_trig2_0", "nComb=0,!BMSok,nfired<3: nfired vs fired BMS-planes",5,0,5,5,0,5);
	trig2H[1]  = new CsHist2D("BMS_trig2_1", "nComb=0,!BMSok,nfired=2: !nfired1 vs !nfired2",5,0,5,5,0,5);
	trig2H[2]  = new CsHist2D("BMS_trig2_2", "SpecMom<cut,nComb=0,!BMSok,nfired=2: !nfired1 vs !nfired2",5,0,5,5,0,5);
	trig2H[3]  = new CsHist2D("BMS_trig2_3", "BMSi vs HitTime-RefTime",150,-15, 15,5,0,5);
	
	//
	//
	pathname =  "/BeamRecons/BMSrecond/RescueAlgo";
	CsHistograms::SetCurrentPath(pathname);
	Resc1D[0]  = new CsHist1D("BMS_resc1_00","Number of 2hit tracks", 100, 0, 100);
	Resc1D[1]  = new CsHist1D("BMS_resc1_01","Time-Chi2", 100, 0, 100);
	Resc1D[2]  = new CsHist1D("BMS_resc1_02","Track-Time", 160, -20, 20);
  }
} 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
void CsBMSrecons::prntclb(ostream &out, double calib[][NBMSHODS]) const
{
      int ipl, i, j;
      cout << "NGODee  "<< nhod<< endl;
      out << resetiosflags(ios::scientific);
      out << setiosflags(ios::fixed|ios::showpoint);
      out << setprecision(4); 
      for(ipl=0; ipl<nhod; ipl++){
            cout <<nhod <<" "<< hodsize[ipl] << endl;
        for(j=0; j<hodsize[ipl]; ) {
           for( i=0; (i<8)&&(j<hodsize[ipl]); i++, j++) out<<calib[j][ipl]<<"  ";
           out << endl;
        }
      }
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: rawevprint(ostream &out) const
//
//   print raw BMS event         
//
{
  int ipl, j, nchnl;
  double time;
// 
  out << resetiosflags(ios::scientific);
  out << setiosflags(ios::fixed|ios::showpoint);
  out << setfill(' ');
  out << setprecision(1) ;
  out << "========================= raw  BMS event ====================" << endl;
//
  list<CsDigit*>::iterator Idig;
  for(ipl=0; ipl<nhod; ipl++) {
     if(!Id[ipl]) continue;
     list<CsDigit*>digit = Id[ipl]->getMyDigits();
     out << hodnames[ipl] << ":  total number of hits: " << digit.size() 
           << "   Hits(channel/time):"<< endl;
     if(!digit.size()) continue;
     j=0;  
     for( Idig = digit.begin(); Idig !=digit.end(); Idig ++ ) {
          CsDigit* digWB = *Idig;
	  nchnl = digWB->getAddress();
    	  time  = digWB->getDatum();
          out << setw(10) << nchnl 
               << setw(7) << time;
          j++; 
          if(j==5) { j=0; out << endl;}
     }
     if(j!=0) out << endl;
  } 
  out << "------------------------- end raw  BMS event ----------------" << endl;
  return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: test(void)
{
      int ipl,i,j,k,l,npart,part[30],ntrackf;
      double tst[4][30];
//
            for(ipl=0; ipl<4; ipl++){
            input_file >> khit16[ipl];
              if(khit16[ipl]) {
	        for(i=0; i<khit16[ipl]; i++) input_file >> ibmz16[i][ipl]; 
	        for(i=0; i<khit16[ipl]; i++) input_file >> tbmz16[i][ipl];
              } 
            }
            nfired=0;
            for(i=0;i<nhod;i++) {if(khit16[i]>0) nfired++;}
//
            for(ipl=0; ipl<4; ipl++){
            cout << "BMS plane, khit: "<< ipl<< "   "<< khit16[ipl]<<endl;
              if(khit16[ipl]) {
	        for(i=0; i<khit16[ipl]; i++) cout << ibmz16[i][ipl]<<" "; 
                cout << endl;
	        for(i=0; i<khit16[ipl]; i++) cout << tbmz16[i][ipl]<<" ";
                cout << endl;
              } 
            }
//          
            input_file >> ntrackf;
            cout<< "BMS:: found  " << ntrackf << endl;
            if(ntrackf) {
               for(i=0; i<ntrackf; i++) {
                input_file  >> part[i] >> tst[0][i] 
                            >> tst[1][i] >> tst[2][i] >> tst[3][i]; 
                cout << part[i] << "  " << tst[0][i] << "  " 
                     << tst[3][i] << "  " << tst[1][i] << "  " << tst[2][i] << endl;; 

               }
            }
   cout << "end BMS input" << endl;
   tstdecode();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: tstdecode(void)
//
//   decode, calibrate and cut BMS data (use just for test)        
//
{
         int ipl, i, ll, nhit, nhit1;
         double tm;
         nfired=0; 
         for(ipl=0; ipl<nhod; ipl++){
            nhit=khit16[ipl];
            nhit1=0;    
            for(i=0; i<nhit; i++) {
                ll=ibmz16[i][ipl]-1;                           //  Changed to move to zero channel count base
                if(ll<=hodsize[ipl]) {
                  tm=tbmz16[i][ipl];
//
                  if(fabs(tm)>BMSpar[0]) continue;         //   time window cut   
                  zhit16[nhit1][ipl]=BMStable[ll][ipl];    //save good hit
                  tbmz16[nhit1][ipl]=tm;
                  ibmz16[nhit1][ipl]=ll;
                  iflc16[nhit1][ipl]=1;
                  nhit1++;
                }
            }
            khit16[ipl]=nhit1;
            if(nhit1) nfired++;
         } 
//
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: readBMSconst(string filename, int set1, int set2)
{
	 char bufor[200];
         ifstream input_file(filename.c_str());
	 if(input_file.fail()){
	    cerr << " Error during opening of input file "<< filename << endl;
	    CsErrLog::Instance()->mes(elFatal,"Error during opening of BMS-Coeff-File!");
            return;
         }
	 else {
	     input_file.getline(bufor, 200);
	     input_file >> p0;
	     
	     int i, j, ll;
	     for( i=0; i<4; i++){
	       input_file >> coeff_avsg16[set1][set2][i];
	     }
	     input_file.getline(bufor, 200);
	     for( i=0; i<4; i++){
	       for( j=0; j<4;j++){
		   input_file>>coeff_evec16[set1][set2][j][i];
	       }
	       input_file.getline(bufor, 200);
	     }
	     for( i=0; i<5; i++){
	       for( j=0; j<3;j++){
	         for(  ll=0; ll<4 ;ll++){
		   input_file >> coeff_coff16[set1][set2][ll][j][i] ;
                 }
		 input_file.getline(bufor, 200);
               }
	     }
	     input_file.close();
	 }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBMSrecons:: prntBMSconst(void)
{
    cout << resetiosflags(ios::scientific);
    cout << setiosflags(ios::fixed|ios::showpoint);
//    cout << header << endl;
    cout << "Pbeam = " << p0 << endl;

    int i, j, ll, a,b;

    for(a=0; a<3; a++){
	for(b=0; b<3; b++){
	    cout<<a<<", "<<b<<endl;

	    cout <<  " AVSG16:" << endl;
	    cout << setprecision(5);
	    for( i=0; i<4; i++){
		cout << setw(16) << coeff_avsg16[a][b][i];
	    }
	    cout << endl;

	    cout <<  " EVEC16:" << endl;
	    for( i=0; i<4; i++){
		for( j=0; j<4;j++){
		    cout << setw(16) << coeff_evec16[a][b][j][i];
		}
		cout << endl;
	    }

	    cout <<  " COFF16:" << endl;
	    for( i=0; i<5; i++){
		for( j=0; j<3;j++){
		    for(  ll=0; ll<4 ;ll++){
			cout << setw(16) << coeff_coff16[a][b][ll][j][i];
		    }
		    cout << endl;
		}
	    }

	}
    }
}

//Return raw BMS hits in selected plane. Takes as input references to two vectors.
//In one hit times are stored. In second decoded Y coordinate of hits are stored.
//Returns true if operation was successfull.
bool CsBMSrecons::getBMSRawHits(int plane, std::vector<double>& time, std::vector<double>& Y)
{
    if(plane<0 || plane>nhod) return false;
    for(int i=0; i<khit16[plane]; i++) {
	time.push_back(tbmz16[i][plane]);
	Y.push_back(-zhit16[i][plane]);
    }
    return true;
}

//Return info about hits that form specified track. Takes track id, and reference to vector.
//In vector ids of hits are stored. Then info about hits can be read from vectors returned by
//CsBMSrecons::getBMSRawHits. Returns true when operation was successfull.
bool CsBMSrecons::getBMStrackHits(int n, std::vector<int>& hits)
{
    if((n<0)||(n>=ntrack)) return false;
    for(int i=0; i<nhod; i++) {
	hits.push_back(TRKindx[n][i]);
    }
    return true;
}

// ***************************************************************************
// *****  setHitPattern: Set BMS bits in argument CsBeam's hit patterns. *****
// ***************************************************************************
bool CsBMSrecons::setHitPattern(CsBeam* beam, int iBMStrack, bool rescue){

   //few variables
   int word, bit;
   unsigned int det_position;

   //get data
   map<CsDetector*,int> det2bit = CsEvent::Instance()->getDetectorMap();
   unsigned int *expected = const_cast<unsigned int*>(beam->getExpectedDetsBitmap());
   unsigned int *fired = const_cast<unsigned int*>(beam->getFiredDetsBitmap());

   //loop over planes
   for(int ipl = 0; ipl < NBMSHODS; ipl++){
      
      //check usage of BM05/BM06
      if( (ipl == 4 && !useBM05) || (ipl == 5 && !useBM06) ) continue;

      //position in the detector map
      det_position = (*(det2bit.find(Id[ipl]))).second;

      word = int(det_position/32);
      bit =      det_position%32;

      //check range
      if( word >= CSTRACK_MAPSIZE ){
         cout << "CsBMSrecons::setHitPattern: Error: array _expectedDets/_firedDets is to small to store hit pattern." << endl;
         cout << "                            The size is " << CSTRACK_MAPSIZE << ", where it needs to be at least " << word+1 << endl;
         return false;
      }

      //set fired "bit" to 0
      fired[word]    = (fired[word] | 1<<bit) ^ 1<<bit;

      //set expected "bit" to 1
      expected[word] = expected[word] | 1<<bit;

      //check if plane fired
      if( !rescue && (TRKindx[iBMStrack][ipl] < 0 || TRKindx[iBMStrack][ipl] >= 9999) ) continue;
      if( rescue && (TRKindx_resc[iBMStrack][ipl] < 0 || TRKindx_resc[iBMStrack][ipl] >= 9999) ) continue;

      //set fired "bit" to 1
      fired[word] = fired[word] | 1<<bit;
   }

   return true;
}

