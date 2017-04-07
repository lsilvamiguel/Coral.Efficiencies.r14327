#include "CsGeom.h"
#include "RecCall.h"
#include "RecOpt.h"
#include  "CsOpt.h"
#include "CsField.h"
#include "CsEvent.h"
#include  "Recon.h"
#include <cstdio>

using namespace std;

extern "C"
void reconini_(int[] , int[], float[], float[], float[], 
               float[], float[], int[], float[],float[],int[]); 

extern "C"
void planes1_(int[], int*,float[], float[],float[], float[], float[], int[], float[],int[]);

extern "C"
void planes2_(int[], int*,float[], float[],float[], float[], float[], int[], float[],int[]);

extern "C"
void planes3_(int*, float[], float[],float[], float[], float[], int[], float[],int[],float[][3]);

void RecCall::ReconIni(){
  
  const int max1=100;
  const int max2=200;
  int plBeforeSM1,plBehindSM1,plBehindRich;
  float ZPla[max1],CosN[max1],SinN[max1],Angel[max1],SpResol[max1],PlEff[max1];
  float ZPla1[max2],CosN1[max2],SinN1[max2],Angel1[max2],SpResol1[max2],PlEff1[max2];
  float ZPla2[max1],CosN2[max1],SinN2[max1],Angel2[max1],SpResol2[max1],PlEff2[max1];
  int  Type[max1],UId[max1],Type1[max2],UId1[max2],Type2[max1],UId2[max1];
   plBeforeSM1=plBehindSM1=plBehindRich=0;
  float geom[2][3],xmag[2] ;
  char* s_mm="MM"; char* s_dc="DC"; char* s_gm="GM"; char* s_st="ST"; char* s_sfi="FI";
  char* s_omega1="PA"; char* s_omega2="PB"; char* s_omega3="PS"; char* s_mw1="MA";
  int DetType;
  

  RecOpt::getRecOptions();   //get RECON options
  list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
  list<CsDetector*>::iterator  ih;
  CsGeom*  g = CsGeom::Instance();
  CsField* f = g->getCsField();
  CsMagInfo* mag =  f->getMagInfo();
  Recon::ref().magType[0]=mag[1].zcm;
  Recon::ref().magType[1]=mag[2].zcm;
  // printf(" MAGNES %f %f \n",Recon::ref().magType[0],Recon::ref().magType[1]);
 
  int schema_number;
  CsOpt* opt=CsOpt::Instance();  opt->getOpt("","reconstruction schema",schema_number);
  if( schema_number == 1) {  Recon::ref().isScheme1=true;}      // only Recon "true==1"
  else if (schema_number == 2) { Recon::ref().isScheme1=false;} // Traffic+Recon "false==0"
  else { printf(" Unknown reconstruction scheme, check in option file \n"); exit(1);}
  printf ("Schema %d\n",Recon::ref().isScheme1);


  if( CsEvent::Instance()->isAMonteCarloEvent() && (RecOpt::McPar[0]==1 || RecOpt::McPar[0]==3) ){
    Recon::ref().isMonteCarlo=true;}
  else {Recon::ref().isMonteCarlo=false;}
  // printf("Monte Carlo  %d\n",Recon::ref().isMonteCarlo);

  for( ih=dets.begin(); ih!=dets.end(); ih++ ) {
    if((*ih)->getZcm() <  0.0 || (*ih)->getZcm() > mag[2].zcm) continue;
    DetType=-1;
    if (strncmp((*ih)->GetTBName().c_str(),s_mm,2)==0) DetType=111; 
    else if (strncmp((*ih)->GetTBName().c_str(),s_dc,2)==0)  DetType=222;
    else if (strncmp((*ih)->GetTBName().c_str(),s_gm,2)==0) { DetType=444;
      geom[0][2]=(float)(((*ih)->getXsiz())/1000);
      geom[1][2]=(float)(((*ih)->getYsiz())/1000);}
    else if (strncmp((*ih)->GetTBName().c_str(),s_st,2)==0) DetType=555;
    else if (strncmp((*ih)->GetTBName().c_str(),s_sfi,2)==0) { DetType=333;
     geom[0][1]=(float)(((*ih)->getXsiz())/1000) ;
     geom[1][1]=(float)(((*ih)->getYsiz())/1000);}
    else if  ( (strncmp((*ih)->GetTBName().c_str(),s_omega1,2)==0) ||
	       (strncmp((*ih)->GetTBName().c_str(),s_omega2,2)==0) ||  
	       (strncmp((*ih)->GetTBName().c_str(),s_omega3,2)==0) ) { DetType=666;
	         geom[0][0]=(float)(((*ih)->getXsiz())/1000) ;
		 geom[1][0]=(float)(((*ih)->getYsiz())/1000);}
    else if (strncmp((*ih)->GetTBName().c_str(),s_mw1,2)==0) ;
    else { printf("Unknown detector %s \n",(*ih)->GetTBName().c_str()); exit(1);}
    if (DetType < 0 ) continue;
      
    if( (*ih)->getZcm() < mag[1].zcm ) { plBeforeSM1++;  // planes before SM1
    ZPla[plBeforeSM1-1]=float((*ih)->getZcm());
    CosN[plBeforeSM1-1]=float((*ih)->getCosAng());
    SinN[plBeforeSM1-1]=float((*ih)->getSinAng());
    Angel[plBeforeSM1-1]=float((*ih)->getAng()); 
    PlEff[plBeforeSM1-1]=float((*ih)->getEff());
    UId[plBeforeSM1-1]=(*ih)->GetID();
    Type[plBeforeSM1-1]=DetType;       
    if((*ih)->hasDrift()){ 
      SpResol[plBeforeSM1-1] = float((*ih)->getVel()*(*ih)->getSpSli());} // for drifts
    else{ SpResol[plBeforeSM1-1] = float(fabs((*ih)->getWirP()/sqrt(12.)));}  // for proportional chambers
    
    }
    else if( (*ih)->getZcm() < 7000.){ plBehindSM1++;  // between SM1 and RICH 
    ZPla1[plBehindSM1-1]=float((*ih)->getZcm());
    CosN1[plBehindSM1-1]=float((*ih)->getCosAng());
    SinN1[plBehindSM1-1]=float((*ih)->getSinAng());
    Angel1[plBehindSM1-1]=float((*ih)->getAng());           
    PlEff1[plBehindSM1-1]=float((*ih)->getEff());
    UId1[plBehindSM1-1]=(*ih)->GetID();
    Type1[plBehindSM1-1]=DetType;
    if((*ih)->hasDrift()){
      SpResol1[plBehindSM1-1] = float((*ih)->getVel()*(*ih)->getSpSli());} // for drifts
    else{ SpResol1[plBehindSM1-1] = float(fabs((*ih)->getWirP()/sqrt(12.)));}   // for proportional chambers
    
    }
    else if ((*ih)->getZcm() > 7000. &&  (*ih)->getZcm() < mag[2].zcm){ plBehindRich++; // beetwen RICH and SM2 
    ZPla2[plBehindRich-1]=float((*ih)->getZcm());
    CosN2[plBehindRich-1]=float((*ih)->getCosAng());
    SinN2[plBehindRich-1]=float((*ih)->getSinAng());
    Angel2[plBehindRich-1]=float((*ih)->getAng()); 
    PlEff2[plBehindRich-1]=float((*ih)->getEff());
    UId2[plBehindRich-1]=(*ih)->GetID();
    Type2[plBehindRich-1]=DetType;
    if((*ih)->hasDrift()){
      SpResol2[plBehindRich-1] = float((*ih)->getVel()*(*ih)->getSpSli()); } // for drift chambers
    else{ SpResol2[plBehindRich-1] =float(fabs((*ih)->getWirP()/sqrt(12.)));} // for proportional chambers  
    } // end if Det in zone
  }  // end of Det loop
  
  xmag[0]=(mag[1].zcm)/1000;  xmag[1]=0.88;

  if (RecOpt::Switch <= 1  ) {  // recon initialization
    reconini_(RecOpt::McPar, RecOpt::inPar, RecOpt::roadPar ,
	      RecOpt::ReProj, RecOpt::roadLengh, RecOpt::chiPar, 
	      RecOpt::ReTarget ,RecOpt::nbPlanes,
	      RecOpt::activeRegion ,xmag,
	      RecOpt::blacklist); 
    
    planes1_(RecOpt::deadPlanes1,&plBeforeSM1,ZPla,Angel,SinN,CosN,SpResol,Type,PlEff,UId);
    planes2_(RecOpt::deadPlanes2,&plBehindSM1,ZPla1,Angel1,SinN1,CosN1,SpResol1,Type1,PlEff1,UId1);
    planes3_(&plBehindRich,ZPla2,Angel2,SinN2,CosN2,SpResol2,Type2,PlEff2,UId2,geom);}
  else if ( RecOpt::Switch >  1 ) { cout << "--->Recon is switched off by RecOpt::Switch option<---" << endl;} 
  }




