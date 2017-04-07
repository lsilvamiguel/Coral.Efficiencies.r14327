// $Id:

/*!
   \file    CsRCIFtoDisplay.cc
   \brief   CORAL Event Display Package.
   \version $Revision: 1.3 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#include "TMath.h"
#include "CsRCIFtoDisplay.h"
#include "CsRichOneDisplay.h"

#include "coral_config.h"

using namespace std;
using namespace CLHEP;

CsRCIFtoDisplay::CsRCIFtoDisplay(){
}

CsRCIFtoDisplay::~CsRCIFtoDisplay(){}

void CsRCIFtoDisplay::Trans(CsRichOneDisplay *event){

  CsEvent * fullevent = CsEvent::Instance();
  CsRCDisplay* display = CsRichOneDisplay::Instance()->getDisplay();

  display->SetRunNumber(fullevent->getRunNumber());  
  display->SetEventNumberInRun(fullevent->getEventNumberInRun());  
  display->SetEventNumberInBurst(fullevent->getEventNumberInBurst());  
  /*
  static int run_num=31815;
  static int eve_num=0;
  static int eve_bur=0;
  eve_num++;
  eve_bur=eve_num%1000;
  display->SetRunNumber(run_num);  
  display->SetEventNumberInRun(eve_num);  
  display->SetEventNumberInBurst(eve_bur);  
  */

}

TRotMatrix *ConvertMatrix(HepMatrix HepRotMat)
{
  Double_t ThetaX=0,ThetaY=0,ThetaZ=0,PhiX=0,PhiY=0,PhiZ=0;
  
    if(HepRotMat.num_col() == 0 || HepRotMat.num_row() == 0){
      cout << "No Rotation Matrics" << endl;
      HepRotMat(1,1) = 1;
      HepRotMat(2,2) = 1;
      HepRotMat(3,3) = 1;
    }    

  ThetaX = TMath::ACos(HepRotMat[0][2])/ TMath::Pi() * 180;
  ThetaY = TMath::ACos(HepRotMat[1][2])/ TMath::Pi() * 180;
  ThetaZ = TMath::ACos(HepRotMat[2][2])/ TMath::Pi() * 180;

  if(HepRotMat[0][0] > 0) 
    PhiX = TMath::ATan(HepRotMat[0][1]/HepRotMat[0][0])/ TMath::Pi() * 180;
  if(HepRotMat[0][0] < 0) 
    PhiX = TMath::ATan(HepRotMat[0][1]/HepRotMat[0][0])/ TMath::Pi() * 180 + 180;
  if(HepRotMat[0][0] == 0 && HepRotMat[0][1]>0 ) 
    PhiX = TMath::ATan(HepRotMat[0][1]/HepRotMat[0][0]);
  if(HepRotMat[0][0] == 0 && HepRotMat[0][1]<0) 
    PhiX = TMath::ATan(HepRotMat[0][1]/HepRotMat[0][0])/ TMath::Pi() * 180 + 180;

  if(HepRotMat[1][0] > 0) 
    PhiY = TMath::ATan(HepRotMat[1][1]/HepRotMat[1][0])/ TMath::Pi() * 180;
  if(HepRotMat[1][0] < 0) 
    PhiY = TMath::ATan(HepRotMat[1][1]/HepRotMat[1][0])/ TMath::Pi() * 180 + 180;
  if(HepRotMat[1][0] == 0 && HepRotMat[1][1]>0 ) 
    PhiY = TMath::ATan(HepRotMat[1][1]/HepRotMat[1][0])/ TMath::Pi() * 180;
  if(HepRotMat[1][0] == 0 && HepRotMat[1][1]<0) 
    PhiY = TMath::ATan(HepRotMat[1][1]/HepRotMat[1][0])/ TMath::Pi() * 180 + 180;
  
  if(HepRotMat[2][0] > 0) 
    PhiZ = TMath::ATan(HepRotMat[2][1]/HepRotMat[2][0])/ TMath::Pi() * 180;
  if(HepRotMat[2][0] < 0) 
    PhiZ = TMath::ATan(HepRotMat[2][1]/HepRotMat[2][0])/ TMath::Pi() * 180 + 180;
  if(HepRotMat[2][0] == 0 && HepRotMat[2][1]>0 ) 
    PhiZ = TMath::ATan(HepRotMat[2][1]/HepRotMat[2][0])/ TMath::Pi() * 180;
  if(HepRotMat[2][0] == 0 && HepRotMat[2][1]<0) 
    PhiZ = TMath::ATan(HepRotMat[2][1]/HepRotMat[2][0])/ TMath::Pi() * 180 + 180;
  return (new TRotMatrix("RotMat","RotMat",ThetaX,ThetaY,PhiY,PhiX,ThetaZ,PhiZ));
}
