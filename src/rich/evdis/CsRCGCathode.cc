// $Id:

/*!
   \file    CsRCGHit.cc
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#include "coral_config.h"
#include "CsRCGCathode.h"

using namespace std;

ClassImp(CsRCGCathode)

 CsRCGCathode::CsRCGCathode(Int_t type,Double_t xcm,Double_t ycm,Double_t zcm,Double_t dxcm,Double_t dycm,Double_t dzcm, Int_t color) : TMarker3DBox(){

 Xcm  = xcm;
 Ycm  = ycm;
 Zcm  = zcm;
 DXcm = dxcm;
 DYcm = dycm;
 DZcm = dzcm;
 SetType(type);
 SetFillStyle(0);
 SetFillColor(color);
 SetPosition(Xcm,Ycm,Zcm);
 SetSize(DXcm,DYcm,DZcm);
 SetLineColor(color);

}
CsRCGCathode::~CsRCGCathode(){}

void CsRCGCathode::SetFactorX(float fact){
  factorX = fact;
  float x,y,z;
  GetPosition(x,y,z);
  SetPosition(x * factorX,y,z);
}

int CsRCGCathode::GetType(){
  return Type;
}

void CsRCGCathode::SetType(int type){
  Type = type;
}

void CsRCGCathode::GetInfo(){
  cout << "Hit  Xcm:" << Xcm 
       << " Ycm:" << Ycm 
       << " Zcm:" << Zcm 
       << endl;  
} 

void CsRCGCathode::ExecuteEvent(Int_t event, Int_t px, Int_t py){
   if (event == kButton1Down){
     GetInfo();
   }
}
