// $Id:

/*!
   \file    CsRCGHit.cc
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#include "coral_config.h"
#include "CsRCGHit.h"

using namespace std;

ClassImp(CsRCGHit)
  CsRCGHit::CsRCGHit(Int_t type,Double_t xcm,Double_t ycm,Double_t zcm,Double_t dxcm,Double_t dycm,Double_t dzcm,Int_t color) : TMarker3DBox(){

 Xcm  = xcm;
 Ycm  = ycm;
 Zcm  = zcm;
 DXcm = dxcm;
 DYcm = dycm;
 DZcm = dzcm;
 SetType(type);
 SetFillStyle(1001);
 SetFillColor(color);
 SetPosition(Xcm,Ycm,Zcm);
 SetSize(DXcm,DYcm,DZcm);
 SetLineColor(color);

}
CsRCGHit::~CsRCGHit(){}

void CsRCGHit::SetFactorX(float fact){
  factorX = fact;
  float x,y,z;
  GetPosition(x,y,z);
  SetPosition(x * factorX,y,z);
}

int CsRCGHit::GetType(){
  return Type;
}

void CsRCGHit::SetType(int type){
  Type = type;
}

void CsRCGHit::GetInfo(){
  cout << "Hit  Xcm:" << Xcm 
       << " Ycm:" << Ycm 
       << " Zcm:" << Zcm 
       << endl;  
} 

void CsRCGHit::ExecuteEvent(Int_t event, Int_t px, Int_t py){
   if (event == kButton1Down){
     GetInfo();
   }
}
