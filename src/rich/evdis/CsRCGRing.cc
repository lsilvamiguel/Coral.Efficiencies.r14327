// $Id:

/*!
   \file    CsRCGRing.cc
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#include "coral_config.h"
#include "CsRCGRing.h"

using namespace std;

ClassImp(CsRCGRing)
CsRCGRing::CsRCGRing(Int_t n, Double_t *ptr , Int_t mcolor, Int_t mtype) : TPolyLine3D(){

  SetLineColor(mcolor);
  SetLineStyle(mtype);
  SetPolyLine( n , ptr );

}
CsRCGRing::~CsRCGRing(){}

int CsRCGRing::GetType(){
   return Type;
}

void CsRCGRing::SetType(int type){
   Type = type;
}

void CsRCGRing::GetInfo(){
} 

void CsRCGRing::ExecuteEvent(Int_t event, Int_t px, Int_t py){
   if (event == kButton1Down){
     GetInfo();
   }
}
