// $Id:

/*!
   \file    CsRCGRingTrack.cc
   \brief   CORAL Event Display Package.
   \version $Revision: 1.4 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#include "coral_config.h"
#include "CsRCGRingTrack.h"

using namespace std;

ClassImp(CsRCGRingTrack)
CsRCGRingTrack::CsRCGRingTrack() : TPolyMarker3D(){
}

CsRCGRingTrack::CsRCGRingTrack(int n,Float_t *x,Marker_t m,Int_t color, Float_t size) : TPolyMarker3D(){

  Option_t t[]="";
  SetPolyMarker(n,x,m,t);
  SetMarkerColor(color);
  SetMarkerSize(size);

}

CsRCGRingTrack::~CsRCGRingTrack(){}


void CsRCGRingTrack::GetInfo(){
  cout << " Particle : " <<  Name 
       << " Charge   : " << Charge 
       << endl;
}

void CsRCGRingTrack::SetCharge(int charge){
Charge = charge;
int color=7; 
 if(Charge == 0) color =  4;
 if(Charge == -1) color = 5;
 if(Charge == 1) color =  3;
 if(Charge >= 2 || Charge <= -2) color =  7;
 SetMarkerColor(color);
}

void CsRCGRingTrack::SetName(const char* name){
  strncpy(Name, name, 80); Name[79] = '\n';
}

void CsRCGRingTrack::ExecuteEvent(Int_t event, Int_t px, Int_t py){
   static Float_t x,y;
   static Float_t xmin,ymin;
   static Float_t xrange,yrange;
   static Float_t pxold, pyold;
   if (event == kButton1Down){
     GetInfo();
   }
}
