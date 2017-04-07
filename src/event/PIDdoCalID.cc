#include"CsEvent.h" 
#include"CsErrLog.h"
#include"CsMCUtils.h"
#include<string>
#include<stdio.h>

#include "Reco/CalorimeterParticle.h"


void PID_doCalID( vector<CsParticle*>& parts ) {

  //***** CUTS & OPTIONS *****
  const int NSigm = 2;

  vector<CsParticle*>::iterator ip;
  map<const CsTrack*,CsParticle*> tp;
  map<Reco::CalorimeterParticle*,CsParticle*> cp;

  for( ip=parts.begin(); ip!=parts.end(); ip++ ) {
    
    const CsTrack* trk = (*ip)->getTrack();
    if( trk != NULL ) tp[ trk ] = (*ip);

    vector<Reco::CalorimeterParticle*> cvp = (*ip)->getCalObjects();
    if( cvp.size() == 1 ) cp[ cvp[0] ] = (*ip);
    
  }

  if( cp.size() == 0 || tp.size() == 0 ) return;
  
  map<const CsTrack*,CsParticle*>::iterator itp;
  map<Reco::CalorimeterParticle*,CsParticle*>::iterator icp;
  
  for( itp=tp.begin(); itp!=tp.end(); itp++ ) {

    const CsTrack* trk = (*itp).first;
    vector<CsHelix>vh = (*itp).first->getHelices();
    if( vh.size() < 2 ) continue;
    CsHelix H1,H2;
    H1 = vh[0];
    H2 = vh[1];
    
    for( icp=cp.begin(); icp!=cp.end(); icp++ ) {
      
      Reco::CalorimeterParticle* cpart = (*icp).first;
      double xCal,yCal,zCal;
      xCal = cpart->GetX();
      yCal = cpart->GetY();
      zCal = cpart->GetZ();

      double xECal,yECal;
      xECal = cpart->GetXerr();
      yECal = cpart->GetYerr();
      
      CsHelix Hc;
      if( fabs( H1.getZ() - zCal ) < fabs( H2.getZ() - zCal ) )
	H1.Extrapolate( zCal, Hc );
      else
	H2.Extrapolate( zCal, Hc );
      
      double xTrk,yTrk;
      xTrk = Hc.getX();
      yTrk = Hc.getY();

      if( ( xTrk - xCal ) < NSigm * xECal &&  
	  ( yTrk - yCal ) < NSigm * yECal ) {
      	
	tp[ trk ] -> addCalobj( cpart );
	
      }

    }

  }

  return;
}

