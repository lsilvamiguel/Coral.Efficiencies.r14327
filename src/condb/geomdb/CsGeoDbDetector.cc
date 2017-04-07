#if __GNUG__ >= 2
#  pragma implementation
#endif
/*!
  \file		CsGeoDbDetector.cc
  \brief	Geometry Detector table interface definition (implementation)
  \author	$Author: miyachi $
  \version	$Revision: 1.1 $
  \date		$Date: 2000/06/07 08:36:42 $
*/

#include "CsGeoDbDetector.h"

//---- CsGeoDbDetector ---------------------------------------------------------

CsGeoDbDetector::CsGeoDbDetector()
{
}

CsGeoDbDetector::~CsGeoDbDetector()
{
}

ostream& operator<<(ostream& os, CsGeoDbDetector& detector){

	string name = detector.name();

	int rotMId(1);

	os	<< setfill( ' ' )	
		<< setiosflags(ios::fixed)
		<< " det"
		<< setw(5)	<< detector.id()
		<< setw(6)	<< detector.name()
		<< setw(3)	<< detector.unit()
		<< setw(5)	<< detector.type()
		<< setprecision(2)	
		<< setw(9)	<< detector.radLen();

	
	vector<double> vDet = detector.param();

	os	<< setw(4)	<< vDet.size();
	for( unsigned int i=0; i < vDet.size(); i++){
		os	<< setprecision(3)	<< setw(8)	<< vDet[i];
	}
	
	os	<< setprecision(4)	
		<< setw(11)	<< detector.center().x()
		<< setw(11)	<< detector.center().y()
		<< setw(11)	<< detector.center().z()
		<< setw(7)	<< rotMId
		<< setw(14)	<< detector.wirD()
		<< setprecision(3)	
		<< setw(9)	<< detector.angle()
		<< setw(6)	<< detector.nWir()
		<< setprecision(4)	
		<< setw(9)	<< detector.wirP()
		<< setprecision(3)	
		<< setw(6)	<< detector.eff()
		<< setw(6)	<< detector.bgk()
		<< setprecision(1)	
		<< setw(7)	<< detector.tGate();
		
	if(detector.hasDrift())
	os	<< setprecision(4)	
		<< setw(11)	<< detector.driftV()
		<< setprecision(1)	
		<< setw(6)	<< detector.T0()
		<< setw(6)	<< detector.dHitRes()
		<< setw(7)	<< detector.sRes()
		<< setw(6)	<< detector.tSlic();


	return os;
}
