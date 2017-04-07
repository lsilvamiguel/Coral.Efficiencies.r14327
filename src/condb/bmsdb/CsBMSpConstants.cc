/*!
   \file    CsBMSpConstants.cc
   \brief   Compass BMS calibration constants persistent object implementation file.
   \author  Yoshiyuki Miyachi
   \version $Revision: 1.1 $
   \date    $Date: 2000/02/28 10:09:15 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsBMSpConstants.ddl"

//---- CsBMSpConstants ---------------------------------------------------------

CsBMSpConstants::CsBMSpConstants(const CsBMSconstants& tConst) {
	tConst.getBMSconst( header_, p0_, avsg16_, evec16_, coff16_);
}

CsBMSpConstants::operator CsBMSconstants () {
	return CsBMSconstants(header_, p0_, avsg16_, evec16_, coff16_);
}

void CsBMSpConstants::printBMSconst(ostream& os){
	os << header_ << endl;
	os << "Pbeam = " << p0_ << endl;
	os <<  " AVSG16:" << endl;

	int i, j, ll;
	for( i=0; i<4; i++){
		os	<< '\t'	<< avsg16_[i]; 
	}
	os << endl;

	os <<  " EVEC16:" << endl;
	for( i=0; i<4; i++){
		for( j=0; j<4;j++){
			os << '\t'	<< evec16_[j][i]; 
		}
		os << endl;
	}

	os <<  " COFF16:" << endl;
	for( i=0; i<5; i++){
		for( j=0; j<3;j++){
			for(  ll=0; ll<4 ;ll++){
				os << '\t'	<< coff16_[ll][j][i]; 
			}
			os << endl;
		}
	}
}

ostream& operator<<(ostream& os, const CsBMSpConstants& pConst){
	pConst.printBMSconst(os);
	return os;
}
