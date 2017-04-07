/*!
  \file		CsMaterial.cc
  \brief	Compass geometry material class implementation file
  \author	$Author: tnagel $
  \version	$Revisopn$
  \date		$Date: 2010/01/28 12:51:25 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include <cmath>
#include <cstdlib>   // for abs()
#include <sstream>
#include "CsErrLog.h"

#include "CsMaterial.h"

using namespace std;

//---- CsMaterial ---------------------------------------------------------
CsMaterial::CsMaterial():
  id_(0), name_("NoName"), radiationLength_(0.0), density_(0.0){
}

CsMaterial::CsMaterial(	const unsigned int& materialId,
					  const string& name,
					  const double& density,
					  const double& radiationLength):
  id_(materialId), name_(name), 
  radiationLength_(radiationLength), density_(density){
}

CsMaterial::CsMaterial(const unsigned int& materialId, const string& name,
  const double& density, const int& nmat,
  vector<double>& a, vector<double>& z, vector<double>& w):
  id_(materialId), name_(name), 
  radiationLength_(0.0), density_(density){
  this->radiationLength( this->rLength( nmat, a, z, w) );
}

CsMaterial::CsMaterial(const CsMaterial& A):
  id_(A.id()), name_(A.name()), 
  radiationLength_(A.radiationLength()), density_(A.density()){
}

CsMaterial::~CsMaterial()
{
}

CsMaterial& CsMaterial::operator=(const CsMaterial& A){
  if(this != &A){
	this->id(A.id());
	this->name(A.name());
	this->radiationLength(A.radiationLength());
	this->density( A.density() );
  }
  return *this;
}

double CsMaterial::coulombCorrection( const double& z ) const {
	const double alpha	= 1.0 / 137.0359895; // fine structure function

	double aZ = ( alpha * z );
	double aZ2 = aZ*aZ;
	double aZ4 = aZ2*aZ2;
	double aZ6 = aZ2*aZ4;
	
	return aZ2*( 1/(1+aZ2) + 0.20206 - 0.0369*aZ2 + 0.0083*aZ4 - 0.0020*aZ6 );
}


double CsMaterial::rLength(	
  const int& nmat, vector<double>& a, vector<double>& z, vector<double>& w) const {

	const double alpha	= 1.0 / 137.0359895; // fine structure function
	const double r_e	= 2.81794092e-13;	// classical electron radius (cm)
	const double n_avo	= 6.0221367e23;

	const int npara = abs(nmat);
	
	vector<double> propotionByWeight( w );

	if( nmat < 0.0 ){
		double A_mol(0.0);
		for( int i=0; i<npara; i++ ) A_mol += a[i]*w[i];

		if( A_mol == 0 ){
		        ostringstream out;
			out << "ERROR in A and Weight arrays. Sum of A is equal to zero." << endl;
			for( int i=0; i<npara; i++ ){
				out << setw(2)	<< i ;
				out	<< ":\t(A:"	<< a[i] 
					<< "\tW:"		<< w[i]
					<< " )"		<< endl;
			}
			CsErrLog::Instance()->mes(elError, out.str());
			return 0.0;
		}	
	
		for( int i=0; i<npara; i++ ) 
		propotionByWeight[i]   = a[i]*w[i]/A_mol ;
	}

	vector<double> p_inv_dX0(propotionByWeight);
	

	for( int i=0; i<npara; i++ ){

		double Z = z[i];
		double logZ = log(Z);

		double C2 = log( 1440.0 ) - 2.0 * logZ / 3.0 ;
		double C3 = log(  183.0 ) - logZ / 3.0 ;

		double cCorr =	this->coulombCorrection( Z );
		double kusai =	C2 / ( C3 - cCorr );

		p_inv_dX0[i] *= 4.0 * alpha * r_e * r_e * n_avo * Z * (Z + kusai);
		p_inv_dX0[i] *= (C3  - cCorr) / a[i];

	}

	double iX(0.0);
	for( int i=0; i<npara; i++ ) iX += p_inv_dX0[i] * this->density() ;

	if( iX != 0 ) return ( 1.0 / iX );
	
	CsErrLog::Instance()->mes(elError, "Inverse of radiation length is equal to zero.");

	return 0.0;
}



istream& operator>>(istream& stream, CsMaterial& material){
  unsigned int        id; 
  string            name;
  double radiationLength;

  stream	>> id >> name >> radiationLength;

  material.id(id);
  material.name(name);
  material.radiationLength(radiationLength);

  return stream;	
}


ostream& operator<<(ostream& stream, const CsMaterial& material){
  stream
	<< "ID:\t"		<< material.id()		<< "\t"
	<< "NAME:\t"		<< material.name() 		<< "\t"
	<< "Density:\t"	<< material.density() 	<< "\t"
	<< "RadLen:\t"		<< material.radiationLength();
  return stream;
}
