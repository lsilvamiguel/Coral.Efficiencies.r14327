#ifndef DigitizerSADCF_h
#define DigitizerSADCF_h

#include "CsDigitizerSADC.h"

class CsCalorimeter;
class TH1D;

////////////////////////////////////////////////////////////////////////////////

class DigitizerSADCF : public DigitizerSADCBase {
	public:
		virtual ~DigitizerSADCF();

		DigitizerSADCF ( CsCalorimeter* parent, Type type , int ic);

		virtual bool FitBase(  const std::vector<uint16> &sample );
		virtual bool Check( void ) const;

	private:
		const double newtonPrecission_;
		const double newtonNIter_;
		const int numberFourierComponents_;
		const int integrationWindow_;
		double parabolaShape_;
		int sampleSize_;
		bool writeTree_;

		std::vector<double> hDer_;
		std::vector<double> hPhase_;
		std::vector<double> hMag_;
// 		TH1D* hEvent_;
// 		TH1D* hMag_;
// 		TH1D* hPhase_;

		void FourierTransform();
		void GetRootNewton(double& root);
		double IntegrateDerivative(const int jmax, const double &dmax);
		void TreatOverflow(const std::vector<uint16> &sample);
		void ManageTree(const double &ampfourier, const double &timefourier);

};

////////////////////////////////////////////////////////////////////////////////

#endif  //         DigitizerSADCF_h
