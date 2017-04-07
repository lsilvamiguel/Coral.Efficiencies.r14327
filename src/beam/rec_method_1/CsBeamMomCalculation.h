#ifndef CsBeamMomCalculation_h 
#define CsBeamMomCalculation_h

#include "CsTypes.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsErrLog.h"
#include "CsOpt.h"
#include "DaqDataDecoding/DaqEvent.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdio.h>

#include "CsBeamReconstructionSupport.h"

#include "beam_coefficients/_fit_160_full.h"
#include "beam_coefficients/_fit_160_BM01.h"
#include "beam_coefficients/_fit_160_BM02.h"
#include "beam_coefficients/_fit_160_BM03.h"
#include "beam_coefficients/_fit_160_BM04.h"
#include "beam_coefficients/_fit_160_BM05.h"
#include "beam_coefficients/_fit_160_BM06.h"

#include "beam_coefficients/_pca_160_full.h"
#include "beam_coefficients/_pca_160_BM01.h"
#include "beam_coefficients/_pca_160_BM02.h"
#include "beam_coefficients/_pca_160_BM03.h"
#include "beam_coefficients/_pca_160_BM04.h"
#include "beam_coefficients/_pca_160_BM05.h"
#include "beam_coefficients/_pca_160_BM06.h"

class CsBMStrack;

//BMS back propagation
class CsBeamMomCalculation{

   private:

      //pointer
      static CsBeamMomCalculation* instance_;

      //constants from file
      double _par_zU;                              //parametrization points
      double _par_zD;

      //pointers to constants
      const int *_pca_gNVariables[7];
      const double *_pca_gEigenVectors[7];
      const double *_pca_gEigenValues[7];
      const double *_pca_gMeanValues[7];
      const double *_pca_gSigmaValues[7];

      const int *_fit_gNVariables[7];
      const int *_fit_gNCoefficients[7];
      const double *_fit_gDMean[7];
      const double *_fit_gXMean[7];
      const double *_fit_gXMin[7];
      const double *_fit_gXMax[7];
      const double *_fit_gCoefficient[7];
      const double *_fit_gCoefficientRMS[7];
      const int *_fit_gPower[7];

      //calculate pca
      void pca_x2p(double *x, double *p, int c);
      
      //calculate momentum
      double fit_p2mom(double *p, int c); 

      //histograms
      int _his_mode;
      void set_histograms();
      
   public:

      //constructor and descructor
      CsBeamMomCalculation();
      ~CsBeamMomCalculation();

      //calculate momentum
      void calculate_momentum(CsBMStrack *track, double *p, double *dp);

      //instance
      static CsBeamMomCalculation* Instance();
};
#endif
