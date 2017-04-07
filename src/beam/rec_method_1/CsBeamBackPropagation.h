#ifndef CsBeamBackPropagation_h 
#define CsBeamBackPropagation_h

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

class CsBMStrack;
class CsBMScluster;

class CsBeamBackPropagation{

   private:

      //pointer
      static CsBeamBackPropagation* instance_;

      //constants from file
      std::string _planeName[6];

      double _spaceResolutionBMS[6];
      double _correctionFactor[6];

      bool _readFromAligFile;

      //constants from functions 
      CsDetector* _planeReference[6];

      //extrapolated positions with errors
      double _x[6];
      double _x_err[6][3];
      double _dx[6];

      double _y[6];
      double _y_err[6][4];
      double _dy[6];

      //extrapolate result to cluster
      double extrapolate_to_cluster(int ipl, CsBMScluster* cluster);

      //histograms
      int _his_mode;
      void set_histograms();

      CsHist2D *his_backprop_01[6];                 //residuals for final BMS tracks
      CsHist2D *his_backprop_02[6];                 //residuals for final BMS tracks vs y, all cases  
      
   public:

      //constructor and desctructor
      CsBeamBackPropagation();
      ~CsBeamBackPropagation();

      //clear class for next event
      void clear();

      //extrapolate
      bool make_back_propagation(CsHelix* BTHelix, double mom);

      //get extrapolated positions
      double get_x(int i){    return _x[i]; }
      double* get_x_err(int i){  return _x_err[i]; }
      double get_dx(int i){   return _dx[i]; }

      double get_y(int i){    return _y[i]; }
      double* get_y_err(int i){  return _y_err[i]; }
      double get_dy(int i){   return _dy[i]; }

      //get coefficients
      double get_correctionFactor(int i){ return _correctionFactor[i];  }
      double get_spaceResolutionBMS(int i){ return _spaceResolutionBMS[i];  }

      //calculate agreement
      double calculate_chi2(CsBMStrack* BMStrack);

      //histogram residuals
      void histogram_residuals(CsHelix* BTHelix, double mom, CsBMStrack* BMStrack);

      //instance
      static CsBeamBackPropagation* Instance();
};

#endif
