#include "CsBeamMomCalculation.h"


//=== CONSTRUCTOR ====================================================================
CsBeamMomCalculation :: CsBeamMomCalculation(){

   //read configuration file 
   read_configuration_key("beam", "histogram_level", _his_mode);

   read_configuration_key("beam_momrec", "mesured_point_U", _par_zU);
   read_configuration_key("beam_momrec", "mesured_point_D", _par_zD);

   //set pointers to constants
   //BM01
   _pca_gNVariables[0]    = &_pca_BM01_gNVariables;
   _pca_gEigenVectors[0]  = _pca_BM01_gEigenVectors;
   _pca_gEigenValues[0]   = _pca_BM01_gEigenValues;
   _pca_gMeanValues[0]    = _pca_BM01_gMeanValues;
   _pca_gSigmaValues[0]   = _pca_BM01_gSigmaValues;

   _fit_gNVariables[0]    = &_fit_BM01_gNVariables;
   _fit_gNCoefficients[0] = &_fit_BM01_gNCoefficients;
   _fit_gDMean[0]         = &_fit_BM01_gDMean;
   _fit_gXMean[0]         = _fit_BM01_gXMean;
   _fit_gXMin[0]          = _fit_BM01_gXMin;
   _fit_gXMax[0]          = _fit_BM01_gXMax;
   _fit_gCoefficient[0]   = _fit_BM01_gCoefficient;
   _fit_gCoefficientRMS[0]= _fit_BM01_gCoefficientRMS;
   _fit_gPower[0]         = _fit_BM01_gPower;

   //BM05
   _pca_gNVariables[1]    = &_pca_BM05_gNVariables;
   _pca_gEigenVectors[1]  = _pca_BM05_gEigenVectors;
   _pca_gEigenValues[1]   = _pca_BM05_gEigenValues;
   _pca_gMeanValues[1]    = _pca_BM05_gMeanValues;
   _pca_gSigmaValues[1]   = _pca_BM05_gSigmaValues;

   _fit_gNVariables[1]    = &_fit_BM05_gNVariables;
   _fit_gNCoefficients[1] = &_fit_BM05_gNCoefficients;
   _fit_gDMean[1]         = &_fit_BM05_gDMean;
   _fit_gXMean[1]         = _fit_BM05_gXMean;
   _fit_gXMin[1]          = _fit_BM05_gXMin;
   _fit_gXMax[1]          = _fit_BM05_gXMax;
   _fit_gCoefficient[1]   = _fit_BM05_gCoefficient;
   _fit_gCoefficientRMS[1]= _fit_BM05_gCoefficientRMS;
   _fit_gPower[1]         = _fit_BM05_gPower;

   //BM02
   _pca_gNVariables[2]    = &_pca_BM02_gNVariables;
   _pca_gEigenVectors[2]  = _pca_BM02_gEigenVectors;
   _pca_gEigenValues[2]   = _pca_BM02_gEigenValues;
   _pca_gMeanValues[2]    = _pca_BM02_gMeanValues;
   _pca_gSigmaValues[2]   = _pca_BM02_gSigmaValues;

   _fit_gNVariables[2]    = &_fit_BM02_gNVariables;
   _fit_gNCoefficients[2] = &_fit_BM02_gNCoefficients;
   _fit_gDMean[2]         = &_fit_BM02_gDMean;
   _fit_gXMean[2]         = _fit_BM02_gXMean;
   _fit_gXMin[2]          = _fit_BM02_gXMin;
   _fit_gXMax[2]          = _fit_BM02_gXMax;
   _fit_gCoefficient[2]   = _fit_BM02_gCoefficient;
   _fit_gCoefficientRMS[2]= _fit_BM02_gCoefficientRMS;
   _fit_gPower[2]         = _fit_BM02_gPower;

   //BM03
   _pca_gNVariables[3]    = &_pca_BM03_gNVariables;
   _pca_gEigenVectors[3]  = _pca_BM03_gEigenVectors;
   _pca_gEigenValues[3]   = _pca_BM03_gEigenValues;
   _pca_gMeanValues[3]    = _pca_BM03_gMeanValues;
   _pca_gSigmaValues[3]   = _pca_BM03_gSigmaValues;

   _fit_gNVariables[3]    = &_fit_BM03_gNVariables;
   _fit_gNCoefficients[3] = &_fit_BM03_gNCoefficients;
   _fit_gDMean[3]         = &_fit_BM03_gDMean;
   _fit_gXMean[3]         = _fit_BM03_gXMean;
   _fit_gXMin[3]          = _fit_BM03_gXMin;
   _fit_gXMax[3]          = _fit_BM03_gXMax;
   _fit_gCoefficient[3]   = _fit_BM03_gCoefficient;
   _fit_gCoefficientRMS[3]= _fit_BM03_gCoefficientRMS;
   _fit_gPower[3]         = _fit_BM03_gPower;

   //BM06
   _pca_gNVariables[4]    = &_pca_BM06_gNVariables;
   _pca_gEigenVectors[4]  = _pca_BM06_gEigenVectors;
   _pca_gEigenValues[4]   = _pca_BM06_gEigenValues;
   _pca_gMeanValues[4]    = _pca_BM06_gMeanValues;
   _pca_gSigmaValues[4]   = _pca_BM06_gSigmaValues;

   _fit_gNVariables[4]    = &_fit_BM06_gNVariables;
   _fit_gNCoefficients[4] = &_fit_BM06_gNCoefficients;
   _fit_gDMean[4]         = &_fit_BM06_gDMean;
   _fit_gXMean[4]         = _fit_BM06_gXMean;
   _fit_gXMin[4]          = _fit_BM06_gXMin;
   _fit_gXMax[4]          = _fit_BM06_gXMax;
   _fit_gCoefficient[4]   = _fit_BM06_gCoefficient;
   _fit_gCoefficientRMS[4]= _fit_BM06_gCoefficientRMS;
   _fit_gPower[4]         = _fit_BM06_gPower;

   //BM04
   _pca_gNVariables[5]    = &_pca_BM04_gNVariables;
   _pca_gEigenVectors[5]  = _pca_BM04_gEigenVectors;
   _pca_gEigenValues[5]   = _pca_BM04_gEigenValues;
   _pca_gMeanValues[5]    = _pca_BM04_gMeanValues;
   _pca_gSigmaValues[5]   = _pca_BM04_gSigmaValues;

   _fit_gNVariables[5]    = &_fit_BM04_gNVariables;
   _fit_gNCoefficients[5] = &_fit_BM04_gNCoefficients;
   _fit_gDMean[5]         = &_fit_BM04_gDMean;
   _fit_gXMean[5]         = _fit_BM04_gXMean;
   _fit_gXMin[5]          = _fit_BM04_gXMin;
   _fit_gXMax[5]          = _fit_BM04_gXMax;
   _fit_gCoefficient[5]   = _fit_BM04_gCoefficient;
   _fit_gCoefficientRMS[5]= _fit_BM04_gCoefficientRMS;
   _fit_gPower[5]         = _fit_BM04_gPower;

   //full
   _pca_gNVariables[6]    = &_pca_full_gNVariables;
   _pca_gEigenVectors[6]  = _pca_full_gEigenVectors;
   _pca_gEigenValues[6]   = _pca_full_gEigenValues;
   _pca_gMeanValues[6]    = _pca_full_gMeanValues;
   _pca_gSigmaValues[6]   = _pca_full_gSigmaValues;

   _fit_gNVariables[6]    = &_fit_full_gNVariables;
   _fit_gNCoefficients[6] = &_fit_full_gNCoefficients;
   _fit_gDMean[6]         = &_fit_full_gDMean;
   _fit_gXMean[6]         = _fit_full_gXMean;
   _fit_gXMin[6]          = _fit_full_gXMin;
   _fit_gXMax[6]          = _fit_full_gXMax;
   _fit_gCoefficient[6]   = _fit_full_gCoefficient;
   _fit_gCoefficientRMS[6]= _fit_full_gCoefficientRMS;
   _fit_gPower[6]         = _fit_full_gPower;


   //set histograms
   //set_histograms();
}

//=== DESTRUCTOR =====================================================================
CsBeamMomCalculation :: ~CsBeamMomCalculation(){
}

//=== POSITION 2 POSITION IN NEW BASE ================================================
void CsBeamMomCalculation :: pca_x2p(double *x, double *p, int c) {
  for (int i = 0; i < (*_pca_gNVariables[c]); i++) {
    p[i] = 0;
    for (int j = 0; j < (*_pca_gNVariables[c]); j++)
      p[i] += (x[j] - _pca_gMeanValues[c][j]) 
        * _pca_gEigenVectors[c][j *  (*_pca_gNVariables[c]) + i] / _pca_gSigmaValues[c][j];

  }
}

//=== POSITION IN NEW BASE 2 MOM =====================================================
double CsBeamMomCalculation :: fit_p2mom(double *x, int c) {
  double returnValue = (*_fit_gDMean[c]);
  int    i = 0, j = 0, k = 0;
  for (i = 0; i < (*_fit_gNCoefficients[c]) ; i++) {
    double term = _fit_gCoefficient[c][i];
    for (j = 0; j < (*_fit_gNVariables[c]); j++) {
      int power = _fit_gPower[c][(*_fit_gNVariables[c]) * i + j]; 
      double p1 = 1, p2 = 0, p3 = 0, r = 0;
      double v =  1 + 2. / (_fit_gXMax[c][j] - _fit_gXMin[c][j]) * (x[j] -_fit_gXMax[c][j]);
      switch(power) {
      case 1: r = 1; break; 
      case 2: r = v; break; 
      default: 
        p2 = v; 
        for (k = 3; k <= power; k++) { 
          p3 = p2 * v;
          p1 = p2; p2 = p3; 
        }
        r = p3;
      }
      term *= r; 
    }
    returnValue += term;
  }
  return returnValue;
}

void CsBeamMomCalculation :: calculate_momentum(CsBMStrack *track, double *p, double *dp){

   //few veriables
   double S, Sz, Szz, Sy, Szy, error2, a, b;

   int i, j;

   int i_segment;
   bool not_empty[2][3];

   double z[2][3];
   double y[2][3];
   double dy[2][3];

   double slope[2];
   double y_max[2];

   int situation[2];

   double x1[4], x2[4];

   //get clusters
   std::vector<CsBMScluster*> clusters = track->get_clusters(); 

   //for upstream and downstream segments
   y_max[0] = y_max[1] = 1.E6;

   for(i_segment = 0; i_segment < 2; i_segment++){

      not_empty[i_segment][0] = false;
      not_empty[i_segment][1] = false;
      not_empty[i_segment][2] = false;

      z[i_segment][0] = 0.;
      z[i_segment][1] = 0.;
      z[i_segment][2] = 0.;

      y[i_segment][0] = 0.;
      y[i_segment][1] = 0.;
      y[i_segment][2] = 0.;

      dy[i_segment][0] = 0.;
      dy[i_segment][1] = 0.;
      dy[i_segment][2] = 0.;

      //downstream
      if(i_segment == 0){
         if( !(clusters[0]->isEmpty()) ){    //BM01
            z[i_segment][0]  = clusters[0]->get_z(); 
            y[i_segment][0]  = clusters[0]->get_y(); 
            dy[i_segment][0] = clusters[0]->get_dy(); 
            not_empty[i_segment][0] = true;
         }
         if( !(clusters[4]->isEmpty()) ){    //BM05
            z[i_segment][1]  = clusters[4]->get_z(); 
            y[i_segment][1]  = clusters[4]->get_y(); 
            dy[i_segment][1] = clusters[4]->get_dy(); 
            not_empty[i_segment][1] = true;
         }
         if( !(clusters[1]->isEmpty()) ){    //BM02
            z[i_segment][2]  = clusters[1]->get_z(); 
            y[i_segment][2]  = clusters[1]->get_y(); 
            dy[i_segment][2] = clusters[1]->get_dy(); 
            not_empty[i_segment][2] = true;
         }
      }

      //upstream
      if(i_segment == 1){
         if( !(clusters[2]->isEmpty()) ){    //BM03
            z[i_segment][0]  = clusters[2]->get_z(); 
            y[i_segment][0]  = clusters[2]->get_y(); 
            dy[i_segment][0] = clusters[2]->get_dy(); 
            not_empty[i_segment][0] = true;
         }
         if( !(clusters[5]->isEmpty()) ){    //BM06
            z[i_segment][1]  = clusters[5]->get_z(); 
            y[i_segment][1]  = clusters[5]->get_y(); 
            dy[i_segment][1] = clusters[5]->get_dy(); 
            not_empty[i_segment][1] = true;
         }
         if( !(clusters[3]->isEmpty()) ){    //BM04
            z[i_segment][2]  = clusters[3]->get_z(); 
            y[i_segment][2]  = clusters[3]->get_y(); 
            dy[i_segment][2] = clusters[3]->get_dy(); 
            not_empty[i_segment][2] = true;
         }
      }

      //which situation we have
      situation[i_segment] = not_empty[i_segment][0] + not_empty[i_segment][1] + not_empty[i_segment][2];

      //calculate slope and y position if possible
      if(situation[i_segment] == 3){

         S = Sz = Szz = Sy = Szy = 0.;

         for(int i_plane = 0; i_plane < 3; i_plane++){

            error2 = pow(dy[i_segment][i_plane], 2);
            S     += 1./error2;
            Sz    += z[i_segment][i_plane]/error2;
            Szz   += pow(z[i_segment][i_plane], 2)/error2;
            Sy    += y[i_segment][i_plane]/error2;
            Szy   += (z[i_segment][i_plane]*y[i_segment][i_plane])/error2;

         }

         a = (S*Szy - Sz*Sy)/(S*Szz - Sz*Sz);
         b = (Szz*Sy - Sz*Szy)/(S*Szz - Sz*Sz);

         //extrapolate to measured points 
         slope[i_segment] = atan(a);
         if(i_segment == 0) y_max[i_segment] = _par_zU*a + b;
         if(i_segment == 1) y_max[i_segment] = _par_zD*a + b;
      }
      
      else if(situation[i_segment] == 2){

         i = j = -1;

         if(not_empty[i_segment][0] && not_empty[i_segment][1]){  i = 0;   j = 1;}
         if(not_empty[i_segment][0] && not_empty[i_segment][2]){  i = 0;   j = 2;}
         if(not_empty[i_segment][1] && not_empty[i_segment][2]){  i = 1;   j = 2;}

         slope[i_segment] = atan((y[i_segment][j]-y[i_segment][i])/(z[i_segment][j]-z[i_segment][i]));

         //extrapolate to measured points 
         if(i_segment == 0){
            y_max[i_segment] = y[i_segment][i] + (_par_zU-z[i_segment][i])*tan(slope[i_segment]); 
         }

         if(i_segment == 1){
            y_max[i_segment] = y[i_segment][i] + (_par_zD-z[i_segment][i])*tan(slope[i_segment]); 
         }
      }

   }

   //calculate momentum
   if(situation[0] > 1 && situation[1] > 1){
      x1[0] = y_max[0];
      x1[1] = 1000*slope[0];
      x1[2] = y_max[1];
      x1[3] = 1000*slope[1];

      pca_x2p(x1, x2, 6);
      *p = 1./fit_p2mom(x2, 6);
   }
   else if(situation[0] == 1 && situation[1] > 1){

      i = -1;

      if(not_empty[0][0]){ i = 0;   }
      if(not_empty[0][1]){ i = 1;   }
      if(not_empty[0][2]){ i = 2;   }

      x1[0] = y[0][i]; 
      x1[1] = y_max[1]; 
      x1[2] = 1000*slope[1];
      x1[3] = 0.; 

      pca_x2p(x1, x2, i);
      *p = 1./fit_p2mom(x2, i);
   }
   else if(situation[0] > 1 && situation[1] == 1){

      i = -1;

      if(not_empty[1][0]){ i = 0;   }
      if(not_empty[1][1]){ i = 1;   }
      if(not_empty[1][2]){ i = 2;   }

      x1[0] = y_max[0]; 
      x1[1] = 1000*slope[0];
      x1[2] = y[1][i]; 
      x1[3] = 0.; 

      pca_x2p(x1, x2, 3+i);
      *p = 1./fit_p2mom(x2, 3+i);
   }
   else{
      *p = 160.;
   }

   return;
}


//=== INSTANCE =======================================================================
CsBeamMomCalculation* CsBeamMomCalculation :: instance_ = NULL;

CsBeamMomCalculation* CsBeamMomCalculation :: Instance() {
   if( instance_ == NULL ){
      instance_ = new CsBeamMomCalculation();
   }
       
   return( instance_ );
}

