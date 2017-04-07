//
// Misc. constants in global scope
//
#ifndef TConstants_h
#define TConstants_h


#include <math.h>

// Numerical arlgorithms' constants

const double TConstants_RKuttaMinStep = 0.0100;  // Ringe Kutta minimal step (100 mic)

// Constants for drawing

const float TConstants_SM1Siz[3]={70, 100, 25};  // 1-st magnet 1/2 sizes
const float TConstants_SM2Siz[3]={200, 100, 50}; // 2-d magnet 1/2 sizes

const float TConstants_A4=210./296.;       //A4 aspect ratio

// Misc
const int TConstants_IDet_max   = 999;     // Detector ID maximum
#ifdef HIGH_NTtrack_MAX
const int TConstants_NTtrack_max = 2000;   // max. number of space track pieces in detector group
#else
const int TConstants_NTtrack_max = 1000;   // max. number of space track pieces in detector group
#endif

//TraFFiC version
const float TConstants_TraFFiC_Version = 1.8;


#endif
