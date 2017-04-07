/*!
  Member function of class Helix
  to add "noise" (mult. scatt.) on the plane ipl
*/

#include "THlx.h"
#include "TSetup.h"
#include <iostream>
#include <math.h>


// add multiple scattering to cov. matrix
// for matherial with radiation length RadLen and thickness x

void THlx::AddNoise(float x, float RadLen)

{
  if(Hpar[5] == 0.) return; // momentum not known. Do nothing.
  
  double len   =  x / RadLen;

  // Lynch and Dahl aproximation for Sigma(Theta_proj) of mult. scatt.
  double SigTheta = 0.0136*fabs(Hpar[5]) * sqrt(len) * (1.+0.038*log(len));
 
  // Noise matrix calculation (NIM A329 (1993) 493-500)
  // Transverse displacement of the track is ignored.
  double p3 = Hpar[3];
  double p4 = Hpar[4];

  double p3p3 = SigTheta*SigTheta * (1 + p3*p3) * (1 + p3*p3 + p4*p4);
  double p4p4 = SigTheta*SigTheta * (1 + p4*p4) * (1 + p3*p3 + p4*p4);
  double p3p4 = SigTheta*SigTheta * p3*p4       * (1 + p3*p3 + p4*p4);

  Hcov[5] += p3p3;
  Hcov[8] += p3p4;
  Hcov[9] += p4p4;

  return;
}

// detector's plane contribution to mult. scatt.

void THlx::AddNoise(int ipl)
{
  const TDetect& d = TSetup::Ref().iPlane2Detect(ipl);
  double x = d.Siz(0)/this->DirCos(1);
  return(AddNoise(x,d.RadLen));
}








