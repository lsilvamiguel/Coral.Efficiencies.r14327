// $Id: THlxExtrap.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!

 Member function of class Helix
 to extrapolate it to plane, perpendicular to the beam
 at position "x".

 Resulting helix is stored in "Hout".

 Numerical Helix extrapolation with simultaneous
 Jacobian calculation in Runge-Kutta propagation function
 is used.
 (thanks to Igor.Gavrilenko@cern.ch)

*/

#include <iostream>
#include <iomanip>
#include "TOpt.h"
#include "TAlgo.h"
#include "TConstants.h"
#include "THlx.h"

using namespace std;

bool THlx::Extrap(double xlast, THlx& Hout) const
{

  if(TOpt::ReMode[10] > 0) { // use old multy-pass helix propagation
    Hout(0)=xlast;
    return(this->NpassExtrap(Hout));
  }

  
  double xfirst = this->Hpar[0];
  double s(0);

  if(fabs(xlast-xfirst) < TConstants_RKuttaMinStep) { // very short step
    Hout = (*this);
    Hout.Hpar[0] = xlast;
    Hout.path=fabs(xlast-xfirst);
    Hout.radLenFr = 0.;
    Hout.eLoss = 0;
    return(true);
  }

  //
  // Undefined momentum. Do straight line extrapolation.
  //

  if(Hpar[5] == 0.){
    Hout.Hpar[0]=xlast;
    double dx=Hout.Hpar[0]-Hpar[0];
    Hout.Hpar[1]=Hpar[1]+Hpar[3]*dx;
    Hout.Hpar[2]=Hpar[2]+Hpar[4]*dx;
    Hout.Hpar[3]=Hpar[3];
    Hout.Hpar[4]=Hpar[4];
    Hout.Hpar[5]=Hpar[5];
    
    // Cov matrix propagation
    Hout.Hcov[0] = Hcov[0] + Hcov[3]*dx + (Hcov[3] + Hcov[5]*dx)*dx;
    Hout.Hcov[1] = Hcov[1] + Hcov[6]*dx + (Hcov[4] + Hcov[8]*dx)*dx;
    Hout.Hcov[2] = Hcov[2] + Hcov[7]*dx + (Hcov[7] + Hcov[9]*dx)*dx;
    Hout.Hcov[3] = Hcov[3] + Hcov[5]*dx;
    Hout.Hcov[4] = Hcov[4] + Hcov[8]*dx;
    Hout.Hcov[5] = Hcov[5];
    Hout.Hcov[6] = Hcov[6] + Hcov[8]*dx;
    Hout.Hcov[7] = Hcov[7] + Hcov[9]*dx;
    Hout.Hcov[8] = Hcov[8];
    Hout.Hcov[9] = Hcov[9];
    Hout.Hcov[10]= Hcov[10]+ Hcov[12]*dx;
    Hout.Hcov[11]= Hcov[11]+ Hcov[13]*dx;
    Hout.Hcov[12]= Hcov[12];
    Hout.Hcov[13]= Hcov[13];
    Hout.Hcov[14]= Hcov[14];
    
    Hout.path = this->Dist(Hout);
    Hout.radLenFr = 0;
    Hout.eLoss = 0;
    return(true);
  }
    
  //
  // Runge-Kutta extrapolation
  //

  // cut on  P < 100 MeV
  if(fabs(1./Hpar[5]) < 0.100) {
    //cout<<"THlx::Extrap() ==> Too low momentum for Runge-Kutta propagation: "
    //<<fabs(1./Hpar[5]) <<" GeV"<<endl;
    return(false);
  }
  


  // Prepare surface description array
  double su[4];
  if(xlast < 0 ) su[0] = -1.;
  else           su[0] =  1.;
  su[1]=su[2]=0.;
  su[3] = fabs(xlast);

  // Prepare track parameters

  double x  = Hpar[0];
  double y  = Hpar[1];
  double z  = Hpar[2];
  double ax = 1./sqrt(Hpar[3]*Hpar[3] + Hpar[4]*Hpar[4] + 1);
  double ay = Hpar[3]*ax;
  double az = Hpar[4]*ax;
  double qP = Hpar[5];
  
  // Prepare Jacobian dP/dH where 
  //
  // P is track parameters vector (size 7), needed for RkutaXX routine,
  // in the direction cosines (ax,ay,az) representation: (x,y,x,ax,ay,az,q/P)  
  //
  // H is helix.
  //

  double P[42] = {
    x,      y,      z,         ax,          ay,         az,     qP,   // Parameters (P) 
    0,      1,      0,          0,           0,           0,     0,   // dP/dHpar[1]			  
    0,      0,      1,          0,           0,           0,     0,   // dP/dHpar[2]
    0,      0,      0,  -ay*ax*ax, ax-ay*ay*ax,   -ay*az*ax,     0,   // dP/dHpar[3]
    0,      0,      0,  -az*ax*ax,   -ay*az*ax, ax-az*az*ax,     0,   // dP/dHpar[4]
    0,      0,      0,          0,           0,           0,     qP   // dP/dHpar[5]*Hpar[5]
  };

  // Do the propagation
  if( !TAlgo::RkutaNoGrad(su,P,s) ) return(false);

  if(P[3] <= 0.) {
    cout<<"THlx::Extrap ==> Trajectory had turned back during extrapolation. P = "
	<<fabs(1./P[6])<<" GeV"<<endl;
    return(false);
  }

  // Output helix parameters

  Hout.Hpar[0] = P[0]; // x
  Hout.Hpar[1] = P[1]; // y
  Hout.Hpar[2] = P[2]; // z
  Hout.Hpar[3] = P[4]/P[3];    // yp = ay/ax
  Hout.Hpar[4] = P[5]/P[3];    // zp = az/ax
  Hout.Hpar[5] = P[6];


  // Calculate helix transformation Jacobian (F = dHout/dHin)  

   
  double p3Ax = -P[4]/(P[3]*P[3]);
  double p3Ay = 1./P[3];
  double p3Az = 0.;
  double p4Ax = -P[5]/(P[3]*P[3]);
  double p4Ay = 0.;
  double p4Az = p3Ay;
  
  
  double f[5][5];

  f[0][0]= P[ 8];    f[1][0]= P[ 9];
  f[0][1]= P[15];    f[1][1]= P[16];
  f[0][2]= P[22];    f[1][2]= P[23];
  f[0][3]= P[29];    f[1][3]= P[30];
  f[0][4]= P[36]/qP; f[1][4]= P[37]/qP;
  
  f[2][0]= p3Ax*P[10]+p3Ay*P[11]+p3Az*P[12];     f[3][0]= p4Ax*P[10]+p4Ay*P[11]+p4Az*P[12];   
  f[2][1]= p3Ax*P[17]+p3Ay*P[18]+p3Az*P[19];     f[3][1]= p4Ax*P[17]+p4Ay*P[18]+p4Az*P[19];
  f[2][2]= p3Ax*P[24]+p3Ay*P[25]+p3Az*P[26];     f[3][2]= p4Ax*P[24]+p4Ay*P[25]+p4Az*P[26];
  f[2][3]= p3Ax*P[31]+p3Ay*P[32]+p3Az*P[33];     f[3][3]= p4Ax*P[31]+p4Ay*P[32]+p4Az*P[33];
  f[2][4]=(p3Ax*P[38]+p3Ay*P[39]+p3Az*P[40])/qP; f[3][4]=(p4Ax*P[38]+p4Ay*P[39]+p4Az*P[40])/qP;
  
  f[4][0]=f[4][1]=f[4][2]=f[4][3]=0; f[4][4]=1;


  // Propagate Cov matrix (F*Cov*F.t)

  double w[5];

  w[0]=Hcov[0] * f[0][0] + Hcov[1] * f[0][1] + Hcov[3] * f[0][2] + Hcov[6] * f[0][3] + Hcov[10]*f[0][4];
  w[1]=Hcov[1] * f[0][0] + Hcov[2] * f[0][1] + Hcov[4] * f[0][2] + Hcov[7] * f[0][3] + Hcov[11]*f[0][4];
  w[2]=Hcov[3] * f[0][0] + Hcov[4] * f[0][1] + Hcov[5] * f[0][2] + Hcov[8] * f[0][3] + Hcov[12]*f[0][4];
  w[3]=Hcov[6] * f[0][0] + Hcov[7] * f[0][1] + Hcov[8] * f[0][2] + Hcov[9] * f[0][3] + Hcov[13]*f[0][4];
  w[4]=Hcov[10]* f[0][0] + Hcov[11]* f[0][1] + Hcov[12]* f[0][2] + Hcov[13]* f[0][3] + Hcov[14]*f[0][4];

  Hout.Hcov[0]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];

  w[0]=Hcov[0] * f[1][0] + Hcov[1] * f[1][1] + Hcov[3] * f[1][2] + Hcov[6] * f[1][3] + Hcov[10]*f[1][4];
  w[1]=Hcov[1] * f[1][0] + Hcov[2] * f[1][1] + Hcov[4] * f[1][2] + Hcov[7] * f[1][3] + Hcov[11]*f[1][4];
  w[2]=Hcov[3] * f[1][0] + Hcov[4] * f[1][1] + Hcov[5] * f[1][2] + Hcov[8] * f[1][3] + Hcov[12]*f[1][4];
  w[3]=Hcov[6] * f[1][0] + Hcov[7] * f[1][1] + Hcov[8] * f[1][2] + Hcov[9] * f[1][3] + Hcov[13]*f[1][4];
  w[4]=Hcov[10]* f[1][0] + Hcov[11]* f[1][1] + Hcov[12]* f[1][2] + Hcov[13]* f[1][3] + Hcov[14]*f[1][4];

  Hout.Hcov[1]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
  Hout.Hcov[2]=w[0]*f[1][0] + w[1]*f[1][1] + w[2]*f[1][2] + w[3]*f[1][3] + w[4]*f[1][4];

  w[0]=Hcov[0] * f[2][0] + Hcov[1] * f[2][1] + Hcov[3] * f[2][2] + Hcov[6] * f[2][3] + Hcov[10]*f[2][4];
  w[1]=Hcov[1] * f[2][0] + Hcov[2] * f[2][1] + Hcov[4] * f[2][2] + Hcov[7] * f[2][3] + Hcov[11]*f[2][4];
  w[2]=Hcov[3] * f[2][0] + Hcov[4] * f[2][1] + Hcov[5] * f[2][2] + Hcov[8] * f[2][3] + Hcov[12]*f[2][4];
  w[3]=Hcov[6] * f[2][0] + Hcov[7] * f[2][1] + Hcov[8] * f[2][2] + Hcov[9] * f[2][3] + Hcov[13]*f[2][4];
  w[4]=Hcov[10]* f[2][0] + Hcov[11]* f[2][1] + Hcov[12]* f[2][2] + Hcov[13]* f[2][3] + Hcov[14]*f[2][4];

  Hout.Hcov[3]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
  Hout.Hcov[4]=w[0]*f[1][0] + w[1]*f[1][1] + w[2]*f[1][2] + w[3]*f[1][3] + w[4]*f[1][4];
  Hout.Hcov[5]=w[0]*f[2][0] + w[1]*f[2][1] + w[2]*f[2][2] + w[3]*f[2][3] + w[4]*f[2][4];

  w[0]=Hcov[0] * f[3][0] + Hcov[1] * f[3][1] + Hcov[3] * f[3][2] + Hcov[6] * f[3][3] + Hcov[10]*f[3][4];
  w[1]=Hcov[1] * f[3][0] + Hcov[2] * f[3][1] + Hcov[4] * f[3][2] + Hcov[7] * f[3][3] + Hcov[11]*f[3][4];
  w[2]=Hcov[3] * f[3][0] + Hcov[4] * f[3][1] + Hcov[5] * f[3][2] + Hcov[8] * f[3][3] + Hcov[12]*f[3][4];
  w[3]=Hcov[6] * f[3][0] + Hcov[7] * f[3][1] + Hcov[8] * f[3][2] + Hcov[9] * f[3][3] + Hcov[13]*f[3][4];
  w[4]=Hcov[10]* f[3][0] + Hcov[11]* f[3][1] + Hcov[12]* f[3][2] + Hcov[13]* f[3][3] + Hcov[14]*f[3][4];

  Hout.Hcov[6]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
  Hout.Hcov[7]=w[0]*f[1][0] + w[1]*f[1][1] + w[2]*f[1][2] + w[3]*f[1][3] + w[4]*f[1][4];
  Hout.Hcov[8]=w[0]*f[2][0] + w[1]*f[2][1] + w[2]*f[2][2] + w[3]*f[2][3] + w[4]*f[2][4];
  Hout.Hcov[9]=w[0]*f[3][0] + w[1]*f[3][1] + w[2]*f[3][2] + w[3]*f[3][3] + w[4]*f[3][4];

    // as dPinv/dpar = 0 if par != Pinv
  Hout.Hcov[10]=Hcov[10]*f[0][0] + Hcov[11]*f[0][1] + Hcov[12]*f[0][2] + Hcov[13]*f[0][3] + Hcov[14]*f[0][4];
  Hout.Hcov[11]=Hcov[10]*f[1][0] + Hcov[11]*f[1][1] + Hcov[12]*f[1][2] + Hcov[13]*f[1][3] + Hcov[14]*f[1][4];
  Hout.Hcov[12]=Hcov[10]*f[2][0] + Hcov[11]*f[2][1] + Hcov[12]*f[2][2] + Hcov[13]*f[2][3] + Hcov[14]*f[2][4];
  Hout.Hcov[13]=Hcov[10]*f[3][0] + Hcov[11]*f[3][1] + Hcov[12]*f[3][2] + Hcov[13]*f[3][3] + Hcov[14]*f[3][4];

  Hout.Hcov[14]=Hcov[14];


  if(Hout.Hcov[0] < 0 || Hout.Hcov[2] < 0 ||
     Hout.Hcov[5] < 0 || Hout.Hcov[9] < 0 ||
     Hout.Hcov[14]< 0){
    cout<<"THlx::Extrapolate() ==> output cov. matrix is wrong "<<endl;
    this->Print("Input helix");
    Hout.Print("Output helix");
    return(false);
  }


  Hout.path     = s;
  Hout.radLenFr = 0;
  Hout.eLoss = 0;
  return(true);
}
