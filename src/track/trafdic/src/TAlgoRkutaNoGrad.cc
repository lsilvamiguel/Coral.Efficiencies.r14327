//                                                                                    //
//------------------------------------------------------------------------------------//
//                                                                                    //
// Runge-Kutta method for tracking a particles through a magnetic field.              //
// Uses Nystroem algorithm (See Handbook Net. Bur. of Standards, procedure 25.5.20    //
//                                                                                    //
// Input parameteres:                                                                 //
//    SU     - plane parameters                                                       //
//    SU [0] - direction cosines normal to surface Ex                                 //
//    SU [1] -          -------                    Ey                                 //
//    SU [2] -          -------                    Ez, Ex*Ex+Ey*Ey+Ez*Ez=1            //
//    SU [3] - distance to surface from (0,0,0) > 0 cm                                //
//    VO     - initial parameters (coordinates(cm), direction cosines,                //
//             charge/momentum (Gev-1) and derivatives this parameters                //
//                                                             Ax*Ax+Ay*Ay+Az*Az=1    //
//      X        Y        Z        Ax       Ay       Az       q/P                     //
//    VO[ 0]   VO[ 1]   VO[ 2]   VO[ 3]   VO[ 4]   VO[ 5]   VO[ 6]                    //
//    dX/dp    dY/dp    dZ/dp    dAx/dp   dAy/dp   dAz/dp   d(q/P)/dp*VO[6]           //
//    VO[ 7]   VO[ 8]   VO[ 9]   VO[10]   VO[11]   VO[12]   VO[13]   d()/dp1          //
//    VO[14]   VO[15]   VO[16]   VO[17]   VO[18]   VO[19]   VO[20]   d()/dp2          //
//    ...........................................................    d()/dpND         //
//                                                                                    //
// Output parameters:                                                                 //
//                                                                                    //
//    VO   -  output parameters and derivatives after propogation in magnetic field   //
//            defined by Mfield (KGauss)                                              //
//    Where a Mfield(R,H) - is interface to magnetic field information                //
//    input  R[ 0],R[ 1],R[ 2] - X     , Y      and Z  of the track                   //
//    output H[ 0],H[ 1],H[ 2] - Hx    , Hy     and Hz of the magnetic field          //
//           H[ 3],H[ 4],H[ 5] - dHx/dx, dHx/dy and dHx/dz                            //
//           H[ 6],H[ 7],H[ 8] - dHy/dx, dHy/dy and dHy/dz                            //
//           H[ 9],H[10],H[11] - dHz/dx, dHz/dy and dHz/dz                            //
//                                                                                    //
//    Path - trajectory length                                                        //
//                                                                                    //
// Authors: R.Brun, M.Hansroul, V.Perevoztchikov (Geant3)                             //
//                                                                                    //
// Modified by I.Gavrilenko for XKALMAN++ 29/06/2000                                  //
//                                                                                    //
// Minor modification by Sergei.Gerassimov@cern.ch for use in Traffic                 //
// Here is assumed, that field is smooth enough i.e. filed gradients are negligible   //
//------------------------------------------------------------------------------------//
//                                                                                    //
#include <iomanip>
#include <iostream>
#include <math.h>
#include "TAlgo.h"
#include "TDisplay.h"

bool TAlgo::RkutaNoGrad (double* SU,double* VO, double& Path) {

  const double EC     = .00014989626;                // 
  const double DLT    = .0002;                       // 
  const double DLT32  = DLT/32.;                     //
  const double Sstop  = .001;                        // Min. step 
  const double pi     = 3.1415927;                   // Pi
  const double P3     = .33333333;                   //
  const double Smax   = 100.;                        // max. step allowed>0 
  const double Wmax   = 5000.;                       // max. way allowed
  const int    ND     = 42, ND1=ND-7;                //
  double*  R          = &VO[0];                      // Start coordinates
  double*  A          = &VO[3];                      // Start directions
  double SA[3]        = {0.,0.,0.};                  // Start directions derivatives 
  double  Pinv        = VO[6]*EC;                    // Invert momentum/2.
  double  Way         = 0.;                          // Total way of the trajectory
  int     Error       = 0;                           // Error of propogation
  //
  // Step estimation until surface
  //
  double Step,An,Dist,Dis,S,Sl=0;
  Dist=SU[3]-R[0]*SU[0]-R[1]*SU[1]-R[2]*SU[2]; An=A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2];
  if (An==0) return(false); Step=Dist/An; 
  if(fabs(Step)>Wmax) {
    std::cout<<"TAlgo::RkutaNoGrad ==> Too long extrapolation requested : "<<Step<<" cm !"<<std::endl;
    std::cout<<"X = "<<R[0]<<" Y = "<<R[1]<<" Z = "<<R[2]
	<<"  COSx = "<<A[0]<<"  COSy = "<<A[1]<<"  COSz = "<<A[2]<<std::endl;
    std::cout<<"Destination  X = "<<SU[0]*SU[3]<<std::endl;
    return(false);
  }
  Step>Smax ? S=Smax : Step<-Smax ? S=-Smax : S=Step;
  //
  // Main cycle of Runge-Kutta method
  //
  //
  while(fabs(Step)>Sstop) {
    double H0[12],H1[12],H2[12],r[3];
    double S3=P3*S, S4=.25*S, PS2=Pinv*S; 
    //
    // First point
    //   
    r[0]=R[0]      ; r[1]=R[1]      ; r[2]=R[2]      ;  TAlgo::Field(r,H0); 
    if(TDisplay::Ptr() != NULL) { // if TDisplay object exists 
      if(TDisplay::Ref().Rkutraj()) TDisplay::Ref().point(r,108); // draw point if requested
    }
    H0[0]*=PS2     ; H0[1]*=PS2     ; H0[2]*=PS2     ; 
    double A0=A[1]*H0[2]-A[2]*H0[1], B0=A[2]*H0[0]-A[0]*H0[2], C0=A[0]*H0[1]-A[1]*H0[0];
    double A2=A[0]+A0              , B2=A[1]+B0              , C2=A[2]+C0              ;
    double A1=A2+A[0]              , B1=B2+A[1]              , C1=C2+A[2]              ;
    //
    // Second point
    //
    r[0]+=A1*S4    ; r[1]+=B1*S4    ; r[2]+=C1*S4    ;   TAlgo::Field(r,H1);    
    if(TDisplay::Ptr() != NULL) { // if TDisplay object exists 
      if(TDisplay::Ref().Rkutraj()) TDisplay::Ref().point(r,108); // draw point if requested
    }
    H1[0]*=PS2     ; H1[1]*=PS2     ; H1[2]*=PS2     ; 
    double A3,B3,C3,A4,B4,C4,A5,B5,C5;
    A3 = B2*H1[2]-C2*H1[1]+A[0]; B3=C2*H1[0]-A2*H1[2]+A[1]; C3=A2*H1[1]-B2*H1[0]+A[2];
    A4 = B3*H1[2]-C3*H1[1]+A[0]; B4=C3*H1[0]-A3*H1[2]+A[1]; C4=A3*H1[1]-B3*H1[0]+A[2];
    A5 = A4-A[0]+A4            ; B5=B4-A[1]+B4            ; C5=C4-A[2]+C4            ;
    //
    // Last point
    //
    r[0]=R[0]+S*A4 ; r[1]=R[1]+S*B4 ; r[2]=R[2]+S*C4 ;  TAlgo::Field(r,H2);
    if(TDisplay::Ptr() != NULL) { // if TDisplay object exists 
      if(TDisplay::Ref().Rkutraj()) TDisplay::Ref().point(r,108); // draw point if requested
    }
    H2[0]*=PS2     ; H2[1]*=PS2     ; H2[2]*=PS2     ; 
    double A6=B5*H2[2]-C5*H2[1], B6=C5*H2[0]-A5*H2[2], C6=A5*H2[1]-B5*H2[0];
    //
    // Test approximation quality on give step and possible step reduction
    //
    double EST;
    if((EST=fabs((A1+A6)-(A3+A4))+fabs((B1+B6)-(B3+B4))+fabs((C1+C6)-(C3+C4)))>DLT) {S*=.5; continue;}
    //
    // Derivatives of track parameters in last point
    //
    for(int i=7; i!=ND; i+=7) {
      double* dR = &VO[i];  double* dA = &VO[i+3]; double dH0,dH1,dH2;
      double dA0   = H0[ 2]*dA[1]-H0[ 1]*dA[2];              // dA0/dp
      double dB0   = H0[ 0]*dA[2]-H0[ 2]*dA[0];              // dB0/dp
      double dC0   = H0[ 1]*dA[0]-H0[ 0]*dA[1];              // dC0/dp
      if(i==ND1) {dA0+=A0; dB0+=B0; dC0+=C0;}
      double dA2   = dA0+dA[0]; 
      double dB2   = dB0+dA[1]; 
      double dC2   = dC0+dA[2];
      double dA3   = dA[0]+dB2*H1[2]-dC2*H1[1];              // dA3/dp
      double dB3   = dA[1]+dC2*H1[0]-dA2*H1[2];              // dB3/dp
      double dC3   = dA[2]+dA2*H1[1]-dB2*H1[0];              // dC3/dp
      if(i==ND1) {dA3+=A3-A[0]; dB3+=B3-A[1]; dC3+=C3-A[2];}
      double dA4   = dA[0]+dB3*H1[2]-dC3*H1[1];              // dA4/dp
      double dB4   = dA[1]+dC3*H1[0]-dA3*H1[2];              // dB4/dp
      double dC4   = dA[2]+dA3*H1[1]-dB3*H1[0];              // dC4/dp
      if(i==ND1) {dA4+=A4-A[0]; dB4+=B4-A[1]; dC4+=C4-A[2];}
      double dA5   = dA4+dA4-dA[0];
      double dB5   = dB4+dB4-dA[1];
      double dC5   = dC4+dC4-dA[2];
      double dA6   = dB5*H2[2]-dC5*H2[1];                    // dA6/dp
      double dB6   = dC5*H2[0]-dA5*H2[2];                    // dB6/dp
      double dC6   = dA5*H2[1]-dB5*H2[0];                    // dC6/dp
      if(i==ND1) {dA6+=A6; dB6+=B6; dC6+=C6;}
      dR[0]+=(dA2+dA3+dA4)*S3; dA[0] = (dA0+dA3+dA3+dA5+dA6)*P3;      
      dR[1]+=(dB2+dB3+dB4)*S3; dA[1] = (dB0+dB3+dB3+dB5+dB6)*P3; 
      dR[2]+=(dC2+dC3+dC4)*S3; dA[2] = (dC0+dC3+dC3+dC5+dC6)*P3;
    }
    if((Way+=fabs(S))>Wmax){ 
      std::cout<<"TAlgo::RkutaNoGrad ==> Trajectory is longer then length limit : "<<Way<<" cm !"
	  <<" P = "<<1./VO[6]<< " GeV"<<std::endl;
      return(false);
    }
    //
    // Track parameters in last point
    //   
    R[0]+=(A2+A3+A4)*S3; A[0]+=(SA[0]=(A0+A3+A3+A5+A6)*P3-A[0]); 
    R[1]+=(B2+B3+B4)*S3; A[1]+=(SA[1]=(B0+B3+B3+B5+B6)*P3-A[1]);
    R[2]+=(C2+C3+C4)*S3; A[2]+=(SA[2]=(C0+C3+C3+C5+C6)*P3-A[2]); Sl=S;
    double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
    A[0]*=CBA; A[1]*=CBA; A[2]*=CBA;
    //
    // Step estimation until surface and test conditions for stop propogation
    //
    Dis = SU[3]-R[0]*SU[0]-R[1]*SU[1]-R[2]*SU[2]; An=A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2];
    if (An==0 || (Dis*Dist>0 && fabs(Dis)>fabs(Dist))) {Error=1; Step=0; break;}   
    Step = Dis/An; Dist=Dis;
    //
    // Possible current step reduction
    //

    if(EST<DLT32                      ) S*=2.;
    if(S*Step<0. || fabs(S)>fabs(Step)) S=Step;

  } //end of main loop

  //
  // Output information preparation for main track parameteres
  //
  if(Sl!=0) Sl=1./Sl;
  A [0]+=(SA[0]*=Sl)*Step; A [1]+=(SA[1]*=Sl)*Step; A [2]+=(SA[2]*=Sl)*Step;
  double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
  VO[0]      = R[0]+Step*(A[0]-.5*Step*SA[0]);                            
  VO[1]      = R[1]+Step*(A[1]-.5*Step*SA[1]);
  VO[2]      = R[2]+Step*(A[2]-.5*Step*SA[2]);
  VO[3]      = A[0]*CBA;
  VO[4]      = A[1]*CBA;
  VO[5]      = A[2]*CBA;
  //
  // Output derivatives of track parameters preparation 
  //
  An = A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2]; An!=0 ? An=1./An : An = 0;
  for(int i=7; i!=ND; i+=7) {
    double* dR = &VO[i];  double* dA = &VO[i+3];
    S = (dR[0]*SU[0]+dR[1]*SU[1]+dR[2]*SU[2])*An;
    dR[0]-=S*A [0];  dR[1]-=S*A [1]; dR[2]-=S*A [2]; 
    dA[0]-=S*SA[0];  dA[1]-=S*SA[1]; dA[2]-=S*SA[2];
  }
  if(Error == 1){
    std::cout<<"TAlgo::RkutaNoGrad ==> Do not getting closer. Path = "
	<<Way<<" cm"<<"  P = "<<1./VO[6]<<" GeV"<<std::endl;
    return(false);
  }

  Path=fabs(Way);
  return(true);
}












