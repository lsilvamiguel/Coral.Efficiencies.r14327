
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Coral.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TAlgo.h"
#include "THlx.h"
#include "TConstants.h"
#include "TEv.h"
#include "TSetup.h"
#include "CA.h"

using namespace std;

void kalman_propagate(ftype R[6], ftype G[6][6], ftype dx){
  R[1]   +=dx*R[3];
  R[2]   +=dx*R[4];
  G[1][1]+=dx*(2.0*G[1][3]+dx*G[3][3]);
  G[1][2]+=dx*(G[2][3]+G[1][4]+dx*G[3][4]);
  G[1][3]+=dx*G[3][3];
  G[1][4]+=dx*G[3][4];
  G[2][2]+=dx*(2.0*G[2][4]+dx*G[4][4]);
  G[2][3]+=dx*G[3][4];
  G[2][4]+=dx*G[4][4];
}
 
ftype kalman_filter(ftype R[6], ftype G[6][6], ftype sina, ftype cosa,
		   ftype measurement, ftype sigma0_2){
  ftype B[6];
  B[1]=cosa*G[1][1]+sina*G[1][2];
  B[2]=cosa*G[1][2]+sina*G[2][2];
  B[3]=cosa*G[1][3]+sina*G[2][3];
  B[4]=cosa*G[1][4]+sina*G[2][4];
  ftype u1 = cosa*R[1]+sina*R[2];
  ftype sigma2 = sqr(cosa)*G[1][1]+2.0*sina*cosa*G[1][2]+sqr(sina)*G[2][2];
  ftype s = sigma0_2+sigma2;
  ftype zeta_k=(measurement-u1);
  for(int j=1;j<=4;j++){
    R[j]+=(B[j]/s)*zeta_k;
    for(int i=1;i<=j;i++)
      G[i][j]-=B[i]*(B[j]/s);
  }
  ftype u2 = cosa*R[1]+sina*R[2];
  return sqr(measurement-u2)/sigma0_2+sqr(u1-u2)/sigma2;
} 

void kanman_ms_correction(ftype R[6], ftype G[6][6], int ipl, ftype pinv){

  const TDetect& d = TSetup::Ref().iPlane2Detect(ipl);

  ftype a = 1/sqrt(1.+ sqr(R[3]) + sqr(R[4]));
  ftype x = d.Siz(0)/a;
  ftype len   =  x / d.RadLen;

  // Lynch and Dahl aproximation for Sigma(Theta_proj) of mult. scatt.
  ftype SigTheta = 0.0136*fabs(pinv) * sqrt(len) * (1.+0.038*log(len));
 
  // Noise matrix calculation (NIM A329 (1993) 493-500)
  // Transverse displacement of the track is ignored.

  ftype h = TOpt::CAOptD[20]*sqr(SigTheta)*(1+sqr(R[3])+sqr(R[4]));
  G[3][3]+=h*(1+sqr(R[3]));
  G[4][4]+=h*(1+sqr(R[4]));
  G[3][4]+=h*R[3]*R[4];
}

void kalman_propagateQK(ftype R[6], ftype G[6][6], ftype x1, ftype x2){

  ftype xfirst = x1;
  ftype xlast = x2;
  
  
  //
  // Undefined momentum or too low momentum for Runge-Kutta propagation. 
  // cut on  P < 100 MeV && P>1000 GeV
  // Do straight line extrapolation.
  //

  if((fabs(R[5])< 0.0001)||(fabs(R[5]) > 10)) goto StraightLine;

  if(fabs(xlast-xfirst) < TConstants_RKuttaMinStep) { // very short step
    goto StraightLine;
  }
   
  //
  // Runge-Kutta extrapolation
  //
  {
    // Prepare surface description array
    double su[] = { ((xlast<0) ?-1. :1.), 0, 0, fabs(xlast)};
  
    // Prepare track parameters
  
    ftype x  = xfirst;
    ftype y  = R[1];
    ftype z  = R[2];
    ftype ax = 1./sqrt(sqr(R[3]) + sqr(R[4]) + 1);
    ftype ay = R[3]*ax;
    ftype az = R[4]*ax;
    ftype qP = R[5];
  
    // Prepare Jacobian dP/dH where 
    //
    // P is track parameters vector (size 7), needed for RkutaBP routine,
    // in the direction cosines (ax,ay,az) representation: (x,y,x,ax,ay,az,q/P)  
    //
    // H is helix.
    //

    double P[42] = {
      x,      y,      z,         ax,          ay,         az,     qP,   // Parameters (P)
      0,      1,      0,          0,           0,           0,     0,   // dP/dR[1]
      0,      0,      1,          0,           0,           0,     0,   // dP/dR[2]
      0,      0,      0,  -ay*ax*ax, ax-ay*ay*ax,   -ay*az*ax,     0,   // dP/dR[3]
      0,      0,      0,  -az*ax*ax,   -ay*az*ax, ax-az*az*ax,     0,   // dP/dR[4]
      0,      0,      0,          0,           0,           0,     qP   // dP/dR[5]*R[5]
    };
    
    // Do the propagation
    double path;
    if( !TAlgo::RkutaNoGrad(su,P,path) ) goto StraightLine;
    
    if(P[3] <= 0.) {
      cout<<"THlx::Extrap ==> Trajectory had turned back during extrapolation. P = "
	  <<fabs(1./P[6])<<" GeV"<<endl;
      goto StraightLine;
    }
    if(fabs(P[0]-xlast) > 0.001) {
      cout<<"THlx::Extrap ==> P[0]!=xlast\n";
      goto StraightLine;
    }
    // Output helix parameters
    
    R[1] = P[1]; // y
    R[2] = P[2]; // z
    R[3] = P[4]/P[3];    // yp = ay/ax
    R[4] = P[5]/P[3];    // zp = az/ax
    R[5] = P[6];

    // Calculate helix transformation Jacobian (F = dHout/dHin)
    
    ftype p3Ax = -P[4]/(P[3]*P[3]);
    ftype p3Ay = 1./P[3];
    ftype p3Az = 0.;
    ftype p4Ax = -P[5]/(P[3]*P[3]);
    ftype p4Ay = 0.;
    ftype p4Az = p3Ay;
        
    ftype f[5][5];
      
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
    
    ftype g[6][6];
    copy_matrix(G,g);
    
    ftype w[5];
    
    w[0]=g[1][1] * f[0][0] + g[1][2] * f[0][1] + g[1][3] * f[0][2] + g[1][4] * f[0][3] + g[1][5]*f[0][4];
    w[1]=g[1][2] * f[0][0] + g[2][2] * f[0][1] + g[2][3] * f[0][2] + g[2][4] * f[0][3] + g[2][5]*f[0][4];
    w[2]=g[1][3] * f[0][0] + g[2][3] * f[0][1] + g[3][3] * f[0][2] + g[3][4] * f[0][3] + g[3][5]*f[0][4];
    w[3]=g[1][4] * f[0][0] + g[2][4] * f[0][1] + g[3][4] * f[0][2] + g[4][4] * f[0][3] + g[4][5]*f[0][4];
    w[4]=g[1][5] * f[0][0] + g[2][5] * f[0][1] + g[3][5] * f[0][2] + g[4][5] * f[0][3] + g[5][5]*f[0][4];

    G[1][1]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
    
    w[0]=g[1][1] * f[1][0] + g[1][2] * f[1][1] + g[1][3] * f[1][2] + g[1][4] * f[1][3] + g[1][5]*f[1][4];
    w[1]=g[1][2] * f[1][0] + g[2][2] * f[1][1] + g[2][3] * f[1][2] + g[2][4] * f[1][3] + g[2][5]*f[1][4];
    w[2]=g[1][3] * f[1][0] + g[2][3] * f[1][1] + g[3][3] * f[1][2] + g[3][4] * f[1][3] + g[3][5]*f[1][4];
    w[3]=g[1][4] * f[1][0] + g[2][4] * f[1][1] + g[3][4] * f[1][2] + g[4][4] * f[1][3] + g[4][5]*f[1][4];
    w[4]=g[1][5] * f[1][0] + g[2][5] * f[1][1] + g[3][5] * f[1][2] + g[4][5] * f[1][3] + g[5][5]*f[1][4];

    G[1][2]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
    G[2][2]=w[0]*f[1][0] + w[1]*f[1][1] + w[2]*f[1][2] + w[3]*f[1][3] + w[4]*f[1][4];

    w[0]=g[1][1] * f[2][0] + g[1][2] * f[2][1] + g[1][3] * f[2][2] + g[1][4] * f[2][3] + g[1][5]*f[2][4];
    w[1]=g[1][2] * f[2][0] + g[2][2] * f[2][1] + g[2][3] * f[2][2] + g[2][4] * f[2][3] + g[2][5]*f[2][4];
    w[2]=g[1][3] * f[2][0] + g[2][3] * f[2][1] + g[3][3] * f[2][2] + g[3][4] * f[2][3] + g[3][5]*f[2][4];
    w[3]=g[1][4] * f[2][0] + g[2][4] * f[2][1] + g[3][4] * f[2][2] + g[4][4] * f[2][3] + g[4][5]*f[2][4];
    w[4]=g[1][5] * f[2][0] + g[2][5] * f[2][1] + g[3][5] * f[2][2] + g[4][5] * f[2][3] + g[5][5]*f[2][4];

    G[1][3]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
    G[2][3]=w[0]*f[1][0] + w[1]*f[1][1] + w[2]*f[1][2] + w[3]*f[1][3] + w[4]*f[1][4];
    G[3][3]=w[0]*f[2][0] + w[1]*f[2][1] + w[2]*f[2][2] + w[3]*f[2][3] + w[4]*f[2][4];
    
    w[0]=g[1][1] * f[3][0] + g[1][2] * f[3][1] + g[1][3] * f[3][2] + g[1][4] * f[3][3] + g[1][5]*f[3][4];
    w[1]=g[1][2] * f[3][0] + g[2][2] * f[3][1] + g[2][3] * f[3][2] + g[2][4] * f[3][3] + g[2][5]*f[3][4];
    w[2]=g[1][3] * f[3][0] + g[2][3] * f[3][1] + g[3][3] * f[3][2] + g[3][4] * f[3][3] + g[3][5]*f[3][4];
    w[3]=g[1][4] * f[3][0] + g[2][4] * f[3][1] + g[3][4] * f[3][2] + g[4][4] * f[3][3] + g[4][5]*f[3][4];
    w[4]=g[1][5] * f[3][0] + g[2][5] * f[3][1] + g[3][5] * f[3][2] + g[4][5] * f[3][3] + g[5][5]*f[3][4];

    G[1][4]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
    G[2][4]=w[0]*f[1][0] + w[1]*f[1][1] + w[2]*f[1][2] + w[3]*f[1][3] + w[4]*f[1][4];
    G[3][4]=w[0]*f[2][0] + w[1]*f[2][1] + w[2]*f[2][2] + w[3]*f[2][3] + w[4]*f[2][4];
    G[4][4]=w[0]*f[3][0] + w[1]*f[3][1] + w[2]*f[3][2] + w[3]*f[3][3] + w[4]*f[3][4];

    // as dPinv/dpar = 0 if par != Pinv
    G[1][5]=g[1][5]*f[0][0] + g[2][5]*f[0][1] + g[3][5]*f[0][2] + g[4][5]*f[0][3] + g[5][5]*f[0][4];
    G[2][5]=g[1][5]*f[1][0] + g[2][5]*f[1][1] + g[3][5]*f[1][2] + g[4][5]*f[1][3] + g[5][5]*f[1][4];
    G[3][5]=g[1][5]*f[2][0] + g[2][5]*f[2][1] + g[3][5]*f[2][2] + g[4][5]*f[2][3] + g[5][5]*f[2][4];
    G[4][5]=g[1][5]*f[3][0] + g[2][5]*f[3][1] + g[3][5]*f[3][2] + g[4][5]*f[3][3] + g[5][5]*f[3][4];

    G[5][5]=g[5][5];

    return;
  }
 StraightLine:
  //cout<<"\nFit with straight line\n";
  ftype dx = x2-x1;
  R[1]   +=dx*R[3];
  R[2]   +=dx*R[4];
  G[1][1]+=dx*(2.0*G[1][3]+dx*G[3][3]);
  G[1][2]+=dx*(G[2][3]+G[1][4]+dx*G[3][4]);
  G[1][3]+=dx*G[3][3];
  G[1][4]+=dx*G[3][4];
  G[2][2]+=dx*(2.0*G[2][4]+dx*G[4][4]);
  G[2][3]+=dx*G[3][4];
  G[2][4]+=dx*G[4][4];
  G[1][5]+=dx*G[3][5];
  G[2][5]+=dx*G[4][5];
  return;
}

ftype kalman_filterQK(ftype R[6], ftype G[6][6], TCAHit *hit){

  ftype Ca = hit->Layer->Ca; 
  ftype Sa = hit->Layer->Sa;
 
  ftype wu = 1./sqr(hit->sigma_u);
  ftype wv = 1./sqr(hit->sigma_v);
  ftype hit_y = Ca*hit->u - Sa*hit->v;
  ftype hit_z = Sa*hit->u + Ca*hit->v;
  ftype w1 = wu*Ca*Ca + wv*Sa*Sa;
  ftype w2 = wv*Ca*Ca + wu*Sa*Sa;
  ftype w3 = wu*Ca*Sa - wv*Ca*Sa;

  ftype a1 = hit_y - R[1];
  ftype a2 = hit_z - R[2];

  ftype d2 = w1*w2 - w3*w3;
  ftype d1 = G[1][1]*G[2][2] - G[1][2]*G[1][2];
  ftype d = 1.0 + G[1][1]*w1 + G[2][2]*w2 + 2.0*G[1][2]*w3 + d2*d1;

  if(d <= 0.) {
    cout<<"THlxUpdate ==> D < 0   : "<<d<<endl;
    return(false);
  }

  d = 1.0/d;

  ftype b11 = G[1][2]*w3 + G[2][2]*w2;
  ftype b21 = G[1][2]*w1 + G[2][2]*w3;
  ftype b31 = d2*(G[1][2]*G[2][3] - G[2][2]*G[1][3] ) - G[1][3] *w1 - G[2][3] *w3;
  ftype b41 = d2*(G[1][2]*G[2][4] - G[2][2]*G[1][4] ) - G[1][4] *w1 - G[2][4] *w3;
  ftype b51 = d2*(G[1][2]*G[2][5]- G[2][2]*G[1][5]) - G[1][5]*w1 - G[2][5]*w3;
  
  ftype b12 = G[1][1]*w3 + G[1][2]*w2;
  ftype b22 = G[1][1]*w1 + G[1][2]*w3;
  ftype b32 = d2*(G[1][2]*G[1][3] - G[1][1]*G[2][3] ) - G[1][3]*w3 - G[2][3]*w2;
  ftype b42 = d2*(G[1][2]*G[1][4] - G[1][1]*G[2][4] ) - G[1][4]*w3 - G[2][4]*w2;
  ftype b52 = d2*(G[1][2]*G[1][5] - G[1][1]*G[2][5])  - G[1][5]*w3 - G[2][5]*w2;

  ftype g[6][6];
  copy_matrix(G,g);

  G[1][1] = d*(g[1][1] + d1*w2);
  G[1][2] = d*(g[1][2] - d1*w3);
  G[1][3] = d*(g[1][3] + b11*g[1][3] - b12*g[2][3]);
  G[1][4] = d*(g[1][4] + b11*g[1][4] - b12*g[2][4]);
  G[1][5] = d*(g[1][5] + b11*g[1][5] - b12*g[2][5]);

  G[2][2] = d*(g[2][2] + d1*w1);
  G[2][3] = d*(g[2][3] - b21*g[1][3] + b22*g[2][3]);
  G[2][4] = d*(g[2][4] - b21*g[1][4] + b22*g[2][4]);
  G[2][5] = d*(g[2][5] - b21*g[1][5] + b22*g[2][5]);

  G[3][3] = g[3][3] + d*(b31*g[1][3] + b32*g[2][3]);
  G[3][4] = g[3][4] + d*(b31*g[1][4] + b32*g[2][4]);
  G[3][5] = g[3][5] + d*(b31*g[1][5] + b32*g[2][5]);

  G[4][4] = g[4][4] + d*(b41*g[1][4] + b42*g[2][4]);
  G[4][5] = g[4][5] + d*(b41*g[1][5] + b42*g[2][5]);

  G[5][5] = g[5][5] + d*(b51*g[1][5] + b52*g[2][5]);

  ftype e1 = w1*a1 + w3*a2;
  ftype e2 = w3*a1 + w2*a2;
  ftype g1 = G[1][1]*e1 + G[1][2]*e2;
  ftype g2 = G[1][2]*e1 + G[2][2]*e2;

  a1 = g1-a1;
  a2 = g2-a2;

  R[1] = R[1] + g1;
  R[2] = R[2] + g2;
  R[3] = R[3] + G[1][3]*e1 + G[2][3]*e2;
  R[4] = R[4] + G[1][4]*e1 + G[2][4]*e2;
  R[5] = R[5] + G[1][5]*e1 + G[2][5]*e2;

  return w1*a1*a1 + a2*(w2*a2 + 2.0*w3*a1) + (g[2][2]*g1*g1 + g2*(g[1][1]*g2 - 2.0*g[1][2]*g1))/d1;
} 

void kanman_ms_correctionQK(ftype R[6], ftype G[6][6], int ipl){

  if(R[5] == 0.) return; // momentum not known

  const TDetect& d = TSetup::Ref().iPlane2Detect(ipl);
  if(R[5] == 0.) return; // momentum not known
  
  ftype x = d.Siz(0)*sqrt(1.+ sqr(R[3]) + sqr(R[4]));
  ftype len   =  x / d.RadLen;

  // Lynch and Dahl aproximation for Sigma(Theta_proj) of mult. scatt.
  ftype SigTheta = 0.0136*fabs(R[5]) * sqrt(len) * (1.+0.038*log(len));

  // Noise matrix calculation (NIM A329 (1993) 493-500)
  // Transverse displacement of the track is ignored.

  ftype h = TOpt::CAOptD[20]*sqr(SigTheta)*(1+sqr(R[3])+sqr(R[4]));

  G[3][3]+=h*(1+sqr(R[3]));
  G[4][4]+=h*(1+sqr(R[4]));
  G[3][4]+=h*R[3]*R[4];
}

void TCATrack::CountPQ(){
  pinv = 4.0*sqrt(sqr((vR[3]+VR[3])/2)+sqr((vR[4]+VR[4])/2))/10;
  ftype Q = 0;
  ftype dx = Vx-vx;
  if(fabs(dx)>0.01){
    ftype Ay = (VR[1]-vR[1])/dx;
    ftype y0 = vR[1];
    ftype z0 = (vR[2]+VR[2])/2;
    ftype x0 = vx;
    for(list<TCAHit*>::iterator i=Hits.begin(); i!=Hits.end(); i++){
      TCAHit *hit = *i;
      ftype u = hit->u;
      ftype C1 = hit->Layer->Ca;
      ftype S1 = hit->Layer->Sa;
      ftype v = z0;
      ftype C2 = 0;
      ftype S2 = 1;
      ftype det = C1*S2 - S1*C2;
      if(fabs(det)<0.001) continue; // bad projection for measure y
      ftype y = (u*S2-v*S1)/det;
      ftype z =-(u*C2-v*C1)/det;
      Q+=(y0+Ay*(hit->x-x0)) - y;
    }
  }
  if(Q<0) Q=-1; else Q=+1;
  vR[5] = VR[5] = Q*pinv;
}

void TCATrack::Fit(ftype x1, ftype x2){

  CountPQ();
  chi2 = 0;

  for(int i=1; i<=4; i++){
    for(int j=1; j<=4; j++){
      vG[i][j]*=TOpt::CAOptD[21];
    }
  }
  //==================== Fit for the last point ======================
  ftype X = vx; // current position of filter

  for(list<TCAHit*>::iterator i=Hits.begin(); i!= Hits.end(); i++){
    TCAHit *hit = *i;
    if(hit==NULL) myexit("Fit : NULL hit in track");
    TLayer &layer = *(hit->Layer);
    ftype sigma = sqr(hit->sigma_u);
    kalman_propagate(vR,vG,layer.X0 - X);
    chi2+=kalman_filter(vR,vG,layer.Sa,layer.Ca,hit->u, sigma);// /hit->u;
    kalman_propagate(vR,vG,layer.size_x/2);
    kanman_ms_correction(vR,vG,layer.global_plane_index,pinv);
    X = layer.X0+layer.size_x/2;
  }

  copy_vector(vR,VR);
  copy_matrix(vG,VG);
  Vx = x2;
  kalman_propagate(VR,VG, Vx-X);
 
  //===================== Fit for the first point ======================

  for(int i=1; i<=4; i++){
    for(int j=1; j<=4; j++){
      vG[i][j]*=TOpt::CAOptD[21]; 
    }
  }
  //chi2 = 0;
  for(list<TCAHit*>::reverse_iterator i=Hits.rbegin(); i!=Hits.rend(); ++i){
    TCAHit *hit = *i;
    if(hit==NULL) myexit("Fit : NULL hit in track");
    TLayer &layer = *(hit->Layer);
    ftype sigma = sqr(hit->sigma_u);
    kanman_ms_correction(vR,vG,layer.global_plane_index,pinv);
    kalman_propagate(vR,vG, layer.X0 - X);
    chi2+=kalman_filter(vR,vG, layer.Sa,layer.Ca,hit->u, sigma);// /hit->u;
    kalman_propagate(vR,vG,-layer.size_x/2);
    X = layer.X0-layer.size_x/2;
  }
  vx = x1;
  kalman_propagate(vR, vG, vx - X);
  chi2/=2*NumberOfHits();
  pinv = fabs(vR[5]);
  
  if((vG[1][1]<0)||(vG[2][2]<0)||(vG[3][3]<0)||(vG[4][4]<0)){
    cout<<"CA: bad cov matrix " <<vG[1][1]<<" "<<vG[2][2]<<" "<<vG[3][3]<<" "<<vG[4][4]<<endl;
  }
}

void TCATrack::Fit(){
  Fit((*Hits.begin())->x, (*Hits.rbegin())->x);
}

void TCATrack::Propagate(ftype x, ftype R[6], ftype G[6][6]){
  ftype dx1 = x-vx;
  ftype dx2 = x-Vx;
  if(fabs(dx1)<fabs(dx2)) count_shift(vG, vR, G, R, dx1);
  else count_shift(VG, VR, G, R, dx2);
}

void TCATrack::FitQK(ftype x1, ftype x2){

  crosspoints.clear();
  CountPQ();
  chi2 = 0;

  for(int i=1; i<=5; i++){
    for(int j=1; j<=5; j++){
      if(i==j) vG[i][j]*=TOpt::CAOptD[22];
      else vG[i][j] = 0;
    }
  }
 
  //==================== Fit for the last point ======================
  ftype X = vx; // current position of filter
   
  for(list<TCAHit*>::iterator i=Hits.begin(); i!= Hits.end(); i++){
    TCAHit *hit = *i;
    if(hit==NULL) myexit("Fit : NULL hit in track");
    TLayer &layer = *(hit->Layer);
    kalman_propagateQK(vR,vG,X, layer.X0);
    X = layer.X0;
    chi2+=kalman_filterQK(vR,vG,hit)/hit->u;
    kalman_propagateQK(vR,vG,X, X+layer.size_x/2);
    kanman_ms_correctionQK(vR,vG,layer.global_plane_index);
    X += layer.size_x/2;
  }

  copy_vector(vR,VR);
  copy_matrix(vG,VG);  
  Vx = X;
  Vx = x2;
  kalman_propagateQK(VR,VG,X,Vx);
 
  //===================== Fit for the first point ======================

  for(int i=1; i<=5; i++){
    for(int j=1; j<=5; j++){
      if(i==j) vG[i][j]*=TOpt::CAOptD[22];
      else vG[i][j] = 0;
    }
  }
  //chi2=0;
  for(list<TCAHit*>::reverse_iterator i=Hits.rbegin(); i!=Hits.rend(); ++i){
    TCAHit *hit = *i;
    if(hit==NULL) myexit("Fit : NULL hit in track");
    TLayer &layer = *(hit->Layer);
    kanman_ms_correctionQK(vR,vG,layer.global_plane_index);
    kalman_propagateQK(vR,vG, X,layer.X0);
    X = layer.X0;
    chi2+=kalman_filterQK(vR,vG, hit)/hit->u;
    ftype u = layer.Ca*vR[1]+layer.Sa*vR[2];
    crosspoints.push_front(u);
    kalman_propagateQK(vR,vG,X,X-layer.size_x/2);
    X -= layer.size_x/2;
  }
  vx = x1;
  kalman_propagateQK(vR, vG,X, vx);
  chi2/=2*NumberOfHits();
  pinv = fabs(vR[5]);
  
  if((vG[1][1]<0)||(vG[2][2]<0)||(vG[3][3]<0)||(vG[4][4]<0)){
    cout<<"CA: bad cov matrix " <<vG[1][1]<<" "<<vG[2][2]<<" "<<vG[3][3]<<" "<<vG[4][4]<<endl;
  }
}

void TCATrack::FitQK(){
  FitQK((*Hits.begin())->x, (*Hits.rbegin())->x);
}

void TCATrack::PropagateQK(ftype x, ftype R[6], ftype G[6][6]){
  ftype dx1 = x-vx;
  ftype dx2 = x-Vx;
  if(fabs(dx1)<fabs(dx2)){
    copy_vector(vR,R);
    copy_matrix(vG,G);
    kalman_propagateQK(R,G, vx, x);
  }
  else{
    copy_vector(VR,R);
    copy_matrix(VG,G);
    kalman_propagateQK(R,G, Vx, x);
  }
}

