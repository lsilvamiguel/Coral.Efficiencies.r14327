// $Id: THlxUpdate.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!

  Member function 
  for Kalman's update this->Helix by Helix H1
  
  \return Hout - updated helix; Chi2 - Chi2 increment

*/

#include <iostream>
#include <stdio.h>
#include <cassert>
#include "THlx.h"
#include "TMtx.h"
#include "THit.h"
#include "TConstants.h"

using namespace std;

bool THlx::Update(const THlx& H1, THlx& Hout, double& Chi2) const

{
  int ier=0;
  double x0,x1;
  TMtx D(5), D0(5), D1(5), V0(5), C0(5,5), W0(5,5), V1(5), C1(5,5), W1(5,5), Vupd(5), Cupd(5,5), Wupd(5,5);
  TMtx Vtmp(5);

  this->Get(x0,V0,C0);
  H1.   Get(x1,V1,C1);

  if(fabs(x1-x0) > 5 * TConstants_RKuttaMinStep) {
    cout<<"THlx::Update ==> Helices are at different X :  X0 = "<<x0<<"   X1 = "<<x1<<endl;
    return(false);
  }

  W0=C0.i(ier);
  if(ier!=0){
    cout<<"THlx::Update ==> Cannot invert C0 matrix. Error = "<<ier<<endl;
    C0.Print("C0");
    assert(false);
  }
  W1=C1.i(ier);
  if(ier!=0){
    cout<<"THlx::Update ==> Cannot invert C1 matrix. Error = "<<ier<<endl;
    C1.Print("C1");
    assert(false);
  }

  Wupd=W0+W1;

  Cupd=Wupd.i(ier);

  if(ier!=0){
    cout<<"THlx::Update ==> Cannot invert Wupd matrix. Error = "<<ier<<endl;
    Wupd.Print("Wupd");
    assert(false);
  }

  D=V1-V0;

  D0 = Cupd*(W1*D);
  Chi2=D0.t()*W0*D0;
  
  D1=D0-D;
  Chi2 += D1.t()*W1*D1;

  Vupd = V0+D0;

  x1=(x0+x1)/2.;
  Hout.Set(x1,Vupd, Cupd);


  return(true);
}

/*!
  Member function   
  for Kalman's Update of helix by coordinate detector measurement 
  
  \return Hout  - updated helix
*/

bool THlx::Update(const THit& h, THlx& Hout) const

{
/*

 //
 // Straight forward method using matrix class (for tests only)
 //
 
  double x0, dChi2;
  int ier;
  TMtx V(5),    C(5,5); 
  TMtx Vm(5),   G(5,5), U(5,5);
  TMtx Vupd(5), Cupd(5,5);
  TMtx I(5,5), K0(5,5), K(5,5), R(5,5);

  I*=0;I(1,1)=I(2,2)=I(3,3)=I(4,4)=I(5,5)=1.;

  this->Get(x0,V,C);

  Vm*=0; // Measurement vector
  Vm(1)=h.U; Vm(2)=h.V; 

  G*=0;  // Mesurement weight matrix
  G(1,1)=1./(h.SigU*h.SigU); 
  G(2,2)=1./(h.SigV*h.SigV); 

  R*=0;  // Rotation to measurement system matrix
  R(1,1)=R(2,2)=R(3,3)=R(4,4)= h.GetDet().Ca;
  R(1,2)=R(3,4)=   h.GetDet().Sa;
  R(2,1)=R(4,3)= - h.GetDet().Sa;
  R(5,5)=1.;

 
  if(C(5,5) < 0) C.Print("Bad cov mtx");

  U = R.t()*G*R;
  K0=(I+C*U);

  K=K0.i(ier);
  if(ier!=0){
    cout<<"THlx::Update ==> Cannot invert K matrix"<<endl;
    G.Print("G");
    R.Print("R");
    U.Print("U");
    C.Print("C");
    K0.Print("K0");
    return(false);
  }

  Cupd=K*C;

//  Vupd=K*V+Cupd*R.t()*G*Vm;

  Vupd = V + Cupd*U*(R.t()*Vm - V);

  Hout.Set(x0,Vupd,Cupd);

  dChi2=
    (Vupd-V).t()*C.i(ier)*(Vupd-V) + (R*Vupd-Vm).t()*G*(R*Vupd-Vm);

*/


//  C' = (1+C*G).i * C
//  X' = X + C'*G*(Xm-X)
// as recomended by I.Gavrilenko (ATLAS)

  double w1 = h.G0;
  double w2 = h.G2;
  double w3 = h.G1;

  double a1 = h.Y - Hpar[1];
  double a2 = h.Z - Hpar[2];

  double d2 = w1*w2 - w3*w3;
  double d1 = Hcov[0]*Hcov[2] - Hcov[1]*Hcov[1];
  double d = 1.0 + Hcov[0]*w1 + Hcov[2]*w2 + 2.0*Hcov[1]*w3 + d2*d1;

  if(d <= 0.) {
    cout<<"THlxUpdate ==> D < 0   : "<<d<<endl;
    return(false);
  }

  d = 1.0/d;

  double b11 = Hcov[1]*w3 + Hcov[2]*w2;
  double b21 = Hcov[1]*w1 + Hcov[2]*w3;
  double b31 = d2*(Hcov[1]*Hcov[4] - Hcov[2]*Hcov[3] ) - Hcov[3] *w1 - Hcov[4] *w3;
  double b41 = d2*(Hcov[1]*Hcov[7] - Hcov[2]*Hcov[6] ) - Hcov[6] *w1 - Hcov[7] *w3;
  double b51 = d2*(Hcov[1]*Hcov[11]- Hcov[2]*Hcov[10]) - Hcov[10]*w1 - Hcov[11]*w3;
  
  double b12 = Hcov[0]*w3 + Hcov[1]*w2;
  double b22 = Hcov[0]*w1 + Hcov[1]*w3;
  double b32 = d2*(Hcov[1]*Hcov[3] - Hcov[0]*Hcov[4] ) - Hcov[3] *w3 - Hcov[4] *w2;
  double b42 = d2*(Hcov[1]*Hcov[6] - Hcov[0]*Hcov[7] ) - Hcov[6] *w3 - Hcov[7] *w2;
  double b52 = d2*(Hcov[1]*Hcov[10]- Hcov[0]*Hcov[11]) - Hcov[10]*w3 - Hcov[11]*w2;

  Hout.Hcov[0] = d*(Hcov[0] + d1*w2);
  Hout.Hcov[1] = d*(Hcov[1] - d1*w3);
  Hout.Hcov[3] = d*(Hcov[3] + b11*Hcov[3] - b12*Hcov[4]);
  Hout.Hcov[6] = d*(Hcov[6] + b11*Hcov[6] - b12*Hcov[7]);
  Hout.Hcov[10]= d*(Hcov[10]+ b11*Hcov[10]- b12*Hcov[11]);
  Hout.Hcov[2] = d*(Hcov[2] + d1*w1);
  Hout.Hcov[4] = d*(Hcov[4] - b21*Hcov[3] +b22*Hcov[4]);
  Hout.Hcov[7] = d*(Hcov[7] - b21*Hcov[6] +b22*Hcov[7]);
  Hout.Hcov[11]= d*(Hcov[11]- b21*Hcov[10]+b22*Hcov[11]);

  Hout.Hcov[5] = Hcov[5] + d*(b31*Hcov[3] +b32*Hcov[4]);
  Hout.Hcov[8] = Hcov[8] + d*(b31*Hcov[6] +b32*Hcov[7]);
  Hout.Hcov[12]= Hcov[12]+ d*(b31*Hcov[10]+b32*Hcov[11]);

  Hout.Hcov[9] = Hcov[9] + d*(b41*Hcov[6] +b42*Hcov[7]);
  Hout.Hcov[13]= Hcov[13]+ d*(b41*Hcov[10]+b42*Hcov[11]);

  Hout.Hcov[14]= Hcov[14]+ d*(b51*Hcov[10]+b52*Hcov[11]);
 
  if(Hout.Hcov[0] < 0 || Hout.Hcov[2] < 0 ||
     Hout.Hcov[5] < 0 || Hout.Hcov[9] < 0 ||
     Hout.Hcov[14]< 0){
    cout<<"THlxUpdate ==> Diagonal element of cov. matrix < 0"<<endl;;
    this->Print("Input helix");
    cout<<"Hit : "<<endl;
    printf("U = %8.2f V = %8.2f G0 = %9.2e G2 = %9.2e G1 = %9.2e \n",
	   h.U, h.V, h.G0, h.G2, h.G1);
    return(false);
  }

  double e1 = w1*a1 + w3*a2;
  double e2 = w3*a1 + w2*a2;
  double g1 = Hout.Hcov[0]*e1 + Hout.Hcov[1]*e2;
  double g2 = Hout.Hcov[1]*e1 + Hout.Hcov[2]*e2;

  a1 = g1-a1;
  a2 = g2-a2;

  double dChi2 = w1*a1*a1 + a2*(w2*a2 + 2.0*w3*a1) +
  (Hcov[2]*g1*g1 + g2*(Hcov[0]*g2 - 2.0*Hcov[1]*g1))/d1;


  Hout.Hpar[0] = Hpar[0];
  Hout.Hpar[1] = Hpar[1] + g1;
  Hout.Hpar[2] = Hpar[2] + g2;
  Hout.Hpar[3] = Hpar[3] + Hout.Hcov[3] *e1 + Hout.Hcov[4] *e2;;
  Hout.Hpar[4] = Hpar[4] + Hout.Hcov[6] *e1 + Hout.Hcov[7] *e2;;
  Hout.Hpar[5] = Hpar[5] + Hout.Hcov[10]*e1 + Hout.Hcov[11]*e2;;

 

  return(true);
}






