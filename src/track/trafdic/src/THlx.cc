// $Id: THlx.cc 14094 2015-11-06 15:28:48Z lsilva $

#include <iostream>
#include <iomanip>
#include <math.h>
#include "THlx.h"
#include "CsHelix.h"
#include "TMtx.h"

using namespace std;

/*! 
  Print helix
  \param str comment to be printed
*/ 
void THlx::Print(const char* str) const
{

  cout.fill(' ');
  cout.setf(ios::showpos);
  cout<<str<<" X = "<<setprecision(7)<<setw(10)<<
    Hpar[0]<<" cm"<<endl;
  cout.setf(ios::scientific);
  cout<<setprecision(3)<<setw(10)<<
    Hpar[1]<<" | "<<Hcov[0]<<endl;
  cout<<setprecision(3)<<setw(10)<<
    Hpar[2]<<" | "<<Hcov[1]<<" "<<Hcov[2]<<endl;
  cout<<setprecision(3)<<setw(10)<<
    Hpar[3]<<" | "<<Hcov[3]<<" "<<Hcov[4]<<" "<<Hcov[5]<<endl;
  cout<<setprecision(3)<<setw(10)<<
    Hpar[4]<<" | "<<Hcov[6]<<" "<<Hcov[7]<<" "<<Hcov[8]<<" "<<Hcov[9]<<endl;
  cout<<setprecision(3)<<setw(10)<<
    Hpar[5]<<" | "<<Hcov[10]<<" "<<Hcov[11]<<" "<<Hcov[12]<<" "<<Hcov[13]<<" "<<Hcov[14]<<endl;

  //Reset some format flags
  cout.setf(ios::fixed,ios::scientific);
  cout.unsetf(ios::showpos);
  cout<<setprecision(7);

  if(Hpar[5] != 0) cout<<"P = "<<1./Hpar[5]<<" +- "<<sqrt(Hcov[14])/(Hpar[5]*Hpar[5])<<endl;; 

}

//! Pack vector and matrix into helix
void THlx::Set(double& x, TMtx& V, TMtx& M)
{
  int i,j;
  Hpar[0]=x;
  for (i=1; i<=5; i++){
    Hpar[i]=V(i);
    for(j=1; j<=i; j++){
      Hcov[(j-1)+((i-1)*i)/2]=M(i,j);
    }
  }
}

//! Upack helix to vector and matrix
void THlx::Get(double& x, TMtx& V, TMtx& M) const
{
  register short int i,j;
  x=Hpar[0];
  for (i=1; i<=5; i++){
    V(i)=Hpar[i];
    for(j=1; j<=i; j++){
      M(i,j)=M(j,i)=Hcov[(j-1)+((i-1)*i)/2];
    }
  }
}

//! Converts CsHelix to "this" helix
void THlx::ImportHelix(const CsHelix& h)
{
  Hpar[0] = h.getZ()/10.;
  Hpar[1] = h.getX()/10.;
  Hpar[2] = h.getY()/10.;
  Hpar[3] = h.getDXDZ();
  Hpar[4] = h.getDYDZ();
  Hpar[5] = h.getCop();

  const double* cov = h.getCov();
  Hcov[0] = cov[0] /100.;
  Hcov[1] = cov[1] /100.; Hcov[2] = cov[2]/100.;
  Hcov[3] = cov[3] /10. ; Hcov[4] = cov[4] /10.; Hcov[5] = cov[5];
  Hcov[6] = cov[6] /10. ; Hcov[7] = cov[7] /10.; Hcov[8] = cov[8];  Hcov[9] = cov[9];
  Hcov[10]= cov[10]/10. ; Hcov[11]= cov[11]/10.; Hcov[12]= cov[12]; Hcov[13]= cov[13]; Hcov[14]=cov[14];
}

//! Converts "this" helix to CsHelix
CsHelix THlx::ExportHelix()
{
  double x,y,z,xp,yp,p;
  double cov[15];

  x = Hpar[1]*10.;
  y = Hpar[2]*10.;
  z = Hpar[0]*10.;
  xp= Hpar[3];
  yp= Hpar[4];
  p = Hpar[5];

  cov[0] = Hcov[0] *100.;
  cov[1] = Hcov[1] *100.; cov[2] = Hcov[2] *100.;
  cov[3] = Hcov[3] *10. ; cov[4] = Hcov[4] *10.; cov[5] = Hcov[5];
  cov[6] = Hcov[6] *10. ; cov[7] = Hcov[7] *10.; cov[8] = Hcov[8];  cov[9] = Hcov[9];
  cov[10]= Hcov[10]*10. ; cov[11]= Hcov[11]*10.; cov[12]= Hcov[12]; cov[13]= Hcov[13]; cov[14]=Hcov[14];

  CsHelix H(x,y,z,xp,yp,p,(const double *)cov);
  return(H);
}
