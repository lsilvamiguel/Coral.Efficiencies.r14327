/*!

 Member function of class Helix
 to extrapolate it to position of
 helix Hout (at the beginning, only fixed coordinate have to be preset). 
 Hout will be overwritten by resulting helix

 Numerical Helix extrapolation with
 numericaly calculated Jacobian, used in cov. matrix propagation.

 If Hcov[0] is set to 0, cov. matrix is not propagated

 Obsolete.

*/

#include <iostream>
#include <iomanip>
#include "TOpt.h"
#include "TAlgo.h"
#include "TConstants.h"
#include "TDisplay.h"
#include "THlx.h"
#include "TMtx.h"

using namespace std;

bool THlx::NpassExtrap(THlx& Hout) const
{

  double xlast = Hout.Hpar[0];
  double xfirst= Hpar[0];

  if(fabs(xlast-xfirst) < TConstants_RKuttaMinStep) { // very short step
    Hout = (*this);
    return(true);
  }

  int i,j,k,ier;
  double dx, dummy[3];
  double f[5][5];
  double w[5];
  double vin[6],vout[6], vecin[6], vecout[6];
  double xx;

  const int nstp = int(fabs(xlast-xfirst)/TOpt::dCut[3]); // N steps for field probing 
  bool helix; // "false" means "use straight line extrapolation"

  // increments for numerical Jacobian calculation
  const double dvin[6]={0.0, 1.e-4, 1.e-4, 1.e-5, 1.e-5, 1.e-5};

  // infinit momentum ? (q/P == 0)
  if(Hpar[5] == 0){ helix=false; goto do_extrap; }

  helix=true; goto do_extrap; ///tmp - always R.-Kutta extrap.
  
  // Probe the field assuming straight line trajectory
  xx = xfirst;
  for(int ii=0; ii <= nstp; ii++){
    dummy[0]=xx; 
    dx=dummy[0]-xfirst;
    dummy[1]=Hpar[1]+Hpar[3]*dx;
    dummy[2]=Hpar[2]+Hpar[4]*dx;
    if(TAlgo::Field(dummy,dummy) > TOpt::dCut[1]) { helix = true; goto do_extrap;} // there is a field on the way
    xx+=(xlast-xfirst)/nstp;
  } 

 do_extrap:

  if(!helix){ // straight line case
    dx=Hout.Hpar[0]-Hpar[0];
    Hout.Hpar[1]=Hpar[1]+Hpar[3]*dx;
    Hout.Hpar[2]=Hpar[2]+Hpar[4]*dx;
    Hout.Hpar[3]=Hpar[3];
    Hout.Hpar[4]=Hpar[4];
    Hout.Hpar[5]=Hpar[5];
    
    if(Hcov[0]==0) goto end; // only track. par. extrap.

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

  } else {    // helix case

    // cut on  P < 100 MeV
    if(fabs(1./Hpar[5]) < 0.100) {
      //cout<<"THlx::Extrapolate() ==> too low momentum for Runge-Kutta propagation : "
      //<<fabs(1./Hpar[5]) <<" GeV"<<endl;
      return(false);
    }

    //Extrapolate state vector by Runge Kutta

    ier=TAlgo::Rkutex(this->Hpar, Hout.Hpar);
    if(ier != 0 ) {
      cout<<"THlx::Extrapolate() ==> (vec.) TAlgo::Rkutex error # "<<ier;
      cout<<"   Momentum = "<<1./Hpar[5]<<"  Gev"<<endl;
      return(false);
    }

    if(Hcov[0]==0) goto end; // only track. par. extrap.
    
    bool flg_sav(false);
    if(TDisplay::Ptr() != NULL){ // if TDisplay object exists
      flg_sav=TDisplay::Ref().Rkutraj();
      TDisplay::Ref().SetRkutraj(false);       // never draw traj. in Jacobian calculation
    }
    //Numerical calculation of Jacobian f[5][5] for helix propagation in magnetic field

    for(i = 1; i <= 2; i++){
      for(j = 1; j <= 5; j++){
	if(i==j){
	  f[j-1][i-1]=1.;
	} else {
	  f[j-1][i-1]=0.;
	}
      }
    }
    for(k = 0; k < 6; k++) vin[k]=Hpar[k];
    for(i = 3; i <= 5; i++){
      vin[i]+=dvin[i];         // increment i-th component
      vout[0]=xlast;
      ier=TAlgo::Rkutex(vin,vout);  // extrapolate
      if(ier!=0) {
	cout<<"THlx::Extrapolate() ==> (cov.) TAlgo::Rkutex error # "<<ier;
	cout<<"   Momentum = "<<1./Hpar[5]<<" Gev"<<endl;
	return(false);
      }
      for(j = 1; j <= 5; j++){
	f[j-1][i-1]=(vout[j]-Hout.Hpar[j])/dvin[i]; // jacobian element dVout(j)/dVin(i)
      }
      vin[i]=Hpar[i];          // Restore i-th component
    }
    f[4][0]=0.; f[4][1]=0.; f[4][2]=0.; f[4][3]=0.; f[4][4]=1.;

    if(TDisplay::Ptr() != NULL){ // if TDisplay object exists
      TDisplay::Ref().SetRkutraj(flg_sav);  // reset prev. status of this flag
    }
/*
    cout<<endl;
    cout<<endl;
    cout<<" dY1/dY0 ="<< f[0][0]<<endl;
    cout<<" dZ1/dY0 ="<< f[1][0]<<endl;
    cout<<" dYp/dY0 ="<< f[2][0]<<endl;
    cout<<" dZp/dY0 ="<< f[3][0]<<endl;
    cout<<" dPi/dY0 ="<< f[4][0]<<endl;
    cout<<endl;
    cout<<" dY1/dZ0 ="<< f[0][1]<<endl;
    cout<<" dZ1/dZ0 ="<< f[1][1]<<endl;
    cout<<" dYp/dZ0 ="<< f[2][1]<<endl;
    cout<<" dZp/dZ0 ="<< f[3][1]<<endl;
    cout<<" dPi/dZ0 ="<< f[4][1]<<endl;
*/

    // cout<<"dY/dPinv = "<<f[0][4]<<" dZ/dPinv = "<<f[1][4]<<endl;

    // Propagate Cov matrix (F*Cov*F.t)

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

    // as dPinv/dx = 0 if x != Pinv
    Hout.Hcov[10]=Hcov[10]*f[0][0] + Hcov[11]*f[0][1] + Hcov[12]*f[0][2] + Hcov[13]*f[0][3] + Hcov[14]*f[0][4];
    Hout.Hcov[11]=Hcov[10]*f[1][0] + Hcov[11]*f[1][1] + Hcov[12]*f[1][2] + Hcov[13]*f[1][3] + Hcov[14]*f[1][4];
    Hout.Hcov[12]=Hcov[10]*f[2][0] + Hcov[11]*f[2][1] + Hcov[12]*f[2][2] + Hcov[13]*f[2][3] + Hcov[14]*f[2][4];
    Hout.Hcov[13]=Hcov[10]*f[3][0] + Hcov[11]*f[3][1] + Hcov[12]*f[3][2] + Hcov[13]*f[3][3] + Hcov[14]*f[3][4];

    Hout.Hcov[14]=Hcov[14];

/*
    // The same with use of TMtx class (for test only)

    TMtx F(5,5), C(5,5), C1(5,5), Dum(5), Par(5);
    double x0;
    this->Get(x0,Dum,C);

    for(i=1; i<=5; i++){
      for(j=1; j<=5; j++){
	F(i,j)=f[i-1][j-1];
      }
    }

    C1 = F*C*F.t(); 

    Hout.Get(x0,Par,C); // for x0, Par

//    cout<<"x0 = "<<x0<<endl;
//    Par.Print("Par");
//    C.Print("C mtx");
//    F.Print("Jacobian F");
//    C1.Print("C1 =F*C*F.t()");

    Hout.Set(x0,Par,C1);
//    Hout.Print("Hout"); cout<<endl;

*/

  } // endif (straight line or helix)

  if(Hout.Hcov[0] < 0 || Hout.Hcov[2] < 0 ||
     Hout.Hcov[5] < 0 || Hout.Hcov[9] < 0 ||
     Hout.Hcov[14]< 0){
    cout<<"THlx::Extrapolate() ==> output cov. matrix is wrong "<<endl;
    this->Print("Input helix");
    Hout.Print("Output helix");
    return(false);
  }

 end:

  /*    
  cout<<"Jacob. old"<<endl;
  cout.setf(ios::showpos);
  cout<<setprecision(3)<<setw(10);
  cout.setf(ios::scientific);
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      cout<<f[i][j]<<"\t ";
    }
    cout<<endl;
  }
  cout.setf(ios::fixed,ios::scientific);
  cout.unsetf(ios::showpos);
  cout<<setprecision(7);
  */

  Hout.path = this->Dist(Hout); //just an estimation of the path length.
  Hout.radLenFr = 0;

  return(true);

}

















