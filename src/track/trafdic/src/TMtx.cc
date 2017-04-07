// $Id: TMtx.cc 13148 2011-12-28 16:55:25Z kbicker $

#include <iostream>
#include <iomanip>
#include "TMtx.h"
#include "TAlgo.h"

using namespace std;

/*!  
  \param str comment to be printed
*/
void TMtx::Print(const char* str) const
{
  if (m==1||n==1) cout<<str<<": vector("<<m<<"x"<<n<<")";
  else            cout<<str<<": matrix("<<m<<"x"<<n<<")";
  cout.setf(ios::scientific);
  cout.setf(ios::showpos);
  for(int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      cout<<setprecision(3)<<setw(13)<<Mtx[i][j];
    }
    cout<<endl;
  }
  //Reset some format flags
  cout.setf(ios::fixed,ios::scientific);
  cout.unsetf(ios::showpos);
  cout<<setprecision(7);
}

//
// Transpose matrix
//

TMtx TMtx::t()
{
  TMtx M(n,m);
  for( int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      M.Mtx[j][i]=Mtx[i][j];
    }
  };
  return(M);
}

//
// Invert symmetric matrix
//

extern "C" void dsinv_(int*, double*, int*, int*); // CERNLIB

TMtx TMtx::invs(int& ierr)
{
  ierr=0;

  if(m!=n){
    cout<<"TMtx::i() ==> Matrix ("<<m<<"x"<<n<<") - not square!"<<endl;
    assert(false);
  }

  double A[n][m];
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++) A[j][i]=Mtx[i][j];
  }

  dsinv_(&m, &A[0][0], &n, &ierr);

  TMtx Res(m,n);

  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++) Res.Mtx[i][j]=A[j][i];
  }

  return(Res);
}
//
// Invert any matrix
//

extern "C" void dinv_(int*, double*, int*, float*,int*);  // CERNLIB

TMtx TMtx::i(int& ierr)
{
  ierr=0;

  if(m!=n){
    cout<<"TMtx::inv() ==> Matrix ("<<m<<"x"<<n<<") - not square!"<<endl;
    assert(false);
  }

  double A[n][m];
  float r[m];
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++) A[j][i]=Mtx[i][j];
  }

  dinv_(&m, &A[0][0], &n, r, &ierr);

  TMtx Res(m,n);

  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++) Res.Mtx[i][j]=A[j][i];
  }

  return(Res);
}

//
// Invert 5x5 matrix
//

TMtx TMtx::i5(int& ierr)
{
  if(m!=5 || n!=5){
    cout<<"TMtx::i5() ==> Matrix ("<<m<<"x"<<n<<") - not (5x5) !"<<endl;
    assert(false);
  }

  ierr=0;
  double C[15], Cinv[15];

  for(int i=0; i<m; i++){
    for(int j=0; j<=i; j++)
      C[j+((i+1)*i)/2]=Mtx[i][j];
  }

  ierr=TAlgo::Inv5(C,Cinv);

  TMtx Res(m,n);

  for(int i=0; i<m; i++){
    for(int j=0; j<=i; j++)
      Res.Mtx[i][j]=Res.Mtx[j][i]=Cinv[j+((i+1)*i)/2];
  }
  return(Res);
}

//
// Overloaded operators:
//
//  "="
//
TMtx& TMtx::operator = (const TMtx &M)
{

  if((M.m!=m)||(M.n!=n)){
    cout<<"TMtx::operator = ()  ==> Matrix dimentions mismatch : "
	<<"("<<m<<"x"<<n<<")=("<<M.m<<"x"<<M.n<<") !"
	<<endl; 
    assert(false);
  };

  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++)
      this->Mtx[i][j]=M.Mtx[i][j];
  };
  return(*this);
}

TMtx::operator double()
{
  if((m != 1)||(n != 1)){
    cout<<"TMtx::operator double() ==> conversion from  M("<<m<<"x"<<n<<") matrix to 'double'"<<endl;
    assert(false);
  };
  return(Mtx[0][0]);
}

//
//  "()"    
//
double& TMtx::operator () (const int i, const int j) const
{

  if((i<=0) || (j<=0) || (i>m)||(j>n)){
    cout<<"TMtx:: operator() ==> try to access element "
      <<i<<","<<j
	<<"  of matrix ("<<m<<"x"<<n<<")"
	  <<endl;
    assert(false);
  };
  return(Mtx[i-1][j-1]);
}

//      
// "*"
//
TMtx TMtx::operator * (const TMtx &M)
{

  if(n != M.m){
    cout<<"TMtx:: operator * "
      <<"("<<m<<"x"<<n<<")*("<<M.m<<"x"<<M.n<<") !"
	<<endl; 
    assert(false);
  };

  TMtx Res(m,M.n);

  for (int i=0; i<m; i++){
    for (int j=0; j<M.n; j++){
      Res.Mtx[i][j]=0;
      for(int k=0; k<n; k++){
	Res.Mtx[i][j] += Mtx[i][k]*M.Mtx[k][j];
      };
    };
  };
  return(Res);
}
//
// "+"
//
TMtx TMtx::operator + (const TMtx &M)
{

  if((n!= M.n)||(m!=M.m)){
    cout<<"TMtx:: operator + "
      <<"("<<m<<"x"<<n<<")+("<<M.m<<"x"<<M.n<<") !"
	<<endl; 
    assert(false);
  };

  TMtx Res(m,n);

  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      Res.Mtx[i][j] = Mtx[i][j] + M.Mtx[i][j];
    };
  };
  return(Res);
}
//
// "-"
//
TMtx TMtx::operator - (const TMtx &M)
{

  if((n!= M.n)||(m!=M.m)){
    cout<<"TMtx:: operator - "
      <<"("<<m<<"x"<<n<<")-("<<M.m<<"x"<<M.n<<") !"
	<<endl; 
    assert(false);
  };

  TMtx Res(m,n);

  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      Res.Mtx[i][j] = Mtx[i][j] - M.Mtx[i][j];
    };
  };
  return(Res);
}

//
//  "+="
//
TMtx& TMtx::operator += (const TMtx &M)
{

  if((n!= M.n)||(m!=M.m)){
    cout<<"TMtx:: operator + "
      <<"("<<m<<"x"<<n<<")+=("<<M.m<<"x"<<M.n<<") !"
	<<endl; 
    assert(false);
  };

  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      Mtx[i][j] += M.Mtx[i][j];
    };
  };
  return(*this);
}

//
//  "*="
//
TMtx& TMtx::operator *= (const double& w)
{
  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      Mtx[i][j] *= w;
    };
  };
  return(*this);
}
