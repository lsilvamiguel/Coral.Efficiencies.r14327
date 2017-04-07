// $Id: THlxHlxChi2.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!

  Member function 
  to calculate Chi2 of difference of "this" helix
  and helix H

  \return Chi2 - Chi2 of difference
  
*/
#include <iostream>
#include "THlx.h"
#include "THit.h"
#include "TMtx.h"
#include "TConstants.h"

using namespace std;

double THlx::HlxChi2(const THlx& H)
{
  TMtx X0(5), X1(5), DX(5), C0(5,5), C1(5,5), C(5,5);
  double x0, x1, chi2;
  int ierr=0;

  this->Get(x0,X0,C0); // unpack this helix to vector and matrix
      H.Get(x1,X1,C1); // unpack helix H to vector and matrix

  if(fabs(x0-x1) > TConstants_RKuttaMinStep){
    cout<<"THlxHlxChi2() ==> Helices are at different X. this->Hpar[0] = "<<Hpar[0]
	<<" H.Hpar[0] = "<<H.Hpar[0]<<endl;
    assert(false);
  }

  DX = X1-X0;
  C  = C0+C1;
  
  //DX.Print("DX");
  //C. Print("C");

  chi2 = DX.t() * (C.i5(ierr)) * DX;

  if(ierr != 0){
    cout<<"THlxHlxChi2() ==> Error in matrix inversion"<<endl;
    (C0+C1).Print();
    assert(false);
  }

  return(chi2);

}





