/*!

  Member function 
  to calculate Chi2 increment for "this" helix
  by measurment h (object of THit class) 

  \return dChi2 - Chi2 increment
  
*/
#include <iostream>
#include "THlx.h"
#include "THit.h"
#include "TMtx.h"

double THlx::HitChi2(const THit& h)
{


  // Chi2 increment with 2x2 matricies

  double wu,wv;
  double det,ci0,ci1,ci2,cn0,cn1,cn2,g0,g1,g2;
  double dm0,dm1, xx0,xx1, xm0,xm1;
  double dChi2;

  // C.i
  det=Hcov[0]*Hcov[2]-Hcov[1]*Hcov[1];
  ci0= Hcov[2]/det;
  ci1=-Hcov[1]/det;
  ci2= Hcov[0]/det;

  // (C.i+G).i
  det=(ci0+h.G0)*(ci2+h.G2)-(ci1+h.G1)*(ci1+h.G1);
  cn0= (ci2+h.G2)/det;
  cn1=-(ci1+h.G1)/det;
  cn2= (ci0+h.G0)/det;

  // m-x0
  dm0=h.Y-Hpar[1];
  dm1=h.Z-Hpar[2];

  // x-x0 = (C.i+G).i * G * (m-x0)
  xx0 = dm0*(cn0*h.G0 + cn1*h.G1) + dm1*(cn0*h.G1 + cn1*h.G2);
  xx1 = dm0*(cn1*h.G0 + cn2*h.G1) + dm1*(cn1*h.G1 + cn2*h.G2);

  // x-m = (C.i+G).i * C.i * -(m-x0)
  xm0 = -( dm0*(ci0*cn0 + ci1*cn1) + dm1*(ci1*cn0 + ci2*cn1) );
  xm1 = -( dm0*(ci0*cn1 + ci1*cn2) + dm1*(ci1*cn1 + ci2*cn2) );

  // dChi2 = (x-x0)*C.i*(x-x0) + (x-m)*G*(x-m)
  dChi2=xx0*( ci0*xx0 +  ci1*xx1) + xx1*( ci1*xx0 +  ci2*xx1)+
        xm0*(h.G0*xm0 + h.G1*xm1) + xm1*(h.G1*xm0 + h.G2*xm1);

  return(dChi2);

}





