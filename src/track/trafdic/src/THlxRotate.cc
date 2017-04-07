
/*!

  Rotate helix to reference system rotated by matrix: 

  \verbatim  
     cos(a)  sin(a)    0      0     0
    -sin(a)  cos(a)    0      0     0
       0      0     cos(a) sin(a)   0
       0      0    -sin(a) cos(a)   0
       0      0        0      0     1
  \endverbatim 

  \param ca   cos(a) 
  \param sa   sin(a)
   
  \return Hout - rotated helix

*/

#include "THlx.h"
#include <iostream>

void THlx::Rotate(double ca, double sa, THlx& Hout)
{
  Hout.Hpar[0]=Hpar[0];

//  H' = R*H

  Hout.Hpar[1] = ca * Hpar[1] + Hpar[2] * sa;
  Hout.Hpar[2] = ca * Hpar[2] - Hpar[1] * sa;
  Hout.Hpar[3] = ca * Hpar[3] + Hpar[4] * sa;
  Hout.Hpar[4] = ca * Hpar[4] - Hpar[3] * sa;
  Hout.Hpar[5] = Hpar[5];

//  C' = R*C*R.t

  Hout.Hcov[0] = ca*(ca * Hcov[0] + Hcov[1] * sa) + sa*(ca * Hcov[1] + Hcov[2] * sa);  
  Hout.Hcov[1] = ca*(ca * Hcov[1] - Hcov[0] * sa) + sa*(ca * Hcov[2] - Hcov[1] * sa);  
  Hout.Hcov[2] =-sa*(ca * Hcov[1] - Hcov[0] * sa) + ca*(ca * Hcov[2] - Hcov[1] * sa);  
  Hout.Hcov[3] = ca*(ca * Hcov[3] + Hcov[6] * sa) + sa*(ca * Hcov[4] + Hcov[7] * sa);  
  Hout.Hcov[4] =-sa*(ca * Hcov[3] + Hcov[6] * sa) + ca*(ca * Hcov[4] + Hcov[7] * sa);  
  Hout.Hcov[5] = ca*(ca * Hcov[5] + Hcov[8] * sa) + sa*(ca * Hcov[8] + Hcov[9] * sa);  
  Hout.Hcov[6] = ca*(ca * Hcov[6] - Hcov[3] * sa) + sa*(ca * Hcov[7] - Hcov[4] * sa);  
  Hout.Hcov[7] =-sa*(ca * Hcov[6] - Hcov[3] * sa) + ca*(ca * Hcov[7] - Hcov[4] * sa);  
  Hout.Hcov[8] = ca*(ca * Hcov[8] - Hcov[5] * sa) + sa*(ca * Hcov[9] - Hcov[8] * sa);  
  Hout.Hcov[9] =-sa*(ca * Hcov[8] - Hcov[5] * sa) + ca*(ca * Hcov[9] - Hcov[8] * sa);  
  Hout.Hcov[10]= ca*Hcov[10] + Hcov[11] * sa;
  Hout.Hcov[11]= ca*Hcov[11] - Hcov[10] * sa;
  Hout.Hcov[12]= ca*Hcov[12] + Hcov[13] * sa;
  Hout.Hcov[13]= ca*Hcov[13] - Hcov[12] * sa;
  Hout.Hcov[14]= Hcov[14];

}
















