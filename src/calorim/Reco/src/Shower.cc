/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Shower.cc,v $
   $Date: 2010/06/18 10:44:20 $
   $Revision: 1.22 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su, Vladimir.Kolosov@cern.ch )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

   Copyright(C): 1999-2000  V.Kolosov,A.Zvyagin

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include <cmath>
#include <iostream>

#include "TMath.h"

#include "Exception.h"
#include "Shower.h"

using namespace std;

namespace Reco {

/////////////////////////////////////////////////////////////////////////////////////

// Lednev's shower profile

/*!
    \brief Calculate expected amplitude fraction in cell
    \param x {cell center distance from particle in X direction in units of radiation length}
    \param y {cell center distance from particle in Y direction in units of radiation length}
    \return {fraction of particle enrgy deposite in cell}
    \callgraph
    \callergraph
    \remarks {
    Hardcoded:
    int N : Number of parameters used in Lednev parameterisation
    double a_[N] : Parameter for Lednev function
    double b_[N] : Parameter for Lednev function
    }
    \todo {Clarify why there are 3 parameter sets used here}
  */


double FLednev(double x, double y)
{
  const int N = 3;
  static const double
    a_[N] = { 0.80 , 0.30 , -0.10 },
    b_[N] = { 0.80 , 0.20 ,  7.6  };  // cm

  double result = 0;

  for( int i=0; i<N; i++ )
  {
    const double &a=a_[i], &b=b_[i];
    result += a * ( atan(x/b) + atan(y/b) + atan(x*y/(b*sqrt(b*b+x*x+y*y))) );
  }

  return result/M_PI/2+0.25;
}


/*!
  \brief Calculate expected amplitude for pcell in a cluster
  \param RadLeng {Radiation length of cell}
  \param E {Total energy of particle}
  \param X {cell center distance from particle in X direction}
  \param Y {cell center distance from particle in X direction}
  \param dX {Possible deviation of X}
  \param dY {Possible deviation of Y}
  \param angelX {Particle angel in X direction}
  \param angelY {Particle angel in Y direction}
  \return {Energy deposite in cell}
  \callgraph
  \callergraph
  \remarks {
  Hardcoded:
  bool turbo : Use very simple shower expectations
  }
*/

double ShowerLednev(double RadLeng,double E,double X, double Y, double dX, double dY, double angleX, double angleY)
{
// GoTo cm
//  bool turbo = true;
  bool turbo = false;
  if(turbo)
  {
    if( (abs(X) < dX/2) &&  (abs(Y) < dY/2) ) return E;
    if( (abs(X) < dX) &&  (abs(Y) < dY) ) return E/4;
    return E/100.;
  }
  else
  {

    // ????????????????
    double f = 10.*RadLeng/23.7;
    return
      E * ( FLednev((X+dX/2)/f,(Y+dY/2)/f) -
            FLednev((X+dX/2)/f,(Y-dY/2)/f) -
            FLednev((X-dX/2)/f,(Y+dY/2)/f) +
            FLednev((X-dX/2)/f,(Y-dY/2)/f) );
  }
}

/////////////////////////////////////////////////////////////////////////////////////

double CoordinateFunction(double RadLeng,double el,double er ,double angl)
{
  double esum = el + er;
  if( esum <= 0 )
    throw Exception("Reco::CoordinateFunction():  bad esum = %g",esum);

  double efrl = el/esum;
//  cout << " Eleft " << el << " Eright " << er << " Esum " << esum << endl;
//  cout << " Coordinate " << acumulative_LeadGlass(efrl) << endl;

//  return acumulative_LeadGlass(efrl)/10.; // coordinate in cm
  bool new_function = true;
  if( new_function )
    return acumulative_LeadGlassAngle(RadLeng, efrl, angl, esum); // coordinate in mm
  else
    return acumulative_LeadGlass(RadLeng,efrl); // coordinate in mm
}

/////////////////////////////////////////////////////////////////////////////////////

double HadronicCoordinateFunction(double NuclLeng,double RadLeng,double el,double er ,double angl)
{
  double esum = el + er;
  if( esum <= 0 )
    throw Exception("Reco::HadronicCoordinateFunction():  bad esum = %g",esum);

  double efrl = el/esum;
  return acumulative_Hadron(NuclLeng,efrl); // coordinate in mm
}

/////////////////////////////////////////////////////////////////////////////////////

double	 acumulative_LeadGlass(double RadLeng,double ar)
{
   double p1 = 1.480;
   double p2 = 1.600;
   double p3 = 0.650;

   double a = ar;
   double s = 1;

   if( ar > 0.5 ) {
    a=1.-ar;
    s=-1.;
    }

   if( a <= 0. ) a=0.0001;
   double y=log(2.*a);
   double x=(p2*y-p1)*y;
   x=x/p3;

//   x = x / 2.8;
   x = x/0.93;
   x = x*RadLeng/23.7;
   return s*x;
}

/////////////////////////////////////////////////////////////////////////////////////

double	 acumulative_LeadGlassAngle(double RadLeng, double ar, double angle, double esum )
{
//    bool debug = true;
   bool debug = false;

   double a = ar;
   double s = 1;

   if( ar > 0.5 ) {
    a=1.-ar;
    s=-1.;
    }


   double x_old = 0.;
   {                       //  evaluate coordinate with old formula
     double p1 = 1.480;
     double p2 = 1.600;
     double p3 = 0.650;

     if( a <= 0. ) a=0.0001;
     double y=log(2.*a);
     double x=(p2*y-p1)*y;
     x=x/p3;

     x = x/0.93;
     x = x*RadLeng/23.7;
     x_old = x;
   }

   double x_new = 0.;
   {                       //  evaluate coordinate with new formula
     double p0 =0.700;
     double p1 =0.250;

     if( a <= 0. ) a=0.0001;
     double y=log(2.*a);
     double x=(p0*y-1.)*y/p1;

     x = x*RadLeng/23.7;
     x_new = x;
   }

   if( debug && esum > 10.)
   {
     cout << " Ar " << ar << " RadLeng " << RadLeng << " angle " << angle << " esum " << esum << endl;
     cout << " Xold " <<  x_old << " deltaX new-old " << x_new-x_old << endl;
   }

//    double x = x_old;
   double x = x_new;
   return s*x;
}

/////////////////////////////////////////////////////////////////////////////////////

double	 acumulative_Hadron(double NuclLeng,double ar)
{
   double p1 = 1.480;
   double p2 = 1.600;
   double p3 = 0.650;

   double a = ar;
   double s = 1;

   if( ar > 0.5 ) {
    a=1.-ar;
    s=-1.;
    }

   if( a <= 0. ) a=0.0001;
   double y=log(2.*a);
   double x=(p2*y-p1)*y;
   x=x/p3;

//   x = x / 2.8;
   x = x/0.93;
   x = x/0.6;
   x = x*NuclLeng/23.7;
   return s*x;
}

/////////////////////////////////////////////////////////////////////////////////////
/*
   unnormalized longitudinal function t**tmax * exp(-t)
*/
double cascade(double xmax,double x)
{
  if(x <= 0) return 0;
  if(xmax <= 0)
  {
    cerr << " Unexpected argument in cascade function xmax=" << xmax << endl;
    return 0;
  }
  double lnt = log(x);
  return exp(lnt*xmax-x);

}

/////////////////////////////////////////////////////////////////////////////////////
/*
   normalized cumulative longitudinal function
*/
double cumulative_cascade(double xmax,double x)
{
  if(x <= 0) return 0;
  if(xmax <= 0)
  {
    cerr << " Unexpected argument in cascade function xmax=" << xmax << endl;
    return 0;
  }

  double norm=0;

  for(int i=0; i<400; i++)
  {
    double xi=(i+0.5)*0.1;
    norm += cascade(xmax,xi)*0.1;
  }

  int i2=(int)(x/0.1);
  if( i2 >= 400) return 1;

  double integral=0;
  for(int i=0; i<i2; i++)
  {
    double xi=(i+0.5)*0.1;
    integral += cascade(xmax,xi)*0.1;
  }

  return integral/norm;

}

/////////////////////////////////////////////////////////////////////////////////////

double	 cascade_Hadron(double E,double xmin,double xmax)
{
  return E*(cumulative_cascade(7,xmax)-cumulative_cascade(7,xmin));
}

/////////////////////////////////////////////////////////////////////////////////////

double ZmidShowerHadronic(double e)
{
  return 1.75;
}

/////////////////////////////////////////////////////////////////////////////////////

double SigmaZmidShowerHadronic(double e)
{
  return 1.;
}

/////////////////////////////////////////////////////////////////////////////////////

double ZmidShowerEM(double e) {
  return ZmidShowerElectron(e, .025);
}

/////////////////////////////////////////////////////////////////////////////////////
/// Weighted mean z position of shower energy deposition.  Obtained by
/// integral from 0 to infinity over t*dE/dt divided by integral from 0 to
/// infinity over dE/dt, with dE/dt taken from PDG 2008 Eq. 27.31, used
/// together with Eq. 27.32 and b assumed constant: b=0.5
/// \param[in] e       shower energy
/// \param[in] e_crit  critical energy

double ZmidShowerElectron(double e, double e_crit)
{
  return log(e/e_crit) + 1.5;
}


/////////////////////////////////////////////////////////////////////////////////////
/// Weighted mean z position of shower energy deposition.  Obtained by
/// integral from 0 to infinity over t*dE/dt divided by integral from 0 to
/// infinity over dE/dt, with dE/dt taken from PDG 2008 Eq. 27.31, used
/// together with Eq. 27.32 and b assumed constant: b=0.5
/// \param[in] e       shower energy
/// \param[in] e_crit  critical energy

double ZmidShowerPhoton(double e, double e_crit)
{
  return log(e/e_crit) + 2.5;
}


/////////////////////////////////////////////////////////////////////////////////////
/// \todo check this
double SigmaZmidShowerEM(double e)
{
  return 2.;
}


/////////////////////////////////////////////////////////////////////////////////////
/// Obtained by integration over Eq. 27.31 of PDG 2008, used together with
/// Eq. 27.32 and with b assumed constant: b=0.5
/// \param[in] e          shower energy
/// \param[in] e_crit     critical energy
/// \param[in] partid     particle id for the shower particle
/// \param[in] thickness  thickness of calorimeter in units of rad.len.

double LongLeakageEM(double e, double e_crit,
                     CalorimeterParticle::ParticleID partid, double thickness)
{
  // assume that shower is fully contained within this amount of radlens
  const double maxthick = 1e30;
  double c;

  if (    partid == CalorimeterParticle::POSITRON
       || partid == CalorimeterParticle::ELECTRON )
    c = 0.5;
  else if ( partid == CalorimeterParticle::GAMMA )
    c = -0.5;
  else
    throw Exception("LongLeakage(): partid %i not implemented.", partid);

  const double p = ( 2. + c + log(e/e_crit) ) / 2.;

  return 1. - TMath::Gamma(p, 0.5*thickness) / TMath::Gamma(p, 0.5*maxthick);
}


/////////////////////////////////////////////////////////////////////////////////////

double	 cascade_EM(double E,double xmin,double xmax)
{
//         er=0.025
// c for electron ?
// c	tmx=alog(e/er)-0.5
// c for gamma?
// 	tmx=alog(e/er)+0.5
//
// 	t=z/rl
// c	write (*,*) ' tmx ',tmx
// 	x=t/2.
// 	xm=tm/2.
// 	if(x.gt.40.) then
// 		azlnt=0.
// 	else
// 		azlnt=exp(xm*alog(x))*exp(-x)
// 	endif

//  double tmax = 7.;
  if(xmin < 0 ) xmin=0;
  if(xmax < 0 ) xmax=0;
  double tmax = log(E/0.025)+0.5;
  return E*(cumulative_cascade(tmax,xmax)-cumulative_cascade(tmax,xmin));
}

/////////////////////////////////////////////////////////////////////////////////////

double	 cascade_MIP(double E,double xmin,double xmax)
{
  double dx=(xmax-xmin)*0.0005;
  if(dx <=0) return 0;
  else return dx;
}

/////////////////////////////////////////////////////////////////////////////////////

double	distmin_LineRectang(double *p01,double *p02,double *p11,double *p12,double *p21,double *p22)
{
  double dmin=distmin_LineLine(p01,p02,p11,p12);
  double d=distmin_LineLine(p01,p02,p12,p21);
  if(d < dmin) dmin=d;
  d=distmin_LineLine(p01,p02,p21,p22);
  if(d < dmin) dmin=d;
  d=distmin_LineLine(p01,p02,p22,p12);
  if(d < dmin) dmin=d;

  return dmin;
}

/////////////////////////////////////////////////////////////////////////////////////

double	distmin_LineLine(double *p11,double *p12,double *p21,double *p22)
{
  double dmin=distmin_PointLine(p11,p21,p22);
  double d=distmin_PointLine(p12,p21,p22);
  if(d < dmin) dmin=d;
  d=distmin_PointLine(p21,p11,p12);
  if(d < dmin) dmin=d;
  d=distmin_PointLine(p22,p11,p12);
  if(d < dmin) dmin=d;
  return dmin;
}

/////////////////////////////////////////////////////////////////////////////////////

double	distmin_PointLine(double *p1,double *p21,double *p22)
{
  double dx = p22[0]-p21[0];
  double dy = p22[1]-p21[1];
  if(fabs(dx) < 0.00001 && fabs(dy) < 0.00001 ) return distmin_PointPoint(p1,p21);
  double lengthsq=dx*dx+dy*dy;
  double alfa = ((p1[0]-p21[0])*dx+(p1[1]-p21[1])*dy)/lengthsq;
  double x_cross[2];
  x_cross[0]=p21[0]+alfa*dx;
  x_cross[1]=p21[1]+alfa*dy;
  if( (x_cross[0]-p21[0])*(p22[0]-x_cross[0]) >= 0 && (x_cross[1]-p21[1])*(p22[1]-x_cross[1]) >= 0 )
  {
    return distmin_PointPoint(p1,x_cross);
  }
  else
  {
    double d1=distmin_PointPoint(p1,p21);
    double d2=distmin_PointPoint(p1,p22);
    if(d2 < d1) return d2;
    else return d1;
  }
}

/////////////////////////////////////////////////////////////////////////////////////

double	distmin_PointPoint(double *p1,double *p2)
{
  double dx = p1[0]-p2[0];
  double dy = p1[1]-p2[1];
  return sqrt( dx*dx+dy*dy );
}

/////////////////////////////////////////////////////////////////////////////////////

double  Shower::aa2f( double ax, double ay )
{
  const static double c0=0.82;
  const static double c1=-2.2;
  const static double c2=-24.30;
  double x = 0.5 - ax;
  double y = 0.5 - ay;
  double x2 = 2.*x;
  double y2 = 2.*y;
  double pnx = x*(1. - x2*x2 );
  double pny = y*(1. - y2*y2 );
  double pxy = c0+c1*(x*x+y*y)+c2*x*x*y*y;
  double aa2 = ax*ay+pnx*pny*pxy;
  return aa2;
}

double  Shower::a2fz( double xc, double yc, double z )
{
  double ax=afz(xc,z);
  double ay=afz(yc,z);
  return aa2f(ax,ay);
}

double Shower::afz(double x, double z)
{
  if( x > 0. )
    return (1.-cumprf(x,z));
  else
    return cumprf(-x,z);
}

double Shower::cumprf(double cx,double cl)
{
  return cumcryex(cx);
}

double Shower::cumcryex( double xc )
{
  double x = fabs( xc );
  double a = acumex(x);
  if( xc > 0 )
    return a;
  else
    return 1.-a;
}

double Shower::acumex( double x )
{
//   const static  double p1 = 1.480;
//   const static  double p2 = 1.600;
//   const static  double p3 = 0.650;
  double a = 0.5*exp((p1_-sqrt(p1_*p1_+4.*p2_*x))/2./p2_);
  return a;
}

double	 Shower::acumulative( double ar)
{
//    double p1 = 1.480;
//    double p2 = 1.600;
//    double p3 = 0.650;

   double a = ar;
   double s = 1;

   if( ar > 0.5 ) {
    a=1.-ar;
    s=-1.;
    }

   if( a <= 0. ) a=0.0001;
   double y=log(2.*a);
   double x=(p2_*y-p1_)*y;
   x=x/p3_;

//   x = x / 2.8;
   x = x/0.93;
   x = x*rad_leng_/23.7;
   return s*x;
}

} // namespace Reco
