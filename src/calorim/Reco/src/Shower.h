/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Shower.h,v $
   $Date: 2010/09/20 10:39:05 $
   $Revision: 1.13 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )
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

#ifndef Shower_Reco________include
#define Shower_Reco________include

#include "Calorimeter.h"   // for Calorimeter::Particle::ParticleID

namespace Reco {

/*! \brief Electromagnetic shower profile.

    \return energy deposition from gamma/electron with energy E in given cell
    \arg E - particle (gamma/electron/positron) energy
    \arg X,Y - particle coordinates
    \arg dX,dY - cell size; cell center has coordinates (0,0)
    \arg angleX,angleY - angle (in radians) between projection of particle track on ZX(ZY) plane and XY plane
*/
// ??? typedef double (*ShowerType)(double RadLeng,double E,double X, double Y, double dX, double dY, double angleX=0, double angleY=0);

double Shower(double RadLeng,double E,double X, double Y, double dX, double dY, double angleX=0, double angleY=0);

/// Lednev's comulative function
double FLednev(double x, double y);

double ShowerLednev( double RadLeng,double E,double X, double Y, double dX, double dY, double angleX=0, double angleY=0);

///  Coordinate as a function of energy deposit left and right from boundary, angle and material
double CoordinateFunction(double RadLeng,double el,double er ,double angl);

///  Coordinate as a function of energy deposit left and right from boundary, angle and material for hadronic shower
double HadronicCoordinateFunction(double NuclLeng,double RadLeng,double el,double er ,double angl);

///  Average Z position of Hadronic Shower calculated in nucler lengths from the cell's front surface
double ZmidShowerHadronic(double e);
double SigmaZmidShowerHadronic(double e);

/// deprecated
double ZmidShowerEM(double e);

///  Average Z position of photon shower calculated in radiation lengths from the cell's front surface
double ZmidShowerPhoton(double e, double e_crit);
///  Average Z position of electron shower calculated in radiation lengths from the cell's front surface
double ZmidShowerElectron(double e, double e_crit);

/// Sigma of EM shower Z position in radiation lengths (always returns 2)
double SigmaZmidShowerEM(double e);

/// Fraction of shower energy that is leaked through the calorimeter.
double LongLeakageEM(double e, double e_crit,
                     CalorimeterParticle::ParticleID partid, double thickness);


/*! \brief Coordinate Function for Lead Glass

    \param ar - energy fraction right from boundary
*/
double	 acumulative_LeadGlass(double RadLeng,double ar);
double	 acumulative_LeadGlassAngle(double RadLeng,double ar, double angl, double esum);
/*! \brief Coordinate Function for Hadron Shower

    \param ar - energy fraction right from boundary
*/
double	 acumulative_Hadron(double NuclLeng,double ar);

/*
   unnormalized longitudinal function t**tmax * exp(-t)
*/
double   cascade(double xmax,double x);
/*
   normalized cumulative longitudinal function
*/
double   cumulative_cascade(double xmax,double x);

/*! \brief  Hadronic Shower energy deposit in the range (MC)
    \param E    - energy in GeV
    \param xmin - min in Nuclear Length units
    \param xmax - max in Nuclear Length units
*/
 double	 cascade_Hadron(double E,double xmin,double xmax);

/*! \brief  Electro-Magnetic Shower energy deposit in the range (MC)
    \param E    - energy in GeV
    \param xmin - min in  Radiation Length units
    \param xmax - max in  Radiation Length units
*/
double	 cascade_EM(double E,double xmin,double xmax);

/*! \brief  MIP energy deposit in the range (MC)
    \param E    - energy in GeV
    \param xmin - min in  density(g/cm**3)*x(cm)  units
    \param xmax - max in  density(g/cm**3)*x(cm)  units
*/
double	 cascade_MIP(double E,double xmin,double xmax);

/*! \brief Few simple plane geometry functions
    \param points p11,p12,p21,p22 (x,y)
*/
double	 distmin_LineLine(double *p11,double *p12,double *p21,double *p22);
double	 distmin_PointLine(double *p1,double *p21,double *p22);
double	 distmin_PointPoint(double *p1,double *p2);
double	 distmin_LineRectang(double *p01,double *p02,double *p11,double *p12,double *p21,double *p22);

// /*! \brief Electro-Magnetic or Hadronic cascades parametrisation used in FMC and reconstruction algorithms
//            We use parametrication in scale variables (PDG):
//            For E-M shower
//            --------------
//                       dE/dt = E0
// */
class Shower
{
  // =========================================================================
  // Constructors and destructor
  // =========================================================================

  public:

    /// Destructor
    virtual            ~Shower   (void) {}

    /// Default constructor
    Shower (void) : p1_ ( 1.480 ), p2_ ( 1.600 ), p3_ ( 0.650 ), rad_leng_ ( 23.7 ) {}

  // ==========================================
  // Methods
  // ==========================================

  protected:

    double aa2f( double ax, double ay );
    double a2fz( double xc, double yc, double z );
    double afz(double x, double z);
    double cumprf(double cx,double cl);
    double cumcryex( double xc );
    double acumex( double x );
    double acumulative( double ar);

  // ==========================================
  //  Attributes, data
  // ==========================================

 private:
    double p1_, p2_, p3_, rad_leng_;
};

} // namespace Reco
#endif // Shower_Reco________include
