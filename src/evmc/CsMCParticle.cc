// $Id: CsMCParticle.cc,v 1.9 2000/07/31 12:08:20 zvyagin Exp $

/*!
   \file    CsMCParticle.cc
   \brief   Compass Monte Carlo particle Class.
   \author  Benigno Gobbo
   \version $Revision: 1.9 $
   \date    $Date: 2000/07/31 12:08:20 $
*/

#include <string.h>
#include <strings.h>
#if defined(__SUNPRO_CC)
# include <ospace/type/prim.h>
#else
# if defined(_AIX) && defined(__cplusplus) && !defined(__GNUC__)
#  include <ospace/type/prim.h>
# endif
#endif

#include "pdgParticles.h"
#include "GeantParticles.h"
#include "CsMCParticle.h"
#include "CsErrLog.h"

CsMCParticle::CsMCParticle() {
  GeantNumber_ = 1;
  bool found = false;
  int pnumber = Gpar2PDGpar( GeantNumber_ );
  for( int i=0; i<nParticles && !found; i++  ) {
    if( PDGpart[i].number == pnumber ) {
      strcpy( name_, PDGpart[0].name );
      charge_      = PDGpart[0].charge;
      mass_        = PDGpart[0].mass;
      massErr_[0]  = PDGpart[0].emass[0];
      massErr_[1]  = PDGpart[0].emass[1];
      width_       = PDGpart[0].width;
      widthErr_[0] = PDGpart[0].ewidth[0];
      widthErr_[1] = PDGpart[0].ewidth[1];
      number_      = PDGpart[0].number;
    }
  }
}


CsMCParticle::CsMCParticle( int gnumber, int isGeantNumber  ) {
    GeantNumber_ = gnumber;
    int pnumber;
    if (isGeantNumber==1) // COMGEANT
        pnumber = Gpar2PDGpar( GeantNumber_ );
    else // TGEANT
    {
        GeantNumber_ = PDGpar2GPar2(gnumber);
        pnumber = gnumber;
    }
    bool found = false;
    if( pnumber < 10000000 ) {
        if( pnumber < 0 ) {
            antiparticle_ = true;
            pnumber = -pnumber;
        }
        else {
            antiparticle_ = false;
        }
        for( int i=0; i<nParticles && !found; i++  ) {
            if( PDGpart[i].number == pnumber ) {
                found = true;
                strcpy( name_, PDGpart[i].name );
                charge_      = PDGpart[i].charge;
                if( antiparticle_ ) charge_ = -charge_;
                mass_        = PDGpart[i].mass;
                massErr_[0]  = PDGpart[i].emass[0];
                massErr_[1]  = PDGpart[i].emass[1];
                width_       = PDGpart[i].width;
                widthErr_[0] = PDGpart[i].ewidth[0];
                widthErr_[1] = PDGpart[i].ewidth[1];
                number_      = PDGpart[i].number;
            }
        }
    }
    else {
        found = true;
        strcpy( name_, GeantPart[GeantNumber_-1].name );
        charge_      = GeantPart[GeantNumber_-1].charge;
        mass_        = GeantPart[GeantNumber_-1].mass;
        massErr_[0]  = 0.0;
        massErr_[1]  = 0.0;
        width_       = GeantPart[GeantNumber_-1].width;
        widthErr_[0] = 0.0;
        widthErr_[1] = 0.0;
        number_      = pnumber;
    }
    if( !found || strlen(name_) == 0 ) {
        name_[0]     = '\0';
        charge_      = 0;
        mass_        = 0.0;
        massErr_[0]  = 0.0;
        massErr_[1]  = 0.0;
        width_       = 0.0;
        widthErr_[0] = 0.0;
        widthErr_[1] = 0.0;
        number_      = 0;
        GeantNumber_ = 0;
//    CsErrLog::Instance()->mes( elWarning, "No particles with this number were found. A dummy particle is being created" );
//    cerr << "CsMCParticle::CsMCParticle():  No particles with this number were found. A dummy particle is being created\n";
    }
}

CsMCParticle::CsMCParticle( const char* pname, int pcharge )
{
  bool found = false;
  for( int i=0; i<nParticles && !found; i++  ) {
    if( strncmp( pname, PDGpart[i].name, strlen(PDGpart[i].name) ) == 0 ) {
      if( PDGpart[i].charge == pcharge ) {
	found = true;
	antiparticle_ = false;
      }
      else if( PDGpart[i].charge == -pcharge ){
	found = true;
	antiparticle_ = true;
      }
      if( found ) {
        strcpy( name_, PDGpart[i].name );
	charge_      = PDGpart[i].charge;
	if( antiparticle_ ) charge_ = -charge_;
	mass_        = PDGpart[i].mass;
	massErr_[0]  = PDGpart[i].emass[0];
	massErr_[1]  = PDGpart[i].emass[1];
	width_       = PDGpart[i].width;
	widthErr_[0] = PDGpart[i].ewidth[0];
	widthErr_[1] = PDGpart[i].ewidth[1];
	number_      = PDGpart[i].number;
	GeantNumber_ = 0;
	for( int j=1; j<92; j++ ) {
	  if( Gpar2PDGpar(j) == number_ ) GeantNumber_ = j;
	}
      }
    }
  }
  if( !found ) {
    name_[0]     = '\0';
    charge_      = 0;
    mass_        = 0.0;
    massErr_[0]  = 0.0;
    massErr_[1]  = 0.0;
    width_       = 0.0;
    widthErr_[0] = 0.0;
    widthErr_[1] = 0.0;
    number_      = 0;
    GeantNumber_ = 0;
//    cerr << "CsMCParticle::CsMCParticle(): No particles with this name and charge were found. A dummy particle is being created\n";
  }
}

CsMCParticle::CsMCParticle( int pnumber, bool antiparticle ) {
  bool found = false;
  for( int i=0; i<nParticles && !found; i++  ) {
    if( pnumber == PDGpart[i].number ) {
      found = true;
      antiparticle_ = antiparticle;
      strcpy( name_, PDGpart[i].name );
      charge_      = PDGpart[i].charge;
      if( antiparticle_ ) charge_ = -charge_;
      mass_        = PDGpart[i].mass;
      massErr_[0]  = PDGpart[i].emass[0];
      massErr_[1]  = PDGpart[i].emass[1];
      width_       = PDGpart[i].width;
      widthErr_[0] = PDGpart[i].ewidth[0];
      widthErr_[1] = PDGpart[i].ewidth[1];
      number_      = PDGpart[i].number;
      GeantNumber_ = 0;
      for( int j=1; j<92; j++ ) {
	if( Gpar2PDGpar(j) == number_ ) GeantNumber_ = j;
      }
    }
  }
  if( !found ) {
    name_[0]     = '\0';
    charge_      = 0;
    mass_        = 0.0;
    massErr_[0]  = 0.0;
    massErr_[1]  = 0.0;
    width_       = 0.0;
    widthErr_[0] = 0.0;
    widthErr_[1] = 0.0;
    number_      = 0;
    GeantNumber_ = 0;
//    cerr << "CsMCParticle::CsMCParticle(): No particles with this name and charge were found. A dummy particle is being created\n";
  }
}


CsMCParticle::CsMCParticle( const CsMCParticle& particle ) {
  strcpy( name_, particle.name_ );
  charge_      = particle.charge_;
  mass_        = particle.mass_;
  massErr_[0]  = particle.massErr_[0];
  massErr_[1]  = particle.massErr_[1];
  width_       = particle.width_;
  widthErr_[0] = particle.widthErr_[0];
  widthErr_[1] = particle.widthErr_[1];
  number_      = particle.number_;
  GeantNumber_ = particle.GeantNumber_;
}


char*  CsMCParticle::getName() {
  return( name_ );
}

int    CsMCParticle::getCharge() {
  return( charge_ );
}

double CsMCParticle::getMass() {
  return( mass_ );
}

double CsMCParticle::getMassErrP() {
  return( massErr_[1] );
}

double CsMCParticle::getMassErrN() {
  return( massErr_[1] );
}

double CsMCParticle::getWidth() {
  return( width_ );
}

double CsMCParticle::getWidthErrP() {
  return( widthErr_[0] );
}

double CsMCParticle::getWidthErrN() {
  return( widthErr_[1] );
}

int CsMCParticle::getNumber() {
  return( number_ );
}

int CsMCParticle::getGeantNumber() {
  return( GeantNumber_ );
}

bool CsMCParticle::operator==( const CsMCParticle& part ) const {
  if( number_ == part.number_ && 
      antiparticle_ == part.antiparticle_ )
    return( true );
  else
    return( false );
}

int CsMCParticle::Gpar2PDGpar( int Gpar ) {

  //                         1          2          3          4          5
  //                     gamma         e+         e-        nue        mu+
  int PDGpar[92] = {0,      22,       -11,        11,        12,       -13,
		    //       6          7          8          9         10
		    //     mu-        pi0        pi+        pi-        K0l
		            13,       111,       211,      -211,       130,
		    //      11         12         13         14         15
		    //      K+         K-          N          P      antiP
		           321,      -321,      2112,      2212,     -2212,
		    //      16         17         18         19         20
		    //     K0s        Eta     Lambda     Sigma+     Sigma0 
		           310,       221,      3122,      3222,      3212,
		    //      21         22         23         24         25
		    //  sigma-        Xi0        Xi-     omega-      antiN
  		          3112,      3322,      3312,       227,     -2112,
		    //      26         27         28         29         30
		    //anLambda   anSigma-   anSigma0   anSigma+    antiXi0
		         -3122,     -3222,     -3212,     -3112,     -3322,
		    //      31         32         33         34         35
		    // antiXi+  antiOmega          -          -         D+
		         -3312,      -227,  10000033,  10000034,       411,
		    //      36         37         38         39         40
		    //      D+         D0     antiD0        Ds+        Ds-
                          -411,       421,      -421,       431,      -431,
		    //      41         42         43         44         45
		    //    K0sd      Lam0d    Sig-ndh       Pnhi   Deuteron  
		           310,      3122,      3112,      2212,  10000045,
		    //      46         47         48         49         50
		    //  Triton      Alpha   Geantino        He3   Cerenkov
		      10000046,  10000047,  10000048,  10000049,  10000050,
		    //      51         52         53         54         55
		    //   Xi-sd     Ome-sd      Sig+2      Sig-2        D+* 
		          3122,       227,      3222,      3112,       413,
		    //      56         57         58         59         60
		    //     D-*        D0*    antiD0*       Ds+*       Ds-*
		          -413,       423,      -423,       433,      -433,
		    //      61         62         63         64         65
		    //       -          -          -          -          -
  		      10000061,  10000062,  10000063,  10000064,  10000065,
		    //      66         67         68         69         70
		    //       -          -          -          -          -
  		      10000066,  10000067,  10000068,  10000069,  10000070,
		    //      71         72         73         74         75
		    //   LamC+       XiC+       XiC0      OmeC0     SigC++
		          4122,      4232,      4132,      4332,      4222,
 		    //      76         77         78         79         80
		    //   SigC+      SigC0          -          -          -
		          4212,      4112,  10000078,  10000079,  10000080,
		    //      81         82         83         84         85
		    // anLamC+          -          -          -          -
                         -4122,  10000082,  10000083,  10000084,  10000085,
		    //      86         87         88         89         90
		    //       -          -          -          -          -
 		      10000086,  10000087,  10000088,  10000089,  10000090,
		    //      91
		    //     Phi
                           333 };
  return( PDGpar[Gpar] );
}

int CsMCParticle::PDGpar2GPar2( int PDGpar)
{
  for (int i =1; i < 92; i ++){
    if (PDGpar == Gpar2PDGpar(i))
      return i;
  }
  return -1;   
}



 
