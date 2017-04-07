// $Id: CsTrafficUtils.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
   \file    CsTrafficUtils.cc
   \brief   Utils for Traffic tracking:
   - genLattice: Generate lattice look-up table == Dico and output to file.
   \author  Yann.Bedfer@cern.ch
   \version $Revision: 13148 $
   \date    $Date: 2011-12-28 17:55:25 +0100 (Wed, 28 Dec 2011) $

*/

#include <stdio.h>
#include "CsErrLog.h"
#include "CsInit.h"
#include "CsGeom.h"
#include "CsTrafficUtils.h"
#include "Traffic.h"
#include "TOpt.h"
#include "TLattice.h"
#include "TSetup.h"

using namespace std;

int GetHitsFromTrack(float *,float *,int,float *);

//Constructor
CsTrafficUtils::CsTrafficUtils() {
  new Traffic; // create Traffic package object
}

//Destructor (empty)
CsTrafficUtils::~CsTrafficUtils() {}

//Method for the generation of the Lattice file == Dico
bool CsTrafficUtils::genLattice() {

  // Instantiate CsInit to access user options.
  CsInit *init = CsInit::Instance();
  unsigned int maxEvents_ = init->getMaxEvents();
  unsigned int firstEvent = init->getSkipEvents();

  // Create lattice grid
  TLattice lat = TLattice(0);
  if (maxEvents_==0) maxEvents_ = lat.getSize(); // #Events = grid size
  else
    CsErrLog::mes(elWarning,"events to read != 0 in options file => Dico may be truncated.");
  if (firstEvent)
    CsErrLog::mes(elWarning,"events to skip != 0 in options file => Dico may be truncated.");

  // Init lattice generation
  lat.GenInit();
  // Write detectors abscissa (so as to document Dico w/ geometry)
  Hits coords;  // As many abscissa as there are hits 
  const TSetup &setup = TSetup::Ref();
  for (int ipl = 0; ipl<(int)setup.vDetect().size(); ipl++) {
    const TDetect &d = setup.vDetect()[ipl];
  }
  double xpl; int ipl; for (ipl = 0, xpl = setup.TargetCenter[0];
			    ipl<(int)setup.vDetect().size(); ipl++) {
    // loop over planes
    double xprv = xpl; xpl = setup.vDetect()[ipl].X(0);
    if (xpl==xprv) continue;         // ...consecutive det's @ same abscissa
    int idx = lat.dico_idx[ipl];
    if (idx==-1) continue;   // ...outside dico
    if (idx==-2) break;      // ...beyond dico
    coords[idx] = xpl;
  }
  lat.GenOutput(coords);

  // loop on events

  for (unsigned int nEvents = firstEvent; nEvents<maxEvents_; nEvents++) {

    if (nEvents%100 == 0 )
      cout<< "Event: "<<nEvents<< endl;

    // Track par's from lattice track server
    float Pxyz[3], Vxyz[3]; int Ipart;
    lat.GenServer(Vxyz,Pxyz,&Ipart);

    // Propagate track and get hits coord's
    int return_val;
    if ((return_val = GetHitsFromTrack(Vxyz,Pxyz,Ipart,coords))!=0) {
      char mes[80];
      if      (return_val==-1) sprintf(mes,"Event %d uncomplete!\n",nEvents);
      else                     sprintf(mes,"Event %d error extrapolating!\n",
				       nEvents);
      CsErrLog::mes(elError,mes);
    }
      

    // OutPut
    lat.GenOutput(coords);
  }
  
  lat.GenClose();         // Close Lattice generation

  return (maxEvents_==lat.getSize() && firstEvent==0);
}

#include "CsMCParticle.h"
#include "TOpt.h"
#include "THlx.h"
#include "TConstants.h"
int GetHitsFromTrack(float *Vxyz,float *Pxyz,int Ipart,float *coords)
{
  const TSetup& setup = TSetup::Ref();
  THlx H0, H1, Hlast;
  double Pmc;

  H0(0) = Vxyz[0]; H0(1) = Vxyz[1]; H0(2) = Vxyz[2];
  H0(3) = Pxyz[1]/Pxyz[0];
  H0(4) = Pxyz[2]/Pxyz[0];
  H0(5) = 1/sqrt(Pxyz[0]*Pxyz[0] +  Pxyz[1]*Pxyz[1] + Pxyz[2]*Pxyz[2]);
  H0(5) *= CsMCParticle(Ipart).getCharge();

  bool error; int ipl; double xpl;
  for (ipl = 0, xpl = setup.TargetCenter[0], error = false;
       ipl<(int)setup.vDetect().size(); ipl++) {

    // *************** LOOP OVER DETECTORS ***************
    const TDetect &d = setup.vDetect(ipl);

    double xprv = xpl; xpl = setup.vDetect()[ipl].X(0);
	
    if (fabs(xpl-xprv)<.0001) continue; // ...consecutive det's @ same abscissa
    const TLattice *lat = TLattice::Ptr();    
    int idx =  lat->dico_idx[ipl];
    if (idx==-1) continue;   // ...outside dico
    if (idx==-2) break;      // ...beyond dico

    if (H0(0)>xpl) {
      printf("\"%s\" %f %f %f\n",d.Name.c_str(),H0(0),xpl,xprv); return -1;
    }

    //       ********** EXTRAPOLATE... **********
    if (!error) {                   // ***** ... USING RKutta *****
      H1 = H0;
      while (fabs(H1(0)-xpl)>TConstants_RKuttaMinStep){ 
	double step = xpl-H0(0)<TOpt::dCut[3] ? xpl-H0(0) : TOpt::dCut[3];
	H1(0) = H0(0)+step;
	if (!H0.Extrapolate(H1)) {
	  error = 1; break;
	}
	H0 = H1;
      }
      if (!error) Hlast = H1; // IF OK: SAVE COORD'S FOR STRAIGHT LINE EXTRAP
    }
    if (error) {         // ***** ... OR STRAIGHT LINE IF RKutta FAILED... *****
      H1(0) = xpl;
      H1(1) = Hlast(1)+(xpl-Hlast(0))*Hlast(3);
      H1(2) = Hlast(2)+(xpl-Hlast(0))*Hlast(4);
    }

    coords[idx] = d.Ca*H1(1)+d.Sa*H1(2);

    H0 = H1;
  } // end of loop over planes

  return 0;
}
