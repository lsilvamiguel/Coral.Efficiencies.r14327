// $Id: PIDdoBeamID.cc,v 1.13 2008/01/08 18:02:27 ybedfer Exp $

/*!
  Flag BEAM tracks.
  Conditioned by option "PID_doBeamID selAllBeams":
  - <=0: Require BMS reco (which is performed by "../beam/CsBeamRecons" and
        yields "CsBeam" objects w/ a "BMSzone" zone.
  - <0 : Reject bad chi2 of BMS<->scifi/Si transport
  - >0 : Accept all beam tracks. The criterion is that track must have its
        first helix upstram of target. The "CsTrack" mothers of "CsBeam's" are
        disregarded (The latter are created by the "../beam" package, w/o
        the scifi/Si tracks there built upon (i.e. their mother) being erased).
         Note: In that case, care has to be taken (in the pieces of software
        producing the beam tracks) that tracks w/ and w/o are distinguishable.
*/

#include"CsEvent.h"
#include"CsOpt.h"
#include"CsErrLog.h"
#include"CsMCUtils.h"
#include"CsGeom.h"
#include"CsBeam.h"
#include<string>
#include<stdio.h>

using namespace std;

void PID_doBeamID( vector<CsParticle*>& parts ) {
  
  static bool first(true); static int selAllBeams(0);  // ***** INITIALISATION
  if (first) {
    CsOpt::Instance()->getOpt("PID_doBeamID","selAllBeams",selAllBeams);
    first = false;
  }


  vector<CsParticle*> badBeams; int nBeams;
  vector<int> ids;
  vector<CsParticle*>::iterator ip;
  for (ip = parts.begin(), nBeams = 0; ip!=parts.end(); ip++ ) {

    // ***** LOOP on ALL PARTICLES: SELECT those w/ BMS TRACKS...
    // ...and store bad chi2, BMS tracks and total beam multiplicity (if needed)

    const CsTrack *cst = (*ip)->getTrack(); if (cst==0) continue;
    const list<CsZone*> &zones = cst->getZones();
    list<CsZone*>::const_iterator iz;
    for (iz = zones.begin(); iz!=zones.end(); iz++) {
      if ((*iz)->getName()=="BMSzone") {

	//                                           *****  STORE BAD CHI2 BEAMS
	const CsBeam *beam; if ((beam = dynamic_cast<const CsBeam*>(cst)) &&
				beam->getChi2CutFlag()) badBeams.push_back(*ip);
	nBeams++;
	ids.push_back(cst->getId());  // ***** STORE CsTrack::id of BEAMS w/ BMS

	(*ip)->setType(CsParticle::SPECIAL);  // ***** FLAG BMS-PARTICLE AS BEAM
      }
    }
  }
  if (nBeams>1) {
    // ***** IF MULTIPLICITY >1, RE-EVALUATE BAD BEAMS' CASE... 
    for (ip = badBeams.begin(); ip!=badBeams.end(); ip++) {
      if (selAllBeams<0) // ...UPON OPTION: EXCLUDE THEM AUTOMATICCALY
	(*ip)->setType(CsParticle::ORDINARY);
      else {             // ...ELSE THEY MIGHT SHARE SAME ID w/ OTHER BEAM
	// (Note that this should never happen (by construction of the beam reco
	// software as of 2011/03). Yet, it does! Cf. event #20983512 in
	// 04W28/cdr09002-37277. And the consequence is that PHAST's
	// "PaEventImport", which is assuming CsTrack::id's to be unique,
	// instantiates only one PaParticle (the 1st encountered one), leaving
	// any primary vertex built on the second one w/o incident PaParticle.
	// Can pose many problems: to quote one, the "BestCoralPrimaryVertex",
	// if affected, ends up being a 2ndary vertex in the output mDST.
	//  I've never observed that two good beams share the same CsTrack::id.
	// If this latter case ever happens, the present patch does not fix it.)
	int id = (*ip)->getTrack()->getId(), i, match;
	for (i=match = 0; i<(int)ids.size(); i++) if (ids[i]==id) match++;
	if (match>1/* =1: no sharing */) (*ip)->setType(CsParticle::ORDINARY);
      }
    }
  }

  if (selAllBeams>0) {              // ***** BMS OPTIONAL *****

    double zTarget = CsGeom::Instance()->getTargetCenter( );
    vector<CsParticle*>::iterator ip;
    for (ip = parts.begin(); ip!=parts.end(); ip++) {
      const CsTrack *cst = (*ip)->getTrack(); if (cst==0) continue;
      const vector<CsHelix> &tpar = cst->getHelices();
      if (tpar.size() && tpar[0].getZ()<zTarget) {

	//                    ***** SELECT TRACKS STARTING UPSTREAM of TARGET...
	const list<CsZone*> &zones = cst->getZones();
	list<CsZone*>::const_iterator iz; int bms;
	for (iz = zones.begin(), bms = 0; iz!=zones.end(); iz++)
	  if((*iz)->getName()=="BMSzone") { bms = 1; break; }
	if (bms) continue;                               // ***** ...AND NOT BMS
	int id = cst->getId();
	for (int i = 0; i<(int)ids.size(); i++)
	  if (ids[i]==id) { bms = 2; break; }     // ***** ...AND NOT BMS MOTHER
	if (bms) continue;

	(*ip)->setType( CsParticle::SPECIAL );
      } 
    } 

  }
  
  
  return;
}








