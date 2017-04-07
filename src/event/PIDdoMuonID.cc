// $Id: PIDdoMuonID.cc,v 1.50 2010/12/06 02:54:49 ybedfer Exp $

/*!
  \file  PIDdoMuonID.cc
  \brief Scattered mu ID (so-called mu'-ID)

  - Enabled upon option "CsKalmanFitting Specials".
  - Purpose:
  - Required:
   i) Obvious things: exclude beam, fringe field, no P, low (< 1 GeV) P.
  ii) Same polarity as beam, which latter is determined on view of SM2 scaling
     factor.
 iii) XX0 > cut.
  iv) W/in acceptance of a trigger hodo system compatible w/ event's trigger
     pattern.
   v) Compatibility of pair of hodo hits, if exist, w/ trigger correlation matrix.
  - Note: there is a basic flaw in the mu'-ID (herein implemented): in its very
   interface! It returns a yes or no answer, based upon, among others, the
   trigger correlation matrix. This, while scanning all the bits of the trigger
   pattern corresponding to hodo-based triggers, if any. Now, let's consider
   a pure CT interaction w/ an accidentally coincident hodo trigger: it will not
   pass the mu'-ID, whereas it's a genuine one. A better interface would be one
   returning a pattern of bits, or one returning at least 2 bits: one inclusive,
   based upon hodo bits, if any, and one pure CT[|highQ2T], not requiring
   consistency w/ the trigger matrix.
*/

#include <stdio.h>
#include <string>
#include <list>
#include <vector>
#include <sstream>
#include <TH2.h>
#include "CsInit.h"
#include "CsErrLog.h"
#include "CsOpt.h"
#include "CsHistograms.h"
#include "CDB.h"
#include "CsRegistry.h"
#include "CsGeom.h"
#include "CsField.h"
#include "CsEvent.h"
#include "CsMCUtils.h"

using namespace std;
using namespace CLHEP;

/*! \brief Statistics output for 
  mu' identification algorithm.
*/
class CsPIDdoMuonIDstat: public CsEndOfJob {

public:
  CsPIDdoMuonIDstat();          //!< constructor
  ~CsPIDdoMuonIDstat();         //!< destructor
  static CsPIDdoMuonIDstat* Instance();
  void Stat( int i, unsigned int TM );
  void StatMW1( int i, unsigned int TM );
  
  //! "End of job" method
  bool end();
    
private:
  static CsPIDdoMuonIDstat* instance_;    //!< The singleton pointer 
  
  int NR, NC;
  float stat[20][11];  // number vs TM
  float statLastEv[20][11];  // number vs TM
  float all[20];
  float allLastEv[20];
};



/*! \brief Container for trigger matrices.
  Used for mu' identification.
*/
class TrigMatrix {

public:
  TrigMatrix( const string &name );
  const string &GetName() const { return name_; }
  void readCalibration(int timePoint);
  bool Accept(int i,int j);

private:
  bool matrix_[32][32];

  CDB *cdb_;
  string  name_;

};



// Interfaces
bool CheckTrigMatrix(const vector<CsCluster*> &clrsF, char cT);
vector<int> IsInside( const vector<CsCluster*> &clrsF, const string &name );


void PID_doMuonID( vector<CsParticle*>& parts ) {

  CsOpt *opt = CsOpt::Instance();
  static int beamCharge(0); static double zT, zSM2;
  static int doMuonID = -1;
  if (!doMuonID) // Early exit if mu'-ID not requested (cf. initialiation infra)
    return;
  else if (doMuonID<0) {// ********** INITIALIZATIONS (cf. also infra)**********
    if (opt->getOpt("CsKalmanFitting","Specials",doMuonID)) {
      if (doMuonID<0 || 2<=doMuonID) CsErrLog::msg(elFatal,__FILE__,__LINE__,
        "Ambiguous argument to \"CsKalmanFitting Specials\" (=%d): "
	"Must be 0 or 1",doMuonID);
    }
    else doMuonID = 0;
    if (!doMuonID) return;

    // ***** BEAM CHARGE: derive it from polarity of SM2 (=last but one magnet)
    CsGeom *g = CsGeom::Instance(); CsField *f = g->getCsField();
    CsMagInfo *mag = f->getMagInfo(); int nMags = f->getNumOfMags();
    if (nMags!=3) CsErrLog::msg(elFatal,__FILE__,__LINE__,
				"# of magnets = %d, while expected = 3",nMags);
    double sm2ScaleF = mag[2].fsc;
    if (!sm2ScaleF) // SM2=0: no telling >0 from <0 muons.
      CsErrLog::mes(elFatal,"SM2 Scaling factor = 0");
    beamCharge = sm2ScaleF>0 ? +1 : -1;
    zSM2 = mag[2].zcm;

    zT = g->getTargetCenter(); // Centre of target
  }

  static bool print(false);
  
  CsPIDdoMuonIDstat* stat = CsPIDdoMuonIDstat::Instance();
  map<CsParticle*,int> statZlast;
  
  //***** CUTS & OPTIONS *****
  bool muMC = false;  // true - MC ID, false - realistic ID.

  list<CsParticle*> partSp; // for MU' candidates
  vector<CsParticle*>::iterator ip;
  partSp.clear();
  
  unsigned int trigger_mask( CsEvent::Instance()->getTriggerMask() );
  stat->Stat(0,trigger_mask);

  for( ip=parts.begin(); ip!=parts.end(); ip++ ) {

    const CsTrack* trk = (*ip)->getTrack();
    if( trk == NULL ) continue;
    vector<CsHelix> tpar = trk->getHelices();

    if( tpar.size() < 2 ) continue;

    if( !muMC ) {    // ********** REALISTIC IDENTIFICATION **********

      //         ***** CUT on NUMBER of RAD LENGTHS => mu ID ***** 
      //  A cut of 30 would reject tracks that do pull the hodo triggers (L/MT),
      // and yet don't amass much XX0, for reasons that are still obscure. Part
      // of the explanation lies in the fact that the thickness of the concrete
      // wall, MUF2, is only 25 X0 (dixit Christian H.), and that some of the
      // muons traverse only that wall, later falling into the central hole of
      // the 2nd, MUF2P, wall. It can aslo be that MUF2 doesn't actually cover
      // completely the acceptance of the hodos (leaving them partially naked)
      // or that it's its descritpion in COMGeant setups, and hence in the MUF2
      // material maps, that's wrong.
      //  A relaxed cut of 20 does not seem to let fake muons through: the
      // distribution of the muons as a f(distance to the axis) remains smooth,
      // meaning w/o any sharp step that would correspond to the fake muons. On
      // the contrary, w/ a cut of 30, this same distribution exhibits a sharp
      // drop (which turns out to be in fact a dip: the tracks traveling on the
      // very edge of the acceptance of the L/M hodos being recuperated by the
      // inner hodo HI5 and the muWall (MUF3) that protects it. The former
      // situation, smoother, is thought to be more reproducible (e.g. by MC)
      // than the latter and is thus retained.
      //  Added lately: I finally opt for 15, despite the fact that I haven't
      // yet checked whether this could let through many more fakes (most
      // probably not). The rationale being: i) enhanced efficiency, ii) the
      // smoother the distribution (the more shallow the dip), the easier it's
      // reproducible in MC, iii) in any case, the cut can be me stricter when
      // the mDSTs produced by coral are analyzed, while the impact of the
      // mu'-ID (which the mu-ID conditions) on the rest of coral is mild: only
      // that in the Kalman filter step of the vertexing, the ID'd mu' will be
      // given precedence over plain tracks.
      if (trk->getXX0()<15) continue;
      (*ip)->setName("mu");
      
      //                  ***** REJECTIONS... *****
      if (tpar.front().getZ()<zT) continue;        // ...Beam track      
      if (tpar[0].getCop()== 0) continue;          // ...Track without momentum
      if (fabs(tpar[0].getCop())>1) continue;      // ...Momentum is very small
      static int DYspecialMuon=-1;
      if (DYspecialMuon==-1) {
	string value;
	CsOpt::Instance()->getOpt("CsKalmanFitting", "DY_MUON_ID_POLARITY", value );
	if (value=="IGNORE") DYspecialMuon=1;
	else                 DYspecialMuon=0;
      }
      if (!DYspecialMuon && 
	  tpar[0].getCop()*beamCharge<0) continue; // ...Wrong polarity
      if (trk->getZones().size()==1) continue;     // ...Fringe-field (w/ mom but only 1 zone)
      if (tpar.front().getZ()>zSM2) continue;      // ...Non bridged, 0xc, track

      stat->Stat(12,trigger_mask);
      
      // ***** NOT CHECKING CONSISTENCY WITH TARGET *****
      // (Checking it:
      //   i) Would not add to the reliability of the mu'-ID, for a mu' can very
      //     well arise from outside the target, since the target pointing that
      //     is built into the trigger logic is only approximate.
      //  ii) Consistency w/ target is anyway expected to be checked at mDST
      //     anlysis level.
      // iii) A fake mu' that would not cross the target would only impact on
      //     the building of a vertex also outside target: not that dramatic.
      // On the contrary allowing mu' outside target will allow to have a better
      // accounting of trigger performances. E.g. one will be able to tell what
      // fraction of triggers, useless ones, are due to bad target pointing.) 

      stat->Stat(13,trigger_mask);
      
      //              *************** CONSISTENCY WITH HODOSCOPE ***************
      if (true) {
	static bool MF(true);    // flag for trigger matrix check
	static int hod_trk(1);
	static double ClrTime(0);
	const map<unsigned int,char>& MASK = CsInit::Instance()->getTCSMasks();
	static list<list<CsDetector*> > hodoSystems; // Typically I, M, L, O
	static unsigned int caloMask = 0, hodoMask = 0;
	const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
	static bool first(true); if (first) {
	  //                     ********** INITIALIZATIONS **********
	
	  //                               ***** OPTIONS *****
	  list<string> names;
	  list<CsDetector*> Dets, dets( CsGeom::Instance()->getDetectors() );
	  list<CsDetector*>::iterator id;
	  list<string>::iterator is;
          
	  int itmp(0);
	  opt->getOpt("Trigger","Print",itmp); if( itmp != 0 ) print = true;
	  opt->getOpt("Trigger","Matix",itmp); if( itmp == 0 ) MF    = false;
	  opt->getOpt("Trigger","tracking",hod_trk);
	  opt->getOpt("Trigger","ClrTime",ClrTime);


	  map<unsigned int,char>::const_iterator in;
	  for (in = MASK.begin(); in!=MASK.end(); in++) {
	    unsigned int mask = (*in).first; char tag = (*in).second;
	    if (tag=='C' || tag=='Q') caloMask |= mask;
	    if (tag=='I' || tag=='M' || tag=='L' || tag=='O' || tag=='G')
	      hodoMask |= mask; 
	  }

	  while (opt->getOptRec("Trigger","det",names)) {
	    for( id=dets.begin(); id!=dets.end(); id++ ) {
	      string name = (*id)->GetTBName();
	      for (is = names.begin(); is!=names.end(); is++) {
	        if (name.find(*is)==0) Dets.push_back( *id );
	      }
	    }
	    if (Dets.size()!=0) {
	      hodoSystems.push_back(Dets);
	      Dets.clear();
	    }
	    
	  }
	  
	  //                            ***** CALIBRATIONS *****
	  CsInit *csInit = CsInit::Instance();
	  if( csInit->useCalibration() && csInit->useCDB() && MF ) MF = true;
	  else                                                     MF = false;
	  
	  first = false;
	} // End of initialisation

	if (print) printf("muonID: trigger_mask=0x%x,  Mom=%.1f\n",
			  trigger_mask,1/tpar[0].getCop());
	
	bool exists(false), goodTrack(false);
	double Utr,Vtr, Uhod,Zhod, VMAX, VMIN, SIGMAU, SIGMAV;
	double *covT;
	HepMatrix covH;
	string name;
	vector<CsCluster*> clrsF;
	char cT='0';
	
	list< list<CsDetector*> >::iterator iS;
	list<CsCluster*>::const_iterator ic;
	list<CsDetector*>::iterator id;
	
	// Check the given trigger is in list of known trigger
	map< unsigned int, char >::const_iterator in;
	for (in = MASK.begin(); in!=MASK.end(); in++) {
	  if (trigger_mask&(*in).first) {
	    exists = true;
	    if (print) cout<<"muonID: current trigger is "<<(*in).second<<endl;
	    break;
	  }
	}
	
	if (!exists &&            // ********** TRIGGER is UNKNOWN... **********
	    tpar.back().getZ()<33000) // => REQUIRE muWall2 RECO 
	  continue;        // (i.e. no attempt @ extrapolating to hodos)

	if ((trigger_mask&caloMask) && !(trigger_mask&hodoMask)) {
	  //                             ********** PURE CALO TRIGGER **********
	  // - Pure calo includes both CT and highQ2T, the latter despite it's
	  //  a mixture of calo and hodo (HQ) and because HQ is very special
	  //  and it would have been difficult to treat it like HI/M/L/O.
	  // - "Pure" means no hodo contribution. Other contributions are OK,
	  //  e.g.: beam trigger.
	  // - No other requirement than XX0, by default. Yet...
	  // - ...The threshold earlier passed may be considered too loose: it's
	  //  been set so as to let through all the hodo triggering muons
	  //  traveling on the rim of the muFilter central hole, which happen to
	  //  collect few X0's, be it that the hodos (it's HL and HM) are
	  //  partially naked in the real life or that it's their COMGeant
	  //  description (ported to material map) that makes them look so from
	  //  coral's point of view.
	  //   One may also consider that the X-check of their consistency w/
	  //  the trigger logic grants the hodo triggering muons some extra
	  //  reliability. And then one may want to compensate the lack of such
	  //  extra reliability in pure-CT by being more demending on XX0.
	  //   If indeed, one can enable the "PID_muID_STRICT_XX0" infra.
#ifdef PID_muID_STRICT_XX0
	  if ((33000<tpar.back().getZ() && tpar.back().getZ()<50000 && 
	       // Require XX0 > THICKNESS of MUF2. I (Y.B.) didn't know a priori
	       // how thick is this. I determined it empirically from RD and
	       // MC(D*) data...
	       trk->getXX0()>65) ||
	      // ...Or special case of tracks through MUF3, where one may want
	      // to be even stricter, rejecting muons traveling on the edge.
	      trk->getXX0()>74)
	    partSp.push_back(*ip);
#else
	  partSp.push_back(*ip);
#endif
	  continue;
	}

	//      ********** CHECK COMPATIBILITY w/ TRIGGER HODOs **********

	// ***** LOOP over HODOSCOPE SYSTEMS *****
	for (iS = hodoSystems.begin(); iS!=hodoSystems.end() && exists; iS++ ) {

	  clrsF.clear();
	  name = (*iS).front()->GetTBName();

	  bool goodSystem = false;

	  // Check whether TriggerMask corresponds to current system
	  map< unsigned int, char >::const_iterator in;
	  for (in = MASK.begin(); in!=MASK.end(); in++)
	    if (name[1]==(*in).second && ((*in).first& trigger_mask)) {
	      goodSystem = true; cT = name[1];
	      if (print) cout<<"muonID: the first detector of group is "<<name<<endl; 
	      break;
	    }
	  if (!goodSystem) // Trigger mask doesn't correspond to current system
	    continue;

	  //                          ***** LOOP over HODOSCOPE PLANES in SYSTEM
	  for (id = (*iS).begin(); id!=(*iS).end(); id++) {
	    CsDetector *hodo = *id; name = hodo->GetTBName();

	    bool goodPlane = false;

            // Check whether hodoscope cluster is included in track
	    if( hod_trk != 0 && tpar.back().getZ() - hodo->getZcm() > -0.1 ) {
	      // Track can contain cluster of fired hodo
	      const list<CsCluster*> &cts = trk->getClusters();
	      for (ic = cts.begin(); ic!=cts.end(); ic++) {
		CsCluster *ct = *ic;
	        const list<CsDetector*> &dets = ct->getDetsList();
	    	const string &nameCl = dets.front()->GetTBName();
		if (name!=nameCl) continue;
		double time(0);
		if (ClrTime!=0) {// ***** REQUEST CLUSTER TIME CONSISTENCY *****
		  if (!ct->getTime(time)) {        // ***** REQUIRE CLUSTER TIME
		    if (print)
		      printf("muonID: %.2fGeV, TRACK hit %s %6.1fmm, no time!\n",
			     1/tpar[0].getCop(),name.c_str(),ct->getU());
		    continue;
		  }
		  if (fabs(time)>ClrTime) {                 // ***** OUT OF TIME
		    if (print)
		      printf("muonID: %.2fGeV, TRACK hit %s %6.1fmm, %4.1fns OUT\n",
			     1/tpar[0].getCop(),name.c_str(),ct->getU(),time);
		    continue;
		  }
		}
		goodPlane = true;
		if (print) 
		  printf("muonID: %.2fGeV, TRACK hit %s %6.1fmm %4.1fns OK\n",
			 1/tpar[0].getCop(),name.c_str(),ct->getU(),time);
		clrsF.push_back(ct);
	      }
	      goodSystem &= goodPlane;
	      if (goodPlane) {  // Current hodo system has valid cluster
		if( statZlast[ *ip ] != 2 ) {   // for statistics
		  if( id == (*iS).begin() ) statZlast[ *ip ] = 1;
		  else                      statZlast[ *ip ] = 2;
		}
	      }
	    }
	  }
	  if (goodSystem && MF) goodSystem &= CheckTrigMatrix(clrsF,cT);
	  if (goodSystem) { goodTrack = true; break; }

	  //   ***** NO FULLY CONSISTENT SET OF ASSOCIATED HODO HITS *****
	  //   *****   RESCUE: CHECK THAT TRACK CROSSES FIRED SLAT   *****
	  for (id = (*iS).begin(), goodSystem = true; id!=(*iS).end(); id++) {
	    //                        ***** LOOP over HODOSCOPE PLANES in SYSTEM
	    CsDetector *hodo = *id; name = hodo->GetTBName();

	    bool goodPlane = false;

	    const list<CsCluster*> &clrs = hodo->getMyClusters();
	    CsHelix Hhod; Zhod = hodo->getZcm();
            if (!tpar.back().Extrapolate(Zhod,Hhod)) // ***** EXTRAP FAILS *****
	      break;

	    covT = Hhod.getCov();
	    if( name[4] == 'X' ) {
	      Utr  = Hhod.getX(); Vtr  = Hhod.getY();
	      VMAX = hodo->getYcm() + hodo->getYsiz()/2;
	      VMIN = hodo->getYcm() - hodo->getYsiz()/2;
	    } else 
	    if( name[4] == 'Y' ) {
	      Utr  = Hhod.getY(); Vtr  = Hhod.getX();
	      VMAX = hodo->getXcm() + hodo->getXsiz()/2;
	      VMIN = hodo->getXcm() - hodo->getXsiz()/2;
	    } else {
	      cout<<"muonID:ERROR: unclear coordinate to check "<<name[4]<<endl;
	      continue;
	    }

	    // ********** LOOP ON HODO HITS TO FIND COMPATIBLE ONES **********
	    for (ic = clrs.begin(); ic!=clrs.end(); ic++) {
	      CsCluster *ct = *ic;
	      if (ClrTime!=0) {  // ***** REQUEST CLUSTER TIME CONSISTENCY *****
	        double time(0);
		if (!ct->getTime(time)) {          // ***** REQUIRE CLUSTER TIME
		  if (print)
		    printf("muonID: %.2fGeV, HODO  hit %s %6.fmm, no time!\n",
			 1/tpar[0].getCop(),name.c_str(),ct->getU());
		  continue;
		}
		if (fabs(time)>ClrTime) {                   // ***** OUT OF TIME
		  if (print)
		    printf("muonID: %.2fGeV, HODO  hit %s %6.1fmm, %.1fns OUT\n",
			   1/tpar[0].getCop(),name.c_str(),ct->getU(),time);
		  continue;
		}
		if (print) 
		  printf("muonID: %.2fGeV, HODO  hit %s %6.1fmm %4.1fns OK\n",
			 1/tpar[0].getCop(),name.c_str(),ct->getU(),time);
	      }
	      
	      Uhod = ct->getU(); covH = ct->getCov();
	      if( name[4] == 'X' ) {
		//SIGMAU = 0.5 * sqrt( 12*covH(1,1) ) + 3*sqrt(covT[0]);
		SIGMAU = sqrt( 12*covH(1,1) );
		SIGMAV = sqrt(covT[2]);
	      } else {
		//SIGMAU = 0.5 * sqrt( 12*covH(1,1) ) + 3*sqrt(covT[2]);
		SIGMAU = sqrt( 12*covH(1,1) );
		SIGMAV = sqrt(covT[0]);
	      }

	      if( fabs(Utr-Uhod) < SIGMAU &&
	          Vtr > VMIN - 3 * SIGMAV &&
		  Vtr < VMAX + 3 * SIGMAV ) {
	        goodPlane = true;
 		if (print) {
		  printf("muonID: %.2fGeV, HODO  hit %s %6.1fmm = "
			 "track%+5.1fmm, dhit = %.1fmm, dtrack = %4.1fmm\n",
			 1/tpar[0].getCop(),name.c_str(),Uhod,Uhod-Utr,
			 SIGMAU,SIGMAV);
		}
		clrsF.push_back(ct);
	      }
	    }
	    goodSystem &= goodPlane;
	    if (!goodSystem) break; // Current hodo system !OK
	  }
	  // TRIGGER MATRIX STUFF
	  if (goodSystem && MF) {
	    if (!CheckTrigMatrix(clrsF,cT)) {
	      goodSystem = false;
	      if( print ) cout<<"muonID: matrix false"<<endl;
	    }
	  }
	  if (goodSystem) {
	    if( print ) cout<<"muonID:"<<cT<<" hodo system OK"<<endl;
	    goodTrack = true; break;  // Current hodo system passed all criteria
	  }
	}

        if (!goodTrack && exists) {
          if (hodoSystems.size()!= 0) {
	    if( print ) cout<<"muonID: next candidate mu' track"<<endl;
	    continue;   // next mu' track
	  }
	  // if no trigger related strings in option file check  
          // for last helix to be downstream of concrete absorber  
	  else if( tpar.back().getZ() < 33000 ) continue;
	}
      } // end of hodoscopes stuff
      else {
        // for last helix to be downstream of concrete absorber 
        if( tpar.back().getZ() < 33000 ) continue;
      }
      
      if( print ) cout<<" ######## SCAT MU ######"<<endl;
      partSp.push_back( *ip );
      //break; // only one mu' for the moment
      
    } else {                               //***** MC IDENTIFICATION
      
      list<CsMCTrack*>trtrs = CsMCUtils::getAssociatedMCTracks(trk,85);
      if(trtrs.size()>0){
	int nhits;
	CsMCTrack* trkMC=CsMCUtils::getAssociatedMCTrack(trk,nhits);
	
	if(trkMC->getParticle()->getGeantNumber()==5) (*ip)->setName("mu");
	else continue;

	if(trkMC->getInVertex()->getGnum()==1 && trkMC->getGnum()==2 )
	  (*ip)->setType( CsParticle::SPECIAL );  //*** IT IS MU'
	
      }
      
    }

  }
  
  if( !muMC ) {                             //***** REALISTIC IDENTIFICATION
    if( partSp.size() > 1 ) {
      list<CsParticle*>::iterator ips;
      for( ips=partSp.begin(); ips!=partSp.end(); ips++ ) {
	(*ips)->setName("mu");
	(*ips)->setType( CsParticle::SPECIAL );
      }
      
    } else if( partSp.size() == 1 ) {
      partSp.front()->setName("mu");
      partSp.front()->setType( CsParticle::SPECIAL );    
    }
    
    ///////// statistics //////////
    if( partSp.size() > 0 ) {
                              stat->Stat(1,trigger_mask);
      if( partSp.size() > 1 ) stat->Stat(2,trigger_mask);
      
      int globZlast = 0;
      list<CsParticle*>::iterator ips;
      for( ips=partSp.begin(); ips!=partSp.end(); ips++ ) {
	if( statZlast[ *ips ] == 1 && globZlast < 1 ) { globZlast = 1; }
	if( statZlast[ *ips ] == 2 )                  { globZlast = 2; break;}
      }
      if(  globZlast == 0   ) stat->Stat(3,trigger_mask);
      if(  globZlast == 1   ) stat->Stat(4,trigger_mask);
      if(  globZlast == 2   ) stat->Stat(5,trigger_mask);
            
    }
    ///////// end of statistics //////////
    
  }
    
  return;
}

/////////////////////////////////////////////////////////////////////
///////////////////////////// ID for CT /////////////////////////////

void muIDinMW1( vector<CsParticle*>& parts )
{
  static int doMuonID = -1;
  if (!doMuonID) // Early exit if mu'-ID not requested (cf. initialiation infra)
    return;

  static bool print(false), hist(false); static int doID(0), beamCharge(0);
  static TH2F *h2[10];
  static bool first(true); if (first) { // ********** INITIALIZATIONS **********
    CsOpt *opt = CsOpt::Instance();
    if (opt->getOpt("CsKalmanFitting","Specials",doMuonID)) {
      if (doMuonID<0 || 2<=doMuonID) CsErrLog::msg(elFatal,__FILE__,__LINE__,
        "Ambiguous argument to \"CsKalmanFitting Specials\" (=%d): "
	"Must be 0 or 1",doMuonID);
    }
    else doMuonID = 0;
    if (!doMuonID) return;

    opt->getOpt("Trigger","Print",print);
    opt->setIntMode("hex"); opt->getOpt("Trigger","MW1_ID",doID);
    opt->setIntMode("dec"); opt->getOpt("Trigger","Hist" ,hist);
    if (hist) {
      CsHistograms::SetCurrentPath("/TRIGGER/");
      h2[0] = new TH2F("h2_CT","CT: MA hits", 9,0,9, 9,0,9 );
      h2[0]->SetXTitle("MA01"); h2[0]->SetYTitle("MA02");      
      CsHistograms::SetCurrentPath("/");
    }  

    // ***** BEAM CHARGE: derive it from polarity of SM2 (=last but one magnet)
    CsGeom *g = CsGeom::Instance(); CsField *f = g->getCsField();
    CsMagInfo *mag = f->getMagInfo(); int nMags = f->getNumOfMags();
    if (!nMags) CsErrLog::mes(elFatal,"No magnet in \"CsField\"");
    double sm2ScaleF = mag[nMags-1].fsc;
    if (!sm2ScaleF) // SM2=0: no telling >0 from <0 muons.
      CsErrLog::mes(elFatal,"SM2 Scaling factor = 0");
    beamCharge = sm2ScaleF>0 ? +1 : -1;

    first = false;
  }

  //   ***** PER EVENT INITIALIZATIONS *****
  int Nmu = 0;
  CsPIDdoMuonIDstat *stat = CsPIDdoMuonIDstat::Instance();
  unsigned int trigger_mask( CsEvent::Instance()->getTriggerMask() );

  if (trigger_mask&doID) {
    
    vector<CsParticle*>::iterator ip;
    for( ip=parts.begin(); ip!=parts.end(); ip++ ) {
      
      const CsTrack* trk = (*ip)->getTrack();
      if( trk == 0 ) continue;
      
      //                  ***** REJECTIONS... *****
      vector<CsHelix> tpar = trk->getHelices();
      if (tpar[0].getCop()== 0) continue;          // ...Track without momentum
      if (fabs(tpar[0].getCop())>1) continue;      // ...Momentum is very small
      static int DYspecialMuon2=-1;
      if (DYspecialMuon2==-1) {
	string value;
	CsOpt::Instance()->getOpt("CsKalmanFitting", "DY_MUON_ID_POLARITY", value );
	if (value=="IGNORE") DYspecialMuon2=1;
	else                 DYspecialMuon2=0;
      }
      if (!DYspecialMuon2 && 
	  tpar[0].getCop()*beamCharge<0) continue; // ...Wrong polarity
      if (tpar.back().getZ()<15000) continue;      // ...Upstream of MA01
      if (trk->getXX0()<30)                        // ...XX0 too small
	// Hits in MAs does not ensure muID per se. We want the track to have
	// actually traversed the mu absorber (and not passed through its hole).
	continue;

      // Not checking consistency with target (cf. comment in "PID_doMuonID").

      //             ***** CHECK # OF HITS in MAs *****
      int N_MA01 = 0, N_MA02 = 0;
      const list<CsCluster*> &clrs = trk->getClusters();
      list<CsCluster*>::const_iterator ic;
      for (ic = clrs.begin(); ic!=clrs.end(); ic++) {
	const list<CsDetector*> &dets = (*ic)->getDetsList();
	const string &nameCl = dets.front()->GetTBName();
	if (nameCl.find("MA01")==0) N_MA01++;
	if (nameCl.find("MA02")==0) N_MA02++;
      }
      if (print && (N_MA01!=0 || N_MA02!=0))
	cout<<"muIDinMW1: #hits(MA01)="<<N_MA01<<", #hits(MA02)="<<N_MA02<<endl;
      if (N_MA01<5 || N_MA02<5) continue;

      //     ***** UPDATE: CsParticle, histo, # of mu's
      (*ip)->setName("mu"); (*ip)->setType(CsParticle::SPECIAL );
      if (hist) h2[0]->Fill(N_MA01,N_MA02);
      Nmu++;
    }
  }
  
  if (Nmu) stat->StatMW1(1,trigger_mask);   // For statistics
}





/////////////////////////////////////////////////////////////
////////////////////// other functions //////////////////////

vector<int> IsInside(const vector<CsCluster*> &clrsF, const string &name) {
  // Returns the indices of the clusters belonging to the argument TBname in the
  // argument vector of clusters.
  vector<int> ID;
  vector<CsCluster*>::const_iterator ic; int i;
  for (ic = clrsF.begin(), i = 0; ic!=clrsF.end(); ic++, i++) {
    const list<CsDetector*> &dets = (*ic)->getDetsList();
    const string &nameCl = dets.front()->GetTBName();
    if (nameCl.find(name)==0) ID.push_back(i);
  }
  return ID;
}




////////////////////////////////////////////////////////////////////////////
///////////////////////////////  TrigMatrix  ///////////////////////////////

TrigMatrix::TrigMatrix( const string &name ) : name_(name) {
  for(int i=0;i<32;i++) for(int j=0;j<32;j++) matrix_[i][j]=true;
  cdb_ = CsInit::Instance()->getDB(); 
}

void TrigMatrix::readCalibration(int timePoint) {
  CDB::Time tp(timePoint,0);
  string data;
  string tmp;
  float calib;
  cdb_->read(GetName(),data,tp);
  if(data.size()!=0) {
    //cout << data << endl<< endl;
    for(unsigned int i=0;i<data.size();i++) if( data[i] != '1' && data[i] != '0' ) {data.replace(i,1,"");i--;}
    if( data.size() != 32*32 ) cout<<"TrigMatrix::readCalibration==> size of "<<name_<<" matrix is not 32*32"<<endl;
    for(unsigned int i=0;i<32;i++) for(unsigned int j=0;j<32;j++) 
      if( 32*i+j < data.size() && data[32*i+j] != '0'  ) matrix_[i][j]=true; else matrix_[i][j]=false;
  }else{
    tm *t = localtime((time_t*)&tp.first);
    cout << GetName() << ", no calibration for local time "
	 << t <<" in CDB"<< endl;
  }
}

inline bool TrigMatrix::Accept(int i,int j) {
  if( i<0 || i>31 ) return false;
  if( j<0 || j>31 ) return false;
  return matrix_[31-i][j];
}

bool CheckTrigMatrix(const vector<CsCluster*> &clrsF, char cT)
{
  // TRIGGER CORRELATION MATRICES
  // - They are retrieved from calibration files.
  // - Ideally those should only be read when the corresponding trigger is
  //  active.
  // - This is indeed the case for the 'G' trigger system, which was not
  //  active in a large fraction of COMPASS data, viz. [2002,2007].
  // - Yet in order to minimize the # of calls to the calibration DB, all of
  //  'IMLO' matrices are read in one go.

  bool ACCEPT(false);

  static TrigMatrix *MIu,*MId, *MMXu,*MMXd,*MMYu,*MMYd, *ML, *MO1,*MO2, *MG;
  static int hist = 0; if (!hist) {
    CsOpt::Instance()->getOpt("Trigger","Hist",hist);
    if (hist<=0) hist = -1;  // Meaning "hist" defined w/ no histogramming.
  }
  static TH2F *hIu[2],*hId[2], *hMXu[2],*hMXd[2],*hMYu[2],*hMYd[2], *hL[2], *hO1[2],*hO2[2], *hG1[2],*hG2[2];

  if (cT=='G') {   //   ********** GEANT HODO SYSTEM **********
    static bool first(true);if (first) {                // ***** INITIALIZATIONS
      first = false;
      MG = new TrigMatrix("HG__Y_M_");
      //                                                      ***** CALIBRATIONS
      CsInit *csInit = CsInit::Instance(); int time = csInit->getCDBUseTime();
      if (time==0) time = csInit->getStartOfRun();
      if (time==0) time = CsEvent::Instance()->getEventTime().secFrEpoch(); //- 2*60*60;//Should be improve! 
      cout<<"muonID: Reading HG trigger matrix calibration\n";
      if (!csInit->getDB()->ConnectDB())
	throw CS::Exception("muonID: can't connect to CDB database");
      MG->readCalibration(time);
      csInit->getDB()->DisconnectDB();
      cout << "muonID: End of reading HG trigger matrix calibration\n";

      if (hist>0) {
	CsHistograms::SetCurrentPath("/TRIGGER/");
	hG1 [0] = new TH2F("hG10" ,"GT: MATRIX" , 32,0,32, 32,0,32 ); hG1 [0]->SetXTitle("HG2"); hG1 [0]->SetYTitle("HG1");
	hG1 [1] = new TH2F("hG11" ,"GT: MATRIX" , 32,0,32, 32,0,32 ); hG1 [1]->SetXTitle("HG2"); hG1 [1]->SetYTitle("HG1");
	hG2 [0] = new TH2F("hG20" ,"GT: MATRIX" , 32,0,32, 32,0,32 ); hG2 [0]->SetXTitle("HG2"); hG2 [0]->SetYTitle("HG1");
	hG2 [1] = new TH2F("hG21" ,"GT: MATRIX" , 32,0,32, 32,0,32 ); hG2 [1]->SetXTitle("HG2"); hG2 [1]->SetYTitle("HG1");
	CsHistograms::SetCurrentPath("/");
      }
    }
    int x1(-1),x2(-1),y1(-1),y2(-1);  
    //                              ***** INDICES of CLUSTERS BELONGING down HGs
    vector<int> iX = IsInside(clrsF,"HG02Y1__");
    vector<int> iY = IsInside(clrsF,"HG01Y1__");
    for (int i = 0; i<(int)iX.size(); i++) {
	for (int j = 0; j<(int)iY.size(); j++) {
	  x1 = clrsF[ iY[j] ]->getDigitsList().front()->getAddress();
	  x2 = clrsF[ iX[i] ]->getDigitsList().front()->getAddress();
	  ACCEPT = MG->Accept(x1,x2);
	  if (hist>0) {
	    hG1[0]->Fill(x2,x1); if (ACCEPT) hG1[1]->Fill(x2,x1);
	  }
	  if (ACCEPT) break;
	}
	if (ACCEPT) break;
      }    
      iX = IsInside(clrsF,"HG02Y2__");
      iY = IsInside(clrsF,"HG01Y1__"); 
      for (int i = 0; i<(int)iX.size(); i++) {
	for (int j = 0; j<(int)iY.size(); j++) {
	  x1 = clrsF[ iY[j] ]->getDigitsList().front()->getAddress();
	  x2 = clrsF[ iX[i] ]->getDigitsList().front()->getAddress();
	  ACCEPT = MG->Accept(x1,x2);
	  if (hist>0) {
	    hG2[0]->Fill(x2,x1); if (ACCEPT) hG2[1]->Fill(x2,x1);
	  }
	  if (ACCEPT) break;
	}
	if (ACCEPT) break;
      }
  }
  else {   //   ******************** IMLO HODO SYSTEMS ********************
    static bool first = true; if (first) {              // ***** INITIALIZATIONS
      first = false;
      MIu  = new TrigMatrix("HI__X_Mu"); MId  = new TrigMatrix("HI__X_Md");
      ML   = new TrigMatrix("HL__X_M_");
      MMXu = new TrigMatrix("HM__X_Mu"); MMXd = new TrigMatrix("HM__Y_Mu");
      MMYu = new TrigMatrix("HM__Y_Mu"); MMYd = new TrigMatrix("HM__Y_Md");
      MO1  = new TrigMatrix("HO__X1M_"); MO2  = new TrigMatrix("HO__X2M_");
      //                                                      ***** CALIBRATIONS
      CsInit *csInit = CsInit::Instance(); int time = csInit->getCDBUseTime();
      if (time==0) time = csInit->getStartOfRun();
      if (time==0) time = CsEvent::Instance()->getEventTime().secFrEpoch(); //- 2*60*60;//Should be improve! 
      cout<<"muonID: Reading IMLO trigger matrix calibrations\n";
      if (!csInit->getDB()->ConnectDB())
	throw CS::Exception("muonID: can't connect to CDB database");
      MIu ->readCalibration(time); MId ->readCalibration(time);
      MMXu->readCalibration(time); MMXd->readCalibration(time);
      MMYu->readCalibration(time); MMYd->readCalibration(time);
      ML  ->readCalibration(time);
      MO1 ->readCalibration(time); MO2 ->readCalibration(time);
      csInit->getDB()->DisconnectDB();
      cout << "muonID: End of reading IMLO trigger matrix calibrations\n";

      if (hist>0) {                                          // ***** HISTOGRAMS
	CsHistograms::SetCurrentPath("/TRIGGER/");    
	hIu [0] = new TH2F("hIu0" ,"IT: MATRIX" , 32,0,32, 32,0,32 ); hIu [0]->SetXTitle("HI5"); hIu [0]->SetYTitle("HI4");
	hIu [1] = new TH2F("hIu1" ,"IT: MATRIX" , 32,0,32, 32,0,32 ); hIu [1]->SetXTitle("HI5"); hIu [1]->SetYTitle("HI4");
	hId [0] = new TH2F("hId0" ,"IT: MATRIX" , 32,0,32, 32,0,32 ); hId [0]->SetXTitle("HI5"); hId [0]->SetYTitle("HI4");
	hId [1] = new TH2F("hId1" ,"IT: MATRIX" , 32,0,32, 32,0,32 ); hId [1]->SetXTitle("HI5"); hId [1]->SetYTitle("HI4");

	hMXu[0] = new TH2F("hMXu0","MXT: MATRIX", 32,0,32, 32,0,32 ); hMXu[0]->SetXTitle("HM5X");hMXu[0]->SetYTitle("HM4X");
	hMXu[1] = new TH2F("hMXu1","MXT: MATRIX", 32,0,32, 32,0,32 ); hMXu[1]->SetXTitle("HM5X");hMXu[1]->SetYTitle("HM4X");
	hMXd[0] = new TH2F("hMXd0","MXT: MATRIX", 32,0,32, 32,0,32 ); hMXd[0]->SetXTitle("HM5X");hMXd[0]->SetYTitle("HM4X");
	hMXd[1] = new TH2F("hMXd1","MXT: MATRIX", 32,0,32, 32,0,32 ); hMXd[1]->SetXTitle("HM5X");hMXd[1]->SetYTitle("HM4X");
	hMYu[0] = new TH2F("hMYu0","MYT: MATRIX", 32,0,32, 32,0,32 ); hMYu[0]->SetXTitle("HM5Y");hMYu[0]->SetYTitle("HM4Y");
	hMYu[1] = new TH2F("hMYu1","MYT: MATRIX", 32,0,32, 32,0,32 ); hMYu[1]->SetXTitle("HM5Y");hMYu[1]->SetYTitle("HM4Y");
	hMYd[0] = new TH2F("hMYd0","MYT: MATRIX", 32,0,32, 32,0,32 ); hMYd[0]->SetXTitle("HM5Y");hMYd[0]->SetYTitle("HM4Y");
	hMYd[1] = new TH2F("hMYd1","MYT: MATRIX", 32,0,32, 32,0,32 ); hMYd[1]->SetXTitle("HM5Y");hMYd[1]->SetYTitle("HM4Y");

	hL [0] = new TH2F("hL0" ,"LT: MATRIX" , 32,0,32, 32,0,32 ); hL [0]->SetXTitle("HL5"); hL [0]->SetYTitle("HL4");
	hL [1] = new TH2F("hL1" ,"LT: MATRIX" , 32,0,32, 32,0,32 ); hL [1]->SetXTitle("HL5"); hL [1]->SetYTitle("HL4");
  
	hO1 [0] = new TH2F("hO10" ,"OT: MATRIX" , 32,0,32, 32,0,32 ); hO1 [0]->SetXTitle("HO4"); hO1 [0]->SetYTitle("HO3");
	hO1 [1] = new TH2F("hO11" ,"OT: MATRIX" , 32,0,32, 32,0,32 ); hO1 [1]->SetXTitle("HO4"); hO1 [1]->SetYTitle("HO3");
	hO2 [0] = new TH2F("hO20" ,"OT: MATRIX" , 32,0,32, 32,0,32 ); hO2 [0]->SetXTitle("HO4"); hO2 [0]->SetYTitle("HO3");
	hO2 [1] = new TH2F("hO21" ,"OT: MATRIX" , 32,0,32, 32,0,32 ); hO2 [1]->SetXTitle("HO4"); hO2 [1]->SetYTitle("HO3");
	CsHistograms::SetCurrentPath("/");
      }
    }
    int x1(-1),x2(-1),y1(-1),y2(-1);  

    if (cT=='I') {  //   ********** INNER HODO SYSTEM **********
      //                            ***** INDICES of CLUSTERS BELONGING down HIs
      vector<int> iX = IsInside(clrsF,"HI05X1_d");
      vector<int> iY = IsInside(clrsF,"HI04X1_d");
      for (int i = 0; i<(int)iX.size(); i++) {
	for (int j = 0; j<(int)iY.size(); j++) {  // ***** LOOP ON CLUSTER PAIRS
	  x1 = clrsF[iY[j]]->getDigitsList().front()->getAddress();
	  x2 = clrsF[iX[i]]->getDigitsList().front()->getAddress();
	  ACCEPT = MId->Accept(x1,x2);   // ***** CHECK WHETHER PAIR w/IN MATRIX
	  if (hist>0) {
	    // Note: Histogramming is bugged. Only 1st pair being histo'd.
	    hId[0]->Fill(x2,x1); if (ACCEPT) hId[1]->Fill(x2,x1);
	  }
	  if (ACCEPT) break;
	}
	if (ACCEPT) break;
      }
      //                                                    ***** SAME W/ up HIs
      // Note: Useless if already ACCEPTED, unless one wants histogramming.
      iX = IsInside(clrsF,"HI05X1_u");
      iY = IsInside(clrsF,"HI04X1_u");
      for (int i = 0; i<(int)iX.size(); i++ ) {
	for (int j = 0; j<(int)iY.size(); j++) {
	  x1 = clrsF[iY[j]]->getDigitsList().front()->getAddress();
	  x2 = clrsF[iX[i]]->getDigitsList().front()->getAddress();
	  ACCEPT = MIu->Accept(x1,x2);
	  if (hist>0) {
	    hIu[0]->Fill(x2,x1); if (ACCEPT) hIu[1]->Fill(x2,x1);
	  }
	  if (ACCEPT) break;
	}
	if (ACCEPT) break;
      }
    }

    else if (cT=='M') {  //   ********** MIDDLE HODO SYSTEM **********
      //                                                       ***** SAME W/ HMs
      vector<int> iX = IsInside(clrsF,"HM05X1_d");
      vector<int> iY = IsInside(clrsF,"HM04X1_d");
      for (int i = 0; i<(int)iX.size(); i++) {
	for (int j = 0; j<(int)iY.size(); j++) {
	  x1 = clrsF[iY[j]]->getDigitsList().front()->getAddress();
	  x2 = clrsF[iX[i]]->getDigitsList().front()->getAddress();
	  ACCEPT = MMXd->Accept(x1,x2);
	  if (hist>0) {
	    hMXd[0]->Fill(x2,x1); if (ACCEPT) hMXd[1]->Fill(x2,x1);
	  }
	  if (ACCEPT) break;
	}
	if (ACCEPT) break;
      }
	    
      bool ACCEPT1(false);   
      iX = IsInside(clrsF,"HM05Y1_d");
      iY = IsInside(clrsF,"HM04Y1_d");
      for (int i = 0; i<(int)iX.size(); i++) {
	for (int j = 0; j<(int)iY.size(); j++) {
	  x1 = clrsF[ iY[j] ]->getDigitsList().front()->getAddress();
	  x2 = clrsF[ iX[i] ]->getDigitsList().front()->getAddress();
	  ACCEPT1 = MMYd->Accept(x1,x2);
	  if (hist>0) {
	    hMYd[0]->Fill(x2,x1); if (ACCEPT1) hMYd[1]->Fill(x2,x1);
	  }
	  if (ACCEPT1) break;
	}
	if (ACCEPT1) break;
      }
      if (ACCEPT && ACCEPT1) ACCEPT = true;
      else                   ACCEPT = false;
	    
      if (!ACCEPT) {     
	iX = IsInside(clrsF,"HM05X1_u");
	iY = IsInside(clrsF,"HM04X1_u");
	for (int i = 0; i<(int)iX.size(); i++) {
	  for (int j = 0; j<(int)iY.size(); j++) {
	    x1 = clrsF[ iY[j] ]->getDigitsList().front()->getAddress();
	    x2 = clrsF[ iX[i] ]->getDigitsList().front()->getAddress();
	    ACCEPT = MMXu->Accept(x1,x2);
	    if (hist>0) {
	      hMXu[0]->Fill(x2,x1); if (ACCEPT) hMXu[1]->Fill(x2,x1);
	    }
	    if (ACCEPT) break;
	  }
	  if (ACCEPT) break;
	}
	ACCEPT1 = false;
	iX = IsInside(clrsF,"HM05Y1_u");
	iY = IsInside(clrsF,"HM04Y1_u");
	for (int i = 0; i<(int)iX.size(); i++) {
	  for (int j = 0; j<(int)iY.size(); j++) {
	    x1 = clrsF[ iY[j] ]->getDigitsList().front()->getAddress();
	    x2 = clrsF[ iX[i] ]->getDigitsList().front()->getAddress();
	    ACCEPT1 = MMYu->Accept(x1,x2);
	    if (hist>0) {
	      hMYu[0]->Fill(x2,x1); if (ACCEPT1) hMYu[1]->Fill(x2,x1);
	    }
	    if (ACCEPT1) break;
	  }
	  if (ACCEPT1) break;
	}
	if (ACCEPT && ACCEPT1) ACCEPT = true;
	else                   ACCEPT = false;
      }
    }

    else if (cT=='L') {  //   ********** LADDER HODO SYSTEM **********
      //                                                       ***** SAME W/ HLs
      vector<int> iX = IsInside(clrsF,"HL05X1_m");
      vector<int> iY = IsInside(clrsF,"HL04X1_m");
      for (int i = 0; i<(int)iX.size(); i++) {
	for (int j = 0; j<(int)iY.size(); j++) {
	  x1 = clrsF[iY[j]]->getDigitsList().front()->getAddress();
	  x2 = clrsF[iX[i]]->getDigitsList().front()->getAddress();
	  ACCEPT = ML->Accept(x1,x2);
	  if (hist>0) {
	    hL[0]->Fill(x2,x1); if (ACCEPT) hL[1]->Fill(x2,x1);
	  }
	  if (ACCEPT) break;
	}
	if (ACCEPT) break;
      }
    }

    else if (cT=='O') {
      //                                                       ***** SAME W/ HOs
      vector<int> iX = IsInside(clrsF,"HO04Y1_m");
      vector<int> iY = IsInside(clrsF,"HO03Y1_m");
      for (int i = 0; i<(int)iX.size(); i++) {
	for (int j = 0; j<(int)iY.size(); j++) {
	  x1 = clrsF[ iY[j] ]->getDigitsList().front()->getAddress();
	  x2 = clrsF[ iX[i] ]->getDigitsList().front()->getAddress();
	  ACCEPT = MO1->Accept(x1,x2);
	  if (hist>0) {
	    hO1[0]->Fill(x2,x1); if (ACCEPT) hO1[1]->Fill(x2,x1);
	  }
	  if (ACCEPT) break;
	}
	if (ACCEPT) break;
      }    
      iX = IsInside(clrsF,"HO04Y2_m");
      iY = IsInside(clrsF,"HO03Y1_m"); 
      for (int i = 0; i<(int)iX.size(); i++) {
	for (int j = 0; j<(int)iY.size(); j++) {
	  x1 = clrsF[ iY[j] ]->getDigitsList().front()->getAddress();
	  x2 = clrsF[ iX[i] ]->getDigitsList().front()->getAddress();
	  ACCEPT = MO2->Accept(x1,x2);
	  if (hist>0) {
	    hO2[0]->Fill(x2,x1); if (ACCEPT) hO2[1]->Fill(x2,x1);
	  }
	  if (ACCEPT) break;
	}
	if (ACCEPT) break;
      }
    }
  }
  return ACCEPT;
}




/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// CsPIDdoMuonIDstat ///////////////////////////////

CsPIDdoMuonIDstat* CsPIDdoMuonIDstat::instance_ = 0;

CsPIDdoMuonIDstat::CsPIDdoMuonIDstat()
{
  NR = 20;
  NC = 11;
  CsRegistry reg;
  reg.EOJRegistration( this );
  instance_=this;
  for(int i=0;i<NR;i++) {
    all[i]       = 0;
    allLastEv[i] = 0;
    for(int j=0;j<NC;j++) {
      stat[i][j]       = 0;
      statLastEv[i][j] = 0;
    }
  }
}

CsPIDdoMuonIDstat::~CsPIDdoMuonIDstat()
{
  instance_=0;
}


CsPIDdoMuonIDstat* CsPIDdoMuonIDstat::Instance() {
  if( instance_ == 0 ) new CsPIDdoMuonIDstat();
  return instance_;
}

void CsPIDdoMuonIDstat::Stat( int ir, unsigned int TM )
{
  if( ir >= NR ) {
    cout<<"ERROR!! CsPIDdoMuonIDstat::Stat: ir="<<ir<<endl;
    exit(1);
  }
  
  for(int i=0;i<NR;i++) for(int j=0;j<NC;j++) {statLastEv[i][j]=0;}
  all[ir]++;
  allLastEv[ir]=1;
  
  unsigned int Mask = 2047;
  TM=(TM&Mask);
  for(int i=0;i<NC;i++) {
    
    //unsigned int tm = (1<<i);
    unsigned int num = (TM&(1<<i));
    bool OR  = (num!=0 );
    bool AND = (num==TM);
    
    if( ir == 0 ) {
      if( OR  ) {stat[0][i]++; statLastEv[0][i]=1; }
      if( AND ) {stat[6][i]++; statLastEv[6][i]=1; }
    } else if( ir == 1 ) {
      if( OR  ) {stat[1][i]++; statLastEv[1][i]=1; }
      if( AND ) {stat[7][i]++; statLastEv[7][i]=1; }
    } else if( ir == 2 ) {
      if( OR  ) {stat[2][i]++; statLastEv[2][i]=1; }
      if( AND ) {stat[8][i]++; statLastEv[8][i]=1; }
    } else if( ir == 3 ) {
      if( OR  ) {stat[3][i]++; statLastEv[3][i]=1; }
      if( AND ) {stat[9][i]++; statLastEv[9][i]=1; }
    } else if( ir == 4 ) {
      if( OR  ) {stat[4][i]++; statLastEv[4][i]=1; }
      if( AND ) {stat[10][i]++;statLastEv[10][i]=1;}
    } else if( ir == 5 ) {
      if( OR  ) {stat[5][i]++; statLastEv[5][i]=1; }
      if( AND ) {stat[11][i]++;statLastEv[11][i]=1;}
    } else if( ir > 11 ) {
      if( AND ) {stat[ir][i]++;statLastEv[ir][i]=1;}
    }
    
  }
}


void CsPIDdoMuonIDstat::StatMW1( int ir, unsigned int TM )
{
  if( ir >= NR ) {
    cout<<"ERROR!! CsPIDdoMuonIDstat::Stat: ir="<<ir<<endl;
    exit(1);
  }
  
  if( allLastEv[ir]==0 ) all[ir]++;
  
  unsigned int Mask = 2047;
  TM=(TM&Mask);
  for(int i=0;i<NC;i++) {
    
    //unsigned int tm = (1<<i);
    unsigned int num = (TM&(1<<i));
    bool OR  = (num!=0 );
    bool AND = (num==TM);
    
    if( ir == 1 ) {
      if( OR  && statLastEv[1][i]==0 ) stat[1][i]++;
      if( AND && statLastEv[7][i]==0 ) stat[7][i]++;
      
      if( OR  && statLastEv[1][i]==1 && statLastEv[2][i]==0 ) stat[2][i]++;
      if( AND && statLastEv[7][i]==1 && statLastEv[8][i]==0 ) stat[8][i]++;
    }
    
  }
}


bool CsPIDdoMuonIDstat::end()
{ 
  //unsigned int Mask = 2047;
  unsigned int Mask = 1823;
  unsigned int N = CsEvent::Instance()->getNumberOfProcessedEvents();
  cout<<setprecision(3)<<endl;
  cout<<"/-----------------------------/ mu' ID statistics \\---------------------------\\"<<endl;
  cout<<"| The number of events (def)             "<<setw(7)<<N<<"       "<<setw(6)<<100<<"%                |"<<endl;
  cout<<"|----------------------+------------------------------------------------------|"<<endl;
  cout<<"| Trigger bin          |";
  cout<<"  all";
  for(int j=0;j<NC;j++) if( (Mask&(1<<j)) != 0 ) cout<<setw(6)<<(j+1);
  cout<<" |"<<endl;
  cout<<"|______________________|______________________________________________________|"<<endl;
  cout<<"|----------------------+------------------------------------------------------|"<<endl;
  
  for(int i=0;i<12;i++) {
         if( i == 0 ) cout<<"| Was called    OR (%) |";
    else if( i == 1 ) cout<<"| mu'           OR (%) |";
    else if( i == 2 ) cout<<"| mu'  N>1      OR (%) |";    
    else if( i == 3 ) cout<<"| mu'  Z < 1    OR (%) |";
    else if( i == 4 ) cout<<"| mu'  Z = 1    OR (%) |";
    else if( i == 5 ) cout<<"| mu'  Z = 2    OR (%) |";
    
    else if( i == 6 ) cout<<"| Was called   AND (%) |     ";
    else if( i == 7 ) cout<<"| mu'          AND (%) |     ";
    else if( i == 8 ) cout<<"| mu'  N>1     AND (%) |     ";
    else if( i == 9 ) cout<<"| mu'  Z < 1   AND (%) |     ";
    else if( i == 10) cout<<"| mu'  Z = 1   AND (%) |     ";
    else if( i == 11) cout<<"| mu'  Z = 2   AND (%) |     ";
    
    else continue;
    
    if( i>=0 && i<=5 ) {
      float val = 100*all[i]/N;
           if( val < 1  ) cout<<setprecision(1);
      else if( val < 10 ) cout<<setprecision(2);
      else                cout<<setprecision(3);
      cout<<setw(5)<<val;
    }
    
    for(int j=0;j<NC;j++) {
      if( (Mask&(1<<j)) != 0 ) {
	float val = 100*stat[i][j]/N;
	     if( val < 1  ) cout<<setprecision(1);
	else if( val < 10 ) cout<<setprecision(2);
	else                cout<<setprecision(3);
	cout<<setw(6)<<val;
      }
    }
    cout<<" |"<<endl;
    if( i == 2 ) cout<<"|----------------------+------------------------------------------------------|"<<endl;
    if( i == 8 ) cout<<"|----------------------+------------------------------------------------------|"<<endl;
    if( i == 5 ) {
      cout<<"|______________________|______________________________________________________|"<<endl;
      cout<<"|----------------------+------------------------------------------------------|"<<endl;
    }
  }
  cout<<"\\______________________|______________________________________________________/"<<endl;
  cout<<endl;

  return true;
}
