// $Id: TEvParaxialPR.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
  \brief Pattern Recognition for ``scattered mu'' (any very forward track in
  fact). Specific to "lattice".
  */

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cassert>
#include "CsErrLog.h"
#include "CsZone.h"
#include "CsCluster.h"
#include "CsDetector.h"
#include "CsGEMDetector.h"
#include "CsEventUtils.h"
#include "TWatches.h"
#include "Traffic.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "TAlgo2.h"
#include "TConstants.h"
#include "TLattice.h"

using namespace std;

void TEv::ParaxialPR()
{
  if (TOpt::ReMode[4]!=0)      // *************** IDEAL PR (MC) ***************
    return;
  const TSetup& setup = TSetup::Ref();
  if (setup.vIplFirst().size()!=5)  // ***** REQUIRE THE 5 STANDARD ZONES *****
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "%d spectrometer zones defined while 5 are expected!",
		  setup.vIplFirst().size());
  Traffic::Ref().Stopwatch.Start(16);

  TAlgo2 algo;
  // igroup = GROUP of DETECTORS, !=ZONE!
  //   - 1st group is all vertical planes in [target,muWall] + planes beyond
  //    wall w/in reach of particles through the hole of the wall. 1st
  //    is search for proj only, once and for all.
  //   - 3rd group is zone [SM2,muWall] + planes beyond w/in reach of hole.
  int igroup;
  // iset = Set, as opposed to group, of detectors.
  //   - It's used by TAlgo methods to discriminate among different groups
  //    spanning a same zone.
  //   - And it's index in PRpar's for given group.
  int iset;
  int ipl;   // ipl  = Plane index

  // Considered, for proj. or space search:
  //    16   FI03:08
  //  + 12   MM01:03
  //  + 44   GM01:11
  //  +  4   PS01
  //  + 21   PA01:11
  //  +  2   HI04
  //    30   ST01:06
  //  + 24   DW
  //  +  9   PB
  //   162
  const int maxpl   = 162, maxhit  = 1620;

   // ******************** PATTERNS of PLANE IDs CONSIDERED ********************
  const int maxplTot = 400, nPlPats = maxplTot/32+1;
  static int plPats[nPlPats];  // Pattern of planes used in proj search
  static int plPat2[nPlPats];  // Pattern of planes used in space search
  static int plPat3[nPlPats];  // Pattern of planes used in pick-up search
  static bool first = true; if (first) {
    first = false;
    // Avoid drift-like detectors, in order to be able to bypass ReLR
    const string FI = "FI", GM = "GM", HI = "HI", PA = "PA";
    const string MM = "MM", PB = "PB";
    const string *TBs[] ={ &FI,   &GM,  &HI,  &PA,  &MM,   &PB};
    const int    i0s[]  ={0xfc, 0x7ff, 0x08, 0x420, 0x4,    0};
    const int    i02[]  ={0xfc, 0x7ff, 0x08, 0x430, 0x7, 0x3f, 0x1, 0x3f,  0x0};
    // In fact, no (yet?) extra plane for pickup w.r.t. track search.
    // Not even HI05:
    //  - That detector has been removed from track search (although it
    //   is definitely relevant for the paraxial case) for fear its inclusion
    //   might skew the segment track too much and degrade its chi2.
    //  - It can (has to) be searched for in a later, pick-up, step
    const int    i03[]  ={0xfc, 0x7ff, 0x08, 0x43f, 0x7, 0x3f, 0x1, 0x3f,  0x0};
    int itb, match;
    for (igroup = 0, memset(plPats,0,nPlPats*sizeof(int)),
	   memset(plPat2,0,nPlPats*sizeof(int)); igroup<4; igroup++) {
      for (ipl = setup.vIplFirst()[igroup]; ipl<=setup.vIplLast()[igroup];
	   ipl++) {
	if (ipl>maxplTot)
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
			"# of planes in zones [0,3] > %d!",maxplTot);
	const TPlane  &p = setup.vPlane(ipl);
	if (p.IFlag==0) continue;                 // PLANE IS OFF...
	const TDetect &d = setup.vDetect(p.IDetRef);
	const string &TB = d.Name.substr(0,2);
	for (itb=match = 0; itb<int(sizeof(TBs)/sizeof(char*)); itb++)
	  if (TB==*TBs[itb]) { match = 1; break; }
	if (!match) continue;
        const char *station = d.Name.c_str()+2; char *end, **endptr = &end;
	int i = strtol(station,endptr,10);
	if (*endptr!=station+2*sizeof(char))
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
			"Cannot decode TBname \"%s\"!\n",d.Name.c_str());
	int j = 1<<i-1;
	if (match && (j&i0s[itb])) plPats[ipl/32] |= 1<<ipl%32;
	if (match && (j&i02[itb])) plPat2[ipl/32] |= 1<<ipl%32;
	if (match && (j&i03[itb])) plPat3[ipl/32] |= 1<<ipl%32;
      }
    }
  }

  int maxproj = setup.vProj().size(); // all proj.

  // working arrays
  int   jpls[1000];       // Index of generic plane in list of active planes
  float xpl [maxpl];               // X of plane
  int   idpl[maxpl];               // Type or index of plane in vDetect
  float res [maxpl];               // Resolution of plane
  float tol [maxpl];               // tolerances for track finding of plane
  int    fh [maxpl], lh[maxpl];    // first and after-the-last hit number on plane
  float cosa[maxpl];               // cos(a)
  float sina[maxpl];               // sin(a) of planes
  float Uhit[maxhit];              // Coord of hit
  int   href[maxhit];              // ref to vHit
  int   ifl [maxhit];              // Flag used hit
  float hinfos[maxhit];

  static int   hits[(TConstants_NTtrack_max+1)*maxpl];  // Found hits

  // ********** ARRAYS TO STORE RESULTS OF TRACK FINDING IN PROJs **********
  int Ntk_prj [4]; int iprojZ = -1;
  static float Y0_prj[4][TConstants_NTtrack_max], Yp_prj[4][TConstants_NTtrack_max];
#ifdef QUADRATIC_PROJ
  static float Y2_prj[TConstants_NTtrack_max];
#else
  float *Y2_prj = 0;
#endif

  float X0 = (float)setup.iPlane2Detect(setup.vIplFirst()[0]).X(0);
  for (int igroup= -1; igroup<3; igroup++) {

    // ******************** LOOP ON ZONES ********************
    // ********** 1ST LOOP = ALL-OUT VERTICAL SEARCH
    // ********** GROUP #3 = ZONE 0x4 + 0x8

    float X0_prv = X0;   int ipl0, ipl1;
    if      (igroup==-1) {
      X0 = (float)setup.iPlane2Detect(setup.vIplFirst()[0]).X(0);
      ipl0 = setup.vIplFirst()[0];      ipl1 = setup.vIplLast()[3];
      iset = 5;
    }
    else if (igroup==2) {
      X0 = (float)setup.iPlane2Detect(setup.vIplFirst()[2]).X(0);
      ipl0 = setup.vIplFirst()[2];      ipl1 = setup.vIplLast()[3];
      iset = 8;
    }
    else {
      X0 = (float)setup.iPlane2Detect(setup.vIplFirst()[igroup]).X(0);
      ipl0 = setup.vIplFirst()[igroup]; ipl1 = setup.vIplLast()[igroup];
      iset = igroup+6;
    }

    int nproj = 0, ndisc = 0;
    int projind[maxproj];   // Array of indices of projections
    int discind[maxproj];   // Array of indices of discarded projections

    int Ntk, nhit, npl, ier;

    // *************** COUNTING PLANES IN proj... ***************
    // *************** ...AND TOTAL #HITS, SINES AND COSINES ... ***************
    nhit = 0;
    float cos_prj[maxproj], sin_prj[maxproj];  // Cos(a) and Sin(a) of proj
    for (int ipr = 0; ipr< maxproj; ipr++) {

      // ********** LOOP ON ALL PROJ'S **********
      int npl = 0, jpr = -1, mhits = 0, nscifis = 0;
      float cos_mean=0, sin_mean=0, w_mean=0;  // Mean cos/sin in proj
      for (ipl = ipl0; ipl<=ipl1; ipl++) {

	// ***** LOOP ON PLANES IN GROUP *****
	const TPlane &p = setup.vPlane(ipl);
	if (!(plPats[ipl/32]&1<<ipl%32)) continue;   // Plane is not considered
	if (npl==0) {
	  int yet, i;
	  for (i=yet = 0; i<nproj; i++) if (p.IProj==projind[i]) {
	    yet = 1; break;
	  }
	  for (i     = 0; i<ndisc; i++) if (p.IProj==discind[i]) {
	    yet = 1; break;
	  }
	  if (!yet) jpr = p.IProj;
	}
	if (p.IProj!=jpr) continue;                // Not current projection
	const TDetect &d = setup.vDetect(p.IDetRef);
	if (nproj==0 && fabs(d.Ca)>.04)  // 1st = Z proj (using a looser cut...
	  continue;  // ...to make for misorientation of individual planes)
	mhits += p.vHitRef().size(); 
	cos_mean += d.Ca*d.Resol*d.Resol; sin_mean += d.Sa*d.Resol*d.Resol;
	w_mean += d.Resol*d.Resol;
	npl++; if (d.IType==22) nscifis++;
      }
      if (npl<TOpt::iPRpar[iset*10] &&      // NOT ENOUGH PLANES IN THIS PROJ...
	  nscifis<2) {                      // ...and not enough scifis
	discind[ndisc++] = jpr; continue;
      }
      cos_prj[nproj] = cos_mean/w_mean; sin_prj[nproj] = sin_mean/w_mean;
      if (fabs(cos_prj[nproj])<.035)        // Y PROJECTION: 2 degrees margin
	iprojZ = nproj;
      projind[nproj++] = jpr; nhit += mhits;
    }  // End of loop over projections
    if (nproj>6)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "# of Projections in zone %d (=%d) > %d",igroup,nproj,6);
    if (iprojZ==-1) CsErrLog::mes(elFatal,"No Y Projection!");
    if (nhit>maxhit) {  // ********** IF TOO MANY HITS => GIVE UP **********
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "# of Hits in group %d > %d",igroup,maxhit);
      Traffic::Ref().Stopwatch.Stop(16); return;
    }

    // ************************************************************************
    // ************************ TRACK FINDING in PROJs ************************
    // ************************************************************************

    for (int ipr = 0; ipr<nproj; ipr++) {

      // *************** LOOP ON PROJ'S ***************
      if (igroup==-1) {
	if (ipr!=iprojZ) continue;    // ***** igroup==-1: Y PROJ ONLY *****
      }
      else if (ipr==iprojZ) continue; // ***** else: ALL BUT Y PROJ'S *****

      int iproj = projind[ipr], nhit=0, npl=0;

      vector<int>::const_iterator ihit;
      for (ipl = ipl0; ipl<=ipl1; ipl++) {
	const TPlane&  p = setup.vPlane(ipl);
	if (!(plPats[ipl/32]&1<<ipl%32)) continue;   // Plane is not considered
	if (iproj!=p.IProj) continue;                // Not current projection

	// ********** LOOP ON PLANES IN PROJ **********
	const TDetect& d = setup.vDetect(p.IDetRef);
	idpl[npl] = d.IType; xpl[npl] =d.X(0); res[npl] = d.Resol;

	fh [npl]=nhit;
	for (ihit = p.vHitRef().begin(); ihit!=p.vHitRef().end(); ihit++) {

	  THit& h = vecHit[*ihit];	  // ***** LOOP ON HITS *****
	  Uhit[nhit]=h.U; href[nhit++]=(*ihit);
	}  // end of loop over hits on the plane
	lh[npl]=nhit;

	if (++npl==maxpl)
	  CsErrLog::msg
	    (elFatal,__FILE__,__LINE__,
	     "# of Planes in Proj (%.3f,%.3f) in Group %d (=%d) > %d",
	     cos_prj[ipr],sin_prj[ipr],igroup,npl,maxpl);
      }   // End of loop on planes
      if (TOpt::Print[5]) {
	printf("Det. group %1u    Proj %2u   Nplanes %3u   NHits %5u  ",
	       igroup,iproj,npl,nhit); cout<<flush;
      }
    
      // ********** RECONSTRUCTION IN ONE PROJ **********

      int Ntk = TConstants_NTtrack_max, ier(0);
      for (int jpl = 0; jpl<npl; jpl++) {           //***** SET TOLERANCES...
	tol[jpl] = res[jpl]*TOpt::dPRpar[iset*10+0];
	// => Add extra tolerance to accomodate curvature, independent
	//   upon detector resolution
	tol[jpl] += TOpt::dPRpar[iset*10+2]*fabs(cos_prj[ipr]) +
	  TOpt::dPRpar[iset*10+3]*fabs(sin_prj[ipr]);
      }

#ifdef QUADRATIC_PROJ
      float *Y2 = igroup==-1 ? Y2_prj : 0;  // In order to accomodate MS
#else
      float *Y2 = 0;
#endif
      ier = algo.FindProjQ(iset,X0,npl,xpl,tol,fh,lh,nhit,Uhit,
			   Ntk,hits, Y0_prj[ipr],Yp_prj[ipr],
			   0,Y2,idpl,href);
      if (ier>1) {
	Traffic::Ref().Stopwatch.Stop(16); return;
      }
      Ntk_prj[ipr] = Ntk;

    } // end of loop over pat. rec.  projections in this group


    if (TOpt::ReMode[5]!=0) {
      Traffic::Ref().Stopwatch.Stop(16); return;   // Space track finding is OFF
    }
    if (igroup==-1) continue;

    // *****************************************************************
    // ******************** TRACK FINDING IN SPACE ********************
    // *****************************************************************

    for (int it = 0; it<Ntk_prj[0]; it++)  // Refer Z tracks to current zone
      Y0_prj[0][it] += (X0-X0_prv)*Yp_prj[0][it];

    nhit = 0; npl = 0; Ntk = TConstants_NTtrack_max;
    vector<int>::const_iterator ihit;
    for (ipl = ipl0; ipl<=ipl1; ipl++) {              // Loop over all planes
      const TPlane&  p = setup.vPlane(ipl);
      if (!(plPat2[ipl/32]&1<<ipl%32)) {              // Plane is not considered
	jpls[ipl-ipl0] = -1; continue;               
      }
      const TDetect& d = setup.vDetect(p.IDetRef);
#if defined SCIFIS_ENHANCED && defined PR_SKIP_NON_ENHANCED
      if (iter<SCIFIS_ENHANCED-1 &&   // If SCIFIS_ENHANCED inactive @ "iter"...
	  (TOpt::iCut[14]&1<<igroup) &&    // ...and enabled for current zone
	  d.IType==22) {
	jpls[ipl-ipl0] = -1; continue;                   // ... Exclude scifis
      }
#endif

      idpl[npl] = p.IDetRef; xpl[npl] = d.X(0);
      cosa[npl] = d.Ca; sina[npl] = d.Sa;
      res[npl]=d.Resol;
      fh [npl]=nhit;
      jpls[ipl-ipl0] = npl;
      for (ihit = p.vHitRef().begin(); ihit!=p.vHitRef().end(); ihit++){

	// loop over hit references on the plane

	THit& h = vecHit[*ihit]; // ref to Thit vector element
	Uhit[nhit]=h.U;
	if (d.IType==22) {
	  double time;
	  if (!h.ptrCl->getTime(time)) hinfos[nhit] = 1.e6;
	  else                         hinfos[nhit] = time;
	}
	href[nhit]=(*ihit);
	if(++nhit==maxhit){
	  if (TOpt::Print[0]) 
	    cout<<"TEv::PrePattern ==> (Space).  More then "<<maxhit<<" hits in group "<<igroup<<endl;
	  Traffic::Ref().Stopwatch.Stop(16); return;
	}
      }  // end of loop over hits on the plane
      lh[npl]=nhit;
      if(++npl == maxpl){
	cout<<"TEv::PrePattern ==> (Space).  More then "<<maxpl<<" planes in group "<<igroup<<endl;
	Traffic::Ref().Stopwatch.Stop(16); return;
      }
    }   // end of loop over planes

    for (int jpl = 0; jpl<npl; jpl++) {           //***** SET TOLERANCES...
      tol[jpl] = res[jpl]*TOpt::dPRpar[igroup*10+0];
      // => Add extra tolerance to accomodate multiscattering, independent
      //   upon detector resolution
      tol[jpl] += TOpt::dPRpar[igroup*10+2]*fabs(cosa[jpl]) +
	TOpt::dPRpar[igroup*10+3]*fabs(sina[jpl]);
    }

#ifdef QUADRATIC_PROJ
    float *Y2 = igroup==1 ? Y2_prj : 0;
#else
    float *Y2 = 0;
#endif
    ier = algo.FindSpace(iset, 0,
			 X0, npl, idpl, cosa, sina, tol, fh, lh, nhit, Uhit, 
			 nproj, Ntk_prj, cos_prj, sin_prj, Y0_prj, Yp_prj, Y2,
			 hinfos,
			 Ntk, hits, ifl);
    if (ier<0)  {
      Traffic::Ref().Stopwatch.Stop(16); return;
    }

    // *************** STORE TRACKS ***************
  
    for (int it = 0, nh = 0; it<Ntk; it++) { // ***** LOOP ON FOUND TRACKS *****
      TTrack t;
      t.Type = 1<<igroup;   // TTrack's bit pattern of zones
      t.Scifi = 1<<igroup;  // TTrack's bit pattern of scifi-only

      // ***** STORE HITS *****
      int ih0 = nh;
      int npats = (ipl1-ipl0+1)/32+1; unsigned int hit_pats[npats];
      bool modified = false; int n_scifis = 0;
      for (int ipat = 0; ipat<npats; ipat++) hit_pats[ipat] = 0;
      while(1){
	int iref = hits[nh++], ref;
	if      (iref==-1) break;                     // next track's hits
	else if (iref<-1)  ref = href[-2-iref];       // alternative placeholder
	else               ref = href[iref];          // normal hit
	if (ref>=int(vecHit.size()))          // Check hit or alternative
	  CsErrLog::mes
	    (elFatal,"space ==> something's wrong with hit references");
	if (iref<-1) continue;                // Skip alternative placeholder
	THit &h = vecHit[ref]; ipl = h.IPlane;
	const TPlane &p = setup.vPlane(ipl);
	const TDetect& d = setup.vDetect(p.IDetRef);
	if (d.IType==22) n_scifis++;
	else t.Scifi = 0;
	t.AddHit(h);
      }
      if (igroup==2 && ipl>=setup.vIplFirst()[3]) t.Type |= 0x8;

      static float QNchi2; // "static" to avoid "warning... uninitialized"
      int QNstatus = 0;
      static int QNmode; static unsigned int QNnh0;
      if (n_scifis>=4) {  // Ensures more than one scifi station
	QNmode = -1; QNnh0 = 5;
      }
      else if (igroup==0) {
	QNmode = 3; QNnh0 = 6;
      }
      else if (igroup<3) {
	QNmode = 5; QNnh0 = 6;
      }
      else QNmode = -1;

      static double incidence; // "static" to avoid "warning... uninitialized"

      int ok = 1;
      static int idebug = 0;
      // If too many (very good) hits have been canceled => give up...
      if (modified &&       // ...but avoid throwing away SCIFIS_ENHANCED
	  (int)t.NHits<TOpt::iPRpar[igroup*10 + 2]) { ok = 0; goto pickup;}

      if (idebug>1) {
	if (it==0) {
	  printf("\nFindSpace %d: %d tracks <\n",igroup,Ntk);
	  for (int iproj = 0; iproj<nproj; iproj++)
	    printf("Proj %d: %d\n",iproj,Ntk_prj[iproj]);
	}
	printf("\nNew track %d\n",it);
	t.DumpHits(0);
      }

      if (TOpt::ReMode[14] && (igroup<3 || igroup==3)) {
	if (t.QuickKF(-1,0)) {          // Straight Backward KF
	  t.QuickKF(1,0);
	  t.IFit |= 0x1;                // Set fit type bit pattern to |Straight

	  if (QNstatus && TOpt::ReMode[19]==0) {
	    if (QNmode!=0) QNmode = 4;
	    if (QNstatus==-1) {   // => Update "QNchi2"
	      t.QNewtonFit(1,QNmode);
	      QNchi2 = t.Chi2tot/(t.NDics-QNnh0); QNstatus = 2;
	    }
	    if (idebug>1) { cout<<"QN\n"; t.DumpHits(1); }

	    // ***** CLEANING *****

	    // Reference to "TLattice". In view of determining whether given plane
	    // has meaningfull guestimates. This is not very satisfying to have a
	    // refernce to "TLattice". Best would have been to have "TTrack" (which
	    // could had guestimates determined by other mean than "TLattice")
	    // return the answer. But for the time being...
	    const TLattice *lat = TLattice::Ptr();    


	    int iref, ih;
	    static int   QNnh; // "static" to avoid "warning... uninitialized"
	    if      (QNchi2<100) {

	      // ***** ALL (reasonably good) CHI2's => do only:
	      // - building hit pattern
	      // - raising alternative ambiguity
	      ok = 2;

	      for (ih = ih0, iref = 0; ih<nh-1; ih++) {
		int ialt = iref;  // Remember previous hit: might be alternative
		iref = hits[ih];
		if (iref<0) continue;          // Bypass alternative placeholder
		THit& h = vecHit[href[iref]];
		ipl = h.IPlane-ipl0;
		if      (lat->dico_idx[h.IPlane]==-1) continue;
		else if (lat->dico_idx[h.IPlane]==-2) break;
		hit_pats[ipl/32] |= 1<<ipl%32;
		if (ialt<0 && ialt!=-2-maxhit      // Alternative exists...
		    && QNchi2<5) {                 // ...and good enough chi2
		  ialt = -2-ialt;
		  if (ialt==maxhit) continue;
		  THit& halt = vecHit[href[ialt]];
		  float guess = t.vGuests[ipl];
		  float residu = fabs(guess-h.U), residualt = fabs(guess-halt.U);
		  CsCluster *ca = halt.ptrCl->getMirrorCluster();
		  if (ca) {
		    float residum = fabs(guess-ca->getU()/10);
		    if (residum<residualt) residualt = residum;
		    else ca = NULL;
		  }
		  if (residualt<residu) {
		    t.SubHit(h);
		    if (ca) {
		      halt.u = ca->getU()/10; halt.ptrCl = ca; halt.Rotate();
		    }
		    t.AddHit(halt);
		    QNstatus = -1; ifl[iref]--; ifl[ialt]++; // Update flags
		    hits[ih] = ialt;      // Replace hit with alternative
		    hits[ih-1] = -2-maxhit;
		  }
		}
	      }
	      QNnh = ih;           // Save "nh" in QN range
	      if (QNstatus==-1) {
		t.QNewtonFit(1,QNmode);
		QNstatus = 2;  QNchi2 = t.Chi2tot/(t.NDics-QNnh0);
		if (idebug==3) { printf("Alt\n"); t.DumpHits(1); }
		t.IFit &= 0x8;       // Flag straight line fit to be updated
	      }
	    }
	    if (1.5<QNchi2 && QNchi2<100) {

	      // ***** BAD CHI2

	      int modified = 0; //TTrack &t_bak = t;
	      for (int kter = 0; kter<2; kter++) {
		static TStation *sPrv; static unsigned int sPat;
		for (ih = ih0, iref = 0, sPrv = 0; ih<QNnh; ih++) {
		  int ialt = iref;// Remember previous hit: might be alternative
		  iref = hits[ih];
		  if (iref<0) continue;        // Bypass alternative placeholder
		  THit& h = vecHit[href[iref]];
		  ipl = h.IPlane-ipl0; int jpl = jpls[ipl];
		  if ((hit_pats[ipl/32]&(1<<ipl%32))==0) {
		    // Must have been erased in a previous kter
		    // Or, if "kter"==0, outside Dico range
		    continue;
		  }
		  float residu = (t.vGuests[ipl]-h.U)/res[jpl];
		  // Isolated?
		  int ipat = ipl/32;
		  unsigned int pat = hit_pats[ipat];
		  const TPlane &p = setup.vPlane(h.IPlane);
		  const TStation *&s = p.Station;
		  if (s!=sPrv) {
		    sPrv = const_cast<TStation*&>(s);
		    int jpat = s->IPat; sPat = pat>>s->JPl;
		    if (jpat<npats && s->JPl!=32 && s->JPl!=0)
		      sPat |= hit_pats[jpat+1]<<(32-s->JPl);
		    sPat &= s->Pat;
		  }
		  bool isolated = int(sPat&s->Pat)==(1<<ipl%32)>>s->JPl;
		  if (fabs(residu)>                   // IF VERY LARGE RESIDUAL
		      TOpt::dPRpar[igroup*10+0]*3 ||
		      kter==1 &&
		      (fabs(residu)>                  // IF LARGE RESIDUAL...
		       TOpt::dPRpar[igroup*10+0]
		       /* Allowance for fit imperfection */ * 1.75 ||
		       fabs(residu)>
		       TOpt::dPRpar[igroup*10+0]*0 && // IF !SO LARGE BUT ALSO
		       (ifl[iref]>1 ||                  // ...OR MULTIHIT
			ialt<0 && -2-ialt!=maxhit ||    // ...OR AMBIGUOUS
			isolated))) {                   // ...OR ISOLATED
		    if (idebug) printf("%d %f %f\n",
				       ipl+ipl0,residu,TOpt::dPRpar[igroup*10]);
		    t.SubHit(h);                           // SUBTRACT HIT
		    vector<float> vguestsp = t.vGuests;
		    if (ialt<0 && -2-ialt!=maxhit) {       // IF ALSO AMBIGUOUS
		      THit& halt = vecHit[href[-2-ialt]];
		      t.AddHit(halt);

		      t.QNewtonFit(1,QNmode);                       // Check
		      if (idebug==3) t.DumpHits(1);
		      float chi2 = t.Chi2tot/(t.NDics-QNnh0);
		      if (chi2<QNchi2*.8) {
			QNstatus = 2;
			ifl[iref]--; QNchi2 = chi2; modified = 1;
			ifl[-2-ialt]++;
			hits[ih] = -2-ialt; hits[ih-1] = -2-maxhit; // hit<-alternative
			continue;
		      }
		      else {
			t.SubHit(halt); QNstatus = 1; 
			t.vGuests = vguestsp;
		      }
		    }
		    if (t.NDics<=QNnh0) {
		      t.AddHit(h); QNstatus = 1;	    
		    }
		    else {
		      t.QNewtonFit(1,QNmode);                        // Check
		      if (idebug==3) t.DumpHits(1);
		      float chi2 = t.Chi2tot;
		      if (chi2<(QNchi2*(t.NDics-QNnh0+1)-residu*residu)*.8) {
			if ((int)t.NHits<TOpt::iPRpar[igroup*10 + 5]) {
			  // Not enough hits any more: Flag track to be released
			  ok = 0; goto next_step;
			}
			QNstatus = 2; hit_pats[ipl/32] ^= 1<<ipl%32;
			// ATTENTION! pourquoi NDics-6?
			ifl[iref]--; QNchi2 = chi2/(t.NDics-6); modified = 1;
			if (QNchi2<2)  goto next_step;
		      }
		      else {                                        // REadd hit
			t.AddHit(h); QNstatus = 1; hit_pats[ipl/32] |= 1<<ipl%32;
			t.vGuests = vguestsp;
		      }
		    }
		  }
		}
	      }
	    next_step:
	      if (modified)
		t.IFit &= 0x8;           // Flag straight line fit to be updated
              if (ok && QNchi2<20.) {
		if (QNstatus!=2) {         // Update fit guestimates...
		  t.QNewtonFit(1,QNmode);
		  QNchi2 = t.Chi2tot/(t.NDics-QNnh0); QNstatus = 2;
		  if (idebug==3) { printf("Clean\n"); t.DumpHits(1); }
		  if (1.75<QNchi2 && QNchi2<5) { 
		    for (ih = ih0, iref = 0; ih<nh-1; ih++) {
		      int ialt = iref;  // Remember previous hit: might be alternative
		      iref = hits[ih];
		      if (iref<0) continue;          // Bypass alternative placeholder
		      THit& h = vecHit[href[iref]];
		      ipl = h.IPlane-ipl0;
		      if ((hit_pats[ipl/32]&(1<<ipl%32))==0) continue;
		      if      (lat->dico_idx[h.IPlane]==-1) continue;
		      else if (lat->dico_idx[h.IPlane]==-2) break;
		      if (ialt<0 && ialt!=-2-maxhit) {   // Alternative exists...
			ialt = -2-ialt;
			if (ialt==maxhit) continue;
			THit& halt = vecHit[href[ialt]];
			float residu    = t.vGuests[ipl]-h.U;
			float residualt = t.vGuests[ipl]-halt.U;
			if (fabs(residualt)<fabs(residu)) {
			  t.SubHit(h); t.AddHit(halt);
			  QNstatus = -1; ifl[iref]--; ifl[ialt]++; // Update flags
			  hits[ih] = ialt;      // Replace hit with alternative
			}
		      }
		    }
		    if (QNstatus!=2) {
		      t.QNewtonFit(1,QNmode);
		      QNchi2 = t.Chi2tot/(t.NDics-QNnh0); QNstatus = 2;
		      if (idebug==3) { printf("Clean\n"); t.DumpHits(1); }
		    }
		  }
		}
              }
	      else ok = 0;
	    }
	    else if (QNchi2>=100) ok = 0;
	    if (idebug) { cout<<"Clean\n"; t.DumpHits(1); }
	  }
	  else {

	    // ********** KF TTracks (no QN) **********

#ifdef CUT_KF_UPON_0
	    // chi2 /=n (and not /=n-5) for fear n<=5
	    if (t.Chi2tot/t.NHits>10) ok = 0;
#endif
	    if (idebug) { cout<<"No QN\n"; t.DumpHits(0); }
	  }
	}
	else ok = 0;               // No KF Fit => Flag track to be released
      }
      else if (igroup>=3) { // Beam and muon tracks: a priori straight track...
	t.QuickKF(-1,0);       // ...base cleaning on straight KF
	float old_chi2 = t.Chi2tot/(t.NHits-4);
	if (old_chi2>10) {                      // ...bad chi2: and try cleaning
	  static THit *worst_hit; float worst_resid = 0;
	  double x0 = t.Hfirst(0), y0 = t.Hfirst(1), z0 = t.Hfirst(2),
	    yp = t.Hfirst(3), zp = t.Hfirst(4);
	  for (int ih = ih0; ih<nh; ih++) {
	    int iref = hits[ih];
	    if (iref<0) continue; // Bypass alternative placeholder
	    THit &h = vecHit[href[iref]];
	    ipl = h.IPlane-ipl0; int jpl = jpls[ipl]; float x = xpl[jpl]-x0;
	    float resid = ((y0+yp*x)*cosa[jpl]+(z0+zp*x)*sina[jpl]-h.U)/res[jpl];
	    //#define PR_DEBUG_STRAIGHT_CLEAN
#ifdef PR_DEBUG_STRAIGHT_CLEAN
	    printf(" %d %.2f",jpl,resid);
#endif
	    if (fabs(resid)>worst_resid) {
	      worst_resid = fabs(resid); worst_hit = &h;
	    }
	  }
#ifdef PR_DEBUG_STRAIGHT_CLEAN
	  printf("\n");
#endif
	  if (worst_resid>TOpt::dPRpar[igroup*10+0]*2) { // Large residu...
	    t.SubHit(*worst_hit); t.QuickKF(-1,0);         // Refit w/o it
	    if (t.Chi2tot/(t.NHits-4)<old_chi2*.5) {       // Big improvement...
	      if ((int)t.NHits<TOpt::iPRpar[igroup*10+5])
		ok = 0; // No more enough hits => Flag track to be released
	    }
	    else {
	      t.AddHit(*worst_hit); t.QuickKF(-1,0);       // ...else restore
	    }
	  }
	}
	t.QuickKF(1,0);                        // Forward straight KF...
	t.IFit |= 0x1; // Set fit type to |Straight
	if (igroup==4 && TOpt::dCut[2]) {      // Beam tracks: standalone...
	  if (t.Chi2tot/(t.NHits-4)>             // ...=> cut bad chi2
	      TOpt::dCut[2]) ok = 0;
	}
      }

    pickup:
      if (ok) {
	if (ok==2 &&                                 // ***** QN Fitted...
	    t.Chi2tot/(t.NDics-4)<TOpt::dCut[60]) {  // ...GOOD CHI2...

	    // *************** PICK-UP HITS ***************

	  // Reference to "TLattice". In view of determining whether given plane
	  // has meaningfull guestimates. This is not very satisfying to have a
	  // reference to "TLattice". Best would have been to have "TTrack" (which
	  // could have guestimates determined by other mean than "TLattice")
	  // return the answer. But for the time being...
	  const TLattice *lat = TLattice::Ptr();    

	  int added = 0;
	  for (ipl = ipl0; ipl<=ipl1; ipl++) {
	    // Use extended set of planes (not only SAS) in the search
	    if (!(plPat3[ipl/32]&1<<ipl%32)) continue;
	    if      (lat->dico_idx[ipl]==-1) continue;
	    else if (lat->dico_idx[ipl]==-2) break;
	    int kpl = ipl-ipl0; if ((hit_pats[ipl/32]&(1<<ipl%32))==0) {
	      const TPlane &p = setup.vPlane(ipl);
	      const TDetect &d = setup.vDetect(p.IDetRef);
	      if (igroup==1) {
		if (980<t.Hlast(0) && t.Hlast(0)<1050) {
		  // Zone #1: Special case of extrapolation over full
		  // length of RICH while track is also defined in ST04
		  // which is not in the dictionary
		  if (750<d.X(0) && d.X(0)<980) continue;
		}
	      }
	      // no hits associated with this track on the plane
	      float yy = t.Hfirst(1)+(d.X(0)-t.Hfirst(0))*t.Hfirst(3);
	      float zz = t.Hfirst(2)+(d.X(0)-t.Hfirst(0))*t.Hfirst(4);
	      if (!d.InActive(yy,zz)) continue;
	      float dUBest = d.Resol*TOpt::dPRpar[igroup*10+0]*
		/* Allowance for fit imperfection */ 1.75;
	      int jpl = jpls[kpl]; if (jpl>=0) { // Plane used in track search
		int ir, irBest;
		for (ir = fh[jpl], irBest = -1; ir<lh[jpl]; ir++) {
		  if (ifl[ir]) continue;
		  THit &h = vecHit[href[ir]];
		  float dU = fabs(t.vGuests[kpl]-h.U); if (dU<dUBest) {
		    dUBest = dU; irBest = ir;
		  }
		}
		if (irBest>=0) {
		  t.AddHit(vecHit[href[irBest]]); added++; ifl[irBest]++;
		}
	      }
	      else {                             // Plane in downstream LAS
		vector<int>::const_iterator ih, ihBest;
		for (ih = p.vHitRef().begin(), ihBest = p.vHitRef().end();
		     ih!=p.vHitRef().end(); ih++) {
		  THit &h = vecHit[*ih];
		  if (h.Status) continue;
		  float dU = fabs(t.vGuests[kpl]-h.U); if (dU<dUBest) {
		    dUBest = dU; ihBest = ih;
		  }
		}
		if (ihBest!=p.vHitRef().end()) {
		  t.AddHit(vecHit[*ihBest]); added++;
		}
	      }
	    }
	  }
	  if (added) {
	    t.QuickKF(-1,0); t.QuickKF(1,0);
	    if (idebug) { cout<<"Pick-up\n"; t.DumpHits(0); }
	  }
	}
	listTrack.push_back(t); // ***** STORE THE TRACK
      }

    }  // End of loop on found tracks

    // We don't flag used hit here (THit::Status=1). We wait for a
    // confirmation by bridging (viz. by "TEv::BridgeSegments" which is
    // expected to be called next (by "CsTrafficPrepattern::doPrepattern"))
    // thatcandidate track is reliable.

  }  // End of loop on det. groups

  Traffic::Ref().Stopwatch.Stop(16);

#ifdef DEBUG
  static int idebug = 0;
  if (idebug) DumpTrackList();
#endif
}
