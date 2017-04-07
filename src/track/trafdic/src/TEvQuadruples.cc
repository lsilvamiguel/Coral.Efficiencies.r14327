// $Id: TEvQuadruples.cc 13548 2012-08-03 13:49:24Z kdupraz $

/*!
  \brief Search for quadruple == 4 mutually consistent hits.
         Optimized for GEM stations, i.e. (Y,Z) + (U,V) w/
	 Z (resp. V) very close in abscissa to Y (resp. U) planes.
  */

//   ipl     - Index of 1st plane of station (which is herein referenced as
//            Y but could have any orientation)
//    hit[]  - hit coordinate in detector reference system (DRS)...
//   yhit[]   ... Y of (U,V) counterpart of (Y,Z) in quad
//   zhit[]   ... Z of ....
//   yref[]  - Ref to vHit of yhit
//   zref[]  - Ref to vHit of zhit
//   ifl[]   - Status flag of THit's: >=0 means corresponding THit is member of
//            a quadruple, the value being = quadruple#.

#include <cstdio>
#include "CsCluster.h"
#include "CsErrLog.h"
#include "CsHistograms.h"
#include "CsGEMDetector.h"
#include "CsPixelGEMDetector.h"
#include "TEv.h"
#include "TOpt.h"
#include "TSetup.h"
using namespace std;

int TEv::Quadruples(int ipl, float *yhit, float *zhit, int *yref, int *zref,
		    int ifl[])
{
  const double delta = TOpt::ReMode[17];

  // Orientations (to be filled upon first encountering corresponding planes)
  float cy, sy, cz, sz, cu, su, cv, sv;
  float szy, syu, szu, syv, szv, svu;

  // Number of hits per plane (which serve in turn to flag 1st encounters))
  int   npy = 0, npz = 0, npu = 0, npv = 0;
  vector<int>::const_iterator ihy, jhy, ihz0, ihz, ihu0, ihu, ihv0, ihv;

  const TSetup &setup = TSetup::Ref();

  const TPlane  &py = setup.vPlane(ipl);
  const TDetect &dy = setup.vDetect(py.IDetRef);
  if (dy.IType!=26 && dy.IType!=28)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "Detector #%d of type %d whereas type 26 or 28 expected!",ipl,dy.IType);

  // IDs of GEM stations, Current ID
  static map<int,int> ids; int id;
  map<int,int>::iterator iid;
  // Parametersisations of amplitude correlation
  const int nGEMs = 14;
  static float ayz0s[nGEMs], ayz1s[nGEMs], ayz2s[nGEMs], auv0s[nGEMs], auv1s[nGEMs], auv2s[nGEMs];
  float ayz0, ayz1, ayz2, auv0, auv1, auv2;

  //#define Quadruples_HISTO
#ifdef Quadruples_HISTO
  static CsHist1D *tGEMYZU[nGEMs], *tGEMYZUV[nGEMs];    // Spatial Correlation
  static CsHist1D *tGEMYZ[nGEMs], *tGEMUV[nGEMs];       // Amplitude Correlation
  static CsHist1D *tGEMok = 0, *tGEMerr = 0;
  static CsHist1D *tGEMt[nGEMs];
#endif

  if ((iid = ids.find(ipl))!=ids.end()) id = (*iid).second;
  else {

    // ******************** INITIALISATION ********************

    id = ids.size(); ids[ipl] = id;
    if (id>=nGEMs) CsErrLog::msg(elFatal,__FILE__,__LINE__,
				 "GEM # >= %d!",nGEMs);

    // ********** I) PARAMETERISATION oF AMPLITUDE CORRELATION **********

    CsGEMDetector *csG = dynamic_cast<CsGEMDetector*>(dy.PtrDet());
    CsPixelGEMDetector *csPG = dynamic_cast<CsPixelGEMDetector*>(dy.PtrDet());
    const float *ampCorr = 0;
    if      (csG)  ampCorr = csG->getAmpCorr();
    else if (csPG) ampCorr = csPG->getAmpCorr();
    else
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"ipl = %d = %s not a CsGEM nor a CsPixelGEM!",ipl,dy.Name.c_str());
    ayz0s[id] = ampCorr[0]; ayz1s[id] = ampCorr[1]; ayz2s[id] = ampCorr[2];
    const TPlane  &pu = setup.vPlane(ipl+2);
    const TDetect &du = setup.vDetect(pu.IDetRef);
    csG = dynamic_cast<CsGEMDetector*>(du.PtrDet());
    csPG = dynamic_cast<CsPixelGEMDetector*>(du.PtrDet());
    if      (csG)  ampCorr = csG->getAmpCorr();
    else if (csPG) ampCorr = csPG->getAmpCorr();
    else
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"ipl = %d = %s not a CsGEM nor a CsPixelGEM!",ipl=1,du.Name.c_str());
    auv0s[id] = ampCorr[0]; auv1s[id] = ampCorr[1]; auv2s[id] = ampCorr[2];

#ifdef Quadruples_HISTO
    // ********** II) HISTOGRAMS **********
    if (TOpt::Hist[6]) {
      if (!tGEMok) {
	tGEMok  = new CsHist1D("tGEMok", "GEM OK"   ,10,-.5,9.5);
	tGEMerr = new CsHist1D("tGEMerr","GEM Error",10,-.5,9.5);
      }
      char name[] = "tGEMYZUVij";
      char title[] = "GM01X1__ Ay-(-00.0+1.000*Az-5.00e-05*Az^{2})     ";
      sprintf(name,"tGEMYZU%d",id+1);
      sprintf(title,"%s U(YZ)-U",dy.Name.c_str());
      tGEMYZU[id] = new CsHist1D(name,title,150,-.5,.5);
      sprintf(name,"tGEMYZUV%d",id+1);
      sprintf(title,"%s V(YZU)-V",dy.Name.c_str());
      tGEMYZUV[id] = new CsHist1D(name,title,120,-.4,.4);
      sprintf(name,"tGEMYZ%d",id+1);
      sprintf(title,"%s Ay-(%5.1f+%5.3f*Az%9.3g*Az^{2})",
	      dy.Name.c_str(),ayz0s[id],ayz1s[id],ayz2s[id]);
      tGEMYZ[id] = new CsHist1D(name,title,250,-500,500);
      sprintf(name,"tGEMUV%d",id+1);
      sprintf(title,"%s Au-(%5.1f+%5.3f*Av%9.3g*Av^{2})",
	      du.Name.c_str(),auv0s[id],auv1s[id],auv2s[id]);
      tGEMUV[id] = new CsHist1D(name,title,250,-500,500);
      sprintf(name,"tGEMt%d",id+1);
      sprintf(title,"%s #sigmat",dy.Name.c_str());
      tGEMt[id] = new CsHist1D(name,title,250,0,60);
    }
#endif
  }

  // Parameterisation of amplitude correlation for current GEM station
  ayz0 = ayz0s[id]; ayz1 = ayz1s[id]; ayz2 = ayz2s[id];
  auv0 = auv0s[id]; auv1 = auv1s[id]; auv2 = auv2s[id];

  // ******************** Y PLANE ********************

  float resol = dy.Resol;
  npy = py.vHitRef().size();
  if (npy ) {
    ihy = py.vHitRef().begin();
    // Flags to store the status, yet quad member or not, of Y hits
    // In principle, the following should work...
    //for (int i = *ihy; i<*ihy+npy; i++) ifl[i] = -1;
    // ...but it may turn out that the ordering of the hits in TPlane:vHitRef
    // differs from that in TEv::vHit. (I don't understand exactly why, but have
    // observed the case when a hit from a pixelGEM of the strip kind delivers
    // indentical twin hits, one on each side.)
    for (jhy = ihy; jhy!=py.vHitRef().end(); jhy++) ifl[*jhy] = -1;
  }
  cy = dy.Ca; sy = dy.Sa;

  // ******************** Z PLANE ********************

  ipl++;
  const TPlane  &pz = setup.vPlane(ipl);
  const TDetect &dz = setup.vDetect(pz.IDetRef);
  if (dz.IType!=26 && dz.IType!=28)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "Detector #%d of type %d whereas type 26 or 28 expected!",ipl,dz.IType);
  npz = pz.vHitRef().size();
  if (npz) {
    ihz0 = pz.vHitRef().begin();
    // Flags to store the status, yet quad member or not, of Z hits
    // In principle, the following should work...
    //for (int i = *ihz0; i<*ihz0+npz; i++) ifl[i] = -1;
    // ...but does not always, cf. the Y plane case supra.
    for (ihz = ihz0; ihz!=pz.vHitRef().end(); ihz++) ifl[*ihz] = -1;
  }
  cz = dz.Ca; sz = dz.Sa;
  szy = sz*cy-cz*sy;

  // ******************** U PLANE ********************

  ipl++;
  const TPlane &pu = setup.vPlane(ipl);
  const TDetect &du = setup.vDetect(pu.IDetRef);
  if (du.IType!=26 && du.IType!=28)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "Detector #%d of type %d whereas type 26 or 28 expected!",ipl,du.IType);
  npu = pu.vHitRef().size();
  if (npu) {
    ihu0 = pu.vHitRef().begin();
    // Flags to store the status, yet quad member or not, of U hits
    // In principle, the following should work...
    //for (int i = *ihu0; i<*ihu0+npu; i++) ifl[i] = -1;
    // ...but does not always, cf. the Y plane case supra.
    for (ihu = ihu0; ihu!=pu.vHitRef().end(); ihu++) ifl[*ihu] = -1;
  }
  cu = du.Ca; su = du.Sa;
  syu = sy*cu-cy*su; szu = sz*cu-cz*su;

  // ******************** V PLANE ********************

  ipl++;
  const TPlane  &pv = setup.vPlane(ipl);
  const TDetect &dv = setup.vDetect(pv.IDetRef);
  if (dv.IType!=26 && dv.IType!=28)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "Detector #%d of type %d whereas type 26 or 28 expected!",ipl,dv.IType);
  npv = pv.vHitRef().size();
  if (npv) {
    ihv0 = pv.vHitRef().begin();
    // Flags to store the status, yet quad member or not, of V hits
    // In principle, the following should work...
    //for (int i = *ihv0; i<*ihv0+npv; i++) ifl[i] = -1;
    // ...but does not always, cf. the Y plane case supra.
    for (ihv = ihv0; ihv!=pv.vHitRef().end(); ihv++) ifl[*ihv] = -1;
  }
  cv = dv.Ca; sv = dv.Sa;
  syv = sy*cv-cy*sv; szv = sz*cv-cz*sv; svu = sv*cu-cv*su;
  double svy = sv*cy-cv*sy, svz = sv*cz-cv*sz;

  // 1 of the PLANES  EMPTY or OFF => Give up (But this takes
  // place still after all "ifl's" have been initialised)
  if (npy*npz*npu*npv==0) return 0;
  if (py.IFlag*pz.IFlag*pu.IFlag*pv.IFlag==0) return 0;

  typedef struct {
    int   jhs[4];
    float Y,V;
    float dV;
  } Quad;
  Quad quads[npy];  // No more quads than Y hits (because unambiguous quads)
  int   nquads = 0;

  // ******************** Y PLANE ********************

  for (; ihy!=py.vHitRef().end(); ihy++) {
    // loop over hit references on Y plane
    THit& hy = vecHit[*ihy]; float Y = hy.U;

    // ******************** Z PLANE ********************


    for (ihz = ihz0; ihz!=pz.vHitRef().end(); ihz++) {
      // loop over hit references on Z plane
      THit& hz = vecHit[*ihz]; float Z = hz.U;

      // ******************** U PLANE ********************

      float Uyz = (Y*szu-Z*syu)/szy;

      for (ihu = ihu0; ihu!=pu.vHitRef().end(); ihu++) {
	// loop over hit references on U plane
	THit& hu = vecHit[*ihu]; float U = hu.U;
	float dU = U-Uyz;
#ifdef Quadruples_HISTO
	if (TOpt::Hist[6]) {  // Fill histograms
	  if (!isMC || hy.IKine>=0 && hz.IKine==hy.IKine && hu.IKine==hy.IKine)
	    tGEMYZU[id]->Fill(dU);
	}
#endif
	if      (dU<-py.Station->URange) continue;
	else if (dU> py.Station->URange) break;

	// ******************** V PLANE ********************

	float Vyz = (Y*szv-Z*syv)/szy;

	for (ihv = ihv0; ihv!=pv.vHitRef().end(); ihv++) {
	  // loop over hit references on V plane
	  THit& hv = vecHit[*ihv]; float V = hv.U;
	  float dV = V-Vyz;
#ifdef Quadruples_HISTO
	  if (TOpt::Hist[6]) {  // Fill histograms
	    if (!isMC || hy.IKine>=0 && hz.IKine==hy.IKine &&
		hu.IKine==hy.IKine && hv.IKine==hy.IKine) {
	      tGEMYZUV[id]->Fill(dV);
	    }
	  }
#endif
	  if      (dV<-py.Station->VRange) continue;
	  else if (dV> py.Station->VRange) break;

	  // ******************** CANDIDATE QUAD ********************

	  // Whether any of the 4 component hits yet quad member
	  int iyzuv, iqalt, nalts, ialts, ambig;
	  for (iyzuv=nalts=ialts=ambig = 0, iqalt = -1; iyzuv<4; iyzuv++) {
	    int iquad;
	    switch (iyzuv) {
	    case 0 : iquad = ifl[*ihy]; break;
	    case 1 : iquad = ifl[*ihz]; break;
	    case 2 : iquad = ifl[*ihu]; break;
	    default:
	    case 3 : iquad = ifl[*ihv];
	    }
	    if      (iquad==-2) {            // Preexisting Ambiguity => give up
	      iqalt = -2; ambig = 1; break;
	    }
	    else if (iquad>=0) {         // Already member of an other quad...
	      if (iquad!=iqalt && iqalt>=0) {   // ...Twice ambiguous => give up
		// ...flag components of 2nd aternative quads as ambiguous (-2)
		for (int jyzuv = 0; jyzuv<4; jyzuv++) {
		  ifl[quads[iquad].jhs[jyzuv]] = -2;
		}
		// ...``delete'' the quads
		quads[iquad].jhs[0] = -1;
		ambig = 1; break;
	      }
	      iqalt = iquad;                     // ...Remember alternative quad
	      nalts++; ialts |= 1<<iyzuv;     // ...Remember alternative coord's
	    }
	  }
	  if (iqalt>=0) {                 // Single ambiguity...
	    if (TOpt::ReMode[17]) { // Try and raise it using GEM correlation
	      if (nalts==3 || nalts==1) {
		// Alternatives have 3 components in common
		int iambi;
		if  (nalts==3) ialts = ~ialts;
		for (iambi = 0; iambi<4; iambi++) { if ((1<<iambi)&ialts) break; }
		if  (nalts==1) iambi = iambi/2 ? 5-iambi : 1-iambi;
		CsCluster *cref, *cl1, *cl2;
		switch (iambi) {
		case 0: cref = hz.ptrCl; cl2 = hy.ptrCl; break;
		case 1: cref = hy.ptrCl; cl2 = hz.ptrCl; break;
		case 2: cref = hv.ptrCl; cl2 = hu.ptrCl; break;
		case 3: cref = hu.ptrCl; cl2 = hv.ptrCl; break;
		default: cref=cl2 = 0;
		}
		cl1 = vecHit[quads[iqalt].jhs[iambi]].ptrCl;
		double ref   = cref->getAllAnalogData()[2], dd;
		double delta1 = cl1->getAllAnalogData()[2];
		double delta2 = cl2->getAllAnalogData()[2];
		switch (iambi) {
		case 0:
		  dd = ayz0+ayz1*delta1+ayz2*delta1*delta1; delta1 = dd - ref;
		  dd = ayz0+ayz1*delta2+ayz2*delta2*delta2; delta2 = dd - ref;
		  break;
		case 1: dd = ayz0+ayz1*ref+ayz2*ref*ref;
		  delta1 -=dd; delta2 -= dd; break;
		case 2:
		  dd = auv0+auv1*delta1+auv2*delta1*delta1; delta1 = dd - ref;
		  dd = auv0+auv1*delta2+auv2*delta2*delta2; delta2 = dd - ref;
		  break;
		case 3: dd = auv0+auv1*ref+auv2*ref*ref;
		  delta1 -=dd; delta2 -= dd; break;
		}
		delta1 = fabs(delta1); delta2 = fabs(delta2);
		if ((fabs(dV)<fabs(quads[iqalt].dV)-resol/2 ||
		     delta2<delta1-delta) &&
		    delta2<2*delta) {
		  // Current quad is much better than old one...
		  // ...flag ambiguous components of old quad as unused (-1)
		  if  (nalts==1) ialts = ~ialts;
		  for (iyzuv= 0; iyzuv<4; iyzuv++) if ((1<<iyzuv)&ialts)
		    ifl[quads[iqalt].jhs[iyzuv]] = -1;
		  // ...replace old quad by current one
		  quads[iqalt].Y = (V*syu-U*syv)/svu;
		  quads[iqalt].V = (Z*svz-Y*svy)/szy;
		  quads[iqalt].jhs[0] = *ihy; quads[iqalt].jhs[1] = *ihz;
		  quads[iqalt].jhs[2] = *ihu; quads[iqalt].jhs[3] = *ihv;
		  quads[iqalt].dV = dV;
		  ifl[*ihy]=ifl[*ihz]=ifl[*ihu]=ifl[*ihv] = iqalt;
		  iqalt = -2;    // To prevent anything from being done infra
		}
		else if (fabs(quads[iqalt].dV)>fabs(dV)-resol/2 &&
			 (delta1>delta2-delta || 
			  delta1>2*delta)) ambig = 1;   // Ambiguity persists
		else  // Old quad is much better than current one: do nothing
		  iqalt = -2;
		}
	      else ambig = 1;  // Complicate case => don't try and raise ambiguity
	    }
	    else ambig = 1;    // Ain't using GEM correlation => ambiguity persists
	  }
	  if (iqalt>=0) {
	    // ...flag components of aternative quad as ambiguous (-2)
	    for (int jyzuv = 0; jyzuv<4; jyzuv++) {
	      ifl[quads[iqalt].jhs[jyzuv]] = -2;
	    }
	    // ...``delete'' the quad
	    quads[iqalt].jhs[0] = -1;
	  }
	  if (ambig) {            // Non raisable ambiguity ...
	    // ...flag all components of current quad as ambiguous (-2)
	    ifl[*ihy]=ifl[*ihz]=ifl[*ihu]=ifl[*ihv] = -2;
	  }
	  else if (iqalt==-1) {

	    // ******************** UNAMBIGUOUS QUAD ********************

	    quads[nquads].Y = (V*syu-U*syv)/svu;
	    quads[nquads].V = (Z*svz-Y*svy)/szy;
	    quads[nquads].jhs[0] = *ihy; quads[nquads].jhs[1] = *ihz;
	    quads[nquads].jhs[2] = *ihu; quads[nquads].jhs[3] = *ihv;
	    quads[nquads].dV = dV;
	    ifl[*ihy]=ifl[*ihz]=ifl[*ihu]=ifl[*ihv] = nquads;
	    nquads++;
	  }

	}   // end loop over v hits
      }   // end loop over u hits
    }   // end loop over z hits
  }   // end loop over y hits

  // ******************** ORDER AND PACK ARRAYS OF HITS ********************

  Quad *qps[nquads]; int mquads = 0;
  for (int i = 0; i<nquads; i++) {
    Quad *qp = &quads[i];
    int jhy = qp->jhs[0];
    if (jhy<0) continue;
    int jhz = qp->jhs[1], jhu = qp->jhs[2], jhv = qp->jhs[3];
    THit &hy = vecHit[jhy], &hz = vecHit[jhz];
    double anay = hy.ptrCl->getAllAnalogData()[2];
    double anaz = hz.ptrCl->getAllAnalogData()[2];
    double deltazy = anaz-ayz0-ayz1*anay-ayz2*anay*anay;
    if (fabs(deltazy)>2.5*delta) {                            // BAD AMP CORR...
      ifl[jhy]=ifl[jhz]=ifl[jhu]=ifl[jhv] = -1;    continue;  // ...=> DISCARD
    }
    THit &hu = vecHit[jhu], &hv = vecHit[jhv];
    double anau = hu.ptrCl->getAllAnalogData()[2];
    double anav = hv.ptrCl->getAllAnalogData()[2];
    double deltavu = anav-auv0-auv1*anau-auv2*anau*anau;
    if (fabs(deltavu)>2.5*delta) {                            // BAD AMP CORR...
      ifl[jhy]=ifl[jhz]=ifl[jhu]=ifl[jhv] = -1;    continue;  // ...=> DISCARD
    }
    double st, st2; int ih;
    double times[] = {hy.Time,hz.Time,hu.Time,hv.Time};
    double sigTs[] = {hy.SigT,hz.SigT,hu.SigT,hv.SigT};
    for (ih = 0, st=st2 = 0; ih<4; ih++) {
      if (sigTs[ih]>100 || sigTs[ih]<0) { st2 = 1.e6; break; }
      double t = times[ih]; st += t; st2 += t*t;
    } 
#ifdef Quadruples_HISTO
    if (TOpt::Hist[6]) tGEMt[id]->Fill(sqrt(st2-st*st/4));
#endif
    if ((st2-st*st/4)>24*24) { // 2-SIGMA(=12ns) CUT ON TIME SPREAD...
      // 2-sigma, i.e. very tight, because associated (X/Y, V/U) hits must be
      // correlated also in time.
      ifl[jhy]=ifl[jhz]=ifl[jhu]=ifl[jhv] = -1;    continue;  // ...=> DISCARD
    }
#ifdef Quadruples_HISTO
    if (TOpt::Hist[6]) {
      tGEMYZ[id]->Fill(deltazy); tGEMUV[id]->Fill(deltavu);
    }
#endif
    if (fabs(qp->dV)>py.Station->VRange/2) continue;
    qps[mquads++] = qp;
  }
  nquads = mquads;
#ifdef Quadruples_HISTO
  if (isMC && TOpt::Hist[6]) {
    // Log errors
    static int n_errors[nGEMs] = {0,0,0,0,0,0,0,0,0,0,0};
    for (int i = 0; i<mquads; i++) {
      Quad *qp = qps[i];
      THit &hy = vecHit[qp->jhs[0]], &hz = vecHit[qp->jhs[1]],
	&hu = vecHit[qp->jhs[2]], &hv = vecHit[qp->jhs[3]];
      int ikine = hy.IKine;
      if (ikine>=0 && hy.IKin2<0) {
	if (hz.IKine!=ikine && (hz.IKin2<0 || hz.IKin2!=ikine) ||
	    hu.IKine!=ikine && (hu.IKin2<0 || hu.IKin2!=ikine) || 
	    hv.IKine!=ikine && (hv.IKin2<0 || hv.IKin2!=ikine)) {
#  define Quadruples_DEBUG
#  ifdef Quadruples_DEBUG
	  printf("** Quadruples: Error\n");
	  printf
	    ("%d GEM %d %d %d %d %d Y %d %d %d Z %d %d %d U %d %d %d V %d %d %d\n",
	     n_errors[id]++,id,hy.IPlane,hz.IPlane,hu.IPlane,hv.IPlane,
	     hy.IHit,hy.IKine,hy.IKin2, hz.IHit,hz.IKine,hz.IKin2,
	     hu.IHit,hu.IKine,hu.IKin2, hv.IHit,hv.IKine,hv.IKin2);
	  double ay = hy.ptrCl->getAllAnalogData()[2];
	  double az = hz.ptrCl->getAllAnalogData()[2];
	  double au = hu.ptrCl->getAllAnalogData()[2];
	  double av = hv.ptrCl->getAllAnalogData()[2];
	  printf("dV %f  %f    dYZ %f dUV %f  %f\n",qp->dV,resol,
		 ay-az,au-av,delta);
#  endif
	  if (TOpt::Hist[6]) tGEMerr->Fill(id);
	}
	else if (TOpt::Hist[6]) tGEMok->Fill(id);  
      }
    }
  }
#endif
  // Order Y hits
  for (int i = 0; i<nquads; i++) {
    for (int j = i+1; j<nquads; j++) {
      Quad *exch = qps[j];
      if (exch->Y<qps[i]->Y) {
	qps[j] = qps[i]; qps[i] = exch;
      }
    }
  }
  // Copy Y hits and reference to output arrays
  for (int i = 0; i<nquads; i++) {
    yhit[i] = qps[i]->Y; yref[i] = qps[i]->jhs[0];
  }
  // Order Z hits
  for (int i = 0; i<nquads; i++) {
    for (int j = i+1; j<nquads; j++) {
      Quad *exch = qps[j];
      if (exch->V<qps[i]->V) {
	qps[j] = qps[i]; qps[i] = exch;
      }
    }
  }
  // Copy Z hits and reference to output arrays
  for (int i = 0; i<nquads; i++) {
    zhit[i] = qps[i]->V; zref[i] = qps[i]->jhs[1];
  }

  return nquads;

}
