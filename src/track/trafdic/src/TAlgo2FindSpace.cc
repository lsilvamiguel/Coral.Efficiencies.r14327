// $Id: TAlgo2FindSpace.cc 14094 2015-11-06 15:28:48Z lsilva $

/*!
  \file TAlgo2FindSpace.cc
  \brief Pattern Recognition of space tracks based on input proj. tracks

  \param \e igr   = Group# (a group being a reconstruction zone (if "igr<5") or a set of zones)
  \param \e X0    = Reference X abscissa

  \param \e npl   = # of detector planes
  \param \e idpl[ipl] = Index of plane "ipl" in TSetup::vDetect
  \param \e cosa[ipl], sina[ipl] = Orientation
  \param \e tols[ipl] = Tolerance on residuals
  \param \e fh  [ipl] = First hit index (in "Uhit")
  \param \e lh  [ipl] = After-the-last hit index ("Uhit"). When no hit: "fh[ipl]==lh[ipl]"

  \param \e nhits = Total # of hits
  \param \e Uhit[ih] = Hit coordinate in detector reference system (DRS)

  \param \e nproj = # of projections where proj. tracks are available
  \param \e Ntrk_prj[ip] = # of tracks in proj. "ip"
  \param \e cos_prj[ip], sin_prj[ip] = Orientation
  \param \e Y0_prj[ip][it] = Intercept of track "it" in proj. "ip"
  \param \e Yp_prj[ip][it] = Slope
  \param \e Y2_prj    [it] = Quadratic term

  \param \e mode  = Processing mode. <=0: Strict requirements, >0: Looser requirements

 Output:
  \param \e Ntk   = # of found tracks

  \param \e hits[] = Array of hit references, separated by "-1": [...,h0_T,h1_T,...,hn_T, -1 ,h0_T+1,...] for hits 0,...,n of successive tracks. Negative references = alternatives.

  \param \e Y0_trk[] = Track intercept at X0
  \param \e Yp_trk[] = Track slope

  \param \e ifl[]  = Used hits

*/

/*
 - The basic algorithm is extremely simple:
  I) Get the proj. tracks found by earlier "FindProjQ" step (cf. "_prj" arrays).
 II) Loop #1 over all combinations of any two of them:
   - Build space tracks and count #hits w/in a road opened about them.
   - Retain those w/ #hits > "iPRpar[2]", evaluate their quality function, Q.
III) Sort combinations according to Q.
 IV) Loop #2 over all retained combinations in descending Q order.
   - Count again #hits, keeping # of already used hits < "iPRpar[1]".
  V) Output "Ntk" space tracks (encoded in "hits") and used hits ("ifl").
     Return a <0 value in case of overflow of buffers.
 - There are a number of complications, though, associated w/:
  I) Low redundancy in VSAT => VSAT ENHANCEMENT.
 II) Ambiguities in drifts => DRIFT SUBDUING.
 - VSAT ENHANCEMENT brings in turn its own complications:
  I) The amount of enhancement depends upon data taking setup.
 II) No reliable Q of VSAT combinations (because have typically no more than 2
    space points, on trajectories that are paraxial => no good chi2) => Need
    retain many of them, if not all => Spiraling combinatorics, w/ many clones.
     => "FindS_DISCARD_SF_CLONES"
     Some more explanations in:
    "http://wwwcompass.cern.ch/twiki/bin/view/DataReconstruction/TrackingSpecial#Beam_Region".
 - Other complications:
  I) Iterative call by "TEv::PrePattern2" w/ varying processing "mode" argument.
 II) Handling of ambiguities in the hit<->track association => Remember
    alternative association in "hits" array (w/ a distinctive negative entry).
III) Track fitting, enabled depending upon zone ("igr") and iteration ("mode"):
    can be either linear or quadratic
 IV) Extra argument "Y2_prj" for quadratic term in Y proj.
     Etc...
*/

#include <iostream>
#include <cstdio>
#include <cmath>
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "TROOT.h"
#include "TMath.h"
#include "Coral.h"
#include "CsErrLog.h"
#include "CsHistograms.h"
#include "TSetup.h"
#include "TOpt.h"
#include "TAlgo2.h"
#include "TConstants.h"
#include "cfortran.h"

using namespace CLHEP;

PROTOCCALLSFSUB6(SORTZV,sortzv,PFLOAT,PINT,INT,INT,INT,INT)
#define SORTZV(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(SORTZV,sortzv,PFLOAT,PINT,INT,INT,INT,INT,A1,A2,A3,A4,A5,A6)

int TAlgo2::FindSpace(int igr, int mode, float X0,
		      int npl,int idpl[],float cosa[],float sina[],float tols[],
		      int fh[],int lh[],
		      int nhits,float Uhit[],
		      int nproj,int Ntk_prj[],float cos_prj[],float sin_prj[], 
		      float Y0_prj[13][TConstants_NTtrack_max],
		      float Yp_prj[13][TConstants_NTtrack_max], float *Y2_prj,
		      float *hinfos,
		      int &Ntk, int hits[], int ifl[])
{
  const TSetup& setup = TSetup::Ref(); int *ihinfos; ihinfos = (int*)hinfos;

  int zone; // ***** DETECTORS GROUP -> SPECTROMETER ZONE *****
  if      (igr<5)  zone = igr;
  else {
    switch (igr) {
      //case 5:  zone = -1; break; // Never occur in SPACE track finding
    case 8:  zone = 5; break;  // Zones 0x4 and 0x8 at once
    case 9:  zone = 0; break;
    default: zone = igr-6;
    }
  }

  unsigned int prj2Projs[nproj]; memset((void*)prj2Projs,0,sizeof(prj2Projs));
  unsigned int zProjs = 0, uv10Projs = 0, UVProjs = 0;// ***** PROJECTIONS *****
  for (int iproj = 0; iproj<(int)setup.vProj().size(); iproj++) {
    unsigned int proj = 1<<iproj; double alpha = setup.vProj()[iproj]/10.;
    if      (fabs(alpha-90)<2)                    zProjs |= proj;
    else if (fabs(alpha-10)<2)                 uv10Projs |= proj;
    else if (18<fabs(alpha) && fabs(alpha)<47)   UVProjs |= proj;
    for (int prj = 0; prj<nproj; prj++) {
      double beta = atan2(sin_prj[prj],cos_prj[prj])*180/M_PI;
      if (fabs(beta-alpha)<2) {
	if (prj2Projs[prj])
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "Group %d: arg proj. %d (%.1f deg.) matches both 0x%x and 0x%x",
			igr,prj,beta,prj2Projs[prj],proj);
	prj2Projs[prj] = proj;
      }
    }
  }
  /*
  for (int prj = 0; prj<nproj; prj++) if (!prj2Projs[prj])
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "Group %d: arg proj. %d does not match any of vecProj",igr,prj);
  */
  const int maxcomb = 20000; // maximum number of combinations
  static float Ysl [maxcomb], Y0 [maxcomb];
  static float Y2  [maxcomb], xYLow[maxcomb], xYUp[maxcomb];
  static float Zsl [maxcomb], Z0 [maxcomb];
  static float Qt  [maxcomb];
  static int  ind  [maxcomb];
  int ipl;

  // ***** SCIFI SPECIFIC
  static int  sfs_only[maxcomb];
  const int nCSfsMx = 4096; // "CSfs" == scifis combinations
  static float timesSfs[nCSfsMx]; static int nScSfs[nCSfsMx], usedSfs[nCSfsMx];
  int nCSfs = 0, cSfsOver = 0;

  // ********** REQUIREMENTS on PR PARAMETERS **********

  int mode_par = mode>0 ? 3 : 0;
  int nHitsMn =                               // ***** MINIMUM # of HITS *****
    TOpt::iPRpar[igr*10+2 + mode_par];   // Depends upon arg "mode"
  if (mode<0) nHitsMn += 1;
  float nHitsMn_Mx = TOpt::iPRpar[igr*10+2]; // Max of min #hits over iterations

  // Special option for DY beamT: dCut[37];
  int beamTinDY = igr==4 && TOpt::iCut[37] ? TOpt::iCut[37] : 0; 

  int scifiType;             // ***** SCIFI/VSAT ENHANCEMENT *****
  bool vSATEnhanced = // Condition vSATEnhanced by zone#, "cf. iCut[16]"
    TOpt::iCut[16]&1<<zone;
  if (igr<2) {// First two groups: the enhancement is delayed so as to limit the
    // combinatorics: very good candidate tracks (fulfilling PR requirements w/o
    // any help from enhancement, and lying typically not in the VSAT region)
    // are found in the first iteration(s), and their hits are removed from the
    // problem. Cf. "iCut[28]".
    // => Condition enhancement by iteration# (retrieved from "mode" arg.)
    bool ok = 1<<(mode>=0 ? mode : 0)&TOpt::iCut[28];
    if (igr==1 && TOpt::iCut[31]) // Special case of 2010-like VSAT setup
      ok |= mode<=0;
    vSATEnhanced &= ok;
  }
  // Disable for low P forward PR: scifis, viz. FI03, plays the role of a SAT (
  // as opposed to VSAT)
  // But keep it enabled for other 2-zone Z-proj. (even "igr==11"!?) 
  vSATEnhanced &= igr!=9;
  // If "vSATEnhanced", flag relevant detector types, here scifi type, w/ 0x200
  // Same if special DY option.
  if (vSATEnhanced || beamTinDY) scifiType = 22|0x200;
  else                           scifiType = -1;

  float incrVSAT = 0; if (vSATEnhanced) switch (zone) {
  case 0:
    if (TOpt::iCut[16]==0x1) // This is made to represent the case of field off
      // single zone (although this definition of "iCut[16]" introduces some
      // ambiguity...): 6 0x1 + 5 (or even, 7 in 2006) 0x2 scifis
      /* */                 incrVSAT = nHitsMn/8.-1;
    else {
      if (TOpt::ReMode[48]) incrVSAT = nHitsMn/4.-1; // 2*MPM's+2*FI's
      else                  incrVSAT = nHitsMn/5.-1;
    }
    break;
  case 1:
    // Case of zone 0x2:
    // - In 2006: #scifis increases from 5 to 7, which could allow us to set
    //  "incrVSAT = nHitsMn/5.-1". Doing so (tested on D*.mc.2006.02.11-13)
    //  does not affect significantly reco perfs though.
    // - In 2010: Still 7 scifis, but also 2 GPs, providing 4 hits. This has 2
    //  effects:
    //   - Easier reconstruction: this time, we can go for "nHitsMn/8.-1".
    //   - Much higher # of combinations in FindProjQ (where one still only
    //    requires 2 hits per proj.) => need for withdrawing as many VSAT hits
    //    as possible before reaching the iteration #'s (typically, mode&0x3,
    //    cf. iCut[28]) where scifi-enhancing takes place. And since the max.
    //    # of VSAT hits is still low (=11), one needs also an enhancement in
    //    1st (mode<=0) iterations.
    //   In order to implement both, special option iCut[31] was introduced.
    //   Note that the behaviour of the routine upon option iCut[31] is not
    //  actually driven by options (iCut[16/28]): it's rather expecting a
    //  specific setting of options in order to work properly. This should  be
    //  changed in the future, but would then require a complete reshuffling.
    // - In 2009 DVCS: intermediate case: 5 scifis + 2 GPs.
    if      (TOpt::iCut[31]>=5) {              // Case 2010
      if (mode<=0) incrVSAT = nHitsMn/9.-1;
      else         incrVSAT = nHitsMn/8.-1;
    }
    else if (TOpt::iCut[31]) {                 // Case 2009 DVCS
      if (mode<=0) incrVSAT = nHitsMn/7.-1;
      else         incrVSAT = nHitsMn/6.-1;
    }
    else           incrVSAT = nHitsMn/4.-1;
    break;
  case 2: incrVSAT = nHitsMn/4.-1; break;
  default: incrVSAT = 0;
  }
  incrVSAT += .001;

  // Safeguard against VSAT too much enhanced
  int nSafe = TOpt::ReMode[48] && zone==0 ? 6 : 5;

  // Stricter requirement for ``long'' track in zone 0x1
  int nHitsMn_FI05 = nHitsMn+TOpt::iCut[5];
  int nHitsMn_GM04 = nHitsMn+(TOpt::iCut[5]>0 ? 1 : 0);
  double xGM04 = 1400-X0, xFI05 = 580-X0, xRICH = zone==1 ? 750-X0 : -1;
  if (igr==4) { // Special case of zone 0x10 (beam): loosen requirement on #hits
    // The rationale is to cope so w/ the large diff in resolution between
    // scifis and Sis, which has at times the routes issued by the proj.
    // search defined by scifis and their combination failing to pickup
    // all of genuine Si hits: the initial search, w/ loosen requirement, will
    // provide, via fitting, a better route, to which the original, restored
    // infra, will be applied.
    nHitsMn -= 1;
  }
  //                                             ***** MINIMUM CHI2 *****
  double chi2_mn = TOpt::dCut[55+zone];
  bool doLinFit  = (TOpt::iCut[18]&1<<zone) &&// Linear fit...
    (mode<=0 ||                  // ...only in earlier iteration(s)...
     (1<<zone&0x18));            // ...or in field-free zones
  bool doQuadFit = (TOpt::iCut[19]&1<<zone) &&// Quad fit...
    mode>1;                                   // ...only in later iterations
  bool doSIFit = // Special case of SI-tracks in zone 0x1...
    // ...Routes about SIs are enlarged (beyond SI's resolution) in order to
    // make for their large diff in resolution w.r.t. other detectors. This has
    // the side-effect of qualifying ghosts which take advantage of the fact to
    // pick up SI hits. => Let's fit those tracks, fitting being here based on
    // the detector's resolution (N.B.: not necessarily the intrinsic one, cf.
    // "dCut[84]"), cf. "wghts" array infra.
    zone==0 && TOpt::ReMode[43] && mode<=0;
  bool doVSATFit; if (!doLinFit) {            // Linear fit in VSAT only
    doVSATFit = (TOpt::iCut[32]&1<<zone) &&
      mode<=0;                   // ...only in earlier iteration(s)
  }
  else doVSATFit = false;
  int nFit;
  if      (doQuadFit)                        nFit = 5;
  else if (doLinFit || doSIFit || doVSATFit) nFit = 4;
  else                                       nFit = 0;
  float  aij[5*npl], wghts[npl];
  double Ajk[15], AjV[5]; static double Vi2;
  float chi2_scale = TOpt::dPRpar[igr*10+0]*TOpt::dPRpar[igr*10+0] *
    /* Empirical factor designed to have greater sensitivity to small variations
       of chi2 for small chi2 */ (TOpt::dCut[62] ? TOpt::dCut[62] : 1);

#define FOUND_vs_EXPECTED 2
  //                                             ***** DETECTOR EFFICIENCY *****
  int nexpected;
  static int nexp0, nexp1; // "static" in order to prevent compiler from warning
  float tracking_eff;
  if (zone==5) { // Special case of paraxial PR: average over zones 0xc
    int i; for (i = 2, tracking_eff = 0; i<4; i++)
      tracking_eff += TOpt::dCut[50+i]    ? TOpt::dCut[50+i]    : .5;
    tracking_eff /= 2;
  }
  else tracking_eff = TOpt::dCut[50+zone] ? TOpt::dCut[50+zone] : .5;

  //                                             ***** SCIFI TIME CUT *****
  float scfDt_mx;
  if (zone==5)  // Special case of paraxial PR...
    scfDt_mx = TOpt::dCut[21+2]; // ...no scifi in zone 0x8: use zone 0x4 alone
  else scfDt_mx = TOpt::dCut[21+zone];
  //                                             ***** DRIFTS SUBDUED? *****
  bool driftsSubdued = 1<<zone&TOpt::iCut[17] ? true : false;

  float xpl [npl];
  if (igr!=9) for (int i = 0; i<nhits; i++) ifl[i]=0;   // FLAGS for HITS

  int typl[npl], prpl[npl];  // ********** DETECTORS TYPES, PROJS **********
#define  FindS_DISCARD_SF_CLONES  // Note: it's in fact ``discard VSAT clones''
#ifdef FindS_DISCARD_SF_CLONES
#  define N_PSFS_MX 16
  int nHsSfs = 0;
  static unsigned int hitsPsSfs[nCSfsMx+1][N_PSFS_MX];
  static int nFoundSfs[nCSfsMx], iQualSfs[nCSfsMx];
  int fh_sfs[npl];  // Array of first scifi hits (analogue to "fh" but specific)
#endif
#ifdef PR_ENHANCE_CENTRAL_GEMs
  float yCentralGEM[npl], zCentralGEM[npl];
#endif
#ifdef PR_GP_WHICHSIDE_INFO
  int ipl28 = -1; double offsets28[4];
  // Upon option or if a single-zone job, disable the making use of the which-
  // side info in stripped GPs, cf. "TEv::PrePattern2".
  bool useGPWhichSide = setup.vIplFirst().size()!=1 && !TOpt::ReMode[45];
#endif
  float wght_scale = 1; for (ipl = 0; ipl<npl; ipl++) {
    const TDetect &d = setup.vDetect(idpl[ipl]);
    if (d.IType!=21) {
      wght_scale = tols[ipl]/setup.vDetect(idpl[ipl]).Resol; break;
    }
  }
  //    ***** TOLERANCES FOR M-PIXELS (as opposed to P-pixels) *****
  // - The uncertainty along v is different from that along u and variable.
  // On 06/21/2012 10:35 AM, THIBAUD Florian wrote:
  // La partie pixellisée des MPs,v2012 est composée d'une partie centrale de
  // 2.56cm (u) x 2.5cm (v) où les pixels ont un pitch en v de 2.5mm, et d'une
  // partie périphérique de 5.12cm x 5cm avec un pitch en v de 6.25mm.
  // - The best way to take this into account would be to derive the tolerance
  //  from the cluster cov. matrix. But in our present scheme, we do not have
  //  this info at hand => let's use built-in values.
  // - Also, in what follows, when it comes to opt for the outer, more favorable
  //  tolerance: we apply a margin of 1./1.6 mm taking into account tracking
  //  uncertainty and misalignment. We also put on a same footing the u and v
  //  axes, hence allowing us to handle all M-pixels in the same way, whatever
  //  be their axes.
  // - These are definitely very rough approximations. And this whole processing
  //  should be revisited at some point, if we want to be more efficient.
  // (N.B.: We ignore P-pixelated and M-pixelated v2011 MPs: not so useful,
  // since these detectors were actually systematically DetNameOff'd.)
  // - Option "dCut[98]" allows to relax the tolerance. This is useful for
  //  processing data where MPs would not have yet ben aligned along their v
  //  dimension. (Let's recall that this v-alignment cannot be obtained w/ the
  //  standard alignment procedure, as opposed to that of pixelated GEM planes:
  //  given that i) the v-granularity of the pixels is very coarse and ii) the
  //  pixelated piece has no fixed position w.r.t. any v-measuring strip piece.)
  float tolMPvI = .25/sqrt(12)*wght_scale, tolMPvO = .625/sqrt(12)*wght_scale;
  if (TOpt::dCut[98]) {
    tolMPvI += TOpt::dCut[98]; tolMPvO += TOpt::dCut[98];
  }
  float minSIRes = TOpt::dCut[84]>.0015 ? TOpt::dCut[84] : .0015;
#define FindS_WO_1ST_SCIFI
//  This patch is to work around a deficiency of the "FindSpace" algorithm in
// dealing w/ time info: the time consistency of a track in progress is checked,
// for each new candidate hit, w.r.t. already associated hits. Therefore, if the
// first associated hit turns out to have a bad timing, all successive well
// timed ones are rejected. In the present (as of 2010/04) implementation of the
// algorithm, only scifis are considered for time info. This patch is meant for
// such a scifi-restricted implementation, and further restricted to cases where
// there are as many as 4 scifis, w/ the rationale being: if a track turns out
// to have few (in fact ==1) scifi hits, the 1st scifi hit is suspicious and the
// search is tried again, based on the same combination of proj. tracks, but
// skipping the 1st scifi plane (in fact, skipping up to 2 scifi planes up to
// the one associated in the 1st place and considered suspicious). The patch
// only works (because of its actual implementation, which is meant to be as
// simple (and hence fast) as possible) in zone 0x10, w/ a setup w/ FI01 at the
// very beginning of the zone.
//  Note: scifis detectors seem to have a sizable rate of firing where the
// response is correct in space, but bad in time, e.g. 07W30/cdr27029-59891
// event #121671115. I haven't ever quantified this, though.
#ifdef FindS_WO_1ST_SCIFI
  int nScifiDets = 0;  // To count the # of scifi planes, and check it's >=4.
#endif
  for (ipl = 0; ipl<npl; ipl++) {
    prpl[ipl] = setup.vPlane(idpl[ipl]).IProj;
    const TDetect &d = setup.vDetect(idpl[ipl]);
    double xi = d.X(0)-X0; xpl[ipl] = xi; int type = d.IType; typl[ipl]= type;
    if (11<=type && type<=18) {
      if (driftsSubdued)        typl[ipl] |= 0x100; // Drifts: subdued? Tag them
#ifdef PR_ST_WHICHSIDE_INFO                        // Use which-side info?...
      if (11<=type && type<=12) typl[ipl] |= 0x400;   // ...Tag central straws
#endif
    }
    if ((vSATEnhanced || 1<igr) && igr!=9) {      // VSAT enhanced?...
      if (type==22 || type==28 || type==29 ||  // ...Tag scifis and pixelGEMs...
	  // ...and pixelated MPs upon option (the option being expected to be
	  // enabled in the case of DVCS 2012 setup where few MPs are at hand).
	  type==32 && TOpt::ReMode[48]) {
	typl[ipl] |= 0x200;
#ifdef FindS_WO_1ST_SCIFI
	if (type==22) nScifiDets++;
#endif
      }
#ifdef PR_ENHANCE_CENTRAL_GEMs               // ...and central GEMs
      // Enhance GEM w/ central zone activated the same way as scifis are.
      // This is conditioned by "PR_ENHANCE_CENTRAL_GEMs", which is expected
      // to be defined by "make ./src/track/lattice", cf. ../../GNUMakefile...
      if (typl[ipl]==26 && (igr==1 || igr==2) &&//...and limited to zones 0x6...
	  !d.DZSize(0)) {             // ...and requiring DZ is indeed activated
	typl[ipl] |= 0x200;
	// Remember also the offsets.
	yCentralGEM[ipl] = d.X(1); zCentralGEM[ipl] = d.X(2);
      }
#endif
    }
    if (type==29 || type==32) typl[ipl] |= 0x800;  // Tag pixel info
#ifdef PR_GP_WHICHSIDE_INFO                        // Use which-side info?...
    if (type==28) {      // Which-side info in stripped pixelGEMS: store offsets
      if (ipl28<0) ipl28 = ipl;
      if (ipl-ipl28<4) {
	// The offset must be that of the associated plane (Y<->X,U<->V). This
	// is taken care of later on. It must be counted along the DRS v axis
	// which is u axis +pi/2 => Therefore for 'X' and 'V' (later to be
	// stored in 'Y' and 'U'), take opposite.
	offsets28[ipl-ipl28] = d.Uorig+d.Pitch*(d.Nwires/2-.5);
	if (d.Name[4]=='X' || d.Name[4]=='V') offsets28[ipl-ipl28] *= -1;
      }
      else if (useGPWhichSide) {// No more than 4 stripped pixelGEMs per zone...
	// ...No else than 4 (at most) in a row, in fact. In cases which do not
	// fit, e.g. single-zone magnets-OFF job, the which-side feature must
	// have been disbaled in "TEv::PrePattern2". Else: fatal.
	const TDetect &first28 = setup.vDetect(idpl[ipl28]);
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
          "PixelGEMs (stripped) in zone 0x%x: #%d(\"%s\") > #%d(\"%s\") + 4",
	  1<<zone,idpl[ipl],d.Name.c_str(),idpl[ipl28],first28.Name.c_str());
      }
    }
#endif

    // ***** "wghts": ALTERNATIVE ARRAY of TOLERANCES *****
    // Purpose: Restore Si resolution in fit and in 2nd loop, i.e.  once the bad
    // definition of the routes by scifis has been cured.
    if (nFit) {
      float wght; // Special cases of Si in zones 0x11, GEMs in zone 0x8 (GM11):
      // Tolerance has been enlarged in PrePattern2 w.r.t. resolution.
      // Restore it. But not entirely in order to make for possible
      // misalignment of Sis/GEM11 (which are difficult to align w.r.t.
      // the rest of the detection).
      float c = cosa[ipl], s = sina[ipl];
      if (igr==3 && TOpt::ReMode[34])
	wght = (typl[ipl]&0xff)==26 ? .0200*wght_scale : tols[ipl];
      else if (igr==4 || igr==0)
	wght = (typl[ipl]&0xff)==21 ? minSIRes*wght_scale : tols[ipl];
      else
	wght = tols[ipl];
      wghts[ipl] = wght;
      int i5 = 5*ipl;
      double sixij = sina[ipl]/wght; aij[i5++] = sixij;
      sixij *= xi; aij[i5++] = sixij;
      double cixij = cosa[ipl]/wght; aij[i5++] = cixij;
      cixij *= xi; aij[i5++] = cixij;
      cixij *= xi; aij[i5]   = cixij;
    }
    else wghts[ipl] = tols[ipl];
#ifdef FindS_DISCARD_SF_CLONES
    if (typl[ipl]&0x200 &&
	(igr!=4 ||// BeamTelescope: no vSAT enhancement => no clones expected...
	 beamTinDY)) { // ...except in the DY setup
      fh_sfs[ipl] = nHsSfs; nHsSfs += lh[ipl]-fh[ipl];
    }
#endif
  } // End loop on planes
#ifdef PR_GP_WHICHSIDE_INFO
  // Preparing for which-side info of stripped pixelGEMs: exchange "offsets28"
  // of associated planes, which are here successive planes.
  for (int i = 0; i<4; i += 2) {
    double o = offsets28[i]; offsets28[i] = offsets28[i+1]; offsets28[i+1] = o;
  }
#endif
#ifdef FindS_DISCARD_SF_CLONES
  int nPsSfs = nHsSfs%32 ? nHsSfs/32+1 : nHsSfs/32;
  nPsSfs = nPsSfs<N_PSFS_MX ? nPsSfs : 0;
#endif
#define TARGET_POINTING_SPACE
#ifdef TARGET_POINTING_SPACE
  float trgtPtngY = zone==0 ? TOpt::dCut[64] : 0;
  float trgtPtngZ = zone<4  ? TOpt::dCut[64] : 0;
#endif

  int nFound, nused;
  float yp, y0, zp, z0, qual, yy, zz;
  int hitsList[npl], nPats = npl/32+1, plPats[nPats];

  //#define DEBUG_SCIFIS_TIMING
#ifdef DEBUG_SCIFIS_TIMING
  static int idebug = 0;
#endif
  int ip1, ncomb;
  //                         Special cases of paraxial and low P forward PRs...
  int ip1Mx = igr>4 ? 0 : nproj;  // ...require 1st proj. (=Zproj)
  int ip2Mx = nproj;
  for (ip1=ncomb = 0; ip1<=ip1Mx; ip1++) { // ***** LOOP ON PROJ'S (1) *****
    int nt1 = Ntk_prj[ip1]; float c1 = cos_prj[ip1], s1 = sin_prj[ip1];
#if defined PR_LOW_P_FORWARD || defined PR_2ZONE_ZPROJ
    // Special 2-zone proj. search in Z
    int ip2_2Z; if (igr==1 && mode==0 &&
		    (TOpt::iPRpar[100] || TOpt::iPRpar[110])) {
      ip2_2Z = TOpt::iPRpar[100]*TOpt::iPRpar[110] ? ip2Mx-2 : ip2Mx-1;
    }
    else ip2_2Z = ip2Mx /* I.e. a value of "ip2" out of reach */;
#endif
    for (int ip2 = ip1+1; ip2<ip2Mx; ip2++) { // ***** LOOP ON PROJ'S (2) *****
      int nt2 = Ntk_prj[ip2]; float c2 = cos_prj[ip2], s2 = sin_prj[ip2];
      float det = c1*s2-c2*s1;
#ifdef  PR_2ZONE_ZPROJ
      if (ip2>=ip2_2Z) {  // Special case of Z PR spanning 2 zones
	if (fabs(c1)<.02) continue;    // Exclude Z*Z combination
      }
#endif
#ifdef SKIP_XWCOMB
      if (igr<2 && mode<=0 && fabs(det)<.3) continue;
#endif
      if (fabs(det)<.12) continue; // Forbid SI combinations  (0,5) or (90,85).
      unsigned int xProjs = ~(prj2Projs[ip1]|prj2Projs[ip2]);

      for (int it1 = 0; it1<nt1; it1++) {// ***** LOOP ON PROJ 1 TRACKS *****
	float y1 = Y0_prj[ip1][it1], yp1 = Yp_prj[ip1][it1];
#ifdef QUADRATIC_PROJ
	float y21 = Y2_prj && fabs(s1)<.05 ? Y2_prj[it1] : 0;
#endif

	for (int it2 = 0; it2<nt2; it2++) { // ***** LOOP ON PROJ 2 TRACKS *****
	  float y2 = Y0_prj[ip2][it2], yp2 = Yp_prj[ip2][it2];
#ifdef QUADRATIC_PROJ
	  float y22 = Y2_prj && fabs(s2)<.05 ? Y2_prj[it2] : 0;
#endif

	  //                  ***** PIVOT LINE *****
	  y0 =  (y1*s2-s1*y2)/det; yp =  (yp1*s2-s1*yp2)/det;
	  z0 = -(y1*c2-c1*y2)/det; zp = -(yp1*c2-c1*yp2)/det;
#ifdef QUADRATIC_PROJ
	  float yq = (y21*s2-s1*y22)/det;
#endif

#ifdef FindS_WO_1ST_SCIFI
	  int iter = 0, ipl1stScifi = -1, nIters = 1; while (iter<nIters) {
	    // In case the timing of the 1st scifi hit is found suspicious (then
	    // "nIters" will have been reset =2, and "ipl1stScifi" =0/1), let's
	    // iterate the hit search, w/ the same combination of proj. tracks.
	    int ipl0 = iter ? ipl1stScifi+1 : 0; iter++;
#else
	    int ipl0 = 0;
#endif
	  nFound = 0; nexpected = 0;

	  //           ***** SCIFI (and VSAT) SPECIFIC *****
	  int nScifis = 0, kScOff = -1; static float st, st2, t_mn, t_mx;
#define SCIFIS_t_PROB
#ifdef SCIFIS_t_PROB
	  static float stProb;
#else
	  static float devt;
#endif
	  int nVSATs;

	  for (ipl = ipl0, memset(plPats,0,nPats*sizeof(int)), nVSATs = 0;
	       ipl<npl; ipl++) {
	    //          ********** LOOP ON PLANES **********
	    // In 2nd iteration, one skips 1st plane. This special feature is
	    // designed to rescue a genuine track which first scifi hit happens
	    // to have a bad timing (which has the consequence that other
	    // scifi hits are rejected because of bad time consistency). It is
	    // expected to be enabled only in case where zone==4 (one restricts
	    // oneself to zone==4 for the time being, in order to stay away from
	    // combinatorics) and zone==4 contains many (>=4) scifi planes and
	    // 1st iteration turned out a combination w/ few (==1) scifi hit.

	    // ***** COUNT HITS ON ALL PLANES IN DEFINED DIRECTION *****

	    yy = y0+yp*xpl[ipl]; zz = z0+zp*xpl[ipl];
#ifdef QUADRATIC_PROJ
	    yy += yq*xpl[ipl]*xpl[ipl];
#endif

	    if (!setup.vDetect(idpl[ipl]).InActive(yy,zz))
	      continue;                      // ***** IN ACTIVE AREA *****
	    nexpected++;

	    if (fh[ipl]==lh[ipl]) continue;
	    float u = yy*cosa[ipl] + zz*sina[ipl];
	    float tol = tols[ipl];
	    if (u<Uhit[fh[ipl]]  -tol ||     // ***** HITS w/in ROUTE *****
		u>Uhit[lh[ipl]-1]+tol) continue;

	    //                               // ***** FIND CLOSEST HIT *****
	    static int k0; int  kl = fh[ipl], kr = lh[ipl]-1; while (kr-kl>1) {
	      k0 = (kl+kr)>>1;
	      if (u<=Uhit[k0]) kr = k0;
	      else             kl = k0;
	    }
	    if (Uhit[kr]-u<u-Uhit[kl]) k0 = kr;
	    else                       k0 = kl;
	    if (fabs(u-Uhit[k0])>tol) continue;

	    int tysp = typl[ipl]; if (tysp&0xe00) { // ***** SPECIAL TYPES *****
	      int type = tysp&0xff;
	      if (tysp==scifiType) {                     // ***** SCIFI SPECIFIC
		// This is done here and not in the reevaluation step because
		// scifi timing may redefine the track.
		float scft = hinfos[k0];
		if (nScifis==0) {
		  nScifis=1; st = scft; st2 = scft*scft; t_mn=t_mx = scft;
		}
		else {
		  if (scft-t_mn>scfDt_mx || t_mx-scft>scfDt_mx) {
		    // Currently selected hit is out of time => try other one?
#ifdef DEBUG_SCIFIS_TIMING
		    if (idebug)
		      printf("%3d %3d %7.3f  %3d %3d\n",ipl,k0,scft,kr,kl);
#endif
		    int k0p, pm; for (pm = 0, k0p = -1; pm<2; pm++) {
		      int kp; if (pm) { kp = k0-1; if (kp<fh[ipl])  continue; }
		      else            { kp = k0+1; if (kp>=lh[ipl]) continue; }
		      if (fabs(u-Uhit[kp])<tol) {
			double scftp = hinfos[kp];
			if (scftp-t_mn<scfDt_mx && t_mx-scftp<scfDt_mx) {
			  if (k0p<0 || fabs(u-Uhit[kp])<fabs(u-Uhit[k0p])) {
			    k0p = kp; scft = scftp;
			  }
			}
		      }
		    }
		    if (k0p>=0) {
		      k0 = k0p;
		      nScifis++; st += scft; st2 += scft*scft;
		      if (scft<t_mn) t_mn = scft; if (scft>t_mx) t_mx = scft;
		    }
		    else if (beamTinDY>1 && kScOff<0)
		      // Special case of beamT in DY:
		      // - BeamT is very crowded => can be that a true, in-time
		      //  hit, is ovverriden by a pileup.
		      // - Redundancy is low => To tell true track from ghost,
		      //  we need very demanding PR requirement.
		      // => Let's allow for one off-time hit.
		      //    But let's not update track timing w/ it.
		      kScOff = k0;
		    else continue;
		  }
		  else {
		    nScifis++; st += scft; st2 += scft*scft;
		    if (scft<t_mn) t_mn = scft; if (scft>t_mx) t_mx = scft;
		  }
		}
#ifdef DEBUG_SCIFIS_TIMING
		if (idebug) {
		  if (nScifis==0) printf("-------- %d\n",ncomb);
		  printf("%3d %3d %7.3f  %3d %3d\n",ipl,k0,scft,kr,kl);
		}
#endif
	      } // End of scifi block
#ifdef PR_ENHANCE_CENTRAL_GEMs
	      else if (tysp==(26|0x200) &&
		       fabs(yy-yCentralGEM[ipl])<4.0 &&
		       fabs(zz-zCentralGEM[ipl])<4.0) nVSATs++;
#endif
	      else if (type==28) {            // ***** PixelGEM SPECIFIC: STRIPS
#ifdef PR_GP_WHICHSIDE_INFO
		int whichSide = ihinfos[k0]>>16; if (whichSide) {
		  float v = zz*cosa[ipl] - yy*sina[ipl] - offsets28[ipl-ipl28];
		  if (whichSide*v<-.5) { // Cut on the distance to the equator
		    // Currently selected hit !compatible w/ which-side info =>
		    // try other one (Note: we only consider here the 2 closest
		    // neighbours, the closest to the ``left'' and the closest
		    // to the ``right'', which is not optimum...)
		    if (kr==kl) continue;
		    if (k0==kr) k0 = kl;
		    else        k0 = kr;
		    whichSide = ihinfos[k0]>>16;
		    if (whichSide && whichSide*v<0) continue;
		    if (fabs(u-Uhit[k0])>tol) continue;
		  }
		}
#endif
		if (tysp==(28|0x200)) nVSATs++;
	      }
#ifdef PR_GP_WHICHSIDE_INFO
	      else if (type==31) {             // ***** PixelMM SPECIFIC: STRIPS
		int whichSide = ihinfos[k0]>>16; if (whichSide) {
		  // For pixelMMs, which are expected to be in zone 0x1, we
		  // assume the offset = 0. Otherwise, same processing as supra.
		  float v = zz*cosa[ipl] - yy*sina[ipl];
		  if (whichSide*v<-.5) { // Cut on the distance to the equator
		    if (kr==kl) continue;
		    if (k0==kr) k0 = kl;
		    else        k0 = kr;
		    whichSide = ihinfos[k0]>>16;
		    if (whichSide && whichSide*v<0) continue;
		    if (fabs(u-Uhit[k0])>tol) continue;
		  }
		}
	      }
#endif
	      else if (tysp&0x800) {       // ***** PixelGEM/MM SPECIFIC: PIXELS
		float v = zz*cosa[ipl] - yy*sina[ipl];
		float tolv; if (type==29) tolv = tol;
		else { // Tolerance for M-pixels: special, variable v-precision
		  // Let's apply the outer tolerance w/ a somewhat arbitrary
		  // margin, cf. explanations supra.
		  if (fabs(yy)>2.4 || fabs(zz)>2.4) tolv = tolMPvO;
		  else                              tolv = tolMPvI;
		}
		if (fabs(v-hinfos[k0])>tolv) {
		  // Currently selected hit not compatible w/ V-wise info =>
		  // try other one (Note: we only consider here the 2 closest
		  // neighbours, the closest to the ``left'' and the closest to
		  // the ``right'', which is not optimum...)
		  if (kr==kl) continue;
		  if (k0==kr) k0 = kl;
		  else        k0 = kr;
		  if (fabs(v-hinfos[k0])>tolv) continue;
		  if (fabs(u-Uhit[k0])>tol) continue;
		}
		if (type==29) { // GEM pixels provide 2D info. MM ones, not.
		  if (!nFound) nexp0 = nexpected; nexpected++;
		  nFound++;
		}
		if (tysp&0x200) { nVSATs++; if (type==29) nVSATs++; }
	      }
#ifdef PR_ST_WHICHSIDE_INFO
	      else if (tysp&0x400) {                    // ***** WHICH-SIDE INFO
		float whichSide = hinfos[k0]; if (whichSide) {
		  float v = zz*cosa[ipl] - yy*sina[ipl];
		  if (whichSide*v<0) {
		    // Currently selected hit !compatible w/ which-side info =>
		    // try other one (Note: we only consider here the 2 closest
		    // neighbours, the closest to the ``left'' and the closest
		    // to the ``right'', which is not optimum...)
		    if (kr==kl) continue;
		    if (k0==kr) k0 = kl;
		    else        k0 = kr;
		    whichSide = hinfos[k0]; if (whichSide && whichSide*v<0)
					      continue;
		    if (fabs(u-Uhit[k0])>tol) continue;
		  }
		}
	      }
#endif
	    } // End of special types
	    if (nFound) nexp1 = nexpected;
	    else        nexp0 = nexpected;

	    nFound++; plPats[ipl/32] |= 1<<(ipl%32); hitsList[ipl] = k0;

	  }  // End of loop over planes
	  // ********** CUT on NUMBER of HITS **********
	  if (nFound+(nVSATs+nScifis)*incrVSAT<nHitsMn) continue;
	  if (nFound<=nSafe && nScifis+nVSATs<4) continue;
	  // In 2012 DVCS: reject combinations w/ nFound = 5 = 2*MPMs+1*FI

#ifdef SCIFIS_t_PROB         // ***** CUT ON TIME CHI2 PROBA *****
	  int disWorstTime;
	  static float tCut;  // "static" to prevent compiler from warning
	  if (nScifis>=2 && // So that, in any case, infra is meaningful...
	      (igr==4 ||
	       zone==0 && nScifis+nVSATs>nFound/2 ||
	       zone>=1 && nScifis+nVSATs>=nFound/2)) { // E.g. 4*FI56+4*GM8
	    // ***** SciFiSi AND SCIFI-(ALMOST)ONLY TRACKS ***** Why here?
	    stProb = TMath::Prob(st2-st*st/nScifis,nScifis);
	    if (stProb<2.5e-2) { // ~=Prob(5,1) => 5 sigma if 1ns resol...
	      //                    ...Prepare for removing worst timed hit
	      float meanTime = st/nScifis; st=st2 = 0;
	      if (t_mx-meanTime>meanTime-t_mn) {
		disWorstTime =  1; tCut = t_mx-.01;
	      }
	      else {
		disWorstTime = -1; tCut = t_mn+.01;
	      }
	      t_mx = -100; t_mn = 100;
	    }
	    else disWorstTime = 0;
	  }
	  else disWorstTime = 0;
#endif
#ifdef FindS_DISCARD_SF_CLONES
	  if (nScifis+nVSATs>=2 && nPsSfs) {  // Init scifi hit pattern
	    for (int jpsf = 0; jpsf<nPsSfs; jpsf++) hitsPsSfs[nCSfs][jpsf] = 0;
	  }
#endif

	  // *************** ReEVALUATE TRACK CANDIDATE ***************
#ifdef TARGET_POINTING_SPACE
	  static float yFirst, zFirst;
#endif
#ifdef PR_USE_GEM_CORR
	  //          ***** GEM (+STRIPPED pixelEGM) SPECIFIC *****
	  int nQuads = 0, nFullQuads = 0, quadStatus = 0;
	  static int prvQuad; // Quadruplet# of previous GEM hit
#endif

	  int nYInfo = 0;  static float xYFirst, xYSecond, xYLast, xYPenUlt;

#define FindS_ENHANCE_XPROJ 1
#ifdef FindS_ENHANCE_XPROJ
	  int nXProjFound = 0; float xProjQual = 0;
#endif
	  static float x; float xPrv = -100;
	  bool first = true; static float xFirst;
#define FindS_RESCUE_VERY_GOOD
#ifdef FindS_RESCUE_VERY_GOOD
	  // Very good tracks (#hits>>,chi2<<) can be outperformed by ghosts
	  // which take advantage of (i) the VSAT bonus and (ii) widely opened
	  // routes. We try here to rescue them: cancelling the bonus and
	  // allowing good chi2 to contribute more than 1 unit to the quality
	  // function. Both steps are risky => therefore we strongly restrict
	  // the rescue's scope. Trying to avoid cases where a ghost would
	  // outperform a genuine track while benefitting from enhanced VSATs
	  // because just below rescue criterion.
	  float bonus;
	  if      (nFound<nHitsMn_Mx+1) bonus = incrVSAT;
	  else if (nFound<nHitsMn_Mx+2) bonus = incrVSAT*3/4;
	  else if (nFound<nHitsMn_Mx+3) bonus = incrVSAT/2;
	  else if (nFound<nHitsMn_Mx+4) bonus = incrVSAT/3;
	  else                          bonus = 0;
#else
#define bonus incrVSAT
#endif
	  int nBonus;
#define FindS_DISCARD_WORST_DIFF
#ifdef FindS_DISCARD_WORST_DIFF
	  // This feature is enabled ony for the scifiSi telescope case. There
	  // the tolerance about Si planes has to be enlarged to make for some
	  // (most?) of the routes being defined by scifis, which have a much
	  // worse resolution. A drawback is that some ghost hits can also be
	  // picked up and ruin the quality of the candidate track. There are
	  // several steps further down the processing that will tend to discard
	  // ghost hits, starting w/ the fit performed w/in "FindSpace" itself (
	  // and this sytematically for zone 0x10, cf. supra) that will redefine
	  // the track's parameters. They will be particularly effective if a
	  // single hit is questionable. Therefore let's remove the worst hit
	  // from the computation of the quality of the track, if it's (i) very
	  // bad and (ii) much (twice) worse than the rest. In the particular
	  // case if the scifi/Si this can give a net (one hit less but better
	  // chi2) to the candidate track.
	  float worstDev = 0, nextWDev = 0; // Worst, next to worst deviations
	  static int iplWorst; // Index of plane of worst deviation
#endif
	  // Min distance between two drift planes that are not associated: the
	  // idea is to give hits less weight when they are associated in a pair
	  // from staggered drifts planes. (For if one falls w/in route, so does
	  // the other.). A lower limit is set by the MB case, where interplane
	  // distance is ~2.8 cm. The upper limit is not much constrained.
	  static float deltax = 3.0;
	  int mFit = nFit; // Fit: linear (==4) or quadratic (==5)
	  if (nVSATs && mFit) mFit = 4;// VSATs: expecting high momenta =>linear
	  if (doVSATFit && nVSATs+nScifis<=TOpt::iCut[31]) mFit = 0;
	  float nAbscissae; unsigned int projs, projPrv; int nProjs;
	  float nAbscissaeRICH; int nProjsRICH, nScifisRICH;
	  int nDrifts;
	  for (ipl = 0, qual = nFound, nAbscissae = 0, projs=projPrv = 0,
		 nProjs=nProjsRICH = 0, nAbscissaeRICH = 0,
		 nScifisRICH=nDrifts=nBonus = 0; ipl<npl; ipl++) {
	    if (!(1<<(ipl%32)&plPats[ipl/32]))
	      continue;         // ********** LOOP on PLANES in TRACK **********
	    int k0 = hitsList[ipl];
	    float c = cosa[ipl], s = sina[ipl]; x = xpl[ipl];
#ifdef QUADRATIC_PROJ
	    yy = yq*x*x+yp*x+y0; zz = zp*x+z0;
#else
	    yy =        yp*x+y0; zz = zp*x+z0;
#endif
	    unsigned int proj = 1<<prpl[ipl];

	    int tysp = typl[ipl], type = tysp&0xff; // ***** DETECTOR TYPE *****
	    if (tysp&0x100) {
	      if (x>xPrv+deltax || proj!=projPrv)
		nAbscissae += 1.01;
	      else if ((zone==1 || zone==3) && (fabs(yy)>65 || fabs(zz)>65))
		nAbscissae += .751;       // Special case of outer region: little redundancy otherwise => give both layers of double layer large weight
	      /*
		else if (type==15)
		nAbscissae += .501;       // DC's 
	      */
	      else
		nAbscissae += .501;       // STs, DWs, MBs (inner region)
	      nDrifts++;
	    }
	    else {
	      if (tysp==scifiType) {      // Scifis
		if (disWorstTime) {
		  float scft = hinfos[k0];
		  if ((disWorstTime>0 && scft>tCut) ||
		      (disWorstTime<0 && scft<tCut)) {
		    nFound--; nScifis--; qual -= 1; continue;
		  }
		  // If "disWorstTime", track timing was reset => re-evaluate 
		  st += scft; st2 += scft*scft;
		  if (scft<t_mn) t_mn = scft;
		  if (scft>t_mx) t_mx = scft;
		}
		nAbscissae += 1.01; nBonus++;
	      }
#ifdef PR_ENHANCE_CENTRAL_GEMs            // Central GEMs
	      else if (tysp==(26|0x200) &&
		       fabs(yy-yCentralGEM[ipl])<4.0 &&
		       fabs(zz-zCentralGEM[ipl])<4.0) {
		nAbscissae += 1.01; nBonus++;
	      }
#endif
	      else if (tysp&0x800) {      // PixelGEM/MM: pixels piece
		nAbscissae += 2.01;
		if (tysp&0x200) { nBonus++; if (type==29) nBonus++; }
		if (!(projs&0x2)) { // If not yet Z-proj.: include it
		  nProjs++; projs |= 0x2;
		}
		float v = zz*c - yy*s; qual -= fabs(v-hinfos[k0])/tols[ipl];
	      }
	      else if (type==28) {        // PixelGEM: strips piece
		nAbscissae += 1.01; if (tysp&0x200) nBonus++;
	      }
	      else nAbscissae += 1.01;
	    }
#ifdef FindS_DISCARD_SF_CLONES
	    if (tysp&0x200) { // This condition may differ from the one above...
	      // ...if !vSATEnhanced or "PR_ENHANCE_CENTRAL_GEMs"
	      if (nPsSfs) {
		int k0sf = fh_sfs[ipl]+k0-fh[ipl];  // Scifi hit index
		int ipsf = k0sf/32;                 // Scifi pattern index
		hitsPsSfs[nCSfs][ipsf] |= 1<<k0sf%32;
	      }
	    }
#endif
#ifdef PR_USE_GEM_CORR
	    // Note: This is not proved to bring any improvement. In addition,
	    // it's inconsisten: what if GEM amplitude correlations are not
	    // enabled by option?
	    if (type==26 || type==28) {
	      int quad = ihinfos[k0]%0x10000;
	      if (x>xPrv+3) { // New GEM station => reset quad status
		if (quadStatus>0) { // Conclude on previous GEM station
		  nQuads++; if (quadStatus==4) nFullQuads++;
		}
		prvQuad = quad; quadStatus = quad>0 ? 1 : 0;
	      }
	      else if (quad>0) {
		if (quad==prvQuad) quadStatus++; // Progressing w/in quad
		else quadStatus = 1;// Entering a quad, but not at its beginning
	      }
	    }
#endif

	    //                      ***** UPDATE: q function / projections *****
	    float u = yy*c + zz*s, dev = fabs(u-Uhit[k0])/tols[ipl];
#ifdef FindS_DISCARD_WORST_DIFF
	    if (igr==4 &&
		// Disregard scifis for simplicity's sake and because the
		// target of "FindS_DISCARD_WORST_DIFF" is rather Si's.
		tysp!=scifiType) {
	      if (dev>worstDev) {
		nextWDev = worstDev; worstDev = dev; iplWorst = ipl;
	      }
	      else if (dev>nextWDev) nextWDev = dev;
	    }
#endif
	    qual -= dev;
	    if (!(proj&projs)) { nProjs++; projs |= proj; }
	    if (zone==1 && x<xRICH) {
	      nProjsRICH = nProjs; nAbscissaeRICH = nAbscissae;
	      if (type==22) nScifisRICH++;
	    }
#ifdef FindS_ENHANCE_XPROJ
	    if (proj&xProjs) { // Contributions of ``cross''-proj.
	      nXProjFound++; xProjQual += fabs(u-Uhit[k0])/tols[ipl];
	    }
#endif

	    if (first) {          // ***** FIRST DETECTOR in TRACK: INIT's *****
	      first = false; xFirst = x;
#ifdef TARGET_POINTING_SPACE
	      yFirst = yy; zFirst = zz;
#endif
	      if (mFit) {
		int j, k, jk, i5 = 5*ipl; float U = Uhit[k0]/wghts[ipl];
		for (j=jk = 0; j<mFit; j++) {
		  for (k = 0; k<=j; k++) Ajk[jk++] = aij[i5+j]*aij[i5+k];
		  AjV[j] = aij[i5+j]*U;
		}
		Vi2 = U*U;
		if (type==29) { // PixelGP (and not PixelMP, where the info...
		  // ...along the v-coordinate is not strongly constrained)
		  // v = u + pi/2, cf. TEv::PrePattern2.
		  i5 = 5*ipl;
		  double bij[5] = { aij[i5+2], aij[i5+3], -aij[i5], -aij[i5+1], -aij[i5+1]*x };
		  float V = hinfos[k0]/wghts[ipl]; for (j=jk = 0; j<mFit; j++) {
		    for (k = 0; k<=j; k++) Ajk[jk++]+= bij[j]*bij[k];
		    AjV[j]+= bij[j]*V;
		  }
		  Vi2+= V*V;
		}
		if (nYInfo==0 && fabs(s)<.8) {
		  xYFirst = x; nYInfo = 1;  // First hit w/ Y info
		}
	      }
	    }
	    else if (mFit) {  // ***** NEXT DETECTORS: PREPARE for FITTING *****
	      int j, k, jk, i5 = 5*ipl; float U = Uhit[k0]/wghts[ipl];
	      for (j=jk = 0; j<mFit; j++) {
		for (k = 0; k<=j; k++) Ajk[jk++]+= aij[i5+j]*aij[i5+k];
		AjV[j]+= aij[i5+j]*U;
	      }
	      Vi2+= U*U;
	      if (type==29) { // PixelGP, cf. supra
		i5 = 5*ipl;
		double bij[5] = { aij[i5+2], aij[i5+3], -aij[i5], -aij[i5+1], -aij[i5+1]*x };
		float V = hinfos[k0]/wghts[ipl]; for (j=jk = 0; j<mFit; j++) {
		  for (k = 0; k<=j; k++) Ajk[jk++]+= bij[j]*bij[k];
		  AjV[j]+= bij[j]*V;
		}
		Vi2+= V*V;
	      }
	      if (fabs(s)<.8) {                  // First/last hit w/ Y info
		if      (nYInfo==0) xYFirst  = x;
		else if (nYInfo==1) xYSecond = x;
		else {
		  xYPenUlt = xPrv; xYLast = x;
		}
		nYInfo++;
	      }
	    }
	    xPrv = x; projPrv = proj;                     // ***** UPDATEs *****
	  }
	  if (nBonus>1) nAbscissae += nBonus*bonus;
	  if (doSIFit &&!(doLinFit || doQuadFit) && // SI-track in zone 0x1...
	      nDrifts) // ...try to avoid fitting low momenta => no drifts
	    mFit = 0;

	  if (nFound<4) continue; // This may(?) happen 'cause of worst time cut

	  if (nAbscissae<nHitsMn) // If drifts subdued: this must be reevaluated
	    continue;

	  if (nDrifts && nDrifts<=3 && // Track has drift hits but very few...
	      // ...which is suspicious for drift stations have many planes...
	      nAbscissae-nDrifts<nHitsMn) {
	    // ...And in addition, it passes PR requirement only thanks to these
	    // (or would hardly do so w/o them, cf. "nAbscissae" increment for
	    // drift planes). Cannot be mistaken for a track travelling in the
	    // SAS/LAS overlap, because few hits overall...
	    if (nScifis || nDrifts<3)
	      continue;//=> Assumed to be fake track picking up ghost drift hits
	  }

	  if (nDrifts==nFound && igr==0 && nAbscissae<nHitsMn+1)
	    // Special case of DC-only track in 2006 setup
	    continue;

	  // Isolated hit: it has a high probability of being a random pickup.
	  // (e.g. case of a ST04 hit separated from the rest by Rich thickness)
	  // => disregard it. At this point, this will only affects the
	  // quality value associated to the track: the cuts on the #hits have
	  // already been performed, the fitting arrays already filled. In
	  // addition the isolated hit will (or rather may, depending on fit)
	  // be reincorporated in the 2nd loop on combinations.
	  if (x-xPrv>150 &&
	      (zone==1 || zone==2)) // Exclude zone #3, because of HI05
	    nAbscissae -= 1.01;
	  // Along the same line: Isolated hits ahead of RICH. Because of the
	  // large separation between up- and down-stream of RICH, the random
	  // pick up probability is particularly high there. And it can concern
	  // more than one hit.
	  if (nProjsRICH==1 &&
	      !nScifis && qual/nFound<.66) {
	      // ...Cf. evt #50343914 of 03P1E wb dump. 
	    nAbscissae -= nAbscissaeRICH; xFirst = xRICH;
	  }
	  if (igr==1 && nFound<16 && mode<=0 && nBonus>6) {
	    // ...Special case of 2010-like VSAT setup: scifis semi-enhanced
	    // upon 1st iter. =>
	    //  - Favour tracks w/ both FI05X and Y.
	    //  - Unfavour isolated hits, not FI05, hits.
	    if      (nProjsRICH==2 && nScifisRICH==2) nAbscissae += 1;
	    else if (nProjsRICH && nProjsRICH<=2 && !nScifisRICH) continue;
	    else if (nProjsRICH==1) {
	      nAbscissae -= nAbscissaeRICH; xFirst = xRICH;
	    }
	  }

#if defined FOUND_vs_EXPECTED && FOUND_vs_EXPECTED == 2
	  //float eff = hasdrift ? tracking_eff/1.2 : tracking_eff;
	  float eff = tracking_eff;
	  if (nFound<eff*(nexp1-nexp0) && nFound<2.5*nHitsMn)
	    continue;
	  // Apply a stricter requirement to ``long'' tracks in zone 0x1
	  // (it's done in 2 steps, w/ looser step introduced to relax the
	  // requirement on tracks through FI05:6+GM05:6, e.g.
	  // mu-track in evt #709 of dstar.2002.04.fz.1 (N.B.: upon earliest
	  // iter's this kind of track has nScifi==0, cf. (SCIFIS_ENHANCED"))
	  if (zone==1 && nScifis==0 && xFirst<xGM04 && xPrv>xGM04 &&
	      (xFirst<xFI05 && nAbscissae<nHitsMn_FI05 ||
	       nAbscissae<nHitsMn_GM04))
	    continue;
#endif

	  bool vsatTrack = igr==4 ||
	    zone==0 && nScifis+nVSATs>nFound/2 ||
	    zone>=1 && nScifis+nVSATs>=nFound/2; // E.g. 4*FI56+4*GM8
	  if (vsatTrack) {
	      // ***** SciFiSi AND SCIFI-(ALMOST)ONLY TRACKS ***** Why here?
#ifdef SCIFIS_t_PROB         // ***** CUT ON TIME CHI2 PROBA *****
	    if (disWorstTime) {
	      if (nScifis>=2) { // So that, in any case, infra is meaningful...
		stProb = TMath::Prob(st2-st*st/nScifis,nScifis);
		if (stProb<2.5e-2) continue;// ~=Prob(5,1)=>5 sigma if 1ns resol
	      }
	      else stProb = 0;
	    }
	    else if (kScOff>=0) {
	      float scft = hinfos[kScOff]; st += scft; st2 = +scft*scft;
	      nScifis++; stProb = TMath::Prob(st2-st*st/nScifis,nScifis);
	    }
#else
	    devt = sqrt(st2*nScifis-st*st)/nScifis; // ?? What if "nScifis==0"
#endif
	    // ********** CUT on SCiFi TIME SPREAD **********
	    if (nScifis && t_mx-t_mn>scfDt_mx) continue;

	    int icsf; // ********** SCIFI COMBI **********
#ifdef FindS_DISCARD_SF_CLONES
	    if (nPsSfs) {
	      int match, iqual = int(qual*10+.5);
	      for (icsf=match = 0; icsf<nCSfs; icsf++) {// Loop on scifis combis
		if (nScifis+nVSATs!=nScSfs[icsf]) continue;
		int ipsf; for (ipsf = 0, match = 1; ipsf<nPsSfs; ipsf++) {
		  if (hitsPsSfs[nCSfs][ipsf]!=hitsPsSfs[icsf][ipsf]) {
		    match = 0; break;
		  }
		}
		if (match) {
		  // Check that the 2 combinations match in total #hits,
		  // including the non-VSAT ones.
		  if (nFound==nFoundSfs[icsf]) break;
		  else match = 0;
		}
	      }
	      if (match && iqual>iQualSfs[icsf] && nScifis+nVSATs>7)
		match = 0;
	      if (match) continue;
	      else { // New scifis combi w/ hit pattern not yet encountered
		if (nCSfs<nCSfsMx) {
		  nFoundSfs[nCSfs] = nFound; iQualSfs[nCSfs] = iqual;
		  icsf = nCSfs; nCSfs++;
		}
		else { icsf = -1; cSfsOver++; }
	      }
	    }
	    else icsf = -1;
#else
	    if (nCSfs<nCSfsMx) {
	      icsf = nCSfs; nCSfs++;
	    }
	    else { icsf = -1; cSfsOver++; }
#endif
	    sfs_only[ncomb] = icsf;
	    if (icsf>=0) {
	      timesSfs[icsf] = st/nScifis; usedSfs[icsf] = 0;
	      nScSfs[icsf] = nScifis+nVSATs;
	    }
	  }
	  else sfs_only[ncomb] = -1;

	  // 2-projection tracks are not reliable. Scifis-tracks are excepted
	  // because they have access to limited redundancy and can easily
	  // fail to collect 3 different projections. They are expected to
	  // gain reliability through bridging.
	  if (nScifis+nVSATs<nFound/2 && nProjs<3) continue;

	  int ifail_qs = -1; static double chi2; static int NDF;
	  if (mFit && nFound>4) {
	    HepVector YZ(5);
	    if (mFit==5 && nFound>8 && nYInfo>=4 &&
		nScifis<4) {   // Skip scifi tracks: they must be straight
	      HepSymMatrix M(5); HepVector V(5);
	      int j, k, jk;
	      for (j=jk = 0; j<5; j++) {
		for (k = 0; k<=j; k++) M[j][k] = Ajk[jk++];
		V[j] = AjV[j];
	      }
	      M.invert(ifail_qs);
	      if (ifail_qs==0) {
		YZ = M*V;
		for (j=jk = 0, chi2 = Vi2; j<5; j++) {
		  double Sj = 0;
		  for (k = 0; k<j; k++) Sj += YZ[j]*YZ[k]*Ajk[jk++];
		  Sj -= YZ[j]*AjV[j];
		  chi2 += 2*Sj + Ajk[jk++]*YZ[j]*YZ[j];
		}
	      }
	      NDF = nFound-5;
	      float xlow = xYSecond-0.5*(xYPenUlt-xYSecond);
	      float xup  = xYSecond+1.5*(xYPenUlt-xYSecond);
	      xYLow[ncomb] = xYFirst<xlow ? xYFirst : xlow;
	      xYUp [ncomb] = xYLast >xup  ? xYLast  : xup;
	    }
	    else if (mFit==4) {
#ifdef FindS_DISCARD_WORST_DIFF
	      if (igr==4 && worstDev>.8 && nextWDev<.4) {
		int ipl = iplWorst, k0 = hitsList[ipl];
		nAbscissae -= 1.01; nFound--; plPats[ipl/32] &= ~(1<<(ipl%32));
		int j, k, jk, i5 = 5*ipl; float U = Uhit[k0]/wghts[ipl];
		for (j=jk = 0; j<mFit; j++) {
		  for (k = 0; k<=j; k++) Ajk[jk++]-= aij[i5+j]*aij[i5+k];
		  AjV[j]-= aij[i5+j]*U;
		}
		Vi2-= U*U;
	      }
#endif
	      HepSymMatrix M(4); HepVector V(4);
	      int j, k, jk;
	      for (j=jk = 0; j<4; j++) {
		for (k = 0; k<=j; k++) M[j][k] = Ajk[jk++];
		V[j] = AjV[j];
	      }
	      M.invert(ifail_qs);
	      if (ifail_qs==0) {
		YZ = M*V; YZ[4] = 0;
		for (j=jk = 0, chi2 = Vi2; j<4; j++) {
		  double Sj = 0;
		  for (k = 0; k<j; k++) Sj += YZ[j]*YZ[k]*Ajk[jk++];
		  Sj -= YZ[j]*AjV[j];
		  chi2 += 2*Sj + Ajk[jk++]*YZ[j]*YZ[j];
		}
	      }
	      NDF = nFound-4;
	    }
	    if (ifail_qs==0) {
	      if (chi2/NDF>chi2_mn) continue;
	      // ***** STORE COMBINATION
	      Ysl[ncomb] = YZ[3]; Y0[ncomb] = YZ[2];
	      if (nFound-NDF==5) Y2[ncomb] = YZ[4];
	      else               Y2[ncomb] = 0;
	      Zsl[ncomb] = YZ[1]; Z0 [ncomb] = YZ[0];
	      qual = 2*nFound*(TMath::Prob(chi2*chi2_scale,NDF));
	      if (igr==4)  // Case of beam track: do put more weight on quality,
		// particularly when, as is the case here, it's determined by
		// a linear fit. This (more weight on quality) may be favorable
		// in any case, but it's not been checked. Let's do it in the
		// case of beam track, where the situation is better under
		// control. In any case, it fixes the case of evt #20975748 in
		// 04/W28/cdr09002-37277.raw...
		qual *= 1.5;
	    }
	  }
	  if (ifail_qs)
	  {                           // ********** STORE COMBINATION **********
	    Ysl[ncomb] = yp; Y0 [ncomb] = y0;
#if   defined QUADRATIC_PROJ
	    Y2 [ncomb] = y2;
#else
	    Y2 [ncomb] = 0;
#endif
	    Zsl[ncomb] = zp; Z0 [ncomb] = z0;
	  }

	  //#define DEBUG_SPACE_FIT
#ifdef DEBUG_SPACE_FIT
	  static int idebug = -1; if (idebug>=0 && zone==idebug) {
	    double chi2p;
	    for (ipl = 0, chi2p = 0; ipl<npl; ipl++) {
	      if ((1<<(ipl%32)&plPats[ipl/32])==0) continue;
	      int k0 = hitsList[ipl];
	      double y = Ysl[ncomb]*xpl[ipl]+Y0[ncomb];
	      double z = Zsl[ncomb]*xpl[ipl]+Z0[ncomb];
	      double u = y*cosa[ipl]+z*sina[ipl];
	      double du = (Uhit[k0]-u)/wghts[ipl]; chi2p += du*du;
	      printf("%2d %4d %.2f  ",ipl,k0,du);
	    }
	    printf("\n    %d  %f  %f  %f\n",nFound,chi2,chi2p,qual);
	  }
#endif

	  // ***** QUALITY ESTIMATOR = Quality factor (less then 1) +...
	  Qt [ncomb] = nAbscissae+qual/nFound;         // ...+ "nAbscissae"
	  //#define FindS_PROMOTE_CHI2
#ifdef FindS_PROMOTE_CHI2
	  float fHits = (nexp1-nexp0)/(float)nHitsMn-1;
	  Qt [ncomb] += qual/nFound*(qual/nFound-.5)*fHits*fHits;
#endif
#ifdef FindS_ENHANCE_XPROJ
	  if ((projs&uv10Projs) && // If 10 degrees UV's alone are relating...
	      (projs&zProjs) && !(projs&UVProjs))  // ...Z to the rest...
	    // ...enhance the contrib of proj. that are neither "ip1" nor "ip2".
	    Qt [ncomb] -= FindS_ENHANCE_XPROJ*xProjQual/nXProjFound;
#endif
	  if (nScifis==nFound)
#ifdef SCIFIS_t_PROB
	    Qt [ncomb] += TOpt::dCut[63]*(stProb-.25); // ...+/- scifi time prob
#else
	    Qt [ncomb] += TOpt::dCut[63]*(1-devt);     // ...+ scifi time rms
#endif

#ifdef TARGET_POINTING_SPACE
	  //Qt [ncomb] += atan2(1,fabs(yp))/M_PI;
	  //Qt [ncomb] += atan2(1,fabs(zp))/M_PI;
	  Qt [ncomb] -= trgtPtngY*fabs(yFirst/(xFirst+X0)-yp);
	  Qt [ncomb] -= trgtPtngZ*fabs(zFirst/(xFirst+X0)-zp);
#endif

	  //#define OLD_FOUND_vs_EXPECTED
#ifdef OLD_FOUND_vs_EXPECTED
	  Qt [ncomb] += ((float)nFound*10)/nexpected;
#else
	  // ***** ``HIT EFFICIENCY`` *****
	  if (beamTinDY)
	    // Consider that all tracks span the full length of the beamT (and
	    // do not enter the zone through the edge). The assumption is
	    // reasonable in the DY setup where the beam is very parallel.
	    Qt [ncomb] += nFound*10./npl;
	  else if (zone) { // ...refrain from applying it to zone 0x1
	    Qt [ncomb] += nFound*6.7/nexpected;
	    if (nexp1<nexpected-nexp0)
	      Qt [ncomb] += nFound*3.3/nexp1;
	    else
	      Qt [ncomb] += nFound*3.3/(nexpected-nexp0);
	  }
	  else
	    Qt [ncomb] += nFound*10./nexpected;
#endif

#ifdef FindS_RESCUE_VERY_GOOD
	  if (qual/nFound>.9 && mode<=1 &&  // Good chi2 (while narrow route)...
	      nFound>nHitsMn_Mx+1 &&        // ...and #hits>>
	      nFound>=nexpected-1)          // ...and missing only one hit
	    Qt [ncomb] += .8+10./nexpected; // ...make(almost) for missing hit
#endif

#ifdef PR_USE_GEM_CORR
	  if (quadStatus>0) { // Conclude on previous GEM station
	    nQuads++; if (quadStatus==4) nFullQuads++;
	  }
	  if (nQuads) Qt [ncomb] += (float)nFullQuads/nQuads -.5;
#endif
#ifdef PR_2ZONE_ZPROJ                      // Special 2-zone proj. search in Z
	  if (ip2>=ip2_2Z) Qt [ncomb] += 2; //  => give it a bonus
#endif
	  if (TOpt::Hist[6]) { // ***** HISTOGRAMS *****
	    static CsHist2D *tYP_FS[5], *tZP_FS[5], *tCH_FS[5], *tQt_FS[5];
	    static bool book = true; if (book) {
	      book = false;
	      CsHistograms::SetCurrentPath("/Traffic/PrePattern");
	      for (int jgr = 0; jgr<5; jgr++) {
		char name[] = "tYP_FSi";
		double xmax = TOpt::dPRpar[1+10*jgr]*1.2;
		if (jgr==0) xmax *= 2;
		sprintf(name,"tYP_FS%d",jgr);
		tYP_FS[jgr] = new CsHist2D(name,name,50,-xmax,xmax,5,-1,3);
		sprintf(name,"tZP_FS%d",jgr);
		tZP_FS[jgr] = new CsHist2D(name,name,50,-xmax,xmax,5,-1,3);
		sprintf(name,"tCH_FS%d",jgr);
		tCH_FS[jgr] = new CsHist2D(name,name,50,0,1,5,-1,3);
		sprintf(name,"tQt_FS%d",jgr);
		tQt_FS[jgr] = new CsHist2D(name,name,50,0,2,5,-1,3);
	      }
	      CsHistograms::SetCurrentPath("/");
	    }
	    if (igr<5) {
	      (tYP_FS[igr])->Fill(yp,(double)mode);
	      (tZP_FS[igr])->Fill(zp,(double)mode);
	      if (ifail_qs==0)
		(tCH_FS[igr])->Fill(chi2/NDF,(double)mode);
	      else
		(tQt_FS[igr])->Fill((nFound-qual)/nFound,(double)mode);
	    }
	  }
	  ind[ncomb] = ncomb+1;
	  if (++ncomb>=maxcomb) {
	    CsErrLog::msg(elError,__FILE__,__LINE__,
	      "Too many combinations %d in zone 0x%x",ncomb,1<<zone);
	    goto sorting;
	  }
#ifdef FindS_WO_1ST_SCIFI
	  if (nScifis==1 && // Few scifi hits: enable "w/o 1st plane" feature...
	      // ...it's restricted to case zone==4 w/ nScifiDets>=4
	      zone==4 && nScifiDets>=4) {
	    for (ipl = 0; ipl<2; ipl++) if (xpl[ipl+1]>xFirst+.1) {
	      ipl1stScifi = ipl; break;
	    }
	    if (ipl1stScifi+1) nIters = 2;
	  }
	  }
#endif
	} // End of loop over tracks on 2-d base proj.
      } // End of loop over tracks on 1-st base proj.

    } // End of loop over 2-d base proj
  } // End of loop over 1-st base proj

 sorting:                       // ********** SORT Qt-WISE **********
  if (cSfsOver)
    CsErrLog::msg(elError,__FILE__,__LINE__,
		  "Too many scifi combinations %d",nCSfsMx+cSfsOver);
  SORTZV(Qt[0],ind[0],ncomb,1,1,0);


  int m; int nt = 0, nh = 0, nHsMn_trk = 0;
  int hh[npl], iflh[npl] /* Flags hits used several times */;
  for (int l = 0; l<ncomb; l++) {
    // ********** LOOP ON COMBINATIONS, starting from Qt max **********

    m = ind[l]-1; int icsf = sfs_only[m];

    // ***** MINIMUM # of HITS, MAXIMUM # of COMMON HITS *****
    int nused_mx = TOpt::iPRpar[igr*10+1]; // And in second iter?!
    if (igr==4) {
      // Special case of zone 0x10: "nHitsMn" has been loosen in the loop
      // over all combinations supra => restore it now that the search route
      // is more precisely defined (by the linear fit performed in 1st loop).
      nHsMn_trk = TOpt::iPRpar[igr*10+2 + mode_par];
    }
    else if (icsf>=0) {
      if      (zone==0) nHsMn_trk = 5;
      else if (zone==1) {
	if (TOpt::iCut[31]) {
	  if (mode<=0) nHsMn_trk = 8;
	  else         nHsMn_trk = 7;
	}
	else           nHsMn_trk = 4;
      }
      else if (zone<=3) nHsMn_trk = 4;
      else              nHsMn_trk = nHitsMn;
      int nSc = nScSfs[icsf];
      if (zone!=1 ||   // Exclude case of scifi_enhanced in zone 0x2,...
	  mode>0) {    // ...upon 1st iters ("mode<=0"), cf. "iCut[31]".
	if (usedSfs[icsf]==0) {
	  usedSfs[icsf] = 1;
	  nused_mx = npl; // Ensures combi is not rejected because of "nused".
	}
	else if (nSc>=5) {
	  nused_mx = 6;
	  // E.g., in zone 0x2, FI06X+U+V+FI05X
	  if (nSc>=7) nused_mx = 7;
	}
	else if (nSc==4 && (1<=zone && zone<=5))
	  nused_mx = 4; // E.g. FI07X+Y  +FI08X
      }
#ifdef DEBUG_SCIFIS_TIMING
      if (idebug) printf("++++++++ %d\n",m);
#endif
    }
    else nHsMn_trk = nHitsMn;

    float nPixels; int nScOffs;
    for (ipl=nFound=nused=nScOffs = 0, nPixels = 0; ipl<npl; ipl++) {
      // ********** HITS on _ALL_ PLANES IN DEFINED ROUTE **********
      if (fh[ipl]==lh[ipl]) continue;

      float xx = xpl[ipl];                  // ***** COORD @ CURRENT PLANE *****
      yy = Y0[m]+Ysl[m]*xx; zz = Z0[m]+Zsl[m]*xx;
      yy += Y2[m]*xx*xx;   // !? what about quadratic proj?
      if (Y2[m] && sina[ipl]<.8 && (xx<xYLow[m] || xYUp[m]<xx)) continue;
      if (!setup.vDetect(idpl[ipl]).InActive(yy,zz)) continue;
      float u = yy*cosa[ipl] + zz*sina[ipl]; // rotate

      //                         ***** ANY HIT on CURRENT PLANE w/in ROUTE *****
      // No longer use enlarged Si tolerances, now that routes are defined
      // by a least square fit (and no longer, possibly, by, bad resolution,
      // scifis)
      float tol = wghts[ipl];
      if (u<Uhit[fh[ipl]]  -tol ||                // ***** HITS w/in ROUTE *****
          u>Uhit[lh[ipl]-1]+tol) continue;

      static int k0; int kl = fh[ipl], kr = lh[ipl]-1;// ***** CLOSEST HIT *****
      while (kr-kl>1) {
	k0 = (kl+kr)>>1;
	if (u<=Uhit[k0]) kr = k0;
	else             kl = k0;
      }
      if (Uhit[kr]-u<u-Uhit[kl]) k0 = kr;
      else                       k0 = kl;
      if (fabs(u-Uhit[k0])>tol) continue;

      int tysp = typl[ipl];
      if (icsf>=0 && tysp==(22|0x200)) {           // ***** SCIFI SPECIFIC *****
#ifdef DEBUG_SCIFIS_TIMING
	if (idebug) {
	  if (nFound==0) printf("======== %d\n",m);
	  printf("%3d %3d %7.3f  %3d %3d\n",ipl,k0,hinfos[k0],kr,kl);
	}
#endif
	// Now that mean time is known, cut @ HALF "scfDt_mx"
	if (fabs(hinfos[k0]-timesSfs[icsf])>scfDt_mx/2) {
	  // Currently selected hit is out of time => try other one?
	  int k0p, pm; for (pm = 0, k0p = -1; pm<2; pm++) {
	    int kp; if (pm) { kp = k0-1; if (kp<fh[ipl])  continue; }
	    else            { kp = k0+1; if (kp>=lh[ipl]) continue; }
	    if (fabs(u-Uhit[kp])<tol) {
	      if (fabs(hinfos[kp]-timesSfs[icsf])>scfDt_mx/2) {
		if (k0p<0 || fabs(u-Uhit[kp])<fabs(u-Uhit[k0p])) k0p = kp;
	      }
	    }
	  }
	  if (k0p>=0) k0 = k0p;
	  else if (beamTinDY>1 && !nScOffs)
	    nScOffs++; // Special case of beamT in DY: cf. supra
	  else continue;
#ifdef DEBUG_SCIFIS_TIMING
	  if (idebug) {
	    printf("%3d %3d %7.3f  %3d %3d\n",ipl,k0,hinfos[k0],kr,kl);
	  }
#endif
	}
	iflh[nFound] = -1;
      }
      else if (tysp&0x800) {    // ***** PIXEL PixelGEM/MM SPECIFIC *****
	float v = zz*cosa[ipl] - yy*sina[ipl];
	float tolv; if (tysp==(29|0x800)) tolv = tol;
	else { // Tolerance for M-pixels: special, variable v-precision
	  // Let's apply the outer tolerance w/ a somewhat arbitrary
	  // margin, cf. explanations supra.
	  if (fabs(yy)>2.4 || fabs(zz)>2.4) tolv = tolMPvO;
	  else                              tolv = tolMPvI;
	}
	//if (fabs(v-hinfos[k0])>tol || ifl[k0]>0) {
	if (fabs(v-hinfos[k0])>tolv) {
	  // Currently selected hit incompatible w/ V-wise info: try other one.
	  // ( - We only consider here the 2 closest neighbours, closest to the
	  //    ``left'' and closest to the ``right'', which is not optimum...
	  //   - We do not require the pixel hit to be free ("ifl" flag), thus
	  //    allowing all combinations of VSATs.)
	  if (kr==kl) continue;
	  if (k0==kr) k0 = kl;
	  else        k0 = kr;
	  if (fabs(v-hinfos[k0])>tolv) continue;
	  if (fabs(u-Uhit[k0])>tol) continue;
	}
	nPixels += tysp==(32|0x800) && TOpt::ReMode[48] ? 1.51 : 2.01;
	iflh[nFound] = -1;
      }
#ifdef PR_ST_WHICHSIDE_INFO
      else if (tysp&0x400) {                       // ***** CHECK WHICHSIDE-INFO
	float whichSide = hinfos[k0]; if (whichSide) {
	  float v = zz*cosa[ipl] - yy*sina[ipl];
	  if (whichSide*v<0) {
	    // Currently selected hit not compatible w/ which-side info
	    if (kr==kl) continue;
	    if (k0==kr) k0 = kl;
	    else        k0 = kr;
	    whichSide = hinfos[k0]; if (whichSide && whichSide*v<0) continue;
	    if (fabs(u-Uhit[k0])>tol) continue;
	  }
	}
	iflh[nFound] = -1;
      }
#endif
#ifdef PR_GP_WHICHSIDE_INFO
      else if ((tysp&0xff)==28) {  // ***** WHICH-SIDE INFO in STRIPPED PixelGEM
	float whichSide = ihinfos[k0]>>16; if (whichSide) {
	  float v = zz*cosa[ipl] - yy*sina[ipl] - offsets28[ipl-ipl28];
	  if (whichSide*v<-.5) { // Cut on the distance to the equator
	    // Currently selected hit not compatible w/ which-side info
	    if (kr==kl) continue;
	    if (k0==kr) k0 = kl;
	    else        k0 = kr;
	    whichSide = ihinfos[k0]>>16; if (whichSide && whichSide*v<0)
	      continue;
	    if (fabs(u-Uhit[k0])>tol) continue;
	  }
	}
	iflh[nFound] = -1;
      }
      else if ((tysp&0xff)==31) {  // ***** WHICH-SIDE INFO in STRIPPED PixelMM
	float whichSide = ihinfos[k0]>>16; if (whichSide) {
	  // For pixelMMs, which are expected to be in zonne 0x1, we
	  // assume the offset = 0. Otherwise, same processing as supra.
	  float v = zz*cosa[ipl] - yy*sina[ipl];
	  if (whichSide*v<-.5) { // Cut on the distance to the equator
	    if (kr==kl) continue;
	    if (k0==kr) k0 = kl;
	    else        k0 = kr;
	    whichSide = ihinfos[k0]>>16; if (whichSide && whichSide*v<0)
	      continue;
	    if (fabs(u-Uhit[k0])>tol) continue;
	  }
	}
	iflh[nFound] = -1;
      }
#endif
      else if (mode>=-1) {     // ***** Determine whether ambiguity => "k_alt"
	int k_alt = -1, k1 = k0+1, k2 = k0-1;
	float delta1 = tol, delta2 = 0.;
	if (k1<lh[ipl] && (delta1 = Uhit[k1]-u)<tol)
	  k_alt = k1;
	if (k2>=fh[ipl] && (delta2 = u-Uhit[k2])<delta1 && delta2<tol)
	  k_alt = k2;

	if (k_alt>=0) {      // ***** Alternative exists
	  if (ifl[k0] &&             // If best is already used...
	      !ifl[k_alt]) {         // ... and alternative is free
	    iflh[nFound] = k0;           // ...Flag hit as ambiguous w/ best
	    k0 = k_alt;                  // ...and use alternative
	  }
	  else iflh[nFound] = k_alt; // ...else flag hit as ambiguous w/ alt
	}
	else iflh[nFound] = -1;
      }

      //if (fabs(y-Uhit[k0])<tol) {  // Why again an "if"?
      hh[nFound++] = k0; // store hit
      if (ifl[k0]>0
	  //#  define FindSpace_ReUSE_HMLO
#ifdef FindSpace_ReUSE_HMLO
	  && typl[ipl]!=42  // HMLO, i.e. crudely segmented hodos
#endif
	  ) nused++;

      //          ***** CUT ON USED HITS
      if (nused>nused_mx) goto next_comb;
    } // End of loop over planes

    //        ********** FINAL CHECK ON #HITS **********
    // It's not assured that it fulfills iPRpar requirements:
    //  - PixelGEMs
    //  - In the case of a fit, the track parameters of the combination
    //   which has fulfilled the requirements have been redefined by the fit.
    // The check is done against "nHsMn_trk" Which is now defined on a track by
    // track basis.
    if (nFound+nPixels<nHsMn_trk) continue;

    //        ********** STORE FOUND HITS **********
    if (mode>=-1 && icsf<0)              // Account for ambiguities
      for (int i = 0; i<nFound; i++) {
	int iref = hh[i], ialt = iflh[i];
	ifl[iref]++;                     // Flag hit as used <-> Increment
	if (ialt>=0) {                   // If alternative...
	  hits[nh++] = -2-ialt;          // => Insert it in hits list (<-1)
	}
	hits[nh++] = iref;
      }
    else
      for (int i = 0; i<nFound; i++) {
	hits[nh++]=hh[i]; ifl[hh[i]]++;  // Flag hit as used <-> Increment
      }
    hits[nh++] = -1;     //end of hit list

    if (++nt>=TConstants_NTtrack_max)  {
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "Too many tracks (=%d) in group %d",nt,igr);
      return -2;
    }

  next_comb:;
  } // End of loop over combinations

  if (TOpt::Print[5])
    printf("  (space) N comb. = %5u   N found tracks = %4u \n", ncomb, nt);

  Ntk = nt;
  return nh;
}
