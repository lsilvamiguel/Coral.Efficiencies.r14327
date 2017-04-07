// $Id: TSetupInit.cc 14094 2015-11-06 15:28:48Z lsilva $

#include <iostream>
#include <list>
#include <algorithm>  
#include "CsInit.h"
#include "CsGeom.h"
#include "CsDetector.h"
#include "CsDriftChamberDetector.h"
#include "CsStrawTubesDetector.h"
#include "CsDWDetector.h"
#include "CsDriftTubeDetector.h"
#include "CsMicroMegaDetector.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TAlgo.h"

#include "cfortran.h"
                       PROTOCCALLSFFUN4(FLOAT,SIMPS,simps,FLOATV,FLOAT,FLOAT,INT)
#define SIMPS(A2,A3,A4,A5)  CCALLSFFUN4(      SIMPS,simps,FLOATV,FLOAT,FLOAT,INT,A2,A3,A4,A5) 

#include "CsHistograms.h"
#include "TH1D.h"

using namespace std;
using namespace CLHEP;

/*!
  Private method for
  initialization of all datamembers of the class TSetup. 
  Called by TSetup constructor.
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TSetupInit":
   i) Added "iPlane" and "associate" member data of class "TPlane".
  ii) Added class "TStation".
 iii) Handle case of pixelGEM/MMs. 
  iv) Determine the hinging points for the bridging over the muFilter.
*/

bool TSetup::Init()
{
  list<CsDetector*> lD = CsGeom::Instance()->getDetectors();
  list<CsDetector*>::iterator id;

  HepMatrix mF, mW;
  HepMatrix rot(3,3,1);

  if(TOpt::dCut[20] <= 0){ 
    cout<<endl;
    cout<<"TSetup::Init() ==> As option dCut[20] <= 0, time measurements will be ignored in tracking"<<endl;
    cout<<endl;
  }

  double totRadLen; int nGEMs;
  for (id = lD.begin(), totRadLen = 0, nGEMs = 0; id!=lD.end(); id++) {

    // ******************** LOOP OVER CORAL DETECTORS ********************

    TDetect d; // create empty detector object

    //                                   ***** DETECTOR HAS TO BE IGNORED? *****
    int siz = sizeof(TOpt::Det2Ignore)/sizeof(string);
    for (int j = 0; j<siz; j++){  // Loop on names of detectors to be ignored
      string &s = TOpt::Det2Ignore[j]; int npos = ((*id)->GetTBName()).find(s);
      if (npos==0 && s.size()>0) {
	setIgnoreIDs.insert((*id)->GetID());
	CsErrLog::msg(elWarning,__FILE__,__LINE__,
 "Det2Ignore: %s is added to 'ignore' list",(*id)->GetTBName().c_str());
	goto next_det;
      }
    }

    mF=(*id)->getRotDRS(); // get detector's frame rotation matrix
    mW=(*id)->getRotWRS(); // get wire ref. system rotation matrix
    
    d.ptrDet= (*id);
    d.iDet  = (*id)->GetID();
    d.iType = (*id)->getType();
    d.name  = (*id)->GetTBName(); 
    for(unsigned int i=0; i < d.name.length(); i++) { // HIGZ does not like undersores
      if(d.name[i] == '_') d.name[i]=' ';
    }
    d.radLen= (*id)->getRdLen();
    if (TOpt::ReMode[8]==1) // Set all rad. lengths to very big value
      d.radLen = 1.E10;
    d.pitch = (*id)->getWirP()/10.;
    d.nwires= (*id)->getNWir();
    d.uorig = (*id)->getWirD()/10.;
    d.x[0]  = (*id)->getZcm ()/10.;
    d.x[1]  = (*id)->getXcm ()/10.;
    d.x[2]  = (*id)->getYcm ()/10.;
    d.siz[0]= (*id)->getZsiz()/10.;
    d.siz[1]= (*id)->getXsiz()/10.;
    d.siz[2]= (*id)->getYsiz()/10.;
    if (d.Name.find("HI")==0 ||
      // - "HI05" plays an important role in the muID: if not attached, the track
      //  will never be ascribed the %X0 associated to the thickness of the
      //  material in front it (so-called "MF3"). On the other hand, its range
      //  is larger than #wires*pitch, because the hodo slabs are larger than
      //  pitch (they overlap). => Exceptionally, set the range = size.
      //   Let's also do the same for "HI04". This time, a missing "HI04" would
      //  not prevent the track from being mu-ID'd. Yet "PID_doMuonID" requires
      //  it for the so-called scattered-mu-ID, which the vertex package requires
      //  to give precedence to the track in the vertex search. In so far as the
      //  inner trigger is specially important for our open charm quest, and
      //  since the range vs. #wires*pitch also holds for HI04, let's make
      //  an other exception and maximize scattered-mu-ID efficiency. (Note: this
      //  is of course at the expense of a loss of purity, which I haven't
      //  estimated...)
	((d.Name.find("DR")==0 || d.Name.find("MA")==0) && d.Name[4]=='X')) {
      // - "DR" and "MA" have a special substructure, and a CsDetector's pitch
      //  and #wires that do not give back the range they span along the u axis.
      d.range = d.siz[1]; d.uedge = -d.siz[1]/2;
    }
    else if (d.Name.find("DR")==0 || d.Name.find("MA")==0) {
      d.range = d.siz[2]; d.uedge = -d.siz[2]/2;
    }
    else {
      d.range = fabs(d.pitch * d.nwires); d.uedge = d.uorig-d.pitch/2.;
    }
    d.hsizY = d.siz[1]/2; d.hsizZ = d.siz[2]/2;

    //            ********** DETECTOR'S TIME RESOLUTION **********
    // - "TResol" is used as safeguard against over-optimistic detector classes:
    //  Time uncertainties better than "TResol" are forced = "TResol".

    d.tresol= -1.; // Default: No measurement of time (e.g. drift detector)
    // Else time resolutions are built-in, infra. (Because det.dat's are not
    // trusted (yet...) on this point.)
    if (TOpt::dCut[20]>0 ||    // When using time explictly requested...
	TOpt::ReMode[11]==2) { // ... or in TrafDic
      if      (d.Name.find("HO")==0) d.tresol =  2.;  // ns
      else if (d.Name.find("HI")==0) d.tresol =  1.5; // ns
      else if (d.Name.find("H" )==0) d.tresol =  1.5; // ns
      else if (d.Name.find("FI")==0) {       // ***** SPECIAL CASE of FI's *****
	const char *stationNum = d.Name.c_str()+2; char *end, **endptr = &end;
	int stNum = strtol(stationNum,endptr,10);
	if (*endptr!=stationNum+2*sizeof(char))
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Detector's name \"%s\" doesn't follow TB scheme\n",d.Name.c_str());
	// Scifi resolution:
	//  - According to Rainer: .5 (J) and .4 (G) ns (while .35 is in MC)
	//  += .1  to be on the safe side (not to change so much w.r.t. previous
	//                                 value, == 1 ns for everybody)
	//  - Not sure that FI15 is of Japanese type(?) 
	if (stNum<5 || stNum==15)    d.tresol = .6; // ns
	else                         d.tresol = .5; // ns
      }
      else if (d.Name.find("SI")==0) {   // Case of SI...
	// ...The best acceptabe time uncertainty is difficult to asses given
	// the fluctuations one observes on the performances of the SI timing
	// in RD, w/, in particular, large offsets.  We decide to have "TResol"
	// set by option ("dCut[75]").
	//  => Let's check that this option has a reasonable value.
	if (TOpt::dCut[75] && TOpt::dCut[75]<1.5)
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Min. acceptable SI time uncertainty \"TraF dCut[75] %f ns\": overoptimistic",
			TOpt::dCut[75]);
	/* */                        d.tresol = TOpt::dCut[75]; // ns
      }
      else if (d.Name.find("MM")==0) { // Time resol. is specific to each plane
	const CsMicroMegaDetector *mM =
	  dynamic_cast<const CsMicroMegaDetector*>(d.PtrDet());
	d.tresol = mM->getTDCResol();
      }
      else if (d.Name.find("MP")==0) d.tresol = 12.; // ns
      else if (d.Name.find("GM")==0) d.tresol = 12.; // ns
      else if (d.Name.find("GP")==0) d.tresol = 12.; // ns
      else if (d.Name.find("PA")==0) d.tresol = 20.; // ns
      else if (d.Name.find("PB")==0) d.tresol = 20.; // ns
      else if (d.Name.find("PS")==0) d.tresol = 20.; // ns
    }

    //            ********** DETECTOR's TYPE **********
    // Overwrite what's in det.dat for some of the TB's.
    // Create a type == 12 for SDC's to distinguish them from Straws
    if (d.Name.find("ST")==0) {
      if (d.Name.find("ST04")==0)    d.iType = 17;   // ST04
      // In 2007/02, I decide to revisit this singling out some of the STs (
      // letting aside the case of ST04, which is specific). Type 13 differs
      // from the rest in the way 2 planes from 2 sub-modules of a same module
      // are allowed to define a segment in a given projection. Why it was
      // decided to grant this possibility to the vertical and outer cases is no
      // longer clear. For the latter, it may be that, in 2002:4, DCs are out of
      // the job and we are left w/ only one module of straws. And this is also
      // the case for any coord at large distance, and specially for the Y coord
      // for which DCs' acceptance is more limited.
      else if (d.Name[7]=='b') {
	if (d.Name[4]!='Y')          d.iType = 11;
	else                         d.iType = 12;
      }
      else {
	if (d.Name[4]!='Y')          d.iType = 13;
	else                         d.iType = 14;
      }
    }
    else if (d.Name.find("DC")==0)   d.iType = 15;
    else if (d.Name.find("DW")==0)   d.iType = 16;
    // Special type for detectors in RICH wall (cf. TEv::PrePattern2)
    else if (d.Name.find("WD")==0 || // (This is DR's name in mc.2006.01.)
	     d.Name.find("DR")==0)                   // DRs  (N.B.: Confusion...
      // ...w/ ST04 is impossible, for the 2 are mutually exclusive.)
                                     d.iType = 17;
    else if (d.Name.find("MB")==0)   d.iType = 18;
    else if (d.Name.find("PS")==0)   d.iType = 2;    // PS
    // Distinct type for pixelGEMs, pixelised and stripped parts...
    else if (d.Name.find("GP")==0) {
      if (d.Name[4]=='P') d.iType = 29; // CsPG pixel core (TBname=="GP..P...")
      //                                   CsPG stripped outskirts (else)
      // (These CsGPs belong to VSAT, as opposed to SAT for standard CsGEMs.)
      else                d.iType = 28;
    }
    else if (d.Name.find("MP")==0) {
      if      (d.Name[4]=='P') d.iType = 32; // CsMP pixel core (TBname=="MP..P...")
      else if (d.Name[4]=='M') d.iType = 32; // CsMP rectangular pixel core (TBname=="MP..M...")
      //                                   CsMP stripped outskirts (else)
      else                     d.iType = 31;
    }
    // A type for MAs, w/ high number, meaning exclude from projection search
    else if (d.Name.find("MA01")==0) d.iType = 39;
    else if (d.Name.find("MA02")==0) d.iType = 40; // Will have lager route
    // Scifis: enforce type code = 22
    else if (d.Name.find("FI")==0)   d.iType = 22;
    // A type for HI04 (It is distinct from all other hodos in that it can
    // be enabled in projection search (cf. TEvPrePattern2)).
    else if (d.Name.find("HI04")==0) d.iType = 43;
    // Special type for cylindrical hodoscopes, so as to be able to easily
    // exclude them (cf. TEv::ImportClusters"): they cannot yet be exploited in
    // TraFDic.
    else if (d.Name[4]=='R')         d.iType = 44;
    // Total. rad length of detectors along beam line (just to print)
    if (d.x[1]==0 && d.x[2]==0 && d.radLen!=0) totRadLen += d.siz[0]/d.radLen;
    if (d.IType==26) nGEMs++;

    //            ********** DETECTOR'S SPACE RESOLUTION **********
    if ((*id)->hasDrift()) {                         // ***** DRIFT DET... *****
      d.kind = 1;  d.resol=(*id)->getVel()/10.*(*id)->getSpSli();
    }
    else {                                // ***** ...ELSE: PITCH/SQRT(12) *****
      // This, again, because det.dat is not trusted yet...
      d.kind = 0; d.resol=fabs(d.pitch)/sqrt(12.);
      if (d.IType==21)             // SIs, for when, e.g., badly aligned:...
	d.resol += TOpt::dCut[84]; // ...position uncertainties correction term
    }

    //            ********** DETECTOR'S ROTATIONS **********
#define TRAFFIC_STRAIGHTEN_HI04X
#ifdef TRAFFIC_STRAIGHTEN_HI04X
    if (d.Name.find("HI04X1")==0) {
      printf("TRAFFIC_STRAIGHTEN_HI04X defined => %s rotation set to unity\n",
	     d.Name.c_str());
      d.cframe = 1; d.sframe = 0; // Frame rotation cos, sin      
      d.ca     = 1; d.sa     = 0; // Wires rotation cos, sin
    }
    else
#endif
      {
	d.cframe =  mF[1][1]; d.sframe = -mF[0][1]; // Frame rotation cos, sin
	d.ca     =  mW[1][1]; d.sa     = -mW[0][1]; // Wires rotation cos, sin
      }
    // center in WRS
    d.xR[0] =  d.x[0];
    d.xR[1] =  d.ca * d.x[1] + d.sa * d.x[2];
    d.xR[2] = -d.sa * d.x[1] + d.ca * d.x[2];
 
    //            ********** DETECTOR'S DEAD ZONE **********
    switch ((*id)->getDZType()) {
    case 0x11:
    case 0x01: 
      d.DZtype = ((*id)->getDZDim()==0.) ?NO :RECTANGULAR; break;
    case 0x15:
    case 0x05: 
      d.DZtype = CIRCULAR; break;
    default:
      d.DZtype = NO;
    }
    if (d.DZtype==RECTANGULAR) d.DZydim  = (*id)->getDZDim()/10;  //  Y 1/2 size
    else                       d.DZydim  = (*id)->getDZDim()/100; // Circular shape => dim = radius squared
    d.DZzdim = (*id)->getDZYdim()/10.;                            // Z 1/2size
    d.DZydrs = (*id)->getDZXdrs()/10.;
    d.DZzdrs = (*id)->getDZYdrs()/10.;
    d.DZymrs = d.cframe*d.DZydrs - d.sframe*d.DZzdrs;  
    d.DZzmrs = d.sframe*d.DZydrs + d.cframe*d.DZzdrs;  
    d.DZymrs += d.x[0];
    d.DZzmrs += d.x[1];
    rot = (*id)->getRotD2DZ();
    d.DZCadrs = rot(1,1);  // (from "asd/lhc++/clhep/manual): "** Note ...
    d.DZSadrs = rot(1,2);  //...that the indexing starts from (1,1) **

    // Check if detector's dead zone had been activated (DZisActive option)
    siz=sizeof(TOpt::DZisActive)/sizeof(string);
    for(int j = 0; j < siz; j++){ // loop over names in "DZisActive" 
      string& s = TOpt::DZisActive[j];
      int npos = d.Name.find(s); // find substring
      if(npos == 0 && s.size()>0){ // found
	d.DZtype = NO;
	cout<<"TSetup::Init() ==> INFO: dead zone in detector "<< d.Name
	    <<" will be considered as active"<<endl;
	break;
      }
    }

    //            ********** DETECTOR'S EMPTY ZONE **********
    // I.e. dead zone when it's empty space.
    // The info (whether massive or empty dead zone) is built-in here.
    // (Note that with the advent of ROOTGeometry, this has become obsolete.)
    if (11<=d.IType && d.IType<=14 /* ST */ || d.IType==17 /* ST04|DR */ ||
	d.IType==18 /* MB */ || 39<=d.IType && d.IType<=40 /* MA */ ||
	41<=d.IType && d.IType<=44 /* Hs, VI */ ||
	d.IType==28 || d.IType==31)// Special case: pixelG/M's stripped piece...
      // ...Their DZ is massive, but is already sensed by tracks when they
      // cross the pixelised central part.
      d.EZtype = d.DZtype;
    else
      d.EZtype = NO;

    vecDetect.push_back(d); // ********** STORE TDetect OBJECT **********

  next_det:;
  } // End of loop over detectors
  cout<<endl;
 
  // Sort vecDetect vector ( see "less" operator)
  sort(vecDetect.begin(), vecDetect.end());
  

  //       *************** FILL vecPlane/vecProj VECTORs ***************
  // (one TPlane per TDetect at the moment)

  int alpha;
  vector<int>::iterator lp;
  vecProj.push_back(0);   // always 1-st proj.
  vecProj.push_back(900); // always 2-d  proj.

  // "SmoothDet" = string used to set a SmoothPos @ the abscissa of TDetect w/
  // matching  name (or average abscissa of group thereof).
  int nSmoothDets = 0; double smoothDetPos = 0;
  map<int,int,less<int> >::const_iterator im;
  for (int i = 0; i<(int)vecDetect.size(); i++) {

    //         *************** LOOP ON TDetects ***************

    TDetect &d = vecDetect[i];
    vecPlane.push_back(TPlane()); TPlane &plane = vecPlane.back();

    //         ***** CHECK IF DETECTOR is OFF in options file BY ID *****
    int siz = sizeof(TOpt::DetOff)/sizeof(int);
    int *p1 = &TOpt::DetOff[0], *p2 = &TOpt::DetOff[siz];
    if (find(p1,p2,d.iDet)!=p2) plane.iFlag = 0;

    //         ***** CHECK IF DETECTOR is OFF in options file (DetNameOff) *****
    siz = sizeof(TOpt::DetNameOff)/sizeof(string);
    for (int j = 0; j<siz; j++) { // Loop over names of switched of detectors
      string &s = TOpt::DetNameOff[j];
      bool found = true;
      if( s.size() > d.Name.size() ) found = false; // Wrong name supplied in opt file
	  // Treat '*' as metechars
      else for( unsigned int name_i=0; name_i < s.size(); name_i++ )
          if( s[name_i] != '*' && s[name_i] != d.Name[name_i] ) {
            found = false; break; // Found difference.
          }
      if(found && s.size()>0){ // found
        plane.iFlag = 0; break;
      }
    }

    //         ***** CHECK IF DETECTOR is to be EXCLUDED FROM FINAL FIT *****
    // By default all detectors (not turned OFF) are included. But one can use
    // the option "Det2Go2Fit" to change this.
    siz = sizeof(TOpt::Det2Go2Fit)/sizeof(string);
    for (int j = 0; j<siz; j++) {   // Loop over names of excluded of detectors
      string &s = TOpt::Det2Go2Fit[j];
      int npos = d.Name.find(s); if (npos==0 && s.size()>0) {
	if (plane.iFlag)  // Exclude turned off detectors
	  plane.iFlag |= 0x2;	// "IFlag&0x2" flags this behaviour.
	break;
      }
    }

    //         ***** SMOOTHING POINT @ DET ABSCISSA... *****
    // ...or abscissa of pair of associated 
    if (!TOpt::SmoothDet.empty() && d.Name.find(TOpt::SmoothDet)==0) {
      smoothDetPos += d.X(0); nSmoothDets++;
    }

    bool isGPP = d.IType==29, isMPM = d.IType == 32;
    // ***** PIXELATED PLANES: CHECK CONSISTENCY W/ OTHER PLANES OF SAME STATION
    if ((isGPP || isMPM) & plane.IFlag) {
      //  I) Case of pixel CsPGP: 0x10|0x20 flag
      // II) Case of pixel CsPMP: 0x40|0x80 flag
      //  In order to be able to easily single out pixel detectors (the actual
      // pixelated ones, tagged w/ a 'P|M' coordinate, as opposed to the
      // stripped pieces of partially pixelated, so-called pixel detectors,
      // which are also CsPG|PM objects, tagged w/ 'X|Y|U|V') in the course of
      // the PR, we flag them w/ a specific "TPlane::iFlag". In so doing, we
      // want to distinguish between XY and UV GPPs, and X|V and Y|U MPMs
      // because PR processes them differently. The ordering of the planes w/in
      // a GP or PM station is fixed so that the index, "jpl", of its stripped
      // twin (i.e. the stripped plane oriented along the pixels' v axis) w.r.t
      // the index, "ipl", of the 'P|M' plane  is:
      //  - 0x10|0x40: "jpl = ipl-2".
      //  - 0x20|0x80: "jpl = ipl+2"
      //  In the PR (in dico fit as far as GPP are concerned, in BackTrack
      // methods for both GPP and MPM), we take advantage of this fixed ordering
      // to retrieve the v coordinate of the pixel piece from its stripped twin.
      //  Also, PR is so built that it expects the TPlanes corresponding to the
      // pixel and stripped pieces of CsPG's to be ordered so that, in a
      // XY+UV station, all 4 stripped pieces come in succession (so as to yield
      // a unique TStation).
      //  => Let's check these orderings.
      // (Note: In coral, the ordering, P1<U<V < Y<X<P2 (resp. M<X|V < Y|U<M),
      // is now enforced in TDetect "<" operator provided P1=U=V < Y=X=P2 (resp.
      // M=X|V < Y|U=M) in input detectors.dat. Before that, we used to rely on
      // the ordering to be artificially encoded in the input detectors.dat,
      // thanks to small differences in Z abscissa. Many detectors.dat's still
      // bear the mark.) 
      string s = d.Name.substr(0,4);         // Station name
      const TDetect *dPrv = i ? &vecDetect[i-1] : 0, *twPrv = dPrv;
      int nDets = (int)vecDetect.size();
      const TDetect *dNxt = i+1<nDets ? &vecDetect[i+1] : 0, *twNxt = dNxt;
      const TDetect *dPPrv = 0, *twPPrv = 0, *dNNxt = 0, *twNNxt = 0;
      if (dPrv) {
	string sPrv = dPrv->Name.substr(0,4); if (sPrv!=s) twPrv = 0;
	else if (dPrv->IType!=28 && dPrv->IType!=31) twPrv =0;
	if (i-1) {
	  dPPrv = &vecDetect[i-2]; string sPPrv = dPPrv->Name.substr(0,4);
	  if (isMPM && sPPrv[1]=='M') sPPrv[1] = 'P';// MM and MP allowed to mix
	  if (sPPrv==s) twPPrv = dPPrv;
	}
      }
      if (dNxt) {
	string sNxt = dNxt->Name.substr(0,4); if (i+2<nDets) {
	  dNNxt = &vecDetect[i+2];
	  string sNNxt = dNNxt->Name.substr(0,4); if (sNxt!=s) twNxt = 0;
	  else if (dNxt->IType!=28 && dNxt->IType!=31) twNxt = 0;
	  if (isMPM && sNNxt[1]=='M') sNNxt[1] = 'P';// MM and MP allowed to mix
	  if (sNNxt==s) twNNxt = dNNxt;
	}
      }
      if (twPrv) {
	// - PixelGP of XY type:  preceding TDetect has type 28 and coord=X
	// - PixelMP of Y|U type: preceding TDetect has type 31 and coord=Y|U
	char coord = dPrv->Name[4];
	if      (isGPP && (coord!='X' || !twPPrv))
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "CsPG Wrong Ordering: %s<%s<%s doesn't follow the P1<U<V < Y<X<P2 scheme",
	    dPPrv?dPPrv->Name.c_str():"??",dPrv->Name.c_str(),d.Name.c_str());
	else if (isMPM && (coord!='Y' && coord!='U'))
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "CsPM Wrong Ordering: %s<%s<%s doesn't follow the M<X|V < Y|U<M scheme",
	    dPPrv?dPPrv->Name.c_str():"??",dPrv->Name.c_str(),d.Name.c_str());
	plane.pixelDet = true;
	if (plane.IFlag) {
	  if (isGPP)       plane.iFlag = 0x10;   // Y   = P - 2
	  else if (twPPrv) plane.iFlag = 0x40;   // X|V = M - 2
	  //else isolated MPM => Let's not flag it.
	}
      }
      else if (twNxt) {
	// - PixelGP of UV type
	// - PixelMP of X or V type
	char coord = dNxt->Name[4];
	if      (isGPP && (coord!='U' || !twNNxt))
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "CsPG Wrong Ordering: %s<%s<%s doesn't follow the P1<U<V < Y<X<P2 scheme",
	    d.Name.c_str(),dNxt->Name.c_str(),dNNxt?dNNxt->Name.c_str():"??");
	else if (isMPM && (coord!='X' && coord!='V')) 
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "CsPM Wrong Ordering: %s<%s<%s doesn't follow the M<X|V < Y|U<M scheme",
	    d.Name.c_str(),dNxt->Name.c_str(),dNNxt?dNNxt->Name.c_str():"??");
	plane.pixelDet = true;
	if (plane.IFlag) {
	  if (isGPP)       plane.iFlag = 0x20;   // V = P + 2
	  else if (twNNxt) plane.iFlag = 0x80;   // Y|U = M + 2
	  //else isolated MPM => Let's not flag it.
	}
      }
      else
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "%s: isolated pixelised %s: #%d while #%d(%s),%d(%s) have types %d and %d",
		      d.Name.c_str(),isGPP?"CsGP":"CsMP",i,
		      i-1,i?vecDetect[i-1].Name.c_str():"?",
		      i+1,i+1<(int)vecDetect.size()?vecDetect[i+1].Name.c_str():"?",
		      i?vecDetect[i-1].IType:-1,
		      i+1<(int)vecDetect.size()?vecDetect[i+1].IType:-1);
    }

    plane.iDetRef = i;    // Store corresponding vecDetect index...
    // ...AND vecPlane index of this TPlane (in view of easily retrieving the
    // index of associated Tplane, cf infra)
    plane.iPlane = vecPlane.size()-1;
    d.iPlane = vecPlane.size()-1;      // ...AND vecPlane index of this TDetect

    im = mapIDvPlane.find(d.iDet);  // map det. ID -> vecPlane index
    if (im!=mapIDvPlane.end())
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "%s: Detector ID(=%d) is not unique !",d.Name.c_str(),d.iDet);
    mapIDvPlane[d.IDet] = d.IPlane; 

    // ********** FILL THE LIST OF ALL EXISTING PROJECTIONS **********
    float phi = atan2(d.sa,d.ca);
    alpha = int(phi*180/M_PI + 0.5 * (phi<0 ? -1 : 1)); // rounding to 1 deg.
    // Patch to avoid SIs misorientations to isolate the detector plane from
    // its counterparts in zone #0 (when in 2003)
    if (d.Name.find("SI")==0) {
      if (fabs(phi*180/M_PI- 0)<=2.0) alpha =  0;
      if (fabs(phi*180/M_PI- 5)<=2.0) alpha =  5;
      if (fabs(phi*180/M_PI-90)<=2.0) alpha = 90;
      if (fabs(phi*180/M_PI+85)<=2.0) alpha =-85;
    }
    // Same thing for FIs
    if (d.Name.find("FI")==0) {
      if (fabs(phi*180/M_PI- 0)<=5.0) alpha =  0;
      if (fabs(phi*180/M_PI-90)<=5.0) alpha = 90;
      if (fabs(phi*180/M_PI-45)<=5.0) alpha = 45;
      if (fabs(phi*180/M_PI+45)<=5.0) alpha =-45;
    }
    // Other patch to avoid isolated detector planes
    if (fabs(double(alpha-90))<=1) alpha = 90;
    if (fabs(double(alpha   ))<=1) alpha = 0;
    if (fabs(double(alpha-45))<=1) alpha = 45;
    if (fabs(double(alpha+45))<=1) alpha =-45;
    if (fabs(double(alpha-30))<=1) alpha = 30;
    if (fabs(double(alpha+30))<=1) alpha =-30;
    // Check GP and GM have either of 0, 90 or +/-45 deg.: this is assumed at
    // some point in the code
    if (d.Name.find("GP")==0 || d.Name.find("GM")==0) {
      if (fabs(alpha- 0)>2 && fabs(alpha-90)>2 &&
	  fabs(alpha-45)>2 && fabs(alpha+45)>2)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "%s (pixel GEM|MM) is at %d deg., while allowed angle are: 0,90,+/-45 deg.\n",
		    d.Name.c_str(),alpha);
    }
    alpha *= 10;
    lp = find(vecProj.begin(), vecProj.end(), alpha);
    if(lp == vecProj.end()) vecProj.push_back(alpha); // add if not yet there
    
    //  ***** STORE SEQUENTIAL NUMBER OF THE PLANE PROJECTION *****
    int ipr; for (ipr = 0; ipr<(int)vecProj.size(); ipr++)
	       if (vecProj[ipr]==alpha) break;
    plane.iProj = ipr;

    if (plane.iFlag!=0) {
      // ********** VARIOUS CHECKS bearing on ACTIVE TPlane's **********

      // Check if original rotation matrix is supported
      mF=d.ptrDet->getRotDRS(); // get detector's frame rotation matrix
#ifdef TRAFFIC_STRAIGHTEN_HI04X
      if (d.Name.find("HI04X1")==0) {
	cout<<"TSetupInit: "<<d.Name<<endl;
	cout<<"WARNING: This detector's been straightened up from track finding"<<endl;
	cout<<"WRS rotation matrix: "<<mW<<endl;
      }
      else
#endif
	{
	  if( mF(3,1) != 0. || mF(3,2) != 0. ||
	      mF(1,3) != 0. || mF(2,3) != 0. ||
	      fabs(mF(3,3)) != 1. ) {
	    cout<<"TSetupInit ==> Detector ID = "<<d.ptrDet->GetID()<<" ("<<d.ptrDet->getName()<<") "
		<<"is not perpendicular to the beam direction (not yet supported)"<<endl;
	    cout<<"WARNING: This detector is excluded from track finding"<<endl;
	    cout<<"DRS rotation matrix: "<<mF<<endl;
	    plane.iFlag = 0; // switch it OFF
	  }
	}

      //                                     ***** MAKE SURE "pitch" BE >0 *****
      if (d.Pitch<0)
	// <0 picth is not supported: we rely on it to be >0 to render easier,
	// and faster, checking that (Y,Z) point is w/ U range (cf. "InActive").
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Detector %s has <0 (=%f) pitch!",d.Name.c_str(),d.Pitch);
    }

  } // End of loop on TDetect's
  // "SmoothDet"...
  if (!TOpt::SmoothDet.empty()) {
    if (nSmoothDets) { // ...found => Set Smooth point @ corresponding abscissa
      TOpt::SmoothPos[0] = smoothDetPos/nSmoothDets;
      CsErrLog::msg(elBasicInfo,__FILE__,__LINE__,
		    "Smooth point set @ average abscissa(=%.1f cm) of detectors matching \"%s\"",
		    smoothDetPos/nSmoothDets,TOpt::SmoothDet.c_str());
    }
    else               // ...not found => Fatal error
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "No detector found matching \"TraF SmoothDet %s\". Check your options file!\"",
		    TOpt::SmoothDet.c_str());
  }

  //          ********** ASSOCIATED PLANE FOR LR AMBIGUITY RAISING **********

  for (int i = 0; i<int(vecPlane.size()); i++) {  // ***** LOOP ON TPlanes *****
    TPlane &pi = vecPlane[i]; TDetect &di = vecDetect[pi.iDetRef];
    CsDetector *ci = di.ptrDet, *associate = NULL;
    if (ci->hasDrift() && (associate = ci->getAssociateDet())) {
      for (int j = 0; j<int(vecPlane.size()); j++) {  // loop over detectors
	TPlane &pj = vecPlane[j]; TDetect &dj = vecDetect[pj.iDetRef];
	if (dj.ptrDet==associate) { pi.Associate = &pj; break; }
      }
    }
  }

  if (TOpt::GEMSpacers.size()) {

    //          ********** GEMs' SPACERS **********

    ifstream f(TOpt::GEMSpacers.c_str(),ios::in);
    if (!f)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Error opening GEMSpacers file \"%s\"",TOpt::GEMSpacers.c_str());
    const int lineSize = 128; char line[lineSize]; string TBname;
    do {
      f.getline(line,lineSize,'\n');
      if (f.eof()) continue;
      if (line[0]=='\0' || line[0]=='/' && line[1]=='/') continue;
      istringstream s(line); s >> TBname; if (!s) continue;
      for (int i = 0; i<int(vecDetect.size()); i++) {
	TDetect &d = vecDetect[i]; if (d.Name.find(TBname)==0) {
	  nGEMs--; 
	  int extrema[6]; for (int j = 0; j<6; j++) {
	    s >> extrema[j]; if (j && extrema[j]<extrema[j-1])
	      CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Error in GEMSpacers file \"%s\" @ \"%s\": %dth extremum < %dth",
		TOpt::GEMSpacers.c_str(),TBname.c_str(),j,j-1);
	  }
	  // We define 2 sets of 3 intervals centered on the 3 spacers,
	  // w/ widthes define by options dCut[72-73].
	  float *bounds= new float[12]; mapID2GEMSpacers[d.IDet] = bounds;
	  for (int set = 0; set<3; set++) {
	    int k = 2*set, l = k+1;
	    float low = extrema[k], up = extrema[l];
	    bounds[k] = d.Uorig+(low-TOpt::dCut[72])*d.Pitch;
	    bounds[l] = d.Uorig+( up+TOpt::dCut[72])*d.Pitch;
	    k += 6; l += 6;
	    bounds[k] = d.Uorig+(low-TOpt::dCut[73])*d.Pitch;
	    bounds[l] = d.Uorig+( up+TOpt::dCut[73])*d.Pitch;
	  }
	  break;
	}
      }
    } while (!f.eof());
    if (nGEMs!=0)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Error in GEMSpacers file \"%s\" only %d GEMs described, out of %d",
	TOpt::GEMSpacers.c_str(),mapID2GEMSpacers.size(),mapID2GEMSpacers.size()+nGEMs);
  }

  //              ********** MAGNET INFO **********
  //   => STORE IT
  //   => DERIVE BEAM CHARGE from it

  CsGeom *g = CsGeom::Instance(); CsField *f = g->getCsField();
  CsMagInfo *mags = f->getMagInfo(); nMags  =  f->getNumOfMags();
  if (NMags!=3) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "# of magnets = %d, while expected = 3",NMags);
  int iCsMag, imag; for (imag=iCsMag = 0; iCsMag<NMags; imag++, iCsMag++) {
    // Loop over magnets in CsMagInfo 
    magType  [imag]     =  mags[iCsMag].mag;
    magCenter[imag][0]  =  mags[iCsMag].zcm/10.;
    magCenter[imag][1]  =  mags[iCsMag].xcm/10.;
    magCenter[imag][2]  =  mags[iCsMag].ycm/10.;
    magScale [imag]     =  mags[iCsMag].fsc;
    magFlag1 [imag]     =  mags[iCsMag].flg1;
    magFlag2 [imag]     = int(mags[iCsMag].flg2);
  }
  // Beam charge: Must be supplied via option.
  // - We here check that the supplied value is compatible w/ the polarity of
  //  the dipole magnets.
  //  (Note: We keep on using option TOpt::iCut, because the polarity retrieved
  //  from "detectors.dat" is no guarantee. When we have access to the polarity
  //  of SM1/2 from the DB, we will require instead the compatibility between
  //  that information and the "detectors.dat", w/o requesting any longer the
  //  user to specify the polarity by option.)
  // - We may still want to allow unorthodox settings, e.g. a low intensity
  //  electron beam w/ <0 fields. Do indeed if beam ID |iCut[30]| = 2.
  Severity severity = abs(TOpt::iCut[30])<2 ? elFatal : elError;
  if (magScale[1]*magScale[2]<0) CsErrLog::msg(severity,__FILE__,__LINE__,
 "Inconsistency in \"detectors.dat\": magnets #1 and #2: scale[1]*scale[2] = %.3f*%.3f <0",magScale[1],magScale[2]);
  if (TOpt::iCut[15]*magScale[1]<0) CsErrLog::msg(severity,__FILE__,__LINE__,
 "Different polarities in \"detectors.dat\" (%c0) and \"TraF iCut[15]\" option (%c0)",magScale[1]>0?'>':'<',TOpt::iCut[15]>0?'>':'<');
  // Check that beam charge has indeed be specified (and != 0) except...
  if (magScale[1] &&// ...for magnets off jobs, where beam charge doesn't matter
      TOpt::iCut[15]==0) CsErrLog::mes(elFatal,
 "Beam charge (option \"TraF iCut[15]\") was not specified");

  if (TOpt::Print[3]>0) {
    cout<<"Magnets: \n";
    for (imag = 0; imag<NMags; imag++){
      cout<<"Mag type = "<<MagType[imag]<<"\t X = "<<MagCenter(imag,0)<<"\t  Scale = "<<MagScale[imag]<<endl;
    }  
    cout<<endl;
  }

  //      ***** MAGNETS: FIELD INTEGRAL, BeamT ANGLE in TRANSVERSE MODE *****
  float fint, xx, a, b, ff[1201]; double xyz[3], bb[3];
  TH1D *hSMs[3]; char hN[] = "hSM1", hT[] = "SM1 Profile";
  CsHistograms::SetCurrentPath("/Traffic");
  for (imag = 0; imag<NMags; imag++) {
    a=b= magCenter[imag][0];
    // Limits of integration? The field integral is used for a 1st estimate of
    // momentum given deviation of downstream w.r.t. upstream segments see
    // "TEv::BridgeSegments[2]", w/ segments being fitted by straight lines.
    // => For SM2, it can be -400 <-> +400 as it had been from the beginning of
    //   TSetup::Init.
    //    For SM1, a good choice should extend from MM02 to ST03/DC05. Therefore,
    //   very crudely -200 <-> +200. At any rate, this prevents the integral
    //   from receiving a strong contribution from the target dipole in the
    //   setups from [2002,2011] where the the distance from SM1 to target is
    //   less than 400.
    //    Let's choose the same for the target (for which the field integral is
    //   also used (see infra) to determine the beamTelescope angle).
    if (imag==2) { a -= 400; b += 400; }
    else         { a -= 200; b += 200; }
    sprintf(hN,"hSM%d",imag); sprintf(hT,"SM%d Profile",imag);
    hSMs[imag] = new TH1D(hN,hT,800,a,b);
    int bin = 1;

    int n(0);
    for (xx = a; xx<=b; xx += 1.0) {
      xyz[0] = xx; xyz[1]=xyz[2] = 0;
      TAlgo::Field(xyz,bb);
      ff[n++]=bb[2];
      hSMs[imag]->SetBinContent(bin++,bb[2]);
    }
    magFieldInt[imag] = SIMPS(ff,a,b,--n);
    if (TOpt::Print[3]>0)
      printf("Magnet %d: X Centre = %6.1fcm     Field integral (Bz in +- %.0f cm) = %8.5f T*m\n",
	     imag,MagCenter(imag,0),b-MagCenter(imag,0),MagFieldInt[imag]/1000);
    if (!imag && MagFlag1[0]>=3) {
      beamTAngle = (MagScale[0]*MagScale[2]>0?-1:1)*
	0.3*fabs(MagFieldInt[0])/1000/TOpt::dCut[4];
      printf("   =>  beamTelescope angle for %.0f GeV = %.2f mrd\n",
	     TOpt::dCut[4],BeamTAngle*1000);
    }
  }

  if (TOpt::Print[3]>0) {  // ***** PRINT-OUT: LIST OF PROJECTIONS *****
    cout<<"List of projections: ";
    for (vector<int>::const_iterator i = vProj().begin(); i!=vProj().end(); i++)
      cout<<*i<<"  "; 
    cout<<endl;
  }

  //              ********** TARGET INFO **********
  //  Target pos.: 2 sources of info: "TraF Target" options and coral.
  //  I) Early COMGEant versions: coral info is not reliable (because derived
  //    not from the specification of the target position but from that of the
  //    target magnet (in "detectors.dat"), which latter might not exist (e.g.
  //    in hadron setup).
  //    => Let's have option's value overwrite coral's.
  // II) COMGeant >= 7.3: Target position (and size) is specified independently
  //    and it is checked that it is indeed read.
  //    => "TraF Target" is deprecated.
  int CGEAversion = atoi(g->getComgeantVers().c_str()+1);
  int CGEArelease = atoi(g->getComgeantVers().c_str()+7);
  if (CGEAversion<7 || CGEAversion==7 && CGEArelease<3) {
    if (TOpt::targetOpt) { // If "TraF Target" has been specified.
      if (fabs(TOpt::Target[0]-g->getTargetCenter()/10)>.01)
	CsErrLog::msg(elError,__FILE__,__LINE__,
 "Target abscissa: TraF option (%fcm) overrides coral's value (%fcm)!",
	  TOpt::Target[0],g->getTargetCenter()/10);
      for (int xyz = 0; xyz<3; xyz++) targetCenter[xyz] = TOpt::Target[xyz];
    }
    else
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "No \"TraF Target\" entry found in options file, while COMGeant,v%d.%d",
	CGEAversion,CGEArelease);
  }
  else {
    if (TOpt::targetOpt)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "COMGeant,v%d.%d>=v7.3 has \"targ\" entries for the target => \"TraF Target\" option is obsolete",
      CGEAversion,CGEArelease);
    targetCenter[0] = g->getTargetCenter()/10;
  }

  //              ********** muWallA INFO **********
  // - Determine the bit pattern of muWallA detectors (MAs and HG02) in the
  //  array of detectors.
  int idet, lastMWA; for (idet = 0, firstMuWallA=lastMuWallA = -1, bMuWallA = 0;
			  idet<(int)vecDetect.size(); idet++) {
    const string &TBname = vecDetect[idet].Name;
    if (TBname.find("MA") && TBname.find("HG02")) continue;
    if (firstMuWallA<0) firstMuWallA = idet; lastMuWallA = idet;
    int jdet = idet-firstMuWallA;
    if (jdet<64) bMuWallA |= static_cast<unsigned long long>(0x1)<<jdet;
  }
  if (firstMuWallA<0) CsErrLog::mes(elFatal,
     "No MuWallA detector (MAs or HG02) found in the setup => Check your \"detectors.dat\"");
  if (lastMuWallA-firstMuWallA>=64) CsErrLog::msg(elFatal,__FILE__,__LINE__,
     "MuWallA detectors (MAs and HG02) span too large an interval(=[%d,%d]) of plane indices",firstMuWallA,lastMWA);



  //     ********** CHECK CORAL'S ZONES DO NOT OVERLAP **********
  // - Which would upset the piece code infra defining TraFFiC groups. (Note
  //  that what really need be avoided are detectors shared by several zones.
  //  Yet, what's checked is the space overlapping, because it's simpler. May
  //  leave some very special probematic cases unaddressed, though.)
  // - Which happens if the user tries and overrides the definition of one zone
  //  (while changing also its name, cf. CsGeom::makeZones), overlooking the
  //  fact that "define zone" is a recursive option.
  list<CsZone*> lZ = g->getZones(); list<CsZone*>::iterator iZ;
  for (iZ = lZ.begin(); iZ!=lZ.end(); iZ++) {
    const CsZone *zi = *iZ; double mini = zi->getZMin(), maxi = zi->getZMax();
    list<CsZone*>::iterator jZ = iZ; for (jZ++; jZ!=lZ.end(); jZ++) {
      const CsZone *zj = *jZ; double minj = zj->getZMin(), maxj = zj->getZMax();
      if (maxj<=mini || minj>=maxi) continue;
      if (minj<mini && maxj>mini+.1 || mini<minj && maxi>minj+.1)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Zones \"%s\"([%.2f,%.2f]) and \"%s\"([%.2f,%.2f]) overlap",
	  zi->getName().c_str(),mini,maxi,zj->getName().c_str(),minj,maxj);
    }
  }
    
  //   *************** ARRANGE DETECTORS into GROUPS ***************
  //   ***************   ACCORDING to CORAL's ZONES  ***************
  // At present: 1 Zone <-> 1 Group
  int igr, ist; for (iZ = lZ.begin(), igr=ist = 0; iZ!=lZ.end(); iZ++, igr++) {

    //     ********** LOOP OVER ZONES **********

    int iplFirst(-1), iplLast(-1), istFirst(-1), istLast(-1), nplanes(0);
    char T = '\0', B = '\0', coord = '\0', cNum = '\0';
    int stNum = -1, nCoords = 0, detType = -1;
    for (int ipl = 0; ipl<int(vecPlane.size()); ipl++) {
      //                      ***** LOOP OVER PLANES ( THEY ARE X ORDERED) *****
      char TPrv = T, BPrv = B, coordPrv = coord, cNumPrv = cNum;
      int stNumPrv = stNum, detTypePrv = detType;
      TPlane &p = vecPlane[ipl]; TDetect &d = vecDetect[p.IDetRef];
      double x = d.X(0); detType = d.IType;
      T = d.Name[0]; B = d.Name[1]; coord = d.Name[4]; cNum = d.Name[5];
      const char *stationNum = d.Name.c_str()+2; char *end, **endptr = &end;
      stNum = strtol(stationNum,endptr,10);
      if (T=='M' && B=='P') {
	B = 'M'; detType = 27; // Mix MMs and MPs in a same TStation
      }
      if (*endptr!=stationNum+2*sizeof(char))
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Detector's name \"%s\" doesn't follow TB scheme\n",d.Name.c_str());
      if ((*iZ)->getZMin()/10<x && x<=(*iZ)->getZMax()/10) { 
	//                      ***** BUILD VECTOR OF STATIONS: vecStation *****
	if (istFirst==-1) istFirst = ist;
	if (coord!=coordPrv || cNum!=cNumPrv) nCoords++;
	if (T!=TPrv || B!=BPrv ||               // NEW TStation: new TB...
	    // ...or, case of pixelGEM w/ GP..P forming a singleton TStation of
	    // type 29 next to TStations GP..X/Y or GP..U/V (of type 28)...
	    detType!=detTypePrv &&              // ...or new type...
	    (detType<11    || 14<detType     || // ...except special case of ST
	     detTypePrv<11 || 14<detTypePrv) ||
	    (stNum!=stNumPrv ||                 // ...or new station number...
	     T=='S' && B=='T' &&                // ...or special case of straws:
	     // Straws can be problematic, for we want to limit the extent of
	     // a TStation to what corresponds, approximately, to a space point
	     // and therefore, in the straws case, to a sub-module. In addition,
	     // there are 3 coord's per sub-module and 6 TPlane's per coord =>
	     // a full module would contain 36 TPlane's: more than TStation's
	     // data member Pat could handle. The 2 sub-modules of a module
	     // cannot easily be distinguished: they bear the same TB, have the
	     // same number and 2 TPlane's from 2 different sub-modules can be
	     // closer to one another than 2 TPlane's w/in a sub-module. =>
	     // We decide to create a new TStation when the count of coords > 3.
	     // (There's a check that total # of TPlanes is <=18, cf. infra.)
	     // (N.B.: DWs stations do not correspond to a space-point w/ this
	     // scheme: they get only 2 coord's. This because XU and YV have
	     // different station#. But, w/ whatever scheme, we cannot group
	     // them in super-sets w/o any intervening PA.)
	     nCoords>3) && (T!='P' || B!='B') &&
	    (T!='D' || B!='R') || // Exclude DR: we want a unique DR TStation
	    T=='P' && B=='B' && nCoords>3) {    // ...or special case of PBs
	  nCoords = 1;
	  vecStation.push_back(TStation());       // Add element to "vecStation"
	  TStation &s = vecStation.back();
	  s.Type = detType; s.IStation = vecStation.size()-1;
	  if      (detType==11) {               // Scifis
	    s.URange=s.VRange = .25;
	  }
	  else if (detType==26 || detType==28) {// Plain CsGEMs (excluding CsPG)
	    static int igem = 0;
	    static string gem_names;
	    if (detType==26) {
	      gem_names += d.Name;
	      if (igem>=11)
		CsErrLog::mes(elFatal,"More than 11 stations of \"GM\" planes: "+gem_names);
	      else {
		s.URange = TOpt::dCut[28+igem]; s.VRange = TOpt::dCut[igem+39];
	      }
	      igem++;
	    }
	    else if (detType==28) {// PixelGEMs, CsGEM & CsPG parts
	      // Recycle the characteristics of the latest (or, if none yet, the
	      // 1st) enountered standard GEM. (This is a preliminary solution,
	      // until the optimum range is determined from data. It allows to
	      // keep from introducing a new bunch of "dCut's".)
	      int ig = igem ? igem-1 : 0;
	      s.URange = TOpt::dCut[28+ig]; s.VRange = TOpt::dCut[ig+39];
	    }
	  }
	  else {                                // Else
	    s.URange = 10*d.Resol; s.VRange = 3*d.Resol;
	  }
	  istLast = ist; nplanes = 0; ist++;
	} // End new TStation
	TStation &s = vecStation.back(); nplanes++;
	s.IPlanes.push_back(ipl);        // Reference TPlane in host TStation
	iplLast = ipl; if (iplFirst==-1) iplFirst = ipl;
	if (s.Type==15 && nplanes>8) // Case DC: require #planes<8, because...
	  // ...this limitation is taken for granted in "TEv::TracksFit2" 
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "DC TStation #%d (\"%s\",...) has more than 8 coordinates!",
			s.IStation,d.Name.c_str());
      }
    } // End of loop over planes
    if (iplFirst!=-1 && iplLast!=-1) {
      vecIplFirst.push_back(iplFirst); vecIplLast.push_back (iplLast);
      mapGroup2Zone[igr] = *iZ; // group->zone map
      // Fill the plane patterns of all the TStation's in current zone
      for (int jst = istFirst; jst<=istLast; jst++) {
	TStation &s = vecStation[jst];
	if ((int)s.IPlanes.size()>18) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "More than 18 planes in station #%d",jst);
	s.IPat = (s.IPlanes[0]-iplFirst)/32;
	s.JPl  = (s.IPlanes[0]-iplFirst)%32;
	int jpl; for (jpl = 0, s.Pat = 0; jpl<(int)s.IPlanes.size(); jpl++) {
	  s.Pat |= 1<<jpl;
	}
      }
    }
    else
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "CsZone %.3f<Z<%.3f mm contains NO detectors",
	(*iZ)->getZMin(),(*iZ)->getZMax());

  } // End loop over zones

  for (ist = 0; ist<(int)vecStation.size(); ist++) {
    //                 ***** LOOP on TStation's *****
    // - Finalize TStation construction: compute longitudinal abscissa X0.
    // - X-reference: add reference to TStation to all its component TPlane's.
    // - Check consistency w/ reconstruction options:
    //   - Option TOpt::ReMode[17], i.e. GEM amplitude correlation =>
    //     GM and GPXYUV stations, in zones 0x16, must have 4 planes.
    // - Debug printout.
    TStation &s = vecStation[ist]; const TDetect *d = 0;
    int ipl; for (ipl = 0, s.X0 = 0; ipl<(int)s.IPlanes.size(); ipl++) {
      int jpl = s.IPlanes[ipl]; TPlane &p = vecPlane[jpl]; p.Station = &s;
      d = &vecDetect[p.IDetRef]; s.X0 += d->X(0);
      if (TOpt::Print[3]>1)
	printf("%s: ID=%4u Station=%2d X=%8.2f Uorig=%8.2f Pitch,sigma=%6.2f,%5.2f (mm) Angle=%5.1f Proj.=%2u\n",
		 d->Name.c_str(),d->IDet,s.IStation,d->X(0),d->Uorig,
		 10*d->Pitch,10*d->Resol,atan2(d->Sa,d->Ca)*180/M_PI,p.IProj);
    }
    s.X0 /= s.IPlanes.size(); // Average TStation longitudinal abscissa
    if (TOpt::Print[3]>1)
      printf("%s    :         Station=%2d X=%8.2f Type=%d\n",
	     d->Name.substr(0,4).c_str(),s.IStation,s.X0,s.Type);
    if (TOpt::ReMode[17] &&            // GEM amplitude correlation enabled
	(s.Type==26 || s.Type==28)) {  // GM or stripped piece of GP 
      // => Determine reconstruction zone
      int igr, jgr = 0, ipl0 = s.IPlanes[0];
      for (igr = 0, jgr = -1; igr<(int)vecIplFirst.size(); igr++) {
	if (ipl0<vecIplFirst[igr] || vecIplLast[igr]<ipl0) continue;
	jgr = igr; break;
      }
      if (jgr<0)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "Detector %s @ Z = %.4f cm does not belong to any reconstruction zone",
		      d->Name.c_str(),d->X(0));
      else if (jgr && igr<4 && // Exclude zone 0x11, cf. TEv::PrePattern2
	       s.IPlanes.size()<4)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "TStation %s of Type %d, a stripped GEM, has %d planes instead of 4, "
	  "while option \"TraF ReMode[17]\" is enabled",
		      d->Name.substr(0,4).c_str(),s.Type,s.IPlanes.size());
    }
  }
  if (vecIplFirst.size()==0)
    CsErrLog::mes(elFatal,
 "# of Zones found to be == 0! Check their definition in options file.");

  // ********** OVERWRITE DETECTORS ATTRIBUTES in SPECIAL CASES **********

  for (int igr = 0; igr<int(vecIplFirst.size()); igr++) {
    for (int ipl = vecIplFirst[igr]; ipl<=vecIplLast[igr]; ipl++) {
      TPlane &p = vecPlane[ipl];
      TDetect &d = vecDetect[p.iDetRef];
      if (d.Name.find("H")==0) {

	// ***** RESET HODOSCOPE RESOLUTION *****
	// So as to have PR routes match (be at least as large as) their
	// pitch (times some margin to correct pitch for overlap)

	double route = d.Resol*TOpt::dPRpar[igr*10]; // PR: select 1st iteration
#define RESET_HODO_RESOL
#ifdef RESET_HODO_RESOL
	if (route<d.Pitch*.8) {
	  CsErrLog::msg(elWarning,__FILE__,__LINE__,
 "%s: Resolution %.3f -> %.3f *= %.3f = PR Route > Pitch = %.3f",d.Name.c_str(),
	    d.resol,d.resol*d.Pitch/route*.8,TOpt::dPRpar[igr*10],d.Pitch);
	    d.resol *= d.Pitch/route*.8;
	}
#endif
      }
    }
  }
  

  if (TOpt::Print[3]>0) { // ********** PRINT to STDOUT **********
    for(int igr = 0; igr < int(vecIplFirst.size()); igr++){
      int nact=0; map<int,int> mPrj; map<int,int>::iterator imPrj;
      if (TOpt::Print[3]>1) cout<<endl<<endl<<
			      "   --- Detector's group # "<<igr<<" ---"<<endl;
      for (int ipl=vecIplFirst[igr]; ipl <= vecIplLast[igr]; ipl++) {
	if(TOpt::Print[3]>1){ //  more detailed printout
	  cout<<"   TraF. plane # "<<ipl<<" : "<<iPlane2Detect(ipl).Name
	      <<" proj. "<<vecPlane[ipl].IProj<<" id "<<iPlane2Detect(ipl).IDet;
	  if(vecPlane[ipl].IFlag == 0) cout<<" is OFF";
	  cout<<endl;
	} 
	if (vecPlane[ipl].IFlag!=0) { // plane is ON
	  nact++; mPrj[vecPlane[ipl].IProj] += 1;
	}
      }
      cout<<endl<<"Detector group # "<<igr
	  <<"  \t First plane : "<<vecIplFirst[igr]
	  <<"  \t Last plane  : "<<vecIplLast [igr]<<" \t ("<<nact<<" are ON)"<<endl;
      cout<<"   (approx. proj. angle) [N planes]: ";
      for(imPrj=mPrj.begin(); imPrj!=mPrj.end(); imPrj++){
	cout<<vecProj[(*imPrj).first]<<" ["<<(*imPrj).second<<"]   ";
      }
      cout<<endl;
    }
  }

  //   *************** MATERIAL MAPS ***************
  if (TOpt::ReMode[20]>0) { // Use of material maps requested
    CsMaterialMap *map = CsGeom::Instance()->getCsMaterialMap();

    usingROOTGeometry = map->usingROOTGeometry();

    if (!usingROOTGeometry)
      {
	if (!map || map->getNofMaps()==0)
	  CsErrLog::mes(elFatal,
			"Extrapolation with material map requested (cf. \"TraF ReMode[20]\"), "
			"while no material map avalable!");
	vecMaterialBorders = map->getZoneBorders();
	if (vecMaterialBorders.size()%2!=0)
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
			"Inconsistency: size of CsMaterialMap::getZoneBorders() = %d not even",
			vecMaterialBorders.size());
	for (int i = 0; i<(int)vecMaterialBorders.size(); i += 2) {
	  float xi = vecMaterialBorders[i], yi = vecMaterialBorders[i+1];
	  for (int j = i; j<(int)vecMaterialBorders.size(); j += 2) {
	    float xj = vecMaterialBorders[j], yj = vecMaterialBorders[j+1];
	    if (xj<xi) {
	      vecMaterialBorders[i]   = xj; vecMaterialBorders[j]   = xi;
	      vecMaterialBorders[i+1] = yj; vecMaterialBorders[j+1] = yi;
	      xi = xj; yi = yj;
	    }
	  }
	}
	for (int i = 0; i<(int)vecMaterialBorders.size(); i += 2) {
	  vecMaterialBorders[i  ] /= 10.; vecMaterialBorders[i+1] /= 10.;
	  printf("Material map #%d defined in [%7.2f,%7.2f] cm\n",
		 i/2,vecMaterialBorders[i],vecMaterialBorders[i+1]);
	}
      }
  }

  //        ***** PREPARE for BRIDGING OVER FILTER *****
  if (vecIplFirst.size()>2 && (vecIplFirst.size()<4 || 6<vecIplFirst.size()))
    // Although the basic principles of TraFFiC are setup independent, there are
    // here and there some restrictive assumptions, particularly in TraFDic.
    // These were devised to fit the case of a setup w/ 5 zones (viz. 0x10:
    // before target, 0x1: target<->SM1, 0x2: SM1<->SM2, 0x4: SM2<->Wall,
    // 0x8: beyonf wall). The case of a single zone setup is also covered (it
    // has to be, in order to process magnets off alignment runs). In order to
    // warn the putative user that would try to run TraFDic on a fancy setup
    // case and make him aware of these restrictions, let's enforce here a
    // requirement on the # of zones. To the 1- and 5-zone setup cases, which
    // are regularly used, I add the cases of
    // - a single zone spectrometer combined w/ a beamTelescope which might
    //  also be useful,
    // - a 4-zone setup where a single zone spans both beamTelescope zone and
    //  zone 0x1, that may be used in the alignment procedure to align the
    //  beamTelescope w.r.t. the spectrometer,
    // - a 6-zone setup that allow to also cover the BMS zone.
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "%d tracking zones, while TraFDic is only able to deal w/ 1, 2, 5 or 6",vecIplFirst.size());
  if (vecIplFirst.size()==5 || vecIplFirst.size()==6) {
    if (TOpt::ReMode[20]>0 &&
	!CsGeom::Instance()->getCsMaterialMap()->usingROOTGeometry()) {
      // If the use of material maps is requested, check that there's one
      // available describing the ECAL2+HCAL2+MUF2 muFilter.
      double anteFilter = vecDetect[vecIplLast[2]].X(0);
      double postFilter = vecDetect[vecIplFirst[3]].X(0);
      int i, muFilterMap;
      for (i = 0, muFilterMap = -1; i<(int)vecMaterialBorders.size()/2; i++)
	if (anteFilter<vecMaterialBorders[2*i] &&
	    vecMaterialBorders[2*i+1]<postFilter) {
	  muFilterMap = i; break;
	}
      if (muFilterMap<0) {
 printf("Searching for the material map describing ECAL2+HCAL2+MUF2...\n");
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "None found among the %d maps supplied in input, fitting in [%.1f,%.1f] cm\n",
		      vecMaterialBorders.size()/2,anteFilter,postFilter);
      }
    }
    //  The bridging over the muFilter (in order to connect 0x4 and 0x8 track
    // segments) is based on hinging points, at which 0x4 and 0x8 extrapolations
    // are required to meet. (Note: The definition of the muFilter here meant is
    // an extended one, including E/HCAL2). The simple hinging point connection
    // is followed by a more rigorous fitting, using ROOTGeometry/matmaps.)
    //  The hinging points are derived from ROOTGeometry/matmaps. There are 4,
    // corresponding to 4 different regions defined in the yz plane. For each,
    // an helix is extrapoled through the material maps.
    ComputeHinges(muFilterHinges);
    MuFilterHinges = muFilterHinges;
  }
  else MuFilterHinges = 0;

  cout<<endl<<"-----------------------------------------"<<endl<<endl;

  return true;
  
}

void TSetup::ComputeHinges(double *hinges, // NULL pointer or array of 3 doubles
			   THlx   *Hs)     // NULL pointer or array of 9 THlx
  const
{
  //  We determine the ``hinging point'' of the bridging over the muFilter for
  // 4 regions ("muFilter" is here the whole of ECAL2+HCAL2+MUF2):
  //  IV) Beyond the reach of ECAL2, where the absorber consists of the muFilter
  //     proper, MF2, and HCAL2.
  //   I) In the central hole of HCAL2, where it consists of ECAL2.
  //  II) In the central hole of MF2, where it consists of ECAL2+HCAL2.
  // III) In between, where  we have all of ECAL2 + HCAL2 + MF2.
  //  The hinging point is the, approximate, point (its abscissa along Z, in
  // fact) where the upstream and downstream track segment meet on average in
  // the simplified view of the process of multi-scattering through the
  // absorber where the track experiences a single kink.
  //  It's defined as the center of gravity of the thickness of absorber,
  // when weighting by the fraction of radiation length.
  //  If there is no ECAL2, only 2 regions: (I) central and (II) external.
  double anteFilter = vecDetect[vecIplLast [2]].X(0);
  double postFilter = vecDetect[vecIplFirst[3]].X(0);
  double E = // We determine the radiation length by propagating an helix
    // through the absorber. => Need assign it a momentum. The precise value
    // of which is not that important.
    TOpt::dCut[4] ? TOpt::dCut[4] : 160;
  double ys[4]; // Boudaries delimiting the 3 regions
  ys[0] = TOpt::Calo[16];                      // HCAL2 hole center: On the axis
  ys[1] = ys[0]+TOpt::Calo[18];                  // Edge of HCAL2 hole
  ys[2] = TOpt::MuonWall[16]+TOpt::MuonWall[18]; // Edge of MF2 hole
  int nRegions; if (TOpt::Calo[30]>1) {    // !=(double)0 => ECAL2 exists...
    nRegions = 4;
    ys[3] = TOpt::Calo[31]+TOpt::Calo[34];   // ...Outer edge of ECAL2
  }
  else {
    nRegions = 3;
    ys[3] = TOpt::Calo[11]+TOpt::Calo[14];   // ...Outer edge of HCAL2
  }
  //#define TSetup_DEBUG
#ifdef TSetup_DEBUG
  printf("Calo[8] = %.2f, ys = %.2f %.2f %.2f\n",
	 TOpt::Calo[8],ys[0],ys[1],ys[2]);
#endif
  for (int region = 0; region<nRegions; region++) {
    THlx Hfirst, Hlast;
    Hfirst(0) = anteFilter; Hfirst(2)=Hfirst(3)=Hfirst(4) = 0; Hfirst(5) = 1/E;
    switch (region) {
      // Y coord of the helix (it's otherwise in midplane and // to the axis
    case 0: Hfirst(1) = (ys[0]+ys[1])/2; break;
    case 1: Hfirst(1) = (ys[1]+ys[2])/2; break;
    case 2: Hfirst(1) = (ys[2]+ys[3])/2; break;
    case 3: Hfirst(1) = ys[3] + /* some arbitrary margin */ 10;
    }
    if (Hs) Hs[3*region] = Hfirst; THlx Htmp(Hfirst);
    int step, nSteps = int((postFilter-anteFilter)*2); // 5mm steps
    double sXX0, sxXX0; for (step = 0, sXX0=sxXX0 = 0; step<=nSteps; step++) {
      if (step==nSteps) Hlast(0) = postFilter;
      else              Hlast(0) = anteFilter+1+step/2.; // 5mm steps
      Htmp.Extrapolate(Hlast,true); double XX0 = Hlast.RadLenFr();
      sXX0 += XX0; sxXX0 += XX0*(Htmp(0)+Hlast(0))/2;
      Htmp = Hlast;
    }
    if (Hs) { Hs[3*region+1] = Hlast; Hs[3*region+1].SetRadLenFr(sXX0); }
    double hinge = sxXX0/sXX0; if (hinges) hinges[region] = hinge;

    Hlast(0) = hinge; Hfirst.Extrapolate(Hlast,true);
    if (Hs) Hs[3*region+2] = Hlast;
  }
}

void TSetup::Update()
// Update field integral for each magnet which scale factor turns out to have changed
{
  const char *magnets[] = {"Target","SM1","SM2"};

  CsGeom *g = CsGeom::Instance(); CsField *f = g->getCsField();
  CsMagInfo *mags = f->getMagInfo();
  // Note: it was already checked that there are 3 magnets, cf. TSetup::Init.
  int iCsMag, imag; for (imag=iCsMag = 0; iCsMag<NMags; imag++, iCsMag++) {
    if (fabs(magScale[imag]-mags[iCsMag].fsc)>1e-4) {
      float oldScale = magScale[imag]; double oldFint = magFieldInt[imag];
      float fint, a, b, xx, ff[1201]; double xyz[3] = {0,0,0}, bb[3];
      a=b = magCenter[imag][0]; a -= 400; b += 400;
      int n; for (xx = a, n = 0; xx<=b; xx += 1.0) {
	xyz[0] = xx; TAlgo::Field(xyz,bb); ff[n++] = bb[2];
      }
      magScale[imag] = mags[iCsMag].fsc; magFieldInt[imag] = SIMPS(ff,a,b,--n);
      printf("\nTraffic ==> %s's scale %.4f -> %.4f => Update field integral %.4f -> %.4f T*m\n",
	     magnets[imag],oldScale,magScale[imag],oldFint/1000,magFieldInt[imag]/1000);
    }
  }
}
