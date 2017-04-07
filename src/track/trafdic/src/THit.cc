// $Id: THit.cc 14069 2015-09-17 20:44:46Z lsilva $

#ifdef __STRICT_ANSI__
#define __GES_STRICT_ANSI__
#undef __STRICT_ANSI__
#endif
#include <math.h>
#ifdef __GES_STRICT_ANSI__
#define __STRICT_ANSI
#endif

#include <iostream>
#include <iterator>
#include "TSetup.h"
#include "THit.h"
#include "TEv.h"

using namespace std;

/*!
  Print THit's, and associated (if any) THitMC's, info.  
*/

/*
  (Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/THit.cc":
  - Print which-side info for split channels of straws (except ST04), MAs and
  stripped pixelGEMs.
*/

void THit::Print(int ifl) const
{
  const TEv &ev = TEv::Ref();
  const THit &h = *this; int ih = h.IHit, ikin  = h.IKine, ikin2 = h.IKin2;
  const TDetect &d = h.DetRef();

  printf("Cluster # %4u  pl. %3u (%s)  X = %8.2f  ",
	 ih,h.IPlane,d.Name.c_str(),d.X(0));

  if (ev.IsMC()) {      // ********** MC: HIT's ORIGIN (PARTICLE<MOTHER or GHOST
    if (ikin<-1)
      printf(" Mirror hit in a drift detector\n");
    else if (ikin==-1 && ikin2==-1)
      printf(" cluster not associated with MC track");

    else {      
      string sPart, sMother;
      if (ikin>=0) {
	int iv = ev.vKine(ikin).IVtxOrig; if (iv<0) {
	  cout<<"\nTHit::Print: Wrong origin vertex index for KINE track\n";
	  return;
	}
	sPart = ev.vKine(ikin).Name();
	int mother = ev.vVtxMC(iv).IOrig; if (mother<0)
	  sMother = iv==0 ? "prim.vtx": " pileup "; 
	else
	  sMother = ev.vKine(mother).Name();
      }
      else {
	sPart = "part?"; sMother = "vertex?";
      }
      cout<<"MC tr. # "<<ikin<<" ("<<sPart<<")   from "<<sMother<<"  IOrig = "<<h.IOrig;
      if (ikin2>=0) {
	int iv = ev.vKine(ikin2).IVtxOrig; if (iv<0) {
	  cout<<"\nTHit::Print: Wrong origin vertex index for KINE track\n";
	  return;
	}
	sPart = ev.vKine(ikin2).Name();
	int mother = ev.vVtxMC(iv).IOrig; if (mother<0)
	  sMother = iv==0 ? "prim.vtx": " pileup "; 
	else 
	  sMother = ev.vKine(mother).Name(); 
	cout<<", tr. # "<<ikin2<<" ("<<sPart<<")   from "<<sMother;
      }
    }
  } // end of isMC()
  cout<<endl;

  if (!h.sTrackID().empty()) {   // ********** ASSOCIATED RECONSTRUCTED TRACK(S)
    cout<<"Hit belongs to reconstructed track(s) with ID = ";
    copy(h.sTrackID().begin(), h.sTrackID().end(), ostream_iterator<int>(cout, "  "));
    cout<<endl;
  }

  //                                            ********** SPATIAL and TIME INFO
  char whichSide = '\0';                                 // ***** WHICH-SIDE?...
  if (d.IType==11 || d.IType==12 ||
      d.IType==39 || d.IType==40) {            // ...in STs (except ST04) or MAs
    double side = h.ptrCl->getDigitsList().front()->getData()[1];
    if      (side>.5)  whichSide = d.Name[4]=='Y'? '-' : '+';
    else if (side<-.5) whichSide = d.Name[4]=='Y'? '+' : '-';
  }
  else if (d.IType==28 || d.IType==31) {          // ...in STRIPPED PixelGEM/MMs
    const vector<double> &info = h.ptrCl->getAllAnalogDataErrors();
    if (info.size()<4) whichSide = '?';
    else {
      // Note that only those strips that cross the pixelised central piece
      // yield an unambiguous info per se. For all others, the cluster might sit
      // close to the border between the 2 sides. The unambiguous strips span
      // the size of the pixelised piece, viz. 32 (raws of pixels) * .1 (pixel
      // pitch), are 32 / .04 (strip pitch) = 80 out of 256. They correspond to
      // channels [88,167]. Note also that, in any case, part of the <0 side may
      // correspond to a >0 V coordinate because of the detector's offset w.r.t.
      // spectrometer axis.
      if      (info[3]> 0.1)     whichSide = '+';
      else if (info[3]<-0.1)     whichSide = '-';
      else /* if == (double)0 */ whichSide = ' ';
    }
  }

  printf("U = %8.2f V = %c%8.2f SigU = %9.2e SigV = %9.2e  ",
	 h.U,whichSide,h.V,h.SigU,h.SigV);
  if   (h.SigT>0) printf("T = %8.2f+/-%5.2f (ns)\n",h.Time,h.SigT);
  else {
    // TraFDic retains time info only for time-measuring, as opposed to drift,
    // detectors, cf. TDetect::TResol in TSetup::Init and TEv::ImportClusters.
    // => Let's get the info from CsCluster
    double t; h.PtrClus()->getTime(t); printf("T = %8.2f+/-??.?? (ns)\n",t);
  }

  if (11<=d.IType && d.IType<=18)     // ********** DRIFT HITS: LR proba, Status
    printf("Mirror %d LR %.1f%% Status %d\n",h.Mirror,h.ptrCl->getLRProb()*100,h.Status);

  int ndig=0;  set<TDigit>::const_iterator id;         // ********** DIGITS INFO
  for(id = h.sDigits().begin(); id != h.sDigits().end(); id++){
    const TDigit& d = (*id); ndig++;
    cout<<"Dig. "<<ndig<<":   Wire # "<<d.IWire<<"  ADC or TDC info : ";
    for(int j=0; j < int(d.vDigInfo.size()); j++){
      cout<<d.vDigInfo[j]<<"  ";
    }
    cout<<endl;
  }
  if(h.PtrClusMirr() != NULL) { // DC hit
    cout<<"Drift distance = "<<10000.*h.DeltaR<<" mic.";
    if(h.DRsignMC != 0) cout<<" sign known from MC is "<<h.DRsignMC;
    cout<<endl;
  }


  if(ifl > 1) {
    cout<<" Rotated hit paramaters: "<<endl;
    cout<<" Y = "<<h.y<<"   Z = "<<h.z<<endl;
    cout<<" G0 = "<<h.g0<<endl;
    cout<<" G1 = "<<h.g1<<"  G2 = "<<h.g2<<endl;
  }

  cout<<endl;
}


//! returns reference to detector for "this" hit
const TDetect& THit::DetRef() const
{
  const TSetup& setup = TSetup::Ref();
  if(IPlane < 0 || IPlane >= int(setup.vPlane().size())){
    cout<<"THit::DetRef ==> IPlane = "<<IPlane<<endl; 
    assert(false); 
  }
  if(setup.vPlane(IPlane).IDetRef < 0 ||
     setup.vPlane(IPlane).IDetRef >= int(setup.vDetect().size())){
    cout<<"THit::DetRef ==> IDetRef = "<<setup.vPlane(IPlane).IDetRef<<endl; 
    assert(false); 
  }
  return(setup.vDetect(setup.vPlane(IPlane).IDetRef));
};



// Rotate measurement to MRS

void THit::Rotate()
{
  const TDetect& d = DetRef();

  if( ::std::isnan(u)    ||  ::std::isnan(v)    ||
      ::std::isnan(sigu) ||  ::std::isnan(sigv) ||
      ::std::isnan(d.Ca) ||  ::std::isnan(d.Sa)) {
    cout<<"THit::Rotate ==> NaN on input !"<<endl;
    cout<<"                 Bad hit info:"<<endl;
    this->Print(1);
    assert(false);
  }
  if( ::std::isinf(u)    ||  ::std::isinf(v)    ||
      ::std::isinf(sigu) ||  ::std::isinf(sigv) ||
      ::std::isinf(d.Ca) ||  ::std::isinf(d.Sa)) {
    cout<<"THit::Rotate ==> Infinity on input!"<<endl;
    cout<<"                 Bad hit info:"<<endl;
    this->Print(1);
    assert(false);
  }

  if(sigu < 0.0001){
    cout<<"THit::Rotate ==> wrong Sig U = "<<sigu
	<<"  from "<<d.Name<<endl;
    assert(false);
  }
  if(sigu > 1.e6){
    cout<<"THit::Rotate ==> non-realistic Sig U = "<<sigu
	<<"  from "<<d.Name<<endl;
  }
  double wu = 1./(sigu*sigu);

  if(sigv < 0.0001){
    cout<<"THit::Rotate ==> wrong Sig V = "<<sigv
	<<"  from "<<d.Name<<endl;
    assert(false);
  }
  if(sigv > 1.e6){
    cout<<"THit::Rotate ==> non-realistic Sig V = "<<sigv
	<<"  from "<<d.Name<<endl;
  }
  double wv = 1./(sigv*sigv);

  // rotation back to MRS (i.e. by transpose R matrix)

  y = d.Ca*u - d.Sa*v;
  z = d.Sa*u + d.Ca*v;

  g0=wu*d.Ca*d.Ca + wv*d.Sa*d.Sa;
  g1=wu*d.Ca*d.Sa - wv*d.Ca*d.Sa;
  g2=wv*d.Ca*d.Ca + wu*d.Sa*d.Sa;

  if(g0 <= 0  || g2 <= 0){
    cout<<"THit::Rotate ==> wrong weight matrix after rotation :"<<endl;
    this->Print(2);
    assert(false);
  }
  if(sqrt(1./g0) < 0.0001 || sqrt(1./g2) < 0.0001){
    cout<<"THit::Rotate ==> non-realistic weight matrix after rotation :"<<endl;
    this->Print(2);
    assert(false);
  }

};


/*!

  Use drift information, according to flag
  State of the object is changed.
  \param idf - 0 - No deltaR used (or not DC); -1 - minus deltaR; +1 - plus deltaR
  2 - special (not yet done)
 
*/
void THit::UseDrift(int idf)

{ 
  if(idf == 0) return; // not a drift detector. Exit.

  if(ptrClMirr == NULL && idf != 0){
    const TSetup& setup = TSetup::Ref();
    const TDetect& d = setup.iPlane2Detect(IPlane);
    cout<<"THit::UseDrift ==> request to use drift info. for non-drift detectors' hit"<<endl
	<<"Hit # "<<IHit<<" on the plane "<<IPlane<<" ("<<d.Name<<")"
	<<" with DeltaR = "<<DeltaR<<" flag "<<idf
	<<endl<<endl;
    return;
  }

  // take more precise  hit coordinate

  if(idf == -1) u -= deltaR; 
  if(idf == +1) u += deltaR; 

  // change error
  if(idf == +1 || idf == -1 ) {
    sigu = siguDC;
  } else {
    sigu = 2*deltaR; // special (to be improved)
  }

  // get rotated to MRS coordinates and weight matrix
  Rotate();
}
