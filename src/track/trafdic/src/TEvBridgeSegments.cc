
#include <iostream>
#include "TOpt.h"
#include "TDisplay.h"
#include "TSetup.h"
#include "TEv.h"
#include "CsHistograms.h"

using namespace std;

/*! 
  \brief Pair of track segments

  Local class of TEv::BridgeSegments function
  used for sorting of track segment pair candidates

 */
class TPairToBridge {
public:
  double Chi2;                //!< track candidate Chi2
  double Pinv;                //!< track candidate q/P
  double ErrP;                //!< track candidate Sigma(q/P)^2
  list<TTrack>::iterator iL;  //!< "upstream"   track segment iterator
  list<TTrack>::iterator iR;  //!< "downstream" track segment iterator
  //! "Less" operator
  bool operator < (const TPairToBridge& tp) const { return (Chi2 < tp.Chi2);} 
};

/*!

  Bridge track segments through the SM1 and SM2 

*/

void TEv::BridgeSegments()
{
  if(TOpt::ReMode[3] > 0) return; // bridging is OFF

  const TSetup& setup     = TSetup::  Ref();
  TDisplay& display = TDisplay::Ref();

  bool hist(false); if(TOpt::Hist[4] > 0) hist=true;
  int print = TOpt::Print[4]; 

  static bool first=true;
  static CsHist1D* h[50];
  bool ret;

  if(first && hist){
    first=false;
    CsHistograms::SetCurrentPath("/Traffic/BridgeSegment");    
    /*
    HBOOK1(TOpt::Hoffs+100,  "Field integral mag 1 (KGauss*cm)",  100,  -2000.,      0.,   0.);
    HBOOK1(TOpt::Hoffs+200,  "Field integral mag 2 (KGauss*cm)",  100,   3500.,      5500.,   0.);
    for(int ii = 1; ii <=2; ii++){  // mag 1,2
      HBOOK1(TOpt::Hoffs+100*ii+1,"(1/P)-(1/P)mc",               100,  -0.0025,   0.0025,   0.);
      HBOOK1(TOpt::Hoffs+100*ii+2,"((1/P)-(1/P)mc)/Sigma",       100,  -50,   50.,   0.);
      HBOOK2(TOpt::Hoffs+100*ii+4,"Delta dip.(true pair)",        50,   0., 100., 50,  -0.05 , 0.05 ,  0.);
      HBOOK2(TOpt::Hoffs+100*ii+5,"Delta Z.(true pair)",          50,   0., 100., 50,  -10.0,  10.0,     0.);
      HBOOK1(TOpt::Hoffs+100*ii+6,"Chi2 (true pair)",            100,   0., 100.,     0.);
      HBOOK1(TOpt::Hoffs+100*ii+7,"Delta dip.(false pair)",      100,  -0.025 , 0.025 ,  0.);
      HBOOK1(TOpt::Hoffs+100*ii+8,"Delta Z.(false pair)",        100,  -5.0, 5.0,   0.);
      HBOOK1(TOpt::Hoffs+100*ii+10,"P of all true pairs.",       310,  0.,   310.,     0.);
      HBOOK1(TOpt::Hoffs+100*ii+11,"P after Z geom. cut.  ",     310,  0.,   310.,     0.);
      HBOOK1(TOpt::Hoffs+100*ii+12,"P after P estim. cut. ",     310,  0.,   310.,     0.);
      HBOOK1(TOpt::Hoffs+100*ii+13,"P after prelim Chi2 cut",    310,  0.,   310.,     0.);
      HBOOK1(TOpt::Hoffs+100*ii+14,"P after last Chi2 cut.",     310,  0.,   310.,     0.);
      HBOOK1(TOpt::Hoffs+100*ii+15,"P of all wrong pairs.",      310,  0.,   310.,     0.);
      HBOOK1(TOpt::Hoffs+100*ii+16,"P of bridged wrong pairs.",  310,  0.,   310.,     0.);
      
      HBOOK1(TOpt::Hoffs+100*ii+17,"Chi2/ndf 1-st iter. true",       100,   0.,    1000.,     0.);
      HBOOK1(TOpt::Hoffs+100*ii+18,"Chi2/ndf 1-st iter. false",      100,   0.,    1000.,     0.);
      
      HBOOK1(TOpt::Hoffs+100*ii+30,"N iterations",               20,    0.0,   20.,   0.);
    }
    */
    h[0]  = new CsHist1D("sm1h0", "SM1. All pairs Delta Y (cm)",     100, -10, 10);
    h[1]  = new CsHist1D("sm1h1", "SM1. All pairs Delta Z (cm)",     100, -10, 10);
    h[2]  = new CsHist1D("sm1h2", "SM1. All pairs Delta Azi (mrad)", 100, -10, 10);
    h[3]  = new CsHist1D("sm1h3", "SM1. All pairs Delta Dip (mrad)", 100, -10, 10);
    h[4]  = new CsHist1D("sm1h4", "SM1. All pairs Delta Yp (mrad)",  100, -10, 10);
    h[5]  = new CsHist1D("sm1h5", "SM1. All pairs Delta Zp (mrad)",  100, -10, 10);
    h[6]  = new CsHist1D("sm1h6", "SM1. All pairs Chi2/ndf",         100,   0, 500);
    h[7]  = new CsHist1D("sm1h7", "SM1. All pairs Chi2/ndf (zoom)",  100,   0, 50);

    h[10] = new CsHist1D("sm2h0", "SM2. All pairs Delta Y (cm)",     100, -2.5, 2.5);
    h[11] = new CsHist1D("sm2h1", "SM2. All pairs Delta Z (cm)",     100, -5, 5);
    h[12] = new CsHist1D("sm2h2", "SM2. All pairs Delta Azi (mrad)", 100, -2.5, 2.5);
    h[13] = new CsHist1D("sm2h3", "SM2. All pairs Delta Dip (mrad)", 100, -5, 5);
    h[14] = new CsHist1D("sm2h4", "SM2. All pairs Delta Yp (mrad)",  100, -2.5, 2.5);
    h[15] = new CsHist1D("sm2h5", "SM2. All pairs Delta Zp (mrad)",  100, -5, 5);
    h[16] = new CsHist1D("sm2h6", "SM2. All pairs Chi2/ndf",         100,   0, 500);
    h[17] = new CsHist1D("sm2h7", "SM2. All pairs Chi2/ndf (zoom)",  100,   0, 50);
    CsHistograms::SetCurrentPath("/");

  }
  
  TPairToBridge tp;
  list<TPairToBridge> lTrackPair1, lTrackPair2;
  list<TPairToBridge>::iterator itp1,itp2;

  for(int imag = 2; imag > 0; imag--){  // loop over magnets tracks are bridged through

    TTrack tr,tl,tx;

    double x0(0), fldint(0);
    bool good_pair;
    bool bad_pair ;
    float pimc(0), pinv, pinv_prev;
    double chi2, chi2_prev;
    int iter;

    list<TTrack>::iterator il,ir;
    for(il=listTrack.begin(); il != listTrack.end(); il++) {  // loop over left track pieces
      if(imag == 1 && (*il).Type != 1) continue;
      if(imag == 2 && (*il).Type != 2 && (*il).Type != 3) continue;
      tl = (*il); // copy left track

      for(ir=listTrack.begin(); ir != listTrack.end(); ir++) {  // loop over right track pieces
	if(il == ir) continue; // the same
	if(imag == 1 && (*ir).Type != 2 && (*ir).Type != 6) continue;
	if(imag == 2 && (*ir).Type != 4) continue;
	tr = (*ir); // copy right track

	good_pair=false; // = true for correct pair of well reconstructed tracks
	bad_pair =false; // = true for wrong   pair of well reconstructed tracks
	if(tl.IKine >=0 && tr.IKine >=0 && tl.IKine == tr.IKine) {
	  good_pair = true;
	  pimc = float(vecKine[tl.IKine].Pinv()); 
	} 

	if(TOpt::ReMode[6] > 0 && (!good_pair) ) continue; // "ideal" bridging

	if(tl.IKine >=0 && tr.IKine >=0 && tl.IKine != tr.IKine) {
	  bad_pair  = true;
	  pimc = float(vecKine[tl.IKine].Pinv()); // take left track piece mom. 
	} 

	if(hist){
	  if(good_pair) {
	    ///HF1(TOpt::Hoffs+imag*100+10,(1./fabs(pimc)),1.);
	    float finteg=(tl.Hfirst(3)-tr.Hlast(3))/(1.E-3 * 0.3 * vecKine[tl.IKine].Pinv());
	    ///HF1(TOpt::Hoffs+imag*100,finteg, 1.);
	  }
	  if(bad_pair) {
	    ///HF1(TOpt::Hoffs+imag*100+15,(1./fabs(pimc)),1.);
	  }
	}

	if(imag == 1) {
	  x0 =   setup.MagCenter(setup.NMags-2, 0);
	  fldint=setup.MagFieldInt[setup.NMags-2];
	}
	if (imag == 2) {
	  x0 =   setup.MagCenter(setup.NMags-1, 0);
	  fldint=setup.MagFieldInt[setup.NMags-1]; 
	}
	
	// Rough track par. comparisson
	double dazi= tr.Hfirst.azi() - tl.Hlast.azi();
	double ddip= tr.Hfirst.dip() - tl.Hlast.dip();
	double dyp= tr.Hfirst(3) - tl.Hlast(3);
	double dzp= tr.Hfirst(4) - tl.Hlast(4);
	double dy  = (tr.Hfirst(1) + (x0-tr.Hfirst(0))*tr.Hfirst(3)) 
	           - (tl.Hlast (1) + (x0-tl.Hlast (0))*tl.Hlast (3));
	double dz  = (tr.Hfirst(2) + (x0-tr.Hfirst(0))*tan(tr.Hfirst.dip())) 
	           - (tl.Hlast (2) + (x0-tl.Hlast (0))*tan(tl.Hlast.dip()));

	if(hist){
	  h[(imag-1)*10 + 0]->Fill(dy);
	  h[(imag-1)*10 + 1]->Fill(dz);
	  h[(imag-1)*10 + 2]->Fill(dazi*1000.);
	  h[(imag-1)*10 + 3]->Fill(ddip*1000.);
	  h[(imag-1)*10 + 4]->Fill(dyp*1000.);
	  h[(imag-1)*10 + 5]->Fill(dzp*1000.);
	}
	
	if(print > 0){
	  cout <<"\n\n BridegeSegments ==> "<<tl.Id<<"("<<tl.IKine<<")   "<<tr.Id<<"("<<tr.IKine<<")";
	  if(good_pair) cout<<"  Pmc = "<<1./pimc;
	  cout<<"  mag = "<<imag<<" FldInt = "<<fldint<<" x0 = "<<x0<<endl;
	  cout<<"Delat Azi (mrad) "<<dazi*1000.<<"   Delta Dip (mrad) = "<<ddip*1000.
	      <<"    Delta Z = "<<dz<<endl;
	}
	if(hist){
	  if(good_pair){ // correct bridge
	    ///HF2(TOpt::Hoffs+100*imag+4,(1./fabs(pimc)), float(fabs(ddip)),1.);
	    ///HF2(TOpt::Hoffs+100*imag+5,(1./fabs(pimc)), float(fabs(dz)),  1.);
	  }
	  if(bad_pair){ // incorrect bridge
	    ///HF1(TOpt::Hoffs+100*imag+7,float(fabs(ddip)),1.);
	    ///HF1(TOpt::Hoffs+100*imag+8,float(fabs(dz)),1.);
	  }
	}

	// Cuts
	if(TOpt::ReMode[6] == 0 ){
	  if(fabs(ddip*1000.) > 5.0)   continue;
	  if(fabs(dz)         > 5.0)   continue;
	}

	if(hist){
	  if(good_pair) {
	    ///HF1(TOpt::Hoffs+imag*100+11,(1./fabs(pimc)),1.);
	  }
	}

	//Estimate momentum

	pinv=(tl.Hlast(3)-tr.Hfirst(3))/(1.E-3 * 0.3 * fldint);
	
	if(print > 0) cout<<"P estim. =  "<<1./pinv<<endl;
	
	// Cut
	if(TOpt::ReMode[6] == 0 ){
	  if(fabs(1./pinv) > 500. || fabs(1./pinv) < 0.1) goto next_pair;
	}

	if(hist){
	  if(good_pair) {
	    ///HF1(TOpt::Hoffs+imag*100+12,(1./fabs(pimc)),1.);
	  }
	}

	// Build merged track
	tx=tl; tx.Append(tr,"tmp"); tx.Chi2tot = 1e6;

	// Initial values
	chi2=1.E10; tx.Hfirst(5) = pinv; tx.Hfirst(5,5) = 1.E-4;

	// Super ideal bridge: get init value from MC
	if (TOpt::ReMode[6] > 1) tx.Hfirst(5) = pimc;


	iter = 0;
	while(iter < 15) { //  <---------- iteration loop

	  pinv_prev = tx.Hfirst(5); chi2_prev = chi2;

	  tx.Hlast(5) = tx.Hfirst(5);  tx.Hlast(5,5) = tx.Hfirst(5,5);  // set momentum and error

          //display.SetRkutraj(true);
	  ret = tx.QuickKF(-1,1);
	  //display.SetRkutraj(false);
	  if( ! ret ) goto next_pair;  // let's try to fit this candidate

	  chi2=tx.Chi2tot/tx.NHits;

	  if(print > 0){
	    cout<<"Iteration # :  "<<iter<<endl;
	    cout<<"Chi2/Nhits = "<<chi2<<endl;
	  }
	  if(hist && iter == 0){
	    h[(imag-1)*10 + 6]->Fill(chi2);
	    h[(imag-1)*10 + 7]->Fill(chi2);
	    if(good_pair){
	      ///HF1(TOpt::Hoffs+100*imag+17,float(chi2),1.);
	    }
	    if(bad_pair){
	      ///HF1(TOpt::Hoffs+100*imag+18,float(chi2),1.);
	    }
	  }

	  // Cuts
	  if(TOpt::ReMode[6] == 0 ) {
	    if(fabs(1./tx.Hfirst(5)) < 0.1) {
	      if(print > 1) cout<<"Too small momentum: "<<1./tx.Hfirst(5)<<endl;
	      goto next_pair;
	    }
	    if(fabs(1./tx.Hfirst(5)) > 500.) {
	      if(print > 1) cout<<"Too big momentum: "<<1./tx.Hfirst(5)<<endl;
	      goto next_pair;
	    }
	    if(chi2 > TOpt::dCut[7]) {
	      if(print > 1) cout<<"Didn't passed preliminary Chi2 cut (dCut[7])"<<endl;
	      goto next_pair;
	    }
	  }

	  if(print > 0){
	    cout<<"P = "<<1./tx.Hfirst(5)<<"\t  Err = "<<sqrt(tx.Hfirst(5,5))<<endl;
	  }

	  // Convergency test
	  int ichi2(0), ichi2_prev(0);
	  ichi2      = int(chi2     *10); // keep only 1 decimal digit
	  ichi2_prev = int(chi2_prev*10); // keep only 1 decimal digit

	  if(ichi2 == ichi2_prev) {
	    if(print > 1) cout<<"---> Converged !"<<endl;
	    break; // converged
	  }
	  
	  iter++;
	  if(TOpt::ReMode[6] > 0 ) break; // only one iteration for "ideal" bridging
	}; // <---------- end of iteration loop 

	if(iter == 15 && TOpt::Print[0] != 0) {
	  cout<<"TEv::BridgeSegments() ==> Maximum number of iterations ("<<iter<<") is reached. Mag = "<<imag;
	  if(good_pair) cout<<"  (good  pair) ";
	  if(bad_pair ) cout<<"  (wrong pair) ";
	  cout<<endl;
	  cout<<"1/P  prev = "<<pinv_prev<<"\t   1/P  last = "<<tx.Hfirst(5)<<" ("<<1./tx.Hfirst(5)<<" Gev)"<<endl;
	  cout<<"Chi2 prev = "<<chi2_prev<<"\t   Chi2 last = "<<chi2<<endl;
	  goto next_pair; // Next track pair.
	}

       	if(tx.Hfirst(5,5) < 0) { // abnormal cov. matrix element
	  cout<<"TEvBridgeSegments() ==> COV(1/P,1/P) < 0 !"<<endl;
	  goto next_pair; // Next track pair.
	}

	if(hist){
	  if(good_pair){
	    ///HF1(TOpt::Hoffs+100*imag+6,float(chi2),1.);
	    ///HF1(TOpt::Hoffs+imag*100+13,(1./fabs(pimc)),1.);
	  }
	}

	// Cut ( final Chi2)
	if(TOpt::ReMode[6] == 0 ) {
	  if(chi2 > TOpt::dCut[9]) {
	    if(print > 1) cout<<"Didn't passed final Chi2 cut (dCut[9])"<<endl;
	    goto next_pair;
	  }
	}
	// Histograms
	if(hist){
	  float pi, sig, dpinv;
	  if(good_pair){ // well reconstructed track
	    int ikin = tr.IKine;
	    pi = float(tx.Hfirst(5));
	    sig = float(sqrt(tx.Hfirst(5,5)));
	    dpinv=fabs(pi)-fabs(pimc);
	    ///HF1(TOpt::Hoffs+imag*100+1,dpinv,1.);
	    ///HF1(TOpt::Hoffs+imag*100+2,dpinv/sig,1.);
	    ///HF1(TOpt::Hoffs+imag*100+14,(1./fabs(pimc)),1.);
	    ///HF1(TOpt::Hoffs+imag*100+30,(iter+0.5),1.);
	  }
	  if(bad_pair){ // if wrong pair was bridged
	    ///HF1(TOpt::Hoffs+imag*100+16,(1./fabs(pimc)),1.);
	  }
	}
	
	// Store track pair candidates

	tp.Chi2 = chi2; 
	tp.iL   = il;
	tp.iR   = ir;
	tp.Pinv = tx.Hfirst(5);
	tp.ErrP = tx.Hfirst(5,5);
	if(imag == 1) lTrackPair1.push_back(tp);
	if(imag == 2) lTrackPair2.push_back(tp);
	

      next_pair:;
      }// end of loop over right track pieces
    }  // end of loop over left  track pieces

    tr.Id = ++TTrack::TrackCounter;// make IDs of temporary tracks unique
    tl.Id = ++TTrack::TrackCounter;// before destructors will be called
    tx.Id = ++TTrack::TrackCounter;// (to prevent removal of it's IDs from sTrackIDs of TKine)
    
  } // end of imag loop
  

  // Sort merged pairs by "Chi2 per hit" of the fit

  lTrackPair1.sort();
  lTrackPair2.sort();

  // Build global tracks

  for(itp1=lTrackPair1.begin(); itp1!= lTrackPair1.end(); itp1++){       // loop over 1-st magnet track pairs
    if((*(*itp1).iL).Type != 1 || (*(*itp1).iL).IMark == -1)  continue; // (IMark = -1 means "already appended")
    if((*(*itp1).iR).Type != 2 || (*(*itp1).iR).IMark == -1 ) continue;
    (*(*itp1).iL).Append((*(*itp1).iR));  // append track type 010 to track type 001
    (*(*itp1).iL).Hlast(5)  =(*itp1).Pinv;// store momentum in last point Helix
    (*(*itp1).iL).Hlast(5,5)=(*itp1).ErrP;
    // look if there is continuation after 2-d magnet
    for(itp2=lTrackPair2.begin(); itp2!= lTrackPair2.end(); itp2++){     // loop over 2-d magnet track pairs
      if((*itp1).iR != (*itp2).iL) continue; // track pairs sould have common track type 010
      (*(*itp1).iL).Append((*(*itp2).iR));   // append track type 100
      (*(*itp1).iL).Hlast(5)  =(*itp2).Pinv; // store momentum in last point Helix
      (*(*itp1).iL).Hlast(5,5)=(*itp2).ErrP;
      break;
    }
  }
  for(itp2=lTrackPair2.begin(); itp2!= lTrackPair2.end(); itp2++){     // loop over 2-d magnet track pairs
    if((*(*itp2).iL).Type != 2 || (*(*itp2).iL).IMark == -1) continue;
    if((*(*itp2).iR).Type != 4 || (*(*itp2).iR).IMark == -1) continue;
    (*(*itp2).iL).Append((*(*itp2).iR));  // append track type 100 to track type 010
    (*(*itp2).iL).Hlast(5)  =(*itp2).Pinv;// store momentum in last point Helix
    (*(*itp2).iL).Hlast(5,5)=(*itp2).ErrP;
  }
  

  // Erase tracks pieces and make backward/forward global fit 
  
  list<TTrack>::iterator it = listTrack.begin();
  while(it != listTrack.end()) { // track loop 
    if((*it).IMark == -1){
      listTrack.erase(it++); // erase appended track piece
    } else {
      if((*it).Type == 3 || (*it).Type == 6 || (*it).Type == 7){ // long track only
	if(
	   ! (*it).QuickKF(-1,1) ||  // global Kalman refit backward
	   ! (*it).QuickKF( 1,1)     // global Kalman refit  forward
	   )
	  {
	    listTrack.erase(it++); // erase track if the fit had been failed
	  }
      }
      it++;
    }
  }// end of track loop

  /*
  it = listTrack.begin();
  while(it != listTrack.end()) { 
      cout<<"after TTrackBridgeSegments  ID = "<<(*it).Id
	<<" Type = "<<(*it).Type<<" Nhit "<<(*it).NHits<<endl;
    it++;
  }
  */
}












