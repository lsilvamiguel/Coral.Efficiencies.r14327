/*!

  Dummy function for alternative track segments bridging
  Activated by setting option ReMode [12] to 1


  Currently, SiFi - MWPC-B bridging attempt.

*/

#include "TEv.h"
#include "TTrackPair.h"
#include "CsHistograms.h"

void TEv::BridgeSegments1()
{
  
  static bool first = true;
  static CsHist1D* h1[10];
  static CsHist2D* h2[10];
  if(first){
    first = false;
    CsHistograms::SetCurrentPath("/Traffic/BridgeSegments1");
    h1[0] = new CsHist1D( "SiFi-MWPC_00",  "DIP MWPC - DIP SiFi  (mrad)" ,  100,  -25.,  25.);
    h1[1] = new CsHist1D( "SiFi-MWPC_01",  "Y MWPC - Y SiFi  (cm)" ,        100,  -10.,  10.);
    h2[0] = new CsHist2D( "SiFi-MWPC_10",  "Y SiFi VS Y MWPC (cm)" ,        100, -25., 25., 100, -25., 25.);
    h2[1] = new CsHist2D( "SiFi-MWPC_11",  "DIP SiFi VS DIP MWPC (mrad)" ,  100, -20., 20., 100, -10., 10.);
    CsHistograms::SetCurrentPath("/");
  }

  TTrackPair tp;
  std::list<TTrackPair> lTP;
  std::list<TTrackPair>::iterator itp;

  THlx Ht;
  std::list<TTrack>::iterator it1,it2;
  for(it1 = listTrack.begin(); it1 != listTrack.end(); it1++){ // loop over MWPC TTracks

    THlx& H1 = (*it1).Hfirst;
    if(H1(0) < 3300.) continue;

    //cuts (see Z VS Y MWPC-B track distrib.)
    //if(H1(1) <  36. || H1(1) > 50.) continue;
    
    for(it2 = listTrack.begin(); it2 != listTrack.end(); it2++){ // loop over SiFi TTracks

      THlx& H2 = (*it2).Hlast;
      if(H2(0) > 350.) continue; 

      //cut
      //(fabs((*it2).MeanTime) > 3.0) continue; // off-time track

      double ddip = (H1.dip()-H2.dip())*1000.;
      h1[0]->Fill(ddip);
      double z1,z2,x0;
      x0 =  H1(0);
      z1 =  H1(2);
      z2 =  H2(2)+(x0-H2(0))*tan(H2.dip());
      double dz = z1-z2;
      h1[1]->Fill(dz);
      h2[0]->Fill(z1, z2);
      h2[1]->Fill(H1.dip()*1000., H2.dip()*1000.);

      // cuts (taken from prev. plots)
      if( fabs(dz)   >  3. )    continue;
      if( fabs(ddip) >  6. )    continue;

      // store track pair
      tp.Chi2 = fabs(dz); // tmp
      tp.iL   = it2;
      tp.iR   = it1;
      lTP.push_back(tp);

    }
  }

  if(lTP.size() == 0) return;

  lTP.sort();
  
  // Append track (only first (best) pair is taken. tmp)
  TTrack& tr = *(lTP.front().iR);
  TTrack& tl = *(lTP.front().iL);
  tl.Append(tr);
  
  std::list<TTrack>::iterator it = listTrack.begin();
  while(it != listTrack.end()) { // track loop 
    if((*it).IMark == -1){  // alredy appended track piece
      listTrack.erase(it++);     // erase; next track
    } else {                // normal track
      if((*it).Type == 9){      // long track only
	(*it).Hfirst(5)=1./200.; // mom. estimation
	if(! (*it).QuickKF( 1,1) ||  // global Kalman refit backward
	   ! (*it).QuickKF(-1,1)) {  // global Kalman refit  forward
	  listTrack.erase(it++); // erase track if the fit had been failed
	} else { // Chi2 cut
	  //if((*it).Chi2tot/((*it).NHits-5.) > 100. ){
	  //listTrack.erase(it++); 
	  //}
	}
      }
      it++; // next track
    }
  }// end of track loop

}
