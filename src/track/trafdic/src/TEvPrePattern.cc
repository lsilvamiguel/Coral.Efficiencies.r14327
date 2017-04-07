/*!
  \brief Preliminary pattern recognition
  Track segments finding in different detector groups

*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include "CsZone.h"
#include "TEv.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TDisplay.h"
#include "TAlgo.h"
#include "TConstants.h"

using namespace std;

// prototypes

#ifdef TIME_CHECK
#include "CsStopwatch.h"

extern bool TiTrack;
extern double find_proj[5];
extern int n_hit[5]; 
#endif


void TEv::PrePattern(CsZone* zone)
  
{
  const TSetup& setup = TSetup::Ref();

  int igroup, ipass; 
  for(igroup = 0; igroup < int(setup.vIplFirst().size()); igroup++){ //loop over all det. groups
    if(setup.Group2Zone(igroup) == zone) break; // it's requested zone
  }
  //
  // "ideal" pattern recognition of the track segments (for MC only)
  //
  if(TOpt::ReMode[4] != 0){
    for(int itk=0; itk < int(vecKine.size()); itk++){ // loop over KINE tracks
      TTrack t; // create new track
      t.Type=1<<igroup;
      vector<int>::const_iterator ivh;
      for(ivh=vecKine[itk].vHitRef().begin(); ivh!=vecKine[itk].vHitRef().end(); ivh++){ // loop over MC track hits
	if(vecHit[(*ivh)].IOrig !=0) continue;  // skip not original hits
	int ipl=vecHit[(*ivh)].IPlane;
	if(setup.vPlane(ipl).IFlag == 0) continue; // this plane had been switched off
	if(ipl > setup.vIplLast()[igroup] || ipl < setup.vIplFirst()[igroup] ) continue; // only current group hits
	t.AddHit(vecHit[(*ivh)]);               // add hit to the track
      } //end of look over MC track hits
      // cout<< "MC track # "<<itk<< " group = "<<igroup<<"  Nhits = "<<t.NHits<<endl;

      if(int(t.NHits) >= TOpt::iPRpar[igroup*10 + 2]){ // ususl cut on Nhits
	t.FindKine(); // find corresponding Kine track (just a test of procedure)
	listTrack.push_back(t); // store the track
      }
    } // end of loop over KINE tracks

    return;
  }

  //--------------------------------------------------------
    
    // special case of only 4 SciFi planes
  int nFI(0), nOther(0);
  for(int ipl=setup.vIplFirst()[igroup]; ipl <= setup.vIplLast()[igroup]; ipl++){ // loop over planes in group
    if(setup.vPlane(ipl).IFlag==0) continue; // skip switched off detector
    string name = setup.iPlane2Detect(ipl).Name;
    if(name.find("FI") == 0) nFI++;
    else                     nOther++;
  }
  if(nFI == 4 && nOther == 0){
    PrePattern1(zone);
    return;
  }

  //--------------------------------------------------------


  // params

  int maxhit  = 3000; 
  int maxpl   = 300;
  int maxntk  = TConstants_NTtrack_max;
  
  int maxproj; // projections to use for pat. rec.
  maxproj = setup.vProj().size(); // all proj.

  // working arrays
  float xpl [maxpl];               // X of plane
  int   idpl [maxpl];              // index of plane in vDetect
  float res [maxpl];               // resolution of plane
  float tol [maxpl];               // tolerances for track finding of plane
  float sigt [maxpl];              // time resolution of plane squared
  float ymax[maxpl];               // max coord on plane
  float ymin[maxpl];               // min coord on plane
  int    fh [maxpl], lh[maxpl];    // first and after-the-last hit number on plane
  int   iprj[maxpl];               // prijection index on the plane
  float cosa[maxpl];               // cos(a)
  float sina[maxpl];               // sin(a) of planes
  float Uhit[maxhit];              // Y of hit
  float Thit[maxhit];              // time of hit
  int   href[maxhit];              // ref to vHit
  int   hits[(maxntk+1)*maxpl];    // array of found  hits on proj. tr.
  float y0[maxntk], yp[maxntk];    // parameters of found proj. tr.
  float tT[maxntk], sT [maxntk];   // time and sigma  of found proj. tr.

  float X0 = float(setup.iPlane2Detect(setup.vIplFirst()[igroup]).X(0)); // reference plane

  for(int ipass = 1; ipass <= TOpt::iPRpar[50]; ipass++) {// loop over pattern recognition passes

    int nproj = 0;
    int projind[maxproj]; // array of indecies of projections in this detector group

    int Ntk,nhit,npl,ier;
  
    // just counting planes in proj
    for(int ipr=0; ipr < maxproj; ipr++){ // loop over all pat. rec. projectons
      int npl=0; 
      for(int ipl=setup.vIplFirst()[igroup]; ipl <= setup.vIplLast()[igroup]; ipl++){ // loop over planes in group
	const TPlane&  p = setup.vPlane(ipl);
	if(p.IFlag == 0) continue;        // plane is off
	if(ipr != p.IProj) continue;      // not current projection
	npl++;
      }
      if(npl == 0) continue; // no planes of projection "ipr" in this group
      if(npl < TOpt::iPRpar[igroup*10+0]) continue; // not enough planes in this projection
      projind[nproj++]=ipr;
    } // end of loop over projections    

    int Ntk_prj [nproj];
    float cos_prj[nproj], sin_prj[nproj]; // to store Cos(a) and Sin(a) of projections
    if(nproj > 10){ 
      cout<<"TEv::PrePattern ==> Number of projection in detector group "<<igroup <<" is "<<nproj<<endl;
      assert(false);
    }
    // arrays to store results of track finding in projections
    float Y0_prj[nproj][maxntk];
    float Yp_prj[nproj][maxntk];
    float tT_prj[nproj][maxntk];
    float sT_prj[nproj][maxntk];
    
    //
    // Fill working arrays for track finding in projections
    //
    
    for(int ipr=0; ipr < nproj; ipr++){  // loop over pat. rec. pojectons 
      int iproj = projind[ipr];
      int nhit=0, npl=0; //reset counters
      vector<int>::const_iterator i;
      float cos_mean=0, sin_mean=0; // to calculate mean Cos(a), Sin(a) over planes in projection
      double worse_res = 0;
      for(int ipl=setup.vIplFirst()[igroup]; ipl <= setup.vIplLast()[igroup]; ipl++){ // loop over planes
	const TPlane&  p = setup.vPlane(ipl);
	if(p.IFlag == 0) continue; // plane is off
	if(iproj != p.IProj) continue;  // not required projection
	idpl[npl] = p.IDetRef;
	const TDetect& d = setup.vDetect(p.IDetRef);
	cos_mean+=d.Ca; sin_mean+=d.Sa;
      
	xpl[npl] =d.X(0);
	ymin[npl]=d.XR(1)-d.Range/2;
	ymax[npl]=d.XR(1)+d.Range/2;
      
	res[npl] = d.Resol;

	if(d.TResol < 0) sigt[npl] = -1; // means "detector do not measure time"
	else             sigt[npl] = d.TResol*d.TResol;

	if(res[npl] > worse_res) worse_res = res[npl];
	fh [npl]=nhit;
	for(i=p.vHitRef().begin(); i != p.vHitRef().end(); i++){      // loop over hit references on the plane
	  THit& h = vecHit[*i]; // ref to Thit vector element
	  if( !h.sTrackID().empty() ) continue; // hit had been already used
	  Uhit[nhit]=h.U; href[nhit]=(*i);
	  Thit[nhit]=h.Time;
	  if(++nhit == maxhit){ 
	    if(TOpt::Print[0] != 0) 
	      cout<<"TEv::PrePattern ==> (Proj).  More then "<<maxhit<<" hits in group "<<igroup<<endl;
	    return;
	  }
	}  // end of loop over hits on the plane
	lh[npl]=nhit;
	if(++npl == maxpl){ 
	  cout<<"TEv::PrePattern ==> (Proj).  More then "<<maxpl<<" planes in group "<<igroup<<endl;
	  assert(false);
	}
      }   // end of loop over planes

      cos_prj[ipr]=cos_mean/npl; sin_prj[ipr]=sin_mean/npl;
      if(TOpt::Print[5]){
	printf("Det. group %1u    Proj %2u   Nplanes %3u   NHits %5u  ", igroup,iproj,npl,nhit); 
	cout<<flush;
      }
    
      //
      // reconstruction in one projection
      //
    
      int ier(0); Ntk = 0;
#ifdef TIME_CHECK
  CsStopwatch stopwatch;	  
  int star=stopwatch.start();
#endif      
      if(false){
	///// some other algorithm ...
      } else {
	for(int i=0; i < npl;  i++) { // set tolerances 
	  //tol[i]=worse_res*TOpt::dPRpar[igroup*10+0]; // worse resolution * coeff.
	  tol[i]=res[i]*TOpt::dPRpar[igroup*10+0];    // plane resolution * coeff.
	} 
	Ntk = maxntk; // Ntk is "bidirectional" parameter		
	//cout<<" recostruction in proj "<<ipr<<endl;
	ier = TAlgo::FindProj(igroup, ipass, X0, npl, xpl, tol, sigt, fh, lh, nhit, Uhit, Thit, Ntk, hits, y0, yp, tT, sT);
      }
#ifdef TIME_CHECK
if(TiTack){
     find_proj[igroup]+=stopwatch.stop(star);
     n_hit[igroup]+=nhit;
}
#endif      
      if(ier > 1) { 
	cout<<"TEv::PrePattern ==> (Proj) error in group "<<igroup<<endl;
	return;
      }

      // store projection track parameters
      Ntk_prj[ipr]=Ntk;
      for(int jtr =0; jtr < Ntk; jtr++){
	Y0_prj[ipr][jtr]=y0[jtr];
	Yp_prj[ipr][jtr]=yp[jtr];
	tT_prj[ipr][jtr]=tT[jtr];
	sT_prj[ipr][jtr]=sT[jtr];
      }
      
      if(TOpt::ReMode[5] != 0){ // save projection tracks
	// Store tracks
	for(int it=0, nh=0; it < Ntk; it++){ // loop over found tracks
	  TTrack t;
	  t.Type=1<<igroup;
	  // store hits
	  while(1){
	    if(hits[nh] == -1) {nh++; break;}
	    if(href[hits[nh]] < 0 || href[hits[nh]] >= int(vecHit.size()) ){
	      cout<<"TEv::PrePattern() proj. ==> something is wrong with hit referencies"<<endl;
	      assert(false);
	    }
	    t.AddHit(vecHit[href[hits[nh]]]); nh++;
	  }
	  listTrack.push_back(t); // store the track
	} // end of loop over found tracks
      } //end of "if ReMode[5] != 0" 

    } // end of loop over pat. rec.  projections in this group


    //
    // Fill arrays for track finding in space
    //
    if(TOpt::ReMode[5] != 0) return;

    vector<int>::const_iterator i;
    double worse_res=0; // will be worse det. res in group
    nhit=0; npl=0; // reset counters
    for(int ipl=setup.vIplFirst()[igroup]; ipl <= setup.vIplLast()[igroup]; ipl++){ //loop over all planes in the group

      const TPlane&  p = setup.vPlane(ipl);
      if(p.IFlag == 0) continue; // plane is off
      const TDetect& d = setup.vDetect(p.IDetRef);

      idpl[npl] = p.IDetRef;
      iprj[npl] = p.IProj;
      cosa[npl]=d.Ca; sina[npl]=d.Sa;

      res[npl]=d.Resol;

      if(d.TResol < 0) sigt[npl] = -1; // means "detector do not measure time"
      else             sigt[npl] = d.TResol*d.TResol;

      if(res[npl] > worse_res) worse_res = res[npl];

      fh [npl]=nhit;
      for(i=p.vHitRef().begin(); i != p.vHitRef().end(); i++){      // loop over hit references on the plane
	THit& h = vecHit[*i]; // ref to Thit vector element
	if( !h.sTrackID().empty() ) continue; // hit had been already used
	Uhit[nhit]=h.U; href[nhit]=(*i);
	Thit[nhit]=h.Time;
	if(++nhit == maxhit){
	  if(TOpt::Print[0] != 0) 
	    cout<<"TEv::PrePattern ==> (Space).  More then "<<maxhit<<" hits in group "<<igroup<<endl;
	  return;
	}
      }  // end of loop over hits on the plane

      lh[npl]=nhit;
      if(++npl == maxpl){
	cout<<"TEv::PrePattern ==> (Space).  More then "<<maxpl<<" planes in group "<<igroup<<endl;
	assert(false);
      }
      //cout<<" igr = "<<igroup<<" npl = "<<npl<<" id pl. "<<d.IDet
      //<<" Nhits = "<<p.vHitRef().size()<<" lh-fh = "<<lh[npl-1]-fh[npl-1]<<endl;

    } // end of loop over planes


    // Fill alignment histograms (if requested)
    TAlgo::Alignment(igroup, npl, idpl, iprj, cosa, sina, tol, fh, lh, nhit, Uhit);


    for(int i=0; i < npl;  i++) { // set tolerances 
      //tol[i]=worse_res*TOpt::dPRpar[igroup*10+0]; // worse resolution * coeff.
      tol[i]=res[i]*TOpt::dPRpar[igroup*10+0];    // plane resolution * coeff.
    } 
    Ntk = maxntk; // Ntk is "bidirectional" parameter
    ier = TAlgo::FindSpace(igroup, ipass, X0, npl, idpl, iprj, cosa, sina, tol, sigt, fh, lh, nhit, Uhit, Thit, 
			 nproj, Ntk_prj, cos_prj, sin_prj, (float**)Y0_prj, (float**)Yp_prj,
			 (float**)tT_prj, (float**)sT_prj, Ntk, hits);
    

    if(ier != 0) {
      cout<<"TEv::PrePattern ==> (Space) error in group "<<igroup<<endl;
      return;
    }


    // Store tracks
  
    for(int it=0, nh=0; it < Ntk; it++){ // loop over found tracks
      TTrack t;
      t.Type  = 1<<igroup;
      t.IMark = ipass; // just for information
      // store hits

      while(1){
	if(hits[nh] == -1) {nh++; break;} // next track's hits
	if(href[hits[nh]] < 0 || href[hits[nh]] >= int(vecHit.size()) ){
	  cout<<"TEv::PrePattern() space ==> something is wrong with hit referencies"<<endl;
	  assert(false);
	}
	t.AddHit(vecHit[href[hits[nh]]]); nh++;
      }
      // cut
      if(t.NHits < 4 ) continue; // it's hard to call it "Track" :-)
      if(t.sProj.size() < 3) {
	cout<<"TEv::PrePattern ==> "<<t.sProj.size()<<"-projections space track !"<<endl;
	t.Print(2);
      }
      if(TOpt::dCut[2] != 0.0) { // Chi2 cut is ON
	// fit by straight line just for Chi2 selection
	t.QuickKF( 1,0); 
	t.QuickKF(-1,0);

	// Chi2 cut
	if(t.Chi2tot/t.NHits > TOpt::dCut[2]) continue; // skip
      }
      // Time cuts inside function
      if( !t.UseHitTime() ) continue; // times are incompatible. Skip the track.

      // store the track
      listTrack.push_back(t); 

    } // end of loop over found tracks
    

  } // end of loop over passes  
  
  // some final operations with found tracks
  list<TTrack>::iterator it;
  for(it =  listTrack.begin(); it != listTrack.end(); it++){
    (*it).FindKine(); // in case of MC data, find corresponding Kine track
  }

};
















