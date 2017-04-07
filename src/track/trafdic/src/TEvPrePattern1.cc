/*!

  Special time-based prepattern 
  for zones with 4 SciFi planes only


*/

#include "TEv.h"
#include "TAlgo.h"
#include "TOpt.h"
#include "CsZone.h"
#include "CsHistograms.h"

using namespace std;

void TEv::PrePattern1(CsZone* zone)
  
{

  const TSetup& setup = TSetup::Ref();

  int igroup = -1; 
  for(igroup = 0; igroup < int(setup.vIplFirst().size()); igroup++){ //loop over all det. groups
    if(setup.Group2Zone(igroup) == zone) break; // reconstruction in this group was requested
  }

  int ipl1 = setup.vIplFirst()[igroup];
  int ipl2 = setup.vIplLast ()[igroup];
  int npl_all = ipl2-ipl1+1;

  int ipl_on[npl_all]; // map(switched on plane # --> vecPlane index) 
  int Nplanes=0;       // N switched on planes
  for(int i=ipl1; i<=ipl2; i++){
    if(setup.vPlane(i).IFlag != 0) ipl_on[Nplanes++]=i;
  }

  static CsHist1D*  SFTrk   [100];
  static CsHist2D*  SFTrk2D [100];

  static bool first=true;
  if(first){
    first=false;
    if(TOpt::Print[0] != 0)
       cout<<"TEvPrePattern1() ==> Time-based pre-pattern is used for zone # "<<igroup<<endl;

    CsHistograms::SetCurrentPath("/Traffic/PrePattern1");
    SFTrk[10]  = new CsHist1D( "sftrk10","delta T of all hit pairs of track candidates",200, -10, 10);
    SFTrk[11]  = new CsHist1D( "sftrk11","mean T of all track candidates",200, -10, 10);
    SFTrk[12]  = new CsHist1D( "sftrk12","mean T of all track candidates",200, -50, 50);

    // Book time correlation/autocorrelation plots
    char id[20];
    for(int ip1 = 0; ip1 < Nplanes; ip1 ++){
      sprintf(id,"sftrk%02u",   20+ip1);
      SFTrk  [20+ip1] = new CsHist1D( id, "delta T of all hit pairs on _the_same_ plane", 200, -20, 20);
      sprintf(id,"sftrk2d%02u", 20+ip1);
      SFTrk2D[20+ip1] = new CsHist2D( id, "delta(T) VS delta(Nfiber)", 96, 0, 96, 100, -20, 20);
    }

    // track params plots
    SFTrk2D[0]  = new CsHist2D( "sftrk2d00","Off-time: Y  VS X  (cm)",    100,  -5,  5, 100,  -5,  5);
    SFTrk2D[1]  = new CsHist2D( "sftrk2d01","Off-time: Y' VS X' (mrad)",  100,  -8,  8, 100,  -8,  8);
    SFTrk2D[2]  = new CsHist2D( "sftrk2d02","Off-time: X' VS X",          100,  -5,  5, 100,  -8,  8);
    SFTrk2D[3]  = new CsHist2D( "sftrk2d03","Off-time: Y' VS Y",          100,  -5,  5, 100,  -8,  8);

    SFTrk2D[10] = new CsHist2D( "sftrk2d10","In-time : Y  VS X  (cm)",    100,  -5,  5, 100,  -5,  5);
    SFTrk2D[11] = new CsHist2D( "sftrk2d11","In-time : Y' VS X' (mrad)",  100,  -8,  8, 100,  -8,  8);
    SFTrk2D[12] = new CsHist2D( "sftrk2d12","In-time : X' VS X",          100,  -5,  5, 100,  -8,  8);
    SFTrk2D[13] = new CsHist2D( "sftrk2d13","In-time : Y' VS Y",          100,  -5,  5, 100,  -8,  8);

    CsHistograms::SetCurrentPath("/");
  }


  // Fill time autocorrelation plots
  for(int ip1 = 0; ip1 < Nplanes; ip1 ++){
    const TPlane&  p1 = setup.vPlane(ipl_on[ip1]);
    vector<int>::const_iterator i1,i2;
    for(i1=p1.vHitRef().begin(); i1 != p1.vHitRef().end(); i1++){ // loop over hit references on the plane p1
      THit& h1 = vecHit[*i1];   // ref to THit vector element
      if(h1.sDigits().empty()) continue;
      const TDigit& d = *(h1.sDigits().begin()); 
      if(d.vDigInfo.empty()) continue;
      float t1 = d.vDigInfo[0]; float w1=d.IWire;
      for(i2 = i1; i2 != p1.vHitRef().end(); i2++){ // loop over hit references on the plane p1
	if( (*i1) == (*i2) ) continue; // same hit
	THit& h2 = vecHit[*i2]; // ref to THit vector element
	if(h2.sDigits().empty()) continue;
	const TDigit& d = *(h2.sDigits().begin()); 
	if(d.vDigInfo.empty()) continue;
	float t2 = d.vDigInfo[0]; float w2=d.IWire;
	SFTrk  [20+ip1]->Fill((t2-t1));
	SFTrk2D[20+ip1]->Fill(fabs(w2-w1), (t2-t1));
      }
    }
  }

  //
  // Buld track candidates like all possible combinations of hits on _all_ planes of the group
  //

  // pack hits and hit times to local arrays
  int Nhits[Nplanes];
  THit** hp[Nplanes];
  float* t [Nplanes];
  for(int ip = 0; ip < Nplanes; ip++){ // loop over plane in group
    int ipl = ipl1+ip; // plane #
    const TPlane&  p = setup.vPlane(ipl_on[ip]);
    Nhits[ip] = p.vHitRef().size(); // store N hits per plane
    hp[ip] = new THit* [Nhits[ip]]; // working array with THit*'s
    t [ip] = new float [Nhits[ip]]; // working array with hit times
    int nhit=0;
    for(int ih = 0; ih < Nhits[ip]; ih++){ // loop over hits on plane
      THit* h = &vecHit[p.vHitRef()[ih]];
      hp[ip][nhit] = h;
      THit& h0 = vecHit[p.vHitRef()[ih]];
      if(h0.SigT < 0) continue; // time not measured
      t [ip][nhit] = h0.Time; 
      // here one can put some cuts
      //if(...) continue;
      nhit++;
    }
    Nhits[ip]=nhit; // redefine as some hits could be skiped
  }

  // Do combinatorics

  list<TTrack> tt;
  int icomb[Nplanes];
  while(1){ // loop over hit combinations (track candidates with Nhits = Nplanes)

    bool ok=true;

    if(! TAlgo::NLoopComb(Nplanes, Nhits, icomb)) break;

    // Take every pair if hits in combination to look on time difference
    float t1,t2,dt;
    for(int ip1 = 0; ip1 < Nplanes-1; ip1 ++){
      t1 =  t[ip1][icomb[ip1]];
      for(int ip2 = ip1+1; ip2 < Nplanes; ip2 ++){
	t2 =  t[ip2][icomb[ip2]];
	dt = (t2-t1);
	SFTrk[10]->Fill(dt);
	// if time difference for _any_ pair of hits of the track candidate
	// is out of cut, goto next track candidate (hit combination)
	if(fabs(dt) > TOpt::dCut[21]) { ok=false; goto exit;}; // Cut

      }
    }

  exit:
    if(ok){
      // create track candidate
      TTrack tr;
      tr.Type=1<<igroup;
      // add hits to track
      for(int ip = 0; ip < Nplanes; ip ++){
	tr.AddHit(*hp[ip][icomb[ip]]);
      }
      // preliminary fit by straight line
      tr.QuickKF(1,0);
      if( !tr.UseHitTime() ) continue; // times are incompatible. Skip the track.
      SFTrk[11]->Fill(tr.MeanTime);
      SFTrk[12]->Fill(tr.MeanTime);
      tt.push_back(tr); // store found track candidate 
    }

  }// end of loop over hit combinations (track candidates)

  // tracks filtering out
  list<TTrack>::iterator it,it1,it2;
  for(it1 = tt.begin(); it1 != tt.end(); it1++){
    for(it2 = it1; it2 != tt.end(); it2++){
    }
  }

  for(it = tt.begin(); it != tt.end(); it++){
    if((*it).IMark != -1) this->listTrack.push_back(*it); // store found track into event
  }



  // study found track properties
  THlx H0;
  for(it = listTrack.begin(); it != listTrack.end(); it++){
    TTrack& tr = (*it);
    if(tr.Type != 1<<igroup) continue; // not my track
    //H0(0)=0; tr.Hlast.Extrapolate(H0);
    H0=tr.Hlast;
    if(fabs(tr.MeanTime) > 3.0) { // off-time tracks (out of +- 3 ns to trigger time)
      SFTrk2D[0] ->Fill(float(        H0(1)), float(        H0(2)));
      SFTrk2D[1] ->Fill(float(1000. * H0(3)), float(1000. * H0(4)));
      SFTrk2D[2] ->Fill(float(        H0(1)), float(1000. * H0(3)));
      SFTrk2D[3] ->Fill(float(        H0(2)), float(1000. * H0(4)));
    } else { // in-time tracks
      SFTrk2D[10]->Fill(float(        H0(1)), float(        H0(2)));
      SFTrk2D[11]->Fill(float(1000. * H0(3)), float(1000. * H0(4)));
      SFTrk2D[12]->Fill(float(        H0(1)), float(1000. * H0(3)));
      SFTrk2D[13]->Fill(float(        H0(2)), float(1000. * H0(4)));
    }
  } 


  // free space of localy reserved arrays
  for(int ip = 0; ip < Nplanes; ip++){ // loop over planes in group
    delete[] hp[ip];
    delete[] t [ip];
  }


};


















