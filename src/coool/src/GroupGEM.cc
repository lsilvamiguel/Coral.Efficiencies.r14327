#include <cmath>
#include "GroupGEM.h"
#include "TThread.h"
#include "GroupPanel.h"
#include "PlaneGEM.h"
#include "TriggerTime.h"
#include "DaqEvent.h"
#include "GeomPlaneGEM.h"

ClassImp(GroupGEM);
ClassImp(GroupGEMTrack);

void GroupGEM::Init() {
  std::string name0 = fName + "_U_vs_V";
  fHistList.push_back(new TH2F(name0.c_str(), name0.c_str(),
			       320, -16, 16, 320, -16, 16));
  fHistList[0]->SetOption("COLZ");
  fHistList[0]->SetStats(false);

  std::string name1 = fName + "_Y_vs_X";
  fHistList.push_back(new TH2F(name1.c_str(), name1.c_str(),
			       320, -16, 16, 320, -16, 16));
  fHistList[1]->SetOption("COLZ");
  fHistList[1]->SetStats(false);

  std::string name2 = fName + "_U_V_AmplCorr";
  fHistList.push_back(new TH2F(name2.c_str(), name2.c_str(),
			       50, 0, 1000, 50, 0, 1000));
  fHistList[2]->SetOption("COLZ");

  std::string name3 = fName + "_Y_X_AmplCorr";
  fHistList.push_back(new TH2F(name3.c_str(), name3.c_str(),
			       50, 0, 1000, 50, 0, 1000));
  fHistList[3]->SetOption("COLZ");

  std::string name4 = fName + "_A2U/A2V";
  fHistList.push_back(new TH1F(name4.c_str(), name4.c_str(),
			       50, -2.45, 2.55));

  std::string name5 = fName + "_A2Y/A2X";
  fHistList.push_back(new TH1F(name5.c_str(), name5.c_str(),
			       50, -2.45, 2.55));
}

void GroupGEM::EndEvent(const CS::DaqEvent &event) {

  if(fPlanes.size()<2) return;
  if (thr_flag) TThread::Lock();

  for (int k=0; k<2; k++) {

    // k: 0-VU 1-XY

    const GeomPlane *xgeom = NULL;
    const GeomPlane *ygeom = NULL;

    for (register unsigned int pl=0; pl<fPlanes.size(); pl++) {
      char planeLetter = fPlanes[pl]->GetName()[4];
      if (planeLetter =='V' && k==0) xgeom = fPlanes[pl]->GetGeometry();
      if (planeLetter =='U' && k==0) ygeom = fPlanes[pl]->GetGeometry();
      if (planeLetter =='X' && k==1) xgeom = fPlanes[pl]->GetGeometry();
      if (planeLetter =='Y' && k==1) ygeom = fPlanes[pl]->GetGeometry();
    }


    if (xgeom && ygeom) {
      const std::set<CCluster1>& xclu = xgeom->GetClusters();
      const std::set<CCluster1>& yclu = ygeom->GetClusters();

      int xsize = xclu.size();
      int ysize = yclu.size();
      float xamp2[xsize];
      float yamp2[ysize];
      float xpos[xsize];
      float ypos[ysize];
      bool xused[xsize], yused[ysize];

      int nx=0;
      for (std::set<CCluster1>::iterator xcl=xclu.begin();
	   xcl!=xclu.end(); xcl++) if (xcl->dt.size()>=3) {
	     xamp2[nx] = xcl->dt[2];
	     xpos[nx] = xcl->pos;
	     xused[nx] = false;
	     nx++;
	   }

      int ny=0;
      for (std::set<CCluster1>::iterator ycl=yclu.begin();
	   ycl!=yclu.end(); ycl++) if (ycl->dt.size()>=3) {
	     yamp2[ny] = ycl->dt[2];
	     ypos[ny] = ycl->pos;
	     yused[ny] = false;
	     ny++;
	   }
      
      int msize = nx*ny;
      int xindex[msize], yindex[msize];
      float distance[msize];
      for (int i=0; i<msize; i++) distance[i]=10000;

      for (int ix=0; ix<nx; ix++)
	for (int iy=0; iy<ny; iy++) {

	  // Distance to yamp2 = xamp2
	  float this_dist = fabs(yamp2[iy]-xamp2[ix])/sqrt(2.);
	  
	  // Maximum deviation (this_dist < slope*average+offset)
	  this_dist = this_dist - 0.05*(xamp2[ix]+yamp2[iy])/2.;
	    
	  // Sort according to this_dist
	  for (int i=0; i<msize; i++){
	    if (this_dist<distance[i]) {
	      for (int ii=msize-2; ii>=i; ii--){
		distance[ii+1] = distance[ii];
		xindex[ii+1] = xindex[ii];
		yindex[ii+1] = yindex[ii];
	      }
	      distance[i]=this_dist;
	      xindex[i]=ix;
	      yindex[i]=iy;
	      break;
	    }
	  }
	}
      
      // Associate clusters accordingly
      for (int i=0; i<msize; i++) {
	
	// Check if within road
	if (distance[i]<80.) {
	  
	  // Clusters already associated
	  if (!xused[xindex[i]] && !yused[yindex[i]]) {
	    fHistList[k  ]->Fill(xpos[xindex[i]], ypos[yindex[i]]);
	    fHistList[k+2]->Fill(xamp2[xindex[i]], yamp2[yindex[i]]);
	    fHistList[k+4]->Fill(yamp2[yindex[i]]/xamp2[xindex[i]]-1.);
	    xused[xindex[i]]=true;
	    yused[yindex[i]]=true;
	  }
	}
	else break;
      }
    }
  }
  
  if (thr_flag) TThread::UnLock();
}

void GroupGEM::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}


void GroupGEMTrack::Init() {

#if USE_TRACK == 1
  // Tracking for the GEMs
  //fTrackMaker = new TrackMaker(GetName(),2,1e-15);
  fTrackMaker = new TrackMaker(GetName(), 2, 0);
  // 2: not space points (chi2)
  // 0.05 cut on chi2 probability

  //sprintf(sciName[0], "FI06X1__");
  //sprintf(sciName[1], "FI07X1__");
  sprintf(sciName[0], "FI05X1__");
  sprintf(sciName[1], "FI05Y1__");

  for(unsigned i=0; i<fPlanes.size(); i++) {
    if ( strcmp(fPlanes[i]->GetName(), "HMSC1") == 0 ) {
      jitterClock = (PlaneTrigHodo*) fPlanes[i];
      std::cout << "Found jitter clock information in "
	   <<jitterClock->GetName()<<std::endl;
    }
    for (int scpl=0; scpl<2; scpl++) {
      if ( strcmp(fPlanes[i]->GetName(), sciName[scpl]) == 0 ) {
	scifiTimePlane[scpl] = (PlaneScifiG*) fPlanes[i];
	std::cout << "Get SciFi1 timing information from "
	     <<scifiTimePlane[scpl]->GetName()<<std::endl;
      }
    }

    const GeomPlane *geom = fPlanes[i]->GetGeometry();

    if (geom) {
      int id = geom->GetID();
      const char *name = geom->GetName();
      int nwires = geom->GetNWires();
      float x = geom->GetX(); 		
      float y = geom->GetY();
      float z = geom->GetZ();
      float dx = geom->GetDX();
      float dy = geom->GetDY();
      float dz = geom->GetDZ();
      float angle = geom->GetAngle();
      float inpitch = geom->GetPitch();
      float ctmin = -500;
      float ctmax =  500;
      fTrackMaker->AddTracker(id,name,nwires,x,y,z,dx,dy,dz,angle,
			      inpitch,inpitch,
			      ctmin,ctmax,0);
    }
  }

  // Clusters per planes allowed
  maxNumberofClusters=3;
  ifstream mNumbCl("maxClusterGEM.num");
  if (mNumbCl) mNumbCl >> maxNumberofClusters;
  mNumbCl.close();
  std::cout << "taking " << maxNumberofClusters << " as maximum number of clusters"
       << std::endl;

  fTrackMaker->Init();
  fTrackMaker->SetNTrackers(6);
  //  fTrackMaker->Dump();

  // Adds tracking histogram to this's histogram list
  std::vector<TH1*>& trackhists=fTrackMaker->GetHistograms();
  for(unsigned i=0; i<trackhists.size();i++)
    fHistList.push_back(trackhists[i]);

  std::string name1 = fName + "_tcsJitter";
  fHistList.push_back(new TH1F(name1.c_str(), name1.c_str(), 160, -40, 40));
  hTCSjitter = fHistList.back();

  std::string name2 = fName + "_Sci1Time";
  fHistList.push_back(new TH1F(name2.c_str(), name2.c_str(),
			       1000, -200, 200));
  hSciTime[0] = fHistList.back();

  std::string name6 = fName + "_Sci2Time";
  fHistList.push_back(new TH1F(name6.c_str(), name6.c_str(),
			       200, -200, 200));
  hSciTime[1] = fHistList.back();

  std::string name3 = fName + "_Sci1_TvcCh";
  fHistList.push_back(new TH2F(name3.c_str(), name3.c_str(),
			       400, -3.25, 196.75, 200, -17000, -14000));
  hSciTimeVsChan[0] = (TH2*) fHistList.back();

  std::string name7 = fName + "_Sci2_TvcCh";
  fHistList.push_back(new TH2F(name7.c_str(), name7.c_str(),
			       400, -3.25, 196.75, 200, -17000, -14000));
  hSciTimeVsChan[1] = (TH2*) fHistList.back();

  std::string name4 = fName + "_GEMA1dA2VsT";
  fHistList.push_back(new TH2F(name4.c_str(), name4.c_str(), 140, -60, 80,
			       84, -0.1, 2.0));
  hGEMRatio12VsSci = (TH2*) fHistList.back();

  std::string name43 = fName + "_GEMTi12VsT";
  fHistList.push_back(new TH2F(name43.c_str(), name43.c_str(),
			       120, -120, 120, 84, -0.1, 2.0));
  hGEMTime12VsSci = (TH2*) fHistList.back();

  std::string name5 = fName + "_GEMTiming12";
  fHistList.push_back(new TH1F(name5.c_str(), name5.c_str(), 500, -250, 250));
  hGEMTime12 = fHistList.back();

  std::string name42 = fName + "_GEMA0dA2VsT";
  fHistList.push_back(new TH2F(name42.c_str(), name42.c_str(), 140, -60, 80,
			       84, -0.1, 2.0));
  hGEMRatio02VsSci = (TH2*) fHistList.back();

  std::string name44 = fName + "_GEMTi02VsT";
  fHistList.push_back(new TH2F(name44.c_str(), name44.c_str(),
			       100, -100, 100, 84, -0.1, 2.0));
  hGEMTime02VsSci = (TH2*) fHistList.back();

  std::string name52 = fName + "_GEMTiming02";
  fHistList.push_back(new TH1F(name52.c_str(),name52.c_str(), 500, -250, 250));
  hGEMTime02 = fHistList.back();

  std::string name53 = fName + "_GEMTimingCorr";
  fHistList.push_back(new TH2F(name53.c_str(),name53.c_str(),
			       100, -250, 250, 100, -250, 250));
  hGEMTimeCorr = (TH2*) fHistList.back();

  std::string name8 = fName + "_Chi22";
  fHistList.push_back(new TH1F(name8.c_str(), name8.c_str(), 1000, 0, 15));
  hChi22 = fHistList.back();

  std::string name9 = fName + "_dT_scifi";
  fHistList.push_back(new TH1F(name9.c_str(), name9.c_str(),
			       500, -10000, 10000));
  hdTSci = fHistList.back();


  std::cout << "trying to read scifi calib" << std::endl;

  for (int scif=0; scif<2; scif++) {
    ifstream readSciFiCalib;
    char scfname[100];
    sprintf(scfname, "calib_link/%s", sciName[scif]);
    readSciFiCalib.open(scfname);

    if (scifiTimePlane[scif] && readSciFiCalib) {
      std::cout << " entering..." << std::endl;
      std::cout << " found for " << scifiTimePlane[scif]->GetName() << std::endl;
      float PosOfWire, PosOfWireBefore;
      int iwire=-1, nummyNum;
      int owire=-2;
      while (owire!=iwire) {
	owire=iwire;
	float scica;
	readSciFiCalib >> iwire >> scica >> nummyNum;
	scifical [scif][iwire] = scica;
	PosOfWire = scifiTimePlane[scif]->GetGeometry()->Wire2Pos(iwire);
	if (iwire==0) pos0[scif] = PosOfWire;
	if (iwire==1) {
	  slop[scif] = 1/(PosOfWire-PosOfWireBefore);
	  std::cout << "Using for channel recovery of "
	       << scifiTimePlane[scif]->GetName()
	       << " (pos-("<<pos0[scif]<<"))*"<<slop[scif]<<std::endl;
	}
	if (iwire>1 &&
	    fabs(iwire-int((PosOfWire-pos0[scif])*slop[scif]+0.5))>1e-8) {
	  std::cout <<"Warning: Nonequidistant wires or not ordered cal file! "
	       <<std::endl;
	}
	std::cout << "identifying ch "<<iwire<<" with pos "
	     << PosOfWire << " and t0 " << scica <<std::endl;
	PosOfWireBefore = PosOfWire;
      }
      nSciFi[scif] = iwire;
    }
    readSciFiCalib.close();
  }

  // Prepare ratio->time conversion
  for (int i=1; i<GEMi; i++) {
    float GEMdt = (GEMcti[i]-GEMcti[i-1]);
    dtdr0[i] = GEMdt/(GEMcr0[i]-GEMcr0[i-1]);
    dtdr1[i] = GEMdt/(GEMcr1[i]-GEMcr1[i-1]);
  }
#endif
}

void GroupGEMTrack::ResetHistograms() {
  Group::ResetHistograms();
#if USE_TRACK == 1
  fTrackMaker->ResetHistos();
#endif
}

void GroupGEMTrack::Tracking() {
#if USE_TRACK == 1

  // reset tracking
  fTrackMaker->Reset();

  // passing clusters to tracking

  for(unsigned i=0; i<fPlanes.size(); i++) {
    if(const GeomPlane *geom = fPlanes[i]->GetGeometry()) {
      const std::set<CCluster1>& clusters = geom->GetClusters();
      for(std::set<CCluster1>::iterator ic=clusters.begin();
	  ic!=clusters.end(); ic++) {
	fTrackMaker->AddClusterWRS(geom->GetID(), ic->pos, ic->size, ic->res);

	// Find out which Tracker is filled: Here
	// we still know cluster timing and can transfer it
	const std::vector<Tracker*>& trs = fTrackMaker->GetTrackers();
	if (ic->dt.size()>4) {
	  for (int ii=0; ii<trs.size(); ii++)
	    if ( geom->GetID() == trs[ii]->GetId() ) {
	      const std::vector<Cluster*>& cli = trs[ii]->GetClusters();
	      cli[cli.size()-1]->fT  = ic->dt[3]; // a1/a2
	      cli[cli.size()-1]->fT2 = ic->dt[4]; // a0/a2
	    }
	}
	else {
	  //	  float thistime=CS::ChipF1::TimeDifference(ic->dt[0], triggerTimeHR);
	  float thistime=ic->dt[0];
	  for (int ii=0; ii<trs.size(); ii++)
	    if ( geom->GetID() == trs[ii]->GetId() ) {
	      const std::vector<Cluster*>& cli = trs[ii]->GetClusters();
	      cli[cli.size()-1]->fT = thistime;
	    }
	  for (int scpl=0; scpl<2; scpl++)
	    if (strcmp(fPlanes[i]->GetName(), sciName[scpl])==0) {
	      int channel = int((ic->pos-pos0[scpl])*slop[scpl]+0.5);
	      hSciTimeVsChan[scpl]->Fill(channel, thistime);
	      hSciTime[scpl]->Fill(thistime-scifical[scpl][channel]);
	    }
	}
      }
    }
  }

  fTrackMaker->Clusterize(maxNumberofClusters);
  fTrackMaker->Candidates();
  fTrackMaker->Residuals();

  const std::vector<Track*>& goodTracks = fTrackMaker->GetAllTracks();

  for (int i=0; i<goodTracks.size(); i++) {
    float lgCh2=log10(goodTracks[i]->GetChi2()+1);
    hChi22->Fill(lgCh2);
    if (goodTracks[i]->GetChi2() < 200) {
      const std::vector<Cluster*>& trClusters = goodTracks[i]->GetClusters();
      float fitime[2] = {-471100, 471100};
      for (int j=0; j<trClusters.size(); j++) {
	if (trClusters[j])
	  for (int scpl=0; scpl<2; scpl++)
	    if ( strcmp(trClusters[j]->GetTracker()->GetName(),
			sciName[scpl]) == 0) {
	      int channel = int((trClusters[j]->fPos-pos0[scpl])*
				slop[scpl]+0.5);
	      fitime[scpl] =  trClusters[j]->fT - scifical[scpl][channel];
	      hSciTime[scpl]->Fill(fitime[scpl]);
	    }
      }
      hdTSci->Fill(fitime[0]-fitime[1]);

      float mft;
      if (scifiTimePlane[0] && scifiTimePlane[1])
	mft = (fitime[0]+fitime[1])/2;
      else {
	if (scifiTimePlane[0]) mft = fitime[0];
	if (scifiTimePlane[1]) mft = fitime[1];
      }

      if (fabs(fitime[0]-fitime[1])<100) {
	
	for (int j=0; j<trClusters.size(); j++) {
	  if (trClusters[j]) {
	    if ( strcmp(trClusters[j]->GetTracker()->GetName(),
			"GM03X1__") == 0) {
	      float fiti_corr = tcsJitter + mft * 0.1289;
	
	      float r0 = trClusters[j]->fT2;
	      float r1 = trClusters[j]->fT;

	      //int i0=0; while (r0<GEMcr0[i0]) i0++;
	      //int i1=0; while (r1<GEMcr1[i1]) i1++;
	
	      //float gm_t0 = (i0==0) ? -4711 :
	      //GEMcti[i0] + dtdr0[i0]*(r0-GEMcr0[i0]);
	      //float gm_t1 = (i1==0) ?  4711 :
	      //GEMcti[i1] + dtdr1[i1]*(r1-GEMcr1[i1]);

	      float rr0, rr1;
	      if      (r0<0.2) rr0 = 0.2;
	      else if (r0>1.1) rr0 = 1.1;
	      else             rr0 = r0;
	      if      (r1<0.3) rr1 = 0.3;
	      else if (r1>1.0) rr1 = 1.0;
	      else             rr1 = r1;

	      float gm_t0 = 30 + 22.2 * log (1.65/rr0-1);   //  14.6....74.0
	      float gm_t1 = 72 + 27.0 * log (1.23/rr1-1);   //  32.3...102.5

	      float gemtime12 = -gm_t1 + fiti_corr;
	      float gemtime02 = -gm_t0 + fiti_corr;

	      hGEMRatio12VsSci->Fill(fiti_corr, r1);
	      hGEMRatio02VsSci->Fill(fiti_corr, r0);
	      hGEMTime12VsSci->Fill(gemtime12, r1);
	      hGEMTime02VsSci->Fill(gemtime02, r0);
	
	      hGEMTime12->Fill(gemtime12);
	      hGEMTime02->Fill(gemtime02);
	      hGEMTimeCorr->Fill(gemtime12, gemtime02);


	    }
	  }
	}
      }
    }
  }

#endif
}

void GroupGEMTrack::EndEvent(const CS::DaqEvent &event) {
#if USE_TRACK == 1

  // Get jitter of TCS clock for APV timing correction
  /*
  float doubles0 = CS::ChipF1::GetTT().GetTimeNorm();
  triggerTimeHR  = CS::ChipF1::GetTT().GetTimeHigh();
  float doubles3=0, doubles4=0;
  float jitter=2000;
  if (!jitterClock) {std::cerr<<"Warning!!! No jitter clock"<<std::endl; return;}
  for (int i=0; i<jitterClock->GetNhits(); i++) {
    int ch=int(jitterClock->GetChanValues()[i]);
    if (ch==30) { //(ch==4||ch==30) { //(ch==29){ //(ch==4||ch==30) {
      float tmp3 = jitterClock->GetTimeValues()[i] + doubles0;
      unsigned short tcs=(unsigned short)tmp3;
      unsigned short trigger=(unsigned short)doubles0;
      if(trigger-tcs<jitter) {
	jitter=trigger-tcs;
	doubles3=tcs;
	doubles4=doubles0-doubles3;
      }
    }
  }
  triggerTime = doubles0;
  tcsJitter=(float((int(doubles4+10707))%200)) * 0.1289;

  hTCSjitter->Fill(tcsJitter);
  */
#endif
}

void GroupGEMTrack::ControlPanel(const TGWindow *p, const TGWindow *main) {
  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}




