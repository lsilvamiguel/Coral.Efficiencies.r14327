#include "GroupMwpc.h"
#include "TThread.h"
//#include "GroupMwpcPanel.h"
#include "ChipF1.h"

#include <stdio.h>
#include <fstream>
#include <vector>
#include "PlaneMwpc.h"
#include "GroupMwpc.h"

#include "GeomPlane.h"
#include "MwpcReconst.h"
#include "TThread.h"
#include "MwpcEventDisplay.h"

using namespace std;


ClassImp(GroupMwpc);

const int GroupMwpc::fGATE = 300;

#define TOTMIN  850
#define MULTMAX  10
#define DIFFXMAX 10

CChamber* chamber[11] =
{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

static PlaneGeometry* config[4] =
{
  &XPlane, &YPlane, &UPlane, &VPlane
};



void GroupMwpc::Init() {
  const char* name = GetName();
  if(strncmp(name,"ST_",3) == 0)
    {
      sscanf(name,"ST_%02d",&fNstation);
      if(fNstation >= 1 && fNstation <= 11)
	{
	  if(fNstation == 1) // A* chamber
	    {
	      config[Y_PLANE] = &YPlane;
	    }
	  else
	    {
	      config[Y_PLANE] = NULL;
	    }

	  char histname[32];
	  chamber[fNstation-1] = new CChamber(config);

	  sprintf(histname,"%s_X_vs_Y",name);
	  fHxvsy = new TH2F(histname, name, 145,-77,77,105,-52,52);
	  AddHistogram(fHxvsy);
	  fHxvsy->GetXaxis()->SetTitle("x-coordinate [cm]");
	  fHxvsy->GetYaxis()->SetTitle("y-coordinate [cm]");

	  sprintf(histname,"%s_error",name);
	  fHrecErr = new TH1F(histname, "Reconstruction error (cm)", 40,0,4);
	  AddHistogram(fHrecErr);

	  sprintf(histname,"%s_evwin",name);
	  fEvDisplay = new EventDisplay(histname);
	  AddHistogram(fEvDisplay);

	  sprintf(histname,"%s_nrec",name);
	  fHnRecEv = new TH1F(histname, "Number of reconstructed tracks", 10,0,10);
	  AddHistogram(fHnRecEv);

	  sprintf(histname,"%s_X_distr",name);
	  fHx = new TH1F(histname, "X distribution", 160,-80,80);
	  AddHistogram(fHx);
	  fHx->GetXaxis()->SetTitle("x-coordinate [cm]");

	  sprintf(histname,"%s_Y_distr",name);
	  fHy = new TH1F(histname, "Y distribution", 233,-60.0,70.0);
	  AddHistogram(fHy);
	  fHy->GetXaxis()->SetTitle("y-coordinate [cm]");

	  sprintf(histname,"%s_XPLANE_EFF",name);
	  fHxEff = new TH1F(histname, "X-plane efficiency", 100,50.0,100.0);
	  AddHistogram(fHxEff);
	  fHxEff->GetXaxis()->SetTitle("Efficiency (%)");

	  sprintf(histname,"%s_UPLANE_EFF",name);
	  fHuEff = new TH1F(histname, "U-plane efficiency", 100,50.0,100.0);
	  AddHistogram(fHuEff);
	  fHuEff->GetXaxis()->SetTitle("Efficiency (%)");

	  sprintf(histname,"%s_VPLANE_EFF",name);
	  fHvEff = new TH1F(histname, "V-plane efficiency", 100,50.0,100.0);
	  AddHistogram(fHvEff);
	  fHvEff->GetXaxis()->SetTitle("Efficiency (%)");
	}
      else
	fNstation = -1;
    }
  else
    fNstation = -1;

#if USE_TRACK == 1
  fTrackMaker = new TrackMaker(GetName(),2,0.1);
  if (fDoTracking) {
  for(unsigned i=0; i<fPlanes.size(); i++) {
    const GeomPlane *geom = fPlanes[i]->GetGeometry();
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
    fTrackMaker->AddTracker(id,name,nwires,x,y,z,dx,dy,dz,angle,inpitch,inpitch,
			    ctmin,ctmax,0);
    }
  }

  fTrackMaker->Init();

  // Adds tracking histogram to this's histogram list
  vector<TH1*>& trackhists=fTrackMaker->GetHistograms();
  for(unsigned i=0; i<trackhists.size();i++)
    fHistList.push_back(trackhists[i]);
#endif
}

void GroupMwpc::ResetHistograms() {
  Group::ResetHistograms();
#if USE_TRACK == 1
  fTrackMaker->ResetHistos();
#endif
}

void GroupMwpc::Tracking() {
  if (thr_flag) TThread::Lock();

#if USE_TRACK == 1
  // lstrack tracking

  // reset tracking
  fTrackMaker->Reset();

  // passing clusters to tracking
  for(unsigned i=0; i<fPlanes.size(); i++) {
    if(const GeomPlane *geom = fPlanes[i]->GetGeometry()) {
      const set<CCluster1>& clusters = geom->GetClusters();
      for(set<CCluster1>::iterator ic=clusters.begin(); ic!=clusters.end(); ic++)
	fTrackMaker->AddClusterWRS(geom->GetID(), ic->pos, ic->size, ic->res);
    }
  }

  fTrackMaker->Clusterize(2);
  fTrackMaker->Candidates();

#else
  // Andrea's tracking
  fNevent++;
  if(fNstation != -1) {
    MwpcEvent* reconstEvent = chamber[fNstation-1]->GetReconPoints();
    fHnRecEv->Fill(reconstEvent->fNpoints);
    if(fNevent%999 == 0) {
      if(chamber[fNstation-1]->CalculateEfficiency()) {
	fHxEff->Fill(chamber[fNstation-1]->GetEfficiency(X_PLANE)*100.);
	fHuEff->Fill(chamber[fNstation-1]->GetEfficiency(U_PLANE)*100.);
	fHvEff->Fill(chamber[fNstation-1]->GetEfficiency(V_PLANE)*100.);
      }
    }
    if(reconstEvent->fNpoints) {
      ((EventDisplay*)fEvDisplay)->Fill(chamber[fNstation-1]->GetHits());
      for(int i = 0; i < reconstEvent->fNpoints; i++) {
	fHxvsy->Fill(reconstEvent->fX[i],reconstEvent->fY[i]);
	fHrecErr->Fill(reconstEvent->fErr[i]);
	fHx->Fill(reconstEvent->fX[i]);
	fHy->Fill(reconstEvent->fY[i]);
      }
    }
    chamber[fNstation-1]->ClearHits();
  }


#endif // USE_TRACK == 1
  if (thr_flag) TThread::UnLock();
}




