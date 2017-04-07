#include "GroupMuonWallA.h"
#include "TThread.h"
#include "GroupMuonWallAPanel.h"
#include "ChipF1.h"

#include "GeomPlane.h"

ClassImp(GroupMuonWallA);

const int GroupMuonWallA::fGATE = 300;
#define TOTMIN  850
#define MULTMAX  10
#define DIFFXMAX 10

GroupMuonWallA::GroupMuonWallA(const char* name)
  : Group(name) {}

//#if USE_TRACK == 1
//    void GroupMuonWallA::Init() {

//      //  fTrackMaker = new TrackMaker(GetName(),2,0.1);
//      fTrackMaker = new TrackMaker(GetName(),2, 100.0);
//      for(unsigned i=0; i<fPlanes.size(); i++) {
//        const GeomPlane *geom = fPlanes[i]->GetGeometry(); 
//        int id = geom->GetID(); 
//        const char *name = geom->GetName(); 
//        int nwires = geom->GetNWires(); 
//        float x = geom->GetX(); 		    
//        float y = geom->GetY(); 
//        float z = geom->GetZ();
//        float dx = geom->GetDX();
//        float dy = geom->GetDY(); 
//        float dz = geom->GetDZ();
//        float angle = geom->GetAngle();
//    //      float inpitch = .036;
//    //      float outpitch = .042; 
//    //      float ctmin = -500; 
//    //      float ctmax =  500; 
//        float inpitch = 1.;
//        float outpitch = 1.; 
//        float ctmin = -2000; 
//        float ctmax = +2000; 

//        //float ctotmin = 600;
//        float ctotmin = -2000.;
  
//        fTrackMaker->AddTracker(id,name,nwires,x,y,z,dx,dy,dz,angle,inpitch,outpitch,
//    			    ctmin,ctmax,ctotmin);
//      }  
  
//      fTrackMaker->Init();

//      // Adds tracking histogram to this's histogram list 
//      vector<TH1*>& trackhists=fTrackMaker->GetHistograms();
//      for(unsigned i=0; i<trackhists.size();i++) 
//        fHistList.push_back(trackhists[i]);
//    }

//    void GroupMuonWallA::ResetHistograms() {
//      Group::ResetHistograms();

//      if (thr_flag) TThread::Lock();
//      fTrackMaker->ResetHistos();
//      if (thr_flag) TThread::UnLock();
//    }

//  #include "PlaneMuonWallA.h"

//    void GroupMuonWallA::Tracking() {

//      if (thr_flag) TThread::Lock();

//      // reset tracking
//      fTrackMaker->Reset();

//      // passing clusters to tracking
//      for(unsigned i=0; i<fPlanes.size(); i++) {
//        if(const GeomPlane *geom = fPlanes[i]->GetGeometry()) {
//          const set<CCluster1>& clusters = geom->GetClusters();
//          for(set<CCluster1>::iterator ic=clusters.begin(); ic!=clusters.end(); ic++) { 
//    	fTrackMaker->AddClusterWRS(geom->GetID(), ic->pos, ic->size, ic->res);
//          }
//        }
//      }

//      for(unsigned i=0; i<fPlanes.size(); i++) {  
//        if(const GeomPlane *geom = fPlanes[i]->GetGeometry()) {
//          Tracker *ct=fTrackMaker->GetTracker(geom->GetID());
//          if(ct) ct->Activate(0);
//          else {
//    	cerr<<"Cannot find Tracker with ID "<<geom->GetID()<<endl;
//    	exit(1);
//          }
//          // TSV      fTrackMaker->Clusterize(1);
//          fTrackMaker->Clusterize(2);
//          fTrackMaker->Candidates();
//          if(ct) ct->Activate(1);    
//        }
//      }
//      if (thr_flag) TThread::UnLock();
//    }
//#endif

void GroupMuonWallA::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new GroupMuonWallAPanel(p, main, 100, 100, this);
}













