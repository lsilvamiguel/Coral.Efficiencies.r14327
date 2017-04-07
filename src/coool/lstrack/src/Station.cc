#include "Station.h"
#include "Tracker.h"

#include <iostream>
#include <strstream>

#include "TBRIK.h"

#define DISPLAY_X_OFFSET 0

ClassImp(Station);

const unsigned int Station::fNdoublets = 2;

Station::Station (const char* name,Doublet* d1, Doublet* d2) : fName(name){    
  fDoublets.push_back(d1);
  fDoublets.push_back(d2);
  
  fPos = (d1->GetPos() + d2->GetPos())/2.;
  fProfile = new TH2F("profile","Profile",50,-20,20,50,-20,20);
}

const std::vector<Point3*>& Station::Coinc3(int plane) {

  const std::vector<Point3*>& p3uv = fDoublets[0]->GetPoints(); 
  //  const std::vector<Point3*>& p3xy = fDoublets[1]->GetPoints(); 
  
  std::vector<Tracker*> planesinxy = fDoublets[1]->GetTrackers();

  for(unsigned i=0; i<p3uv.size(); i++) {
    Track track(p3uv[i]->fX,
		p3uv[i]->fY,
		p3uv[i]->fZ,0,0,0);
    if(planesinxy[0]->CheckCandidate(&track)) {
      fPoints.push_back(new Point3(fPos,p3uv[i]->fY,p3uv[i]->fZ));
      fProfile->Fill(p3uv[i]->fY,p3uv[i]->fZ);
    }
  }
    

//    TClonesArray *p3uv=fDoublets[0]->GetPoints();
//    TClonesArray *p3xy=fDoublets[1]->GetPoints();
  
//    if(p3uv->GetSize()&&p3xy->GetSize()) {
    
//      TIter nextuv(p3uv);
//      while(PointN *pointuv=(PointN*)nextuv()) {
//        Track track((pointuv->GetValues())[0],
//  		  (pointuv->GetValues())[1],
//  		  (pointuv->GetValues())[2],0,0,0);
//        std::vector<Tracker*> pinxy = fDoublets[1]->GetTrackers();
//        if(pinxy[0]->CheckCandidate(&track,1)) {
//  	fPoints.push_back(new Point3(fPos,
//  				     (pointuv->GetValues())[1],
//  				     (pointuv->GetValues())[2]));
//  	fProfile->Fill((pointuv->GetValues())[1],(pointuv->GetValues())[2]);
//        }
      
//      //    if(pinxy[1]->CheckCandidate(&track,5)) {}
//      }
    
//      //    TIter nextxy(p3xy);      
//      //    while(PointN *pointxy=(PointN*)nextxy()) {
//      //      new(tracks[fNcomb++]) Track(point0->GetValues(),point1->GetValues(),0.);
//      //    }
//    }
  return fPoints;
}

void Station::Draw(TNode *worldnode) {

  worldnode->cd();
  TBRIK *mark;
  TNode *node;
  for(unsigned int i=0; i<fPoints.size();i++) {
    mark=new TBRIK("mark","mark","full",.1,.1,.1);
    mark->SetLineColor(2);
    node=new TNode("node","node",mark,fPoints[i]->fX-DISPLAY_X_OFFSET,fPoints[i]->fY,fPoints[i]->fZ);  
  }
}

void Station::Reset() {
  
  for(unsigned i=0; i<fPoints.size(); i++)
    delete fPoints[i];
  
  fPoints.clear();
}
