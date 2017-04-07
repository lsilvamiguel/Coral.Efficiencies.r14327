#include "Doublet.h"

#include <iostream>

#include "TBRIK.h"

#define DISPLAY_X_OFFSET 0

ClassImp(Doublet);

Doublet::Doublet(const char* name, Tracker* p1, Tracker* p2) 
  : fName(name) {
  //, fPoints(new TClonesArray("PointN",10)), fNpoints(0) {

  std::cout<<"Doublet "<<name<<" "<<p1->GetName()<<" and "<<p2->GetName()<<std::endl;
  fTrackers.push_back(p1);
  fTrackers.push_back(p2);

  fPos = (p1->GetX() + p2->GetX())/2.;
}


Doublet::~Doublet() {
  Reset();
}


//  void Doublet::AddTracker(Tracker* plane) {

//    std::cerr<<"DON'T USE THAT !!! "<<std::endl;

//    if(fTrackers.size()<fCnofplanes) {
//      std::cout<<plane->GetName()<<" added to Doublet"<<std::endl;
//      fTrackers.push_back(plane);
//    }
//    else std::cerr<<"error in Doublet::AddTracker : number of planes = "
//  	   <<fCnofplanes<<std::endl;
//  }

const std::vector<Point3*>& Doublet::Point() {

  //Get the clusters and various parameters for each plane
  //Float_t *c0=fTrackers[0]->GetClusters();
  //Float_t *c1=fTrackers[1]->GetClusters();

  const std::vector<Cluster*>& c0 = fTrackers[0]->GetClusters();
  const std::vector<Cluster*>& c1 = fTrackers[1]->GetClusters();

  Float_t sinth0=fTrackers[0]->GetSin();
  Float_t costh0=fTrackers[0]->GetCos();
  Float_t y0=fTrackers[0]->GetY();
  Float_t z0=fTrackers[0]->GetZ();

  Float_t sinth1=fTrackers[1]->GetSin();
  Float_t costh1=fTrackers[1]->GetCos(); 
  Float_t y1=fTrackers[1]->GetY();
  Float_t z1=fTrackers[1]->GetZ();

  //  TClonesArray &points=*fPoints;
  Float_t point[3];
  point[0]=(fTrackers[0]->GetX() + fTrackers[1]->GetX())/2.;

  for(size_t i0=0; i0 < c0.size(); i0++) {
    for(size_t i1=0; i1 < c1.size(); i1++) {

      Float_t a0=c0[i0]->fPos;
      Float_t a1=c1[i1]->fPos;

      point[1]=((z0*sinth0*sinth1+a0*sinth1+sinth1*costh0*y0)-
		(z1*sinth0*sinth1+a1*sinth0+sinth0*costh1*y1))
	/(-sinth0*costh1+sinth1*costh0);
      
      point[2]=0;
      if(sinth0)
	point[2]=z0+a0/sinth0-(costh0/sinth0)*(point[1]-y0);
      else if(sinth1)
	point[2]=z1+a1/sinth1-(costh1/sinth1)*(point[1]-y1);
      else std::cerr<<"Doublet::Point() : 2 X planes in Doublet !"<<std::endl;

      if( InActiveZone(point[1],point[2]) )
	//new(points[fNpoints++]) PointN(3,point);
	fPoints.push_back(new Point3(point[0],point[1],point[2]));
    }
  }
  return fPoints;
}

Bool_t Doublet::InActiveZone(Float_t y, Float_t z) {
 
//    if(z<50 && z>-50) return true;
//    else return false;
  return true;
}


void Doublet::Reset() {

//    for(Int_t i=0;i<fNpoints;i++) 
//      fPoints[i][0]=fPoints[i][1]=0;
//    
  //fNpoints=0;
  for(unsigned i=0; i!=fPoints.size();i++)
    delete fPoints[i];
  fPoints.clear();
}

void Doublet::Draw(TNode *worldnode) {

  worldnode->cd();

  TBRIK *mark;
  TNode *node;
  for(unsigned i=0; i<fPoints.size(); i++) {
    mark=new TBRIK("mark","mark","full",.5,.5,.5);
    mark->SetLineColor(2);
    
    node=new TNode("node","node",mark,fPoints[i]->fX-DISPLAY_X_OFFSET,
		   fPoints[i]->fY,fPoints[i]->fZ); 
  }

  //  TIter next(fPoints);
  //    while(PointN* point=(PointN*)next()){
  //      Float_t *xyz=point->GetValues();
  //      mark=new TBRIK("mark","mark","full",.5,.5,.5);
  //      mark->SetLineColor(2);
  
  //      node=new TNode("node","node",mark,xyz[0]-DISPLAY_X_OFFSET,xyz[1],xyz[2]);
  //    }
}

void Doublet::Dump() {

  std::cout<<"Doublet"<<std::endl;
  for(unsigned int i=0; i<fTrackers.size();i++)
    fTrackers[i]->Dump();
}
