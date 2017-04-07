#include "TrackMaker.h"

#include <iostream>
#include <iomanip>
#include <strstream>
#include <regex.h>

#include "TLine.h"
#include "TString.h"
#include "TText.h"
#include "TF1.h"
#include "TEllipse.h"
#include "TMath.h"
#include "TPad.h"
#include "TProfile.h"

#define DISPLAY_X_OFFSET 0
#define MAXID 10000

ClassImp(TrackMaker);

const unsigned int TrackMaker::fCnofDoublets = 2;


TrackMaker::TrackMaker(const char* name, Int_t mode, float minpchi2) 
  : fName(name),fMinSize(8),fMinpchi2(minpchi2), fMaxchi2(0),
    fR2min(0),fR2max(0),
    fYpmin(0),fYpmax(0),
    fZpmin(0),fZpmax(0),fXmin(0xffff),fXmax(-0xffff),fXRecon(0),
    fChi2Sum(0), fChi2N(0), 
    fEventMask(0), 
    fMode(mode), fNcombOk(0),
    fCurCombSize(0), 
    gGeometry(new TGeometry()) {
  
}

TrackMaker::TrackMaker(const char* name,Int_t mode, float minpchi2,
		       ifstream &in, TFile *file) 
  :
  fName(name),fFile(file),
  fMinSize(8),fMinpchi2(minpchi2),fMaxchi2(0),
  fR2min(0),fR2max(0),fYpmin(0),fYpmax(0),
  fZpmin(0),fZpmax(0),fXmin(0xffff),fXmax(-0xffff),fXRecon(0),
  fChi2Sum(0), fChi2N(0), 
  fEventMask(0),
  fMode(mode),fNcombOk(0),fCurCombSize(0), 
  gGeometry(new TGeometry()) {
    

  if(mode != LSF) {
    std::cerr<<"TrackMaker::TrackMaker : Use mode 2 (Least Squares Fit). Other modes not implemented or not debugged..."<<std::endl;
    exit(0);
  }

  if(!in) {
    std::cerr<<"WARNING in TrackMaker::TrackMaker. cannot read from stream, using default setup"<<std::endl;
    Default();
  }
  else {
    const int linesize=1000;
    char s[linesize];
    int pos=0;
    while(1) { 
      in.seekg(pos);
      in.get(s,linesize);
      if(!in) {
	std::cerr<<"End of file"<<std::endl;
	break;
      }
      pos = (int)in.tellg() + 1;
      
      if(std::string(s).empty()) continue;
      std::istrstream lin(s);  
      
      std::string tag;
      lin>>tag;
      if(tag == "//") continue;
      if(tag == "TM") {
	// TrackMaker options
	
	float ypmin,ypmax,zpmin,zpmax;
	lin>>ypmin>>ypmax>>zpmin>>zpmax;
	if(!lin) 
	  std::cerr<<"Bad TM options"<<std::endl;
	else {
	  SetYpRange(ypmin,ypmax);
	  SetZpRange(zpmin,zpmax);
	}
      }
      else {
	// plane declaration

	lin.seekg(0);
	int active;
	int id, nchan;
	char name[20];
	float x,y,z,angle,ipitch,opitch, ctnsigma; 
	float ctmin,ctmax,ctotmin,crmin,crmax ;
	lin>>active>>id>>name>>nchan>>x>>y>>z>>angle
	   >>ipitch>>opitch
	   >>ctmin>>ctmax>>ctotmin
	   >>crmin>>crmax;
	if(!lin) { 
	  std::cerr<<"Bad Format ! Using default setup"<<std::endl;
	  Default();
	}
	else {
	  Tracker *p=AddTracker(id,name,nchan,x,y,z,1,20,20,angle,ipitch,opitch);
	  p->Activate(active);
	  p->SetCTimeRange(ctmin,ctmax);
	  p->SetCheckRange(crmin,crmax);
	  p->SetCTotMin(ctotmin);
	}
      }
    }
  }
  Init();
}

void TrackMaker::TimeCuts(int n) {
  for(unsigned i=0; i<fTrackers.size(); i++) {
    fTrackers[i]->TimeCuts(n, fFile);
  }
}

void TrackMaker::TimeCuts(float min, float max) {
  for(unsigned i=0; i<fTrackers.size(); i++) {
    fTrackers[i]->SetCTimeRange(min, max);
  }
}

void TrackMaker::Default() {
  std::cerr<<"No default TrackMaker configuration defined. Exiting."<<std::endl;
  exit(0);
}

TrackMaker::~TrackMaker() {
  if(fMode == LSF) {
    delete fConstVec;
  }
}

void TrackMaker::Init() {    

#ifdef __FASTSKIP__
  fastSkip = new FastSkip(fTrackers.size());
#endif
  
  fX = new float[fMaxTracks];
  fY = new float[fMaxTracks];
  fZ = new float[fMaxTracks];
  fYp = new float[fMaxTracks];
  fZp = new float[fMaxTracks];
  fChi2 = new float[fMaxTracks];
  
#ifdef __FASTSKIP__
  fCurClusterInds = new int[fTrackers.size()];
  fCurClusterMult = new int[fTrackers.size()];
#endif
  
  float maxoffset=0;
  float halfsize=0;
  for(unsigned i=0; i<fTrackers.size(); i++) {
    
    // determine zone boundaries 
    if(fTrackers[i]->GetX() < fXmin) fXmin=fTrackers[i]->GetX();
    if(fTrackers[i]->GetX() > fXmax) fXmax=fTrackers[i]->GetX();
    if(fTrackers[i]->GetHalfSize() > halfsize) 
      halfsize = fTrackers[i]->GetHalfSize();
    float offset=fabs(fTrackers[i]->GetY());
    if(offset > maxoffset) 
      maxoffset = offset;
    offset=fabs(fTrackers[i]->GetZ());
    if(offset > maxoffset) 
      maxoffset = offset;
    
				   
    // DEBUG all histograms from all trackers are put in TrackMaker
    std::vector<TH1*>& hists = fTrackers[i]->GetHistograms();
    for(unsigned j=0; j<hists.size(); j++) 
      fHists.push_back(hists[j]);
  }

  // reconstruction plane is taken to be at the center of the zone
  fXRecon = (fXmin+fXmax)/2.;

  // histograms boundaries
  float tmax = maxoffset + halfsize;
  float maxp = 0.05;

  // book histograms
  
  fWorldBrik=new TBRIK("fWorldBrik","world volume","void"
		       ,0.1,0.1,0.1);  
  fWorldNode = new TNode("fWorldNode","world node",fWorldBrik);


  fHzvsy = new TH2F((fName+"_zvsy").c_str(),"z vs y @ x=0",500,-tmax,tmax,
		    500,-tmax,tmax);
  fHists.push_back(fHzvsy);
  fHyp = new TH1F((fName+"_yp").c_str(),"yp  @ x=0",100,-maxp,maxp);
  fHists.push_back(fHyp);
  fHypvsy = new TH2F((fName+"_ypvsy").c_str(),"yp vs y @ x=0",100,-tmax,tmax,
		     100,-maxp,maxp);
  fHists.push_back(fHypvsy);
  fHzp = new TH1F((fName+"_zp").c_str(),"zp",100,-maxp,maxp);
  fHists.push_back(fHzp);
  fHzpvsz = new TH2F((fName+"_zpvsz").c_str(),"zp vs z @ x=0",
		     100,-tmax,tmax,100,-maxp,maxp);
  fHists.push_back(fHzpvsz);
  fHchi2 = new TH1F((fName+"_chi2").c_str(),"Chi2",1000,0,200);  
  fHists.push_back(fHchi2);
  fHlogchi2 = new TH1F((fName+"_logChi2").c_str(),"log10(Chi2)",480,-4,8);  
  fHists.push_back(fHlogchi2);
  fHprobchi2 = new TH1F((fName+"_probchi2").c_str(),"Chi2 probablity",
			500,0,1);  
  fHists.push_back(fHprobchi2);
  
  fHNtracks = new TH1F((fName+"_ntracks").c_str(),"Tracks per event",
		       100,-0.5,99.5);  
  fHists.push_back(fHNtracks);
  
  // allocates space for vector of constants
  fConstVec = new TVector(4);
}

void TrackMaker::OutputTree(TTree *tree) {
  
  tree->Branch("ntracks",&fNcombOk,"ntracks/I",32000);
  tree->SetBranchStatus("ntracks",1);

  tree->Branch("x",fX,"x[ntracks]/F",32000);
  tree->SetBranchStatus("x",1);

  tree->Branch("y",fY,"y[ntracks]/F",32000);
  tree->SetBranchStatus("y",1);

  tree->Branch("z",fZ,"z[ntracks]/F",32000);
  tree->SetBranchStatus("z",1);

  tree->Branch("yp",fYp,"yp[ntracks]/F",32000);
  tree->SetBranchStatus("yp",1);

  tree->Branch("zp",fZp,"zp[ntracks]/F",32000);
  tree->SetBranchStatus("zp",1);

  tree->Branch("chi2",fChi2,"chi2[ntracks]/F",32000);  
  tree->SetBranchStatus("chi2",1);

  tree->Branch("mask",&fEventMask,"mask/I",32000);
  tree->SetBranchStatus("mask",1);
}


Tracker* TrackMaker::AddTracker(Int_t id,const char *name, 
				Int_t nwires, Float_t x, 
				Float_t y, Float_t z,
				float dx, float dy, float dz,
				Float_t angle, Float_t inpitch, Float_t outpitch) {
  
  //fTrackers[id] = new Tracker(name,nwires,x,y,z,angle,pitch);
  //fCheckTrackers[id] = fTrackers[id];
  Tracker *p=new Tracker(id,name,nwires,x,y,z,dx,dy,dz,angle,inpitch,outpitch);
  AddTracker(p);
  return p;
}

Tracker* TrackMaker::AddTracker(Int_t id,const char *name, 
				Int_t nwires, Float_t x, 
				Float_t y, Float_t z, 
				float dx, float dy, float dz,
				Float_t angle, Float_t inpitch, 
				Float_t outpitch, float ctmin, float ctmax, 
				float ctotmin) {
  
  Tracker *p = AddTracker(id,name,nwires,x,y,z,dx,dy,dz,angle,inpitch,outpitch);
  p->SetCTimeRange(ctmin,ctmax);
  p->SetCTotMin(ctotmin);
  return p;
}

void TrackMaker::AddTracker(const Tracker *tracker) {
  if(!tracker) return;
  Tracker *mycopy = new Tracker(*tracker);
  AddTracker(mycopy);
}

void TrackMaker::AddTracker(Tracker *tracker) {
  
  int id = tracker->GetId();
  if(id > MAXID) {
    std::cerr<<"detector ID is bigger than "<<MAXID<<"... sorry."<<std::endl;
    return;
  }
  else if(id<0) {
    std::cerr<<"negative detector ID ... sorry."<<std::endl;
    return;
  }
  else if(id<fTrackersID.size())
    fTrackersID[id] = tracker;
  else {
    if(id>fTrackersID.size())
      fTrackersID.resize(id,0);
    fTrackersID.push_back(tracker);
  }

  fTrackers.push_back(tracker);
  fCheckTrackers.push_back(tracker);  
}

Doublet* TrackMaker::CreateDoublet(const char* name, unsigned int id1, unsigned int id2) {

  Doublet *d=0;

  // Looking for both planes int the plane vector
  Tracker *p1=0;
  Tracker *p2=0;
  for(unsigned int i=0; i<fTrackers.size(); i++) {
    if(fTrackers[i]->GetId()==id1)
      p1=fTrackers[i];
    else if(fTrackers[i]->GetId()==id2)
      p2=fTrackers[i];
  }
  
  if(p1 && p2) {
    d=new Doublet(name,p1,p2);
    fDoublets.push_back(d);
  }
  else {
    std::cerr<<"Cannot create doublet : "<<id1<<" or "<<id2
	<<" do not exist in :"<<std::endl;
    for(unsigned int i=0; i<fTrackers.size(); i++) 
      std::cout<<fTrackers[i]->GetId()<<"\t";
    std::cout<<std::endl;
  }
  return d;
}

Station* TrackMaker::CreateStation(const char* name, Doublet* d1, Doublet *d2) {
  Station *s=new Station(name,d1,d2);
  fStations.push_back(s);
  return s;
}

void TrackMaker::ResetHistos() {

  fChi2Sum=0;
  fChi2N=0;

  // resets histograms for all planes
  for(unsigned i=0; i<fTrackers.size(); i++) 
    fTrackers[i]->ResetHistos();
  
  for(unsigned i=0; i<fHists.size(); i++) 
    fHists[i]->Reset();
}

void TrackMaker::Reset() {

  //resetting combinations
  fNcombOk=0;
  fNcomb=0;
  fCurComb.clear();

  //resetting track arrays
  for(unsigned int i=0;i<fTracks.size();i++) 
    delete fTracks[i];
  fTracks.clear();
  fTracksOk.clear();
  
  for(unsigned int i=0;i<fTrackers.size();i++) 
    fTrackers[i]->ResetHits();
  for(unsigned int i=0;i<fDoublets.size();i++) 
    fDoublets[i]->Reset();
  for(unsigned int i=0;i<fStations.size();i++) 
    fStations[i]->Reset();



//    for (iLp=fTrackers.begin(); iLp!=fTrackers.end(); iLp++) {
//      (*iLp).second->ResetHits();
//    }  
//    for (iLs=fDoublets.begin(); iLs!=fDoublets.end(); iLs++) {
//      (*iLs).second->Reset();
//    } 
}


void TrackMaker::AddDigit(unsigned int id, Int_t channel, int t, int tot) {
  
  //OPTIMIZE
//    for(unsigned int i=0;i<fTrackers.size();i++) {
//      if(id==fTrackers[i]->GetId()) {
//        fTrackers[i]->AddDigit(channel,t,tot);
//        break;
//      }
//      if(i==fTrackers.size())
//        std::cerr<<"No detector with ID "<<id<<std::endl;
//    } 
  
  if(!fTrackersID[id])
    std::cerr<<"No detector with ID "<<id<<std::endl;
  else 
    fTrackersID[id]->AddDigit(channel,t,tot);
}

void TrackMaker::AddClusterWRS(unsigned int id, float pos, int size, float res) {
  
  //OPTIMIZE put this stuff in a map !!!
  unsigned int i=0;
  for(;i<fTrackers.size();i++) {
    if(id==fTrackers[i]->GetId()) {
      fTrackers[i]->AddClusterWRS(pos,size,res);
      break;
    }
  } 
//    if(i==fTrackers.size())
//      std::cerr<<"No detector with ID "<<id<<std::endl;
}

void TrackMaker::AddClusterWRS(unsigned int id, float pos, float time,
			       int size, float res) {
  
  //OPTIMIZE put this stuff in a map !!!
  unsigned int i=0;
  for(;i<fTrackers.size();i++) {
    if(id==fTrackers[i]->GetId()) {
      fTrackers[i]->AddClusterWRS(pos,time,size,res);
      break;
    }
  } 
//    if(i==fTrackers.size())
//      std::cerr<<"No detector with ID "<<id<<std::endl;
}

void TrackMaker::AddCluster(unsigned int id, float channel, float t, float tot,
			    float size) {
  
  //OPTIMIZE put this stuff in a map !!!
  for(unsigned int i=0;i<fTrackers.size();i++) {
    if(id==fTrackers[i]->GetId()) {
      fTrackers[i]->AddCluster(channel,t,tot,size);
      break;
    }
//      if(i==fTrackers.size())
//        std::cerr<<"No detector with ID "<<id<<std::endl;
  } 
}


Int_t TrackMaker::Clusterize(int nclustperplane) {

  fBadEvent=false;

  int nclusters=0;

  // clusterize each plane
  for(unsigned int i=0;i<fTrackers.size();i++) {
    int nc=fTrackers[i]->Clusterize();
 
    if(nclustperplane && nc > nclustperplane && fTrackers[i]->IsActive())
      fBadEvent=true; 
    
    nclusters+=nc;
  } 

  return nclusters;
}

Int_t TrackMaker::Candidates() {

  if(fBadEvent) return 0;

  switch(fMode) {
    
  case SP: 
    {
      const std::vector<Point3*>& p1s=fStations[0]->Coinc3(0);
      const std::vector<Point3*>& p3s=fStations[1]->Coinc3(0);
      
      //TClonesArray& tracks=*fTracks;
      
      if(p1s.empty() || p3s.empty()) return 0;
      else {
	
	for(unsigned i1=0;i1<p1s.size();i1++) {
	  for(unsigned i3=0;i3<p3s.size();i3++) {
	    Track *t=new Track(p1s[i1]->fX,p1s[i1]->fY,p1s[i1]->fZ,
			       p3s[i3]->fX,p3s[i3]->fY,p3s[i3]->fZ,0.);
	    fTracks.push_back(t);
	    fHyp->Fill(t->GetYp());
	    fHzp->Fill(t->GetZp());
	    t->Move(0.);
	    fHypvsy->Fill(t->GetY(),t->GetYp());
	    fHzpvsz->Fill(t->GetZ(),t->GetZp());	
	  }
	} 
      }
      break;
    }
    
  case LSF : 
    {

#ifdef __FASTSKIP__
      int combMult=1;
      for (int i=0; i<fTrackers.size(); i++){
	if ( fTrackers[i]->IsActive() ) {
	  fCurClusterMult[i] = fTrackers[i]->GetClusters().size();
	  // std::cout << fCurClusterMult[i]<<"*";
	  combMult *= fCurClusterMult[i];
	  const std::vector<Cluster*>& clusters = fTrackers[i]->GetClusters();
	  if (clusters.size()>0) {
	    double pos = clusters[0]->fPos;
	    bool isOK=true;
	    for (int j=1; j<clusters.size(); j++){
	      double npos = clusters[j]->fPos;
	      if (npos<pos) isOK=false;
	      pos=npos;
	    }
	    if (!isOK)
	      std::cout << " Cluster pos not increasing for tracker "
		   << fTrackers[i]->GetName() << std::endl;
	  }
	}
	else fCurClusterMult[i]=0;
      }
      //      std::cout << " = Combination multiplicity: " << combMult << std::endl;
      
      fastSkip->Reset(fCurClusterMult);
#endif

      NextTracker();
      fHNtracks->Fill(fTracks.size());
      
      if (fTracks.size()<2)
	Efficiency();
      
      // fastSkip->PrintStatistics();
      
    }
  }
  return fNcombOk;
}

int TrackMaker::NextTracker() {

  static int iplane = -1;

  //std::cout <<iplane<<"#";
  
  iplane++;
  
  if(static_cast<unsigned>(iplane) == fTrackers.size()) {
    
#ifdef __FASTSKIP__
    if (!fastSkip->CheckCluster(fCurClusterInds))
#endif
      {
	//std::cout<<"Cluster pattern:";
	//for (int i=0; i<fTrackers.size(); i++) {
	//  std::cout<< fCurClusterInds[i]<<","<<fCurComb[i];
	//  if (fCurComb[i]) std::cout << ","<<fCurComb[i]->GetTracker()->GetName();
	//  std::cout<<"|";
	//}
	//std::cout << std::endl;
	LsfTrack();
      }
    
    iplane--;
    // std::cout<<"last plane : "<<fTrackers[iplane]->GetName()<<std::endl;
    return 0;
  }
  
  const std::vector<Cluster*>& clusters = fTrackers[iplane]->GetClusters();
  if(clusters.empty() || !(fTrackers[iplane]->IsActive())) {
    fCurComb.push_back(0x0);    
#ifdef __FASTSKIP__
    fCurClusterInds[iplane]=0;
#endif
    int n=NextTracker();
    //    std::cout<<"plane "<<fTrackers[iplane]->GetName()<<" empty "<<n<<std::endl;
    fCurComb.pop_back();
    iplane--;
    return n;
  }
  
  for(unsigned ic=0; ic<clusters.size(); ic++) {
    fCurComb.push_back(clusters[ic]);
    if (fTrackers[iplane]->IsActive()) fCurCombSize++;
#ifdef __FASTSKIP__
    fCurClusterInds[iplane]=ic+1;
#endif
    NextTracker();
    //      if(!NextTracker()) {
    //        LsfTrack();
    //      }
    fCurComb.pop_back();
    if (fTrackers[iplane]->IsActive()) fCurCombSize--;
  } 
  
  //std::cout<<fTrackers[iplane]->GetName()<<" : "<<clusters.size()<<" clusters processed"<<std::endl;
  iplane--;
  return clusters.size();
}

void TrackMaker::LsfTrack() {

  //int config=0;
  //std::cout<<"Combination : "<<fCurComb.size()<<" clusters   ";
  //for(unsigned i=0; i<fCurComb.size(); i++) {
  //  if(fCurComb[i]) {
      // this stuff is not yet used. 
      //config = config | (fCurComb[i]->fZone << i);
      //std::cout<<fCurComb[i]->fPos<<" ";
  //  }
    //  else 
      //std::cout<<"X ";
  //}
  //std::cout<<std::endl;

  // Get cov matrix associated to this combination 
  // std::cout<<"Associated CovMatrix label : "<<hex<<config<<std::endl;

  MakeTrack();
}


void TrackMaker::MakeTrack() {
  
  fConstVec->Zero();

  if(fCurCombSize<fMinSize) return;
  fNcomb++;

  for(unsigned i=0; i<fCurComb.size(); i++) {
    if(fCurComb[i] && fCurComb[i]->GetTracker()->IsActive()) {
      Tracker *plane = fCurComb[i]->GetTracker();//fTracker;
      float sigma2 = fCurComb[i]->fRes*fCurComb[i]->fRes; 
                     // plane->GetPitch()*plane->GetPitch()/12.;
      float sm = (plane->GetA() + fCurComb[i]->fPos)/sigma2;
      (*fConstVec)(0) += sm*plane->GetCos();
      (*fConstVec)(1) += sm*plane->GetSin();
      (*fConstVec)(2) += sm*plane->GetXc();
      (*fConstVec)(3) += sm*plane->GetXs();
    }
  }
  
  //  static TMatrix invmatrix(4,4);
  //  static bool firstTime=true;
  //  if (firstTime) {
  TMatrix matrix(4,4);
  for(unsigned i=0; i<fCurComb.size(); i++) {
    if(fCurComb[i] && fCurComb[i]->GetTracker()->IsActive()) {
      Tracker *plane = fCurComb[i]->GetTracker();//fTracker;
      float sigma2 = fCurComb[i]->fRes*fCurComb[i]->fRes; //plane->GetPitch()*plane->GetPitch()/12.;
      matrix(0,0) += plane->GetC2()/sigma2;
      matrix(0,1) += plane->GetSc()/sigma2;
      matrix(0,2) += plane->GetXc2()/sigma2;
      matrix(0,3) += plane->GetXsc()/sigma2;
      
      matrix(1,0) += plane->GetSc()/sigma2;
      matrix(1,1) += plane->GetS2()/sigma2;
      matrix(1,2) += plane->GetXsc()/sigma2;
      matrix(1,3) += plane->GetXs2()/sigma2;
      
      matrix(2,0) += plane->GetXc2()/sigma2;
      matrix(2,1) += plane->GetXsc()/sigma2;
      matrix(2,2) += plane->GetX2c2()/sigma2;
      matrix(2,3) += plane->GetX2sc()/sigma2;
      
      matrix(3,0) += plane->GetXsc()/sigma2;
      matrix(3,1) += plane->GetXs2()/sigma2;
      matrix(3,2) += plane->GetX2sc()/sigma2;
      matrix(3,3) += plane->GetX2s2()/sigma2;
    }
  }
  
  static bool firstTime=true;
  if (firstTime) {
    std::cout<<"Summed element matrix:" <<std::endl;
    matrix.Print();
    TMatrix iinvmatrix(4,4);
    iinvmatrix = matrix.Invert();
    std::cout<<"Inverted matrix:"<<std::endl;
    iinvmatrix.Print();
    
    for(unsigned i=0; i<fCurComb.size(); i++) {
      if(fCurComb[i] && fCurComb[i]->GetTracker()->IsActive()) {
	Tracker *plane = fCurComb[i]->GetTracker();//fTracker;
	std::cout << "Tracker "<<plane->GetName()
	     <<" at "<<plane->GetX()<<" cm"<<std::endl;
	
	// Get track parameters at detector position
	double x=plane->GetX();
	double trans_elements[4][4] = {{1,0,x,0},
				       {0,1,0,x},
				       {0,0,1,0},
				       {0,0,0,1} };
	TMatrix transp_translate(4,4);
	TMatrix translate(4,4);
	for (int ii=0; ii<4; ii++)
	  for (int jj=0; jj<4; jj++) {
	    translate       (ii,jj)=trans_elements[ii][jj];
	    transp_translate(jj,ii)=trans_elements[ii][jj];
	  }
	//translate.Print();
	//transp_translate.Print();
	TMatrix transCov = translate*iinvmatrix*transp_translate;
	transCov.Print();
      }
    }
    firstTime=false;
  }
  
  TMatrix invmatrix(TMatrix::kInverted,matrix);

  (*fConstVec)*=invmatrix;

  // invmatrix.Print();

  //  (fConstVec->operator*=)(*fCovMat);

  float y0 = (*fConstVec)(0);
  float z0 = (*fConstVec)(1);
  float tany = (*fConstVec)(2);
  float tanz = (*fConstVec)(3);
  //  std::cout << "Used matrix, result:" << y0<<","<<tany<<std::endl;

  // chi2 
  float chi2 = 0;
  for(unsigned i=0; i<fCurComb.size(); i++) {
    if(fCurComb[i] && fCurComb[i]->GetTracker()->IsActive()) {
      Tracker *plane = fCurComb[i]->GetTracker();//fTracker;
      float sigma2 = fCurComb[i]->fRes*fCurComb[i]->fRes;//plane->GetPitch()*plane->GetPitch()/12.;
      float sqnum = fCurComb[i]->fPos - plane->GetCos()*(y0+
							 tany*plane->GetX()) 
	- plane->GetSin()*(z0+tanz*plane->GetX()) + plane->GetA();      
      chi2 += sqnum*sqnum/sigma2;
      
#ifdef __FASTSKIP__
      int aClInd = abs(fCurClusterInds[i]);
      if (sqnum<0) fCurClusterInds[i] = - aClInd;
      else         fCurClusterInds[i] =   aClInd;
#endif

    }
  }
  int ndf = fCurCombSize-4;
  chi2 = chi2/static_cast<float>(ndf);
  float probchi2 = TMath::Prob(chi2,ndf);
  fHchi2->Fill(chi2);
  fHlogchi2->Fill(log10(chi2));
  fHprobchi2->Fill(probchi2);
  
  // cut on chi2
  if ( (fMinpchi2 && probchi2<fMinpchi2) ||
       (fMaxchi2 && chi2>fMaxchi2)          ) {
#ifdef __FASTSKIP__
    fastSkip->AddSkip(fCurClusterInds);
#endif
    return;
  }
  
  Track *track = new Track(0,y0,z0,tany,tanz,0);
  
  fChi2Sum += chi2;
  fChi2N++;

  // histograms filling 

  track->Move(fXRecon);
  float y=track->GetY();
  float z=track->GetZ();

  fHzvsy -> Fill(y,z);
  fHyp -> Fill(tany);
  fHzp -> Fill(tanz);
  fHypvsy -> Fill(y,tany);
  fHzpvsz -> Fill(z,tanz);

  
  // surface cuts
  float r2=y*y+z*z;
  if(fR2min && fR2max && (r2<fR2min || r2>fR2max)) 
    return; 
  
  //angle cuts 
  if(fYpmin && fYpmax && (tany<fYpmin || tany>fYpmax))
    return; 
  if(fZpmin && fZpmax && (tanz<fZpmin || tanz>fZpmax))
    return;

  // Track is OK
  
  //std::cout << "*** Good track ***" << std::endl;
  //std::cout << y<<","<< tany<<","<< z<<","<< tanz<<", tgrackers:"<< std::endl;
  //for(unsigned i=0; i<fCurComb.size(); i++) {
  //  if(fCurComb[i] && fCurComb[i]->GetTracker()->IsActive()) {
  //    std::cout << fCurComb[i]->GetTracker()->GetName();
  //  }
  // }
  //std::cout << std::endl;
  
  track->SetChi2(chi2);
  track->SetChi2Prob(probchi2);

  // tree variables
  if(fNcombOk<fMaxTracks) {
    fX[fNcombOk] = fXRecon;
    fY[fNcombOk] = y;
    fZ[fNcombOk] = z;
    fYp[fNcombOk] = tany;
    fZp[fNcombOk] = tanz;
    fChi2[fNcombOk] = chi2;
  }
  fNcombOk++;
  
  std::vector<Cluster*> fEffCurComb = fCurComb;  
  
  // Efficiency
  for(unsigned i=0; i<fTrackers.size(); i++) {
    if(!fTrackers[i]->IsActive()) {
      int EffClusterID = fTrackers[i]->Efficiency(track, 0);
      //std::cout << "EffClusterID is" << EffClusterID << std::endl;
      if(EffClusterID<1)
  	fEventMask=0;
      else {
  	fEventMask=1;
  	fEffCurComb[i] = fTrackers[i]->GetClusters()[EffClusterID-1];
      }
    }
  }

  track->SetClusters(fEffCurComb);
  
  fTracks.push_back(track); 

}

void TrackMaker::Efficiency() {
  for (unsigned j=0; j<fTracks.size(); j++) {
    for(unsigned i=0; i<fTrackers.size(); i++) {
      if(!fTrackers[i]->IsActive()) {
	fTrackers[i]->Efficiency(fTracks[j], 1);
      }
    }
  }
}

Int_t TrackMaker::CheckCandidates(Float_t checkwidth, unsigned int okthr) {

  
  if(okthr > fCheckTrackers.size()) {
    std::cout<<"TrackMaker::CheckCandidates : asking for "<<okthr
	<<"compatible planes, but number of planes for checking = "
	<<fCheckTrackers.size()<<std::endl;
    return 0;
  } 

  //TClonesArray& tracksok=*fTracksOk;
  //TIter nt(fTracks);
  //  map<Int_t, Tracker*>::iterator np;

  //Int_t fNcombOk=0;  

  for(unsigned j=0; j<fTracks.size(); j++) {
    unsigned int nok=0;
 
    for(unsigned int i=0; i<fCheckTrackers.size();i++) {
      if(fCheckTrackers[i]->CheckCandidate(fTracks[j])>0)
	nok++;
    }

    if (nok>=okthr) {
      fTracksOk.push_back(fTracks[j]);
    }    
  }

//    while(Track* track=(Track*) nt()) {

//      unsigned int nok=0;

//      for(unsigned int i=0; i<fCheckTrackers.size();i++) {

//        Float_t xp = fCheckTrackers[i]->GetX();
//        track->Move(xp);
//        if(fCheckTrackers[i]->CheckCandidate(track->GetY(),track->GetZ(),checkwidth))
//  	nok++;
//      }
//      if (nok>=okthr) {
//        tracksok[fNcombOk++]=track;
//      }
//    }
  return fNcombOk;
}


void TrackMaker::Residuals() {

  for(unsigned i=0; i<fTracks.size(); i++) {
    for(unsigned j=0; j<fTrackers.size(); j++) {
      fTracks[i]->Move(fTrackers[j]->GetX());
      fTrackers[j]->Residuals(fTracks[i]->GetY(),fTracks[i]->GetZ());
    }
  }
}

void TrackMaker::PlotProbchi2() {
  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  fHprobchi2->Draw();
  TLine *l = new TLine(fMinpchi2,0,fMinpchi2,fHprobchi2 ->GetMaximum());
  l->Draw();
}


void TrackMaker::PlotYp() {
  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  fHyp->Draw();
  TLine *l = new TLine(fYpmin,0,fYpmin, fHyp->GetMaximum());
  l->Draw();
  TLine *l2 = new TLine(fYpmax,0,fYpmax, fHyp->GetMaximum());
  l2->Draw();
}

void TrackMaker::PlotZp() {
  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  fHzp->Draw();
  TLine *l = new TLine(fZpmin,0,fZpmin, fHzp->GetMaximum());
  l->Draw();
  TLine *l2 = new TLine(fZpmax,0,fZpmax, fHzp->GetMaximum());
  l2->Draw();
}

void TrackMaker::PlotResVsPos() {

  if(fTrackers.size()<6) return;

  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();

  mpad->Divide(3,2);

  for(unsigned i=0; i<fTrackers.size(); i++) {
    mpad->cd(i+1); 
    TProfile *prof = fTrackers[i]->fHresVSpos;
    prof->Fit("pol1");
    TF1* func=prof->GetFunction("pol1");
    func->SetLineColor(5);
    float p0 = func->GetParameter(0);
    float p1 = func->GetParameter(1);
    char text[50];
    sprintf(text,"%f%s%f", p1," *X + ",p0);
    
    float xmin = prof->GetXaxis()->GetXmin();
    float xmax = prof->GetXaxis()->GetXmax();
    float ymin = prof->GetYaxis()->GetXmin();
    float ymax = prof->GetYaxis()->GetXmax();
    TText *t = new TText(xmin + (xmax-xmin)/10., ymin + (ymax-ymin)/10.,
			 text);
    t->Draw();
  }
  mpad->cd();
}


void TrackMaker::PlotRes() {
  if(fTrackers.size()<6) return;

  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  mpad->Divide(3,2);
  for(unsigned i=0; i<fTrackers.size(); i++) {
    mpad->cd(i+1); 
    fTrackers[i]->PlotRes();
  }  
  mpad->cd();
}

void TrackMaker::PlotCTimes() {
  if(fTrackers.size()<6) return;

  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  mpad->Divide(4,3);
  for(unsigned i=0; i<fTrackers.size(); i++) {
    mpad->cd(i+1); 
    fTrackers[i]->PlotCTime();
  }  
  mpad->cd();
}

void TrackMaker::PlotProfile() {
  if(fTrackers.size()<6) return;

  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  fHzvsy->Draw();
  TEllipse *inner=new TEllipse(0,0,sqrt(fR2min));
  TEllipse *outer=new TEllipse(0,0,sqrt(fR2max));
  
  inner->SetLineColor(4);
  outer->SetLineColor(4);
  inner->Draw();
  outer->Draw();
}


void TrackMaker::Draw(bool planes, bool stations, bool doublets, bool tracks) {

  gGeometry->RecursiveRemove(fWorldNode);
  fWorldNode = new TNode("fWorldNode","world node",fWorldBrik); 

  if(planes)
    for(unsigned int i=0;i<fTrackers.size();i++) 
      fTrackers[i]->Draw(fWorldNode);

  
  if(doublets)
    for(unsigned int i=0;i<fDoublets.size();i++) 
      fDoublets[i]->Draw(fWorldNode);
 
  if(stations)
    for(unsigned int i=0;i<fStations.size();i++) 
      fStations[i]->Draw(fWorldNode);
 
  fWorldNode->Draw();
  fWorldNode->cd();

  if(tracks) {
    for(unsigned int i=0;i<fTracks.size();i++) 
      fTracks[i]->Draw(fXmin-100,fXmax+100,DISPLAY_X_OFFSET);
  }

}

std::vector<Tracker*>&  TrackMaker::GetTrackers(const char* selection) {
  
  fTrackersSel.clear();
  for( unsigned int iInf = 0; iInf < fTrackers.size(); iInf++ ) {
    if(string_match(fTrackers[iInf]->GetName(),selection))
      fTrackersSel.push_back( fTrackers[iInf] );  
  }
  
  return fTrackersSel;
}

void TrackMaker::OutConfig(const char *filename) {
  
  ostream *out;

  if(!strcmp("",filename))
    out = &std::cout;
  else 
    out = new ofstream(filename);
  
  for(unsigned i=0; i<fTrackers.size();i++) {
    (*out)<<std::setw(4)<<(int)fTrackers[i]->IsActive()<<"\t"
	  <<std::setw(4)<<fTrackers[i]->GetId()<<"\t"
	  <<std::setw(10)<<fTrackers[i]->GetName()<<"\t"
	  <<std::setw(5)<<fTrackers[i]->GetNwires()<<"\t"
	  <<std::setw(10)<<fTrackers[i]->GetX()<<"\t"
	  <<std::setw(10)<<fTrackers[i]->GetY()<<"\t"
	  <<std::setw(10)<<fTrackers[i]->GetZ()<<"\t"
	  <<std::setw(5)<<fTrackers[i]->GetAngle()<<"\t"
	  <<std::setw(10)<<fTrackers[i]->GetIp()<<"\t"
	  <<std::setw(10)<<fTrackers[i]->GetOp()<<"\t"
	  <<std::setw(5)<<fTrackers[i]->GetTMin()<<"\t"
	  <<std::setw(5)<<fTrackers[i]->GetTMax()<<"\t"
	  <<fTrackers[i]->GetTotMin()<<"\t"
	  <<fTrackers[i]->GetCheckMin()<<"\t"
	  <<fTrackers[i]->GetCheckMax()<<std::endl;
  }
  (*out)<<std::endl;
  (*out)<<"TM"<<"\t"<<fYpmin<<"\t"<<fYpmax<<"\t"
	<<fZpmin<<"\t"<<fZpmax<<std::endl;
  
}

void TrackMaker::Dump() {
  
  std::cout<<"**********        TrackMaker "<<fName<<"        ************"<<std::endl;

  std::cout<<"Cuts_________________________________________"<<std::endl;
  std::cout<<"\tmin combination size : "<<fMinSize<<std::endl;
  std::cout<<"\tmax chi2         : "<<fMaxchi2<<std::endl;
  std::cout<<"\tmin chi2 prob    : "<<fMinpchi2<<std::endl;
  std::cout<<"\tmin/max  yp      : "<<fYpmin<<" "<<fYpmax<<std::endl;
  std::cout<<"\tmin/max  zp      : "<<fZpmin<<" "<<fZpmax<<std::endl;

  std::cout<<"Trackers_____________________________________"<<std::endl;
  for(unsigned int i=0;i<fTrackers.size();i++) 
    fTrackers[i]->Dump();
}

bool TrackMaker::string_match(const char *str, const char* pattern) {
  int    status;
  regex_t    re;
  
  if( regcomp(&re, pattern, REG_EXTENDED|REG_NOSUB) != 0 )
    return false;
  
  status = regexec(&re, str, (size_t) 0, NULL, 0);
  regfree(&re);
  
  if (status != 0)
    return false;
  
  return true;
}

#ifdef __FASTSKIP__

FastSkip::FastSkip(int ntr) {
  ntrackers = ntr;
  max_index = new int[ntrackers];
  multiplik = new int[ntrackers];
  jmin      = new short[ntrackers];
  jmax      = new short[ntrackers];
  fastCount = new int[5];
#ifdef __FASTSKIP_TIMER__
  fBenchLoopTracker = new TStopwatch();
  fBenchLoopTracker->Stop();
  fBenchCheckCluster = new TStopwatch();
  fBenchCheckCluster->Stop();
#endif
  shMultiplik = new int[ntrackers];
  skipShort = NULL;
}

void
FastSkip::Reset(int *cluMult){
  combMult=1;
  lat=0;  // last active tracker
  for (int i=ntrackers-1; i>=0; i--) {
    max_index[i] = cluMult[i];    
    if (!lat) if (cluMult[i]!=0) lat=i;
    multiplik[i] = combMult;
    if (cluMult[i]) {
      combMult*=cluMult[i];
    }
  }

  shMult=1;
  for (int i=lat-1; i>=0; i--) {
    shMultiplik[i] = shMult;
    if (cluMult[i]) shMult = shMult*cluMult[i] + 1;
  }
  delete skipShort;
  skipShort = new unsigned short[shMult*2];
  for (int i=0; i<shMult; i++) {
    skipShort[2*i  ] = 0;
    skipShort[2*i+1] = max_index[lat]+1;  
  }
  
  for (int i=0; i<5; i++) fastCount[i]=0;

#ifdef __FASTSKIP_TIMER__
  fBenchLoopTracker->Reset();
  fBenchLoopTracker->Stop();
  fBenchCheckCluster->Reset();
  fBenchCheckCluster->Stop();
#endif
}

void
FastSkip::PrintStatistics(){
  std::cout << "AddSkip:                  " << fastCount[0]<<std::endl;
  std::cout << "Checked Clusters:         " << fastCount[1]<<std::endl;
  std::cout << "Last tracker comparisons: " << fastCount[2]<<std::endl;
  std::cout << "change table updated:     " << fastCount[3]<<std::endl;
#ifdef __FASTSKIP_TIMER__
  std::cout << "time in LoopTracker: "<<fBenchLoopTracker->RealTime()<<std::endl;
  std::cout << "time in CheckCluster: "<<fBenchCheckCluster->RealTime()<<std::endl;
#endif
}

bool 
FastSkip::CheckCluster(int *cl) {
#ifdef __FASTSKIP_TIMER__
  fBenchCheckCluster->Start(false);
#endif
  fastCount[1]++;

  int this_index=0;

  for (int i=0; i<lat; i++) {
    int c=abs(cl[i]);
    if (0)  std::cout << c << "|";
    if (c) this_index += (c-1) * shMultiplik[i] + 1;
  }

  bool valShort = !(abs(cl[lat])>skipShort[2*this_index] &&
		    abs(cl[lat])<skipShort[2*this_index+1]);

#ifdef __FASTSKIP_TIMER__
  fBenchCheckCluster->Stop();
#endif
  return valShort;
};


bool
FastSkip::ShortLT(int ti, int i, int *cl) {
  
  if (cl[i]==0) {
    if (i<lat-1) return ShortLT(ti, i+1, cl);
    else return true;
  }
  
  short cll=cl[lat];

  bool goChange = !((cll<0 && -cll<=skipShort[2*ti  ]) ||
		    (cll>0 &&  cll>=skipShort[2*ti+1]));

  if (!goChange) return false;

  bool haveChange=false;
  for (int j=jmin[i]; j<=jmax[i]; j++) {
    bool jChange=false;
    if (skipPast) {
      if ( j < -cl[i] ) continue;
      if ( j >  cl[i] ) skipPast=false;
    }
    int ni = ti + (j-1) * shMultiplik[i] + 1;
    if (i==lat-1) {
      if (cl[i+1]<0 && -cl[i+1] > skipShort[2*ni  ])
	{skipShort[2*ni  ] = -cl[i+1]; jChange=true;}
      if (cl[i+1]>0 &&  cl[i+1] < skipShort[2*ni+1]) 
	{skipShort[2*ni+1] =  cl[i+1]; jChange=true;}
      fastCount[2]++;
    }
    else {
      jChange=ShortLT(ni, i+1, cl);
    }
    if (jChange) haveChange=true;  
  }
  
  if (!skipPast && haveChange) {
    int jjmin=max_index[i];
    int jjmax=0;
    for (int j=1; j<=max_index[i]; j++) {
      int ni = ti + (j-1) * shMultiplik[i] + 1;
      if (jjmin>skipShort[2*ni  ]) jjmin=skipShort[2*ni  ];
      if (jjmax<skipShort[2*ni+1]) jjmax=skipShort[2*ni+1];
    }
    skipShort[2*ti  ] = jjmin;
    skipShort[2*ti+1] = jjmax;
    
    if (haveChange!=goChange) std::cout << "Tracker "<<i<< " Alarm!! for cell " 
				   << ti << std::endl;
    fastCount[3]++;
  }
  
  return haveChange;
}

void
FastSkip::AddSkip(int *cl) {
  //
  if (0) {
    std::cout <<"*** ----------------->      Adding: |";
    for (int i=0; i<lat; i++) std::cout << cl[i] << "|";
    std::cout<<cl[lat]<<std::endl;
  }
  //
#ifdef __FASTSKIP_TIMER__
  fBenchLoopTracker->Start(false);
#endif
  fastCount[0]++;
  
  for (int i=0; i<lat+1; i++) {
    jmin[i] = (cl[i]>0) ?        cl[i] :          1; 
    jmax[i] = (cl[i]>0) ? max_index[i] : abs(cl[i]);
  }

  skipPast=true;
  ShortLT(0, 0, cl);
  
#ifdef __FASTSKIP_TIMER__
  fBenchLoopTracker->Stop();
#endif
}

#endif





