#include "Tracker.h"

#include <iostream>
#include <strstream>

#include "Fits.h"
#include "TLine.h"
#include "TText.h"
#include "TF1.h"
#include "TEllipse.h"
#include "TMath.h"
#include "TPad.h"
#include "TBRIK.h"

#define DISPLAY_X_OFFSET 0

#define WIRE_DY 0.02
#define WIRE_DX 1. 

ClassImp(Tracker);

const float Tracker::sqrt12 = sqrt(12);

Tracker::Tracker(int id,const char *name, Int_t nwires, float x, 
		 float y, float z, float dx, float dy, float dz,
		 float angle, float inpitch, float outpitch) 
  : fClustering(true),fNeff(0), fNtot(0), fNbackg(0), fIsActive(true),
    fNsigmat(0),
    fId(id), fName(name), fX(x), fY(y), fZ(z),
    fYbck(y), fZbck(z),
    fDx(dx), fDy(dy), fDz(dz), fAngle(angle), fR2Min(3*3), fNwires(nwires),
    fPitch((inpitch+outpitch)/2.), fIp(inpitch), fOp(outpitch),
    fPrevChan(-1), fCTMin(-65000), fCTMax(65000), fCTotMin(0), fHeffcor(0){

  // a few parameters --------------------------------------------

  fPitchs[0] = fIp;
  fHalfsize = fNwires*fIp/2.;

  Init();

  std::cout<<"Tracker : "<<fName<<" created "<<id<<" "<<nwires<<" "
      <<x<<" "<<y<<" "<<z<<" "<<angle<<" "<<fIp<<" "<<fOp<<std::endl;
}


Tracker::Tracker(const Tracker& t) 
  : fClustering(t.fClustering),fNeff(t.fNeff), fNtot(t.fNtot), 
  fNbackg(t.fNbackg), fIsActive(t.fIsActive), fNsigmat(t.fNsigmat),
  fId(t.fId), fName(t.fName), fX(t.fX), fY(t.fY), fZ(t.fZ),
  fYbck(t.fYbck), fZbck(t.fZbck),  
  fDx(t.fDx), fDy(t.fDy), fDz(t.fDz), fAngle(t.fAngle),fHalfsize(t.fHalfsize), 
  fR2Min(t.fR2Min), fNwires(t.fNwires),
  fPitch(t.fPitch), fIp(t.fIp), fOp(t.fOp), fPitchs(t.fPitchs),
  fPrevChan(t.fPrevChan), fCTMin(t.fCTMin), fCTMax(t.fCTMax), 
  fCTotMin(t.fCTotMin), fHeffcor(t.fHeffcor)
{
  Init();

  std::cout<<"Tracker : "<<fName<<" created "<<fId<<" "<<fNwires<<" "
      <<fX<<" "<<fY<<" "<<fZ<<" "<<fAngle<<" "<<fIp<<" "<<fOp<<std::endl;
}

void Tracker::Init() {
  fCheckWidth = 9*GetPitch()/sqrt12;
  fCheckMin = -fCheckWidth/2.;
  fCheckMax = fCheckWidth/2.;

  TrackingParameters();

  // histograms ------------------------------------------------
  char ht[50];
  char hn[20];

  sprintf(hn,"%s%s",fName.c_str(),"_hcpos");
  sprintf(ht,"%s %s",fName.c_str()," clusters distribution");
  fHcpos=new TH1F(hn,ht,100,-fHalfsize,fHalfsize);
  fHists.push_back(fHcpos);

  sprintf(hn,"%s%s",fName.c_str(),"_hdcpos");
  sprintf(ht,"%s %s",fName.c_str()," Distance between Clusters");
  fHdcpos=new TH1F(hn,ht,500,-1,1);
  fHists.push_back(fHdcpos);

  sprintf(hn,"%s%s",fName.c_str(),"_hres");
  sprintf(ht,"%s %s",fName.c_str()," Residuals");
  fHres=new TH1F(hn,ht,100,-0.2,0.2);
  fHists.push_back(fHres);

  sprintf(hn,"%s%s",fName.c_str(),"_hresl");
  sprintf(ht,"%s %s",fName.c_str()," Residuals large");
  fHresl=new TH1F(hn,ht,100,-4,4);
  fHists.push_back(fHresl);
  
  sprintf(hn,"%s%s",fName.c_str(),"_hresz");
  sprintf(ht,"%s %s",fName.c_str()," Residuals zoom");
  fHresz=new TH1F(hn,ht,200,-0.05,0.05);
  fHists.push_back(fHresz);
  
  sprintf(hn,"%s%s",fName.c_str(),"_hrespos");
  sprintf(ht,"%s %s",fName.c_str()," Residuals vs. Position");
  fHrespos=new TH2F(hn,ht,200,-0.1,0.1,50,-3,3);
  fHists.push_back(fHrespos);
  
  sprintf(hn,"%s%s",fName.c_str(),"_hresVSpos");
  sprintf(ht,"%s %s",fName.c_str()," Residuals");
  fHresVSpos=new TProfile(hn,ht,50,-20,20);
  fHresVSpos->SetMinimum(-0.5);
  fHresVSpos->SetMaximum(0.5);
  fHists.push_back(fHresVSpos);
  
  sprintf(hn,"%s%s",fName.c_str(),"_hctime");
  sprintf(ht,"%s %s",fName.c_str()," Cluster Time");
  fHctime=new TH1F(hn,ht,100,-1000,1000);
  fHists.push_back(fHctime);
  
  sprintf(hn,"%s%s",fName.c_str(),"_hctot");
  sprintf(ht,"%s %s",fName.c_str()," Cluster Tot");
  fHctot=new TH1F(hn,ht,100,0,4000);
  fHists.push_back(fHctot);
  
  sprintf(hn,"%s%s",fName.c_str(),"_hchits");
  sprintf(ht,"%s %s",fName.c_str()," Number of Clusters");
  fHchits=new TH1F(hn,ht,25,0,25);
  fHists.push_back(fHchits);
  
  sprintf(hn,"%s%s",fName.c_str(),"_heff");
  sprintf(ht,"%s %s",fName.c_str()," Seen Tracks");
  fHeff=new TH2F(hn,ht,50,-20,20,50,-20,20);
  fHists.push_back(fHeff);

  sprintf(hn,"%s%s",fName.c_str(),"_heffcor");
  sprintf(ht,"%s %s",fName.c_str()," Corrected Efficiency Profile");
  fHeffcor=new TH2F(hn,ht,50,-20,20,50,-20,20);
  fHists.push_back(fHeffcor);
  fHeffcor->SetStats(kFALSE);

  sprintf(hn,"%s%s",fName.c_str(),"_tcks");
  sprintf(ht,"%s %s",fName.c_str()," All Tracks");
  fHtot=new TH2F(hn,ht,50,-20,20,50,-20,20);
  fHists.push_back(fHtot);

  sprintf(hn,"%s%s",fName.c_str(),"_bg");
  sprintf(ht,"%s %s",fName.c_str()," Background");
  fHbackg=new TH2F(hn,ht,50,-20,20,50,-20,20);
  fHists.push_back(fHbackg);

  //
  // 1-d fine resolution efficiency histograms
  //
  
  sprintf(hn,"%s%s",fName.c_str(),"_h1deff");
  sprintf(ht,"%s %s",fName.c_str()," 1D Seen Tracks");
  fH1deff=new TH1F(hn, ht, 1000, -20, 20);
  fHists.push_back(fH1deff);
  
  sprintf(hn,"%s%s",fName.c_str(),"_h1deffcor");
  sprintf(ht,"%s %s",fName.c_str()," 1D Corrected Efficiency Profile");
  fH1deffcor=new TH1F(hn, ht, 1000, -20, 20);
  fHists.push_back(fH1deffcor);
  fH1deffcor->SetStats(kFALSE);
  
  sprintf(hn,"%s%s",fName.c_str(),"_1dtcks");
  sprintf(ht,"%s %s",fName.c_str()," 1D All Tracks");
  fH1dtot=new TH1F(hn, ht, 1000, -20, 20);
  fHists.push_back(fH1dtot);
  
  sprintf(hn,"%s%s",fName.c_str(),"_1dbg");
  sprintf(ht,"%s %s",fName.c_str()," 1D Background");
  fH1dbackg=new TH1F(hn, ht, 1000, -20, 20);
  fHists.push_back(fH1dbackg);
  
  sprintf(hn,"%s%s",fName.c_str(),"_tres");
  sprintf(ht,"%s %s",fName.c_str()," Time Resolution");
  fHtres=new TH1F(hn,ht,60,-60,60);
  fHists.push_back(fHtres);
   
  fHitPat =  new int[fNwires/32+1];
  fHitToTs =  new int[fNwires];
  fHitTimes =  new int[fNwires];
  fNwireshit=0;

  fRotMatrix=new TRotMatrix("fRotMatrix","Rotation matrix",90,0,90-fAngle,90,
			    fAngle,-90);
  
}

Tracker::~Tracker() {
  ResetHits();
}

void Tracker::SetX(float x) {
  fX=x;
  TrackingParameters();
}

void Tracker::SetY(float y) {
  fY=y;
  TrackingParameters();
}

void Tracker::SetZ(float z) {
  fZ=z;
  TrackingParameters();
}

void Tracker::MovePerp(float offset) {
  fYbck = fY;
  fZbck = fZ;

  fY += offset * GetCos();
  fZ += offset * GetSin();
  TrackingParameters();
}

void Tracker::MoveBack() {
  fY = fYbck;
  fZ = fZbck;
  TrackingParameters();
}

void Tracker::TrackingParameters() {

  fSin=sin(fAngle*TMath::Pi()/180.);
  fCos=cos(fAngle*TMath::Pi()/180.);
  
  fXc = fX*fCos;
  fC2 = fCos*fCos;
  fXc2 = fX*fC2;
  fX2c2 = fX*fXc2;

  fXs = fX*fSin;
  fS2 = fSin*fSin;
  fXs2 = fX*fS2;
  fX2s2 = fX*fXs2;
  
  fSc = fSin*fCos;
  fXsc = fX*fSin*fCos;
  fX2sc = fX*fXsc;

  fA = fZ*fSin + fY*fCos;
  //fA = fY;
  fAc = fA*fCos;
  fXac = fX*fAc;
  fAs = fA*fSin;
  fXas = fX*fAs;
}

void Tracker::AddSubDetector(int id, const char *name, Int_t nwires, float x, 
			     float y, float z, float dx, float dy, float dz,
			     float angle, float inpitch, float outpitch) {

  typedef std::map<int,double>::iterator IP;
  std::cout<<"Adding subdetector"<<std::endl;
  
  if(fName!=name)
    throw "Tracker::AddSubDetector : different names !";
  if(fabs(angle - fAngle) > 0.001)
    throw "Tracker::AddSubDetector : different angles !";
  
  // CONTINUITY TEST (in WRS) 
  double tolerance=fPitchs[0]/10.;


  // old half size of the active zone (WRS)
//    double halfxs = 0;
//    int lastw=0;
//    double lastp=0.;
//    for(IP i=fPitchs.begin(); i!=fPitchs.end(); i++) {
//      halfxs += (i->first-lastw)*lastp;
//      lastw = i->first;
//      lastp = i->second;
//    }
//    halfxs += (fNwires - lastw)*lastp; 
//    halfxs = halfxs/2.;

  // vector from old detector center to subdetector center in WRS
  float cscy =  GetCos()*(y-fY) + GetSin()*(z-fZ);
  float cscz = -GetSin()*(y-fY) + GetCos()*(z-fZ);

  double dist=sqrt((fabs(cscy)-(fHalfsize+nwires*inpitch/2.))
		   *(fabs(cscy)-(fHalfsize+nwires*inpitch/2.))
		   + cscz*cscz);  

  if(fabs(dist)>tolerance)
    throw "Tracker::AddSubDetector : Continuity test failed";

  // storing subdet pitch, calculating center pos in WRS

  float cmwrsy =  GetCos()*fY + GetSin()*fZ;

  if(cscy>0) {
    fPitchs[fNwires] = inpitch;
    cmwrsy += nwires*inpitch/2.;
  }
  else {
    std::map<int, double> newpitch;
    newpitch[0] = inpitch;
    
    for(IP i=fPitchs.begin(); i!=fPitchs.end(); i++) {
      newpitch[i->first + nwires] = i->second;
    }
    fPitchs = newpitch;
    cmwrsy -= nwires*inpitch/2.;
  }
  fNwires += nwires;

  //new size
  fDy += dy;
  fHalfsize += nwires*inpitch/2.;

  // new center pos in MRS
  float cmwrsz = -GetSin()*fY + GetCos()*fZ;
  fY = GetCos()*cmwrsy - GetSin()*cmwrsz;
  fZ = GetSin()*cmwrsy + GetCos()*cmwrsz;
} 
 

void Tracker::TimeCuts(int n, TFile *file) {

//    string histname("MM/");
//    histname += GetName();
//    histname += "/";
//    histname += GetName();    
//    histname += "_ct";
//    TH1F *h=(TH1F*)file->Get(histname.c_str()); 
//    if(h) {
//      float sigma; float mean;
//      FitCTime(h,  sigma, mean);
//      //SetCTimeRange(mean-n*sigma,mean+n*sigma);
//      SetCTimeRange(0,500);
//      fNsigmat = n;
//    }
//    else {
//      cerr<<"Tracker::TimeCuts() Cannot find histogram "<<histname;
//    }
//    std::cout<<"Time Cuts : "<<GetTMin()<<" "<<GetTMax()<<std::endl;
  SetCTimeRange(-500,500);
}

void Tracker::Activate(bool active) {
  fIsActive=active;
}

void Tracker::ResetHistos() {

  fNtot=0;
  fNeff=0;
  fNbackg=0;
  

  for(unsigned i=0; i<fHists.size(); i++)
    fHists[i]->Reset();
}

void Tracker::AddDigit(Int_t chan, int time, int tot) {
  
  //const int fHIT_TMIN = -10000-1000;
  //const int fHIT_TMAX = -10000+1500;
  const int fMT_T0    = -10000;
  
  //if (fHIT_TMIN<time && time<fHIT_TMAX) {
  if (chan!=fPrevChan) {
    // NEW WIRE => SET PATTERN for CLUSTER SEARCH and tot and Time 
    fPrevChan = chan;
    fHitToTs[chan] = tot; fHitTimes[chan] = time;
    fHitPat[chan/32] |= 1<<chan%32;
  }
  else if (fabs(time-fMT_T0)<fabs(fHitTimes[chan]-fMT_T0)) {
    // RETAIN HIT CLOSEST TO T0 => UPDATE ToT and Time
    fHitToTs[chan] = tot; fHitTimes[chan] = time;
  } 
}

void Tracker::AddClusterWRS(float pos, int size, float res) {

  //  fClustering = false;
  fClustering = true; //jaf 17.8.03

  fHcpos->Fill(pos);
  //fHctime->Fill(time);
  //fHctot->Fill(tot);

  //if(IsActive()) {
  //  if(time<fCTMin || time > fCTMax) return;
  //  if(tot<fCTotMin) return;
  //}

  fClusters.push_back(new Cluster(pos, 0, 0, size, res, this));
}

void Tracker::AddClusterWRS(float pos, float time, int size, float res) {

  // fClustering = false;
  fClustering = true;  //jaf 17.8.03

  fHcpos->Fill(pos);
  fHctime->Fill(time);
  //fHctot->Fill(tot);

  //if(IsActive()) {
  //  if(time<fCTMin || time > fCTMax) return;
  //  if(tot<fCTotMin) return;
  //}

  fClusters.push_back(new Cluster(pos, time, 0, size, res, this));
}


void Tracker::AddCluster(float chan, float time, float tot, float size) {

  fClustering = false;

  float pos = C2x(chan);

  fHcpos->Fill(pos);
  fHctime->Fill(time);
  fHctot->Fill(tot);

  if(IsActive()) {
    if(time<fCTMin || time > fCTMax) return;
    if(tot<fCTotMin) return;
  }

  float res = 0;
  //int zone = 0;
  if(255.5<chan && chan<767.5) {
    res = fIp/sqrt12;
  }
  else {
    res = fOp/sqrt12;
    //zone = 1;
  }
  fClusters.push_back(new Cluster(pos, time, tot, size, res,this));
}

void Tracker::ResetHits() {
  
  fNwireshit=0;  
  for(int i=0; i<fNwires/32+1; i++)
    fHitPat[i]=0;

  for(unsigned i=0; i!=fClusters.size();i++)
    delete fClusters[i];
  fClusters.clear();
  
} 


int Tracker::Clusterize() {

  if(!fClustering) {
    // Plane::AddCluster was used
    
    std::vector<Cluster*> ctmp;

    for(unsigned i=0; i<fClusters.size(); i++) {
      if(!i) 
	ctmp.push_back(fClusters[i]); 
      else { 
	if(ctmp.back()->fPos != fClusters[i]->fPos) 
	  ctmp.push_back(fClusters[i]);
	else {
	  // multiple cluster (2 clusters at the same pos)
	  delete fClusters[i];
	}
      }
    }
    fClusters = ctmp;
    //    return fClusters.size();
  }
  fHchits->Fill(fClusters.size());
  return fClusters.size();  //jaf 17.8.03
  // return 0;
}


Int_t Tracker::CheckCandidate(float y,float z,
			      float checkmin, float checkmax,
			      bool fillHistos) {

  if(!InActiveZone(y,z)) return -1;
  
  y=y-fY;  //referential associated to this plane
  z=z-fZ;

  //distance to the center of the "wire" passing through (y,z)
  float a=GetCos()*y + GetSin()*z;
  float b=-GetSin()*y + GetCos()*z;
  
  //find a cluster close enough in this event
  
  if (fillHistos) {
    fNtot++;
    fHtot->Fill(a,b);
    fH1dtot->Fill(a);
  }
  
  float maxbg=2*checkmax - checkmin;
  float minbg=checkmax;
  int found = 0;
  bool foundbg = false;
  for(unsigned i=0;i<fClusters.size();i++) {

    if (fillHistos) {
      fHres->Fill(fClusters[i]->fPos-a);
      fHresl->Fill(fClusters[i]->fPos-a);
      fHresz->Fill(fClusters[i]->fPos-a);
      fHrespos->Fill(fClusters[i]->fPos-a, fClusters[i]->fPos);      
      fHresVSpos->Fill(fClusters[i]->fPos, fClusters[i]->fPos-a);      
    }

    if(!found &&
       fClusters[i]->fPos-a < checkmax && 
       fClusters[i]->fPos-a > checkmin) { 
      found=i+1;
      if (fillHistos) {
	fHeff->Fill(a,b);
	fH1deff->Fill(a);
	fNeff++;   
      }
    } 
    if(!foundbg && 
       fClusters[i]->fPos-a < maxbg && 
       fClusters[i]->fPos-a > minbg) {
      foundbg=true;
      if (fillHistos) {
	fHbackg->Fill(a,b);
	fH1dbackg->Fill(a);
	fNbackg++;
      }
    }
  }
  return found;
}

Int_t Tracker::CheckCandidate(Track *track,float checkmin, float checkmax,
			      bool fillHistos) {
  
  track->Move(fX);
  return CheckCandidate(track->GetY(),track->GetZ(),checkmin,checkmax,
			fillHistos);
}

Int_t Tracker::CheckCandidate(Track *track, bool fillHistos) {
  
  track->Move(fX);
  return CheckCandidate(track->GetY(),track->GetZ(),fCheckMin,fCheckMax,
			fillHistos);
}

void Tracker::PlotEffProfile() {
  
  if(IsActive()) {
    std::cerr<<"WARNING : this tracker is used in tracking !"<<std::endl;
  }

  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());
  
  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  EffProfile()->Draw("colz");
  WriteEfficiency();
}

TH2F* Tracker::EffProfile() {
  TH1F denom1d(*fH1dtot);
  denom1d.Add(fH1dbackg,-1);
  fH1deffcor->Reset();
  fH1deffcor->Add(fH1deff,1);
  fH1deffcor->Add(fH1dbackg,-1);
  fH1deffcor->Divide(&denom1d);
  fH1deffcor->SetMinimum(0);
  
  TH2F denom(*fHtot);
  denom.Add(fHbackg,-1);
  fHeffcor->Reset();
  fHeffcor->Add(fHeff,1);
  fHeffcor->Add(fHbackg,-1);
  fHeffcor->Divide(&denom);
  fHeffcor->SetMinimum(0);
  
  return fHeffcor;
}

void Tracker::WriteEfficiency() {
  char eff[5];
  sprintf(eff,"%3.1f",GetEfficiency()*100);
  TText *teff = new TText(15,-15,eff);
  teff->Draw();
}

float Tracker::GetEfficiency() const {
  if (fNtot && fNbackg<fNtot) {
    float ntotnb=static_cast<float>(fNtot-fNbackg);
    return (fNeff-fNbackg)/ntotnb;
  }
  else return -1;   
}

float Tracker::GetBGEfficiency() const {
  if (fNtot) {
    return fNbackg/static_cast<float>(fNtot);
  }
  else return -1;
}

float Tracker::GetDEfficiency() const {
  if(fNbackg<fNeff) {
    float ea = fNeff/static_cast<float>(fNtot);
    float pb = fNbackg/static_cast<float>(fNtot);

    return 
      (sqrt(ea*(1-ea))/(ea-pb) + 
       sqrt(pb*(1-pb))*(1-ea)/((ea-pb)*(1-pb)))/sqrt(fNtot);
  }
  else return -1;
}


int Tracker::Efficiency(Track* track, bool fillHistos) {
  
  if(fIsActive) 
    std::cerr<<"WARNING ! (Tracker::Efficiency) "<<GetName()<<" active in tracking"<<std::endl;
  
  return CheckCandidate(track, fCheckMin, fCheckMax, fillHistos); 
}

bool Tracker::InActiveZone(Track* track) {
 
  track->Move(fX);
  return InActiveZone(track->GetY(),track->GetZ());
}

bool Tracker::InActiveZone(float y,float z) {
  y=y-fY;  //referential associated to this plane
  z=z-fZ;
  
  float a=GetCos()*y + GetSin()*z;  
  float b=-GetSin()*y + GetCos()*z;  
  
  if((a*a+b*b) < fR2Min || a>fDy || b>fDz)
    return false;
  else 
    return true;
}

float Tracker::GetResidualOfCluster(float pos, float y, float z) {
  y=y-fY;
  z=z-fZ;
  float a= GetCos()*y + GetSin()*z;
  //float b=-GetSin()*y + GetCos()*z;
  return pos-a;
}

void Tracker::Residuals(float y,float z) {

  y=y-fY;  //referential associated to this plane
  z=z-fZ;
  
  //distance to the center of the "wire" passing through (y,z)
  float a= GetCos()*y + GetSin()*z;
  float b=-GetSin()*y + GetCos()*z;
//    for(Int_t i=0;i<fNclusters;i++) {
//      fResiduals->Fill(fClAbscs[i]-a);
//    }  
  for(size_t i=0;i<fClusters.size();i++) {
    fHres->Fill(fClusters[i]->fPos-a);
    fHresl->Fill(fClusters[i]->fPos-a);
    fHresz->Fill(fClusters[i]->fPos-a);
    fHrespos->Fill(fClusters[i]->fPos-a,b);    
    fHresVSpos->Fill(fClusters[i]->fPos,fClusters[i]->fPos-a);    
  }  
}


void Tracker::Dump() {

  std::cout<<"Tracker "<<fName<<" x: "<<fX<<" y: "<<fY<<" z: "<<fZ
      <<" ang: "<<fAngle<<" pitch: "<<fPitch<<std::endl;

  if(!fClusters.empty()) {
    std::cout<<fClusters.size()<<"\tClusters found : ";
    //    for(Int_t icl=0;icl<fNclusters;icl++){
    //      std::cout<<fClAbscs[icl]<<"\t";
    //    }
    for(size_t i=0;i<fClusters.size();i++) {
      std::cout<<fClusters[i]->fPos<<"\t";
    } 
    std::cout <<std::endl;
  }
}


void Tracker::Draw(TNode *worldnode) {


  worldnode->cd();
  TBRIK *det = new TBRIK("det","det","void",fDx,fDy,fDz);
  TNode *detnode = new TNode("detnode","detnode",det,fX-DISPLAY_X_OFFSET,fY,fZ,fRotMatrix);
  detnode->cd();

  TBRIK *wire;
  TNode *node;
  for(size_t i=0;i<fClusters.size();i++) {

    wire = new TBRIK("BRIK","BRIK","full",WIRE_DX,fPitch,fDz);
    wire -> SetLineColor(4);
//      node = new TNode("NODE","NODE",wire,0.,C2x(fClAbscs[i]),0.);
    node = new TNode("NODE","NODE",wire,0.,fClusters[i]->fPos,0.);
  }
}

void Tracker::PlotRes() {
  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  
  fHres->Draw();
  TLine *l = new TLine(fCheckMin,0,fCheckMin, fHres->GetMaximum());
  l->Draw();
  TLine *l2 = new TLine(fCheckMax,0,fCheckMax, fHres->GetMaximum());
  l2->Draw();
}

void Tracker::DrawResCuts() {

  TLine *l = new TLine(fCheckMin,0,fCheckMin, fHres->GetMaximum());
  l->Draw();
  TLine *l2 = new TLine(fCheckMax ,0,fCheckMax, fHres->GetMaximum());
  l2->Draw();
}

void Tracker::PlotCTime() {
  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  
  fHctime->Draw();
  TLine *l = new TLine(fCTMin,0,fCTMin, fHctime->GetMaximum());
  l->Draw();
  TLine *l2 = new TLine(fCTMax,0,fCTMax, fHctime->GetMaximum());
  l2->Draw();
}

void Tracker::PlotCTot() {
  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  
  fHctot->Draw();
  TLine *l = new TLine(fCTotMin,0,fCTotMin, fHctot->GetMaximum());
  l->Draw();
}

void Tracker::PlotCPos() {
  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    std::cerr<<"create a pad first..."<<std::endl;
    return;
  }
  mpad->Clear();
  
  fHcpos->Draw();
  TLine *l = new TLine(-sqrt(fR2Min) + GetOffset(),0,-sqrt(fR2Min) + GetOffset(), fHcpos->GetMaximum());
  l->Draw();
  TLine *l2 = new TLine(sqrt(fR2Min) + GetOffset(),0,sqrt(fR2Min) + GetOffset(), fHcpos->GetMaximum());
  l2->Draw();
}

float Tracker::C2x(float wire)  // Channel number to Position X
{
  //  float x;

//    if(fOp) { // new mumegas
//      if      (255.5<c && c<767.5) x = (c-fNwires/2.-0.5) * fIp;
//      else if (767.5<=c)           x =  256 * fIp + (c-fNwires/2.-0.5-256) * fOp;
//      else                         x = -256 * fIp + (c-fNwires/2.+0.5+256) * fOp;
//    }
//    else { // old mumegas
//      x = (c-575.5) * fIp;
//    }

  int lastch=0;
  double lastp=0.;
  double dist=0.;
  double firstpitch=0.;
  for(std::map<int,double>::iterator ipitch=fPitchs.begin(); ipitch!=fPitchs.end(); ipitch++ ) {
    if(!firstpitch) firstpitch=ipitch->second;
    if(wire < ipitch->first) break;
    dist += (ipitch->first-lastch)*lastp;
    lastp = ipitch->second;
    lastch = ipitch->first;
  }
  dist += ((wire-lastch+0.5)*lastp - fHalfsize);

  return dist;
}

