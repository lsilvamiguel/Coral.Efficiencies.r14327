#include "CsMW1Pid.h"
#include "Reco/CalorimeterParticle.h"

using namespace std;
using namespace CLHEP;

CsMW1Pid::CsMW1Pid()
{
  // get detectors
  Coral* coral        = Coral::Instance();
  //CsOpt* opt = CsOpt::Instance();
  list<CsDetector*> dets = coral->getDetectors();
  list<CsDetector*>::iterator ItDet;

  int mw1i = 0;

  for(ItDet = dets.begin(); ItDet != dets.end(); ++ItDet)
  {
    if(!strncmp((*ItDet)->GetTBName().c_str(),"MA",2)) {
      mw1s.insert(make_pair(*ItDet,mw1i++));
      cout<<"Inserted "<<(*ItDet)->GetTBName()<<" with id "<<mw1i-1<<endl;
    }
  }

  for(int i = 0; i < 16; i++) {
    CsDetector* det = NULL;
    for(map<CsDetector*,int>::iterator di=mw1s.begin();di!=mw1s.end();++di) {
      if((*di).second == i) {
	det = (*di).first;
	break;
      }
    }
    if(!det) continue;
    string tit = string("h")+det->GetTBName()+"_Resi";
    string title = string(det->GetTBName()) + " residuals";
    hMW1Resi[i] = new TH1F(tit.c_str(), title.c_str(), 200, -50, 50);
  }
}


bool CsMW1Pid::doPid(const CsTrack* trk, int& nplanes, float& E)
{
  //cout << "CsMW1Pid::doPid() called." << endl;
  nplanes = -1;

  E = -1.;

  // Track selection criteria
  if(trk->getHelices().front().getCop() <= 0.) {
    return false;
  }

  // Exclude beam tracks
  if(trk->getHelices().front().getZ() < 0) return false;

  bool in_active_area = true;
  for(map<CsDetector*,int>::iterator di=mw1s.begin();di!=mw1s.end();++di) {
    THlx H;
    H.ImportHelix(trk->getHelices().front());
    THlx H1;
    
    H1(0) = (*di).first->getZcm()/10.;
    H.Extrapolate(H1);

    if(!(*di).first->inActiveArea(H1(1)*10.,H1(2)*10.)) {
      in_active_area = false;
      break;
    }
  }

  if(!in_active_area) return false;

  bool fired[2][8];
  for(unsigned int i = 0; i < 2; i++) {
    for(unsigned int j = 0; j < 8; j++) {
      fired[i][j] = false;
    }
  }
  for(map<CsDetector*,int>::iterator di=mw1s.begin();di!=mw1s.end();++di) {
    list<CsCluster*> clusters = (*di).first->getMyClusters();
    CsCluster* closer_cls = NULL;
    double minres = 1000000;
    double Utrack, sigmaU;
    THlx H;
    H.ImportHelix(trk->getHelices().front());
    THlx H1;
    
//     cout<<"CsMW1Pid::doPid(): plane "<<(*di).first->GetTBName()
// 	<<".getZcm() = "<<(*di).first->getZcm()<<endl;
    
    H1(0) = (*di).first->getZcm()/10.;
    H.Extrapolate(H1);
    
    HepMatrix wM = (*di).first->getRotWRS();
    
    Utrack = (wM(1,1)*H1(1) + wM(2,1)*H1(2) + wM(3,1)*H1(0))*10.;
    sigmaU = sqrt(wM(1,1)*H1(1,1) + wM(2,1)*H1(2,2));
//     cout<<"CsMW1Pid::doPid(): plane "<<(*di).first->GetTBName()
// 	<<": wM(1,1)*H1(1,1)="<<wM(1,1)<<"*"<<H1(1,1)
// 	<<"; wM(2,1)*H1(2,2)="<<wM(2,1)<<"*"<<H1(2,2)
// 	<<"; sigmaU=sqrt("<<wM(1,1)*H1(1,1) + wM(2,1)*H1(2,2)<<")"<<endl;
    for(list<CsCluster*>::iterator ci = clusters.begin();
	ci != clusters.end(); ++ci) {
      double  cU = (*ci)->getU();
      //double  cV = (*ci)->getV();
//       double  cW = (*ci)->getW();
//       cout<<"CsMW1Pid::doPid(): plane "<<(*di).first->GetTBName()
// 	  <<"; cW = "<<cW<<endl;
      
      double resi = ( cU - Utrack );
      if(fabs(resi) < fabs(minres)) {
	minres = resi;
	closer_cls = *ci;
      }
    }
    if((*di).second < 16) {
      hMW1Resi[(*di).second]->Fill(minres);
      float max_dU = sigmaU*3.>10. ? sigmaU*3. : 10.;
      if(fabs(minres) < max_dU) {
	int i = (*di).second / 8;
	int j = (*di).second % 8;
	fired[i][j] = true;
      }
    } else {
      cout<<(*di).first->GetTBName()<<" id = "<<(*di).second<<endl;
    }
  }

  int count[2] = {0,0};
  for(unsigned int i = 0; i < 2; i++) {
    for(unsigned int j = 0; j < 8; j++) {
      if(fired[i][j]) count[i]++;
    }
  }

  if(count[0] < 4) return false;

  //if(count[0] < count[1]) nplanes = count[0];
  //else 
  nplanes = count[1];

  CsEvent* event = CsEvent::Instance();
  vector< CsParticle * > pv = event->getParticles();
  vector< CsParticle * >::iterator pi;
  CsParticle * part = NULL;
  for(pi = pv.begin(); pi != pv.end(); ++pi) {
    if((*pi)->getTrack() == trk) {
      part = (*pi);
      break;
    }
  }

  if(!part) return true;

  vector< Reco::CalorimeterParticle * > cv = part->getCalObjects();
  vector< Reco::CalorimeterParticle * >::iterator ci;
  for(ci = cv.begin(); ci != cv.end(); ++ci) {
    //cout<<"CaloHistos::Fill(): E("<<(*ci)->GetCalorimeterName()
    //    <<") = "<<(*ci)->GetE()<<endl;
    if((*ci)->GetCalorimeterName()[0] != 'H') continue;
    switch((*ci)->GetCalorimeterName()[3]) {
    case '1':
      E = (*ci)->GetE();
      //cout<<"hHCAL1_E_vs_p["<<i<<"]->Fill("<<p<<","<<(*ci)->GetE()<<")"
      //    <<endl;
      break;
    case '2':
      //hHCAL2_E_vs_p[i]->Fill(p,(*ci)->GetE());
      break;
    }
  }

  return true;
}

