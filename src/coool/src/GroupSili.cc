#include <cmath>
#include "GroupSili.h"
#include "TThread.h"
#include "TStyle.h"
#include "GroupPanel.h"
#include "PlaneSili.h"
#include "TriggerTime.h"
#include "DaqEvent.h"
#include "GeomPlaneSili.h"

ClassImp(GroupSili);
ClassImp(GroupSiliTrack);

void GroupSili::Init() {
  // wish from M. Wiesmann for this histogram: swap sign of V vertical axis
  // in order to have equivalent, by 5 deg inclined beam profile for UV
  // as for XY
  std::string name0 = fName + "_-V_vs_U";
  fHistList.push_back(new TH2F(name0.c_str(), name0.c_str(),
			       320, -4.5, 4.5, 320, -4.5, 4.5));
  fHistList[0]->SetStats(false);
  fHistList[0]->GetXaxis()->SetTitle("u [cm]");
  fHistList[0]->GetYaxis()->SetTitle("-v [cm]");

  std::string name1 = fName + "_Y_vs_X";
  fHistList.push_back(new TH2F(name1.c_str(), name1.c_str(),
			       320, -4.5, 4.5, 320, -4.5, 4.5));

  fHistList[1]->SetStats(false);

  std::string name2 = fName + "_V_U_AmplCorr";
  fHistList.push_back(new TH2F(name2.c_str(), name2.c_str(),
			       250, -0.5, 499.5, 250, -0.5, 499.5));

  std::string name3 = fName + "_X_Y_AmplCorr";
  fHistList.push_back(new TH2F(name3.c_str(), name3.c_str(),
			       250, -0.5, 499.5, 250, -0.5, 499.5));


  
}

void GroupSili::EndEvent(const CS::DaqEvent &event) {

  if (thr_flag) TThread::Lock();


  for (int k=0; k<2; k++) {
    // k: 0-UV 1-XY

    const GeomPlane *xgeom = NULL;
    const GeomPlane *ygeom = NULL;

    for (register unsigned int pl=0; pl<fPlanes.size(); pl++) {
      char planeLetter = fPlanes[pl]->GetName()[4];
      if (planeLetter =='U' && k==0) xgeom = fPlanes[pl]->GetGeometry();
      if (planeLetter =='V' && k==0) ygeom = fPlanes[pl]->GetGeometry();
      if (planeLetter =='X' && k==1) xgeom = fPlanes[pl]->GetGeometry();
      if (planeLetter =='Y' && k==1) ygeom = fPlanes[pl]->GetGeometry();
    }

    if (xgeom && ygeom) {
      const std::set<CCluster1>& xclu = xgeom->GetClusters();
      const std::set<CCluster1>& yclu = ygeom->GetClusters();

      int xsize = xclu.size();
      int ysize = yclu.size();
      float xamp[xsize], xamp1[xsize];
      float yamp[ysize], yamp1[ysize];
      float xpos[xsize];
      float ypos[ysize];

      int ix=0;
      for (std::set<CCluster1>::iterator xcl=xclu.begin();
	   xcl!=xclu.end(); xcl++) if (xcl->dt.size()>=3) {
	     xamp [ix] = xcl->dt[0] + xcl->dt[1] + xcl->dt[2];
	     xamp1[ix] = xcl->dt[1] / xcl->dt[2];
	     xpos[ix] = xcl->pos;
	     ix++;
	   }

      int iy=0;
      for (std::set<CCluster1>::iterator ycl=yclu.begin();
	   ycl!=yclu.end(); ycl++) if (ycl->dt.size()>=3) {
	     yamp [iy] = ycl->dt[0] + ycl->dt[1] + ycl->dt[2];
	     yamp1[iy] = ycl->dt[1] / ycl->dt[2];
	     ypos[iy] = ycl->pos;
	     iy++;
	   }
      
      int msize = xsize<ysize?xsize:ysize;
      int xindex[msize], yindex[msize];
      float distances[msize];
      for (int i=0; i<msize; i++) distances[i]=10000;

      for (int iix=0; iix<ix; iix++)
	for (int iiy=0; iiy<iy; iiy++) {
	  //float this_dist = fabs(float(xamp[iix]-yamp[iiy])/
	  //		 float(xamp[iix]+yamp[iiy]));
	  float this_dist = fabs(float(xamp1[iix]-yamp1[iiy])/
				 float(xamp1[iix]+yamp1[iiy])) + 
	    0.1 * fabs(float(xamp[iix]-yamp[iiy])/
		       float(xamp[iix]+yamp[iiy]));;
	  for (int i=msize-1; i>=0; i--) {
	    if (this_dist<0.5 && this_dist<distances[i]) {
	      if (i!=msize-1) {
		distances[i+1]=distances[i];
		xindex   [i+1]=xindex   [i];
		yindex   [i+1]=yindex   [i];
	      }
	      distances[i]=this_dist;
	      xindex   [i]=iix;
	      yindex   [i]=iiy;
	    }
	  }
	}
      for (int i=0; i<msize; i++) if (distances[i]!=10000) {
	if (k==0) fHistList[k]->Fill(xpos[xindex[i]], -ypos[yindex[i]]);
	else      fHistList[k]->Fill(xpos[xindex[i]],  ypos[yindex[i]]);
	fHistList[k+2]->Fill(xamp[xindex[i]], yamp[yindex[i]]);
      }
    }
  }

  if (thr_flag) TThread::UnLock();
}

void GroupSili::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}


void GroupSiliTrack::Init() {

#if USE_TRACK == 1


  maxNumberofClusters=6;   // Clusters per planes allowed
  scifiTimeDiffCut=2;
  scifiTriggerCut=2;
  siliconTimeCutMode=0;
  siliconActiveTimeCut=0;
  chi2Cut=5;
  int minNumberofPlanes=6;
  double chi2probCut=1e-5;
  deactivatePlane=0;
  char tps[10][500];

  char parFileName[500];
  sprintf(parFileName, "/afs/e18.ph.tum.de/compass/analysis/primakoff/timecalib/analysis_coool/configs/TrackPar/Tracking_%s.par", this->GetName());  
  ifstream trackpar(parFileName);
  if (!trackpar){ trackpar.open("Tracking.par");
  std::cout << "NO Track par File in " << parFileName << " found! USE DEFAULT VALUES!!! " << std::endl; 
  }
  if (trackpar) {
    trackpar >> maxNumberofClusters >> tps[0];// >> std::endl;
    trackpar >> minNumberofPlanes   >> tps[1];// >> std::endl;
    trackpar >> chi2probCut         >> tps[2];// >> std::endl;
    trackpar >> chi2Cut             >> tps[3];// >> std::endl;
    trackpar >> scifiTimeDiffCut    >> tps[4];// >> std::endl;
    trackpar >> scifiTriggerCut     >> tps[5];// >> std::endl;
    trackpar >> siliconTimeCutMode  >> tps[6];// >> std::endl;
    trackpar >>  siliconActiveTimeCut >>tps[7];
    trackpar >> deactivatePlane     >> tps[8];// >> std::endl;
  }
  trackpar.close();

  std::cout << "Read from " << parFileName <<std::endl;
  std::cout << maxNumberofClusters << " is maxNumberofClusters" << std::endl;
  std::cout << minNumberofPlanes   << " is minNumberofPlanes"   << std::endl;
  std::cout << chi2probCut         << " is chi2probCut"         << std::endl;
  std::cout << chi2Cut             << " is chi2Cut"             << std::endl;
  std::cout << scifiTimeDiffCut    << " is scifiTimeDiffCut"    << std::endl;
  std::cout << scifiTriggerCut     << " is scifiTriggerCut"     << std::endl;
  std::cout << siliconTimeCutMode  << " is siliconTimeCutMode"  << std::endl;
  std::cout << siliconActiveTimeCut << " is  siliconActiveTimeCut" << std::endl;
 
 if (deactivatePlane==0)
    std::cout << " All planes active" << std::endl;
  else {
    deactivatedPlaneName = new char[strlen(tps[8])+1];
    sprintf(deactivatedPlaneName, "%s", tps[8]);
    std::cout << "Plane " <<   deactivatedPlaneName << " will be deactivated." 
	 << std::endl;
    deactivatedPlane = new Tracker*[deactivatePlane];
    for (int i=0; i<deactivatePlane; i++) deactivatedPlane[i]=NULL;
  }
  // Tracking for the Silis
  fTrackMaker = new TrackMaker(GetName(), 2, chi2probCut);
  fTrackMaker->SetChi2(chi2Cut);
  // 2: not space points (chi2)
  
  // std::cout << " #### fPlanes.size()  "<<  fPlanes.size()<< std::endl;
  for(unsigned i=0; i<fPlanes.size(); i++) {
    if ( strncmp(fPlanes[i]->GetName(), "SI", 2) == 0)
      silName.push_back(fPlanes[i]->GetName());

    if ( strncmp(fPlanes[i]->GetName(), "FI", 2) == 0){
      sciName.push_back(fPlanes[i]->GetName());
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
			      inpitch,inpitch, ctmin,ctmax,0);
    }
  }
  if (deactivatePlane) {
    int dap=0;
    for(unsigned i=0; i<fPlanes.size(); i++) {  
      for (unsigned s=0; s<strlen(deactivatedPlaneName)-4; s++) {
	if (strncmp(fPlanes[i]->GetName(), deactivatedPlaneName+s, 5)==0) {
	  if(const GeomPlane *geom = fPlanes[i]->GetGeometry()) {
	    if (dap==deactivatePlane) {
	      std::cout << "Too many deactivate Plane names..." << std::endl;
	      break;
	    }
	    deactivatedPlane[dap]=fTrackMaker->GetTracker(geom->GetID());
	    if (deactivatedPlane[dap]) {
	      std::cout << "Plane "<<deactivatedPlane[dap]->GetName()
		   << " deactivated."<<std::endl;
	      deactivatedPlane[dap]->Activate(0);
	      deactivatedPlane[dap]->SetRMin(0);
	      // include error on SciFi track (>> than Silicon resolution)
	      if (strncmp(deactivatedPlaneName, "SI", 2) == 0)
		deactivatedPlane[dap]->SetCheckRange(-0.06, 0.06);
	    }
	    dap++;
	  }
	}
      }
    }
  }
  
  fTrackMaker->Init();
  fTrackMaker->SetNTrackers(minNumberofPlanes);
  // Adds tracking histogram to this's histogram list
  std::vector<TH1*>& trackhists=fTrackMaker->GetHistograms();
  for(unsigned i=0; i<trackhists.size();i++){
    float halfWidth=3.5;
    if (strncmp(trackhists[i]->GetName()+9, "heff", 4) == 0 ||
	strncmp(trackhists[i]->GetName()+9, "tcks", 4) == 0 ||
	strncmp(trackhists[i]->GetName()+9, "bg",   2) == 0 ) {
      trackhists[i]->GetXaxis()->SetLimits(-halfWidth, halfWidth);
      trackhists[i]->GetYaxis()->SetLimits(-halfWidth, halfWidth);
      trackhists[i]->SetOption("colz");
      gStyle->SetPalette(1);
    }
    if (strncmp(trackhists[i]->GetName()+9, "h1deff", 6) == 0 ||
	strncmp(trackhists[i]->GetName()+9, "1dtcks", 6) == 0 ||
	strncmp(trackhists[i]->GetName()+9, "1dbg",   4) == 0 ) {
      trackhists[i]->GetXaxis()->SetLimits(-halfWidth, halfWidth);
      trackhists[i]->GetYaxis()->SetLimits(-halfWidth, halfWidth);
    }
    fHistList.push_back(trackhists[i]);
  }
  
  hSiliR0VsSciTCS = new TH2*[silName.size()];
  hSiliR1VsSciTCS = new TH2*[silName.size()];
  hSiliR0VsR1 = new TH2*[silName.size()];
  char hsiln[3][silName.size()][50];
  
  for (unsigned i=0; i< silName.size(); i++) {
    sprintf(hsiln[0][i],"%sSiliR0VsSci", silName[i]);
    fHistList.push_back(new TH2F(hsiln[0][i], hsiln[0][i], 
				 420, -60, 150, 128, -0.1, 3.1));
    hSiliR0VsSciTCS[i] = (TH2*) fHistList.back();
    hSiliR0VsSciTCS[i]->SetOption("col");
    
    sprintf(hsiln[1][i],"%sSiliR1VsSci", silName[i]);
    fHistList.push_back(new TH2F(hsiln[1][i], hsiln[1][i], 
				 420, -60, 150, 128, -0.1, 3.1));
    hSiliR1VsSciTCS[i] = (TH2*) fHistList.back();
    hSiliR1VsSciTCS[i]->SetOption("col");

    sprintf(hsiln[1][i],"%sSiliR0VsR1", silName[i]);
    fHistList.push_back(new TH2F(hsiln[1][i], hsiln[1][i], 
				 150, -0.1, 3.0, 150, -0.1, 3.0));
    hSiliR0VsR1[i] = (TH2*) fHistList.back();
    hSiliR0VsR1[i]->SetOption("col");
  }
  
  std::string name7 = fName + "Beam Profile";
  fHistList.push_back(new TH2F(name7.c_str(), name7.c_str(), 
			       640, -8, 8, 640, -8, 8));
  hHitMap = (TH2*) fHistList.back();
  hHitMap->SetOption("col");
  
  
  std::string name8 = fName + "Beam Divergence Profile";
  fHistList.push_back(new TH2F(name8.c_str(), name8.c_str(), 
			       500, -0.005, 0.005, 500, -0.005, 0.005));
  hHitDMap = (TH2*) fHistList.back();
  hHitDMap->SetOption("col");
  
  std::string name9 = fName + "_dT_scifi";
  fHistList.push_back(new TH1F(name9.c_str(), name9.c_str(), 100, -20, 20));
  hdTSci = fHistList.back();
  
  //std::string name8 = fName + "_Chi22";
  //fHistList.push_back(new TH1F(name8.c_str(), name8.c_str(), 
  //                             1000, 0, 100000));
  //hChi22 = fHistList.back();

#endif
}

void GroupSiliTrack::ResetHistograms() {
#if USE_TRACK == 1
  Group::ResetHistograms();
  fTrackMaker->ResetHistos();
#endif
}

void GroupSiliTrack::Tracking() {
#if USE_TRACK == 1
if (thr_flag) TThread::Lock();

  // reset tracking
  fTrackMaker->Reset();
  
  // passing clusters to tracking
  // time cuts applied
  for(unsigned i=0; i<fPlanes.size(); i++) {
    if(const GeomPlane *geom = fPlanes[i]->GetGeometry()) {
      int isNoFI = strncmp(fPlanes[i]->GetName(), "FI", 2);
      int isNoSI = strncmp(fPlanes[i]->GetName(), "SI", 2);
      if (isNoFI&&isNoSI) continue;
      bool SIactive = fTrackMaker->
	GetTracker(fPlanes[i]->GetGeometry()->GetID())->IsActive();
      const std::set<CCluster1>& clusters = geom->GetClusters();
      for(std::set<CCluster1>::iterator ic=clusters.begin();
	  ic!=clusters.end(); ic++) {
	
	/*	std::cout << "#### if Bedingung " << " isNoFI " << isNoFI
		  << "  fabs(ic->dt[0]) " <<   fabs(ic->dt[0]) 
		  << "  isNoSI " << isNoSI
		  << "  SIactive " << SIactive
		  << " fabs(ic->dt[7]) " << fabs(ic->dt[7])
		  << std::endl;*/

	if ( ( isNoFI || fabs(ic->dt[0])<scifiTriggerCut )
	     && 
	     ( isNoSI ||
	       // do tracking only with in-time, but take all inactive plane
	       // clusters:
		(SIactive==false
		 && (fabs(ic->dt[7])<float(siliconTimeCutMode)))
		||
		(SIactive==true
		 && (fabs(ic->dt[7])<float(siliconActiveTimeCut) ) ) ) 
	     //(ic->dt[3]<1.15 && fabs(ic->dt[4]-1.0)<0.5) )   )
	     ) {

	  fTrackMaker->AddClusterWRS(geom->GetID(), 
				     ic->pos, ic->size, ic->res);
	  
	  // Find out which Tracker is filled: Here
	  // we still know cluster timing and can transfer it
	  const std::vector<Tracker*>& trs = fTrackMaker->GetTrackers();
	  
	  
	  /*  std::cerr<<"Our trackers:"<<std::endl;
	  for (int ii=0; ii<trs.size(); ii++) {
	    std::cerr <<trs[ii]<<" ";
	  }
	  std::cerr<<std::endl;
	  */
	  
	  if (strncmp(fPlanes[i]->GetName(), "SI", 2)==0) {
	    tcsphase = ((PlaneSili*)fPlanes[i])->GetTCSphase();
	    for (unsigned int ii=0; ii<trs.size(); ii++){
	      if (unsigned (geom->GetID()) == trs[ii]->GetId() ) {
		const std::vector<Cluster*>& cli = trs[ii]->GetClusters();
		cli[cli.size()-1]->fTot = ic->dt[2]; // Amplitude a2
		cli[cli.size()-1]->fT   = ic->dt[3]; // a0/a2
		cli[cli.size()-1]->fT2  = ic->dt[4]; // a1/a2
		cli[cli.size()-1]->fT_APV_0 = ic->dt[5]; // t(a0/a2)
		cli[cli.size()-1]->fT_APV_1 = ic->dt[6]; // t(a1/a2)
		cli[cli.size()-1]->fT_APV_c = ic->dt[7]; // weighted mean
	      }
	    }
	  }
	  
	  if (strncmp(fPlanes[i]->GetName(), "FI", 2)==0){
	    for (unsigned int ii=0; ii<trs.size(); ii++)
	      if (unsigned (geom->GetID()) == trs[ii]->GetId() ) {
		const std::vector<Cluster*>& cli = trs[ii]->GetClusters();
		cli[cli.size()-1]->fT=ic->dt[0];
	      }
	  }
	}
      }
    }
  }
  
  fTrackMaker->Clusterize(maxNumberofClusters);
  fTrackMaker->Candidates();
  fTrackMaker->Residuals();
  
  // Candidates() 
  //     calls NextTracker() 
  //              calls MakeTrack() 
  //                       calls fTrackers[i]->Efficiency(fTracks.back(), 0);
  //                                             (no histogram filling)
  //     and
  //     calls if (fTracks.size()<2) 
  //                 TrackMaker::Efficiency();
  //                       calls fTrackers[i]->Efficiency(fTracks[j], 1); 
  //                                             (histogram filling)

  for (int dap=0; dap<deactivatePlane; dap++)
    if (deactivatedPlane[dap]) 
      deactivatedPlane[dap]->EffProfile();
 
   
  const std::vector<Track*>& goodTracks = fTrackMaker->GetAllTracks();
  
  //  std::cout << "Tracks:" << goodTracks.size() << std::endl;

  ntrack=goodTracks.size();
  if (ntrack>=1) {
    crossing_point[0]=goodTracks[0]->GetY();
    crossing_point[1]=goodTracks[0]->GetZ();
    crossing_point[2]=goodTracks[0]->GetX();
    slopes[0]=goodTracks[0]->GetYp();
    slopes[1]=goodTracks[0]->GetZp();
    if (ntrack>=2) {
      crossing_point2[0]=goodTracks[1]->GetY();
      crossing_point2[1]=goodTracks[1]->GetZ();
      crossing_point2[2]=goodTracks[1]->GetX();
      slopes2[0]=goodTracks[1]->GetYp();
      slopes2[1]=goodTracks[1]->GetZp();
    }
    else {
      crossing_point2[0]=-10000;
      crossing_point2[1]=-10000;
      crossing_point2[2]=-10000;
      slopes2[0]=-10000;
      slopes2[1]=-10000;
    }
  }
  else {
    crossing_point[0]=-10000;
    crossing_point[1]=-10000;
    crossing_point[2]=-10000;
    slopes[0]=-10000;
    slopes[1]=-10000;
    crossing_point2[0]=-10000;
    crossing_point2[1]=-10000;
    crossing_point2[2]=-10000;
    slopes2[0]=-10000;
    slopes2[1]=-10000;
  }
  
  unsigned NBSIL = silName.size(); // number of silicon planes
  
  if(firstread){ 
    silMult = new int[NBSIL];    // cluster multiplicity[12] changed to[20]
    silResi = new float[NBSIL];  // cluster residual in current track
    silPosi = new float[NBSIL];  // cluster position
    silAmpl = new float[NBSIL];  // cluster amplitude
    silRat0 = new float[NBSIL];  // cluster r0
    silRat1 = new float[NBSIL];  // cluster r1
    silT0   = new float[NBSIL];    // time from r0
    silT1   = new float[NBSIL];    // time from r1
    silTa   = new float[NBSIL];    //
       
    tsit = new TTree(this->GetName(),"Silicon Tracking"); 
    for(unsigned i=0; i< NBSIL; i++) {
      TString sn(silName[i]);
      tsit->Branch(sn+"mult", &silMult[i], sn+"/I");
      tsit->Branch(sn+"resi", &silResi[i], sn+"/F");
      tsit->Branch(sn+"posi", &silPosi[i], sn+"/F");
      tsit->Branch(sn+"ampl", &silAmpl[i], sn+"/F");
      tsit->Branch(sn+"rat0", &silRat0[i], sn+"/F");
      tsit->Branch(sn+"rat1", &silRat1[i], sn+"/F");
      tsit->Branch(sn+"t0",   &silT0[i]  , sn+"/F");
      tsit->Branch(sn+"t1",   &silT1[i]  , sn+"/F");
      tsit->Branch(sn+"ta",   &silTa[i]  , sn+"/F");
    }
    tsit->Branch("timefib",     &TimeFiber  , "timefib/F");
    tsit->Branch("timephase",   &TimePhase  , "timepha/F");
    tsit->Branch("chi2",  &chi2,  "chi2/F");
    tsit->Branch("chi2p", &chi2p, "chi2p/F");
    firstread=false;
  }
  
  //  std::cout << "####### goodTracks.size() "<< goodTracks.size() << std::endl;
  for (unsigned int i=0; i< goodTracks.size(); i++) {
    //    hChi22->Fill(goodTracks[i]->GetChi2());
    //    std::cout << "####### IN good track loop " << std::endl;
    //   std::cout << "goodTracks[i]->GetChi2() "<< goodTracks[i]->GetChi2() 
    //	       <<  " <  " << chi2Cut <<  " ? "<< std::endl;
    // std::cout << " goodTracks[i]->GetChi2Prob() " << goodTracks[i]->GetChi2Prob() << std::endl;

    if (goodTracks[i]->GetChi2() < chi2Cut) {
      chi2  = goodTracks[i]->GetChi2();
      chi2p = goodTracks[i]->GetChi2Prob();
      const std::vector<Cluster*> trClusters = goodTracks[i]->GetClusters();
      // std::vector<Cluster*> trClusters(goodTracks[i]->GetClusters());
      /*
      std::cerr << "There is good track with "<<trClusters.size()<<"clusters"<<std::endl;
      std::cerr << "addresses:";
      for (int ii=0; ii<trClusters.size(); ii++) {
	std::cerr<<trClusters[ii]<<" ";
	if (trClusters[ii]) std::cerr<<"("<<trClusters[ii]->GetTracker()<<")";
      }
      std::cerr << std::endl;
      */

      if (trClusters.size()!=0) {
	float fitime[sciName.size()];
	float minfitime= 100000;
	float maxfitime=-100000;
	for (unsigned int j=0; j<trClusters.size(); j++) {
	  Cluster *jcl = trClusters[j];
	  if (jcl) {
	    Tracker *jtr = jcl->GetTracker();
	    if (jtr) {
	      if (jtr->GetName()) {
		for (unsigned int k=0; k< sciName.size(); k++) {
		  if ( strcmp(jtr->GetName(), sciName[k]) == 0) {
		    fitime[k] = jcl->fT;
		    if (fitime[k]>maxfitime) maxfitime=fitime[k];
		    if (fitime[k]<minfitime) minfitime=fitime[k];
		  }
		}
	      }
	    }
	  }
	}
	hdTSci->Fill(maxfitime-minfitime);
	//	std::cout << " ##### maxfitime-minfitime " 
	//	     << fabs(maxfitime-minfitime) << " < " 
	//	     << scifiTimeDiffCut << " ? " << std::endl;
	if (fabs(maxfitime-minfitime)<scifiTimeDiffCut) {
	  float fitime = (maxfitime+minfitime)/2;
	  TimeFiber  =  fitime;
	  TimePhase  =  tcsphase;
	  double oldX = goodTracks[i]->GetX();
	  for (unsigned int j=0; j<trClusters.size(); j++) {
	    if (trClusters[j]) {
	      Cluster *jcl = trClusters[j];
	      Tracker *jtr = jcl->GetTracker();
	      for (unsigned int k=0; k<silName.size(); k++) {

		if ( strcmp(trClusters[j]->GetTracker()->GetName(),
			    silName[k]) == 0) {
		  //	  std::cout << "###fitime " <<  fitime << " tcsphase " << tcsphase << std::endl; 
		  float fiti_corr = -fitime + tcsphase;
		  float r1 = trClusters[j]->fT;
		  float r2 = trClusters[j]->fT2;
		  hSiliR0VsSciTCS[k]->Fill(fiti_corr, r1);
		  hSiliR1VsSciTCS[k]->Fill(fiti_corr, r2);
                  hSiliR0VsR1[k]->Fill(r2, r1);
		  silMult[k] = (int)trClusters[j]->fSize;
		  silPosi[k] = trClusters[j]->fPos;
		  silAmpl[k] = trClusters[j]->fTot;
		  silRat0[k] = r1;
		  silRat1[k] = r2;
		  silT0[k]   =  trClusters[j]->fT_APV_0;
		  silT1[k]   =  trClusters[j]->fT_APV_1;
		  silTa[k]   =  trClusters[j]->fT_APV_c;
		  goodTracks[i]->Move(jtr->GetX());
		  silResi[k] = jtr->
		    GetResidualOfCluster(trClusters[j]->fPos, 
					 goodTracks[i]->GetY(),
					 goodTracks[i]->GetZ());
		  //  std::cout <<silMult[k]<<","<<silPosi[k]<<","
		  //   <<silResi[k]<<std::endl;
		}
	      }
	    }
	  }
	  goodTracks[i]->Move(oldX);
	  hHitMap->Fill(goodTracks[i]->GetY(), goodTracks[i]->GetZ());
	  hHitDMap->Fill(goodTracks[i]->GetYp(), goodTracks[i]->GetZp());
	}
      }
      tsit->Fill();
    }
  }
  
  if (thr_flag) TThread::UnLock();
  
#endif
}

void GroupSiliTrack::EndEvent(const CS::DaqEvent &event) {
#if USE_TRACK == 1
  //std::cout << "will process event" << event.GetEventNumberInRun() << " in Run "
  //     <<  event.GetRunNumber() << std::endl;
#endif
}

void GroupSiliTrack::ControlPanel(const TGWindow *p, const TGWindow *main) {
  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}











