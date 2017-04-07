
#include "GroupECAL1.h"

#define RiAPV_A2_CUT 12

ClassImp(GroupECAL1);


GroupECAL1::GroupECAL1(const char* name): Group(name), fPlanesOK(false), fPlaneECAL1(0),
    fPlaneFEM(0), fHch2Dall(0), fHchamp2Dall(0),
    fHavgamp2Dall(0), fHsigamp2Dall(0), fHavgampCM2Dall(0), fHsigampCM2Dall(0),
    fHevtamp2Dall(0), fHevtampcut2Dall(0), fHapvvsa2meanall(0), hLasAmplrenF(1056)
    {}


void GroupECAL1::Init() {

  bool noupdate=true;

  int nbPlanes = fPlanes.size();
  if (nbPlanes == 0) return;
  
  RefLedMin=0.01;
  RefLedMax=3.;

  for (register int ii = 0; ii < nbPlanes ; ii++) {
    const PlaneFEM* planeFEM = dynamic_cast<const PlaneFEM*>(fPlanes[ii]);
    const PlaneECAL1* planeECAL1 = dynamic_cast<const PlaneECAL1*>(fPlanes[ii]);
    if (planeFEM ) {
      fPlaneFEM = planeFEM;
    }
    if (planeECAL1 ) {
      fPlaneECAL1 = planeECAL1;
    }
  }

  if ( fPlaneFEM == 0 || fPlaneECAL1 == 0 ) return;

  fPlanesOK = true;
  int nchan = fPlaneECAL1->GetNchannels();
  int nrow = fPlaneECAL1->GetNrows();
  int ncol = fPlaneECAL1->GetNcols();

  string name = fName + "_Las/Reference";
  string title = fName + "Las ren FEM / Beam offspill reference(%)";
  href= new TProfile(name.c_str(), title.c_str(), nchan,0.,nchan,"s");
  AddHistogram(href);

  name = fName + "_Las/RMS_renFEM_vs_chnl#";
  title = fName + "Las/RMS ren FEM vs channel#";
  hLedProf=new TProfile(name.c_str(),title.c_str(),nchan, 0., nchan,"s");
  AddHistogram(hLedProf);

  name = fName + "_MnLas_renFEM_vs_chnl#";
  title = fName + " Mean Las ren FEM vs channel#";
  hmnled=new TH1F_Ref(name.c_str(),title.c_str(), 
                              nchan, 0., nchan, fRateCounter, noupdate);
  AddHistogram(hmnled);
  
  if (fExpertHistos) {
    for ( int icol=0; icol<ncol; icol++) {
    	for( int irow=0; irow<nrow; irow++ ) {
	int ch=irow + icol*nrow;
 	char b1[222],b2[222];
	
	sprintf(b1,"%s_%d_%d_Las_Ampl_renFEM",fName.c_str(),icol,irow);
	sprintf(b2,"%s %d %d Las Amplitude ren FEM",fName.c_str(),icol,irow);
	hLasAmplrenF[ch]=new TH1F_Ref(b1,b2,1000, 0., 2., fRateCounter);
	AddHistogram(hLasAmplrenF[ch]);
	
	}
    }
  }
  
  OpenReference();
  if (fReferenceDirectory) {
  ((TH1F_Ref*)hmnled)->SetReference(fReferenceDirectory);
  	if (fExpertHistos) {
	  for ( int icol=0; icol<ncol; icol++) {
	    for( int irow=0; irow<nrow; irow++ ) {
	      int ch=irow + icol*nrow;
  		((TH1F_Ref*)hLasAmplrenF[ch])->SetReference(fReferenceDirectory);
	    }
	  }
  	}
  }  
  
}


void GroupECAL1::EndEvent(const CS::DaqEvent &event) {

  if (fPlanesOK == false) return;

  if ( fPlaneFEM == 0 || fPlaneECAL1 == 0 ) return;

  int nchan = fPlaneECAL1->GetNchannels();
  const double * pulseECAL1 = fPlaneECAL1->GetPulseAmp();
  const vector<double>& pulseFEM = fPlaneFEM->GetFEMAmp();

  if (thr_flag) TThread::Lock();

  if (pulseFEM[2]) { 
  if (hmnled!=NULL) hmnled->Reset();
  const float* myref=((TH1F_Ref*)hmnled)->GetRef();
    for (register int i=0; i<nchan; i++) {
    hLedProf->Fill(i, pulseECAL1[i] / pulseFEM[2]);
    hmnled->Fill(i, hLedProf->GetBinContent(i));
    if(myref!=NULL) {
    double Ref=myref[i];
    double Diff=100.;
    if(Ref>RefLedMin&&Ref<RefLedMax) Diff=(hLedProf->GetBinContent(i)/Ref)*100.;
       if(Diff<2.) Diff=2.;  
    href->Fill(i, Diff);
    }
  if (fExpertHistos) {
    hLasAmplrenF[i]->Fill(pulseECAL1[i] / pulseFEM[2]);
    }
  }
  }

  if (thr_flag) TThread::UnLock();
}




