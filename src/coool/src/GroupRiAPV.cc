#include "GroupRiAPV.h"
#include "PlaneRiAPV.h"

#define RiAPV_A2_CUT 12

ClassImp(GroupRiAPV);


GroupRiAPV::GroupRiAPV(const char* name): Group(name), fNbPlanes(0), fHch2Dall(0), fHchamp2Dall(0),
    fHavgamp2Dall(0), fHsigamp2Dall(0), fHavgampCM2Dall(0), fHsigampCM2Dall(0),
    fHevtamp2Dall(0), fHevtampcut2Dall(0), fHapvvsa2meanall(0)
    {}


void GroupRiAPV::Init() {

//   addPdRef(0, PdRefType(0, 3, -1));
//   addPdRef(1, PdRefType(0, 2, -1));
//   addPdRef(2, PdRefType(1, 3, -1));
//   addPdRef(4, PdRefType(2, 3, -1));
//   addPdRef(6, PdRefType(3, 3, -1));
//   addPdRef(7, PdRefType(3, 2, -1));  // bad PD numbering here...
//   addPdRef(8, PdRefType(0, 1, -1));
//   addPdRef(9, PdRefType(0, 0, -1));
//   addPdRef(11, PdRefType(1, 0, -1));
//   addPdRef(13, PdRefType(2, 0, -1));
//   addPdRef(14, PdRefType(3, 1, -1));
//   addPdRef(15, PdRefType(3, 0, -1));

  addPdRef(0, PdRefType(3, 3, -1));
  addPdRef(1, PdRefType(3, 2, -1));
  addPdRef(2, PdRefType(2, 3, -1));
  addPdRef(4, PdRefType(1, 3, -1));
  addPdRef(6, PdRefType(0, 3, -1));
  addPdRef(7, PdRefType(0, 2, -1));
  addPdRef(8, PdRefType(3, 1, -1));
  addPdRef(9, PdRefType(3, 0, -1));
  addPdRef(11, PdRefType(2, 0, -1));
  addPdRef(13, PdRefType(1, 0, -1));
  addPdRef(14, PdRefType(0, 1, -1));
  addPdRef(15, PdRefType(0, 0, -1));

  char idpd[3];

  fNbPlanes = fPlanes.size();
  fNbPlanesUsed = 0;

  for (register int ii = 0; ii < fNbPlanes; ii++) {
    const char* pname = fPlanes[ii]->GetName();
    const PlaneRiAPV* riplane = dynamic_cast<const PlaneRiAPV*>(fPlanes[ii]);
    if (riplane == 0) {
      std::cerr<<"GroupRiAPV::Init: unknown Plane "<<pname<<" in group RiAPV"<<std::endl;
      return;
    }
    if (sscanf(pname, "RA01P%2c", idpd) == 0) {
      std::cerr<<"GroupRiAPV::Init: can't recognize PlaneRiAPV name "<<pname<<", ignoring it"<<std::endl;
      return;
    }
    idpd[2]=0;
    unsigned int pdnb = strtoul(idpd, 0, 10);
//     if (pdnb<0 || pdnb>15 || pdnb==3 || pdnb==5 || pdnb==10 || pdnb==12)  // unsigned ints are never < 0
    if (pdnb>15 || pdnb==3 || pdnb==5 || pdnb==10 || pdnb==12) {
      std::cerr<<"GroupRiAPV::Init: not valid PlaneRiAPV name "<<pname<<" (bad pd number), ignoring it"<<std::endl;
      return;
    }

    pdRef[pdnb].plindex = ii;
    fNbPlanesUsed++;
  }

  std::string name;

  // 2D Hit profile
  name = fName + "_2D_profile_all";
  fHch2Dall=new TH2F(name.c_str(),name.c_str(), 72*4, 0., 72.*4, 72*4, 0., 72.*4);
  fHch2Dall->SetOption("colz");
  fHch2Dall->GetXaxis()->SetNdivisions(12,kFALSE);
  fHch2Dall->GetYaxis()->SetNdivisions(4,kFALSE);
  AddHistogram(fHch2Dall);

  // 2D Hit profile
  name = fName + "_2D_amp_all";
  fHchamp2Dall=new TH2F(name.c_str(),name.c_str(), 72*4, 0., 72.*4, 72*4, 0., 72.*4);
  fHchamp2Dall->SetOption("colz");
  fHchamp2Dall->GetXaxis()->SetNdivisions(12,kFALSE);
  fHchamp2Dall->GetYaxis()->SetNdivisions(4,kFALSE);
  AddHistogram(fHchamp2Dall);

  if (fExpertHistos) {
    // 2D a2 average profile
    name = fName + "_2D_a2_avg_all";
    fHavgamp2Dall=new TH2D(name.c_str(),name.c_str(), 72*4, 0., 72.*4, 72*4, 0., 72.*4);
    fHavgamp2Dall->SetOption("col");
    AddHistogram(fHavgamp2Dall);

    // 2D a2 sigma profile
    name = fName + "_2D_a2_sig_all";
    fHsigamp2Dall=new TH2D(name.c_str(),name.c_str(), 72*4, 0., 72.*4, 72*4, 0., 72.*4);
    fHsigamp2Dall->SetOption("col");
    AddHistogram(fHsigamp2Dall);

    // 2D a2 CM average profile
    name = fName + "_2D_a2CM_avg_all";
    fHavgampCM2Dall=new TH2D(name.c_str(),name.c_str(), 72*4, 0., 72.*4, 72*4, 0., 72.*4);
    fHavgampCM2Dall->SetOption("col");
    AddHistogram(fHavgampCM2Dall);

    // 2D a2 CM sigma profile
    name = fName + "_2D_a2CM_sig_all";
    fHsigampCM2Dall=new TH2D(name.c_str(),name.c_str(), 72*4, 0., 72.*4, 72*4, 0., 72.*4);
    fHsigampCM2Dall->SetOption("col");
    AddHistogram(fHsigampCM2Dall);

    // 2D plot where bin value is the a2 mean for each apv (shown in x,y)
    name = fName + "_2D_a2_mean_vs_apv_all";
    fHapvvsa2meanall=new TH2F(name.c_str(),name.c_str(), 12*4, 0., 12.*4, 4*4, 0., 4.*4);
    fHapvvsa2meanall->SetOption("colz");
    AddHistogram(fHapvvsa2meanall);
  }

  // 2D event with a2 amplitude
  name = fName + "_2D_a2CM_event_all";
  fHevtamp2Dall=new TH2D(name.c_str(),name.c_str(), 72*4, 0., 72.*4, 72*4, 0., 72.*4);
  fHevtamp2Dall->SetOption("col");
  AddHistogram(fHevtamp2Dall);

  // 2D event with a2 amplitude with cut
  name = fName + "_2D_a2CM_cut_event_all";
  fHevtampcut2Dall=new TH2D(name.c_str(),name.c_str(), 72*4, 0., 72.*4, 72*4, 0., 72.*4);
  fHevtampcut2Dall->SetOption("col");
  AddHistogram(fHevtampcut2Dall);
}


void GroupRiAPV::addPdRef(unsigned int i, const PdRefType& rf) {
  pdRef.insert(make_pair(i,rf));
}


void GroupRiAPV::EndEvent(const CS::DaqEvent &event) {

  if (fNbPlanesUsed<1) return;
  if (thr_flag) TThread::Lock();

  typedef std::map<unsigned int, PdRefType>::iterator MI;
  for (MI pditer = pdRef.begin(); pditer != pdRef.end(); pditer++) {
    if (pditer->second.plindex == -1) continue;
    PdRefType& pd = pditer->second;
    const PlaneRiAPV* plane = dynamic_cast<const PlaneRiAPV*>(fPlanes[pd.plindex]);
    register int xpos = pd.xpos, ypos = pd.ypos;


    // Get digits (single strip data)
    for (register int idg = 0; idg < plane->fVch->GetNvalues(); idg++) {

//      register int channel = (int) plane->fVch->GetValues()[idg];
//       register int a0 = (int) plane->fVa0->GetValues()[idg];
//       register int a1 = (int) plane->fVa1->GetValues()[idg];
      register int a2 = (int) plane->fVa2->GetValues()[idg];
//       register int t = (int) plane->fVt->GetValues()[idg];
//       float a12 = (a2==0) ? 0 : float(a1)/float(a2);
//       float a02 = (a2==0) ? 0 : float(a0)/float(a2);
      register int xpx = (int) plane->fVxpx->GetValues()[idg];
      register int ypx = (int) plane->fVypx->GetValues()[idg];

      if (a2 > RiAPV_A2_CUT) fHch2Dall->Fill(xpx + 72*xpos, ypx + 72*ypos);
      fHchamp2Dall->Fill(xpx + 72*xpos, ypx + 72*ypos, a2);
    }

    if((fRateCounter&0x3f)==0) {

      register TH2D* avgamp2Dpl = plane->fHavgamp2D;
      if (fHavgamp2Dall && avgamp2Dpl) for (register int ix = 1; ix <= avgamp2Dpl->GetNbinsX(); ix++) {
        for (register int iy = 1; iy <= avgamp2Dpl->GetNbinsY(); iy++) {
          register double val = avgamp2Dpl->GetBinContent(ix, iy);
          fHavgamp2Dall->SetBinContent(ix+72*xpos, iy+72*ypos, val);
        }
      }

      register TH2D* sigamp2Dpl = plane->fHsigamp2D;
      if (fHsigamp2Dall && sigamp2Dpl) for (register int ix = 1; ix <= sigamp2Dpl->GetNbinsX(); ix++) {
        for (register int iy = 1; iy <= sigamp2Dpl->GetNbinsY(); iy++) {
          register double val = sigamp2Dpl->GetBinContent(ix, iy);
          fHsigamp2Dall->SetBinContent(ix+72*xpos, iy+72*ypos, val);
        }
      }

      register TH2D* avgampCM2Dpl = plane->fHavgampCM2D;
      if (fHavgampCM2Dall && avgampCM2Dpl) for (register int ix = 1; ix <= avgampCM2Dpl->GetNbinsX(); ix++) {
        for (register int iy = 1; iy <= avgampCM2Dpl->GetNbinsY(); iy++) {
          register double val = avgampCM2Dpl->GetBinContent(ix, iy);
          fHavgampCM2Dall->SetBinContent(ix+72*xpos, iy+72*ypos, val);
        }
      }

      register TH2D* sigampCM2Dpl = plane->fHsigampCM2D;
      if (fHsigampCM2Dall && sigampCM2Dpl) for (register int ix = 1; ix <= sigampCM2Dpl->GetNbinsX(); ix++) {
        for (register int iy = 1; iy <= sigampCM2Dpl->GetNbinsY(); iy++) {
          register double val = sigampCM2Dpl->GetBinContent(ix, iy);
          fHsigampCM2Dall->SetBinContent(ix+72*xpos, iy+72*ypos, val);
        }
      }

      register TH2D* evtamp2Dpl = plane->fHevtamp2D;
      if (fHevtamp2Dall && evtamp2Dpl) for (register int ix = 1; ix <= evtamp2Dpl->GetNbinsX(); ix++) {
        for (register int iy = 1; iy <= evtamp2Dpl->GetNbinsY(); iy++) {
          register double val = evtamp2Dpl->GetBinContent(ix, iy);
          fHevtamp2Dall->SetBinContent(ix+72*xpos, iy+72*ypos, iy, val);
        }
      }

      register TH2D* evtampcut2Dpl = plane->fHevtampcut2D;
      if (fHevtampcut2Dall && evtampcut2Dpl) for (register int ix = 1; ix <= evtampcut2Dpl->GetNbinsX(); ix++) {
        for (register int iy = 1; iy <= evtampcut2Dpl->GetNbinsY(); iy++) {
          register double val = evtampcut2Dpl->GetBinContent(ix, iy);
          fHevtampcut2Dall->SetBinContent(ix+72*xpos, iy+72*ypos, val);
        }
      }

      register TH2F* apvvsa2mean = plane->fHapvvsa2mean;
      if (fHapvvsa2meanall && apvvsa2mean) for (register int ix = 1; ix <= apvvsa2mean->GetNbinsX(); ix++) {
        for (register int iy = 1; iy <= apvvsa2mean->GetNbinsY(); iy++) {
          register double val = apvvsa2mean->GetBinContent(ix, iy);
          fHapvvsa2meanall->SetBinContent(ix+12*xpos, iy+4*ypos, val);
        }
      }

    }
  }

  if (thr_flag) TThread::UnLock();
}




