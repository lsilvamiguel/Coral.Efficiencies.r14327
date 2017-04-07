#include "GroupRICH_MAPMT.h"
#include "PlaneRICH_MAPMT.h"
#include "TStyle.h"

/* ----------------------------------------------------------------------
   Maintainer & author: Frank Nerling 
   New Class for monitoring  RICH_MAPMT: (x, y, time), created April 2006
   Here: Group of all 4 MAPMT Quadrants  u
   Updated Version: 11.09.2006 (time cut fitted to 2007 data May 2007)
   ---------------------------------------------------------------------- */

ClassImp(GroupRICH_MAPMT)

void GroupRICH_MAPMT::Init() {

  gStyle -> SetPalette(1);

  fNbPlanes = fPlanes.size();
  std::string name;
  std::string strtmp;
  fHch2D_cc_b = 0;

  // 2D Hit profile
  name = fName + "_2D_hitmapPhotonView";
  fHch2D_cc = new TH2F(name.c_str(),name.c_str(),96,0.,96,96,0.,96);
  fHch2D_cc -> SetOption("colz"); 
  fHch2D_cc -> SetContour(99);
  fHch2D_cc -> GetXaxis() -> SetNdivisions(-424); 
  fHch2D_cc -> GetYaxis() -> SetNdivisions(-424); 
  fHch2D_cc -> GetXaxis()-> SetLabelSize(0.03);
  fHch2D_cc -> GetYaxis()-> SetLabelSize(0.03);
  fHch2D_cc -> GetXaxis() -> SetTitle("channel#");
  fHch2D_cc -> GetYaxis() -> SetTitle("channel#");
  fHch2D_cc -> GetYaxis() -> SetTitleOffset(1.2);
  fHch2D_cc -> GetZaxis() -> SetTitle("number of hits");
  fHch2D_cc -> SetStats(kFALSE); 
  AddHistogram(fHch2D_cc);

  // 2D Hit profile - ElectronicView
  name = fName + "_2D_hitmapElectrView";
  fHch2D_cc_b = new TH2F(name.c_str(),name.c_str(),96,0.,96,96,0.,96);
  fHch2D_cc_b -> SetOption("colz"); 
  fHch2D_cc_b -> SetContour(99);
  fHch2D_cc_b -> GetXaxis() -> SetNdivisions(-424); 
  fHch2D_cc_b -> GetYaxis() -> SetNdivisions(-424); 
  fHch2D_cc_b -> GetXaxis()-> SetLabelSize(0.03);
  fHch2D_cc_b -> GetYaxis()-> SetLabelSize(0.03);
  fHch2D_cc_b -> GetXaxis() -> SetTitle("channel#");
  fHch2D_cc_b -> GetYaxis() -> SetTitle("channel#");
  fHch2D_cc_b -> GetYaxis() -> SetTitleOffset(1.2);
  fHch2D_cc_b -> GetZaxis() -> SetTitle("number of hits");
  fHch2D_cc_b -> SetStats(kFALSE); 
  AddHistogram(fHch2D_cc_b);

  // Single events: 2D sum of times histogram
  name = fName + "_2D_SingleEvtPhotView";
  fHchamp2D = new TH2F(name.c_str(),name.c_str(),96,0.,96,96,0.,96);
  fHchamp2D -> SetOption("colz"); 
  fHchamp2D -> GetXaxis() -> SetNdivisions(-424); 
  fHchamp2D -> GetYaxis() -> SetNdivisions(-424); 
  fHchamp2D -> GetXaxis()-> SetLabelSize(0.03);
  fHchamp2D -> GetYaxis()-> SetLabelSize(0.03);
  fHchamp2D -> GetXaxis() -> SetTitle("channel#");
  fHchamp2D -> GetYaxis() -> SetTitle("channel#");
  fHchamp2D -> GetYaxis() -> SetTitleOffset(1.2);
  fHchamp2D -> GetZaxis() -> SetTitle("time (ns)");
  fHchamp2D -> SetStats(kFALSE); 
  //  fHchamp2D -> SetMinimum(1672); // apply 20ns time-cut, 2006 data
  //  fHchamp2D -> SetMaximum(1692); 
  fHchamp2D -> SetMinimum(1546); // apply 10ns time-cut, 2007 data
  fHchamp2D -> SetMaximum(1556); 
//  fHchamp2D -> SetMinimum(1695); // apply 10ns time-cut
//  fHchamp2D -> SetMaximum(1705); 
  AddHistogram(fHchamp2D);
}


void GroupRICH_MAPMT::EndEvent(const CS::DaqEvent &event) {

  if (fNbPlanes<1) return;
  if (thr_flag) TThread::Lock();

  if((fRateCounter%50)){ 
    fHchamp2D -> Reset();
    fHchamp2D -> SetContour(99);
    }


  for (register int ip = 0; ip < fNbPlanes; ip++) {
    const PlaneRICH_MAPMT* plane = dynamic_cast<const PlaneRICH_MAPMT*>(fPlanes[ip]);
    // Get digits (single strip data)
    for (register int idg = 0; idg < plane -> fVrow -> GetNvalues(); idg++) {
      register int row = (int) plane -> fVrow -> GetValues()[idg];
      register int col = (int) plane -> fVcol -> GetValues()[idg];
      register int amp = (int) plane -> fVamp -> GetValues()[idg];
      if (ip == 0) {
         row += 48;
         col += 48;
      }
      else if (ip == 1){
         row += 48;
      }
      else if(ip == 2){
         col += 48;
      }
      else if (ip>3) std::cout << "wrong plane No (quad > 4) "<< std::endl;

      fHch2D_cc -> Fill(col,row);
      if((fRateCounter%50)){ 
      fHchamp2D -> Fill(col,row, -amp); 
      }
    }
      const register TH2F* ptrHch2D_cc_b = plane -> fHrc1b;
      if (fHch2D_cc_b) for (register int ix = 1; ix <= ptrHch2D_cc_b -> GetNbinsX(); ix++) {
        for (register int iy = 1; iy <= ptrHch2D_cc_b -> GetNbinsY(); iy++) {
          register double val = ptrHch2D_cc_b -> GetBinContent(ix, iy);
 	  switch (ip){
	   case 0: fHch2D_cc_b -> SetBinContent(ix,iy+48, val); 
	           break;
	   case 1: fHch2D_cc_b -> SetBinContent(ix+48,iy+48, val);
	           break;
	   case 2: fHch2D_cc_b -> SetBinContent(ix,iy, val);
	           break;
	   case 3: fHch2D_cc_b -> SetBinContent(ix+48,iy, val);
	           break;
	   default: std::cout << "wrong plane No (quad > 4) "<< std::endl;
	           break;
          }      
        }
      }
  }
  if (thr_flag) TThread::UnLock();
}






