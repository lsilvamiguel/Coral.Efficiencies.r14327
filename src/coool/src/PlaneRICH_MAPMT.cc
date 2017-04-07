#include "TStyle.h"
#include "PlaneRICH_MAPMT.h"
#include "TF1.h"
#include <sstream>
#include <iostream>

/* ---------------------------------------------------------------------
   Maintainer & author: Frank Nerling 
   New Class for monitoring RICH_MAPMT: (x, y, time), created April 2006
   Updated Version: 11.09.2006 
   --------------------------------------------------------------------- */

using namespace std;

ClassImp(PlaneRICH_MAPMT);


PlaneRICH_MAPMT::PlaneRICH_MAPMT(const char *detname, int ncol, int nrow,
				 int center, int width)
: Plane2V(detname,nrow, ncol, center, width){

}


void PlaneRICH_MAPMT::Init(TTree* tree) {

  Plane2V::Init(tree);
  gStyle -> SetPalette(1);

  std::string strtmp;

// times distribution
  strtmp = fName + "_times";
  fHa -> SetNameTitle(strtmp.c_str(),strtmp.c_str()); 
  fHa -> GetXaxis()->SetTitle("time (ns)");
  ((TH1F_Ref*)fHa)->SetReference(fReferenceDirectory);

  strtmp = fName + "_times_vs_chan#";
  fHavsadr -> SetNameTitle(strtmp.c_str(),strtmp.c_str()); 
  fHavsadr -> GetYaxis()->SetTitle("time (ns)");
  fHavsadr -> GetYaxis()->SetTitleOffset(1.2);
  fHavsadr -> GetXaxis()->SetTitle("channel#");

// further histos:
  strtmp = fName + "_times_vs_chan#cc";
  fHtime_vs_chan = new TH2F(strtmp.c_str(),"times vs. channel#",2304,0,2304,2000,-1800,-1600);
  fHtime_vs_chan -> SetOption("colz");
  fHtime_vs_chan -> GetXaxis()->SetTitle("channel#");
  fHtime_vs_chan -> GetYaxis()->SetTitle("time (ns)");
  fHtime_vs_chan -> GetYaxis()->SetTitleOffset(1.2);
  AddHistogram(fHtime_vs_chan);

// hitmap photon view
  fHrc1 -> GetXaxis() -> SetNdivisions(-412); 
  fHrc1 -> GetYaxis() -> SetNdivisions(-412); 
  fHrc1 -> GetXaxis() -> SetTitle("channel#");
  fHrc1 -> GetYaxis() -> SetTitle("channel#");
  fHrc1 -> GetYaxis() -> SetTitleOffset(1.2);
  fHrc1 -> GetZaxis() -> SetTitle("number of hits");
  fHrc1 -> SetOption("colz"); 
  fHrc1 -> SetContour(99);
  strtmp = fName + "_PhotonView";
  fHrc1 -> SetNameTitle(strtmp.c_str(),strtmp.c_str()); 

// hitmap electronic view
  fHrc1b = (TH2F*) fHrc1 -> Clone();
  strtmp = fName + "_ElectrView";
  fHrc1b -> SetNameTitle(strtmp.c_str(),strtmp.c_str()); 
  fHrc1b -> GetXaxis() -> SetNdivisions(-412); 
  fHrc1b -> GetYaxis() -> SetNdivisions(-412); 
  fHrc1b -> GetXaxis() -> SetTitle("channel#");
  fHrc1b -> GetYaxis() -> SetTitle("channel#");
  fHrc1b -> GetYaxis() -> SetTitleOffset(1.2);
  fHrc1b -> GetZaxis() -> SetTitle("number of hits");
//  fHrc1b -> GetXaxis()-> SetLabelSize(0.03);
//  fHrc1b -> GetYaxis()-> SetLabelSize(0.03);
  AddHistogram(fHrc1b);

// singleEvtPhotView
  fHrca -> GetXaxis() -> SetNdivisions(-412); 
  fHrca -> GetYaxis() -> SetNdivisions(-412); 
  fHrca -> GetXaxis() -> SetTitle("channel#");
  fHrca -> GetYaxis() -> SetTitle("channel#");
  fHrca -> GetYaxis() -> SetTitleOffset(1.2);
  fHrca -> GetZaxis() -> SetTitle("time (ns)");
  fHrca -> SetOption("colz"); 
  fHrca -> SetStats(kFALSE); 
  //  fHrca -> SetMinimum(1672); // apply 20ns time-cut, 2006 data
  //  fHrca -> SetMaximum(1692); 
//  fHrca -> SetMinimum(1650); // apply 100ns time-cut, 2006 data
//  fHrca -> SetMaximum(1750); 
//  fHrca -> SetMinimum(1690); // apply 20ns time-cut
//  fHrca -> SetMaximum(1710); 
//  fHrca -> SetMinimum(1695); // apply 10ns time-cut
//  fHrca -> SetMaximum(1705); 
//-- 2007
  fHrca -> SetMinimum(1546); // apply 10ns time-cut, 2007 data
  fHrca -> SetMaximum(1556); 
  strtmp = fName + "_SingleEvtPhotonView";
  fHrca -> SetNameTitle(strtmp.c_str(),strtmp.c_str()); 
}


PlaneRICH_MAPMT::~PlaneRICH_MAPMT() {
}


void PlaneRICH_MAPMT::StoreDigit(CS::Chip::Digit* digit) {

  Plane2V::StoreDigit (digit);

  int quad_num = atoi( digit -> GetDetID().GetName().substr(5,2).c_str() );

  if( quad_num < 0 || quad_num > 12 ){
  printf("PlaneRICH_MAPMT::StoreDigit(): Bad cathode number%d.\n",quad_num);
  return;
  }

  std::vector<float> data = digit -> GetNtupleData();

  if((data[0] >= 0 && data[0] <= 48) && (data[1] >= 0 && data[1] <= 48)) {
  }
  else     printf("PlaneRICH_MAPMT::StoreDigit(): Bad Channel No %d %d \n",(int) data[0],(int)data[1]);
}


void PlaneRICH_MAPMT::EndEvent(const CS::DaqEvent &event) {

  if (thr_flag) TThread::Lock();

// Single event:
  if((fRateCounter%50)){ 
   fHrca -> Reset();
   fHrca -> SetContour(99);
  }
  for (int i = 0; i < fNhits; i++) {
    int row = fRow[i];
    int col = fCol[i];
    int amp = fAmp[i];
    int chan;

/* Conversion from coral to elctronic view - implemented & checked F.N. 03.08.2006
   (re-rotation of each MAPMT by 180 deg & re-swapping at vertical at centre of Xaxis) */

    int PMTx, PMTy;
    int Xpmt, Ypmt;
    int row_b, col_b;  

// Redo the optical mirroring:
    PMTy = (int) row / 4;  		
    PMTx = (int) col / 4;  		
    assert( PMTy < 12  && PMTy >= 0 && PMTx < 12  && PMTx >= 0);

    Ypmt = row - 4*PMTy;	
    Xpmt = col - 4*PMTx;	
    assert( Ypmt < 4  && Ypmt >= 0 && Xpmt < 4  && Xpmt >= 0);

    Ypmt = 3 - Ypmt;		
    Xpmt = 3 - Xpmt;		
    assert( Ypmt < 4  && Ypmt >= 0 && Xpmt < 4  && Xpmt >= 0);
    
    row_b = Ypmt + 4*PMTy;	
    col_b = Xpmt + 4*PMTx;
    assert( row_b < 48  && row_b >= 0 && col_b < 48  && col_b >= 0);

// Swapping of Xcoordiante -> to go to coral system of reference (photonview):
    col_b = 47 - col_b;
    assert( col_b < 48  && col_b >= 0);
// End of conversion from coral to elctronic view - implemented & checked F.N. 03.08.2006

    if(fVrow -> Test(row) &&
       fVcol -> Test(col) &&
       fVamp -> Test(amp) &&
       fVrow -> Test(row_b) &&
       fVcol -> Test(col_b)) {    

      fHa -> Fill(amp);
      fHrc1 -> Fill(col,row);
      fHrc1b -> Fill(col_b,row_b);
// Single event:
      if((fRateCounter%50)){ fHrca -> Fill(col,row, -amp); }
      
      chan = row * 48 + col;
      assert( chan < 2304  && chan >= 0 );

      fHtime_vs_chan -> Fill(chan,amp);
      fHavsadr -> Fill(chan,amp);

      fVrow -> Store(row);
      fVcol -> Store(col);
      fVamp -> Store(amp);
       
      fNhitsKept++;
    }
  }
  fHhit -> Fill(fNhitsKept);
  if (thr_flag) TThread::UnLock();
}








