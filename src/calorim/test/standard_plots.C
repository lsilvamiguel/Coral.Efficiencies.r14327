// ROOT macro to create standard plots for all calorimeters
// requires TRAFDIC histogram file produced with the option:
//   CsBuildPart	Hist		5
// or for higher resolution, use eg.:
//   CsBuildPart	Hist		9
//
// based on CALposisitions.cxx by Alex Austregesilo
//
// Doesn't work with CINT interpretation mode, but works well with CINT
// compilation mode.  See USAGE message below!
// The ROOT bugs which are the cause of that restriction are supposed to be
// fixed in v5.28, cf. https://savannah.cern.ch/bugs/?64311 and
// https://savannah.cern.ch/bugs/?72472


#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "TCanvas.h"
#include "TCutG.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TString.h"
#include "TStyle.h"
#include "TROOT.h"

using namespace std;

#define DEBUG 0

void standard_plots() {
  cerr << endl;
  cerr << "USAGE (from shell):   root standard_plots.C+\\(\\\"trafdic.root\\\"\\)" << endl;
  cerr << endl;
}


// Calculate rebin factor which leads to target bincount larger or equal to
// "target" and which is a divisor of the existing bincount.
unsigned rebin_factor(unsigned bincount, unsigned target) {
  unsigned r = bincount / target;
  for ( ; bincount % r != 0 && r > 1; r-- ) ;
  return r;
}


class myRow {
public:
  myRow(void) {};
  myRow(TString name_, vector<double> drangex_, vector<double> drangey_) :
    name(name_),
    drangex(drangex_),
    drangey(drangey_),
    fittype(""),
    profile(0),
    bincountx(0),
    bincounty(0),
    xtitle(""),
    ytitle(""),
    mirror(false),
    marker(6)
  {};
  
  TString        name;
  vector<double> drangex;   // histogram display range (x)
  vector<double> drangey;   // histogram display range (y)
  string         fittype;
  vector<double> fitrange;
  vector<double> sigmapreset;
  int            profile;
  int            bincountx; // target bincount (x)
  int            bincounty; // target bincount (y)
  TString        xtitle;
  TString        ytitle;
  bool           mirror;
  int            marker;
};

class myCanvas {
public:
  myCanvas(TString cname, TString ctitle, int ndet=4);

  // add row of histograms to display (with drawing ranges)
  // if a histogram is mirrored, the ranges apply to the mirrored histogram
  // if a profile plot is requested, the ranges define the cut region for computing the profile
  void AddRow(TString name_, double* drangex_=NULL, double* drangey_=NULL);

  // add fit to previous row
  void SetFit(string fittype_, double* fitrange_=NULL, double* sigmapreset_=NULL);

  // add axis titles to previous row
  void SetAxisTitle(TString xtitle="", TString ytitle="");

  // set marker type
  void SetMarker(int marker_);

  // display previous row as a profile instead of a histogram
  // apply re-bin down to ~bincount bins
  void SetProfile(int bincount=0);

  // specify target bincount (X).  rebinning will be done, trying to get as
  // close as possible to the target bincount.
  void SetBinCountX(int bincount);

  // specify target bincount (Y).  rebinning will be done, trying to get as
  // close as possible to the target bincount.
  void SetBinCountY(int bincount);

  // exchange X and Y
  void SetMirror(void);

  // show the previously configured canvas
  void Display(void);
private:
  TCanvas*         canvas;
  TString          name;
  vector< myRow >  rows;
  int              NDet;
};

myCanvas::myCanvas(TString cname, TString ctitle, int ndet) : name(cname), NDet(ndet) {
  canvas = new TCanvas(cname, ctitle, 1200, 900);
}

void
myCanvas::AddRow(TString name_, double* drangex_, double* drangey_) {
  vector<double> drx, dry;
  for (int i=0; i < NDet; i++) {
    if (drangex_)
      drx.push_back(drangex_[i]);
    if (drangey_)
      dry.push_back(drangey_[i]);
  }

  assert( ! (drx.size() == 0 && drangex_ != NULL) );
  assert( ! (drx.size() != 0 && drangex_ == NULL) );
  assert( ! (dry.size() == 0 && drangey_ != NULL) );
  assert( ! (dry.size() != 0 && drangey_ == NULL) );

  rows.push_back(myRow(name_, drx, dry));
}

void
myCanvas::SetFit(string fittype_, double* fitrange_, double* sigmapreset_) {
  vector<double> fr, sp;

  for (int i=0; i < NDet; i++) {
    if (fitrange_)
      fr.push_back(fitrange_[i]);
    if (sigmapreset_)
      sp.push_back(sigmapreset_[i]);
  }

  rows.back().fittype     = fittype_;
  rows.back().fitrange    = fr;
  rows.back().sigmapreset = sp;
}

void
myCanvas::SetAxisTitle(TString xtitle, TString ytitle) {
  rows.back().xtitle = xtitle;
  rows.back().ytitle = ytitle;
}

void
myCanvas::SetMarker(int marker_) {
  rows.back().marker = marker_;
}

void
myCanvas::SetMirror(void) {
  rows.back().mirror = true;
}

void
myCanvas::SetProfile(int bincount) {
  rows.back().profile = bincount;
}

void
myCanvas::SetBinCountX(int bincount) {
  rows.back().bincountx = bincount;
}

void
myCanvas::SetBinCountY(int bincount) {
  rows.back().bincounty = bincount;
}

void
myCanvas::Display(void) {
  canvas->Divide(NDet, rows.size());
  
  // List of all detectors
  TString modlist[] = {"EC01", "EC02", "HC01", "HC02"};


  // Loop over detectors
  TString unique(name);
  for (int m=0; m<NDet; m++) {
    
    // Module name
    TString &module = modlist[m];
    
    // make directory name
    char hdir[200];
    sprintf(hdir, "Calorimeter_%sP1__", (const char *)module);
    
    // change directory
    gDirectory->cd(hdir);

    int i=0;
    for (vector<myRow>::const_iterator row=rows.begin(); row != rows.end(); row++, i++) {
      
      canvas->cd(m+1 + NDet*i);

      // Get histogram
      TH1D       *h1d = NULL;
      TH2D       *h2d = NULL;
      TProfile   *p1d = NULL;
      TProfile2D *p2d = NULL;
      TProfile3D *p3d = NULL;
      TH2        *h2  = NULL;
      TH3        *h3  = NULL;

      TH1 *h  = NULL;
      TH1 *h1 = dynamic_cast<TH1*>( gROOT->FindObject(row->name) );
      if ( !h1 ) continue;
      if (( h3 = dynamic_cast<TH3*>( h1 ) )) {
	p3d = dynamic_cast<TProfile3D*>( h3 );
	assert( p3d );
	h = new TProfile3D( *p3d );
      } else if (( h2 = dynamic_cast<TH2*>( h1 ) )) {
	if (( p2d = dynamic_cast<TProfile2D*>( h2 ) )) {
	  h = new TProfile2D( *p2d );
	} else if (( h2d = dynamic_cast<TH2D*>( h2 ) )) {
	  h = new TH2D( *h2d );
	}
      } else if (( h1d = dynamic_cast<TH1D*>( h1 ) )) {
        if (( p1d = dynamic_cast<TProfile*>( h1d ) )) {
          h = new TProfile( *p1d );
        } else {
          h = new TH1D( *h1d );
        }
      }
      
      if (h == NULL) {
	cerr << "histogram " << i << ":" << row->name << " for "<< module <<" not found" << endl;
	continue;
      }

      if (row->mirror) {
	TH2D *h2 = dynamic_cast<TH2D*>(h);
	assert(h2);
        if (DEBUG)
          cout << "Mirroring: x=" << h2->GetNbinsX() << ", y=" << h2->GetNbinsY() << endl;
        h = new TH2D((TString(h2->GetName())+"_mirror").Data(), h2->GetTitle(),
                     h2->GetNbinsY(), h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax(),
                     h2->GetNbinsX(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
        for (int x=1; x <= h2->GetNbinsX(); x++)
          for (int y=1; y <= h2->GetNbinsY(); y++)
            h->SetBinContent(y, x, h2->GetBinContent(x, y));
      }

      // set ranges
      double xmin = h->GetXaxis()->GetXmin();
      double xmax = h->GetXaxis()->GetXmax();
      double ymin = h->GetYaxis()->GetXmin();
      double ymax = h->GetYaxis()->GetXmax();
      if (DEBUG)
        cout << module << " " << row->name << endl;
      if (row->drangex.size()) {
        xmin = -row->drangex.at(m);
        xmax =  row->drangex.at(m);
        h->GetXaxis()->SetRangeUser(xmin, xmax);
        if (DEBUG)
          cout << "xrange: " << xmin << " ... " << xmax << " (reduced)" << endl;
      } else {
        if (DEBUG)
          cout << "xrange: " << xmin << " ... " << xmax << endl;
      }
      if (row->drangey.size()) {
        ymin = -row->drangey.at(m);
        ymax =  row->drangey.at(m);
        h->GetYaxis()->SetRangeUser(ymin, ymax);
        if (DEBUG)
          cout << "yrange: " << ymin << " ... " << ymax << " (reduced)" << endl;
      } else {
        if (DEBUG)
          cout << "yrange: " << ymin << " ... " << ymax << endl;
      }

      
      // 2d profile from TProfile3D
      if ( h3 ) {
	h2 = p2d = h3->Project3DProfile(TString("yx")+unique); unique += "x";
      }

      // calculate profile
      if (row->profile) {
        double xpoly[4], ypoly[4];
        // bottom left
        xpoly[0] = xmin;
        ypoly[0] = ymin;
        // top left
        xpoly[1] = xmin;
        ypoly[1] = ymax;
        // top right
        xpoly[2] = xmax;
        ypoly[2] = ymax;
        // bottom right
        xpoly[3] = xmax;
        ypoly[3] = ymin;
        TCutG *cutg = new TCutG("cutg", 4, xpoly, ypoly);
        if (DEBUG)
          cutg->Print();

	TH2 *h2 = dynamic_cast<TH2*>(h);
	assert(h2);

	// rebin to a number of bins >= the bincount-goal given in row->profile
	if ( row->profile != 0 ) {
	  unsigned r = rebin_factor( h2->GetNbinsY(), row->profile );
	  if ( r > 1 )
	    h2->RebinY(r);
	}
	h = h2->ProfileY(unique, 1, -1, "[cutg]"); unique += "x";

        // set ranges for profile
        h->GetXaxis()->SetRangeUser(ymin, ymax);
        h->GetYaxis()->SetRangeUser(xmin/2., xmax/2.);

        // plot profile
	h->Draw();

      } else { 
	// no profile calculation

	// rebinning (1D)
	if ( row->bincountx && !h2 ) {
	    // 1D case
	    unsigned r = rebin_factor( h->GetNbinsX(), row->bincountx );
	    if ( r > 1 )
	      h->Rebin(r);
	  }
	
	// rebinning 2D X
	if ( row->bincountx && h2 ) {
	  unsigned r = rebin_factor( h->GetNbinsX(), row->bincountx );
	  if ( r > 1 )
	    h2->RebinX(r);
	}

	// rebinning 2D Y
	if ( row->bincounty && h2 ) {
	  unsigned r = rebin_factor( h->GetNbinsY(), row->bincounty );
	  if ( r > 1 )
	    h2->RebinY(r);
	}
	
	if ( p2d ) {

	  p2d->Draw("colz");
	
	} else {
	  
	  // plot histogram
	  h->SetMarkerStyle(row->marker);
	  h->Draw();
	}
      }


      // fit range
      double fmin, fmax;
      if (row->fitrange.size() == 0) {
	fmin = h->GetXaxis()->GetXmin();
	fmax = h->GetXaxis()->GetXmax();
      } else {
	fmin = - row->fitrange.at(m);
	fmax =   row->fitrange.at(m);
      }

      if (row->fittype == "linear") {
	TF1 *f = new TF1("f", "pol1(0)");
	f->SetLineColor(kRed);
	h->Fit(f, "qrfm", "", fmin, fmax);

	cout << module << setw(20) << row->name
             << "  slope: " << setprecision(3) << setw(5)
             << f->GetParameter(1) * 100. << "%" << endl;
      }

      if (row->fittype == "gauss") {
	
	TF1 *f = new TF1("f", "gaus(0)+pol1(3)");
	
	// Gauss Constant
	f->SetParameter(0, .9*h->GetMaximum());
	f->SetParLimits(0, .001*h->GetMaximum(), h->GetMaximum());
	// Mean
	f->SetParameter(1, 0);
	// Sigma
	f->SetParameter(2, row->sigmapreset.at(m)); 
	f->SetParLimits(2, .01, 1000);
	// Backgr. Constant
	f->SetParameter(3, .1*h->GetMaximum());
	f->SetParLimits(3, -h->GetMaximum(), h->GetMaximum());
	// Backgr. Slope
	f->SetParameter(4, 0);
	
	f->SetLineColor(kRed);
	h->Fit(f, "qrm", "", fmin, fmax);
        double A     = f->GetParameter(0);
	double mean  = f->GetParameter(1);
	double sigma = f->GetParameter(2);
        double a     = f->GetParameter(3);
        double b     = f->GetParameter(4);
        double area  = A / sigma / sqrt(2*TMath::Pi());
        double sb    = area / (a*6*sigma);
	
	cout << module << setw(20) << row->name << "  mean: "
             << setprecision(3) << setw(5) << mean
             << ", sigma: " << setw(5) << sigma
             << ", S/B: " << setprecision(4) << setw(5) << sb << endl;
      }

      if (row->xtitle != "")
        h->GetXaxis()->SetTitle( row->xtitle );
      if (row->ytitle != "")
        h->GetYaxis()->SetTitle( row->ytitle );

    }
    
    gDirectory->cd("..");
  }
}


void
standard_plots(TString file, TString detector="all") {

  TH1::SetDefaultSumw2(kTRUE);

  // Graphics setup
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(10);
  gStyle->SetOptFit(1);
  gStyle->SetMarkerStyle(8);  

  // Open ROOT file
  cout << "Opening file: " << file << endl;
  TFile *fi = new TFile(file);

  // go to head director
  gDirectory->cd("CsBPMonitor");

  
  // definition of the individual canvases
  {
    myCanvas *c = new myCanvas("c0", "times + energy");
    {
      double drangex[] = { 25,  25,  25,  25};
      double frange[]  = { 20,  20,  20,  20};
      double spreset[] = {  2,   2,   2,   2};
      c->AddRow("Time", drangex, NULL);
      c->SetFit("gauss", frange, spreset);
    }
    {
      double drangex[] = { .5,  .5,  .5,  .5};
      double frange[]  = { .3,  .3,  .3,  .3};
      double spreset[] = { .1,  .1,  .1,  .1};
      c->AddRow("deltaEs", drangex, NULL);
      c->SetFit("gauss", frange, spreset);
    }
    c->Display();
  }

  {
    myCanvas *c = new myCanvas("c1", "X position");
    {
      double drangex[] = {100, 100, 200, 200};
      double frange[]  = { 50,  12, 100, 100};
      double spreset[] = { 15,   5,  50,  50};
      c->AddRow("deltaX", drangex, NULL);
      c->SetFit("gauss", frange, spreset);
      c->AddRow("deltaXEcut", drangex, NULL);
      c->SetFit("gauss", frange, spreset);
    }
    c->Display();
  }

  {
    myCanvas *c = new myCanvas("c2", "Y positions");
    {
      double drangex[] = {100, 100, 200, 200};
      double frange[]  = { 50,  12, 100, 100};
      double spreset[] = { 15,   5,  50,  50};
      c->AddRow("deltaY", drangex, NULL);
      c->SetFit("gauss", frange, spreset);
      c->AddRow("deltaYEcut", drangex, NULL);
      c->SetFit("gauss", frange, spreset);
    }
    c->Display();
  }

  {
    myCanvas *c = new myCanvas("c3", "X pitch");
    {
      double drange[] = { 200,   30,  400,  400};
      double frange[] = {1800, 1100, 1900, 2000};
      c->AddRow("deltaXcaloX", NULL, drange);
      c->SetMirror();
      c->SetAxisTitle("track X [mm]", "#DeltaX [mm]");

      c->AddRow("deltaXcaloX", drange, NULL);
      c->SetAxisTitle("track X [mm]", "#DeltaX [mm]");
      c->SetProfile(30);
      c->SetFit("linear", frange);
    }
    c->Display();
  }

  {
    myCanvas *c = new myCanvas("c4", "Y pitch");
    {
      double drange[] = { 200,  50,  400, 400};
      double frange[] = {1000, 600, 1200, 900};
      c->AddRow("deltaYcaloY", NULL, drange);
      c->SetMirror();
      c->SetAxisTitle("track Y [mm]", "#DeltaY [mm]");

      c->AddRow("deltaYcaloY", drange, NULL);
      c->SetAxisTitle("track Y [mm]", "#DeltaY [mm]");
      c->SetProfile(30);
      c->SetFit("linear", frange);
    }
    c->Display();
  }

  {
    myCanvas *c = new myCanvas("c5", "Cluster <--> Track association");
    {
      c->AddRow("effCaloAss");
      c->AddRow("effTrackAss");
    }
    c->Display();
  }

  {
    myCanvas *c = new myCanvas("c6", "Cluster <--> Track association");
    {
      c->AddRow("effCaloAss2D");
      //      c->SetBinCountX(64);
      //      c->SetBinCountY(48);
      c->AddRow("effTrackAss2D");
      //      c->SetBinCountX(64);
      //      c->SetBinCountY(48);
    }
    c->Display();
  }

  for (int i=0; i<8; i++) {
    stringstream s; s << i;
    myCanvas *c = new myCanvas(s.str().c_str(), TString("Cluster <--> Track association (") + s.str() + ")");
    {
      c->AddRow(TString("effTrackAss3D")+s.str());
      //      c->SetBinCountX(64);
      //      c->SetBinCountY(48);
      c->AddRow(TString("effTrackAss3D")+s.str());
      //      c->SetBinCountX(64);
      //      c->SetBinCountY(48);
    }
    c->Display();
  }
}


int main(int argc, char **argv) {
  assert(argc >= 1);

  standard_plots(TString(argv[1]));
  
  return 0;
}
