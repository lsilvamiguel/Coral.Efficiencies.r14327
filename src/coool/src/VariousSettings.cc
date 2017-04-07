#include <unistd.h>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>
#include <iostream>
#include "TCanvas.h"
#include "TPaveStats.h"
#include "../expat/xmlparse/xmlparse.h"

#include "VariousSettings.h"
#include "TStyle.h"

#include <iostream>
#include <streambuf>
#ifndef __CINT_
  #include <dirent.h>
  #include <fnmatch.h>
#endif

ClassImp(VariousSettings)

const char* ORDERLY_SHIFT_PLOTS_DIR    = "/online/daq/shift/orderly_shift_plots";
const char* SHIFT_DIR                  = "/online/daq/shift";
const char* ORDERLY_SHIFT_PLOTS_FNAME  = "orderly_shift_plots.ps";

const char* VariousSettings::fPersoFile = ".cooolrc";
int VariousSettings::fMaxItemInList = 25;


namespace {  // private namespace for save functions

Bool_t SaveGuiSettings(ofstream& guisetfile)
{
   // Save the settings of the monitoring GUI
   // For the moment it saves principally the names of the histos in
   // all the canvas windows
   // NB: the file format is still preliminary !!!!

  register TObject* cvo;
  register TSeqCollection *allCanvas;
  TDatime temps;

  temps.Set();
  guisetfile << "<guisetup version=\"0.4\" date=\"" << temps.AsSQLString();
  guisetfile << "\">\n";
  allCanvas = gROOT->GetListOfCanvases();
  TIter nextcv(allCanvas);

  while ((cvo = (TObject*) nextcv())) {
    register TCanvas* cv;
    if ((cv = (TCanvas*)(cvo))) {
      std::string cvname(cv->GetName());
      if (cvname == "fCanvasWindow") continue;  // to not record content of main window

      register TPaveStats* pv;
      pv = (TPaveStats*) cv->FindObject("stats");
      register TObject* p;
      TList *subpads=cv->GetListOfPrimitives();
      TIter next(subpads);
      guisetfile << "  <canvas name=\"" << cv->GetName() << "\" geometry=\"";
      guisetfile << cv->GetWindowWidth() << "x" << cv->GetWindowHeight();
      guisetfile << "+" << cv->GetWindowTopX() << "+" << cv->GetWindowTopY() << "\"";
      guisetfile << " logXYZ=\"" << cv->GetLogx() << cv->GetLogy() << cv->GetLogz() <<"\"";
      if (pv && pv->GetOptStat()) {
        guisetfile << " ostat=\"" << pv->GetOptStat() <<"\"";
        guisetfile << " ostatfmt=\"" << pv->GetStatFormat() <<"\"";
        guisetfile << " ostatXY=\"" << pv->GetX1NDC() << " " << pv->GetY1NDC();
        guisetfile << " " << pv->GetX2NDC() << " " << pv->GetY2NDC() <<"\"";
      }
      else guisetfile << " ostat=\"\"";
      guisetfile << ">\n";

      // loop on pads inside the canvas
      while ((p=(TObject*)next())) {
	register TPad *pad=dynamic_cast<TPad*>(p);
	if(pad!=NULL) {
	  register TObject* q;
          Double_t xlow, ylow, xup, yup;
          pad->GetPadPar(xlow, ylow, xup, yup);
          pv = (TPaveStats*) pad->FindObject("stats");
	  TIter nextpd(pad->GetListOfPrimitives());
          guisetfile << "    <pad name=\"" << pad->GetName() << "\" bounds=\"";
          guisetfile << xlow << " " << ylow << " ";
          guisetfile << xup << " " << yup << "\"";
          guisetfile << " logXYZ=\"" << pad->GetLogx() << pad->GetLogy() << pad->GetLogz() <<"\"";
          if (pv && pv->GetOptStat()) {
            guisetfile << " ostat=\"" << pv->GetOptStat() <<"\"";
            guisetfile << " ostatfmt=\"" << pv->GetStatFormat() <<"\"";
            guisetfile << " ostatXY=\"" << pv->GetX1NDC() << " " << pv->GetY1NDC();
            guisetfile << " " << pv->GetX2NDC() << " " << pv->GetY2NDC() <<"\"";
          }
          else guisetfile << " ostat=\"\"";
          guisetfile << ">\n";
	  while ((q = (TObject*) nextpd())) {
	    register RefHist *leref=dynamic_cast<RefHist*>(q);
            if (leref) continue;  // no ref plot there
	    if (q->IsA()->InheritsFrom("TH1")) {  // it is an histo
              register TH1 *leh1 = dynamic_cast<TH1*>(q);
              register TAxis* ax;
              ax = leh1->GetXaxis();
              guisetfile << "      <histo name=\"" << leh1->GetName() << "\"";
              guisetfile << " Xbinrange=\"" << ax->GetFirst() << " " << ax->GetLast() << "\"";
              if (leh1->GetDimension() > 1) {
                ax = leh1->GetYaxis();
                guisetfile << " Ybinrange=\"" << ax->GetFirst() << " " << ax->GetLast() << "\"";
              }
              if (leh1->GetDimension() > 2) {
                ax = leh1->GetZaxis();
                guisetfile << " Zbinrange=\"" << ax->GetFirst() << " " << ax->GetLast() << "\"";
              }
              guisetfile << " options=\"" << leh1->GetOption() << "\"/>\n";
	    }     // NB: we don't loop recursively on subsubpads (may be later...)
	  }
          guisetfile << "    </pad>\n";
	}

	// histos in the canvas directly
	register RefHist *leref=dynamic_cast<RefHist*>(p);
	if (!leref) {
	  register TH1 *leh1=dynamic_cast<TH1*>(p);
          register TAxis* ax;
	  if (leh1!=NULL) {
            ax = leh1->GetXaxis();
            guisetfile << "      <histo name=\"" << leh1->GetName() << "\"";
            guisetfile << " Xbinrange=\"" << ax->GetFirst() << " " << ax->GetLast() << "\"";
            if (leh1->GetDimension() > 1) {
              ax = leh1->GetYaxis();
              guisetfile << " Ybinrange=\"" << ax->GetFirst() << " " << ax->GetLast() << "\"";
            }
            if (leh1->GetDimension() > 2) {
              ax = leh1->GetZaxis();
              guisetfile << " Zbinrange=\"" << ax->GetFirst() << " " << ax->GetLast() << "\"";
            }
            guisetfile << " options=\"" << leh1->GetOption() << "\"/>\n";
	  }
        }
      }
      guisetfile << "  </canvas>\n";
    }
  }

  guisetfile << "</guisetup>\n\n";
  return true;
}


Bool_t SavePlanesSettings(ofstream& guisetfile, Monitor* monitor)
{
   // Save the settings of the planes (histos size and cuts)

  TDatime temps;

  temps.Set();
  guisetfile << "<planesetup version=\"0.1\" date=\"" << temps.AsSQLString();
  guisetfile << "\">\n";

  typedef std::map<std::string,Plane*>::iterator LI;
  typedef std::vector<Variable*>::iterator VI;
  typedef std::vector<TH1*>::iterator HI;

  std::map<std::string,Plane*>& detmap=monitor->GetDetMap();
  for(LI i=detmap.begin();i!=detmap.end();i++) {
    register Plane* pln;
    pln = i->second;
    guisetfile << "  <plane name=\"" << pln->GetName() << "\" type=\"";
    guisetfile << pln->IsA()->GetName() << "\">\n";

    std::vector<Variable*>& varsvec = pln->GetVariables();
    for (VI j = varsvec.begin(); j != varsvec.end(); j++) {
      register Variable* var;
      var = *j;
      guisetfile << "    <variable name=\"" << var->GetName();
      guisetfile << "\" min=\"" << var->GetMin();
      guisetfile << "\" max=\"" << var->GetMax() << "\"/>\n";
    }

    std::vector<TH1*>& histosvec = pln->GetHistoList();
    for (HI k = histosvec.begin(); k != histosvec.end(); k++) {
      register TH1* leh1;
      register TAxis* ax;
      leh1 = *k;
      ax = leh1->GetXaxis();
      guisetfile << "    <histo name=\"" << leh1->GetName();
      guisetfile << "\" type=\"" << leh1->IsA()->GetName() << "\"";
      guisetfile << " Xbinminmax=\"" << ax->GetNbins() << " ";
      guisetfile << ax->GetXmin() << " ";
      guisetfile << ax->GetXmax() << "\"";
      if (leh1->GetDimension() > 1) {
        ax = leh1->GetYaxis();
        guisetfile << " Ybinminmax=\"" << ax->GetNbins() << " ";
        guisetfile << ax->GetXmin() << " ";
        guisetfile << ax->GetXmax() << "\"";
      }
      if (leh1->GetDimension() > 2) {
        ax = leh1->GetZaxis();
        guisetfile << " Zbinminmax=\"" << ax->GetNbins() << " ";
        guisetfile << ax->GetXmin() << " ";
        guisetfile << ax->GetXmax() << "\"";
      }
      guisetfile << "/>\n";
    }

    guisetfile << "  </plane>\n";
  }


  guisetfile << "</planesetup>\n\n";
  return true;
}


Bool_t SaveGroupsSettings(ofstream& guisetfile, Monitor* monitor)
{
   // Save the settings of the planes (histos size and cuts)

  TDatime temps;

  temps.Set();
  guisetfile << "<groupsetup version=\"0.1\" date=\"" << temps.AsSQLString();
  guisetfile << "\">\n";

  typedef std::map<std::string,Group*>::iterator LI;
//   typedef vector<Variable*>::iterator VI;
  typedef std::vector<TH1*>::iterator HI;

  std::map<std::string,Group*>& groupmap=monitor->GetGroupMap();
  for(LI i=groupmap.begin();i!=groupmap.end();i++) {
    register Group* grp;
    grp = i->second;
    guisetfile << "  <group name=\"" << grp->GetName() << "\" type=\"";
    guisetfile << grp->IsA()->GetName() << "\">\n";

//     vector<Variable*>& varsvec = pln->GetVariables();
//     for (VI j = varsvec.begin(); j != varsvec.end(); j++) {
//       register Variable* var;
//       var = *j;
//       guisetfile << "    <variable name=\"" << var->GetName();
//       guisetfile << "\" min=\"" << var->GetMin();
//       guisetfile << "\" max=\"" << var->GetMax() << "\"/>\n";
//     }

    std::vector<TH1*>& histosvec = grp->GetHistoList();
    for (HI k = histosvec.begin(); k != histosvec.end(); k++) {
      register TH1* leh1;
      register TAxis* ax;
      leh1 = *k;
      ax = leh1->GetXaxis();
      guisetfile << "    <histo name=\"" << leh1->GetName();
      guisetfile << "\" type=\"" << leh1->IsA()->GetName() << "\"";
      guisetfile << " Xbinminmax=\"" << ax->GetNbins() << " ";
      guisetfile << ax->GetXmin() << " ";
      guisetfile << ax->GetXmax() << "\"";
      if (leh1->GetDimension() > 1) {
        ax = leh1->GetYaxis();
        guisetfile << " Ybinminmax=\"" << ax->GetNbins() << " ";
        guisetfile << ax->GetXmin() << " ";
        guisetfile << ax->GetXmax() << "\"";
      }
      if (leh1->GetDimension() > 2) {
        ax = leh1->GetZaxis();
        guisetfile << " Zbinminmax=\"" << ax->GetNbins() << " ";
        guisetfile << ax->GetXmin() << " ";
        guisetfile << ax->GetXmax() << "\"";
      }
      guisetfile << "/>\n";
    }

    guisetfile << "  </group>\n";
  }

  guisetfile << "</groupsetup>\n\n";
  return true;
}

} // end of private namespace



Bool_t VariousSettings::SaveConfigSettings(const char* filename, int type_save)
{
   // Save all the settings of rtree

  if (!filename) return false;
  ofstream guisetfile(filename);
  if (!guisetfile)  {
#ifndef __CINT__
    CS::Exception("Can not write GUI settings file.\n");
#endif
    return false;
  }

  guisetfile << "<general_setup>\n\n";
  if (type_save == 0 || type_save == 1) SaveGuiSettings(guisetfile);
  if (type_save == 0 || type_save == 2) {
    SavePlanesSettings(guisetfile, fMonitor);
    SaveGroupsSettings(guisetfile, fMonitor);
  }
  guisetfile << "\n</general_setup>\n";
  guisetfile.close();

  if (type_save == 0) {
    fConfigFileName = filename;
  }
  AddNameToList(filename);
  return true;
}



// ------------------------------------------------------------------------

// Read-out code

namespace {  // private namespace for this file

class UserData
{
  public:
            UserData(void)
                { Clear(); depth = 0; psfile = 0; pspage = 0; monitor = 0;
                  filename = 0; datestr = ""; runnb = 0; evtnb = 0; }
    void    Clear(void)
                { canvasname=""; padname=""; histoname=""; options="";
                  planename = ""; planetype = ""; plane = 0;
                  version=1; canvas = 0; }
    void    ClearCanvas(void)
                { canvasname=""; canvasgeom="570x332+300+300"; canvaslog="000";
                  canvas = 0; padnb = 0; skipcanvas = false; ClearPad(); }
    void    ClearPad(void)
                { padname=""; padbounds="0.01 0.01 0.99 0.99"; padlog="000";
                  pad = 0; ClearHisto(); ClearPave(); ClearTitPave();}
    void    ClearHisto(void)
                { histoname="";
                  histoXbinminmax=""; histoYbinminmax="";
                  histoZbinminmax="";
                  histoXbinrange=""; histoYbinrange=""; histoZbinrange="";
                  histooptions=""; }
    void    ClearPave()
                { paveostat = ""; paveostatfmt = ""; paveostatXY = ""; }
    void    ClearTitPave()
                { pavetitle = ""; pavetitleXY = ""; }
    void    ClearPlane(void)
                { planename = ""; planetype = ""; plane = 0; }
    void    ClearGroup(void)
                { groupname = ""; grouptype = ""; group = 0; }
    void    ClearVar(void)
                { varname = ""; varmin = "0"; varmax = "100"; }
    Monitor* monitor;
    int     read_mode;
    int     type_read;
    std::string  canvasname;
    std::string  canvasgeom;
    std::string  canvaslog;
    TCanvas* canvas;
    bool    skipcanvas;
    int     padnb;
    std::string  padname;
    std::string  padbounds;
    std::string  padlog;
    TPad*   pad;
    std::string  paveostat;
    std::string  paveostatfmt;
    std::string  paveostatXY;
    std::string  pavetitle;
    std::string  pavetitleXY;
    std::string  histoname;
    std::string  histoXbinminmax;
    std::string  histoYbinminmax;
    std::string  histoZbinminmax;
    std::string  histoXbinrange;
    std::string  histoYbinrange;
    std::string  histoZbinrange;
    std::string  histooptions;
    std::string  planename;
    std::string  planetype;
    Plane*  plane;
    std::string  groupname;
    std::string  grouptype;
    Group*  group;
    std::string  varname;
    std::string  varmin;
    std::string  varmax;
    int     version;
    int     depth;
    std::string  options;
    TPostScript* psfile;
    int     pspage;
    const char*  filename;
    std::string  datestr;
    unsigned int runnb;
    int     evtnb;
};


void decodePlanesSettings(UserData& d, std::string& stname, const char** atts)
{
  // Decoding of the plane header
  if (stname == "plane") {
    register int i;
    d.ClearPlane();
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "name") d.planename = val;
      if (opt == "type") d.planetype = val;
    }
    if (d.planename != "") {
      typedef std::map<std::string,Plane*>::iterator DMI;
      DMI dmiit = (d.monitor->GetDetMap()).find(d.planename);
      if (dmiit != (d.monitor->GetDetMap()).end()) { d.plane = dmiit->second; }
        else { d.plane = 0; }
    }
    if (d.plane) d.plane->ResetHistograms();
  }

  // Decoding of the Variable entry
  if (d.plane && (stname == "variable")) {
    register int i;
    Variable* var = 0;
    d.ClearVar();
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "name") d.varname = val;
      if (opt == "min") d.varmin = val;
      if (opt == "max") d.varmax = val;
    }
    typedef std::vector<Variable*>::iterator VI;
    std::vector<Variable*>& varsvec = d.plane->GetVariables();
    if (d.varname != "") for (VI j=varsvec.begin(); j!=varsvec.end(); j++) {
      if ((*j)->GetName() == d.varname) var = *j;
    }
    if (var == 0) {
      std::cerr << "Unknown Variable " << d.varname;
      std::cerr << " in Plane " << d.planename << std::endl;
    }
    else if ((d.varmin != "") && (d.varmax != "")) {
      int vmin, vmax;
      sscanf(d.varmin.c_str(), "%d", &vmin);
      sscanf(d.varmax.c_str(), "%d", &vmax);
      var->SetRange(vmin, vmax);
    }
  }

  // Decoding of the histo entry
  if (d.plane && (stname == "histo")) {
    register int i;
    d.ClearHisto();
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "name") d.histoname = val;
      if (opt == "Xbinminmax") d.histoXbinminmax = val;
      if (opt == "Ybinminmax") d.histoYbinminmax = val;
      if (opt == "Zbinminmax") d.histoZbinminmax = val;
    }
    register TH1* histo = 0;
    typedef std::vector<TH1*>::iterator HI;
    std::vector<TH1*>& histosvec = d.plane->GetHistoList();
    if (d.histoname != "") for (HI j=histosvec.begin(); j!=histosvec.end(); j++) {
      if ((*j)->GetName() == d.histoname) histo = *j;
    }
    if (histo == 0) {
      std::cerr << "Histo " << d.histoname << " not found in Plane " << d.planename <<std::endl;
    }
    else {
      int nbdim = histo->GetDimension();
      switch (nbdim) {
        case 1:
          if (d.histoXbinminmax != "") {
            int xbin;
            float xmin, xmax;
            sscanf(d.histoXbinminmax.c_str(), "%d %f %f", &xbin, &xmin, &xmax);
            histo->SetBins(xbin,xmin,xmax);
          }
        case 2:
          if (d.histoXbinminmax != "" && d.histoYbinminmax != "") {
            int xbin, ybin;
            float xmin, xmax, ymin, ymax;
            sscanf(d.histoXbinminmax.c_str(), "%d %f %f", &xbin, &xmin, &xmax);
            sscanf(d.histoYbinminmax.c_str(), "%d %f %f", &ybin, &ymin, &ymax);
            histo->SetBins(xbin, xmin, xmax, ybin, ymin, ymax);
          }
        case 3:
          if (d.histoXbinminmax != "" && d.histoYbinminmax != "" && d.histoZbinminmax != "") {
            int xbin, ybin, zbin;
            float xmin, xmax, ymin, ymax, zmin, zmax;
            sscanf(d.histoXbinminmax.c_str(), "%d %f %f", &xbin, &xmin, &xmax);
            sscanf(d.histoYbinminmax.c_str(), "%d %f %f", &ybin, &ymin, &ymax);
            sscanf(d.histoZbinminmax.c_str(), "%d %f %f", &zbin, &zmin, &zmax);
            histo->SetBins(xbin, xmin, xmax, ybin, ymin, ymax, zbin, zmin, zmax);
          }
      }
      histo->Reset();
    }
  }
}


void decodeGroupsSettings(UserData& d, std::string& stname, const char** atts)
{
  // Decoding of the group header
  if (stname == "group") {
    register int i;
    d.ClearGroup();
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "name") d.groupname = val;
      if (opt == "type") d.grouptype = val;
    }
    if (d.groupname != "") {
      typedef std::map<std::string,Group*>::iterator DGI;
      DGI dgiit = (d.monitor->GetGroupMap()).find(d.groupname);
      if (dgiit != (d.monitor->GetGroupMap()).end()) { d.group = dgiit->second; }
        else { d.group = 0; }
    }
  }

  // Decoding of the Variable entry
//   if (d.group && (stname == "variable")) {
//     register int i;
//     Variable* var = 0;
//     d.ClearVar();
//     for (i=0; atts[i]!=0; i+=2 ) {
//       std::string opt(atts[i]), val(atts[i+1]);
//       if (opt == "name") d.varname = val;
//       if (opt == "min") d.varmin = val;
//       if (opt == "max") d.varmax = val;
//     }
//     typedef vector<Variable*>::iterator VI;
//     vector<Variable*>& varsvec = d.plane->GetVariables();
//     if (d.varname != "") for (VI j=varsvec.begin(); j!=varsvec.end(); j++) {
//       if ((*j)->GetName() == d.varname) var = *j;
//     }
//     if (var == 0) {
//       std::cerr << "Unknown Variable " << d.varname;
//       std::cerr << " in Plane " << d.planename << std::endl;
//     }
//     else if ((d.varmin != "") && (d.varmax != "")) {
//       int vmin, vmax;
//       sscanf(d.varmin.c_str(), "%d", &vmin);
//       sscanf(d.varmax.c_str(), "%d", &vmax);
//       var->SetRange(vmin, vmax);
//     }
//   }

  // Decoding of the histo entry
  if (d.group && (stname == "histo")) {
    register int i;
    d.ClearHisto();
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "name") d.histoname = val;
      if (opt == "Xbinminmax") d.histoXbinminmax = val;
      if (opt == "Ybinminmax") d.histoYbinminmax = val;
      if (opt == "Zbinminmax") d.histoZbinminmax = val;
    }
    register TH1* histo = 0;
    typedef std::vector<TH1*>::iterator HI;
    std::vector<TH1*>& histosvec = d.group->GetHistoList();
    if (d.histoname != "") for (HI j=histosvec.begin(); j!=histosvec.end(); j++) {
      if ((*j)->GetName() == d.histoname) histo = *j;
    }
    if (histo == 0) {
      std::cerr << "Histo " << d.histoname << " not found in Group " << d.groupname <<std::endl;
    }
    else {
      int nbdim = histo->GetDimension();
      switch (nbdim) {
        case 1:
          if (d.histoXbinminmax != "") {
            int xbin;
            float xmin, xmax;
            sscanf(d.histoXbinminmax.c_str(), "%d %f %f", &xbin, &xmin, &xmax);
            histo->SetBins(xbin,xmin,xmax);
          }
        case 2:
          if (d.histoXbinminmax != "" && d.histoYbinminmax != "") {
            int xbin, ybin;
            float xmin, xmax, ymin, ymax;
            sscanf(d.histoXbinminmax.c_str(), "%d %f %f", &xbin, &xmin, &xmax);
            sscanf(d.histoYbinminmax.c_str(), "%d %f %f", &ybin, &ymin, &ymax);
            histo->SetBins(xbin, xmin, xmax, ybin, ymin, ymax);
          }
        case 3:
          if (d.histoXbinminmax != "" && d.histoYbinminmax != "" && d.histoZbinminmax != "") {
            int xbin, ybin, zbin;
            float xmin, xmax, ymin, ymax, zmin, zmax;
            sscanf(d.histoXbinminmax.c_str(), "%d %f %f", &xbin, &xmin, &xmax);
            sscanf(d.histoYbinminmax.c_str(), "%d %f %f", &ybin, &ymin, &ymax);
            sscanf(d.histoZbinminmax.c_str(), "%d %f %f", &zbin, &zmin, &zmax);
            histo->SetBins(xbin, xmin, xmax, ybin, ymin, ymax, zbin, zmin, zmax);
          }
      }
      histo->Reset();
    }
  }
}


void decodeGuiSettings(UserData& d, std::string& stname, const char** atts)
{
  // Decoding of the canvas header
  if (stname == "canvas") {
    d.ClearCanvas();
    d.depth++;
    register int i;
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "name") d.canvasname = val;
      if (opt == "geometry") d.canvasgeom = val;
      if (opt == "logXYZ") d.canvaslog = val;
      if (opt == "ostat") d.paveostat = val;
      if (opt == "ostatfmt") d.paveostatfmt = val;
      if (opt == "ostatXY") d.paveostatXY = val;
      if (opt == "title") d.pavetitle = val;
      if (opt == "titleXY") d.pavetitleXY = val;
    }
    if (d.canvasname == "") {
      std::cerr << "XML:decodeGuiSettings(): a canvas has no name\n";
      return;
    }

    int wx, wy, px, py;
    sscanf(d.canvasgeom.c_str(), "%dx%d+%d+%d", &wx, &wy, &px, &py);
    register TCanvas *cv = 0;
    std::string vcanvasname = d.canvasname;
    if (d.canvasname == "fChanDisplay") { // we don't touch the panel window with wire time spectrums
      d.skipcanvas = true;
      return;
    }
    if (d.psfile && (d.canvasname == "fCanvasWindow")) { // we don't print the main window
      d.skipcanvas = true;
      return;
    }
    if (d.psfile) d.canvasname = std::string("PStmpcanvasname_") + d.canvasname;
    cv = dynamic_cast<TCanvas*>(gROOT->FindObject(d.canvasname.c_str()));
    if (cv) {
      cv->Clear();
    } else {
      if (d.psfile) cv = new TCanvas("", d.canvasname.c_str(), 290, 200);
      else cv = new TCanvas(d.canvasname.c_str(), d.canvasname.c_str(), px, py, wx, wy);
    }
    if (d.psfile) cv->SetBatch(kTRUE);  // useless as canvas has no name, but...
    d.canvas = cv;
    d.pad = (TPad*) cv;
    if (d.psfile) {
      char strtmp[500];
      const char *lastname = d.filename;
      if ((lastname = strrchr(d.filename, '/'))) lastname++;
        else lastname = d.filename;
      d.psfile->On();
      d.psfile->NewPage();
      d.pspage++;
      sprintf(strtmp, "Run %d %s: %s - Page %d (%d evts)  %s", d.runnb, vcanvasname.c_str(), lastname, d.pspage, d.evtnb, d.datestr.c_str());
      d.psfile->PrintStr(" 0.865 1.02 Zone@");
      d.psfile->Range(29, 19);
      d.psfile->SetTextAlign(23);
      d.psfile->SetTextSize(0.025);
      if (d.filename) d.psfile->Text(0.5, 1, strtmp);
      d.psfile->Off();
    }
    cv->cd();
    register const char* logXYZ = d.canvaslog.c_str();
    if (logXYZ[0] == '1') cv->SetLogx(1); else cv->SetLogx(0);
    if (logXYZ[1] == '1') cv->SetLogy(1); else cv->SetLogy(0);
    if (logXYZ[2] == '1') cv->SetLogz(1); else cv->SetLogz(0);
    if (d.canvasname == "fCanvasWindow") return;
    cv->SetWindowSize(wx, wy);
    cv->SetWindowPosition(px, py);
    cv->Modified();
    if (!d.psfile) return;
    cv->Draw();
    return;
  }


  // Decoding of the pad header
  if (stname == "pad" && d.read_mode == 1) {
    d.ClearPad();
    d.depth++;
    if (d.skipcanvas) return;
    d.padnb++;
    if (d.canvas == 0) {
      std::cerr << "XML:decodeGuiSettings(): a pad is defined not within a canvas\n";
      return;
    }
    register int i;
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "name") d.padname = val;
      if (opt == "bounds") d.padbounds = val;
      if (opt == "logXYZ") d.padlog = val;
      if (opt == "ostat") d.paveostat = val;
      if (opt == "ostatfmt") d.paveostatfmt = val;
      if (opt == "ostatXY") d.paveostatXY = val;
      if (opt == "title") d.pavetitle = val;
      if (opt == "titleXY") d.pavetitleXY = val;
    }
    float lx, ly, ux, uy;
    sscanf(d.padbounds.c_str(), "%f %f %f %f", &lx, &ly, &ux, &uy);
//     fprintf(stderr, "pad padnb %d: lx %f ly %f ux %f uy %f\n", d.padnb, lx, ly, ux ,uy);
    if (d.psfile) d.padname = std::string("PStmppadname_") + d.padname;
    register TPad *pad = 0;
    pad = dynamic_cast<TPad*>(gROOT->FindObject(d.padname.c_str()));
    if (pad) { pad->Clear(); SafeDelete(pad); }
    d.canvas->cd();
    if (d.psfile) { ly *= 0.975; uy *= 0.975; }  // to give some space to the title
    pad = new TPad(d.padname.c_str(), d.padname.c_str(), lx, ly, ux, uy);
    register const char* logXYZ = d.padlog.c_str();
    if (logXYZ[0] == '1') pad->SetLogx(1); else pad->SetLogx(0);
    if (logXYZ[1] == '1') pad->SetLogy(1); else pad->SetLogy(0);
    if (logXYZ[2] == '1') pad->SetLogz(1); else pad->SetLogz(0);
    pad->SetNumber(d.padnb);
    d.pad = pad;
    pad->Modified();
    pad->Draw();
    if (!d.psfile) return;
    pad->Update();
    pad->Draw();
    return;
  }

  // Decoding the histo entry
  if (stname == "histo") {
    d.ClearHisto();
    if (d.skipcanvas) return;
    if (d.pad == 0) {
      CS::Exception("XML:decodeGuiSettings(): an histo is defined not within a canvas or a pad\n");
      std::cerr << "XML:decodeGuiSettings(): an histo is defined not within a canvas or a pad\n";
      return;
    }
    register int i;
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "name") d.histoname = val;
      if (opt == "Xbinrange") d.histoXbinrange = val;
      if (opt == "Ybinrange") d.histoYbinrange = val;
      if (opt == "Zbinrange") d.histoZbinrange = val;
      if (opt == "options") d.histooptions = val;
    }
    register TH1* histo;
    TFile* rootfile = d.monitor->GetRootFile();
    if (rootfile) {
      rootfile->cd();
      histo = dynamic_cast<TH1*>(rootfile->FindObjectAny(d.histoname.c_str()));
      if (!histo) {
        histo = d.monitor->GetHisto(d.histoname.c_str());
      }
    } else {
      histo = dynamic_cast<TH1*>(gROOT->FindObject(d.histoname.c_str()));
    }
    if (histo == 0) {
      CS::Exception("XML:decodeGuiSettings(): the histo %s is not known\n", d.histoname.c_str());
      fprintf(stderr, "XML:decodeGuiSettings(): the histo %s is not known\n", d.histoname.c_str());
      return;
    }
    d.pad->cd();
    register TAxis* ax;
    int first, last;
    if (d.histoXbinrange != "") {
      ax = histo->GetXaxis();
      sscanf(d.histoXbinrange.c_str(), "%d %d", &first, &last);
      ax->SetRange(first, last);
    }
    if ((histo->GetDimension() > 1) && (d.histoYbinrange != "")) {
      ax = histo->GetYaxis();
      sscanf(d.histoYbinrange.c_str(), "%d %d", &first, &last);
      ax->SetRange(first, last);
    }
    if ((histo->GetDimension() > 2) && (d.histoZbinrange != "")) {
      ax = histo->GetZaxis();
      sscanf(d.histoZbinrange.c_str(), "%d %d", &first, &last);
      ax->SetRange(first, last);
    }
    histo->Draw(d.histooptions.c_str());
    d.pad->Modified();
    d.pad->Update();
//     if (d.psfile) d.canvas->Update();
    return;
  }
}




// XML_StartElementHandler
void startElement(void *userData, const char *name, const char **atts)
{
  UserData &d = *reinterpret_cast<UserData*>(userData);
  std::string stname(name);

  for( register int i=0; atts[i]!=0; i+=2 ) {
    if (atts[i+1]==0 )
    {
      std::cerr << "XML:startElement(): an attribute is zero!!";
      return;
    }
  }

  // Decoding of the guisetup header
  if (stname == "guisetup") {
    d.Clear();
    d.depth=0;
    if ((d.type_read == 0) || (d.type_read == 1)) d.read_mode = 1;
      else d.read_mode = 0;
    register int i;
    register double vers = 0;
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "version") vers = atof(atts[i+1]);
    }
    if (vers == 0) {
      std::cerr << "XML:startElement(): no version number in this file\n";
    }
    if (vers < 0.2) {
      std::cerr << "XML:startElement(): this file may not be read (bad version number)\n";
    }
    return;
  }

  // Decoding of the planesetup header
  if (stname == "planesetup") {
    d.Clear();
    d.depth=0;
    if ((d.type_read == 0) || (d.type_read == 2)) d.read_mode = 2;
      else d.read_mode = 0;
    register int i;
    register double vers = 0;
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "version") vers = atof(atts[i+1]);
    }
    if (vers == 0) {
      std::cerr << "XML:startElement(): no version number in this file\n";
    }
    if (vers < 0.1) {
      std::cerr << "XML:startElement(): this file may not be read (bad version number)\n";
    }
    return;
  }

  // Decoding of the groupsetup header
  if (stname == "groupsetup") {
    d.Clear();
    d.depth=0;
    if ((d.type_read == 0) || (d.type_read == 2)) d.read_mode = 3;
      else d.read_mode = 0;
    register int i;
    register double vers = 0;
    for (i=0; atts[i]!=0; i+=2 ) {
      std::string opt(atts[i]), val(atts[i+1]);
      if (opt == "version") vers = atof(atts[i+1]);
    }
    if (vers == 0) {
      std::cerr << "XML:startElement(): no version number in this file\n";
    }
    if (vers < 0.1) {
      std::cerr << "XML:startElement(): this file may not be read (bad version number)\n";
    }
    return;
  }

  if (d.read_mode == 1) decodeGuiSettings(d, stname, atts);
  if (d.read_mode == 2) decodePlanesSettings(d, stname, atts);
  if (d.read_mode == 3) decodeGroupsSettings(d, stname, atts);
}


/* XML_EndElementHandler */
void endElement(void *userData, const char *name)
{
  UserData &d = *reinterpret_cast<UserData*>(userData);
  std::string stname(name);
  d.depth--;

  if (stname == "guisetup") d.read_mode = 0;
  if (stname == "planesetup") d.read_mode = 0;
  if (stname == "groupsetup") d.read_mode = 0;

  if ((d.read_mode == 1) && ((stname == "pad") || (stname == "canvas"))) {
    if (d.skipcanvas) return;

    d.pad->Update();
    register TPaveStats *pv;
    pv = (TPaveStats*) d.pad->FindObject("stats");
    if (pv && (d.paveostat != "")) {
      int ostat;
      if (sscanf(d.paveostat.c_str(), "%d", &ostat)) pv->SetOptStat(ostat);
      pv->SetStatFormat(d.paveostatfmt.c_str());
      Double_t osx1, osy1, osx2, osy2;
      if (sscanf(d.paveostatXY.c_str(), "%lf %lf %lf %lf", &osx1, &osy1, &osx2, &osy2) == 4) {
        pv->SetX1NDC(osx1);
        pv->SetY1NDC(osy1);
        pv->SetX2NDC(osx2);
        pv->SetY2NDC(osy2);
        pv->ConvertNDCtoPad();
      }
      d.ClearPave();
      d.pad->Modified();
      if (d.psfile) d.canvas->Update();
    }
    
    TPaveText* pvtit = (TPaveText*) d.pad->FindObject("title");
    if (pvtit && (d.pavetitle != "")) {
      Double_t titx1, tity1, titx2, tity2;
      if (sscanf(d.pavetitleXY.c_str(), "%lf %lf %lf %lf", 
		 &titx1, &tity1, &titx2, &tity2) == 4) {
        pvtit->SetX1NDC(titx1);
        pvtit->SetY1NDC(tity1);
        pvtit->SetX2NDC(titx2);
        pvtit->SetY2NDC(tity2);
        pvtit->ConvertNDCtoPad();
      }
      d.ClearTitPave();
      d.pad->Modified();
      if (d.psfile) d.canvas->Update();
    }

    if (d.psfile) {
      if (stname == "pad") {
        d.psfile->On();
        d.pad->Paint();
        d.psfile->Off();
      }
      if (stname == "canvas") {
        if (d.pad == d.canvas) {
          d.psfile->On();
          d.canvas->Paint();
          d.psfile->Off();
        }
        d.canvas->Close();
        SafeDelete(d.canvas);
        d.ClearCanvas();
      }
    }
  }
}


/* XML_DefaultHandler */
void default_handler(void *userData,const XML_Char *s,int len)
{
  char buf[500];
  for (register int i=0; i<len; i++) buf[i]=s[i];
  buf[len]=0;
  return;
}

/* XML_CharacterDataHandler */
void characterdata_handler(void *userData,const XML_Char *s,int len)
{
  char buf[500];
  for (register int i=0; i<len; i++) buf[i]=s[i];
  buf[len]=0;
  return;
}

/* XML_ProcessingInstructionHandler */
void instruction_handler(void *userData,const XML_Char *target, const XML_Char *data)
{
  std::cerr << "instruction_handler  target= " << target << " data= " << data <<std::endl;
  return;
}

/* XML_CommentHandler */
void comment_handler(void *userData, const XML_Char *data)
{
  std::cerr << "Comment: " << data <<std::endl;
}
}  // namespace (private namespace for this file)



Bool_t VariousSettings::LoadConfigFile(std::string filename, int type_read, Bool_t save_name)
{
  XML_Parser parser = XML_ParserCreate(NULL);
  UserData d;

  const unsigned buf_size=1000000;
  char *buf = new char[buf_size];
  XML_SetUserData                       (parser, &d);
  XML_SetElementHandler                 (parser, startElement, endElement);
//   XML_SetCharacterDataHandler           (parser, characterdata_handler);
  XML_SetProcessingInstructionHandler   (parser, instruction_handler);
  XML_SetDefaultHandler                 (parser, default_handler);
  XML_SetCommentHandler                 (parser, comment_handler);

  ifstream in(filename.c_str());
  if (!in) {
#ifndef __CINT__
    std::cerr << "VariousSettings::LoadConfigFile: Can not read GUI settings file " <<filename<<std::endl;;
#endif
    return false;
  }

  d.monitor = fMonitor;
  d.read_mode = 0;
  d.type_read = type_read;
  d.psfile = 0;
  try
  {
    int done=0;
    if (fMonitor->GetThreadState()) if (thr_flag) TThread::Lock();
    do
    {
      in.read(buf,buf_size);
      if( ((unsigned int)in.gcount())==buf_size )
#ifndef __CINT__
        throw CS::Exception("MainFrame::LoadConfigFile(): too short buffer.");
#endif

      done = ((unsigned int)in.gcount()) < buf_size;
//  std::cerr << "appel XML_Parse: in.gcount()=" << in.gcount() << " done=" <<done <<std::endl;
      if( !XML_Parse(parser, buf, in.gcount(), done) )
#ifndef __CINT__
        CS::Exception("MainFrame::LoadConfigFile(): Error: %s  at line %d (buf: %s)\n",
		      XML_ErrorString(XML_GetErrorCode(parser)),
		      XML_GetCurrentLineNumber(parser), buf).Print();
#endif
    } while (!done);
    if (fMonitor->GetThreadState()) if (thr_flag) TThread::UnLock();
  }
  catch(...)
  {
    if (fMonitor->GetThreadState()) if (thr_flag) TThread::UnLock();
    delete [] buf;
    throw;
  }

  delete [] buf;
  in.close();
  XML_ParserFree(parser);

  if (save_name) {
    fConfigFileName = filename;
    AddNameToList(fConfigFileName);
  }
  return true;
}


Bool_t VariousSettings::LoadConfigSettings(std::string filename, int type_read, Bool_t save_name)
{
  if (filename == "") {
    std::cerr << "VariousSettings::LoadConfigSettings: no file name given\n";
    return kFALSE;
  }
  struct stat info;
  if( 0!=stat(filename.c_str(),&info) ) {
    std::cerr << "VariousSettings::LoadConfigSettings: Can't stat file "<<filename<<std::endl;
    return kFALSE;
  }
  short dirflag;
  if (S_ISDIR(info.st_mode)) dirflag = 1;
  else if (S_ISREG(info.st_mode)) dirflag = 0;
  else {
    std::cerr << "VariousSettings::LoadConfigSettings: "<<filename;
    std::cerr << " not a file or a directory\n";
    return kFALSE;
  }

  if (dirflag) {
    struct dirent *filedir;
    std::string filedirname;
    DIR* dirfd;
    std::cerr << "Reading all config files in "<<filename<<" directory:\n";
    dirfd = opendir(filename.c_str());
    while ((filedir = readdir(dirfd))) {
      register char* strcfg = strstr(filedir->d_name,".cfg");
      register char* strxml = strstr(filedir->d_name,".xml");
      if (!((strcfg && (strlen(strcfg) == 4)) || (strxml && (strlen(strxml) == 4)))) continue;
      filedirname = filename + "/" + filedir->d_name;
      std::cerr << "  - reading "<<filedir->d_name<<" -> ";
      if (LoadConfigFile(filedirname, type_read, kFALSE)) std::cerr << "success\n";
        else std::cerr << "failed\n";
    }
    closedir(dirfd);
  }
  else {
    Bool_t retcfgfile = LoadConfigFile(filename, type_read, save_name);
    if (retcfgfile) std::cerr << "Config file "<<filename<<" read\n";
    else std::cerr<<"Config file "<<filename<<" can't be read\n";
    return retcfgfile;
  }

  return kFALSE;
}


Bool_t VariousSettings::GuiSettingsToPS(const char* filename, TPostScript* psfile)
{
  XML_Parser parser = XML_ParserCreate(NULL);
  UserData d;

  const unsigned buf_size=1000000;
  char *buf = new char[buf_size];
  XML_SetUserData                       (parser, &d);
  XML_SetElementHandler                 (parser, startElement, endElement);
//   XML_SetCharacterDataHandler           (parser, characterdata_handler);
  XML_SetProcessingInstructionHandler   (parser, instruction_handler);
  XML_SetDefaultHandler                 (parser, default_handler);
  XML_SetCommentHandler                 (parser, comment_handler);

  ifstream in(filename);
  if (!in) {
#ifndef __CINT__
  std::cerr << "VariousSettings::GuiSettingsToPS: Can not read GUI settings file "<<filename<<std::endl;
#endif
    return false;
  }
  if (!psfile) {
    std::cerr << "VariousSettings::GuiSettingsToPS: No TPostScript handler defined\n";
    return false;
  }

  d.monitor = fMonitor;
  d.read_mode = 0;
  d.type_read = 1;
  d.psfile = psfile;
  d.filename = filename;
  time_t temps = time(0);
  d.datestr = ctime(&temps);
  d.runnb = fMonitor->GetRunNumber();
  d.evtnb = fMonitor->GetEventNumber();
  try
  {
    int done=0; 
    if (fMonitor->GetThreadState()) if (thr_flag) TThread::Lock();
    do
    {
      in.read(buf,buf_size);
      if( ((unsigned int)in.gcount())==buf_size )
#ifndef __CINT__
        throw CS::Exception("MainFrame::GuiSettingsToPS(): too short buffer.");
#endif

      done = ((unsigned int)in.gcount()) < buf_size;
      if( !XML_Parse(parser, buf, in.gcount(), done) )
#ifndef __CINT__
        CS::Exception("MainFrame::LoadSettings(): Error: %s  at line %d (buf: %s)\n",
		      XML_ErrorString(XML_GetErrorCode(parser)),
		      XML_GetCurrentLineNumber(parser), buf).Print();
#endif
    } while (!done);
    if (fMonitor->GetThreadState()) if (thr_flag) TThread::UnLock();
  }
  catch(...)
  {
    std::cerr<<"VariousSettings::GuiSettingsToPS : unknown exception !"<<std::endl;
    if (fMonitor->GetThreadState()) if (thr_flag) TThread::UnLock();
    delete [] buf;
    throw;
  }

  delete [] buf;
  in.close();
  XML_ParserFree(parser);
  return true;
}


#ifndef __CINT_
// filter function for directory scaning of .cfg files
int localFilterDotcfg (const struct dirent *dent) {
  return ! fnmatch("*.cfg", dent->d_name, FNM_CASEFOLD);
}


// comparison function to sort directory using inodes numbers
// int localSortByInode (const void *a, const void *b)  // no more accepted by gcc
int localSortByInode (const struct dirent **d1, const struct dirent **d2) {
  // the use of d_ino element of the dirent structure is not authorized by POSIX, bug may come from here
  // but I'm too lazy to stat the files just to get back the inode so... (DN 27/6/2011)
//   const struct dirent **d1 = (const struct dirent **) a;
//   const struct dirent **d2 = (const struct dirent **) b;
  if ((*d1)->d_ino > (*d2)->d_ino) return 1;
  if ((*d1)->d_ino < (*d2)->d_ino) return -1;
  return 0;
}
#endif


Bool_t VariousSettings::GuiSettingsCreatePS(std::string filename, const TList* filelist) {
  struct stat info;
  short dirflag;
  TVirtualPad* savepad = gPad;
  TPostScript* psfile = 0;

  if (filename != "") {
    if ( 0!=stat(filename.c_str(),&info) ) {
      std::cerr << "VariousSettings::GuiSettingsCreatePS: Can't stat file "<<filename<<std::endl;
      return kFALSE;
    }

    if (S_ISDIR(info.st_mode)) dirflag = 1;
    else if (S_ISREG(info.st_mode)) dirflag = 0;
    else {
      std::cerr << "VariousSettings::GuiSettingsCreatePS: "<<filename;
      std::cerr << " not a file or a directory\n";
      return kFALSE;
    }

    std::string psname = filename;
    unsigned int possl = filename.rfind("/");  // to remove all path, the file is created in the current dir
    if (possl != std::string::npos) psname=filename.substr(possl+1);
    unsigned int pospt = psname.rfind(".");
    if (pospt != std::string::npos) psname.erase(pospt); // to remove the last .something suffix
    psname+=".ps";
    psfile = new TPostScript(psname.c_str(), 115); // Landscape orientation, color
    psfile->Range(29, 19);
    psfile->Off();

    if (dirflag) {
#ifndef __CINT_
      struct dirent **namelist;
      std::string filedirname;
      std::cerr << "Create PS file "<<psname<<" from all config files in "<<filename<<" directory...\n";
      register int nbfiles = scandir(filename.c_str(), &namelist, &localFilterDotcfg, &localSortByInode);
      if (nbfiles < 0)
        perror("VariousSettings::GuiSettingsCreatePS: scandir error");
      else {
        for (register int ii=0; ii<nbfiles; ii++) {
          register char* strcfg = strstr(namelist[ii]->d_name,".cfg");  // just in case...
          if (!(strcfg && (strlen(strcfg) == 4))) continue;
          filedirname = filename + "/" + namelist[ii]->d_name;
          std::cerr << "  ...including "<<namelist[ii]->d_name<<std::endl;
          GuiSettingsToPS(filedirname.c_str(), psfile);
          free (namelist [ii]);
        }
        free (namelist);
      }
#endif
    }
    else GuiSettingsToPS(filename.c_str(), psfile);


  } else {  // filename == ""


    if (!filelist) {
      std::cerr << "VariousSettings::GuiSettingsCreatePS: no file name given\n";
      return kFALSE;
    } else {

      register int nbfiles = filelist->GetEntries();
      if (nbfiles < 0)
        perror("VariousSettings::GuiSettingsCreatePS: problem in file list");
      else {
        std::string psname = fPwd+"/config_file_list.ps";
        std::cerr << "Create PS file "<<psname<<" from all config files in file list...\n";
        psfile = new TPostScript(psname.c_str(), 115); // Landscape orientation, color
        psfile->Range(29, 19);
        psfile->Off();
        for (register int ifile = 0; ifile < nbfiles; ifile++) {
          const char* iname = filelist->At(ifile)->GetName();
          const char* strcfg = strstr(iname,".cfg");  // just in case...
          if (!(strcfg && (strlen(strcfg) == 4))) continue;
          std::cerr << "  ...including "<<iname<<std::endl;
          GuiSettingsToPS(iname, psfile);
        }
      }

    }
  }

  gPad = 0;
  psfile->On();
  psfile->Close();
  SafeDelete(psfile);
  if(savepad) savepad->cd();
  return kTRUE;
}


Bool_t VariousSettings::GuiSettingsCreateOrderlyShiftPlots() {
  //std::string filename = std::string(ORDERLY_SHIFT_PLOTS_DIR);
  std::string filename = std::string(SHIFT_DIR)+"/" + 
    std::string(ORDERLY_SHIFT_PLOTS_DIR);
  struct stat info;
  if( 0!=stat(filename.c_str(),&info) ) {
    std::cerr << "VariousSettings::GuiSettingsCreateOrderlyShiftPlots: Can't stat file "<<filename<<std::endl;
    return kFALSE;
  }

  short dirflag = 1;

  std::string psname = std::string(SHIFT_DIR)+"/"+ORDERLY_SHIFT_PLOTS_FNAME;
  cout << "PS file name: " << psname << endl;
  TVirtualPad* savepad = gPad;
  TPostScript* psfile = new TPostScript(psname.c_str(), 115); // Landscape orientation, color
  psfile->Range(29, 19);
  psfile->Off();

  if (dirflag) {
#ifndef __CINT_
    struct dirent **namelist;
    std::string filedirname;
    std::cerr << "Create PS file "<<psname<<" from all config files in "<<filename<<" directory...\n";
    register int nbfiles = scandir(filename.c_str(), &namelist, &localFilterDotcfg, &localSortByInode);
    if (nbfiles < 0)
      perror("VariousSettings::GuiSettingsCreatePS: scandir error");
    else {
      for (register int ii=0; ii<nbfiles; ii++) {
        register char* strcfg = strstr(namelist[ii]->d_name,".cfg");  // just in case...
        if (!(strcfg && (strlen(strcfg) == 4))) continue;
        filedirname = filename + "/" + namelist[ii]->d_name;
        std::cerr << "  ...including "<<namelist[ii]->d_name<<std::endl;
        GuiSettingsToPS(filedirname.c_str(), psfile);
        free (namelist [ii]);
      }
      free (namelist);
    }
#endif
  }
  else GuiSettingsToPS(filename.c_str(), psfile);

  gPad = 0;
  psfile->On();
  psfile->Close();
  SafeDelete(psfile);
  if(savepad) savepad->cd();
  return kTRUE;
}



// ---------------------------------------------------------------------

// personal parameters .cooolrc file handling


bool VariousSettings::ReadPersoParam() {
  std::string file=getenv("HOME");
  file+="/"; std::string file2 = file;
  file+= fPersoFile;
  register bool fflag;
  fflag = ReadPersoParam(file.c_str());
  if (fflag == false) {
    file2 += ".CompMon"; // old name for .cooolrc
    fflag = ReadPersoParam(file2.c_str());
  }
  return fflag;
}


bool VariousSettings::ReadPersoParam(const char * filename) {

  bool filelist_mode = false;

  ifstream in(filename);

  if(! in) {
    std::cerr<<"VariousSettings::ReadPersoParam : Cannot open personal parameters file "<<filename<<std::endl;
    return false;
  }
  while(1) {
    std::string tname;
    int flag=0;
    in >> tname;
    if(in.fail()) break;
    if (tname == "Config_files_list:") { filelist_mode = true; continue;}
    if (tname == "End_config_files_list") { filelist_mode = false; continue;}
    if (filelist_mode) {
      fConfigFileNameList.push_back(tname);
      continue;
    }
    else {
      in >> flag;
      std::map<std::string,bool>::iterator im = fMonitor->GetDetInTree().find(tname);
      if(im!=fMonitor->GetDetInTree().end())
	im->second = flag;
    }
  }
  in.close();

  return true;
}


bool VariousSettings::WritePersoParam(const char * filename) {

  std::string file=getenv("HOME");
  file+="/";
  file+=filename;

  ofstream out(file.c_str(), std::ios::out | std::ios::trunc);

  if(! out) {
    std::cerr<<"Monitor::WritePersoParam : Cannot write personal parameters file"<<std::endl;
    return false;
  }

  typedef std::map<std::string,bool>::iterator IM;
  for(IM i=fMonitor->GetDetInTree().begin();i!=fMonitor->GetDetInTree().end();i++)
    out<<i->first<<"\t"<<static_cast<int>(i->second)<<std::endl; 

  out << "\n\nConfig_files_list:\n\n";
  typedef std::vector<std::string>::iterator LS;
  for(LS j = fConfigFileNameList.begin(); j != fConfigFileNameList.end(); j++) {
    out << (*j) <<std::endl;
  }
  out << "\nEnd_config_files_list\n";
  out.close();

  return true;
}


void VariousSettings::AddNameToList(std::string name) {

  typedef std::vector<std::string>::iterator LS;
  LS i = fConfigFileNameList.begin();
  int nb_elem = 0;
  while (i != fConfigFileNameList.end()) {
    bool fglast = false;
    if ((*i) != name) { i++; nb_elem++; continue; }
    LS i2 = i;
    i++;
    if (i == fConfigFileNameList.end()) {
      fglast = true;
    }
    fConfigFileNameList.erase(i2);
    if (fglast) break;
  }
  if (nb_elem >= fMaxItemInList) fConfigFileNameList.pop_back();
  fConfigFileNameList.insert(fConfigFileNameList.begin(), name);
}









