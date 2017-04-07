/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionTest.cc,v $
   $Date: 2011/01/31 20:35:45 $
   $Revision: 1.47 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     Denis     Murashev  ( Denis.Mourachev@cern.ch, Murashev@sirius.ihep.su )

   Copyright(C): 1999-2001  V.Kolosov, A.Zvyagin, D.Murashev

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

// --- Standard C/C++ library ---
#include <cassert>
#include <cmath>
#include <iostream>

// --- Histo library ---
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TMath.h"

// --- Internal files ---
#include "Calorimeter.h"

#include "CalorimeterHist.h"
#include "CellDataRaw.h"
#include "CellType.h"
#include "OneParticleResponse.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::ProcessOnDuty(void)
{
  FillInternalCorrelations();
  FillExternalCorrelations();
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

CalorimeterHist::CalorimeterHist(Calorimeter *calo) :
calorimeter_(calo)
{
  raw_histos_booked=false;
  led_histos_booked=false;
  ped_histos_booked=false;
  test_histos_booked=false;
  root_dir        = NULL;
  stat_dir        = NULL;
  test_dir        = NULL;
  cells_dir       = NULL;
  cells_time_dir  = NULL;
  cells_led_dir   = NULL;
  cells_ped_dir   = NULL;
  noise_dir       = NULL;
  cells_calib_dir = NULL;

  h1_ngen = NULL;
  h1_ngen_Out = NULL;
  h1_nrec = NULL;
  h2_ngen_nrec = NULL;

  h1_Xgen_Out = NULL;
  h1_Ygen_Out = NULL;
  h1_Egen_Out = NULL;
  h2_XgYg_Out = NULL;

  h1_Xgen = NULL;
  h1_Ygen = NULL;
  h1_Rgen = NULL;
  h1_Egen = NULL;
  h2_XgYg = NULL;

  h1_Xgen_Norec = NULL;
  h1_Ygen_Norec = NULL;
  h1_Egen_Norec = NULL;
  h2_XgYg_Norec = NULL;

  h1_Xrec = NULL;
  h1_Yrec = NULL;
  h1_Rrec = NULL;
  h1_Erec = NULL;
  h1_Sumrec = NULL;
  h2_XrYr = NULL;

  h1_Xdelta = NULL;
  h1_Ydelta = NULL;
  h2_EgEr = NULL;
  h2_TgTr = NULL;
  h2_XgXr = NULL;
  h2_YgYr = NULL;
  h2_XYdelta = NULL;
  h1_Edelta = NULL;
  h1_EdeltaN = NULL;

  h1_XdeltaHadron = NULL;
  h1_YdeltaHadron = NULL;
  h1_XdeltaHadron40 = NULL;
  h1_YdeltaHadron40 = NULL;
  h2_EgErHadron = NULL;
  h1_XdeltaMuon = NULL;
  h1_YdeltaMuon = NULL;
  h2_EgErMuon = NULL;

  h2_EdeltaEg = NULL;
  h2_EdeltaN_Eg = NULL;
  h2_XdeltaXg = NULL;
  h2_YdeltaYg = NULL;
  h2_XdeltaYdelta = NULL;

 // For Cells Test
  h2_XgYgCellTest = NULL;
  h2_XrYrCellTest = NULL;
  h2_XrYrMissedCell = NULL;
  h2_XgXrCellTest = NULL;
  h2_YgYrCellTest = NULL;

 // For Near Gamma Test
  h2_Dist2 = NULL;
  h2_Dist2emin = NULL;
  h2_Dist2Cell = NULL;

 // For total hited cells and total energy at a time
  h1_ESumADC = NULL;
  h1_ETotal = NULL;
  h2_NHitedEtotal = NULL;
  h2_XhitYhit = NULL;
  h2_XhitYhitE = NULL;;

  h1_CalibFactorBase = NULL;
  h1_CalibFactorLED = NULL;
  p1_CalibFactorEdep = NULL;
  p1_CalibFactorTis = NULL;

  for(int i=0; i<(int)calorimeter_->NCells(); i++) h1_ErCell.push_back(NULL);
  for(int i=0; i<(int)calorimeter_->NCells(); i++) h1_ADCmCell.push_back(NULL);

  h1_Time = NULL;
  h2_TimeE = NULL;

  for(int i=0; i<(int)calorimeter_->NCells(); i++) h1_TimeCell.push_back(NULL);

  for(int i=0; i<(int)calorimeter_->NCells(); i++) h1_RAWCell.push_back(NULL);
  for(int i=0; i<(int)calorimeter_->NCells(); i++) h1_LEDCell.push_back(NULL);
  for(int i=0; i<(int)calorimeter_->NCells(); i++) h1_PEDCell.push_back(NULL);

  for(int i=0; i<(int)calorimeter_->NCells(); i++) h1_CalibCell.push_back(NULL);
  for(int i=0; i<(int)calorimeter_->NCells(); i++) h2_CalibInSpill.push_back(NULL);

  for(int i=0; i<(int)calorimeter_->NCells(); i++)
  {
    h1_NOISECell.push_back(NULL);
    h1_NOISEGamCell.push_back(NULL);
  }

  h2_XYAll = NULL;
  for(int i=0; i<(int)calorimeter_->NCells(); i++) h2_XdeltaYdeltaCell.push_back(NULL);
 // Statistic Histo
  h1_calib_new = NULL;
  h1_calib_old = NULL;
  h1_calib_prim = NULL;
  h1_calib_new_over_calib_prim = NULL;
  h1_calib_new_over_calib_old = NULL;
  h1_calib_old_over_calib_prim = NULL;
  h1_led_new = NULL;
  h1_led_old = NULL;
  h1_ped_new = NULL;
  h1_ped_old = NULL;
  h1_time_new = NULL;
  h1_time_old = NULL;
  h1_ecut_low_old = NULL;
  h1_ecut_high_old = NULL;
  h1_ecut_low_new = NULL;
  h1_ecut_high_new = NULL;
  h1_ecut_sparse_mode = NULL;
  h1_ecut_min = NULL;
}

////////////////////////////////////////////////////////////////////////////////

CalorimeterHist::~CalorimeterHist(void)
{
  std::cerr << " We dont expect call to CalorimeterHist destructor !!! " << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void CalorimeterHist::Reset(void)
{
  if(h1_ngen != NULL) h1_ngen->Reset();
  if(h1_ngen_Out != NULL) h1_ngen_Out->Reset();
  if(h2_ngen_nrec != NULL) h2_ngen_nrec->Reset();

  if(h1_Xgen_Out != NULL) h1_Xgen_Out->Reset();
  if(h1_Ygen_Out != NULL) h1_Ygen_Out->Reset();
  if(h1_Egen_Out != NULL) h1_Egen_Out->Reset();
  if(h2_XgYg_Out != NULL) h2_XgYg_Out->Reset();

  if(h1_Xgen != NULL) h1_Xgen->Reset();
  if(h1_Ygen != NULL) h1_Ygen->Reset();
  if(h1_Rgen != NULL) h1_Rgen->Reset();
  if(h1_Egen != NULL) h1_Egen->Reset();
  if(h2_XgYg != NULL) h2_XgYg->Reset();

  if(h1_Xgen_Norec != NULL) h1_Xgen_Norec->Reset();
  if(h1_Ygen_Norec != NULL) h1_Ygen_Norec->Reset();
  if(h1_Egen_Norec != NULL) h1_Egen_Norec->Reset();
  if(h2_XgYg_Norec != NULL) h2_XgYg_Norec->Reset();

  if(h1_nrec != NULL) h1_nrec->Reset();
  if(h1_Xrec != NULL) h1_Xrec->Reset();
  if(h1_Yrec != NULL) h1_Yrec->Reset();
  if(h1_Rrec != NULL) h1_Rrec->Reset();
  if(h1_Erec != NULL) h1_Erec->Reset();
  if(h1_Sumrec != NULL) h1_Sumrec->Reset();
  if(h2_XrYr != NULL) h2_XrYr->Reset();

  if(h1_Xdelta != NULL) h1_Xdelta->Reset();
  if(h1_Ydelta != NULL) h1_Ydelta->Reset();
  if(h2_EgEr != NULL) h2_EgEr->Reset();
  if(h2_TgTr != NULL) h2_TgTr->Reset();
  if(h2_XgXr != NULL) h2_XgXr->Reset();
  if(h2_YgYr != NULL) h2_YgYr->Reset();
  if(h2_XYdelta != NULL) h2_XYdelta->Reset();
  if(h1_Edelta != NULL) h1_Edelta->Reset();
  if(h1_EdeltaN != NULL) h1_EdeltaN->Reset();

  if(h1_XdeltaHadron != NULL) h1_XdeltaHadron->Reset();
  if(h1_YdeltaHadron != NULL) h1_YdeltaHadron->Reset();
  if(h1_XdeltaHadron40 != NULL) h1_XdeltaHadron40->Reset();
  if(h1_YdeltaHadron40 != NULL) h1_YdeltaHadron40->Reset();
  if(h2_EgErHadron != NULL) h2_EgErHadron->Reset();
  if(h1_XdeltaMuon != NULL) h1_XdeltaMuon->Reset();
  if(h1_YdeltaMuon != NULL) h1_YdeltaMuon->Reset();
  if(h2_EgErMuon != NULL) h2_EgErMuon->Reset();

  if(h2_EdeltaEg != NULL) h2_EdeltaEg->Reset();
  if(h2_EdeltaN_Eg != NULL) h2_EdeltaN_Eg->Reset();
  if(h2_XdeltaXg != NULL) h2_XdeltaXg->Reset();
  if(h2_YdeltaYg != NULL) h2_YdeltaYg->Reset();
  if(h2_XdeltaYdelta != NULL) h2_XdeltaYdelta->Reset();

 // For Cells Test
  if(h2_XgYgCellTest != NULL) h2_XgYgCellTest->Reset();
  if(h2_XrYrCellTest != NULL) h2_XrYrCellTest->Reset();
  if(h2_XrYrMissedCell != NULL) h2_XrYrMissedCell->Reset();
  if(h2_XgXrCellTest != NULL) h2_XgXrCellTest->Reset();
  if(h2_YgYrCellTest != NULL) h2_YgYrCellTest->Reset();
 // For Near Gamma Test
  if(h2_Dist2 != NULL) h2_Dist2->Reset();
  if(h2_Dist2emin != NULL) h2_Dist2emin->Reset();
  if(h2_Dist2Cell != NULL) h2_Dist2Cell->Reset();

 // For total hited cells and total energy at a time
  if(h1_ESumADC != NULL) h1_ESumADC->Reset();
  if(h1_ETotal != NULL) h1_ETotal->Reset();
  if(h2_NHitedEtotal != NULL) h2_NHitedEtotal->Reset();
  if(h2_XhitYhit != NULL) h2_XhitYhit->Reset();
  if(h2_XhitYhitE != NULL) h2_XhitYhitE->Reset();

  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h1_ErCell[i] != NULL) h1_ErCell[i]->Reset();
  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h1_ADCmCell[i] != NULL) h1_ADCmCell[i]->Reset();

  if(h1_Time != NULL) h1_Time->Reset();
  if(h2_TimeE != NULL) h2_TimeE->Reset();

  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h1_TimeCell[i] != NULL) h1_TimeCell[i]->Reset();

  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h1_RAWCell[i] != NULL) h1_RAWCell[i]->Reset();
  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h1_LEDCell[i] != NULL) h1_LEDCell[i]->Reset();
  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h1_PEDCell[i] != NULL) h1_PEDCell[i]->Reset();

  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h1_CalibCell[i] != NULL) h1_CalibCell[i]->Reset();

  if(h2_XYAll != NULL) h2_XYAll->Reset();

  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h2_XdeltaYdeltaCell[i] != NULL) h2_XdeltaYdeltaCell[i]->Reset();

//  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h1_NOISECell[i] != NULL)  h1_NOISECell[i]->Reset();
//  for(int i=0; i<(int)calorimeter_->NCells(); i++) if(h1_NOISEGamCell[i] != NULL) h1_NOISEGamCell[i]->Reset();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ReconstructionTest(void)
{
  bool debug_ReconstructionTest=false;
//  bool debug_ReconstructionTest=true;
  bool print_missed_cell_warning = false;

  debug_ReconstructionTest=options.debug_reconstruction_test;

  if(debug_ReconstructionTest) cout << " ReconstructionTest " << GetName() << endl;
  bool ok=false;
  if( !options.fill_histos )
    return;
  bool need_cells_timing_histo= options.enable_time_calibration;

  if( gDirectory==NULL )
  {
    cerr << "Calorimeter::ReconstructionTest():  ROOT file was not opend.\n"
         << "                                    I do not want to work.\n";
    return;
  }

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.

// create test histos directory sructure
  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  if( h.root_dir==NULL )
  {
//     TDirectory *dir_check;
    char dir_name[111];
  if(debug_ReconstructionTest) cout << " ReconstructionTest gDirectory->cd(/) " << endl;
    ok = gDirectory->cd("/");
  if(debug_ReconstructionTest) cout << " ReconstructionTest gDirectory->cd(/) OK " << ok << endl;
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
  }
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.root_dir->cd() " << endl;
  ok = h.root_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.root_dir->cd() OK " << ok <<  endl;
  assert(ok);

  float cellx=40;
  float celly=40;
  if( NCells() > 0 )
  {
//    cellx = cells[0].GetCellType().GetSizeX();
//    celly = cells[0].GetCellType().GetSizeY();
    cellx = GetMinCellSizeX();
    celly = GetMinCellSizeY();
  }
  else
  {
    cerr << " Calorimeter::ReconstructionTest ERROR: Calorimeter " << GetName() << " is empty !!! " << endl;
    return;
  }

// TODO should be optional

  double e_cut = 2.;
  bool only_reco_hist = options.fill_only_reco_histos_in_reco_test;
  if(!h.test_histos_booked)
  {
    char hist_name[132];
    char dir_name[111];
// #warning "WARNING not the best place for booking and filling "
// storing some statistic about cells
    sprintf(dir_name,"%s_stat",GetName().c_str());
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.stat_dir->cd() " << endl;
    h.stat_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok=h.stat_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.stat_dir->cd() OK " << ok <<endl;
    assert(ok);

    float e_peak = options.calibration_energy;
    if( e_peak <= 0. ) e_peak = 10.;

    double maxdev_cf_joujou = 0.20;
    double ecut_hist_max = 20.;

    sprintf(hist_name,"%s:: Peak Position of %7.2f GeV in ADCs CALIB OLD ",GetName().c_str(),e_peak);
    h.h1_calib_old = myROOT_utils::TH1D_checked("CalibOLD",hist_name,4096,0,4096);

    sprintf(hist_name,"%s:: CALIBRATION (OLD-PRIM)/PRIM ",GetName().c_str());
    h.h1_calib_old_over_calib_prim = myROOT_utils::TH1D_checked("CalibOLDoverPRIM",hist_name,800,-1.,3.);

    sprintf(hist_name,"%s:: Peak Position of %7.2f GeV in ADCs CALIB NEW ",GetName().c_str(),e_peak);
    h.h1_calib_new = myROOT_utils::TH1D_checked("CalibNEW",hist_name,4096,0,4096);

    sprintf(hist_name,"%s:: CALIBRATION (NEW-PRIM)/PRIM ",GetName().c_str());
    h.h1_calib_new_over_calib_prim = myROOT_utils::TH1D_checked("CalibNEWoverPRIM",hist_name,800,-1.,3.);

    sprintf(hist_name,"%s:: CALIBRATION (NEW-OLD)/OLD ",GetName().c_str());
    h.h1_calib_new_over_calib_old = myROOT_utils::TH1D_checked("CalibNEWoverOLD",hist_name,800,-1.,3.);

    sprintf(hist_name,"%s:: CALIBRATION (MONITOR-OLD)/OLD ",GetName().c_str());
    h.h1_calib_monitor_over_calib_old = myROOT_utils::TH1D_checked("CalibMONITORoverOLD",hist_name,800,-1.,3.);

    sprintf(hist_name,"%s:: LED (MONITOR-OLD)/OLD ",GetName().c_str());
    h.h1_led_monitor_over_led_old = myROOT_utils::TH1D_checked("LEDMONITORoverOLD",hist_name,800,-1.,3.);

    double xmi = GetXmin() , xma = GetXmax() , ymi = GetYmin() , yma = GetYmax();
    sprintf(hist_name,"%s XY for OUT of CALIBRATION (MONITOR-OLD)/OLD ",GetName().c_str());
    h.h2_XYout_calib_monitor_over_old = myROOT_utils::TH2D_checked("CalibXYMONITORoverOLD",hist_name,100,xmi,xma,100,ymi,yma);

    sprintf(hist_name,"%s:: Peak Position of %7.2f GeV in ADCs CALIB PRIM ",GetName().c_str(),e_peak);
    h.h1_calib_prim = myROOT_utils::TH1D_checked("CalibPRIM",hist_name,4096,0,4096);

    sprintf(hist_name,"%s:: LED OLD ",GetName().c_str());
    h.h1_led_old = myROOT_utils::TH1D_checked("LEDOLD",hist_name,4096,0,4096);

    sprintf(hist_name,"%s:: LED NEW ",GetName().c_str());
    h.h1_led_new = myROOT_utils::TH1D_checked("LEDNEW",hist_name,4096,0,4096);

    sprintf(hist_name,"%s:: PED OLD ",GetName().c_str());
    h.h1_ped_old = myROOT_utils::TH1D_checked("PEDOLD",hist_name,4096,0,4096);

    sprintf(hist_name,"%s:: PED NEW ",GetName().c_str());
    h.h1_ped_new = myROOT_utils::TH1D_checked("PEDNEW",hist_name,4096,0,4096);

    sprintf(hist_name,"%s:: TIME OLD ",GetName().c_str());
    h.h1_time_old = myROOT_utils::TH1D_checked("TIMEOLD",hist_name,1000,-500.,500.);

    sprintf(hist_name,"%s:: TIME NEW ",GetName().c_str());
    h.h1_time_new = myROOT_utils::TH1D_checked("TIMENEW",hist_name,1000,-500.,500.);

    sprintf(hist_name,"%s:: Ecut Cell Low OLD ",GetName().c_str());
    h.h1_ecut_low_old = myROOT_utils::TH1D_checked("ECUT_CELL_LOW_OLD",hist_name,1000,0,ecut_hist_max);

    sprintf(hist_name,"%s:: Ecut Gamma High OLD ",GetName().c_str());
    h.h1_ecut_high_old = myROOT_utils::TH1D_checked("ECUT_GAMMA_HIGH_OLD",hist_name,1000,0,ecut_hist_max);

    sprintf(hist_name,"%s:: Ecut Cell Low NEW ",GetName().c_str());
    h.h1_ecut_low_new = myROOT_utils::TH1D_checked("ECUT_CELL_LOW_ NEW",hist_name,1000,0,ecut_hist_max);

    sprintf(hist_name,"%s:: Ecut Gamma High NEW ",GetName().c_str());
    h.h1_ecut_high_new = myROOT_utils::TH1D_checked("ECUT_GAMMA_HIGH_ NEW",hist_name,1000,0,ecut_hist_max);

    sprintf(hist_name,"%s:: FE Ecut (sparse mode)  ",GetName().c_str());
    h.h1_ecut_sparse_mode = myROOT_utils::TH1D_checked("ECUT_SPARSE_MODE",hist_name,1000,0,ecut_hist_max);

    sprintf(hist_name,"%s:: Ecut min(sparse mode?)  ",GetName().c_str());
    h.h1_ecut_min = myROOT_utils::TH1D_checked("ECUT_MIN",hist_name,1000,0,ecut_hist_max);

    bool debug_status = debug_ReconstructionTest;
//     debug_ReconstructionTest = true;
    if(debug_ReconstructionTest) cout << " Start BOOKING subset histo " << GetName() << " subsets " << subsets_list2store_calib_histo_.size() << endl;
    for (size_t isub = 0; isub <subsets_list2store_calib_histo_.size(); isub++)
    {
      size_t issub = subsets_list2store_calib_histo_[isub];
      char name[132];
      sprintf(name,"CalibOLD_%s",GetSubSets()[issub].GetName().c_str());
      sprintf(hist_name,"%s::SubSet::%s Peak Position of %7.2f GeV in ADCs CALIB OLD ",GetName().c_str(), GetSubSets()[issub].GetName().c_str(), e_peak);
      TH1D *h1 = myROOT_utils::TH1D_checked(name,hist_name,4096,0,4096);
      if(debug_ReconstructionTest) cout << " BOOK " << hist_name << endl;
      h.h1_calib_old_subset.push_back(h1);

      sprintf(name,"CalibMONITORoverOLD_%s",GetSubSets()[issub].GetName().c_str());
      sprintf(hist_name,"%s::SubSet::%s CALIBRATION (MONITOR-OLD)/OLD ",GetName().c_str(), GetSubSets()[issub].GetName().c_str());
      h1 = myROOT_utils::TH1D_checked(name,hist_name,800,-1.,3.);
      h.h1_calib_monitor_over_old_subset.push_back(h1);

      double xmis = GetSubSets()[issub].GetXmin() , xmas = GetSubSets()[issub].GetXmax() ,
          ymis = GetSubSets()[issub].GetYmin() , ymas = GetSubSets()[issub].GetYmax();

      sprintf(name,"CalibXYMONITORoverOLD_%s",GetSubSets()[issub].GetName().c_str());
      sprintf(hist_name,"%s::SubSet::%s XY for OUT of CALIBRATION (MONITOR-OLD)/OLD ",GetName().c_str(), GetSubSets()[issub].GetName().c_str());
      TH2D *h2 = myROOT_utils::TH2D_checked(name,hist_name,100,xmis,xmas,100,ymis,ymas);
      h.h2_XYout_calib_monitor_over_old_subset.push_back(h2);
    }

    if(debug_ReconstructionTest) cout << " Start FILLING subset histo " << endl;

    for (size_t icell = 0; icell < NCells(); icell++)
    {
      double xcell = cells[icell].GetX();
      double ycell = cells[icell].GetY();
      double cf_old = cells_info[CALIB][OLD][icell].GetMean();
      double cf_prim = cells_info[CALIB][PRIM][icell].GetMean();
      double cf_monitor = cells_info[CALIB][MONITOR][icell].GetMean();
      h.h1_calib_old->Fill(  e_peak/cf_old);
      h.h1_calib_prim->Fill( e_peak/cf_prim);
      if(cf_prim > 0. && cf_old != cf_prim ) h.h1_calib_old_over_calib_prim->Fill( cf_old/cf_prim -1. );
      if( cf_old == cf_prim ) h.h1_calib_old_over_calib_prim->Fill(-2.);
      if(cf_old > 0.)
      {
        double frmcf =  cf_monitor/cf_old -1.;
        h.h1_calib_monitor_over_calib_old->Fill( frmcf );
//        if(  frmcf - 0.14 < - 0.2 || frmcf - 0.14 > 0.2 ) cout <<" Out CF of range (-0.14) " << frmcf-0.14 << GetCellName(icell) << endl;
        if(  (frmcf - 0.14) < - maxdev_cf_joujou || (frmcf - 0.14) > maxdev_cf_joujou )
        {
          h.h2_XYout_calib_monitor_over_old->Fill( xcell, ycell );
        }
      }

      h.h1_time_old->Fill( cells_info[TIME][OLD][icell].GetMean() );

      h.h1_led_old->Fill(cells_info[LED][OLD][icell].GetMean());
      h.h1_ped_old->Fill(cells_info[PED][OLD][icell].GetMean());
      h.h1_ecut_low_old->Fill(energy_cut_bad_cells_old[icell]);
      h.h1_ecut_high_old->Fill(energy_gamma_cut_bad_cells_old[icell]);
      h.h1_ecut_sparse_mode->Fill(energy_cut_sparse_mode[icell]);
      h.h1_ecut_min->Fill( GetMinEnergyCutInCell(icell) );

      if( subsets_list2store_calib_histo_.size() !=  h.h1_calib_old_subset.size() )
      {
         cerr << " ne kleit " << subsets_list2store_calib_histo_.size() << "  hsz "<<  h.h1_calib_old_subset.size() << endl;
	 exit(1);
      }
      for (size_t isub = 0; isub <subsets_list2store_calib_histo_.size(); isub++)
      {
        size_t issub = subsets_list2store_calib_histo_[isub];
        if( GetSubSets()[issub].IsMyCell(icell) )
	{
	  if( h.h1_calib_old_subset[isub] == NULL )
	  {
            cerr << " ne kleit H=NULL " << isub << "  issub "<<  issub << endl;
	    exit(1);
	  }
	  h.h1_calib_old_subset[isub]->Fill(  e_peak/cf_old);
          if(cf_old > 0.)
          {
            double frmcf =  cf_monitor/cf_old -1.;
            h.h1_calib_monitor_over_old_subset[isub]->Fill( frmcf );
            if(  (frmcf - 0.14) < - maxdev_cf_joujou || (frmcf - 0.14) > maxdev_cf_joujou )
            {
              h.h2_XYout_calib_monitor_over_old_subset[isub]->Fill( xcell, ycell );
            }
          }
        }
      }

    }

    debug_ReconstructionTest = debug_status;

// end storing some statistic about cells
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.root_dir->cd() " << endl;
    ok=h.root_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.root_dir->cd() OK " << ok << endl;
    assert(ok);

    sprintf(dir_name,"%s_test",GetName().c_str());
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.test_dir->cd() " << endl;
    h.test_dir=myROOT_utils::TDirectory_checked(dir_name);
    ok=h.test_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.test_dir->cd() OK " << ok << endl;
    assert(ok);

// Create subdirectory  for the cells
    sprintf(dir_name,"%s_Cells",GetName().c_str());
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.cells_dir->cd() " << endl;
    h.cells_dir=myROOT_utils::TDirectory_checked(dir_name);
    ok=h.cells_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.cells_dir->cd() OK " << ok << endl;
    assert(ok);

// Timing in cells
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.root_dir->cd() " << endl;
    ok=h.root_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.root_dir->cd() OK " << ok << endl;
    assert(ok);

// Create subdirectory  for the cells
    sprintf(dir_name,"%s_CellsTime",GetName().c_str());
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.cells_time_dir->cd() " << endl;
    h.cells_time_dir=myROOT_utils::TDirectory_checked(dir_name);
    ok=h.cells_time_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.cells_time_dir->cd() OK " << ok << endl;
    assert(ok);

// Create subdirectory  for the XY in cells
    ok=h.root_dir->cd();
    assert(ok);
    sprintf(dir_name,"%s_CellsXY",GetName().c_str());
    h.cells_XY_dir=myROOT_utils::TDirectory_checked(dir_name);

  if(debug_ReconstructionTest) cout << " ReconstructionTest h.test_dir->cd() OK "  << endl;
    ok=h.test_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.test_dir->cd() OK " << ok << endl;
    assert(ok);


    int   n_gam_max = 100;
    int   n_xbin1 = 200, n_ybin1 = 200, n_ebin1 = options.hist_energy_nbins_1d;
    int   n_xbin2 =  50, n_ybin2 =  50, n_ebin2 = options.hist_energy_nbins_2d;
    float xmin = (*this).GetXmin() , xmax = (*this).GetXmax() ,
          ymin = (*this).GetYmin() , ymax = (*this).GetYmax() ,
          xmingate=xmin-(xmax-xmin)/4. , xmaxgate=xmax+(xmax-xmin)/4.,
          ymingate=ymin-(ymax-ymin)/4. , ymaxgate=ymax+(ymax-ymin)/4.,
          emin = 0, emax = options.hist_energy_max;
    int cellb = 50;
//    Timing information
    int   n_time_bins = options.hist_time_nbins_1d;
    int   n_time2 = options.hist_time_nbins_2d;
//     float time_min = -50 , time_max = 50;
    float time_min = options.hist_time_min , time_max = options.hist_time_max;
// Testing and setting the good emax
//    if (hint_particles.size() > 0) {
//      emax = 1;
//      for( vector<CalorimeterParticle>::iterator it=hint_particles.begin(); it!=hint_particles.end(); it++ )
//	if ((float)it->GetE() > emax) emax = (float)it->GetE();
//      if (emax > 25) emax *= 2;
//      else (emax *= 3);
//    }
// End of emax
    sprintf(hist_name,"%s:: Number of reconstructed particles",GetName().c_str());
    h.h1_nrec = myROOT_utils::TH1D_checked("Nrec",hist_name,n_gam_max,0,n_gam_max);
    if( !only_reco_hist )
    {
      sprintf(hist_name,"%s:: Number of generated particles Missed ",GetName().c_str());
      h.h1_ngen_Out = myROOT_utils::TH1D_checked("NgenOut",hist_name,n_gam_max,0,n_gam_max);
      sprintf(hist_name,"%s:: Number of generated particles",GetName().c_str());
      h.h1_ngen = myROOT_utils::TH1D_checked("Ngen",hist_name,n_gam_max,0,n_gam_max);
      sprintf(hist_name,"%s:: Ngenerated vs Nreconstructed",GetName().c_str());
      h.h2_ngen_nrec = myROOT_utils::TH2D_checked("Ngen_Nrec",hist_name,n_gam_max,0,n_gam_max,n_gam_max,0,n_gam_max);

//  Histgramms for generated particles missed calorimeter
      sprintf(hist_name,"%s:: X distribution of generated particles Missed",GetName().c_str());
      h.h1_Xgen_Out = myROOT_utils::TH1D_checked("XgenOut",hist_name,n_xbin1,xmingate,xmaxgate);
      sprintf(hist_name,"%s:: Y distribution of generated particles Missed",GetName().c_str());
      h.h1_Ygen_Out = myROOT_utils::TH1D_checked("YgenOut",hist_name,n_ybin1,ymingate,ymaxgate);
      sprintf(hist_name,"%s:: Energy distribution of generated particles Missed",GetName().c_str());
      h.h1_Egen_Out = myROOT_utils::TH1D_checked("EgenOut",hist_name,n_ebin1,emin,emax);
      sprintf(hist_name,"%s:: X_gen vs Y_gen Missed",GetName().c_str());
      h.h2_XgYg_Out = myROOT_utils::TH2D_checked("XgYgOut",hist_name,n_xbin2,xmingate,xmaxgate,n_ybin2,ymingate,ymaxgate);

//  Histgramms for generated particles hited calorimeter
      sprintf(hist_name,"%s:: X distribution of generated particles",GetName().c_str());
      h.h1_Xgen = myROOT_utils::TH1D_checked("Xgen",hist_name,n_xbin1,xmin,xmax);
      sprintf(hist_name,"%s:: Y distribution of generated particles",GetName().c_str());
      h.h1_Ygen = myROOT_utils::TH1D_checked("Ygen",hist_name,n_ybin1,ymin,ymax);
      sprintf(hist_name,"%s:: Energy distribution of generated particles",GetName().c_str());
      h.h1_Egen = myROOT_utils::TH1D_checked("Egen",hist_name,n_ebin1,emin,emax);
      sprintf(hist_name,"%s:: X_gen vs Y_gen",GetName().c_str());
      h.h2_XgYg = myROOT_utils::TH2D_checked("XgYg",hist_name,n_xbin2,xmin,xmax,n_ybin2,ymin,ymax);

      sprintf(hist_name,"%s:: X generated not reconstructed",GetName().c_str());
      h.h1_Xgen_Norec = myROOT_utils::TH1D_checked("Xgen_Norec",hist_name,n_xbin1,xmingate,xmaxgate);
      sprintf(hist_name,"%s:: Y generated not reconstructed",GetName().c_str());
      h.h1_Ygen_Norec = myROOT_utils::TH1D_checked("Ygen_Norec",hist_name,n_ybin1,ymingate,ymaxgate);
      sprintf(hist_name,"%s:: R Distribution of generated particles",GetName().c_str());
      h.h1_Rgen = myROOT_utils::TH1D_checked("Rgen",hist_name,0,0.,std::sqrt(xmax*xmax+ymax*ymax));
      sprintf(hist_name,"%s:: Energy generated not reconstructed",GetName().c_str());
      h.h1_Egen_Norec = myROOT_utils::TH1D_checked("Egen_Norec",hist_name,n_ebin1,emin,emax);
      sprintf(hist_name,"%s:: X_gen vs Y_gen not reconstructed",GetName().c_str());
      h.h2_XgYg_Norec = myROOT_utils::TH2D_checked("XgYg_Norec",hist_name,n_xbin2,xmingate,xmaxgate,n_ybin2,ymingate,ymaxgate);
    }

    sprintf(hist_name,"%s:: X Distribution of reconstructed particles",GetName().c_str());
    h.h1_Xrec = myROOT_utils::TH1D_checked("Xrec",hist_name,n_xbin1,xmin,xmax);
    sprintf(hist_name,"%s:: Y Distribution of reconstructed particles",GetName().c_str());
    h.h1_Yrec = myROOT_utils::TH1D_checked("Yrec",hist_name,n_ybin1,ymin,ymax);
    sprintf(hist_name,"%s:: R Distribution of reconstructed particles",GetName().c_str());
    h.h1_Rrec = myROOT_utils::TH1D_checked("Rrec",hist_name,500,0.,std::sqrt(xmax*xmax+ymax*ymax));
    sprintf(hist_name,"%s:: Energy distribution of reconstructed particles",GetName().c_str());
    h.h1_Erec = myROOT_utils::TH1D_checked("Erec",hist_name,n_ebin1,emin,emax);
    sprintf(hist_name,"%s:: Summ of energies of reconstructed particles",GetName().c_str());
    h.h1_Sumrec = myROOT_utils::TH1D_checked("Sumrec",hist_name,n_ebin1,emin,emax);
    int iecut = (int) (1000.*e_cut);
    sprintf(hist_name,"%s:: Summ of energies with Ecut > %d MeV of reconstructed particles",GetName().c_str(),iecut);
    h.h1_SumrecCut = myROOT_utils::TH1D_checked("SumrecCut",hist_name,n_ebin1,emin,emax);
    sprintf(hist_name,"%s:: X_rec vs Y_rec",GetName().c_str());
    h.h2_XrYr = myROOT_utils::TH2D_checked("XrYr",hist_name,n_xbin2,xmin,xmax,n_ybin2,ymin,ymax);
    sprintf(hist_name,"%s:: X_rec vs Y_rec missed Cell? ",GetName().c_str());
    h.h2_XrYrMissedCell = myROOT_utils::TH2D_checked("XrYrMissedCell",hist_name,n_xbin2,xmin,xmax,n_ybin2,ymin,ymax);


    if( !only_reco_hist )
    {

      sprintf(hist_name,"%s:: X_gen - X_rec",GetName().c_str());
      h.h1_Xdelta = myROOT_utils::TH1D_checked("Xdelta",hist_name,100,-2*cellx,2*cellx);
      sprintf(hist_name,"%s:: Y_gen - Y_rec",GetName().c_str());
      h.h1_Ydelta = myROOT_utils::TH1D_checked("Ydelta",hist_name,100,-2*celly,2*celly);

      sprintf(hist_name,"%s:: E_gen vs E_rec",GetName().c_str());
      h.h2_EgEr = myROOT_utils::TH2D_checked("EgEr",hist_name,n_ebin2,emin,emax,n_ebin2,emin,emax);
      sprintf(hist_name,"%s:: Time_gen vs Time_rec",GetName().c_str());
      h.h2_TgTr = myROOT_utils::TH2D_checked("TgTr",hist_name,n_time2,time_min,time_max,n_time2,time_min,time_max);
      sprintf(hist_name,"%s:: X_gen vs X_rec",GetName().c_str());
      h.h2_XgXr = myROOT_utils::TH2D_checked("XgXr",hist_name,4*n_xbin2,xmin,xmax,4*n_xbin2,xmin,xmax);
      sprintf(hist_name,"%s:: Y_gen vs Y_rec",GetName().c_str());
      h.h2_YgYr = myROOT_utils::TH2D_checked("YgYr",hist_name,n_ybin2,ymin,ymax,n_ybin2,ymin,ymax);

      sprintf(hist_name,"%s:: X_gen-X_rec vs Y_gen-Y_rec",GetName().c_str());
      h.h2_XYdelta = myROOT_utils::TH2D_checked("XYdelta",hist_name,50,-2*cellx,2*cellx,50,-2*celly,2*celly);
      sprintf(hist_name,"%s:: E_gen - E_rec",GetName().c_str());
      h.h1_Edelta = myROOT_utils::TH1D_checked("Edelta",hist_name,200,-1,1);

      sprintf(hist_name,"%s:: X_gen - X_rec Hadrons",GetName().c_str());
      h.h1_XdeltaHadron = myROOT_utils::TH1D_checked("XdeltaHadron",hist_name,100,-2*cellx,2*cellx);
      sprintf(hist_name,"%s:: Y_gen - Y_rec Hadrons",GetName().c_str());
      h.h1_YdeltaHadron = myROOT_utils::TH1D_checked("YdeltaHadron",hist_name,100,-2*celly,2*celly);
      sprintf(hist_name,"%s:: X_gen - X_rec Hadrons E >40 Gev",GetName().c_str());
      h.h1_XdeltaHadron40 = myROOT_utils::TH1D_checked("XdeltaHadron40",hist_name,100,-2*cellx,2*cellx);
      sprintf(hist_name,"%s:: Y_gen - Y_rec Hadrons E >40 Gev",GetName().c_str());
      h.h1_YdeltaHadron40 = myROOT_utils::TH1D_checked("YdeltaHadron40",hist_name,100,-2*celly,2*celly);
      sprintf(hist_name,"%s:: E_gen vs E_rec Hadrons",GetName().c_str());
      h.h2_EgErHadron = myROOT_utils::TH2D_checked("EgErHadron",hist_name,n_ebin2,emin,emax,n_ebin2,emin,emax);
      sprintf(hist_name,"%s:: X_gen - X_rec Hadrons",GetName().c_str());
      h.h1_XdeltaMuon = myROOT_utils::TH1D_checked("XdeltaMuon",hist_name,100,-2*cellx,2*cellx);
      sprintf(hist_name,"%s:: Y_gen - Y_rec Muons",GetName().c_str());
      h.h1_YdeltaMuon = myROOT_utils::TH1D_checked("YdeltaMuon",hist_name,100,-2*celly,2*celly);
      sprintf(hist_name,"%s:: E_gen vs E_rec Muons",GetName().c_str());
      h.h2_EgErMuon = myROOT_utils::TH2D_checked("EgErMuon",hist_name,n_ebin2,emin,emax,n_ebin2,emin,emax);

      sprintf(hist_name,"%s:: (E_gen - E_rec) vs E_gen",GetName().c_str());
      h.h2_EdeltaEg = myROOT_utils::TH2D_checked("EdeltaEg",hist_name,50,-1,1,n_ebin1,emin,emax);
      sprintf(hist_name,"%s:: (E_gen - E_rec)/E_gen vs E_gen",GetName().c_str());
      h.h2_EdeltaN_Eg = myROOT_utils::TH2D_checked("EdeltaN_Eg",hist_name,50,-0.5,0.5,n_ebin2,emin,emax);
      sprintf(hist_name,"%s:: (E_gen - E_rec)/E_gen ",GetName().c_str());
      h.h1_EdeltaN = myROOT_utils::TH1D_checked("EdeltaN",hist_name,200,-0.5,0.5);
      sprintf(hist_name,"%s:: (X_gen - X_rec) vs X_gen",GetName().c_str());
      h.h2_XdeltaXg = myROOT_utils::TH2D_checked("XdeltaXg",hist_name,100,-2*cellx,2*cellx,n_xbin2,xmin,xmax);
      sprintf(hist_name,"%s:: (Y_gen - Y_rec) vs Y_gen",GetName().c_str());
      h.h2_YdeltaYg = myROOT_utils::TH2D_checked("YdeltaYg",hist_name,100,-2*celly,2*celly,n_ybin1,ymin,ymax);
      sprintf(hist_name,"%s:: (X_gen - X_rec) vs (Y_gen - Y_rec)",GetName().c_str());
      h.h2_XdeltaYdelta = myROOT_utils::TH2D_checked("XdeltaYdelta",hist_name,100,-2*cellx,2*cellx,100,-2*celly,2*celly);

// Cell Test (added by DM)
      sprintf(hist_name,"%s:: XY gen in Cell",GetName().c_str());
      h.h2_XgYgCellTest = myROOT_utils::TH2D_checked("XgYg_in_Cell", hist_name, cellb, -cellx, cellx, cellb, -celly, celly);
      sprintf(hist_name,"%s:: X gen vs X rec in Cell",GetName().c_str());
      h.h2_XgXrCellTest = myROOT_utils::TH2D_checked("XgXr_in_Cell", hist_name, cellb, -cellx, cellx, cellb, -celly, celly);
      sprintf(hist_name,"%s:: Y gen vs Y rec in Cell",GetName().c_str());
      h.h2_YgYrCellTest = myROOT_utils::TH2D_checked("YgYr_in_Cell", hist_name, cellb, -cellx, cellx, cellb, -celly, celly);

      sprintf(hist_name,"%s:: Xrec - Xgen as function Xgen  ",GetName().c_str());
      h.p1_DxXgCell = myROOT_utils::TProfile_checked("DxXg_in_Cell ", hist_name, cellb, -cellx, cellx );
      sprintf(hist_name,"%s:: Yrec - Ygen as function Ygen   ",GetName().c_str());
      h.p1_DyYgCell = myROOT_utils::TProfile_checked("DyYg_in_Cell ", hist_name, cellb, -cellx, cellx );
    }

    sprintf(hist_name,"%s:: XY rec in Cell",GetName().c_str());
    h.h2_XrYrCellTest = myROOT_utils::TH2D_checked("XrYr_in_Cell", hist_name, cellb, -cellx, cellx, cellb, -celly, celly);
//    Two near gamma problem
    sprintf(hist_name,"%s::  Hi2/Ndf as energy function  ",GetName().c_str());
    h.p1_Hi2Ndf = myROOT_utils::TProfile_checked("Hi2Ndf", hist_name, 50, emin,emax );
    sprintf(hist_name,"%s::  Probability as energy function ",GetName().c_str());
    h.p1_Hi2Prob = myROOT_utils::TProfile_checked("Hi2Prob", hist_name, 50, emin,emax );

    sprintf(hist_name,"%s::  Probability ",GetName().c_str());
    h.h1_Prob = myROOT_utils::TH1D_checked("Prob", hist_name, 100, 0.,1. );
    sprintf(hist_name,"%s::  Probability for 2x2 zone",GetName().c_str());
    h.h1_Prob_4 = myROOT_utils::TH1D_checked("Prob_4", hist_name, 100, 0.,1. );
    sprintf(hist_name,"%s::  Probability for 3x3 minus 2x2 zone",GetName().c_str());
    h.h1_Prob_9m4 = myROOT_utils::TH1D_checked("Prob_9m4", hist_name, 100, 0.,1. );
    sprintf(hist_name,"%s::  Probability for 5x5 minus 3x3 zone",GetName().c_str());
    h.h1_Prob_25m9 = myROOT_utils::TH1D_checked("Prob_25m9", hist_name, 100, 0.,1. );

    sprintf(hist_name,"%s::  Egamma vs Probability ",GetName().c_str());
    h.h2_ProbE = myROOT_utils::TH2D_checked("ProbE", hist_name, 100, 0.,1., 50, emin,emax );
    sprintf(hist_name,"%s::  Egamma vs Probability for 2x2 zone ",GetName().c_str());
    h.h2_ProbE_4 = myROOT_utils::TH2D_checked("ProbE_4", hist_name, 100, 0.,1., 50, emin,emax );
    sprintf(hist_name,"%s::  Egamma vs Probability for 3x3 minus 2x2 zone ",GetName().c_str());
    h.h2_ProbE_9m4 = myROOT_utils::TH2D_checked("ProbE_9m4", hist_name, 100, 0.,1., 50, emin,emax );
    sprintf(hist_name,"%s::  Egamma vs Probability for 5x5 minus 3x3 zone ",GetName().c_str());
    h.h2_ProbE_25m9 = myROOT_utils::TH2D_checked("ProbE_25m9", hist_name, 100, 0.,1., 50, emin,emax );

    sprintf(hist_name,"%s::  (Er-Et)/SDe as energy fraction function for 2x2 zone ",GetName().c_str());
    h.h2_DENorm_4 = myROOT_utils::TH2D_checked("DENorm_4", hist_name, 200,-10., 10., 50, 0, 1);
    sprintf(hist_name,"%s::  (Er-Et)/SDe as energy fraction function for 3x3 minus 2x2 zone ",GetName().c_str());
    h.h2_DENorm_9m4 = myROOT_utils::TH2D_checked("DENorm_9m4", hist_name, 200,-10., 10., 50, 0, 1);
    sprintf(hist_name,"%s::  (Er-Et)/SDe as energy fraction function for 5x5 minus 3x3 zone ",GetName().c_str());
    h.h2_DENorm_25m9 = myROOT_utils::TH2D_checked("DENorm_25m9", hist_name, 200,-10., 10., 50, 0, 1);

    sprintf(hist_name,"%s::  (Er-Et)/SDe for 2x2 zone ",GetName().c_str());
    h.h1_Pool_4 = myROOT_utils::TH1D_checked("Pool_4", hist_name, 200,-10., 10.);
    sprintf(hist_name,"%s::  (Er-Et)/SDe for 3x3 minus 2x2 zone ",GetName().c_str());
    h.h1_Pool_9m4 = myROOT_utils::TH1D_checked("Pool_9m4", hist_name, 200,-10., 10.);
    sprintf(hist_name,"%s::  (Er-Et)/SDe for 5x5 minus 3x3 zone ",GetName().c_str());
    h.h1_Pool_25m9 = myROOT_utils::TH1D_checked("Pool_25m9", hist_name, 200,-10., 10.);

    sprintf(hist_name,"%s::  Egamma vs Hi2 ",GetName().c_str());
    h.h2_Hi2E = myROOT_utils::TH2D_checked("Hi2E", hist_name, 100, 0.,100., 50, emin,emax );
    sprintf(hist_name,"%s::  Hi2 as energy function if ngen.ne.nrec  ",GetName().c_str());
    h.p1_Hi2Cells = myROOT_utils::TProfile_checked("Hi2Cells", hist_name, 50, 0, 1 );
    sprintf(hist_name,"%s::  (Er-Et)/SDe as energy fraction function ",GetName().c_str());
    h.p1_DENorm = myROOT_utils::TProfile_checked("DENorm", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  absDe/SDe as energy fraction function ",GetName().c_str());
    h.p1_DEabsNorm = myROOT_utils::TProfile_checked("DEabsNorm", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  (Eexp-Ecalc)/Ecalc as energy fraction function ",GetName().c_str());
    h.p1_DEcalc = myROOT_utils::TProfile_checked("DEcalc", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  Himax as energy fraction function if ngen.eq.nrec ",GetName().c_str());
    h.p1_HiMaxOK  = myROOT_utils::TProfile_checked("HimaxOK", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  Hi-Himax vs Himax  if ngen.eq.nrec",GetName().c_str());
    h.h2_HiHiMaxOK = myROOT_utils::TH2D_checked("HiHimaxOK", hist_name, 250, 0., 500., 250, 0., 500.);
    sprintf(hist_name,"%s::  Himax as energy fraction function if ngen.ne.nrec ",GetName().c_str());
    h.p1_HiMaxBad  = myROOT_utils::TProfile_checked("HimaxBad", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  Hi-Himax vs Himax if ngen.ne.nrec ",GetName().c_str());
    h.h2_HiHiMaxBad = myROOT_utils::TH2D_checked("HiHimaxBad", hist_name, 250, 0., 500., 250, 0., 500.);

    sprintf(hist_name,"%s::  Dist vs emin/emax 2 gamma e1+e2>5 GeV",GetName().c_str());
    h.h2_Dist2 = myROOT_utils::TH2D_checked("Dist2", hist_name, 200, 0, 400, 50, 0, 1);
    sprintf(hist_name,"%s::  Dist vs emin e1+e2>5 GeV",GetName().c_str());
    h.h2_Dist2emin = myROOT_utils::TH2D_checked("Dist2emin", hist_name, 200, 0, 400, 50, 0, 10);
    sprintf(hist_name,"%s::  DistXY 2  ",GetName().c_str());
    h.h2_Dist2Cell = myROOT_utils::TH2D_checked("Dist2Cell", hist_name, 100, -400, 400, 100, -400, 400);

//    Timing information

    sprintf(hist_name, "%s:: Timing in all cells",GetName().c_str());
    h.h1_Time = myROOT_utils::TH1D_checked("Time", hist_name, n_time_bins, time_min, time_max);
    sprintf(hist_name, "%s:: Time vs E  ",GetName().c_str());
    h.h2_TimeE = myROOT_utils::TH2D_checked("TimeE", hist_name,n_time_bins, time_min, time_max,n_ebin2,emin,emax);

//    For total hited cells and total energy at a time
    sprintf(hist_name, "%s:: Summ of Energy in all cells",GetName().c_str());
    h.h1_ETotal = myROOT_utils::TH1D_checked("Etotal", hist_name, n_ebin1, emin, emax);
    sprintf(hist_name, "%s:: Summ Energy of ADC",GetName().c_str());
    h.h1_ESumADC = myROOT_utils::TH1D_checked("ESumADC", hist_name, 1000, 0, 10000);
    sprintf(hist_name, "%s:: Number Hited Cells vs E total ",GetName().c_str());
    h.h2_NHitedEtotal = myROOT_utils::TH2D_checked("NHitedEtotal", hist_name, NCells(), 0, NCells(),50,0,100);
    sprintf(hist_name,"%s:: X_hit vs Y_hit cells",GetName().c_str());
    h.h2_XhitYhit = myROOT_utils::TH2D_checked("XhitYhit",hist_name,n_xbin2,xmin,xmax,n_ybin2,ymin,ymax);
    sprintf(hist_name,"%s:: X_hit vs Y_hit cells Ewighted",GetName().c_str());
    h.h2_XhitYhitE = myROOT_utils::TH2D_checked("XhitYhitE",hist_name,n_xbin2,xmin,xmax,n_ybin2,ymin,ymax);

//   Collect general statistic on calibrations
    sprintf(hist_name, "%s:: Calibration Factor Base ",GetName().c_str());
    h.h1_CalibFactorBase = myROOT_utils::TH1D_checked("CalibFactorBase", hist_name, 10000, 0, 1.);
    sprintf(hist_name, "%s:: Calibration Factor LED ",GetName().c_str());
    h.h1_CalibFactorLED = myROOT_utils::TH1D_checked("CalibFactorLED", hist_name, 10000, 0, 10.);
    sprintf(hist_name, "%s:: Calibration Factor Energy dependent ",GetName().c_str());
    h.p1_CalibFactorEdep = myROOT_utils::TProfile_checked("CalibFactorEdep", hist_name, 2000, 0, 200.);
    sprintf(hist_name, "%s:: Calibration Factor Time in spill dependent ",GetName().c_str());
    h.p1_CalibFactorTis = myROOT_utils::TProfile_checked("CalibFactorTis", hist_name, 200, 0, 200.);

//   Book subset histo
    for (size_t isub = 0; isub <subsets_list2store_calib_histo_.size(); isub++) // just use this list
    {
      size_t issub = subsets_list2store_calib_histo_[isub];
      char hname[132];
      sprintf(hname,"Time_%s",GetSubSets()[issub].GetName().c_str());
      sprintf(hist_name, "SubSet::%s Timing  ", GetSubSets()[issub].GetName().c_str() );
      TH1D *h1 = myROOT_utils::TH1D_checked(hname, hist_name, n_time_bins, time_min, time_max);
      h.h1_Time_subset.push_back(h1);
    }



  if(debug_ReconstructionTest) cout << " ReconstructionTest h.cells_dir->cd()  " << endl;
    ok=h.cells_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.cells_dir->cd() OK " << ok << endl;
    assert(ok);

  if(debug_ReconstructionTest)
  {
    cout << " ReconstructionTest BOOKING plenty of CELLs histos " <<  endl;
    cout << " h.h1_ErCell.size() " <<h.h1_ErCell.size() << endl;
    cout << " h.h1_ADCmCell.size() " <<h.h1_ADCmCell.size() << endl;
  }


    for( size_t i=0; i< NCells(); i++ )
    {
//      name << "Cell" << GetCellName(i);
      if( debug_ReconstructionTest ) cout << " Book for cell " << i << "CellName " << GetCellName(i) << endl;
      char name[132];
      sprintf(name,"Cell%s",GetCellName(i).c_str());
      sprintf(hist_name," Reconstructed Particle Energy in cell %s ",GetCellName(i).c_str());
//      h.h1_ErCell[i] = myROOT_utils::TH1D_checked(name.str(),hist_name,500,emin,emax);
      h.h1_ErCell[i] = myROOT_utils::TH1D_checked(name,hist_name,500,emin,emax);
      char name1[132];
//      name1 << "ADC_Cell" << GetCellName(i);
      sprintf(name1,"ADC_Cell%s",GetCellName(i).c_str());
      sprintf(hist_name," ADC Energy in main cell %s ",GetCellName(i).c_str());
      h.h1_ADCmCell[i] = myROOT_utils::TH1D_checked(name1,hist_name,4000,0,4000);
    }



  if(debug_ReconstructionTest)
  {
    cout << " ReconstructionTest BOOKING more CELLs histos " <<  endl;
    cout << " options.need_cells_XY_histo " << options.need_cells_XY_histo <<  endl;
  }
    ok=h.cells_XY_dir->cd();
    assert(ok);
    int nbinsx = (int)(((xmax-xmin)*1.0001)/cellx);
    int nbinsy = (int)(((ymax-ymin)*1.0001)/celly);

    h.h2_XYAll = myROOT_utils::TH2D_checked("XYAll"," X Y hits in cells with E gt 0.7 GeV ",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
//     if( options.need_cells_XY_histo)
//     {
//       for( size_t i=0; i< NCells(); i++ )
//       {
//         char name[132];
//         sprintf(name,"dXdYforCell%s",GetCellName(i).c_str());
//         sprintf(hist_name," (X - Xcorr ) vs (Y - Ycorr) in cell %s ",GetCellName(i).c_str());
// 	double xc = cells[i].GetX();
// 	double yc = cells[i].GetY();
//
//         h.h2_XdeltaYdeltaCell[i] = myROOT_utils::TH2D_checked(name,hist_name,nbinsx,xmin-xc,xmax-xc,nbinsy,ymin-yc,ymax-yc);
//       }
//     }

  if(debug_ReconstructionTest) cout << " ReconstructionTest END of BOOKING plenty of CELLs histos MAY BE HERE !!!!!!" << ok << endl;

  if(debug_ReconstructionTest) cout << " ReconstructionTest h.cells_time_dir->cd()  " << endl;
    ok=h.cells_time_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.cells_time_dir->cd()  " << endl;
    assert(ok);
    if( need_cells_timing_histo)
    {
      for( size_t i=0; i< NCells(); i++ )
      {
        char name[132];
        sprintf(name,"Time%s",GetCellName(i).c_str());
        sprintf(hist_name," Time in cell %s ",GetCellName(i).c_str());
        h.h1_TimeCell[i] = myROOT_utils::TH1D_checked(name,hist_name,n_time_bins,time_min,time_max);
      }
    }

    h.test_histos_booked=true;
  if(debug_ReconstructionTest) cout << " ReconstructionTest BOOKING OK  " << endl;
  }
//  end booking

  if(options.test_histo_raw_info) TestHistoRawInfo();


  if(debug_ReconstructionTest) cout << " ReconstructionTest START FILLING  " << endl;
  ok=h.test_dir->cd();
  if(debug_ReconstructionTest) cout << " ReconstructionTest h.test_dir->cd() OK " << ok << endl;
  assert(ok);

  vector<CalorimeterParticle> hint_particles_in,hint_particles_out;

  for( vector<CalorimeterParticle>::iterator it=hint_particles.begin(); it!=hint_particles.end(); it++ )
  {
     if( it->GetMainCells().size() > 0 )
     {
        hint_particles_in.push_back( *it );
     }
     else
     {
        hint_particles_out.push_back( *it );
     }
  }

  if(debug_ReconstructionTest) cout << " ReconstructionTest only_reco_hist " << only_reco_hist << endl;
  size_t ngen_out = hint_particles_out.size();

  if(debug_ReconstructionTest) cout << " ReconstructionTest ngen_out " << ngen_out << endl;
  if( !only_reco_hist )
  {
    h.h1_ngen_Out->Fill(ngen_out);

    if(ngen_out > 0 )
    {
      for( size_t it=0; it != ngen_out; it++ )
      {
        h.h1_Xgen_Out->Fill(hint_particles_out[it].GetX());
        h.h1_Ygen_Out->Fill(hint_particles_out[it].GetY());
        h.h1_Egen_Out->Fill(hint_particles_out[it].GetE());
        h.h2_XgYg_Out->Fill((float)hint_particles_out[it].GetX(),(float)hint_particles_out[it].GetY());
      }
    }
  }

  size_t ngen = hint_particles_in.size();
  size_t nrec = reco_particles.size();
  if(debug_ReconstructionTest) cout << " ReconstructionTest ngen " << ngen << " nrec " << nrec << endl;

  h.h1_nrec->Fill(nrec);

  if( !only_reco_hist )
  {
    h.h1_ngen->Fill(ngen);
    h.h2_ngen_nrec -> Fill ( (float)ngen , (float)nrec );
  }

  if( (int)nrec > options.ngam_max4histo ) return;

  if(debug_ReconstructionTest) cout << " ReconstructionTest  Fillig histogramms for generated particles " << endl;
  // Fillig histogramms for generated particles
  if( !only_reco_hist )
  {
    if(ngen > 0 )
    {
      for( size_t it=0; it != ngen; it++ )
      {
        double xg = hint_particles_in[it].GetX();
        double yg = hint_particles_in[it].GetY();
        double rg=std::sqrt(xg*xg+yg*yg);
        h.h1_Xgen->Fill(xg);
        h.h1_Ygen->Fill(yg);
        h.h1_Rgen->Fill(rg);
        h.h1_Egen->Fill(hint_particles_in[it].GetE());
        h.h2_XgYg->Fill((float)xg,(float)yg);
      }
    }
  }
  // Fillig histogramms for generated particles if reconstructed are absent
  if(nrec == 0)
  {
    if(ngen > 0 )
    {
       for( size_t it=0; it != ngen; it++ )
       {
          double xg = hint_particles_in[it].GetX();
          double yg = hint_particles_in[it].GetY();
          h.h1_Xgen_Norec->Fill((float)xg);
          h.h1_Ygen_Norec->Fill((float)yg);
          h.h1_Egen_Norec->Fill(hint_particles_in[it].GetE());
          h.h2_XgYg_Norec->Fill((float)xg,(float)yg);
       }
    }
  }
  if(debug_ReconstructionTest)
    cout << " ReconstructionTest  Fillig histogramms for reconstructed particles nrec=" << nrec << endl;
  // Fillig histogramms for reconstructed particles
  double summ_rec = 0.;
  double summ_rec_cut = 0.;
  if(nrec > 0 )
  {
    for( size_t it=0; it != nrec; it++ )  // cycle over reconstructed particles
    {
      const vector<size_t> &hitedcells =reco_particles[it].GetMainCells();
      if(reco_particles[it].HasTime())
      {
//          cout << " E " << reco_particles[it].GetE() << " Time " << reco_particles[it].GetTime() << endl;
         double time = reco_particles[it].GetTime();
         h.h1_Time->Fill(time);
         h.h2_TimeE->Fill(time,reco_particles[it].GetE());
         for (size_t isub = 0; isub <subsets_list2store_calib_histo_.size(); isub++)
         {
           size_t issub = subsets_list2store_calib_histo_[isub];
           if( hitedcells.size() <= 0 ) continue;
           int icell = hitedcells[0];
           if( GetSubSets()[issub].IsMyCell(icell) )
	   {
	     if( h.h1_Time_subset[isub] == NULL )
	     {
               cerr << GetName() << " h1_Time_subset " << isub << "  issub "<<  issub << " no histo! " << endl;
	       exit(1);
	     }
	     h.h1_Time_subset[isub]->Fill( time);
           }
         }
      }
      if(debug_ReconstructionTest)
         cout << " it=" << it << endl;
      double xr = reco_particles[it].GetX();
      double yr = reco_particles[it].GetY();
      if(debug_ReconstructionTest)
         cout << " it=" << it << " xr=" << xr << " yr=" << yr << " e=" << reco_particles[it].GetE() << endl;
      double rr = std::sqrt(xr*xr+yr*yr);
      h.h1_Xrec->Fill(xr);
      h.h1_Yrec->Fill(yr);
      h.h1_Rrec->Fill(rr);
      h.h1_Erec->Fill(reco_particles[it].GetE());
      summ_rec += reco_particles[it].GetE();
      if( reco_particles[it].GetE() > e_cut ) summ_rec_cut += reco_particles[it].GetE();
      h.h2_XrYr->Fill((float)xr,(float)yr,1);

      if(debug_ReconstructionTest)
        cout << " Cycle over hitted cells " << reco_particles[it].GetMainCells().size() << endl;
      for(size_t i=0; i< hitedcells.size(); i++)
      {
        size_t incell=hitedcells[i];
        if(reco_particles[it].HasTime())
        {
          double time = reco_particles[it].GetTime();
          if(need_cells_timing_histo)
          {
            ok=( h.h1_TimeCell[incell] != NULL);
            assert( ok);
            h.h1_TimeCell[incell]->Fill(time);
          }
        }

        ok=( h.h1_ErCell[incell] != NULL);
        assert(ok);
        h.h1_ErCell[incell]->Fill(reco_particles[it].GetE());

        for( size_t j=0; j!=signals_.size(); j++ )
        {
          if(incell == signals_[j].GetCellIdx())
          {
            ok=( h.h1_ADCmCell[incell] != NULL);
            assert( ok);
            if(signals_[j].HasAmplitude() ) h.h1_ADCmCell[incell]->Fill(signals_[j].GetAmplitude());
            break;
          }
        }
      }
      if(debug_ReconstructionTest)
        cout << " Cycle over hitted cells ok " << endl;

    }  // end of cycle over reconstructed particles
  } // nrec > 0

  if(debug_ReconstructionTest) cout << " ReconstructionTest Check of internal reconstruction info  still nrec=" << nrec <<endl;
   bool hist_internal_info = options.fill_internal_histos;
//    bool internal_histo_specific_demand = false;
//    if( hist_internal_info  )
//    {
//      if(debug_ReconstructionTest) cout << " ReconstructionTest Hist internal info" << endl;
//      if( nrec == 1 )
//      {
//        double emin_specific_demand = 8.5;
//        double emax_specific_demand = 11.5;
//        if( reco_particles[0].GetE() >= emin_specific_demand &&
//            reco_particles[0].GetE() < emax_specific_demand )
//        {
//          internal_histo_specific_demand = true;
//        }
//      }
//
//    }

   if( options.fill_internal_histos )
   {
     if(debug_ReconstructionTest) cout << " ReconstructionTest Now realy goto HistInternalInfo() " << endl;
//     cerr << " Some bug in HistInternalInfo() " << endl;
     HistInternalInfo();
     if(debug_ReconstructionTest) cout << " ReconstructionTest HistInternalInfo() OK " << endl;
   }

  if(debug_ReconstructionTest) cout << " ReconstructionTest plane Hist Sumrec and SumrecCut " << endl;
  h.h1_Sumrec->Fill(summ_rec);
  h.h1_SumrecCut->Fill(summ_rec_cut);
  if(debug_ReconstructionTest) cout << " ReconstructionTest Hist calib_data_store signals_.size() = " <<
                                                                                   signals_.size() << endl;
  // Use cells_XY_histo to search correlations
  for(size_t i=0; i< signals_.size(); i++)
  {
    size_t incell=signals_[i].GetCellIdx();
    double e = signals_[i].GetEnergy();
    if( e < 0.7 ) continue;
    if( incell >= NCells() )
    {
       cerr << " ERROR detected in ReconstructionTest WRONG adress in cell = " << incell << endl;
    }
    double x = cells[incell].GetX();
    double y = cells[incell].GetY();
    if(h.h2_XYAll != NULL) h.h2_XYAll->Fill(x,y);
  }
  if(debug_ReconstructionTest) cout << " ReconstructionTest Hist calib_data_store OK " << endl;

//   if( options.need_cells_XY_histo)
//   {
//     for( size_t it=0; it != nrec; it++ )
//     {
//       size_t incell0 = reco_particles[it].GetMainCells()[0];
//       double x0 = cells[incell0].GetX();
//       double y0 = cells[incell0].GetY();
//       for(size_t i=0; i< calib_data_store.size(); i++)
//       {
//         size_t incell=calib_data_store[i].GetAddress();
//         double e = calib_data_store[i].GetAmplitude();
// 	if( incell == incell0 ) continue;
// 	if( e < 0.7 ) continue;
//         double x = cells[incell].GetX();
//         double y = cells[incell].GetY();
//         if(h.h2_XdeltaYdeltaCell[incell0] != NULL) h.h2_XdeltaYdeltaCell[incell0]->Fill(x-x0,y-y0);
//       }
//     }
//   }
  if(debug_ReconstructionTest) cout << " ReconstructionTest Start  Linking generated particles to reconstructed particles ngen = " << ngen << endl;
  // Linking generated particles to reconstructed particles
  if( !only_reco_hist )
  {
    if( nrec > 0 && ngen > 0 )
    {
      for( size_t itgen=0; itgen != ngen; itgen++ )
      {
         double xgen0 = hint_particles_in[itgen].GetX();
         double ygen0 = hint_particles_in[itgen].GetY();
         double zgen = hint_particles_in[itgen].GetZ();
         double egen = hint_particles_in[itgen].GetE();
         double axgen = hint_particles_in[itgen].GetAngleX();
         double aygen = hint_particles_in[itgen].GetAngleY();
         for( size_t itrec=0; itrec != nrec; itrec++ )
         {
            double xrec = reco_particles[itrec].GetX();
            double yrec = reco_particles[itrec].GetY();
            double zrec = reco_particles[itrec].GetZ();
            double erec = reco_particles[itrec].GetE();
//           if( fabs(ygen0-yrec) < 4*celly && fabs(xgen0-xrec) < 4*cellx )
// 	  {
//             cout << GetName() << " er= " << erec << " xr= " << xrec << " yr= " << yrec << " zr= " << zrec << endl;
//             cout << " eg= " << egen << " xg= " << xgen0 << " yg= " << ygen0 << " zg= " << zgen <<
// 	                                                " ax= " << axgen << " ay= " << aygen << endl;
//             cout <<  " xgnew= " << xgen0 + axgen*(zrec - zgen) << " ygnew= " << ygen0 + aygen*(zrec - zgen) << endl;
// 	  }
  	    double xgen = xgen0 + axgen*(zrec - zgen);
	    double ygen = ygen0 + aygen*(zrec - zgen);

            if( fabs(ygen-yrec) < 2*celly ) h.h1_Xdelta->Fill(xgen-xrec);
            if( fabs(xgen-xrec) < 2*cellx ) h.h1_Ydelta->Fill(ygen-yrec);
            h.h2_XYdelta->Fill((float)(xgen-xrec),(float)(ygen-yrec));

//  move into  hint_gate            h.h2_EgEr->Fill(egen,erec,1);
            if( fabs(ygen-yrec) < 2*celly ) h.h2_XgXr->Fill(xgen,xrec,1);
            if( fabs(xgen-xrec) < 2*cellx ) h.h2_YgYr->Fill(ygen,yrec,1);

            h.h2_XYdelta->Fill((float)(xgen-xrec),(float)(ygen-yrec));
            h.h1_Edelta->Fill(egen-erec);

            h.h2_EdeltaEg->Fill((float)(egen-erec),(float)egen);
            h.h2_EdeltaN_Eg->Fill((float)(1.-erec/egen),(float)egen);
            h.h1_EdeltaN->Fill((float)(1.-erec/egen));
            if( fabs(ygen-yrec) < 2*celly ) h.h2_XdeltaXg->Fill((float)(xgen-xrec),(float)xgen);
            if( fabs(xgen-xrec) < 2*cellx ) h.h2_YdeltaYg->Fill((float)(ygen-yrec),(float)ygen);

            h.h2_XdeltaYdelta->Fill((float)(xgen-xrec),(float)(ygen-yrec));
//           if( options.need_cells_XY_histo)
//           {
//             vector<size_t> hitedcells =reco_particles[itrec].GetMainCells();
//             for(size_t i=0; i< hitedcells.size(); i++)
//             {
//               size_t incell=hitedcells[i];
//               if(h.h2_XdeltaYdeltaCell[incell] != NULL) h.h2_XdeltaYdeltaCell[incell]->Fill((float)(xgen-xrec),(float)(ygen-yrec));
//             }
//           }

// #warning "WARNING!! ReconstructionTest:: hardcoded constant in particles linking "
            double hint_gate =1.2*cellx;
            if( abs(xgen-xrec) < hint_gate/2 && abs(ygen-yrec) < hint_gate/2 )
            {
              int pid    = hint_particles_in[itgen].GetID();
              vector<size_t> hitedcells =hint_particles_in[itgen].GetMainCells();
              h.h2_EgEr->Fill(egen,erec,1);
              double trec = reco_particles[itrec].GetTime();
              double tgen = hint_particles_in[itgen].GetTime();
              h.h2_TgTr->Fill(tgen,trec,1);
//            cout << " Linked with Particle ID " << pid << " in cell " << incell << endl;
//            cout << " Erec = " << erec << " Egen = " << egen << endl;
//            cout << " Xrec = " << xrec << " Xgen = " << xgen << endl;
//            cout << " Yrec = " << yrec << " Ygen = " << ygen << endl;
              double sigma_egen = std::sqrt(0.65*0.65*egen + 0.07*egen*egen);
//            double sigma_egen = 0.65*std::sqrt(egen) + 0.07*egen;
              if( erec > (egen - 2*sigma_egen)  ) // hadron
              {
                h.h1_XdeltaHadron->Fill((float)(xgen-xrec));
                h.h1_YdeltaHadron->Fill((float)(ygen-yrec));
                if( egen>40)
                {
                  h.h1_XdeltaHadron40->Fill((float)(xgen-xrec));
                  h.h1_YdeltaHadron40->Fill((float)(ygen-yrec));
                }
                h.h2_EgErHadron->Fill(egen,erec,1);
              }
              else    // muon
              {
                h.h1_XdeltaMuon->Fill((float)(xgen-xrec));
                h.h1_YdeltaMuon->Fill((float)(ygen-yrec));
                h.h2_EgErMuon->Fill(egen,erec,1);
              }

              for(size_t i=0; i< hitedcells.size(); i++)
              {
                size_t incell=hitedcells[i];

//            cout << " Linked with Particle ID " << pid << " in cell " << incell << endl;
//            cout << " Erec = " << erec << " Egen = " << egen << endl;
//            cout << " Xrec = " << xrec << " Xgen = " << xgen << endl;
//            cout << " Yrec = " << yrec << " Ygen = " << ygen << endl;
            //  Link with charged particle
                if( pid == 9 )
                {
                 // Energy around
                   double e_around = 0.;
                   for (vector<size_t>::const_iterator it=cells[incell].GetNeighbors().begin()+1;
                        it!=cells[incell].GetNeighbors().end(); it++) {
                      for (vector<CellDataRaw>::const_iterator itSig=signals_.begin(); itSig!=signals_.end(); itSig++)
                        if ( itSig->GetCellIdx()==(*it) )
                          e_around += itSig->GetEnergy();
                   }
//                  cout <<  " Energy in cell "  << real_amplitudes[incell].first
//                                            <<          " and around "      <<  e_around << endl;
                }
              }
            } // linked
         } // reco cycle
      if(debug_ReconstructionTest) cout << " ReconstructionTest OPYAT' zaputalsya VDYM !!! " << endl;
      }// generated cycle
    } // check ngen>0 && nrec>0
  } // check only_reco_hist option

  if(debug_ReconstructionTest) cout << " ReconstructionTest Propustil dve skobki gde ya ne znaju !! " << endl;

// Filling XYCellTest
  if(debug_ReconstructionTest) cout << " Filling XYCellTest yakoby!!  " << endl;
  double xr, yr, zr, x, y;
  int n;
  if(debug_ReconstructionTest) cout << " ngen =  " << ngen << " hint_particles_in.size() = " <<
                                                                       hint_particles_in.size() << endl;

  if( ngen != hint_particles_in.size() )
  {
     cerr << " ERROR in ReconstructionTest Polnyi OUT!! ngen =  " << ngen << " hint_particles_in.size() = " <<
                                                                       hint_particles_in.size() << endl;
     exit(1);
  }

  if( !only_reco_hist )
  {
    if(ngen > 0 )
    {
      for( size_t it=0; it != ngen; it++ )
      {
        xr = hint_particles_in[it].GetX();
        yr = hint_particles_in[it].GetY();
        zr = hint_particles_in[it].GetZ();

        vector <double> xhint;
        xhint.push_back(xr);
        xhint.push_back(yr);
        xhint.push_back(zr);
        n = FindCell (xhint);
//      n = FindCell (xr, yr, zr);
        x = xr - cells[n].GetX();
        y = yr - cells[n].GetY();
        h.h2_XgYgCellTest->Fill( (float) x, (float) y );
      }
    }
  }
  if(debug_ReconstructionTest) cout << " Calorimeter::ReconstructionTest nrec= " << nrec << endl;

  if(nrec > 0 )
  {
    if( nrec != reco_particles.size() )
    {
      cerr << " ERROR in ReconstructionTest Polnyi OUT!! nrec =  " << nrec << " reco_particles.size() = " <<
                                                                       reco_particles.size() << endl;
      exit(1);
    }
    for( size_t it=0; it != nrec; it++ )
    {
      xr = reco_particles[it].GetX();
      yr = reco_particles[it].GetY();
      zr = reco_particles[it].GetZ();

      vector <double> vxr;
      vxr.push_back(xr);
      vxr.push_back(yr);
      vxr.push_back(zr);

      n = FindCell (vxr);
//      n = FindCell (xr, yr, zr);
      if( n < 0 || n >= (int)NCells() )
      {
        h.h2_XrYrMissedCell->Fill( (float) xr, (float) yr );
        if( print_missed_cell_warning )
        {
          cerr << " WARNING!! in " << GetName() <<
             " Calorimeter::ReconstructionTest Non existig cell " << n << endl;
          reco_particles[it].Print();
          cerr << " Try to continue " << endl;
        }
        continue;
//         exit(1);
      }
      x = xr - cells[n].GetX();
      y = yr - cells[n].GetY();
      h.h2_XrYrCellTest->Fill( (float) x, (float) y );
    }
  }

  if(debug_ReconstructionTest) cout << " Calorimeter::ReconstructionTest cycle nrecXngen " <<endl;
  if( !only_reco_hist )
  {
  double xg, yg, zg, axg, ayg;
//  if(ngen == 1 && nrec == 1)
    for( size_t itr=0; itr != nrec; itr++ )
    {
      for( size_t itg=0; itg != ngen; itg++ )
      {
        xr = reco_particles[itr].GetX();
        yr = reco_particles[itr].GetY();
        zr = reco_particles[itr].GetZ();

        xg = hint_particles_in[itg].GetX();
        yg = hint_particles_in[itg].GetY();
        zg = hint_particles_in[itg].GetZ();
        axg = hint_particles_in[itg].GetAngleX();
        ayg = hint_particles_in[itg].GetAngleY();
//	xg = xg + axg*(zr - zg);
//	yg = yg + ayg*(zr - zg);

	vector <double> vxg;
	vxg.push_back(xg);
	vxg.push_back(yg);
	vxg.push_back(zg);

        n = FindCell ( vxg );
//        n = FindCell (xg, yg, zg);
	if( n <= 0 ) continue;
        xg = xg - cells[n].GetX();
        yg = yg - cells[n].GetY();

        // The same cell position should subtracted n = FindCell (xr, yr, zr);
        xr = xr - cells[n].GetX();
        yr = yr - cells[n].GetY();

	if( fabs(xg-xr) < 2.*cells[n].GetCellType().GetSizeX() &&
	      fabs(yg-yr) < 2.*cells[n].GetCellType().GetSizeY()  )
	{
          h.h2_XgXrCellTest->Fill( (float) xg, (float) xr );
          h.h2_YgYrCellTest->Fill( (float) yg, (float) yr );

          h.p1_DxXgCell->Fill( xg, xr-xg );
          h.p1_DyYgCell->Fill( yg, yr-yg );
        }
      }
    }
  } // check only_reco_hist option

  if(debug_ReconstructionTest) cout << " Calorimeter::ReconstructionTest Two near gamma problem " <<  endl;

//    Two near gamma problem
// # warning " Calorimeter::ReconstructionTest:: This cut is only for GAMS valid !!! "
  int n_near_hole=0;
  bool near_hole=false;
  int cell_near_hole=-1;
  for( size_t it=0; it != nrec; it++ )
  {
    double x=reco_particles[it].GetX();
    double y=reco_particles[it].GetY();
    double e=reco_particles[it].GetE();
    if(fabs(x) < 40. && fabs(y) < 40.) n_near_hole++;
    if(e > 5 && fabs(x) < 40. && fabs(y) < 40.)
    {
      near_hole=true;
      if(cell_near_hole == -1) cell_near_hole=(int)it;
      else if(cell_near_hole > 0 && e > reco_particles[cell_near_hole].GetE()) cell_near_hole=(int)it;
    }
  }
  if( nrec == 3 && n_near_hole < 2 )
  {
    for( size_t it=0; it != nrec; it++ )
    {
      double x=reco_particles[it].GetX();
      double y=reco_particles[it].GetY();
      double e=reco_particles[it].GetE();
//       if(cell_near_hole != -1 && it != cell_near_hole)
//       {
        double dist = std::sqrt (x*x+y*y);
        h.h2_Dist2emin->Fill( float(dist), e,1);
//       }
//       if(cell_near_hole != -1 && it == cell_near_hole)
//       {
        h.h2_Dist2Cell->Fill( float(x), float(y),1);
//       }
    }

    for( size_t it1=0; it1 != nrec-1; it1++ )
    {
      double x1=reco_particles[it1].GetX();
      double y1=reco_particles[it1].GetY();
      double e1=reco_particles[it1].GetE();
      for( size_t it2=it1+1; it2 != nrec; it2++ )
      {
         double x2=reco_particles[it2].GetX();
         double y2=reco_particles[it2].GetY();
         double e2=reco_particles[it2].GetE();
         double dist = std::sqrt ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

         double fr=e1/e2;
         double emax = e2;
         double emin = e1;
         double xn = x1;
         double yn = y1;
         double xm = x2;
         double ym = y2;
         if(e1 > e2)
         {
           fr=e2/e1;
           emax = e1;
           emin = e2;
           xn = x2;
           yn = y2;
           xm = x1;
           ym = y1;
         }

         if(emax >5) h.h2_Dist2->Fill( float(dist), fr,1);
//          if(emax >5) h.h2_Dist2emin->Fill( float(dist), emin,1);

//          if(reco_particles[it1].GetMainCells().size() <=0 ||
//             reco_particles[it1].GetMainCells().size() >  1 ) break;
//          if(reco_particles[it2].GetMainCells().size() <= 0 ||
//             reco_particles[it2].GetMainCells().size() >  1 ) break;
//          int n1 = reco_particles[it1].GetMainCells()[0];
//          int n2 = reco_particles[it2].GetMainCells()[0];
//          double xcell1 = cells[n1].GetX();
//          double ycell1 = cells[n1].GetY();
//          double xcell2 = cells[n2].GetX();
//          double ycell2 = cells[n2].GetY();
//         if(emin < 1) h.h2_Dist2Cell->Fill( float(xn), float(yn),1);
//          if(emin < 1) h.h2_Dist2Cell->Fill( float(xm), float(ym),1);

//          if(e2 > 5 && e1 < 0.6)
//          {
//            h.h2_Dist2Cell->Fill( float(xcell2-xcell1), float(ycell2-ycell1),1);
//          }
//          if(e1 > 5 && e2 < 0.6)
//          {
//            h.h2_Dist2Cell->Fill( float(xcell1-xcell2), float(ycell1-ycell2),1);
//          }


      }
    }
  }

  // Finish   Two near gamma problem
  if(debug_ReconstructionTest) cout << " Calorimeter::ReconstructionTest Two near gamma problem OK " <<  endl;

  if(debug_ReconstructionTest) cout << " Calorimeter::ReconstructionTest Filling total histograms " <<  endl;
  // Filling "total" histograms
  int NTotal = 0;
  double ETotal = 0.0;
  double ESumADC = 0.0;

  for( size_t it1=0; it1!=signals_.size(); it1++ )
  {
     double e=signals_[it1].GetEnergy();
     int icell = signals_[it1].GetCellIdx();
     ETotal += e;
     if (e >= 0.01) NTotal ++;
     double xcell = cells[icell].GetX();
     double ycell = cells[icell].GetY();
     h.h2_XhitYhit->Fill(xcell,ycell);
     h.h2_XhitYhitE->Fill(xcell,ycell,e);
     if(signals_[it1].HasAmplitude() ) ESumADC += signals_[it1].GetAmplitude();
     h.h1_CalibFactorBase->Fill( GetCalibrationFactorBase( icell) );
     h.h1_CalibFactorLED->Fill( GetCalibrationFactorLED( icell) );
     h.p1_CalibFactorEdep->Fill( e, GetCalibrationFactorEdep( icell, e) );
//     h.p1_CalibFactorTis->Fill( tis, GetCalibrationFactorEdep( icell, tis) );
  }
  h.h1_ETotal->Fill (ETotal);
  h.h2_NHitedEtotal->Fill (NTotal,ETotal);
  h.h1_ESumADC->Fill (ESumADC);

  if(debug_ReconstructionTest) cout << " NOW ReconstructionTest dir_save->cd(); " << GetName() << endl;

  dir_save->cd();

  if(debug_ReconstructionTest) cout << " ReconstructionTest OK " << GetName() << endl;
}

////////////////////////////////////////////////////////////////////////////////

void  Calorimeter::TestHistoRawInfo( void )
{
}

////////////////////////////////////////////////////////////////////////////////

bool  Calorimeter::CheckBoundaryX4MakeProfiles( double xb, const vector < double > &boundary,
                                                                 const std::vector<size_t> &long_list ) const
{
  bool ok = true;
  double tolerance = options.tolerance_for_nearby_cells;
  for( int ib=0; ib< (int)boundary.size(); ib++ )
  {
    if( fabs(xb-boundary[ib] ) < tolerance )
    {
      ok = false;
      break;
    }
  }
  if( ok == false ) return ok;
  if( ok )
  {
    int nl=0;
    int nr=0;
    for( int ibc=0; ibc< (int)long_list.size(); ibc++ )
    {
      size_t cell = long_list[ibc];
      double xcell = cells[cell].GetX();
//      double ycell = cells[cell].GetY();
//      double zcell = cells[cell].GetZ();
      double sxcell = cells[cell].GetCellType().GetSizeX();
//      double sycell = cells[cell].GetCellType().GetSizeY();

      double xl = xcell-sxcell/2.;
      double xr = xcell+sxcell/2.;
//      double yl = ycell-sycell/2.;
//      double yr = ycell+sycell/2.;

      if( xl+tolerance < xb && xr-tolerance > xb )
      {
        ok = false;
        break;
      }
      if( xcell < xb ) nl++;
      if( xcell > xb ) nr++;
    }
    if( ok && (nl == 0 || nr == 0) ) ok= false;
  }
  return ok;
}

bool  Calorimeter::CheckBoundaryY4MakeProfiles( double xb, const vector < double > &boundary,
                                                                 const std::vector<size_t> &long_list ) const
{
  bool ok = true;
  double tolerance = options.tolerance_for_nearby_cells;
  for( int ib=0; ib< (int)boundary.size(); ib++ )
  {
    if( fabs(xb-boundary[ib] ) < tolerance )
    {
      ok = false;
      break;
    }
  }
  if( ok == false ) return ok;
  if( ok )
  {
    int nl=0;
    int nr=0;
    for( int ibc=0; ibc< (int)long_list.size(); ibc++ )
    {
      size_t cell = long_list[ibc];
//      double xcell = cells[cell].GetX();
      double ycell = cells[cell].GetY();
//      double zcell = cells[cell].GetZ();
//      double sxcell = cells[cell].GetCellType().GetSizeX();
      double sycell = cells[cell].GetCellType().GetSizeY();

//      double xl = xcell-sxcell/2.;
//      double xr = xcell+sxcell/2.;
      double yl = ycell-sycell/2.;
      double yr = ycell+sycell/2.;

      if( yl+tolerance < xb && yr-tolerance > xb )
      {
        ok = false;
        break;
      }
      if( ycell < xb ) nl++;
      if( ycell > xb ) nr++;
    }
    if( ok && (nl == 0 || nr == 0) ) ok= false;
  }
  return ok;
}

void  Calorimeter::MakeProfiles( void )
{
//   bool debug = true;
  bool debug = false;
  if( !options.fill_histos )
    return;
  if( !options.make_profiles )
    return;
  if( profiles_hist_ == NULL )
  {
    profiles_hist_ = new ProfilesHisto();
    if( profiles_hist_ == NULL )
    {
      cerr << " ERROR in Calorimeter::ProfilesHisto: No memory for histograms !!" << endl;
      exit(1);
    }
    ProfilesHisto &h = *profiles_hist_;
    if( debug ) cout << " Calorimeter::ProfilesHisto Start Booking " << endl;
    if(!gDirectory->cd("/")) assert(false);
    char hist_name[132];
    char dir_name[111];

    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    TDirectory *top_dir = myROOT_utils::TDirectory_checked(dir_name);
    bool ok = top_dir->cd();
    if( !ok ) assert(false);
    sprintf(dir_name,"%s_ProfilesHisto",GetName().c_str());
    TDirectory *profiles_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok = profiles_dir->cd();
    if( !ok ) assert(false);

    double dmaxx = 5.*GetMaxCellSizeX();
    double dminx = -dmaxx;
    double dmaxy = 5.*GetMaxCellSizeY();
    double dminy = -dmaxy;


    sprintf(hist_name,"%s:: Electrons Profile X ",GetName().c_str());
    h.ProfileX_Electrons = myROOT_utils::TProfile_checked("ElectronsProfileX", hist_name, 200, dminx,dmaxx );

    sprintf(hist_name,"%s:: Electrons Profile Y ",GetName().c_str());
    h.ProfileY_Electrons = myROOT_utils::TProfile_checked("ElectronsProfileY", hist_name, 200, dminy,dmaxy );

    sprintf(hist_name,"%s:: Muons Profile X ",GetName().c_str());
    h.ProfileX_Muons = myROOT_utils::TProfile_checked("MuonsProfileX", hist_name, 200, dminx,dmaxx );

    sprintf(hist_name,"%s:: Muons Profile Y ",GetName().c_str());
    h.ProfileY_Muons = myROOT_utils::TProfile_checked("MuonsProfileY", hist_name, 200, dminy,dmaxy );

    sprintf(hist_name,"%s:: Pions Profile X ",GetName().c_str());
    h.ProfileX_Pions = myROOT_utils::TProfile_checked("PionsProfileX", hist_name, 200, dminx,dmaxx );

    sprintf(hist_name,"%s:: Pions Profile Y ",GetName().c_str());
    h.ProfileY_Pions = myROOT_utils::TProfile_checked("PionsProfileY", hist_name, 200, dminy,dmaxy );

    for ( int ien=0; ien < (int)h.ENERGY_BINS_SIZE; ien++ )
      for ( int iang=0; iang < (int)h.BIG_ANGLE_BINS_SIZE; iang++ )
      {
        char name[132];
        sprintf(name,"ElectronsProfileX_En_%d_BigAng_%d",int(1000.1*h.energy_range[ien]),int(1000.1*h.big_angle_range[iang]));
        sprintf(hist_name," %s:: Electrons Profile X BigAngle(mrad) %d -- %d Energy(MeV) %d -- %d ",
            GetName().c_str(), int(1000.*h.big_angle_range[iang]),int(1000.*h.big_angle_range[iang+1]),
            int(1000.*h.energy_range[ien]),int(1000.*h.energy_range[ien+1]) );
        h.ProfileX_BigAngle_Energy_Electrons[iang][ien] = myROOT_utils::TProfile_checked(name, hist_name, 200, dminx,dmaxx );

        sprintf(name,"ElectronsProfileX_En_%d_SmallAng_%d",int(1000.1*h.energy_range[ien]),int(1000.1*h.small_angle_range[iang]));
        sprintf(hist_name," %s:: Electrons Profile X SmallAngle(mrad) %d -- %d Energy(MeV) %d -- %d ",
            GetName().c_str(), int(1000.*h.small_angle_range[iang]),int(1000.*h.small_angle_range[iang+1]),
            int(1000.*h.energy_range[ien]),int(1000.*h.energy_range[ien+1]) );
        h.ProfileX_SmallAngle_Energy_Electrons[iang][ien] = myROOT_utils::TProfile_checked(name, hist_name, 200, dminx,dmaxx );

        sprintf(name,"ElectronsProfileY_En_%d_BigAng_%d",int(1000.1*h.energy_range[ien]),int(1000.1*h.big_angle_range[iang]));
        sprintf(hist_name," %s:: Electrons Profile Y BigAngle(mrad) %d -- %d Energy(MeV) %d -- %d ",
            GetName().c_str(), int(1000.*h.big_angle_range[iang]),int(1000.*h.big_angle_range[iang+1]),
            int(1000.*h.energy_range[ien]),int(1000.*h.energy_range[ien+1]) );
        h.ProfileY_BigAngle_Energy_Electrons[iang][ien] = myROOT_utils::TProfile_checked(name, hist_name, 200, dminx,dmaxx );

        sprintf(name,"ElectronsProfileY_En_%d_SmallAng_%d",int(1000.1*h.energy_range[ien]),int(1000.1*h.small_angle_range[iang]));
        sprintf(hist_name," %s:: Electrons Profile Y SmallAngle(mrad) %d -- %d Energy(MeV) %d -- %d ",
            GetName().c_str(), int(1000.*h.small_angle_range[iang]+0.1),int(1000.*h.small_angle_range[iang+1]+0.1),
            int(1000.*h.energy_range[ien]),int(1000.*h.energy_range[ien+1]) );
        h.ProfileY_SmallAngle_Energy_Electrons[iang][ien] = myROOT_utils::TProfile_checked(name, hist_name, 200, dminx,dmaxx );
      }

    if( debug ) cout << " Calorimeter::ProfilesHisto Booking  OK " << endl;
  }
  ProfilesHisto &h = *profiles_hist_;

  if( debug ) cout << " Calorimeter::MakeProfiles " << GetName() << " debuging " << endl;
  for( vector<CalorimeterParticle>::const_iterator phint= hint_particles.begin(); phint!= hint_particles.end(); phint++ )
  {
    int mcell = -1;
    if( (phint)->GetMainCells().empty() ) continue; // the particle finaly missed calorimeter
    mcell = (phint)->GetMainCells()[0];
    if( mcell < 0 )
    {
      cerr << " Calorimeter::MakeProfiles " << GetName() << " cell does not exist ...?? continue " << endl;
      continue;
    }

    double xhint = (phint)->GetX();
    double yhint = (phint)->GetY();
    double zhint = (phint)->GetZ();
    CalorimeterParticle::ParticleID id_hint = (phint)->GetID();
    double ehint = (phint)->GetE();
    double axhint = (phint)->GetAngleX();
    double ayhint = (phint)->GetAngleY();

    int index_e = h.IndexEnergy(ehint);
    int index_big_ax = h.IndexBigAngle(axhint);
    int index_small_ax = h.IndexSmallAngle(axhint);
    int index_big_ay = h.IndexBigAngle(ayhint);
    int index_small_ay = h.IndexSmallAngle(ayhint);

    double distmin = 1000000000.;
    double distminx = 1000000000.;
    double distminy = 1000000000.;
    for( vector<CalorimeterParticle>::const_iterator pohint= hint_particles.begin(); pohint!= hint_particles.end(); pohint++ )
    {
      if( phint == pohint ) continue;
      double zohint = (pohint)->GetZ();
      double xohint = (pohint)->GetX()+(pohint)->GetAngleX()*(zhint - zohint);
      double yohint = (pohint)->GetY()+(pohint)->GetAngleY()*(zhint - zohint);
      double dx = xhint-xohint;
      double dy = yhint-yohint;
      double dist = sqrt( dx*dx + dy*dy );
      if( dist < distmin )
      {
        distmin = dist;
        distminx = dx;
        distminy = dy;
      }
    }

    if( debug ) cout << " distmin " << distmin << endl;

//    double xmcell = cells[mcell].GetX();
//    double ymcell = cells[mcell].GetY();
//    double zmcell = cells[mcell].GetZ();
    double sxmcell = cells[mcell].GetCellType().GetSizeX();
    double symcell = cells[mcell].GetCellType().GetSizeY();

    if( debug ) cout << " cell size x " << sxmcell << " cell size y " << symcell << endl;
    if( debug ) cout << " sdistmin " << distmin/sxmcell << " distminx " << fabs(distminx)/sxmcell << " distminy " << fabs(distminy)/symcell << endl;

//    if( distmin/sxmcell < 5. ) continue;
    if( debug ) cout << " OK " << !(fabs(distminx)/sxmcell > 5. || fabs(distminy)/symcell > 5.) << endl;
    if( !(fabs(distminx)/sxmcell > 5. || fabs(distminy)/symcell > 5.) ) continue;
//    if( !WellInActiveArea( xmcell, ymcell, zmcell ) ) continue;


    OneParticleResponse phint_responce ((*phint),this);
    const std::vector<size_t> short_list = phint_responce.GetCells();
    phint_responce.ExtendList();                       // dont care about list size
    phint_responce.CalculateExpectedResponse();
    const std::vector<size_t> long_list = phint_responce.GetCells();
    if( debug )
    {
      cout << " Main cell " << mcell <<" List size " << long_list.size() << endl;
      for ( int i=0; i< (int)long_list.size(); i++ )
      {
        if( i < (int)short_list.size() ) cout << " List short" << short_list[i];
        cout << " long " << long_list[i] << endl;
      }
    }
    vector < double > xboundary;
    vector < pair < double,double > > exboundary;
    vector < double > yboundary;
    vector < pair < double,double > > eyboundary;
//    double tolerance = options.tolerance_for_nearby_cells;
    for( int ic=0; ic< (int)long_list.size(); ic++ )
    {
      int cell = long_list[ic];
      double xcell = cells[cell].GetX();
      double ycell = cells[cell].GetY();
//      double zcell = cells[cell].GetZ();
      double sxcell = cells[cell].GetCellType().GetSizeX();
      double sycell = cells[cell].GetCellType().GetSizeY();

      double xl = xcell-sxcell/2.;
      double xr = xcell+sxcell/2.;
      double yl = ycell-sycell/2.;
      double yr = ycell+sycell/2.;

      if( CheckBoundaryX4MakeProfiles(xl, xboundary, long_list ) ) xboundary.push_back(xl);
      if( CheckBoundaryX4MakeProfiles(xr, xboundary, long_list ) ) xboundary.push_back(xr);
      if( CheckBoundaryY4MakeProfiles(yl, yboundary, long_list ) ) yboundary.push_back(yl);
      if( CheckBoundaryY4MakeProfiles(yr, yboundary, long_list ) ) yboundary.push_back(yr);

    }
    if( debug ) cout << " Boundary X size " << xboundary.size() << endl;
    if( debug ) cout << " Boundary Y size " << yboundary.size() << endl;
    double etot=0.;
    for (size_t ic=0; ic<long_list.size(); ic++)
    {
      size_t cell = long_list[ic];
      for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++)
        if (it->GetCellIdx()==cell)
          etot += it->GetEnergy();
    }
    if( debug ) cout << " Etotal " << etot << endl;
    for (size_t ib=0; ib<xboundary.size(); ib++)
    {
      double el=0.;
      double er=0.;
      double xb = xboundary[ib];
      double xp = xhint - xb;
      for (size_t ic=0; ic<long_list.size(); ic++)
      {
        size_t cell = long_list[ic];
        double xcell = cells[cell].GetX();
        for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++)
          if (it->GetCellIdx()==cell) {
            if (xcell < xb)
              el += it->GetEnergy();
            if (xcell > xb)
              er += it->GetEnergy();
          }
      }
      if( etot > 0. )
      {
        if( id_hint == CalorimeterParticle::POSITRON || id_hint == CalorimeterParticle::ELECTRON )
        {
          h.ProfileX_Electrons->Fill(xp,er/etot);
          h.ProfileX_BigAngle_Energy_Electrons[index_big_ax][index_e]->Fill(xp,er/etot) ;
          h.ProfileX_SmallAngle_Energy_Electrons[index_small_ax][index_e]->Fill(xp,er/etot) ;
        }
        else if( id_hint == CalorimeterParticle::MUON_PLUS || id_hint == CalorimeterParticle::MUON_MINUS )
        {
          h.ProfileX_Muons->Fill(xp,er/etot);
        }
        else if( id_hint == CalorimeterParticle::PION_PLUS || id_hint == CalorimeterParticle::PION_MINUS )
        {
          h.ProfileX_Pions->Fill(xp,er/etot);
        }
      }
    }

    for (size_t ib=0; ib<yboundary.size(); ib++)
    {
      double el=0.;
      double er=0.;
      double yb = yboundary[ib];
      double yp = yhint - yb;
      for (size_t ic=0; ic<long_list.size(); ic++)
      {
        size_t cell = long_list[ic];
        double ycell = cells[cell].GetY();
        for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++)
          if (it->GetCellIdx()==cell) {
            if (ycell < yb)
              el += it->GetEnergy();
            if (ycell > yb)
              er += it->GetEnergy();
          }
      }
      if( etot > 0. )
      {
        if( id_hint == CalorimeterParticle::POSITRON || id_hint == CalorimeterParticle::ELECTRON )
        {
          h.ProfileY_Electrons->Fill(yp,er/etot);
          h.ProfileY_BigAngle_Energy_Electrons[index_big_ay][index_e]->Fill(yp,er/etot) ;
          h.ProfileY_SmallAngle_Energy_Electrons[index_small_ay][index_e]->Fill(yp,er/etot) ;
        }
        else if( id_hint == CalorimeterParticle::MUON_PLUS || id_hint == CalorimeterParticle::MUON_MINUS )
        {
          h.ProfileY_Muons->Fill(yp,er/etot);
        }
        else if( id_hint == CalorimeterParticle::PION_PLUS || id_hint == CalorimeterParticle::PION_MINUS )
        {
          h.ProfileY_Pions->Fill(yp,er/etot);
        }
      }
    }
  }
  if( debug ) cout << " Calorimeter::MakeProfiles " << GetName() << " debuging OK " << endl;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::HistInternalInfo(void)
{
//  bool debug = true;
  bool debug = false;
  if(debug) cout << " HistInternalInfo Started " << GetName() << endl;
  if( !options.fill_histos )
    return;
  if( !options.fill_internal_histos )
    return;
  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
  if( internal_hist_ == NULL )
  {

    if(debug) cout << " Make new HistInternalInfo " << endl;
    internal_hist_ = new InternalHisto(this);
    if( debug ) cout << " Really Make new HistInternalInfo internal_hist_ = " << internal_hist_ << endl;
    if( internal_hist_ == NULL )
    {
      cerr << " ERROR in Calorimeter::HistInternalInfo: No memory for histograms !!" << endl;
      exit(1);
    }
    InternalHisto *ih = internal_hist_;
    if( debug ) cout << " Calorimeter::HistInternalInfo Start Booking " << endl;
    if( debug ) cout << " Calorimeter::HistInternalInfo Yeee Start Booking " << endl;
    if(!gDirectory->cd("/")) assert(false);
    char hist_name[132];
    char dir_name[111];

    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    TDirectory *top_dir = myROOT_utils::TDirectory_checked(dir_name);
    if( debug ) cout << "  Try cd to " << dir_name << endl;
    bool ok = top_dir->cd();
    if( !ok ) assert(false);
    sprintf(dir_name,"%s_test",GetName().c_str());
    TDirectory *stat_dir = myROOT_utils::TDirectory_checked(dir_name);
    if( debug ) cout << " Try cd to " << dir_name << endl;
    ok = stat_dir->cd();
    if( !ok ) assert(false);
    sprintf(dir_name,"%s_internal_info",GetName().c_str());
    TDirectory *internal_info_dir = myROOT_utils::TDirectory_checked(dir_name);
    if( debug ) cout << " Try cd to " << dir_name << endl;
    ok = internal_info_dir->cd();
    if( !ok ) assert(false);

    double emin = 0, emax = options.hist_energy_max;

    if( debug ) cout << " Try to book h.p1_Hi2Ndf " << ih->p1_Hi2Ndf << endl;

    sprintf(hist_name,"%s::  Hi2/Ndf as energy function ",GetName().c_str());
    ih->p1_Hi2Ndf = myROOT_utils::TProfile_checked("Hi2Ndf", hist_name, 50, emin,emax );
    sprintf(hist_name,"%s::  Probability as energy function ",GetName().c_str());
    ih->p1_Hi2Prob = myROOT_utils::TProfile_checked("Hi2Prob", hist_name, 50, emin,emax );
    sprintf(hist_name,"%s::  Hi2 in cell as energy fraction function   ",GetName().c_str());
    ih->p1_Hi2Cells = myROOT_utils::TProfile_checked("Hi2Cells", hist_name, 50,  0, 1 );
    sprintf(hist_name,"%s::  De/SDe as energy fraction function ",GetName().c_str());
    ih->p1_DENorm = myROOT_utils::TProfile_checked("DENorm", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  absDe/SDe as energy fraction function ",GetName().c_str());
    ih->p1_DEabsNorm = myROOT_utils::TProfile_checked("DEabsNorm", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  (Eexp-Ecalc)/Ecalc as energy fraction function ",GetName().c_str());
    ih->p1_DEcalc = myROOT_utils::TProfile_checked("DEcalc", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  Himax as energy fraction function if ngen.eq.nrec ",GetName().c_str());
    ih->p1_HiMaxOK  = myROOT_utils::TProfile_checked("HimaxOK", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  Hi-Himax vs Himax  if ngen.eq.nrec",GetName().c_str());
    ih->h2_HiHiMaxOK = myROOT_utils::TH2D_checked("HiHimaxOK", hist_name, 250, 0., 500., 250, 0., 500.);
    sprintf(hist_name,"%s::  Himax as energy fraction function if ngen.ne.nrec ",GetName().c_str());
    ih->p1_HiMaxBad  = myROOT_utils::TProfile_checked("HimaxBad", hist_name, 50, 0, 1);
    sprintf(hist_name,"%s::  Hi-Himax vs Himax if ngen.ne.nrec ",GetName().c_str());
    ih->h2_HiHiMaxBad = myROOT_utils::TH2D_checked("HiHimaxBad", hist_name, 250, 0., 500., 250, 0., 500.);

//     TProfile *p1_ClusterSizeOfInterest;
//     TProfile *p1_ClusterSize_absEcut;
//     TProfile *p1_ClusterSize_relEcut;
    sprintf(hist_name,"%s::  Cluster size Of potential Interest ",GetName().c_str());
    ih->p1_ClusterSizeOfInterest = myROOT_utils::TProfile_checked("ClusterSizeOfInterest", hist_name, 50, emin,emax );
    sprintf(hist_name,"%s::  Cluster size with absolute Ecut 50 MeV ",GetName().c_str());
    ih->p1_ClusterSize_absEcut = myROOT_utils::TProfile_checked("ClusterSize_absEcut", hist_name, 50, emin,emax );
    sprintf(hist_name,"%s::  Cluster size  with relative Ecut 2 Sigma ",GetName().c_str());
    ih->p1_ClusterSize_relEcut = myROOT_utils::TProfile_checked("ClusterSize_relEcut", hist_name, 50, emin,emax );

    if( debug ) cout << " Calorimeter::HistInternalInfo Booking  OK " << endl;
  }

  dir_save->cd();

  if(debug) cout << " HistInternalInfo OK " << GetName() << endl;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillInternalCorrelations(void)
{
//  bool debug = true;
//  bool debug = false;
  if( !options.fill_internal_correlations )
    return;

  if( internal_correlations_all_info_.size() == 0 )
  {
    stat_total_internal_correlations_ = 0;
    internal_correlations_info_ = new std::vector< StatInfo >[ NCells() ];
    if( internal_correlations_info_ == NULL )
    {
      std::cerr << " ERROR in Calorimeter::FillInternalCorrelations: memory problems?? " << endl;
      exit(1);
    }
    for( size_t i=0; i< NCells(); i++ )
    {
      internal_correlations_all_info_.push_back(StatInfo(0.,0.,0.));
      for( size_t j=0; j< NCells(); j++ )
      {
        internal_correlations_info_[j].push_back(StatInfo(0.,0.,0.));
      }
    }
  }

  stat_total_internal_correlations_++;

  if( internal_correlations_hist_ == NULL )
  {
    internal_correlations_hist_ = new CorrelationsHisto(this);
    if( internal_correlations_hist_ == NULL )
    {
      cerr << " ERROR in Calorimeter::FillInternalCorrelations: No memory for histograms !!" << endl;
      exit(1);
    }
    CorrelationsHisto &h = *internal_correlations_hist_;  // create a short name
    if(!gDirectory->cd("/")) assert(false);
    char hist_name[132];
    char dir_name[111];

    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    TDirectory *top_dir = myROOT_utils::TDirectory_checked(dir_name);
    bool ok = top_dir->cd();
    if( !ok ) assert(false);
    sprintf(dir_name,"Internal_correlations %s",GetName().c_str());
    TDirectory *internal_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok = internal_dir->cd();
    if( !ok ) assert(false);
    float xmin = (*this).GetXmin() , xmax = (*this).GetXmax() ,
          ymin = (*this).GetYmin() , ymax = (*this).GetYmax();
    double cellx = GetMinCellSizeX();
    double celly = GetMinCellSizeY();

    int nbinsx = (int)(((xmax-xmin)*1.0001)/cellx);
    int nbinsy = (int)(((ymax-ymin)*1.0001)/celly);

    h.h2_XYAll = myROOT_utils::TH2D_checked("XYAll"," X Y hits in cells with E gt 0.7 GeV ",nbinsx,xmin,xmax,nbinsy,ymin,ymax);


    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"dXdYforCell%s",GetCellName(i).c_str());
      sprintf(hist_name," (X - Xcorr ) vs (Y - Ycorr) in cell %s ",GetCellName(i).c_str());
      double xc = cells[i].GetX();
      double yc = cells[i].GetY();

      TH2D *h2 = myROOT_utils::TH2D_checked(name,hist_name,nbinsx,xmin-xc,xmax-xc,nbinsy,ymin-yc,ymax-yc);
      h.h2_XdeltaYdeltaCell.push_back(h2);

      sprintf(name,"CorrXYforCell%s",GetCellName(i).c_str());
      sprintf(hist_name," Correlations (X - Xcorr ) vs (Y - Ycorr) in cell %s ",GetCellName(i).c_str());
      h2 = myROOT_utils::TH2D_checked(name,hist_name,nbinsx,xmin-xc,xmax-xc,nbinsy,ymin-yc,ymax-yc);
      h.h2_XYCellCorrelations.push_back(h2);
    }
  }
  CorrelationsHisto &h= *internal_correlations_hist_;  // create a short name

  double emin4correlations = 0.7;

  for(size_t i1=0; i1< signals_.size(); i1++)
  {
    size_t incell0=signals_[i1].GetCellIdx();
    double e0 = signals_[i1].GetEnergy();
    if( e0 < emin4correlations ) continue;
    double x0 = cells[incell0].GetX();
    double y0 = cells[incell0].GetY();
    if(h.h2_XYAll != NULL) h.h2_XYAll->Fill(x0,y0);

    internal_correlations_all_info_[incell0].Add(e0);
    for(size_t i=0; i< signals_.size(); i++)
    {
      size_t incell=signals_[i].GetCellIdx();
      double e = signals_[i].GetEnergy();
      if( e < emin4correlations ) continue;
      double x = cells[incell].GetX();
      double y = cells[incell].GetY();
      if( incell != incell0 )
      {
        if(h.h2_XdeltaYdeltaCell[incell0] != NULL) h.h2_XdeltaYdeltaCell[incell0]->Fill(x-x0,y-y0);
      }
      internal_correlations_info_[incell0][incell].Add(e);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectInternalCorrelations(void) const
{
//  bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::InspectInternalCorrelations for " << GetName() << endl;
  if( !options.fill_internal_correlations )
    return 0;
  if( internal_correlations_hist_ == NULL ) return 0; // Histograms were not initialized and filled
  CorrelationsHisto &h= *internal_correlations_hist_;  // create a short name
  if( !h.h2_XYCellCorrelations.size() == NCells() )
  {
    std::cerr << " Calorimeter::InspectInternalCorrelations Internal ERROR h2_XYCellCorrelations.size() = " <<
                  h.h2_XYCellCorrelations.size()  << std::endl;
  }

//   for( size_t i=0; i< NCells(); i++ )
//   {
//     double entries = h.h2_XdeltaYdeltaCell[i]->GetEntries();
//     TAxis *axis = h.h2_XdeltaYdeltaCell[i]->GetXaxis();
// //    Int_t ncx   = axis->GetNbins();
//     Int_t ncx   = h.h2_XdeltaYdeltaCell[i]->GetNbinsX();
//     Int_t ncy   = h.h2_XdeltaYdeltaCell[i]->GetNbinsY();
// //    Int_t ncx   = axis->GetNbinsX();
// //    Int_t ncy   = axis->GetNbinsY();
//     double summ = 0;
//     Int_t binx,biny;
//     for (binx=ncx+1;binx>=1;binx--)
//       for (biny=ncy+1;biny>=1;biny--)
//       {
//         double content = h.h2_XdeltaYdeltaCell[i]->GetBinContent(binx,biny);
// 	summ += content;
// //        Axis_t ecut = h.h2_XdeltaYdeltaCell[i]->GetBinCenter(bin);
//       }
//     double summ1 = 0;
//     for( size_t j=0; j< NCells(); j++ )
//     {
//       double stat = internal_correlations_info_[i][j].GetEntries();
//       summ1 += stat;
//     }
//     cout << " Cell " << i << " Entries  " << entries << " Summ " << summ << " Summ1 " << summ1 << endl;
//   }
  if( stat_total_internal_correlations_ <= 0 ) return -1; // no entries

  vector < double > prob_noise( NCells() );
  for( size_t i=0; i< NCells(); i++ )
  {
    prob_noise[i] = internal_correlations_all_info_[i].GetEntries()/stat_total_internal_correlations_;
  }

  for( size_t i=0; i< NCells(); i++ )
  {

    double entries = internal_correlations_info_[i][i].GetEntries();
    double x0 = cells[i].GetX();
    double y0 = cells[i].GetY();
    for( size_t j=0; j< NCells(); j++ )
    {
      double stat_noise = prob_noise[j]*entries;
      double sigma_stat_noise = std::sqrt( stat_noise );
      double stat = internal_correlations_info_[i][j].GetEntries() - stat_noise;
      double x = cells[j].GetX();
      double y = cells[j].GetY();
      if( stat > 2.*sigma_stat_noise && i != j )
      {
        h.h2_XYCellCorrelations[i]-> Fill( x-x0, y-y0, stat );
      }
    }
  }

  return 0;

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillExternalCorrelations(void)
{
//  bool debug = true;
//  bool debug = false;
  if( !options.fill_external_correlations )
    return;
  if( external_correlations_all_info_.size() == 0 )
  {
    stat_total_external_correlations_ = 0;
    external_correlations_info_ = new std::vector< StatInfo >[ NCells() ];
    external_correlations_hist_ = new CorrelationsHisto(this);
// Start booking
    CorrelationsHisto &h = *external_correlations_hist_;  // create a short name
    if(!gDirectory->cd("/")) assert(false);
    char hist_name[132];
    char dir_name[111];

    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    TDirectory *top_dir = myROOT_utils::TDirectory_checked(dir_name);
    bool ok = top_dir->cd();
    if( !ok ) assert(false);
    sprintf(dir_name,"External_correlations %s",GetName().c_str());
    TDirectory *internal_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok = internal_dir->cd();
    if( !ok ) assert(false);
    float xmin = (*this).GetXmin() , xmax = (*this).GetXmax() ,
          ymin = (*this).GetYmin() , ymax = (*this).GetYmax();
    double cellx = GetMinCellSizeX();
    double celly = GetMinCellSizeY();

    int nbinsx = (int)(((xmax-xmin)*1.0001)/cellx);
    int nbinsy = (int)(((ymax-ymin)*1.0001)/celly);

    h.h2_XYAll = myROOT_utils::TH2D_checked("XYAll"," X Y hits in cells with E gt 0.7 GeV ",nbinsx,xmin,xmax,nbinsy,ymin,ymax);


    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"dXdYforCell%s",GetCellName(i).c_str());
      sprintf(hist_name," (X - Xcorr ) vs (Y - Ycorr) in cell %s ",GetCellName(i).c_str());
      double xc = cells[i].GetX();
      double yc = cells[i].GetY();

      TH2D *h2 = myROOT_utils::TH2D_checked(name,hist_name,nbinsx,xmin-xc,xmax-xc,nbinsy,ymin-yc,ymax-yc);
      h.h2_XdeltaYdeltaCell.push_back(h2);

      sprintf(name,"CorrXYforCell%s",GetCellName(i).c_str());
      sprintf(hist_name," Correlations (X - Xcorr ) vs (Y - Ycorr) in cell %s ",GetCellName(i).c_str());
      h2 = myROOT_utils::TH2D_checked(name,hist_name,nbinsx,xmin-xc,xmax-xc,nbinsy,ymin-yc,ymax-yc);
      h.h2_XYCellCorrelations.push_back(h2);
    }
// Finish booking

    if( external_correlations_info_ == NULL )
    {
      std::cerr << " ERROR in Calorimeter::FillExternalCorrelations: memory problems?? " << endl;
      exit(1);
    }
    for( size_t i=0; i< NCells(); i++ )
    {
      external_correlations_all_info_.push_back(StatInfo(0.,0.,0.));
      for( size_t j=0; j< NCells(); j++ )
      {
        external_correlations_info_[j].push_back(StatInfo(0.,0.,0.));
      }
    }
  }
// Finish initialisation
  CorrelationsHisto &h= *external_correlations_hist_;  // create a short name

  stat_total_external_correlations_++;

  double emin4correlations = 0.7;

  for(size_t itg=0; itg< hint_particles.size(); itg++)
  {
//    double xg = hint_particles[itg].GetX();
//    double yg = hint_particles[itg].GetY();
//    double zg = hint_particles[itg].GetZ();
//     axg = hint_particles[itg].GetAngleX();
//     ayg = hint_particles[itg].GetAngleY();

    int incell0 = FindCell (hint_particles[itg]);

    if( incell0 < 0 ) continue;
    double e0 = hint_particles[itg].GetE();
    double x0 = cells[incell0].GetX();
    double y0 = cells[incell0].GetY();
//     if( e0 < emin4correlations ) continue;

    external_correlations_all_info_[incell0].Add(e0);
    for(size_t i=0; i<signals_.size(); i++)
    {
      size_t incell=signals_[i].GetCellIdx();
      double e = signals_[i].GetEnergy();
      double x = cells[incell].GetX();
      double y = cells[incell].GetY();
      if( e < emin4correlations ) continue;

      if(h.h2_XdeltaYdeltaCell[incell0] != NULL) h.h2_XdeltaYdeltaCell[incell0]->Fill(x-x0,y-y0);

      external_correlations_info_[incell0][incell].Add(e);
    }

  }


}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectExternalCorrelations(void) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::InspectExternalCorrelations for " << GetName() << endl;
  if( !options.fill_external_correlations )
    return 0;

  if( external_correlations_hist_ == NULL ) return 0; // Histograms were not initialized and filled
  CorrelationsHisto &h= *external_correlations_hist_;  // create a short name
  if( !h.h2_XYCellCorrelations.size() == NCells() )
  {
    std::cerr << " Calorimeter::InspectExternalCorrelations Internal ERROR h2_XYCellCorrelations.size() = " <<
                  h.h2_XYCellCorrelations.size()  << std::endl;
  }

  if( stat_total_external_correlations_ <= 0 ) return -1; // no entries

  vector < double > prob_noise( NCells() );
  for( size_t i=0; i< NCells(); i++ )
  {
    prob_noise[i] = external_correlations_all_info_[i].GetEntries()/stat_total_external_correlations_;
  }

  for( size_t i=0; i< NCells(); i++ )
  {

    double entries = external_correlations_info_[i][i].GetEntries();
    double x0 = cells[i].GetX();
    double y0 = cells[i].GetY();
    for( size_t j=0; j< NCells(); j++ )
    {
      double stat_noise = prob_noise[j]*entries;
      double sigma_stat_noise = std::sqrt( stat_noise );
      double stat = external_correlations_info_[i][j].GetEntries() - stat_noise;
      double x = cells[j].GetX();
      double y = cells[j].GetY();
      if( stat > 2.*sigma_stat_noise && i != j )
      {
        h.h2_XYCellCorrelations[i]-> Fill( x-x0, y-y0, stat );
      }
    }
  }

  return 0;


}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::ShowEvent          (void)
{
//   cout << " Start with visualisation of " << GetName() << endl;
  bool ok=false;
  if( !options.fill_histos )
    return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);
  sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }
//   sprintf(dir_name,"%s_test",GetName().c_str());
//   if( !gDirectory->cd(dir_name) )
//   {
//     cerr << "Can not change directory to " << dir_name << endl;
//     dir_save->cd();
//     return;
//   }

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c) {
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);
  }

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  if (h.c == NULL) {
//   create canvas
//     cout << " create canvas " << endl;
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    h.c->Clear();
    h.c->Divide(3,2);
    h.c->cd(1);
    if(h.h1_nrec != NULL) h.h1_nrec->Draw();
    h.c->cd(2);
    if(h.h1_Erec != NULL) h.h1_Erec->Draw();
    h.c->cd(3);
    if(h.h2_XrYr != NULL) h.h2_XrYr->Draw("BOX");
    h.c->cd(4);
    if(h.h1_ETotal != NULL) h.h1_ETotal->Draw();
    h.c->cd(5);
    if(h.h1_ESumADC != NULL) h.h1_ESumADC->Draw();
    h.c->cd(6);
    if(h.h1_Time != NULL) h.h1_Time->Draw();
  }
  else {
    h.c->Clear();
    h.c->Divide(3,2);
    h.c->cd(1);
    if(h.h1_nrec != NULL) h.h1_nrec->Draw();
    h.c->cd(2);
    if(h.h1_Erec != NULL) h.h1_Erec->Draw();
    h.c->cd(3);
    if(h.h2_XrYr != NULL) h.h2_XrYr->Draw("BOX");
    h.c->cd(4);
    if(h.h1_ETotal != NULL) h.h1_ETotal->Draw();
    h.c->cd(5);
    if(h.h1_ESumADC != NULL) h.h1_ESumADC->Draw();
    h.c->cd(6);
    if(h.h1_Time != NULL) h.h1_Time->Draw();
    h.c->Update ();
    h.c->Show ();
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillHistoRAW(void)
{
//   bool debug = true;
  bool debug = false;

  bool ok=false;
  if( !options.fill_histos )
    return;
  if( !options.fill_raw_histos )
    return;

  if( debug ) cout << " Calorimeter::FillHistoRAW " << GetName() << " debug " << endl;

  if( gDirectory==NULL )
  {
// #warning Some bug in ReconstructionTest:: Check(gDirectory==NULL) does not help if ROOT file was not opend
    cerr << "Calorimeter::ReconstructionTest():  ROOT file was not opend.\n"
         << "                                    I do not want to work.\n";
    return;
  }

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
// create directory structure for RAW info
  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if( h.root_dir==NULL )
  {
    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
  }
  ok=(h.root_dir->cd());
  assert(ok);

  if(!h.raw_histos_booked)
  {
    char dir_name[111];
// Create subdirectory  for the cells
    sprintf(dir_name,"%s_RAW_Cells",GetName().c_str());
    h.cells_raw_dir = myROOT_utils::TDirectory_checked(dir_name);

    ok=(h.cells_raw_dir->cd());
    assert(ok);

    char hist_name[132];
    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"RAW%s",GetCellName(i).c_str());
      sprintf(hist_name," RAW data in cell %s ",GetCellName(i).c_str());
      h.h1_RAWCell[i] = myROOT_utils::TH1D_checked(name,hist_name,4000,0.,4000.);

    }
    h.raw_histos_booked=true;
  }

  if( debug ) cout << " signals_.size() = " << signals_.size() << endl;

  for( vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++ )
  {
    if( it->GetCellIdx() >= NCells() )
    {
      cerr << " Bad raw data in " << GetName() << " icell = " << it->GetCellIdx() << endl;
      continue;
    }
    if( h.h1_RAWCell[it->GetCellIdx()] != NULL ) h.h1_RAWCell[it->GetCellIdx()]->Fill(it->GetAmplitude());
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillHistoLED(void)
{
  bool ok=false;
  if( !options.fill_histos )
    return;

  if( gDirectory==NULL )
  {
// #warning Some bug in ReconstructionTest:: Check(gDirectory==NULL) does not help if ROOT file was not opend
    cerr << "Calorimeter::FillHistoLED():  ROOT file was not opend.\n"
         << "                                    I do not want to work.\n";
    return;
  }

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
//  TDirectory *dir_check = NULL;
// create directory structure for LEDs
  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if( h.root_dir==NULL )
  {
    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
  }
  ok=(h.root_dir->cd());
  assert(ok);

  if(!h.led_histos_booked)
  {
    char dir_name[111];
// Create subdirectory  for the cells
    sprintf(dir_name,"%s_LED_Cells",GetName().c_str());
    h.cells_led_dir = myROOT_utils::TDirectory_checked(dir_name);

    ok=(h.cells_led_dir->cd());
    assert(ok);

    char hist_name[132];
    int nbins = options.led_histo_size;
    double max = options.hist_led_max;
    for( size_t i=0; i< NCells(); i++ )
    {
//kill      strstream name;
      char name[132];
//kill      name.form("LED%s",GetCellName(i).c_str());
      sprintf(name,"LED%s",GetCellName(i).c_str());
//      cout << " Now in LED Cells for BOOKING " << GetName() << " i=" << i << " name " <<  name.str() <<   endl;
      sprintf(hist_name," LED in cell %s ",GetCellName(i).c_str());
//kill      h.h1_LEDCell[i] = myROOT_utils::TH1D_checked(name.str(),hist_name,500,0.,5000.);
//      h.h1_LEDCell[i] = myROOT_utils::TH1D_checked(name,hist_name,500,0.,5000.);
      h.h1_LEDCell[i] = myROOT_utils::TH1D_checked(name,hist_name,nbins,0.,max);

    }
    h.led_histos_booked=true;
  }


  for( vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++ )
  {
    h.h1_LEDCell[it->GetCellIdx()]->Fill(it->GetAmplitude());
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillHistoPED(void)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " DEBUG  FillHistoPED for " << GetName() << endl;

  bool ok=false;
  if( !options.fill_histos )
    return;

  if( gDirectory==NULL )
  {
    cerr << "Calorimeter::FillHistoPED():  ROOT file was not opend.\n"
         << "                                    I do not want to work.\n";
    return;
  }

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
//  TDirectory *dir_check = NULL;
// create directory structure for LEDs
  if( debug ) cout << " calo_hist " << calo_hist << endl;
  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if( h.root_dir==NULL )
  {
    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
    if( debug ) cout << " Make top dir OK " << endl;
  }
  ok=(h.root_dir->cd());
  assert(ok);

  if(!h.ped_histos_booked)
  {
    if( debug ) cout << " BOOK PED HISTOh. h1_PEDCell.size() " << h.h1_PEDCell.size() << endl;
    char dir_name[111];
// Create subdirectory  for the cells
    sprintf(dir_name,"%s_PED_Cells",GetName().c_str());
    h.cells_ped_dir = myROOT_utils::TDirectory_checked(dir_name);

    ok=(h.cells_ped_dir->cd());
    assert(ok);

    char hist_name[132];
    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"PED%s",GetCellName(i).c_str());
      sprintf(hist_name," PED in cell %s ",GetCellName(i).c_str());
      h.h1_PEDCell[i] = myROOT_utils::TH1D_checked(name,hist_name,1000,0.,1000.);

    }
    h.ped_histos_booked=true;
    if( debug ) cout << " BOOK PED HISTO OK !!!!" << endl;
  }

  if( debug ) cout << " Fill PED HISTO " << endl;

  for( vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++ )
  {
    h.h1_PEDCell[it->GetCellIdx()]->Fill(it->GetAmplitude());
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::DrawCellHisto          (int icell)
{
  bool ok=false;
//   cout << " Start with visualisation of " << GetName() << endl;
  if( icell <0 || icell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << endl;
     return;
  }

  if( !options.fill_histos )
    return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);
  sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }
  sprintf(dir_name,"%s_test",GetName().c_str());

//   if( !gDirectory->cd(dir_name) )
//   {
//     cerr << "Can not change directory to " << dir_name << endl;
//     dir_save->cd();
//     return;
//   }

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c) {
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);
  }

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  if (h.c == NULL) {
//   create canvas
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    h.c->Clear();
    h.c->Divide(2,2);
    h.c->cd(1);
    if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
    h.c->cd(2);
    if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    h.c->cd(3);
    if(h.h1_ADCmCell[icell] != NULL) h.h1_ADCmCell[icell]->Draw();
    h.c->cd(4);
    if(h.h1_RAWCell[icell] != NULL) h.h1_RAWCell[icell]->Draw();
//    if(options.enable_time_calibration) if(h.h1_TimeCell[icell] != NULL) h.h1_TimeCell[icell]->Draw();
//    if(h.h1_ETotal != NULL) h.h1_ETotal->Draw();
  }
  else {
    h.c->Clear();
    h.c->Divide(2,2);
    h.c->cd(1);
    if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
    h.c->cd(2);
    if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    h.c->cd(3);
    if(h.h1_ADCmCell[icell] != NULL) h.h1_ADCmCell[icell]->Draw();
    h.c->cd(4);
    if(h.h1_RAWCell[icell] != NULL) h.h1_RAWCell[icell]->Draw();
//    if(options.enable_time_calibration) if(h.h1_TimeCell[icell] != NULL) h.h1_TimeCell[icell]->Draw();
//    if(h.h1_ETotal != NULL) h.h1_ETotal->Draw();
    h.c->Update ();
    h.c->Show ();
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::DrawMonitorCellHisto          (int icell)
{
  bool ok=false;
  cout << " Start with visualisation of MONITORING " << GetName() << endl;
  if( icell <0 || icell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << endl;
     return;
  }

  if( !options.fill_histos )
    return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);
  sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }
  sprintf(dir_name,"%s_test",GetName().c_str());

//   if( !gDirectory->cd(dir_name) )
//   {
//     cerr << "Can not change directory to " << dir_name << endl;
//     dir_save->cd();
//     return;
//   }

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c) {
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);
  }

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  if (h.c == NULL) {
//   create canvas
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    h.c->Clear();
    h.c->Divide(2);
    h.c->cd(1);
    if(cells_store[CALIB][MONITOR] != NULL ) cells_store[CALIB][MONITOR][icell].DrawHisto();
    h.c->cd(2);
    if(cells_store[LED][OLD] != NULL ) cells_store[LED][OLD][icell].DrawHisto();
  }
  else {
    h.c->Clear();
    h.c->Divide(2);
    h.c->cd(1);
    if(cells_store[CALIB][MONITOR] != NULL ) cells_store[CALIB][MONITOR][icell].DrawHisto();
    h.c->cd(2);
    if(cells_store[LED][OLD] != NULL ) cells_store[LED][OLD][icell].DrawHisto();
    h.c->Update ();
    h.c->Show ();
  }

  dir_save->cd();
}

/////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::DrawMonitorCellHistoInSpill          (int icell)
{
  bool ok=false;
  cout << " Start with visualisation of InSpill MONITORING " << GetName() << endl;
  if( icell <0 || icell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << endl;
     return;
  }

  if( p1_led_in_spills_.size()== 0 )
  {
    cout << " NO MonitorInSpill Histos  " << endl;
    return;
  }


//   char dir_name[111];
//   TDirectory *dir_save = gDirectory;
//   ok=(gDirectory->cd("/"));
//   assert(ok);
//   sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
//   if( !gDirectory->cd(dir_name) )
//   {
//     cerr << "Can not change directory to " << dir_name << endl;
//     dir_save->cd();
//     return;
//   }
//   sprintf(dir_name,"MonitorLED");

//   if( !gDirectory->cd(dir_name) )
//   {
//     cerr << "Can not change directory to " << dir_name << endl;
//     dir_save->cd();
//     return;
//   }

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != c_led_in_spills_) {
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);
  }

  if (curcan == NULL)
  {
    // Canvas was killed probably
    if (c_led_in_spills_ != NULL) c_led_in_spills_ = NULL;
  }

// end of check


  if (c_led_in_spills_ == NULL) {
//   create canvas
    char canvas_name[132];
    sprintf(canvas_name,"MonitorLED_%s_c",GetName().c_str());
    c_led_in_spills_ = new TCanvas(canvas_name,canvas_name);
    c_led_in_spills_->Clear();
    c_led_in_spills_->Divide(2);
    c_led_in_spills_->cd(1);
    if(p1_led_in_spills_[icell] != NULL ) p1_led_in_spills_[icell]->Draw();
    c_led_in_spills_->cd(2);
    if(cells_store[LED][OLD] != NULL ) cells_store[LED][OLD][icell].DrawHisto();
  }
  else {
    c_led_in_spills_->Clear();
    c_led_in_spills_->Divide(2);
    c_led_in_spills_->cd(1);
    if(p1_led_in_spills_[icell] != NULL ) p1_led_in_spills_[icell]->Draw();
    c_led_in_spills_->cd(2);
    if(cells_store[LED][OLD] != NULL ) cells_store[LED][OLD][icell].DrawHisto();
    c_led_in_spills_->Update ();
    c_led_in_spills_->Show ();
  }

//   dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::DrawCellsHistoLED       (int from_cell,int how_much,int step)
{
  bool ok=false;
  int to_cell = from_cell + (how_much -1)*step;
  if( from_cell <0 || from_cell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << from_cell << endl;
     return;
  }

  if( to_cell <0 || to_cell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << to_cell << endl;
     return;
  }

  if( !options.fill_histos )
    return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);
  sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }
  sprintf(dir_name,"%s_test",GetName().c_str());

  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }

  int nhist = how_much;
  int nx = (int)std::sqrt((double)nhist);
  int ny=nx;
  while( nhist-nx*ny > 0 ) nx += 1;

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c) {
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);
  }

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  if (h.c == NULL) {
//   create canvas
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    h.c->Clear();
    h.c->Divide(nx,ny);
    int ic = 1;
    for( int icell=from_cell; icell < to_cell+1; icell +=step )
    {
      if(icell>= (int)NCells()) break;
      h.c->cd(ic);
      ic++;
      if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    }
  }
  else {
    h.c->Clear();
    h.c->Divide(nx,ny);
    int ic = 1;
    for( int icell=from_cell; icell < to_cell+1; icell +=step )
    {
      if(icell>= (int)NCells()) break;
      h.c->cd(ic);
      ic++;
//      if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
      if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    }
    h.c->Update ();
    h.c->Show ();
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::DrawCellsHistoGAM       (int from_cell,int how_much,int step)
{

  bool ok=false;
  int to_cell = from_cell + (how_much -1)*step;
  if( from_cell <0 || from_cell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << from_cell << endl;
     return;
  }

  if( to_cell <0 || to_cell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << to_cell << endl;
     return;
  }

  if( !options.fill_histos )
    return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);
  sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }
//   sprintf(dir_name,"%s_test",GetName().c_str());
//
//   if( !gDirectory->cd(dir_name) )
//   {
//     cerr << "Can not change directory to " << dir_name << endl;
//     dir_save->cd();
//     return;
//   }

  int nhist = how_much;
  int nx = (int)std::sqrt((double)nhist);
  int ny=nx;
  while( nhist-nx*ny > 0 ) nx += 1;

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c) {
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);
  }

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  if (h.c == NULL) {
//   create canvas
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    h.c->Clear();
    h.c->Divide(nx,ny);
    int ic = 1;
    for( int icell=from_cell; icell!=to_cell+1;  icell +=step )
    {
      if(icell>= (int)NCells()) break;
      h.c->cd(ic);
      ic++;
      if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
//      if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    }
  }
  else {
    h.c->Clear();
    h.c->Divide(nx,ny);
    int ic = 1;
    for( int icell=from_cell; icell < to_cell+1; icell +=step )
    {
      if(icell>= (int)NCells()) break;
      h.c->cd(ic);
      ic++;
      if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
//      if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    }
    h.c->Update ();
    h.c->Show ();
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::FitCellsHistoGAM       (int from_cell,int how_much,int step)
{
  bool ok=false;

  int to_cell = from_cell + (how_much -1)*step;
  if( from_cell <0 || from_cell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << from_cell << endl;
     return;
  }

  if( to_cell <0 || to_cell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << to_cell << endl;
     return;
  }

  if( !options.fill_histos )
    return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);
  sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }
//   sprintf(dir_name,"%s_test",GetName().c_str());
//
//   if( !gDirectory->cd(dir_name) )
//   {
//     cerr << "Can not change directory to " << dir_name << endl;
//     dir_save->cd();
//     return;
//   }

  int nhist = how_much;
  int nx = (int)std::sqrt((double)nhist);
  int ny=nx;
  while( nhist-nx*ny > 0 ) nx += 1;

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c) {
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);
  }

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  if (h.c == NULL) {
//   create canvas
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    h.c->Clear();
    h.c->Divide(nx,ny);
    int ic = 1;
    for( int icell=from_cell; icell < to_cell+1; icell +=step )
    {
      if(icell>= (int)NCells()) break;
      h.c->cd(ic);
      ic++;
      if(h.h1_ErCell[icell] != NULL)
      {
//         double vmin = h.h1_ErCell[icell]->GetXaxis()->GetXmin();
//         double vmax = h.h1_ErCell[icell]->GetXaxis()->GetXmax();
//        TF1 *func = new TF1("fit","gaus",vmin,vmax);
        h.h1_ErCell[icell]->Fit("fit","R0");
        h.h1_ErCell[icell]->Draw();
      }
//      if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
//      if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    }
  }
  else {
    h.c->Clear();
    h.c->Divide(nx,ny);
    int ic = 1;
    for( int icell=from_cell; icell < to_cell+1; icell +=step )
    {
      if(icell>= (int)NCells()) break;
      h.c->cd(ic);
      ic++;
      if(h.h1_ErCell[icell] != NULL)
      {
//         double vmin = h.h1_ErCell[icell]->GetXaxis()->GetXmin();
//         double vmax = h.h1_ErCell[icell]->GetXaxis()->GetXmax();
//        TF1 *func = new TF1("fit","gaus",vmin,vmax);
        h.h1_ErCell[icell]->Fit("fit","R0");
        h.h1_ErCell[icell]->Draw();
      }
//      if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
//      if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    }
    h.c->Update ();
    h.c->Show ();
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::FitCellsHistoLED       (int from_cell,int how_much,int step)
{
  bool ok=false;

  int to_cell = from_cell + (how_much -1)*step;
  if( from_cell <0 || from_cell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << from_cell << endl;
     return;
  }

  if( to_cell <0 || to_cell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << to_cell << endl;
     return;
  }

  if( !options.fill_histos )
    return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);
  sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }
//   sprintf(dir_name,"%s_test",GetName().c_str());
//
//   if( !gDirectory->cd(dir_name) )
//   {
//     cerr << "Can not change directory to " << dir_name << endl;
//     dir_save->cd();
//     return;
//   }

  int nhist = how_much;
  int nx = (int)std::sqrt((double)nhist);
  int ny=nx;
  while( nhist-nx*ny > 0 ) nx += 1;

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c) {
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);
  }

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  if (h.c == NULL) {
//   create canvas
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    h.c->Clear();
    h.c->Divide(nx,ny);
    int ic = 1;
    for( int icell=from_cell; icell < to_cell+1; icell +=step )
    {
      if(icell>= (int)NCells()) break;
      h.c->cd(ic);
      ic++;
      if(h.h1_LEDCell[icell] != NULL)
      {
//         double vmin = h.h1_LEDCell[icell]->GetXaxis()->GetXmin();
//         double vmax = h.h1_LEDCell[icell]->GetXaxis()->GetXmax();
//        TF1 *func = new TF1("fit","gaus",vmin,vmax);
        h.h1_LEDCell[icell]->Fit("fit","R0");
        h.h1_LEDCell[icell]->Draw();
      }
//      if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
//      if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    }
  }
  else {
    h.c->Clear();
    h.c->Divide(nx,ny);
    int ic = 1;
    for( int icell=from_cell; icell < to_cell+1; icell +=step )
    {
      if(icell>= (int)NCells()) break;
      h.c->cd(ic);
      ic++;
      if(h.h1_LEDCell[icell] != NULL)
      {
        double vmin = h.h1_LEDCell[icell]->GetXaxis()->GetXmin();
        double vmax = h.h1_LEDCell[icell]->GetXaxis()->GetXmax();
        TF1 *func = new TF1("fit","gaus",vmin,vmax);
        h.h1_LEDCell[icell]->Fit("fit","R0");
        h.h1_LEDCell[icell]->Draw();

        cout << " Fit LED in cell " << icell << " Mean " <<func->GetParameter(1) <<
                 " Sigma " << func->GetParameter(2) << " Stat " << func->GetParameter(0) << endl;
      }
//      if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
//      if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    }
    h.c->Update ();
    h.c->Show ();
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void    Calorimeter::DrawFitCellsHisto (CellInfoType info_type,const char *opt,
                                               int from_cell,int how_much,int step,
                                               double fmin, double fmax)
{
  bool ok=false;

  bool use_cell_name_to_finde_histo = false;
  string opts;
//   cout << " Calorimeter::DrawCellsHistoXY " << info_type << endl;
  if(step<=0) step=1;
  if(how_much<=0) how_much=1;
  int to_cell = from_cell + step*(how_much-1);

  int nhist = how_much;
  int nx = (int)std::sqrt((double)nhist);
  int ny=nx;
  while( nhist-nx*ny > 0 ) nx += 1;

//   if( !options.fill_histos )
//     return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if (info_type == CALIB && h.h1_CalibCell.size() != NCells() ) return;
  if (info_type == LED && h.h1_LEDCell.size() != NCells() ) return;
  if (info_type == PED && h.h1_PEDCell.size() != NCells() ) return;
  if (info_type == TIME && h.h1_TimeCell.size() != NCells() ) return;
  if (info_type == NOISE && h.h1_NOISECell.size() != NCells() ) return;

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);

  if (info_type == CALIB) sprintf(dir_name,"Calorimeter_%s/%s_test/%s_Cells",GetName().c_str(),GetName().c_str(),GetName().c_str());
  else if(info_type == LED) sprintf(dir_name,"Calorimeter_%s/%s_LED_Cells",GetName().c_str(),GetName().c_str());
  else if(info_type == PED) sprintf(dir_name,"Calorimeter_%s/%s_PED_Cells",GetName().c_str(),GetName().c_str());
  else if(info_type == TIME) sprintf(dir_name,"Calorimeter_%s/%s_CellsTime",GetName().c_str(),GetName().c_str());
  else if(info_type == NOISE) sprintf(dir_name,"Calorimeter_%s/%s_CellsNoise",GetName().c_str(),GetName().c_str());
  else return;

  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c)
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  bool canvas_new=false;
  if (h.c == NULL) //   create new canvas
  {
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    canvas_new=true;
  }

  if(0!=strcmp(opt,"FIT_ON_OFF" ) )   // In case Change Fit validity dont touch canvas
  {
    h.c->Clear();
    h.c->Divide(nx,ny);
  }
  int ic = 1;
  int fit_status=0;
  bool fit_valid=true;
    for( int icell=from_cell; icell < to_cell+1; icell +=step )// cycle over cells
    {
      if(0!=strcmp(opt,"FIT_ON_OFF" ) )
      {
        h.c->cd(ic);
        ic++;
      }
      if(icell < 0 || icell >= (int)NCells() )
      {
        cerr << " FATAL!!  DrawFitCellsHistoXY INTERNAL ERROR !!! Wrong cell " << icell << endl;
        exit(1);
      }
      TH1D *h1 = NULL;
      if( use_cell_name_to_finde_histo )
      {
        char cell_name[132];
// option??        if      (info_type == CALIB) sprintf(cell_name,"Cell_%d",icell);
        if      (info_type == CALIB) sprintf(cell_name,"Calib_%d",icell);
        else if (info_type == LED)   sprintf(cell_name,"LED_%d",icell);
        else if (info_type == PED)   sprintf(cell_name,"PED_%d",icell);
        else if (info_type == TIME)  sprintf(cell_name,"Time_%d",icell);
        else if (info_type == NOISE)  sprintf(cell_name,"Noise_%d",icell);
        else    return;
        h1 = (TH1D*) gDirectory->Get(cell_name);
      }
      else
      {
// option??        if      (info_type == CALIB) h1=h.h1_ErCell[icell];
        if      (info_type == CALIB) h1=h.h1_CalibCell[icell];
        else if (info_type == LED)   h1=h.h1_LEDCell[icell];
        else if (info_type == PED)   h1=h.h1_PEDCell[icell];
        else if (info_type == TIME)  h1=h.h1_TimeCell[icell];
        else if (info_type == NOISE)  h1=h.h1_NOISECell[icell];
        else    return;
      }


      if(h1 == NULL) continue;

      if( 0==strcmp(opt,"DRAW") )          // Draw
      {
        if((fmax-fmin)> 0.000001)
        {
          h1->GetXaxis()->SetRange((int)fmin ,(int)fmax);
        }
        h1->Draw();
      }
      else if(  0==strcmp(opt,"FIT" ) )      // Fit
      {
        double vmin = h1->GetXaxis()->GetXmin();
        double vmax = h1->GetXaxis()->GetXmax();
        TF1 *func = new TF1("fit","gaus",vmin,vmax);
//         TF1 *func = new TF1("fit","gaus");

//        h1->Fit("fit","R0");

        h1->Fit("fit","LR0");

//         cout << " Fit LED in cell " << icell << " Mean " <<func->GetParameter(1) <<
//                " Sigma " << func->GetParameter(2) << " Stat " << func->GetParameter(0) << endl;
        fit_info[info_type][icell].fit_parameters.clear();
        fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(func->GetParameter(0),0));
        fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(func->GetParameter(1),0));
        fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(func->GetParameter(2),0));
//         h1->GetXaxis()->SetRange(vmin,vmax);

// //          h1->Integral(bin1,bin2);
// //          c1->SetLogy();
// //          c1->SetLogy(0);
        const Int_t kNotDraw = 1<<9;
//         h1->GetXaxis()->SetRange(binmin,binmax);
        if(h1->GetFunction("fit") != NULL) h1->GetFunction("fit")->ResetBit(kNotDraw);
        h1->Draw();
      }
      else if(0==strcmp(opt,"FIT_ON_OFF" ) )   // Change Fit validity
      {
        if( icell < 0 || icell >= (int)NCells() ) continue;
        if( fit_status == 0 ) // Check enable or disable fit results
        {
          if( fit_info[info_type][icell].fit_ok==true )
          {
            fit_status=-1;
            fit_valid=false;
          }
          else
          {
            fit_status= 1;
            fit_valid=true;
          }
        }

        if( fit_valid ) // Enable Fit
        {
          fit_info[info_type][icell].fit_ok=true;
          if( fit_info[info_type][icell].fit_parameters.empty() )
          {
             cerr << " ERROR !!! :: FIT CAN'T BE VALID WITHOUT FIT PARAMETERS: " <<
                                                 GetName() << " cell# " << icell << endl;
             fit_info[info_type][icell].fit_ok=false;
          }
        }
        else
        {
          fit_info[info_type][icell].fit_ok=false;
        }

	if(info_type == CALIB )
	{
          if( CellIsBad( icell, OLD ) )
	  {
	    if( fit_info[info_type][icell].fit_ok == true )
	    {
	      cout << " WARNING!! Disable fit result for the BAD cell " << GetCellName(icell ) <<" in " << GetName() << endl;
	      fit_info[info_type][icell].fit_ok=false;
	    }
	  }
	}

        continue;
      }
      else
      {
         cerr << " ERROR !!! :: " << GetName() <<
                 " DrawFitCellsHistoXY  UNKNOUN option " << opt << endl;
      }
    }                                         // end of cycle over cells

  if(!canvas_new)
  {
    if(0!=strcmp(opt,"FIT_ON_OFF" ) )   // In case Change Fit validity dont touch canvas
    {
      h.c->Update ();
      h.c->Show ();
    }
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::DrawFitCellsHistoXY       (CellInfoType info_type,const char *opt,
                                             int xmin,int xhow_much,int xstep,
                                             int ymin,int yhow_much,int ystep,
                                             double fmin, double fmax)
{
  bool debug = false;
//   bool debug = true;
  if( !XYRegularGrid() )
  {
    cerr << " ERROR!! call DrawFitCellsHistoXY for " << GetName() << " But XY not initialized " << endl;
    return;
  }
  bool ok=false;
  string opts;
  if( debug )  cout << " Calorimeter::DrawCellsHistoXY " << GetName() <<
                       " info " << info_type << " CalorimeterHist " << calo_hist << endl;
  if(xstep<=0) xstep=1;
  if(xhow_much<=0) xhow_much=1;
  int xmax = xmin + xstep*(xhow_much-1);

  if(ystep<=0) ystep=1;
  if(yhow_much<=0) yhow_much=1;
  int ymax = ymin + ystep*(yhow_much-1);

//   if( !options.fill_histos )
//     return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if( debug )
  {
    if (info_type == CALIB ) cout << " Calib size " << h.h1_CalibCell.size() << endl;
    if (info_type == LED )  cout << " LED size " << h.h1_LEDCell.size()  << endl;
    if (info_type == PED )  cout << " PED size " << h.h1_PEDCell.size()  << endl;
    if (info_type == TIME )  cout << " Time size " << h.h1_TimeCell.size()  << endl;
    if (info_type == NOISE ) cout << " NOISE size " <<  h.h1_NOISECell.size()  << endl;
  }
  if (info_type == CALIB && h.h1_CalibCell.size() != NCells() ) return;
  if (info_type == LED && h.h1_LEDCell.size() != NCells() ) return;
  if (info_type == PED && h.h1_PEDCell.size() != NCells() ) return;
  if (info_type == TIME && h.h1_TimeCell.size() != NCells() ) return;
  if (info_type == NOISE && h.h1_NOISECell.size() != NCells() ) return;

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);

  bool valid_info_type = true;
  if (info_type == CALIB) sprintf(dir_name,"Calorimeter_%s/CellsCalibration_%s",GetName().c_str(),GetName().c_str());
  else if(info_type == LED) sprintf(dir_name,"Calorimeter_%s/%s_LED_Cells",GetName().c_str(),GetName().c_str());
  else if(info_type == PED) sprintf(dir_name,"Calorimeter_%s/%s_PED_Cells",GetName().c_str(),GetName().c_str());
  else if(info_type == TIME) sprintf(dir_name,"Calorimeter_%s/%s_CellsTime",GetName().c_str(),GetName().c_str());
  else if(info_type == NOISE) sprintf(dir_name,"Calorimeter_%s/%s_CellsNoise",GetName().c_str(),GetName().c_str());
  else {
    sprintf(dir_name,"Calorimeter_%s/%s_test/%s_Cells",GetName().c_str(),GetName().c_str(),GetName().c_str());
    valid_info_type = false;
  }

  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c)
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  bool canvas_new=false;
  if (h.c == NULL) //   create new canvas
  {
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    canvas_new=true;
  }

  if(0!=strcmp(opt,"FIT_ON_OFF" ) )   // In case Change Fit validity dont touch canvas
  {
    h.c->Clear();
    h.c->Divide(xhow_much,yhow_much);
  }
  int ic = 1;
  int fit_status=0;
  bool fit_valid=true;
  for( int iy=ymax; iy >= ymin;  iy -=ystep )  // cycle over cells
    for( int ix=xmin; ix <= xmax;  ix +=xstep )
    {
      int icell=GetCellOfColumnRow(ix,iy);
      if(0!=strcmp(opt,"FIT_ON_OFF" ) )
      {
        h.c->cd(ic);
        ic++;
      }
      if( icell < 0 ) continue;

      char cell_name[132];
//       if      (info_type == CALIB) sprintf(cell_name,"CalibX_%d_Y_%d",ix,iy);
//       else if (info_type == LED)   sprintf(cell_name,"LEDX_%d_Y_%d",ix,iy);
//       else if (info_type == PED)   sprintf(cell_name,"PEDX_%d_Y_%d",ix,iy);
//       else if (info_type == TIME)  sprintf(cell_name,"TimeX_%d_Y_%d",ix,iy);
//       else if (info_type == NOISE)  sprintf(cell_name,"NOISEX_%d_Y_%d",ix,iy);
//       else    sprintf(cell_name,"CellX_%d_Y_%d",ix,iy);
      if      (info_type == CALIB) sprintf(cell_name,"Calib%s",GetCellName(icell).c_str());
      else if (info_type == LED)   sprintf(cell_name,"LED%s",GetCellName(icell).c_str());
      else if (info_type == PED)   sprintf(cell_name,"PED%s",GetCellName(icell).c_str());
      else if (info_type == TIME)  sprintf(cell_name,"Time%s",GetCellName(icell).c_str());
      else if (info_type == NOISE)  sprintf(cell_name,"NOISE%s",GetCellName(icell).c_str());
      else    sprintf(cell_name,"Cell%s",GetCellName(icell).c_str());
      TH1D *h1 = (TH1D*) gDirectory->Get(cell_name);
      if(h1 == NULL) continue;
      if(icell < 0 || icell >= (int)NCells() )
      {
        cerr << " FATAL!!  DrawFitCellsHistoXY INTERNAL ERROR !!! Wrong cell " << icell << endl;
        exit(1);
      }
      double scale = 0;
      double maximum = h1->GetXaxis()->GetXmax();
      double minimum = h1->GetXaxis()->GetXmin();
//      cout << " Maximum" << maximum  << " Minimum " << minimum  << " Nbins " << h1->GetNbinsX() << endl;
      if( h1->GetNbinsX() > 0) scale = (maximum-minimum)/h1->GetNbinsX();

      if( 0==strcmp(opt,"DRAW") )          // Draw
      {
        if((fmax-fmin)> 0.000001 && info_type != TIME)
        {
          h1->GetXaxis()->SetRange( h1->GetXaxis()->FindBin(fmin) , h1->GetXaxis()->FindBin(fmax) );
        }
        h1->Draw();
      }
      else if(  0==strcmp(opt,"FIT" ) )      // Fit
      {
        double vmin = minimum;
        double vmax = maximum;
//        cout << " vmin " << vmin << " vmax " << vmax  << endl;
//         if(fmax-fmin> 0.001)
//         {
//           vmin=fmin;
//           vmax=fmax;
//         }

        TF1 *func = new TF1("fit","gaus",vmin,vmax);
//         TF1 *func = new TF1("fit","gaus");
        h1->Fit("fit","R0");
//         cout << " Fit LED in cell " << icell << " Mean " <<func->GetParameter(1) <<
//                " Sigma " << func->GetParameter(2) << " Stat " << func->GetParameter(0) << endl;
        double amplitude = func->GetParameter(0);
        double mean      = func->GetParameter(1);
        double sigma     = func->GetParameter(2);
        double integral  = amplitude*M_PI*sigma/scale;
        if( valid_info_type )
        {
          fit_info[info_type][icell].fit_parameters.clear();
          fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(integral,0));
          fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(mean,0));
          fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(sigma,0));
        }
//         h1->GetXaxis()->SetRange(vmin,vmax);

// //          h1->Integral(bin1,bin2);
// //          c1->SetLogy();
// //          c1->SetLogy(0);
        const Int_t kNotDraw = 1<<9;
//         h1->GetXaxis()->SetRange(binmin,binmax);
        if(h1->GetFunction("fit") != NULL) h1->GetFunction("fit")->ResetBit(kNotDraw);
        h1->Draw();
      }
      else if(0==strcmp(opt,"FIT_ON_OFF" ) && valid_info_type )   // Change Fit validity
      {
        cout << " Change of fit validity by FIT_ON_OFF button cell " << icell <<" is Bad  = " << CellIsBad(icell,OLD) << endl;
        if( icell < 0 || icell >= (int)NCells() ) continue;
        if( fit_status == 0 ) // Check enable or disable fit results
        {
          if( fit_info[info_type][icell].fit_ok==true )
          {
            fit_status=-1;
            fit_valid=false;
          }
          else
          {
            fit_status= 1;
            fit_valid=true;
          }
        }

        if( fit_valid ) // Enable Fit
        {
          fit_info[info_type][icell].fit_ok=true;
          if( fit_info[info_type][icell].fit_parameters.empty() )
          {
             cerr << " ERROR !!! :: FIT CAN'T BE VALID WITHOUT FIT PARAMETERS:" <<
                         GetName() << " X=" << ix << " Y=" << iy << endl;
             fit_info[info_type][icell].fit_ok=false;
          }
        }
        else
        {
          fit_info[info_type][icell].fit_ok=false;
        }

	if(info_type == CALIB )
	{
          if( CellIsBad( icell, OLD ) )
	  {
	    if( fit_info[info_type][icell].fit_ok == true )
	    {
	      cout << " WARNING!! Disable fit result for the BAD cell " << GetCellName(icell ) <<" in " << GetName() << endl;
	      fit_info[info_type][icell].fit_ok=false;
	    }
	  }
	}

        continue;
      }
      else
      {
         cerr << " ERROR !!! :: " << GetName() <<
                 " DrawFitCellsHistoXY  UNKNOUN option " << opt << endl;
      }
    }                                         // end of cycle over cells

  if(!canvas_new)
  {
    if(0!=strcmp(opt,"FIT_ON_OFF" ) )   // In case Change Fit validity dont touch canvas
    {
      h.c->Update ();
      h.c->Show ();
    }
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::DrawFitCellsHistoXYSubset  ( int sub_set, CellInfoType info_type,const char *opt,
                                             int xmin,int xhow_much,int xstep,
                                             int ymin,int yhow_much,int ystep,
                                             double fmin, double fmax)
{
  bool debug = false;
//   bool debug = true;
  if( debug ) cout << " Calorimeter::DrawFitCellsHistoXYSubset " << GetName() << " debug subset " << sub_set << " total " << GetSubSets().size() <<endl;
  if( sub_set < 0 || sub_set >= (int)GetSubSets().size() ) return;



  if( !GetSubSets()[sub_set].XYInitialized() )
  {
    cerr << " ERROR!! call DrawFitCellsHistoXYSubset for " << GetSubSets()[sub_set].GetName() << " But XY not initialized " << endl;
    return;
  }
  bool ok=false;
  string opts;
//   cout << " Calorimeter::DrawCellsHistoXY " << info_type << endl;
  if(xstep<=0) xstep=1;
  if(xhow_much<=0) xhow_much=1;
  int xmax = xmin + xstep*(xhow_much-1);

  if(ystep<=0) ystep=1;
  if(yhow_much<=0) yhow_much=1;
  int ymax = ymin + ystep*(yhow_much-1);


//   if( !options.fill_histos )
//     return;
  if( debug ) cout << " Calorimeter::DrawFitCellsHistoXYSubset " << GetName() << " check CalorimeterHist " << calo_hist << endl;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  if( debug ) cout << " Calorimeter::DrawFitCellsHistoXYSubset " << GetName() << " CalorimeterHist " << calo_hist << " OK " << endl;
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if( debug )
  {
    if      (info_type == CALIB  ) cout << " CALIB Nhisto " << h.h1_CalibCell.size() << endl;
    else if (info_type == LED  ) cout << " LED Nhisto " << h.h1_LEDCell.size() << endl;
    else if (info_type == PED  ) cout << " PED Nhisto " << h.h1_PEDCell.size() << endl;
    else if (info_type == TIME  ) cout << " TIME Nhisto " << h.h1_TimeCell.size() << endl;
    else if (info_type == NOISE  ) cout << " NOISE Nhisto " << h.h1_NOISECell.size() << endl;
    else  cout << " Info Type unknown " << endl;
  }
  if( debug ) cout << " Calorimeter::DrawFitCellsHistoXYSubset " << GetName() << " OK NCells() " << NCells() << endl;

  if (info_type == CALIB && h.h1_CalibCell.size() != NCells() ) return;
  if (info_type == LED && h.h1_LEDCell.size() != NCells() ) return;
  if (info_type == PED && h.h1_PEDCell.size() != NCells() ) return;
  if (info_type == TIME && h.h1_TimeCell.size() != NCells() ) return;
  if (info_type == NOISE && h.h1_NOISECell.size() != NCells() ) return;

  if( debug ) cout << " Calorimeter::DrawFitCellsHistoXYSubset " << GetName() << " size check OK " << endl;

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);

  bool valid_info_type = true;
  if (info_type == CALIB) sprintf(dir_name,"Calorimeter_%s/CellsCalibration_%s",GetName().c_str(),GetName().c_str());
  else if(info_type == LED) sprintf(dir_name,"Calorimeter_%s/%s_LED_Cells",GetName().c_str(),GetName().c_str());
  else if(info_type == PED) sprintf(dir_name,"Calorimeter_%s/%s_PED_Cells",GetName().c_str(),GetName().c_str());
  else if(info_type == TIME) sprintf(dir_name,"Calorimeter_%s/%s_CellsTime",GetName().c_str(),GetName().c_str());
  else if(info_type == NOISE) sprintf(dir_name,"Calorimeter_%s/%s_CellsNoise",GetName().c_str(),GetName().c_str());
  else if(info_type == RAW) sprintf(dir_name,"Calorimeter_%s/%s_test/%s_Cells",
                                                            GetName().c_str(),GetName().c_str(),GetName().c_str());
  else {
    sprintf(dir_name,"Calorimeter_%s/%s_test/%s_Cells",GetName().c_str(),GetName().c_str(),GetName().c_str());
    valid_info_type = false;
  }

  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to " << dir_name << endl;
    dir_save->cd();
    return;
  }
  if( debug ) cout << " Calorimeter::DrawFitCellsHistoXYSubset " << GetName() << " directory " << dir_name << " OK " <<  endl;

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c)
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  bool canvas_new=false;
  if (h.c == NULL) //   create new canvas
  {
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    canvas_new=true;
  }
  if( debug ) cout << " Calorimeter::DrawFitCellsHistoXYSubset " << GetName() << " opt " << opt << endl;

  if(0!=strcmp(opt,"FIT_ON_OFF" ) )   // In case Change Fit validity dont touch canvas
  {
    h.c->Clear();
    h.c->Divide(xhow_much,yhow_much);
  }
  int ic = 1;
  int fit_status=0;
  bool fit_valid=true;
  for( int iy=ymax; iy >= ymin;  iy -=ystep )  // cycle over cells
    for( int ix=xmin; ix <= xmax;  ix +=xstep )
    {
      int icell = GetSubSets()[sub_set].GetCellOfColumnRow(ix,iy);

      if(0!=strcmp(opt,"FIT_ON_OFF" ) )
      {
        h.c->cd(ic);
        ic++;
      }
      if(icell < 0 || icell >= (int)NCells() )
      {
        cerr << " calorimeter " << GetName() << " DrawFitCellsHistoXYSubset " <<  GetSubSets()[sub_set].GetName() <<
	                                                                 " INTERNAL ERROR !!! Wrong cell " << icell << endl;
        continue;
      }

      char cell_name[132];
//       if      (info_type == CALIB) sprintf(cell_name,"CalibX_%d_Y_%d",ix,iy);
//       else if (info_type == LED)   sprintf(cell_name,"LEDX_%d_Y_%d",ix,iy);
//       else if (info_type == PED)   sprintf(cell_name,"PEDX_%d_Y_%d",ix,iy);
//       else if (info_type == TIME)  sprintf(cell_name,"TimeX_%d_Y_%d",ix,iy);
//       else if (info_type == NOISE)  sprintf(cell_name,"NOISEX_%d_Y_%d",ix,iy);
//       else    sprintf(cell_name,"CellX_%d_Y_%d",ix,iy);
      if      (info_type == CALIB) sprintf(cell_name,"Calib%s",GetCellName(icell).c_str());
      else if (info_type == LED)   sprintf(cell_name,"LED%s",GetCellName(icell).c_str());
      else if (info_type == PED)   sprintf(cell_name,"PED%s",GetCellName(icell).c_str());
      else if (info_type == TIME)  sprintf(cell_name,"Time%s",GetCellName(icell).c_str());
      else if (info_type == NOISE)  sprintf(cell_name,"NOISE%s",GetCellName(icell).c_str());
      else if (info_type == RAW)  sprintf(cell_name,"Cell%s",GetCellName(icell).c_str());
      else    sprintf(cell_name,"Cell%s",GetCellName(icell).c_str());

      TH1D *h1 = (TH1D*) gDirectory->Get(cell_name);


      if(h1 == NULL)
      {
        cerr <<  " Error !! DrawFitCellsHistoXYSubset Histo name mismatch?? " << cell_name << endl;
        continue;
      }
      if(icell < 0 || icell >= (int)NCells() )
      {
        cerr << " FATAL!!  DrawFitCellsHistoXY INTERNAL ERROR !!! Wrong cell " << icell << endl;
        exit(1);
      }
      double scale = 0;
      double maximum = h1->GetXaxis()->GetXmax();
      double minimum = h1->GetXaxis()->GetXmin();
//      cout << " Maximum" << maximum  << " Minimum " << minimum  << " Nbins " << h1->GetNbinsX() << endl;
      if( h1->GetNbinsX() > 0) scale = (maximum-minimum)/h1->GetNbinsX();

      if( 0==strcmp(opt,"DRAW") )          // Draw
      {
        if((fmax-fmin)> 0.000001 && info_type != TIME)
        {
          h1->GetXaxis()->SetRange( h1->GetXaxis()->FindBin(fmin) , h1->GetXaxis()->FindBin(fmax) );
        }
        h1->Draw();
      }
      else if(  0==strcmp(opt,"FIT" ) )      // Fit
      {
        double vmin = minimum;
        double vmax = maximum;
//        cout << " vmin " << vmin << " vmax " << vmax  << endl;
//         if(fmax-fmin> 0.001)
//         {
//           vmin=fmin;
//           vmax=fmax;
//         }

        TF1 *func = new TF1("fit","gaus",vmin,vmax);
//         TF1 *func = new TF1("fit","gaus");
        h1->Fit("fit","R0");
//         cout << " Fit LED in cell " << icell << " Mean " <<func->GetParameter(1) <<
//                " Sigma " << func->GetParameter(2) << " Stat " << func->GetParameter(0) << endl;
        double amplitude = func->GetParameter(0);
        double mean      = func->GetParameter(1);
        double sigma     = func->GetParameter(2);
        double integral  = amplitude*M_PI*sigma/scale;
        if( valid_info_type )
        {
          fit_info[info_type][icell].fit_parameters.clear();
          fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(integral,0));
          fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(mean,0));
          fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(sigma,0));
        }
//         h1->GetXaxis()->SetRange(vmin,vmax);

// //          h1->Integral(bin1,bin2);
// //          c1->SetLogy();
// //          c1->SetLogy(0);
        const Int_t kNotDraw = 1<<9;
//         h1->GetXaxis()->SetRange(binmin,binmax);
        if(h1->GetFunction("fit") != NULL) h1->GetFunction("fit")->ResetBit(kNotDraw);
        h1->Draw();
      }
      else if(0==strcmp(opt,"FIT_ON_OFF" ) && valid_info_type )   // Change Fit validity
      {
        cout << " try FIT_ON_OFF " << fit_status << " info_type " << info_type << " cell " << icell << endl;

        if( icell < 0 || icell >= (int)NCells() ) continue;
        if( fit_status == 0 ) // Check enable or disable fit results
        {
          if( fit_info[info_type][icell].fit_ok==true )
          {
            fit_status=-1;
            fit_valid=false;
          }
          else
          {
            fit_status= 1;
            fit_valid=true;
          }
        }

        if( fit_valid ) // Enable Fit
        {
          fit_info[info_type][icell].fit_ok=true;
          if( fit_info[info_type][icell].fit_parameters.empty() )
          {
             cerr << " ERROR !!! :: FIT CAN'T BE VALID WITHOUT FIT PARAMETERS:" <<
                         GetName() << " X=" << ix << " Y=" << iy << endl;
             fit_info[info_type][icell].fit_ok=false;
          }
        }
        else
        {
          fit_info[info_type][icell].fit_ok=false;
        }
        continue;
      }
      else
      {
         cerr << " ERROR !!! :: " << GetName() <<
                 " DrawFitCellsHistoXY  UNKNOUN option " << opt << endl;
      }
    }                                         // end of cycle over cells

  if(!canvas_new)
  {
    if(0!=strcmp(opt,"FIT_ON_OFF" ) )   // In case Change Fit validity dont touch canvas
    {
      h.c->Update ();
      h.c->Show ();
    }
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void    Calorimeter::DrawMonitorHisto (CellInfoType info_type,const char *opt,
                                               int from_cell,int how_much,int step,
                                               double fmin, double fmax)
{
  bool ok=false;

  bool use_cell_name_to_finde_histo = false;
  string opts;
//   cout << " Calorimeter::DrawCellsHistoXY " << info_type << endl;
  if(step<=0) step=1;
  if(how_much<=0) how_much=1;
  int to_cell = from_cell + step*(how_much-1);

  int nhist = how_much;
  int nx = (int)std::sqrt((float)nhist);
  int ny=nx;
  while( nhist-nx*ny > 0 ) nx += 1;

  if( !options.fill_histos )
    return;

  if( calo_hist==NULL )
  {
    cout << " NO Test Histos  " << endl;
    return;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  char dir_name[111];
  TDirectory *dir_save = gDirectory;
  ok=(gDirectory->cd("/"));
  assert(ok);

//   if (info_type == CALIB) sprintf(dir_name,"Calorimeter_%s/%s_test/%s_Cells",GetName().c_str(),GetName().c_str(),GetName().c_str());
//   else if(info_type == LED) sprintf(dir_name,"Calorimeter_%s/%s_LED_Cells",GetName().c_str(),GetName().c_str());
//   else if(info_type == PED) sprintf(dir_name,"Calorimeter_%s/%s_PED_Cells",GetName().c_str(),GetName().c_str());
//   else if(info_type == TIME) sprintf(dir_name,"Calorimeter_%s/%s_CellsTime",GetName().c_str(),GetName().c_str());
//   else return;

  if( !gDirectory->cd(dir_name) )
  {
    cerr << "Can not change directory to "<< dir_name << endl;
    dir_save->cd();
    return;
  }

//  Checking if canvas in ROOT list
  TSeqCollection *canvases;
  canvases = gROOT->GetListOfCanvases ();
  TCanvas *curcan = (TCanvas *) canvases->First ();
  while (curcan != NULL && curcan != h.c)
    curcan = (TCanvas *) canvases->After((TCanvas *) curcan);

  if (curcan == NULL)
    if (h.c != NULL) h.c = NULL;
// end of check

  bool canvas_new=false;
  if (h.c == NULL) //   create new canvas
  {
    char canvas_name[132];
    sprintf(canvas_name,"%s_c",GetName().c_str());
    h.c = new TCanvas(canvas_name,canvas_name);
    canvas_new=true;
  }

  if(0!=strcmp(opt,"FIT_ON_OFF" ) )   // In case Change Fit validity dont touch canvas
  {
    h.c->Clear();
    h.c->Divide(nx,ny);
  }
  int ic = 1;
    for( int icell=from_cell; icell < to_cell+1; icell +=step )// cycle over cells
    {
      if(0!=strcmp(opt,"FIT_ON_OFF" ) )
      {
        h.c->cd(ic);
        ic++;
      }
      if(icell < 0 || icell >= (int)NCells() )
      {
        cerr << " FATAL!!  DrawFitCellsHistoXY INTERNAL ERROR !!! Wrong cell " <<
                                                                         icell << endl;
        exit(1);
      }
      TH1D *h1 = NULL;
      if( use_cell_name_to_finde_histo )
      {
        char cell_name[132];
        if      (info_type == CALIB) sprintf(cell_name,"Cell_%d",icell);
        else if (info_type == LED)   sprintf(cell_name,"LED_%d",icell);
        else if (info_type == PED)   sprintf(cell_name,"PED_%d",icell);
        else if (info_type == TIME)  sprintf(cell_name,"Time_%d",icell);
        else    return;
        h1 = (TH1D*) gDirectory->Get(cell_name);
      }
      else
      {
        if      (info_type == CALIB) h1=h.h1_ErCell[icell];
        else if (info_type == LED)   h1=h.h1_LEDCell[icell];
        else if (info_type == PED)   h1=h.h1_PEDCell[icell];
        else if (info_type == TIME)  h1=h.h1_TimeCell[icell];
        else    return;
      }


      if(h1 == NULL) continue;

      if( 0==strcmp(opt,"DRAW") )          // Draw
      {
        if((fmax-fmin)> 0.000001)
        {
          h1->GetXaxis()->SetRange((int)fmin ,(int)fmax);
        }
        h1->Draw();
      }
      else
      {
         cerr << " ERROR !!! :: " << GetName() <<
              " DrawFitCellsHistoXY  UNKNOUN option " << opt << endl;
      }
    }                                         // end of cycle over cells

  if(!canvas_new)
  {
      h.c->Update ();
      h.c->Show ();
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FitAllCellsHisto( CellInfoType info_type, const char *opt )
{
  bool debug = true;
//   bool debug = false;

  if( debug ) cout << " Calorimeter::FitAllCellsHisto " << GetName() << " info_type = " << info_type << endl;

  if( calo_hist == NULL ) return;
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  bool no_histo = false;

  if (info_type == CALIB && h.h1_CalibCell.size() != NCells() ) no_histo = true;
  if (info_type == LED && h.h1_LEDCell.size() != NCells() ) no_histo = true;
  if (info_type == PED && h.h1_PEDCell.size() != NCells() ) no_histo = true;
  if (info_type == TIME && h.h1_TimeCell.size() != NCells() ) no_histo = true;
  if (info_type == NOISE && h.h1_NOISECell.size() != NCells() ) no_histo = true;

  if( debug )
  {
    if( no_histo )
      cout << " Calorimeter::FitAllCellsHisto " << GetName() << " we have no histo for info_type = " << info_type << endl;
    else
      cout << " Calorimeter::FitAllCellsHisto " << GetName() << " we have ! histo for info_type = " << info_type <<" go ahead !! " << endl;
  }

  if( no_histo ) return;

  for( int icell = 0; icell < (int)NCells(); icell++ )
  {
    FitCellHisto( icell, info_type, opt);
  }
  if( debug ) cout << " Calorimeter::FitAllCellsHisto " << GetName() << " OK " << endl;

}

////////////////////////////////////////////////////////////////////////////////

class SFitGauss
{
      public:

        /// Destructor
                           ~SFitGauss            (void) {}

        /// Default constructor
                            SFitGauss            (TH1D * h1,double entries_min,
                                                            double stat_max_bin_min,
                                                            double sigma_rel_min,
                                                            double sigma_rel_max,
                                                            double sigma_abs_min,
                                                            double sigma_abs_max,
                                                            double mean_min,
                                                            double mean_max,
                                                            double stat_usefull_min,
                                                            double stat_usefull_max  ) :

							    h1_(h1),
                                                            entries_min_(entries_min),
                                                            stat_max_bin_min_(stat_max_bin_min),
                                                            sigma_rel_min_(sigma_rel_min),
                                                            sigma_rel_max_(sigma_rel_max),
                                                            sigma_abs_min_(sigma_abs_min),
                                                            sigma_abs_max_(sigma_abs_max),
                                                            mean_min_(mean_min),
                                                            mean_max_(mean_max),
                                                            stat_usefull_min_(stat_usefull_min),
                                                            stat_usefull_max_(stat_usefull_max)
							                                       {
											          if( h1_==NULL ) exit(1);
											          Clear();
                                                                                                  bin_width_ = h1_->GetBinWidth(0);
                                                                                                  entries_ = h1_->GetEntries();
                                                                                                  eff_entries_  = h1_->GetEffectiveEntries();
                                                                                                  max_  = h1_->GetMaximum();
                                                                                                  max_bin_  = h1_->GetMaximumBin();
                                                                                                  x_max_bin_  = h1_->GetBinCenter(max_bin_);
                                                                                                  n_bin_  = h1_->GetNbinsX();

												  chi2_=-1.;
                                                                                                  amplitude_=-1.;
                                                                                                  mean_=-1.;
                                                                                                  sigma_=-1.;
                                                                                                  integral_=-1.;
											       }
        /// Clear
        void            Clear                   (void) {}

        /// Fit
        bool            Fit                     (void) {
//							  return FitNormal();
						         return FitNeNormal();
	                                               }
        /// Fit
        bool            FitNeNormal                     (void) {
	                                                  bool fitok = false;
							 // int flag_failure = 0;

							  double frcut =0.10;
							  double xmax = h1_->GetBinCenter(n_bin_)+h1_->GetBinWidth(n_bin_)/2.;
							  vector <int > intervals = GetLogIntervals(xmax, 0.,frcut );

                                                          if( entries_ > entries_min_ )
                                                          {
                                                            fitok = true;

	                                                  }

							  if( fitok )
							  {
							    cout << " xmax value " << xmax << endl;
							    cout << " Intervals " << frcut  << " N intervals " << intervals.size() <<" : " << endl;
							    for( int i =0; i< (int)intervals.size(); i++)
							    {
							      if( i!= 0 ) cout << " --";
							      cout << " " << intervals[i];
							    }
							    cout << endl;

							    cout << " Statistic " << endl;
							    for( int i = 0; i< (int)intervals.size()-1; i++)
							    {
							      int bin_max = intervals[i];
							      int bin_min = intervals[i+1];
							      cout << "  " << GetStatInBins( bin_min, bin_max );

							    }
							    cout << endl;
							  }
                                                          return fitok;
	                                               }
        /// Fit
        bool            FitNormal               (void) {
	                                                  bool fitok = false;
							  int flag_failure = 0;

                                                          if( entries_ > entries_min_ )
                                                          {
	                                                    int max_bin = GetMaxBin();
							    double x_max = h1_->GetBinCenter(max_bin);							  							    double max = h1_->GetBinContent(max_bin);

// 	                                                    double rms = GetFWHMaxBin( max_bin );
	                                                    double rms = x_max*0.1;
	                                                    pair < double,double> mrms = GetMRMSinRange(x_max - 4.*rms, x_max + 4.*rms );

							    if( mrms.second > 0. ) rms = mrms.second;

	                                                    double vmin =  x_max - 5.*rms;
	                                                    double vmax =  x_max + 5.*rms;
                                                            TF1 *ffff = new TF1("ffff","gaus",vmin,vmax);

							    ffff->SetParameters( max, x_max , rms);

							    h1_->Fit("ffff","RQ0N");
                                                            chi2_ = ffff->GetChisquare();
                                                            amplitude_ = ffff->GetParameter(0);
                                                            mean_          = ffff->GetParameter(1);
                                                            sigma_         = ffff->GetParameter(2);
                                                            integral_  = amplitude_*M_PI*sigma_/bin_width_;
                                                            //delete ffff;
							    fitok = true;
                                                            if( mean_ < mean_min_ )
							    {
							      flag_failure += 1;
							      fitok=false;
							    }
                                                            if( mean_ > mean_max_ )
							    {
							      flag_failure += 2;
							      fitok=false;
							    }
                                                            double stat_usefull = integral_/entries_;
                                                            if( stat_usefull < stat_usefull_min_ )
							    {
							      flag_failure += 4;
							      fitok=false;
							    }
                                                            if( stat_usefull > stat_usefull_max_ )
							    {
							      flag_failure += 8;
							      fitok=false;
							    }

							    if( fitok )
							    {
							      h1_->Fit("ffff","RQ0");
                                                              chi2_ = ffff->GetChisquare();
                                                              amplitude_ = ffff->GetParameter(0);
                                                              mean_          = ffff->GetParameter(1);
                                                              sigma_         = ffff->GetParameter(2);
                                                              integral_  = amplitude_*M_PI*sigma_/bin_width_;
                                                              const Int_t kNotDraw = 1<<9;
							      if(h1_->GetFunction("ffff") != NULL) h1_->GetFunction("ffff")->ResetBit(kNotDraw);
	                                                    }
							    else
							    {
                                                               if( max > 0 ) cout << " Flag failure " << flag_failure << endl;
							    }
	                                                 }
                                                         return fitok;
	                                               }

         std::vector <  int >  GetLogIntervals( double from_x, double to_x, double frac_x )
	                                               {
						         std::vector < int > intervals;
							 if( from_x > 0 && to_x >= 0 && frac_x && frac_x < 1. )
							 {
							   double xmax = from_x;
							   while(1)
							   {
							     double xmin = xmax - xmax*frac_x;
							     int bin_max =(int)(h1_->FindBin(xmax));
							     int bin_min =(int)(h1_->FindBin(xmin));
							     if( bin_max <= bin_min )
							     {
							       intervals.push_back(bin_min);
							       break;
							     }
							     intervals.push_back(bin_max);
							     xmax = h1_->GetBinCenter(bin_min)-1.1*h1_->GetBinWidth(bin_min)/2.;
							   }
							 }
							 return intervals;
						       }

         /// \retun bin at maximum
        double            GetStatInBins                ( int bin_min, int bin_max ) {
	                                                  double stat = 0.;
                                                          for( int ibin=bin_min; ibin < bin_max; ibin++ )
	                                                  {
							    stat  += (double)h1_->GetBinContent( (Int_t)ibin );
	                                                  }
							  return stat;
	                                               }
         /// \retun bin at maximum
        int            GetMaxBin                (void) {
                                                          int bin_min = (int)(h1_->FindBin(mean_min_))+1;
	                                                  int max_bin = bin_min;
                                                          double max = (double) h1_->GetBinContent(  (Int_t)bin_min );
                                                          for( int ibin=bin_min; ibin < n_bin_; ibin++ )
	                                                  {
                                                            if( (double)h1_->GetBinContent( (Int_t)ibin ) > max )
	                                                    {
                                                               max = (double)h1_->GetBinContent( ibin );
                                                               max_bin = ibin;
	                                                    }
	                                                  }
							  return max_bin;
	                                               }

         double         GetFWHMaxBin   ( int max_bin ) {
       	                                                  int bin_right = max_bin;
       	                                                  int bin_left = max_bin;
                                                          double max = h1_->GetBinContent(max_bin);
	                                                  while( bin_left > 0 )
	                                                  {
                                                            if( h1_->GetBinContent(bin_left) < max/2. ) break;
	                                                    bin_left--;
	                                                  }

	                                                  while( bin_right < h1_->GetNbinsX() )
	                                                  {
                                                            if( h1_->GetBinContent(bin_right) < max/2. ) break;
	                                                    bin_right++;
	                                                  }

	                                                  return double(bin_right-bin_left+1)*bin_width_;
	                                               }

         std::pair<double,double>    GetMRMSinRange( double xmin, double xmax ) {
                                                          int bin_left =(int)(h1_->FindBin(xmin));
       	                                                  int bin_right =(int)(h1_->FindBin(xmax));
							  if( bin_left > bin_right ) return std::pair<double,double> (0,0);
							  if( bin_left ==  bin_right ) bin_right= bin_left+1;
							  double sum = 0.;
							  double ssum = 0.;
							  for (int ib=bin_left; ib<= bin_right; ib++ )
							  {
							    double x = h1_->GetBinCenter(ib);
                                                            double s = x*h1_->GetBinContent(ib);
							    sum += s;
							    ssum += s*s;
							  }
							  double entr = bin_right-bin_left+1;
							  double ms = sum/entr;
							  double rms = ssum/entr - ms*ms;
							  if( rms > 0 ) rms = sqrt(rms);
							  return std::pair<double,double> (ms,rms);

	                                               }

        void                           Print ( void ) const { std::cout <<" Mean " << mean_ <<" Sigma " << sigma_ <<
	                                                                                       " Integral " << integral_ << std::endl; }

     private:

  TH1D * h1_;
  double entries_min_;
  double stat_max_bin_min_;
  double sigma_rel_min_;
  double sigma_rel_max_;
  double sigma_abs_min_;
  double sigma_abs_max_;
  double mean_min_;
  double mean_max_;
  double stat_usefull_min_;
  double stat_usefull_max_;

  double bin_width_;
  double entries_;
  double eff_entries_;
  double max_;
  int max_bin_;
  double x_max_bin_;
  int n_bin_;

     public:
   Double_t chi2_;
   double amplitude_;
   double mean_;
   double sigma_;
   double integral_;

};

void Calorimeter::FitCellHisto( int icell, CellInfoType info_type, const char *opt )
{
//   bool debug = true;
  bool debug = false;
//   if( GetName() == "HC02P1__" ) debug = true;
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if( debug )
  {
    cout << " FitCellHisto " << GetName() << " cell " << GetCellName(icell) <<" info " <<  info_type << endl;
  }

  double entries_min=10;
  double stat_max_bin_min=5;
  double sigma_rel_min=0.0001;
  double sigma_rel_max=0.2;
  double sigma_abs_min=0.1;
  double sigma_abs_max=100;
  double mean_min=-1.0e+07;
  double mean_max=1.0e+07;
  double stat_usefull_min=0.05;
  double stat_usefull_max=5.;

//   if( icell < 0 || icell >= NCells() ) return;

//   if (info_type == CALIB && h.h1_CalibCell.size() != NCells() ) return;
//   if (info_type == LED && h.h1_LEDCell.size() != NCells() ) return;
//   if (info_type == PED && h.h1_PEDCell.size() != NCells() ) return;
//   if (info_type == TIME && h.h1_TimeCell.size() != NCells() ) return;
//   if (info_type == NOISE && h.h1_NOISECell.size() != NCells() ) return;

  if (info_type == CALIB )
  {
  }
  else if (info_type == LED  )
  {
    entries_min=10;
    stat_max_bin_min=5;
    sigma_rel_min=0.0001;
    sigma_rel_max=0.3;
    sigma_abs_min=0.1;
    sigma_abs_max=200.;
    mean_min=0.;
    mean_max=5000.;
    stat_usefull_min=0.1;
    stat_usefull_max=5.;
  }
  else if (info_type == PED  )
  {
  }
  else if (info_type == TIME  )
  {
  }
  else if (info_type == NOISE  )
  {
  }

  TH1D *h1 = NULL;

  bool valid_info_type = false;
  if	  (info_type == CALIB) h1=h.h1_CalibCell[icell];
  else if (info_type == LED)   h1=h.h1_LEDCell[icell];
  else if (info_type == PED)   h1=h.h1_PEDCell[icell];
  else if (info_type == TIME)  h1=h.h1_TimeCell[icell];
  else if (info_type == NOISE)  h1=h.h1_NOISECell[icell];
  else    return;

  if( h1==NULL ) return;

  valid_info_type = true;
  bool fitok=true;

 SFitGauss sfg( h1,entries_min,
                   stat_max_bin_min,
                   sigma_rel_min,
                   sigma_rel_max,
                   sigma_abs_min,
                   sigma_abs_max,
                   mean_min,
                   mean_max,
                   stat_usefull_min,
                   stat_usefull_max  );

 fitok = sfg.Fit();

 if( valid_info_type )
 {
   fit_info[info_type][icell].fit_parameters.clear();
   fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(sfg.integral_,0));
   fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(sfg.mean_,0));
   fit_info[info_type][icell].fit_parameters.push_back(pair<double,double>(sfg.sigma_,0));
 }

 if( fitok )
 {
   if( debug ) cout << " FitCellHisto Fit OK! cell " << GetCellName(icell) <<
                           " Stat " <<fit_info[info_type][icell].fit_parameters[0].first <<
                           " Mean " << fit_info[info_type][icell].fit_parameters[1].first <<
                           " Sigma " <<fit_info[info_type][icell].fit_parameters[2].first << endl;
 }

 fit_info[info_type][icell].fit_ok=fitok;

}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::ResetHisto      (void)
{
  if(calo_hist !=NULL) calo_hist->Reset();
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::PrintPS      (void)
{
  char ps_file_name[111];
  sprintf(ps_file_name,"%s_canvas.ps",GetName().c_str());
  if(calo_hist !=NULL) calo_hist->c->Print(ps_file_name);
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::TrigGroupTest(void)
{
//  bool debug = true;
  bool debug = false;
  if( trigger_groups_index_.size() <= 0) return;
  if( debug ) cout << " Calorimeter::TrigGroupTest for " << GetName() << " is called " << endl;
  bool ok=false;
  if( !options.fill_histos )
    return;
  if( debug ) cout << " Calorimeter::TrigGroupTest here we have trigger groups " << endl;

  if( gDirectory==NULL )
  {
    cerr << "Calorimeter::TrigGroupTest():  ROOT file was not opend.\n"
         << "                                    I do not want to work.\n";
    return;
  }

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.

// create test histos directory sructure
  if( trig_group_hist_==NULL ) trig_group_hist_ = new TrigGroupHist;
  TrigGroupHist &h=*trig_group_hist_;  // create a short name for calo_hist

  if( h.root_dir==NULL )
  {
    if( debug ) cout << " Calorimeter::TrigGroupTest start booking " << endl;
//     TDirectory *dir_check;
    char dir_name[111];
    ok = gDirectory->cd("/");
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s_TriggerGroups",GetName().c_str());
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
  }
  if( trigger_groups_index_.size() <= 0) return;
  ok = h.root_dir->cd();
  assert(ok);

  if(!h.test_histos_booked)
  {
    double emin = 0, emax = options.hist_energy_max;
    int   n_xbin2 =  100, n_ybin2 =  100;
    //int n_ebin2 =  100;
    double xmin = (*this).GetXmin() , xmax = (*this).GetXmax();
    double ymin = (*this).GetYmin() , ymax = (*this).GetYmax();
    cout << " Calorimeter::TrigGroupTest Ngroups=" << trigger_groups_index_.size() << " in " << GetName() << endl;
    for( int igr=0; igr < (int)trigger_groups_index_.size(); igr++ )
    {
      h.root_dir->cd();
      int sub_set_index = trigger_groups_index_[igr];
      SubSet &subset = calorimeter_sub_sets[sub_set_index];
      char hist_name[132];
//       char dir_name[111];
//       sprintf(dir_name,"%s",subset.GetName().c_str());
//
//       TDirectory* cur_dir = myROOT_utils::TDirectory_checked(dir_name);
//      if( cur_dir == NULL ) continue;
      char name[132];
      sprintf(name,"summrec_in_group%s",subset.GetName().c_str());
      sprintf(hist_name," Summ E of the group %s ",subset.GetName().c_str());
      TH1D * h1 = myROOT_utils::TH1D_checked(name, hist_name,1000, emin, emax );
      h.h1_energy_.push_back(h1);

      sprintf(name,"summax_in_group%s",subset.GetName().c_str());
      sprintf(hist_name," Max Summ E of the group %s ",subset.GetName().c_str());
      h1 = myROOT_utils::TH1D_checked(name, hist_name,1000, emin, emax );
      h.h1_max_energy_.push_back(h1);

      sprintf(name,"xyE_in_group%s",subset.GetName().c_str());
      sprintf(hist_name," Y vs X in group %s weighted with E  ",subset.GetName().c_str());
      TH2D * h2 = myROOT_utils::TH2D_checked(name, hist_name,n_xbin2,xmin,xmax,n_ybin2,ymin,ymax );
      h.h2_XYe_.push_back(h2);

      sprintf(name,"EE_in_group%s",subset.GetName().c_str());
      sprintf(hist_name," Egr vs Ecal in group %s  ",subset.GetName().c_str());
      h2 = myROOT_utils::TH2D_checked(name, hist_name, 200, emin, emax , 200, emin, emax );
      h.h2_Ee_.push_back(h2);

      sprintf(name,"E_in_group%s",subset.GetName().c_str());
      sprintf(hist_name," E measured in group %s  ",subset.GetName().c_str());
      h1 = myROOT_utils::TH1D_checked(name, hist_name, 1000, emin, emax );
      h.h1_E_.push_back(h1);

      sprintf(name,"FindeMe_in_group%s",subset.GetName().c_str());
      sprintf(hist_name," Finde Me in group %s  ",subset.GetName().c_str());
      TProfile *p1 = myROOT_utils::TProfile_checked(name, hist_name, trigger_groups_index_.size(), 0.-0.5, (double)trigger_groups_index_.size()-0.5 );
      h.p1_FindMe_.push_back(p1);
    }
    h.test_histos_booked=true;
  }

  double esummax = 0;
  int max_group = -1;
  for( int igr=0; igr < (int)trigger_groups_index_.size(); igr++ )
  {
    int sub_set_index = trigger_groups_index_[igr];
    SubSet &subset = calorimeter_sub_sets[sub_set_index];
    double esum = 0;
    for ( int i=0; i < (int)subset.GetCells().size(); i++)
    {
      size_t icell = subset.GetCells()[i];
      double x = cells[icell].GetX();
      double y = cells[icell].GetY();

      // search for this cell in the list of signals
      for (vector<CellDataRaw>::const_iterator it=signals_.end(); it!=signals_.begin(); it++)
        if ( it->GetCellIdx()==icell )
        {
          double e = it->GetEnergy();
          if( e > options.cell_energy_threshold )
          {
            esum +=e;
            h.h2_XYe_[igr]->Fill(x,y,e);
          }
        }
    }
    h.h1_energy_[igr]->Fill(esum);
    if( esum > esummax )
    {
      esummax = esum;
      max_group = igr;
    }
    h.h2_Ee_[igr]->Fill( esum, subset.GetEnergy() );
    h.h1_E_[igr]->Fill( subset.GetEnergy());
    if( esum > 0. )
    {
      for( int igr1=0; igr1 < (int)trigger_groups_index_.size(); igr1++ )
      {
        int sub_set_index1 = trigger_groups_index_[igr1];
        SubSet &subset1 = calorimeter_sub_sets[sub_set_index1];
        h.p1_FindMe_[igr]->Fill(igr1, subset1.GetEnergy()/esum );
      }
    }
  }
  if( max_group >= 0 ) h.h1_max_energy_[max_group]->Fill(esummax);

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

TrigGroupHist::TrigGroupHist(void)
{
  test_histos_booked=false;
  root_dir = NULL;
  h1_energy_.clear();
  h1_max_energy_.clear();
  h2_XYe_.clear();
  c = NULL;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillParticlesCaloAssociationHist ( void )
{
  bool debug = false;
//   bool debug = true;
  bool debug_gamma = false;
  debug = options.debug_reconstruction_combined;

//   bool debug = true;
//   bool debug_gamma = true;
  if( debug ) cout << " Calorimeter::FillParticlesCaloAssociationHist level_histos = " <<
                      options.level_particle_calo_association_histo <<" h= " << make_particles_histo_ << endl;
  if( options.level_particle_calo_association_histo <= 0 ) return;
  int levelBPHistos = options.level_particle_calo_association_histo;
  bool is_ECAL2 = (  GetName() == "EC02P1__" );
//  bool is_ECAL1 = (  GetName() == "EC01P1__" );

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.

  if( make_particles_histo_ == NULL )
  {
    double emin = 0., emax = options.hist_energy_max;

    bool ok = gDirectory->cd("/");
    if(!ok) cerr << " Can not change to top dir !! " << endl;
    char dir_name[111];
    sprintf(dir_name,"MakeParticlesHisto_%s",GetName().c_str());
    if( debug ) cout << " new MakeParticlesHisto " << endl;
    make_particles_histo_ = new MakeParticlesHisto( dir_name );
    MakeParticlesHisto &h = *make_particles_histo_;  // create a short name for histo store
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
    if(!h.root_dir->cd()) assert(false);
//     for( std::vector <Reco::Calorimeter*>::const_iterator it = calorimeters_.begin(); it != calorimeters_.end(); it++)
    if( ! h.booked_ok_ )
    {
        if(!h.root_dir->cd()) assert(false);
//         sprintf(dir_name,"Calorimeter_%s",  GetName().c_str());
//         h.calo_dirs = myROOT_utils::TDirectory_checked(dir_name) );
//         if(!h.calo_dirs.back()->cd()) assert(false);
	if( debug ) cout << " MakeParticlesHisto Book for " <<  GetName() << endl;

        double eMax = emax;
//         double hiMax = 20.;
        double hiMax = 50.;
//         double xMin =  GetXmin()+ GetPositionX();
//         double xMax =  GetXmax()+ GetPositionX();
//         double yMin =  GetYmin()+ GetPositionY();
//         double yMax =  GetYmax()+ GetPositionY();
        double xMin =  GetXmin();
        double xMax =  GetXmax();
        double yMin =  GetYmin();
        double yMax =  GetYmax();

        double xMaxCell =  GetMaxCellSizeX();
        double yMaxCell =  GetMaxCellSizeY();

        char name[128];
        char title[100];

        sprintf(name,"caloE");
        sprintf(title," ECalorimeter clusters in  (%s)", GetName().c_str());
        h.caloE = myROOT_utils::TH1D_checked (name,title,10000,0,eMax);
        sprintf(name,"caloSumE");
        sprintf(title," E Summ of Calorimeter clusters in  (%s)", GetName().c_str());
        h.caloSumE = myROOT_utils::TH1D_checked (name,title,1000,0,eMax);
        sprintf(name,"deltaX");
        sprintf(title," XCalorimeter - XTrack in  (%s)", GetName().c_str());
        h.deltaX = myROOT_utils::TH1D_checked (name,title,1000,-500.,500.);

        sprintf(name,"deltaXelectron_InGoodZone");
        sprintf(title," XCalorimeter - X electron InGoodZone in  (%s)", GetName().c_str());
        h.deltaXelectron_InGoodZone = myROOT_utils::TH1D_checked (name,title,2000,-250.,250.);
        sprintf(name,"deltaYelectron_InGoodZone");
        sprintf(title," YCalorimeter - Y electron InGoodZone in  (%s)", GetName().c_str());
        h.deltaYelectron_InGoodZone = myROOT_utils::TH1D_checked (name,title,2000,-250.,250.);
        sprintf(name,"deltaEsElectron_InGoodZone");
        sprintf(title," ECalorimeter/Eelectron-1 InGoodZone  in  (%s)", GetName().c_str());
        h.deltaEsElectron_InGoodZone = myROOT_utils::TH1D_checked (name,title,1000,-0.5,0.5);

        sprintf(name,"deltaXelectron");
        sprintf(title," XCalorimeter - X electron in  (%s)", GetName().c_str());
        h.deltaXelectron = myROOT_utils::TH1D_checked (name,title,2000,-250.,250.);

        sprintf(name,"deltaXelectron_NearBadCell");
        sprintf(title," XCalorimeter - X electron NearBadCell in  (%s)", GetName().c_str());
        h.deltaXelectron_NearBadCell = myROOT_utils::TH1D_checked (name,title,2000,-250.,250.);

        sprintf(name,"deltaXpion");
        sprintf(title," XCalorimeter - X pion in  (%s)", GetName().c_str());
        h.deltaXpion = myROOT_utils::TH1D_checked (name,title,1000,-500.,500.);
        sprintf(name,"deltaXmuon");
        sprintf(title," XCalorimeter - X muon in  (%s)", GetName().c_str());
        h.deltaXmuon = myROOT_utils::TH1D_checked (name,title,1000,-500.,500.);
        sprintf(name,"deltaXwEas");
        sprintf(title," XCalorimeter - XTrack completely associated in  (%s)", GetName().c_str());
        h.deltaXwEas = myROOT_utils::TH1D_checked (name,title,500,-100.,100.);
        sprintf(name,"deltaXS");
        sprintf(title," (XCalorimeter - XTrack)/SigmaXCalo in  (%s)", GetName().c_str());
        h.deltaXS = myROOT_utils::TH1D_checked (name,title,200,-10.,10.);
        sprintf(name,"deltaYwEas");
        sprintf(title," YCalorimeter - YTrack completely associated in  (%s)", GetName().c_str());
        h.deltaYwEas = myROOT_utils::TH1D_checked (name,title,500,-100.,100.);
        sprintf(name,"deltaY");
        sprintf(title," YCalorimeter - YTrack in  (%s)", GetName().c_str());
        h.deltaY = myROOT_utils::TH1D_checked (name,title,1000,-500.,500.);

        sprintf(name,"deltaYelectron_NearBadCell");
        sprintf(title," YCalorimeter - Y electron NearBadCell in  (%s)", GetName().c_str());
        h.deltaYelectron_NearBadCell = myROOT_utils::TH1D_checked (name,title,2000,-250.,250.);

        sprintf(name,"deltaYelectron");
        sprintf(title," YCalorimeter - Y electron in  (%s)", GetName().c_str());
        h.deltaYelectron = myROOT_utils::TH1D_checked (name,title,2000,-250.,250.);

        sprintf(name,"deltaYpion");
        sprintf(title," YCalorimeter - Y pion in  (%s)", GetName().c_str());
        h.deltaYpion = myROOT_utils::TH1D_checked (name,title,1000,-500.,500.);
        sprintf(name,"deltaYmuon");
        sprintf(title," YCalorimeter - Y muon in  (%s)", GetName().c_str());
        h.deltaYmuon = myROOT_utils::TH1D_checked (name,title,1000,-500.,500.);
        sprintf(name,"deltaYS");
        sprintf(title," (YCalorimeter - YTrack)/SigmaYCalo in  (%s)", GetName().c_str());
        h.deltaYS = myROOT_utils::TH1D_checked (name,title,200,-10.,10.);
        sprintf(name,"deltaEs");
        sprintf(title," ECalorimeter/Etrack - 1  in  (%s)", GetName().c_str());
        h.deltaEs = myROOT_utils::TH1D_checked (name,title,1000,-0.5,0.5);

        sprintf(name,"deltaEsElectron");
        sprintf(title," ECalorimeter-Eelectron  in  (%s)", GetName().c_str());
        h.deltaEselectron = myROOT_utils::TH1D_checked (name,title,1000,-0.5,0.5);

        sprintf(name,"deltaEsElectron_NearBadCell");
        sprintf(title," ECalorimeter/Eelectron-1 NearBadCell  in  (%s)", GetName().c_str());
        h.deltaEsElectron_NearBadCell = myROOT_utils::TH1D_checked (name,title,1000,-0.5,0.5);


        sprintf(name,"deltaEsPion");
        sprintf(title," ECalorimeter-Epion  in  (%s)", GetName().c_str());
        h.deltaEspion = myROOT_utils::TH1D_checked (name,title,1000,-0.5,0.5);
        sprintf(name,"deltaEsMuon");
        sprintf(title," ECalorimeter-Emuon  in  (%s)", GetName().c_str());
        h.deltaEsmuon = myROOT_utils::TH1D_checked (name,title,1000,-0.5,0.5);

        sprintf(name,"neutralE");
        sprintf(title," E neutrals  in  (%s)", GetName().c_str());
        h.neutralE = myROOT_utils::TH1D_checked (name,title,1000,0.,eMax);
        sprintf(name,"gammaE");
        sprintf(title," E gammas  in  (%s)", GetName().c_str());
        h.gammaE = myROOT_utils::TH1D_checked (name,title,1000,0.,eMax);
        sprintf(name,"electronE");
        sprintf(title," E electrons-positrons  in  (%s)", GetName().c_str());
        h.electronE = myROOT_utils::TH1D_checked (name,title,1000,0.,eMax);
        sprintf(name,"pionE");
        sprintf(title," E pions  in  (%s)", GetName().c_str());
        h.pionE = myROOT_utils::TH1D_checked (name,title,1000,0.,eMax);
        sprintf(name,"muonE");
        sprintf(title," E muons  in  (%s)", GetName().c_str());
        h.muonE = myROOT_utils::TH1D_checked (name,title,1000,0.,eMax);
        if(levelBPHistos > 1)
	{
//           int fact1Dim = levelBPHistos - 1;
          int fact2Dim = levelBPHistos - 1;
          sprintf(name,"caloXtrackX");
          sprintf(title," XCalorimeter vs XTrack in  (%s)", GetName().c_str());
          h.caloXtrackX = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,xMin,xMax);
          sprintf(name,"caloYtrackY");
          sprintf(title," YCalorimeter vs YTrack in  (%s)", GetName().c_str());
          h.caloYtrackY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,yMin,yMax,100*fact2Dim,yMin,yMax);

  //  caloXCellX   caloYCellY

          sprintf(name,"caloXtrackXcell");
          sprintf(title," XCalorimeter vs XTrack in the cell of (%s)", GetName().c_str());
          h.caloXtrackXcell = myROOT_utils::TH2D_checked (name,title,50,-2.*xMaxCell,2.*xMaxCell,50,-2.*xMaxCell,2.*xMaxCell);

          sprintf(name,"caloYtrackYcell");
          sprintf(title," YCalorimeter vs YTrack in the cell of  (%s)", GetName().c_str());
          h.caloYtrackYcell = myROOT_utils::TH2D_checked (name,title,50,-2.*yMaxCell,2.*yMaxCell,50,-2.*yMaxCell,2.*yMaxCell);

          sprintf(name,"caloXmuonXcell");
          sprintf(title," XCalorimeter vs Xmuon in the cell of (%s)", GetName().c_str());
          h.caloXmuonXcell = myROOT_utils::TH2D_checked (name,title,50,-2.*xMaxCell,2.*xMaxCell,50,-2.*xMaxCell,2.*xMaxCell);

          sprintf(name,"caloYmuonYcell");
          sprintf(title," YCalorimeter vs Ymuon in the cell of  (%s)", GetName().c_str());
          h.caloYmuonYcell = myROOT_utils::TH2D_checked (name,title,50,-2.*yMaxCell,2.*yMaxCell,50,-2.*yMaxCell,2.*yMaxCell);

          sprintf(name,"caloXpionXcell");
          sprintf(title," XCalorimeter vs Xpion in the cell of (%s)", GetName().c_str());
          h.caloXpionXcell = myROOT_utils::TH2D_checked (name,title,50,-2.*xMaxCell,2.*xMaxCell,50,-2.*xMaxCell,2.*xMaxCell);

          sprintf(name,"caloYpionYcell");
          sprintf(title," YCalorimeter vs Ypion in the cell of  (%s)", GetName().c_str());
          h.caloYpionYcell = myROOT_utils::TH2D_checked (name,title,50,-2.*yMaxCell,2.*yMaxCell,50,-2.*yMaxCell,2.*yMaxCell);

          sprintf(name,"caloXelectronXcell");
          sprintf(title," XCalorimeter vs Xelectron in the cell of (%s)", GetName().c_str());
          h.caloXelectronXcell = myROOT_utils::TH2D_checked (name,title,50,-2.*xMaxCell,2.*xMaxCell,50,-2.*xMaxCell,2.*xMaxCell);

          sprintf(name,"caloYelectronYcell");
          sprintf(title," YCalorimeter vs Yelectron in the cell of  (%s)", GetName().c_str());
          h.caloYelectronYcell = myROOT_utils::TH2D_checked (name,title,50,-2.*yMaxCell,2.*yMaxCell,50,-2.*yMaxCell,2.*yMaxCell);

  //   deltaXCellX   deltaYCellY

          sprintf(name,"deltaXtrackXcell");
          sprintf(title," XTrack vs (XCalorimeter - XTrack) in the cell of (%s)", GetName().c_str());
          h.deltaXtrackXcell = myROOT_utils::TH2D_checked (name,title,100,-xMaxCell,xMaxCell,50,-2.*xMaxCell,2.*xMaxCell);

          sprintf(name,"deltaYtrackYcell");
          sprintf(title," YTrack vs (YCalorimeter - YTrack) in the cell of  (%s)", GetName().c_str());
          h.deltaYtrackYcell = myROOT_utils::TH2D_checked (name,title,100,-yMaxCell,yMaxCell,50,-2.*yMaxCell,2.*yMaxCell);

          sprintf(name,"deltaXmuonXcell");
          sprintf(title," XTrack vs (XCalorimeter - XTrack) in the cell of (%s)", GetName().c_str());
          h.deltaXmuonXcell = myROOT_utils::TH2D_checked (name,title,100,-xMaxCell,xMaxCell,50,-2.*xMaxCell,2.*xMaxCell);

          sprintf(name,"deltaYmuonYcell");
          sprintf(title," YTrack vs (YCalorimeter - YTrack) in the cell of  (%s)", GetName().c_str());
          h.deltaYmuonYcell = myROOT_utils::TH2D_checked (name,title,100,-yMaxCell,yMaxCell,50,-2.*yMaxCell,2.*yMaxCell);

          sprintf(name,"deltaXpionXcell");
          sprintf(title," XTrack vs (XCalorimeter - XTrack) in the cell of (%s)", GetName().c_str());
          h.deltaXpionXcell = myROOT_utils::TH2D_checked (name,title,100,-xMaxCell,xMaxCell,50,-2.*xMaxCell,2.*xMaxCell);

          sprintf(name,"deltaYpionYcell");
          sprintf(title," YTrack vs (YCalorimeter - YTrack) in the cell of  (%s)", GetName().c_str());
          h.deltaYpionYcell = myROOT_utils::TH2D_checked (name,title,100,-yMaxCell,yMaxCell,50,-2.*yMaxCell,2.*yMaxCell);

          sprintf(name,"deltaXelectronXcell");
          sprintf(title," XTrack vs (XCalorimeter - XTrack) in the cell of (%s)", GetName().c_str());
          h.deltaXelectronXcell = myROOT_utils::TH2D_checked (name,title,100,-xMaxCell,xMaxCell,50,-2.*xMaxCell,2.*xMaxCell);

          sprintf(name,"deltaYelectronYcell");
          sprintf(title," YTrack vs (YCalorimeter - YTrack) in the cell of  (%s)", GetName().c_str());
          h.deltaYelectronYcell = myROOT_utils::TH2D_checked (name,title,100,-yMaxCell,yMaxCell,50,-2.*yMaxCell,2.*yMaxCell);

          sprintf(name,"deltaXelectronXcell_selected");
          sprintf(title," XTrack vs (XCalorimeter - XTrack)_selected in the cell of (%s)", GetName().c_str());
          h.deltaXelectronXcell_selected = myROOT_utils::TH2D_checked (name,title,100,-xMaxCell,xMaxCell,50,-2.*xMaxCell,2.*xMaxCell);

          sprintf(name,"deltaYelectronYcell_selected");
          sprintf(title," YTrack vs (YCalorimeter - YTrack)_selected in the cell of  (%s)", GetName().c_str());
          h.deltaYelectronYcell_selected = myROOT_utils::TH2D_checked (name,title,100,-yMaxCell,yMaxCell,50,-2.*yMaxCell,2.*yMaxCell);

     //

          sprintf(name,"trackXY");
          sprintf(title," Y vs X Track in  (%s)", GetName().c_str());
          h.trackXY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"caloXY");
          sprintf(title," Y vs X CaloCluster in  (%s)", GetName().c_str());
          h.caloXY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);

          sprintf(name,"deltaXdeltaY");
          sprintf(title," (YCalorimeter - YTrack) vs (XCalorimeter - XTrack) in  (%s)", GetName().c_str());
          h.deltaXdeltaY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-500.,500.,100*fact2Dim,-500.,500.);

          sprintf(name,"deltaXdeltaYall");
          sprintf(title," (YCalorimeter - YTrack) vs (XCalorimeter - XTrack) for all in  (%s)", GetName().c_str());
          h.deltaXdeltaYall = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2500.,2500.,100*fact2Dim,-2500.,2500.);

          sprintf(name,"deltaXtrackX");
          sprintf(title," XTrack vs (XCalorimeter - XTrack) in  (%s)", GetName().c_str());
          h.deltaXtrackX = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*xMaxCell,2.5*xMaxCell,100*fact2Dim,xMin,xMax);

          sprintf(name,"MeanDeltaXelectronX");
          sprintf(title," X electron vs Mean(XCalorimeter - X electron) in  (%s)", GetName().c_str());
          h.MeanDeltaXelectronX = myROOT_utils::TProfile_checked (name,title,100*fact2Dim,xMin,xMax);

          sprintf(name,"deltaXelectronX");
          sprintf(title," X electron vs (XCalorimeter - X electron) in  (%s)", GetName().c_str());
          h.deltaXelectronX = myROOT_utils::TH2D_checked (name,title,200*fact2Dim,-2.5*xMaxCell,2.5*xMaxCell,100*fact2Dim,xMin,xMax);

          sprintf(name,"MeanDeltaXelectronX_selected");
          sprintf(title," X electron vs Mean(XCalorimeter - X electron)_selected in  (%s)", GetName().c_str());
          h.MeanDeltaXelectronX_selected = myROOT_utils::TProfile_checked (name,title,100*fact2Dim,xMin,xMax);

          sprintf(name,"deltaXelectronX_selected");
          sprintf(title," X electron vs (XCalorimeter - X electron)_selected in  (%s)", GetName().c_str());
          h.deltaXelectronX_selected = myROOT_utils::TH2D_checked (name,title,200*fact2Dim,-2.5*xMaxCell,2.5*xMaxCell,100*fact2Dim,xMin,xMax);

          sprintf(name,"deltaXpionX");
          sprintf(title," X pion vs (XCalorimeter - X pion) in  (%s)", GetName().c_str());
          h.deltaXpionX = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*xMaxCell,2.5*xMaxCell,100*fact2Dim,xMin,xMax);
          sprintf(name,"deltaXmuonX");
          sprintf(title," X muon vs (XCalorimeter - X muon) in  (%s)", GetName().c_str());
          h.deltaXmuonX = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*xMaxCell,xMaxCell,100*fact2Dim,xMin,xMax);

          sprintf(name,"deltaXelectronE");
          sprintf(title," Ecalo electron vs (XCalorimeter - X electron) in  (%s)", GetName().c_str());
          h.deltaXelectronE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*xMaxCell,2.5*xMaxCell,100,0.,eMax);
          sprintf(name,"deltaXpionE");
          sprintf(title," Ecalo pion vs (XCalorimeter - X pion) in  (%s)", GetName().c_str());
          h.deltaXpionE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*xMaxCell,2.5*xMaxCell,100,0.,eMax);
          sprintf(name,"deltaXmuonE");
          sprintf(title," Ecalo muon vs (XCalorimeter - X muon) in  (%s)", GetName().c_str());
          h.deltaXmuonE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*xMaxCell,2.5*xMaxCell,100,0.,eMax);

          sprintf(name,"deltaYelectronE");
          sprintf(title," Ecalo electron vs (YCalorimeter - Y electron) in  (%s)", GetName().c_str());
          h.deltaYelectronE = myROOT_utils::TH2D_checked (name,title,200*fact2Dim,-2.5*yMaxCell,2.5*yMaxCell,100,0.,eMax);

          sprintf(name,"deltaYpionE");
          sprintf(title," Ecalo pion vs (YCalorimeter - Y pion) in  (%s)", GetName().c_str());
          h.deltaYpionE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*yMaxCell,2.5*yMaxCell,100,0.,eMax);

          sprintf(name,"deltaYmuonE");
          sprintf(title," Ecalo muon vs (YCalorimeter - Y muon) in  (%s)", GetName().c_str());
          h.deltaYmuonE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*yMaxCell,2.5*yMaxCell,100,0.,eMax);

          sprintf(name,"deltaYtrackY");
          sprintf(title," YTrack vs (YCalorimeter - YTrack) in  (%s)", GetName().c_str());
          h.deltaYtrackY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*yMaxCell,2.5*yMaxCell,100*fact2Dim,yMin,yMax);

          sprintf(name,"deltaYelectronY");
          sprintf(title," Y electron vs (YCalorimeter - Y electron) in  (%s)", GetName().c_str());
          h.deltaYelectronY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*yMaxCell,2.5*yMaxCell,100*fact2Dim,yMin,yMax);

          sprintf(name,"MeanDeltaYelectronY");
          sprintf(title," Y electron vs Mean(YCalorimeter - Y electron) in  (%s)", GetName().c_str());
          h.MeanDeltaYelectronY = myROOT_utils::TProfile_checked (name,title,100*fact2Dim,yMin,yMax);

          sprintf(name,"deltaYelectronY_selected");
          sprintf(title," Y electron vs (YCalorimeter - Y electron)_selected in  (%s)", GetName().c_str());
          h.deltaYelectronY_selected = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*yMaxCell,2.5*yMaxCell,100*fact2Dim,yMin,yMax);

          sprintf(name,"MeanDeltaYelectronY_selected");
          sprintf(title," Y electron vs Mean(YCalorimeter - Y electron)_selected in  (%s)", GetName().c_str());
          h.MeanDeltaYelectronY_selected = myROOT_utils::TProfile_checked (name,title,100*fact2Dim,yMin,yMax);

          sprintf(name,"deltaYpionY");
          sprintf(title," Y pion vs (YCalorimeter - Y pion) in  (%s)", GetName().c_str());
          h.deltaYpionY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*yMaxCell,2.5*yMaxCell,100*fact2Dim,yMin,yMax);
          sprintf(name,"deltaYmuonY");
          sprintf(title," Y muon vs (YCalorimeter - Y muon) in  (%s)", GetName().c_str());
          h.deltaYmuonY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-2.5*yMaxCell,2.5*yMaxCell,100*fact2Dim,yMin,yMax);

          sprintf(name,"deltaSqXtrackSqX");
          sprintf(title," covTrackX vs (XCalorimeter - XTrack)**2 in  (%s)", GetName().c_str());
          h.deltaSqXtrackSqX = myROOT_utils::TH2D_checked (name,title,100,0.,50.,100,0.,50.);
          sprintf(name,"deltaSqXelectronSqX");
          sprintf(title," covElectronX vs (XCalorimeter - XElectron)**2 in  (%s)", GetName().c_str());
          h.deltaSqXelectronSqX = myROOT_utils::TH2D_checked (name,title,100,0.,50.,100,0.,50.);
          sprintf(name,"deltaSqXpionSqX");
          sprintf(title," covPionX vs (XCalorimeter - Xpion)**2 in  (%s)", GetName().c_str());
          h.deltaSqXpionSqX = myROOT_utils::TH2D_checked (name,title,100,0.,50.,100,0.,50.);
          sprintf(name,"deltaSqXmuonSqX");
          sprintf(title," covMuonX vs (XCalorimeter - Xmuon)**2 in  (%s)", GetName().c_str());
          h.deltaSqXmuonSqX = myROOT_utils::TH2D_checked (name,title,100,0.,50.,100,0.,50.);

          sprintf(name,"deltaSqYtrackSqY");
          sprintf(title," covTrackY vs (YCalorimeter - YTrack)**2 in  (%s)", GetName().c_str());
          h.deltaSqYtrackSqY = myROOT_utils::TH2D_checked (name,title,100,0.,50.,100,0.,50.);
          sprintf(name,"deltaSqYelectronSqY");
          sprintf(title," covElectronY vs (YCalorimeter - YElectron)**2 in  (%s)", GetName().c_str());
          h.deltaSqYelectronSqY = myROOT_utils::TH2D_checked (name,title,100,0.,50.,100,0.,50.);
          sprintf(name,"deltaSqYpionSqY");
          sprintf(title," covPionY vs (YCalorimeter - Ypion)**2 in  (%s)", GetName().c_str());
          h.deltaSqYpionSqY = myROOT_utils::TH2D_checked (name,title,100,0.,50.,100,0.,50.);
          sprintf(name,"deltaSqYmuonSqY");
          sprintf(title," covMuonY vs (YCalorimeter - Ymuon)**2 in  (%s)", GetName().c_str());
          h.deltaSqYmuonSqY = myROOT_utils::TH2D_checked (name,title,100,0.,50.,100,0.,50.);

           sprintf(name,"deltaXtrackXwEas");
          sprintf(title," XTrack vs (XCalorimeter - XTrack) in  (%s) associated with E ", GetName().c_str());
          h.deltaXtrackXwEas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-500.,500.,100*fact2Dim,xMin,xMax);
          sprintf(name,"deltaYtrackYwEas");
          sprintf(title," YTrack vs (YCalorimeter - YTrack) in  (%s) associated with E ", GetName().c_str());
          h.deltaYtrackYwEas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-500.,500.,100*fact2Dim,yMin,yMax);
          sprintf(name,"trackXYas");
          sprintf(title," YTrack vs XTrack associated in  (%s)", GetName().c_str());
          h.trackXYas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);

          sprintf(name,"electronXYas");
          sprintf(title," Yelectron vs Xelectron associated in  (%s)", GetName().c_str());
          h.electronXYas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"pionXYas");
          sprintf(title," Ypion vs Xpion associated in  (%s)", GetName().c_str());
          h.pionXYas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"muonXYas");
          sprintf(title," Ymuon vs Xmuon associated in  (%s)", GetName().c_str());
          h.muonXYas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);


          sprintf(name,"trackXYnotas");
          sprintf(title," YTrack vs XTrack hit but not associated in  (%s)", GetName().c_str());
          h.trackXYnotas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);

          sprintf(name,"electronXYnotas");
          sprintf(title," Yelectron vs Xelectron not associated in  (%s)", GetName().c_str());
          h.electronXYnotas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"pionXYnotas");
          sprintf(title," Ypion vs Xpion not associated in  (%s)", GetName().c_str());
          h.pionXYnotas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"muonXYnotas");
          sprintf(title," Ymuon vs Xmuon not associated in  (%s)", GetName().c_str());
          h.muonXYnotas = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);

          sprintf(name,"trackXYas_E");
          sprintf(title," YTrack vs XTrack( with momentum) associated in  (%s)", GetName().c_str());
          h.trackXYas_E = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"trackXYnotas_E");
          sprintf(title," YTrack vs XTrack( with momentum) not associated in  (%s)", GetName().c_str());
          h.trackXYnotas_E = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"trackXYas_noE");
          sprintf(title," YTrack vs XTrack( without momentum) associated in  (%s)", GetName().c_str());
          h.trackXYas_noE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"trackXYnotas_noE");
          sprintf(title," YTrack vs XTrack( without momentum) not associated in  (%s)", GetName().c_str());
          h.trackXYnotas_noE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"caloEtrackE");
          sprintf(title," ETrack vs ECalorimeter associated in  (%s)", GetName().c_str());
          h.caloEtrackE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,0,eMax,100*fact2Dim,0,eMax);
          sprintf(name,"deltaEtrackE");
          sprintf(title," ECalorimeter-ETrack vs Etrack associated in  (%s)", GetName().c_str());
          h.deltaEtrackE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-50.,50.,100*fact2Dim,0,eMax);
          sprintf(name,"deltaEsEtrack");
          sprintf(title," ETrack vs ECalorimeter - Etrack associated in  (%s)", GetName().c_str());
          h.deltaEsEtrack = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-0.5,0.5,100*fact2Dim,0,eMax);

          sprintf(name,"deltaEsEelectron");
          sprintf(title," Eelectron vs ECalorimeter - Eelectron associated in  (%s)", GetName().c_str());
          h.deltaEsEelectron = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-0.5,0.5,100*fact2Dim,0,eMax);

          sprintf(name,"MeanDeltaEsElectronY");
          sprintf(title," Mean ECalorimeter/Eelectron-1 as function of Y electron in  (%s)", GetName().c_str());
          h.MeanDeltaEsElectronY = myROOT_utils::TProfile_checked (name,title,50*fact2Dim,yMin,yMax);

          sprintf(name,"MeanDeltaEsElectronX");
          sprintf(title," Mean ECalorimeter/Eelectron-1 as function of X electron in  (%s)", GetName().c_str());
          h.MeanDeltaEsElectronX = myROOT_utils::TProfile_checked (name,title,50*fact2Dim,yMin,yMax);
//           if(  XYRegularGrid() )
//           {
//              TProfile** arr = new (TProfile*) [ GetNRows()];
//              h.MeanDeltaEsElectronXRow.push_back(arr);
//              for( int ii=0; ii< GetNRows(); ii++)
//              {
//                sprintf(name,"MeanDeltaEsElectronXRow_%d",ii);
//                sprintf(title," Mean ECalorimeter/Eelectron-1 as function of X electron in row %d  (%s)",ii, GetName().c_str());
//                arr[ii] = myROOT_utils::TProfile_checked (name,title,50*fact2Dim,xMin,xMax);
//              }
//
//              arr = new (TProfile*) [ GetNColumns()];
//              h.MeanDeltaEsElectronYColumn.push_back(arr);
//              for( int ii=0; ii< GetNColumns(); ii++)
//              {
//                sprintf(name,"MeanDeltaEsElectronYColumn_%d",ii);
//                sprintf(title," Mean ECalorimeter/Eelectron-1 as function of Y electron in column %d  (%s)",ii, GetName().c_str());
//                arr[ii] = myROOT_utils::TProfile_checked (name,title,50*fact2Dim,yMin,yMax);
//              }
//
//           }

          sprintf(name,"deltaEsEpion");
          sprintf(title," Epion vs ECalorimeter - Epion associated in  (%s)", GetName().c_str());
          h.deltaEsEpion = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-0.5,0.5,100*fact2Dim,0,eMax);
          sprintf(name,"deltaEsEmuon");
          sprintf(title," Emuon vs ECalorimeter - Emuon associated in  (%s)", GetName().c_str());
          h.deltaEsEmuon = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-0.5,0.5,100*fact2Dim,0,eMax);

          sprintf(name,"caloTtrackT");
          sprintf(title," TimeTrack vs TimeCalorimeter associated in  (%s)", GetName().c_str());
          h.caloTtrackT = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,-50,50,100*fact2Dim,-50,50);
          sprintf(name,"neutralXY");
          sprintf(title," Y vs X of neutral patricles in  (%s)", GetName().c_str());
          h.neutralXY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"gammaXY");
          sprintf(title," Y vs X of gamma patricles in  (%s)", GetName().c_str());
          h.gammaXY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"electronXY");
          sprintf(title," Y vs X of electron-positron patricles in  (%s)", GetName().c_str());
          h.electronXY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"pionXY");
          sprintf(title," Y vs X of pion patricles in  (%s)", GetName().c_str());
          h.pionXY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);
          sprintf(name,"muonXY");
          sprintf(title," Y vs X of muon patricles in  (%s)", GetName().c_str());
          h.muonXY = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,xMin,xMax,100*fact2Dim,yMin,yMax);

          sprintf(name,"caloEelectronE");
          sprintf(title," E electrons-positrons vs ECalorimeter associated in  (%s)", GetName().c_str());
          h.caloEelectronE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,0,eMax,100*fact2Dim,0,eMax);
          sprintf(name,"caloEpionE");
          sprintf(title," E pions vs ECalorimeter associated in  (%s)", GetName().c_str());
          h.caloEpionE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,0,eMax,100*fact2Dim,0,eMax);
          sprintf(name,"caloEmuonE");
          sprintf(title," E muons vs ECalorimeter associated in  (%s)", GetName().c_str());
          h.caloEmuonE = myROOT_utils::TH2D_checked (name,title,100*fact2Dim,0,eMax,100*fact2Dim,0,eMax);

          sprintf(name,"distGamPiE");
          sprintf(title," E pion deposit vs pion-gamma dist in ECalorimeter(%s)", GetName().c_str());
          h.distGamPiE = myROOT_utils::TH2D_checked (name,title,100,0.,2000.,50*fact2Dim,0,100.);

          sprintf(name,"distGamGamE");
          sprintf(title," E gamma deposit vs gamma-gamma dist in ECalorimeter(%s)", GetName().c_str());
          h.distGamGamE = myROOT_utils::TH2D_checked (name,title,100,0.,2000.,100*fact2Dim,0,200.);

          sprintf(name,"distXYGamGam");
          sprintf(title," distx gamma-gamma deposit vs disty gamma-gamma in ECalorimeter(%s)", GetName().c_str());
          h.distXYGamGam = myROOT_utils::TH2D_checked (name,title,100,-1000.,1000.,100,-500.,500.);

          sprintf(name,"distXYElectronGam");
          sprintf(title," distx Electron-gamma deposit vs disty Electron-gamma in ECalorimeter(%s)", GetName().c_str());
          h.distXYElectronGam = myROOT_utils::TH2D_checked (name,title,100,-1000.,1000.,100,-500.,500.);

          sprintf(name,"distXYPiGam");
          sprintf(title," distx pion-gamma deposit vs disty pion-gamma in ECalorimeter(%s)", GetName().c_str());
          h.distXYPiGam = myROOT_utils::TH2D_checked (name,title,100,-1000.,1000.,100,-500.,500.);

          sprintf(name,"muonEmE4low");
          sprintf(title," E out of 4 vs E4 out of Main for muons with low energy(below 6) deposit in   (%s)", GetName().c_str());
          h.muonEmE4low = myROOT_utils::TH2D_checked (name,title,120,-0.1,1.1,120,-0.1,1.1);
          sprintf(name,"muonEmE4high");
          sprintf(title," E out of 4 vs E4 out of Main for muons with low energy(above 6) deposit in   (%s)", GetName().c_str());
          h.muonEmE4high = myROOT_utils::TH2D_checked (name,title,120,-0.1,1.1,120,-0.1,1.1);

          sprintf(name,"pionEmE4low");
          sprintf(title," E out of 4 vs E4 out of Main for pions with low energy(below 6) deposit in   (%s)", GetName().c_str());
          h.pionEmE4low = myROOT_utils::TH2D_checked (name,title,120,-0.1,1.1,120,-0.1,1.1);
          sprintf(name,"pionEmE4high");
          sprintf(title," E out of 4 vs E4 out of Main for pions with low energy(above 6) deposit in   (%s)", GetName().c_str());
          h.pionEmE4high = myROOT_utils::TH2D_checked (name,title,120,-0.1,1.1,120,-0.1,1.1);

          sprintf(name,"electronEmE4low");
          sprintf(title," E out of 4 vs E4 out of Main for electrons with low energy(below 6) deposit in   (%s)", GetName().c_str());
          h.electronEmE4low = myROOT_utils::TH2D_checked (name,title,120,-0.1,1.1,120,-0.1,1.1);
          sprintf(name,"electronEmE4high");
          sprintf(title," E out of 4 vs E4 out of Main for electrons with low energy(above 6) deposit in   (%s)", GetName().c_str());
          h.electronEmE4high = myROOT_utils::TH2D_checked (name,title,120,-0.1,1.1,120,-0.1,1.1);

          sprintf(name,"gammaEmE4low");
          sprintf(title," E out of 4 vs E4 out of Main for gammas with low energy(below 6) deposit in   (%s)", GetName().c_str());
          h.gammaEmE4low = myROOT_utils::TH2D_checked (name,title,120,-0.1,1.1,120,-0.1,1.1);
          sprintf(name,"gammaEmE4high");
          sprintf(title," E out of 4 vs E4 out of Main for gammas with low energy(above 6) deposit in   (%s)", GetName().c_str());
          h.gammaEmE4high = myROOT_utils::TH2D_checked (name,title,120,-0.1,1.1,120,-0.1,1.1);

       }

       sprintf(name,"MGamGam");
       sprintf(title," Mass gamma-gamma in ECalorimeter(%s)", GetName().c_str());
       h.MGamGam = myROOT_utils::TH1D_checked (name,title,1000,0.,2.);

       sprintf(name,"MGamGam_NearBadCell");
       sprintf(title," Mass gamma-gamma Near Bad Cell in ECalorimeter(%s)", GetName().c_str());
       h.MGamGam_NearBadCell = myROOT_utils::TH1D_checked (name,title,1000,0.,2.);

       sprintf(name,"Ngamma");
       sprintf(title," N gamma patricles in  (%s)", GetName().c_str());
       h.Ngamma = myROOT_utils::TH1D_checked (name,title,100,0.,100.);
       sprintf(name,"NKl");
       sprintf(title," N Kl patricles in  (%s)", GetName().c_str());
       h.NKl = myROOT_utils::TH1D_checked (name,title,100,0.,100.);
       sprintf(name,"NPiplus");
       sprintf(title," N Piplus patricles in  (%s)", GetName().c_str());
       h.NPiplus = myROOT_utils::TH1D_checked (name,title,100,0.,100.);
       sprintf(name,"NPiminus");
       sprintf(title," N Piminus patricles in  (%s)", GetName().c_str());
       h.NPiminus = myROOT_utils::TH1D_checked (name,title,100,0.,100.);
       sprintf(name,"NMuplus");
       sprintf(title," N Muplus patricles in  (%s)", GetName().c_str());
       h.NMuplus = myROOT_utils::TH1D_checked (name,title,100,0.,100.);
       sprintf(name,"NMuminus");
       sprintf(title," N Muminus patricles in  (%s)", GetName().c_str());
       h.NMuminus = myROOT_utils::TH1D_checked (name,title,100,0.,100.);
       sprintf(name,"NElectrons");
       sprintf(title," N electrons in  (%s)", GetName().c_str());
       h.NElectrons = myROOT_utils::TH1D_checked (name,title,100,0.,100.);
       sprintf(name,"NPositrons");
       sprintf(title," N positrons in  (%s)", GetName().c_str());
       h.NPositrons = myROOT_utils::TH1D_checked (name,title,100,0.,100.);

        sprintf(name,"effTrackAss");
        sprintf(title," Probability of Track Calorimeter association in  (%s)", GetName().c_str());
        h.effTrackAss = myROOT_utils::TProfile_checked (name,title,50,0.,50.);
        sprintf(name,"effPiplusAss");
        sprintf(title," Probability of Piplus Calorimeter association in  (%s)", GetName().c_str());
        h.effPiplusAss = myROOT_utils::TProfile_checked (name,title,50,0.,50.);
        sprintf(name,"effPiminusAss");
        sprintf(title," Probability of Piminus Calorimeter association in  (%s)", GetName().c_str());
        h.effPiminusAss = myROOT_utils::TProfile_checked (name,title,50,0.,50.);
        sprintf(name,"effMuplusAss");
        sprintf(title," Probability of Muplus Calorimeter association in  (%s)", GetName().c_str());
        h.effMuplusAss = myROOT_utils::TProfile_checked (name,title,50,0.,50.);
        sprintf(name,"effMuminusAss");
        sprintf(title," Probability of Muminus Calorimeter association in  (%s)", GetName().c_str());
        h.effMuminusAss = myROOT_utils::TProfile_checked (name,title,50,0.,50.);
        sprintf(name,"effElectronAss");
        sprintf(title," Probability of Electron Calorimeter association in  (%s)", GetName().c_str());
        h.effElectronAss = myROOT_utils::TProfile_checked (name,title,50,0.,50.);
        sprintf(name,"effPositronAss");
        sprintf(title," Probability of Positron Calorimeter association in  (%s)", GetName().c_str());
        h.effPositronAss = myROOT_utils::TProfile_checked (name,title,50,0.,50.);

        sprintf(name,"effCaloAss");
        sprintf(title," Probability of Calorimeter cluster association with Track  in  (%s)", GetName().c_str());
        h.effCaloAss = myROOT_utils::TProfile_checked (name,title,50,0.,50.);

        sprintf(name,"ZlastTrackAss");
        sprintf(title," Probability of Track Calorimeter association in  (%s)", GetName().c_str());
        h.ZlastTrackAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastPiminusAss");
        sprintf(title," Probability of Piminus Calorimeter association in  (%s)", GetName().c_str());
        h.ZlastPiminusAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastPiplusAss");
        sprintf(title," Probability of Piplus Calorimeter association in  (%s)", GetName().c_str());
        h.ZlastPiplusAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastMuminusAss");
        sprintf(title," Probability of Muminus Calorimeter association in  (%s)", GetName().c_str());
        h.ZlastMuminusAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastMuplusAss");
        sprintf(title," Probability of Muplus Calorimeter association in  (%s)", GetName().c_str());
        h.ZlastMuplusAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastElectronAss");
        sprintf(title," Probability of Electron Calorimeter association in  (%s)", GetName().c_str());
        h.ZlastElectronAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastPositronAss");
        sprintf(title," Probability of Positron Calorimeter association in  (%s)", GetName().c_str());
        h.ZlastPositronAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);

        sprintf(name,"ZlastTrackNotAss");
        sprintf(title," Probability of Track Calorimeter Notassociation in  (%s)", GetName().c_str());
        h.ZlastTrackNotAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastPiminusNotAss");
        sprintf(title," Probability of Piminus Calorimeter Notassociation in  (%s)", GetName().c_str());
        h.ZlastPiminusNotAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastPiplusNotAss");
        sprintf(title," Probability of Piplus Calorimeter Notassociation in  (%s)", GetName().c_str());
        h.ZlastPiplusNotAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastMuminusNotAss");
        sprintf(title," Probability of Muminus Calorimeter Notassociation in  (%s)", GetName().c_str());
        h.ZlastMuminusNotAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastMuplusNotAss");
        sprintf(title," Probability of Muplus Calorimeter Notassociation in  (%s)", GetName().c_str());
        h.ZlastMuplusNotAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastElectronNotAss");
        sprintf(title," Probability of Electron Calorimeter Notassociation in  (%s)", GetName().c_str());
        h.ZlastElectronNotAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);
        sprintf(name,"ZlastPositronNotAss");
        sprintf(title," Probability of Positron Calorimeter Notassociation in  (%s)", GetName().c_str());
        h.ZlastPositronNotAss = myROOT_utils::TH1D_checked (name,title,5000,-1000.,49000.);

//    Two near gamma problem
        sprintf(name,"%s::  Hi2 as energy function Electrons ",GetName().c_str());
        h.electronHi2OK = myROOT_utils::TH2D_checked("ElectronHi2OK", name, 50, emin,emax, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 as energy function Electrons  ",GetName().c_str());
        h.electronHi2Bad = myROOT_utils::TH2D_checked("ElectronHi2Bad", name, 50, emin,emax, 250, 0., hiMax );
        sprintf(name,"%s::  Ereal/Etot as coordinate function Electrons ",GetName().c_str());
        h.electronERNorm = myROOT_utils::TProfile_checked("ElectronERNorm", name, 100, 0, 100.);
        sprintf(name,"%s::  Eexpect/Etot as coordinate  function Electrons ",GetName().c_str());
        h.electronETNorm = myROOT_utils::TProfile_checked("ElectronETNorm", name, 100, 0, 100.);
        sprintf(name,"%s::  De/SDe as energy fraction function Electrons ",GetName().c_str());
        h.electronDENorm = myROOT_utils::TProfile_checked("ElectronDENorm", name, 50, 0, 1);
        sprintf(name,"%s::  absDe/SDe as energy fraction function Electrons ",GetName().c_str());
        h.electronDEabsNorm = myROOT_utils::TProfile_checked("ElectronDEabsNorm", name, 50, 0, 1);
        sprintf(name,"%s::  (Eexp-Ecalc)/Ecalc as energy fraction function Electrons ",GetName().c_str());
        h.electronDEcalc = myROOT_utils::TProfile_checked("ElectronDEcalc", name, 50, 0, 1);
        sprintf(name,"%s::  Himax as energy fraction function Electrons ",GetName().c_str());
        h.electronHiMaxOK  = myROOT_utils::TH2D_checked("ElectronHimaxOK", name, 50, 0, 1, 250, 0., hiMax);
        sprintf(name,"%s::  Hi-Himax vs Himax  Electrons",GetName().c_str());
        h.electronHiHiMaxOK = myROOT_utils::TH2D_checked("ElectronHiHimaxOK", name, 250, 0., hiMax, 250, 0., hiMax);
        sprintf(name,"%s::  Himax as energy fraction function Electrons ",GetName().c_str());
        h.electronHiMaxBad  = myROOT_utils::TH2D_checked("ElectronHimaxBad", name, 50, 0, 1, 250, 0., hiMax);
        sprintf(name,"%s::  Hi-Himax vs Himax Electrons ",GetName().c_str());
        h.electronHiHiMaxBad = myROOT_utils::TH2D_checked("ElectronHiHimaxBad", name, 250, 0., hiMax, 250, 0., hiMax);

        sprintf(name,"%s::  Hi2 Gammas OK  ",GetName().c_str());
        h.gammaHi2OK = myROOT_utils::TH1D_checked("GammaHi2OK", name, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 Gammas Bad ",GetName().c_str());
        h.gammaHi2Bad = myROOT_utils::TH1D_checked("GammaHi2Bad", name, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 as energy function Gammas ",GetName().c_str());
        h.gammaHi2EOK = myROOT_utils::TH2D_checked("GammaHi2EOK", name, 50, emin,emax, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 as energy function Gammas  ",GetName().c_str());
        h.gammaHi2EBad = myROOT_utils::TH2D_checked("GammaHi2EBad", name, 50, emin,emax, 250, 0., hiMax );

        sprintf(name,"%s::  Hi2 Gammas_1gEvt OK  ",GetName().c_str());
        h.gammaHi2OK_1gEvt = myROOT_utils::TH1D_checked("GammaHi2OK_1gEvt", name, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 Gammas_1gEvt Bad ",GetName().c_str());
        h.gammaHi2Bad_1gEvt = myROOT_utils::TH1D_checked("GammaHi2Bad_1gEvt", name, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 as energy function Gammas_1gEvt ",GetName().c_str());
        h.gammaHi2EOK_1gEvt = myROOT_utils::TH2D_checked("GammaHi2EOK_1gEvt", name, 50, emin,emax, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 as energy function Gammas_1gEvt  ",GetName().c_str());
        h.gammaHi2EBad_1gEvt = myROOT_utils::TH2D_checked("GammaHi2EBad_1gEvt", name, 50, emin,emax, 250, 0., hiMax );

        sprintf(name,"%s::  Hi2 Gammas_2gEvt OK  ",GetName().c_str());
        h.gammaHi2OK_2gEvt = myROOT_utils::TH1D_checked("GammaHi2OK_2gEvt", name, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 Gammas_2gEvt Bad ",GetName().c_str());
        h.gammaHi2Bad_2gEvt = myROOT_utils::TH1D_checked("GammaHi2Bad_2gEvt", name, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 as energy function Gammas_2gEvt ",GetName().c_str());
        h.gammaHi2EOK_2gEvt = myROOT_utils::TH2D_checked("GammaHi2EOK_2gEvt", name, 50, emin,emax, 250, 0., hiMax );
        sprintf(name,"%s::  Hi2 as energy function Gammas_2gEvt  ",GetName().c_str());
        h.gammaHi2EBad_2gEvt = myROOT_utils::TH2D_checked("GammaHi2EBad_2gEvt", name, 50, emin,emax, 250, 0., hiMax );


        sprintf(name,"%s::  Ereal/Etot as coordinate function Gammas ",GetName().c_str());
        h.gammaERNorm = myROOT_utils::TProfile_checked("GammaERNorm", name, 100, emin, emax );
        sprintf(name,"%s::  Eexpect/Etot as coordinate  function Gammas ",GetName().c_str());
        h.gammaETNorm = myROOT_utils::TProfile_checked("GammaETNorm", name, 100, 0, 100.);
        sprintf(name,"%s::  De/SDe as energy fraction function Gammas ",GetName().c_str());
        h.gammaDENorm = myROOT_utils::TProfile_checked("GammaDENorm", name, 500, 0, 1);
        sprintf(name,"%s::  absDe/SDe as energy fraction function Gammas ",GetName().c_str());
        h.gammaDEabsNorm = myROOT_utils::TProfile_checked("GammaDEabsNorm", name, 250, 0, 1);
        sprintf(name,"%s::  (Eexp-Ecalc)/Ecalc as energy fraction function Gammas ",GetName().c_str());
        h.gammaDEcalc = myROOT_utils::TProfile_checked("GammaDEcalc", name, 500, 0, 1);

        sprintf(name,"%s::  Himax as energy fraction function Gammas OK ",GetName().c_str());
        h.gammaHiMaxOK  = myROOT_utils::TH2D_checked("GammaHimaxOK", name, 500, 0, 1, 250, 0., hiMax);
        sprintf(name,"%s::  Hi-Himax vs Himax  Gammas OK ",GetName().c_str());
        h.gammaHiHiMaxOK = myROOT_utils::TH2D_checked("GammaHiHimaxOK", name, 250, 0., hiMax, 250, 0., hiMax);
        sprintf(name,"%s::  Himax as energy fraction function Gammas Bad ",GetName().c_str());
        h.gammaHiMaxBad  = myROOT_utils::TH2D_checked("GammaHimaxBad", name, 500, 0, 1, 250, 0., hiMax);
        sprintf(name,"%s::  Hi-Himax vs Himax Gammas Bad ",GetName().c_str());
        h.gammaHiHiMaxBad = myROOT_utils::TH2D_checked("GammaHiHimaxBad", name, 250, 0., hiMax, 250, 0., hiMax);

        sprintf(name,"%s::  HiSign as energy fraction function Gammas OK ",GetName().c_str());
        h.gammaHiSignOK  = myROOT_utils::TH2D_checked("GammaHiSignOK", name, 500, 0, 1, 250, -hiMax, hiMax);
        sprintf(name,"%s::  Hi vs HiSign  Gammas OK ",GetName().c_str());
        h.gammaHiHiSignOK = myROOT_utils::TH2D_checked("GammaHiHiSignOK", name, 250, 0., hiMax, 250, -hiMax, hiMax);
        sprintf(name,"%s::  HiSign as energy fraction function Gammas Bad ",GetName().c_str());
        h.gammaHiSignBad  = myROOT_utils::TH2D_checked("GammaHiSignBad", name, 50, 0, 1, 250, -hiMax, hiMax);
        sprintf(name,"%s::  Hi vs HiSign Gammas Bad ",GetName().c_str());
        h.gammaHiHiSignBad = myROOT_utils::TH2D_checked("GammaHiHiSignBad", name, 250, 0., hiMax, 250, -hiMax, hiMax);

        sprintf(name,"%s::  Himax as energy fraction function Gammas_1gEvt OK ",GetName().c_str());
        h.gammaHiMaxOK_1gEvt  = myROOT_utils::TH2D_checked("GammaHimaxOK_1gEvt", name, 50, 0, 1, 250, 0., hiMax);
        sprintf(name,"%s::  Hi-Himax vs Himax  Gammas_1gEvt OK ",GetName().c_str());
        h.gammaHiHiMaxOK_1gEvt = myROOT_utils::TH2D_checked("GammaHiHimaxOK_1gEvt", name, 250, 0., hiMax, 250, 0., hiMax);
        sprintf(name,"%s::  Himax as energy fraction function Gammas_1gEvt Bad ",GetName().c_str());
        h.gammaHiMaxBad_1gEvt  = myROOT_utils::TH2D_checked("GammaHimaxBad_1gEvt", name, 50, 0, 1, 250, 0., hiMax);
        sprintf(name,"%s::  Hi-Himax vs Himax Gammas_1gEvt Bad ",GetName().c_str());
        h.gammaHiHiMaxBad_1gEvt = myROOT_utils::TH2D_checked("GammaHiHimaxBad_1gEvt", name, 250, 0., hiMax, 250, 0., hiMax);

        sprintf(name,"%s::  HiSign as energy fraction function Gammas_1gEvt OK ",GetName().c_str());
        h.gammaHiSignOK_1gEvt  = myROOT_utils::TH2D_checked("GammaHiSignOK_1gEvt", name, 50, 0, 1, 250, -hiMax, hiMax);
        sprintf(name,"%s::  Hi vs HiSign  Gammas_1gEvt OK ",GetName().c_str());
        h.gammaHiHiSignOK_1gEvt = myROOT_utils::TH2D_checked("GammaHiHiSignOK_1gEvt", name, 250, 0., hiMax, 250, -hiMax, hiMax);
        sprintf(name,"%s::  HiSign as energy fraction function Gammas_1gEvt Bad ",GetName().c_str());
        h.gammaHiSignBad_1gEvt  = myROOT_utils::TH2D_checked("GammaHiSignBad_1gEvt", name, 50, 0, 1, 250, -hiMax, hiMax);
        sprintf(name,"%s::  Hi vs HiSign Gammas_1gEvt Bad ",GetName().c_str());
        h.gammaHiHiSignBad_1gEvt = myROOT_utils::TH2D_checked("GammaHiHiSignBad_1gEvt", name, 250, 0., hiMax, 250, -hiMax, hiMax);

        sprintf(name,"%s::  Himax as energy fraction function Gammas_2gEvt OK ",GetName().c_str());
        h.gammaHiMaxOK_2gEvt  = myROOT_utils::TH2D_checked("GammaHimaxOK_2gEvt", name, 50, 0, 1, 250, 0., hiMax);
        sprintf(name,"%s::  Hi-Himax vs Himax  Gammas_2gEvt OK ",GetName().c_str());
        h.gammaHiHiMaxOK_2gEvt = myROOT_utils::TH2D_checked("GammaHiHimaxOK_2gEvt", name, 250, 0., hiMax, 250, 0., hiMax);
        sprintf(name,"%s::  Himax as energy fraction function Gammas_2gEvt Bad ",GetName().c_str());
        h.gammaHiMaxBad_2gEvt  = myROOT_utils::TH2D_checked("GammaHimaxBad_2gEvt", name, 50, 0, 1, 250, 0., hiMax);
        sprintf(name,"%s::  Hi-Himax vs Himax Gammas_2gEvt Bad ",GetName().c_str());
        h.gammaHiHiMaxBad_2gEvt = myROOT_utils::TH2D_checked("GammaHiHimaxBad_2gEvt", name, 250, 0., hiMax, 250, 0., hiMax);

        sprintf(name,"%s::  HiSign as energy fraction function Gammas_2gEvt OK ",GetName().c_str());
        h.gammaHiSignOK_2gEvt  = myROOT_utils::TH2D_checked("GammaHiSignOK_2gEvt", name, 50, 0, 1, 250, -hiMax, hiMax);
        sprintf(name,"%s::  Hi vs HiSign  Gammas_2gEvt OK ",GetName().c_str());
        h.gammaHiHiSignOK_2gEvt = myROOT_utils::TH2D_checked("GammaHiHiSignOK_2gEvt", name, 250, 0., hiMax, 250, -hiMax, hiMax);
        sprintf(name,"%s::  HiSign as energy fraction function Gammas_2gEvt Bad ",GetName().c_str());
        h.gammaHiSignBad_2gEvt  = myROOT_utils::TH2D_checked("GammaHiSignBad_2gEvt", name, 50, 0, 1, 250, -hiMax, hiMax);
        sprintf(name,"%s::  Hi vs HiSign Gammas_2gEvt Bad ",GetName().c_str());
        h.gammaHiHiSignBad_2gEvt = myROOT_utils::TH2D_checked("GammaHiHiSignBad_2gEvt", name, 250, 0., hiMax, 250, -hiMax, hiMax);

        sprintf(name,"%s::  Dist vs emin/emax 2 gamma e1+e2>5 GeV",GetName().c_str());
        h.h2_Dist2 = myROOT_utils::TH2D_checked("Dist2", name, 200, 0, 400, 50, 0, 1);
        sprintf(name,"%s::  Dist vs emin e1+e2>5 GeV",GetName().c_str());
        h.h2_Dist2emin = myROOT_utils::TH2D_checked("Dist2emin", name, 200, 0, 400, 50, 0, 10);
        sprintf(name,"%s::  DistXY 2  ",GetName().c_str());
        h.h2_Dist2Cell = myROOT_utils::TH2D_checked("Dist2Cell", name, 100, -400, 400, 100, -400, 400);

      if( debug ) cout << " MakeParticlesHisto Book for " <<  GetName() << " OK " <<endl;
    }
    h.booked_ok_ = true;

  }
  MakeParticlesHisto &h = *make_particles_histo_;  // create a short name for histo store
  if( !h.booked_ok_ )
  {
    cerr << " FATAL ERROR!! Calorimeter::FillParticlesCaloAssociationHist booking problems ... memory?? " << endl;
    exit(1);
  }

  if( debug ) cout << " Calorimeter::FillParticlesCaloAssociationHist Ntot = " <<
                                                    hint_particles.size() << endl;

  if( debug ) cout << " Nparticles " <<   hint_particles.size() << endl;

  if( (int)reco_particles.size() > options.ngam_max4histo ) return;


  double DistNCellSize = 1.5;
  vector < vector < int > > hpart_tr_ass( hint_particles.size() );
  vector < vector < int > > rpart_tr_ass( reco_particles.size() );

  for(size_t ip =0; ip<hint_particles.size(); ip++)                           //   loop over Particles
  {
    const CalorimeterParticle &p = hint_particles[ip];
    if( debug ) cout << " Particle " << ip << endl;
    if( p.GetID() == CalorimeterParticle::GAMMA || p.GetID() == CalorimeterParticle::KAON_0_LONG ) continue;
      double xtrk0 = p.GetX();
      double ytrk0 = p.GetY();
      double zcal0 = p.GetZ();
      double momentum = p.GetE();
      if( fabs(momentum - 9999.) < 0.0001 ) momentum = 0.; // Momentum not defined
      double track_time = p.GetTime();

      bool track_associated = false;

//       for( int icc=0; icc < calo_reco_particles_[it_calo].size(); icc++)                   // loop over all clusters
      for(size_t icc=0; icc<reco_particles.size(); icc++)                   // loop over all clusters
      {
         const CalorimeterParticle &pe = reco_particles[icc];
//         const CalorimeterParticle &pe =calo_reco_particles_[it_calo][icc];
//         int cell_main = pe.GetMainCells()[0];
//         int x =  GetColumnOfCell(cell_main);
//         int y =  GetRowOfCell(cell_main);
//         double ec = pe.GetE();
//
//         double eg = pe.GetE();
        double xcalc = pe.GetX();
        double ycalc = pe.GetY();
//        double zcalc = pe.GetZ();
        if(levelBPHistos > 1)
        {
            h.deltaXdeltaYall->Fill(xcalc-xtrk0,ycalc-ytrk0);
        }
      }


//
//       const vector < pair < pair <int,int>, TrackPar> >  &calorim_associated = p.GetCaloAssociated();
//       if( debug ) cout << " calorim_associated size " << calorim_associated.size() << endl;
//       for( int ic=0; ic < calorim_associated.size(); ic++)                   // loop over associated clusters
//       {
      for(size_t icc=0; icc<reco_particles.size(); icc++)                   // loop over all clusters
      {
         const CalorimeterParticle &pe = reco_particles[icc];
//         if( it_calo != calorim_associated[ic].first.first ) continue;
//         track_associated = true;
//         int icc = calorim_associated[ic].first.second;
//         const CalorimeterParticle &pe =calo_reco_particles_[it_calo][icc];
        double eg = pe.GetE();
        double xcalc = pe.GetX();
        double ycalc = pe.GetY();
        double zcalc = pe.GetZ();

        double xmain = -100000.;
        double ymain = -100000.;
        double emain = 0.;
        double e4 = 0.;
        double esumm = 0.;

        bool cluster_near_bad_cell = false;
        int cell_main = pe.GetMainCells()[0];
        double cell_size_x = cells[cell_main].GetCellType().GetSizeX();
        double cell_size_y = cells[cell_main].GetCellType().GetSizeY();
        xmain =  GetCells()[cell_main].GetX();
        ymain =  GetCells()[cell_main].GetY();
        if(  BadCellsAround(cell_main).size() > 0 ) cluster_near_bad_cell=true;

        if( fabs( xcalc-xtrk0 ) > DistNCellSize*cell_size_x ) continue;
        if( fabs( ycalc-ytrk0 ) > DistNCellSize*cell_size_y ) continue;

        track_associated = true;

        rpart_tr_ass[icc].push_back(ip);
        hpart_tr_ass[ip].push_back(icc);

//         cout << GetName() << " dX " << xcalc-xtrk0 << " dY " << ycalc-ytrk0 << endl;

        double excal = pe.GetXerr();
        double eycal = pe.GetYerr();
        double xtrk = xtrk0;
        double ytrk = ytrk0;
//         const double *covtrk = trackpar_at_destination.GetCov();
        if(fabs(zcal0-zcalc) > 0.1 ) // OK precision first!
        {
//            TrackPar trackpar_precise;
// //           ExtrapolateExternal(trackpar_at_destination,zcalc, trackpar_precise, false);
//            Extrapolate(trackpar_at_destination,zcalc, trackpar_precise, false);
//            xtrk = trackpar_precise.GetX();
//            ytrk = trackpar_precise.GetY();
//            covtrk = trackpar_precise.GetCov();
        }

        double dxcelltrk = xtrk-xmain;
        double dycelltrk = ytrk-ymain;
        double dxcellcalc = xcalc-xmain;
        double dycellcalc = ycalc-ymain;
        double dxct = xcalc-xtrk;
        double dyct = ycalc-ytrk;
//         cout << " ecalo " << pe.GetE() << endl;
//         cout << " xtrk " << xtrk << " xcalc " << xcalc << " xmain " << xmain << endl;
//         cout << " ytrk " << ytrk << " ycalc " << ycalc << " ymain " << ymain << endl;
//         cout << " dxcelltrk " << dxcelltrk << " dxcellcalc " << dxcellcalc << endl;
//         cout << " dycelltrk " << dycelltrk << " dycellcalc " << dycellcalc << endl;


        h.deltaX->Fill(xcalc-xtrk);
        h.deltaXS->Fill((xcalc-xtrk)/excal);
        h.deltaY->Fill(ycalc-ytrk);
        h.deltaYS->Fill((ycalc-ytrk)/eycal);
        if( p.GetID() == CalorimeterParticle::ELECTRON || p.GetID() == CalorimeterParticle::POSITRON )
        {
          if( debug ) cout << " Filling histo for electrons " << endl;
          bool good_zone4electrons_selected = GoodZone4Electrons( xtrk, ytrk, eg );

          if( good_zone4electrons_selected )
          {
            h.deltaXelectron_InGoodZone->Fill(xcalc-xtrk);
            h.deltaYelectron_InGoodZone->Fill(ycalc-ytrk);
            h.deltaEsElectron_InGoodZone->Fill(eg/momentum-1.);
          }


          h.electronE->Fill(eg);

          h.deltaXelectron->Fill(xcalc-xtrk);
          h.deltaYelectron->Fill(ycalc-ytrk);

          if( cluster_near_bad_cell )
          {
            h.deltaXelectron_NearBadCell->Fill(xcalc-xtrk);
            h.deltaYelectron_NearBadCell->Fill(ycalc-ytrk);
          }
          if(levelBPHistos > 1)
          {
            h.deltaXelectronX->Fill(xcalc-xtrk,xtrk);
            h.deltaYelectronY->Fill(ycalc-ytrk,ytrk);
            h.MeanDeltaXelectronX->Fill(xtrk,xcalc-xtrk);
            h.MeanDeltaYelectronY->Fill(ytrk,ycalc-ytrk);

            h.caloXelectronXcell->Fill(dxcellcalc,dxcelltrk);
            h.caloYelectronYcell->Fill(dycellcalc,dycelltrk);
            h.deltaXelectronXcell->Fill(dxct,dxcelltrk);
            h.deltaYelectronYcell->Fill(dyct,dycelltrk);
//             h.deltaSqXelectronSqX->Fill((xcalc-xtrk)*(xcalc-xtrk),covtrk[0]);
//             h.deltaSqYelectronSqY->Fill((ycalc-ytrk)*(ycalc-ytrk),covtrk[1]);
            h.deltaXelectronE->Fill(xcalc-xtrk,eg);
            h.deltaYelectronE->Fill(ycalc-ytrk,eg);
            if( good_zone4electrons_selected )
            {
              h.deltaXelectronX_selected->Fill(xcalc-xtrk,xtrk);
              h.deltaYelectronY_selected->Fill(ycalc-ytrk,ytrk);
              h.MeanDeltaXelectronX_selected->Fill(xtrk,xcalc-xtrk);
              h.MeanDeltaYelectronY_selected->Fill(ytrk,ycalc-ytrk);
              h.deltaXelectronXcell_selected->Fill(dxct,dxcelltrk);
              h.deltaYelectronYcell_selected->Fill(dyct,dycelltrk);
            }
          }
        }

        if( p.GetID() == CalorimeterParticle::PION_MINUS || p.GetID() == CalorimeterParticle::PION_PLUS )
        {
          if( debug ) cout << " Filling histo for pions " << endl;
          h.pionE->Fill(eg);
          h.deltaXpion->Fill(xcalc-xtrk);
          h.deltaYpion->Fill(ycalc-ytrk);
          if(levelBPHistos > 1)
          {
            h.deltaXpionX->Fill(xcalc-xtrk,xtrk);
            h.deltaYpionY->Fill(ycalc-ytrk,ytrk);
            h.caloXpionXcell->Fill(dxcellcalc,dxcelltrk);
            h.caloYpionYcell->Fill(dycellcalc,dycelltrk);
            h.deltaXpionXcell->Fill(dxct,dxcelltrk);
            h.deltaYpionYcell->Fill(dyct,dycelltrk);
//             h.deltaSqXpionSqX->Fill((xcalc-xtrk)*(xcalc-xtrk),covtrk[0]);
//             h.deltaSqYpionSqY->Fill((ycalc-ytrk)*(ycalc-ytrk),covtrk[2]);
            h.deltaXpionE->Fill(xcalc-xtrk,eg);
            h.deltaYpionE->Fill(ycalc-ytrk,eg);
          }
        }

        if( p.GetID() == CalorimeterParticle::MUON_MINUS || p.GetID() == CalorimeterParticle::MUON_PLUS )
        {
          if( debug ) cout << " Filling histo for muons " << endl;
          h.muonE->Fill(eg);
          h.deltaXmuon->Fill(xcalc-xtrk);
          h.deltaYmuon->Fill(ycalc-ytrk);
          if(levelBPHistos > 1)
          {
            h.deltaXmuonX->Fill(xcalc-xtrk,xtrk);
            h.deltaYmuonY->Fill(ycalc-ytrk,ytrk);
            h.caloXmuonXcell->Fill(dxcellcalc,dxcelltrk);
            h.caloYmuonYcell->Fill(dycellcalc,dycelltrk);
            h.deltaXmuonXcell->Fill(dxct,dxcelltrk);
            h.deltaYmuonYcell->Fill(dyct,dycelltrk);
//             h.deltaSqXmuonSqX->Fill((xcalc-xtrk)*(xcalc-xtrk),covtrk[0]);
//             h.deltaSqYmuonSqY->Fill((ycalc-ytrk)*(ycalc-ytrk),covtrk[2]);
            h.deltaXmuonE->Fill(xcalc-xtrk,eg);
            h.deltaYmuonE->Fill(ycalc-ytrk,eg);
          }
        }
        if(levelBPHistos > 1)
        {
          h.caloXtrackX->Fill(xcalc,xtrk);
          h.caloYtrackY->Fill(ycalc,ytrk);
          h.caloXtrackXcell->Fill(dxcellcalc,dxcelltrk);
          h.caloYtrackYcell->Fill(dycellcalc,dycelltrk);
          h.deltaXtrackXcell->Fill(dxct,dxcelltrk);
          h.deltaYtrackYcell->Fill(dyct,dycelltrk);
          h.deltaXtrackX->Fill(xcalc-xtrk,xtrk);
          h.deltaYtrackY->Fill(ycalc-ytrk,ytrk);
          h.deltaXdeltaY->Fill(xcalc-xtrk,ycalc-ytrk);
//           h.deltaSqXtrackSqX->Fill((xcalc-xtrk)*(xcalc-xtrk),covtrk[0]);
//           h.deltaSqYtrackSqY->Fill((ycalc-ytrk)*(ycalc-ytrk),covtrk[2]);

          h.caloEtrackE->Fill(eg,momentum);
          h.deltaEtrackE->Fill(eg-momentum,momentum);

          if( p.GetID() == CalorimeterParticle::ELECTRON || p.GetID() == CalorimeterParticle::POSITRON )
          {
            if( debug ) cout << " Filling 2D histo for electrons " << endl;
            h.caloEelectronE->Fill(eg,momentum);
            if(eg < 6 ) h.electronEmE4low->Fill(e4-emain,1.-e4);
            if(eg >= 6 ) h.electronEmE4high->Fill(e4-emain,1.-e4);
          }

          if( p.GetID() == CalorimeterParticle::PION_MINUS || p.GetID() == CalorimeterParticle::PION_PLUS )
          {
            if( debug ) cout << " Filling 2D histo for pions " << endl;
            h.caloEpionE->Fill(eg,momentum);
            if(eg < 6 ) h.pionEmE4low->Fill(e4-emain,1.-e4);
            if(eg >= 6 ) h.pionEmE4high->Fill(e4-emain,1.-e4);
          }

          if( p.GetID() == CalorimeterParticle::MUON_MINUS || p.GetID() == CalorimeterParticle::MUON_PLUS )
          {
            if( debug ) cout << " Filling 2D histo for muons " << endl;
            h.caloEmuonE->Fill(eg,momentum);
            if(eg < 6 ) h.muonEmE4low->Fill(e4-emain,1.-e4);
            if(eg >= 6 ) h.muonEmE4high->Fill(e4-emain,1.-e4);
          }

        }

//        if( momentum > 5.0 && fabs(eg - momentum) < 10.*sqrt(momentum)*0.07 ) // associated with energy ( calorimeter default particle )
        if( momentum > 2.0  ) // associated with energy ( calorimeter default particle )
        {
          if( debug ) cout << " Filling histo for tracks with momentum " << endl;
          h.deltaEs->Fill( eg/momentum - 1. );
          h.deltaXwEas->Fill(xcalc-xtrk);
          h.deltaYwEas->Fill(ycalc-ytrk);
          if( p.GetID() == CalorimeterParticle::ELECTRON || p.GetID() == CalorimeterParticle::POSITRON )
          {
            if( debug ) cout << " Filling histo for electrons with momentum " << endl;
            h.deltaEselectron->Fill( eg/momentum - 1. );
            if(levelBPHistos > 1)
            {
              h.MeanDeltaEsElectronX->Fill( xtrk, eg/momentum - 1. );
              h.MeanDeltaEsElectronY->Fill( ytrk, eg/momentum - 1. );
//             if(  XYRegularGrid() )
//             {
//               int ncolumn =  GetColumnOfCell( cell_main );
//               int nrow =  GetRowOfCell( cell_main );
//               h.MeanDeltaEsElectronXRow[it_calo][nrow]->Fill( xtrk, eg/momentum - 1. );
//               h.MeanDeltaEsElectronYColumn[it_calo][ncolumn]->Fill( ytrk, eg/momentum - 1. );
//             }
            }

            if( cluster_near_bad_cell )
            {
              h.deltaEsElectron_NearBadCell->Fill( eg/momentum - 1.);
            }
          }

          if( p.GetID() == CalorimeterParticle::PION_MINUS || p.GetID() == CalorimeterParticle::PION_PLUS )
          {
            if( debug ) cout << " Filling histo for pions with momentum " << endl;
            h.deltaEspion->Fill( eg/momentum - 1. );
          }

          if( p.GetID() == CalorimeterParticle::MUON_MINUS || p.GetID() == CalorimeterParticle::MUON_PLUS )
          {
            if( debug ) cout << " Filling histo for muons with momentum " << endl;
            h.deltaEsmuon->Fill( eg/momentum - 1. );
          }
          if(levelBPHistos > 1)
          {
            if( debug ) cout << " Filling 2D histo for tracks with momentum " << endl;
            h.deltaEsEtrack->Fill( eg/momentum - 1., momentum );
            h.deltaXtrackXwEas->Fill(xcalc-xtrk,xtrk);
            h.deltaYtrackYwEas->Fill(ycalc-ytrk,ytrk);
            if( p.GetID() == CalorimeterParticle::ELECTRON || p.GetID() == CalorimeterParticle::POSITRON )
            {
              if( debug ) cout << " Filling 2D histo for electrons with momentum " << endl;
              h.deltaEsEelectron->Fill( eg/momentum - 1.,momentum );
            }

            if( p.GetID() == CalorimeterParticle::PION_MINUS || p.GetID() == CalorimeterParticle::PION_PLUS )
            {
              if( debug ) cout << " Filling 2D histo for pions with momentum " << endl;
              h.deltaEsEpion->Fill( eg/momentum - 1.,momentum );
            }

            if( p.GetID() == CalorimeterParticle::MUON_MINUS || p.GetID() == CalorimeterParticle::MUON_PLUS )
            {
              if( debug ) cout << " Filling 2D histo for muons with momentum " << endl;
              h.deltaEsEmuon->Fill( eg/momentum - 1.,momentum );
            }
          }
        }
        if(pe.HasTime()  )
        {
          if( debug ) cout << " Filling timing histo  " << endl;
          if(levelBPHistos > 1) h.caloTtrackT->Fill(pe.GetTime(),track_time );
        }
        if( debug ) cout << " GOto next cluster  " << endl;

      }                                                                  // end of loop over associated clusters


      if( debug ) cout << " Filling XY associate/not histo  " << endl;
      if(levelBPHistos > 1)
      {
        h.trackXY->Fill(xtrk0,ytrk0);

        if( p.GetID() == CalorimeterParticle::ELECTRON || p.GetID() == CalorimeterParticle::POSITRON )
        {
           if( debug ) cout << " Filling XY associate/not histo for electrons " << endl;
           h.electronXY->Fill(xtrk0,ytrk0);
           if(track_associated)
             h.electronXYas->Fill(xtrk0,ytrk0);
           else
             h.electronXYnotas->Fill(xtrk0,ytrk0);
        }
        if( p.GetID() == CalorimeterParticle::PION_PLUS || p.GetID() == CalorimeterParticle::PION_MINUS )
        {
           if( debug ) cout << " Filling XY associate/not histo for pions " << endl;
           h.pionXY->Fill(xtrk0,ytrk0);
           if(track_associated)
             h.pionXYas->Fill(xtrk0,ytrk0);
           else
             h.pionXYnotas->Fill(xtrk0,ytrk0);
        }
        if( p.GetID() == CalorimeterParticle::MUON_PLUS || p.GetID() == CalorimeterParticle::MUON_MINUS )
        {
           if( debug ) cout << " Filling XY associate/not histo for muons " << endl;
           if( debug ) cout << " muon xtrk0 " << xtrk0 << " ytrk0 " << ytrk0 << endl;
           h.muonXY->Fill(xtrk0,ytrk0);
           if( debug ) cout << " xtrk0 " << xtrk0 << " ytrk0 " << ytrk0 << endl;
           if( debug ) cout << "  track_associated " << track_associated << endl;

           if(track_associated)
	   {
             if( debug ) cout << "  track associated " << endl;
             h.muonXYas->Fill(xtrk0,ytrk0);
	   }
           else
	   {
             if( debug ) cout << "  track not associated " << endl;
             h.muonXYnotas->Fill(xtrk0,ytrk0);
	   }
        }
      }

      if( debug ) cout << " Continue Filling XY associate/not histo for all tracks " << endl;
      if(track_associated)
      {
        if(levelBPHistos > 1) h.trackXYas->Fill(xtrk0,ytrk0);
        if(momentum > 0.01 )  {
          if(levelBPHistos > 1) h.trackXYas_E->Fill(xtrk0,ytrk0);
        }
        else  {
          if(levelBPHistos > 1) h.trackXYas_noE->Fill(xtrk0,ytrk0);
        }
      }
      else
      {
        if(levelBPHistos > 1) h.trackXYnotas->Fill(xtrk0,ytrk0);
        if(momentum > 0.01 )  {
          if(levelBPHistos > 1) h.trackXYnotas_E->Fill(xtrk0,ytrk0);
        }
        else  {
          if(levelBPHistos > 1) h.trackXYnotas_noE->Fill(xtrk0,ytrk0);
        }
      }

      if( debug ) cout << " Filling association efficiency for all tracks " << endl;
      h.effTrackAss->Fill(momentum,track_associated);
      if( p.GetID() == CalorimeterParticle::ELECTRON )
      {
        h.effElectronAss->Fill(momentum,track_associated);
//         if( track_associated )
//           h.ZlastElectronAss->Fill(zlast);
//         else
//           h.ZlastElectronNotAss->Fill(zlast);
      }
      else if( p.GetID() == CalorimeterParticle::POSITRON )
      {
        h.effPositronAss->Fill(momentum,track_associated);
//         if( track_associated )
//           h.ZlastPositronAss->Fill(zlast);
//         else
//           h.ZlastPositronNotAss->Fill(zlast);
      }
      else if( p.GetID() == CalorimeterParticle::PION_PLUS )
      {
        h.effPiplusAss->Fill(momentum,track_associated);
//         if( track_associated )
//           h.ZlastPiplusAss->Fill(zlast);
//         else
//           h.ZlastPiplusNotAss->Fill(zlast);
      }
      else if( p.GetID() == CalorimeterParticle::PION_MINUS )
      {
        h.effPiminusAss->Fill(momentum,track_associated);
//         if( track_associated )
//           h.ZlastPiminusAss->Fill(zlast);
//         else
//           h.ZlastPiminusNotAss->Fill(zlast);
      }
      else if( p.GetID() == CalorimeterParticle::MUON_PLUS )
      {
        h.effMuplusAss->Fill(momentum,track_associated);
//         if( track_associated )
//           h.ZlastMuplusAss->Fill(zlast);
//         else
//           h.ZlastMuplusNotAss->Fill(zlast);
      }
      else if( p.GetID() == CalorimeterParticle::MUON_MINUS )
      {
        h.effMuminusAss->Fill(momentum,track_associated);
//         if( track_associated )
//           h.ZlastMuminusAss->Fill(zlast);
//         else
//           h.ZlastMuminusNotAss->Fill(zlast);
      }

    if( debug ) cout << " Go to next particle " << endl;
  }                                                                      // end of loop over Particles
//
//
// Pure Calorimeters Histograms
  if( debug ) cout << " Pure Calorimeters Histograms " << endl;

//  bool associate[reco_particles.size()];
  double sume = 0.;
  bool development = false;
  if( is_ECAL2 ) development = false;
  for(size_t ir =0; ir<reco_particles.size(); ir++ )
  {
    const CalorimeterParticle &pe = reco_particles[ir];
//    associate[ir]=false;
    double xcalc = pe.GetX();
    double ycalc = pe.GetY();
    double e = pe.GetE();

    sume += e;
    h.caloE->Fill( e );
    if(levelBPHistos > 1)
      h.caloXY->Fill(xcalc,ycalc);

  }

  unsigned int countGamma(0);
  unsigned int countKl(0);
  unsigned int countPiplus(0);
  unsigned int countPiminus(0);
  unsigned int countMuplus(0);
  unsigned int countMuminus(0);
  unsigned int countElectrons(0);
  unsigned int countPositrons(0);

  // loop over particles
  for (size_t ipp=0; ipp<hint_particles.size(); ipp++) {
    const CalorimeterParticle& pip = hint_particles[ipp];

    if      ( pip.GetID() == CalorimeterParticle::GAMMA       ) countGamma++;
    else if ( pip.GetID() == CalorimeterParticle::POSITRON    ) countPositrons++;
    else if ( pip.GetID() == CalorimeterParticle::ELECTRON    ) countElectrons++;
    else if ( pip.GetID() == CalorimeterParticle::MUON_PLUS   ) countMuplus++;
    else if ( pip.GetID() == CalorimeterParticle::MUON_MINUS  ) countMuminus++;
    else if ( pip.GetID() == CalorimeterParticle::KAON_0_LONG ) countKl++;
    else if ( pip.GetID() == CalorimeterParticle::PION_PLUS   ) countPiplus++;
    else if ( pip.GetID() == CalorimeterParticle::PION_MINUS  ) countPiminus++;
  }

  if( debug ) cout << " Fill Muliplicity Histo " << endl;
  h.caloSumE->Fill  ( sume );
  h.Ngamma->Fill    ( countGamma );
  h.NKl->Fill       ( countKl );
  h.NPiplus->Fill   ( countPiplus );
  h.NPiminus->Fill  ( countPiminus );
  h.NMuplus->Fill   ( countMuplus );
  h.NMuminus->Fill  ( countMuminus );
  h.NElectrons->Fill( countElectrons );
  h.NPositrons->Fill( countPositrons );

  if( debug_gamma ) cout << " NKl = " << countKl << endl;
  for (size_t ip=0; ip<hint_particles.size(); ip++)                           //   loop over Particles
  {
    const CalorimeterParticle &p = hint_particles[ip];
    if( debug ) cout << " Particle " << ip << endl;
    if( !( p.GetID() == CalorimeterParticle::GAMMA || p.GetID() == CalorimeterParticle::KAON_0_LONG ) ) continue;
    if( (p.GetID() == CalorimeterParticle::GAMMA) && (GetType() == Calorimeter::Hadronic) ) continue;
    if( (p.GetID() == CalorimeterParticle::KAON_0_LONG) && (GetType() == Calorimeter::ElectroMagnetic) ) continue;

    double xtrk0 = p.GetX();
    double ytrk0 = p.GetY();
    double zcal0 = p.GetZ();
    double momentum = p.GetE();
//    double track_time = p.GetTime();
    if(p.GetID() == CalorimeterParticle::KAON_0_LONG )
    {
      h.neutralE->Fill(momentum);
      if(levelBPHistos > 1)
      {
        h.neutralXY->Fill( xtrk0, ytrk0 );
      }
    }

    if(p.GetID() == CalorimeterParticle::GAMMA )
    {
      h.gammaE->Fill(momentum);
      if(levelBPHistos > 1)
      {
        h.gammaXY->Fill( xtrk0, ytrk0 );
      }
    }
  }

  if( debug ) cout << " Gamma-gamma Ngamma = " << countGamma << endl;
  if( debug ) cout << " NElectrons = "         << countElectrons << endl;
  if( debug ) cout << " NPositrons = "         << countPositrons << endl;
  if( debug ) cout << " NPiplus = "            << countPiplus << endl;
  if( debug ) cout << " NPiminus = "           << countPiminus << endl;

  unsigned int n_gamma4mgg(0);

  // count gammas over a threshold
  for (size_t ipp=0; ipp<hint_particles.size(); ipp++) {
    const CalorimeterParticle& pip = hint_particles[ipp];

    if ( pip.GetID() == CalorimeterParticle::GAMMA )
      if (pip.GetE() > options.ecut4mgg_histo)
        n_gamma4mgg++;
  }
  // loop over particles
  for (size_t ig1=0; ig1<hint_particles.size(); ig1++) {
    if ( debug ) cout << " ig1 = " << ig1 << endl;
    const CalorimeterParticle& g1 = hint_particles[ig1];

    // select gammas
    if ( !g1.GetID()==CalorimeterParticle::GAMMA )
      continue;

    double eg1 = g1.GetE();
    double xg1 = g1.GetX();
    double yg1 = g1.GetY();

    // loop over particles a second time to build combinations with
    // the first gamma
    for (size_t ig2=ig1+1; ig2<hint_particles.size(); ig2++) {
      if( debug ) cout << " ig2 = " << ig2 << endl;
      const CalorimeterParticle& g2 = hint_particles[ig2];

      double eg2 = g2.GetE();
      double xg2 = g2.GetX();
      double yg2 = g2.GetY();

      double distx = xg1 - xg2;
      double disty = yg1 - yg2;
      double dist = sqrt( distx*distx + disty*disty );

      if ( g2.GetID()==CalorimeterParticle::GAMMA ) {
        TLorentzVector pgg = g1.GetMomentum() + g2.GetMomentum();
        double mgg = pgg.M();
        if( eg1 > options.ecut4mgg_histo && eg2 > options.ecut4mgg_histo && n_gamma4mgg == 2 ) {
          if(levelBPHistos > 1) {
            h.distGamGamE->Fill(dist,eg1+eg2);
            h.distXYGamGam->Fill(distx,disty);
          }

          h.MGamGam->Fill(mgg);
        }

        if( g1.GetMainCells().size() <= 0 ||  g2.GetMainCells().size() <= 0 ) {
          if( g1.GetMainCells().size() <= 0 ) cerr << " Calorimeter " << GetName()
                                                    << " main cells not defined in hint gamma " << ig1 << endl;
          if( g2.GetMainCells().size() <= 0 ) cerr << " Calorimeter " << GetName()
                                                    << " main cells not defined in hint gamma " << ig2 << endl;
          continue;
        }

        if( debug ) cout << " g1 mains " << g1.GetMainCells().size() << " g2 mains " << g2.GetMainCells().size() << endl;
        int main_gam1 = g1.GetMainCells()[0];
        int main_gam2 = g2.GetMainCells()[0];
        if( debug ) cout << " main_gam1 = " << main_gam1 <<  " main_gam2 = " << main_gam2 << endl;
        if( BadCellsAround(main_gam1).size() > 0 || BadCellsAround(main_gam2).size() > 0 ) {
          h.MGamGam_NearBadCell->Fill(mgg);
        }
      }

      if(levelBPHistos > 1) {
        if ( g2.GetID()==CalorimeterParticle::ELECTRON )
          h.distXYElectronGam->Fill(-distx,-disty);

        if ( g2.GetID()==CalorimeterParticle::POSITRON )
          h.distXYElectronGam->Fill(-distx,-disty);

        if ( g2.GetID()==CalorimeterParticle::PION_PLUS ) {
          h.distGamPiE->Fill(dist,g2.GetE());
          h.distXYPiGam->Fill(-distx,-disty);
        }

        if ( g2.GetID()==CalorimeterParticle::PION_MINUS ) {
          h.distGamPiE->Fill(dist,g2.GetE());
          h.distXYPiGam->Fill(-distx,-disty);
        }
      }
    }
  }

  for(size_t icc=0; icc<reco_particles.size(); icc++)                   // loop over all clusters
  {
    const CalorimeterParticle &pe = reco_particles[icc];
    double eg = pe.GetE();
    if( rpart_tr_ass[icc].size() > 0 )
      h.effCaloAss->Fill(eg,1.);
    else
      h.effCaloAss->Fill(eg,0.);
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillFitInfoHist ( void )
{
  bool debug = false;
//   bool debug = true;

  if( debug ) cout << " Calorimeter::FillFitInfoHist options.fill_fit_info_histo = " <<
                               options.fill_fit_info_histo <<" h= " << fit_info_histo_ << endl;
  if( debug )	cout << " Calorimeter::FillFitInfoHist exit for debug " << endl;
  if( !options.fill_fit_info_histo ) return;

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
  TDirectory *dir_checked = NULL;

  if( fit_info_histo_ == NULL )
  {
    fit_info_histo_ = new FitInfoHisto();
    FitInfoHisto &h = *fit_info_histo_;
    bool ok = false;
    if( h.root_dir==NULL )
    {
      char dir_name[111];
      if(debug) cout << " FillFitInfoHist gDirectory->cd(/) " << endl;
      ok = gDirectory->cd("/");
      if(debug) cout << " FillFitInfoHist gDirectory->cd(/) OK " << ok << endl;
      assert(ok);
      sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
      dir_checked = myROOT_utils::TDirectory_checked(dir_name);
      ok = dir_checked->cd();
      assert(ok);
      if(debug) cout << " FillFitInfoHist gDirectory->cd("<< string(dir_name) << " ) OK " << ok << endl;
      sprintf(dir_name,"FitInfo");
      h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
    }
    if(debug) cout << " FillFitInfoHist h.root_dir->cd() " << h.root_dir << endl;
    ok = h.root_dir->cd();
    if(debug) cout << " FillFitInfoHist h.root_dir->cd() OK " << ok <<  endl;
    assert(ok);
    if( NCells() <= 0 )
    {
      cerr << " Calorimeter::FillFitInfoHist ERROR: Calorimeter " << GetName() << " is empty !!! " << endl;
      return;
    }

    double e_peak = options.calibration_energy;
    if(!h.booked_ok_)
    {
      char hist_name[132];
      char dir_name[111];

      sprintf(hist_name,"%s:: CALIB Fitted Peak Position in GeV ",GetName().c_str());
      h.MeanCALIB = myROOT_utils::TH1D_checked("FitCalibMean",hist_name,5000,0,2.*e_peak);
      sprintf(hist_name,"%s:: CALIB Fitted Sigma in GeV ",GetName().c_str());
      h.SigmaCALIB = myROOT_utils::TH1D_checked("FitCalibSigma",hist_name,5000,0,e_peak/2.);
      sprintf(hist_name,"%s:: CALIB mean Norm ",GetName().c_str());
      h.MeanNCALIB = myROOT_utils::TH1D_checked("FitCalibMeanN",hist_name,5000,-2.5,2.5);
      sprintf(hist_name,"%s:: CALIB Sigma Norm ",GetName().c_str());
      h.SigmaNCALIB = myROOT_utils::TH1D_checked("FitCalibSigmaN",hist_name,5000,-2.5,2.5);

      sprintf(hist_name,"%s:: CALIB Fitted Stat ",GetName().c_str());
      h.StatCALIB = myROOT_utils::TH1D_checked("FitCalibStat",hist_name,5000,0,100000.);
      if( debug ) cout << " MeanCALIB " << h.MeanCALIB <<" SigmaCALIB " << h.SigmaCALIB <<" StatCALIB " << h.StatCALIB  << endl;

    }
    if(debug) cout << " Start FILLING histo " << endl;

    for (size_t icell = 0; icell < NCells(); icell++)
    {
      if( debug ) cout << " cell " << icell << " " << GetCellName(icell) << endl;
      const FitInfo &fiti = GetFitInfo( CALIB, icell);
      if( fiti.fit_ok )
      {
        if( debug ) cout << " Fit OK " << fiti.fit_parameters.size() << endl;
        if( fiti.fit_parameters.size() >=3 )
	{
	  double stat = fiti.fit_parameters[0].first;
	  double mean = fiti.fit_parameters[1].first;
	  double sigma = fiti.fit_parameters[2].first;
          if( debug ) cout << " mean " << mean <<" sigma " << sigma <<" stat " << stat << endl;
          h.MeanCALIB->Fill( mean );
          h.SigmaCALIB->Fill( sigma );
          h.StatCALIB->Fill( stat );
          h.MeanNCALIB->Fill( mean/e_peak-1. );
          h.SigmaNCALIB->Fill( sigma/mean );
        }
      }
      else
      {
        if( debug ) cout << " Fit Not OK " << endl;
      }
    }

  }
  dir_save->cd();

//   if( debug )
//   {
//     cout << " Calorimeter::FillFitInfoHist exit for debug " << endl;
//     exit(0);
//   }

}

////////////////////////////////////////////////////////////////////////////////
} // namespace Reco
