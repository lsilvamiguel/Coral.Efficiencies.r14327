#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include "DaqDataDecoding/ChipADC.h"
#include "DaqDataDecoding/ChipSADC.h"
#include "DaqDataDecoding/ChipF1.h"
#include "CsEvent.h"

#include "CsECAL0.h"

#include "Reco/CalorimeterHist.h"
#include "Reco/CellDataRaw.h"
#include "Reco/Exception.h"
#include "Reco/Shower.h"

#include "CsTimingInfo.h"

using namespace std;
using namespace Reco;

////////////////////////////////////////////////////////////////////////////////

CsECAL0::CsECAL0(const string &name, const string &geom_file)
 : CsCalorimeter(name, geom_file)
{
}

////////////////////////////////////////////////////////////////////////////////

void  CsECAL0::Initialize( void )
{
  CsCalorimeter::Initialize();
  for( size_t it=0; it!=NCells(); it++ )
  {
    sadc_readout_[it] = false;
  }

  InitXY();
  if( GetNColumns() != 33 )
  {
    cout << " FATAL bug in calorimeter " << GetName() << " geometry: Ncols=" << GetNColumns() << endl;
    exit(1);
  }
  if( GetNRows() != 27 )
  {
    cout << " FATAL bug in calorimeter " << GetName() << " geometry: Nrows= " << GetNRows() << endl;
    exit(1);
  }

  if( options.print_general_info ) PrintGeneralInfo();
}

////////////////////////////////////////////////////////////////////////////////

void  CsECAL0::InitOptions( void )
{
  CsCalorimeter::InitOptions();

  options.ecut4mgg_histo                = 5.;
  options.default_calibration           = 0.020;
  options.readout_sparsified_           = true;
  options.delta_spars_mode_             = 3.;

  // MC options
  options.mc_default_calibration     = 0.020;
  options.mc_make_real_digitization  = true;
  options.mc_make_fiadc_digitization = true;
  options.mc_fiadc_sparce_delta      = 3.;

  options.SetMiscDoubleOption(1, 3.6 );
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL0::SetDefaultSADCDigitizationParameters( void )
{
  sadc_format_version_=1;
  sadc_decode_version_=1;
  sadc_delta_ = 10;
  sadc_ped_min_ = 0;
  sadc_ped_max_ = sadc_ped_min_+4;
  sadc_signal_min_ = 5;
  sadc_signal_max_ = sadc_signal_min_+20;
  sadc_front_ = 7;
  sadc_front_gate_ = 3;
  sadc_max_position_ = sadc_front_ + sadc_front_gate_;
  coeff_convert2max_ = 0.12;
  sadc_clock_ = 12.86;
  make_tcs_corrections_ = true;
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL0::FillHistoLED4SADC(void)
{
//  bool debug = true;
 bool debug = false;
 bool ok=false;
  if( debug ) cout << " CsECAL0::FillHistoLED4SADC options.fill_histos = " << options.fill_histos << endl;
  if( !options.fill_histos )
    return;

  static bool first_time = true;
  static TH1D* h1_LEDPED4SADC[5000];
  static TH1D* h1_LEDCell4SADC[5000];
  static TProfile* p1_LEDsCell4SADC[5000];
  static TProfile* p1_LEDDiff4SADC[5000];
  static TH1D* h1_LEDsTime4SADC[5000];
  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
  if( first_time )
  {
    if( debug ) cout << " CsECAL0::FillHistoLED4SADC BOOKING " << endl;
    first_time = false;
    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    TDirectory *dir = myROOT_utils::TDirectory_checked(dir_name);
    ok=(dir->cd());
    assert(ok);
    sprintf(dir_name,"%s_LED_Cells4SADC",GetName().c_str());
    dir = myROOT_utils::TDirectory_checked(dir_name);
    ok=(dir->cd());
    assert(ok);
    char hist_name[132];
    for( size_t i=0; i< NCells(); i++ )
    {
//      if( !sadc_readout_[i] ) continue;
      char name[132];
      sprintf(name,"LEDPED4SADC%s",GetCellName(i).c_str());
      sprintf(hist_name," PED in LED event in cell %s ",GetCellName(i).c_str());
      h1_LEDPED4SADC[i] = myROOT_utils::TH1D_checked(name,hist_name,500,0.,500.);
      sprintf(name,"LED4SADC%s",GetCellName(i).c_str());
      sprintf(hist_name," LED in cell %s ",GetCellName(i).c_str());
      h1_LEDCell4SADC[i] = myROOT_utils::TH1D_checked(name,hist_name,128,0.,128.);
      sprintf(name,"LEDsub4SADC%s",GetCellName(i).c_str());
      sprintf(hist_name," LED with pedestal subtracted in cell %s ",GetCellName(i).c_str());
      p1_LEDsCell4SADC[i] = myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.);
      sprintf(name,"LEDDiff4SADC%s",GetCellName(i).c_str());
      sprintf(hist_name," LED Differential signal in cell %s ",GetCellName(i).c_str());
      p1_LEDDiff4SADC[i] = myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.);
      sprintf(name,"LEDTime4SADC%s",GetCellName(i).c_str());
      sprintf(hist_name," LED Time signal in cell %s ",GetCellName(i).c_str());
      h1_LEDsTime4SADC[i] = myROOT_utils::TH1D_checked(name,hist_name,500,-50.,50.);
    }
    if( debug ) cout << " CsECAL0::FillHistoLED4SADC BOOKING  OK " << endl;
  }

  if( debug ) cout << " CsECAL0::FillHistoLED4SADC Start filling sadc_samples_.size() = " << sadc_samples_.size() << endl;
//   for( size_t it=0; it < sadc_samples_.size(); it++ )
  for( map< int, const std::vector<CS::uint16> * >::const_iterator its=sadc_samples_.begin(); its != sadc_samples_.end(); its++ )
  {
    if( debug ) cout << " New sample " <<endl;
//     int ic = sadc_samples_[it].first;
    int ic = its->first;
//     const vector<uint16> &sample = *sadc_samples_[it].second;
    const vector<uint16> &sample = *(its->second);
//    vector<uint16> sample = *sadc_samples_[it].second;
    if( sadc_format_version_ == 1 )
    {
//       for( unsigned i=0; i < sample.size(); i++ )
//       {
//         sample[i] = 1023-sample[i];
//       }
    }

    double ped =0.;

    for( int is=sadc_ped_min_; is <= sadc_ped_max_; is++)
    {
      ped += sample[is];
    }
    ped /= (sadc_ped_max_-sadc_ped_min_+1);
    h1_LEDPED4SADC[ic]->Fill(ped);

    for( int is=0; is < (int)(sample.size()-1); is++)
    {
      h1_LEDCell4SADC[ic]->Fill(is, sample[is]);
      p1_LEDsCell4SADC[ic]->Fill(is, double(sample[is])-ped );
      if(is > 0 ) p1_LEDDiff4SADC[ic]->Fill(is, double(sample[is] - sample[is-1]) );
//       cout << " " << sample[is]-cells_info[PED][OLD][icell].GetMean();
//       summ +=  sample[is];
    }
//     summ /=  double(sample.size()-1);
//  look for corresponding time
    for( size_t itt=0; itt < signals_.size(); itt++ )
    {
      if( signals_[itt].GetCellIdx() == (unsigned)ic )
        h1_LEDsTime4SADC[ic]->Fill( signals_[itt].GetTime() );
    }
    if( debug ) cout << " sample OK " << endl;
  }
  if( debug ) cout << " CsECAL0::FillHistoLED4SADC  OK " << endl;

  dir_save->cd();
  if( debug ) cout << " CsECAL0::FillHistoLED4SADC relly OK " << endl;
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL0::CalibrateAtElectronBeam(void)
{
//   bool debug = true;
  bool debug = false;
//  bool only_one_particle = false;
  bool only_one_particle = true;
  if( debug ) cout << " CsECAL0::CalibrateAtElectronBeam Try to calibrate at electron beam " << endl;
// Disable nparticle=1  if(reco_particles.size()!=1) return;
  if( only_one_particle && reco_particles.size()!=1 ) return;
  int ngam = 0;
  for( int i=0; i < (int)reco_particles.size(); i++ )
    if( reco_particles[i].GetE() > options.particle_energy_threshold ) ngam++;

  if( debug ) cout << " ngam= " << ngam << endl;

  int nb_min_part = 1;
  int nb_max_part = 10;
  if( ngam < nb_min_part) return;
  if( ngam > nb_max_part ) return;
  double e_electron = 40.;
  double emin_electron_cut = 5.;
  double emax_electron_cut = 200.;
  if(options.calibration_energy >0) e_electron = options.calibration_energy;
  if(options.calibration_energy_min >= 0) emin_electron_cut = options.calibration_energy_min;
  if(options.calibration_energy_max >  0) emax_electron_cut = options.calibration_energy_max;


  int indx = 0;
  if( reco_particles.size() > 1 )
    for (size_t i=0; i<reco_particles.size(); i++)
      if( reco_particles[i].GetE() > reco_particles[indx].GetE() ) indx = i;

  double e = reco_particles[indx].GetE(); // This is particle with maximum energy
  if( !(e > options.calibration_energy_min && e < options.calibration_energy_max) ) return;

  double factor =  e/e_electron;
  size_t icell = reco_particles[indx].GetClusterData()[0].first;
  cells_info[CALIB][NEW][icell].Add( factor , 1 );

  if( !options.fill_histos ) return;
  if( debug ) cout << " Try to Book/Fill calibration histos at electron beam " << endl;

  if( gDirectory==NULL )
  {
    return;
  }

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.

// create test histos directory sructure
  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  bool ok = false;
  if( h.root_dir==NULL )
  {
    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    h.root_dir=myROOT_utils::TDirectory_checked(dir_name);
  }
  ok=(h.root_dir->cd());
  assert(ok);
  if( h.cells_calib_dir==NULL )
  {
    char dir_name[111];
    sprintf(dir_name,"CellsCalibration_%s",GetName().c_str());
    h.cells_calib_dir=myROOT_utils::TDirectory_checked(dir_name);
    ok=(h.cells_calib_dir->cd());
    assert(ok);

    char hist_name[132];

    //    double emin = 0, emax = options.hist_energy_max;
    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"Calib%s",GetCellName(i).c_str());
      sprintf(hist_name," Electron peak in cell %s ",GetCellName(i).c_str());
//      h.h1_CalibCell[i] = myROOT_utils::TH1D_checked(name,hist_name,100,emin,emax);
      h.h1_CalibCell[i] = myROOT_utils::TH1D_checked(name,hist_name,options.calib_histo_size,0.,options.hist_energy_max);

    }
    options.calib_histo_units_in_gev=true;
  }
  if(!h.cells_calib_dir->cd()) assert(false);

  if(h.h1_CalibCell[icell]!=NULL) h.h1_CalibCell[icell]->Fill(e);
  else cout << " h.h1_CalibCell["<<icell<<"]=NULL"<<endl;

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL0::CalibrateAtElectronBeamSpecial(void)
{
//  bool debug = true;
  bool debug = false;
  if( debug ) cout << " CsECAL0::CalibrateAtElectronBeamSpecial Try to calibrate at electron beam " << endl;
// Disable nparticle=1  if(reco_particles.size()!=1) return;
  int ngam = 0;
  for (size_t i=0; i<reco_particles.size(); i++)
    if( reco_particles[i].GetE() > options.particle_energy_threshold ) ngam++;

  if( debug ) cout << " ngam= " << ngam << endl;

  int nb_min_part = 2;
  int nb_max_part = 2;
  if( ngam < nb_min_part) return;
  if( ngam > nb_max_part ) return;
  double e_electron = 40.;
  double emin_electron_cut = 5.;
  double emax_electron_cut = 200.;
  if(options.calibration_energy >0) e_electron = options.calibration_energy;
  if(options.calibration_energy_min >= 0) emin_electron_cut = options.calibration_energy_min;
  if(options.calibration_energy_max >  0) emax_electron_cut = options.calibration_energy_max;


  size_t indx0 = 0;
  if( reco_particles.size() > 1 )
    for(size_t i=0; i<reco_particles.size(); i++)
      if( reco_particles[i].GetE() > reco_particles[indx0].GetE() ) indx0 = i;

  size_t indx1 = indx0;
  if( reco_particles.size() > 1 )
  {
    for(size_t i=0; i<reco_particles.size(); i++)
    {
      if( i != indx0 ) // skip particle with max energy
      {

        if( indx1 == indx0 )
          indx1 = i;
        else
          if(reco_particles[i].GetE() > reco_particles[indx1].GetE()) indx1 = i;
      }
    }
  }

  if(indx1 == indx0) return;

  double e = reco_particles[indx0].GetE() +
                reco_particles[indx1].GetE(); // This Energy summ of two particles with maximum energy
  if( !(e > options.calibration_energy_min && e < options.calibration_energy_max) ) return;

  double factor =  e/e_electron;
  size_t icell0 = reco_particles[indx0].GetClusterData()[0].first;
  size_t icell1 = reco_particles[indx1].GetClusterData()[0].first;
  double disty = reco_particles[indx0].GetY() - reco_particles[indx1].GetY();
// disty selection
  if( fabs(disty) > 50. ) return;

  cells_info[CALIB][NEW][icell0].Add( factor , 1 );
  cells_info[CALIB][NEW][icell1].Add( factor , 1 );

  if( !options.fill_histos ) return;
  if( debug )
    cout << " CsECAL0::CalibrateAtElectronBeamSpecial Try to Book/Fill calibration histos at electron beam " << endl;

  if( gDirectory==NULL )
  {
    return;
  }

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.

// create test histos directory sructure
  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  bool ok=false;
  if( h.root_dir==NULL )
  {
    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    h.root_dir=myROOT_utils::TDirectory_checked(dir_name);
  }
  ok=(h.root_dir->cd());
  assert(ok);
  if( h.cells_calib_dir==NULL )
  {
    char dir_name[111];
    sprintf(dir_name,"CellsCalibration_%s",GetName().c_str());
    h.cells_calib_dir=myROOT_utils::TDirectory_checked(dir_name);
    ok=(h.cells_calib_dir->cd());
    assert(ok);

    char hist_name[132];

    double emin = 0, emax = options.hist_energy_max;
    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"Calib%s",GetCellName(i).c_str());
      sprintf(hist_name," Electron peak of 2 gammas in cell %s ",GetCellName(i).c_str());
      h.h1_CalibCell[i] = myROOT_utils::TH1D_checked(name,hist_name,100,emin,emax);

    }
    options.calib_histo_units_in_gev=true;
  }
  if(!h.cells_calib_dir->cd()) assert(false);

  if(h.h1_CalibCell[icell0]!=NULL) h.h1_CalibCell[icell0]->Fill(e);
  else cout << " h.h1_CalibCell["<<icell0<<"]=NULL"<<endl;

  if(h.h1_CalibCell[icell1]!=NULL) h.h1_CalibCell[icell1]->Fill(e);
  else cout << " h.h1_CalibCell["<<icell1<<"]=NULL"<<endl;

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void   CsECAL0::DrawCellHisto          (int icell)
{
//   cout << " Start with visualisation of " << GetName() << endl;
  bool ok = false;
  if( icell <0 || icell >= (int)NCells() )
  {
     cerr << " wrong cell number in DrawCellHisto " << endl;
     return;
  }

  if( !options.fill_histos )
    return;

  CalorimeterHist &h= *calo_hist;  // create a short name for calo_hist
//   CalorimeterCalibHist &hc=  *calib_hist;  // create a short name for calib_hist

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

//   sprintf(dir_name,"CellsCalibration_%s",GetName().c_str());
//   if( !gDirectory->cd(dir_name) )
//   {
//     cerr.form("Can not change directory to %s\n",dir_name);
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

  if (h.c == NULL)
  {
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
    if(h.h1_CalibCell[icell] != NULL) h.h1_CalibCell[icell]->Draw();
    h.c->cd(4);
    if(h.h1_ADCmCell[icell] != NULL) h.h1_ADCmCell[icell]->Draw();
  }
  else
  {
    h.c->Clear();
    h.c->Divide(2,2);
    h.c->cd(1);
    if(h.h1_ErCell[icell] != NULL) h.h1_ErCell[icell]->Draw();
    h.c->cd(2);
    if(h.h1_LEDCell[icell] != NULL) h.h1_LEDCell[icell]->Draw();
    h.c->cd(3);
    if(h.h1_CalibCell[icell] != NULL) h.h1_CalibCell[icell]->Draw();
    h.c->cd(4);
    if(h.h1_ADCmCell[icell] != NULL) h.h1_ADCmCell[icell]->Draw();
    h.c->Update ();
    h.c->Show ();
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

bool   CsECAL0::GoodZone4Electrons          (double x, double y, double e)  const
{
  if( fabs(y) > 300. ) return false;
  if( fabs(x) < 250. || fabs(x) > 1100. ) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
#include "CsCalorimeterHist.h"

void CsECAL0::EndOfJob( void )
{
  Reco::Calorimeter::EndOfJob();
  cout << " CsECAL0 now we need to write SADC profiles, wait a bit .. " << endl;
  if( test_histo_sadc_led_ != NULL )
  {
    TestHistoSADC_LED &h = *test_histo_sadc_led_;  // create a short name
    cout << "SADC profile LED for " << GetName() << endl;
    for( size_t i=0; i< subsets_list2store_sadc_histo_.size()+1; i++ )
    {
      if( i > 0 )
      {
        int isub = subsets_list2store_sadc_histo_[i-1];
        cout << "SubSet_" << GetSubSets()[isub].GetName() << endl;
      }

      cout << " ShapeNM10SADC " << endl;
      TProfile *p = h.p1_ShapeNM10SADC[i];
      TAxis *axis = p->GetXaxis();
      Int_t ncx   = axis->GetNbins();
      Int_t bin;
      for (bin=1; bin < ncx+1; bin++)
      {
        Double_t v = p->GetBinContent(bin);
        char c[100];
        if( v == 0. )
          sprintf(c," %2.0f,",v);
        else
          sprintf(c," %5.3f,",v);
        string ss(c);
        cout << ss;
      }
      cout << endl;

      cout << " ShapeNM100SADC " << endl;
      p = h.p1_ShapeNM100SADC[i];
      axis = p->GetXaxis();
      ncx   = axis->GetNbins();
      for (bin=1; bin < ncx+1; bin++)
      {
        Double_t v = p->GetBinContent(bin);
        char c[100];
        if( v == 0. )
          sprintf(c," %2.0f,",v);
        else
          sprintf(c," %5.3f,",v);
        string ss(c);
        cout << ss;
      }
      cout << endl;

      cout << " ShapeNM500SADC " << endl;
      p = h.p1_ShapeNM500SADC[i];
      axis = p->GetXaxis();
      ncx   = axis->GetNbins();
      for (bin=1; bin < ncx+1; bin++)
      {
        Double_t v = p->GetBinContent(bin);
        char c[100];
        if( v == 0. )
          sprintf(c," %2.0f,",v);
        else
          sprintf(c," %5.3f,",v);
        string ss(c);
        cout << ss;
      }
      cout << endl;

      cout << " ShapeNM1000SADC " << endl;
      p = h.p1_ShapeNM1000SADC[i];
      axis = p->GetXaxis();
      ncx   = axis->GetNbins();
      for (bin=1; bin < ncx+1; bin++)
      {
        Double_t v = p->GetBinContent(bin);
        char c[100];
        if( v == 0. )
          sprintf(c," %2.0f,",v);
        else
          sprintf(c," %5.3f,",v);
        string ss(c);
        cout << ss;
      }
      cout << endl;

    }
  }

  if( test_histo_sadc_more_ != NULL )
  {
    TestHistoSADCMore &h = *test_histo_sadc_more_;  // create a short name
    cout << "SADC profile RealSignals for " << GetName() << endl;
    for( size_t i=0; i< subsets_list2store_sadc_histo_.size()+1; i++ )
    {
      if( i > 0 )
      {
        int isub = subsets_list2store_sadc_histo_[i-1];
        cout << "SubSet_" << GetSubSets()[isub].GetName() << endl;
      }

      cout << " ShapeNM10SADC " << endl;
      TProfile *p = h.p1_ShapeNM10SADC[i];
      TAxis *axis = p->GetXaxis();
      Int_t ncx   = axis->GetNbins();
      Int_t bin;
      for (bin=1; bin < ncx+1; bin++)
      {
        Double_t v = p->GetBinContent(bin);
        char c[100];
        if( v == 0. )
          sprintf(c," %2.0f,",v);
        else
          sprintf(c," %5.3f,",v);
        string ss(c);
        cout << ss;
      }
      cout << endl;

      cout << " ShapeNM100SADC " << endl;
      p = h.p1_ShapeNM100SADC[i];
      axis = p->GetXaxis();
      ncx   = axis->GetNbins();
      for (bin=1; bin < ncx+1; bin++)
      {
        Double_t v = p->GetBinContent(bin);
        char c[100];
        if( v == 0. )
          sprintf(c," %2.0f,",v);
        else
          sprintf(c," %5.3f,",v);
        string ss(c);
        cout << ss;
      }
      cout << endl;

      cout << " ShapeNM500SADC " << endl;
      p = h.p1_ShapeNM500SADC[i];
      axis = p->GetXaxis();
      ncx   = axis->GetNbins();
      for (bin=1; bin < ncx+1; bin++)
      {
        Double_t v = p->GetBinContent(bin);
        char c[100];
        if( v == 0. )
          sprintf(c," %2.0f,",v);
        else
          sprintf(c," %5.3f,",v);
        string ss(c);
        cout << ss;
      }
      cout << endl;

      cout << " ShapeNM1000SADC " << endl;
      p = h.p1_ShapeNM1000SADC[i];
      axis = p->GetXaxis();
      ncx   = axis->GetNbins();
      for (bin=1; bin < ncx+1; bin++)
      {
        Double_t v = p->GetBinContent(bin);
        char c[100];
        if( v == 0. )
          sprintf(c," %2.0f,",v);
        else
          sprintf(c," %5.3f,",v);
        string ss(c);
        cout << ss;
      }
      cout << endl;

    }
  }
}

////////////////////////////////////////////////////////////////////////////////

int  CsECAL0::StoreSADCDigit              (const int& icell, const CS::ChipSADC::Digit& digit, const bool& isLed) {
  if( CsCalorimeter::StoreSADCDigit(icell, digit, isLed) ) {

    CS::ChipSADC::DataID data_id = digit.GetDataID();
    int16 srcid = data_id.GetSourceID();
    int16 port = data_id.GetPort();
    int16 chip =  data_id.GetChip();
    int16 channel =  digit.GetChannel();

    if( channel > 15 )
      printf("\tsrcid = %i , port = %i , chip = %i , channel = %i\n", srcid, port, chip ,channel);

    float time_cor = 0.;
    if( !skip_cs_time_corrections_ )
    {
      try
      {
        time_cor = CsTimingInfo::Instance()->GetEC02Shift(srcid, port, chip, channel);
      }
      catch( std::exception &e )
      {
        cerr << " CsECAL0::StoreSADCDigit CsTimingInfo->GetEC02Shift exit4debug  Exception:\n" << e.what() << endl;
        exit(1);
      }
      catch( ... )
      {
        cerr << " From CsECAL0::StoreSADCDigit CsTimingInfo->GetEC02Shift exit4debug Unknown exception: " << endl;
        exit(1);
      }
    }
    signals_.back().SetTime(signals_.back().GetTime()+time_cor);
    return 1;

  }

  return 0;


}
