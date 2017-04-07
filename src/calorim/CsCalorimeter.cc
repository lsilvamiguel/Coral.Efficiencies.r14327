/*!
   \file    CsCalorimeter.cc
   \brief   Base class for COMPASS calorimeters
   \version $Revision: 1.139 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \author  Denis Murashev
   \date    $Date: 2011/02/16 15:54:16 $
*/

#include <fstream>
#include <iostream>
#include <limits>    // for NaN
#include <sstream>
#include <string>

#include "CsInit.h"
#include "CsEvent.h"
#include "CsCalorimeter.h"
#include "CsCalorimeterGUI.h"
#include "CsOpt.h"
#include "CsHistograms.h"
#include "DaqDataDecoding/ChipADC.h"
#include "DaqDataDecoding/ChipF1.h"
#include "CDB.h"
#include "CsDigitizerSADC.h"
#include "DigitizerSADCN.h"
#include "DigitizerSADCF.h"
#include "CsTimingInfo.h"
#include "coral_config.h"
#include "CsGeant3.h"


#include "Reco/CellDataRaw.h"
#include "Reco/Exception.h"
#include "Reco/ShowerProfileLednev.h"
#define CDB_LINE_MAX 132

using namespace std;
using namespace Reco;
using CS::DetID;

CsCalorimeter::CsCalorimeter(const string &det_name, const string &geom_file) :
 Reco::Calorimeter(det_name),
 CsDet(DetID(det_name),det_name),
 skip_decoding_(false),
 skip_led_decoding_(false),
 skip_cs_time_corrections_(false),
 new_sadc_decoding_(false),
 make_sadc_cluster_filtering_(false),
 make_sadc_histo_(false),
 correct_leds_by_fem_signal_(false),
 mantab_sadc_(NULL),
 tcs_phase_(numeric_limits<double>::quiet_NaN()),
 time_in_spill_(numeric_limits<double>::quiet_NaN()),
 raw_amp_cut_delta_(0.),
 raw_time_cut_min_(-10000.),
 raw_time_cut_max_(10000.),
 sadc_shape_filter_apply_(false),
 sadc_shape_filter_slc_use_(false),
 sadc_shape_filter_line_fit_use_(true),
 test_histo_sadc_(NULL),
 test_histo_sadc_more_(NULL),
 test_histo_sadc_led_(NULL),
 dig_sadc_(NULL),
 dig_led_sadc_(NULL),
 hists_book_level(0),
 clean_cs_digits_(false),
 add_sadcinfo_to_cs_digits_(false),
 add_all_sadcinfo_to_cs_digits_(false),
 add_timeinfo_to_cs_digits_(false),
 skip_make_cs_digits_(false),
 hist_MC_hitE(NULL),
 hist_MC_hitN(NULL),
 hist_MC_hitT(NULL),
 hist_MC_hitNcellMCTrackID(NULL),
 cs_calo_hist(NULL),
 response_mc(MCSmeared)
{
  bool debug = false;
  if( debug ) cout << " CsCalorimeter::CsCalorimeter debug " << GetName() << endl;

  options.update_ecut_using_prob_cell_noise     = true;
  subsets_list2store_sadc_histo_.clear();

  ReadGeom(det_name, geom_file);

  InitMemory();

  if( CsInit::Instance()->IsAMonteCarloJob()) options.is_real_data = false;

  trigger_groups.clear();

  CsHistograms::SetCurrentPath("/Calorimeter");
  SetHistogramsBaseDir(gDirectory);

  if( debug ) cout << " CsCalorimeter::CsCalorimeter debug OK " << GetName() << endl;
}

///////////////////////////////////////////////////////////////////////////////

void  CsCalorimeter::Initialize( void )
{
  bool debug = false;
  if( debug ) cout << " CsCalorimeter::Initialize debug " << GetName() << endl;

  if( debug ) cout << " InitOptions " << GetName() << endl;
  InitOptions();
  if( debug ) cout << " ReadOptions " << GetName() << endl;
  ReadOptions();
  if( debug ) cout << " SetOptionsOK " << GetName() << endl;
  options.SetOK();

  if( debug ) cout << " UpdateAfterOptionsSettings " << GetName() << endl;
  UpdateAfterOptionsSettings();

  if( debug ) cout << " Setup readout type info " << GetName() << endl;
  for(size_t it=0; it!=NCells(); it++)
  {
    read_out_fe_amp_.push_back(UNKNOWN_FEA);
  }

  if( debug ) cout << " Reco::Calorimeter::Init " << GetName() << endl;
  Reco::Calorimeter::Init();
  if( debug ) cout << " Reco::Calorimeter::InitXY " << GetName() << endl;
  Reco::Calorimeter::InitXY();

  if( debug ) cout << " InitDefaultReadOut " << GetName() << endl;
  InitDefaultReadOut();

  if( debug ) cout << " PrintGeneralInfo " << GetName() << endl;
  if( GetOptions().print_general_info ) PrintGeneralInfo();

  if( debug ) cout << " CsCalorimeter::Initialize debug OK " << GetName() << endl;

}

////////////////////////////////////////////////////////////////////////////////

void  CsCalorimeter::InitOptions( void )
{
  Reco::Calorimeter::InitOptions();
  options.SetMiscDoubleOption(1, 2.8 );
}

////////////////////////////////////////////////////////////////////////////////

void  CsCalorimeter::ReadOptions( void )
{
  string o;
  if( CsOpt::Instance()->getOpt( GetName(), "MC response", o ) )
  {
    if( o=="MCExact" ) {
      response_mc = MCExact;
      options.mc_smear_response   = false;
    } else
    if( o=="MCSmeared" ) {
      response_mc = MCSmeared;
      options.mc_smear_response   = true;
      options.mc_smear_constant   = true;
      options.mc_smear_stochastic = true;
      options.mc_smear_readout    = true;
    } else
    {
      response_mc = MCSmeared;
      options.mc_smear_response   = true;
      options.mc_smear_constant   = true;
      options.mc_smear_stochastic = true;
      options.mc_smear_readout    = true;
      cerr << "Bad MC response option for " << getName() << ", I will use MCSmeared " << endl;
    }
  }
  if( CsOpt::Instance()->getOpt( GetName(), "correct_leds_by_fem_signal" ) )
  {
    correct_leds_by_fem_signal_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " correct_leds_by_fem_signal option was set !  " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "clean_cs_digits" ) )
  {
    clean_cs_digits_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " clean_cs_digits option was set!!!! This defenitely make a crash for standard CORAL !!! " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "add_sadcinfo_to_cs_digits" ) )
  {
    add_sadcinfo_to_cs_digits_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " add_sadcinfo_to_cs_digits option was set ! " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "add_all_sadcinfo_to_cs_digits" ) )
  {
    add_all_sadcinfo_to_cs_digits_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " add_all_sadcinfo_to_cs_digits option was set ! " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "add_timeinfo_to_cs_digits" ) )
  {
    add_timeinfo_to_cs_digits_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " add_timeinfo_to_cs_digits option was set ! " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "skip_make_cs_digits" ) )
  {
    skip_make_cs_digits_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " skip_make_cs_digits_cs_digits option was set ! " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "make_sadc_histo" ) )
  {
    make_sadc_histo_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " make_sadc_histo option was set!!!! " << endl;
  }
 if( CsOpt::Instance()->getOpt( GetName(), "make_sadc_cluster_filtering" ) )
  {
    make_sadc_cluster_filtering_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " make_sadc_cluster_filtering_ option was set!!!! " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "new_sadc_decoding" ) )
  {
    new_sadc_decoding_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " make_sadc_cluster_filtering_ option was set!!!! " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "skip_decoding" ) )
  {
    skip_decoding_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " skip_decoding option was set!!!! " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "skip_led_decoding" ) )
  {
    skip_led_decoding_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " skip_led_decoding option was set!!!! " << endl;
  }
  if( CsOpt::Instance()->getOpt( GetName(), "skip_cs_time_corrections" ) )
  {
    skip_cs_time_corrections_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " skip_cs_time_corrections option was set!!!! " << endl;
  }

  // read the shower profile for the fortran reconstruction
  if (CsOpt::Instance()->getOpt(GetName(), "FortranShowerProfA", options.fortran_showerprof_a)) {
    if (options.fortran_showerprof_a.size()!=4) {
      size_t size=options.fortran_showerprof_a.size();
      options.fortran_showerprof_a.clear();
      std::cout << "CsCalorimter::ReadOptions: Size of shower profile A array is "
                << size << ", not 4!" << std::endl;
      exit(1);
      // TODO: this should be a fatal exception, but it is caught somewhere
      throw Reco::Exception("CsCalorimter::ReadOptions: Size of shower profile A array is %zu, not 4!",
                            size);
    }
  }
  if (CsOpt::Instance()->getOpt(GetName(), "FortranShowerProfB", options.fortran_showerprof_b)) {
    if (options.fortran_showerprof_b.size()!=4) {
      size_t size=options.fortran_showerprof_b.size();
      options.fortran_showerprof_b.clear();
      std::cout << "CsCalorimter::ReadOptions: Size of shower profile B array is "
                << size << ", not 4!" << std::endl;
      exit(1);
      // TODO: this should be a fatal exception, but it is caught somewhere
      throw Reco::Exception("CsCalorimter::ReadOptions: Size of shower profile B array is %zu, not 4!",
                            size);
    }
  }
  if (options.fortran_showerprof_a.size()!=options.fortran_showerprof_b.size()) {
    size_t sizeA=options.fortran_showerprof_a.size();
    size_t sizeB=options.fortran_showerprof_b.size();
    options.fortran_showerprof_a.clear();
    options.fortran_showerprof_b.clear();
    std::cout << "CsCalorimter::ReadOptions: Size of shower profile arrays does not agree (a: "
                << sizeA << ", b: " << sizeB << ")" << std::endl;
    exit(1);
    // TODO: this should be a fatal exception, but it is caught somewhere
    throw Reco::Exception("CsCalorimter::ReadOptions: Size of shower profile arrays does not agree (a: %zu, b: %zu)",
                          sizeA, sizeB);
  }

  std::list<int> msadc_srcids;
  if( CsOpt::Instance()->getOpt( GetName(), "MSADC SrcIDs:", msadc_srcids ) )
  {
    cout << "WARNING!! In CsCalorimeter " << getName() <<
              " The following SrcIDs were by force assigned to MSADC: ";
    for( std::list<int>::const_iterator it = msadc_srcids.begin(); it!=msadc_srcids.end(); it++ )
    {
      msadc_srcids_.push_back(*it);
      cout <<" "<< msadc_srcids_.back();
    }
    cout << endl;
  }

  // Try to read BookHistograms option from CORAL configuration file.
  if( !CsOpt::Instance()->getOpt( GetName(), "BookHistograms", hists_book_level ) )
    hists_book_level=0; // We failed do this. Let's use the default value.

  // Book histograms only if hists_book_level>0
  if( hists_book_level>0 )
  {
    // Change histograms directory to the "/DetectorName"
    CsHistograms::SetCurrentPath(string("/")+GetName());
    // Create a histogram.
    hist_MC_hitE = new CsHist1F("MC_hitE","Energy from MC hits",100,0,100);
    hist_MC_hitN = new CsHist1F("MC_hitNcell","Energy from MC hits in Cells",NCells(),0,NCells());
    hist_MC_hitT = new CsHist1F("MC_hitT","Time from MC hits",100,0,100);
    if( hists_book_level>1 ) hist_MC_hitNcellMCTrackID = new CsHist2F("MC_hitNcellMCTrackID"," MC track ID from MC hits in Cells",NCells(),0,NCells(),50,0,50.);
    // Go back to the top directory
    CsHistograms::SetCurrentPath("/");
    options.SetRecoOptions("hist=yes");
  }
  else
    options.SetRecoOptions("hist=no");

  //  Well ... bad, I know ...
  if( GetName() == "EC01P1__" ) options.tolerance_for_nearby_cells = 15.0;
  if( GetName() == "EC02P1__" ) options.tolerance_for_nearby_cells = 2.5;

  // Obtain RecoOptions from CORAL.  First apply last "RecoOption" line, then
  // apply all "MoreRecoOptions" lines.
  list<string> reco_opt_list;
  CsOpt::Instance()->getOpt( GetName(), "RecoOptions",     reco_opt_list, false );
  CsOpt::Instance()->getOpt( GetName(), "MoreRecoOptions", reco_opt_list, true );
  for( list<string>::const_iterator it=reco_opt_list.begin(); it!=reco_opt_list.end(); it++ )
    options.SetRecoOptions(*it);

  // read some calibrations the are mainly specific to the ReconstructionCombined
  std::vector<double> paramD;
  if (CsOpt::Instance()->getOpt(GetName(), "ParamEnergyError", paramD)) {
    if (paramD.size()!=3)
      CsErrLog::msg(elFatal, __FILE__, __LINE__,
                    "CsCalorimeter::ReadOptions: parameterization of energy error for calorimeter %s is wrong, expected 3 parameters, found %zu.",
                    getName(), paramD.size());
    options.param_energy_error = paramD;
  }

  if (CsOpt::Instance()->getOpt(GetName(), "ParamTimeError", paramD)) {
    if (paramD.size()!=3)
      CsErrLog::msg(elFatal, __FILE__, __LINE__,
                    "CsCalorimeter::ReadOptions: parameterization of time error for calorimeter %s is wrong, expected 3 parameters, found %zu.",
                    getName(), paramD.size());
    options.param_time_error = paramD;
  }

  std::vector<int> paramI;
  if (CsOpt::Instance()->getOpt(GetName(), "RecoCombined NrAllowedShowers", paramI))
    options.combined_allowed_showers = paramI;
}


////////////////////////////////////////////////////////////////////////////////

void  CsCalorimeter::InitDefaultReadOut( void )
{
  string value;

  // Set default to MaxAdvanced
  options.digitizer_method = CsDigitizerSADC::MaxAdvanced;
  options.digitizer_method_led = CsDigitizerSADC::MaxAdvanced;

  if( CsOpt::Instance()->getOpt( GetName(), "DIGITIZER_METHOD", value ) )
    {
      if( value == "MaxAdvanced" )
        options.digitizer_method = CsDigitizerSADC::MaxAdvanced;
      else if( value == "MaxSimple" )
        options.digitizer_method = CsDigitizerSADC::MaxSimple;
      else if( value == "SummRange" )
        options.digitizer_method = CsDigitizerSADC::SummRange;
      else if( value == "ShapeFit" )
        options.digitizer_method = CsDigitizerSADC::ShapeFit;
      else if( value == "FitPulse" )
        options.digitizer_method = CsDigitizerSADC::FitPulse;
      else
        throw Reco::Exception("CsCalorimeter::InitDefaultReadOut: Bad option DIGITIZER_METHOD = \"%s\" \n", value.c_str());
    }



  // Set default to be the same as for physiks pulses
  options.digitizer_method_led = options.digitizer_method;
   if( CsOpt::Instance()->getOpt( GetName(), "DIGITIZER_METHOD_LED", value ) )
    {
      if( value == "MaxAdvanced" )
        options.digitizer_method_led = CsDigitizerSADC::MaxAdvanced;
      else if( value == "MaxSimple" )
        options.digitizer_method_led = CsDigitizerSADC::MaxSimple;
      else if( value == "SummRange" )
        options.digitizer_method_led = CsDigitizerSADC::SummRange;
      else if( value == "ShapeFit" )
        options.digitizer_method_led = CsDigitizerSADC::ShapeFit;
      else if( value == "FitPulse" )
        options.digitizer_method = CsDigitizerSADC::FitPulse;
      else
        throw Reco::Exception("CsCalorimeter::InitDefaultReadOut: Bad option DIGITIZER_METHOD_LED = \"%s\" \n", value.c_str());
    }

}

////////////////////////////////////////////////////////////////////////////////

void  CsCalorimeter::InitReadOut( void )
{
  bool debug = false;
  CsCalorimeter::GetSADCDigitizationOptions(); // Move this call to implementation classes
  InitArraySADC();

  if( debug ) cout << " CsCalorimeter:: create SADC digitizers " <<  endl;
  if( GetDigitizersSADC().size() != NCells() )
  {
    cerr <<" CsCalorimeter::InitReadOut " << GetName() << " internal error " << GetDigitizersSADC().size() << endl;
    exit(1);
  }

  ManagerShapeTableSADC *mtab = NULL;
  if( sadc_decode_version_ == 2 ) //Create Table Manager
  {
    mtab = new ManagerShapeTableSADC(this);
    mantab_sadc_= mtab;
  }


  for ( unsigned ic=0; ic<NCells(); ic++ )
  {
    if(read_out_fe_amp_[ic] == CsCalorimeter::SADC_FEA ||read_out_fe_amp_[ic] == CsCalorimeter::MSADC_FEA )
    {
      CsDigitizerSADC::Type type_sadc = CsDigitizerSADC::SADC;
      if( read_out_fe_amp_[ic] == CsCalorimeter::MSADC_FEA)
        type_sadc = CsDigitizerSADC::MSADC;
      switch(sadc_decode_version_) {
        case 1: {
          CsDigitizerSADC::Shape rd_table = CsDigitizerSADC::RD_ECAL2;
          CsDigitizerSADC *d = new CsDigitizerSADC( this, (CsDigitizerSADC::Method)options.digitizer_method, type_sadc, rd_table );
          d->use_shape_filter_line_fit_=sadc_shape_filter_line_fit_use_;
          AddDigitizerSADC(d,ic);
          break;
        }
        case 2: { // We'll use hardcoded values for sadc_decode_version_ for a while it is not yet stable code
          DigitizerSADCN *d = new DigitizerSADCN( this, type_sadc);
          d->mtab_ = mtab;
          AddDigitizerSADC(d,ic);
          break;
        }
        case 3: {
          DigitizerSADCF *d = new DigitizerSADCF( this, type_sadc, ic);
          AddDigitizerSADC(d,ic);
          break;
        }
        default:
          std::cerr <<" InitReadOutSADC " << GetName() <<" unknown sadc_decode_version = " << sadc_decode_version_ << std::endl;
          exit(1);
          break;
      }
      if( !GetDigitizersSADC()[ic] )
      {
        cerr << " Memory problems in CsCalorimeter::InitReadOut " << GetName() << " ?? Cell " << GetCellName(ic) << endl;
        exit(1);
      }

      {
        CsDigitizerSADC::Shape led_table = CsDigitizerSADC::LED_ECAL2;
        CsDigitizerSADC *d = new CsDigitizerSADC( this, (CsDigitizerSADC::Method)options.digitizer_method_led, type_sadc, led_table );
        if( !d )
        {
          cerr << " Memory problems in CsCalorimeter::InitReadOut " << GetName() << " ?? " << endl;
          exit(1);
        }
        d->option_store_stat_info_ = true;
        AddLedDigitizerSADC(d,ic);
      }


    }
    else if(read_out_fe_amp_[ic] == CsCalorimeter::FIADC_FEA )
    {
//        This is FIADC, it seems no actions needed for the moment
    }
    else if(read_out_fe_amp_[ic] == CsCalorimeter::UNKNOWN_FEA  )
    {
      cerr << " Warning in calorimeter " << GetName() <<" !! : The cell " <<
                  GetCellName(ic) << " does not attached to any readout !!! " << endl;
    }
    else
    {
      cerr << " Internal error : Should never be here " << endl;
    }
  }

  if( GetDigitizersSADC().size() != NCells() )
  {
    cerr <<" CsCalorimeter::InitReadOut " << GetName() << " internal error " << GetDigitizersSADC().size() << endl;
    exit(1);
  }
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::ReadGeom(const string &det_name, const string &geom_file) {
  bool debug = false;
  string com_geant_det; // non empty (3bytes string) if calorimeter was found.
  ifstream f(geom_file.c_str());

  if( !f.is_open() )
    throw Reco::Exception("CsCalorimeter::ReadGeom():  can not open file \"%s\"",geom_file.c_str());

  // Loop over all lines of input file.
  map < double, pair <int,int> > rad_len_stats;
  for( size_t line_n=1; true; line_n++ ) {
    string line;
    getline(f,line);

    if( f.eof() )
      break;

    istringstream s(line.c_str());

    string opt;
    s >> opt;

    if( opt == "calo" )
    {
      string TBname, name;
      double pos_xyz[3]; // in cm

      if( !(s >> TBname >> name >> pos_xyz[2] >> pos_xyz[0] >> pos_xyz[1]) )
        throw Reco::Exception("CsCalorimeter::ReadGeom(): Bad format in file %s:  \"%s\"",
                               geom_file.c_str(),line.c_str());

      if( TBname==det_name )
      {
        if( !com_geant_det.empty() )
          throw Reco::Exception("CsCalorimeter::ReadGeom(): there are two "
                                "definitions for calorimeter %s", det_name.c_str());

        com_geant_det=string(name,0,3);

        SetPosition(pos_xyz[0]*10., pos_xyz[1]*10., pos_xyz[2]*10.);
      }
    }
    else
    if( opt == "cmtx" )
    {
      // Calorimeter matrix definition

      string name;
      size_t ID, type, nmod, n_rows, n_cols, offset[2];
      double x_pos, y_pos, z_pos;
      double x_size, y_size, z_size;
      double x_step, y_step;
      double tgate, Tresh, GeV_to_ADC;
      double rad_len, abs_len, StochT, ConstT, active_material;
      size_t nholes;

      if( !(
      s >> ID
        >> name
        >> type
        >> nmod
// #warning Amount of cells in configuration file is (Y,X) and NOT (X,Y)
        >> n_rows >> n_cols
// #warning Parameter "offset" is not used.
        >> offset[0] >> offset[1]
        >> z_size >> x_size >> y_size
        >> z_pos  >> x_pos  >> y_pos
        >> x_step >> y_step
        >> tgate >> Tresh >> GeV_to_ADC
        >> rad_len >> abs_len >> StochT >> ConstT >> active_material >> nholes )
       )
        throw Reco::Exception("CsCalorimeter::ReadGeom(): Bad format in file %s:  \"%s\"",
                               geom_file.c_str(), line.c_str());

      // critical energy is only contained in new versions of detectors.dat,
      // therefore we default to a reasonable (average) value
      double crit_energy = 11.;
      s >> crit_energy;
      crit_energy /= 1000.;  // convert to GeV

      string hw_type;
      s >> hw_type;

      string com_geant_det_alt = com_geant_det;
      com_geant_det_alt[1] = com_geant_det[2];
      com_geant_det_alt[2] = com_geant_det[1];

      // discard lines belonging to different calorimeters
      if( (name.compare(0,3,com_geant_det) ) && (name.compare(0,3,com_geant_det_alt) ) )
        continue;

      if( nholes!=0 )
        throw Reco::Exception("CsCalorimeter::ReadGeom(): nholes!=0 is not "
                              "supported yet. File \"%s\"", geom_file.c_str());

      // Convert cm to mm
      x_size  *= 10.;
      y_size  *= 10.;
      z_size  *= 10.;
      x_pos   *= 10.;
      y_pos   *= 10.;
      z_pos   *= 10.;
      x_step  *= 10.;
      y_step  *= 10.;
      rad_len *= 10.;
      abs_len *= 10.;

      rad_len_stats[rad_len].first++;
      rad_len_stats[rad_len].second += nmod;

      if ( n_rows*n_cols != nmod )
        throw Reco::Exception("CsCalorimeter::ReadGeom(): "
			      "nrow*ncol!=nmod in file %s:  \"%s\"",
			      geom_file.c_str(), line.c_str());

      double ReadOutT = 0.01;  // quick fix for ReadOutTerm

      char cell_type_name[500];
      sprintf(cell_type_name, "Cell_%zu_%s", cells_type.size(), com_geant_det.c_str());
      CellType new_cell_type = Reco::CellType(cell_type_name,
					      x_size, y_size, z_size,
					      x_step, y_step,
					      rad_len, abs_len, crit_energy,
					      StochT, ConstT, ReadOutT, active_material,
					      0.,  // density is unused
					      hw_type);
      InsertCellType( new_cell_type );
      Reco::CellType *cell_type = &cells_type.back();

      ReadShowerProfile(*cell_type);

      matrixes.push_back( CellsMatrix("matrix", *cell_type, n_rows, n_cols,
                                      x_pos, y_pos, z_pos) );

#ifdef calorim_debug
      matrixes.back().Print(cout,det_name+"("+com_geant_det+"):  ");
#endif

      for( size_t y=0; y<n_rows; y++ )
        for( size_t x=0; x<n_cols; x++ )
        {
          const double pos_x = x*cell_type->GetStepX() + x_pos;
          const double pos_y = y*cell_type->GetStepY() + y_pos;
          const double pos_z = 0                       + z_pos;

          if( !matrixes.back().IsInHole(x,y) )
            cells.push_back( Reco::Cell(*cell_type,true,pos_x,pos_y,pos_z) );
        }

      comgeant_first_cell_n.push_back(ID);
    }
  }

  if( debug )
  {
    for ( map< double,pair<int,int> >::const_iterator rl=rad_len_stats.begin(); rl!=rad_len_stats.end(); rl++)
      cout << GetName() << ": " << rl->second.first << " matrices with "
           << rl->second.second << " blocks with rad. len. " << rl->first << "mm" << endl;
  }

  if( com_geant_det.empty() )
    throw Reco::Exception("CsCalorimeter::ReadGeom():  calorimeter \"%s\" was not found in file \"%s\"",
                           GetName().c_str(),geom_file.c_str());

}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::ReadShowerProfile(Reco::CellType& cellType) {
    // use Reco::ShowerProfile... classes?
    std::string tmpSP;
    if (CsOpt::Instance()->getOpt(GetName(), "ShowerProfile", tmpSP)) {
        if (tmpSP=="Lednev") {
            std::string hwtype = Reco::CellType::GetHwTypeString(cellType.GetHwType());
            std::string key    = "ShowerProfileParams " + hwtype;

            // Lednev shower profile requires 4 parameters per cell type
            std::vector<double> params;
            if (CsOpt::Instance()->getOpt(GetName(), key, params)) {
                if (params.size()==4 || params.size()==6) {
                    ShowerProfile* showerProfile = new ShowerProfileLednev(this, cellType, params);
                    cellType.SetShowerProfile(showerProfile);
                } else
                    CsErrLog::msg(elFatal, __FILE__, __LINE__,
                                  "%s: shower profile for cell type \"%s\" needs to have 4 or 6 parameters, not %zu.",
                                  GetName().c_str(), hwtype.c_str(), params.size());

            } else {
                if (options.combined_reconstruction) // strictly speaking this is only really fatal in case
                                                     // of using ReconstructionCombined
                    CsErrLog::msg(elFatal, __FILE__, __LINE__,
                                  "%s: cannot find parameters for shower profile of cell type \"%s\"",
                                  GetName().c_str(), hwtype.c_str());
                else
                    CsErrLog::msg(elInfo, __FILE__, __LINE__,
                                  "%s: cannot find parameters for shower profile of cell type \"%s\"",
                                  GetName().c_str(), hwtype.c_str());
            }
        } else {
            if (options.combined_reconstruction) // strictly speaking this is only really fatal in case
                                                 // of using ReconstructionCombined
                CsErrLog::msg(elFatal, __FILE__, __LINE__,
                              "%s: unknown shower profile \"%s\"",
                              GetName().c_str(), tmpSP.c_str());
            else
                CsErrLog::msg(elInfo, __FILE__, __LINE__,
                              "%s: unknown shower profile \"%s\"",
                              GetName().c_str(), tmpSP.c_str());
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::ProcLED( void )
{
  Calorimeter::ProcLED();
}

////////////////////////////////////////////////////////////////////////////////

bool CsCalorimeter::IsItMSADCSrcID(int src_id) const
{
  bool yes_it_is = false;
  for( size_t it=0; it<msadc_srcids_.size(); it++ )
  {
    if( src_id == msadc_srcids_[it] )
    {
      yes_it_is = true;
      break;
    }
  }
  return yes_it_is;
}

////////////////////////////////////////////////////////////////////////////////

double CsCalorimeter::GetFFTwidth(const int icell) const {
  if(FFTwidth_.size() <= (unsigned int)icell) {
    std::cerr << " GetFFTwidth " << GetName() <<" FFTwidth calibratrion has not enough entries" << std::endl;
    std::cerr << " expected at least " << icell << " found "<< FFTwidth_.size() << std::endl;
    exit(1);
  }
  return FFTwidth_.at(icell);
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::SetDefaultSADCDigitizationParameters( void )
{
// SADC digitisation parameters
  sadc_format_version_=1;
  sadc_decode_version_=1;
  sadc_delta_ = 1;
  sadc_ped_min_ = 0;
  sadc_ped_max_ = sadc_ped_min_+4;
  sadc_signal_min_ = 5;
  sadc_signal_max_ = sadc_signal_min_+20;
  sadc_front_ = 7;
  sadc_front_gate_ = 3;
  coeff_convert2max_ = 0.12;
  sadc_clock_ = 12.86;
  sadc_max_position_  = sadc_front_ + sadc_front_gate_ ;
  make_tcs_corrections_ = true;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::SetDaqDataDecodingInfoGeneral( const multimap<CS::Chip::DataID,CS::Chip::Digit*> &daqmap )
{
  for( multimap<CS::Chip::DataID,CS::Chip::Digit*>::const_iterator it=daqmap.begin();
       it!=daqmap.end(); it++ )
  {
    SetDaqDataDecodingInfo( *(it->second) );
  }
  InitReadOut();
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::SetDaqDataDecodingInfo( const CS::Chip::Digit &d  )
{
  bool debug = false;
  const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit*>( &d );
  if(ds != NULL )
  {
    if( GetName() == ds->GetDetID().GetName() )
    {
      int32 x_digit = ds->GetX();
      int32 y_digit = ds->GetY();
      int icell = GetCellOfColumnRow( x_digit, y_digit);
      if( icell >= 0 )
      {
        if( debug )
        {
          cout <<" Yees !! it is defenitely MY DIGIT !! For cell " << icell << endl;
          cout << " It was FE " << read_out_fe_amp_[icell] << " and now its is " << CsCalorimeter::SADC_FEA << endl;
        }

        CS::ChipSADC::DataID id( ds->GetDataID() );
        int src_id = id.GetSourceID();
        if( read_out_fe_amp_[icell]!=CsCalorimeter::UNKNOWN_FEA )
        {
          cerr << " Warning in calorimeter " << GetName() <<"!! : The cell " <<
                      GetCellName(icell) << " is attached to more than one readout!!!" << endl;
        }
        if( IsItMSADCSrcID(src_id) )
        {
          read_out_fe_amp_[icell]=CsCalorimeter::MSADC_FEA;
        }
        else
        {
          read_out_fe_amp_[icell]=CsCalorimeter::SADC_FEA;
        }
      }
    }
    return;
  }

  const CS::ChipADC::Digit *df = dynamic_cast<const CS::ChipADC::Digit*>( &d );
  if(df != NULL )
  {
    if( GetName() == df->GetDetID().GetName() )
    {
      int32 x_digit = df->GetX();
      int32 y_digit = df->GetY();
      int icell = GetCellOfColumnRow( x_digit, y_digit);
      if( icell >= 0 )
      {
        if( debug )
        {
          cout <<" Yees !! it is defenitely MY DIGIT !! For cell " << icell << endl;
          cout << " It was FE " << read_out_fe_amp_[icell] << " and now its is " << CsCalorimeter::FIADC_FEA << endl;
        }
        if( read_out_fe_amp_[icell]!=CsCalorimeter::UNKNOWN_FEA )
        {
          cerr << " Warning in calorimeter " << GetName() <<"!! : The cell " <<
                      GetCellName(icell) << " is attached to more than one readout!!!" << endl;
        }
        read_out_fe_amp_[icell]=CsCalorimeter::FIADC_FEA;
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::GetSADCDigitizationOptions( void)
{
  bool debug = false;

  SetDefaultSADCDigitizationParameters();
  string value;
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_DECODE_VERSION", value ) )
  {
    sadc_decode_version_ = atoi(value.c_str());
    if( debug )
      cout << GetName() << " SADC_DECODE_VERSION" << sadc_decode_version_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_MAKE_TCS_CORRECTIONS", value ) )
  {
    if( value == "NO" || value == "no" || value == "No" )
    {
      make_tcs_corrections_ = false;
    }
    else
    {
      make_tcs_corrections_ = true;
    }

    if( debug )
    {
      if( make_tcs_corrections_ )
         cout << " GetName() " << " MAKE_TCS_CORRECTIONS ON " << endl;
      else
         cout << " GetName() " << " MAKE_TCS_CORRECTIONS OFF " << endl;
    }
  }
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_DELTA", value ) )
  {
    sadc_delta_ = atoi(value.c_str());
    if( debug ) cout << " GetName() " << "SADC_DELTA" << sadc_delta_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_PED_MIN", value ) )
  {
    sadc_ped_min_ = atoi(value.c_str());
    if( debug ) cout << " GetName() " << "SADC_PED_MIN" << sadc_ped_min_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_PED_MAX", value ) )
  {
    sadc_ped_max_ = atoi(value.c_str());
    if( debug ) cout << " GetName() " << "SADC_PED_MAX" << sadc_ped_max_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_SIGNAL_MIN", value ) )
  {
    sadc_signal_min_ = atoi(value.c_str());
    if( debug ) cout << " GetName() " << "SADC_SIGNAL_MIN" << sadc_signal_min_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_SIGNAL_MAX", value ) )
  {
    sadc_signal_max_ = atoi(value.c_str());
    if( debug ) cout << " GetName() " << "SADC_SIGNAL_MAX" << sadc_signal_max_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_FRONT", value ) )
  {
    sadc_front_ = atoi(value.c_str());
    if( debug ) cout << " GetName() " << "SADC_FRONT" << sadc_front_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_FRONT_GATE", value ) )
  {
    sadc_front_gate_ = atoi(value.c_str());
    if( debug ) cout << " GetName() " << "SADC_FRONT_GATE" << sadc_front_gate_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "SADC_CONVERT2MAX", value ) )
  {
    coeff_convert2max_ = atof(value.c_str());
    if( debug ) cout << " GetName() " << "SADC_CONVERT2MAX" << coeff_convert2max_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "RAW_AMP_CUT_DELTA", value ) )
  {
    raw_amp_cut_delta_ = atof(value.c_str());
    if( debug ) cout << " GetName() " << "RAW_AMP_CUT_DELTA" << raw_amp_cut_delta_ << endl;
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( GetName(), "RAW_TIME_CUT_MIN", value ) )
  {
    raw_time_cut_min_ = atof(value.c_str());
    if( debug ) cout << " GetName() " << "RAW_TIME_CUT_MIN" << raw_time_cut_min_ << endl;
    value.clear();
  }

  if( CsOpt::Instance()->getOpt( GetName(), "SADC_SHAPE_FILTER_ON", value ) )
  {
    if( value == "YES" || value == "yes" )
    {
      sadc_shape_filter_apply_ = true;
    }
    else
    {
      sadc_shape_filter_apply_ = false;
    }
  }

  if( CsOpt::Instance()->getOpt( GetName(), "SADC_SHAPE_FILTER_SLC_USE", value ) )
  {
    if( value == "YES" || value == "yes" )
    {
      sadc_shape_filter_slc_use_ = true;
    }
    else
    {
      sadc_shape_filter_slc_use_ = false;
    }
  }

  if( CsOpt::Instance()->getOpt( GetName(), "SADC_SHAPE_FILTER_LINE_FIT_USE", value ) )
  {
    if( value == "YES" || value == "yes" )
    {
      sadc_shape_filter_line_fit_use_ = true;
    }
    else
    {
      sadc_shape_filter_line_fit_use_ = false;
    }
  }

  if( CsOpt::Instance()->getOpt( GetName(), "RAW_TIME_CUT_MAX", value ) )
  {
    raw_time_cut_max_ = atof(value.c_str());
    if( debug ) cout << " GetName() " << "RAW_TIME_CUT_MAX" << raw_time_cut_max_ << endl;
    value.clear();
  }

  double shape1000[] = {   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
                            0.,    0.01,0.046,0.174, 0.526, 0.885,    1.,0.9345, 0.806,0.6585,
                        0.5353,0.4348,0.3556,0.2951,0.2482,0.2109,0.1809,0.1575,0.1395,0.1253,
                        0.1093,0.0919, 0.073, 0.061, 0.051, 0.045, 0.040, 0.035, 0.030, 0.025,
                         0.020, 0.015, 0.010, 0.008, 0.006, 0.004, 0.002, 0.001, 0.0007,0.0005};

  for( size_t i=0; i<50; i++ )
  {
    sadc_shape1000_.push_back(shape1000[i]);
  }
}

////////////////////////////////////////////////////////////////////////////////

bool CsCalorimeter::AddMCHit(int detector_number,const void *data)
{
//  bool debug = true;
  bool debug = false;
  bool result = false;

  if( detector_number!=4 )
    return false;                     // No, this is not a calorimter hit...

  const CalorimeterMCData &d = *(const CalorimeterMCData*)(data);

  if( debug ) cout << "cell_id = " << d.cell_id << "  E=" << d.dE << "  dT=" << d.dT <<"  Track_id=" << d.track_id <<"\n";

  list<CellsMatrix>::iterator mtx = matrixes.begin();

  for( size_t i=0, cell_n=0; i<comgeant_first_cell_n.size(); i++,mtx++ )
  {
    assert(mtx!=matrixes.end());
    size_t first_cell = comgeant_first_cell_n[i];

    if( d.cell_id>first_cell && (d.cell_id-first_cell)<=mtx->Size() )
    {
      // OK, this is hit from THIS calorimeter and THIS matrix.

      cell_n += d.cell_id-first_cell-1;         // This is cell number in the range [0,NCells())
      if( cell_n>=NCells() )
      {
        cerr << "  ERROR in " << GetName() <<
                "  CsCalorimeter::AddMCHit(): d.cell_id= " << d.cell_id <<
                "  first_cell= " << first_cell <<
                "  cell_n= " << cell_n <<
                "  NCells= " << NCells() << endl;
        throw Reco::Exception("CsCalorimeter::AddMCHit(): d.cell_id=%d  first_cell=%d  cell_n=%d   NCells=%d",
                               d.cell_id,first_cell,cell_n,NCells());
      }

      #if calorim_debug>1
        cout << "HIT: Calorimeter " << getName() << " matrix " << i << " cell_n=" <<
                       cell_n << " dE=" << d.dE << endl;
      #endif
      if( debug )
      {
        cout << "HIT: Calorimeter " << getName() << " matrix " << i << " cell_n=" <<
                       cell_n << " dE=" << d.dE << endl;
      }
      if( hist_MC_hitN!=NULL )
        hist_MC_hitN->Fill(cell_n,d.dE);
      if( hist_MC_hitT!=NULL )
        hist_MC_hitT->Fill(d.dT,d.dE);
      if( hist_MC_hitNcellMCTrackID !=NULL )
        hist_MC_hitNcellMCTrackID->Fill(float(cell_n),float(d.track_id),d.dE);

      mc_input_.push_back( CellDataRaw(cell_n,
                                       d.dE/cells[cell_n].GetCellType().GetActiveMaterial(),
                                       d.dT) );

      result = true;
      break;
    }
    else
      cell_n += mtx->Size();
  }

  // Fill histogram.
  if( hist_MC_hitE!=NULL ) {
    double sum(0.);
    for (vector<CellDataRaw>::const_iterator it=mc_input_.begin(); it!=mc_input_.end(); it++)
      sum += it->GetEnergy();
    hist_MC_hitE->Fill(sum);
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::SetEventIDInfo( void )
{
//   bool debug = true;
  bool debug = false;
  int    run_num       = CsEvent::Instance()->getRunNumber();
  int    ev_in_run     = CsEvent::Instance()->getEventNumberInRun();
  int    spill_num     = CsEvent::Instance()->getBurstNumber();
  int    ev_in_spill   = CsEvent::Instance()->getEventNumberInBurst();
  double time_in_spill = CsEvent::Instance()->getTimeInSpill();

  // keep lowest 16 bits of trigger mask only, the other bits may contain
  // information from online filter or TCS phase measurement
  unsigned int trigger_mask = CsEvent::Instance()->getTriggerMask() & 0xffff;
  bool is_led_event = (trigger_mask==4096);

// TODO: time_t and ev_type need proper programming
  CsTime time = CsEvent::Instance()->getEventTime();
  time_t timet = CsEvent::Instance()->getEventTime().secFrEpoch();
  int ev_type = 1;

  UpdateEventID( ev_type, run_num, spill_num, ev_in_spill, ev_in_run, trigger_mask, timet, time_in_spill, is_led_event );
  if( debug ) cout <<" time_t evid.GetTime " << GetEventID().GetTime() << " double evid.GetTimeInSpilll = " << GetEventID().GetTimeInSpill() << endl;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::Clear(void)
{
  bool debug = false;
  Calorimeter::Clear();
  if( clean_cs_digits_ )
  {
    if( debug ) cout << " CsCalorimeter::Clear " << GetName() << " CleanCsDigits  " << endl;
    CleanCsDigits();
  }
  else
  {
    CsDet::Clear();
  }
  for(vector<TrigGroup>::iterator it=trigger_groups.begin(); it!=trigger_groups.end(); it++) it->Clear();


// Clean all DigitizersSADC
  ClearDigitizersSADC();

  sadc_samples_.clear();
  tcs_phase_     = numeric_limits<double>::quiet_NaN();
  time_in_spill_ = numeric_limits<double>::quiet_NaN();

}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::CleanCsDigits(void)
{
  for( std::list<CsDigit*>::iterator it = myDigits_.begin(); it != myDigits_.end(); it++ )
  {
    delete (*it);
  }
  CsDet::Clear();
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::DrawCellHisto(int icell)
{
  Reco::Calorimeter::DrawCellHisto( icell );
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::readMCCalibration(time_t timePoint)
{
  bool debug = false;
  if( debug ) cout << " CsCalorimeter::readMCCalibration " << GetName() << " mccdb_ " << mccdb_ << endl;

  // in case of combined reconstruction read some additional stuff from a calibration file
  // (this has to be before the check whether a MC data base was defined, because here
  // it is only read from a file)
  if (GetOptions().combined_reconstruction && GetOptions().use_combined_calibration) {
    std::string filename;
    if (CsOpt::Instance()->getOpt(GetName(), "RecoCombined McCalib", filename)) {
      if (filename=="")
        CsErrLog::msg(elFatal, __FILE__, __LINE__,
                      "CsCalorimeter::ReadOptions: calorimeter %s: reading of calibration file for ReconstructionCombined requested, but no file given.",
                      getName());
    }

    options.combined_calibration="";
    FILE* calibfd;
    char tmpbuf[8192];
    register int nbread = 0;
    calibfd = fopen(filename.c_str(), "r");
    if ( !calibfd ) {
      CsErrLog::msg(elFatal, __FILE__, __LINE__,
                    "CsCalorimeter::ReadOptions: calorimeter %s: calib file %s not found",
                    getName(), filename.c_str());
    }
    while (!feof(calibfd)) {
      nbread = fread(tmpbuf, 1, 8192, calibfd);
      if (!nbread) { break;}
      options.combined_calibration.append(tmpbuf, nbread);
    }
    fclose(calibfd);

    if (options.combined_calibration == "") {
      CsErrLog::msg(elFatal, __FILE__, __LINE__,
                    "CsCalorimeter::ReadOptions: calorimeter %s: empty string from calib file %s",
                    getName(), filename.c_str());
    }
  }

  if( mccdb_ == NULL ) return;
  CDB::Time tp(timePoint,0);
  tm *t = localtime(&tp.first);
// Calibrations reading
  try
  {
    if( debug ) cout << " CsCalorimeter::readMCCalibration try to read _CALIB " << endl;
    string s("");
    mccdb_->read(getName(), s, tp, "_CALIB");
    if (s == "") throw Exception("empty string from calib file");
    int err = 0;
    err = InputCalibInfo(MC,s);
    if( debug ) cout << " CsCalorimeter::readMCCalibration read _CALIB " << " error_code =" << err << endl;
  }
  catch( const std::exception &e )
  {
    cerr << "CsCalorimeter::readMCCalibration():_CALIB " << getName() << " error in reading for time point ";
    cerr << t << ": " << e.what() << endl;
  }
// Bad cells reading
  try
  {
    if( debug ) cout << " CsCalorimeter::readMCCalibration try to read _BadCells " << endl;
    string s("");
    mccdb_->read(getName(), s, tp, "_MC_BadCells");
    if (s == "") throw Exception("empty string from calib file");
    int err = InputBadCellsInfo(MC,s);
    if (debug) cout << " CsCalorimeter::readMCCalibration read _MC_BadCells " << " error_code =" << err << s <<endl;
  }
  catch( const std::exception &e )
  {
    cerr << "CsCalorimeter::readMCCalibration():_MC_BadCells " << getName() << " error in reading for time point ";
    cerr << t << ": " << e.what() << endl;
  }
// Bad cells reading for reconstruction
  try
  {
    if( debug ) cout << " CsCalorimeter::readMCCalibration try to read _BadCells " << endl;
    string s("");
    mccdb_->read(getName(), s, tp, "_BadCells");
    if (s == "") throw Exception("empty string from calib file");
    int err = InputBadCellsInfo(OLD,s);
    if (debug) cout << " CsCalorimeter::readMCCalibration read _BadCells " << " error_code =" << err << s <<endl;
  }
  catch( const std::exception &e )
  {
    cerr << "CsCalorimeter::readMCCalibration():_BadCells " << getName() << " error in reading for time point ";
    cerr << t << ": " << e.what() << endl;
  }
// energy dependent corrections
  if ( options.edep_corr != "" ) {
    try
    {
      if( debug ) cout << " CsCalorimeter::readCalibration try to read " << options.edep_corr << endl;
      string s("");
      mccdb_->read(getName(), s, tp, options.edep_corr);
      if (s == "") throw Exception("empty string from calib file");
      InputEdepCorr(s);
      if (debug) cout << " CsCalorimeter::readMCCalibration read " << options.edep_corr << "  " << s <<endl;
    }
    catch( const std::exception &e)
    {
      cerr << "CsCalorimeter::readMCCalibration():" << options.edep_corr << " " << getName() << " error in reading for time point ";
      cerr << t << ": " << e.what() << endl;
    }
  }
  if( debug ) cout << " CsCalorimeter::readMCCalibration fini " << GetName() << endl;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::readCalibration(time_t timePoint)
{
  bool debug = false;
  if( debug ) cout << " CsCalorimeter::readCalibration " << GetName() << endl;
  string tag;
  CDB::Time tp(timePoint,0);
  tm *t = localtime(&tp.first);

  bool db_debug = ( getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_CALO")!=0 );

  if (CsInit::Instance()->useConditionsDB())
  {
    cout << " CsCalorimeter::readCalibration for " << getName()
	 << " :: Sorry, ConditionsDB Calibrations for calos "
	 << "are not implemented for the moment. " << endl;
    UpdateInternalAfterNewSettings();
    return;
  }

  // MeV per ADC count
  tag = "_CALIB";
  try {
    if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
    string s("");
    cdb_->read(getName(), s, tp, "_CALIB");
    if (s == "") throw Exception("empty string from calib file");
    int err = InputCalibInfo(OLD,s);
    if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
			<< " error_code =" << err << s << endl;
  }
  catch( const std::exception &e ) {
    cerr << "CsCalorimeter::readCalibration(): " << tag << " " << getName() << " error in reading for time point ";
    cerr << t << ": " << e.what() << endl;
  }

  // frequency-based noise suppression
  if ( options.use_ProbEnT ) {
    tag = "_ProbEnT";
    if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
    string s("");
    cdb_->read(getName(), s, tp, tag);
    if (s == "") throw Exception("empty string from calib file");
    int err = InputProbEnT(s);
    if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
			<< " error_code =" << err << s << endl;
  }

  // Time Calibrations (SADC t0) reading
  tag = "_TimeCALIB";
  try {
    if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
    string s("");

    string filename;
    if( CsOpt::Instance()->getOpt( GetName(), "TimeCALIBFile",  filename) ) {
      cout<<"------------------------------------------------------------------\n";
      cout<<"--Warning "<<GetName()<<" Time calibration file specified manually\n";
      cout<<"------------------------------------------------------------------\n";
      ifstream infile(filename.c_str());
      string line;
      while(getline(infile,line)) {
        s +=line;
        s +='\n';
      }
    }
    else
      cdb_->read(getName(), s, tp, tag);
    if (s == "") throw Exception("empty string from calib file");
    int err = InputTimeCalibInfo(OLD, s);
    if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
			<< " error_code =" << err << s << endl;
  }
  catch( const std::exception &e ) {
    cerr << "CsCalorimeter::readCalibration(): " << tag << " " << getName() << " error in reading for time point ";
    cerr << t << ": " << e.what() << endl;
  }

  //ECAL2 Time Jump corrections
  tag = "_TimeCorr";
  try {
      if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
      string s("");
      if(strcmp(getName(),"EC02P1__")==0) {
        string filename;
        if( CsOpt::Instance()->getOpt( GetName(), "TimeCorrFile",  filename) ) {
          cout<<"------------------------------------------------------------------\n";
          cout<<"--Warning "<<GetName()<<" Time correction file specified manually\n";
          cout<<"------------------------------------------------------------------\n";
          ifstream infile(filename.c_str());
          string line;
          while(getline(infile,line)){
            s +=line;
            s +='\n';
          }
        }
        else
          cdb_->read(getName(), s, tp, tag);
        if (s == "") throw Exception("empty string from calib file");
        int err = CsTimingInfo::Instance()->ReadEC02ShiftFile(s);
        if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
          << " error_code =" << err << s << endl;
      }
  }
  catch( const std::exception &e ) {
    cerr << "CsCalorimeter::readCalibration(): " << tag << " " << getName() << " error in reading for time point ";
    cerr << t << ": " << e.what() << endl;
  }

  if(sadc_decode_version_ == 3) {
    tag = "_FFTWidth";
    try {
      if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
      string s("");

      string filename;
      if( CsOpt::Instance()->getOpt( GetName(), "FFTWidthFile",  filename) ) {
        cout<<"------------------------------------------------------------------\n";
        cout<<"--Warning "<<GetName()<<" Fourier width calibration file specified manually\n";
        cout<<"------------------------------------------------------------------\n";
        ifstream infile(filename.c_str());
        string line;
        while(getline(infile,line)) {
          s +=line;
          s +='\n';
        }
      }
      else
        cdb_->read(getName(), s, tp, tag);
      if (s == "") throw Exception("empty string from calib file");
      int err = InputWidthCalibInfo(OLD, s);
      if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
        << " error_code =" << err << s << endl;
    }
    catch( const std::exception &e ) {
      cerr << "CsCalorimeter::readCalibration(): " << tag << " " << getName() << " error in reading for time point ";
      cerr << t << ": " << e.what() << endl;
    }
  }


  if ( options.use_led_ref_calib_correction ) {
    tag = "_LED";
    if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
    string s;
    cdb_->read(getName(), s, tp, tag);
    if (s == "") {
      stringstream ss;
      ss << "CsCalorimeter::readCalibration(): calorimeter " << getName() << " tag = "
	 << tag << " empty string from calib file in reading for time point " << t;
      throw Exception(ss.str().c_str());
    }

    int err;
    if( !XYRegularGrid() )
      err = InputLEDInfo(OLD, s);
    else
      err = InputLEDInfoXY(OLD, s);

    if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
			<< " error_code =" << err << s << endl;
  }

  if ( options.use_led_ref_calib_correction ) {
    tag = "_LED_PRIM";
    if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
    string s("");
    cdb_->read(getName(), s, tp, tag);
    if (s == "") {
      stringstream ss;
      ss << "CsCalorimeter::readCalibration(): calorimeter " << getName() << " tag = "
	 << tag << " empty string from calib file in reading for time point " << t;
      throw Exception(ss.str().c_str());
    }

    int err;
    if( !XYRegularGrid() )
      err = InputLEDInfo(PRIM, s);
    else
      err = InputLEDInfoXY(PRIM, s);

    if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
			<< " error_code =" << err << s << endl;
  }

  if ( options.use_led_ref_calib_correction_in_spills && !XYRegularGrid() ) {
    tag = "_LEDinSpills";
    if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
    string s;
    cdb_->read(getName(), s, tp, tag);
    if (s == "") {
      stringstream ss;
      ss << "CsCalorimeter::readCalibration(): calorimeter " << getName() << " tag = "
	 << tag << " empty string from calib file in reading for time point " << t;
      throw Exception(ss.str().c_str());
    }
    int err = InputLEDInfoInSpills(OLD, s);
    if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
			<< " error_code =" << err << s << endl;
  }

  if ( options.use_led_ref_calib_correction_in_spills && XYRegularGrid() ) {
    tag = "_LEDinSpillsXY";
    if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
    string s;
    cdb_->read(getName(), s, tp, tag);
    if (s == "") {
      tag = "_LEDinSpills";
      if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
      string s;
      cdb_->read(getName(), s, tp, tag);
      if (s == "") {
	stringstream ss;
	ss << "CsCalorimeter::readCalibration(): calorimeter " << getName() << " tag = "
	   << tag << " empty string from calib file in reading for time point " << t;
	throw Exception(ss.str().c_str());
      }
      int err = InputLEDInfoInSpills(OLD, s);
      if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
			  << " error_code =" << err << s << endl;
    } else {
      int err = InputLEDInfoInSpillsXY(OLD, s);
      if( db_debug ) cout << " CsCalorimeter::readCalibration read " << tag
			  << " error_code =" << err << s << endl;
    }
  }

  // energy dependent corrections
  if ( options.edep_corr != "" ) {
    tag = options.edep_corr;
    if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
    string s;
    string filename;
    if( CsOpt::Instance()->getOpt( GetName(), "edepFile",  filename) ) {
      cout<<"------------------------------------------------------------------\n";
      cout<<"--Warning "<<GetName()<<" Edep calibration file specified manually\n";
      cout<<"--to"<< filename <<"\n";
      cout<<"------------------------------------------------------------------\n";
      ifstream infile(filename.c_str());
      string line;
      while(getline(infile,line)) {
        s +=line;
        s +='\n';
      }
    }
    else
      cdb_->read(getName(), s, tp, tag);
    if (s == "") {
      stringstream ss;
      ss << "CsCalorimeter::readCalibration(): calorimeter " << getName() << " tag = "
    	 << tag << " empty string from calib file in reading for time point " << t;
      throw Exception(ss.str().c_str());
    }
    InputEdepCorr(s);
  }

  // time in spill dependent corrections
  if ( options.tisdep_corr != "" ) {
    tag = options.tisdep_corr;
    if( debug ) cout << " CsCalorimeter::readCalibration try to read " << tag << endl;
    string s;
    string filename;
    if( CsOpt::Instance()->getOpt( GetName(), "tdepFile",  filename) ) {
      cout<<"------------------------------------------------------------------\n";
      cout<<"--Warning "<<GetName()<<" Tdep calibration file specified manually\n";
      cout<<"--to"<< filename <<"\n";
      cout<<"------------------------------------------------------------------\n";
      ifstream infile(filename.c_str());
      string line;
      while(getline(infile,line)) {
        s +=line;
        s +='\n';
      }
    }
    else
      cdb_->read(getName(), s, tp, tag);
    if (s == "") {
      stringstream ss;
      ss << "CsCalorimeter::readCalibration(): calorimeter " << getName() << " tag = "
    	 << tag << " empty string from calib file in reading for time point " << t;
      throw Exception(ss.str().c_str());
    }
    InputTiSdepCorr(s);
  }

  // in case of combined reconstruction read some additional stuff from the calibration database
  // ReconstructionCombined
  if (GetOptions().combined_reconstruction && GetOptions().use_combined_calibration) {
    tag = "_RecoCombined";
    if (debug)
      std::cout << "CsCalorimeter::readCalibration: try to read " << tag << std::endl;

    cdb_->read(getName(), options.combined_calibration, tp, tag);
    // fatal error as we request to read this
    if (options.combined_calibration == "") {
      char date[50];
      if (strftime(date, sizeof(date), "%c", t)==0) {
        strcpy(date, "<unknown date>");
      }
      CsErrLog::msg(elFatal, __FILE__, __LINE__,
                    "CsCalorimeter::readCalibration: calorimeter %s, tag %s: empty string from calib file in reading for time point %s",
                    getName(), tag.c_str(), date);
    }
  }

  const std::vector< std::string >&calib_tags = GetInputInfoNames();
    for( unsigned itc=0; itc < calib_tags.size(); itc++) {
      if( calib_tags[itc] == "_TrigGroupTimeCALIB" ) {           // Trigger Groups Time Calibrations reading
        try {
          if( debug ) cout << " CsCalorimeter::readCalibration try to read " << calib_tags[itc] << endl;
          string s("");
          cdb_->read(getName(), s, tp, calib_tags[itc]);
          if (s == "") throw Exception("empty string from calib file");
          int err = InputTrGrTimeCalibInfo(s);
          if( db_debug ) cout << " CsCalorimeter::readCalibration read " << calib_tags[itc]
			      << " error_code =" << err << s << endl;
        }
        catch( const std::exception &e ) {
          cerr << "CsCalorimeter::readCalibration(): " << calib_tags[itc] << " " << getName() << " error in reading for time point ";
          cerr << t << ": " << e.what() << endl;
        }
      }
      else if( calib_tags[itc] ==  "_TrigGroupCALIB" ) { // Trigger Groups energy calibrations (ECAL1, HCAL1/2 have summation of several cells into trigger groups)
        try {
          if( debug ) cout << " CsCalorimeter::readCalibration try to read " << calib_tags[itc] << endl;
          string s("");
          cdb_->read(getName(), s, tp, calib_tags[itc]);
          if (s == "") throw Exception("empty string from calib file");
          int err = InputTrGrCalibInfo(s);
          if( db_debug ) cout << " CsCalorimeter::readCalibration read " << calib_tags[itc]
      			<< " error_code =" << err << s << endl;
        }
        catch( const std::exception &e ) {
          cerr << "CsCalorimeter::readCalibration(): " << calib_tags[itc] << " " << getName() << " error in reading for time point ";
          cerr << t << ": " << e.what() << endl;
        }
      }
      else if( calib_tags[itc] ==  "_BadCells" ) {  // bad/dead cells (amplitude should be put to zero for these cells)
        try {
          if( debug ) cout << " CsCalorimeter::readCalibration try to read " << calib_tags[itc] << endl;
          string s("");
          cdb_->read(getName(), s, tp, calib_tags[itc]);
          if (s == "") throw Exception("empty string from calib file");
          int err = InputBadCellsInfo(OLD,s);
          if( db_debug ) cout << " CsCalorimeter::readCalibration read " << calib_tags[itc]
	      		<< " error_code =" << err << s << endl;
        }
        catch( const std::exception &e ) {
          cerr << "CsCalorimeter::readCalibration(): " << calib_tags[itc] << " " << getName() << " error in reading for time point ";
          cerr << t << ": " << e.what() << endl;
        }
      }
  // SADC pedestals:  After some discussion with Stefan Huber it seems likely
  // that the reading of pedestals is optional and that there is some code
  // which determines pedestals event-by-event.  But Sebastian Uhl adds that
  // in his opinion, these baselines are required to correct for overflows
  // (cells in saturation), which however wasn't necessary for 2009 anymore.
      else if( calib_tags[itc] ==  "_SADCInfo" ) {
        try {
          if( debug ) cout << " CsCalorimeter::readCalibration try to read " << calib_tags[itc] << endl;
          string s("");
          cdb_->read(getName(), s, tp, calib_tags[itc]);
          if (s == "") throw Exception("empty string from calib file");
          int err = InputSADCInfo(s);
          if( db_debug ) cout << " CsCalorimeter::readCalibration read " << calib_tags[itc]
	      		<< " error_code =" << err << s << endl;
        }
        catch( const std::exception &e ) {
          cerr << "CsCalorimeter::readCalibration(): " << calib_tags[itc] << " " << getName() << " error in reading for time point ";
          cerr << t << ": " << e.what() << endl;
        }
      }
      else if( calib_tags[itc] ==  "_SADCInfo" ) {
        try {
          if( debug ) cout << " CsCalorimeter::readCalibration try to read " << calib_tags[itc] << endl;
          string s("");
          cdb_->read(getName(), s, tp, calib_tags[itc]);
          if (s == "") throw Exception("empty string from calib file");
          int err = InputSADCInfo(s);
          if( db_debug ) cout << " CsCalorimeter::readCalibration read " << calib_tags[itc]
	      		<< " error_code =" << err << s << endl;
        }
        catch( const std::exception &e ) {
          cerr << "CsCalorimeter::readCalibration(): " << calib_tags[itc] << " " << getName() << " error in reading for time point ";
          cerr << t << ": " << e.what() << endl;
        }
      }
      else if( calib_tags[itc] ==  "_ShapeTableSADC" ) {
        try {
          debug = true;
          if( debug ) cout << " CsCalorimeter::readCalibration try to read " << calib_tags[itc] << endl;
          string s("");
          cdb_->read(getName(), s, tp, calib_tags[itc]);
          if (s == "") throw Exception("empty string from calib file");
          int err = InputShapeTableSADC(s);
          if( db_debug ) cout << " CsCalorimeter::readCalibration read " << calib_tags[itc]
	      		<< " error_code =" << err << s << endl;
        }
        catch( const std::exception &e ) {
          cerr << "CsCalorimeter::readCalibration(): " << calib_tags[itc] << " " << getName() << " error in reading for time point ";
          cerr << t << ": " << e.what() << endl;
        }
      }
      else if( calib_tags[itc] ==  "_TimeFWHMSADC" ) {
        try {
          if( debug ) cout << " CsCalorimeter::readCalibration try to read " << calib_tags[itc] << endl;
          string s("");
          cdb_->read(getName(), s, tp, calib_tags[itc]);
          if (s == "") throw Exception("empty string from calib file");
          int err = InputTimeFWHMSADC(s);
          if( db_debug ) cout << " CsCalorimeter::readCalibration read " << calib_tags[itc]
	      		<< " error_code =" << err << s << endl;
        }
        catch( const std::exception &e ) {
          cerr << "CsCalorimeter::readCalibration(): " << calib_tags[itc] << " " << getName() << " error in reading for time point ";
          cerr << t << ": " << e.what() << endl;
        }
      }
      else if( calib_tags[itc] ==  "_CALIB" || calib_tags[itc] ==  "_TimeCALIB" || calib_tags[itc] ==  "_TimeCorr" || calib_tags[itc] == "_ProbEnT"
                  || calib_tags[itc] ==  "_LED_PRIM"  || calib_tags[itc] ==  "_LED"   || calib_tags[itc] ==  "_LEDinSpills"  || calib_tags[itc] ==  "_LEDinSpillsXY"  ) {
          cerr <<" Warning! You dont need(not considered) to put in your Reco options: " << GetName() <<" MoreRecoOptions  SET_INPUT_INFO="<< calib_tags[itc] << endl;
      }

    }   // end of cycle over calib_tags


  UpdateInternalAfterNewSettings();
}

////////////////////////////////////////////////////////////////////////////////
// This function should be Called from Reco::Calorimeter::UpdateInternalAfterNewSettings as a virtual function
void  CsCalorimeter::UpdateFrontEndDependentSettings( void )
{
  bool debug = false;
  if( debug ) cout <<" CsCalorimeter::UpdateFrontEndDependentSettings debug " << GetName() << endl;
  Reco::Calorimeter::UpdateFrontEndDependentSettings();
  if( sadc_decode_version_ == 2 ) // Initialize DigitizerSADCN
  {
    for( size_t ic=0; ic<NCells(); ic++ )
    {
      bool digok = reinterpret_cast<DigitizerSADCN *>(GetDigitizerSADC2Modify(ic))-> Check();  // No need in dynamic_cast I hope
      if( !digok)
      {
        cerr << GetName() << " CsCalorimeter::UpdateFrontEndDependentSettings : Some fatal problems for DigitizerSADCN initialization were detected cell ic  " << ic <<" " <<GetCellName(ic) << endl;
        exit(1);
      }
    }
  }
  if( debug ) cout <<" CsCalorimeter::UpdateFrontEndDependentSettings debug OK " << GetName() << endl;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::WriteCalib(void)
{
  tm tt[2];  // Start,Finish
          // Start time
          tt[0].tm_year = 2000 - 1900;     // year 2000
          tt[0].tm_mon  = 4    - 1;        // month number 4 - April
          tt[0].tm_mday = 25;              // day 25
          tt[0].tm_hour = 0;
          tt[0].tm_min  = 0;
          tt[0].tm_sec  = 0;

          // Finish time
          tt[1].tm_year = 2000 - 1900;
          tt[1].tm_mon  = 4    - 1;
          tt[1].tm_mday = 25;
          tt[1].tm_hour = 23;
          tt[1].tm_min  = 59;
          tt[1].tm_sec  = 59;
  try
  {
    string s("");
    cout << " Try WriteToDataBase " << GetName() << endl;
    int err = OutputCalibInfo(s);
    CsDet::WriteToDataBase(string(getName())+"_CALIB", s,tt[0],tt[1]);
    s.clear();
    err = OutputTimeCalibInfo(s);
    CsDet::WriteToDataBase(string(getName())+"_TimeCALIB", s,tt[0],tt[1]);
    s.clear();
    err = OutputTrGrCalibInfo(s);
    CsDet::WriteToDataBase(string(getName())+"_TrigGroupCALIB", s,tt[0],tt[1]);
    s.clear();
    err = OutputTrGrTimeCalibInfo(s);
    CsDet::WriteToDataBase(string(getName())+"_TrigGroupTimeCALIB", s,tt[0],tt[1]);
    cout << " OK WriteToDataBase " << GetName() << endl;
//     return;
  }
  catch( const std::exception &e )
  {
    cerr << e.what() << "\n";
  }
  catch( ... )
  {
    cout << " Return FAIL " << endl;
    cerr << "CsCalorimeter::WriteCalib(): " << getName() << " error in writing for time point ";
    cerr<<tt[0].tm_mday<<"."<<tt[0].tm_mon+1<<"."<<tt[0].tm_year+1900<<" "
      <<tt[0].tm_hour<<":"<<tt[0].tm_min<<":"<<tt[0].tm_sec<<endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputCalibInfo(size_t when, const std::string &s)
{
//  bool debug = true;
  bool debug = false;
  if( debug ) cout << " CsCalorimeter::InputCalibInfo " << GetName() << " when = " << when << endl;
  return Reco::Calorimeter::InputCalibInfo( when,  s);
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputAnyCalibInfo(const std::string &tag, std::string &s,const std::string &comment) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " CsCalorimeter::OutputAnyCalibInfo " << GetName() << " tag = " << tag << endl;

  if( tag == "_SADCInfo" )
    return OutputSADCInfo(s);
  else if( tag == "_TrigGroupCALIB" )
    return OutputTrGrCalibInfo(s);
  else if( tag == "_TrigGroupTimeCALIB" )
    return OutputTrGrTimeCalibInfo(s);
  else if( tag == "_LED" )
    return OutputLEDInfo(s,comment);
  else if( tag == "_TimeCALIB4LED" ) {
    if( XYRegularGrid() )
      return OutputTimeCalibInfoXY4LED(s);
    else
      return OutputTimeCalibInfo4LED(s);
  }
  else
    return Calorimeter::OutputAnyCalibInfo(tag, s);
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputAnyCalibInfo(const std::string &tag,const std::string &s)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " CsCalorimeter::InputAnyCalibInfo " << GetName() << endl;

  if( tag == "_SADCInfo" )
  {
    return InputSADCInfo(s);
  }

  return Calorimeter::InputAnyCalibInfo(tag, s);
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputSADCInfo( const string &s)
{
  bool debug = false;
  if( debug )
  {
    cout << " Calorimeter::InputSADCInfo " << GetName() << " debug " << endl;
    cout << s;
  }
  istringstream is(s.c_str());
  char calorim_name[132],dummy[132],cellname[132];
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  sscanf(str.c_str(),"%s %s ",calorim_name,dummy);
  getline(is,str);

  getline(is,str);

  int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    // there is one type of calibration with 9 rows of data, and one with 7
    // call the one with 7 the old one
    bool oldformat(false);

    int icell,idsadc, iledev;
    float ped, dped, convert, unknown_entry_a, unknown_entry_b;

    if ( sscanf(str.c_str()," %d %d %d %g %g %g %g %g %s \n", &icell, &idsadc, &iledev,  &ped,  &dped,  &convert, &unknown_entry_a, &unknown_entry_b, cellname)!=9 )
    {
      oldformat = true;

      int ret = sscanf(str.c_str()," %d %d %d %g %g %g %s \n", &icell, &idsadc, &iledev,  &ped,  &dped,  &convert, cellname );
      assert(ret==7);
    }

    int jcell = FindCellByName( cellname );
    if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      cerr << " Unexpected Cell name " << cellname << " in CsCalorimeter::InputSADCInfo " <<
                                                                        GetName() << "  " << endl;
      cerr << " Input strig: " << endl;
      cerr << str << endl;
      cerr << " It well might be that you use wrong calibration file for this detector " << GetName() << endl;
      stringstream ss;
      ss << " Unexpected Cell name " << GetName() << " CsCalorimeter::InputSADCInfo " << str;
      throw Exception(ss.str().c_str());
    }

    if( icell != jcell )
    {
      icell = jcell;
      not_matched_cells_cnt ++;
    }
    if( GetLedDigitizersSADC()[icell] != NULL )
    {
      double ped_odd = ped - dped;
      double ped_even = ped + dped;
      double dpedd = dped;
      double cfcv = convert;

      GetDigitizersSADC()[icell]->SetCalibInfo( ped_odd, ped_even, dpedd, cfcv);
      GetLedDigitizersSADC()[icell]->SetCalibInfo( ped_odd, ped_even, dpedd, cfcv);
    }
  }
  if( not_matched_cells_cnt > 0 )
  {
    cerr << " WARNING!!! CsCalorimeter::InputSADCInfo " << GetName() <<
               " Not matching in cells id was detected " << not_matched_cells_cnt << " times " << endl;
    cerr << " You use wrong calibaration file or calibrations were produced with different geometry descriptor !!! " << endl;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputSADCInfo( string &s ) const
{
  bool debug = false;
  char o[500000];
  o[0]=0;
  int nsadc = 0;
  for( size_t it=0; it < NCells(); it++ )
  {
    if( GetLedDigitizerSADC(it) != NULL ) nsadc++;
  }
  sprintf(o,"### \n");
  sprintf(o+strlen(o)," %s SADC Info for %d out of Ncells=%zu \n", GetName().c_str(), nsadc, NCells());
  sprintf(o+strlen(o),"### \n");
  sprintf(o+strlen(o)," CELL#   MSADC? NOTLED PED  DPED    CFMAX2SUMM       SIGMA    SIGMA_CM    NAME  \n" );
  for( size_t it=0; it < NCells(); it++ )
  {
    if( GetLedDigitizerSADC(it) == NULL ) continue;
    const CsDigitizerSADC * digl = dynamic_cast <const CsDigitizerSADC *>(GetLedDigitizerSADC(it));
    if( digl == NULL )
    {
      cerr <<" Fatal error " << GetName() << " CsCalorimeter::OutputSADCInfo not yet implemented for DigitizerSADCN( " << endl;
      exit(1);
    }

    int notledev = 0;
    int is_MSADC = (int)digl->GetType();
    double stat_led = digl->stat_ped_odd_.GetEntries();
    double ped_odd = digl->stat_ped_odd_.GetMean();
    double ped_even = digl->stat_ped_even_.GetMean();
    double sig_odd = digl->stat_ped_odd_.GetSigma();
    double sig_even = digl->stat_ped_even_.GetSigma();

    float cfmax2summ = 1.;

    float sig = (sig_odd + sig_even)/2.;
    float ped = (ped_odd + ped_even)/2.;
    float dped = (ped_even - ped_odd)/2.;
    float sig_ped = digl->stat_ped_.GetSigma();

    sprintf(o+strlen(o)," %6zu %4d  %4d %6.1f %6.1f %12.6f  %12.4f %12.4f %s \n",it, is_MSADC, notledev, ped, dped,
                                                    cfmax2summ, sig, sig_ped, GetCellName(it).c_str());

  }
  string ss(o);
  s = ss;
  if( debug ) cout << s;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputTimeFWHMSADC( const string &s)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " CsCalorimeter::InputTimeFWHMSADC debug " << endl;

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret >= 2 );

//  Read 2 lines of comments
  getline(is,str);
  getline(is,str);
//  Read 1 line of comment and time offset
  getline(is,str);
  getline(is,str);
//  Read calibration
  uint not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell;
    float t,st,w;
    float mean_hfwn;
    float sigma_hfwn;
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %d %g %g %g %g %g %s", &icell, &w, &t, &st, &mean_hfwn, &sigma_hfwn, cellname);
    assert( ret == 7 );
    int jcell = FindCellByName( cellname );
    if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      cerr << " Unexpected Cell name " << cellname << " in Calorimeter::InputTimeFWHMSADC " <<
                                                        GetName() <<  endl;
      cerr << " Input string: " << endl;
      cerr << str << endl;
      cerr << " It well might be that you are use wrong calibration file for this detector " << GetName() << endl;
      stringstream ss;
      cerr << " Unexpected Cell name " << GetName() << "CsCalorimeter::InputTimeFWHMSADC " << str << endl;
      exit(1);;
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if(icell >=0 && icell < (int)NCells() )
    {
      if( debug )
      {
        cout << "  cell " << icell << " w " << w << " t " << t <<" st " << st << endl;
      }
      if( w == 0. ) w=1.;
      DigitizerSADCBase * dgsadcany = GetDigitizerSADC2Modify ( icell );
      DigitizerSADCN * dgsadc = dynamic_cast<DigitizerSADCN *>(dgsadcany);
      if( dgsadc != NULL )
      {
        dgsadc->calib_fwhm_tab_=t;
        dgsadc->calib_fwhm_sigma_tab_=st;
        dgsadc->calib_hfwn_tab_=mean_hfwn;
        dgsadc->calib_hfwn_sigma_tab_=sigma_hfwn;
      }
      else
      {
        cerr <<" Missuse of calibrations dedicated to DigitizerSADCN which is not initialized! " << endl;
      }
    }
  }

  if( debug )
  {
    cout << " CsCalorimeter::TimeFWHMSADC stop for debug " << GetName() << endl;
    exit(0);
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputShapeTableSADC( const string &s)
{
  if( mantab_sadc_ != NULL )
    return mantab_sadc_->InputShapeTableSADC(s);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputCalibInfo( std::string& s ) const
{
  if( !XYRegularGrid() )
    return Calorimeter::OutputCalibInfo(s);
  else
    return Calorimeter::OutputCalibInfoXY(s);
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputLEDInfo( size_t when, const std::string &s)
{
  if( !XYRegularGrid() )
    return Calorimeter::InputLEDInfo( when, s);
  else
    return Calorimeter::InputLEDInfoXY( when, s);
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputLEDInfo( std::string& s, const std::string &comment) const
{
// Overwrite the comment
  string comment_new("### No FEM corrections ");
  if( !comment.empty() ) comment_new = comment;
  if( !XYRegularGrid() )
    return Calorimeter::OutputLEDInfo(s, comment_new);
  else
    return Calorimeter::OutputLEDInfoXY(s, comment_new);
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputEcutCellsInfo(const std::string &s)
{
  if( !XYRegularGrid() )
    return Calorimeter::InputEcutCellsInfo(s);
  else
    return Calorimeter::InputEcutCellsInfoXY(s);
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputEcutCellsInfo( std::string& s ) const
{
  if( !XYRegularGrid() )
    return Calorimeter::OutputEcutCellsInfo(s);
  else
    return Calorimeter::OutputEcutCellsInfoXY(s);
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputTimeCalibInfo( size_t when, const std::string &s)
{
  if( !XYRegularGrid() )
    return Calorimeter::InputTimeCalibInfo( when, s);
  else
    return Calorimeter::InputTimeCalibInfoXY( when, s);
}

////////////////////////////////////////////////////////////////////////////////
int CsCalorimeter::InputWidthCalibInfo(size_t when, const string &s) {
  istringstream is(s);
  char calorim_name[256],dummy[256],cellname[256];

  string str;
  getline(is,str);
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
  getline(is,str);

  int icell;
  float t,st,w;
  while( getline(is,str) )
  {
    ret = sscanf(str.c_str(), " %d %g %g %g %s", &icell, &w, &t, &st, cellname);
    assert( ret == 5 );
    FFTwidth_.push_back(w);
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputTimeCalibInfo( std::string& s ) const
{
  if( !XYRegularGrid() )
    return Calorimeter::OutputTimeCalibInfo(s);
  else
    return Calorimeter::OutputTimeCalibInfoXY(s);
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputTimeCalibInfo4LED( string &s ) const
{
// Check if all cells attached to SADC readout
  for( size_t it=0; it < NCells(); it++ )
  {
    if( GetLedDigitizerSADC(it) == NULL ) return 1;
  }

  char o[500000];
  o[0]=0;
  sprintf(o," %s Calorimeter Time Calibration File for LEDs \n",GetName().c_str());
  sprintf(o+strlen(o)," Cell#    Statistic   Time     SigmaTime  CellName    \n");
  sprintf(o+strlen(o),"                      (ns)       (ns)                 \n");

  float t_offset=0.;
  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    if( GetLedDigitizerSADC(icell) == NULL ) continue;
    const CsDigitizerSADC * digl = dynamic_cast <const CsDigitizerSADC *>(GetLedDigitizerSADC(icell));
    if( digl == NULL )
    {
      cerr <<" Fatal error " << GetName() << " CsCalorimeter::OutputTimeCalibInfo4LED not yet implemented for DigitizerSADCN " << endl;
      exit(1);
    }
    const Reco::StatInfo &ledtime = digl->stat_time_;
    double entries = ledtime.GetEntries();
    if( entries > 10 )
    {
      double t=ledtime.GetMean();
      t_offset += t;
    }
  }

  t_offset = t_offset/float(NCells());
  sprintf(o+strlen(o)," Time  Offset  \n");
  sprintf(o+strlen(o)," %11.2f \n",t_offset );

  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    if( GetLedDigitizerSADC(icell) == NULL ) continue;
    const CsDigitizerSADC * digl = dynamic_cast <const CsDigitizerSADC *>(GetLedDigitizerSADC(icell));
    if( digl == NULL )
    {
      cerr <<" Fatal error " << GetName() << " CsCalorimeter::OutputTimeCalibInfo4LED not yet implemented for DigitizerSADCN " << endl;
      exit(1);
    }
    const Reco::StatInfo &ledtime = digl->stat_time_;
    double entries = ledtime.GetEntries();
    if( entries > 10 )
    {
      double t=ledtime.GetMean();
      double st=ledtime.GetSigma();
      sprintf(o+strlen(o)," %4zu %11.2f %11.2f %11.2f %s \n",
                           icell,(float)entries,(float)t-t_offset,(float)st, GetCellName(icell).c_str() );
    }
    else
    {
      double t=0.;
      double st=1.;
      sprintf(o+strlen(o)," %4zu %11.2f %11.2f %11.2f %s \n",
                         icell,1.,(float)t-t_offset,(float)st, GetCellName(icell).c_str());
    }
  }
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputTimeCalibInfoXY4LED( string& s ) const
{
  assert( XYRegularGrid() );
// Check if all cells attached to SADC readout
  for( size_t it=0; it < NCells(); it++ )
  {
    if( GetLedDigitizerSADC(it) == NULL ) return 1;
  }

  char o[500000];
  o[0] = 0;
  sprintf(o," %s COMPASS Calorimeter Time Calibration File for LEDs \n",GetName().c_str());
  sprintf(o+strlen(o),"   X     Y        Time     SigmaTime  Statistic    \n");
  float t_offset=0.;
  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    if( GetLedDigitizerSADC(icell) == NULL ) continue;
    const CsDigitizerSADC * digl = dynamic_cast <const CsDigitizerSADC *>(GetLedDigitizerSADC(icell));
    if( digl == NULL )
    {
      cerr <<" Fatal error " << GetName() << " CsCalorimeter::OutputTimeCalibInfoXY4LED not yet implemented for DigitizerSADCN " << endl;
      exit(1);
    }
    const Reco::StatInfo &ledtime = digl->stat_time_;
    double entries = ledtime.GetEntries();
    if( entries > 10 )
    {
      double t=ledtime.GetMean();
      t_offset += t;
    }
    else
    {
      continue;
    }
  }

  t_offset = t_offset/float(NCells());
  sprintf(o+strlen(o)," Time_Offset %11.2f  \n",t_offset);

  for( int x=0; x != GetNColumns(); x++ )
    for( int y=0; y != GetNRows(); y++ )
    {
      int icell = GetCellOfColumnRow(x,y);
      if( icell >= 0 )
      {
        if( GetLedDigitizerSADC(icell) == NULL ) continue;
        const CsDigitizerSADC * digl = dynamic_cast <const CsDigitizerSADC *>(GetLedDigitizerSADC(icell));
        if( digl == NULL )
        {
          cerr <<" Fatal error " << GetName() << " CsCalorimeter::OutputTimeCalibInfoXY4LED not yet implemented for DigitizerSADCN " << endl;
          exit(1);
        }
        const Reco::StatInfo &ledtime = digl->stat_time_;
        double entries = ledtime.GetEntries();
        if( entries > 10 )
        {
          double t=ledtime.GetMean();
          double st=ledtime.GetSigma();
          sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,(float)t-t_offset,(float)st,entries);
        }
        else
        {
          double t=0.;
          double st=1.;
          sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,(float)t-t_offset,(float)st,1.);
        }
      }
      else
      {
        sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,-1.,-1.,-1.);
      }
    }
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputTrGrCalibInfo(const string &s)
{
  istringstream is(s.c_str());
  char calorim_name[132],dummy[132];
  string str;

//  Read line of comments
  getline(is,str);
  sscanf(str.c_str(),"%s %s ",calorim_name,dummy);

  getline(is,str);
  float energy_calib=0.;
  sscanf(str.c_str()," %s %g",dummy,&energy_calib);

//  Read 2 lines of comments
  getline(is,str);
  getline(is,str);

//  Read calibration
  while( getline(is,str) )
  {
    int l,x,y;
    float w,m,s,c,sc;
    sscanf(str.c_str()," %d %d %d %g %g %g %g %g \n",&l,&x,&y,&m,&s,&w,&c,&sc);
    int trig_gr = GetTrigGroupNumber(l,x,y);
    if( trig_gr >= 0 ) group_info[TrigGroup::CALIB][OLD][trig_gr]=Reco::StatInfo(w,c/1000.,sc/1000.);
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputTrGrCalibInfo(string& s) const
{
  if( NGroups() <= 0) goto out;
  char o[500000];
  o[0] = 0;
  sprintf(o," %s COMPASS Calorimeter TriggerGroups Energy Calibration File \n",GetName().c_str());
  sprintf(o+strlen(o)," Calibration_Energy(GeV)    %7.2f \n",options.calibration_energy );
  sprintf(o+strlen(o)," Layer  X     Y      MeanPeak   SigmaPeak  Statistic     Coeff    SigmaCoeff \n");
  sprintf(o+strlen(o),"                   (ADCchanels)                     (MeV/ADCchanel)          \n");
  for(size_t it=0; it!=NGroups(); it++)
  {
    int gr_layer = trigger_groups[it].trig_group_data.GetLayer();
    int gr_x     = trigger_groups[it].trig_group_data.GetX();
    int gr_y     = trigger_groups[it].trig_group_data.GetY();
    double entries = group_info[TrigGroup::CALIB][NEW][it].GetEntries();
    if( entries > 10 )
    {
      double cfnew=group_info[TrigGroup::CALIB][NEW][it].GetMean()*
                          group_info[TrigGroup::CALIB][OLD][it].GetMean();
      double scfnew=group_info[TrigGroup::CALIB][NEW][it].GetSigma()*
                           group_info[TrigGroup::CALIB][OLD][it].GetMean();
      double peak_adc=0;
      double sigma_peak_adc=0;
      if(options.calibration_energy > 0)
      {
        if(cfnew > 0 )
        {
          peak_adc = options.calibration_energy/cfnew;
          sigma_peak_adc = options.calibration_energy*scfnew/cfnew;
        }
      }
      sprintf(o+strlen(o)," %4d %4d %4d %11.2f %11.2f %11.2f %11.2f %11.2f \n",gr_layer,gr_x,gr_y,
                                peak_adc,sigma_peak_adc,entries,1000.*cfnew,1000.*scfnew);
    }
    else
    {
      double old_entries = group_info[TrigGroup::CALIB][OLD][it].GetEntries();
      double cfnew=group_info[TrigGroup::CALIB][OLD][it].GetMean();
      double scfnew=group_info[TrigGroup::CALIB][OLD][it].GetSigma();
      double peak_adc=0;
      double sigma_peak_adc=0;
      if(options.calibration_energy > 0)
      {
        if(cfnew > 0 )
        {
          peak_adc = options.calibration_energy/cfnew;
          sigma_peak_adc = options.calibration_energy*scfnew/cfnew;
        }
      }
      sprintf(o+strlen(o)," %4d %4d %4d %11.2f %11.2f %11.2f %11.2f %11.2f \n",gr_layer,gr_x,gr_y,
                              peak_adc,sigma_peak_adc,old_entries,1000.*cfnew,1000.*scfnew);
    }

  }
out:
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::InputTrGrTimeCalibInfo(const string &s)
{
  istringstream is(s.c_str());
  char calorim_name[132],dummy[132];
  string str;

//  Read line of comments
  getline(is,str);
  sscanf(str.c_str(),"%s %s ",calorim_name,dummy);

  getline(is,str);
  float time_average=0.;
  sscanf(str.c_str()," %s %g",dummy,&time_average);

// Set time_average as default
  for(size_t it=0; it!=NGroups(); it++)
  {
    group_info[TrigGroup::TIME1][OLD][it]=Reco::StatInfo(1,time_average,1.);
    group_info[TrigGroup::TIME2][OLD][it]=Reco::StatInfo(1,time_average,1.);
  }

//  Read 2 lines of comments
  getline(is,str);
  getline(is,str);

//  Read calibration
  while( getline(is,str) )
  {
    int l,x,y;
    float t1,st1,t2,st2;
    sscanf(str.c_str()," %d %d %d %g %g %g %g \n",&l,&x,&y,&t1,&st1,&t2,&st2);
    int trig_gr = GetTrigGroupNumber(l,x,y);
    if( trig_gr >= 0 )
    {
      group_info[TrigGroup::TIME1][OLD][trig_gr]=Reco::StatInfo(1,t1+time_average,st1);
      group_info[TrigGroup::TIME2][OLD][trig_gr]=Reco::StatInfo(1,t2+time_average,st2);
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::OutputTrGrTimeCalibInfo( string& s ) const
{
  char o[500000];
  o[0] = 0;
  if( NGroups() > 0)
  {
// First calculate average time
  double time_average_new=0;
  double time_average=0;
  double entries=0;
  for(size_t it=0; it!=NGroups(); it++)
  {
    double t1   = group_info[TrigGroup::TIME1][OLD][it].GetMean();
    double st1   = group_info[TrigGroup::TIME1][OLD][it].GetSigma();
    double entries1 = group_info[TrigGroup::TIME1][NEW][it].GetEntries();
    if( entries1 > 10 )
    {
      t1=group_info[TrigGroup::TIME1][NEW][it].GetMean()+
                          group_info[TrigGroup::TIME1][OLD][it].GetMean();
      st1=group_info[TrigGroup::TIME1][NEW][it].GetSigma();
      time_average_new += entries1*t1;
      entries += entries1;
    }

    double t2   = group_info[TrigGroup::TIME2][OLD][it].GetMean();
    double st2   = group_info[TrigGroup::TIME2][OLD][it].GetSigma();
    double entries2 = group_info[TrigGroup::TIME2][NEW][it].GetEntries();
    if( entries2 > 10 )
    {
      t2=group_info[TrigGroup::TIME2][NEW][it].GetMean()+
                          group_info[TrigGroup::TIME2][OLD][it].GetMean();
      st2=group_info[TrigGroup::TIME2][NEW][it].GetSigma();
      time_average_new += entries2*t2;
      entries += entries2;
    }
    time_average +=t1;
    time_average +=t2;
  }
  time_average /= 2*NGroups();
  if(entries >0 ) time_average = time_average_new/entries;
  sprintf(o," %s COMPASS Calorimeter TriggerGroups Time Calibration File \n",GetName().c_str());
  sprintf(o," Average_Time_constant(ns)    %7.2f \n",time_average );
  sprintf(o," Layer  X     Y         T1      SigmaT1        T2       SigmaT2 \n");
  sprintf(o,"                       (ns)      (ns)                           \n");
  for(size_t it=0; it!=NGroups(); it++)
  {
    int gr_layer = trigger_groups[it].trig_group_data.GetLayer();
    int gr_x     = trigger_groups[it].trig_group_data.GetX();
    int gr_y     = trigger_groups[it].trig_group_data.GetY();
    double t1   = group_info[TrigGroup::TIME1][OLD][it].GetMean();
    double st1   = group_info[TrigGroup::TIME1][OLD][it].GetSigma();
    double entries1 = group_info[TrigGroup::TIME1][NEW][it].GetEntries();
    if( entries1 > 10 )
    {
      t1=group_info[TrigGroup::TIME1][NEW][it].GetMean()+
                          group_info[TrigGroup::TIME1][OLD][it].GetMean();
      st1=group_info[TrigGroup::TIME1][NEW][it].GetSigma();

    }
    double t2   = group_info[TrigGroup::TIME2][OLD][it].GetMean();
    double st2   = group_info[TrigGroup::TIME2][OLD][it].GetSigma();
    double entries2 = group_info[TrigGroup::TIME2][NEW][it].GetEntries();
    if( entries2 > 10 )
    {
      t2=group_info[TrigGroup::TIME2][NEW][it].GetMean()+
                          group_info[TrigGroup::TIME2][OLD][it].GetMean();
      st2=group_info[TrigGroup::TIME2][NEW][it].GetSigma();
    }
    sprintf(o," %4d %4d %4d %11.2f %11.2f %11.2f %11.2f \n",gr_layer,gr_x,gr_y,
                                          t1-time_average,st1,t2-time_average,st2  );
  }
  }
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::MakeMCDigitization(void)
{
  bool debug = false;

  if( debug ) cout << " CsCalorimeter::MakeMCDigitization " << GetName() << " debug " << endl;

  if( debug ) cout << " call Reco::Calorimeter::MakeMCDigitization " << endl;
  if( debug ) cout << " Before MakeMCDigitizatio signals_ " << signals_.size() << endl;
  Reco::Calorimeter::MakeMCDigitization();
  if( debug ) cout << " After MakeMCDigitizatio signals_ " << signals_.size() << endl;

  if( debug ) cout << " CsCalorimeter::MakeMCDigitization " << GetName() << " debug  OK " << endl;

}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::makeMCDecoding(void)
{
  // fill the timing of the trigger groups
  // TODO: smearing of the time?
  for( size_t it=0; it!=NGroups(); it++ )
  {
    double gr_sumE = 0.;
    double gr_time = 0.;

    for (vector<CellDataRaw>::const_iterator itMC=mc_input_.begin(); itMC!=mc_input_.end(); itMC++) {
      for (vector<size_t>::const_iterator it2=trigger_groups[it].GetCells().begin();
           it2!=trigger_groups[it].GetCells().end(); it2++) {
        if ( itMC->GetCellIdx()==(*it2) ) {
          gr_sumE += itMC->GetEnergy();
          gr_time += itMC->GetEnergy() * itMC->GetTime();
        }
      }
    }

    if(gr_sumE > 0. )
    {
      trigger_groups[it].trig_group_data.SetTime1(gr_time/gr_sumE);
      trigger_groups[it].trig_group_data.SetTime2(gr_time/gr_sumE);
    }
  }

  MakeMCDigitization();

  // TODO: setting the error on the time with a proper value
  //       at the moment this is done as for real data, the
  //       error is simply set to the SADC clock
  for (vector<CellDataRaw>::iterator itSig=signals_.begin(); itSig!=signals_.end(); itSig++)
    itSig->SetTimeErr(12.86);

  // Calculate summ ADC in groups
  for( size_t it=0; it!=NGroups(); it++ )
  {
    double gr_sumE = 0.;
    double gr_time = 0.;
    for (vector<CellDataRaw>::const_iterator itSig=signals_.begin(); itSig!=signals_.end(); itSig++) {
      for(vector<size_t>::const_iterator it2=trigger_groups[it].GetCells().begin();
          it2!=trigger_groups[it].GetCells().end(); it2++ ) {
         if( itSig->GetCellIdx()==(*it2) )
           gr_sumE += itSig->GetEnergy();
      }
    }

    trigger_groups[it].trig_group_data.SetSumADC(gr_sumE);
    trigger_groups[it].trig_group_data.SetAmpl(gr_sumE);
  }

  MakeCsDigits();
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::MakeCsDigits( void )
{
  bool debug = false;
  if( skip_make_cs_digits_ ) return;

  if( debug ) cout << " CsCalorimeter::MakeCsDigits " << GetName() <<
                       " signals_ size=" << signals_.size() <<
                                       " N CsDigist=" << myDigits_.size() << endl;
  for( size_t it=0; it!=signals_.size(); it++ )
  {
    if( add_timeinfo_to_cs_digits_ )
    {
      double data[3];
      int icell = signals_[it].GetCellIdx();
      data[0]=signals_[it].GetEnergy();
      data[1]=signals_[it].GetTime();
      if( debug ) cout << " it " << it << " cell " << icell << " data " << data[0] << endl;

      CsDigit *d = new CsDigit( *this, icell, data, 2 );
      myDigits_.push_back(d);

    }
    else if( add_sadcinfo_to_cs_digits_ || add_all_sadcinfo_to_cs_digits_ )
    {
      double data[1002];
      int icell = signals_[it].GetCellIdx();
      data[0]=signals_[it].GetEnergy();
      data[1]=signals_[it].GetTime();
      if( debug ) cout << " it " << it << " cell " << icell << " data " << data[0] << endl;
      const vector<CS::uint16> &sample = *GetSADCSample(icell);
      if( sample.size() < 1000 )
      {
        for( int i=0; i < (int)sample.size(); i++ )
        {
          data[2+i]=sample[i];
        }
      }
      CsDigit *d = new CsDigit( *this, icell, data, 2+sample.size() );
      myDigits_.push_back(d);

    }
    else if( CsEvent::Instance()->isAMonteCarloEvent() )
    {
      int icell = signals_[it].GetCellIdx();
      double data[2];
      data[0]=signals_[it].GetEnergy();
      data[1]=signals_[it].GetTime();

      CsDigit *d = new CsDigit(*this, icell, data, 2);
      myDigits_.push_back(d);
    }
    else
    {
      double data[3];
      int icell = signals_[it].GetCellIdx();
      data[0]=signals_[it].GetEnergy();
      if( debug ) cout << " it " << it << " cell " << icell << " data " << data[0] << endl;

      CsDigit *d = new CsDigit( *this, icell, data, 1 );
      myDigits_.push_back(d);
    }
  }

  if( add_all_sadcinfo_to_cs_digits_ )
  {
    for(size_t icell=0; icell<NCells(); icell++ )
    {
      const vector<CS::uint16> *asample = GetSADCSample(icell);
      if( asample == NULL ) continue;
      size_t it =0;
      for( it=0; it!=signals_.size(); it++ )
      {
        if( icell == signals_[it].GetCellIdx()) break;
      }
      if( it != signals_.size() ) continue;

      const vector<CS::uint16> &sample = *asample;
      double data[1002];
      data[0]=0.;
      data[1]=-1000000.;
      if( sample.size() < 1000 )
      {
        for( int i=0; i < (int)sample.size(); i++ )
        {
          data[2+i]=sample[i];
        }
      }
      CsDigit *d = new CsDigit( *this, icell, data, 2+sample.size() );
      myDigits_.push_back(d);
    }
  }

  if( debug ) cout << " CsCalorimeter::MakeCsDigits on output N CsDigist=" << myDigits_.size() <<endl;
}

////////////////////////////////////////////////////////////////////////////////

bool CsCalorimeter::Reconstruction(void)
{
  bool debug = false;
  if( debug ) cout << " This is CsCalorimeter  Reconstruction for " << GetName() <<
                                               " nrec= " << reco_particles.size() << endl;
  if( debug ) cout << "  signals_.size() = " << signals_.size() << endl;


  // do the reconstruction
  bool ok=Reco::Calorimeter::Reconstruction();
  if( debug )
  {
    cout << " Reco::Calorimeter::Reconstruction  OK " << endl;
    cout << " Nrec " << reco_particles.size() << endl;
    for( unsigned i=0; i < reco_particles.size(); i++ )
    {
      cout << i << " E " << reco_particles[i].GetE() << endl;
    }
  }

  {
// Add SADC timing
    bool exit_if_abnormal_energy_detected = false;
    for( size_t ip = 0; ip < reco_particles.size(); ip++ )
    {
      const vector<size_t> &hitedcells = reco_particles[ip].GetMainCells();
      if( debug ) cout << ip << " CsCalorimeter::Reconstruction  main cells vector size = " << hitedcells.size() << endl;
      if( reco_particles[ip].GetE() > 1000. )
      {
        cerr << " CsCalorimeter::Reconstruction abnormal energy " << reco_particles[ip].GetE() <<" was detected " << endl;
        if( exit_if_abnormal_energy_detected )
        {
          cerr << " Exit for debug " << endl;
          exit(1);
        }
      }

      // Set SADC Shape flag
      if(!hitedcells.empty())
      {
        size_t icell = hitedcells[0];
	const DigitizerSADCBase * dig = GetDigitizerSADC ( icell );
	if( dig != NULL )
	{
          bool good_sadc_shape = !dig->IsNoiseByShape();
	  reco_particles[ip].SetMiscInfo(Reco::CalorimeterParticle::SHAPE_SADC_OK, double(good_sadc_shape) );
	}
      }

      // do not set the time for particles which already haven one
      // for example from the shower fit
      if (reco_particles[ip].HasTime())
        continue;

      if(!hitedcells.empty())
      {
        const size_t icell = hitedcells[0]; // main cell
        for (size_t i=0; i<signals_.size(); i++)
          if (signals_[i].GetCellIdx() == icell) {
            if (signals_[i].HasTime()) {
              // in principle one could check above that also a time error is set,
              // however, an error on the time should always be given if a time is
              // set, so trigger a warning if this is not the case.
              const double time = signals_[i].GetTime();
              const double timeErr = signals_[i].GetTimeErr();

              reco_particles[ip].SetTime(time, timeErr);
              if (std::abs(time) > 10000.) {
                const double ecell = signals_[i].GetEnergy();
                cerr << " CsCalorimeter::Reconstruction Abnormal time detected E " << reco_particles[ip].GetE()
                     << " E cell " << ecell
                     << " Time " << time <<" tcs_phase " << tcs_phase_ <<  endl;
              }
            }
            break;
          }
      }
      else
      {
        cerr << " Very, very BAD! Csalorimeter::Reconstruction main cells vector is empty !!!" << ip << endl;
      }
    }
  }

//   cout << " CsDigits in " << GetName() <<  " = " << getMyDigits().size() << endl;

  return ok;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::CreateGUI()
{
  if( gui!=NULL )
  {
    CsErrLog::msg(elFatal, __FILE__, __LINE__,
                  "CsCalorimeter::CreateGUI: attempted to create a second GUI for %s.",
                  getName());
    return;
  }

  gui = new CsCalorimeterGUI(*this,NULL,GetName().c_str());
}

////////////////////////////////////////////////////////////////////////////////

double vvalue( int index, const vector<CS::uint16>& v )
{
  if( index >= 0 && index < (int)v.size() )
  {
    return (double)v[index];
  }
  else
  {
    return 0.;
  }
}

////////////////////////////////////////////////////////////////////////////////

double vvalue( int index, const vector<double>& v )
{
  if( index >= 0 && index < (int)v.size() )
  {
    return (double)v[index];
  }
  else
  {
    return 0.;
  }
}

////////////////////////////////////////////////////////////////////////////////

const vector<CS::uint16>* CsCalorimeter::GetSADCSample ( int icell ) const
{
  map< int, const std::vector<CS::uint16> *  >::const_iterator e = sadc_samples_.find(icell);
  if( e!= sadc_samples_.end() ) return e->second;
  return NULL;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::PrintSADCSettings ( void ) const
{
  cout << " SADC settings in " << GetName();
  cout << " sadc_format_version_ " << sadc_format_version_;
  cout << " sadc_decode_version_ " << sadc_decode_version_;
  cout << " sadc_front_ " << sadc_front_;
  cout << " sadc_front_gate_ " << sadc_front_gate_;
  cout << " sadc_ped_min_ " << sadc_ped_min_;
  cout << " sadc_ped_max_ " << sadc_ped_max_;
  cout << " sadc_signal_min_ " << sadc_signal_min_;
  cout << " sadc_signal_max_ " << sadc_signal_max_;
  cout << " sadc_delta_ " << sadc_delta_;
  cout << " sadc_clock_ " << sadc_clock_;
  cout << " coeff_convert2max_ " << coeff_convert2max_;
  cout << " sadc_max_position_ " << sadc_max_position_;
  cout << " make_tcs_corrections_ " << make_tcs_corrections_;
  cout << endl;
}

////////////////////////////////////////////////////////////////////////////////

vector<CsTrigGroupData> CsCalorimeter::GetTrigGroupData  (double energy_threshold) const
{
  vector<CsTrigGroupData> trig_group_data;
  for( size_t it=0; it!=NGroups(); it++ )
    if(trigger_groups[it].trig_group_data.GetSumADC() >= energy_threshold )
      trig_group_data.push_back(trigger_groups[it].trig_group_data);
  return trig_group_data;
}

////////////////////////////////////////////////////////////////////////////////

vector<CsTrigGroupData> CsCalorimeter::GetTrigGroupDataXY (double energy_threshold,
                                                               double x,double y) const
{
  vector<CsTrigGroupData> trig_group_data;

  int icell = FindCellInternal(x-GetPositionX(),y-GetPositionY(),0);
//  int icell = FindCell(x-GetPositionX(),y-GetPositionY(),0);
  if(icell < 0) return trig_group_data;
  for( size_t it=0; it!=NGroups(); it++ )
  {
     for( vector<size_t>::const_iterator it1=trigger_groups[it].GetCells().begin();
                                         it1!=trigger_groups[it].GetCells().end(); it1++ )
     {
       if(icell == (int)*it1)
       {
          if(trigger_groups[it].trig_group_data.GetSumADC() >= energy_threshold )
            trig_group_data.push_back(trigger_groups[it].trig_group_data);
          break;
       }
     }
  }
  return trig_group_data;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::CalibrateTrigGroupTime(void)
{
  for( size_t it=0; it!=NGroups(); it++ )
  {
    if( trigger_groups[it].trig_group_data.IsTime1Valid() )
    {
      double time = trigger_groups[it].trig_group_data.GetTime1();
      if(fabs(time) < 5) group_info[TrigGroup::TIME1][NEW][it].Add(time,1.);
    }
    if( trigger_groups[it].trig_group_data.IsTime2Valid() )
    {
      double time = trigger_groups[it].trig_group_data.GetTime2();
      if(fabs(time) < 5) group_info[TrigGroup::TIME2][NEW][it].Add(time,1.);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

int CsCalorimeter::GetTrigGroupNumber(int layer,int x,int y) const
{
  for(size_t it=0; it!=NGroups(); it++)
  {
    if(trigger_groups[it].trig_group_data.GetLayer() == layer &&
               trigger_groups[it].trig_group_data.GetX() == x &&
               trigger_groups[it].trig_group_data.GetY() == y) return it;
  }
  return -1;
}

////////////////////////////////////////////////////////////////////////////////

CsCalorimeter::TrigGroup::TrigGroup   (Calorimeter *c, const string &tg_name, const vector<size_t> &tg,
                                                                   int g_l,int g_x,int g_y) :
  Reco::Calorimeter::SubSet(c,tg_name,tg),
  trig_group_data( g_l, g_x, g_y)
{}

////////////////////////////////////////////////////////////////////////////////

CsCalorimeter::TrigGroup::TrigGroup ( const  CsCalorimeter::TrigGroup &tg) :
  Reco::Calorimeter::SubSet( tg.calorimeter,tg.name,tg.GetCells()),
  trig_group_data (tg.trig_group_data)
{}

////////////////////////////////////////////////////////////////////////////////

 CsCalorimeter::TrigGroup & CsCalorimeter::TrigGroup::operator = (const  CsCalorimeter::TrigGroup &tg)
{
  if( &tg!=this )
  {
    trig_group_data = tg.trig_group_data;
    name            = tg.name;
    sub_set_cells   = tg.sub_set_cells;
    info_           = tg.info_;
    finfo_          = tg.finfo_;
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

 CsTrigGroupData & CsTrigGroupData::operator = (const  CsTrigGroupData &tg)
{
  if( &tg!=this )
  {
    fAmpl           = tg.fAmpl;
    fTime1          = tg.fTime1;
    fTime2          = tg.fTime2;
    fSumADC         = tg.fSumADC;
    time1_is_valid  = tg.time1_is_valid;
    time2_is_valid  = tg.time2_is_valid;
    fId_layer       = tg.fId_layer;
    fId_x           = tg.fId_x;
    fId_y           = tg.fId_y;
  }
  return *this;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DigitizerSADCBase * CsCalorimeter::GetDigitizerSADC2Modify ( size_t icell )
{
  if(dig_sadc_!=NULL && GetDigitizersSADC().size() == NCells() && icell < NCells() ) return GetDigitizersSADC()[icell];
  return NULL;
}

///////////////////////////////////////////////////////////////////////////////////////////////

const DigitizerSADCBase * CsCalorimeter::GetDigitizerSADC ( size_t icell ) const
{
  if(dig_sadc_!=NULL && GetDigitizersSADCs().size() == NCells() && icell < NCells() ) return GetDigitizersSADCs()[icell];
  return NULL;
}

///////////////////////////////////////////////////////////////////////////////////////////////

const DigitizerSADCBase * CsCalorimeter::GetLedDigitizerSADC ( size_t icell ) const
{
  if(dig_led_sadc_!=NULL && GetLedDigitizersSADCs().size() == NCells() && icell < NCells() ) return GetLedDigitizersSADCs()[icell];
  return NULL;
}

///////////////////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::AddDigitizerSADC( DigitizerSADCBase * dig, size_t icell )
{
  if(GetDigitizersSADC().size() == NCells() && icell < NCells() )
  {
    GetDigitizersSADC()[icell] = dig;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::AddLedDigitizerSADC( DigitizerSADCBase * dig, size_t icell )
{
  if( GetLedDigitizersSADC().size() == NCells() && icell < NCells() )
  {
    GetLedDigitizersSADC()[icell] = dig;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////

std::vector<DigitizerSADCBase *> &CsCalorimeter::GetDigitizersSADC( void )
{
  return dig_sadc_->vsadc;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<DigitizerSADCBase *> &CsCalorimeter::GetLedDigitizersSADC( void )
{
  return dig_led_sadc_->vsadc;
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<DigitizerSADCBase *> &CsCalorimeter::GetDigitizersSADCs( void ) const {return (dig_sadc_->vsadc);}

////////////////////////////////////////////////////////////////////////////////

const std::vector<DigitizerSADCBase *> &CsCalorimeter::GetLedDigitizersSADCs( void ) const {return (dig_led_sadc_->vsadc);}

////////////////////////////////////////////////////////////////////////////////

void  CsCalorimeter::InitArraySADC( void )
{
  bool debug = false;
  if( debug ) std::cout <<" CsCalorimeter::InitArraySADC " << sadc_decode_version_ << std::endl;
  dig_sadc_ = new DigitizerStore(NCells());
  dig_led_sadc_ = new DigitizerStore(NCells());
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::ClearDigitizersSADC (void)
{
  if(dig_sadc_ != NULL)
    dig_sadc_->Clear();
}

////////////////////////////////////////////////////////////////////////////////

void DigitizerStore::Clear (void)
{
  for( size_t ic=0; ic< vsadc.size(); ic++)
  {
    if(vsadc[ic] != NULL ) vsadc[ic]->Clear();
  }
}

////////////////////////////////////////////////////////////////////////////////
