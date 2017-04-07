/*!
   \file    CsECAL1.cc
   \brief   Base class for COMPASS calorimeters
   \version $Revision: 1.48 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \date    $Date: 2011/02/01 22:05:52 $
*/

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "CsOpt.h"
#include "CsErrLog.h"

#include "DaqDataDecoding/ChipADC.h"
#include "DaqDataDecoding/ChipSADC.h"
#include "DaqDataDecoding/ChipF1.h"

#include "CsECAL1.h"
#include "CsCalorimeterGUI.h"
#include "CsDigitizerSADC.h"

#include "Reco/CellDataRaw.h"
#include "Reco/Exception.h"

using namespace std;
using namespace Reco;

////////////////////////////////////////////////////////////////////////////////

CsECAL1::CsECAL1(const string &name, const string &geom_file)
 : CsCalorimeter(name, geom_file),
  index_gams_(0),
  index_maintz_top_(0),
  index_maintz_bottom_(0),
  index_olga_saleve_(0),
  index_olga_jura_(0),
  new_style_sadc_timing_(false)
{
  bool debug = false;
  if( debug ) cout << " CsECAL1::CsECAL1 debug " << endl;

  ClearSubSet();

  double gams_size = 38.;
  double maintz_size = 75.;
  double olga_size = 143.;

// GAMS cells first
  {
    vector < size_t > cells_list;
    for( size_t it=0; it!=NCells(); it++ )
    {
      double sx = cells[it].GetCellType().GetSizeX();
      if( fabs( sx - gams_size ) < 5. )
      {
        cells_list.push_back(it);
      }
    }

    const char* name = "EC01P1_GAMS";
    AddSubSet( SubSet( this, name,cells_list),true );
    cell_sizex_gams_ = cells[cells_list[0]].GetCellType().GetSizeX();
    cell_sizey_gams_ = cells[cells_list[0]].GetCellType().GetSizeY();

  }
  index_gams_=GetSubSets().size()-1;
// Maintz top cells
  {
    vector < size_t > cells_list;
    for( size_t it=0; it!=NCells(); it++ )
    {
      double sx = cells[it].GetCellType().GetSizeX();
      if( fabs( sx - maintz_size ) < 5. && cells[it].GetY() > 0. )
      {
        cells_list.push_back(it);
      }
    }

    const char* name = "EC01P1_MaintzTop";
    AddSubSet( SubSet( this, name,cells_list),true );
    cell_sizex_maintz_ = cells[cells_list[0]].GetCellType().GetSizeX();
    cell_sizey_maintz_ = cells[cells_list[0]].GetCellType().GetSizeY();
  }
  index_maintz_top_=GetSubSets().size()-1;
// Maintz bottom cells
  {
    vector < size_t > cells_list;
    for( size_t it=0; it!=NCells(); it++ )
    {
      double sx = cells[it].GetCellType().GetSizeX();
      if( fabs( sx - maintz_size ) < 5. && cells[it].GetY() < 0. )
      {
        cells_list.push_back(it);
      }
    }

    const char* name = "EC01P1_MaintzBottom";
    AddSubSet( SubSet( this, name,cells_list),true );
    cell_sizex_maintz_ = cells[cells_list[0]].GetCellType().GetSizeX();
    cell_sizey_maintz_ = cells[cells_list[0]].GetCellType().GetSizeY();
  }
  index_maintz_bottom_=GetSubSets().size()-1;
// Olga Saleve cells
  {
    vector < size_t > cells_list;
    for( size_t it=0; it!=NCells(); it++ )
    {
      double sx = cells[it].GetCellType().GetSizeX();
      if( fabs( sx - olga_size ) < 5.&& cells[it].GetX() < 0.)
      {
        cells_list.push_back(it);
      }
    }

    const char* name = "EC01P1_OlgaSaleve";
    AddSubSet( SubSet( this, name,cells_list),true );
    cell_sizex_olga_ = cells[cells_list[0]].GetCellType().GetSizeX();
    cell_sizey_olga_ = cells[cells_list[0]].GetCellType().GetSizeY();
  }
  index_olga_saleve_=GetSubSets().size()-1;
// Olga Jura cells
  {
    vector < size_t > cells_list;
    for( size_t it=0; it!=NCells(); it++ )
    {
      double sx = cells[it].GetCellType().GetSizeX();
      if( fabs( sx - olga_size ) < 5.&& cells[it].GetX() > 0. )
      {
        cells_list.push_back(it);
      }
    }

    const char* name = "EC01P1_OlgaJura";
    AddSubSet( SubSet( this, name,cells_list),true );
    cell_sizex_olga_ = cells[cells_list[0]].GetCellType().GetSizeX();
    cell_sizey_olga_ = cells[cells_list[0]].GetCellType().GetSizeY();
  }
  index_olga_jura_=GetSubSets().size()-1;

  if( !calorimeter_sub_sets[index_gams_].InitXY( 2.) )
  {
    cerr << " Error in " << GetName() << " subset " << calorimeter_sub_sets[index_gams_].GetName() << "  not XY structure !! " << endl;
    exit(1);
  }
  else
  {
    if( debug )
    {
      cout << " subset " << calorimeter_sub_sets[index_gams_].GetName() << "  XY structure !! " << endl;
      cout << " ncol " <<calorimeter_sub_sets[index_gams_].GetNColumns() << " nrows " << calorimeter_sub_sets[index_gams_].GetNRows() <<
            " xstep " << calorimeter_sub_sets[index_gams_].GetXStep() << " ystep " << calorimeter_sub_sets[index_gams_].GetYStep() << endl;
    }
  }

  if( !calorimeter_sub_sets[index_maintz_top_].InitXY( 6.) )
  {
    cerr << " Error in " << GetName() << " subset " << calorimeter_sub_sets[index_maintz_top_].GetName() << "  not XY structure !! " << endl;
    exit(1);
  }
  else
  {
    if( debug )
    {
      cout << " subset " << calorimeter_sub_sets[index_maintz_top_].GetName() << "  XY structure !! " << endl;
      cout << " ncol " <<calorimeter_sub_sets[index_maintz_top_].GetNColumns() << " nrows " << calorimeter_sub_sets[index_maintz_top_].GetNRows() <<
            " xstep " << calorimeter_sub_sets[index_maintz_top_].GetXStep() << " ystep " << calorimeter_sub_sets[index_maintz_top_].GetYStep() << endl;
    }
  }

  if( !calorimeter_sub_sets[index_maintz_bottom_].InitXY( 6.) )
  {
    cerr << " Error in " << GetName() << " subset " << calorimeter_sub_sets[index_maintz_bottom_].GetName() << "  not XY structure !! " << endl;
    exit(1);
  }
  else
  {
    if( debug )
    {
      cout << " subset " << calorimeter_sub_sets[index_maintz_bottom_].GetName() << "  XY structure !! " << endl;
      cout << " ncol " <<calorimeter_sub_sets[index_maintz_bottom_].GetNColumns() << " nrows " << calorimeter_sub_sets[index_maintz_bottom_].GetNRows() <<
            " xstep " << calorimeter_sub_sets[index_maintz_bottom_].GetXStep() << " ystep " << calorimeter_sub_sets[index_maintz_bottom_].GetYStep() << endl;
    }
  }

  if( !calorimeter_sub_sets[index_olga_saleve_].InitXY( 5.) )
  {
    cerr << " Error in " << GetName() << " subset " << calorimeter_sub_sets[index_olga_saleve_].GetName() << "  not XY structure !! " << endl;
    exit(1);
  }
  else
  {
    if( debug )
    {
      cout << " subset " << calorimeter_sub_sets[index_olga_saleve_].GetName() << "  XY structure !! " << endl;
      cout << " ncol " <<calorimeter_sub_sets[index_olga_saleve_].GetNColumns() << " nrows " << calorimeter_sub_sets[index_olga_saleve_].GetNRows() <<
            " xstep " << calorimeter_sub_sets[index_olga_saleve_].GetXStep() << " ystep " << calorimeter_sub_sets[index_olga_saleve_].GetYStep() << endl;
    }
  }

  if( !calorimeter_sub_sets[index_olga_jura_].InitXY( 5.) )
  {
    cerr << " Error in " << GetName() << " subset " << calorimeter_sub_sets[index_olga_jura_].GetName() << "  not XY structure !! " << endl;
    exit(1);
  }
  else
  {
    if( debug )
    {
      cout << " subset " << calorimeter_sub_sets[index_olga_jura_].GetName() << "  XY structure !! " << endl;
      cout << " ncol " <<calorimeter_sub_sets[index_olga_jura_].GetNColumns() << " nrows " << calorimeter_sub_sets[index_olga_jura_].GetNRows() <<
            " xstep " << calorimeter_sub_sets[index_olga_jura_].GetXStep() << " ystep " << calorimeter_sub_sets[index_olga_jura_].GetYStep() << endl;
    }
  }

  if( debug ) cout << " CsECAL1::CsECAL1 debug OK " << endl;
}

////////////////////////////////////////////////////////////////////////////////

void  CsECAL1::Initialize( void )
{
  CsCalorimeter::Initialize();

  if( options.mixed_blocks ) {
    InitXY();
  }
  else
  {
    cerr <<" Wrong settings!! ECAL1 has a mixed blocks structure " << endl;
    options.mixed_blocks = true;
    InitXY();
  }

  double gams_cells_factor = 0.610/0.7;
  double maintz_cells_factor = 0.480/0.7;
  double olga_cells_factor = 0.750/0.7;


  for ( size_t i=0; i<calorimeter_sub_sets[index_gams_].Size(); i++ )
  {
    size_t icell = calorimeter_sub_sets[index_gams_].GetCells()[i];
    individual_calib_factor_[icell] = gams_cells_factor;
  }

  for ( size_t i=0; i<calorimeter_sub_sets[index_maintz_top_].Size(); i++ )
  {
    size_t icell = calorimeter_sub_sets[index_maintz_top_].GetCells()[i];
    individual_calib_factor_[icell] = maintz_cells_factor;
  }

  for ( size_t i=0; i<calorimeter_sub_sets[index_maintz_bottom_].Size(); i++ )
  {
    size_t icell = calorimeter_sub_sets[index_maintz_bottom_].GetCells()[i];
    individual_calib_factor_[icell] = maintz_cells_factor;
  }

  for ( size_t i=0; i<calorimeter_sub_sets[index_olga_saleve_].Size(); i++ )
  {
    size_t icell = calorimeter_sub_sets[index_olga_saleve_].GetCells()[i];
    individual_calib_factor_[icell] = olga_cells_factor;
  }

  for ( size_t i=0; i<calorimeter_sub_sets[index_olga_jura_].Size(); i++ )
  {
    size_t icell = calorimeter_sub_sets[index_olga_jura_].GetCells()[i];
    individual_calib_factor_[icell] = olga_cells_factor;
  }

  subsets_list2store_sadc_histo_.push_back(index_gams_);
  subsets_list2store_sadc_histo_.push_back(index_maintz_top_);
  subsets_list2store_sadc_histo_.push_back(index_maintz_bottom_);
  subsets_list2store_sadc_histo_.push_back(index_olga_saleve_);
  subsets_list2store_sadc_histo_.push_back(index_olga_jura_);

  subsets_list2store_calib_histo_.push_back(index_gams_);
  subsets_list2store_calib_histo_.push_back(index_maintz_top_);
  subsets_list2store_calib_histo_.push_back(index_maintz_bottom_);
  subsets_list2store_calib_histo_.push_back(index_olga_saleve_);
  subsets_list2store_calib_histo_.push_back(index_olga_jura_);
  if( options.print_general_info ) PrintGeneralInfo();
}

////////////////////////////////////////////////////////////////////////////////

void  CsECAL1::InitOptions( void )
{
  CsCalorimeter::InitOptions();

  options.mixed_blocks                  = true;
  options.tolerance_for_nearby_cells    = 15.0;
  options.use_preshower_correction      = false;
  options.preshower_length              = 3.14;
  options.ecut4mgg_histo                = 1.5;
  options.default_calibration           = 0.060;
  options.readout_sparsified_           = true;
  options.delta_spars_mode_             = 3.;

  // MC options
  options.mc_default_calibration     = 0.060;
  options.mc_make_real_digitization  = true;
  options.mc_make_fiadc_digitization = true;
  options.mc_fiadc_sparce_delta      = 3.;

  options.SetMiscDoubleOption(1, 5.6 );
  if( CsOpt::Instance()->getOpt( GetName(), "new_style_sadc_timing" ) )
  {
    new_style_sadc_timing_ = true;
    cout << "WARNING!! In CsCalorimeter " << getName() << " make_sadc_histo option was set!!!! " << endl;
  }
  InitPreshowerCorrections();
}

////////////////////////////////////////////////////////////////////////////////

void  CsECAL1::InitReadOut( void )
{
  bool debug = false;
// new SADC stuff
  CsCalorimeter::GetSADCDigitizationOptions(); // Move this call to implementation classes
  if( debug ) std::cout <<" CsECAL1::InitReadOut sadc_decode_version_ " << sadc_decode_version_ << std::endl;
  if( sadc_decode_version_ != 1 ) return CsCalorimeter::InitReadOut();  // Jump to the base class implementation
  InitArraySADC();

// More SADC development

  if( debug ) cout << " CsECAL1:: create SADC digitazers " <<  endl;
  for ( unsigned ic=0; ic<NCells(); ic++ )
  {
    CsDigitizerSADC::Shape rd_table = CsDigitizerSADC::RD_ECAL1_GAMS;
    CsDigitizerSADC::Shape led_table = CsDigitizerSADC::LED_ECAL1_GAMS;
//    index_olga_jura_
    if( calorimeter_sub_sets[index_olga_jura_].IsMyCell(ic) )
    {
      rd_table = CsDigitizerSADC::RD_ECAL1_OLGA;
      led_table = CsDigitizerSADC::LED_ECAL1_OLGA;
    }
    else if( calorimeter_sub_sets[index_olga_saleve_].IsMyCell(ic) )
    {
      rd_table = CsDigitizerSADC::RD_ECAL1_OLGA;
      led_table = CsDigitizerSADC::LED_ECAL1_OLGA;
    }
    else if( calorimeter_sub_sets[index_maintz_bottom_].IsMyCell(ic) )
    {
      rd_table = CsDigitizerSADC::RD_ECAL1_MAINZ;
      led_table = CsDigitizerSADC::LED_ECAL1_MAINZ;
    }
    else if( calorimeter_sub_sets[index_maintz_top_].IsMyCell(ic) )
    {
      rd_table = CsDigitizerSADC::RD_ECAL1_MAINZ;
      led_table = CsDigitizerSADC::LED_ECAL1_MAINZ;
    }
    else
    {
      rd_table = CsDigitizerSADC::RD_ECAL1_GAMS;
      led_table = CsDigitizerSADC::LED_ECAL1_GAMS;
    }

    CsDigitizerSADC::Type type_sadc = CsDigitizerSADC::SADC;

    CsDigitizerSADC *d = new CsDigitizerSADC( this, (CsDigitizerSADC::Method)options.digitizer_method, type_sadc, rd_table );
    if( d == NULL )
    {
      cerr << " Memory problems in CsECAL1::CsECAL1:: ?? " << endl;
      exit(1);
    }
    AddDigitizerSADC(d,ic);

    CsDigitizerSADC *dl= new CsDigitizerSADC( this , (CsDigitizerSADC::Method)options.digitizer_method_led, type_sadc, led_table );
    if( dl == NULL )
    {
      cerr << " Memory problems in CsECAL1::CsECAL1:: ?? " << endl;
      exit(1);
    }
    dl->option_store_stat_info_ = true;
    AddLedDigitizerSADC(dl,ic);
  }
  if( debug ) cout << " CsECAL1:: create SADC digitazers OK " << endl;

// end new SADC stuff
  if( debug ) cout << " CsECAL1::CsECAL1 goto SADCDigitization settings " << endl;
  GetSADCDigitizationOptions();

// Add more to SADC digitization
  assert( GetDigitizersSADC().size() == NCells() );
  assert( GetLedDigitizersSADC().size() == NCells() );
  for( unsigned ic = 0; ic < NCells(); ic++ )
  {
    assert( GetDigitizersSADC()[ic] != NULL );
    assert( GetLedDigitizersSADC()[ic] != NULL );
    CsDigitizerSADC * dig =  dynamic_cast<CsDigitizerSADC *>(GetDigitizersSADC()[ic]);
    CsDigitizerSADC * digl =  dynamic_cast<CsDigitizerSADC *>(GetLedDigitizersSADC()[ic]);
    if( dig == NULL )
    {
       cerr <<" CsECAL1::InitReadOut CsDigitizerSADC = NULL some internal error in GetDigitizersSADC() initialization " << GetCellName(ic) << endl;
       exit(1);
    }
    if( digl == NULL )
    {
       cerr <<" CsECAL1::InitReadOut CsDigitizerSADC = NULL some internal error in GetLedDigitizersSADC() initialization " << GetCellName(ic) <<  endl;
       exit(1);
    }

    if( debug ) cout << " Start SADC settings for the cell " << ic <<  endl;
    dig->use_shape_filter_line_fit_=sadc_shape_filter_line_fit_use_;

    vector < double > vpar;
    switch( dig->GetMethod() ) {
    case CsDigitizerSADC::MaxSimple:
      vpar.push_back(1.);
      vpar.push_back(sadc_ped_min_);
      vpar.push_back(sadc_ped_max_);
      break;
    case CsDigitizerSADC::SummRange:
      vpar.push_back(coeff_convert2max_);
      vpar.push_back(sadc_ped_min_);
      vpar.push_back(sadc_ped_max_);
      vpar.push_back(sadc_signal_min_);
      vpar.push_back(sadc_signal_max_);
      vpar.push_back(sadc_front_);
      vpar.push_back(sadc_front_gate_);
      break;
    default:
      break;
    }

    if( !vpar.empty() ) {

      if( debug )
        cout << " Load dig_ settings " << dig << endl;

      bool ok = dig->LoadSettings(vpar);

      if( !ok )
        cerr << " Bad SADC settings for the cell " << ic << " Calorimeter " << GetName() << " take actions pls. ! " << endl;

      if( debug )
        cout << " Load dig_led settings " << digl << endl;

      ok = digl->LoadSettings(vpar);

      if( !ok )
        cerr << " Bad SADC settings for the cell " << ic << " Calorimeter " << GetName() << " take actions pls. ! " << endl;
    }
  }

  if( sadc_shape_filter_slc_use_ )
  {
// Configure usage of SLC ( Severin Lois Charpignon Filter )
    for( unsigned ic = 0; ic < NCells(); ic++ )
    {
      CsDigitizerSADC * dig =  dynamic_cast<CsDigitizerSADC *>(GetDigitizersSADC()[ic]);
      if( dig == NULL ) continue;
      if( calorimeter_sub_sets[index_olga_jura_].IsMyCell(ic) )
      {
      }
      else if( calorimeter_sub_sets[index_olga_saleve_].IsMyCell(ic) )
      {
      }
      else if( calorimeter_sub_sets[index_maintz_bottom_].IsMyCell(ic) )
      {
      }
      else if( calorimeter_sub_sets[index_maintz_top_].IsMyCell(ic) )
      {
      }
      else
      {
        // Set up shape_filter flag
        dig->use_shape_filter_slc_=true;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

bool CsECAL1::Reconstruction(void)
{
// ?? SetEventIDInfo() Must be done at Coral Level for any CsCalorimeter ... TODO check
// get event number, trigger info, ... for processing
  SetEventIDInfo();
  return CsCalorimeter::Reconstruction();
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1::ScaleCalib2006( size_t when )
{
  cout << " WARNING!!!!! CsECAL1::ScaleCalib2006 was called CALIBRATIONS HAS CHANGED!!!!  " << endl;
  double gams_size = 38.;
  double maintz_size = 75.;
  double olga_size = 143.;
  int nproblems = 0;
  for( size_t it=0; it!=NCells(); it++ )
  {
    double sx = cells[it].GetCellType().GetSizeX();

    double entries = cells_info[CALIB][when][it].GetEntries();
    double fac=cells_info[Reco::Calorimeter::CALIB][when][it].GetMean();
    double sfac=cells_info[Reco::Calorimeter::CALIB][when][it].GetSigma();
    if( fabs( sx - gams_size ) < 5. )
    {
      cells_info[CALIB][when][it]=StatInfo(entries,fac/1.1,sfac);
    }
    else if(  fabs( sx - maintz_size ) < 5. )
    {
      cells_info[CALIB][when][it]=StatInfo(entries,fac/1.15,sfac);
    }
    else if(  fabs( sx - olga_size ) < 5. )
    {
      cells_info[CALIB][when][it]=StatInfo(entries,fac/1.1,sfac);
    }
    else
    {
      nproblems++;
    }
  }

  if( nproblems != 0 )
  {
    cerr << " ECAL1 internal problem " << endl;
    exit(1);
  }
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1::ScaleCalib2007( size_t when )
{
  cout << " WARNING!!!!! CsECAL1::ScaleCalib2007 was called CALIBRATIONS HAS CHANGED!!!!  " << endl;
  double gams_size = 38.;
  double maintz_size = 75.;
  double olga_size = 143.;
  int nproblems = 0;
  for( size_t it=0; it!=NCells(); it++ )
  {
    double sx = cells[it].GetCellType().GetSizeX();

    double entries = cells_info[CALIB][when][it].GetEntries();
    double fac=cells_info[Reco::Calorimeter::CALIB][when][it].GetMean();
    double sfac=cells_info[Reco::Calorimeter::CALIB][when][it].GetSigma();
    if( fabs( sx - gams_size ) < 5. )
    {
      cells_info[CALIB][when][it]=StatInfo(entries,fac*0.8447,sfac);
    }
    else if(  fabs( sx - maintz_size ) < 5. )
    {
      cells_info[CALIB][when][it]=StatInfo(entries,fac*0.6509,sfac);
    }
    else if(  fabs( sx - olga_size ) < 5. )
    {
      cells_info[CALIB][when][it]=StatInfo(entries,fac*0.7618,sfac);
    }
    else
    {
      nproblems++;
    }
  }

  if( nproblems != 0 )
  {
    cerr << " ECAL1 internal problem " << endl;
    exit(1);
  }
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1::SetDefaultSADCDigitizationParameters( void )
{
  sadc_format_version_=1;
  sadc_decode_version_=1;
  sadc_delta_ = 1;
  sadc_ped_min_ = 0;
  sadc_ped_max_ = 9;
  sadc_signal_min_ = 10;
  sadc_signal_max_ = 30;
  sadc_front_ = 13;
  sadc_front_gate_ = 3;
  coeff_convert2max_ = 0.12;
  sadc_clock_ = 12.86;
  sadc_max_position_  = sadc_front_ + sadc_front_gate_;
  make_tcs_corrections_ = true;
}

////////////////////////////////////////////////////////////////////////////////

int CsECAL1::GetGAMSCellOfColumnRow( int x, int y) const
{
  return GetSubSets()[index_gams_].GetCellOfColumnRow(x,y);
}

////////////////////////////////////////////////////////////////////////////////

int CsECAL1::GetMaintzCellOfColumnRow( int x, int y) const
{
  int nrows = GetSubSets()[index_maintz_top_].GetNRows();
  if( y >= nrows )
  {
    return GetSubSets()[index_maintz_top_].GetCellOfColumnRow(x,y-nrows);
  }
  else
  {
    return GetSubSets()[index_maintz_bottom_].GetCellOfColumnRow(x,y);
  }
}

////////////////////////////////////////////////////////////////////////////////

int CsECAL1::GetOlgaCellOfColumnRow( int x, int y) const
{
  bool debug = false;
  if(debug ) cout << "CsECAL1::GetOlgaCellOfColumnRow x " << x << " y " << y << endl;
  int ncolomns = GetSubSets()[index_olga_jura_].GetNColumns();
  if( x >= ncolomns )
  {
    return GetSubSets()[index_olga_jura_].GetCellOfColumnRow(x-ncolomns,y);
  }
  else
  {
    return GetSubSets()[index_olga_saleve_].GetCellOfColumnRow(x,y);
  }
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1::SetDaqDataDecodingInfo( const CS::Chip::Digit &d  )
{
  bool debug = false;
  const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit*>( &d );
  if(ds != NULL )
  {
    std::string gname("EC01P00");
    std::string mname("EC01P01");
    std::string oname("EC01P02");
    std::string femname("EC01FEM");
    if( gname == ds->GetDetID().GetName() )
    {
      int32 x_digit = ds->GetX();
      int32 y_digit = ds->GetY();
      int icell = GetGAMSCellOfColumnRow(x_digit,y_digit);
      if( debug )
      {
        cout <<" Yees !! it is defenitely ECAL1 GAMS DIGIT !! For cell " << icell << endl;
      }
      if( icell >= 0 && icell < (int)NCells() )
          read_out_fe_amp_[icell]=CsCalorimeter::SADC_FEA;
      else
        cerr << " ERROR in CsECAL1::SetDaqDataDecodingInfo digit not identified for " << ds->GetDetID().GetName() <<
                                                        " invalid cell " << icell << " x = " << x_digit << " y = " << y_digit << endl;
    }
    else if( mname == ds->GetDetID().GetName() )
    {
      int32 x_digit = ds->GetX();
      int32 y_digit = ds->GetY();
      int icell = GetMaintzCellOfColumnRow(x_digit,y_digit);
      if( debug )
      {
        cout <<" Yees !! it is defenitely ECAL1 MAINZ DIGIT !! For cell " << icell << endl;
      }
      if( icell >= 0 && icell < (int)NCells() )
        read_out_fe_amp_[icell]=CsCalorimeter::SADC_FEA;
      else
        cerr << " ERROR in CsECAL1::SetDaqDataDecodingInfo digit not identified for " << ds->GetDetID().GetName() <<
                                                        " invalid cell " << icell << " x = " << x_digit << " y = " << y_digit << endl;
    }
    else if( oname == ds->GetDetID().GetName() )
    {
      int32 x_digit = ds->GetX();
      int32 y_digit = ds->GetY();
      int icell = GetOlgaCellOfColumnRow(x_digit,y_digit);
      if( debug )
      {
        cout <<" Yees !! it is defenitely ECAL1 OLGA DIGIT !! For cell " << icell << endl;
      }
      if( icell >= 0 && icell < (int)NCells() )
        read_out_fe_amp_[icell]=CsCalorimeter::SADC_FEA;
      else
        cerr << " ERROR in CsECAL1::SetDaqDataDecodingInfo digit not identified for " << ds->GetDetID().GetName() <<
                                                        " invalid cell " << icell << " x = " << x_digit << " y = " << y_digit << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1::DecodeChipDigits(const CS::Chip::Digits &digits)
{
  bool debug = false;
  if( debug ) cout << " debug CsECAL1::DecodeChipDigits skip_decoding_ = " << skip_decoding_ << endl;
  if( skip_decoding_ ) return;
  if( debug ) cout << " debug CsECAL1::DecodeChipDigits " << GetTBName() << " total digits size " << digits.size() << endl;

  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type

  std::string gname("EC01P00");
  std::pair<m_it,m_it> m_range = digits.equal_range( gname );
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " decode GAMS digit " << endl;
    const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
    int32 x_digit = ds->GetX();
    int32 y_digit = ds->GetY();
    int icell = GetGAMSCellOfColumnRow(x_digit,y_digit);
    if( icell < 0 )
    {
      cerr << " ERROR CsECAL1::DecodeChipDigits in GAMS wrong SADC digit x=" << x_digit << " y=" << y_digit << endl;
      ds->Print();
      exit(1);
    }
    if( debug ) cout << " Store GAMS digit x=" << x_digit << " y=" << y_digit << endl;
    StoreSADCDigit( icell, *ds , false );
  }

  std::string mname("EC01P01");
  m_range = digits.equal_range( mname );
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " decode Maintz digit " << endl;
    const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
    int32 x_digit = ds->GetX();
    int32 y_digit = ds->GetY();
    int icell = GetMaintzCellOfColumnRow(x_digit,y_digit);
    if( icell < 0 )
    {
      cerr << " ERROR CsECAL1::DecodeChipDigits in Maintz wrong SADC digit x=" << x_digit << " y=" << y_digit << endl;
      ds->Print();
      exit(1);
    }
    if( debug ) cout << " Store GAMS digit x=" << x_digit << " y=" << y_digit << endl;
    StoreSADCDigit( icell, *ds , false );
  }

  std::string oname("EC01P02");
  m_range = digits.equal_range( oname );
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " decode Olga digit " << endl;
    const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
    int32 x_digit = ds->GetX();
    int32 y_digit = ds->GetY();
    int icell = GetOlgaCellOfColumnRow(x_digit,y_digit);
    if( icell < 0 )
    {
      cerr << " ERROR CsECAL1::DecodeChipDigits in Olga wrong SADC digit x=" << x_digit << " y=" << y_digit << endl;
      ds->Print();
      exit(1);
    }
    if( debug ) cout << " Store Olga digit x=" << x_digit << " y=" << y_digit << endl;
    StoreSADCDigit( icell, *ds , false );
  }

  MakeCsDigits();

  if( make_sadc_histo_ ) DecodingTestSADC();
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1::DecodeChipDigitsLEDEvent(const CS::Chip::Digits &digits)
{
  bool debug = false;
  bool debug_new_fem = false;
  if( debug ) cout << " debug CsECAL1::DecodeChipDigits skip_led_decoding_ = " << skip_led_decoding_ << endl;
  if( skip_led_decoding_ ) return;
  if( debug ) cout << " debug CsECAL1::DecodeChipDigits " << GetTBName() << " total digits size " << digits.size() << endl;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  std::string gname("EC01P00");
  std::pair<m_it,m_it> m_range = digits.equal_range( gname );
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " decode GAMS digit " << endl;
    const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
    int32 x_digit = ds->GetX();
    int32 y_digit = ds->GetY();
    int icell = GetGAMSCellOfColumnRow(x_digit,y_digit);
    if( icell < 0 )
    {
      cerr << " ERROR CsECAL1::DecodeChipDigits in GAMS wrong SADC digit x=" << x_digit << " y=" << y_digit << endl;
      ds->Print();
      exit(1);
    }
    if( debug ) cout << " Store GAMS digit x=" << x_digit << " y=" << y_digit << endl;
    StoreSADCDigit( icell, *ds , true);
  }

  std::string mname("EC01P01");
  m_range = digits.equal_range( mname );
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " decode Maintz digit " << endl;
    const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
    int32 x_digit = ds->GetX();
    int32 y_digit = ds->GetY();
    int icell = GetMaintzCellOfColumnRow(x_digit,y_digit);
    if( icell < 0 )
    {
      cerr << " ERROR CsECAL1::DecodeChipDigits in Maintz wrong SADC digit x=" << x_digit << " y=" << y_digit << endl;
      ds->Print();
      exit(1);
    }
    if( debug ) cout << " Store GAMS digit x=" << x_digit << " y=" << y_digit << endl;
    StoreSADCDigit( icell, *ds , true);
  }

  std::string oname("EC01P02");
  m_range = digits.equal_range( oname );
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " decode Olga digit " << endl;
    const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
    int32 x_digit = ds->GetX();
    int32 y_digit = ds->GetY();
    int icell = GetOlgaCellOfColumnRow(x_digit,y_digit);
    if( icell < 0 )
    {
      cerr << " ERROR CsECAL1::DecodeChipDigits in Olga wrong SADC digit x=" << x_digit << " y=" << y_digit << endl;
      ds->Print();
      exit(1);
    }
    if( debug ) cout << " Store Olga digit x=" << x_digit << " y=" << y_digit << endl;
    StoreSADCDigit( icell, *ds , true);
  }

  MakeCsDigits();
  StoreRawDataStatistic();
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1::CreateGUI()
{
  if( gui!=NULL )
  {
    CsErrLog::msg(elFatal, __FILE__, __LINE__,
                  "CsECAL1::CreateGUI: attempted to create a second GUI for %s.",
                  getName());
    return;
  }

 gui = new CsECAL1GUI(*this,NULL,GetName().c_str());
}

////////////////////////////////////////////////////////////////////////////////

string CsECAL1::GetCellName            (int icell) const
{
  char cell_name[132];
  double sx = cells[icell].GetCellType().GetSizeX();
  if( fabs( sx - cell_sizex_gams_) < 5. )
  {
    int x = GetSubSets()[index_gams_].GetColumnOfCell( icell );
    int y = GetSubSets()[index_gams_].GetRowOfCell( icell );
    if( x < 0 || y < 0 )
    {
       cerr << " CsECAL1::GetCellName internal error cell GAMS ? " << icell << " x " << x << " y " << endl;
       exit(1);
    }
    sprintf(cell_name,"GAMS_%d_%d",x,y);
  }
  else if( fabs( sx - cell_sizex_maintz_) < 5. )
  {
    if( cells[icell].GetY() > 0 )
    {
      int x = GetSubSets()[index_maintz_top_].GetColumnOfCell( icell );
      int y = GetSubSets()[index_maintz_top_].GetRowOfCell( icell );
      if( x < 0 || y < 0 )
      {
        cerr << " CsECAL1::GetCellName internal error cell MaintzTop ? " << icell << " x " << x << " y " << endl;
        exit(1);
      }
      sprintf(cell_name,"MaintzTop_%d_%d",x,y);
    }
    else
    {
      int x = GetSubSets()[index_maintz_bottom_].GetColumnOfCell( icell );
      int y = GetSubSets()[index_maintz_bottom_].GetRowOfCell( icell );
      if( x < 0 || y < 0 )
      {
        cerr << " CsECAL1::GetCellName internal error cell MaintzBottom ? " << icell << " x " << x << " y " << endl;
        exit(1);
      }
      sprintf(cell_name,"MaintzBottom_%d_%d",x, y);
    }
  }
  else
  {
    if( cells[icell].GetX() > 0 )
    {
      int x = GetSubSets()[index_olga_jura_].GetColumnOfCell( icell );
      int y = GetSubSets()[index_olga_jura_].GetRowOfCell( icell );
      if( x < 0 || y < 0 )
      {
        cerr << " CsECAL1::GetCellName internal error cell OlgaJura ? " << icell << " x " << x << " y " << endl;
        exit(1);
      }
      sprintf(cell_name,"OlgaJura_%d_%d",x, y);
    }
    else
    {
      int x = GetSubSets()[index_olga_saleve_].GetColumnOfCell( icell );
      int y = GetSubSets()[index_olga_saleve_].GetRowOfCell( icell );
      if( x < 0 || y < 0 )
      {
        cerr << " CsECAL1::GetCellName internal error cell OlgaSaleve ? " << icell << " x " << x << " y " << endl;
        exit(1);
      }
      sprintf(cell_name,"OlgaSaleve_%d_%d",x, y);
    }
  }

  return cell_name;

}

////////////////////////////////////////////////////////////////////////////////

int CsECAL1::OutputCalibInfo( string &s ) const
{
  bool debug = false;

  if( debug ) cout << " Calorimeter::OutputCalibInfo for " <<
                   (this)->GetName() << " this = " << this << endl;

  char o[500000];
  o[0]=0;
  if( debug ) cout << " At start strlen(o) = " << strlen(o) << endl;

  sprintf(o,"### \n");
  sprintf(o+strlen(o),"Calibration Calorimeter:  %s  Ncells=  %zu \n",GetName().c_str(),NCells());
  sprintf(o+strlen(o),"### \n");

  for( size_t isb=0; isb < GetSubSets().size(); isb++ )
  {
  for( size_t indx=0; indx < GetSubSets()[isb].GetCells().size(); indx++ )
  {
     size_t icell = GetSubSets()[isb].GetCells()[indx];
     int x = GetSubSets()[isb].GetColumnOfCell(icell);
     int y = GetSubSets()[isb].GetRowOfCell(icell);
     if( isb == index_gams_ )
     {
       sprintf(o+strlen(o)," G %d %d ",x,y);
     }
     else if( isb == index_maintz_top_ )
     {
       sprintf(o+strlen(o)," M %d %d ",x,y+GetSubSets()[index_maintz_top_].GetNRows());
     }
     else if( isb == index_maintz_bottom_ )
     {
       sprintf(o+strlen(o)," M %d %d ",x,y);
     }
     else if( isb == index_olga_saleve_ )
     {
       sprintf(o+strlen(o)," O %d %d ",x,y);
     }
     else if( isb == index_olga_jura_ )
     {
       sprintf(o+strlen(o)," O %2d %2d ",x+GetSubSets()[index_olga_jura_].GetNColumns(),y);
     }
     else
     {
       cerr << " ECAL1 internal error " << endl;
       exit(1);
     }

     double entries = cells_info[CALIB][NEW][icell].GetEntries();
     double fac=cells_info[Reco::Calorimeter::CALIB][NEW][icell].GetMean();
     double sfac=cells_info[Reco::Calorimeter::CALIB][NEW][icell].GetSigma();

     bool tolerance_check = true;
     if(fac>0.) tolerance_check = (fabs( 1.-1./fac ) > options.tolerate2keepold_calib );

     if( fac <= 0 )
     {
       if( debug ) cout << " ERROR OutputCalibInfo " << GetName() << " Ignore Bad calibration " << fac <<
                                      " in cell " << GetCellName(icell) << " use old one! " << endl;
     }

     if( entries > options.new_calib_stat_min && fac > 0 && tolerance_check )
     {
       if( debug )
       {
         cout << " Stat OK in cell " << GetCellName(icell) <<
                 " cells_info[CALIB][NEW][icell].GetMean() " << cells_info[CALIB][NEW][icell].GetMean() <<
                 " cells_info[CALIB][OLD][icell].GetMean() " << cells_info[CALIB][OLD][icell].GetMean() << endl;

       }

       double cfnew= cells_info[Reco::Calorimeter::CALIB][OLD][icell].GetMean()/fac;
       double scfnew= sfac/fac * cfnew;

       if( debug )
       {
         printf(" %8.2f %8.2f %8.2f %s %zu \n",
                    1000.*cfnew,1000.*scfnew,cells_info[CALIB][NEW][icell].GetEntries(), GetCellName(icell).c_str(), icell );
       }
       sprintf(o+strlen(o)," %8.2f %8.2f %8.2f %s %zu \n",
                    1000.*cfnew,1000.*scfnew,cells_info[CALIB][NEW][icell].GetEntries(), GetCellName(icell).c_str(), icell );
     }
     else
     {

       if( debug )
       {
         cout << " save OLD " << endl;
         printf(" %8.2f %8.2f %8.2f %s %zu \n",
          1000.*cells_info[CALIB][OLD][icell].GetMean(),
          1000.*cells_info[CALIB][OLD][icell].GetSigma(),
            1., GetCellName(icell).c_str(), icell);
       }
        sprintf(o+strlen(o)," %8.2f %8.2f %8.2f %s %zu \n",
          1000.*cells_info[CALIB][OLD][icell].GetMean(),
          1000.*cells_info[CALIB][OLD][icell].GetSigma(),
            1., GetCellName(icell).c_str(), icell);
    }
  }
  }
  if( debug ) cout << " At finish strlen(o) = " << strlen(o) << endl;
  string ss(o);
  s = ss;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int CsECAL1::InputCalibInfo(size_t when, const string &s )
{
//   bool debug = true;
  bool debug = false;
  bool error = 0;
  istringstream is(s.c_str());
  char calorim_name[132],dummy[132];
  string str;
  if( !(when == OLD || when == PRIM || when == MC) )
  {
    cerr << " CsECAL1::InputCalibInfo supposed to be OLD, PRIM or MC " << when << endl;
    return 2;
  }


  //  Read line of comments
    getline(is,str);
    getline(is,str);
    if( debug ) cout << str << endl;
    sscanf(str.c_str(),"%s %s ",calorim_name,dummy);
    getline(is,str);

    if( debug )
    {
      cout << " Calorimeter " << calorim_name <<" with " << dummy  << endl;
    }

    int not_matched_cells_cnt = 0;
  //  Read calibration
    while( getline(is,str) )
    {
      if( debug ) cout << str << endl;
      char subset_name[132];
      char cell_name[132];
      int icell,x,y;
      float w,c,sc;
      if( !sscanf(str.c_str()," %s %d %d %g %g %g %s %d \n",subset_name,&x,&y,&c,&sc,&w,cell_name,&icell) )
      {
        cerr << " FRORMAT ERROR for InputCalibInfo in line " << str << endl;
        error = 1;
        break;
      }
      int jcell = FindCellByName( cell_name );
      if(debug) cout << " Cell name " << cell_name << " jcell " << jcell << endl;
      if( !(jcell >=0 && jcell < (int)NCells()) )
      {
        cerr << " Unexpected Cell name " << cell_name << " in CsECAL1::InputCalibInfo " <<
                                                               GetName() << "  " << when  << endl;
        cerr << " Input strig: " << endl;
        cerr << str << endl;
        cerr << " It well might be that you use wrong calibration file for this detector " << GetName() << endl;
        stringstream ss;
        ss << " Unexpected Cell name " << GetName() << " CsECAL1::InputCalibInfo " << str;
        throw Exception(ss.str().c_str());
      }

      if( icell != jcell )
      {
        icell = jcell;
        not_matched_cells_cnt ++;
      }
      if( debug )
      {
        cout << subset_name << " x " << x << " y " << y << " c " << c/1000. << " sc " << sc/1000. << " stat " << w << " icell " << icell << endl;
      }

      if(icell >=0 && icell < (int)NCells() )
        SetCellInfo(CALIB, when,icell,Reco::StatInfo(w,c/1000.,sc/1000.));
    }

    if( not_matched_cells_cnt > 0 )
    {
      cerr << " WARNING!!! CsECAL1::InputCalibInfo " << GetName() <<
               " Not matching in cells id was detected " << not_matched_cells_cnt << " times " << endl;
      cerr << " You use wrong calibaration file or calibrations were produced with different geometry descriptor !!! " << endl;
    }

  return error;

}

////////////////////////////////////////////////////////////////////////////////
#include "CsCalorimeterHist.h"

void CsECAL1::EndOfJob( void )
{
  Reco::Calorimeter::EndOfJob();
  cout << " CsECAL1 now we need to write SADC profiles, wait a bit .. " << endl;
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

int CsECAL1::InputJOUJOUInfo(const string &s)
{
  bool debug = false;
  istringstream is(s.c_str());
  string str;
  if( debug ) cout <<  " InputJOUJOUInfo " << GetName() << endl;

  vector < pair<int,int> > gamsxy;
  for( int ix=0; ix < 44; ix++ )
  {
    for( int iy=0; iy < 24; iy++ )
    {
      if( iy > 3 && iy<20 && ix > 7 && ix<36 ) continue;

      gamsxy.push_back( pair<int,int>(ix,iy) );
    }
  }
  if( debug ) cout <<" cells in GAMS part " << gamsxy.size() << endl;

//  Read calibration
  float p0, p1, p2;
  int i0;
  float v0, v1, v2, v3, v4, v5, v6;
  getline(is,str);
  sscanf(str.c_str()," %g %g %g %d %g %g %g %g %g %g %g",&p0,&p1,&p2,&i0,&v0,&v1,&v2,&v3,&v4,&v5,&v6);
  if( debug ) cout << " Electron beam energy ?? " << p2 << endl;
  bool gams = true;
  bool mainz = false;
  bool olga = false;
  while( getline(is,str) )
  {
    int ncell, stat0, stat1;
    float pedpos, elpos, ledpos;
    sscanf(str.c_str()," %d %g %g %g %d %d",&ncell,&pedpos,&elpos,&ledpos,&stat0,&stat1);
    if( debug ) cout << " " << str << " " << endl;
    int ix_gams = -1;
    int iy_gams = -1;
    int ix_mainz = -1;
    int iy_mainz = -1;
    int ix_olga = -1;
    int iy_olga = -1;
    int icell = -1;
    if( gams )
    {
      if( ncell >= (int)gamsxy.size() )
      {
         cerr << " Internal ERROR ncell = " <<  ncell <<" But GAMS size = " << gamsxy.size() << endl;
      }
      ix_gams = gamsxy[ncell].first;
      iy_gams = gamsxy[ncell].second;

      icell = GetSubSets()[index_gams_].GetCellOfColumnRow( ix_gams, iy_gams );
    }
    else if( mainz)
    {
      ix_mainz = (ncell - (int)gamsxy.size() )/26;
      iy_mainz = (ncell - (int)gamsxy.size() ) - ix_mainz*26;
      if( iy_mainz > 12 )
        icell = GetSubSets()[index_maintz_top_].GetCellOfColumnRow( ix_mainz, iy_mainz -13 );
      else
        icell = GetSubSets()[index_maintz_bottom_].GetCellOfColumnRow( ix_mainz, iy_mainz );

    }
    else if( olga)
    {
      ix_olga = (ncell- (int)gamsxy.size() - 26*22 )/20;
      iy_olga = (ncell- (int)gamsxy.size() - 26*22 ) - ix_olga*20;
      if( ix_olga > 7 )
        icell = GetSubSets()[index_olga_jura_].GetCellOfColumnRow( ix_olga-8, iy_olga );
      else
        icell = GetSubSets()[index_olga_saleve_].GetCellOfColumnRow( ix_olga, iy_olga );
    }

    if( debug )
    {
      if( gams )
      {
        cout << " ncell " << ncell << " icell " << icell << " ix_gams " << ix_gams << " iy_gams " << iy_gams << endl;
      }
      else if( mainz)
      {
        cout << " ncell " << ncell << " icell " << icell << " ix_mainz " << ix_mainz << " iy_mainz " << iy_mainz << endl;
      }
      else if( olga)
      {
        cout << " ncell " << ncell << " icell " << icell  << " ix_olga " << ix_olga << " iy_olga " << iy_olga << endl;
      }
      else
      {
        cout << " ncell " << ncell << " icell " << icell  << " NONE  Should never be here" << endl;
        exit(1);
      }
    }

    if(icell < 0 )
    {
      cerr << " Internal ERROR No cell ncell = " << ncell << endl;
    }


    double w=1.;;
    double sc=0.001;
    double c=15./elpos;
    SetCellInfo(CALIB, MONITOR,icell,Reco::StatInfo(w,c,sc));

    c=ledpos;
    SetCellInfo(LED, MONITOR,icell,Reco::StatInfo(w,c,sc));

    if( ncell == (int)gamsxy.size() -1 )
    {
      cout <<" Last GAMS cell " << ncell << " switch to Mainz " << endl;
      gams = false;
      mainz = true;
      olga = false;
    }
    if( ncell == (int)gamsxy.size() -1+26*22 )
    {
      cout <<" Last Maintz cell " << ncell << " switch to Olga " << endl;
      gams = false;
      mainz = false;
      olga = true;
    }
  }
  if( debug )
  {
    cout <<  " InputJOUJOUInfo end to debug " << GetName() << endl;
    exit(0);
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

