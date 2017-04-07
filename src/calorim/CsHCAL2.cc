/*!
   \file    CsHCAL2.cc
   \brief   Implementation class for HCAL2 COMPASS calorimeter
   \version $Revision: 1.38 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \author  Denis Murashev
   \date    $Date: 2011/02/01 22:05:52 $
*/

#include <cmath>

#include "coral_config.h"

#include "CsEvent.h"
#include "DaqDataDecoding/DaqEvent.h"
#include "DaqDataDecoding/Scaler.h"
#include "CsHCAL2.h"
#include "CsEvent.h"

#include "Reco/CalorimeterHist.h"
#include "Reco/CellDataRaw.h"

using namespace std;
using namespace Reco;

////////////////////////////////////////////////////////////////////////////////

CsHCAL2::CsHCAL2(const string &name, const string &geom_file) :
 CsCalorimeter(name, geom_file),
 offset_2x2_trig(0),
 offset_layer1_trig(0),
 offset_layer2_trig(0),
 offset_layer3_trig(0),
 offset_layer4_trig(0)
{
}

////////////////////////////////////////////////////////////////////////////////

void  CsHCAL2::InitOptions( void )
{
  CsCalorimeter::InitOptions();

  options.calo_type                     = Hadronic;
  options.particle_default              = CalorimeterParticle::PION_PLUS;
  options.default_calibration           = 0.100;
  options.readout_sparsified_           = true;
  options.delta_spars_mode_             = 3.;

  // MC options
  options.mc_default_calibration     = 0.100;
  options.mc_make_real_digitization  = true;
  options.mc_make_fiadc_digitization = true;
  options.mc_fiadc_sparce_delta      = 3.;
}

////////////////////////////////////////////////////////////////////////////////

void  CsHCAL2::Initialize( void )
{
  CsCalorimeter::Initialize();

  InitXY();

// Not sure. # warning CsHCAL2:: Temporary trigger groups setting. Should be removed by CsCalorimeter options from det.dat  file.
// Set trigger groups in HCAL2
  offset_2x2_trig = 0;
  for( int yg=0; yg!=  GetNRows()/2; yg++ ) for( int xg=0; xg!= GetNColumns() /2; xg++ )
  {
    vector<size_t> g;
    int x = xg*2;
    int y = yg*2;
    int icell = GetCellOfColumnRow(x,y+1);
    if(icell >= 0 ) g.push_back(icell);
    icell = GetCellOfColumnRow(x+1,y+1);
    if(icell >= 0 ) g.push_back(icell);
    icell = GetCellOfColumnRow(x,y);
    if(icell >= 0 ) g.push_back(icell);
    icell = GetCellOfColumnRow(x+1,y);
    if(icell >= 0 ) g.push_back(icell);
    char tg_name[132];
    sprintf(tg_name," 2x2_TrGr_%d_%d ",xg,yg);
    TrigGroup tg( this, tg_name,g,3,xg,yg);
    if(g.size() > 0 ) trigger_groups.push_back(tg);
//    cout << " Trig group " << trigger_groups.size()-1 << " Amount of
//    elements " << trigger_groups.back().NElements() << " GroupName " <<
//    trigger_groups.back().GetName() << endl;
  }
// Trigger groups setting
  offset_layer1_trig = NGroups();
  for( int yg=0; yg!=2; yg++ ) for( int xg=0; xg!=5; xg++ )
  {
//    cout << " xg " << xg << " yg " << yg << endl;
    vector<size_t> g;
    int x = xg*4;
    int y = yg*4+2;
    for( int yc=0; yc!=4; yc++ ) for( int xc=0; xc!=4; xc++ )
    {
      int xcc = x+xc;
      int ycc = y+yc;
      int icell = GetCellOfColumnRow(xcc,ycc);
      if(icell >= 0 ) g.push_back(icell);
    }
    char tg_name[132];
    int ibox = 5*yg+xg;
    if(ibox > 4 ) {ibox -= 5;} else {ibox += 5;}
    ibox=ibox+1;
    if(ibox > 5 ) ibox=ibox+1;
    sprintf(tg_name," 4x4_Layer1_%d_%d Box%d",xg,yg,ibox);
    TrigGroup tg( this, tg_name,g,1,xg,yg);
    if(g.size() > 0 ) trigger_groups.push_back(tg);
//    cout << " Trig group " << trigger_groups.size()-1 << " N elements " << trigger_groups.back().NElements() <<
//                                                              " GroupName " << trigger_groups.back().GetName() << endl;
  }

  offset_layer2_trig = NGroups();
  for( int yg=0; yg!=2; yg++ ) for( int xg=0; xg!=5; xg++ )
  {
//    cout << " xg " << xg << " yg " << yg << endl;
    vector<size_t> g;
    int x = xg*4+2;
    int y = yg*4;
    for( int yc=0; yc!=4; yc++ ) for( int xc=0; xc!=4; xc++ )
    {
      int xcc = x+xc;
      int ycc = y+yc;
      int icell = GetCellOfColumnRow(xcc,ycc);
      if(icell >= 0 ) g.push_back(icell);
    }
    char tg_name[132];
    int ibox = 5*yg+xg;
    if(ibox > 4 ) {ibox -= 5;} else {ibox += 5;}
    ibox=ibox+1;
    sprintf(tg_name," 4x4_Layer2_%d_%d Box%d",xg,yg,ibox);
    TrigGroup tg( this, tg_name,g,2,xg,yg);
    if(g.size() > 0 ) trigger_groups.push_back(tg);
//    cout << " Trig group " << trigger_groups.size()-1 << " N elements " << trigger_groups.back().NElements() <<
//                                                              " GroupName " << trigger_groups.back().GetName() << endl;
  }

  offset_layer3_trig = NGroups();
  for( int yg=0; yg!=2; yg++ ) for( int xg=0; xg!=5; xg++ )
  {
//    cout << " xg " << xg << " yg " << yg << endl;
    vector<size_t> g;
    int x = xg*4+2;
    int y = yg*4+2;
    for( int yc=0; yc!=4; yc++ ) for( int xc=0; xc!=4; xc++ )
    {
      int xcc = x+xc;
      int ycc = y+yc;
      int icell = GetCellOfColumnRow(xcc,ycc);
      if(icell >= 0 ) g.push_back(icell);
    }
    char tg_name[132];
    int ibox = 5*yg+xg;
    if(ibox > 4 ) {ibox -= 5;} else {ibox += 5;}
    ibox=ibox+1;
    sprintf(tg_name," 4x4_Layer3_%d_%d Box%d",xg,yg,ibox);
    TrigGroup tg( this, tg_name,g,3,xg,yg);
    if(g.size() > 0 ) trigger_groups.push_back(tg);
//    cout << " Trig group " << trigger_groups.size()-1 << " N elements " << trigger_groups.back().NElements() <<
//                                                              " GroupName " << trigger_groups.back().GetName() << endl;
  }

  offset_layer4_trig = NGroups();
  for( int yg=0; yg!=2; yg++ ) for( int xg=0; xg!=5; xg++ )
  {
//    cout << " xg " << xg << " yg " << yg << endl;
    vector<size_t> g;
    int x = xg*4;
    int y = yg*4;
    for( int yc=0; yc!=4; yc++ ) for( int xc=0; xc!=4; xc++ )
    {
      int xcc = x+xc;
      int ycc = y+yc;
      int icell = GetCellOfColumnRow(xcc,ycc);
      if(icell >= 0 ) g.push_back(icell);
    }
    char tg_name[132];
    int ibox = 5*yg+xg;
    if(ibox > 4 ) {ibox -= 5;} else {ibox += 5;}
    ibox=ibox+1;
    sprintf(tg_name," 4x4_Layer4_%d_%d Box%d",xg,yg,ibox);
    TrigGroup tg( this, tg_name,g,4,xg,yg);
    if(g.size() > 0 ) trigger_groups.push_back(tg);
//    cout << " Trig group " << trigger_groups.size()-1 << " N elements " << trigger_groups.back().NElements() <<
//                                                              " GroupName " << trigger_groups.back().GetName() << endl;
  }

  Reco::StatInfo zero = Reco::StatInfo(0.,0.,0.);
  double time_group_default_calibration = -295;
  double energy_group_default_calibration = 0.15;
  for(size_t it=0; it!=NGroups(); it++)
  {
    group_info[TrigGroup::CALIB][OLD].push_back(Reco::StatInfo(1,energy_group_default_calibration,1.));
    group_info[TrigGroup::CALIB][NEW].push_back(zero);
    group_info[TrigGroup::TIME1][OLD].push_back(Reco::StatInfo(1,time_group_default_calibration,1.));
    group_info[TrigGroup::TIME1][NEW].push_back(zero);
    group_info[TrigGroup::TIME2][OLD].push_back(Reco::StatInfo(1,time_group_default_calibration,1.));
    group_info[TrigGroup::TIME2][NEW].push_back(zero);

    if( it> 0 ) group_info[TrigGroup::CALIB][OLD].back() = Reco::StatInfo(1,energy_group_default_calibration/4.,1.);
  }

  double sizez0 = cells[0].GetCellType().GetSizeZ();
  for( size_t it=0; it!=NCells(); it++ )
  {
    double sizez = cells[it].GetCellType().GetSizeZ();
    individual_calib_factor_[it] = sizez/sizez0;
//     cout << it << " zc " << cells[it].GetZ() << " factor " << individual_calib_factor_[it] << endl;
  }

//  if( options.print_general_info ) PrintGeneralInfo();

}

////////////////////////////////////////////////////////////////////////////////

void CsHCAL2::DecodeChipDigit(const CS::Chip::Digit &digit)
{
//   bool debug_sadc = true;
  bool debug_sadc = false;
  const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  const CS::ChipADC::Digit *d = dynamic_cast<const CS::ChipADC::Digit *>(&digit);
  const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit *>(&digit);
  if( (dt==NULL) && (d==NULL) && (ds==NULL) )
  {
    cout << "CsHCAL2::DecodeChipDigit(): Wrong digit type"  << endl;
    throw CS::Exception("CsCalorimeter::DecodeChipDigit(): Wrong digit type");
  }
  if( d!=NULL )
  {
    int32 x_digit = d->GetX();
    int32 y_digit = d->GetY();
    int32 e_digit = d->GetAmplitude();
#if calorim_debug>1
    printf(" DecodeChipDigit: Calorimeter %s x= %d y= %d E=%d \n",getName(),x_digit,y_digit,e_digit);
#endif
    int icell = GetCellOfColumnRow(x_digit,y_digit);
    if( icell >= 0 )
    {
//  Amplitudes from calorimeter cells
      signals_.push_back( CellDataRaw(icell) );
      signals_.back().SetAmplitude(e_digit);
      signals_.back().SetEnergy(e_digit * cells_info[CALIB][OLD][icell].GetMean());
//  Summ ADC energy in groups
      for( size_t it=0; it!=NGroups(); it++ )
      {
        double e=signals_.back().GetEnergy();
        int icell = signals_.back().GetCellIdx();
        for( vector<size_t>::const_iterator it2=trigger_groups[it].GetCells().begin();
                                            it2!=trigger_groups[it].GetCells().end(); it2++ )
          if( icell == (int)*(it2) ) trigger_groups[it].trig_group_data.AddADC(e);
      }

    }
    else
    {
//  Check if it is trigger sum amplitudes
// old maps       if( y_digit == 0 && x_digit >= 100 && x_digit <= 155 )
       if( y_digit >= 300 && y_digit <= 304 && x_digit >= 0 )
       {
//   This is trigger sum of 2x2 cells matrixes
//          printf("DecodeChipDigit:  2x2 TrigGroup = %d ESum=%d OF CALORIMETER %s \n",
//                                                                        x_digit-100,e_digit,getName());
//  Check are this trigger groups initialised ?
          if(offset_2x2_trig >=0 )
          {
            int itrig =  x_digit+11*(y_digit-300)+offset_2x2_trig;
            if( itrig >= (int)trigger_groups.size() )
            {
               printf(" Error:: DecodeChipDigit: Calorimeter %s TriggerGroup %d is OUT of RANGE %zu \n",getName(),itrig,trigger_groups.size());
            }
            else
            {
               trigger_groups[itrig].trig_group_data.SetAmpl( e_digit*group_info[TrigGroup::CALIB][OLD][itrig].GetMean() );
            }
          }
       }
//       else if( y_digit == 0 && x_digit >= 200 && x_digit <= 210 )
       else if( y_digit >= 100 && y_digit <= 101 && x_digit >= 0)
       {
//   This is trigger sum of 4x4 cells matrixes
//          printf("DecodeChipDigit:  4x4 TrigGroupLayer1 BOX = %d ESum=%d OF CALORIMETER %s \n",
//                                                                        x_digit-200,e_digit,getName());
//  Check are this trigger groups initialised ?
          if(offset_layer1_trig >=0 )
          {
            int itrig =  x_digit+5*(y_digit-100)+offset_layer1_trig;
            if( itrig >= (int)trigger_groups.size() )
            {
               printf(" Error:: DecodeChipDigit: Calorimeter %s TriggerGroup %d is OUT of RANGE %zu \n",getName(),itrig,trigger_groups.size());
            }
            else
            {
               trigger_groups[itrig].trig_group_data.SetAmpl( e_digit*group_info[TrigGroup::CALIB][OLD][itrig].GetMean() );
            }
          }
       }
//        else if( y_digit == 0 && x_digit >= 300 && x_digit <= 309 )
       else if( y_digit >= 200 && y_digit <= 201 && x_digit >= 0 )
       {
//   This is trigger sum of 4x4 cells matrixes
//          printf("DecodeChipDigit:  4x4 TrigGroupLayer2 BOX = %d ESum=%d OF CALORIMETER %s \n",
//                                                                        x_digit-300,e_digit,getName());
//  Check are this trigger groups initialised ?
          if(offset_layer2_trig >=0 )
          {
            int itrig =  x_digit+5*(y_digit-200)+offset_layer2_trig;
            if( itrig >= (int)trigger_groups.size() )
            {
               printf(" Error:: DecodeChipDigit: Calorimeter %s TriggerGroup %d is OUT of RANGE %zu \n",getName(),itrig,trigger_groups.size());
            }
            else
            {
               trigger_groups[itrig].trig_group_data.SetAmpl( e_digit*group_info[TrigGroup::CALIB][OLD][itrig].GetMean() );
            }
          }
       }
       else
       {
       #if calorim_debug>1
         if( icell < 0 ) printf("DecodeChipDigit: digit x= %d y= %d E=%d OUT OF CALORIMETER %s \n",
                                                                        x_digit,y_digit,e_digit,getName());
       #endif
       }
    }
//     cout << "CsHCAL2::DecodeChipDigit(): data size "  << raw_data_store.size() << endl;
  }
  else if( dt!=NULL )
  {
    //    const float F1bin = 0.1289231; // ns
    int32 ch_pos_digit = dt->GetChannelPos();
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();

  double time=dt->GetTimeDecoded();
       #if calorim_debug>1
      printf("CsHCAL2::DecodeChipDigit: Calorimeter %s Ch %d X %d Y %d ChPos %d Time %d \n",
                                                getName(),ch_digit,x_digit,y_digit,ch_pos_digit,t_digit);
      #endif
//        printf("CsHCAL2::DecodeChipDigit: Calorimeter %s Ch %d ChPos %d Time %d \n",
//                                                                  getName(),ch_digit,ch_pos_digit,t_digit);
    if( ch_pos_digit == 1 && y_digit >= 100 && y_digit <= 101 )
    {
      if(offset_layer1_trig >=0 )
      {
         int itrig =  x_digit+5*(y_digit-100)+offset_layer1_trig;
         trigger_groups[itrig].trig_group_data.SetTime1( time - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
//        printf("CsHCAL2::DecodeChipDigit: Calorimeter %s layer 1 Ch %d  Time %d \n",
//                                                                  getName(),itrig-offset_layer1_trig,t_digit);
      }
    }
    if( ch_pos_digit == 2 && y_digit >= 100 && y_digit <= 101 )
    {
      if(offset_layer1_trig >=0 )
      {
         int itrig =  x_digit+5*(y_digit-100)+offset_layer1_trig;
         trigger_groups[itrig].trig_group_data.SetTime2( time - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
//        printf("CsHCAL2::DecodeChipDigit: Calorimeter %s layer 1 Ch %d  Time %d \n",
//                                                                  getName(),itrig-offset_layer1_trig,t_digit);
      }
    }
    if( ch_pos_digit == 1 && y_digit >= 200 && y_digit <= 201 )
    {
      if(offset_layer2_trig >=0 )
      {
         int itrig =  x_digit+5*(y_digit-200)+offset_layer2_trig;
         trigger_groups[itrig].trig_group_data.SetTime1( time - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
//        printf("CsHCAL2::DecodeChipDigit: Calorimeter %s layer 2 Ch %d  Time %d \n",
//                                                                  getName(),itrig-offset_layer2_trig,t_digit);
      }
    }
    if( ch_pos_digit == 2 && y_digit >= 200 && y_digit <= 201 )
    {
      if(offset_layer2_trig >=0 )
      {
         int itrig =  x_digit+5*(y_digit-200)+offset_layer2_trig;
         trigger_groups[itrig].trig_group_data.SetTime2( time - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
//        printf("CsHCAL2::DecodeChipDigit: Calorimeter %s layer 2 Ch %d  Time %d \n",
//                                                                  getName(),itrig-offset_layer2_trig,t_digit);
      }
    }
//     if(ch_digit >=0 && ch_digit <= (int)NGroups() )
//     {
//         int itrig =  ch_digit;
//         trigger_groups[itrig].SetTime( t_digit );
//     }


  }
  else if( ds!=NULL ) // ChipSADC
  {
    if( debug_sadc ) cout << " SADC digit in " << GetName() << endl;
    DecodeChipSADCDigit(*ds);
  }

}

////////////////////////////////////////////////////////////////////////////////

void CsHCAL2::DecodeTriggerGroups( void )
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " CsHCAL2::DecodeTriggerGroups " << GetName() << " debug " << endl;
  if( debug ) cout << " CsHCAL2::DecodeTriggerGroups " << GetName() << " debug " << endl;
  if( debug ) cout << " CsHCAL2::DecodeTriggerGroups " << GetName() << " debug " << endl;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  if( debug ) cout << " std::multimap declaration OK " << endl;
  std::string tgname("HC02P1T1");
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  std::pair<m_it,m_it> m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  if( debug ) cout << " std::CsEvent::Instance()->getChipDigits() .equal_range( " << tgname << " OK " << endl;

  // loop on all found digits
  int cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( dt == NULL ) continue;
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer1_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer1_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME1][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime1( time - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P1T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( dt == NULL ) continue;
     if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer1_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer1_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME2][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime2( time - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P2T1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P2T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P3T1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P3T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P4T1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( dt == NULL ) continue;
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer2_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer2_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME1][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime1( time - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P4T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( dt == NULL ) continue;
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer2_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer2_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME2][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime2( time - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

//   m_range = CsEvent::Instance()->getChipDigits().equal_range(GetName());
//   cnt = 0;
//   for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
//     cout << " HCAL2 for "<< GetName() <<" Didit #" << cnt++ << " found " << endl;

  if( debug ) cout << " CsHCAL2::DecodeTriggerGroups Finished OK " << endl;

}

////////////////////////////////////////////////////////////////////////////////

void CsHCAL2::DecodeTriggerGroups(const CS::Chip::Digits &dddigits )
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " CsHCAL2::DecodeTriggerGroups " << endl;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  if( debug ) cout << " std::multimap declaration OK " << endl;
  std::string tgname("HC02P1T1");
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
//  std::pair<m_it,m_it> m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  std::pair<m_it,m_it> m_range = dddigits.equal_range(tgname);
  if( debug ) cout << " dddigits.equal_range( " << tgname << " OK " << endl;

  // loop on all found digits
  int cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer1_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer1_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME1][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime1( time - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P1T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
//  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  m_range = dddigits.equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
     if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer1_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer1_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME2][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime2( time - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P2T1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
//  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  m_range = dddigits.equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P2T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
//  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  m_range = dddigits.equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P3T1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
//  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  m_range = dddigits.equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P3T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
//  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  m_range = dddigits.equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P4T1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
//  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  m_range = dddigits.equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer2_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer2_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME1][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime1( time - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC02P4T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
//  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  m_range = dddigits.equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer2_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer2_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME2][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime2( time - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

//   m_range = CsEvent::Instance()->getChipDigits().equal_range(GetName());
//   cnt = 0;
//   for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
//     cout << " HCAL2 for "<< GetName() <<" Didit #" << cnt++ << " found " << endl;


// Amplitudes info
  tgname = "HC02P1S0";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = dddigits.equal_range(tgname);
  if( debug ) cout << " dddigits.equal_range( " << tgname << " OK " << endl;

  // loop on all found digits
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipADC::Digit *d = dynamic_cast<const CS::ChipADC::Digit *>(d_it->second);
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = d->GetX();
    int32 y_digit = d->GetY();
    int32 e_digit = d->GetAmplitude();
    if(offset_2x2_trig >=0 )
    {
      int itrig =  x_digit+11*y_digit+offset_2x2_trig;
      if( itrig >= (int)trigger_groups.size() )
      {
    	 printf(" Error:: DecodeChipDigit: Calorimeter %s TriggerGroup %d is OUT of RANGE %zu \n",getName(),itrig,trigger_groups.size());
      }
      else
      {
         if( itrig >= (int)group_info[TrigGroup::CALIB][OLD].size() ) continue;
    	 trigger_groups[itrig].trig_group_data.SetAmpl( e_digit*group_info[TrigGroup::CALIB][OLD][itrig].GetMean() );
      }
    }
  }
  if( debug ) cout << " Loop over digits " << cnt << "  for " << tgname << "  OK " << endl;
//   This is trigger sum of 4x4 cells matrixes
//          printf("DecodeChipDigit:  4x4 TrigGroupLayer1 BOX = %d ESum=%d OF CALORIMETER %s \n",
//                                                                        x_digit-200,e_digit,getName());

  tgname = "HC02P1S1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = dddigits.equal_range(tgname);
  if( debug ) cout << " dddigits.equal_range( " << tgname << " OK " << endl;

  // loop on all found digits
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipADC::Digit *d = dynamic_cast<const CS::ChipADC::Digit *>(d_it->second);
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = d->GetX();
    int32 y_digit = d->GetY();
    int32 e_digit = d->GetAmplitude();
//  Check are this trigger groups initialised ?
    if(offset_layer1_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer1_trig;
      if( itrig >= (int)trigger_groups.size() )
      {
    	 printf(" Error:: DecodeChipDigit: Calorimeter %s TriggerGroup %d is OUT of RANGE %zu \n",getName(),itrig,trigger_groups.size());
      }
      else
      {
         if( itrig >= (int)group_info[TrigGroup::CALIB][OLD].size() ) continue;
    	 trigger_groups[itrig].trig_group_data.SetAmpl( e_digit*group_info[TrigGroup::CALIB][OLD][itrig].GetMean() );
      }
    }
  }
  if( debug ) cout << " Loop over digits " << cnt << "  for " << tgname << "  OK " << endl;

//   This is trigger sum of 4x4 cells matrixes
//          printf("DecodeChipDigit:  4x4 TrigGroupLayer2 BOX = %d ESum=%d OF CALORIMETER %s \n",
//                                                                        x_digit-300,e_digit,getName());
  tgname = "HC02P1S2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = dddigits.equal_range(tgname);
  if( debug ) cout << " dddigits.equal_range( " << tgname << " OK " << endl;

  // loop on all found digits
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipADC::Digit *d = dynamic_cast<const CS::ChipADC::Digit *>(d_it->second);
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = d->GetX();
    int32 y_digit = d->GetY();
    int32 e_digit = d->GetAmplitude();
    if(offset_layer2_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer2_trig;
      if( itrig >= (int)trigger_groups.size() )
      {
    	 printf(" Error:: DecodeChipDigit: Calorimeter %s TriggerGroup %d is OUT of RANGE %zu \n",getName(),itrig,trigger_groups.size());
      }
      else
      {
         if( itrig >= (int)group_info[TrigGroup::CALIB][OLD].size() ) continue;
    	 trigger_groups[itrig].trig_group_data.SetAmpl( e_digit*group_info[TrigGroup::CALIB][OLD][itrig].GetMean() );
      }
    }
  }
  if( debug ) cout << " Loop over digits " << cnt << " for " << tgname << "  OK " << endl;

//   This is trigger sum of 4x4 cells matrixes
//          printf("DecodeChipDigit:  4x4 TrigGroupLayer2 BOX = %d ESum=%d OF CALORIMETER %s \n",
//                                                                        x_digit-300,e_digit,getName());
  tgname = "HC02P1S3";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = dddigits.equal_range(tgname);
  if( debug ) cout << " dddigits.equal_range( " << tgname << " OK " << endl;
  // loop on all found digits
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipADC::Digit *d = dynamic_cast<const CS::ChipADC::Digit *>(d_it->second);
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = d->GetX();
    int32 y_digit = d->GetY();
    int32 e_digit = d->GetAmplitude();
    if(offset_layer3_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer3_trig;
      if( itrig >= (int)trigger_groups.size() )
      {
    	 printf(" Error:: DecodeChipDigit: Calorimeter %s TriggerGroup %d is OUT of RANGE %zu \n",getName(),itrig,trigger_groups.size());
      }
      else
      {
         if( itrig >= (int)group_info[TrigGroup::CALIB][OLD].size() ) continue;
    	 trigger_groups[itrig].trig_group_data.SetAmpl( e_digit*group_info[TrigGroup::CALIB][OLD][itrig].GetMean() );
      }
    }
  }
  if( debug ) cout << " Loop over digits " << cnt << " for " << tgname << "  OK " << endl;

//   This is trigger sum of 4x4 cells matrixes
//          printf("DecodeChipDigit:  4x4 TrigGroupLayer2 BOX = %d ESum=%d OF CALORIMETER %s \n",
//                                                                        x_digit-300,e_digit,getName());
  tgname = "HC02P1S4";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = dddigits.equal_range(tgname);
  if( debug ) cout << " dddigits.equal_range( " << tgname << " OK " << endl;
  // loop on all found digits
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipADC::Digit *d = dynamic_cast<const CS::ChipADC::Digit *>(d_it->second);
    if( debug ) cout << " HCAL2 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = d->GetX();
    int32 y_digit = d->GetY();
    int32 e_digit = d->GetAmplitude();
    if(offset_layer4_trig >=0 )
    {
      int itrig =  x_digit+5*y_digit+offset_layer4_trig;
      if( itrig >= (int)trigger_groups.size() )
      {
    	 printf(" Error:: DecodeChipDigit: Calorimeter %s TriggerGroup %d is OUT of RANGE %zu \n",getName(),itrig,trigger_groups.size());
      }
      else
      {
         if( itrig >= (int)group_info[TrigGroup::CALIB][OLD].size() ) continue;
    	 trigger_groups[itrig].trig_group_data.SetAmpl( e_digit*group_info[TrigGroup::CALIB][OLD][itrig].GetMean() );
      }
    }
  }
  if( debug ) cout << " Loop over digits " << cnt << " for " << tgname << "  OK " << endl;


  if( debug ) cout << " CsHCAL2::DecodeTriggerGroups Finished OK " << endl;

}
