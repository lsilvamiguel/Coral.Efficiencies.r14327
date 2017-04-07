/*!
   \file    CsHCAL1.cc
   \brief   Implementation class for HCAL1 COMPASS calorimeter
   \version $Revision: 1.37 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \author  Denis Murashev
   \date    $Date: 2011/02/01 22:05:52 $
*/
#include "coral_config.h"

#include "CsOpt.h"
#include "CsErrLog.h"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TF1.h"

#include "CsHCAL1.h"
#include "CsEvent.h"

#include "Reco/CalorimeterHist.h"
#include "Reco/CellDataRaw.h"

using namespace std;
using namespace Reco;

////////////////////////////////////////////////////////////////////////////////

CsHCAL1::CsHCAL1(const string &name, const string &geom_file)
 : CsCalorimeter(name, geom_file)
{
//   bool debug = true;

  InitXY();

  offset_layer1_trig = NGroups();
  for( int yg=4; yg!=-1; yg-- ) for( int xg=0; xg!=6; xg++ )
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
//    sprintf(tg_name," 4x4_Layer1_%d_%d Number%d",xg,4-yg,trigger_groups.size());
    sprintf(tg_name," 4x4_Layer1_%d_%d Number%zu",xg,yg,trigger_groups.size());
    TrigGroup tg(this,tg_name,g,1,xg,yg);
    if(g.size() > 0 ) trigger_groups.push_back(tg);
//    cout << " Trig group " << trigger_groups.size()-1 <<
//                                " x " << trigger_groups.back().trig_group_data.GetX() <<
//                                " y " << trigger_groups.back().trig_group_data.GetY() <<  endl;
//    cout << " Trig group " << trigger_groups.size()-1 << " N elements " << trigger_groups.back().NElements() <<
//                                                              " GroupName " << trigger_groups.back().GetName() << endl;
  }

  offset_layer2_trig = NGroups();
  for( int yg=3; yg!=-1; yg-- ) for( int xg=0; xg!=7; xg++ )
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
//    sprintf(tg_name," 4x4_Layer2_%d_%d Number%d",xg,3-yg,trigger_groups.size()-offset_layer2_trig);
    sprintf(tg_name," 4x4_Layer2_%d_%d Number%zu",xg,yg,trigger_groups.size()-offset_layer2_trig);
    TrigGroup tg(this, tg_name,g,2,xg,yg);
    if(g.size() > 0 ) trigger_groups.push_back(tg);
//    cout << " Trig group " << trigger_groups.size()-1 << " N elements " << trigger_groups.back().NElements() <<
//                                                              " GroupName " << trigger_groups.back().GetName() << endl;
  }

  Reco::StatInfo zero = Reco::StatInfo(0.,0.,0.);
  double time_group_default_calibration = -935;
  double energy_group_default_calibration = 0.05;
  for(size_t it=0; it!=NGroups(); it++)
  {
    group_info[TrigGroup::CALIB][OLD].push_back(Reco::StatInfo(1,energy_group_default_calibration,1.));
    group_info[TrigGroup::CALIB][NEW].push_back(zero);
    group_info[TrigGroup::TIME1][OLD].push_back(Reco::StatInfo(1,time_group_default_calibration,1.));
    group_info[TrigGroup::TIME1][NEW].push_back(zero);
    group_info[TrigGroup::TIME2][OLD].push_back(Reco::StatInfo(1,time_group_default_calibration,1.));
    group_info[TrigGroup::TIME2][NEW].push_back(zero);
  }

  if( options.print_general_info ) PrintGeneralInfo();

}

////////////////////////////////////////////////////////////////////////////////

void  CsHCAL1::InitOptions( void )
{
  CsCalorimeter::InitOptions();

  options.calo_type                     = Hadronic;
  options.particle_default              = CalorimeterParticle::PION_PLUS;
  options.default_calibration           = 0.060;
  options.readout_sparsified_           = true;
  options.delta_spars_mode_             = 3.;

  // MC options
  options.mc_default_calibration     = 0.060;
  options.mc_make_real_digitization  = true;
  options.mc_make_fiadc_digitization = true;
  options.mc_fiadc_sparce_delta      = 3.;

}

////////////////////////////////////////////////////////////////////////////////

void  CsHCAL1::Initialize( void )
{
  CsCalorimeter::Initialize();
//   if( fem_ != NULL ) fem_->Init( this );
}

////////////////////////////////////////////////////////////////////////////////

void CsHCAL1::DecodeChipDigit(const CS::Chip::Digit &digit)
{
  bool debug_sadc = false;
//   bool debug_sadc = true;
  const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  const CS::ChipADC::Digit *d = dynamic_cast<const CS::ChipADC::Digit *>(&digit);
  const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit *>(&digit);
  if( (dt==NULL) && (d==NULL) && (ds==NULL) )
  {
    cout << "CsHCAL1::DecodeChipDigit(): Wrong digit type"  << endl;
    throw CS::Exception("CsCalorimeter::DecodeChipDigit(): Wrong digit type");
  }
  if( d != NULL )
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
       if( trigger_groups.size() == 0 ) return;

       if( y_digit >= 100 && y_digit <= 104 && x_digit >= 0 )
       {
//   This is trigger sum of 4x4 cells matrixes
//          printf("DecodeChipDigit:  4x4 TrigGroupLayer1 BOX = %d ESum=%d OF CALORIMETER %s \n",
//                                                                        x_digit-200,e_digit,getName());
//  Check are this trigger groups initialised ?
          if(offset_layer1_trig >=0 )
          {
            y_digit = 4-(y_digit-100);
            int itrig =  x_digit+6*y_digit+offset_layer1_trig;
            if(itrig >= 16 ) itrig -=2; // vo-kak blin!!!
            if( itrig >= (int)trigger_groups.size() )
            {
               printf(" Error:: DecodeChipDigit: Calorimeter %s TriggerGroup %d is OUT of RANGE %zu \n",
                                                                   getName(),itrig,trigger_groups.size());
            }
            else
            {
               trigger_groups[itrig].trig_group_data.SetAmpl( e_digit*group_info[TrigGroup::CALIB][OLD][itrig].GetMean() );
            }
          }
       }
       else if( y_digit >= 200 && y_digit <= 203 && x_digit >= 0 )
       {
//   This is trigger sum of 4x4 cells matrixes
//          printf("DecodeChipDigit:  4x4 TrigGroupLayer2 BOX = %d ESum=%d OF CALORIMETER %s \n",
//                                                                        x_digit-300,e_digit,getName());
//  Check are this trigger groups initialised ?
          if(offset_layer2_trig >=0 )
          {
            y_digit = 3-(y_digit-200);
            int itrig =  x_digit+7*y_digit+offset_layer2_trig;
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
//     cout << "CsHCAL1::DecodeChipDigit(): data size "  << raw_data_store.size() << endl;


  }
  else if( dt!=NULL )
  {
    if( trigger_groups.size() == 0 ) return;
    int32 ch_pos_digit = dt->GetChannelPos();
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();

//    double t0 = CsEvent::Instance()->getTriggerTime();
//    double time = CS::ChipF1::TimeDifference(dt->GetAmplitude(),CsEvent::Instance()->getTriggerTime());
  double time=dt->GetTimeDecoded();
      #if calorim_debug>1
      printf("CsHCAL1::DecodeChipDigit: Calorimeter %s Ch %d X %d Y %d ChPos %d Time %d \n",
                                                getName(),ch_digit,x_digit,y_digit,ch_pos_digit,t_digit);
      #endif

    if(ch_pos_digit == 1 && y_digit >= 100 && y_digit <= 104 )
    {
      if(offset_layer1_trig >=0 )
      {
//          cout << "Layer 1 T1 xdigit " << x_digit  << " ydigit " << y_digit-100 ;
         y_digit = 4-(y_digit-100);
         int itrig =  x_digit+6*y_digit+offset_layer1_trig;
            if(itrig >= 16 ) itrig -=2; // vo-kak blin!!!
//         cout << " L " << trigger_groups[itrig].trig_group_data.GetLayer()
//              << " x " << trigger_groups[itrig].trig_group_data.GetX()
//              << " y " <<  trigger_groups[itrig].trig_group_data.GetY() << endl;
         trigger_groups[itrig].trig_group_data.SetTime1( time - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
//        printf("CsHCAL1::DecodeChipDigit: Calorimeter %s layer 1 Ch %d  Time %d \n",
//                                                                  getName(),itrig-offset_layer1_trig,t_digit);
      }
    }
    if( ch_pos_digit == 2 && y_digit >= 100 && y_digit <= 104 )
    {
      if(offset_layer1_trig >=0 )
      {
//          cout << "Layer 1 T2 xdigit " << x_digit  << " ydigit " << y_digit-100 ;
          y_digit = 4-(y_digit-100);
         int itrig =  x_digit+6*y_digit+offset_layer1_trig;
            if(itrig >= 16 ) itrig -=2; // vo-kak blin!!!
//         cout << " L " << trigger_groups[itrig].trig_group_data.GetLayer()
//              << " x " << trigger_groups[itrig].trig_group_data.GetX()
//              << " y " <<  trigger_groups[itrig].trig_group_data.GetY() << endl;
         trigger_groups[itrig].trig_group_data.SetTime2( time  - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
//        printf("CsHCAL1::DecodeChipDigit: Calorimeter %s layer 1 Ch %d  Time %d \n",
//                                                                  getName(),itrig-offset_layer1_trig,t_digit);
      }
    }
    if( ch_pos_digit == 1 && y_digit >= 200 && y_digit <= 203 && x_digit >= 0 )
    {
      if(offset_layer2_trig >=0 )
      {
//          cout << "Layer 2 T1 xdigit " << x_digit  << " ydigit " << y_digit-200 ;
         y_digit = 3-(y_digit-200);
         int itrig =  x_digit+7*y_digit+offset_layer2_trig;
//          cout << " L " << trigger_groups[itrig].trig_group_data.GetLayer()
//              << " x " << trigger_groups[itrig].trig_group_data.GetX()
//              << " y " <<  trigger_groups[itrig].trig_group_data.GetY() << endl;

         trigger_groups[itrig].trig_group_data.SetTime1( time  - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
//        printf("CsHCAL1::DecodeChipDigit: Calorimeter %s layer 2 Ch %d  Time %d \n",
//                                                                  getName(),itrig-offset_layer2_trig,t_digit);
      }
    }
    if( ch_pos_digit == 2 && y_digit >= 200 && y_digit <= 203 && x_digit >= 0 )
    {
      if(offset_layer2_trig >=0 )
      {
//          cout << "Layer 2 T2 xdigit " << x_digit  << " ydigit " << y_digit-200 ;
         y_digit = 3-(y_digit-200);
         int itrig =  x_digit+7*y_digit+offset_layer2_trig;
//          cout << " L " << trigger_groups[itrig].trig_group_data.GetLayer()
//              << " x " << trigger_groups[itrig].trig_group_data.GetX()
//              << " y " <<  trigger_groups[itrig].trig_group_data.GetY() << endl;
         trigger_groups[itrig].trig_group_data.SetTime2( time - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
//        printf("CsHCAL1::DecodeChipDigit: Calorimeter %s layer 2 Ch %d  Time %d \n",
//                                                                  getName(),itrig-offset_layer2_trig,t_digit);
      }
    }


  }
  else if( ds!=NULL ) // ChipSADC
  {
    if( debug_sadc ) cout << " SADC digit in " << GetName() << endl;
    DecodeChipSADCDigit(*ds);
  }
}

////////////////////////////////////////////////////////////////////////////////

void CsHCAL1::DecodeTriggerGroups( void )
{
  bool debug = true;
//   bool debug = false;
  if( debug ) cout << " CsHCAL1::DecodeTriggerGroups " << endl;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  if( debug ) cout << " std::multimap declaration OK " << endl;
  std::string tgname("HC02P1T1");
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  std::pair<m_it,m_it> m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  if( debug ) cout << " std::CsEvent::Instance()->getChipDigits() .equal_range( " << tgname << " OK " << endl;
  int cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( dt == NULL ) continue;
    if( debug ) cout << " HCAL1 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer1_trig >=0 )
    {
      y_digit = 4-y_digit;
      int itrig =  x_digit+6*y_digit+offset_layer1_trig;
      if(itrig >= 16 ) itrig -=2; // vo-kak blin!!!
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME1][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime1( time - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC01P1T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( dt == NULL ) continue;
    if( debug ) cout << " HCAL1 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer1_trig >=0 )
    {
      y_digit = 4-y_digit;
      int itrig =  x_digit+6*y_digit+offset_layer1_trig;
      if(itrig >= 16 ) itrig -=2; // vo-kak blin!!!
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME2][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime2( time  - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC01P2T1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( dt == NULL ) continue;
    if( debug ) cout << " HCAL1 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer2_trig >=0 )
    {
      y_digit = 3-y_digit;
      int itrig =  x_digit+7*y_digit+offset_layer2_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME1][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime1( time  - group_info[TrigGroup::TIME1][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC01P2T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    const CS::ChipF1::Digit *dt = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
    if( dt == NULL ) continue;
    if( debug ) cout << " HCAL1 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
    int32 x_digit = dt->GetX();
    int32 y_digit = dt->GetY();
    double time=dt->GetTimeDecoded();
    if(offset_layer2_trig >=0 )
    {
      y_digit = 3-y_digit;
      int itrig =  x_digit+7*y_digit+offset_layer2_trig;
      if(itrig < 0 || itrig >= (int)trigger_groups.size() ) continue;
      if( itrig >= (int)group_info[TrigGroup::TIME2][OLD].size() ) continue;
      trigger_groups[itrig].trig_group_data.SetTime2( time - group_info[TrigGroup::TIME2][OLD][itrig].GetMean());
    }
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC01P3T1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL1 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC01P3T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL1 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC01P4T1";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL1 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;

  tgname = "HC01P4T2";
  if( debug ) cout << " Try to Loop over digits for " << tgname << endl;
  m_range = CsEvent::Instance()->getChipDigits().equal_range(tgname);
  cnt = 0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    if( debug ) cout << " HCAL1 for "<< tgname <<" Didit #" << cnt++ << " found " << endl;
  }
  if( debug ) cout << " Loop over digits for " << tgname << "  OK " << endl;
  if( debug ) cout << " CsHCAL1::DecodeTriggerGroups Finished OK " << endl;

}
