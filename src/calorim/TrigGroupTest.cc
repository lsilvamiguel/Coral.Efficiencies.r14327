/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/TrigGroupTest.cc,v $
   $Date: 2010/06/29 15:47:03 $
   $Revision: 1.21 $
   -------------------------------------------------------------------------
   \file    TrigGroupTest.cc
   \brief   Test function for trigger and timing stuff of the calorimeter
   \version $  $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \author  Denis Murashev
   \date    $ $

*/

// --- Standard C/C++ library ---
#include <cassert>
#include <cmath>
#include <iostream>

#include "coral_config.h"
// --- Histo library ---
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "Reco/myROOT_utils.h"

// --- Internal files ---

#include "CsEvent.h"

#include "CsCalorimeterHist.h"
#include "CsCalorimeter.h"

#include "Reco/CellDataRaw.h"

using namespace std;

void CsCalorimeter::TrigGroupTest(void)
{
  Reco::Calorimeter::TrigGroupTest(); // create a parallel stuff

  bool ok=false;
  if( !options.fill_histos )
    return;

//  cout << " Now in TrigGroupTest for " << GetName() << endl;

  if( gDirectory==NULL )
  {
    cerr << "CsCalorimeter::TrigGroupTest():  ROOT file was not opend.\n"
         << "                                 Please fix this .\n";
    return;
  }
  if(NGroups() <=0 ) return;
  TDirectory *dir_save=gDirectory; // Save current ROOT directory.

  if( cs_calo_hist==NULL )
  {
//  cout << " Now in TrigGroupTest for BOOKING " << GetName() << endl;
    cs_calo_hist = new CsCalorimeterHist;
    CsCalorimeterHist &h=*cs_calo_hist;  // create a short name for calo_hist

    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"%s_TriGroup_test",GetName().c_str());
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
    assert( NULL != h.root_dir);
    ok=(h.root_dir->cd());
    assert(ok);


    char hist_name[132];

    sprintf(hist_name,"%s:: Amount of Trigger Groups with Hits ",GetName().c_str());
    h.h1_NHits = myROOT_utils::TH1D_checked("NHits",hist_name,50,0,50);
    sprintf(hist_name,"%s:: Trigger Groups with Hits ",GetName().c_str());
    h.h1_HitedGroups = myROOT_utils::TH1D_checked("HitedGroups",hist_name,NGroups(),0,NGroups());
    sprintf(hist_name,"%s:: Groups vs Groups ",GetName().c_str());
    h.h2_GroupVsGroup = myROOT_utils::TH2D_checked("GroupVsGroup",hist_name,NGroups(),0,NGroups(),NGroups(),0,NGroups());
    sprintf(hist_name,"%s:: Groups vs Cells X ",GetName().c_str());
    h.h2_GroupVsCellX = myROOT_utils::TH2D_checked("GroupVsCellX",hist_name,NGroups(),0,NGroups(),100,(*this).GetXmin(),(*this).GetXmax());
    sprintf(hist_name,"%s:: Groups vs Cells Y ",GetName().c_str());
    h.h2_GroupVsCellY = myROOT_utils::TH2D_checked("GroupVsCellY",hist_name,NGroups(),0,NGroups(),100,(*this).GetYmin(),(*this).GetYmax());
    sprintf(hist_name,"%s:: T1Groups vs EGroups ",GetName().c_str());
    h.h2_T1GroupVsEGroup = myROOT_utils::TH2D_checked("T1GroupVsEGroup",hist_name,NGroups(),0,NGroups(),NGroups(),0,NGroups());
    sprintf(hist_name,"%s:: T2Groups vs EGroups ",GetName().c_str());
    h.h2_T2GroupVsEGroup = myROOT_utils::TH2D_checked("T2GroupVsEGroup",hist_name,NGroups(),0,NGroups(),NGroups(),0,NGroups());
    sprintf(hist_name,"%s:: TGroups vs logic_cell E_wighted ",GetName().c_str());
    h.h2_TGroupVsNCellE = myROOT_utils::TH2D_checked("TGroupVsNCellE",hist_name,NGroups(),0,NGroups(),NCells(),0,NCells());

// Create subdirectory
    ok=(h.root_dir->cd());
    assert(ok);
    sprintf(dir_name,"%s_GeometryTest",GetName().c_str());
    h.geom_test_dir = myROOT_utils::TDirectory_checked(dir_name);
    assert( NULL!=h.geom_test_dir);
    ok=(h.geom_test_dir->cd());
    assert(ok);

    for( size_t i=0; i< NGroups(); i++ )
    {
      char name[132];
      sprintf(name,"Geom%3.3zu",i);
      sprintf(hist_name,"%s::%s:: Geometry of Group %zu ",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h2_Geom[i] = myROOT_utils::TH2D_checked(name,hist_name,100,(*this).GetXmin(),(*this).GetXmax(),100,(*this).GetYmin(),(*this).GetYmax());
    }

    for( size_t i=0; i< NGroups(); i++ )
    {
      char name[132];
      sprintf(name,"SearchADC%3.3zu",i);
      sprintf(hist_name,"%s::%s::Search cell of ADC  Group %zu ",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h2_SearchADC[i] = myROOT_utils::TH2D_checked(name,hist_name,100,(*this).GetXmin(),(*this).GetXmax(),100,(*this).GetYmin(),(*this).GetYmax());
    }
    for( size_t i=0; i< NGroups(); i++ )
    {
      char name[132];
      sprintf(name,"SearchT1%3.3zu",i);
      sprintf(hist_name,"%s::%s::Search cell of T1  Group %zu ",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h2_SearchT1[i] = myROOT_utils::TH2D_checked(name,hist_name,100,(*this).GetXmin(),(*this).GetXmax(),100,(*this).GetYmin(),(*this).GetYmax());
    }
    for( size_t i=0; i< NGroups(); i++ )
    {
      char name[132];
      sprintf(name,"SearchT2%3.3zu",i);
      sprintf(hist_name,"%s::%s::Search cell of T2  Group %zu ",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h2_SearchT2[i] = myROOT_utils::TH2D_checked(name,hist_name,100,(*this).GetXmin(),(*this).GetXmax(),100,(*this).GetYmin(),(*this).GetYmax());
    }


    ok=(h.root_dir->cd());
    assert(ok);
// Create subdirectory
    sprintf(dir_name,"%s_AllGroups",GetName().c_str());
    h.all_groups_dir=myROOT_utils::TDirectory_checked(dir_name);
    assert( NULL!=h.all_groups_dir );
    ok=(h.all_groups_dir->cd());
    assert(ok);

    for( size_t i=0; i< NGroups(); i++ )
    {
      char name[132];
      sprintf(name,"Egr%3.3zu",i);
//      cout << " Now in TrigGroupTest for BOOKING " << GetName() << " i=" << i << " name " <<  name.str() <<   endl;
      sprintf(hist_name,"%s::%s:: Energy in Group %zu ",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h1_Egroups[i] = myROOT_utils::TH1D_checked(name,hist_name,100,0,50);
    }

    ok=(h.root_dir->cd());
    assert(ok);
// Create another subdirectory
    sprintf(dir_name,"AllGroupsCorrelations");
    h.all_gcorr_dir=myROOT_utils::TDirectory_checked(dir_name);
    assert( NULL!=h.all_gcorr_dir);
    ok=(h.all_gcorr_dir->cd());
    assert(ok);

    for( size_t i=0; i< NGroups(); i++ )
    {
// # warning : TrigGroupTest: small problem with names
      char name[132],name1[132],name2[132];
      sprintf(name,"Egr%3.3zu",i);
//      cout << " Now in TrigGroupTest for BOOKING " << GetName() << " i=" << i << " name " <<  name.str() <<   endl;
      sprintf(hist_name,"%s::%s:: TrigADCEnergy vs Esumm in Group %zu ",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h2_EtgEcm[i] = myROOT_utils::TH2D_checked(name,hist_name,100,0,100,100,0,100);

      sprintf(name1,"EgrT1%3.3zu",i);
      sprintf(hist_name,"%s::%s:: TrigADCEnergy vs Esumm in Group %zu if T1 on",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h2_EtgEcmT1[i] = myROOT_utils::TH2D_checked(name1,hist_name,100,0,100,100,0,100);

      sprintf(name2,"EgrT2%3.3zu",i);
      sprintf(hist_name,"%s::%s:: TrigADCEnergy vs Esumm in Group %zu if T2 on",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h2_EtgEcmT2[i] = myROOT_utils::TH2D_checked(name2,hist_name,100,0,100,100,0,100);
    }
    ok=(h.root_dir->cd());
    assert(ok);

// Create subdirectory  for timing
    sprintf(dir_name,"%s_Timing",GetName().c_str());
    h.timing_dir=myROOT_utils::TDirectory_checked(dir_name);
    assert( NULL!=h.timing_dir);
    ok=(h.timing_dir->cd());
    assert(ok);

    int   n_time_bins = 200;
    float time_min = -100 , time_max = 100;

    for( size_t i=0; i< NGroups(); i++ )
    {
      char name1[132];
      sprintf(name1,"T1g%3.3zu",i);
      sprintf(hist_name,"%s::%s:: Timing in Group %zu ",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h2_Timing1[i] = myROOT_utils::TH2D_checked(name1,hist_name,n_time_bins,time_min,time_max,50,0.,100.);
      char name2[132];
      sprintf(name2,"T2g%3.3zu",i);
      sprintf(hist_name,"%s::%s:: Timing in Group %zu ",GetName().c_str(),trigger_groups[i].GetName().c_str(),i);
      h.h2_Timing2[i] = myROOT_utils::TH2D_checked(name2,hist_name,n_time_bins,time_min,time_max,50,0.,100.);
    }

    ok=(h.root_dir->cd());
    assert(ok);

// Create subdirectory  for the cells
    sprintf(dir_name,"%s_Cells",GetName().c_str());
    h.cells_dir=myROOT_utils::TDirectory_checked(dir_name);
    assert( NULL!=h.cells_dir);
    ok=(h.cells_dir->cd());
    assert(ok);

    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"Cell%s",GetCellName(i).c_str());
//      cout << " Now in TrigGroupTest for BOOKING " << GetName() << " i=" << i << " name " <<  name.str() <<   endl;
      sprintf(hist_name," Energy in cell %s ",GetCellName(i).c_str());
      h.h1_ECell[i] = myROOT_utils::TH1D_checked(name,hist_name,500,0.,50.);
    }

    ok=(h.root_dir->cd());
    assert(ok);

  }

////////////////////////////////////////  Filling histograms  ////////////////////////////////////////////////////////////////
//  if(calib_data_store.size() > 30) return;

//  cout << " Now in TrigGroupTest for Filling " << GetName() << endl;

  CsCalorimeterHist &h=*cs_calo_hist;  // create a short name for calo_hist
  ok=(h.root_dir->cd());
  assert(ok);
  for( size_t it1=0; it1!=signals_.size(); it1++ )
  {
     double e=signals_[it1].GetEnergy();
     size_t icell = signals_[it1].GetCellIdx();
     h.h1_ECell[icell]->Fill(e);
  }
  for( size_t it=0; it!=NGroups(); it++ )
  {
     for( vector<size_t>::const_iterator it1=trigger_groups[it].GetCells().begin(); it1!=trigger_groups[it].GetCells().end(); it1++ )
     {
             double x = cells[*it1].GetX();
             double y = cells[*it1].GetY();
             h.h2_Geom[it]->Fill(x,y);
     }
  }

  int nhit_groups =0;
  for( size_t it=0; it!=NGroups(); it++ )
  {
    double Egr = trigger_groups[it].trig_group_data.GetAmpl();
    double Egr_sum = trigger_groups[it].trig_group_data.GetSumADC();
    double T1 = trigger_groups[it].trig_group_data.GetTime1();
    double T2 = trigger_groups[it].trig_group_data.GetTime2();
    h.h2_EtgEcm[it]->Fill(Egr,Egr_sum);
    if(T1!=0) h.h2_EtgEcmT1[it]->Fill(Egr,Egr_sum);
    if(T2!=0) h.h2_EtgEcmT2[it]->Fill(Egr,Egr_sum);
    if( Egr > 0. )
    {
      h.h1_HitedGroups->Fill(it);
      h.h1_Egroups[it]->Fill(Egr);
//      h.h2_EtgEcm[it]->Fill(Egr,Egr_sum);
      nhit_groups++;
// Group -Group correlation
      for( size_t it1=0; it1!=NGroups(); it1++ )
      {
        if(it1 >it )
        {
          double Egr1 = trigger_groups[it1].trig_group_data.GetAmpl();
          if( Egr1 > 0. ) h.h2_GroupVsGroup->Fill(it,it1);
        }
      }
// Group - CellsXY correlation
      for( size_t it1=0; it1!=signals_.size(); it1++ )
      {
          double e=signals_[it1].GetEnergy();
          if(e > 0.)
          {
             int icell = signals_[it1].GetCellIdx();
             double x = cells[icell].GetX();
             double y = cells[icell].GetY();
             h.h2_GroupVsCellX->Fill(it,x);
             h.h2_GroupVsCellY->Fill(it,y);
          }
      }
// Search for cells correlated with  Group ADC
      double EgrSearcADC=5.;
      double egam=2.5;
//      cout << GetName() << " group " << it << " Egr " << Egr << " Ng " << reco_particles.size() << endl;
      if( Egr > EgrSearcADC && reco_particles.size() > 0 && reco_particles.size() <= 3 )
      {
//         for( size_t it1=0; it1!=calib_data_store.size(); it1++ )
//         {
//           double e=calib_data_store[it1].GetAmplitude();
//           if(e > 0.)
//           {
//              int icell = calib_data_store[it1].GetAddress();
//              double x = cells[icell].GetX();
//              double y = cells[icell].GetY();
//              h.h2_SearchADC[it]->Fill(x,y,e);
//           }
//         }

        for( size_t ng=0; ng < reco_particles.size(); ng++ )
        {
          double e=  reco_particles[ng].GetE();
          const vector<size_t> &hitedcells =reco_particles[ng].GetMainCells();
          if(e > egam)
          {
             int icell = hitedcells[0];
             double x = cells[icell].GetX();
             double y = cells[icell].GetY();
             h.h2_SearchADC[it]->Fill(x,y,e);
          }
        }
      }
    }

    if( T1 != 0. )
    {
// Group -Group Time-Amplitude correlation
      for( size_t it1=0; it1!=NGroups(); it1++ )
      {
         double Egr1 = trigger_groups[it1].trig_group_data.GetAmpl();
         h.h2_T1GroupVsEGroup->Fill(it,it1,Egr1);

      }
// Search for cells correlated with  Group T1
      for( size_t it1=0; it1!=signals_.size(); it1++ )
      {
          double e=signals_[it1].GetEnergy();
          if(e > 0.)
          {
             int icell = signals_[it1].GetCellIdx();
             double x = cells[icell].GetX();
             double y = cells[icell].GetY();
             h.h2_SearchT1[it]->Fill(x,y,e);
          }
      }
    }

    if( T2 != 0. )
    {
// Group -Group Time-Amplitude correlation
      for( size_t it1=0; it1!=NGroups(); it1++ )
      {
         double Egr1 = trigger_groups[it1].trig_group_data.GetAmpl();
         h.h2_T2GroupVsEGroup->Fill(it,it1,Egr1);

      }
      for( size_t it1=0; it1!=signals_.size(); it1++ )
      {
          double e=signals_[it1].GetEnergy();
          if(e > 0.)
          {
             int icell = signals_[it1].GetCellIdx();
             h.h2_TGroupVsNCellE->Fill(it,icell,e);
          }
      }
// Search for cells correlated with  Group T2
      for( size_t it1=0; it1!=signals_.size(); it1++ )
      {
          double e=signals_[it1].GetEnergy();
          if(e > 0.)
          {
             int icell = signals_[it1].GetCellIdx();
             double x = cells[icell].GetX();
             double y = cells[icell].GetY();
             h.h2_SearchT2[it]->Fill(x,y,e);
          }
      }
     }
  }

  h.h1_NHits->Fill(nhit_groups);

  for( size_t it=0; it!=NGroups(); it++ )
  {
    double time1 = trigger_groups[it].trig_group_data.GetTime1();
    double time2 = trigger_groups[it].trig_group_data.GetTime2();
    double Egr = trigger_groups[it].trig_group_data.GetAmpl();
    if(time1 != 0 ) h.h2_Timing1[it]->Fill(time1,Egr);
    if(time2 != 0 ) h.h2_Timing2[it]->Fill(time2,Egr);
  }

  dir_save->cd();
//  cout << " Now in TrigGroupTest GetOut " << GetName() << endl;
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::DecodingTestSADC(void)
{
//  cout << " CsCalorimeter::DecodingTestSADC" << GetName() << endl;
  bool ok=false;
  if( !options.fill_histos )
    return;

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
  if( test_histo_sadc_ == NULL )
  {
    test_histo_sadc_ = new TestHistoSADC();
    TestHistoSADC &h = *test_histo_sadc_;  // create a short name

    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    TDirectory *dir = myROOT_utils::TDirectory_checked(dir_name);
    ok=(dir->cd());
    assert(ok);
    sprintf(dir_name,"%s_DecodingTestSADC",GetName().c_str());
    dir = myROOT_utils::TDirectory_checked(dir_name);
    for( size_t i=0; i< subsets_list2store_sadc_histo_.size()+1; i++ )
    {
      ok=(dir->cd());
      assert(ok);
      if( i > 0 )
      {
        int isub = subsets_list2store_sadc_histo_[i-1];
        sprintf(dir_name,"SubSet_%s",GetSubSets()[isub].GetName().c_str());
        TDirectory *subdir = myROOT_utils::TDirectory_checked(dir_name);
        ok=(subdir->cd());
        assert(ok);
      }

      int max_amp = 4000;
      char hist_name[132];
      char name[132];
      sprintf(name,"NSamplesSADC");
      sprintf(hist_name," Amount of samples ");
      h.h1_NSamplesSADC.push_back( myROOT_utils::TH1D_checked(name,hist_name,100,0.,100.) );
      sprintf(name,"SizeSamplesSADC");
      sprintf(hist_name," Sample size ");
      h.h1_SizeSamplesSADC.push_back( myROOT_utils::TH1D_checked(name,hist_name,100,0.,100.) );
      sprintf(name,"ChMaxAmp");
      sprintf(hist_name," Max amplitude vs Max channel ");
      h.h2_ChMaxAmp.push_back( myROOT_utils::TH2D_checked(name,hist_name,100,0.,100.,100,0.,max_amp) );
      sprintf(name,"SummMax");
      sprintf(hist_name," Max amplitude vs Sample summ ");
      h.h2_SummMax.push_back( myROOT_utils::TH2D_checked(name,hist_name,100,0.,10000.,100,0.,max_amp) );
      sprintf(name,"SignalMax");
      sprintf(hist_name," Max amplitude vs Sample signal ");
      h.h2_SignalMax.push_back( myROOT_utils::TH2D_checked(name,hist_name,100,0.,10000.,100,0.,max_amp) );
      sprintf(name,"PedestalsSADC");
      sprintf(hist_name," Pedestals of SADC ");
      h.h1_Ped.push_back( myROOT_utils::TH1D_checked(name,hist_name,1000,0.,1000.) );
      sprintf(name,"SummSADC");
      sprintf(hist_name," Summ of SADC pedestal subtracted");
      h.h1_Summ.push_back( myROOT_utils::TH1D_checked(name,hist_name,10000,0.,10000.) );
      sprintf(name,"PedestalsSADC");
      sprintf(hist_name," Signal of SADC pedestal subtracted ");
      h.h1_Signal.push_back( myROOT_utils::TH1D_checked(name,hist_name,10000,0.,10000.) );
      sprintf(name,"Shape10SADC");
      sprintf(hist_name," Shape of SADC signals summ gt 10 lt 100 ");
      h.p1_Shape10SADC.push_back( myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.) );
      sprintf(name,"Shape100SADC");
      sprintf(hist_name," Shape of SADC signals summ gt 100 lt 500");
      h.p1_Shape100SADC.push_back( myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.) );
       sprintf(name,"Shape500SADC");
      sprintf(hist_name," Shape of SADC signals summ gt 500 lt 100 ");
      h.p1_Shape500SADC.push_back( myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.) );
      sprintf(name,"Shape1000SADC");
      sprintf(hist_name," Shape of SADC signals summ gt 1000 ");
      h.p1_Shape1000SADC.push_back( myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.) );

      sprintf(name,"ShapeNM10SADC");
      sprintf(hist_name," ShapeNM  of SADC signals summ gt 10 lt 100 ");
      h.p1_ShapeNM10SADC.push_back( myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.) );
      sprintf(name,"ShapeNM100SADC");
      sprintf(hist_name," ShapeNM of SADC signals summ gt 100 lt 500");
      h.p1_ShapeNM100SADC.push_back( myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.) );
       sprintf(name,"ShapeNM500SADC");
      sprintf(hist_name," ShapeNM of SADC signals summ gt 500 lt 100 ");
      h.p1_ShapeNM500SADC.push_back( myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.) );
      sprintf(name,"ShapeNM1000SADC");
      sprintf(hist_name," ShapeNM of SADC signals summ gt 1000 ");
      h.p1_ShapeNM1000SADC.push_back( myROOT_utils::TProfile_checked(name,hist_name,128,0.,128.) );

      sprintf(name,"TimeSADC");
      sprintf(hist_name," Time for SADC ");
      h.h1_TimeSADC.push_back( myROOT_utils::TH1D_checked(name,hist_name,1000,-50.,50.) );
      sprintf(name,"TimeSADCvsE");
      sprintf(hist_name," E vs Time for SADC ");
      h.h2_TimeSADCvsE.push_back( myROOT_utils::TH2D_checked(name,hist_name,50,-50.,50.,50,0.,50.) );
      sprintf(name,"TimeSADCvsTCS");
      sprintf(hist_name," TCS corr vs Time for SADC ");
      h.h2_TimeSADCvsTCS.push_back( myROOT_utils::TH2D_checked(name,hist_name,50,-50.,50.,50,-50.,50.) );
    }
  }

//   cout << " Calorimeter " << GetName() << " N SubSets " <<  GetSubSets().size() << endl;
//   for ( int isub=0; isub < GetSubSets().size(); isub++ )
//   {
//     cout << " SubSet " <<  GetSubSets()[isub].GetName() << endl;
//   }

  TestHistoSADC &h = *test_histo_sadc_;  // create a short name

  h.h1_NSamplesSADC[0]->Fill( (double)(sadc_samples_.size()) );

//  for( size_t it=0; it < sadc_samples_.size(); it++ )
  for( map< int, const vector<CS::uint16> * >::const_iterator its=sadc_samples_.begin(); its != sadc_samples_.end(); its++ )
  {
//    int ic = sadc_samples_[it].first;
    int ic = its->first;
    vector < int >ssl;
    ssl.push_back(0);
    for( size_t is = 0; is < subsets_list2store_sadc_histo_.size(); is++ )
    {
      int isub = subsets_list2store_sadc_histo_[is];
      if(  GetSubSets()[isub].IsMyCell(ic) ) ssl.push_back(is+1);
    }

//    const vector<CS::uint16> &sample = *sadc_samples_[it].second;
    const vector<CS::uint16> &sample = *(its->second);
    double ped =0.;

    for( size_t iis=0; iis<ssl.size(); iis++)
    {
      h.h1_SizeSamplesSADC[ ssl[iis] ]->Fill( (double)( sample.size() ) );
    }

    for( int is=sadc_ped_min_; is <= sadc_ped_max_; is++)
    {
      ped += sample[is];
    }
    ped /= (sadc_ped_max_-sadc_ped_min_+1);
    for( size_t iis=0; iis<ssl.size(); iis++)
    {
      h.h1_Ped[ ssl[iis] ]->Fill(ped);
    }

    double amp_max = -1000.;
    double summ = 0.;
    int ch_max = -1;

    for( int is=0; is < (int)(sample.size()-1); is++)
    {
      double amp = double(sample[is])-ped;
      if( amp > amp_max )
      {
        amp_max = amp;
        ch_max = is;
      }
      summ += amp;
    }

    for( size_t iis=0; iis<ssl.size(); iis++)
    {
      h.h2_ChMaxAmp[ ssl[iis] ] ->Fill(double(ch_max),amp_max);

      h.h1_Summ[ ssl[iis] ]->Fill(summ);
      h.h2_SummMax[ ssl[iis] ]->Fill(summ,amp_max);
    }

    double signal = 0.;
    for( int is=sadc_signal_min_; is <= sadc_signal_max_; is++)
    {
      signal += sample[is] - ped;
    }

    for( size_t iis=0; iis<ssl.size(); iis++)
    {
      h.h1_Signal[ ssl[iis] ]->Fill(signal);
      h.h2_SignalMax[ ssl[iis] ]->Fill(signal,amp_max);
    }


    if( summ > 10. && summ <= 100. )
    {
      for( int is=0; is < (int)(sample.size()-1); is++)
      {
        for( size_t iis=0; iis<ssl.size(); iis++)
        {
//          h.p1_Shape100SADC[ ssl[iis] ]->Fill( double(is), (double(sample[is])-ped)/summ );
          h.p1_Shape10SADC[ ssl[iis] ]->Fill( double(is), (double(sample[is])-ped)/amp_max );
          h.p1_ShapeNM10SADC[ ssl[iis] ]->Fill( double(is-ch_max+64), (double(sample[is])-ped)/amp_max );
        }
      }
    }
    if( summ > 100. && summ <= 500. )
    {
      for( int is=0; is < (int)(sample.size()-1); is++)
      {
        for( size_t iis=0; iis<ssl.size(); iis++)
        {
//          h.p1_Shape100SADC[ ssl[iis] ]->Fill( double(is), (double(sample[is])-ped)/summ );
          h.p1_Shape100SADC[ ssl[iis] ]->Fill( double(is), (double(sample[is])-ped)/amp_max );
          h.p1_ShapeNM100SADC[ ssl[iis] ]->Fill( double(is-ch_max+64), (double(sample[is])-ped)/amp_max );
        }
      }
    }
    if( summ > 500. && summ <= 1000. )
    {
      for( int is=0; is < (int)(sample.size()-1); is++)
      {
        for( size_t iis=0; iis<ssl.size(); iis++)
        {
//           h.p1_Shape500SADC[ ssl[iis] ]->Fill( double(is), (double(sample[is])-ped)/summ );
          h.p1_Shape500SADC[ ssl[iis] ]->Fill( double(is), (double(sample[is])-ped)/amp_max );
          h.p1_ShapeNM500SADC[ ssl[iis] ]->Fill( double(is-ch_max+64), (double(sample[is])-ped)/amp_max );
        }
      }
    }
    if( summ > 1000. )
    {
      for( int is=0; is < (int)(sample.size()-1); is++)
      {
        for( size_t iis=0; iis<ssl.size(); iis++)
        {
 //          h.p1_Shape1000SADC[ ssl[iis] ]->Fill( double(is), (double(sample[is])-ped)/summ );
          h.p1_Shape1000SADC[ ssl[iis] ]->Fill( double(is), (double(sample[is])-ped)/amp_max );
          h.p1_ShapeNM1000SADC[ ssl[iis] ]->Fill( double(is-ch_max+64), (double(sample[is])-ped)/amp_max );
        }
      }
    }

  }

  double TCS_T0 = 40.; // maximum of TCS phase distr.
  double tcs_cor = CsEvent::Instance()->getTCSPhaseTime()-TCS_T0;
  for( size_t it=0; it < signals_.size(); it++ )
  {
    int ic = signals_[it].GetCellIdx();
    vector < int >ssl;
    ssl.push_back(0);
    for( size_t is = 0; is < subsets_list2store_sadc_histo_.size(); is++ )
    {
      int isub = subsets_list2store_sadc_histo_[is];
      if(  GetSubSets()[isub].IsMyCell(ic) ) ssl.push_back(is+1);
    }

    if(signals_[it].GetTimeErr() < 40.  )
    {
      for( size_t iis=0; iis<ssl.size(); iis++)
      {
        h.h2_TimeSADCvsE[ ssl[iis] ]->Fill( signals_[it].GetTime(), signals_[it].GetEnergy() );
      }
      if( signals_[it].GetEnergy() > 1. )
      {
        for( size_t iis=0; iis<ssl.size(); iis++)
        {
          h.h1_TimeSADC[ ssl[iis]]->Fill( signals_[it].GetTime() );
          h.h2_TimeSADCvsTCS[ ssl[iis]]->Fill( signals_[it].GetTime() - tcs_cor, tcs_cor );
        }
      }
    }
  }

  dir_save->cd();
}

////////////////////////////////////////////////////////////////////////////////

void CsCalorimeter::TestHistoRawInfo(void)
{
//   bool debug = true;
//   if( debug ) cout << " CsCalorimeter::TestHistoRawInfo " << GetName() << " debug " << endl;
//   cout << " The program DecodingTestSADCMore breaks with seg.viol. in HCALs " << endl;
//   DecodingTestSADCMore();
}

////////////////////////////////////////////////////////////////////////////////
