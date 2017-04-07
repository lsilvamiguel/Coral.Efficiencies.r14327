///////////////////////////////////////////////////////////////////////////////////////////////
/////  Implementation DevMonitorFEM
///////////////////////////////////////////////////////////////////////////////////////////////

#include "Reco/CellDataRaw.h"
#include "Reco/StatInfo.h"

#include "CsOpt.h"
#include "CsCalorimeter.h"
#include "CsDigitizerSADC.h"

#include "DevMonitorFEM.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////

DevMonitorFEM::DevMonitorFEM( void )
{
  name_ ="";
  fem_is_configured_ = false; 
  nFEMs_ = 0;
//  fem_norm_id_ = 0; 
  fem_range_min_ = 0.; 
  fem_range_max_ = 0.;
  c_ = NULL;
  Clear();
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::Clear( void )
{
  if( !fem_is_configured_ ) return;
//  Clear fem_signal
//  fem_measured_ = false;
  for( unsigned i=0; i< fem_signal_.size(); i++ )
  {
    fem_signal_[i]=std::pair<bool,double>(false,0.);
  }
  fem_is_good4correction_=false;   // ?? Actually it is a final decision after all processing. This ia a result after decoding and led raw data processing. Should we recover raw signals??
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::Init( CsCalorimeter *c, const std::string &name )
{
  bool debug = false;
  if( c == NULL ) return;
  c_ = c;
  if( debug ) cout << " After Init fem_norm_id_ size " << fem_norm_id_.size() << endl;
  string value;
  double fem_ref_default = 500.;

// Setting the name while it is already second parameter of Init function. Not necessary complicated if not buggy.
// Correct FEM_NAME is mandatory in case it is used to set proper options in the options files But this is not the case?
  if( CsOpt::Instance()->getOpt( c_->GetName(), "FEM_NAME", value ) )
  {
    if( debug ) cout <<" DevMonitorFEM::Init Set FEM_NAME " <<  value  << endl;
    SetName(value);
    value.clear();
  }
// So this should e always true once we created DevMonitorFEM object. We need to get rid of this useless flag
  fem_is_configured_ = true;

  std::list<int> femnorm_ids;
  if( CsOpt::Instance()->getOpt( c_->GetName(), "FEM_NORM_IDs:", femnorm_ids ) )
  {
    for( std::list<int>::const_iterator it = femnorm_ids.begin(); it!=femnorm_ids.end(); it++ )
    {
      fem_norm_id_.push_back(*it);
    }      
  }
  else
  {
    if( CsOpt::Instance()->getOpt( c_->GetName(), "FEM_NORM_ID", value ) )
    {
      int fem_norm_id = atoi(value.c_str());
      fem_norm_id_.push_back(fem_norm_id);
      value.clear();
    }
  }


  if( CsOpt::Instance()->getOpt( c_->GetName(), "FEM_AMOUNT", value ) )
  {
    nFEMs_ = atoi(value.c_str());
    value.clear();
  }

// This settings is not optimal in case of several fems for normalization Something to think about.
  if( CsOpt::Instance()->getOpt( c_->GetName(), "FEM_RANGE_MIN", value ) )
  {
    fem_range_min_ = atof(value.c_str());
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( c_->GetName(), "FEM_REF_DEFAULT", value ) )
  {
    fem_ref_default = atof(value.c_str());
    value.clear();
  }
  if( CsOpt::Instance()->getOpt( c_->GetName(), "FEM_RANGE_MAX", value ) )
  {
    fem_range_max_ = atof(value.c_str());
    value.clear();
  }
// So options are known.
// Using known information we trying to set-up FEM configuration covering
// all years before&including 2012
// As we are already inside FEM object we assume it was created on purpose 
// and hence that FEM correctins are switched on.
  if( c_->GetName() == "EC01P1__" )
  {
    ConfigureECAL1FEM();
//      if( fem_norm_id_.size() == 1 ) // Only one FEM => One laser in configuration flashing all cells. No need in source object.
//      {
//         sources_.push_back(NULL);        
//         for( size_t i=0; i< nFEMs_; i++ )
//         {
//           map_fem2source_.push_back(0);
//         }
//      }
//      else if( fem_norm_id_.size() == 2 ) //  ECAL1 2012 Lasers configuration only. That case will be cross-checked later with mapping files.
//      {
//         SourceLaserLED *s = new SourceLaserLED( c_,"ECAL1_ShashlykLaser");
//         s->Initialize();
//         sources_.push_back(s);
//         s = new SourceLaserLED( c_,"ECAL1_MainLaser");
//         s->Initialize();
//         sources_.push_back(s);
//         // We shall use fem_norm_id_  vector sorted in the corresponding sources order.  This is done to minimize configuration parameters. So re-shufle the order.
//         // We need probably better solution. ...  
//         for( size_t i=0; i< nFEMs_; i++ )
//         {
//           if( i > 7 && i <=10 ) 
//             map_fem2source_.push_back(0);
//           else if ( i >= 0 && i <= 7 )
//             map_fem2source_.push_back(1);
//         }
//         cout << " map_fem2source_.size() = " << map_fem2source_.size() << endl;
//         std::vector <int> fem_norm_id_sv = fem_norm_id_;
// 
// // Reseting the meaning
//         for( size_t i=0; i< fem_norm_id_sv.size(); i++ )
//         {
//           fem_norm_id_[i]= -1;  // Reset the meaning
//         }
// // Setting up in right order. Exploating some extra knowledge which probabbly better to configure. Ok for start-up.
//         for( size_t i=0; i< fem_norm_id_sv.size(); i++ )
//         {
//           int isource =  map_fem2source_ [fem_norm_id_sv[i]];
//           fem_norm_id_[ isource] = fem_norm_id_sv[i];
//         }
// // Checkingng the result
//         for( size_t i=0; i< fem_norm_id_sv.size(); i++ )
//         {
//           if( fem_norm_id_[i]== -1)
//           {
//             cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" Duplicate and/or NaN assignment of normalization FEM Bad source = " << i <<" Exit " << endl;
//             exit(1);
//           }
//         }
//         if( fem_norm_id_.size() != sources_.size() )
//         {
//            cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<
//                                                               " Not cosistent Nsources = " << sources_.size() <<" NormFems = " << fem_norm_id_.size() << endl;
//            exit(1);
//         } 
//         InitMapCell2FEM();   
//      }
//      else if( fem_norm_id_.size() == 0 )  // We assume that we must define FEM for normalization if FEM obect is created. It is probably too strict. We shall see later.
//      {
//         cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" No FEMs defined for normalization n = " << fem_norm_id_.size() <<" Exit " << endl;
//         exit(1);
//      }
//      else
//      {
//         cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" Defined too many FEMs for corrections n = " << fem_norm_id_.size() <<" Exit " << endl;
//         exit(1);
//      }
  }
  else if( c_->GetName() == "EC00P1__" )
  {
    ConfigureECAL0FEM();
  }
  else
  {
     if( fem_norm_id_.size() == 1 ) // Only one FEM => One laser in configuration flashing all cells. No need in source object.
     {
        sources_.push_back(NULL);
        for( size_t i=0; i< nFEMs_; i++ )
        {
          map_fem2source_.push_back(0);
        }
     }
     else if( fem_norm_id_.size() == 0 )  // We assume that we must define FEM for normalization if FEM obect is created. It is probably too strict. We shall see later.
     {
        cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" No FEMs defined for normalization n = " << fem_norm_id_.size() <<" Exit " << endl;
        exit(1);
     }
     else  
     {
        cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" Defined too many FEMs for normalization n = " << fem_norm_id_.size() <<" Exit " << endl;
        exit(1);
     }
  }  

// Consistency check
  if( fem_is_configured_ )
  {
    if( c_ == NULL )
    {
      cerr <<" FEM configuration fatal error CsCalorimeter was not specified " << endl;
      exit(1);
    } 

    if( nFEMs_ <= 0 )
    {
      cerr <<" FEM configuration fatal error no FEMs in " << c_->GetName() << endl;
      exit(1);
    }
    
    if( fem_signal_.size() > 0 )
    {
      cerr <<" FEM configuration fatal error: FEMs MUST be configured in DevMonitorFEM::Init function for " << c_->GetName() << endl;
      exit(1);
    }
    InitMemory(fem_ref_default);
  }


}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::ConfigureECAL1FEM( void )
{
     if( fem_norm_id_.size() == 1 ) // Only one FEM => One laser in configuration flashing all cells. No need in source object.
     {
        sources_.push_back(NULL);        
        for( size_t i=0; i< nFEMs_; i++ )
        {
          map_fem2source_.push_back(0);
        }
     }
     else if( fem_norm_id_.size() == 2 ) //  ECAL1 2012 Lasers configuration only. That case will be cross-checked later with mapping files.
     {
        SourceLaserLED *s = new SourceLaserLED( c_,"ECAL1_ShashlykLaser");
        s->Initialize();
        sources_.push_back(s);
        s = new SourceLaserLED( c_,"ECAL1_MainLaser");
        s->Initialize();
        sources_.push_back(s);
        // We shall use fem_norm_id_  vector sorted in the corresponding sources order.  This is done to minimize configuration parameters. So re-shufle the order.
        // We need probably better solution. ...  
        for( size_t i=0; i< nFEMs_; i++ )
        {
          if( i > 7 && i <=10 ) 
            map_fem2source_.push_back(0);
          else if ( i >= 0 && i <= 7 )
            map_fem2source_.push_back(1);
        }
        cout << " map_fem2source_.size() = " << map_fem2source_.size() << endl;
        std::vector <int> fem_norm_id_sv = fem_norm_id_;

// Reseting the meaning
        for( size_t i=0; i< fem_norm_id_sv.size(); i++ )
        {
          fem_norm_id_[i]= -1;  // Reset the meaning
        }
// Setting up in right order. Exploating some extra knowledge which probabbly better to configure. Ok for start-up.
        for( size_t i=0; i< fem_norm_id_sv.size(); i++ )
        {
          int isource =  map_fem2source_ [fem_norm_id_sv[i]];
          fem_norm_id_[ isource] = fem_norm_id_sv[i];
        }
// Checkingng the result
        for( size_t i=0; i< fem_norm_id_sv.size(); i++ )
        {
          if( fem_norm_id_[i]== -1)
          {
            cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" Duplicate and/or NaN assignment of normalization FEM Bad source = " << i <<" Exit " << endl;
            exit(1);
          }
        }
        if( fem_norm_id_.size() != sources_.size() )
        {
           cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<
                                                              " Not cosistent Nsources = " << sources_.size() <<" NormFems = " << fem_norm_id_.size() << endl;
           exit(1);
        } 
        InitMapCell2FEM();   
     }
     else if( fem_norm_id_.size() == 0 )  // We assume that we must define FEM for normalization if FEM obect is created. It is probably too strict. We shall see later.
     {
        cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" No FEMs defined for normalization n = " << fem_norm_id_.size() <<" Exit " << endl;
        exit(1);
     }
     else
     {
        cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" Defined too many FEMs for corrections n = " << fem_norm_id_.size() <<" Exit " << endl;
        exit(1);
     }
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::ConfigureECAL0FEM( void )
{
   bool debug = false;
   if( debug ) cout <<" DevMonitorFEM::ConfigureECAL0FEM debug We should go step by step " << endl;
   if( fem_norm_id_.size() == 1 ) // Only one FEM => One laser in configuration flashing all cells. No need in source object.
   {
     sources_.push_back(NULL);        
     for( size_t i=0; i< nFEMs_; i++ )
     {
       map_fem2source_.push_back(0);
     }
   }
   else if( fem_norm_id_.size() == 4 )
   {
//      cerr <<"  DevMonitorFEM::ConfigureECAL0FEM: ...Well You probably first need to implememnt FemsMapping to ECAL0 cells ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" Defined too many FEMs for corrections n = " << fem_norm_id_.size() <<" Exit " << endl;
//      exit(1);
        SourceLaserLED *s = new SourceLaserLED( c_,"ECAL0_LED0");
        s->Initialize();
        sources_.push_back(s);
        s = new SourceLaserLED( c_,"ECAL0_LED1");
        s->Initialize();
        sources_.push_back(s);
        s = new SourceLaserLED( c_,"ECAL0_LED2");
        s->Initialize();
        sources_.push_back(s);
        s = new SourceLaserLED( c_,"ECAL0_LED3");
        s->Initialize();
        sources_.push_back(s);
        for( size_t i=0; i< nFEMs_; i++ )
        {
           map_fem2source_.push_back(i);
        }

        for( size_t ic=0; ic < c_->NCells(); ic++ )
        {
          for( size_t ifm=0; ifm< fem_norm_id_.size(); ifm++ )
          {
            if( sources_[ifm]->IsMyCell(ic) )
               map_cell2fem_[ic]=fem_norm_id_[ifm];
          }
        }
   }
   else if( fem_norm_id_.size() == 0 )  // We assume that we must define FEM for normalization if FEM obect is created. It is probably too strict. We shall see later.
   {
     cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" No FEMs defined for normalization n = " << fem_norm_id_.size() <<" Exit " << endl;
     exit(1);
   }
   else  
   {
     cerr <<" CONFIGURATION FATAL ERROR ! DevMonitorFEM::Init " <<GetName() <<" for " << c_->GetName() <<" No FEMs defined for normalization n = " << fem_norm_id_.size() <<" Exit " << endl;
     exit(1);
   }
}

////////////////////////////////////////////////////////////////////////////////

int DevMonitorFEM::GetFEMIndexForCellNormalization( int cell_id ) const
{
  if( fem_norm_id_.size() == 0 )
  {
    return -1;                                 // None is setup for monitoring
  }  
  else if( fem_norm_id_.size() == 1 )
  { 
    return fem_norm_id_[0];             // No choice, only one fem is specified
  }  
  else if( fem_norm_id_.size() == 2 )
  {
    if(  map_cell2fem_.empty() )
    {
       cerr <<" DevMonitorFEM::GetFEMIndexForCellNormalization : FATAL !We haveambiguity configuration problem: " <<
                                         fem_norm_id_.size() <<" FEMs declared for normalization but cells configuration is not defined " << endl;
       exit(1);                                  
    }  
    std::map<int, int>::const_iterator it;
    if( (it= map_cell2fem_.find(cell_id))==map_cell2fem_.end() )
     return -1;
     
    return it->second; 
  }
  else if( fem_norm_id_.size() == 4 )
  {
    if(  map_cell2fem_.empty() )
    {
       cerr <<" DevMonitorFEM::GetFEMIndexForCellNormalization : FATAL !We haveambiguity configuration problem: " <<
                                         fem_norm_id_.size() <<" FEMs declared for normalization but cells configuration is not defined " << endl;
       exit(1);                                  
    }  
    std::map<int, int>::const_iterator it;
    if( (it= map_cell2fem_.find(cell_id))==map_cell2fem_.end() )
    {
      cerr <<" Normalization FEM for cell " << cell_id << " Not found " << endl;
      return -1;
    } 
//    cout <<" Normalization FEM for cell " << cell_id << " is " << it->second << endl;
    return it->second; 
  }
  else 
  {
       cerr <<" DevMonitorFEM::GetFEMIndexForCellNormalization : FATAL !We have ambiguity configuration problem: " <<
                                         fem_norm_id_.size() <<" FEMs declared for normalization but cells configuration is not defined " << endl;
       exit(1);                                  
  }
}

////////////////////////////////////////////////////////////////////////////////

std::pair<bool,double>  DevMonitorFEM::GetCorrectionFactor( int cell_id ) const
{
//   bool debug = true;
  bool debug = false;
//  if( c_->GetName() == "EC01P1__") debug = true;
  std::pair < bool,double> fem_corr (false,1.);
  if( fem_is_configured_ )
  {
    if( debug ) cout << " DevMonitorFEM::GetCorrectionFactor Now try to investigate in " << c_->GetName() <<" how to correct_leds_by_fem_signal " << endl;
//    cout <<" DevMonitorFEM::GetCorrectionFactor( " << cell_id <<") is dummy for the moment. Please have look fem_corr = " << fem_corr << endl;
    int  fem_norm_id = GetFEMIndexForCellNormalization( cell_id );
//    cout <<" DevMonitorFEM::GetCorrectionFactor( " << cell_id <<") is dummy for the moment. Please have look fem_norm_id = " << fem_norm_id << endl;
    if( fem_norm_id < 0 ) return fem_corr;
    if( debug ) cout <<" " <<c_->GetName() <<" fem_measured_ " << fem_signal_[fem_norm_id].first  << endl;
    if( fem_signal_[fem_norm_id].first )
    {
      int iref = fem_norm_id;
      if( debug ) cout <<" Fem ref " << iref << " Ref value " << fem_ref_[iref].GetMean() << " Actual " << fem_signal_[iref].second << endl;
      if((int)fem_signal_.size() > iref )
      {
        double fem_norm = fem_ref_[iref].GetMean();
        if( fem_norm > fem_range_min_ &&  fem_norm < fem_range_max_ && fem_norm > 0. )
        {
          fem_corr =std::pair < bool,double> ( true, fem_signal_[iref].second/fem_norm);
        }
        else
        {
          cerr << " WARNING !!! " << GetName() << " correct_leds_by_fem_signal bad range ref = " << fem_norm << endl;
        }
      }
      else
      {
        cerr << " correct_leds_by_fem_signal internal error iref = " << iref <<
                                 " but fem_signal_.size() " << fem_signal_.size() << endl;
        exit(1);                         
      }

    }
    if( debug )
    { 
      cout << " fem_corr = " << fem_corr.second;
      if( fem_corr.first )
        cout <<" OK " << endl;
      else  
        cout <<" NOT OK " << endl;
    }  
  }
  return fem_corr;
}

////////////////////////////////////////////////////////////////////////////////

int DevMonitorFEM::InputFEMInfo( size_t when, const string &s)
{
//   bool debug = true;
  bool debug = false;
  if( debug )
  { 
    cout << " DevMonitorFEM::InputFEMInfo " << c_->GetName() << " debug " << endl;
//    cout << s;
  }
  istringstream is(s.c_str());
  char calorim_name[132],dummy[132];
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  sscanf(str.c_str(),"%s %s ",calorim_name,dummy);  
  getline(is,str);

  while( getline(is,str) )
  {
    int ifem;
    float fm,sfm,efm,fmn,sfmn,efmn;

    if( debug ) cout << str;
    
    sscanf(str.c_str()," %d %g %g %g %g %g %g \n", &ifem, &fm, &sfm,  &efm,  &fmn,  &sfmn, &efmn );
    if( ifem < 0 )
    {
       cerr << " Error in DevMonitorFEM::InputFEMInfo " << GetName() << " bad FEM index " << ifem << endl;
       continue;
    }
    if( ifem >=  (int)fem_old_.size() || ifem >=  (int)fem_ref_.size() )
    {
       if( fem_old_.size() == 0 || fem_ref_.size() == 0 )
       {
         cerr << " Error in DevMonitorFEM::InputFEMInfo " << GetName() << " Array for FEM reading  was not initialsed " << endl;
         return 1;
       }

       cerr << " Error in DevMonitorFEM::InputFEMInfo " << GetName() << " bad FEM index " << ifem << endl;
       continue;
    }
    
    
    
    if( when == Reco::Calorimeter::OLD )
    {
      if( debug ) cout << " ifem " <<  ifem  << " set to OLD  entries " << efm << " mean " << fm << " sigma " << sfm << endl;
      fem_old_[ifem] = Reco::StatInfo( efm, fm, sfm );
    }
    else if( when == Reco::Calorimeter::PRIM )
    {
      if( debug ) cout << " ifem " <<  ifem  << " set to PRIM  entries " << efm << " mean " << fm << " sigma " << sfm << endl;
      fem_ref_[ifem] = Reco::StatInfo( efm, fm, sfm );
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int DevMonitorFEM::OutputFEMInfo( string &s ) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) 
  {
    cout << " OtputFEMInfo for  " << c_->GetName() << " start with debug " << endl;
  }
  char o[500000];
  o[0]=0;
  sprintf(o,"### \n");
  sprintf(o+strlen(o),"FEMs Calorimeter: %s N_FEMs=%zu \n",GetName().c_str(),fem_new_.size());
  sprintf(o+strlen(o)," N_FEM  Mean   Sigma    Stat    MeanN   SigmaN  StatN  \n");

  for( size_t ifem=0; ifem!=fem_new_.size(); ifem++ )
  {
    double fm = fem_new_[ifem].GetMean();
    double sfm = fem_new_[ifem].GetSigma();
    double efm = fem_new_[ifem].GetEntries();
    double fmn = fem_norm_new_[ifem].GetMean();
    double sfmn = fem_norm_new_[ifem].GetSigma();
    double efmn = fem_norm_new_[ifem].GetEntries();
    sprintf(o+strlen(o)," %zu  %7.2f  %7.2f %7.2f  %7.3f  %7.4f %7.2f\n",ifem,fm,sfm,efm,fmn,sfmn,efmn);
  }
  string ss(o);
  s = ss;
  if( debug ) cout << s;
  
//   if( debug ) 
//   {
//     cout << " OtputFEMInfo for  " << GetName() << " stop for debug " << endl;
//     exit(0);
//   }
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////

int DevMonitorFEM::OutputFEMInfoInSpills( string &s ) const
{

  bool debug = false;
//   bool debug = true;
  if( debug )
  {
    cout << " DevMonitorFEM::OutputFEMInfoInSpills " << GetName() << " debug " << endl;
  }
  unsigned NFEM = nFEMs_;

  if( debug )
  {
    cout <<" NFEM " << NFEM <<  " NledEv " << led_events_id_.size() << " FemStoreSize " << store_femn_all_[0].size() << endl;
  }

  if( store_femn_all_[0].size() == 0 ) return -1;
  if( led_events_id_.size() == 0 ) return -1;
  if( store_femn_all_[0].size() != led_events_id_.size() ) 
  {
    cerr << " ERROR in DevMonitorFEM::OutputFEMInfoInSpills internal problem " << GetName() <<
                  " store_femn_all_[0].size() = " << store_femn_all_[0].size() <<
                                  " led_events_id_.size() = " << led_events_id_.size() << endl;
    exit(1);                              
  }
  
  double ndata[led_events_id_.size()];
  double ngood[led_events_id_.size()];
  int run_min = 1000000;
  int run_max = -1000000;
  for( unsigned il=0; il < led_events_id_.size(); il++ ) 
  {
    ndata[il] = 0.;
    ngood[il] = 0.;
    if( led_events_id_[il].GetRunNumber() < run_min ) run_min = led_events_id_[il].GetRunNumber();
    if( led_events_id_[il].GetRunNumber() > run_max ) run_max = led_events_id_[il].GetRunNumber();
  }
  int nruns = run_max - run_min + 1 ;
  if( nruns <= 0 || nruns >100000 )
  {
    cerr << " ERROR in DevMonitorFEM::OutputFEMInfoInSpills internal problem " << GetName() <<" nrun_min = " << run_min <<
                                                                          " nrun_max = " << run_max << endl;
    exit(1);                                                                      
  }
  vector < Reco::StatInfo > femn_in_spill[nruns][NFEM];
  int spill_min[nruns];
  int spill_max[nruns];
  for( int irun=0 ; irun < nruns; irun++)
  {
    spill_min[irun] = 1000000;
    spill_max[irun] = -1000000;
  }
  
  for( unsigned il=0; il < led_events_id_.size(); il++ ) 
  {
    int run = led_events_id_[il].GetRunNumber();
    int irun = run-run_min;
    if( led_events_id_[il].GetBurstNumber() < spill_min[irun] ) spill_min[irun] = led_events_id_[il].GetBurstNumber();
    if( led_events_id_[il].GetBurstNumber() > spill_max[irun] ) spill_max[irun] = led_events_id_[il].GetBurstNumber();
    if( debug ) cout << " Spill " << led_events_id_[il].GetBurstNumber() <<
                 " Event in spill " << led_events_id_[il].GetEventNumberInBurst() << endl;
  }

  for( int irun=0 ; irun < nruns; irun++)
  {
    for( int ifm=0 ; ifm < (int)NFEM; ifm++ )
    {
      for( int isp=0 ; isp < spill_max[irun]; isp++)
      {
        femn_in_spill[irun][ifm].push_back( Reco::StatInfo() );
      }
    }
  }

  vector< pair <int,int> > event_min_max[nruns];
  for( int irun=0 ; irun < nruns; irun++)
  {
    for( int isp=0 ; isp < spill_max[irun]; isp++)
    {
      event_min_max[irun].push_back( pair <int,int>( 10000000,-10000000 ) );
    }
  }

  for( unsigned il=0; il < led_events_id_.size(); il++ ) 
  {
    int run = led_events_id_[il].GetRunNumber();
    int irun = run-run_min;
    int spill = led_events_id_[il].GetBurstNumber();
    int isp = spill- spill_min[irun];
    pair < int, int > &min_max = event_min_max[irun][isp];
    if( debug ) cout << " event " << led_events_id_[il].GetEventNumberInBurst() <<
                                  " min " << min_max.first <<
                                  " max " << min_max.second << endl;
    if( min_max.first > led_events_id_[il].GetEventNumberInBurst() ) min_max.first = led_events_id_[il].GetEventNumberInBurst();
    if( min_max.second < led_events_id_[il].GetEventNumberInBurst() ) min_max.second = led_events_id_[il].GetEventNumberInBurst();
    if( debug ) cout << " event " << led_events_id_[il].GetEventNumberInBurst() <<
                                  " New min " << min_max.first <<
                                  " New max " << min_max.second << endl;
  }

    if( debug ) cout << " led_events_id_.size() " << led_events_id_.size() <<
                          " store_fem_all_[0].size() = " << store_fem_all_[0].size() <<
                                  " store_femn_all_.size() = " << store_femn_all_[0].size() << endl;
    Reco::StatInfo fem_no_cuts[NFEM];
    Reco::StatInfo femn_thr[NFEM];
    Reco::StatInfo fem_strict[NFEM];
    Reco::StatInfo femn_strict[NFEM];

    for( unsigned ifm=0; ifm < NFEM; ifm++ )
    {
      fem_no_cuts[ifm].Clear();
      femn_thr[ifm].Clear();
      fem_strict[ifm].Clear();
      femn_strict[ifm].Clear();
    }
    double thr_norm = 0.01;
    double fem_range_min = 30.;
    double fem_range_max = 900.;
    for( size_t ifm=0; ifm < NFEM; ifm++ )
    {
      if( debug ) cout <<" store_fem_all_ ["<< ifm <<"].size() = " << 
                     store_fem_all_[ifm].size() << " X-check" << led_events_id_.size() << endl;
      if( store_femn_all_[ifm].size() != led_events_id_.size() ) continue;
      for( unsigned il=0; il < led_events_id_.size(); il++ ) 
      {
        fem_no_cuts[ifm].Add( store_femn_all_[ifm][il] );
        if( store_femn_all_[ifm][il] > thr_norm )
        {
          femn_thr[ifm].Add( store_femn_all_[ifm][il] );
        }
        if( store_fem_all_[ifm][il] > fem_range_min && store_fem_all_[ifm][il] < fem_range_max)
        {
          fem_strict[ifm].Add( store_fem_all_[ifm][il] );
        }
      }
    }
  for( unsigned ifm=0; ifm < NFEM; ifm++ )
  {
    femn_strict[ifm].Clear();
    for( unsigned il=0; il < led_events_id_.size(); il++ ) 
    {
      if( store_femn_all_[ifm].size() != led_events_id_.size() ) continue;
      int irun = led_events_id_[il].GetRunNumber()-run_min;
      int ispill = led_events_id_[il].GetBurstNumber() - 1;
      if( store_femn_all_[ifm][il] > thr_norm )
      {
        if( fabs( store_femn_all_[ifm][il] - femn_thr[ifm].GetMean() ) <  10.*femn_thr[ifm].GetSigma()+0.0001  )
        {
          femn_strict[ifm].Add( store_femn_all_[ifm][il] );
          femn_in_spill[irun][ifm][ispill].Add( store_femn_all_[ifm][il] );
        }
      }
    }
  }

    
  char o[5000000];
  o[0]=0;
  sprintf(o,"### \n");
  sprintf(o+strlen(o),"  OUTPUT FEM info in spills : %s NFEMS=%d \n",GetName().c_str(),NFEM);
  sprintf(o+strlen(o),"  Cell      MeanPeak   SigmaPeak  Statistic  NormFlag  CellName \n");

  for( size_t ifm=0; ifm!=NFEM; ifm++ )
  {
     if( store_femn_all_[ifm].size() != led_events_id_.size() ) continue;
     int normflag = 1;      
       sprintf(o+strlen(o)," %zu %7.3f %7.3f %7.3f %5d   %s \n",ifm,
             femn_strict[ifm].GetMean(),
             femn_strict[ifm].GetSigma(),
             femn_strict[ifm].GetEntries(),
             normflag,
	     "femname" );
      for( int irun=0 ; irun < nruns; irun++)
      {
//        cout <<" run " << run_min+irun << endl;
        sprintf(o+strlen(o)," run %d spills %d \n",run_min+irun, spill_max[irun]);
        int icnt20 = 0;
        for( int isp=0 ; isp < spill_max[irun]; isp++)
        {
          if( femn_in_spill[irun][ifm][isp].GetEntries() > 0 )
          {
        
            double ddd = 1000.*(femn_in_spill[irun][ifm][isp].GetMean()/femn_thr[ifm].GetMean()-1.);
//            cout <<" " << int(ddd);
            sprintf(o+strlen(o)," %d",int(ddd));
          }
          else
          {
//            cout <<" nan " << endl;
//             sprintf(o+strlen(o)," %d",-1001);
            sprintf(o+strlen(o)," nan ");
          }
          icnt20++;
          if( icnt20 == 20 )
          {
            icnt20 = 0;
            sprintf(o+strlen(o),"\n");
//            cout << endl;   
          }
        }
        if( icnt20 != 0 ) sprintf(o+strlen(o),"\n");
//        cout << endl;
      }
  }
  string ss(o);
//  s += ss;
  s = ss;
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::MonitorInSpill( void )
{
  bool debug = false;
//   if( c_->GetName() == "EC01P1__")   debug = true;
  if( !fem_is_configured_ ) return;
  if( !c_->GetOptions().store_leds_all  ) return;
  led_events_id_.push_back(  c_->GetEventID() );
  for( int ifm =0; ifm < (int)fem_signal_.size(); ifm++ )
  {
    if( fem_signal_[ifm].first )
      store_fem_all_[ifm].push_back(fem_signal_[ifm].second);
    else   
      store_fem_all_[ifm].push_back(-1. );                            // We should skip this -1 at the end carefully
    int is = map_fem2source_[ifm];
    int fem_norm_id = fem_norm_id_[is];
    if( fem_signal_[fem_norm_id].first )
    {
      double fem_norm = fem_signal_[fem_norm_id].second;
      if(  fem_norm < 0.001 || fem_norm > 100000. )
      {
        cerr <<" Seems internal error ? DevMonitorFEM::MonitorInSpil  Fem signal used for normalization has inproper value ." << fem_norm << "  .. Let it be fatal.. We'll see .. " <<endl;
        exit(1);
      }
      if( fem_signal_[ifm].first )
        store_femn_all_[ifm].push_back(fem_signal_[ifm].second/fem_norm);
      else   
        store_femn_all_[ifm].push_back(-1. );                            // We should skip this -1 at the end carefully
    }
    else
    {
      store_femn_all_[ifm].push_back(-1. );
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::InitReadOut( void )
{
  bool debug = true;
  if(  c_->GetName() == "EC01P1__" )
  {
    if( debug ) cout <<" DevMonitorFEM::InitReadOut " << c_->GetName() << " nFEMs_ = " << nFEMs_ << endl;
    for( int nfem=0; nfem < (int)nFEMs_; nfem++ )
    {
      CsDigitizerSADC::Type type_fem_sadc = CsDigitizerSADC::SADC;
      CsDigitizerSADC::Shape fem_table = CsDigitizerSADC::LED_ECAL1_GAMS;
      CsDigitizerSADC *sd = new CsDigitizerSADC( c_ , CsDigitizerSADC::MaxSimple, type_fem_sadc, fem_table );
      dig_max_fem_.push_back(sd);
      vector < double> vpar;
      vpar.push_back(1.);
      double fem_ped_min = 0.;
      double fem_ped_max = 2.;
      vpar.push_back(fem_ped_min);
      vpar.push_back(fem_ped_max);
      bool ok = dig_max_fem_.back()->LoadSettings(vpar);
      if( !ok) cerr << " Some problem in loading SADC settings to FEM ?? " << c_->GetName() << endl;
    }
  }
  else if(  c_->GetName() == "HC01P1__" )
  {
  // Clearly not correct settings for HCAL1 FEM. Presently I have no idea abot the signal shape.
  //  Anyway we are using methods ignoring any shape dependence.
    for( int nfem=0; nfem < (int)nFEMs_; nfem++ )
    {
      CsDigitizerSADC::Type type_fem_sadc = CsDigitizerSADC::SADC;
      CsDigitizerSADC::Shape fem_table = CsDigitizerSADC::LED_ECAL1_GAMS;
      CsDigitizerSADC *sd = new CsDigitizerSADC( c_ , CsDigitizerSADC::MaxSimple, type_fem_sadc, fem_table );
      dig_max_fem_.push_back(sd);
      vector < double> vpar;
      vpar.push_back(1.);
      double fem_ped_min = 0.;
      double fem_ped_max = 2.;
      vpar.push_back(fem_ped_min);
      vpar.push_back(fem_ped_max);
      bool ok = dig_max_fem_.back()->LoadSettings(vpar);
      if( !ok) cerr << " Some problem in loading SADC settings to FEM ?? " << c_->GetName() << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::InitMapCell2FEM( void )
{
//   bool debug = true;
  bool debug = false;
  if( fem_norm_id_.size() <=1 ) return;
  if( c_ != NULL )
  {
    if( c_ ->GetName() == "EC01P1__" )
    {
// Configuration consistency check must be done elsewere
      if( fem_norm_id_.size() != 2 ) return;
      if( sources_.size() != 2 ) return;       
      for( size_t ic=0; ic < c_->NCells(); ic++ )
      {
          if( sources_[0]->IsMyCell(ic) )
             map_cell2fem_[ic]=fem_norm_id_[0];
          else if( sources_[1]->IsMyCell(ic) ) 
            map_cell2fem_[ic]=fem_norm_id_[1];
          else
          {
            cerr <<" Warning !!!  DevMonitorFEM::InitMapCell2FEM " << GetName() << " For " << c_ ->GetName() <<
                                                " Cell " << ic << " "  << c_->GetCellName(ic)  << " Does not assigned to any LED/Laser source Will not be corrected !!! " << endl;
          }  
       } // seting up map elements
    }  // EC01P1__ case
  } 
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::InitMemory( double fem_ref_default )
{
  for( unsigned nfem=0; nfem < nFEMs_; nfem++ )
  {
    fem_signal_.push_back(  std::pair<bool,double>( false,0.) );
    fem_norm_new_.push_back(Reco::StatInfo());
    fem_new_.push_back(Reco::StatInfo());
    fem_old_.push_back(Reco::StatInfo());
    fem_ref_.push_back(Reco::StatInfo(1.,fem_ref_default,1.));
  }
  assert( fem_signal_.size() == nFEMs_);
  assert( fem_norm_new_.size() == nFEMs_);
  assert( fem_new_.size() == nFEMs_);
  assert( fem_old_.size() == nFEMs_);
  assert( fem_ref_.size() == nFEMs_);
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::DecodeChipDigitsECAL0( const CS::Chip::Digits& digits )
{
//   bool debug = true;
  bool debug = false;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  std::pair<m_it,m_it> m_range;

  {
    if( debug ) cout <<" DevMonitorFEM::DecodeChipDigits My name is " <<   GetName() << endl;

    std::string pindiodname( GetName().c_str()); // "EC00FEM_"
    m_range = digits.equal_range( pindiodname );

    for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
    {
      if( debug ) cout << " decode PIN diod digit " << endl;
      const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
      const CS::ChipSADC::DigitCalib *dsc = dynamic_cast<const CS::ChipSADC::DigitCalib*>( ds );
      if( debug ) cout << " dsc = " << dsc << " ds " << ds << endl;


      int fem_index = ds->GetChannel();
      if( debug ) cout << " Channel " <<  fem_index   << endl;
      if( dsc != NULL )
      {
         const std::vector<CS::ChipSADC::DigitCalib::Box> &boxes = dsc->GetBoxes();
//          if( boxes.size() > 0 )
//          {
//             cout <<" Clearly we have some Boxes !!! " <<   boxes.size() <<" for RefDet " << dsc->GetRefDet()  << " Chan# " <<dsc->GetChannel() << endl;
//          }
      }
      else
      { 
        int32 x_digit = ds->GetX();
        int32 y_digit = ds->GetY();
        if( debug ) cout << "  DevMonitorFEM::DecodeChipDigitsECAL0 decode FEM digit x=" << x_digit << " y=" << y_digit << " Chan# " <<ds->GetChannel() <<  endl;
        if( x_digit > 1 ) continue;
        if( y_digit > 1 ) continue;
        int32 femnum = 2*y_digit+x_digit;
// Some messy swaps        
        if( femnum == 2 )
          femnum = 3;
        else if(  femnum == 3 ) 
          femnum = 2;

        if( femnum >= (int)nFEMs_ ) continue;
        if( debug ) cout <<" DevMonitorFEM::DecodeChipDigitsECAL0 changing fem_index from " << fem_index << " to " << femnum << endl;
        fem_index = femnum;
      }
      
      if( debug ) cout <<" dig_max_fem_.size() " << dig_max_fem_.size() << endl;

      if( (int)dig_max_fem_.size() <= fem_index ) continue;                       // This is a bug in configuration Break the execution?

      const vector<CS::uint16>& sample = ds->GetSamples();
//       vector<CS::uint16> sample_inv;
      if( debug ) cout <<" EC00FEM " << fem_index <<" Samples: ";
      for( int is=0; is < (int)sample.size(); is++)
      {
        if( debug ) cout <<"  " <<sample[is];
//         sample_inv.push_back(1023-sample[is]);
      }
     if( debug ) cout << endl;

//  What is the meaning  of the Fit parameters like     fem_index  and ds->GetDetID().GetName() ????

//      dig_max_fem_[fem_index]->Fit( sample_inv , fem_index, ds->GetDetID().GetName());
      dig_max_fem_[fem_index]->Fit( sample , fem_index);

      double signal_fem_max = dig_max_fem_[fem_index]->GetSignal(); 
      double time_fem_max = dig_max_fem_[fem_index]->GetTime();

      SetSignal( fem_index, signal_fem_max);
   }
// // New fem imlementation 
// // Correct laser amplitudes in calorimeter
// //    bool all_fems_ok = c_->ApplyCorrectionsMonitorFEM();
// //     if( all_fems_ok ) TurnedOutGood4Corrections();            // Set up the finall result of FEM corrections. As a matter of fact statement after all corrections.
//                                                                                            // Basically fem_finally_turned_out_to_be_good4correction_ Not perfect cross settings
//     CalculateInternalCrossNormalization();                         // Just monitoring function. Must be displaced from the decoding sequence.
  }
}
////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::DecodeChipDigitsECAL1( const CS::Chip::Digits& digits )
{
//   bool debug = true;
  bool debug = false;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  std::pair<m_it,m_it> m_range;

  {
    if( debug ) cout <<" DevMonitorFEM::DecodeChipDigits My name is " <<   GetName() << endl;
    std::string pindiodname( GetName().c_str()); // "EC01FEM"
    m_range = digits.equal_range( pindiodname );
    for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
    {
      if( debug ) cout << " decode PIN diod digit " << endl;
      const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
      const CS::ChipSADC::DigitCalib *dsc = dynamic_cast<const CS::ChipSADC::DigitCalib*>( ds );
      int fem_index = ds->GetChannel();
      if( dsc != NULL )
      {
         const std::vector<CS::ChipSADC::DigitCalib::Box> &boxes = dsc->GetBoxes();
//          if( boxes.size() > 0 )
//          {
//             cout <<" Clearly we have some Boxes !!! " <<   boxes.size() <<" for RefDet " << dsc->GetRefDet()  << " Chan# " <<dsc->GetChannel() << endl;
//          }
      }
      else
      { 
        int32 x_digit = ds->GetX();
        int32 y_digit = ds->GetY();
        if( debug ) cout << "  CsECAL1::DecodeChipDigits decode FEM digit x=" << x_digit << " y=" << y_digit << " Chan# " <<ds->GetChannel() <<  endl;
        if( x_digit != 0 ) continue;
        if( y_digit < 0 || y_digit >= (int)nFEMs_ ) continue;
        fem_index = y_digit;
      }

      if( (int)dig_max_fem_.size() <= fem_index ) continue;                       // This is a bug in configuration Break the execution?

      const vector<CS::uint16>& sample = ds->GetSamples();
      vector<CS::uint16> sample_inv;
      for( int is=0; is < (int)sample.size(); is++)
      {
        sample_inv.push_back(1023-sample[is]);
      }

//  What is the meaning  of the Fit parameters like     fem_index  and ds->GetDetID().GetName() ????

      dig_max_fem_[fem_index]->Fit( sample_inv , fem_index);

      double signal_fem_max = dig_max_fem_[fem_index]->GetSignal(); 
      double time_fem_max = dig_max_fem_[fem_index]->GetTime();

      SetSignal( fem_index, signal_fem_max);
   }
// New fem imlementation 
// Correct laser amplitudes in calorimeter
//    bool all_fems_ok = c_->ApplyCorrectionsMonitorFEM();
//     if( all_fems_ok ) TurnedOutGood4Corrections();            // Set up the finall result of FEM corrections. As a matter of fact statement after all corrections.
                                                                                           // Basically fem_finally_turned_out_to_be_good4correction_ Not perfect cross settings
    CalculateInternalCrossNormalization();                         // Just monitoring function. Must be displaced from the decoding sequence.
  }
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::DecodeChipDigitsHCAL1( const CS::Chip::Digits& digits )
{
//   bool debug = true;
  bool debug = false;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  std::pair<m_it,m_it> m_range;

  {
    std::string pindiodname( c_->GetName().c_str()); // 
    m_range = digits.equal_range( pindiodname );
    for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
    {
//     if( debug ) cout << " decode PIN diod digit " << endl;
       const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit *>(d_it->second); 
       if( ds == NULL ) assert(false);
       int32 x_digit = ds->GetX();
       int32 y_digit = ds->GetY();
       if( x_digit == 27 && y_digit == 0 )
       { 
         if( dig_max_fem_.size() != 1 ) assert(false);
         const std::vector<CS::uint16>& sample = ds->GetSamples();
         dig_max_fem_[0]->Fit( sample , y_digit);
         double signal_fem_max = dig_max_fem_[0]->GetSignal();
         double time_fem_max = dig_max_fem_[0]->GetTime();
//          cout <<" This HCAL1 pindiod in LED event Try to store  it " << dig_max_fem_.size() << " signal " << signal_fem_max << endl;
         SetSignal( 0,signal_fem_max);
         break;
       } 
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::DecodeChipDigits( const CS::Chip::Digits& digits )
{
//   bool debug = true;
  bool debug = false;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  std::pair<m_it,m_it> m_range;
  if( c_->GetName() == "EC01P1__")
  {
    DecodeChipDigitsECAL1( digits );
//     if( debug ) cout <<" DevMonitorFEM::DecodeChipDigits My name is " <<   GetName() << endl;
//     std::string pindiodname( GetName().c_str()); // "EC01FEM"
//     m_range = digits.equal_range( pindiodname );
//     for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
//     {
//       if( debug ) cout << " decode PIN diod digit " << endl;
//       const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
//       const CS::ChipSADC::DigitCalib *dsc = dynamic_cast<const CS::ChipSADC::DigitCalib*>( ds );
//       int fem_index = ds->GetChannel();
//       if( dsc != NULL )
//       {
//          const std::vector<CS::ChipSADC::DigitCalib::Box> &boxes = dsc->GetBoxes();
// //          if( boxes.size() > 0 )
// //          {
// //             cout <<" Clearly we have some Boxes !!! " <<   boxes.size() <<" for RefDet " << dsc->GetRefDet()  << " Chan# " <<dsc->GetChannel() << endl;
// //          }
//       }
//       else
//       { 
//         int32 x_digit = ds->GetX();
//         int32 y_digit = ds->GetY();
//         if( debug ) cout << "  CsECAL1::DecodeChipDigits decode FEM digit x=" << x_digit << " y=" << y_digit << " Chan# " <<dsc->GetChannel() <<  endl;
//         if( x_digit != 0 ) continue;
//         if( y_digit < 0 || y_digit >= (int)nFEMs_ ) continue;
//         fem_index = y_digit;
//       }
// 
//       if( (int)dig_max_fem_.size() <= fem_index ) continue;                       // This is a bug in configuration Break the execution?
// 
//       const vector<CS::uint16>& sample = ds->GetSamples();
//       vector<CS::uint16> sample_inv;
//       for( int is=0; is < (int)sample.size(); is++)
//       {
//         sample_inv.push_back(1023-sample[is]);
//       }
// 
// //  What is the meaning  of the Fit parameters like     fem_index  and ds->GetDetID().GetName() ????
// 
//       dig_max_fem_[fem_index]->Fit( sample_inv , fem_index, ds->GetDetID().GetName());
// 
//       double signal_fem_max = dig_max_fem_[fem_index]->GetSignal(); 
//       double time_fem_max = dig_max_fem_[fem_index]->GetTime();
// 
//       SetSignal( fem_index, signal_fem_max);
//    }
// // New fem imlementation 
// // Correct laser amplitudes in calorimeter
// //    bool all_fems_ok = c_->ApplyCorrectionsMonitorFEM();
// //     if( all_fems_ok ) TurnedOutGood4Corrections();            // Set up the finall result of FEM corrections. As a matter of fact statement after all corrections.
//                                                                                            // Basically fem_finally_turned_out_to_be_good4correction_ Not perfect cross settings
//     CalculateInternalCrossNormalization();                         // Just monitoring function. Must be displaced from the decoding sequence.
  }
  else if( c_->GetName() == "HC01P1__")
  {
    DecodeChipDigitsHCAL1( digits );
//     std::string pindiodname( c_->GetName().c_str()); // 
//     m_range = digits.equal_range( pindiodname );
//     for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
//     {
// //     if( debug ) cout << " decode PIN diod digit " << endl;
//        const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit *>(d_it->second); 
//        if( ds == NULL ) assert(false);
//        int32 x_digit = ds->GetX();
//        int32 y_digit = ds->GetY();
//        if( x_digit == 27 && y_digit == 0 )
//        { 
//          if( dig_max_fem_.size() != 1 ) assert(false);
//          const std::vector<CS::uint16>& sample = ds->GetSamples();
//          dig_max_fem_[0]->Fit( sample , y_digit, ds->GetDetID().GetName());
//          double signal_fem_max = dig_max_fem_[0]->GetSignal();
//          double time_fem_max = dig_max_fem_[0]->GetTime();
// //          cout <<" This HCAL1 pindiod in LED event Try to store  it " << dig_max_fem_.size() << " signal " << signal_fem_max << endl;
//          SetSignal( 0,signal_fem_max);
//          break;
//        } 
//     }
  }
  else if( c_->GetName() == "EC00P1__")
  {
    DecodeChipDigitsECAL0( digits );
  }
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::SetSignal( int fem_index, double signal)
{
  if( fem_index <0 || fem_index >= (int)nFEMs_)
  {
    cerr <<" DevMonitorFEM::SetSignal " << c_->GetName() <<" bad fem index " << fem_index << endl;
    return;
  } 
  fem_signal_[fem_index] = std::pair<bool,double> ( true,signal);
  fem_new_[fem_index].Add( signal );
//  fem_measured_= true;
}

////////////////////////////////////////////////////////////////////////////////

void DevMonitorFEM::CalculateInternalCrossNormalization( void )
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " DevMonitorFEM::CalculateInternalCrossNormalization debug in " << GetName() <<"  fem_is_configured_ " <<  fem_is_configured_ << endl;
  if( fem_is_configured_ )
  {
     for( size_t in=0; in< fem_norm_id_.size(); in++ )
     {
       int fem_norm_id = fem_norm_id_[in];
       if( ! fem_signal_[fem_norm_id].first ) continue;  // Normalization FEM amplitude is not available in this event
 // apply fem normalisation  
      int index_FEM_to_normalise = fem_norm_id;
      if( debug ) cout << " Calculate FEM noralized statistic using " << index_FEM_to_normalise << endl;
      double fem_range_min = fem_range_min_;
      double fem_range_max = fem_range_max_;
      double fem_norm = fem_signal_[ index_FEM_to_normalise].second;
      if( !(fem_norm > fem_range_min &&  fem_norm < fem_range_max) ) continue; // defined range is not specific for this FEM
      if( !(fem_norm > 0.001 &&  fem_norm < 100000.) ) continue;                         // still some precautions
      for( int ifm =0; ifm < (int)fem_signal_.size(); ifm++ )
      {
        if( !fem_signal_[ifm].first ) continue;                                                           // FEM was not measured in this event
        if( (int)in != map_fem2source_[ifm] ) continue;                                            // Normalization only FEMs which mapped to the same source as norm FEM
        fem_norm_new_[ifm].Add(  fem_signal_[ifm].second/fem_norm );
      }
    }
    if( debug ) cout << " Calculate FEM normalized statistic OK " << endl;
  }
  if( debug ) cout << " DevMonitorFEM::CalculateInternalCrossNormalization OK " << endl;
}

////////////////////////////////////////////////////////////////////////////////

std::string   DevMonitorFEM::GetSummaryComment( void ) const
{
  if( !fem_is_configured_ ) return std::string("### No FEM corrections ");
  if( fem_ref_.size() == 0 ) return std::string("### No FEM corrections ");
//   int iref = fem_norm_id_;                  
//   if( iref < 0 || iref >= (int)fem_ref_.size() ) return std::string("### No FEM corrections ");
  if( fem_norm_id_.size() == 0 ) return std::string("### No FEM corrections ");
  char o[50000];
  o[0]=0;
  int iref = fem_norm_id_[0];
  sprintf(o,"### FEM %s configured NORM_ID = %d FEM_REF = %7.3f FEM_ACTUAL = %7.3f ",
                 GetName().c_str(),fem_norm_id_[0],fem_ref_[iref].GetMean(),fem_new_[iref].GetMean());
  for( size_t i=1; i < fem_norm_id_.size(); i++)
  {
    iref = fem_norm_id_[i];
    if( (int)fem_ref_.size() <= iref )
    {
      sprintf(o+strlen(o)," MORE FEM CONFIGURATION PROBLEMS DETECTED");
      cerr <<" fem_ref_.size() " << fem_ref_.size() <<" Not fit to fem_id " << iref << endl;
      break;
    }
    sprintf(o+strlen(o)," NORM_ID = %d FEM_REF = %7.3f FEM_ACTUAL = %7.3f",iref,fem_ref_[iref].GetMean(),fem_new_[iref].GetMean());
  }             
  string ss(o); 
  return ss;
}

////////////////////////////////////////////////////////////////////////////////

SourceLaserLED::SourceLaserLED ( CsCalorimeter *c, const std::string &name ) : name_(name),c_(c)
{
}

////////////////////////////////////////////////////////////////////////////////

void SourceLaserLED:: InitCellsMap ( void )
{
  bool debug = false;
  if( c_ != NULL )
  {
    if( c_ ->GetName() == "EC01P1__" )
    {
      if( !(GetName() == "ECAL1_ShashlykLaser" || GetName() == "ECAL1_MainLaser" ) )
      {
         cerr <<" FATAL CONFIGURATION ERROR in SourceLaserLED:: InitCellsMap Unknown source : " << GetName() << endl;
         exit(1);
      }
      for( size_t ic=0; ic < c_->NCells(); ic++ )
      {
        double rl =  c_->GetCells()[ic].GetCellType().GetRadiationLength();
        bool is_shashlyk = ( fabs(rl - 17.65) < 0.1 );
        if( debug ) cout << " Cell " << ic <<" RadLength " << rl <<" is Shashlyk " <<  is_shashlyk << endl;
        if( GetName() == "ECAL1_ShashlykLaser" )
        {
          if(is_shashlyk) 
             map_cell2source_[ic]=true;
          else   
             map_cell2source_[ic]=false;
        }  
        else if( GetName() == "ECAL1_MainLaser" )
        {   
          if( !is_shashlyk) 
            map_cell2source_[ic]=true;
          else   
             map_cell2source_[ic]=false;
        }  
      }
    }
    else if( c_ ->GetName() == "EC00P1__" ) 
    {
        size_t xmi = 0;
        size_t xma = 0;
        size_t ymi = 0;
        size_t yma = 0;
        if( GetName() == "ECAL0_LED0" )
        {
          // new box
           xmi = 3;
           xma = 14;
           ymi = 3;
           yma = 11;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
          // new box
           xmi = 0;
           xma = 2;
           ymi = 6;
           yma = 14;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
          // new box
           xmi = 3;
           xma = 17;
           ymi = 0;
           yma = 2;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
        }
        else if( GetName() == "ECAL0_LED1" )
        {
          // new box
           xmi = 18;
           xma = 29;
           ymi = 0;
           yma = 2;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
          // new box
           xmi = 15;
           xma = 29;
           ymi = 3;
           yma = 11;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
          // new box
           xmi = 30;
           xma = 32;
           ymi = 6;
           yma = 8;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
        }
        else if( GetName() == "ECAL0_LED2" )
        {
          // new box
           xmi = 30;
           xma = 32;
           ymi = 9;
           yma = 20;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
          // new box
           xmi = 18;
           xma = 29;
           ymi = 12;
           yma = 23;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
          // new box
           xmi = 21;
           xma = 29;
           ymi = 24;
           yma = 26;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
        }
        else if( GetName() == "ECAL0_LED3" )
        {
          // new box
           xmi = 3;
           xma = 28;
           ymi = 24;
           yma = 26;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
          // new box
           xmi = 3;
           xma = 17;
           ymi = 12;
           yma = 23;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
          // new box
           xmi = 0;
           xma = 2;
           ymi = 15;
           yma = 20;
           for( size_t ix=xmi; ix<xma+1; ix++ )
             for( size_t iy=ymi; iy<yma+1; iy++ )
             {
                int icell =c_-> GetCellOfColumnRow(ix,iy);
                if( icell < 0 ) continue;
                map_cell2source_[icell]=true;
             }
        }
    } // EC00P1__
  }
}

////////////////////////////////////////////////////////////////////////////////

void SourceLaserLED::Initialize ( void )
{
  InitCellsMap();
}

////////////////////////////////////////////////////////////////////////////////
    
bool SourceLaserLED::IsMyCell(  int cell_id ) const
{
  if( map_cell2source_.empty() ) return true; // By default all cells attached to the led/laser source. Or the source is flashing the space?
  std::map<int, bool>::const_iterator it;
  if ( (it= map_cell2source_.find(cell_id))==map_cell2source_.end() )
    return false;
  return it->second;
}

////////////////////////////////////////////////////////////////////////////////

