
#include <iomanip>
#include <sstream>

#include "TDirectory.h"
#include "TProfile.h"

#include "Reco/CellDataRaw.h"
#include "Reco/myROOT_utils.h"

#include "StoreLED.h"

#include "HistoInSpillLED.h"

#define CDB_LINE_MAX 132

using namespace Reco;
using namespace std;
 
const double StoreCaloEvents::spill_gate   = 22.7;
const double StoreCaloEvents::spill_tcycle = 46.792;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StoreCaloEvents:: StoreCaloEvents ( Reco::Calorimeter *c ) : c_(c),run_min_(10000000),run_max_(0),spill_min_(10000000),spill_max_(0),dummy_spill_info_(-1,-1),ncells_(0)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  StoreCaloEvents::AddEvent (  CaloEvent  * ev) 
{   
   if( ev == NULL ) return;
   size_t run = ev->event_id_.GetRunNumber();
   size_t spill = ev->event_id_.GetBurstNumber();
   std::pair< size_t, size_t > indx( run, spill);
   events_[indx].push_back(ev);

   if(run_min_ > run ) run_min_= run;
   if(run_max_ < run ) run_max_ = run;
   SetSpillsInfo( run, spill );

   vev_id_.push_back(ev->event_id_);
   if( vev_id_.size() == 1 ) ncells_ =ev->data_.size();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  StoreCaloEvents::SetSpillsInfo ( size_t run, size_t spill ) 
{   
   if(spill_min_ > spill ) spill_min_= spill;
   if( spill_max_ < spill ) spill_max_= spill;
   std::pair<  int, int> smima = GetSpillsInfo( run );
   if( smima .first == -1 )
   {
     spill_min_max_[run] = std::pair< int, int> (spill,spill);
   }
   else
   {
     if(smima.first > (int)spill ) spill_min_max_[run].first = spill;
     if(smima.second < (int)spill ) spill_min_max_[run].second = spill;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const std::pair<  int, int> & StoreCaloEvents::GetSpillsInfo( size_t run ) const
{
  std::map <  size_t, std::pair< int, int > > ::const_iterator it = spill_min_max_.find(run);
  if( it !=   spill_min_max_.end() )
    return it->second;
  else 
   return dummy_spill_info_; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const StatInfo &  StoreCaloEvents::GetValue( size_t run, size_t spill, size_t cell ) const
{
  std::pair< size_t, size_t > indx( run, spill);
  std::map <  std::pair< size_t, size_t >,  std::vector<StatInfo> > ::const_iterator it = leds_in_spill.find(indx);
  if( it !=  leds_in_spill.end() )
  {
    if( cell < it->second.size() )
      return it->second[cell];
    else  
      return zero; 
  } 
  else
   return zero; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const StatInfo &  StoreCaloEvents::GetStrictAverageValue ( size_t cell ) const
{
  if( cell < led_strict.size() )
    return led_strict[cell];
  else
    return zero;  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void StoreCaloEvents::Process ( void )
{
//   bool debug = true;
  bool debug = false;
  
  if( debug ) cout << " StoreCaloEvents::Process debug NEvents = " << vev_id_.size() << endl;
  if( debug ) cout <<" Run min " << run_min_ <<" Run max " << run_max_<< " Spill min " << spill_min_ <<  " Spill max " << spill_max_ << endl;
  if( vev_id_.size() == 0 ) return;
  size_t ncells = NCells();
  if( debug ) cout <<" ncells " << ncells << endl;
  if( ncells <= 0 )
  {
    cerr <<"StoreCaloEvents::Process Internall error ncells " << ncells << endl;
    exit(1);
  }
  if( Processed() )
  {
    cerr <<"Usage ERROR: StoreCaloEvents::Process must be done only once " << endl;
    exit(1);
  }
  int ntimeok = 0;
  int ntimeneok = 0;

  for( unsigned iev=0; iev< vev_id_.size(); iev++ )
  {
    size_t run = vev_id_[iev].GetRunNumber();
    size_t spill = vev_id_[iev].GetBurstNumber();
    std::pair< size_t, size_t > indx( run, spill);
    double time_in_spill = vev_id_[iev].GetTimeInSpill();
    if( time_in_spill < 0.00001 )
      ntimeneok++;
    else  
      ntimeok++;
    const time_t  &tt = vev_id_[iev].GetTime();
    int ev_cnt = statinspill_[indx].ev_cnt;
    if( ev_cnt == 0 )
    {
      statinspill_[indx].tmin = tt;
      statinspill_[indx].tmax = tt;
      statinspill_[indx].evmin = vev_id_[iev].GetEventNumberInBurst();
      statinspill_[indx].evmax = vev_id_[iev].GetEventNumberInBurst();
      statinspill_[indx].ev_cnt++;
    }
    else
    {
      if( statinspill_[indx].tmin > tt )
      {
        statinspill_[indx].tmin = tt;
        statinspill_[indx].evmin = vev_id_[iev].GetEventNumberInBurst();
      }
      if( statinspill_[indx].tmax < tt )
      {
        statinspill_[indx].tmax = tt;
        statinspill_[indx].evmax = vev_id_[iev].GetEventNumberInBurst();
      }
      statinspill_[indx].ev_cnt++;
    }
//    cout <<" ev_cnt  " << ev_cnt << endl;
//    if( debug ) cout << " iev " << iev <<" run " << run << " spill " << spill <<" time_in_spill " << time_in_spill << " time " << tt << endl;
  }

  for( unsigned icell=0; icell < ncells; icell++ )
  {
    led_no_cuts.push_back( StatInfo() );
    led_thr.push_back( StatInfo() );
    led_strict.push_back( StatInfo() );
  }  
  leds_in_spill.clear();

  for( unsigned irun=run_min_; irun <= run_max_; irun++ )
  {
    for( unsigned ispill=spill_min_; ispill <= spill_max_; ispill++ )
    {
       std::pair< size_t, size_t > indx( irun, ispill);
       vector<StatInfo> st( ncells );
       leds_in_spill[indx] =st;
    }
  }  

  double thr = 20.;
  for( unsigned icell=0; icell < ncells; icell++ )
  {
     for( unsigned irun=run_min_; irun <= run_max_; irun++ )
     {
       for( unsigned ispill=spill_min_; ispill <= spill_max_; ispill++ )
       {
         std::pair< size_t, size_t > indx( irun, ispill);
         std::vector<  CaloEvent  * > &v = events_[indx];
         if( v.size() == 0 ) continue;
         for( unsigned il=0; il < v.size(); il++ ) 
         {
           double amp =  v[il]->data_[icell];
           led_no_cuts[icell].Add( amp );
           if( amp > thr )
           {
             led_thr[icell].Add( amp );
           }
         }
      }
    }

     for( unsigned irun=run_min_; irun <= run_max_; irun++ )
     {
       for( unsigned ispill=spill_min_; ispill <= spill_max_; ispill++ )
       {
         std::pair< size_t, size_t > indx( irun, ispill);
         std::vector< CaloEvent  * > &v = events_[indx];
         if( v.size() == 0 ) continue;
         for( unsigned il=0; il < v.size(); il++ ) 
         {
           double amp =  v[il]->data_[icell];
           if( amp > thr && fabs( amp- led_thr[icell].GetMean() ) <  10.*led_thr[icell].GetSigma()  )
           {
             led_strict[icell].Add( amp );
             if(leds_in_spill[indx].size() != ncells )
             {
                cerr <<"StoreCaloEvents::Process Internall error " << endl;
                exit(1);
             }
             leds_in_spill[indx][icell].Add( amp );
           }
         }
      }
    }

  } // new cell

  if( debug ) cout << " StoreCaloEventse::Process debug OK " << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StoreLED::StoreLED ( Reco::Calorimeter *c )  :   c_(c), store_leds_all_(NULL), store_leds_in_spill_all_(NULL), histo_in_spill_led_(NULL)
{
//   bool debug = true; 
  bool debug = false; 
  if( c_ == NULL )
  {
    cerr <<" Error in StoreLED::StoreLED: calorimeter is NULL " << endl;
    assert( false);
  }
  if( store_leds_all_ == NULL ) store_leds_all_ = new StoreCaloEvents( c_ );
  if( store_leds_in_spill_all_ == NULL ) store_leds_in_spill_all_ = new StoreCaloEvents( c_ );
  if( debug ) cout <<" StoreLED:: StoreLED  debug " << c_->GetName() <<  "OK " << endl;
} 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void StoreLED::AddEvent ( bool off_spill )
{

  CaloEvent * led_ev = new CaloEvent( NCells(), c_->GetEventID() );
  for (vector<CellDataRaw>::const_iterator it=c_->GetSignals().begin(); it!=c_->GetSignals().end(); it++ )
  {
    size_t icell = it->GetCellIdx();
    double amp_led = it->GetAmplitude();
    led_ev->data_[icell] = amp_led;
  }
  if( off_spill )
    store_leds_all_->AddEvent( led_ev);
  else  
    store_leds_in_spill_all_->AddEvent( led_ev);

} 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HistoInSpillLED::HistoInSpillLED ( void ) : h2fitpar_good4monitoring(NULL),h2fitpar_bad4monitoring(NULL)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void StoreLED::FillHistoInSpillLED ( size_t icell, double framp, double frspill, bool allok )
{
  if(histo_in_spill_led_ == NULL )
  {
    histo_in_spill_led_ = new HistoInSpillLED();
    if( histo_in_spill_led_ == NULL ) assert(false);
    char dir_name[111];
    TDirectory *dir_save = gDirectory;
    bool ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",c_->GetName().c_str());
    TDirectory *root_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok =root_dir->cd();
    assert(ok);
    sprintf(dir_name,"%s_InSpillLED_Cells",c_->GetName().c_str());
    TDirectory *dirlis = myROOT_utils::TDirectory_checked(dir_name);
    ok = dirlis->cd();
    assert(ok);

    HistoInSpillLED &h = *histo_in_spill_led_;
    size_t nbins = 20;
    double spmin = 0.;
    double spmax = 1.;
    h.h2fitpar_good4monitoring = myROOT_utils::TH2D_checked( "h2fitpar_good4monitoring"," m1 vs m0   for  events good4monitring ",1000,-0.5,0.5,1000,-0.5,0.5);
    h.h2fitpar_bad4monitoring = myROOT_utils::TH2D_checked( "h2fitpar_bad4monitoring"," m1 vs m0   for  events bad4monitring ",1000,-0.5,0.5,1000,-0.5,0.5);
    for( size_t icell=0; icell< NCells(); icell++ )
    {
      char name[132];
      char hist_name[132];
      sprintf(name,"InSpillLED%s",c_->GetCellName(icell).c_str());
      sprintf(hist_name," LED in cell %s ",c_->GetCellName(icell).c_str());
      TProfile *p = myROOT_utils::TProfile_checked(name,hist_name,nbins,spmin,spmax);
      h.p1_InSpillLEDCell.push_back(p);
    }  
    if( h.p1_InSpillLEDCell.size() != NCells() ) assert(false);
    dir_save->cd();
  }
  HistoInSpillLED &h = *histo_in_spill_led_;
  h.p1_InSpillLEDCell[icell]->Fill(0.5+frspill,framp);
}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void StoreLED::FitHistoInSpillLED( void )
{
  if(histo_in_spill_led_ == NULL ) return;
  HistoInSpillLED &h = *histo_in_spill_led_;
  if(h.p1_InSpillLEDCell.size() != NCells() ) return;
  for( size_t icell=0; icell< NCells(); icell++ )
  {
    if( h.p1_InSpillLEDCell[icell] == NULL ) assert(false);        
    double m0 = ledinspill0_[icell].GetMean();
    double m1 = ledinspill1_[icell].GetMean();
    double m2 = ledinspill2_[icell].GetMean();
    if( ledinspill0_[icell].GetEntries() > 0. )
    {
      if( store_leds_all_->led_thr[icell].GetMean() > 30. )
        h.h2fitpar_good4monitoring->Fill(m0, m1);
      else  
        h.h2fitpar_bad4monitoring->Fill(m0, m1);

      TF1 *fl = new TF1("fl","[0]+(x-0.5)*[1]",-0.1,1.1);
      fl->SetParameter(0,m0 );
      fl->SetParameter(1,m1 );
// In the run time it turned out that Add function Not implemented for TProfile  ...    h.p1_InSpillLEDCell[icell]->Add(fl);
    }  
// Apply polinomilal fit on the profile
//    h.p1_InSpillLEDCell[icell]
  }
}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void StoreLED::Process ( void )
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::ProcessStoredLED " << c_->GetName() << " debug  " <<  endl;
  store_leds_all_->Process();
  store_leds_in_spill_all_->Process();
     
  ledinspill0_.clear();
  ledinspill1_.clear();
  ledinspill2_.clear();
  StatInfo zero;

  std::vector<double>  stat(NCells(),0.);  
  std::vector<double>  ai(NCells(),0.);  
  std::vector<double> aiti(NCells(),0.);  
  std::vector<double> aiti2(NCells(),0.);  
  std::vector<double> ti(NCells(),0.);  
  std::vector<double> ti2(NCells(),0.);  
  std::vector<double> ti3(NCells(),0.);  
  std::vector<double> ti4(NCells(),0.);  
     
  if(  store_leds_all_->Processed() &&  store_leds_in_spill_all_->Processed() )
  {
    StoreCaloEvents &os = *store_leds_all_;
    StoreCaloEvents &is = *store_leds_in_spill_all_;
    size_t ncells = os.NCells();
     for( unsigned icell=0; icell < NCells(); icell++ )
     {
       ledinspill0_.push_back(zero);
       ledinspill1_.push_back(zero);
       ledinspill2_.push_back(zero);
     }
    size_t run_min = os.run_min_;
    size_t run_max = os.run_max_;
    for( unsigned irun=run_min; irun <= run_max; irun++ )
    {
      const std::pair<  int, int> & smima = os.GetSpillsInfo(irun);
      if( smima.first < 0 || smima.second < 0 ) continue;
      size_t spill_min = (size_t) smima.first;
      size_t spill_max = (size_t) smima.second;
      if( debug ) cout <<" run " << irun << endl;
      for( size_t ispill=spill_min; ispill <= spill_max; ispill++ )
      {
        std::pair< size_t, size_t > indx( irun, ispill);
       std::vector< CaloEvent  * > &evinspill = is.events_[indx];
        if( debug ) cout <<" spill " << ispill <<" evmin out spill " << os.statinspill_[indx].evmin <<" Events in spill " <<  evinspill.size() << endl;
        size_t evmaxinspill = os.statinspill_[indx].evmin;
        
        if( evmaxinspill  <= 0 ) continue;
        if( evmaxinspill  <= 100000. ) continue;   // Consider only high intencity spills. This is temporary check. Must be configured properly and improved.

        if( evinspill.size() == 0 ) continue;
        for( unsigned iev=0; iev < evinspill.size(); iev++ ) 
        {
          double frspill = double( evinspill[iev]->event_id_.GetEventNumberInBurst() )/double(evmaxinspill)-0.5;
          if( debug ) cout << iev << " ev # burst " << evinspill[iev]->event_id_.GetEventNumberInBurst() << " SpillFr " << frspill << endl;
          for( unsigned icell=0; icell < ncells; icell++ )
          {
            double ampos = os.led_thr[icell].GetMean();
            double sampos = os.led_thr[icell].GetSigma();
            if(  ampos < 20. ) continue;
            double amp =  evinspill[iev]->data_[icell];
            double framp = amp/ampos-1.;
            if( debug ) cout <<icell <<" frspill " << frspill <<" ampos " << ampos <<" sig proc.  " << 100.*sampos/ampos <<" framp " << framp << endl;
            bool goodled = false;
            if( ledtis0_.size() == NCells() )
            {
              if( ampos > 30. )
              {
                double m0 = ledtis0_[icell];
                double m1 = ledtis1_[icell];
                double m2 = ledtis2_[icell];
                double cledtis = 1./(1.+m0 + m1*frspill + m2*frspill*frspill);
                double ampc = amp*cledtis;
                double frampc = ampc/ampos-1.;
// test corrections                framp = frampc;
                goodled = true;  
                FillHistoInSpillLED(icell,framp,frspill,goodled);
              }
            }
            
            if( ! goodled ) FillHistoInSpillLED(icell,framp,frspill,goodled);

//             ledinspill0_[icell].Add(framp);
//             ledinspill1_[icell].Add(framp*frspill);
//             ledinspill2_[icell].Add(framp*frspill*frspill);

            stat[icell] += 1.;  
            ai[icell] += framp;  
            aiti[icell] += framp*frspill;  
            aiti2[icell] += framp*frspill*frspill;  
            ti[icell] += frspill;  
            ti2[icell] += frspill*frspill;  
            ti3[icell] += frspill*frspill*frspill;  
            ti4[icell] += frspill*frspill*frspill*frspill;  
          }  // end of cells cycle
        }  // end of events in spill cycle 
      } // end of spills cycle 
    }// end of runs cycle  

    int fit_corrections_power = 1;
    if(  fit_corrections_power ==  2 )
    {
      for( unsigned icell=0; icell < ncells; icell++ )
      {
        if( stat[icell] <= 0 ) continue;
        ai[icell] /= stat[icell];  
        aiti[icell] /= stat[icell]; 
        aiti2[icell] /= stat[icell]; 
        ti[icell] /= stat[icell];  
        ti2[icell] /= stat[icell];  
        ti3[icell] /= stat[icell];  
        ti4[icell] /= stat[icell];  
      
        double at = aiti[icell]-ai[icell]*ti[icell];
        double at2 = aiti2[icell]-ai[icell]*ti2[icell];
        double s2 =  ti2[icell] - ti[icell]*ti[icell];
        double s3 =  ti3[icell] - ti[icell]*ti2[icell];
        double s4 =  ti4[icell] - ti2[icell]*ti2[icell];
        double det = s4*s2 - s3*s3;
        if( det < 0.0001 )
        {
          cout <<" cell " << icell <<" " << c_->GetCellName(icell) << " det " << det << endl; 
          continue;
        }  
        double m2 = (at2*s2-at*s3)/det;
        double m1 = (at*s4-at2*s3)/det;
        double m0 = ai[icell]-m1*ti[icell]-m2*ti2[icell];
        ledinspill0_[icell].Add(m0);
        ledinspill1_[icell].Add(m1);
        ledinspill2_[icell].Add(m2);
      }
    }
    else if(  fit_corrections_power ==  0 )
    {
      for( unsigned icell=0; icell < ncells; icell++ )
      {
        if( stat[icell] <= 0 ) continue;
        ai[icell] /= stat[icell];  
        aiti[icell] /= stat[icell]; 
        aiti2[icell] /= stat[icell]; 
        ti[icell] /= stat[icell];  
        ti2[icell] /= stat[icell];  
        ti3[icell] /= stat[icell];  
        ti4[icell] /= stat[icell];  

        double m2 = 0.;
        double m1 = 0.;
        double m0 = ai[icell];
        ledinspill0_[icell].Add(m0);
        ledinspill1_[icell].Add(m1);
        ledinspill2_[icell].Add(m2);
      
      }
    }
    else if(  fit_corrections_power ==  1 )
    {
      for( unsigned icell=0; icell < ncells; icell++ )
      {
        if( stat[icell] <= 0 ) continue;
        ai[icell] /= stat[icell];  
        aiti[icell] /= stat[icell]; 
        aiti2[icell] /= stat[icell]; 
        ti[icell] /= stat[icell];  
        ti2[icell] /= stat[icell];  
        ti3[icell] /= stat[icell];  
        ti4[icell] /= stat[icell];  
        double at = aiti[icell]-ai[icell]*ti[icell];
        double at2 = aiti2[icell]-ai[icell]*ti2[icell];
        double s2 =  ti2[icell] - ti[icell]*ti[icell];
        double s3 =  ti3[icell] - ti[icell]*ti2[icell];
        double s4 =  ti4[icell] - ti2[icell]*ti2[icell];

        double m2 = 0.;
        double m1 = at/s2;
        double m0 = ai[icell]-m1*ti[icell];
        ledinspill0_[icell].Add(m0);
        ledinspill1_[icell].Add(m1);
        ledinspill2_[icell].Add(m2);
      }
    } 
    
    FitHistoInSpillLED();

  } 
}

////////////////////////////////////////////////////////////////////////////////
// Real corrections reading will be implemented later
//
int StoreLED::InputTimeInSpillLED(const string &s)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) 
  {
     cout <<" StoreLED::InputTimeInSpillLED " << c_->GetName() << " debug " << endl;
  }

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  getline(is,str);

//  Read calibration
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell,stat;
    float m0,m1,m2;

    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    int ret = sscanf(str.c_str(), " %d %g %g %g %d %s", &icell, &m0, &m1, &m2, &stat, cellname);
    assert( ret == 6 );
    int jcell = c_->FindCellByName( cellname );

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << c_->GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }
      
//     if(icell >=0 && icell < (int)NCells() )
//       c_->SetCorrTimeInSpillLED(icell,m0/10000.,m1/10000.,m2/10000., double(stat) );
  }
  if( debug ) 
  {
     cout <<" StoreLED::InputTimeInSpillLED " << c_->GetName() << " debug OK " << endl;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int StoreLED::OutputTimeInSpillLED( string &s ) const
{
 bool debug = false;
//   bool debug = true;
//  if( c_->GetName() == "EC02P1__") debug = true;
  if( debug )
  {
    cout << " StoreLED::OutputTimeInSpillLED " << c_->GetName() << " debug " << endl;
  }
// check that processed store_leds_all_
  if( !( store_leds_all_->Processed() &&  store_leds_in_spill_all_->Processed())  )
  {
    cout << " ERROR in StoreLED::OutputTimeInSpillLED " << c_->GetName() << " StoreCaloEvents not Processed " << endl;
    return -1;
  }  
  if( ledinspill0_.size() != NCells() || ledinspill1_.size() != NCells()  || ledinspill2_.size() != NCells()  )
  {
    cerr  << " ERROR in StoreLED::OutputTimeInSpillLED " << c_->GetName() << "  ledinspill size not consistent " << endl;
    cerr << " NCells " << NCells() << " but " << ledinspill0_.size() << " " << ledinspill1_.size() << " " << ledinspill2_.size() << " " << endl;
    exit( -1);
  }  
  std::ostringstream o;
  o << "### " << std::endl;
  o << "LEDs Calorimeter: " << c_->GetName() << " Ncells=" << NCells() << " " << std::endl;
  o << "  Cell    InSpill:   Mean0   Mean1  Mean2   Statistic   CellName " << std::endl;

  if( debug )
  {
     cout <<" LEDs Calorimeter: " << c_->GetName() <<" Ncells= " << NCells() << endl;
     cout <<"    Cell    InSpill:   Mean0   Mean1  Mean2   Statistic   CellName  " << endl;
  }

  for( size_t icell=0; icell < NCells(); icell++ )
  {
     double entries = ledinspill0_[icell].GetEntries();
     int stat = (int)entries;
     double m0 = (10000.*ledinspill0_[icell].GetMean());
     double m1 = (10000.*ledinspill1_[icell].GetMean());
     double m2 = (10000.*ledinspill2_[icell].GetMean());
       if( debug ) cout <<"Cell  " <<icell << " m0 " << m0 << " m1 " << m1 <<
                                                                                " m2 " << m2 <<" Stat " << entries <<" " << c_->GetCellName(icell) <<endl;

       o << " " << icell << " "
         << std::fixed << std::setw(7) << std::setprecision(2) << m0 << "  "
         << std::fixed << std::setw(7) << std::setprecision(2) << m1 << "  "
         << std::fixed << std::setw(7) << std::setprecision(2) << m2 << "  "
         << std::setw(10) << stat << "   "
         << c_->GetCellName(icell) << " " << std::endl;

  }
  s = o.str();
  if( debug )
  {
    cout << " StoreLED::OutputTimeInSpillLED " << c_->GetName() << "  debug finished OK " << endl;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int StoreLED::OutputLEDInfoInSpills( string &s, const string &comment ) const
{
 bool debug = false;
//   bool debug = true;
  
  if( debug )
  {
    cout << " StoreLED::OutputLEDInfoInSpills " << c_->GetName() << " debug " << endl;
  }
//  if( !options.store_leds_all ) return -1;
  if( store_leds_all_ == NULL )  return -1;
// check that processed store_leds_all_
  if( !store_leds_all_->Processed() ) 
  {
    cout << " ERROR in StoreLED::OutputLEDInfoInSpills " << c_->GetName() << " StoreCaloEvents store_leds_all_ was not Processed " << endl;
    return -1;
  }
  if( debug )
  {
    cout << " StoreLED::OutputLEDInfoInSpills " << c_->GetName() << " Start the output  " << endl;
  }
  std::ostringstream o;
  if( comment == "") { 
    o << "### " << std::endl;
  } else {
    o << comment.c_str() << std::endl;
  }
  o << "LEDs Calorimeter: " << c_->GetName() << " Ncells=" << c_->NCells() << " " << std::endl;
  o << "  Cell      MeanPeak   SigmaPeak  Statistic  NormFlag  CellName " << std::endl;

  if( debug )
  {
     cout <<" LEDs Calorimeter: " << c_->GetName() <<" Ncells= " << c_->NCells() << endl;
     cout <<"  Cell      MeanPeak   SigmaPeak  Statistic  NormFlag  CellName   " << endl;
  }

  for( size_t icell=0; icell < c_->NCells(); icell++ )
  {
     int normflag = 1;
     const StatInfo &led_strict = store_leds_all_->GetStrictAverageValue( icell );    
//     double entries = led_strict.GetEntries();
       if( debug ) cout <<"Cell  " <<icell << " Mean " << led_strict.GetMean()<< " Sig " << led_strict.GetSigma() <<
                                                                                " Stat " << led_strict.GetEntries() <<" FL " << normflag <<" " << c_->GetCellName(icell) <<endl;

       o << " " << icell << " "
         << std::fixed << std::setw(7) << std::setprecision(2) << led_strict.GetMean() << " "
         << std::fixed << std::setw(7) << std::setprecision(2) << led_strict.GetSigma() << " "
         << std::fixed << std::setw(7) << std::setprecision(2) << led_strict.GetEntries() << " "
         << std::setw(5) << normflag << "   "
         << c_->GetCellName(icell) << " " << std::endl;

     for( unsigned irun=store_leds_all_->run_min_; irun <= store_leds_all_->run_max_; irun++ )
     {
       if( debug ) cout <<" run " << irun <<" spills " << store_leds_all_->spill_max_-1 << endl;
       o << " run " << irun << " spills " << store_leds_all_->spill_max_ << " " << std::endl;
       int icnt20 = 0;
       for( unsigned ispill= 1; ispill <= store_leds_all_->spill_max_; ispill++ )
       {
         const StatInfo &leds_in_spill =  store_leds_all_->GetValue( irun, ispill, icell );

          if( leds_in_spill.GetEntries() > 0 && led_strict.GetMean() > 0. )
          {
        
            double ddd = 1000.*(leds_in_spill.GetMean()/led_strict.GetMean()-1.);
//            cout <<" " << int(ddd);
            o << " " << int(ddd);
          }
          else
          {
//            cout <<" nan " << endl;
            o << " nan ";
          }
          icnt20++;
          if( icnt20 == 20 )
          {
            icnt20 = 0;
            o << std::endl;
//            cout << endl;   
          }
        }
        if( icnt20 != 0 ) o << std::endl;
//        cout << endl;
      }
  }
  s = o.str();
  if( debug )
  {
    cout << " StoreLED::OutputLEDInfoInSpills " << c_->GetName() << "  debug finished OK " << endl;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int StoreLED::OutputLEDInfoInSpillsXY( string &s, const string &comment ) const
{
  bool debug = false;
//   bool debug = true;
  
  if( debug )
  {
    cout << " StoreLED::OutputLEDInfoInSpills " << c_->GetName() << " debug " << endl;
  }
//   if( !options.store_leds_all ) return -1;
// check that processed store_leds_all_
  if( !store_leds_all_->Processed() ) 
  {
    cout << " ERROR in StoreLED::OutputLEDInfoInSpills " << c_->GetName() << " StoreCaloEvents store_leds_all_ was not Processed " << endl;
    return -1;
  }
  std::ostringstream o;
  o << "### " << std::endl;
  o << "LEDs Calorimeter: " << c_->GetName() << " Ncells=" << c_->NCells() << " " << std::endl;
  o << "  X  Y  Cell     MeanPeak   SigmaPeak  Statistic  NormFlag  CellName " << std::endl;

  assert( c_->XYRegularGrid() );

  for( size_t icell=0; icell!=c_->NCells(); icell++ )
  {
     int normflag = 1;
     int x = -1;     
     int y = -1;
     if( c_->XYRegularGrid() )
     {     
       x = c_->GetColumnOfCell(icell);     
       y = c_->GetRowOfCell(icell);
     }     
     const StatInfo &led_strict = store_leds_all_->GetStrictAverageValue( icell );    
//     double entries = led_strict.GetEntries();
       o << " " << x << " "
         << y << " "
         << icell << " "
         << std::fixed << std::setw(7) << std::setprecision(2) << led_strict.GetMean() << " "
         << std::fixed << std::setw(7) << std::setprecision(2) << led_strict.GetSigma() << " "
         << std::fixed << std::setw(7) << std::setprecision(2) << led_strict.GetEntries() << " "
         << std::setw(5) << normflag << "   "
         << c_->GetCellName(icell) << " " << std::endl;
     for( unsigned irun=store_leds_all_->run_min_; irun <= store_leds_all_->run_max_; irun++ )
     {
       if( debug ) cout <<" run " << irun <<" spills " << store_leds_all_->spill_max_-1 << endl;
       o << " run " << irun << " spills " << store_leds_all_->spill_max_ << " " << std::endl;
       int icnt20 = 0;
       for( unsigned ispill= 1; ispill <= store_leds_all_->spill_max_; ispill++ )
       {
         const StatInfo &leds_in_spill =  store_leds_all_->GetValue( irun, ispill, icell );
          if( leds_in_spill.GetEntries() > 0 && led_strict.GetMean() > 0. )
          {
        
            double ddd = 1000.*(leds_in_spill.GetMean()/led_strict.GetMean()-1.);
//            cout <<" " << int(ddd);
            o << " " << int(ddd);
          }
          else
          {
//            cout <<" nan " << endl;
            o << " nan ";
          }
          icnt20++;
          if( icnt20 == 20 )
          {
            icnt20 = 0;
//            cout << endl;   
            o << std::endl;
          }
        }
        if( icnt20 != 0 ) o << std::endl;
//        cout << endl;
      }
  }
  s = o.str();
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

