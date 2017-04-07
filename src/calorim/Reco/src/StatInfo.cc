/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/StatInfo.cc,v $
   $Date: 2010/10/05 18:10:14 $
   $Revision: 1.13 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

   Copyright(C): 2002  V.Kolosov,A.Zvyagin

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
#include <iostream>
#include <cstdio>
#include <cmath>
#include <sstream>

// --- Internal files ----
#include "StatInfo.h"
#include "Exception.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

StatInfo::StatInfo(void)
{
  Clear();
}

////////////////////////////////////////////////////////////////////////////////

StatInfo::StatInfo(double weights_sum, double mean, double sigma)
{
  momenta[0] = weights_sum;
  momenta[1] = weights_sum * mean;
  momenta[2] = weights_sum * (sigma*sigma+mean*mean);
}

////////////////////////////////////////////////////////////////////////////////

void StatInfo::Clear(void)
{
  for(size_t i=0; i<momenta_size; i++ )
    momenta[i] = 0;
}

////////////////////////////////////////////////////////////////////////////////

void StatInfo::Set(double weights_sum, double mean, double sigma)
{
  momenta[0] = weights_sum;
  momenta[1] = weights_sum * mean;
  momenta[2] = weights_sum * (sigma*sigma+mean*mean);
}

////////////////////////////////////////////////////////////////////////////////

double StatInfo::GetMean(void) const
{
  if(momenta[0]!=0)
  {
    return momenta[1]/momenta[0];
  }
  else
  {
    return momenta[1];
  }
}

////////////////////////////////////////////////////////////////////////////////

void StatInfo::Add(double a, double weight)
{
  double mom=weight;
  for(size_t i=0; i<momenta_size; i++ )
  {
    momenta[i] += mom;
    mom *= a;
  }
}

////////////////////////////////////////////////////////////////////////////////

void StatInfo::Translate(double c)
{
  momenta[2] += 2.*momenta[1]*c+momenta[0]*c*c;
  momenta[1] += momenta[0]*c;
}

////////////////////////////////////////////////////////////////////////////////

void StatInfo::Scale(double c)
{
  momenta[1] *= c;
  momenta[2] *= c*c;
}

////////////////////////////////////////////////////////////////////////////////

double StatInfo::GetSigma(void) const
{
  if(momenta[0]==0) return momenta[2];
  const double mean = GetMean();
  double f = 1.;
  if( GetEntries() > 1. ) f= GetEntries()/(GetEntries()-1.);
  return sqrt( f*fabs(momenta[2]/GetEntries() - mean*mean) );  // 'abs' used to overcome precison problems
}

////////////////////////////////////////////////////////////////////////////////

void StatInfo::Result(double &n,double &mean,double &sigma) const
{
  n = StatInfo::momenta[0];
  if( n==0 )
  {
     mean  =0;
     sigma =0;
     return;
  }
  mean  = GetMean();
  sigma = GetSigma();
}

////////////////////////////////////////////////////////////////////////////////

ostream &operator << (ostream &o, const StatInfo &c)
{
  char s[200];
  sprintf(s,"mean=%9.3e   sigma=%9.3e   N=mom0=%9.3e",
          c.GetMean(),c.GetSigma(),c.GetEntries());
  o << s;
  return o;
}

////////////////////////////////////////////////////////////////////////////////

istream & operator >> (istream &in, StatInfo &c)
{
  c.Clear();
  string s;
  char dummy[111];

  getline(in,s);

  size_t len;
  sscanf(s.c_str(),"%s %s %c %zu",dummy,dummy,dummy,&len);

  return in;
}

////////////////////////////////////////////////////////////////////////////////

StatInfo &StatInfo::operator = (const StatInfo &c)
{
  if( &c!=this )
  {
    momenta[0]  = c.momenta[0];
    momenta[1]  = c.momenta[1];
    momenta[2]  = c.momenta[2];
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

StatInfo &StatInfo::operator += (const StatInfo &c)
{
  momenta[0]  += c.momenta[0];
  momenta[1]  += c.momenta[1];
  momenta[2]  += c.momenta[2];
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

StatInfoStore::StatInfoStore(void) :
  run_start_(0),
  run_finish_(-1),
  monitor_hist_(0),
  min_range_(0),
  min_stat_(1),
  min_disp_add_(0),
  min_disp_frac_(0),
  data_updated_(false)
{
}

////////////////////////////////////////////////////////////////////////////////

StatInfoStore::StatInfoStore(int run_start, int run_finish) :
  run_start_(run_start),
  run_finish_(run_finish),
  monitor_hist_(0),
  min_range_(0),
  min_stat_(1),
  min_disp_add_(0),
  min_disp_frac_(0),
  data_updated_(false)
{
  Init();
}

////////////////////////////////////////////////////////////////////////////////

StatInfoStore::StatInfoStore(int run_start, vector<StatInfo> &stat_info) :
  run_start_(run_start),
  stat_info_(stat_info),
  monitor_hist_(0),
  min_range_(0),
  min_stat_(1),
  min_disp_add_(0),
  min_disp_frac_(0),
  data_updated_(false)
{
  run_finish_ = run_start + stat_info.size() - 1;
}

////////////////////////////////////////////////////////////////////////////////

StatInfoStore &StatInfoStore::operator = (const StatInfoStore &c)
{
  if( &c!=this )
  {
    run_start_    = c.run_start_;
    run_finish_   = c.run_finish_;
    stat_info_    = c.stat_info_;
    points_       = c.points_;
    for( int i=0; i < 100; i++)
      parameters_[i] = c.parameters_[i];
    monitor_hist_ = 0;  // to avoid memory leaks
    min_range_    = c.min_range_;
    min_stat_     = c.min_stat_;
    min_disp_add_  = c.min_disp_add_;
    min_disp_frac_ = c.min_disp_frac_;
    data_updated_ = true;
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

StatInfoStore &StatInfoStore::operator += (const StatInfoStore &c)
{
  data_updated_ = true;
  if( run_start_ == c.run_start_ && run_finish_ == c.run_finish_ )
  {
    for( int i = 0; i < Size(); i++)
      stat_info_[i] += c.stat_info_[i];
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::Init( int run_start, int run_finish )
{
  stat_info_.clear();
  points_.clear();
  run_start_ = run_start;
  run_finish_ = run_finish;
  if( run_finish_ >= run_start_ )
  {
    points_.insert( points_.end(), 0 );
    points_.insert( points_.end(), Size() );
    points_.sort();
  }
  Init();
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::Init( void )
{
  for( int i=0; i < Size(); i++ )
    stat_info_.push_back(StatInfo());
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::Clear(void)
{
  data_updated_ = true;
  for( int i=0; i < Size(); i++)
    stat_info_[i].Clear();
  points_.clear();
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::SetIntervals(list<int> &points)
{
  points_.clear();
  if( Size() <= 0 ) return;
  for( list<int>::const_iterator it=points.begin(); it != points.end(); it++)
  {
     if( (*it) >= 0 && (*it) <= Size() )  points_.insert( points_.end(), *it );
  }
  points_.sort();
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::SetInfo(int run_start, vector<StatInfo> &stat_info)
{
  data_updated_ = true;
  if( monitor_hist_ != 0 )
  {
    if( run_start_ != run_start || (int)stat_info.size() != Size() )
    {
      cerr << " ERROR in StatInfoStore::SetInfo  new DATA but old Histos !!!! " << endl;
      assert(false);
    }
  }
  run_start_= run_start;
  run_finish_ = run_start + stat_info.size() - 1;
  stat_info_ = stat_info;
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::SetInfo(int run, StatInfo &stat_info)
{
  data_updated_ = true;
  if( run < run_start_ || run > run_finish_ ) return;
  int bin = run - run_start_;
  if( (unsigned)bin >= stat_info_.size() ) assert(false);
  stat_info_[bin] = stat_info;
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::Add(int run, StatInfo &stat_info)
{
  data_updated_ = true;
  if( run < run_start_ || run > run_finish_ ) return;
  int bin = run - run_start_;
  if( (unsigned)bin >= stat_info_.size() ) assert(false);
  stat_info_[bin] += stat_info;
}

/////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::Add(int run, double data)
{
  data_updated_ = true;
  if( run < run_start_ || run > run_finish_ ) return;
  int bin = run - run_start_;
  if( (unsigned)bin >= stat_info_.size() ) assert(false);
  stat_info_[bin].Add(data);
}

///////////////////////////////////////////////////////////////////////////////

StatInfo StatInfoStore::CalculateAverage(void)
{
  StatInfo average;
  for( size_t i=0; i < stat_info_.size(); i++ )
  {
    if(stat_info_[i].GetEntries() >= min_stat_ )
      average += stat_info_[i];
  }
  return average;
}

////////////////////////////////////////////////////////////////////////////////

StatInfo StatInfoStore::CalculateAverageInInterval( int run )
{
  if(points_.empty())
  {
    cerr << " ERROR:: StatInfoStore seems not initialsed points_.empty() !! " << endl;
    if( RunInRange( run ) )
    {
      return CalculateAverage();
    }
    else
    {
      StatInfo empty;
      return empty;
    }
  }

  list<int>::const_iterator it;
  list<int>::const_iterator it_previous;
  int bin = run - run_start_;
  int range_min = 0;
  int range_max = -1;
  for( it=points_.begin(); it != points_.end(); it++)
  {
     if( (*it) > bin )
     {
       range_max = (*it);
       break;
     }
     it_previous = it;
  }
  if( it != points_.begin() )
  {
    range_min = (*it_previous);
    return CalculateAverage( range_min + run_start_, range_max + run_start_ - 1 );
  }
  else
  {
    StatInfo empty;
    return empty;
  }

}

////////////////////////////////////////////////////////////////////////////////

StatInfo StatInfoStore::CalculateAverage( int run_min, int run_max )
{
  StatInfo average;
  int bin_min = run_min - run_start_;
  int bin_max = run_max - run_start_;

  if( bin_max < bin_min ) return average;
  for( int i = bin_min; i <= bin_max; i++ )
  {
    if( i >= 0 && i < Size() )
    {
      if( stat_info_[i].GetEntries() >= min_stat_ )
      {
        average += stat_info_[i];
      }
    }
  }
  return average;
}

////////////////////////////////////////////////////////////////////////////////

double StatInfoStore::EstimateSigma(int from_bin, int to_bin )
{
  StatInfo delta;
  int min = from_bin;
  int max = to_bin;
  if(min >= max) return 0;
  if( min < 0 ) min=0;
  if( max > (int)stat_info_.size()-1 ) max=stat_info_.size()-1;

  for( int i=min; i < max; i++ )
  {
    if(stat_info_[i].GetEntries() >= min_stat_  && stat_info_[i+1].GetEntries() >= min_stat_)
    {
      double d = stat_info_[i].GetMean() - stat_info_[i+1].GetMean();
      delta.Add( d * d );
    }
  }
  return sqrt( delta.GetMean()/2. );
}

////////////////////////////////////////////////////////////////////////////////

StatInfo StatInfoStore::MakeRelativeNormalizationToAverage(void)
{
  StatInfo average = CalculateAverage();
  if( average.GetMean() == 0 ) assert(false);
  for( size_t i=0; i!=stat_info_.size(); i++ )
  {
    stat_info_[i].Translate(-average.GetMean());
    stat_info_[i].Scale(1./average.GetMean());
  }
  return average;
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::FindIntervals(int max_intervals,double hi_break,double hi_penalty )
{
  list<int> lpoints;
  lpoints.insert(lpoints.end(),0);
  lpoints.insert(lpoints.end(),stat_info_.size());
  double params[200];
  double hi0 = FitInIntervals( lpoints, params );
//   cout << " hi0 " << hi0;
  list<int> lpoints1 = lpoints;
  list<int> lpoints2 = lpoints;
  double hi,hi1;
  hi = hi0;
  for(int i=0; i< max_intervals; i++)
  {
    if(lpoints1.size() > (unsigned)maximum_intervals) break;
    hi1 = AddNewPoint( lpoints1, params, hi_penalty );
    if( hi - hi1 < hi_break ) break;
    hi = hi1;
//     cout << " hi" << i+1 << " " << hi1;
//    hi2 = AddNewInterval( lpoints, params, hi_penalty );
  }
//   cout << endl;
  points_ = lpoints1;
  for(int i=0; i<200; i++)
    parameters_[i] = params[i];
}

////////////////////////////////////////////////////////////////////////////////

double StatInfoStore::FitInInterval(int from_bin,int to_bin,double *parameters)
{
  static const double eps = 1.e-10;
//  double sigma_estimated = EstimateSigma( from_bin, to_bin);
  parameters[0]=0.;
  parameters[1]=0.;
  int min = from_bin;
  int max = to_bin;
  if(min >= max) return 0;
  if( min < 0 ) min=0;
  if( max > (int)stat_info_.size()-1 ) max=stat_info_.size()-1;

  double sum_weight=0.;
  double sum_x=0.;
  double sum_y=0.;
  int nbin = max-min+1;
  if(nbin < 1) return 0;
  if(nbin == 1)
  {
    parameters[0]=stat_info_[min].GetMean();
    return 0;
  }
  vector < double > ww(nbin);
  vector < double > x (nbin);
  vector < double > y (nbin);

  for( int it=0; it < nbin; it++)
  {
    double w = 0;
    if(stat_info_[min+it].GetEntries() < min_stat_ )
    {
      y[it] = 0;
      w=0;  // Set weight to zero; Will be ignored in fit
    }
    else
    {
      y[it] = stat_info_[min+it].GetMean();
      double add_disp = min_disp_add_ + min_disp_frac_ * y[it] * y[it]; // add some constants to dispersion to be more robust
//      add_disp = 0.;
//      w = stat_info_[min+it].GetSigma()/sqrt(stat_info_[min+it].GetEntries())+0.2*sigma_estimated;
      w = stat_info_[min+it].GetSigma()/sqrt(stat_info_[min+it].GetEntries());
      w = 1./(w*w+add_disp);
    }

    x[it] = it;
    ww[it] = w;
    sum_weight += w;
    sum_x += x[it]*w;
    sum_y += y[it]*w;
  }

  if(sum_weight < eps) return 0;

  double average_x = sum_x/sum_weight;
  double average_y = sum_y/sum_weight;

  double sum_xx=0.;
  double sum_xy=0.;
  double sum_yy=0.;

  for( int it=0; it < nbin; it++)
  {
    double xx = x[it]-average_x;
    double yy = y[it]-average_y;
    sum_xx += xx*xx*ww[it];
    sum_xy += xx*yy*ww[it];
    sum_yy += yy*yy*ww[it];
  }

  double a = sum_xy/sum_xx;
  double b = a*( (double)(nbin-1)/2. - average_x ) + average_y;
  parameters[0]=b;
  parameters[1]=a;
  return a*a*sum_xx - 2.*a*sum_xy + sum_yy;
}

////////////////////////////////////////////////////////////////////////////////

double StatInfoStore::FitInIntervals( list<int> &points, double *parameters )
{
  if(points.size() < 2) return 0;
  double hi = 0;
  list<int>::const_iterator it = points.begin();
  for( int i = 0; i < (int)points.size()-1; i++)
  {
    int from_point = *it++;
    int to_point = *it;
    hi += FitInInterval( from_point, to_point - 1, &parameters[2*i]);
  }
  return hi;
}

////////////////////////////////////////////////////////////////////////////////

double StatInfoStore::AddNewPoint( list<int> &points, double *parameters, double penalty )
{
  if( points.size() < 2 ) assert(false);
  int min = *points.begin();
  int max = *( --points.end() );
  double hi_best = FitInIntervals( points, parameters );
  int ibest = -1;
  for( int i = min; i < max; i++)
  {
    list<int>::const_iterator it;
    double hi_penalty = 0;
    for( it = points.begin(); it != points.end(); it++)
    {
      if( fabs( double(*it-i) ) < min_range_ ) hi_penalty += penalty*(min_range_ - fabs( double(*it-i) ));
    }
    if( it != points.end()) continue;

    list<int> points_plus = points;
    points_plus.insert( points_plus.end(), i );
    points_plus.sort();
    double hi = FitInIntervals( points_plus, parameters ) + hi_penalty;

    if(hi < hi_best)
    {
      hi_best = hi;
      ibest = i;
    }
  }

  if( ibest == -1)
  {
//     cerr << " StatInfoStore::AddNewPoint just no need to add new point to "<< points.size() <<
//                                                                " hi_best " << hi_best  << endl;
    return hi_best;  // possibly statistic is empty
  }

  points.insert( points.end(), ibest );
  points.sort();
  return FitInIntervals( points, parameters );
}

////////////////////////////////////////////////////////////////////////////////

double StatInfoStore::AddNewInterval( list<int> &points, double *parameters, double penalty )
{
  if( points.size() < 2 )  assert(false);
  int min = *points.begin();
  int max = *(--points.end());
  double hi_best = FitInIntervals( points, parameters );
  int ibest = -1;
  int jbest = -1;
  for( int i = min; i < max-1; i++)
  {
    list<int>::const_iterator it;
    for( it = points.begin(); it != points.end(); it++)
    {
      if( fabs( double(*it-i) ) < min_range_ ) break;
    }
    if( it != points.end()) continue;
    for( int j = i+1; j < max; j++)
    {
      double hi_penalty = 0;
      for( it = points.begin(); it != points.end(); it++)
      {
        if( fabs( double(*it-i) ) < min_range_ ) hi_penalty += penalty*(min_range_ - fabs( double(*it-i) ));
        if( fabs( double(*it-j) ) < min_range_ ) hi_penalty += penalty*(min_range_ - fabs( double(*it-j) ));
      }
      list<int> points_plus = points;
      points_plus.insert(points_plus.end(),i);
      points_plus.insert(points_plus.end(),j);
      points_plus.sort();
      double hi = FitInIntervals( points_plus, parameters ) + hi_penalty;
      if(hi < hi_best)
      {
        hi_best = hi;
        ibest = i;
        jbest = j;
      }
    }
  }
  if( ibest == -1 || jbest == -1 )
  {
//     cerr << " StatInfoStore::AddNewInterval just no need to add new point to "<< points.size() <<
//                                                                " hi_best " << hi_best  << endl;
    return hi_best;  // possibly statistic is empty
  }
  points.insert(points.end(),ibest);
  points.insert(points.end(),jbest);
  points.sort();
  return FitInIntervals( points, parameters );
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::BookHisto(const string &name, const string &histo_path, int histo_level )
{
  if( Size() <= 0 ) return; // StatInfoStore has no data
  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
  if( monitor_hist_ == 0 )
  {
    monitor_hist_ = new Monitor_Histo;
    Monitor_Histo* h = monitor_hist_;
    if(!gDirectory->cd("/")) assert(false);
    if(histo_path.empty())
    {
      h->root_dir = gDirectory;
    }
    else
    {
      h->root_dir=myROOT_utils::TDirectory_checked(histo_path.c_str());
    }
    if(!h->root_dir->cd()) assert(false);

    char hist_name[132];
    char hname[132];
    if( histo_level >= 0 )
    {
      sprintf(hname,"%s",name.c_str());
      sprintf(hist_name," Monitor value %s ", name.c_str() );
      h->h1_monitor = myROOT_utils::TH1D_checked(hname,hist_name,run_finish_-run_start_+1,run_start_-0.5,run_finish_+0.5);
    }
    if( histo_level > 0 )
    {
      sprintf(hname,"Norm_%s",name.c_str());
      sprintf(hist_name," Monitor normalised %s ", name.c_str() );
      h->h1_monitor_norm = myROOT_utils::TH1D_checked(hname,hist_name,run_finish_-run_start_+1,run_start_-0.5,run_finish_+0.5);
    }
    if( histo_level > 0 )
    {
      sprintf(hname,"Entries_%s",name.c_str());
      sprintf(hist_name," Monitor entries %s ", name.c_str() );
      h->h1_monitor_entries = myROOT_utils::TH1D_checked(hname,hist_name,run_finish_-run_start_+1,run_start_-0.5,run_finish_+0.5);
    }
    if( histo_level > 0 )
    {
      sprintf(hname,"Sigma_%s",name.c_str());
      sprintf(hist_name," Monitor sigma %s ", name.c_str() );
      h->h1_monitor_sigma = myROOT_utils::TH1D_checked(hname,hist_name,run_finish_-run_start_+1,run_start_-0.5,run_finish_+0.5);
    }
  }
  if( !dir_save->cd() ) assert(false);
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::FillHisto(void)
{
//   bool debug = true;
  bool debug = false;
  if( monitor_hist_ == 0 ) return;
  data_updated_ = false;
  double sigma_estimated = EstimateSigma( 0, stat_info_.size()-1);
  Monitor_Histo* h = monitor_hist_;
  StatInfo average = CalculateAverage();
  if( debug ) cout << " StatInfoStore::FillHisto size =" << stat_info_.size() << endl;
  for( size_t i=0; i!=stat_info_.size(); i++ )
  {
    if( debug ) cout << " i =" << i << endl;
    if( stat_info_[i].GetEntries() > min_stat_ )
    {
      double mean = stat_info_[i].GetMean();
      if( debug ) cout << " mean =" << mean << " entries " << stat_info_[i].GetEntries() << " Sigma " << stat_info_[i].GetSigma() << endl;
      double sigma = stat_info_[i].GetSigma()/sqrt(stat_info_[i].GetEntries()) + 0.2*sigma_estimated;
      double add_disp = min_disp_add_+ min_disp_frac_*(mean*mean); // add some constants to dispersion to be more robust
      sigma = sqrt( sigma*sigma + add_disp );
      if(h->h1_monitor != 0)
        h->h1_monitor->SetBinContent( i, mean);
      if(h->h1_monitor != 0)
        h->h1_monitor->SetBinError( i, sigma);
      double mean_average = average.GetMean();
      if( debug ) cout << " mean_average =" << mean_average << " h1_monitor_norm = " << h->h1_monitor_norm << endl;
      if(h->h1_monitor_norm != 0)
        h->h1_monitor_norm->SetBinContent( i, mean/mean_average - 1.);

      if(h->h1_monitor_norm != 0)
        h->h1_monitor_norm->SetBinError( i, sigma/mean_average );

      if( debug ) cout << " entries =" << stat_info_[i].GetEntries() << endl;
      if(h->h1_monitor_entries != 0)
        h->h1_monitor_entries->SetBinContent( i, stat_info_[i].GetEntries());

      if( debug ) cout << " sigma =" << stat_info_[i].GetEntries() << endl;
      if(h->h1_monitor_sigma != 0)
        h->h1_monitor_sigma->SetBinContent( i, stat_info_[i].GetSigma() );
      if( debug ) cout << " sigma/sqrt(N) =" << endl;
      if(h->h1_monitor_sigma != 0)
        h->h1_monitor_sigma->SetBinError( i, stat_info_[i].GetSigma()/sqrt(stat_info_[i].GetEntries()));
    }
    if( debug ) cout << " i =" << i << " OK " <<  endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::DrawHisto(void)
{
  if( Size() <= 0 ) return;
  if( monitor_hist_ == 0 ) return;
  Monitor_Histo* h = monitor_hist_;
  if( data_updated_ ) FillHisto();
  if( h->h1_monitor != NULL ) h->h1_monitor->Draw();
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::AddFitToHisto(void)
{
  if( monitor_hist_ == 0 ) return;
  if( points_.size() < 2 ) return;
  Monitor_Histo* h = monitor_hist_;
  list<int>::const_iterator list_it = points_.begin();
  for( size_t i=0; i < points_.size()-1; i++)
  {
    double min =run_start_+(*list_it++)-1;
    double max =run_start_+(*list_it)-1;
//kill    strstream hname;
    char hname[132];
//kill    hname.form("f_%d",i);
    sprintf(hname,"f_%zu",i);
//kill    TF1 f = TF1(hname.str(),"[0]+[1]*x",min,max);
    TF1 f = TF1(hname,"[0]+[1]*x",min,max);
    double x_middle = (max+min)/2.;
    double b_inner = parameters_[i*2];
    double a = parameters_[i*2+1];
    double b = -a*x_middle + b_inner;
//    f.FixParameter(0,b);
//    f.FixParameter(1,a);
    f.SetParameters(b,a);
//kill    if( i == 0 ) h->h1_monitor->Fit(hname.str(),"QOR");
//kill    else h->h1_monitor->Fit(hname.str(),"QOR+");
    if( i == 0 ) h->h1_monitor->Fit(hname,"QOR");
    else h->h1_monitor->Fit(hname,"QOR+");
    const Int_t kNotDraw = 1<<9;
//kill    h->h1_monitor->GetFunction(hname.str())->ResetBit(kNotDraw);
    h->h1_monitor->GetFunction(hname)->ResetBit(kNotDraw);
  }
}

////////////////////////////////////////////////////////////////////////////////

double StatInfoStore::CalculateFitValue( int run )
{
  if(points_.empty()) return 0;
  list<int>::const_iterator it;
  list<int>::const_iterator it_previous;
  int bin = run - run_start_;
//  int range_min = 0;
  int range_max = -1;
  int interval = -1;
  for( it=points_.begin(); it != points_.end(); it++)
  {
     if( (*it) > bin )
     {
       range_max = (*it);
       break;
     }
     it_previous = it;
     interval++;
  }
  if( it != points_.begin() )
  {
    double min = (*it_previous);
    double max = (*it);
    double x_middle = (max+min)/2.;

    double b_inner = parameters_[interval*2];
    double a = parameters_[interval*2+1];
//    cout << " interval " << interval << " a=" << a << " b=" << b_inner;
    double value = a*( bin - x_middle ) + b_inner;
//    cout << " Fit value " << value << " Measured value " << stat_info_[bin].GetMean();
    return value;
  }
  else
  {
    return 0;
  }
}

////////////////////////////////////////////////////////////////////////////////

StatInfo StatInfoStore::GetBinValue( int run ) const
{
  int bin = run - run_start_;
  return stat_info_[bin];
}

////////////////////////////////////////////////////////////////////////////////

void StatInfoStore::Print( void )
{
  cout << " StatInfoStore::PrintOut for the Range: run_min = " << run_start_ << " -- run_max= " << run_finish_ << endl;
  cout << " Total of " << NIntervals() << " of intervals " << endl;
  list<int>::const_iterator it = points_.begin();
  for(int i=0; i < NIntervals(); i++ )
  {
    cout << " From " << (*it++);
    cout << " to " << (*it)-1 << endl;
  }

  for( int run = run_start_; run <= run_finish_; run++ )
  {
    int bin = run - run_start_;
    cout << " For run " << run << " Fit value=" << CalculateFitValue( run ) << " Measured value " << stat_info_[bin].GetMean() << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

ostream & operator << (ostream &o,const vector<Reco::StatInfo> &v)
{
  o << "Cells amount = " << v.size() << endl;
  o << " Cell      Mean     Sigma       WeightsSum " << endl;
  for( size_t i=0; i<v.size(); i++ )
  {
    const Reco::StatInfo &c=v[i];
    o << i << c.GetMean() << c.GetSigma() << c.GetEntries();
  }
  return o;
}

////////////////////////////////////////////////////////////////////////////////

istream & operator >> (istream &in,vector<Reco::StatInfo> &v)
{
  v.clear();
  string s;
  char dummy[111];

  getline(in,s);

  size_t len;
  sscanf(s.c_str(),"%s %s %c %zu",dummy,dummy,dummy,&len);

  getline(in,s);

  while( getline(in,s) )
  {
    unsigned i;
    float weight,mean,sigma;
    sscanf(s.c_str(),"%u %g %g %g",&i,&mean,&sigma,&weight);
    if( i!=v.size() )
      throw Reco::Exception("istream & operator >> (istream &in,vector<Reco::StatInfo> &v):  wrong cell number %d. I expect %d.",i,v.size());
    v.push_back(Reco::StatInfo(weight,mean,sigma));

  }

  if( v.size()!=len )
    throw Reco::Exception("istream & operator >> (istream &in,vector<Reco::StatInfo> &v):  bad size: I expect %d but I have read %d",
                     len,v.size());
  return in;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco

////////////////////////////////////////////////////////////////////////////////

ostream & operator << (ostream &o,const vector<Reco::StatInfo> &v)
{
  o << " StatInfo size = " << v.size() << endl;
  o << " index      Mean     Sigma       WeightsSum " << endl;
  for( size_t i=0; i<v.size(); i++ )
  {
    const Reco::StatInfo &c=v[i];
    o << i << c.GetMean() << c.GetSigma() << c.GetEntries();
  }
  return o;
}

////////////////////////////////////////////////////////////////////////////////

istream & operator >> (istream &in,vector<Reco::StatInfo> &v)
{
  v.clear();
  string s;
  char dummy[111];

  getline(in,s);

  size_t len;
  sscanf(s.c_str(),"%s %s %c %zu",dummy,dummy,dummy,&len);

  getline(in,s);

  while( getline(in,s) )
  {
    unsigned i;
    float weight,mean,sigma;
    sscanf(s.c_str(),"%u %g %g %g",&i,&mean,&sigma,&weight);
    if( i!=v.size() )
      throw Reco::Exception("istream & operator >> (istream &in,vector<Reco::StatInfo> &v):  wrong cell number %d. I expect %d.",i,v.size());
    v.push_back(Reco::StatInfo(weight,mean,sigma));

  }

  if( v.size()!=len )
    throw Reco::Exception("istream & operator >> (istream &in,vector<Reco::StatInfo> &v):  bad size: I expect %d but I have read %d",
                     len,v.size());
  return in;
}

////////////////////////////////////////////////////////////////////////////////
