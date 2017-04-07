// --- Standard C/C++ library ---
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

// --- Internal files ----
#include "Calorimeter.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

std::string TimeConvert2string( const time_t  &time )
{
  char timechar[132];
  time_t  ltime = time;
  size_t nt = strftime(timechar, 132, "%a %b %d %T %Y", gmtime( &ltime));
  return std::string (timechar);
}

////////////////////////////////////////////////////////////////////////////////

pair<bool,vector<double> > Calorimeter::AlignmentInTime::GetPosition( const time_t  &time ) const
{
//   bool debug = true;
  bool debug = false;
  vector<double> position;
  bool ok= false;
// Configuration attempt
  size_t ip_xrel = 0;
  size_t ip_yrel = 1;
  size_t ip_xabs = 2;
  size_t ip_yabs = 3;
//  The default is the relative counters usage
//   size_t ip_x = ip_xrel;
//   size_t ip_y = ip_yrel;
  size_t ip_x = ip_xabs;
  size_t ip_y = ip_yabs;


  position.push_back(0.);
  position.push_back(0.);
  if( debug )
  {
    string timestr = TimeConvert2string( time );
    cout << " Calorimeter::AlignmentInTime::GetPosition time = " << timestr <<" MapSize " << map_time2position_.size() << endl;
  }
  if( map_time2position_.size() <= 0 ) return pair<bool,vector<double> >(ok,position);
  if( debug )
  {
    cout << " Compare from time " << TimeConvert2string(map_time2position_[0].first) <<
                          " to " << TimeConvert2string(map_time2position_.back().first) << endl;
  }
  for( size_t it=0; it < map_time2position_.size(); it++ )
  {
    if( map_time2position_[it].first >= time )
    {
      if( debug )
      {
        cout <<" At it = " << it <<" Got it ! " <<  TimeConvert2string( map_time2position_[it].first ) << endl;
        if( it > 0 ) cout <<" Previous " << it-1 << TimeConvert2string( map_time2position_[it-1].first ) << endl;
      }
      assert( map_time2position_[it].second.size() == 4 );
      if( it == 0 )
      {
        const pair < time_t, vector < double > > &p1 =  map_time2position_[it];
        position[0]= p1.second[ip_x];
        position[1]= p1.second[ip_y];
        ok = true;
        return pair<bool,vector<double> >(ok,position);
      }
      else
      {
        const pair < time_t, vector < double > > &p1 =  map_time2position_[it];
        const pair < time_t, vector < double > > &p0 =  map_time2position_[it-1];

        double dpx = p1.second[ip_x]-p0.second[ip_x];
        double dpy = p1.second[ip_y]-p0.second[ip_y];
        time_t dte = time - p0.first;
        time_t dt = p1.first - p0.first;
        double drt = 0.;
        if( dt > 0 ) drt = double(dte)/double(dt);
        dpx = dpx*drt;
        dpy = dpy*drt;
        if( debug )
        {
          cout << "  Time = " << TimeConvert2string(p0.first) <<" < " <<
                                      TimeConvert2string(time) << " < " <<
                                             TimeConvert2string(p1.first) << endl;
          cout << " Position is found (mm) x_rel=" << p1.second[0] << " y_rel=" << p1.second[1] <<
                                                  " x_abs=" << p1.second[2] << " y_abs=" << p1.second[3] << endl;
          cout << " Previous position  x_rel=" << p0.second[0] << " y_rel=" << p0.second[1] <<
                                                  " x_abs=" << p0.second[2] << " y_abs=" << p0.second[3] << endl;
          cout << " it map = " << it <<" rel time " << dte <<" time gate " << dt << " dpx = " << dpx <<" dpy " << dpy  <<  endl;
          cout << " Calculated position  x_rel=" << p1.second[0]-dpx << " y_rel=" << p1.second[1]-dpy <<
                                                  " x_abs=" << p1.second[2]-dpx << " y_abs=" << p1.second[3]-dpy << endl;
        }
//         position[0]= p1.second[0]-dpx;
//         position[1]= p1.second[1]-dpy;
        position[0]= p0.second[ip_x]+dpx;
        position[1]= p0.second[ip_y]+dpy;
        ok = true;
        return pair<bool,vector<double> >(ok,position);
      }
    }
  }
// pick-up last element
  const vector < double > &v = map_time2position_.back().second;
  if( v.size() == 4 )
  {
    if( debug ) cout << " Position is found (mm) x_rel=" << v[0] << " y_rel=" << v[1] <<
                                                " x_abs=" << v[2] << " y_abs=" << v[3] << endl;
    position[0]= v[ip_x];
    position[1]= v[ip_y];
    ok = true;
    return pair<bool,vector<double> >(ok,position);
  }
  else
  {
    cerr << " Error in UpdatePositionInTime wrong vector size = " << v.size() << endl;
    assert(false);
  }
  cerr << " Internal bug in UpdatePositionInTime should never be here " << endl;
  assert(false);
}

/////////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::UpdatePositionInTime( const time_t  &time )
{
//   bool debug = true;
  bool debug = false;

//   if( GetName() == "EC02P1__") debug = true;

  if( debug )
  {
//     char timechar[132];
//     time_t  ltime = time;
//     size_t nt = strftime(timechar, 132, "%a %b %d %T %Y", gmtime( &ltime));
//     string timestr(timechar);
    string timestr = TimeConvert2string( time );
    cout << " Try to UpdatePositionInTime for calorimeter " << GetName() << endl;
    cout << " Event time = " << timestr << " x_at_startrun = " << position_atstartofrun_[0] <<
                                                        " y_at_startrun = " << position_atstartofrun_[1] << endl;
  }

  pair<bool,vector<double> > position_in_time = alignment_in_time_.GetPosition( time );
  if( !position_in_time.first ) return false;
  position[0]= position_in_time.second[0];
  position[1]= position_in_time.second[1];
  position_correction_[0]=position[0]-position_atstartofrun_[0];
  position_correction_[1]=position[1]-position_atstartofrun_[1];
  if( debug )
  {
    cout << " x_updated = " << position[0] <<" y_updated = " << position[1] << endl;
  }
  return true;

//   if( map_time2position_.size() <= 0 ) return false;
//   if( debug )
//   {
//     cout << " Compare from time " << map_time2position_[0].first << " to " << map_time2position_.back().first << endl;
//   }
//   for( size_t it=0; it < map_time2position_.size(); it++ )
//   {
//     if( map_time2position_[it].first >= time )
//     {
//       assert( map_time2position_[it].second.size() == 4 );
//       if( it == 0 )
//       {
//         const pair < time_t, vector < double > > &p1 =  map_time2position_[it];
//         position[0]= p1.second[0];
//         position[1]= p1.second[1];
//         position_correction_[0]=position[0]-position_atstartofrun_[0];
//         position_correction_[1]=position[1]-position_atstartofrun_[1];
//         return true;
//       }
//       else
//       {
//         const pair < time_t, vector < double > > &p1 =  map_time2position_[it];
//         const pair < time_t, vector < double > > &p0 =  map_time2position_[it-1];
//
//         double dpx = p1.second[0]-p0.second[0];
//         double dpy = p1.second[1]-p0.second[1];
//         time_t dte = time - p0.first;
//         time_t dt = p1.first - p0.first;
//         double drt = 0.;
//         if( dt > 0 ) drt = double(dte)/double(dt);
//         dpx = dpx*drt;
//         dpy = dpy*drt;
//         if( debug )
//         {
//           cout << "  Time = " << p0.first <<" < " << time << " < " << p1.first << endl;
//           cout << " Position is found (mm) x_rel=" << p1.second[0] << " y_rel=" << p1.second[1] <<
//                                                   " x_abs=" << p1.second[2] << " y_abs=" << p1.second[3] << endl;
//           cout << " Previous position  x_rel=" << p0.second[0] << " y_rel=" << p0.second[1] <<
//                                                   " x_abs=" << p0.second[2] << " y_abs=" << p0.second[3] << endl;
//           cout << " it map = " << it <<" rel time " << dte <<" time gate " << dt << " dpx = " << dpx <<" dpy " << dpy  <<  endl;
//           cout << " Calculated position  x_rel=" << p1.second[0]-dpx << " y_rel=" << p1.second[1]-dpy <<
//                                                   " x_abs=" << p1.second[2]-dpx << " y_abs=" << p1.second[3]-dpy << endl;
//         }
//         position[0]= p1.second[0]-dpx;
//         position[1]= p1.second[1]-dpy;
//         position_correction_[0]=position[0]-position_atstartofrun_[0];
//         position_correction_[1]=position[1]-position_atstartofrun_[1];
//         return true;
//       }
//     }
//   }
// // pick-up last element
//   const vector < double > &v = map_time2position_.back().second;
//   if( v.size() == 4 )
//   {
//     if( debug ) cout << " Position is found (mm) x_rel=" << v[0] << " y_rel=" << v[1] <<
//                                                 " x_abs=" << v[2] << " y_abs=" << v[3] << endl;
//     position[0]= v[0];
//     position[1]= v[1];
//     position_correction_[0]=position[0]-position_atstartofrun_[0];
//     position_correction_[1]=position[1]-position_atstartofrun_[1];
//     return true;
//   }
//   else
//   {
//     cerr << " Error in UpdatePositionInTime wrong vectro size = " << v.size() << endl;
//     return false;
//   }
//   cerr << " Internal bug in UpdatePositionInTime should never be here " << endl;
//   return false;
}

/////////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::UpdatePositionInBurst( int burst_number, const time_t  &time )
{
//   bool debug = true;
  bool debug = false;

//   if( GetName() == "EC02P1__") debug = true;

  if( debug )
  {
    char timechar[132];
    time_t  ltime = time;
    size_t nt = strftime(timechar, 132, "%a %b %d %T %Y", gmtime( &ltime));
    string timestr(timechar);
    cout << " UpdatePositionInBurst " << GetName() <<" Time " << timestr <<
                       " burst_number " << burst_number << " size " <<  map_burst2position_.size() << endl;

    cout << " current  x " << position[0] << " start run  x " << position_atstartofrun_[0] << endl;
    cout << " current  y " << position[1] << " start run  y " << position_atstartofrun_[1] << endl;
    cout << " current  z " << position[2] << " start run  z " << position_atstartofrun_[2] << endl;
  }

  if( map_burst2position_.empty() )
  {
//     cerr << " Error! UpdatePositionInBurst " << GetName() <<
//             "  map_burst2position_ is empty ! Please consider  UpdatePositionInTime instead " << endl;
    return UpdatePositionInTime( time );
  }
  std::map< int, std::vector< double> >::const_iterator end = map_burst2position_.end();
  std::map< int, std::vector< double> >::const_iterator e =  map_burst2position_.find(burst_number);
  if( e!=end )
  {
    const std::vector< double> &v = e->second;
    if( v.size() == 3 )
    {
      if( debug ) cout << " correct x =  " << v[0] << " y= " << v[1] << endl;
      position[0]=v[0];
      position[1]=v[1];
      position[2]=v[3];
      position_correction_[0]=position[0]-position_atstartofrun_[0];
      position_correction_[1]=position[1]-position_atstartofrun_[1];
      position_correction_[2]=position[2]-position_atstartofrun_[2];
   }
    else if ( v.size() == 2 )
    {
      if( debug ) cout << " correct x =  " << v[0] << " y= " << v[1] << endl;

      position[0]=v[0];
      position[1]=v[1];
      position_correction_[0]=position[0]-position_atstartofrun_[0];
      position_correction_[1]=position[1]-position_atstartofrun_[1];

    }
    else
    {
      if( debug ) cout << " PositionInBurst not found  " << endl;
      return false;
    }
    return true;
  }
  else
  {
    if( debug ) cout << " PositionInBurst not found  " << endl;
    return false;
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StorePositionInfo( const vector< pair< StatInfo, StatInfo> > &diff_in_bursts )
{
  bool debug = false;
  positionX_in_burst_new_.clear();
  positionY_in_burst_new_.clear();

  for( size_t n=0; n< diff_in_bursts.size(); n++ )
  {
    std::map< int, std::vector< double> >::const_iterator end = map_burst2position_.end();
    std::map< int, std::vector< double> >::const_iterator e =  map_burst2position_.find(n+1);
    double positionx = position[0];
    double positiony = position[1];
    double positionz = position[2];
    if( e!=end )
    {
      const std::vector< double> &v = e->second;
      if( v.size() == 3 )
      {
        if( debug ) cout << " correct x =  " << v[0] << " y= " << v[1] << endl;
       positionx=v[0];
       positiony=v[1];
       positionz=v[3];
//       position_correction_[0]=position[0]-position_atstartofrun_[0];
//       position_correction_[1]=position[1]-position_atstartofrun_[1];
//       position_correction_[2]=position[2]-position_atstartofrun_[2];
      }
      else if ( v.size() == 2 )
      {
        if( debug ) cout << " correct x =  " << v[0] << " y= " << v[1] << endl;

       positionx=v[0];
       positiony=v[1];
//       position_correction_[0]=position[0]-position_atstartofrun_[0];
//       position_correction_[1]=position[1]-position_atstartofrun_[1];

      }
    }


    double statx = diff_in_bursts[n].first.GetEntries();
    double meanx = diff_in_bursts[n].first.GetMean();
    double sigmax = diff_in_bursts[n].first.GetSigma();

    double staty = diff_in_bursts[n].second.GetEntries();
    double meany = diff_in_bursts[n].second.GetMean();
    double sigmay = diff_in_bursts[n].second.GetSigma();

    positionX_in_burst_new_.push_back( StatInfo(statx, positionx - meanx, sigmax) );
    positionY_in_burst_new_.push_back( StatInfo(staty, positiony - meany, sigmay) );
  }

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StorePositionInfo( int burst_number, const time_t  &time, double x_minus_xref, double y_minus_yref, double w )
{
  if( burst_number <= 0 )
  {
    cerr << " Wrong burst number " << burst_number << endl;
    return;
  }

  if( (int)positionX_in_burst_new_.size() < burst_number )
  {
    size_t nadd = burst_number - positionX_in_burst_new_.size();
    for( size_t n=0; n< nadd; n++ )
    {
      positionX_in_burst_new_.push_back( StatInfo(0.,0.,0.) );
      positionY_in_burst_new_.push_back( StatInfo(0.,0.,0.) );
    }
  }

  positionX_in_burst_new_[burst_number-1].Add( position[0]-x_minus_xref, w);
  positionY_in_burst_new_[burst_number-1].Add( position[1]-y_minus_yref, w);
}

////////////////////////////////////////////////////////////////////////////////

pair<bool,vector<double> > Calorimeter::GetPosition( const time_t  &time ) const
{
  return alignment_in_time_.GetPosition( time );
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::CorrectPosition( void )
{
  position[0] += position_correction_[0];
  position[1] += position_correction_[1];
  position[2] += position_correction_[2];
}

////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco

