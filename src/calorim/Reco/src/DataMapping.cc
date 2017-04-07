/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/DataMapping.cc,v $
   $Date: 2010/09/21 10:27:57 $
   $Revision: 1.14 $
   -------------------------------------------------------------------------
*/
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "Calorimeter.h"

#include "CellDataRaw.h"
#include "Exception.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ReadMap(const string &file_name)
{
  ifstream f(file_name.c_str());
  if( !f.is_open() )
    throw Exception("Calorimeter::ReadMap(): can not open file \"%s\"",file_name.c_str());

//  chan_to_cell.clear();

  string all;

  while( !f.eof() )
  {
    string line;
    getline(f,line);
    if( line.length()<=0 || line[0]=='#' )
      continue;
    all += line + ' ';
  }

  // Put all symbols to lower case..
  for( size_t i=0; i<all.length(); i++ )
    all[i]=tolower(all[i]);

  vector< pair<size_t,size_t> > range[2]; // channels, cells
  vector< int > step[2];                  // step in the range

  size_t what=0;  // 0 for channels and 1 for cells
  bool range_flag=false;
  bool step_flag=false;
  int offset[2];
  offset[0]=0;
  offset[1]=0;
  bool offset_flag[2];
  offset_flag[0]=false;
  offset_flag[1]=false;

  char *all_str = strdup( all.c_str() );
  const char delim[]=" 	";  // SPACE and TAB
  for( char *s=strtok(all_str,delim); s!=NULL; s=strtok(NULL,delim) )
  {
//     cout << s << endl;
    if( 0==strcmp(s,"channels") )
    {
      what=0;
      continue;
    }

    if( 0==strcmp(s,"channels_offset") )
    {
      offset_flag[0]=true;
      continue;
    }

    if( 0==strcmp(s,"cells") )
    {
      what=1;
      continue;
    }

    if( 0==strcmp(s,"cells_offset") )
    {
      offset_flag[1]=true;
      continue;
    }

    if( 0==strcmp(s,"step") )
    {
      step_flag=true;
      continue;
    }

    if( step_flag )
    {
      char *result=s;
      step[what].back()=strtol(s,&result,10);
      if( *result!=0 )
        throw Exception("Calorimeter::ReadMap(): file \"%s\".  Error in format.",file_name.c_str());
      step_flag=false;
      continue;
    }

    if( offset_flag[0] )
    {
      char *result=s;
      offset[0]=strtol(s,&result,10);
      if( *result!=0 )
        throw Exception("Calorimeter::ReadMap(): file \"%s\".  Error in format.",file_name.c_str());
      offset_flag[0]=false;
      continue;
    }

    if( offset_flag[1] )
    {
      char *result=s;
      offset[1]=strtol(s,&result,10);
      if( *result!=0 )
        throw Exception("Calorimeter::ReadMap(): file \"%s\".  Error in format.",file_name.c_str());
      offset_flag[1]=false;
      continue;
    }

    if( *s=='-' )
    {
      if( range_flag )
        throw Exception("Calorimeter::ReadMap(): file \"%s\".   Double range flag.",file_name.c_str());

      range_flag = true;
      s++;
      if( *s==0 )
        continue;  // go to next token
     // continue in the case of   -12333
    }

    char *p=strchr(s,'-');        // situation:   12-34  or     12- 34
    if( p!=NULL )
      *p=0;

    char *result=s;

    if( !range_flag )
    {
      int v =strtol(s,&result,10)+offset[what];
      range[what].push_back( pair<size_t,size_t>(v,v) );
      step[what].push_back( 0 );
      if( *result!=0 )
        throw Exception("Calorimeter::ReadMap(): file \"%s\".  Error in format.",file_name.c_str());
      if( p!=NULL )
      {
        range_flag=true;
        s = p+1;
        if( *s==0 )
          continue;   // this is situation     12- 34
        // this is situation 12-34
      }
    }

    if( range_flag )
    {
      assert(range[what].size()>0);
      range[what].back().second=strtol(s,&result,10)+offset[what];
      step[what].back()=1;
      if( *result!=0 )
        throw Exception("Calorimeter::ReadMap(): file \"%s\".  Error in format.",file_name.c_str());
      range_flag=false;
    }
  }
  free( all_str );

  for( size_t i=0; i<2; i++ )
  {
    for( size_t it=0; it!=range[i].size(); it++ )
      if( step[i][it]!=0 )
      {
        int r=range[i][it].second-range[i][it].first;
//         if(r == 0 )
//           throw Exception("Calorimeter::ReadMap(): file \"%s\"WRONG RANGE %4d - %4d or STEP %4d\n",
//                               file_name.c_str(),range[i][it].first,range[i][it].second,step[i][it]);
        int n=r/step[i][it];
        if(n<0)
          throw Exception("Calorimeter::ReadMap(): file \"%s\" WRONG RANGE %4d - %4d or STEP %4d\n",
                               file_name.c_str(),range[i][it].first,range[i][it].second,step[i][it]);
        if(r-n*step[i][it] != 0)
          throw Exception("Calorimeter::ReadMap(): file \"%s\" WRONG RANGE %4d - %4d or STEP %4d\n",
                               file_name.c_str(),range[i][it].first,range[i][it].second,step[i][it]);
      }
      else
      {
//        if(range[i][it].first != range[i][it].second )
          throw Exception("Calorimeter::ReadMap(): file \"%s\" WRONG RANGE %4d - %4d or STEP %4d\n",
                               file_name.c_str(),range[i][it].first,range[i][it].second,step[i][it]);
      }
  }

//   for( size_t i=0; i<2; i++ )
//   {
//     cout << "This is " << (i==0?"channels":"cells") << " list\n";
//     for( size_t it=0; it!=range[i].size(); it++ )
//       if( step[i][it]!=0 )
//         printf("%4d - %4d step %4d\n",range[i][it].first,range[i][it].second,step[i][it]);
//       else
//         printf("%4d\n",range[i][it].first);
//   }

  vector< size_t > map_cells,map_channels; // channels, cells
  for( size_t it=0; it!=range[0].size(); it++ )
  {
    if(step[0][it] == 0 )
    {
      map_channels.push_back(range[0][it].first);
    }
    else
    {
      for( size_t it1=range[0][it].first; it1!=range[0][it].second; it1 += step[0][it])
      {
        map_channels.push_back(it1);
      }
      map_channels.push_back(range[0][it].second);
    }
  }

  for( size_t it=0; it!=range[1].size(); it++ )
  {
    if(step[1][it] == 0 )
    {
      map_cells.push_back(range[1][it].first);
    }
    else
    {
      for( size_t it1=range[1][it].first; it1!=range[1][it].second; it1 += step[1][it])
      {
        map_cells.push_back(it1);
      }
      map_cells.push_back(range[1][it].second);
    }
  }

  if( map_channels.size() != map_cells.size() )
      throw Exception("Calorimeter::ReadMap(): file \"%s\" number of Channels %4d But number of Cells %4d \n",
                             file_name.c_str(),map_channels.size(),map_cells.size());

  for( size_t it=0; it < NCells(); it++ )
  {

  }

  for( size_t it=0; it < NCells(); it++ )
  {
    size_t n=0;
    for( size_t  it1=0; it1!=map_cells.size(); it1++)
    {
      if( it == map_cells[it1])
      {
        map___chan_to_cell.insert( pair<size_t,size_t>(map_channels[it1],it) );
        n++;
      }
    }
    if(n!=1) cout << " wrong Mapping " << n << endl;
    if( n==0 )
      throw Exception("Calorimeter::ReadMap(): file \"%s\" Mapping of Cell %4d is undeclared \n",
                                                                           file_name.c_str(),it);
    if( n>1 )
      throw Exception("Calorimeter::ReadMap(): file \"%s\" Mapping of Cell %4d is multiply declared \n",
                                                                           file_name.c_str(),it);
  }

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SetMap(const map<size_t,size_t> &chan_to_cell)
{
  size_t n=0; // counter for bad elements of the map.

  for( map<size_t,size_t>::const_iterator it=chan_to_cell.begin(); it!=chan_to_cell.end(); it++ )
    if( it->second>=NCells() )
    {
      n++;
      cerr << "Calorimeter::SetMap(): There are " << NCells() <<
              " cells, but in your map there is channel=" << it->first <<
              " with cell=" << it->second << endl;
    }

  if( n>0 )
    throw Exception("Calorimeter::SetMap(): there are %d bad elements in the map.",n);
  else
    map___chan_to_cell = chan_to_cell;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::Remap(const vector<CellDataRaw>& data_in, vector<CellDataRaw>& data_out) const
{
  if( &data_in==&data_out )
    throw Exception("Calorimter::Remap():  data_in and data_out are the same containers.");

  size_t n=0; // counter for converted cells
  map<size_t,size_t>::const_iterator end=map___chan_to_cell.end();

  for( vector<CellDataRaw>::const_iterator it=data_in.begin(); it!=data_in.end(); it++ )
  {
    map<size_t,size_t>::const_iterator e = map___chan_to_cell.find(it->GetCellIdx());
    if( e!=end )
    {
      if( !(e->second < NCells()) ) assert(false);
      data_out.push_back(CellDataRaw(e->second, *it));
    }
    else
    {
      n++;
      if( options.print_unmaped_data )
        cerr << "Calorimter::Remap():  unknown address " << it->GetCellIdx() << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::RemapCalibInternal( void )
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::RemapInternal " << GetName() << " map2remap_.size() = " << map2remap_.size() << endl;
  if( map2remap_.size() == 0 ) return;

  map<int,int>::const_iterator end = map2remap_.end();
  vector<CellDataRaw> signals_new;

  int nremap = 0;

  for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++) {
    map<int,int>::const_iterator e = map2remap_.find( it->GetCellIdx() );
    if( e!=end )
    {
      if( e->second < 0 || e->second > (int)NCells() ) continue;
      nremap++;
      signals_new.push_back( CellDataRaw(e->second, *it) );
    }
    else
    {
      signals_new.push_back( *it );
    }
  }
  if( debug ) cout << " Calorimeter::RemapInternal " << GetName() <<
                       " remmaped " << nremap << " out of " << signals_new.size() <<  endl;
  signals_ = signals_new;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::RemapCalibCellInfo(CellInfoType what, CalibTimeType when)
{
  if( map2remap_.size() == 0 ) return;
  if( what < 0 || what > (int)CellInfoTypeSize ) return;
  if( when < 0 || when > (int)CalibTimeTypeSize ) return;

  list < pair <int, StatInfo> > s;
  for( map<int,int>::const_iterator it=map2remap_.begin(); it != map2remap_.end(); it++ )
  {
    if( it->first >= 0 &&  it->second >= 0 )
      s.push_back( pair <int, StatInfo> (it->first, GetCellInfo(what, when, it->second) ) );
  }
  for( list< pair <int, StatInfo> >::const_iterator it=s.begin(); it != s.end(); it++ )
  {
     SetCellInfo(what, when, it->first, it->second );
  }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco
