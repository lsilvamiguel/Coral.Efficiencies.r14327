
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iostream>

#include "StoreDeltaADC.h"


#define CDB_LINE_MAX 132

using namespace std;
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StoreDeltaADC::StoreDeltaADC ( void ) 
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int StoreDeltaADC::GetDelta ( const std::pair<size_t,size_t> &src_id ) const
{
  std::map <  std::pair<size_t,size_t>, size_t > ::const_iterator it = srcid_delta_.find(src_id);
  if( it !=  srcid_delta_.end() )
    return    (int)it->second;
  else
    return -1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool StoreDeltaADC::UpdateDelta ( const std::pair<size_t,size_t> &src_id, size_t delta )
{
  std::map <  std::pair<size_t,size_t>, size_t > ::iterator it = srcid_delta_.find(src_id);
  if( it !=  srcid_delta_.end() )
  {
    if( it->second > delta )
    { 
      it->second = delta;
      return true;
    }
    return false;  
  }  
  else
  {
   srcid_delta_[src_id]=delta ;
   return true;
  }  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair <int,int>  StoreDeltaADC::MiniMaxSADC(  const std::vector<unsigned short>& samples )
{
  int min = 10000; 
  int max = -10000;
  for( size_t is=0; is<samples.size(); is++ )
  {
    if( samples[is] < min ) min = samples[is];
    if( samples[is] > max ) max = samples[is];
  }
  return std::pair<int,int>(min,max);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int StoreDeltaADC::InputDeltaADC(const std::string &s)
{
// //   bool debug = true;
//   bool debug = false;
//   if( debug ) 
//   {
//      cout <<" StoreDelatADC::InputDeltaADCD " << c_->GetName() << " debug " << endl;
//   }
// 
//   istringstream is(s);
//   char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
//   string str;
// 
// //  Read line of comments
//   getline(is,str);
//   getline(is,str);
//   getline(is,str);
// 
// //  Read calibration
//   unsigned int not_matched_cells_cnt = 0;
//   while( getline(is,str) )
//   {
//     int icell,stat;
//     float m0,m1,m2;
// 
//     assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
//     int ret = sscanf(str.c_str(), " %d %g %g %g %d %s", &icell, &m0, &m1, &m2, &stat, cellname);
//     assert( ret == 6 );
//     int jcell = c_->FindCellByName( cellname );
// 
//     if( icell != jcell ) {
//       icell = jcell;
//       not_matched_cells_cnt++;
//       if (not_matched_cells_cnt==1) {
// 	cerr << "Notice: " <<  __func__ << " " << c_->GetName() << " "
// 	     << "cell id and cell name do not match.  (The calibration file "
// 	     << "was produced with a different geometry description, which "
// 	     << "might be an indication that you're using the wrong file.)"
// 	     << endl;
//       }
//     }
//       
//     if(icell >=0 && icell < (int)NCells() )
//       c_->SetCorrTimeInSpillLED(icell,m0/10000.,m1/10000.,m2/10000., double(stat) );
//   }
//   if( debug ) 
//   {
//      cout <<" StoreDelatADC::InputDeltaADCD " << c_->GetName() << " debug OK " << endl;
//   }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int StoreDeltaADC::OutputDeltaADC( std::string &s ) const
{
//  bool debug = false;
  bool debug = true;
  int nsources = 0;
  char o[5000000];
  o[0]=0;
  sprintf(o,"### \n");
  sprintf(o+strlen(o),"  SrcID  ChipID  Delta \n");
  if( debug ) std::cout <<" SrcID ChipID Delta " << std::endl;

  std::map <  std::pair<size_t,size_t>, size_t > ::const_iterator it;
  for( it = srcid_delta_.begin(); it != srcid_delta_.end(); it++ )
  {
     sprintf(o+strlen(o)," %zu  %zu %zu  \n",it->first.first, it->first.second, it->second );
     if( debug ) std::cout <<" " << nsources << " " <<it->first.first << " " << it->first.second << " " <<it->second << std::endl;
     nsources++;
  }
  string ss(o);
//  s += ss;
  s = ss;
  if( debug )
  {
    std::cout << " StoreDeltaADC::OutputDeltaADC " << nsources << "  debug finished OK " << std::endl;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

