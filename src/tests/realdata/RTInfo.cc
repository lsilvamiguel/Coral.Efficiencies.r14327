// $Id: RTInfo.cc,v 1.3 2008/04/10 13:43:20 rgazda Exp $

/*!
   \file    RTInfo.cc
   \brief   to store detector associated RT relation informations
   \author  Hugo Pereira
   \version $Revision: 1.3 $
   \date    $Date: 2008/04/10 13:43:20 $
*/

#include "RTInfo.h"
#include "TMath.h"
#include "Fit.h"

#include <iostream>
#include <fstream>

using namespace std;

ClassImp( RTInfo )

//_________________________________________________________
RTInfo::RTInfo( string TBName )
{
  rtGrid_.clear();
  t0Par_.clear();
  rtPar_.clear();
  valid_start_="YYYY-MM-DD-hh:mm:ss";
  valid_stop_ ="YYYY-MM-DD-hh:mm:ss";
}

//_________________________________________________________________________________________
void RTInfo::_AddRTGridPoint( const double t, const double r, const double res  )
{
  RTGridPoint gp;
  gp.t=t;
  gp.r=r;
  gp.res=res;
  
  rtGrid_.push_back( gp );
  return;
}

//_________________________________________________________________________________________
bool RTInfo::SetValidity( const char* start, const char* stop )
{

  //=== Check validity_start format
  if( _CheckValidityFormat( start ) ) valid_start_ =string(start);
  else cout << "RTInfo::SetValidity - ERROR: wrong format for valid_start_ \n";  

  //=== Check validity_stop format
  if( _CheckValidityFormat( stop ) ) valid_stop_ =string(stop);
  else cout << "RTInfo::SetValidity - ERROR: wrong format for valid_stop_ \n";

  //=== dump and return
  cout << "RTInfo::SetValidity - INFO: validity for \"" 
    << TBName_ << "\" starts  at \""
    << valid_start_ << "\", stops  at \""
    << valid_stop_ << "\"\n";
  return true;
} 


//_________________________________________________________________________________________
void RTInfo::DumpT0Par( void )
{
  for( unsigned int iP = 0; iP < t0Par_.size(); iP++ )
  cout << " T0Par[" << iP <<"] = "<<t0Par_[iP]<< endl;  
}

//_________________________________________________________________________________________
void RTInfo::DumpRTPar( void )
{
  for( unsigned int iP = 0; iP < rtPar_.size(); iP++ )
  cout << " RTPar[" << iP <<"] = "<<rtPar_[iP]<< endl;  
}

//_________________________________________________________________________________________
void RTInfo::DumpRTGrid( void )
{
  for( unsigned int iP = 0; iP < rtGrid_.size(); iP++ )
  cout << " rtGrid[" << iP <<"].t   = "<< rtGrid_[iP].t
       << " rtGrid[" << iP <<"].r   = "<< rtGrid_[iP].r
       << " rtGrid[" << iP <<"].res = "<< rtGrid_[iP].res
       << endl;  
}


//_________________________________________________________________________________________
bool RTInfo::WriteRTParToDB( const bool useMyT0, const double myT0 )
{
  //=== check TBName and RT Parameters ===
  if( !( TBName_.size() && HasRTPar() ) ) {
    cout << "RTInfo::WriteRTParToDBFile - ERROR: either TBName or rtPar_ not set." << endl;
    return false;
  }

  //=== define and check outfile ===
  string outFile = TBName_+"~~start-"+valid_start_+"~~finish-"+valid_stop_;
  ofstream out( outFile.c_str(), ios::out );
  if( !out ) {
    cout << "RTInfo::WriteRTParToDBFile - ERROR: Cannot write to file \"" << outFile << "\"." << endl;
    return false;
  }
  
  cout << "RTInfo::WriteRTParToDBFile - INFO: DBFile is \"" << outFile << "\"." << endl;
  
  
  //=== Set RTRelations parameters
  vector<double> rtPar;
  string rtType;
  
  #ifdef _FIT1_
    // polynom fit 
    rtType = "RTFit1";
    for( unsigned int iP=0; iP < rtPar_.size(); iP++ ) rtPar.push_back( rtPar_[iP] );
  #else
    // hyperbolic tangeant fit
    // first parameter ( detector pitch ) is not written.
    rtType = "RTFit2";
    for( unsigned int iP=1; iP < rtPar_.size(); iP++ ) rtPar.push_back( rtPar_[iP] );
  #endif
  
  //=== write to file ===
  char* text = new char[200];
  
  if( useMyT0 ) { 
    // myT0 is used for t0 in the DB. it is written first. 
    // RTRelation first parameter is then properly corrected.
    // This because of redoundancy in the T0 an RTRelation parametrisation definitions.
    snprintf(text, 200, "%9.1f\n", myT0 );
    out<<text;
    out << rtType << " ";
    for( unsigned int iP=0; iP < rtPar.size(); iP++ ) {
	snprintf(text, 200, "%12.4e", (iP==0) ? (rtPar[iP]-myT0) : rtPar[iP] );
	out<<text;
    }
  } else {
    // First RTRelation fit Parameter is used as T0 and written first to the DB
    // it is then set to 0 when writting the RTRelation fit to the DB
    // This because of redoundancy in the T0 an RTRelation parametrisation definitions.
    snprintf(text, 200, "%9.1f\n", rtPar[0] );
    out<<text;
    out << rtType << "  ";
    for( unsigned int iP=0; iP < rtPar.size(); iP++ ) {
	snprintf(text, 200, "%12.4e", (iP==0) ? 0.0:rtPar[iP] );
	out<<text;
    }
  }
  
  SafeDelete(text);
  out.close();
  return true;
}

//_________________________________________________________________________________________
bool RTInfo::WriteRTGridToDB( const bool useMyT0, const double myT0 )
{
  if( !( HasTBName() && HasRTGrid() ) ) {
    cout << "RTInfo::WriteRTGridToDBFile - ERROR: either TBName or rtGrid_ not set." << endl;
    return false;
  }
  
  //=== define and check outfile ===
  string outFile=TBName_+"~~start-"+valid_start_+"~~finish-"+valid_stop_;
  ofstream out( outFile.c_str(), ios::out );
  if( !out ) {
    cout << "RTInfo::WriteRTGridToDBFile - INFO: Cannot write to file \""
         << outFile << "\"." 
         << endl;
    return false;
  }
  
  cout << "RTInfo::WriteRTGridToDBFile - INFO: DBFile is \""
       << outFile << "\"." 
       << endl;
  
  //=== write to file ===
  char* text = new char[200];
  
  if( useMyT0 ) {
    // myT0 is used for t0 in the DB. it is written first. 
    // it is then subtracted to all grid point time values,
    // This because of redoundancy in the T0 an RTRelation parametrisation definitions.
    snprintf(text, 200, "%9.1f\n", myT0 );
    out<<text;
    out << "RTGrid\n";
    for( unsigned int iG=0; iG < rtGrid_.size(); iG++ ) {
	snprintf(text, 200, "%9.4e %9.4e %9.4e\n",
		rtGrid_[iG].t - myT0, 
		rtGrid_[iG].r, 
		rtGrid_[iG].res 
		);
	out<<text;
    }
  } else {
    // First grid point time value is used as T0 and written first to the DB.
    // it is then subtracted to all grid point time values,
    // This because of redoundancy in the T0 an RTRelation parametrisation definitions.
    snprintf(text, 200, "%9.1f\n", rtGrid_[0].t );
    out<<text;
    out << "RTGrid\n";
    for( unsigned int iG=0; iG < rtGrid_.size(); iG++ ) {
	snprintf(text, 200, "%9.4e %9.4e %9.4e\n",
		rtGrid_[iG].t - rtGrid_[0].t, 
		rtGrid_[iG].r, 
		rtGrid_[iG].res 
		);
	out<<text;
    }
  }
  out.close();
  return true;
}

//_______________________________________________________________________________
double RTInfo::GetRfromT( const double t, bool &error )
{
  error = false;
  double r = 0;
  
  //=== RESULT FROM THE FIT ===
  if( HasRTPar() ) {
  
  
    #ifdef _FIT1_
    //=== polynomial version ===
    double tC = t-rtPar_[0];  // offset corrected time
    double v = 0;
    r = 0;
    for( unsigned int i = 1; i < rtPar_.size(); i++ ) {
      v += rtPar_[i]*i*pow( tC, (int)i-1);
      r += rtPar_[i]*pow( tC, (int)i);
    }
    if( v < 0 ) error = true;
  
    #else
    cout << "Fit Procedure for _FIT2_ not implemented.\n";
    error = true;
    r=0;
    #endif
    
  //=== RESULT FROM RAW GRID ===
  } else if( HasRTGrid() ) {
    unsigned int ig = 0;
	  double r0 = -1, r1 = -1;
	  double t0 = -1, t1 = -1;
  	while( ig < rtGrid_.size() && rtGrid_[ig].t < t ){
	    r0 = rtGrid_[ig].r;
	  	t0 = rtGrid_[ig].t;
	  	ig++;
	  }
	
	  if( ig == 0 ){ r = rtGrid_.front().r; error = true; }
	  else if( ig == rtGrid_.size() ){ r = rtGrid_.back().r; error=true; }
	  else {
	  	r1 = rtGrid_[ig].r;
	  	t1 = rtGrid_[ig].t;
	  	r = (t0 == t1 ) ? r1:( r0 + (t-t0)*(r1-r0)/(t1-t0) );
	  }

  //=== NOTHING KNOWN ABOUT RT RELATION ===
  } else {
    cout << "RTInfo::GetRfromT. Nothing known about rt relation.";
    r = 0;
    error = true;
  }

  //=== some checks ===
  if( r < 0 ) {
    r =0;
    error = true;
  } 
  
  return r;
}

//_______________________________________________________________________________
bool RTInfo::_CheckValidityFormat( const char* stamp )
{
  if( strlen( stamp ) < strlen("YYYY-MM-DD-hh:mm:ss") ) return false;
  if( stamp[4] != '-' ||
    stamp[7]  != '-' ||
    stamp[10] != '-' ||
    stamp[13] != ':' ||
    stamp[16] != ':' ) return false;
  return true;
}
