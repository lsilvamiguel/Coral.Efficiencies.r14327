// $Id: Utils.cc,v 1.6 2010/03/16 08:33:04 suhl Exp $

/*!
  \file    Utils.cc
  \brief   some utilities used everywhere
  \author  Hugo Pereira
  \version $Revision: 1.6 $
  \date    $Date: 2010/03/16 08:33:04 $
*/
#include "Utils.h"
#include "fdstream.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>

#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;

//_________________________________
ClassImp(Utils)

void Utils::DumpMatrix( TMatrix m, std::ostream &out )
{
    char* text = new char[200];
  for( int i=0; i< m.GetNrows(); i++ ) {
    for( int j=0; j< m.GetNcols(); j++ ) 
    snprintf(text, 200, "  %10f", m(i,j) );
    out<<text;
    out << endl;
  }
  out << endl;
  
  return;
}

//_________________________________
TMatrix Utils::InvertMatrix( TMatrix m )
{
  TMatrix out( 2,2 );
  if( !( m.GetNrows() == 2 && m.GetNcols() == 2 ) ) {
    cerr << "Utils::InvertMatrix - Matrix is not 2x2. Aborted.\n";
    return out;
  }
  
  double det = m(0,0)*m(1,1)-m(0,1)*m(1,0);
  if( fabs( det ) < 1e-6 ) {
    cerr << "Utils::InvertMatrix - Determinant to small. Aborted.\n";
    return out;
  }
  
  out(0,0) = m(1,1)/det;  out(0,1) = -m(0,1)/det;
  out(1,0) = -m(1,0)/det; out(1,1) = m(0,0)/det;  
  return out;
}

//________________________________________
int Utils::GetColorCode( const char* TBName )
{
  string TB = string( TBName ).substr(0,2);
  if( TB == "SI" ) return 2;
  if( TB == "FI" ) return 1;
  if( TB == "MM" ) return 4;
  if( TB == "GM" ) return 4;
  if( TB == "DC" ) return 2;
  if( TB == "ST" ) return 6;
  if( TB == "DW" ) return 5;
  if( TB == "PA" || TB == "PB" || TB == "PS" ) return 3;
  if( TB == "MA" || TB == "MB" ) return 1;
  if( TB == "HI" || TB == "HO" || TB == "HL" || TB == "HM" ) return 2;

  return 1;
}

#ifndef __CINT__ 
//________________________________________
TH1* Utils::Get( vector< string > files, const char* name )
{
  TH1* histo = 0;
  for( unsigned int i=0; i<files.size(); i++ ) {
    TFile* f = new TFile( files[i].c_str() );
    if( !f->IsOpen() ) continue;
    TH1* h = (TH1*) gROOT->FindObject(name);
    if( h ) {
      if( !i ) histo = (TH1*) h->Clone();
      else histo->Add( h );
    } else cout << "Utils::Get - histogram \"" << name << "\" not found in file \"" << files[i] << ".\n";
  }
  
  if( !histo ) cout << "Utils::Get - no match for \"" << name << "\".\n";
  return histo;
}
    
//________________________________________
TH2* Utils::Get2D( vector< string > files, const char* name )
{
  TH2* histo = 0;
  for( unsigned int i=0; i<files.size(); i++ ) {
    TFile* f = new TFile( files[i].c_str() );
    if( !f->IsOpen() ) continue;
    TH2* h = (TH2*) gROOT->FindObject(name);
    if( h ) {
      if( !i ) histo = (TH2*) h->Clone();
      else histo->Add( h );
    } else cout << "Utils::Get2D - histogram \"" << name << "\" not found in file \"" << files[i] << ".\n";
  }
  
  if( !histo ) cout << "Utils::Get2D - no match for \"" << name << "\".\n";
  return histo;
}
#endif
//______________________________________________________
double Utils::HDiv(TH1* h1, TH1* h2, TH1* h3) 
{

  unsigned int n1 = h1->GetNbinsX();
  unsigned int n2 = h2->GetNbinsX();
  unsigned int n3 = h3->GetNbinsX();

  if(!(n1 == n2 && n2 == n3)){
    cout << "Utils::HDiv - ERROR: Different number of bins.\n";
    cout << "   " << n1 << ", " << n2 << ", " << n3 << endl;
    return 0;
  } 
  return HDiv( h1, h2, h3, 1, n1 ); 
}
    
//______________________________________________________
double Utils::HDiv(TH2* h1, TH2* h2, TH2* h3) 
{
  
  unsigned int nx1 = h1->GetNbinsX();
  unsigned int nx2 = h2->GetNbinsX();
  unsigned int nx3 = h3->GetNbinsX();
  
  unsigned int ny1 = h1->GetNbinsY();
  unsigned int ny2 = h2->GetNbinsY();
  unsigned int ny3 = h3->GetNbinsY();

  if(!( nx1 == nx2 && 
     nx2 == nx3 && 
     ny1==ny2 && 
     ny2==ny3 )){
    cout << "Utils::HIntegrate - ERROR: Different number of bins.\n";
    cout << "   " << nx1 << ", " << nx2 << ", " << nx3 << endl;
    cout << "   " << ny1 << ", " << ny2 << ", " << ny3 << endl;
    return 0;
  } 

  double b1,b2,b3,e3;
  for(unsigned int i = 1; i < nx1+1; i++)
  for(unsigned int j = 1; j < ny1+1; j++) {   
    //loop over bins
    b1= (double) h1->GetBinContent(i,j);
    b2= (double) h2->GetBinContent(i,j);

    if(b2 != 0){
      b3 = b1/b2;
      e3=sqrt(b3*(1-b3))/sqrt(b2);
    } else {
      b3 = 0.; 
      e3 = 0.;
    }
    h3->SetBinContent(i,j,b3);
    h3->SetBinError  (i,j,e3);
  }
  return ( (double) h1->GetEntries() / (double) h2->GetEntries() );
} 

//______________________________________________________
double Utils::HDiv(TH1* h1, TH1* h2, TH1* h3, unsigned int i1, unsigned int i2) 
{
  
  unsigned int n1 = h1->GetNbinsX();
  unsigned int n2 = h2->GetNbinsX();
  unsigned int n3 = h3->GetNbinsX();

  if(!(n1 == n2 && n2 == n3)){
    cout << "Utils::HDiv - ERROR: Different number of bins.\n"; 
    cout << "   " << n1 << ", " << n2 << ", " << n3 << endl;
    return 0;
  } 

  double b1,b2,b3,e3;
  for(unsigned int i = 1; i < n1+1; i++) {
    if( i >= i1 && i< i2+1 ) {
   
      //loop over bins
      b1= (double) h1->GetBinContent(i);
      b2= (double) h2->GetBinContent(i);

      if(b2 != 0){
        b3 = b1/b2;
        e3=sqrt(b3*(1-b3))/sqrt(b2);
      } else {
        b3 = 0.; 
        e3 = 0.;
      }
      h3->SetBinContent(i,b3);
      h3->SetBinError  (i,e3);
    } else {
      h3->SetBinContent(i,0);
      h3->SetBinError  (i,0);
    }
  }
  return ( (double) h1->GetEntries() / (double) h2->GetEntries() );
}


//________________________________________________________________________________
TCanvas* Utils::DivideCanvas( TCanvas *cv, int nPads )
{  
  if( !cv ) return 0;
  
  if( nPads<=1 ) return cv;
  int sqnP = (unsigned int) (sqrt( (double)nPads ));
  int nC = sqnP;
  int nL = sqnP;
  
  while( nC*nL < nPads ) 
  if( nC < nL ) nC++;
  else nL++;
  
  cv->Divide( nC, nL, 0.005, 0.005, 0 );
  return cv;
}


#ifndef __CINT__
//____________________________________________________________________
string Utils::CopyToLocal( string filename )
{
  size_t start = filename.rfind("/"); 
  if( start == string::npos ) return string( filename );
  if( start == filename.size()-1 ) return string( filename );
  string out = filename.substr( start+1, filename.size() );
  string command = string( "cp " ) + filename+ " "+ out;
  cout << "Utils::CopyToLocal - INFO: " << filename << " -> " << out << endl;
  system( command.c_str() );
  return out;
} 
  
//____________________________________________________________________
string Utils::GetTimeStamp( string format )
{
  string command = string("date +")+format;
  FILE* fp;
  fp=popen(command.c_str(),"r");
  boost::fdistream in(fileno(fp));
  string out;
  in >> out;
  return out;
}

//____________________________________________________________________
vector<string> Utils::GetFiles( string fileSelection )
{
  vector<string> out;
  string command = string("ls -1 ") + fileSelection;
  FILE* fp;
  fp=popen(command.c_str(),"r");
  boost::fdistream in(fileno(fp));
 
  if( !in ) {
    cout << "Utils::GetFiles - ERROR: cannot execute command \"" << command << "\".\n";
    return out;
  } 
  
  char line[128];
  while( !in.eof() ) {
    in.getline( line, sizeof( line ) );
    
    istringstream s(line);
    string name; s >> name;
    if( !name.size() ) continue;
    
    if( access( name.c_str(), R_OK ) ){
      cout << "Utils::GetFiles - INFO: cannot access file \"" << name << "\".\n";
      continue;
    }
    
    out.push_back( name );
  }
    
  return out;
}
	  
//_________________________________
string Utils::MakeBackup( string filename )
{
  if( access( filename.c_str(), R_OK ) ) return filename;
  string backup;
  unsigned int ver=0;
  do{
    char* buf = new char[128];
    sprintf( buf, "%s.%i", filename.c_str(), ver );
    backup = string( buf );
    delete buf;
    ver++;
  } while  ( !access( backup.c_str(), R_OK ) );
  char* buf = new char[256];
  sprintf( buf,"cp %s %s", filename.c_str(), backup.c_str() );
  system( buf );
  delete buf;
  cout << "Utils::MakeBackup - file \""
    << backup 
    << "\" created.\n";
  return backup;
}

//_________________________________
bool Utils::Match( string s1, string s2 )
{
  
  //=== check sizes
  if( s2.size() > s1.size() ) return false;
  
  //=== check letter/letter
  bool accept = true;
  for( unsigned int i=0; i < s2.size(); i++ ) 
  if( s2[i] != '*' && s2[i] != s1[i] ) {
    accept = false;
    break;
  }
  
  return accept;
}

//_________________________________
string Utils::RemovePath( string filename )
{
  size_t found = filename.rfind("/");
  if( found==string::npos ) return filename;
  if( found==filename.size()-1 ) return "";
  return filename.substr( found+1, filename.size()-found-1);
}
   
#endif //!<__CINT__
