// $Id: Utils.cc,v 1.17 2002/09/16 17:23:30 hpereira Exp $
#include "Utils.h"
#include <unistd.h>
#include <procbuf.h>
#include <stdio.h>
#include <iostream>

#include <TPad.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TMatrix.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

//____________________________________________________________________
ClassImp( Utils )
char* Utils::GetTimeStamp( void )
{
  string buf( "date +%d_%m_%Y.%H_%M_%S" );
  static char stampS[56];
  
  procbuf inBuf;	
  inBuf.open( buf.c_str(), ios::in );
	istream in( &inBuf );
  string out;
  in >> out;
  sprintf( stampS, "%s", out.c_str() );
  return stampS;
}
	  
//______________________________________________________________________
string Utils::Convert(string sIn, string c1, string c2) 
{
	string sOut="";
	for( unsigned int i = 0; i < sIn.size(); i++ )
	if( sIn.substr( i, c1.size() ) == c1 ) {
	  sOut += c2;
		i+=c1.size()-1;
	} else sOut += sIn.substr( i, 1 );
	
	return sOut;
}
 
//______________________________________________________
double Utils::HDiv(TH1* h1, TH1* h2, TH1D* h3) 
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
double Utils::HDiv(TH1* h1, TH1* h2, TH1D* h3, unsigned int i1, unsigned int i2) 
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
    
//______________________________________________________
double Utils::HDiv(TH2* h1, TH2* h2, TH2D* h3) 
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
void Utils::HIntegrate(TH1S* h1, TH1S* h2) 
{
  
  unsigned int n1 = h1->GetNbinsX();
  unsigned int n2 = h2->GetNbinsX();
  if( n1 != n2 ){
    cout << "Utils::HIntegrate - ERROR: Different number of bins.\n"; 
    cout << "   " << n1 << " vs " << n2 << endl;
    return;
  }  

  double sum = 0;
  for( unsigned int i = 1; i < n1+1; i++ ) 
  h2->SetBinContent( i, ( sum+=h1->GetBinContent(i) ) );
  
  return;
}
       
//______________________________________________________
void Utils::HDerive(TH1S* h1, TH1S* h2) 
{
  
  unsigned int n1 = h1->GetNbinsX();
  unsigned int n2 = h2->GetNbinsX();
  if( n1 != n2 ){
    cout << "Utils::HIntegrate - ERROR: Different number of bins.\n"; 
    cout << "   " << n1 << " vs " << n2 << endl;
    return;
  }  

  double sum = 0;
  h2->SetBinContent( 1, 0 );
  for( unsigned int i = 2; i < n1+1; i++ ) 
  h2->SetBinContent( i, h1->GetBinContent(i)- h1->GetBinContent(i-1) );
  
  return;
}
      
//______________________________________________________
Utils::WStatus Utils::wStat_ = NO;
Utils::WStatus Utils::Writable( string fileName) 
{
  //=== check if ALWAYS was already replied
  if( wStat_ == ALWAYS_W ) return YES;
  
  //=== check if filename is new ===
  if( access( fileName.c_str(), R_OK ) ) return YES; 

  //=== check if NEVER was already replied
  if( wStat_ == NEVER_W ) {
    cout << "Utils::Writable - File \"" << fileName << "\" already exist.\n";
    return NO;
  }
  
  //=== 
  cout << "Utils::Writable - File \"" << fileName << "\" already exist. Overwrite\n";
  cout << "Utils::Writable - ([y] yes | [n] no | [a] always | [N] never | [e] exit) ?";
  
  while( 1 ) {
    char Rep = getchar();
    if( Rep == 'y' ) return YES;
    else if( Rep == 'n' ) { getchar(); return NO; }
    else if( Rep == 'e' ) { getchar(); return EXIT; }
    else if( Rep == 'a' ) { getchar(); wStat_ = ALWAYS_W; return YES; }
    else if( Rep == 'N' ) { getchar(); wStat_ = NEVER_W;  return NO; }
    else cout << "Utils::Writable - ERROR: bad key. Overwrite ?"; 
  }
  
  //=== This should not happen
  return NO;
}

//_________________________________
void Utils::DumpMatrix( TMatrix m, ostream &out=cout )
{
  for( int i=0; i< m.GetNrows(); i++ ) {
    for( int j=0; j< m.GetNcols(); j++ ) 
    out.form("  %10f", m(i,j) );
    out << endl;
  }
  out << endl;
  
  return;
}

//_______________________________________________________________________________
void Utils::SetVisibility( TCanvas *tcv )
{
  if( !tcv ) return;
  string Name( tcv->GetName() );
  
  TPad* pad;
  if ((pad = (TPad*)gROOT->FindObject(Name.c_str()))==0) {
    cout.form( "Utils::SetVisibility -ERROR: No \"%s\" TPad!\n",Name.c_str());
    return;
  }
  
  char* PadName = new char[128];
  tcv->cd();
  unsigned int nPads = 0;
  do { 
    nPads++; sprintf(PadName,"%s_%i",Name.c_str(),nPads );
  } while( gROOT->FindObject(PadName) );
  
  nPads--;
  cout << "Utils::SetVisibility -INFO: nPads= " << nPads << ".\n";
  double lsize, csiz;
  if (nPads>6) { lsize = 0.06; csiz = .065; }
  else { lsize = 0.05; csiz = .05; }
  
  for( unsigned int iPad = 0; iPad <= nPads; iPad++ ) {
    tcv->cd(iPad);
    TIter next(gPad->GetListOfPrimitives());
    
    TH1 *histo = 0; 
    TObject *obj; 
    int dim;
    
    while (obj = next())
    if (strncmp(obj->ClassName(),"TH",2)==0) {
  	  histo = (TH1*)obj;
  	  if (strncmp(obj->ClassName(),"TH1",3)==0) dim = 1;
  	  else dim = 2;
  	  break;
    }    
      
    if (!histo) {
      cout.form( "Utils::SetVisibility -ERROR: No histo in TPad \"%s\".\n",gPad->GetName());
      continue;
    }
    
    double xmax = histo->GetXaxis()->GetXmax(), xmin = histo->GetXaxis()->GetXmin();
    double ymax, ymin;
    if (dim==2) {
      ymax = histo->GetYaxis()->GetXmax(); 
      ymin = histo->GetYaxis()->GetXmin();
    } else { 
      ymax = histo->GetMaximum(); 
      ymin = 0;  
    }
    
    double dx = xmax-xmin, dy = ymax-ymin;
    printf("%d %s %f %f\n",iPad,histo->GetName(),xmin,xmax);
    
    // EDIT AXIS
    histo->GetXaxis()->SetNdivisions(505);
    histo->GetXaxis()->SetLabelSize(lsize);
    histo->GetYaxis()->SetNdivisions(505);
    histo->GetYaxis()->SetLabelSize(lsize);
  
    // EDIT TITLE
    TPaveText *title = (TPaveText*)gPad->GetPrimitive("title");
    Coord_t x1 = title->GetX1NDC();
    double x2 = title->GetX2NDC();
    dx = x2-x1;
    x2 += .25*dx; title->SetX2NDC(x2);
    Coord_t y1 = title->GetY1NDC();
    double y2 = title->GetY2NDC();
    dy = y2-y1;
    y1 -= dy; title->SetY1NDC(y1);
    title->SetTextSize(csiz);
    title->Draw(); gPad->Update();
    
    // EDIT PaveStats
    TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    x1 = st->GetX1NDC(), x2 = st->GetX2NDC(), dx = x2-x1;
    if (dim==1) x1 -= .15*dx;
    else        x1 -= .5*dx;
    st->SetX1NDC(x1);
    y1 = st->GetY1NDC(), y2 = st->GetY2NDC(), dy = y2-y1;
    if   (nPads>8) y1 -= 1.2*dy;
    else if (nPads>5) y1 -= .9*dy;
    else              y1 -= .5*dy;
    st->SetY1NDC(y1);
    st->SetTextSize(csiz);
    st->Draw(); gPad->Update();
  
    if (dim==2) {
      // EDIT PaveStats Fit
      TPaveText *st = (TPaveStats*)gPad->GetPrimitive("sRes");
      if (st) {
  	    x1 = st->GetX1NDC(), x2 = st->GetX2NDC(), dx = x2-x1;
    	  x1 -= .1*dx; st->SetX1NDC(x1);
    	  y1 = st->GetY1NDC(), y2 = st->GetY2NDC(), dy = y2-y1;
    	  y1 -= 1.8*dy; st->SetY1NDC(y1);
    	  y2 -= 1.3*dy; st->SetY2NDC(y2);
    	  st->SetTextSize(lsize);
    	  st->Draw(); gPad->Update();
      }
    }
  
    // Check next pad exists
    sprintf(PadName,"%s_%d",Name.c_str(),iPad+1);
    tcv->cd(); // Go back to main pad
  } 
}
   
   
