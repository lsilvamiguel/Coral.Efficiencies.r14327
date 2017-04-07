// $Id: CheckTracks.cc,v 1.34 2011/02/03 14:11:27 ybedfer Exp $

/*!
  \file    CheckTracks.cc
  \brief   To check coral output alignment tree
  \author  Hugo Pereira
  \version $Revision: 1.34 $
  \date    $Date: 2011/02/03 14:11:27 $
  */

#include "CheckTracks.h"
#include "Macro.h"
#include "Tracks.h"
#include "DetFileManager.h"
#include "DetectorInfo.h"
#include "Utils.h"
#include "TH2Fit.h"
#include "Fit.h"
#include "Opt.h"

#include <unistd.h>
#include <stdio.h>

#include <vector>
#include <sstream>

#include <TChain.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TH3.h>
#include <TPaveText.h>
#include <TText.h>
#include <TPaveStats.h>
#include <TMatrix.h>
#include <TGraphErrors.h>
#include <TCut.h>
#include <TProfile.h>

using namespace std;

//________________________________________
ClassImp(CheckTracks)
CheckTracks::CheckTracks(const char* trackfileselection, const char* detectorfile, bool magnets_on ):
  TObject(),
  isBatch_( false ),
  glCut_(""),
  color_(1),
  ps_(0),
  cv_(0),
  df_(0),
  tracks_(0),
  detectorFile_( "" )
{
  if( strlen( detectorfile ) ) LoadDetectorFile( detectorfile );
  if( strlen( trackfileselection ) ) LoadTracks( trackfileselection, magnets_on ); 
  gROOT->SetStyle("Plain");
  //  gStyle->SetOptDate();
  //  gStyle->SetOptFile();
}

//________________________________________
DetFileManager* CheckTracks::LoadDetectorFile( const char* name )
{
    if( df_ ) SafeDelete(df_);
    df_ = new DetFileManager( name );
    detectorFile_ = string(name);
    return df_;
}

//________________________________________
Tracks* CheckTracks::LoadTracks( const char* fileselection, bool magnets_on )
{
    if( tracks_ ) SafeDelete(tracks_);
    cout << "CheckTracks::LoadTracks - magnets: " << ((magnets_on)?"on":"off") << ".\n";
    tracks_ = new Tracks( magnets_on );  
    tracks_->AddToChain( fileselection );
    trackFileSelection_.clear();
    trackFiles_.clear();

    trackFileSelection_.push_back(fileselection);
    vector< string > files = Utils::GetFiles( fileselection );
    for( unsigned i = 0; i < files.size(); i++ ) trackFiles_.push_back( files[i] );

    return tracks_;
}
//________________________________________
Tracks* CheckTracks::AddToTracks( const char* fileselection )
{
    if( !tracks_ ) {
	cout << "CheckTracks::AddToTracks - Tracks object is not defined.\n";
	return 0;
    }
    tracks_->AddToChain( fileselection, isBatch_ );

    trackFileSelection_.push_back(fileselection);
    vector< string > files = Utils::GetFiles( fileselection );
    for( unsigned i = 0; i < files.size(); i++ ) trackFiles_.push_back( files[i] );
    return tracks_;
}

//________________________________________
TCanvas* CheckTracks::MakeCanvas( int width, int height )
{
    if( cv_ ) { cv_->Close(); SafeDelete(cv_); }
    cv_ = new TCanvas("CheckTracks","CheckTracks", width, height );
    return cv_;
}

//________________________________________
TPostScript* CheckTracks::OpenPS( const char* name, int type )
{
    if( ps_ ) { ps_->On(); ps_->Close(); SafeDelete(ps_); }
    ps_ = new TPostScript( name, type );
    ps_->Range(20, 29);
    ps_->Zone();
    ps_->Off();
    return ps_;
}

//________________________________________
void CheckTracks::ClosePS( void )
{   
    if( !ps_ ) return;
    _NewPage();
    ps_->On(); 
    ps_->Close(); 
    SafeDelete(ps_);
    ps_ = 0;
}    

//________________________________________________________________________________
void CheckTracks::MakeFrontPage( void )  
{ 
    if( !ps_ ) return;
    if( !cv_ ) MakeCanvas();
    _NewPage();
    ps_->On();
    ps_->SetTextSize(0.025);
    ps_->SetTextAlign(21);
    double step = 1.0/( trackFileSelection_.size()+3);
    for( unsigned int i=0; i < trackFileSelection_.size(); i++ )
	ps_->TextNDC(0.5, 1.0-step*(i+1), Utils::RemovePath(trackFileSelection_[i]).c_str());

    int i=trackFileSelection_.size();
    ps_->TextNDC(0.5, 1.0-step*(i+1), Utils::RemovePath(detectorFile_).c_str() ); i++;
    ps_->TextNDC(0.5, 1.0-step*(i+1), Utils::GetTimeStamp("%d/%m/%y-%H:%M:%S").c_str() );
    ps_->Off();
    return;
}

//________________________________________
void CheckTracks::Draw( 
	const char* var, 
	const char* cut, 
	const char* opt, 
	int nEnt,
	const char* XTitle,
	const char* YTitle )
{
    if( !tracks_) { cout << "CheckTracks::Draw - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();

    if( string( opt ).find("same")==string::npos ) {
	cv_->Clear();
	color_=1;
    } else color_++;

    Utils::DivideCanvas( cv_, 1 );
    cv_->cd(1);

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();

    char *buf = new char[64];
    sprintf( buf, "htemp_%i", color_ );
    string var_H = string(var)+">>"+buf;

    TH1 *h;
    TCut tCut( cut );
    if( glCut_.size() ) tCut+=TCut( glCut_.c_str() );

    if( (h = (TH1*) gDirectory->Get(buf)) ) SafeDelete(h);
    tracks_->GetChain()->Draw( var_H.c_str(), tCut, "goff", nEnt );
    h = (TH1*) gDirectory->Get(buf);

    cout << "CheckTracks::Draw - var \"" 
	<< var << "\": " 
	<< int((h)?h->GetEntries():-1)
				    << " entries.\n";

    if( h && h->GetEntries() ) {
	if( XTitle ) h->SetXTitle( XTitle );
	if( YTitle ) h->SetYTitle( YTitle );
	h->SetLineColor( color_ );
	h->Draw(opt);
    }
    SafeDelete( buf );

    _UpdateFromPad(gPad);
    _UpdateFromCanvas();
    if( string( opt ).find("same")==string::npos ) _MakeTitle( "Draw - %s:%s", var, (const char*) tCut );
    return;
}  

//________________________________________
void CheckTracks::DrawDet(  
	const char* var, 
	const char* detselection, 
	const char* cut, 
	const char* opt, 
	int nEnt,
	const char* XTitle,
	const char* YTitle )
{
    if( !df_ )     { cout << "CheckTracks::DrawDet - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawDet - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetChain()->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawDet - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();
    for( unsigned int i=0; i<nPads; i++ ) {

	//! select main detector
	DetectorInfo *det = dets[i]->GetMain();

	//! Set Cut
	char *buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut(buf);
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! Set histogram title
	string var_H = string(var)+">>"+det->TBName_;
	TH1 *h;
	if( (h = (TH1*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var_H.c_str(), detCut, "goff", nEnt );

	//! Retrieve and Draw histogram
	if( nPads ) cv_->cd(i+1);
	h = (TH1*) gDirectory->Get(det->TBName_.c_str());

	cout << "CheckTracks::DrawDet - var \"" 
	    << var 
	    << "\" det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h)?h->GetEntries():-1)
					<< " entries.\n";

	if( h && h->GetEntries() ) {
	    if( XTitle ) h->SetXTitle( XTitle );
	    if( YTitle ) h->SetYTitle( YTitle );
	    h->Draw(opt); 
	}

	_UpdateFromPad( gPad );


	//! loop over sub detectors
	for( unsigned int id=0; id< dets[i]->sub_.size(); id++ ) { 

	    DetectorInfo *det = dets[i]->sub_[id];

	    // skip main detector as already drawn
	    if( det == dets[i]->GetMain() ) continue;

	    //! for Pixel Detector: create corresponding dvVect names
	    string pvar;
	    string pcut;
	    if (det->TBName_.substr(0, 2)=="GP" ||
		det->TBName_.substr(0, 2)=="MP") {
	      pvar = "T_dvVect";
	      pcut = "abs(T_dvVect)<1";
	    }
	    
	    //! Set Cut
	    char *buf = new char[64];
	    sprintf( buf, "T_detVect==%i",det->id_ );
	    TCut detCut;
	    if (det->TBName_.substr(0, 2)=="GP" ||
		det->TBName_.substr(0, 2)=="MP") // Pixel det exception
	      detCut = TCut(pcut.c_str()) + TCut(buf);
	    else detCut = TCut(cut) + TCut(buf);
	    
	    if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );

	    //! Set histogram title
	    sprintf( buf, "%s_%i", det->TBName_.c_str(), id );
	    string name_H( buf );
	    string var_H;
	    if (det->TBName_.substr(0, 2)=="GP" ||
		det->TBName_.substr(0, 2)=="MP") // Pixel det exception
	      var_H = string(pvar)+">>"+name_H;
	    else var_H = string(var)+">>"+name_H;
	    SafeDelete(buf); 
	    
	    TH1 *h;
	    if( (h = (TH1*) gDirectory->Get(name_H.c_str())) ) SafeDelete(h);
	    tracks_->GetChain()->Draw( var_H.c_str(), detCut, "goff", nEnt );
	    
	    //! Retrieve and Draw histogram
	    if( nPads ) cv_->cd(i+1);
	    h = (TH1*) gDirectory->Get(name_H.c_str());
	    
	    cout << "CheckTracks::DrawDet - var \"" 
		 << var 
		 << "\" det \"" 
		 << det->TBName_.c_str() 
		 << "\": " 
		 << int((h)?h->GetEntries():-1)
		 << " entries.\n";
	    if( h && h->GetEntries() ) { 
	      if( XTitle ) h->SetXTitle( XTitle );
	      if( YTitle ) h->SetYTitle( YTitle );
	      h->SetLineColor( id+1 );
	      h->Draw("same");
	      _UpdateFromPad( gPad );
	    }
	} // loop over subdetectors      
	
    }
    
    _UpdateFromCanvas();
    if( string( opt ).find("same")==string::npos ) _MakeTitle("DrawDet \"%s\" - %s:%s", detselection, var, cut );
    
}

//________________________________________
void CheckTracks::DrawDU( const char* detselection, const char* cut, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::DrawDU - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawDU - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawDU - " << dets.size() << " detector(s) selected.\n";
    vector< TH3* > hV;

    //! Create all histograms
    for( unsigned int i=0; i<dets.size(); i++ ) {

	//! Generate var
	char *buf = new char[256];
	DetectorInfo *det = dets[i]->GetMain();
	TMatrix r = det->rotM_;
	if( tracks_->MagnetsOn() ) 
	    sprintf( buf, 
		    "T_duVect:((%f)*T_xVect+(%f)*T_yVect):((%f)*T_xVect+(%f)*T_yVect)>>%s",
		    r(0,0), 
		    r(1,0),
		    r(0,1),
		    r(0,0),
		    det->TBName_.c_str() );
	else 
	    sprintf( buf, 
		    "T_duVect:((%f)*(T_xLx+T_txLx*(%f-T_zLx))+(%f)*(T_yLx+T_tyLx*(%f-T_zLx)))"
		    ":((%f)*(T_xLx+T_txLx*(%f-T_zLx))+(%f)*(T_yLx+T_tyLx*(%f-T_zLx)))>>%s",
		    r(0,0), det->zcm_,
		    r(1,0), det->zcm_,
		    r(0,1), det->zcm_,
		    r(0,0), det->zcm_,
		    det->TBName_.c_str() ); 
	string var( buf );
	SafeDelete(buf);  

	//! Set Cut
	buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut( buf );
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! Set histogram title
	TH3 *h;
	if( (h = (TH3*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );
	h = (TH3*) gDirectory->Get(det->TBName_.c_str());
	if( h ) {
	    h->SetXTitle( "[mm]" );
	    h->SetYTitle( "[mm]" );
	    h->SetZTitle( "[mm]" );
	}
	hV.push_back( h );

	if( !h ) continue;

	cout << "CheckTracks::DrawDU - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h)?h->GetEntries():-1)
					<< " entries.\n";

    }

    //! Project/draw histograms
    //! DU Plots
    unsigned int nPads = dets.size();
    cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );
    if( ps_ ) _NewPage();

    for( unsigned int i=0; i<hV.size(); i++ ) {

	if( nPads ) cv_->cd(i+1);
	if( hV[i] && hV[i]->GetEntries() ) hV[i]->Project3D("z")->Draw(); 
	_UpdateFromPad( gPad );       
    }

    _UpdateFromCanvas( );
    if( ps_ ) _MakeTitle( "DrawDU - (%s) %s", detselection, cut );

    //! DU vs U Plots
    cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );
    if( ps_ ) _NewPage();

    for( unsigned int i=0; i<hV.size(); i++ ) {

	if( nPads ) cv_->cd(i+1);

	if( hV[i] && hV[i]->GetEntries() )  {
	    hV[i]->Project3D("zy")->Draw("box"); 
	    gPad->SetLogz();
	}

	_UpdateFromPad( gPad );       
    }

    _UpdateFromCanvas( );
    if( ps_ ) _MakeTitle( "DrawDUvsU - (%s) %s", detselection, cut );

    //! DU vs V Plots
    cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );
    if( ps_ ) _NewPage();

    for( unsigned int i=0; i<hV.size(); i++ ) {

	if( nPads ) cv_->cd(i+1);

	if( hV[i] && hV[i]->GetEntries() )  {
	    hV[i]->Project3D("zx")->Draw("box"); 
	    gPad->SetLogz();
	}

	_UpdateFromPad( gPad );       
    }

    _UpdateFromCanvas( );
    if( ps_ ) _MakeTitle( "DrawDUvsV - (%s) %s", detselection, cut );
    return;
}

//__________________________________________________________________
void CheckTracks::DrawDU_LR( const char* detselection, const char* cut, 
	int nEnt,
	double offsetN, 
	double offsetP )
{
    if( !df_ ) { cout << "CheckTracks::DrawDU_LR - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawDU_LR - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawDU_LR - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );
    if( ps_ ) _NewPage();

    //! Check TBranch T_rVect is present
    if( !tracks_->GetChain()->GetBranch( "T_rVect"   ) ) {
	cout << "CheckTracks::DrawDU_LR - Could not find branch T_rVect. aborted.\n";
	return;
    }

    for( unsigned int i=0; i<dets.size(); i++ ) {
	DetectorInfo *det = dets[i]->GetMain();
	if( nPads ) cv_->cd(i+1);
	if( det->type_ != 11 ) {         
	    cout << "CheckTracks::DrawDU_LR - det \"" 
		<< det->TBName_.c_str() 
		<< "\": not a drift-like detector.\n";
	    continue;
	} 

	//! Set detector cut
	char *buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut(buf);
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! Set neg R histogram title, var and cut
	string title_N = det->TBName_+"_N";
	buf = new char[64];
	if( offsetN ) sprintf( buf, "T_duVect+%f", offsetN );
	else sprintf( buf, "T_duVect" );
	string var_N = string(buf) + ">>" + title_N;
	TCut detCut_N = detCut + TCut("T_rVect<0");
	SafeDelete(buf);    

	TH1 *h_N;
	if( (h_N = (TH1*) gDirectory->Get(title_N.c_str())) ) SafeDelete(h_N);
	tracks_->GetChain()->Draw( var_N.c_str(), detCut_N, "goff", nEnt );
	h_N = (TH1*) gDirectory->Get(title_N.c_str());

	//! Set pos R histogram title, var and cut
	string title_P = det->TBName_+"_P";
	buf = new char[64];
	if( offsetP ) sprintf( buf, "T_duVect+%f", offsetP );
	else sprintf( buf, "T_duVect" );
	string var_P = string(buf) + ">>" + title_P;
	TCut detCut_P = detCut + TCut("T_rVect>0");
	SafeDelete(buf);    

	TH1 *h_P;
	if( (h_P = (TH1*) gDirectory->Get(title_P.c_str())) ) SafeDelete(h_P);
	tracks_->GetChain()->Draw( var_P.c_str(), detCut_P, "goff", nEnt );
	h_P = (TH1*) gDirectory->Get(title_P.c_str());

	if( (!(h_N && h_N->GetEntries())) && (!(h_P && h_P->GetEntries())) ) continue;

	//! Change histogram colors and axis title
	if( h_N ) {
	    h_N->SetXTitle( "[mm]" );
	    h_N->SetLineColor( 4 );
	}

	if( h_P ) {
	    h_P->SetXTitle( "[mm]" );
	    h_P->SetLineColor( 2 );
	}

	//! Get Histo with max value
	TH1* h1, *h2;
	if( h_N && ((!h_P) || h_P->GetMaximum()<h_N->GetMaximum() ) ) { h1 = h_N; h2 = h_P; }
	else{ h1 = h_P; h2 = h_N; }

	//! Dump number of entries
	cout << "CheckTracks::DrawDU_LR - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h_N)?h_N->GetEntries():-1)
					    << "/" 
						<< int((h_P)?h_P->GetEntries():-1)
										<< " entries.\n";

	//!Draw histograms
	if( h1 ) {
	    h1->Draw();
	    if( h2 ) h2->Draw("same");
	} else if( h2 ) h2->Draw();
	_UpdateFromPad( gPad );
    }

    _UpdateFromCanvas();
    _MakeTitle("DrawDU_LR \"%s\" - %s", detselection, cut );

    return;
}

//________________________________________
void CheckTracks::DrawDUvsU( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::DrawDUvsU - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawDUvsU - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawDUvsU - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	if( nPads ) cv_->cd(i+1);

	//! Generate var
	char *buf = new char[128];
	DetectorInfo *det = dets[i]->GetMain();
	TMatrix r = det->rotM_;
	if( tracks_->MagnetsOn() ) 
	    sprintf( buf, 
		    "T_duVect:((%f)*T_xVect+(%f)*T_yVect)>>%s",
		    r(0,0), 
		    r(1,0),
		    det->TBName_.c_str() );
	else 
	    sprintf( buf, 
		    "T_duVect:((%f)*(T_xLx+T_txLx*(%f-T_zLx))+(%f)*(T_yLx+T_tyLx*(%f-T_zLx)))>>%s",
		    r(0,0), det->zcm_,
		    r(1,0), det->zcm_,
		    det->TBName_.c_str() ); 
	string var( buf );
	SafeDelete(buf);  

	//! Set Cut
	buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut( buf );
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! erase old histogram with same name, if needed, create new one
	TH2 *h;
	if( (h = (TH2*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	//! Retrieve Histogram
	h = (TH2*) gDirectory->Get(det->TBName_.c_str());
	if(!h) {
	    cout << "CheckTracks::DrawDUvsU - det \"" 
		<< det->TBName_.c_str() 
		<< "\": -1 entries"<<endl;
	    return;
	}
	double mean = h->GetMean();
	double rms = h->GetRMS();
	SafeDelete(h);

	//! Create profile. Then fill it
	string title = var.substr(0,var.find(">>"))+string("{")+string(detCut)+string("}");
	TProfile* hh = new TProfile(det->TBName_.c_str(),title.c_str(),50,mean-3*rms,mean+3*rms);
	Double_t* w = new Double_t[tracks_->GetChain()->GetSelectedRows()];
	for(Long64_t i = 0; i<tracks_->GetChain()->GetSelectedRows(); i++) w[i]=1;
	hh->FillN(tracks_->GetChain()->GetSelectedRows(),tracks_->GetChain()->GetV2(),tracks_->GetChain()->GetV1(),w,1);

	cout << "CheckTracks::DrawDUvsU - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((hh)?hh->GetEntries():-1)
					<< " entries.\n";


	if( hh )
	    if( hh->GetEntries() )  {
		hh->Draw(opt);
		hh->SetXTitle("[mm]");
		hh->SetYTitle("[mm]");
		gPad->SetLogz();
	    }

	_UpdateFromPad( gPad );       


    }

    _UpdateFromCanvas( );
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "DrawDUvsU - (%s) %s", detselection, cut );

}

//________________________________________
void CheckTracks::PrintDUvsU( const char* detselection, float du, float theta, float umin,float umax,const char* cut, const char* opt, int nEnt ){

    if( !df_ ) { cout << "CheckTracks:PrintDUvsU - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::PrintDUvsU - no tracks loaded.\n"; return; }
    //  if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    vector< DetectorInfo* > alldets= df_->GetDetInfo();
    cout << "CheckTracks::PrintDUvsU - " << dets.size() << " detector(s) selected.\n";

    //  unsigned int nPads = dets.size();
    //if( string( opt ).find("same")==string::npos ) cv_->Clear();
    //if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	tracks_->Init();

	//! Generate var

	DetectorInfo* det = dets[i]->GetMain();


	TMatrix r = det->rotM_;



	//char* buf = new char[128];
	//sprintf( buf, "T_detVect==%i",det->id_ );
	//TCut detCut = TCut(cut) + TCut( buf );
	//if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	//SafeDelete(buf);       

	// tracks_->AddCut(buf);


	long count=0;
	double x=0.,x2=0.,y=0.,xy=0.,delta=0.; // line sub parameters;
	double A,B,dA,dB;
	double sigmay2=0;

	while(tracks_->GetNextEntry()){



	    bool  iffound=0;
	    int sid=0; //selected id

	    for(int j=0; j<tracks_->T_nDets; j++){
		if(tracks_->T_detVect[j]==det->id_&&fabs(tracks_->T_duVect[j])<du&&
			((!tracks_->MagnetsOn()&&sqrt(pow(tracks_->T_txLx,2.0)+pow(tracks_->T_tyLx,2.0))<theta)||
			 (tracks_->MagnetsOn()&&sqrt(pow(tracks_->T_txVect[j],2.0)+pow(tracks_->T_tyVect[j],2.0))<theta))&&
			tracks_->T_uVect[j]>umin&&tracks_->T_uVect[j]<umax){ 

		    iffound=1;
		    count+=1;
		    sid=j;
		    j=tracks_->T_nDets+1;
		}   
	    }






	    if(iffound){

		x+= tracks_->T_uVect[sid];
		x2+=tracks_->T_uVect[sid]*tracks_->T_uVect[sid];
		y+= tracks_->T_duVect[sid];
		xy+=tracks_->T_uVect[sid]*tracks_->T_duVect[sid];
	    }
	}  //while loop


	delta=count*x2-x*x;
	A=(x2*y-x*xy)/delta;
	B=(count*xy-x*y)/delta;
	//    dB=0.1*sqrt(count/delta);
	//    cout<<A<<"  "<<B<<"<<"   "<<count<<endl;



	tracks_->Init();   //second loop needed to compute errors.

	while(tracks_->GetNextEntry()){


	    for(int j=0; j<tracks_->T_nDets; j++){
		if(tracks_->T_detVect[j]==det->id_&&fabs(tracks_->T_duVect[j])<du&&        
			((!tracks_->MagnetsOn()&&sqrt(pow(tracks_->T_txLx,2.0)+pow(tracks_->T_tyLx,2.0))<theta)||
			 (tracks_->MagnetsOn()&&sqrt(pow(tracks_->T_txVect[j],2.0)+pow(tracks_->T_tyVect[j],2.0))<theta))&&
			tracks_->T_uVect[j]>umin&&tracks_->T_uVect[j]<umax){ 

		    sigmay2+=pow(tracks_->T_duVect[j]-A-tracks_->T_uVect[j]*B,2.0)/(count-2);

		    j=tracks_->T_nDets+1;
		}   
	    }
	}   //while loop2


	dB=sqrt(sigmay2)*sqrt(count/delta);
	dA=sqrt(sigmay2)*sqrt(x2/delta);
	cout<<"------------------------------------"<<endl;
	cout<<det->TBName_.c_str()<<"number of tracks:"<<count<<endl; 

	cout<<"const:"<<A<<"+/-"<<dA<<"<<[mm]"<<endl<<endl;
	cout<<"SLOPE in 10-4:"<<10000*B<<"+/-"<<10000*dB<<endl;
	//        cout<<"------------------------------------"<<endl;
    }  //det loop
}

//________________________________________
void CheckTracks::DrawDUvsV( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::DrawDUvsV - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawDUvsV - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawDUvsV - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	if( nPads ) cv_->cd(i+1);

	//! Generate var
	char *buf = new char[128];
	DetectorInfo* det = dets[i]->GetMain();
	TMatrix r = det->rotM_;
	if( tracks_->MagnetsOn() )
	    sprintf( buf, 
		    "T_duVect:((%f)*T_xVect+(%f)*T_yVect)>>%s",
		    r(0,1), 
		    r(0,0),
		    det->TBName_.c_str() );
	else 
	    sprintf( buf, 
		    "T_duVect:((%f)*(T_xLx+T_txLx*(%f-T_zLx))+(%f)*(T_yLx+T_tyLx*(%f-T_zLx)))>>%s",
		    r(0,1), det->zcm_,
		    r(0,0), det->zcm_,
		    det->TBName_.c_str() ); 
	string var(buf);
	SafeDelete(buf);   

	//! Set Cut
	buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut(buf);
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! erase old histogram with same name, if needed, create new one
	TH2 *h;
	if( (h = (TH2*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	//! Retrieve Histogram
	h = (TH2*) gDirectory->Get(det->TBName_.c_str());
	if(!h) {
	    cout << "CheckTracks::DrawDUvsV - det \"" 
		<< det->TBName_.c_str() 
		<< "\": -1 entries"<<endl;
	    return;
	}
	double mean = h->GetMean();
	double rms = h->GetRMS();
	SafeDelete(h);

	//! Create profile. Then fill it
	string title = var.substr(0,var.find(">>"))+string("{")+string(detCut)+string("}");
	TProfile* hh = new TProfile(det->TBName_.c_str(),title.c_str(),50,mean-3*rms,mean+3*rms);
	Double_t* w = new Double_t[tracks_->GetChain()->GetSelectedRows()];
	for(Long64_t i = 0; i<tracks_->GetChain()->GetSelectedRows(); i++) w[i]=1;
	hh->FillN(tracks_->GetChain()->GetSelectedRows(),tracks_->GetChain()->GetV2(),tracks_->GetChain()->GetV1(),w,1);

	cout << "CheckTracks::DrawDUvsV - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((hh)?hh->GetEntries():-1)
					<< " entries.\n";

	if( hh )
	    if( hh->GetEntries() )  {
		hh->Draw(opt);
		hh->SetXTitle("[mm]");
		hh->SetYTitle("[mm]");
		gPad->SetLogz();
	    }

	_UpdateFromPad( gPad );    
    }

    _UpdateFromCanvas();
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "DrawDUvsV - (%s) %s", detselection, cut );

}

//_________________________________________________________________________________
void CheckTracks::DrawDUvsP( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::DrawDUvsP - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawDUvsP - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawDUvsP - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	//! Generate var
	char *buf = new char[128];
	DetectorInfo* det = dets[i]->GetMain();
	TMatrix r = det->rotM_;
	if( tracks_->MagnetsOn() )
	    sprintf( buf, 
		    "T_duVect:1/T_cop>>%s",
		    det->TBName_.c_str() );
	else {
	  SafeDelete(buf);
	  return;
	}

	string var(buf);
	SafeDelete(buf);   

	//! Set Cut
	buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut(buf);
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! erase old histogram with same name, if needed, create new one
	TH2 *h;
	if( (h = (TH2*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	cout<<detCut<<endl; //temp
	//! Retrieve and Draw Histogram
	if( nPads ) cv_->cd(i+1);
	h = (TH2*) gDirectory->Get(det->TBName_.c_str());

	cout << "CheckTracks::DrawDUvsP - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h)?h->GetEntries():-1)
					<< " entries.\n";

	if( h && h->GetEntries() )  {
	    h->SetXTitle("[GeV/c]");
	    h->SetYTitle("[mm]");
	    h->Draw(opt); 
	    gPad->SetLogz();

	}
	_UpdateFromPad( gPad );    
    }

    _UpdateFromCanvas();
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "DrawDUvsP - (%s) %s", detselection, cut );

}

//____________________________________________________________________________
void CheckTracks::PrintDUvsP( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks:PrintDUvsP - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::PrintDUvsP - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    vector< DetectorInfo* > alldets= df_->GetDetInfo();
    cout << "CheckTracks::PrintDUvsP - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	tracks_->Init();

	//! Generate var
	char *buf = new char[128];
	DetectorInfo* det = dets[i]->GetMain();


	TMatrix r = det->rotM_;
	if( tracks_->MagnetsOn() ) 
	    sprintf( buf, 
		    "1/T_cop>100");
	else {
	  SafeDelete(buf);   
	  return;
	}

	string var(buf);
	SafeDelete(buf);   


	long count=0;
	float momarray[]={-150.,-35.,-20.,-12.,-8.,-5., -3., 0., 3., 5., 8., 12., 20., 35., 50., 65.,
	    80., 95., 110., 125., 135., 145., 155., 165., 175., 185.,250.};

	// float  meanmom[sizeof(momarray)/sizeof(float)-1]= {}; // what is it for? (jj)
	double  meanDU[sizeof(momarray)/sizeof(float)-1]= {};
	double sigmaDU[sizeof(momarray)/sizeof(float)-1]= {};
	int      numDU[sizeof(momarray)/sizeof(float)-1]= {}; 

	while(tracks_->GetNextEntry()){


	    //       vector< DetectorInfo* > detintrack;
	    //       for(int i=0;i<tracks_->T_nDets;i++)
	    //           detintrack.push_back(df_->GetDetInfo(tracks_->T_detVect[i]));

	    //      if(!tracks_->AcceptEntry(alldets,NPARTRCK)) continue;
	    bool  iffound=0;
	    int sid=0; 

	    for(int j=0; j<tracks_->T_nDets; j++){
		if(tracks_->T_detVect[j]==det->id_){ 
		    iffound=1;
		    count+=1;
		    sid=j;
		    j=tracks_->T_nDets+1;
		}   
	    }



	    if(iffound){ 

		for(unsigned int j=0;j<sizeof(momarray)/sizeof(float)-1;j++){
		    if(1.0/tracks_->T_cop>momarray[j]&&1.0/tracks_->T_cop<momarray[j+1]&&
			    fabs(tracks_->T_duVect[sid])<10*det->res_){

			numDU[j]+=1;
			meanDU[j]+=tracks_->T_duVect[sid];
			sigmaDU[j]+=tracks_->T_duVect[sid]*tracks_->T_duVect[sid];

		    }}}
	}  //while loop

	//       cout<<count<<endl;

	cout<<"------------------------------------"<<endl;
	cout<<det->TBName_.c_str()<<"number of tracks:"<<count<<endl; 
	for(unsigned int j=0;j<sizeof(momarray)/sizeof(float)-1;j++){

	    meanDU[j]=meanDU[j]/numDU[j];
	    sigmaDU[j]=sqrt(sigmaDU[j]/numDU[j]-meanDU[j]*meanDU[j])/sqrt((double)numDU[j]);
	    cout<<momarray[j]<<':'<<momarray[j+1]<<"   "<<1000*meanDU[j]<<"+/-"<<1000*sigmaDU[j]<<"  "<<numDU[j]<<endl;
	    //       cout<<1000*meanDU[j]<<", "<<endl;
	}



    }  //det loop



}

//_______________________________________________________________________________________________
void CheckTracks::DrawDUvsT( const char* detselection, const char* cut,  const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::DrawDUvsT - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawDUvsT - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawDUvsT - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	if( nPads ) cv_->cd(i+1);

	//! Generate var
	char *buf = new char[128];
	DetectorInfo *det = dets[i]->GetMain();
	TMatrix r = det->rotM_;
	if( tracks_->MagnetsOn() ) 
	    sprintf( buf, 
		    "T_duVect:((%f)*T_txVect+(%f)*T_tyVect)>>%s",
		    r(0,0), 
		    r(1,0),
		    det->TBName_.c_str() );
	else 
	    sprintf( buf, 
		    "T_duVect:((%f)*T_txLx+(%f)*T_tyLx)>>%s",
		    r(0,0), 
		    r(1,0),
		    det->TBName_.c_str() );

	string var(buf);
	SafeDelete(buf);  

	//! Set Cut
	buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut( buf );
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! erase old histogram with same name, if needed, create new one
	TH2 *h;
	if( (h = (TH2*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	//! Retrieve the 1D Histogram
	h = (TH2*) gDirectory->Get(det->TBName_.c_str());
	if(!h) {
	    cout << "CheckTracks::DrawDUvsT - det \"" 
		<< det->TBName_.c_str() 
		<< "\": -1 entries"<<endl;
	    return;
	}
	double mean = h->GetMean();
	double rms = h->GetRMS();
	SafeDelete(h);

	//! Create profile. Then fill it
	string title = var.substr(0,var.find(">>"))+string("{")+string(detCut)+string("}");
	TProfile* hh = new TProfile(det->TBName_.c_str(),title.c_str(),50,mean-3*rms,mean+3*rms);
	Double_t* w = new Double_t[tracks_->GetChain()->GetSelectedRows()];
	for(Long64_t i = 0; i<tracks_->GetChain()->GetSelectedRows(); i++) w[i]=1;
	hh->FillN(tracks_->GetChain()->GetSelectedRows(),tracks_->GetChain()->GetV2(),tracks_->GetChain()->GetV1(),w,1);

	cout << "CheckTracks::DrawDUvsT - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((hh)?hh->GetEntries():-1)
					  << " entries."<<endl;

	//! Draw the plot
	if( hh )
	    if( hh->GetEntries() ){
		hh->Draw(opt);
		hh->SetXTitle("[tg theta]");
		hh->SetYTitle("[mm]");
		gPad->SetLogz();
	    }

	_UpdateFromPad( gPad );      


    }

    _UpdateFromCanvas( );
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "DrawDUvsT - (%s) %s", detselection, cut );

}

//________________________________________________________________________
void CheckTracks::DrawDUvsT90( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::DrawDUvsT90 - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawDUvsT90 - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawDUvsT90 - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	if( nPads ) cv_->cd(i+1);

	//! Generate var
	char *buf = new char[128];
	DetectorInfo *det = dets[i]->GetMain();
	TMatrix r = det->rotM_;
	if( tracks_->MagnetsOn() ) 
	    sprintf( buf, 
		    "T_duVect:((%f)*T_txVect+(%f)*T_tyVect)>>%s",
		    r(0,1), 
		    r(0,0),
		    det->TBName_.c_str() );
	else 
	    sprintf( buf, 
		    "T_duVect:((%f)*T_txLx+(%f)*T_tyLx)>>%s",
		    r(0,1), 
		    r(0,0),
		    det->TBName_.c_str() );

	string var(buf);
	SafeDelete(buf);  

	//! Set Cut
	buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut( buf );
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! erase old histogram with same name, if needed, create new one
	TH2 *h;
	if( (h = (TH2*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	//! Retrieve and Draw Histogram
	h = (TH2*) gDirectory->Get(det->TBName_.c_str());
	if(!h) {
	    cout << "CheckTracks::DrawDUvsT - det \"" 
		<< det->TBName_.c_str() 
		<< "\": -1 entries"<<endl;
	    return;
	}
	double mean = h->GetMean();
	double rms = h->GetRMS();
	SafeDelete(h);

	//! Create profile. Then fill it
	string title = var.substr(0,var.find(">>"))+string("{")+string(detCut)+string("}");
	TProfile* hh = new TProfile(det->TBName_.c_str(),title.c_str(),100,mean-3*rms,mean+3*rms);
	Double_t* w = new Double_t[tracks_->GetChain()->GetSelectedRows()];
	for(Long64_t i = 0; i<tracks_->GetChain()->GetSelectedRows(); i++) w[i]=1;
	hh->FillN(tracks_->GetChain()->GetSelectedRows(),tracks_->GetChain()->GetV2(),tracks_->GetChain()->GetV1(),w,1);

	cout << "CheckTracks::DrawDUvsT90 - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((hh)?hh->GetEntries():-1)
					<< " entries.\n";

	if( hh && hh->GetEntries() )  {
	    hh->Draw(opt);
	    hh->SetXTitle("[tg theta]");
	    hh->SetYTitle("[mm]");

	    gPad->SetLogz();
	}

	_UpdateFromPad( gPad );       

    }

    _UpdateFromCanvas( );
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "DrawDUvsT90 - (%s) %s", detselection, cut );

}

//________________________________________________________________________________
void CheckTracks::DrawDUvsZ( const char* detselection, const char* cut, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::DrawDUvsZ - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawDUvsZ - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();

    cv_->Clear();
    _NewPage();

    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawDUvsZ - " << dets.size() << " detector(s) selected.\n";

    //! Create the TGraph for all points
    TGraphErrors *tG_all = _AddTGE();
    string title = string(detselection)+"::"+string(cut);
    tG_all->SetTitle( title.c_str() );
    unsigned int iTG_all=0;

    //! Create TGraph for detector types
    vector< TGraphErrors* > tG;
    vector< int> iTG;

    const int iFI = 0; tG.push_back( _AddTGE( 1 ) ); iTG.push_back( 0 );
    const int iGM = 1; tG.push_back( _AddTGE( 2 ) ); iTG.push_back( 0 );
    const int iDC = 2; tG.push_back( _AddTGE( 3 ) ); iTG.push_back( 0 );
    const int iMM = 3; tG.push_back( _AddTGE( 4 ) ); iTG.push_back( 0 );
    const int iSI = 4; tG.push_back( _AddTGE( 5 ) ); iTG.push_back( 0 );

    for( unsigned int i=0; i<dets.size(); i++ ) {
	DetectorInfo* det = dets[i]->GetMain();
	string name =  det->TBName_+"_"+"du";
	string title = name + ":" + string(cut);  
	string select = string("T_duVect>>")+name;    

	char* buf = new char[256];
	sprintf( buf, "T_detVect==%i", det->id_ );
	TCut detCut = TCut( buf )+cut;
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);


	TH1S *h;
	if( (h = (TH1S*) gDirectory->Get(name.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw(select.c_str(), detCut, "goff", nEnt );
	h = (TH1S*) gDirectory->Get(name.c_str());
	if( !h ) continue;

	cout << "CheckTracks::DrawDUvsZ - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h)?h->GetEntries():-1)
					<< " entries.\n";

	tG_all->SetPoint( iTG_all, det->zcm_, h->GetMean() );
	tG_all->SetPointError(iTG_all, 0, h->GetRMS()/sqrt(h->GetEntries()) );
	iTG_all++;

	//! Fill subTGraphs
	string TB = det->TBName_.substr(0, 2);
	if( TB == "FI" ) {
	    tG[iFI]->SetPoint( iTG[iFI], det->zcm_, h->GetMean() );
	    tG[iFI]->SetPointError(iTG[iFI], 0, h->GetRMS()/sqrt(h->GetEntries()) );
	    iTG[iFI]++;
	} else if( TB == "GM" ) {
	    tG[iGM]->SetPoint( iTG[iGM], det->zcm_, h->GetMean() );
	    tG[iGM]->SetPointError(iTG[iGM], 0, h->GetRMS()/sqrt(h->GetEntries()) );
	    iTG[iGM]++;
	} else if( TB == "MM" ) {
	    tG[iMM]->SetPoint( iTG[iMM], det->zcm_, h->GetMean() );
	    tG[iMM]->SetPointError(iTG[iMM], 0, h->GetRMS()/sqrt(h->GetEntries()) );
	    iTG[iMM]++;
	} else if( TB == "DC" ) {
	    tG[iDC]->SetPoint( iTG[iDC], det->zcm_, h->GetMean() );
	    tG[iDC]->SetPointError(iTG[iDC], 0, h->GetRMS()/sqrt(h->GetEntries()) );
	    iTG[iDC]++;
	} else if( TB == "SI" ) {
	    tG[iSI]->SetPoint( iTG[iSI], det->zcm_, h->GetMean() );
	    tG[iSI]->SetPointError(iTG[iSI], 0, h->GetRMS()/sqrt(h->GetEntries()) );
	    iTG[iSI]++;
	}  

    }
    tG_all->GetXaxis()->SetTitle( "[mm]" );
    tG_all->GetYaxis()->SetTitle( "[mm]" );
    tG_all->Draw( "AP" );
    for( unsigned int i=0; i<tG.size(); i++ ) 
	if( iTG[i] ) {
	    tG[i]->GetXaxis()->SetTitle( "[mm]" );
	    tG[i]->GetYaxis()->SetTitle( "[mm]" );
	    tG[i]->Draw("P");
	}
    _UpdateFromCanvas();
    _MakeTitle( "DrawDUvsZ - (%s) %s", detselection, cut );
}    

//________________________________________
void CheckTracks::DrawTrackMlt( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::DrawTrackMlt - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawTrackMlt - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetChain()->GetEntries();

    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawTrackMlt - " << dets.size() << " detector(s) selected.\n";
    int nDets = dets.size();
    TH1S *h = new TH1S("TrakMlt", "TrakMlt", nDets+1, -0.5, nDets+0.5 );

    //! dump detectors
    for( int id=0; id<nDets; id++ )
	cout << "CheckTracks::DrawTrackMl - \"" << dets[id]->GetMain()->TBName_ << "\" ("
	    << dets[id]->GetMain()->id_ << ") added .\n";

    for( int i=0; i<nEnt; i++ ) {

	//! Try to retrieve entry i
	if( !tracks_->GetEntry(i) ) return;

	//! Check if entry match the cut
	if( !tracks_->GetChain()->Draw("T_evt",cut, "goff", 1, i ) ) continue;

	//! Loop over dets in the track
	int n=0;
	for( int jd=0; jd<nDets; jd++ ) 
	    for( int id=0; id<tracks_->T_nDets; id++ ) 
		if( tracks_->T_detVect[id] == dets[jd]->GetMain()->id_ ) {
		    n++;
		    break; //!< Detectors are only counted once/track
		}
	h->Fill( n );
    }  

    Utils::DivideCanvas( cv_, 1 );
    cv_->cd(1);
    h->Draw( opt );

    _UpdateFromPad(gPad);
    _UpdateFromCanvas();
    if( string( opt ).find("same")==string::npos ) _MakeTitle("DrawTrackMlt \"%s\" - %s", detselection, cut );

}

//________________________________________
void CheckTracks::DrawRT( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::DrawRT - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::DrawRT - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::DrawRT - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {
	DetectorInfo* det = dets[i]->GetMain();

	//! Check Detector name
	string TB = det->TBName_.substr(0,2);
	if( !( TB=="DC" || TB=="DW" || TB=="MB" || TB=="ST" ) ) {
	    cout << "CheckTracks::DrawRT - det \"" 
		<< det->TBName_.c_str() 
		<< "\": not a drift-like detector.\n";
	    continue;
	}

	//! Generate detector Cut
	char* buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut(buf);
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! Generate var [u_track-u_wire = (u_cluster-u_wire)-(u_cluster-u_track) = T_rVect-T_duVect]
	string var = "T_rVect-T_duVect:T_tVect>>"+det->TBName_;
	TH1 *h;
	if( (h = (TH1*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	//! Retrieve and Draw Histogram
	if( nPads ) cv_->cd(i+1);
	h = (TH1*) gDirectory->Get(det->TBName_.c_str());

	cout << "CheckTracks::DrawRT - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h)?h->GetEntries():-1)
					<< " entries.\n";

	if( h && h->GetEntries() ) {
	    h->SetXTitle( "[ns]" );
	    h->SetYTitle( "[mm]" );
	    h->Draw(opt); 
	}

	_UpdateFromPad( gPad );    
    }

    _UpdateFromCanvas();
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "DrawDUvsV - (%s) %s", detselection, cut );

}


//________________________________________
void CheckTracks::FitDUvsU( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::FitDUvsU - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::FitDUvsU - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::FitDUvsU - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {
	if( nPads ) cv_->cd(i+1);

	//! Generate var
	char *buf = new char[128];
	DetectorInfo *det = dets[i]->GetMain();
	TMatrix r = det->rotM_;
	if( tracks_->MagnetsOn() ) 
	    sprintf( buf, 
		    "T_duVect:((%f)*T_xVect+(%f)*T_yVect)>>%s",
		    r(0,0), 
		    r(1,0),
		    det->TBName_.c_str() );
	else 
	    sprintf( buf, 
		    "T_duVect:((%f)*(T_xLx+T_txLx*(%f-T_zLx))+(%f)*(T_yLx+T_tyLx*(%f-T_zLx)))>>%s",
		    r(0,0), det->zcm_,
		    r(1,0), det->zcm_,
		    det->TBName_.c_str() ); 
	string var( buf );
	SafeDelete(buf);  

	//! Set Cut
	buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut( buf );
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! erase old histogram with same name, if needed, create new one
	TH2 *h;
	if( (h = (TH2*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	//! Retrieve and Draw Histogram
	h = (TH2*) gDirectory->Get(det->TBName_.c_str());

	cout << "CheckTracks::FitDUvsU - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h)?h->GetEntries():-1)
					<< " entries - ";

	if( h && h->GetEntries() )  {
	    h->SetXTitle( "[mm]" );
	    h->SetYTitle( "[mm]" );
	    h->Draw(opt); 
	    gPad->SetLogz();

	    // do the fit
	    double uMin =  h->GetXaxis()->GetBinCenter( 1 );
	    double uMax =  h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
	    char f_name[256];
	    sprintf(f_name,"P1_%i",i);
	    TF1* f = new TF1(f_name, Fit::P1, uMin, uMax, 2);
	    f->SetParameter( 0, 0 );
	    f->SetParameter( 1, 0 );
	    f->SetLineColor( 2 );
	    f->SetLineWidth( 2 );
	    f->SetFillStyle( 0 );
	    TH2Fit* fit = new TH2Fit( f, 2 );
	    fit->Fit( h, uMin, uMax );
	    printf( "dU=%.3g+%.3g*U\n", f->GetParameter(0), f->GetParameter(1) );
	    _PutText( gPad, 0.15, 0.15, 0.06, "dU=%.3g+%.3g*U", f->GetParameter(0), f->GetParameter(1) );
	    f->Draw("same");
	} else cout << endl;

	_UpdateFromPad( gPad );       


    }

    _UpdateFromCanvas( );
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "FitDUvsU - (%s) %s", detselection, cut );

}

//________________________________________
void CheckTracks::FitDUvsV( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ )     { cout << "CheckTracks::FitDUvsV - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::FitDUvsV - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::FitDUvsV - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	//! Generate var
	char *buf = new char[128];
	DetectorInfo* det = dets[i]->GetMain();
	TMatrix r = det->rotM_;
	if( tracks_->MagnetsOn() )
	    sprintf( buf, 
		    "T_duVect:((%f)*T_xVect+(%f)*T_yVect)>>%s",
		    r(0,1), 
		    r(0,0),
		    det->TBName_.c_str() );
	else 
	    sprintf( buf, 
		    "T_duVect:((%f)*(T_xLx+T_txLx*(%f-T_zLx))+(%f)*(T_yLx+T_tyLx*(%f-T_zLx)))>>%s",
		    r(0,1), det->zcm_,
		    r(0,0), det->zcm_,
		    det->TBName_.c_str() ); 
	string var(buf);
	SafeDelete(buf);   

	//! Set Cut
	buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut(buf);
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! erase old histogram with same name, if needed, create new one
	TH2 *h;
	if( (h = (TH2*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	//! Retrieve and Draw Histogram
	if( nPads ) cv_->cd(i+1);
	h = (TH2*) gDirectory->Get(det->TBName_.c_str());

	cout << "CheckTracks::FitDUvsV - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h)?h->GetEntries():-1)
					<< " entries - ";

	if( h && h->GetEntries() )  {
	    h->SetXTitle( "[mm]" );
	    h->SetYTitle( "[mm]" );
	    h->Draw(opt); 
	    gPad->SetLogz();

	    // do the fit
	    double uMin =  h->GetXaxis()->GetBinCenter( 1 );
	    double uMax =  h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
	    char f_name[256];
	    sprintf(f_name,"P1_%i",i);
	    TF1* f = new TF1(f_name, Fit::P1, uMin, uMax, 2);
	    f->SetParameter( 0, 0 );
	    f->SetParameter( 1, 0 );
	    f->SetLineColor( 2 );
	    f->SetLineWidth( 2 );
	    f->SetFillStyle( 0 );
	    TH2Fit* fit = new TH2Fit( f, 2 );
	    fit->Fit( h, uMin, uMax );
	    printf( "dU=%.3g+%.3g*V\n", f->GetParameter(0), f->GetParameter(1) );
	    _PutText( gPad, 0.15, 0.15, 0.06, "dU=%.3g+%.3g*V", f->GetParameter(0), f->GetParameter(1) );
	    f->Draw("same");

	} else cout << endl;
	_UpdateFromPad( gPad );    
    }

    _UpdateFromCanvas();
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "FitDUvsV - (%s) %s", detselection, cut );

}

//________________________________________
void CheckTracks::FitDUvsT( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::FitDUvsT - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::FitDUvsT - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::FitDUvsT - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	//! Generate var
	char *buf = new char[128];
	DetectorInfo *det = dets[i]->GetMain();
	TMatrix r = det->rotM_;

	if( tracks_->MagnetsOn() ) 
	    snprintf( buf, 127,
		    "T_duVect:((%f)*T_txVect+(%f)*T_tyVect)>>%s",
		    r(0,0), 
		    r(1,0),
		    det->TBName_.c_str() );
	else 
	    snprintf( buf, 127,
		    "T_duVect:((%f)*T_txLx+(%f)*T_tyLx)>>%s",
		    r(0,0), 
		    r(1,0),
		    det->TBName_.c_str() );

	string var( buf );
	SafeDelete(buf);  

	//! Set Cut
	buf = new char[64];
	snprintf( buf, 64, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut( buf );
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);

	//! erase old histogram with same name, if needed, create new one
	TH2 *h;
	if( (h = (TH2*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	//! Retrieve and Draw Histogram
	if( nPads ) cv_->cd(i+1);
	h = (TH2*) gDirectory->Get(det->TBName_.c_str());

	cout << "CheckTracks::FitDUvsT - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h)?h->GetEntries():-1)
					<< " entries - ";

	if( h && h->GetEntries() )  {
	    h->SetXTitle( "tg theta" );
	    h->SetYTitle( "[mm]" );
	    h->Draw(opt); 
	    gPad->SetLogz();

	    // do the fit
	    double uMin =  h->GetXaxis()->GetBinCenter( 1 );
	    double uMax =  h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
	    buf = new char[64];
	    snprintf(buf, 63, "P1_%i",i);
	    TF1* f = new TF1(buf, Fit::P1, uMin, uMax, 2);
	    SafeDelete(buf);
	    f->SetParameter( 0, 0 );
	    f->SetParameter( 1, 0 );
	    f->SetLineColor( 2 );
	    f->SetLineWidth( 2 );
	    f->SetFillStyle( 0 );
	    TH2Fit* fit = new TH2Fit( f, 2 );
	    fit->Fit( h, uMin, uMax );
	    //      cout.form( "dU=%.3g+%.3g*T\n", f->GetParameter(0), f->GetParameter(1) );
	    //      _PutText( gPad, 0.15, 0.15, 0.06, "dU=%.3g+%.3g*T", f->GetParameter(0), f->GetParameter(1) );
	    buf = new char[128];
	    snprintf(buf, 127, "dU=%.3g+%.3g*T", f->GetParameter(0), f->GetParameter(1));
	    cout<<buf<<endl;
	    SafeDelete(buf);
	    _PutText( gPad, 0.15, 0.15, 0.06, "dZ=%.3g",f->GetParameter(1) );

	    //       cout.form( "dZ=%.3g+-%.3g\n", f->GetParameter(1), f->GetParError(1) );
	    //      _PutText( gPad, 0.15, 0.15, 0.06, "dZ=%.3g+-%.3g", f->GetParameter(1), f->GetParError(1) );
	    f->Draw("same");

	} else cout << endl;

	_UpdateFromPad( gPad );       


    }

    _UpdateFromCanvas( );
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "FitDUvsT - (%s) %s", detselection, cut );

}

//________________________________________
void CheckTracks::FitDUvsT90( const char* detselection, const char* cut, const char* opt, int nEnt )
{
    if( !df_ ) { cout << "CheckTracks::FitDUvsT90 - no detectorfile loaded.\n"; return; }
    if( !tracks_ ) { cout << "CheckTracks::FitDUvsT90 - no tracks loaded.\n"; return; }
    if( !cv_ ) MakeCanvas();
    if( !nEnt ) nEnt = (int) tracks_->GetEntries();
    vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
    cout << "CheckTracks::FitDUvsT90 - " << dets.size() << " detector(s) selected.\n";

    unsigned int nPads = dets.size();
    if( string( opt ).find("same")==string::npos ) cv_->Clear();
    if( nPads ) Utils::DivideCanvas( cv_, nPads );

    if( ps_ && string( opt ).find("same")==string::npos ) _NewPage();

    for( unsigned int i=0; i<dets.size(); i++ ) {

	//! Generate var
	char *buf = new char[128];
	DetectorInfo *det = dets[i]->GetMain();
	TMatrix r = det->rotM_;

	if( tracks_->MagnetsOn() ) 
	    sprintf( buf, 
		    "T_duVect:((%f)*T_txVect+(%f)*T_tyVect)>>%s",
		    r(0,1), 
		    r(0,0),
		    det->TBName_.c_str() );
	else 
	    sprintf( buf, 
		    "T_duVect:((%f)*T_txLx+(%f)*T_tyLx)>>%s",
		    r(0,1), 
		    r(0,0),
		    det->TBName_.c_str() );

	string var( buf );
	SafeDelete(buf);  

	//! Set Cut
	buf = new char[64];
	sprintf( buf, "T_detVect==%i",det->id_ );
	TCut detCut = TCut(cut) + TCut( buf );
	if( glCut_.size() ) detCut+=TCut( glCut_.c_str() );
	SafeDelete(buf);    

	//! erase old histogram with same name, if needed, create new one
	TH2 *h;
	if( (h = (TH2*) gDirectory->Get(det->TBName_.c_str())) ) SafeDelete(h);
	tracks_->GetChain()->Draw( var.c_str(), detCut, "goff", nEnt );

	//! Retrieve and Draw Histogram
	if( nPads ) cv_->cd(i+1);
	h = (TH2*) gDirectory->Get(det->TBName_.c_str());

	cout << "CheckTracks::FitDUvsT90 - det \"" 
	    << det->TBName_.c_str() 
	    << "\": " 
	    << int((h)?h->GetEntries():-1)
					<< " entries - ";

	if( h && h->GetEntries() )  {
	    h->SetXTitle( "tg theta" );
	    h->SetYTitle( "[mm]" );
	    h->Draw(opt); 
	    gPad->SetLogz();

	    // do the fit
	    double uMin =  h->GetXaxis()->GetBinCenter( 1 );
	    double uMax =  h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
	    buf = new char[64];
	    snprintf(buf, 63, "P1_%i",i);
	    TF1* f = new TF1(buf, Fit::P1, uMin, uMax, 2);
	    SafeDelete(buf);
	    f->SetParameter( 0, 0 );
	    f->SetParameter( 1, 0 );
	    f->SetLineColor( 2 );
	    f->SetLineWidth( 2 );
	    f->SetFillStyle( 0 );
	    TH2Fit* fit = new TH2Fit( f, 2 );
	    fit->Fit( h, uMin, uMax );
	    buf = new char[128];
	    snprintf(buf, 127, "dU=%.3g+%.3g*T", f->GetParameter(0), f->GetParameter(1));
	    cout<<buf<<endl;
	    _PutText( gPad, 0.15, 0.15, 0.06, "dU=%.3g+%.3g*T", f->GetParameter(0), f->GetParameter(1) );
	    f->Draw("same");

	} else cout << endl;

	_UpdateFromPad( gPad );       


    }

    _UpdateFromCanvas( );
    if( ps_ && string( opt ).find("same")==string::npos ) _MakeTitle( "FitDUvsT90 - (%s) %s", detselection, cut );

}

//________________________________________
void CheckTracks::DrawFromOpt( char* file )
{

    //!
    cout << "File is " << file << endl;
    if( access( file, R_OK ) != 0 ) return;
    int argc=2;
    char* argv[] = { "CheckTracks", file, 0 }; 
    Opt *opt = Opt::Instance( argc, argv );

    //! Load detector Table  
    string detFile;
    if( opt->getOpt( "detector", "table", detFile ) ) opt->expand( detFile );
    else {
	cout << "CheckTracks::DrawFromOpt - no detector table.\n";
	return;
    } 

    LoadDetectorFile( detFile.c_str() );
    fflush( stdout );

    //! Load Tracks
    if( tracks_ ) SafeDelete( tracks_ ); 

    bool magnets_on = opt->getOpt("main", "magnets on" );
    tracks_ = new Tracks( magnets_on );

    string fileselection;
    while( opt->getOptRec( "input","tree", fileselection ) ) {
	tracks_->AddToChain( fileselection.c_str(), isBatch_ );
	trackFileSelection_.push_back( fileselection );
	vector< string > files = Utils::GetFiles( fileselection );
	for( unsigned i = 0; i < files.size(); i++ ) trackFiles_.push_back( files[i] );
    }

    if( !tracks_->GetEntries() ) {
	cout << "CheckTracks::DrawFromOpt - No tracks found.\n";
	return;
    }
    fflush( stdout );

    //! Make canvas
    vector<int> geometry;
    if( !( opt->getOpt( "CheckTracks", "geometry", geometry) && geometry.size() >= 2 ) ) {
	geometry.push_back( 700 );
	geometry.push_back( 700 );
    }

    MakeCanvas( geometry[0], geometry[1] );

    //! Make PSFile
    string psFile;
    if( opt->getOpt( "CheckTracks", "psfile", psFile ) ) opt->expand( psFile );
    else psFile = "checkTracks.ps";

    Utils::MakeBackup( psFile );
    OpenPS( psFile.c_str() );
    MakeFrontPage();

    //! Read global cut
    if( !opt->getOpt( "CheckTracks", "globalCut", glCut_ ) ) glCut_ = "";

    //! Dump options
    cout << "===========================================\n";
    cout << "CheckTracks::DrawFromOpt - Magnets        : "   << ((magnets_on)?"on":"off") << ".\n";
    cout << "CheckTracks::DrawFromOpt - detector table : \"" << detFile << "\".\n";
    cout << "CheckTracks::DrawFromOpt - postscript file: \"" << psFile << "\".\n";
    cout << "CheckTracks::DrawFromOpt - geometry       : "   << geometry[0] << "x" << geometry[1] << ".\n";
    cout << "CheckTracks::DrawFromOpt - Global cut     : "   << ((glCut_.size())?glCut_.c_str():"none") << ".\n";
    cout << "CheckTracks::DrawFromOpt - is a batch job : "   << ((isBatch_)? "yes":"no") << ".\n";
    cout << "===========================================\n";
    fflush( stdout );

    //! Make Single plots
    list<string> dOpt;
    while( opt->getOptRec( "CheckTracks","draw", dOpt ) ) 
    {
	unsigned int i=0;
	string var;
	string cut("");
	string opt("");
	int nEnt = 0;
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: var=(*It);  break;
		case 1: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 2: cut=(*It);  break;
		case 3: opt=(*It);  break; 
		default: break;
	    }
	}
	if( opt=="0" ) opt="";
	if( cut=="0" ) cut="";

	string XTitle = (var == "T_duVect" )?"[mm]":"";
	Draw( var.c_str(), cut.c_str(), opt.c_str(), nEnt, XTitle.c_str() );
	fflush( stdout );
    } 

    //! Make detector plots
    while( opt->getOptRec( "CheckTracks","drawDet", dOpt ) ) 
    {
	unsigned int i=0;
	string var;
	string cut("");
	string opt("");
	int nEnt = 0;
	string dets("");
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: var=(*It);  break;
		case 1: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 2: cut=(*It);  break;
		case 3: opt=(*It);  break;
		default: dets+=" "+(*It); break;
	    }
	}
	if( opt=="0" ) opt="";
	if( cut=="0" ) cut="";
	string XTitle = (var == "T_duVect" )?"[mm]":"";
	DrawDet( var.c_str(), dets.c_str(), cut.c_str(), opt.c_str(), nEnt, XTitle.c_str() );
	fflush( stdout );
    }

    //! Make DU Correlation plots
    while( opt->getOptRec( "CheckTracks","drawDU", dOpt ) ) 
    {
	string cut("");
	string dets("");
	int nEnt = 0;
	unsigned int i=0;
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 1: cut=(*It);  break;
		default: dets+=" "+(*It); break;
	    }
	}
	if( cut=="0" ) cut="";
	DrawDU( dets.c_str(), cut.c_str(), nEnt );
	fflush( stdout );
    }

    //! Make DU_LR plots
    while( opt->getOptRec( "CheckTracks","drawDU_LR", dOpt ) ) 
    {
	string cut("");
	string dets("");
	int nEnt = 0;
	unsigned int i=0;
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 1: cut=(*It);  break;
		default: dets+=" "+(*It); break;
	    }
	}
	if( cut=="0" ) cut="";
	DrawDU_LR( dets.c_str(), cut.c_str(), nEnt );
	fflush( stdout );
    }

    //! Make DUvsU Correlation plots
    while( opt->getOptRec( "CheckTracks","drawUCor", dOpt ) ) 
    {
	string cut("");
	string opt("");
	string dets("");
	int nEnt = 0;
	unsigned int i=0;
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 1: cut=(*It);  break;
		case 2: opt=(*It);  break; 
		default: dets+=" "+(*It); break;
	    }
	}
	if( opt=="0" ) opt="";
	if( cut=="0" ) cut="";
	DrawDUvsU( dets.c_str(), cut.c_str(), opt.c_str(), nEnt );
	fflush( stdout );
    }

    //! Fit DUvsU Correlation plots
    while( opt->getOptRec( "CheckTracks","fitUCor", dOpt ) ) 
    {
	string cut("");
	string opt("");
	string dets("");
	int nEnt = 0;
	unsigned int i=0;
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 1: cut=(*It);  break;
		case 2: opt=(*It);  break; 
		default: dets+=" "+(*It); break;
	    }
	}
	if( opt=="0" ) opt="";
	if( cut=="0" ) cut="";
	FitDUvsU( dets.c_str(), cut.c_str(), opt.c_str(), nEnt );
	fflush( stdout );
    }

    //! Make DUvsV Correlation plots
    while( opt->getOptRec( "CheckTracks","drawVCor", dOpt ) ) 
    {
	unsigned int i=0;
	string cut("");
	string opt("");
	int nEnt = 0;
	string dets("");
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 1: cut=(*It);  break;
		case 2: opt=(*It);  break; 
		default: dets+=" "+(*It); break;
	    }
	}
	if( opt=="0" ) opt="";
	if( cut=="0" ) cut="";
	DrawDUvsV( dets.c_str(), cut.c_str(), opt.c_str(), nEnt );
	fflush( stdout );
    }

    //! Fit DUvsV Correlation plots
    while( opt->getOptRec( "CheckTracks","fitVCor", dOpt ) ) 
    {
	unsigned int i=0;
	string cut("");
	string opt("");
	int nEnt = 0;
	string dets("");
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 1: cut=(*It);  break;
		case 2: opt=(*It);  break; 
		default: dets+=" "+(*It); break;
	    }
	}
	if( opt=="0" ) opt="";
	if( cut=="0" ) cut="";
	FitDUvsV( dets.c_str(), cut.c_str(), opt.c_str(), nEnt );
	fflush( stdout );
    }

    //! Make DUvsT Correlation plots
    while( opt->getOptRec( "CheckTracks","drawTCor", dOpt ) )   
    {
	unsigned int i=0;
	string cut("");
	string opt("");
	int nEnt = 0;
	string dets("");
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 1: cut=(*It);  break;
		case 2: opt=(*It);  break;
		default: dets+=" "+(*It); break;
	    }
	}
	if( opt=="0" ) opt="";
	if( cut=="0" ) cut="";
	DrawDUvsT( dets.c_str(), cut.c_str(), opt.c_str(), nEnt );
	fflush( stdout );
    }

    //! Fit DUvsT Correlation plots
    while( opt->getOptRec( "CheckTracks","fitTCor", dOpt ) ) 
    {
	unsigned int i=0;
	string cut("");
	string opt("");
	int nEnt = 0;
	string dets("");
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: sscanf( (*It).c_str(), "%i", &nEnt );  break;
		case 1: cut=(*It);  break;
		case 2: opt=(*It);  break; 
		default: dets+=" "+(*It); break;
	    }
	}
	if( opt=="0" ) opt="";
	if( cut=="0" ) cut="";
	FitDUvsT( dets.c_str(), cut.c_str(), opt.c_str(), nEnt );
	fflush( stdout );
    }

    //! Make Track Multiplicity plots
    while( opt->getOptRec( "CheckTracks","drawTrMlt", dOpt ) ) 
    {
	unsigned int i=0;
	string cut("");
	string opt("");
	int nEnt = 0;
	string dets("");
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: s >> nEnt;  break;
		case 1: cut=(*It);  break;
		case 2: opt=(*It);  break; 
		default: dets+=" "+(*It); break;
	    }
	}
	if( cut=="0" ) cut="";
	if( opt=="0" ) opt="";
	DrawTrackMlt( dets.c_str(), cut.c_str(), opt.c_str(), nEnt );
	fflush( stdout );
    }

    //! Make DUvsZ Correlation plots
    while( opt->getOptRec( "CheckTracks","drawZCor", dOpt ) )
    {
	unsigned int i=0;
	string cut("");
	int nEnt = 0;
	string dets("");
	for( list<string>::iterator It = dOpt.begin(); It != dOpt.end(); It++, i++ ) {
	    istringstream s( (*It).c_str(), ios::in );
	    switch( i ) {
		case 0: s >> nEnt;  break;
		case 1: cut=(*It);  break;
		default: dets+=" "+(*It); break;
	    }
	}
	if( cut=="0" ) cut="";
	DrawDUvsZ( dets.c_str(), cut.c_str(), nEnt );
	fflush( stdout );
    }

    //! Make bit pattern plots
    list<string> zones;
    while( opt->getOptRec("CheckTracks","bitPattern", zones ) ) {

	int zoneMask=0;
	if( !zones.size() ) zoneMask = 0xf;
	else for( list<string>::iterator I=zones.begin(); I!=zones.end(); I++ ) {
	    int zone;
	    istringstream s( (*I).c_str(), ios::in ); s >> zone;
	    zoneMask += (1<<zone) ;
	}

	if( !cv_ ) MakeCanvas();
	cv_->Clear();
	_NewPage();
	Macro::BitPattern( trackFiles_, cv_, zoneMask );
	_UpdateFromCanvas();
	_MakeTitle( "BitPattern - (0x%x)", zoneMask );

    }

    _NewPage();
    ClosePS();
    cout << "CheckTracks::DrawFromOpt - Done.\n";
    return;

}

//_________________________________________________________
TGraphErrors* CheckTracks::_AddTGE( int color )
{
    TGraphErrors *tg = new TGraphErrors();
    tg->SetMarkerColor(color);
    tg->SetLineColor(color);
    tg->SetLineWidth(1);
    tg->SetMarkerStyle(20);
    tg->SetMarkerSize(1);
    return tg;
}

//________________________________________________________________________________
void CheckTracks::_NewPage( void )  
{ 
    if( !ps_ ) return;
    ps_->On();
    ps_->NewPage();
    ps_->Off();
}

//________________________________________________________________________________
void CheckTracks::_UpdateFromCanvas( void )  
{ 
    if( !( ps_ && cv_ ) ) return;
    ps_->On();
    cv_->Update();
    ps_->Off();
}

//________________________________________________________________________________
void CheckTracks::_UpdateFromPad( TVirtualPad* pad )  
{ 
    if( !( ps_ && pad ) ) return;
    ps_->On();
    pad->Update();
    ps_->Off();
}

//________________________________________________________________________________
void CheckTracks::_PutText( TVirtualPad* pad, double x, double y, double size, const char* format, ... )  
{ 
    if( !pad ) return;
    char* buf = new char[128];
    va_list p;
    va_start(p,format);
    vsprintf(buf, format, p);
    va_end(p);

    TText *text = new TText( );
    text->SetTextColor( 2 );
    text->SetTextSize( size );

    if( ps_ ) ps_->On();
    text->DrawTextNDC( x, y, buf );
    if( ps_ ) ps_->Off();

    SafeDelete( buf );
}  


//________________________________________________________________________________
void CheckTracks::_MakeTitle( const char* format, ... )  
{ 
    char* buf = new char[128];
    va_list p;
    va_start(p,format);
    vsprintf(buf, format, p);
    va_end(p);

    if( !(ps_&&cv_) ) {
      SafeDelete(buf);
      return;
    }
    cv_->cd();
    TPaveText *text = new TPaveText(0.2,0.95,0.8,1.0,"NDC");
    text->SetBorderSize( 0 );
    text->SetLineWidth( 0 );
    text->SetFillStyle( 0 );
    text->AddText( buf );
    ps_->On();
    text->Paint();
    ps_->Off();
    SafeDelete(buf);
}    
