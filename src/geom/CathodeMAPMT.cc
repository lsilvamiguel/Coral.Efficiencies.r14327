/*!
   \file    
   \brief   
   \author  
   \version 
   \date    
*/

//////////////////// RICH1MAPMT FUNCTIONS ///////////////////////////////////////////

#include "CsOpt.h"
#include "CathodeMAPMT.h"
#include "CsMCRICH1Hit.h"
#include "RayTrace.h"

using namespace std;



//////////////////////////////////////////////////////////////////////////////////

CathodeMAPMT::CathodeMAPMT (int id ) : CsRICH1UpGrade::CathodePlane( id ) ,
				       use_quantum_efficiency_corrections_(false),
				       use_obsolete_lens_numbering_(false),
				       use_old_zebra_format_(false),
				       rescale_wavelength_(false),
				       peak_eff_(1.)
{
  read_out_type_ = MAPMT;
  pad_size_ = 12.;
  InitCsOpt();
}

//////////////////////////////////////////////////////////////////////////////////

void CathodeMAPMT::InitCsOpt( void )
{
  CsRICH1UpGrade::CathodePlane::InitCsOpt();

  CsOpt* opt = CsOpt::Instance();

  float f =  peak_eff_;
  opt->CsOpt::getOpt( "MCSettingsCathodeMAPMT", "PeakEff", f );
  if(f>0. && f <=1.)    peak_eff_=f;
  else cerr << "CathodeMAPMT::InitCsOpt: peak_eff_ not between 0 and 1. Setting it to 1." << endl;

  string str;

  if ( opt->CsOpt::getOpt( "MCSettingsCathodeMAPMT", "UseQuantEffCorrections", str ) )
    {
      if( sizeof(str) > 0 )
	{
	  if( str[0] == 'Y' )
	    {
	      cout << "CathodeMAPMT: Using efficiency correction with peak efficiency set to  " << peak_eff_ << endl;
	      use_quantum_efficiency_corrections_ = true;
	    }
	}
    }

  if ( opt->CsOpt::getOpt( "MCSettingsCathodeMAPMT", "UseObsoleteLensNumbering", str ) )
    {
      if( sizeof(str) > 0 )
	{
	  if( str[0] == 'Y' )
	    {
	      cout << "CathodeMAPMT: Using obsolete lens numbering scheme!" << endl;
	      use_obsolete_lens_numbering_ = true;
	    }
	}
    }

  if ( opt->CsOpt::getOpt( "MCSettingsCathodeMAPMT", "UseOldZebraFormat", str ) )
    {
      if( sizeof(str) > 0 )
	{
	  if( str[0] == 'Y' )
	    {
	      cout << "CathodeMAPMT: Using old zebra format (i.e. without lenses)." << endl;
	      use_old_zebra_format_ = true;
	    }
	}
    }

  if ( opt->CsOpt::getOpt( "MCSettingsCathodeMAPMT", "RescaleWavelength", str ) )
    {
      if( sizeof(str) > 0 )
	{
	  if( str[0] == 'Y' )
	    {
	      cout << "CathodeMAPMT: Rescaling wavelength to PMT sensitive domain!" << endl;
	      cout << "              This overrides the efficiency correction!" << endl;
	      rescale_wavelength_ = true;
	    }
	}
    }
}

//////////////////////////////////////////////////////////////////////////////////

void CathodeMAPMT::Clear( void )
{
  CsRICH1UpGrade::CathodePlane::Clear();
  pseudopads_.clear(); 
}

//////////////////////////////////////////////////////////////////////////////////

void CathodeMAPMT::MakeMCResponse( void )
{
  //  bool debug = true;
  bool debug = false;
  if( debug ) cout << " Well Start CathodeMAPMT::MakeMCResponse in debug " << endl;

  static OptSyst Telescope; 
  
  SortMCHits();
  for ( unsigned ih=0; ih< NMCHits(); ih++ )
    {
      const CsMCRICH1Hit &hit=  *mchits_[ih];
//      double x = hit.getYdetDRS();
//      double y = hit.getZdetDRS();
// convert x and y to cm
      double x = hit.getYdetDRS()/10.;
      double y = hit.getZdetDRS()/10.;
      double t = hit.getDTime();
      double e = hit.getPhotEnergy();
//      double xmrs = hit.getX();
//      double ymrs = hit.getY();
//      double zmrs = hit.getZ();
//      double xref = hit.getXref();
//      double yref = hit.getYref();
//      double zref = hit.getZref();
// reshuffle coordinates to correspond to COMGEANT MRS
// Everything in mm, but this should be ok, as we are only
// interested in the direction...
      double ymrs = hit.getX();
      double zmrs = hit.getY();
      double xmrs = hit.getZ();
      double yref = hit.getXref();
      double zref = hit.getYref();
      double xref = hit.getZref();
      int dn = hit.getLens();
      int dnc = hit.getCathode();
      if( debug ) cout << " Try to make MAPT signal at  x " << x << "  y " << y << " e " << e ;
      if( debug ) cout << "   lens  " << dn << " cath  " << dnc ;
      if( debug ) cout << "   xmrs " << xmrs << "  ymrs " << ymrs << " zmrs " << zmrs;
      if( debug ) cout << "   xref " << xref << "  yref " << yref << " zref " << zref << endl;


//       cerr <<  x << "  " << y << "  " << e ;
//       cerr << " " << dn << "  " << dnc ;
//       cerr << "  " << xmrs << "   " << ymrs << "  " << zmrs;
//       cerr << " " << xref << "  " << yref << " " << zref << endl ;

//       if(dn < 17) cerr << "Not on central cathode? dn = " << dn ;
      
      LightRay photon;
      int icath = -1;
      int lens = -1;
      int ix_lens, iy_lens, ix_pin, iy_pin;
            
      ix_lens = -1;
      iy_lens = -1;
      ix_pin  = -1;
      iy_pin  = -1;      

      if(use_old_zebra_format_){
	icath=dnc;
	lens=Telescope.FindFirstLens(xref,  yref,  zref,  xmrs,  ymrs,  zmrs, x, y, e, photon);
	x=photon.getval(1)/10.;
	y=photon.getval(2)/10.;
	if(icath > 8) y=-y;
	if(icath < 9) x=-x;
	if ( debug ){
	  cout << "----> Photon parameters for Raytrace: " ;
	  photon.print();
	  cout << "----> Cathode, lens number and coordinates on first lens: " << icath << " " << lens << " " << x << " " << y << endl; 
	}
      }
      
      else{
	// Set photon parameters in the local lens system
	icath=Telescope.TransPhotonDir( xref,  yref,  zref,  xmrs,  ymrs,  zmrs, x, y, e, dn, photon);
	if(icath < 0) continue; 
	//  if(icath != dnc ) cerr << "Different opinion about Cathode number... icath = " << icath << "  dnc = " << dnc << endl;
	
	//  simple (=new) scheme 
	if(icath ==  4) lens = dn -  16;
	if(icath ==  6) lens = dn - 160;
	if(icath == 11) lens = dn - 304;
	if(icath == 13) lens = dn - 448;
      }
      
      if(lens < 1 || lens > 144){
	cerr << "CathodeMAPMT::MakeMCResponse: Lens number outside allowed range:  " << lens << endl;
	continue;
      }      
      
      // scheme used by Luis in the first versions of the zebra file
      if(use_obsolete_lens_numbering_){
	if(icath ==  4){
	  iy_lens = 11 - (lens-1)/12;
	  ix_lens = (lens-1) - (11-iy_lens)*12;
	}
	if(icath ==  6){ 
	  iy_lens = 11 - (lens-1)/12;
	  ix_lens = 11 - (lens-1) + (11-iy_lens)*12;
	}
	if(icath == 11){ 
	  iy_lens = (lens-1)/12;
	  ix_lens = (lens-1) - iy_lens*12;		
	}
	if(icath == 13){
	  iy_lens = (lens-1)/12;
	  ix_lens = 11 - (lens-1) + iy_lens*12;		
	}
      }
      else{
	iy_lens = (lens-1)/12;
	ix_lens = (lens-1) - iy_lens*12;
      }

      
      // find the pseudo-pad from the lens coordinates
      double ysize = 47.7/40.;
      double xsize = 44.8/40.;
      int ix=(int)floor((x+2.*xsize)/xsize);
      int iy=(int)floor((y+2.*ysize)/ysize);
      int ix_pseudo = ix_lens*4 + ix ; 
      int iy_pseudo = iy_lens*4 + iy ; 
      
      // rescale the wavelength to the sensitive domain, if used on old zebra files. Ignore efficiency correction!
      if( rescale_wavelength_ ){
	if(debug) cout << "----> Original wavelength    " << photon.getval(7) << endl; 
	photon.SetWave(9.5*photon.getval(7)-1320);
	if(debug) cout << "----> Rescaled wavelength to " << photon.getval(7) << endl; 
      }
      else{
	// set to 1 if photon is outside PMT efficiency region
	if( use_quantum_efficiency_corrections_) {
	  int notdet = 0;
	  notdet = Telescope.CheckEff(photon,peak_eff_);
	  if(notdet) {
	    //        cerr << "Wavelength not in sensitive PMT region" << endl;
	    continue;
	  }
	}    
      }
     

      // The following should be the same for all options! 
      
      // If we get here, we finally trace the photon through the telescope
      int ipin  = -1;
      ipin=Telescope.RayTrace(photon);
      //      cerr << "Photon traced to pin " << ipin << endl;
      if(ipin < 0) continue;
      
      iy_pin = ipin/4;
      ix_pin = ipin - iy_pin*4;
      
      // account for optics
      ix_pin = 3-ix_pin;
      iy_pin = 3-iy_pin;
      
      // account for different orientation of ref. frames
      if(icath ==  4 || icath ==  6) ix_pin = 3 - ix_pin;
      if(icath == 11 || icath == 13) iy_pin = 3 - iy_pin;
      
      int ix_pad = ix_lens*4+ix_pin;
      int iy_pad = iy_lens*4+iy_pin;
      
      //      cerr << "Successfully determined pseudo-pads " << ix_pseudo << " " << iy_pseudo << endl; 
      
      double amp = 10. ;
      
      if( debug ) cout << "----> Pseudo-pads: " << ix_pseudo << " " << iy_pseudo << " " << ix << " " << iy << endl;
      if( debug ) cout << "====> Pads       : " << ix_pad << " " << iy_pad << endl;
      
      // I'm not sure which value id_ will have, therefore I set it "by hand"
      pseudopads_.push_back( CsRICH1UpGrade::CathodePAD( icath-1, ix_pseudo, iy_pseudo , amp, t) );
      pads_.push_back( CsRICH1UpGrade::CathodePAD( icath-1, ix_pad, iy_pad , amp, t) );
      
      
    }
  if( debug ) cout << " CathodeMAPMT::MakeMCResponse in debug OK " << endl;
}

//////////////////////////////////////////////////////////////////////////////////

void    CathodeMAPMT::FillMCDecodingHisto( void )
{
  if( !flags_.histo_mc_decoding_booked_ )
    {
      flags_.histo_mc_decoding_booked_=true;
      char path[132],hist_name[132];
      sprintf(path,"/CsRICH1UpGrade/MCHits/Cathode_%d",id_);
      CsHistograms::SetCurrentPath(path);

      histo_.h1_Time = new  CsHist1D("Time"," Time ",500,-50.,50.);

      sprintf(hist_name," Cathode %d:  Y vs X Pads ",id_+1);
      histo_.h2_padXY= new  CsHist2D("padXY",hist_name,48,0.,48.,48,0.,48.);

      sprintf(hist_name," Cathode %d:  Y vs X Pseudo pads ",id_+1);
      histo_.h2_pseudoXY= new  CsHist2D("pseudoXY",hist_name,48,0.,48.,48,0.,48.);

      sprintf(hist_name," Cathode %d: pad - pseudo in X ",id_+1);
      histo_.h1_diffX = new  CsHist1D("diffX",hist_name,10,-5.,5.);

      sprintf(hist_name," Cathode %d: pad - pseudo in Y ",id_+1);
      histo_.h1_diffY = new  CsHist1D("diffY",hist_name,10,-5.,5.);
    }
  
  for ( unsigned i=0; i< pads_.size(); i++ )
    {
      int ix = pads_[i].ix_;
      int iy = pads_[i].iy_;
      double time = pads_[i].time_;
      histo_.h1_Time->Fill( time );
      histo_.h2_padXY->Fill( double(ix), double(iy) );
      int ixp = pseudopads_[i].ix_;
      int iyp = pseudopads_[i].iy_;
      histo_.h2_pseudoXY->Fill( double(ixp), double(iyp) );
      histo_.h1_diffX->Fill( double(ix-ixp) );
      histo_.h1_diffY->Fill( double(iy-iyp) );
    }
}

//////////////////////////////////////////////////////////////////////////////////
