/*!
   \file    
   \brief   
   \author  Kolosov Vladimir            Vladimir.Kolosov@cern.ch  
   \version 
   \date    
*/

////////////////// CathodeAPV Geometry ///////////////////////

/*
           ==<==  FE numbering dir 0-11,
             /\   APV numbering dir 0-3  
    ====================================                     
    I  ===<==  ===<==  ===<==  ===<==  I
    I  I    I  I    I  I    I  I    I  I
    I  I  0\/  I  2\/  I  4\/  I  6\/  I
    I  I    I  I    I  I    I  I    I  I
    I  ======  ======  ======  ======  I
    I                                  I             
    I  ======  ======  ======  ======  I                 Y
    I  I    I  I    I  I    I  I    I  I              /\
    I /\  1 I /\  3 I /\  5 I /\  7 I  I   Coral RF    |      Z   
    I  I    I  I    I  I    I  I    I  I               |     /
    I  ==>===  ==>===  ==>===  ==>===  I               |   /
    ====================================               | /
 Jura                                  Saleve    X<----/
    ====================================
    I  ===<==  ===<==  ===<==  ===<==  I
    I  I    I  I    I  I    I  I    I  I
    I  I  8\/  I 10\/  I 12\/  I 14\/  I
    I  I    I  I    I  I    I  I    I  I
    I  ======  ======  ======  ======  I
    I                                  I
    I  ======  ======  ======  ======  I
    I  I    I  I    I  I    I  I    I  I
    I /\  9 I /\ 11 I /\ 13 I /\ 15 I  I
    I  I    I  I    I  I    I  I    I  I
    I  ==>===  ==>===  ==>===  ==>===  I
    ====================================
*/

////////////////// CathodeAPV FUNCTIONS //////////////////////////////////////////
#include "CsOpt.h"
#include "CathodeAPV.h"
#include "CsMCRICH1Hit.h"
#include <cassert>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////

const double ElectronicAPV::ballistic_deficit_ = 0.43;

//////////////////////////////////////////////////////////////////////////////////

CathodeAPV::CathodeAPV ( int id ) : 
 CsRICH1UpGrade::CathodePlane( id ), 
 use_quantum_efficiency_corrections_(false),
 peak_eff_(1.),
 rd_read_out_dirx_(1),
 rd_read_out_diry_(1)
{
//   bool debug = false;
  bool debug = true;
  read_out_type_ = APV;
  nwires_=333;
  wire0_position_ = -290.; // Outside first pad?
//  wire0_position_ = -286.;
  wire_pitch_size_ = 4.;
  
//   electron_drift_velocity_ = 1./30.;
  electron_drift_velocity_ = 1./10.;
  multiplication_gain_ = 15000.;
  InitCsOpt();
//   apv_.PrintSettings();

  SetMCReadoutParmeters();

//   ElectronicAPV apv;
//   apv.TestSignal();

// Set readout directions for RD decoding
  if( id_<0 || id_ > 15 ) assert(false);
  SetRDreadoutDirs( 1,1 );

//  SetRDreadoutDirs( -1,-1 );

//   int xcor = -1;
//   int ycor = -1;
//   switch (id_)
//   {
//     case 15:
//      SetRDreadoutDirs( -1*xcor, 1*ycor );
//      break;
//     case 14:
//      SetRDreadoutDirs(  1*xcor,-1*ycor );
//      break;
//     case 13:
//      SetRDreadoutDirs( -1*xcor, 1*ycor );
//      break;
//     case 12:
//      SetRDreadoutDirs(  1*xcor,-1*ycor );
//      break;
//     case 11:
//      SetRDreadoutDirs( -1*xcor, 1*ycor );
//      break;
//     case 10:
//      SetRDreadoutDirs(  1*xcor,-1*ycor );
//      break;
//     case 9:
//      SetRDreadoutDirs( -1*xcor, 1*ycor );
//      break;
//     case 8:
//      SetRDreadoutDirs(  1*xcor,-1*ycor );
//      break;
//     case 7:
//      SetRDreadoutDirs( -1*xcor, 1*ycor );
//      break;
//     case 6:
//      SetRDreadoutDirs(  1*xcor,-1*ycor );
//      break;
//     case 5:
//      SetRDreadoutDirs( -1*xcor, 1*ycor );
//      break;
//     case 4:
//      SetRDreadoutDirs(  1*xcor,-1*ycor );
//      break;
//     case 3:
//      SetRDreadoutDirs( -1*xcor, 1*ycor );
//      break;
//     case 2:
//      SetRDreadoutDirs(  1*xcor,-1*ycor );
//      break;
//     case 1:
//      SetRDreadoutDirs( -1*xcor, 1*ycor );
//      break;
//     case 0:
//      SetRDreadoutDirs(  1*xcor,-1*ycor );
//      break;
//   }
}

//////////////////////////////////////////////////////////////////////////////////

void CathodeAPV::SetMCReadoutParmeters( void )
{
  for( unsigned ixpad = 0; ixpad < nxpads_; ixpad++)
  {
    double r = 0.; 
    for( unsigned iypad = 0; iypad < nypads_; iypad++)
    {
      double amp_unit = apv_.amp_unit_;
      r = CsRICH1UpGrade::Randm();
      amp_unit += 0.1*apv_.amp_unit_*(r-0.5);
      double noise = apv_.noise_;
      r = CsRICH1UpGrade::Randm();
      noise += 0.2*apv_.noise_*(r-0.5);
      double ped = apv_.pedestal_;
      r = CsRICH1UpGrade::Randm();
      ped += amp_unit*r;
      ch_apv_[ixpad][iypad].Set( amp_unit, noise, apv_.signal_cut_in_sigmas_, ped );
    } 
  } 
}

//////////////////////////////////////////////////////////////////////////////////

void CathodeAPV::InitCsOpt( void )
{
   CsRICH1UpGrade::CathodePlane::InitCsOpt();
   CsOpt* opt = CsOpt::Instance();
   float f = multiplication_gain_;
   opt->CsOpt::getOpt( "MCSettingsCathodeAPV", "MultiGain", f );
   multiplication_gain_=f;

   f =  peak_eff_;
   opt->CsOpt::getOpt( "MCSettingsCathodeAPV", "PeakEff", f );
   if(f>0. && f <=1.)    peak_eff_=f;
   else cerr << "CathodeAPV::InitCsOpt: peak_eff_ not between 0 and 1. Setting it to 1." << endl;

   string str;
   if ( opt->CsOpt::getOpt( "MCSettingsCathodeAPV", "UseQuantEffCorrections", str ) )
   {
     if( sizeof(str) > 0 )
     {
       if( str[0] == 'Y' )
       {
         use_quantum_efficiency_corrections_ = true;
       }
     }
   }
}

//////////////////////////////////////////////////////////////////////////////////

void    CathodeAPV::Clear( void )
{

  CsRICH1UpGrade::CathodePlane::Clear(); 

  mc_avalanches_.clear();
  for( unsigned ix=0; ix< nxpads_; ix++)
    for( unsigned iy=0; iy< nypads_; iy++)
    {
      signals_[ix][iy].clear();
    }
  data_apv_.clear();
}

//////////////////////////////////////////////////////////////////////////////////

void    CathodeAPV::FillMCDecodingHisto( void )
{
//   if( !flags_.histo_mc_decoding_booked_ )
//   {
//     flags_.histo_mc_decoding_booked_=true;
//     char path[132],hist_name[132];
//     sprintf(path,"/CsRICH1UpGrade/MCHits/Cathode_%d",id_);
//     CsHistograms::SetCurrentPath(path);
//     histo_.h1_A2 = new  CsHist1D("A2"," APV A2 ",500,0.,500.);
//     histo_.h1_rA1A2 = new  CsHist1D("rA1A2"," APV Ratio A1/A2 ",100,0.,2.);
//     histo_.h1_rA1A2cut = new  CsHist1D("rA1A2cut"," APV Ratio A1/A2 if A2>10 ",100,0.,2.);
//     histo_.h2_A2A1 = new  CsHist2D("A2A1"," APV A1 vs A2 ",50,0.,200.,50,0.,200.);
//     histo_.h2_A2A0 = new  CsHist2D("A2A0"," APV A0 vs A2 ",50,0.,200.,50,0.,200.);
//   }
// 
//   for( unsigned it=0; it< NPadsFired(); it++ )
//   {
//     CsRICH1UpGrade::CathodePAD &pad = GetPads()[it];
//     int ixpad =  pad.ix_;
//     int iypad =  pad.iy_;
//     double amp =  pad.amp_;
//     double time =  pad.time_;
//     assert ( pad.el_digits_.size() == 3 );
//     double a0 = (double)pad.el_digits_[0];
//     double a1 = (double)pad.el_digits_[1];
//     double a2 = (double)pad.el_digits_[2];
//     histo_.h1_A2->Fill(a2);
//     if( a2 != 0 ) histo_.h1_rA1A2->Fill(a1/a2);
//     if( a2 > 10. ) histo_.h1_rA1A2cut->Fill(a1/a2);
//     histo_.h2_A2A1->Fill(a2,a1);
//     histo_.h2_A2A0->Fill(a2,a0);
//   }
  static vector< double > st;
  if( !flags_.histo_mc_decoding_booked_ )
  {
    flags_.histo_mc_decoding_booked_=true;
    char path[132],hist_name[132];
    sprintf(path,"/CsRICH1UpGrade/MCHits/Cathode_%d",id_);
    CsHistograms::SetCurrentPath(path);
    histo_.h1_A2 = new  CsHist1D("A2"," APV A2 ",500,0.,500.);
    histo_.h1_rA1A2 = new  CsHist1D("rA1A2"," APV Ratio A1/A2 ",100,0.,2.);
    histo_.h1_rA1A2cut = new  CsHist1D("rA1A2cut"," APV Ratio A1/A2 if A2>10 ",100,0.,2.);
    histo_.h2_A2A1 = new  CsHist2D("A2A1"," APV A1 vs A2 ",255,0.,255.,255,0.,255.);
    histo_.h2_A2A0 = new  CsHist2D("A2A0"," APV A0 vs A2 ",64,0.,63.,255,0.,255.);
    histo_.h2_A2A1_A0corr = new  CsHist2D("A2A1_A0corr"," APV A1-A0 vs A2-A0 ",300,-100.,200.,300,-100.,200.);
    histo_.h2_Chisq_Afit = new  CsHist2D("Chisq_Afit"," APV Afit vs Chisq ",200,0.,100.,200,0.,200.);
    st = apv_.ReadTestSignal();
    st[0]=st[0]/st[2];
    st[1]=st[1]/st[2];
    st[2]=1.;
  }

  for( unsigned it=0; it< NPadsFired(); it++ )
  {
    CsRICH1UpGrade::CathodePAD &pad = GetPads()[it];
    int ixpad =  pad.ix_;
    int iypad =  pad.iy_;
    double amp =  pad.amp_;
    double time =  pad.time_;
    assert ( pad.el_digits_.size() == 3 );
    double a0 = (double)pad.el_digits_[0];
    double a1 = (double)pad.el_digits_[1];
    double a2 = (double)pad.el_digits_[2];
    histo_.h1_A2->Fill(a2);
    if( a2 != 0 ) histo_.h1_rA1A2->Fill(a1/a2);
    if( a2 > 10. ) histo_.h1_rA1A2cut->Fill(a1/a2);
    histo_.h2_A2A1->Fill(a2,a1);
    histo_.h2_A2A0->Fill(a2,a0);
    histo_.h2_A2A1_A0corr->Fill(a2-a0,a1-a0);
 //  Fit with shape a0,a1,a2   
    {
//       vector< double > st = apv_.ReadTestSignal();
// //       cout << " Ideal signal: a0t " <<  st[0] << " a1t " <<  st[1] << " a2t " <<  st[2] << endl;
// //       cout << " Real signal: a0 " << a0 << " a1 " <<  a1 << " a2 " <<  a2 << endl;
//       st[0]=st[0]/st[2];
//       st[1]=st[1]/st[2];
//       st[2]=1.;
      double sumtsq =  st[0]*st[0]+st[1]*st[1]+st[2]*st[2];
      double sumata =  st[0]*a0 + st[1]*a1 + st[2]*a2;
      double afit = sumata/sumtsq;
      if( afit < 0 ) afit = 0.;
      double da0 = (a0 - afit*st[0]);
      double da1 = (a1 - afit*st[1]);
      double da2 = (a2 - afit*st[2]);
      double chisq = da0*da0 + da1*da1  + da2*da2;
      double disp = ch_apv_[ixpad][iypad].noise_/ch_apv_[ixpad][iypad].amp_unit_;
      disp *= disp;
      chisq = chisq/disp/2.;
//       cout << " Ch_APV  x = " << ixpad << " y= " << iypad << " noise " << ch_apv_[ixpad][iypad].noise_ << " disp " << disp << endl;
//       cout << " chisq " << chisq << " afit " << afit << endl;
      histo_.h2_Chisq_Afit->Fill(chisq,afit);
    }

  }
  
  
}

//////////////////////////////////////////////////////////////////////////////////

void    CathodeAPV::MakeMCResponse( void )
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " CathodeAPV::MakeMCResponse in " << id_ << endl;
  SortMCHits();
  for ( unsigned ih=0; ih<mchits_.size(); ih++ )
  {
    bool makeit = MakeAvalanche( mchits_[ih] );
  }
  MakeElectronicMCDigits();
}

//////////////////////////////////////////////////////////////////////////////////

bool    CathodeAPV::MakeAvalanche( const CsMCRICH1Hit *hit )
{
//   bool debug = true;
  bool debug = false;

  double x = hit->getYdetDRS();
  double y = hit->getZdetDRS();
  double t = hit->getDTime();
  double e = hit->getPhotEnergy();

  if( debug )
  {
    cout << " Try to make Avalanche at  x " << x << "  y " << y << " e " << e << endl;
  }
  if( use_quantum_efficiency_corrections_ )
  {
    double eff = CathodeQuantumEfficiency( e, x, y);
    if( eff< 0. || eff > 1. ) 
    {
       cerr << " CathodeAPV::MakeAvalanche:: Warning wrong efficiency value e " << e << " eff " << eff << endl;
    }
    eff *= efficiency_reduction_factor_;
    double r = CsRICH1UpGrade::Randm();
    if( debug )
    {
      cout << " QuantEff " << eff;
      if( r > eff ) 
        cout << " No photo-electron " << endl;
      else 
        cout << " photo electron created " << endl;
    }
    if( r > eff ) return false;
  }

  int iwire = (int)(( x - wire0_position_ )/wire_pitch_size_);
  double xwire = wire0_position_ + iwire*wire_pitch_size_;
  double dx = x - xwire;
  if( dx > wire_pitch_size_/2. )
  {
    iwire++;
    xwire += wire_pitch_size_;
    dx -= wire_pitch_size_;    
  }
  pair < int, int> pad_indexes = GetPadIndexes( xwire, y );
  if( pad_indexes.first  < 0 ) return false;

  int ixpad = pad_indexes.first; 
  int iypad = pad_indexes.second;
  double drift_time = fabs(dx)/electron_drift_velocity_;
  Avalanche ava( xwire, y, 1., t + drift_time, multiplication_gain_ );
  mc_avalanches_.push_back( ava );
  
  if( debug ) cout << " xwire " << xwire << " y " << y <<  endl;
  
  double xpad =  pad0_x_position_ + ixpad*pad_size_;
  double ypad =  pad0_y_position_ + iypad*pad_size_;
  if( debug ) cout << " xpad " << xpad << " ypad " << ypad <<  endl;
  double dxpad = xwire - xpad;
  double dypad = y - ypad;
  
  if( debug ) cout << " Ava " << mc_avalanches_.back().amp_ << " dxpad " << dxpad <<  " dypad " << dypad << endl;
  for( int ix=-2; ix<=2; ix++)  
    for( int iy=-2; iy<=2; iy++)
    {
      double epad = mc_avalanches_.back().amp_*PadResponse(dxpad,dypad,ix,iy);
      int ixpads = ixpad + ix;
      int iypads = iypad + iy;
      if( ixpads >= 0 && ixpads < (int)nxpads_ && iypads >= 0 && iypads < (int)nypads_ )
      {
        signals_[ixpads][iypads].push_back( CathodeAPV::Signal(epad,mc_avalanches_.back().time_,hit) );
      }
    }


  if( debug )
  {
    cout << " Hit Wire N = " << iwire  << "  dx " << dx << 
            "  drift time " << drift_time << "  Amplitude " << mc_avalanches_.back().amp_ << endl;
    if( iwire < 0 || fabs(dx) > wire_pitch_size_/2.+0.05 )
    {
      cerr << " CathodeAPV::MakeAvalanche bug in geometry!! " << 
                " Wire N = " << iwire  << "  dx " << dx << "  drift time " << drift_time << endl;
      exit(1);          
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////////

void    CathodeAPV::MakeElectronicMCDigits( void )
{
// //   bool debug = true;
//   bool debug = false;
//   if( debug ) cout << " CathodeAPV::MakeElectronicMCDigits id " <<  id_ << endl;
//   
//   for( unsigned ixpad = 0; ixpad < nxpads_; ixpad++)
//   {  
//     for( unsigned iypad = 0; iypad < nypads_; iypad++)
//     {
//       apv_.Clear();
//       for( unsigned is=0; is< signals_[ixpad][iypad].size(); is++ )
//       {
//         apv_.AddSignal( signals_[ixpad][iypad][is].amp_, signals_[ixpad][iypad][is].time_ );
//       }
//       double amp = 0.;
//       vector < int > amps = apv_.Read( ch_apv_[ixpad][iypad] );
//       if( amps.size() == 0 )
//       {
// //         if( debug ) cout << " No signal " << amps.size()  << endl;
//         continue;
//       }
// 
//       
//       if( amps.size() == 3 )
//       {
//         amp = (double)amps[2];
// //         if( debug ) cout << " Cathode " << id_ << " APV Ch0 " << amps[0] << " Ch1 " << amps[1]<< " Ch2 " << amps[2] << endl;
//       }
//       else
//       {
// //         if( debug ) cout << " APV Digit wrong sample size " << amps.size()  << endl;
//          cerr << " APV Digit wrong sample size " << amps.size()  << endl;
//       }  
//       double time = 0.;
//       
//       if( amp > 0 )
//       { 
// //         if( debug ) cout << " add new APV digit " << endl;
//         pads_.push_back( CsRICH1UpGrade::CathodePAD( id_, ixpad, iypad, amp, time) );
//         int data_id = SetPadAdr( ixpad,iypad );
//         data_apv_.push_back( ElectronicAPV::DataMCAPV( data_id, amps) );
//         pads_.back().el_digits_ = amps;
// //         if( debug ) cout << " add new APV digit OK " << endl;
//       }
//     }
//   }
// 
//   if( debug ) cout << " CathodeAPV::MakeElectronicMCDigits id " <<  id_ <<" OK data size " <<  pads_.size() << endl;
//   bool use_chisq_cut_apv_digits = true;
  bool use_chisq_cut_apv_digits = false;
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " CathodeAPV::MakeElectronicMCDigits id " <<  id_ << endl;
  
  for( unsigned ixpad = 0; ixpad < nxpads_; ixpad++)
  {  
    for( unsigned iypad = 0; iypad < nypads_; iypad++)
    {
      apv_.Clear();
      for( unsigned is=0; is< signals_[ixpad][iypad].size(); is++ )
      {
        apv_.AddSignal( signals_[ixpad][iypad][is].amp_, signals_[ixpad][iypad][is].time_ );
      }
      double amp = 0.;
      vector < int > amps = apv_.Read( ch_apv_[ixpad][iypad] );
      if( amps.size() == 0 )
      {
//         if( debug ) cout << " No signal " << amps.size()  << endl;
        continue;
      }

      
      if( amps.size() == 3 )
      {
        amp = (double)amps[2];
//         if( debug ) cout << " Cathode " << id_ << " APV Ch0 " << amps[0] << " Ch1 " << amps[1]<< " Ch2 " << amps[2] << endl;
      }
      else
      {
//         if( debug ) cout << " APV Digit wrong sample size " << amps.size()  << endl;
         cerr << " APV Digit wrong sample size " << amps.size()  << endl;
      }  
      double time = 0.;
      static vector< double > st;
      if( st.size() == 0 )
      {
        st = apv_.ReadTestSignal();
        st[0]=st[0]/st[2];
        st[1]=st[1]/st[2];
        st[2]=1.;
      }
      
      if( amp > 0 )
      {
//  APV MC chi^2  development       
        double sumtsq =  st[0]*st[0]+st[1]*st[1]+st[2]*st[2];
        double sumata =  st[0]*amps[0] + st[1]*amps[1] + st[2]*amps[2];
        double afit = sumata/sumtsq;
        double da0 = (amps[0] - afit*st[0]);
        double da1 = (amps[1] - afit*st[1]);
        double da2 = (amps[2] - afit*st[2]);
        double chisq = da0*da0 + da1*da1  + da2*da2;
// ???      chisq = (chisq/2.2)/(3.-1.);
        double disp = ch_apv_[ixpad][iypad].noise_/ch_apv_[ixpad][iypad].amp_unit_;
        disp *= disp;
        chisq = chisq/disp/2.;
        bool chisq_ok=false;
        if( chisq < 5. && afit > 10.*chisq && afit > 10. ) chisq_ok=true;
        bool digit_ok = true;
        if( use_chisq_cut_apv_digits ) digit_ok = chisq_ok;
       
//         if( debug ) cout << " add new APV digit " << endl;
        if( digit_ok )
        {
          pads_.push_back( CsRICH1UpGrade::CathodePAD( id_, ixpad, iypad, amp, time) );
          int data_id = SetPadAdr( ixpad,iypad );
          data_apv_.push_back( ElectronicAPV::DataMCAPV( data_id, amps) );
          pads_.back().el_digits_ = amps;
//           if( debug ) cout << " add new APV digit OK " << endl;
        }
      }
      else
      {
         cout << " APV a2=" << amp << endl;  
      }
    }
  }

  if( debug ) cout << " CathodeAPV::MakeElectronicMCDigits id " <<  id_ <<" OK data size " <<  pads_.size() << endl;
  
  
  
}

//////////////////////////////////////////////////////////////////////////////////

double  CathodeAPV::CathodeQuantumEfficiency(  double e, double x, double y )
{  
  return CathodeQuantumEfficiencyTable( e )/peak_eff_;
}

//////////////////////////////////////////////////////////////////////////////////

// c*********** find proper energy bin (jen) 
// 
// 	 enrg=5.5
// 	 rj=0.
// 	 jen=1
// 
//  20      if(eph.ge.(enrg+rj/10.) .and.eph.lt.(enrg+(rj+1.)/10.) ) then
//           goto 25
//          else
// 
//           if(jen.lt.24) then
//            jen=jen+1
//            rj=rj+1.
//            goto 20
//           endif
// 
// c!!! 25      endif
//          endif
//  25      continue
// 
// 
//          enrg=enrg+rj/10.
// 
// c********** CsI QE test
//  
//           qe=csiqe(jen)+(csiqe(jen+1)-csiqe(jen))*(eph-enrg)/0.1
//           xqe = qe*ridfact(icath)
//---   CsI quantum efficiency ( Eph 5.5-7.9 eV )   (def.=0.)
    const static double CsI_Quantum_Efficiency_Table[] = {
    0.00025,  0.00035,  0.0005,  0.005,  0.0075,
    0.0125,   0.01775,  0.0235,  0.04,   0.064,  
    0.094,    0.116,    0.1585,  0.192,  0.207,  
    0.23,     0.2525,   0.27,    0.2845, 0.3035,
    0.32,     0.34,     0.345,   0.36,   0.37            };

double  CathodeAPV::CathodeQuantumEfficiencyTable(  double e )
{
//   return 1.;
  if( e < 5.5 ) return 0.;
  int ie = (int)((e - 5.5)/0.1);
  if( ie >= 24 ) return 1.; //?????
  double elow = 5.5 + 0.1*double(ie);
  double ediff = (e-elow)/0.1;
  double qe = CsI_Quantum_Efficiency_Table[ie] +
                (CsI_Quantum_Efficiency_Table[ie+1] - CsI_Quantum_Efficiency_Table[ie])*ediff;
  return qe;              
}

//////////////////////////////////////////////////////////////////////////////////

CathodeAPV::Avalanche::Avalanche(  double x, double y, double e, double time, double multiplication )
{
  x_ = x;
  y_ = y;
  time_ = time;
  if( multiplication > 0. )
  { 
    amp_ = e*CsRICH1UpGrade::RanExp(multiplication);
  }
  else
  {
    amp_ = 0.;
  }  
}

//////////////////////////////////////////////////////////////////////////////////

void ElectronicAPV::ChannelData::Set(double amp_unit, double noise, double signal_cut, double  ped )
{
  amp_unit_ = fabs( amp_unit );
  noise_ = fabs( noise );
  signal_cut_in_sigmas_ = fabs( signal_cut );
  pedestal_ = ped;  
  signal_cut_in_units_ = (int)(signal_cut_in_sigmas_*noise_/amp_unit);
  int_pedestal_ = (int)(pedestal_/amp_unit_);
  if( pedestal_ - int_pedestal_ >= 0.5 ) int_pedestal_++;

}

//////////////////////////////////////////////////////////////////////////////////

ElectronicAPV::Shaper::Shaper( const std::vector < double > &shape, double ballistic_deficit,
                                                         double time_tik,unsigned zero_time_ch) :
  shape_( shape ),
  ballistic_deficit_( ballistic_deficit ),
  time_tik_( time_tik ),
  zero_time_ch_( zero_time_ch )                                       
{
  double max = -10e+10;
  int max_ch = -1;
  for( unsigned i=0; i < shape_.size(); i++ )
  {
    if( shape_[i] > max )
    { 
      max = shape_[i];
      max_ch = i;
    }
  }

  assert( max>0.);

  double factor = ballistic_deficit_/max;
  for( unsigned i=0; i < shape_.size(); i++ )
  {
    shape_[i] *= factor;
  }

}

//////////////////////////////////////////////////////////////////////////////////

double ElectronicAPV::Shaper::GetSignal( double amp, double time, double time_measuring ) const
{
  double dt = time_measuring - time;
  if( dt < 0 ) return 0.;
  dt += time_tik_/2.;
  int ch = zero_time_ch_ + int(dt/time_tik_);
  if( ch < 0 || ch >= (int)shape_.size() ) return 0.;
  return amp*shape_[ch];
}

//////////////////////////////////////////////////////////////////////////////////

ElectronicAPV::ElectronicAPV( void )
{
  add_noise_ = true;
  calculate_signal_only_in_measuring_poins_ = true;

  noise_ = 680.;
  pedestal_ = 5000.;
  noise_scale_factor_ = 1.;

  amp_unit_ = 320.;    
  time_tik_ = 25.;
  zero_time_ch_ = 999;

  signal_cut_in_sigmas_=3.;
  
  read_ch0_ = (int)zero_time_ch_ - 1;
  read_ch1_ = (int)zero_time_ch_ + 9;
  read_ch2_ = (int)zero_time_ch_ + 21;
  
  if( SIZE_ != 2000 )
  {
    cerr << " ERROR Initialisation ElectronicAPV: wrong size!! " << SIZE_ << endl;
    exit(1);
  }   

  vector < double > shape;
  for( unsigned i=0; i< 250; i++ )
  {
    double t = i*time_tik_;
    shape.push_back( ElectronicAPV::ShaperSignal0(t) );
  }
  shaper_ = new Shaper( shape, ballistic_deficit_, time_tik_, 21 );

  InitCsOpt();
  Clear();
}

//////////////////////////////////////////////////////////////////////////////////

void ElectronicAPV::InitCsOpt( void )
{
   CsOpt* opt = CsOpt::Instance();
   float f = amp_unit_;
   opt->CsOpt::getOpt( "MCSettingsAPV", "UnitScale", f );
   amp_unit_=f;

   f=noise_;
   opt->CsOpt::getOpt( "MCSettingsAPV", "PedestalsNoise", f );
   noise_=f;

   f=noise_scale_factor_;
   opt->CsOpt::getOpt( "MCSettingsAPV", "PedestalsNoiseScaleFactor", f );
   noise_scale_factor_=f;

   f = signal_cut_in_sigmas_;
   opt->CsOpt::getOpt( "MCSettingsAPV", "SignalCut", f );
   signal_cut_in_sigmas_ = f;

   std::string value = "";
   if( opt->CsOpt::getOpt( "MCSettingsAPV", "AddNoise", value ) )
   {
     if( value == "No" || value == "NO" || value == "N"  || value == "no"  || value == "n" ) add_noise_ = false;
   }
   
   if( opt->CsOpt::getOpt( "MCSettingsAPV", "CalculateInMeasuringPoints", value ) )
   {
     if( value == "No" || value == "NO" || value == "N"  || value == "no"  || value == "n" ) 
        calculate_signal_only_in_measuring_poins_ = false;
   }
}

//////////////////////////////////////////////////////////////////////////////////

void ElectronicAPV::PrintSettings( void )
{
   cout << " APV settings print out " << endl;
   cout << " APV amp unit(e) = " << amp_unit_ << endl;
   cout << " APV ped noise(e) = " << noise_ << endl;
   cout << " APV signal cut in ped sigma noise = " << signal_cut_in_sigmas_ << endl;

}

//////////////////////////////////////////////////////////////////////////////////

double ElectronicAPV::ShaperSignal0( double tns )
{
  double t = tns/25.;
  double tau1 = 20.;
  double tau2 = 19.9999;
  double tau3 = 19.;
  double tau4 = 73.;
  double shift = 20.;
//  double shift = 0.;
  if( t-shift < 0. ) return 0.;
  double amp = (tau1/(tau1-tau2))*(exp(-(t-shift)/tau1)-exp(-(t-shift)/tau2));
  return amp;
}

//////////////////////////////////////////////////////////////////////////////////

void ElectronicAPV::TestSignal( void )
{
// Date: Tue, 13 Sep 2005 10:23:00 +0200
// From: Igor Konorov <igor.konorov@cern.ch>
// To: Vladimir Kolosov <Vladimir.Kolosov@cern.ch>
// Subject: rc-cr shaper
// 
// 1,45*(1-EXP(-B39/tau3))*EXP(-B39/tau4)
// (tau1/(tau1-tau2))*(EXP(-(B43-shift)/tau1)-EXP(-(B43-shift)/tau2))
// 
// tau1 = 20
// tau2 = 19.9999
// tau3 = 19
// tau4 = 73

  double tau1 = 20.;
  double tau2 = 19.9999;
  double tau3 = 19.;
  double tau4 = 73.;
  double shift = 20.;

  for( int it=0; it<220; it++ )
  {
    double t=0.5+(double)it;
    
    double amp=0.;
    if( t-shift > 0 ) amp = (tau1/(tau1-tau2))*(exp(-(t-shift)/tau1)-exp(-(t-shift)/tau2));
    double amp1 = ElectronicAPV::ShaperSignal0( (t-shift)*25. );
    cout << " t= " << t*25. << " amp=" << amp <<  " amp1=" << amp1 <<endl;
  } 
}

//////////////////////////////////////////////////////////////////////////////////

void ElectronicAPV::Clear( void )
{
  for( unsigned i=0; i < SIZE_; i++ )
  {
    data_[i]=0.;
  }
}

//////////////////////////////////////////////////////////////////////////////////

// void ElectronicAPV::AddSignal( double amp, double time )
// {
//   if( calculate_signal_only_in_measuring_poins_ )
//   {
//     double time_measure = ( read_ch2_- (int)zero_time_ch_ )*time_tik_;
//     data_[read_ch2_] += shaper_->GetSignal( amp, time, time_measure );
//     time_measure = ( read_ch1_- (int)zero_time_ch_ )*time_tik_;
//     data_[read_ch1_] += shaper_->GetSignal( amp, time, time_measure );
//     time_measure = ( read_ch0_- (int)zero_time_ch_ )*time_tik_;
//     data_[read_ch0_] += shaper_->GetSignal( amp, time, time_measure );
//   }
//   else
//   {
//     for( unsigned i=0; i < (int)SIZE_; i++ )
//     {
//        double time_measure = ( i - (int)zero_time_ch_ )*time_tik_;
//        data_[i] += shaper_->GetSignal( amp, time, time_measure );
//     }
//   }
// }
// //////////////////////////////////////////////////////////////////////////////////

void ElectronicAPV::AddSignal( double amp, double time )
{
  if( calculate_signal_only_in_measuring_poins_ )
  {
    double time_measure2 = ( (int)read_ch2_- (int)zero_time_ch_ )*time_tik_;
    double a2 = shaper_->GetSignal( amp, time, time_measure2 );
    data_[read_ch2_] += a2;
    double time_measure1 = ( (int)read_ch1_- (int)zero_time_ch_ )*time_tik_;
    double a1 = shaper_->GetSignal( amp, time, time_measure1 );
    data_[read_ch1_] += a1;
    double time_measure0 = ( (int)read_ch0_- (int)zero_time_ch_ )*time_tik_;
    double a0 = shaper_->GetSignal( amp, time, time_measure0 );
    data_[read_ch0_] += a0;
//     cout << " amp " << amp << " time " << time << " a2 " << a2 << " a1 " << a1 << " a0 " << a0 << endl;
  }
  else
  {
    for( unsigned i=0; i < (int)SIZE_; i++ )
    {
       double time_measure = ( i - (int)zero_time_ch_ )*time_tik_;
       data_[i] += shaper_->GetSignal( amp, time, time_measure );
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////

vector<double> ElectronicAPV::ReadTestSignal( void )
{
  vector <double> out;
  double time = 0.;
  double amp = 1.;
  double time_measure0 = ( (int)read_ch0_- (int)zero_time_ch_ )*time_tik_;
  double a0 = shaper_->GetSignal( amp, time, time_measure0 );
  out.push_back(a0);
  double time_measure1 = ( (int)read_ch1_- (int)zero_time_ch_ )*time_tik_;
  double a1 = shaper_->GetSignal( amp, time, time_measure1 );
  out.push_back(a1);
  double time_measure2 = ( (int)read_ch2_- (int)zero_time_ch_ )*time_tik_;
  double a2 = shaper_->GetSignal( amp, time, time_measure2 );
  out.push_back(a2);
  return out;
}

//////////////////////////////////////////////////////////////////////////////////

vector<int> ElectronicAPV::Read( const ChannelData &ch )
{
  bool debug = false;
//   bool debug = true;
  bool zero_supression = true;
  double s2 = data_[read_ch2_] + ch.pedestal_;
  double s1 = data_[read_ch1_] + ch.pedestal_;
  double s0 = data_[read_ch0_] + ch.pedestal_;

// Add noise not in proper place just to speed up the process
  if( add_noise_ )
  {
    pair < double, double> rg = CsRICH1UpGrade::RanGau();
    double sig = ch.noise_*rg.first;
    s2 += sig*noise_scale_factor_;
    rg = CsRICH1UpGrade::RanGau();
    sig = ch.noise_*rg.first;
    s1 += sig*noise_scale_factor_;
    rg = CsRICH1UpGrade::RanGau();
    sig = ch.noise_*rg.first;
    s0 += sig*noise_scale_factor_;
  }
  double signal2_in_units = s2/ch.amp_unit_;
  double signal1_in_units = s1/ch.amp_unit_;
  double signal0_in_units = s0/ch.amp_unit_;

  int signal2 = (int)(signal2_in_units);
  int signal1 = (int)(signal1_in_units);
  int signal0 = (int)(signal0_in_units);

  if( debug ) cout << " signal_in_units " << signal2_in_units << " signal2 " <<  signal2 << endl;
  vector <int > output;
// Subtract pedestals

  signal2 -= ch.int_pedestal_;
  signal1 -= ch.int_pedestal_;
  signal0 -= ch.int_pedestal_;
  if( debug ) cout << "  signal2 - pedestal " <<  signal2 << endl;

  output.push_back( signal0 );
  output.push_back( signal1 );
  output.push_back( signal2 );

  if( zero_supression )
  {
    if( debug ) cout << " ch.signal_cut_in_units_  " << ch.signal_cut_in_units_  << endl;
    if( signal2 < ch.signal_cut_in_units_ ) output.clear();
  }
// remove last bit
  if( output.size() == 3 )
  {
//     output[0] = 2*(output[0]/2);
//     output[1] = 2*(output[1]/2);
//     output[2] = 2*(output[2]/2);
    ElectronicAPV::Encode(output[0],output[1],output[2],output);
//     cout << " s2 " << signal2 << " s1 " <<  signal1 << " s0 " <<  signal0 << endl;
//     cout << " data2 " << data_[read_ch2_] << " data1 " << data_[read_ch1_] << " data0 " << data_[read_ch0_] << endl;
  }

  return output;
  
}

//////////////////////////////////////////////////////////////////////////////////

void ElectronicAPV::Encode( int a0, int a1, int a2, vector <int> &out )
{
  bool new_encoding = true;
//   if( new_encoding )
//   {
//     if(a0 > 511 )
//     {
//       int a = (a0 - 512)/8+512;
//       if( a > 1022 ) a = 1022;
//       out[0] = a;      
//     }
//     else if( a0 <= 0)
//     {
//       out[0] = 0;      
//     }
//     else
//     {
//       out[0] = a0;      
//     }
// 
//     if(a1 > 511 )
//     {
//       int a = (a1 - 512)/8+512;
//       if( a > 1023 ) a = 1023;
//       out[1] = a;      
//     }
//     else if( a1 <= 0)
//     {
//       out[1] = 0;      
//     }
//     else
//     {
//       out[1] = a1;      
//     }
// 
//     if(a2 > 1023 )
//     {
//       int a = (a2 - 1024)/16+1024;
//       if( a > 2047 ) a = 2047;
//       out[2] = a;      
//     }
//     else if( a2 <= 0)
//     {
//       out[2] = 1;      
//     }
//     else
//     {
//       out[2] = a2;      
//     }
//   }
  if( new_encoding )
  {
    if(a0 > 63 )
    {
      int a = (a0 - 64)/8+64;
      if( a > 126 ) a = 126;
      out[0] = a;      
    }
    else if( a0 <= 0)
    {
      out[0] = 0;      
    }
    else
    {
      out[0] = a0;      
    }

    if(a1 > 255 )
    {
      int a = (a1 - 256)/8+256;
      if( a > 511 ) a = 511;
      out[1] = a;      
    }
    else if( a1 <= 0)
    {
      out[1] = 0;      
    }
    else
    {
      out[1] = a1;      
    }

    if(a2 > 255 )
    {
      int a = (a2 - 256)/8+256;
      if( a > 511 ) a = 511;
      out[2] = a;      
    }
    else if( a2 <= 0)
    {
      out[2] = 1;      
    }
    else
    {
      out[2] = a2;      
    }
  }
  else
  {
    out[0] = 2*(a0/2);
    out[1] = 2*(a1/2);
    out[2] = 2*(a2/2);
  }
}


// //////////////////////////////////////////////////////////////////////////////////
// 
// vector<int> ElectronicAPV::Read( const ChannelData &ch )
// {
//   bool debug = false;
// //   bool debug = true;
//   bool zero_supression = true;
//   double s2 = data_[read_ch2_] + ch.pedestal_;
//   double s1 = data_[read_ch1_] + ch.pedestal_;
//   double s0 = data_[read_ch0_] + ch.pedestal_;
// 
// // Add noise not in proper place just to seed up the process
//   if( add_noise_ )
//   {
//     pair < double, double> rg = CsRICH1UpGrade::RanGau();
//     double sig = ch.noise_*rg.first;
//     s2 += sig;
//     rg = CsRICH1UpGrade::RanGau();
//     sig = ch.noise_*rg.first;
//     s1 += sig;
//     rg = CsRICH1UpGrade::RanGau();
//     sig = ch.noise_*rg.first;
//     s0 += sig;
//   }
//   double signal2_in_units = s2/ch.amp_unit_;
//   double signal1_in_units = s1/ch.amp_unit_;
//   double signal0_in_units = s0/ch.amp_unit_;
// 
//   int signal2 = (int)(signal2_in_units);
//   int signal1 = (int)(signal1_in_units);
//   int signal0 = (int)(signal0_in_units);
// 
//   if( debug ) cout << " signal_in_units " << signal2_in_units << " signal2 " <<  signal2 << endl;
//   vector <int > output;
// // Subtract pedestals
// 
//   signal2 -= ch.int_pedestal_;
//   signal1 -= ch.int_pedestal_;
//   signal0 -= ch.int_pedestal_;
//   if( debug ) cout << "  signal2 - pedestal " <<  signal2 << endl;
// 
//   output.push_back( signal0 );
//   output.push_back( signal1 );
//   output.push_back( signal2 );
// 
//   if( zero_supression )
//   {
//     if( debug ) cout << " ch.signal_cut_in_units_  " << ch.signal_cut_in_units_  << endl;
//     if( signal2 < ch.signal_cut_in_units_ ) output.clear();
//   }
// // remove last bit
//   if( output.size() == 3 )
//   {
//     output[0] = 2*(output[0]/2);
//     output[1] = 2*(output[1]/2);
//     output[2] = 2*(output[2]/2);
//   }
// 
//   return output;
//   
// }
// 
//////////////////////////////////////////////////////////////////////////////////

const static double padtable[1600] = 
{
//   1  1
     2.07309E-07,
     1.22228E-05,
     1.22252E-05,
     2.07193E-07,
   0.,
     5.10901E-05,
     4.59625E-02,
     4.59607E-02,
     5.11046E-05,
   0.,
     1.48618E-04,
    0.453040,
    0.452964,
     1.48583E-04,
   0.,
     5.34265E-06,
     8.18364E-04,
     8.18553E-04,
     5.34310E-06,
   0.,
   0.,
     1.25631E-07,
     1.25186E-07,
   0.,
   0.,
//   1  2
     1.69018E-07,
     1.14939E-05,
     1.29448E-05,
     2.52946E-07,
   0.,
     3.97201E-05,
     4.08603E-02,
     5.10754E-02,
     6.56747E-05,
   0.,
     1.13704E-04,
    0.382067,
    0.523911,
     1.94183E-04,
     9.82488E-12,
     4.30081E-06,
     7.56939E-04,
     8.80099E-04,
     6.62448E-06,
   0.,
   0.,
     1.17561E-07,
     1.33426E-07,
   0.,
   0.,
//   1  3
     1.36688E-07,
     1.07648E-05,
     1.36517E-05,
     3.06743E-07,
   0.,
     3.08302E-05,
     3.58627E-02,
     5.60329E-02,
     8.43394E-05,
     3.17406E-10,
     8.69975E-05,
    0.315252,
    0.590723,
     2.53522E-04,
     5.76280E-09,
     3.45207E-06,
     6.96059E-04,
     9.40389E-04,
     8.18863E-06,
   0.,
   0.,
     1.09295E-07,
     1.41492E-07,
   0.,
   0.,
//   1  4
     1.09416E-07,
     1.00366E-05,
     1.43432E-05,
     3.69395E-07,
   0.,
     2.39372E-05,
     3.11280E-02,
     6.07500E-02,
     1.08191E-04,
     2.71822E-09,
     6.65154E-05,
    0.255740,
    0.650180,
     3.31070E-04,
     1.70820E-08,
     2.76202E-06,
     6.35643E-04,
     9.98822E-04,
     1.00998E-05,
   0.,
   0.,
     1.01695E-07,
     1.49571E-07,
   0.,
   0.,
//   1  5
     8.66259E-08,
     9.31725E-06,
     1.50112E-05,
     4.42718E-07,
   0.,
     1.85492E-05,
     2.67121E-02,
     6.51379E-02,
     1.38650E-04,
     7.83387E-09,
     5.08532E-05,
    0.204725,
    0.701111,
     4.32138E-04,
     3.34434E-08,
     2.20603E-06,
     5.77384E-04,
     1.05568E-03,
     1.24212E-05,
   0.,
   0.,
     9.36791E-08,
     1.57187E-07,
     8.23665E-11,
   0.,
//   1  6
     6.81823E-08,
     8.62114E-06,
     1.56503E-05,
     5.28109E-07,
   0.,
     1.43703E-05,
     2.26982E-02,
     6.91295E-02,
     1.77576E-04,
     1.53828E-08,
     3.88684E-05,
    0.162364,
    0.743341,
     5.63781E-04,
     5.65586E-08,
     1.75719E-06,
     5.21307E-04,
     1.10966E-03,
     1.52314E-05,
   0.,
   0.,
     8.62608E-08,
     1.64308E-07,
     3.34795E-10,
   0.,
//   1  7
     5.25577E-08,
     7.93738E-06,
     1.62519E-05,
     6.27453E-07,
   0.,
     1.11229E-05,
     1.91099E-02,
     7.26735E-02,
     2.27054E-04,
     2.68759E-08,
     2.97062E-05,
    0.127989,
    0.777551,
     7.35403E-04,
     8.77827E-08,
     1.39573E-06,
     4.67794E-04,
     1.15973E-03,
     1.86378E-05,
     2.87185E-11,
   0.,
     7.86226E-08,
     1.71504E-07,
     8.21121E-10,
   0.,
//   1  8
     4.00851E-08,
     7.27411E-06,
     1.68091E-05,
     7.41382E-07,
   0.,
     8.60103E-06,
     1.59555E-02,
     7.57554E-02,
     2.89919E-04,
     4.29727E-08,
     2.26833E-05,
     1.00364E-01,
    0.804972,
     9.58865E-04,
     1.29333E-07,
     1.10455E-06,
     4.17253E-04,
     1.20676E-03,
     2.27197E-05,
     7.58727E-10,
   0.,
     7.15072E-08,
     1.77776E-07,
     1.61267E-09,
   0.,
//   1  9
     2.99971E-08,
     6.64455E-06,
     1.73164E-05,
     8.71501E-07,
   0.,
     6.64603E-06,
     1.32199E-02,
     7.84268E-02,
     3.69824E-04,
     6.47326E-08,
     1.73217E-05,
     7.84181E-02,
    0.826619,
     1.24957E-03,
     1.83774E-07,
     8.71599E-07,
     3.69807E-04,
     1.24919E-03,
     2.76316E-05,
     2.70778E-09,
   0.,
     6.47314E-08,
     1.83525E-07,
     2.54000E-09,
   0.,
//   1  10
     2.16724E-08,
     6.03994E-06,
     1.77797E-05,
     1.02068E-06,
   0.,
     5.12681E-06,
     1.08840E-02,
     8.06783E-02,
     4.71068E-04,
     9.37763E-08,
     1.32188E-05,
     6.11011E-02,
    0.843546,
     1.62763E-03,
     2.55167E-07,
     6.84868E-07,
     3.26056E-04,
     1.28741E-03,
     3.34770E-05,
     5.91333E-09,
   0.,
     5.81860E-08,
     1.88980E-07,
     3.89537E-09,
   0.,
//   1  11
     1.54700E-08,
     5.46807E-06,
     1.81960E-05,
     1.19093E-06,
   0.,
     3.95290E-06,
     8.90607E-03,
     8.25186E-02,
     5.99139E-04,
     1.32538E-07,
     1.00887E-05,
     4.75120E-02,
    0.856657,
     2.11921E-03,
     3.48962E-07,
     5.35910E-07,
     2.85937E-04,
     1.32120E-03,
     4.04550E-05,
     1.05488E-08,
   0.,
     5.16721E-08,
     1.93599E-07,
     5.45619E-09,
   0.,
//   1  12
     1.03075E-08,
     4.92675E-06,
     1.85525E-05,
     1.38332E-06,
   0.,
     3.04280E-06,
     7.24991E-03,
     8.39921E-02,
     7.61154E-04,
     1.83786E-07,
     7.69148E-06,
     3.68867E-02,
    0.866668,
     2.75783E-03,
     4.71818E-07,
     4.17235E-07,
     2.49228E-04,
     1.34911E-03,
     4.87200E-05,
     1.74564E-08,
   0.,
     4.59420E-08,
     1.97525E-07,
     7.28410E-09,
   0.,
//   1  13
     6.52299E-09,
     4.42156E-06,
     1.88348E-05,
     1.59976E-06,
   0.,
     2.33833E-06,
     5.87287E-03,
     8.51915E-02,
     9.63827E-04,
     2.51050E-07,
     5.85823E-06,
     2.85877E-02,
    0.874115,
     3.58654E-03,
     6.32525E-07,
     3.23319E-07,
     2.16304E-04,
     1.37326E-03,
     5.84369E-05,
     2.72422E-08,
   0.,
     4.02424E-08,
     2.00917E-07,
     9.54209E-09,
   0.,
//   1  14
     3.84566E-09,
     3.95034E-06,
     1.90594E-05,
     1.84311E-06,
   0.,
     1.79344E-06,
     4.73553E-03,
     8.60449E-02,
     1.21922E-03,
     3.38754E-07,
     4.46285E-06,
     2.21320E-02,
    0.879528,
     4.66110E-03,
     8.43947E-07,
     2.48512E-07,
     1.86557E-04,
     1.39060E-03,
     6.98384E-05,
     3.99342E-08,
   0.,
     3.52181E-08,
     2.03795E-07,
     1.20659E-08,
   0.,
//   1  15
     1.89687E-09,
     3.51692E-06,
     1.92324E-05,
     2.11536E-06,
   0.,
     1.37314E-06,
     3.80606E-03,
     8.66585E-02,
     1.53896E-03,
     4.53163E-07,
     3.39539E-06,
     1.71087E-02,
    0.883153,
     6.05579E-03,
     1.12130E-06,
     1.89498E-07,
     1.60275E-04,
     1.40387E-03,
     8.31225E-05,
     5.66077E-08,
   0.,
     3.00290E-08,
     2.05407E-07,
     1.49361E-08,
   0.,
//   1  16
     6.91815E-10,
     3.11670E-06,
     1.93271E-05,
     2.41712E-06,
   0.,
     1.04732E-06,
     3.04850E-03,
     8.70267E-02,
     1.93724E-03,
     6.02622E-07,
     2.58005E-06,
     1.32159E-02,
    0.885233,
     7.85923E-03,
     1.48310E-06,
     1.43424E-07,
     1.36882E-04,
     1.41240E-03,
     9.85977E-05,
     7.84080E-08,
   0.,
     2.57464E-08,
     2.06348E-07,
     1.82928E-08,
   0.,
//   1  17
     1.11113E-10,
     2.75149E-06,
     1.93691E-05,
     2.74998E-06,
     1.02798E-10,
     7.96596E-07,
     2.43479E-03,
     8.71642E-02,
     2.43489E-03,
     7.97027E-07,
     1.95980E-06,
     1.01999E-02,
    0.885889,
     1.01984E-02,
     1.95891E-06,
     1.07065E-07,
     1.16525E-04,
     1.41501E-03,
     1.16443E-04,
     1.07122E-07,
   0.,
     2.19579E-08,
     2.06652E-07,
     2.18480E-08,
   0.,
//   1  18
   0.,
     2.41667E-06,
     1.93249E-05,
     3.11660E-06,
     6.85273E-10,
     6.02715E-07,
     1.93783E-03,
     8.69985E-02,
     3.04780E-03,
     1.04802E-06,
     1.48300E-06,
     7.85541E-03,
    0.885273,
     1.32097E-02,
     2.58072E-06,
     7.84878E-08,
     9.86300E-05,
     1.41144E-03,
     1.36845E-04,
     1.43196E-07,
   0.,
     1.82088E-08,
     2.06859E-07,
     2.58126E-08,
   0.,
//   1  19
   0.,
     2.11543E-06,
     1.92292E-05,
     3.51683E-06,
     1.93609E-09,
     4.53585E-07,
     1.53837E-03,
     8.66852E-02,
     3.80581E-03,
     1.37284E-06,
     1.12106E-06,
     6.05437E-03,
    0.883133,
     1.71045E-02,
     3.39529E-06,
     5.67260E-08,
     8.31632E-05,
     1.40417E-03,
     1.60238E-04,
     1.89796E-07,
   0.,
     1.49125E-08,
     2.05258E-07,
     3.02564E-08,
   0.,
//   1  20
   0.,
     1.84304E-06,
     1.90668E-05,
     3.95155E-06,
     3.77989E-09,
     3.38589E-07,
     1.21948E-03,
     8.60701E-02,
     4.73760E-03,
     1.79358E-06,
     8.44156E-07,
     4.66137E-03,
    0.879498,
     2.21336E-02,
     4.46349E-06,
     3.99053E-08,
     6.98433E-05,
     1.39093E-03,
     1.86631E-04,
     2.48591E-07,
   0.,
     1.19790E-08,
     2.03689E-07,
     3.52159E-08,
   0.,
//   1  21
   0.,
     1.60015E-06,
     1.88337E-05,
     4.42109E-06,
     6.58926E-09,
     2.51019E-07,
     9.63996E-04,
     8.52088E-02,
     5.87190E-03,
     2.33753E-06,
     6.32881E-07,
     3.58580E-03,
    0.874107,
     2.85807E-02,
     5.85981E-06,
     2.70561E-08,
     5.84372E-05,
     1.37299E-03,
     2.16233E-04,
     3.23222E-07,
   0.,
     9.54209E-09,
     2.00899E-07,
     4.03658E-08,
   0.,
//   1  22
   0.,
     1.38328E-06,
     1.85432E-05,
     4.92810E-06,
     1.03718E-08,
     1.83809E-07,
     7.60809E-04,
     8.39986E-02,
     7.24935E-03,
     3.04329E-06,
     4.71779E-07,
     2.75808E-03,
    0.866659,
     3.68895E-02,
     7.69072E-06,
     1.76395E-08,
     4.87023E-05,
     1.34961E-03,
     2.49326E-04,
     4.17093E-07,
   0.,
     7.33346E-09,
     1.97419E-07,
     4.57316E-08,
   0.,
//   1  23
   0.,
     1.19053E-06,
     1.81928E-05,
     5.46757E-06,
     1.54175E-08,
     1.32674E-07,
     5.99054E-04,
     8.25087E-02,
     8.90515E-03,
     3.95257E-06,
     3.48774E-07,
     2.11883E-03,
    0.856654,
     4.75271E-02,
     1.00866E-05,
     1.05848E-08,
     4.04460E-05,
     1.32087E-03,
     2.85912E-04,
     5.36022E-07,
   0.,
     5.43997E-09,
     1.93646E-07,
     5.19842E-08,
   0.,
//   1  24
   0.,
     1.02092E-06,
     1.77870E-05,
     6.04013E-06,
     2.18007E-08,
     9.39385E-08,
     4.71104E-04,
     8.06454E-02,
     1.08810E-02,
     5.12787E-06,
     2.55001E-07,
     1.62753E-03,
    0.843590,
     6.10935E-02,
     1.32193E-05,
     5.77758E-09,
     3.34794E-05,
     1.28755E-03,
     3.26131E-04,
     6.85038E-07,
   0.,
     3.79013E-09,
     1.88903E-07,
     5.78221E-08,
   0.,
//   1  25
   0.,
     8.71167E-07,
     1.73247E-05,
     6.64252E-06,
     2.98903E-08,
     6.47264E-08,
     3.70003E-04,
     7.84279E-02,
     1.32166E-02,
     6.64527E-06,
     1.83835E-07,
     1.24921E-03,
    0.826614,
     7.84248E-02,
     1.73138E-05,
     2.51731E-09,
     2.76235E-05,
     1.24956E-03,
     3.70007E-04,
     8.71282E-07,
   0.,
     2.59402E-09,
     1.83495E-07,
     6.45586E-08,
   0.,
//   1  26
   0.,
     7.40510E-07,
     1.68025E-05,
     7.27605E-06,
     4.01264E-08,
     4.28884E-08,
     2.89871E-04,
     7.57692E-02,
     1.59596E-02,
     8.59729E-06,
     1.29255E-07,
     9.58419E-04,
    0.804998,
     1.00322E-01,
     2.26719E-05,
     7.48394E-10,
     2.27172E-05,
     1.20606E-03,
     4.17301E-04,
     1.10403E-06,
   0.,
     1.61201E-09,
     1.77766E-07,
     7.16540E-08,
   0.,
//   1  27
   0.,
     6.27147E-07,
     1.62523E-05,
     7.93557E-06,
     5.26787E-08,
     2.72378E-08,
     2.27031E-04,
     7.26624E-02,
     1.91038E-02,
     1.11232E-05,
     8.76299E-08,
     7.35772E-04,
    0.777596,
    0.127961,
     2.96945E-05,
     2.11627E-11,
     1.86316E-05,
     1.16031E-03,
     4.67823E-04,
     1.39593E-06,
   0.,
     8.41215E-10,
     1.71077E-07,
     7.85802E-08,
   0.,
//   1  28
   0.,
     5.28038E-07,
     1.56475E-05,
     8.62076E-06,
     6.82719E-08,
     1.54686E-08,
     1.77590E-04,
     6.91121E-02,
     2.27019E-02,
     1.43740E-05,
     5.65840E-08,
     5.64032E-04,
    0.743392,
    0.162326,
     3.88758E-05,
   0.,
     1.52304E-05,
     1.10979E-03,
     5.21057E-04,
     1.75824E-06,
   0.,
     3.58233E-10,
     1.64363E-07,
     8.61930E-08,
   0.,
//   1  29
   0.,
     4.42961E-07,
     1.50148E-05,
     9.32188E-06,
     8.68880E-08,
     7.78258E-09,
     1.38682E-04,
     6.51392E-02,
     2.67179E-02,
     1.85546E-05,
     3.34386E-08,
     4.32048E-04,
    0.701187,
    0.204642,
     5.08483E-05,
   0.,
     1.24189E-05,
     1.05597E-03,
     5.77255E-04,
     2.20714E-06,
   0.,
     8.76614E-11,
     1.56897E-07,
     9.38072E-08,
   0.,
//   1  30
   0.,
     3.69487E-07,
     1.43526E-05,
     1.00363E-05,
     1.09456E-07,
     2.76833E-09,
     1.08208E-04,
     6.07509E-02,
     3.11314E-02,
     2.39222E-05,
     1.70468E-08,
     3.31138E-04,
    0.650185,
    0.255730,
     6.65187E-05,
   0.,
     1.00978E-05,
     9.98903E-04,
     6.35840E-04,
     2.76362E-06,
   0.,
   0.,
     1.49411E-07,
     1.01609E-07,
   0.,
//   1  31
   0.,
     3.06548E-07,
     1.36598E-05,
     1.07642E-05,
     1.36430E-07,
     3.21941E-10,
     8.43882E-05,
     5.60423E-02,
     3.58700E-02,
     3.08396E-05,
     5.67250E-09,
     2.53586E-04,
    0.590657,
    0.315302,
     8.69749E-05,
   0.,
     8.18498E-06,
     9.40435E-04,
     6.95730E-04,
     3.45026E-06,
   0.,
   0.,
     1.41290E-07,
     1.09234E-07,
   0.,
//   1  32
   0.,
     2.53044E-07,
     1.29464E-05,
     1.14998E-05,
     1.69084E-07,
   0.,
     6.57113E-05,
     5.10613E-02,
     4.08388E-02,
     3.97117E-05,
     8.31438E-12,
     1.94198E-04,
    0.523915,
    0.382098,
     1.13746E-04,
   0.,
     6.62525E-06,
     8.80258E-04,
     7.57216E-04,
     4.30028E-06,
   0.,
   0.,
     1.33515E-07,
     1.17411E-07,
   0.,
//   2  1
   0.,
     1.25453E-07,
     1.25707E-07,
   0.,
   0.,
     5.34311E-06,
     8.18512E-04,
     8.18151E-04,
     5.34119E-06,
   0.,
     1.48654E-04,
    0.452985,
    0.453052,
     1.48582E-04,
   0.,
     5.10660E-05,
     4.59493E-02,
     4.59418E-02,
     5.10987E-05,
   0.,
     2.07489E-07,
     1.22195E-05,
     1.22227E-05,
     2.07496E-07,
   0.,
//   2  2
   0.,
     1.17280E-07,
     1.33407E-07,
   0.,
   0.,
     4.29997E-06,
     7.57301E-04,
     8.79973E-04,
     6.62243E-06,
   0.,
     1.13768E-04,
    0.382082,
    0.523915,
     1.94209E-04,
     7.55812E-12,
     3.97213E-05,
     4.08404E-02,
     5.10757E-02,
     6.56884E-05,
   0.,
     1.68986E-07,
     1.14980E-05,
     1.29463E-05,
     2.52989E-07,
   0.,
//   2  3
   0.,
     1.09401E-07,
     1.41631E-07,
   0.,
   0.,
     3.45122E-06,
     6.95724E-04,
     9.40436E-04,
     8.18980E-06,
   0.,
     8.70220E-05,
    0.315230,
    0.590705,
     2.53684E-04,
     5.82473E-09,
     3.08563E-05,
     3.58795E-02,
     5.60571E-02,
     8.43785E-05,
     3.12951E-10,
     1.36261E-07,
     1.07712E-05,
     1.36594E-05,
     3.06877E-07,
   0.,
//   2  4
   0.,
     1.01399E-07,
     1.49219E-07,
   0.,
   0.,
     2.76347E-06,
     6.35972E-04,
     9.99341E-04,
     1.00986E-05,
   0.,
     6.65340E-05,
    0.255783,
    0.650114,
     3.31217E-04,
     1.69156E-08,
     2.39296E-05,
     3.11369E-02,
     6.07623E-02,
     1.08244E-04,
     2.81516E-09,
     1.09351E-07,
     1.00408E-05,
     1.43503E-05,
     3.70043E-07,
   0.,
//   2  5
   0.,
     9.37147E-08,
     1.56814E-07,
     7.48212E-11,
   0.,
     2.20592E-06,
     5.77468E-04,
     1.05603E-03,
     1.24259E-05,
   0.,
     5.08607E-05,
    0.204720,
    0.701098,
     4.32221E-04,
     3.36664E-08,
     1.85547E-05,
     2.67109E-02,
     6.51580E-02,
     1.38676E-04,
     7.77724E-09,
     8.68622E-08,
     9.32296E-06,
     1.50094E-05,
     4.42672E-07,
   0.,
//   2  6
   0.,
     8.61632E-08,
     1.63765E-07,
     3.37363E-10,
   0.,
     1.75738E-06,
     5.21133E-04,
     1.10966E-03,
     1.52353E-05,
   0.,
     3.88467E-05,
    0.162403,
    0.743295,
     5.63487E-04,
     5.63897E-08,
     1.43606E-05,
     2.26956E-02,
     6.91384E-02,
     1.77483E-04,
     1.54080E-08,
     6.79778E-08,
     8.61863E-06,
     1.56449E-05,
     5.28305E-07,
   0.,
//   2  7
   0.,
     7.85203E-08,
     1.71438E-07,
     8.65727E-10,
   0.,
     1.39595E-06,
     4.67498E-04,
     1.16015E-03,
     1.86357E-05,
     2.64518E-11,
     2.96986E-05,
    0.127973,
    0.777579,
     7.35421E-04,
     8.75605E-08,
     1.11291E-05,
     1.91021E-02,
     7.26695E-02,
     2.26997E-04,
     2.71194E-08,
     5.27794E-08,
     7.93844E-06,
     1.62508E-05,
     6.26855E-07,
   0.,
//   2  8
   0.,
     7.12709E-08,
     1.77797E-07,
     1.56281E-09,
   0.,
     1.10439E-06,
     4.16979E-04,
     1.20683E-03,
     2.27366E-05,
     7.26239E-10,
     2.26817E-05,
     1.00372E-01,
    0.804948,
     9.58842E-04,
     1.29306E-07,
     8.59753E-06,
     1.59598E-02,
     7.57670E-02,
     2.89942E-04,
     4.29334E-08,
     4.02665E-08,
     7.27632E-06,
     1.68093E-05,
     7.40898E-07,
   0.,
//   2  9
   0.,
     6.45394E-08,
     1.83287E-07,
     2.62977E-09,
   0.,
     8.71721E-07,
     3.69998E-04,
     1.24952E-03,
     2.76392E-05,
     2.61276E-09,
     1.73206E-05,
     7.84187E-02,
    0.826632,
     1.24953E-03,
     1.83765E-07,
     6.64541E-06,
     1.32251E-02,
     7.84068E-02,
     3.69938E-04,
     6.46951E-08,
     3.01390E-08,
     6.64449E-06,
     1.73218E-05,
     8.71764E-07,
   0.,
//   2  10
   0.,
     5.79545E-08,
     1.89250E-07,
     3.88016E-09,
   0.,
     6.84940E-07,
     3.25997E-04,
     1.28781E-03,
     3.34872E-05,
     5.80324E-09,
     1.32153E-05,
     6.11157E-02,
    0.843568,
     1.62739E-03,
     2.55124E-07,
     5.12676E-06,
     1.08817E-02,
     8.06443E-02,
     4.70975E-04,
     9.38347E-08,
     2.19257E-08,
     6.03766E-06,
     1.77820E-05,
     1.02056E-06,
   0.,
//   2  11
   0.,
     5.18486E-08,
     1.93343E-07,
     5.54304E-09,
   0.,
     5.35929E-07,
     2.85886E-04,
     1.32077E-03,
     4.04514E-05,
     1.06149E-08,
     1.00847E-05,
     4.75101E-02,
    0.856646,
     2.12013E-03,
     3.48984E-07,
     3.95200E-06,
     8.90383E-03,
     8.25342E-02,
     5.99082E-04,
     1.32386E-07,
     1.54094E-08,
     5.46555E-06,
     1.81870E-05,
     1.19049E-06,
   0.,
//   2  12
   0.,
     4.60251E-08,
     1.97574E-07,
     7.30388E-09,
   0.,
     4.17239E-07,
     2.49314E-04,
     1.34913E-03,
     4.86918E-05,
     1.76552E-08,
     7.68934E-06,
     3.68761E-02,
    0.866677,
     2.75688E-03,
     4.71863E-07,
     3.04243E-06,
     7.24883E-03,
     8.39966E-02,
     7.60375E-04,
     1.83285E-07,
     1.03732E-08,
     4.92437E-06,
     1.85378E-05,
     1.38253E-06,
   0.,
//   2  13
   0.,
     4.02903E-08,
     2.00939E-07,
     9.59548E-09,
   0.,
     3.23234E-07,
     2.16135E-04,
     1.37274E-03,
     5.84321E-05,
     2.71806E-08,
     5.85931E-06,
     2.85911E-02,
    0.874143,
     3.58694E-03,
     6.32687E-07,
     2.33746E-06,
     5.87084E-03,
     8.51619E-02,
     9.63983E-04,
     2.50808E-07,
     6.51953E-09,
     4.42232E-06,
     1.88290E-05,
     1.60058E-06,
   0.,
//   2  14
   0.,
     3.50790E-08,
     2.03748E-07,
     1.20198E-08,
   0.,
     2.48695E-07,
     1.86569E-04,
     1.39106E-03,
     6.98497E-05,
     3.98633E-08,
     4.46380E-06,
     2.21269E-02,
    0.879522,
     4.66153E-03,
     8.43880E-07,
     1.79348E-06,
     4.73746E-03,
     8.60517E-02,
     1.21976E-03,
     3.38512E-07,
     3.77890E-09,
     3.95097E-06,
     1.90673E-05,
     1.84314E-06,
   0.,
//   2  15
   0.,
     3.05035E-08,
     2.05605E-07,
     1.48776E-08,
   0.,
     1.89764E-07,
     1.60229E-04,
     1.40395E-03,
     8.31671E-05,
     5.65697E-08,
     3.39701E-06,
     1.71141E-02,
    0.883144,
     6.05498E-03,
     1.12044E-06,
     1.37290E-06,
     3.80720E-03,
     8.66629E-02,
     1.53828E-03,
     4.53117E-07,
     1.88153E-09,
     3.51719E-06,
     1.92312E-05,
     2.11496E-06,
   0.,
//   2  16
   0.,
     2.57347E-08,
     2.06632E-07,
     1.82829E-08,
   0.,
     1.43147E-07,
     1.36928E-04,
     1.41191E-03,
     9.85834E-05,
     7.87859E-08,
     2.58114E-06,
     1.32166E-02,
    0.885248,
     7.85745E-03,
     1.48334E-06,
     1.04770E-06,
     3.04861E-03,
     8.70123E-02,
     1.93819E-03,
     6.02661E-07,
     6.88875E-10,
     3.11686E-06,
     1.93317E-05,
     2.41744E-06,
   0.,
//   2  17
   0.,
     2.18264E-08,
     2.07106E-07,
     2.19488E-08,
   0.,
     1.07024E-07,
     1.16467E-04,
     1.41440E-03,
     1.16478E-04,
     1.07067E-07,
     1.95837E-06,
     1.01985E-02,
    0.885936,
     1.01941E-02,
     1.95864E-06,
     7.96234E-07,
     2.43421E-03,
     8.71255E-02,
     2.43339E-03,
     7.96493E-07,
     1.13357E-10,
     2.75097E-06,
     1.93655E-05,
     2.75251E-06,
     9.97538E-11,
//   2  18
   0.,
     1.81939E-08,
     2.06700E-07,
     2.60097E-08,
   0.,
     7.85139E-08,
     9.86600E-05,
     1.41213E-03,
     1.36922E-04,
     1.43448E-07,
     1.48306E-06,
     7.86090E-03,
    0.885242,
     1.32120E-02,
     2.58097E-06,
     6.02971E-07,
     1.93797E-03,
     8.70199E-02,
     3.04906E-03,
     1.04759E-06,
   0.,
     2.41663E-06,
     1.93271E-05,
     3.11618E-06,
     6.58246E-10,
//   2  19
   0.,
     1.50434E-08,
     2.05326E-07,
     3.03305E-08,
   0.,
     5.68382E-08,
     8.31597E-05,
     1.40358E-03,
     1.60215E-04,
     1.89747E-07,
     1.12027E-06,
     6.05422E-03,
    0.883155,
     1.71027E-02,
     3.39444E-06,
     4.53440E-07,
     1.53858E-03,
     8.66640E-02,
     3.80658E-03,
     1.37257E-06,
   0.,
     2.11511E-06,
     1.92262E-05,
     3.51767E-06,
     1.88440E-09,
//   2  20
   0.,
     1.19963E-08,
     2.04065E-07,
     3.51804E-08,
   0.,
     3.97036E-08,
     6.98673E-05,
     1.39097E-03,
     1.86655E-04,
     2.48666E-07,
     8.44151E-07,
     4.66294E-03,
    0.879527,
     2.21267E-02,
     4.46245E-06,
     3.38955E-07,
     1.21952E-03,
     8.60451E-02,
     4.73799E-03,
     1.79413E-06,
   0.,
     1.84331E-06,
     1.90685E-05,
     3.95335E-06,
     3.74619E-09,
//   2  21
   0.,
     9.53336E-09,
     2.00685E-07,
     4.02460E-08,
   0.,
     2.70395E-08,
     5.84377E-05,
     1.37306E-03,
     2.16244E-04,
     3.23279E-07,
     6.33069E-07,
     3.58668E-03,
    0.874136,
     2.85973E-02,
     5.86083E-06,
     2.50615E-07,
     9.63903E-04,
     8.51609E-02,
     5.87268E-03,
     2.33912E-06,
   0.,
     1.59931E-06,
     1.88338E-05,
     4.42236E-06,
     6.56378E-09,
//   2  22
   0.,
     7.35659E-09,
     1.97833E-07,
     4.59386E-08,
   0.,
     1.77123E-08,
     4.86993E-05,
     1.34958E-03,
     2.49303E-04,
     4.17313E-07,
     4.71905E-07,
     2.75766E-03,
    0.866677,
     3.68875E-02,
     7.68906E-06,
     1.83849E-07,
     7.60767E-04,
     8.39844E-02,
     7.24782E-03,
     3.04137E-06,
   0.,
     1.38284E-06,
     1.85388E-05,
     4.92573E-06,
     1.02492E-08,
//   2  23
   0.,
     5.46121E-09,
     1.93461E-07,
     5.18285E-08,
   0.,
     1.06311E-08,
     4.04570E-05,
     1.32117E-03,
     2.86004E-04,
     5.36277E-07,
     3.48994E-07,
     2.11889E-03,
    0.856655,
     4.75224E-02,
     1.00873E-05,
     1.32818E-07,
     5.99236E-04,
     8.25123E-02,
     8.90428E-03,
     3.95358E-06,
   0.,
     1.19068E-06,
     1.81950E-05,
     5.46713E-06,
     1.53724E-08,
//   2  24
   0.,
     3.87935E-09,
     1.89014E-07,
     5.79048E-08,
   0.,
     5.82847E-09,
     3.34944E-05,
     1.28790E-03,
     3.26283E-04,
     6.85175E-07,
     2.55065E-07,
     1.62726E-03,
    0.843511,
     6.11102E-02,
     1.32176E-05,
     9.39722E-08,
     4.71046E-04,
     8.07051E-02,
     1.08830E-02,
     5.12771E-06,
   0.,
     1.02100E-06,
     1.77855E-05,
     6.04114E-06,
     2.16855E-08,
//   2  25
   0.,
     2.60900E-09,
     1.83423E-07,
     6.46448E-08,
   0.,
     2.59577E-09,
     2.76283E-05,
     1.24964E-03,
     3.69920E-04,
     8.71279E-07,
     1.83605E-07,
     1.24989E-03,
    0.826615,
     7.84063E-02,
     1.73227E-05,
     6.46656E-08,
     3.69860E-04,
     7.84361E-02,
     1.32250E-02,
     6.64354E-06,
   0.,
     8.71351E-07,
     1.73193E-05,
     6.64155E-06,
     2.99721E-08,
//   2  26
   0.,
     1.55326E-09,
     1.77564E-07,
     7.14699E-08,
   0.,
     7.40926E-10,
     2.27218E-05,
     1.20659E-03,
     4.17186E-04,
     1.10463E-06,
     1.29214E-07,
     9.58632E-04,
    0.804998,
     1.00357E-01,
     2.26815E-05,
     4.30005E-08,
     2.89930E-04,
     7.57487E-02,
     1.59437E-02,
     8.59879E-06,
   0.,
     7.40890E-07,
     1.68046E-05,
     7.27489E-06,
     4.02102E-08,
//   2  27
   0.,
     8.51741E-10,
     1.71288E-07,
     7.87609E-08,
   0.,
     2.79631E-11,
     1.86282E-05,
     1.16044E-03,
     4.67797E-04,
     1.39492E-06,
     8.74583E-08,
     7.35659E-04,
    0.777587,
    0.127973,
     2.96946E-05,
     2.69031E-08,
     2.27043E-04,
     7.26552E-02,
     1.91079E-02,
     1.11191E-05,
   0.,
     6.27287E-07,
     1.62485E-05,
     7.93824E-06,
     5.25773E-08,
//   2  28
   0.,
     3.75621E-10,
     1.64363E-07,
     8.60112E-08,
   0.,
   0.,
     1.52364E-05,
     1.10971E-03,
     5.21250E-04,
     1.75662E-06,
     5.63411E-08,
     5.63837E-04,
    0.743260,
    0.162449,
     3.88750E-05,
     1.55266E-08,
     1.77513E-04,
     6.91165E-02,
     2.27062E-02,
     1.43735E-05,
   0.,
     5.28021E-07,
     1.56485E-05,
     8.61861E-06,
     6.80741E-08,
//   2  29
   0.,
     6.95295E-11,
     1.57083E-07,
     9.36318E-08,
   0.,
   0.,
     1.24231E-05,
     1.05605E-03,
     5.77325E-04,
     2.20773E-06,
     3.36273E-08,
     4.32241E-04,
    0.701175,
    0.204639,
     5.08595E-05,
     7.77181E-09,
     1.38678E-04,
     6.51548E-02,
     2.67174E-02,
     1.85573E-05,
   0.,
     4.42845E-07,
     1.50132E-05,
     9.32344E-06,
     8.70011E-08,
//   2  30
   0.,
   0.,
     1.49722E-07,
     1.01553E-07,
   0.,
   0.,
     1.00896E-05,
     9.98878E-04,
     6.35958E-04,
     2.76314E-06,
     1.69832E-08,
     3.31018E-04,
    0.650211,
    0.255699,
     6.65118E-05,
     2.71127E-09,
     1.08206E-04,
     6.07527E-02,
     3.11348E-02,
     2.39345E-05,
   0.,
     3.69670E-07,
     1.43488E-05,
     1.00377E-05,
     1.09365E-07,
//   2  31
   0.,
   0.,
     1.41432E-07,
     1.09307E-07,
   0.,
   0.,
     8.18754E-06,
     9.40655E-04,
     6.95788E-04,
     3.45011E-06,
     5.61449E-09,
     2.53523E-04,
    0.590683,
    0.315311,
     8.69709E-05,
     3.24175E-10,
     8.43415E-05,
     5.60156E-02,
     3.58614E-02,
     3.08446E-05,
   0.,
     3.06340E-07,
     1.36590E-05,
     1.07659E-05,
     1.36255E-07,
//   2  32
   0.,
   0.,
     1.33380E-07,
     1.17617E-07,
   0.,
   0.,
     6.61937E-06,
     8.79899E-04,
     7.57072E-04,
     4.29904E-06,
     9.82329E-12,
     1.94122E-04,
    0.523902,
    0.382120,
     1.13697E-04,
   0.,
     6.56572E-05,
     5.10326E-02,
     4.08594E-02,
     3.97049E-05,
   0.,
     2.53036E-07,
     1.29468E-05,
     1.14962E-05,
     1.69125E-07
 };

//////////////////////////////////////////////////////////////////////////////////

double  CathodeAPV::PadResponse( double dxpad, double dypad, int ixs, int iys )
{
  int iy = (int)(16+dypad/0.25);
  if( iy < 0 || iy > 31 )
  {
    cerr << " PadResponse iy not perfect " << iy << " dypad " << dypad << endl;
    iy=0;
  }
  int ix = (int)(1+dxpad/8.);
  if( ix < 0 || ix > 1 )
  {
    cerr << " PadResponse ix not perfect " << ix << " dxpad " << dxpad << endl;
    ix =0;
  }
  if( ixs < -2 || ixs > 2 || iys < -2 || iys > 2 ) return 0.;

  int is = (ixs+2)*5+(iys+2);
  
  int indx_new = ix*800+iy*25+is;
  assert( indx_new < 1600 );
  assert( indx_new >= 0 );
  
  return padtable[indx_new];
}

//////////////////////////////////////////////////////////////////////////////////


