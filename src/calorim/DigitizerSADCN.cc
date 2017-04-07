// --- Standard C/C++ library ---
# include <sstream>
# include <cassert>

#include "DigitizerSADCN.h"
#include "CsCalorimeter.h"

using namespace std;
//using namespace MN;

#define CDB_LINE_MAX 132

//namespace MN {

const size_t ManagerShapeTableSADC::NTABS  =  25;
const double ManagerShapeTableSADC::HFWmin =   0.2;
const double ManagerShapeTableSADC::HFWmax =   0.45;
const double ManagerShapeTableSADC::Xmin   =  -1.5;
const double ManagerShapeTableSADC::Xmax   =   4.5;
const double ManagerShapeTableSADC::NXbins = 600;

////////////////////////////////////////////////////////////////////////////////

bool DigitizerSADCN::Check( void ) const
{
// Nothing to Init just check
  if( GetParent() == NULL   )
  {
    cerr <<" DigitizerSADCN::Check(): GetParent() == NULL which is fatal in present code implementation " << endl;
    exit(1);
  }
  const CsCalorimeter * c = GetParent();
  if( c == NULL )
  {
    cerr <<" DigitizerSADCN::Init(): GetParent()  is not Calorimeter which is fatal in present code implementation " << endl;
    exit(1);
  }

//  mtab_ = c->GetManShapeTableSADC();
  if( mtab_ == NULL )
  {
    cerr <<" DigitizerSADCN::Check():  GetManShapeTableSADC() == NULL which is fatal in present code implementation " << endl;
    exit(1);
  }
// TODO  mtab_->GetProfileMaps

  bool oktables = mtab_->CheckTables();
//  cout <<" oktables " << oktables << endl;
  if( !oktables )
  {
    cerr <<" DigitizerSADCN::Check(): ShapeTable Manager was not properly initialized from calibration file: No way to perform digitization !! " << endl;
//    exit(1);
    return oktables;
  }

  if( calib_fwhm_tab_ <= 1. )
  {
    cerr <<"  DigitizerSADCN::Check(): calib_fwhm_tab_ out of range " << calib_fwhm_tab_ << endl;
    return false;
  }

 if(  calib_hfwn_tab_ <= 0.01 )
  {
    cerr <<"  DigitizerSADCN::Check(): calib_hfwn_tab_ out of range " << calib_hfwn_tab_ << endl;
    return false;
  }

 return true;
}

////////////////////////////////////////////////////////////////////////////////

void DigitizerSADCN::Init( void )
{
  if( mtab_ != NULL ) return;
  if( GetParent() == NULL   )
  {
    cerr <<" DigitizerSADCN::Init(): GetParent() == NULL which is fatal in present code implementation " << endl;
    exit(1);
  }
  const CsCalorimeter * c = GetParent();
  if( c == NULL )
  {
    cerr <<" DigitizerSADCN::Init(): GetParent()  is not Calorimeter which is fatal in present code implementation " << endl;
    exit(1);
  }
//  mtab_ = c->GetManShapeTableSADC();
  if( mtab_ == NULL )
  {
    cerr <<" DigitizerSADCN::Init():  GetManShapeTableSADC() == NULL which is fatal in present code implementation " << endl;
    exit(1);
  }
// TODO  mtab_->GetProfileMaps

}

////////////////////////////////////////////////////////////////////////////////

double DigitizerSADCN::GetTabShapeValue( double time ) const
{
// Input parameter time is a time in ns relative to forward pulse front which is positioned at zero value for our tables
// Special scaling due to average tabels usage
  double x = time/calib_fwhm_tab_;
  return mtab_->Value( calib_hfwn_tab_, x );
}

////////////////////////////////////////////////////////////////////////////////

std::pair <double,double> DigitizerSADCN::GetTabShapeWValue( double wtime ) const
{
// Input parameter time is a time in ns relative to forward pulse front which is positioned at zero value for our tables
// Special scaling due to average tabels usage
  double x = wtime/calib_fwhm_tab_;
  std::pair <double,double> ret = mtab_->WValue( calib_hfwn_tab_, x );
  ret.second *= calib_fwhm_tab_;  // Need to rescale time position for plato begining
  return ret;
}

////////////////////////////////////////////////////////////////////////////////

bool DigitizerSADCN::FitBase(  const std::vector<uint16> &sample )
{
  Init(); // ?? Not clearwhatto Init ??? TODO toberemoved
  FitMaxAdvancedN(sample);
  return ResultIsReady();
}

////////////////////////////////////////////////////////////////////////////////

void DigitizerSADCN::SetDefaultCalibrations ( void )
{
// Harcoded table default parameters to start with. Parameters must be cell dependent.
  double calib_fwhm = 100.;
  double calib_fwhm_sigma = 10.;
  double calib_hfwn = 0.33;
  double calib_hfwn_sigma = 0.022;
  double parFWHMn_HFW_fwhmn_c = 0.;
  double parFWHMn_HFW_fwhmn_r = 3.;
//  double parFWHMn_HFW_hfw_c = 0.3;
  double parFWHMn_HFW_hfw_c = 0.;
  double parFWHMn_HFW_hfw_r = 0.15;
  if( GetParentName() == "HC02P1__" )
  {
    calib_fwhm = 112.5;
    calib_fwhm_sigma = 15.;
    calib_hfwn = 0.3076;
    calib_hfwn_sigma = 0.02717;
  }
  else if( GetParentName() == "HC01P1__" )
  {
    calib_fwhm = 104.5;
    calib_fwhm_sigma = 10.;
    calib_hfwn = 0.3024;
    calib_hfwn_sigma = 0.0166;
  }
  else if( GetParentName() == "EC01P1__" )
  {
    calib_fwhm = 95.0;
    calib_fwhm_sigma = 15.;
    calib_hfwn = 0.3141;
    calib_hfwn_sigma = 0.026;

    parFWHMn_HFW_fwhmn_c = -0.7;
    parFWHMn_HFW_fwhmn_r = 3.;
//    parFWHMn_HFW_hfw_c = 0.27;
    parFWHMn_HFW_hfw_c = 0.;
    parFWHMn_HFW_hfw_r = 0.15;
  }
  else if( GetParentName() == "EC02P1__" )
  {
    calib_fwhm = 87.5;
    calib_fwhm_sigma = 10.;
    calib_hfwn = 0.3241;
    calib_hfwn_sigma = 0.022;

    parFWHMn_HFW_fwhmn_c = -0.5;
    parFWHMn_HFW_fwhmn_r = 3.;
//    parFWHMn_HFW_hfw_c = 0.27;
    parFWHMn_HFW_hfw_c = 0.;
    parFWHMn_HFW_hfw_r = 0.15;
  }
  else if( GetParentName() == "EC00P1__" )
  {
    calib_fwhm = 79.5;
    calib_fwhm_sigma = 10.;
    calib_hfwn = 0.367;
    calib_hfwn_sigma = 0.01373;
  }

  calib_fwhm_tab_ = calib_fwhm;
  calib_fwhm_sigma_tab_ = calib_fwhm_sigma;
  calib_hfwn_tab_ = calib_hfwn;
  calib_hfwn_sigma_tab_ = calib_hfwn_sigma;
  parFWHMn_HFW_fwhmn_c_ = parFWHMn_HFW_fwhmn_c;
  parFWHMn_HFW_fwhmn_r_ = parFWHMn_HFW_fwhmn_r;
  parFWHMn_HFW_hfw_c_ = parFWHMn_HFW_hfw_c;
  parFWHMn_HFW_hfw_r_ = parFWHMn_HFW_hfw_r;
}

////////////////////////////////////////////////////////////////////////////////

std::vector <double>  fit_parabola( const std::vector<double> &arg, const std::vector<double> &value )
{
  bool debug = false;
  std::vector <double> result;
  result.push_back(0.);
  result.push_back(0.);
  result.push_back(0.);
  if( arg.size() <3 || arg.size() != value.size() ) return result; // TODO break execution

  double s0=0.;
  double s1=0.;
  double s2=0.;
  double s3=0.;
  double s4=0.;
  double sy0=0.;
  double sy1=0.;
  double sy2=0.;

  for( int is=0; is < (int)arg.size(); is++)
  {
    double x = arg[is];
    double xx = x*x;
    double xxx = xx*x;
    double xxxx = xx*xx;
    s0 += 1.;
    s1 += x;
    s2 += xx;
    s3 += xxx;
    s4 += xxxx;
    double y = value[is];
    sy0 += y;
    sy1 += y*x;
    sy2 += y*xx;
  }

  s1 /= s0;
  s2 /= s0;
  s3 /= s0;
  s4 /= s0;
  sy0 /= s0;
  sy1 /= s0;
  sy2 /= s0;

  double det = s4*s2 + 2.*s3*s2*s1 - s2*s2*s2 - s4*s1*s1 -s3*s3;
  double a = (( sy2-sy0*s2)*(s2-s1*s1) - ( sy1-sy0*s1)*(s3-s1*s2))/det;
  double b =  ( sy1 - a*s3 - sy0*s1 + a*s2*s1 )/(s2-s1*s1);
  double c =  sy0 - a*s2 - s1*b;
  double xc = -b/(2*a);
  double yc = c- a*xc*xc;
  if( debug ) cout  << " xc " << xc << endl;
  result[0] = xc;
  result[1] = sqrt( fabs(1./a));
  result[2] = yc;
  return result;
}

void DigitizerSADCN::FitMaxAdvancedN ( const std::vector<uint16> &sample )
{
  // Control debuging
  bool debug = debug_;
  double signal_cut2print4debug = 2.;

  if( debug )
    cout << " DigitizerSADC::FitMaxAdvanced debug " << endl;

  // Reset
  Clear();

  if( debug )
    cout << " DigitizerSADC::FitMaxAdvanced Clear() OK " << endl;
  if( debug )
    PrintSettings();

  // Do basline subtraction and base shape analysis
  GetCcSample( sample );
//  GetCSample( sample );


  if( debug ) cout << " DigitizerSADC::FitMaxAdvanced overflow_amplitude_overflow_amplitude_overflow_amplitude_ " << overflow_amplitude_ << endl;

// Overflow not detected corrected due to OLD shape references

  // Detect overflows
// bool debug_over = true;
  bool debug_over = false;
//  if( debug_over) cout << " OverDebug " << GetParentName() <<" overflow_amplitude_ " << overflow_amplitude_ << endl;

  vector <int> overflow_samples;
  unsigned int new_over(0);
  double time_over = 0.;
  double amp_over = 0.;

// Search for overflow and fill some primaty information

  {

    for( int is=0; is < (int)sample.size(); is++)
      {
        if( sample[is] >  overflow_amplitude_ )
          {
            if( sample[is] >=  overflow_amplitude_)
              {
                if( !new_over )
                  {
                    new_over++;
                    overflow_samples.push_back(is);
                  }
                else if( overflow_samples.back() < is-1 )
                  {
                    new_over++;
                    overflow_samples.push_back(is);
                  }
                else
                  {
                    overflow_samples.push_back(is);
                  }
              }
          }
      }
      if( debug_over && new_over )
      {
        cout << " OverDebug " << GetParentName() <<" new_over " << new_over << " nover " << overflow_samples.size() << endl;
        if( overflow_samples.size() >= 9 )  PrintSample( sample );
      }


    // Treat overflows
    if( new_over  )
    {
      if( new_over > 1 )
      {
        cerr << " Two samples? And both in overflow ?? Ne byvaet !.. " << endl;
        PrintSample( sample );
      }
      double over_plato_width = double(overflow_samples.size())*GetSADCClock();
      int start_plato = overflow_samples[0];
      int end_plato = overflow_samples.back();
      std::pair <double,double> tabhint = GetTabShapeWValue( over_plato_width );
      if( debug_over ) cout <<" Multiplicator " << tabhint.first <<" Plato starts at time " << tabhint.second << endl;
// Za vremya plato prinimaen NACHALO bina soverfolow !!!
      time_over = double(start_plato-0.5)*GetSADCClock();
      double amp_over_actual = csample_[start_plato];    // It must be equal: ( overflow_amplitude_-GetDynamicPed()  ) which will be nice to test
// As we now in test mode we'll replace   amp_over_actual  by:
// TODO check the approx equality for real overflows
      amp_over_actual  = overflow_amplitude_-GetDynamicPed();

      amp_over = amp_over_actual/tabhint.first;
      time_over -= tabhint.second;


      if( debug_over ) cout <<" amp_over " << amp_over  <<" time_over " << time_over << endl;
// Try development with Pulse Shape Overflows
 double timeerr_over = 0.;
 double chisq_over = 0.;
 int maxgate_over = 5;
 bool make_psh_fit_over = true;
 if( make_psh_fit_over )
 {  // begin neshapefit scope

     bool printsam_over = false;
//      bool printsam_over = true;
//      if( signal> 200. ) printsam = true;
//    if( GetParentName() == "EC02P1__" && signal> 1000. ) printsam = true;

//     if( printsam_over ) cout <<" parent " << GetParent()  << "  "<< GetParentName() << " fwhm " << calib_fwhm_tab_ <<  " hfwn " << calib_hfwn_tab_ <<endl;

//     if( GetParent() != NULL   )
//     {
//       const CalorimeterMN * c = dynamic_cast<const CalorimeterMN *> (GetParent());
//       if( c != NULL )
//       {

// This must be hidden
//      const ManagerShapeTableSADC *mtab = c->GetManShapeTableSADC();
//       const ManagerShapeTableSADC *mtab = mtab_;
//       if( mtab != NULL )
//       {
//   All exists

        if( printsam_over ) cout << GetParentName() <<" sample  amp_over " << amp_over <<" time_over " << time_over << endl;

        std::vector<double> pp;
        std::vector<double> py;
        std::vector<double> t;
        std::vector<double> ss;
        size_t nfitp = 6;
// We estimate time ucertainty in case of overflows as +/- 6 ns from width estimator. To be sure we apply for the fit +/-12ns this range can be squeezed later( now +/- 6ns)
        double stept = 1.0; // can be amplitude/energy dependent

        for( size_t np=0; np<nfitp; np++)
        {
          pp.push_back(0.);
          py.push_back(0.);
          ss.push_back(0.);
          t.push_back( stept*(-3.+double(np))  );
        }

        double timepar = time_over;
        double signalpar = amp_over;
        double yy =0.;
        int min_sample = start_plato - maxgate_over;
        int max_sample = end_plato + maxgate_over;
// This value(amp_over_norm)  gives an overflow level cut line on the shape curve. We relay on the amplitude estimation (signalpar) and hence this value is a matter of iterations.
//  There are more than one approach:
//  1) skip in the fit bins with overflow. It is OK as we can use amplitude estimator without change. Disadvatage is that if the shape value below overflow we will not add penalty for ChiSq.
//  2) Modify shape by overflow flat top, but then there is strong dependence from heght parameter.
//    variable (amp_over_norm) can be used in 2) approach ( see upper)
        double amp_over_norm = amp_over_actual/signalpar;
        if( debug ) cout <<" amp_over_norm " << amp_over_norm << endl;
        int add_bins2fit = 3;
        for( int is=min_sample; is <= max_sample; is++)
        {
          if( is <= 0 ) continue;
          if( is >= (int)sample.size() ) break;
          if( is < start_plato - add_bins2fit) continue;
          if( is < end_plato + add_bins2fit) continue;
          if( is >= start_plato && is <= end_plato ) continue; // exclude overflows from the fit
          if( sample[is] >=  overflow_amplitude_) continue;

          double timex = double(is)*sadc_clock_-timepar;
// vtab and dy here is just for printing. Must be commented!
          double vtab = GetTabShapeValue( timex );
          double dy = csample_[is]-signalpar*vtab;

          double y = csample_[is];
          double sigy = 1.+0.01*fabs(y);
          double w = 1./sigy/sigy;
          yy += y*y*w;
          for( size_t np=0; np<nfitp; np++)
          {
             double p = GetTabShapeValue( timex + t[np] );
             pp[np] += p*p*w;
             py[np] += p*y*w;
          }
          if( printsam_over ) cout <<" timex " << timex  << " V " <<   csample_[is] <<"( " << signalpar*vtab <<") "<<" dy " << dy;
        }

        if( printsam_over ) cout << endl;
        if( printsam_over ) cout << "  sp5 scan  signalpar = " << signalpar << endl;
        for( size_t np=0; np<nfitp; np++)
        {
          double sp = py[np]/pp[np];
          ss[np] = yy-2.*sp*py[np]+sp*sp*pp[np];
          if( printsam_over ) cout << "  "<< sp <<"( " << ss[np]  << " ) ";
        }
        if( printsam_over ) cout << endl;
        std::vector<double> tcorr = fit_parabola(t,ss);
//  time1 keep unchanged
//  correct time TODO check  tcorr.second which is Time fit sigma
// ??        signalpar

// Ne po fitu i hardcoded tuftu gonju. nado schtat' v minimume.
// nado vse proschitat' i sravnit'. Pust' dolgo.

// Ignore firt results
         amp_over = py[3]/pp[3];

//         if(  fabs( tcorr[1] ) < 20.  )
//         {
//           amp_over = py[3]/pp[3];
//           time_over += tcorr[0];
//           timeerr_over = tcorr[1];
//           chisq_over = tcorr[2];
//         }
//         else // No corrections
//         {
//           chisq_over = ss[4];
//         }


        if( printsam_over ) cout <<" amp_over " << amp_over << " tcorr " << tcorr[0] <<" terr_over " << timeerr_over <<" chisq_over " << chisq_over << endl;
//       }   //   Fini with precautions that All exists
//       }
//     }
 }  // end neshapefit scope

//         if( debug )
//           {
//             cout << " NOverflow =  " << overflow_samples.size() << " Overflow? " << endl;
//             for( int is=0; is < (int)sample.size(); is++)
//               {
//                 cout <<" " << sample[is];
//               }
//             cout << endl;
//
//             cout << " Csample " << endl;
//             for ( int is=0; is < (int)sample.size(); is++ )
//               {
//                 cout << " " << (int)(10.*csample_[is]);
//               }
//             cout << endl;
//           }
//
//         vector <int> over;
//         unsigned int new_over(0);
//
//
//         // Calculate time
//         int plato_size = overflow_samples.back() - overflow_samples.front() + 1;
//
//         if( debug )
//           cout << " Overflow detected from sample["<< overflow_samples.front() <<"]="<<sample[overflow_samples.front()] <<
//             "to sample["<< overflow_samples.back() <<"]="<<sample[overflow_samples.back()] <<
//             " with plato size = "<< plato_size << endl;
//
//         double value_ower = SADC_shape_value_ower[shape_table_][plato_size];
//         double index_ower = SADC_shape_index_ower[shape_table_][plato_size];
//         double time = overflow_samples.front() + SADC_shape_max_position[shape_table_]
//           - SADC_shape_max_from_front[shape_table_] - index_ower;
//
//         if( debug )
//           cout << " SADC_shape_max_position = " << SADC_shape_max_position[shape_table_] <<
//             " SADC_shape_max_from_front = " << SADC_shape_max_from_front[shape_table_] <<
//             " from index_ower to front = " << SADC_shape_max_position[shape_table_]-SADC_shape_max_from_front[shape_table_]-index_ower << endl
//                << " Time ower in clocks before TCS correction = " << time << endl;
//
//         // Convert to ns
//         time *= sadc_clock_;
//
//         if( debug )
//           cout << " Over Value_ower " << value_ower << " index_ower " << index_ower << endl;
//
//         double ped = (ped_odd_old_+ped_even_old_)/2.;
//         double csignal = overflow_amplitude_ - ped;
//
//         if( debug )
//           cout <<" Csample value at first over point(calculated!) " << csignal << endl;
//
//         double signal = csignal/value_ower;
//
//         // Prepare for precise corrections, need more investigations, skip for a while Fri Feb 13 13:34:05 CET 2009
//         if( debug )
//           cout << " Value " << signal << endl
//                << "********************" << endl
//                << " Value + ped " << signal + ped <<  endl
//                << "--------------------" << endl;
//
//         // Set return structure
//         result_.ampl = signal;
//         result_.time = time;
//         result_.ready = 1;
//
//         // Return after overflow
// //        return result_;
//         return;

      }  // End  Method FitOverflow


  }

  double time = 0.;

  if( debug ) cout << " sample " << sample[0] <<" " << sample[1] <<" " << sample[2] <<" " << sample[3] <<" " << sample[4] <<" "<< endl;


  double signal = result_.ampl;
  int max_position = result_.max_pos;

// Commented as I need to clarify origin and usage of this data member: front_position_
//  int front_position = front_position_;


  // calculate chi square
//   int ishape0 = (int)(max_position - SADC_shape_max_position[shape_table_]);
  double max = signal;
  if( debug ) cout << " max " << max << endl;
//  double chisq = 0.;
  const double ped_ = result_.base + result_.base_diff / 2.;
// Drop any Shape usage, drop any   chisq usage!!
//   for ( int is=0; is < (int)sample.size(); is++ )
//     {
//       double f = Vshape( is-ishape0, max, ped_ );
//       double r = csample_[is];
//       double dvt = DVshape( is-ishape0, max, ped_ );
//       dvt = dvt*dvt;
//       double d =  r - f;
//       chisq += d*d/(1.+dvt);
//     }

//   if( debug )
//     {
//       cout << " Csample " << endl;
//       for ( int is=0; is < (int)sample.size(); is++ )
//         {
//           cout << " " << (int)(10.*csample_[is]);
//         }
//       cout << endl;
//       cout << " Vshape " << endl;
//       for ( int is=0; is < (int)sample.size(); is++ )
//         {
//           cout << " " << (int)(10.*Vshape( is-ishape0, max, ped_ ));
//         }
//       cout << endl;
//       cout << " Initial chisq " << chisq << " prob " << TMath::Prob( chisq, 31)  << endl;
//       }

// // Drop StraightLine usage, drop any   chisq of NiseHypot. usage!!
//   { // Calculate streight line chi-square to compare to shape
//     double chisq_line_new = 0.;
//     double time_line = 0.;
//     double signal_line = 0.;
//     {
//
//       double sam_s = 0.;
//       double sam_ss = 0.;
//       double sam_sf = 0.;
//       double sam_f = 0.;
//       double sam_ff = 0.;
//
//       double sam_n = (double)csample_.size();
//       for ( unsigned int i=0; i < csample_.size(); i++ )
//         {
//           double s = (double)(csample_[i]);
//           double f = (double)(i);
//           sam_s += s;
//           sam_ss += s*s;
//           sam_sf += s*f;
//           sam_f += f;
//           sam_ff += f*f;
//         }
//       sam_s /= sam_n;
//       sam_ss /= sam_n;
//       sam_sf /= sam_n;
//       sam_f /= sam_n;
//       sam_ff /= sam_n;
//
//       double a_new = ( sam_sf - sam_f*sam_s )/( sam_ff - sam_f*sam_f );
//       double ped_new = sam_s -  a_new*sam_f;
//       for ( int i=0; i < (int)csample_.size(); i++ )
//         {
//           double d =  (double)((double)(csample_[i])-ped_new ) - a_new*double(i);
//           chisq_line_new += d*d;
//         }
//       if( a_new >= 0. )
//         {
//           time_line = (double)(csample_.size());
//           signal_line = a_new*(double)(csample_.size());
//           time_line *= sadc_clock_;
//           //       double tcs_cor = CsEvent::Instance()->getTCSPhaseTime()-TCS_T0;
//           //       if( MAKE_TCS_CORRECTIONS ) time_line += tcs_cor;
//         }
//       else
//         {
//           time_line = -2.;
//           //      signal_line = ped_new + a_new*time_line;
//           signal_line = -a_new*(double)(csample_.size());
//           time_line *= sadc_clock_;
//           //       double tcs_cor = CsEvent::Instance()->getTCSPhaseTime()-TCS_T0;
//           //       if( MAKE_TCS_CORRECTIONS ) time_line += tcs_cor;
//         }
//       if( debug )   cout << " chisq_line_new " << chisq_line_new << endl;
//     }
//
//     if( chisq_line_new < chisq )
//       {
//         result_.time = time_line;
//         result_.ampl = signal_line;
//         result_.ready = 1;
//         if( debug )
//           {
//             cout << " Staight line fit looks better then Shape fit " << endl;
//             cout << " Considering this amlitude as Garbaje we need to return something Thinking ... " << endl;
//             cout << " Signal out of range! MaxAdvanced signal " << GetSignal() <<
//               " pedestal " << GetPed() << " time(~ in clocks) " << (int)(GetTime()/sadc_clock_) <<
//               " max " << max_position <<
//               " Xcheck " <<  sample[max_position]-ped_ <<" Sample size " << sample.size() << endl;
//           }
// //        return result_;
//         return;
//       }
//   }

//   if( debug ) cout << " MaxAdvanced signal " << signal << " max_position " << max_position <<
//                 " front_position " << front_position <<" diff " << max_position-front_position-0.5 << endl;
//
//   double expected_front_position = max_position - SADC_shape_max_from_front[shape_table_];
//
//   if( debug ) cout << " front_position " << front_position  <<
//                 " expected front_position " << max_position - SADC_shape_max_from_front[shape_table_] <<
//                 " sadc_front_gate " << SADC_shape_sadc_front_gate[shape_table_] << endl;
//
//   // Start of signal fits expected start of signal
//   if( fabs(float(expected_front_position - front_position)) < SADC_shape_sadc_front_gate[shape_table_] )
//     {
//
// /// \todo  TODO: remove old time calculation
//       // Not used at the moment
//       int front_min = front_position - SADC_shape_sadc_front_gate[shape_table_] - 1;
//       int front_max = front_position + SADC_shape_sadc_front_gate[shape_table_] + 1;
//
//     }

  // Start of signal does not match expectations
/// \todo  TODO calculate time for mismatch between expected and real signal time
//  Skip this method as refereing to some Shape parameters
//  time = CalcTime( HalfMax, sample, (unsigned int)max_position);

  std::pair<double,double> dtime2 = CalcTime2_halfMax( sample, (unsigned int)max_position);
  double time1 = dtime2.first;
  time = time1;  // Replace time calculated by the old method by new time calculation

  double time2 = dtime2.second;
// Dont care
  double max_time = double(max_position)*sadc_clock_;
//  cout << " time " << time <<  " time1 " << time1 << " time2 " << time2 <<" max time" << max_time << " FWHM " << time2-time1 << " HFW( HalfFrontWidth )  " << max_time-time1<< endl;
  if( !isfinite(time) )
    {
      cerr <<" Unexpected!! NAN is NOT(yet?) used in Any of  DigitizerSADCN calculations " << endl;
      exit(1);
//       time = expected_front_position;
//       // GoTo ns
//       time *= sadc_clock_;
    }


  double sxx=0.;
  double sy=0.;
//      int maxrng = 2;
      int maxrng = 1;
  double fmax = csample_[max_position];
  for( int is=max_position-maxrng; is <= max_position+maxrng; is++)
    {
      if( is <= 0 ) continue;
      if( is >= (int)sample.size() ) break;
      double dx = is-max_position;
      sxx += dx*dx;
      double dy = csample_[is]-fmax;
      sy += dy;
    }

  double a = 0.;
  if( sxx > 0 ) a = sy/sxx;

  if( debug ) cout << " MaxAdvanced a in max " << a << " max[" <<max_position<<"] = "<< csample_[max_position] << endl;
  if( a < 0 )
    {
      double ss = 0.;
      double sx = 0.;
      double sy = 0.;
      double syx = 0.;
      double sxx = 0.;
      double sxxx = 0.;

      // Recalculate signal
      for( int is=max_position-maxrng; is <= max_position+maxrng; is++)
        {
          if( is <= 0 ) continue;
          if( is >= (int)sample.size() ) break;
          double dx = is-max_position;
          double dx2 = dx*dx;
          double dy = csample_[is]-fmax;
          if( debug ) cout << " dx " << dx <<" dy " << dy << endl;
          ss += 1.;
          sx += dx;
          sy += dy;
          syx += dy*dx;
          sxx += dx2;
          sxxx += dx2*dx;
        }

      if( ss <= 0 )
        {
          cerr << "  EROOR EROOR EROOR EROOR MaxAdvanced internal error ss = " << ss << endl;
        }
      else
        {
          sx /= ss;
          sy /= ss;
          syx /= ss;
          sxx /= ss;
          sxxx /= ss;
          double b =  ( syx - a*sxxx - sy*sx + a*sxx*sx )/(sxx-sx*sx);
          double c =  sy - a*sxx - sx*b;
          double xmax = -b/(2*a);
          double ymax = c- a*xmax*xmax;
          if( debug ) cout << " MaxAdvanced xmax = " << xmax << " ymax = " << ymax << endl;
          if( debug )
            {
              cout << " X check max_position = " << max_position << " fmax " << fmax <<
                " X check1 " << csample_[max_position] << " X check2 signal " << signal << endl;
            }

          if( xmax < -2 ||  xmax > 2 || ymax < -1. || ymax > 0.3*signal || a > - 0.19 )
            {
              if ( debug ) {
                cout << " Not expected: MaxAdvanced xmax = " << xmax << " ymax = " << ymax <<" a= " << a << endl;
                cout << " X check max_position = " << max_position << " fmax " << fmax <<
                  " X check1 " << csample_[max_position] << " X check2 signal " << signal << endl;
              }
            }
          else
            {
              if( debug )
                {
                  cout << " Expected: MaxAdvanced xmax = " << xmax << " ymax = " << ymax <<" a= " << a <<  endl;
                  cout << " X check max_position = " << max_position << " fmax " << fmax <<
                    " X check1 " << csample_[max_position] << " X check2 signal " << signal <<" xmax " << xmax << endl;
                  for( int is=max_position-2; is <= max_position+2; is++)
                  {
                    if( is <= 0 ) continue;
                    if( is >= (int)sample.size() ) break;
                    cout <<"  v [ " << is - max_position << " ] " << csample_[is];
                  }
                  cout << endl;
                }
//                 if(signal > 100.)
//                 {
//                   cout << " Expected: MaxAdvanced xmax = " << xmax << " ymax = " << ymax <<" a= " << a <<  endl;
//                   cout << " X check max_position = " << max_position << " fmax " << fmax <<
//                     " X check1 " << csample_[max_position] << " X check2 signal " << signal <<" xmax " << xmax << endl;
//                   for( int is=max_position-2; is <= max_position+2; is++)
//                   {
//                     if( is <= 0 ) continue;
//                     if( is >= (int)sample.size() ) break;
//                     cout <<"  v [ " << is - max_position << " ] " << csample_[is];
//                   }
//                   cout << endl;
//                 }
                signal += ymax;
                max_time += xmax*sadc_clock_;
              }
        }
    }

// The Shape is NOT yet ready here!

// Try development with Pulse Shape
 double timeerr = 0.;
 double chisq = 0.;
 int maxgate = 5;
 {  // begin neshapefit scope

// Prepare input information to use tables
//  calib_fwhm_tab_
//  calib_hfwn_tab_

      bool printsam = false;
//      if( signal> 200. ) printsam = true;
//    if( GetParentName() == "EC02P1__" && signal> 1000. ) printsam = true;

      if( printsam ) cout <<" parent " << GetParent()  << "  "<< GetParentName() << " fwhm " << calib_fwhm_tab_ <<  " hfwn " << calib_hfwn_tab_ <<endl;
//     if( GetParent() != NULL   )
//     {
//       const CalorimeterMN * c = dynamic_cast<const CalorimeterMN *> (GetParent());
//       if( c != NULL )
//       {

// This must be hidden
//      const ManagerShapeTableSADC *mtab = c->GetManShapeTableSADC();
      const ManagerShapeTableSADC *mtab = mtab_;
//       if( mtab != NULL )
//       {
//   All exists

        if( printsam ) cout << GetParentName() <<" sample  signal " << signal <<" time " << time << endl;

        std::vector<double> pp;
        std::vector<double> py;
        std::vector<double> t;
        std::vector<double> ss;
        size_t nfitp = 6;
        double stept = 0.5; // can be amplitude/energy dependent
        for( size_t np=0; np<nfitp; np++)
        {
          pp.push_back(0.);
          py.push_back(0.);
          ss.push_back(0.);
          t.push_back( stept*(-3.+double(np))  );
        }

        double timepar = time1;
        double signalpar = signal;
        double yy =0.;

        for( int is=max_position-maxgate; is <= max_position+maxgate; is++)
        {
          if( is <= 0 ) continue;
          if( is >= (int)sample.size() ) break;
          double timex = double(is)*sadc_clock_-timepar;
//          double vtab = GetTabShapeValue( timex );  // we aim to relace by this call

          double x = timex/calib_fwhm_tab_;
          double vtab = mtab->Value( calib_hfwn_tab_, x );
          double dy = csample_[is]-signalpar*vtab;
          double y = csample_[is];
          double sigy = 1.+0.01*fabs(y);
          double w = 1./sigy/sigy;
          yy += y*y*w;
          for( size_t np=0; np<nfitp; np++)
          {
             double xp = x + t[np]/calib_fwhm_tab_;
             double p = mtab->Value( calib_hfwn_tab_, xp );
             pp[np] += p*p*w;
             py[np] += p*y*w;
          }
          if( printsam ) cout <<" x " << x  << " V " <<   csample_[is] <<"( " << signal*vtab <<") "<<" dy " << dy;
        }

        if( printsam ) cout << endl;
        if( printsam ) cout << "  sp5 scan  signal = " << signal << endl;
        for( size_t np=0; np<nfitp; np++)
        {
          double sp = py[np]/pp[np];
          ss[np] = yy-2.*sp*py[np]+sp*sp*pp[np];
          if( printsam ) cout << "  "<< sp <<"( " << ss[np]  << " ) ";
        }
        if( printsam ) cout << endl;
        std::vector<double> tcorr = fit_parabola(t,ss);
//  time1 keep unchanged
//  correct time TODO check  tcorr.second which is Time fit sigma
// ??        signalpar

// Ne po fitu i hardcoded tuftu gonju. nado schtat' v minimume.
// nado vse proschitat' i sravnit'. Pust' dolgo.


        if(  fabs( tcorr[1] ) < 20.  )
        {
          signal = py[3]/pp[3];
          time += tcorr[0];
          timeerr = tcorr[1];
          chisq = tcorr[2];
        }
        else // No corrections
        {
          chisq = ss[4];
        }


        if( printsam ) cout <<" tcorr " << tcorr[0] <<" terr " << timeerr <<" chisq " << chisq << endl;
//       }   //   Fini with precautions that All exists
//       }
//     }
 }  // end neshapefit scope

 if( debug_over && new_over  )
 {
    cout <<" signal " << signal << " time " << time << endl;
 }
 // Prepare result
//  cout <<" preliminary signal " << result_.ampl  << " final signal " << signal <<endl;
  result_.ampl = signal;
  result_.time = time;
  result_.timeFWHM = time2-time1;
  result_.timeHFW = max_time-time1; // should be smeared by TCSphase
  result_.time1 = time1;
  result_.time2 = time2;
  result_.ready = 1;
  result_.chi2 = chisq;
  result_.ndf = maxgate*2 + 1 - 2;
  if( new_over  )
  {
     result_.nover = new_over;
     result_.ampl_overflow = amp_over;
     result_.time_overflow = time_over;

     if( debug_over ) cout <<" OverflowDetected " << OverflowDetected() << endl;
  }
  else
  {
     result_.nover = new_over;
     result_.ampl_overflow = amp_over;
     result_.time_overflow = time_over;
  }
// Is the Shape ready here? Yes, it is. result_ must be filled in. TODO:: For Many key variables( like FWHM and HFW) we must be refer to result itself but not digitizer!
  is_noise_by_shape_anal_   = BadShape();

  if( debug && signal > signal_cut2print4debug )
    cout << " MaxAdvanced signal " << GetSignal() <<
      " pedestal " << GetPed() << " time(~ in clocks) " << (int)(time/sadc_clock_) <<" max " << max_position <<
      " Xcheck " <<  sample[max_position]-ped_ <<" Sample size " << sample.size() << endl;
  if( debug ) cout << " MaxAdvanced time " << GetTime() << endl;

//  return result_;
  return;
}

/////////////////////////////////////////////////////////

std::pair<double, double> DigitizerSADCN::CalcTime2_halfMax(const vector<short unsigned int> &sample, const unsigned int max_position)
 {
//    bool debug = true;
   bool debug = false;
// No assumptions on SADC shape
  if( debug )
  {
    cout <<" CalcTime2_halfMax " << max_position  << endl;
    for( int is =0;  is < (int)sample.size(); is++)
    {
      if( is == (int)max_position ) cout <<"  !";
      cout <<"  "<<  csample_[is];
      if( is == (int)max_position ) cout <<"  !";
    }
    cout << endl;
  }


// !!!  Otkuda eto -0.5 ????? Revision!!!!!!! Ili tak nado?
// Vrode-by eto tak nado. Esli max u nas samyi levyi bin, to my predpolagaem, chto levee
// naxoditsya 0, ni i time togda dolzhen byt' -0.5
// Analogichno s zadnim frontom time budet +0.5

  double time = max_position-0.5;

  double fmax = csample_[max_position];
  for( int is = max_position-1;  is >= -1; is--)
  {
    if( is >= 0 )
    {
      if( csample_[is] <= fmax/2. )
      {
        time = is + (fmax/2.- csample_[is])/(csample_[is+1]- csample_[is]);
        break;
      }
    }
    else
    {
      if( debug ) cout << " time " << time <<" expected " << -0.5 << endl;
      break;
    }
  }

  double time2 = max_position+0.5;
  for( int is = max_position+1;  is <= (int)sample.size(); is++)
  {
    if( is < (int)sample.size() )
    {
      if( csample_[is] <= fmax/2. )
      {
// error!        time2 = is + (fmax/2.- csample_[is])/(csample_[is-1] - csample_[is]);
        time2 = is - (fmax/2.- csample_[is])/(csample_[is-1] - csample_[is]);
        break;
      }
    }
    else
    {
      if( debug ) cout << " time2 " << time2 <<" expected " << (int)sample.size()-0.5 << endl;
      break;
    }
  }
  if( debug ) cout <<" Time2 ( bins) " << time << " max_position " << max_position << " time2 " << time2 <<" sample size " << sample.size() << endl;
//  cout << " What we get so far? max_position " << max_position << " result_.ampl = " << result_.ampl  << " sample max position " << csample_[max_position] << endl;
// Try to estimate time_max from parabola fit for +/- 1(2?) bins around maximum
// Ne zdes' Sorry
// end of time_max detection

  time *= sadc_clock_;
  time2 *= sadc_clock_;

  return std::pair< double, double> (time,time2);
}

//////////////////////////////////////////////////////////

bool   DigitizerSADCN::BadShape ( void ) const
{
  return !Fit2dcutFWHMnHFW( 1.2 );
}

//////////////////////////////////////////////////////////

double   DigitizerSADCN::GetRadiusFWHMnHFW   ( void ) const
{
  double x =( GetFWHMTimeNorm() -parFWHMn_HFW_fwhmn_c_ )/parFWHMn_HFW_fwhmn_r_;
  double y =( GetHFWTimeNorm() -parFWHMn_HFW_hfw_c_ )/parFWHMn_HFW_hfw_r_;
  double r = sqrt( x*x + y*y);
  return r;
}

//////////////////////////////////////////////////////////

bool   DigitizerSADCN::Fit2dcutFWHMnHFW   ( double rcut ) const
{
  double r = GetRadiusFWHMnHFW();
  return (r<rcut);
}

//////////////////////////////////////////////////////////

double   DigitizerSADCN::GetSigmaFWHMTime ( void ) const
{
//  return calib_fwhm_sigma_tab_;

//  cout <<" result_.ampl " << result_.ampl << endl;
  double e = 0.05*result_.ampl;
  double factor = ( 0.02*0.02 + 0.07*0.07/e);
  factor = sqrt(factor);
  return calib_fwhm_tab_*factor;
}

//////////////////////////////////////////////////////////

double   DigitizerSADCN::GetFWHMTimeNorm ( void ) const
{
  double fwhmn =  (GetFWHMTime()-calib_fwhm_tab_)/GetSigmaFWHMTime();
  return fwhmn;
}

//////////////////////////////////////////////////////////

double   DigitizerSADCN::GetHFWTimeNorm    ( void ) const
{
  double hfwn = GetHFWTime()/calib_fwhm_tab_-calib_hfwn_tab_;
  return hfwn;
}

////////////////////////////////////////////////////////////////////////////////

double ManagerShapeTableSADC::Value( double hfw, double x ) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout <<" x " << x << " vtab_.size() " << vtab_.size() << endl;
  if(vtab_.size() == 0 ) return 0.;
  const std::map<double,double> * tab = GetTab( hfw );
  if( debug ) cout <<" tab * " << tab << endl;
  if(  tab != NULL) return ValueTab(x, tab);
  return 0.;
}

////////////////////////////////////////////////////////////////////////////////

std::pair< double, double> ManagerShapeTableSADC::WValue( double hfw, double x ) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout <<" x " << x << " vtab_.size() " << vtab_.size() << endl;
  if(vtab_.size() == 0 ) return std::pair< double, double>(0.,0.);
  const std::map<double,std::pair<double,double> > * tab = GetWTab( hfw );
  if( debug ) cout <<" tab * " << tab << endl;
  if(  tab != NULL) return WValueTab(x, tab);
  return std::pair< double, double>(0.,0.);
}

////////////////////////////////////////////////////////////////////////////////

const std::map<double,double> * ManagerShapeTableSADC::GetTab( double v ) const
{
//   if( vtab_.size() == 0 ) return NULL;
//   if( vtab_.size() == 1) return vtab_.begin()->second;
//   std::map<double, std::map<double,double> *>::const_iterator hi = vtab_.lower_bound(v);
//  // variable is larger than highest data point?
//   if ( hi == vtab_.end() )
//     return vtab_.rbegin()->second;
//
//  // variable is equal or smaller than the lowest data point?
//   if ( hi == vtab_.begin() )
//     return vtab_.begin()->second;
//   std::map<double, std::map<double,double> *>::const_iterator lo = hi;
//   lo--;
//   return lo->second;
// Return tab directly
  if( vtab_.size() == 0 ) return NULL;
  if( vtab_.size() == 1) return vtab_.begin()->second.pmap_;
  std::map<double, ProfileMaps >::const_iterator hi = vtab_.lower_bound(v);
 // variable is larger than highest data point?
  if ( hi == vtab_.end() )
    return vtab_.rbegin()->second.pmap_;

 // variable is equal or smaller than the lowest data point?
  if ( hi == vtab_.begin() )
    return vtab_.begin()->second.pmap_;
  std::map<double, ProfileMaps>::const_iterator lo = hi;
  lo--;
  return lo->second.pmap_;
}
////////////////////////////////////////////////////////////////////////////////

const std::map<double, std::pair<double,double > > *  ManagerShapeTableSADC::GetWTab( double v ) const
{
  if( vtab_.size() == 0 ) return NULL;
  if( vtab_.size() == 1) return vtab_.begin()->second.wmap_;
  std::map<double, ProfileMaps >::const_iterator hi = vtab_.lower_bound(v);
 // variable is larger than highest data point?
  if ( hi == vtab_.end() )
    return vtab_.rbegin()->second.wmap_;

 // variable is equal or smaller than the lowest data point?
  if ( hi == vtab_.begin() )
    return vtab_.begin()->second.wmap_;
  std::map<double, ProfileMaps>::const_iterator lo = hi;
  lo--;
  return lo->second.wmap_;
}

// ////////////////////////////////////////////////////////////////////////////////
//
// const std::map <double,std::pair<double,double > > * ManagerShapeTableSADC::GetTab( double v ) const
// {
//   if( vtab_.size() == 0 ) return NULL;
//   if( vtab_.size() == 1) return vtab_.begin()->second.wmap_;
//   std::map<double, ProfileMaps >::const_iterator hi = vtab_.lower_bound(v);
//  // variable is larger than highest data point?
//   if ( hi == vtab_.end() )
//     return vtab_.rbegin()->second.wmap_;
//
//  // variable is equal or smaller than the lowest data point?
//   if ( hi == vtab_.begin() )
//     return vtab_.begin()->second.pmap_;
//   std::map<double, ProfileMaps>::const_iterator lo = hi;
//   lo--;
//   return lo->second.wmap_;
// }

////////////////////////////////////////////////////////////////////////////////
// Static functions

double ManagerShapeTableSADC::ValueTab( double v, const std::map<double,double> * tab ) const
{
  if( tab->size() == 0 ) return 0.;
  if( tab->size() == 1) return tab->begin()->second;
  std::map<double, double>::const_iterator hi = tab->lower_bound(v);
 // variable is larger than highest data point?
  if ( hi == tab->end() )
    return tab->rbegin()->second;

 // variable is equal or smaller than the lowest data point?
  if ( hi == tab->begin() )
    return tab->begin()->second;
  std::map<double, double>::const_iterator lo = hi;
  lo--;
 // linear interpolation
  return lo->second
        + (hi->second - lo->second) * (v - lo->first) / (hi->first - lo->first);
}

////////////////////////////////////////////////////////////////////////////////

std::pair<double,double>  ManagerShapeTableSADC::WValueTab( double v, const std::map<double, std::pair<double,double> > * tab ) const
{
  if( tab->size() == 0 ) return std::pair<double,double>(0.,0.);
  if( tab->size() == 1) return tab->begin()->second;
  std::map<double, std::pair<double,double> >::const_iterator hi = tab->lower_bound(v);
 // variable is larger than highest data point?
  if ( hi == tab->end() )
    return tab->rbegin()->second;

 // variable is equal or smaller than the lowest data point?
  if ( hi == tab->begin() )
    return tab->begin()->second;
  std::map<double, std::pair<double,double>  >::const_iterator lo = hi;
  lo--;
 // linear interpolation
  const std::pair<double,double> &lo2 = lo->second;
  const std::pair<double,double> &hi2 = hi->second;
  double ff = (v - lo->first) / (hi->first - lo->first);
  std::pair <double,double> ret = lo2;
  ret.first += (hi2.first - lo2.first)*ff;
  ret.second += (hi2.second - lo2.second)*ff;
  return ret;
}

///////////////////////////////////////////////////////////////////////////////

bool ManagerShapeTableSADC::CheckTables( void ) const
{
  cout <<" ManagerShapeTableSADC::CheckTables just Checking the size for the moment NTABLES = " << vtab_.size() << endl;
  return ( vtab_.size() > 0 );
}

////////////////////////////////////////////////////////////////////////////////

bool ManagerShapeTableSADC::AddTable( double hfwmin, const std::vector <int> &profile )
{
  bool debug = false;
  bool print_warnings = false;
  std::map<double,double> * tab = new std::map<double,double> ;
  int max=0;
  int itmax = -1;
  for( size_t it=0; it<profile.size(); it++)
  {
    if(profile[it] > max )
    {
      max = profile[it];
      itmax = it;
    }
  }

  if( max < 9000 )
  {
    cerr <<" ManagerShapeTableSADC::AddTable hfwmin " << hfwmin <<" Not expected max = " << max <<" : drop table creation " << endl;
    exit(1);
  }

  double xvalmax =   Xmin + double(itmax)*stepX_;
  for( size_t it=0; it<profile.size(); it++)
  {
     double x = Xmin + double(it)*stepX_;
     double v = double(profile[it])/double(max);
     if( v == 0. ) continue;
//     cout <<" x " << x << " v " << v << endl;
     if( !tab->insert(std::make_pair(x,v)).second )
     {
       cerr <<" ManagerShapeTableSADC::AddTable hfwmin " << hfwmin <<" insertion failed at x = " << x <<" : drop table creation " << endl;
       exit(1);
     }
  }


  std::map<double,std::pair<double, double> > * wtab = new std::map<double,std::pair<double, double> > ;
// Create Table of widths
  for( size_t iv=0; iv< 100; iv++)
  {
     double tagyval =(0.5+double(iv))/100.;
     double xval = Xmin;
     size_t ix=0;
     double ymin = 0.;
     double ymax = 0.;
     double xmin = Xmin;
     double xmax = Xmin;
     bool minok = false;
     bool maxok = false;
     for( ix=0; ix<10000; ix++)
     {
        xval = Xmin+0.005* double(ix);
        double y =  ValueTab( xval,  tab);
//        cout <<" xval " << xval << "   Y= " << y << endl;
        if( minok )
        {
          if( y < tagyval )
          {
            if( xval < xvalmax ) continue;

            if( maxok )
            {
              cerr <<" Create Table of widths internal error " << endl;
              exit(1);
            }
            else
            {
              maxok = true;
              ymax = y;
              xmax = xval;
              break;
            }
          }
        }
        else
        {
          if( y > tagyval )
          {
            minok = true;
            ymin = y;
            xmin = xval;
          }
        }
     }


     if( debug) cout <<" tagval " << tagyval << " xmin " << xmin <<" ymin " << ymin <<" xmax " << xmax <<" ymax " << ymax << " Width " << xmax - xmin << endl;
     double width = xmax - xmin;
     std::pair<double, double> val( tagyval, xmin);
     std::pair<double, std::pair<double, double> > widval( width,val);

     if( wtab->insert(widval).second )
     {
       if( print_warnings) cerr <<" Warning SADC width map filling width "<< width <<" was inserted already " << endl;
     }
  }

  ProfileMaps tabs(tab,wtab);

  if( !vtab_.insert( std::make_pair(hfwmin,tabs) ).second )
  {
     cerr <<" ManagerShapeTableSADC::AddTable hfwmin " << hfwmin <<" insertion vtab failed: drop table creation " << endl;
     exit(1);
  }

//  cout << " exit4debug " << endl;
//  exit (0);
  return true;
}

////////////////////////////////////////////////////////////////////////////////

int ManagerShapeTableSADC::InputShapeTableSADC(const string &s)
{
//   bool debug = true;
  bool debug = false;
  istringstream is(s);
  string str;

//  Read line of comments
  getline(is,str);
  int ret;
//  Read input table by table
// Get shapetable name together with calo name
  while( getline(is,str) )
  {
    if( debug )  cout <<" New Table " << str << endl;
// Get table parameters
    getline(is,str);
    float hfwbinmin, hfwbinmax, entries, xmin, xmax;
    int ncx;
    ret = sscanf(str.c_str()," %f  %f  %f  %f  %f  %d ", &hfwbinmin, &hfwbinmax, &entries, &xmin, &xmax, &ncx);
    if( debug )  cout <<" hfwbinmin " << hfwbinmin << " hfwbinmax " << hfwbinmax <<" entries " << entries <<" xmin " << xmin <<" xmax " << xmax  <<" ret " << ret << endl;
    assert( ret == 6 );
    if( ncx != 600 )
    {
      cerr <<" ManagerShapeTableSADC::InputShapeTableSADC format error not expected ncx = " << ncx << endl;
      exit(1);
    }
    double step = (xmax-xmin)/double(ncx);
    if( debug ) cout <<" step " << step << endl;
//    ret = sscanf(str.c_str(),"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f%f %f %f %f %f ", &);
    assert( ret == 6 );
    std::vector<int> profile(600);
    if( debug ) cout <<" profile " << endl;
    int max=0;
    for(size_t il=0; il<30; il++)
    {
      getline(is,str);
      istringstream iis(str);
      for(size_t it=0; it<20; it++)
      {
        iis >> profile[it+20*il];
        if( debug ) cout <<" [ " <<it+20*il <<"] " << profile[it+20*il];
        if(profile[it+20*il] > max ) max =  profile[it+20*il];
      }
      if( debug ) cout << endl;
    }
// Read sigma separator
    getline(is,str);
// Read sigma
    std::vector<int> sigma(600);
    if( debug ) cout <<" sigma " << endl;
    for(size_t il=0; il<30; il++)
    {
      getline(is,str);
      istringstream iis(str);
      for(size_t it=0; it<20; it++)
      {
        iis >> sigma[it+20*il];
        if( debug ) cout <<" [ " <<it+20*il <<"] " << sigma[it+20*il];
      }
      if( debug ) cout << endl;
    }
// Read stat separator
    getline(is,str);
// Read stat
    std::vector<int> stat(600);
    if( debug ) cout <<" stat  "<< endl;
    for(size_t il=0; il<30; il++)
    {
      getline(is,str);
      istringstream iis(str);
      for(size_t it=0; it<20; it++)
      {
        iis >> stat[it+20*il];
        if( debug ) cout <<" [ " <<it+20*il <<"] " << stat[it+20*il];
      }
      if( debug ) cout << endl;
    }
// Store table in persistent object
    bool ok = AddTable( hfwbinmin,  profile );
    if( debug ) cout <<" AddTable ok " << ok << endl;
  }
  bool oktables = CheckTables();
  if( debug ) cout <<" No more tables oktables " << oktables << endl;

  return 0;
}

//////////////////////////////////////////////////////////
// Destructor

DigitizerSADCN::~DigitizerSADCN  (void) {
}

//////////////////////////////////////////////////////////

//} // namespace MN

////////////////////////////////////////////////////////////////////////////////

