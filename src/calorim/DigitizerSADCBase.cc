#include "DigitizerSADCBase.h"
#include "CsCalorimeter.h"
#include "CsDigitizerSADC.h"

#include <algorithm>
#include <numeric>

const double DigitizerSADCBase::sadc_clock_ = 12.86;
const double DigitizerSADCBase::TCS_T0      = 40.;

////////////////////////////////////////////////////////////////////////////////
DigitizerSADCBase::DigitizerSADCBase ( CsCalorimeter* parent,  Type type ) :
                  parent_(parent), type_(type),
                  debug_(false),
                  overflow_amplitude_(1023), is_noise_by_shape_anal_(false),
                  sadc_ped_max_(4),
                  mycell_(-1) {
                    // Set overflow acording to SADC type: SADC 10Bit, MSADC 12Bit
                    overflow_amplitude_ = type_ == CsDigitizerSADC::MSADC ? 4095 : 1023;
}

////////////////////////////////////////////////////////////////////////////////
std:: string DigitizerSADCBase::GetParentName( void ) const {
   if( parent_ == NULL ) return ""; return  parent_->GetName();
}

////////////////////////////////////////////////////////////////////////////////
void DigitizerSADCBase::PrintSettings ( void ) const {
  std::cout << " SADC::  Type " << type_;
  std::cout << " sadc_clock_ " << sadc_clock_;
  std::cout << " sadc_ped_max_ " << sadc_ped_max_;
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
void DigitizerSADCBase::ClearResult() {
  memset( &result_, 0, sizeof(result_t) );
}

////////////////////////////////////////////////////////////////////////////////
unsigned int DigitizerSADCBase::DetectOverflow ( const std::vector<uint16> &sample ) {
  unsigned int new_over(0);
  for( int is = 0; is < (int)sample.size(); is++) {
    if( sample[is] >=  overflow_amplitude_) {
      if( !new_over ) {
        new_over++;
        overflow_samples_.push_back(is);
      }
      else if( overflow_samples_.back() < is-1 ) {
        new_over++;
        overflow_samples_.push_back(is);
      }
      else {
        overflow_samples_.push_back(is);
      }
    }
  }
  return new_over;
}

////////////////////////////////////////////////////////////////////////////////
/*!
  \brief Generate corrected samples and make basic shape analysis
*/
void DigitizerSADCBase::GetCcSample ( const std::vector<uint16> &sample ) {
  bool debug = false;
  if( debug ) std::cout << " sample " << sample[0] <<
                 " " << sample[1] <<" " << sample[2] <<" " << sample[3] <<" " << sample[4] <<" "<< std::endl;

  ped_odd_old_ = 0.;
  ped_even_old_ =0.;
  double stat = 0.;
  for( int is=0; is <= sadc_ped_max_; is+=2 ) {
      ped_odd_old_ += sample[is];
      ped_even_old_ += sample[is+1];
      stat++;
  }
  if( stat > 0 ) {
      ped_odd_old_ /= stat;
      ped_even_old_ /= stat;
  }

  // Subtract pedestals
  for( int is=0; is < (int)sample.size(); is+=2) {
      csample_.push_back( double(sample[is]) - ped_odd_old_ );
      csample_.push_back( double(sample[is + 1]) - ped_even_old_ );
      if( debug ) std::cout << (int)(10.*csample_.back() ) <<" ";
  }

  // Find min, max, summ, slope
  double summ = std::accumulate(csample_.begin(), csample_.end(), 0.);
  int max_position = std::max_element(csample_.begin(), csample_.end()) - csample_.begin();
  int min_position = std::min_element(csample_.begin(), csample_.end()) - csample_.begin();
  double signal = csample_[max_position];
  double min = csample_[min_position];
  // We need this min_position to calculate delta correctly!

  result_.base = GetDynamicPed();
  result_.base_diff = (ped_odd_old_ - ped_even_old_);

  // Store some results, could be intermediate
  result_.ampl = signal;

  //  double time = max_position - SADC_shape_max_from_front[shape_table_];
  // We'd better calculate true time here??
  double time = max_position - 3.;
  time *= sadc_clock_;

  result_.time = time;
  // More technical info to output
  result_.dbase = result_.base_diff/2.;
  result_.max_pos = max_position;

  if( debug ) std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
void DigitizerSADCBase::PrintSample ( const std::vector<uint16> &sample ) {
  std::cout <<" Sample ";
  for( int i=0; i< (int)sample.size(); i++) {
    std::cout << sample[i] <<" ";
  }
  std::cout << std::endl;
}

void DigitizerSADCBase::PrintType (){
  std::cout<<"Base\n";
}
