#ifndef inside_Calorimeter_h
#  error Only for inclusion from inside Calorimeter.h!
#endif


class SummaryInspectLED
{
 public:
  SummaryInspectLED( void ) { Clear(); }

  void Clear ( void )      {
    init_ = false;
    result_reported_ = false;
    stat_max_ = 0;
    no_entries_ = 0;
    low_statistic_ = 0;
    small_value_ = 0;
    cells_ok_ = 0;
  }

  void Print( void ) const {
    if( init_ )
      {
	if( result_reported_ )
	  std::cout << " Reported result ";
	else
	  std::cout << " Not Reported result ";

	std::cout << " SummaryInspectLED STATISTIC Total=" << stat_max_ <<
	  " Dead leds " << no_entries_ << " Low statistic " << low_statistic_ <<
	  " Small LEDs " << small_value_ <<
	  " NCells OK " << cells_ok_ << std::endl;
      }
    else
      {
	std::cout << " SummaryInspectLED Not Init " << std::endl;
      }
  }
 public:
  bool init_;
  bool result_reported_;
  int stat_max_;
  int no_entries_;
  int low_statistic_;
  int small_value_;
  int cells_ok_;
};
