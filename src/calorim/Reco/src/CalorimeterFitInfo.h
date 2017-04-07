#ifndef inside_Calorimeter_h
#  error Only for inclusion from inside Calorimeter.h!
#endif


/*! \brief  Collection of fit parameters and validation flags
  used to store fit results for each cell
*/
class FitInfo
{
 public:
  /// Default constructor
  FitInfo       (void) : fit_ok(false) {}
  /// Copy constructor
  FitInfo  (const FitInfo &fi) {*this=fi;}
  // =========================================================================
  // Operators
  // =========================================================================
 public:
  /// Assignment operator
  FitInfo         &operator =  (const FitInfo &fi){
    if( &fi!=this )
      {
	fit_ok= fi.fit_ok;
	fit_parameters = fi.fit_parameters;
      } return *this;}
  // ============================================================================
  // Attributes, data
  // ============================================================================
 public:
  bool                           fit_ok;
  std::vector<std::pair<double,double> >   fit_parameters;
};

